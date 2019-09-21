'''
by Evan Mayer

Perform a frequency-switched measurement of spectra made from RTL-SDR samples.
Based on sample code from roger-'s pyrtlsdr library.
Integrate over a specified length of time and save off to a text file.
'''

# Import plotting necessities
import argparse
import numpy as np
from scipy.signal import welch
import time

from rtlsdr import RtlSdr


def run_fswitch( NFFT, gain, rate, fc, fthrow, t_int ):
    '''
    Inputs:
    NFFT:    Number of elements to sample from the SDR IQ timeseries: powers of 2 are most efficient
    gain:    Requested SDR gain (dB)
    rate:    SDR sample rate, intrinsically tied to bandwidth in SDRs (Hz)
    fc:      base center frequency (Hz)
    fthrow:  frequency to switch to, samples subtracted from spectrum of fc
             (Hz)
    t_int:   total integration time (s)

    Returns:
    freqs_fold:  Frequencies of the resulting spectrum, centered at fc (Hz), numpy array
    p_fold:   Folded frequency-switched spectrum (dB/Hz) numpy array,
                 following Herschel docs at http://herschel.esac.esa.int/hcss-doc-15.0/load/hifi_um/html/hcb_pfsw.html
    '''

    sdr = RtlSdr()
    sdr.rs = rate # Rate of Sampling (intrinsically tied to bandwidth with SDR dongles)
    sdr.fc = fc
    sdr.gain = gain
    print('  sample rate: %0.6f MHz' % (sdr.rs/1e6))
    print('  center frequency %0.6f MHz' % (sdr.fc/1e6))
    print('  gain: %d dB' % sdr.gain)

    # Set up arrays to store total power calculated from I-Q samples
    p_xx_on_tot = np.zeros(NFFT)
    p_xx_off_tot = np.zeros(NFFT)
    cnt_on = 0
    cnt_off = 0

    # Set the baseline time
    start_time = time.time()
    print('Integration began at {}'.format(time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime(start_time))))
    # Time integration loop
    # TODO: This could absolutely be sped up using a smart buffer, but this is python after all
    # If you really need speed, just use rtl_power_fftw, that's good C++ with a buffer

    # Since we essentially sample as long as we want on each frequency,
    # Estimate the power spectrum by Bartlett's method:
    while time.time()-start_time < t_int:
        # Switch to the center frequency
        sdr.fc = fc
        # Collect on-frequency samples
        iq_on = np.zeros(NFFT, dtype=complex)
        iq_on += sdr.read_samples(NFFT)
        if iq_on.any():
            cnt_on += 1

        # Switch to the switched frequency
        sdr.fc = fc + fthrow
        # Collect off-frequency samples
        iq_off = np.zeros(NFFT, dtype=complex)
        iq_off += sdr.read_samples(NFFT)
        if iq_off.any():
            cnt_off += 1
        
        #p_xx_on, freqs_on = psd(iq_on, NFFT=NFFT, Fs=rate, scale_by_freq=False)
        #p_xx_off, freqs_off = psd(iq_off, NFFT=NFFT, Fs=rate, scale_by_freq=False)

        # Following https://en.wikipedia.org/wiki/Bartlett%27s_method: 
        # Use scipy.signal.welch to compute 1 periodogram for each freq hop.
        # For non-overlapping intervals, which we have because we are sampling
        # the timeseries as it comes in, the welch() method is equivalent to
        # Bartlett's method.
        # We therefore have an N=NFFT-point data segment split up into K=1 non-
        # overlapping segments, of length M=NFFT.
        # This means we can call welch() on each set of samples from the SDR,
        # accumulate them, and average later by the number of hops on each freq
        # to reduce the noise while still following Barlett's method, and
        # without keeping huge arrays of iq samples around in RAM.
        freqs_on, p_xx_on = welch(iq_on, fs=rate, nperseg=NFFT, noverlap=0, scaling='spectrum', return_onesided=False)
        freqs_off, p_xx_off = welch(iq_off, fs=rate, nperseg=NFFT, noverlap=0, scaling='spectrum', return_onesided=False)
        p_xx_on_tot += p_xx_on
        p_xx_off_tot += p_xx_off
    
    end_time = time.time()
    print('Integration ended at {} after {} seconds.'.format(time.strftime('%a, %d %b %Y %H:%M:%S'), end_time-start_time))
    print('{} spectra were measured at {}.'.format(cnt_on, fc))
    print('{} spectra were measured at {}.'.format(cnt_off, fc+fthrow))

    # Unfortunately, welch() with return_onesided=False does a sloppy job
    # of returning the arrays in what we'd consider the "right" order,
    # so we have to swap the first and last halves to avoid artifacts
    # in the plot.
    half_len = len(freqs_on)//2
    # Swap frequencies:
    tmp_first = freqs_on[:half_len].copy() 
    tmp_last = freqs_on[half_len:].copy()
    freqs_on[:half_len] = tmp_last
    freqs_on[half_len:] = tmp_first

    tmp_first = freqs_off[:half_len].copy()
    tmp_last = freqs_off[half_len:].copy()
    freqs_off[:half_len] = tmp_last
    freqs_off[half_len:] = tmp_first

    # Swap powers:
    tmp_first = p_xx_on_tot[:half_len].copy()
    tmp_last = p_xx_on_tot[half_len:].copy()
    p_xx_on_tot[:half_len] = tmp_last
    p_xx_on_tot[half_len:] = tmp_first

    tmp_first = p_xx_off_tot[:half_len].copy()
    tmp_last = p_xx_off_tot[half_len:].copy()
    p_xx_off_tot[:half_len] = tmp_last
    p_xx_off_tot[half_len:] = tmp_first

    # Compute the average power spectrum based on the number of spectra read
    p_avg_on = 10.*np.log10(p_xx_on_tot / cnt_on)
    p_avg_off = 10.*np.log10(p_xx_off_tot / cnt_off)
    # Compute the difference spectrum
    p_diff = p_avg_on - p_avg_off

    # Shift frequency spectra back to the intended range
    freqs_on = freqs_on + fc
    freqs_off = freqs_off + fc + fthrow

    save_spectrum('on_avg.txt', freqs_on, p_avg_on)
    save_spectrum('off_avg.txt', freqs_off, p_avg_off)

    # Shift-and-add to fold
    # Magnitude of frequency shift in bin space:
    epsilon = (freqs_on[1] - freqs_on[0]) / 2.
    fc_idx = np.where(np.abs(freqs_on - fc) < epsilon)[0][0]
    fthrow_idx = np.where(np.abs(freqs_on - (fc-fthrow)) < epsilon)[0][0]
    bin_throw = np.abs(fthrow_idx - fc_idx)
    # Folding procedure is a shift-and-add-negate, then average
    p_fold = (p_diff[bin_throw:-1] - p_diff[0:NFFT-bin_throw-1]) / 2.
    # Get the folded, upper segment of freqs to return it,
    # throwing away the unfolded part of the spectrum. Is there a way around this?
    freqs_fold = freqs_on[bin_throw:-1]

    # nice and tidy
    sdr.close()

    return freqs_fold, p_fold


def save_spectrum(filename, freqs, p_xx):
    '''
    Save the results of integration to a file.
    '''
    header='\n\n\n\n\n'
    np.savetxt(filename, np.column_stack((freqs, p_xx)), delimiter=' ', header=header)
    print('Results were written to {}.'.format(filename))

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gain', dest='gain', default=30.0, action='store', type=float, help='Requested SDR gain (dB), default 30.0')
    parser.add_argument('-r', '--rate', dest='rate', default=2.56e6, action='store', type=float, help='Sample rate (Hz), default 2.56e6 Hz')
    parser.add_argument('-c', '--fc', dest='fc', default=1.4202e9, action='store', type=float, help='Requested base frequency (Hz), default = 1.420e9 Hz')
    parser.add_argument('-d', '--dfthrow', dest='fthrow', default=4.0e5, action='store', type=float, help='Delta to fc applied on every switching integration cycle (Hz), default 4.0e5 Hz')
    parser.add_argument('-t', '--time', dest='t_int', default=1, action='store', type=float, help='Total amount of time to integrate frequency switched measurements (s) default 1 s')
    parser.add_argument('-o', '--output', dest='output', default='rtl_freq_switch_int.tsv', action='store', type=str, help='Filename to save integrated spectrum to. Default is tab delimited, rtl_freq_switch_int.tsv.')

    args = parser.parse_args()

    NFFT=2**12

    freqs, p_xx_avg = run_fswitch( NFFT, args.gain, args.rate, args.fc, args.fthrow, args.t_int )
    save_spectrum( args.output, freqs, p_xx_avg )

