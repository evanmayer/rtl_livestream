'''
by Evan Mayer

Perform a frequency-switched measurement of spectra made from RTL-SDR samples.
Based on sample code from roger-'s pyrtlsdr library.
Integrate over a specified length of time and save off to a text file.
'''

# Import plotting necessities
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.mlab import psd
import time

from rtlsdr import RtlSdr


def run_fswitch( NFFT, gain, rate, fc, fthrow, fswitch, t_int ):
    '''
    Inputs:
    NFFT:    Pass-through for matplotlib plot psd: powers of 2 are most efficient
    gain:    Requested SDR gain (dB)
    rate:    SDR sample rate, intrinsically tied to bandwidth in SDRs (Hz)
    fc:      base center frequency (Hz)
    fthrow:  frequency to switch to, samples subtracted from spectrum of fc
             (Hz)
    fswitch: switching rate (Hz)
    t_int:   total integration time (s)

    Returns:
    freqs:   Frequencies of the resulting spectrum, centered at fc (Hz), numpy array
    p_xx_tot: ( Power spectral density on - power spectral density off ) / 2 (dB/Hz) numpy array
    '''

    sdr = RtlSdr()
    sdr.rs = rate # Rate of Sampling (intrinsically tied to bandwidth with SDR dongles)
    sdr.fc = fc
    sdr.gain = gain
    print('  sample rate: %0.6f MHz' % (sdr.rs/1e6))
    print('  center frequency %0.6f MHz' % (sdr.fc/1e6))
    print('  gain: %d dB' % sdr.gain)

    # The amount of time to integrate on each frequency
    delta_t = 1./(2.*fswitch)

    # The number of iq samples to read on each call to the RtlSdr
    NUM_READ_SAMPLES = NFFT #* 1000

    # Set up arrays to store integrated power from psd
    p_xx_on_tot = np.zeros(NFFT)
    p_xx_off_tot = np.zeros(NFFT)
    
    # Set the baseline time
    start_time = time.time()
    print('Integration began at {}'.format(time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime(start_time))))

    # Time integration loop
    while time.time()-start_time < t_int:
        # The start of the on-frequency cycles
        start_time_on = time.time()
        # Integrate on the target frequency
        iq_on = []
        while time.time()-start_time_on < delta_t:
            # Collect on-frequency samples
            iq_on.extend(sdr.read_samples(NUM_READ_SAMPLES))
        # Take the PSD of the on samples
        p_xx_on, freqs_on = psd(iq_on, NFFT=NFFT, Fs=rate/1e6)
        p_xx_on_tot += p_xx_on                                 

        # Switch to the switched frequency
        sdr.fc = fc + fthrow

        start_time_off = time.time()
        # Integrate on the switched frequency
        iq_off = []
        while time.time()-start_time_off < delta_t:
            # Collect off-frequency samples
            iq_off.extend(sdr.read_samples(NUM_READ_SAMPLES))
        # And the off samples
        p_xx_off, freqs_off = psd(iq_off, NFFT=NFFT, Fs=rate/1e6)
        p_xx_off_tot += p_xx_off      
            
        if not p_xx_on_tot.any() or not p_xx_off_tot.any():
            print('Warning: One of the frequency-switched spectra was empty. Consider lowering the switching frequency in order for RTL-SDR sampling to catch up.')
    end_time = time.time()
    print('Integration ended at {} after {} seconds.'.format(time.strftime('%a, %d %b %Y %H:%M:%S'), end_time-start_time))

    # Add folded samples to the current and update
    p_xx_tot = (p_xx_on_tot - p_xx_off_tot)/2.
    # Shift frequency spectrum back to the intended range
    freqs_on = freqs_on*1e6 + fc

    # nice and tidy
    sdr.close()

    return freqs_on, p_xx_tot


def save_spectrum(filename, freqs, p_xx):
    '''
    Save the results of integration to a file.
    '''
    header='\n\n\n\n\n'
    np.savetxt(filename, np.column_stack((freqs, p_xx)), delimiter=' ',)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gain', dest='gain', default=30.0, action='store', type=float, help='Requested SDR gain (dB), default 30.0')
    parser.add_argument('-r', '--rate', dest='rate', default=2.56e6, action='store', type=float, help='Sample rate (Hz), default 2.56e6 Hz')
    parser.add_argument('-c', '--fc', dest='fc', default=1.4202e9, action='store', type=float, help='Requested base frequency (Hz), default = 1.420e9 Hz')
    parser.add_argument('-d', '--dfthrow', dest='fthrow', default=4.0e5, action='store', type=float, help='Delta to fc applied on every switching integration cycle (Hz), default 4.0e5 Hz')
    parser.add_argument('-s', '--fswitch', dest='fswitch', default=0.2, action='store', type=float, help='Frequency at which SDR tuning frequency is switched from fc to fc+fthrow. One cycle is from fc to fc+fthrow bck to fc (Hz), default 0.2 Hz')
    parser.add_argument('-t', '--time', dest='t_int', default=1, action='store', type=float, help='Total amount of time to integrate frequency switched measurements (s) default 1 s')
    parser.add_argument('-o', '--output', dest='output', default='rtl_freq_switch_int.tsv', action='store', type=str, help='Filename to save integrated spectrum to. Default is tab delimited, rtl_freq_switch_int.tsv.')

    args = parser.parse_args()

    NFFT=2048

    freqs, p_xx_tot = run_fswitch( NFFT, args.gain, args.rate, args.fc, args.fthrow, args.fswitch , args.t_int )
    save_spectrum( args.output, freqs, p_xx_tot )

