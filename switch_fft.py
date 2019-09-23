'''
by Evan Mayer

Perform a noise source-switched measuremnt from RTL-SDR samples.
I use the nooelec H1 Barebones sawbird filter/LNA combo.
Uses the rtl_biast utility to switch the RTL-SDR blog v3 dongle's bias-tee
on and off, but these lines can be easily removed.
Depends on pyrtlsdr from roger-.
Requires the gpiozero module, which is installed by default in Raspbian.
Integrate over a specified length of time. Begin with a measurement of the
receiver switched noise source, then revisit it every 5 minutes afterward.
'''

import argparse
import gpiozero
import os
from rtlsdr import RtlSdr
import numpy as np
from scipy.signal import welch
import subprocess
import time


def run_integrate( NFFT, gain, rate, fc, t_int ):
    '''
    Inputs:
    NFFT:    Number of elements to sample from the SDR IQ timeseries: powers of 2 are most efficient
    gain:    Requested SDR gain (dB)
    rate:    SDR sample rate, intrinsically tied to bandwidth in SDRs (Hz)
    fc:      base center frequency (Hz)
    t_int:   total integration time (s)
                                                                                                                           
    Returns:
    freqs_on: Frequencies of the resulting spectrum, centered at fc (Hz), numpy array
    p_avg_on: Average power spectral density in dB/Hz
    '''

    # Ensure the bias tee is turned off, then turn it on.                       
    basepath = os.path.expanduser('~/scratch/')
    cmd = [os.path.join(basepath, 'rtl_biast', 'build', 'src', 'rtl_biast'), '-d 0', '-b 0']
    print('Disabling rtl_biast.')
    ret = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, capture_output=False)
    print('Enabling rtl_biast.')
    cmd[-1] = '-b 1'
    ret = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, capture_output=False)
                                                                                                                           
    print('Initializing rtl-sdr with pyrtlsdr.')
    sdr = RtlSdr()
    sdr.rs = rate
    sdr.fc = fc
    sdr.gain = gain
    print('  sample rate: %0.6f MHz' % (sdr.rs/1e6))
    print('  center frequency %0.6f MHz' % (sdr.fc/1e6))
    print('  gain: %d dB' % sdr.gain)

    # Set up arrays to store total power calculated from I-Q samples
    p_xx_on_tot = np.zeros(NFFT)
    cnt_on = 0
                                                                                                                           
    # Set the baseline time
    start_time = time.time()
    last_noise_time = start_time
    time_since_noise = 0.0;
    print('Integration began at {}'.format(time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime(start_time))))
    # Time integration loop
    # TODO: This could absolutely be sped up using a smart buffer, but this is python after all
    # If you really need speed, just use rtl_power_fftw, that's good C++ with a buffer
    # The last time I tried this, subprocess.run() had a lot of overhead for
    # spawning and waiting for child processes to end at the rates desired for
    # frequency switching with rtl-power-fftw.
                                                                                                                           
    # Since we essentially sample as long as we want on each frequency,
    # Estimate the power spectrum by Bartlett's method:
    iq_on = np.zeros(NFFT, dtype=complex)
    while time.time()-start_time < t_int:
        # Collect on-frequency samples
        iq_on = sdr.read_samples(NFFT)
        cnt_on += 1
                                                                                                                           
        # Following https://en.wikipedia.org/wiki/Bartlett%27s_method: 
        # Use scipy.signal.welch to compute 1 periodogram for each source.
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
        p_xx_on_tot += p_xx_on
    
    end_time = time.time()
    print('Integration ended at {} after {} seconds.'.format(time.strftime('%a, %d %b %Y %H:%M:%S'), end_time-start_time))
    print('{} spectra were measured at {}.'.format(cnt_on, fc))
                                                                                                                           
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
                                                                                                                           
    # Swap powers:
    tmp_first = p_xx_on_tot[:half_len].copy()
    tmp_last = p_xx_on_tot[half_len:].copy()
    p_xx_on_tot[:half_len] = tmp_last
    p_xx_on_tot[half_len:] = tmp_first

    # Interpolate across central point to remove DC spike
    p_xx_on_tot[half_len] = (p_xx_on_tot[half_len - 1] + p_xx_on_tot[half_len + 1]) / 2
                                                                                                                           
    # Compute the average power spectrum based on the number of measurements
    p_avg_on = 10.*np.log10(p_xx_on_tot / cnt_on )
                                                                                                                           
    # Shift frequency spectra back to the intended range
    freqs_on = freqs_on + fc
                                                                                                                           
    # nice and tidy
    sdr.close()

    print('Disabling rtl_biast.')
    cmd[-1] = '-b 0'
    ret = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, capture_output=False)
    
    return freqs_on, p_avg_on
    
    
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
    parser.add_argument('-t', '--time', dest='t_int', default=1, action='store', type=float, help='Total amount of time to integrate frequency switched measurements (s) default 1 s')
    parser.add_argument('-o', '--output', dest='output', default='switch_fft.tsv', action='store', type=str, help='Filename to save integrated spectrum to. Default is tab delimited, switch_fft.tsv.')

    args = parser.parse_args()

    NFFT=2**8

    CTRL_PIN = 'GPIO17'
    SENS_PIN = 'GPIO27'
    noise_ctrl = gpiozero.DigitalOutputDevice(CTRL_PIN) 
    noise_sens = gpiozero.DigitalInputDevice(SENS_PIN)
    # Take a 30s noise source measurement and save to timestamped file.
    print('Enabling noise source: \nSwitching GPIO pin {} high.'.format(
        CTRL_PIN))
    # Switch CTRL pin and check Vcc to confirm
    noise_ctrl.on()
    assert 1 == noise_ctrl.value, '{} did not return HIGH, noise source not switched on.'.format(CTRL_PIN)

    freqs_noise, p_xx_avg_noise = run_integrate( NFFT, args.gain, args.rate, args.fc, 30.0 )
    save_spectrum('noise_{}.txt'.format(int(time.time())),
                  freqs_noise,
                  p_xx_avg_noise)

    print('Disabling noise source: \nSwitching GPIO pin {} low.'.format(
        CTRL_PIN))
    noise_ctrl.off()
    assert 0 == noise_ctrl.value, '{} did not return LOW, noise source not switched off.'.format(CTRL_PIN)

    freqs, p_xx_avg = run_integrate( NFFT, args.gain, args.rate, args.fc, args.t_int )
    save_spectrum(args.output,
                  freqs, 
                  p_xx_avg)

