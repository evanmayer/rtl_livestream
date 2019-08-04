'''
by Evan Mayer

Perform a frequency-switched measurement of spectra made from RTL-SDR samples.
Based on sample code from roger-'s pyrtlsdr library.
'''

# Import plotting necessities
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.mlab import psd
import time

from rtlsdr import RtlSdr


def run_fswitch( NFFT, gain, rate, fc, fthrow, fswitch ):
    '''
    NFFT:    Pass-trough for matplotlib plot psd: powers of 2 are most efficient
    gain:    Requested SDR gain (dB)
    rate:    SDR sample rate, intrinsically tied to bandwidth in SDRs (Hz_
    fc:      base center frequency (Hz)
    fthrow:  frequency to switch to, samples subtracted from spectrum of fc
             (Hz)
    fswitch: switching rate (Hz)
    '''

    # Initial plotting
    fig, ax = plt.subplots()
    plt.grid()
    ax.set_ylim(-1e-4, 1e-4)

    sdr = RtlSdr()
    sdr.rs = rate # Rate of Sampling (intrinsically tied to bandwidth with SDR dongles)
    sdr.fc = fc
    sdr.gain = gain
    print('  sample rate: %0.6f MHz' % (sdr.rs/1e6))
    print('  center frequency %0.6f MHz' % (sdr.fc/1e6))
    print('  gain: %d dB' % sdr.gain)

    # The amount of time to integrate on each frequency
    delta_t = 1./(2.*fswitch)

    iq = sdr.read_samples(NFFT) # get initial data from sdr
    p_xx, freqs, psdlines = plt.psd(iq, NFFT=NFFT, Fc=sdr.fc/1e6, Fs=sdr.rs/1e6, return_line=True, animated=True)
    line = psdlines[0]


    '''
    update the the amplitude y-vals in the only psd line instance being plotted, [0]
    inputs "data", a generator from func data_gen
    returns size 1 tuple containing the updated line object
    '''
    def update(data):
        line.set_ydata(data)
        return line,


    '''
    wraps read_samples, the pyrtlsdr data getter, to pass new y-vals into
    animation below

    Integrates frequency-switched data.

    inputs: sdr, the pyrtlsdr object
    outputs: psd_y, a generator that gives the ydata of the psd plot.
    '''
    def data_gen():
        # Ensure set to base frequency at start of every cycle
        sdr.fc = fc

        # Set the time: the start of the integration cycles
        start_time = time.time()

        p_xx_on_tot = np.zeros(NFFT)
        # Integrate on the target frequency
        iq_on = []
        while time.time()-start_time < delta_t:
            # Collect on-frequency samples
            iq_on.extend(sdr.read_samples(NFFT))
        # Take the PSD of the on samples
        p_xx_on, freqs_on = psd(iq_on, NFFT=NFFT, Fs=sdr.rs/1e6)
        p_xx_on_tot += p_xx_on

        # Switch to the switched frequency
        sdr.fc = fc + fthrow

        start_time = time.time()

        p_xx_off_tot = np.zeros(NFFT)
        # Integrate on the switched frequency
        iq_off = []
        while time.time()-start_time < delta_t:
            # Collect off-frequency samples
            iq_off.extend(sdr.read_samples(NFFT))
        # And the off samples
        p_xx_off, freqs_off = psd(iq_off, NFFT=NFFT, Fs=sdr.rs/1e6)
        p_xx_off_tot += p_xx_off

        if not p_xx_on_tot.any() or not p_xx_off_tot.any():
            print('Warning: One of the frequency-switched spectra was empty. Consider lowering the switching frequency in order for RTL-SDR sampling to catch up.')
        # Add folded samples to the current and update
        p_xx_tot = (p_xx_on_tot - p_xx_off_tot)/2.

        yield p_xx_tot


    # play animation
    ani = animation.FuncAnimation(fig, update, data_gen, interval=0, blit=True)
    plt.show()

    # nice and tidy
    sdr.close()

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gain', dest='gain', default=30.0, action='store', type=float, help='Requested SDR gain (dB)')
    parser.add_argument('-r', '--rate', dest='rate', default=2.56e6, action='store', type=float, help='Sample rate (Hz)')
    parser.add_argument('-c', '--fc', dest='fc', default=1.4202e9, action='store', type=float, help='Requested base frequency (Hz)')
    parser.add_argument('-t', '--fthrow', dest='fthrow', default=0.0004e9, action='store', type=float, help='Delta to fc applied on every switching integration cycle. (Hz)')
    parser.add_argument('-s', '--fswitch', dest='fswitch', default=0.2, action='store', type=float, help='Frequency at which SDR tuning frequency is switched from fc to fc+fthrow. One cycle is from fc to fc+fthrow bck to fc. (Hz)')
    args = parser.parse_args()

    NFFT=1024

    run_fswitch( NFFT, args.gain, args.rate, args.fc, args.fthrow, args.fswitch )
