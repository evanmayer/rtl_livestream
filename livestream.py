'''
by Evan Mayer
for Steward Observatory Marrone Lab
8/30/2017

Livestream a frequency spectrum at one center frequency

Uses matplotlib animation API and pyrtlsdr library.

This is a combination of a matplotlib example code found at
http://matplotlib.org/examples/animation/random_data.html
and roger-'s waterfall plotter demo found at 
https://github.com/roger-/pyrtlsdr.
'''
from __future__ import print_function

# Import plotting necessities
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.mlab import psd

from rtlsdr import RtlSdr

# Use psd parameters and var names from roger-'s demo_waterfall.py file
NFFT = 512
NUM_SAMPLES_PER_SCAN = NFFT*16

# Choose spectrum window by altering these values.
sdr = RtlSdr()
sdr.rs = 2.4e6 # Rate of Sampling (intrinsically tied to bandwidth with SDR dongles)
sdr.fc = 89.1e6 # Frequency Center
sdr.gain = 10
print('  sample rate: %0.6f MHz' % (sdr.rs/1e6))
print('  center frequency %0.6f MHz' % (sdr.fc/1e6))
print('  gain: %d dB' % sdr.gain)

# Initial plotting
fig, ax = plt.subplots()
iq = sdr.read_samples(NUM_SAMPLES_PER_SCAN) # get initial data from sdr
p_xx, freqs, psdlines = plt.psd(iq, NFFT=NFFT, Fc = sdr.fc/1e6, Fs=sdr.rs/1e6, 
    return_line = True)
ax.set_ylim(-60, -10) # y-vals in dB/Hz

'''
update the the amplitude y-vals in the only psd line instance being plotted, [0]
inputs "data", a generator from func data_gen
returns size 1 tuple containing first element of lines object "psdlines"
'''
def update(data):
    psdlines[0].set_ydata(data)
    return psdlines[0],

'''
wraps read_samples, the pyrtlsdr data getter, to pass new y-vals into 
animation below
inputs: none
outputs: psd_y, a generator that gives the ydata of the psd plot.
'''
def data_gen():
        iq = sdr.read_samples(NUM_SAMPLES_PER_SCAN)
        p_xx, freqs, psdlines = plt.psd(iq, NFFT=NFFT, Fc = sdr.fc/1e6, Fs=sdr.rs/1e6, 
            return_line = True, animated = True)
        psd_y = psdlines[0].get_ydata()
        yield psd_y

# play animation
ani = animation.FuncAnimation(fig, update, data_gen, interval=0, blit = True)
plt.show()

# nice and tidy
sdr.close()
