# rtl_livestream
Shows an animated 1D plot of the PSD from an RTL-SDR dongle, like slices of a waterfall plot.
Uses matplotlib's FuncAnimation library to update the plot.

Depends on:
numpy
matplotlib
pyrtlsdr by roger-: https://github.com/roger-/pyrtlsdr

Currently dumb; you must edit the python source to change the SDR's window size and center frequency.
I may do the right thing later and add command line args. For now, initialize with:
$ python livestream.py

This is a derivative work (8/30/2017) of roger-'s pyrtlsdr code demo_waterfall.py, which is provided under a GPL v3.0 license.
The same will apply to this.
https://www.gnu.org/licenses/gpl-3.0.en.html
