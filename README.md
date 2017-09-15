# rtl_livestream
Shows a 1D plot of the PSD from an RTL-SDR dongle, like slices of a waterfall plot.
Uses matplotlib's FuncAnimation library to update the plot.

Currently dumb; you must edit the python source to change the SDR's window size and center frequency.
I may do the right thing later and add command line args. For now, initialize with:
$ python livestream.py
