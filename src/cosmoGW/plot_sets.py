"""
plot_sets.py is a Python routine that includes some of the settings used in the
generation of plots for the cosmoGW project (see https://github.com/AlbertoRoper/cosmoGW).

Author: Alberto Roper Pol
Created: 01/01/2021
Updated: 01/11/2023 (release of the cosmoGW code)
"""

import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=20)
plt.rcParams.update({'axes.labelsize': 'x-large',
                     'axes.titlesize': 'x-large',
                     'xtick.labelsize': 'x-large',
                     'ytick.labelsize': 'x-large',
                     'legend.fontsize': 'x-large'})

def axes_lines(ax=[], both=True):

    """
    Function that includes some default settings for the lines in the loglog
    plots.
    """

    if ax == []: ax = plt.gca()
    ax.tick_params(axis='y', direction='in', length=12)
    ax.tick_params(axis='x', direction='in', length=12, pad=10)
    ax.tick_params(which='minor', direction='in', length=6)
    if both:
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
