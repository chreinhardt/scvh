#!/usr/bin/env python2
"""
Plot the EOS tables.
"""
from matplotlib import *
from matplotlib.pyplot import *
import numpy

def main():
    """
    Setup the plot.
    """

    # Set a font
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 10.0

    # Legend
    # mpl.rcParams['legend.handlelength']  = 2.9
    rcParams['legend.handlelength']  = 0.5
    rcParams['legend.frameon']       = False
    rcParams['legend.numpoints']     = 1
    rcParams['legend.scatterpoints'] = 1

    # Adjust axes line width
    rcParams['axes.linewidth']   = 0.5

    # Adjust ticks
    rcParams['xtick.major.size'] = 4
    rcParams['xtick.minor.size'] = 2
    rcParams['ytick.major.size'] = 4
    rcParams['ytick.minor.size'] = 2

    # Adjust Font Size
    rcParams['xtick.labelsize']  = 'x-small'
    rcParams['ytick.labelsize']  = 'x-small'
    rcParams['axes.labelsize']   = 'small'

    # Adjust Line Width and Marker Size
    rcParams['lines.markersize']  = 5
    rcParams['lines.linewidth']  = 0.5

    # Restore classic font used for math
    rcParams['mathtext.fontset'] = 'cm'
    rcParams['mathtext.rm']      = 'serif'

    # Set Up Figure, Single Column MNRAS
    fig = gcf()
    ax = gca()
    fig, ax = subplots(1,1)
    fig.set_size_inches(8.27*0.39,8.27*(6./8.)*0.39)

    """
    Entropy of the SCVH EOS table for H.
    """
    nRho = 201
    nT   = 100

    data = numpy.loadtxt("scvheoscalcentropy_gsl.txt")

    rho    = data[:, 0] 
    T      = data[:, 1]
    s      = data[:, 2]
    s_calc = data[:, 3]

    rho = rho[0:nRho]
    T   = T[0:numpy.size(T):nRho]

    rho_min = numpy.min(rho)
    rho_max = numpy.max(rho)
    T_min   = numpy.min(T)
    T_max   = numpy.max(T)
    
    print "rho_min=", rho_min
    print "rho_max=", rho_max
    print "T_min  =", T_min
    print "T_max  =", T_max

    print "rho=", rho
    print "T  =", T
    print

    index_T = numpy.array([19, 39, 59, 79, 99])

    # Split into arrays of constant T
    s_array = numpy.split(s, nT)
    s_calc_array = numpy.split(s_calc, nT)
    
    """
    Plot s(rho, T).
    """
    # Plot all isotherms
    for i in range(0, nT):
        loglog(rho, s_array[i], '-', color='blue')
        loglog(rho, s_calc_array[i], '--', color='green')

    title("SCVH EOS for H")
    xlabel("Log Density [g cm$^{-3}$]")
    ylabel("Log Entropy [erg g$^{-1}$ K$^{-1}$]")

    savefig('scvheoscalcentropy_gsl.png', dpi=300, bbox_inches='tight')

    fig = gcf()
    fig.clear()

    # Plot only a few isotherms
    for i in index_T:
        loglog(rho, s_array[i], '-', color='blue')
        loglog(rho, s_calc_array[i], '--', color='green')
        print T[index_T]

    title("SCVH EOS for H")
    xlabel("Log Density [g cm$^{-3}$]")
    ylabel("Log Entropy [erg g$^{-1}$ K$^{-1}$]")

    savefig('scvheoscalcentropy_gsl.png', dpi=300, bbox_inches='tight')

    show()

    exit(0)

if __name__ == '__main__':
    main()

