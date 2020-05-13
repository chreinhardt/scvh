#!/usr/bin/env python2
"""
Plot s(rho, T) for different isotherms.
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
    rcParams['axes.titlesize']   = 'small'

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
    SCVH EOS table for H.
    """
    data = numpy.loadtxt("scvh_h_dt_cgs.csv", delimiter=",", skiprows=1) 
    #data = numpy.loadtxt("scvh_he_dt_cgs.csv", delimiter=",", skiprows=1)
    #data = numpy.loadtxt("scvh_hhe_y0.275_dt_cgs.csv", delimiter=",", skiprows=1)

    logT_table    = data[:, 0] 
    logrho_table  = data[:, 1]
    logP_table    = data[:, 2]
    logu_table    = data[:, 3]
    logs_table    = data[:, 4]

    # All SCVH EOS tables obtained from Ravit have the same size
    nRho = 201
    nT   = 100

    logrho_table = logrho_table[0:nRho]
    logT_table = logT_table[0:numpy.size(logT_table):nRho]

    logrho_min = numpy.min(logrho_table)
    logrho_max = numpy.max(logrho_table)
    logT_min   = numpy.min(logT_table)
    logT_max   = numpy.max(logT_table)
    
    dlogrho = logrho_table[1:]-logrho_table[:-1]
    dlogT = logT_table[1:]-logT_table[:-1]

    print "logrho_min=", logrho_min
    print "logrho_max=", logrho_max
    print "logT_min  =", logT_min
    print "logT_max  =", logT_max

    print "logrho_table=", logrho_table
    print "logT_table  =", logT_table
    print

    index = numpy.min(numpy.where(logT_table>=2.1))
    print "i=", index, "logT=", logT_table[index]

    # Split into arrays of constant T
    logP_array = numpy.split(logP_table, nT)
    logu_array = numpy.split(logu_table, nT)
    logs_array = numpy.split(logs_table, nT)
    
    """
    Entropy.
    """
    data = numpy.loadtxt("testscvheosentropy.txt")
    
    num_isotherms = numpy.size(data[0,1:])
    
    logrho = data[:,0]
    logs   = data[:,1:num_isotherms+1]

    print "num_isotherms=", num_isotherms

    # Plot all isotherms
    for i in range(0, num_isotherms):
        if i<nT:
            plot(logrho_table, logs_array[i], '-')
        plot(logrho, logs[:, i], '--o')

    title("SCVH EOS for H")
    xlabel("Log Density")
    ylabel("Log Entropy")

    savefig('testscvheosentropy.png', dpi=300, bbox_inches='tight')
    
    show()

    fig = gcf()
    fig.clear()

    # Now plot the isotherms below logT=2.0
    for i in range(0, index+1):
        plot(logrho_table, logs_array[i], '-')
        plot(logrho, logs[:, i], '-o')

    for i in range(index+1, index+20):
        plot(logrho_table, logs_array[i], '--')
        plot(logrho, logs[:, i], '--*')

    # Zoom in to the weird region
    xlim(-1.72, -1.62)
    ylim(8.5, 8.6) 

    title("SCVH EOS for H")
    xlabel("Log Density")
    ylabel("Log Entropy")

    savefig('testscvheosentropy_zoom.png', dpi=300, bbox_inches='tight')

    show()
    exit(0)

if __name__ == '__main__':
    main()

