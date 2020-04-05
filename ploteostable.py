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
    SCVH EOS table for H.
    """
    nRho = 201
    nT   = 100

    """
    Pressure.
    """
    data = numpy.loadtxt("testscvheos_h_p_rho.txt")

    logrho = data[:,0]
    logP   = data[:,1:nT+1]

    logrho_min = numpy.min(logrho)
    logrho_max = numpy.max(logrho)

    print "logrho_min=", logrho_min
    print "logrho_max=", logrho_max
    print "logrho=", logrho
    print

    # Plot all isotherms
    for i in range(0, nT):
        plot(logrho, logP[:,i], '-')

    """
    # Now plot the isotherms below logT=2.0
    for i in range(0, index+1):
        plot(logrho_table, logP_array[i], '-', linewidth=1)
    
    for i in range(index+1, index+20):
        plot(logrho_table, logP_array[i], '--', linewidth=1)
    """

    # Zoom in to the weird region
    #xlim(-3.5, -3.0)
    #ylim(5.8, 7.0) 

    xlabel("Log Density")
    ylabel("Log Pressure")

    savefig('ploteostable_h_pofrho.png', dpi=300, bbox_inches='tight')

    #show()
    fig = gcf()
    fig.clear()

    """
    Internal energy.
    """
    data = numpy.loadtxt("testscvheos_h_u_rho.txt")

    logrho = data[:,0]
    logu   = data[:,1:nT+1]

    logrho_min = numpy.min(logrho)
    logrho_max = numpy.max(logrho)

    print "logrho_min=", logrho_min
    print "logrho_max=", logrho_max
    print "logrho=", logrho
    print

    # Plot all isotherms
    for i in range(0, nT):
        plot(logrho, logu[:,i], '-')

    """
    # Now plot the isotherms below logT=2.0
    for i in range(0, index+1):
        plot(logrho_table, logP_array[i], '-', linewidth=1)
    
    for i in range(index+1, index+20):
        plot(logrho_table, logP_array[i], '--', linewidth=1)
    """

    # Zoom in to the weird region
    #xlim(-3.5, -3.0)
    #ylim(5.8, 7.0) 

    xlabel("Log Density")
    ylabel("Log Internal Energy")

    savefig('ploteostable_h_uofrho.png', dpi=300, bbox_inches='tight')

    #show()
    fig = gcf()
    fig.clear()

    """
    Entropy.
    """
    data = numpy.loadtxt("testscvheos_h_s_rho.txt")

    logrho = data[:,0]
    logs   = data[:,1:nT+1]

    logrho_min = numpy.min(logrho)
    logrho_max = numpy.max(logrho)

    print "logrho_min=", logrho_min
    print "logrho_max=", logrho_max
    print "logrho=", logrho
    print

    # Plot all isotherms
    for i in range(0, nT):
        plot(logrho, logs[:,i], '-')

    """
    # Now plot the isotherms below logT=2.0
    for i in range(0, index+1):
        plot(logrho_table, logP_array[i], '-', linewidth=1)
    
    for i in range(index+1, index+20):
        plot(logrho_table, logP_array[i], '--', linewidth=1)
    """

    # Zoom in to the weird region
    #xlim(-3.5, -3.0)
    #ylim(5.8, 7.0) 

    xlabel("Log Density")
    ylabel("Log Entropy")

    savefig('ploteostable_h_sofrho.png', dpi=300, bbox_inches='tight')

    show()
    fig = gcf()
    fig.clear()

    exit(0)

if __name__ == '__main__':
    main()

