#!/usr/bin/env python2
"""
Test if the EOS table is properly read.
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
    data = numpy.loadtxt("scvh_h_dt_cgs.csv", delimiter=",", skiprows=1) 
    #data = numpy.loadtxt("scvh_he_dt_cgs.csv", delimiter=",", skiprows=1)
    #data = numpy.loadtxt("scvh_hhe_y0.275_dt_cgs.csv", delimiter=",", skiprows=1)

    logT_table    = data[:, 0] 
    logrho_table  = data[:, 1]
    logP_table    = data[:, 2]
    logu_table    = data[:, 3]
    logs_table    = data[:, 4]

    # All EOS tables obtained from Ravit have the same size
    nRho = 201
    nT   = 100

    print "nRho=", nRho
    print "nT=", nT

    logrho_table = logrho_table[0:nRho]
    logT_table = logT_table[0:numpy.size(logT_table):nRho]

    logrho_min = numpy.min(logrho_table)
    logrho_max = numpy.max(logrho_table)
    logT_min   = numpy.min(logT_table)
    logT_max   = numpy.max(logT_table)

    print "logrho_min=", logrho_min
    print "logrho_max=", logrho_max
    print "logT_min  =", logT_min
    print "logT_max  =", logT_max

    print "logrho_table=", logrho_table
    print "logT_table  =", logT_table
    print

    # Split into arrays of constant T
    logP_array = numpy.split(logP_table, nT)
    logu_array = numpy.split(logu_table, nT)
    logs_array = numpy.split(logs_table, nT)
    
    eps = 1e-6

    """
    Pressure.
    """
    # Read P(rho) for different T 
    data = numpy.loadtxt("testscvheos_h_p_rho.txt")

    logrho = data[:,0]
    logP   = data[:,1:nT+1]

    # Verify that the two tables agree
    index = numpy.where(numpy.abs(logP.transpose()-logP_array) > eps)

    if (numpy.size(index) > 0):
        print "Pressure differs:", index
        print "eps=", eps
        print "Check failed."
        exit(1)

    # Plot all isotherms
    for i in range(0, nT):
        plot(logrho_table, logP_array[i], '-')
        plot(logrho, logP[:,i], '--')

    xlabel("Log Density")
    ylabel("Log Pressure")

    savefig('testscvheosreadtable_h_pofrho.png', dpi=300, bbox_inches='tight')

    # Plot the difference in pressure between the two tables 
    imshow(numpy.abs(logP.transpose()-logP_array))
    colorbar()

    title("Pressure")
    xlabel("Density")
    ylabel("Temperature")

    savefig('testscvheosreadtable_h_pressure.png', dpi=300, bbox_inches='tight')

    fig = gcf()
    fig.clear()

    """
    Internal energy.
    """
    # Read u(rho) for different T 
    data = numpy.loadtxt("testscvheos_h_u_rho.txt")

    logrho = data[:,0]
    logu   = data[:,1:nT+1]
    
    # Verify that the two tables agree
    index = numpy.where(numpy.abs(logu.transpose()-logu_array) > eps)

    if (numpy.size(index) > 0):
        print "Internal energy differs:", index
        print "eps=", eps
        print "Check failed."
        exit(1)

    # Plot all isotherms
    for i in range(0, nT):
        plot(logrho_table, logu_array[i], '-')
        plot(logrho, logu[:,i], '--')

    xlabel("Log Density")
    ylabel("Log Internal Energy")

    savefig('testscvheosreadtable_uofrho.png', dpi=300, bbox_inches='tight')

    fig = gcf()
    fig.clear()

    # Plot the difference in internal energy between the two tables 
    imshow(numpy.abs(logu.transpose()-logu_array))
    colorbar()

    title("Internal energy")
    xlabel("Density")
    ylabel("Temperature")

    savefig('testscvheosreadtable_h_intenergy.png', dpi=300, bbox_inches='tight')

    fig = gcf()
    fig.clear()

    """
    Entropy.
    """
    # Read s(rho) for different T 
    data = numpy.loadtxt("testscvheos_h_s_rho.txt")

    logrho = data[:,0]
    logs   = data[:,1:nT+1]

    # Verify that the two tables agree
    index = numpy.where(numpy.abs(logs.transpose()-logs_array) > eps)

    if (numpy.size(index) > 0):
        print "Entropy differs:", index
        print "eps=", eps
        print "Check failed."
        exit(1)

    # Plot all isotherms
    for i in range(0, nT):
        plot(logrho_table, logs_array[i], '-')
        plot(logrho, logs[:,i], '--')

    xlabel("Log Density")
    ylabel("Log Entropy")

    savefig('testscvheosreadtable_sofrho.png', dpi=300, bbox_inches='tight')

    fig = gcf()
    fig.clear()

    # Plot the difference in internal energy between the two tables 
    imshow(numpy.abs(logs.transpose()-logs_array))
    colorbar()

    title("Entropy")
    xlabel("Density")
    ylabel("Temperature")

    savefig('testscvheosreadtable_h_entropy.png', dpi=300, bbox_inches='tight')

    exit(0)

if __name__ == '__main__':
    main()

