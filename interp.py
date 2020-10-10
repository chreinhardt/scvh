#!/usr/bin/env python2
"""
Plot the original EOS.
"""
from matplotlib import *
from matplotlib.pyplot import *
import numpy
from scipy import interpolate 

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

    print numpy.where(logrho_table == logrho_table[0])

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

    index = numpy.min(numpy.where(logT_table>=2.0))
    print "i=", index, "logT=", logT_table[index]

    # Split into arrays of constant T
    logP_array = numpy.split(logP_table, nT)
    logu_array = numpy.split(logu_table, nT)
    logs_array = numpy.split(logs_table, nT)

    # Generate an array with shape (nT, nRho)
    logP_array_new = numpy.reshape(logP_table, (-1, nRho)) 
    logu_array_new = numpy.reshape(logu_table, (-1, nRho)) 
    logs_array_new = numpy.reshape(logs_table, (-1, nRho)) 

    # Returns a function that does interpolate bilinear 
    interpP = interpolate.interp2d(logrho_table, logT_table, logP_array_new, kind='linear')    
    interpU = interpolate.interp2d(logrho_table, logT_table, logu_array_new, kind='linear')    
    interpS = interpolate.interp2d(logrho_table, logT_table, logs_array_new, kind='linear')    

    logT_int = numpy.linspace(logT_min, logT_max, nT)
    logrho_int = numpy.linspace(logrho_min, logrho_max, nRho)

    """
    Pressure.
    """
    logP_int = interpP(logrho_int, logT_int)

    # Plot all isotherms from the original EOS table
    for i in range(0, nT):
        plot(logrho_table, logP_array[i], '-')

    # Plot interpolated values
    for i in range(0, nT):
        plot(logrho_int, logP_int[i, :], '.')

    title("SCVH EOS for H")
    xlabel("Log Density")
    ylabel("Log Pressure")

    savefig('interp_pressure.png', dpi=300, bbox_inches='tight')
    
    #show()

    fig = gcf()
    fig.clear()

    """
    Internal energy.
    """
    logu_int = interpU(logrho_int, logT_int)

    # Plot all isotherms from the original EOS table
    for i in range(0, nT):
        plot(logrho_table, logu_array[i], '-')

    # Plot interpolated values
    for i in range(0, nT):
        plot(logrho_int, logu_int[i, :], '.')

    xlabel("Log Density")
    ylabel("Log Internal Energy")

    savefig('interp_intenergy.png', dpi=300, bbox_inches='tight')

    #show()
    fig = gcf()
    fig.clear()

    """
    Entropy.
    """
    logs_int = interpS(logrho_int, logT_int)

    # Plot all isotherms from the original EOS table
    for i in range(0, nT):
        plot(logrho_table, logs_array[i], '-')

    # Plot interpolated values
    for i in range(0, nT):
        plot(logrho_int, logs_int[i, :], '.')

    xlabel("Log Density")
    ylabel("Log Entropy")

    savefig('interp_entropy.png', dpi=300, bbox_inches='tight')
    
    show()
    exit(1)
    fig = gcf()
    fig.clear()

    # Plot SCVH only
    print "Plottting from index", index+2, " logT=", logT_table[index+2]
    for i in range(index+2,  nT):
        plot(logrho_table, logs_array[i], '-')

    xlabel("Log Density")
    ylabel("Log Entropy")

    savefig('ploteostable_original_h_sofrho_scvheosonly.png', dpi=300, bbox_inches='tight')

    show()
    exit(0)

if __name__ == '__main__':
    main()

