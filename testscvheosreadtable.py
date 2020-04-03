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
    
    logT_table    = data[:, 0] 
    logrho_table  = data[:, 1]
    logP_table    = data[:, 2]
    logu_table    = data[:, 3]
    logu_table    = data[:, 4]

    print numpy.where(logrho_table == logrho_table[0])
    print numpy.where(logT_table == logT_table[0])

    exit(1)
    nRho = 201
    nT   = 42
    
    rho_min = numpy.min(rho_table)
    rho_max = numpy.max(rho_table)
    T_min   = numpy.min(T_table)
    T_max   = numpy.max(T_table)

    print "rho_min=", rho_min
    print "rho_max=", rho_max
    print "T_min  =", T_min
    print "T_max  =", T_max

    rho_table = rho_table[0:nRho]
    T_table = T_table[0:numpy.size(T_table):nRho]

    print "rho_table=", rho_table
    print "T_table  =", T_table
    print
    
    # Split into arrays of constant T
    P_array = numpy.split(P_table, nT)
    u_array = numpy.split(u_table, nT)
    
    """
    Pressure.
    """
    # Read P(rho) for different T 
    data = numpy.loadtxt("testreos3_h_p_rho.txt")

    rho = data[:,0]
    P   = data[:,1:nT+1]

    # Plot all isotherms
    for i in range(0, nT):
        loglog(rho_table, P_array[i], '-', linewidth=1)
        loglog(rho, P[:,i], '--', linewidth=1)

    xlabel("Density [g cm$^{-3}$]")
    ylabel("Pressure [GPa]")

    savefig('testreos3readtable_pofrho.png', dpi=300, bbox_inches='tight')

    show()
    fig = gcf()
    fig.clear()
    exit(1)


    """
    Internal energy.
    """
    # Read u(rho) for different T 
    data = numpy.loadtxt("testreos3_h_u_rho.txt")

    rho = data[:,0]
    u   = data[:,1:nT+1]

    # Plot all isotherms
    for i in range(0, nT):
        semilogx(rho_table, u_array[i], '-', linewidth=1)
        semilogx(rho, u[:,i], '--', linewidth=1)

    xlabel("Density [g cm$^{-3}$]")
    ylabel("Internal energy [kJ g${-1}$]")

    savefig('testreos3readtable_uofrho.png', dpi=300, bbox_inches='tight')

    show()
    fig = gcf()
    fig.clear()
    exit(1)

    # Read u(T) for different rho
    data = numpy.loadtxt("testreos3_h_u_T.txt")

    T = data[:,0]
    u = data[:,1:nRho+1]

    # Plot all curves of constant rho
    for i in range(0, nRho):
        semilogx(T_table, u_array[i], '-', linewidth=1)
        semilogx(T, u[:,i], '--', linewidth=1)

    xlabel("Temperature [K]")
    ylabel("Internal energy [kJ g${-1}$]")

    savefig('testreos3u_uofT.png', dpi=300, bbox_inches='tight')

    show()
    exit(0)

if __name__ == '__main__':
    main()

