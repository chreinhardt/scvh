#!/usr/bin/env python2
"""
Test the derivatives.
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
    Plot dP/drho(rho, T).
    """
    data = numpy.loadtxt("testscvheosderivs_dpdrhoofrhot.txt")

    print numpy.where(data <= 0.0)

    imshow(data)
    colorbar()

    xlabel("nT")
    ylabel("nRho")

    savefig('testscvheosderivs_dpdrhoofrhot.png', dpi=300, bbox_inches='tight')
    #show()

    fig = gcf()
    fig.clear()

    """
    Plot dP/dT(rho, T).
    """
    data = numpy.loadtxt("testscvheosderivs_dpdtofrhot.txt")

    print numpy.where(data <= 0.0)

    imshow(data)
    colorbar()

    xlabel("nT")
    ylabel("nRho")

    savefig('testscvheosderivs_dpdtofrhot.png', dpi=300, bbox_inches='tight')
    
    #show()
    exit(0)

if __name__ == '__main__':
    main()

