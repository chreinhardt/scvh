#!/usr/bin/env python3
from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import os

# Fix issues with gui
mpl.use('Agg')

def read_eos_table_dt(filename, delimiter=None):
    data = np.loadtxt(filename, delimiter=delimiter) 

    logT_table   = data[:, 0]
    logrho_table = data[:, 1]
    logP_table   = data[:, 2]
    logu_table   = data[:, 3]
    logs_table   = data[:, 4]

    # All the SCvH EOS tables are tabulated along isotherms
    logT = np.unique(logT_table)
    nT = np.size(logT)

    # For the mixture tables the number of grid points in rho the same for each isotherm
    nRho = int(np.size(logT_table)/nT)

    print("nT = {:} nRho = {:}".format(nT, nRho))

    logrho = logrho_table[0:nRho]

    logrho_min = np.min(logrho)
    logrho_max = np.max(logrho)
    logT_min   = np.min(logT)
    logT_max   = np.max(logT)

    print("logrho_min = {:}".format(logrho_min))
    print("logrho_max = {:}".format(logrho_max))
    print("logT_min   = {:}".format(logT_min))
    print("logT_max   = {:}".format(logT_max))
    print()

    # Split into arrays of constant T
    logP_array = np.split(logP_table, nT)
    logu_array = np.split(logu_table, nT)
    logs_array = np.split(logs_table, nT)

    # Generate 2d arrays
    logP = np.vstack(logP_array)
    logu = np.vstack(logu_array)
    logs = np.vstack(logs_array)

    eos_table_dt = {
        "nT":         nT,
        "nRho":       nRho,
        "logT":       logT,
        "logrho":     logrho,
        "logrho_min": logrho_min,
        "logrho_max": logrho_max,
        "logT_min":   logT_min,
        "logT_max":   logT_max,
        "logP":       logP,
        "logu":       logu,
        "logs":       logs,
    }
    
    return eos_table_dt

def main():

    scvh_ext_dt = read_eos_table_dt("scvh_extended_dt_hydrogen_722_helium_278.data")

    data = np.loadtxt("logrho_limit_pt.txt")

    logT_limit = data[:,0]
    logrho_limit_min = data[:,1]
    logrho_limit_max = data[:,2]

    # Plot the data
    fig, ax = plt.subplots(1, 1)

    # Isentropes
    ax.contour(scvh_ext_dt['logrho'], scvh_ext_dt['logT'], scvh_ext_dt['logs'], levels=200)

    # Limit of the original table in P-T
    ax.plot(logrho_limit_min, logT_limit, color='red')
    ax.plot(logrho_limit_max, logT_limit, color='red')

    ax.set_xlim(scvh_ext_dt['logrho_min'], scvh_ext_dt['logrho_max'])
    ax.set_ylim(scvh_ext_dt['logT_min'], scvh_ext_dt['logT_max'])

    plt.savefig("plot_isentropes_dt.png", dpi=100, bbox_inches='tight')

    exit(0)

if __name__ == '__main__':
    main()

