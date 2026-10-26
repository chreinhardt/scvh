#!/usr/bin/env python3
"""
Plot the entropy of the original and the fixed eos table.
"""
from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


# Fix issues with gui
mpl.use('Agg')


def read_eos_table_dt(filename, delimiter=None):
    """ Read an EOS table that is tabulated in rho and T. """
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


def write_eos_table(output_file, logrho_axis, logT_axis, logP_array, logu_array, logs_array, delimiter=" ", comments="#", input_file=""):
    """ Write an EOS table to an output file. """
    nRho = np.size(logrho_axis)
    nT = np.size(logT_axis)

    # Store the different isotherms in 1d arrays
    logrho = np.tile(logrho_axis, nT)
    logT = np.repeat(logT_axis, nRho)

    logP = logP_array.flatten(order='C')
    logu = logu_array.flatten(order='C')
    logs = logs_array.flatten(order='C')
 
    header = " nT = {:} nRho= {:} (input file: {:})\n"\
             " {:} {:} {:} {:} {:}".format(nT, nRho, input_file, "logT [K]", "logRho [g/cc]", "logP [barye]", "logE [erg/g]", "logS [erg/g/K]")

    np.savetxt(output_file, np.column_stack([logT, logrho, logP, logu, logs]), header=header, fmt='%.8e', delimiter=delimiter, comments=comments)


def main():
    scvh_old = read_eos_table_dt("scvh_extended_dt_hydrogen_722_helium_278.data")
    scvh_new = read_eos_table_dt("scvh_extended_dt_hydrogen_722_helium_278.data_fixed_entropy")

    data = np.loadtxt("logrho_limit_pt.txt")

    logT_limit = data[:,0]
    logrho_limit_min = data[:,1]
    logrho_limit_max = data[:,2]

   # Check if the old and the new table agree
    if np.any(np.abs((scvh_old['logrho']-scvh_new['logrho'])/scvh_old['logrho']) > 1e-8):
        print("logrho axis differ.")

    if np.any(np.abs((scvh_old['logT']-scvh_new['logT'])/scvh_old['logT']) > 1e-8):
        print("logT axis differ.")

    if np.any(np.abs((scvh_old['logP']-scvh_new['logP'])/scvh_old['logP']) > 1e-8):
        print("logP axis differ.")

    if np.any(np.abs((scvh_old['logu']-scvh_new['logu'])/scvh_old['logu']) > 1e-8):
        print("logu axis differ.")

    if np.any(np.abs((scvh_old['logs']-scvh_new['logs'])/scvh_old['logs']) > 1e-8):
        print("logs axis differ.")

    # Plot isentropes
    fig, ax = plt.subplots(1, 1)

    levels = np.linspace(np.amin(scvh_old['logs']), np.amax(scvh_old['logs']), 201)
    levels = np.sort(levels)

    ax.contour(scvh_new['logrho'], scvh_new['logT'], scvh_new['logs'], levels=levels)
    ax.contour(scvh_old['logrho'], scvh_old['logT'], scvh_old['logs'], levels=levels, linestyles='dashed')

    # Show limits of the original table
    ax.plot(logrho_limit_max, logT_limit, color='red')
    ax.plot(logrho_limit_min, logT_limit, color='red')

    ax.set_xlim(scvh_old['logrho_min'], scvh_old['logrho_max'])
    ax.set_ylim(scvh_old['logT_min'], scvh_old['logT_max'])

    ax.set(xlabel="log(rho) [g cm$^{-3}$]", ylabel="log(T) [K]")

    # Show where logs is NaN
    index_nan = np.where(np.isnan(scvh_new['logs']) == True)
    
    xx, yy = np.meshgrid(scvh_new['logrho'], scvh_new['logT'])
    ax.scatter(xx[index_nan], yy[index_nan])
        
    plt.savefig("fix_entropy_dt.png", dpi=200, bbox_inches='tight')
    exit(0)

if __name__ == '__main__':
    main()

