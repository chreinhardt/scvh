#!/usr/bin/env python3
"""
The entropy in the region, where the EOS was extended to low temperatures (logT < 2.1) is
inconsistent with the 1. law of thermodynamics. This mean that if T(rho, s=const.) is calculated
from

du = P/rho^2*drho + T*ds

and using a root finder

s(rho=const, T) - s = 0

the curves are inconsistent above logT = 2.1 due to a kink in the entropy. In order to obtain more
accurate (but still inconsistent) values the entropy is calculated from a diatomic ideal gas.

Note that the kink in the entropy is not due to a different reference entropy s0 but it seems that
each isochore (rho=const.) has a different reference value so simply shifting the entropy value so
that they match at logT = 2.1 does not work.
"""
from __future__ import print_function
import numpy as np


def entropy_ideal_gas(rho, T, rho0, T0, s0, mu=2.3, gamma=5.0/3.0):
    # s0 = s(rho0, T0) has to be consistent
    kB = 1.38e-16  # erg/K
    m_H = 1.67e-24 # g
    return kB/(mu*m_H)*(1.0/(gamma-1.0)*np.log(T/T0) - np.log(rho/rho0)) + s0


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

    scvh_ext_dt = read_eos_table_dt("scvh_extended_dt_hydrogen_722_helium_278.data")

    # At logT = 2.1 the SCvH EOS table starts
    index_start_scvh = np.min(np.where(scvh_ext_dt['logT'] >= 2.1)[0])

    print("The SCvH EOS starts from index: {:}".format(index_start_scvh))

    # Extend the entropy to low T with an ideal gas
    logs_shifted = np.copy(scvh_ext_dt['logs'])

    verbose = False
 
    for index_logrho0 in range(scvh_ext_dt['nRho']):
        # Calculate reference point s0 for each isochore
        index_logT0 = index_start_scvh

        logrho0 = scvh_ext_dt['logrho'][index_logrho0]
        logT0 = scvh_ext_dt['logT'][index_logT0]

        logs0 = scvh_ext_dt['logs'][index_logT0, index_logrho0]

        if verbose is True:
            print("Reference point: logs0 = {:} logrho0 = {:} (iRho= {:}) logT0 = {:} (iT= {:})".format(logs0, logrho0, index_logrho0, logT0, index_logT0))

        T0 = 10.0**logT0
        rho0 = 10**logrho0
        s0 = 10.0**logs0

        T = 10.0**scvh_ext_dt['logT'][0:index_logT0]
        rho = 10.0**scvh_ext_dt['logrho'][index_logrho0]
    
        logs_of_rho = np.log10(entropy_ideal_gas(rho, T, rho0, T0, s0, mu=2.3, gamma=1.4))
        
        if np.any(np.isnan(logs_of_rho)):
            print("Invalid values found, skip logrho= {:}".format(np.log10(rho)))
        else:
            logs_shifted[0:index_logT0,index_logrho0] = np.log10(entropy_ideal_gas(rho, T, rho0, T0, s0, mu=2.3, gamma=1.4))


    # Write corrected eos table
    write_eos_table("scvh_extended_dt_hydrogen_722_helium_278.data_fixed_entropy", scvh_ext_dt['logrho'], scvh_ext_dt['logT'], scvh_ext_dt['logP'], scvh_ext_dt['logu'], logs_shifted, input_file="scvh_extended_dt_hydrogen_722_helium_278.data")
   
   # Check if the old and the new table agree
    scvh_old = read_eos_table_dt("scvh_extended_dt_hydrogen_722_helium_278.data")
    scvh_new = read_eos_table_dt("scvh_extended_dt_hydrogen_722_helium_278.data_fixed_entropy")

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

    exit(0)

if __name__ == '__main__':
    main()

