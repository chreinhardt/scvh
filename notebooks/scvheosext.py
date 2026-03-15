#!/home/creinh/python_venv/python3
""" 
This file provides functions to read, write and modify the modified tables of the SCvH EOS
extended to lower pressures and temperatures.
"""
#import matplotlib.pyplot as plt
#import matplotlib as mpl
from scipy import interpolate as interp
from scipy import optimize
import numpy as np

def read_original_table_pt(input_file):
    """ Read the original eos tables in P-T where each isotherm has a different number of grid points. """
    data = np.loadtxt(input_file) 

    logT_table   = data[:, 0] 
    logP_table   = data[:, 1]
    frac_H2      = data[:, 2]
    frac_H       = data[:, 3]
    logrho_table = data[:, 4]
    logu_table   = data[:, 5]
    logs_table   = data[:, 6]

    # All the SCvH EOS tables are tabulated along isotherms
    logT_table_axis = np.unique(logT_table)
    nT = np.size(logT_table_axis)

    print("Number of isotherms: nT = {:}".format(nT))

    # The number of grid points in P are different for each isotherm
    logP_table_axis = list()
    logrho_isotherm = list()
    logu_isotherm = list()
    logs_isotherm = list()

    for logT in logT_table_axis:
        logP_table_axis.append(logP_table[np.where(logT_table == logT)])
        logrho_isotherm.append(logrho_table[np.where(logT_table == logT)])
        logu_isotherm.append(logu_table[np.where(logT_table == logT)])
        logs_isotherm.append(logs_table[np.where(logT_table == logT)])

    logT_min = np.min(logT_table_axis)
    logT_max = np.max(logT_table_axis)

    print("logT_min = {:}".format(logT_min))
    print("logT_max = {:}".format(logT_max))
    print()

    # Determine the number of grid points in P for each isotherm
    nP = list()
    for logP in logP_table_axis:
        nP.append(np.size(logP))

    # Check if the distance between each data point in logP is constant
    dlogP_global = []
    for indexT in range(nT):
        # dLogP from each grid point
        dlogP_isotherm = logP_table_axis[indexT][1:]- logP_table_axis[indexT][:-1]

        nP_isotherm =nP[indexT]
        dlogP_global.append((logP_table_axis[indexT][-1] - logP_table_axis[indexT][0])/(nP_isotherm-1))

        
        diff = np.abs((dlogP_isotherm - dlogP_global[indexT])/dlogP_isotherm)

        if np.size(np.where(diff > 1e-10)) != 0:
            print("isotherm {:}: dlogP not constant for all grid points".format(indexT))
            return None

    # Check that dlogP is the same for all isotherms
    dlogP = np.array(dlogP_global)
    diff = np.abs((dlogP[1:]-dlogP[:-1])/dlogP[1:])
    print(diff)
    
    if np.size(np.where(diff > 1e-10)) != 0:
        print("dlogP not constant for all isotherms")
        return None    
    else:
        dlogP = np.min(dlogP)

    eos_table = {
        'logT': logT_table_axis,
        'logP': logP_table_axis,
        'logrho': logrho_isotherm,
        'logu': logu_isotherm,
        'logs': logs_isotherm,
        'nT': nT,
        'nP': nP,
        'dlogP': dlogP,
    }

    return eos_table

def read_eos_table_rhot(filename, delimiter=None):
    """ Read a rectangular EOS table in rho and T """
    data = np.loadtxt(filename, delimiter=delimiter) 

    logT_table   = data[:, 0]
    logrho_table = data[:, 1]
    logP_table   = data[:, 2]
    logu_table   = data[:, 3]
    logs_table   = data[:, 4]

    # All the SCvH EOS tables are tabulated along isotherms
    logT = np.unique(logT_table)
    nT = np.size(logT)

    num_entry = np.size(logT_table)

    # Check if the table is rectangular
    if num_entry % nT != 0:
        print("Table not rectangular.")
        return None

    # For the mixture tables the number of grid points in rho the same for each isotherm
    nRho = int(num_entry/nT)
    
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


def read_eos_table_pt(filename, delimiter=None):
    """ Read a rectangular EOS table in P and T """
    data = np.loadtxt(filename, delimiter=delimiter) 

    logT_table   = data[:, 0] 
    logP_table   = data[:, 1]
    logrho_table = data[:, 2]
    logu_table   = data[:, 3]
    logs_table   = data[:, 4]

    # All the SCvH EOS tables are tabulated along isotherms
    logT = np.unique(logT_table)
    nT = np.size(logT)

    num_entry = np.size(logT_table)

    # Check if the table is rectangular
    if num_entry % nT != 0:
        print("Table not rectangular.")
        print("num_entry= {:} nT= {:}".format(num_entry, nT))
        print(logT)
        return None

    # For the mixture tables the number of grid points in rho the same for each isotherm
    nP = int(num_entry/nT)
    
    print("nT = {:} nP = {:}".format(nT, nP))
    
    logP = logP_table[0:nP]
    
    logP_min = np.min(logP)
    logP_max = np.max(logP)
    logT_min   = np.min(logT)
    logT_max   = np.max(logT)

    print("logP_min = {:}".format(logP_min))
    print("logP_max = {:}".format(logP_max))
    print("logT_min   = {:}".format(logT_min))
    print("logT_max   = {:}".format(logT_max))
    print()

    # Split into arrays of constant T
    logrho_array = np.split(logrho_table, nT)
    logu_array = np.split(logu_table, nT)
    logs_array = np.split(logs_table, nT)

    # Generate 2d arrays
    logrho = np.vstack(logrho_array)
    logu = np.vstack(logu_array)
    logs = np.vstack(logs_array)
      
    eos_table_pt = {
        "nT":       nT,
        "nP":       nP,
        "logT":     logT,
        "logP":     logP,
        "logP_min": logP_min,
        "logP_max": logP_max,
        "logT_min": logT_min,
        "logT_max": logT_max,
        "logrho":   logrho,
        "logu":     logu,
        "logs":     logs,
    }
    
    return eos_table_pt


def write_eos_table_rhot(filename, eos_table, delimiter=" ", comments="#"):
    """ Write a rectangular EOS table to an output file """
    logrho = np.tile(eos_table['logrho'], eos_table['nT'])
    logT = np.repeat(eos_table['logT'], eos_table['nRho'])

    logP = eos_table['logP'].flatten(order='C')
    logu = eos_table['logu'].flatten(order='C')
    logs = eos_table['logs'].flatten(order='C')   

    header = "nT = {:} nRho= {:}\n"\
             "{:},{:},{:},{:},{:}".format(eos_table['nT'], eos_table['nRho'], "logT [K]", "logRho [g/cc]", "logP [barye]", "logE [erg/g]", "logS [erg/g/K]")
    
    np.savetxt(filename, np.column_stack([logT, logrho, logP, logu, logs]), header=header, fmt='%15.8e', delimiter=delimiter, comments=comments)


def generate_rect_table_pt(eos_table):
    """ Generate a new table that has the same number of grid points in logP for all isotherms. """
    nP = np.max(eos_table['nP'])
    print("nP = {:}".format(nP))

    # Determine the min and max of logP for all isotherms
    logP_min = np.min(np.concatenate(([logP for logP in eos_table['logP']])))
    logP_max = np.max(np.concatenate(([logP for logP in eos_table['logP']])))

    print("logP_min= {:} logP_max= {:}".format(logP_min, logP_max))

    logP_int = logP_min + (logP_max-logP_min)/(nP-1)*np.arange(0, nP)

    logrho_array = []
    logu_array = []
    logs_array = []

    # For each isotherm generate an array with nP grid points in logP
    for indexT in range(eos_table['nT']):
        logrho_interp = interp.interp1d(eos_table['logP'][indexT], eos_table['logrho'][indexT], kind='linear', fill_value='extrapolate')
        logu_interp = interp.interp1d(eos_table['logP'][indexT], eos_table['logu'][indexT], kind='linear', fill_value='extrapolate')
        logs_interp = interp.interp1d(eos_table['logP'][indexT], eos_table['logs'][indexT], kind='linear', fill_value='extrapolate')

        logrho_array.append(logrho_interp(logP_int))
        logu_array.append(logu_interp(logP_int))
        logs_array.append(logs_interp(logP_int))

    eos_table_rect = {
        'nT': eos_table['nT'],
        'nP': nP,
        'logT': eos_table['logT'],
        'logP': logP_int,
        'logrho': np.array(logrho_array),
        'logu': np.array(logu_array),
        'logs': np.array(logs_array),
    }

    return eos_table_rect


def func_logrho_pt(logP, *my_args):
    """ Calculate logrho(logT, logP) - logrho for the root finder. """
    logT, logrho, logrho_tp  = my_args
    return logrho_tp((logT, logP))-logrho


def inv_eos_table_rhot(eos_table_rect, logrho_min, logrho_max, nRho, logP_min=-1e2, logP_max=1e2, rtol=None):
    """ Invert a rectangular eos table from rho(P, T) to P(rho, T) """
    logrho_axis = np.linspace(logrho_min, logrho_max, nRho)

    # Initialize interpolators (extrapolate if input is out of bounds)
    logrho_tp = interp.RegularGridInterpolator((eos_table_rect['logT'], eos_table_rect['logP']), eos_table_rect['logrho'], method='linear', fill_value=None, bounds_error=False)
    logu_tp = interp.RegularGridInterpolator((eos_table_rect['logT'], eos_table_rect['logP']), eos_table_rect['logu'], method='linear', fill_value=None, bounds_error=False)
    logs_tp = interp.RegularGridInterpolator((eos_table_rect['logT'], eos_table_rect['logP']), eos_table_rect['logs'], method='linear', fill_value=None, bounds_error=False)

    logP_array = []
    logu_array = []
    logs_array = []
  
    for index_T in range(eos_table_rect['nT']):
        # This is shitty but it seems that root_scalar does not work with an array of inputs
        logP_rhot = np.zeros(nRho)
        logu_rhot = np.zeros(nRho)
        logs_rhot = np.zeros(nRho)

        for index_rho in range(nRho):
            logrho = logrho_axis[index_rho]
            logT = eos_table_rect['logT'][index_T]
            
            # Calculate logP(logrho, logT)
            try:
                sol = optimize.root_scalar(func_logrho_pt, args=(logT, logrho, logrho_tp), bracket=[logP_min, logP_max], method='brentq', rtol=rtol)
            except ValueError as e:
                print("Inversion failed ({:}): logT= {:} logrho= {:}".format(e, logT, logrho))
                return None

            #sol = optimize.root_scalar(fun, args=(logT, logrho, logrho_tp), bracket=[logP_min, logP_max], method='brentq')
            logP = sol.root
    
            logP_rhot[index_rho] = logP
            logu_rhot[index_rho] = logu_tp((logT, logP))
            logs_rhot[index_rho] = logs_tp((logT, logP))

        logP_array.append(logP_rhot)
        logu_array.append(logu_rhot)
        logs_array.append(logs_rhot)

    # Return the inverted EOS table
    eos_table = {
        'nT': eos_table_rect['nT'],
        'nRho': nRho,
        'logT': eos_table_rect['logT'],
        'logrho': logrho_axis,
        'logP': np.array(logP_array),
        'logu': np.array(logu_array),
        'logs': np.array(logs_array),  
    }

    return eos_table

