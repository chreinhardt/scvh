#!/home/creinh/python_venv/bin/python3
""" Generate a rectangular eos table in P-T and the invert it to rho-T. """
from scipy import interpolate as interp
from scipy import optimize
import numpy as np
import argparse

import scvheosext

def main(): 
    #input_file = "../eos-tables-pt/hydrogen_scvh_extended.data"
    #output_file = "../h_ext_dt.data"

    # Command line parameters
    parser = argparse.ArgumentParser(description="Invert the original eos tables from P-T to rho-T.")

    # Required positional argument
    parser.add_argument('input_file', type=str, help="Input eos table (in P-T, not rectangular)")
    parser.add_argument('output_file', type=str, help="Output eos table")
    parser.add_argument('--logrho_min', type=float, help="Minimum density (default: -20.0)", default=-20.0)
    parser.add_argument('--logrho_max', type=float, help="Maximum density (default: 2.95)", default=2.95)
    parser.add_argument('--nrho', type=int, help="Number of grid points in logrho (default: 460)", default=460)

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    logrho_min_ref = -20.0
    logrho_max_ref = 2.95
    nrho_ref = 460

    logrho_min = args.logrho_min
    logrho_max = args.logrho_max
    nrho = args.nrho

    print("Invert eos table: {:}".format(input_file))
    print("nrho= {:} logrho_min= {:} logrho_max= {:}".format(nrho, logrho_min, logrho_max))

    # Load original (not rectangular) eos table
    eos_table_pt = scvheosext.read_original_table_pt(input_file)

    # Generate rectangular table extrapolating all isotherms
    eos_table_rect_pt = scvheosext.generate_rect_table_pt(eos_table_pt)

    # Invert eos table
    eos_table_rhot = scvheosext.inv_eos_table_rhot(eos_table_rect_pt, logrho_min, logrho_max, nrho)

    scvheosext.write_eos_table_rhot(output_file, eos_table_rhot)

    exit(0)

if __name__ == '__main__':
    main()
