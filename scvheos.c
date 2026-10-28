/*
 * Copyright (c) 2020 Christian Reinhardt
 *
 * This file provides an interface for the SCvH equation of state (Saumon et al. 1995) for Hydrogen
 * and Helium that was extended to lower temperatures by Vazan et al. (2013).
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 26.10.2026 (1.3.0)
 */

/* getline() is POSIX.1-2008; declare it under -std=c99 too. */
#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "scvheos.h"

/*
 * 10^x. The tables are base-10 logarithms (log10 is hard-wired in every
 * forward conversion), so the exponential is fixed too; exp(x ln 10) is
 * cheaper than pow(10, x) and accurate to an ulp or two.
 */
static double scvheosExp10(double x) {
    const double LN10 = 2.302585092994045684;   /* M_LN10 is not ISO C */
    return exp(x*LN10);
}

/*
 * Lookup failure handling, the EOSlib convention (see scvheos.h).
 *
 * A failed lookup never exits. The function counts the failure in
 * Mat->nLookupFail, prints a message if the library was compiled with
 * -DEOSLIB_VERBOSE (SCVHEOS_WARN) and returns SCVHEOS_FAIL_VALUE (in log
 * space from the log-space functions), or falls back to a table edge where
 * the contract says so (T(rho, u) -> T_min, T(rho, s) -> T_min/T_max).
 * SCVHEOS_FAIL_RETURN(Mat, ...) does the three steps in one line.
 *
 * The linear-space wrappers must test the log value BEFORE exponentiating,
 * because pow(10, -1e50) == 0 is a perfectly valid-looking number and the
 * failure would be lost. scvheosFailed(x) (scvheos.h) is true for the
 * sentinel and for NaN, in log or linear space. Initialisation errors
 * (missing or defective table) are fatal.
 */
#define SCVHEOS_FAIL_RETURN(Mat, ...) do { (Mat)->nLookupFail++; SCVHEOS_WARN(__VA_ARGS__); return SCVHEOS_FAIL_VALUE; } while (0)

/*
 * Interpolation.
 *
 * The tables are bilinear in (logrho, logT) and everything is done by the
 * interpBilinear core on Mat->grid: forward lookups, the closed-form
 * inversion of the piecewise linear line obtained by fixing one coordinate
 * (T(rho, u), T(rho, s), rho(P, T)), the partials on the bilinear patch and
 * the node derivatives of the sound-speed table. The core is agnostic about
 * the logarithms: this file passes log10 values in and exponentiates on the
 * way out.
 *
 * The inversions use the strict mode: a target outside [f(T_min), f(T_max)]
 * of the column (or [f(rho_min), f(rho_max)] of the row) fails with
 * INTERP_BELOW/INTERP_ABOVE, any target inside gives a root. On a
 * non-monotonic column or row (a table defect, reported once by
 * scvheosCheckMonotonic at init) that is the root the bisection lands on, not
 * necessarily the physical one, and a target attained only in a dip below the
 * first or a rise above the last node fails like one outside. The value
 * tolerance absorbs a target that left the column by a 10^x <-> log10 round
 * trip.
 */
static const INTERPINVERTOPT scvheosInvertOpt = { 0, SCVHEOS_EDGE_TOL, 0.0, 0.0 };   /* INTERP_INVERT_STRICT with this library's tolerance */

/* Name of an inversion status, for the diagnostics. */
static const char *scvheosInvertStatus(int iret) {
    switch (iret) {
        case INTERP_BELOW:        return "below the table";
        case INTERP_ABOVE:        return "above the table";
        default:                  return "not attained";
    }
}

/*
 * Locate (logrho, logT) on the grid with the edge tolerance. Returns
 * INTERP_SUCCESS and the cell, or the core's status code (the caller prints).
 */
static int scvheosLocateLog(const SCVHEOSMAT *Mat, double logrho, double logT, INTERPCELL *cell) {
    return interpLocateTol(&Mat->grid, logrho, logT, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, cell);
}

/*
 * Check the monotonicity the inversions rely on, once per table at init:
 * u(T) and s(T) increasing along every rho column, P(rho) increasing along
 * every T row. Bilinear interpolation between two monotonic node columns is
 * monotonic, so checking the node columns covers every logrho. Reports the
 * number of defective columns/rows and the first few (with -DEOSLIB_VERBOSE only:
 * every thread initialises its own instance, so the report repeated per thread). The inversions do not
 * detect them (that would cost O(n) per call): there they return a root that
 * need not be the physical one, so a suspicious T(rho, u), T(rho, s) or
 * rho(P, T) can be traced to the table with this report.
 */
static void scvheosCheckMonotonic(SCVHEOSMAT *Mat) {
    const int nT = Mat->nT, nRho = Mat->nRho, nReport = 5;
    int nBadU = 0, nBadS = 0, nBadP = 0;

    for (int j = 0; j < nRho; j++) {
        int badU = interpNodeLineMonotonic(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, j);
        int badS = interpNodeLineMonotonic(&Mat->grid, Mat->dLogSArray, INTERP_ALONG_Y, j);
        if (badU) { if (nBadU < nReport) SCVHEOS_WARN("scvheosInitMaterial: iMat= %i: u(T) not increasing at jRho= %i (logrho= %g), e.g. iT= %i.\n", Mat->iMat, j, Mat->dLogRhoAxis[j], badU); nBadU++; }
        if (badS) { if (nBadS < nReport) SCVHEOS_WARN("scvheosInitMaterial: iMat= %i: s(T) not increasing at jRho= %i (logrho= %g), e.g. iT= %i.\n", Mat->iMat, j, Mat->dLogRhoAxis[j], badS); nBadS++; }
    }
    for (int i = 0; i < nT; i++) {
        int badP = interpNodeLineMonotonic(&Mat->grid, Mat->dLogPArray, INTERP_ALONG_X, i);
        if (badP) { if (nBadP < nReport) SCVHEOS_WARN("scvheosInitMaterial: iMat= %i: P(rho) not increasing at iT= %i (logT= %g), e.g. jRho= %i.\n", Mat->iMat, i, Mat->dLogTAxis[i], badP); nBadP++; }
    }
    if (nBadU || nBadS || nBadP) {
        SCVHEOS_WARN("scvheosInitMaterial: iMat= %i: non-monotonic columns: u(T) %i of %i, s(T) %i of %i, P(rho) %i of %i; "
                "the inversions T(rho, u), T(rho, s) and rho(P, T) may return another root there.\n",
                Mat->iMat, nBadU, nRho, nBadS, nRho, nBadP, nT);
    }
}

/*
 * Initalize a material.
 *
 * Initialize the material data, read the EOS tables and convert them to code units if desired.
 */
SCVHEOSMAT *scvheosInitMaterial(int iMat, double dKpcUnit, double dMsolUnit) {
    /* Gasoline/pkdgrav3 unit constants. Must match the hydro code and the other EOS libraries exactly. */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    SCVHEOSMAT *Mat;
    int nRho;
    int nT;
    char inFile[256];
    int nSkip;

    /* 
     * Allocate memory and initialize the data.
     *
     * Memory for the EOS tables is allocated in scvheosReadTable. The number of data points nRho
     * and nT are set in the function if the table was read correctly.
     */
    /*
     * Unit system (the same rule in EOSlib and every sublibrary): both units > 0 and finite are code
     * units, both <= 0 is cgs; anything else (exactly one of them <= 0, NaN, inf) is a caller error.
     */
    if (!(((dKpcUnit > 0.0) && (dMsolUnit > 0.0) && isfinite(dKpcUnit) && isfinite(dMsolUnit)) || ((dKpcUnit <= 0.0) && (dMsolUnit <= 0.0)))) {
        fprintf(stderr, "scvheosInitMaterial: invalid unit system dKpcUnit= %g dMsolUnit= %g (give both > 0 for code units or both <= 0 for cgs).\n", dKpcUnit, dMsolUnit);
        assert(0);
        exit(1);
    }

    Mat = (SCVHEOSMAT *) calloc(1, sizeof(SCVHEOSMAT));
    assert(Mat != NULL);

    Mat->iMat = iMat;
    Mat->dKpcUnit = dKpcUnit;
    Mat->dMsolUnit = dMsolUnit;
    Mat->nRho = 0;
    Mat->nT = 0;

    /* Reference density, a convention for a gas; see SCVHEOS_RHO0_CGS. */
    Mat->rho0 = SCVHEOS_RHO0_CGS;

    /*
     * Load the EOS table.
     */
    switch(iMat)
    {
        case SCVHEOS_H:
            /*
             * Hydrogen.
             */
            strcpy(inFile, "scvh_h_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            strcpy(Mat->MatString, "SCvH EOS H (Saumon et al. 1995).");
            break;
        case SCVHEOS_HE:
            /*
             * Helium.
             */
            strcpy(inFile, "scvh_he_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            strcpy(Mat->MatString, "SCvH EOS He (Saumon et al. 1995).");
            break;
        case SCVHEOS_HHE:
            /*
             * Hydrogen / Helium mixture in Solar abundance.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            strcpy(Mat->MatString, "SCvH EOS H-He Y=0.275 (Saumon et al. 1995).");
            break;
        case SCVHEOS_HHE_LOWRHOT:
            /*
             * Hydrogen / Helium mixture in Solar abundance limited to low rho and T.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs_lowrhot.txt");
            nRho = 80;
            nT = 48;
            nSkip = 2;
            strcpy(Mat->MatString, "SCvH EOS H-He Y=0.275 (Saumon et al. 1995, Vazan et al. 2013).");
            break;
        case SCVHEOS_HHE_EXT:
            /*
             * Hydrogen / Helium mixture (X=0.722, Y=0.278), full extended EOS table
             * (the low-rho/T cut this id once loaded is no longer used).
             */
            strcpy(inFile, "scvh_extended_dt_hydrogen_722_helium_278.data");
            nRho = 460;
            nT = 76;
            nSkip = 1;

            strcpy(Mat->MatString, "SCvH EOS H-He X=0.722 Y=0.278 (Saumon et al. 1995, Vazan et al. 2013).");
            break;
        case SCVHEOS_HE_EXT:
            /*
             * Helium EOS table extended to lower pressures and temperatures .
             */
            strcpy(inFile, "helium_scvh_extended_dt.data");
            nRho = 460;
            nT = 76;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS He extended to lower pressures and temperatures (Saumon et al. 1995, Vazan et al. 2013).");
            break;
        case SCVHEOS_HHE_H2O_Z_1:
            /*
             * H-He and H2O mixture obtained from Simon.
             */
            strcpy(inFile, "scvh_h2o_mixture_table_Z0.100.data");
            nRho = 128;
            nT = 64;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS H-He and H2O mixture Z=0.1.");
            break;
        case SCVHEOS_HHE_H2O_Z_05:
            /*
             * H-He and H2O mixture obtained from Simon.
             */
            strcpy(inFile, "scvh_h2o_mixture_table_Z0.050.data");
            nRho = 128;
            nT = 64;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS H-He and H2O mixture Z=0.05.");
            break;
        case SCVHEOS_HHE_H2O_Z_012:
            /*
             * H-He and H2O mixture obtained from Simon.
             */
            strcpy(inFile, "scvh_h2o_mixture_table_Z0.012.data");
            nRho = 128;
            nT = 64;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS H-He and H2O mixture Z=0.012.");
            break;
        default:
            /* Unknown material. Do not fall through: Mat is freed here. */
            fprintf(stderr, "scvheosInitMaterial: iMat= %i is not an SCvH EOS material.\n", iMat);
            scvheosFinalizeMaterial(Mat);
            assert(0);
            exit(1);
    }

    /*
     * Allocate memory and read the EOS table.
     */
    if (scvheosReadTable(Mat, inFile, nRho, nT, nSkip) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "scvheosInitMaterial: Could not read EOS table %s (see above).\n", inFile);
        scvheosFinalizeMaterial(Mat);
        assert(0);
        exit(1);
    }

    /* The view the interpBilinear core works on: x = logrho (outer), y = logT (inner). */
    Mat->grid.nx = Mat->nRho;
    Mat->grid.ny = Mat->nT;
    Mat->grid.x = Mat->dLogRhoAxis;
    Mat->grid.y = Mat->dLogTAxis;

    /* The lookup limits are the table bounds; no extrapolation. */
    Mat->LogRhoMin = Mat->dLogRhoAxis[0];
    Mat->LogRhoMax = Mat->dLogRhoAxis[Mat->nRho-1];
    Mat->LogTMin = Mat->dLogTAxis[0];
    Mat->LogTMax = Mat->dLogTAxis[Mat->nT-1];

    /*
     * rho0 must be a density the EOS can be evaluated at (EOSSolveBC starts
     * its lower bracket at 0.999*rho0). If the table does not cover the
     * conventional value, fall back to the geometric centre of the rho axis.
     * Still in cgs here.
     */
    {
        double logrho0 = log10(Mat->rho0);
        if (!(logrho0 > Mat->LogRhoMin) || !(logrho0 < Mat->LogRhoMax)) {
            logrho0 = 0.5*(Mat->LogRhoMin + Mat->LogRhoMax);
            fprintf(stderr, "scvheosInitMaterial: iMat= %i: reference density %g g/cc outside the "
                    "table [%g, %g] g/cc, using %g g/cc instead.\n", iMat, Mat->rho0,
                    scvheosExp10(Mat->LogRhoMin), scvheosExp10(Mat->LogRhoMax),
                    scvheosExp10(logrho0));
            Mat->rho0 = scvheosExp10(logrho0);
        }
    }

    /*
     * Convert from cgs to code units.
     */
    if ((Mat->dKpcUnit > 0.0) && (Mat->dMsolUnit > 0.0)) {
        /* code energy per unit mass --> erg per g */
        Mat->dErgPerGmUnit = GCGS*Mat->dMsolUnit*MSOLG/(Mat->dKpcUnit*KPCCM);
        /* code density --> g per cc */
        Mat->dGmPerCcUnit = (Mat->dMsolUnit*MSOLG)/pow(Mat->dKpcUnit*KPCCM,3.0);

        for (int i=0; i<Mat->nRho; i++) {
            Mat->dLogRhoAxis[i] -= log10(Mat->dGmPerCcUnit);
        }
        
        /*
         * The entropy table is NOT converted: following the EOSlib convention
         * (see EOSSofRhoU/EOSSofRhoT in EOSlib.h) the specific entropy is
         * always in cgs, erg g^-1 K^-1, whatever the units of rho, u and P.
         * Only entropy differences within one material are meaningful anyway.
         */
        for (int i=0; i<Mat->nT; i++) {
            for (int j=0; j<Mat->nRho; j++) {
                Mat->dLogPArray[j*Mat->nT+i] -= log10(Mat->dErgPerGmUnit*Mat->dGmPerCcUnit);
                Mat->dLogUArray[j*Mat->nT+i] -= log10(Mat->dErgPerGmUnit);
            }
        }
 
        /* Convert table limits and reference density. */
        Mat->LogRhoMin -= log10(Mat->dGmPerCcUnit);
        Mat->LogRhoMax -= log10(Mat->dGmPerCcUnit);
        Mat->rho0 /= Mat->dGmPerCcUnit;
    } else {
        /* No unit system given: everything stays in cgs. */
        Mat->dErgPerGmUnit = 1.0;
        Mat->dGmPerCcUnit = 1.0;
    }

    /* The sound speed table is built from dLogPArray/dLogUArray directly; it only needs the tables in code units. */
    if (scvheosGenerateSoundSpeedTable(Mat) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "scvheosInitMaterial: Could not generate table for sound speed.\n");
        scvheosFinalizeMaterial(Mat);
        assert(0);
        exit(1);
    }

    scvheosCheckMonotonic(Mat);

    return Mat;
}

/*
 * Free memory. Safe to call on a partially initialized material (calloc
 * zeroes all pointers, free(NULL) is fine).
 */
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat) {
    if (Mat != NULL) {
        free(Mat->dLogRhoAxis);
        free(Mat->dLogTAxis);

        free(Mat->dLogUArray);
        free(Mat->dLogPArray);
        free(Mat->dLogSArray);
        free(Mat->dLogCArray);

        free(Mat);
    }

    return SCVHEOS_SUCCESS;
}

/*
 * Read an EOS table from a file.
 *
 * The function also allocates memory of the tables and sets nRho and nT.
 *
 * SCVHEOS: 
 *
 * H:    nRho = 201, nT = 100, logarithmic spacing (base 10)
 * He:   nRho = 201, nT = 100, logarithmic spacing (base 10)
 * H/He: nRho = 201, nT = 100, logarithmic spacing (base 10)
 */
int scvheosReadTable(SCVHEOSMAT *Mat, const char *chInFile, int nRho, int nT, int nSkip) {
    FILE *fp = NULL;
    char *chLine = NULL;
    size_t nCharMax = 256;
    int iRet;
    int iStatus = SCVHEOS_FAIL;
    int i, j;

    if (Mat == NULL)  {
        return SCVHEOS_FAIL;
    }

    if ((Mat->dLogRhoAxis != NULL) || (Mat->dLogTAxis != NULL) || (Mat->dLogUArray != NULL) ||
            (Mat->dLogPArray != NULL) || (Mat->dLogSArray != NULL)) {
        return SCVHEOS_FAIL;
    }

    /* A cell needs two nodes on each axis (interpBilinear's convention; reos3ReadTable checks the same). */
    if ((nRho < 2) || (nT < 2)) {
        fprintf(stderr, "scvheosReadTable: nRho= %i nT= %i (need at least 2 each).\n", nRho, nT);
        return SCVHEOS_FAIL;
    }

    /*
     * The table arrays stay allocated on failure: scvheosFinalizeMaterial
     * frees them. Only the local resources (fp, chLine) are released here,
     * on every path, via the cleanup label.
     */
    Mat->dLogRhoAxis = (double *) calloc(nRho, sizeof(double));
    Mat->dLogTAxis = (double *) calloc(nT, sizeof(double));

    /* The tables are flat rho-major arrays, A[j*nT + i]. */
    Mat->dLogUArray = (double *) calloc(nT*nRho, sizeof(double));
    Mat->dLogPArray = (double *) calloc(nT*nRho, sizeof(double));
    Mat->dLogSArray = (double *) calloc(nT*nRho, sizeof(double));

    if ((Mat->dLogRhoAxis == NULL) || (Mat->dLogTAxis == NULL) || (Mat->dLogUArray == NULL) ||
            (Mat->dLogPArray == NULL) || (Mat->dLogSArray == NULL)) {
        goto cleanup;
    }

    chLine = (char *) calloc(nCharMax, sizeof(char));
    if (chLine == NULL) {
        goto cleanup;
    }

    /* Open the file and skip the first few lines. */
    fp = fopen(chInFile, "r");
    if (fp == NULL) {
        goto cleanup;
    }

    for (i=0; i<nSkip; i++) {
        if (getline(&chLine, &nCharMax, fp) == -1) {
            goto cleanup;
        }
    }

    /*
     * The file lists one (logT, logrho, logP, logu, logs) row per node, T
     * outermost. The axes are taken from the first T row (for rho) and the
     * first rho column (for T) and every later row must repeat them exactly:
     * a mismatch means the file is not a rectangular grid in the assumed
     * order, which the interpolator would silently misread.
     */
    for (i=0; i<nT; i++) {
        for (j=0; j<nRho; j++) {
            double logT, logrho, logP, logu, logs;

            if (getline(&chLine, &nCharMax, fp) == -1) {
                fprintf(stderr, "scvheosReadTable: %s: unexpected end of file at iT= %i jRho= %i "
                        "(expected nT= %i x nRho= %i nodes).\n", chInFile, i, j, nT, nRho);
                goto cleanup;
            }

            iRet = sscanf(chLine, "%lf %lf %lf %lf %lf", &logT, &logrho, &logP, &logu, &logs);

            /* Check if the number of matches is correct. */
            if (iRet != 5) {
                fprintf(stderr, "scvheosReadTable: %s: expected 5 columns at iT= %i jRho= %i, got %i.\n",
                        chInFile, i, j, iRet);
                goto cleanup;
            }

            /* The file stores logarithms; nan or inf would poison the interpolator. */
            if (!isfinite(logT) || !isfinite(logrho) || !isfinite(logP) || !isfinite(logu) || !isfinite(logs)) {
                fprintf(stderr, "scvheosReadTable: %s: non-finite entry at iT= %i jRho= %i "
                        "(logT= %g logrho= %g logP= %g logu= %g logs= %g).\n",
                        chInFile, i, j, logT, logrho, logP, logu, logs);
                goto cleanup;
            }

            /* Axes: set from the first row/column, verify afterwards. */
            if (i == 0) {
                Mat->dLogRhoAxis[j] = logrho;
            } else if (logrho != Mat->dLogRhoAxis[j]) {
                fprintf(stderr, "scvheosReadTable: %s: logrho= %.15g at iT= %i jRho= %i differs from "
                        "the rho axis (%.15g): not a rectangular grid.\n",
                        chInFile, logrho, i, j, Mat->dLogRhoAxis[j]);
                goto cleanup;
            }

            if (j == 0) {
                Mat->dLogTAxis[i] = logT;
            } else if (logT != Mat->dLogTAxis[i]) {
                fprintf(stderr, "scvheosReadTable: %s: logT= %.15g at iT= %i jRho= %i differs from "
                        "the T axis (%.15g): not a rectangular grid.\n",
                        chInFile, logT, i, j, Mat->dLogTAxis[i]);
                goto cleanup;
            }

            /* We store the data as A[T][rho]. */
            Mat->dLogPArray[j*nT+i] = logP;
            Mat->dLogUArray[j*nT+i] = logu;
            Mat->dLogSArray[j*nT+i] = logs;
        }
    }

    /* The interpolation requires strictly increasing axes. */
    for (j=1; j<nRho; j++) {
        if (!(Mat->dLogRhoAxis[j] > Mat->dLogRhoAxis[j-1])) {
            fprintf(stderr, "scvheosReadTable: %s: rho axis not strictly increasing at jRho= %i "
                    "(%.15g <= %.15g).\n", chInFile, j, Mat->dLogRhoAxis[j], Mat->dLogRhoAxis[j-1]);
            goto cleanup;
        }
    }
    for (i=1; i<nT; i++) {
        if (!(Mat->dLogTAxis[i] > Mat->dLogTAxis[i-1])) {
            fprintf(stderr, "scvheosReadTable: %s: T axis not strictly increasing at iT= %i "
                    "(%.15g <= %.15g).\n", chInFile, i, Mat->dLogTAxis[i], Mat->dLogTAxis[i-1]);
            goto cleanup;
        }
    }

    Mat->nRho = nRho;
    Mat->nT = nT;
    iStatus = SCVHEOS_SUCCESS;

cleanup:
    if (fp != NULL) fclose(fp);
    free(chLine);   /* getline may have reallocated it; free(NULL) is fine */

    return iStatus;
}

/*
 * Generate an array that contains the sound speed at each EOS table data point.
 *
 * The sound speed is calculated from
 *
 * cs^2 = dP/drho + T/(rho^2*C_v)*(dP/dT)^2
 *
 * where
 *
 * c_v = du/dT
 *
 * is the specific heat capacity at constant volume. The derivatives are the
 * central differences of the log tables on the grid (interpNodeDeriv: second
 * order, symmetric about the node, one-sided at the ends), taken directly
 * from dLogPArray and dLogUArray:
 * dP/drho = (P/rho) dlogP/dlogrho, dP/dT = (P/T) dlogP/dlogT, c_v = (u/T) dlogu/dlogT.
 *
 * Where the table is defective (c_v <= 0 or a negative cs^2, which happens
 * e.g. in the H-He/H2O Z=0.05 table) cs^2 is replaced by P/rho at the node:
 * this has the right units whatever the table is in, and for a gas it is a
 * lower bound on cs^2 = gamma*P/rho, so the hydro code errs towards a short
 * timestep. A zero or -inf must never be stored, because bilinear
 * interpolation would spread it into the neighbouring cells. The number of
 * replaced nodes is reported on stderr.
 */
int scvheosGenerateSoundSpeedTable(SCVHEOSMAT *Mat) {
    int nT, nRho;
    double *rhoNode, *TNode;
    int nReplaced = 0;
    const int nReport = 10;   /* print at most this many defective nodes */

    /* Check before the first dereference (the sizes were read in the initializers before). */
    if (Mat == NULL)
        return SCVHEOS_FAIL;
    nT = Mat->nT;
    nRho = Mat->nRho;

    if ((Mat->dLogRhoAxis == NULL) || (Mat->dLogTAxis == NULL) || (Mat->dLogPArray == NULL) || (Mat->dLogUArray == NULL))
        return SCVHEOS_FAIL;

    if (Mat->dLogCArray != NULL)
        return SCVHEOS_FAIL;

    Mat->dLogCArray = (double *) calloc((size_t)nT*nRho, sizeof(double));
    rhoNode = (double *) malloc((size_t)nRho*sizeof(double));
    TNode = (double *) malloc((size_t)nT*sizeof(double));

    if ((Mat->dLogCArray == NULL) || (rhoNode == NULL) || (TNode == NULL)) {
        free(rhoNode);
        free(TNode);
        return SCVHEOS_FAIL;
    }

    for (int j = 0; j < nRho; j++) rhoNode[j] = scvheosExp10(Mat->dLogRhoAxis[j]);
    for (int i = 0; i < nT; i++) TNode[i] = scvheosExp10(Mat->dLogTAxis[i]);

    for (int j = 0; j < nRho; j++) {
        const double *logP = Mat->dLogPArray + (size_t)j*nT;   /* column at rho_j, consecutive in T */
        const double *logu = Mat->dLogUArray + (size_t)j*nT;
        const double rho = rhoNode[j];

        for (int i = 0; i < nT; i++) {
            const double T = TNode[i];
            const double P = scvheosExp10(logP[i]);
            const double u = scvheosExp10(logu[i]);
            /* Along rho the nodes are nT apart, along T they are consecutive. */
            const double dlPdlr = interpNodeDeriv(Mat->dLogRhoAxis, nRho, Mat->dLogPArray + i, (size_t)nT, j);
            const double dlPdlT = interpNodeDeriv(Mat->dLogTAxis, nT, logP, 1, i);
            const double dludlT = interpNodeDeriv(Mat->dLogTAxis, nT, logu, 1, i);
            const double dPdrho = P/rho*dlPdlr;
            const double dPdT = P/T*dlPdlT;
            const double cv = u/T*dludlT;
            double cs2 = dPdrho + T/(rho*rho*cv)*dPdT*dPdT;

            /* Defective node: non-finite (cv = 0), or not positive. */
            if (!isfinite(cs2) || !(cs2 > 0.0)) {
                if (nReplaced < nReport) {
                    SCVHEOS_WARN("scvheosGenerateSoundSpeedTable: iMat= %i iT= %i jRho= %i "
                            "logrho= %g logT= %g dPdRho= %g dPdT= %g cv= %g cs2= %g, "
                            "replaced by P/rho= %g\n", Mat->iMat, i, j,
                            Mat->dLogRhoAxis[j], Mat->dLogTAxis[i], dPdrho, dPdT, cv, cs2, P/rho);
                }
                nReplaced++;
                cs2 = P/rho;
            }

            /* P/rho > 0 at every node (the table stores log P), so this is finite. */
            assert(isfinite(cs2) && cs2 > 0.0);

            Mat->dLogCArray[(size_t)j*nT + i] = 0.5*log10(cs2);
        }
    }

    free(rhoNode);
    free(TNode);

    if (nReplaced > 0) {
        SCVHEOS_WARN("scvheosGenerateSoundSpeedTable: iMat= %i: cs^2 replaced by P/rho at %i of %i nodes.\n",
                Mat->iMat, nReplaced, nT*nRho);
    }

    return SCVHEOS_SUCCESS;
}

/*
 * Calculate logP(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (counted in nLookupFail) if (logrho, logT) is
 * not finite or outside of the table (closed intervals, see SCVHEOS_EDGE_TOL).
 */
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    INTERPCELL cell;

    if (scvheosLocateLog(Mat, logrho, logT, &cell) != INTERP_SUCCESS) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogPofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
    }

    return interpCellValue(&Mat->grid, Mat->dLogPArray, &cell);
}

/*
 * Calculate logU(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (counted in nLookupFail) if (logrho, logT) is
 * not finite or outside of the table.
 */
double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    INTERPCELL cell;

    if (scvheosLocateLog(Mat, logrho, logT, &cell) != INTERP_SUCCESS) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogUofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
    }

    return interpCellValue(&Mat->grid, Mat->dLogUArray, &cell);
}

/*
 * Calculate logs(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (counted in nLookupFail) if (logrho, logT) is
 * not finite or outside of the table.
 */
double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    INTERPCELL cell;

    if (scvheosLocateLog(Mat, logrho, logT, &cell) != INTERP_SUCCESS) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogSofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
    }

    return interpCellValue(&Mat->grid, Mat->dLogSArray, &cell);
}

/*
 * Calculate logcs(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (counted in nLookupFail) if (logrho, logT) is
 * not finite or outside of the table.
 */
double scvheosLogCsofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    INTERPCELL cell;

    if (scvheosLocateLog(Mat, logrho, logT, &cell) != INTERP_SUCCESS) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogCsofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
    }

    return interpCellValue(&Mat->grid, Mat->dLogCArray, &cell);
}


/*
 * Calculate the pressure P(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logP;

    logrho = log10(rho);
    logT = log10(T);

    /* Interpolate in the table. */
    logP = scvheosLogPofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logP)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logP);
}

/*
 * Calculate the internal energy u(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logu;

    logrho = log10(rho);
    logT = log10(T);

    logu = scvheosLogUofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logu)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logu);
}

/*
 * Calculate the specific entropy s(rho, T).
 *
 * Returned in cgs (erg g^-1 K^-1) regardless of the code units of rho, u
 * and P (EOSlib convention). Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosSofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logs;

    logrho = log10(rho);
    logT = log10(T);

    /* Interpolate in the table. */
    logs = scvheosLogSofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logs)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logs);
}

/*
 * Calculate the sound speed cs(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logcs;

    logrho = log10(rho);
    logT = log10(T);

    logcs = scvheosLogCsofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logcs)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logcs);
}

/*
 * Locate the state (logrho, logu) in the table: the cell (ix along rho, iy
 * along T) with its local coordinates, from the inversion of the u column,
 * plus logT. Shared by scvheosLogTofLogRhoLogU, scvheosPCTofRhoU and the
 * (rho, u) derivatives.
 *
 * Returns 0 on success, -1 (counted in nLookupFail, message under
 * EOSLIB_VERBOSE) if logrho is not finite (rho <= 0) or outside the table.
 *
 * If logrho is in the table but logu cannot be inverted (not finite, i.e.
 * u <= 0, or outside the column [u(T_min), u(T_max)]) the state falls back
 * to T_min, the cold curve of the table: cell iy = 0, v = 0, logT = LogTMin.
 * This is the EOSlib contract for T(rho, u) (same as ANEOSTofRhoU); the
 * fallback is counted in nLookupFail and reported under EOSLIB_VERBOSE, and
 * the function returns 0, so P, cs and the derivatives are then those of
 * (rho, T_min). On a non-monotonic column a logu inside [u(T_min), u(T_max)]
 * gives a root, possibly not the physical one (see "Interpolation" above).
 */
static int scvheosLocateLogRhoLogU(SCVHEOSMAT *Mat, const char *chFunc, double logrho, double logu,
                                   INTERPCELL *cell, double *plogT) {
    int iret;

    if (!isfinite(logrho)) {
        Mat->nLookupFail++;
        SCVHEOS_WARN("%s: invalid input logrho= %g logu= %g.\n", chFunc, logrho, logu);
        return -1;
    }

    /* rho first: outside the table nothing can be evaluated. */
    if (interpLocateAxisTol(Mat->dLogRhoAxis, Mat->nRho, logrho, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, &cell->ix, &cell->u) != INTERP_SUCCESS) {
        Mat->nLookupFail++;
        SCVHEOS_WARN("%s: logrho= %15.7E outside of the EOS table (logrho in [%g, %g]).\n", chFunc, logrho, Mat->LogRhoMin, Mat->LogRhoMax);
        return -1;
    }

    iret = isfinite(logu)
         ? interpInvertLine(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, cell->ix, cell->u, logu, &scvheosInvertOpt, &cell->iy, &cell->v)
         : INTERP_INVALID;

    if (iret != INTERP_SUCCESS) {
        /* Fallback, by contract: the lowest tabulated temperature. */
        Mat->nLookupFail++;
        SCVHEOS_WARN("%s: logu= %g %s at logrho= %g "
                "(column u(T_min)= %g u(T_max)= %g), falling back to logT_min= %g.\n", chFunc, logu,
                (iret == INTERP_INVALID) ? "invalid" : scvheosInvertStatus(iret), logrho,
                interpLineValue(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, cell->ix, cell->u, 0),
                interpLineValue(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, cell->ix, cell->u, Mat->nT - 1),
                Mat->LogTMin);
        cell->iy = 0;
        cell->v = 0.0;
    }

    *plogT = interpLineCoord(&Mat->grid, INTERP_ALONG_Y, cell->iy, cell->v);
    return 0;
}

/* Same in linear space, for scvheosPCTofRhoU and the (rho, u) derivatives. */
static int scvheosLocateRhoU(SCVHEOSMAT *Mat, const char *chFunc, double rho, double u,
                             INTERPCELL *cell, double *plogT) {
    return scvheosLocateLogRhoLogU(Mat, chFunc, log10(rho), log10(u), cell, plogT);
}

/*
 * Calculate the pressure P(rho, u), the sound speed cs(rho, u) and the
 * temperature T(rho, u) from one cell lookup.
 *
 * The inversion of the u column gives the T segment; P and cs are then the
 * bilinear values of their tables in that same cell, with no further search.
 * This is the hydro-code entry point (EOSPCTofRhoU) and mirrors
 * ANEOSPCTofRhoU.
 *
 * If rho is outside the table returns SCVHEOS_FAIL_VALUE and sets *c to it
 * as well, and *T to T_min (T(rho, u) never fails, see scvheosTofRhoU); if
 * u cannot be inverted the state falls back to (rho, T_min) (see
 * scvheosLocateLogRhoLogU).
 */
double scvheosPCTofRhoU(SCVHEOSMAT *Mat, double rho, double u, double *c, double *T) {
    INTERPCELL cell;
    double logT;

    *c = SCVHEOS_FAIL_VALUE;
    *T = scvheosExp10(Mat->LogTMin);   /* T(rho, u) never fails: T_min by contract */

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &cell, &logT) != 0) return SCVHEOS_FAIL_VALUE;

    *T = scvheosExp10(logT);
    *c = scvheosExp10(interpCellValue(&Mat->grid, Mat->dLogCArray, &cell));
    return scvheosExp10(interpCellValue(&Mat->grid, Mat->dLogPArray, &cell));
}

/*
 * Calculate the pressure P(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if rho is outside the table; P(rho, T_min) if
 * u cannot be inverted (see scvheosLocateLogRhoLogU).
 */
double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double c, T;
    return scvheosPCTofRhoU(Mat, rho, u, &c, &T);
}

/*
 * Calculate the sound speed cs(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if rho is outside the table; cs(rho, T_min) if
 * u cannot be inverted (see scvheosLocateLogRhoLogU).
 */
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double c, T;
    scvheosPCTofRhoU(Mat, rho, u, &c, &T);
    return c;
}

/*
 * Calculate the specific entropy s(rho, u) in cgs (EOSlib convention, see
 * scvheosSofRhoT), evaluated in the cell of the u-column inversion like
 * scvheosPCTofRhoU.
 *
 * Returns SCVHEOS_FAIL_VALUE if rho is outside the table; s(rho, T_min) if
 * u cannot be inverted. Either is counted once in nLookupFail (EOSlib built
 * s(rho, u) from scvheosTofRhoU and scvheosSofRhoT, which counted a density
 * outside the table twice: the T fallback and the failed s lookup).
 */
double scvheosSofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    INTERPCELL cell;
    double logT;

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &cell, &logT) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosExp10(interpCellValue(&Mat->grid, Mat->dLogSArray, &cell));
}

/*
 * Calculate the temperature logT(logrho, logu) by direct inversion of the
 * u column at logrho.
 *
 * Never returns the sentinel (same contract as ANEOSTofRhoU): if logrho is
 * not finite or outside the table, or logu is not finite or outside the
 * column [u(T_min), u(T_max)], returns logT_min and counts the fallback (see
 * scvheosLocateLogRhoLogU). On a non-monotonic column it returns a root, not
 * necessarily the physical one.
 */
double scvheosLogTofLogRhoLogU(SCVHEOSMAT *Mat, double logrho, double logu) {
    INTERPCELL cell;
    double logT;

    if (scvheosLocateLogRhoLogU(Mat, __func__, logrho, logu, &cell, &logT) != 0) return Mat->LogTMin;
    return logT;
}

/*
 * Calculate the temperature logT(logrho, logs) by direct inversion of the
 * s column at logrho (interpInvertY, strict mode).
 *
 * logs below s(T_min) or above s(T_max) of the column returns logT_min or
 * logT_max respectively, counted in Mat->nLookupFail (ANEOSmaterial
 * convention). Returns SCVHEOS_FAIL_VALUE if the input is not finite or
 * logrho is outside the table. On a non-monotonic s column a logs inside
 * [s(T_min), s(T_max)] gives a root, not necessarily the physical one.
 */
double scvheosLogTofLogRhoLogS(SCVHEOSMAT *Mat, double logrho, double logs) {
    INTERPCELL cell;
    double logT;
    int iret;

    if (!isfinite(logrho) || !isfinite(logs)) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogTofLogRhoLogS: invalid input logrho= %g logs= %g.\n", logrho, logs);
    }

    iret = interpInvertYTol(&Mat->grid, Mat->dLogSArray, logrho, logs, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, &scvheosInvertOpt, &cell, &logT);

    if (iret == INTERP_OUTSIDE_XMIN || iret == INTERP_OUTSIDE_XMAX) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogTofLogRhoLogS: logrho= %15.7E outside of the EOS table.\n", logrho);
    }

    /*
     * logs outside the column [s(T_min), s(T_max)]: clamp to the table edge,
     * the convention ANEOSisentropicU uses (cold curve u(rho2, T_min) when S
     * is not attained at rho2). Unlike there the fallback is counted in
     * Mat->nLookupFail and reported with SCVHEOS_WARN.
     */
    if (iret == INTERP_BELOW || iret == INTERP_ABOVE) {
        int ix = 0; double w = 0.0;
        interpLocateAxisTol(Mat->dLogRhoAxis, Mat->nRho, logrho, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, &ix, &w);
        Mat->nLookupFail++;
        if (iret == INTERP_BELOW) {
            SCVHEOS_WARN("scvheosLogTofLogRhoLogS: logs= %g below s(T_min)= %g at logrho= %g, "
                         "falling back to logT_min= %g\n", logs,
                         interpLineValue(&Mat->grid, Mat->dLogSArray, INTERP_ALONG_Y, ix, w, 0), logrho, Mat->LogTMin);
            return Mat->LogTMin;
        }
        SCVHEOS_WARN("scvheosLogTofLogRhoLogS: logs= %g above s(T_max)= %g at logrho= %g, "
                     "falling back to logT_max= %g\n", logs,
                     interpLineValue(&Mat->grid, Mat->dLogSArray, INTERP_ALONG_Y, ix, w, Mat->nT - 1), logrho, Mat->LogTMax);
        return Mat->LogTMax;
    }

    /* Every other status (INTERP_INVALID) needs a NaN input, rejected above. */
    if (iret != INTERP_SUCCESS) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogTofLogRhoLogS: inversion failed (status %d) at logs= %g logrho= %g.\n", iret, logs, logrho);
    }

    return logT;
}

/*
 * Calculate the density logrho(logP, logT) by direct inversion of the P row
 * at logT (interpInvertX, strict mode).
 *
 * Returns SCVHEOS_FAIL_VALUE (counted in nLookupFail) if the input is not
 * finite, logT is outside the table, or logP is outside the row
 * [P(rho_min), P(rho_max)]. Where P(rho) is not monotonic the inversion
 * returns a root, not necessarily the physical one.
 */
double scvheosLogRhoofLogPLogT(SCVHEOSMAT *Mat, double logP, double logT) {
    INTERPCELL cell;
    double logrho;
    int iret;

    if (!isfinite(logP) || !isfinite(logT)) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogRhoofLogPLogT: invalid input logP= %g logT= %g.\n", logP, logT);
    }

    iret = interpInvertXTol(&Mat->grid, Mat->dLogPArray, logT, logP, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, &scvheosInvertOpt, &cell, &logrho);

    if (iret == INTERP_OUTSIDE_YMIN || iret == INTERP_OUTSIDE_YMAX) {
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogRhoofLogPLogT: logT= %15.7E outside of the EOS table.\n", logT);
    }
    if (iret != INTERP_SUCCESS) {
        int iy = 0; double v = 0.0;
        interpLocateAxisTol(Mat->dLogTAxis, Mat->nT, logT, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, &iy, &v);
        SCVHEOS_FAIL_RETURN(Mat, "scvheosLogRhoofLogPLogT: logP= %g %s at logT= %g "
                "(row P(rho_min)= %g P(rho_max)= %g, logrho in [%g, %g]).\n", logP,
                scvheosInvertStatus(iret), logT,
                interpLineValue(&Mat->grid, Mat->dLogPArray, INTERP_ALONG_X, iy, v, 0),
                interpLineValue(&Mat->grid, Mat->dLogPArray, INTERP_ALONG_X, iy, v, Mat->nRho - 1),
                Mat->LogRhoMin, Mat->LogRhoMax);
    }

    return logrho;
}

/*
 * Calculate the temperature T(rho, u). Never fails: returns T_min if rho is
 * outside the table or u cannot be inverted (see scvheosLogTofLogRhoLogU).
 */
double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    return scvheosExp10(scvheosLogTofLogRhoLogU(Mat, log10(rho), log10(u)));
}

/*
 * Calculate the temperature T(rho, s). s is in cgs (erg g^-1 K^-1), as
 * returned by scvheosSofRhoT.
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion fails.
 */
double scvheosTofRhoS(SCVHEOSMAT *Mat, double rho, double s) {
    double logrho;
    double logs;
    double logT;

    logrho = log10(rho);
    logs = log10(s);

    logT = scvheosLogTofLogRhoLogS(Mat, logrho, logs);
    if (scvheosFailed(logT)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logT);
}

/*
 * Calculate the density rho(P, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion fails (callers such as
 * EOSWoolfsonCoeff test for rho <= 0).
 */
double scvheosRhoofPT(SCVHEOSMAT *Mat, double P, double T) {
    double logP;
    double logT;
    double logrho;

    logP = log10(P);
    logT = log10(T);

    logrho = scvheosLogRhoofLogPLogT(Mat, logP, logT);
    if (scvheosFailed(logrho)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logrho);
}

/*
 * Partial derivatives of the tabulated log quantities with respect to logrho
 * (bWrtRho = TRUE) or logT (bWrtRho = FALSE), evaluated on the bilinear
 * patch that contains (logrho, logT). Only defined inside the table: the
 * lookup limits are the table bounds and the bilinear patch has no
 * extrapolation. At a node the cell to the right/above is used (the last
 * cell at the upper table edge). Returns SCVHEOS_FAIL_VALUE (counted in
 * nLookupFail) outside the table.
 */
static double scvheosdLogdLog(SCVHEOSMAT *Mat, const double *Array,
                              double logrho, double logT, int bWrtRho, const char *chFunc) {
    INTERPCELL cell;
    double dlr, dlT;

    if (scvheosLocateLog(Mat, logrho, logT, &cell) != INTERP_SUCCESS) {
        SCVHEOS_FAIL_RETURN(Mat, "%s: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", chFunc, logrho, logT);
    }

    interpCellGradient(&Mat->grid, Array, &cell, &dlr, &dlT);
    return bWrtRho ? dlr : dlT;
}

/*
 * Calculate the derivative dlogP/dlogrho(logrho, logT).
 */
double scvheosdLogPdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->dLogPArray, logrho, logT, TRUE, __func__);
}

/*
 * Calculate the derivative dlogP/dlogT(logrho, logT).
 */
double scvheosdLogPdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->dLogPArray, logrho, logT, FALSE, __func__);
}

/*
 * Calculate the derivative dlogu/dlogrho(logrho, logT).
 */
double scvheosdLogUdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->dLogUArray, logrho, logT, TRUE, __func__);
}

/*
 * Calculate the derivative dlogu/dlogT(logrho, logT).
 */
double scvheosdLogUdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->dLogUArray, logrho, logT, FALSE, __func__);
}

/*
 * Calculate the derivative dlogs/dlogrho(logrho, logT).
 */
double scvheosdLogSdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->dLogSArray, logrho, logT, TRUE, __func__);
}

/*
 * Calculate the derivative dlogs/dlogT(logrho, logT).
 */
double scvheosdLogSdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->dLogSArray, logrho, logT, FALSE, __func__);
}

/*
 * Locate (rho, T) in the table for the (rho, T) derivatives. Returns 0, or
 * -1 (counted in nLookupFail, message under EOSLIB_VERBOSE) if the state is
 * not finite or outside the table.
 */
static int scvheosLocateRhoT(SCVHEOSMAT *Mat, const char *chFunc, double rho, double T, INTERPCELL *cell) {
    double logrho, logT;

    if (!(rho > 0.0) || !(T > 0.0) || !isfinite(rho) || !isfinite(T)) {
        Mat->nLookupFail++;
        SCVHEOS_WARN("%s: invalid input rho= %g T= %g.\n", chFunc, rho, T);
        return -1;
    }
    logrho = log10(rho);
    logT = log10(T);
    if (scvheosLocateLog(Mat, logrho, logT, cell) != INTERP_SUCCESS) {
        Mat->nLookupFail++;
        SCVHEOS_WARN("%s: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", chFunc, logrho, logT);
        return -1;
    }
    return 0;
}

/*
 * (f/x) dlogf/dlogx on the patch: the linear-space derivative of one table
 * in the located cell with respect to rho (bWrtRho) or T, where x is that
 * coordinate's linear value. The (rho, T) derivatives below return
 * SCVHEOS_FAIL_VALUE (counted in nLookupFail) outside the table.
 */
static double scvheosdFdX(SCVHEOSMAT *Mat, const double *Array, const INTERPCELL *cell, double x, int bWrtRho) {
    double dlr, dlT;
    interpCellGradient(&Mat->grid, Array, cell, &dlr, &dlT);
    return scvheosExp10(interpCellValue(&Mat->grid, Array, cell))/x*(bWrtRho ? dlr : dlT);
}

/*
 * Calculate the derivative dP/drho(rho, T) = (P/rho) dlogP/dlogrho on the
 * bilinear patch.
 */
double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    INTERPCELL cell;
    if (scvheosLocateRhoT(Mat, __func__, rho, T, &cell) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosdFdX(Mat, Mat->dLogPArray, &cell, rho, TRUE);
}

/*
 * Calculate the derivative dP/dT(rho, T) = (P/T) dlogP/dlogT.
 */
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    INTERPCELL cell;
    if (scvheosLocateRhoT(Mat, __func__, rho, T, &cell) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosdFdX(Mat, Mat->dLogPArray, &cell, T, FALSE);
}

/*
 * Calculate the derivative du/drho(rho, T) = (u/rho) dlogu/dlogrho.
 */
double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    INTERPCELL cell;
    if (scvheosLocateRhoT(Mat, __func__, rho, T, &cell) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosdFdX(Mat, Mat->dLogUArray, &cell, rho, TRUE);
}

/*
 * Calculate the derivative du/dT(rho, T) = (u/T) dlogu/dlogT, the specific
 * heat capacity at constant volume.
 */
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    INTERPCELL cell;
    if (scvheosLocateRhoT(Mat, __func__, rho, T, &cell) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosdFdX(Mat, Mat->dLogUArray, &cell, T, FALSE);
}

/*
 * Calculate the derivative ds/drho(rho, T) = (s/rho) dlogs/dlogrho (s in cgs).
 */
double scvheosdSdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    INTERPCELL cell;
    if (scvheosLocateRhoT(Mat, __func__, rho, T, &cell) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosdFdX(Mat, Mat->dLogSArray, &cell, rho, TRUE);
}

/*
 * Calculate the derivative ds/dT(rho, T) = (s/T) dlogs/dlogT (s in cgs).
 */
double scvheosdSdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    INTERPCELL cell;
    if (scvheosLocateRhoT(Mat, __func__, rho, T, &cell) != 0) return SCVHEOS_FAIL_VALUE;
    return scvheosdFdX(Mat, Mat->dLogSArray, &cell, T, FALSE);
}

/*
 * The four partials of logP and logu on the patch, for the (rho, u) derivatives.
 */
static void scvheosCellPartials(const SCVHEOSMAT *Mat, const INTERPCELL *cell,
                                double *dlPdlr, double *dlPdlT, double *dludlr, double *dludlT) {
    interpCellGradient(&Mat->grid, Mat->dLogPArray, cell, dlPdlr, dlPdlT);
    interpCellGradient(&Mat->grid, Mat->dLogUArray, cell, dludlr, dludlT);
}

/*
 * Calculate dP/drho at constant u, analytically from the bilinear patch:
 *
 *   dP/drho|_u = dP/drho|_T - dP/dT|_rho * (du/drho|_T) / (du/dT|_rho)
 *
 * With the tables in log form, dP/drho|_T = (P/rho) dlogP/dlogrho etc., so
 *
 *   dP/drho|_u = (P/rho) [dlogP/dlogrho - dlogP/dlogT * dlogu/dlogrho / dlogu/dlogT].
 *
 * Returns SCVHEOS_FAIL_VALUE if rho is outside the table or du/dT|rho <= 0
 * in the cell (the heat capacity must be positive; a non-positive value is a
 * table artefact). If u cannot be inverted the derivative is that of
 * (rho, T_min) (see scvheosLocateLogRhoLogU).
 */
double scvheosdPdRhoofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    INTERPCELL cell;
    double logT, P;
    double dlPdlr, dlPdlT, dludlr, dludlT;
    const long nLookupFail0 = Mat->nLookupFail;

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &cell, &logT) != 0) return SCVHEOS_FAIL_VALUE;

    scvheosCellPartials(Mat, &cell, &dlPdlr, &dlPdlT, &dludlr, &dludlT);

    if (!(dludlT > 0.0)) {
        /* Counted once per call: not again after a T fallback (the cell at T_min can have dlogu/dlogT <= 0). */
        Mat->nLookupFail = nLookupFail0;
        SCVHEOS_FAIL_RETURN(Mat, "scvheosdPdRhoofRhoU: dlogu/dlogT= %g <= 0 at rho= %g u= %g (cell jRho= %i iT= %i).\n",
                dludlT, rho, u, cell.ix, cell.iy);
    }

    P = scvheosExp10(interpCellValue(&Mat->grid, Mat->dLogPArray, &cell));

    return P/rho*(dlPdlr - dlPdlT*dludlr/dludlT);
}

/*
 * Calculate dP/du at constant rho, analytically from the bilinear patch:
 *
 *   dP/du|_rho = dP/dT|_rho / du/dT|_rho = (P/u) dlogP/dlogT / dlogu/dlogT.
 *
 * Same failure conditions as scvheosdPdRhoofRhoU.
 */
double scvheosdPdUofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    INTERPCELL cell;
    double logT, P;
    double dlPdlr, dlPdlT, dludlr, dludlT;
    const long nLookupFail0 = Mat->nLookupFail;

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &cell, &logT) != 0) return SCVHEOS_FAIL_VALUE;

    scvheosCellPartials(Mat, &cell, &dlPdlr, &dlPdlT, &dludlr, &dludlT);

    if (!(dludlT > 0.0)) {
        /* Counted once per call: not again after a T fallback (the cell at T_min can have dlogu/dlogT <= 0). */
        Mat->nLookupFail = nLookupFail0;
        SCVHEOS_FAIL_RETURN(Mat, "scvheosdPdUofRhoU: dlogu/dlogT= %g <= 0 at rho= %g u= %g (cell jRho= %i iT= %i).\n",
                dludlT, rho, u, cell.ix, cell.iy);
    }

    P = scvheosExp10(interpCellValue(&Mat->grid, Mat->dLogPArray, &cell));

    return P/u*dlPdlT/dludlT;
}

/*
 * Calculate u2 so that (rho1, u1) and (rho2, u2) are on the same isentrope.
 *
 * If the entropy S(rho1, u1) is not attained at rho2 the state is clamped to
 * the table edge: u(rho2, T_min) (the cold curve) when it would fall below
 * the table, u(rho2, T_max) above it. This follows ANEOSisentropicU. The start
 * point is clamped the same way: (rho1, u1) above the u column starts from
 * s(rho1, T_max), the top row, and below the cold curve (or u1 <= 0, NaN) from
 * s(rho1, T_min), the T(rho, u) fallback. Up to 1.3.0 both started from
 * s(rho1, T_min), so a state above the table followed the coldest isentrope
 * and lost all its thermal energy; the energy above the top row is lost
 * either way. Returns SCVHEOS_FAIL_VALUE if rho1 or rho2 is outside the
 * table. A call counts at most once in nLookupFail, however many of these
 * steps fail or fall back (up to 1.3.0 rho1 outside the table counted twice).
 */
double scvheosIsentropicU(SCVHEOSMAT *Mat, double rho1, double u1, double rho2) {
    const long nLookupFail0 = Mat->nLookupFail;
    INTERPCELL cell;
    double logT1, T2;
    double S = SCVHEOS_FAIL_VALUE;
    double u2 = SCVHEOS_FAIL_VALUE;

    /* s(rho1, u1) in the cell of (rho1, u1), as scvheosSofRhoU: T fallback and density check in one lookup. */
    if (scvheosLocateRhoU(Mat, __func__, rho1, u1, &cell, &logT1) == 0) {
        /*
         * u1 not inverted (counted, the only count of a successful locate): the locate fell back
         * to T_min. Above the column (log u1 > log u(rho1, T_max)) the isentrope starts from the
         * top row instead. The strict inversion fails only outside the end nodes, so this is
         * above both of them. u1 <= 0 or NaN gives a NaN logu and stays at T_min.
         */
        if (Mat->nLookupFail > nLookupFail0 &&
            log10(u1) > interpLineValue(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, cell.ix, cell.u, Mat->nT - 1)) {
            cell.iy = Mat->nT - 2;
            cell.v = 1.0;
        }
        S = scvheosExp10(interpCellValue(&Mat->grid, Mat->dLogSArray, &cell));
    }
    if (!scvheosFailed(S)) {
        T2 = scvheosTofRhoS(Mat, rho2, S);
        if (!scvheosFailed(T2)) {
            u2 = scvheosUofRhoT(Mat, rho2, T2);
            if (scvheosFailed(u2)) u2 = SCVHEOS_FAIL_VALUE;
        }
    }

    if (Mat->nLookupFail > nLookupFail0) Mat->nLookupFail = nLookupFail0 + 1;
    return u2;
}

/*
 * The cold curve of the table: u(rho, T_min), evaluated on the first T row (no 10^x round trip of
 * T_min).
 */
double scvheosUCold(SCVHEOSMAT *Mat, double rho) {
    double logu = scvheosLogUofLogRhoLogT(Mat, log10(rho), Mat->LogTMin);
    if (scvheosFailed(logu)) return SCVHEOS_FAIL_VALUE;
    return scvheosExp10(logu);
}

int scvheosPrintMat(SCVHEOSMAT *Mat, FILE *fp) {
    if (!fp) return SCVHEOS_FAIL;
	fprintf(fp,"# Material: %i (%s)\n", Mat->iMat, Mat->MatString);
	fprintf(fp,"# Reference density rho0 (cgs units): %g\n", Mat->rho0*Mat->dGmPerCcUnit);
	fprintf(fp,"# Table size: nRho = %d, nT = %d\n", Mat->nRho, Mat->nT);
    fprintf(fp,"# Table boundaries (cgs units):\n");
    fprintf(fp,"# LogRhoMin=%15.7E LogRhoMax=%15.7E\n", Mat->dLogRhoAxis[0]+log10(Mat->dGmPerCcUnit), Mat->dLogRhoAxis[Mat->nRho-1]+log10(Mat->dGmPerCcUnit));
    fprintf(fp,"# LogTMin=%15.7E LogTMax=%15.7E\n", Mat->dLogTAxis[0], Mat->dLogTAxis[Mat->nT-1]);
    fprintf(fp,"# Failed lookups and fallbacks (nLookupFail): %ld\n", Mat->nLookupFail);

    return SCVHEOS_SUCCESS;
}

/*
 * Check if (logrho, logT) is inside of the eos table.
 */
int scvheosCheckTableBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    if ((logrho >= Mat->dLogRhoAxis[0]) && (logrho <= Mat->dLogRhoAxis[Mat->nRho-1]) && (logT >= Mat->dLogTAxis[0]) &&
        (logT <= Mat->dLogTAxis[Mat->nT-1])) {
        return TRUE;
    }

    return FALSE;
}

/*
 * Check if LogRhoMin <= logrho <= LogRhoMax and LogTMin <= logT <= LogTMax.
 */
int scvheosCheckBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    if ((logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
        return FALSE;
    }

    return TRUE;
}

/*
 * Check if a value (rho, u) is inside of the EOS table.
 *
 * Returns SCVHEOS_SUCCESS, or one of SCVHEOS_OUTSIDE_RHOMIN, _RHOMAX, _TMIN,
 * _TMAX. The energy bounds are u(rho, T_min) and u(rho, T_max) of the table
 * column, which assumes u(T) is monotonic at fixed rho (the same assumption
 * scvheosLogTofLogRhoLogU makes). The bounds are compared with the same
 * roundoff tolerance the inversion applies to its target, so a state that
 * passes this test can be inverted.
 *
 * rho <= 0 reports RHOMIN and u <= 0 reports TMIN: log10 of a non-positive
 * number is NaN or -inf, and the comparisons below would otherwise let the
 * state through. Non-finite input is reported the same way.
 */
int scvheosIsInExtrapLimit(SCVHEOSMAT *Mat, double rho, double u) {
    double logrho;
    double logu;

    if (!(rho > 0.0) || !isfinite(rho)) return SCVHEOS_OUTSIDE_RHOMIN;
    if (!(u > 0.0) || !isfinite(u)) return SCVHEOS_OUTSIDE_TMIN;

    logrho = log10(rho);
    logu = log10(u);

    /* Same roundoff tolerance at the edges as the lookups themselves. */
    {
        int ix = 0; double w = 0.0;
        int iret = interpLocateAxisTol(Mat->dLogRhoAxis, Mat->nRho, logrho, SCVHEOS_EDGE_TOL, SCVHEOS_EDGE_TOL_ABS, &ix, &w);
        if (iret == INTERP_BELOW) return SCVHEOS_OUTSIDE_RHOMIN;
        if (iret != INTERP_SUCCESS) return SCVHEOS_OUTSIDE_RHOMAX;
        {
            const double umin = interpLineValue(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, ix, w, 0);
            const double umax = interpLineValue(&Mat->grid, Mat->dLogUArray, INTERP_ALONG_Y, ix, w, Mat->nT - 1);
            const double eps = SCVHEOS_EDGE_TOL*(1.0 + fabs(logu));
            if (logu < umin - eps) return SCVHEOS_OUTSIDE_TMIN;
            if (logu > umax + eps) return SCVHEOS_OUTSIDE_TMAX;
        }
    }

    return SCVHEOS_SUCCESS;
}