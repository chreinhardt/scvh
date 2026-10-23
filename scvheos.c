/*
 * Copyright (c) 2020 Christian Reinhardt
 *
 * This file provides an interface for the SCvH equation of state (Saumon et al. 1995) for Hydrogen
 * and Helium that was extended to lower temperatures by Vazan et al. (2013).
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 23.10.2026
 */

/* getline() is POSIX.1-2008; declare it under -std=c99 too. */
#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp2d.h>
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
 * Range check with roundoff tolerance at the ends. A state on the table edge
 * converted through 10^x and log10 lands within a few ulps of the bound (and
 * the edge-roundoff direction depends on the libm), so it is accepted and
 * clamped onto the bound rather than rejected. NaN and +-inf fail. Returns 1
 * if *px lies in [lo, hi] after clamping, 0 otherwise.
 */
static int scvheosInRange(double *px, double lo, double hi) {
    const double eps = 1e-13*(1.0 + fabs(*px));
    if (!(*px >= lo)) { if (*px >= lo - eps) *px = lo; else return 0; }
    if (!(*px <= hi)) { if (*px <= hi + eps) *px = hi; else return 0; }
    return 1;
}

/*
 * Lookup failure handling.
 *
 * The log-space functions return SCVHEOS_FAIL_VALUE when a lookup fails and
 * SCVH_NOEXIT is defined, otherwise they exit. The linear-space wrappers must
 * test the log value BEFORE exponentiating, because pow(10, -1e50) == 0 is a
 * perfectly valid-looking number and the failure would be lost.
 *
 * scvheosFailed(x) is true for the sentinel and for NaN, in log or linear space.
 */
static int scvheosFailed(double x) {
    return !(x > -1e40);
}

#ifdef SCVH_NOEXIT
#define SCVHEOS_FAIL_OR_EXIT() return SCVHEOS_FAIL_VALUE
#else
#define SCVHEOS_FAIL_OR_EXIT() exit(1)
#endif

/*
 * Non-fatal warnings for lookups that fall back to a table edge. Off by
 * default: a gas particle expanding into vacuum hits the T_min edge at
 * every step, and a message per particle per step would flood stderr.
 * Mat->nLookupFail counts them regardless; compile with -DSCVHEOS_VERBOSE
 * to also see each one.
 */
#ifdef SCVHEOS_VERBOSE
#define SCVHEOS_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define SCVHEOS_WARN(...) ((void)0)
#endif

/*
 * Direct inversion of the bilinear interpolant.
 *
 * The tables are bilinear in (logT, logrho). At fixed logrho (between the
 * rho nodes j and j+1, weight w) the interpolant is PIECEWISE LINEAR in
 * logT with node values
 *
 *     f_i = (1-w) A[j][i] + w A[j+1][i],   i = 0 .. nT-1,
 *
 * and likewise at fixed logT it is piecewise linear in logrho. So f = target
 * is solved exactly: locate the segment [i, i+1] with f_i <= target <= f_i+1
 * by bisection on the node index (log2(n) steps, two array reads each) and
 * invert the linear segment in closed form. This replaces the GSL root
 * finder (allocation + ~37 interpolations per call) and has no tolerance.
 *
 * Both require the same thing the root finder assumed: f_i increases with
 * i along the searched axis. Where a column is not monotonic the bisection
 * lands on a segment that does not contain the target and that is reported.
 */

/*
 * Locate x in a strictly increasing axis of n >= 2 nodes. Returns j with
 * axis[j] <= x <= axis[j+1], 0 <= j <= n-2, and the weight
 * w = (x - axis[j]) / (axis[j+1] - axis[j]) in [0, 1]. At x == axis[n-1]
 * the last cell is returned with w = 1, the convention GSL uses, so the
 * two interpolants agree on the whole closed interval.
 * The caller guarantees axis[0] <= x <= axis[n-1].
 */
static int scvheosFindCell(const double *axis, int n, double x, double *w) {
    int lo = 0;
    int hi = n - 1;

    while (hi - lo > 1) {
        int mid = lo + (hi - lo)/2;
        if (axis[mid] <= x) lo = mid; else hi = mid;
    }
    *w = (x - axis[lo])/(axis[lo+1] - axis[lo]);
    return lo;
}

/*
 * Solve f(x) = target along one axis of the table, where the node values
 * are f_i = (1-w) a0[i*stride] + w a1[i*stride] for the n nodes axis[0..n-1].
 *
 * For an inversion along T at fixed logrho: axis = dLogTAxis, n = nT,
 * a0 = Array + j*nT, a1 = a0 + nT, stride = 1 (consecutive T nodes).
 * For an inversion along rho at fixed logT: axis = dLogRhoAxis, n = nRho,
 * a0 = Array + i, a1 = a0 + 1, stride = nT (consecutive rho nodes).
 *
 * Returns  0 and *px on success, with the segment index *plo (0 <= lo <= n-2)
 *            and the position *pt in [0, 1] within it if the pointers are
 *            not NULL, so the caller can evaluate the other tables in the
 *            same cell without a second search,
 *         -1 if target < f_0       (below the table along this axis),
 *         +1 if target > f_{n-1}   (above the table along this axis),
 *         -2 if the located segment does not bracket the target, which
 *            happens only when f is not monotonically increasing.
 * On a flat segment (f_i == f_i+1 == target) the lower node is returned.
 */
static int scvheosInvertColumn(const double *axis, int n, const double *a0, const double *a1,
                               int stride, double w, double target, double *px,
                               int *plo, double *pt) {
    const double w0 = 1.0 - w;
    int lo = 0;
    int hi = n - 1;
    double flo = w0*a0[0] + w*a1[0];
    double fhi = w0*a0[(size_t)(n-1)*stride] + w*a1[(size_t)(n-1)*stride];
    double df, t;

    /*
     * A state on the table edge round-trips through 10^x and log10 with a
     * few ulps of error, so a target can land just outside [f_0, f_n-1]
     * (e.g. u(rho, T_max) fed back into T(rho, u)). Tolerate that at the
     * roundoff level, far below any physical scale, instead of failing.
     */
    {
        const double eps = 1e-13*(1.0 + fabs(target));
        if (!(target >= flo)) { if (target >= flo - eps) target = flo; else return -1; }
        if (!(target <= fhi)) { if (target <= fhi + eps) target = fhi; else return +1; }
    }

    while (hi - lo > 1) {
        int mid = lo + (hi - lo)/2;
        double fm = w0*a0[(size_t)mid*stride] + w*a1[(size_t)mid*stride];
        if (fm <= target) { lo = mid; flo = fm; } else { hi = mid; fhi = fm; }
    }

    /* On a monotonic column flo <= target <= fhi holds here by construction. */
    if (!(flo <= target) || !(target <= fhi)) return -2;

    df = fhi - flo;
    t = (df > 0.0) ? (target - flo)/df : 0.0;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    *px = axis[lo] + t*(axis[hi] - axis[lo]);
    if (plo) *plo = lo;
    if (pt) *pt = t;
    return 0;
}

/*
 * Bilinear value of a table in the cell (j, i) at weights w (rho) and v (T):
 * the same formula gsl_interp2d_bilinear evaluates, on the raw array.
 */
static double scvheosCellValue(const double *Array, int nT, int j, int i, double w, double v) {
    const double *c0 = Array + (size_t)j*nT + i;
    const double *c1 = c0 + nT;
    return (1.0 - w)*((1.0 - v)*c0[0] + v*c0[1]) + w*((1.0 - v)*c1[0] + v*c1[1]);
}

/*
 * Partial derivative of one table on the bilinear patch (j, i) at (w, t),
 * with respect to logrho (bWrtRho = TRUE) or logT (bWrtRho = FALSE).
 * On a bilinear patch the derivative along one axis is linear in the other
 * coordinate and jumps at the cell boundaries; it is the derivative
 * gsl_interp2d_eval_deriv_* would return, without the accelerator lookups.
 */
static double scvheosCellPartial(const SCVHEOSMAT *Mat, const double *Array, int j, int i,
                                 double w, double t, int bWrtRho) {
    const int nT = Mat->nT;
    const double *A0 = Array + (size_t)j*nT + i, *A1 = A0 + nT;

    if (bWrtRho) {
        return ((1.0 - t)*(A1[0] - A0[0]) + t*(A1[1] - A0[1]))/(Mat->dLogRhoAxis[j+1] - Mat->dLogRhoAxis[j]);
    }
    return ((1.0 - w)*(A0[1] - A0[0]) + w*(A1[1] - A1[0]))/(Mat->dLogTAxis[i+1] - Mat->dLogTAxis[i]);
}

/*
 * The four partials of logP and logu on the patch, for the (rho, u) derivatives.
 */
static void scvheosCellPartials(const SCVHEOSMAT *Mat, int j, int i, double w, double t,
                                double *dlPdlr, double *dlPdlT, double *dludlr, double *dludlT) {
    *dlPdlr = scvheosCellPartial(Mat, Mat->dLogPArray, j, i, w, t, TRUE);
    *dlPdlT = scvheosCellPartial(Mat, Mat->dLogPArray, j, i, w, t, FALSE);
    *dludlr = scvheosCellPartial(Mat, Mat->dLogUArray, j, i, w, t, TRUE);
    *dludlT = scvheosCellPartial(Mat, Mat->dLogUArray, j, i, w, t, FALSE);
}

/*
 * Check the monotonicity the inversions rely on, once per table at init:
 * u(T) and s(T) increasing along every rho column, P(rho) increasing along
 * every T row. Bilinear interpolation between two monotonic node columns is
 * monotonic, so checking the node columns covers every logrho. Reports the
 * number of defective columns/rows and the first few, so a later
 * "non-monotonic column" failure of an inversion can be traced to the table.
 */
static void scvheosCheckMonotonic(SCVHEOSMAT *Mat) {
    const int nT = Mat->nT, nRho = Mat->nRho, nReport = 5;
    int nBadU = 0, nBadS = 0, nBadP = 0;

    for (int j = 0; j < nRho; j++) {
        const double *u = Mat->dLogUArray + (size_t)j*nT;
        const double *s = Mat->dLogSArray + (size_t)j*nT;
        int badU = 0, badS = 0;
        for (int i = 1; i < nT; i++) {
            if (!(u[i] > u[i-1])) badU = i;
            if (!(s[i] > s[i-1])) badS = i;
        }
        if (badU) { if (nBadU < nReport) fprintf(stderr, "scvheosInitMaterial: iMat= %i: u(T) not increasing at jRho= %i (logrho= %g), e.g. iT= %i.\n", Mat->iMat, j, Mat->dLogRhoAxis[j], badU); nBadU++; }
        if (badS) { if (nBadS < nReport) fprintf(stderr, "scvheosInitMaterial: iMat= %i: s(T) not increasing at jRho= %i (logrho= %g), e.g. iT= %i.\n", Mat->iMat, j, Mat->dLogRhoAxis[j], badS); nBadS++; }
    }
    for (int i = 0; i < nT; i++) {
        int badP = 0;
        for (int j = 1; j < nRho; j++) {
            if (!(Mat->dLogPArray[(size_t)j*nT + i] > Mat->dLogPArray[(size_t)(j-1)*nT + i])) badP = j;
        }
        if (badP) { if (nBadP < nReport) fprintf(stderr, "scvheosInitMaterial: iMat= %i: P(rho) not increasing at iT= %i (logT= %g), e.g. jRho= %i.\n", Mat->iMat, i, Mat->dLogTAxis[i], badP); nBadP++; }
    }
    if (nBadU || nBadS || nBadP) {
        fprintf(stderr, "scvheosInitMaterial: iMat= %i: non-monotonic columns: u(T) %i of %i, s(T) %i of %i, P(rho) %i of %i; "
                "inversions fail (u, P) or clamp (s) where a target falls into such a column.\n",
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
            exit(1);
    }

    /*
     * Allocate memory and read the EOS table.
     */
    if (scvheosReadTable(Mat, inFile, nRho, nT, nSkip) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "scvheosInitMaterial: Could not read EOS table %s (see above).\n", inFile);
        scvheosFinalizeMaterial(Mat);
        exit(1);
    }

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

    /*
     * Initialize the GSL interpolator. Bilinear is not a choice: the direct
     * inversions (scvheosInvertColumn), the cell evaluation
     * (scvheosCellValue) and the analytic derivatives (scvheosCellPartials)
     * are written for the bilinear patch and would silently disagree with a
     * bicubic interpolant.
     */
    const gsl_interp2d_type *InterpType = gsl_interp2d_bilinear;

    /* One accelerator pair serves all four tables: they share the axes. */
    Mat->xAcc = gsl_interp_accel_alloc();
    Mat->yAcc = gsl_interp_accel_alloc();

    /* x corresponts to T and y corresponds to rho. */
    Mat->InterpLogP = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);
    Mat->InterpLogU = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);
    Mat->InterpLogS = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);
    Mat->InterpLogCs = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);

    gsl_interp2d_init(Mat->InterpLogP, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogPArray, Mat->nT, Mat->nRho);
    gsl_interp2d_init(Mat->InterpLogU, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogUArray, Mat->nT, Mat->nRho);
    gsl_interp2d_init(Mat->InterpLogS, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogSArray, Mat->nT, Mat->nRho);

    /* The sound speed table is built from dLogPArray/dLogUArray directly; it only needs the tables in code units. */
    if (scvheosGenerateSoundSpeedTable(Mat) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "scvheosInitMaterial: Could not generate table for sound speed.\n");
        scvheosFinalizeMaterial(Mat);
        exit(1);
    }

    gsl_interp2d_init(Mat->InterpLogCs, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogCArray, Mat->nT, Mat->nRho);

    scvheosCheckMonotonic(Mat);

    return Mat;
}

/*
 * Free memory. Safe to call on a partially initialized material (calloc
 * zeroes all pointers, and the GSL free functions accept NULL).
 */
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat) {
    if (Mat != NULL) {
        gsl_interp2d_free(Mat->InterpLogP);
        gsl_interp2d_free(Mat->InterpLogU);
        gsl_interp2d_free(Mat->InterpLogS);
        gsl_interp2d_free(Mat->InterpLogCs);

        gsl_interp_accel_free(Mat->xAcc);
        gsl_interp_accel_free(Mat->yAcc);

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

    /*
     * The table arrays stay allocated on failure: scvheosFinalizeMaterial
     * frees them. Only the local resources (fp, chLine) are released here,
     * on every path, via the cleanup label.
     */
    Mat->dLogRhoAxis = (double *) calloc(nRho, sizeof(double));
    Mat->dLogTAxis = (double *) calloc(nT, sizeof(double));

    /* The GSL interpolation functions require the tables are 1D arrays. */
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

    /* GSL requires strictly increasing axes. */
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
 * Derivative of f at node k of a strictly increasing axis x[0..n-1] by the
 * second-order central difference on a (possibly non-uniform) grid,
 *
 *   f'(x_k) = [h0^2 f_k+1 - h1^2 f_k-1 + (h1^2 - h0^2) f_k] / (h0 h1 (h0 + h1)),
 *   h0 = x_k - x_k-1,  h1 = x_k+1 - x_k,
 *
 * and the one-sided difference at the two ends. f is accessed with a stride.
 */
static double scvheosNodeDeriv(const double *x, int n, const double *f, int stride, int k) {
    if (k == 0)     return (f[stride] - f[0])/(x[1] - x[0]);
    if (k == n - 1) return (f[(size_t)k*stride] - f[(size_t)(k-1)*stride])/(x[k] - x[k-1]);
    {
        const double h0 = x[k] - x[k-1], h1 = x[k+1] - x[k];
        const double fm = f[(size_t)(k-1)*stride], f0 = f[(size_t)k*stride], fp = f[(size_t)(k+1)*stride];
        return (h0*h0*fp - h1*h1*fm + (h1*h1 - h0*h0)*f0)/(h0*h1*(h0 + h1));
    }
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
 * central differences of the log tables on the grid (second order, symmetric
 * about the node), taken directly from dLogPArray and dLogUArray:
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
    const int nT = Mat->nT, nRho = Mat->nRho;
    double *rhoNode, *TNode;
    int nReplaced = 0;
    const int nReport = 10;   /* print at most this many defective nodes */

    if (Mat == NULL)
        return SCVHEOS_FAIL;

    if ((Mat->dLogRhoAxis == NULL) || (Mat->dLogTAxis == NULL) || (Mat->dLogPArray == NULL) || (Mat->dLogUArray == NULL))
        return SCVHEOS_FAIL;

    if (Mat->dLogCArray != NULL)
        return SCVHEOS_FAIL;

    /* The GSL interpolation functions require the tables are 1D arrays. */
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
            const double dlPdlr = scvheosNodeDeriv(Mat->dLogRhoAxis, nRho, Mat->dLogPArray + i, nT, j);
            const double dlPdlT = scvheosNodeDeriv(Mat->dLogTAxis, nT, logP, 1, i);
            const double dludlT = scvheosNodeDeriv(Mat->dLogTAxis, nT, logu, 1, i);
            const double dPdrho = P/rho*dlPdlr;
            const double dPdT = P/T*dlPdlT;
            const double cv = u/T*dludlT;
            double cs2 = dPdrho + T/(rho*rho*cv)*dPdT*dPdT;

            /* Defective node: non-finite (cv = 0), or not positive. */
            if (!isfinite(cs2) || !(cs2 > 0.0)) {
                if (nReplaced < nReport) {
                    fprintf(stderr, "scvheosGenerateSoundSpeedTable: iMat= %i iT= %i jRho= %i "
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
        fprintf(stderr, "scvheosGenerateSoundSpeedTable: iMat= %i: cs^2 replaced by P/rho at %i of %i nodes.\n",
                Mat->iMat, nReplaced, nT*nRho);
    }

    return SCVHEOS_SUCCESS;
}

// Functions that have to be implemented or added
// We also need derivatives to calculate the sound speed (maybe do this once and make a table)?


/*
 * Calculate logP(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logP;

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax) || !scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogPofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogP, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogPArray, logT, logrho,
			    Mat->xAcc, Mat->yAcc, &logP) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogPofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
 
    return logP;
}

/*
 * Calculate logU(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logu;

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax) || !scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogUofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogU, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogUArray, logT, logrho,
			    Mat->xAcc, Mat->yAcc, &logu) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogUofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
    
    return logu;
}

/*
 * Calculate logs(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logs;

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax) || !scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogSofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogS, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogSArray, logT, logrho,
			    Mat->xAcc, Mat->yAcc, &logs) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogSofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
    
    return logs;
}

/*
 * Calculate logcs(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogCsofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logcs;

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax) || !scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogCsofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogCs, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogCArray, logT, logrho,
			    Mat->xAcc, Mat->yAcc, &logcs) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogCsofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
 
    return logcs;
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
 * Locate the state (rho, u) in the table: the rho cell j with weight w, and
 * the T cell i with position t from the inversion of the u column, plus logT.
 * Shared by scvheosPCTofRhoU and the (rho, u) derivatives.
 *
 * Returns 0 on success, -1 (after a message naming chFunc) if the input is
 * non-finite or non-positive, rho is outside the table, or u is outside the
 * column or in a non-monotonic column. The caller decides between the
 * sentinel and exit.
 */
static int scvheosLocateRhoU(SCVHEOSMAT *Mat, const char *chFunc, double rho, double u,
                             int *pj, int *pi, double *pw, double *pt, double *plogT) {
    const int nT = Mat->nT;
    const double *a0;
    double logrho, logu, w;
    int j, iret;

    if (!(rho > 0.0) || !(u > 0.0) || !isfinite(rho) || !isfinite(u)) {
        fprintf(stderr, "%s: invalid input rho= %g u= %g.\n", chFunc, rho, u);
        return -1;
    }

    logrho = log10(rho);
    logu = log10(u);

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax)) {
        fprintf(stderr, "%s: logrho= %15.7E outside of the EOS table.\n", chFunc, logrho);
        return -1;
    }

    j = scvheosFindCell(Mat->dLogRhoAxis, Mat->nRho, logrho, &w);
    a0 = Mat->dLogUArray + (size_t)j*nT;

    iret = scvheosInvertColumn(Mat->dLogTAxis, nT, a0, a0 + nT, 1, w, logu, plogT, pi, pt);

    if (iret != 0) {
        const double w0 = 1.0 - w;
        fprintf(stderr, "%s: logu= %g %s at logrho= %g "
                "(column u(T_min)= %g u(T_max)= %g).\n", chFunc, logu,
                (iret == -1) ? "below the table" : (iret == +1) ? "above the table" : "in a non-monotonic column",
                logrho, w0*a0[0] + w*a0[nT], w0*a0[nT-1] + w*a0[2*nT-1]);
        return -1;
    }

    *pj = j;
    *pw = w;
    return 0;
}

/*
 * Calculate the pressure P(rho, u), the sound speed cs(rho, u) and the
 * temperature T(rho, u) from one cell lookup.
 *
 * The inversion of the u column gives the T segment; P and cs are then the
 * bilinear values of their tables in that same cell, with no further search
 * and no detour through GSL. This is the hydro-code entry point (EOSPCTofRhoU)
 * and mirrors ANEOSPCTofRhoU.
 *
 * On failure (see scvheosLocateRhoU) returns SCVHEOS_FAIL_VALUE and sets *c
 * and *T to it as well (or exits without SCVH_NOEXIT).
 */
double scvheosPCTofRhoU(SCVHEOSMAT *Mat, double rho, double u, double *c, double *T) {
    double logT, w, t;
    int j, i;

    *c = SCVHEOS_FAIL_VALUE;
    *T = SCVHEOS_FAIL_VALUE;

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &j, &i, &w, &t, &logT) != 0) {
        SCVHEOS_FAIL_OR_EXIT();
    }

    *T = scvheosExp10(logT);
    *c = scvheosExp10(scvheosCellValue(Mat->dLogCArray, Mat->nT, j, i, w, t));
    return scvheosExp10(scvheosCellValue(Mat->dLogPArray, Mat->nT, j, i, w, t));
}

/*
 * Calculate the pressure P(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion or the lookup fails.
 */
double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double c, T;
    return scvheosPCTofRhoU(Mat, rho, u, &c, &T);
}

/*
 * Calculate the sound speed cs(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion or the lookup fails.
 */
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double c, T;
    scvheosPCTofRhoU(Mat, rho, u, &c, &T);
    return c;
}

/*
 * Calculate the temperature logT(logrho, logu) by direct inversion of the
 * u column at logrho (see scvheosInvertColumn).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if the input is
 * not finite, logrho is outside the table, logu is outside the column
 * [u(T_min), u(T_max)], or the column is not monotonic at logu.
 */
double scvheosLogTofLogRhoLogU(SCVHEOSMAT *Mat, double logrho, double logu) {
    const int nT = Mat->nT;
    const double *a0;
    double w;
    double logT;
    int j, iret;

    /* log10 of a non-positive rho or u is NaN/-inf: nothing to invert. */
    if (!isfinite(logrho) || !isfinite(logu)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogU: invalid input logrho= %g logu= %g.\n", logrho, logu);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogU: logrho= %15.7E outside of the EOS table.\n", logrho);
        SCVHEOS_FAIL_OR_EXIT();
    }

    j = scvheosFindCell(Mat->dLogRhoAxis, Mat->nRho, logrho, &w);
    a0 = Mat->dLogUArray + (size_t)j*nT;

    iret = scvheosInvertColumn(Mat->dLogTAxis, nT, a0, a0 + nT, 1, w, logu, &logT, NULL, NULL);

    if (iret != 0) {
        const double w0 = 1.0 - w;
        fprintf(stderr, "scvheosLogTofLogRhoLogU: logu= %g %s at logrho= %g "
                "(column u(T_min)= %g u(T_max)= %g, logT in [%g, %g]).\n", logu,
                (iret == -1) ? "below the table" : (iret == +1) ? "above the table" : "in a non-monotonic column",
                logrho, w0*a0[0] + w*a0[nT], w0*a0[nT-1] + w*a0[2*nT-1], Mat->LogTMin, Mat->LogTMax);
        SCVHEOS_FAIL_OR_EXIT();
    }


    return logT;
}

/*
 * Calculate the temperature logT(logrho, logs) by direct inversion of the
 * s column at logrho (see scvheosInvertColumn).
 *
 * logs below s(T_min) or above s(T_max) of the column returns logT_min or
 * logT_max respectively, counted in Mat->nLookupFail (ANEOSmaterial
 * convention). Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if
 * the input is not finite, logrho is outside the table, or the column is
 * not monotonic at logs.
 */
double scvheosLogTofLogRhoLogS(SCVHEOSMAT *Mat, double logrho, double logs) {
    const int nT = Mat->nT;
    const double *a0;
    double w;
    double logT;
    int j, iret;

    if (!isfinite(logrho) || !isfinite(logs)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogS: invalid input logrho= %g logs= %g.\n", logrho, logs);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogS: logrho= %15.7E outside of the EOS table.\n", logrho);
        SCVHEOS_FAIL_OR_EXIT();
    }

    j = scvheosFindCell(Mat->dLogRhoAxis, Mat->nRho, logrho, &w);
    a0 = Mat->dLogSArray + (size_t)j*nT;

    iret = scvheosInvertColumn(Mat->dLogTAxis, nT, a0, a0 + nT, 1, w, logs, &logT, NULL, NULL);

    /*
     * logs outside the column [s(T_min), s(T_max)]: clamp to the table edge,
     * the convention ANEOSisentropicU uses (cold curve u(rho2, T_min) when S
     * is not attained at rho2). Unlike there the fallback is counted in
     * Mat->nLookupFail and reported with SCVHEOS_WARN rather than silent.
     */
    if (iret == -1) {
        Mat->nLookupFail++;
        SCVHEOS_WARN("scvheosLogTofLogRhoLogS: logs= %g below s(T_min)= %g at logrho= %g, "
                     "falling back to logT_min= %g\n", logs, (1.0 - w)*a0[0] + w*a0[nT], logrho, Mat->LogTMin);
        return Mat->LogTMin;
    }

    if (iret == +1) {
        Mat->nLookupFail++;
        SCVHEOS_WARN("scvheosLogTofLogRhoLogS: logs= %g above s(T_max)= %g at logrho= %g, "
                     "falling back to logT_max= %g\n", logs, (1.0 - w)*a0[nT-1] + w*a0[2*nT-1], logrho, Mat->LogTMax);
        return Mat->LogTMax;
    }

    if (iret != 0) {
        fprintf(stderr, "scvheosLogTofLogRhoLogS: s column not monotonic at logs= %g logrho= %g.\n", logs, logrho);
        SCVHEOS_FAIL_OR_EXIT();
    }


    return logT;
}

/*
 * Calculate the density logrho(logP, logT) by direct inversion of the P row
 * at logT (see scvheosInvertColumn; the search runs along rho with stride nT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if the input is
 * not finite, logT is outside the table, logP is outside the row
 * [P(rho_min), P(rho_max)], or P(rho) is not monotonic at logP.
 */
double scvheosLogRhoofLogPLogT(SCVHEOSMAT *Mat, double logP, double logT) {
    const int nT = Mat->nT;
    const double *a0;
    double v;
    double logrho;
    int i, iret;

    if (!isfinite(logP) || !isfinite(logT)) {
        fprintf(stderr, "scvheosLogRhoofLogPLogT: invalid input logP= %g logT= %g.\n", logP, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (!scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
        fprintf(stderr, "scvheosLogRhoofLogPLogT: logT= %15.7E outside of the EOS table.\n", logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    i = scvheosFindCell(Mat->dLogTAxis, nT, logT, &v);
    a0 = Mat->dLogPArray + i;

    iret = scvheosInvertColumn(Mat->dLogRhoAxis, Mat->nRho, a0, a0 + 1, nT, v, logP, &logrho, NULL, NULL);

    if (iret != 0) {
        const double v0 = 1.0 - v;
        const size_t last = (size_t)(Mat->nRho - 1)*nT;
        fprintf(stderr, "scvheosLogRhoofLogPLogT: logP= %g %s at logT= %g "
                "(row P(rho_min)= %g P(rho_max)= %g, logrho in [%g, %g]).\n", logP,
                (iret == -1) ? "below the table" : (iret == +1) ? "above the table" : "in a non-monotonic row",
                logT, v0*a0[0] + v*a0[1], v0*a0[last] + v*a0[last+1], Mat->LogRhoMin, Mat->LogRhoMax);
        SCVHEOS_FAIL_OR_EXIT();
    }


    return logrho;
}

/*
 * Calculate the temperature T(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion fails.
 */
double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double logrho;
    double logu;
    double logT;

    logrho = log10(rho);
    logu = log10(u);

    logT = scvheosLogTofLogRhoLogU(Mat, logrho, logu);
    if (scvheosFailed(logT)) return SCVHEOS_FAIL_VALUE;

    return scvheosExp10(logT);
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
 * cell at the upper table edge), the convention scvheosFindCell and GSL share.
 */
static double scvheosdLogdLog(SCVHEOSMAT *Mat, const double *Array,
                              double logrho, double logT, int bWrtRho, const char *chFunc) {
    double w, t;
    int j, i;

    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax) || !scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
	    fprintf(stderr, "%s: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", chFunc, logrho, logT);
        exit(1);
    }

    j = scvheosFindCell(Mat->dLogRhoAxis, Mat->nRho, logrho, &w);
    i = scvheosFindCell(Mat->dLogTAxis, Mat->nT, logT, &t);

    return scvheosCellPartial(Mat, Array, j, i, w, t, bWrtRho);
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
 * -1 after a message if the state is not finite or outside the table.
 */
static int scvheosLocateRhoT(SCVHEOSMAT *Mat, const char *chFunc, double rho, double T,
                             int *pj, int *pi, double *pw, double *pt) {
    double logrho, logT;

    if (!(rho > 0.0) || !(T > 0.0) || !isfinite(rho) || !isfinite(T)) {
        fprintf(stderr, "%s: invalid input rho= %g T= %g.\n", chFunc, rho, T);
        return -1;
    }
    logrho = log10(rho);
    logT = log10(T);
    if (!scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax) || !scvheosInRange(&logT, Mat->LogTMin, Mat->LogTMax)) {
        fprintf(stderr, "%s: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", chFunc, logrho, logT);
        return -1;
    }
    *pj = scvheosFindCell(Mat->dLogRhoAxis, Mat->nRho, logrho, pw);
    *pi = scvheosFindCell(Mat->dLogTAxis, Mat->nT, logT, pt);
    return 0;
}

/*
 * Calculate the derivative dP/drho(rho, T) = (P/rho) dlogP/dlogrho on the
 * bilinear patch. Returns SCVHEOS_FAIL_VALUE (or exits) outside the table.
 */
double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double w, t;
    int j, i;

    if (scvheosLocateRhoT(Mat, __func__, rho, T, &j, &i, &w, &t) != 0) SCVHEOS_FAIL_OR_EXIT();
    return scvheosExp10(scvheosCellValue(Mat->dLogPArray, Mat->nT, j, i, w, t))/rho
           *scvheosCellPartial(Mat, Mat->dLogPArray, j, i, w, t, TRUE);
}

/*
 * Calculate the derivative dP/dT(rho, T) = (P/T) dlogP/dlogT.
 */
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double w, t;
    int j, i;

    if (scvheosLocateRhoT(Mat, __func__, rho, T, &j, &i, &w, &t) != 0) SCVHEOS_FAIL_OR_EXIT();
    return scvheosExp10(scvheosCellValue(Mat->dLogPArray, Mat->nT, j, i, w, t))/T
           *scvheosCellPartial(Mat, Mat->dLogPArray, j, i, w, t, FALSE);
}

/*
 * Calculate the derivative du/drho(rho, T) = (u/rho) dlogu/dlogrho.
 */
double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double w, t;
    int j, i;

    if (scvheosLocateRhoT(Mat, __func__, rho, T, &j, &i, &w, &t) != 0) SCVHEOS_FAIL_OR_EXIT();
    return scvheosExp10(scvheosCellValue(Mat->dLogUArray, Mat->nT, j, i, w, t))/rho
           *scvheosCellPartial(Mat, Mat->dLogUArray, j, i, w, t, TRUE);
}

/*
 * Calculate the derivative du/dT(rho, T) = (u/T) dlogu/dlogT, the specific
 * heat capacity at constant volume.
 */
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double w, t;
    int j, i;

    if (scvheosLocateRhoT(Mat, __func__, rho, T, &j, &i, &w, &t) != 0) SCVHEOS_FAIL_OR_EXIT();
    return scvheosExp10(scvheosCellValue(Mat->dLogUArray, Mat->nT, j, i, w, t))/T
           *scvheosCellPartial(Mat, Mat->dLogUArray, j, i, w, t, FALSE);
}

/*
 * Calculate the derivative ds/drho(rho, T) = (s/rho) dlogs/dlogrho (s in cgs).
 */
double scvheosdSdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double w, t;
    int j, i;

    if (scvheosLocateRhoT(Mat, __func__, rho, T, &j, &i, &w, &t) != 0) SCVHEOS_FAIL_OR_EXIT();
    return scvheosExp10(scvheosCellValue(Mat->dLogSArray, Mat->nT, j, i, w, t))/rho
           *scvheosCellPartial(Mat, Mat->dLogSArray, j, i, w, t, TRUE);
}

/*
 * Calculate the derivative ds/dT(rho, T) = (s/T) dlogs/dlogT (s in cgs).
 */
double scvheosdSdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double w, t;
    int j, i;

    if (scvheosLocateRhoT(Mat, __func__, rho, T, &j, &i, &w, &t) != 0) SCVHEOS_FAIL_OR_EXIT();
    return scvheosExp10(scvheosCellValue(Mat->dLogSArray, Mat->nT, j, i, w, t))/T
           *scvheosCellPartial(Mat, Mat->dLogSArray, j, i, w, t, FALSE);
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
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (rho, u)
 * cannot be located or du/dT|rho <= 0 in the cell (the heat capacity must be
 * positive; a non-positive value is a table artefact).
 */
double scvheosdPdRhoofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double logT, w, t, P;
    double dlPdlr, dlPdlT, dludlr, dludlT;
    int j, i;

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &j, &i, &w, &t, &logT) != 0) {
        SCVHEOS_FAIL_OR_EXIT();
    }

    scvheosCellPartials(Mat, j, i, w, t, &dlPdlr, &dlPdlT, &dludlr, &dludlT);

    if (!(dludlT > 0.0)) {
        fprintf(stderr, "scvheosdPdRhoofRhoU: dlogu/dlogT= %g <= 0 at rho= %g u= %g (cell jRho= %i iT= %i).\n",
                dludlT, rho, u, j, i);
        SCVHEOS_FAIL_OR_EXIT();
    }

    P = scvheosExp10(scvheosCellValue(Mat->dLogPArray, Mat->nT, j, i, w, t));

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
    double logT, w, t, P;
    double dlPdlr, dlPdlT, dludlr, dludlT;
    int j, i;

    if (scvheosLocateRhoU(Mat, __func__, rho, u, &j, &i, &w, &t, &logT) != 0) {
        SCVHEOS_FAIL_OR_EXIT();
    }

    scvheosCellPartials(Mat, j, i, w, t, &dlPdlr, &dlPdlT, &dludlr, &dludlT);

    if (!(dludlT > 0.0)) {
        fprintf(stderr, "scvheosdPdUofRhoU: dlogu/dlogT= %g <= 0 at rho= %g u= %g (cell jRho= %i iT= %i).\n",
                dludlT, rho, u, j, i);
        SCVHEOS_FAIL_OR_EXIT();
    }

    P = scvheosExp10(scvheosCellValue(Mat->dLogPArray, Mat->nT, j, i, w, t));

    return P/u*dlPdlT/dludlT;
}

/*
 * Calculate u2 so that (rho1, u1) and (rho2, u2) are on the same isentrope.
 *
 * If the entropy S(rho1, u1) is not attained at rho2 the state is clamped to
 * the table edge: u(rho2, T_min) (the cold curve) when it would fall below
 * the table, u(rho2, T_max) above it. This follows ANEOSisentropicU; the
 * fallback is counted in Mat->nLookupFail. Returns SCVHEOS_FAIL_VALUE if
 * rho1 or rho2 is outside the table or (rho1, u1) cannot be inverted.
 */
double scvheosIsentropicU(SCVHEOSMAT *Mat, double rho1, double u1, double rho2) {
    double T1;
    double T2;
    double S;
    double u2;

    T1 = scvheosTofRhoU(Mat, rho1, u1);
    if (scvheosFailed(T1)) return SCVHEOS_FAIL_VALUE;

    S = scvheosSofRhoT(Mat, rho1, T1);
    if (scvheosFailed(S)) return SCVHEOS_FAIL_VALUE;

    T2 = scvheosTofRhoS(Mat, rho2, S);
    if (scvheosFailed(T2)) return SCVHEOS_FAIL_VALUE;
   
    u2 = scvheosUofRhoT(Mat, rho2, T2);
    if (scvheosFailed(u2)) return SCVHEOS_FAIL_VALUE;

    return u2;
}


int scvheosPrintMat(SCVHEOSMAT *Mat, FILE *fp) {
    if (!fp) return SCVHEOS_FAIL;
	fprintf(fp,"# Material: %i (%s)\n", Mat->iMat, Mat->MatString);
	fprintf(fp,"# Reference density rho0: %g\n", Mat->rho0);
	fprintf(fp,"# Table size: nRho = %d, nT = %d\n", Mat->nRho, Mat->nT);
    fprintf(fp,"# Table boundaries (cgs units):\n");
    fprintf(fp,"# LogRhoMin=%15.7E LogRhoMax=%15.7E\n", Mat->dLogRhoAxis[0]+log10(Mat->dGmPerCcUnit), Mat->dLogRhoAxis[Mat->nRho-1]+log10(Mat->dGmPerCcUnit));
    fprintf(fp,"# LogTMin=%15.7E LogTMax=%15.7E\n", Mat->dLogTAxis[0], Mat->dLogTAxis[Mat->nT-1]);

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
 * scvheosLogTofLogRhoLogU makes). A state that passes this test can be
 * inverted: both bracket ends are inside the table and the root is bracketed.
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
    if (logrho < Mat->LogRhoMin && !scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax)) return SCVHEOS_OUTSIDE_RHOMIN;
    if (logrho > Mat->LogRhoMax && !scvheosInRange(&logrho, Mat->LogRhoMin, Mat->LogRhoMax)) return SCVHEOS_OUTSIDE_RHOMAX;
    {
        const double eps = 1e-13*(1.0 + fabs(logu));
        if (logu < scvheosLogUofLogRhoLogT(Mat, logrho, Mat->LogTMin) - eps) return SCVHEOS_OUTSIDE_TMIN;
        if (logu > scvheosLogUofLogRhoLogT(Mat, logrho, Mat->LogTMax) + eps) return SCVHEOS_OUTSIDE_TMAX;
    }

    return SCVHEOS_SUCCESS;
}