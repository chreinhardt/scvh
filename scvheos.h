/*
 * Header file for scvheos.c.
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 23.10.2026
 */
#ifndef SCVHEOS_HINCLUDED
#define SCVHEOS_HINCLUDED

#include <stdio.h>   /* FILE in scvheosPrintMat */
#include <gsl/gsl_interp2d.h>

/*
 * Version history
 *   1.1.0  failure sentinel SCVHEOS_FAIL_VALUE in linear space; live T_max
 *          check and u <= 0 rejection in scvheosIsInExtrapLimit; cs^2 floor;
 *          entropy in cgs (EOSlib convention); isentrope clamped at the table
 *          edge like ANEOSisentropicU, counted in nLookupFail; table
 *          validation (finite, rectangular, increasing axes, monotonicity
 *          report); direct closed-form inversion of the bilinear interpolant
 *          replaces the GSL root finders; scvheosPCTofRhoU; analytic
 *          derivatives on the bilinear patch; sound-speed table from central
 *          differences of the arrays; one accelerator pair; removed
 *          dLogBase, dGasConst, dSecUnit, InterpType, the *_GSL_Params
 *          structs and *_GSL_rootfinder callbacks; SCVHEOS_FAIL = -99;
 *          SCVHEOS_HHE_EXT for id 114 (old name kept as alias).
 */
#define SCVHEOS_VERSION_TEXT      "1.1.0"
#define SCVHEOS_VERSION_MAJOR     1
#define SCVHEOS_VERSION_MINOR     1
#define SCVHEOS_VERSION_PATCH     0

/*
 * Status codes. SCVHEOS_FAIL is a generic failure (initialization, I/O);
 * it must not collide with the SCVHEOS_OUTSIDE_* codes below, which
 * scvheosIsInExtrapLimit returns and EOSlib maps one by one.
 */
#define SCVHEOS_SUCCESS 0
#define SCVHEOS_FAIL  -99

/*
 * Value returned by the EOS functions when a lookup or an inversion fails
 * and the library is compiled with SCVH_NOEXIT (without it the functions
 * exit). It is returned both by the log-space functions (in log space) and
 * by the linear-space wrappers (rho, u, T, P, s, cs), so callers can test
 * for x < -1e40 in either case. A failed inversion (target outside the
 * table column, or a non-monotonic column) is treated the same way as any
 * other failed lookup.
 */
#define SCVHEOS_FAIL_VALUE (-1e50)

/*
 * Reference density rho0 (cgs). A gas has no natural reference density;
 * this value is a convention and is used for
 *  - ballic: the density the model starts from,
 *  - EOSSolveBC: the lower bracket for rho2 starts at 0.999*rho0 and is
 *    shrunk by 0.9 at most EOS_SOLVEBC_MAX_EXPAND times, so interface
 *    densities down to ~3e-8 g/cc are reachable from this value,
 *  - EOSWoolfsonCoeff (obsolete): the 0.8*rho0 cutoff.
 * scvheosInitMaterial moves it to the geometric centre of the rho axis if
 * it lies outside the table, since EOSSolveBC evaluates the EOS there.
 */
#define SCVHEOS_RHO0_CGS 1e-3

/*
 * Material ids (EOSlib reserves 110-130 for SCvH).
 *
 * 114 loads the full extended H-He table (460 x 76), not the low-rho/T cut
 * it was once limited to; SCVHEOS_HHE_EXT is the accurate name and
 * SCVHEOS_HHE_EXT_LOWRHOT is kept as a synonym for existing parameter files.
 * 115 is reserved but has no table yet.
 */
#define SCVHEOS_H               110
#define SCVHEOS_HE              111
#define SCVHEOS_HHE             112
#define SCVHEOS_HHE_LOWRHOT     113
#define SCVHEOS_HHE_EXT         114
#define SCVHEOS_HHE_EXT_LOWRHOT SCVHEOS_HHE_EXT   /* legacy name */
#define SCVHEOS_H_EXT           115               /* reserved, no table */
#define SCVHEOS_HE_EXT          116
#define SCVHEOS_HHE_H2O_Z_1     117
#define SCVHEOS_HHE_H2O_Z_05    118
#define SCVHEOS_HHE_H2O_Z_012   119


/*
 * Return codes of scvheosIsInExtrapLimit (besides SCVHEOS_SUCCESS).
 */
#define SCVHEOS_OUTSIDE_RHOMIN -1
#define SCVHEOS_OUTSIDE_RHOMAX -2
#define SCVHEOS_OUTSIDE_TMIN   -3
#define SCVHEOS_OUTSIDE_TMAX   -4

/*
 * TRUE and FALSE. Other EOS headers (reos3.h) define them too, so only
 * define them if nobody else has; the values must be 1 and 0 everywhere.
 */
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

/*
 * One SCVHEOSMAT is NOT safe to share between threads: the gsl_interp_accel
 * caches are mutated by every lookup. pkdgrav3 holds one material per
 * thread, which is what this library assumes.
 */
typedef struct scvheosMat {
    int iMat;
    int nRho;
    int nT;
    char MatString[256];
    double rho0;

    /* Lookup limits: the table bounds, no extrapolation. */
    double LogRhoMin;
    double LogRhoMax;
    double LogTMin;
    double LogTMax;

    /* Convert units to code units. */
    double dKpcUnit;
    double dMsolUnit;
    double dErgPerGmUnit;
    double dGmPerCcUnit;

    /* GSL 2D interpolators (bilinear; see scvheosInitMaterial). */
    gsl_interp2d *InterpLogP;
    gsl_interp2d *InterpLogU;
    gsl_interp2d *InterpLogS;
    gsl_interp2d *InterpLogCs;

    /* Accelerator, shared by the four interpolators (same axes). */
    gsl_interp_accel *xAcc;
    gsl_interp_accel *yAcc;

    /* The EOS tables are base-10 logarithms, laid out A[j*nT + i] (j: rho, i: T). */
    double *dLogRhoAxis;
    double *dLogTAxis;
    double *dLogUArray;
    double *dLogPArray;
    double *dLogSArray;
    double *dLogCArray;

    /*
     * Number of lookups that had to fall back to a table edge (currently:
     * scvheosLogTofLogRhoLogS clamping to T_min/T_max). Same role as
     * nLookupFail in ANEOSmaterial. Only incremented, never reset.
     */
    long nLookupFail;
} SCVHEOSMAT;

/* Functions to initialize and finalize materials. */
SCVHEOSMAT *scvheosInitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat);
int scvheosReadTable(SCVHEOSMAT *Mat, const char *chInFile, int nRho, int nT, int nSkip);
int scvheosGenerateSoundSpeedTable(SCVHEOSMAT *Mat);

/* Calculate EOS values. */
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosLogCsofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosSofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u);
/* P, cs and T from one cell lookup; the hydro-code entry point. */
double scvheosPCTofRhoU(SCVHEOSMAT *Mat, double rho, double u, double *c, double *T);

/* Calculate derivatives. */
double scvheosdLogPdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosdLogPdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosdLogUdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosdLogUdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosdLogSdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosdLogSdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosLogTofLogRhoLogU(SCVHEOSMAT *Mat, double logrho, double logu);
double scvheosLogTofLogRhoLogS(SCVHEOSMAT *Mat, double logrho, double logs);
double scvheosLogRhoofLogPLogT(SCVHEOSMAT *Mat, double logP, double logT);

double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosTofRhoS(SCVHEOSMAT *Mat, double rho, double s);
double scvheosRhoofPT(SCVHEOSMAT *Mat, double P, double T);

/* Derivatives on the bilinear patch: df/dx = f/x * dlogf/dlogx. */
double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdSdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdSdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdPdRhoofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosdPdUofRhoU(SCVHEOSMAT *Mat, double rho, double u);

double scvheosIsentropicU(SCVHEOSMAT *Mat, double rho1, double u1, double rho2);

int scvheosPrintMat(SCVHEOSMAT *Mat, FILE *fp);

int scvheosCheckTableBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
int scvheosCheckBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
int scvheosIsInExtrapLimit(SCVHEOSMAT *Mat, double rho, double u);
#endif