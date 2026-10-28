/*
 * Header file for scvheos.c.
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 26.10.2026
 */
#ifndef SCVHEOS_HINCLUDED
#define SCVHEOS_HINCLUDED

#include <stdio.h>   /* FILE in scvheosPrintMat */
#include "interpBilinear/interpBilinear.h"   /* git subtree of the interpBilinear repository */

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
 *   1.2.0  interpolation through the header-only interpBilinear core (git
 *          subtree in interpBilinear/): the four forward lookups, the
 *          inversions, the partials, the sound-speed derivatives and the
 *          monotonicity report all use it; GSL is gone entirely (no
 *          gsl_interp2d, no accelerators, no GSL headers in this header).
 *          The edge tolerance is SCVHEOS_EDGE_TOL relative to the table
 *          bound instead of 1e-13*(1 + |x|); SCVHEOSMAT is read-only after
 *          initialisation apart from nLookupFail. Same inversion (strict
 *          mode) and same failure contract as 1.1.
 *   1.3.0  EOSlib failure contract (see SCVHEOS_FAIL_VALUE): lookups never
 *          exit, SCVH_NOEXIT and SCVHEOS_VERBOSE are gone, messages only
 *          with -DEOSLIB_VERBOSE (SCVHEOS_WARN), every failure counted in
 *          nLookupFail, T(rho, u) falls back to T_min like ANEOSTofRhoU,
 *          scvheosFailed() in the header; the log-space derivatives no
 *          longer exit outside the table. scvheosGenerateSoundSpeedTable checks
 *          Mat == NULL before the first dereference. An inconsistent unit pair
 *          (exactly one of dKpcUnit, dMsolUnit <= 0, NaN, inf) is fatal; it
 *          silently meant cgs. Requires interpBilinear 1.0.3 with
 *          INTERP_ABI_VERSION 1 (SCVHEOS_INTERP_*, prefixed: the unprefixed
 *          INTERP_*_REQUIRED clashed with the other vendoring libraries).
 *          New scvheosUCold, so that EOSlib no longer reads LogTMin.
 *          scvheosPrintMat prints nLookupFail.
 *          New scvheosSofRhoU (one cell lookup, like scvheosPCTofRhoU): EOSlib
 *          built s(rho, u) from T(rho, u) and s(rho, T), which counted a density
 *          outside the table twice. scvheosIsentropicU and the (rho, u) derivatives
 *          count a call at most once in nLookupFail (the isentrope counted rho1
 *          outside the table twice, the derivatives a T fallback into a column with
 *          dlogu/dlogT <= 0 twice).
 *          Edge tolerance on the log axes SCVHEOS_EDGE_TOL*|bound| +
 *          SCVHEOS_EDGE_TOL_ABS (a log bound at 0 had none); uses the *Tol
 *          functions of interpBilinear (to be required when its version is bumped).
 *          scvheosReadTable rejects nRho < 2 or nT < 2 (no cell), like reos3ReadTable.
 *          scvheosPrintMat prints rho0 in cgs with the unit, like the table boundaries and
 *          reos3PrintMat (it printed the code-unit value unlabelled).
 *          Documentation of the inversions corrected: a non-monotonic column or row does not
 *          fail (the strict inversion of interpBilinear cannot detect it), it gives a root that
 *          need not be the physical one; only a target outside the end nodes falls back, clamps
 *          or fails. No change in behaviour.
 *          scvheosIsentropicU starts an isentrope from (rho1, u1) above the u column at
 *          s(rho1, T_max), the top row, like the clamp at rho2; it started from s(rho1, T_min)
 *          (the T(rho, u) fallback), so u2 was the cold curve at rho2 and the thermal energy
 *          was lost. Below the cold curve it still starts from s(rho1, T_min).
 *          No output at initialisation unless compiled with -DEOSLIB_VERBOSE (SCVHEOS_WARN): the
 *          monotonicity report and the sound speed repairs were printed on stderr unconditionally,
 *          with one instance per thread 64 times. Errors before a fatal exit are still printed.
 *          scvheosCheckTableBoundsLogRhoLogT and scvheosCheckBoundsLogRhoLogT use the edge
 *          tolerance of the lookups (they are TRUE exactly where scvheosLocateLog succeeds); they
 *          compared without a tolerance and rejected states on the edge that every lookup accepts.
 *          The note that SCVHEOS_RHO0_CGS lies outside the table (rho0 moved to the centre of the rho
 *          axis) is SCVHEOS_WARN as well; it was the last unconditional fprintf of the init that is
 *          not followed by a fatal exit.
 */
#define SCVHEOS_VERSION_TEXT      "1.3.0"
#define SCVHEOS_VERSION_MAJOR     1
#define SCVHEOS_VERSION_MINOR     3
#define SCVHEOS_VERSION_PATCH     0

/*
 * interpBilinear: the minimum version and the layout of its types this library needs. The names
 * carry the library prefix because ANEOSmaterial, scvheos and reos3 all vendor the core and their
 * headers meet in one translation unit (EOSlib), where only the first copy of the core is
 * included (see INTERP_ABI_VERSION in interpBilinear.h). Checked at compile time.
 */
#define SCVHEOS_INTERP_VERSION_MAJOR_REQUIRED 1
#define SCVHEOS_INTERP_VERSION_MINOR_REQUIRED 0
#define SCVHEOS_INTERP_VERSION_PATCH_REQUIRED 3
#define SCVHEOS_INTERP_ABI_VERSION 1
#if (INTERP_VERSION_MAJOR * 10000 + INTERP_VERSION_MINOR * 100 + INTERP_VERSION_PATCH) < \
    (SCVHEOS_INTERP_VERSION_MAJOR_REQUIRED * 10000 + SCVHEOS_INTERP_VERSION_MINOR_REQUIRED * 100 + SCVHEOS_INTERP_VERSION_PATCH_REQUIRED)
#error "scvheos: the interpBilinear in effect is older than the required version (in EOSlib: the copy of the first library header included)"
#endif
#if !defined(INTERP_ABI_VERSION) || (INTERP_ABI_VERSION != SCVHEOS_INTERP_ABI_VERSION)
#error "scvheos: the interpBilinear in effect has another type layout (INTERP_ABI_VERSION) than this library; update the vendored copies together"
#endif

/*
 * Roundoff tolerance at the table edges, relative to the bound: a state on
 * the edge converted through 10^x and log10 lands within a few ulps of the
 * bound (in which direction depends on the libm), so it is accepted and
 * clamped onto the bound rather than rejected. The same tolerance is used
 * on the targets of the inversions (u, s, P in log space).
 *
 * The axes are log10 rho and log10 T, where the round-trip error of
 * log10(10^x) is absolute (~1e-16) as well as relative to |x|. The band at a
 * bound b is therefore SCVHEOS_EDGE_TOL*|b| + SCVHEOS_EDGE_TOL_ABS (the *Tol
 * functions of interpBilinear): with the relative part alone a bound at or
 * near 0, which the code-unit shift of log rho can produce, had no tolerance.
 */
#define SCVHEOS_EDGE_TOL 1e-13
#define SCVHEOS_EDGE_TOL_ABS 1e-13

/*
 * Status codes. SCVHEOS_FAIL is a generic failure (initialization, I/O);
 * it must not collide with the SCVHEOS_OUTSIDE_* codes below, which
 * scvheosIsInExtrapLimit returns and EOSlib maps one by one.
 */
#define SCVHEOS_SUCCESS 0
#define SCVHEOS_FAIL  -99

/*
 * Failure contract (the EOSlib convention, shared by all EOS libraries).
 *
 * A failed lookup never exits. The EOS functions return SCVHEOS_FAIL_VALUE
 * when the state is not positive and finite or outside the table, when an
 * inversion has no solution (rho(P, T) outside the row) or when a
 * derivative is undefined; it is returned both by the log-space functions
 * (in log space) and by the linear-space wrappers (rho, u, T, P, s, cs), and
 * scvheosPCTofRhoU sets *c to it as well. Two inversions fall back to a
 * table edge instead:
 *  - T(rho, u) never returns the sentinel (same contract as ANEOSTofRhoU):
 *    rho outside the table, or u below the cold curve u(rho, T_min), above
 *    u(rho, T_max) or not positive gives T_min, the lowest tabulated
 *    temperature (scvheosPCTofRhoU sets *T to it too); at a valid rho all
 *    (rho, u) functions then evaluate at (rho, T_min);
 *  - T(rho, s) at a valid rho: s below s(rho, T_min) gives T_min, above
 *    s(rho, T_max) gives T_max (the ANEOSisentropicU convention for the
 *    isentrope).
 * On a non-monotonic column or row (a table defect, reported at
 * initialisation with -DEOSLIB_VERBOSE) the inversions T(rho, u), T(rho, s) and rho(P, T) return a
 * root, not necessarily the physical one; they cannot detect it.
 *
 * Every failure and every fallback is counted in Mat->nLookupFail, once per
 * call (a call that fails or falls back counts 1, however many of its steps
 * do), and reported on stderr only when the library is compiled with -DEOSLIB_VERBOSE
 * (SCVHEOS_WARN). Test a returned value with scvheosFailed(x) (true for the
 * sentinel and for NaN, in log or linear space); EOSlib exposes the same
 * test as EOSFailed(). Errors during initialisation (missing or defective
 * table, unknown material) are fatal.
 */
#define SCVHEOS_FAIL_VALUE (-1e50)

static inline int scvheosFailed(double x) {
    return !(x > -1e40);
}

/*
 * Diagnostic messages of the lookups, only with -DEOSLIB_VERBOSE. The
 * silent variant still type-checks the arguments but compiles to nothing.
 */
#ifdef EOSLIB_VERBOSE
#define SCVHEOS_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define SCVHEOS_WARN(...) do { if (0) fprintf(stderr, __VA_ARGS__); } while (0)
#endif

/*
 * Reference density rho0 (cgs). A gas has no natural reference density;
 * this value is a convention and is used for
 *  - ballic: the density the model starts from,
 *  - EOSSolveBC: the lower bracket for rho2 starts at 0.999*rho0 and is
 *    shrunk by 0.9 at most EOS_SOLVEBC_MAX_EXPAND times, so interface
 *    densities down to ~3e-8 g/cc are reachable from this value,
 *  - EOSWoolfsonCoeff (obsolete): the cutoff below which no interface
 *    correction is done, 0.8*rho0 = 8e-4 g/cc. EOSlib applies it to the
 *    gases too; rho0forInterfaceCorrection (0 for SCvH) is not read.
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
 * After scvheosInitMaterial a SCVHEOSMAT is read-only except for the
 * nLookupFail counter, which the clamping lookups increment without any
 * locking. pkdgrav3 holds one material per thread anyway.
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

    /* The EOS tables are base-10 logarithms, laid out A[j*nT + i] (j: rho, i: T). */
    double *dLogRhoAxis;
    double *dLogTAxis;
    INTERPGRID grid;   /* {nRho, nT, dLogRhoAxis, dLogTAxis}: the view the interpBilinear core works on */
    double *dLogUArray;
    double *dLogPArray;
    double *dLogSArray;
    double *dLogCArray;

    /*
     * Number of failed lookups and edge fallbacks since initialisation (see
     * SCVHEOS_FAIL_VALUE). Same role as nLookupFail in ANEOSmaterial. Only
     * incremented by the library, never reset.
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
/* s(rho, u) in cgs from one cell lookup; SCVHEOS_FAIL_VALUE if rho is outside the table,
 * s(rho, T_min) if u cannot be inverted (counted once). */
double scvheosSofRhoU(SCVHEOSMAT *Mat, double rho, double u);
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

/* u(rho, T_min), the cold curve of the table (code units); SCVHEOS_FAIL_VALUE (counted) outside the table. */
double scvheosUCold(SCVHEOSMAT *Mat, double rho);

int scvheosPrintMat(SCVHEOSMAT *Mat, FILE *fp);

/*
 * TRUE if (logrho, logT) is inside the table, exactly where the (logrho, logT) lookups find it (the
 * same edge tolerance SCVHEOS_EDGE_TOL), FALSE otherwise or for NaN; not counted. The lookup limits
 * LogRhoMin..LogTMax are the table bounds, so both functions are the same test.
 */
int scvheosCheckTableBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
int scvheosCheckBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
int scvheosIsInExtrapLimit(SCVHEOSMAT *Mat, double rho, double u);
#endif