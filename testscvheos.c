/*  This file is part of scvheos.
 *  Copyright (c) 2020-2026 Thomas Meier & Christian Reinhardt
 *
 *  scvheos is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  scvheos is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with scvheos.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Test suite for the scvheos library.
 *
 * Every test is an invariant that must hold for any valid table, so no
 * reference values are hard-coded. Usage:
 *
 *   testscvheos [iMat] [dKpcUnit dMsolUnit]
 *
 * Defaults: iMat = SCVHEOS_HHE_EXT (114), Earth-like code units. Pass
 * negative units to run in cgs. The table file for iMat must be in the
 * working directory. Exit status is the number of failed checks.
 *
 * The tables are bilinear in (logT, logrho). The inversions are closed-form
 * on that interpolant, so round trips are tested to roundoff, not to a
 * solver tolerance. Non-monotonic columns (u(T), s(T)) or rows (P(rho)) are
 * table defects: they are reported by the first test and skipped where the
 * library's documented response (failure or clamp) is not an invariant.
 *
 * The failure contract (SCVHEOS_FAIL_VALUE instead of exit) only exists when
 * the library is built with -DSCVH_NOEXIT; those checks are compiled in only
 * then (see the Makefile target, which builds its own scvheos.noexit.o).
 * The library prints a diagnostic for every failed lookup; those messages are
 * expected in the last test.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "scvheos.h"

/* ---------------------------------------------------------------- harness */

static int nChecks = 0;
static int nFails = 0;
static const char *currentTest = "";

#define CHECK(cond, ...) do { \
    nChecks++; \
    if (!(cond)) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n"); \
    } \
} while (0)

#define CHECK_REL(a, b, tol, ...) do { \
    double _a = (a), _b = (b); \
    double _den = fmax(fabs(_a), fabs(_b)); \
    double _err = (_den > 0.0) ? fabs(_a - _b) / _den : fabs(_a - _b); \
    nChecks++; \
    if (!(_err <= (tol))) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, " (got %.15e vs %.15e, rel err %.3e > %.1e)\n", _a, _b, _err, (double)(tol)); \
    } \
} while (0)

/* Absolute tolerance, for quantities in log10 space. */
#define CHECK_ABS(a, b, tol, ...) do { \
    double _a = (a), _b = (b); \
    double _err = fabs(_a - _b); \
    nChecks++; \
    if (!(_err <= (tol))) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, " (got %.15e vs %.15e, abs err %.3e > %.1e)\n", _a, _b, _err, (double)(tol)); \
    } \
} while (0)

/* Relative tolerance with a floor: for derivatives that are legitimately zero
 * (an ideal-gas region has du/drho|T = 0), where a pure relative check has no
 * scale. The error is measured against max(|a|, |b|, scale). */
#define CHECK_NEAR(a, b, tol, scale, ...) do { \
    double _a = (a), _b = (b); \
    double _den = fmax(fmax(fabs(_a), fabs(_b)), (scale)); \
    double _err = fabs(_a - _b) / _den; \
    nChecks++; \
    if (!(_err <= (tol))) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, " (got %.15e vs %.15e, err %.3e > %.1e relative to %.3e)\n", _a, _b, _err, (double)(tol), _den); \
    } \
} while (0)

#define TEST(name) do { currentTest = name; printf("  %-44s", name); fflush(stdout); } while (0)
#define END_TEST() do { printf("%s\n", nFails == failsBefore ? "ok" : "FAILED"); } while (0)
/* Informational line inside a running test; re-indents so END_TEST aligns. */
#define NOTE(...) do { printf("\n    note: "); printf(__VA_ARGS__); printf("\n  %-44s", ""); } while (0)

static int isSentinel(double v) { return v < -1e40; }

/* Table access: A[j*nT + i], j along rho, i along T. */
#define LOGP(m, j, i) ((m)->dLogPArray[(size_t)(j) * (m)->nT + (i)])
#define LOGU(m, j, i) ((m)->dLogUArray[(size_t)(j) * (m)->nT + (i)])
#define LOGS(m, j, i) ((m)->dLogSArray[(size_t)(j) * (m)->nT + (i)])
#define LOGC(m, j, i) ((m)->dLogCArray[(size_t)(j) * (m)->nT + (i)])

static double tenpow(double x) { return exp(x * 2.302585092994045684); }

/* Sample points: interior of the table, log-spaced in rho and T. */
/* 100 for a quick run (seconds); 1000 for a full sweep (a minute). */
#define NSAMPLE_RHO 100
#define NSAMPLE_T   100
#define SAMPLE_OFFSET 0.0137   /* keeps log-spaced samples off the table nodes */

static double sampleLogRho(SCVHEOSMAT *m, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5 + SAMPLE_OFFSET) / NSAMPLE_RHO;
    return m->LogRhoMin + f * (m->LogRhoMax - m->LogRhoMin);
}

static double sampleLogT(SCVHEOSMAT *m, int k)
{
    double f = 0.05 + 0.90 * (k + 0.5 + SAMPLE_OFFSET) / NSAMPLE_T;
    return m->LogTMin + f * (m->LogTMax - m->LogTMin);
}

static double sampleRho(SCVHEOSMAT *m, int k) { return tenpow(sampleLogRho(m, k)); }
static double sampleT(SCVHEOSMAT *m, int k)   { return tenpow(sampleLogT(m, k)); }

/* Cell containing x on a strictly increasing axis: axis[j] <= x < axis[j+1]. */
static int findCell(const double *axis, int n, double x)
{
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) { int mid = (lo + hi) / 2; if (axis[mid] <= x) lo = mid; else hi = mid; }
    return lo;
}

/* 1 if f(T) is strictly increasing along the whole columns j and j+1 of a table. */
static int monotonicInColumn(SCVHEOSMAT *m, const double *A, int j)
{
    for (int i = 0; i < m->nT - 1; i++) {
        if (!(A[(size_t)j * m->nT + i + 1] > A[(size_t)j * m->nT + i])) return 0;
        if (!(A[(size_t)(j+1) * m->nT + i + 1] > A[(size_t)(j+1) * m->nT + i])) return 0;
    }
    return 1;
}

/* 1 if P(rho) is strictly increasing along the whole rows i and i+1. */
static int pressureMonotonicInRow(SCVHEOSMAT *m, int i)
{
    for (int j = 0; j < m->nRho - 1; j++) {
        if (!(LOGP(m, j+1, i) > LOGP(m, j, i))) return 0;
        if (!(LOGP(m, j+1, i+1) > LOGP(m, j, i+1))) return 0;
    }
    return 1;
}

/* The bilinear value of a table at (logrho, logT), computed independently. */
static double bilinear(SCVHEOSMAT *m, const double *A, double lr, double lt)
{
    int j = findCell(m->dLogRhoAxis, m->nRho, lr), i = findCell(m->dLogTAxis, m->nT, lt);
    double w = (lr - m->dLogRhoAxis[j]) / (m->dLogRhoAxis[j+1] - m->dLogRhoAxis[j]);
    double v = (lt - m->dLogTAxis[i]) / (m->dLogTAxis[i+1] - m->dLogTAxis[i]);
    const double *c0 = A + (size_t)j * m->nT + i, *c1 = c0 + m->nT;
    return (1.0 - w) * ((1.0 - v) * c0[0] + v * c0[1]) + w * ((1.0 - v) * c1[0] + v * c1[1]);
}

/* ------------------------------------------------------------------ tests */

static void testInit(SCVHEOSMAT *m, int iMat, double dKpcUnit, double dMsolUnit)
{
    int failsBefore = nFails;
    TEST("initialisation and table sanity");
    CHECK(m != NULL, "material is NULL");
    CHECK(m->iMat == iMat, "iMat %d != %d", m->iMat, iMat);
    CHECK(m->nRho >= 2 && m->nT >= 2, "table too small: %d x %d", m->nRho, m->nT);
    CHECK(strlen(m->MatString) > 0, "MatString empty");
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after init", m->nLookupFail);

    for (int j = 1; j < m->nRho; j++)
        CHECK(m->dLogRhoAxis[j] > m->dLogRhoAxis[j-1], "rho axis not increasing at %d", j);
    for (int i = 1; i < m->nT; i++)
        CHECK(m->dLogTAxis[i] > m->dLogTAxis[i-1], "T axis not increasing at %d", i);

    /* The lookup limits are the table bounds. */
    CHECK(m->LogRhoMin == m->dLogRhoAxis[0] && m->LogRhoMax == m->dLogRhoAxis[m->nRho-1], "rho limits are not the table bounds");
    CHECK(m->LogTMin == m->dLogTAxis[0] && m->LogTMax == m->dLogTAxis[m->nT-1], "T limits are not the table bounds");
    CHECK(m->grid.nx == m->nRho && m->grid.ny == m->nT && m->grid.x == m->dLogRhoAxis && m->grid.y == m->dLogTAxis, "interpBilinear grid view inconsistent");

    /* rho0 must be a density the EOS can be evaluated at (EOSSolveBC starts there). */
    CHECK(m->rho0 > 0.0 && log10(m->rho0) > m->LogRhoMin && log10(m->rho0) < m->LogRhoMax, "rho0 = %g outside the table", m->rho0);

    /* Unit factors: identity in cgs, consistent otherwise. */
    if (dKpcUnit > 0.0 && dMsolUnit > 0.0) {
        CHECK(m->dErgPerGmUnit > 0.0 && m->dGmPerCcUnit > 0.0 && isfinite(m->dErgPerGmUnit), "unit factors");
    } else {
        CHECK(m->dErgPerGmUnit == 1.0 && m->dGmPerCcUnit == 1.0, "cgs unit factors are not 1");
    }

    /* Tables must be finite everywhere, including the generated cs table. */
    int bad = 0;
    for (int j = 0; j < m->nRho; j++)
        for (int i = 0; i < m->nT; i++)
            if (!isfinite(LOGP(m, j, i)) || !isfinite(LOGU(m, j, i)) || !isfinite(LOGS(m, j, i)) || !isfinite(LOGC(m, j, i))) bad++;
    CHECK(bad == 0, "%d non-finite table entries", bad);

    /* Thermodynamic consistency of the nodes: u(T), s(T) increasing along rho
     * columns, P(rho) increasing along T rows. Violations are table defects
     * (reported, not failed); the tests below skip them. */
    int nUBad = 0, nSBad = 0, nPBad = 0;
    for (int j = 0; j < m->nRho; j++)
        for (int i = 0; i < m->nT - 1; i++) {
            if (!(LOGU(m, j, i+1) > LOGU(m, j, i))) nUBad++;
            if (!(LOGS(m, j, i+1) > LOGS(m, j, i))) nSBad++;
        }
    for (int i = 0; i < m->nT; i++)
        for (int j = 0; j < m->nRho - 1; j++)
            if (!(LOGP(m, j+1, i) > LOGP(m, j, i))) nPBad++;
    if (nUBad || nSBad || nPBad)
        NOTE("%d cells with du/dT <= 0, %d with ds/dT <= 0, %d with dP/drho <= 0 (table defects)", nUBad, nSBad, nPBad);
    if (nUBad || nSBad) {
        /* Envelope of the defective region: per density decade, the T range
         * (lowest to highest cell top) in which some cell has du/dT <= 0 or
         * ds/dT <= 0, so one can see which simulation states are affected. */
        printf("\n    note: T range of cells with du/dT <= 0 (u) or ds/dT <= 0 (s), per density decade (cgs):");
        int jStart = 0;
        while (jStart < m->nRho) {
            double decadeLo = floor(m->dLogRhoAxis[jStart] + log10(m->dGmPerCcUnit));
            double uLo = INFINITY, uHi = -INFINITY, sLo = INFINITY, sHi = -INFINITY;
            int jEnd = jStart;
            while (jEnd < m->nRho && m->dLogRhoAxis[jEnd] + log10(m->dGmPerCcUnit) < decadeLo + 1.0) {
                for (int i = 0; i < m->nT - 1; i++) {
                    if (!(LOGU(m, jEnd, i+1) > LOGU(m, jEnd, i))) { uLo = fmin(uLo, m->dLogTAxis[i]); uHi = fmax(uHi, m->dLogTAxis[i+1]); }
                    if (!(LOGS(m, jEnd, i+1) > LOGS(m, jEnd, i))) { sLo = fmin(sLo, m->dLogTAxis[i]); sHi = fmax(sHi, m->dLogTAxis[i+1]); }
                }
                jEnd++;
            }
            if (uHi > -INFINITY || sHi > -INFINITY) {
                printf("\n          logrho in [%g, %g):", decadeLo, decadeLo + 1.0);
                if (uHi > -INFINITY) printf("  u: T in [%.3g, %.3g] K", tenpow(uLo), tenpow(uHi));
                if (sHi > -INFINITY) printf("  s: T in [%.3g, %.3g] K", tenpow(sLo), tenpow(sHi));
            }
            jStart = jEnd;
        }
        printf("\n  %-44s", "");
    }
    END_TEST();
}

static void testForwardLookups(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("forward lookups are finite and match the bilinear interpolant");
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double lr = sampleLogRho(m, a), lt = sampleLogT(m, b);
        double rho = tenpow(lr), T = tenpow(lt);
        double u = scvheosUofRhoT(m, rho, T);
        double P = scvheosPofRhoT(m, rho, T);
        double c = scvheosCsofRhoT(m, rho, T);
        double S = scvheosSofRhoT(m, rho, T);
        CHECK(isfinite(u) && u > 0.0, "u(rho=%g,T=%g) = %g", rho, T, u);
        CHECK(isfinite(P) && P > 0.0, "P(rho=%g,T=%g) = %g", rho, T, P);
        CHECK(isfinite(c) && c > 0.0, "c(rho=%g,T=%g) = %g", rho, T, c);
        CHECK(isfinite(S) && S > 0.0, "S(rho=%g,T=%g) = %g", rho, T, S);
        /* The core must agree with an independent bilinear evaluation. */
        CHECK_ABS(scvheosLogPofLogRhoLogT(m, lr, lt), bilinear(m, m->dLogPArray, lr, lt), 1e-12, "logP interpolant at (%g,%g)", lr, lt);
        CHECK_ABS(scvheosLogUofLogRhoLogT(m, lr, lt), bilinear(m, m->dLogUArray, lr, lt), 1e-12, "logu interpolant at (%g,%g)", lr, lt);
        CHECK_ABS(scvheosLogSofLogRhoLogT(m, lr, lt), bilinear(m, m->dLogSArray, lr, lt), 1e-12, "logs interpolant at (%g,%g)", lr, lt);
        CHECK_ABS(scvheosLogCsofLogRhoLogT(m, lr, lt), bilinear(m, m->dLogCArray, lr, lt), 1e-12, "logcs interpolant at (%g,%g)", lr, lt);
    }
    /* Nodes must be reproduced exactly, including the edges: the lookups
     * tolerate the 10^x <-> log10 roundoff at the table edge. The round trip
     * loses |log10 f| ln10 eps per operation, which for u ~ 1e-50 in code
     * units (cold, dense corner) is 1e-13, hence the tolerance. */
    for (int jj = 0; jj <= 8; jj++)
        for (int ii = 0; ii <= 8; ii++) {
            int j = (int)((long)jj * (m->nRho - 1) / 8), i = (int)((long)ii * (m->nT - 1) / 8);
            double rho = tenpow(m->dLogRhoAxis[j]), T = tenpow(m->dLogTAxis[i]);
            CHECK_REL(scvheosUofRhoT(m, rho, T), tenpow(LOGU(m, j, i)), 1e-12, "u at node (%d,%d)", j, i);
            CHECK_REL(scvheosPofRhoT(m, rho, T), tenpow(LOGP(m, j, i)), 1e-12, "P at node (%d,%d)", j, i);
            CHECK_REL(scvheosSofRhoT(m, rho, T), tenpow(LOGS(m, j, i)), 1e-12, "s at node (%d,%d)", j, i);
            CHECK_REL(scvheosCsofRhoT(m, rho, T), tenpow(LOGC(m, j, i)), 1e-12, "cs at node (%d,%d)", j, i);
        }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after in-table lookups", m->nLookupFail);
    END_TEST();
}

static void testTemperatureInversion(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("T(rho,u) and T(rho,s) invert the forward lookups");
    int nSkippedU = 0, nSkippedS = 0;
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double lr = sampleLogRho(m, a), lt = sampleLogT(m, b);
        double rho = tenpow(lr), T = tenpow(lt);
        int j = findCell(m->dLogRhoAxis, m->nRho, lr);

        /* u(T): the closed-form inversion is exact on a monotonic column. */
        if (monotonicInColumn(m, m->dLogUArray, j)) {
            double u = scvheosUofRhoT(m, rho, T);
            double T2 = scvheosTofRhoU(m, rho, u);
            CHECK_REL(T2, T, 1e-12, "T(rho,u) roundtrip at rho=%g T=%g", rho, T);
            CHECK_REL(scvheosPofRhoU(m, rho, u), scvheosPofRhoT(m, rho, T), 1e-12, "P(rho,u) vs P(rho,T) at rho=%g T=%g", rho, T);
            CHECK_REL(scvheosCsofRhoU(m, rho, u), scvheosCsofRhoT(m, rho, T), 1e-12, "c(rho,u) vs c(rho,T) at rho=%g T=%g", rho, T);
        } else nSkippedU++;

        /* s(T): same, through the entropy column. */
        if (monotonicInColumn(m, m->dLogSArray, j)) {
            /* Looser than for u: s(T) has nearly flat stretches, where the
             * roundoff of s maps to a larger error in T. */
            double S = scvheosSofRhoT(m, rho, T);
            double T2 = scvheosTofRhoS(m, rho, S);
            CHECK_REL(T2, T, 1e-9, "T(rho,s) roundtrip at rho=%g T=%g", rho, T);
        } else nSkippedS++;
    }
    /* Exactly at the column ends the inversion must return the table edge. */
    for (int a = 0; a < NSAMPLE_RHO; a += NSAMPLE_RHO / 10 + 1) {
        double lr = sampleLogRho(m, a);
        int j = findCell(m->dLogRhoAxis, m->nRho, lr);
        if (!monotonicInColumn(m, m->dLogUArray, j)) continue;
        double uMin = scvheosLogUofLogRhoLogT(m, lr, m->LogTMin), uMax = scvheosLogUofLogRhoLogT(m, lr, m->LogTMax);
        CHECK_ABS(scvheosLogTofLogRhoLogU(m, lr, uMin), m->LogTMin, 1e-12, "u(T_min) inverts to T_min at logrho=%g", lr);
        CHECK_ABS(scvheosLogTofLogRhoLogU(m, lr, uMax), m->LogTMax, 1e-12, "u(T_max) inverts to T_max at logrho=%g", lr);
    }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after inversions", m->nLookupFail);
    if (nSkippedU || nSkippedS)
        NOTE("%d T(rho,u) and %d T(rho,s) round trips skipped (non-monotonic column, table defect)", nSkippedU, nSkippedS);
    END_TEST();
}

static void testDensityInversion(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("rho(P,T) inverts P(rho,T)");
    int nSkipped = 0;
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double lr = sampleLogRho(m, a), lt = sampleLogT(m, b);
        int i = findCell(m->dLogTAxis, m->nT, lt);
        if (!pressureMonotonicInRow(m, i)) { nSkipped++; continue; }
        double rho = tenpow(lr), T = tenpow(lt);
        double P = scvheosPofRhoT(m, rho, T);
        double rho2 = scvheosRhoofPT(m, P, T);
        CHECK_REL(rho2, rho, 1e-12, "rho(P,T) roundtrip at rho=%g T=%g", rho, T);
    }
    /* Exactly at a node the inversion must return that node. */
    {
        int j = m->nRho / 2, i = m->nT / 2;
        if (pressureMonotonicInRow(m, i)) {
            double rho = scvheosRhoofPT(m, tenpow(LOGP(m, j, i)), tenpow(m->dLogTAxis[i]));
            CHECK_REL(rho, tenpow(m->dLogRhoAxis[j]), 1e-12, "rho(P,T) at node (%d,%d)", j, i);
        }
    }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after density inversions", m->nLookupFail);
    if (nSkipped) NOTE("%d rho(P,T) round trips skipped (P(rho) not increasing in the row, table defect)", nSkipped);
    END_TEST();
}

static void testFusedLookup(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("scvheosPCTofRhoU matches the separate calls");
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double lr = sampleLogRho(m, a);
        int j = findCell(m->dLogRhoAxis, m->nRho, lr);
        if (!monotonicInColumn(m, m->dLogUArray, j)) continue;
        double rho = tenpow(lr), T = sampleT(m, b);
        double u = scvheosUofRhoT(m, rho, T);
        double c, Tf;
        double P = scvheosPCTofRhoU(m, rho, u, &c, &Tf);
        CHECK(!isSentinel(P) && !isSentinel(c) && !isSentinel(Tf), "fused lookup failed at rho=%g T=%g", rho, T);
        CHECK_REL(P, scvheosPofRhoT(m, rho, T), 1e-12, "fused P at rho=%g T=%g", rho, T);
        CHECK_REL(c, scvheosCsofRhoT(m, rho, T), 1e-12, "fused c at rho=%g T=%g", rho, T);
        CHECK_REL(Tf, T, 1e-12, "fused T at rho=%g T=%g", rho, T);
        CHECK_REL(scvheosPofRhoU(m, rho, u), P, 0.0, "PofRhoU vs fused");
        CHECK_REL(scvheosCsofRhoU(m, rho, u), c, 0.0, "CsofRhoU vs fused");
    }
    END_TEST();
}

static void testDerivatives(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("analytic derivatives match the interpolant");
    /* Inside a cell the interpolant is bilinear in (logT, logrho): the log
     * partials are exactly the central differences of the log lookups, and the
     * linear-space derivatives follow from df/dx = (f/x) dlogf/dlogx. */
    for (int j = 1; j < m->nRho - 2; j += (m->nRho - 3) / 8 + 1)
        for (int i = 1; i < m->nT - 2; i += (m->nT - 3) / 8 + 1) {
            double lr = 0.5 * (m->dLogRhoAxis[j] + m->dLogRhoAxis[j+1]);
            double lt = 0.5 * (m->dLogTAxis[i] + m->dLogTAxis[i+1]);
            double hr = 0.01 * (m->dLogRhoAxis[j+1] - m->dLogRhoAxis[j]);
            double ht = 0.01 * (m->dLogTAxis[i+1] - m->dLogTAxis[i]);
            double rho = tenpow(lr), T = tenpow(lt);

            double fd;
            fd = (scvheosLogPofLogRhoLogT(m, lr + hr, lt) - scvheosLogPofLogRhoLogT(m, lr - hr, lt)) / (2 * hr);
            CHECK_NEAR(scvheosdLogPdLogRhoofLogRhoLogT(m, lr, lt), fd, 1e-9, 1.0, "dlogP/dlogrho in cell (%d,%d)", j, i);
            fd = (scvheosLogPofLogRhoLogT(m, lr, lt + ht) - scvheosLogPofLogRhoLogT(m, lr, lt - ht)) / (2 * ht);
            CHECK_NEAR(scvheosdLogPdLogTofLogRhoLogT(m, lr, lt), fd, 1e-9, 1.0, "dlogP/dlogT in cell (%d,%d)", j, i);
            fd = (scvheosLogUofLogRhoLogT(m, lr + hr, lt) - scvheosLogUofLogRhoLogT(m, lr - hr, lt)) / (2 * hr);
            CHECK_NEAR(scvheosdLogUdLogRhoofLogRhoLogT(m, lr, lt), fd, 1e-9, 1.0, "dlogu/dlogrho in cell (%d,%d)", j, i);
            fd = (scvheosLogUofLogRhoLogT(m, lr, lt + ht) - scvheosLogUofLogRhoLogT(m, lr, lt - ht)) / (2 * ht);
            CHECK_NEAR(scvheosdLogUdLogTofLogRhoLogT(m, lr, lt), fd, 1e-9, 1.0, "dlogu/dlogT in cell (%d,%d)", j, i);
            fd = (scvheosLogSofLogRhoLogT(m, lr + hr, lt) - scvheosLogSofLogRhoLogT(m, lr - hr, lt)) / (2 * hr);
            CHECK_NEAR(scvheosdLogSdLogRhoofLogRhoLogT(m, lr, lt), fd, 1e-9, 1.0, "dlogs/dlogrho in cell (%d,%d)", j, i);
            fd = (scvheosLogSofLogRhoLogT(m, lr, lt + ht) - scvheosLogSofLogRhoLogT(m, lr, lt - ht)) / (2 * ht);
            CHECK_NEAR(scvheosdLogSdLogTofLogRhoLogT(m, lr, lt), fd, 1e-9, 1.0, "dlogs/dlogT in cell (%d,%d)", j, i);

            /* Linear-space partials: P and u are exponential in log rho / log T
             * inside the cell, so the central difference has an O(h^2) error.
             * Scale: f/x, the prefactor of the log derivative. */
            double rp = tenpow(lr + 1e-3 * hr), rm = tenpow(lr - 1e-3 * hr);
            double Tp = tenpow(lt + 1e-3 * ht), Tm = tenpow(lt - 1e-3 * ht);
            double P = scvheosPofRhoT(m, rho, T), u = scvheosUofRhoT(m, rho, T), S = scvheosSofRhoT(m, rho, T);
            fd = (scvheosPofRhoT(m, rp, T) - scvheosPofRhoT(m, rm, T)) / (rp - rm);
            CHECK_NEAR(scvheosdPdRhoofRhoT(m, rho, T), fd, 1e-7, P / rho, "dP/drho|T in cell (%d,%d)", j, i);
            fd = (scvheosPofRhoT(m, rho, Tp) - scvheosPofRhoT(m, rho, Tm)) / (Tp - Tm);
            CHECK_NEAR(scvheosdPdTofRhoT(m, rho, T), fd, 1e-7, P / T, "dP/dT|rho in cell (%d,%d)", j, i);
            fd = (scvheosUofRhoT(m, rp, T) - scvheosUofRhoT(m, rm, T)) / (rp - rm);
            CHECK_NEAR(scvheosdUdRhoofRhoT(m, rho, T), fd, 1e-7, u / rho, "du/drho|T in cell (%d,%d)", j, i);
            fd = (scvheosUofRhoT(m, rho, Tp) - scvheosUofRhoT(m, rho, Tm)) / (Tp - Tm);
            CHECK_NEAR(scvheosdUdTofRhoT(m, rho, T), fd, 1e-7, u / T, "du/dT|rho in cell (%d,%d)", j, i);
            fd = (scvheosSofRhoT(m, rp, T) - scvheosSofRhoT(m, rm, T)) / (rp - rm);
            CHECK_NEAR(scvheosdSdRhoofRhoT(m, rho, T), fd, 1e-7, S / rho, "ds/drho|T in cell (%d,%d)", j, i);
            fd = (scvheosSofRhoT(m, rho, Tp) - scvheosSofRhoT(m, rho, Tm)) / (Tp - Tm);
            CHECK_NEAR(scvheosdSdTofRhoT(m, rho, T), fd, 1e-7, S / T, "ds/dT|rho in cell (%d,%d)", j, i);

            /* (rho, u) derivatives against finite differences of P(rho, u). The
             * stencil must stay in this cell, which the small step guarantees
             * on a monotonic column; a defective cell gives the sentinel. */
            if (!monotonicInColumn(m, m->dLogUArray, j)) continue;
            double up = scvheosUofRhoT(m, rho, Tp), um = scvheosUofRhoT(m, rho, Tm);
            fd = (scvheosPofRhoU(m, rho, up) - scvheosPofRhoU(m, rho, um)) / (up - um);
            CHECK_NEAR(scvheosdPdUofRhoU(m, rho, u), fd, 1e-7, P / u, "dP/du|rho in cell (%d,%d)", j, i);
            fd = (scvheosPofRhoU(m, rp, u) - scvheosPofRhoU(m, rm, u)) / (rp - rm);
            CHECK_NEAR(scvheosdPdRhoofRhoU(m, rho, u), fd, 1e-7, P / rho, "dP/drho|u in cell (%d,%d)", j, i);
        }
    END_TEST();
}

static void testSoundSpeedTable(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("sound-speed table from the nodes");
    /* Interior nodes: recompute cs^2 = dP/drho + T/(rho^2 cv) (dP/dT)^2 with
     * second-order central differences of the log tables, independently of
     * the library, and compare. Nodes where that gives cs^2 <= 0 (table
     * defect) carry the documented floor P/rho instead. */
    int nFloor = 0;
    for (int j = 1; j < m->nRho - 1; j++)
        for (int i = 1; i < m->nT - 1; i++) {
            double h0, h1, fm, f0, fp;
            double rho = tenpow(m->dLogRhoAxis[j]), T = tenpow(m->dLogTAxis[i]);
            double P = tenpow(LOGP(m, j, i)), u = tenpow(LOGU(m, j, i));
            h0 = m->dLogRhoAxis[j] - m->dLogRhoAxis[j-1]; h1 = m->dLogRhoAxis[j+1] - m->dLogRhoAxis[j];
            fm = LOGP(m, j-1, i); f0 = LOGP(m, j, i); fp = LOGP(m, j+1, i);
            double dlPdlr = (h0*h0*fp - h1*h1*fm + (h1*h1 - h0*h0)*f0) / (h0*h1*(h0 + h1));
            h0 = m->dLogTAxis[i] - m->dLogTAxis[i-1]; h1 = m->dLogTAxis[i+1] - m->dLogTAxis[i];
            fm = LOGP(m, j, i-1); f0 = LOGP(m, j, i); fp = LOGP(m, j, i+1);
            double dlPdlT = (h0*h0*fp - h1*h1*fm + (h1*h1 - h0*h0)*f0) / (h0*h1*(h0 + h1));
            fm = LOGU(m, j, i-1); f0 = LOGU(m, j, i); fp = LOGU(m, j, i+1);
            double dludlT = (h0*h0*fp - h1*h1*fm + (h1*h1 - h0*h0)*f0) / (h0*h1*(h0 + h1));
            double cv = u / T * dludlT;
            double cs2 = P / rho * dlPdlr + T / (rho * rho * cv) * (P / T * dlPdlT) * (P / T * dlPdlT);
            double csTable = tenpow(LOGC(m, j, i));
            if (isfinite(cs2) && cs2 > 0.0) {
                CHECK_REL(csTable, sqrt(cs2), 1e-10, "cs at node (%d,%d)", j, i);
            } else {
                nFloor++;
                CHECK_REL(csTable, sqrt(P / rho), 1e-10, "cs floor at defective node (%d,%d)", j, i);
            }
        }
    if (nFloor) NOTE("%d interior nodes carry the cs^2 = P/rho floor (table defect)", nFloor);
    END_TEST();
}

static void testIsentrope(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("isentropic evolution conserves tabulated entropy");
    int nSkipped = 0;
    long nf0 = m->nLookupFail;
    for (int a = 1; a < NSAMPLE_RHO - 1; a++) for (int b = 2; b < NSAMPLE_T - 2; b++) {
        double lr = sampleLogRho(m, a), T1 = sampleT(m, b);
        double rho1 = tenpow(lr);
        int j1 = findCell(m->dLogRhoAxis, m->nRho, lr);
        if (!monotonicInColumn(m, m->dLogUArray, j1) || !monotonicInColumn(m, m->dLogSArray, j1)) { nSkipped++; continue; }
        double u1 = scvheosUofRhoT(m, rho1, T1);
        double S1 = scvheosSofRhoT(m, rho1, T1);

        /* Identity: same rho reproduces u. The s(T) round trip amplifies
         * roundoff on flat stretches of the entropy column. */
        CHECK_REL(scvheosIsentropicU(m, rho1, u1, rho1), u1, 1e-9, "isentrope to same rho at rho=%g T=%g", rho1, T1);

        /* Small compression and expansion inside the table: S is conserved
         * wherever the isentrope stays inside the s column at rho2. */
        for (int dir = -1; dir <= 1; dir += 2) {
            double lr2 = lr + dir * 0.5 * (m->dLogRhoAxis[j1+1] - m->dLogRhoAxis[j1]);
            if (lr2 < m->LogRhoMin || lr2 > m->LogRhoMax) continue;
            int j2 = findCell(m->dLogRhoAxis, m->nRho, lr2);
            if (!monotonicInColumn(m, m->dLogUArray, j2) || !monotonicInColumn(m, m->dLogSArray, j2)) { nSkipped++; continue; }
            double rho2 = tenpow(lr2);
            double sLo = scvheosSofRhoT(m, rho2, tenpow(m->LogTMin)), sHi = scvheosSofRhoT(m, rho2, tenpow(m->LogTMax));
            CHECK(!isSentinel(sLo) && !isSentinel(sHi), "s at the column ends failed at rho=%g", rho2);
            if (!(S1 > sLo && S1 < sHi)) continue;          /* leaves the table: clamp by contract, tested below */
            long nf = m->nLookupFail;
            double u2 = scvheosIsentropicU(m, rho1, u1, rho2);
            CHECK(!isSentinel(u2), "isentropicU failed rho=%g->%g", rho1, rho2);
            CHECK_REL(scvheosSofRhoT(m, rho2, scvheosTofRhoU(m, rho2, u2)), S1, 1e-10, "S conserved rho=%g->%g T=%g", rho1, rho2, T1);
            CHECK(m->nLookupFail == nf, "in-table isentrope counted as a fallback");
        }
    }
    /* Leaving the table along the isentrope: clamp to the cold curve (below)
     * or the top row (above), counted in nLookupFail, like ANEOSisentropicU. */
    {
        double lr1 = sampleLogRho(m, NSAMPLE_RHO - 1), lr2 = sampleLogRho(m, 0);   /* strong expansion */
        int j1 = findCell(m->dLogRhoAxis, m->nRho, lr1), j2 = findCell(m->dLogRhoAxis, m->nRho, lr2);
        if (monotonicInColumn(m, m->dLogUArray, j1) && monotonicInColumn(m, m->dLogSArray, j1) &&
            monotonicInColumn(m, m->dLogUArray, j2) && monotonicInColumn(m, m->dLogSArray, j2)) {
            double rho1 = tenpow(lr1), rho2 = tenpow(lr2);
            double T1 = tenpow(m->LogTMin + 0.01 * (m->LogTMax - m->LogTMin));
            double u1 = scvheosUofRhoT(m, rho1, T1), S1 = scvheosSofRhoT(m, rho1, T1);
            double sLo = scvheosSofRhoT(m, rho2, tenpow(m->LogTMin));
            if (S1 < sLo) {
                long nf = m->nLookupFail;
                double u2 = scvheosIsentropicU(m, rho1, u1, rho2);
                CHECK_REL(u2, scvheosUofRhoT(m, rho2, tenpow(m->LogTMin)), 1e-12, "isentrope below the table clamps to the cold curve");
                CHECK(m->nLookupFail == nf + 1, "clamp not counted in nLookupFail (%ld -> %ld)", nf, m->nLookupFail);
            }
        }
    }
    m->nLookupFail = nf0;
    if (nSkipped) NOTE("%d isentrope checks skipped (non-monotonic column, table defect)", nSkipped);
    END_TEST();
}

static void testExtrapLimit(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("scvheosIsInExtrapLimit codes");
    double rhoOk = sampleRho(m, NSAMPLE_RHO / 2);
    double uOk = scvheosUofRhoT(m, rhoOk, sampleT(m, NSAMPLE_T / 2));
    double rhoLo = 0.5 * tenpow(m->LogRhoMin), rhoHi = 2.0 * tenpow(m->LogRhoMax);
    double uLo = 0.5 * scvheosUofRhoT(m, rhoOk, tenpow(m->LogTMin));
    double uHi = 2.0 * scvheosUofRhoT(m, rhoOk, tenpow(m->LogTMax));

    CHECK(uLo > 0.0 && uHi > 0.0, "u at the column ends failed");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, uOk) == SCVHEOS_SUCCESS, "in table");
    /* Exactly on the edges, through the linear-space round trip. */
    CHECK(scvheosIsInExtrapLimit(m, tenpow(m->LogRhoMin), uOk) == SCVHEOS_SUCCESS, "on rhoMin");
    CHECK(scvheosIsInExtrapLimit(m, tenpow(m->LogRhoMax), scvheosUofRhoT(m, tenpow(m->LogRhoMax), sampleT(m, NSAMPLE_T / 2))) == SCVHEOS_SUCCESS, "on rhoMax");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, 2.0 * uLo) == SCVHEOS_SUCCESS, "on the cold curve");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, 0.5 * uHi) == SCVHEOS_SUCCESS, "on the top row");
    CHECK(scvheosIsInExtrapLimit(m, rhoLo, uOk) == SCVHEOS_OUTSIDE_RHOMIN, "below rhoMin");
    CHECK(scvheosIsInExtrapLimit(m, rhoHi, uOk) == SCVHEOS_OUTSIDE_RHOMAX, "above rhoMax");
    CHECK(scvheosIsInExtrapLimit(m, -1.0, uOk) == SCVHEOS_OUTSIDE_RHOMIN, "negative rho");
    CHECK(scvheosIsInExtrapLimit(m, 0.0, uOk) == SCVHEOS_OUTSIDE_RHOMIN, "zero rho");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, uLo) == SCVHEOS_OUTSIDE_TMIN, "below cold curve");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, uHi) == SCVHEOS_OUTSIDE_TMAX, "above top row");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, -1.0) == SCVHEOS_OUTSIDE_TMIN, "negative u");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, 0.0) == SCVHEOS_OUTSIDE_TMIN, "zero u");
    CHECK(scvheosIsInExtrapLimit(m, rhoOk, NAN) == SCVHEOS_OUTSIDE_TMIN, "NaN u");
    /* Consistency: anything the status accepts must evaluate without a sentinel. */
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double lr = sampleLogRho(m, a);
        int j = findCell(m->dLogRhoAxis, m->nRho, lr);
        if (!monotonicInColumn(m, m->dLogUArray, j)) continue;
        double rho = tenpow(lr), u = scvheosUofRhoT(m, rho, sampleT(m, b));
        if (scvheosIsInExtrapLimit(m, rho, u) == SCVHEOS_SUCCESS) {
            double c, T;
            double P = scvheosPCTofRhoU(m, rho, u, &c, &T);
            CHECK(!isSentinel(P) && !isSentinel(c) && !isSentinel(T), "in-table state produced a sentinel at rho=%g", rho);
        } else {
            CHECK(0, "in-table state rejected at rho=%g u=%g", rho, u);
        }
    }
    END_TEST();
}

#ifdef SCVH_NOEXIT
static void testFailureContract(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("out-of-table inputs return SCVHEOS_FAIL_VALUE");
    long nf0 = m->nLookupFail;
    double rhoOk = sampleRho(m, NSAMPLE_RHO / 2), TOk = sampleT(m, NSAMPLE_T / 2);
    double rhoLo = 0.5 * tenpow(m->LogRhoMin), rhoHi = 2.0 * tenpow(m->LogRhoMax);
    double TLo = 0.5 * tenpow(m->LogTMin), THi = 2.0 * tenpow(m->LogTMax);
    double uOk = scvheosUofRhoT(m, rhoOk, TOk);
    double uLo = 0.5 * scvheosUofRhoT(m, rhoOk, tenpow(m->LogTMin));
    double uHi = 2.0 * scvheosUofRhoT(m, rhoOk, tenpow(m->LogTMax));
    double c, T, P;

    CHECK(uLo > 0.0 && uHi > 0.0, "u at the column ends failed");

    /* Forward lookups. */
    CHECK(scvheosUofRhoT(m, rhoLo, TOk) == SCVHEOS_FAIL_VALUE, "UofRhoT below rhoMin");
    CHECK(scvheosUofRhoT(m, rhoHi, TOk) == SCVHEOS_FAIL_VALUE, "UofRhoT above rhoMax");
    CHECK(scvheosUofRhoT(m, rhoOk, TLo) == SCVHEOS_FAIL_VALUE, "UofRhoT below TMin");
    CHECK(scvheosUofRhoT(m, rhoOk, THi) == SCVHEOS_FAIL_VALUE, "UofRhoT above TMax");
    CHECK(scvheosPofRhoT(m, rhoHi, TOk) == SCVHEOS_FAIL_VALUE, "PofRhoT above rhoMax");
    CHECK(scvheosSofRhoT(m, rhoHi, TOk) == SCVHEOS_FAIL_VALUE, "SofRhoT above rhoMax");
    CHECK(scvheosCsofRhoT(m, rhoHi, TOk) == SCVHEOS_FAIL_VALUE, "CsofRhoT above rhoMax");
    CHECK(scvheosUofRhoT(m, -1.0, TOk) == SCVHEOS_FAIL_VALUE, "UofRhoT negative rho");
    CHECK(scvheosUofRhoT(m, rhoOk, -1.0) == SCVHEOS_FAIL_VALUE, "UofRhoT negative T");

    /* Inversions: u and P fail outside the column; s clamps and counts. */
    CHECK(scvheosTofRhoU(m, rhoOk, uLo) == SCVHEOS_FAIL_VALUE, "TofRhoU below cold curve");
    CHECK(scvheosTofRhoU(m, rhoOk, uHi) == SCVHEOS_FAIL_VALUE, "TofRhoU above top row");
    CHECK(scvheosTofRhoU(m, rhoHi, uOk) == SCVHEOS_FAIL_VALUE, "TofRhoU above rhoMax");
    CHECK(scvheosTofRhoU(m, rhoOk, -1.0) == SCVHEOS_FAIL_VALUE, "TofRhoU negative u");
    CHECK(scvheosPofRhoU(m, rhoOk, uHi) == SCVHEOS_FAIL_VALUE, "PofRhoU above top row");
    CHECK(scvheosCsofRhoU(m, rhoOk, uLo) == SCVHEOS_FAIL_VALUE, "CsofRhoU below cold curve");
    CHECK(scvheosRhoofPT(m, 1e-3 * scvheosPofRhoT(m, tenpow(m->LogRhoMin), TOk), TOk) == SCVHEOS_FAIL_VALUE, "RhoofPT below the row");
    CHECK(scvheosRhoofPT(m, 1e3 * scvheosPofRhoT(m, tenpow(m->LogRhoMax), TOk), TOk) == SCVHEOS_FAIL_VALUE, "RhoofPT above the row");
    CHECK(scvheosRhoofPT(m, 1.0, THi) == SCVHEOS_FAIL_VALUE, "RhoofPT above TMax");
    CHECK(m->nLookupFail == nf0, "failures counted as fallbacks (%ld -> %ld)", nf0, m->nLookupFail);
    {
        double sLo = 0.5 * scvheosSofRhoT(m, rhoOk, tenpow(m->LogTMin)), sHi = 2.0 * scvheosSofRhoT(m, rhoOk, tenpow(m->LogTMax));
        CHECK(sLo > 0.0 && sHi > 0.0, "s at the column ends failed");
        CHECK_REL(scvheosTofRhoS(m, rhoOk, sLo), tenpow(m->LogTMin), 1e-12, "TofRhoS clamps to TMin");
        CHECK_REL(scvheosTofRhoS(m, rhoOk, sHi), tenpow(m->LogTMax), 1e-12, "TofRhoS clamps to TMax");
        CHECK(m->nLookupFail == nf0 + 2, "clamps not counted (%ld -> %ld)", nf0, m->nLookupFail);
        CHECK(scvheosTofRhoS(m, rhoHi, sLo) == SCVHEOS_FAIL_VALUE, "TofRhoS above rhoMax");
    }

    /* Fused lookup sets all three outputs. */
    c = T = 1.0; P = scvheosPCTofRhoU(m, rhoHi, uOk, &c, &T);
    CHECK(P == SCVHEOS_FAIL_VALUE && c == SCVHEOS_FAIL_VALUE && T == SCVHEOS_FAIL_VALUE, "PCTofRhoU above rhoMax: P=%g c=%g T=%g", P, c, T);
    c = T = 1.0; P = scvheosPCTofRhoU(m, rhoOk, uHi, &c, &T);
    CHECK(P == SCVHEOS_FAIL_VALUE && c == SCVHEOS_FAIL_VALUE && T == SCVHEOS_FAIL_VALUE, "PCTofRhoU above top row: P=%g c=%g T=%g", P, c, T);
    c = T = 1.0; P = scvheosPCTofRhoU(m, rhoOk, -1.0, &c, &T);
    CHECK(P == SCVHEOS_FAIL_VALUE && c == SCVHEOS_FAIL_VALUE && T == SCVHEOS_FAIL_VALUE, "PCTofRhoU negative u: P=%g c=%g T=%g", P, c, T);

    /* Derivatives and the isentrope. */
    CHECK(scvheosdPdRhoofRhoU(m, rhoOk, uHi) == SCVHEOS_FAIL_VALUE, "dPdRhoofRhoU above top row");
    CHECK(scvheosdPdUofRhoU(m, rhoHi, uOk) == SCVHEOS_FAIL_VALUE, "dPdUofRhoU above rhoMax");
    CHECK(scvheosdPdRhoofRhoT(m, rhoOk, THi) == SCVHEOS_FAIL_VALUE, "dPdRhoofRhoT above TMax");
    CHECK(scvheosdPdTofRhoT(m, rhoLo, TOk) == SCVHEOS_FAIL_VALUE, "dPdTofRhoT below rhoMin");
    CHECK(scvheosIsentropicU(m, rhoOk, uOk, rhoHi) == SCVHEOS_FAIL_VALUE, "IsentropicU to rho above rhoMax");
    CHECK(scvheosIsentropicU(m, rhoHi, uOk, rhoOk) == SCVHEOS_FAIL_VALUE, "IsentropicU from rho above rhoMax");
    END_TEST();
}
#endif

/* ------------------------------------------------------------------- main */

int main(int argc, char *argv[])
{
    int iMat = SCVHEOS_HHE_EXT;
    double dKpcUnit = 2.06701e-13;   /* Earth radius */
    double dMsolUnit = 4.80438e-08;  /* Earth mass */

    if (argc > 1) iMat = atoi(argv[1]);
    if (argc > 3) { dKpcUnit = atof(argv[2]); dMsolUnit = atof(argv[3]); }

    printf("scvheos %s (interpBilinear %s) test suite: iMat = %d, dKpcUnit = %g, dMsolUnit = %g\n",
           SCVHEOS_VERSION_TEXT, INTERP_VERSION_TEXT, iMat, dKpcUnit, dMsolUnit);

    SCVHEOSMAT *m = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("  %s: nRho = %d, nT = %d\n", m->MatString, m->nRho, m->nT);

    testInit(m, iMat, dKpcUnit, dMsolUnit);
    testForwardLookups(m);
    testTemperatureInversion(m);
    testDensityInversion(m);
    testFusedLookup(m);
    testDerivatives(m);
    testSoundSpeedTable(m);
    testIsentrope(m);
    testExtrapLimit(m);
#ifdef SCVH_NOEXIT
    testFailureContract(m);   /* last: deliberately triggers failure messages */
#else
    printf("  (failure-contract test skipped: build with -DSCVH_NOEXIT)\n");
#endif

    scvheosFinalizeMaterial(m);

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}