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
 * The failure contract (SCVHEOS_FAIL_VALUE, the T_min fallback of T(rho, u),
 * the T_min/T_max clamps of T(rho, s), nLookupFail) is tested last; with
 * -DEOSLIB_VERBOSE the library prints a diagnostic for every failed lookup,
 * so messages are expected there.
 */
/* fork() and waitpid() for the fatal paths (POSIX). */
#define _POSIX_C_SOURCE 200809L
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>   /* open in scvheosInitQuiet */
#include <sys/wait.h>
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

static int isSentinel(double v) { return scvheosFailed(v); }

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
    /* Argument checks: a NULL material is rejected before any dereference, and an
     * existing table is not regenerated (or leaked). */
    CHECK(scvheosGenerateSoundSpeedTable(NULL) == SCVHEOS_FAIL, "NULL material must be rejected");
    {
        double *c0 = m->dLogCArray;
        CHECK(scvheosGenerateSoundSpeedTable(m) == SCVHEOS_FAIL && m->dLogCArray == c0, "existing sound-speed table must not be regenerated");
    }
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
    /*
     * A start point above the u column: the isentrope starts from the top row s(rho1, T_max), the
     * same isentrope as from u(rho1, T_max) itself, counted once (up to scvheos 1.3.0 it started
     * from the T_min fallback and u2 was the cold curve at rho2). Below the cold curve (and for
     * u1 <= 0) it starts from s(rho1, T_min), the isentrope from the cold curve.
     */
    {
        static const double fac[5] = {1.0 / 5.0, 1.0 / 1.2, 1.0, 1.2, 5.0};
        const double Tmax = tenpow(m->LogTMax), Tmin = tenpow(m->LogTMin);
        int nTested = 0;
        for (int a = 1; a < NSAMPLE_RHO - 1; a += 7) {
            double lr1 = sampleLogRho(m, a), rho1 = tenpow(lr1);
            int j1 = findCell(m->dLogRhoAxis, m->nRho, lr1);
            if (!monotonicInColumn(m, m->dLogUArray, j1) || !monotonicInColumn(m, m->dLogSArray, j1)) continue;
            double uTop = scvheosUofRhoT(m, rho1, Tmax), uCold = scvheosUofRhoT(m, rho1, Tmin);
            const double u1s[3] = { 2.0 * uTop, 0.5 * uCold, -1.0 };
            const double uRef1[3] = { uTop, uCold, uCold };
            CHECK(scvheosIsInExtrapLimit(m, rho1, u1s[0]) == SCVHEOS_OUTSIDE_TMAX, "u above the column at rho=%g", rho1);
            for (int f = 0; f < 5; f++) {
                double rho2 = fac[f] * rho1;
                if (log10(rho2) < m->LogRhoMin || log10(rho2) > m->LogRhoMax) continue;
                for (int q = 0; q < 3; q++) {
                    double uRef = scvheosIsentropicU(m, rho1, uRef1[q], rho2);
                    long nf = m->nLookupFail;
                    double u2 = scvheosIsentropicU(m, rho1, u1s[q], rho2);
                    CHECK(m->nLookupFail == nf + 1, "isentrope from u1=%g outside the column counted %ld times (rho=%g->%g)", u1s[q], m->nLookupFail - nf, rho1, rho2);
                    CHECK(!isSentinel(u2), "isentrope from u1=%g outside the column failed (rho=%g->%g)", u1s[q], rho1, rho2);
                    CHECK_REL(u2, uRef, 1e-12, "isentrope from u1=%g = isentrope from the column end (rho=%g->%g)", u1s[q], rho1, rho2);
                    if (f == 2) CHECK_REL(u2, uRef1[q], 1e-9, "isentrope from u1=%g to the same rho is the column end (rho=%g)", u1s[q], rho1);
                }
                CHECK(scvheosIsentropicU(m, rho1, u1s[0], rho2) > scvheosUofRhoT(m, rho2, Tmin) * (1.0 + 1e-12),
                      "isentrope from above the column on the cold curve (rho=%g->%g)", rho1, rho2);
                nTested++;
            }
        }
        CHECK(nTested > 0, "no isentrope from outside the column was tested");
        NOTE("%d isentropes started outside the u column", nTested);
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

/* The cold curve accessor EOSlib uses: u(rho, T_min) on the first T row, the sentinel (counted) outside. */
static void testUCold(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("scvheosUCold");
    long nf0 = m->nLookupFail;
    for (int a = 0; a < NSAMPLE_RHO; a++) {
        double lr = sampleLogRho(m, a), rho = tenpow(lr);
        /* The reference goes through 10^x <-> log10 round trips (also of T_min), each losing
         * |log10 u| ln10 eps; u reaches 1e-42 in code units in the real tables (cold dense corner).
         * scvheosUCold reads the first row directly and is the more exact of the two. */
        double uref = scvheosUofRhoT(m, rho, tenpow(m->LogTMin));
        double tol = 1e-14 + 16.0 * DBL_EPSILON * 2.302585092994046 * (1.0 + fabs(log10(uref)));
        CHECK_REL(scvheosUCold(m, rho), uref, tol, "UCold = u(rho, T_min) at rho=%g", rho);
        CHECK_ABS(log10(scvheosUCold(m, rho)), bilinear(m, m->dLogUArray, lr, m->LogTMin), 1e-12, "UCold on the first row at logrho=%g", lr);
    }
    CHECK_REL(scvheosUCold(m, tenpow(m->LogRhoMax)), tenpow(LOGU(m, m->nRho - 1, 0)), 1e-12, "UCold at the last node");
    CHECK(m->nLookupFail == nf0, "in-table UCold counted");
    CHECK(scvheosFailed(scvheosUCold(m, 2.0 * tenpow(m->LogRhoMax))), "UCold above rhoMax");
    CHECK(scvheosFailed(scvheosUCold(m, -1.0)), "UCold for rho < 0");
    CHECK(m->nLookupFail == nf0 + 2, "UCold failures counted (%ld -> %ld)", nf0, m->nLookupFail);
    m->nLookupFail = nf0;
    END_TEST();
}

/*
 * scvheosReadTable needs two nodes on each axis (one cell): nRho or nT below 2 is rejected before
 * the file is read, like reos3ReadTable (up to scvheos 1.3.0 it was accepted, and the grid of a
 * single node has no cell). scvheosPrintMat prints rho0 in cgs with the unit (it printed the
 * code-unit value unlabelled).
 */
static void testReadTableSizeAndPrint(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("table size check, rho0 in PrintMat");
    static const int sizes[3][2] = { {1, 76}, {460, 1}, {0, 0} };
    for (int k = 0; k < 3; k++) {
        SCVHEOSMAT *t = (SCVHEOSMAT *) calloc(1, sizeof(SCVHEOSMAT));
        /* The expected message goes to /dev/null. */
        fflush(stderr);
        int err = dup(2), devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) { dup2(devnull, 2); close(devnull); }
        int r = scvheosReadTable(t, "scvh_extended_dt_hydrogen_722_helium_278.data", sizes[k][0], sizes[k][1], 1);
        fflush(stderr);
        if (err >= 0) { dup2(err, 2); close(err); }
        CHECK(r == SCVHEOS_FAIL, "scvheosReadTable with nRho= %d nT= %d returned %d", sizes[k][0], sizes[k][1], r);
        CHECK(t->dLogRhoAxis == NULL && t->dLogTAxis == NULL && t->dLogUArray == NULL, "nothing allocated for nRho= %d nT= %d", sizes[k][0], sizes[k][1]);
        free(t);
    }
    {
        FILE *fp = tmpfile();
        char line[512];
        double rho0 = -1.0;
        CHECK(fp != NULL, "tmpfile");
        if (fp != NULL) {
            scvheosPrintMat(m, fp);
            rewind(fp);
            while (fgets(line, sizeof(line), fp) != NULL) {
                const char *p = strstr(line, "rho0 (cgs units):");
                if (p != NULL) rho0 = atof(p + strlen("rho0 (cgs units):"));
            }
            fclose(fp);
        }
        CHECK_REL(rho0, m->rho0 * m->dGmPerCcUnit, 1e-5, "scvheosPrintMat: rho0 (cgs units)");
    }
    END_TEST();
}

static void testFailureContract(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("failure contract: sentinel, T_min fallback, count");
    long nf0 = m->nLookupFail, nf;
    double rhoOk = sampleRho(m, NSAMPLE_RHO / 2), TOk = sampleT(m, NSAMPLE_T / 2);
    double rhoLo = 0.5 * tenpow(m->LogRhoMin), rhoHi = 2.0 * tenpow(m->LogRhoMax);
    double TLo = 0.5 * tenpow(m->LogTMin), THi = 2.0 * tenpow(m->LogTMax);
    double TMin = tenpow(m->LogTMin), TMax = tenpow(m->LogTMax);
    double uOk = scvheosUofRhoT(m, rhoOk, TOk);
    double uLo = 0.5 * scvheosUofRhoT(m, rhoOk, TMin);
    double uHi = 2.0 * scvheosUofRhoT(m, rhoOk, TMax);
    double c, T, P;

    CHECK(uLo > 0.0 && uHi > 0.0, "u at the column ends failed");

    /* Forward lookups: sentinel, each counted once. */
    nf = m->nLookupFail;
    CHECK(scvheosFailed(scvheosUofRhoT(m, rhoLo, TOk)), "UofRhoT below rhoMin");
    CHECK(scvheosFailed(scvheosUofRhoT(m, rhoHi, TOk)), "UofRhoT above rhoMax");
    CHECK(scvheosFailed(scvheosUofRhoT(m, rhoOk, TLo)), "UofRhoT below TMin");
    CHECK(scvheosFailed(scvheosUofRhoT(m, rhoOk, THi)), "UofRhoT above TMax");
    CHECK(scvheosFailed(scvheosPofRhoT(m, rhoHi, TOk)), "PofRhoT above rhoMax");
    CHECK(scvheosFailed(scvheosSofRhoT(m, rhoHi, TOk)), "SofRhoT above rhoMax");
    CHECK(scvheosFailed(scvheosCsofRhoT(m, rhoHi, TOk)), "CsofRhoT above rhoMax");
    CHECK(scvheosFailed(scvheosUofRhoT(m, -1.0, TOk)), "UofRhoT negative rho");
    CHECK(scvheosFailed(scvheosUofRhoT(m, rhoOk, -1.0)), "UofRhoT negative T");
    CHECK(m->nLookupFail == nf + 9, "forward failures not counted (%ld -> %ld)", nf, m->nLookupFail);

    /* T(rho, u) at a valid rho never fails: it falls back to T_min and counts. */
    nf = m->nLookupFail;
    CHECK_REL(scvheosTofRhoU(m, rhoOk, uLo), TMin, 1e-12, "TofRhoU below cold curve -> T_min");
    CHECK_REL(scvheosTofRhoU(m, rhoOk, uHi), TMin, 1e-12, "TofRhoU above top row -> T_min");
    CHECK_REL(scvheosTofRhoU(m, rhoOk, -1.0), TMin, 1e-12, "TofRhoU negative u -> T_min");
    CHECK_REL(scvheosTofRhoU(m, rhoOk, NAN), TMin, 1e-12, "TofRhoU NaN u -> T_min");
    CHECK(m->nLookupFail == nf + 4, "fallbacks not counted (%ld -> %ld)", nf, m->nLookupFail);
    /* ... and the (rho, u) functions then evaluate at (rho, T_min). */
    CHECK_REL(scvheosPofRhoU(m, rhoOk, uHi), scvheosPofRhoT(m, rhoOk, TMin), 1e-12, "PofRhoU above top row = P(rho, T_min)");
    CHECK_REL(scvheosCsofRhoU(m, rhoOk, uLo), scvheosCsofRhoT(m, rhoOk, TMin), 1e-12, "CsofRhoU below cold curve = cs(rho, T_min)");
    c = T = 1.0; P = scvheosPCTofRhoU(m, rhoOk, uHi, &c, &T);
    CHECK_REL(T, TMin, 1e-12, "PCTofRhoU above top row: T");
    CHECK_REL(P, scvheosPofRhoT(m, rhoOk, TMin), 1e-12, "PCTofRhoU above top row: P");
    CHECK_REL(c, scvheosCsofRhoT(m, rhoOk, TMin), 1e-12, "PCTofRhoU above top row: c");
    c = T = 1.0; P = scvheosPCTofRhoU(m, rhoOk, -1.0, &c, &T);
    CHECK_REL(T, TMin, 1e-12, "PCTofRhoU negative u: T");
    CHECK_REL(scvheosIsentropicU(m, rhoOk, uLo, rhoOk), scvheosUofRhoT(m, rhoOk, TMin), 1e-9, "IsentropicU from below the cold curve = u(rho, T_min)");
    /* rho outside the table: sentinel for P and c, T_min for T (never fails). */
    CHECK_REL(scvheosTofRhoU(m, rhoHi, uOk), TMin, 1e-12, "TofRhoU above rhoMax -> T_min");
    CHECK_REL(scvheosTofRhoU(m, rhoLo, uOk), TMin, 1e-12, "TofRhoU below rhoMin -> T_min");
    CHECK_REL(scvheosTofRhoU(m, -1.0, uOk), TMin, 1e-12, "TofRhoU negative rho -> T_min");
    c = T = 1.0; P = scvheosPCTofRhoU(m, rhoHi, uOk, &c, &T);
    CHECK(scvheosFailed(P) && scvheosFailed(c) && T == TMin, "PCTofRhoU above rhoMax: P=%g c=%g T=%g", P, c, T);

    /* rho(P, T). */
    CHECK(scvheosFailed(scvheosRhoofPT(m, 1e-3 * scvheosPofRhoT(m, tenpow(m->LogRhoMin), TOk), TOk)), "RhoofPT below the row");
    CHECK(scvheosFailed(scvheosRhoofPT(m, 1e3 * scvheosPofRhoT(m, tenpow(m->LogRhoMax), TOk), TOk)), "RhoofPT above the row");
    CHECK(scvheosFailed(scvheosRhoofPT(m, 1.0, THi)), "RhoofPT above TMax");
    CHECK(scvheosFailed(scvheosRhoofPT(m, -1.0, TOk)), "RhoofPT negative P");

    /* T(rho, s) clamps to the edges and counts. */
    {
        double sLo = 0.5 * scvheosSofRhoT(m, rhoOk, TMin), sHi = 2.0 * scvheosSofRhoT(m, rhoOk, TMax);
        CHECK(sLo > 0.0 && sHi > 0.0, "s at the column ends failed");
        nf = m->nLookupFail;
        CHECK_REL(scvheosTofRhoS(m, rhoOk, sLo), TMin, 1e-12, "TofRhoS clamps to TMin");
        CHECK_REL(scvheosTofRhoS(m, rhoOk, sHi), TMax, 1e-12, "TofRhoS clamps to TMax");
        CHECK(m->nLookupFail == nf + 2, "clamps not counted (%ld -> %ld)", nf, m->nLookupFail);
        CHECK(scvheosFailed(scvheosTofRhoS(m, rhoHi, sLo)), "TofRhoS above rhoMax");
    }

    /* Derivatives and the isentrope. */
    CHECK(scvheosFailed(scvheosdPdUofRhoU(m, rhoHi, uOk)), "dPdUofRhoU above rhoMax");
    CHECK(scvheosFailed(scvheosdPdRhoofRhoT(m, rhoOk, THi)), "dPdRhoofRhoT above TMax");
    CHECK(scvheosFailed(scvheosdPdTofRhoT(m, rhoLo, TOk)), "dPdTofRhoT below rhoMin");
    CHECK(scvheosFailed(scvheosdLogPdLogRhoofLogRhoLogT(m, m->LogRhoMax + 1.0, log10(TOk))), "dLogPdLogRho above rhoMax");
    CHECK(scvheosFailed(scvheosIsentropicU(m, rhoOk, uOk, rhoHi)), "IsentropicU to rho above rhoMax");
    CHECK(scvheosFailed(scvheosIsentropicU(m, rhoHi, uOk, rhoOk)), "IsentropicU from rho above rhoMax");

    CHECK(m->nLookupFail > nf0, "nothing counted (%ld -> %ld)", nf0, m->nLookupFail);
    m->nLookupFail = nf0;
    END_TEST();
}

/*
 * Initialise a material with the unit pair (dKpcUnit, dMsolUnit) in a child process. Returns 0 if
 * the library accepted it, 1 if it stopped the program (assert(0), or exit(1) with -DNDEBUG), -1
 * for anything else.
 */
static int initUnitsRejected(int iMat, double dKpcUnit, double dMsolUnit)
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);   /* the expected diagnostic */
        SCVHEOSMAT *m = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
        scvheosFinalizeMaterial(m);
        _exit(0);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    if ((WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) || (WIFEXITED(status) && WEXITSTATUS(status) == 1)) return 1;
    return -1;
}

/*
 * One rule for the unit pair in EOSlib and every sublibrary: both units > 0 (and finite) are code
 * units, both <= 0 is cgs, anything else is fatal. (Before, exactly one unit <= 0 made igeos divide
 * by zero, Tillotson assert, and ANEOS, SCvH and REOS3 silently use cgs.)
 */
static void testUnitPair(int iMat)
{
    int failsBefore = nFails;
    TEST("unit pair: code units, cgs or fatal");
    static const struct { double dKpcUnit, dMsolUnit; int bRejected; const char *what; } cases[] = {
        {2.06701e-13, 4.80438e-08, 0, "both > 0: code units"},
        {0.0, 0.0, 0, "both 0: cgs"},
        {-1.0, -1.0, 0, "both < 0: cgs"},
        {0.0, -1.0, 0, "0 and < 0: cgs"},
        {2.06701e-13, 0.0, 1, "dMsolUnit = 0 only"},
        {0.0, 4.80438e-08, 1, "dKpcUnit = 0 only"},
        {-1.0, 4.80438e-08, 1, "dKpcUnit < 0 only"},
        {2.06701e-13, NAN, 1, "dMsolUnit NaN"},
        {INFINITY, 4.80438e-08, 1, "dKpcUnit inf"},
        {NAN, NAN, 1, "both NaN"},
    };
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        int r = initUnitsRejected(iMat, cases[k].dKpcUnit, cases[k].dMsolUnit);
        CHECK(r == cases[k].bRejected, "%s: %s (got %d)", cases[k].what, cases[k].bRejected ? "must be fatal" : "must be accepted", r);
    }
    END_TEST();
}

/* ------------------------------------------------------------------- main */

/*
 * The border strips of the table: the first and last density cell over all temperature cells, and
 * the first and last temperature cell over all density cells. The interior tests sample 5..95% of
 * the log axes and never reach them (an error confined to an edge cell, like the last-cell bug of
 * ANEOSRhoofPU up to ANEOSmaterial 1.2.0, would go unnoticed). Per point at off-grid fractions of
 * the cell, in log space:
 *   - P, u, cs (and s) are 10^(bilinear value of the log table in the cell);
 *   - P, cs, T of the fused lookup equal the separate (rho, u) calls bit for bit;
 *   - T(rho, u) is a root of u(rho, T) = u, and T itself where the u column is monotonic (tolerance
 *     from the cell's dlogu/dlogT); on a non-monotonic column T(rho, u) may be another root, or a
 *     counted T_min fallback where u lies outside [u(T_min), u(T_max)] of the column (a dip of
 *     u(T) below the cold row);
 *   - rho(P, T) reproduces P where the P row is monotonic;
 *   - no other in-table lookup counts a failure.
 */
static void testEdgeCells(SCVHEOSMAT *m)
{
    int failsBefore = nFails;
    TEST("edge cells: first and last rho and T cells");
    static const double frac[3] = { 0.31, 0.53, 0.77 };
    const double eps = DBL_EPSILON, ln10 = 2.302585092994045684;
    long nPoints = 0, nNonMono = 0, nFallback = 0, nOtherRoot = 0;

    for (int strip = 0; strip < 4; strip++) {
        const int nAlong = (strip < 2) ? m->nT - 1 : m->nRho - 1;
        const int stride = nAlong / 300 + 1;
        for (int k = 0; k < nAlong; k += stride) {
            const int j = (strip == 0) ? 0 : (strip == 1) ? m->nRho - 2 : k;   /* rho cell */
            const int i = (strip == 2) ? 0 : (strip == 3) ? m->nT - 2 : k;     /* T cell */
            const int bUMono = monotonicInColumn(m, m->dLogUArray, j);
            for (int fx = 0; fx < 3; fx++) for (int fy = 0; fy < 3; fy++) {
                const double lr = m->dLogRhoAxis[j] + frac[fx] * (m->dLogRhoAxis[j+1] - m->dLogRhoAxis[j]);
                const double lt = m->dLogTAxis[i] + frac[fy] * (m->dLogTAxis[i+1] - m->dLogTAxis[i]);
                const double rho = tenpow(lr), T = tenpow(lt);
                /* The roundoff of a value through 10^x and log10: |log10 f| ln10 eps per operation. */
                const double lscale = fabs(lr) + fabs(lt) + 1.0;
                long nf = m->nLookupFail;
                nPoints++;

                double P = scvheosPofRhoT(m, rho, T), u = scvheosUofRhoT(m, rho, T), cs = scvheosCsofRhoT(m, rho, T);
                const double lP = bilinear(m, m->dLogPArray, lr, lt), lu = bilinear(m, m->dLogUArray, lr, lt), lc = bilinear(m, m->dLogCArray, lr, lt);
                CHECK(!isSentinel(P) && fabs(P / tenpow(lP) - 1.0) <= 64 * eps * ln10 * (lscale + fabs(lP)), "P in edge cell (%d,%d): %.17g vs %.17g", j, i, P, tenpow(lP));
                CHECK(!isSentinel(u) && fabs(u / tenpow(lu) - 1.0) <= 64 * eps * ln10 * (lscale + fabs(lu)), "u in edge cell (%d,%d): %.17g vs %.17g", j, i, u, tenpow(lu));
                CHECK(!isSentinel(cs) && fabs(cs / tenpow(lc) - 1.0) <= 64 * eps * ln10 * (lscale + fabs(lc)), "cs in edge cell (%d,%d)", j, i);
                {
                    double S = scvheosSofRhoT(m, rho, T); const double ls = bilinear(m, m->dLogSArray, lr, lt);
                    CHECK(!isSentinel(S) && fabs(S / tenpow(ls) - 1.0) <= 64 * eps * ln10 * (lscale + fabs(ls)), "s in edge cell (%d,%d)", j, i);
                }
                CHECK(m->nLookupFail == nf, "forward lookups in edge cell (%d,%d) counted", j, i);

                /* The fused lookup and the separate calls. */
                nf = m->nLookupFail;
                double cf, Tf, Pf = scvheosPCTofRhoU(m, rho, u, &cf, &Tf);
                const long nfscvheosPCTofRhoU = m->nLookupFail - nf;
                CHECK(Pf == scvheosPofRhoU(m, rho, u) && cf == scvheosCsofRhoU(m, rho, u) && Tf == scvheosTofRhoU(m, rho, u),
                      "scvheosPCTofRhoU vs separate calls in edge cell (%d,%d): P %.17g/%.17g T %.17g/%.17g", j, i, Pf, scvheosPofRhoU(m, rho, u), Tf, scvheosTofRhoU(m, rho, u));

                /* T(rho, u): a root, T itself on a monotonic column, else possibly a counted fallback. */
                if (nfscvheosPCTofRhoU == 0) {
                    const double u2 = scvheosUofRhoT(m, rho, Tf);
                    CHECK(fabs(u2 / u - 1.0) <= 256 * eps * ln10 * (lscale + fabs(lu)), "u(rho, T(rho, u)) in edge cell (%d,%d): T = %g, T(rho, u) = %g", j, i, T, Tf);
                    if (bUMono) {
                        /* dlogT = dlogu / (dlogu/dlogT of the cell). */
                        const double w = frac[fx];
                        const double slope = fabs(((1.0 - w) * (m->dLogUArray[(size_t)j * m->nT + i + 1] - m->dLogUArray[(size_t)j * m->nT + i])
                                                   + w * (m->dLogUArray[(size_t)(j+1) * m->nT + i + 1] - m->dLogUArray[(size_t)(j+1) * m->nT + i]))
                                                  / (m->dLogTAxis[i+1] - m->dLogTAxis[i]));
                        CHECK(fabs(Tf / T - 1.0) <= 256 * eps * ln10 * (lscale + fabs(lu)) * (1.0 + 1.0 / slope), "T(rho, u) in edge cell (%d,%d): %.17g vs %.17g", j, i, Tf, T);
                    } else if (fabs(Tf / T - 1.0) > 1e-6) {
                        nOtherRoot++;
                    }
                } else {
                    CHECK(!bUMono, "T(rho, u) fell back on a monotonic column in edge cell (%d,%d)", j, i);
                    CHECK(nfscvheosPCTofRhoU == 1, "fallback counted %ld times in edge cell (%d,%d)", nfscvheosPCTofRhoU, j, i);
                    nFallback++;
                }
                if (!bUMono) nNonMono++;

                /* rho(P, T) on a monotonic row. */
                if (pressureMonotonicInRow(m, i)) {
                    nf = m->nLookupFail;
                    double rP = scvheosRhoofPT(m, P, T);
                    CHECK(!isSentinel(rP) && m->nLookupFail == nf, "rho(P,T) failed in edge cell (%d,%d)", j, i);
                    if (!isSentinel(rP)) CHECK(fabs(scvheosPofRhoT(m, rP, T) / P - 1.0) <= 1e-10, "P(rho(P,T),T) in edge cell (%d,%d)", j, i);
                }
                m->nLookupFail = nf;
            }
        }
    }
    if (nNonMono)
        NOTE("%ld of %ld edge-cell states on a non-monotonic u column (table defect): %ld T_min fallbacks, %ld other roots",
             nNonMono, nPoints, nFallback, nOtherRoot);
    END_TEST();
}

/*
 * Edge tolerance at a log bound near 0. The code-unit shift of log rho is arbitrary, so a unit pair can
 * put the end of the log rho axis at (or within roundoff of) 0, where a tolerance relative to the
 * bound alone is zero. A state on the table edge, converted from cgs to code units as a hydro code
 * holds it, lands a few ulps off the bound (absolute ~1e-16 after log10) and failed up to scvheos
 * 1.3.0 for about half the unit pairs (92 of 201 below). The band is SCVHEOS_EDGE_TOL*|bound| +
 * SCVHEOS_EDGE_TOL_ABS. Both ends of the rho axis; a density clearly outside still fails.
 */
/* scvheosInitMaterial reports its table checks on stderr; quiet for the 402 materials below. */
static SCVHEOSMAT *scvheosInitQuiet(int iMat, double dKpcUnit, double dMsolUnit)
{
    fflush(stdout); fflush(stderr);
    int out = dup(1), err = dup(2), devnull = open("/dev/null", O_WRONLY);
    if (devnull >= 0) { dup2(devnull, 1); dup2(devnull, 2); close(devnull); }
    SCVHEOSMAT *m = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fflush(stdout); fflush(stderr);
    if (out >= 0) { dup2(out, 1); close(out); }
    if (err >= 0) { dup2(err, 2); close(err); }
    return m;
}

static void testLogBoundAtZero(int iMat, double dKpcUnit, double dMsolUnit)
{
    int failsBefore = nFails;
    TEST("edge tolerance at a log rho bound of 0");
    SCVHEOSMAT *m0 = scvheosInitQuiet(iMat, dKpcUnit, dMsolUnit);
    const double g0 = m0->dGmPerCcUnit;
    const double logRhoEdgeCGS[2] = { m0->dLogRhoAxis[0] + log10(g0), m0->dLogRhoAxis[m0->nRho-1] + log10(g0) };
    const double T = pow(10.0, 0.5 * (m0->dLogTAxis[0] + m0->dLogTAxis[m0->nT-1]));
    int nFailed[2] = {0, 0}, nTried = 0, nDefect = 0;
    scvheosFinalizeMaterial(m0);

    for (int end = 0; end < 2; end++) {
        for (int k = -100; k <= 100; k++) {
            /* dGmPerCcUnit is proportional to dMsolUnit: this puts the axis end at ~0 in code units. */
            double M = dMsolUnit * pow(10.0, logRhoEdgeCGS[end]) / g0 * (1.0 + k * 1e-15);
            SCVHEOSMAT *m = scvheosInitQuiet(iMat, dKpcUnit, M);
            const double bound = end ? m->dLogRhoAxis[m->nRho-1] : m->dLogRhoAxis[0];
            const double rho = pow(10.0, logRhoEdgeCGS[end]) / m->dGmPerCcUnit;   /* the edge, from cgs */
            CHECK(fabs(bound) < 1e-12, "axis end %g not near 0 (k = %d)", bound, k);
            /* The forward lookups at the edge decide whether the edge state is accepted. */
            double P = scvheosPofRhoT(m, rho, T);
            double s = scvheosSofRhoT(m, rho, T);
            double u = scvheosUofRhoT(m, rho, T);
            if (isSentinel(P) || isSentinel(s) || isSentinel(u)) {
                nFailed[end]++;
            } else {
                /*
                 * T(rho, u) on the edge column. On a monotonic column it is a root of u(rho, T) = u and
                 * involves no fallback. The real SCvH table has columns at the rho edges where u(T) is
                 * not monotonic: there the strict inversion may find another root than T (still a root),
                 * or u lies outside [u(T_min), u(T_max)] of the column (a dip below the cold row) and T
                 * falls back to T_min (counted); both are the documented behaviour for a table defect,
                 * not an edge rejection.
                 */
                const int ixEdge = end ? m->nRho - 1 : 0;
                const int bMonotonic = (interpNodeLineMonotonic(&m->grid, m->dLogUArray, INTERP_ALONG_Y, ixEdge) == 0);
                long nf = m->nLookupFail;
                double Tu = scvheosTofRhoU(m, rho, u);
                if (m->nLookupFail == nf) {
                    CHECK_REL(scvheosUofRhoT(m, rho, Tu), u, 1e-9, "u(rho, T(rho, u)) at the edge (end %d, k = %d, T = %g, T(rho, u) = %g)", end, k, T, Tu);
                    if (bMonotonic) CHECK_REL(Tu, T, 1e-9, "T(rho, u) on a monotonic edge column (end %d, k = %d)", end, k);
                } else {
                    CHECK(!bMonotonic, "T(rho, u) fell back on a monotonic edge column (end %d, k = %d)", end, k);
                    nDefect++;
                }
            }
            /* Clearly outside (1e-10 in log rho) still fails. */
            double rhoOut = rho * (end ? 1.0 + 2.3e-10 : 1.0 - 2.3e-10);
            CHECK(isSentinel(scvheosPofRhoT(m, rhoOut, T)), "P(rho, T) %s the table accepted (end %d, k = %d)",
                  end ? "above" : "below", end, k);
            nTried++;
            scvheosFinalizeMaterial(m);
        }
    }
    CHECK(nFailed[0] == 0 && nFailed[1] == 0, "edge state at a log bound of 0 rejected for %d (lower end) and %d (upper end) of %d unit pairs",
          nFailed[0], nFailed[1], nTried / 2);
    if (nDefect) NOTE("%d T(rho, u) at the edge fell back to T_min (non-monotonic u(T) edge column, table defect)", nDefect);
    END_TEST();
}

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
    testUCold(m);
    testReadTableSizeAndPrint(m);
    testFailureContract(m);   /* last: triggers failure messages under EOSLIB_VERBOSE */
    testUnitPair(iMat);
    testEdgeCells(m);
    testLogBoundAtZero(iMat, dKpcUnit, dMsolUnit);

    scvheosFinalizeMaterial(m);

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}