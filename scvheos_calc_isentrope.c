/*
 * Calculate u(rho) for constant entropy for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  01.11.2022
 * Modified:
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "scvheos.h"

/* Isentrope. */
struct ISENTROPE {
    SCVHEOSMAT *Mat;
    int n;
    double *rho;
    double *u;
    double *T;
    double *s;
};

/* Parameters for calculating the derivatives. */
struct deriv_params_type {
    SCVHEOSMAT *Mat;
    int iStatus;
};

/*
 * This function calculates derivative
 *
 * du/drho = P/rho^2
 *
 * x:   Integration variable rho
 * y[]: Array containing u(rho)
 * f[]: Array that returns du/drho
 */
int derivs(double x, const double y[], double f[], void *params) {
    double u = y[0];
    double rho = x;
    struct deriv_params_type *deriv_params = (struct deriv_params_type *) params;
    SCVHEOSMAT *Mat = deriv_params->Mat;
    int iRet;

    /* Be sure to reset iStatus. */
    deriv_params->iStatus = GSL_SUCCESS;

    iRet = scvheosIsInExtrapLimit(Mat, rho, u);

    if (iRet != SCVHEOS_SUCCESS) {
        // CR
        fprintf(stderr, "derivs: iRet = %i\n", iRet);

        switch(iRet) {
            case SCVHEOS_OUTSIDE_RHOMIN:
                fprintf(stderr, "rho= %15.7E logrho= %15.7E smaller than LogRhoMin= %15.7E\n",
                        rho, log10(rho), Mat->LogRhoMin);
                break;
            case SCVHEOS_OUTSIDE_RHOMAX:
                fprintf(stderr, "rho= %15.7E logrho= %15.7E larger than LogRhoMax= %15.7E\n",
                        rho, log10(rho), Mat->LogRhoMax);
                break;
            case SCVHEOS_OUTSIDE_TMIN:
                fprintf(stderr, "u= %15.7E logu= %15.7E smaller than LoguMin= %15.7E\n",
                        u, log10(u), scvheosLogUofLogRhoLogT(Mat, log10(rho), Mat->LogTMin));
                break;
            case SCVHEOS_OUTSIDE_TMAX:
                fprintf(stderr, "u= %15.7E logu= %15.7E larger than LoguMax= %15.7E\n",
                        u, log10(u), scvheosLogUofLogRhoLogT(Mat, log10(rho), Mat->LogTMax));
                break;
        }

        /* Return the error code from EOSlib. */
        deriv_params->iStatus = iRet;
        return GSL_EBADFUNC;
    }

    f[0] = scvheosPofRhoU(Mat, rho, u)/(rho*rho);

    return GSL_SUCCESS;
}

/*
 * Check if a value (rho, T) is inside of the extrapolation limits.
 */
int scvheosIsInExtrapLimitRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho = log10(rho);
    double logT = log10(T);

    if (logrho < Mat->LogRhoMin) return SCVHEOS_OUTSIDE_RHOMIN;
    if (logrho > Mat->LogRhoMax) return SCVHEOS_OUTSIDE_RHOMAX;
    if (logT < Mat->LogTMin) return SCVHEOS_OUTSIDE_TMIN;
    if (logT > Mat->LogTMax) return SCVHEOS_OUTSIDE_TMAX;

    return SCVHEOS_SUCCESS;
}

struct ISENTROPE *SolveIsentrope(SCVHEOSMAT *Mat, double *rhoaxis, double Ti, int nRho) {
    double rho;
    double u;
    double T;
    struct ISENTROPE *isentrope;
    /* ODE system. */
    size_t ndim = 1;
    struct deriv_params_type deriv_params = {Mat, GSL_SUCCESS};
    gsl_odeiv2_system sys = {derivs, NULL, ndim, &deriv_params};
    /* Driver. */
    const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;
    double eps_abs = 0.0;
    double eps_rel = 1e-8;
    double h;
    gsl_odeiv2_driver *driver;
    int iRet;
    
    assert(rhoaxis != NULL);
    assert(nRho > 1);

    /* Assume that the points are equally spaced. */
    h = (rhoaxis[1]-rhoaxis[0])*1e-6;

    driver = gsl_odeiv2_driver_alloc_y_new(&sys, step_type, h, eps_abs, eps_rel);

    isentrope = (struct ISENTROPE *) calloc(1, sizeof(struct ISENTROPE));
    assert(isentrope != NULL);

    isentrope->rho = (double *) calloc(nRho, sizeof(double));
    assert(isentrope->rho != NULL);

    isentrope->u = (double *) calloc(nRho, sizeof(double));
    assert(isentrope->u != NULL);

    isentrope->T = (double *) calloc(nRho, sizeof(double));
    assert(isentrope->T != NULL);

    isentrope->s = (double *) calloc(nRho, sizeof(double));
    assert(isentrope->s != NULL);

    /* Initial values. */
    rho = rhoaxis[0];
    T = Ti;
    
    iRet = scvheosIsInExtrapLimitRhoT(Mat, rho, T);
    
    switch(iRet) {
        case SCVHEOS_OUTSIDE_RHOMIN:
            fprintf(stderr, "rho= %15.7E logrho= %15.7E smaller than LogRhoMin= %15.7E\n",
                    rho, log10(rho), Mat->LogRhoMin);
            break;
        case SCVHEOS_OUTSIDE_RHOMAX:
            fprintf(stderr, "rho= %15.7E logrho= %15.7E larger than LogRhoMax= %15.7E\n",
                    rho, log10(rho), Mat->LogRhoMax);
            break;
        case SCVHEOS_OUTSIDE_TMIN:
            fprintf(stderr, "T= %15.7E logT= %15.7E smaller than LogTMin= %15.7E\n",
                    T, log10(T), Mat->LogTMin);
            break;
        case SCVHEOS_OUTSIDE_TMAX:
            fprintf(stderr, "T= %15.7E logT= %15.7E larger than LogTMax= %15.7E\n",
                    T, log10(T), Mat->LogTMax);
            break;
    }

    // CR
    fprintf(stderr, "Initial values:\n");
    fprintf(stderr, "rho=%15.7E T=%15.7E\n", rho, T);

    u = scvheosUofRhoT(Mat, rho, T);

    isentrope->Mat = Mat;
    isentrope->n = nRho;

    isentrope->rho[0] = rho;
    isentrope->u[0] = u;
    isentrope->T[0] = T;
    isentrope->s[0] = scvheosSofRhoT(Mat, rho, T);
 
    for (int i=1; i<nRho; i++) {
        int status = gsl_odeiv2_driver_apply(driver, &rho, rhoaxis[i], &u);

        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Error: %i (GSL) %i (EOSlib)\n", status, deriv_params.iStatus);

            free(isentrope->rho); 
            free(isentrope->u); 
            free(isentrope->T); 
            free(isentrope->s); 
            free(isentrope);

            return NULL;
        }
 
        T = scvheosTofRhoU(Mat, rho, u);

        isentrope->rho[i] = rho;
        isentrope->u[i] = u;
        isentrope->T[i] = T;
        isentrope->s[i] = scvheosSofRhoT(Mat, rho, T);
    }

    return isentrope;
}

int PrintIsentrope(struct ISENTROPE *isentrope) {
    fprintf(stdout, "#%14s%15s%15s%15s\n", "rho", "u", "T", "s");

    for (int i=0; i<isentrope->n; i++) {
        fprintf(stdout, "%15.7E%15.7E%15.7E%15.7E\n", isentrope->rho[i], isentrope->u[i], isentrope->T[i], isentrope->s[i]);
    }

    return 0;
}

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE_LOWRHOT;
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    struct ISENTROPE *isentrope;
    double rhoi;
    double Ti;
    double rhof;
    double *rhoaxis;
    int nRho;

    dKpcUnit = 0.0;
    dMsolUnit = 0.0;
    //dKpcUnit = 4.84821E-09;
    //dMsolUnit = 9.53869E-04;


    if (argc != 5) {
        fprintf(stderr, "Usage: scvheos_calc_isentrope <rhoi> <Ti> <rhof> <nRho>\n");
        exit(1);
    }

    rhoi = atof(argv[1]);
    Ti = atof(argv[2]);
    rhof = atof(argv[3]);
    nRho = atoi(argv[4]);

    assert(rhoi > 0.0);
    assert(Ti > 0.0);
    assert(rhof > 0.0);
    assert(nRho > 1);

    rhoaxis = (double *) calloc(nRho, sizeof(double));
    assert(rhoaxis != NULL);
    
    /* Generate logarithmic axis. */
    for (int i=0; i<nRho; i++) {
        rhoaxis[i] = rhoi*pow(10.0, (log10(rhof)-log10(rhoi))/(nRho-1)*i);
    }

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    /* Solve ODE. */
    isentrope = SolveIsentrope(Mat, rhoaxis, Ti, nRho);

    PrintIsentrope(isentrope);

    scvheosFinalizeMaterial(Mat);

    return 0;
}
