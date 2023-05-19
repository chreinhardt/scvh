/*
 * A simple program to test the SCVH EOS library.
 *
 * Calculate T(rho) along different isentropes.
 *
 * Author: Christian Reinhardt
 * Created: 04.04.2023
 * Modified:
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "scvheos.h"

/*
 * This function calculates derivatives for the GSL ode solver.
 */
int derivs(double rho, const double u[], double dudrho[], void *params) {
    SCVHEOSMAT *Mat = (SCVHEOSMAT *) params;
    double P;

    /* Calculate pressure. */
    P = scvheosPofRhoU(Mat, rho, u[0]);

    if (P < 0.0) {
        return GSL_EBADFUNC;
    }

    dudrho[0] = P/(rho*rho);

    return GSL_SUCCESS;
}

int main(int argc, char **argv) {
    // SCvH material
    SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE_EXT_LOWRHOT;
    //int iMat = SCVHEOS_HHE;
    //double dKpcUnit = 0.0;
	//double dMsolUnit = 0.0;
    double dKpcUnit = 4.8482100E-09;
	double dMsolUnit = 9.5386900E-04;
    double *logrhoAxis;
    //double *logs;
    double LogRhoMin;
    double LogRhoMax;
    double LogTMin;
    int nRho;
    /* ODE system. */
    size_t ndim = 1;
    //gsl_odeiv2_system sys = {derivs, NULL, ndim, Mat};
    gsl_odeiv2_system sys;
    /* Driver. */
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
    double eps_abs = 0.0;
    double eps_rel = 1e-8;
    double h = 1e-10;
    //gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, T, h, eps_abs, eps_rel);
    gsl_odeiv2_driver *d;
    FILE *fp;

    fprintf(stderr, "SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Doing %s interpolation for logP(logrho, logT)\n", gsl_interp2d_name(Mat->InterpLogP));
    fprintf(stderr, "\n");
    
    /* Number of data points along one isentrope. */
    nRho = 1000;

#if 0    
    LogRhoMin = Mat->dLogRhoAxis[0];
    LogRhoMax = Mat->dLogRhoAxis[Mat->nRho-1];

    LogTMin = Mat->dLogTAxis[0];
#endif

#if 0
    LogRhoMin = log10(4.5800000E-12);
    LogRhoMax = log10(5.0300000E-08);

    LogTMin = log10(4.9000000E+01);
#endif

#if 0
    /* 3MJ model of Ravit starting from the surface. */
    LogRhoMin = log10(4.6564849E-03);
    LogRhoMax = log10(9.1895024E-01);
    LogTMin = log10(2.2900000E+01);
#endif

#if 0
    /* 3MJ model of Ravit starting from the core. */
    LogRhoMin = log10(9.1895024E-01);
    LogRhoMax = log10(4.6564849E-03);
    LogTMin = log10(2.0700000E+02);
#endif

    /* 10MJ (late) model of Ravit starting from the surface. */
    LogRhoMin = log10(8.0782957E-03);
    LogRhoMax = log10(8.8720148E+01);
    LogTMin = log10(4.9000000E+01);

    /* Fix sign of the step size if nescessary. */
    if (LogRhoMin > LogRhoMax) h *= (-1);

    /* Generate rho axis and calculate entropy for each isentrope. */
    logrhoAxis = (double *) calloc(nRho, sizeof(double));
    //logs = (double *) calloc(nT, sizeof(double));


    for (int i=0; i<nRho; i++) {
        logrhoAxis[i] = LogRhoMin + i*(LogRhoMax-LogRhoMin)/(nRho-1);
    }

    /* Initialize entropy and internal energy. */
    double logs = scvheosLogSofLogRhoLogT(Mat, logrhoAxis[0], LogTMin);
    
    double rho = pow(10.0, logrhoAxis[0]);
    double u = scvheosUofRhoT(Mat, rho, pow(10.0, LogTMin));

    //gsl_odeiv2_system sys = {derivs, NULL, ndim, Mat};
    sys.function = derivs;
    sys.jacobian = NULL;
    sys.dimension = ndim;
    sys.params = Mat;

    d = gsl_odeiv2_driver_alloc_y_new(&sys, T, h, eps_abs, eps_rel);

    /* Calculate T(rho, s=const) for different entropies. */
    fp = fopen("calc_isentrope_logrhologt.txt", "w");

    fprintf(fp, "# logT(logrho, logs) (nRho = %i)\n", nRho);
    fprintf(fp, "# logs=%15.7E\n", logs);
    fprintf(fp, "# dKpcUnit=%15.7E dMsolUnit=%15.7E\n", Mat->dKpcUnit, Mat->dMsolUnit);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogS));
    fprintf(fp, "#%14s%15s%15s%15s%15s%15s\n", "logrho", "logT", "logs", "logs_int", "err", "logT_PdV");

    fprintf(fp, "%15.7E%15.7E%15.7E%15.7E%15.7E%15.7E\n", logrhoAxis[0], LogTMin, logs, logs, 0.0, log10(scvheosTofRhoU(Mat, rho, u)));
 
    for (int i=1; i<nRho; i++) {
        /* Calculate T(rho, s). */
        double logT = scvheosLogTofLogRhoLogS(Mat, logrhoAxis[i], logs);
        double logs_int = scvheosLogSofLogRhoLogT(Mat, logrhoAxis[i], logT); 

        /* Calculate T from the first law of thermodynamics. */
        int status = gsl_odeiv2_driver_apply(d, &rho, pow(10.0, logrhoAxis[i]), &u);

        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Error: status= %i\n", status);
        }

        fprintf(fp, "%15.7E%15.7E%15.7E%15.7E%15.7E%15.7E\n", logrhoAxis[i], logT, logs, logs_int, fabs((logs-logs_int)/logs), log10(scvheosTofRhoU(Mat, rho, u)));
    }

    fclose(fp);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    //free(logrhoAxis);
    //free(logTAxis);

    return 0;
}
