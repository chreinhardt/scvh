/*
 * A simple program to test the SCVH EOS library.
 *
 * Test the inversion routines T(rho, u), T(rho, s) and rho(P, T).
 *
 * Author: Christian Reinhardt
 * Created: 16.08.2022
 * Modified: 18.08.2022
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
    // SCvH material
    SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE_LOWRHOT;
    //int iMat = SCVHEOS_HHE;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double *logrhoAxis;
    double *logTAxis;
    int nRho;
    int nT;
    FILE *fp;

    fprintf(stderr, "SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Doing %s interpolation for logP(logrho, logT)\n", gsl_interp2d_name(Mat->InterpLogP));
    fprintf(stderr, "\n");

    nRho = Mat->nRho*10;
    nT = Mat->nT*10;

    /* Generate rho and T axis. */
    logrhoAxis = (double *) calloc(nRho, sizeof(double));
    logTAxis = (double *) calloc(nT, sizeof(double));

    for (int i=0; i<nRho; i++) {
        logrhoAxis[i] = Mat->dLogRhoAxis[0] + i*(Mat->dLogRhoAxis[Mat->nRho-1]-Mat->dLogRhoAxis[0])/(nRho-1);
    }

    for (int i=0; i<nT; i++) {
        logTAxis[i] = Mat->dLogTAxis[0] + i*(Mat->dLogTAxis[Mat->nT-1]-Mat->dLogTAxis[0])/(nT-1);
    }

#if 0
    nRho = (Mat->nRho-1)*2+1;
    nT = (Mat->nT-1)*2+1;

    /* Generate rho and T axis. */
    logrhoAxis = (double *) calloc(nRho, sizeof(double));
    logTAxis = (double *) calloc(nT, sizeof(double));

    for (int i=0; i<Mat->nRho-1; i++) {
        for (int j=0; j<2; j++) {
            logrhoAxis[i*2+j] = Mat->dLogRhoAxis[i] + j*(Mat->dLogRhoAxis[i+1]-Mat->dLogRhoAxis[i])/2;
        }
    }

    logrhoAxis[nRho-1] = Mat->dLogRhoAxis[Mat->nRho-1];

    for (int i=0; i<Mat->nT-1; i++) {
        for (int j=0; j<2; j++) {
            logTAxis[i*2+j] = Mat->dLogTAxis[i] + j*(Mat->dLogTAxis[i+1]-Mat->dLogTAxis[i])/2;
        }
    }

    logTAxis[nT-1] = Mat->dLogTAxis[Mat->nT-1];
#endif

    /* Write both axis to a file. */
    fp = fopen("testscvheosinv_rhoaxis.txt", "w");

    for (int i=0; i<nRho; i++) {
        fprintf(fp, "%15.7E\n", logrhoAxis[i]);
    }

    fclose(fp);

    fp = fopen("testscvheosinv_taxis.txt", "w");

    for (int i=0; i<nT; i++) {
        fprintf(fp, "%15.7E\n", logTAxis[i]);
    }

    fclose(fp);

    /* Calculate the relative error on logT(logrho, logu). */
    fprintf(stderr, "Calculate logT(logrho, logu) on a logrho x logT grid.\n");
    fp = fopen("testscvheosinv_logtoflogrhologu.txt", "w");

    fprintf(fp, "# logT(logrho, logu) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            double logu = scvheosLogUofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]);
            double logT = scvheosLogTofLogRhoLogU(Mat, logrhoAxis[j], logu);
            fprintf(fp, "%15.7E", fabs((logT-logTAxis[i])/logT));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Calculate the relative error on logT(logrho, logs). */
    fprintf(stderr, "Calculate logT(logrho, logs) on a logrho x logT grid.\n");
    fp = fopen("testscvheosinv_logtoflogrhologs.txt", "w");

    fprintf(fp, "# logT(logrho, logs) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogS));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            double logs = scvheosLogSofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]);
            double logT = scvheosLogTofLogRhoLogS(Mat, logrhoAxis[j], logs);
            fprintf(fp, "%15.7E", fabs((logT-logTAxis[i])/logT));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Calculate the relative error on logrho(logP, logT). */
    fprintf(stderr, "Calculate logrho(logP, logT) on a logrho x logT grid.\n");
    fp = fopen("testscvheosinv_logrhooflogplogt.txt", "w");

    fprintf(fp, "# logrho(logP, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            double logP = scvheosLogPofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]);
            double logrho = scvheosLogRhoofLogPLogT(Mat, logP, logTAxis[i]);
            fprintf(fp, "%15.7E", fabs((logrho-logrhoAxis[j])/logrho));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    //free(logrhoAxis);
    //free(logTAxis);

    return 0;
}
