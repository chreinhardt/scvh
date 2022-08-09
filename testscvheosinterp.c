/*
 * A simple program to test the SCVH EOS library.
 *
 * Test if interpolations works correctly.
 *
 * Author: Christian Reinhardt
 * Created: 09.08.2022
 * Modified:
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

    /* Write both axis to a file. */
    fp = fopen("testscvheosinterp_rhoaxis.txt", "w");

    for (int i=0; i<nRho; i++) {
        fprintf(fp, "%15.7E\n", logrhoAxis[i]);
    }

    fclose(fp);

    fp = fopen("testscvheosinterp_taxis.txt", "w");

    for (int i=0; i<nT; i++) {
        fprintf(fp, "%15.7E\n", logTAxis[i]);
    }

    fclose(fp);

    /* logP(logrho, logT). */
    fprintf(stderr, "Interpolating logP(logrho, logT).\n");
    fp = fopen("testscvheosinterp_logpress.txt", "w");

    fprintf(fp, "# Pressure logP(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosLogPofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* logu(logrho, logT). */
    fprintf(stderr, "Interpolating logU(logrho, logT).\n");
    fp = fopen("testscvheosinterp_intenergy.txt", "w");

    fprintf(fp, "# Internal energy logU(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosLogUofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    free(logrhoAxis);
    free(logTAxis);

    return 0;
}
