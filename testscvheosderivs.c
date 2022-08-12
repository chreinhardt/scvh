/*
 * A simple program to test the SCVH EOS library.
 *
 * Test if the derivatives are correct.
 *
 * Author: Christian Reinhardt
 * Created: 10.08.2022
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

#if 0
    // CR: Try calculating derivatives on the grid points of the eos table
    nRho = Mat->nRho;
    nT = Mat->nT;
    logrhoAxis = Mat->dLogRhoAxis;
    logTAxis = Mat->dLogTAxis;
#endif

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

    /* dlogP/dlogrho(logrho, logT). */
    fprintf(stderr, "Calculate dlogP/dlogrho(logrho, logT).\n");
    fp = fopen("testscvheosinterp_dlogpdlogrho.txt", "w");

    fprintf(fp, "# dlogP/dlogrho(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosdLogPdLogRhoofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* dlogP/dlogT(logrho, logT). */
    fprintf(stderr, "Calculate dlogP/dlogT(logrho, logT).\n");
    fp = fopen("testscvheosinterp_dlogpdlogt.txt", "w");

    fprintf(fp, "# dlogP/dlogT(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosdLogPdLogTofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* dlogu/dlogrho(logrho, logT). */
    fprintf(stderr, "Calculate dlogu/dlogrho(logrho, logT).\n");
    fp = fopen("testscvheosinterp_dlogudlogrho.txt", "w");

    fprintf(fp, "# dlogu/dlogrho(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosdLogUdLogRhoofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* dlogu/dlogT(logrho, logT). */
    fprintf(stderr, "Calculate dlogu/dlogT(logrho, logT).\n");
    fp = fopen("testscvheosinterp_dlogudlogt.txt", "w");

    fprintf(fp, "# dlogu/dlogT(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosdLogUdLogTofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* dlogs/dlogrho(logrho, logT). */
    fprintf(stderr, "Calculate dlogs/dlogrho(logrho, logT).\n");
    fp = fopen("testscvheosinterp_dlogsdlogrho.txt", "w");

    fprintf(fp, "# dlogs/dlogrho(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogS));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosdLogSdLogRhoofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* dlogs/dlogT(logrho, logT). */
    fprintf(stderr, "Calculate dlogs/dlogT(logrho, logT).\n");
    fp = fopen("testscvheosinterp_dlogsdlogt.txt", "w");

    fprintf(fp, "# dlogs/dlogT(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogS));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosdLogSdLogTofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
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
