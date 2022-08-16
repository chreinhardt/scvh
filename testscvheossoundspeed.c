/*
 * A simple program to test the SCVH EOS library.
 *
 * Test if the calculations of the sound speed works correctly.
 *
 * Author: Christian Reinhardt
 * Created: 15.08.2022
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
    fprintf(stderr, "\n");

    /* Calculate the sound speed on the grid points of the EOS table. */
    fprintf(stderr, "Calculate logcs(logrho, logT).\n");
    fp = fopen("testscvheossoundspeed_table.txt", "w");

    fprintf(fp, "# Sound speed logcs(logrho, logT) of the EOS table (nRho = %i nT= %i)\n", Mat->nRho, Mat->nT);
    fprintf(fp, "# Interpolator P: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));
    fprintf(fp, "# Interpolator U: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<Mat->nT; i++) {
        for (int j=0; j<Mat->nRho; j++) {
            fprintf(fp, "%15.7E", gsl_interp2d_get(Mat->InterpLogCs, Mat->dLogCArray, i, j));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

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
    fp = fopen("testscvheossoundspeed_rhoaxis.txt", "w");

    for (int i=0; i<nRho; i++) {
        fprintf(fp, "%15.7E\n", logrhoAxis[i]);
    }

    fclose(fp);

    fp = fopen("testscvheossoundspeed_taxis.txt", "w");

    for (int i=0; i<nT; i++) {
        fprintf(fp, "%15.7E\n", logTAxis[i]);
    }

    fclose(fp);

    /* Calculate logcs(logrho, logT) for interpolated values. */
    fprintf(stderr, "Calculate logcs(logrho, logT).\n");
    fp = fopen("testscvheossoundspeed.txt", "w");

    fprintf(fp, "# Sound speed logcs(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator P: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));
    fprintf(fp, "# Interpolator U: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            double cs = scvheosLogCsofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]);

            if (isfinite(cs)) {
            } else {
                printf("i=%i j=%i: logrho=%15.7E logT=%15.7E logcs= %15.7E\n", i, j, cs, Mat->dLogRhoAxis[j], Mat->dLogTAxis[i]);
                exit(1);
            }
            fprintf(fp, "%15.7E", cs);
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    //if (logrhoAxis) free(logrhoAxis);
    //if (logTAxis) free(logTAxis);

    return 0;
}
