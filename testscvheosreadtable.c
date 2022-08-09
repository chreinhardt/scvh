/*
 * A simple program to test the SCVH EOS library.
 *
 * Test if the tables are read correctly.
 *
 * Author: Christian Reinhardt
 * Created: 03.04.2020
 * Modified:
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
    SCVHEOSMAT *Mat;
    //int iMat = SCVHEOS_HHE_LOWRHOT;
    int iMat = SCVHEOS_HHE;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    FILE *fp;
    int i, j;

    printf("SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("Done.\n");

    fprintf(stderr, "logT=");
    for (i=0; i<Mat->nT; i++) {
        fprintf(stderr, " %g", Mat->dLogTAxis[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "logRho=");
    for (j=0; j<Mat->nRho; j++) {
        fprintf(stderr, " %g", Mat->dLogRhoAxis[j]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    printf("Writing EOS table to file\n");

    /*
     * Print P(rho) for different T.
     */
    fp = fopen("testscvheos_h_p_rho.txt", "w");

    for (j=0; j<Mat->nRho; j++) {
        fprintf(fp, "%15.7E", Mat->dLogRhoAxis[j]);
        for (i=0; i<Mat->nT; i++) {
            fprintf(fp, "%15.7E", Mat->dLogPArray[j*Mat->nT+i]);
        }
        
        fprintf(fp, "\n");
    }
    fclose(fp);

    /*
     * Print u(rho) for different T.
     */
    fp = fopen("testscvheos_h_u_rho.txt", "w");

    for (j=0; j<Mat->nRho; j++) {
        fprintf(fp, "%15.7E", Mat->dLogRhoAxis[j]);
        for (i=0; i<Mat->nT; i++) {
            fprintf(fp, "%15.7E", Mat->dLogUArray[j*Mat->nT+i]);
        }
        
        fprintf(fp, "\n");
    }
    fclose(fp);

    /*
     * Print s(rho) for different T.
     */
    fp = fopen("testscvheos_h_s_rho.txt", "w");

    for (j=0; j<Mat->nRho; j++) {
        fprintf(fp, "%15.7E", Mat->dLogRhoAxis[j]);
        for (i=0; i<Mat->nT; i++) {
            fprintf(fp, "%15.7E", Mat->dLogSArray[j*Mat->nT+i]);
        }
        
        fprintf(fp, "\n");
    }
    fclose(fp);

    /*
     * Print the EOS table.
     */
    fp = fopen("testscvheosreadtable.txt", "w");

    for (i=0; i<Mat->nT; i++) {
        for (j=0; j<Mat->nRho; j++) {
            fprintf(fp, "%.8f ", Mat->dLogTAxis[i]);
            fprintf(fp, "%.8f ", Mat->dLogRhoAxis[j]);
            fprintf(fp, "%.8f ", gsl_interp2d_get(Mat->InterpLogP, Mat->dLogPArray, i, j));
            fprintf(fp, "%.8f ", gsl_interp2d_get(Mat->InterpLogU, Mat->dLogUArray, i, j));
            fprintf(fp, "%.8f", gsl_interp2d_get(Mat->InterpLogS, Mat->dLogSArray, i, j));
            fprintf(fp, "\n");
        }
    }
    fclose(fp);


    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);
    printf("Done.\n");

    return 0;
}
