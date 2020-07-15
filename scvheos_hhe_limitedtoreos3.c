/*
 * Print the pressure and internal energy of a H/He mixture in the region where SCVH EOS and REOS3
 * are defined.
 *
 * Author:   Christian Reinhardt
 * Created:  15.07.2020
 * Modified:
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double rho, T;
    /* Start from log(T) = 2.01 */
    int iTMin = 19;
    /* Start from log(rho) = -7.45 */
    int iRhoMin = 11;
    FILE *fp;
    int i, j;

    printf("SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("Done.\n");

    /* Print P(rho, T) on the grid points. */
    fp = fopen("scvheos_hhe_pressure.txt", "w");

    for (j=iRhoMin; j<Mat->nRho; j++) {
        rho = pow(10.0, Mat->dLogRhoAxis[j]);
        fprintf(fp, "%15.7E", rho);
        for (i=iTMin; i<Mat->nT-1; i++) {
            T = pow(10.0,  Mat->dLogTAxis[i]);
            fprintf(fp, "%15.7E", pow(10.0, Mat->dLogPArray[i][j]));
            //fprintf(fp, "%15.7E", scvheosPofRhoT(Mat, rho, T));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Print u(rho, T) on the grid points. */
    fp = fopen("scvheos_hhe_intenergy.txt", "w");

    for (j=iRhoMin; j<Mat->nRho; j++) {
        rho = pow(10.0, Mat->dLogRhoAxis[j]);
        fprintf(fp, "%15.7E", rho);
        for (i=iTMin; i<Mat->nT-1; i++) {
            T = pow(10.0,  Mat->dLogTAxis[i]);
            fprintf(fp, "%15.7E", pow(10.0, Mat->dLogUArray[i][j]));
            //fprintf(fp, "%15.7E", scvheosUofRhoT(Mat, rho, T));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);

    return 0;
}
