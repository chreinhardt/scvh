/*
 * A simple program to test the SCVH EOS library.
 *
 * Author:   Christian Reinhardt
 * Created:  06.04.2020
 * Modified:
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_H;
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double logrho, logT;
    double rho, T;
    FILE *fp;
    int i, j;

    dKpcUnit = 0.0;
    dMsolUnit = 0.0;

    printf("SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("Done.\n");

    /* Interpolate in rho for differen isotherms. */
    fp = fopen("testscvheospressure.txt", "w");

    for (j=1; j<Mat->nRho-1; j++) {
        logrho = 0.5*(Mat->dLogRhoAxis[j+1] + Mat->dLogRhoAxis[j]);
        fprintf(fp, "%15.7E", logrho);

        for (i=0; i<Mat->nT-1; i++) {
            //logT =  Mat->dLogTAxis[i];
            logT = 0.5*(Mat->dLogTAxis[i+1]+Mat->dLogTAxis[i]);
            fprintf(fp, "%15.7E", scvheosLogPofRhoT(Mat, logrho, logT));
        }
        
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Calculate the pressure at the grid points that are within the range of REOS3. */
    fp = fopen("testscvheospressure_grid.txt", "w");

    for (j=10; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);

        for (i=15; i<Mat->nT-1; i++) {
            T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
            fprintf(fp, "%15.7E", scvheosPofRhoT(Mat, rho, T));
        }
        
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);
    printf("Done.\n");

    return 0;
}
