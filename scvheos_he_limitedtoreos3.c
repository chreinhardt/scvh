/*
 * Print the pressure and internal energy of helium in the region where SCVH EOS and REOS3 are
 * defined.
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
    int iMat = SCVHEOS_HE;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    /* Start from log(T) = 2.01 */
    int iTMin = 19;
    /* Start from log(rho) = -7.45 */
    int iRhoMin = 11;
    FILE *fp;
    int i, j;

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    
    fprintf(stderr, "logrho_min= %15.7E logrho_max= %15.7E\n", Mat->dLogRhoAxis[iRhoMin],
            Mat->dLogRhoAxis[Mat->nRho-1]);

    fprintf(stderr, "logT_min=   %15.7E logT_max=   %15.7E\n", Mat->dLogTAxis[iTMin],
            Mat->dLogTAxis[Mat->nT-1]);

    /* Print P(rho, T) on the grid points. */
    fp = fopen("scvheos_he_pressure.txt", "w");

    fprintf(fp, "# logrho logP (iMat= %i, dKpcUnit= %15.7E dMsolUnit= %15.7E)\n", Mat->iMat, dKpcUnit,
            dMsolUnit);
    for (j=iRhoMin; j<Mat->nRho; j++) {
        fprintf(fp, "%15.7E", Mat->dLogRhoAxis[j]);
        for (i=iTMin; i<Mat->nT; i++) {
            fprintf(fp, "%15.7E", Mat->dLogPArray[i][j]);
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Print u(rho, T) on the grid points. */
    fp = fopen("scvheos_he_intenergy.txt", "w");

    fprintf(fp, "# logrho logu (iMat= %i, dKpcUnit= %15.7E dMsolUnit= %15.7E)\n", Mat->iMat, dKpcUnit,
            dMsolUnit);
    for (j=iRhoMin; j<Mat->nRho; j++) {
        fprintf(fp, "%15.7E", Mat->dLogRhoAxis[j]);
        for (i=iTMin; i<Mat->nT; i++) {
            fprintf(fp, "%15.7E", Mat->dLogUArray[i][j]);
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);

    return 0;
}
