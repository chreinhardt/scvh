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
	REOS3MAT *Mat;
    int iMat = SCVHEOS_H;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    FILE *fp;
    int i, j;

    printf("REOS3: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("Done.\n");

    printf("Writing EOS table to file\n");

    /*
     * Print P(rho) for different T.
     */
    fp = fopen("testscvheos_h_p_rho.txt", "w");

    for (j=0; j<Mat->nRho; j++) {
        fprintf(fp, "%15.7E", Mat->dLogRhoAxis[j]);
        for (i=0; i<Mat->nT; i++) {
            fprintf(fp, "%15.7E", Mat->dLogPArray[i][j]);
        }
        
        fprintf(fp, "\n");
    }
    fclose(fp);

    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);
    printf("Done.\n");

    return 0;
}
