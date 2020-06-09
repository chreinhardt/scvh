/*
 * A simple program to test the SCVH EOS library.
 *
 * Author:   Christian Reinhardt
 * Created:  08.06.2020
 * Modified: 09.06.2020
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
    double rho, T;
    FILE *fp;
    int i, j;

    /* Use the original units. */
    dKpcUnit = 0.0;
    dMsolUnit = 0.0;

    printf("SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("Done.\n");
    
	/* Write the (rho, T) within the range of REOS3 to file. */
    fp = fopen("testscvheosderivs_rhoT_grid.txt", "w");

    for (j=10; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);
		fprintf(fp, "%15.7E", rho);
    	fprintf(fp, "\n");
    }
	fprintf(fp, "\n");
	for (i=15; i<Mat->nT-1; i++) {
		T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
		fprintf(fp, "%15.7E", T);
		fprintf(fp, "\n");
	}

    fclose(fp);

    /* Calculate dP/drho(rho, T) on the grid points of the EOS table. */	
    fp = fopen("testscvheosderivs_dpdrhoofrhot.txt", "w");

    for (j=0; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase, Mat->dLogRhoAxis[j]);

        for (i=0; i<Mat->nT-1; i++) {
            T = pow(Mat->dLogBase, Mat->dLogTAxis[i]);

            //fprintf(stderr, "rho= %g T= %g dPdRho= %g\n", rho, T, scvheosdPdRhoofRhoT(Mat, rho, T));
            fprintf(fp, "%15.7E", scvheosdPdRhoofRhoT(Mat, rho, T));
        }
        
        fprintf(fp, "\n");
    }
	
    fclose(fp);

    /* Do the same but only for the grid points that are within the range of REOS3. */
    fp = fopen("testscvheosderivs_dpdrhoofrhot_reos3.txt", "w");

    for (j=10; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);

        for (i=15; i<Mat->nT-1; i++) {
			T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
            fprintf(fp, "%15.7E", scvheosPofRhoT(Mat, rho, T));
        }
        fprintf(fp, "\n");
    }

    /* Calculate dP/dT(rho, T) on the grid points of the EOS table. */	
    fp = fopen("testscvheosderivs_dpdtofrhot.txt", "w");

    for (j=0; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase, Mat->dLogRhoAxis[j]);

        for (i=0; i<Mat->nT-1; i++) {
            T = pow(Mat->dLogBase, Mat->dLogTAxis[i]);

            fprintf(fp, "%15.7E", scvheosdPdTofRhoT(Mat, rho, T));
        }
        
        fprintf(fp, "\n");
    }
	
    fclose(fp);

    /* Do the same but only for the grid points that are within the range of REOS3. */
    fp = fopen("testscvheosderivs_dpdtofrhot_reos3.txt", "w");

    for (j=10; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);

        for (i=15; i<Mat->nT-1; i++) {
			T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
            fprintf(fp, "%15.7E", scvheosdPdTofRhoT(Mat, rho, T));
        }
        fprintf(fp, "\n");
    }

    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);
    printf("Done.\n");

    return 0;
}
