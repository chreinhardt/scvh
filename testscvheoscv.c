/*
 * A simple program to test the SCVH EOS library.
 * 
 * Calculate the specific heat capacity from dUdT at constant rho.
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
    //double dKpcUnit = 2.06701e-13;
	//double dMsolUnit = 4.80438e-08;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double rho, T;
    double cv;
    FILE *fp;
    int i, j;

    printf("SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("Done.\n");

    /* Calculate the specific heat capacity at the grid points that are within the range of REOS3. */
    fp = fopen("testscvheoscv_grid.txt", "w");

    for (j=10; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);

        for (i=15; i<Mat->nT-1; i++) {
            T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
            
            // if ((i==0) && (j==0)) printf("rho0= %g T= %g\n", rho, T);
            //cv = scvheosdUdTofRhoT(Mat, rho, T);
            cv = T*scvheosdSdTofRhoT(Mat, rho, T);
            
            cv = (pow(Mat->dLogBase, Mat->dLogUArray[i+1][j])-pow(Mat->dLogBase, Mat->dLogUArray[i][j]))/(pow(Mat->dLogBase, Mat->dLogTAxis[i+1])-pow(Mat->dLogBase, Mat->dLogTAxis[i]));
            if (cv > 0.0) {
            } else {
                //fprintf(stderr, "rho= %15.7E T= %15.7E cv= %15.7E\n", rho, T, cv);
                fprintf(stderr, "rho= %15.7g T= %15.7g cv= %15.7g\n", rho, T, cv);
            }

            fprintf(fp, "%15.7E", cv);
        }
        
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Calculate the difference in the specific heat capacity if it is calculated from u or s. */
    fp = fopen("testscvheoscv_grid_diff.txt", "w");

    for (j=0; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);

        for (i=0; i<Mat->nT-1; i++) {
            T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
            cv = T*scvheosdSdTofRhoT(Mat, rho, T);
    
            if (fabs(cv-scvheosdUdTofRhoT(Mat, rho, T))/cv < 1e-2) {
                fprintf(fp, "%3i", 1);
            } else if (fabs(cv-scvheosdUdTofRhoT(Mat, rho, T))/cv < 1e-1) {
                fprintf(fp, "%3i", 2);
            } else {
                fprintf(fp, "%3i", 3);
            }

            //fprintf(fp, "%15.7E", (cv-scvheosdUdTofRhoT(Mat, rho, T))/cv);
        }
        
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* Mark where cv < 0. */
    fp = fopen("testscvheoscv_grid_neg.txt", "w");

    for (j=10; j<Mat->nRho-1; j++) {
        rho = pow(Mat->dLogBase,  Mat->dLogRhoAxis[j]);

        for (i=15; i<Mat->nT-1; i++) {
            T = pow(Mat->dLogBase,  Mat->dLogTAxis[i]);
            //cv = T*scvheosdSdTofRhoT(Mat, rho, T);
            cv = scvheosdUdTofRhoT(Mat, rho, T);
    
            if (cv > 0.0) {
                fprintf(fp, "%3i", 1);
            } else {
                fprintf(fp, "%3i", 0);
            }

            //fprintf(fp, "%15.7E", (cv-scvheosdUdTofRhoT(Mat, rho, T))/cv);
        }
        
        fprintf(fp, "\n");
    }
    printf("Free memory\n");
    scvheosFinalizeMaterial(Mat);
    printf("Done.\n");

    return 0;
}
