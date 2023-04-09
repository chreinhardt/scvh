/*
 * A simple program to test the SCVH EOS library.
 *
 * Calculate T(rho) along different isentropes.
 *
 * Author: Christian Reinhardt
 * Created: 04.04.2023
 * Modified:
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
    // SCvH material
    SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE_EXT_LOWRHOT;
    //int iMat = SCVHEOS_HHE;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double *logrhoAxis;
    //double *logs;
    double LogRhoMin;
    double LogRhoMax;
    double LogTMin;
    int nRho;
    FILE *fp;

    fprintf(stderr, "SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Doing %s interpolation for logP(logrho, logT)\n", gsl_interp2d_name(Mat->InterpLogP));
    fprintf(stderr, "\n");
    
    /* Number of data points along one isentrope. */
    nRho = 100;
    
    LogRhoMin = Mat->dLogRhoAxis[0];
    LogRhoMax = Mat->dLogRhoAxis[Mat->nRho-1];

    LogTMin = Mat->dLogTAxis[0];

    /* Generate rho axis and calculate entropy for each isentrope. */
    logrhoAxis = (double *) calloc(nRho, sizeof(double));
    //logs = (double *) calloc(nT, sizeof(double));


    for (int i=0; i<nRho; i++) {
        logrhoAxis[i] = LogRhoMin + i*(LogRhoMax-LogRhoMin)/(nRho-1);
    }

    double logs = scvheosLogSofLogRhoLogT(Mat, logrhoAxis[0], LogTMin);
    
    /* Calculate T(rho, s=const) for different entropies. */
    fp = fopen("calc_isentrope_logrhologt.txt", "w");

    fprintf(fp, "# logT(logrho, logs) (nRho = %i)\n", nRho);
    fprintf(fp, "# logs=%15.7E\n", logs);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogS));

    for (int i=0; i<nRho; i++) {
        double logT = scvheosLogTofLogRhoLogS(Mat, logrhoAxis[i], logs);
        
        fprintf(fp, "%15.7E%15.7E%15.7E\n", logrhoAxis[i], logT, logs);
    }

    fclose(fp);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    //free(logrhoAxis);
    //free(logTAxis);

    return 0;
}
