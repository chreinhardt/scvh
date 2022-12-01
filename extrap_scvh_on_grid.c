/*
 * Extrapolate the SCvH EOS table for H-He on a different grid to check if details how the table is
 * extended make a significant difference.
 *
 * Author: Christian Reinhardt
 * Created: 30.11.2022
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
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double *logrho;
    double *logT;
    double logrho_min = -14.0;
    double logrho_max = -4.0;
    double logT_min = 1.06;
    double logT_max = 3.46;
    int nRho = 201;
    int nT = 31;
    FILE *fp;

    fprintf(stderr, "SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    /* Generate rho and T axis. */
    logrho = (double *) calloc(nRho, sizeof(double));
    logT = (double *) calloc(nT, sizeof(double));

    for (int i=0; i<nRho; i++) {
        logrho[i] = logrho_min + i*(logrho_max-logrho_min)/(nRho-1);
    }

    for (int i=0; i<nT; i++) {
        logT[i] = logT_min + i*(logT_max-logT_min)/(nT-1);
    }

    /* Write an EOS table for the extrapolated values. */
    fp = fopen("extrap_scvh_on_grid.txt", "w");

    fprintf(fp, "# logT [K] logRho [g/cc] logP [barye] logE [erg/g] logS [erg/g/K]\n");

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            double logP = scvheosLogPofLogRhoLogT(Mat, logrho[j], logT[i]);
            double logu = scvheosLogUofLogRhoLogT(Mat, logrho[j], logT[i]);
            double logs = scvheosLogSofLogRhoLogT(Mat, logrho[j], logT[i]);
            
            fprintf(fp, "%.8e %.8e %.8e %.8e %.8e\n", logT[i], logrho[j], logP, logu, logs);
        }
    }

    fclose(fp);
#if 0


    /* Write both axis to a file. */
    fp = fopen("extrap_scvh_on_grid_rhoaxis.txt", "w");

    for (int i=0; i<nRho; i++) {
        fprintf(fp, "%15.7E\n", logrho[i]);
    }

    fclose(fp);

    fp = fopen("extrap_scvh_on_grid_taxis.txt", "w");

    for (int i=0; i<nT; i++) {
        fprintf(fp, "%15.7E\n", logT[i]);
    }

    fclose(fp);

    /* logP(logrho, logT). */
    fprintf(stderr, "Interpolating logP(logrho, logT).\n");
    fp = fopen("testscvheosinterp_logpress.txt", "w");

    fprintf(fp, "# Pressure logP(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogP));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosLogPofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);

    /* logu(logrho, logT). */
    fprintf(stderr, "Interpolating logU(logrho, logT).\n");
    fp = fopen("testscvheosinterp_logintenergy.txt", "w");

    fprintf(fp, "# Internal energy logU(logrho, logT) (nRho = %i nT= %i)\n", nRho, nT);
    fprintf(fp, "# Interpolator: %s (GSL)\n", gsl_interp2d_name(Mat->InterpLogU));

    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            fprintf(fp, "%15.7E", scvheosLogUofLogRhoLogT(Mat, logrhoAxis[j], logTAxis[i]));
        } 
        fprintf(fp, "\n");
    }

    fclose(fp);
#endif

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    free(logrho);
    free(logT);

    return 0;
}
