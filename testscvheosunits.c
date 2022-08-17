/*
 * A simple program to test the SCVH EOS library.
 *
 * Test if the values are correctly converted to code units.
 *
 * Author: Christian Reinhardt
 * Created: 17.08.2022
 * Modified:
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
    // SCvH material
    SCVHEOSMAT *Mat1;
    SCVHEOSMAT *Mat2;
    int iMat = SCVHEOS_HHE_LOWRHOT;
    //int iMat = SCVHEOS_HHE;
    /* L_unit = 1 RE, v_unit = 1 km/s. */
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double *logrhoAxis;
    double *logTAxis;
    int nRho;
    int nT;

    fprintf(stderr, "SCVHEOS: Initializing materials %i\n", iMat); 
    Mat1 = scvheosInitMaterial(iMat, 0.0, 0.0);
    Mat2 = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "\n");
    
    assert(Mat1->iMat == Mat2->iMat);

    nRho = Mat1->nRho*10;
    nT = Mat1->nT*10;

    /* Generate rho and T axis. */
    logrhoAxis = (double *) calloc(nRho, sizeof(double));
    logTAxis = (double *) calloc(nT, sizeof(double));

    for (int i=0; i<nRho; i++) {
        logrhoAxis[i] = Mat1->dLogRhoAxis[0] + i*(Mat1->dLogRhoAxis[Mat1->nRho-1]-Mat1->dLogRhoAxis[0])/(nRho-1);
    }

    for (int i=0; i<nT; i++) {
        logTAxis[i] = Mat1->dLogTAxis[0] + i*(Mat1->dLogTAxis[Mat1->nT-1]-Mat1->dLogTAxis[0])/(nT-1);
    }

    /* Check if the values agree. */
    for (int i=0; i<nT; i++) {
        for (int j=0; j<nRho; j++) {
            double rel_err = 1e-10;
            double rho_cgs = pow(10.0, logrhoAxis[j]);
            double T = pow(10.0, logTAxis[i]);

            double P_cgs = scvheosPofRhoT(Mat1, rho_cgs, T);
            double P_codeunits = scvheosPofRhoT(Mat2, rho_cgs/Mat2->dGmPerCcUnit, T);
            double delta_P = (P_cgs-P_codeunits*Mat2->dErgPerGmUnit*Mat2->dGmPerCcUnit)/P_cgs;

            double u_cgs = scvheosUofRhoT(Mat1, rho_cgs, T);
            double u_codeunits = scvheosUofRhoT(Mat2, rho_cgs/Mat2->dGmPerCcUnit, T);
            double delta_u = (u_cgs-u_codeunits*Mat2->dErgPerGmUnit)/u_cgs;

            double s_cgs = scvheosSofRhoT(Mat1, rho_cgs, T);
            double s_codeunits = scvheosSofRhoT(Mat2, rho_cgs/Mat2->dGmPerCcUnit, T);
            double delta_s = (s_cgs-s_codeunits*Mat2->dErgPerGmUnit)/s_cgs;
            
            assert(delta_P <= rel_err);
            assert(delta_u <= rel_err);
            assert(delta_s <= rel_err);
        } 
    }

    /* Free memory. */
    scvheosFinalizeMaterial(Mat1);
    scvheosFinalizeMaterial(Mat2);
    
    //free(logrhoAxis);
    //free(logTAxis);

    return 0;
}
