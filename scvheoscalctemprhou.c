/*
 * Calculate T(rho, u) for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  18.08.2022
 * Modified:
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE_LOWRHOT;
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double rho;
    double u;
    double T;

    //dKpcUnit = 0.0;
    //dMsolUnit = 0.0;
    dKpcUnit = 4.84821E-09;
    dMsolUnit = 9.53869E-04;


    if (argc != 3) {
        fprintf(stderr, "Usage: scvheoscalctemprhou <rho> <u>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    u = atof(argv[2]);

    assert(rho > 0.0);
    assert(u > 0.0);

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    T = scvheosTofRhoU(Mat, rho, u);

    printf("T=%15.7E K (rho= %15.7E g/cm^3 u= %15.7E erg/g)\n", T, rho*Mat->dGmPerCcUnit, u*Mat->dErgPerGmUnit);

    scvheosFinalizeMaterial(Mat);

    return 0;
}
