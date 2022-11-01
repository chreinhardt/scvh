/*
 * Calculate u(rho, T) for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  14.10.2022
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
        fprintf(stderr, "Usage: scvheoscalcurhotemp <rho> <T>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    T = atof(argv[2]);

    assert(rho > 0.0);
    assert(T > 0.0);

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    u = scvheosUofRhoT(Mat, rho, T);

    printf("u=%15.7E erg/g (rho= %15.7E g/cm^3 T= %15.7E K)\n", u*Mat->dErgPerGmUnit,rho*Mat->dGmPerCcUnit, T);

    scvheosFinalizeMaterial(Mat);

    return 0;
}
