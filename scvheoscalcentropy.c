/*
 * Calculate s(rho, T) for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  19.06.2020
 * Modified:
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_HHE_EXT_LOWRHOT;
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    double rho, T;

    dKpcUnit = 0.0;
    dMsolUnit = 0.0;

#if 0
    if (argc != 3) {
        fprintf(stderr, "Usage: scvheoscalcentropy <rho> <T>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    T = atof(argv[2]);

    assert(rho > 0.0);
    assert(T > 0.0);

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Done.\n");

    printf("rho= %15.7E g/cm^3 T= %15.7E K s= %15.7E erg/g/K\n", rho, T, scvheosSofRhoT(Mat, rho, T));
#endif

//    double logrho = -11.339134521996131;
//    double logT = 1.6901960800285136;
    double logrho = -7.298432014944073;
    double logT = 3.2227164711475833;


    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Done.\n");

    double logs = scvheosLogSofLogRhoLogT(Mat, logrho, logT);


    printf("logrho= %15.7E g/cm^3 logT= %15.7E K logs= %15.7E erg/g/K\n", logrho, logT, logs);
    scvheosFinalizeMaterial(Mat);

    return 0;
}
