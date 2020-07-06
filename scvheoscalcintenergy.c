/*
 * Calculate u(rho, T) for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  06.07.2020
 * Modified:
 */
#include <stdlib.h>
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

    dKpcUnit = 0.0;
    dMsolUnit = 0.0;

    if (argc != 3) {
        fprintf(stderr, "Usage: scvheoscalcintenergy <rho> <T>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    T = atof(argv[2]);

    assert(rho > 0.0);
    assert(T > 0.0);

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    printf("rho= %15.7E g/cm^3 T= %15.7E K u= %15.7E erg/g\n", rho, T, scvheosUofRhoT(Mat, rho, T));

    scvheosFinalizeMaterial(Mat);

    return 0;
}
