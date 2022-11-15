/*
 * Calculate P(rho, T) for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  08.11.2022
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
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double rho;
    double T;
    double P;

#if 0
    dKpcUnit = 2.06701e-13;
	dMsolUnit = 4.80438e-08;

    dKpcUnit = 4.84821E-09;
    dMsolUnit = 9.53869E-04;
#endif

    if (argc != 3) {
        fprintf(stderr, "Usage: scvheos_calc_pofrhot <rho> <T>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    T = atof(argv[2]);

    assert(rho > 0.0);
    assert(T > 0.0);

    fprintf(stderr, "SCVH EOS: Initializing material %i (dKpcUnit=%15.7E dMsolUnit=%15.7E)\n",
            iMat, dKpcUnit, dMsolUnit); 

    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    P = scvheosPofRhoT(Mat, rho, T);

    //printf("P=%15.7E erg/cm^3 (rho= %15.7E g/cm^3 T= %15.7E K)\n", P, rho, T);
    printf("P=%15.7E code units (rho= %15.7E code units T= %15.7E K)\n", P, rho, T);
    printf("P=%15.7E erg/cm^3 (rho= %15.7E g/cm^3 T= %15.7E K)\n", P*(Mat->dErgPerGmUnit/Mat->dGmPerCcUnit), rho*Mat->dGmPerCcUnit, T);

    scvheosFinalizeMaterial(Mat);

    return 0;
}
