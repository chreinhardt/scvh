/*
 * Calculate rho(P, T) for the SCVH EOS.
 *
 * Author:   Christian Reinhardt
 * Created:  15.01.2025
 * Modified:
 */
#include <stdlib.h>
//#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat =  0.0;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    double rho = 0.0;
    double T = 0.0;
    double P = 0.0;

    if (argc != 4) {
        fprintf(stderr, "Usage: scvheos_calc_rhoofpt <P> <T> <iMat>\n");
        exit(1);
    }

    P = atof(argv[1]);
    T = atof(argv[2]); 
    iMat = atoi(argv[3]);

    assert(P > 0.0);
    assert(T > 0.0);

    fprintf(stderr, "SCVH EOS: Initializing material %i (dKpcUnit=%15.7E dMsolUnit=%15.7E)\n",
            iMat, dKpcUnit, dMsolUnit); 

    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    rho = scvheosRhoofPT(Mat, P, T);

    //printf("P=%15.7E erg/cm^3 (rho= %15.7E g/cm^3 T= %15.7E K)\n", P, rho, T);
    printf("rho=%15.7E code units (P= %15.7E code units T= %15.7E K)\n", rho, P, T);
    printf("rho=%15.7E erg/cm^3   (P= %15.7E g/cm^3 T= %15.7E K)\n", rho*Mat->dGmPerCcUnit, P*(Mat->dErgPerGmUnit/Mat->dGmPerCcUnit), T);

    scvheosFinalizeMaterial(Mat);

    return 0;
}
