/*
 * Print material data.
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
    SCVHEOSMAT *Mat;
    int iMat;
    //int iMat = SCVHEOS_HHE_LOWRHOT;
    //int iMat = SCVHEOS_HHE;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;

    if (argc != 2) {
        fprintf(stderr, "Usage: scvheosprintmat <iMat>\n");
        exit(1);
    }

    iMat = atoi(argv[1]);

    fprintf(stderr, "SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "\n");
    
    scvheosPrintMat(Mat, stderr);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);
    
    //free(logrhoAxis);
    //free(logTAxis);

    return 0;
}
