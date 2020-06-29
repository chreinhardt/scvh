/*
 * Calculate s(rho, T) for the SCVH EOS on the grid points of REOS3.
 *
 * Author:   Christian Reinhardt
 * Created:  24.06.2020
 * Modified: 29.06.2020
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "scvheos.h"

int main(int argc, char **argv) {
	SCVHEOSMAT *Mat;
    int iMat = SCVHEOS_H;
    double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
    int const nRho = 87; 
    double rho[nRho];
    int const nT = 33;
    double T[nT];
    int i;
    int j;
     
    /* rho axis of REOS3 within the range of SCVH EOS. */
    rho[0] = 3.00000000e-08;
    rho[1] = 4.00000000e-08;
    rho[2] = 5.00000000e-08;
    rho[3] = 6.00000000e-08;
    rho[4] = 7.00000000e-08;
    rho[5] = 8.00000000e-08;
    rho[6] = 9.00000000e-08;
    rho[7] = 1.00000000e-07;
    rho[8] = 2.00000000e-07;
    rho[9] = 3.00000000e-07;
    rho[10] = 4.00000000e-07;
    rho[11] = 5.00000000e-07;
    rho[12] = 6.00000000e-07;
    rho[13] = 7.00000000e-07;
    rho[14] = 8.00000000e-07;
    rho[15] = 9.00000000e-07;
    rho[16] = 1.00000000e-06;
    rho[17] = 2.00000000e-06;
    rho[18] = 3.00000000e-06;
    rho[19] = 4.00000000e-06;
    rho[20] = 5.00000000e-06;
    rho[21] = 6.00000000e-06;
    rho[22] = 7.00000000e-06;
    rho[23] = 8.00000000e-06;
    rho[24] = 9.00000000e-06;
    rho[25] = 1.00000000e-05;
    rho[26] = 2.00000000e-05;
    rho[27] = 3.00000000e-05;
    rho[28] = 4.00000000e-05;
    rho[29] = 5.00000000e-05;
    rho[30] = 6.00000000e-05;
    rho[31] = 7.00000000e-05;
    rho[32] = 8.00000000e-05;
    rho[33] = 9.00000000e-05;
    rho[34] = 1.00000000e-04;
    rho[35] = 2.00000000e-04;
    rho[36] = 3.00000000e-04;
    rho[37] = 4.00000000e-04;
    rho[38] = 5.00000000e-04;
    rho[39] = 6.00000000e-04;
    rho[40] = 7.00000000e-04;
    rho[41] = 8.00000000e-04;
    rho[42] = 9.00000000e-04;
    rho[43] = 1.00000000e-03;
    rho[44] = 2.00000000e-03;
    rho[45] = 3.00000000e-03;
    rho[46] = 4.00000000e-03;
    rho[47] = 5.00000000e-03;
    rho[48] = 6.00000000e-03;
    rho[49] = 7.00000000e-03;
    rho[50] = 8.00000000e-03;
    rho[51] = 9.00000000e-03;
    rho[52] = 1.00000000e-02;
    rho[53] = 2.00000000e-02;
    rho[54] = 3.00000000e-02;
    rho[55] = 4.00000000e-02;
    rho[56] = 5.00000000e-02;
    rho[57] = 6.00000000e-02;
    rho[58] = 7.00000000e-02;
    rho[59] = 8.00000000e-02;
    rho[60] = 9.00000000e-02;
    rho[61] = 1.00000000e-01;
    rho[62] = 2.00000000e-01;
    rho[63] = 3.00000000e-01;
    rho[64] = 4.00000000e-01;
    rho[65] = 5.00000000e-01;
    rho[66] = 6.00000000e-01;
    rho[67] = 7.00000000e-01;
    rho[68] = 8.00000000e-01;
    rho[69] = 9.00000000e-01;
    rho[70] = 1.00000000e+00;
    rho[71] = 1.20000000e+00;
    rho[72] = 1.40000000e+00;
    rho[73] = 1.60000000e+00;
    rho[74] = 1.80000000e+00;
    rho[75] = 2.00000000e+00;
    rho[76] = 2.20000000e+00;
    rho[77] = 3.00000000e+00;
    rho[78] = 4.00000000e+00;
    rho[79] = 5.00000000e+00;
    rho[80] = 9.00000000e+00;
    rho[81] = 1.20000000e+01;
    rho[82] = 1.50000000e+01;
    rho[83] = 2.00000000e+01;
    rho[84] = 3.00000000e+01;
    rho[85] = 5.00000000e+01;
    rho[86] = 7.00000000e+01;
    //rho[87] = 1.00000000e+02;

    /* T axis of REOS3 within the range of SCVH EOS. */
    T[0] = 6.00000000e+01;
    T[1] = 1.00000000e+02;
    T[2] = 2.00000000e+02;
    T[3] = 3.00000000e+02;
    T[4] = 5.00000000e+02;
    T[5] = 7.00000000e+02;
    T[6] = 1.00000000e+03;
    T[7] = 2.00000000e+03;
    T[8] = 3.00000000e+03;
    T[9] = 4.00000000e+03;
    T[10] = 5.00000000e+03;
    T[11] = 8.00000000e+03;
    T[12] = 1.00000000e+04;
    T[13] = 1.20000000e+04;
    T[14] = 1.50000000e+04;
    T[15] = 2.00000000e+04;
    T[16] = 3.00000000e+04;
    T[17] = 5.00000000e+04;
    T[18] = 7.50000000e+04;
    T[19] = 1.00000000e+05;
    T[20] = 1.50000000e+05;
    T[21] = 2.00000000e+05;
    T[22] = 2.50000000e+05;
    T[23] = 3.00000000e+05;
    T[24] = 3.50000000e+05;
    T[25] = 4.00000000e+05;
    T[26] = 4.50000000e+05;
    T[27] = 5.00000000e+05;
    T[28] = 6.00000000e+05;
    T[29] = 7.00000000e+05;
    T[30] = 8.00000000e+05;
    T[31] = 9.00000000e+05;
    T[32] = 1.00000000e+06;

    fprintf(stderr, "SCVH EOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    fprintf(stderr, "Done.\n");
    
#if 0 
    /* Print T axis. */
    printf("T= \n");
    for (i=0; i<nT; i++) {
        printf("%15.7E\n", T[i]);
    }
    printf("\n");
    
    /* Print rho axis. */
    printf("rho= \n");
    for (j=0; j<nRho; j++) {
        printf("%15.7E\n", rho[j]);
    }
    printf("\n");
#endif

    /* Print the entropy. */
    for (i=0; i<nT; i++) {
        for (j=0; j<nRho; j++) {
            printf("%15.7E %15.7E %15.7E\n", rho[j], T[i], scvheosSofRhoT(Mat, rho[j], T[i]));
        }
    }

    scvheosFinalizeMaterial(Mat);

    return 0;
}
