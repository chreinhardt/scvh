/*
 * This program calculates the pressure and entropy for Ravit's models in different units using the
 * SCvH EOS for H-He.
 *
 * Author:   Christian Reinhardt
 * Created:  09.11.2022
 * Modified: 
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "scvheos.h"

struct Model {
    double *R;
    double *rho;
    double *P;
    double *T;
    double *M;
    double *s;
    int nTable;
};

struct Model *ReadModel(char *chFile) {
    struct Model *model;
    FILE *fp;
    char *chLine;
    size_t nCharMax = 256;
    int iRet;
    int nTableMax = 100000;
    int i;

    /* Allocate memory. */
    model = (struct Model *) calloc(1, sizeof(struct Model));
    assert(model != NULL);

    model->R = (double *) calloc(nTableMax, sizeof(double));
    assert(model->R != NULL);
    model->rho = (double *) calloc(nTableMax, sizeof(double));
    assert(model->rho != NULL);
    model->P = (double *) calloc(nTableMax, sizeof(double));
    assert(model->P != NULL);
    model->T = (double *) calloc(nTableMax, sizeof(double));
    assert(model->T != NULL);
    model->M = (double *) calloc(nTableMax, sizeof(double));
    assert(model->M != NULL);
    model->s = (double *) calloc(nTableMax, sizeof(double));
    assert(model->s != NULL);
    model->nTable = 0;

    chLine = (char *) calloc(nCharMax, sizeof(char));

    /* Open the file. */
    fp = fopen(chFile, "r");
    assert(fp != NULL);

    i = 0;
    while (getline(&chLine, &nCharMax, fp) != -1) {
        /* Check if its a comment. */
        if (strchr(chLine, '#') != NULL) continue;

        iRet = sscanf(chLine, "%lf %lf %lf %lf %lf", &model->R[i], &model->rho[i],
                      &model->P[i], &model->T[i],  &model->M[i]);

        /* Check if the number of matches is correct. */
        assert(iRet == 5);
        
        /* Check that the values are sensible. */
        assert(model->R[i] >= 0.0);
        assert(model->rho[i] >= 0.0);
        assert(model->T[i] >= 0.0);

        i++;

        /* Initialize entropy to a weird value. */
        model->s[i] = -1e50;

        assert(i < nTableMax);
    }

    model->nTable = i;

    fprintf(stderr, "ReadModel: read %i lines.\n", model->nTable);

    return model;
}

int main(int argc, char **argv) {
    // SCvH EOS library
    SCVHEOSMAT *Mat;
    double dKpcUnit;
    double dMsolUnit;
    int iMat = 113;
    // Units
    const double MSOLG = 1.99e33;
    const double KPCCM = 3.085678e21;
    double dLUnit;
    double dMUnit;
    // model
    struct Model *model;

#if 0
    /* L_unit = 1 RE, v_unit = 1 km/s. */
    dKpcUnit = 2.06701e-13;
    dMsolUnit = 4.80438e-08;

    /* L_unit = 1 AU, M_unit = 1 MJ. */
    dKpcUnit = 4.84821E-09;
    dMsolUnit = 9.53869E-04;
#endif
    /* cgs */
    dKpcUnit = 0;
    dMsolUnit = 0;

	if (argc != 2) {
		fprintf(stderr,"Usage: calc_ravit_model_codeunits <model.dat>\n");
		exit(1);
	}

    /* Read equilibrium model in cgs. */
    model = ReadModel(argv[1]);

    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);
    
    dLUnit = Mat->dKpcUnit*KPCCM;
    dMUnit = Mat->dMsolUnit*MSOLG;

    /* Calculate entropy. */
    for (int i=0; i<model->nTable; i++) {
        model->s[i] = scvheosSofRhoT(Mat, model->rho[i]/Mat->dGmPerCcUnit, model->T[i]);
        /* Convert to cgs. */
        model->s[i] *= Mat->dErgPerGmUnit;
    }

#if 0	
    /* Print the model. */
    fprintf(stdout, "#%13s%11s%11s%11s%11s\n", "R [cm]", "rho [g/cm^3]", "P [cgs]", "T [K]", "M [g]");

    for (int i=0; i<model->nTable; i++) {
        fprintf(stdout, "%13.3E%11.3E%11.3E%11.3E%11.3E\n", model->R[i], model->rho[i],
                model->P[i], model->T[i], model->M[i]);
    }
#endif
    /* Print the model. */
    fprintf(stdout, "#%14s%15s", "R [cm]", "R");
    fprintf(stdout, "%15s%15s", "rho [cgs]", "rho");
    fprintf(stdout, "%15s%15s", "P [cgs]", "P");
    fprintf(stdout, "%15s", "T [K]");
    fprintf(stdout, "%15s%15s", "M [g]", "M");
    fprintf(stdout, "%15s%15s", "s [cgs]", "s");
    fprintf(stdout, "\n");

    for (int i=0; i<model->nTable; i++) {
        fprintf(stdout, "%15.7E%15.7E", model->R[i], model->R[i]/dLUnit);
        fprintf(stdout, "%15.7E%15.7E", model->rho[i], model->rho[i]/Mat->dGmPerCcUnit);
        fprintf(stdout, "%15.7E%15.7E", model->P[i], model->P[i]/(Mat->dErgPerGmUnit*Mat->dGmPerCcUnit));
        fprintf(stdout, "%15.7E", model->T[i]);
        fprintf(stdout, "%15.7E%15.7E", model->M[i], model->M[i]/dMUnit);
        fprintf(stdout, "%15.7E%15.7E", model->s[i], model->s[i]/Mat->dErgPerGmUnit); 
        fprintf(stdout, "\n");
    }

    /* Free memory. */ 
    scvheosFinalizeMaterial(Mat);	

    free(model);

    return 0;
}
