/*
 * Author:   Christian Reinhardt
 * Created:  31.10.2022
 * Modified: 
 *
 * This program calculates the entropy of a model generated with ballic.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "scvheos.h"

struct BallicModel {
    double *R;
    double *rho;
    double *u;
    double *M;
    int *iMat;
    double *P;
    double *T;
    double *s;
    int nTable;
};

struct BallicModel *ReadBallicModel(char *chFile) {
    struct BallicModel *model;
    FILE *fp;
    char *chLine;
    size_t nCharMax = 256;
    int iRet;
    int nTableMax = 100000;
    int i;

    /* Allocate memory. */
    model = (struct BallicModel *) calloc(1, sizeof(struct BallicModel));
    assert(model != NULL);

    model->R = (double *) calloc(nTableMax, sizeof(double));
    assert(model->R != NULL);
    model->rho = (double *) calloc(nTableMax, sizeof(double));
    assert(model->rho != NULL);
    model->u = (double *) calloc(nTableMax, sizeof(double));
    assert(model->u != NULL);
    model->M = (double *) calloc(nTableMax, sizeof(double));
    assert(model->M != NULL);
    model->iMat = (int *) calloc(nTableMax, sizeof(int));
    assert(model->iMat != NULL);
    model->P = (double *) calloc(nTableMax, sizeof(double));
    assert(model->P != NULL);
    model->T = (double *) calloc(nTableMax, sizeof(double));
    assert(model->T != NULL);
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

        iRet = sscanf(chLine, "%lf %lf %lf %lf %d %lf %lf", &model->R[i], &model->rho[i],
                      &model->M[i], &model->u[i],  &model->iMat[i], &model->P[i], &model->T[i]);

        /* Check if the number of matches is correct. */
        assert(iRet == 7);
        
        /* Check that the values are sensible. */
        assert(model->R[i] >= 0.0);
        assert(model->rho[i] >= 0.0);

        i++;

        /* Initialize entropy to a weird value. */
        model->s[i] = -1e50;

        assert(i < nTableMax);
    }

    model->nTable = i;

    fprintf(stderr, "ReadBallicModel: read %i lines.\n", model->nTable);

    return model;
}

int main(int argc, char **argv) {
    // SCvH EOS library
    SCVHEOSMAT *Mat;
    double dKpcUnit;
    double dMsolUnit;
    int iMat;
    // ballic model
    struct BallicModel *model;

    /* L_unit = 1 RE, v_unit = 1 km/s. */
    dKpcUnit = 2.06701e-13;
    dMsolUnit = 4.80438e-08;

    /* L_unit = 1 AU, M_unit = 1 MJ. */
    dKpcUnit = 4.84821E-09;
    dMsolUnit = 9.53869E-04;


	if (argc != 3) {
		fprintf(stderr,"Usage: calc_model_entropy <ballic.model> <iMat>\n");
		exit(1);
	}

    /* Read equilibrium model. */
    model = ReadBallicModel(argv[1]);

    iMat = atoi(argv[2]);

    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    for (int i=0; i<model->nTable; i++) {
        model->s[i] = scvheosSofRhoT(Mat, model->rho[i], model->T[i]);
    }
	
    /* Print the model. */
    fprintf(stdout, "#%14s%15s%15s%15s%4s%15s%15s%15s\n", "R", "rho", "M", "u", "mat", "P", "T", "s");

    for (int i=0; i<model->nTable; i++) {
        fprintf(stdout, "%15.7E%15.7E%15.7E%15.7E%4i%15.7E%15.7E%15.7E\n", model->R[i], model->rho[i],
                model->M[i], model->u[i], model->iMat[i], model->P[i], model->T[i], model->s[i]);
    }

    /* Free memory. */ 
    scvheosFinalizeMaterial(Mat);	

    free(model);

    return 0;
}
