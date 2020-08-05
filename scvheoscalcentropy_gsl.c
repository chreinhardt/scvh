/*
 * Calculate the entropy for SCVH EOS from Miguel et al (2016).
 *
 * Author:   Christian Reinhardt
 * Date:     03.08.2020
 * Modified:  
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "scvheos.h"

struct rhoFuncEntropyParams {
    SCVHEOSMAT *Mat;
    double T0;
};

struct TFuncEntropyParams {
    SCVHEOSMAT *Mat;
    double rho;
};

/* 
 * Calculate f_rho(rho) = P(T0, rho)/T0/rho^2. 
 */
double rhoFuncEntropy(double rho, void *params) {
    struct rhoFuncEntropyParams *p;
    SCVHEOSMAT *Mat;
    double T0;

    p = (struct rhoFuncEntropyParams *) params;
    
    Mat = p->Mat;
    T0 = p->T0;

    return scvheosPofRhoT(Mat, rho, T0)/T0/(rho*rho);
}

/*
 * Calculate f_T(T) = u(rho, T)/T^2.
 */
double TFuncEntropy(double T, void *params) {
    struct TFuncEntropyParams *p;
    SCVHEOSMAT *Mat;
    double rho;

    p = (struct TFuncEntropyParams *) params;
    
    Mat = p->Mat;
    rho = p->rho;

    return scvheosUofRhoT(Mat, rho, T)/(T*T);
}

double SolveIntegralRho(SCVHEOSMAT *Mat, double rho0, double rho1, double T0) {
    /* GSL Integrator. */
    gsl_integration_workspace *w;
    gsl_function F;
    struct rhoFuncEntropyParams Params;
    /* Other variables. */
    double result;
    double error;
    const double err_abs = 0.0;
    const double err_rel = 1e-4;

    /* Initialize the parameters. */
    Params.Mat = Mat;
    Params.T0 = T0;

    /* Initialize the function to be integrated. */
    F.function = &rhoFuncEntropy;
    F.params = &Params;

    /* Initialize integrator. */
    w = gsl_integration_workspace_alloc(1000000);

    //gsl_integration_qags(&F, rho0, rho1, err_abs, err_rel, 1000000, w, &result, &error);
    gsl_integration_qag(&F, rho0, rho1, err_abs, err_rel, 1000000, GSL_INTEG_GAUSS61, w, &result, &error);

    gsl_integration_workspace_free(w);
    
    return result;
}


double SolveIntegralT(SCVHEOSMAT *Mat, double T0, double T1, double rho) {
    /* GSL Integrator. */
    gsl_integration_workspace *w;
    gsl_function F;
    struct TFuncEntropyParams Params;
    /* Other variables. */
    double result;
    double error;
    const double err_abs = 0.0;
    const double err_rel = 1e-4;

    /* Initialize the parameters. */
    Params.Mat = Mat;
    Params.rho = rho;

    /* Initialize the function to be integrated. */
    F.function = &TFuncEntropy;
    F.params = &Params;

    /* Initialize integrator. */
    w = gsl_integration_workspace_alloc(1000000);

    //gsl_integration_qags(&F, T0, T1, err_abs, err_rel, 1000000, w, &result, &error);
    gsl_integration_qag(&F, T0, T1, err_abs, err_rel, 1000000, GSL_INTEG_GAUSS61, w, &result, &error);

    gsl_integration_workspace_free(w);
    
    return result;
}

double SolveIntegral(SCVHEOSMAT *Mat, double rho, double T, double rho0, double T0, double s0) {
    double s;

    s = scvheosUofRhoT(Mat, rho, T)/T - SolveIntegralRho(Mat, rho0, rho, T0) + SolveIntegralT(Mat, T0, T, rho) + s0;

    s -= scvheosUofRhoT(Mat, rho0, T0);
    return s; 
}

#if 0
/*
 * Calculate the entropy at each grid point of the EOS table.
 */
int scvheosCalcEntropy(SCVHEOSMAT *Mat, double *s, double rho0, double T0, double s0) {
    int i, j;

    fprintf(stderr, "rho0= %g T0= %g s0= %g\n", rho0, T0, s0);
    fprintf(stderr, "rho= %g T= %g\n",  Mat->rhoAxis[0], Mat->TAxis[0]);
    /* Calculate the entropy of the first grid point. */
    s[0] = SolveIntegral(Mat, Mat->rhoAxis[0], Mat->TAxis[0], rho0, T0, s0);
    
    fprintf(stderr, "rho= %g T= %g s= %g\n",  Mat->rhoAxis[0], Mat->TAxis[0], s0);

    /* Calculate the entropy for all temperatures at rho=rho_min. */
    for (i=1; i<Mat->nT; i++) {
        s[0*Mat->nT+i] = SolveIntegral(Mat, Mat->rhoAxis[0], Mat->TAxis[i], Mat->rhoAxis[0], Mat->TAxis[i-1], s[0*Mat->nT+(i-1)]);
        //s[0*Mat->nT+i] += reos3UofRhoT(Mat,Mat->rhoAxis[0], Mat->TAxis[i-1]);
        fprintf(stderr, "i=%i: s0= %15.7E s1= %15.7E\n", i, s[0*Mat->nT+(i-1)], s[0*Mat->nT+i]);
    }

    /* Now calculate the entropy for each isotherm. */
    for (i=0; i<Mat->nT; i++) {
        for (j=1; j<Mat->nRho; j++) {
            s[j*Mat->nT+i] = SolveIntegral(Mat, Mat->rhoAxis[j], Mat->TAxis[i], Mat->rhoAxis[j-1], Mat->TAxis[i], s[(j-1)*Mat->nT+i]);
            //s[j*Mat->nT+i] += reos3UofRhoT(Mat,Mat->rhoAxis[j-1], Mat->TAxis[i]);
        }
    }

    return REOS3_SUCCESS;
}
#endif

int main(int argc, char **argv) {
    // REOS3 material
    SCVHEOSMAT *Mat;
    //int iMat = REOS3_H_F;
    int iMat = SCVHEOS_H;
    double dKpcUnit = 0.0;
    double dMsolUnit = 0.0;
    double rho;
    double T;
    double *s;
    double rho0;
    double T0;
    double s0;
    FILE *fp;
    int i, j;

    fprintf(stderr, "SCVHEOS: Initializing material %i\n", iMat); 
    Mat = scvheosInitMaterial(iMat, dKpcUnit, dMsolUnit);

    /* Do not use the standard GSL error handler. */
    gsl_set_error_handler_off();

    /* Allocate memory as a 1D array like in GSL. */
    s = (double *) calloc(Mat->nT*Mat->nRho, sizeof(double));
    assert(s != NULL);
 
    /* Reference point from Miguel et al. (2016). */
    T0 = 29500.0;
    rho0 = 0.036;

    /* Reference point from Miguel et al. (2018) as suggested by Yamila. */
    T0 = 150.0;
    rho0 = 0.3;
    s0 = 0.0;

    /* Reference point from Nadine Nettelmann. */
    T0 = 1000.0;
    rho0 = 0.01;
    //s0 = 1.03783e+11;
    
    rho = rho0;
    T = T0;
    
    s0 = 1.044139746e+11;

    printf("rho= %g [g/cm^3] T= %g [K] s= %g [erg/g/K] (rho0= %g [g/cm^3] T0= %g [K] s0= %g [erg/g/K])\n", rho, T, SolveIntegral(Mat, rho, T, rho0, T0, s0), rho0, T0, s0);

    //exit(1); 

    /*
     * Calculate the entropy on each grid point from
     *
     * s(T_i, rho_j) = u(T_i, rho_j)/T_i - int(P(T0, rho)/T0*1/(rho^2)drho) + int(u(T, rho_j)/T^2 dT + s0
     *
     * where the first integral is from rho0 to rho_j and the second from T0 to T_i.
     */
    fprintf(stderr, "rho0= %15.7E T0= %15.7E s0= %15.7E\n", rho0, T0, s0);

#if 0
    for (j=0; j<Mat->nRho; j++) {
        for (i=0; i<Mat->nT; i++) {
            rho = Mat->rhoAxis[j];
            T = Mat->TAxis[i];
            //fprintf(stderr, "rho= %g T= %g\n", rho, T);
            s[j*Mat->nT+i] = SolveIntegral(Mat, rho, T, rho0, T0, s0);
        }
    }
#endif

#if 0
    if (scvheosCalcEntropy(Mat, s, rho0, T0, s0) != REOS3_SUCCESS) {
        fprintf(stderr, "reos3CalcEntropy: Could not calculate the entropy.\n");
        exit(1);
    }
#endif

    fp = fopen("scvheoscalcentropy_gsl.txt", "w");


    /* Print the entropy table in the same format as the EOS tables. */
    for (i=0; i<Mat->nT; i++) {
        for (j=0; j<Mat->nRho; j++) {
            //fprintf(fp, "%15.7E%15.7E%15.7E\n", Mat->rhoAxis[j], Mat->TAxis[i], s[j*Mat->nT+i]);
            rho = pow(10.0, Mat->dLogRhoAxis[j]);
            T = pow(10.0, Mat->dLogTAxis[i]);
        
            fprintf(fp, "%15.7E%15.7E%15.7E%15.7E\n", rho, T, pow(10.0, Mat->dLogSArray[i][j]), SolveIntegral(Mat, pow(10.0, Mat->dLogRhoAxis[j]), pow(10.0, Mat->dLogTAxis[i]), rho0, T0, s0));
        }
    }
    fclose(fp);

    /* Free memory. */
    scvheosFinalizeMaterial(Mat);

    return 0;
}

