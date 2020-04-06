/*
 * Copyright (c) 2020 Christian Reinhardt
 *
 * This file proves an interface for the SCvH equation of state (Saumon et al. 1995) for Hydrogen
 * and Helium.
 *
 * Parts of the code are from Thomas Meier's ANEOS library (www.github.com...).
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 05.04.2020
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "scvheos.h"
#include "interpBilinear.h"

/*
 * Initalize a material.
 *
 * Initialize the material data, read the EOS tables and convert them to code units if desired.
 */
SCVHEOSMAT *scvheosInitMaterial(int iMat, double dKpcUnit, double dMsolUnit) {
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    const double NA = 6.022e23;          /* Avogadro's number */
    SCVHEOSMAT *Mat;
    int nRho;
    int nT;
    char inFile[256];
    int i, j;

    /* 
     * Allocate memory and initialize the data.
     *
     * Memory for the EOS tables is allocated in scvheosReadTable. The number of data points nRho
     * and nT are set in the function if the table was read correctly.
     */
    Mat = (SCVHEOSMAT *) calloc(1, sizeof(SCVHEOSMAT));
    assert(Mat != NULL);

    Mat->iMat = iMat;
    Mat->dKpcUnit = dKpcUnit;
    Mat->dMsolUnit = dMsolUnit;
    Mat->nRho = 0;
    Mat->nT = 0;
    Mat->dLogBase = 10.0;

    /*
     * Load the EOS table.
     */
    switch(iMat)
    {
        case SCVHEOS_H:
            /*
             * Hydrogen.
             */
            strcpy(inFile, "scvh_h_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            break;
        case SCVHEOS_HE:
            /*
             * Helium.
             */
            strcpy(inFile, "scvh_he_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            break;
        case SCVHEOS_HHE:
            /*
             * Hydrogen / Helium mixture in Solar abundance.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            break;
        default:
            /* Unknown material */
            scvheosFinalizeMaterial(Mat);
            assert(0);
    }

    /*
     * Allocate memory and read the EOS table.
     */
    if (scvheosReadTable(Mat, inFile, nRho, nT) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "SCVH EOS: Could not open EOS table %s.\n", inFile);
        scvheosFinalizeMaterial(Mat);
        exit(1);
    }

    /*
     * Convert from cgs to code units.
     */
    if ((dKpcUnit > 0.0) && (dMsolUnit > 0.0))
    {
        Mat->dGasConst = Mat->dKpcUnit*KPCCM*KBOLTZ
            /MHYDR/GCGS/Mat->dMsolUnit/MSOLG;
        /* code energy per unit mass --> erg per g */
        Mat->dErgPerGmUnit = GCGS*Mat->dMsolUnit*MSOLG/(Mat->dKpcUnit*KPCCM);
        /* code density --> g per cc */
        Mat->dGmPerCcUnit = (Mat->dMsolUnit*MSOLG)/pow(Mat->dKpcUnit*KPCCM,3.0);
        /* code time --> seconds */
        Mat->dSecUnit = sqrt(1/(Mat->dGmPerCcUnit*GCGS));
    }
    return Mat;
}

/*
 * Free memory.
 */
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat) {
    int i;

    if (Mat != NULL) {
        if (Mat->dLogRhoAxis != NULL) free(Mat->dLogRhoAxis);
        if (Mat->dLogTAxis != NULL) free(Mat->dLogTAxis);

        if (Mat->dLogUArray != NULL) {
            for (i=0; i<Mat->nT; i++) {
                if (Mat->dLogUArray[i] != NULL) free(Mat->dLogUArray[i]);
            }
            free(Mat->dLogUArray);
        }

        if (Mat->dLogPArray != NULL) {
            for (i=0; i<Mat->nT; i++) {
                if (Mat->dLogPArray[i] != NULL) free(Mat->dLogPArray[i]);
            }
            free(Mat->dLogPArray);
        }

        if (Mat->dLogSArray != NULL) {
            for (i=0; i<Mat->nT; i++) {
                if (Mat->dLogSArray[i] != NULL) free(Mat->dLogSArray[i]);
            }
            free(Mat->dLogSArray);
        }

        if (Mat->dLogCArray != NULL) {
            for (i=0; i<Mat->nT; i++) {
                if (Mat->dLogCArray[i] != NULL) free(Mat->dLogCArray[i]);
            }
            free(Mat->dLogCArray);
        }

        free(Mat);
    }

    return SCVHEOS_SUCCESS;
}

/*
 * Read an EOS table from a file.
 *
 * The function also allocates memory of the tables and sets nRho and nT.
 *
 * SCVHEOS: 
 *
 * H:    nRho = 201, nT = 100, logarithmic spacing (base 10)
 * He:   nRho = 201, nT = 100, logarithmic spacing (base 10)
 * H/He: nRho = 201, nT = 100, logarithmic spacing (base 10)
 */
int scvheosReadTable(SCVHEOSMAT *Mat, char *chInFile,  int nRho, int nT) {
    FILE *fp;
    char *chLine;
    size_t nCharMax = 256;
    int nSkip = 1;
    int iRet;
    int i, j;

    if (Mat == NULL)  {
        return SCVHEOS_FAIL;
    }

    if ((Mat->dLogRhoAxis != NULL) || (Mat->dLogTAxis != NULL) || (Mat->dLogUArray != NULL) ||
            (Mat->dLogPArray != NULL) || (Mat->dLogSArray != NULL)) {
        return SCVHEOS_FAIL;
    }

    Mat->dLogRhoAxis = (double *) calloc(nRho, sizeof(double));
    Mat->dLogTAxis = (double *) calloc(nT, sizeof(double));

    /* Thomas' code requires the tables to be of the form u[T][rho] and P[T][rho]. */
    Mat->dLogUArray = (double **) calloc(nT, sizeof(double*));
    Mat->dLogPArray = (double **) calloc(nT, sizeof(double*));
    Mat->dLogSArray = (double **) calloc(nT, sizeof(double*));

    for (i=0; i<nT; i++) {
        Mat->dLogUArray[i] = (double *) calloc(nRho, sizeof(double));
        Mat->dLogPArray[i] = (double *) calloc(nRho, sizeof(double));
        Mat->dLogSArray[i] = (double *) calloc(nRho, sizeof(double));
    }

    if ((Mat->dLogRhoAxis == NULL) || (Mat->dLogTAxis == NULL) || (Mat->dLogUArray == NULL) ||
            (Mat->dLogPArray == NULL) || (Mat->dLogSArray == NULL)) {
        return SCVHEOS_FAIL;
    }

    chLine = (char *) calloc(nCharMax, sizeof(char));

    /* Open the file and skip the first few lines. */
    fp = fopen(chInFile, "r");
    if (fp == NULL) {
        return SCVHEOS_FAIL;
    }

    for (i=0; i<nSkip; i++) {
        if (getline(&chLine, &nCharMax, fp) == -1) {
            return SCVHEOS_FAIL;
        }
    }

    for (i=0; i<nT; i++) {
        for (j=0; j<nRho; j++) {
            if (getline(&chLine, &nCharMax, fp) == -1) {
                return SCVHEOS_FAIL;
            }

            iRet = sscanf(chLine, "%lf %lf %lf %lf %lf", &Mat->dLogTAxis[i], &Mat->dLogRhoAxis[j], &Mat->dLogPArray[i][j], &Mat->dLogUArray[i][j], &Mat->dLogSArray[i][j]);

            /* Check if the number of matches is correct. */
            if (iRet != 5) {
                return SCVHEOS_FAIL;
            }

            if ((pow(Mat->dLogBase, Mat->dLogRhoAxis[j]) < 0.0) ||
                (pow(Mat->dLogBase, Mat->dLogTAxis[i]) < 0.0) ||
                (pow(Mat->dLogBase, Mat->dLogUArray[i][j]) < 0.0) ||
                (pow(Mat->dLogBase, Mat->dLogSArray[i][j]) < 0.0)) {
                return SCVHEOS_FAIL;
            }
        }
    }

    fclose(fp);
    free(chLine);

    Mat->nRho = nRho;
    Mat->nT = nT;

    return SCVHEOS_SUCCESS;
}

#if 0
/*
 * Generate an array that contains the sound speed at each EOS table data point.
 * 
 * The sound speed is calculated from
 *
 * cs^2 = dP/drho + T/(rho^2*C_v)*(dP/dT)^2
 *
 * where
 *
 * c_v = du/dT
 *
 * is the specific heat capacity at constant volume. The derivatives are obtained
 * numerically.
 */
int scvheosGenerateSoundSpeedTable(SCVHEOSMAT *Mat) {
    double P;
    double dPdrho;
    double dPdT;
    double cv;
    double cs2;
    int i, j;

    if (Mat == NULL)
        return SCVHEOS_FAIL;
    
    if ((Mat->dLogRhoAxis == NULL) || (Mat->dLogTAxis == NULL))
        return SCVHEOS_FAIL;

    if (Mat->dLogCArray != NULL)
        return SCVHEOS_FAIL;

    Mat->dLogArray = (double **) calloc(Mat->nT, sizeof(double*));

    if (Mat->dLogArray == NULL)
        return SCVHEOS_FAIL;

    /* nRho and nT are set in scvheosReadTable(). */
    for (i=0; i<Mat->nT; i++) {
        Mat->dLogCArray[i] = (double *) calloc(Mat->nRho, sizeof(double));
        if (Mat->dLogArray[i] == NULL)
            return SCVHEOS_FAIL;

    }

    for (i=0; i<Mat->nT; i++) {
        for (j=0; j<Mat->nRho; j++) {
            P = scvheosPofRhoT(Mat, Mat->dLogRhoAxis[j], Mat->dLogTAxis[j]);
            dPdrho = scvheosdPdrhoofRhoT(Mat, Mat->dLogRhoAxis[j], Mat->dLogTAxis[j]);
            dPdT = scvheosdPdTofRhoT(Mat, Mat->dLogRhoAxis[j], Mat->dLogTAxis[j]);
            cv = scvheosdUdrhoofRhoT(Mat, Mat->dLogRhoAxis[j], Mat->dLogTAxis[j]);

            cs2 = dPdrho + Mat->dLogTAxis[i]/(Mat->rhoAxis[j]*Mat->rhoAxis[j]*cv)*dPdT;

            assert(cs2 > 0.0);
            Mat->cArray[i][j] = sqrt(cs2);
        }
    }
}
#endif
// Functions that have to be implemented or added
// We also need derivatives to calculate the sound speed (maybe do this once and make a table)?


/*
 * Calculate logP(logrho, logT).
 */
double scvheosLogPofRhoT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logP;

    logP = interpolateValueBilinear(logrho, logT, Mat->nT, Mat->nRho, Mat->dLogRhoAxis, Mat->dLogTAxis, Mat->dLogPArray);
    
    return logP;
}

/*
 * Calculate the pressure P(rho, T).
 */
double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double P;

    logrho = log10(rho);
    logT = log10(T);

    /* Interpolate in the table. */
    P = pow(Mat->dLogBase, scvheosLogPofRhoT(Mat, logrho, logT));

    if (P < 0.0) {
        fprintf(stderr, "scvheosPofRhoT: Negative pressure for rho=%15.7E T=%15.7E (P=%15.7E)\n", rho, T, P);
    }

    return P;
}

#if 0
/*
 * Calculate the pressure P(rho, u).
 */
double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double T;
    double logP;
    double P;

    T = scvheosTofRhoU(Mat, rho, u);
    logP = interpolateValueBilinear(rho, T, Mat->nT, Mat->nRho, Mat->rhoAxis, Mat->TAxis, Mat->pArray);

    P = pow(Mat->dLogBase, logP);

    if (P < 0.0) {
        fprintf(stderr, "scvheosPofRhoU: Negative pressure for rho=%15.7E T=%15.7E (P=%15.7E)\n", rho, T, P);
    }

    return P;
}
#endif

/*
 * Calculate logU(logrho, logT).
 */
double scvheosLogUofRhoT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logu;

    logu = interpolateValueBilinear(logrho, logT, Mat->nT, Mat->nRho, Mat->dLogRhoAxis,
                                    Mat->dLogTAxis, Mat->dLogUArray);

    return logu;
}

/*
 * Calculate the internal energy u(rho, T).
 */
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double u;

    logrho = log10(rho);
    logT = log10(T);
    
    u = pow(Mat->dLogBase, scvheosLogUofRhoT(Mat, logrho, logT));

    if (u < 0.0) {
        fprintf(stderr, "scvheosUofRhoT: Negative internal energy for rho=%15.7E T=%15.7E (u=%15.7E)\n", rho, T, u);
    }

    return u;
}

#if 0
/*
 * Calculate the sound speed cs(rho, T).
 */
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double cs;

    cs = interpolateValueBilinear(rho, T, Mat->nT, Mat->nRho, Mat->rhoAxis, Mat->TAxis, Mat->cArray);

    if (cs > 0.0) {
    } else {
        fprintf(stderr, "scvheosCsofRhoT: Failed for rho=%15.7E T=%15.7E (cs=%15.7E)\n", rho, T, cs);
    }
    return cs;
}

/*
 * Calculate the sound speed cs(rho, u).
 */
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double T;
    double cs;

    T = scvheosTofRhoU(Mat, rho, u);
    cs = interpolateValueBilinear(rho, T, Mat->nT, Mat->nRho, Mat->rhoAxis, Mat->TAxis, Mat->cArray);

    if (cs > 0.0) {
    } else {
        fprintf(stderr, "scvheosCsofRhoT: Failed for rho=%15.7E T=%15.7E (cs=%15.7E)\n", rho, T, cs);
    }
    return cs;
}
#endif

#if 0
/*
 * Calculate T(rho, u).
 */
double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double logT;
    double T;

    logT = backwardInterpolateTemperatureBilinear(rho, u, Mat->nT, Mat->nRho, Mat->rhoAxis, Mat->TAxis, Mat->uArray);

    if (T < 0.0) {
        fprintf(stderr, "scvheosTofRhoU: Failed for rho=%15.7E u=%15.7E (T=%15.7E)\n", rho, u, T);
    }
    return T;
}

/*
 * Calculate dPdrho(rho, T) for constant T.
 */
double scvheosdPdrhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {

}

/*
 * Calculate dPdT(rho, T) for constant rho.
 */
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {


}

/*
 * Calculate dudrho(rho, T) for constant T.
 */
double scvheosdUdrhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {

}
#endif
#if 0
/*
 * Below are the interpolation functions from Thomas' ANEOS code.
 *
 * interpolateValueBilinear: interpolate a quantity z(rho, T) in rho and T
 *
 * backwardInterpolateTemperatureBilinear: determine T(rho, z)
 *
 * backwardInterpolateDensityBilinear: determine rho(T, z)
 *
 * Currently all functions do bilinear interpolation, for a more smooth EOS as SCVHEOS a higher order
 * interpolator could work too.
 */

/*
 * Internal function to backward interpolate the temperature using bilinear interpolation
 */
double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // searching the rho interval containing the rho value
    int i=0;
    for (int k=0; k<nRho; k++)
    {
        if (rhoAxis[k]>rho)
        {
            i = k-1;
            break;
        }
    }

    // selecting the rectangles that may contain the correct value
    int *indices = (int *)malloc((nT-1) * sizeof(int));
    int qq =0;
    for (int k=0; k<nT-1; k++)
    {
        // calculate minimum of the four points
        double minimum = fmin(fmin(fmin(zArray[k][i],zArray[k+1][i]),zArray[k][i+1]),zArray[k+1][i+1]);
        // calculate maximum of the four points
        double maximum = fmax(fmax(fmax(zArray[k][i],zArray[k+1][i]),zArray[k][i+1]),zArray[k+1][i+1]);

        if (z <= maximum && z >= minimum)
        {
            qq=qq+1;
            indices[k] = 1;
        } else {
            indices[k] = 0;
        }
    }

    // calculating the inverted bilinear interpolation for each of the selected rectangles
    // until the value is found
    double T = -1e50;

    for (int j=0; j<(nT-1); j++)
    {
        if (indices[j]==0)
        {
            continue;
        }
        double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);

        double f00=zArray[j][i];
        double f01=zArray[j+1][i];
        double f10=zArray[j][i+1];
        double f11=zArray[j+1][i+1];

        double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));

        if (y>= 0 && y<=1)
        {
            T = (TAxis[j+1]-TAxis[j])*y+TAxis[j];
            break;
        }

    }

    free(indices);
    return T;
}

/*
 * Internal function to backward interpolate the density using bilinear interpolation
 */
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // searching the T interval containing the T value
    int i=0;
    for (int k=0; k<nT; k++)
    {
        if (TAxis[k]>T)
        {
            i = k-1;
            break;
        }
    }

    // selecting the rectangles that may contain the correct value
    int *indices = (int *)malloc((nRho-1) * sizeof(int));
    int qq =0;
    for (int k=0; k<nRho-1; k++)
    {
        // calculate minimum of the four points
        double minimum = fmin(fmin(fmin(zArray[i][k],zArray[i+1][k]),zArray[i][k+1]),zArray[i+1][k+1]);
        // calculate maximum of the four points
        double maximum = fmax(fmax(fmax(zArray[i][k],zArray[i+1][k]),zArray[i][k+1]),zArray[i+1][k+1]);

        if (z <= maximum && z >= minimum)
        {
            qq=qq+1;;
            indices[k] = 1;
        } else {
            indices[k] = 0;
        }
    }

    // calculating the inverted bilinear interpolation for each of the selected rectangles
    // until the value is found
    double rho = -1e50;

    for (int j=0; j<(nRho-1); j++)
    {
        if (indices[j]==0)
        {
            continue;
        }
        double y=(T-TAxis[i])/(TAxis[i+1]-TAxis[i]);

        double f00=zArray[i][j];
        double f01=zArray[i+1][j];
        double f10=zArray[i][j+1];
        double f11=zArray[i+1][j+1];

        double x = -(z - f01*y + f00*(y - 1))/(y*(f01 - f11) - (f00 - f10)*(y - 1));

        if (x>= 0 && x<=1)
        {
            rho = (rhoAxis[j+1]-rhoAxis[j])*x+rhoAxis[j];
            break;
        }
    }

    free(indices);
    return rho;
}

/*
 * Internal function to interpolate a value using bilinear interpolation
 */
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // searching for grid rectangle containing the point
    // could be calculated if grid is guarantied to be logarithmic
    // as we do not assume that, we search for the grid rectangle
    int i=0;
    for (int k=0; k<nRho; k++)
    {
        if (rhoAxis[k]>rho)
        {
            i = k-1;
            break;
        }
    }
    int j=0;
    for (int k=0; k<nT; k++)
    {
        if (TAxis[k]>T)
        {
            j = k-1;
            break;
        }
    }

    // point not in grid
    if (i < 0 || j < 0)
    {
        fprintf(stderr, "interpolateValueBilinear: rho=%15.7E T=%15.7E is not in the grid (i=%i j=%i)\n", rho, T, i, j);
        return -1e50;
    }

    // CR: Check if rho < rho_max and T<T_max

    // scaling grid rectangle to unit square
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);

    // calculating bilinear interpolation
    double f00=zArray[j][i];
    double f01=zArray[j+1][i];
    double f10=zArray[j][i+1];
    double f11=zArray[j+1][i+1];

    double z = -1e50;

    z = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
    return z;
}

/*
 * Internal function to interpolate a value using bilinear interpolation in logarithmic quantities
 */
double interpolateValueBilinearLog(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    double dLogRho[2];
    double dLogT[2];

    // searching for grid rectangle containing the point
    // could be calculated if grid is guarantied to be logarithmic
    // as we do not assume that, we search for the grid rectangle
    int i=0;
    for (int k=0; k<nRho; k++)
    {
        if (rhoAxis[k]>rho)
        {
            i = k-1;
            break;
        }
    }
    int j=0;
    for (int k=0; k<nT; k++)
    {
        if (TAxis[k]>T)
        {
            j = k-1;
            break;
        }
    }

    // point not in grid
    if (i < 0 || j < 0)
    {
        fprintf(stderr, "interpolateValueBilinear: rho=%15.7E T=%15.7E is not in the grid (i=%i j=%i)\n", rho, T, i, j);
        return -1e50;
    }

    // CR: Check if rho < rho_max and T<T_max
    
    dLogRho[0] = log10(rhoAxis[i]);
    dLogRho[1] = log10(rhoAxis[i+1]);
    dLogT[0] = log10(TAxis[j]);
    dLogT[1] = log10(TAxis[j+1]);

    // scaling grid rectangle to unit square
    double x=(log10(rho)-dLogRho[0])/(dLogRho[1]-dLogRho[0]);
    double y=(log10(T)-dLogT[0])/(dLogT[1]-dLogT[0]);

    // calculating bilinear interpolation
    double f00=log10(zArray[j][i]);
    double f01=log10(zArray[j+1][i]);
    double f10=log10(zArray[j][i+1]);
    double f11=log10(zArray[j+1][i+1]);

    double z = -1e50;

    z = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
    return pow(10.0, z);
    //return z;
}
#endif
