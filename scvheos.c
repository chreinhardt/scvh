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

/*
 * Calculate logs(logrho, logT).
 */
double scvheosLogSofRhoT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logs;

    logs = interpolateValueBilinear(logrho, logT, Mat->nT, Mat->nRho, Mat->dLogRhoAxis, Mat->dLogTAxis, Mat->dLogSArray);
    
    return logs;
}

/*
 * Calculate the entropy s(rho, T).
 */
double scvheosSofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double s;

    logrho = log10(rho);
    logT = log10(T);

    /* Interpolate in the table. */
    s = pow(Mat->dLogBase, scvheosLogSofRhoT(Mat, logrho, logT));

    return s;
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
