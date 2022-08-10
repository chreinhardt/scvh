/*
 * Copyright (c) 2020 Christian Reinhardt
 *
 * This file proves an interface for the SCvH equation of state (Saumon et al. 1995) for Hydrogen
 * and Helium that was extended to lower temperatures by Vazan et al. (2013).
 *
 * Parts of the code are from Thomas Meier's ANEOS library (www.github.com...).
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 08.08.2022
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include "scvheos.h"

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
    int nSkip;
    /* GSL interpolator type */
    const gsl_interp2d_type *InterpType;
    int bInterpBilinear = 1;

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
            nSkip = 1;
            break;
        case SCVHEOS_HE:
            /*
             * Helium.
             */
            strcpy(inFile, "scvh_he_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            break;
        case SCVHEOS_HHE:
            /*
             * Hydrogen / Helium mixture in Solar abundance.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            break;
        case SCVHEOS_HHE_LOWRHOT:
            /*
             * Hydrogen / Helium mixture in Solar abundance limited to low rho and T.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs_lowrhot.txt");
            nRho = 80;
            nT = 49;
            nSkip = 2;
            break;
        default:
            /* Unknown material */
            scvheosFinalizeMaterial(Mat);
            assert(0);
    }

    /*
     * Allocate memory and read the EOS table.
     */
    if (scvheosReadTable(Mat, inFile, nRho, nT, nSkip) != SCVHEOS_SUCCESS) {
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

    /* Define limits of the table (or extrapolation). */
    Mat->LogRhoMin = -15.0;
    Mat->LogRhoMax = 3.0;
    Mat->LogTMin = 1.0;
    Mat->LogTMax = 1e8;

    /*
     * Initialize the GSL interpolator.
     */
    if (bInterpBilinear) {
        InterpType = gsl_interp2d_bilinear;
    } else {
        InterpType = gsl_interp2d_bicubic;
    }

    Mat->xAccP = gsl_interp_accel_alloc();
    Mat->yAccP = gsl_interp_accel_alloc();
    Mat->xAccU = gsl_interp_accel_alloc();
    Mat->yAccU = gsl_interp_accel_alloc();
    Mat->xAccS = gsl_interp_accel_alloc();
    Mat->yAccS = gsl_interp_accel_alloc();
    Mat->xAccCs = gsl_interp_accel_alloc();
    Mat->yAccCs = gsl_interp_accel_alloc();

    /* x corresponts to T and y corresponds to rho. */
    Mat->InterpLogP = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);
    Mat->InterpLogU = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);
    Mat->InterpLogS = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);
    Mat->InterpLogCs = gsl_interp2d_alloc(InterpType, Mat->nT, Mat->nRho);

    gsl_interp2d_init(Mat->InterpLogP, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogPArray, Mat->nT, Mat->nRho);
    gsl_interp2d_init(Mat->InterpLogU, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogUArray, Mat->nT, Mat->nRho);
    gsl_interp2d_init(Mat->InterpLogS, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogSArray, Mat->nT, Mat->nRho);

#if 0
    /* The sound speed can only be calculated after the interpolation in P and u is initialized. */
    if (scvheosGenerateSoundSpeedTable(Mat) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "scvheosInitMaterial: Could not generate table for sound speed.\n");
        exit(1);
    }
#endif
    gsl_interp2d_init(Mat->InterpLogCs, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogCArray, Mat->nT, Mat->nRho);

    return Mat;
}

/*
 * Free memory.
 */
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat) {
    if (Mat != NULL) {
        if (Mat->dLogRhoAxis != NULL) free(Mat->dLogRhoAxis);
        if (Mat->dLogTAxis != NULL) free(Mat->dLogTAxis);

        if (Mat->dLogUArray != NULL) free(Mat->dLogUArray);
        if (Mat->dLogPArray != NULL) free(Mat->dLogPArray); 
        if (Mat->dLogSArray != NULL) free(Mat->dLogSArray);
        if (Mat->dLogCArray != NULL) free(Mat->dLogCArray);

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
int scvheosReadTable(SCVHEOSMAT *Mat, char *chInFile,  int nRho, int nT, int nSkip) {
    FILE *fp;
    char *chLine;
    size_t nCharMax = 256;
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

    /* The GSL interpolation functions require the tables are 1D arrays. */
    Mat->dLogUArray = (double *) calloc(nT*nRho, sizeof(double));
    Mat->dLogPArray = (double *) calloc(nT*nRho, sizeof(double));
    Mat->dLogSArray = (double *) calloc(nT*nRho, sizeof(double));

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

            /* We store the data as A[T][rho]. */
            iRet = sscanf(chLine, "%lf %lf %lf %lf %lf", &Mat->dLogTAxis[i], &Mat->dLogRhoAxis[j],
                                                         &Mat->dLogPArray[j*nT+i],
                                                         &Mat->dLogUArray[j*nT+i],
                                                         &Mat->dLogSArray[j*nT+i]);

            /* Check if the number of matches is correct. */
            if (iRet != 5) {
                return SCVHEOS_FAIL;
            }

            if ((pow(Mat->dLogBase, Mat->dLogRhoAxis[j]) < 0.0) ||
                (pow(Mat->dLogBase, Mat->dLogTAxis[i]) < 0.0) ||
                (pow(Mat->dLogBase, Mat->dLogUArray[j*nT+i]) < 0.0) ||
                (pow(Mat->dLogBase, Mat->dLogSArray[j*nT+i]) < 0.0)) {
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
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logP;

    if ((logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogPofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        exit(1);
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogP, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogPArray, logT, logrho,
			    Mat->xAccP, Mat->yAccP, &logP) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogPofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
	    exit(1);
    }
 
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
    P = pow(Mat->dLogBase, scvheosLogPofLogRhoLogT(Mat, logrho, logT));

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
double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logu;

    if ((logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogUofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        exit(1);
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogU, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogUArray, logT, logrho,
			    Mat->xAccU, Mat->yAccU, &logu) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogUofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
	    exit(1);
    }
    
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
    
    u = pow(Mat->dLogBase, scvheosLogUofLogRhoLogT(Mat, logrho, logT));

    return u;
}

/*
 * Calculate logs(logrho, logT).
 */
double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logs;

    if ((logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogSofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        exit(1);
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogS, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogSArray, logT, logrho,
			    Mat->xAccS, Mat->yAccS, &logs) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogSofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
	    exit(1);
    }
    
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
    s = pow(Mat->dLogBase, scvheosLogSofLogRhoLogT(Mat, logrho, logT));

    return s;
}

/*
 * Calculate the derivative dlogP/dlogrho(logrho, logT).
 */
double scvheosdLogPdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    /* Finite difference. */
    double h = 1e-5*logrho;
    double dLogPdLogRho;

    if (!scvheosCheckBoundsLogRhoLogT(Mat, logrho, logT)) {
	    fprintf(stderr, "scvheosdPdRhoofRhoT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        exit(1);
    }

    /* If (rho, T) is inside of the EOS table use GSL. */
    if (scvheosCheckTableBoundsLogRhoLogT(Mat, logrho, logT)) {
        if (gsl_interp2d_eval_deriv_y_e(Mat->InterpLogP, Mat->dLogTAxis, Mat->dLogRhoAxis,
            Mat->dLogPArray, logT, logrho, Mat->xAccP, Mat->yAccP, &dLogPdLogRho) == GSL_EDOM) {
            fprintf(stderr, "scvheosdLogPdLogRhoofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
            exit(1);
        }
        return dLogPdLogRho;
    }
    
    if ((logrho-h > Mat->LogRhoMin) && (logrho+h < Mat->LogRhoMax)) {
        /* Central difference. */
        dLogPdLogRho = (scvheosLogPofLogRhoLogT(Mat, logrho+h, logT)-scvheosLogPofLogRhoLogT(Mat, logrho-h, logT))/(2.0*h);
    } else if (logrho-h > Mat->LogRhoMin) {
        /* Backward finite difference. */
        dLogPdLogRho = (scvheosLogPofLogRhoLogT(Mat, logrho, logT) - scvheosLogPofLogRhoLogT(Mat, logrho-h, logT))/h;
    } else if (logrho+h < Mat->LogRhoMax) {
        /* Forward finite difference. */
        dLogPdLogRho = (scvheosLogPofLogRhoLogT(Mat, logrho+h, logT) - scvheosLogPofLogRhoLogT(Mat, logrho, logT))/h;
    } else {
        /* Both points are problematic so h is reduced. */
        h *= 1e-4;
        dLogPdLogRho = (scvheosLogPofLogRhoLogT(Mat, logrho+h, logT) - scvheosLogPofLogRhoLogT(Mat, logrho-h, logT))/(2.0*h);
    }

    return dLogPdLogRho;
}

/*
 * Calculate the derivative dPdRho(rho, T).
 */
double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    /* Finite difference. */
    double h = 1e-5*rho;
    double dPdRho;

    dPdRho = (scvheosPofRhoT(Mat, rho+h, T) - scvheosPofRhoT(Mat, rho-h, T))/(2.0*h);
    return dPdRho;
}

/*
 * Calculate the derivative dPdT(rho, T).
 */
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    /* Finite difference. */
    double h = 1e-5*T;
    double dPdT;

    dPdT = (scvheosPofRhoT(Mat, rho, T+h) - scvheosPofRhoT(Mat, rho, T-h))/(2.0*h);
    return dPdT;
}

/*
 * Calculate the derivative dUdRho(rho, T).
 */
double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    /* Finite difference. */
    double h = 1e-5*rho;
    double dUdRho;

    dUdRho = (scvheosUofRhoT(Mat, rho+h, T)-scvheosUofRhoT(Mat, rho-h, T))/(2.0*h);
    return dUdRho;
}

/*
 * Calculate the derivative dUdT(rho, T).
 */
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    /* Finite difference. */
    double h = 1e-5*T;
    double dUdT;

    dUdT = (scvheosUofRhoT(Mat, rho, T+h)-scvheosUofRhoT(Mat, rho, T-h))/(2.0*h);
    return dUdT;
}

/*
 * Calculate the derivative dSdRho(rho, T).
 */
double scvheosdSdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    /* Finite difference. */
    double h = 1e-5*rho;
    double dSdRho;

    dSdRho = (scvheosSofRhoT(Mat, rho+h, T)-scvheosSofRhoT(Mat, rho-h, T))/(2.0*h);
    return dSdRho;
}
/*
 * Calculate the derivative dSdT(rho, T).
 */
double scvheosdSdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    /* Finite difference. */
    double h = 1e-5*T;
    double dSdT;

    dSdT = (scvheosSofRhoT(Mat, rho, T+h)-scvheosSofRhoT(Mat, rho, T-h))/(2.0*h);
    return dSdT;
}

/*
 * Check if (logrho, logT) is inside of the eos table.
 */
int scvheosCheckTableBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    if ((logrho >= Mat->dLogRhoAxis[0]) && (logrho <= Mat->dLogRhoAxis[Mat->nRho-1]) && (logT >= Mat->dLogTAxis[0]) &&
        (logT <= Mat->dLogTAxis[Mat->nT-1])) {
        return TRUE;
    }

    return FALSE;
}

/*
 * Check if LogRhoMin <= logrho <= LogRhoMax and LogTMin <= logT <= LogTMax.
 */
int scvheosCheckBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    if ((logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
        return FALSE;
    }

    return TRUE;
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
