/*
 * Copyright (c) 2020 Christian Reinhardt
 *
 * This file proves an interface for the SCvH equation of state (Saumon et al. 1995) for Hydrogen
 * and Helium that was extended to lower temperatures by Vazan et al. (2013).
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 16.04.2023
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_roots.h>
#include "scvheos.h"

/*
 * Lookup failure handling.
 *
 * The log-space functions return SCVHEOS_FAIL_VALUE when a lookup fails and
 * SCVH_NOEXIT is defined, otherwise they exit. The linear-space wrappers must
 * test the log value BEFORE exponentiating, because pow(10, -1e50) == 0 is a
 * perfectly valid-looking number and the failure would be lost.
 *
 * scvheosFailed(x) is true for the sentinel and for NaN, in log or linear space.
 */
static int scvheosFailed(double x) {
    return !(x > -1e40);
}

#ifdef SCVH_NOEXIT
#define SCVHEOS_FAIL_OR_EXIT() return SCVHEOS_FAIL_VALUE
#else
#define SCVHEOS_FAIL_OR_EXIT() exit(1)
#endif

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

    /* Set the reference density for ballic (careful: this fails for clumps!!). */ 
    Mat->rho0 = 1e-3; 

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
            strcpy(Mat->MatString, "Material not implemented yet.");
            break;
        case SCVHEOS_HE:
            /*
             * Helium.
             */
            strcpy(inFile, "scvh_he_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            strcpy(Mat->MatString, "Material not implemented yet.");
            break;
        case SCVHEOS_HHE:
            /*
             * Hydrogen / Helium mixture in Solar abundance.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs.txt");
            nRho = 201;
            nT = 100;
            nSkip = 1;
            strcpy(Mat->MatString, "Material not implemented yet.");
            break;
        case SCVHEOS_HHE_LOWRHOT:
            /*
             * Hydrogen / Helium mixture in Solar abundance limited to low rho and T.
             */
            strcpy(inFile, "scvh_hhe_y0.275_dt_cgs_lowrhot.txt");
            nRho = 80;
            nT = 48;
            nSkip = 2;
            strcpy(Mat->MatString, "SCvH EOS H-He Y=0.275 (Saumon et al. 1995, Vazan et al. 2013).");
            break;
        case SCVHEOS_HHE_EXT_LOWRHOT:
            /*
             * Hydrogen / Helium mixture (X=0.722, Y=0.278) based on the extended EOS tables
             * limited to low rho and T.
             */
#if 0
            strcpy(inFile, "scvh_extended_dt_hydrogen_722_helium_278_lowrhot.txt");
            nRho = 201;
            nT = 31;
            nSkip = 2;
#endif
            strcpy(inFile, "scvh_extended_dt_hydrogen_722_helium_278.data");
            nRho = 460;
            nT = 76;
            nSkip = 1;

            strcpy(Mat->MatString, "SCvH EOS H-He X=0.722 Y=0.278 (Saumon et al. 1995, Vazan et al. 2013).");
            break;
        case SCVHEOS_HE_EXT:
            /*
             * Helium EOS table extended to lower pressures and temperatures .
             */
            strcpy(inFile, "helium_scvh_extended_dt.data");
            nRho = 460;
            nT = 76;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS He extended to lower pressures and temperatures (Saumon et al. 1995, Vazan et al. 2013).");
            break;
        case SCVHEOS_HHE_H2O_Z_1:
            /*
             * H-He and H2O mixture obtained from Simon.
             */
            strcpy(inFile, "scvh_h2o_mixture_table_Z0.100.data");
            nRho = 128;
            nT = 64;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS H-He and H2O mixture Z=0.1.");
            break;
        case SCVHEOS_HHE_H2O_Z_05:
            /*
             * H-He and H2O mixture obtained from Simon.
             */
            strcpy(inFile, "scvh_h2o_mixture_table_Z0.050.data");
            nRho = 128;
            nT = 64;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS H-He and H2O mixture Z=0.05.");
            break;
        case SCVHEOS_HHE_H2O_Z_012:
            /*
             * H-He and H2O mixture obtained from Simon.
             */
            strcpy(inFile, "scvh_h2o_mixture_table_Z0.012.data");
            nRho = 128;
            nT = 64;
            nSkip = 2;

            strcpy(Mat->MatString, "SCvH EOS H-He and H2O mixture Z=0.012.");
            break;
        default:
            /* Unknown material */
            scvheosFinalizeMaterial(Mat);
            assert(0);
    }

    /* Currently only SCVHEOS_HHE_EXT_LOWRHOT works. */
    //assert(iMat == SCVHEOS_HHE_EXT_LOWRHOT);

    /*
     * Allocate memory and read the EOS table.
     */
    if (scvheosReadTable(Mat, inFile, nRho, nT, nSkip) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "SCVH EOS: Could not open EOS table %s.\n", inFile);
        scvheosFinalizeMaterial(Mat);
        exit(1);
    }

    /* Define limits of the table (or extrapolation). */
#if 0
    Mat->LogRhoMin = -25.0;
    Mat->LogRhoMax = 3.0;
    Mat->LogTMin = -2;
    Mat->LogTMax = 6;
#endif
#if 0
    //Mat->LogRhoMin = -10.5;
    Mat->LogRhoMin = -18.0;
    Mat->LogRhoMax = -0.9;
    Mat->LogTMin = -2.0;
    Mat->LogTMax = 4.7;
#endif
#if 0
    Mat->LogRhoMin = Mat->dLogRhoAxis[0];
    Mat->LogRhoMax = Mat->dLogRhoAxis[Mat->nRho-1];
    Mat->LogTMin = Mat->dLogTAxis[0];
    //Mat->LogTMax = Mat->dLogTAxis[Mat->nT-1];
    Mat->LogTMax = 3.6;
#endif
#if 0
    Mat->LogRhoMin = -25.0;
    Mat->LogRhoMax = 1.0;
    Mat->LogTMin = -1;
    Mat->LogTMax = 4;
#endif
    Mat->LogRhoMin = Mat->dLogRhoAxis[0];
    Mat->LogRhoMax = Mat->dLogRhoAxis[Mat->nRho-1];
    Mat->LogTMin = Mat->dLogTAxis[0];
    Mat->LogTMax = Mat->dLogTAxis[Mat->nT-1];

    /*
     * Convert from cgs to code units.
     */
    if ((Mat->dKpcUnit > 0.0) && (Mat->dMsolUnit > 0.0)) {
        Mat->dGasConst = Mat->dKpcUnit*KPCCM*KBOLTZ
            /MHYDR/GCGS/Mat->dMsolUnit/MSOLG;
        /* code energy per unit mass --> erg per g */
        Mat->dErgPerGmUnit = GCGS*Mat->dMsolUnit*MSOLG/(Mat->dKpcUnit*KPCCM);
        /* code density --> g per cc */
        Mat->dGmPerCcUnit = (Mat->dMsolUnit*MSOLG)/pow(Mat->dKpcUnit*KPCCM,3.0);
        /* code time --> seconds */
        Mat->dSecUnit = sqrt(1/(Mat->dGmPerCcUnit*GCGS));

        for (int i=0; i<Mat->nRho; i++) {
            Mat->dLogRhoAxis[i] -= log10(Mat->dGmPerCcUnit);
        }
        
        for (int i=0; i<Mat->nT; i++) {
            for (int j=0; j<Mat->nRho; j++) {
                Mat->dLogPArray[j*Mat->nT+i] -= log10(Mat->dErgPerGmUnit*Mat->dGmPerCcUnit);
                Mat->dLogUArray[j*Mat->nT+i] -= log10(Mat->dErgPerGmUnit);
                Mat->dLogSArray[j*Mat->nT+i] -= log10(Mat->dErgPerGmUnit);
            }
        }
 
        /* Convert table limits and reference density. */
        Mat->LogRhoMin -= log10(Mat->dGmPerCcUnit);
        Mat->LogRhoMax -= log10(Mat->dGmPerCcUnit);
        Mat->rho0 /= Mat->dGmPerCcUnit;
    } else {
        /* Prevent problems if dKpcUnit or dMsolUnit are not set. */
        Mat->dGasConst = 0.0;
        Mat->dErgPerGmUnit = 1.0;
        Mat->dGmPerCcUnit = 1.0;
        Mat->dSecUnit = 1.0;
    }

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

    /* The sound speed can only be calculated after the interpolation in P and u is initialized. */
    if (scvheosGenerateSoundSpeedTable(Mat) != SCVHEOS_SUCCESS) {
        fprintf(stderr, "scvheosInitMaterial: Could not generate table for sound speed.\n");
        exit(1);
    }

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
 *
 * Where the table is defective (c_v <= 0 or a negative cs^2, which happens
 * e.g. in the H-He/H2O Z=0.05 table) cs^2 is replaced by P/rho at the node:
 * this has the right units whatever the table is in, and for a gas it is a
 * lower bound on cs^2 = gamma*P/rho, so the hydro code errs towards a short
 * timestep. A zero or -inf must never be stored, because bilinear
 * interpolation would spread it into the neighbouring cells. The number of
 * replaced nodes is reported on stderr.
 */
int scvheosGenerateSoundSpeedTable(SCVHEOSMAT *Mat) {
    double rho;
    double T;
    double P;
    double dPdrho;
    double dPdT;
    double cv;
    double cs2;
    int nReplaced = 0;
    const int nReport = 10;   /* print at most this many defective nodes */

    if (Mat == NULL)
        return SCVHEOS_FAIL;
    
    if ((Mat->dLogRhoAxis == NULL) || (Mat->dLogTAxis == NULL))
        return SCVHEOS_FAIL;

    if (Mat->dLogCArray != NULL)
        return SCVHEOS_FAIL;

    /* The GSL interpolation functions require the tables are 1D arrays. */
    Mat->dLogCArray = (double *) calloc(Mat->nT*Mat->nRho, sizeof(double));

    if (Mat->dLogCArray == NULL)
        return SCVHEOS_FAIL;

    for (int i=0; i<Mat->nT; i++) {
        for (int j=0; j<Mat->nRho; j++) {
            rho = pow(Mat->dLogBase, Mat->dLogRhoAxis[j]);
            T = pow(Mat->dLogBase, Mat->dLogTAxis[i]);
            P = pow(Mat->dLogBase, Mat->dLogPArray[j*Mat->nT+i]);

            dPdrho = scvheosdPdRhoofRhoT(Mat, rho, T);
            dPdT = scvheosdPdTofRhoT(Mat, rho, T);
            cv = scvheosdUdTofRhoT(Mat, rho, T);

            cs2 = dPdrho + T/(rho*rho*cv)*dPdT*dPdT;

            /* Defective node: non-finite (cv = 0), or not positive. */
            if (!isfinite(cs2) || !(cs2 > 0.0)) {
                if (nReplaced < nReport) {
                    fprintf(stderr, "scvheosGenerateSoundSpeedTable: iMat= %i iT= %i jRho= %i "
                            "logrho= %g logT= %g dPdRho= %g dPdT= %g cv= %g cs2= %g, "
                            "replaced by P/rho= %g\n", Mat->iMat, i, j,
                            Mat->dLogRhoAxis[j], Mat->dLogTAxis[i], dPdrho, dPdT, cv, cs2, P/rho);
                }
                nReplaced++;
                cs2 = P/rho;
            }

            /* P/rho > 0 at every node (the table stores log P), so this is finite. */
            assert(isfinite(cs2) && cs2 > 0.0);

            Mat->dLogCArray[j*Mat->nT+i] = 0.5*log10(cs2);
        }
    }

    if (nReplaced > 0) {
        fprintf(stderr, "scvheosGenerateSoundSpeedTable: iMat= %i: cs^2 replaced by P/rho at %i of %i nodes.\n",
                Mat->iMat, nReplaced, Mat->nT*Mat->nRho);
    }

    return SCVHEOS_SUCCESS;
}

// Functions that have to be implemented or added
// We also need derivatives to calculate the sound speed (maybe do this once and make a table)?


/*
 * Calculate logP(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logP;

    if (!isfinite(logrho) || !isfinite(logT) ||
        (logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogPofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogP, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogPArray, logT, logrho,
			    Mat->xAccP, Mat->yAccP, &logP) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogPofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
 
    return logP;
}

/*
 * Calculate logU(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logu;

    if (!isfinite(logrho) || !isfinite(logT) ||
        (logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogUofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogU, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogUArray, logT, logrho,
			    Mat->xAccU, Mat->yAccU, &logu) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogUofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
    
    return logu;
}

/*
 * Calculate logs(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logs;

    if (!isfinite(logrho) || !isfinite(logT) ||
        (logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogSofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogS, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogSArray, logT, logrho,
			    Mat->xAccS, Mat->yAccS, &logs) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogSofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
    
    return logs;
}

/*
 * Calculate logcs(logrho, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if (logrho, logT)
 * is not finite or outside of the table.
 */
double scvheosLogCsofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    double logcs;

    if (!isfinite(logrho) || !isfinite(logT) ||
        (logrho < Mat->LogRhoMin) || (logrho > Mat->LogRhoMax) || (logT < Mat->LogTMin) || (logT > Mat->LogTMax)) {
	    fprintf(stderr, "scvheosLogCsofLogRhoLogT: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    if (gsl_interp2d_eval_e_extrap(Mat->InterpLogCs, Mat->dLogTAxis, Mat->dLogRhoAxis, Mat->dLogCArray, logT, logrho,
			    Mat->xAccCs, Mat->yAccCs, &logcs) != GSL_SUCCESS) {
	    fprintf(stderr, "scvheosLogCsofLogRhoLogT: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", logrho, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }
 
    return logcs;
}


/*
 * Calculate the pressure P(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logP;

    logrho = log10(rho);
    logT = log10(T);

    /* Interpolate in the table. */
    logP = scvheosLogPofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logP)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logP);
}

/*
 * Calculate the internal energy u(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logu;

    logrho = log10(rho);
    logT = log10(T);

    logu = scvheosLogUofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logu)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logu);
}

/*
 * Calculate the entropy s(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosSofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logs;

    logrho = log10(rho);
    logT = log10(T);

    /* Interpolate in the table. */
    logs = scvheosLogSofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logs)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logs);
}

/*
 * Calculate the sound speed cs(rho, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the lookup fails.
 */
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double logcs;

    logrho = log10(rho);
    logT = log10(T);

    logcs = scvheosLogCsofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logcs)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logcs);
}

/*
 * Calculate the pressure P(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion or the lookup fails.
 */
double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double logrho;
    double logu;
    double logT;
    double logP;
    double P;

    logrho = log10(rho);
    logu = log10(u);

    /* Calculate logT. */
    logT = scvheosLogTofLogRhoLogU(Mat, logrho, logu);
    if (scvheosFailed(logT)) return SCVHEOS_FAIL_VALUE;

    /* Interpolate in the table. */
    logP = scvheosLogPofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logP)) return SCVHEOS_FAIL_VALUE;

    P = pow(Mat->dLogBase, logP);

    // CR
    assert(P > 0.0);
    return P;
}

/*
 * Calculate the sound speed cs(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion or the lookup fails.
 */
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double logrho;
    double logu;
    double logT;
    double logcs;

    logrho = log10(rho);
    logu = log10(u);

    /* Calculate logT. */
    logT = scvheosLogTofLogRhoLogU(Mat, logrho, logu);
    if (scvheosFailed(logT)) return SCVHEOS_FAIL_VALUE;

    logcs = scvheosLogCsofLogRhoLogT(Mat, logrho, logT);
    if (scvheosFailed(logcs)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logcs);
}

/*
 * Calculate the temperature logT(logrho, logu).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if the input is
 * not finite or the root cannot be bracketed.
 */
double scvheosLogTofLogRhoLogU(SCVHEOSMAT *Mat, double logrho, double logu) {
    /* GSL root finder */
    gsl_root_fsolver *Solver;
    const gsl_root_fsolver_type *SolverType;
    gsl_function F;
    struct LogUofLogRhoLogT_GSL_Params Params;
    const double err_abs = 0.0;
    const double err_rel = 1e-10;
    int status;
    int max_iter = 1000;
    double logT_min, logT_max;
    double logT = 0.0;

    /* log10 of a non-positive rho or u is NaN/-inf: nothing to invert. */
    if (!isfinite(logrho) || !isfinite(logu)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogU: invalid input logrho= %g logu= %g.\n", logrho, logu);
        SCVHEOS_FAIL_OR_EXIT();
    }

    /* Initialize the parameters. */
    Params.Mat = Mat;
    Params.logrho = logrho;
    Params.logu = logu;

    /* Initialize the function used for root finding. */
    F.function = &LogUofLogRhoLogT_GSL_rootfinder;
    F.params = &Params;

    /* Initialize the root finder. */
    SolverType = gsl_root_fsolver_brent;
    SolverType = gsl_root_fsolver_bisection;
    Solver = gsl_root_fsolver_alloc(SolverType);
    assert(Solver != NULL);

    /* Set minimum and maximum temperature. */ 
    logT_min = Mat->LogTMin;
    logT_max = Mat->LogTMax;

    /* Check if logu < logu(logrho, logT_min) or logu > logu(logrho, logT_max) and set a minimum or maximum value. */
#if 0
    if (logu < scvheosLogUofLogRhoLogT(Mat, logrho, logT_min)) return logT_min;
    if (logu > scvheosLogUofLogRhoLogT(Mat, logrho, logT_max)) return logT_max;
#endif
    /*
    /// CR: Careful, this assumes that u(T) is monotonic AND u(T_min) < u(T_max)
    if (logu < scvheosLogUofLogRhoLogT(Mat, logrho, logT_min)) assert(logu >= scvheosLogUofLogRhoLogT(Mat, logrho, logT_min));
    if (logu > scvheosLogUofLogRhoLogT(Mat, logrho, logT_max)) assert(logu <= scvheosLogUofLogRhoLogT(Mat, logrho, logT_max));
    */

    /* Make sure the root is bracketed.*/
    if (!(LogUofLogRhoLogT_GSL_rootfinder(logT_min, &Params)*LogUofLogRhoLogT_GSL_rootfinder(logT_max, &Params) <= 0.0)) {
        // CR: 04.10.2022
        fprintf(stderr, "scvheosLogTofLogRhoLogU:\n");

        fprintf(stderr, "logrho= %g ", logrho);
        fprintf(stderr, "logu= %g\n", logu);
        fprintf(stderr, "LogTMin= %g ", Mat->LogTMin);
        fprintf(stderr, "LogTMax= %g\n", Mat->LogTMax);
        fprintf(stderr, "loguMin= %g\n", scvheosLogUofLogRhoLogT(Mat, logrho, Mat->LogTMin));
        fprintf(stderr, "loguMax= %g\n", scvheosLogUofLogRhoLogT(Mat, logrho, Mat->LogTMax));
        fprintf(stderr, "loguMax= %g\n", LogUofLogRhoLogT_GSL_rootfinder(logT_min, &Params)*LogUofLogRhoLogT_GSL_rootfinder(logT_max, &Params));

        fprintf(stderr, "Could not bracket root.\n");
        gsl_root_fsolver_free(Solver);
        SCVHEOS_FAIL_OR_EXIT();
    }

    gsl_root_fsolver_set(Solver, &F, logT_min, logT_max);

    for (int i=0; i<max_iter; i++) {
        /* Do one iteration of the root solver. */
        status = gsl_root_fsolver_iterate(Solver);

        /* Estimate of the root. */
        logT = gsl_root_fsolver_root(Solver);

        /* Current interval that brackets the root. */
        logT_min = gsl_root_fsolver_x_lower(Solver);
        logT_max = gsl_root_fsolver_x_upper(Solver);

        /* Test for convergence. */
        status = gsl_root_test_interval(logT_min, logT_max, err_abs, err_rel);

#if 0
        if (status == GSL_SUCCESS)
            fprintf(stderr, "Converged: x= %g\n", x);
#endif
        if (status != GSL_CONTINUE) break;
    }

    //if (status != GSL_SUCCESS) logT = -1.0;
    // Add stricter tests
    assert(status == GSL_SUCCESS); 
    double logu_int = scvheosLogUofLogRhoLogT(Mat, logrho, logT);

    if (fabs(logu-logu_int)/logu < 1e-3) {
    } else {
        if (fabs(logu) > 1e-6) {
            fprintf(stderr, "scvheosLogTofLogRhoLogU: logrho=%15.7E logu=%15.7E logT=%15.7E logu_int=%15.7E\n", logrho, logu, logT, logu_int);
            assert(fabs(logu-logu_int)/logu < 1e-3);
        }
    }

    //assert(fabs(logu-logu_int)/logu < 1e-3);

    gsl_root_fsolver_free(Solver);

    return logT;
}

/*
 * Calculate the temperature logT(logrho, logs).
 *
 * Returns SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if the input is
 * not finite or the root cannot be bracketed.
 */
double scvheosLogTofLogRhoLogS(SCVHEOSMAT *Mat, double logrho, double logs) {
    /* GSL root finder */
    gsl_root_fsolver *Solver;
    const gsl_root_fsolver_type *SolverType;
    gsl_function F;
    struct LogSofLogRhoLogT_GSL_Params Params;
    const double err_abs = 0.0;
    const double err_rel = 1e-10;
    int status;
    int max_iter = 1000;
    double logT_min, logT_max;
    double logT = 0.0;

    if (!isfinite(logrho) || !isfinite(logs)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogS: invalid input logrho= %g logs= %g.\n", logrho, logs);
        SCVHEOS_FAIL_OR_EXIT();
    }

    /* Initialize the parameters. */
    Params.Mat = Mat;
    Params.logrho = logrho;
    Params.logs = logs;

    /* Initialize the function used for root finding. */
    F.function = &LogSofLogRhoLogT_GSL_rootfinder;
    F.params = &Params;

    /* Initialize the root finder. */
    SolverType = gsl_root_fsolver_brent;
    SolverType = gsl_root_fsolver_bisection;
    Solver = gsl_root_fsolver_alloc(SolverType);
    assert(Solver != NULL);

    /* Set minimum and maximum temperature. */ 
    logT_min = Mat->LogTMin;
    logT_max = Mat->LogTMax;

    /* Check if logs < logs(logrho, logT_min) or logs > logs(logrho, logT_max) and set a minimum or maximum value. */
    if (logs < scvheosLogSofLogRhoLogT(Mat, logrho, logT_min)) {
        gsl_root_fsolver_free(Solver);
        return logT_min;
    }

    if (logs > scvheosLogSofLogRhoLogT(Mat, logrho, logT_max)) {
        gsl_root_fsolver_free(Solver);
        return logT_max;
    }

    /* Make sure the root is bracketed.*/
    if (!(LogSofLogRhoLogT_GSL_rootfinder(logT_min, &Params)*LogSofLogRhoLogT_GSL_rootfinder(logT_max, &Params) <= 0.0)) {
        fprintf(stderr, "scvheosLogTofLogRhoLogS: Could not bracket root (logrho= %g logs= %g).\n", logrho, logs);
        gsl_root_fsolver_free(Solver);
        SCVHEOS_FAIL_OR_EXIT();
    }

    gsl_root_fsolver_set(Solver, &F, logT_min, logT_max);

    for (int i=0; i<max_iter; i++) {
        /* Do one iteration of the root solver. */
        status = gsl_root_fsolver_iterate(Solver);

        /* Estimate of the root. */
        logT = gsl_root_fsolver_root(Solver);

        /* Current interval that brackets the root. */
        logT_min = gsl_root_fsolver_x_lower(Solver);
        logT_max = gsl_root_fsolver_x_upper(Solver);

        /* Test for convergence. */
        status = gsl_root_test_interval(logT_min, logT_max, err_abs, err_rel);

#if 0
        if (status == GSL_SUCCESS)
            fprintf(stderr, "Converged: x= %g\n", x);
#endif
        if (status != GSL_CONTINUE) break;
    }

    //if (status != GSL_SUCCESS) logT = -1.0;
    assert(status == GSL_SUCCESS); 
    double logs_int = scvheosLogSofLogRhoLogT(Mat, logrho, logT);
    assert(fabs(logs-logs_int)/logs < 1e-3);

    gsl_root_fsolver_free(Solver);

    return logT;
}

/*
 * Calculate the density logrho(logP, logT).
 *
 * Returns SCVHEOS_FAIL_VALUE if the root finder does not converge, and
 * SCVHEOS_FAIL_VALUE (or exits without SCVH_NOEXIT) if the input is not
 * finite or the root cannot be bracketed.
 */
double scvheosLogRhoofLogPLogT(SCVHEOSMAT *Mat, double logP, double logT) {
    /* GSL root finder */
    gsl_root_fsolver *Solver;
    const gsl_root_fsolver_type *SolverType;
    gsl_function F;
    struct LogPofLogRhoLogT_GSL_Params Params;
    const double err_abs = 0.0;
    const double err_rel = 1e-10;
    int status;
    int max_iter = 1000;
    double logrho_min, logrho_max;
    double logrho = 0.0;

    if (!isfinite(logP) || !isfinite(logT)) {
        fprintf(stderr, "scvheosLogRhoofLogPLogT: invalid input logP= %g logT= %g.\n", logP, logT);
        SCVHEOS_FAIL_OR_EXIT();
    }

    /* Initialize the parameters. */
    Params.Mat = Mat;
    Params.logT = logT;
    Params.logP = logP;

    /* Initialize the function used for root finding. */
    F.function = &LogPofLogRhoLogT_GSL_rootfinder;
    F.params = &Params;

    /* Initialize the root finder. */
    SolverType = gsl_root_fsolver_brent;
    SolverType = gsl_root_fsolver_bisection;
    Solver = gsl_root_fsolver_alloc(SolverType);
    assert(Solver != NULL);

    /* Set minimum and maximum density. */ 
    logrho_min = Mat->LogRhoMin;
    logrho_max = Mat->LogRhoMax;

    /* Make sure the root is bracketed.*/
    if (!(LogPofLogRhoLogT_GSL_rootfinder(logrho_min, &Params)*LogPofLogRhoLogT_GSL_rootfinder(logrho_max, &Params) <= 0.0)) {
        fprintf(stderr, "scvheosLogRhoofLogPLogT: Could not bracket root (logP= %g logT= %g).\n", logP, logT);
        gsl_root_fsolver_free(Solver);
        SCVHEOS_FAIL_OR_EXIT();
    }

    gsl_root_fsolver_set(Solver, &F, logrho_min, logrho_max);

    for (int i=0; i<max_iter; i++) {
        /* Do one iteration of the root solver. */
        status = gsl_root_fsolver_iterate(Solver);

        /* Estimate of the root. */
        logrho = gsl_root_fsolver_root(Solver);

        /* Current interval that brackets the root. */
        logrho_min = gsl_root_fsolver_x_lower(Solver);
        logrho_max = gsl_root_fsolver_x_upper(Solver);

        /* Test for convergence. */
        status = gsl_root_test_interval(logrho_min, logrho_max, err_abs, err_rel);

#if 0
        if (status == GSL_SUCCESS)
            fprintf(stderr, "Converged: x= %g\n", x);
#endif
        if (status != GSL_CONTINUE) break;
    }

    /*
     * No convergence: return the failure sentinel. A value such as -1.0 would
     * become a valid-looking density (0.1) after exponentiation.
     */
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "scvheosLogRhoofLogPLogT: root finder did not converge (logP= %g logT= %g).\n", logP, logT);
        logrho = SCVHEOS_FAIL_VALUE;
    }

    gsl_root_fsolver_free(Solver);

    return logrho;
}

/*
 * Calculate the temperature T(rho, u).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion fails.
 */
double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double logrho;
    double logu;
    double logT;

    logrho = log10(rho);
    logu = log10(u);

    logT = scvheosLogTofLogRhoLogU(Mat, logrho, logu);
    if (scvheosFailed(logT)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logT);
}

/*
 * Calculate the temperature T(rho, s).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion fails.
 */
double scvheosTofRhoS(SCVHEOSMAT *Mat, double rho, double s) {
    double logrho;
    double logs;
    double logT;

    logrho = log10(rho);
    logs = log10(s);

    logT = scvheosLogTofLogRhoLogS(Mat, logrho, logs);
    if (scvheosFailed(logT)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logT);
}

/*
 * Calculate the density rho(P, T).
 *
 * Returns SCVHEOS_FAIL_VALUE if the inversion fails (callers such as
 * EOSWoolfsonCoeff test for rho <= 0).
 */
double scvheosRhoofPT(SCVHEOSMAT *Mat, double P, double T) {
    double logP;
    double logT;
    double logrho;

    logP = log10(P);
    logT = log10(T);

    logrho = scvheosLogRhoofLogPLogT(Mat, logP, logT);
    if (scvheosFailed(logrho)) return SCVHEOS_FAIL_VALUE;

    return pow(Mat->dLogBase, logrho);
}

/*
 * Partial derivatives of the tabulated log quantities with respect to logrho
 * (bWrtRho = TRUE) or logT (bWrtRho = FALSE), from the GSL interpolator.
 *
 * The derivative is only defined inside the table. The extrapolation limits
 * LogRhoMin/Max and LogTMin/Max coincide with the table bounds, so the
 * previous finite-difference fallback for the extrapolated region was dead
 * code (and used a step proportional to logrho/logT, which vanishes at
 * logrho = 0 or logT = 0). If extrapolation is ever re-enabled a fallback
 * with an absolute step in log space (a fraction of the grid spacing) has
 * to be added here, since GSL has no extrapolating derivative.
 *
 * With bilinear interpolation the derivative is piecewise constant and
 * jumps at the cell boundaries; at a node GSL evaluates the cell to the
 * right/above (the last cell at the upper table edge).
 */
static double scvheosdLogdLog(SCVHEOSMAT *Mat, const gsl_interp2d *Interp, const double *Array,
                              gsl_interp_accel *xAcc, gsl_interp_accel *yAcc,
                              double logrho, double logT, int bWrtRho, const char *chFunc) {
    double dLogdLog;
    int status;

    if (!scvheosCheckTableBoundsLogRhoLogT(Mat, logrho, logT)) {
	    fprintf(stderr, "%s: logrho= %15.7E logT= %15.7E outside of the EOS table.\n", chFunc, logrho, logT);
        exit(1);
    }

    /* x corresponds to T and y to rho. */
    if (bWrtRho) {
        status = gsl_interp2d_eval_deriv_y_e(Interp, Mat->dLogTAxis, Mat->dLogRhoAxis, Array,
                                             logT, logrho, xAcc, yAcc, &dLogdLog);
    } else {
        status = gsl_interp2d_eval_deriv_x_e(Interp, Mat->dLogTAxis, Mat->dLogRhoAxis, Array,
                                             logT, logrho, xAcc, yAcc, &dLogdLog);
    }

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "%s: Interpolation failed (logrho= %15.7E logT= %15.7E).\n", chFunc, logrho, logT);
        exit(1);
    }

    return dLogdLog;
}

/*
 * Calculate the derivative dlogP/dlogrho(logrho, logT).
 */
double scvheosdLogPdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->InterpLogP, Mat->dLogPArray, Mat->xAccP, Mat->yAccP,
                           logrho, logT, TRUE, __func__);
}

/*
 * Calculate the derivative dlogP/dlogT(logrho, logT).
 */
double scvheosdLogPdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->InterpLogP, Mat->dLogPArray, Mat->xAccP, Mat->yAccP,
                           logrho, logT, FALSE, __func__);
}

/*
 * Calculate the derivative dlogu/dlogrho(logrho, logT).
 */
double scvheosdLogUdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->InterpLogU, Mat->dLogUArray, Mat->xAccU, Mat->yAccU,
                           logrho, logT, TRUE, __func__);
}

/*
 * Calculate the derivative dlogu/dlogT(logrho, logT).
 */
double scvheosdLogUdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->InterpLogU, Mat->dLogUArray, Mat->xAccU, Mat->yAccU,
                           logrho, logT, FALSE, __func__);
}

/*
 * Calculate the derivative dlogs/dlogrho(logrho, logT).
 */
double scvheosdLogSdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->InterpLogS, Mat->dLogSArray, Mat->xAccS, Mat->yAccS,
                           logrho, logT, TRUE, __func__);
}

/*
 * Calculate the derivative dlogs/dlogT(logrho, logT).
 */
double scvheosdLogSdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT) {
    return scvheosdLogdLog(Mat, Mat->InterpLogS, Mat->dLogSArray, Mat->xAccS, Mat->yAccS,
                           logrho, logT, FALSE, __func__);
}

/*
 * Calculate the derivative dP/drho(rho, T).
 */
double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double P;
    double dLogPdLogRho;
    
    logrho = log10(rho);
    logT = log10(T);

    // CR: Here it is probably a bit more efficient if we calculate P=10**log(P)
    P = scvheosPofRhoT(Mat, rho, T);
    dLogPdLogRho = scvheosdLogPdLogRhoofLogRhoLogT(Mat, logrho, logT);
    
    /* dP/dx = P/x*dln(P)/dln(rho). */
    return P/rho*dLogPdLogRho;
}

/*
 * Calculate the derivative dP/dT(rho, T).
 */
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double P;
    double dLogPdLogT;
    
    logrho = log10(rho);
    logT = log10(T);

    P = scvheosPofRhoT(Mat, rho, T);
    dLogPdLogT = scvheosdLogPdLogTofLogRhoLogT(Mat, logrho, logT);
    
    return P/T*dLogPdLogT;
}

/*
 * Calculate the derivative du/drho(rho, T).
 */
double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double u;
    double dLogUdLogRho;
    
    logrho = log10(rho);
    logT = log10(T);

    u = scvheosUofRhoT(Mat, rho, T);
    dLogUdLogRho = scvheosdLogUdLogRhoofLogRhoLogT(Mat, logrho, logT);
    
    return u/rho*dLogUdLogRho;
}

/*
 * Calculate the derivative du/dT(rho, T).
 */
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T) {
    double logrho;
    double logT;
    double u;
    double dLogUdLogT;
    
    logrho = log10(rho);
    logT = log10(T);

    u = scvheosUofRhoT(Mat, rho, T);
    dLogUdLogT = scvheosdLogUdLogTofLogRhoLogT(Mat, logrho, logT);
    
    return u/T*dLogUdLogT;
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
 * Calculate dPdrho(rho, u).
 */
double scvheosdPdRhoofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double h = 1e-6*rho;
    double dPdRho;
    double UMinRight;
    double UMinLeft;
    double TMin;

    /* Limit of the EOS table defines the minimum temperature. */
    TMin = pow(Mat->dLogBase, Mat->dLogTAxis[0]);

    /* Check if (rho-h, u) and (rho+h, u) are inside of the EOS table. */
    UMinLeft = scvheosUofRhoT(Mat, rho-h, TMin);
    UMinRight = scvheosUofRhoT(Mat, rho+h, TMin);

    if ((UMinLeft <= u) && (UMinRight <= u)) {
        /* Both left and right value are smaller than u. */
        dPdRho = (scvheosPofRhoU(Mat, rho+h, u) - scvheosPofRhoU(Mat, rho-h, u))/(2.0*h);
    } else if (UMinLeft <= u) {
        /* Backward finite difference. */
        dPdRho = (scvheosPofRhoU(Mat, rho, u) - scvheosPofRhoU(Mat, rho-h, u))/h;
    } else if (UMinRight <= u) {
        /* Forward finite difference. */
        dPdRho = (scvheosPofRhoU(Mat, rho+h, u) - scvheosPofRhoU(Mat, rho, u))/h;
    } else {
        /* Both points are problematic so h is reduced. */
        h *= 1e-4;
        dPdRho = (scvheosPofRhoU(Mat, rho+h, u) - scvheosPofRhoU(Mat, rho-h, u))/(2.0*h);
    }
    return dPdRho;
}

/*
 * Calculate dPdu(rho, u).
 */
double scvheosdPdUofRhoU(SCVHEOSMAT *Mat, double rho, double u) {
    double h = fabs(1e-5*u);
    double dPdU;

    dPdU = (scvheosPofRhoU(Mat, rho, u+h) - scvheosPofRhoU(Mat, rho, u))/h;

    return dPdU;
}

/*
 * Calculate u2 so that (rho1, u1) and (rho2, u2) are on the same isentrope.
 *
 * Returns SCVHEOS_FAIL_VALUE if any of the lookups fails.
 */
double scvheosIsentropicU(SCVHEOSMAT *Mat, double rho1, double u1, double rho2) {
    double T1;
    double T2;
    double S;
    double u2;

    T1 = scvheosTofRhoU(Mat, rho1, u1);
    if (scvheosFailed(T1)) return SCVHEOS_FAIL_VALUE;

    S = scvheosSofRhoT(Mat, rho1, T1);
    if (scvheosFailed(S)) return SCVHEOS_FAIL_VALUE;

    T2 = scvheosTofRhoS(Mat, rho2, S);
    if (scvheosFailed(T2)) return SCVHEOS_FAIL_VALUE;
   
    u2 = scvheosUofRhoT(Mat, rho2, T2);
    if (scvheosFailed(u2)) return SCVHEOS_FAIL_VALUE;

    return u2;
}


int scvheosPrintMat(SCVHEOSMAT *Mat, FILE *fp) {
    if (!fp) return SCVHEOS_FAIL;
	fprintf(fp,"# Material: %i (%s)\n", Mat->iMat, Mat->MatString);
	fprintf(fp,"# Reference density rho0: %g\n", Mat->rho0);
	fprintf(fp,"# Table size: nRho = %d, nT = %d\n", Mat->nRho, Mat->nT);
    fprintf(fp,"# Table boundaries (cgs units):\n");
    fprintf(fp,"# LogRhoMin=%15.7E LogRhoMax=%15.7E\n", Mat->dLogRhoAxis[0]+log10(Mat->dGmPerCcUnit), Mat->dLogRhoAxis[Mat->nRho-1]+log10(Mat->dGmPerCcUnit));
    fprintf(fp,"# LogTMin=%15.7E LogTMax=%15.7E\n", Mat->dLogTAxis[0], Mat->dLogTAxis[Mat->nT-1]);

    return SCVHEOS_SUCCESS;
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

/*
 * Check if a value (rho, u) is inside of the EOS table.
 *
 * Returns SCVHEOS_SUCCESS, or one of SCVHEOS_OUTSIDE_RHOMIN, _RHOMAX, _TMIN,
 * _TMAX. The energy bounds are u(rho, T_min) and u(rho, T_max) of the table
 * column, which assumes u(T) is monotonic at fixed rho (the same assumption
 * scvheosLogTofLogRhoLogU makes). A state that passes this test can be
 * inverted: both bracket ends are inside the table and the root is bracketed.
 *
 * rho <= 0 reports RHOMIN and u <= 0 reports TMIN: log10 of a non-positive
 * number is NaN or -inf, and the comparisons below would otherwise let the
 * state through. Non-finite input is reported the same way.
 */
int scvheosIsInExtrapLimit(SCVHEOSMAT *Mat, double rho, double u) {
    double logrho;
    double logu;

    if (!(rho > 0.0) || !isfinite(rho)) return SCVHEOS_OUTSIDE_RHOMIN;
    if (!(u > 0.0) || !isfinite(u)) return SCVHEOS_OUTSIDE_TMIN;

    logrho = log10(rho);
    logu = log10(u);

    if (logrho < Mat->LogRhoMin) return SCVHEOS_OUTSIDE_RHOMIN;
    if (logrho > Mat->LogRhoMax) return SCVHEOS_OUTSIDE_RHOMAX;
    if (logu < scvheosLogUofLogRhoLogT(Mat, logrho, Mat->LogTMin)) return SCVHEOS_OUTSIDE_TMIN;
    if (logu > scvheosLogUofLogRhoLogT(Mat, logrho, Mat->LogTMax)) return SCVHEOS_OUTSIDE_TMAX;

    return SCVHEOS_SUCCESS;
}

/* Functions required by the GSL root finder. */
double LogUofLogRhoLogT_GSL_rootfinder(double logT, void *params) {
    SCVHEOSMAT *Mat;
    double logrho;
    double logu;

    struct LogUofLogRhoLogT_GSL_Params *p;

    p = (struct LogUofLogRhoLogT_GSL_Params *) params;

    Mat = p->Mat;
    logrho = p->logrho;
    logu = p->logu;

    return (scvheosLogUofLogRhoLogT(Mat, logrho, logT)-logu);
}

double LogSofLogRhoLogT_GSL_rootfinder(double logT, void *params) {
    SCVHEOSMAT *Mat;
    double logrho;
    double logs;

    struct LogSofLogRhoLogT_GSL_Params *p;

    p = (struct LogSofLogRhoLogT_GSL_Params *) params;

    Mat = p->Mat;
    logrho = p->logrho;
    logs = p->logs;

    return (scvheosLogSofLogRhoLogT(Mat, logrho, logT)-logs);
}

double LogPofLogRhoLogT_GSL_rootfinder(double logrho, void *params) {
    SCVHEOSMAT *Mat;
    double logT;
    double logP;

    struct LogPofLogRhoLogT_GSL_Params *p;

    p = (struct LogPofLogRhoLogT_GSL_Params *) params;

    Mat = p->Mat;
    logT = p->logT;
    logP = p->logP;

    return (scvheosLogPofLogRhoLogT(Mat, logrho, logT)-logP);
}