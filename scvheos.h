/*
 * Header file for scvheos.c.
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 05.04.2020
 */
#ifndef SCVHEOS_HINCLUDED
#define SCVHEOS_HINCLUDED

#include <gsl/gsl_interp2d.h>

/*
 * Version.
 */
#define SCVHEOS_VERSION_TEXT      "1.0.0"
#define SCVHEOS_VERSION_MAJOR     1
#define SCVHEOS_VERSION_MINOR     0
#define SCVHEOS_VERSION_PATCH     0

/*
 * Error codes.
 */
#define SCVHEOS_SUCCESS 0
#define SCVHEOS_FAIL   -1

/*
 * Material ids.
 */
#define SCVHEOS_H           110
#define SCVHEOS_HE          111
#define SCVHEOS_HHE         112
#define SCVHEOS_HHE_LOWRHOT 113

/*
 * Define TRUE and FALSE.
 */
#define FALSE 0
#define TRUE  1

typedef struct scvheosMat {
    int iMat;
    int nRho;
    int nT;

    /* Define limits of the table (or extrapolation). */
    double LogRhoMin;
    double LogRhoMax;
    double LogTMin;
    double LogTMax;

    /* Convert units to code units. */
    double dKpcUnit;
    double dMsolUnit;
    double dGasConst;
	double dErgPerGmUnit;
	double dGmPerCcUnit;
	double dSecUnit;

    /* GSL 2D interpolator */
    const gsl_interp2d_type *InterpType;
    gsl_interp2d *InterpLogP;
    gsl_interp2d *InterpLogU;
    gsl_interp2d *InterpLogS;
    gsl_interp2d *InterpLogCs;

    /* Accelerator */
    gsl_interp_accel *xAccP;
    gsl_interp_accel *yAccP;
    gsl_interp_accel *xAccU;
    gsl_interp_accel *yAccU;
    gsl_interp_accel *xAccS;
    gsl_interp_accel *yAccS;
    gsl_interp_accel *xAccCs;
    gsl_interp_accel *yAccCs;

    /* The EOS tables are logarithmic.*/
    double *dLogRhoAxis;
    double *dLogTAxis;
    double *dLogUArray;
    double *dLogPArray;
    double *dLogSArray;
    double *dLogCArray;
    
    /* Base of the logarith. */
    double dLogBase;
} SCVHEOSMAT;

/* Functions to initialize and finalize materials. */
SCVHEOSMAT *scvheosInitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat);
int scvheosReadTable(SCVHEOSMAT *Mat, char *chInFile, int nRho, int nT, int nSkip);

/* Calculate EOS values. */
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T);
//double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u);

double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosSofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdSdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdSdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);


int scvheosCheckTableBoundsLogRhoLogT(REOS3MAT *Mat, double logrho, double logT);
int scvheosCheckBoundsLogRhoLogT(REOS3MAT *Mat, double logrho, double logT);

# if 0
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u);

double scvheosdPdrhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdUdrhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
#endif
#endif

