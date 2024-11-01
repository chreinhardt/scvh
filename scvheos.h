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
#define SCVHEOS_VERSION_TEXT      "1.0.1"
#define SCVHEOS_VERSION_MAJOR     1
#define SCVHEOS_VERSION_MINOR     0
#define SCVHEOS_VERSION_PATCH     1

/*
 * Error codes.
 */
#define SCVHEOS_SUCCESS 0
#define SCVHEOS_FAIL   -1

/*
 * Material ids.
 */
#define SCVHEOS_H               110
#define SCVHEOS_HE              111
#define SCVHEOS_HHE             112
#define SCVHEOS_HHE_LOWRHOT     113
#define SCVHEOS_HHE_EXT_LOWRHOT 114

/*
 * Define error codes.
 */
#define SCVHEOS_OUTSIDE_RHOMIN -1
#define SCVHEOS_OUTSIDE_RHOMAX -2
#define SCVHEOS_OUTSIDE_TMIN   -3
#define SCVHEOS_OUTSIDE_TMAX   -4

/*
 * Define TRUE and FALSE.
 */
#define FALSE 0
#define TRUE  1

typedef struct scvheosMat {
    int iMat;
    int nRho;
    int nT;
    char MatString[256];
    double rho0;

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

/* Data structure needed for the GSL root finder. */
struct LogUofLogRhoLogT_GSL_Params {
    SCVHEOSMAT *Mat;
    double logrho;
    double logu;
};

struct LogSofLogRhoLogT_GSL_Params {
    SCVHEOSMAT *Mat;
    double logrho;
    double logs;
};

struct LogPofLogRhoLogT_GSL_Params {
    SCVHEOSMAT *Mat;
    double logT;
    double logP;
};

/* Functions to initialize and finalize materials. */
SCVHEOSMAT *scvheosInitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat);
int scvheosReadTable(SCVHEOSMAT *Mat, char *chInFile, int nRho, int nT, int nSkip);
int scvheosGenerateSoundSpeedTable(SCVHEOSMAT *Mat);

/* Calculate EOS values. */
double scvheosLogPofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosLogUofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosLogSofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosLogCsofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosSofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u);

/* Calculate derivatives. */
double scvheosdLogPdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosdLogPdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosdLogUdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosdLogUdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosdLogSdLogRhoofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosdLogSdLogTofLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);

double scvheosLogTofLogRhoLogU(SCVHEOSMAT *Mat, double logrho, double logu);
double scvheosLogTofLogRhoLogS(SCVHEOSMAT *Mat, double logrho, double logs);
double scvheosLogRhoofLogPLogT(SCVHEOSMAT *Mat, double logP, double logT);

double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosRhoofPT(SCVHEOSMAT *Mat, double P, double T);

/* Use df/dx = f/x*ln(f)/ln(x). */
double scvheosdPdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdUdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdUdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdSdRhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdSdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);

double scvheosdPdRhoofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosdPdUofRhoU(SCVHEOSMAT *Mat, double rho, double u);

double scvheosIsentropicU(SCVHEOSMAT *Mat, double rho1, double u1, double rho2);

int scvheosPrintMat(SCVHEOSMAT *Mat, FILE *fp);

int scvheosCheckTableBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
int scvheosCheckBoundsLogRhoLogT(SCVHEOSMAT *Mat, double logrho, double logT);
int scvheosIsInExtrapLimit(SCVHEOSMAT *Mat, double rho, double u);

/* Functions required by the GSL root finder. */
double LogUofLogRhoLogT_GSL_rootfinder(double logT, void *params);
double LogSofLogRhoLogT_GSL_rootfinder(double logT, void *params);
double LogPofLogRhoLogT_GSL_rootfinder(double logrho, void *params);
#endif

