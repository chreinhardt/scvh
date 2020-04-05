/*
 * Header file for scvheos.c.
 *
 * Author:   Christian Reinhardt
 * Created:  02.04.2020
 * Modified: 05.04.2020
 */
#ifndef SCVHEOS_HINCLUDED
#define SCVHEOS_HINCLUDED

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
#define SCVHEOS_H    110
#define SCVHEOS_HE   111
#define SCVHEOS_HHE  112

/*
 * Define TRUE and FALSE.
 */
#define FALSE 0
#define TRUE  1

typedef struct scvheosMat {
    int iMat;
    int nRho;
    int nT;

    /* Convert units to code units. */
    double dKpcUnit;
    double dMsolUnit;
    double dGasConst;
	double dErgPerGmUnit;
	double dGmPerCcUnit;
	double dSecUnit;

    /* The EOS tables are logarithmic.*/
    double *dLogRhoAxis;
    double *dLogTAxis;
    double **dLogUArray;
    double **dLogPArray;
    double **dLogSArray;
    double **dLogCArray;
    
    /* Base of the logarith. */
    double dLogBase;
} SCVHEOSMAT;

/* Functions to initialize and finalize materials. */
SCVHEOSMAT *scvheosInitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
int scvheosFinalizeMaterial(SCVHEOSMAT *Mat);
int scvheosReadTable(SCVHEOSMAT *Mat, char *chInFile, int nRho, int nT);

/* Calculate EOS values. */
double scvheosLogPofRhoT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosPofRhoT(SCVHEOSMAT *Mat, double rho, double T);
//double scvheosPofRhoU(SCVHEOSMAT *Mat, double rho, double u);

double scvheosLogUofRhoT(SCVHEOSMAT *Mat, double logrho, double logT);
double scvheosUofRhoT(SCVHEOSMAT *Mat, double rho, double T);

# if 0
double scvheosCsofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosCsofRhoU(SCVHEOSMAT *Mat, double rho, double u);
double scvheosTofRhoU(SCVHEOSMAT *Mat, double rho, double u);

double scvheosdPdrhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdPdTofRhoT(SCVHEOSMAT *Mat, double rho, double T);
double scvheosdUdrhoofRhoT(SCVHEOSMAT *Mat, double rho, double T);
#endif
#endif

