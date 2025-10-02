/* IONOSPHERE_util.h:  Header file for IONOSPHERE utility subroutines. */

/*==== s y s t e m  r o u t i n e s  &  d e f i n i t i o n s ====*/

/* Use the preprocessor to include the math, I/O, and
   system baseline libraries and subroutines */

#include <stdio.h>
#include <math.h>

/*==== r o u t i n e s ===========================================*/

/* Use the preprocessor to define various application routines. */

static float sqrarg;
#define IONOSPHERE_SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define IONOSPHERE_DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define IONOSPHERE_DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define IONOSPHERE_DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define IONOSPHERE_FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define IONOSPHERE_FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define IONOSPHERE_LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define IONOSPHERE_LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IONOSPHERE_IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IONOSPHERE_IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define IONOSPHERE_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void IONOSPHERE_error();
float *IONOSPHERE_vector();
float **IONOSPHERE_matrix();
float **IONOSPHERE_submatrix();
float **IONOSPHERE_convert_matrix();
float ***IONOSPHERE_f3tensor();
double *IONOSPHERE_dvector();
double **IONOSPHERE_dmatrix();
int *IONOSPHERE_ivector();
int **IONOSPHERE_imatrix();
unsigned char *IONOSPHERE_cvector();
unsigned long *IONOSPHERE_lvector();
void IONOSPHERE_free_vector();
void IONOSPHERE_free_dvector();
void IONOSPHERE_free_ivector();
void IONOSPHERE_free_cvector();
void IONOSPHERE_free_lvector();
void IONOSPHERE_free_matrix();
void IONOSPHERE_free_submatrix();
void IONOSPHERE_free_convert_matrix();
void IONOSPHERE_free_dmatrix();
void IONOSPHERE_free_imatrix();
void IONOSPHERE_free_f3tensor();

void IONOSPHERE_read_line();
void IONOSPHERE_concatenate();
void IONOSPHERE_copy();
int IONOSPHERE_compare();

