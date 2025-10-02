/* CAUTION: This is the traditional K&R C (only) version of the Numerical
   Recipes utility file NR_util.h. */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define NR_SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define NR_DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define NR_DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define NR_DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define NR_FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define NR_FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define NR_LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define NR_LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define NR_IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define NR_IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define NR_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void NR_error();
float *NR_vector();
float **NR_matrix();
float **NR_submatrix();
float **NR_convert_matrix();
float ***NR_f3tensor();
double *NR_dvector();
double **NR_dmatrix();
int *NR_ivector();
int **NR_imatrix();
unsigned char *NR_cvector();
unsigned long *NR_lvector();
void NR_free_vector();
void NR_free_dvector();
void NR_free_ivector();
void NR_free_cvector();
void NR_free_lvector();
void NR_free_matrix();
void NR_free_submatrix();
void NR_free_convert_matrix();
void NR_free_dmatrix();
void NR_free_imatrix();
void NR_free_f3tensor();

#endif /* _NR_UTILS_H_ */
