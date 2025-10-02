/* Use the preprocessor to include IONOSPHERE header file IONOSPHERE_util.h. */

#include "IONOSPHERE_util.h"
#define IONOSPHERE_END 1
#define FREE_ARG char*

/*******************************************************************/

void IONOSPHERE_error(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
        void exit();

        fprintf(stderr,"\nNumerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

float *IONOSPHERE_vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+IONOSPHERE_END)*sizeof(float)));
        if (!v) IONOSPHERE_error("allocation failure in vector()");
        return v-nl+IONOSPHERE_END;
}

int *IONOSPHERE_ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+IONOSPHERE_END)*sizeof(int)));
        if (!v) IONOSPHERE_error("allocation failure in ivector()");
        return v-nl+IONOSPHERE_END;
}

unsigned char *IONOSPHERE_cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
        unsigned char *v;

        v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+IONOSPHERE_END)*sizeof(unsigned char)));
        if (!v) IONOSPHERE_error("allocation failure in cvector()");
        return v-nl+IONOSPHERE_END;
}

unsigned long *IONOSPHERE_lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
        unsigned long *v;

        v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+IONOSPHERE_END)*sizeof(long)));
        if (!v) IONOSPHERE_error("allocation failure in lvector()");
        return v-nl+IONOSPHERE_END;
}

double *IONOSPHERE_dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((unsigned int) ((nh-nl+1+IONOSPHERE_END)*sizeof(double)));
        if (!v) IONOSPHERE_error("allocation failure in dvector()");
        return v-nl+IONOSPHERE_END;
}

float **IONOSPHERE_matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((unsigned int)((nrow+IONOSPHERE_END)*sizeof(float*)));
        if (!m) IONOSPHERE_error("allocation failure 1 in matrix()");
        m += IONOSPHERE_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+IONOSPHERE_END)*sizeof(float)));
        if (!m[nrl]) IONOSPHERE_error("allocation failure 2 in matrix()");
        m[nrl] += IONOSPHERE_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

double **IONOSPHERE_dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((unsigned int)((nrow+IONOSPHERE_END)*sizeof(double*)));
        if (!m) IONOSPHERE_error("allocation failure 1 in matrix()");
        m += IONOSPHERE_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+IONOSPHERE_END)*sizeof(double)));
        if (!m[nrl]) IONOSPHERE_error("allocation failure 2 in matrix()");
        m[nrl] += IONOSPHERE_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

int **IONOSPHERE_imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;

        /* allocate pointers to rows */
        m=(int **) malloc((unsigned int)((nrow+IONOSPHERE_END)*sizeof(int*)));
        if (!m) IONOSPHERE_error("allocation failure 1 in matrix()");
        m += IONOSPHERE_END;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+IONOSPHERE_END)*sizeof(int)));
        if (!m[nrl]) IONOSPHERE_error("allocation failure 2 in matrix()");
        m[nrl] += IONOSPHERE_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

float **IONOSPHERE_submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
        long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
        float **m;

        /* allocate array of pointers to rows */
        m=(float **) malloc((unsigned int) ((nrow+IONOSPHERE_END)*sizeof(float*)));
        if (!m) IONOSPHERE_error("allocation failure in submatrix()");
        m += IONOSPHERE_END;
        m -= newrl;

        /* set pointers to rows */
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

float **IONOSPHERE_convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
long nch,ncl,nrh,nrl;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((unsigned int) ((nrow+IONOSPHERE_END)*sizeof(float*)));
        if (!m) IONOSPHERE_error("allocation failure in convert_matrix()");
        m += IONOSPHERE_END;
        m -= nrl;

        /* set pointers to rows */
        m[nrl]=a-ncl;
        for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
        /* return pointer to array of pointers to rows */
        return m;
}

float ***IONOSPHERE_f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
        float ***t;

        /* allocate pointers to pointers to rows */
        t=(float ***) malloc((unsigned int)((nrow+IONOSPHERE_END)*sizeof(float**)));
        if (!t) IONOSPHERE_error("allocation failure 1 in f3tensor()");
        t += IONOSPHERE_END;
        t -= nrl;

        /* allocate pointers to rows and set pointers to them */
        t[nrl]=(float **) malloc((unsigned int)((nrow*ncol+IONOSPHERE_END)*sizeof(float*)));
        if (!t[nrl]) IONOSPHERE_error("allocation failure 2 in f3tensor()");
        t[nrl] += IONOSPHERE_END;
        t[nrl] -= ncl;

        /* allocate rows and set pointers to them */
        t[nrl][ncl]=(float *) malloc((unsigned int)((nrow*ncol*ndep+IONOSPHERE_END)*sizeof(float)));
        if (!t[nrl][ncl]) IONOSPHERE_error("allocation failure 3 in f3tensor()");
        t[nrl][ncl] += IONOSPHERE_END;
        t[nrl][ncl] -= ndl;

        for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
        for(i=nrl+1;i<=nrh;i++) {
                t[i]=t[i-1]+ncol;
                t[i][ncl]=t[i-1][ncl]+ncol*ndep;
                for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

        /* return pointer to array of pointers to rows */
        return t;
}

void IONOSPHERE_free_vector(v,nl,nh)
float *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-IONOSPHERE_END));
}

void IONOSPHERE_free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-IONOSPHERE_END));
}

void IONOSPHERE_free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
        free((FREE_ARG) (v+nl-IONOSPHERE_END));
}

void IONOSPHERE_free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
        free((FREE_ARG) (v+nl-IONOSPHERE_END));
}

void IONOSPHERE_free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-IONOSPHERE_END));
}

void IONOSPHERE_free_matrix(m,nrl,nrh,ncl,nch)
float **m;
long nch,ncl,nrh,nrl;
/* free a float matrix allocated by matrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-IONOSPHERE_END));
        free((FREE_ARG) (m+nrl-IONOSPHERE_END));
}

void IONOSPHERE_free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-IONOSPHERE_END));
        free((FREE_ARG) (m+nrl-IONOSPHERE_END));
}

void IONOSPHERE_free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-IONOSPHERE_END));
        free((FREE_ARG) (m+nrl-IONOSPHERE_END));
}

void IONOSPHERE_free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
        free((FREE_ARG) (b+nrl-IONOSPHERE_END));
}

void IONOSPHERE_free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
        free((FREE_ARG) (b+nrl-IONOSPHERE_END));
}

void IONOSPHERE_free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
float ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a float f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-IONOSPHERE_END));
        free((FREE_ARG) (t[nrl]+ncl-IONOSPHERE_END));
        free((FREE_ARG) (t+nrl-IONOSPHERE_END));
}

/*******************************************************************

 IONOSPHERE_read_line: This subroutine reads in a line of text from the
                       standard input file or terminal.

 Variable description:

 (call statement parameter list)

 buffer              Line of text.

 Begin subroutine IONOSPHERE_read_line. */

void IONOSPHERE_read_line(buffer)

/* Variable declaration. */

char buffer[];
{

  /* Local variable declaration. */

  int i;
  char single_char;

  /* Get the character string using getchar. */

  i=0;
  do
    {
       single_char = getchar();
       buffer[i] = single_char;
       ++i;
    }
  while(single_char != '\n');

  buffer[i-1] = '\0';

  return;

/* End subroutine IONOSPHERE_read_line. */

}

/*******************************************************************

 IONOSPHERE_concatenate: This subroutine concatenates two character strings.

 Variable description:

 (call statement parameter list)

 string1, string2          Input character strings.

 result                    Result of concatenation.

 Begin subroutine IONOSPHERE_concatenate. */

void IONOSPHERE_concatenate(string1, string2, result)

/* Variable declaration. */

char string1[], string2[], result[];
{

  /* Local variable declaration. */

  int i, j;

  /* Copy string1 to result. */

  for ( i = 0; string1[i] != '\0'; ++i )
    result[i] = string1[i];

  /* Copy string2 to result. */

  for ( j = 0; string2[j] != '\0'; ++j )
    result[i+j] = string2[j];

  /* Terminate the concatenated string with a null. */

  result[i+j] = '\0';

  return;

/* End subroutine IONOSPHERE_concatenate. */

}

/*******************************************************************

 IONOSPHERE_copy: This subroutine copies one character string to another.

 Variable description:

 (call statement parameter list)

 string_from, string_to          Character strings.

 Begin subroutine IONOSPHERE_copy. */

void IONOSPHERE_copy(string_from, string_to)

/* Variable declaration. */

char string_from[], string_to[];
{

  /* Local variable declaration. */

  int i;

  /* Copy string_from to string_to. */

  for ( i = 0; string_from[i] != '\0'; ++i )
    string_to[i] = string_from[i];

  string_to[i] = '\0';

  return;

/* End subroutine IONOSPHERE_copy. */

}

/*******************************************************************

 IONOSPHERE_compare: This function compares two character strings
                     and determines whether they are equal or not,
                     returning 1 (TRUE) if in fact the two strings
                     are identical and 0 (FALSE) if they are not.

 Variable description:

 (call statement parameter list)

 string1, string2          Input character strings.

 Begin subroutine IONOSPHERE_compare. */

int IONOSPHERE_compare(string1, string2)

/* Variable declaration. */

char string1[], string2[];
{

  /* Local variable declaration. */

  int i=0, answer;

  /* Compare string1 and string2 */

  while ( string1[i] == string2[i] &&
          string1[i] != '\0' && string2[i] != '\0' )
     ++i;

  if (string1[i] == '\0' && string2[i] == '\0') {
     answer = 1;
  } else {
     answer = 0;
  } /* endif */

  return(answer);

/* End function IONOSPHERE_compare. */

}
