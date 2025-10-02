/* IONOSPHERE.h:  Header file for IONOSPHERE subroutines. */

/*==== d a t a  a n d  c o n s t a n t s =========================*/

/* Use the preprocessor to define various constants and data. */

#define IONOSPHERE_numberofnodes 17

#define IONOSPHERE_Radius_Earth 6378000.00
#define IONOSPHERE_Height_Earth 400000.00
#define IONOSPHERE_Density_Earth 1.33808e-20
#define IONOSPHERE_SoundSpeed_Earth 50000.00

#define IONOSPHERE_MU 1.256637e-06
#define IONOSPHERE_PI 3.141592654
#define IONOSPHERE_Theta_0 0.00001

#define IONOSPHERE_NPRE 4
#define IONOSPHERE_NPOST 4
#define IONOSPHERE_NEXACT 10
#define IONOSPHERE_NGMAX 15
#define IONOSPHERE_OMEGA 1.00
#define IONOSPHERE_TOLER 5.0e-05

/* Use the preprocessor to define various application routines. */

/*==== r o u t i n e s ===========================================*/

/* Use the preprocessor to define various application routines. */

void IONOSPHERE_conductance();
void IONOSPHERE_multigrid();
void IONOSPHERE_fgrid();
void IONOSPHERE_cgrid();
void IONOSPHERE_restrict();
void IONOSPHERE_prolong();
void IONOSPHERE_smoother();
void IONOSPHERE_residual();
void IONOSPHERE_update();
void IONOSPHERE_matcopy();
void IONOSPHERE_matzero();
void IONOSPHERE_magfac();
void IONOSPHERE_magvel();
void IONOSPHERE_current();

