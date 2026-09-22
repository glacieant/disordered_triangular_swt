#ifndef LIBRARY_INCLUSION
#define LIBRARY_INCLUSION

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <fftw3.h>
#include "mt64.h"

#endif

// Lattice constant

#define latconst 1.0

// Global zero tolerance

#define gtol (1.0e-8)*1.0

// pi
#define pi 3.14159265358979323846

// Spin magnitude S
#define S 0.5

// Lattice co-ordination number
#define ZCO 6

#ifndef COMPLEX_STRUCT
#define COMPLEX_STRUCT

typedef struct _complex{

	double re;
	double im;

} complex;

#endif

// Three-component vector
#ifndef VECTOR_STRUCT
#define VECTOR_STRUCT

typedef struct _vect{
	double x;
	double y;
	double z;
} vect;

#endif

// Three-component complex vector

#ifndef COMPLEX_VECTOR_STRUCT
#define COMPLEX_VECTOR_STRUCT

typedef struct _cmplxvect{
	complex x;
	complex y;
	complex z;
} cmplxvect;

#endif

// 3D solid angle
#ifndef ANGLE3D_STRUCT
#define ANGLE3D_STRUCT

typedef struct _angle3D{
	double th;
	double ph;
} angle3D;

#endif

// Lattice point structure

#ifndef LATTICE_STRUCT
#define LATTICE_STRUCT

typedef struct _lattpoint{

	int *neighbour;
	int *nn_nbr;
	vect pos;
	angle3D ang;

} lattpoint;

#endif



//Complex algebra routines//


complex cmplxadd(complex A,complex B);
complex cmplxmult(complex A,complex B);
complex realmult(double A,complex B);
complex cmplxconj(complex A);
double cmplxmod(complex A);
double cmplxmodsqr(complex A);

//Complex square matrix printer//

void printzmat(int SIZE, complex *MAT);
void printmat(int SIZE, double *MAT);


//External linear algebra routines//

extern void zheev_(char* jobz, char* uplo, int* n,
		complex* a, int* lda,double* w,
		complex* work, int* lwork, double* rwork, 
		int* info );

extern void zgemm_(char *transa, char *transb, int *m, int *n, 
		int *k,complex *alpha,complex *A,int *lda, 
		complex *B,int *ldb,complex *beta,complex *C,
		int *ldc);

extern void dgemm_(char *transa, char *transb, int *m, int *n, 
		int *k,double *alpha,double  *A,int *lda, 
		double *B,int *ldb,double *beta,double *C,
		int *ldc);

extern void dgemv_(char *trans,int *m, int *n, 
		double *alpha,double  *A,int *lda, 
		double *X,int *incx,double *beta,double *Y,
		int *incy);

extern void zgeev_(char *jobvl,char *jobvr,int *n,complex *a,
		int *lda,complex *w, complex *vl, int *ldvl,
		complex *vr, int *ldvr,complex *work,int *lwork,
		double *rwork, int *info);

extern void dsyev_(char *jobz,char *uplo,int *n,double *a,
		int *lda,double *w,double *work, 
		int *lwork, int *info);

extern void dgesvd_(char *jobu,char *jobvt,int *m,int *n,
		double *a,int *lda,double *sval,double *u,
		int *ldu, double *vt, int *ldvt,
		double *work, int *lwork, int *info);

//Classical algorithm and lattice routines//

double Fxx(int i, int j, double **J,lattpoint *lsite);
double Fyy(int i, int j, double **J,lattpoint *lsite);
double Fzz(int i, int j, double **J,lattpoint *lsite);
double Fxy(int i, int j, double **J,lattpoint *lsite);
complex A(int i, int j, double **J,lattpoint *lsite);
complex B(int i, int j, double **J,lattpoint *lsite); 


double Fxz(int i, int j, double **J,lattpoint *lsite);
double Fyz(int i, int j, double **J,lattpoint *lsite);
complex dthj_Fz(int i, int j, double **J,lattpoint *lsite);
complex dphj_Fz(int i, int j, double **J,lattpoint *lsite);
complex dthi_Fz(int i, int j, double **J,lattpoint *lsite);
complex dphi_Fz(int i, int j, double **J,lattpoint *lsite);

void gen_ham(int MDIM,lattpoint *lsite,double **J,complex *M); 

void shuffle(int diter,int PARX,int ORDX,int *array,int n);

void classic_algo(int diter,int PARX,lattpoint *lsite,double **J,
		int N_lattpoint,int iter_NUMG);

void create_trilattice(int N_lattpoint, lattpoint *lsite);
void create_auxlattice(int N_lattpoint,int **sublattice, int **eltriangle);

void gen_angles(int N_lattpoint,double *thetalist,
		double *philist,lattpoint *lsite);
void gen_couplings(int N_lattpoint,lattpoint *lsite,
		double ALPHA,double DELTA,double **J);
int subcycle(int i);

void copylattice(int N_lattpoint,lattpoint *lat_SRC,lattpoint *lat_DEST);
void copycoupling(int N_lattpoint,double **J_SRC,double **J_DEST);
FILE *openfile(const char *name, const char *mode);
void readconstant(FILE *fyl, const char *KEY, char *Value);
int readint(FILE *fyl, const char *KEY);
double readdouble(FILE *fyl, const char *KEY);

//Random list generator//

void randgen_cpl(int N_lattpoint,int NMX,int diter,
		double delta,double *randlist);

void randgen_ang(int N_lattpoint,int NMX,int par,
		double ANGDISPAR,double *thetalist,double *philist);

//Numerical Bogoliubov routine//

void nubog(int NL,double BTOL,complex *M,complex *T,double *OMEGA);

//Vector algebra routines//

double vectnorm(vect R);
double vectnormsqr(vect R);
vect vectmult(double c,vect R);
vect unitvect(vect R);
vect vectsum(vect R1,vect R2);
double vectdot(vect R1,vect R2);
vect vectcross(vect R1,vect R2);

//Complex vector algebra routine//

cmplxvect cmplxconjvect(cmplxvect X);
cmplxvect cmplxvectcross(cmplxvect R1,cmplxvect R2);
cmplxvect cmplxvectsum(cmplxvect X,cmplxvect Y);
double cmplxvectnorm(cmplxvect X);
double cmplxvectnormsqr(cmplxvect X);
cmplxvect re2cmplxvect(vect X);
cmplxvect realvectmult(double C,cmplxvect X);
complex cmplxvectdot(cmplxvect X,cmplxvect Y);
cmplxvect unitcmplxvect(cmplxvect R);

//Observables//

double SYSAVG(double X,int NSMPLE);
double SYSERR(double X,double XSQR,int NSMPLE);

void ground_en(int N_lattpoint,lattpoint *lsite,double **J,
		double *OMEGA,double *EN,double *EN_var);

void cl_ground_en(int N_lattpoint,lattpoint *lsite,double **J,
		double *cl_EN,double *cl_EN_var);

void stag_spin(int N_lattpoint, lattpoint *lsite, complex **alpha,
		double *OMEGA,double *SPIN,double *SPIN_var,
		double *spin_grid,double *spin_grid_var);

void qsubmag(int N_lattpoint,lattpoint *lsite,
		int **sublattice,double *dthph,
		complex *T,double *OMEGA,complex **alpha,
		double *SBLMAG,double *SBLCHI,
		double *SBLMAG_VAR,double *SBLCHI_VAR);

void cl_sublatmag(int N_lattpoint,lattpoint *lsite,
		int **sublattice,double *CLSBL_ANG,
		double *CLSBL_CHI,double *CLSBL_ANG_VAR,
		double *CLSBL_CHI_VAR);

void tri_cl_sublatmag(int N_lattpoint,lattpoint *lsite, 
		int **eltriangle,double *tri_CLSBL_ANG, 
		double *tri_CLSCAL_CHI,double *tri_CLSBL_CHI,
		double *tri_CLSBL_ANG_VAR, 
		double *tri_CLSCAL_CHI_VAR,
		double *tri_CLSBL_CHI_VAR);

void tri_qsubmag_ex(int N_lattpoint,
		lattpoint *lsite,
		int **eltriangle,double *dthph,
		complex *T,double *OMEGA,
		complex **alpha,
		double *tri_SBLMAG_ex,
		double *tri_SCALCHI_ex,
		double *tri_SBLCHI_ex,
		double *tri_SBLMAG_ex_VAR,
		double *tri_SCALCHI_ex_VAR,
		double *tri_SBLCHI_ex_VAR);

void im_chi(double womega,vect q,double GW,
		int N_lattpoint,lattpoint *lsite,
		double *dthph,complex *T,double *OMEGA,
		double *SUSC,double *SUSC_VAR,
		double *CUSC,double *CSUSC_VAR);

void CHI_MAT(int N_lattpoint,lattpoint *lsite,
		double *dthph,complex *T,double *OMEGA,
		double ***CHIQ,double **SZERO);

void eigenview(int EIGENCOUNT,int N_lattpoint,
		complex *T,double *OMEGA,double *eigen_grid,double *eigen_grid_var);

void IPR_CALC(double IPRMIN,double IPRMAX,int IPRGRID,double GW,
		double BTOL,int N_lattpoint,complex *T,
		double *OMEGA,double *IPR,double *IPR_var);

void histo_angles(int N_lattpoint,lattpoint *lsite,double *dthph,
		double *cangle,double *cangle_var,
		double *qangle,double *qangle_var);

//Utilities for the main code//

void set_boundary(double XQMIN, double YQMIN,double XQMAX,
		double YQMAX, double WMIN, double WMAX,
		int FGRID, int XMGRID, int YMGRID,
		int L,
		int *XMBASE,int *YMBASE,int *WRANGE, int *XMRANGE,
		int *YMRANGE, int *dx, int *dy);


void calc_strfact(int STRSWITCH,double GW, lattpoint *lsite,
		double *dthph,complex *T,double *OMEGA,
		double WMIN, double WMAX, double DELTA,
		double ***IMCHI,double ***IMCHI_var, 
		int HCUT, int VCUT,
		int L,int WRANGE,int FGRID,
		int XMBASE, int YMBASE,
		int XMRANGE, int YMRANGE, int dx, int dy,
		int RANGEMIN, int RANGEMAX, int drmx);

void print_strfact(FILE *outgrid, int STRSWITCH,
		double WMIN, double WMAX, 
		double DELTA, double ALPHA,
		double ***IMCHI,double ***IMCHI_var, 
		int HCUT, int VCUT,
		int L,int WRANGE,int FGRID,
		int XMBASE, int YMBASE,
		int XMRANGE, int YMRANGE, int dx, int dy,
		int RANGEMIN, int RANGEMAX, int drmx,
		int iter_DISD);

void mf_params(int N_lattpoint,lattpoint *lsite,
		complex *T,complex **alpha,complex **beta);

void mf_angles(int N_lattpoint,lattpoint *lsite,double **J,
		complex **alpha,complex **beta,double BTOL,double *dthph);



void cl_moment(int N_lattpoint, lattpoint *lsite, 
        double *CMOM_X,
        double *CMOM_Y,
        double *CMOM_Z,
        double *CMOM_var);

void q_moment(int N_lattpoint, lattpoint *lsite, 
        complex **alpha,double *OMEGA,
        double *dthph,
        double *QMOM_X,
        double *QMOM_Y,
        double *QMOM_Z,
        double *QMOM_var);
