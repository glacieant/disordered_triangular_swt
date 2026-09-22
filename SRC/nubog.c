#include "trilattice_sw.h"



//**** Numerical Bogoliubov routine *****//

// Linearized matrix array

/*

   NSYS = (2L)*(2L), Hamiltonian matrix size
   NL = 2L, matrix order
   M, filled Hamiltonian matrix
   T, diagonalizing matrix
   OMEGA, eigenvalues

*/

void nubog(int NL,double BTOL,complex *M,complex *T,double *OMEGA) {

	int i,j,k,l;

	// Hamiltonian matrix size

	int NSYS;
	NSYS=NL*NL;

	// S.M matrix, S = diag(1 ... 1 -1 ... -1) of size 2L

	complex *SM;
	SM=(complex *)malloc(NSYS*sizeof(complex));

	for (j=0;j<NL;j++) {

		for (i=0;i<NL/2;i++) {

			SM[i+j*NL]=M[i+j*NL];

		}
		for (i=NL/2;i<NL;i++) {

			SM[i+j*NL]=realmult(-1.0,M[i+j*NL]);

		}

	}

	// Right eigenvectors Z and eigenvalues LAMBDA of SM

	complex *Z;
	Z=(complex *)malloc(NSYS*sizeof(complex));

	int info1,lwork1;

	complex* work1;
	complex wkopt1;	

	double *rwork1;
	rwork1=(double *)malloc(2*NL*sizeof(double));

	complex *eigenval, *vl;
	eigenval=(complex *)malloc(NL*sizeof(complex));
	vl=(complex *)malloc(NSYS*sizeof(complex));


	lwork1=-1;
	zgeev_("N","V",&NL,SM,&NL,eigenval,vl,&NL,Z,&NL,&wkopt1,
			&lwork1,rwork1,&info1);
	lwork1 = (int)wkopt1.re;
	work1 = (complex *)malloc(lwork1*sizeof(complex));
	zgeev_("N","V",&NL,SM,&NL,eigenval,vl,&NL,Z,&NL,work1,
			&lwork1,rwork1,&info1);

	free(work1);
	free(rwork1);
	free(vl);

	complex *SZ;
	SZ=(complex *)malloc(NSYS*sizeof(complex));

	for (j=0;j<NL;j++) {

		for (i=0;i<NL/2;i++) {

			SZ[i+j*NL]=Z[i+j*NL];

		}
		for (i=NL/2;i<NL;i++) {

			SZ[i+j*NL]=realmult(-1.0,Z[i+j*NL]);

		}

	}

	complex *ELL;
	ELL=(complex *)malloc(NSYS*sizeof(complex));

	complex alpha={1.0,0.0};
	complex beta={0.0,0.0};

	zgemm_("C","N",&NL,&NL,&NL,&alpha,Z,&NL,SZ,&NL,&beta,ELL,&NL);

	int lwork2,info2;

	complex *work2;
	double *rwork2;
	complex wkopt2;

	double *reigenval;
	reigenval=(double *)malloc(NL*sizeof(double));

	rwork2=(double *)malloc((3*NL-2)*sizeof(double));

	lwork2 = -1;
	zheev_("V","U",&NL,ELL,&NL,reigenval,&wkopt2,&lwork2,rwork2,&info2);
	lwork2 = (int)wkopt2.re;
	work2 = (complex *)malloc(lwork2*sizeof(complex));
	zheev_("V","U",&NL,ELL,&NL,reigenval,work2,&lwork2,rwork2,&info2);
	free(work2);
	free(rwork2);


	// Bogoliubov transformation matrix
	complex *ELLD;
	ELLD=(complex *)malloc(NSYS*sizeof(complex));

#pragma omp parallel private(i,j)
	{
#pragma omp for
	for (i=0;i<NL;i++) {

		for(j=0;j<NL;j++) {

			ELLD[i+j*NL]=realmult(1.0/sqrt(fabs(reigenval[j])),ELL[i+j*NL]);
		}

	}
	}

	zgemm_("N","N",&NL,&NL,&NL,&alpha,Z,&NL,ELLD,&NL,&beta,T,&NL);

	// Hamiltonian eigenmodes, OMEGA = U^H.LAMBDA.U

#pragma omp parallel private(i,j)
	{
#pragma omp for
	for (i=0;i<NL;i++) {

		OMEGA[i]=0.0;
		for (j=0;j<NL;j++) {

			OMEGA[i]+=cmplxmult(eigenval[j],cmplxmult(cmplxconj(ELL[j+i*NL]),ELL[j+i*NL])).re;

		}
	}
	}

	// Zero modes with divergent norm
	
#pragma omp parallel private(i,j)
	{
#pragma omp for
	for (i=0;i<NL/2;i++) {

		if (fabs(OMEGA[i])<BTOL) {

			for (j=0;j<NL;j++) {

				T[j+i*NL].re=0.0;
				T[j+i*NL].im=0.0;
		
			}
		}

	}
	}

	free(ELLD);
	free(eigenval);
	free(reigenval);
	free(SM);
	free(Z);
	free(SZ);
	free(ELL);

}
