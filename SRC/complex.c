#include "trilattice_sw.h"

//****Complex algebra routines***//


complex cmplxadd(complex A,complex B) {

	complex C;
	C.re=A.re+B.re;
	C.im=A.im+B.im;

	return C;

}

complex cmplxmult(complex A,complex B) {

	complex C;
	C.re=A.re*B.re-A.im*B.im;
	C.im=A.re*B.im+A.im*B.re;

	return C;

}


complex realmult(double A,complex B) {

	complex C;
	C.re=A*B.re;
	C.im=A*B.im;

	return C;

}

complex cmplxconj(complex A) {

	complex C;
	C.re=A.re;
	C.im=-A.im;

	return C;
}

double cmplxmod(complex A) {

	return sqrt(A.re*A.re+A.im*A.im);

}

double cmplxmodsqr(complex A) {

	return A.re*A.re+A.im*A.im;

}
//*****Complex square matrix printer****//

void printzmat(int SIZE, complex *MAT) {

	int i, j;

	int DIM=sqrt(SIZE);

	FILE *fyl;

	char ofname[200];
	snprintf(ofname,200*sizeof(char),"MAT_%d",(int)time(NULL));

	fyl=fopen(ofname,"a");

	fprintf(fyl,"\n");

	for (i=0;i<DIM;i++) {

		for (j=0;j<DIM;j++) {

			fprintf(fyl,"(%lf,%lf)\t",MAT[i+j*DIM].re,MAT[i+j*DIM].im);

		}
		fprintf(fyl,"\n");
	}

	fclose(fyl);

}

