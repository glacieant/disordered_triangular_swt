#include "trilattice_sw.h"

//****The vector algebra routines*****//


//Norm of a vector
double vectnorm(vect R) {

	double norm_val;
	norm_val=sqrt(R.x*R.x+R.y*R.y+R.z*R.z);
	return norm_val;

}

//Norm of a vector
double vectnormsqr(vect R) {

	double norm_val;
	norm_val=R.x*R.x+R.y*R.y+R.z*R.z;
	return norm_val;

}
//Scalar multiplication to a vector
vect vectmult(double c,vect R) {

	vect Rprime;
	Rprime.x=c*R.x;
	Rprime.y=c*R.y;
	Rprime.z=c*R.z;
	return Rprime;

}

//Unit vector from a cartesian vector
vect unitvect(vect R) {

	vect u={0.0,0.0,0.0};
	double X=vectnorm(R);
	if (X>=gtol) {
		u=vectmult(1.0/X,R);
	}
	return u;
}

//Summing up vectors
vect vectsum(vect R1,vect R2) {

	vect R;
	R.x=R1.x+R2.x;
	R.y=R1.y+R2.y;
	R.z=R1.z+R2.z;
	return R;

}

//Dot product of vector
double vectdot(vect R1,vect R2) {

	double dotprod;
	dotprod=R1.x*R2.x+R1.y*R2.y+R1.z*R2.z;
	return dotprod;

}

//Cross product of a vector
vect vectcross(vect R1, vect R2) {

	vect crossprod;
	crossprod.x=R1.y*R2.z-R1.z*R2.y;
	crossprod.y=-(R1.x*R2.z-R1.z*R2.x);
	crossprod.z=R1.x*R2.y-R1.y*R2.x;

	return crossprod;

}

cmplxvect cmplxconjvect(cmplxvect X) {

	cmplxvect Y;
	Y=X;
	Y.x.im=-Y.x.im;
	Y.y.im=-Y.y.im;
	Y.z.im=-Y.z.im;

	return Y;

}

cmplxvect cmplxvectsum(cmplxvect X,cmplxvect Y) {

	cmplxvect Z;

	Z.x=cmplxadd(X.x,Y.x);
	Z.y=cmplxadd(X.y,Y.y);
	Z.z=cmplxadd(X.z,Y.z);

	return Z;

}

complex cmplxvectdot(cmplxvect X,cmplxvect Y) {

	complex Z={0.0,0.0};

	Z=cmplxadd(Z,cmplxmult(X.x,Y.x));
	Z=cmplxadd(Z,cmplxmult(X.y,Y.y));
	Z=cmplxadd(Z,cmplxmult(X.z,Y.z));

	return Z;

}

double cmplxvectnorm(cmplxvect X) {

	double NORM;
	NORM=cmplxmodsqr(X.x)
		+cmplxmodsqr(X.y)+cmplxmodsqr(X.z);

	NORM=sqrt(NORM);
	return NORM;

}

double cmplxvectnormsqr(cmplxvect X) {

	double NORM;
	NORM=cmplxmodsqr(X.x)
		+cmplxmodsqr(X.y)+cmplxmodsqr(X.z);

	return NORM;

}
cmplxvect realvectmult(double C,cmplxvect X) {

	cmplxvect Y;

	Y.x=realmult(C,X.x);
	Y.y=realmult(C,X.y);
	Y.z=realmult(C,X.z);

	return Y;
}

cmplxvect cmplxvectcross(cmplxvect R1,cmplxvect R2) {

	cmplxvect compcross;

	compcross.x=cmplxadd(cmplxmult(R1.y,R2.z),
			realmult(-1.0,cmplxmult(R1.z,R2.y)));
	compcross.y=realmult(-1.0,cmplxadd(cmplxmult(R1.x,R2.z),
				realmult(-1.0,cmplxmult(R1.z,R2.x))));
	compcross.z=cmplxadd(cmplxmult(R1.x,R2.y),
			realmult(-1.0,cmplxmult(R1.y,R2.x)));

	return compcross;
}

cmplxvect unitcmplxvect(cmplxvect R) {

	cmplxvect u={{0.0,0.0},{0.0,0.0},{0.0,0.0}};
	double X=cmplxvectnorm(R);
	if (X>=gtol) {
		u=realvectmult(1.0/X,R);
	}
	return u;
}

cmplxvect re2cmplxvect(vect X) {

	cmplxvect Y={{0.0,0.0},{0.0,0.0},{0.0,0.0}};

	Y.x.re=X.x;
	Y.y.re=X.y;
	Y.z.re=X.z;

	return Y;

}
