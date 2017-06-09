#include "trilattice_sw.h"


// Random array shuffler

void shuffle(int diter,int PARX,int ORDX,int *array,int n) {    

	int i,j,t;

	struct timeval now;
	gettimeofday(&now,NULL);
	long long tv[4]={now.tv_sec,diter,PARX,ORDX}, tv_leng=4;
	init_by_array64(tv,tv_leng);

	if (n > 1) {

		for (i=0;i<(n-1);i++) {
			j = i + (genrand64_int64() >> 11 ) / (9007199254740991 / (n - i) + 1);
			t = array[j];
			array[j] = array[i];
			array[i] = t;
		}

	}
}


// Performing the classical algorithm

void classic_algo(int diter,int PARX,lattpoint *lsite,
		double **J,int N_lattpoint,int iter_NUMG) {


	int i,j,k,l,p,q;

	// Allocating the spins and the local fields

	vect nullvect={0.0,0.0,0.0};

	vect *spinvect;
	vect *spinvectpr;
	vect localvect;

	spinvect=(vect *)malloc(N_lattpoint*sizeof(vect));
	spinvectpr=(vect *)malloc(N_lattpoint*sizeof(vect));

	// Array of spins to be flipped

	int *FSPIN;
	FSPIN=(int *)malloc(N_lattpoint*sizeof(int));

	for (i=0;i<N_lattpoint;i++) {
		FSPIN[i]=i;
	}

	// Defining tolerances

	double hyper_norm; // This is just a non-zero garbage value
	double single_tol=gtol;
	double tot_tol=single_tol;

	//  The algorithm to find the classical ground state	

	// Defining the spins and the local fields


	for (i=0;i<N_lattpoint;i++) {

		spinvect[i].x=(sin(lsite[i].ang.th)*cos(lsite[i].ang.ph));
		spinvect[i].y=(sin(lsite[i].ang.th)*sin(lsite[i].ang.ph));
		spinvect[i].z=(cos(lsite[i].ang.th));
		spinvectpr[i]=spinvect[i];

	}

	double NRM, NRMTEMP;

	for (i=0;i<iter_NUMG;i++) {
		hyper_norm=0.0;
		NRM=0.0;	
		for (l=0;l<N_lattpoint;l++) {
			j=FSPIN[l];
			localvect=nullvect;
			for (k=0;k<ZCO;k++) {
				localvect=vectsum(localvect,vectmult(J[j][lsite[j].neighbour[k]],
							spinvect[lsite[j].neighbour[k]]));
				localvect=vectsum(localvect,vectmult(J[j][lsite[j].nn_nbr[k]],
							spinvect[lsite[j].nn_nbr[k]]));

			}
			spinvect[j]=vectmult(-1.0,unitvect(localvect));
			NRMTEMP=vectnorm(vectsum(vectmult(-1.0,spinvectpr[j]),spinvect[j]));
			NRMTEMP=sqrt(NRMTEMP);
			if (NRMTEMP > NRM) {
				NRM=NRMTEMP;
			}
			hyper_norm+=NRMTEMP*NRMTEMP;
			spinvectpr[j]=spinvect[j];
		}
		hyper_norm=sqrt(hyper_norm);
		shuffle(diter,PARX,i,FSPIN,N_lattpoint);
		if (hyper_norm <= tot_tol && NRM <= tot_tol) {
			break;
		}
	}

    printf("%d\n",i);

	// Getting the theta's and phi's back;


	for (i=0;i<N_lattpoint;i++) {

		if (spinvect[i].z >= 0) {

			lsite[i].ang.th=acos(spinvect[i].z);

		}

		else if (spinvect[i].z < 0) {

			lsite[i].ang.th=pi-acos(fabs(spinvect[i].z));

		}

		if (fabs(spinvect[i].x)<=gtol && spinvect[i].y >= 0.0) {

			lsite[i].ang.ph=pi/2.0;

		}

		else if (fabs(spinvect[i].x)<=gtol && spinvect[i].y < 0.0) {

			lsite[i].ang.ph=3.0*pi/2.0;

		}

		else if (spinvect[i].x > gtol && spinvect[i].y >=0.0) {

			lsite[i].ang.ph=atan(spinvect[i].y/spinvect[i].x);

		}

		else if (spinvect[i].x < -gtol && spinvect[i].y >=0.0) {

			lsite[i].ang.ph=pi-atan(spinvect[i].y/fabs(spinvect[i].x));

		}

		else if (spinvect[i].x < -gtol && spinvect[i].y < 0.0) {

			lsite[i].ang.ph=pi+atan(spinvect[i].y/spinvect[i].x);
		}

		else if (spinvect[i].x > gtol && spinvect[i].y < 0.0) {

			lsite[i].ang.ph=2.0*pi-atan(fabs(spinvect[i].y)/spinvect[i].x);
		}

	}	

	free(FSPIN);
	free(spinvect);
	free(spinvectpr);

}
