#include "trilattice_sw.h"

void gen_couplings(int N_lattpoint,lattpoint *lsite,double ALPHA,double *randlist,double **J) {

	int i,j,k,l,p,q;

	// Assigning the random couplings

	for (i=0;i<N_lattpoint;i++) {
		for (j=0;j<N_lattpoint;j++) {
			J[i][j]=-103.0;
		}
	}

	j=0;

	for (i=0;i<N_lattpoint;i++) {
		for (k=0;k<ZCO;k++) {
			if (J[i][lsite[i].neighbour[k]]<-51.5) {
				J[i][lsite[i].neighbour[k]]=randlist[j];
				J[lsite[i].neighbour[k]][i]=J[i][lsite[i].neighbour[k]];
				j++;
			}
			if (J[i][lsite[i].nn_nbr[k]]<-51.5) {
				J[i][lsite[i].nn_nbr[k]]=ALPHA*randlist[j];
				J[lsite[i].nn_nbr[k]][i]=J[i][lsite[i].nn_nbr[k]];
				j++;
			}

		}

	}

	for (i=0;i<N_lattpoint;i++) {
		for (j=0;j<N_lattpoint;j++) {
			if (J[i][j]<-51.5 || i==j) {
				J[i][j]=0.0;
			}
		}
	}



}


void gen_angles(int N_lattpoint,double *thetalist,double *philist,lattpoint *lsite) {

	//Assigning random orientation of the spin vector

	int i,j,k,l,p,q;

	for (i=0;i<N_lattpoint;i++) {


		lsite[i].ang.th=thetalist[i];
		lsite[i].ang.ph=philist[i];


	}



}

void randgen(int N_lattpoint,int par,double delta,double ANGDISPAR,double *randlist,double *thetalist,double *philist) {

	int i,j;

	int cuplnum=N_lattpoint*ZCO;
	int NL=sqrt(N_lattpoint);

	double fleet;
	// Generating the random number array to be assigned to the couplings

	struct timeval now;
	gettimeofday(&now,NULL);
	long long t[2]={par,now.tv_sec}, t_leng=2;
	init_by_array64(t,t_leng);

//	srand(1342569);

	for (i=0;i<cuplnum;i++) {		
		randlist[i]=(1.0-delta+2.0*delta*((double)rand()/RAND_MAX));
	}

	for (i=0;i<N_lattpoint;i++) {

		philist[i]=2.0*pi*((double)rand()/RAND_MAX);
		thetalist[i]=(pi/2.0)*((double)rand()/RAND_MAX);

	}

	// Generating the random angle array to be assigned to the solid angles
/*
	for (i=0;i<angnum;i++) {

		if (i==0) {
			fleet=0.0;
			philist[i]=fleet-ANGDISPAR+ANGDISPAR*((double)rand()/RAND_MAX);
			thetalist[i]=pi/2.0-ANGDISPAR+ANGDISPAR*((double)rand()/RAND_MAX);
		}

		else if (i!=0 && i%NL==0) {
			fleet=philist[i-1];
			philist[i]=fleet-ANGDISPAR+ANGDISPAR*((double)rand()/RAND_MAX);
			thetalist[i]=pi/2.0-ANGDISPAR+ANGDISPAR*((double)rand()/RAND_MAX);
		}

		else if (i%NL!=0) {
			fleet=philist[i-1]+2.0*pi/3.0;
			if (fleet>=2.0*pi) {
				fleet=fleet-2.0*pi;
			}
			else if (fleet<0) {
				fleet=fleet+2.0*pi;
			}
			philist[i]=fleet-ANGDISPAR+ANGDISPAR*((double)rand()/RAND_MAX);
			thetalist[i]=pi/2.0-ANGDISPAR+ANGDISPAR*((double)rand()/RAND_MAX);
		}

	}
*/

}
