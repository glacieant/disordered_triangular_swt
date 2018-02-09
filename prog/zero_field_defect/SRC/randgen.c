#include "trilattice_sw.h"

void gen_couplings(int N_lattpoint,lattpoint *lsite,double ALPHA,double DELTA,double **J) {

	int i,j,k,l,p,q;

    int L = sqrt(N_lattpoint);

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
				J[i][lsite[i].neighbour[k]]=1.0;
				J[lsite[i].neighbour[k]][i]=J[i][lsite[i].neighbour[k]];
				j++;
			}
			if (J[i][lsite[i].nn_nbr[k]]<-51.5) {
				J[i][lsite[i].nn_nbr[k]]=ALPHA;
				J[lsite[i].nn_nbr[k]][i]=J[i][lsite[i].nn_nbr[k]];
				j++;
			}
		}

	}

    int imp = (L*(L-1))/2-1;
	J[imp][lsite[imp].neighbour[3]]=1.0-DELTA;
	J[lsite[imp].neighbour[3]][imp]=J[imp][lsite[imp].neighbour[3]];
	
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

	int i;

	for (i=0;i<N_lattpoint;i++) {


		lsite[i].ang.th=thetalist[i];
		lsite[i].ang.ph=philist[i];


	}



}

void randgen_ang(int N_lattpoint,int NMX,int par,
		double ANGDISPAR,double *thetalist,double *philist) {

	int i,j;

	int NL=sqrt(N_lattpoint);

	double fleet;
	
	struct timeval now;
	gettimeofday(&now,NULL);
	long long t[3]={par,now.tv_sec,NMX}, t_leng=3;
	init_by_array64(t,t_leng);

	// Generating the random angle array to be assigned to the solid angles

	for (i=0;i<N_lattpoint;i++) {
        
		if (i==0) {
			fleet=0.0;
			philist[i]=fleet-ANGDISPAR+ANGDISPAR*genrand64_real1();
			thetalist[i]=pi/2.0;//-ANGDISPAR+ANGDISPAR*genrand64_real1();
		}

		else if (i!=0 && i%NL==0) {
			fleet=philist[i-1];
			philist[i]=fleet-ANGDISPAR+ANGDISPAR*genrand64_real1();
			thetalist[i]=pi/2.0;//-ANGDISPAR+ANGDISPAR*genrand64_real1();
		}

		else if (i!=0 && i%NL!=0) {
			fleet=philist[i-1]+2.0*pi/3.0;
			if (fleet>=2.0*pi) {
				fleet=fleet-2.0*pi;
			}
			else if (fleet<0) {
				fleet=fleet+2.0*pi;
			}
			philist[i]=fleet-ANGDISPAR+ANGDISPAR*genrand64_real1();
			thetalist[i]=pi/2.0;//-ANGDISPAR+ANGDISPAR*genrand64_real1();
		}

	}

}



void randgen_cpl(int N_lattpoint,int NMX,int diter,
		double delta,double *randlist) {

	int i,j;

	int N_j=N_lattpoint*ZCO;

	double fleet;
	// Generating the random number array to be assigned to the couplings

	struct timeval now;
	gettimeofday(&now,NULL);
	long long t[3]={diter,now.tv_sec,NMX}, t_leng=3;
	init_by_array64(t,t_leng);

	for (i=0;i<N_j;i++) {		
		randlist[i]=(1.0-delta+2.0*delta*genrand64_real1());
	}

}

