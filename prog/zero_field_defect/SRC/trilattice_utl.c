#include "trilattice_sw.h"


void set_boundary(double XQMIN, double YQMIN,double XQMAX,
		double YQMAX, double WMIN, double WMAX,
		int FGRID, int XMGRID, int YMGRID,
		int L,
		int *XMBASE,int *YMBASE,int *WRANGE, int *XMRANGE,
		int *YMRANGE, int *dx, int *dy) {


	// Setting up the frequency grid

	if (FGRID==1 || WMIN==WMAX) {

		(*WRANGE)=1;

	}

	else {

		(*WRANGE)=FGRID;
	}

	// Setting up the momentum range

	if (XQMIN>=0.0 && XQMIN <=1.0) {

		(*XMBASE)=(int)((XQMIN)*L);

	}

	else {
		printf("Correct XQMIN Input");
		exit(0);
	}

	if (YQMIN>=0.0 && YQMIN <=1.0) {

		(*YMBASE)=(int)((YQMIN)*L);

	}

	else {
		printf("Correct YQMIN Input");
		exit(0);
	}


	if (XQMAX>=0.0 && XQMAX <=1.0) {

		(*XMRANGE)=(int)((XQMAX)*L);
	}

	else {
		printf("Correct XQMAX Input");
		exit(0);
	}

	if (YQMAX>=0.0 && YQMAX<=1.0) {

		(*YMRANGE)=(int)((YQMAX)*L);
	}

	else {
		printf("Correct YQMAX Input");
		exit(0);
	}

	// Coarse graining of the grid


	if (((*XMRANGE)-(*XMBASE))>XMGRID) {

		(*dx)=((*XMRANGE)-(*XMBASE))/XMGRID;

	}
	else {

		(*dx)=1;

	}


	if (((*YMRANGE)-(*YMBASE))>YMGRID) {

		(*dy)=((*YMRANGE)-(*YMBASE))/YMGRID;

	}
	else {

		(*dy)=1;

	}


}

void calc_strfact(int STRSWITCH,double GW, lattpoint *lsite,
		double *dthph,complex *T,double *OMEGA,
		double WMIN, double WMAX, double DELTA,
		double ***IMCHI,double ***IMCHI_var, 
		int HCUT, int VCUT,
		int L,int WRANGE,int FGRID,
		int XMBASE, int YMBASE,
		int XMRANGE, int YMRANGE, int dx, int dy,
		int RANGEMIN, int RANGEMAX, int drmx) {

	int i,j,k,l,m,n;

	int N_lattpoint=L*L;

	cmplxvect *SMAG;
	SMAG=(cmplxvect *)calloc(N_lattpoint,sizeof(cmplxvect));

	double ***CHIQ;
	CHIQ=(double ***)malloc(N_lattpoint*sizeof(double **));
	double **SZERO;
	SZERO=(double **)malloc(L*sizeof(double *));

	for (i=0;i<L;i++) {

		SZERO[i]=(double *)calloc(L,sizeof(double));

	}

	for (i=0;i<N_lattpoint;i++) {

		CHIQ[i]=(double **)malloc(L*sizeof(double *));

		for (j=0;j<L;j++) {

			CHIQ[i][j]=(double *)calloc(L,sizeof(double));

		}

	}

	cmplxvect MAGZ;
	double MAGQ;
	double womega;
	double gweight;

	complex QZERO={0.0,0.0};

	if (XMRANGE < L) {
		XMRANGE++;
	}

	if (YMRANGE < L) {
		YMRANGE++;
	}

	double TEMPIMCHI=0.0;

	if(STRSWITCH==1) {

		CHI_MAT(N_lattpoint,lsite,dthph,T,OMEGA,CHIQ,SZERO);

#pragma omp parallel private(i,j,k,l,MAGQ,gweight,womega,TEMPIMCHI)
		{
#pragma omp for
			for (j=YMBASE;j<YMRANGE;j++) {

				for (k=XMBASE;k<XMRANGE;k++) {

					for (i=0;i<WRANGE;i++) {

						MAGQ=0.0;

						womega=WMIN+((WMAX-WMIN)/FGRID)*i;

						for (l=0;l<N_lattpoint;l++) {

							//Note that the excited spectra is 2*OMEGA[l]
							gweight=exp(-pow(womega-fabs(2.0*OMEGA[l]),2.0)/(GW*GW))/(sqrt(2.0*pi)*GW);
							MAGQ+=gweight*(CHIQ[l][j][k]);

						}

						TEMPIMCHI=MAGQ;

						TEMPIMCHI+=exp(-pow(womega,2.0)/(GW*GW))/(sqrt(2.0*pi)*GW)
							*(SZERO[j][k]);
						IMCHI[i][j][k]+=TEMPIMCHI;
						IMCHI_var[i][j][k]+=TEMPIMCHI*TEMPIMCHI;


					}
				}
			}
		}
	}

	else if(STRSWITCH==2) {

		if (STRSWITCH ==2 && RANGEMAX*HCUT < L && RANGEMAX*VCUT< L) {
			RANGEMAX++;
		}


		CHI_MAT(N_lattpoint,lsite,dthph,T,OMEGA,CHIQ,SZERO);

#pragma omp parallel private(i,k,l,MAGQ,gweight,womega,TEMPIMCHI)
		{
#pragma omp for
			for (k=RANGEMIN;k<RANGEMAX;k++) {

				for (i=0;i<WRANGE;i++) {

					MAGQ=0.0;

					womega=WMIN+((WMAX-WMIN)/FGRID)*i;

					for (l=0;l<N_lattpoint;l++) {

						//Note that the excited spectra is 2*OMEGA[l]
						gweight=exp(-pow(womega-fabs(2.0*OMEGA[l]),2.0)/(GW*GW))/(sqrt(2.0*pi)*GW);
						MAGQ+=gweight*(CHIQ[l][k*VCUT][k*HCUT]);

					}

					TEMPIMCHI=MAGQ;

					TEMPIMCHI+=exp(-pow(womega,2.0)/(GW*GW))/(sqrt(2.0*pi)*GW)
						*(SZERO[k*VCUT][k*HCUT]);
					IMCHI[i][k][0]+=TEMPIMCHI;
					IMCHI_var[i][k][0]+=TEMPIMCHI*TEMPIMCHI;


				}
			}
		}
	}

	for (i=0;i<L;i++) {

		free(SZERO[i]);

	}
	free(SZERO);

	for (i=0;i<N_lattpoint;i++) {


		for (j=0;j<L;j++) {

			free(CHIQ[i][j]);

		}

		free(CHIQ[i]);
	}

	free(CHIQ);

	free(SMAG);
}


void print_strfact(FILE *outgrid, int STRSWITCH,
		double WMIN, double WMAX, 
		double DELTA, double ALPHA,
		double ***IMCHI,double ***IMCHI_var, 
		int HCUT, int VCUT,
		int L,int WRANGE,int FGRID,
		int XMBASE, int YMBASE,
		int XMRANGE, int YMRANGE, int dx, int dy,
		int RANGEMIN,int RANGEMAX, int drmx,
		int iter_DISD) {

	int i,j,k,l,m,n;
	vect b1={2.0*pi,-2.0*pi/sqrt(3.0),0.0};
	vect b2={0.0,4.0*pi/(sqrt(3.0)),0.0};

	vect FLOATVEC={0.0,0.0,0.0};

	// Momentum is printed in the units of (1/lattice constant)

	if (XMRANGE < L) {
		XMRANGE++;
	}

	if (YMRANGE < L) {
		YMRANGE++;
	}

	if(STRSWITCH==1) {

		for (i=0;i<WRANGE;i++) {

			for (j=YMBASE;j<YMRANGE;j=j+drmx) {

				for (k=XMBASE;k<XMRANGE;k=k+drmx) {

					FLOATVEC=vectsum(vectmult((1.0*k)/L,b1),
							vectmult((1.0*j)/L,b2));
					fprintf(outgrid,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n"
							,L,DELTA,WMIN+((WMAX-WMIN)/FGRID)*i,FLOATVEC.x,
							FLOATVEC.y,IMCHI[i][j][k]/iter_DISD,
							sqrt(fabs(IMCHI_var[i][j][k]/iter_DISD
							-pow(IMCHI[i][j][k]/iter_DISD,2.0)))/sqrt(iter_DISD),
							iter_DISD);
				}
			}

		}
	}

	else if(STRSWITCH==2) {

		if (RANGEMAX*HCUT < L && RANGEMAX*VCUT< L) {

			RANGEMAX++;

		}

		for (i=0;i<WRANGE;i++) {

			for (k=RANGEMIN;k<RANGEMAX;k=k+drmx) {


				FLOATVEC=vectsum(vectmult((1.0*k*HCUT)/L,b1),
						vectmult((1.0*k*VCUT)/L,b2));
				fprintf(outgrid,"%d\t%lf\t%lf\t%lf\t%lf\t%.12f\t%.12f\t%lf\t%d\n"
						,L,DELTA,WMIN+((WMAX-WMIN)/FGRID)*i,FLOATVEC.x,
						FLOATVEC.y,IMCHI[i][k][0]/iter_DISD,
						sqrt(fabs(IMCHI_var[i][k][0]/iter_DISD
						-pow(IMCHI[i][k][0]/iter_DISD,2.0)))/sqrt(iter_DISD),
						ALPHA,
						iter_DISD);


			}
		}

	}



}


void copylattice(int N_lattpoint,lattpoint *lat_SRC,lattpoint *lat_DEST) {


	int i,j;
#pragma omp parallel private(i,j)
	{
#pragma omp for
	for (i=0;i<N_lattpoint;i++) {

		for (j=0;j<ZCO;j++) {
			lat_DEST[i].neighbour[j]=lat_SRC[i].neighbour[j];
			lat_DEST[i].nn_nbr[j]=lat_SRC[i].nn_nbr[j];
		}
		lat_DEST[i].pos=lat_SRC[i].pos;
		lat_DEST[i].ang=lat_SRC[i].ang;

	}
	}
}



void copycoupling(int N_lattpoint,double **J_SRC,double **J_DEST) {


	int i,j;
#pragma omp parallel private(i,j)
	{
#pragma omp for
	for (i=0;i<N_lattpoint;i++) {
		for (j=0;j<N_lattpoint;j++) {
			J_DEST[i][j]=J_SRC[i][j];
		}
	}
	}
}


