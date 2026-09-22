#include "trilattice_sw.h"

int main(int argc,char *argv[]) {


	// Argument variables
	int L=atoi(argv[1]); // System size
	int iter_DISD=atoi(argv[2]); // Disorder iterations
	double DELTA=atof(argv[3]); // Disorder width
	double ALPHA=atof(argv[4]); // NNN coupling strength
	double BTOL=atof(argv[5]); // Bogoliubov tolerance
	int NMX=atoi(argv[6]); // Run index
	double XQMIN=atof(argv[7]); // Minimum X momentum
	double XQMAX=atof(argv[8]); // Maximum X momentum
	double YQMIN=atof(argv[9]); // Minimum Y momentum
	double YQMAX=atof(argv[10]); // Maximum Y momentum
	int VCUT=atoi(argv[11]); // Vertical cut ratio
	int HCUT=atoi(argv[12]); // Horizontal cut ratio
	int PAR_X=atoi(argv[13]); // Parallel classical samples

	struct timeval begin,end;
	long secs_used, micros_used;

	// Other input parameters

	FILE *fyl;
	fyl=openfile("../constants.txt","r");

	int iter_NUMG=(int)readdouble(fyl,"CITER"); // Classical iterations
	int LSTSWITCH=readint(fyl,"LSTX"); // List observable switch
	int STRSWITCH=readint(fyl,"STRX"); // Structure factor switch
	double GW=readdouble(fyl,"GWIDTH"); // Delta-function width
	double ANGDISPAR=readdouble(fyl,"ANGDWIDTH")*pi;
	int EIGENCOUNT=readint(fyl,"ECOUNT");

	// Momentum and frequency of the structure-factor scaling
	vect Q_F={atof(argv[14])*pi,atof(argv[15])*pi,0.0};
	Q_F.x=Q_F.x/latconst;
	Q_F.y=Q_F.y/latconst;
	double W_F=readdouble(fyl,"WW");

	// Frequency and momentum ranges

	double IPRMIN=readdouble(fyl,"IMIN");
	double IPRMAX=readdouble(fyl,"IMAX");
	double WMIN=readdouble(fyl,"MINF");
	double WMAX=readdouble(fyl,"MAXF");

	// Momentum and frequency grids

	int XMGRID=readint(fyl,"XMOMNUM");
	int YMGRID=readint(fyl,"YMOMNUM");
	int FGRID=readint(fyl,"FREQNUM");
	int IPRGRID=readint(fyl,"IGRID");

	fclose(fyl);

	// Start momenta

	int XMBASE;
	int YMBASE;

	// Frequency and momentum loop ranges

	int XMRANGE;
	int YMRANGE;
	int WRANGE;

	int dx, dy;

	set_boundary(XQMIN,YQMIN,XQMAX,
			YQMAX,WMIN,WMAX,FGRID,XMGRID,YMGRID,
			L,
			&XMBASE,&YMBASE,&WRANGE,&XMRANGE,
			&YMRANGE,&dx,&dy); 

	// BZ cut parameters

	int nx, ny;
	int RANGEMIN, RANGEMAX, drmx;

	if (STRSWITCH==2) {

		if (VCUT==0 && HCUT!=0) {

			RANGEMIN=XMBASE/HCUT;
			RANGEMAX=XMRANGE/HCUT;
			drmx=dx;

		}
		else if (HCUT==0 && VCUT!=0) {

			RANGEMIN=YMBASE/VCUT;
			RANGEMAX=YMRANGE/VCUT;
			drmx=dy;

		}

		else {

			nx=(int)(XMRANGE/HCUT);
			ny=(int)(YMRANGE/VCUT);
			RANGEMAX= (nx >= ny)? ny : nx;

			nx=XMBASE/HCUT;
			ny=YMBASE/VCUT;
			RANGEMIN= (nx >= ny)? nx : ny;

			drmx= (dy >= dx)? dx : dy;

		}

	}

	// Output

	// List output
	char ofname[200];
	snprintf(ofname,200*sizeof(char),"../DATA/RAW/DATA_L-%d_ORD-%d.txt"
			,L,NMX);

	FILE *outf;
	outf=fopen(ofname,"w+");

	// Grid output
	char gridfname[200];
	snprintf(gridfname,200*sizeof(char),"../DATA/RAW/GRID_L-%d_ORD-%d.txt"
			,L,NMX);

	FILE *outgrid;

	outgrid=fopen(gridfname,"w+");

	// Eigenfunction
	char eigenfname[200];
	snprintf(eigenfname,200*sizeof(char),"../DATA/RAW/EIGEN_L-%d_ORD-%d.txt"
			,L,NMX);

	FILE *eigengrid;
	eigengrid=fopen(eigenfname,"w+");

	// Inverse Participation Ratio
	char iprfname[200];
	snprintf(iprfname,200*sizeof(char),"../DATA/RAW/IPR_L-%d_ORD-%d.txt"
			,L,NMX);

	FILE *iprlist;
	iprlist=fopen(iprfname,"w+");

	// Bond-angles
	char anglefname[200];
	snprintf(anglefname,200*sizeof(char),"../DATA/RAW/ANGLE_L-%d_ORD-%d.txt"
			,L,NMX);

	FILE *anglelist;
	anglelist=fopen(anglefname,"w+");

	// Lattice size

	int N_lattpoint=L*L;

	//Number of couplings

	int N_j=N_lattpoint*ZCO;

	int i,j,k,l,diter;
	int p,q;


	// Lattice points

	lattpoint *lsite;
	lsite=(lattpoint *)malloc(N_lattpoint*sizeof(lattpoint));

	for (i=0;i<N_lattpoint;i++) {

		lsite[i].neighbour=(int *)malloc(ZCO*sizeof(int));
		lsite[i].nn_nbr=(int *)malloc(ZCO*sizeof(int));

	}

	// Temporary lattice sites, parallel comparison

	lattpoint **lsite_tmp;
	lsite_tmp=(lattpoint **)malloc(PAR_X*sizeof(lattpoint *));

	for (i=0;i<PAR_X;i++) {
		lsite_tmp[i]=(lattpoint *)malloc(N_lattpoint*sizeof(lattpoint)); 
		for (j=0;j<N_lattpoint;j++) {

			lsite_tmp[i][j].neighbour=(int *)malloc(ZCO*sizeof(int));
			lsite_tmp[i][j].nn_nbr=(int *)malloc(ZCO*sizeof(int));


		}
	}

	double *ENX;
	ENX=(double *)malloc(PAR_X*sizeof(double));
	double *ENX_VAR;
	ENX_VAR=(double *)malloc(PAR_X*sizeof(double));
	double ENX_MIN;

	// Sublattice arrays

	int **sublattice;
	sublattice=(int **)malloc((ZCO/2)*sizeof(int *));

	for (i=0;i<(ZCO/2);i++) {
		sublattice[i]=(int *)malloc(N_lattpoint*sizeof(int));
	}

	// Elementary upward triangles

	int **eltriangle;
	eltriangle=(int **)malloc((L-1)*(L-1)*sizeof(int *));

	for (i=0;i<(L-1)*(L-1);i++) {
		eltriangle[i]=(int *)malloc((ZCO/2)*sizeof(int));
	}


	// Random solid-angle arrays

	double **thetalist;
	double **philist;
	thetalist=(double **)malloc(PAR_X*sizeof(double *));
	philist=(double **)malloc(PAR_X*sizeof(double *));
	for (i=0;i<PAR_X;i++) {
		thetalist[i]=(double *)malloc(N_lattpoint*sizeof(double));
		philist[i]=(double *)malloc(N_lattpoint*sizeof(double));
	}

	// Couplings

	double **J;
	J=(double **)malloc(N_lattpoint*sizeof(double *));
	for (i=0;i<N_lattpoint;i++) {
		J[i]=(double *)malloc(N_lattpoint*sizeof(double));
	}


	// Bogoliubov matrices

	complex *M;
	complex *T;
	double *OMEGA;

	int MDIM;
	MDIM=2*N_lattpoint;


	M=(complex *)malloc(MDIM*MDIM*sizeof(complex));
	T=(complex *)malloc(MDIM*MDIM*sizeof(complex));
	OMEGA=(double *)malloc(MDIM*sizeof(double));

	// Dynamical structure-factor matrix

	double ***IMCHI;
	double ***IMCHI_var;

	if (STRSWITCH==1) {

		IMCHI=(double ***)malloc(WRANGE*sizeof(double **));
		IMCHI_var=(double ***)malloc(WRANGE*sizeof(double **));

		for (i=0;i<WRANGE;i++) {

			IMCHI[i]=(double **)malloc((YMRANGE+1)*sizeof(double *));
			IMCHI_var[i]=(double **)malloc((YMRANGE+1)*sizeof(double *));

			for (j=0;j<=YMRANGE;j++) {

				IMCHI[i][j]=(double *)calloc((XMRANGE+1),sizeof(double));
				IMCHI_var[i][j]=(double *)calloc((XMRANGE+1),sizeof(double));

			}
		}

	}

	else if (STRSWITCH ==2) {

		IMCHI=(double ***)malloc(WRANGE*sizeof(double **));
		IMCHI_var=(double ***)malloc(WRANGE*sizeof(double **));

		for (i=0;i<WRANGE;i++) {

			IMCHI[i]=(double **)malloc((RANGEMAX+1)*sizeof(double *));
			IMCHI_var[i]=(double **)malloc((RANGEMAX+1)*sizeof(double *));

			for (j=0;j<=RANGEMAX;j++) {

				IMCHI[i][j]=(double *)calloc((1),sizeof(double));
				IMCHI_var[i][j]=(double *)calloc((1),sizeof(double));

			}
		}
	}

	// IPR parameters
	double *IPR;
	double *IPR_var;

	if (LSTSWITCH != 0) {

		IPR=(double *)calloc(IPRGRID,sizeof(double));
		IPR_var=(double *)calloc(IPRGRID,sizeof(double));

	}

	// Bond-angle parameters
	double *cangle;
	double *cangle_var;
	double *qangle;
	double *qangle_var;

	if (LSTSWITCH != 0) {

		cangle=(double *)calloc((N_j/2),sizeof(double));
		cangle_var=(double *)calloc((N_j/2),sizeof(double));
		qangle=(double *)calloc((N_j/2),sizeof(double));
		qangle_var=(double *)calloc((N_j/2),sizeof(double));

	}

	// Spin value

	double TEMPSPIN=0.0;
	double SPIN=0.0;
	double SPIN_var=0.0;

	// Ground-state energy

	double EN=0.0;
	double ENTEMP=0.0;
	double EN_var=0.0;

	// Classical ground-state energy

	double cl_EN=0.0;
	double cl_ENTEMP=0.0;
	double cl_EN_var=0.0;

	// Sublattice magnetization

	double SBL_MAG=0.0;
	double SBL_MAGTEMP=0.0;
	double SBL_MAGTEMP1=0.0;
	double SBL_MAGTEMP2=0.0;
	double SBL_MAG_var=0.0;
	double SBL_CHI=0.0;
	double SBL_CHITEMP=0.0;
	vect SBL_CHITEMP1;
	cmplxvect SBL_CHITEMP2;
	double SBL_CHI_var=0.0;

	double tri_SBL_MAG=0.0;
	double tri_SBL_MAGTEMP=0.0;
	double tri_SBL_MAG_var=0.0;
	double tri_SBL_CHI=0.0;
	double tri_SBL_CHITEMP=0.0;
	cmplxvect tri_SBL_CHITEMP1;
	double tri_SBL_CHI_var=0.0;
	double tri_SCALAR_CHI=0.0;
	double tri_SCALAR_CHITEMP=0.0;
	double tri_SCALAR_CHI_var=0.0;

	// Classical sublattice magnetization

	double cl_SBL_MAG=0.0;
	double cl_SBL_MAGTEMP=0.0;
	double cl_SBL_MAG_var=0.0;
	double cl_SBL_CHI=0.0;
	double cl_SBL_CHITEMP=0.0;
	double cl_SBL_CHI_var=0.0;

	double tri_cl_SBL_MAG=0.0;
	double tri_cl_SBL_MAGTEMP=0.0;
	double tri_cl_SBL_MAG_var=0.0;
	double tri_cl_SBL_CHI=0.0;
	double tri_cl_SBL_CHITEMP=0.0;
	vect tri_cl_SBL_CHITEMP1;
	double tri_cl_SBL_CHI_var=0.0;
	double tri_cl_SCALAR_CHI=0.0;
	double tri_cl_SCALAR_CHITEMP=0.0;
	double tri_cl_SCALAR_CHI_var=0.0;

	// Susceptibility scaling

	double SUSC=0.0;
	double SUSC_VAR=0.0;
	double CSUSC=0.0;
	double CSUSC_VAR=0.0;

	// Eigenvector map

	double *eigen_grid;
	eigen_grid=(double *)calloc(N_lattpoint,sizeof(double));
	double *eigen_grid_var;
	eigen_grid_var=(double *)calloc(N_lattpoint,sizeof(double));

	// Local spin map

	double *spin_grid;
	spin_grid=(double *)calloc(N_lattpoint,sizeof(double));
	double *spin_grid_var;
	spin_grid_var=(double *)calloc(N_lattpoint,sizeof(double));


	// Mean-field angle fluctuation

	complex **alpha;
	complex **beta;
	alpha=(complex **)malloc(N_lattpoint*sizeof(complex *));
	beta=(complex **)malloc(N_lattpoint*sizeof(complex *));
	for (i=0;i<N_lattpoint;i++) {
		alpha[i]=(complex *)malloc(N_lattpoint*sizeof(complex));
		beta[i]=(complex *)malloc(N_lattpoint*sizeof(complex));
	}

	double *dthph;
	dthph=(double *)malloc(MDIM*sizeof(double));

	// Disorder iterations

	gettimeofday(&begin,NULL);

	// Lattice blocks

	create_auxlattice(N_lattpoint,sublattice,eltriangle);

	for (i=0;i<PAR_X;i++) {

		create_trilattice(N_lattpoint,lsite_tmp[i]);

	}
	// Classical uniform moment

	double CMOM_X=0.0;
	double CMOM_Y=0.0;
	double CMOM_Z=0.0;
	double CMOM_var=0.0;

    // Quantum uniform moment

	double QMOM_X=0.0;
	double QMOM_Y=0.0;
	double QMOM_Z=0.0;
	double QMOM_var=0.0;

	for (diter=0;diter<iter_DISD;diter++) {

		// Random couplings

		gen_couplings(N_lattpoint,lsite_tmp[0],
				ALPHA,DELTA,J); // Common random Hamiltonian


		for (i=0;i<PAR_X;i++) {

			// Random solid-angle array

			randgen_ang(N_lattpoint,NMX,(i+diter*PAR_X),
					ANGDISPAR,
					thetalist[i],philist[i]);


		}

#pragma omp parallel private(i)
		{
#pragma omp for
			for (i=0;i<PAR_X;i++) {

				ENX[i]=0.0;
				ENX_VAR[i]=0.0;

				// Random spin orientation

				gen_angles(N_lattpoint,thetalist[i],
						philist[i],lsite_tmp[i]);

				// Classical ground state

				classic_algo(diter,i,lsite_tmp[i],J,
						N_lattpoint,iter_NUMG);	

				// Energy

				cl_ground_en(N_lattpoint,lsite_tmp[i],
						J,&ENX[i],&ENX_VAR[i]);


			}
		}

		// Lowest-energy configuration

		j=0;
		ENX_MIN=0.0;
		for (i=0;i<PAR_X;i++) {
			if (ENX[i]<=ENX_MIN) {
				j=i;
				ENX_MIN=ENX[i];
			}

		}
        
		copylattice(N_lattpoint,lsite_tmp[j],lsite);

		// Spin-wave Hamiltonian matrix

		gen_ham(MDIM,lsite,J,M);

		// Bogoliubov transformation

		nubog(MDIM,BTOL,M,T,OMEGA);	

		// Mean-field parameters

		mf_params(N_lattpoint,lsite,T,alpha,beta);

		// Angle corrections

		mf_angles(N_lattpoint,lsite,J,alpha,beta,BTOL,dthph);

		// List variables

		if (LSTSWITCH !=0) {

			// Magnetization

			stag_spin(N_lattpoint,lsite,alpha,
					OMEGA,&SPIN,&SPIN_var,
					spin_grid,spin_grid_var);

            // Energy

			ground_en(N_lattpoint,lsite,J,
					OMEGA,&EN,&EN_var);

			// Classical energy

			cl_ground_en(N_lattpoint,lsite,
					J,&cl_EN,&cl_EN_var);

			// Structure-factor scaling

			im_chi(W_F,Q_F,GW,
					N_lattpoint,lsite,
					dthph,T,OMEGA,
					&SUSC,&SUSC_VAR,
					&CSUSC,&CSUSC_VAR);

			// Sublattice magnetization

			qsubmag(N_lattpoint,lsite,
					sublattice,dthph,
					T,OMEGA,alpha,&SBL_MAG,
					&SBL_CHI,&SBL_MAG_var,
					&SBL_CHI_var);

			cl_sublatmag(N_lattpoint,lsite,
					sublattice,
					&cl_SBL_MAG,
					&cl_SBL_CHI,
					&cl_SBL_MAG_var,
					&cl_SBL_CHI_var);

			// Elementary-triangle magnetization

			tri_cl_sublatmag(N_lattpoint,
					lsite,eltriangle,
					&tri_cl_SBL_MAG,
					&tri_cl_SCALAR_CHI,
					&tri_cl_SBL_CHI,
					&tri_cl_SBL_MAG_var,
					&tri_cl_SCALAR_CHI_var,
					&tri_cl_SBL_CHI_var);

			tri_qsubmag_ex(N_lattpoint,
					lsite,eltriangle,
					dthph,T,OMEGA,
					alpha,
					&tri_SBL_MAG,
					&tri_SCALAR_CHI,
					&tri_SBL_CHI,
					&tri_SBL_MAG_var,
					&tri_SCALAR_CHI_var,
					&tri_SBL_CHI_var);

			// Eigenvectors on the lattice

			eigenview(EIGENCOUNT,N_lattpoint,T,
					OMEGA,eigen_grid,eigen_grid_var);

			// IPR

			IPR_CALC(IPRMIN,IPRMAX,IPRGRID,GW,
					BTOL,N_lattpoint,T,OMEGA,IPR,IPR_var);

			// Bond angles

			histo_angles(N_lattpoint,lsite,dthph,
					cangle,cangle_var,qangle,qangle_var);

			cl_moment(N_lattpoint,lsite,
					&CMOM_X,&CMOM_Y,&CMOM_Z,
					&CMOM_var);

			q_moment(N_lattpoint,lsite,alpha,
					OMEGA,dthph,
					&QMOM_X,&QMOM_Y,&QMOM_Z,
					&QMOM_var);

		}

		// Structure factor

		if (STRSWITCH != 0) {

			calc_strfact(STRSWITCH,GW,lsite,
					dthph,T,OMEGA,
					WMIN, WMAX, DELTA,
					IMCHI,IMCHI_var, 
					HCUT,VCUT,
					L,WRANGE,FGRID,
					XMBASE,YMBASE,
					XMRANGE,YMRANGE,dx,dy,
					RANGEMIN,RANGEMAX,drmx);

		}

	}

	gettimeofday(&end,NULL);
	secs_used =(end.tv_sec-begin.tv_sec);
	micros_used =((secs_used*1000000)+end.tv_usec)-begin.tv_usec;
    printf("%d iterations for L = %d done in %lf seconds\n",iter_DISD,L,(double)micros_used/pow(10.0,6.0));

	int LAT_SMPLE=N_lattpoint*iter_DISD;
	int TRI_SMPLE=(L-1)*(L-1)*iter_DISD;

	if (LSTSWITCH != 0) {
		fprintf(outf,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%.12f\t%.12f\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%.12f\t%.12f\t%d"
				"\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\n",
				L, // 1
				DELTA, // 2
				SYSAVG(EN,iter_DISD), // 3
				SYSERR(EN,EN_var,iter_DISD), // 4
				SYSAVG(SPIN,LAT_SMPLE), // 5
				SYSERR(SPIN,SPIN_var,LAT_SMPLE), // 6
				SYSAVG(SBL_MAG,iter_DISD), // 7
				SYSERR(SBL_MAG,SBL_MAG_var,iter_DISD), // 8
				SYSAVG(SBL_CHI,iter_DISD), // 9
				SYSERR(SBL_CHI,SBL_CHI_var,iter_DISD), // 10
				SYSAVG(cl_SBL_MAG,iter_DISD), // 11
				SYSERR(cl_SBL_MAG,cl_SBL_MAG_var,iter_DISD), // 12
				SYSAVG(cl_SBL_CHI,iter_DISD), // 13
				SYSERR(cl_SBL_CHI,cl_SBL_CHI_var,iter_DISD), // 14
				SYSAVG(tri_SBL_MAG,TRI_SMPLE), // 15
				SYSERR(tri_SBL_MAG,tri_SBL_MAG_var,TRI_SMPLE), // 16
				SYSAVG(tri_SBL_CHI,iter_DISD), // 17
				SYSERR(tri_SBL_CHI,tri_SBL_CHI_var,iter_DISD), // 18
				SYSAVG(tri_cl_SBL_MAG,TRI_SMPLE), // 19
				SYSERR(tri_cl_SBL_MAG,tri_cl_SBL_MAG_var,TRI_SMPLE), // 20
				SYSAVG(tri_cl_SBL_CHI,iter_DISD), // 21
				SYSERR(tri_cl_SBL_CHI,tri_cl_SBL_CHI_var,iter_DISD), // 22
				SYSAVG(cl_EN,iter_DISD), // 23
				SYSERR(cl_EN,cl_EN_var,iter_DISD), // 24
				SYSAVG(SUSC,iter_DISD), // 25
				SYSERR(SUSC,SUSC_VAR,iter_DISD), // 26
				ALPHA, // 27
				SYSAVG(tri_SCALAR_CHI,TRI_SMPLE), // 28
				SYSERR(tri_SCALAR_CHI,tri_SCALAR_CHI_var,TRI_SMPLE), // 29
				SYSAVG(tri_cl_SCALAR_CHI,TRI_SMPLE), // 30
				SYSERR(tri_cl_SCALAR_CHI,tri_cl_SCALAR_CHI_var,TRI_SMPLE), // 31
				SYSAVG(CSUSC,iter_DISD), // 32
				SYSERR(CSUSC,CSUSC_VAR,iter_DISD), // 33
				iter_DISD, //34
				SYSAVG(CMOM_X,iter_DISD), // 35
				SYSAVG(CMOM_Y,iter_DISD), // 36
				SYSAVG(CMOM_Z,iter_DISD), // 37
				SYSAVG(QMOM_X,iter_DISD), // 38
				SYSAVG(QMOM_Y,iter_DISD), // 39
				SYSAVG(QMOM_Z,iter_DISD) // 40
				);

		for (i=0;i<IPRGRID;i++) {

			fprintf(iprlist,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
					L,DELTA,
					IPRMIN+((IPRMAX-IPRMIN)*i)/IPRGRID,
					SYSAVG(IPR[i],iter_DISD),
					SYSERR(IPR[i],IPR_var[i],iter_DISD),
					ALPHA,
					iter_DISD);

		}

		for (i=0;i<(N_j/2);i++) {

			fprintf(anglelist,"%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%d\n",
					L,
					DELTA,
					ALPHA,
					i,
					SYSAVG(cangle[i],iter_DISD),
					SYSERR(cangle[i],cangle_var[i],iter_DISD),
					SYSAVG(qangle[i],iter_DISD),
					SYSERR(qangle[i],qangle_var[i],iter_DISD),
					iter_DISD
			       );

		}


		for (i=0;i<N_lattpoint;i++) {

			fprintf(eigengrid,"%d\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\n",
					L,DELTA,ALPHA,
					(i%L),(i/L),
					SYSAVG(eigen_grid[i],iter_DISD),
					SYSERR(eigen_grid[i],eigen_grid_var[i],iter_DISD),
					SYSAVG(spin_grid[i],iter_DISD),
					SYSERR(spin_grid[i],spin_grid_var[i],iter_DISD),
					iter_DISD
			       );

		}


	}

	if (LSTSWITCH != 0) {

		free(IPR);
		free(IPR_var);
		free(cangle);
		free(cangle_var);
		free(qangle);
		free(qangle_var);

	}

	fclose(outf);
	fclose(eigengrid);
	fclose(iprlist);
	fclose(anglelist);

	if (LSTSWITCH == 0) {

		remove(ofname);
		remove(eigenfname);
		remove(iprfname);

	}

	if (STRSWITCH != 0) {


		print_strfact(outgrid,STRSWITCH,
				WMIN, WMAX,
				DELTA,ALPHA,
				IMCHI,IMCHI_var, 
				HCUT,VCUT,
				L,WRANGE,FGRID,
				XMBASE,YMBASE,
				XMRANGE,YMRANGE,dx,dy,
				RANGEMIN,RANGEMAX,drmx,
				iter_DISD); 


	}

	fclose(outgrid);

	if(STRSWITCH == 0) {

		remove(gridfname);

	}

	if (STRSWITCH==1) {

		for (i=0;i<WRANGE;i++) {

			for (j=0;j<=YMRANGE;j++)	{

				free(IMCHI[i][j]);
				free(IMCHI_var[i][j]);

			}

			free(IMCHI[i]);
			free(IMCHI_var[i]);



		}


		free(IMCHI);
		free(IMCHI_var);
	}

	else if (STRSWITCH==2) {

		for (i=0;i<WRANGE;i++) {

			for (j=0;j<=RANGEMAX;j++)	{

				free(IMCHI[i][j]);
				free(IMCHI_var[i][j]);

			}

			free(IMCHI[i]);
			free(IMCHI_var[i]);



		}


		free(IMCHI);
		free(IMCHI_var);



	}

	for (i=0;i<N_lattpoint;i++) {
		free(J[i]);
		free(alpha[i]);
		free(beta[i]);
	}
	free(J);
	free(alpha);
	free(beta);

	for (i=0;i<PAR_X;i++) {
		free(thetalist[i]);
		free(philist[i]);
	}
	free(thetalist);
	free(philist);
	free(dthph);
	for (i=0;i<(ZCO/2);i++) {
		free(sublattice[i]);
	}
	free(sublattice);
	for (i=0;i<(L-1)*(L-1);i++) {
		free(eltriangle[i]);
	}
	free(eltriangle);

	for (i=0;i<N_lattpoint;i++) {
		free(lsite[i].neighbour);
		free(lsite[i].nn_nbr);
	}

	for (i=0;i<PAR_X;i++) {

		for (j=0;j<N_lattpoint;j++) {
			free(lsite_tmp[i][j].neighbour);
			free(lsite_tmp[i][j].nn_nbr);
		}
		free(lsite_tmp[i]);
	}

	free(lsite);
	free(lsite_tmp);
	free(ENX);
	free(ENX_VAR);
	free(M);
	free(T);
	free(OMEGA);
	free(eigen_grid);
	free(eigen_grid_var);
	free(spin_grid);
	free(spin_grid_var);

	return 0;
}
