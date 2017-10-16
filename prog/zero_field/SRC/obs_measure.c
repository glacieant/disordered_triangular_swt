#include "trilattice_sw.h"

double SYSAVG(double X,int NSMPLE) {

	return (X/NSMPLE);

}

double SYSERR(double X,double XSQR,int NSMPLE) {

	double S1=(XSQR/NSMPLE);
	double S2=pow(SYSAVG(X,NSMPLE),2.0);
	return sqrt(fabs(S1-S2))*((S1 > S2) ? 1.0 : -1.0);

}


void ground_en(int N_lattpoint,lattpoint *lsite,double **J,
		double *OMEGA,double *EN,double *EN_var) {

	int i,j, k;

	int N_j=N_lattpoint*ZCO/2;

	double H0=0;

	for (i=0;i<N_lattpoint;i++) {

		for (k=0;k<ZCO;k++) {

			j=lsite[i].neighbour[k];
			H0+=0.5*S*(S+1)*Fzz(i,j,J,lsite);

			j=lsite[i].nn_nbr[k];
			H0+=0.5*S*(S+1)*Fzz(i,j,J,lsite);

		}


	}

	for (i=0;i<N_lattpoint;i++) {

		H0+=fabs(OMEGA[i]);
	}

	H0=H0/N_j;
	(*EN)+=H0;
	(*EN_var)+=H0*H0;

}

void cl_ground_en(int N_lattpoint,lattpoint *lsite,double **J,
		double *cl_EN,double *cl_EN_var) {

	int i,j, k;

	int N_j=N_lattpoint*ZCO/2;

	double H0=0;

	for (i=0;i<N_lattpoint;i++) {

		for (k=0;k<ZCO;k++) {

			j=lsite[i].neighbour[k];
			H0+=0.5*S*S*Fzz(i,j,J,lsite);

			j=lsite[i].nn_nbr[k];
			H0+=0.5*S*S*Fzz(i,j,J,lsite);

		}

	}

	H0=H0/N_j;
	(*cl_EN)+=H0;
	(*cl_EN_var)+=H0*H0;

}

void stag_spin(int N_lattpoint, lattpoint *lsite, complex **alpha,
		double *OMEGA,double *SPIN,double *SPIN_var,
		double *spin_grid,double *spin_grid_var) {

	int i,j;

	int DIM=2*N_lattpoint;

	double VAL=0.0;
	double VAL2=0.0;
	double VAL2_VAR=0.0;

	for (i=0;i<N_lattpoint;i++) {

		VAL=S-alpha[i][i].re;
		spin_grid[i]+=VAL;
		spin_grid_var[i]+=VAL*VAL;
		VAL2+=VAL;
		VAL2_VAR+=VAL*VAL;

	}

	VAL2=VAL2;
	VAL2_VAR=VAL2_VAR;
	(*SPIN)+=VAL2;
	(*SPIN_var)+=VAL2_VAR;

}


void cl_sublatmag(int N_lattpoint,lattpoint *lsite,
		int **sublattice,double *CLSBL_ANG,
		double *CLSBL_CHI,double *CLSBL_ANG_VAR,
		double *CLSBL_CHI_VAR) {

	int i,j,k,l,m,n,p,q;

	int DIM=2*N_lattpoint;
	int subnum=N_lattpoint/3;

	double VAL1=0.0;
	double VAL2=0.0;
	double VAL3=0.0;

	vect *spin;
	spin=(vect *)calloc(3,sizeof(vect));
	vect *spintemp;
	spintemp=(vect *)calloc(3,sizeof(vect));


	for (i=0;i<subnum;i++) {

		VAL1=S;
		VAL2=S;
		VAL3=S;

		l=sublattice[0][i];
		m=sublattice[1][i];
		n=sublattice[2][i];

		spintemp[0].x=-sin(lsite[l].ang.th)*cos(lsite[l].ang.ph)*VAL1;
		spintemp[0].y=sin(lsite[l].ang.th)*sin(lsite[l].ang.ph)*VAL1;
		spintemp[0].z=cos(lsite[l].ang.th)*VAL1;

		spin[0]=vectsum(spin[0],spintemp[0]);

		spintemp[1].x=-sin(lsite[m].ang.th)*cos(lsite[m].ang.ph)*VAL2;
		spintemp[1].y=sin(lsite[m].ang.th)*sin(lsite[m].ang.ph)*VAL2;
		spintemp[1].z=cos(lsite[m].ang.th)*VAL2;

		spin[1]=vectsum(spin[1],spintemp[1]);

		spintemp[2].x=-sin(lsite[n].ang.th)*cos(lsite[n].ang.ph)*VAL3;
		spintemp[2].y=sin(lsite[n].ang.th)*sin(lsite[n].ang.ph)*VAL3;
		spintemp[2].z=cos(lsite[n].ang.th)*VAL3;

		spin[2]=vectsum(spin[2],spintemp[2]);

	}

	double VAL=0.0;

	for (p=0;p<3;p++) {

		VAL+=vectdot(spin[p],spin[p]);

	}

	vect CHIR={0.0,0.0,0.0};

	CHIR=vectsum(CHIR,vectcross(spin[0],spin[1]));
	CHIR=vectsum(CHIR,vectcross(spin[1],spin[2]));
	CHIR=vectsum(CHIR,vectcross(spin[2],spin[1]));

	double CHI_VAL=0.0;
	double SBL_VAL=0.0;

	CHI_VAL=(72.0/(sqrt(3.0)))*
		vectnorm(CHIR)/(N_lattpoint*N_lattpoint);
	SBL_VAL=(12.0*VAL)/(N_lattpoint*N_lattpoint);

	(*CLSBL_ANG)+=SBL_VAL;
	(*CLSBL_CHI)+=CHI_VAL;
	(*CLSBL_ANG_VAR)+=SBL_VAL*SBL_VAL;
	(*CLSBL_CHI_VAR)+=CHI_VAL*CHI_VAL;

	free(spin);
	free(spintemp);


}


void qsubmag(int N_lattpoint,lattpoint *lsite,
		int **sublattice,double *dthph,complex *T,double *OMEGA,
		complex **alpha,double *SBLMAG,double *SBLCHI,
		double *SBLMAG_VAR,double *SBLCHI_VAR) {

	int i,j,k,l,m,n,p,nu;
	int ix;

	int DIM=2*N_lattpoint;
	int sblnum=N_lattpoint/3;

	complex Y1={0.0,0.0};
	complex Y2={0.0,0.0};
	complex Y3={0.0,0.0};
	complex VAL2={0.0,0.0};
	complex VAL3={0.0,0.0};

	complex QZERO={0.0,0.0};
	cmplxvect QZEROVECT={{0.0,0.0},{0.0,0.0},{0.0,0.0}};

	cmplxvect *SMAG;
	SMAG=(cmplxvect *)malloc((ZCO/2)*N_lattpoint*sizeof(cmplxvect));
	double SMAGEX=0.0;
	cmplxvect SCHIEX={{0.0,0.0},{0.0,0.0},{0.0,0.0}};
	cmplxvect SCHIEX_TMP={{0.0,0.0},{0.0,0.0},{0.0,0.0}};
	double SCHIEX_x_re=0.0,SCHIEX_x_im=0.0, \
		SCHIEX_y_re=0.0,SCHIEX_y_im=0.0, \
		SCHIEX_z_re=0.0,SCHIEX_z_im=0.0;

#pragma omp parallel private(n,nu,ix,l,Y1,Y2,Y3, \
		VAL3,VAL2,SCHIEX_TMP) reduction(+:SMAGEX,SCHIEX_x_re,SCHIEX_x_im, \
			SCHIEX_y_re,SCHIEX_y_im, \
			SCHIEX_z_re,SCHIEX_z_im)
	{
#pragma omp for
	for (n=0;n<N_lattpoint;n++) {

		for (nu=0;nu<(ZCO/2);nu++) {

			SMAG[nu+n*(ZCO/2)]=QZEROVECT;

			for (ix=0;ix<sblnum;ix++) {

				l=sublattice[nu][ix];

				Y1.re=cos(lsite[l].ang.th)*cos(lsite[l].ang.ph);
				Y1.im=-sin(lsite[l].ang.ph);
				Y2.re=-cos(lsite[l].ang.th)*sin(lsite[l].ang.ph);
				Y2.im=-cos(lsite[l].ang.ph);
				Y3.re=sin(lsite[l].ang.th);
				Y3.im=0.0;

				VAL3=cmplxconj(T[l+n*DIM]);	
				VAL2=cmplxconj(T[l+N_lattpoint+n*DIM]);	

				SMAG[nu+n*(ZCO/2)].x=cmplxadd(SMAG[nu+n*(ZCO/2)].x,realmult(sqrt(S/2.0),
							cmplxadd(cmplxmult(Y1,VAL2),cmplxmult(cmplxconj(Y1),VAL3))));
				SMAG[nu+n*(ZCO/2)].y=cmplxadd(SMAG[nu+n*(ZCO/2)].y,realmult(sqrt(S/2.0),
							cmplxadd(cmplxmult(Y2,VAL2),cmplxmult(cmplxconj(Y2),VAL3))));
				SMAG[nu+n*(ZCO/2)].z=cmplxadd(SMAG[nu+n*(ZCO/2)].z,realmult(sqrt(S/2.0),
							cmplxadd(cmplxmult(Y3,VAL2),cmplxmult(cmplxconj(Y3),VAL3))));

			}

			SMAGEX+=cmplxvectnormsqr(SMAG[nu+n*(ZCO/2)]);

		}

		SCHIEX_TMP=cmplxvectcross(SMAG[0+n*(ZCO/2)],cmplxconjvect(SMAG[1+n*(ZCO/2)]));
		SCHIEX_TMP=cmplxvectsum(SCHIEX_TMP,
				cmplxvectcross(SMAG[1+n*(ZCO/2)],cmplxconjvect(SMAG[2+n*(ZCO/2)])));
		SCHIEX_TMP=cmplxvectsum(SCHIEX_TMP,
				cmplxvectcross(SMAG[2+n*(ZCO/2)],cmplxconjvect(SMAG[0+n*(ZCO/2)])));

		SCHIEX_x_re+=SCHIEX_TMP.x.re;
		SCHIEX_y_re+=SCHIEX_TMP.y.re;
		SCHIEX_z_re+=SCHIEX_TMP.z.re;
		SCHIEX_x_im+=SCHIEX_TMP.x.im;
		SCHIEX_y_im+=SCHIEX_TMP.y.im;
		SCHIEX_z_im+=SCHIEX_TMP.z.im;


	}
	}

	free(SMAG);

	SCHIEX.x.re=SCHIEX_x_re;
	SCHIEX.y.re=SCHIEX_y_re;
	SCHIEX.z.re=SCHIEX_z_re;
	SCHIEX.x.im=SCHIEX_x_im;
	SCHIEX.y.im=SCHIEX_y_im;
	SCHIEX.z.im=SCHIEX_z_im;

	double SVAL1=0.0;
	double SVAL2=0.0;
	double SVAL3=0.0;

	vect *spin;
	spin=(vect *)calloc(3,sizeof(vect));

	double spin0_x=0.0,spin0_y=0.0,spin0_z=0.0,
	       spin1_x=0.0,spin1_y=0.0,spin1_z=0.0,
	       spin2_x=0.0,spin2_y=0.0,spin2_z=0.0;

#pragma omp parallel private(i,SVAL1,SVAL2,SVAL3, \
		l,m,n) reduction(+:spin0_x,spin0_y,spin0_z, \
			spin1_x,spin1_y,spin1_z, \
			spin2_x,spin2_y,spin2_z)
	{
#pragma omp for
	for (i=0;i<sblnum;i++) {

		l=sublattice[0][i];
		m=sublattice[1][i];
		n=sublattice[2][i];

		SVAL1=S-alpha[l][l].re;
		SVAL2=S-alpha[m][m].re;
		SVAL3=S-alpha[n][n].re;

		spin0_x+=-sin(lsite[l].ang.th)*cos(lsite[l].ang.ph)*SVAL1;
		spin0_x+=-(cos(lsite[l].ang.th)*cos(lsite[l].ang.ph)*dthph[l]-
				sin(lsite[l].ang.th)*sin(lsite[l].ang.ph)*dthph[l+N_lattpoint]);

		spin0_y+=sin(lsite[l].ang.th)*sin(lsite[l].ang.ph)*SVAL1;
		spin0_y+=(cos(lsite[l].ang.th)*sin(lsite[l].ang.ph)*dthph[l]+
				sin(lsite[l].ang.th)*cos(lsite[l].ang.ph)*dthph[l+N_lattpoint]);

		spin0_z+=cos(lsite[l].ang.th)*SVAL1;
		spin0_z+=-sin(lsite[l].ang.th)*dthph[l];

		spin1_x+=-sin(lsite[m].ang.th)*cos(lsite[m].ang.ph)*SVAL2;
		spin1_x+=-(cos(lsite[m].ang.th)*cos(lsite[m].ang.ph)*dthph[m]-
				sin(lsite[m].ang.th)*sin(lsite[m].ang.ph)*dthph[m+N_lattpoint]);

		spin1_y+=sin(lsite[m].ang.th)*sin(lsite[m].ang.ph)*SVAL2;
		spin1_y+=(cos(lsite[m].ang.th)*sin(lsite[m].ang.ph)*dthph[m]+
				sin(lsite[m].ang.th)*cos(lsite[m].ang.ph)*dthph[m+N_lattpoint]);

		spin1_z+=cos(lsite[m].ang.th)*SVAL2;
		spin1_z+=-sin(lsite[m].ang.th)*dthph[m];

		spin2_x+=-sin(lsite[n].ang.th)*cos(lsite[n].ang.ph)*SVAL3;
		spin2_x+=-(cos(lsite[n].ang.th)*cos(lsite[n].ang.ph)*dthph[n]-
				sin(lsite[n].ang.th)*sin(lsite[n].ang.ph)*dthph[n+N_lattpoint]);

		spin2_y+=sin(lsite[n].ang.th)*sin(lsite[n].ang.ph)*SVAL3;
		spin2_y+=(cos(lsite[n].ang.th)*sin(lsite[n].ang.ph)*dthph[n]+
				sin(lsite[n].ang.th)*cos(lsite[n].ang.ph)*dthph[n+N_lattpoint]);

		spin2_z+=cos(lsite[n].ang.th)*SVAL3;
		spin2_z+=-sin(lsite[n].ang.th)*dthph[n];

	}
	}

	spin[0]=(vect){spin0_x,spin0_y,spin0_z};
	spin[1]=(vect){spin1_x,spin1_y,spin1_z};
	spin[2]=(vect){spin2_x,spin2_y,spin2_z};

	double VAL=0.0;

	for (p=0;p<3;p++) {

		VAL+=vectdot(spin[p],spin[p]);

	}

	vect chir={0.0,0.0,0.0};

	chir=vectsum(chir,vectcross(spin[0],spin[1]));
	chir=vectsum(chir,vectcross(spin[1],spin[2]));
	chir=vectsum(chir,vectcross(spin[2],spin[1]));

	free(spin);

	double SBL_MZ;
	cmplxvect SBL_CHIZ;
	SBL_MZ=12.0*(VAL+SMAGEX)/(N_lattpoint*N_lattpoint);
	SBL_CHIZ=re2cmplxvect(chir);
	SBL_CHIZ=cmplxvectsum(SBL_CHIZ,SCHIEX);
	SBL_CHIZ=realvectmult(72.0/(sqrt(3.0)
				*N_lattpoint*N_lattpoint),
			SBL_CHIZ);

	(*SBLMAG)+=SBL_MZ;
	(*SBLCHI)+=cmplxvectnorm(SBL_CHIZ);
	(*SBLMAG_VAR)+=SBL_MZ*SBL_MZ;
	(*SBLCHI_VAR)+=cmplxvectnormsqr(SBL_CHIZ);


}


void tri_cl_sublatmag(int N_lattpoint,lattpoint *lsite, 
		int **eltriangle,double *tri_CLSBL_ANG, 
		double *tri_CLSCAL_CHI,double *tri_CLSBL_CHI,
		double *tri_CLSBL_ANG_VAR, 
		double *tri_CLSCAL_CHI_VAR,
		double *tri_CLSBL_CHI_VAR) {

	int i,j,k,l,m,n,p,q;

	int L=sqrt(N_lattpoint);
	int DIM=2*N_lattpoint;
	int trinum=(L-1)*(L-1);

	vect *spin;

	spin=(vect *)calloc((ZCO/2),sizeof(vect));

	double VAL=0.0;
	double VAL_TMP=0.0;
	double VAL_VAR=0.0;
	vect ZEROVECT={0.0,0.0,0.0};
	vect CHIR={0.0,0.0,0.0};
	vect CHIRTEMP={0.0,0.0,0.0};
	double SCALAR_CHIR_TMP=0.0;
	double SCALAR_CHIR=0.0;
	double SCALAR_CHIR_VAR=0.0;

	for (i=0;i<trinum;i++) {

		l=eltriangle[i][0];
		m=eltriangle[i][1];
		n=eltriangle[i][2];

		spin[0].x=-S*sin(lsite[l].ang.th)*cos(lsite[l].ang.ph);
		spin[0].y=S*sin(lsite[l].ang.th)*sin(lsite[l].ang.ph);
		spin[0].z=S*cos(lsite[l].ang.th);

		spin[1].x=-S*sin(lsite[m].ang.th)*cos(lsite[m].ang.ph);
		spin[1].y=S*sin(lsite[m].ang.th)*sin(lsite[m].ang.ph);
		spin[1].z=S*cos(lsite[m].ang.th);

		spin[2].x=-S*sin(lsite[n].ang.th)*cos(lsite[n].ang.ph);
		spin[2].y=S*sin(lsite[n].ang.th)*sin(lsite[n].ang.ph);
		spin[2].z=S*cos(lsite[n].ang.th);

		VAL_TMP=vectnormsqr(
				vectsum(spin[0],
					vectsum(spin[1],spin[2])));
		VAL+=(1.0/3.0)*VAL_TMP;
		VAL_VAR+=VAL_TMP*VAL_TMP;

		CHIRTEMP=vectcross(spin[0],spin[1]);
		CHIRTEMP=vectsum(CHIRTEMP,vectcross(spin[1],spin[2]));
		CHIRTEMP=vectsum(CHIRTEMP,vectcross(spin[2],spin[0]));

		CHIR=vectsum(CHIR,unitvect(CHIRTEMP));

		SCALAR_CHIR_TMP=pow(8.0*vectdot(spin[0],vectcross(spin[1],spin[2])),2.0);
		SCALAR_CHIR+=SCALAR_CHIR_TMP;
		SCALAR_CHIR_VAR+=SCALAR_CHIR_TMP*SCALAR_CHIR_TMP;

	}

	double TSBLCHI;

	TSBLCHI=(1.0/trinum)*vectnorm(CHIR);

	(*tri_CLSBL_ANG)+=VAL;
	(*tri_CLSBL_CHI)+=TSBLCHI;
	(*tri_CLSCAL_CHI)+=SCALAR_CHIR;
	(*tri_CLSBL_ANG_VAR)+=VAL_VAR;
	(*tri_CLSBL_CHI_VAR)+=TSBLCHI*TSBLCHI;
	(*tri_CLSCAL_CHI_VAR)+=SCALAR_CHIR_VAR;

	free(spin);


}


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
		double *tri_SBLCHI_ex_VAR) {

	int i,j,k,l,m,n,nu;
	int IXP,JXP,KXP;
	int ix;
	int l1,l2,l3;

	int L=sqrt(N_lattpoint);
	int DIM=2*N_lattpoint;
	int trinum=(L-1)*(L-1);

	complex Y11={0.0,0.0};
	complex Y21={0.0,0.0};
	complex Y31={0.0,0.0};
	complex Y12={0.0,0.0};
	complex Y22={0.0,0.0};
	complex Y32={0.0,0.0};
	complex Y13={0.0,0.0};
	complex Y23={0.0,0.0};
	complex Y33={0.0,0.0};

	complex VAL2={0.0,0.0};
	complex VAL3={0.0,0.0};

	double ZVAL1,ZVAL2,ZVAL3;
	complex ZVTMP1,ZVTMP2,ZVTMP3;

	complex QZERO={0.0,0.0};
	cmplxvect QZEROVECT={{0.0,0.0},{0.0,0.0},{0.0,0.0}};

	cmplxvect SMAG1,SMAG2,SMAG3;
	cmplxvect STEMP1, STEMP2, STEMP3;

	int QDIM=N_lattpoint*trinum;
	cmplxvect *QSPIN1;
	QSPIN1=(cmplxvect *)malloc(QDIM*sizeof(cmplxvect));
	cmplxvect *QSPIN2;
	QSPIN2=(cmplxvect *)malloc(QDIM*sizeof(cmplxvect));
	cmplxvect *QSPIN3;
	QSPIN3=(cmplxvect *)malloc(QDIM*sizeof(cmplxvect));

	double *SMAGTEMP;
	SMAGTEMP=(double *)malloc(trinum*sizeof(double));
	cmplxvect *SCHITEMP;
	SCHITEMP=(cmplxvect *)malloc(trinum*sizeof(cmplxvect));
	complex *SCALCHITEMP;
	SCALCHITEMP=(complex *)malloc(trinum*sizeof(complex));

#pragma omp parallel private(ix,i,j,k,IXP,JXP,KXP,l1,l2,l3, \
		Y11,Y21,Y31,Y12,Y22,Y32,Y13,Y23,Y33, \
		VAL3,VAL2, \
		ZVAL1,ZVAL2,ZVAL3, \
		ZVTMP1,ZVTMP2,ZVTMP3, \
		STEMP1,STEMP2,STEMP3, \
		SMAG1,SMAG2,SMAG3)
	{
#pragma omp for
		for (ix=0;ix<trinum;ix++) {

			SMAGTEMP[ix]=0.0;
			SCALCHITEMP[ix]=QZERO;
			SCHITEMP[ix]=QZEROVECT;

			l1=eltriangle[ix][0];
			Y11.re=cos(lsite[l1].ang.th)*cos(lsite[l1].ang.ph);
			Y11.im=-sin(lsite[l1].ang.ph);
			Y21.re=-cos(lsite[l1].ang.th)*sin(lsite[l1].ang.ph);
			Y21.im=-cos(lsite[l1].ang.ph);
			Y31.re=sin(lsite[l1].ang.th);
			Y31.im=0.0;

			l2=eltriangle[ix][1];
			Y12.re=cos(lsite[l2].ang.th)*cos(lsite[l2].ang.ph);
			Y12.im=-sin(lsite[l2].ang.ph);
			Y22.re=-cos(lsite[l2].ang.th)*sin(lsite[l2].ang.ph);
			Y22.im=-cos(lsite[l2].ang.ph);
			Y32.re=sin(lsite[l2].ang.th);
			Y32.im=0.0;

			l3=eltriangle[ix][2];
			Y13.re=cos(lsite[l3].ang.th)*cos(lsite[l3].ang.ph);
			Y13.im=-sin(lsite[l3].ang.ph);
			Y23.re=-cos(lsite[l3].ang.th)*sin(lsite[l3].ang.ph);
			Y23.im=-cos(lsite[l3].ang.ph);
			Y33.re=sin(lsite[l3].ang.th);
			Y33.im=0.0;

			for (i=0;i<N_lattpoint;i++) {

				IXP=ix+i*trinum;

				VAL3=cmplxconj(T[l1+i*DIM]);	
				VAL2=cmplxconj(T[l1+N_lattpoint+i*DIM]);

				SMAG1.x=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y11,VAL2),cmplxmult(cmplxconj(Y11),VAL3)));
				SMAG1.y=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y21,VAL2),cmplxmult(cmplxconj(Y21),VAL3)));
				SMAG1.z=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y31,VAL2),cmplxmult(cmplxconj(Y31),VAL3)));

				QSPIN1[IXP]=SMAG1;

				VAL3=cmplxconj(T[l2+i*DIM]);	
				VAL2=cmplxconj(T[l2+N_lattpoint+i*DIM]);

				SMAG2.x=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y12,VAL2),cmplxmult(cmplxconj(Y12),VAL3)));
				SMAG2.y=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y22,VAL2),cmplxmult(cmplxconj(Y22),VAL3)));
				SMAG2.z=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y32,VAL2),cmplxmult(cmplxconj(Y32),VAL3)));


				QSPIN2[IXP]=SMAG2;

				VAL3=cmplxconj(T[l3+i*DIM]);	
				VAL2=cmplxconj(T[l3+N_lattpoint+i*DIM]);	

				SMAG3.x=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y13,VAL2),cmplxmult(cmplxconj(Y13),VAL3)));
				SMAG3.y=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y23,VAL2),cmplxmult(cmplxconj(Y23),VAL3)));
				SMAG3.z=realmult(sqrt(S/2.0),
						cmplxadd(cmplxmult(Y33,VAL2),cmplxmult(cmplxconj(Y33),VAL3)));

				QSPIN3[IXP]=SMAG3;

				SMAGTEMP[ix]+=cmplxvectnormsqr(
						cmplxvectsum(SMAG1,
							cmplxvectsum(SMAG2,SMAG3)));

				SCHITEMP[ix]=cmplxvectsum(SCHITEMP[ix],
						cmplxvectcross(SMAG1,cmplxconjvect(SMAG2)));
				SCHITEMP[ix]=cmplxvectsum(SCHITEMP[ix],
						cmplxvectcross(SMAG2,cmplxconjvect(SMAG3)));
				SCHITEMP[ix]=cmplxvectsum(SCHITEMP[ix],
						cmplxvectcross(SMAG3,cmplxconjvect(SMAG1)));

			}

			SMAG1=QZEROVECT;
			SMAG2=QZEROVECT;
			SMAG3=QZEROVECT;

			ZVAL1=S-alpha[l1][l1].re;
			ZVAL2=S-alpha[l2][l2].re;
			ZVAL3=S-alpha[l3][l3].re;
			
			SMAG1.x.re=-sin(lsite[l1].ang.th)*cos(lsite[l1].ang.ph)*ZVAL1;
			SMAG1.x.re+=-(cos(lsite[l1].ang.th)*cos(lsite[l1].ang.ph)*dthph[l1]-
					sin(lsite[l1].ang.th)*sin(lsite[l1].ang.ph)*dthph[l1+N_lattpoint]);

			SMAG1.y.re=sin(lsite[l1].ang.th)*sin(lsite[l1].ang.ph)*ZVAL1;
			SMAG1.y.re+=(cos(lsite[l1].ang.th)*sin(lsite[l1].ang.ph)*dthph[l1]+
					sin(lsite[l1].ang.th)*cos(lsite[l1].ang.ph)*dthph[l1+N_lattpoint]);

			SMAG1.z.re=cos(lsite[l1].ang.th)*ZVAL1;
			SMAG1.z.re+=-sin(lsite[l1].ang.th)*dthph[l1];

			SMAG2.x.re=-sin(lsite[l2].ang.th)*cos(lsite[l2].ang.ph)*ZVAL2;
			SMAG2.x.re+=-(cos(lsite[l2].ang.th)*cos(lsite[l2].ang.ph)*dthph[l2]-
					sin(lsite[l2].ang.th)*sin(lsite[l2].ang.ph)*dthph[l2+N_lattpoint]);

			SMAG2.y.re=sin(lsite[l2].ang.th)*sin(lsite[l2].ang.ph)*ZVAL2;
			SMAG2.y.re+=(cos(lsite[l2].ang.th)*sin(lsite[l2].ang.ph)*dthph[l2]+
					sin(lsite[l2].ang.th)*cos(lsite[l2].ang.ph)*dthph[l2+N_lattpoint]);

			SMAG2.z.re=cos(lsite[l2].ang.th)*ZVAL2;
			SMAG2.z.re+=-sin(lsite[l2].ang.th)*dthph[l2];

			SMAG3.x.re=-sin(lsite[l3].ang.th)*cos(lsite[l3].ang.ph)*ZVAL3;
			SMAG3.x.re+=-(cos(lsite[l3].ang.th)*cos(lsite[l3].ang.ph)*dthph[l3]-
					sin(lsite[l3].ang.th)*sin(lsite[l3].ang.ph)*dthph[l3+N_lattpoint]);

			SMAG3.y.re=sin(lsite[l3].ang.th)*sin(lsite[l3].ang.ph)*ZVAL3;
			SMAG3.y.re+=(cos(lsite[l3].ang.th)*sin(lsite[l3].ang.ph)*dthph[l3]+
					sin(lsite[l3].ang.th)*cos(lsite[l3].ang.ph)*dthph[l3+N_lattpoint]);

			SMAG3.z.re=cos(lsite[l3].ang.th)*ZVAL3;
			SMAG3.z.re+=-sin(lsite[l3].ang.th)*dthph[l3];

			SMAGTEMP[ix]+=cmplxvectnormsqr(
					cmplxvectsum(SMAG1,
						cmplxvectsum(SMAG2,SMAG3)));

			SCHITEMP[ix]=cmplxvectsum(SCHITEMP[ix],
					cmplxvectcross(SMAG1,cmplxconjvect(SMAG2)));
			SCHITEMP[ix]=cmplxvectsum(SCHITEMP[ix],
					cmplxvectcross(SMAG2,cmplxconjvect(SMAG3)));
			SCHITEMP[ix]=cmplxvectsum(SCHITEMP[ix],
					cmplxvectcross(SMAG3,cmplxconjvect(SMAG1)));

			SCALCHITEMP[ix]=cmplxadd(SCALCHITEMP[ix],
					cmplxvectdot(SMAG1,cmplxvectcross(SMAG2,SMAG3)));
			SCALCHITEMP[ix]=cmplxadd(SCALCHITEMP[ix],
					cmplxvectdot(SMAG2,cmplxvectcross(SMAG3,SMAG1)));
			SCALCHITEMP[ix]=cmplxadd(SCALCHITEMP[ix],
					cmplxvectdot(SMAG3,cmplxvectcross(SMAG1,SMAG2)));

		}
	}

	cmplxvect SCHIEX=QZEROVECT;

	double SMAGEX_TMP=0.0;
	double SMAGEX=0.0;
	double SMAGEX_VAR=0.0;
	double SCALCHIEX_TMP=0.0;
	double SCALCHIEX=0.0;
	double SCALCHIEX_VAR=0.0;

	for (ix=0;ix<trinum;ix++) {

		SMAGEX_TMP=(1.0/3.0)*SMAGTEMP[ix];
		SMAGEX+=SMAGEX_TMP;
		SMAGEX_VAR+=SMAGEX_TMP*SMAGEX_TMP;

		SCALCHIEX_TMP=cmplxmodsqr(realmult(8.0/3.0,SCALCHITEMP[ix]));
		SCALCHIEX+=SCALCHIEX_TMP;
		SCALCHIEX_VAR+=SCALCHIEX_TMP*SCALCHIEX_TMP;

		SCHIEX=cmplxvectsum(SCHIEX,
				unitcmplxvect(SCHITEMP[ix]));

	}

	double TSBLCHI;

	TSBLCHI=(1.0/trinum)*cmplxvectnorm(SCHIEX);

	(*tri_SBLMAG_ex)+=SMAGEX;
	(*tri_SBLCHI_ex)+=TSBLCHI;
	(*tri_SCALCHI_ex)+=SCALCHIEX;
	(*tri_SBLMAG_ex_VAR)+=SMAGEX_VAR;
	(*tri_SBLCHI_ex_VAR)+=TSBLCHI*TSBLCHI;
	(*tri_SCALCHI_ex_VAR)+=SCALCHIEX_VAR;

	free(SMAGTEMP);
	free(SCHITEMP);
	free(SCALCHITEMP);
	free(QSPIN1);
	free(QSPIN2);
	free(QSPIN3);
}



void im_chi(double womega,vect q,double GW,
		int N_lattpoint,lattpoint *lsite,
		double *dthph,complex *T,double *OMEGA,
		double *SUSC,double *SUSC_VAR,
		double *CSUSC,double *CSUSC_VAR) {

	int i,j,k,l,nu;

	int DIM=2*N_lattpoint;
	int L=sqrt(N_lattpoint);

	complex Y1={0.0,0.0};
	complex Y2={0.0,0.0};
	complex Y3={0.0,0.0};
	complex Y4={0.0,0.0};
	complex Y5={0.0,0.0};
	complex Y6={0.0,0.0};

	double VAL1=0.0;
	complex VAL2={0.0,0.0};
	complex VAL3={0.0,0.0};
	complex Z={0.0,0.0};

	double phase=0.0;

	complex QZERO={0.0,0.0};
	cmplxvect QZEROVECT={{0.0,0.0},{0.0,0.0},{0.0,0.0}};

	cmplxvect *SZEROMAT;
	SZEROMAT=(cmplxvect *)malloc(N_lattpoint*sizeof(cmplxvect));
	cmplxvect *SMAG;
	SMAG=(cmplxvect *)malloc(N_lattpoint*sizeof(cmplxvect));
	cmplxvect SZERO;
	cmplxvect *CZEROMAT;
	CZEROMAT=(cmplxvect *)malloc(N_lattpoint*sizeof(cmplxvect));
	cmplxvect CZERO;

#pragma omp parallel private(i,j,VAL1,VAL2,VAL3, \
		phase,Z,Y1,Y2,Y3,Y4,Y5,Y6)
	{
#pragma omp for
		for (j=0;j<N_lattpoint;j++) {

			SMAG[j]=QZEROVECT;
			SZEROMAT[j]=QZEROVECT;
			CZEROMAT[j]=QZEROVECT;

			for (i=0;i<N_lattpoint;i++) {


				VAL1=1.0*S/N_lattpoint;

				phase=vectdot(q,lsite[i].pos);
				Z.re=cos(phase);
				Z.im=-sin(phase);

				Y1.re=cos(lsite[i].ang.th)*cos(lsite[i].ang.ph);
				Y1.im=-sin(lsite[i].ang.ph);
				Y2.re=-cos(lsite[i].ang.th)*sin(lsite[i].ang.ph);
				Y2.im=-cos(lsite[i].ang.ph);
				Y3.re=sin(lsite[i].ang.th);
				Y3.im=0.0;

				VAL1-=pow(cmplxmod(T[i+j*DIM]),2.0);
				VAL3=cmplxconj(T[i+j*DIM]);	
				VAL2=cmplxconj(T[i+N_lattpoint+j*DIM]);	

				SMAG[j].x=cmplxadd(SMAG[j].x,cmplxmult(Z,
							cmplxadd(cmplxmult(Y1,VAL2),cmplxmult(cmplxconj(Y1),VAL3))));
				SMAG[j].y=cmplxadd(SMAG[j].y,cmplxmult(Z,
							cmplxadd(cmplxmult(Y2,VAL2),cmplxmult(cmplxconj(Y2),VAL3))));
				SMAG[j].z=cmplxadd(SMAG[j].z,cmplxmult(Z,
							cmplxadd(cmplxmult(Y3,VAL2),cmplxmult(cmplxconj(Y3),VAL3))));

				Y4.re=-sin(lsite[i].ang.th)*cos(lsite[i].ang.ph)*VAL1;
				Y4.im=0.0;

				Y5.re=sin(lsite[i].ang.th)*sin(lsite[i].ang.ph)*VAL1;
				Y5.im=0.0;

				Y6.re=cos(lsite[i].ang.th)*VAL1;
				Y6.im=0.0;

				Y4=cmplxmult(Z,Y4);
				Y5=cmplxmult(Z,Y5);
				Y6=cmplxmult(Z,Y6);
				SZEROMAT[j].x=cmplxadd(SZEROMAT[j].x,Y4);
				SZEROMAT[j].y=cmplxadd(SZEROMAT[j].y,Y5);
				SZEROMAT[j].z=cmplxadd(SZEROMAT[j].z,Y6);

				Y4.re=-sin(lsite[i].ang.th)*cos(lsite[i].ang.ph)*(1.0*S/N_lattpoint);
				Y4.im=0.0;

				Y5.re=sin(lsite[i].ang.th)*sin(lsite[i].ang.ph)*(1.0*S/N_lattpoint);
				Y5.im=0.0;

				Y6.re=cos(lsite[i].ang.th)*(1.0*S/N_lattpoint);
				Y6.im=0.0;

				Y4=cmplxmult(Z,Y4);
				Y5=cmplxmult(Z,Y5);
				Y6=cmplxmult(Z,Y6);
				CZEROMAT[j].x=cmplxadd(CZEROMAT[j].x,Y4);
				CZEROMAT[j].y=cmplxadd(CZEROMAT[j].y,Y5);
				CZEROMAT[j].z=cmplxadd(CZEROMAT[j].z,Y6);

			}

		}
	}

	SZERO=QZEROVECT;
	CZERO=QZEROVECT;

	for (j=0;j<N_lattpoint;j++) {

		phase=vectdot(q,lsite[j].pos);
		Z.re=cos(phase);
		Z.im=-sin(phase);
		Y4.re=-(cos(lsite[j].ang.th)*cos(lsite[j].ang.ph)*dthph[j]-
				sin(lsite[j].ang.th)*sin(lsite[j].ang.ph)*dthph[j+N_lattpoint]);
		Y4.im=0.0;
		Y5.re=(cos(lsite[j].ang.th)*sin(lsite[j].ang.ph)*dthph[j]+
				sin(lsite[j].ang.th)*cos(lsite[j].ang.ph)*dthph[j+N_lattpoint]);
		Y5.im=0.0;
		Y6.re=-sin(lsite[j].ang.th)*dthph[j];
		Y6.im=0.0;
		Y4=cmplxmult(Z,Y4);
		Y5=cmplxmult(Z,Y5);
		Y6=cmplxmult(Z,Y6);
		SZEROMAT[j].x=cmplxadd(SZEROMAT[j].x,Y4);
		SZEROMAT[j].y=cmplxadd(SZEROMAT[j].y,Y5);
		SZEROMAT[j].z=cmplxadd(SZEROMAT[j].z,Y6);

		SZERO=cmplxvectsum(SZERO,SZEROMAT[j]);
		CZERO=cmplxvectsum(CZERO,CZEROMAT[j]);

	}

	double MAGQ=0.0;
	double gweight=0.0;
	double MAGC=0.0;


	for (l=0;l<N_lattpoint;l++) {

		//Noting that the excited spectra is 2*OMEGA[l]
		gweight=exp(-pow(womega-fabs(2.0*OMEGA[l]),2.0)/(GW*GW))/(sqrt(2.0*pi)*GW);
		MAGQ+=gweight*
			((S/2.0)/(N_lattpoint*N_lattpoint))*(cmplxmodsqr(SMAG[l].x)+cmplxmodsqr(SMAG[l].y)+cmplxmodsqr(SMAG[l].z));
	}


	MAGQ+=exp(-pow(womega,2.0)/(GW*GW))/(sqrt(2.0*pi)*GW)*
		((S/2.0)/(N_lattpoint*N_lattpoint))*(cmplxmodsqr(SZERO.x)+cmplxmodsqr(SZERO.y)+cmplxmodsqr(SZERO.z));

	MAGC=exp(-pow(womega,2.0)/(GW*GW))/(sqrt(2.0*pi)*GW)*
		((S/2.0)/(N_lattpoint*N_lattpoint))*(cmplxmodsqr(CZERO.x)+cmplxmodsqr(CZERO.y)+cmplxmodsqr(CZERO.z));

	free(SZEROMAT);
	free(SMAG);
	free(CZEROMAT);

	(*SUSC)+=MAGQ;
	(*SUSC_VAR)+=MAGQ*MAGQ;
	(*CSUSC)+=MAGC;
	(*CSUSC_VAR)+=MAGC*MAGC;

}




void CHI_MAT(int N_lattpoint,lattpoint *lsite,
		double *dthph,complex *T,double *OMEGA,
		double ***CHIQ,double **SZERO) {

	int i,j,k,l,nu,ix;

	int m, n;

	int DIM=2*N_lattpoint;
	int L=sqrt(N_lattpoint);

	complex Y1={0.0,0.0};
	complex Y2={0.0,0.0};
	complex Y3={0.0,0.0};

	double VAL1=0.0;
	complex VAL2={0.0,0.0};
	complex VAL3={0.0,0.0};

	complex SMAGX;
	complex SMAGY;
	complex SMAGZ;

	fftw_complex *inX;
	fftw_complex *outX;
	fftw_plan plan_forwardX;

	fftw_complex *inY;
	fftw_complex *outY;
	fftw_plan plan_forwardY;

	fftw_complex *inZ;
	fftw_complex *outZ;
	fftw_plan plan_forwardZ;

	inX = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	inY = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	inZ = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	outX = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	outY = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	outZ = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);

	fftw_complex *INSZX;
	fftw_complex *OSZX;
	fftw_plan plan_forwardSZX;

	fftw_complex *INSZY;
	fftw_complex *OSZY;
	fftw_plan plan_forwardSZY;

	fftw_complex *INSZZ;
	fftw_complex *OSZZ;
	fftw_plan plan_forwardSZZ;

	INSZX = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	INSZY = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	INSZZ = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	OSZX = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	OSZY = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);
	OSZZ = fftw_malloc(sizeof(fftw_complex)*N_lattpoint);

	plan_forwardX = fftw_plan_dft_2d (L, L, inX, outX,
			FFTW_FORWARD, FFTW_ESTIMATE);
	plan_forwardY = fftw_plan_dft_2d (L, L, inY, outY,
			FFTW_FORWARD, FFTW_ESTIMATE);
	plan_forwardZ = fftw_plan_dft_2d (L, L, inZ, outZ,
			FFTW_FORWARD, FFTW_ESTIMATE);

	for (i=0;i<N_lattpoint;i++) {

		INSZX[i][0]=0.0;
		INSZX[i][1]=0.0;
		INSZY[i][0]=0.0;
		INSZY[i][1]=0.0;
		INSZZ[i][0]=0.0;
		INSZZ[i][1]=0.0;

	}

	for (j=0;j<N_lattpoint;j++) {

		for (i=0;i<N_lattpoint;i++) {

			VAL1=1.0*S/N_lattpoint;

			Y1.re=cos(lsite[i].ang.th)*cos(lsite[i].ang.ph);
			Y1.im=-sin(lsite[i].ang.ph);
			Y2.re=-cos(lsite[i].ang.th)*sin(lsite[i].ang.ph);
			Y2.im=-cos(lsite[i].ang.ph);
			Y3.re=sin(lsite[i].ang.th);
			Y3.im=0.0;

			VAL1-=pow(cmplxmod(T[i+j*DIM]),2.0);
			VAL3=cmplxconj(T[i+j*DIM]);	
			VAL2=cmplxconj(T[i+N_lattpoint+j*DIM]);	

			SMAGX=cmplxadd(cmplxmult(Y1,VAL2),cmplxmult(cmplxconj(Y1),VAL3));
			SMAGY=cmplxadd(cmplxmult(Y2,VAL2),cmplxmult(cmplxconj(Y2),VAL3));
			SMAGZ=cmplxadd(cmplxmult(Y3,VAL2),cmplxmult(cmplxconj(Y3),VAL3));

			inX[i][0]=SMAGX.re;
			inX[i][1]=SMAGX.im;
			inY[i][0]=SMAGY.re;
			inY[i][1]=SMAGY.im;
			inZ[i][0]=SMAGZ.re;
			inZ[i][1]=SMAGZ.im;

			INSZX[i][0]+=-sin(lsite[i].ang.th)*cos(lsite[i].ang.ph)*VAL1;
			INSZX[i][1]+=0.0;

			INSZY[i][0]+=sin(lsite[i].ang.th)*sin(lsite[i].ang.ph)*VAL1;
			INSZY[i][1]+=0.0;

			INSZZ[i][0]+=cos(lsite[i].ang.th)*VAL1;
			INSZZ[i][1]+=0.0;

		}

		INSZX[j][0]+=-(cos(lsite[j].ang.th)*cos(lsite[j].ang.ph)*dthph[j]-
				sin(lsite[j].ang.th)*sin(lsite[j].ang.ph)*dthph[j+N_lattpoint]);
		INSZY[j][0]+=(cos(lsite[j].ang.th)*sin(lsite[j].ang.ph)*dthph[j]+
				sin(lsite[j].ang.th)*cos(lsite[j].ang.ph)*dthph[j+N_lattpoint]);
		INSZZ[j][0]+=-sin(lsite[j].ang.th)*dthph[j];

		fftw_execute ( plan_forwardX );

		fftw_execute ( plan_forwardY );

		fftw_execute ( plan_forwardZ );


		for (i=0;i<N_lattpoint;i++) {

			n=i%(L);
			m=i/(L);

			CHIQ[j][m][n]=((S/2.0)/(N_lattpoint*N_lattpoint))*
				(pow(outX[i][0],2.0)+pow(outX[i][1],2.0)
				 +pow(outY[i][0],2.0)+pow(outY[i][1],2.0)
				 +pow(outZ[i][0],2.0)+pow(outZ[i][1],2.0));


		}


	}

	fftw_destroy_plan ( plan_forwardX);
	fftw_destroy_plan ( plan_forwardY);
	fftw_destroy_plan ( plan_forwardZ);

	fftw_free(inX);
	fftw_free(inY);
	fftw_free(inZ);
	fftw_free(outX);
	fftw_free(outY);
	fftw_free(outZ);

	plan_forwardSZX = fftw_plan_dft_2d (L, L, INSZX, OSZX,
			FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute ( plan_forwardSZX );

	plan_forwardSZY = fftw_plan_dft_2d (L, L, INSZY, OSZY,
			FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute ( plan_forwardSZY );

	plan_forwardSZZ = fftw_plan_dft_2d (L, L, INSZZ, OSZZ,
			FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute ( plan_forwardSZZ );

	fftw_destroy_plan ( plan_forwardSZX);
	fftw_destroy_plan ( plan_forwardSZY);
	fftw_destroy_plan ( plan_forwardSZZ);

	for (i=0;i<N_lattpoint;i++) {

		n=i%(L);
		m=i/(L);

		SZERO[m][n]=((S/2.0)/(N_lattpoint*N_lattpoint))*
			(pow(OSZX[i][0],2.0)+pow(OSZX[i][1],2.0)
			 +pow(OSZY[i][0],2.0)+pow(OSZY[i][1],2.0)
			 +pow(OSZZ[i][0],2.0)+pow(OSZZ[i][1],2.0));



	}

	fftw_free(INSZX);
	fftw_free(INSZY);
	fftw_free(INSZZ);
	fftw_free(OSZX);
	fftw_free(OSZY);
	fftw_free(OSZZ);

	fftw_cleanup();

}



void eigenview(int EIGENCOUNT,int N_lattpoint,complex *T,double *OMEGA,double *eigen_grid,double *eigen_grid_var) {


	int i,j;
	int DIM=2*N_lattpoint;
	double VAL;

	for (i=0;i<N_lattpoint;i++) {

		VAL=cmplxmodsqr(T[i+(N_lattpoint-EIGENCOUNT)*DIM]);

		eigen_grid[i]+=VAL;
		eigen_grid_var[i]+=VAL*VAL;


	}


}

void IPR_CALC(double IPRMIN,double IPRMAX,int IPRGRID,double GW,double BTOL,int N_lattpoint,complex *T,double *OMEGA,double *IPR,double *IPR_var) {



	int i,j,k;
	int DIM=2*N_lattpoint;
	double gweight;
	double womega;
	double IVAL;
	double NUVAL;
	double DEVAL;
	double VNORM;
	double VAL;
	double IPRVAL;
	double MAXW;
	double MINW;
	MAXW=0.0;
	MINW=1.0/gtol;
	for (i=0;i<N_lattpoint;i++) {
		if (MAXW<fabs(OMEGA[i])) {
			MAXW=fabs(OMEGA[i]);
		}
		if (MINW>fabs(OMEGA[i]) && fabs(OMEGA[i])>BTOL) {
			MINW=fabs(OMEGA[i]);
		}
	}
#pragma omp parallel private(i,j,k,womega,DEVAL,NUVAL,\
		IVAL,VNORM,gweight,VAL,IPRVAL)
	{
#pragma omp for		
		for (i=0;i<IPRGRID;i++) {
			womega=IPRMIN+((IPRMAX-IPRMIN)/IPRGRID)*i;
			DEVAL=0.0;
			NUVAL=0.0;
			IPRVAL=0.0;
			if (womega>=MINW && womega<=2.0*MAXW) {
				for (k=0;k<N_lattpoint;k++) {
					if (fabs(OMEGA[k])>BTOL) {
						gweight=exp(-pow(womega-fabs(2.0*OMEGA[k]),2.0)/(GW*GW))/(sqrt(2.0*pi)*GW);
						IVAL=0.0;
						VNORM=0.0;
						for (j=0;j<N_lattpoint;j++) {
							VAL=cmplxmodsqr(T[j+k*DIM]);
							VNORM+=VAL;
							IVAL+=VAL*VAL;
						}
						NUVAL+=gweight*IVAL/(VNORM*VNORM);
						DEVAL+=gweight;
					}
				}
				IPRVAL=NUVAL/DEVAL;
			}
			IPR[i]+=IPRVAL;
			IPR_var[i]+=IPRVAL*IPRVAL;

		}
	}

}



void histo_angles(int N_lattpoint,lattpoint *lsite,double *dthph,
		double *cangle,double *cangle_var,
		double *qangle,double *qangle_var) {

	int i,j,k,nj;

	int *connect;
	connect=(int *)calloc(N_lattpoint*N_lattpoint,sizeof(int));

	nj=0;

	double Omega,th_i,th_j,ph_i,ph_j;
	double dQomega;

	for (i=0;i<N_lattpoint;i++) {

		th_i=lsite[i].ang.th;
		ph_i=lsite[i].ang.ph;
		for (k=0;k<ZCO;k++) {

			j=lsite[i].neighbour[k];
			th_j=lsite[j].ang.th;
			ph_j=lsite[j].ang.ph;

			if (connect[i+j*N_lattpoint]<1) {
				
				Omega=cos(th_i)*cos(th_j)+
					sin(th_i)*sin(th_j)*cos(ph_i-ph_j);
				cangle[nj]+=Omega;
				cangle_var[nj]+=Omega*Omega;
				dQomega=(1.0/S)*(
					cos(ph_i-ph_j)*(sin(th_i)*cos(th_j)*dthph[j]
						+cos(th_i)*sin(th_j)*dthph[i])
					-sin(ph_i-ph_j)*sin(th_i)*sin(th_j)*(dthph[i+N_lattpoint]
						-dthph[j+N_lattpoint])
					-(sin(th_i)*cos(th_j)*dthph[i]+cos(th_i)*sin(th_j)*dthph[j]));
				Omega+=dQomega;
				qangle[nj]+=Omega;
				qangle_var[nj]+=Omega*Omega;
				nj++;
				connect[i+j*N_lattpoint]=2;
				connect[j+i*N_lattpoint]=2;

			}

		}

	}

	free(connect);

}

