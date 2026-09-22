#include "trilattice_sw.h"

// Matrix components

double Fxx(int i, int j, double **J,lattpoint *lsite) {

	double VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL=J[i][j]*(cos(th_i)*cos(th_j)*cos(ph_i)*cos(ph_j)
			+cos(th_i)*cos(th_j)*sin(ph_i)*sin(ph_j)
			+sin(th_i)*sin(th_j));

	return VAL;

}

double Fxz(int i, int j, double **J,lattpoint *lsite) {

	double VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL=J[i][j]*(-cos(th_i)*sin(th_j)*cos(ph_i)*cos(ph_j)
			-cos(th_i)*sin(th_j)*sin(ph_i)*sin(ph_j)
			+sin(th_i)*cos(th_j));

	return VAL;

}

double Fyy(int i, int j, double **J,lattpoint *lsite) {

	double VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL=J[i][j]*(cos(ph_i)*cos(ph_j)+sin(ph_i)*sin(ph_j));

	return VAL;

}

double Fyz(int i, int j, double **J,lattpoint *lsite) {

	double VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL=-J[i][j]*(sin(th_j)*sin(ph_i-ph_j));

	return VAL;

}

double Fzz(int i, int j, double **J,lattpoint *lsite) {

	double VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL=J[i][j]*(sin(th_i)*sin(th_j)*cos(ph_i)*cos(ph_j)
			+sin(th_i)*sin(th_j)*sin(ph_i)*sin(ph_j)
			+cos(th_i)*cos(th_j));

	return VAL;

}

double Fxy(int i, int j, double **J,lattpoint *lsite) {

	double VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL=J[i][j]*(cos(th_i)*cos(ph_i)*sin(ph_j)
			-cos(th_i)*sin(ph_i)*cos(ph_j));

	return VAL;

}

complex dthj_Fz(int i, int j, double **J,lattpoint *lsite) {

	complex VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL.re=J[i][j]*(-cos(th_i)*cos(th_j)*cos(ph_i)*cos(ph_j)
			-cos(th_i)*cos(th_j)*sin(ph_i)*sin(ph_j)
			-sin(th_i)*sin(th_j));
	VAL.im=J[i][j]*(cos(th_j)*sin(ph_i-ph_j));

	return VAL;

}

complex dphj_Fz(int i, int j, double **J,lattpoint *lsite) {

	complex VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL.re=J[i][j]*(cos(th_i)*cos(th_j)*cos(ph_i)*sin(ph_j)
			-cos(th_i)*cos(th_j)*sin(ph_i)*cos(ph_j));
	VAL.im=J[i][j]*(-sin(th_j)*cos(ph_i-ph_j));

	return VAL;

}

complex dthi_Fz(int i, int j, double **J,lattpoint *lsite) {

	complex VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL.re=J[i][j]*(+sin(th_i)*sin(th_j)*cos(ph_i)*cos(ph_j)
			+sin(th_i)*sin(th_j)*sin(ph_i)*sin(ph_j)
			+cos(th_i)*cos(th_j));
	VAL.im=0.0;

	return VAL;

}

complex dphi_Fz(int i, int j, double **J,lattpoint *lsite) {

	complex VAL;

	double th_i=lsite[i].ang.th;
	double th_j=lsite[j].ang.th;
	double ph_i=lsite[i].ang.ph;
	double ph_j=lsite[j].ang.ph;


	VAL.re=J[i][j]*(cos(th_i)*cos(th_j)*sin(ph_i)*cos(ph_j)
			-cos(th_i)*cos(th_j)*cos(ph_i)*sin(ph_j));
	VAL.im=J[i][j]*(sin(th_j)*cos(ph_i-ph_j));

	return VAL;

}


complex A(int i, int j, double **J,lattpoint *lsite) {

	complex VAL;

	int k,l,nbr;

	VAL.re=0.0;
	VAL.im=0.0;

	if (i==j) {

		for (nbr=0;nbr<ZCO;nbr++) {

			k=lsite[i].neighbour[nbr];
			VAL.re-=0.5*S*Fzz(i,k,J,lsite);
			k=lsite[i].nn_nbr[nbr];
			VAL.re-=0.5*S*Fzz(i,k,J,lsite);


		}	

	}

	else {

		VAL.re=0.25*S*(Fxx(i,j,J,lsite)+Fyy(i,j,J,lsite));
		VAL.im=-0.25*S*(Fxy(i,j,J,lsite)-Fxy(j,i,J,lsite));

	}	

	return VAL;

}


complex B(int i, int j, double **J,lattpoint *lsite) {

	complex VAL;

	int k,l,nbr;

	VAL.re=0.0;
	VAL.im=0.0;


	VAL.re=0.25*S*(Fxx(i,j,J,lsite)-Fyy(i,j,J,lsite));
	VAL.im=0.25*S*(Fxy(i,j,J,lsite)+Fxy(j,i,J,lsite));

	return VAL;

}



// Spin-wave Hamiltonian matrix



void gen_ham(int MDIM,lattpoint *lsite, double **J,complex *M) {

	int N_lp=MDIM/2;

	int i,j,k,l,p,q;

#pragma omp parallel private(p,q,i,j)
	{
#pragma omp for
		for (p=0;p<MDIM;p++) {

			for (q=0;q<MDIM;q++) {

				M[p+q*MDIM].re=0.0;
				M[p+q*MDIM].im=0.0;

				if (p<N_lp && q<N_lp) {

					M[p+q*MDIM]=A(p,q,J,lsite);
				} 

				else if (p>=N_lp && q>=N_lp) {

					i=p-N_lp;
					j=q-N_lp;

					M[p+q*MDIM]=A(j,i,J,lsite);

				} 

				else if (p<N_lp && q>=N_lp) {

					j=q-N_lp;

					M[p+q*MDIM]=B(p,j,J,lsite);

				} 

				else if (p>=N_lp && q<N_lp) {

					i=p-N_lp;

					M[p+q*MDIM]=cmplxconj(B(i,q,J,lsite));

				} 

			}

		}
	}

}


// Mean-field parameters


void mf_params(int N_lattpoint,lattpoint *lsite,complex *T,
	       complex **alpha,complex **beta) {

	int i,j,k,p,q;

	int MDIM=N_lattpoint*2;
	complex QZERO={0.0,0.0};

#pragma omp parallel private(p,q,i,j,k)
	{
#pragma omp for
		for (i=0;i<N_lattpoint;i++) {

			alpha[i][i]=QZERO;
			beta[i][i]=QZERO;

			for (k=0;k<N_lattpoint;k++) {

				alpha[i][i].re+=cmplxmodsqr(T[i+k*MDIM]);
				beta[i][i]=cmplxadd(beta[i][i],
						cmplxmult(cmplxconj(T[i+N_lattpoint+k*MDIM]),
							T[i+k*MDIM]));

			}
			for (j=0;j<ZCO;j++) {

				p=lsite[i].neighbour[j];
				q=lsite[i].nn_nbr[j];

				alpha[i][p]=QZERO;
				beta[i][p]=QZERO;
				alpha[i][q]=QZERO;
				beta[i][q]=QZERO;

				for (k=0;k<N_lattpoint;k++) {

					alpha[i][p]=cmplxadd(alpha[i][p],
							cmplxmult(cmplxconj(T[i+k*MDIM]),
								T[p+k*MDIM]));
					beta[i][p]=cmplxadd(beta[i][p],
							cmplxmult(cmplxconj(T[i+N_lattpoint+k*MDIM]),
								T[p+k*MDIM]));

					alpha[i][q]=cmplxadd(alpha[i][q],
							cmplxmult(cmplxconj(T[i+k*MDIM]),
								T[q+k*MDIM]));
					beta[i][q]=cmplxadd(beta[i][q],
							cmplxmult(cmplxconj(T[i+N_lattpoint+k*MDIM]),
								T[q+k*MDIM]));

				}

			}

		}

	}

}


void mf_angles(int N_lattpoint,lattpoint *lsite,double **J,
		complex **alpha,complex **beta,double BTOL,double *dthph) {


	int i,j,k,p,q;

	int MDIM=N_lattpoint*2;
	complex QZERO={0.0,0.0};
	complex ii={0.0,1.0};
	complex x;

	double *mat;
	mat=(double *)calloc(MDIM*MDIM,sizeof(double));
	double *vec;
	vec=(double *)calloc(MDIM,sizeof(double));

	// Matrix-vector linear system

#pragma omp parallel private(i,j,p,q,x)
	{
#pragma omp for
		for (i=0;i<N_lattpoint;i++) {

			for (j=0;j<ZCO;j++) {

				p=lsite[i].neighbour[j];
				q=lsite[i].nn_nbr[j];

				// NN contribution

				// Right-hand vector

				x=cmplxadd(alpha[p][p],realmult(1.0/4.0,
							cmplxadd(realmult(2.0,alpha[i][i]),
								cmplxconj(beta[i][i]))));
				vec[i]+=-1.0*Fxz(i,p,J,lsite)*x.re;
				vec[i+N_lattpoint]+=-1.0*Fxz(i,p,J,lsite)*x.im;

				x=cmplxadd(alpha[i][p],cmplxconj(beta[i][p]));
				vec[i]+=-1.0*Fxz(p,i,J,lsite)*x.re;
				vec[i+N_lattpoint]+=-1.0*Fxz(p,i,J,lsite)*x.im;

				x=cmplxadd(alpha[p][p],realmult(1.0/4.0,
							cmplxadd(realmult(2.0,alpha[i][i]),
								realmult(-1.0,cmplxconj(beta[i][i])))));
				x=cmplxmult(ii,x);
				vec[i]+=Fyz(i,p,J,lsite)*x.re;
				vec[i+N_lattpoint]+=Fyz(i,p,J,lsite)*x.im;

				x=cmplxadd(alpha[i][p],realmult(-1.0,cmplxconj(beta[i][p])));
				x=cmplxmult(ii,x);
				vec[i]+=Fyz(p,i,J,lsite)*x.re;
				vec[i+N_lattpoint]+=Fyz(p,i,J,lsite)*x.im;

				// Matrix elements

				// Neighbouring spin-angle derivative

				x=dthj_Fz(i,p,J,lsite);
				mat[i+p*MDIM]=x.re;
				mat[i+N_lattpoint+p*MDIM]=x.im;


				x=dphj_Fz(i,p,J,lsite);
				mat[i+(p+N_lattpoint)*MDIM]=x.re;
				mat[i+N_lattpoint+(p+N_lattpoint)*MDIM]=x.im;

				// Self derivative

				x=dthi_Fz(i,p,J,lsite);
				mat[i+i*MDIM]+=x.re;
				mat[i+N_lattpoint+i*MDIM]+=x.im;

				x=dphi_Fz(i,p,J,lsite);
				mat[i+(i+N_lattpoint)*MDIM]+=x.re;
				mat[i+N_lattpoint+(i+N_lattpoint)*MDIM]+=x.im;


				// NNN contribution

				// Right-hand vector

				x=cmplxadd(alpha[q][q],realmult(1.0/4.0,
							cmplxadd(realmult(2.0,alpha[i][i]),
								cmplxconj(beta[i][i]))));
				vec[i]+=-1.0*Fxz(i,q,J,lsite)*x.re;
				vec[i+N_lattpoint]+=-1.0*Fxz(i,q,J,lsite)*x.im;

				x=cmplxadd(alpha[i][q],cmplxconj(beta[i][q]));
				vec[i]+=-1.0*Fxz(q,i,J,lsite)*x.re;
				vec[i+N_lattpoint]+=-1.0*Fxz(q,i,J,lsite)*x.im;

				x=cmplxadd(alpha[q][q],realmult(1.0/4.0,
							cmplxadd(realmult(2.0,alpha[i][i]),
								realmult(-1.0,cmplxconj(beta[i][i])))));
				x=cmplxmult(ii,x);
				vec[i]+=Fyz(i,q,J,lsite)*x.re;
				vec[i+N_lattpoint]+=Fyz(i,q,J,lsite)*x.im;

				x=cmplxadd(alpha[i][q],realmult(-1.0,cmplxconj(beta[i][q])));
				x=cmplxmult(ii,x);
				vec[i]+=Fyz(q,i,J,lsite)*x.re;
				vec[i+N_lattpoint]+=Fyz(q,i,J,lsite)*x.im;

				// Matrix elements

				// Neighbouring spin-angle derivative

				x=dthj_Fz(i,q,J,lsite);
				mat[i+q*MDIM]=x.re;
				mat[i+N_lattpoint+q*MDIM]=x.im;


				x=dphj_Fz(i,q,J,lsite);
				mat[i+(q+N_lattpoint)*MDIM]=x.re;
				mat[i+N_lattpoint+(q+N_lattpoint)*MDIM]=x.im;

				// Self derivative

				x=dthi_Fz(i,q,J,lsite);
				mat[i+i*MDIM]+=x.re;
				mat[i+N_lattpoint+i*MDIM]+=x.im;

				x=dphi_Fz(i,q,J,lsite);
				mat[i+(i+N_lattpoint)*MDIM]+=x.re;
				mat[i+N_lattpoint+(i+N_lattpoint)*MDIM]+=x.im;



			}

			// Left-hand matrix

			vec[i]=-vec[i];
			vec[i+N_lattpoint]=-vec[i+N_lattpoint];

		}
	}

	// Singular value decomposition

	double *sval;
	double *u;
	double *vt;
	sval=(double *)malloc(MDIM*sizeof(double));
	u=(double *)malloc(MDIM*MDIM*sizeof(double));
	vt=(double *)malloc(MDIM*MDIM*sizeof(double));

	int lwork, info;
	double wkopt;
	double *work;
	lwork=-1;
	dgesvd_("A","A",&MDIM,&MDIM,mat,&MDIM,sval,
			u,&MDIM,vt,&MDIM,&wkopt,&lwork,&info);
	lwork=(int)wkopt;
	work=(double *)malloc(lwork*sizeof(double));
	dgesvd_("A","A",&MDIM,&MDIM,mat,&MDIM,sval,
			u,&MDIM,vt,&MDIM,work,&lwork,&info);
	free(work);

	for (i=0;i<MDIM;i++) {
		if (fabs(sval[i])>BTOL) {
			sval[i]=1.0/sval[i];	
		}
		else {
			sval[i]=0.0;	
		}
	}

	double *ut;
	ut=(double *)malloc(MDIM*MDIM*sizeof(double));

#pragma omp parallel private(i,j)
	{
#pragma omp for
	for (i=0;i<MDIM;i++) {
		for (j=0;j<MDIM;j++) {
			ut[i+j*MDIM]=sval[i]*u[j+i*MDIM];
		}
	}
	}

	double a=1.0;
	double b=0.0;
	dgemm_("T","N",&MDIM,&MDIM,&MDIM,&a,vt,&MDIM,ut,&MDIM,&b,mat,&MDIM);

	int incx=1;
	int incy=1;
	dgemv_("N",&MDIM,&MDIM,&a,mat,&MDIM,vec,&incx,&b,dthph,&incy);

	free(sval);
	free(u);
	free(ut);
	free(vt);
	free(mat);
	free(vec);


}


