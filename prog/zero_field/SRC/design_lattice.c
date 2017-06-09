#include "trilattice_sw.h"


// Takes as an argument one of the sublattice type 0, 1 or 2 and
// returns the next sublattice type in cyclic order

int subcycle(int i) {

	int j;

	if (i==0) {
		j=1;
	}

	else if (i==1) {
		j=2;
	}

	else if (i==2) {
		j=0;
	}

	return j;

}

void create_trilattice(int N_lattpoint, lattpoint *lsite) {

	int i,j,k,l,p,q,NUM;

	int L=sqrt(N_lattpoint);

	//Assigning the lattice geometry with 
	//the nearest and the next nearest neigbour

	for (i=0;i<N_lattpoint;i++) {

		//Automated assignment of 6 neighbours and 6
		//nearest neighbours for each lattice point

		if (floor((1.0*i-1.0)/L)==floor((1.0*i)/L)) {
			lsite[i].neighbour[0]=i-1;
		}
		else if (floor((1.0*i-1.0)/L)!=floor((1.0*i)/L)) {
			lsite[i].neighbour[0]=i-1+L;
		}
		if (floor((1.0*i+1.0)/L)==floor((1.0*i)/L)) {
			lsite[i].neighbour[1]=i+1;
		}
		else if (floor((1.0*i+1.0)/L)!=floor((1.0*i)/L)) {
			lsite[i].neighbour[1]=i+1-L;
		}

		if ((i+L)>=N_lattpoint) {
			NUM=i+L-L*L;		

			lsite[i].neighbour[2]=NUM;
			if (floor((1.0*NUM-1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].neighbour[3]=NUM-1;
			}
			else if (floor((1.0*NUM-1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].neighbour[3]=NUM-1+L;
			}
	
			if (floor((1.0*NUM+1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[0]=NUM+1;
			}
			else if (floor((1.0*NUM+1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[0]=NUM+1-L;
			}
			if (floor((1.0*NUM-2.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[1]=NUM-2;
			}
			else if (floor((1.0*NUM-2.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[1]=NUM-2+L;
			}

		}
		else if ((i+L)<N_lattpoint) {
			NUM=i+L;		
			lsite[i].neighbour[2]=NUM;
			if (floor((1.0*NUM-1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].neighbour[3]=NUM-1;
			}
			else if (floor((1.0*NUM-1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].neighbour[3]=NUM-1+L;
			}

			if (floor((1.0*NUM+1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[0]=NUM+1;
			}
			else if (floor((1.0*NUM+1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[0]=NUM+1-L;
			}
			if (floor((1.0*NUM-2.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[1]=NUM-2;
			}
			else if (floor((1.0*NUM-2.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[1]=NUM-2+L;
			}

		}
		if ((i-L)<0) {
			NUM=i-L+L*L;		
			lsite[i].neighbour[4]=NUM;
			if (floor((1.0*NUM+1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].neighbour[5]=NUM+1;
			}
			else if (floor((1.0*NUM+1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].neighbour[5]=NUM+1-L;
			}	

			if (floor((1.0*NUM+2.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[2]=NUM+2;
			}
			else if (floor((1.0*NUM+2.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[2]=NUM+2-L;
			}	
			if (floor((1.0*NUM-1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[3]=NUM-1;
			}
			else if (floor((1.0*NUM-1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[3]=NUM-1+L;
			}	

		}
		else if ((i-L)>=0) {
			NUM=i-L;		
			lsite[i].neighbour[4]=NUM;
			if (floor((1.0*NUM+1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].neighbour[5]=NUM+1;
			}
			else if (floor((1.0*NUM+1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].neighbour[5]=NUM+1-L;
			}

			if (floor((1.0*NUM+2.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[2]=NUM+2;
			}
			else if (floor((1.0*NUM+2.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[2]=NUM+2-L;
			}	
			if (floor((1.0*NUM-1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[3]=NUM-1;
			}
			else if (floor((1.0*NUM-1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[3]=NUM-1+L;
			}	

		}

		if ((i+2*L)>=N_lattpoint) {
			NUM=i+2*L-L*L;		
			if (floor((1.0*NUM-1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[4]=NUM-1;
			}
			else if (floor((1.0*NUM-3.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[4]=NUM-1+L;
			}	
		}
		else if ((i+2*L)<N_lattpoint) {
			NUM=i+2*L;		
			if (floor((1.0*NUM-1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[4]=NUM-1;
			}
			else if (floor((1.0*NUM-3.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[4]=NUM-1+L;
			}	
		}
		if ((i-2*L)<0) {
			NUM=i-2*L+L*L;		
			if (floor((1.0*NUM+1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[5]=NUM+1;
			}
			else if (floor((1.0*NUM+1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[5]=NUM+1-L;
			}	
		}
		else if ((i-2*L)>=0) {
			NUM=i-2*L;		
			if (floor((1.0*NUM+1.0)/L)==floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[5]=NUM+1;
			}
			else if (floor((1.0*NUM+1.0)/L)!=floor((1.0*NUM)/L)) {
				lsite[i].nn_nbr[5]=NUM+1-L;
			}	
		}

		//Assigning co-ordinates of the lattice

		lsite[i].pos.x=(i/L)*latconst*1.0/2.0+(i%L)*latconst;
		lsite[i].pos.y=floor((1.0*i)/L)*latconst*sqrt(3.0)/2.0;
		lsite[i].pos.z=0.0;

	}

}

void create_auxlattice(int N_lattpoint,int **sublattice, int **eltriangle) {

	int i,j,k,l,p,q,NUM;

	int L=sqrt(N_lattpoint);

	int CYCL=0;
	int *SUBNUM;

	SUBNUM=(int *)malloc(3*sizeof(int));

	SUBNUM[0]=0;
	SUBNUM[1]=0;
	SUBNUM[2]=0;

	for (i=0;i<N_lattpoint;i++) {

		//Adding the lattice point to appropriate sublattice directory		

		sublattice[CYCL][SUBNUM[CYCL]]=i;
		SUBNUM[CYCL]+=1;

		if (((i+1)%L)!=0) {
			CYCL=subcycle(CYCL);
		}


	}

	for (i=0;i<(L-1)*(L-1);i++) {

		j=i+i/(L-1);
		eltriangle[i][0]=j;
		eltriangle[i][1]=j+1;
		eltriangle[i][2]=j+L;


	}

	free(SUBNUM);
}

