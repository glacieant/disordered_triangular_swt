#cython: boundscheck=False, wraparound=False, nonecheck=False
#cython: cdivision=True, profile=True

### This is a program with additional    ###
### algorithms required for computation. ###

# the zero temperature monte carlo algorithm
def create(long NSYS,long ZCO,
        double S,long[:,:] nbr,
        double HFIELD,double[:,:] J,
        double[:,:] ham):

    # declaring C datatypes for loop variables
    cdef long i, j, k

    # looping over matrix dimension
    for i in range(0,NSYS):

        #diagonal entry
        ham[i,i] = HFIELD

        #looping over neighbours
        for k in range(0,ZCO):

            #neighbour index
            j = nbr[i,k]

            #update to the diagonal value
            ham[i,i] -= J[i,j]*S

            #setting off diagonal value
            ham[i,j] = 0.5*J[i,j]*S



