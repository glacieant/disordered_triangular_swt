#cython: boundscheck=False, wraparound=False, nonecheck=False
#cython: cdivision=True, profile=True

### This is a program with additional    ###
### algorithms required for computation. ###

# the zero temperature monte carlo algorithm
def classic_zmc(long CLNUM,long NSYS,long[:,:] nbr,double[:,:] J,
        double[:,:,:] M,double TOL):

    # declaring C datatypes for loop variables
    cdef long i, j, k, p
    cdef long MAX_ITER = CLNUM
    cdef double SQTOL = 0.0
    cdef double TMPTOL = 0.0
    cdef double ITERTOL = 0.0
    cdef double VNORM
    cdef double LOCVEC0,LOCVEC1,LOCVEC2

    for i in range(0,MAX_ITER):

        SQTOL = 0.0
        TMPTOL = 0.0
        for j in range(0,NSYS):

            LOCVEC0 = 0.0
            LOCVEC1 = 0.0
            LOCVEC2 = 0.0
            for p in range(0,12):
                k = nbr[j,p]
                LOCVEC0 += J[j,k]*M[1,k,0]
                LOCVEC1 += J[j,k]*M[1,k,1]
                #LOCVEC2 += J[j,k]*M[1,k,2]

            VNORM = (LOCVEC0**2
                    +LOCVEC1**2
                    #+LOCVEC2**2
                    )**0.5
            if VNORM >= 10.0**(-6.0):
                M[1,j,0] = -LOCVEC0/VNORM
                M[1,j,1] = -LOCVEC1/VNORM
                #M[1,j,2] = -LOCVEC2/VNORM
            
            else:
                M[1,j,0] = 1.0
                M[1,j,1] = 0.0
                #M[1,j,2] = 0.0

            ITERTOL = ((M[1,j,0]-M[0,j,0])**2
                    +(M[1,j,1]-M[0,j,1])**2
                    #+(M[1,j,2]-M[0,j,2])**2
                    )**0.5

            if ITERTOL>TMPTOL:
                TMPTOL = ITERTOL
            SQTOL += ITERTOL**2
            M[0,j,0] = M[1,j,0]
            M[0,j,1] = M[1,j,1]
            #M[0,j,2] = M[1,j,2]

        SQTOL = (SQTOL**0.5/(NSYS**0.5))

        if SQTOL <= TOL and TMPTOL <= TOL:
            break


