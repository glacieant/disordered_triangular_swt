#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import sys
import numpy as np
import gen_ham
import randgen
import lattice_map as lmap

# seeding the random number generators
np.random_intel.seed(0)

# defining the main function block

def main(const):

    ## setting the calculation parameters ##

    LSYS = const.LSYS # system length
    NSYS = LSYS*LSYS # system size
    ZCO = const.ZCO # the count of nearest and next nearest neighbours
    ITERDISD = const.ITERDISD # disorder iteration number
    HFIELD = const.HFIELD # magnitude of the external field
    DELTA = const.DELTA # disorder magnitude array
    DSIZE = len(DELTA)
    ALPHA = const.ALPHA # J'/J ratio array
    ASIZE = len(ALPHA)
    GTOL = const.GTOL # a global tolerance for the code
    S = const.S # spin magnitude value
    DNMR = const.DNMR # denominator for batch execution
    MINF = const.MINF # minimum frequency for structure factor
    MAXF = const.MAXF # maximum frequency for structure factor
    FGRID = const.FGRID # frequency grid for structure factor
    GWDTH = const.GWDTH # width of the delta function 

    ## creating the nearest and next nearest neighbour map for ##
    ## the triangular lattice                                  ##

    nbr = np.zeros((NSYS,ZCO),dtype=np.int)
    lmap.design(LSYS,nbr)

    # initiating empty parameter arrays of fixed shape

    # hamiltonian,coupling, eigenvalues and eigenvectors
    ham = np.zeros((NSYS,NSYS),dtype=np.float64)
    J = np.zeros((NSYS,NSYS),dtype=np.float64)
    eigvec = np.zeros((NSYS,NSYS),dtype=np.float64)
    eigval = np.zeros(NSYS,dtype=np.float64)

    # dynamic structure factor
    evq = np.zeros((NSYS,LSYS,LSYS),dtype=np.complex)
    dynstrc = np.zeros((DSIZE,ASIZE,FGRID,LSYS,LSYS),dtype=np.float64)
    dynstrc_var = np.zeros((DSIZE,ASIZE,FGRID,LSYS,LSYS),dtype=np.float64)

    # outer loop over disorder and NNN coupling strength
    for DNUM in range(0,DSIZE):

        for ANUM in range(0,ASIZE):
            
            
            # iteration loop for disorder
            for i in range(0,ITERDISD):

                # randomised couplings
                randgen.coupling(NSYS,ZCO,
                        nbr,DELTA[DNUM],
                        ALPHA[ANUM],J)


                # hamiltonian matrix
                gen_ham.create(NSYS,ZCO,S,
                        nbr,HFIELD,J,ham)

                # diagonalising the matrix
                eigval, eigvec = np.linalg.eigh(ham)

                # fourier transformed vectors
                for j in range(0,NSYS):

                    # fourier transforming the unitary matrix
                    evq[j] = np.fft.fft2(np.reshape(
                        eigvec[:,j],(LSYS,LSYS)),norm='ortho')

                # loop over eigenvalues
                for j in range(0,NSYS):
                    
                    # loop over frequencies
                    for k in range(0,FGRID):
                        
                        # guassian approx. of delta fn.
                        FREQ = MINF + ((MAXF-MINF)*k)/(FGRID-1)
                        DWDTH = np.exp((-(FREQ-eigval[j])**2)/
                                GWDTH)/(GWDTH*np.pi)**0.5
                        strc = S*(DWDTH*evq[j]*evq[j].conj()).real
                        
                        dynstrc[DNUM,ANUM,k] += strc
                        dynstrc_var[DNUM,ANUM,k] += strc**2



        # the final output
        np.savez_compressed("out/strfc.npz",
                LSYS=LSYS,
                DELTA=DELTA,
                ALPHA=ALPHA,
                ITERDISD=ITERDISD,
                MINF=MINF,
                MAXF=MAXF,
                FGRID=FGRID,
                eigval=eigval,
                eigvec=eigvec,
                dynstrc=dynstrc,
                dynstrc_var=dynstrc_var)





