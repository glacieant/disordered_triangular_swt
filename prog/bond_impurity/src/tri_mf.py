#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import sys
import numpy as np
from numpy import random_intel
import lattice_map as lmap
import algo
import gen_ham
import calc_mfparam as mfp
import tool

# seeding the random number generators
np.random_intel.seed(0)

# defining the main function block

def main(const):

    ## setting the calculation parameters ##

    LSYS = const.LSYS # system length
    NSYS = LSYS*LSYS # system size
    MDIM = NSYS*2 # dimension of mean-field matrix
    ZCO = const.ZCO # the count of nearest and next nearest neighbours
    ITERDISD = const.ITERDISD # disorder iteration number
    ANGVAR = const.ANGVAR # fluctuation of inital angle
    BONDVAR = const.BONDVAR # fluctuation of initial bonds
    BOOTNUM = const.BOOTNUM # bootstrap sample number
    SIMSZE = const.SIMSZE # simulation iteration number
    DELTA = const.DELTA # magntitude of disorder
    ALPHA = const.ALPHA # magnititude of J'/J
    MFTOL = const.MFTOL # tolerance for the mean-field calculation
    GTOL = const.GTOL # a global tolerance for the code
    TEMP = const.TEMP # the temperature for modulating the fermi function
    XPAR = const.XPAR # the fudge parameter for decoupling
    UPAR = const.UPAR # the update ration for iterations
    CLNUM = const.CLNUM # the maximum number of iteration for the classical routine
    DNMR = const.DNMR # denominator for batch execution

    ## creating the nearest and next nearest neighbour map for ##
    ## the triangular lattice                                  ##

    nbr = np.zeros((NSYS,ZCO),dtype=np.int64)
    lmap.lattice_map(LSYS,nbr)

    # initiating empty parameter arrays of fixed shape

    # bond params
    bond = np.zeros((2,NSYS,NSYS),dtype=np.complex128)
    chi = np.zeros((NSYS,NSYS),dtype=np.complex128)
    
    # field params
    M = np.zeros((2,NSYS,3),dtype=np.float64)
    B = np.zeros((NSYS,3),dtype=np.float64)

    # occupancy and chemical potential
    fnum = np.zeros(NSYS,dtype=np.float64)
    lmult = np.zeros(NSYS,dtype=np.float64)

    # the hamiltonian matrix and coupling
    J = np.zeros((NSYS,NSYS),dtype=np.float64)
    ham = np.zeros((MDIM,MDIM),dtype=np.complex128)

    # iteration loop for disorder
    for i in range(0,ITERDISD):

        # fixing the coupling matrix
        tool.init_cpl(NSYS,nbr,
                ZCO,XPAR,
                DELTA,ALPHA,J)

        # looping over bootstrapped initialisations
        for g in range(0,BOOTNUM): 

            # initiating parameters
            tool.init_param(NSYS,nbr,
                    ZCO,XPAR,
                    DELTA,ALPHA,
                    ANGVAR,BONDVAR,
                    J,bond,chi,
                    M,B,fnum,lmult)
            

            # classical algorithm to get the magnetic
            # ground state
            # algo.classic_zmc(CLNUM,NSYS,nbr,J,M,GTOL)

            # rubbish tolerance
            TOL = 100.0
            # rubbish reference energies
            EN = 100.0
            EN_M = 0.0
            # iterations for mean-field solution convergence
            # with number conservation for fermions
            ###
            for j in range(0,SIMSZE):

                # checking the tolerance
                if TOL > MFTOL:

                    M[0] = M[1].copy()
                    bond[0] = bond[1].copy()

                    # creating the mean-field hamiltonian
                    gen_ham.mfmatrix(NSYS,chi,
                            B,lmult,ham)

                    # diagonalising the hamiltonian
                    en_mf, v_mf = np.linalg.eigh(ham)

                    # calculating the mean-field parameters
                    mfp.paramset(NSYS,j,nbr,
                            XPAR,UPAR,TEMP,
                            J,en_mf,v_mf,
                            bond,chi,
                            M,B,
                            fnum,lmult)
                    
                    # differences in the update step
                    # norm difference
                    # individual difference
                    TOL = np.amax([np.linalg.norm(bond[1]-bond[0]),
                        np.linalg.norm(M[1]-M[0]),
                        np.linalg.norm(fnum[1])])/LSYS
                    
                else:
                    break



            ## output data ##

            # saving the final output data
            np.savez_compressed("out/data/FNL"+
                    "_L_"+str(LSYS)+
                    "_DLT_"+str(DELTA)+
                    "_ALP_"+str(ALPHA)+
                    "_DISD_"+str(i)+
                    "_BOOT_"+str(g)+
                    ".npz",
                    T=TEMP,
                    J=J,
                    inum=j,
                    err=TOL,
                    ensys=en_mf,
                    vsys=v_mf,
                    bond=bond[1],
                    spin=M[1],
                    lmult=lmult)

