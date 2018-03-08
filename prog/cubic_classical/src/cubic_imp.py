#!/usr/bin/env python

### This is the main program for the classical energy ###
### minimization in a J1-J4 cubic lattice.            ###

import sys
import numpy as np
import algo
import tool

# seeding the random number generators
np.random.seed(0)

# defining the main function block

def main(const):

    ## setting the calculation parameters ##

    LSYS = const.LSYS # system length
    NSYS = LSYS**3 # system size
    ITERDISD = const.ITERDISD # disorder iteration number
    ANGVAR = const.ANGVAR # fluctuation of inital angle
    BOOTNUM = const.BOOTNUM # bootstrap sample number
    DELTA = const.DELTA # magntitude of disorder
    ALPHA = const.ALPHA # magnititude of J'/J
    GTOL = const.GTOL # a global tolerance for the code
    CLNUM = const.CLNUM # the maximum number of iteration for the classical routine
    DNMR = const.DNMR # denominator for batch execution

    # field params
    M = np.zeros((2,NSYS,3),dtype=np.float64)
    M0 = np.zeros((2,NSYS,3),dtype=np.float64)

    # the hamiltonian matrix and coupling
    J = np.zeros((NSYS,NSYS),dtype=np.float64)
    
    # iteration loop for disorder
    for i in range(0,ITERDISD):

        # fixing the coupling matrix
        tool.init_cpl(LSYS,DELTA,ALPHA,J)
        
        # looping over bootstrapped initialisations
        for g in range(0,BOOTNUM): 

            # initiating parameters
            tool.init_param(LSYS,DELTA,ALPHA,ANGVAR,M)
            tool.init_param(LSYS,DELTA,ALPHA,ANGVAR,M0)
            
            # classical algorithm to get the magnetic
            # ground state of the pure and the bond 
            # impurity system

            algo.classic_zmc(CLNUM,LSYS,J,M0,GTOL)

            """

            M[1,:] = M0[1,:]
            
            impsite = (((LSYS+1)/2)*(LSYS**2) 
                    + (LSYS*((LSYS+1)/2-1)+LSYS/2))

            J[impsite-1,impsite] = 0.0
            J[impsite,impsite-1] = 0.0
            """

            algo.classic_zmc(CLNUM,LSYS,J,M,GTOL)

            ## output data ##

            # saving the final output data
            np.savez_compressed("out/data/FNL"+
                    "_L_"+str(LSYS)+
                    "_DLT_"+str(DELTA)+
                    "_ALP_"+str(ALPHA)+
                    "_DISD_"+str(i)+
                    "_BOOT_"+str(g)+
                    ".npz",
                    J=J,
                    spin0=M0[1],
                    spin=M[1],
                    )

