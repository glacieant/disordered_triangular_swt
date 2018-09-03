#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import sys
import numpy as np
import lattice_map as lmap
import algo
import tool

# defining the main function block

def main(const):

    ## setting the calculation parameters ##

    LSYS = const.LSYS # system length
    NSYS = LSYS*LSYS # system size
    ZCO = const.ZCO # the count of nearest and next nearest neighbours
    ITERDISD = const.ITERDISD # disorder iteration number
    ANGVAR = const.ANGVAR # fluctuation of inital angle
    BOOTNUM = const.BOOTNUM # bootstrap sample number
    DELTA = const.DELTA # magntitude of disorder
    ALPHA = const.ALPHA # magnititude of J'/J
    GTOL = const.GTOL # a global tolerance for the code
    CLNUM = const.CLNUM # the maximum number of iteration for the classical routine
    DNMR = const.DNMR # denominator for batch execution

    # seeding the random number generators
    # np.random.seed([LSYS,DNMR])

    ## creating the nearest and next nearest neighbour map for ##
    ## the triangular lattice                                  ##

    nbr = np.zeros((NSYS,ZCO),dtype=np.int64)
    lmap.lattice_map(LSYS,nbr)
    
    # field params
    M = np.zeros((2,NSYS,2),dtype=np.float64)
    M0 = np.zeros((2,NSYS,2),dtype=np.float64)
    
    # the hamiltonian matrix and coupling
    J = np.zeros((NSYS,NSYS),dtype=np.float64)
    
    # iteration loop for disorder
    for i in range(0,ITERDISD):

        # fixing the coupling matrix
        tool.init_cpl(NSYS,nbr,
                    ZCO,DELTA,ALPHA,J)
        
        ENMIN = -1000000.0
        # looping over bootstrapped initialisations
        for g in range(0,BOOTNUM):

            # initiating parameters
            if g%2 == 0:
                tool.init_param_ord(NSYS,nbr,
                    ZCO,
                    DELTA,ALPHA,
                    ANGVAR,M0)
            else:
                tool.init_param_rand(NSYS,nbr,
                        ZCO,
                        DELTA,ALPHA,
                        ANGVAR,M0)
           
            # classical algorithm to get the magnetic
            # ground state of impurity system

            algo.classic_zmc(CLNUM,NSYS,nbr,J,M0,GTOL)

            ENTEMP = tool.en_calc(NSYS,nbr,J,M0[1])

            if ENTEMP > ENMIN:

                ENMIN = ENTEMP

                np.copyto(M,M0)


        ## output data ##

        # saving the final output data
        np.savez_compressed("out/data/FNL"+
                "_L_"+str(LSYS)+
                "_DLT_"+str(DELTA)+
                "_ALP_"+str(ALPHA)+
                "_DISD_"+str(i)+
                "_BOOT_"+str(0)+
                "_DNM_"+str(DNMR)+
                ".npz",
                J=J,
                spin=M[1],
                )

