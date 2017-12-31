#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import sys
import numpy as np
#from numpy import random_intel
import lattice_map as lmap
import algo
import tool

# seeding the random number generators
np.random.seed(0)

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

    ## creating the nearest and next nearest neighbour map for ##
    ## the triangular lattice                                  ##

    nbr = np.zeros((NSYS,ZCO),dtype=np.int64)
    lmap.lattice_map(LSYS,nbr)
    
    # field params
    M = np.zeros((2,NSYS,3),dtype=np.float64)
    M0 = np.zeros((2,NSYS,3),dtype=np.float64)

    # the hamiltonian matrix and coupling
    J = np.zeros((NSYS,NSYS),dtype=np.float64)
    
    # relative angle between pure and impure system
    theta = np.zeros((NSYS,NSYS),dtype=np.float64)

    # iteration loop for disorder
    for i in range(0,ITERDISD):
        # looping over bootstrapped initialisations
        for g in range(0,BOOTNUM): 

            # fixing the coupling matrix
            tool.init_cpl(NSYS,nbr,
                    ZCO,DELTA,ALPHA,J)

            # initiating parameters
            tool.init_param(NSYS,nbr,
                    ZCO,
                    DELTA,ALPHA,
                    ANGVAR,M)
            
            # classical algorithm to get the magnetic
            # ground state of impurity system

            # introducing sharp anomalous dipole density
            # wall
            
            HSYS = LSYS*(LSYS+1)/2
            for site in range(0,HSYS):
                for p in [3]:
                    J[site,nbr[site,p]] = 0.25
                    J[nbr[site,p],site] = 0.25
            for site in range(HSYS,NSYS):
                for p in [1]:
                    J[site,nbr[site,p]] = 0.75
                    J[nbr[site,p],site] = 0.75

            algo.classic_zmc(CLNUM,NSYS,nbr,J,M,GTOL)

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
                    theta=theta,
                    spin0=M0[1],
                    spin=M[1],
                    )

