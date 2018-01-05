#!/usr/bin/env python

### This is an assortment of programs, ###
### useful for the main computation    ###

import numpy as np

# the finite temperature fermi function
def fermi(E,T):

    # with a factor of 2 for degenracy
    return 1.0/(np.exp(E/T)+1.0)

def init_cpl(NSYS,nbr,
        ZCO,XPAR,
        DELTA,ALPHA,J):

    # zeroing out the couplings

    J[:] = 0.0
    
    if DELTA != 0.0:
        JBANK = np.random.lognormal(1.0,DELTA,size=(NSYS*12))
    else:
        JBANK = np.full((NSYS*12),1.0)
    
    #JBANK = np.full((NSYS*12),1.0)
    q = 0
    for i in range(0,NSYS):
        for m in range(0,6):
            j = nbr[i,m]
            if j > i:
                J[i][j] = JBANK[q]
                J[j][i] = J[i][j]
                q += 1
        for m in range(6,12):
            j = nbr[i,m]
            if j > i:
                J[i][j] = ALPHA*JBANK[q]
                J[j][i] = J[i][j]
                q += 1


def init_param(NSYS,nbr,
        ZCO,XPAR,
        DELTA,ALPHA,
        ANGVAR,BONDVAR,
        J,bond,chi,
        M,B,fnum,lmult):
    
    # zeroing out the parameter matrices
    
    chi[:] = 0.0
    bond[:] = 0.0
    M[:] = 0.0
    B[:] = 0.0
    lmult[:] = 0.0
    fnum[:] = 0.0
    
    # setting up the bond matrices
    for i in range(0,NSYS):

        for j in nbr[i]:

            if j > i:

                bond[1,i,j] = (np.random.uniform(-BONDVAR,BONDVAR)
                        + np.random.uniform(-BONDVAR,BONDVAR)*1.0j)
                bond[1,j,i] = bond[1,i,j].conj()
                chi[i,j] = 0.5*(1.0-XPAR)*J[i,j]*bond[1,i,j]
                chi[j,i] = chi[i,j].conj()
    
    M[1] = np.random.uniform(-ANGVAR,ANGVAR,size=(NSYS,3))
    #M[1] = SPIN/np.einsum("ij,ij->i",SPIN,SPIN)**0.5
    
    for i in range(0,NSYS):
        
        for j in nbr[i]:
            
            B[i] += 0.5*XPAR*J[i,j]*M[1,j]

