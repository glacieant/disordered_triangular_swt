#!/usr/bin/env python

### This is an assortment of programs, ###
### useful for the main computation    ###

import numpy as np

# the finite temperature fermi function
def fermi(E,T):

    # with a factor of 2 for degenracy
    return 1.0/(np.exp(E/T)+1.0)

def init_cpl(NSYS,nbr,
        ZCO,DELTA,
        ALPHA,J):

    # zeroing out the couplings

    J[:] = 0.0
    
    for i in range(0,NSYS):
        for k in range(0,6):
            j = nbr[i,k]
            if j > i:
                J[i,j] = 1.0
                J[j,i] = 1.0


def init_param(NSYS,nbr,
        ZCO,
        DELTA,ALPHA,
        ANGVAR,M):
    
    # zeroing out the parameter matrices
    
    
    M[:] = 0.0
    
    M[1] = np.random.uniform(-ANGVAR,ANGVAR,size=(NSYS,3))
    M[1] = np.einsum('i,ij->ij',np.einsum('ij,ij->i',M[1],M[1])**(-0.5),M[1])

    """
    # setting up the spin and field matrices
    ROT = np.array([[np.cos(2.0*np.pi/3.0),-np.sin(2.0*np.pi/3.0),0],
        [np.sin(2.0*np.pi/3.0),np.cos(2.0*np.pi/3.0),0],[0,0,1]])

    LSYS = int(NSYS**0.5)
    for i in range(0,NSYS):
    
        if i==0:
            fleet=0.0
            M[0,i] = np.array([np.cos(fleet),np.sin(fleet),0])
            if ANGVAR != 0.0:
                phi = np.random_intel.uniform(-ANGVAR,ANGVAR)
            else:
                phi =0.0
            M[1,i] = np.einsum('ab,b->a',[[np.cos(phi),-np.sin(phi),0],
                [np.sin(phi),np.cos(phi),0],[0,0,1]],
                M[0,i])
        elif i!=0 and i%LSYS!=0:
            M[0,i] = np.einsum('ab,b->a',ROT,M[0,i-1])
            if ANGVAR != 0.0:
                phi = np.random_intel.uniform(-ANGVAR,ANGVAR)
            else:
                phi =0.0
            M[1,i] = np.einsum('ab,b->a',[[np.cos(phi),-np.sin(phi),0],
                [np.sin(phi),np.cos(phi),0],[0,0,1]],
                M[0,i])
        elif i!=0 and i%LSYS==0:
            M[0,i] = M[0,i-1]
            if ANGVAR != 0.0:
                phi = np.random_intel.uniform(-ANGVAR,ANGVAR)
            else:
                phi =0.0
            M[1,i] = np.einsum('ab,b->a',[[np.cos(phi),-np.sin(phi),0],
                [np.sin(phi),np.cos(phi),0],[0,0,1]],
                M[0,i])


    """
