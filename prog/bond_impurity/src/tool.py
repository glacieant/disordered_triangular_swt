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

    if DELTA != 0.0:
        JBANK = np.random.lognormal(1.0,DELTA,size=(NSYS*12))
    else:
        JBANK = np.full((NSYS*12),1.0)

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
        ZCO,
        #XPAR,
        DELTA,ALPHA,
        ANGVAR,
        #BONDVAR,
        #J,
        #bond,chi,
        M
        #,B
        #,fnum,lmult
        ):
    
    # zeroing out the parameter matrices
    
    #chi[:] = 0.0
    #bond[:] = 0.0
    M[:] = 0.0
    #B[:] = 0.0
    #lmult[:] = 0.0
    #fnum[:] = 0.0

    """
    # setting up the bond matrices
    for i in range(0,NSYS):

        for j in nbr[i]:

            if j > i:

                bond[1,i,j] = (np.random.uniform(-BONDVAR,BONDVAR)
                        + np.random.uniform(-BONDVAR,BONDVAR)*1.0j)
                bond[1,j,i] = bond[1,i,j].conj()
                chi[i,j] = 0.5*(1.0-XPAR)*J[i,j]*bond[1,i,j]
                chi[j,i] = chi[i,j].conj()
    """

    #M[1] = np.random.uniform(-ANGVAR,ANGVAR,size=(NSYS,3))

    
    # setting up the spin and field matrices
    ROT = np.array([[np.cos(2.0*np.pi/3.0),-np.sin(2.0*np.pi/3.0),0.0],
        [np.sin(2.0*np.pi/3.0),np.cos(2.0*np.pi/3.0),0.0],
        [0.0,0.0,1.0]])

    LSYS = int(NSYS**0.5)
    for i in range(0,NSYS):
    
        if i==0:
            fleet=0.0
            M[0,i] = np.array([np.cos(fleet),np.sin(fleet),0.0])
            if ANGVAR != 0.0:
                phi = np.random_intel.uniform(-ANGVAR,ANGVAR)
            else:
                phi =0.0
            M[1,i] = np.einsum('ab,b->a',[[np.cos(phi),-np.sin(phi),0.0],
                [np.sin(phi),np.cos(phi),0.0],
                [0.0,0.0,1.0]],
                M[0,i])
        elif i%LSYS!=0:
            M[0,i] = np.einsum('ab,b->a',ROT,M[0,i-1])
            if ANGVAR != 0.0:
                phi = np.random_intel.uniform(-ANGVAR,ANGVAR)
            else:
                phi =0.0
            M[1,i] = np.einsum('ab,b->a',[[np.cos(phi),-np.sin(phi),0.0],
                [np.sin(phi),np.cos(phi),0.0],
                [0.0,0.0,1.0]],
                M[0,i])
        elif i!=0 and i%LSYS==0:
            M[0,i] = M[0,i-1]
            if ANGVAR != 0.0:
                phi = np.random_intel.uniform(-ANGVAR,ANGVAR)
            else:
                phi =0.0
            M[1,i] = np.einsum('ab,b->a',[[np.cos(phi),-np.sin(phi),0.0],
                [np.sin(phi),np.cos(phi),0.0],
                [0.0,0.0,1.0]],
                M[0,i])

    """
    for i in range(0,NSYS):
        
        for j in nbr[i]:
            
            B[i] += (0.5*
                    #XPAR*
                    J[i,j]*M[1,j])

    """
