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
        JBANK = np.random.uniform(1.0-DELTA,1.0+DELTA,size=(NSYS*12))
    else:
        JBANK = np.full((NSYS*6),1.0)
    q = 0
    for i in range(0,NSYS):
        for m in range(0,6):
            j = nbr[i,m]
            if j > i:
                J[i,j] = JBANK[q]
                J[j,i] = J[i,j]
                q += 1
        for m in range(6,12):
            j = nbr[i,m]
            if j > i:
                J[i,j] = ALPHA*JBANK[q]
                J[j,i] = J[i,j]
                q += 1

def en_calc(NSYS,nbr,J,M):

    EN = 0.0

    for i in range(0,NSYS):

        for j in nbr[i]:

            EN += 0.5*J[i,j]*np.dot(M[i],M[j])

    EN = EN/(3*NSYS)

    return EN

   
def init_param_rand(NSYS,nbr,
        ZCO,
        DELTA,ALPHA,
        ANGVAR,M):
    
    # zeroing out the parameter matrices
    
    
    M[:] = 0.0
   
    
    M[1] = np.random.uniform(-ANGVAR,ANGVAR,size=(NSYS,2))
    M[1] = np.einsum('i,ij->ij',np.einsum('ij,ij->i',M[1],M[1])**(-0.5),M[1])
    
def init_param_ord(NSYS,nbr,
        ZCO,
        DELTA,ALPHA,
        ANGVAR,M):
    
    # zeroing out the parameter matrices
    
    M[:] = 0.0
    
    # setting up the spin and field matrices
    ROT = np.array([[np.cos(2.0*np.pi/3.0),-np.sin(2.0*np.pi/3.0)],
        [np.sin(2.0*np.pi/3.0),np.cos(2.0*np.pi/3.0)]])

    LSYS = int(NSYS**0.5)

    VAR = ANGVAR

    for i in range(0,NSYS):
    
        if i==0:
            fleet=0.0
            M[0,i] = np.array([np.cos(fleet),np.sin(fleet)])
            if VAR != 0.0:
                ph = VAR*np.random.uniform(-1.0,1.0)
            else:
                ph = 0.0
            M[1,i] = np.einsum('ab,b->a',[
                [np.cos(ph),-np.sin(ph)],
                [np.sin(ph),np.cos(ph)]
                ],
                M[0,i])
        elif i!=0 and i%LSYS!=0:
            M[0,i] = np.einsum('ab,b->a',ROT,M[0,i-1])
            if VAR != 0.0:
                ph = VAR*np.random.uniform(-1.0,1.0)
            else:
                ph = 0.0
            M[1,i] = np.einsum('ab,b->a',[
                [np.cos(ph),-np.sin(ph)],
                [np.sin(ph),np.cos(ph)]
                ],
                M[0,i])
        elif i!=0 and i%LSYS==0:
            M[0,i] = M[0,i-1]
            if VAR != 0.0:
                ph = VAR*np.random.uniform(-1.0,1.0)
            else:
                ph = 0.0
            M[1,i] = np.einsum('ab,b->a',[
                [np.cos(ph),-np.sin(ph)],
                [np.sin(ph),np.cos(ph)]
                ],
                M[0,i])

