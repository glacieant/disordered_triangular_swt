#!/usr/bin/env python

### This is an assortment of programs, ###
### useful for the main computation    ###

import numpy as np

# the finite temperature fermi function
def fermi(E,T):

    # with a factor of 2 for degenracy
    return 1.0/(np.exp(E/T)+1.0)

def init_cpl(LSYS,DELTA,ALPHA,J):

    NSYS = LSYS**3
    
    # zeroing out the couplings
    J[:] = 0.0
   
    if DELTA != 0.0:
        JBANK = np.random.uniform(1.0-DELTA,1.0+DELTA,size=(NSYS*12))
    else:
        JBANK = np.full((NSYS*12),1.0)
    
    x = 0
    for i in range(0,LSYS):

        for j in range(0,LSYS):
        
            for k in range(0,LSYS):

                p = i + j*LSYS +k*(LSYS**2)

                for NN in [-1,1]:
                    if (i+NN) >= 0 and (i+NN) < LSYS:
                        q = (i + NN) + j*LSYS + k*(LSYS**2)
                    elif (i+NN) >= 0 and (i+NN) >= LSYS:
                        q = (i - LSYS + NN) + j*LSYS + k*(LSYS**2)
                    elif (i+NN) < 0:
                        q = (i + LSYS + NN) + j*LSYS + k*(LSYS**2)

                    if p < q:
                        J[p,q] = JBANK[x]
                        J[q,p] = J[p,q]
                        x += 1

                    if (j+NN) >= 0 and (j+NN) < LSYS:
                        q = i + (j + NN)*LSYS + k*(LSYS**2)
                    elif (j+NN) >= 0 and (j+NN) >= LSYS:
                        q = i + (j - LSYS + NN)*LSYS + k*(LSYS**2)
                    elif (j+NN) < 0:
                        q = i+ (j + LSYS + NN)*LSYS + k*(LSYS**2)
                    
                    if p < q:
                        J[p,q] = JBANK[x]
                        J[q,p] = J[p,q]
                        x += 1

                    if (k+NN) >= 0 and (k+NN) < LSYS:
                        q = i + j*LSYS +  (k + NN)*(LSYS**2)
                    elif (k+NN) >= 0 and (k+NN) >= LSYS:
                        q = i + j*LSYS + (k - LSYS + NN)*(LSYS**2)
                    elif (k+NN) < 0:
                        q = i+ j*LSYS + (k + LSYS + NN)*(LSYS**2)
                    
                    if p < q:
                        J[p,q] = JBANK[x]
                        J[q,p] = J[p,q]
                        x += 1

                for NN in [-2,2]:
                    if (i+NN) >= 0 and (i+NN) < LSYS:
                        q = (i + NN) + j*LSYS +k*(LSYS**2)
                    elif (i+NN) >= 0 and (i+NN) >= LSYS:
                        q = (i - LSYS + NN) + j*LSYS +k*(LSYS**2)
                    elif (i+NN) < 0:
                        q = (i + LSYS + NN) + j*LSYS +k*(LSYS**2)
                    
                    if p < q:
                        J[p,q] = ALPHA*JBANK[x]
                        J[q,p] = J[p,q]
                        x += 1

                    if (j+NN) >= 0 and (j+NN) < LSYS:
                        q = i + (j + NN)*LSYS +k*(LSYS**2)
                    elif (j+NN) >= 0 and (j+NN) >= LSYS:
                        q = i + (j - LSYS + NN)*LSYS +k*(LSYS**2)
                    elif (j+NN) < 0:
                        q = i+ (j + LSYS + NN)*LSYS +k*(LSYS**2)
                    
                    if p < q:
                        J[p,q] = ALPHA*JBANK[x]
                        J[q,p] = J[p,q]
                        x += 1

                    if (k+NN) >= 0 and (k+NN) < LSYS:
                        q = i + j*LSYS +  (k + NN)*(LSYS**2)
                    elif (k+NN) >= 0 and (k+NN) >= LSYS:
                        q = i + j*LSYS + (k - LSYS + NN)*(LSYS**2)
                    elif (k+NN) < 0:
                        q = i+ j*LSYS + (k + LSYS + NN)*(LSYS**2)
                    
                    if p < q:
                        J[p,q] = ALPHA*JBANK[x]
                        J[q,p] = J[p,q]
                        x += 1

def init_param(LSYS,DELTA,ALPHA,ANGVAR,M):
    
    # zeroing out the parameter matrices
    
    M[:] = 0.0

    NSYS = LSYS**3
   
    """
    #M[1] = np.random.uniform(-ANGVAR,ANGVAR,size=(NSYS,3))
    #M[1] = np.einsum('i,ij->ij',np.einsum('ij,ij->i',M[1],M[1])**(-0.5),M[1])

    
    for i in range(0,NSYS):

        #M[1,i] = np.array([1.0,0.0,0.0])
        RND = np.random.uniform(0.0,ANGVAR)
        RX = i+(i/LSYS)%2 + (i/(LSYS**2))%2
        M[1,i] = ((-1.0)**RX)*np.array([np.sqrt(1.0-RND**2),RND,0.0])
    """
    # setting up the spin and field matrices

    M[1,0] = np.array([1.0,0.0,0.0])
    for i in range(1,NSYS):
        
        if (i%LSYS == 0):
            if ((i%(LSYS*LSYS))==0.0):
                phx = 2.0*np.pi/3.0
            else:
                phx = -2.0*np.pi/3.0
        else:
            phx = 2.0*np.pi/3.0

        ROT = np.array([[np.cos(phx),-np.sin(phx),0],
            [np.sin(phx),np.cos(phx),0],[0,0,1]])
        M[1,i] = np.einsum('ab,b',ROT,M[1,i-1])
