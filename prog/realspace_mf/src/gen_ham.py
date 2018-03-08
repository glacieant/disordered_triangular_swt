#!/usr/bin/env python

### This is the program to create the meanfield hamiltonian ###
### from the mean-field parameters ###
import numpy as np

def mfmatrix(NSYS,XPAR,nbr,J,bond,M,lmult,ham):

    ham[:] = 0.0

    chi = (1-XPAR)*J*bond[1]/2

    ham[:NSYS,:NSYS] = ham[NSYS:,NSYS:] = -chi.conj()

    for i in range(0,NSYS):

        H = np.array([0.0,0.0,0.0])

        for k in nbr[i]:
            H += XPAR*J[i,k]*M[1,k]

        # the contribution from the local fields
        
        hdotsigma = 0.5*np.array([[H[2],H[0]-H[1]*1.0j],
            [H[0]+H[1]*1.0j,-H[2]]])

        ham[i,i] += hdotsigma[0,0]
        ham[i,i+NSYS] += hdotsigma[0,1]
        ham[i+NSYS,i] += hdotsigma[1,0]
        ham[i+NSYS,i+NSYS] += hdotsigma[1,1]

        # the contribution from the lagrange multiplier

        ham[i,i] -= lmult[1,i]
        ham[i+NSYS,i+NSYS] -= lmult[1,i]


