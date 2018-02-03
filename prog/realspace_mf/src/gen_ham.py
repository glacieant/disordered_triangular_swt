#!/usr/bin/env python

### This is the program to create the meanfield hamiltonian ###
### from the mean-field parameters ###
import numpy as np

def mfmatrix(NSYS,chi,B,lmult,ham):

    ham[:] = 0.0+0.0j

    ham[:NSYS,:NSYS] = ham[NSYS:,NSYS:] = -chi.conj()

    for i in range(0,NSYS):

        # the contribution from the local fields
        
        offblock = 0.5*np.array([[B[i][2],B[i][0]-B[i][1]*1.0j],
            [B[i][0]+B[i][1]*1.0j,-B[i][2]]])

        ham[i,i] += offblock[0,0]
        ham[i,i+NSYS] += offblock[0,1]
        ham[i+NSYS,i] += offblock[1,0]
        ham[i+NSYS,i+NSYS] += offblock[1,1]

        # the contribution from the lagrange multiplier

        ham[i,i] += lmult[i]
        ham[i+NSYS,i+NSYS] += lmult[i]


