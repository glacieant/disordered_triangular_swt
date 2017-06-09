#!/usr/bin/env python

### program to generate random couplings ###

import numpy as np
from numpy import random_intel

def coupling(NSYS,ZCO,nbr,
        DELTA,ALPHA,J):

    # randomised couplings
    if DELTA != 0.0:
        JSTAK =  np.random_intel.uniform(
                1.0-DELTA,1.0+DELTA,size=(NSYS,ZCO))
    else:
        JSTAK = np.full((NSYS,ZCO),1.0)

    # scaling the next nearest neigbour
    # couplings compared to the original
    # couplings
    JSTAK[:,ZCO/2:ZCO] *= ALPHA

    # loop over the lattice sites to
    # to assign random couplings
    for i in range(0,NSYS):

        # loop over neighbouring sites
        for k in range(0,ZCO):

            j = nbr[i,k]
            J[i][j] = JSTAK[i][k]

    J = (J + np.transpose(J))/2.0

