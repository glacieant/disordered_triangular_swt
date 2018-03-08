#!/usr/bin/env python

### This is the subroutine to calculate the various mean ###
### field parameters from the eigenvectors of the mean   ###
### hamiltonian

import numpy as np
import tool


def paramset(NSYS,nbr,UPAR,TEMP,
        eigen,vec,bond,M,lmult):

    # the fermi function array for temperature T
    farray = tool.fermi(eigen,TEMP)

    vec = np.einsum('ij,j->ij',vec,farray)
    
    # iteration over the lattice sites
    for i in range(0,NSYS):

        # getting the flavour number
        fdagf = (np.vdot(vec[i,:],vec[i,:]) + 
                np.vdot(vec[i+NSYS,:],vec[i+NSYS,:]) - 1.0).real

        # setting up the lagrange multiplier
        # lx = np.random_intel.uniform(10.0,20.0)
        lx = 5
        lmult[1,i] -= fdagf/lx

        # the flavour variable matrix
        ada = np.vdot(vec[i,:],vec[i,:])
        adb = np.vdot(vec[i,:],vec[i+NSYS,:])
        bda = adb.conj()
        bdb = np.vdot(vec[i+NSYS,:],vec[i+NSYS,:])
        
        # getting the magnetic moments
        MAG = 0.5*np.array([adb+bda,
            (-adb+bda)*1.0j,ada-bdb]).real
        rx = np.random_intel.uniform(UPAR[0][0],
                UPAR[0][1])
        rx = 0.5
        M[1,i] = ((1.0-rx**2)**0.5)*M[1,i]+ rx*MAG

    tx = np.random_intel.uniform(UPAR[1][0],UPAR[1][1])
    #tx  = 0.5
    for i in range(0,NSYS):

        """
        # iteration over the neighbours
        j = nbr[i,0]

        xbond = (np.vdot(vec[i,:],vec[j,:])
                +np.vdot(vec[i+NSYS,:],vec[j+NSYS,:]))
        xbond = 1.0
        #tx = np.random_intel.uniform(UPAR[1][0],UPAR[1][1])
        parabond = (1.0-tx)*bond[1,i,j] + tx*xbond
        bond[1,i,j] = parabond
        bond[1,j,i] = bond[1,i,j].conj()
        """
        for m in range(0,12):

            j = nbr[i,m]

            if j > i:

                # computing the bond operators
                xbond = (np.vdot(vec[i,:],farray*vec[j,:])
                        +np.vdot(vec[i+NSYS,:],
                            farray*vec[j+NSYS,:]))
                #tx = np.random_intel.uniform(UPAR[1][0],
                #        UPAR[1][1])
                #tx = 0.1
                bond[1,i,j] = (1.0-tx)*bond[1,i,j] + tx*xbond
                #bond[1,i,j] = xbond
                bond[1,j,i] = bond[1,i,j].conj()
