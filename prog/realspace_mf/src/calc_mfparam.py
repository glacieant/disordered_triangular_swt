#!/usr/bin/env python

### This is the subroutine to calculate the various mean ###
### field parameters from the eigenvectors of the mean   ###
### hamiltonian

import numpy as np
import tool


def paramset(NSYS,INUM,nbr,XPAR,
        UPAR,T,J,eigen,vec,
        bond,chi,
        M,B,
        fnum,lmult):

    # the fermi function array for temperature T
    farray = tool.fermi(eigen,T)

    # iteration over the lattice sites
    for i in range(0,NSYS):

        # getting the flavour number
        fnum[i] = (np.vdot(vec[i,:],
            farray*vec[i,:]) + 
            np.vdot(vec[i+NSYS,:],
                farray*vec[i+NSYS,:]) - 1.0).real

        # setting up the lagrange multiplier
        lmult[i] += fnum[i]/20.0

        # the flavour variable matrix
        ada = np.vdot(vec[i,:],farray*vec[i,:])
        adb = np.vdot(vec[i,:],farray*vec[i+NSYS,:])
        bda = adb.conj()
        bdb = np.vdot(vec[i+NSYS,:],farray*vec[i+NSYS,:])
        
        # getting the magnetic moments
        MAG = 0.5*np.array([adb+bda,
            (-adb+bda)*1.0j,ada-bdb]).real
        rx = np.random_intel.uniform(UPAR[0][0],
                UPAR[0][1])
        M[1,i] = (1.0-rx)*M[1,i]+ rx*MAG

    for i in range(0,NSYS):

        # initiating Zero magnetic fields
        BMAG = np.array([0.0,0.0,0.0])
        # iteration over the neighbours
        for j in nbr[i]:

            # computing the bond operators
            if j > i:

                xbond = (np.vdot(vec[i,:],farray*vec[j,:])
                        +np.vdot(vec[i+NSYS,:],
                            farray*vec[j+NSYS,:]))
                tx = np.random_intel.uniform(UPAR[1][0],
                        UPAR[1][1])
                bond[1,i,j] = (1.0-tx)*bond[1,i,j] + tx*xbond
                bond[1,j,i] = bond[1,i,j].conj()
                # updating hopping matrix
                chi[i,j] = 0.5*(1.0-XPAR)*J[i,j]*bond[1,i,j]
                chi[j,i] = chi[i,j].conj()

            # computing the magnetic fields
            BMAG += XPAR*J[i,j]*M[1,j]

        # Updating B fields
        B[i] = BMAG


