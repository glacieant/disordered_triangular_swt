#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
import sys
import subprocess
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# the font styleset
from matplotlib import rcParams
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.family'] = 'serif'

## picturing the output data

# file name pattern

fpat = re.compile('([^/]*)X([^/]*)X([^/]*)X([^/]*)X([^/]*).dat')

# importing lattice geometry system
import lattice_map as lmap
ZCO = 12

X = np.zeros(1)
Y = np.zeros(1)
Z = np.zeros(1)

a = 1.0

# three diffrent translation vector

avec = np.array([[-0.5,((3.0)**0.5)/2],
    [0.5,((3.0)**0.5)/2],
    [1.0,0.0],
    [-1.5,((3.0)**0.5)/2],
    [0.0,(3.0**0.5)],
    [1.5,((3.0)**0.5)/2]])*a

for fname in glob.iglob('*.dat'):

    match = fpat.match(fname)
    LS = match.group(1)
    CL_EN = match.group(2)
    Q_EN = match.group(3)
    DELTA = match.group(4)
    ALPHA = match.group(5)


    L = int(LS)
    N = L**2

    nbr = np.zeros((N,ZCO),dtype=np.int)
    lmap.lattice_map(L,nbr)

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                X[i,j]=(i+j*(1.0/2.0))*a
                Y[i,j]=j*((3.0)**0.5/2.0)*a

    # plotting the configuration

    fig, ax = plt.subplots()

    # defining the plane for the spins

    TH, PH = np.loadtxt(fname,usecols=(0,1),unpack=True)

    S_X = np.sin(TH)*np.cos(PH)
    S_Y = np.sin(TH)*np.sin(PH)
    S_Z = np.cos(TH)

    spin = np.column_stack((S_X,S_Y,S_Z))

    e1 = spin[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    U = np.einsum('ij,j->i',spin,e1).reshape((L,L)).transpose()
    V = np.einsum('ij,j->i',spin,e2).reshape((L,L)).transpose()

    for i in range(0,N):

        for p in range(1,4):

            j = nbr[i,p]
            q = p -1

            cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
            cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]

            # plotting the lattice
            line = plt.Line2D(cor_x,cor_y,
                    color='gray',
                    alpha=0.25,
                    ls='solid',
                    lw=0.3)
            ax.add_line(line)

    # picturing the spin orientation
    Q=ax.quiver(X,Y,U,V,
            color='black',pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            )
    plt.axis([X.min()-2.0*a,X.max()+2.0*a,
        Y.min()-2.0*a,Y.max()+2.0*a])

    plt.suptitle("Spin Configuration for L="+LS,
            x=0.5, y=0.99,fontsize=16)
    plt.title(r"$E_{\mathrm{cl}}/L^{2}$ = "+CL_EN
            +"       "
            +r"$E_{\mathrm{Q}}/L^{2}$ = "+Q_EN
            ,x=0.485,fontsize=12)
    plt.xlabel(r'$\Delta$ ='+DELTA
            +"      "
            +r'$\alpha$ ='+ALPHA
    ,fontsize=12)
    fig.tight_layout(pad=1.6,h_pad=1.6,w_pad=1.6)
    fig.savefig(LS+"_"+DELTA+"_"+ALPHA+"_"+CL_EN+".pdf")
    plt.close('all')
