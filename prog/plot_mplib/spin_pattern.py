#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
folder = "single_impurity"
#folder = "zero_field_classical"
os.chdir("../"+folder+"/out/data")
import sys
import subprocess
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.figure as figure
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# the font styleset
#from matplotlib import rcParams
#rcParams['font.serif'] = ['Times New Roman']
#rcParams['font.family'] = 'serif'

## picturing the output data

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*)_DISD_([^/]*)_BOOT_([^/]*).npz')

if not os.path.isdir("../plot"):
    os.mkdir("../plot")

HASH = set([])

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    STR_L = match.group(1)
    STR_DELTA = match.group(2)
    STR_ALPHA = match.group(3)
    HSH = STR_L+STR_DELTA+STR_ALPHA
    HASH.add(HSH)

HASH = list(HASH)

# the master arrays of data
HSHNUM = len(HASH)
hshchar = [0]*HSHNUM
sfc = [0]*HSHNUM
n_sfc = np.zeros(HSHNUM)
en = np.zeros(HSHNUM)

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

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    STR_L = match.group(1)
    L = int(STR_L)
    STR_DELTA = match.group(2)
    STR_ALPHA = match.group(3)

    HSHCHK = STR_L+STR_DELTA+STR_ALPHA
    INDX = HASH.index(HSHCHK)
    hshchar[INDX] = [STR_L,STR_DELTA,STR_ALPHA]

    FNDATA = np.load(fname)

    IDISD = int(match.group(4))
    BTNUM = int(match.group(5))
   
    J = FNDATA['J']
    spin = FNDATA['spin']
    spin0 = FNDATA['spin0']

    N = L**2

    nbr = np.zeros((N,ZCO),dtype=np.int)
    sub = np.zeros(N,dtype=np.int)
    lmap.lattice_map(L,nbr)
    lmap.sublattice_map(L,sub)

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                X[i,j]=(i+j*(1.0/2.0))*a
                Y[i,j]=j*((3.0)**0.5/2.0)*a

    # plotting the configuration

    #fig, ax = plt.subplots()
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])
    
    # defining the plane for the spins

    R1 = np.array([[-1.0,0.0],[0.0,1.0]])
    phx = -np.pi/6
    R2 = np.array([[np.cos(phx),-np.sin(phx)],[np.sin(phx),np.cos(phx)]])
    R = np.dot(R1,R2)

    e1 = spin[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    UP = np.einsum('ij,j->i',spin,e1)
    VP = np.einsum('ij,j->i',spin,e2)

    U = np.reshape(R[0,0]*UP + R[0,1]*VP,(L,L)).transpose()
    V = np.reshape(R[1,0]*UP + R[1,1]*VP,(L,L)).transpose()

    e1 = spin0[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin0[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    UP0 = np.einsum('ij,j->i',spin0,e1)
    VP0 = np.einsum('ij,j->i',spin0,e2)

    U0 = np.reshape(R[0,0]*UP0 + R[0,1]*VP0,(L,L)).transpose()
    V0 = np.reshape(R[1,0]*UP0 + R[1,1]*VP0,(L,L)).transpose()

    SX = np.zeros((L,L),dtype=np.float)
    SY = np.zeros((L,L),dtype=np.float)

    SU = np.zeros((L,L),dtype=np.float)
    SV = np.zeros((L,L),dtype=np.float)
    LSU = np.zeros((L,L),dtype=np.float)
    LSV = np.zeros((L,L),dtype=np.float)
    RSU = np.zeros((L,L),dtype=np.float)
    RSV = np.zeros((L,L),dtype=np.float)
    SU0 = np.zeros((L,L),dtype=np.float)
    SV0 = np.zeros((L,L),dtype=np.float)

    impsite = (L*(L-1))/2-1
    impcor = np.array([X[impsite%L,impsite/L]+0.5*a,Y[impsite%L,impsite/L]])
    for i in range(0,N):

        x = i%L
        y = i/L
        
        """
        cor = np.array([X[x,y],Y[x,y]])
        dixt = cor-impcor
        dx = (5)*a
        dy = (5)*(np.sqrt(3.0)*a/2)

        if ((np.abs(dixt[0])/dx+np.abs(dixt[1])/dy)<=a):
        """
        
        SX[x,y] = X[x,y]
        SY[x,y] = Y[x,y]

        SU0[x,y] = U0[x,y]
        SV0[x,y] = V0[x,y]

        if sub[i] == 0:
            LSU[x,y] = U[x,y]
            LSV[x,y] = V[x,y]
        elif sub[i] == 1:
            RSU[x,y] = U[x,y]
            RSV[x,y] = V[x,y]
        else:
            SU[x,y] = U[x,y]
            SV[x,y] = V[x,y]

    # picturing the spin orientation
    pivot="mid"
    width=0.002
    headwidth=5
    headlength=6

    #LT = ax.scatter(SX,SY,
    #        color='k',
    #        zorder=3)

    LQ=ax.quiver(SX,SY,LSU,LSV,
            color='red',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            #angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )
    RQ=ax.quiver(SX,SY,RSU,RSV,
            color='limegreen',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            #angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )
    Q=ax.quiver(SX,SY,SU,SV,
            color='royalblue',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            #angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )

    Q0=ax.quiver(SX,SY,SU0,SV0,
            color='black',
            linewidth=0.5,
            facecolor='none',
            #edgecolor='k',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            #angles='xy',
            scale=1,
            scale_units='xy',
            zorder=1
            )
    
    ellipse = Ellipse(xy=(impcor[0],impcor[1]),width=1.75*a,height=0.4*a,
            edgecolor='navajowhite',facecolor='navajowhite',
            zorder=0)
    ax.add_patch(ellipse)
    plt.axis([X.min()- 2*a,X.max() + 2*a,
        Y.min()- 2*a,Y.max() + 2*a])
    ax.axis('off')

    fig.savefig("../plot/esbc_"+
            STR_L+
            STR_DELTA+
            STR_ALPHA+
            '_DNMR_'+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf",
            #bbox_inches='tight'
            )
    plt.close('all')

for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('pdftk ../plot/esbc_'+
                L+
                DLT+
                ALP+
                '_DNMR_'+
                '* '+
                'cat output ../plot/CONFIG_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf',shell=True)

# removing split files
subprocess.call('rm ../plot/esbc_*',shell=True)


