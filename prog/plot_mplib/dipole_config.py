#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
folder = "single_impurity"
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

a = 0.5

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

    #fig, ax = plt.subplots()
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])
    
    # defining the plane for the spins

    R11 = -np.sqrt(3.0)/2.0
    R12 = 1.0/2.0

    e1 = spin[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    UP = 0.5*np.einsum('ij,j->i',spin,e1)
    VP = 0.5*np.einsum('ij,j->i',spin,e2)

    U = np.reshape(R11*UP + R12*VP,(L,L)).transpose()
    V = np.reshape(-R12*UP + R11*VP,(L,L)).transpose()

    e1 = spin0[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin0[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    UP0 = 0.5*np.einsum('ij,j->i',spin0,e1)
    VP0 = 0.5*np.einsum('ij,j->i',spin0,e2)

    U0 = np.reshape(R11*UP0 + R12*VP0,(L,L)).transpose()
    V0 = np.reshape(-R12*UP0 + R11*VP0,(L,L)).transpose()

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

    """
    VMIN = 0.0
    VMAX = 2.0
    VGRID = 20
    
    JMIN = np.amin(J)
    JMAX = np.amax(J)

    JWDTH = (JMAX-JMIN)/VGRID

    JLW = np.zeros((N,N))
    
    for i in range(0,N):

        for j in nbr[i]:

            JLW[i,j] = (VMIN + 
                    (np.floor((J[i,j]-JMIN)/JWDTH)
                        /VGRID)*(VMAX-VMIN)
                    )
    """
    impsite = (N+1)/2-1
    impcor = np.array([X[impsite%L,impsite/L]+0.5*a,Y[impsite%L,impsite/L]])
    for i in range(0,N):

        x = i%L
        y = i/L

        cor = np.array([X[x,y],Y[x,y]])
        dixt = cor-impcor

        d = np.linalg.norm(dixt)

        if d < 4*a:

            SX[x,y] = X[x,y]
            SY[x,y] = Y[x,y]

            SU0[x,y] = U0[x,y]
            SV0[x,y] = V0[x,y]

            if dixt[0] < 0.0:
                LSU[x,y] = U[x,y]
                LSV[x,y] = V[x,y]
            elif dixt[0] > 0.0:
                RSU[x,y] = U[x,y]
                RSV[x,y] = V[x,y]
            else:
                SU[x,y] = U[x,y]
                SV[x,y] = V[x,y]




            # plotting the lattice
            #line = plt.Line2D(cor_x,cor_y,
            #        color='gray',
            #        alpha=0.1,
            #        ls='solid',
            #        lw=0.25)
            #ax.add_line(line)

            # plotting the couplings
            #line = plt.Line2D(cor_x,cor_y,
            #        color='red',
            #        alpha=0.5,
            #        ls='solid',
            #        lw=JLW[i,j])
            #ax.add_line(line)

        # plotting additonal bonds
        """
        for p in range(0,3):

            j = nbr[i,6+p]
            q = p +3

            cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
            cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]
        """
            # plotting the couplings
            #line = plt.Line2D(cor_x,cor_y,
            #        color='red',
            #        alpha=0.5,
            #        ls='solid',
            #        lw=JLW[i,j])
            #ax.add_line(line)

    """
    # energy etimate
    ENX = 0.0
    for i in range(0,N):
        for k in range(0,6):
            j = nbr[i,k]
            ENX += J[i,j]*np.dot(spin[i],spin[j])

    ENX = ENX/N

    ax.text(-1.25*a,(3.5*L/4.5)*a,
            r'$E/N$ = '+str("%.4f" % ENX)
            ,
            fontsize=12)

    # labeling scale of coupling and bonds
    ax.text(-1.25*a,(3.0*L/4.5)*a,
            r'$J_{\mathrm{max}}$ = '+str("%.4f" % JMAX)
            ,
            fontsize=12)
    """
    # picturing the spin orientation
    LQ=ax.quiver(SX,SY,LSU,LSV,
            color='red',pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )
    RQ=ax.quiver(SX,SY,RSU,RSV,
            color='blue',pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )
    Q=ax.quiver(SX,SY,SU,SV,
            color='black',pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )

    Q0=ax.quiver(SX,SY,SU0,SV0,
            color='gray',
            edgecolor='gray',
            #linestyle='dashed',
            #linewidth=2,
            #facecolor='none',
            #width=0.0001,
            #headwidth=100,
            #headlength=100,
            #hatch='ooo',
            pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            alpha=0.5,
            zorder=1
            )
    
    ellipse = Ellipse(xy=(impcor[0],impcor[1]),width=1.25*a,height=0.25*a,
            edgecolor='palegreen',facecolor='palegreen',
            zorder=0)
    ax.add_patch(ellipse)
    plt.axis([impcor[0]-4.5*a,impcor[0]+4.5*a,
        impcor[1]-4.5*a,impcor[1]+4.5*a])
    ax.axis('off')
    
    # adding label patch
    """
    patch = [patches.Ellipse(xy={0,0},width=1.25,height=0.25,color='palegreen',alpha=0.5)]
    label = [r'Dipolar perturbation']

    legend=ax.legend(patch,label,
            loc='lower right',shadow=True)

    plt.suptitle(r"Spin Configuration", 
            x=0.5, y=0.99, fontsize=16)
    plt.title(r'Anomalous dipole density',x=0.485,fontsize=12)

    
    plt.suptitle("Spin Configuration" 
            , x=0.5, y=0.99, fontsize=16)
    plt.title(r"$\Delta$ = "+STR_DELTA+" , "
            +r"$\alpha$ = "+STR_ALPHA+" , "
            +"Sample. = "+str(IDISD)+" , "
            +"Init. config. = "+str(BTNUM),x=0.485,fontsize=12)
    """
    #fig.tight_layout(pad=1.6,h_pad=1.6,w_pad=1.6)
    fig.savefig("../plot/esbc_"+
            STR_L+
            STR_DELTA+
            STR_ALPHA+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf",
            bbox_inches='tight'
            )
    plt.close('all')

for L in zip(*hshchar)[0]:
    for DLT in zip(*hshchar)[1]:
        for ALP in zip(*hshchar)[2]:
            subprocess.call('pdftk ../plot/esbc_'+
                L+
                DLT+
                ALP+
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


