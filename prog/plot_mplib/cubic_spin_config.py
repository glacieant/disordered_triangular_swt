#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
folder = "cubic_classical"
os.chdir("../"+folder+"/out/data_bk")
import sys
import subprocess
import re
import glob
import numpy as np
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import matplotlib.figure as figure
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
    td_spin = FNDATA['spin']
    td_spin0 = FNDATA['spin0']

    N = L**3

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                X[i,j]= i*a
                Y[i,j]= j*a
    # defining the plane for the spins
    
    e1 = td_spin0[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,td_spin0[1]))
        if np.linalg.norm(e2) > 10.0**(-5):
            e2 = e2/np.linalg.norm(e2)
        else:
            e2 = np.array([e1[1],-e1[0],0.0])
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])


    for zix in range((L+1)/2,(L+1)/2+1):
    #for zix in range(0,L):

        # plotting the configuration

        w,h = figure.figaspect(1)
        fig = plt.figure(figsize=(w,h))
        ax = fig.add_axes([0,0,1,1])

        ix = (zix)*(L**2)
        jx = ix + L**2

        spin = td_spin[ix:jx]

        U = (3.0*a/4)*np.einsum('ij,j->i',spin,e1).reshape((L,L)).transpose()
        V = (3.0*a/4)*np.einsum('ij,j->i',spin,e2).reshape((L,L)).transpose()

        spin0 = td_spin0[ix:jx]

        #spin0[15] = [0,0,0]

        U0 = (3.0*a/4)*np.einsum('ij,j->i',spin0,e1).reshape((L,L)).transpose()
        V0 = (3.0*a/4)*np.einsum('ij,j->i',spin0,e2).reshape((L,L)).transpose()

        RSU = np.zeros((L,L),dtype=np.float)
        RSV = np.zeros((L,L),dtype=np.float)
        GSU = np.zeros((L,L),dtype=np.float)
        GSV = np.zeros((L,L),dtype=np.float)
        BSU = np.zeros((L,L),dtype=np.float)
        BSV = np.zeros((L,L),dtype=np.float)
        
        SU0 = np.zeros((L,L),dtype=np.float)
        SV0 = np.zeros((L,L),dtype=np.float)

        SX = np.zeros((L,L),dtype=np.float)
        SY = np.zeros((L,L),dtype=np.float)
        
        sub = np.zeros(L*L,dtype=np.int)
        lmap.sublattice_map(L,sub)

        for i in range(0,L*L):

            x = i%L
            y = i/L

            if ( x>=a and x <= 7*a ) \
                    and ( y>=a and y <= 7*a):
                
                SX[x,y] = X[x,y]
                SY[x,y] = Y[x,y]
                SU0[x,y] = U0[x,y]
                SV0[x,y] = V0[x,y]
                if sub[i] == 0:
                    RSU[x,y] = U[x,y]
                    RSV[x,y] = V[x,y]
                elif sub[i] == 1:
                    GSU[x,y] = U[x,y]
                    GSV[x,y] = V[x,y]
                else:
                    BSU[x,y] = U[x,y]
                    BSV[x,y] = V[x,y]
                    
        for i in range(1,L-1):

            for j in range(1,L-1):

                for NN in [-1,1]:

                    cor_x = [X[i,j],X[i,j]+NN*a/2]
                    cor_y = [Y[i,j],Y[i,j]]

                    p = i + j*L + ix
                    if i + NN >= 0:
                        q = (i+NN)%L + j*L + ix
                    elif i + NN < 0:
                        q = (i+NN+L) + j*L + ix

                    # plotting the couplings
                    line = plt.Line2D(cor_x,cor_y,
                            color='black',
                            alpha=1,
                            ls='solid',
                            lw=2*J[p,q],
                            zorder = 0
                            )
                    ax.add_line(line)

                    cor_x = [X[i,j],X[i,j]]
                    cor_y = [Y[i,j],Y[i,j]+NN*a/2]

                    p = i + j*L + ix
                    if j + NN >= 0:
                        q = i + ((j+NN)%L)*L + ix
                    elif j + NN < 0:
                        q = i + (j+NN+L)*L + ix

                    # plotting the couplings
                    line = plt.Line2D(cor_x,cor_y,
                            color='black',
                            alpha=1.0,
                            ls='solid',
                            lw=2*J[p,q],
                            zorder = 1
                            )
                    ax.add_line(line)

        # picturing the spin orientation
        pivot="mid"
        width=0.008
        headwidth=5
        headlength=6

        LT = ax.scatter(SX,SY,
                color='k',
                zorder=3)

        RQ=ax.quiver(SX,SY,RSU,RSV,
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
        GQ=ax.quiver(SX,SY,GSU,GSV,
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
        BQ=ax.quiver(SX,SY,BSU,BSV,
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
                linewidth=1,
                facecolor='none',
                edgecolor='k',
                width=width,
                headwidth=headwidth,
                headlength=headlength,
                pivot=pivot,
                #angles='xy',
                scale=1,
                scale_units='xy',
                zorder=1)

        plt.axis(([a/2,15*a/2,a/2,15*a/2]))
     
        ax.axis('off')

        fig.savefig("../plot/cubic_spin_config"+
                "_L_"+
                STR_L+
                "_DLT_"+
                STR_DELTA+
                "_ALP_"+
                STR_ALPHA+
                "_DNM_"+
                str("%06d" % IDISD)+
                str("%06d" % BTNUM)+
                ".pdf"
                ,bbox_inches='tight'
                )
        plt.close('all')

for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite '+
                    '-dPDSETTINGS=/prepress -sOutputFile='+
                '../plot/CUBIC_SPIN_CONFIG_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/cubic_spin_config'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '* ',shell=True)
            subprocess.call('rm ../plot/cubic_spin_config'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '*',shell=True)

