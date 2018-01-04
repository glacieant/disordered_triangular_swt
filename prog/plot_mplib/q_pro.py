#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = "sharp_wall"
#folder = "single_impurity"
folder = "zero_field_classical"
os.chdir("../"+folder+"/out/data")
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

a = 0.5

# three diffrent translation vector

"""
avec = np.array([[-0.5,((3.0)**0.5)/2],
    [0.5,((3.0)**0.5)/2],
    [1.0,0.0],
    [-1.5,((3.0)**0.5)/2],
    [0.0,(3.0**0.5)],
    [1.5,((3.0)**0.5)/2]])*a
"""

evec = np.array([
    [-((3.0)**0.5)/2,0.5],
    [0.0,1.0],
    [((3.0)**0.5)/2,0.5],
    [((3.0)**0.5)/2,-0.5],
    [0.0,-1.0],
    [-((3.0)**0.5)/2,-0.5],
    ])

avec = evec*a/2


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

    ax = a*(3.0**0.5)/2.0
    ay = a/2.0

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                #X[i,j]=(i+j*(1.0/2.0))*a +2*a
                #Y[i,j]=j*((3.0)**0.5/2.0)*a + 2*a
                X[i,j]= (i+j)*ax+3*ax
                Y[i,j]= (j-i)*ay
    """
    # defining the plane for the spins
    e1 = spin[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    U = 0.5*np.einsum('ij,j->i',spin,e1).reshape((L,L)).transpose()
    V = 0.5*np.einsum('ij,j->i',spin,e2).reshape((L,L)).transpose()
    """

    US = np.zeros(N,dtype=np.float)
    VS = np.zeros(N,dtype=np.float)
    
    for i in range(0,N):
        p = 0
        q = 2
        s1 = np.arccos(np.dot(spin[i],spin[nbr[i,p]]))
        s2 = np.arccos(np.dot(spin[i],spin[nbr[i,q]]))
        US[i] = s1 + s2
        VS[i] = (s2 - s1)/np.sqrt(3)
    
    U = np.reshape(US,(L,L)).transpose()#%(2.0*np.pi)#)/(2.0*np.pi)
    V = np.reshape(VS,(L,L)).transpose()#%(2.0*np.pi)#)/(2.0*np.pi)

    #print U
    #print "====Shit====="
    #print V
    

    # fixing colormap for line plotting

    cmstyle = cm.plasma
    VMIN = np.pi/2
    VMAX = 3.0*np.pi/2

    Norm = Normalize(vmin=VMIN,vmax=VMAX,clip=False)
    scalarMap = cm.ScalarMappable(norm=Norm,cmap=cmstyle)
    
    # getting some colorbar

    # plotting the configuration

    #fig, ax = plt.subplots()
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1.0/(3.0)**(0.5)])
    ax1 = fig.add_axes([0.6,0.05,0.35,0.02])

    cbar = colorbar.ColorbarBase(ax1,cmap=cmstyle,
            norm=Norm,
            ticks=(VMIN,
                VMIN+(VMAX-VMIN)/4,
                VMIN+2*(VMAX-VMIN)/4,
                VMIN+3*(VMAX-VMIN)/4,
                VMAX),
            orientation='horizontal') 
    cbar.set_label(r'$|\vec{Q}_i|$',
            fontsize=12,
            labelpad=-45)
    cbar.set_ticklabels([r'$\pi/2$',
        r'$3\pi/4$',
        r'$\pi$',
        r'$5\pi/4$',
        r'$3\pi/2$',
        ])
   
    """
    for i in range(0,N):

        for p in range(0,6):

            j = nbr[i,p]

            cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[p,0]]
            cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[p,1]]

            #cval = scalarMap.to_rgba(
            #        np.arccos(np.dot(spin[i],spin[j])))

            #cval = scalarMap.to_rgba(i)
            # plotting the couplings
            line = plt.Line2D(cor_x,cor_y,
                    color='gray',
                    alpha=1,
                    ls='solid',
                    #lw=JLW[i,j],
                    lw=J[i,j],
                    zorder = 0
                    )
            ax.add_line(line)
    """
    # picturing the spin orientation
    UVNORM = np.sqrt(U*U+V*V)
    eU = a*U/np.sqrt(U*U+V*V)
    eV = a*V/np.sqrt(U*U+V*V)

    Q=ax.quiver(X,Y,eU,eV,UVNORM,
            cmap=cmstyle,
            norm=Norm,
            pivot='tail',
            angles='xy',
            scale=1,
            scale_units='xy',
            width=0.002,
            headwidth=4,
            headlength=5,
            alpha=1.0,
            zorder = 1.0
            )

    ax.axis('off')

    fig.savefig("../plot/qdom_"+
            STR_L+
            STR_DELTA+
            STR_ALPHA+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf"
            ,bbox_inches='tight'
            )
    plt.close('all')

for L in zip(*hshchar)[0]:
    for DLT in zip(*hshchar)[1]:
        for ALP in zip(*hshchar)[2]:
            subprocess.call('pdftk ../plot/qdom_'+
                L+
                DLT+
                ALP+
                '* '+
                'cat output ../plot/QDOMAIN_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf',shell=True)

# removing split files
subprocess.call('rm ../plot/qdom_*',shell=True)


