#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = "sharp_wall"
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
avec = np.array([
    [-((3.0)**0.5)/2,0.5],
    [0.0,1.0],
    [((3.0)**0.5)/2,0.5],
    [((3.0)**0.5)/2,-0.5],
    [0.0,-1.0],
    [-((3.0)**0.5)/2,-0.5],
    ])*a/2


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

    VMIN = 0.0
    VMAX = 4.0
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
    
    # fixing colormap for line plotting

    Norm = Normalize(vmin=0,vmax=np.pi,clip=False)
    scalarMap = cm.ScalarMappable(norm=Norm,cmap=cm.viridis)
    
    # getting some colorbar

    # plotting the configuration

    #fig, ax = plt.subplots()
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])
    ax1 = fig.add_axes([0.6,0.05,0.35,0.02])

    cbar = colorbar.ColorbarBase(ax1,cmap=cm.viridis,
            norm=Norm,ticks=(0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),
            orientation='horizontal') 
    cbar.set_label(r'$\cos^{-1}(\vec{S}_i\cdot\vec{S}_{i+1})$',
            fontsize=12,
            labelpad=-45)
    cbar.set_ticklabels([r'$0$',
        r'$\pi/4$',
        r'$\pi/2$',
        r'$3\pi/4$',
        r'$\pi$',
        ])
   
    for i in range(0,N):

        for p in range(0,6):

            j = nbr[i,p]

            cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[p,0]]
            cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[p,1]]

            """

            # plotting the lattice
            line = plt.Line2D(cor_x,cor_y,
                    color='gray',
                    alpha=0.1,
                    ls='solid',
                    lw=0.25)
            ax.add_line(line)

            """

            cval = scalarMap.to_rgba(
                    np.arccos(np.dot(spin[i],spin[j])))

            #cval = scalarMap.to_rgba(i)
            # plotting the couplings
            line = plt.Line2D(cor_x,cor_y,
                    color=cval,
                    alpha=1,
                    ls='solid',
                    #lw=JLW[i,j],
                    lw=2*J[i,j],
                    zorder = 0
                    )
            ax.add_line(line)

        """

        # plotting additonal bonds
        for p in range(0,3):

            j = nbr[i,6+p]
            q = p +3

            cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
            cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]

            # plotting the couplings
            line = plt.Line2D(cor_x,cor_y,
                    color='red',
                    alpha=0.5,
                    ls='solid',
                    lw=JLW[i,j])
            ax.add_line(line)
        """
    # plotting the colorbar
    #fig.colorbar(CS3,shrink=0.5)
    #sm = ax1.cm.ScalarMappable(cmap=cm.viridis,norm=Norm)
    #sm._A = []
    #cbar = fig.colorbar(sm,shrink=0.5,ticks=(0,np.pi/4,np.pi/2,3*np.pi/4,np.pi))
   
    
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
    Q=ax.quiver(X,Y,U,V,
            color='red',
            pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            width=0.002,
            headwidth=3,
            headlength=5,
            alpha=1.0,
            zorder = 1.0
            )

    """
    Q0=ax.quiver(X,Y,U0,V0,
            color='blue',pivot='mid',
            angles='xy',
            scale=1,
            scale_units='xy',
            )
    """
    #plt.axis([X.min()-2.0*a,X.max()+2.0*a,
    #    Y.min()-2.0*a,Y.max()+2.0*a])
    ax.axis('off')

    """
    # adding label patch
    patch = [mpatches.Patch(color='red',alpha=0.5)]
    label = [r'Coupling, $J$']

    legend=ax.legend(patch,label,
            loc='lower right',shadow=True)

    """
    #plt.suptitle(r"Spin Configuration", 
    #        x=0.5, y=0.99, fontsize=16)
    #plt.title(r'Anamalous dipole density domains',x=0.625,fontsize=12)

    """
    plt.suptitle("Spin Configuration" 
            , x=0.5, y=0.99, fontsize=16)
    plt.title(r"$\Delta$ = "+STR_DELTA+" , "
            +r"$\alpha$ = "+STR_ALPHA+" , "
            +"Sample. = "+str(IDISD)+" , "
            +"Init. config. = "+str(BTNUM),x=0.485,fontsize=12)
    """
    #ax.set_aspect('equal')
    #fig.tight_layout(pad=0.5,h_pad=0.5,w_pad=0.5)
    fig.savefig("../plot/qdom_"+
            STR_L+
            STR_DELTA+
            STR_ALPHA+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf"
            #,bbox_inches='tight'
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


