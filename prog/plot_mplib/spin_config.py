#!/opt/intel/intelpython2/bin/python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = "sharp_wall"
folder = "single_impurity"
#folder = "zero_field_classical"
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
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')

# the font styleset
from matplotlib import rcParams
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.family'] = 'serif'

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

"""
avec = np.array([[-0.5,((3.0)**0.5)/2],
    [0.5,((3.0)**0.5)/2],
    [1.0,0.0],
    [-1.5,((3.0)**0.5)/2],
    [0.0,(3.0**0.5)],
    [1.5,((3.0)**0.5)/2]])*a
"""

avec = np.array([
    [-1.0,0.0],
    [-0.5,(3.0**0.5)/2],
    [0.5,(3.0**0.5)/2],
    [1.0,0.0],
    [0.5,-(3.0**0.5)/2],
    [-0.5,-(3.0**0.5)/2],
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
    sub = np.zeros(N,dtype=np.int)
    lmap.lattice_map(L,nbr)
    lmap.sublattice_map(L,sub)

    ax = a*(3.0**0.5)/2.0
    ay = a/2.0

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                X[i,j]=(i+j*(1.0/2.0))*a +2*a
                Y[i,j]=j*((3.0)**0.5/2.0)*a + 2*a
                #X[i,j]= (i+j)*ax+3*ax
                #Y[i,j]= (j-i)*ay

    # plotting the configuration

    #fig, ax = plt.subplots()
    w,h = figure.figaspect(1)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1.0/(3.0)**(0.5)])
    #ax = fig.add_axes([0,0,1,1])
    ax1 = fig.add_axes([0.6,0.05,0.35,0.02])
    # defining the plane for the spins

    #R1 = np.array([[-1.0,0.0],[0.0,1.0]])
    #phx = -np.pi/6
    #R2 = np.array([[np.cos(phx),-np.sin(phx)],[np.sin(phx),np.cos(phx)]])
    #R = np.dot(R1,R2)

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

    #U = np.reshape(R[0,0]*UP + R[0,1]*VP,(L,L)).transpose()
    #V = np.reshape(R[1,0]*UP + R[1,1]*VP,(L,L)).transpose()

    SX = np.zeros((L,L),dtype=np.float)
    SY = np.zeros((L,L),dtype=np.float)

    SU = np.zeros((L,L),dtype=np.float)
    SV = np.zeros((L,L),dtype=np.float)
    LSU = np.zeros((L,L),dtype=np.float)
    LSV = np.zeros((L,L),dtype=np.float)
    RSU = np.zeros((L,L),dtype=np.float)
    RSV = np.zeros((L,L),dtype=np.float)

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
    
    # fixing colormap for line plotting

    cmstyle = cm.Greys

    Norm = Normalize(vmin=VMIN,vmax=VMAX,clip=False)
    scalarMap = cm.ScalarMappable(norm=Norm,cmap=cmstyle)
    
    # getting some colorbar

    cbar = colorbar.ColorbarBase(ax1,cmap=cmstyle,
            norm=Norm,
            ticks=(VMIN,
                VMIN+(VMAX-VMIN)/4,
                VMIN+2*(VMAX-VMIN)/4,
                VMIN+3*(VMAX-VMIN)/4,
                VMAX),
            orientation='horizontal') 
    cbar.set_label(r'$J$',
            fontsize=12,
            labelpad=-45)
    """
    cbar.set_ticklabels([r'$\pi/2$',
        r'$3\pi/4$',
        r'$\pi$',
        r'$5\pi/4$',
        r'$3\pi/2$',
        ])
    """
    for i in range(0,N):

        x = i/L
        y = i%L

        cor = np.array([X[x,y],Y[x,y]])

        SX[x,y] = X[x,y]
        SY[x,y] = Y[x,y]

        #SU0[x,y] = U0[x,y]
        #SV0[x,y] = V0[x,y]
        
        if sub[i]==0:
            LSU[x,y] = U[x,y]
            LSV[x,y] = V[x,y]
        elif sub[i]==1:
            RSU[x,y] = U[x,y]
            RSV[x,y] = V[x,y]
        else:
            SU[x,y] = U[x,y]
            SV[x,y] = V[x,y]


        for p in range(0,6):

            j = nbr[i,p]

            cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[p,0]]
            cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[p,1]]

            cval = scalarMap.to_rgba(JLW[i,j])

            # plotting the couplings
            line = plt.Line2D(cor_x,cor_y,
                    color=cval,
                    alpha=1,
                    ls='solid',
                    #lw=0.5*J[i,j],
                    zorder = 0
                    )
            ax.add_line(line)

    # picturing the spin orientation
    pivot="mid"
    width=0.002
    headwidth=4
    headlength=5

    LQ=ax.quiver(SX,SY,LSU,LSV,
            color='red',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )
    RQ=ax.quiver(SX,SY,RSU,RSV,
            color='green',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )
    Q=ax.quiver(SX,SY,SU,SV,
            color='blue',
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            pivot=pivot,
            angles='xy',
            scale=1,
            scale_units='xy',
            zorder=2
            )

    #plt.axis([X.min()-2.0*a,X.max()+2.0*a,
    #    Y.min()-2.0*a,Y.max()+2.0*a])
    ax.axis('off')

    fig.savefig("../plot/spin_config"+
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
                '../plot/SPIN_CONFIG_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/spin_config'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '* ',shell=True)
            subprocess.call('rm ../plot/spin_config'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '*',shell=True)

