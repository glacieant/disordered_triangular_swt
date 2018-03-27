#!/opt/intel/intelpython2/bin/python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = "sharp_wall"
#folder = "single_impurity"
folder = "zero_field_classical"
#folder = "zero_field_dilution"
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
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.optimize import root

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

a = 6*np.pi
#a = np.pi/64

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
"""
dvec = np.array([
    [-1.0,0.0],
    [-0.5,(3.0**0.5)/2],
    [0.5,(3.0**0.5)/2],
    [1.0,0.0],
    [0.5,-(3.0**0.5)/2],
    [-0.5,-(3.0**0.5)/2],
    ])
"""
dvec = np.array([
    [1.0,0.0],
    [0.5,(3.0**0.5)/2],
    [-0.5,(3.0**0.5)/2]
    ])
# defining the optimizing function

"""
def qval(Q,D):
    val = ((np.cos(np.dot(Q,dvec[0]))-D[0])**2
            +(np.cos(np.dot(Q,dvec[1]))-D[1])**2
            +(np.cos(np.dot(Q,dvec[2]))-D[2])**2
            +(np.cos(np.dot(Q,dvec[3]))-D[3])**2
            +(np.cos(np.dot(Q,dvec[4]))-D[4])**2
            +(np.cos(np.dot(Q,dvec[5]))-D[5])**2
            )
    return val
"""

def qval(Q,D):
    val = ((np.cos(np.dot(Q,dvec[0]))-D[0])**2
            +(np.cos(np.dot(Q,dvec[1]))-D[1])**2
            +(np.cos(np.dot(Q,dvec[2]))-D[2])**2
            )
    return val


# reference 120 degree Q        
RHS = np.zeros(6,dtype=np.float)
for j in range(0,6):
    RHS[j] = -0.5

UZ, VZ = minimize(qval,
        [4.0*np.pi/3,0.0],
        args=(RHS),
        method='Powell'
        ).x

# gaussian envelope

def gauss(x,sigma):

    return ((1.0/(sigma*np.sqrt(2*np.pi)))
            *np.exp(-(x**2)/(2*(sigma**2))))

# gaussian parameters

SIGMA = 10

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

    N = L**2

    nbr = np.zeros((N,ZCO),dtype=np.int)
    lmap.lattice_map(L,nbr)

    LTR = L-1
    NTR = LTR**2

    elt = np.zeros((NTR,3),dtype=np.int)
    lmap.eltriangle(L,elt)

    ax = a*(3.0**0.5)/2.0
    ay = 0.5*a

    if len(X) != NTR:
        
        # laying out the lattice skeleton
        X = np.zeros((LTR,LTR))
        Y = np.zeros((LTR,LTR))
        for i in range(0,LTR):
            for j in range(0,LTR):
                # the rhombus lattice
                # X[i,j] = (i+j)*ax+3*ax
                # Y[i,j] = (j-i)*ay
                # the slanted lattice
                X[i,j] = i*a + j*ay
                Y[i,j] = j*ax 
   
    chir = np.zeros((NTR,3),dtype=np.float) 

    for i in range(0,NTR):
        chir[i] = np.cross(spin[elt[i,0]],
                spin[elt[i,1]])
        chir[i] += np.cross(spin[elt[i,1]],
                spin[elt[i,2]])
        chir[i] += np.cross(spin[elt[i,2]],
                spin[elt[i,0]])

    e1 = chir[0]
    """
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,chir[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
    """
    e1 = e1/np.linalg.norm(e1)
    e2 = np.cross(e1,np.cross(e1,[0.0,0.0,1.0]))

    U = np.einsum('ij,j->i',chir,e1).reshape((LTR,LTR)).transpose()
    V = np.einsum('ij,j->i',chir,e2).reshape((LTR,LTR)).transpose()

    eta = 0.0000001
    #eU = U/np.sqrt(U*U+V*V+eta)
    #eV = V/np.sqrt(U*U+V*V+eta)
    UVNORM = np.sqrt(U*U+V*V)
    eU = U
    eV = V
    UVA = np.arctan2(eV+eta,eU+eta)
   

    # fixing colormap for line plotting

    cmstyle = cm.hsv
    VMIN = -np.pi
    VMAX = np.pi

    Norm = Normalize(vmin=VMIN,vmax=VMAX,clip=False)
    scalarMap = cm.ScalarMappable(norm=Norm,cmap=cmstyle)
    
    # getting some colorbar

    # plotting the configuration

    #fig, ax = plt.subplots()
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1.0/(3.0)**(0.5)])
    #ax = fig.add_axes([0,0,1,1])
    ax1 = fig.add_axes([0.6,0.05,0.35,0.02])

    cbar = colorbar.ColorbarBase(ax1,cmap=cmstyle,
            norm=Norm,
            ticks=(VMIN,
                VMIN+(VMAX-VMIN)/4,
                VMIN+2*(VMAX-VMIN)/4,
                VMIN+3*(VMAX-VMIN)/4,
                VMAX),
            orientation='horizontal') 
   
    
    #cbar.set_label(r'$|\vec{Q}_i-\vec{Q}_0|$',
    #        fontsize=12,
    #        labelpad=-45)
    cbar.set_label(r'$\tan^{-1}\left(\frac{\Delta Q^y}{\Delta Q^x}\right)$',
            fontsize=12,
            labelpad=-45)

    cbar.set_ticklabels([r'$-\pi$',
        r'$-\pi/2$',
        r'$0$',
        r'$\pi/2$',
        r'$\pi$',
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

    Q=ax.quiver(X,Y,eU,eV,UVA,
            cmap=cmstyle,
            norm=Norm,
            pivot='mid',
            angles='xy',
            scale=0.1,
            scale_units='xy',
            width=0.002,
            headwidth=4,
            headlength=5,
            alpha=1.0,
            zorder = 1.0
            )
    """
    DELTA = float(STR_DELTA)
    ax.text(X.min(),Y.max(),
            r'$\frac{\delta J}{J}$='
            +str("%.1f" % DELTA),
            fontsize=20
            )
    """
    QMAX = np.amax(UVNORM)
    ax.text(X.min(),Y.max(),
            r'$\Delta Q_{\textrm{max}}$='
            +str("%.4f" % QMAX),
            fontsize=20
            )
    ax.axis('off')

    fig.savefig("../plot/chirdom"+
            "_L_"+
            STR_L+
            "_DLT_"+
            STR_DELTA+
            "_ALP_"+
            STR_ALPHA+
            "_DNM_"+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf",
            bbox_inches='tight',
            transparent=True
            )
    plt.close('all')


for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite '+
                    '-dPDSETTINGS=/prepress -sOutputFile='+
                '../plot/CHIRDOMAIN_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/chirdom'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '* ',shell=True)
            subprocess.call('rm ../plot/chirdom'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '*',shell=True)

