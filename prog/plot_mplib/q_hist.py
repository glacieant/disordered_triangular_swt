#!/opt/intel/intelpython2/bin/python

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

a = np.pi

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

dvec = np.array([
    [-1.0,0.0],
    [-0.5,(3.0**0.5)/2],
    [0.5,(3.0**0.5)/2],
    [1.0,0.0],
    [0.5,-(3.0**0.5)/2],
    [-0.5,-(3.0**0.5)/2],
    ])

# defining the optimizing function

def qval(Q,D):
    val = ((np.cos(np.dot(Q,dvec[0]))-D[0])**2
            +(np.cos(np.dot(Q,dvec[1]))-D[1])**2
            +(np.cos(np.dot(Q,dvec[2]))-D[2])**2
            +(np.cos(np.dot(Q,dvec[3]))-D[3])**2
            +(np.cos(np.dot(Q,dvec[4]))-D[4])**2
            +(np.cos(np.dot(Q,dvec[5]))-D[5])**2
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

SIGMA = 4

UX = np.zeros(1,dtype=np.float)
VX = np.zeros(1,dtype=np.float)

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
    ay = 0.5*a

    if len(X) != N:
        
        # laying out the lattice skeleton
        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                # the rhombus lattice
                # X[i,j] = (i+j)*ax+3*ax
                # Y[i,j] = (j-i)*ay
                # the slanted lattice
                X[i,j] = i*a + j*ay
                Y[i,j] = j*ax 
   
    US = np.zeros(N,dtype=np.float)
    VS = np.zeros(N,dtype=np.float)

    #ERR = []
    QSOLX = 4.0*np.pi/3.0
    QSOLY = 0.0
    for i in range(0,N):
        RHS = np.zeros(6,dtype=np.float)
        for j in range(0,6):
            RHS[j] = np.dot(spin[i],spin[nbr[i,j]])
        
        SOL = minimize(qval,[QSOLX,QSOLY],
                args=(RHS),
                #bounds=((0,2*np.pi),(0,2*np.pi)),
                method='Powell',
                tol=1e-6
                )

        US[i], VS[i] = SOL.x

        #QSOLX = US[i]
        #QSOLY = VS[i]
        
    U = np.reshape(US,(L,L)).transpose()
    V = np.reshape(VS,(L,L)).transpose()

    # reprocessing to make data non-local
    Qx = np.zeros((L,L),dtype=np.float)
    Qy = np.zeros((L,L),dtype=np.float)
    for i in range(0,L):
        for j in range(0,L):
            NUM = 0
            for k in range(-SIGMA,SIGMA+1):
                for m in range(-SIGMA,SIGMA+1):
                    x = (i + k)%L
                    y = (j + m)%L
                    DIST = np.sqrt((X[x,y]-X[i,j])**2 
                            + (Y[x,y]-Y[i,j])**2)
                    if (DIST <= SIGMA*a) or (DIST > np.sqrt(2)*SIGMA*a):
                        Qx[i,j] += U[x,y]*gauss(DIST,SIGMA*a)
                        Qy[i,j] += V[x,y]*gauss(DIST,SIGMA*a)
                        NUM += gauss(DIST,SIGMA*a)
                        #NUM += 1
           
            Qx[i,j] *= 1.0/NUM
            Qy[i,j] *= 1.0/NUM
    
    U = Qx.flatten()
    V = Qy.flatten()

    U = (U - UZ)
    V = (V - VZ)

    UX = np.append(UX,U)
    VX = np.append(VX,V)

# plotting the histogram

w,h = figure.figaspect(1.0)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

# the histogram of the dara
ax.hist(UX,bins=50,normed=1,facecolor='royalblue')

plt.title(r'Histogram of $\Delta Q_x$')
fig.savefig("../plot/QXHIST.pdf"
        ,bbox_inches='tight'
        )
plt.close('all')

# plotting the histogram

w,h = figure.figaspect(1.0)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

# the histogram of the dara
ax.hist(VX,bins=50,normed=1,facecolor='royalblue')

plt.title(r'Histogram of $\Delta Q_y$')
fig.savefig("../plot/QYHIST.pdf"
        ,bbox_inches='tight'
        )
plt.close('all')

