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
from scipy.stats import norm
from scipy.optimize import curve_fit

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
UX = [0]*HSHNUM
VX = [0]*HSHNUM
DELTA = [0]*HSHNUM

# importing lattice geometry system
import lattice_map as lmap
ZCO = 12

a = np.pi

X = np.zeros(1)
Y = np.zeros(1)
Z = np.zeros(1)

dvec = np.array([
    [1.0,0.0],
    [0.5,(3.0**0.5)/2],
    [-0.5,(3.0**0.5)/2]
    ])
# defining the optimizing function

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
    spin0 = FNDATA['spin0']

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
                X[i,j] = (i+j)*ax+3*ax
                Y[i,j] = (j-i)*ay
                # the slanted lattice
                # X[i,j] = i*a + j*ay
                # Y[i,j] = j*ax 
   
    US = np.zeros(NTR,dtype=np.float)
    VS = np.zeros(NTR,dtype=np.float)

    #ERR = []
    QSOLX = 4.0*np.pi/3.0
    QSOLY = 0.0
    RHS = np.zeros(3,dtype=np.float)
    for i in range(0,NTR):
        RHS[0] = np.dot(spin[elt[i,0]],spin[elt[i,1]])
        RHS[1] = np.dot(spin[elt[i,0]],spin[elt[i,2]])
        RHS[2] = np.dot(spin[elt[i,1]],spin[elt[i,2]])
        SOL = minimize(qval,[QSOLX,QSOLY],
                args=(RHS),
                #bounds=((0,2*np.pi),(0,2*np.pi)),
                method='Powell',
                tol=1e-6
                )

        US[i], VS[i] = SOL.x
        #QSOLX = US[i]
        #QSOLY = VS[i]

    U = np.reshape(US,(LTR,LTR)).transpose()
    V = np.reshape(VS,(LTR,LTR)).transpose()
    
    # reprocessing to make data non-local
    Qx = np.zeros((LTR,LTR),dtype=np.float)
    Qy = np.zeros((LTR,LTR),dtype=np.float)
    for i in range(0,LTR):
        for j in range(0,LTR):
            NUM = 0
            for k in range(-SIGMA,SIGMA+1):
                for m in range(-SIGMA,SIGMA+1):
                    x = (i + k)%LTR
                    y = (j + m)%LTR
                    DIST = np.linalg.norm(np.array(k*dvec[0]+m*dvec[1]))
                    if (DIST <= 2*SIGMA*a):
                        Qx[i,j] += U[x,y]*gauss(DIST,SIGMA*a)
                        Qy[i,j] += V[x,y]*gauss(DIST,SIGMA*a)
                        NUM += gauss(DIST,SIGMA*a)
           
            Qx[i,j] *= 1.0/NUM
            Qy[i,j] *= 1.0/NUM
    
    U = Qx.flatten()
    V = Qy.flatten()

    U = (U - UZ)
    V = (V - VZ)

    UX[INDX] = np.append(UX[INDX],U)
    VX[INDX] = np.append(VX[INDX],V)
    DELTA[INDX] = float(STR_DELTA)

np.savez_compressed("../plot/QHIST.npz",
        DELTA=DELTA,
        UX=UX,
        VX=VX)


QSIGMA = np.zeros(HSHNUM,dtype=np.float)

for i in range(0,HSHNUM):

    STR_L = hshchar[i][0]
    STR_DELTA = hshchar[i][1]
    STR_ALPHA = hshchar[i][2]
    L = int(STR_L)
    N = L**2

    # plotting the histogram

    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])

    # best fit of data
    (mu, sigma) = norm.fit(UX[i])

    DELTA[i] = float(STR_DELTA)
    QSIGMA[i] = sigma

    # the histogram of the data
    ax.hist(UX[i],bins=50,normed=1,facecolor='royalblue')

    plt.title(r'Histogram of $\Delta Q_x$')
    fig.savefig("../plot/QXHIST_L_"+
            "_L_"+
            STR_L+
            "_DLT_"+
            STR_DELTA+
            "_ALP_"+
            STR_ALPHA+
            "_DNM_"+
            ".pdf",
            bbox_inches='tight',
            transparent=True
            )
    plt.close('all')

# fit function
def power(x,a):

    return a*(x**2)

# plotting the histogram

w,h = figure.figaspect(1.0)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

# the histogram of the data

ax.loglog(DELTA,QSIGMA,
        ls='None',
        marker='s',
        ms=12,
        mew=2,
        mfc='None',
        mec='red',
        basex=10,basey=10,
        zorder=1
        )

popt, pcov = curve_fit(power,DELTA,QSIGMA)

DELTAX = np.linspace(np.amin(DELTA),np.amax(DELTA),num=200)

ax.loglog(DELTAX,power(DELTAX, *popt),
        color='royalblue',
        linestyle='--',
        lw=4,
        #label=r'Power law fit, $\delta\theta\sim1/r$',
        #label=r'Exponential fit fit, $\delta\theta\sim e^{-r}$'
        #basex=10,basey=10,
        zorder=2
        #+' $n$ = '
        #+str("%.4f" % n)
        )
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
#plt.legend(loc='best')
#plt.suptitle(r"$\delta\theta$(r) vs $r$", 
#        x=0.5, fontsize=16)
plt.ylabel(r'$|\delta Q_x|$',fontsize=30)
plt.xlabel(r'$\delta J/J$',fontsize=30)
plt.xticks(rotation='vertical')
plt.tick_params(which='both',width=2,labelsize=30)
plt.tick_params(which='major',length=20)
plt.tick_params(which='minor',length=10)

fig.savefig("../plot/QXSIGMA_VS_DELTA"+
        ".pdf",
        bbox_inches='tight',
        transparent=True
        )
plt.close('all')


