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

FNDATA = np.load("../plot/QHIST.npz")

DELTA = FNDATA['DELTA']
UX = FNDATA['UX']
VX = FNDATA['VX']

HSHNUM = len(DELTA) 

QSIGMA = np.zeros(HSHNUM,dtype=np.float)

for i in range(0,HSHNUM):
    
    """
    STR_L = hshchar[i][0]
    STR_DELTA = hshchar[i][1]
    STR_ALPHA = hshchar[i][2]
    L = int(STR_L)
    N = L**2
    """
    
    DLT = DELTA[i]

    # plotting the histogram

    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])

    # best fit of data
    (mu, sigma) = norm.fit(UX[i])

    QSIGMA[i] = sigma

    # the histogram of the data
    ax.hist(UX[i],bins=50,normed=1,facecolor='royalblue')

    plt.title(r'Histogram of $\Delta Q_x$')
    """
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
    """
    fig.savefig("../plot/QXHIST"+
            "_DLT_"+
            str("%.2f" % DLT)+
            "_DNM_"+
            ".pdf",
            bbox_inches='tight',
            transparent=True
            )
    plt.close('all')

# fit function
def power(x,a):

    return a*x

# plotting the histogram

w,h = figure.figaspect(1.0)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

# the histogram of the data

ax.plot(DELTA,QSIGMA,
        ls='None',
        marker='s',
        ms=12,
        mew=2,
        mfc='None',
        mec='red',
        #basex=10,basey=10,
        zorder=1
        )

popt, pcov = curve_fit(power,DELTA[:3],QSIGMA[:3])

DELTAX = np.linspace(np.amin(DELTA),np.amax(DELTA),num=200)

ax.plot(DELTAX,power(DELTAX, *popt),
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


