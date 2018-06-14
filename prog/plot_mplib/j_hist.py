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

## picturing the output data

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*)_DISD_([^/]*)_BOOT_([^/]*)_DNM_([^/]*).npz')

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
INDXCOUNT = [0]*HSHNUM
DELTA = [0]*HSHNUM

import lattice_map as lmap

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

    N = L**2

    nbr = np.zeros((N,12),dtype=np.int)
    lmap.lattice_map(L,nbr)

    JX = np.zeros(N*3)

    o = 0

    for i in range(0,N):

        for k in range(0,6):

            j = nbr[i,k]

            if j > i:
                JX[o] = J[i,j]
                o += 1
    
    U = JX

    if INDXCOUNT[INDX] == 0:
        UX[INDX] = U
    else:
        UX[INDX] = np.append(UX[INDX],U)
    INDXCOUNT[INDX] += 1
    print UX[INDX].min()
    DELTA[INDX] = float(STR_DELTA)

np.savez_compressed("../plot/JHIST.npz",
        DELTA=DELTA,
        UX=UX)

JSIGMA = np.zeros(HSHNUM,dtype=np.float)

for i in range(0,HSHNUM):
    
    DLT = DELTA[i]

    # plotting the histogram

    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])

    # best fit of data
    (mu, sigma) = norm.fit(UX[i])

    JSIGMA[i] = sigma

    # the histogram of the data
    ax.hist(UX[i],bins=50,normed=1,facecolor='royalblue')

    plt.title(r'Histogram of $J$')

    fig.savefig("../plot/JHIST"+
            "_DLT_"+
            str("%.2f" % DLT)+
            "_DNM_"+
            ".pdf",
            bbox_inches='tight',
            transparent=True
            )
    plt.close('all')
