#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = raw_input("Name of the target folder? ")
folder = "single_impurity"
os.chdir("../"+folder+"/out/data")
import sys
import subprocess
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.figure as figure
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit


## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# the font styleset
#from matplotlib import rcParams
#rcParams['font.serif'] = ['Times New Roman']
#rcParams['font.family'] = 'serif'

# setting MKL configuration to allways optimise 
os.putenv("MKL_DYNAMIC","FALSE")

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*)_DISD_([^/]*)_BOOT_([^/]*).npz')

if not os.path.isdir("../plot"):
    os.mkdir("../plot")

DHASH = set([])
LHASH = set([])

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))

    HSH = str(DELTA)
    DHASH.add(HSH)

    HSH = str(L)
    LHASH.add(HSH)

DHASH = list(DHASH)
DHSHNUM = len(DHASH)
LHASH = list(LHASH)
LHSHNUM = len(LHASH)

DLTBNK = np.zeros((DHSHNUM,LHSHNUM),dtype=np.float)

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))
    ALPHA = float(match.group(3))

    DHSHCHK = str(DELTA)
    LHSHCHK = str(L)
    DINDX = DHASH.index(DHSHCHK)
    LINDX = LHASH.index(LHSHCHK)

    FNDATA = np.load(fname)

    spin = FNDATA['spin']
    
    totspin = np.sum(spin,axis=0)
    DLTBNK[DINDX,LINDX] = np.linalg.norm(totspin)

# fit function
def power(x,a,b):

    return a*x + b

def invd(x,a,b):

    return a/x + b

def expf(x,a,b):

    return a*np.exp(-b*x)


DLT = np.zeros(DHSHNUM)
TSPIN = np.zeros(DHSHNUM)

for i in range(0,DHSHNUM):
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])
   
    d = float(DHASH[i])
    l = np.array(LHASH).astype(np.float)
    dtheta = DLTBNK[i]
    
    """
    ax.loglog(l,dtheta,'bx',
            ms=10,
            mew=2,
            basex=10,basey=10
            )
    """
    ax.plot(l,dtheta,'bx',
            ms=10,
            mew=2,
            )

    
    popt, pcov = curve_fit(power,l,dtheta)

    lp = np.linspace(np.amin(l),np.amax(l),num=200)
    """    
    ax.loglog(lp,power(lp, *popt),
            'r--',
            lw=2,
            label=r'\Delta='+str("%.4f" % d),
            basex=10,basey=10
            )
    """
    ax.plot(lp,power(lp, *popt),
            'r--',
            lw=2,
            label=r'\Delta='+str("%.4f" % d),
            )


    tspin = power(10*(10),*popt)

    DLT[i] = d
    TSPIN[i] = tspin
   
    plt.title(r'$\Delta$='+str("%.4f" % d),fontsize=20)
    plt.ylabel(r'$|\vec{S}_{\textrm{tot}}|$',fontsize=20)
    plt.xlabel(r'$L$',fontsize=20)
    plt.tick_params(which='both',width=2,labelsize=16)
    plt.tick_params(which='major',length=8)
    plt.tick_params(which='minor',length=4)

    #plt.ylim([10.0**(-2),0.5])
    plt.ylim([np.amin(dtheta),np.amax(dtheta)])
    
    fig.savefig("../plot/fract-"
            +"dlt_"
            +str("%.4f" % d)
            +".pdf",
            bbox_inches='tight'
            )

    plt.close('all')

subprocess.call('pdftk ../plot/fract-dlt_*'
        +' cat output ../plot/FRACT_VS_L.pdf',
        shell=True)

subprocess.call('rm ../plot/fract-dlt_*',
        shell=True)

ARG = np.argsort(DLT)
TSPIN = TSPIN[ARG]
DLT = np.sort(DLT)

w,h = figure.figaspect(1.0)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

ax.plot(DLT,TSPIN,'bx',
        ms=10,
        mew=2
        )

popt, pcov = curve_fit(power,DLT,TSPIN)

dp = np.linspace(np.amin(DLT),np.amax(DLT),num=200)

ax.plot(dp,power(dp, *popt),
        'r--',
        lw=2,
        label=r'\Delta='+str("%.4f" % d),
        )

plt.ylim([10.0**(-2),0.5])
plt.ylabel(r'$|\vec{S}_{\textrm{tot}}|$',fontsize=20)
plt.xlabel(r'$\Delta$',fontsize=20)
plt.tick_params(which='both',width=2,labelsize=16)
plt.tick_params(which='major',length=8)
plt.tick_params(which='minor',length=4)

fig.savefig("../plot/DELTA_VS_FRACT.pdf",
        bbox_inches='tight'
        )
plt.close('all')

