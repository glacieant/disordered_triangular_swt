#!/usr/bin/env python

import os
import datetime
import multiprocessing as mp
#folder = "xy_disorder"
folder = "xy_gaussian"
#folder = "hbrg_gaussian"
os.chdir("../"+folder+"/out/data")
#os.chdir("../"+folder+"/out/data_18.07.18")
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
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

# TU colors

tu_dunkelblau = '#07284A'
tu_grau = '#5F6967'
tu_blau = '#0067A5'
tu_cyan = '#00A3DA'
tu_dunkelgruen = '#008644'
tu_gruen = '#5EB245'
tu_rot = '#DE4A39'
tu_dunkelrot = '#BD252C'


## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

# setting MKL configuration to allways optimise 
os.putenv("MKL_DYNAMIC","FALSE")

# importing lattice geometry system

import lattice_map as lmap

## picturing the output data

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*)_DISD_([^/]*)_BOOT_([^/]*)_DNM_([^/]*).npz')

if not os.path.isdir("../plot"):
    os.mkdir("../plot")

HASH = set([])
LHASH = set([])
DHASH =set([])
AHASH = set([])

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))
    ALPHA = float(match.group(3))

    HSH = str(L)+str(DELTA)+str(ALPHA)
    HASH.add(HSH)
    LHASH.add(str(L))
    DHASH.add(str(DELTA))
    AHASH.add(str(ALPHA))

HASH = list(HASH)

# the master arrays of data
HSHNUM = len(HASH)
hshchar = [0]*HSHNUM
SPC = [0]*HSHNUM
EFC = [0]*HSHNUM
EFC_ERR = [0]*HSHNUM
N_EFC = np.zeros(HSHNUM)

# some arrays for the scaling analysis

TLRAY = 1.0/np.array(list(LHASH)).astype(float)
LRAY = np.sort(TLRAY)
LORD = np.argsort(TLRAY)

TDRAY = np.array(list(DHASH)).astype(float)
DRAY = np.sort(TDRAY)
DORD = np.argsort(TDRAY)

LNUM = len(LHASH)
DNUM = len(DHASH)
ANUM = len(AHASH)

ENRGY = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
SCOM = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
ENRGY_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)

# Defining model function to be used to fit 1d data

def linear(x,a,b):

    return a-b*x


# reciprocal lattice vectors
b1 = np.array([2.0*np.pi,-2.0*np.pi/3.0**0.5])
b2 = np.array([0.0,4.0*np.pi/3.0**0.5])

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))
    ALPHA = float(match.group(3))

    HSHCHK = str(L)+str(DELTA)+str(ALPHA)
    INDX = HASH.index(HSHCHK)
    hshchar[INDX] = [str(L),str(DELTA),str(ALPHA)]

    FNDATA = np.load(fname)

    spin = FNDATA['spin']
    J = FNDATA['J']

    N = L**2

    NBR = np.zeros((N,6),dtype=np.int)

    lmap.lattice_map(L,NBR)

    ENCL = 0.0

    for i in range(0,N):

        for j in NBR[i]:

            ENCL += 0.5*J[i,j]*np.dot(spin[i],spin[j])

    #ENCL = 0.5*np.einsum("ij,ik,jk",J,spin,spin)/N

    print ENCL

    N_EFC[INDX] += 1

    SPC[INDX] = spin[0].size
    EFC[INDX] += ENCL
    EFC_ERR[INDX] += ENCL*ENCL

LHASH = list(LHASH)
DHASH = list(DHASH)
AHASH = list(AHASH)

for i in range(0,HSHNUM):

    STR_L = hshchar[i][0]
    STR_DELTA = hshchar[i][1]
    STR_ALPHA = hshchar[i][2]

    EFX = EFC[i]/N_EFC[i]
    EFX_ERR = np.sqrt((EFC_ERR[i]/N_EFC[i]-SFX**2)/N_EFC[i])

    LN = LHASH.index(STR_L)
    DN = DHASH.index(STR_DELTA)
    AN = AHASH.index(STR_ALPHA)

    SCOM[LN,DN,AN] = SPC[i]
    ENRGY[LN,DN,AN] = EFX
    ENRGY_ERR[LN,DN,AN] = EFX_ERR


for AN in range(0,ANUM):

    w,h = figure.figaspect(1.0)
    afig = plt.figure(figsize=(w,h))
    ax = afig.add_axes([0.26,0.15,0.685,0.8])

    for LN in range(0,LNUM):

        Y = ENRGY[LN,:,AN]
        Y_ERR = ENERGY_ERR[LN,:,AN]
        X = np.array(DHASH).astype(float)

        line = ax.errorbar(X,
                Y,
                yerr=Y_ERR,
                lw=2,
                marker='.',
                ms=10,
                label=r'$L = $'+LHASH[LN]
                )

        DHLIST = [float(DHASH[i]) for i in range(0,DNUM)]
        handles1, labels1 = ax.get_legend_handles_labels()

        tup = sorted(zip(DHLIST,handles1,labels1))
        DHLIST, handles1, labels1 = zip(*tup)

        legend1 = ax.legend(handles1,labels1,
                loc='best',
                fontsize=12,
                markerscale=1,
                facecolor='w',
                edgecolor='k',
                framealpha=1)

    if SCOM[0,0,AN] == 3:
        ax.set_title('Antiferromagnetic XY')
    else:
        ax.set_title('Antiferromagnetic Heisenberg')

    ax.set_ylabel(r'$E$',fontsize=30)
    ax.set_xlabel(r'$\delta J/J$',fontsize=30)
        
    ax.tick_params(which='both',width=2,
            labelsize=30,direction='in',
            bottom=True,top=True,
            left=True,right=True)
    ax.tick_params(which='major',length=20)
    ax.tick_params(which='minor',length=10)
    
    afig.savefig("../plot/E-VS-DELTA"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )
    
    plt.close('all')

