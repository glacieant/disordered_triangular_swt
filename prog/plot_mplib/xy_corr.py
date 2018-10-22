#!/usr/bin/env python

import os
import datetime
import multiprocessing as mp
folder = "xy_gaussian"
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
SFC = [0]*HSHNUM
SFC_ERR = [0]*HSHNUM
RFC = [0]*HSHNUM
N_SFC = np.zeros(HSHNUM)

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

ITNUM = np.zeros((LNUM,DNUM,ANUM),dtype=np.int)

# Defining model function to be used to fit 1d data

def gaussian(x,a,b,c):

    return a*np.exp(-(x-b)**2/(2.0*c**2))

def parabola(x,a,b,c):

    return a*(x**2)+b*x+c

def skew_gauss(x,alpha,a,b,c):

    return a*skewnorm.pdf(x,alpha,b,np.sqrt(2.0)*c)

def linear(x,a,b):

    return a + b*x

def invpar(x,a,b):

    return a + b/(x**2)


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

    N = L**2

    sbltc = np.zeros(N)
    lmap.sublattice_map(L,sbltc)

    SFM = []
    RFM = []

    for i in range(0,L/2):
        if sbltc[i] == 0:
            SFM.append(np.dot(spin[0],spin[i]))
            RFM.append(i)

    N_SFC[INDX] += 1

    NPSFM = np.array(SFM)
    NPRFM = np.array(RFM)
    
    SFC[INDX] += NPSFM
    SFC_ERR[INDX] += NPSFM*NPSFM
    RFC[INDX] += NPRFM

LHASH = list(LHASH)
DHASH = list(DHASH)
AHASH = list(AHASH)

for i in range(0,HSHNUM):

    STR_L = hshchar[i][0]
    STR_DELTA = hshchar[i][1]
    STR_ALPHA = hshchar[i][2]
    L = int(STR_L)
    N = L**2

    SFX = SFC[i]/N_SFC[i]
    SFX_ERR = np.sqrt((SFC_ERR[i]/N_SFC[i]-SFX**2)/N_SFC[i])
    RFX = RFC[i]/N_SFC[i]

    SFC[i] = SFX
    SFC_ERR[i] = SFX_ERR
    RFC[i] = RFX

    LN = LHASH.index(STR_L)
    DN = DHASH.index(STR_DELTA)
    AN = AHASH.index(STR_ALPHA)

    ITNUM[LN,DN,AN] = i

CLR = [tu_rot,tu_blau,tu_gruen]

for LN in range(0,LNUM):

    for AN in range(0,ANUM):

        w,h = figure.figaspect(1.0)
        afig = plt.figure(figsize=(w,h))
        ax = afig.add_axes([0.26,0.15,0.685,0.8])
        
        for DN in range(0,DNUM):

            ix = ITNUM[LN,DN,AN]

            line1 = ax.errorbar(1.0/RFC[ix][1:],
                    SFC[ix][1:],
                    yerr=SFC_ERR[ix][1:],
                    lw=2,
                    marker='.',
                    ms=10,
                    label=r'$\Delta=$'+DHASH[DN],
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

        #time = str(datetime.datetime.now()) 
        #itmin = str(np.amin(ITNUM[:,:,AN])) 
        
        #ax.set_title('DATE = '
        #        +time
        #        +' MINIMUM REALIZATION = '
        #        +itmin
        #        )

        ax.set_title('Antiferromagnetic Heisenberg')
        #ax.set_title('Ferromagnetic Heisenberg')

        ax.set_xlim(left=0.0)
        #ax.set_ylim(bottom=0.0)

        ax.set_ylabel(r'$C(r)$',fontsize=30)
        ax.set_xlabel(r'$1/r$',fontsize=30)
            
        ax.tick_params(which='both',width=2,
                labelsize=30,direction='in',
                bottom=True,top=True,
                left=True,right=True)
        ax.tick_params(which='major',length=20)
        ax.tick_params(which='minor',length=10)
        
        afig.savefig("../plot/CORR-HBRG"
                +"_L_"
                +str("%d" % int(LHASH[LN]))
                +"_ALPHA_"
                +str("%.4f" % float(AHASH[AN]))
                +".pdf",
                bbox_inches='tight'
                )
        
        plt.close('all')

