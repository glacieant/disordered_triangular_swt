#!/usr/bin/env python

import os
import datetime
import multiprocessing as mp
folder = "xy_disorder"
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
sfc = [0]*HSHNUM
sfc_var = [0]*HSHNUM
n_sfc = np.zeros(HSHNUM)

# some arrays for the scaling analysis

TLRAY = 1.0/np.array(list(LHASH)).astype(float)
INVLRAY = np.sort(TLRAY)
LORD = np.argsort(TLRAY)

TDRAY = np.array(list(DHASH)).astype(float)
DRAY = np.sort(TDRAY)
DORD = np.argsort(TDRAY)

LNUM = len(LHASH)
DNUM = len(DHASH)
ANUM = len(AHASH)

SFC_MAX = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
SFC_MAX_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)

SFC_DELTA = np.zeros((DNUM,ANUM),dtype=np.float)
DELTA_ARR = np.zeros((DNUM,ANUM),dtype=np.float)


# Defining model function to be used to fit 1d data

def linear(x,a,b):

    return a + b*x

def res_log(invL,c1,c2,c3):

    return c1 + c2*invL + c3*invL*np.log(invL)

def power_law(x,a,b):

    return a*(x**b)

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

    
    spin_X = spin[:,0].reshape((L,L))
    spin_Y = spin[:,1].reshape((L,L))
    SFC_X = np.abs(np.fft.fft2(spin_X,norm='ortho'))**2
    SFC_Y = np.abs(np.fft.fft2(spin_Y,norm='ortho'))**2

    SFC = (SFC_X + SFC_Y)
    
    # Only considering Q = (4*pi/3,0)
    KMAX_X = L/3 
    KMAX_Y = 2*L/3

    SFM = SFC[KMAX_X,KMAX_Y]/N 

    #SFM = np.linalg.norm(np.sum(spin,axis=0))/N

    n_sfc[INDX] += 1
    
    sfc[INDX] += SFM
    sfc_var[INDX] += SFM*SFM


LHASH = list(LHASH)
DHASH = list(DHASH)
AHASH = list(AHASH)

for i in range(0,HSHNUM):

    STR_L = hshchar[i][0]
    STR_DELTA = hshchar[i][1]
    STR_ALPHA = hshchar[i][2]
    L = int(STR_L)
    N = L**2

    SFX = sfc[i]/n_sfc[i]
    SFX_ERR = np.sqrt((sfc_var[i]/n_sfc[i]-SFX**2)/n_sfc[i])

    # scaling data

    LN = LHASH.index(STR_L)
    DN = DHASH.index(STR_DELTA)
    AN = AHASH.index(STR_ALPHA)

    SFC_MAX[LN,DN,AN] = SFX
    SFC_MAX_ERR[LN,DN,AN] = SFX_ERR


for AN in range(0,ANUM):

    w,h = figure.figaspect(1.0)
    afig = plt.figure(figsize=(w,h))
    ax = afig.add_axes([0.26,0.15,0.685,0.8])

    for DN in range(0,DNUM):

        INVNRAY = INVLRAY**2

        popt, pcov = curve_fit(power_law,
                INVLRAY,
                SFC_MAX[:,DN,AN][LORD]
                )

        YDATA = SFC_MAX[:,DN,AN][LORD]
        YDATA_ERR = SFC_MAX_ERR[:,DN,AN][LORD]
        XDATA = INVLRAY

        ax.errorbar(XDATA,
                YDATA,
                yerr=YDATA_ERR,
                ls='None',
                marker='.',
                ms=10,
                )

        popt, pcov = curve_fit(linear,XDATA,YDATA)

        PDATA = np.linspace(XDATA.min(),
                XDATA.max(),
                200)
        color = ax.get_lines()[-1].get_color()
        ax.plot(PDATA,linear(PDATA, *popt),
                ls='--',
                lw=2,
                color=color,
                label=r'$\Delta=$'+DHASH[DN]+r' $\alpha=$'+str(popt[1])
                )

    DHLIST = [float(DHASH[i]) for i in range(0,DNUM)]
    handles, labels = ax.get_legend_handles_labels()

    tup = sorted(zip(DHLIST,handles,labels))
    DHLIST, handles, labels = zip(*tup)
    
    """
    legend = ax.legend(handles,labels,
            loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)
            
    """

    #ax.set_xlim(left=0.0,right=0.25)
    #ax.set_ylim(bottom=0.0,top=0.01)
 
    #ax.set_title(r'$\alpha=$'+str("%.2f" % float(AHASH[AN])))
    ax.set_ylabel(r'$S(Q)/L^2$',fontsize=20)
    ax.set_xlabel(r'$1/L$',fontsize=20)
    #ax.set_xlabel(r'$1/L$',fontsize=20)
    #ax.set_xlabel(r'$|(1/L)\log(1/L)|$',fontsize=20)

    afig.savefig("../plot/SFC_DELTA"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )

    plt.close('all')

