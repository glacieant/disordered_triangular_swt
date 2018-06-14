#!/usr/bin/env python

import os
#folder = "sharp_wall"
#folder = "single_impurity"
folder = "zero_field_classical"
os.chdir("../"+folder+"/out/data_20.04.18")
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

## picturing the output data

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*)_DISD_([^/]*)_BOOT_([^/]*).npz')

if not os.path.isdir("../plot"):
    os.mkdir("../plot")

HASH = set([])

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))
    ALPHA = float(match.group(3))

    HSH = str(L)+str(DELTA)+str(ALPHA)
    HASH.add(HSH)

HASH = list(HASH)

# the master arrays of data
HSHNUM = len(HASH)
hshchar = [0]*HSHNUM
sfc = [0]*HSHNUM
sfc_var = [0]*HSHNUM
n_sfc = np.zeros(HSHNUM)

# Defining model function to be used to fit 1d data

def gauss(x,x0,sigma):

    return ((1.0/(sigma*np.sqrt(2*np.pi)))
            *np.exp(-((x-x0)**2)/(2*(sigma**2))))


# reciprocal lattice vectors
b1 = np.array([2.0*np.pi,-2.0*np.pi/3.0**0.5])
b2 = np.array([0.0,4.0*np.pi/3.0**0.5])

import lattice_map as lmap

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

    LTR = L-1
    NTR = LTR**2

    elt = np.zeros((NTR,3),dtype=np.int)
    lmap.eltriangle(L,elt)

    tsr_spin = spin[elt]

    epsilon = np.zeros((3,3,3),dtype=np.float)
    epsilon[0,1,2] = epsilon[1,2,0] = epsilon[2,0,1] = 1.0
    epsilon[0,2,1] = epsilon[2,1,0] = epsilon[1,0,2] = -1.0

    alt = np.zeros((3,3),dtype=np.float)
    alt[0,1] = alt[1,2] = alt[2,0] = 1.0
    
    chir = 2.0*np.einsum('iab,icd,ac,kbd->ik',
            tsr_spin,tsr_spin,alt,epsilon)/(3.0*np.sqrt(3.0))

    chir_X = chir[:,0].reshape((LTR,LTR))
    chir_Y = chir[:,1].reshape((LTR,LTR))
    chir_Z = chir[:,2].reshape((LTR,LTR))
    fchir_X = np.abs(np.fft.fft2(chir_X,norm='ortho'))**2
    fchir_Y = np.abs(np.fft.fft2(chir_Y,norm='ortho'))**2
    fchir_Z = np.abs(np.fft.fft2(chir_Z,norm='ortho'))**2

    FCHIR = fchir_X + fchir_Y + fchir_Z

    KMAX_X = 0
    KMAX_Y = 0

    SFM = FCHIR[:,KMAX_Y] 
    
    n_sfc[INDX] += 1
    
    sfc[INDX] += SFM
    #sfc_var[INDX] += SFM*SFM


for i in range(0,HSHNUM):

    STR_L = hshchar[i][0]
    STR_DELTA = hshchar[i][1]
    STR_ALPHA = hshchar[i][2]

    SFX = sfc[i]/n_sfc[i]
    
    # fitting a gaussian to extract a pattern along 
    # the high symmetry point GAMMA

    fig, ax = plt.subplots()

    # data for gaussian fitting
    YDATA = SFX
    XDATA = np.arange(YDATA.size)
    ax.plot(XDATA,YDATA,'b.',label='Data')

    # finding the gauss fit

    #popt, pcov = curve_fit(gauss,XDATA,YDATA)

    #a = popt[0]
    #b = popt[1]
    #c = popt[2]

    #PDATA = np.linspace(XDATA.min(),XDATA.max(),
    #        200)
    #ax.plot(PDATA,gauss(PDATA, *popt),'k-')

    plt.legend(loc='best')

    fig.savefig("../plot/XI_HIST"+
            "_L_"+
            STR_L+
            "_DLT_"+
            STR_DELTA+
            "_ALP_"+
            STR_ALPHA+
            ".pdf",
            bbox_inches='tight',
            transparent=True
            )
    plt.close('all')

