#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
os.chdir("../out/data")
import sys
import subprocess
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

## The tex style commands
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')

# the font styleset
from matplotlib import rcParams
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.family'] = 'serif'

# setting MKL configuration to allways optimise 
os.putenv("MKL_DYNAMIC","FALSE")

# the finite temperature fermi function
def fermi(E,T):

    return 1.0/(np.exp(E/T)+1.0)

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
schir = [0]*HSHNUM
vchir = [0]*HSHNUM
n_vchir = np.zeros(HSHNUM)
en = np.zeros(HSHNUM)

# defining a tolerance limit for calculations
zero = 10.0**(-8)

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))
    ALPHA = float(match.group(3))

    HSHCHK = str(L)+str(DELTA)+str(ALPHA)
    INDX = HASH.index(HSHCHK)
    hshchar[INDX] = [str(L),str(DELTA),str(ALPHA)]

    FNDATA = np.load(fname)

    TOL = FNDATA['err'][()]

    if TOL < 101:

        n_vchir[INDX] += 1

        T = FNDATA['T'][()]
        ensys = FNDATA['ensys']
        lmult = FNDATA['lmult']
        spin = FNDATA['spin']

        N = L**2

        # creating elementary triangle indices
        TNUM = int((L-1)**2)
        eltri = np.zeros((TNUM,3),dtype=np.int32)

        for i in range(0,TNUM):
            j = i + i/(L-1)
            eltri[i,0] = j
            eltri[i,1] = j+1
            eltri[i,2] = j+L

        en[INDX] += 0.5*(np.sum(fermi(ensys,T)*ensys)
                -np.sum(lmult))/N

        # calculating vector chirality
        
        vec = np.array([0.0,0.0,0.0])
        scal = 0.0
        for i in range(0,TNUM):
            vec += np.cross(spin[eltri[i,0],:],spin[eltri[i,1],:])
            vec += np.cross(spin[eltri[i,1],:],spin[eltri[i,2],:])
            vec += np.cross(spin[eltri[i,2],:],spin[eltri[i,0],:])
            scal += np.absolute(np.dot(spin[eltri[i,0],:],
                    np.cross(spin[eltri[i,1],:],
                        spin[eltri[i,2],:])))
        schir[INDX] += (1.0/TNUM)*scal
        vchir[INDX] += (1.0/(TNUM*3*0.5*3.0**0.5))*np.linalg.norm(vec)


for i in range(0,HSHNUM):

    if n_vchir[i] > 0:

        STR_L = hshchar[i][0]
        STR_DELTA = hshchar[i][1]
        STR_ALPHA = hshchar[i][2]
        L = int(STR_L)
        DLT = float(STR_DELTA)
        ALP = float(STR_ALPHA)
        schir[i]=(1.0/n_vchir[i])*schir[i]
        vchir[i]=(1.0/n_vchir[i])*vchir[i]

        print L,DLT,ALP,schir[i],vchir[i]


