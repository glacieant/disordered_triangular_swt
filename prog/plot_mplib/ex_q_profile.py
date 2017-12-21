#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
folder = "sharp_wall"
os.chdir("../"+folder+"/out/data")
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
sfc = [0]*HSHNUM
n_sfc = np.zeros(HSHNUM)
en = np.zeros(HSHNUM)

# importing lattice geometry system
import lattice_map as lmap
ZCO = 12

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

    n_sfc[INDX] += 1

    spin = FNDATA['spin']

    N = L**2
    nbr = np.zeros((N,ZCO),dtype=np.int)
    lmap.lattice_map(L,nbr)

    qnn = np.zeros(N,dtype=np.float64)

    for i in range(0,N):

        qnn[i] = 0.0

        for j in range(0,6):

            qnn[i] += np.dot(spin[i],spin[nbr[i,j]])

        qnn[i] = qnn[i]/6

    sfc[INDX] += qnn


# reciprocal lattice vectors
b1 = np.array([2.0*np.pi,-2.0*np.pi/3.0**0.5])
b2 = np.array([0.0,4.0*np.pi/3.0**0.5])

for i in range(0,HSHNUM):

    if n_sfc[i] > 0:

        fig, ax = plt.subplots()


        STR_L = hshchar[i][0]
        STR_DELTA = hshchar[i][1]
        STR_ALPHA = hshchar[i][2]
        L = int(STR_L)
        N = L**2

        sfc[i] *= (1.0/n_sfc[i])
        sfc[i] = np.arccos(sfc[i])
        en[i] *= (1.0/n_sfc[i])

        # setting up the clipping of the density plot
        MAXCLIP = np.pi
        MINCLIP = 0.0

        cax = ax.imshow(sfc[i].reshape(L,L),
                norm=Normalize(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                interpolation='nearest',
                extent=(0.0,1.0,0.0,1.0),
                cmap=cm.Accent,aspect='auto')
        cbar = fig.colorbar(cax,shrink=0.5,format='%.0e')
        plt.suptitle(r"Q Profile, $\cos^{-1}(1/6\sum_{\hat{r}}\vec{S}_i\cdot\vec{S}_{i+\hat{r}})$", 
                x=0.5,y=1.0, fontsize=16)
        """
        plt.title(r"$\Delta$ = "+STR_DELTA+" , "
                +r"$\alpha$ = "+STR_ALPHA+" , "
                +"|E|/$L^{2}$ = "+str("%f" % en[i]),
                x=0.6,fontsize=12)
        """     
        smax = np.amax(sfc[i])/np.pi
        smin = np.amin(sfc[i])/np.pi
        plt.figtext(.8,.85,
                r'$Q^{\mathrm{max}}$ = '+str("%f" % smax)
                +r' $ \pi $',
                fontsize=12)
        plt.figtext(.8,.75,
                r'$Q^{\mathrm{min}}$ = '+str("%f" % smin)
                +r' $ \pi $',
                fontsize=12)
        plt.xlabel('L = '+STR_L,fontsize=16)
        fig.tight_layout(pad=2.5,h_pad=2.5,w_pad=2.5)
        fig.savefig("../plot/qprof_"+
                STR_L+
                STR_DELTA+
                STR_ALPHA+
                ".pdf"
                )

        plt.close('all')

for L in zip(*hshchar)[0]:
    subprocess.call('pdftk ../plot/qprof_'+
            L+
            '* '+
            'cat output ../plot/QPROFILE_L_'+
            L+
            '.pdf',shell=True)

# removing split files
subprocess.call('rm ../plot/qprof_*',shell=True)
