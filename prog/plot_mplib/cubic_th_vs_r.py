#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
folder = "cubic_classical"
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
sfc = [0]*HSHNUM
n_sfc = np.zeros(HSHNUM)
en = np.zeros(HSHNUM)

X = np.zeros(1)
Y = np.zeros(1)
Z = np.zeros(1)

a = 1.0

# fit function
def power(x,a,b):

    return a - b*x

def invd(x,a):

    return a/x**2

def expf(x,a,b):

    return a*np.exp(-b*x)


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
    td_spin = FNDATA['spin']
    td_spin0 = FNDATA['spin0']

    N = L**3

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                X[i,j]= i*a
                Y[i,j]= j*a


    for zix in range((L+1)/2,(L+1)/2+1):
 
        ix = (zix)*(L**2)
        jx = ix + L**2

        spin = td_spin[ix:jx]
        spin0 = td_spin0[ix:jx]
        
        theta = np.arccos(np.einsum('ij,ij->i',spin,spin0))
        theta *= np.sign(np.cross(spin,spin0)[:,2])

        # plotting the configuration

        w,h = figure.figaspect(1)
        fig = plt.figure(figsize=(w,h))
        ax = fig.add_axes([0,0,1,1])
        
        shift = 0
        indx = (L*((L+1)/2-1)+L/2) + shift
        inde = indx + (L+1)/2-1 - shift
        l = np.array(range(0,(inde-indx)),
                dtype=np.float64)
        l += shift+0.5

        dtheta = np.abs(theta[indx:inde])

        indx0 = indx
        inde0 = inde
        l0 = np.array(range(0,(inde0-indx0)),
                dtype=np.float64)
        l0 += 0.5

        dtheta0 = np.abs(theta[indx0:inde0])

        ax.loglog(l0,dtheta0,
                ls='None',
                marker='s',
                ms=12,
                mew=2,
                mfc='None',
                mec='red',
                basex=10,basey=10,
                zorder=1
                )
        
        popt, pcov = curve_fit(invd,l,
                dtheta
                )

        #plotting the fit
        lp = np.linspace(0.5,(L+1.0)/2.0-1.0+0.5,num=200)
        ax.loglog(lp,invd(lp, *popt),
                color='royalblue',
                linestyle='--',
                lw=4,
                #label=r'Power law fit, $\delta\theta\sim1/r$',
                basex=10,basey=10,
                zorder=2
                )

        #plt.xlim([0.5,10])
        #plt.ylim([10.0**(-3),1.0])
        plt.ylabel(r'$\delta\theta(r)$',fontsize=30)
        plt.xlabel(r'$r$',fontsize=30)
        plt.tick_params(which='both',width=2,labelsize=30)
        plt.tick_params(which='major',length=20)
        plt.tick_params(which='minor',length=10)

        fig.savefig("../plot/cubic_th_r_"+
                STR_L+
                STR_DELTA+
                STR_ALPHA+
                "DNMR"+
                str("%06d" % IDISD)+
                str("%06d" % BTNUM)+
                str("%06d" % zix)+
                ".pdf"
                ,bbox_inches='tight'
                )
        plt.close('all')

for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('pdftk ../plot/cubic_th_r_'+
                L+
                DLT+
                ALP+
                "DNMR"+
                '* '+
                'cat output ../plot/CUBIC_TH_VS_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf',shell=True)
            subprocess.call('rm ../plot/cubic_th_r_'+
                L+
                DLT+
                ALP+
                "DNMR"+
                '* ',shell=True)

