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
from scipy.optimize import curve_fit


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

    #TOL = FNDATA['err'][()]

    #if TOL < 101:

    n_sfc[INDX] += 1

    theta = FNDATA['theta']
    sfc[INDX] += theta

# fit function
def power(x,a,b):

    return a*x**(-b)

def invd(x,a):

    return a/x


for i in range(0,HSHNUM):

    if n_sfc[i] > 0:

        fig, ax = plt.subplots()


        STR_L = hshchar[i][0]
        STR_DELTA = hshchar[i][1]
        STR_ALPHA = hshchar[i][2]
        L = int(STR_L)
        N = L**2

        sfc[i] *= (1.0/(n_sfc[i]*N))
        
        indx = (N+1)/2
        inde = indx+(L+1)/2-1
        l = np.arange(1,(L+1)/2-1)
    
        ax.plot(l,sfc[i][indx+1:inde],'ro')
        
        #fitting data
        #popt, pcov = curve_fit(power,l,
        #        sfc[i][indx+1:inde],
        #        p0=[sfc[i][indx+1]/l[0],1.05])
        popt, pcov = curve_fit(invd,l,
                sfc[i][indx+1:inde],
                p0=[sfc[i][indx+1]/l[0]])

        #n = popt[0]
        n = 1
        #plotting the fit
        ax.plot(l,invd(l, *popt),'r-',
                label=r'Power law fit, $\theta\sim1/r^n$, with'
                +' $n$ = '
                +str("%.4f" % n))
        plt.legend(loc='best')
        plt.suptitle(r"$\theta$(r) vs $r$", 
                x=0.5, fontsize=16)
        plt.title(r"$\Delta$ = "+STR_DELTA+" , "
                +r"$\alpha$ = "+STR_ALPHA+" , "
                +"|E|/$L^{2}$ = "+str("%f" % en[i]),
                x=0.6,fontsize=12)
        plt.xlabel('r = '+STR_L,fontsize=16)
        fig.tight_layout(pad=2.5,h_pad=2.5,w_pad=2.5)
        fig.savefig("../plot/tht_"+
                STR_L+
                STR_DELTA+
                STR_ALPHA+
                ".pdf"
                )

        plt.close('all')

for L in zip(*hshchar)[0]:
    subprocess.call('pdftk ../plot/tht_'+
            L+
            '* '+
            'cat output ../plot/THT_L_'+
            L+
            '.pdf',shell=True)

# removing split files
subprocess.call('rm ../plot/tht_*',shell=True)
