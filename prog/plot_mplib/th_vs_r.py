#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = raw_input("Name of the target folder? ")
folder = "single_impurity"
os.chdir("../"+folder+"/out/data_bk")
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

    spin = FNDATA['spin']
    spin0 = FNDATA['spin0']

    
    # defining the plane for the spins
    
    N = L**2
    indx = N - 1
    
    """
    e1 = spin[indx]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    U = np.einsum('ij,j->i',spin,e1)
    V = np.einsum('ij,j->i',spin,e2)
    
    e1 = spin0[indx]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin0[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    U0 = np.einsum('ij,j->i',spin0,e1)
    V0 = np.einsum('ij,j->i',spin0,e2)

    theta = np.arccos(U*U0 + V*V0)
    signm = np.sign(np.arcsin(U0*V - U*V0))
    theta *= signm
    
    
    
    # rotating spins to preferred direction
    ROT = np.zeros((3,3),dtype=np.float64)
    for i in range(0,3):
        for j in range(0,3):
            ROT[i,j] = spin0[0,i]*spin[0,j]

    r_spin = np.einsum("ij,kj->ki",ROT,spin)
    theta = np.arccos(np.einsum('ij,ij->i',r_spin,spin0))
    #theta *= np.sign(np.cross(r_spin,spin0)[:,2])
    #theta = FNDATA['theta']
    """    
    theta = np.arccos(np.einsum('ij,ij->i',spin,spin0))
    theta *= np.sign(np.cross(spin,spin0)[:,2])
    sfc[INDX] += theta

# fit function
def power(x,a,b):

    return a - b*x

def invd(x,a):

    return a/x

def expf(x,a,b):

    return a*np.exp(-b*x)

for i in range(0,HSHNUM):

    if n_sfc[i] > 0:

        #fig, ax = plt.subplots()

        STR_L = hshchar[i][0]
        STR_DELTA = hshchar[i][1]
        STR_ALPHA = hshchar[i][2]
        L = int(STR_L)
        N = L**2
        
        sfc[i] *= 1.0/n_sfc[i]
        """
        cax = ax.imshow(np.reshape(sfc[i],(L,L)),
                origin='lower',
                norm=Normalize(vmin=-0.4,vmax=0.4,clip=False),
                interpolation='none',
                cmap=cm.viridis,
                aspect='auto')
        cbar = fig.colorbar(cax, shrink=0.5)
        plt.suptitle(r'Numerical Simulation',x=0.5,fontsize=16)
        ax.set_aspect('equal')
        fig.tight_layout(pad=2.5,h_pad=2.5,w_pad=2.5)
        fig.savefig("../plot/dnsty_"+
                STR_L+
                STR_DELTA+
                STR_ALPHA+
                ".pdf"
                )
        
        plt.close('all')
        """
        #fig, ax = plt.subplots()
        w,h = figure.figaspect(1.0)
        fig = plt.figure(figsize=(w,h))
        ax = fig.add_axes([0,0,1,1])
        
        shift = 2
        indx = (L*(L-1))/2+shift
        inde = indx+(L+1)/2-1
        l = np.array(range(0,(inde-indx)),
                dtype=np.float64)
        l += shift+0.5

        dtheta = np.abs(sfc[i][indx:inde])

        indx0 = (L*(L-1))/2
        inde0 = indx0 + 50 #+(L+1)/2-1
        l0 = np.array(range(0,(inde0-indx0)),
                dtype=np.float64)
        l0 += 0.5

        dtheta0 = np.abs(sfc[i][indx0:inde0])

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
        
        #fitting data
        #popt, pcov = curve_fit(expf,l,
        #        sfc[i][indx:inde],
        #        p0=[sfc[i][indx]*np.exp(l[0]),1.05])
       
        #print l[1:20]
        #print dtheta[1:20]
        popt, pcov = curve_fit(invd,l,
                dtheta,
                #bounds=(-2.0,2.0),
                #method='dogbox',
                #jac='cs'
                )

        #print popt

        #plotting the fit
        lp = np.linspace(0.5,(L+1.0)/2.0-1.0+0.5,num=200)
        ax.loglog(lp,invd(lp, *popt),
                color='royalblue',
                linestyle='--',
                lw=4,
                label=r'Power law fit, $\delta\theta\sim1/r$',
                #label=r'Exponential fit fit, $\delta\theta\sim e^{-r}$'
                basex=10,basey=10,
                zorder=2
                #+' $n$ = '
                #+str("%.4f" % n)
                )
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        #plt.legend(loc='best')
        #plt.suptitle(r"$\delta\theta$(r) vs $r$", 
        #        x=0.5, fontsize=16)
        plt.xlim([0.5,50])
        plt.ylim([10.0**(-3),1.0])
        plt.ylabel(r'$\delta\theta(r)$',fontsize=30)
        plt.xlabel(r'$r$',fontsize=30)
        plt.tick_params(which='both',width=2,labelsize=30)
        plt.tick_params(which='major',length=20)
        plt.tick_params(which='minor',length=10)
        #fig.tight_layout(pad=2.5,h_pad=2.5,w_pad=2.5)
        fig.savefig("../plot/tht_"+
                STR_L+
                STR_DELTA+
                STR_ALPHA+
                ".pdf",
                bbox_inches='tight'
                )
        
        
        plt.close('all')

for L in set(zip(*hshchar)[0]):
    subprocess.call('pdftk ../plot/tht_'+
            L+
            '* '+
            'cat output ../plot/THT_L_'+
            L+
            '.pdf',shell=True)

# removing split files
subprocess.call('rm ../plot/tht_*',shell=True)
#subprocess.call('rm ../plot/dnsty_*',shell=True)
