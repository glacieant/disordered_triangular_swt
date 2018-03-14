#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
#folder = "sharp_wall"
#folder = "single_impurity"
folder = "zero_field_classical"
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
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

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

# reciprocal lattice vectors
b1 = np.array([2.0*np.pi,-2.0*np.pi/3.0**0.5])
b2 = np.array([0.0,4.0*np.pi/3.0**0.5])

# Defining model function to be used to fit 1d data

def gaussian(x,a,b,c):
    
    return a*np.exp(-(x-b)**2/(2.0*c**2))

def parabola(x,a,b,c):
    
    return a*(x**2)+b*x+c


def skew_gauss(x,alpha,a,b,c):

    return a*skewnorm.pdf(x,alpha,b,np.sqrt(2.0)*c)


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

    spin_X = spin[:,0].reshape((L,L))
    spin_Y = spin[:,1].reshape((L,L))
    spin_Z = spin[:,2].reshape((L,L))
    SFC_X = np.abs(np.fft.fft2(spin_X,norm='ortho'))**2
    SFC_Y = np.abs(np.fft.fft2(spin_Y,norm='ortho'))**2
    SFC_Z = np.abs(np.fft.fft2(spin_Z,norm='ortho'))**2

    SFC = SFC_X + SFC_Y + SFC_Z

    sfc[INDX] += SFC

for i in range(0,HSHNUM):

    if n_sfc[i] > 0:

        fig, ax = plt.subplots()

        STR_L = hshchar[i][0]
        STR_DELTA = hshchar[i][1]
        STR_ALPHA = hshchar[i][2]
        L = int(STR_L)
        N = L**2

        sfc[i] *= (1.0/(n_sfc[i]*N))
        en[i] *= (1.0/n_sfc[i])

        # processing the BZ cut
        HCUT = 1
        VCUT = 2
        XMAXMO = 1.0
        YMAXMO = 1.0

        if VCUT == 0 and HCUT != 0:
            RMAX = int((XMAXMO*L)/HCUT)
        elif VCUT != 0 and HCUT == 0:
            RMAX = int((YMAXMO*L)/VCUT)
        else:
            NX = int((XMAXMO*L)/HCUT)
            NY = int((YMAXMO*L)/VCUT)
            if NX >= NY:
                RMAX = NY
            else:
                RMAX = NX
            
        SFC_CUT = np.zeros(RMAX)

        for k in range(0,RMAX):
            SFC_CUT[k] = sfc[i][k*VCUT,k*HCUT]

        #fitting a gaussian to extract a pattern along the high symmetry points

        # GAMMA point

        fig, ax = plt.subplots()

        # data for gaussian fitting
        YDATA = SFC_CUT
        XDATA = np.linspace(0.0,XMAXMO*(2.0*np.pi),YDATA.size)
        ax.plot(XDATA,YDATA,'b.',label='Data')

        # finding the parabolic fit

        MAX_N = np.argmax(YDATA)
        Y_VAL = YDATA[MAX_N-1:MAX_N+2]
        X_VAL = XDATA[MAX_N-1:MAX_N+2]

        popt, pcov = curve_fit(parabola,X_VAL,Y_VAL)

        sigma_p = -1/(2*popt[0])

        PDATA = np.linspace(X_VAL.min(),
                X_VAL.max(),
                200)
        #sigma = np.absolute(popt[3])
        ax.plot(PDATA,parabola(PDATA, *popt),'k-',
        #ax.plot(XDATA,skew_gauss(XDATA, *popt),'r-',
                label=
                #r'Gaussian fit ~ $e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ ,'
                #label=r'Skew-normal fit around'+'\n'+r'($\frac{x-\mu}{\sqrt{2}\sigma}$) with,'
                #+r' $\sigma$ = '
                #+str("%.4f" % sigma)
                #+"\n"
                r'Parabolic fit, $ax^2+bx+c$,'
                +r' $\sigma_p=-\frac{1}{2a}=$'
                +str("%.4f" % sigma_p)
                )

        # finding the fit
        popt, pcov = curve_fit(gaussian, XDATA, YDATA)
        #popt, pcov = curve_fit(skew_gauss, XDATA, YDATA,p0=[-1.25,np.amax(YDATA),-1.25,0.05])
        sigma = np.absolute(popt[2])
        LXDATA = np.linspace(0.0,XMAXMO*(2.0*np.pi),200)
        #sigma = np.absolute(popt[3])
        ax.plot(LXDATA,gaussian(LXDATA, *popt),'r-',
        #ax.plot(XDATA,skew_gauss(XDATA, *popt),'r-',
                label=
                r'Gaussian fit ~ $e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ ,'
                #label=r'Skew-normal fit around'+'\n'+r'($\frac{x-\mu}{\sqrt{2}\sigma}$) with,'
                +r' $\sigma$ = '
                +str("%.4f" % sigma)
                #+"\n"
                #r'Parabolic fit, $ax^2+bx+c$,'
                #+r' $\sigma_p=-\frac{1}{2a}=$'
                #+str("%.4f" % sigma_p)
                )
        #plt.ylim([0,SFC_MAX])
        #plt.ylabel(r'$S_{\Gamma}/L^2$', fontsize=18)
        #plt.xlabel(r'$\omega$/J', fontsize=18)
        plt.legend(loc='best')
        #plt.title(r'$S_{\Gamma}(\omega)$ at $\Delta$='
        #        +str("%.3f" % DELTA[d])+r', $\alpha$='+str("%.3f" % ALPHA[a]), 
        #        fontsize=20,y=1.09)
 
        fig.savefig("../plot/XI_STRFC"+
                "_L_"+
                STR_L+
                "_DLT_"+
                STR_DELTA+
                "_ALP_"+
                STR_ALPHA+
                "_DNM_"+
                ".pdf",
                bbox_inches='tight',
                transparent=True
                )
        plt.close('all')

