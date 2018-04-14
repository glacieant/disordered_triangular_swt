#!/usr/bin/env python

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
import matplotlib.figure as figure
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# setting MKL configuration to allways optimise 
os.putenv("MKL_DYNAMIC","FALSE")

## picturing the output data

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*)_DISD_([^/]*)_BOOT_([^/]*).npz')

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
sfc_nxt = [0]*HSHNUM
sfc_nxt_var = [0]*HSHNUM
n_sfc = np.zeros(HSHNUM)

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

SFC_MAX = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
XI_EST = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)

SFC_MAX_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
XI_EST_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)

XI_DELTA = np.zeros((DNUM,ANUM),dtype=np.float)
DELTA_ARR = np.zeros((DNUM,ANUM),dtype=np.float)

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

    spin_X = spin[:,0].reshape((L,L))
    spin_Y = spin[:,1].reshape((L,L))
    spin_Z = spin[:,2].reshape((L,L))
    SFC_X = np.abs(np.fft.fft2(spin_X,norm='ortho'))**2
    SFC_Y = np.abs(np.fft.fft2(spin_Y,norm='ortho'))**2
    SFC_Z = np.abs(np.fft.fft2(spin_Z,norm='ortho'))**2

    SFC = SFC_X + SFC_Y + SFC_Z
    
    # Only considering Q = (4*pi/3,0)
    KMAX_X = L/3 
    KMAX_Y = 2*L/3

    SFM = SFC[KMAX_X,KMAX_Y]/N    
    NSFM = SFC[KMAX_X,KMAX_Y+1]/N
    
    n_sfc[INDX] += 1
    
    sfc[INDX] += SFM
    sfc_var[INDX] += SFM*SFM
    sfc_nxt[INDX] += NSFM
    sfc_nxt_var[INDX] += NSFM*NSFM


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
    NSFX = sfc_nxt[i]/n_sfc[i]
    SFX_ERR = np.sqrt((sfc_var[i]/n_sfc[i]-SFX**2)/n_sfc[i])
    NSFX_ERR = np.sqrt((sfc_nxt_var[i]/n_sfc[i]-NSFX**2)/n_sfc[i])

    # scaling data

    LN = LHASH.index(STR_L)
    DN = DHASH.index(STR_DELTA)
    AN = AHASH.index(STR_ALPHA)

    Q = np.linalg.norm(b1)/L

    SFC_MAX[LN,DN,AN] = SFX
    XI_EST[LN,DN,AN] = np.sqrt(SFX/NSFX-1)/Q
    SFC_MAX_ERR[LN,DN,AN] = SFX_ERR
    XI_EST_ERR[LN,DN,AN] = (1.0/(2*Q*np.sqrt(SFX/NSFX-1)))*(
            SFX_ERR/NSFX + SFX*NSFX_ERR/(NSFX**2)
            )

    """
    # fitting a gaussian to extract a pattern along 
    # the high symmetry points

    # GAMMA point

    fig, ax = plt.subplots()

    # data for gaussian fitting
    YDATA = SFC_CUT
    XDATA = np.linspace(0.0,(2.0*np.pi),YDATA.size)
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
    ax.plot(PDATA,parabola(PDATA, *popt),'k-',
            label=
            r'Parabolic fit, $ax^2+bx+c$,'
            +r' $\sigma_p=-\frac{1}{2a}=$'
            +str("%.4f" % sigma_p)
            )

    # finding the fit
    popt, pcov = curve_fit(gaussian, XDATA, YDATA)
    sigma = np.absolute(popt[2])
    LXDATA = np.linspace(0.0,(2.0*np.pi),200)
    ax.plot(LXDATA,gaussian(LXDATA, *popt),'r-',
            label=
            r'Gaussian fit ~ $e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ ,'
            +r' $\sigma$ = '
            +str("%.4f" % sigma)
            )
    plt.legend(loc='best')

    fig.savefig("../plot/XI_STRFC"+
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

    """

for AN in range(0,ANUM):

    w,h = figure.figaspect(1.0)
    afig = plt.figure(figsize=(w,h))
    ax = afig.add_axes([0.26,0.15,0.685,0.8])
    bfig = plt.figure(figsize=(w,h))
    bx = bfig.add_axes([0.26,0.15,0.685,0.8])
    cfig = plt.figure(figsize=(w,h))
    cx = cfig.add_axes([0.26,0.15,0.685,0.8])
    dfig = plt.figure(figsize=(w,h))
    dx = dfig.add_axes([0.26,0.15,0.685,0.8])


    for DN in range(0,DNUM):

        line1 = ax.errorbar(LRAY,
                SFC_MAX[:,DN,AN][LORD],
                yerr=SFC_MAX_ERR[:,DN,AN][LORD],
                lw=2,
                marker='.',
                ms=10,
                label=r'$\delta J/J=$'+str("%.2f" % float(DHASH[DN]))
                )

        line2 = bx.errorbar(LRAY,
                XI_EST[:,DN,AN][LORD]*LRAY,
                yerr=XI_EST_ERR[:,DN,AN][LORD]*LRAY,
                lw=2,
                marker='.',
                ms=10,
                label=r'$\delta J/J=$'+str("%.2f" % float(DHASH[DN]))
                )


        popt, pcov = curve_fit(linear,LRAY[0:5],
                XI_EST[:,DN,AN][LORD][0:5])

        XI_DELTA[DN,AN] = popt[0]
        DELTA_ARR[DN,AN] = float(DHASH[DN])

    DHLIST = [float(DHASH[i]) for i in range(0,DNUM)]
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = bx.get_legend_handles_labels()

    tup = sorted(zip(DHLIST,handles1,labels1,handles2,labels2))
    DHLIST, handles1, labels1, handles2, labels2 = zip(*tup)

    legend1 = ax.legend(handles1,labels1,
            loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)

    legend2 = bx.legend(handles1,labels1,
            loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)


    #ax.set_title(r'$\alpha=$'+str("%.2f" % float(AHASH[AN])))
    #bx.set_title(r'$\alpha=$'+str("%.2f" % float(AHASH[AN])))

    ax.set_ylabel(r'$S(Q)/L$',fontsize=20)
    bx.set_ylabel(r'$\xi/L$',fontsize=20)
    ax.set_xlabel(r'$1/L$',fontsize=20)
    bx.set_xlabel(r'$1/L$',fontsize=20)

    afig.savefig("../plot/CHI_Q"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )

    bfig.savefig("../plot/XI_Q"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )

    # plotting the cumulant figure

    for LN in range(0,LNUM):

        LV = int(LHASH[LN])
        line3 = dx.errorbar(DRAY,
                XI_EST[LN,:,AN][DORD]/LV,
                yerr=XI_EST_ERR[LN,:,AN][DORD]/LV,
                #color=CLR[DN],
                lw=2,
                marker='.',
                ms=10,
                label=r'$L=$'+LHASH[LN]
                )

    LHLIST = [float(LHASH[i]) for i in range(0,LNUM)]
    handles3, labels3 = dx.get_legend_handles_labels()

    tup = sorted(zip(LHLIST,handles3,labels3))
    LHLIST, handles3, labels3 = zip(*tup)

    legend3 = dx.legend(handles3,labels3,
            loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)

    dx.set_ylabel(r'$\xi/L$',fontsize=20)
    dx.set_xlabel(r'$\delta J/J$',fontsize=20)

    dfig.savefig("../plot/XI_BINDER"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )


    # fitting the corr. length vs Delta exponential

    YDATA = np.log(XI_DELTA[:,AN])
    XDATA = DELTA_ARR[:,AN]
    cx.plot(XDATA,YDATA,'b.',label='Data')

    popt, pcov = curve_fit(invpar,XDATA,YDATA)

    PDATA = np.linspace(XDATA.min(),
            XDATA.max(),
            200)
    cx.plot(PDATA,invpar(PDATA, *popt),'k-',
            label=
            r'$\log\xi = a + b/(\delta J/J)^2$,'
            +'\n'
            +r' $a=$'
            +str("%.4f" % popt[0])
            +r' $b=$'
            +str("%.4f" % popt[1])
            )

    legend3 = cx.legend(loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)

    cx.set_title(r'$\alpha=$'+str("%.2f" % float(AHASH[AN])))
    cx.set_ylabel(r'$\log\xi$',fontsize=20)
    cx.set_xlabel(r'$\delta J/J$',fontsize=20)

    cfig.savefig("../plot/XI_DELTA"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )

    plt.close('all')

