#!/usr/bin/env python

import os
import datetime
import multiprocessing as mp
#folder = "sharp_wall"
#folder = "single_impurity"
folder = "xy_gaussian"
#folder = "test_chir_2"
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
sfc = [0]*HSHNUM
sfc_var = [0]*HSHNUM
sfc_nxt = [0]*HSHNUM
sfc_nxt_var = [0]*HSHNUM
cfc = [0]*HSHNUM
cfc_var = [0]*HSHNUM
csfc = [0]*HSHNUM
csfc_var = [0]*HSHNUM
csfc_nxt = [0]*HSHNUM
csfc_nxt_var = [0]*HSHNUM
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
CFC = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
CSFC = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
CXI_EST = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
ITNUM = np.zeros((LNUM,DNUM,ANUM),dtype=np.int)

SFC_MAX_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
XI_EST_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
CFC_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
CSFC_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)
CXI_EST_ERR = np.zeros((LNUM,DNUM,ANUM),dtype=np.float)

XI_DELTA = np.zeros((DNUM,ANUM),dtype=np.float)
CXI_DELTA = np.zeros((DNUM,ANUM),dtype=np.float)
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
    SFC_X = np.abs(np.fft.fft2(spin_X,norm='ortho'))**2
    SFC_Y = np.abs(np.fft.fft2(spin_Y,norm='ortho'))**2

    SFC = (SFC_X + SFC_Y)
    
    # Only considering Q = (4*pi/3,0)
    KMAX_X = L/3 
    KMAX_Y = 2*L/3

    SFM = SFC[KMAX_X,KMAX_Y]/N    
    NSFM = (SFC[KMAX_X+1,KMAX_Y]+SFC[KMAX_X-1,KMAX_Y]
            +SFC[KMAX_X,KMAX_Y+1]+SFC[KMAX_X,KMAX_Y-1]
            )/N/4

    n_sfc[INDX] += 1
    
    sfc[INDX] += SFM
    sfc_var[INDX] += SFM*SFM
    sfc_nxt[INDX] += NSFM
    sfc_nxt_var[INDX] += NSFM*NSFM

    LTR = L-1
    NTR = LTR**2

    elt = np.zeros((NTR,3),dtype=np.int)
    lmap.eltriangle(L,elt)

    tsr_spin = spin[elt]

    epsilon = np.zeros((2,2),dtype=np.float)
    epsilon[0,1] =  1.0
    epsilon[1,0] = -1.0

    alt = np.zeros((3,3),dtype=np.float)
    alt[0,1] = alt[1,2] = alt[2,0] = 1.0
    
    chir = 2.0*np.einsum('iab,icd,ac,bd->i',
            tsr_spin,tsr_spin,alt,epsilon)/(3.0*np.sqrt(3.0))

    chir = chir.reshape((LTR,LTR))
    fchir = np.abs(np.fft.fft2(chir,norm='ortho'))**2

    fchir = fchir

    KMAX_X = 0
    KMAX_Y = 0

    CFM = fchir[KMAX_X,KMAX_Y]/NTR    
    NCFM = (fchir[KMAX_X+1,KMAX_Y]+fchir[KMAX_X-1,KMAX_Y]
            +fchir[KMAX_X,KMAX_Y+1]+fchir[KMAX_X,KMAX_Y-1]
            )/NTR/4

    CVL = np.linalg.norm(np.sum(chir,axis=0))/(3*np.sqrt(3)*NTR/2)

    cfc[INDX] += CVL
    cfc_var[INDX] += CVL*CVL

    csfc[INDX] += CFM
    csfc_var[INDX] += CFM*CFM
    csfc_nxt[INDX] += NCFM
    csfc_nxt_var[INDX] += NCFM*NCFM


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

    CFX = cfc[i]/n_sfc[i]
    CFX_ERR = np.sqrt((cfc_var[i]/n_sfc[i]-CFX**2)/n_sfc[i])
 
    CSFX = csfc[i]/n_sfc[i]
    CSFX_ERR = np.sqrt((csfc_var[i]/n_sfc[i]-CSFX**2)/n_sfc[i])
    NCSFX = csfc_nxt[i]/n_sfc[i]
    NCSFX_ERR = np.sqrt((csfc_nxt_var[i]/n_sfc[i]-NCSFX**2)/n_sfc[i])
       
    # scaling data

    LN = LHASH.index(STR_L)
    DN = DHASH.index(STR_DELTA)
    AN = AHASH.index(STR_ALPHA)

    #Q = np.linalg.norm(b1)/L
    #QC = np.linalg.norm(b1)/LTR

    #eric's formula

    Q = np.sqrt(16*(np.sin(np.pi/L)**2)/3)
    QC = np.sqrt(16*(np.sin(np.pi/LTR)**2)/3)

    ITNUM[LN,DN,AN] = n_sfc[i]
    SFC_MAX[LN,DN,AN] = SFX
    CFC[LN,DN,AN] = CFX
    CSFC[LN,DN,AN] = CSFX

    XI_EST[LN,DN,AN] = np.sqrt(np.abs(SFX/NSFX-1))/Q
    CXI_EST[LN,DN,AN] = np.sqrt(np.abs(CSFX/NCSFX-1))/QC
    
    SFC_MAX_ERR[LN,DN,AN] = SFX_ERR
    CFC_ERR[LN,DN,AN] = CFX_ERR
    CSFC_ERR[LN,DN,AN] = CSFX_ERR
    XI_EST_ERR[LN,DN,AN] = (1.0/(2*Q*np.sqrt(np.abs(SFX/NSFX-1))))*(
            SFX_ERR/NSFX + SFX*NSFX_ERR/(NSFX**2)
            )
    CXI_EST_ERR[LN,DN,AN] = (1.0/(2*QC*np.sqrt(np.abs(CSFX/NCSFX-1))))*(
            CSFX_ERR/NCSFX + CSFX*NCSFX_ERR/(NCSFX**2)
            )

CLR = [tu_rot,tu_blau,tu_gruen]

for AN in range(0,ANUM):

    w,h = figure.figaspect(1.0)
    afig = plt.figure(figsize=(w,h))
    ax = afig.add_axes([0.26,0.15,0.685,0.8])
    a2fig = plt.figure(figsize=(w,h))
    a2x = a2fig.add_axes([0.26,0.15,0.685,0.8])
    a3fig = plt.figure(figsize=(w,h))
    a3x = a3fig.add_axes([0.26,0.15,0.685,0.8])
    
    for DN in range(0,DNUM):

        line1 = ax.errorbar(LRAY,
                SFC_MAX[:,DN,AN][LORD],
                yerr=SFC_MAX_ERR[:,DN,AN][LORD],
                lw=2,
                marker='.',
                ms=10,
                label=r'$\Delta=$'+DHASH[DN],
                )

        line3 = a2x.errorbar(LRAY,
                CFC[:,DN,AN][LORD],
                yerr=CFC_ERR[:,DN,AN][LORD],
                lw=2,
                marker='.',
                ms=10,
                label=r'$\delta J/J=$'+str("%.2f" % float(DHASH[DN]))
                )

        
        line4 = a3x.errorbar(LRAY,
                CSFC[:,DN,AN][LORD],
                yerr=CSFC_ERR[:,DN,AN][LORD],
                lw=2,
                marker='.',
                ms=10,
                label=r'$\delta J/J=$'+str("%.2f" % float(DHASH[DN]))
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
    legend3 = a2x.legend(handles1,labels1,
            loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)

    legend4 = a3x.legend(handles1,labels1,
            loc='best',
            fontsize=12,
            markerscale=1,
            facecolor='w',
            edgecolor='k',
            framealpha=1)

    time = str(datetime.datetime.now()) 
    itmin = str(np.amin(ITNUM[:,:,AN])) 
    
    ax.set_title('DATE = '
            +time
            +' MINIMUM REALIZATION = '
            +itmin
            )

    #ax.set_xlim([0.0,0.06])
    #a2x.set_xlim([0.0,0.06])
    #a3x.set_xlim([0.0,0.06])
    #bx.set_xlim([0.0,0.06])
    #b2x.set_xlim([0.0,0.06])
   
    ax.set_xlim(left=0.0)
    a2x.set_xlim(left=0.0)
    a3x.set_xlim(left=0.0)
  
    #ax.set_ylim(bottom=0.45,top=0.5)
    a2x.set_ylim(bottom=0.0)
    a3x.set_ylim(bottom=0.0)

    ax.set_ylabel(r'$S(Q)$',fontsize=30)
    a2x.set_ylabel(r'$|\mathbf{\chi}|$',fontsize=20)
    a3x.set_ylabel(r'$S_\chi(0)/L$',fontsize=20)
    ax.set_xlabel(r'$1/L$',fontsize=30)
    a2x.set_xlabel(r'$1/L$',fontsize=20)
    a3x.set_xlabel(r'$1/L$',fontsize=20)
        
    ax.tick_params(which='both',width=2,
            labelsize=30,direction='in',
            bottom=True,top=True,
            left=True,right=True)
    ax.tick_params(which='major',length=20)
    ax.tick_params(which='minor',length=10)
    
    a2x.tick_params(which='both',width=2,
            labelsize=30,direction='in',
            bottom=True,top=True,
            left=True,right=True)
    a2x.tick_params(which='major',length=20)
    a2x.tick_params(which='minor',length=10)

    a3x.tick_params(which='both',width=2,
            labelsize=30,direction='in',
            bottom=True,top=True,
            left=True,right=True)
    a3x.tick_params(which='major',length=20)
    a3x.tick_params(which='minor',length=10)


    afig.savefig("../plot/S_Q-XY"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )
    a2fig.savefig("../plot/CHIR-XY"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )

    a3fig.savefig("../plot/CS_Q-XY"
            +"_ALPHA_"
            +str("%.4f" % float(AHASH[AN]))
            +".pdf",
            bbox_inches='tight'
            )
    
    plt.close('all')

