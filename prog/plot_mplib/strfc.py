#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
os.chdir("../sharp_wall/out/data")
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

# defining a tolerance limit for calculations
zero = 10.0**(-8)

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

        # setting up the clipping of the density plot
        MAXCLIP = 1.0
        MINCLIP = 10**(-10)
        # safely bottoming out for log scale
        sfc[i] += MINCLIP/(10.0**5)
 
        #cfx = np.concatenate((np.flipud(sfc[i]),sfc[i]),axis=0)
        #cfx = np.concatenate((bfx,bfx),axis=0)
       
        # shearing the date into the BZ
        btx = 1000
        dtx = btx/2
        bfx = np.zeros((2*btx,2*btx),dtype=np.float64)
        dfx = np.zeros((2*dtx,2*dtx),dtype=np.float64)
        
        for m in range(-btx,btx):
            for n in range(-btx,btx):
                G1 = (1.0*m/btx)
                G2 = (1.0*n/btx)
                p = int(2*L*G1)
                q = int(2*L*(G1+G2*3.0**0.5)/2)
                if p >= 0 and p < L and q >= 0 and q < L:
                    bfx[n+btx,m+btx] = sfc[i][p,q]
                
                if p >= L and p < 2*L and q >= L and q < 2*L:
                    bfx[n+btx,m+btx] = sfc[i][p-L,q-L]
                if p >= 0 and p < L and q >= L and q < 2*L:
                    bfx[n+btx,m+btx] = sfc[i][p,q-L]
                if p >= L and p < 2*L and q >= 0 and q < L:
                    bfx[n+btx,m+btx] = sfc[i][p-L,q]
                
        #post processing to fill out the brillouin zone
        for x in range(-dtx,dtx):
            for y in range(-dtx,dtx):
                if x < 0:
                    if (1.0*y/x) <= 1.0/(3.0**0.5): 
                        # rotating -120 degrees
                        xp = int((1.0/2.0)*(-x+y*3.0**0.5))
                        yp = int((1.0/2.0)*(-x*3.0**0.5-y))
                        if xp >= -btx and xp < btx and yp >= -btx and yp < btx:
                            dfx[y+dtx,x+dtx] = bfx[yp+btx,xp+btx]
                    elif (1.0*y/x) > 1.0/(3.0**0.5): 
                        # rotating +120 degrees
                        xp = int((1.0/2.0)*(-x-y*3.0**0.5))
                        yp = int((1.0/2.0)*(x*3.0**0.5-y))
                        if xp >= -btx and xp < btx and yp >= -btx and yp < btx:
                            dfx[y+dtx,x+dtx] = bfx[yp+btx,xp+btx]
                elif x > 0:
                    if (1.0*y/x) <= -1.0/(3.0**0.5): 
                        # rotating +120 degrees
                        xp = int((1.0/2.0)*(-x-y*3.0**0.5))
                        yp = int((1.0/2.0)*(x*3.0**0.5-y))
                        if xp >= -btx and xp < btx and yp >= -btx and yp < btx:
                            dfx[y+dtx,x+dtx] = bfx[yp+btx,xp+btx]
                    else:
                        # copying original
                        dfx[y+dtx,x+dtx] = bfx[y+btx,x+btx]
                
                elif x == 0:
                    if y < 0:
                        # rotating +120 degrees
                        xp = int((1.0/2.0)*(-x-y*3.0**0.5))
                        yp = int((1.0/2.0)*(x*3.0**0.5-y))
                        if xp >= -btx and xp < btx and yp >= -btx and yp < btx:
                            dfx[y+dtx,x+dtx] = bfx[yp+btx,xp+btx]
                    else:
                        # copying original
                        dfx[y+dtx,x+dtx] = bfx[y+btx,x+btx]
        
        cax = ax.imshow(dfx,
                origin='lower',
                #norm=Normalize(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                norm=LogNorm(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                interpolation='nearest',
                #extent=(0.0,2.0*np.pi,0.0,2.0*np.pi),
                cmap=cm.viridis,aspect='auto')
        cbar = fig.colorbar(cax,shrink=0.5,format='%.0e')
        
        plt.suptitle(r"Static Structure Factor, $\chi (q)/L^{2}$", 
                x=0.5, fontsize=16)
        plt.title(r"$\Delta$ = "+STR_DELTA+" , "
                +r"$\alpha$ = "+STR_ALPHA+" , "
                +"|E|/$L^{2}$ = "+str("%f" % en[i]),
                x=0.6,fontsize=12)
        xticks = [0,dtx,2*dtx]
        xtickslabel = [r'$-2\pi$',0,r'$2\pi$']
        yticks = [0,dtx,2*dtx]
        ytickslabel = [r'$-2\pi$',0,r'$2\pi$']
        plt.xticks(xticks,xtickslabel)
        plt.yticks(yticks,ytickslabel)
        plt.text(dtx+dtx/20,dtx+dtx/25,r'$\Gamma$',color='black',fontsize=16)
        plt.text(2*dtx/3+dtx+dtx/25,dtx+dtx/25,r'$K$',color='black',fontsize=16)
        plt.text(
                int((np.pi/3.0**0.5)*dtx/(2.0*np.pi)+dtx+dtx/25),
                int((np.pi/3.0**0.5)*dtx/(2.0*np.pi)+dtx++dtx/25),
                r'$M$',color='black',fontsize=16)
        smax = np.amax(sfc[i])
        #smax = np.amax(dfx)
        ax.set_aspect('equal')
        plt.figtext(.8,.85,
                r'$\chi^{\mathrm{max}}/L^{2}$ = '+str("%.2e" % smax),
                fontsize=12)
        plt.xlabel('L = '+STR_L,fontsize=16)
        fig.tight_layout(pad=2.5,h_pad=2.5,w_pad=2.5)
        fig.savefig("../plot/stsc_"+
                STR_L+
                STR_DELTA+
                STR_ALPHA+
                ".pdf"
                )

        plt.close('all')

for L in zip(*hshchar)[0]:
    subprocess.call('pdftk ../plot/stsc_'+
            L+
            '* '+
            'cat output ../plot/STRFC_L_'+
            L+
            '.pdf',shell=True)

# removing split files
subprocess.call('rm ../plot/stsc_*',shell=True)
