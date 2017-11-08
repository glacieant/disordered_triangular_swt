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

    #T = FNDATA['T'][()]
    #ensys = FNDATA['ensys']
    #lmult = FNDATA['lmult']
    #vsys = FNDATA['vsys']
    spin = FNDATA['spin']

    N = L**2

    #en[INDX] += 0.5*(np.sum(fermi(ensys,T)*ensys)
    #        -np.sum(lmult))/N
   
    # reshaping the eigenvectors 
    #farray = fermi(ensys,T)
    #hsys = vsys[:,np.where(farray>zero)]
    #[NVEC,d,NEN] = np.array(hsys.shape)
    #if NVEC != 2*N:
    #    print "eigenvector size wrong!"
    #    quit()
    #vsys = np.reshape(hsys,(2*N,NEN))
    #U = np.zeros((2,NEN,L,L),dtype=np.complex64)
    #U[0] = np.reshape(np.einsum('ij->ji',vsys[:N,:]),(NEN,L,L))
    #U[1] = np.reshape(np.einsum('ij->ji',vsys[N:,:]),(NEN,L,L))
    #UDAGU = np.einsum('amij,bnij->mnabij',U.conj(),U)
    #GAMMA = np.fft.fft2(UDAGU.real,norm='ortho')
    #GAMMA += 1.0j*np.fft.fft2(UDAGU.imag,norm='ortho')
    
    #SFC = 2*np.einsum("mmabij,nnabij->ij",GAMMA,GAMMA.conj())
    #SFC = -2*np.einsum("mnabij,mnabij->ij",GAMMA,GAMMA.conj())
    #SFC -= np.einsum("mmaaij,nnbbij->ij",GAMMA,GAMMA.conj())
    #SFC += np.einsum("mnaaij,mnbbij->ij",GAMMA,GAMMA.conj())
    SFC_X = np.abs(np.fft.fft2(spin[:,0].reshape(L,L),norm='ortho'))**2
    SFC_Y = np.abs(np.fft.fft2(spin[:,1].reshape(L,L),norm='ortho'))**2
    SFC_Z = np.abs(np.fft.fft2(spin[:,2].reshape(L,L),norm='ortho'))**2
    sfc[INDX] += SFC_X + SFC_Y + SFC_Z


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

        sfc[i] *= (1.0/(n_sfc[i]*N))
        #sfc[i] *= (1.0/n_sfc[i])
        en[i] *= (1.0/n_sfc[i])

        # setting up the clipping of the density plot
        MAXCLIP = np.amax(sfc[i])
        #MAXCLIP = 10.0**(-2)
        MINCLIP = MAXCLIP/10.0**(5)
        #MINCLIP = np.amin(sfc[i])
        # safely bottoming out for log scale
        sfc[i] += MINCLIP/(10.0**5)

        #plt.plot(sfc[i][0])
        #plt.show()
        #raw_input()

        cax = ax.imshow(sfc[i],
                #norm=Normalize(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                norm=LogNorm(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                interpolation='nearest',
                extent=(0.0,1.0,0.0,1.0),
                cmap=cm.jet,aspect='auto')
        cbar = fig.colorbar(cax,shrink=0.5,format='%.0e')
        """        
        plt.suptitle(r"Static Structure Factor, $\chi (q)$", 
                x=0.5, fontsize=16)
        plt.title('Two bond impurities in The Middle',x=0.6,fontsize=12)
        #plt.title(r'Anomalous dipole density along $\hat{e}_3=(-1/2,\sqrt{3}/2)$',x=0.6,fontsize=12)


        """
        plt.suptitle(r"Static Structure Factor, $\chi (q)/L^{2}$", 
                x=0.5, fontsize=16)
        plt.title(r"$\Delta$ = "+STR_DELTA+" , "
                +r"$\alpha$ = "+STR_ALPHA+" , "
                +"|E|/$L^{2}$ = "+str("%f" % en[i]),
                x=0.6,fontsize=12)
                
        smax = np.amax(sfc[i])
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
