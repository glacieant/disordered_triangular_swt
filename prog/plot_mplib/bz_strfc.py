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

## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# the font styleset
# from matplotlib import rcParams
# rcParams['font.serif'] = ['Times New Roman']
# rcParams['font.family'] = 'serif'

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
pm_sfc = [0]*HSHNUM
mp_sfc = [0]*HSHNUM
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
    sfc[INDX] += SFC_X+SFC_Y+SFC_Z
    """
    PM_SFC_X = np.fft.fft(
            np.fft.ifft(spin[:,0].reshape(L,L),norm='ortho')
            ,axis=0,norm='ortho').real
    PM_SFC_Y = np.fft.fft(
            np.fft.ifft(spin[:,1].reshape(L,L),norm='ortho')
            ,axis=0,norm='ortho').real
    PM_SFC_Z = np.fft.fft(
            np.fft.ifft(spin[:,2].reshape(L,L),norm='ortho')
            ,axis=0,norm='ortho').real
    pm_sfc[INDX] += PM_SFC_X+PM_SFC_Y+PM_SFC_Z
    MP_SFC_X = np.fft.ifft(
            np.fft.fft(spin[:,0].reshape(L,L),norm='ortho')
            ,axis=0,norm='ortho').real
    MP_SFC_Y = np.fft.ifft(
            np.fft.fft(spin[:,1].reshape(L,L),norm='ortho')
            ,axis=0,norm='ortho').real
    MP_SFC_Z = np.fft.ifft(
            np.fft.fft(spin[:,2].reshape(L,L),norm='ortho')
            ,axis=0,norm='ortho').real
    pm_sfc[INDX] += MP_SFC_X+MP_SFC_Y+MP_SFC_Z
    """

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
        en[i] *= (1.0/n_sfc[i])

        # setting up the clipping of the density plot
        MAXCLIP = np.amax(sfc[i])
        #MAXCLIP = 10.0**(-2)
        #MINCLIP = 10.0**(-4)
        MINCLIP = np.amin(sfc[i])
        # safely bottoming out for log scale
        #sfc[i] += MINCLIP/(10.0**5)
        
        # shearing the date into the BZ
        LX = 200
        sfx = np.zeros((LX,LX),dtype=np.float64)
        
        for m in range(0,LX):
            for n in range(0,LX):
                G1 = (1.0*m/LX)
                G2 = (1.0*n/LX)
                p = int(L*G1)
                q = int(L*(G1+G2*3.0**0.5)/2)
                if p >= 0 and p < L:
                    if q >= 0 and q < L:
                        sfx[n,m] = sfc[i][p,q]

       
        """
        #plt.plot(sfc[i][0])
        #plt.show()
        #raw_input()
        # creating the momentum space mesh
        DX = 100
        sfc_q = np.zeros((DX,DX),dtype=np.float64)
        for X in range(0,DX):
            for Y in range(0,DX):
                qx = 4*np.pi*X/DX-2*np.pi
                qy = 4*np.pi*Y/DX-2*np.pi
                if qx < 0:
                    qx = 2*np.pi + qx
                if qy < 0:
                    qy = 2*np.pi + qy
                m1 = int(qx*L/(2.0*np.pi))
                m2 = int((qx+qy*3.0**0.5)*L/(2.0*np.pi))
                if abs(m1) < L and abs(m2) < L:
                    if m1 > 0 and m2 > 0:
                        sfc_q[X][Y] = sfc[i][m1][m2]
                    elif m1 > 0 and m2 < 0:
                        sfc_q[X][Y] = pm_sfc[i][m1][abs(m2)]
                    elif m1 < 0 and m2 > 0:
                        sfc_q[X][Y] = pm_sfc[i][abs(m1)][m2]
                    elif m1 < 0 and m2 < 0:
                        sfc_q[X][Y] = pm_sfc[i][abs(m1)][abs(m2)]

        """
        cax = ax.imshow(sfc[i],
                norm=Normalize(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                #norm=LogNorm(vmin=MINCLIP,vmax=MAXCLIP,clip=False),
                interpolation='nearest',
                #extent=(0.0,1.0,0.0,1.0),
                cmap=cm.jet,aspect='auto')

        print sfc[i][8,16]
        cbar = fig.colorbar(cax,shrink=0.5,format='%.0e')

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
        fig.savefig("../plot/strfc"+
                "_L_"+
                STR_L+
                "_DLT_"+
                STR_DELTA+
                "_ALP_"+
                STR_ALPHA+
                "_DNM"+
                ".pdf",
                bbox_inches='tight',
                transparent=True
                )
        plt.close('all')


for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite '+
                    '-dPDSETTINGS=/prepress -sOutputFile='+
                '../plot/STRFC_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/strfc'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM"+
                '* ',shell=True)
            subprocess.call('rm ../plot/strfc'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM"+
                '*',shell=True)

