#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
folder = "realspace_mf"
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
#from matplotlib import rcParams
#rcParams['font.serif'] = ['Times New Roman']
#rcParams['font.family'] = 'serif'

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

# importing lattice geometry system
import lattice_map as lmap
ZCO = 12

X = np.zeros(1)
Y = np.zeros(1)
Z = np.zeros(1)

a = 0.5

# three diffrent translation vector

avec = np.array([[-0.5,((3.0)**0.5)/2],
    [0.5,((3.0)**0.5)/2],
    [1.0,0.0],
    [-1.5,((3.0)**0.5)/2],
    [0.0,(3.0**0.5)],
    [1.5,((3.0)**0.5)/2]])*a

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

    TOL = FNDATA['err'][()]

    if TOL < 101:
        
        IDISD = int(match.group(4))
        BTNUM = int(match.group(5))
       
        T = FNDATA['T'][()]
        J = FNDATA['J']
        inum = FNDATA['inum']
        err=FNDATA['err']
        ensys = FNDATA['ensys']
        lmult = FNDATA['lmult']
        bond = np.absolute(FNDATA['bond'])
        spin = FNDATA['spin']

        N = L**2

        nbr = np.zeros((N,ZCO),dtype=np.int)
        lmap.lattice_map(L,nbr)

        # laying out the lattice skeleton
        if len(X) != N:

            X = np.zeros((L,L))
            Y = np.zeros((L,L))
            for i in range(0,L):
                for j in range(0,L):
                    X[i,j]=(i+j*(1.0/2.0))*a
                    Y[i,j]=j*((3.0)**0.5/2.0)*a

        # plotting the configuration

        fig, ax = plt.subplots()

        # defining the plane for the spins

        S_SQR = np.einsum('ij->i',spin**2)**(0.5)
        SMAX = np.amax(S_SQR)

        e1 = spin[0]
        if np.linalg.norm(e1) > 10.0**(-5):
            e1 = e1/np.linalg.norm(e1)
            e2 = np.cross(e1,np.cross(e1,spin[1]))
            e2 = e2/np.linalg.norm(e2)
        else:
            e1 = np.array([1.0,0.0,0.0])
            e2 = np.array([0.0,1.0,0.0])

        U = np.einsum('ij,j->i',spin,e1).reshape((L,L)).transpose()
        V = np.einsum('ij,j->i',spin,e2).reshape((L,L)).transpose()

        VMIN = 0.0
        VMAX = 2.0
        VGRID = 20
        
        BMIN = np.amin(bond)
        BMAX = np.amax(bond)
        JMIN = np.amin(J)
        JMAX = np.amax(J)

        BWDTH = (BMAX-BMIN)/VGRID
        JWDTH = (JMAX-JMIN)/VGRID

        JLW = np.zeros((N,N))
        BLW = np.zeros((N,N))

        if (BMAX > 10.0**(-4)):

            for i in range(0,N):

                for j in nbr[i]:

                    JLW[i,j] = (VMIN + 
                            (np.floor((J[i,j]-JMIN)/JWDTH)
                                /VGRID)*(VMAX-VMIN)
                            )
                    BLW[i,j] = (VMIN + 
                            (np.floor((bond[i,j]-BMIN)/BWDTH)
                                /VGRID)*(VMAX-VMIN)
                            )
                

            for i in range(0,N):

                for p in range(1,4):

                    j = nbr[i,p]
                    q = p -1

                    cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
                    cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]

                    # plotting the lattice
                    line = plt.Line2D(cor_x,cor_y,
                            color='gray',
                            alpha=0.25,
                            ls='solid',
                            lw=0.3)
                    ax.add_line(line)

                    # plotting the couplings
                    line = plt.Line2D(cor_x,cor_y,
                            color='red',
                            alpha=0.5,
                            ls='solid',
                            lw=JLW[i,j])
                    ax.add_line(line)
                    
                    # plotting the bonds
                    line = plt.Line2D(cor_x,cor_y,
                            color='blue',
                            alpha=0.5,
                            ls='dashed',
                            dashes=(1.0,1.0),
                            lw=BLW[i,j])
                    ax.add_line(line)
                    
                
                # plotting additonal bonds
                for p in range(0,3):

                    j = nbr[i,6+p]
                    q = p +3

                    cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
                    cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]

                    # plotting the couplings
                    line = plt.Line2D(cor_x,cor_y,
                            color='red',
                            alpha=0.5,
                            ls='solid',
                            lw=JLW[i,j])
                    ax.add_line(line)

                    line = plt.Line2D(cor_x,cor_y,
                            color='blue',
                            alpha=0.5,
                            ls='dashed',
                            dashes=(1.0,1.0),
                            lw=BLW[i,j])
                    ax.add_line(line)

        else:


            for i in range(0,N):

                for j in nbr[i]:

                    JLW[i,j] = (VMIN + 
                            (np.floor((J[i,j]-JMIN)/JWDTH)
                                /VGRID)*(VMAX-VMIN)
                            )

            for i in range(0,N):

                for p in range(1,4):

                    j = nbr[i,p]
                    q = p -1

                    cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
                    cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]

                    # plotting the lattice
                    line = plt.Line2D(cor_x,cor_y,
                            color='gray',
                            alpha=0.1,
                            ls='solid',
                            lw=0.25)
                    ax.add_line(line)

                    # plotting the couplings
                    line = plt.Line2D(cor_x,cor_y,
                            color='red',
                            alpha=0.5,
                            ls='solid',
                            lw=JLW[i,j])
                    ax.add_line(line)

                # plotting additonal bonds
                for p in range(0,3):

                    j = nbr[i,6+p]
                    q = p +3

                    cor_x = [X[i%L,i/L],X[i%L,i/L]+avec[q,0]]
                    cor_y = [Y[i%L,i/L],Y[i%L,i/L]+avec[q,1]]

                    # plotting the couplings
                    line = plt.Line2D(cor_x,cor_y,
                            color='red',
                            alpha=0.5,
                            ls='solid',
                            lw=JLW[i,j])
                    ax.add_line(line)


        # labeling scale of coupling and bonds
        ax.text(-1.25*a,(3.0*L/4.5)*a,r'$b_{\mathrm{max}}$ = '+str("%.4f" % BMAX)
                +'\n'+r'$J_{\mathrm{max}}$ = '+str("%.4f" % JMAX)
                +'\n'+r'$S_{\mathrm{max}}$ = '+str("%.4f" % SMAX),
                fontsize=12)

        # picturing the spin orientation
        Q=ax.quiver(X,Y,U,V,
                color='black',pivot='mid',
                angles='xy',
                scale=1,
                scale_units='xy',
                )
        plt.axis([X.min()-2.0*a,X.max()+2.0*a,
            Y.min()-2.0*a,Y.max()+2.0*a])

        # adding label patch
        patch = [mpatches.Patch(color='red',alpha=0.5),
                mpatches.Patch(color='blue',alpha=1.0)]
        label = [r'Coupling, $J$',r'Bond parameter, $b$']

        legend=ax.legend(patch,label,
                loc='lower right',shadow=True)

        EN = 0.5*(np.sum(fermi(ensys,T)*ensys)
                -np.sum(lmult))/N

        #plt.suptitle("Spin & Bond Configuration", x=0.5, y=0.99, fontsize=16)
        plt.title(r"$\Delta$ = "+STR_DELTA+" , "
                +r"$\alpha$ = "+STR_ALPHA
                #+" , "
                #+"E/$L^{2}$ = "+str("%.5f" % EN)
                #+" , "
                #+"Sample. = "+str(IDISD)+" , "
                #+"Init. config. = "+str(BTNUM)
                ,x=0.485,fontsize=15)
        
        plt.xlabel(r'Iteration count = '+
                str("%d" % inum)+
                r', Error = '+
                str("%.4e" % err))
        
        #ax.axis('off')
        #fig.tight_layout(pad=1.6,h_pad=1.6,w_pad=1.6)

        fig.savefig("../plot/esbc"+
                "_L_"+
                STR_L+
                "_DLT_"+
                STR_DELTA+
                "_ALP_"+
                STR_ALPHA+
                "_DNM_"+
                str("%06d" % IDISD)+
                str("%06d" % BTNUM)+
                ".pdf",
                bbox_inches='tight',
                )
        plt.close('all')


for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite '+
                    '-dPDSETTINGS=/prepress -sOutputFile='+
                '../plot/MF_CONFIG_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/esbc'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '* ',shell=True)
            subprocess.call('rm ../plot/esbc'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '*',shell=True)

