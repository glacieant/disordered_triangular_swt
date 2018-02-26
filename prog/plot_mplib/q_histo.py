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
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import matplotlib.figure as figure
import matplotlib.mlab as mlab
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm


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

# importing lattice geometry system
import lattice_map as lmap
ZCO = 12

X = np.zeros(1)
Y = np.zeros(1)
Z = np.zeros(1)

a = 1.0

# three diffrent translation vector

"""
avec = np.array([[-0.5,((3.0)**0.5)/2],
    [0.5,((3.0)**0.5)/2],
    [1.0,0.0],
    [-1.5,((3.0)**0.5)/2],
    [0.0,(3.0**0.5)],
    [1.5,((3.0)**0.5)/2]])*a
"""

evec = np.array([
    [-((3.0)**0.5)/2,0.5],
    [0.0,1.0],
    [((3.0)**0.5)/2,0.5],
    [((3.0)**0.5)/2,-0.5],
    [0.0,-1.0],
    [-((3.0)**0.5)/2,-0.5],
    ])

avec = evec*a/2

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    STR_L = match.group(1)
    L = int(STR_L)
    STR_DELTA = match.group(2)
    STR_ALPHA = match.group(3)

    FNDATA = np.load(fname)

    J = FNDATA['J']
    spin = FNDATA['spin']
    spin0 = FNDATA['spin0']

    N = L**2

    nbr = np.zeros((N,ZCO),dtype=np.int)
    lmap.lattice_map(L,nbr)

    if np.abs(np.float(STR_DELTA)) < 10.0**(-4) \
            and np.abs(np.float(STR_ALPHA)) < 10.0**(-4) :
        
        US = np.zeros(N,dtype=np.float)
        VS = np.zeros(N,dtype=np.float)
        
        for i in range(0,N):
            p = 0
            pp = 3
            q = 2
            qq = 5
            s1 = (np.arccos(np.dot(spin[i],spin[nbr[i,p]]))+
                    np.arccos(np.dot(spin[i],spin[nbr[i,pp]])))/2
            s2 = (np.arccos(np.dot(spin[i],spin[nbr[i,q]]))+
                    np.arccos(np.dot(spin[i],spin[nbr[i,qq]])))/2
            US[i] = s1
            VS[i] = (2.0*s2 - US[i])/np.sqrt(3)

        UZ = US[i]
        VZ = VS[i]

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
    spin = FNDATA['spin']
    spin0 = FNDATA['spin0']

    N = L**2

    nbr = np.zeros((N,ZCO),dtype=np.int)
    lmap.lattice_map(L,nbr)

    ax = a*(3.0**0.5)/2.0
    ay = a/2.0

    # laying out the lattice skeleton
    if len(X) != N:

        X = np.zeros((L,L))
        Y = np.zeros((L,L))
        for i in range(0,L):
            for j in range(0,L):
                #X[i,j]=(i+j*(1.0/2.0))*a +2*a
                #Y[i,j]=j*((3.0)**0.5/2.0)*a + 2*a
                X[i,j]= (i+j)*ax+3*ax
                Y[i,j]= (j-i)*ay
    """
    # defining the plane for the spins
    e1 = spin[0]
    if np.linalg.norm(e1) > 10.0**(-5):
        e1 = e1/np.linalg.norm(e1)
        e2 = np.cross(e1,np.cross(e1,spin[1]))
        e2 = e2/np.linalg.norm(e2)
    else:
        e1 = np.array([1.0,0.0,0.0])
        e2 = np.array([0.0,1.0,0.0])

    U = 0.5*np.einsum('ij,j->i',spin,e1).reshape((L,L)).transpose()
    V = 0.5*np.einsum('ij,j->i',spin,e2).reshape((L,L)).transpose()
    """

    US = np.zeros(N,dtype=np.float)
    VS = np.zeros(N,dtype=np.float)
    
    for i in range(0,N):
        p = 0
        pp = 3
        q = 2
        qq = 5
        s1 = (np.arccos(np.dot(spin[i],spin[nbr[i,p]]))+
                np.arccos(np.dot(spin[i],spin[nbr[i,pp]])))/2
        s2 = (np.arccos(np.dot(spin[i],spin[nbr[i,q]]))+
                np.arccos(np.dot(spin[i],spin[nbr[i,qq]])))/2
        US[i] = s1
        VS[i] = (2.0*s2 - US[i])/np.sqrt(3)

    U = US[i] - UZ
    V = VS[i] - VZ

    # plotting the histogram

    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])

    # the histogram of the dara
    ax.hist(U,bins=50,normed=1,facecolor='royalblue')

    plt.title(r'Histogram of $\Delta Q_x$')
    fig.savefig("../plot/qhist_x"+
            "_L_"+
            STR_L+
            "_DLT_"+
            STR_DELTA+
            "_ALP_"+
            STR_ALPHA+
            "_DNM_"+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf"
            ,bbox_inches='tight'
            )
    plt.close('all')

    # plotting the histogram

    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0,0,1,1])

    # the histogram of the dara
    ax.hist(V,bins=50,normed=1,facecolor='royalblue')

    plt.title(r'Histogram of $\Delta Q_y$')
    fig.savefig("../plot/qhist_y"+
            "_L_"+
            STR_L+
            "_DLT_"+
            STR_DELTA+
            "_ALP_"+
            STR_ALPHA+
            "_DNM_"+
            str("%06d" % IDISD)+
            str("%06d" % BTNUM)+
            ".pdf"
            ,bbox_inches='tight'
            )
    plt.close('all')

for L in set(zip(*hshchar)[0]):
    for DLT in set(zip(*hshchar)[1]):
        for ALP in set(zip(*hshchar)[2]):
            subprocess.call('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite '+
                    '-dPDSETTINGS=/prepress -sOutputFile='+
                '../plot/QXHIST_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/qhist_x'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '* ',shell=True)
            subprocess.call('rm ../plot/qhist_x'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '*',shell=True)
            subprocess.call('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite '+
                    '-dPDSETTINGS=/prepress -sOutputFile='+
                '../plot/QYHIST_L_'+
                L+
                '_DELTA_'+
                DLT+
                '_ALPHA_'+
                ALP+
                '.pdf '+
                '../plot/qhist_y'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '* ',shell=True)
            subprocess.call('rm ../plot/qhist_y'+
                "_L_"+
                L+
                "_DLT_"+
                DLT+
                "_ALP_"+
                ALP+
                "_DNM_"+
                '*',shell=True)

