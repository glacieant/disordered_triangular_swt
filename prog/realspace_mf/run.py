#!/usr/bin/env python

##################################################
# This is the wrapper for the computation of the #
# triangular lattice large-N mean-field program  #
##################################################

# system size
LSYS = [24] 
# co-ordination number of the lattice
ZCO = 12
# disorder iteration number 
ITERDISD = [1]
# initial angle fluctuation
ANGVAR = [0.5]
# initial bond fluctuation
BONDVAR = [0.5]
# number of bootstrapping
BOOTNUM = [1]
# maximum number of iteration for 
# the mean-field simulation
SIMSZE = [1000] 
# tolerance for the mean-field routine
MFTOL = [10.0**(-5)]
# a global tolerance value
GTOL = 10.0**(-10) 
# the disorder amplitude
DELTA = [0.0] 
# the ratio between nearest and next 
# nearest couplings
ALPHA = [0.125] 
# tempering the fermi function
TEMP = [0.01] 
# the fudge parameter
XPAR = [0.5] 
# update ratio for magnetisation and bonds
UPAR = [[0.3,0.7],[0.3,0.7]]
# maximum iteration for classical algorithm
CLNUM = 10**6

# checking the lengths of simulation arrays #

if len(LSYS) == len(ITERDISD) \
        == len(ANGVAR) == len(BONDVAR) \
        == len(BOOTNUM) == len(SIMSZE) \
        == len(MFTOL) == len(TEMP) \
        == len(XPAR):
            pass
else:
    print "System size/simulation array length wrong"
    quit()

# importing processing tools #

import os
import subprocess
import collections

# checking if the directory is intact #

if not os.path.exists("src"):
    print "Executables could not be found"
    quit()
if not os.path.exists("src/tri_mf.py"):
    print "Main executables could not be found"
    quit()
if not os.path.exists("out"):
    os.makedirs("out")
if not os.path.exists("out/data"):
    os.makedirs("out/data")

# Setting MKL_DYNAMIC configuration

os.putenv("MKL_DYNAMIC","FALSE")

# Importing the main script

import src.tri_mf as tmf

# fixing the simulation parameter tuple #

ivar = collections.namedtuple('ivar',
        'LSYS ZCO ITERDISD \
                ANGVAR BONDVAR BOOTNUM \
                SIMSZE MFTOL GTOL \
                DELTA ALPHA TEMP \
                XPAR UPAR CLNUM DNMR')

isize = len(LSYS)
dsize = len(DELTA)
asize = len(ALPHA)

# denominator descriptor
DNMR = 0
for i in range(0,isize):

    for j in range(0,dsize):

        for k in range(0,asize):

            const = ivar(LSYS = LSYS[i],
                    ZCO = ZCO,
                    ITERDISD = ITERDISD[i],
                    ANGVAR = ANGVAR[i],
                    BONDVAR = BONDVAR[i],
                    BOOTNUM = BOOTNUM[i],
                    SIMSZE = SIMSZE[i],
                    MFTOL = MFTOL[i],
                    GTOL = GTOL,
                    DELTA = DELTA[j],
                    ALPHA = ALPHA[k],
                    TEMP = TEMP[i],
                    XPAR = XPAR[i],
                    UPAR = UPAR,
                    CLNUM = CLNUM,
                    DNMR = DNMR)

            tmf.main(const)
            DNMR += 1
