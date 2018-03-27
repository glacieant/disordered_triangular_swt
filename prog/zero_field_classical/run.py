#!/usr/bin/env python
import multiprocessing as mp

##################################################
# This is the wrapper for the computation of the #
# triangular lattice large-N mean-field program  #
##################################################

# system size
LSYS = [30,36] 
# co-ordination number of the lattice
ZCO = 12
# disorder iteration number 
ITERDISD = [10,10]
# initial angle fluctuation
ANGVAR = [0.5,0.5]
# number of bootstrapping
BOOTNUM = [10,10]
# a global tolerance value
GTOL = 10.0**(-10) 
# the disorder amplitude
DELTA = [0.8] 
# the ratio between nearest and next 
# nearest couplings
ALPHA = [0.0] 
# maximum iteration for classical algorithm
CLNUM = 10**7

# checking the lengths of simulation arrays #

if len(LSYS) == len(ITERDISD) \
        == len(ANGVAR) == len(BOOTNUM):
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
if not os.path.exists("src/tri_imp.py"):
    print "Main executables could not be found"
    quit()
if not os.path.exists("out"):
    os.makedirs("out")
if not os.path.exists("out/data"):
    os.makedirs("out/data")

# Setting MKL_DYNAMIC configuration

os.putenv("MKL_DYNAMIC","FALSE")

# Importing the main script

import src.tri_imp as tmi

# fixing the simulation parameter tuple #

ivar = collections.namedtuple('ivar',
        'LSYS ZCO ITERDISD \
                ANGVAR BOOTNUM \
                GTOL DELTA ALPHA \
                CLNUM DNMR')

# generate paramter range

isize = len(LSYS)
dsize = len(DELTA)
asize = len(ALPHA)

# the parameter function

def const(DNMR):

    carray = range(len(DNMR))
    for DNM in DNMR:
        i = (DNM%(isize*dsize))%isize 
        j = (DNM%(isize*dsize))/isize 
        k = DNM/(isize*dsize)

        carray[DNM] = ivar(LSYS = LSYS[i],
                ZCO = ZCO,
                ITERDISD = ITERDISD[i],
                ANGVAR = ANGVAR[i],
                BOOTNUM = BOOTNUM[i],
                GTOL = GTOL,
                DELTA = DELTA[j],
                ALPHA = ALPHA[k],
                CLNUM = CLNUM,
                DNMR = DNM)

    return carray

DNX = range(isize*dsize*asize)

NUM_CORES = mp.cpu_count()

# batch processing
pool = mp.Pool(processes=NUM_CORES)
pool.map(tmi.main,const(DNX))
