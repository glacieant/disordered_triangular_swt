#!/usr/bin/env python
import multiprocessing as mp

##################################################
# This is the wrapper for the computation of the #
# triangular lattice large-N mean-field program  #
##################################################

# system size
LSYS = [48,54,60] 
# co-ordination number of the lattice
ZCO = 12
# disorder iteration number 
ITERDISD = [50,50,50]
# initial angle fluctuation
ANGVAR = [1.0,1.0,1.0]
# number of bootstrapping
BOOTNUM = [50,50,50]
# a global tolerance value
GTOL = 10.0**(-6) 
# bimodal strength
PROB = [0.5] 
# the disorder amplitude
DELTA = [0.2,0.6,0.99] 
# the ratio between nearest and next 
# nearest couplings
ALPHA = [0.0] 
# maximum iteration for classical algorithm
CLNUM = 10**6

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

import src.tri_dis as tmd

# fixing the simulation parameter tuple #

ivar = collections.namedtuple('ivar',
        'LSYS ZCO ITERDISD \
                ANGVAR BOOTNUM \
                GTOL PROB DELTA ALPHA \
                CLNUM DNMR')

# generate paramter range

isize = len(LSYS)
psize = len(PROB)
dsize = len(DELTA)
asize = len(ALPHA)
tdisd = sum(ITERDISD)

# the parameter function

def const(DNMR):

    carray = range(len(DNMR))
    for DNM in DNMR:
        PRB = (DNM%(asize*dsize*psize))%(dsize*psize)%psize
        DLT = (DNM%(asize*dsize*psize))%(dsize*psize)/psize
        ALP = (DNM%(asize*dsize*psize))/(dsize*psize)
        SIM = (DNM/(asize*dsize*psize))

        for LI in range(isize):

            SIMX = sum(ITERDISD[:(LI+1)])
            if SIMX-SIM > 0:
                  
                carray[DNM] = ivar(LSYS = LSYS[LI],
                        ZCO = ZCO,
                        ITERDISD = 1,
                        ANGVAR = ANGVAR[LI],
                        BOOTNUM = BOOTNUM[LI],
                        GTOL = GTOL,
                        PROB = PROB[PRB],
                        DELTA = DELTA[DLT],
                        ALPHA = ALPHA[ALP],
                        CLNUM = CLNUM,
                        DNMR = DNM)

                break

    return carray

DNX = range(psize*dsize*asize*tdisd)

# batch processing
pool = mp.Pool()
pool.map(tmd.main,const(DNX))
