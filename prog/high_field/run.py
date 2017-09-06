#!/usr/bin/env python

##################################################
# This is the wrapper for the computation of the #
# triangular lattice large-h SWT approximation   #
##################################################

# system size
LSYS = [60] 
# co-ordination number of the lattice
ZCO = 12 
# disorder iteration number 
ITERDISD = [20]
# magnitude of the external field
HFIELD = [0]
# the disorder amplitude
DELTA = [0.0,0.2,0.55,0.65,0.8,0.9,0.99] 
# the ratio between nearest and next 
# nearest couplings
ALPHA = [0.0,0.025,0.05,0.075,0.1,0.125] 
# a global tolerance value
GTOL = 10.0**(-8.0) 
# the spin magnitude
S = 0.5
# the minimum frequency
MINF = -6.0
# the maximum frequency
MAXF = 0.0
# the frequency grid 
FGRID = 400
# width of the delta function
GWDTH = 0.01

# checking the lengths of simulation arrays #

if len(LSYS) == len(ITERDISD) \
        == len(HFIELD):
            pass
else:
    print "LSYS/ITERDISD/HFIELD input wrong"
    quit()

# importing processing tools #

import os
import subprocess
import collections

# checking if the directory is intact #

if not os.path.exists("src"):
    print "Executables could not be found."
    quit()
if not os.path.exists("src/tri_hfswt.py"):
    print "Main executables could not be found."
    quit()
if not os.path.exists("out"):
    os.makedirs("out")

# Setting MKL_DYNAMIC configuration

os.putenv("MKL_DYNAMIC","FALSE")

# Importing the main script

import src.tri_hfswt as thf

# fixing the simulation parameter tuple #

ivar = collections.namedtuple('ivar',
        'LSYS ZCO ITERDISD \
                HFIELD DELTA \
                ALPHA GTOL \
                S DNMR \
                MINF MAXF \
                FGRID GWDTH')

isize = len(LSYS)

# denominator descriptor
DNMR = 0
for i in range(0,isize):

    const = ivar(LSYS = LSYS[i],
            ZCO = ZCO,
            ITERDISD = ITERDISD[i],
            HFIELD = HFIELD[i],
            DELTA = DELTA,
            ALPHA = ALPHA,
            GTOL = GTOL,
            S = S,
            DNMR = DNMR,
            MINF = MINF,
            MAXF = MAXF,
            FGRID = FGRID,
            GWDTH = GWDTH)

    thf.main(const)
    DNMR += 1
