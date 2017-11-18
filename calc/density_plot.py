#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
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
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit


## The tex style commands
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')

# the font styleset
from matplotlib import rcParams
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.family'] = 'serif'


STR = 'sup_6-dipole'
dSTR = np.loadtxt(STR+'.dat')

fig, ax = plt.subplots()

cax = ax.imshow(dSTR,
        origin='lower',
        norm=Normalize(vmin=-0.4,vmax=0.4,clip=False),
        interpolation='nearest',
        cmap=cm.viridis,
        aspect='auto')

cbar = fig.colorbar(cax,shrink=0.5)

plt.suptitle(r'Superposition of Dipoles',x=0.5,fontsize=16)
ax.set_aspect('equal')
fig.tight_layout(pad=2.5,h_pad=2.5,w_pad=2.5)
fig.savefig(STR+'.pdf')
plt.close('all')
