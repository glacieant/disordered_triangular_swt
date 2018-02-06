#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as figure
from matplotlib.colors import Normalize


## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# plotting the configuration

w,h = figure.figaspect(0.5)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

def energy(a,b,c,Q):

    EX = ((1+a)*np.cos(Q.real)+
            (1+b)*np.cos(0.5*Q.real+0.5*Q.imag*(3.0**0.5))
            +(1+c)*np.cos(-0.5*Q.real+0.5*Q.imag*(3.0**0.5))
            )

    return EX

NX = 500

EN_ARRAY = np.zeros((NX,NX),dtype=np.float)

QX = np.arange(-2*np.pi,2*np.pi+4*np.pi/NX,4*np.pi/NX)
QY = np.arange(-2*np.pi,2*np.pi+4*np.pi/NX,4*np.pi/NX)
Q = QX[:,None] + QY[None,:]*1.0j

a = 0
b = 0
c = 0

EN_ARRAY = energy(a,b,c,Q)

MINCLIP = -1.5
MAXCLIP = 0.0

cax = ax.imshow(EN_ARRAY,
        norm=Normalize(vmin=MINCLIP,vmax=MAXCLIP,clip=True),
        interpolation='nearest'
        )

fig.savefig("anisotropy_Q.pdf"
        ,bbox_inches='tight'
        )
plt.close('all')

