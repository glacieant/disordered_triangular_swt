#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as figure
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm


## The tex style commands
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')

#loading data
sx, sy, sz = np.loadtxt('out',usecols=(0,1,2),unpack=True)
sx0, sy0, sz0 = np.loadtxt('out0',usecols=(0,1,2),unpack=True)

## picturing the output data

X = np.zeros(1)
Y = np.zeros(1)
Z = np.zeros(1)

a = 1.0

# three diffrent translation vector

avec = np.array([
    [-((3.0)**0.5)/2,0.5],
    [0.0,1.0],
    [((3.0)**0.5)/2,0.5],
    [((3.0)**0.5)/2,-0.5],
    [0.0,-1.0],
    [-((3.0)**0.5)/2,-0.5],
    ])*a/2

L = 24
N = L*L

# laying out the lattice skeleton
ax = a*(3.0**0.5)/2.0
ay = a/2.0

X = np.zeros((L,L))
Y = np.zeros((L,L))
for i in range(0,L):
    for j in range(0,L):
        #X[i,j]=(i+j*(1.0/2.0))*a +2*a
        #Y[i,j]=j*((3.0)**0.5/2.0)*a + 2*a
        X[i,j]= (i+j)*ax+3*ax
        Y[i,j]= (j-i)*ay

# getting the spin vector

spin = np.column_stack((sx,sy,sz))
spin0 = np.column_stack((sx0,sy0,sz0))

e1 = spin0[0]
if np.linalg.norm(e1) > 10.0**(-5):
    e1 = e1/np.linalg.norm(e1)
    e2 = np.cross(e1,np.cross(e1,spin0[1]))
    e2 = e2/np.linalg.norm(e2)
else:
    e1 = np.array([1.0,0.0,0.0])
    e2 = np.array([0.0,1.0,0.0])

U0 = np.einsum('ij,j->i',spin0,e1).reshape((L,L)).transpose()
V0 = np.einsum('ij,j->i',spin0,e2).reshape((L,L)).transpose()

U = np.einsum('ij,j->i',spin,e1).reshape((L,L)).transpose()
V = np.einsum('ij,j->i',spin,e2).reshape((L,L)).transpose()

# plotting the configuration

w,h = figure.figaspect(1)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1.0/(3.0)**(0.5)])

# picturing the spin orientation
pivot="mid"
width=0.002
headwidth=4
headlength=5

Q=ax.quiver(X,Y,U0,V0,
        color='gray',
        width=width,
        headwidth=headwidth,
        headlength=headlength,
        pivot=pivot,
        angles='xy',
        scale=1,
        scale_units='xy',
        zorder=1,
        alpha=0.5
        )

Q=ax.quiver(X,Y,U,V,
        color='blue',
        width=width,
        headwidth=headwidth,
        headlength=headlength,
        pivot=pivot,
        angles='xy',
        scale=1,
        scale_units='xy',
        zorder=2
        )
ax.axis('off')

fig.savefig("spin_layout.pdf",bbox_inches='tight'
            )
plt.close('all')

