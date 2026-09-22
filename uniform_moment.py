#!/usr/bin/env python

# A matplotlib plotting script to plot a simple x-y plot with yerrorbars
import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# TU colors

tu_dunkelblau = '#07284A'
tu_grau = '#5F6967'
tu_blau = '#0067A5'
tu_cyan = '#00A3DA'
tu_dunkelgruen = '#008644'
tu_gruen = '#5EB245'
tu_rot = '#DE4A39'
tu_dunkelrot = '#BD252C'

## mathtext style
plt.rcParams['mathtext.fontset']='cm'
plt.rc('font',family='serif')

DELTA_FIT = int(sys.argv[1])
S = float([ln.split()[2] for ln in open('../INC/trilattice_sw.h') if ln.startswith('#define S ')][0])

# Loading data
l, delta, alpx, s_x, s_y, s_z = np.loadtxt('LISTDATA.txt', usecols=(0,1,26,34,35,36), unpack=True)
qs_x, qs_y, qs_z = np.loadtxt('LISTDATA.txt', usecols=(37,38,39), unpack=True)

# Getting unique data arrays

data_len = len(l)
L=np.unique(l)
LSET=len(L)
ALPX=np.unique(alpx)
ASET=len(ALPX)
DELTA=np.unique(delta)
DSET=len(DELTA)

#### Principle data array ####

data_C = np.zeros((LSET,DSET,ASET),dtype=float)
data_Q = np.zeros((LSET,DSET,ASET),dtype=float)
data_T = np.zeros((LSET,DSET,ASET),dtype=float)

##### Defining fitting functions ####

def power(x,a):

    return np.abs(a)*x

def power_2(x,a,b):

    return a + b*x

def power_3(x,a,b,c):

    return a + b*x + c*(x**2)

def invd(x,a,b):

    return a/x + b

###### looping through the dataset ##############

for g in range(0,data_len):

    ln = np.where(L==l[g])[0].item(0)
    dn = np.where(DELTA==delta[g])[0].item(0)
    an = np.where(ALPX==alpx[g])[0].item(0)

    s = np.array([s_x[g],s_y[g],s_z[g]])
    qs = np.array([qs_x[g],qs_y[g],qs_z[g]])


    data_C[ln,dn,an] = np.linalg.norm(s)
    data_Q[ln,dn,an] = np.linalg.norm(qs)
    data_T[ln,dn,an] = np.linalg.norm(S*s+qs)

##### Figure generation ##########

subprocess.call('mkdir -p PLOT',shell=True)

def style(ax,xlabel,ylabel):

    ax.tick_params(which='major',width=2,
                   labelsize=22,direction='in',
                   length=10,
                   bottom=True,top=False,
                   left=True,right=False)
    ax.set_xlabel(xlabel,fontsize=22)
    ax.set_ylabel(ylabel,fontsize=22)
    ax.set_box_aspect(1)

def fitlegend(ax):

    legend=ax.legend(loc='best',
                fontsize=18,
                facecolor='w',
                edgecolor='k',
                framealpha=1)

DLALPHA_C = np.zeros((DSET,ASET),dtype=float)
DLALPHA_Q = np.zeros((DSET,ASET),dtype=float)
DLALPHA_T = np.zeros((DSET,ASET),dtype=float)

fyl = open('UMOM.txt','w+')

print("ALPHA", "DELTA", "CL", "QQ", file=fyl)

fyl2C = open('L_SCALE_C.txt','w+')

print("ALPHA", "DELTA", "a0", "a1", "a2", file=fyl2C)

fyl2Q = open('L_SCALE_Q.txt','w+')

print("ALPHA", "DELTA", "a0", "a1", "a2", file=fyl2Q)

fyl2T = open('L_SCALE_TOT.txt','w+')

print("ALPHA", "DELTA", "a0", "a1", "a2", file=fyl2T)

fyl3C = open('DELTA_SCALE_C.txt','w+')

print("ALPHA", "m0/Delta", file=fyl3C)

fyl3Q = open('DELTA_SCALE_Q.txt','w+')

print("ALPHA", "m1/Delta", file=fyl3Q)

fyl3T = open('DELTA_SCALE_TOT.txt','w+')

print("ALPHA", "mtot/Delta", file=fyl3T)

FITLABEL = r'$a^{(0)}+a^{(1)}/L+a^{(2)}/L^2$'

for an in range(0,ASET):
    for dn in range(0,DSET):

        invL = L**(-1.0)
        invLP = np.linspace(0.0,np.amax(invL),num=200)

        for data, name, ylabel, store, tfyl in [
                (data_C, "PLOT/UMOM_C-", r'$\delta m^{(0)}(\delta J/J = -%g)$' % DELTA[dn], DLALPHA_C, fyl2C),
                (data_Q, "PLOT/UMOM_Q", r'$\delta m^{(1)}(\delta J/J = -%g)$' % DELTA[dn], DLALPHA_Q, fyl2Q),
                (data_T, "PLOT/UMOM_TOT-", r'$m_{tot}(\delta J/J = -%g)$' % DELTA[dn], DLALPHA_T, fyl2T)]:

            afig, ax = plt.subplots()

            dmx = data[:,dn,an]

            ax.plot(invL,dmx,
                    ls='None',
                    marker='s',
                    ms=12,
                    mew=2,
                    mfc='None',
                    mec=tu_rot,
                    zorder=1
                    )

            popt, pcov = curve_fit(power_3,invL,dmx)

            store[dn,an] = popt[0]

            print(ALPX[an], DELTA[dn], \
                    str("%.4f" % popt[0]), str("%.4f" % popt[1]), \
                    str("%.4f" % popt[2]), file=tfyl)

            ax.plot(invLP,power_3(invLP, *popt),
                    color=tu_blau,
                    linestyle='--',
                    lw=4,
                    label=FITLABEL,
                    zorder=2
                    )
            ax.set_xlim(left=0.0)
            style(ax,r'$1/L$',ylabel)
            fitlegend(ax)

            afig.savefig(name
                    +"DELTA_"
                    +str("%.4f" % DELTA[dn])
                    +"_ALPHA_"
                    +str("%.4f" % ALPX[an])
                    +".pdf",
                    bbox_inches='tight',transparent=True
                    )

            plt.close('all')

        print(ALPX[an], DELTA[dn], \
                str("%.4f" % DLALPHA_C[dn,an]), str("%.4f" % DLALPHA_Q[dn,an]), file=fyl)

    DELTAX = np.linspace(np.amin(DELTA),np.amax(DELTA),num=200)

    for vals, name, ylabel, tfyl in [
            (DLALPHA_C[:,an], "PLOT/UMOM_C-DELTA_VS_ALPHA", r'$m_{tot}/S$', fyl3C),
            (DLALPHA_Q[:,an], "PLOT/UMOM_Q-DELTA_VS_ALPHA", r'$\delta m^{(1)}$', fyl3Q),
            (DLALPHA_C[:,an]-DLALPHA_Q[:,an], "PLOT/UMOM_DELTA_VS_ALPHA", r'$m_{tot}/S$', None),
            (DLALPHA_T[:,an], "PLOT/UMOM_TOT-DELTA_VS_ALPHA", r'$m_{tot}$', fyl3T)]:

        afig, ax = plt.subplots()

        ax.plot(DELTA,vals,
                ls='None',
                marker='s',
                ms=12,
                mew=2,
                mfc='None',
                mec=tu_rot,
                zorder=1
                )

        popt, pcov = curve_fit(power,DELTA[:DELTA_FIT],vals[:DELTA_FIT])

        if tfyl is not None:
            print(ALPX[an], str("%.4f" % popt[0]), file=tfyl)

        ax.plot(DELTAX,power(DELTAX, *popt),
                color=tu_blau,
                linestyle='--',
                lw=4,
                zorder=2
                )

        style(ax,r'$-\delta J/J$',ylabel)

        afig.savefig(name
                +".pdf",
                bbox_inches='tight',transparent=True
                )

        plt.close('all')

for f in [fyl, fyl2C, fyl2Q, fyl2T, fyl3C, fyl3Q, fyl3T]:
    f.close()
