#!/usr/bin/env python

# A matplotlib plotting script to plot a simple x-y plot with yerrorbars
import os
os.chdir("../zero_field_defect/DATA")
import sys
import subprocess
import re
import glob
import numpy as np
import matplotlib.figure as figure
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# The tex style commands
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

############## the dataset range ################

CLR=['red','royalblue','lawngreen']

# Loading data
l, delta, alpx, s_x, s_y, s_z = np.loadtxt('LISTDATA.dat', usecols=(0,1,2,4,5,6), unpack=True)
qs_x, qs_y, qs_z = np.loadtxt('LISTDATA.dat', usecols=(7,8,9), unpack=True)

# Getting unique data arrays

data_len = len(l)
L=np.unique(l)
LSET=len(L)
ALPX=np.unique(alpx)
ASET=len(ALPX)
DELTA=np.unique(delta)
DSET=len(DELTA)

#### Principle data array ####

data_C = np.zeros((LSET,DSET,ASET),dtype=np.float)
data_Q = np.zeros((LSET,DSET,ASET),dtype=np.float)
data = np.zeros((LSET,DSET,ASET),dtype=np.float)

S = 0.5

##### Defining fitting functions ####

def power(x,a):

    return a*x

def power_2(x,a,b,c):

    return a + b*x

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
    data[ln,dn,an] = S*np.linalg.norm(s+(1.0/S)*qs)

##### Figure generation ##########

subprocess.call('mkdir -p PLOT',shell=True)

DLALPHA = np.zeros((DSET,ASET),dtype=np.float)

fyl = open('UMOM.dat','w+')

print >> fyl, "ALPHA", "DELTA", "CL", "QQ"

fyl2 = open('L_SCALE.dat','w+')

print >> fyl2, "ALPHA", "DELTA", "m0", "m1", "m2"

fyl3 = open('DELTA_SCALE.dat','w+')

print >> fyl3, "ALPHA", "m0", "m1", "m2"

for an in range(0,ASET):
    for dn in range(0,DSET):
        
        w,h = figure.figaspect(1.0)
        fig = plt.figure(figsize=(w,h))
        ax = fig.add_axes([0.26,0.15,0.685,0.8])

        dmx_C = data_C[:,dn,an] 
        dmx_Q = data_Q[:,dn,an] 
        dmx = data[:,dn,an]

        invL = L**(-1.0)

        ax.plot(invL,dmx,
                ls='None',
                marker='s',
                ms=12,
                mew=2,
                mfc='None',
                mec='red',
                zorder=1
                )

        popt_C, pcov_C = curve_fit(power_2,invL,dmx_C) 
        popt_Q, pcov_Q = curve_fit(power_2,invL,dmx_Q) 
        popt, pcov = curve_fit(power_2,invL,dmx)

        DLALPHA[dn,an] = popt[0]

        print >> fyl, ALPX[an], DELTA[dn], \
                str("%.4f" % popt_C[0]), str("%.4f" % popt_Q[0])
        
        print >> fyl2, ALPX[an], DELTA[dn], \
                str("%.4f" % popt[0]), str("%.4f" % popt[1]), \
                str("%.4f" % popt[2])

        invLP = np.linspace(0.0,0.1,num=200)
        
        ax.plot(invLP,power_2(invLP, *popt),
                color='royalblue',
                linestyle='--',
                lw=4,
                label=r'$a^{(0)}_L+a^{(1)}_L/L+a^{(2)}_L/L^2$',
                zorder=2
                )
        ax.xaxis.set_ticks(np.arange(0.0,0.12,0.02))
        plt.xlim([0.0,0.1])
        plt.ylabel(r'$m_{\textrm{imp}}(\delta J = - J)$'
                ,fontsize=30)
        plt.xlabel(r'$1/L$',fontsize=30)
        plt.tick_params(which='both',width=2,labelsize=30)
        plt.tick_params(which='major',length=20)
        plt.tick_params(which='minor',length=10)

        legend=ax.legend(loc='best',
                    fontsize=20,
                    markerscale=2,
                    shadow=True)
        #plt.ylim([10.0**(-2),0.5])
        #plt.xlim([0.0,0.0])
        #plt.ylim([np.amin(dmx),np.amax(dmx)])
        
        fig.savefig("PLOT/UMOM-"
                +"DELTA_"
                +str("%.4f" % DELTA[dn])
                +"_ALPHA_"
                +str("%.4f" % ALPX[an])
                +".pdf",
                #bbox_inches='tight'
                )
        
        plt.close('all')
    
    w,h = figure.figaspect(1.0)
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_axes([0.26,0.15,0.685,0.8])

    ax.plot(DELTA,DLALPHA[:,an],
            ls='None',
            marker='s',
            ms=12,
            mew=2,
            mfc='None',
            mec='red',
            zorder=1
            )
    popt, pcov = curve_fit(power,DELTA,DLALPHA[:,an])

    DELTAX = np.linspace(np.amin(DELTA),np.amax(DELTA),num=200)
    
    print >> fyl3, ALPX[an], str("%.4f" % popt[0])

    ax.plot(DELTAX,power(DELTAX, *popt),
            color='royalblue',
            linestyle='--',
            lw=4,
            zorder=2
            )
    
    plt.ylabel(r'$m_{\textrm{imp}}$',fontsize=30)
    plt.xlabel(r'$\delta J$',fontsize=30)

    plt.tick_params(which='both',width=2,labelsize=30)
    plt.tick_params(which='major',length=20)
    plt.tick_params(which='minor',length=10)

    """
    legend=ax.legend(loc='best',
                fontsize=20,
                markerscale=2,
                shadow=True)
    """
    fig.savefig("PLOT/UMOM-DELTA_VS_ALPHA"
            +".pdf",
            #bbox_inches='tight'
            )
    
    plt.close('all')

