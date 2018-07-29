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

# TU colors

tu_dunkelblau = '#07284A'
tu_grau = '#5F6967'
tu_blau = '#0067A5'
tu_cyan = '#00A3DA'
tu_dunkelgruen = '#008644'
tu_gruen = '#5EB245'
tu_rot = '#DE4A39'
tu_dunkelrot = '#BD252C'

## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

# setting MKL configuration to allways optimise 
os.putenv("MKL_DYNAMIC","FALSE")

CLR = [tu_rot,tu_blau,tu_gruen]

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

##### Figure generation ##########

subprocess.call('mkdir -p PLOT',shell=True)

DLALPHA_C = np.zeros((DSET,ASET),dtype=np.float)
DLALPHA_Q = np.zeros((DSET,ASET),dtype=np.float)

fyl = open('UMOM.dat','w+')

print >> fyl, "ALPHA", "DELTA", "CL", "QQ"

fyl2C = open('L_SCALE_C.dat','w+')

print >> fyl2C, "ALPHA", "DELTA", "a0", "a1", "a2"

fyl2Q = open('L_SCALE_Q.dat','w+')

print >> fyl2Q, "ALPHA", "DELTA", "a0", "a1", "a2"

fyl3C = open('DELTA_SCALE_C.dat','w+')

print >> fyl3C, "ALPHA", "m0/Delta"

fyl3Q = open('DELTA_SCALE_Q.dat','w+')

print >> fyl3Q, "ALPHA", "m1/Delta"


for an in range(0,ASET):
    for dn in range(0,DSET):
        
        w,h = figure.figaspect(0.5)
        afig = plt.figure(figsize=(w,h))
        ax = afig.add_axes([0.26,0.15,0.685,0.8])
        bfig = plt.figure(figsize=(w,h))
        bx = bfig.add_axes([0.26,0.15,0.685,0.8])

        dmx_C = data_C[:,dn,an] 
        dmx_Q = data_Q[:,dn,an] 

        invL = L**(-1.0)

        ax.plot(invL,dmx_C,
                ls='None',
                marker='s',
                ms=12,
                mew=2,
                mfc='None',
                mec=tu_rot,
                zorder=1
                )

        bx.plot(invL,dmx_Q,
                ls='None',
                marker='s',
                ms=12,
                mew=2,
                mfc='None',
                mec=tu_rot,
                zorder=1
                )


        popt_C, pcov_C = curve_fit(power_3,invL,dmx_C) 
        popt_Q, pcov_Q = curve_fit(power_3,invL,dmx_Q) 

        DLALPHA_C[dn,an] = popt_C[0]
        DLALPHA_Q[dn,an] = popt_Q[0]

        print >> fyl, ALPX[an], DELTA[dn], \
                str("%.4f" % popt_C[0]), str("%.4f" % popt_Q[0])
        
        print >> fyl2C, ALPX[an], DELTA[dn], \
                str("%.4f" % popt_C[0]), str("%.4f" % popt_C[1]), \
                str("%.4f" % popt_C[2])

        print >> fyl2Q, ALPX[an], DELTA[dn], \
                str("%.4f" % popt_Q[0]), str("%.4f" % popt_Q[1]), \
                str("%.4f" % popt_Q[2])

        invLP = np.linspace(0.0,0.1,num=200)
        
        ax.plot(invLP,power_3(invLP, *popt_C),
                color=tu_blau,
                linestyle='--',
                lw=4,
                label=r'$a^{(0)}+a^{(1)}/L+a^{(2)}/L^2$',
                zorder=2
                )
        ax.set_ylabel(r'$\delta m^{(0)}(\delta J = - J)$'
                ,fontsize=30)
        ax.set_xlim([0.0,0.1])
        ax.set_xlabel(r'$1/L$',fontsize=30)
        ax.tick_params(which='both',width=2,labelsize=30,direction='in')
        ax.tick_params(which='major',length=20)
        ax.tick_params(which='minor',length=10)
        
        legend=ax.legend(loc='best',
                    fontsize=12,
                    markerscale=1,
                    facecolor='w',
                    edgecolor='k',
                    framealpha=1)
        #plt.ylim([10.0**(-2),0.5])
        #plt.xlim([0.0,0.0])
        #plt.ylim([np.amin(dmx),np.amax(dmx)])
        
        afig.savefig("PLOT/UMOM_C-"
                +"DELTA_"
                +str("%.4f" % DELTA[dn])
                +"_ALPHA_"
                +str("%.4f" % ALPX[an])
                +".pdf",
                bbox_inches='tight'
                )
         
        bx.plot(invLP,power_3(invLP, *popt_Q),
                color=tu_blau,
                linestyle='--',
                lw=4,
                label=r'$a^{(0)}+a^{(1)}/L+a^{(2)}/L^2$',
                zorder=2
                )
        
        bx.set_xlim([0.0,0.1])
        bx.set_ylabel(r'$\delta m^{(1)}(\delta J = - J)$'
                ,fontsize=30)
        bx.set_xlabel(r'$1/L$',fontsize=30)
        bx.tick_params(which='both',width=2,labelsize=30,direction='in')
        bx.tick_params(which='major',length=20)
        bx.tick_params(which='minor',length=10)

        legend=bx.legend(loc='best',
                    fontsize=12,
                    markerscale=1,
                    facecolor='w',
                    edgecolor='k',
                    framealpha=1)
        #plt.ylim([10.0**(-2),0.5])
        #plt.xlim([0.0,0.0])
        #plt.ylim([np.amin(dmx),np.amax(dmx)])
        
        bfig.savefig("PLOT/UMOM_Q"
                +"DELTA_"
                +str("%.4f" % DELTA[dn])
                +"_ALPHA_"
                +str("%.4f" % ALPX[an])
                +".pdf",
                bbox_inches='tight'
                )
        

        plt.close('all')
    
    w,h = figure.figaspect(0.5)
    afig = plt.figure(figsize=(w,h))
    ax = afig.add_axes([0.26,0.15,0.685,0.8])
    bfig = plt.figure(figsize=(w,h))
    bx = bfig.add_axes([0.26,0.15,0.685,0.8])
    cfig = plt.figure(figsize=(w,h))
    cx = cfig.add_axes([0.26,0.15,0.685,0.8])

    ax.plot(DELTA,DLALPHA_C[:,an],
            ls='None',
            marker='s',
            ms=12,
            mew=2,
            mfc='None',
            mec=tu_rot,
            zorder=1
            )

    bx.plot(DELTA,DLALPHA_Q[:,an],
            ls='None',
            marker='s',
            ms=12,
            mew=2,
            mfc='None',
            mec=tu_rot,
            zorder=1
            )

    cx.plot(DELTA,DLALPHA_C[:,an]-DLALPHA_Q[:,an],
            ls='None',
            marker='s',
            ms=12,
            mew=2,
            mfc='None',
            mec=tu_rot,
            zorder=1
            )


    popt_C, pcov_C = curve_fit(power,DELTA[:3],DLALPHA_C[:3,an])
    popt_Q, pcov_Q = curve_fit(power,DELTA[:3],DLALPHA_Q[:3,an])
    popt, pcov = curve_fit(power,DELTA[:3],DLALPHA_C[:3,an]-DLALPHA_Q[:3,an])

    DELTAX = np.linspace(np.amin(DELTA),np.amax(DELTA),num=200)
    
    print >> fyl3C, ALPX[an], str("%.4f" % popt_C[0])
    print >> fyl3Q, ALPX[an], str("%.4f" % popt_Q[0])

    ax.plot(DELTAX,power(DELTAX, *popt_C),
            color=tu_blau,
            linestyle='--',
            lw=4,
            zorder=2
            )
    
    ax.set_ylabel(r'$m_{tot}/S$',fontsize=30)
    ax.set_xlabel(r'$\delta J/J$',fontsize=30)
    ax.tick_params(which='both',width=2,
            labelsize=30,direction='in',
            bottom=True,top=True,
            left=True,right=True)
    ax.tick_params(which='major',length=20)
    ax.tick_params(which='minor',length=10)
    ax.xaxis.set_major_locator(plt.FixedLocator(locs=[0.0,0.4,0.8]))

    """
    legend=ax.legend(loc='best',
                fontsize=20,
                markerscale=2,
                shadow=True)
    """
    afig.savefig("PLOT/UMOM_C-DELTA_VS_ALPHA"
            +".pdf",
            bbox_inches='tight'
            )
 
    bx.plot(DELTAX,power(DELTAX, *popt_Q),
            color=tu_blau,
            linestyle='--',
            lw=4,
            zorder=2
            )
    
    bx.set_ylabel(r'$\delta m^{(1)}$',fontsize=30)
    bx.set_xlabel(r'$\delta J/J$',fontsize=30)
    bx.tick_params(which='both',width=2,labelsize=30,direction='in')
    bx.tick_params(which='major',length=20)
    bx.tick_params(which='minor',length=10)
    bx.xaxis.set_major_locator(plt.FixedLocator(locs=[0.0,0.4,0.8]))
    #bx.yaxis.set_major_locator(plt.FixedLocator(locs=[0.0,0.2,0.4]))

    """
    legend=ax.legend(loc='best',
                fontsize=20,
                markerscale=2,
                shadow=True)
    """
    bfig.savefig("PLOT/UMOM_Q-DELTA_VS_ALPHA"
            +".pdf",
            bbox_inches='tight'
            )
    
  
    cx.plot(DELTAX,power(DELTAX, *popt),
            color=tu_blau,
            linestyle='--',
            lw=4,
            zorder=2
            )
    
    cx.set_ylabel(r'$m_{tot}/S$',fontsize=30)
    cx.set_xlabel(r'$\delta J/J$',fontsize=30)
    cx.tick_params(which='both',width=2,
            labelsize=30,direction='in',
            bottom=True,top=True,
            left=True,right=True)
    cx.tick_params(which='major',length=20)
    cx.tick_params(which='minor',length=10)
    cx.xaxis.set_major_locator(plt.FixedLocator(locs=[0.0,0.4,0.8]))
    #bx.yaxis.set_major_locator(plt.FixedLocator(locs=[0.0,0.2,0.4]))

    """
    legend=ax.legend(loc='best',
                fontsize=20,
                markerscale=2,
                shadow=True)
    """
    cfig.savefig("PLOT/UMOM_DELTA_VS_ALPHA"
            +".pdf",
            bbox_inches='tight'
            )
  
    plt.close('all')

