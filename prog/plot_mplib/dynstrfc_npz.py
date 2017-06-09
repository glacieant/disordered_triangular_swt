#!/usr/bin/env python

## The import commands
import os
os.chdir("../high_field/out/")
from math import log
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

## The tex style commands
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

## loading the dataset

DXDATA = np.load("strfc.npz")

L = DXDATA['LSYS'][()]
DELTA = DXDATA['DELTA']
ALPHA = DXDATA['ALPHA']
ITERDISD = DXDATA['ITERDISD'][()]
MINF = DXDATA['MINF'][()]
MAXF = DXDATA['MAXF'][()]
FGRID = DXDATA['FGRID'][()]
FWDTH = DXDATA['FWDTH'][()]
SFC = DXDATA['dynstrc']

## processing data for momentum space layout

for d in range(0,len(DELTA)):

    for a in range(0,len(ALPHA)):

        # processing the BZ cut
        # the reciprocal vectors are
        # (2*pi/a,-2*pi/sqrt(3)/a) and
        # (0,4*pi/sqrt(3)/a)
        HCUT = [2,1]
        VCUT = [1,0]
        XMAXMO = [2.0/3.0,1.0/3.0]
        YMAXMO = [1.0,1.0]
        RMAX = [0.0,0.0]

        for i in range(0,2):

            if VCUT[i] == 0 and HCUT[i] != 0:
                RMAX[i] = int((XMAXMO[i]*L)/HCUT[i])
            elif VCUT[i] != 0 and HCUT[i] == 0:
                RMAX[i] = int((YMAXMO[i]*L)/VCUT[i])
            else:
                NX = int((XMAXMO[i]*L)/HCUT[i])
                NY = int((YMAXMO[i]*L)/VCUT[i])
                if NX >= NY:
                    RMAX[i] = NY
                else:
                    RMAX[i] = NX
                
        SFC_CUT = np.zeros((FGRID,np.sum(RMAX)))

        for w in range(0,FGRID):
            for k in range(0,RMAX[0]):
                SFC_CUT[w,k] = SFC[d,a,w,(RMAX[0]-k-1)*VCUT[0],(RMAX[0]-k-1)*HCUT[0]]
            for k in range(0,RMAX[1]):
                SFC_CUT[w,k+RMAX[0]] = SFC[d,a,w,k*VCUT[1],k*HCUT[1]]

        # performing the disorder average
        SFC_CUT *= (1.0/ITERDISD)


        # Setting up the clipping of the density plot 
        MAXCLIP=np.amax(SFC_CUT)
        MINCLIP=MAXCLIP/(10**4)
        mxc=np.log10(MAXCLIP)
        mnc=np.log10(MINCLIP)
        CLIPTICK=np.logspace(mnc,mxc,num=4,endpoint=True,base=10.0)

        #Making plot with vertical colorbar
        fig, ax = plt.subplots()

        #Imshow plotting the data
        #cax=ax.imshow(np.flipud(zmat), norm=Normalize(vmin=MINCLIP, vmax=MAXCLIP, clip=False),
        cax=ax.imshow(np.flipud(SFC_CUT), norm=LogNorm(vmin=MINCLIP, vmax=MAXCLIP, clip=False),
                interpolation='nearest',
                extent=(0.0, 1.0, MINF, MAXF),
                cmap=cm.afmhot ,aspect='auto')
        cbar=fig.colorbar(cax,shrink=0.5,
        #        ticks=CLIPTICK,
                format='%.0e')
        xtickpos=[0.0,0.5,1.0]
        xticklabel=[r'M',r'$\Gamma$',r'K']
        plt.xticks(xtickpos,xticklabel)
        plt.xlabel(r'\textbf{q}', fontsize=18)
        plt.ylabel(r'$\omega$/J', fontsize=18)
        plt.title(r'S(\textbf{q},$\omega$) at $\Delta$='
                +str("%.3f" % DELTA[d])+r', $\alpha$='+str("%.3f" % ALPHA[a]), 
                fontsize=20,y=1.09)
        fig.tight_layout(pad=1.0,h_pad=1.0,w_pad=1.0)
        fig.savefig('GRID-w-k-'+str("%.3f" % DELTA[d])+'-'+str("%.3f" % ALPHA[a])+'.pdf',format='pdf')
        plt.close('all')
