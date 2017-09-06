#!/usr/bin/env python

## The import commands
import os
os.chdir("../high_field/out/")
import subprocess
from math import log
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.stats import skewnorm

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
SFC = DXDATA['dynstrc']

# post processing

# disorder averaging
SFC *= SFC/(ITERDISD*L**2)
# leveling the bottom        
SFC += 10.0**(-10)
# the max value
SFC_MAX = 0.05

# Defining model function to be used to fit 1d data

def gaussian(x,a,b,c):
    
    return a*np.exp(-(x-b)**2/(2.0*c**2))

def skew_gauss(x,alpha,a,b,c):

    return a*skewnorm.pdf(x,alpha,b,np.sqrt(2.0)*c)

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

        # Setting up the clipping of the density plot 
        MAXCLIP=SFC_MAX
        MINCLIP=MAXCLIP/(10.0**5)

        #Making plot with vertical colorbar
        fig, ax = plt.subplots()

        #Imshow plotting the data
        cax=ax.imshow(np.flipud(SFC_CUT), norm=LogNorm(vmin=MINCLIP, vmax=MAXCLIP, clip=True),
                interpolation='nearest',
                extent=(0.0, 1.0, MINF, MAXF),
                cmap=cm.inferno,aspect='auto')
        cbar=fig.colorbar(cax,shrink=0.5,format='%.0e')
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

        #fitting a gaussian to extract a pattern along the high symmetry points

        # GAMMA point

        fig, ax = plt.subplots()

        # data for gaussian fitting
        YDATA = SFC_CUT[:,RMAX[0]-1]
        XDATA = np.linspace(MINF,MAXF,YDATA.size)
        ax.plot(XDATA,YDATA,'b.',label='Data')
        
        # finding the fit
        #popt, pcov = curve_fit(gaussian, XDATA, YDATA,p0=[np.amax(YDATA),-1.5,0.1])
        popt, pcov = curve_fit(skew_gauss, XDATA, YDATA,p0=[-1.25,np.amax(YDATA),-1.25,0.05])
        #sigma = np.absolute(popt[2])
        sigma = np.absolute(popt[3])
        #ax.plot(XDATA,gaussian(XDATA, *popt),'r-',
        ax.plot(XDATA,skew_gauss(XDATA, *popt),'r-',
                #label=r'Gaussian fit ~ $e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ ,'
                label=r'Skew-normal fit around'+'\n'+r'($\frac{x-\mu}{\sqrt{2}\sigma}$) with,'
                +' $\sigma$ = '
                +str("%.4f" % sigma))
        plt.ylim([0,SFC_MAX])
        plt.ylabel(r'$S_{\Gamma}/L^2$', fontsize=18)
        plt.xlabel(r'$\omega$/J', fontsize=18)
        plt.legend(loc='best')
        plt.title(r'$S_{\Gamma}(\omega)$ at $\Delta$='
                +str("%.3f" % DELTA[d])+r', $\alpha$='+str("%.3f" % ALPHA[a]), 
                fontsize=20,y=1.09)
        fig.tight_layout(pad=1.0,h_pad=1.0,w_pad=1.0)
        fig.savefig('GAUSS-GAMMA-'+str("%.3f" % DELTA[d])+'-'+str("%.3f" % ALPHA[a])+'.pdf',format='pdf')
        plt.close('all')
       
        
        # M point

        fig, ax = plt.subplots()

        # data for gaussian fitting
        YDATA = SFC_CUT[:,0]
        XDATA = np.linspace(MINF,MAXF,YDATA.size)
        ax.plot(XDATA,YDATA,'b.',label='Data')
        
        # finding the fit
        #popt, pcov = curve_fit(gaussian, XDATA, YDATA,p0=[np.amax(YDATA),-1.5,0.1])
        popt, pcov = curve_fit(skew_gauss, XDATA, YDATA,p0=[1.25,np.amax(YDATA),-3.85,0.1])
        #sigma = np.absolute(popt[2])
        sigma = np.absolute(popt[3])
        #ax.plot(XDATA,gaussian(XDATA, *popt),'r-',
        ax.plot(XDATA,skew_gauss(XDATA, *popt),'r-',
                #label=r'Gaussian fit ~ $e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ ,'
                label=r'Skew-normal fit around'+'\n'+r'($\frac{x-\mu}{\sqrt{2}\sigma}$) with,'
                +' $\sigma$ = '
                +str("%.4f" % sigma))
        plt.ylim([0,SFC_MAX])
        plt.ylabel(r'$S_{\textbf{M}}/L^2$', fontsize=18)
        plt.xlabel(r'$\omega$/J', fontsize=18)
        plt.legend(loc='best')
        plt.title(r'$S_{\textbf{M}}$ at $\Delta$='
                +str("%.3f" % DELTA[d])+r', $\alpha$='+str("%.3f" % ALPHA[a]), 
                fontsize=20,y=1.09)
        fig.tight_layout(pad=1.0,h_pad=1.0,w_pad=1.0)
        fig.savefig('GAUSS-M-'+str("%.3f" % DELTA[d])+'-'+str("%.3f" % ALPHA[a])+'.pdf',format='pdf')
        plt.close('all')



# merging data, cleaning directory
subprocess.call('pdfunite GRID-* strfc.pdf',shell=True)
subprocess.call('rm GRID-*',shell=True)
subprocess.call('pdfunite GAUSS-GAMMA* fwhm-GAMMA.pdf',shell=True)
subprocess.call('pdfunite GAUSS-M* fwhm-M.pdf',shell=True)
subprocess.call('rm GAUSS-*',shell=True)


