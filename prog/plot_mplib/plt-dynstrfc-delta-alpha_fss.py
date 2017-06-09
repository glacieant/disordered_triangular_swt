# The import commands
import os
os.chdir("../high_field/out/")
from math import log
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# The tex style commands
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

L, DISD, ALPX, OMEGA, KX, KY, SFC, SFC_VAR  = np.loadtxt('GRX', usecols=(0,1,7,2,3,4,5,6), unpack=True)

#Setting up multiple dataset

DELTA=np.unique(DISD)
DSET=len(DELTA)
ALPHA=np.unique(ALPX)
ASET=len(ALPHA)
WFREQ=len(np.unique(OMEGA))
XARR = KX+1j*KY
XFREQ=len(np.unique(XARR))

#Positioning the Gamma properly

kyf=len(np.unique(KY))-1
kxf=XFREQ-kyf
gmpos=(1.0*kyf)/(kxf+kyf)

#The data array
DRANGE=WFREQ*XFREQ;

# Setting up the clipping of the density plot 
MAXCLIP=np.amax(SFC)
MINCLIP=MAXCLIP/(10**4)
mxc=np.log10(MAXCLIP)
mnc=np.log10(MINCLIP)
CLIPTICK=np.logspace(mnc,mxc,num=4,endpoint=True,base=10.0)

for i in range(0,DSET):

    for j in range(0,ASET):

        stride=i*(DRANGE*ASET)+j*DRANGE

        #Making plot with vertical colorbar
        fig, ax = plt.subplots()

        #Formatting the data

        xp=KX[stride:stride+DRANGE]
        xpp=KY[stride:stride+DRANGE]
        yp=OMEGA[stride:stride+DRANGE]
        zp=SFC[stride:stride+DRANGE]+10**-13
        
        #Creating a rectangular grid

        #The width of the pixel
        zmat=zp.reshape((WFREQ,XFREQ))

        #Imshow plotting the data
        #cax=ax.imshow(np.flipud(zmat), norm=Normalize(vmin=MINCLIP, vmax=MAXCLIP, clip=False),
        cax=ax.imshow(np.flipud(zmat), norm=LogNorm(vmin=MINCLIP, vmax=MAXCLIP, clip=False),
                interpolation='nearest',
                extent=(0.0, 1.0, np.amin(yp), np.amax(yp)),
                cmap=cm.afmhot ,aspect='auto')
        cbar=fig.colorbar(cax,shrink=0.5,
        #        ticks=CLIPTICK,
                format='%.0e')
        xtickpos=[0.0,gmpos,1.0]
        xticklabel=[r'M',r'$\Gamma$',r'K']
        plt.xticks(xtickpos,xticklabel)
        plt.xlabel(r'\textbf{q}', fontsize=18)
        plt.ylabel(r'$\omega$/J', fontsize=18)
        plt.title(r'S(\textbf{q},$\omega$) at $\Delta$='+str("%.3f" % DELTA[i])+r', $\alpha$='+str("%.3f" % ALPHA[j]), fontsize=20,y=1.09)
        fig.tight_layout(pad=1.0,h_pad=1.0,w_pad=1.0)
        fig.savefig('GRID-w-k-'+str("%.3f" % DELTA[i])+'-'+str("%.3f" % ALPHA[j])+'.pdf',format='pdf')
        plt.close('all')
