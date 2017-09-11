#!/usr/bin/env python

### This is the main program for mean field calculation of ###
### the disordered triangular lattice heisenberg model     ###

import os
os.chdir("../out/")
import sys
import subprocess
import re
import glob
import numpy as np
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

## Unlocking all the cpu ##

os.putenv("MKL_DYNAMIC","FALSE")

## The tex style commands
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')

# the font styleset
from matplotlib import rcParams
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.family'] = 'serif'

# the finite temperature fermi function
def fermi(E,T):
    
    return 1.0/(np.exp(E/T)+1.0)

def invfermi(E,T):
    
    return 1.0/(np.exp(-E/T)+1.0)

def susvc(SUSC,MTRX,N):

    for j in range(0,N):

        SUSC += (np.einsum('ijk,ijk',MTRX[i],
            np.einsum('aij,ij->aij',MTRX[j].conj(),DNM)))

## picturing the output data

# file name pattern

fpat = re.compile('FNL_L_([^/]*)_DLT_([^/]*)_ALP_([^/]*).npz')

if not os.path.isdir("plot"):
    os.mkdir("plot")

for fname in glob.iglob('*.npz'):

    match = fpat.match(fname)
    L = int(match.group(1))
    DELTA = float(match.group(2))
    ALPHA = float(match.group(3))

    FNDATA = np.load(fname)

    ITERDISD = FNDATA['ITERDISD'][()]
    BOOTNUM = FNDATA['BOOTNUM'][()]
    SIMSZE = FNDATA['SIMSZE'][()]
    T = FNDATA['T'][()]
    ensys = FNDATA['ensys']
    vsys = FNDATA['vsys']
    bond = FNDATA['bond']
    spin = FNDATA['spin']

    # calculating the onsite susceptibility

    for IDISD in range(0,ITERDISD):

        for BTNUM in range(0,BOOTNUM):

            N = L**2
            MDIM = 2*N

            SUSC = np.zeros(N,dtype=np.complex)
            MTRX = np.zeros((N,3,MDIM,MDIM),dtype=np.complex)
            DNM = (np.absolute(ensys[IDISD,BTNUM])[:, None] +
                    np.absolute(ensys[IDISD,BTNUM])[None, :] +
                    0.0001)**(-1.0)

            for i in range(0,N):

                evec = np.array([vsys[IDISD,BTNUM,i,:]*fermi(ensys[IDISD,BTNUM],T),
                    vsys[IDISD,BTNUM,i+N,:]*fermi(ensys[IDISD,BTNUM],T)]).conj()
                hvec = np.array([vsys[IDISD,BTNUM,i,:]*invfermi(ensys[IDISD,BTNUM],T),
                    vsys[IDISD,BTNUM,i+N,:]*invfermi(ensys[IDISD,BTNUM],T)])
                
                MTRX[i,0] = 0.5*DNM*(np.einsum('i,j->ij',
                    evec[0],hvec[1])+np.einsum('i,j->ij',
                    evec[1],hvec[0]))

                MTRX[i,1] = 0.5j*DNM*(-np.einsum('i,j->ij',
                    evec[0],hvec[1])+np.einsum('i,j->ij',
                    evec[1],hvec[0]))

                MTRX[i,2] = 0.5*DNM*(np.einsum('i,j->ij',
                    evec[0],hvec[0])-np.einsum('i,j->ij',
                        evec[1],hvec[1]))
            
            #Parallel(n_jobs=4)(delayed(susvc)(SUSC[i],MTRX,N) for i in range(0,N))

            for i in range(0,N):

                for j in range(0,N):

                    SUSC[i] += np.tensordot(MTRX[i],MTRX[j].conj(),axes=((0,1,2),(0,1,2)))
            """
            for i in range(0,N):

                for j in range(0,N):

                    SUSC[i] += (np.einsum('ijk,ijk',MTRX[i],
                            np.einsum('aij,ij->aij',MTRX[j].conj(),DNM)))

            """
            fig, ax = plt.subplots()

            sfc = SUSC.reshape((L,L)).real
 
            MAXCLIP = 10.0
            MINCLIP = MAXCLIP/(10.0**10)
           
            cax = ax.imshow(np.absolute(sfc),
                    norm=LogNorm(vmin=MINCLIP,vmax=MAXCLIP,clip=True),
                    interpolation='nearest',
                    extent=(0.0,1.0,0.0,1.0),
                    cmap=cm.inferno,aspect='auto')
            cbar = fig.colorbar(cax,shrink=0.5,format='%.0e')

            EN = np.sum(fermi(ensys[IDISD,BTNUM],T)*ensys[IDISD,BTNUM])/N

            plt.suptitle(r"Local Susceptibility, $\chi_{ij}'(\omega = 0)$", 
                    x=0.5, y=0.99, fontsize=16)
            plt.title(r"$\Delta$ = "+str(DELTA)+" , "
                    +r"$\alpha$ = "+str(ALPHA)+" , "
                    +"|E|/$L^{2}$ = "+str("%.3f" % EN)+" , "
                    +"Sample. = "+str(IDISD)+" , "
                    +"Init. config. = "+str(BTNUM),x=0.485,fontsize=12)
            smax = sfc.max()/N
            plt.figtext(.8,.85,
                    r"$\chi_{ij}^{'\mathrm{max}}/L^{2}$ = "+str("%.2e" % smax),
                    fontsize=12)
            fig.tight_layout(pad=1.6,h_pad=1.6,w_pad=1.6)
            fig.savefig("plot/susc_"+
                    str("%04d" % L)+
                    str("%.3f" % DELTA)+
                    str("%.3f" % ALPHA)+
                    str("%04d" % IDISD)+
                    str("%04d" % BTNUM)+
                    ".pdf"
                    )

            plt.close('all')


    subprocess.call('pdfunite plot/susc_'+
            str("%04d" % L)+
            str("%.3f" % DELTA)+
            str("%.3f" % ALPHA)+
            '* '+
            'plot/LSSC'+
            '_L_'+str(L)+
            '_DLT_'+str(DELTA)+
            '_ALP_'+str(ALPHA)+
            '.pdf',shell=True)

    # removing split files
    subprocess.call('rm plot/susc_*',shell=True)
