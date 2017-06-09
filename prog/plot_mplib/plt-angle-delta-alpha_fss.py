# A matplotlib plotting script to plot a simple x-y plot with yerrorbars
import os
os.chdir("../ALPHA_DELTA_PHASE/wth_angle_corr_j2/DATA/")
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import container
from matplotlib import lines
from matplotlib.legend_handler import HandlerErrorbar

# The tex style commands
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

############## the dataset range ################

CLR=['black','royalblue','dimgray','gold','limegreen','red','deeppink']

#loading data
L, delta, alpx, data, data_corr = np.loadtxt('ANGLEDATA.dat', usecols=(0,1,2,4,6), unpack=True)

ALPX=np.unique(alpx)
NSET=len(ALPX)
LX=np.unique(L)
LSET=len(LX)
DELTA=np.unique(delta)
DSET=len(DELTA)

###### looping through the dataset ##############

#Extracting data with yerrorbar every N'th row starting with the M'th row

STRIDE=0

for l in range(0,LSET):

    SYS=(LX[l]**2)*3

    for d in range(0,DSET):

        for a in range(0,NSET):

            #Making a plot with legends

            fig, ax = plt.subplots()

            M=STRIDE
            N=M+SYS
            STRIDE=N
            xp=data[M:N]
            xpcorr=data_corr[M:N]
            #Plotting with predefined labels
            ax.hist(xp,bins=np.linspace(-1.0,1.0,200),alpha=0.5,color=CLR[0])
            ax.hist(xpcorr,bins=np.linspace(-1.0,1.0,200),alpha=0.5,color=CLR[3])

            patch = []
            label = []
            draw = mpatches.Patch(color=CLR[0])
            patch.append(draw)
            char = r'$\theta_i^C$ , classical angles'
            label.append(char)
            draw = mpatches.Patch(color=CLR[3])
            patch.append(draw)
            char = r'$\theta_i^Q$ , quantum corrected angles'
            label.append(char)

            legend=ax.legend(patch,label,
                    loc='best',shadow=True,
                    borderpad=0.5, labelspacing=0.5,
                    handlelength=1.0,handletextpad=0.5,
                    columnspacing=0.5,ncol=1)

            plt.xlim(-1.0,1.0)
            plt.title(r'Histogram of bond angle cosines, L='+str("%d"%LX[l])+',$\Delta=$'
                    +str("%.3f"%DELTA[d])+r',$\alpha=$'+str("%.3f"%ALPX[a]),fontsize=15, y=1.15)

            #Setting up figure lay out
            fig.tight_layout(pad=1.0,h_pad=1.0,w_pad=1.0)
            #plt.show()
            fig.savefig("hist_l_="+str("%d"%LX[l])+"_delta="+str("%.3f"%DELTA[d])+"_alpha="+str("%.3f"%ALPX[a])+".pdf")
            plt.close('all')


