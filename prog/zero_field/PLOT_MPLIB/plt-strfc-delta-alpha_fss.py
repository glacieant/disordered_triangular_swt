# A matplotlib plotting script to plot a simple x-y plot with yerrorbars
import os
os.chdir("../")
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import container
from matplotlib import lines
from matplotlib.legend_handler import HandlerErrorbar

# The tex style commands
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

############## the dataset range ################

CLR=['red','black','blue','brown','grey','sandybrown','pink','cyan']

IPX = 0

##### Quantum data #######

#loading data
l, delta, alpx, strfc, strfc_var, sample_size  = np.loadtxt('DATA/LISTDATA.dat', usecols=(0,1,26,24,25,33), unpack=True)

strfccl, strfccl_var  = np.loadtxt('DATA/LISTDATA.dat', usecols=(31,32), unpack=True)

##### Getting the standard error #####

strfc_var=(strfc_var/sample_size)**(0.5)

ALPX=np.unique(alpx)
NSET=len(ALPX)
L=np.unique(l)
LSET=len(L)
DELTA=np.unique(delta)
DSET=len(DELTA)

#Finite size scaling

eta=1.0

#Here reshaping the x axis such that x_i is replaced with 1.0/x_i
l=l**(-eta)

###### looping through the dataset ##############

#Extracting data with yerrorbar every N'th row starting with the M'th row

for d in range(0,DSET):
 
    #Making a plot with legends
    fig, ax = plt.subplots()
  
    for i in range(0,NSET):

        M=i+d*NSET
        N=NSET*DSET
        xp=l[M::N]
        yp=strfc[M::N]
        yerrorbarp=strfc_var[M::N]
        #Plotting with predefined labels
        line=ax.errorbar(xp, yp, yerr=yerrorbarp, xerr=None,
                        ecolor=CLR[i], elinewidth=2, marker='.', ms=10)

        #Setting the line width and style
        line[0].set_linestyle('solid')
        line[0].set_linewidth(2)
        line[0].set_color(CLR[i])

    for i in range(0,NSET):

        M=i+d*NSET
        N=NSET*DSET
        xpc=l[M::N]
        ypc=strfccl[M::N]
        yerrorbarpc=strfccl_var[M::N]
        #Plotting with predefined labels
        line=ax.errorbar(xpc, ypc, yerr=yerrorbarpc, xerr=None,
                ecolor=CLR[i], elinewidth=2, marker='.', ms=10)

        #Setting the line width and style
        line[0].set_linestyle('dashed')
        line[0].set_linewidth(2)
        line[0].set_color(CLR[i])


    #Setting up the legend and its position
    patch = []
    label = []
    for i in range(0,NSET):
        draw = mpatches.Patch(color=CLR[i])
        patch.append(draw)
        char = r'$\alpha $='+ str("%.3f" % ALPX[i])
        label.append(char)

    legend=ax.legend(patch,label,
            loc='best',shadow=True,
            borderpad=0.5, labelspacing=0.5,
            handlelength=1.0,handletextpad=0.0,
            columnspacing=0.5,ncol=2)
    #The frame is matplotlib.patches.Rectangel instance surrounding the legend
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    # Setting the frame fontzise and line width
    for i in legend.get_texts():
        i.set_fontsize(6)
    for i in legend.get_patches():
        i.set_width(10.0)
        i.set_height(5.0)
    cpatch=plt.Line2D((0,1),(0,0),color='k',ls='dashed',lw=2)
    clabel='Classical'
    qpatch=plt.Line2D((0,1),(0,0),color='k',ls='solid',lw=2)
    qlabel='Quantum'
    fig.legend([cpatch,qpatch],[clabel,qlabel],ncol=2,
            loc='lower center',fontsize=12)
    plt.xlabel(('$L'+'^{-'+str(eta)+'}$'),fontsize=18, labelpad=10)
    
    plt.suptitle(r'Structure factor $\chi(Q)/L^2$ at $(4\pi/3,0)$',fontsize=15, y=1.15)
    plt.title(r'$\Delta=$'+str("%.3f" % DELTA[d]),fontsize=12)

    #Setting up the tick lines and label formatting
    xmin=0.0
    xmax=np.amax(l)*1.1
    ymin=0.0
    ymax=0.6
    tickfreq=4
    xtickspace=(xmax-xmin)/tickfreq
    ytickspace=(ymax-ymin)/(tickfreq*4)
    plt.xticks(np.arange(xmin,xmax+xtickspace, xtickspace))
    plt.yticks(np.arange(ymin,ymax+ytickspace, ytickspace))
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.xaxis.set_tick_params(width=2,length=4,colors='r',right='off',labelsize=18)
    ax.yaxis.set_tick_params(width=2,length=4,colors='r',top='off',labelsize=10)
    #Setting up figure lay out
    fig.tight_layout(pad=1.0,h_pad=1.0,w_pad=1.0)

    #plt.show()
    fig.savefig("DATA/RAW/plt_strfc_delta="+str("%.3f"%DELTA[d])+str("%d"%IPX)+".pdf")
    IPX += 1
    plt.close('all')

if not os.path.exists("DATA/PLOT"):
    os.makedirs("DATA/PLOT")

# merging files and cleaning directories
subprocess.call('pdfunite DATA/RAW/plt_strfc_* DATA/PLOT/strfc.pdf',shell=True)
subprocess.call('rm DATA/RAW/plt_strfc_*',shell=True)
