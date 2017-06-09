# A matplotlib plotting script to plot a simple x-y plot with yerrorbars
import os
os.chdir("../ALPHA_DELTA_PHASE/swt-mft-j1-j2-alpha-delta-job_1/plot/")
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

##### Quantum data #######

#loading data
l, delta, alpx, sblcorr, sblcorr_var = np.loadtxt('LISTDATA.dat', usecols=(0,1,26,6,7), unpack=True)

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

##### Classical data #####

#loading data
sblcorrcl, sblcorrcl_var = np.loadtxt('LISTDATA.dat', usecols=(10,11), unpack=True)

###### looping through the dataset ##############

#Extracting data with yerrorbar every N'th row starting with the M'th row

for d in range(0,DSET):
 
    #Making a plot with legends
    fig, ax = plt.subplots()
  
    for i in range(0,NSET):

        M=i+d*NSET
	N=NSET*DSET
	xp=l[M::N]
        yp=sblcorr[M::N]
        yerrorbarp=sblcorr_var[M::N]
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
        ypc=sblcorrcl[M::N]
        yerrorbarpc=sblcorrcl_var[M::N]
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
            char = r'$\alpha$='+ str("%.3f" % ALPX[i])
            label.append(char)

    legend=ax.legend(patch,label,
                    loc='best',shadow=True,
                    borderpad=0.5, labelspacing=0.5,
                    handlelength=1.0,handletextpad=0.0,
                    columnspacing=0.5,ncol=2)
    #The frame is matplotlib.patches.Rectangel instance surrounding the legend
    frame =legend.get_frame()
    frame.set_facecolor('0.90')
    #Setting the frame fontzise and line width
    for i in legend.get_texts():
            i.set_fontsize(6)
    for i in legend.get_patches():
            i.set_width(10.0)
            i.set_height(5.0)
    cpatch=plt.Line2D((0,1),(0,0),color='k',ls='dashed',lw=2)
    clabel='Classical'
    qpatch=plt.Line2D((0,1),(0,0),color='k',ls='solid',lw=2)
    qlabel='Quantum'
    fig.legend([cpatch,qpatch],[clabel,qlabel],ncol=2,bbox_to_anchor=(0.5,0.9),
            loc='upper center',fontsize=12)
    plt.xlabel(('$L'+'^{-'+str(eta)+'}$'),fontsize=18, labelpad=10)
    plt.title(r'Sublattice correlation, $\Delta=$'+str("%.3f" % DELTA[d]),fontsize=15, y=1.15)
    #Setting up the tick lines and label formatting
    xmin=0.0
    xmax=np.amax(l)*1.1
    ymin=0.0
    ymax=np.amax(sblcorr[d*NSET::DSET])*1.05
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
    fig.savefig("sblcorr_delta="+str("%.3f" % DELTA[d])+".pdf")
    plt.close('all')

