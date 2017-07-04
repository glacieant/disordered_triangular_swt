#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# the font styleset
from matplotlib import rcParams
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.family'] = 'serif'

fig, ax = plt.subplots()

# LRO data
strfc = np.array([[0.105,0.0],[0.09,0.2],[0.08,0.3],
    [0.055,0.5],[0.05,0.55],[0.035,0.65],
    [0.03,0.7],[0.02,0.75],[0.015,0.8],[0.01,0.85],[0.005,0.9],[0.0,0.99]])

# LRO plot
ax.plot(strfc[:,0],strfc[:,1],'bs')

# defining parabolic fitting function
def sfc_b(x,a,b,c,d,e):
    return a + b*x + c*(x**2) + d*(x**3) + e*(x**4)

# LRO fitting and plotting
N = len(strfc[:,0])
# more weight for the end and mid points
sigma = np.ones(N)
sigma[[0,N/2,-1]] = 0.01
popt, pcov = curve_fit(sfc_b,strfc[:,0],strfc[:,1],
        sigma=sigma)
ax.plot(strfc[:,0],sfc_b(strfc[:,0], *popt),'b-')
ax.text(0.025,0.3,'Long range AFM order',
        fontsize=12,
        bbox={'facecolor':'white',
            'alpha':1.0,
            'edgecolor':'none',
            'pad':5})
ax.fill_between(np.append(strfc[:,0],0.0),0,
        np.append(sfc_b(strfc[:,0], *popt),0.99),facecolor='skyblue')

# local moment data

moment_1 = np.array([[0.105,0.0],[0.095,0.2],[0.09,0.3],
    [0.08,0.5],[0.085,0.55],[0.1,0.6]])

moment_2 = np.array([[0.14,0.0],[0.125,0.2],[0.12,0.3],
    [0.11,0.5],[0.105,0.55],[0.1,0.6]])

# local moment plot

ax.plot(moment_1[:,0],moment_1[:,1],'ro')
ax.plot(moment_2[:,0],moment_2[:,1],'ro')
ax.text(0.0935,0.3,'Non-magnetic \n phase',
        fontsize=12,
        bbox={'facecolor':'white',
            'alpha':1.0,
            'edgecolor':'none',
            'pad':5})
ax.fill_betweenx(moment_1[:,1],
        moment_1[:,0],
        moment_2[:,0],facecolor='indianred')

ax.text(0.095,0.8,'Disordered paramagnet',
        fontsize=12,
        bbox={'facecolor':'white',
            'alpha':1.0,
            'edgecolor':'none',
            'pad':5})

plt.xlabel(r'$\alpha$',fontsize=18)
plt.ylabel(r'$\Delta$',fontsize=18)
plt.title(
'Phase diagram of disordered $J_1-J_2$ triangular lattice'+'\n Heisenberg model',
y=1.05,fontsize=18)

x_major_ticks = np.arange(0.0,0.16,0.02)
x_minor_ticks = np.arange(0.0,0.16,0.005)
y_major_ticks = np.arange(0.0,0.99,0.2)
y_minor_ticks = np.arange(0.0,0.99,0.05)

ax.set_xticks(x_major_ticks)
ax.set_xticks(x_minor_ticks,minor=True)
ax.set_yticks(y_major_ticks)
ax.set_yticks(y_minor_ticks,minor=True)

ax.set_xlim([0.0,0.14])
ax.set_ylim([0.0,0.99])

ax.grid(which='minor',alpha=0.4)
ax.grid(which='major',alpha=0.6)

fig.tight_layout(pad=1.6,h_pad=1.6,w_pad=1.6)
fig.savefig("phase_diagram.pdf")
