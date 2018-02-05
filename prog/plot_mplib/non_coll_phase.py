#!/usr/bin/env python

import os
folder = "zero_field_classical"
import matplotlib.pyplot as plt
import matplotlib.figure as figure
from matplotlib.patches import FancyArrow
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle


## The tex style commands
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

# plotting the configuration

w,h = figure.figaspect(0.5)
fig = plt.figure(figsize=(w,h))
ax = fig.add_axes([0,0,1,1])

# first phase diagram

ffs = 34

arr1 = FancyArrow(0.25,0.75,0.9-0.25,0.0,
        ec='cornflowerblue',
        fc='cornflowerblue',
        width=0.025,
        head_width=0.08,
        head_length=0.05,
        zorder=1
        )
circ1 = Ellipse((0.25,0.75),0.04,0.08,
        ec='chocolate',
        fc='chocolate',
        zorder=2
        )

ax.text(0.9,0.625,r'$\delta J/J$',
        fontsize=30)

ax.text(0.245,0.625,r'$0$',
        fontsize=ffs)

ax.text(-0.05,0.75-0.025/2,r'(a) classical',
        fontsize=ffs)

ax.text(0.2125,0.85,r'LRO',
        color='chocolate',
        fontsize=ffs)

ax.text(0.5,0.85,r'spin glass',
        color='cornflowerblue',
        fontsize=ffs)

ax.add_artist(arr1)
ax.add_artist(circ1)

# second phase diagram

rect2 = Rectangle((0.25,0.25),
        (0.9-0.25)/2,0.025,
        ec='cornflowerblue',
        fc='cornflowerblue',
        zorder=1
        )

rect3 = Rectangle(((0.9-0.25)/2+0.25,0.225),
        0.001,0.075,
        ec='black',
        fc='black',
        zorder=2
        )

arr2 = FancyArrow((0.9-0.25)/2+0.25+0.001,
        0.25+0.025/2,0.9-((0.9-0.25)/2+0.25+0.001),0.0,
        ec='seagreen',
        fc='seagreen',
        width=0.025,
        head_width=0.08,
        head_length=0.05,
        zorder=1
        )

circ2 = Ellipse((0.25,0.2625),0.04,0.08,
        ec='chocolate',
        fc='chocolate',
        zorder=2
        )

ax.text(0.525,0.125,
        r'$\left(\delta J/J\right)_{c}$',
        fontsize=ffs)

ax.text(0.9,0.125,r'$\delta J/J$',
        fontsize=ffs)

ax.text(0.245,0.125,r'$0$',
        fontsize=ffs)

ax.text(-0.05,0.25,r'(b) quantum',
        fontsize=ffs)

ax.text(0.2125,0.35,r'LRO',
        color='chocolate',
        fontsize=ffs)

ax.text(0.34,0.35,r'spin glass',
        color='cornflowerblue',
        fontsize=ffs)

ax.text(0.615,0.35,r'random singlet',
        color='seagreen',
        fontsize=ffs)

ax.add_artist(rect2)
ax.add_artist(rect3)
ax.add_artist(circ2)
ax.add_artist(arr2)

ax.axis('off')
fig.savefig("fig_pd.pdf"
        ,bbox_inches='tight'
        )
plt.close('all')

