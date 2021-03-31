#This code is for making a plot of mass loss
#Created by Malia Jenks

import numpy as np
import pylab
from math import sqrt

runname = '32node_77_cgs_spitzer'
datfile = runname + '/relax.dat'
time, mass, = np.loadtxt(datfile,skiprows=1,usecols=(0,1),unpack=True)

golden = (1.0 + sqrt(5.0)) / 2.0
 
figprops = dict(figsize=(8., 8./golden), dpi=128)
adjustprops = dict(left=0.1, bottom=0.12, right=0.97, top=0.93, wspace=0.3, hspace=0.3)
 
fig1 = pylab.figure(**figprops)
fig1.subplots_adjust(**adjustprops)
 
ax = fig1.add_subplot(111)

p1, = ax.plot(time, mass, linewidth=2.0)
 
ax.set_xlabel("Time (s)", fontsize = 20)
ax.set_ylabel("Mass (g)", fontsize = 20)
name = 'plots/massloss_' + runname
pngname = name + ".png"
fig1.savefig(pngname)
 
pylab.close
