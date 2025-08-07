#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import FuncFormatter
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.cm as cm 
from matplotlib.font_manager import FontProperties
import pylab as pl
from scipy.ndimage import zoom

mpl.style.use('classic')
# Data for plotting
data  =np.loadtxt('./data2D.dat')
mubdata = [0,50,100,150,200,250,300,350,400,450,500,550,575,600,625,650,675,700]
Tdata= np.loadtxt('./T.dat')
data_interp = zoom(data, 5, order=3) 
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)

xnew, ynew = np.meshgrid(mubdata, Tdata)
#print(xnew.shape)
#print(ynew.shape)
#ax1.contour(xnew, ynew, data.T, levels=5, colors='black')
ax1.imshow(data_interp.T,cmap='plasma',origin='lower',aspect='auto',interpolation='gaussian', vmin=0.1, vmax=0.41)
#ax1.contour(data.T, levels=[0.39], colors='black')

ax1.set_ylabel(r'$\beta$', fontsize=13, color='black')
ax1.set_xlabel(r'$\mathrm{ln}\,t$', fontsize=13, color='black')
#ax1.legend(loc=(0.01,0.01),fontsize='6.5',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([0,5*20,0,5*8000])
#ax1.text(-3.5, 0.403, r'0.402', fontsize=12, color='k')
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.16, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./Tscaling_hot.pdf")