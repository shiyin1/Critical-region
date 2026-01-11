#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.style.use('classic')

# Data for plotting
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = X**2 + Y**2
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111, projection='3d')
ax1.plot_surface(X, Y, Z)
#ax1.text(221,1.2,r'$\mu=290\,\mathrm{MeV}$',fontsize=10)
#ax1.set_xscale('log')
#ax1.axis([6*10**-6,0.03,0.3,0.45])

ax1.set_xlabel(r'$t$', fontsize=14, color='black')
ax1.set_ylabel(r'$\mathrm{Leading\,\,order\,\,Fit}\,\,\beta$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize=9.5,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.145, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("beta_LO.pdf")
