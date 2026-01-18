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
kappa_LPA=np.loadtxt('./kappa_LPA.dat')*197.33**2
k_LPA=np.loadtxt('./k_LPA.dat')
kappa_LPAp=np.loadtxt('./kappa_LPAp.dat')*197.33**2
k_LPAp=np.loadtxt('./k_LPAp.dat')
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(k_LPA,kappa_LPA/(118.7192459*k_LPA),'-',c='#8589B0',linewidth=2.5,alpha=1.,label=r'$LPA,\,\,\mu_B=0$',zorder=0)
ax1.plot(k_LPAp,kappa_LPAp/(110.0184745*k_LPAp),'-',c='#91B09F',linewidth=2.5,alpha=1.,label=r'$LPA^\prime,\,\,\mu_B=0$',zorder=2)
ax1.plot([5e-3,5e-3],[0,1],dashes=[1,1],c='r',linewidth=1.,alpha=0.5,zorder=1)
ax1.plot([1e-2,1e-2],[0,1],dashes=[1,1],c='m',linewidth=1.,alpha=0.5,zorder=3)
#ax1.fill_between([1e-4,3e-2],[0,0],[1,1],color='gray',edgecolor='none',alpha=0.5)
#ax1.text(221,1.2,r'$\mu=290\,\mathrm{MeV}$',fontsize=10)
ax1.set_xscale('log')
ax1.axis([10**-3,500,0,0.1])

ax1.set_xlabel(r'$k\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\frac{\bar\rho}{T_c\,k}$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize=9.5,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.17, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("kappa_k.pdf")
