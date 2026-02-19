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
#T=np.loadtxt('Tem1/buffer/TMeV.dat')
LOdata=np.loadtxt('./LObeta.dat')
NLOdata=np.loadtxt('./NLObeta.dat')
LOdatamub550=np.loadtxt('./LObetamub550.dat')
NLOdatamub550=np.loadtxt('./NLObetamub550.dat')

beta=0.4005422694019118605889736521756
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(np.exp(LOdata[:,0]),LOdata[:,1],'-',c='#313772',linewidth=2.5,alpha=1.,label=r'$\mathrm{LO}\,\,\mu_B=0$')
ax1.plot(NLOdata[:,0],NLOdata[:,1],'-',c='#b7282e',linewidth=2.5,alpha=1.,label=r'$\mathrm{NLO}\,\,\mu_B=0$')
ax1.plot(np.exp(LOdatamub550[:,0]),LOdatamub550[:,1],dashes=[2,1],c='#478ecc',linewidth=2.5,alpha=1.,label=r'$\mathrm{LO}\,\,\mu_B=550\,\mathrm{MeV}$')
ax1.plot(NLOdatamub550[:,0],NLOdatamub550[:,1],dashes=[2,1],c='#dc917b',linewidth=2.5,alpha=1.,label=r'$\mathrm{NLO}\,\,\mu_B=550\,\mathrm{MeV}$')

ax1.plot([10**-12,10**4],[beta,beta],dashes=[2,1],color='gray',linewidth=1.,alpha=0.6)
ax1.fill_between([10**-12,10**4],np.array([beta,beta])*0.99,np.array([beta,beta])*1.01,color='gray',edgecolor='none',alpha=0.3)
ax1.fill_betweenx([0,1],[4*10**-6,4*10**-6],[0.01/700,0.01/700],color='red',edgecolor='none',alpha=0.3,zorder=10)
ax1.set_xscale('log')
ax1.text(2*10**-5,0.407,r'$\mathrm{LPA}$',fontsize=15)
ax1.axis([6*10**-6,0.03,0.37,0.44])

ax1.set_xlabel(r'$-t$', fontsize=14, color='black')
ax1.set_ylabel(r'$\mathrm{LO\,\,&\,\,NLO\,\,Fit}\,\,\beta$', fontsize=14, color='black')

ax1.legend(loc=(0.1,0.01),fontsize=9.5,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.145, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("beta_NLO_LPA.pdf")
