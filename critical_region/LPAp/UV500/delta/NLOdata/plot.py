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
LOdata=np.loadtxt('./LOdelta.dat')
NLOdata=np.loadtxt('./NLOdelta.dat')
LOdatamub500=np.loadtxt('./LOmub500.dat')
NLOdatamub500=np.loadtxt('./NLOdeltamub500.dat')

#NLOdata[0:23,1]=NLOdata[24,1]
NLOdata[25:43,1]=NLOdata[44,1]
NLOdata[45:48,1]=NLOdata[49,1]

NLOdatamub500[25:43,1]=NLOdatamub500[44,1]
NLOdatamub500[45:48,1]=NLOdatamub500[49,1]
delta=4.79
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(LOdata[:,0],LOdata[:,1],'-',c='#313772',linewidth=2.5,alpha=1.,label=r'$\mathrm{LO}\,\,\mu_B=0$')
ax1.plot(NLOdata[:,0],NLOdata[:,1],'-',c='#b7282e',linewidth=2.5,alpha=1.,label=r'$\mathrm{NLO}\,\,\mu_B=0$')
ax1.plot(LOdatamub500[:,0],LOdatamub500[:,1],dashes=[2,1],c='#478ecc',linewidth=2.5,alpha=1.,label=r'$\mathrm{LO}\,\,\mu_B=500\,\mathrm{MeV}$')
ax1.plot(NLOdatamub500[:,0],NLOdatamub500[:,1],dashes=[2,1],c='#dc917b',linewidth=2.5,alpha=1.,label=r'$\mathrm{NLO}\,\,\mu_B=500\,\mathrm{MeV}$')

ax1.plot([10**-4,10**4],[delta,delta],dashes=[2,1],color='gray',linewidth=1.,alpha=0.6)
ax1.fill_between([10**-4,10**4],np.array([delta,delta])*0.99,np.array([delta,delta])*1.01,color='gray',edgecolor='none',alpha=0.3)
ax1.fill_betweenx([0,10],[10**-6,10**-6],[0.005,0.005],color='red',edgecolor='none',alpha=0.3,zorder=10)
ax1.text(10**-2,5.3,r'$\mathrm{LPA}^\prime$',fontsize=15)
ax1.set_xscale('log')
ax1.axis([2*10**-4,300,2.,8.5])

ax1.set_xlabel(r'$m_\pi\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\mathrm{LO\,\,&\,\,NLO\,\,Fit}\,\,\delta$', fontsize=14, color='black')

ax1.legend(loc=(0.21, 0.05),fontsize=9.5,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.12, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("delta_NLO_LPAp.pdf")
