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
beta0p36=np.loadtxt('./beta0p36_fit.dat')
beta0p37=np.loadtxt('./beta0p37_fit.dat')
beta0p38=np.loadtxt('./beta0p38_fit.dat')
nu0p01=np.loadtxt('./nu0p01_fit.dat')
nu0p02=np.loadtxt('./nu0p02_fit.dat')
nu0p03=np.loadtxt('./nu0p03_fit.dat')
pb=np.loadtxt('./pb.dat')

pb2nd=np.loadtxt('./pb2nd.dat')
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(pb[:,0],nu0p01,'-',c='#206FB6',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;correlation\;length}=\pm 1\%\nu$')
ax1.plot(pb[:,0],nu0p02,'-',c='#6BADD7',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;correlation\;length}=\pm 2\%\nu$')
ax1.plot(pb[:,0],nu0p03,'-',c='#C5DAEE',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;correlation\;length}=\pm 3\%\nu$')
ax1.plot(pb[:,0],beta0p36,'-',c='#FDDFD0',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;order\;parameter}=0.36$')
ax1.plot(pb[:,0],beta0p37,'-',c='#FC9171',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;order\;parameter}=0.37$')
ax1.plot(pb[:,0],beta0p38,'-',c='#EE3B2A',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;order\;parameter}=0.38$')
# ax1.fill_between(pb[:,0],pb[:,1],beta0p36,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],beta0p37,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],beta0p38,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],nu0p01,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],nu0p02,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],nu0p03,color='gray',edgecolor='none',alpha=0.3)
ax1.fill_betweenx([0,200],[600,600],[1000,1000],color='gray',edgecolor='none',alpha=0.3,zorder=10)
ax1.plot(pb2nd[:,0],pb2nd[:,1]-0.5,'-',c='k',linewidth=1.5,alpha=1.,label=r'$\mathrm{2nd-Phase\;boundary}$')
#ax1.set_xscale('log')
ax1.text(300,125,r'$\mathrm{LPA}$',fontsize=15)
ax1.axis([0,700,50,135])
ax1.set_xlabel(r'$\mu_B\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$T\,[\mathrm{MeV}]$', fontsize=14, color='black')

ax1.legend(loc=(0.01,0.01),fontsize=7.9,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.145, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("beta_phase.pdf")
