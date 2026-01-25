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
delta0p01=np.loadtxt('./delta0p01.dat')
delta0p02=np.loadtxt('./delta0p02.dat')
delta0p03=np.loadtxt('./delta0p03.dat')
mub=[0,100,200,300,400,450,500,530]
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(mub,delta0p01,'-',c='#206FB6',linewidth=1.5,alpha=1.,label=r'$\mathrm{Slope\;of\;}m_\pi=\pm 1\%\delta$')
ax1.plot(mub,delta0p02,'-',c='#6BADD7',linewidth=1.5,alpha=1.,label=r'$\mathrm{Slope\;of\;}m_\pi=\pm 2\%\delta$')
ax1.plot(mub,delta0p03,'-',c='#C5DAEE',linewidth=1.5,alpha=1.,label=r'$\mathrm{Slope\;of\;}m_\pi=\pm 3\%\delta$')
# ax1.fill_between(pb[:,0],pb[:,1],beta0p36,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],beta0p37,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],beta0p38,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],nu0p01,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],nu0p02,color='gray',edgecolor='none',alpha=0.3)
# ax1.fill_between(pb[:,0],pb[:,1],nu0p03,color='gray',edgecolor='none',alpha=0.3)
#ax1.fill_betweenx([0,200],[600,600],[1000,1000],color='gray',edgecolor='none',alpha=0.3,zorder=10)
#ax1.plot(pb2nd[:,0],pb2nd[:,1]-0.5,'-',c='k',linewidth=1.5,alpha=1.,label=r'$\mathrm{2nd-Phase\;boundary}$')
ax1.set_yscale('log')
ax1.text(270,30,r'$\mathrm{LPA}^\prime$',fontsize=15)
ax1.axis([0,600,10**-2,100])
ax1.set_xlabel(r'$\mu_B\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$m_\pi\,[\mathrm{MeV}]$', fontsize=14, color='black')

ax1.legend(loc=(0.01,0.1),fontsize=10,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.15, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("delta_phase_LPAp.pdf")
