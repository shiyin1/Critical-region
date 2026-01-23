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
mub0data=np.loadtxt('./mub0.dat')
mub300data=np.loadtxt('./mub300.dat')
mub400data=np.loadtxt('./mub400.dat')
mub500data=np.loadtxt('./mub500.dat')
mub550data=np.loadtxt('./mub550.dat')
mub575data=np.loadtxt('./mub575.dat')
mub600data=np.loadtxt('./mub600.dat')
delta=5.
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(mub0data[:,0],mub0data[:,1],'-',c='#F4DDB6',linewidth=2.5,alpha=1.,label=r'$\mu_B=0$')
ax1.plot(mub300data[:,0],mub300data[:,1],'-',c='#F5CABC',linewidth=2.5,alpha=1.,label=r'$\mu_B=300\,\mathrm{MeV}$')
ax1.plot(mub400data[:,0],mub400data[:,1],'-',c='#F4B4BE',linewidth=2.5,alpha=1.,label=r'$\mu_B=400\,\mathrm{MeV}$')
ax1.plot(mub500data[:,0],mub500data[:,1],'-',c='#CEA2B5',linewidth=2.5,alpha=1.,label=r'$\mu_B=500\,\mathrm{MeV}$')
ax1.plot(mub550data[:,0],mub550data[:,1],'-',c='#9C9CC0',linewidth=2.5,alpha=1.,label=r'$\mu_B=550\,\mathrm{MeV}$')
ax1.plot(mub575data[:,0],mub575data[:,1],'-',c='#6F8BB8',linewidth=2.5,alpha=1.,label=r'$\mu_B=575\,\mathrm{MeV}$')
ax1.plot(mub600data[:,0],mub600data[:,1],'-',c='#59649D',linewidth=2.5,alpha=1.,label=r'$\mu_B=600\,\mathrm{MeV}$')
ax1.plot([10**-4,10**4],[delta,delta],dashes=[2,1],color='gray',linewidth=1.,alpha=0.6)
ax1.fill_between([10**-4,10**4],np.array([delta,delta])*0.99,np.array([delta,delta])*1.01,color='gray',edgecolor='none',alpha=0.3)
ax1.fill_betweenx([0,10],[10**-6,10**-6],[0.002,0.002],color='red',edgecolor='none',alpha=0.3,zorder=10)
ax1.text(3*10**-3,5.2,r'$\mathrm{LPA}$',fontsize=15)
ax1.set_xscale('log')
ax1.axis([3*10**-4,500,3.0,6.1])

ax1.set_xlabel(r'$m_\pi\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\mathrm{Leading\,\,order\,\,Fit}\,\,\delta$', fontsize=14, color='black')

ax1.legend(loc=(0.12,-0.01),fontsize=9.5,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.145, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("delta_LO.pdf")
