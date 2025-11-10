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
mub50data=np.loadtxt('./mub50.dat')
mub100data=np.loadtxt('./mub100.dat')
mub150data=np.loadtxt('./mub150.dat')
mub200data=np.loadtxt('./mub200.dat')
mub250data=np.loadtxt('./mub250.dat')
mub300data=np.loadtxt('./mub300.dat')
mub350data=np.loadtxt('./mub350.dat')
mub400data=np.loadtxt('./mub400.dat')
mub450data=np.loadtxt('./mub450.dat')
mub500data=np.loadtxt('./mub500.dat')
mub550data=np.loadtxt('./mub550.dat')
mub600data=np.loadtxt('./mub600.dat')
mub650data=np.loadtxt('./mub650.dat')
mub700data=np.loadtxt('./mub700.dat')
beta=0.4022
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(np.exp(mub0data[:,0]),mub0data[:,1],'-',c='#F4DDB6',linewidth=2.5,alpha=1.,label=r'$\mu_B=0$')
# ax1.plot(mub50data[:,0],mub50data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu=320\,\mathrm{MeV}$')
# ax1.plot(mub100data[:,0],mub100data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu_B=100\,\mathrm{MeV}$')
# ax1.plot(mub150data[:,0],mub150data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu=320\,\mathrm{MeV}$')
# ax1.plot(mub200data[:,0],mub200data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu_B=200\,\mathrm{MeV}$')
# ax1.plot(mub250data[:,0],mub250data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu=320\,\mathrm{MeV}$')
ax1.plot(np.exp(mub300data[:,0]),mub300data[:,1],'-',c='#F5CABC',linewidth=2.5,alpha=1.,label=r'$\mu_B=300\,\mathrm{MeV}$')
# ax1.plot(mub350data[:,0],mub350data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu=320\,\mathrm{MeV}$')
ax1.plot(np.exp(mub400data[:,0]),mub400data[:,1],'-',c='#F4B4BE',linewidth=2.5,alpha=1.,label=r'$\mu_B=400\,\mathrm{MeV}$')
# ax1.plot(mub450data[:,0],mub450data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu=320\,\mathrm{MeV}$')
ax1.plot(np.exp(mub500data[:,0]),mub500data[:,1],'-',c='#CEA2B5',linewidth=2.5,alpha=1.,label=r'$\mu_B=500\,\mathrm{MeV}$')
# ax1.plot(mub550data[:,0],mub550data[:,1],'-',linewidth=2.5,alpha=0.6,label=r'$\mu=320\,\mathrm{MeV}$')
ax1.plot(np.exp(mub600data[:,0]),mub600data[:,1],'-',c='#9C9CC0',linewidth=2.5,alpha=1.,label=r'$\mu_B=600\,\mathrm{MeV}$')
ax1.plot(np.exp(mub650data[:,0]),mub650data[:,1],'-',c='#6F8BB8',linewidth=2.5,alpha=1.,label=r'$\mu_B=650\,\mathrm{MeV}$')
ax1.plot(np.exp(mub700data[:,0]),mub700data[:,1],'-',c='#59649D',linewidth=2.5,alpha=1.,label=r'$\mu_B=700\,\mathrm{MeV}$')
ax1.plot([10**-12,10**4],[beta,beta],dashes=[2,1],color='gray',linewidth=1.,alpha=0.6)
ax1.fill_between([10**-12,10**4],np.array([beta,beta])*0.99,np.array([beta,beta])*1.01,color='gray',edgecolor='none',alpha=0.3)
#ax1.text(221,1.2,r'$\mu=290\,\mathrm{MeV}$',fontsize=10)
ax1.set_xscale('log')
ax1.axis([10**-5,0.1,0.33,0.42])

ax1.set_xlabel(r'$t$', fontsize=14, color='black')
ax1.set_ylabel(r'$\mathrm{Leading\,\,order\,\,Fit}\,\,\beta$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize=9.5,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,scatterpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.14, left=0.145, right=0.95, hspace=0.35,wspace=0.35)

fig.savefig("beta_LO.pdf")
