#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl
#from scipy.interpolate import spline
from matplotlib import cm
from matplotlib import axes
from matplotlib.font_manager import FontProperties
import pylab as pl

mpl.style.use('classic')
# Data for plotting
msigma1=np.loadtxt('./msigma1.dat')
msigma2=np.loadtxt('./msigma2.dat')
msigma3=np.loadtxt('./msigma3.dat')
msigma4=np.loadtxt('./msigma4.dat')
mub1=np.loadtxt('./mub1.dat')
mub2=np.loadtxt('./mub2.dat')
mub3=np.loadtxt('./mub3.dat')
mub4=np.loadtxt('./mub4.dat')
T=np.loadtxt('./T.dat')
T=T-0.3

msigma2 = np.nan_to_num(msigma2, nan=100.0) 
msigma3 = np.nan_to_num(msigma3, nan=100.0) 
msigma4 = np.nan_to_num(msigma4, nan=100.0) 

beta0p36=np.loadtxt('./beta0p36.dat')
beta0p37=np.loadtxt('./beta0p37.dat')
beta0p38=np.loadtxt('./beta0p38.dat')
mubbeta=[0,50,100,150,200,250,300,350,400,450,500,550,575,600]

nu0p01=np.loadtxt('./nu0p01.dat')
nu0p02=np.loadtxt('./nu0p02.dat')
nu0p03=np.loadtxt('./nu0p03.dat')
mubnu=[0,50,100,150,200,250,300,350,400,450,500,550,575,600]

pb=np.loadtxt('./pb.dat')
pb2nd=np.loadtxt('./pb2nd.dat')
background=np.zeros([2101,3201])
# Create figure
fig=plt.figure(figsize=(4.57, 3))
#fig=plt.figure()
####################################################################################################
# Create figure
#fig=plt.figure()
ax2=fig.add_subplot(111)
im1=ax2.pcolormesh(mub1, T, 1./msigma1.T, cmap='Reds', vmin=0, vmax=1)
im2=ax2.pcolormesh(mub2, T, 1./msigma2.T, cmap='Reds', vmin=0, vmax=1)
im3=ax2.pcolormesh(mub3, T, 1./msigma3.T, cmap='Reds', vmin=0, vmax=1)
im4=ax2.pcolormesh(mub4, T, 1./msigma4.T, cmap='Reds', vmin=0, vmax=1)
vnorm = mpl.colors.Normalize(vmin=0, vmax=1)
plt.rcParams['font.size'] = 7
cbar=plt.colorbar(im1,fraction=0.031, pad=0.01,norm=vnorm)
cbar.set_label(r'$\xi\,\,[\mathrm{MeV}^{-1}]$', rotation=270,fontsize=10, labelpad=15)
# ax2.plot(mubnu,nu0p01,'-',c='#206FB6',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;correlation\;length}=\pm 1\%\nu$')
# ax2.plot(mubnu,nu0p02,'-',c='#6BADD7',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;correlation\;length}=\pm 2\%\nu$')
# ax2.plot(mubnu,nu0p03,'-',c='#C5DAEE',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;correlation\;length}=\pm 3\%\nu$')
# ax2.plot(mubbeta,beta0p36,'-',c='#FDDFD0',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;order\;parameter}=0.36$')
# ax2.plot(mubbeta,beta0p37,'-',c='#FC9171',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;order\;parameter}=0.37$')
# ax2.plot(mubbeta,beta0p38,'-',c='#EE3B2A',linewidth=.5,alpha=1.,label=r'$\mathrm{Slope\;of\;order\;parameter}=0.38$')
# ax2.plot(pb2nd[:,0],pb2nd[:,1]-0.5,'-',c='k',linewidth=1.5,alpha=1.,label=r'$\mathrm{2nd-Phase\;boundary}$')
# ax2.fill_betweenx([0,200],[600,600],[1000,1000],color='gray',edgecolor='none',alpha=0.3,zorder=10)

#im2=ax2.imshow(background-2, cmap=plt.get_cmap('binary'),interpolation='nearest',vmin=-10,vmax=10,zorder=1,alpha=0.8)#plt.cm.hot_r)
#ax2.invert_yaxis()
#And,=ax2.plot(mubdata,FOT2*10-400,color='b',linewidth=1.5,alpha=0.6,label=r'freezeout: Andronic et al.',zorder=3)
#star7,=ax2.plot(mubdata,FOT3*10-400,dashes=[4,1,2,1],color='r',linewidth=1.5,alpha=0.6,label=r'freezeout: STAR Fit I',zorder=5)
#star4,=ax2.plot(mubdata,FOT*10-400,dashes=[5,2],color='g',linewidth=1.5,alpha=0.6,label=r'freezeout: STAR Fit II',zorder=4)
#CEP,=ax2.plot(603*4,107*10-400,marker='o',color='r',linewidth=0,markersize=6,markeredgewidth=1.5,label=r'CEP',zorder=3)

plt.axis([0.,594.,50.,130.])
ax2.text(490, 120, r'$\mathrm{LPA}$',fontsize=15, color='k')

#plt.yticks([100,600,1100,1600,2100],[50,100,150,200,250])
#plt.xticks([0,400,800,1200,1600,2000,2400,2800,3200],[0,100,200,300,400,500,600,700,800])
ax2.set_xlabel(r'$\mu_B\,[\mathrm{MeV}]$', fontsize=12, color='black')
ax2.set_ylabel(r'$T\,[\mathrm{MeV}]$', fontsize=12, color='black')
#ax1.set_title(r'$R^B_{42}$',loc='right')
for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(8)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(8)

#ax2.legend(loc=0,fontsize=8,frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)

fig.subplots_adjust(top=0.90, bottom=0.15, left=0.12, right=0.89, hspace=0.1, wspace=0.1)

fig.savefig("beta_nu_phase_LPA_v2.pdf",dpi=300)
