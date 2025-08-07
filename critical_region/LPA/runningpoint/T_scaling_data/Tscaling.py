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

mpl.style.use('classic')
# Data for plotting
mub0  =np.loadtxt('./mub0.dat')
mub50 =np.loadtxt('./mub50.dat')
mub100=np.loadtxt('./mub100.dat')
mub150=np.loadtxt('./mub150.dat')
mub200=np.loadtxt('./mub200.dat')
mub250=np.loadtxt('./mub250.dat')
mub300=np.loadtxt('./mub300.dat')
mub350=np.loadtxt('./mub350.dat')
mub400=np.loadtxt('./mub400.dat')
mub450=np.loadtxt('./mub450.dat')
mub500=np.loadtxt('./mub500.dat')
mub550=np.loadtxt('./mub550.dat')
mub600=np.loadtxt('./mub600.dat')
mub650=np.loadtxt('./mub650.dat')
mub700=np.loadtxt('./mub700.dat')

ps=np.arange(1, 2000, 10)
#print(mu)
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(mub0[:,0]  ,mub0[:,1],  color=[0.00,0.20,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=0$') 
ax1.plot(mub50[:,0] ,mub50[:,1], color=[0.05,0.25,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=50\,\mathrm{MeV} $') 
ax1.plot(mub100[:,0],mub100[:,1],color=[0.10,0.30,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=100\,\mathrm{MeV} $') 
ax1.plot(mub150[:,0],mub150[:,1],color=[0.15,0.35,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=150\,\mathrm{MeV} $') 
ax1.plot(mub200[:,0],mub200[:,1],color=[0.20,0.40,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=200\,\mathrm{MeV} $') 
ax1.plot(mub250[:,0],mub250[:,1],color=[0.25,0.45,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=250\,\mathrm{MeV} $') 
ax1.plot(mub300[:,0],mub300[:,1],color=[0.30,0.50,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=300\,\mathrm{MeV} $') 
ax1.plot(mub350[:,0],mub350[:,1],color=[0.35,0.55,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=350\,\mathrm{MeV} $') 
ax1.plot(mub400[:,0],mub400[:,1],color=[0.40,0.60,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=400\,\mathrm{MeV} $') 
ax1.plot(mub450[:,0],mub450[:,1],color=[0.45,0.65,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=450\,\mathrm{MeV} $') 
ax1.plot(mub500[:,0],mub500[:,1],color=[0.50,0.70,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=500\,\mathrm{MeV} $') 
ax1.plot(mub550[:,0],mub550[:,1],color=[0.55,0.75,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=550\,\mathrm{MeV} $') 
ax1.plot(mub600[:,0],mub600[:,1],color=[0.60,0.80,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=600\,\mathrm{MeV} $') 
ax1.plot(mub650[:,0],mub650[:,1],color=[0.65,0.85,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=650\,\mathrm{MeV} $') 
ax1.plot(mub700[:,0],mub700[:,1],color=[0.70,0.90,0.5],linewidth=2.,alpha=0.7,label=r'$ \mu_B=700\,\mathrm{MeV} $') 

ax1.plot([-13,10],[0.402,0.402],color='k',dashes=[0.5,0.5],linewidth=1.,alpha=0.7) 

ax1.set_ylabel(r'$\beta$', fontsize=13, color='black')
ax1.set_xlabel(r'$\mathrm{ln}\,t$', fontsize=13, color='black')
ax1.legend(loc=(0.01,0.01),fontsize='6.5',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1,scatterpoints=1)
ax1.axis([-11.5,-1.8,0.32,0.41])
ax1.text(-3.5, 0.403, r'0.402', fontsize=12, color='k')
for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.16, right=0.95, hspace=0.33,wspace=0.2)

fig.savefig("./Tscaling.pdf")


