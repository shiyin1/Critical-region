#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D 

mpl.style.use('classic')

# Data for plotting
pb=np.loadtxt('./pb.dat')
mub=np.loadtxt('./mub.dat')
beta=np.loadtxt('./betaregion.dat')

pb2=np.loadtxt('./pb2.dat')
mub2=np.loadtxt('./mub2.dat')
delta=np.loadtxt('./deltaregion.dat')

# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111, projection='3d')
l=30
ax1.plot(pb[:,0],pb[:,1]*0-3*l,pb[:,1],'-',c='b',linewidth=1.,alpha=1.,label=r'$\mu_B=0$')

ax1.fill_between(mub, mub*0-3*l, beta ,mub,mub*0-3*l, pb[:,1], color='skyblue', alpha=0.5, edgecolor='none',zorder=0)  
ax1.fill_between(mub2,delta*0-3*l, pb2[:,1],mub2, np.log(delta)*l, pb2[:,1] , color='#F4DDB6', alpha=0.5, edgecolor='none',zorder=0)  
ax1.quiver(0, 0-3*l, pb[0,1], 0, 6*l, 0, arrow_length_ratio=0.1, color='b', linewidth=1.)
#黑色主坐标轴（真实物理轴）
ax1.quiver(0, -3*l, 0, 700, 0, 0, arrow_length_ratio=0.05, color='k', linewidth=1.5)
ax1.quiver(0, -3*l, 0, 0, 8*l, 0, arrow_length_ratio=0.1, color='k', linewidth=1.5)
ax1.quiver(0, -3*l, 0, 0, 0, 160, arrow_length_ratio=0.15, color='k', linewidth=1.5)

ax1.set_xlim(0, 700)
ax1.set_ylim(-3*l, 5*l)
ax1.set_zlim(0, 160)
#ax1.set_box_aspect((1, 1, 1))
ax1.view_init(elev=20, azim=-45)
#ax1.view_init(elev=20, azim=-160)

# 背景轴设置
ax1.xaxis.pane.set_visible(False)
ax1.yaxis.pane.set_visible(False)
ax1.zaxis.pane.set_visible(False)
ax1.invert_yaxis()
ax1.set_xlabel(r'$\mu$', fontsize=14, color='black')
ax1.set_ylabel(r'$m_\pi$', fontsize=14, color='black')
ax1.set_zlabel(r'$T$', fontsize=14, color='black')
# ax1.set_yscale(    'symlog',
#     linthresh=1e-3,   # 线性区间半宽
#     linscale=1.0)
for axis in [ax1.xaxis, ax1.yaxis, ax1.zaxis]:
    axis.line.set_color((0.7, 0.7, 0.7))
    axis.line.set_linewidth(0.8)

ax1.tick_params(colors='0', width=0.8, labelsize=8)

#fig.subplots_adjust(top=0.9, bottom=0.14, left=0.1, right=0.95, hspace=0.35,wspace=0.35)
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax1.set_position([0.0, 0.05, 0.95, 0.95])

fig.savefig("phase_3D.pdf")
