# -*- coding: utf-8 -*-
"""
Created on Thu May 14 17:37:17 2020

@author: Jay
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def f(x,y):
    return np.exp(-(x*x + y*y))

def f_ft(kx,ky):
    return 0.5*np.exp(-(kx*kx + ky*ky)/4)

xmin = -20.0
xmax = 20.0
nx_points = 64
dx = (xmax - xmin)/(nx_points-1)

ymin = -20.0
ymax = 20.0
ny_points = 64
dy = (ymax - ymin)/(ny_points-1)

x = np.linspace(xmin,xmax,nx_points)
y = np.linspace(ymin,ymax,ny_points)
X,Y = np.meshgrid(x,y)

sampled_data = f(X,Y)
'''
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, sampled_data, 50, cmap='binary')
'''

nft = np.fft.fft2(sampled_data)
kx = 2*np.pi*np.fft.fftfreq(nx_points,dx)
ky = 2*np.pi*np.fft.fftfreq(ny_points,dy)
Kx,Ky = np.meshgrid(kx,ky)
factor = np.exp(-1j*Kx*xmin)*np.exp(-1j*Ky*ymin)

aft = dx*dy*factor*nft/(2*np.pi)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(Kx, Ky, np.real(aft), rstride=1, cstride=1,cmap='viridis', edgecolor='none')
#ax.plot_surface(Kx, Ky, f_ft(Kx,Ky),  rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#plt.title("Analytical 2dFT of gaussian")
plt.title("Numerical 2dFT of gaussian")

ax.view_init(60, 60)
fig