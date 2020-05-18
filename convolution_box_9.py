# -*- coding: utf-8 -*-
"""
Created on Wed May 13 22:14:39 2020

@author: Jay
"""

import numpy as np
from matplotlib.pyplot import plot, show, legend, title

def f(x):
    if(x<=1 and x>=-1):
        return 1
    else: return 0

xmin = -5.0
xmax = 5.0
n_points = 256

dx = (xmax - xmin)/(n_points-1)

x = np.linspace(xmin,xmax,n_points)
sampled_data = [f(xi) for xi in x]
sampled_data = np.pad(sampled_data,(0,n_points),mode='constant',constant_values=(0, 0)) 
    
nft = np.fft.fft(sampled_data, norm='ortho')
cnv = np.fft.ifft(nft*nft, norm='ortho')
acnv = dx*np.sqrt(2*n_points)*cnv
#acnv = np.fft.fftshift(acnv)

x_a = int(n_points/2)
x_b = int(3*n_points/2)

f_box = [f(xi) for xi in x]
plot(x,f_box,'r',markersize='7',label='Box')
plot(x,np.real(acnv[x_a:x_b]),'g',label='Convolution')
legend()
title("Convolution of box function with itself")
show()