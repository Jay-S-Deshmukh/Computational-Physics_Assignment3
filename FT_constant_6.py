# -*- coding: utf-8 -*-
"""
Created on Mon May 11 21:25:17 2020

@author: Jay
"""

import numpy as np
from matplotlib.pyplot import plot, show, title, legend

def f(x):
    return 1

xmin = -50.0
xmax = 50.0
n_points = 256

dx = (xmax - xmin)/(n_points-1)

x = np.fft.fftshift(np.linspace(xmin,xmax,n_points))
sampled_data = [f(xi) for xi in x]
    
nft = np.fft.fft(sampled_data, norm='ortho')
karr = np.fft.fftfreq(n_points,dx)
karr = 2*np.pi*karr
factor = np.exp(-1j*karr*xmin)

aft = dx*np.sqrt(n_points/(2*np.pi))*factor*nft

plot(karr,np.abs(aft),'r^',label='Numerical')
legend()
title("FT of constant")
show()