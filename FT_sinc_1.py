# -*- coding: utf-8 -*-
"""
Created on Mon May 11 12:51:05 2020

@author: Jay
"""

import numpy as np
from matplotlib.pyplot import plot, show, legend, title
#import csv

def f(x):
    if(x==0):
        return 1
    return np.sin(x)/x

def f_ft(k):
    if(k<=1 and k>=-1):
        return np.sqrt(np.pi/2)
    else: return 0

xmin = -30.0
xmax = 30.0
n_points = 64
dx = (xmax - xmin)/(n_points-1)

x = np.linspace(xmin,xmax,n_points)
sampled_data = [f(xi) for xi in x] 
  
nft = np.fft.fft(sampled_data, norm='ortho')
k = np.fft.fftfreq(n_points,dx)
k = 2*np.pi*k
factor = np.exp(-1j*k*xmin)

aft = dx*np.sqrt(n_points/(2*np.pi))*factor*nft

aft = np.fft.fftshift(aft)
k = np.fft.fftshift(k)

'''x = []
y = []

with open('sinc_fftw.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))'''

ft_analytical = [f_ft(ki) for ki in k]
plot(k,np.real(aft),'r^',label='Numerical')
plot(k,ft_analytical,'b',label='Analytical')
#plot(k,np.real(aft),'r^',label='Numpy', markersize = 10)
#plot(x,y,'g^',label='FFTW3')
legend()
title("FT of sin(x)/x")
show()