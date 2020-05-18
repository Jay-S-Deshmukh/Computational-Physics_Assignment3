# -*- coding: utf-8 -*-
"""
Created on Sat May  9 21:13:57 2020

@author: Dell
"""

import numpy as np
from matplotlib.pyplot import plot, show, title, legend
import csv

def f(x):
    return np.exp(-x*x)

def f_ft(k):
    return np.exp(-k*k/4)/np.sqrt(2)

xmin = -20.0
xmax = 20.0
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

x = []
y = []

with open('gauss_fftw.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))


ft_analytical = [f_ft(ki) for ki in k]
plot(x,y,'r^',label='FFTW',markersize='6')
#plot(k,aft,'r.',label='Numpy',markersize='9')
plot(k,ft_analytical,'b',label='Analytical',markersize='6')
legend()
title("FT of e^(-x*x)")