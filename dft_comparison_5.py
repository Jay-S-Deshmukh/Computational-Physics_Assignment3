# -*- coding: utf-8 -*-
"""
Created on Mon May 11 13:26:17 2020

@author: Jay
"""
import numpy as np
from matplotlib.pyplot import plot, show, legend, title
import time
pi = np.pi

def f(x):
    return 0.0000009989738464355469*x*np.log2(x)

def dft_manual(x):
    
    n = np.size(x)
    w = np.zeros(n)
    
    for i in range(n):
        for j in range(n):
            w[i] = w[i] + x[i]*np.exp(-1j*2*pi*i*j/n)
            
    return w/np.sqrt(n)

def dft_numpy(x):
    
    return np.fft.fft(x)

manual_time = np.zeros(100)
fft_time = np.zeros(100)
n_array = np.arange(1,101)

for n in n_array:
    
    x = np.random.randint(100,size=n)
    
    start1 = time.time()
    w1=dft_manual(x)
    end1 = time.time()
    manual_time[n-1] = end1-start1
    
    start2 = time.time()
    w2=dft_numpy(x)
    end2 = time.time()
    
    fft_time[n-1] = end2-start2
    
plot(n_array,manual_time,'b^',label='Manual')
plot(n_array,fft_time,'r^',label='FFT')
#plot(n_array,f(n_array),label='nlogn')
title("Time taken vs Array size")
legend()
show()
    
    