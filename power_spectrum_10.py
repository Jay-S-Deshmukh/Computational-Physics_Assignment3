# -*- coding: utf-8 -*-
"""
Created on Sat May 16 19:46:01 2020

@author: Jay
"""
import numpy as np
import matplotlib.pyplot as plt

f = open("data.txt","r").read().split('\n')

data = []
for x in f:    
    data.append(float(x))

X = np.array(data)
N = X.size

X_k = np.fft.fft(X)  
P = (1/N)*np.abs(X_k)**2
k = 2*np.pi*np.fft.fftfreq(N,1)

k_bins = 10
kmax = np.amax(k)
kmin = np.amin(k)
dk = (kmax - kmin)/(k_bins)

P_binned_vals = np.zeros((k_bins,2))
for i in range(k_bins):
    for ki in k:
        if ki >= (kmin + i*dk) and ki <= (kmin + (i+1)*dk):
            P_binned_vals[i][0] = P_binned_vals[i][0] + P[i]
            P_binned_vals[i][1] = P_binned_vals[i][1] + 1
            
P_binned = np.zeros(N)
for i in range(N):
    j = int(np.floor(5*i/N))
    P_binned[i] = P_binned_vals[j][0]/P_binned_vals[j][1]

idx = np.argsort(k)
plt.plot(k[idx],P[idx],'b',label='Unbinned')
plt.plot(k[idx],P_binned[idx],'r',label='Binned')
plt.legend()
plt.title("Power Spectrum")
plt.show()

plt.figure(figsize=(15,4))
plt.plot(X)
plt.title("Data")
plt.show()
plt.figure(figsize=(15,4))
plt.plot(X_k)
plt.title("DFT")
plt.show()