#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculating the SNR of a supplied gravitational waveform (from the other
scripts).

@author: woutervanzeist
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

"""
    Inputs
"""

#columns: time, frequency, amp
inputarray=np.genfromtxt("100 Mpc grid of waveforms/1.4_0.8.dat")
#closest to GW150914 system is 1.4_0.8 to 1.5_0.9

noisearray=np.genfromtxt("May 2019 ASD curves/o1.txt")		#O1 LIGO noise
#columns: frequency, amp

d=440                                               #distance (Mpc)
#default/GW150914: 440
#GW151012: 1.2_0.6, 1080   GW151226: 0.9_0.6, 450

findchirp_array = \
    np.genfromtxt("traces/MakeWaves-3/mergerwave_1.40000_0.800000.dat")
#calculated assuming 1 Mpc, hence scaling:
for i in range(findchirp_array.shape[0]):
    findchirp_array[i,3] /= d

"""
    Calculation
"""

H0=70                                   #Hubble's constant (km/s/Mpc)
c=2.99792458e8                                      #speed of light (m/s)
#z = (H0/(c/1000))*d                     #simplified redshift (valid for low z)
                                        #factor of 1000 is for m -> km 
z = 1.0832e-12*d**3 - 1.7022e-8*d**2 + 0.00021614*d #more precise (from JJ)

#first part of SNR calculation that can be outside loop because it does not
#involve the noisearray

for i in range(inputarray.shape[0]):                #implementing redshift
    inputarray[i,1] = inputarray[i,1] / (1+z)
    inputarray[i,2] = inputarray[i,2] / (d/100)     #distance adjustment
    #100 because waveforms are all generated assuming 100 MPc
    #this bit also needs to be outside all loops because its self-reference
    #means that otherwise the values will get smaller each iteration
    
for i in range(findchirp_array.shape[0]):
    findchirp_array[i,1] /= (1+z)   #redshift adjustment for trace

freqmax=np.amax(inputarray[:,1])
freqmin=np.amin(inputarray[:,1])
#used solely by noiseratio calculation but defined here for speed

#FINDCHIRP-inspired empirical approximation of Fourier transform
fourieramp = np.zeros((inputarray.shape[0]))
for i in range(len(fourieramp)):
    fourieramp[i] = (inputarray[i,2]) / (inputarray[i,1]**(11/6))
                    #*fourier_freq**(5/6)          #* 1.15 (trace scaling)
                    #11/6 in this version with adjusted traces
                    
#scaling to match Fourier-transformed waveform and FINDCHIRP trace at 10 Hz
fourier_10Hz = np.searchsorted(inputarray[:,1],10)
trace_10Hz = np.searchsorted(findchirp_array[:,1],10)
f_t_ratio = fourieramp[fourier_10Hz] / findchirp_array[trace_10Hz,3]
for i in range(len(fourieramp)):
    fourieramp[i] /= f_t_ratio #scale so amplitudes are equal at 10 Hz

smoothinput = interp1d(inputarray[:,1], fourieramp[:], kind='cubic')
#previously:
#smoothinput = interp1d(inputarray[:,1], inputarray[:,2], kind='cubic')
#interpolating input frequency-amp curve so it can be calculated for the other
#grid of frequencies used by noisearray

#rest of SNR calculation is put in loop

def noiseratio(i):                                  #h/S part of SNR equation
    if noisearray[i,0] > freqmax or noisearray[i,0] < freqmin:
            ratio = 0
            #just set ratio to zero outside the actual frequency range of the
            #gravitational waveform
            #not doing this causes searchsorted to break when searching for
            #frequencies higher than amax(inputarray)
    else:
        ratio = (smoothinput(noisearray[i,0]) / noisearray[i,1])**2
        #[i,0] is the frequency we're comparing these at
        #squaring is to convert ASD to PSD and amp to h*(f)h(f)
        #f to balance ASD units (from noise and from df) is no longer
        #needed when we have FINDCHIRP-based transform
    return ratio

#test
#a = np.zeros(noisearray.shape[0])
#for i in range(len(a)):
#    a[i] = noiseratio(i)

#plt.plot(range(len(a)),a)
#plt.plot(range(2500),a[0:2500])
#plt.plot(noisearray[:,0],a)

def df(i):                                          #derivative as differences
    #unusual averages at ends (not that these will be relevant here)
    if i == 0:
        deriv = noisearray[1,0] - noisearray[0,0]
    elif i == noisearray.shape[0] - 1:              #final value in array
        deriv = noisearray[i,0] - noisearray[i-1,0]
        #regular averages, based on those used in merger part fhatdot calculation
    else:
        deriv = 0.5*(noisearray[i+1,0] - noisearray[i-1,0])
    return deriv

#test
#a = np.zeros(2500)
#for i in range(len(a)):
#    a[i] = df(i)
    
#plt.plot(range(len(a)),a)

sumpart = np.zeros(noisearray.shape[0])             #slices of SNR integral
for i in range(len(sumpart)):
    sumpart[i] = noiseratio(i)*df(i)

#final part of Barrett et al. calculation
ind_SNR = np.sqrt(4*sum(sumpart))

print(ind_SNR)

#example graph
plt.figure(1)
plt.plot(noisearray[:,0],noisearray[:,1])
plt.plot(inputarray[:,1],fourieramp[:])
#plt.plot(findchirp_array[:,1],findchirp_array[:,2])
plt.plot(findchirp_array[:,1],findchirp_array[:,3])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Strain/√(f (Hz))')
plt.title('Fourier approx. vs. O1 noise vs. FINDCHIRP trace')
#plt.axis((0,300,10**-24,10**-21))
plt.axis((10,1000,10**-24,10**-21))