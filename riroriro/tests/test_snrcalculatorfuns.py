#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the snrcalculatorfuns module.
"""

import numpy as np
import pytest
import riroriro.snrcalculatorfuns as snr

def test_snrcalculatorfuns_errors():
    """
    Testing improper inputs to snrcalculatorfuns functions.
    """
    
    with pytest.raises(AssertionError):
        snr.polynomial_redshift('foo')
        
    with pytest.raises(AssertionError):
        snr.polynomial_redshift(-0.5)
        
    with pytest.raises(AssertionError):
        snr.redshift_distance_adjustment([0,1,2],0.2,0.2)
        
    print('test_snrcalculatorfuns_errors completed successfully')
        
def test_snrcalculatorfuns_numerical():
    """
    Testing whether the SNR-calculating functions behave as expected for a
    given waveform data file.
    """
    
    import urllib.request
    
    #we load in the three files of example data; ideally np.genfromtxt would be
    #used for this but these files are stored on the riroriro_tutorials github
    #repository rather than local files
    
    #loading in inputarray
    inputarray_raw = urllib.request.urlopen('https://raw.githubusercontent.com'
        '/wvanzeist/riroriro_tutorials/main/example_1.4_0.8.dat')
    inputarray_lines = [line.decode("utf-8") for line in inputarray_raw]
    inputarray = np.empty((len(inputarray_lines),3))
    for i in range(len(inputarray)):
        line = inputarray_lines[i].split("\t")
        for j in range(3):
            inputarray[i,j] = float(line[j])
            
    #loading in findchirp_array
    fcarray_raw = urllib.request.urlopen('https://raw.githubusercontent.com'
        '/wvanzeist/riroriro_tutorials/main/mergerwave_1.40000_0.800000.dat')
    fcarray_lines = [line.decode("utf-8") for line in fcarray_raw]
    findchirp_array = np.empty((len(fcarray_lines),4))
    for i in range(len(findchirp_array)):
        line = fcarray_lines[i].split("\t")
        for j in range(4):
            findchirp_array[i,j] = float(line[j])
            
    #loading in noisearray
    noisearray_raw = urllib.request.urlopen('https://raw.githubusercontent.com'
        '/wvanzeist/riroriro_tutorials/main/noise_spectra/o3_l1.txt')
    noisearray_lines = [line.decode("utf-8") for line in noisearray_raw]
    noisearray = np.empty((len(noisearray_lines),2))
    for i in range(len(noisearray)):
        line = noisearray_lines[i].split("\t")
        for j in range(2):
            noisearray[i,j] = float(line[j])
    
    #now, we perform the tests
    
    d = 100.0
    
    z = snr.polynomial_redshift(d)
    adj_inputarray = snr.redshift_distance_adjustment(inputarray,d,z)
    freqmax, freqmin = snr.frequency_limits(adj_inputarray)
    fourieramp = snr.findchirp_fourier(adj_inputarray,findchirp_array,d,z)
    
    noise_freq_amp = snr.amplitude_interpolation(adj_inputarray,fourieramp,\
                                                 noisearray,freqmax,freqmin)
    liv_SNR = snr.individual_detector_SNR(noisearray,noise_freq_amp)
    #specifically for LIGO Livingston noise spectrum
    
    assert np.isclose(liv_SNR,284.81739686496985), ('The SNR value is not as '
        'expected. Please make sure you have set the inputs correctly.')
    
    print('test_snrcalculatorfuns_numerical completed successfully')