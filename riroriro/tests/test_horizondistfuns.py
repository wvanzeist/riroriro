#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the horizondistfuns module.
"""

import numpy as np
import pytest
import riroriro.horizondistfuns as hor

def test_horizondistfuns_errors():
    """
    Testing improper inputs to horizondistfuns functions.
    """
    
    with pytest.raises(AssertionError):
        hor.horizon_distance_calculation(np.empty((3)),np.empty((3)),\
                                         np.empty((3)),'quad')
            
    print('test_horizondistfuns_errors completed successfully')

def test_horizondistfuns_numerical():
    """
    Testing whether functions in horizondistfuns behave as expected
    numerically.
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
    
    horizon_distance = hor.horizon_distance_calculation(inputarray,\
        findchirp_array,[noisearray],'mean')
    #only a single detector, so choice of method is irrelevant
    
    assert np.isclose(horizon_distance,2340.33203125), ('The horizon distance '
        'is not as expected. Please make sure you have set the inputs '
        'correctly.')
    
    print('test_horizondistfuns_numerical completed successfully')
