#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parts of the procedure for calculating the SNR of a gravitational waveform
(from gwexporter or otherwise), collected into modular functions.
"""

import numpy as np

def polynomial_redshift(d):
    """
    Polynomial approximation of calculating redshift corresponding to a given
    distance.
    
    Parameters
    ----------
    d: float
        A luminosity distance, in Mpc.
        
    Results
    -------
    z: float
        The redshift corresponding to the input distance.
    """
    
    #input type checking
    assert type(d) == float, 'd should be a float.'
    
    #sanity check: distance should not be negative
    assert d >= 0, 'The distance should be a positive number.'
    
    #polynomial approximation of redshift conversion
    z = 1.0832e-12*d**3 - 1.7022e-8*d**2 + 0.00021614*d
    
    return z

def redshift_distance_adjustment(inputarray,d,z):
    """
    Adjusts the frequencies and amplitudes in the input gravitational waveform
    to account for the effects of distance/redshift.
    
    Parameters
    ----------
    inputarray: numpy.ndarray
        The time, frequency and amplitude data of the gravitational waveform,
        in the format used by waveform_exporter() in gwexporter.
    d: float
        The luminosity distance to the merging binary, in Mpc.
    z: float
        The redshift corresponding to the input distance.
        
    Results
    -------
    adjustedarray: numpy.ndarray
        inputarray, but with the frequency and amplitudes adjusted.
    """
    
    #input type checking
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    assert type(d) == float, 'd should be a float.'
    assert type(z) == float, 'z should be a float.'
    
    adjustedarray = np.zeros(inputarray.shape)
    
    for i in range(inputarray.shape[0]):
        adjustedarray[i,0] = inputarray[i,0]
        adjustedarray[i,1] = inputarray[i,1] / (1+z)    #frequency redshifting
        adjustedarray[i,2] = inputarray[i,2] / (d/100)  #distance adjustment
        
    return adjustedarray

def frequency_limits(inputarray):
    """
    Calculates the upper and lower limits of the frequency of the gravitational
    waveform in inputarray, which are used by the subsequent noiseratio
    calculation.
    
    Parameters
    ----------
    inputarray: numpy.ndarray
        The time, frequency and amplitude data of the gravitational waveform;
        should have been adjusted by redshift_distance_adjustment().
        
    Results
    -------
    (freqmax,freqmin): tuple of floats
        The upper and lower limits of the waveform signal frequency,
        respectively.
    """
    
    #input type checking
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    
    freqmax=np.amax(inputarray[:,1])
    freqmin=np.amin(inputarray[:,1])
    
    return (freqmax,freqmin)