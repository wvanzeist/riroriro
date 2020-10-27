#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parts of the procedure for calculating the horizon distance of a gravitational
waveform (from gwexporter or otherwise), collected into modular functions.
"""

import numpy as np

import riroriro.snrcalculatorfuns as snr

def compact_SNR_calculation(inputarray,findchirp_array,noisearray_list,method,\
                            d):
    """
    Runs through all of the functions of snrcalculatorfuns to obtain a SNR from
    an individual detector. This function is mainly included not to be called
    directly, but rather by horizon_distance_calculation().
    
    Parameters
    ----------
    inputarray: numpy.ndarray
        The time, frequency and amplitude data of the gravitational waveform,
        in the format used by waveform_exporter() in gwexporter.
    findchirp_array: numpy.ndarray
        The array output by FINDCHIRP. The second column is frequency, the
        fourth is (Fourier-transformed) strain amplitude, the other columns
        are irrelevant. A grid of sample findchirp_arrays can be found at
        https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
    noisearray_list: list of numpy.ndarrays
        Each item in this list should be an array describing the noise spectrum
        of a detector; in each noise spectrum, it is assumed that frequency
        values are in the first column and ASD noise levels in the second.
    method: str
        If 'quad', returns the quadrature SNR across the detectors in
        noisearray_list. If 'mean', returns the mean of the SNRs with each
        individual detector (simulating one random detector in operation). If
        only one detector is included in noisearray_list, these methods are
        equivalent.
    d: float
        The luminosity distance to the merging binary, in Mpc.
        
    Returns
    -------
    final_SNR: float
        The SNR of the simulated gravitational waveform, for the detectors in
        noisearray and assuming optimal alignment.
    """
    
    #no input type checks here because the functions of snrcalculatorfuns (and
    #horizon_distance_calculation) already include them
    
    z = snr.polynomial_redshift(d)
    inputarray_2 = snr.redshift_distance_adjustment(inputarray,d,z)
    freqmax, freqmin = snr.frequency_limits(inputarray_2)
    fourieramp = snr.findchirp_fourier(inputarray_2,findchirp_array,d,z)
    
    ind_SNR = np.empty((len(noisearray_list)))
    for i in range(len(noisearray_list)):
        noisearray = noisearray_list[i]
        noise_freq_amp = snr.amplitude_interpolation(inputarray_2,fourieramp, \
                                                 noisearray,freqmax,freqmin)
        ind_SNR[i] = snr.individual_detector_SNR(noisearray,noise_freq_amp)
    
    if method == 'quad':                #quadrature SNR
        final_SNR = np.sqrt(sum(i*i for i in ind_SNR)) #sqrt of sum of squares
        final_SNR = float(final_SNR)    #becomes numpy.float64 otherwise
    elif method == 'mean':              #mean-of-individual SNR
        final_SNR = sum(ind_SNR)/len(ind_SNR)
    else:
        raise ValueError('method must be either \'quad\' or \'mean\'.')
        
    return final_SNR

def horizon_distance_calculation(inputarray,findchirp_array,noisearray_list,\
                                 method):
    """
    Calculates the horizon distance (maximum distance at which something can
    be observed) given optimal alignment for a given merger.
    
    Parameters
    ----------
    inputarray: numpy.ndarray
        The time, frequency and amplitude data of the gravitational waveform,
        in the format used by waveform_exporter() in gwexporter.
    findchirp_array: numpy.ndarray
        The array output by FINDCHIRP. The second column is frequency, the
        fourth is (Fourier-transformed) strain amplitude, the other columns
        are irrelevant. A grid of sample findchirp_arrays can be found at
        https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
    noisearray_list: list of numpy.ndarrays
        Each item in this list should be an array describing the noise spectrum
        of a detector; in each noise spectrum, it is assumed that frequency
        values are in the first column and ASD noise levels in the second.
    method: str
        If 'quad', uses the quadrature SNR across the detectors in
        noisearray_list. If 'mean', uses the mean of the SNRs with each
        individual detector (simulating one random detector in operation). If
        only one detector is included in noisearray_list, these methods are
        equivalent.
        
    Returns
    -------
    horizon_dist: float
        The horizon distance of the given merger, for the given detector(s).
    """
    
    #input type checking ('method' addressed further on)
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    assert type(findchirp_array) == np.ndarray, ('findchirp_array should be an'
                                                 ' array.')
    assert type(noisearray_list) == list, ('noisearray_list should be a list. '
         'Even if using a single detector, please still wrap the noisearray '
         'with [].')
    
    #initialising calculation thresholds
    horizon_SNR = 8                         #minimum observable SNR
    SNR_accuracy = 0.01                 #how close to SNR we call satisfactory
    dist = 10.0                             #low starting distance in Mpc
    check_SNR = []                          #used for comparing to thresholds
    
    #initial order of magnitude estimating loop
    while True:                             #loop terminated from inside
        check_SNR = compact_SNR_calculation(inputarray,findchirp_array,\
                                            noisearray_list,method,dist)
        if check_SNR < horizon_SNR:
            break
        dist *= 10   #keep increasing until we cross horizon for the first time
        
    #now, we set two bounds we'll keep adjusting until we get close enough to
    #the horizon SNR
    upper_bound = dist
    lower_bound = dist/10
    
    #we check halfway between the bounds, adjusting the bounds each time to
    #keep the horizon distance between them and tighten the interval, until we
    #get close enough to the horizon SNR
    while True:                             #loop terminated from inside
        test_dist = (upper_bound + lower_bound)/2
        check_SNR = compact_SNR_calculation(inputarray,findchirp_array,\
                                            noisearray_list,method,test_dist)
        if abs(check_SNR - horizon_SNR) <= SNR_accuracy:
            break                           #end if sufficiently close
        elif check_SNR < horizon_SNR:
            upper_bound = test_dist         #horizon distance must be smaller
        elif check_SNR > horizon_SNR:
            lower_bound = test_dist         #horizon distance must be greater
            
    #if SNR is close enough to threshold, we call the current distance the
    #horizon distance
    horizon_dist = test_dist
    
    return horizon_dist
