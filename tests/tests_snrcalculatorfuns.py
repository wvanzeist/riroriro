#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the snrcalculatorfuns module.
"""

import numpy as np
import pytest
import snrcalculatorfuns as snr

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
        
def test_snrcalculatorfuns_numerical(inputarray,findchirp_array,noisearray):
    """
    Testing whether the SNR-calculating functions behave as expected for a
    given waveform data file.
    For inputarray and findchirp_array, please use the files
    'example_1.4_0.8.dat' and
    'findchirp_traces/mergerwave_1.40000_0.800000.dat', respectively, which are
    provided in the riroriro_tutorials directory.
    For noisearray, please use the file 'o3_l1.txt' from
    https://dcc.ligo.org/LIGO-T1500293/public
    """
    
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
                                                    'expected.')