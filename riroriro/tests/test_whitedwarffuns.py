#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the whitedwarffuns module.
"""

import numpy as np
import pytest
from random import random
import whitedwarffuns as wdf #put riroriro.whitedwarffuns when publishing

def test_whitedwarffuns_errors():
    """
    Testing improper inputs to whitedwarffuns functions.
    """
    
    with pytest.raises(AssertionError):
        wdf.wd_amplitude(1.0,1.0,[])
        
    with pytest.raises(AssertionError):
        wdf.characteristic_strain(1.0,1.0,[])
        
    with pytest.raises(AssertionError):
        wdf.wd_polarisations(1.0,1.0,1.0,1.0,1.0,1.0,1.0)
        
    with pytest.raises(AssertionError):
        wdf.wd_inclination([1.0],cosiota=1.0)
        
    with pytest.raises(TypeError):
        wdf.wd_inclination(1.0,1.0)
        
    with pytest.raises(AssertionError):
        wdf.wd_inclination_polarisations(1.0,[1.0],cosiota=1.0)
        
    with pytest.raises(TypeError):
        wdf.wd_inclination_polarisations([1.0],[1.0])
        
    with pytest.raises(TypeError):
        wdf.wd_inclination_polarisations([1.0],[1.0],foo=1.0)
        
    with pytest.raises(TypeError):
        wdf.wd_inclination_polarisations([1.0],[1.0],cosiota=1.0,iota=1.0)
        
    with pytest.raises(AssertionError):
        wdf.wd_inclination_polarisations([1.0],[1.0],iota=4.0)
        
    with pytest.raises(AssertionError):
        wdf.wd_inclination_polarisations([1.0],[1.0],cosiota=2.0)
        
    with pytest.raises(AssertionError):
        wdf.lisa_rotation([1.0],[1.0],1.0)
        
    with pytest.raises(AssertionError):
        wdf.wd_binary_vectors(1.0,1.0,1.0,[])
        
    with pytest.raises(AssertionError):
        wdf.lisa_rotation(1.0,[1.0],[1.0])
    
    with pytest.raises(AssertionError):
        wdf.lisa_phase_modulation([1.0],[1.0],[1.0],[1.0])
    
    with pytest.raises(AssertionError):
        wdf.lisa_frequency_modulation([1.0],1.0,[1.0],1.0)
        
    with pytest.raises(AssertionError):
        wdf.lisa_amplitude_modulation([1.0],[1.0],[1.0],[1.0])
        
    with pytest.raises(AssertionError):
        wdf.lisa_detector_response([1.0],[1.0],[1.0],[1.0],[1.0])
        
    print('test_whitedwarffuns_errors completed successfully')
    
def test_whitedwarffuns_numerical():
    """
    Testing whether functions in whitedwarffuns behave as expected numerically.
    """
    
    amp = wdf.wd_amplitude(0.5,0.01,100.0) #arbitrary "realistic" parameters
    assert np.isclose(amp,1.7377730758133952e-20), ('The instantaneous '
        'amplitude is not as expected.')
    
    #L should be a unit vector for any valid input orientation
    #(the same is not the case for P as it is a cross product of unit vectors)
    theta = random()*np.pi
    phi = random()*2*np.pi
    iota = random()*np.pi
    chi = random()*np.pi
    L, P = wdf.wd_binary_vectors(theta,phi,iota,chi)
    L_size = np.sqrt(sum(i*i for i in L))
    assert np.isclose(L_size,1), ('L is not a unit vector for this '
        'orientation: %f, %f, %f, %f.' % (theta,phi,iota,chi))
    
    #testing P for a semi-arbitrary valid input
    L_2, P_2 = wdf.wd_binary_vectors(0.7,1.1,0.5,0.3)
    assert np.isclose(P_2[0],-0.4573372405242637), ('The x-component of P is '
                                                    'not as expected.')
    
    #testing the beam pattern for a semi-arbitrary input
    Fplus, Fcross = wdf.instantaneous_beam_pattern(1.7,0.3,0.9)
    assert np.isclose(Fplus,-0.02446701323291403), 'Fplus is not as expected.'
    assert np.isclose(Fcross,0.4250762606055975), 'Fcross is not as expected.'
    
    print('test_whitedwarffuns_numerical completed successfully')
