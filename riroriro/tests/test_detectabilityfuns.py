#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the detectabilityfuns module.
"""

import numpy as np
import pytest
import riroriro.detectabilityfuns as det

def test_detectabilityfuns_errors():
    """
    Testing improper inputs to detectabilityfuns functions.
    """
    
    with pytest.raises(AssertionError):
        det.cdf_generator(123.4)
    
    with pytest.raises(AssertionError):
        det.specific_orientation_SNR('foo',0.1,0.1,0.1,20.0)
        
    with pytest.raises(AssertionError):
        det.specific_orientation_SNR(2.0,0.1,0.1,0.1,20.0)
        
    with pytest.raises(AssertionError):
        det.specific_orientation_SNR(0.1,7.0,0.1,0.1,20.0)
        
    with pytest.raises(AssertionError):
        det.specific_orientation_SNR(0.1,0.1,2.0,0.1,20.0)
        
    with pytest.raises(AssertionError):
        det.specific_orientation_SNR(0.1,0.1,0.1,4.0,20.0)
    
    with pytest.raises(ValueError):
        det.specific_orientation_SNR(0.1,0.1,0.1,0.1,20.0,'foo')
        
    print('test_detectabilityfuns_errors completed successfully')

def test_detectabilityfuns_numerical():
    """
    Testing whether functions in detectabilityfuns behave as expected
    numerically.
    """
    
    input_SNR = 20.0 #semi-arbitrary choice
    
    Theta_CDF, min_CDF, max_CDF = det.cdf_generator()
    detectability = det.detectability_calculator(Theta_CDF,min_CDF,max_CDF,\
                                                 input_SNR)
    assert 0.358 <= detectability <= 0.362, ('The detectability %f is outside '
        'the expected range of 0.358 to 0.362.' % (detectability))
    #This is actually expected to be slightly different each time because the
    #projection function CDF is dynamically generated via random variables.
    #The included range somewhat overstates the expected variability to avoid
    #false positives of error detection.
    
    adjusted_SNR = det.specific_orientation_SNR(50.0,50.0,50.0,50.0,20.0,'deg')
                                    #semi-arbitrary parameters
    assert np.isclose(adjusted_SNR,9.011032429277753), ('The orientation-'
        'adjusted SNR is not as expected.')
    
    print('test_detectabilityfuns_numerical completed successfully')