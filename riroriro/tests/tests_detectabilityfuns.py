#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the detectabilityfuns module.
"""

import pytest
import riroriro.detectabilityfuns as det

def test_detectabilityfuns_errors():
    """
    Testing improper inputs to detectabilityfuns functions.
    """
    
    with pytest.raises(AssertionError):
        det.cdf_generator(123.4)

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