#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the horizondistfuns module.
"""

import numpy as np
import pytest
import horizondistfuns as hor

def test_horizondistfuns_errors():
    """
    Testing improper inputs to horizondistfuns functions.
    """
    
    with pytest.raises(AssertionError):
        hor.horizon_distance_calculation(np.empty((3)),np.empty((3)),\
                                         np.empty((3)),'quad')

def test_horizondistfuns_numerical(inputarray,findchirp_array,noisearray):
    """
    Testing whether functions in horizondistfuns behave as expected
    numerically.
    For inputarray and findchirp_array, please use the files
    'example_1.4_0.8.dat' and
    'findchirp_traces/mergerwave_1.40000_0.800000.dat', respectively, which are
    provided in the riroriro_tutorials directory.
    For noisearray, please use the file 'o3_l1.txt' from
    https://dcc.ligo.org/LIGO-T1500293/public
    """
    
    horizon_distance = hor.horizon_distance_calculation(inputarray,\
        findchirp_array,[noisearray],'mean')
    #only a single detector, so choice of method is irrelevant
    
    assert np.isclose(horizon_distance,2340.33203125), ('The horizon distance '
                                                        'is not as expected.')