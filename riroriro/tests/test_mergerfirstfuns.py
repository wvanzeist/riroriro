#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the mergerfirstfuns module.
"""

import numpy as np
import pytest
import riroriro.mergerfirstfuns as me1

def test_mergerfirstfuns_errors():
    """
    Testing improper inputs to mergerfirstfuns functions.
    """
    
    with pytest.raises(AssertionError):
        me1.quasi_normal_modes('foo')
    
    with pytest.raises(AssertionError):
        me1.quasi_normal_modes(0.75)
        
    with pytest.raises(AssertionError):
        me1.quasi_normal_modes(-0.1)
        
    print('test_mergerfirstfuns_errors completed successfully')

def test_mergerfirstfuns_numerical():
    """
    Testing whether functions in mergerfirstfuns behave as expected
    numerically.
    """
    
    #generic system, somewhat like GW150914
    M = 60.0 #(Msun)
    eta = 0.24
    
    sfin, wqnm = me1.quasi_normal_modes(eta)
    alpha, b, C, kappa = me1.gIRS_coefficients(eta,sfin)
    fhat, m_omega = me1.merger_freq_calculation(wqnm,b,C,kappa)
    fhatdot = me1.fhat_differentiation(fhat)
    m_time = me1.merger_time_conversion(M)
    
    assert len(fhat) == 201, 'The length of fhat is not as expected.'
    
    assert np.isclose(fhatdot[200],-1.6137233253473387e-07), ('The final value'
        ' of fhatdot is not as expected.')
    
    assert np.isclose(m_time[200],0.029538), ('The final value of '
        'm_time is not as expected.')
    
    print('test_mergerfirstfuns_numerical completed successfully')