#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the inspiralfuns module.
"""

import numpy as np
import pytest
import riroriro.inspiralfuns as ins

def test_inspiralfuns_errors():
    """
    Testing improper inputs to inspiralfuns functions.
    """
    
    with pytest.raises(AssertionError):
        ins.get_M_and_eta(m1='foo',m2='bar')
        
    with pytest.raises(TypeError):
        ins.get_M_and_eta(m1=30.0,m2=15.0,q=0.5)
        
    with pytest.raises(TypeError):
        ins.get_M_and_eta(m1=30.0)
        
    with pytest.raises(ValueError):
        ins.get_M_and_eta(logMc=30.0,q=1.5)
        
    with pytest.raises(AssertionError):
        ins.startx(30.0,'foo')
        
    with pytest.raises(ValueError):
        ins.endx(0.1,'not_BH_or_NS')
        
    print('test_inspiralfuns_errors completed successfully')

def test_inspiralfuns_numerical():
    """
    Testing whether functions in inspiralfuns behave as expected numerically.
    """
    
    #The same parameters, function calls etc. as in the tutorial
    #example_GW.ipynb are used.
    
    #GW150914-like
    logMc=1.4
    q=0.8
    
    #defaults
    flow=10.0           #(Hz)
    merger_type='BH'
    D=100.0             #(Mpc)
    
    M, eta = ins.get_M_and_eta(logMc=logMc,q=q)
    start_x = ins.startx(M,flow)
    end_x = ins.endx(eta,merger_type)
    x, xtimes, dt = ins.PN_parameter_integration(start_x,end_x,M,eta)
    
    assert len(x) == 2133591, 'The length of x is not as expected.'
    
    #realtimes = ins.inspiral_time_conversion(xtimes,M)
    i_phase, omega, freq = ins.inspiral_phase_freq_integration(x,dt,M)
    r, rdot = ins.radius_calculation(x,M,eta)
    A1, A2 = ins.a1_a2_calculation(r,rdot,omega,D,M,eta)
    Aorth, Adiag = ins.inspiral_strain_polarisations(A1,A2,i_phase)
    i_amp = ins.inspiral_strain_amplitude(Aorth,Adiag)
    
    assert np.isclose(i_amp[2133590],2.0433359836624698e-20), ('The final '
       'value of i_amp is not as expected.')
    
    print('test_inspiralfuns_numerical completed successfully')