#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for Riroriro functions.
"""

import numpy as np
import pytest

import inspiralfuns as ins
import mergerfirstfuns as me1
import matchingfuns as mat
import mergersecondfuns as me2
#import gwexporter as gwe
import snrcalculatorfuns as snr
import horizondistfuns as hor
import detectabilityfuns as det

def test_inspiralfuns():
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
        
def test_mergerfirstfuns():
    """
    Testing improper inputs to mergerfirstfuns functions.
    """
    
    with pytest.raises(AssertionError):
        me1.quasi_normal_modes('foo')
    
    with pytest.raises(AssertionError):
        me1.quasi_normal_modes(0.75)
        
    with pytest.raises(AssertionError):
        me1.quasi_normal_modes(-0.1)
        
def test_matchingfuns():
    """
    Testing improper inputs to matchingfuns functions.
    """
    
    with pytest.raises(AssertionError):
        mat.min_switch_ind_finder('foo',[0,2],[0,2],[0,2])
        
def test_mergersecondfuns():
    """
    Testing improper inputs to mergersecondfuns functions.
    """
    
    with pytest.raises(AssertionError):
        me2.merger_phase_calculation(1.4,1,[0,2],[0.2])

def test_snrcalculatorfuns():
    """
    Testing improper inputs to snrcalculatorfuns functions.
    """
    
    with pytest.raises(AssertionError):
        snr.polynomial_redshift('foo')
        
    with pytest.raises(AssertionError):
        snr.polynomial_redshift(-0.5)
        
    with pytest.raises(AssertionError):
        snr.redshift_distance_adjustment([0,1,2],0.2,0.2)

def test_horizondistfuns():
    """
    Testing improper inputs to horizondistfuns functions.
    """
    
    with pytest.raises(AssertionError):
        hor.horizon_distance_calculation(np.empty((3)),np.empty((3)),\
                                         np.empty((3)),'quad')

def test_detectabilityfuns():
    """
    Testing improper inputs to detectabilityfuns functions.
    """
    
    with pytest.raises(AssertionError):
        det.cdf_generator(123.4)

def test_all():
    """
    Runs all tests.
    """
    
    test_inspiralfuns()
    test_mergerfirstfuns()
    test_matchingfuns()
    test_mergersecondfuns()
    test_snrcalculatorfuns()
    test_horizondistfuns()
    test_detectabilityfuns()