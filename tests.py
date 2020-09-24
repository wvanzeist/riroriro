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
import gwexporter as gwe
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

def test_all_errors():
    """
    Runs all tests of error-raising.
    """
    
    test_inspiralfuns()
    test_mergerfirstfuns()
    test_matchingfuns()
    test_mergersecondfuns()
    test_snrcalculatorfuns()
    test_horizondistfuns()
    test_detectabilityfuns()
    
def test_gw_numerical():
    """
    Testing whether the GW-simulating functions behave as expected for a
    specific example system.
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
    x, xtimes, dt = ins.x_integration(start_x,end_x,M,eta)
    realtimes = ins.inspiral_time_conversion(xtimes,M)
    i_phase, omega, freq = ins.inspiral_phase_freq_integration(x,dt,M)
    r, rdot = ins.radius_calculation(x,M,eta)
    A1, A2 = ins.a1_a2_calculation(r,rdot,omega,D,M,eta)
    Aorth, Adiag = ins.inspiral_strain_polarisations(A1,A2,i_phase)
    i_amp = ins.inspiral_strain_amplitude(Aorth,Adiag)
    
    i_time = ins.list_size_reducer(100,realtimes)
    i_omega = ins.list_size_reducer(100,omega)
    i_phase = ins.list_size_reducer(100,i_phase)
    i_amp = ins.list_size_reducer(100,i_amp)
    
    sfin, wqnm = me1.quasi_normal_modes(eta)
    alpha, b, C, kappa = me1.gIRS_coefficients(eta,sfin)
    fhat, m_omega = me1.merger_freq_calculation(wqnm,b,C,kappa)
    fhatdot = me1.fhat_differentiation(fhat)
    m_time = me1.merger_time_conversion(M)
    
    min_switch_ind = mat.min_switch_ind_finder(i_time,i_omega,m_time,m_omega)
    final_i_index = mat.final_i_index_finder(min_switch_ind,i_omega,m_omega)
    time_offset = mat.time_offset_finder(min_switch_ind,final_i_index,i_time,\
                                         m_time)
    i_m_time, i_m_omega = mat.time_frequency_stitching(min_switch_ind,\
                   final_i_index,time_offset,i_time,i_omega,m_time,m_omega)
    i_m_freq = mat.frequency_SI_units(i_m_omega,M)
    
    #m_phase = me2.merger_phase_calculation(min_switch_ind,final_i_index,\
    #                                       i_phase,m_omega)
    #i_m_phase = me2.phase_stitching(final_i_index,i_phase,m_phase)
    m_amp = me2.merger_strain_amplitude(min_switch_ind,final_i_index,alpha,\
                                        i_amp,m_omega,fhat,fhatdot)
    i_m_amp = me2.amplitude_stitching(final_i_index,i_amp,m_amp)
    
    try:
        gwarray = gwe.waveform_arrayer(i_m_time,i_m_freq,i_m_amp)
    except AttributeError:      #different name used in fourier branch
        gwarray = gwe.waveform_arrayer_3col(i_m_time,i_m_freq,i_m_amp)
    
    assert len(gwarray) == 16114, ('The length of the data array is not as '
                                   'expected.')
    
    assert np.sum(np.isclose(gwarray[16113,:],[6.237875956808586864e+00,\
        2.933614258599319555e+02,1.483560630863543344e-23])) == 3, ('The '
        'values in the final row of the data array are not as expected.')
                                                                    
def test_snr_numerical(inputarray,findchirp_array,noisearray):
    """
    Testing whether the SNR-calculating functions behave as expected for a
    given waveform data file.
    For inputarray and findchirp_array, please use the files
    'example_1.4_0.8.dat' and 'example_1.4_0.8_findchirp.dat', respectively,
    which are provided in the riroriro_tutorials directory.
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
    

def test_all_numerical(inputarray,findchirp_array,noisearray):
    """
    Runs all numerical tests of calculations.
    For inputarray and findchirp_array, please use the files
    'example_1.4_0.8.dat' and 'example_1.4_0.8_findchirp.dat', respectively,
    which are provided in the riroriro_tutorials directory.
    For noisearray, please use the file 'o3_l1.txt' from
    https://dcc.ligo.org/LIGO-T1500293/public
    """
    
    test_gw_numerical()
    test_snr_numerical(inputarray,findchirp_array,noisearray)
    
def test_all(inputarray,findchirp_array,noisearray):
    """
    Runs all tests.
    For inputarray and findchirp_array, please use the files
    'example_1.4_0.8.dat' and 'example_1.4_0.8_findchirp.dat', respectively,
    which are provided in the riroriro_tutorials directory.
    For noisearray, please use the file 'o3_l1.txt' from
    https://dcc.ligo.org/LIGO-T1500293/public
    """
    
    test_all_errors()
    print("Error testing completed.")
    test_all_numerical(inputarray,findchirp_array,noisearray)
    print("Numerical testing completed.")
