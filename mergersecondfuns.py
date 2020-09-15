#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The post-matching parts of the procedure for simulating the merger/ringdown
portions of gravitational waves from binary black holes, collected into modular
functions.
"""

import numpy as np

def merger_phase_calculation(min_switch_ind,final_i_index,i_phase,m_omega):
    """
    Calculation of the orbital phase for the merger/ringdown portion, based
    on Buskirk et al. (2019) equation 21.
    
    Parameters
    ----------
    min_switch_ind: int
        The index in the merger/ringdown data where the switch from inspiral to
        merger/ringdown should occur, from min_switch_ind_finder() in
        matchingfuns.
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder() in matchingfuns.
    i_phase: list of floats
        Values of orbital phase at each timestep for the inspiral portion, from
        inspiral_phase_freq_integration() in inspiralfuns.
    m_omega: list of floats
        Values of angular frequency over time for the merger/ringdown portion,
        from merger_freq_calculation() in mergerfirstfuns.
        
    Returns
    -------
    m_phase: list of floats
        Values of orbital phase over time for the merger/ringdown portion.
    """
    
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_phase) == list, 'i_phase should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    m_phase = np.empty((201 - min_switch_ind))
    #starts at min_switch_ind instead of 0 because the part of the
    #merger/ringdown data before that gets cut off by the matching
    m_phase[0] = i_phase[final_i_index]
    #matching phase at start of merger/ringdown to phase at end of inspiral
    for i in range(min_switch_ind + 1,201):
        m_phase[i - min_switch_ind] = m_phase[i - min_switch_ind - 1] + \
            m_omega[i]
                                        #Euler integration Buskirk of eq. 21
    
    #output type conversion
    m_phase = list(m_phase)
    
    return m_phase

def phase_stitching(final_i_index,i_phase,m_phase):
    """
    Stitching together the inspiral and merger/ringdown portions of the phase
    lists to give a combined list with the correct matching.
    
    Parameters
    ----------
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder() in matchingfuns.
    i_phase: list of floats
        Values of orbital phase at each timestep for the inspiral portion, from
        inspiral_phase_freq_integration() in inspiralfuns.
    m_phase: list of floats
        Values of orbital phase over time for the merger/ringdown portion, from
        merger_phase_calculation().
        
    Returns
    -------
    i_m_phase: list of floats
        Values of orbital phase over time for the entire duration of the
        gravitational waveform.
    """
    
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_phase) == list, 'i_phase should be a list.'
    assert type(m_phase) == list, 'm_phase should be a list.'
    
    i_m_phase = np.concatenate((i_phase[:final_i_index],m_phase))
    
    #output type conversion
    i_m_phase = list(i_m_phase)
    
    return i_m_phase

def merger_strain_amplitude(min_switch_ind,final_i_index,alpha,i_amp,m_omega,\
                            fhat,fhatdot):
    """
    Calculating the amplitude of strain for the merger/ringdown portion, based
    on Buskirk et al. (2019) equation 16.
    
    Parameters
    ----------
    min_switch_ind: int
        The index in the merger/ringdown data where the switch from inspiral to
        merger/ringdown should occur, from min_switch_ind_finder() in
        matchingfuns.
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder() in matchingfuns.
    alpha: float
        A gIRS coefficient, from gIRS_coefficients() in mergerfirstfuns.
    i_amp: list of floats
        The values of the amplitude of the GW strain over time for the inspiral
        portion, from inspiral_strain_amplitude() in inspiralfuns.
    m_omega: list of floats
        Values of angular frequency over time for the merger/ringdown portion,
        from merger_freq_calculation() in mergerfirstfuns.
    fhat: list of floats
        Values of a sort of frequency parameter called fhat (f^) over time,
        from merger_freq_calculation() in mergerfirstfuns.
    fhatdot: list of floats
        Values of the time-derivative of fhat over time, from
        fhat_differentiation() in mergerfirstfuns.
        
    Returns
    -------
    m_amp: list of floats
        The values of the amplitude of the GW strain over time for the
        merger/ringdown portion.
    """
    
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(alpha) == float, 'alpha should be a float.'
    assert type(i_amp) == list, 'i_amp should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    assert type(fhat) == list, 'fhat should be a list.'
    assert type(fhatdot) == list, 'fhatdot should be a list.'
    
    m_amp = np.empty((201 - min_switch_ind))
    for i in range(min_switch_ind,201):         #initial unscaled calculation
        m_amp[i - min_switch_ind] = (1/(2e0*m_omega[i])) * ((abs(fhatdot[i]))/\
             (1 + alpha*(fhat[i]**2 - fhat[i]**4)))**(1/2)
                                                #Buskirk eq. 16 (note 2e0)           
        
    matching_amplitude = i_amp[final_i_index]   #amplitude at end of inspiral
    scaling_ratio = matching_amplitude / m_amp[0]
    for i in range(len(m_amp)):
        m_amp[i] = scaling_ratio * m_amp[i]
        #rescaling for continuity between inspiral and merger/ringdown
        
    #output type conversion
    m_amp = list(m_amp)
    
    return m_amp

def amplitude_stitching(final_i_index,i_amp,m_amp):
    """
    Stitching together the inspiral and merger/ringdown portions of the
    amplitude lists to give a combined list with the correct matching.
    
    Parameters
    ----------
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder() in matchingfuns.
    i_amp: list of floats
        The values of the amplitude of the GW strain over time for the inspiral
        portion, from inspiral_strain_amplitude() in inspiralfuns.
    m_amp: list of floats
        The values of the amplitude of the GW strain over time for the
        merger/ringdown portion, from merger_strain_amplitude().
        
    Returns
    -------
    i_m_amp: list of floats
        The values of the amplitude of the GW strain over time for the entire
        duration of the gravitational waveform.
    """
    
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_amp) == list, 'i_amp should be a list.'
    assert type(m_amp) == list, 'm_amp should be a list.'
    
    i_m_amp = np.concatenate((i_amp[:final_i_index],m_amp))
    
    #output type conversion
    i_m_amp = list(i_m_amp)
    
    return i_m_amp