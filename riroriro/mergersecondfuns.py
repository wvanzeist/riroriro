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

def merger_polarisations(final_i_index,m_amp,m_phase,i_Aorth):
    """
    Calculating the values of the two polarisations of strain for the merger.
    
    Parameters
    ----------
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder() in matchingfuns.
    m_amp: list of floats
        The values of the amplitude of the GW strain over time for the
        merger/ringdown portion, from merger_strain_amplitude().
    m_phase: list of floats
        Values of orbital phase over time for the merger/ringdown portion, from
        merger_phase_calculation().
    i_Aorth: list of floats
        The values of the orthogonal/plus polarisation of strain over time for
        the inspiral portion, from inspiral_strain_polarisations() in
        inspiralfuns.
        
    Returns
    -------
    [m_Aorth,m_Adiag]: list of lists of floats
        The first list is the values of the orthogonal/plus polarisation of
        strain over time, the second list is the diagonal/cross polarisation.
    """
    
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(m_amp) == list, 'm_amp should be a list.'
    assert type(m_phase) == list, 'm_phase should be a list.'
    assert type(i_Aorth) == list, 'i_Aorth should be a list.' 
    
    m_Aorth = np.empty((len(m_amp)))            #orthogonal/plus polarisation
    m_Adiag = np.empty((len(m_amp)))            #diagonal/cross polarisation
    
    #adjusting phase to ensure continuity between inspiral and merger
    #polarisations at the switching point
    comparison_phase_1 = (1/2) * np.arccos(i_Aorth[final_i_index] / m_amp[0])
    #phase so that m_Aorth = i_Aorth at the switching point
    #this ensures continuity in amplitude, but for continuity in gradient of
    #amplitude the required value could be either this or π - this:
    comparison_phase_2 = np.pi - comparison_phase_1
    #now a function is required that tests both and sees which one has the
    #correct sign of the gradient
    #the following "sign" constants indicate whether the function is
    #increasing or decreasing at the switching point
    sign_i_Aorth = np.sign(i_Aorth[final_i_index] - i_Aorth[final_i_index-1])
    sign_comparison = np.sign(np.cos(2*(m_phase[1] - m_phase[0] + \
        comparison_phase_1)) - np.cos(2*comparison_phase_1))
    if sign_i_Aorth == sign_comparison: #signs match, phase_1 is correct
        comparison_phase = comparison_phase_1
    else:                               #signs don't match, phase_2 is correct
        comparison_phase = comparison_phase_2
    
    phase_difference = m_phase[0] - comparison_phase
    adjusted_m_phase = np.empty((len(m_phase)))
    for i in range(len(m_amp)):    #adjusting phase for polarisation continuity
        adjusted_m_phase[i] = m_phase[i] - phase_difference
    
    for i in range(len(m_amp)):
        m_Aorth[i] = m_amp[i] * np.cos(2*adjusted_m_phase[i])
        m_Adiag[i] = m_amp[i] * np.sin(2*adjusted_m_phase[i])
    
    #output type conversion
    m_Aorth = list(m_Aorth)
    m_Adiag = list(m_Adiag)
    
    return [m_Aorth,m_Adiag]

def polarisation_stitching(final_i_index,i_Aorth,i_Adiag,m_Aorth,m_Adiag):
    """
    Stitching together the inspiral and merger/ringdown portions of the
    polarisation lists to give combined lists with the correct matching.
    
    Parameters
    ----------
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder() in matchingfuns.
    i_Aorth: list of floats
        The values of the orthogonal/plus polarisation of strain over time for
        the inspiral portion, from inspiral_strain_polarisations() in
        inspiralfuns.
    i_Adiag: list of floats
        The values of the diagonal/cross polarisation of strain over time for
        the inspiral portion, from inspiral_strain_polarisations() in
        inspiralfuns.
    m_Aorth: list of floats
        The values of the orthogonal/plus polarisation of strain over time for
        the merger/ringdown portion, from merger_polarisations().
    m_Adiag: list of floats
        The values of the diagonal/cross polarisation of strain over time for
        the merger/ringdown portion, from merger_polarisations().
        
    Returns
    -------
    [i_m_Aorth,i_m_Adiag]: list of lists of floats
        The first list is the combined orthogonal/plus polarisation values, the
        second list is the combined diagonal/cross polarisation values.
    """
    
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_Aorth) == list, 'i_Aorth should be a list.'
    assert type(i_Adiag) == list, 'i_Adiag should be a list.'
    assert type(m_Aorth) == list, 'm_Aorth should be a list.'
    assert type(m_Adiag) == list, 'm_Adiag should be a list.'
    
    i_m_Aorth = np.concatenate((i_Aorth[:final_i_index],m_Aorth))
    i_m_Adiag = np.concatenate((i_Adiag[:final_i_index],m_Adiag))
    
    #output type conversion
    i_m_Aorth = list(i_m_Aorth)
    i_m_Adiag = list(i_m_Adiag)
    
    return [i_m_Aorth,i_m_Adiag]
