#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parts of the procedure for matching together the inspiral and merger/ringdown
portions of the gravitational waveforms of binary black holes, collected into
modular functions.
"""

import numpy as np

def MQdiff(i,i_time,i_omega,m_time,m_omega):
    """
    A function that calculates a "matching quantity" that is subsequently used
    to determine the best offset and switching point for matching together the
    inspiral and merger/ringdown portions (the precise method was revised
    several times). This function should usually not be called directly, but
    rather by min_switch_ind_finder().
    
    Parameters
    ----------
    i: int
        An index in the range of merger/ringdown data points.
    i_time: list of floats
        Real time values for the inspiral portion, from
        inspiral_time_conversion() in inspiralfuns.
    i_omega: list of floats
        Values of angular frequency over time for the inspiral portion, from
        inspiral_phase_freq_integration() in inspiralfuns.
    m_time: list of floats
        Real time values for the merger/ringdown portion, from
        merger_time_conversion() in mergerfirstfuns.
    m_omega: list of floats
        Values of angular frequency over time for the merger/ringdown portion,
        from merger_freq_calculation() in mergerfirstfuns.
        
    Returns
    -------
    df_diff: float (or nan)
        A value describing the difference in gradient of frequency between the
        inspiral and merger/ringdown portions for points of matching frequency.
    """
    
    #input type checking
    assert type(i) == int, 'i should be an int.'
    assert type(i_time) == list, 'i_time should be a list.'
    assert type(i_omega) == list, 'i_omega should be a list.'
    assert type(m_time) == list, 'm_time should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    #making sure i is in the range of indices used for the merger/ringdown
    assert 0 <= i <= 200, 'i should be in range(201).'
    
    #with this MQ method, we iterate for each frequency value in the merger
    #part, look for the closest frequency in the inspiral part, measure the df
    #difference between these points in the two waveforms, and then select the
    #frequency with the minimum df difference as the switching point, adjusting
    #the waveform times appropriately so the combined waveform is continuous in
    #f and t
    
    try:
        closest_index = np.searchsorted(i_omega, m_omega[i], side='right')
        #with this method, we use searchsorted in the frequency domain instead
        #of the time domain
        #it does assume frequency increases monotonously, which it should
        df_diff = abs((m_omega[i] - m_omega[i-1])/(m_time[i] - m_time[i-1]) - \
                      (i_omega[closest_index] - i_omega[closest_index - 1])/ \
                      (i_time[closest_index] - i_time[closest_index - 1]))
    except IndexError:
        df_diff = np.nan
        #to get rid of known searchsorted errors for some irrelevant
        #frequencies
        
    return df_diff

def min_switch_ind_finder(i_time,i_omega,m_time,m_omega):
    """
    Finds the index in the merger/ringdown data where the switch from inspiral
    to merger/ringdown should occur, as part of the matching process.
    
    Parameters
    ----------
    i_time: list of floats
        Real time values for the inspiral portion, from
        inspiral_time_conversion() in inspiralfuns.
    i_omega: list of floats
        Values of angular frequency over time for the inspiral portion, from
        inspiral_phase_freq_integration() in inspiralfuns.
    m_time: list of floats
        Real time values for the merger/ringdown portion, from
        merger_time_conversion() in mergerfirstfuns.
    m_omega: list of floats
        Values of angular frequency over time for the merger/ringdown portion,
        from merger_freq_calculation() in mergerfirstfuns.
        
    Returns
    -------
    min_switch_ind: int
        The index in the merger/ringdown data where the switch from inspiral to
        merger/ringdown should occur.
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'
        
    minMQ = 1000                                    #deliberately overly high
    MQ = []                                         #holds MQ values
    min_switch_ind = []
    
    #(method description is in MQdiff() code)
    
    for i in range(1,len(m_time)):
        #MQ calculation at each point in merger waveform
        #EXCEPT 0, because then the [i-1] in MQdiff rolls over to the last
        #value in the array (and the first index shouldn't be the switch point
        #anyway)
        MQ = MQdiff(i,i_time,i_omega,m_time,m_omega)
        if MQ < minMQ:
            minMQ = MQ                              #update minimum
            min_switch_ind = i                      #time of new minimum
            
    return min_switch_ind

def final_i_index_finder(min_switch_ind,i_omega,m_omega):
    """
    Finds what the last index in the inspiral data before the switch to the
    merger/ringdown should be, as part of the matching process.
    
    Parameters
    ----------
    min_switch_ind: int
        The index in the merger/ringdown data where the switch from inspiral to
        merger/ringdown should occur, from min_switch_ind_finder().
    i_omega: list of floats
        Values of angular frequency over time for the inspiral portion, from
        inspiral_phase_freq_integration() in inspiralfuns.
    m_omega: list of floats
        Values of angular frequency over time for the merger/ringdown portion,
        from merger_freq_calculation() in mergerfirstfuns.
        
    Returns
    -------
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown.
    """
    
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(i_omega) == list, 'i_omega should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    final_i_index = np.searchsorted(i_omega, m_omega[min_switch_ind], \
                                    side='right')
    
    #output type conversion
    final_i_index = int(final_i_index)    
    
    return final_i_index

def time_offset_finder(min_switch_ind,final_i_index,i_time,m_time):
    """
    Calculates what the offset between the time values of the inspiral and
    merger/ringdown portions should be to match them together.
    
    Parameters
    ----------
    min_switch_ind: int
        The index in the merger/ringdown data where the switch from inspiral to
        merger/ringdown should occur, from min_switch_ind_finder().
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder().
    i_time: list of floats
        Real time values for the inspiral portion, from
        inspiral_time_conversion() in inspiralfuns.
    m_time: list of floats
        Real time values for the merger/ringdown portion, from
        merger_time_conversion() in mergerfirstfuns.
        
    Returns
    -------
    time_offset: float
        The offset between the time values of the inspiral and merger/ringdown
        portions.
    """
    
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_time) == list, 'i_time should be a list.'
    assert type(m_time) == list, 'm_time should be a list.'
    
    time_offset = i_time[final_i_index] - m_time[min_switch_ind]
    
    #output type conversion
    time_offset = float(time_offset)
    
    return time_offset

def time_frequency_stitching(min_switch_ind,final_i_index,time_offset,i_time,\
                             i_omega,m_time,m_omega):
    """
    Stitches together the inspiral and merger/ringdown portions of the time and
    angular frequency lists to give combined lists for these with the correct
    matching.
    
    Parameters
    ----------
    min_switch_ind: int
        The index in the merger/ringdown data where the switch from inspiral to
        merger/ringdown should occur, from min_switch_ind_finder().
    final_i_index: int
        The last index in the inspiral data before the switch to the merger/
        ringdown, from final_i_index_finder().
    time_offset: float
        The offset between the time values of the inspiral and merger/ringdown
        portions, from time_offset_finder().
    i_time: list of floats
        Real time values for the inspiral portion, from
        inspiral_time_conversion() in inspiralfuns.
    i_omega: list of floats
        Values of angular frequency over time for the inspiral portion, from
        inspiral_phase_freq_integration() in inspiralfuns.
    m_time: list of floats
        Real time values for the merger/ringdown portion, from
        merger_time_conversion() in mergerfirstfuns.
    m_omega: list of floats
        Values of angular frequency over time for the merger/ringdown portion,
        from merger_freq_calculation() in mergerfirstfuns.
        
    Returns
    -------
    [i_m_time,i_m_omega]: list of lists of floats
        The first list is the combined time values, the second list is the
        combined angular frequency values.
    """
    
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(time_offset) == float, 'time_offset should be a float.'
    assert type(i_time) == list, 'i_time should be a list.'
    assert type(i_omega) == list, 'i_omega should be a list.'
    assert type(m_time) == list, 'm_time should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    min_offset_m_time = np.empty((len(m_time)))
    for i in range(len(m_time)):                    #offsetting to match i_time
        min_offset_m_time[i] = m_time[i] + time_offset
        
    #now we stitch the inspiral and merger frequency waveforms together    
    i_m_omega = []                                  #combined omega
    i_m_time = []                                   #combined time
    for i in range(final_i_index):                  #inspiral segment
        i_m_omega.append(i_omega[i])
        i_m_time.append(i_time[i])
        
    for i in range(min_switch_ind,len(m_time)):     #merger segment
        i_m_omega.append(m_omega[i])
        i_m_time.append(min_offset_m_time[i])
        
    return [i_m_time,i_m_omega]

def frequency_SI_units(i_m_omega,M):
    """
    The angular frequency in geometric units translated to ordinary/temporal
    frequency in SI units (Hz). Useful for plotting and also required for the
    SNR calculator.
    
    Parameters
    ----------
    i_m_omega: list of floats
        Values of angular frequency over time for the entire duration of the
        gravitational waveform, from time_frequency_stitching().
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta() in
        inspiralfuns.
    
    Returns
    -------
    i_m_freq: list of floats
        Values of frequency in Hz for the entire duration of the gravitational
        waveform.
    """
    
    #input type checking
    assert type(i_m_omega) == list, 'i_m_omega should be a list.'
    assert type(M) == float, 'M should be a float.'
    
    pi=np.pi
    Msuns=4.923e-6                                  #geometric unit conversion
    
    i_m_freq = np.empty((len(i_m_omega)))
    for i in range(len(i_m_omega)):
        i_m_freq[i] = i_m_omega[i] / (M*Msuns*pi)
        
    #output type conversion
    i_m_freq = list(i_m_freq)
        
    return i_m_freq