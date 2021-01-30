************
matchingfuns
************

This is the documentation for the matchingfuns module, which consists of parts of the procedure for matching together the inspiral and merger/ringdown portions of the gravitational waveforms of binary black holes, collected into modular functions.

MQdiff
======

``MQdiff(i,i_time,i_omega,m_time,m_omega)``

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

min_switch_ind_finder
=====================

``min_switch_ind_finder(i_time,i_omega,m_time,m_omega)``

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
    
final_i_index_finder
====================

``final_i_index_finder(min_switch_ind,i_omega,m_omega)``

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

time_offset_finder
==================

``time_offset_finder(min_switch_ind,final_i_index,i_time,m_time)``

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
    
time_frequency_stitching
========================

``time_frequency_stitching(min_switch_ind,final_i_index,time_offset,i_time,i_omega,m_time,m_omega)``

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
    
frequency_SI_units
==================

``frequency_SI_units(i_m_omega,M)``

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
