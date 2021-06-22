****************
mergersecondfuns
****************

This is the documentation for the mergersecondfuns module, which consists of the post-matching parts of the procedure for simulating the merger/ringdown portions of gravitational waves from binary black holes, collected into modular functions.

merger_phase_calculation
========================

``merger_phase_calculation(min_switch_ind,final_i_index,i_phase,m_omega)``

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

phase_stitching
===============

``phase_stitching(final_i_index,i_phase,m_phase)``

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

merger_strain_amplitude
=======================

``merger_strain_amplitude(min_switch_ind,final_i_index,alpha,i_amp,m_omega,fhat,fhatdot)``

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

amplitude_stitching
===================

``amplitude_stitching(final_i_index,i_amp,m_amp)``

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

merger_polarisations
====================

``merger_polarisations(final_i_index,m_amp,m_phase,i_Aorth)``

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

polarisation_stitching
======================

``polarisation_stitching(final_i_index,i_Aorth,i_Adiag,m_Aorth,m_Adiag)``

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
