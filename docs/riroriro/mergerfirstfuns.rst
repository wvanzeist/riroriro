***************
mergerfirstfuns
***************

This is the documentation for the mergerfirstfuns module, which consists of the pre-matching parts of the procedure for simulating the merger/ringdown portions of gravitational waves from binary black holes, collected into modular functions.

quasi_normal_modes
==================

``quasi_normal_modes(eta)``

Calculation of the final spin and quasi-normal mode factor used in the
calculation of angular frequency for the merger/ringdown waveform, based on
Buskirk et al. (2019) equations 19 and 20.

Parameters
----------
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta() in inspiralfuns.

Returns
-------
(sfin,wqnm): tuple of floats
    The first constant is the final spin, the second is the quasi-normal
    mode factor used in subsequent calculations.

gIRS_coefficients
=================

``gIRS_coefficients(eta,sfin)``

Calculation of several gIRS (generic implicit rotating source)-related
coefficients used in the calculation of angular frequency for the merger/
ringdown waveform, based on Buskirk et al. (2019) Appendix C.

Parameters
----------
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta() in inspiralfuns.
sfin: float
    Final spin value, from quasi_normal_modes().
    
Returns
-------
(alpha,b,C,kappa): tuple of floats
    Four gIRS-related constants used in subsequent calculations.
    (NOTE: alpha is not used by anything in mergerfirstfuns but *is* used
    in mergersecondfuns.)
    
merger_freq_calculation
=======================

``merger_freq_calculation(wqnm,b,C,kappa)``

Calculation of orbital angular frequency for the merger/ringdown portion,
based on Buskirk et al. (2019) equations 17 and 18.

Parameters
----------
wqnm: float
    Quasi-normal mode factor, from quasi_normal_modes().
b: float
    A gIRS coefficient, from gIRS_coefficients().
C: float
    A gIRS coefficient, from gIRS_coefficients().
kappa: float
    A gIRS coefficient, from gIRS_coefficients().

Returns
-------
[fhat,m_omega]: list of lists of floats
    First list is the values over time of a sort of frequency parameter
    called fhat (f^), second list is the angular frequency.
    
fhat_differentiation
====================

``fhat_differentiation(fhat)``

Calculation of derivative of fhat used by amplitude calculations in
mergersecondfuns.

Parameters
----------
fhat: list of floats
    Values of a sort of frequency parameter called fhat (f^) over time,
    from merger_freq_calculation().
    
Returns
-------
fhatdot: list of floats
    Values of the time-derivative of fhat over time.
    
merger_time_conversion
======================

``merger_time_conversion(M)``

Calculating times in real units corresponding to the times in geometric
units used by other merger/ringdown functions.

Parameters
----------
M: float
    Total mass of the binary, can be obtained from get_M_and_eta() in
    inspiralfuns.
    
Returns
-------
m_time: list of floats
    The list of timesteps used by other merger/ringdown functions, but in
    seconds instead of geometric units.
