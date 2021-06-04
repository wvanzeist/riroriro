**************
whitedwarffuns
**************

This is the documentation for the whitedwarffuns module, which consists of parts of the procedure for simulating the gravitational waveform of a white dwarf binary and evaluating LISA's response to such a waveform, collected into modular functions.

wd_amplitude
============

``wd_amplitude(Mc,freq,d)``

Calculates the instantaneous strain amplitude of a white dwarf binary,
based on Shah et al. (2012) equation 3/Królak et al. (2004) equation 17.

Parameters
----------
Mc: float
    The chirp mass of the binary, in solar masses.
freq: float
    The gravitational wave frequency of the binary, in Hz.
d: float
    The distance to the binary, in pc.
    
Returns
-------
amp: float
    The strain amplitude of the gravitational wave signal emitted by the
    binary (unitless).

characteristic_strain
=====================

``characteristic_strain(amp,freq,T_obs)``

Calculates the expected characteristic strain of a monochromatic
(non-chirping) white dwarf binary observed for a given amount of time,
based on Kupfer et al. (2018) equation 8/Moore et al. (2015) equation 41.

Parameters
----------
amp: float
    The strain amplitude of the gravitational wave signal emitted by the
    binary (unitless), from wd_amplitude().
freq: float
    The gravitational wave frequency of the binary, in Hz.
T_obs: float
    The length of time for which the binary is observed, in sec.
    
Returns
-------
hc: float
    The characteristic strain

wd_polarisations
================

``wd_polarisations(Mc,freq,d,T_sim,init_phase=0.0,sampling_freq=None,chirp=False)``

Calculates the values of the two polarisation of strain for a non-chirping
or weakly chirping white dwarf binary, based on Shah et al. (2012)
equations 1,2,4 and Królak et al. (2004) equation 21. Here we treat the
inclination as being optimal; effects of inclination are applied by a
subsequent function.

Parameters
----------
Mc: float
    The chirp mass of the binary, in solar masses.
freq: float
    The gravitational wave frequency of the binary (at t=0), in Hz.
d: float
    The distance to the binary, in pc.
T_sim: float
    The length of time to simulate the binary's waveform for, in sec.
init_phase: float
    The orbital phase of the binary at t=0, in rad. Default: 0.0.
sampling_freq: float or None
    The frequency at which measurements of the strain are to be taken. If
    not specified (set as None) the default sampling frequency is four
    times the Nyquist frequency.
chirp: bool
    If False, binary is simulated as monochromatic (non-chirping). If True,
    a weakly chirping binary is simulated via Shah et al. (2012) equation
    4. Default: False.
    
Returns
-------
[times,h_orth,h_diag]: list of lists of floats
    The first list is the times at which the strain has been calculated,
    the second is the values of the orthogonal/plus polarisation of strain
    at each timestep and the third is the diagonal/cross polarisation.
    
wd_inclination
==============

``wd_inclination(amp,**kwargs)``

Applies the effects of the binary's inclination of the strain amplitude of
a white dwarf binary, based on Królak et al. (2004) equation 16.

Parameters
----------
amp: float
    The strain amplitude of the gravitational wave signal emitted by the
    binary (unitless), from wd_amplitude().
kwargs: The inclination of the binary, expressed as either:
    iota: float
        The inclination angle, in radians.
    cosiota: float
        The cosine of the inclination angle.
        
Returns
-------
[amp_orth,amp_diag]: list of floats
    The first value is the adjusted strain amplitude for the
    orthogonal/plus polarisation, the second is for the diagonal/cross
    polarisation.
    
wd_inclination_polarisations
============================

``wd_inclination_polarisations(horth,hdiag,**kwargs)``

Similar to wd_inclination(), but instead of applying the inclination to a
single amplitude value, this applies it to the values of the two
polarisations of strain over time, based on Shah et al. (2012) equations
1,2 and Królak et al. (2004) equation 16.

Parameters
----------
horth: list of floats
    The orthogonal/plus polarisation of strain over time, from
    wd_polarisations().
hdiag: list of floats
    The diagonal/cross polarisation of strain over time, from
    wd_polarisations().
kwargs: The inclination of the binary, expressed as either:
    iota: float
        The inclination angle, in radians.
    cosiota: float
        The cosine of the inclination angle.

Returns
-------
[adj_horth,adj_hdiag]: list of lists of floats
    The first list is the adjusted orthogonal/plus values, the second is
    the adjusted diagonal/cross values.
    
lisa_rotation
=============

``lisa_rotation(times,eta_0=0.0,xi_0=0.0)``

Generates two time series of angles describing LISA's rotation over time,
including its rotation around the Sun and about its own axis, based on
Cutler (1998) equation 3.3 and Królak et al. (2004) pages 3 and 4.

Parameters
----------
times: list of floats
    The times at which strain has been calculated, from wd_polarisations().
eta_0: float
    The orbital phase of LISA around the Sun at t=0, in rad. Default: 0.0.
xi_0: float
    The orbital phase of LISA around its own axis at t=0, in rad. Default:
    0.0.
    
Returns
-------
[eta,xi]: list of lists of floats
    The first list is the orbital phases/angles of LISA around the Sun over
    time, the second is the same for LISA's rotation around its own axis.
    
wd_binary_vectors
=================

``wd_binary_vectors(theta_s,phi_s,iota,chi)``

Defines two vectors that describe properties of a binary in a stationary
frame of reference, which are used for calculating the polarisation angle
in lisa_angle_conversion(). They are defined here for speed, as these
calculations do not need to be in the loop of lisa_angle_conversion().

Parameters
----------
theta_s: float
    The ecliptic latitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary).
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
phi_s: float
    The ecliptic longitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary).
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
iota: float
    The inclination angle of the binary, in radians.
    This angle is invariant between the stationary and rotating frames of
    reference.
chi: float
    An angle that relates to the ascending node of the binary, in radians.
    It seems to be equivalent to the angle of the ascending node ± π/2.
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
    
Returns
-------
[L,P]: list of lists of floats
    The first list is the angular momentum vector, the second is the vector
    of the cross product of the line of sight and the angular momentum
    (giving the principal direction of orthogonal/plus polarisation).
    
lisa_angle_conversion
=====================

``lisa_angle_conversion(theta_s,phi_s,iota,eta,xi,L,P)``

Converts three angles (that are used for calculating LISA's response to a
signal) from a stationary frame of reference to one that follows LISA's
rotation, based on Cutler (1998) equation 3.4 and Królak et al. (2004)
equation 6 with zeta = -π/6 (for theta and phi) and a derivation from
Apostolatos et al. (1994) equation 5 (for psi).

Parameters
----------
theta_s: float
    The ecliptic latitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary).
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
phi_s: float
    The ecliptic longitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary).
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
iota: float
    The inclination angle of the binary, in radians.
    This angle is invariant between the stationary and rotating frames of
    reference.
eta: list of floats
    The orbital phases/angles of LISA around the Sun over time, from
    lisa_rotation().
    Not to be confused with the symmetric mass ratio, which also has the
    symbol eta.
xi: list of floats
    The orbital phases/angles of LISA around its own axis over time, from
    lisa_rotation().
L: list of floats
    The three-dimensional vector of the angular momentum of the binary,
    from wd_binary_vectors().
P: list of floats
    The three-dimensional vector of the principal direction of orthogonal/
    plus polarisation (equal to the vector cross product of the line of
    sight and the angular momentum), from wd_binary_vectors().
    
Returns
-------
[theta,phi,psi]: list of lists of floats
    The first list is the values of latitude with respect to LISA's
    rotating frame of reference at each timestep, the second is the
    corresponding longitudes and the third is the polarisation angles.
    
instantaneous_beam_pattern
==========================

``instantaneous_beam_pattern(theta_d,phi_d,psi_d)``

Outputs coefficients, known as the detector beam-pattern coefficients,
describing LISA's response to a binary's GW signal depending on the
binary's orientation/alignment as specified by three angles, based on
Cutler (1998) equation 3.12.
This function uses the alignment at a specific instant for the detector's
(rotating) frame of reference.

(These particular coefficients are also applicable to LIGO/Virgo and thus
also appear in the detectabilityfuns module.)

Parameters
----------
theta_d: float
    The relative latitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary). Ranges
    from 0 to π rad (180 deg).
    This angle should be with respect to the detector's (rotating) frame
    of reference.
phi_d: float
    The relative longitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary). Ranges
    from 0 to 2π rad (360 deg).
    This angle should be with respect to the detector's (rotating) frame
    of reference.
psi_d: float
    The polarisation angle of the binary. Ranges from 0 to π (180 deg).
    This angle should be with respect to the detector's (rotating) frame
    of reference.

Results
-------
[Fplus,Fcross]: list of floats
    The first quantity is the beam-pattern coefficient for
    orthogonal/plus polarisation, the second is for diagonal/cross
    polarisation.
    
lisa_beam_pattern
=================

``lisa_beam_pattern(theta,phi,psi)``

Evaluates instantaneous_beam_pattern() for each timestep in the white dwarf
binary simulation, giving the beam pattern coefficients over time.

Parameters
----------
theta: list of floats
    The source's latitude with respect to LISA's rotating frame of
    reference at each timestep, from lisa_angle_conversion().
phi: list of floats
    The source's longitude with respect to LISA's rotating frame of
    reference at each timestep, from lisa_angle_conversion().
psi: list of floats
    The source's polarisation angle in LISA's rotating frame of reference
    at each timestep, from lisa_angle_conversion().
    
Returns
-------
[Fplus_t,Fcross_t]: list of lists of floats
    The first list is the beam-pattern coefficients for orthogonal/plus
    polarisation at each time step, the second is the diagonal/cross
    polarisation coefficients.
    
lisa_phase_modulation
=====================

``lisa_phase_modulation(Fplus_t,Fcross_t,amp_orth,amp_diag)``

Calculates coefficients describing the phase modulation that affects LISA's
response to an incoming GW signal, based on Cutler (1998) equation 3.15b/
Cornish et al. (2003) equation 5.

Parameters
----------
Fplus_t: list of floats
    The detector beam-pattern coefficients for orthogonal/plus polarisation
    over time, from lisa_beam_pattern().
Fcross_t: list of floats
    The detector beam-pattern coefficients for diagonal/cross polarisation
    over time, from lisa_beam_pattern().
amp_orth: float
    The inclination-adjusted orthogonal/plus polarisation of the strain
    amplitude, from wd_inclination().
amp_diag: float
    The inclination-adjusted diagonal/cross polarisation of the strain
    amplitude, from wd_inclination().
    
Returns
-------
varphi_p: list of floats
    The phase modulation coefficients at each timestep.
    
lisa_frequency_modulation
=========================

``lisa_frequency_modulation(times,freq,theta_s,phi_s)``

Calculates coefficients describing the frequency (Doppler) modulation that
affects LISA's response to an incoming GW signal, based on Cutler (1998)
equation 3.15c/Cornish et al. (2003) equation 4.

Parameters
----------
times: list of floats
    The times at which strain has been calculated, from wd_polarisations().
freq: float
    The gravitational wave frequency of the binary (at t=0), in Hz.
theta_s: float
    The ecliptic latitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary).
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
phi_s: float
    The ecliptic longitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary).
    This angle should be given with respect to a stationary frame of
    reference (celestial reference).
    
Returns
-------
varphi_d: list of floats
    The frequency (Doppler) modulation coefficients at each timestep.
    
lisa_amplitude_modulation
=========================

``lisa_amplitude_modulation(Fplus_t,Fcross_t,amp_orth,amp_diag)``

Applies the amplitude modulation that affects LISA's response to the strain
amplitude of an incoming GW signal, based on Cutler (1998) equation 3.15a/
Cornish et al. (2003) equation 3.

Parameters
----------
Fplus_t: list of floats
    The detector beam-pattern coefficients for orthogonal/plus polarisation
    over time, from lisa_beam_pattern().
Fcross_t: list of floats
    The detector beam-pattern coefficients for diagonal/cross polarisation
    over time, from lisa_beam_pattern().
amp_orth: float
    The inclination-adjusted orthogonal/plus polarisation of the strain
    amplitude, from wd_inclination().
amp_diag: float
    The inclination-adjusted diagonal/cross polarisation of the strain
    amplitude, from wd_inclination().

Returns
-------
A_mod: list of floats
    The strain amplitude of the gravitational wave signal with modulation
    applied, at each timestep.
    
lisa_detector_response
======================

``lisa_detector_response(times,A_mod,varphi_d,varphi_p,freq,init_phase=0.0)``

For an incoming GW signal, calculates LISA's detector response (the strain
observed in the detector) to that signal, based on Cornish et al. (2003)
equations 1,2 (also Shah et al. (2012) equations 6,8, though the form given
there does not correctly match that of Cornish).
This form is valid if the change of frequency of the binary over the
duration of the observation is much smaller than the frequency itself;
otherwise, the more complicated (and time-consuming to evaluate) Cutler
(1998) equation 3.14 would need to be used.

Parameters
----------
times: list of floats
    The times at which strain has been calculated, from wd_polarisations().
A_mod: list of floats
    The modulation-adjusted strain amplitude of the gravitational wave
    signal over time, from lisa_amplitude_modulation().
varphi_d: list of floats
    The frequency (Doppler) modulation coefficients at each timestep, from
    lisa_frequency_modulation().
varphi_p: list of floats
    The phase modulation coefficients at each timestep, from
    lisa_phase_modulation().
freq: float
    The gravitational wave frequency of the binary (at t=0), in Hz.
init_phase: float
    The phase angle of the signal at t=0. Default: 0.0.
    
Returns
-------
A_lisa: list of floats
    The strain amplitude of the GW signal over time as it is detected by
    LISA, account for the various time-varying factors affecting LISA's
    response.
