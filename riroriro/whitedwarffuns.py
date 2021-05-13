#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for evaluating the gravitational waves emitted by non-chirping or
weakly chirping white dwarf binaries (for LISA).
"""

import numpy as np

def wd_amplitude(Mc,freq,d):
    """
    Calculates the instantaneous strain amplitude of a white dwarf binary,
    based on Shah et al. (2012) equation 3/Królak et al. (2004) equation 17.
    
    Parameters
    ----------
    Mc: float
        The chirp mass of the binary, in solar masses.
    freq: float
        The orbital frequency of the binary, in Hz.
    d: float
        The distance to the binary, in pc.
        
    Returns
    -------
    amp: float
        The strain amplitude of the gravitational wave signal emitted by the
        binary (unitless).
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
        
    #basic constants
    pi=np.pi
    c=2.99792458e8 #m s^-1
    G=6.674e-11 #m^3 kg^-1 s^-2
    
    #factors for conversion to SI units for convenience
    Msun_to_kg=1.9885e30
    pc_to_m=3.0857e16
    
    amp = (4*(G*Mc*Msun_to_kg)**(5/3) * (pi*freq)**(2/3)) / (c**4 * d*pc_to_m)
    
    return amp

def characteristic_strain(amp,freq,T_obs):
    """
    Calculates the expected characteristic strain of a monochromatic
    (non-chirping) white dwarf binary observed for a given amount of time,
    based on Kupfer et al. (2018) equation 8/Moore et al. (2015) equation 41.
    
    Parameters
    ----------
    amp: float
        The strain amplitude of the gravitational wave signal emitted by the
        binary (unitless), from wd_amplitude().
    freq: float
        The orbital frequency of the binary, in Hz.
    T_obs: float
        The length of time for which the binary is observed, in sec.
        
    Returns
    -------
    hc: float
        The characteristic strain
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    hc = np.sqrt(freq * T_obs) * amp
    
    return hc

def wd_polarisations(Mc,freq,d,T_sim,init_phase=0.0,chirp=False):
    """
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
        The orbital frequency of the binary (at t=0), in Hz.
    d: float
        The distance to the binary, in pc.
    T_sim: float
        The length of time to simulate the binary's waveform for, in sec.
    init_phase: float
        The orbital phase of the binary at t=0, in rad. Default: 0.0.
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
    """
    
    #input type checking
    assert type(Mc) == float, 'Mc should be a float.'
    assert type(freq) == float, 'freq should be a float.'
    assert type(d) == float, 'd should be a float.'
    assert type(T_sim) == float, 'T_sim should be a float.'
    assert type(init_phase) == float, 'init_phase should be a float.'
    assert type(chirp) == bool, 'chirp should be a bool.'
    
    #basic constants
    pi=np.pi
    c=2.99792458e8 #m s^-1
    G=6.674e-11 #m^3 kg^-1 s^-2
    
    #factors for conversion to SI units for convenience
    Msun_to_kg=1.9885e30
    
    amp = wd_amplitude(Mc,freq,d)                       #Shah eq. 3
    
    if chirp == True:
        freq_dot = (96/5) * (freq/(Mc*Msun_to_kg)) * (pi*freq * \
            Mc*Msun_to_kg)**(8/3) * (G/c**3)**(5/3)
        #Shah eq. 4 with dimensional correction term based on Królak eq. 21
    else:
        freq_dot = 0
        
    #the freq_dot equation does not need to be integrated to obtain freq as
    #long as freq >> freq_dot
    
    no_of_samples = int(np.ceil(T_sim * 8*freq))
    #sampling at t=0 Nyquist frequency times 4
    times = np.linspace(0,T_sim,no_of_samples) #times to calculate strain at
    
    h_orth = np.empty((no_of_samples))          #orthogonal/plus polarisation
    h_diag = np.empty((no_of_samples))          #diagonal/cross polarisation
    
    for i in range(no_of_samples):
        h_orth[i] = amp * np.cos(2*pi*freq*times[i] + pi*freq_dot*times[i]**2 \
                              + init_phase)             #Shah eq. 1
        h_diag[i] = amp * np.sin(2*pi*freq*times[i] + pi*freq_dot*times[i]**2 \
                              + init_phase)             #Shah eq. 2
            
    #output type conversion
    times = list(times)
    h_orth = list(h_orth)
    h_diag = list(h_diag)
    
    return [times,h_orth,h_diag]

def wd_inclination(horth,hdiag,**kwargs):
    """
    Applies the effects of the binary's inclination on the strain values of a
    white dwarf binary, based on Shah et al. (2012) equations 1,2 and Królak et
    al. (2004) equation 16.

    Parameters
    ----------
    horth: list of floats
        The orthogonal/plus polarisation of strain over time, from
        wd_polarisations().
    hdiag: list of floats
        The diagonal/cross polarisation of strain over time, from
        wd_polarisations().
    **kwargs: The inclination of the binary, expressed as either:
        iota: float
            The inclination angle, in radians.
        cosiota: float
            The cosine of the inclination angle.

    Returns
    -------
    [adj_horth,adj_hdiag]: list of lists of floats
        The first list is the adjusted orthogonal/plus values, the second is
        the adjusted diagonal/cross values.
    """
    
    #input type checking
    assert type(horth) == list, 'horth should be a list.'
    assert type(hdiag) == list, 'hdiag should be a list.'
    for each_variable in kwargs.values():
        assert type(each_variable) == float, ('The inclination should be a '
                                              'float.')
    
    #making sure the user specifies at least one of iota and cosiota
    if not ('iota' in kwargs or 'cosiota' in kwargs):
        raise TypeError('Please specify iota or cosiota.')
        
    #making sure the user doesn't specify both iota and cosiota
    if 'iota' in kwargs and 'cosiota' in kwargs:
        raise TypeError('Please specify either iota or cosiota, not both.')
    
    if 'iota' in kwargs:
        #checking iota is within the expected range
        assert 0 <= kwargs['iota'] <= np.pi, ('iota should be between 0 and π '
                                              'rad.')
        cosiota = np.cos(kwargs['iota'])
    elif 'cosiota' in kwargs:
        #checking cosiota is within the expected range
        assert -1 <= kwargs['cosiota'] <= 1, ('cosiota should be between -1 '
                                              'and 1.')
        cosiota = kwargs['cosiota']
        
    #applying inclination to strain
    adj_horth = np.empty((len(horth)))
    adj_hdiag = np.empty((len(hdiag)))
    for i in range(len(horth)):
        adj_horth[i] = 0.5*(1 + cosiota**2) * horth[i]  #part of Shah eq. 1
        adj_hdiag[i] = cosiota * hdiag[i]               #part of Shah eq. 2
    
    #output type conversion
    adj_horth = list(adj_horth)
    adj_hdiag = list(adj_hdiag)
    
    return [adj_horth,adj_hdiag]

def lisa_rotation(times,eta_0=0.0,xi_0=0.0):
    """
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
    """
    
    #input type checking
    assert type(times) == list, 'times should be a list.'
    assert type(eta_0) == float, 'eta_0 should be a float.'
    assert type(xi_0) == float, 'xi_0 should be a float.'
    
    #basic constants
    pi=np.pi
    year_to_sec=3.154e7
    
    eta = np.empty((len(times)))
    xi = np.empty((len(times)))
    for i in range(len(times)):
        eta[i] = 2*pi*times[i]/year_to_sec + eta_0
        xi[i] = -2*pi*times[i]/year_to_sec + xi_0           #note sign
    
    #output type conversion
    eta = list(eta)
    xi = list(xi)
    
    return [eta,xi]

def wd_binary_vectors(theta_s,phi_s,iota,chi):
    """
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
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    #components of L
    Lx = (np.cos(chi) + np.cos(phi_s)**2*np.sin(theta_s)**2*(1 - np.cos(chi)))\
        *np.cos(phi_s)*np.sin(theta_s + iota) + (0.5*np.sin(2*phi_s)* \
        np.sin(theta_s)**2*(1 - np.cos(chi)) - np.cos(theta_s)*np.sin(chi)) \
        *np.sin(phi_s)*np.sin(theta_s + iota) + (0.5*np.cos(phi_s)* \
        np.sin(2*theta_s)*(1 - np.cos(chi)) + np.sin(phi_s)*np.sin(theta_s)* \
        np.sin(chi))*np.cos(theta_s + iota)
    Ly = (0.5*np.sin(2*phi_s)*np.sin(theta_s)**2*(1 - np.cos(chi)) + \
        np.cos(theta_s)*np.sin(chi))*np.cos(phi_s)*np.sin(theta_s + iota) + \
        (np.cos(chi) + np.sin(phi_s)**2*np.sin(theta_s)**2*(1 - np.cos(chi))) \
        *np.sin(phi_s)*np.sin(theta_s + iota) + (0.5*np.sin(phi_s)* \
        np.sin(2*theta_s)*(1 - np.cos(chi)) - np.cos(phi_s)*np.sin(theta_s)* \
        np.sin(chi))*np.cos(theta_s + iota)
    Lz = (0.5*np.cos(phi_s)*np.sin(2*theta_s)*(1 - np.cos(chi)) - \
        np.sin(phi_s)*np.sin(theta_s)*np.sin(chi))*np.cos(phi_s)* \
        np.sin(theta_s + iota) + (0.5*np.sin(phi_s)*np.sin(2*theta_s)*(1 - \
        np.cos(chi)) + np.cos(phi_s)*np.sin(theta_s)*np.sin(chi)) \
        *np.sin(phi_s)*np.sin(theta_s + iota) + (np.cos(chi) + \
        np.cos(theta_s)**2*(1 - np.cos(chi)))*np.cos(theta_s + iota)
    
    #components of P
    Px = np.sin(phi_s)*np.sin(theta_s)*Lz - np.cos(theta_s)*Ly
    Py = np.cos(theta_s)*Lx - np.cos(phi_s)*np.sin(theta_s)*Lz
    Pz = np.cos(phi_s)*np.sin(theta_s)*Ly - np.sin(phi_s)*np.sin(theta_s)*Lx
    
    #collecting components into lists representing vectors
    L = [Lx,Ly,Lz]
    P = [Px,Py,Pz]
    
    return [L,P]

def lisa_angle_conversion(theta_s,phi_s,iota,eta,xi,L,P):
    """
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
    """
    
    #input type checking
    assert type(theta_s) == float, 'theta_s should be a float.'
    assert type(phi_s) == float, 'phi_s should be a float.'
    assert type(iota) == float, 'iota should be a float.'
    assert type(eta) == list, 'eta should be a list.'
    assert type(xi) == list, 'xi should be a list.'
    assert type(L) == list, 'L should be a list.'
    assert type(P) == list, 'P should be a list.'
    
    theta = np.empty((len(eta)))
    phi = np.empty((len(eta)))
    psi = np.empty((len(eta)))
    
    for i in range(len(eta)):
        theta[i] = np.arccos(0.5*(-np.sqrt(3)*np.cos(eta[i])*np.cos(phi_s)* \
            np.sin(theta_s) - np.sqrt(3)*np.sin(eta[i])*np.sin(phi_s)* \
            np.sin(theta_s) + np.cos(theta_s)))
        phi[i] = np.arctan2((-np.sin(eta[i])*np.sin(xi[i]) + \
            0.5*np.cos(eta[i])*np.cos(xi[i]))*np.cos(phi_s)*np.sin(theta_s) + \
            (np.cos(eta[i])*np.sin(xi[i]) + 0.5*np.sin(eta[i])*np.cos(xi[i])) \
            *np.sin(phi_s)*np.sin(theta_s) + 0.5*np.sqrt(3)*np.cos(xi[i])* \
            np.cos(theta_s), \
            (np.sin(eta[i])*np.cos(xi[i]) + 0.5*np.cos(eta[i])*np.sin(xi[i])) \
            *np.cos(phi_s)*np.sin(theta_s) + (-np.cos(eta[i])*np.cos(xi[i]) + \
            0.5*np.sin(eta[i])*np.sin(xi[i]))*np.sin(phi_s)*np.sin(theta_s) + \
            0.5*np.sqrt(3)*np.sin(xi[i])*np.cos(theta_s)) #"y", then "x"
        #check handedness
        psi[i] = np.arctan2(0.5*(-np.sqrt(3)*np.cos(eta[i])*L[0] - np.sqrt(3) \
            *np.sin(eta[i])*L[1] + L[2]) - np.cos(iota)*np.cos(theta[i]), \
            0.5*(-np.sqrt(3)*np.cos(eta[i])*P[0] - np.sqrt(3)*np.sin(eta[i]) \
            *P[1] + P[2]))
    
    #output type conversion
    theta = list(theta)
    phi = list(phi)
    psi = list(psi)
    
    return [theta,phi,psi]

def instantaneous_beam_pattern(theta_d,phi_d,psi_d):
    """
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
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    Fplus = 0.5*(1 + np.cos(theta_d)**2)*np.cos(2*phi_d)*np.cos(2*psi_d) - \
        np.cos(theta_d)*np.sin(2*phi_d)*np.sin(2*psi_d) #Cutler eq. 3.12a
    Fcross = 0.5*(1 + np.cos(theta_d)**2)*np.cos(2*phi_d)*np.sin(2*psi_d) + \
        np.cos(theta_d)*np.sin(2*phi_d)*np.cos(2*psi_d) #Cutler eq. 3.12b
    
    #output type conversion
    Fplus = float(Fplus)
    Fcross = float(Fcross)
    
    return [Fplus,Fcross]

def lisa_beam_pattern(theta,phi,psi):
    """
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
        The source's polarisation angle in LISA's rotating frame of reference,
        from lisa_angle_conversion().
        
    Returns
    -------
    [Fplus_t,Fcross_t]: list of lists of floats
        The first list is the beam-pattern coefficients for orthogonal/plus
        polarisation at each time step, the second is the diagonal/cross
        polarisation coefficients.
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'
        
    Fplus_t = np.empty((len(theta)))
    Fcross_t = np.empty((len(theta)))
    
    #instantaneous_beam_pattern() for each timestep
    for i in range(len(theta)):
        Fplus_t[i], Fcross_t[i] = instantaneous_beam_pattern(theta[i],phi[i], \
                                                             psi[i])
    
    #output type conversion
    Fplus_t = list(Fplus_t)
    Fcross_t = list(Fcross_t)
    
    return [Fplus_t,Fcross_t]

def lisa_phase_modulation(adj_horth,adj_hdiag,Fplus_t,Fcross_t):
    """
    Calculates coefficients describing the phase modulation that affects LISA's
    response to an incoming GW signal, based on Cutler (1998) equation 3.15b/
    Cornish et al. (2003) equation 5.
    
    Parameters
    ----------
    adj_horth: list of floats
        The inclination-adjusted orthogonal/plus polarisation of strain over
        time, from wd_inclination().
    adj_hdiag: list of floats
        The inclination-adjusted diagonal/cross polarisation of strain over
        time, from wd_inclination().
    Fplus_t: list of floats
        The detector beam-pattern coefficients for orthogonal/plus polarisation
        over time, from lisa_beam_pattern().
    Fcross_t: list of floats
        The detector beam-pattern coefficients for diagonal/cross polarisation
        over time, from lisa_beam_pattern().
        
    Returns
    -------
    varphi_p: list of floats
        The phase modulation coefficients at each timestep.
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'
        
    varphi_p = np.empty((len(adj_horth)))
    for i in range(len(varphi_p)):
        varphi_p[i] = -np.arctan2(adj_hdiag[i]*Fcross_t[i], \
                                  adj_horth[i]*Fplus_t[i])
            
    #output type conversion
    varphi_p = list(varphi_p)
    
    return varphi_p

def lisa_frequency_modulation(times,freq,theta_s,phi_s):
    """
    Calculates coefficients describing the frequency (Doppler) modulation that
    affects LISA's response to an incoming GW signal, based on Cutler (1998)
    equation 3.15c/Cornish et al. (2003) equation 4.
    
    Parameters
    ----------
    times: list of floats
        The times at which strain has been calculated, from wd_polarisations().
    freq: float
        The orbital frequency of the binary (at t=0), in Hz.
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
    """
    
    #input type checking
    assert type(times) == list, 'times should be a list.'
    assert type(freq) == float, 'freq should be a float.'
    assert type(theta_s) == float, 'theta_s should be a float.'
    assert type(phi_s) == float, 'phi_s should be a float.'
    
    #basic constants
    pi=np.pi
    c=2.99792458e8 #m s^-1
    year_to_sec=3.154e7
    AU_to_m=1.495978707e11
    
    #Cornish eq. 4 calculation
    varphi_d = np.empty((len(times)))
    for i in range(len(times)):
        varphi_d[i] = 2*pi*freq*(AU_to_m/c)*np.sin(theta_s)*np.cos(2*pi* \
            times[i]/year_to_sec - phi_s)
    
    #output type conversion
    varphi_d = list(varphi_d)
    
    return varphi_d

def lisa_amplitude_modulation(adj_horth,adj_hdiag,Fplus_t,Fcross_t):
    """
    Applies the amplitude modulation that affects LISA's response to the strain
    amplitude of an incoming GW signal, based on Cutler (1998) equation 3.15a/
    Cornish et al. (2003) equation 3.
    
    Parameters
    ----------
    adj_horth: list of floats
        The inclination-adjusted orthogonal/plus polarisation of strain over
        time, from wd_inclination().
    adj_hdiag: list of floats
        The inclination-adjusted diagonal/cross polarisation of strain over
        time, from wd_inclination().
    Fplus_t: list of floats
        The detector beam-pattern coefficients for orthogonal/plus polarisation
        over time, from lisa_beam_pattern().
    Fcross_t: list of floats
        The detector beam-pattern coefficients for diagonal/cross polarisation
        over time, from lisa_beam_pattern().
    
    Returns
    -------
    A_mod: list of floats
        The strain amplitude of the gravitational wave signal with modulation
        applied, at each timestep.
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'
        
    A_mod = np.empty((len(adj_horth)))
    for i in range(len(A_mod)):
        A_mod[i] = np.sqrt((adj_horth[i]*Fplus_t[i])**2 + (adj_hdiag[i]* \
                                                           Fcross_t[i])**2)
    
    #output type conversion
    A_mod = list(A_mod)
    
    return A_mod

def lisa_detector_response(times,A_mod,varphi_d,varphi_p,freq,init_phase=0.0):
    """
    For an incoming GW signal, calculates LISA's detector response (the strain
    observed in the detector) to that signal, based on Cornish et al. (2003)
    equations 1,2/Shah et al. (2012) equations 6,8.
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
        The orbital frequency of the binary (at t=0), in Hz.
    init_phase: float
        The phase angle of the signal at t=0. Default: 0.0.
        
    Returns
    -------
    A_lisa: list of floats
        The strain amplitude of the GW signal over time as it is detected by
        LISA, account for the various time-varying factors affecting LISA's
        response.
    """
    
    #input type checking
    assert type(times) == list, 'times should be a list.'
    assert type(A_mod) == list, 'A_mod should be a list.'
    assert type(varphi_d) == list, 'varphi_d should be a list.'
    assert type(varphi_p) == list, 'varphi_p should be a list.'
    assert type(freq) == float, 'freq should be a float.'
    assert type(init_phase) == float, 'init_phase should be a float.'
    
    A_lisa = np.empty((len(times)))
    for i in range(len(times)):
        A_lisa[i] = A_mod[i]*np.cos(2*np.pi*freq*times[i] + init_phase + \
                                    varphi_d[i] + varphi_p[i])
            
    #output type conversion
    A_lisa = list(A_lisa)
    
    return A_lisa
