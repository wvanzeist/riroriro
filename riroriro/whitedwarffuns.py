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
    white dwarf binary, based on Shah et al. (2012) equations 1,2.

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
        assert 0 <= kwargs['iota'] <= np.pi/2, ('iota should be between 0 and '
                                                'π/2 rad.')
        cosiota = np.cos(kwargs['iota'])
    
    if 'cosiota' in kwargs:
        #checking cosiota is within the expected range
        assert 0 <= kwargs['cosiota'] <= 1, ('cosiota should be between 0 and '
                                             '1.')
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
        The ecliptic latitude, one of the angles describing the direction of
        the line of sight to the gravitational wave source relative to the axes
        of the detector’s arms (sky-location coordinates of the binary). Ranges
        from 0 to π/2 rad (90 deg).
        This angle should be with respect to the detector's (rotating) frame
        of reference.
    phi_d: float
        The ecliptic longitude, one of the angles describing the direction of
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
