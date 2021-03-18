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
        h_orth[i] = amp * np.cos(2*pi*freq*times[i] + pi*freq_dot*times**2 + \
                              init_phase)               #Shah eq. 1
        h_diag[i] = amp * np.sin(2*pi*freq*times[i] + pi*freq_dot*times**2 + \
                              init_phase)               #Shah eq. 1
            
    #output type conversion
    times = list(times)
    h_orth = list(h_orth)
    h_diag = list(h_diag)
    
    return [times,h_orth,h_diag]
