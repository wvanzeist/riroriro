#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The pre-matching parts of the procedure for simulating the merger/ringdown
portions of gravitational waves from binary black holes, collected into modular
functions.
"""

import numpy as np

def quasi_normal_modes(eta):
    """
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
    """
    
    #input type checking
    assert type(eta) == float, 'eta should be a float.'
    
    #basic sanity check for eta; should not be necessary if inspiralfuns is
    #working correctly but included just to be safe since this is the first
    #function in a new file
    assert eta <= 0.25 and eta > 0, ('eta should be positive and no larger '
                                     'than 0.25.')
    
    sfin = 2*np.sqrt(3)*eta - (390/79)*eta**2 + (2379/287)*eta**3 - \
        (4621/276)*eta**4               #final spin; based on Buskirk eq. 20
    sfin = float(sfin)                  #becomes numpy.float64 otherwise
    wqnm = 1 - 0.63*(1 - sfin)**0.3     #based on Buskirk eq. 19
    
    return (sfin,wqnm)

def gIRS_coefficients(eta,sfin):
    """
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
    """
    
    #input type checking
    assert type(eta) == float, 'eta should be a float.'
    assert type(sfin) == float, 'sfin should be a float.'
    
    #constants from Buskirk Appendix C
    Q = 2/((1 - sfin)**0.45)
    alpha = Q**-2 * (16313/562 + (21345/124)*eta)
    b = 16014/979 - (29132/1343)*eta**2
    C = 206/903 + (180/1141)*np.sqrt(eta) + (424/1205)*eta**2*(1/np.log(eta))
    kappa = 713/1056 - (23/193)*eta
    
    #output type conversion
    C = float(C)
    
    return (alpha,b,C,kappa)

def merger_freq_calculation(wqnm,b,C,kappa):
    """
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
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    #time in geometric units    
    time = np.empty((201))    #can't set it directly as range because of typing
    for i in range(201):
        time[i] = range(-100,101)[i]
        
    fhat = np.empty((201))
    m_omega = np.empty((201))
    for i in range(201):
        fhat[i] = (C/2) * (1 + 1/kappa)**(1 + kappa) * (1 - (1 + \
            (1/kappa)*np.exp(-2*time[i]*(1/b)))**-kappa)
        #based on Buskirk eq. 17
        m_omega[i] = 0.5 * wqnm * (1 - fhat[i])
        #based on Buskirk eq. 18, with 0.5 to standardise weird angle
        #definitions in Buskirk paper
        
    #output type conversion
    fhat = list(fhat)
    m_omega = list(m_omega)
        
    return [fhat,m_omega]

def fhat_differentiation(fhat):
    """
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
    """
    
    #input type checking
    assert type(fhat) == list, 'fhat should be a list.'
    
    #making sure fhat has the correct length, to avoid IndexErrors
    assert len(fhat) == 201, 'fhat should have length 201.'
    
    fhatdot = np.empty((201))
    #unusual averages at ends
    fhatdot[0] = fhat[1] - fhat[0]
    fhatdot[200] = fhat[200] - fhat[199]
    for i in range(1,200):              #derivative approximated as differences
        fhatdot[i] = 0.5*(fhat[i+1] - fhat[i-1])
        
    #output type conversion
    fhatdot = list(fhatdot)
    
    return fhatdot

def merger_time_conversion(M):
    """
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
    """
    
    #input type checking
    assert type(M) == float, 'M should be a float.'
    
    #geometric unit conversion
    Msuns=4.923e-6
    
    time = range(-100,101)                          #time in geometric units
    m_time = np.empty((201))                        #time in s (initialisation)
    for i in range(201):
        m_time[i] = time[i]*M*Msuns
        
    #output type conversion
    m_time = list(m_time)
    
    return m_time