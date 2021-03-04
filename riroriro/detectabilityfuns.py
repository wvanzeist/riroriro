#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parts of the procedure for calculating the detectability fraction of a merger
given its optimal-alignment SNR, collected into modular functions.
"""

import numpy as np
from random import random
from scipy.interpolate import interp1d

def cdf_generator(N=10**6):
    """
    Generates the cumulative distribution function (CDF) of the projection
    function Theta, for use with detectability_calculator(), based on Finn
    (1996), Belczynski et al. (2013), Belczynski et al. (2014).
    
    Parameters
    ----------
    N: int
        The number of random samples of Theta you want to take to build the
        CDF. Default: 10**6.
        
    Returns
    -------
    Theta_CDF: function
        The CDF of the projection function Theta.
    min_CDF: float
        The lower boundary of the range over which Theta_CDF is defined.
    max_CDF: float
        The upper boundary of the range over which Theta_CDF is defined.
    """
    
    #input type checking
    assert type(N) == int, 'N should be an int.'
    
    pi=np.pi
    
    Theta_dist=np.empty((N))
    for i in range(N):
        psi=random()*pi
        phi=random()*2*pi
        costh=random()
        cosiota=random()
        #sinth = np.sqrt(1 - costh**2)
        #siniota = np.sqrt(1 - sinth**2)
        Fplus = 0.5*(1 + costh**2)*np.cos(2*phi)*np.cos(2*psi) - \
            costh*np.sin(2*phi)*np.sin(2*psi)
        Fcross = 0.5*(1 + costh**2)*np.cos(2*phi)*np.sin(2*psi) + \
            costh*np.sin(2*phi)*np.cos(2*psi)
        Theta = 0.5*np.sqrt(Fplus**2*(1 + cosiota**2)**2 + \
            4*Fcross**2*cosiota**2)
        Theta_dist[i] = Theta
        
    Theta_dist = np.sort(Theta_dist)                    #ordering
    Theta_CDF = interp1d(Theta_dist, np.linspace(0,1,N), kind='cubic')
    #this takes a value of Theta and returns the proportion of values of Theta
    #lower than that value of Theta
    
    #limits of the range over which Theta_CDF is defined
    min_CDF = Theta_dist[0]
    max_CDF = Theta_dist[N-1]
    
    #output type conversion
    min_CDF = float(min_CDF)
    max_CDF = float(max_CDF)
    
    return Theta_CDF, min_CDF, max_CDF

def detectability_calculator(Theta_CDF,min_CDF,max_CDF,SNR_in):
    """
    Given the optimal-alignment SNR of a merger, this function returns the
    fraction of arbitrary orientations in which the merger would be expected to
    be observable (i.e. have a SNR above 8).
    
    Parameters
    ----------
    Theta_CDF: function
        The CDF of the projection function Theta, from cdf_generator().
    min_CDF: float
        The lower boundary of the range over which Theta_CDF is defined, from
        cdf_generator().
    max_CDF: float
        The upper boundary of the range over which Theta_CDF is defined, from
        cdf_generator().
    SNR_in: float
        The optimal-alignment SNR of the merger in question, can be obtained
        from snrcalculatorfuns.
        
    Returns
    -------
    det: float
        The detectability fraction of the merger.
    """
    
    #input type checking
    assert type(min_CDF) == float, 'min_CDF should be a float.'
    assert type(max_CDF) == float, 'max_CDF should be a float.'
    assert type(SNR_in) == float, 'SNR_in should be a float.'
    
    if SNR_in <= (8 / max_CDF):
        det = 0                             #always undetectable
        #set separately so Theta_CDF doesn't get input above 1
        #the max_CDF (just below 1) has to be included because the function
        #Theta_dist is not defined for the small range between max_CDF and 1
        #and this can occasionally cause problems if just 8 is used
    elif SNR_in >= (8 / min_CDF):
        det = 1                             #always detectable
        #similarly, Theta_CDF is not defined between 0 and min_CDF
    else:
        threshold = 8 / SNR_in
        #how much the SNR would need to be lowered to hit the threshold SNR=8
        undet = Theta_CDF(threshold)
        #the proportion of Theta values that would lower SNR below 8
        det = 1 - undet             #proportion that would keep SNR above 8
    
    #output type conversion
    det = float(det)
    
    return det

def specific_orientation_SNR(theta,phi,iota,psi,SNR_in,angle_unit='rad'):
    """
    Given the optimal-alignment SNR of a merger, this function returns the SNR
    that would result if the detector and binary had a specific orientation/
    alignment, specified by four angles as in Finn (1996), Belczynski et al.
    (2013), Belczynski et al. (2014).
    
    Parameters
    ----------
    theta: float
        One of the angles descriving the direction of the line of sight to the
        gravitational wave source relative to the axes of the detector’s arms
        (sky-location coordinates). Ranges from 0 to π/2 rad (90 deg).
    phi: float
        One of the angles descriving the direction of the line of sight to the
        gravitational wave source relative to the axes of the detector’s arms
        (sky-location coordinates). Ranges from 0 to 2π rad (360 deg).
    iota: float
        The inclination angle of the binary. Ranges from 0 to π/2 rad (90 deg).
    psi: float
        The polarisation angle of the binary. Ranges from 0 to π (180 deg).
    SNR_in: float
        The optimal-alignment SNR of the merger in question, can be obtained
        from snrcalculatorfuns.
    angle_unit: str
        Specifies whether the input angles are given in 'rad' or 'deg'; the
        default is 'rad'.
    
    Returns
    -------
    SNR_out: float
        The SNR of the merger in question at the specific orientation given by
        the input angles.
    """
    
    #input type checking
    assert type(theta) == float, 'theta should be a float.'
    assert type(phi) == float, 'phi should be a float.'
    assert type(iota) == float, 'iota should be a float.'
    assert type(psi) == float, 'psi should be a float.'
    assert type(SNR_in) == float, 'SNR_in should be a float.'
    
    pi=np.pi
    
    #converting from deg to rad, if necessary
    if angle_unit == 'deg':
        theta *= 2*pi/360
        phi *= 2*pi/360
        iota *= 2*pi/360
        psi *= 2*pi/360
    elif angle_unit == 'rad':
        pass
    else:
        raise ValueError('angle_unit must be either \'rad\' or \'deg\'.')
    
    #checking input angles are within the expected ranges
    assert 0 <= theta <= pi/2, ('theta should be between 0 and π/2 rad (90 '
                                'deg).')
    assert 0 <= phi <= 2*pi, 'phi should be between 0 and 2π rad (360 deg).'
    assert 0 <= iota <= pi/2, 'iota should be between 0 and π/2 rad (90 deg).'
    assert 0 <= psi <= pi, 'psi should be between 0 and π rad (180 deg).'
    
    #calculating projection function
    Fplus = 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.cos(2*psi) - \
        np.cos(theta)*np.sin(2*phi)*np.sin(2*psi)
    Fcross = 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.sin(2*psi) + \
        np.cos(theta)*np.sin(2*phi)*np.cos(2*psi)
    Theta_proj = 0.5*np.sqrt(Fplus**2*(1 + np.cos(iota)**2)**2 + \
        4*Fcross**2*np.cos(iota)**2)
        
    SNR_out = SNR_in * Theta_proj               #applying projection function
    
    return SNR_out
