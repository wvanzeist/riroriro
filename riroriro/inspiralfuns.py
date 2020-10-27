#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parts of the procedure for simulating the inspiral portions of gravitational
waves, collected into modular functions.
"""

import numpy as np

def get_M_and_eta(**kwargs):
    """
    Gives total mass (M) and symmetric mass ratio (eta) from either m1 and m2
    OR logMc and q; M and eta are used by many functions in the GW synthesis.
    
    Parameters
    ----------
    First method
    m1: float
        Mass of one object in binary in solar masses.
    m2: float
        Mass of other object in binary in solar masses.
    Second method
    logMc: float
        log10(the chirp mass of the binary in solar masses).
    q: float
        The mass ratio of the objects in the binary.
    
    Returns
    -------
    (M,eta): tuple of floats
        The total mass of the binary, followed by the symmetric mass ratio.
    """
    
    #making sure that all of the inputs are floats
    for each_variable in kwargs.values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    #making sure that the user doesn't put arguments from both methods at once
    if ('m1' in kwargs or 'm2' in kwargs) and ('logMc' in kwargs or 'q' in \
       kwargs):
        raise TypeError('Please don\'t mix the two input methods.')
        
    #making sure that both of the arguments for a method are specified
    if (('m1' in kwargs) != ('m2' in kwargs)) or (('logMc' in kwargs) != ('q' \
        in kwargs)):
        raise TypeError('Please specify both arguments required by a method.')
        
    if 'q' in kwargs and kwargs['q'] > 1:    #making sure that q is not above 1
        raise ValueError('The mass ratio needs to be less than or equal to 1, '
                         'i.e. the smaller mass divided by the larger.')
        
    if 'm1' in kwargs:                       #first method
        M = kwargs['m1'] + kwargs['m2']
        eta = (kwargs['m1']*kwargs['m2'])/M**2    
    elif 'q' in kwargs:                      #second method
        eta = kwargs['q'] / ((1+kwargs['q'])**2)
        M = 10**kwargs['logMc'] / (eta**(3/5))
        
    return (M,eta)

def startx(M,flow):
    """
    Gives starting value/lower boundary for integration of post-Newtonian
    parameter, based on Buskirk et al. (2019) equation 22.
    
    Parameters
    ----------
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
    flow: float
        Lower cutoff frequency for the gravitational waveform, which we usually
        set to be 10 Hz.
        
    Returns
    -------
    value: float
        The starting value for the post-Newtonian integration.
    """
    
    #input type checking
    assert type(M) == float, 'M should be a float.'
    assert type(flow) == float, 'flow should be a float.'
    
    #basic constants
    pi=np.pi
    c=2.99792458e8
    #geometric unit conversion
    Msunkg=1.9891e30
    
    value = (6.673e-11 * M*Msunkg*pi*flow*c**-3)**(2/3)     #Buskirk eq. 22
    
    return value

def endx(eta,merger_type):
    """
    Gives ending value/upper boundary for integration of post-Newtonian
    parameter, based on Buskirk et al. (2019) equation 23.
    
    Parameters
    ----------
    eta: float
        Symmetric mass ratio of the binary, can be obtained from
        get_M_and_eta().
    merger_type: string
        'BH' for a BH-BH merger, 'NS' for a BH-NS or NS-NS merger
        
    Returns
    -------
    value: float
        The ending value for the post-Newtonian integration.
    """
    
    #input type checking
    assert type(eta) == float, 'eta should be a float.'
    
    if merger_type == 'BH':
        value = (1/3)*(1 + (7/18)*eta)
        #Buskirk eq. 23, with 1/6 -> 1/3 to ensure overlap with merger portion
        #for the matching script
    elif merger_type == 'NS':
        value = (1/6)*(1 + (7/18)*eta)                  #Buskirk eq. 23
        #1/6 because 3RSch is the ISCO for BH-NS and approx. the touching point
        #for NS-NS, so the end condition for both
    else:
        raise ValueError('merger_type must be either \'BH\' for BH-BH or '
                         '\'NS\' for BH-NS or NS-NS.')
    return value

def PNderiv(x,M,eta):
    """
    Encodes the differential equation for the post-Newtonian parameter that is
    integrated by x_integration(), based on Huerta et al. (2017) and Buskirk et
    al. (2019). This function should usually not be called directly, but rather
    by x_integration().
    
    Parameters
    ----------
    x: float
        The post-Newtonian parameter, the variable being integrated over.
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
    eta: float
        Symmetric mass ratio of the binary, can be obtained from
        get_M_and_eta().
    
    Returns
    -------
    Mdxdt: float
        The value of M * (dx/dt) for the input x, as given by the differential
        equation.
    """
    
    #basic constants
    pi=np.pi
    #coefficients from Huerta article
    a0 = 12.8*eta
    a2 = eta*(-16*(924*eta + 743))/420
    a3 = 51.2*pi*eta
    a4 = (1903104*eta**2 + 3934368*eta + 1091296)*(eta/45360)
    #Buskirk also includes a8 to a12 (6PN)
    a8 = 170.799 - 742.551*eta + 370.173*eta**2 - 43.4703*eta**3 - \
        0.0249486*eta**4 + (14.143 - 150.692*eta)*np.log(x)
    a9 = 1047.25 - 2280.56*eta + 923.756*eta**2 + 22.7462*eta**3 - \
        102.446*np.log(x)
    a10 = 714.739 - 1936.48*eta + 3058.95*eta**2 - 514.288*eta**3 + \
        29.5523*eta**4 - 0.185941*eta**5 + (-3.00846 + 1019.71*eta + \
        1146.13*eta**2)*np.log(x)
    a11 = 3622.99 - 11498.7*eta + 12973.5*eta**2 - 1623*eta**3 + \
        25.5499*eta**4 + (83.1435 - 1893.65*eta)*np.log(x)
    a12 = 11583.1 - 45878.3*eta + 33371.8*eta**2 - 7650.04*eta**3 + \
        648.748*eta**4 - 14.5589*eta**5 - 0.0925075*eta**6 + (-1155.61 + \
        7001.79*eta - 2135.6*eta**2 - 2411.92*eta**3)*np.log(x) + \
        33.2307*np.log(x**2)
    Mdxdt = a0*x**5 + a2*x**6 + a3*x**6.5 + a4*x**7 + ((1/M)*(64/5)*eta*x**5 \
        * (1 + a8*x**4 + a9*x**4.5 + a10*x**5 + a11*x**5.5 + a12*x**6))
    return Mdxdt

def PN_parameter_integration(start,end,M,eta):
    """
    Integrates the PNderiv() differential equation for the post-Newtonian
    parameter, x.
    
    Parameters
    ----------
    start: float
        The starting value/lower boundary of the integration, from startx().
    end: float
        The ending value/upper boundary of the integration, from endx().
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
    eta: float
        Symmetric mass ratio of the binary, can be obtained from
        get_M_and_eta().
        
    Returns
    -------
    [x,xtimes,dt]: list of lists of floats
        First list is the series of values of the post-Newtonian parameter x
        that has been integrated, second list is the time corresponding to each
        value of x (data point), third list is the timestep between each pair
        of data points.
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    dt=[]                                       #holds the timestep
    xtimes=[0]                                  #times when x is evaluated
    i=0                                         #loop counter
    #Euler integration with variable timestep
    x=[start]                                   #post-Newtonian parameter
    while x[i] <= end:
        dt.append((10**-6*x[i])/PNderiv(x[i],M,eta))
                                                #timestep with 10^-6 threshold
        x.append(x[i] + 10**-6*x[i])            #Euler step
        xtimes.append(xtimes[i] + dt[i])        #increment time-storer
        i += 1                                  #increment counter
    
    return [x,xtimes,dt]

def inspiral_time_conversion(xtimes,M):
    """
    Converting times in geometric units from x_integration() to times in real
    units.
    
    Parameters
    ----------
    xtimes: list of floats
        Times in geometric units of data points in the integration of the post-
        Newtonian parameter, from PN_parameter_integration().
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
        
    Returns
    -------
    realtimes: list of floats
        xtimes, but in seconds instead of geometric units.
    """
    
    #input type checking
    assert type(xtimes) == list, 'xtimes should be a list.'
    assert type(M) == float, 'M should be a float.'
    
    #geometric unit conversion
    Msuns=4.923e-6
    
    realtimes = np.empty((len(xtimes)))         #initialisation of list
    for i in range(len(xtimes)):
        realtimes[i] = xtimes[i]*M*Msuns
    
    #output type conversion
    realtimes = list(realtimes)
    
    return realtimes

def inspiral_phase_freq_integration(x,dt,M):
    """
    Integration of orbital phase and angular frequency for the inspiral, using
    the post-Newtonian parameter, based on Buskirk et al. (2019) equation 7.
    
    Parameters
    ----------
    x: list of floats
        Values of the post-Newtonian parameter over time, from
        PN_parameter_integration().
    dt: list of floats
        Timesteps in geometric units between each value of xtimes, from
        PN_parameter_integration().
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
        
    Returns
    -------
    [i_phase,omega,freq]: list of lists of floats
        First list is the values of orbital phase at each timestep, second list
        is the angular frequency, third list is the frequency of the GW signal.
    """
    
    #input type checking
    assert type(x) == list, 'x should be a list.'
    assert type(dt) == list, 'dt should be a list.'
    assert type(M) == float, 'M should be a float.'
    
    #basic constants
    pi=np.pi
    #geometric unit conversion
    Msuns=4.923e-6
    
    i_phase = np.empty((len(x)))                #orbital phase
    omega = np.empty((len(x)))                  #angular frequency
    freq = np.empty((len(x)))                   #frequency
    
    #initial values
    i_phase[0]=0
    omega[0]=x[0]**1.5
    freq[0]=omega[0]/(M*Msuns)
    
    #for the phase integration, we are tied to the times at which we evaluated
    #x
    for i in range(len(x)-1):                   #-1 because always i+1 terms
        omega[i+1] = x[i+1]**1.5
        i_phase[i+1] = i_phase[i] + omega[i+1]*dt[i]
        freq[i+1] = omega[i+1] / (M*Msuns*pi)
    
    #output type conversion
    i_phase = list(i_phase)
    omega = list(omega)
    freq = list(freq)
    
    return [i_phase,omega,freq]

def radius_calculation(x,M,eta):
    """
    Calculation of orbital radius (and time-derivative of radius) for the
    binary for each timestep during the inspiral, based on Buskirk et al.
    (2019).
    
    Parameters
    ----------
    x: list of floats
        Values of the post-Newtonian parameter over time, from
        PN_parameter_integration().
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
    eta: float
        Symmetric mass ratio of the binary, can be obtained from
        get_M_and_eta().
        
    Returns
    -------
    [r,rdot]: list of lists of floats
        First list is the values of the orbital radius (in geometric units) at
        each timestep, second list is the time-derivative of the radius (used
        by strain calculations).
    """
    
    #input type checking
    assert type(x) == list, 'x should be a list.'
    assert type(M) == float, 'M should be a float.'
    assert type(eta) == float, 'eta should be a float.'
    
    #post-Newtonian correction coefficients from Buskirk
    r0pn = 1
    r1pn = -1 + 0.333333*eta
    r2pn = 4.75*eta + 0.111111*eta**2
    r3pn = -7.51822*eta - 3.08333*eta**2 + 0.0246914*eta**3
    
    r = np.empty((len(x)))                      #orbital radius (geometric u.)
    rdot = np.empty((len(x)))                   #derivative of radius
    
    for i in range(len(x)):                     #Buskirk radius equations
        r[i] = r0pn*(1/x[i]) + r1pn + r2pn*x[i] + r3pn*x[i]**2
        rdot[i] = PNderiv(x[i],M,eta) * (-2*r0pn*x[i]**-2 + r2pn + 2*r3pn*x[i])
        
    #output type conversion
    r = list(r)
    rdot = list(rdot)
    
    return [r,rdot]

def a1_a2_calculation(r,rdot,omega,D,M,eta):
    """
    Calculation of A1 and A2, two coefficients used in the calculation of
    strain polarisations, based on Buskirk et al. (2019) equation 9.
    
    Parameters
    ----------
    r: list of floats
        Values of the orbital radius over time, from radius_calculation().
    rdot: list of floats
        Values of the time-derivative of the radius, from radius_calculation().
    omega: list of floats
        Values of the angular frequency over time, from
        inspiral_phase_freq_integration().
    D: float
        Distance from the detector to the binary, in Mpc. IMPORTANT: if you
        want to feed the strain values into the SNR calculator, use the default
        distance of 100 Mpc here and instead set the distance when using the
        SNR functions.
    M: float
        Total mass of the binary, can be obtained from get_M_and_eta().
    eta: float
        Symmetric mass ratio of the binary, can be obtained from
        get_M_and_eta().
        
    Returns
    -------
    [A1,A2]: list of lists of floats
        The first list is the values  of the A1 parameter used in strain
        calculation over time, the second list is the A2 parameter.
    """
    
    #input type checking
    assert type(r) == list, 'r should be a list.'
    assert type(rdot) == list, 'rdot should be a list.'
    assert type(omega) == list, 'omega should be a list.'
    assert type(D) == float, 'D should be a float.'
    assert type(M) == float, 'M should be a float.'
    assert type(eta) == float, 'eta should be a float.'
    
    Dkm = D * 3.086e19                          #conversion from Mpc to km
    
    A1 = np.empty((len(r)))
    A2 = np.empty((len(r)))
    
    for i in range(len(r)):                     #based on Buskirk eq. 9
        A1[i] = (-2*M*eta*(1/Dkm))*(rdot[i]**2 + (r[i]*omega[i])**2 + 1/r[i])
        A2[i] = (-2*M*eta*(1/Dkm))*(2*r[i]*rdot[i]*omega[i])
    
    #output type conversion
    A1 = list(A1)
    A2 = list(A2)
    
    return [A1,A2]

def inspiral_strain_polarisations(A1,A2,i_phase):
    """
    Calculating the values of the two polarisations of strain for the inspiral,
    using the coefficients from a1_a2_calculation().
    
    Parameters
    ----------
    A1: list of floats
        Values of the first strain coefficient over time, from
        a1_a2_calculation().
    A2: list of floats
        Values of the second strain coefficient over time, from
        a1_a2_calculation().
    i_phase: list of floats
        Values of the orbital phase at each timestep, from
        inspiral_phase_freq_integration().
        
    Returns
    -------
    [Aorth,Adiag]: list of lists of floats
        The first list is the values of the orthogonal/plus polarisation of
        strain over time, the second list is the diagonal/cross polarisation.
    """
    
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'
        
    Aorth = np.empty((len(i_phase)))            #orthogonal/plus polarisation
    Adiag = np.empty((len(i_phase)))            #diagonal/cross polarisation
        
    for i in range(len(i_phase)):
        Aorth[i] = A1[i]*np.cos(2*i_phase[i]) + A2[i]*np.sin(2*i_phase[i])
        Adiag[i] = A1[i]*np.sin(2*i_phase[i]) - A2[i]*np.cos(2*i_phase[i])
        
    #output type conversion
    Aorth = list(Aorth)
    Adiag = list(Adiag)
    
    return [Aorth,Adiag]

def inspiral_strain_amplitude(Aorth,Adiag):
    """
    Calculating the amplitude of the strain from the polarisations.
    
    Parameters
    ----------
    Aorth: list of floats
        The values of the orthogonal/plus polarisation of strain over time,
        from inspiral_strain_polarisations().
    Adiag: list of floats
        The values of the diagonal/cross polarisation of strain over time, from
        inspiral_strain_polarisations().
        
    Returns
    -------
    i_amp: list of floats
        The values of the amplitude of the GW strain over time (unitless).
    """
    
    #input type checking
    assert type(Aorth) == list, 'Aorth should be a list.'
    assert type(Adiag) == list, 'Adiag should be a list.'
    
    i_amp = np.empty((len(Aorth)))
    for i in range(len(Aorth)):
        i_amp[i] = np.sqrt(Aorth[i]**2 + Adiag[i]**2)
        
    #output type conversion
    i_amp = list(i_amp)
    
    return i_amp

def list_size_reducer(reduction_factor,your_list):
    """
    Optional function to reduce the size of the lists output by the inspiral
    functions (not the merger lists, as those are much shorter), in order to
    reduce filesize to conserve storage space.
    NOTES:
    The typical reduction factor we have used in our research using this code
    is 100.
    The inspiral lists used by the matching/merger portions are realtimes,
    omega, i_phase and i_amp so if you reduce one of these you should reduce
    all of them.
    
    Parameters
    ----------
    reduction_factor: int
        The factor you want to reduce the list length by.
    your_list: list
        The list you want to reduce.
        
    Returns
    -------
    reduced_list: list
        your_list, in reduced form.
    """
    
    #input type checking
    assert type(reduction_factor) == int, 'reduction_factor should be an int.'
    assert type(your_list) == list, ('The thing to be reduced needs to be a '
                                     'list.')
    
    #create new list with every nth point of your_list
    reduced_list = [your_list[0]]
    for i in range(reduction_factor,len(your_list),reduction_factor):
        reduced_list.append(your_list[i])
        
    return reduced_list