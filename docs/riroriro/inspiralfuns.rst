************
inspiralfuns
************

This is the documentation for the inspiralfuns module, which consists of parts of the procedure for simulating the inspiral portions of gravitational waves, collected into modular functions.

get_M_and_eta
=============

``get_M_and_eta(**kwargs)``

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

startx
======

``startx(M,flow)``

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

endx
====

``endx(eta,merger_type)``

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

PNderiv
=======

``PNderiv(x,M,eta)``

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

PN_parameter_integration
========================

``PN_parameter_integration(start,end,M,eta)``

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
    
inspiral_time_conversion
========================

``inspiral_time_conversion(xtimes,M)``

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
    
inspiral_phase_freq_integration
===============================

``inspiral_phase_freq_integration(x,dt,M)``

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
    
radius_calculation
==================

``radius_calculation(x,M,eta)``

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

a1_a2_calculation
=================

``a1_a2_calculation(r,rdot,omega,D,M,eta)``

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

inspiral_strain_polarisations
=============================

``inspiral_strain_polarisations(A1,A2,i_phase)``

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
    
inspiral_strain_amplitude
=========================

``inspiral_strain_amplitude(Aorth,Adiag)``

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

list_size_reducer
=================

``list_size_reducer(reduction_factor,your_list)``

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
