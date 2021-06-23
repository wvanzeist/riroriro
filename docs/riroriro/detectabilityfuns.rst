*****************
detectabilityfuns
*****************

This is the documentation for the detectabilityfuns module, which consists of parts of the procedure for calculating the detectability fraction of a merger given its optimal-alignment SNR, collected into modular functions.

cdf_generator
=============

``cdf_generator(N=10**6)``

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

detectability_calculator
========================

``detectability_calculator(Theta_CDF,min_CDF,max_CDF,SNR_in)``

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

specific_orientation_SNR
========================

``specific_orientation_SNR(theta,phi,iota,psi,SNR_in,angle_unit='rad')``

Given the optimal-alignment SNR of a merger, this function returns the SNR
that would result if the detector and binary had a specific orientation/
alignment, specified by four angles as in Finn (1996), Belczynski et al.
(2013), Belczynski et al. (2014).

Parameters
----------
theta: float
    The relative latitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary). Ranges
    from 0 to π rad (180 deg).
phi: float
    The relative longitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary). Ranges
    from 0 to 2π rad (360 deg).
iota: float
    The inclination angle of the binary. Ranges from 0 to π rad (180 deg).
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
