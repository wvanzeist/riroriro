*****************
detectabilityfuns
*****************

This is the documentation for the detectabilityfuns module, which consists of parts of the procedure for calculating the detectability fraction of a merger given its optimal-alignment SNR, collected into modular functions.

cdf_generator
=============

cdf_generator(N=10**6)

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

detectability_calculator(Theta_CDF,min_CDF,max_CDF,SNR_in)

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
