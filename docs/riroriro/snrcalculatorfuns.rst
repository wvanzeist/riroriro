*****************
snrcalculatorfuns
*****************

This is the documentation for the snrcalculatorfuns module, which consists of parts of the procedure for calculating the SNR of a gravitational waveform (from gwexporter or otherwise), collected into modular functions.

polynomial_redshift
===================

``polynomial_redshift(d)``

Polynomial approximation of calculating redshift corresponding to a given
distance.

Parameters
----------
d: float
    A luminosity distance, in Mpc.
    
Returns
-------
z: float
    The redshift corresponding to the input distance.
    
redshift_distance_adjustment
============================

``redshift_distance_adjustment(inputarray,d,z)``

Adjusts the frequencies and amplitudes in the input gravitational waveform
to account for the effects of distance/redshift.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform,
    in the format used by waveform_exporter() in gwexporter.
d: float
    The luminosity distance to the merging binary, in Mpc.
z: float
    The redshift corresponding to the input distance.
    
Returns
-------
adjustedarray: numpy.ndarray
    inputarray, but with the frequency and amplitudes adjusted.
    
frequency_limits
================

``frequency_limits(inputarray)``

Calculates the upper and lower limits of the frequency of the gravitational
waveform in inputarray, which are used by amplitude_interpolation().

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform;
    should have been adjusted by redshift_distance_adjustment().
    
Returns
-------
(freqmax,freqmin): tuple of floats
    The upper and lower limits of the waveform signal frequency,
    respectively.
    
findchirp_fourier
=================

``findchirp_fourier(inputarray,findchirp_array,d,z)``

Approximation of a Fourier transform on the gravitational waveform data,
using the frequency spectrum output by the simpler model FINDCHIRP (Allen
et al., 2012) for calibration.
NOTE: May in the future be replaced by something fft-based.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform;
    should have been adjusted by redshift_distance_adjustment().
findchirp_array: numpy.ndarray
    The array output by FINDCHIRP. The second column is frequency, the
    fourth is (Fourier-transformed) strain amplitude, the other columns
    are irrelevant. A grid of sample findchirp_arrays can be found at
    https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
d: float
    The luminosity distance to the merging binary, in Mpc.
z: float
    The redshift corresponding to the input distance.
    
Returns
-------
fourieramp: list
    Fourier-transformed/calibrated amplitudes at each frequency value in
    inputarray.
    
amplitude_interpolation
=======================

``amplitude_interpolation(inputarray,fourieramp,noisearray,freqmax,freqmin)``

The simulated gravitational waveform data and the detector noise spectrum
are assumed to have amplitude data at different sets of frequencies, so
this function uses scipy's interp1d to calculate the waveform amplitude
values at the frequencies used by the detector data.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform;
    should have been adjusted by redshift_distance_adjustment().
fourieramp: list
    Fourier-transformed/calibrated amplitudes at each frequency value in
    inputarray, from findchirp_fourier().
noisearray: numpy.ndarray
    Data on the noise spectrum of the detector; it is assumed that
    frequency values are in the first column and ASD noise levels in the
    second.
freqmax: float
    The upper limit of the waveform signal frequency, from
    frequency_limits().
freqmin: float
    The lower limit of the waveform signal frequency, from
    frequency_limits().

Returns
-------
noise_freq_amp: list
    Waveform amplitudes as in fourieramp, but over the set of frequencies
    in noisearray rather than those in inputarray.
    
individual_detector_SNR
=======================

``individual_detector_SNR(noisearray,noise_freq_amp)``

Calculates the single-detector optimal-alignment SNR by comparing the
waveform frequency spectrum and detector noise spectrum using the method of
Barrett et al. (2018).

Parameters
----------
noisearray: numpy.ndarray
    Data on the noise spectrum of the detector; it is assumed that
    frequency values are in the first column and ASD noise levels in the
    second.
noise_freq_amp: list
    Amplitudes of the simulated gravitational waveform, over the set of
    frequencies of noisearray, from amplitude_interpolation().
    
Returns
-------
ind_SNR: float
    The SNR of the simulated gravitational waveform, for the detector in
    noisearray and assuming optimal alignment.
