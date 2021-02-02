***************
horizondistfuns
***************

This is the documentation for the horizondistfuns module, which consists of parts of the procedure for calculating the horizon distance of a gravitational waveform (from gwexporter or otherwise), collected into modular functions.

compact_SNR_calculation
=======================

``compact_SNR_calculation(inputarray,findchirp_array,noisearray_list,method,d)``

Runs through all of the functions of snrcalculatorfuns to obtain a SNR from
an individual detector. This function is mainly included not to be called
directly, but rather by horizon_distance_calculation().

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform,
    in the format used by waveform_exporter() in gwexporter.
findchirp_array: numpy.ndarray
    The array output by FINDCHIRP. The second column is frequency, the
    fourth is (Fourier-transformed) strain amplitude, the other columns
    are irrelevant. A grid of sample findchirp_arrays can be found at
    https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
noisearray_list: list of numpy.ndarrays
    Each item in this list should be an array describing the noise spectrum
    of a detector; in each noise spectrum, it is assumed that frequency
    values are in the first column and ASD noise levels in the second.
method: str
    If 'quad', returns the quadrature SNR across the detectors in
    noisearray_list. If 'mean', returns the mean of the SNRs with each
    individual detector (simulating one random detector in operation). If
    only one detector is included in noisearray_list, these methods are
    equivalent.
d: float
    The luminosity distance to the merging binary, in Mpc.
    
Returns
-------
final_SNR: float
    The SNR of the simulated gravitational waveform, for the detectors in
    noisearray and assuming optimal alignment.
    
horizon_distance_calculation
============================

``horizon_distance_calculation(inputarray,findchirp_array,noisearray_list,method)``

Calculates the horizon distance (maximum distance at which something can
be observed) given optimal alignment for a given merger.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform,
    in the format used by waveform_exporter() in gwexporter.
findchirp_array: numpy.ndarray
    The array output by FINDCHIRP. The second column is frequency, the
    fourth is (Fourier-transformed) strain amplitude, the other columns
    are irrelevant. A grid of sample findchirp_arrays can be found at
    https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
noisearray_list: list of numpy.ndarrays
    Each item in this list should be an array describing the noise spectrum
    of a detector; in each noise spectrum, it is assumed that frequency
    values are in the first column and ASD noise levels in the second.
method: str
    If 'quad', uses the quadrature SNR across the detectors in
    noisearray_list. If 'mean', uses the mean of the SNRs with each
    individual detector (simulating one random detector in operation). If
    only one detector is included in noisearray_list, these methods are
    equivalent.
    
Returns
-------
horizon_dist: float
    The horizon distance of the given merger, for the given detector(s).
