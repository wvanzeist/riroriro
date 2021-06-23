#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parts of the procedure for calculating the SNR of a gravitational waveform
(from gwexporter or otherwise), collected into modular functions.
"""

import numpy as np
from scipy.interpolate import interp1d

def polynomial_redshift(d):
    """
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
    """
    
    #input type checking
    assert type(d) == float, 'd should be a float.'
    
    #sanity check: distance should not be negative
    assert d >= 0, 'The distance should be a positive number.'
    
    #polynomial approximation of redshift conversion
    z = 1.0832e-12*d**3 - 1.7022e-8*d**2 + 0.00021614*d
    
    return z

def redshift_distance_adjustment(inputarray,d,z):
    """
    Adjusts the frequencies and amplitudes in the input gravitational waveform
    to account for the effects of distance/redshift.
    
    Parameters
    ----------
    inputarray: numpy.ndarray
        The time, frequency and amplitude data of the gravitational waveform,
        in the format used by either waveform_exporter_3col() or
        waveform_exporter_4col() in gwexporter.
    d: float
        The luminosity distance to the merging binary, in Mpc.
    z: float
        The redshift corresponding to the input distance.
        
    Returns
    -------
    adjustedarray: numpy.ndarray
        inputarray, but with the frequency and amplitudes adjusted.
    """
    
    #input type checking
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    assert type(d) == float, 'd should be a float.'
    assert type(z) == float, 'z should be a float.'
    
    adjustedarray = np.zeros(inputarray.shape)
    
    for i in range(inputarray.shape[0]):
        adjustedarray[i,0] = inputarray[i,0]
        adjustedarray[i,1] = inputarray[i,1] / (1+z)    #frequency redshifting
        adjustedarray[i,2] = inputarray[i,2] / (d/100)  #distance adjustment
        if inputarray.shape[1] == 4:                #with both polarisations
            adjustedarray[i,3] = inputarray[i,3] / (d/100)
        
    return adjustedarray

def frequency_limits(inputarray):
    """
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
    """
    
    #input type checking
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    
    freqmax=np.amax(inputarray[:,1])
    freqmin=np.amin(inputarray[:,1])
    
    #output type conversion
    freqmax=float(freqmax)
    freqmin=float(freqmin)
    
    return (freqmax,freqmin)

def findchirp_fourier(inputarray,findchirp_array,d,z):
    """
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
    """
    
    #input type checking
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    assert type(findchirp_array) == np.ndarray, ('findchirp_array should be an'
                                                 ' array.')
    assert type(d) == float, 'd should be a float.'
    assert type(z) == float, 'z should be a float.'
    
    #findchirparray scaling and redshifting
    adj_findchirp_array = np.empty((findchirp_array.shape))
    for i in range(findchirp_array.shape[0]):
        adj_findchirp_array[i,0] = findchirp_array[i,0]
        adj_findchirp_array[i,1] = findchirp_array[i,1]/(1+z)
        #redshift adjustment for trace
        adj_findchirp_array[i,2] = findchirp_array[i,2]
        adj_findchirp_array[i,3] = findchirp_array[i,3]/d
        #calculated assuming 1 Mpc, hence scaling
    
    #FINDCHIRP-inspired empirical approximation of Fourier transform
    fourieramp = np.empty((inputarray.shape[0]))
    for i in range(len(fourieramp)):
        fourieramp[i] = (inputarray[i,2]) / (inputarray[i,1]**(11/6))
                        #*fourier_freq**(5/6)          #* 1.15 (trace scaling)
                        #11/6 in this version with adjusted traces
                        
    #scaling to match Fourier-transformed waveform and FINDCHIRP trace at 10 Hz
    fourier_10Hz = np.searchsorted(inputarray[:,1],10)
    trace_10Hz = np.searchsorted(adj_findchirp_array[:,1],10)
    f_t_ratio = fourieramp[fourier_10Hz] / adj_findchirp_array[trace_10Hz,3]
    for i in range(len(fourieramp)):
        fourieramp[i] /= f_t_ratio #scale so amplitudes are equal at 10 Hz
        
    #output type conversion
    fourieramp = list(fourieramp)
    
    return fourieramp

def proper_fourier(inputarray,freqmax):
    """
    Fourier transform of the gravitational wave amplitudes to the frequencies,
    using a fft function.
    NOTE: This function is intended to replace findchirp_fourier(). It requires
    the version of inputarray that has the strain polarisations.
    
    Parameters
    ----------
    inputarray: numpy.ndarray
        The time, frequency and strain polarisation data of the gravitational
        waveform; should have been adjusted by redshift_distance_adjustment().
    freqmax: float
        The upper limit of the waveform signal frequency, from
        frequency_limits().
        
    Returns
    -------
    fourieramp: list
        Fourier-transformed/calibrated amplitudes at each frequency value in
        inputarray.
    """
    
    #import scipy.signal as signal
    #from astropy.timeseries import LombScargle
    from scipy.interpolate import interp1d
    
    #input type checking
    assert type(inputarray) == np.ndarray, 'inputarray should be an array.'
    assert type(freqmax) == float, 'freqmax should be a float.'
    
    """ #Lomb-Scargle methods
    #complex amplitude (combining polarisations) to be fft input
    fourierinput = np.empty((len(inputarray)),dtype=complex)
    for i in range(len(fourierinput)):
        fourierinput[i] = inputarray[i,2] - 1j*inputarray[i,3]
    
    #we cannot use a normal fft here because the measurements are not evenly
    #spaced in time; we use the alternative Lomb-Scargle method instead
    #fourieramp = signal.lombscargle(inputarray[:,0],fourierinput, \
    #                                inputarray[:,1])
    fourieramp = LombScargle(inputarray[:,0], \
                             fourierinput).power(inputarray[:,1])
    """
    
    #FFT functions require amplitude measurements evenly spaced in time, so we
    #have to interpolate the amplitude polarisation values
    Aorth_inter = interp1d(inputarray[:,0], inputarray[:,2], kind='cubic')
    Adiag_inter = interp1d(inputarray[:,0], inputarray[:,3], kind='cubic')
    #For the array we feed to the FFT, we want to balance accuracy with
    #CPU time and memory; the minimum sampling rate to prevent aliasing is the
    #Nyquist rate, twice the highest frequency in the signal
    nyquist_period = 1/(4*freqmax)      #period corresponding to Nyquist freq.
    #(except freq. is doubled to reduce FFT glitching at high frequencies)
    #This period won't be an exact divisor of the signal's duration, so we look
    #for the next shorter period that *is*.
    signal_duration = inputarray[len(inputarray)-1,0]
    nyquist_time_ratio = signal_duration / nyquist_period
    no_of_samples = int(np.ceil(nyquist_time_ratio))
    #higher ratio, lower period (becomes float if not set as int)
    sampling_period = signal_duration / no_of_samples
    
    #now we use amp_time_inter to generate measurements at this period
    sampling_times = np.linspace(0,signal_duration,no_of_samples)
    even_fourierinput = np.empty((no_of_samples),dtype=complex)
    for i in range(no_of_samples):
        even_fourierinput[i] = Aorth_inter(sampling_times[i]) - \
            1j*Adiag_inter(sampling_times[i])
    
    raw_fourieramp = np.fft.fft(even_fourierinput)      #Fourier transform
    raw_fourierfreqs = np.fft.fftfreq(no_of_samples,d=sampling_period)
    
    #removing irrelevant negative frequencies from the output
    cplx_fourieramp = raw_fourieramp[:int(len(raw_fourieramp)/2)]
    fourierfreqs = raw_fourierfreqs[:int(len(raw_fourierfreqs)/2)]
    
    #converting FFT strain values from complex to absolute/real amplitude
    fourieramp = np.empty((len(cplx_fourieramp)))
    for i in range(len(fourieramp)):
        fourieramp[i] = abs(cplx_fourieramp[i])
    
    #output type conversion
    fourieramp = list(fourieramp)
    fourierfreqs = list(fourierfreqs)
    
    #return fourieramp    
    #return even_fourierinput, fourieramp         #testing
    return [fourieramp,fourierfreqs]

def amplitude_interpolation(freqinput,fourieramp,noisearray,freqmax,freqmin):
    """
    The simulated gravitational waveform data and the detector noise spectrum
    are assumed to have amplitude data at different sets of frequencies, so
    this function uses scipy's interp1d to calculate the waveform amplitude
    values at the frequencies used by the detector data.
    
    Parameters
    ----------
    freqinput: either inputarray if using findchirp_fourier, or fourier_freqs
        if using proper_fourier.
        inputarray: numpy.ndarray
            The time, frequency and amplitude data of the gravitational
            waveform; should have been adjusted by
            redshift_distance_adjustment().
        fourierfreqs: list
            The frequency values corresponding to each Fourier-transformed
            amplitude measurement, from findchirp_fourier.
    fourieramp: list
        Fourier-transformed/calibrated amplitudes at each frequency value in
        inputarray, from findchirp_fourier() or proper_fourier.
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
    """
    
    #input type checking
    assert (type(freqinput) == np.ndarray) or (type(freqinput) == list), \
        ('freqinput should be an array (if using findchirp_fourier) or a list '
        '(if using proper_fourier).')
    assert type(fourieramp) == list, 'fourieramp should be a list.'
    assert type(noisearray) == np.ndarray, 'noisearray should be an array.'
    assert type(freqmax) == float, 'freqmax should be a float.'
    assert type(freqmin) == float, 'freqmin should be a float.'
    
    if type(freqinput) == np.ndarray:   #inputarray, from findchirp_fourier
        smoothinput = interp1d(freqinput[:,1], fourieramp[:], kind='cubic')
    else:                               #fourierfreqs, from proper_fourier
        smoothinput = interp1d(freqinput[:], fourieramp[:], kind='cubic')
    #interpolating input frequency-amp curve so it can be calculated for the
    #other grid of frequencies used by noisearray
    
    noise_freq_amp = np.empty((len(noisearray)))
    
    for i in range(len(noise_freq_amp)):
        if noisearray[i,0] > freqmax or noisearray[i,0] < freqmin:
            noise_freq_amp[i] = 0
            #the smoothinput function is only defined between freqmin and
            #freqmax; this just sets the amplitude to zero outside the actual
            #simulated range of the waveform
        else:
            noise_freq_amp[i] = smoothinput(noisearray[i,0])
            
    #output type conversion
    noise_freq_amp = list(noise_freq_amp)
    
    return noise_freq_amp

def individual_detector_SNR(noisearray,noise_freq_amp):
    """
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
    """
    
    #input type checking
    assert type(noisearray) == np.ndarray, 'noisearray should be an array.'
    assert type(noise_freq_amp) == list, 'noise_freq_amp should be a list.'
    
    noiseratio = np.empty((len(noise_freq_amp)))
    
    #h/S part of SNR equation
    for i in range(len(noiseratio)):
        noiseratio[i] = (noise_freq_amp[i] / noisearray[i,1])**2
        #squaring is to convert ASD to PSD and amp to h*(f)h(f)
        #f to balance ASD units (from noise and from df) is no longer
        #needed when we have FINDCHIRP-based transform
        #(the freqmax-freqmin stuff that was here in the original SNR
        #calculator script is moved to the noise_freq_amp calculation in this
        #version)
        
    SNR_df = np.empty((len(noise_freq_amp)))
    
    #df in SNR equation
    for i in range(len(SNR_df)):                #derivative as differences
        #unusual averages at ends (not that these will be relevant here)
        if i == 0:
            SNR_df[i] = noisearray[1,0] - noisearray[0,0]
        elif i == noisearray.shape[0] - 1:              #final value in array
            SNR_df[i] = noisearray[i,0] - noisearray[i-1,0]
            #regular averages, based on those used in fhatdot calculation
        else:
            SNR_df[i] = 0.5*(noisearray[i+1,0] - noisearray[i-1,0])
    
    #integrating over frequency as Riemann sum
    sumpart = np.empty(noisearray.shape[0])             #slices of SNR integral
    for i in range(len(sumpart)):
        sumpart[i] = noiseratio[i]*SNR_df[i]
    
    #final part of Barrett et al. calculation
    ind_SNR = np.sqrt(4*sum(sumpart))
    
    return ind_SNR
