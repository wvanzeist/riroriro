**********
gwexporter
**********

This is the documentation for the gwexporter module, which consists of functions to output the gravitational waveforms from inspiralfuns, mergerfirstfuns, matchingfuns and mergersecondfuns into a file.

waveform_exporter
=================

``waveform_exporter(time,freq,amp,path)``

Function to export a simulated gravitational waveform into a file.

Parameters
----------
time: list of floats
    The time at each data point. For a BH-BH merger, use i_m_time from
    time_frequency_stitching in matchingfuns. For a BH-NS or NS-NS merger,
    use i_time (realtimes) from inspiral_time_conversion in inspiralfuns.
freq: list of floats
    The frequency of the GW signal at each data point. For a BH-BH merger,
    use i_m_freq from frequency_SI_units in matchingfuns. For a BH-NS or
    NS-NS merger, use i_freq (freq) from inspiral_phase_freq_integration
    in inspiralfuns.
amp: list of floats
    The amplitude of the GW strain at each data point. For a BH-BH merger,
    use i_m_amp from amplitude_stitching in mergersecondfuns. For a BH-NS
    or NS-NS merger, use i_amp from inspiral_strain_amplitude in
    inspiralfuns.
path: str
    The file path to the location/document where you want to save the
    simulated gravitational waveform data.
    
Returns
-------
An output file containing an array wherein the first column is the time,
the second is the frequency and the third is the amplitude.

waveform_arrayer
================

``waveform_arrayer(time,freq,amp)``

Function to collate important data of the simulated gravitational waveform
(for SNR calculation etc.) into a single array for ease of storage. Similar
to waveform_exporter(), but instead of outputting the area into a file,
this function outputs the data inline in Python, in a numpy.ndarray.

Parameters
----------
time: list of floats
    The time at each data point. For a BH-BH merger, use i_m_time from
    time_frequency_stitching in matchingfuns. For a BH-NS or NS-NS merger,
    use i_time (realtimes) from inspiral_time_conversion in inspiralfuns.
freq: list of floats
    The frequency of the GW signal at each data point. For a BH-BH merger,
    use i_m_freq from frequency_SI_units in matchingfuns. For a BH-NS or
    NS-NS merger, use i_freq (freq) from inspiral_phase_freq_integration
    in inspiralfuns.
amp: list of floats
    The amplitude of the GW strain at each data point. For a BH-BH merger,
    use i_m_amp from amplitude_stitching in mergersecondfuns. For a BH-NS
    or NS-NS merger, use i_amp from inspiral_strain_amplitude in
    inspiralfuns.
    
Returns
-------
exportarray: numpy.ndarray
    An array containing the important data of the simulated gravitational
    waveform: the first column is the time, the second is the frequency and
    the third is the amplitude.
