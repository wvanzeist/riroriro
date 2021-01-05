---
title: 'Riroriro: Simulating gravitational waves and calculating SNRs in Python'
tags:
  - Python
  - astronomy
  - gravitational waves
  - black holes
  - neutron stars
authors:
  - name: Wouter G. J. van Zeist
    affiliation: 1
  - name: Héloïse F. Stevance
    affiliation: 1
  - name: J. J. Eldridge
    affiliation: 1
affiliations:
  - name: Department of Physics, University of Auckland, New Zealand
    index: 1
date: 2021
bibliography: paper.bib
---

# Summary

`Riroriro` is a Python package to simulate the gravitational waveforms of binary mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors.

The gravitational waveform simulation of `Riroriro` is based upon the methods of @buskirk2019, a paper which describes a computational implementation of an earlier theoretical gravitational waveform model by @huerta2017, using post-Newtonian expansions and an approximation called the implicit rotating source to simplify the Einstein field equations and simulate gravitational waves. `Riroriro`'s calculation of signal-to-noise ratios (SNR) of gravitational wave events is based on the methods of @barrett2018, with the simpler gravitational wave model `Findchirp` [@findchirp] being used for comparison and calibration in these calculations.

# Statement of Need

Gravitational waves have long been an area of research in astronomy, and particularly since the first observation of gravitational waves was announced in 2016 [@gw150914discovery] this area has seen a lot of research activity. Observations of gravitational waves from binary mergers can provide unique information about their progenitors and stellar populations, especially when combined with electromagnetic observations in the field called multi-messenger astronomy. A major factor in the successful detection and analysis of gravitational wave signals is the creation of simulations of such signals which observed data can be compared to. Because of this, multiple gravitational wave models have been created over the years.

The gravitational wave observatories LIGO/Virgo use their own models as templates, and various other research groups have also created models. For example, prior to @huerta2017 and @buskirk2019 there were other simulations of gravitational waves using similar methods of post-Newtonian expansions and/or limited numerical relativity simulations [@apostolatos1995; @buonanno2003; @ajith2007], though these had less mathematical accuracy than the later models. There are also models using alternative methods to the common post-Newtonian expansions, such as representing the two-body system of merging binaries as an effective one-body system [@buonanno1999], solving full numerical relativity equations but only at limits of the binary interaction and extrapolating for the rest of the waveform [@baker2002a; @baker2002b] or using numerically greedy algorithms to predict the shape of a wide range of waveforms based on a small number of given waveforms [@field2014].

However, for the most part these articles do not go into much detail about the computational implementation of their simulations and have not published their full database or its source code. Some have published it but in a language that is not free to use, like the original model of @buskirk2019 which is written in Mathematica. There are some open-source software packages in free-to-use programming languages that relate to gravitational waves, like `Bilby` [@bilby] and `PyCBC` [@pycbc] in Python, but these are more focused on performing analysis and parameter estimation on given signals rather than the forward modelling of signals performed by `Riroriro` and the models mentioned in the previous paragraph.

`Riroriro` combines areas covered by previous models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility in Python, a programming language that is ubiquitous in astronomy. Aside from being a research tool, `Riroriro` is also designed to be easy to use and modify, and it can also be used as an educational tool for students learning about gravitational waves.

# Features

Features of `Riroriro` include:

- Simulating the gravitational waveform signal from a binary merger of two black holes, two neutron stars or a black hole and a neutron star and outputting the data of this signal in terms of frequency and strain amplitude.
- Using a gravitational wave output and given a detector noise spectrum (such spectra are made publicly available by LIGO), calculating the signal-to-noise ratio (SNR) of the signal at a given distance assuming optimal alignment.
- Calculating the horizon distance (maximum distance at which an event could be observed) for a gravitational wave model and a given detector.
- Given the optimal-alignment SNR of an event, evaluating its detectability, the probability that the event would be detected with a SNR above the commonly used threshold of 8, if the alignment would be arbitrary.

In addition, we have created Jupyter Notebook tutorials to help users get started with `Riroriro` ([see tutorials](https://github.com/wvanzeist/riroriro_tutorials)).

# Research

`Riroriro` has been used for research in conjunction with `BPASS`, a suite of computer programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions and predicts their electromagnetic spectral emission [@bpass1; @bpass2]. There is also a Python interface for `BPASS` called `Hoki` [@hoki]. This research took rates of formation of merging systems from `BPASS` and then evaluated the detectability of the gravitational wave signals from those systems using `Riroriro`. This was done to obtain predictions of the rates at which gravitational waves of different types would be expected to be observed, which can then be directly compared to those events found by the LIGO/Virgo gravitational wave observatories [@massdistribution].

# Acknowledgments

HFS and JJE acknowledge support from the University of Auckland and also the Royal Society of New Zealand Te Apārangi under the Marsden Fund​.

# References