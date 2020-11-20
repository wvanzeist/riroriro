---
title: 'Riroriro: A code for simulating gravitational waves and calculating SNRs in Python'
tags:
  - Python
  - astronomy
  - gravitational waves
  - black holes
  - neutron stars
authors:
  - name: Wouter G. J. van Zeist
    affiliation: 1
  - name: J. J. Eldridge
    affiliation: 1
affiliations:
  - name: Department of Physics, University of Auckland, New Zealand
    index: 1
date: 2020
bibliography: paper.bib
---

# Summary

`Riroriro` is a Python package containing several modules to simulate the gravitational waveforms of binary mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors.

The gravitational waveform simulation of `Riroriro` is based upon the methods of @buskirk2019, a paper which describes a computational implementation of an earlier theoretical gravitational waveform model by @huerta2017, using post-Newtonian expansions and an approximation called the implicit rotating source to simplify the Einstein field equations and simulate gravitational waves. `Riroriro`'s calculation of signal-to-noise ratios (SNR) of gravitational wave events is based on the methods of @barrett2018, with the simpler gravitational wave model `Findchirp` [@findchirp] being used for comparison and calibration in these calculations.

# Motivations

Since the detection of gravitational waves was first announced in 2016, these have become a major area of research activity in astrophysics. A major factor in the successful detection and analysis of gravitational wave signals is the creation of simulations of such signals which observed data can be compared to. Because of this, many different gravitational wave models have been created over the years.

The gravitational wave observatories LIGO/Virgo use their own models as templates, and various other research groups have also created models. For example, prior to @huerta2017 and @buskirk2019 there were several other simulations of gravitational waves using similar methods of post-Newtonian expansions and/or limited numerical relativity simulations [@apostolatos1995; @buonanno2003; @ajith2007], though these had less mathematical accuracy than the later models. There are also models using alternative methods to the common post-Newtonian expansions, such as representing the two-body system of merging binaries as an effective one-body system [@buonanno1999], solving full numerical relativity equations but only at limits of the binary interaction and extrapolating for the rest of the waveform [@baker2002a; @baker2002b] or using numerically greedy algorithms to predict the shape of a wide range of waveforms based on a small number of given waveforms [@field2014]. However, for the most part these articles do not go into much detail about the computational implementation of their simulations and have not published their full database or its source code, or have published it but in a language that is not free to use, like the original model of @buskirk2019 which is written in `Mathematica`.

We intend for `Riroriro` to have value by combining together areas covered by previous models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility, and avoiding the deficiencies of older models by publishing our code open-source in the free-to-use language `Python`. The code of `Riroriro` is also published on Github to allow others to raise issues and propose  new features, and the code split into modularised functions in order to allow for individual sections of the algorithm to be easily edited or substituted if, for example, improvements are made in the underlying science.

# Features

Features of `Riroriro` include:

- Simulating the gravitational waveform signal from a binary merger of two black holes, two neutron stars or a black hole and a neutron star and outputting the data of this signal in terms of frequency and strain amplitude.
- Using such a gravitational wave output, or any other gravitational wave model output in the appropriate form, and given a detector noise spectrum (such spectra are made publicly available by LIGO), calculating the signal-to-noise ratio (SNR) of the signal at a given distance assuming optimal alignment.
- Calculating the horizon distance (maximum distance at which an event could be observed) for a gravitational wave model and a given detector.
- Given the optimal-alignment SNR of an event, evaluating its detectability, the probability that the event would be detected with a SNR above the commonly used threshold of 8, if the alignment would be arbitrary.

In addition, we have created Jupyter Notebook tutorials to help users get started with `Riroriro` ([see tutorials](https://github.com/wvanzeist/riroriro_tutorials)).

# Research

`Riroriro` has been used for research in conjunction with `BPASS`, a a suite of computer programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions [@bpass1, @bpass2], which also has a `Python` interface called `Hoki` [@hoki]. This research took rates of formation of merging systems from `BPASS` and then evaluated the detectability of the gravitational wave signals from those systems using `Riroriro` to obtain predictions of the rates at which gravitational waves of different types would be expected to be observed by the gravitational wave detectors LIGO and Virgo, which can be compared to their actual observations [ref for in-prep mass distribution paper?].

# Acknowledgments

# References
