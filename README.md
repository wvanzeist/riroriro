# riroriro

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4588070.svg)](https://doi.org/10.5281/zenodo.4588070)

**Riroriro** is a set of Python modules containing functions to simulate the gravitational waveforms of mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors. Riroriro combines areas covered by previous gravitational wave models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility in Python, a programming language that is ubiquitous in astronomy. Aside from being a research tool, Riroriro is also designed to be easy to use and modify, and it can also be used as an educational tool for students learning about gravitational waves.

The modules “inspiralfuns”, “mergerfirstfuns”, “matchingfuns”, “mergersecondfuns” and “gwexporter”, in that order, can be used to simulate the strain amplitude and frequency of a merger gravitational waveform. The module “snrcalculatorfuns” can compare such a simulated waveform to a detector noise spectrum to calculate a signal-to-noise ratio (SNR) for that signal for that detector. The module “horizondistfuns” calculates the horizon distance of a merger given its waveform, and the module “detectabilityfuns” evaluates the detectability of a merger given its SNR.

Riroriro is installable via pip:

    pip install riroriro

More information on the pip installation can be found here: https://pypi.org/project/riroriro/

Tutorials for Riroriro can be found here: https://github.com/wvanzeist/riroriro_tutorials

Full documentation of each of the functions of Riroriro can be found here: https://wvanzeist.github.io/

Riroriro is one of several Python packages associated with **BPASS** (Binary Population And Spectral Synthesis), a suite of programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions. Each of these associated packages are named after native animals of New Zealand. The riroriro (*Gerygone igata*, also known as the grey warbler) is a small bird that can be recognised by its distinctive melodious call but is rarely seen, similarly to how black hole binary mergers are detected by their gravitational wave signals rather than visually.

The central website of BPASS, which also contains links to related programs, can be found here: https://bpass.auckland.ac.nz
