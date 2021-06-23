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

## Paper

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02968/status.svg)](https://doi.org/10.21105/joss.02968)

A paper describing Riroriro has been published in the Journal of Open Source Software. If you use Riroriro in your work, please cite this paper! https://doi.org/10.21105/joss.02968

    @ARTICLE{2021JOSS....6.2968V,
           author = {{van Zeist}, Wouter G.~J. and {Stevance}, H{\'e}lo{\"i}se F. and {Eldridge}, J.~J.},
            title = "{Riroriro: Simulating gravitational waves and evaluating their detectability in Python}",
          journal = {The Journal of Open Source Software},
         keywords = {Python, neutron stars, astronomy, gravitational waves, black holes, General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
             year = 2021,
            month = mar,
           volume = {6},
           number = {59},
              eid = {2968},
            pages = {2968},
              doi = {10.21105/joss.02968},
    archivePrefix = {arXiv},
           eprint = {2103.06943},
     primaryClass = {gr-qc},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2021JOSS....6.2968V},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
