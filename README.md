# Extended-Hill-type-Muscle-Model-for-LS-DYNA

The repository has moved to DaRUS - the Data Repository of the University of Stuttgart: [![DOI](https://img.shields.io/badge/DOI-10.18419%2Fdarus--1144-red.svg)](https://doi.org/10.18419/darus-1144)

### Changelog 

[moved repo] - 2022-04-25
* moved to repository to [![DOI](https://img.shields.io/badge/DOI-10.18419%2Fdarus--1144-red.svg)](https://doi.org/10.18419/darus-1144)

[v2.0] - 2019-05-20
* total muscle length is set to the part length (enables routing via part_averaged)
* added internal muscle-level controller (EQ-point theory and reflex)
* ringbuffer using history variables
* damping methods removed
* improved numerical stability, by allowing negative l_ce for extreme slack in the muscles
* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3142387.svg)](https://doi.org/10.5281/zenodo.3142387)

[v1.0] - 2017-07-12
* initial muscle model implementation in LS-Dyna 
* find detailed description in C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt: Implementation and Validation of the Extended Hill-type Muscle Model with Robust Routing Capabilities in LS-DYNA for Active Human Body Models, Biomedical Engineering Online, 2017.
* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.826209.svg)](https://doi.org/10.5281/zenodo.826209)
