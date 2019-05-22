# Extended-Hill-type-Muscle-Model-for-LS-DYNA
This is an open-source muscle model for the commercial finite element analysis software LS-DYNA.

For a detailed description of the model see:

> C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt: Implementation and Validation of the Extended Hill-type Muscle Model with Robust Routing Capabilities in LS-DYNA for Active Human Body Models, Biomedical Engineering Online, 2017.

> O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr: Development of an internal physiological muscle controller within an open‐source Hill‐type material model in LS‐DYNA, Proceedings in Applied Mathematics and Mechanics, Munich, 2018.

> O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr: Integrated Physiologically Motivated Controller for the Open-Source Extended Hill-type Muscle Model in LS-DYNA. Proceedings of IRCOBI Conference, Athens, 2018.

> 
> **Abstract**:
> 
> An extended four element Hill-type muscle model with serial damping, eccentric force-velocity relation, integrated Ca2+ dependent activation dynamics and routing capability is implemented into the general-purpose finite element (FE) simulation software LS-DYNA as a user material for truss elements.
This material model is verified and validated with three different sets of mammalian experimental data, taken from the literature.
It is compared to the *MAT_MUSCLE (*MAT_156) Hill-type muscle model already existing in LS-DYNA , which is currently used in finite element Human Body Models (HBMs).
An application example with an arm model extracted from the FE ViVA OpenHBM is given, taking into account physiological muscles paths.
The simulation results show better material model accuracy, calculation robustness and improved muscle routing capability compared to *MAT_156.
The FORTRAN source code for the user material subroutine `dyn21.f` and the muscle parameters for all simulations, conducted in the study, are given as supplementary material under a open source license.
This enables a quick application of the proposed material model in LS-DYNA, especially in Active Human Body Models (AHBMs).

If you use this material model for scientific purposes, please cite the original research articles.

## Content

1. Source Code
	* commented version of the `dyn21.f`-file
	* patchfile
2. Examples for Release v2.0 will be available soon
	* head-fall test setup of the ViVA HBM with neck musculature using the Extended-Hill-type-Muscle-Model. Run the `main.k`-file.

The commented version of the source code is a stripped version of the full `dyn21.f`-file.
You have to copy the source code in the appropriate sections of the `dyn21.f` File provided by your LS-DYNA distributor.
As an alternative you can do this automatically using the patchfile and the well-known unix command [patch](https://linux.die.net/man/1/patch).

In both cases you have to manually add a root-finding algorithm, see line 505 of the commented version of the source code. We used the ZEROIN function from Sec. 7.2 of:

> Forsythe, G.E.; Malcolm, M.A.; Moler, C.B.: Computer Methods for Mathematical Computations. Prentice Hall Professional Technical Reference, 1977.

### LS-DYNA Version

The code was written for LS-DYNA R8.1.0, more precisely the Usermat package:

`ls-dyna_smp_d_r810_x64_redhat511_ifort131`

was used. You can also check the beginning of the `dyn21.f`-file, which should be

> c
c $Id: dyn21.F 96707 2015-03-18 22:29:10Z ubasu $
> c

With modifications it should compile also with other versions of LS-DYNA.

### Changelog 

[v2.0] - 2019-05-20
* total muscle length is set to the part length (enables routing via part_averaged)
* added internal muscle-level controller (EQ-point theory and reflex)
* ringbuffer using history variables
* damping methods removed
* improved numerical stability, by allowing negative l_ce for extreme slack in the muscles

[v1.0] - 2017-07-12
* initial muscle model implementation in LS-Dyna 
* find detailed description in C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt: Implementation and Validation of the Extended Hill-type Muscle Model with Robust Routing Capabilities in LS-DYNA for Active Human Body Models, Biomedical Engineering Online, 2017.
* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.826209.svg)](https://doi.org/10.5281/zenodo.826209)
