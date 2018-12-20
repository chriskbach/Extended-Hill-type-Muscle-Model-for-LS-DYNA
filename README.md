# Extended-Hill-type-Muscle-Model-for-LS-DYNA
This is an open-source muscle model for the commercial finite element analysis software LS-DYNA.

For a detailed description of the model see:

> C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt: Implementation and Validation of the Extended Hill-type Muscle Model with Robust Routing Capabilities in LS-DYNA for Active Human Body Models, Biomedical Engineering Online, 2017.
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

If you use this material model for scientific purposes, please cite the original research article.

## Content

1. Source Code
	* commented version of the `dyn21.f`-file
	* patchfile
2. Example
	* isometric contraction of the cat soleus muscle. Run the `main.k`-file.
3. Simulation Data
	* simulation results as shown in the research paper. The units are [N,m,s].

The commented version of the source code is a stripped version of the full `dyn21.f`-file.
You have to copy the source code in the appropriate sections of the `dyn21.f` File provided by your LS-DYNA distributor.
As an alternative you can do this automatically using the patchfile and the well-known unix command [patch](https://linux.die.net/man/1/patch).

In both cases you have to manually add a root-finding algorithm, see line 505 of the commented version of the source code. We used the ZEROIN function from Sec. 7.2 of:

> Forsythe, G.E.; Malcolm, M.A.; Moler, C.B.: Computer Methods for Mathematical Computations. Prentice Hall Professional Technical Reference, 1977.

### LS-DYNA Version

The code was written for LS-DYNA R7.1.2, more precisely the Usermat package:

`ls-dyna_mpp_d_r7_1_2_95028_x64_redhat54_ifort131_sse2_platformmpi.usermat.tar.gz`

was used. You can also check the beginning of the `dyn21.f`-file, which should be

> c
> c $Id: dyn21.F 94310 2014-12-05 01:04:21Z tsay $
> c

With modifications it should compile also with other versions of LS-DYNA.

### Changelog 

[v2.0] - 2019-01-01 
* total muscle length is set to part length (enables routing via part_averaged)
* added muscle-level controllers (EQ-point theroy and reflex)

[v1.0] - 2017-07-12
* initial muscle model implementation in LS-Dyna 
* find detailed description in C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt: Implementation and Validation of the Extended Hill-type Muscle Model with Robust Routing Capabilities in LS-DYNA for Active Human Body Models, Biomedical Engineering Online, 2017.
* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.826209.svg)](https://doi.org/10.5281/zenodo.826209)
