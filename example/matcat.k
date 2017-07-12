$
$ This is supplementary material for the publication:
$
$ C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt:
$ Implementation and Validation of the Extended Hill-type Muscle Model with
$ Robust Routing Capabilities in LS-DYNA for Active Human Body Models,
$ Biomedical Engineering Online, 2017.
$
$ If you use this material model for scientific purposes, please cite
$ the original research article.
$
*KEYWORD
*MAT_USER_DEFINED_MATERIAL_MODELS
$ Hill-type Muscle for a cat soleus moerl12
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
         1  1.000E-6        41        28        14         0         3         4
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0         0         0         0                    
$#  ActOpt   STIM_ID        q0    tauq/c betaq/eta         k         m
         2         3    1.0E-4  1.373E-4    5.27E4       2.9     22.54
$#    Fmax    lCEopt     dWdes   nuCEdes     dWasc   nuCEasc     Arel0     Brel0      
      10.0     0.053      0.35       1.5      0.35       3.0      0.07       0.2
$#    Secc      Fecc     LPEE0     nuPEE      FPEE     lSEE0  dUSEEnll    duSEEl
       2.0       1.5       0.9       2.5       2.0     0.060    0.0425     0.017
$#  dFSEE0   Damping    Param1    Param2
       4.0       3.0       0.3      0.01    
*END
