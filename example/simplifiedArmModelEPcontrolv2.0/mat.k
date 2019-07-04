$
$ This is supplementary material for the publications:
$
$ O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr: 
$ Development of an internal physiological muscle controller within an 
$ open‐source Hill‐type material model in LS‐DYNA, Proceedings in Applied 
$ Mathematics and Mechanics, Munich, 2018.
$
$ O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr: 
$ Integrated Physiologically Motivated Controller for the Open-Source Extended 
$ Hill-type Muscle Model in LS-DYNA. Proceedings of IRCOBI Conference, Athens, 2018.
$
$ If you use this material model for scientific purposes, please cite
$ the original research articles.
$
*KEYWORD
$-------------------------------------------------------------------------------
$    Muscle Material Cards
$-------------------------------------------------------------------------------
*MAT_USER_DEFINED_MATERIAL_MODELS
$ MEF muscle material
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
    9000011.00000E-6        41        40      1200         0        31        32
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0         0         0         0                    
$#  ActOpt   STIM_ID        q0    tauq/c betaq/eta         k         m   lOffset
       2.0       6.0     0.0051.37000E-4   52700.0 2.9000001      11.3    
$#    Fmax    lCEopt     dWdes   nuCEdes     dWasc   nuCEasc     Arel0     Brel0 
     420.0     0.102     0.390       1.5     0.390       3.0       0.5       5.0
$#    Secc      Fecc     LPEE0     nuPEE      FPEE     lSEE0  dUSEEnll    duSEEl
       2.0       1.5      0.95       2.5       2.0    0.1899    0.0425     0.017
$#  dFSEE0   Damping    Param1    Param2    Output     dtOut     IBULK        IG
     568.0       3.00.30000001      0.01         2      0.10     0.0051.37000E-4
$#  Regler    Lambda        kp        kd     delay    Switch    thresh   
  &control        16         3         1         0        10      0.05     
*MAT_USER_DEFINED_MATERIAL_MODELS
$ MEE muscle material
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
    9000021.00000E-6        41        40      1200         0        31        32
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0         0         0         0                    
$#  ActOpt   STIM_ID        q0    tauq/c betaq/eta         k         m   lOffset
       2.0       8.0     0.0051.37000E-4   52700.0 2.9000001      11.3    
$#    Fmax    lCEopt     dWdes   nuCEdes     dWasc   nuCEasc     Arel0     Brel0 
    1550.0     0.102     0.390       1.5     0.390       3.0       0.5       5.0
$#    Secc      Fecc     LPEE0     nuPEE      FPEE     lSEE0  dUSEEnll    duSEEl
       2.0       1.5      0.95       2.5       2.0    0.2048    0.0425     0.017
$#  dFSEE0   Damping    Param1    Param2    Output     dtOut     IBULK        IG
     620.0       3.00.30000001      0.01         2      0.10     0.0051.37000E-4
$#  Regler    Lambda        kp        kd     delay    Switch    thresh   
  &control        18         3         1         0        10      0.05 
*MAT_USER_DEFINED_MATERIAL_MODELS
$ BEF muscle material
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
    9000031.00000E-6        41        40      1200         0        31        32
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0         0         0         0                    
$#  ActOpt   STIM_ID        q0    tauq/c betaq/eta         k         m   lOffset
       2.0       7.0     0.0051.37000E-4   52700.0 2.9000001      11.3    
$#    Fmax    lCEopt     dWdes   nuCEdes     dWasc   nuCEasc     Arel0     Brel0 
     420.0     0.051     0.390       1.5     0.390       3.0       0.5       5.0
$#    Secc      Fecc     LPEE0     nuPEE      FPEE     lSEE0  dUSEEnll    duSEEl
       2.0       1.5      0.95       2.5       2.0    0.0948    0.0425     0.017
$#  dFSEE0   Damping    Param1    Param2    Output     dtOut     IBULK        IG
 165.60001       3.00.30000001      0.01         2      0.10     0.0051.37000E-4
$#  Regler    Lambda        kp        kd     delay    Switch    thresh   
  &control        17         3         1         0        10      0.05 
*MAT_USER_DEFINED_MATERIAL_MODELS
$ BEE muscle material
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
    9000041.00000E-6        41        40      1200         0        31        32
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0         0         0         0                    
$#  ActOpt   STIM_ID        q0    tauq/c betaq/eta         k         m   lOffset
       2.0       9.0     0.0051.37000E-4   52700.0 2.9000001      11.3    
$#    Fmax    lCEopt     dWdes   nuCEdes     dWasc   nuCEasc     Arel0     Brel0 
     603.0     0.146     0.390       1.5     0.390       3.0       0.5       5.0
$#    Secc      Fecc     LPEE0     nuPEE      FPEE     lSEE0  dUSEEnll    duSEEl
       2.0       1.5      0.95       2.5       2.0    0.2493   0.0425     0.017
$#  dFSEE0   Damping    Param1    Param2    Output     dtOut     IBULK        IG
     241.2       3.00.30000001      0.01         2      0.10     0.0051.37000E-4
$#  Regler    Lambda        kp        kd     delay    Switch    thresh   
  &control        19         3         1         0        10      0.05 
*MAT_USER_DEFINED_MATERIAL_MODELS
$ BB muscle material
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
    9000051.00000E-6        41        40      1200         0        31        32
$#   ivect     ifail    itherm    ihyper      ieos      lmca    unused    unused
         0         0         0         0         0         0                    
$#  ActOpt   STIM_ID        q0    tauq/c betaq/eta         k         m   lOffset
       2.0      10.0     0.0051.37000E-4   52700.0 2.9000001      11.3     
$#    Fmax    lCEopt     dWdes   nuCEdes     dWasc   nuCEasc     Arel0     Brel0 
     414.0     0.168     0.390       1.5     0.390       3.0       0.5       5.0
$#    Secc      Fecc     LPEE0     nuPEE      FPEE     lSEE0  dUSEEnll    duSEEl
       2.0       1.5      0.95       2.5       2.0    0.2863    0.0425     0.017
$#  dFSEE0   Damping    Param1    Param2    Output     dtOut     IBULK        IG
     241.2       3.00.30000001      0.01         2      0.10     0.0051.37000E-4
$#  Regler    Lambda        kp        kd     delay    Switch    thresh   
  &control        20         3         1         0        10      0.05 
$-------------------------------------------------------------------------------
$    Other material definitions
$-------------------------------------------------------------------------------
*MAT_RIGID
$#     mid        ro         e        pr         n    couple         m     alias
    100001    1500.01.60000E100.30000001       0.0       0.0       0.0          
$#     cmo      con1      con2    
       1.0         7         7
$#lco or a1        a2        a3        v1        v2        v3  
       0.0       0.0       0.0       0.0       0.0       0.0
*MAT_RIGID
$#     mid        ro         e        pr         n    couple         m     alias
    200001    1500.01.60000E100.30000001       0.0       0.0       0.0          
$#     cmo      con1      con2    
       1.0         7         7
$#lco or a1        a2        a3        v1        v2        v3  
       0.0       0.0       0.0       0.0       0.0       0.0
*MAT_RIGID
$#     mid        ro         e        pr         n    couple         m     alias
    300001    1500.01.60000E100.30000001       0.0       0.0       0.0          
$#     cmo      con1      con2    
       1.0         3         4
$#lco ora1        a2        a3        v1        v2        v3  
       0.0       0.0       0.0       0.0       0.0       0.0
*END
