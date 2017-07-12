$
$ This is supplementary material for the publication:
$
$ C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr, S. Schmitt:
$ Implementation and Validation of the Extended Hill-type Muscle Model with
$ Robust Routing Capabilities in LS-DYNA for Active Human Body Models,
$ Biomedical Engineering Online, 2017.
$
$ If you use this material model for scientific purposes, please cite the original
$ research article.
$
*KEYWORD
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R dtout   0.001                                                                 
*TITLE
$#                                                                         title
isometric Muscle contraction
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas      
       1.5         0       0.0       0.01.000000E8
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       0.00.49999998         0       0.0       0.0         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     
       0.0         0         0                           0.0
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
$                                                                              $
$                             DATABASE (OUTPUT)                                $
$                                                                              $
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt     
&dtout             0         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt     
&dtout             0         0         1
*DATABASE_RCFORC
$#      dt    binary      lcur     ioopt     
&dtout             0         0         1
*DATABASE_SPCFORC
$#      dt    binary      lcur     ioopt     
&dtout             0         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
     0.001         0         0         0         0
$#   ioopt     
         0
*DATABASE_HISTORY_NODE
$#     id1       id2       id3       id4       id5       id6       id7       id8
         1         2         0         0         0         0         0         0
*PART
$#                                                                         title
Beam2
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2         1         1         0         0         0         0         0
*SECTION_BEAM
$#   secid    elform      shrf   qr/irid       cst     scoor       nsm   
         1         3       1.0         2         0       0.0       0.0
$#       a     rampt    stress      
       1.0       0.0       0.0
*INCLUDE
matcat.k
$
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         3         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
                 0.1                 0.0
                 0.1                 1.0
                 1.1                 1.0
                 1.1                 0.0
                 1.5                 0.0
*ELEMENT_BEAM
$#   eid     pid      n1      n2      n3     rt1     rr1     rt2     rr2   local
       1       2       1       2       0       0       0       0       0       2
*END
