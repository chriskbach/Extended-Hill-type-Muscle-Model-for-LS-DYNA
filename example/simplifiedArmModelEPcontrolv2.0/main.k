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
*TITLE
$#                                                                         title
Simple Arm Model with EHTM and Equilibrium Point Control
*PARAMETER
$#    name       val$
$-------------------------------------------------------------------------------
$    Define stepsize, endtime and lsdyna standard output
$-------------------------------------------------------------------------------
R  stepsiz    5.0E-5
R   endtim       2.5
R    dtout      0.05
$-------------------------------------------------------------------------------
$    Define controller type
$          0: no controller 
$          1: lambda controller 
$          2: hybrid controller
$-------------------------------------------------------------------------------
R  control         2
$-------------------------------------------------------------------------------
$    Define constant stimulation levels at t = 0
$-------------------------------------------------------------------------------
R     mef1    0.0010                                                                
R     mee1    0.3220                                                               
R     bee1    0.1845                                                                
R     bef1      0.12 
R      bb1    0.1845
$-------------------------------------------------------------------------------
$    Define constant stimulation levels at t = 1.5 $-------------------------------------------------------------------------------
R     mef2    0.3689                                                                
R     mee2   0.30605                                                               
R     bee2   0.29658                                                                
R     bef2       0.4   
R      bb2       0.4
$-------------------------------------------------------------------------------
$    Define constant lambda values at t = 0
$-------------------------------------------------------------------------------
R   lb1mef 9.012E-02
R   lb1mee 6.000E-02
R   lb1bef  6.77E-02
R   lb1bee 7.067E-02
R    lb1bb 7.608E-02
$-------------------------------------------------------------------------------
$    Define constant lambda values at t = 1.5
$-------------------------------------------------------------------------------
R   lb2mef 5.012E-02
R   lb2mee 8.654E-02
R   lb2bef  5.77E-02
R   lb2bee 1.067E-01
R    lb2bb 5.608E-02
$-------------------------------------------------------------------------------
$    Part Shoulder
$-------------------------------------------------------------------------------
*PART
$#                                                                         title
shoulder
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    100000    100000    100001         0         0         0         0         0
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
    100000         2       1.0         2       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
      0.05      0.05      0.05      0.05       0.0       0.0       0.0         0
$-------------------------------------------------------------------------------
$    Part Upper Arm
$-------------------------------------------------------------------------------
*PART_INERTIA
$#                                                                         title
upper arm
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    200000    200000    200001         0         0         0         0         0
$#      xc        yc        zc        tm      ircs    nodeid      
       0.0       0.0       0.0 2.0999999         0    200818
$#     ixx       ixy       ixz       iyy       iyz       izz   
1.00000E-5-5.601E-13-8.1924E-81.00000E-53.08422E-8     0.024
$#     vtx       vty       vtz       vrx       vry       vrz   
       0.0       0.0       0.0       0.0       0.0       0.0
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
    200000         2       1.0         2       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
      0.05      0.05      0.05      0.05       0.0       0.0       0.0         0
$-------------------------------------------------------------------------------
$    Part Lower Arm
$-------------------------------------------------------------------------------
*PART_INERTIA
$#                                                                         title
lower arm
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    300000    300000    300001         0         0         0         0         0
$#      xc        yc        zc        tm      ircs    nodeid      
       0.0       0.0       0.0      1.65         0    300612
$#     ixx       ixy       ixz       iyy       iyz       izz   
1.00000E-5-5.645E-13-8.5340E-81.00000E-53.21300E-8     0.025
$#     vtx       vty       vtz       vrx       vry       vrz   
       0.0       0.0       0.0       0.0       0.0       0.0
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
    300000         2       1.0         2       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
      0.05      0.05      0.05      0.05       0.0       0.0       0.0         0
$-------------------------------------------------------------------------------
$    Parts Muscles 
$-------------------------------------------------------------------------------
*SECTION_BEAM_TITLE
Muskel_Section
$#   secid    elform      shrf   qr/irid       cst     scoor       nsm   
    900000         3       1.0         2         0       0.0       0.0
$#       a     rampt    stress      
      0.05       0.0       0.0
*PART_AVERAGED
$#                                                                         title
MEF
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    900001    900000    900001         0         0         0         0         0
*PART_AVERAGED
$#                                                                         title
MEE
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    900002    900000    900002         0         0         0         0         0
*PART_AVERAGED
$#                                                                         title
BEF
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    900003    900000    900003         0         0         0         0         0
*PART_AVERAGED
$#                                                                         title
BEE
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    900004    900000    900004         0         0         0         0         0
*PART_AVERAGED
$#                                                                         title
BB
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
    900005    900000    900005         0         0         0         0         0
$-------------------------------------------------------------------------------
$    Joint definition
$-------------------------------------------------------------------------------
*CONSTRAINED_JOINT_REVOLUTE_ID
$#     jid                                                                 title
         1Ober-Unterarm
$#      n1        n2        n3        n4        n5        n6       rps      damp
    301319    201679    300263    200335         0         0       1.0     0.001
*CONSTRAINED_JOINT_REVOLUTE_ID
         2Schulter-Oberarm
    100000    200000    301723    201344         0         0       5.0     0.005
$-------------------------------------------------------------------------------
$    Include files
$-------------------------------------------------------------------------------
*INCLUDE
mat.k
*INCLUDE
mesh.k
*INCLUDE
control.k
*INCLUDE
curves.k
*INCLUDE 
muscle_elements.k
*END
