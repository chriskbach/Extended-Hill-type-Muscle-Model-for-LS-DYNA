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
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas      
   &endtim         0       0.0       0.01.000000E8
$
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       0.0       0.9         0       0.0       0.0         1         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     
       0.0         0         0                           0.0
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
    &dtout         0         0         0         0
$#   ioopt     
         0
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
    &dtout         2         0         1       0.0         0
*CONTROL_RIGID
$#     lmf      jntf    orthmd     partm    sparse    metalf    plotel     rbsms
         0         1         0         0         0         0         0         0
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     
         0         0         0       100
*CONTROL_UNITS
$#  length      time      mass      temp    
       0.0       0.0       0.0       0.0
*DATABASE_HISTORY_NODE
$#     id1       id2       id3       id4       id5       id6       id7       id8
    201676    201679   1000021         0         0         0         0         0
*END
