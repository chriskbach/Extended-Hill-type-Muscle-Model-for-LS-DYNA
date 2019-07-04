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
*Parameter
$#   prmr1      val1     
$-------------------------------------------------------------------------------
$    Define Position of Routing Points
$-------------------------------------------------------------------------------
R p_posBBx  0.06487
R p_posBBy  0.11645
R i_posBBx  0.2369 
R i_posBBy  0.2889 
R p_posMFx  0.21379
R p_posMFy  0.37498
R i_posMFx  0.18599
R i_posMFy  0.18599
R p_posMEx  0.26201
R p_posMEy  0.23017
R i_posMEx  0.25942
R i_posMEy  0.2499 
R p_posBFx  0.20065
R p_posBFy  0.22   
R i_posBFx  0.2369 
R i_posBFy  0.2619 
R p_posBEx  0.26201
R p_posBEy  0.23017
R i_posBEx  0.25942
R i_posBEy  0.2499 

R    x_MEF    0.2260
$-------------------------------------------------------------------------------
$    Muscle beam element definition
$-------------------------------------------------------------------------------
$Muscle Beams
*ELEMENT_BEAM
$#   eid     pid      n1      n2      n3     rt1     rr1     rt2     rr2   local
$MEF muscle beams
  900001  900001 1000004 1500000 
  900011  900001 1001173 1000003
  900012  900001 1000003 1000004  
$MEE muscle beams
  900002  900002 2500001 2000001 
  900021  900002 2000001 2000002      
  900022  900002 2000002 2000003
  900023  900002 2000003 2500000
$BRAC muscle beams
  900003  900003 3500003 3000004
  900031  900003 3000004 3000005
  900032  900003 3000005 3500010
$BEE muscle beams 
  900004  900004 1001101 4000000
  900041  900004 4000000 4000001      
  900042  900004 4000001 4000002     
  900044  900004 4000002 4500000
$BB muscle beams
  900005  900005 5500010 5000001
  900051  900005 5000001 5000005     
  900052  900005 5000005  301546 
$-------------------------------------------------------------------------------
$    Attachment and Routing Nodes as Constrained Extra Nodes 
$-------------------------------------------------------------------------------
*CONSTRAINED_EXTRA_NODES_NODE
$     pid        nid
$MEF
   300000    1500010
   300000    1000003
   300000    1000004
   200000    1000005
   200000    1500000
$MEE
   200000    2500001
   200000    2000002
   200000    2000003
   300000    2500000
$BRAC 
   200000    3500003
   200000    3000005
   200000    3000004
   300000    3500010
$BEE   
   200000    4000000
   200000    4000001
   200000    4000002
   300000    4500000
$BB
   200000    5000001
   300000    5500010
   200000    5000005
$-------------------------------------------------------------------------------
$    Attachment and Routing Node Definition
$-------------------------------------------------------------------------------
*NODE
$$ ---------------------- MEF -------------------------
$BeamNodes Cable
 1000003       &p_posMFx       &p_posMFy             0.0       0       0
$Insertion Node 
 1500000       &i_posMFx       &i_posMFy             0.0       0       0
 1500010          &x_MEF          0.4789             0.0       0       0
$
 1000004        0.211473        0.359231             0.0       0       0
 1000005         0.19294        0.233238             0.0       0       0
$$ ---------------------- MEE -------------------------
$BeamNodes Cable
 2000001          0.1800          0.1500             0.0       0       0
 2000002          0.2410          0.2110             0.0       0       0
 2000003       &p_posMEx       &p_posMEy             0.0       0       0
$Origin
 2500001          0.0800          0.0500             0.0       0       0
$Insertion Node
 2500000       &i_posMEx       &i_posMEy             0.0       0       0       
$$ ---------------------- BB -------------------------
$BeamNodes Cable
 5000005       &p_posBBx       &p_posBBy             0.0       0       0
 5000001       0.122213         0.173933             0.0       0       0
$Insertion Nodes
 5500010       &i_posBBx       &i_posBBy             0.0       0       0  
$
$$ ---------------------- BEF -------------------------
$BeamNodes Cable
 3500003          0.1100          0.1300             0.0       0       0
 3000004          0.1700          0.1900             0.0       0       0
 3000005       &p_posBFx       &p_posBFy             0.0       0       0
$Insertion Nodes
 3500010       &i_posBFx       &i_posBFy             0.0       0       0  
$
$$ ---------------------- BEE -------------------------
$BeamNodes Cable
 4000000          0.1355          0.1135             0.0       0       0
 4000001           0.241           0.217             0.0       0       0
 4000002       &p_posBEx       &p_posBEy             0.0       0       0
$Insertion Node
 4500000       &i_posBEx       &i_posBEy             0.0       0       0   
*END
