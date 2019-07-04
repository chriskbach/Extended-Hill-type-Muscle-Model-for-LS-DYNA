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
$    Define Stimulation Curves Stim = Stim(t)
$-------------------------------------------------------------------------------
*DEFINE_CURVE_TITLE
Define Timestep
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0            &stepsiz
                 5.0            &stepsiz
*DEFINE_CURVE_TITLE
MEF stimulation
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         6         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0               &mef1
                 1.5               &mef1
                 1.5               &mef2                
                 5.0               &mef2                
*DEFINE_CURVE_TITLE
BEF stimulation
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         7         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0               &bef1
                 1.5               &bef1
                 1.5               &bef2                
                 5.0               &bef2                
*DEFINE_CURVE_TITLE
MEE stimulation
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         8         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0               &mee1
                 1.5               &mee1
                 1.5               &mee2                
                 5.0               &mee2                
*DEFINE_CURVE_TITLE
BEE stimulation
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         9         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0               &bee1
                 1.5               &bee1 
                 1.5               &bee2                
                 5.0               &bee2        
*DEFINE_CURVE_TITLE
BB stimulation
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
        10         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                &bb1
                 1.5                &bb1 
                 1.5                &bb2                
                 5.0                &bb2  
$-------------------------------------------------------------------------------
$    Define target length curves lambda = lambda(t)
$-------------------------------------------------------------------------------                 
*DEFINE_CURVE_TITLE
MEF lambda 
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
        16         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0             &lb1mef
                 1.5             &lb1mef
                 1.5             &lb2mef            
                 5.0             &lb2mef               
*DEFINE_CURVE_TITLE
BEF lambda 
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
        17         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0             &lb1bef
                 1.5             &lb1bef
                 1.5             &lb2bef            
                 5.0             &lb2bef                    
*DEFINE_CURVE_TITLE
MEE lambda 
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
        18         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0             &lb1mee
                 1.5             &lb1mee
                 1.5             &lb2mee           
                 5.0             &lb2mee                
*DEFINE_CURVE_TITLE
BEE lambda 
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
        19         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0             &lb1bee
                 1.5             &lb1bee
                 1.5             &lb2bee            
                 5.0             &lb2bee       
*DEFINE_CURVE_TITLE
BB lambda 
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
        20         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0              &lb1bb
                 1.5              &lb1bb
                 1.5              &lb2bb            
                 5.0              &lb2bb   
*END
