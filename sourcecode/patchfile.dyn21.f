2270,2271c2270,2273
<    41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
<      1   temper,failel,crv,a(lcma),0,elsiz,idele)
---
>    41   call umat41 (lft,cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
>      1   temper,failel,crv,a(lcma),0,elsiz,idele,r_mem(dm_x),
>      2   r_mem(dm_v),r_mem(dm_a),i,d1,d2,d3,NHISVAR,
>      3   no_hsvs)
2414,2415c2416,2469
<       subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
<      1 temper,failel,crv,cma,qmat,elsiz,idele)
---
> C
> C
> C % If you use this material model for scientific purposes, please cite
> C % the original research articles:
> C % 1) C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr,
> C % S. Schmitt: Implementation and Validation of the Extended Hill-type
> C % Muscle Model with Robust Routing Capabilities in LS-DYNA for Active
> C % Human Body Models, Biomedical Engineering Online, 2017.
> C % 
> C % 2) O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr:
> C % Development of an internal physiological muscle controller within an 
> C % open‐source Hill‐type material model in LS‐DYNA, Proceedings in Applied 
> C % Mathematics and Mechanics, Munich, 2018.
> C % 
> C % 3) O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr:
> C % Integrated Physiologically Motivated Controller for the Open-Source 
> C % Extended Hill-type Muscle Model in LS-DYNA. Proceedings of IRCOBI 
> C % Conference, Athens, 2018.
> C % 
> C % Copyright (c) 2019 on modifications and LS-DYNA implementation belongs to
> C % O. Martynenko, F. Kempter, C. Kleinbach, S. Schmitt and J. Fehr.
> C % Copyright (c) 2017 on modifications and LS-DYNA implementation belongs to
> C % C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr and
> C % S. Schmitt. Original Copyright (c) 2014 belongs to D. Haeufle, M. Guenther,
> C % A. Bayer and S. Schmitt. 
> C % All rights reserved. 
> C % Redistribution and use in source and binary forms, with or without
> C % modification, are permitted provided that the following conditions are
> C % met:
> C %
> C %  1 Redistributions of source code must retain the above copyright notice,
> C %    this list of conditions and the following disclaimer. 
> C %  2 Redistributions in binary form must reproduce the above copyright
> C %    notice, this list of conditions and the following disclaimer in the
> C %    documentation and/or other materials provided with the distribution.
> C %  3 Neither the name of the owner nor the names of its contributors may be
> C %    used to endorse or promote products derived from this software without
> C %    specific prior written permission.
> C %
> C % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
> C % IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
> C % THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
> C % PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE
> C % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
> C % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
> C % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
> C % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
> C % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
> C % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
> C % THE POSSIBILITY OF SUCH DAMAGE.
> 
>       subroutine umat41 (lft,cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
>      1 temper,failel,crv,cma,qmat,elsiz,idele,x,v,a,i,d1,d2,d3,NHISVAR,
>      2 no_hsvs)
2424,2499c2478,2604
< c     isotropic elastic material (sample user subroutine)
< c
< c     Variables
< c
< c     cm(1)=first material constant, here young's modulus
< c     cm(2)=second material constant, here poisson's ratio
< c        .
< c        .
< c        .
< c     cm(n)=nth material constant
< c
< c     eps(1)=local x  strain increment
< c     eps(2)=local y  strain increment
< c     eps(3)=local z  strain increment
< c     eps(4)=local xy strain increment
< c     eps(5)=local yz strain increment
< c     eps(6)=local zx strain increment
< c
< c     sig(1)=local x  stress
< c     sig(2)=local y  stress
< c     sig(3)=local z  stress
< c     sig(4)=local xy stress
< c     sig(5)=local yz stress
< c     sig(6)=local zx stress
< c
< c     hsv(1)=1st history variable
< c     hsv(2)=2nd history variable
< c        .
< c        .
< c        .
< c        .
< c     hsv(n)=nth history variable
< c
< c     dt1=current time step size
< c     capa=reduction factor for transverse shear
< c     etype:
< c        eq."solid" for solid elements
< c        eq."sph" for smoothed particle hydrodynamics
< c        eq."sld2d" for shell forms 13 (2D solids - plane strain)
< c        eq."sldax" for shell forms 14, and 15 (2D solids - axisymmetric)
< c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
< c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
< c        eq."tshel" for thick shell forms 3 and 5
< c        eq."hbeam" for beam element forms 1 and 11
< c        eq."tbeam" for beam element form 3 (truss)
< c        eq."dbeam" for beam element form 6 (discrete)
< c        eq."beam " for all other beam elements
< c
< c     tt=current problem time.
< c
< c     temper=current temperature
< c
< c     failel=flag for failure, set to .true. to fail an integration point,
< c            if .true. on input the integration point has failed earlier
< c
< c     crv=array representation of curves in keyword deck
< c
< c     cma=additional memory for material data defined by LMCA at 
< c       6th field of 2nd crad of *DATA_USER_DEFINED
< c
< c     elsiz=characteristic element size
< c
< c     idele=element id
< c
< c     All transformations into the element local system are
< c     performed prior to entering this subroutine.  Transformations
< c     back to the global system are performed after exiting this
< c     routine.
< c
< c     All history variables are initialized to zero in the input
< c     phase. Initialization of history variables to nonzero values
< c     may be done during the first call to this subroutine for each
< c     element.
< c
< c     Energy calculations for the dyna3d energy balance are done
< c     outside this subroutine.
---
> c     Extended Hill-type muscle model with a contractile element
> c     a parallel elastic element and serial elastic as well as a
> c     serial damping element
> c
> c     The theory for this muscle model is described in:
> c     C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr,
> c     S. Schmitt: Implementation and Validation of the Extended Hill-type
> c     Muscle Model with Robust Routing Capabilities in LS-DYNA for Active
> c     Human Body Models, Biomedical Engineering Online, 2017.
> c
> c     Variables (activation dynamic)
> c     cm(1)=Activation Option (EQ.0. Activation Values see STIM ID, EQ.1 Zajac, EQ.2 Hatze)
> c     cm(2)=GT0: LE.0.0 constant values for STIM or Activation GT.0.0 curve id for STIM or Activation
> c     cm(3)=q0 minimum value of q
> c     cm(4)=tau_q time constant of rising activation LT0 curve id of tau_q over time
> c     cm(4)=c
> c     cm(5)=beta_q ratio between tau_q and time constant of falling activation
> c     cm(5)=eta
> c     cm(6)=k
> c     cm(7)=m
> c     cm(8)=muscle length offset (not needed, use instead *Part_Averaged for routing)
> c
> c     Variables (isometric force)
> c
> c     cm(9)=F_max maximum isometric force
> c     cm(10)=l_CE_opt optimal fibre length
> c     cm(11)=dW_des
> c     cm(12)=nu_CE_des
> c     cm(13)=dW_asc
> c     cm(14)=nu_CE_asc
> c
> c     Variables (Hill Parameter, concentric)
> c
> c     cm(15)=A_rel_0
> c     cm(16)=B_rel_0
> c
> c     Variables (van Soest Parameter, eccentric)
> c
> c     cm(17)=S_ecc
> c     cm(18)=F_ecc
> c
> c     Variables (Parallel elastic element)
> c
> c     cm(19)=L_PEE_0
> c     cm(20)=nu_PEE
> c     cm(21)=F_PEE
> c
> c     Variables (Seriell elastic element)
> c
> c     cm(22)=l_SEE_0
> c     cm(23)=dU_SEE_nllfindkParams
> c     cm(24)=dU_SEE_l
> c     cm(25)=dF_SEE_0
> c
> c     Variables (Damping element)
> c
> c     cm(26)= --- (former damping method)
> c     cm(27)=D_SE
> c     cm(28)=R_SE
> c
> c     Variables (Output definition; musout.(partid))
> c
> c     cm(29)=output method (EQ.0.  no output 
> c                           EQ.1.  basic output (idpart, tt, hsv(2:10)) 
> c                           EQ.2.  advanced output (basic output plus dot_l_CE, dot_l_MTC, lCEdelay, dotlCEdelay))
> c
> c     cm(30)=timestep of outputfile
> c
> c     Variables for Controller
> c
> c     cm(33)=Activation Method (EQ.1. lambda_controller
> c                               EQ.2. hybrid_controller
> c                               EQ.3. reflexive controller)
> c
> c     cm(34)=target l_CE
> c     cm(35)=kp
> c     cm(36)=kd
> c     cm(37)=delay of lCEdelay / dotlCEdelay
> c     cm(38)=time till swap from alpha to lambda / t_PreSim for reflexive controller
> c     cm(39)=threshold for reflex controller (e.g. 0.10 for a 10% strain threshold)
> c
> c [...]
> c
> c     d1  - strain rate/increment in x  direction, local x for shells
> c     d2  - strain rate/increment in y  direction, local y for shells
> c     d3  - strain rate/increment in z  direction, local z for shells
> c ------------------------------------------------------------
> c ------ history variables - overview ------------------------
> c ------------------------------------------------------------
> c         hsv(1)    = sig(1)
> c         hsv(2)    = STIM
> c         hsv(3)    = q
> c         hsv(4)    = F_MTC
> c         hsv(5)    = F_CE
> c         hsv(6)    = F_PEE
> c         hsv(7)    = F_SEE
> c         hsv(8)    = F_SDE
> c         hsv(9)    = elleng (l_MTC)
> c         hsv(10)   = l_CE
> c         hsv(11)   = dot_l_MTC
> c         hsv(12)   = dot_l_CE
> c         hsv(13)   = counter_output
> c         hsv(14)   = F_isom
> c         hsv(15)   = gam_rel
> c         hsv(16)   = l_MTC_0 (initial element length)
> c ------------------------------------------------------------
> c ------ only needed if controller are used ------------------
> c ------------------------------------------------------------
> c         hsv(20)   = lCEdelay
> c         hsv(21)   = dotlCEdelay
> c         hsv(22)   = l_CE_ref (for reflex controller)
> c         hsv(23)   = strain (for reflex controller)
> c         hsv(24)   = STIM_reflex_prev (STIM value of reflex controller in previous timestep)
> c ------------------------------------------------------------
> c ------ delay buffer-----------------------------------------
> c ------------------------------------------------------------
> c         hsv(30)   = buffersize
> c         hsv(31)   = idx_begin_lCEbuffer (hsv(142:145) seems to be used internally by lsdyna, suggestion: use hsv(31)=150)
> c         hsv(32)   = indexr1...index1 of lce-ringbuffer (not eq index of hsv)
> c         hsv(33)   = indexr1...index2 of lce-ringbuffer (not eq index of hsv)
> c         hsv(34)   = lasttime
> c         hsv(36)   = begindotlCE
> c         hsv(37)   = indexdotr1...index1 of dotlce-ringbuffer (not eq index of hsv)    
> c         hsv(38)   = indexdotr2...index2 of dotlce-ringbuffer (not eq index of hsv)  
> c         hsv(39)   = dotlasttime   
> c         hsv(hsv(31):hsv(31)+hsv(30)-1) = ringbuffer_l_CE
> c         hsv(hsv(36):hsv(36)+hsv(30)-1) = ringbuffer_dot_l_CE
2500a2606
> c [...]
2504c2610,2624
<       dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
---
> c
>       common/prescloc/voltot(nlq)
>       common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
>      + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
>      + grvity,idirgv,nodspc,nspcor,nusa
> c$omp threadprivate (/aux33loc/)
> c$omp threadprivate (/prescloc/)
> c
> c############################################
> c#   initializing the variables             #
> c############################################
> c
>       dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
>       dimension x(3,*),v(3,*),a(3,*),qmat(3,3)
>       dimension d1(*),d2(*),d3(*)
2507,2513c2627,2856
<       integer idele
< c
<       if (ncycle.eq.1) then
<         if (cm(16).ne.1234567) then
<           call usermsg('mat41')
<         endif
<       endif
---
>       real*8 l_PEE0, K_PEE, d_SE_max, l_SEE_nll, v_SEE, K_SEE_nl
>       real*8 F_isom, F_PEE, F_SEE, l_SEE, L_A_rel, A_rel, q,dq
>       real*8 B_rel, L_B_rel, Q_Brel, D0, C2, C1, C0, dot_l_CE
>       real*8 O1, O2, O3, F_CE, F_CE_init, F_SDE, F_MTC
>       real*8 dot_l_MTC,Q_Arel, K_SEE_l, F_SUM, AREA,l_CE,q_old
>       real*8 tol,l_MTC_0,STIM,dSTIM,tau,dtau,epsilon
>       real*8 delay, gam_rel, rho_act, lCEdelay, dotlCEdelay
>       real*8 lambda,dlambda
>       integer curve_id, var,k,j,idpart
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
>       common/soundloc/sndspd(nlq),sndsp(nlq),diagm(nlq),sarea(nlq)
>       common/aux33loc/
>      1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
>       common/aux14loc/
>      1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
>      2 sig5(nlq),sig6(nlq),epx1(nlq),epx2(nlq),aux(nlq,14),dig1(nlq),
>      3 sig_pass(nlq),sig_actv(nlq),act_levl(nlq),out_leng(nlq),
>      4 eps_rate(nlq),sig_svs(nlq),sig_sde(nlq)
> c############################################
> c#   calc muscle parameters                 #
> c############################################
>        l_PEE0=cm(19)*cm(10)
>        d_SE_max=cm(27)*(cm(9)*cm(15))/(cm(10)*cm(16))
>        l_SEE_nll=(1.0+cm(23))*cm(22)
>        v_SEE=cm(23)/cm(24)
>        K_SEE_nl=cm(25)/(cm(23)*cm(22))**v_SEE
>        K_SEE_l=cm(25)/(cm(24)*cm(22))
>        K_PEE=cm(21)*(cm(9)/(cm(10)*(cm(11)+1.0-cm(19)))**cm(20))
> c############################################
> c#   calc length of muscle element          #
> c############################################ 
>        idpart=lqfmiv(mxt(lft))
> c muscle length is part length plus offset
>        elleng=sqrt(sarea(i))+cm(8)
>        if (ncycle.LT.1) then 
>         hsv(16)=elleng
>         return
>        endif       
> c     
>        if (ncycle.LT.2) then
> c       
>          epsilon = 1.0
>    99    epsilon = 0.5*epsilon
>          if ( 1.0 + 0.5*epsilon .GT. 1.0 ) GOTO 99
> c
>          hsv(3)=cm(3)
>          hsv(13)=0.0
>          hsv(15)=0.0
>          hsv(16)=elleng
> c     
>        end if
> c     
> c############################################
> c#   calc delayed lCE/dotlCE                #
> c############################################       
> c ncycle.EQ.3, because for ncycle<=2: dt1==0 in case of part_averaged
>        if (ncycle.LT.3.and.cm(33).ge.1.0.and.cm(37).GT.0.0) then
> c             delay = cm(37)
> c beginlCE (150 because of lsdyna messing with our data at lower index)
>               hsv(31) = 150
> c total buffersize based on NHISVAR  (CEILING = "abrunden")         
>               hsv(30) = ceiling((no_hsvs-hsv(31)-4)/2);
> c NHISVAR is the amount of maximum hsv (defined in compile process)
> c no_hsvs is the user defined amount of hsvs (!<=NHISVAR)
> c
> c indexr1  
>               hsv(32) = 0
> c indexr2
>               hsv(33) = 1
> c lasttime
>               hsv(34) = 0
> c begindotlCE
>               hsv(36) = hsv(31)+hsv(30)+2
> c indexdotr1              
>               hsv(37) = 0
> c indexdotr1
>               hsv(38) = 1
> c dotlasttime
>               hsv(39) = 0
>        end if
> c
>        if (ncycle.GE.3.and.cm(33).ge.1.0.and.cm(37).GT.0.0) then
> c all controllers with delay.gt.0
>             lCEdelay = delaylCE(ncycle,hsv,cm,idpart,tt)
>             dotlCEdelay=hsv(12)
>             if (cm(33).le.2.0) then
> c lambda/hybrid controller require delaydotlCE
>                 dotlCEdelay = delaydotlCE(ncycle,hsv,cm,idpart,tt)    
>             endif
>        else if (ncycle.GE.3.and.cm(33).ge.1.0.and.cm(37).EQ.0.0) then
> c all controllers without delay
>             dotlCEdelay=hsv(12)
>             lCEdelay=hsv(10)
>        end if
>        hsv(20) = lCEdelay
>        hsv(21) = dotlCEdelay
> c       
> c############################################
> c#   calc activation                        #
> c############################################
> c
>        if (cm(1).EQ.0) then
> c      Use Activation Values directly
>           if (cm(2).LE.0.0) then
>              q=abs(cm(2))
> c      constant activation level
>           else
>              call crvval(crv,abs(cm(2)),tt,q,dq)
> c      activation level from curve
>           end if
>           hsv(3)=q
>        else if (cm(1).EQ.1.0) then
> c      Calc Activation with Zajac
> c      1. get STIM
>           if (cm(2).LE.0.0) then
>              STIM=abs(cm(2))
> c      constant STIM level
>           else
>              call crvval(crv,cm(2),tt,STIM,dSTIM)
> c      STIM level from curve
>           end if
> c      2. get tau
>           if (cm(4).LT.0.0) then
>              call crvval(crv,abs(cm(4)),tt,tau,dtau)
>           else
>              tau=cm(4)
>           end if
> c      3. calc new q
>           q=hsv(3)+calc_dq(hsv(3),cm,STIM,tau)*dt1
>           hsv(3)=q
>        else if (cm(1).EQ.2.0) then
> c      Calc Activation with Hatze
> c      1. get STIM
>           if (cm(2).LE.0.0.and.cm(33).EQ.0.0) then
>              STIM=abs(cm(2))
> c      constant STIM level
>           else if (cm(2).GT.0.0.and.cm(33).EQ.0.0) then
> c      STIM level from curve
>              call crvval(crv,cm(2),tt,STIM,dSTIM)
>           else if (cm(33).EQ.1.0) then
> c------HERE COMES THE LAMBDA CONTROLLER PART
>           if (cm(34).LE.0.0) then
>             lambda=-cm(34)
>           else
>             call crvval(crv,cm(34),tt,lambda,dlambda)
>           end if
>              if (tt.LT.cm(38)) then
>              call crvval(crv,cm(2),tt,STIM,dSTIM)
>              else
>              STIM=STIM_lambda(cm,hsv,lCEdelay,dotlCEdelay,tt,lambda)
>              end if
> c------HERE COMES THE HYBRID CONTROLLER PART
>           else if (cm(2).LE.0.0.and.cm(33).EQ.2.0) then
>              STIM=abs(cm(2))
>              if (cm(34).LE.0.0) then
>                 lambda=-cm(34)
>              else
>                 call crvval(crv,cm(34),tt,lambda,dlambda)
>              end if
>              STIM=STIM_hybrid(cm,hsv,STIM,lambda,dlambda,
>      1 lCEdelay, dotlCEdelay,tt)
>           else if (cm(2).GT.0.0.and.cm(33).EQ.2.0) then
>              call crvval(crv,cm(2),tt,STIM,dSTIM)
>              if (cm(34).LE.0.0) then
>                 lambda=-cm(34)
>              else
>                  call crvval(crv,cm(34),tt,lambda,dlambda)
>              end if
>              STIM=STIM_hybrid(cm,hsv,STIM,lambda,dlambda,
>      1 lCEdelay, dotlCEdelay,tt)
> c------HERE COMES THE REFLEXIVE CONTROLLER PART     
>           else if (cm(33).EQ.3.0.and.ncycle.GE.3) then
>              STIM=STIM_reflex(ncycle,cm,hsv,tt)
>           end if
>           hsv(2)=STIM
> c         calc new gamma
>           gam_rel=hsv(15)+calc_dgam(cm(7),STIM,hsv(15))*dt1
>           hsv(15)=gam_rel
> c         calc new rho
>           rho_act=cm(4)*cm(5)*(cm(6)-1.0)/(cm(6)-(hsv(10)/cm(10)))
>      1    *(hsv(10)/cm(10))
> c         calc new q
>           q=(cm(3)+(rho_act*gam_rel)**3)/(1.0+(rho_act*gam_rel)**3)
>           hsv(3)=q
>        end if
> c
> c############################################
> c#   calc initial muscle force equilibrium  #
> c############################################
>        if (ncycle.LT.2) then
>          dot_l_MTC=0.0
> c        q=q0 for initialization of l_CE
>          q=cm(3)
> c        Force resulting from the force equilibrium between PEE or CE and SEE or
> c        SDE at starting point (dot_l_CE = 0 ~> F_SDE = 0):
>        if (ncycle.GE.1) then
>          l_CE=ZEROIN(elleng,q,cm)
>          call hill_model(l_CE,elleng,dot_l_MTC,q,cm,
>      1   F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,dt1,hsv,dot_l_CE,tt,i)
>        end if
>        end if
> c
>        if (etype.eq.'tbeam') then
> c
>        if (ncycle.GE.2) then
> c
> c------INITIALIZATION
>         l_CE=hsv(10)
> c
> c  calc dot_l_MTC
> c
> c calculating from strain,timestep and part length
>         dot_l_MTC=eps(i)/dt1*sqrt(sarea(i))
>         hsv(11) = dot_l_MTC
> c
> c############################################
> c#   update force/length/velocity           #
> c############################################
> c
>         call hill_model(l_CE,elleng,dot_l_MTC,q,cm,
>      1  F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,dt1,hsv,dot_l_CE,tt,i)
> c
>        end if
> c
> c  calc stresses from force
> c
>        if (etype.eq.'tbeam') then        
>         sig(1)=F_MTC/(voltot(i)/sqrt(sarea(i)))
> c  crosssection area is volume divided by part length
2515c2858,2860
< c     compute shear modulus, g
---
>         sig(2)=0.0
>         sig(3)=0.0
>        end if
2517,2518c2862
<       g2 =abs(cm(1))/(1.+cm(2))
<       g  =.5*g2
---
> c  write history variables
2520,2538c2864,2881
<       if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
<      1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
<      2     etype.eq.'sph  '.or.etype.eq.'sldax') then
<         if (cm(16).eq.1234567) then
<           call mitfail3d(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
<         else
<           if (.not.failel) then
<           davg=(-eps(1)-eps(2)-eps(3))/3.
<           p=-davg*abs(cm(1))/(1.-2.*cm(2))
<           sig(1)=sig(1)+p+g2*(eps(1)+davg)
<           sig(2)=sig(2)+p+g2*(eps(2)+davg)
<           sig(3)=sig(3)+p+g2*(eps(3)+davg)
<           sig(4)=sig(4)+g*eps(4)
<           sig(5)=sig(5)+g*eps(5)
<           sig(6)=sig(6)+g*eps(6)
<           if (cm(1).lt.0.) then            
<             if (sig(1).gt.cm(5)) failel=.true.
<           endif
<           endif
---
>        hsv(1) = sig(1)
>        hsv(4) = F_MTC
>        hsv(5) = F_CE
>        hsv(6) = F_PEE
>        hsv(7) = F_SEE
>        hsv(8) = F_SDE
>        hsv(9) = elleng
>        hsv(10)= l_CE
>        hsv(12)= dot_l_CE
>        hsv(14)= F_isom
> c  output is written only if output mode is set to 1...normal or 2...expert   
>         if((cm(29).eq.1).or.(cm(29).eq.2)) then   
> c  output is written only if time is GE than dt_out*counter_output    
>          if (tt.GE.(hsv(13)*cm(30))) then
>              call output(cm,idpart,hsv,tt)
> c  counter_output is adjusted
>              hsv(13)=hsv(13)+1
>          end if
2541,2596d2883
<       else if (etype.eq.'shell') then
<         if (cm(16).eq.1234567) then
<           call mitfailure(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
<         else
<           if (.not.failel) then
<           gc    =capa*g
<           q1    =abs(cm(1))*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
<           q3    =1./(q1+g2)
<           eps(3)=-q1*(eps(1)+eps(2))*q3
<           davg  =(-eps(1)-eps(2)-eps(3))/3.
<           p     =-davg*abs(cm(1))/(1.-2.*cm(2))
<           sig(1)=sig(1)+p+g2*(eps(1)+davg)
<           sig(2)=sig(2)+p+g2*(eps(2)+davg)
<           sig(3)=0.0
<           sig(4)=sig(4)+g *eps(4)
<           sig(5)=sig(5)+gc*eps(5)
<           sig(6)=sig(6)+gc*eps(6)
<           if (cm(1).lt.0.) then
<             if (sig(1).gt.cm(5)) failel=.true.
<           endif
<           endif
<         end if
<       elseif (etype.eq.'beam ' ) then
<           q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
<           q3    =q1+2.0*g
<           gc    =capa*g
<           deti  =1./(q3*q3-q1*q1)
<           c22i  = q3*deti
<           c23i  =-q1*deti
<           fac   =(c22i+c23i)*q1
<           eps(2)=-eps(1)*fac-sig(2)*c22i-sig(3)*c23i
<           eps(3)=-eps(1)*fac-sig(2)*c23i-sig(3)*c22i
<           davg  =(-eps(1)-eps(2)-eps(3))/3.
<           p     =-davg*cm(1)/(1.-2.*cm(2))
<           sig(1)=sig(1)+p+g2*(eps(1)+davg)
<           sig(2)=0.0
<           sig(3)=0.0
<           sig(4)=sig(4)+gc*eps(4)
<           sig(5)=0.0
<           sig(6)=sig(6)+gc*eps(6)
< c
<       elseif (etype.eq.'tbeam') then
<         q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
<         q3    =q1+2.0*g
<         deti  =1./(q3*q3-q1*q1)
<         c22i  = q3*deti
<         c23i  =-q1*deti
<         fac   =(c22i+c23i)*q1
<         eps(2)=-eps(1)*fac
<         eps(3)=-eps(1)*fac
<         davg  =(-eps(1)-eps(2)-eps(3))/3.
<         p     =-davg*cm(1)/(1.-2.*cm(2))
<         sig(1)=sig(1)+p+g2*(eps(1)+davg)
<         sig(2)=0.0
<         sig(3)=0.0
< c
2610a2898,3413
> c END OF SUBROUTINE UMAT41
> c
> c############################################
> c#   hill-model subroutine                  #
> c############################################
> c
>       subroutine hill_model(l_CE,elleng,dot_l_MTC,q,cm,
>      1 F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,dt1,hsv,dot_l_CE,tt,i)
> c
>       dimension cm(*),hsv(*)
>       real*8 l_PEE0, K_PEE, d_SE_max, l_SEE_nll, v_SEE, K_SEE_nl
>       real*8 l_CE, F_isom, F_PEE, F_SEE, l_SEE, L_A_rel, A_rel, q
>       real*8 B_rel, Q_Brel, D0, C2, C1, C0, dot_l_CE
>       real*8 O1, O2, O3, F_CE, F_CE_init, F_SDE, F_MTC
>       real*8 dot_l_MTC, elleng, K_SEE_l, Q_Arel, F_SUM 
>       real*8 l_MTC_0, dt1, tt, epsilon, dot_l_SDE
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
> c
> c  isometric force
> c
>          F_isom=calc_F_isom(l_CE,cm)
> c
> c  force of the parallel elastic element PEE
> c
>          F_PEE=calc_F_PEE(l_CE,cm)
> c
> c  force of the serial elastic element SEE
> c
>          F_SEE=calc_F_SEE(l_CE,cm,elleng)
> c
> c  Hill parameters concentric contraction
> c
>          if (l_CE.LT.cm(10)) then
>            L_A_rel=1.0
>          else 
>            L_A_rel=F_isom
>          end if
>          Q_Arel=(0.25)*(1.0+3.0*q)
>          A_rel=cm(15)*L_A_rel*Q_Arel
>          Q_Brel=(1.0/7.0)*(3.0+4.0*q)
>          B_rel=cm(16)*Q_Brel
> c
> c  calculate CE contraction velocity
> c
>          call calc_damping(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,D0,C2,C1,C0)
> c         
>          if ((C1**2.0-4.0*C2*C0).LT.0.0) then
>            dot_l_CE=0.0
>          else if ((D0.EQ.0.0).and.(C2.EQ.0.0)) then
> c        no damping
>            dot_l_CE=C0/C1
>          else
>            dot_l_CE=(-C1-sqrt(C1**2.0-4.0*C2*C0))/(2.0*C2)
>          end if
> c
> c  in case of an eccentric contraction
> c
>          if (dot_l_CE.GT.0.0) then
>             B_rel=(q*F_isom*(1.0-cm(18))/(q*F_isom+A_rel)
>      1 *B_rel/cm(17))
>             A_rel=-cm(18)*q*F_isom
> c
> c  calculate CE eccentric velocity
> c
>             call calc_damping(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,D0,C2,C1,C0)
> c            
>             if ((C1**2.0-4.0*C2*C0).LT.0.0) then
>               dot_l_CE=0.0
>             else if ((D0.EQ.0.0).and.(C2.EQ.0.0)) then
> c           no damping
>               dot_l_CE=-C0/C1
>             else
>               dot_l_CE=(-C1+sqrt(C1**2.0-4.0*C2*C0))/(2.0*C2)
>             end if
>             if (dot_l_CE.LT.0.0) then
> c            call usermsg('eccentric case error in muscle')
>             end if
>          end if
> c
> c  numerical restriction
> c
>          if (l_CE.LT.(0.01*cm(10))) then
>             dot_l_CE=0.0
>          end if
>          if (l_CE.GT.(1.99*cm(10))) then
>            l_SEE=elleng-l_CE
>            dot_l_CE=dot_l_MTC/(1.0+(K_PEE*cm(20)*(l_CE-l_PEE0)**
>      1     (cm(20)-1.0))/(K_SEE_nl*v_SEE*(l_SEE-cm(22))**(v_SEE-1.0)))
>          end if
> c
> c  contractile element force
> c
>          F_CE=cm(9)*(((q*F_isom+A_rel)/(1.0-dot_l_CE/(cm(10)*B_rel)))
>      1   -A_rel)
> c
> c  force of the serial damping element
> c
>          F_SDE=d_SE_max*((1.0-cm(28))*((F_CE+F_PEE)/cm(9))+cm(28))
>      1   *(dot_l_MTC-dot_l_CE)
>          F_MTC=F_SEE+F_SDE
> c
> c  calc l_CE (integrate dot_l_CE with finite difference methode)
> c
>          l_CE=l_CE+dot_l_CE*dt1
>          hsv(10)=l_CE
>          hsv(12)=dot_l_CE
> c
>          return
>          end
> c END OF SUBROUTINE HILL_MODEL
> c      
> c
> c############################################
> c#   ZEROIN function                        #
> c############################################
> c
>       real*8 function ZEROIN(elleng,act,cm)
> c
> c     Put a root-finding algorithm here. We used the ZEROIN function
> c     from Sec. 7.2 of:
> c
> c     Forsythe, G.E.; Malcolm, M.A.; Moler, C.B.: Computer Methods for
> c     Mathematical Computations. Prentice Hall Professional Technical 
> c     Reference, 1977
> c
> c     as interval we used [0,elleng], the function to evaluate F(X)
> c     is calc_F_sum(X,elleng,act,cm).
> 
> c
> c [...]
> c 
> c END OF FUNCTION ZEROIN
> c
> c############################################
> c#   calc_F_SUM function                    #
> c############################################
> c
>       real*8 function calc_F_SUM(l_CE,elleng,q,cm)
>          dimension cm(*)
>          real*8 F_isom, F_PEE, F_SEE, F_CE_init 
>          real*8 l_CE, elleng, q
> c
> c  calc coefficients
> c
>          F_isom=calc_F_isom(l_CE,cm)
>          F_PEE=calc_F_PEE(l_CE,cm)
>          F_SEE=calc_F_SEE(l_CE,cm,elleng) 
>          F_CE_init=cm(9)*q*F_isom
>          calc_F_SUM=F_SEE-F_CE_init-F_PEE
>       return
>       end
> c
> c END OF FUNCTION CALC_F_SUM
> c
> c############################################
> c#   calc_F_isom function                   #
> c############################################
> c
> c
>       real*8 function calc_F_isom(l_CE,cm)
> c  Isometric Force
> c
>          real*8 l_CE
>          dimension cm(*)
>          if (l_CE.GE.cm(10)) then
>            calc_F_isom=exp(-(abs(((l_CE/cm(10))-1.0)/cm(11)))**cm(12))
>          else
>            calc_F_isom=exp(-(abs(((l_CE/cm(10))-1.0)/cm(13)))**cm(14))
>          end if
>       return
>       end
> c
> c END OF FUNCTION CALC_F_ISOM
> c
> c############################################
> c#   calc_F_PEE function                    #
> c############################################
> c
>       real*8 function calc_F_PEE(l_CE,cm)
> c  Force of the parallel elastic element PEE
> c
>          dimension cm(*)
>          real*8 l_CE, l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,l_MTC_0,epsilon   
>          common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,l_MTC_0,epsilon
> c
>          if (l_CE.GE.l_PEE0) then
>            calc_F_PEE=K_PEE*((l_CE-l_PEE0)**cm(20))
>          else
>            calc_F_PEE=0.0
>          end if
>       return
>       end
> c
> c END OF FUNCTION CALC_F_PEE
> c
> c############################################
> c#   calc_F_SEE function                    #
> c############################################
> c
> c
>       real*8 function calc_F_SEE(l_CE,cm,elleng)
> c   5    |    5    |    5    |    5    |    5    |    5    |    5    |    5    |
> c  Force of the serial elastic element SEE
> c
>          dimension cm(*)
>          real*8 l_CE,K_SEE_nl,K_SEE_l,l_SEE,l_SEE_nll
>          real*8 v_SEE,elleng,d_SE_max,l_PEE0,K_PEE,epsilon
>          common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
> c
>          l_SEE=elleng-l_CE
>          if ((l_SEE.LT.l_SEE_nll).and.(l_SEE.GT.cm(22))) then
>            calc_F_SEE=K_SEE_nl*((l_SEE-cm(22))**v_SEE)
>          else if (l_SEE.GE.l_SEE_nll) then
>            calc_F_SEE=cm(25)+K_SEE_l*(l_SEE-l_SEE_nll)
>          else
>            calc_F_SEE=0.0
>            
>          end if
>       return
>       end
> c
> c END OF FUNCTION CALC_F_SEE
> c
> c############################################
> c#   calc_dgam function                     #
> c############################################
> c
>       real*8 function calc_dgam(m,STIM,gam_rel)
>          real*8 m,STIM,gam_rel
>          calc_dgam=m*(STIM-gam_rel)
>       return
>       end
> c
> c END OF FUNCTION CALC_DGAM
> c
> c############################################
> c#   calc_dq function                       #
> c############################################
> c
>       real*8 function calc_dq(q,cm,STIM,tau_q)
> c  calculates activation with the differential equation
> c
>          dimension cm(*)
>          real*8 q,STIM,tau_q
> c
>          calc_dq=1/tau_q*(STIM-STIM*(1.0-cm(5))*(q-cm(3))-cm(5)
>      1   *(q-cm(3)))
>       return
>       end
> c
> c END OF FUNCTION CALC_DQ
> c
> c############################################
> c#   calc_damping subroutine                #
> c############################################
> c
>       subroutine calc_damping(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,D0,C2,C1,C0)
>         dimension cm(*)
>         real*8 q,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1   dot_l_MTC,d_SE_max,D0,C2,C1,C0,epsilon
>         common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1   ,K_SEE_l,epsilon
> c    Force-dependent serial (SE)-Damping
> c
>          D0=cm(10)*B_rel*d_SE_max*(cm(28)+(1.0-cm(28))
>      1   *(q*F_isom+F_PEE/cm(9)))
>          C2=d_SE_max*(cm(28)-(A_rel-F_PEE/cm(9))*(1.0-cm(28)))
>          C1=-C2*dot_l_MTC-D0-F_SEE+F_PEE-cm(9)*A_rel
>          C0=D0*dot_l_MTC+cm(10)*B_rel*(F_SEE-F_PEE-cm(9)*q*F_isom)
> c
>       return
>       end
> c
> c END OF SUBROUTINE CALC_DAMPING
> c
> c############################################
> c#   Lambda Controller function             #
> c############################################
>       real*8 function STIM_lambda(cm,hsv,lCEdelay,dotlCEdelay,tt,lambda)
>          dimension cm(*),hsv(*)
>          real*8 lCEdelay,dotlCEdelay,lambda,tt     
> c         calc routine only if muscle is longer than the target            
>           if (lCEdelay.GT.lambda) then
>           STIM_lambda=(cm(35)*(lCEdelay-lambda)+
>      1       cm(36)*(dotlCEdelay-dlambda))/cm(10)
>           if (STIM_lambda.GT.1.0) then
>             STIM_lambda=1.0
>           else if (STIM_lambda.LT.0.0) then
>             STIM_lambda=0.0
>           else
>           continue
>           end if
>           else 
>             STIM_lambda=0.0
>           end if 
>       return
>       end
> c
> c END OF FUNCTION STIM_LABMDA
> c
> c############################################
> c#   Hybrid Controller function             #
> c############################################
>       real*8 function STIM_hybrid(cm,hsv,STIM_open,lambda,dlambda,
>      1 lCEdelay,dotlCEdelay,tt)
>          dimension cm(*),hsv(*)
>          real*8 STIM_open, lambda, dlambda, tt
>          real*8 lCEdelay, dotlCEdelay
>           if (lCEdelay.GT.lambda) then
>           STIM_hybrid=(STIM_open+(cm(35)*(lCEdelay-lambda)
>      1      +cm(36)*(dotlCEdelay-dlambda))/cm(10))
>           else 
>             STIM_hybrid=STIM_open
>           end if
>           if (STIM_hybrid.GT.1.0) then
>             STIM_hybrid=1.0
>           else if (STIM_hybrid.LT.0.0) then
>             STIM_hybrid=0.0
>           end if
>       return
>       end
> c
> c END OF FUNCTION STIM_HYBRID
> c
> c
> c############################################
> c#   Reflexive Controller function          #
> c############################################
>       real*8 function STIM_reflex(ncycle,cm,hsv,tt)
>          dimension cm(*),hsv(*)
>          real*8 tt
>          integer ncycle
> c        initialization of strain (hsv(23))
>          hsv(23) = 0
> c        initialization of l_CE_ref = l_CE(t=0) OR
> c        if time < t_PreSim then l_CE_ref = l_CE --> l_CE_ref = l_CE(t=t_PreSim)
>          if (ncycle.LE.3.OR.tt.LE.cm(38)) then
>              hsv(22)= hsv(10)         
>          endif
> c   hsv(22) (=l_CE_ref) = l_CE(t=t_PreSim) 
> c   if time > t_PreSim + delay: calc strain (strain = 0 for t<t_PreSim + delay)
>          if(tt.GT.(cm(38)+cm(37))) then
> c   hsv(23) (=strain) = (lCEdelay-l_CE_ref)/l_CE_ref
>           hsv(23) = (hsv(20)-hsv(22))/hsv(22)       
> c   if strain bigger than threshold --> activate reflex
>           if (hsv(23).GT.cm(39)) then
>            STIM_reflex = 1
> c   if reflex was active and strain bigger than zero --> activate reflex 
>           elseif (hsv(24).EQ.1.and.hsv(23).GE.0) then
>            STIM_reflex = 1
> c   otherwise --> deactivate reflex           
>           else
>            STIM_reflex = 0
>           end if    
>          end if
>          hsv(24) = STIM_reflex
>       return
>       end
> c
> c END OF FUNCTION STIM_reflex
> c
> c############################################
> c#   output subroutine                      #
> c############################################
> c
>       subroutine output(cm,idpart,hsv,tt)
> c
>         integer idpart
>         real*8 tt
>         dimension cm(*),hsv(*)
>         character*80 fname,fname1
> c       create filename musout.[PARTID]
>         write(fname1,'(I10.10)') idpart
>         fname = 'musout.'//fname1
>         if (cm(29).eq.1.or.
>      1     (cm(33).NE.1.and.cm(33).NE.2.and.cm(33).NE.3.)) then
> c         write simple output if output option == 1 or controller card is NOT set        
>           if (hsv(13).eq.0.0) then
> c          open/create the file under the filename created above
>            OPEN(86,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN')
> c          write header at the beginning
>            write (86,'('' output for muscle (PartID):'',I10)')idpart
>            write (86,'(''         time     stim_tot            q''
>      1       ''        f_mtc         f_ce        f_pee        f_see''
>      2       ''        f_sde        l_mtc         l_ce'')')
>            write (86,'(10ES13.5E2)')tt,hsv(2:10)
>            CLOSE(86)
>           else
> c          open/create the file under the filename created above
>            OPEN(86,FILE=fname,ACCESS='APPEND',FORM='FORMATTED',
>      1      STATUS='UNKNOWN')
> c          write normal data without header from now on
>            write (86,'(10ES13.5E2)')tt,hsv(2:10)
>            CLOSE(86)
>           end if
>         else 
> c         write extended output if output option == 2 and controller card is set
>           if (hsv(13).eq.0.0) then
> c          open/create the file under the filename created above
>            OPEN(86,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN')
> c          write header at the beginning
>            write (86,'(''advanced output for muscle (PartID):'',I10)')
>      1      idpart
>            write (86,'(''         time     stim_tot            q''
>      1       ''        f_mtc         f_ce        f_pee        f_see''
>      2       ''        f_sde        l_mtc         l_ce''
>      3       ''    dot_l_mtc     dot_l_ce     del_l_ce del_dot_l_ce'')')
>            write (86,'(14ES13.5E2)')tt,hsv(2:10),hsv(11),hsv(12),
>      1       hsv(20),hsv(21)
>            CLOSE(86)
>           else
> c          open/create the file under the filename created above
>            OPEN(86,FILE=fname,ACCESS='APPEND',FORM='FORMATTED',
>      1      STATUS='UNKNOWN')
> c          write normal data without header from now on
>            write (86,'(14ES13.5E2)')tt,hsv(2:10),hsv(11),hsv(12),
>      1       hsv(20),hsv(21)
>            CLOSE(86)
>           end if
>         end if  
>         return
>       end
> c      
> c END OF SUBROUTINE OUTPUT
> c
> c############################################
> c#   delaylCE function                      #
> c############################################
>       real*8 function delaylCE(ncycle,hsv,cm, idpart,tt)
>         real*8 tt, lasttime, dt
>         integer idpart, ncycle
>         dimension cm(*),hsv(*) 
>         integer index1,index2,indexr1, indexr2, iter
> c param declaration
>          lasttime = hsv(34);
> c dt = delay/buffersize
>          dt       = cm(37)/hsv(30);
>          indexr1  = hsv(32);
>          indexr2  = hsv(33);
>          if (ncycle.LE.3) then
> c  initialization
> c  set ring buffer values to 0.   
>           do iter=1,(hsv(30))
>            hsv(hsv(31)+iter-1) = 0;
>           enddo
>          elseif (tt.gt.(lasttime + dt)) then
>           lasttime = lasttime + dt
> c   mod ensures jump from end to the beginning of the ringbuffer
>           indexr1 = mod((hsv(32)+1),(hsv(30)))
>           indexr2 = mod((hsv(33)+1),(hsv(30)))
> c   save current lce
>           index1 = hsv(31) + indexr1
>           hsv(index1) = hsv(10);
>          end if
> c   define absolute index position of current delayedlCE value     
>          index2 = hsv(31) + indexr2
> c   return delayedlCE
>          delaylCE = hsv(index2)       
>          hsv(34) = lasttime
>          hsv(32) = indexr1
>          hsv(33) = indexr2
>          return
>       end
> c
> c END OF FUNCTION DELAYCE
> c
> 
> c############################################
> c#   delayldotCE function                   #
> c############################################
>       real*8 function delaydotlCE(ncycle,hsv,cm, idpart,tt)
>         real*8 tt,lasttimed,dt
>         integer idpart, ncycle
>         dimension cm(*),hsv(*) 
>         integer indexd1,indexd2,indexdr1,indexdr2,iterd
> c   ---------initialization ------
> c   param declaration
>          lasttimed = hsv(39)
>          dt        = cm(37)/hsv(30)
>          indexdr1  = hsv(37)
>          indexdr2  = hsv(38)
>          if (ncycle.LE.3) then
> c   set ring buffer values to 0.   
>           do iterd=1,(hsv(30))
>            hsv(hsv(36)+iterd-1) = 0
>           enddo
>          elseif (tt.gt.(lasttimed + dt)) then
>           lasttimed = lasttimed + dt
> c   mod ensures jump from end to the beginning of the ringbuffer
>           indexdr1 = mod((hsv(37)+1),(hsv(30)))
>           indexdr2 = mod((hsv(38)+1),(hsv(30)))
>           indexd1 = hsv(36) + indexdr1
> c   save current dot_lce
>           hsv(indexd1) = hsv(12)
>          end if
> c   define absolute index position of current delayeddotlCE value 
>          indexd2 = hsv(36) + indexdr2
> c   return delayeddotlCE
>          delaydotlCE = hsv(indexd2)
>          hsv(39) = lasttimed
>          hsv(37) = indexdr1
>          hsv(38) = indexdr2
>         return
>       end
> c
> c END OF FUNCTION DELAYDOTLCE
> c
> c END OF UMAT41
> c
