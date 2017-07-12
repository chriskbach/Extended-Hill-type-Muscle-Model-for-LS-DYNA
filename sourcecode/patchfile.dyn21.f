2254c2254,2255
<      1   temper,failel,crv,a(lcma),0,elsiz,idele)
---
>      1   temper,failel,crv,a(lcma),0,elsiz,idele,r_mem(dm_x),
>      2   r_mem(dm_v),r_mem(dm_a),i,d1,d2,d3)
2396a2398,2435
> C
> C % If you use this material model for scientific purposes, please cite
> C % the original research article:
> C % C. Kleinbach, O. Martynenko, J. Promies, D.F.B. Haeufle, J. Fehr,
> C % S. Schmitt: Implementation and Validation of the Extended Hill-type
> C % Muscle Model with Robust Routing Capabilities in LS-DYNA for Active
> C % Human Body Models, Biomedical Engineering Online, 2017.
> C %
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
2398c2437,2439
<      1 temper,failel,crv,cma,qmat,elsiz,idele)
---
>      1 temper,failel,crv,cma,qmat,elsiz,idele,x,v,a,i,d1,d2,d3)
> c
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
2407c2448,2497
< c     isotropic elastic material (sample user subroutine)
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
> c
> c     cm(3)=q0 minimum value of q
> c     cm(4)=tau_q time constant of rising activation LT0 curve id of tau_q over time
> c     cm(4)=c
> c     cm(5)=beta_q ratio between tau_q and time constant of falling activation
> c     cm(5)=eta
> c     cm(6)=k
> c     cm(7)=m
> c     cm(8)=additional constant muscle length
> c
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
2409c2499,2502
< c     Variables
---
> c     cm(22)=l_SEE_0
> c     cm(23)=dU_SEE_nllfindkParams
> c     cm(24)=dU_SEE_l
> c     cm(25)=dF_SEE_0
2411,2416c2504,2530
< c     cm(1)=first material constant, here young's modulus
< c     cm(2)=second material constant, here poisson's ratio
< c        .
< c        .
< c        .
< c     cm(n)=nth material constant
---
> c     Variables (Damping element)
> c
> c     cm(26)=damping methode
> c     cm(27)=d_PE
> c     cm(27)=d_SE
> c     cm(27)=D_SE
> c     cm(28)=R_SE
> c
> c
> c     Variable (Output definition; fort.(idele))
> c
> c     cm(29)=output method (EQ.0.  no output 
> c                           EQ.1.  basic output (idele, tt, ncycle, q) 
> c                           EQ.2.  basic+Forces (i...q ,F_MTC, F_SEE, F_SDE, F_CE, F_PEE, F_isom)
> c                           EQ.-1. basic+Lengths (i...q, elleng, l_CE, dot_l_MTC, dot_l_CE)
> c                           EQ.-2. basic+Forces+Lengths (i...q, F..., elleng...)
> c
> c     cm(30)=timestep of outputfile 
> c
> c     Variables for Controller
> c
> c     cm(33)=Activation Method
> c     cm(34)=target l_CE
> c     cm(35)=kp
> c     cm(36)=kd
> c     cm(37)=delay
> c     cm(38)=time till swap from alpha to lambda
2419c2533
< c     eps(2)=local y  strain increment
---
> c     eps(2)=local y  strain incrementhsv(14)
2422c2536
< c     eps(5)=local yz strain increment
---
> c     eps(5)=local yz strain incrementhsv(14)
2432,2438c2546,2564
< c     hsv(1)=1st history variable
< c     hsv(2)=2nd history variable
< c        .
< c        .
< c        .
< c        .
< c     hsv(n)=nth history variable
---
> c     d1  - strain rate/increment in x  direction, local x for shells
> c     d2  - strain rate/increment in y  direction, local y for shells
> c     d3  - strain rate/increment in z  direction, local z for shells
> c
> c         hsv(1)=sig(1)
> c         hsv(2)=F_MTC
> c         hsv(3)=F_SEE
> c         hsv(4)=F_SDE
> c         hsv(5)=F_CE
> c         hsv(6)=F_PEE
> c         hsv(7)=F_isom
> c         hsv(8)=elleng
> c         hsv(9)=l_CE
> c         hsv(10)=dot_l_CE
> c         hsv(11)=l_MTC_0
> c         hsv(12)=q_old
> c         hsv(13)=size required for ring buffer for lCE_delay
> c	  hsv(14)=gam_rel
> c         hsv(15)=counter_output
2446c2572
< c        eq."sldax" for shell forms 13 (2D solids - plane strain)
---
> c        eq."sldax" for shell forms 13 (2D solids -hsv(14) plane strain)
2464c2590
< c     cma=additional memory for material data defined by LMCA at 
---
> c     cma=additional memory for material data definhsv(14)ed by LMCA at 
2487c2613,2630
<       dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
---
> c
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
>       common/aux33loc/ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),
>      1 ix6(nlq),ix7(nlq),ix8(nlq),mxt(nlq),ix9(nlq),ix10(nlq)
>       common/prescloc/voltot(nlq)
>       common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
>      + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
>      + grvity,idirgv,nodspc,nspcor,nusa
> c$omp threadprivate (/aux33loc/)
> c$omp threadprivate (/prescloc/)
> c
> c############################################
> c#   Initialisierung der Vars               #
> c############################################
> c
>       dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
>       dimension x(3,*),v(3,*),a(3,*),qmat(3,3)
>       dimension d1(*),d2(*),d3(*)
2489a2633,2825
>       real*8 A1,B1,dot_h,dot_h_0,ent_rate,mech_eff
>       real*8 l_PEE0, K_PEE, d_SE_max, l_SEE_nll, v_SEE, K_SEE_nl
>       real*8 F_isom, F_PEE, F_SEE, l_SEE, L_A_rel, A_rel, q,dq
>       real*8 B_rel, L_B_rel, Q_Brel, D0, C2, C1, C0, dot_l_CE
>       real*8 v_max, P, O1, O2, O3, F_CE, F_CE_init, F_SDE, F_MTC
>       real*8 dot_l_MTC,Q_Arel, K_SEE_l, F_SUM, AREA,l_CE,q_old
>       real*8 n1,n2,n3,v1,v2,dv1,dv2,tol,l_MTC_0,l_CE_temp,STIM,dSTIM
>       real*8 tau,dtau,epsilon, wirk,E_PEE,E_SEE,E_CE,E_MTC,E_SDE,E_ges
>       real*8 delay, gam_rel, rho_act, lCEdelay, dotlCEdelay
>       real*8 lambda,dlambda
>       real*8 :: ZEROIN
>       real*8 :: calc_F_SUM
>       integer curve_id, var,k,j
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
> c############################################
> c#   Calc Muscle Parameters                 #
> c############################################
>          l_PEE0=cm(19)*cm(10)
>          if (cm(26).EQ.3.0) then
>          d_SE_max=cm(27)*(cm(9)*cm(15))/(cm(10)*cm(16))
>          else
>          d_SE_max=0.0
>          end if
>          l_SEE_nll=(1.0+cm(23))*cm(22)
>          v_SEE=cm(23)/cm(24)
>          K_SEE_nl=cm(25)/(cm(23)*cm(22))**v_SEE
>          K_SEE_l=cm(25)/(cm(24)*cm(22))
>          K_PEE=cm(21)*(cm(9)/(cm(10)*(cm(11)+1.0-cm(19)))**cm(20))
> c############################################
> c#      calc Length of Truss Element        #
> c############################################ 
>         elleng=sqrt((x(1,ix1(i))-x(1,ix2(i)))**2
>      1  +(x(2,ix1(i))-x(2,ix2(i)))**2
>      2  +(x(3,ix1(i))-x(3,ix2(i)))**2)+cm(8)
> c muscle length is element length plus offset
> c     
>        if (ncycle.LT.1) then
> c       
>          epsilon = 1.0
>    99    epsilon = 0.5*epsilon
>          if ( 1.0 + 0.5*epsilon .GT. 1.0 ) GOTO 99
> c
>        hsv(11)=elleng
> c     
>        hsv(12)=cm(3)
>        hsv(14)=0.0
>        hsv(15)=0
> c
>        end if
> c############################################
> c#   Calc delayed lCE/dotlCE                #
> c############################################       
> c
> c       if (ncycle.EQ.1.and.cm(33).ge.1.0.and.cm(37).GT.0.0) then
> c              hsv(13)=anint(cm(37)/dt1)
> c       end if
> cc
> c       if (ncycle.GE.1.and.cm(33).ge.1.0.and.cm(37).GT.0.0) then
> c       call delaylCE(i,ncycle,hsv(9),numelb,
> c     1 int(hsv(13)),cm(37),lCEdelay,tt)
> c       call delaydotlCE(i,ncycle,hsv(10),numelb,
> c     1 int(hsv(13)),cm(37),dotlCEdelay,tt)
> c       else if (ncycle.GE.1.and.cm(33).ge.1.0.and.cm(37).EQ.0.0) then
> c       dotlCEdelay=hsv(10)
> c       lCEdelay=hsv(9)
> c       end if
> c       
> c
> c############################################
> c#   Calc Activation                        #
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
>           hsv(12)=q
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
> c      variable tau? not documented!
>           else
>              tau=cm(4)
>           end if
> c      3. calc new q
>           q=hsv(12)+calc_dq(hsv(12),cm,STIM,tau)*dt1
>           hsv(12)=q
>        else if (cm(1).EQ.2.0) then
> c      Calc Activation with Hatze
> c      1. get STIM
>           if (cm(2).LE.0.0.and.cm(33).EQ.0.0) then
>              STIM=abs(cm(2))
> c      constant STIM level
>           else if (cm(2).GT.0.0.and.cm(33).EQ.0.0) then
> c      STIM level from curve
>              call crvval(crv,cm(2),tt,STIM,dSTIM)
> c          else if (cm(33).EQ.1.0) then
> cc      HERE COMES THE LAMBDA CONTROLLER PART
> c          if (cm(34).LE.0.0) then
> c            lambda=-cm(34)
> c          else
> c            call crvval(crv,cm(34),tt,lambda,dlambda)
> c          end if
> c             if (tt.LT.cm(38)) then
> c             call crvval(crv,cm(2),tt,STIM,dSTIM)
> c             else
> c             STIM=STIM_lambda(cm,hsv,lCEdelay,dotlCEdelay,tt,lambda)
> c             end if
> cc      HERE COMES THE HYBRID CONTROLLER PART
> c          else if (cm(2).LE.0.0.and.cm(33).EQ.2.0) then
> c             STIM=abs(cm(2))
> c             if (cm(34).LE.0.0) then
> c                lambda=-cm(34)
> c             else
> c                call crvval(crv,cm(34),tt,lambda,dlambda)
> c             end if
> c             STIM=STIM_hybrid(cm,hsv,STIM,lambda,dlambda,
> c     1 lCEdelay, dotlCEdelay,tt)
> c          else if (cm(2).GT.0.0.and.cm(33).EQ.2.0) then
> c             call crvval(crv,cm(2),tt,STIM,dSTIM)
> c             if (cm(34).LE.0.0) then
> c                lambda=-cm(34)
> c             else
> c                 call crvval(crv,cm(34),tt,lambda,dlambda)
> c             end if
> c             STIM=STIM_hybrid(cm,hsv,STIM,lambda,dlambda,
> c     1 lCEdelay, dotlCEdelay,tt)
>           end if
> c         calc new gamma
>           gam_rel=hsv(14)+calc_dgam(cm(7),STIM,hsv(14))*dt1
>           hsv(14)=gam_rel
> c         calc new rho
>           rho_act=cm(4)*cm(5)*(cm(6)-1.0)/(cm(6)-(hsv(9)/cm(10)))
>      1    *(hsv(9)/cm(10))
> c         calc new q
>           q=(cm(3)+(rho_act*gam_rel)**3)/(1.0+(rho_act*gam_rel)**3)
>           hsv(12)=q
>        end if
> c
> c      
> c
> c############################################
> c#   Calc Initial Muscle Force Equilibrium  #
> c############################################
>        if (ncycle.LT.1) then
>          dot_l_MTC=0.0
> c         
> c        Force resulting from the force equilibrium between PEE or CE and SEE or
> c        SDE at starting point (dot_l_CE = 0 ~> F_SDE = 0):
>          l_CE=ZEROIN(elleng,q,cm)
>          call hill_model(l_CE,elleng,dot_l_MTC,q,cm,
>      1   F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,dt1,hsv,dot_l_CE,tt,i)
>        end if
> c
>        if (etype.eq.'tbeam'.or.'hbeam') then
> c
>        if (ncycle.GE.1) then
> c
> c  INITIALIZATION
> c
>         l_CE=hsv(9)
> c
> c  calc dot_l_MTC
> c
>         n1=x(1,ix1(i))-x(1,ix2(i))
>         n2=x(2,ix1(i))-x(2,ix2(i))
>         n3=x(3,ix1(i))-x(3,ix2(i))
>         v1=v(1,ix1(i))*n1+v(2,ix1(i))*n2+v(3,ix1(i))*n3
>         v2=v(1,ix2(i))*n1+v(2,ix2(i))*n2+v(3,ix2(i))*n3
> c		scaling with element(truss) length
>         dv1=v1/(elleng-cm(8))
>         dv2=v2/(elleng-cm(8))
>         dot_l_MTC=dv1-dv2
> c
> c############################################
> c#   Update Force/Length/Velo               #
> c############################################
2491,2495c2827,2828
<       if (ncycle.eq.1) then
<         if (cm(16).ne.1234567) then
<           call usermsg('mat41')
<         endif
<       endif
---
>         call hill_model(l_CE,elleng,dot_l_MTC,q,cm,
>      1  F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,dt1,hsv,dot_l_CE,tt,i)
2497c2830
< c     compute shear modulus, g
---
>         end if
2499,2500d2831
<       g2 =abs(cm(1))/(1.+cm(2))
<       g  =.5*g2
2502,2521d2832
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
<         end if
2523,2575c2834,2838
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
---
> c  Calc Stresses from Force
> c
>         if (etype.eq.'tbeam') then        
>         sig(1)=F_MTC/(voltot(i)/(elleng-cm(8)))
> c		crosssection area is volume divided by element(truss) length
2577a2841,2861
>         end if
> c
> c  Write History Variables
> c
>          hsv(1)=sig(1)
>          hsv(2)=F_MTC
>          hsv(3)=F_SEE
>          hsv(4)=F_SDE
>          hsv(5)=F_CE
>          hsv(6)=F_PEE
>          hsv(7)=F_isom
>          hsv(8)=elleng
>          hsv(9)=l_CE
>          hsv(10)=dot_l_CE
>          if (tt.GE.(hsv(15)*cm(30))) then
> c  output is written only if time is GE than dt_out*counter_output
>          call output(cm,idele,F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,
>      1 elleng,l_CE,dot_l_MTC,dot_l_CE,ncycle,tt,dt1,q)
> c  counter_output is adjusted
>          hsv(15)=hsv(15)+1
>          end if
2592a2877,3336
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
>       subroutine hill_model(l_CE,elleng,dot_l_MTC,q,cm,
>      1 F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom,dt1,hsv,dot_l_CE,tt,i)
> c
>       dimension cm(*),hsv(*)
>       dimension qmat(3,3)
>       integer i
>       real*8 A1,B1,dot_h,dot_h_0,ent_rate,mech_eff
>       real*8 l_PEE0, K_PEE, d_SE_max, l_SEE_nll, v_SEE, K_SEE_nl
>       real*8 l_CE, F_isom, F_PEE, F_SEE, l_SEE, L_A_rel, A_rel, q
>       real*8 B_rel,Q_Brel,D0,C2,C1,C0,dot_l_CE
>       real*8 v_max, P, O1, O2, O3, F_CE, F_CE_init, F_SDE, F_MTC
>       real*8 dot_l_MTC, elleng, K_SEE_l, Q_Arel, F_SUM, AREA
>       real*8 l_MTC_0,dt1,tt,epsilon, Power, dot_l_SDE
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
> c
> c
>       include 'iounits.inc'
> c
> c  Isometric Force
> c
> c
>          F_isom=calc_F_isom(l_CE,cm)
> c
> c  Force of the parallel elastic element PEE
> c
>          F_PEE=calc_F_PEE(l_CE,cm)
> c
> c  Force of the serial elastic element SEE
> c
>          F_SEE=calc_F_SEE(l_CE,cm,elleng)
> c
> c  Hill Parameters concentric contraction
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
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
>          call calc_damping(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,D0,C2,C1,C0)
> c         
>          if ((C1**2.0-4.0*C2*C0).LT.0.0) then
>            dot_l_CE=0.0
>          else if ((D0.EQ.0.0).and.(C2.EQ.0.0)) then
>            dot_l_CE=C0/C1
>          else
>            dot_l_CE=(-C1-sqrt(C1**2.0-4.0*C2*C0))/(2.0*C2)
>          end if
> c
> c  in case of a negative CE-force (numerical restriction)
> c
>          v_max=B_rel/A_rel*cm(10)*q*F_isom
> c
> c   Future Task: analyse why dot_l_MTC is LT -vmax in quick_release applications
> c
> c         if (dot_l_MTC.LT.(-v_max)) then
> c           P=1e3*cm(9)*cm(10)*cm(16)/cm(15)
> c           O1=-1.0/(P*cm(9)*B_rel*cm(10))
> c           O2=1.0/(P*cm(9))-A_rel/(B_rel*cm(10))
> c           O3=-q*F_isom
> c           if ((O2**2.0-4.0*O1*O3).LT.0.0) then
> c              dot_l_CE=0.0
> c           else
> c              dot_l_CE=(-O2-sqrt(O2**2.0-4.0*O1*O3))/(2.0*O1)E
> c           end if
> c         end if
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
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
>             call calc_damping(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,D0,C2,C1,C0)
> c            
>             if ((C1**2.0-4.0*C2*C0).LT.0.0) then
>               dot_l_CE=0.0
>             else if ((D0.EQ.0.0).and.(C2.EQ.0.0)) then
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
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
>          if (l_CE.LT.(0.01*cm(10))) then
>             dot_l_CE=0.0
>          end if
>          if (l_CE.GT.(1.99*cm(10))) then
>            l_SEE=abs(elleng-l_CE)
>            dot_l_CE=dot_l_MTC/(1.0+(K_PEE*cm(20)*(l_CE-l_PEE0)**
>      1     (cm(20)-1.0))/(K_SEE_nl*v_SEE*(l_SEE-cm(22))**(v_SEE-1.0)))
>          end if
> c
> c  contractile element force
> c
>          F_CE=cm(9)*(((q*F_isom+A_rel)/(1.0-dot_l_CE/(cm(10)*B_rel)))
>      1   -A_rel)
> c
> c  Force of the serial damping element
> c
>          F_SDE=d_SE_max*((1.0-cm(28))*((F_CE+F_PEE)/cm(9))+cm(28))
>      1   *(dot_l_MTC-dot_l_CE)
>          F_MTC=F_SEE+F_SDE
> c
> c  Calc l_CE (integrate dot_l_CE with finite difference methode)
> c
>          l_CE=l_CE+dot_l_CE*dt1
>          hsv(9)=l_CE
>          hsv(10)=dot_l_CE
> c
> c
>          return
>          end
> c
> c############################################
> c#   ZEROIN subroutine                      #
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
> c
> [...]
> c END OF SUBROUTINE ZEROIN
> c
> c
>       real*8 function calc_F_SUM(l_CE,elleng,q,cm)
>       dimension cm(*)
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    | 2  5    |
>          real*8 l_PEE0, K_PEE, v_SEE, K_SEE_nl, K_SEE_l,l_SEE_nll,l_CE
>          real*8 F_isom, F_PEE, F_SEE, F_CE_init, elleng, q,l_SEE_nl
>          real*8 d_SE_max,l_MTC_0,epsilon
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
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
> c
>       real*8 function calc_F_isom(l_CE,cm)
> c   567 |    5    |    5    |    5    |    5    |    5    |    5    |    5    |
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
> c
>       real*8 function calc_F_PEE(l_CE,cm)
> c   5    |    5    |    5    |    5    |    5    |    5    |    5    |    5    |
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
> c
>       real*8 function calc_F_SEE(l_CE,cm,elleng)
> c   5    |    5    |    5    |    5    |    5    |    5    |    5    |    5    |
> c  Force of the serial elastic element SEE
> c
>          dimension cm(*)
>          real*8 l_CE,K_SEE_nl,K_SEE_l,l_SEE,l_SEE_nll
>          real*8 v_SEE,elleng,d_SE_max,l_PEE0,K_PEE,epsilon
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
> c
>          l_SEE=abs(elleng-l_CE)
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
>       real*8 function calc_dgam(m,STIM,gam_rel)
>          real*8 m,STIM,gam_rel
>          calc_dgam=m*(STIM-gam_rel)
>       return
>       end
> c
>       real*8 function calc_dq(q,cm,STIM,tau_q)
> c   5    |    5    |    5    |    5    |    5    |    5    |    5    |    5    |
> c  Calculates activation with the Differentialequation
> c
>          dimension cm(*)
>          real*8 q,STIM,tau_q
> c
>          calc_dq=1/tau_q*(STIM-STIM*(1.0-cm(5))*(q-cm(3))-cm(5)
>      1   *(q-cm(3)))
>       return
>       end
> c
>       subroutine calc_damping(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,D0,C2,C1,C0)
>       dimension cm(*)
>       real*8 q,A_rel,B_rel,F_SEE,F_PEE,F_isom,
>      1 dot_l_MTC,d_SE_max,D0,C2,C1,C0,epsilon
>       common /coeff/ l_PEE0,K_PEE,d_SE_max,l_SEE_nll,v_SEE,K_SEE_nl
>      1 ,K_SEE_l,epsilon
>       if (cm(26).EQ.1.0) then
> c      
> c    constant parallel (PE)-Damping
> c
>          D0=0.0
>          C2=cm(27)
>          C1=-(cm(27)*cm(10)*B_rel+F_SEE-F_PEE+cm(9)*A_rel)
>          C0=cm(10)*B_rel*(F_SEE-F_PEE-cm(9)*q*F_isom)
> c
>       else if (cm(26).EQ.2.0) then
> c      
> c    constant serial (SE)-Damping
> c
>          D0=0.0
>          C2=cm(27)
>          C1=-(cm(27)*(dot_l_MTC+cm(10)*B_rel)+F_SEE-F_PEE+cm(9)*A_rel)
>          C0=cm(10)*B_rel*(cm(27)*dot_l_MTC+F_SEE-F_PEE-cm(9)*q*F_isom)
> c
>       else if (cm(26).EQ.3.0) then
> c      
> c    Force-dependent serial (SE)-Damping
> c
>          D0=cm(10)*B_rel*d_SE_max*(cm(28)+(1.0-cm(28))
>      1   *(q*F_isom+F_PEE/cm(9)))
>          C2=d_SE_max*(cm(28)-(A_rel-F_PEE/cm(9))*(1.0-cm(28)))
>          C1=-C2*dot_l_MTC-D0-F_SEE+F_PEE-cm(9)*A_rel
>          C0=D0*dot_l_MTC+cm(10)*B_rel*(F_SEE-F_PEE-cm(9)*q*F_isom)
> c
>       else
> c      
> c    no Damping
> c      
>          D0=0.0
>          C2=0.0
>          C1=-(F_SEE-F_PEE+cm(9)*A_rel)
>          C0=cm(10)*B_rel*(F_PEE-F_SEE+cm(9)*q*F_isom)
>       end if
>       return
>       end
> c
> c      real*8 function STIM_lambda(cm,hsv,lCEdelay,dotlCEdelay,tt,lambda)
> cc############################################
> cc#            Lambda Controller             #
> cc############################################
> c         dimension cm(*),hsv(*)
> c         real*8 lCEdelay,dotlCEdelay,lambda
> c          STIM_lambda=(cm(35)*(lambda-
> c     1 lCEdelay)+cm(36)*(-dotlCEdelay))/cm(10)
> c          if (STIM_lambda.GT.1.0) then
> c          STIM_lambda=1.0
> c          else if (STIM_lambda.LT.0.0) then
> c          STIM_lambda=0.0
> c          else
> c          continue
> c          end if
> c      return
> c      end
> c
> c      real*8 function STIM_hybrid(cm,hsv,STIM_open,lambda,dlambda,
> c     1 lCEdelay,dotlCEdelay,tt)
> cc############################################
> cc#            Hybrid Controller             #
> cc############################################
> c         dimension cm(*),hsv(*)
> c         real*8 STIM_open, lambda,dlambda
> c         real*8 lCEdelay, dotlCEdelay
> c          STIM_hybrid=(STIM_open+(cm(35)*(lambda-lCEdelay)
> c     1 +cm(36)*(dlambda-dotlCEdelay))/cm(10))
> c          if (STIM_hybrid.GT.1.0) then
> c          STIM_hybrid=1.0
> c          else if (STIM_hybrid.LT.0.0) then
> c          STIM_hybrid=0.0
> c          else
> cc          write(115,*)STIM_hybrid
> c          end if
> c      return
> c      end
> c
>       subroutine output(cm,idele,F_MTC,F_SEE,F_SDE,F_CE,F_PEE,F_isom
>      1 ,elleng,l_CE,dot_l_MTC,dot_l_CE,nintcy,tt,dt1,q)
>       integer idele,nintcy
>       real*8 F_MTC,F_SDE,F_CE,F_PEE,elleng,l_CE,tt
>      1 ,dot_l_MTC,dot_l_CE,q
>       dimension cm(*) 
>       if (abs(cm(29)).GT.0) then
>        write(idele,*)idele
>        write(idele,*)tt
>        write(idele,*)nintcy
>        write(idele,*)q
>       if (abs(cm(29)).GE.2) then
>        write(idele,*)F_MTC
>        write(idele,*)F_SEE
>        write(idele,*)F_SDE
>        write(idele,*)F_CE
>        write(idele,*)F_PEE
>        write(idele,*)F_isom
>       end if
>       if (cm(29).LT.0) then
>        write(idele,*)elleng
>        write(idele,*)l_CE
>        write(idele,*)dot_l_MTC
>        write(idele,*)dot_l_CE
>       end if
>       end if
>       return
>       end
> c      
> c
> c      subroutine delaylCE(i,ncycle,safe,numelb,anzahl,delay
> c     1 ,delayval,tt)
> c      real*8, allocatable, dimension(:,:) :: histlCE
> c      real*8 tt,safe
> c      integer run,jump, anzahl,numelb,curlCE,i,dessteplCE
> c      save histlCE,curlCE,dessteplCE
> c      if (.not. allocated(histlCE)) allocate(histlCE(numelb,anzahl))
> c      if(ncycle.eq.1.and.i.eq.1) then
> c      curlCE=0
> c      end if
> cc
> c      if(ncycle.eq.1.and.i.eq.1) then
> c        run=1
> c        jump=1
> c        DO WHILE(run.LE.numelb)
> c            DO WHILE(jump.LE.anzahl)
> c                histlCE(run,jump)=0.0
> c                jump=jump+1
> c            END DO
> c            run=run+1
> c            jump=1
> c        END DO
> c      end if
> cc
> c      if(i.eq.1) then
> c        curlCE=curlCE+1
> c        if(curlCE.gt.anzahl) then
> c          curlCE=1
> c        end if
> c      end if
> cc
> c      histlCE(i,curlCE)=safe
> cc
> c      if (tt.lt.delay) then
> c         delayval=histlCE(i,1)
> c         dessteplCE=1
> c      else
> c       dessteplCE=curlCE+1
> c       if (dessteplCE.gt.anzahl) then
> c       dessteplCE=1.0
> c       end if
> c       delayval=histlCE(i,dessteplCE)
> c      end if
> cc     
> c      return
> c      end
> c
> c      subroutine delaydotlCE(i,ncycle,safe,numelb,anzahl,delay
> c     1 ,delayval,tt)
> c      real*8, allocatable, dimension(:,:) :: histdotlCE
> c      real*8 tt,safe
> c      integer run,jump, anzahl,numelb,curdotlCE,i,desstepdotlCE
> c      save histdotlCE,curdotlCE,desstepdotlCE
> c      if (.not. allocated(histdotlCE)) 
> c     1 allocate(histdotlCE(numelb,anzahl))
> c      if(ncycle.eq.1.and.i.eq.1) then
> c      curdotlCE=0
> c      end if
> cc
> c      if(ncycle.eq.1.and.i.eq.1) then
> c        run=1
> c        jump=1
> c        DO WHILE(run.LE.numelb)
> c            DO WHILE(jump.LE.anzahl)
> c                histdotlCE(run,jump)=0.0
> c                jump=jump+1
> c            END DO
> c            run=run+1
> c            jump=1
> c        END DO
> c      end if
> cc
> c      if(i.eq.1) then
> c        curdotlCE=curdotlCE+1
> c        if(curdotlCE.gt.anzahl) then
> c          curdotlCE=1
> c        end if
> c      end if
> cc
> c      histdotlCE(i,curdotlCE)=safe
> cc
> c      if (tt.lt.delay) then
> c         delayval=histdotlCE(i,1)
> c         desstepdotlCE=1
> c      else
> c       desstepdotlCE=curdotlCE+1
> c       if (desstepdotlCE.gt.anzahl) then
> c       desstepdotlCE=1.0
> c       end if
> c       delayval=histdotlCE(i,desstepdotlCE)
> c      end if
> cc     
> c      return
> c      end
> c
> c
