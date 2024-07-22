










      subroutine equation_of_state(txt_,time_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 365 - last update: 26-01-23
!______________________________________________________________________
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal
      use module_parallele
      implicit none
      character*17 txt_
      integer time_ ! after, now or before


!...............................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         04/11/03: passage d'un argument dans CALL DENSITY
!                   + nouveau cas destiné à calculer une densité pour l'équation
!                   de la tke qui permette d'équilibrer exactement les transferts
!                   d'energie du courant moyen (mélange vertical T & S) vers le
!                   courant turbulent (terme de flottabilité équation de tke)
!         11/12/04: Premier pas vers une equation d'etat non lineaire dependante
!                   de la pression.
!         28/01/05: ajout blindage. Une puissance 1.5 ne s'applique pas à un
!                   nbre negatif
!         10/05/06: fonctions compatibles avec double precision
!         11/12/07: correction facteur d'echelle dans cas ichoix=1
! 2009.3  01-10-09: Utilisation des nouveaux facteurs d'echelle verticale
!         11-10-09: tem_c et sal_c remplacent thz_z et shz_z
! 2010.11 21-07-10  debug
!         22-07-10  routine density renommee equation_of_state
! 2010.14 16-12-10  ajout EOS Wright (1997) et McDougall et al. (2003)
! S.26    01-04-14  possibilite d'appel a deep_convection
!         16-05-14  ajout routine equation_of_state_potloc_jmfwg
!                   ajout routine equation_of_state_potzref_jmfwg
!         17-05-14  routine equation_of_state_potzref_jmfwg appelee depuis
!                   read_ogcm_fields avec aiguillage selon variable entrante
!                   tem_t sal_t ou temobc_t salobc_t
!         22-05-14  suite point precedent prevoyant le cas ou temofl n'est pas allocate
!         19-06-14  Possibilite d'appeler la routine de densite potentielle a niveau
!                   de reference "z local" depuis module_offline ou read_ogcm.
!                   Ajout d'une routine "potentielle z local" rapide uniquement
!                   depuis le gradient de pression et uniquement pour la tke (pas
!                   compatible avec module_offline et read_ogcm)
!         20-06-14  ajout subroutine equation_of_state_potzref2
!         18-10-14  Ajout densite potentielle algo deduit de McDougall JPO 1987
!         20-10-14  possibilite de calcul densite potentielle zref=0 meme si eos_tkezref<0
!         14-09-15  - ajout de routines generiques equation_of_state_prs_anyv 
!                   et equation_of_state_pot_anyv
!                   - remplacement de const1 et cie par cst1_ et cie
!                   - equation_of_state_full_a1_a7
!         30-09-15 marqueurs de disponibilite pour tableaux generiques anyv3d
!         16-02-17 verification du contenu de anyv3d
! v253    05-05-19 suppression de call deep_convection
! v273    07-02-20 Appliquer A l'EOS des bornes de validites sur T et S
! v328    23-02-22 .or.txt_=='grp') then !23-02-22
! v358    29-10-22 timeobc_ pour temobc et salobc et bornage sur time_ (<=2)
! v365    26-01-23 - prs_=max(x0-0.5*dsig_t(i,j,k)*h_w(i,j),0.)
!                  - boucles de calculs etendues 0:imax+1,0:jmax+1 pour servir au tableau rho_t dans drifter
!...............................................................................

      if(time_/=1.and.time_/=2)   &
      stop ' STOP Bad time_ argument in subroutine equation_of_state'

! Compute the density anomaly (rh-rho) in kg/m3 at the pressure of reference
! rho is a reference density adjusted to the initial state in routine initial_state_eq.F90

      if(txt_=='potential density') then !-potential-potential-potential->

       if(eos_author==0)call equation_of_state_linear(time_)

       if(eos_author==3) then
! z reference different de z=0m:
        if(eos_pgfzref> 0.)call equation_of_state_potzref_jmfwg('pgf',time_)
! z reference z=0m:
        if(eos_pgfzref==0.)call equation_of_state_potential_jmfwg(time_)
! z reference "local depth": !19-06-14
!       if(eos_pgfzref< 0.)call equation_of_state_potloc_jmfwg('pgf',time_)
       endif

!     if(eos_author==1)call equation_of_state_potential_mellor91(time_)
!     if(eos_author==2)call equation_of_state_potential_w97(time_)

!      if(initial_main_status==1) then !111111111111>
!        if(convect_yn==3)call deep_convection !01-04-14 !commentE le 05-05-19
!      endif                           !111111111111>

      return
      endif                                 !-potential-potential-potential->


! Compute the part of the density that is dependent on pressure
! rho is a reference density adjusted to the initial state in routine initial_state_eq.F90
      if(txt_=='compression terms') then !-add-pressure-term-add-pressure-term->

      if(eos_comprs/=1)then ! compression in not computed if eos_comprs/=1
       call equation_of_state_nopressure
       return
      endif

      if(eos_author==1)call equation_of_state_pressure_mellor91(time_)
      if(eos_author==2)call equation_of_state_pressure_w97(time_)
      if(eos_author==3)call equation_of_state_pressure_jmfwg(time_)

      return
      endif                                 !-add-pressure-term-add-pressure-term->


      if(txt_=='at a single point') then !-single-single-single->

! compute the potential density around d single omain-averaged (T,S) in order
! to compute the coef of a linear equation of state:
       call equation_of_state_0d

       return
      endif                                 !-single-single-single->

      if(txt_=='at initial state ') then !-initial-initial->

      if(eos_author==0) then !>>>>>>
         if(eos_comprs==1)stop ' Err eos_author==0 & eos_comprs==1'
         call equation_of_state_potential_jmfwg(time_)
      endif                  !>>>>>>

      if(eos_author==3) then

         if(eos_comprs==1) then !19-06-14
! Si PGF Marsaleix 2011 la constante rho se calcule sur la base de la densite potentielle a z=0m
          call equation_of_state_potential_jmfwg(time_)
         endif
         if(eos_comprs==0) then !19-06-14
! Si PGF Marsaleix 2009 sans effet de compression la constante rho se calcule sur la base
! de la densite potentielle. Le niveau de reference depend du choix dans notebook_eqs
           if(eos_pgfzref==0.)call equation_of_state_potential_jmfwg(time_)
           if(eos_pgfzref> 0.)call equation_of_state_potzref_jmfwg('pgf',time_)
           if(eos_pgfzref< 0.)call equation_of_state_potloc_jmfwg('pgf',time_)
         endif

      endif

! Le PGF est calcule de telle sorte que la constante rho se base uniquement
! sur la densite potentielle. Voila pourquoi les lignes suivantes sont commentees.

!     if(eos_author==1)call equation_of_state_potential_mellor91(time_)
!     if(eos_author==2)call equation_of_state_potential_w97(time_)
!     if(eos_author==3)call equation_of_state_potential_jmfwg(time_)

!     if(eos_comprs/=1)then ! compression in not computed if eos_comprs/=1
!      call equation_of_state_nopressure
!      return
!     endif

!     if(eos_author==1)call equation_of_state_pressure_mellor91(time_)
!     if(eos_author==2)call equation_of_state_pressure_w97(time_)
!     if(eos_author==3)call equation_of_state_pressure_jmfwg(time_)

      return
      endif                                 !-initial-initial->

      stop ' STOP: bad txt_ argument in subroutine equation_of_state'

      end subroutine equation_of_state

!....................................................................................

      subroutine equation_of_state_linear(time_)
      use module_principal
      implicit none
      integer time_ ! after, now or before
      double precision cst1_,cst2_,cst3_

! A linear equation of state with alp_t, alp_s, t0, S0, rho adjusted to the initial state
! in routine initial_state_eq.F90

      cst1_=-rho*alp_t
      cst2_= rho*alp_s
      cst3_=-cst1_*t0-cst2_*s0
      do k=1,kmax
      do j=0,jmax+1 ; do i=0,imax+1 !26-01-23
      rhp_t(i,j,k)=cst1_*tem_t(i,j,k,time_)       &
                  +cst2_*sal_t(i,j,k,time_)       &
                  +cst3_
      enddo          ; enddo
      enddo

      end subroutine equation_of_state_linear

!....................................................................................

      subroutine equation_of_state_potential_mellor91(time_)
      use module_principal
      implicit none
      integer time_ ! after, now or before
      stop ' PASSER En POTLOC'

! Compute (density - rho in kg/m3) at the pressure of reference using the method described in
! Mellor (1991) Journal of Atmospheric and Oceanic Technology. pp 609-611

! rho is a reference density adjusted to the initial state in routine initial_state_eq.F90

      const5=-rho+999.842594 ! 1er coef de XA - rho (reference)

      do j=1,jmax
      do i=1,imax
      if(  mask_t(i,j,kmax+1).eq.1)then !>>>>>>>>>>>>>>>>>>>>
      do k=kmin_w(i,j),kmax

      x2=tem_t(i,j,k,time_)**2
      x3=tem_t(i,j,k,time_)**3
      x4=tem_t(i,j,k,time_)**4

! Calcul de la densite potentielle:

      rhp_t(i,j,k)=                                   &
! XA-RH0:
       const5+                                        &
       6.793952e-2*tem_t(i,j,k,time_)-             &
       9.095290e-3*x2+                                &
       1.001685e-4*x3-                                &
       1.120083e-6*x4+                                &
       6.536332e-9*tem_t(i,j,k,time_)**5           &

! +XB*S:
       +(                                             &
         8.24493e-1-                                  &
       4.0899e-3*tem_t(i,j,k,time_)+               &
       7.6438e-5*x2-                                  &
       8.2467e-7*x3+                                  &
       5.3875e-9*x4                                   &
        )*sal_t(i,j,k,time_)                       &

! +XC*S**1.5:
       +(                                             &
         -5.72466e-3+                                 &
       1.0227e-4*tem_t(i,j,k,time_)-               &
       1.6546e-6*x2                                   &
        )*max(sal_t(i,j,k,time_),zero)**1.5        &    !28/01/05

! +XD*S**2:
       +4.8314e-4*sal_t(i,j,k,time_)**2

      enddo
      endif                         !>>>>>>>>>>>>>>>>>>>>

      enddo
      enddo

      end subroutine equation_of_state_potential_mellor91

!....................................................................................

      subroutine equation_of_state_pressure_mellor91(time_)
      use module_principal
      implicit none
      integer time_ ! after, now or before

! Compute rh - rhp (rh is the full density, rhp is the potential density)  i.e. the
! part of the density which is only due to the compression effect

! We use the method described in:
! Mellor (1991) Journal of Atmospheric and Oceanic Technology. pp 609-611

! anyv3d(i,j,k,0) will be, temporarily, the part of the density related to pressure
! (see also presgrad.F90 in which potential density and "pressure density" are used
! in different ways).

      const3=+1.e4
      do j=1,jmax
      do i=1,imax
      if(  mask_t(i,j,kmax+1).eq.1)then !>>>>>>>>>>>>>>>>>>>>
      do k=kmin_w(i,j),kmax

      const2=-depth_t(i,j,k)/                                     &
           ( 1449.2+1.34*(sal_t(i,j,k,time_)-35.)              &
                    +4.55*tem_t(i,j,k,time_)                   &
                   -0.045*tem_t(i,j,k,time_)**2                &  !21-07-10
               -0.00821*depth_t(i,j,k)                            &
                +15.d-9*depth_t(i,j,k)**2  )**2

      anyv3d(i,j,k,0)=const3*const2*(1.-0.2*const2)


      enddo
      endif                         !>>>>>>>>>>>>>>>>>>>>

      enddo
      enddo

      end subroutine equation_of_state_pressure_mellor91

!....................................................................................

      subroutine equation_of_state_potential_opa(time_)
      use module_principal
      implicit none
!.....
      double precision zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze,zbw,zb &
       ,zd,zc,zaw,za,zb1,za1,zkw,zk0
!.....
      integer time_ ! after, now or before

! Compute (density - rho in kg/m3) at the pressure of reference as in opa

! rho is a reference density adjusted to the initial state in routine initial_state_eq.F90

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      zt = tem_t(i,j,k,time_)
      zs = sal_t(i,j,k,time_)

      ! square root salinity
      zsr=sqrt(abs(  sal_t(i,j,k,time_) ))

      ! compute volumic mass pure water at atm pressure
      zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
                -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
      ! seawater volumic mass atm pressure
      zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
             -4.0899e-3 ) *zt+0.824493
      zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
      zr4= 4.8314e-4

      ! potential volumic mass (reference to the surface)
      rhp_t(i,j,k) = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1 - rho


      enddo
      enddo
      enddo

      end subroutine equation_of_state_potential_opa

!....................................................................................

      subroutine equation_of_state_pressure_opa(time_)
      use module_principal
      implicit none
!.....
      double precision zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze,zbw,zb &
       ,zd,zc,zaw,za,zb1,za1,zkw,zk0
!.....
      integer time_ ! after, now or before

! Compute rh - rhp (rh is the full density, rhp is the potential density)  i.e. the
! part of the density which is only due to the compression effect

! We use the opa scheme

! - rho is a reference density adjusted to the initial state in routine initial_state_eq.F90
! - rhp_t+rho is the potential density (kg/m3)
! - anyv3d(i,j,k,0) will be, temporarily, the part of the density related to pressure
! (see also presgrad.F90 in which potential density and "pressure density" are used
! in different ways).

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      zt = tem_t(i,j,k,time_)
      zs = sal_t(i,j,k,time_)

      ! depth
      zh =-depth_t(i,j,k)

      ! square root salinity
      zsr=sqrt(abs(  sal_t(i,j,k,time_) ))

      ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
      zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
      zb = zbw + ze * zs

      zd = -2.042967e-2
      zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
      zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
      za = ( zd*zsr + zc ) *zs + zaw


      zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
      za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
      zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
      zk0= ( zb1*zsr + za1 )*zs + zkw

      anyv3d(i,j,k,0)=(rhp_t(i,j,k)+rho)/(1.0-zh/(zk0-zh*(za-zh*zb))) &! full density
                     -(rhp_t(i,j,k)+rho)                               !-potential density

      enddo
      enddo
      enddo

      end subroutine equation_of_state_pressure_opa

!..............................................................................

      subroutine equation_of_state_0d
      use module_principal
      implicit none


      if(sal_0d.lt.0.)then !--------->                                 !28/01/05
       write(6,*)
       write(6,*)'equation d''etat non lineaire ne fonctionne pas pour'
       write(6,*)'salinités négatives donc '
       stop 'density.f cas 6'
      endif                !--------->                                 !28/01/05

      const5=-rho+999.842594 ! 1er coef de XA - constante

      x2=tem_0d**2
      x3=tem_0d**3
      x4=tem_0d**4

! Calcul de la densite potentielle:

      rho_0d=                                                           &
! XA-RH0:
       const5+                                                          &
       6.793952e-2*tem_0d-                                              &
       9.095290e-3*x2+                                                  &
       1.001685e-4*x3-                                                  &
       1.120083e-6*x4+                                                  &
       6.536332e-9*tem_0d**5                                            &

! +XB*S:
       +(                                                               &
         8.24493e-1-                                                    &
       4.0899e-3*tem_0d+                                                &
       7.6438e-5*x2-                                                    &
       8.2467e-7*x3+                                                    &
       5.3875e-9*x4                                                     &
        )*sal_0d                                                        &

! +XC*S**1.5:
       +(                                                               &
         -5.72466e-3+                                                   &
       1.0227e-4*tem_0d-                                                &
       1.6546e-6*x2                                                     &
        )*sal_0d**1.5                                                   &

! +XD*S**2:
       +4.8314e-4*sal_0d**2

      end subroutine equation_of_state_0d

!.....................................................................

      subroutine equation_of_state_nopressure
      use module_principal
      implicit none

! If no compression, density = potential density:
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
      anyv3d(i,j,k,0)=0.
      enddo
      enddo
      enddo

      end subroutine equation_of_state_nopressure

!.....................................................................

      subroutine equation_of_state_coef
      use module_principal
      implicit none

! Compute eos_ref, a coefficient of a linear relation approximating the
! compression effect in order to remove a z-reference profil from the
! 3d density field before the computation of the pressure gradient force.
! See also routines equation_of_state_pressure_
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

       sum1=sum1+mask_t(i,j,k)*depth_t(i,j,k)
       sum2=sum2+mask_t(i,j,k)*(anyv3d(i,j,k,0)-rhp_t(i,j,k))

      enddo
      enddo
      enddo
      if(sum1==0.)stop 'STOP sum1=0 in routine equation_of_state_coef'
!     eos_coef=sum2/sum1

      end subroutine equation_of_state_coef

!................................................................................

      subroutine equation_of_state_potential_jmfwg(time_)
      use module_principal
      use module_parallele
      implicit none
      double precision tem_,tem2_,tem3_,tem4_    &
                      ,sal_,sal2_
      integer time_ ! after, now or before

! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3


      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=tem_t(i,j,k,time_)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

      sal_=max(sal_t(i,j,k,time_),zero)
      sal2_=sal_*sal_

       rhp_t(i,j,k)=                                          &
          (c1_jmfwg+c2_jmfwg *tem_+c3_jmfwg*tem2_             &
                   +c4_jmfwg *tem3_                           &
                   +c5_jmfwg *sal_                            &
                   +c6_jmfwg *sal_*tem_                       &
                   +c7_jmfwg *sal2_ )                         &
           /(1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_           &
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_           &
                   +c18_jmfwg*sal_                            &
                   +c19_jmfwg*tem_*sal_                       &
                   +c20_jmfwg*sal_*tem3_                      &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_))  &
           -rho

      enddo
      enddo
      enddo

      end subroutine equation_of_state_potential_jmfwg

!................................................................................

      subroutine equation_of_state_pressure_jmfwg(time_)
      use module_principal
      implicit none
      double precision tem_,tem2_,tem3_,tem4_    &
                      ,sal_,sal2_,b0_,a0_,prs_,prs2_ &
                      ,prs3_
      integer time_ ! after, now or before


! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3


!     do k=1,kmax
      do j=1,jmax
      do i=1,imax

      x0=0.
      do k=kmax,1,-1

      x0=x0+h_w(i,j)*dsig_t(i,j,k)
!     prs_=-depth_t(i,j,k)
      prs_=x0-0.5*dsig_t(i,j,k)*h_w(i,j)
      prs2_= prs_*prs_
      prs3_=prs2_*prs_
       

      tem_=tem_t(i,j,k,time_)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

      sal_=max(sal_t(i,j,k,time_),zero)
      sal2_=sal_*sal_


      b0_=                                                      &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                              &
                   +c5_jmfwg*sal_                               &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_          & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                              &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


      anyv3d(i,j,k,0)=                                              &

      ( b0_                                                      &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prs_       &
      +(c11_jmfwg+c12_jmfwg*tem2_)*                prs2_  )   &

      /( a0_                                                     &
        +c23_jmfwg*         prs_                                 &
        +c24_jmfwg*tem3_*prs2_                                &
        +c25_jmfwg*tem_ *prs3_   )                            &

       -b0_/a0_


      enddo
      enddo
      enddo

      end subroutine equation_of_state_pressure_jmfwg

!................................................................................

      subroutine equation_of_state_pressure_w97(time_)
      use module_principal
      implicit none
      double precision       &
       l_,a0_,p0_,rho_,prs_,tem_,sal_,tem2_    &
      ,tem3_
!                      tem_,tem2_,tem3_,tem4_    &
!                     ,sal_,sal2_,b0_,a0_,prs_,prs2_ &
!                     ,prs3_,
      integer time_ ! after, now or before

! Wright (JAOT 1997) equation of state:
! coefficients are those of its TABLE 1 column: "Extended formula fit over the reduced range"


      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=tem_t(i,j,k,time_)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_

!     sal_=max(sal_t(i,j,k,time_),zero)
      sal_=sal_t(i,j,k,time_)

! Pression en Pa:
      prs_=-rho*grav*depth_t(i,j,k)

! Coef alpha0:
      a0_=a0_w97+a1_w97*tem_+a2_w97*sal_
! Coef p0:
      p0_=b0_w97+b1_w97*tem_+b2_w97*tem2_+b3_w97*tem3_   &
                   +b4_w97*sal_+b5_w97*tem_*sal_
! Coef lambda:
      l_= c0_w97+c1_w97*tem_+c2_w97*tem2_+c3_w97*tem3_   &
                   +c4_w97*sal_+c5_w97*tem_*sal_

! rh(T,S,p)-rh(T,S,0):
      anyv3d(i,j,k,0)=(prs_+p0_)/(l_+a0_*(prs_+p0_)) &
                              -p0_ /(l_+a0_*         p0_)

      enddo
      enddo
      enddo


      end subroutine equation_of_state_pressure_w97

!................................................................................

      subroutine equation_of_state_potential_w97(time_)
      use module_principal
      implicit none
      double precision       &
       l_,a0_,p0_,rho_,prs_,tem_,sal_,tem2_    &
      ,tem3_
      integer time_ ! after, now or before

! Wright (JAOT 1997) equation of state:
! coefficients are those of its TABLE 1 column: "Extended formula fit over the reduced range"


      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=tem_t(i,j,k,time_)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_

!     sal_=max(sal_t(i,j,k,time_),zero)
      sal_=sal_t(i,j,k,time_)

! Pression en Pa:
      prs_=-rho*grav*depth_t(i,j,k)

! Coef alpha0:
      a0_=a0_w97+a1_w97*tem_+a2_w97*sal_
! Coef p0:
      p0_=b0_w97+b1_w97*tem_+b2_w97*tem2_+b3_w97*tem3_   &
                   +b4_w97*sal_+b5_w97*tem_*sal_
! Coef lambda:
      l_= c0_w97+c1_w97*tem_+c2_w97*tem2_+c3_w97*tem3_   &
                   +c4_w97*sal_+c5_w97*tem_*sal_

! potential density anomaly:
      rhp_t(i,j,k)=p0_ /(l_+a0_*p0_) - rho

      enddo
      enddo
      enddo


      end subroutine equation_of_state_potential_w97

!..........................................................................
! Quelques lignes pour passer de la temperature à la temperature
! potentielle:
!.....Pression hydrostatique
!---->Constantes
!     P0=1.01325
!     P=P0+RHO*G*(-ZW)/1E5
! où ZW est la profondeur en metre (prof3D_z par ex)
! soit la temperature et la salinitie TPC et SPC
! la temperature potentielle est TPC+TPC1+TPC2+TPC3 avec:
!       TPC1=-P*(3.6504E-4+
!    1   8.3198E-5*TPC-
!    2   5.4065E-7*TPC**2+
!    3   4.0274E-9*TPC**3)
!     TPC2=-P*(SPC-35.)*
!    1     (1.7439E-5-
!    2     2.9778E-7*TPC)-
!    3     P**2*(8.9309E-7-
!    4     3.1628E-8*TPC+
!    5     2.1987E-10*TPC**2)
!     TPC3=4.1057E-9*(SPC-35.)*P**2-
!    1     P**3*(-1.6056E-10+
!    2     5.0484E-12*TPC)
!..........................................................................



!..........................................................................
! Quelques lignes pour passer de la temperature à la temperature
! potentielle:
!.....Pression hydrostatique
!---->Constantes
!     P0=1.01325
!     P=P0+RHO*G*(-ZW)/1E5
! où ZW est la profondeur en metre (prof3D_z par ex)
! soit la temperature et la salinitie TPC et SPC
! la temperature potentielle est TPC+TPC1+TPC2+TPC3 avec:
!       TPC1=-P*(3.6504E-4+
!    1   8.3198E-5*TPC-
!    2   5.4065E-7*TPC**2+
!    3   4.0274E-9*TPC**3)
!     TPC2=-P*(SPC-35.)*
!    1     (1.7439E-5-
!    2     2.9778E-7*TPC)-
!    3     P**2*(8.9309E-7-
!    4     3.1628E-8*TPC+
!    5     2.1987E-10*TPC**2)
!     TPC3=4.1057E-9*(SPC-35.)*P**2-
!    1     P**3*(-1.6056E-10+
!    2     5.0484E-12*TPC)
!..........................................................................


!................................................................................

      subroutine equation_of_state_potloc_jmfwg(txt_,timeobc_)    !16-05-14
      use module_principal
      use module_parallele
      implicit none
      double precision tem_,tem2_,tem3_,tem4_    &
                      ,sal_,sal2_,b0_,a0_        &
                      ,prsup_,prsup2_,prsup3_    &
                      ,prdwn_,prdwn2_,prdwn3_
!                     ,rhsup_,rhdwn_
      character*3 txt_
      integer timeobc_, time_ !29-10-22

      time_=min(timeobc_,2) ! parce que possibilite timeobc_=4 le 29-10-22

! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3

! Calcule la densite complete "en double" c.a.d. sur deux niveaux de pression
! (au dessus au dessous) afin d'oter par differentiation l'effet de compression
! et de fournir une densite potentielle a reference locale au lieu de la
! densite potentielle referencee a une pression constante (par ex 0, ou autre...)
! Complements d'explications dans:
! https://docs.google.com/document/d/1DPvKOZ3rV-FQfEaYsX-bTDwVmDg3wB4kbPcdS1TJqoc/edit

      ksecu=0
         if(txt_=='pgf'     &
        .or.txt_=='tke'     & !19-06-14
        .or.txt_=='grp') then !23-02-22
       k1=1 ; k2=0 ; ksecu=1 ! select tem_t(time_)  & sal_t(time_)
      endif
      if(txt_=='obc') then
       k1=0 ; k2=1 ; ksecu=1 ! select temobc_t(time_) & salobc_t(time_)
      endif
      if(txt_=='ofl') then
       k1=0 ; k2=0 ; ksecu=1 ! select temofl_t & salofl_t
      endif
      if(ksecu==0)stop 'txt_ inconnu equation_of_state_potloc_jmfwg'

      if(txt_=='ofl') then !oooooo> !18-10-14

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=temofl_t(i,j,k,1)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

      sal_=max(salofl_t(i,j,k,1),zero)
      sal2_=sal_*sal_

! prsup = Pression facette superieure
      prsup_ =-depth_w(i,j,k+1)
      prsup2_= prsup_*prsup_
      prsup3_=prsup2_*prsup_

! prdwn = Pression facette inferieure
      prdwn_ =-depth_w(i,j,k)
      prdwn2_= prdwn_*prdwn_
      prdwn3_=prdwn2_*prdwn_

      b0_=                                                   &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                           &
                   +c5_jmfwg*sal_                            &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_             & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                           &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


!     Densite facette superieure:
      anyv3d(i,j,k,3)=                                       &

      ( b0_                                                  &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prsup_       &
      +(c11_jmfwg+c12_jmfwg*tem2_)*             prsup2_  )   &

      /( a0_                                                 &
        +c23_jmfwg*      prsup_                              &
        +c24_jmfwg*tem3_*prsup2_                             &
        +c25_jmfwg*tem_ *prsup3_   )

!     Densite facette inferieure:
      anyv3d(i,j,k,2)=                                       &

      ( b0_                                                  &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prdwn_       &
      +(c11_jmfwg+c12_jmfwg*tem2_)*             prdwn2_  )   &

      /( a0_                                                 &
        +c23_jmfwg*      prdwn_                              &
        +c24_jmfwg*tem3_*prdwn2_                             &
        +c25_jmfwg*tem_ *prdwn3_   )

      enddo
      enddo
      enddo

      else               !oooooo>

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=k1*   tem_t(i,j,k,time_)      &
          +k2*temobc_t(i,j,k,timeobc_)  !29-10-22     
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

      sal_=max( k1*   sal_t(i,j,k,time_)      &
               +k2*salobc_t(i,j,k,timeobc_),0.) !29-10-22
      sal2_=sal_*sal_

! prsup = Pression facette superieure
      prsup_ =-depth_w(i,j,k+1)
      prsup2_= prsup_*prsup_
      prsup3_=prsup2_*prsup_

! prdwn = Pression facette inferieure
      prdwn_ =-depth_w(i,j,k)
      prdwn2_= prdwn_*prdwn_
      prdwn3_=prdwn2_*prdwn_

      b0_=                                                   &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                           &
                   +c5_jmfwg*sal_                            &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_             & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                           &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


!     Densite facette superieure:
      anyv3d(i,j,k,3)=                                       &

      ( b0_                                                  &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prsup_       &
      +(c11_jmfwg+c12_jmfwg*tem2_)*             prsup2_  )   &

      /( a0_                                                 &
        +c23_jmfwg*      prsup_                              &
        +c24_jmfwg*tem3_*prsup2_                             &
        +c25_jmfwg*tem_ *prsup3_   )

!     Densite facette inferieure:
      anyv3d(i,j,k,2)=                                       &

      ( b0_                                                  &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prdwn_       &
      +(c11_jmfwg+c12_jmfwg*tem2_)*             prdwn2_  )   &

      /( a0_                                                 &
        +c23_jmfwg*      prdwn_                              &
        +c24_jmfwg*tem3_*prdwn2_                             &
        +c25_jmfwg*tem_ *prdwn3_   )

      enddo
      enddo
      enddo

      endif              !oooooo>

!..................

      do j=1,jmax
      do i=1,imax
      xy_t(i,j,1)=0.
      enddo
      enddo
      if(txt_=='pgf'.or.txt_=='tke') then !vvvvvvvvv>

! Cas txt_='pgf' ou 'ofl' densite potentielle dans rhp_t(i,j,k)
       do k=kmax,1,-1 ; do j=1,jmax ; do i=1,imax

! x1 = delta_rho_p = Increment local de densite due a la pression seule:
        x1=anyv3d(i,j,k,3)-anyv3d(i,j,k,2)
! Cumul vertical de delta_rho_p:
        xy_t(i,j,1)=xy_t(i,j,1)+x1

! Anomalie de densite potentielle:
           rhp_t(i,j,k)=0.5*(anyv3d(i,j,k,3)+anyv3d(i,j,k,2))   & !   densite totale au centre de la cellule
                       +xy_t(i,j,1)-0.5*x1                      & ! - cumul de l'effet de compression (note: methode rectangle)
                       -rho                                       ! - constante (de l'ordre de 1020, 1030, ...)

       enddo ; enddo ; enddo

      else                                !vvvvvvvvv>

! Cas txt_='obc' ou 'ofl' densite potentielle dans anyv3d(i,j,k,1)
       do k=kmax,1,-1 ; do j=1,jmax ; do i=1,imax

! x1 = delta_rho_p = Increment local de densite due a la pression seule:
        x1=anyv3d(i,j,k,3)-anyv3d(i,j,k,2)
! Cumul vertical de delta_rho_p:
        xy_t(i,j,1)=xy_t(i,j,1)+x1

! Anomalie de densite potentielle:
        anyv3d(i,j,k,1)=0.5*(anyv3d(i,j,k,3)+anyv3d(i,j,k,2))   & !   densite totale au centre de la cellule
                       +xy_t(i,j,1)-0.5*x1                      & ! - cumul de l'effet de compression (note: methode rectangle)
                       -rho                                       ! - constante (de l'ordre de 1020, 1030, ...)

       enddo ; enddo ; enddo

      endif                !vvvvvvvvv>


      end subroutine equation_of_state_potloc_jmfwg

!................................................................................

      subroutine equation_of_state_potloc_w97(time_) !16-05-14
      use module_principal
      implicit none
      double precision       &
       l_,a0_,p0_,rho_,prsup_,tem_,sal_,tem2_    &
      ,tem3_,prdwn_
      integer time_ ! after, now or before

! Wright (JAOT 1997) equation of state:
! coefficients are those of its TABLE 1 column: "Extended formula fit over the reduced range"


      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=tem_t(i,j,k,time_)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_

!     sal_=max(sal_t(i,j,k,time_),zero)
      sal_=sal_t(i,j,k,time_)

! Pression en Pa:
!     prsup_=-rho*grav*depth_t(i,j,k)
      prsup_=-rho*grav*depth_w(i,j,k+1)
      prdwn_=-rho*grav*depth_w(i,j,k  )

! Coef alpha0:
      a0_=a0_w97+a1_w97*tem_+a2_w97*sal_
! Coef p0:
      p0_=b0_w97+b1_w97*tem_+b2_w97*tem2_+b3_w97*tem3_   &
                   +b4_w97*sal_+b5_w97*tem_*sal_
! Coef lambda:
      l_= c0_w97+c1_w97*tem_+c2_w97*tem2_+c3_w97*tem3_   &
                   +c4_w97*sal_+c5_w97*tem_*sal_

! rh(T,S,prsup):
      anyv3d(i,j,k,2)=(prsup_+p0_)/(l_+a0_*(prsup_+p0_))

! rh(T,S,prdwn):
      anyv3d(i,j,k,1)=(prdwn_+p0_)/(l_+a0_*(prdwn_+p0_))

      enddo
      enddo
      enddo


!..................

      do j=1,jmax
      do i=1,imax
      xy_t(i,j,1)=0.
      enddo
      enddo
      do k=kmax,1,-1
      do j=1,jmax
      do i=1,imax

! x1 = delta_rho_p = Increment local de densite due a la pression seule:
      x1=anyv3d(i,j,k,2)-anyv3d(i,j,k,1)

! Cumul vertical de delta_rho_p:
      xy_t(i,j,1)=xy_t(i,j,1)+x1

! Anomalie de densite potentielle:
      rhp_t(i,j,k)=0.5*(anyv3d(i,j,k,2)+anyv3d(i,j,k,1))   & !   densite totale au centre de la cellule
                  +xy_t(i,j,1)-0.5*x1                      & ! - cumul de l'effet de compression (note: methode rectangle)
                  -rho                                       ! - constante (de l'ordre de 1020, 1030, ...)
      enddo
      enddo
      enddo

      end subroutine equation_of_state_potloc_w97

!................................................................................

      subroutine equation_of_state_potzref_jmfwg(txt_,time_)
      use module_principal
      implicit none
      double precision tem_,tem2_,tem3_,tem4_    &
                      ,sal_,sal2_,b0_,a0_,prs_,prs2_ &
                      ,prs3_,zref_
      integer time_     ! after, now or before
      character*3 txt_

! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3

      zref_=max(eos_tkezref,zero) ! Si eos_tkezref<0 (zref local) alors zref=0 !20-10-14
                                  ! utile pour ecriture densite potentielle dans fichiers netcdf

      if(    txt_=='pgf' &
         .or.txt_=='tke') then !vvvvvvvvvvvvvvvvvvvvvv> !17-05-14

      if(txt_=='pgf')zref_=eos_pgfzref !20-06-14

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=tem_t(i,j,k,time_)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

!     prs_=-depth_t(i,j,k)
      prs_=zref_
      prs2_= prs_*prs_
      prs3_=prs2_*prs_

      sal_=max(sal_t(i,j,k,time_),zero)
      sal2_=sal_*sal_


      b0_=                                                   &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                           &
                   +c5_jmfwg*sal_                            &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_             & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                           &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


      rhp_t(i,j,k)=                                           &

      ( b0_                                                   &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prs_          &
      +(c11_jmfwg+c12_jmfwg*tem2_)*                prs2_  )   &

      /( a0_                                                  &
        +c23_jmfwg*      prs_                                 &
        +c24_jmfwg*tem3_*prs2_                                &
        +c25_jmfwg*tem_ *prs3_   )                            &

       -rho

      enddo
      enddo
      enddo
      return
      endif                !vvvvvvvvvvvvvvvvvvvvvv>

      if(txt_=='obc'.or.    &
         txt_=='ofl'.or.    &
         txt_=='grp') then !oooooooooooooooooooooo>

      if(txt_=='obc')then
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,2)=temobc_t(i,j,k,time_) !22-05-14
        anyv3d(i,j,k,3)=salobc_t(i,j,k,time_)
       enddo       ; enddo       ; enddo
      endif
      if(txt_=='ofl')then
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,2)=temofl_t(i,j,k,1) !22-05-14
        anyv3d(i,j,k,3)=salofl_t(i,j,k,1)
       enddo       ; enddo       ; enddo
      endif
      if(txt_=='grp')then                 !20-10-14
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,2)=tem_t(i,j,k,1) 
        anyv3d(i,j,k,3)=sal_t(i,j,k,1)
       enddo       ; enddo       ; enddo
      endif

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_=anyv3d(i,j,k,2) !22-05-14
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

!     prs_=-depth_t(i,j,k)
      prs_=zref_
      prs2_= prs_*prs_
      prs3_=prs2_*prs_

      sal_=anyv3d(i,j,k,3) !22-05-14
      sal_=max(sal_,zero)
      sal2_=sal_*sal_


      b0_=                                                   &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                           &
                   +c5_jmfwg*sal_                            &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_             & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                           &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


!     rhp_t(i,j,k)=                                           &
      anyv3d(i,j,k,1)=                                        &

      ( b0_                                                   &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prs_          &
      +(c11_jmfwg+c12_jmfwg*tem2_)*                prs2_  )   &

      /( a0_                                                  &
        +c23_jmfwg*         prs_                              &
        +c24_jmfwg*tem3_*prs2_                                &
        +c25_jmfwg*tem_ *prs3_   )                            &

       -rho

      enddo
      enddo
      enddo
      return
      endif                !oooooooooooooooooooooo>

      stop 'txt non defini subroutine equation_of_state_potzref_jmfwg'

      end subroutine equation_of_state_potzref_jmfwg

!................................................................................

      subroutine equation_of_state_potloc2 !19-06-14
      use module_principal
      implicit none
      double precision prsup_,prsup2_,prsup3_,prdwn_,prdwn2_,prdwn3_

! SAME AS equation_of_state_potzref_jmfwg but the EOS T,S dependent coeficients
! are already stored in "anyv3d" arrays, following the pressure gradient routine
! pressure_gradient_jmfwg06_hyb

! Note de l'auteur: j'ai verifie:
! 1 - que le gradient vertical de rhp etait localement
! egal au gradient vertical de rhp calculee "a une profondeur de reference" quand
! z est egal a la profondeur de reference en question.
! 2 - que le gradient vertical est nul si le gradient vertical de T et S sont nuls

! On peut a la place l'algo de McDougall JPO 1987 qui est equivalent:
!     call equation_of_state_rhpmdjpo87(1) !18-10-14
!     return

      i=imax/2 ; j=jmax/2 ; k=kmax/2 !16-02-17
      if(checkanyv3d(1)/=anyv3d(i,j,k,1))stop ' Err 1 eosA anyv3d '
      if(checkanyv3d(2)/=anyv3d(i,j,k,2))stop ' Err 2 eosA anyv3d '
      if(checkanyv3d(3)/=anyv3d(i,j,k,3))stop ' Err 3 eosA anyv3d '
      if(checkanyv3d(4)/=anyv3d(i,j,k,4))stop ' Err 4 eosA anyv3d '
      if(checkanyv3d(6)/=anyv3d(i,j,k,6))stop ' Err 6 eosA anyv3d '
      if(checkanyv3d(7)/=anyv3d(i,j,k,7))stop ' Err 7 eosA anyv3d '

! Densite potentielle "z local"
      do j=0,jmax+1 ; do i=0,imax+1 !26-01-23
       xy_t(i,j,1)=0.
      enddo ; enddo

      do k=kmax,1,-1
      do j=0,jmax+1 ; do i=0,imax+1 !26-01-23

! prsup = Pression facette superieure
      prsup_ =-depth_w(i,j,k+1)
      prsup2_= prsup_*prsup_
      prsup3_=prsup2_*prsup_

! prdwn = Pression facette inferieure
      prdwn_ =-depth_w(i,j,k)
      prdwn2_= prdwn_*prdwn_
      prdwn3_=prdwn2_*prdwn_

! Densite pression facette superieure:
          x2=( anyv3d(i,j,k,1)             &
              +anyv3d(i,j,k,2)*prsup_      &
              +anyv3d(i,j,k,3)*prsup2_ )   &
            /( anyv3d(i,j,k,4)             &
                    +c23_jmfwg*prsup_      &
              +anyv3d(i,j,k,6)*prsup2_     &
              +anyv3d(i,j,k,7)*prsup3_ )

! Densite pression facette inferieure:
          x1=( anyv3d(i,j,k,1)             &
              +anyv3d(i,j,k,2)*prdwn_      &
              +anyv3d(i,j,k,3)*prdwn2_ )   &
            /( anyv3d(i,j,k,4)             &
                    +c23_jmfwg*prdwn_      &
              +anyv3d(i,j,k,6)*prdwn2_     &
              +anyv3d(i,j,k,7)*prdwn3_ )

! x2-x1 = delta_rho_p = Increment local de densite due a la pression seule:

! Cumul vertical de delta_rho_p:
      xy_t(i,j,1)=xy_t(i,j,1)+x2-x1

! Anomalie de densite potentielle:
      rhp_t(i,j,k)=0.5*(x1+x2)               & !   densite totale au centre de la cellule
                  +xy_t(i,j,1)-0.5*(x2-x1)   & ! - cumul de l'effet de compression (note: methode rectangle)
                  -rho                         ! - constante (de l'ordre de 1020, 1030, ...)

      enddo ; enddo
      enddo

      end subroutine equation_of_state_potloc2

!................................................................................

      subroutine equation_of_state_rhp_a1_a7(time_)
      use module_principal
      implicit none
      double precision tem_,sal_
      integer time_


      do k=kmax,1,-1
      do j=0,jmax+1 ; do i=0,imax+1 !26-01-23

!     tem_=tem_t(i,j,k,time_)
!     sal_=max(sal_t(i,j,k,time_),zero)
! Appliquer A l'EOS des bornes de validites sur T et S            !07-02-20
      tem_=min(max(tem_t(i,j,k,time_),-3.),50.) !07-02-20
      sal_=min(max(sal_t(i,j,k,time_), 0.),50.) !07-02-20

! Compute density rho_    from (sal_   ,tem_   ,prs_   ) according
! to Jackett et al 2006 (JAOT):
!     tem2_= tem_*tem_
!     tem3_=tem2_*tem_
!     tem4_=tem3_*tem_
!     sal2_= sal_*sal_

      anyv3d(i,j,k,1)= c1_jmfwg                 &
                      +c2_jmfwg*tem_            &
                      +c3_jmfwg*tem_**2         &
                      +c4_jmfwg*tem_**3         &
                      +c5_jmfwg*sal_            &
                      +c6_jmfwg*sal_*tem_       &
                      +c7_jmfwg*sal_**2

      anyv3d(i,j,k,2)= c8_jmfwg                 &
                      +c9_jmfwg*tem_**2         &
                     +c10_jmfwg*sal_

      anyv3d(i,j,k,3)=c11_jmfwg                 &
                     +c12_jmfwg*tem_**2

      anyv3d(i,j,k,4)=1.                        &
                     +c14_jmfwg*tem_            &
                     +c15_jmfwg*tem_**2         &
                     +c16_jmfwg*tem_**3         &
                     +c17_jmfwg*tem_**4         &
                     +c18_jmfwg*sal_            &
                     +c19_jmfwg*tem_*sal_       &
                     +c20_jmfwg*sal_*tem_**3    &
                     +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem_**2)

!     anyv3d(i,j,k,5)=c23_jmfwg

      anyv3d(i,j,k,6)=c24_jmfwg*tem_**3

      anyv3d(i,j,k,7)=c25_jmfwg*tem_

! Potential density anomaly:
      rhp_t(i,j,k)=anyv3d(i,j,k,1)/anyv3d(i,j,k,4) - rho

      enddo ; enddo
      enddo

      i=imax/2 ; j=jmax/2 ; k=kmax/2 !16-02-17
      checkanyv3d(1)=anyv3d(i,j,k,1)
      checkanyv3d(2)=anyv3d(i,j,k,2)
      checkanyv3d(3)=anyv3d(i,j,k,3)
      checkanyv3d(4)=anyv3d(i,j,k,4)
      checkanyv3d(6)=anyv3d(i,j,k,6)
      checkanyv3d(7)=anyv3d(i,j,k,7)

      end subroutine equation_of_state_rhp_a1_a7

!................................................................................

      subroutine equation_of_state_potzref2 !20-06-14
      use module_principal
      implicit none
      double precision prs_,prs2_,prs3_

! SAME AS equation_of_state_potzref_jmfwg but the EOS T,S dependent coeficients
! are already stored in "anyv3d" arrays, following the pressure gradient routine
! pressure_gradient_jmfwg06_hyb

      i=imax/2 ; j=jmax/2 ; k=kmax/2 !16-02-17
      if(checkanyv3d(1)/=anyv3d(i,j,k,1))stop ' Err 1 eosA anyv3d '
      if(checkanyv3d(2)/=anyv3d(i,j,k,2))stop ' Err 2 eosA anyv3d '
      if(checkanyv3d(3)/=anyv3d(i,j,k,3))stop ' Err 3 eosA anyv3d '
      if(checkanyv3d(4)/=anyv3d(i,j,k,4))stop ' Err 4 eosA anyv3d '
      if(checkanyv3d(6)/=anyv3d(i,j,k,6))stop ' Err 6 eosA anyv3d '
      if(checkanyv3d(7)/=anyv3d(i,j,k,7))stop ' Err 7 eosA anyv3d '


! McDougall JPO 1987:
      if(eos_tkezref<0.) &
      stop ' eos_tkezref<0 in equation_of_state_potzref2'

      prs_ =eos_tkezref
      prs2_= prs_*prs_
      prs3_=prs2_*prs_

      do k=kmax,1,-1
      do j=0,jmax+1 ; do i=0,imax+1 !26-01-23

! Densite pression facette superieure:
       rhp_t(i,j,k)= ( anyv3d(i,j,k,1)              &
                      +anyv3d(i,j,k,2)*prs_         &
                      +anyv3d(i,j,k,3)*prs2_ )      &
                    /( anyv3d(i,j,k,4)              &
                            +c23_jmfwg*prs_         &
                      +anyv3d(i,j,k,6)*prs2_        &
                      +anyv3d(i,j,k,7)*prs3_ )-rho

      enddo ; enddo
      enddo


      end subroutine equation_of_state_potzref2

!................................................................................

      subroutine equation_of_state_rhpmdjpo87(t_) !18-10-14
      use module_principal
      implicit none
      double precision tem_,sal_,prs_,alphaobeta_,beta_
      integer t_

! Si on commente ces lignes rhp(kmax) reste sur la derniere valeur calculee
! au moment du calcul du PGF
!     do j=1,jmax
!     do i=1,imax
!      rhp_t(i,j,kmax)=0.
!     enddo
!     enddo

! Neutral Surfaces, Trevor J. McDougall, JPO 1987:
! http://dx.doi.org/10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

      do k=kmax  ,2,-1
      do j=1,jmax
      do i=1,imax

      tem_=0.5*(tem_t(i,j,k,t_)+tem_t(i,j,k-1,t_))         ! potential temperature at depth_w
      sal_=0.5*(sal_t(i,j,k,t_)+sal_t(i,j,k-1,t_)) - 35.0  ! salinity anomaly (s-35) at depth_w
      prs_=-depth_w(i,j,k)                                 ! pressure in db approximated to -depth in meters

! Test values : beta_=0.72088e-3 alphaobeta_=0.34763
!               tem_=10. ; sal_=40.-35. ; prs_=4000.

             alphaobeta_ = ( ( ( - 0.255019e-07 *tem_+ 0.298357e-05 )*tem_  &   ! ratio alpha/beta
                     &                               - 0.203814e-03 )*tem_  &
                     &                               + 0.170907e-01 )*tem_  &
                     &   +         0.665157e-01                     &
                     &   +     ( - 0.678662e-05 *sal_               &
                     &           - 0.846960e-04 *tem_+ 0.378110e-02 )*sal_  &
                     &   +   ( ( - 0.302285e-13 *prs_               &
                     &           - 0.251520e-11 *sal_               &
                     &           + 0.512857e-12 *tem_*tem_ ) *prs_  &
                     &           - 0.164759e-06 *sal_               &
                     &        +(   0.791325e-08 *tem_- 0.933746e-06 )*tem_  &
                     &                               + 0.380374e-04 )*prs_
                     !
                     beta_  = ( ( -0.415613e-09 *tem_+ 0.555579e-07 )*tem_     &   ! beta
                     &                               - 0.301985e-05 )*tem_     &
                     &   +       0.785567e-03                       &
                     &   + (     0.515032e-08 *sal_                 &
                     &         + 0.788212e-08 *tem_- 0.356603e-06 )*sal_    &
                     &   + ( (   0.121551e-17 *prs_                 &
                     &         - 0.602281e-15 *sal_                 &
                     &         - 0.175379e-14 *tem_+ 0.176621e-12 )*prs_    &
                     &                             + 0.408195e-10  *sal_    &
                     &     + ( - 0.213127e-11 *tem_+ 0.192867e-09 )*tem_    &
                     &                             - 0.121555e-07 )*prs_

!#ifdef bidon
!        alphaobeta_ =0.665157e-1         +0.170907e-1*tem_       &
!                    -0.203814e-3*tem_**2 +0.298357e-5*tem_**3    &
!                    -0.255019e-7*tem_**4                         &
!             +sal_*( 0.378110e-2         -0.846960e-4*tem_       &
!                    -0.164759e-6*prs_    -0.251520e-11*prs_**2)  &
!        +(sal_**2)*(-0.678662e-5)                                &
!             +prs_*( 0.380374e-4         -0.933746e-6*tem_       &
!                    +0.791325e-8*tem_**2)                        &
!        +0.512857e-12*(prs_*tem_)**2-0.302285e-13*prs_**3
!      write(6,*)' 0.34763 ',alphaobeta_
!      stop 'toto'
!#endif


        rhp_t(i,j,k-1)=                        &
        rhp_t(i,j,k  )-                        &
        rho*beta_*( alphaobeta_*( tem_t(i,j,k-1,t_)-tem_t(i,j,k,t_) )   &
                              - ( sal_t(i,j,k-1,t_)-sal_t(i,j,k,t_) ) )

      enddo
      enddo
      enddo

! Test values : beta_=0.72088e-3 alphaobeta_=0.34763 when
!     write(6,*)' 0.72088e-3 ',beta_
!     write(6,*)' 0.34763    ',alphaobeta_

      end subroutine equation_of_state_rhpmdjpo87

!................................................................................
!#ifdef bidon
      subroutine equation_of_state_pot_anyv(istr_,iend_,jstr_,jend_)
      use module_principal
      use module_parallele
      implicit none
      double precision tem_,tem2_,tem3_,tem4_    &
                      ,sal_,sal2_,cst1_,cst2_,cst3_
      integer time_                     & ! after, now or before
             ,istr_,iend_,jstr_,jend_

! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3

! verifier la disponibilite de anyv3d(:,:,:,id_rhp) !30-09-15
      if(anyv3d(-1,-1,0,id_rhp)==-9999.) &
      stop 'Err anyv3d id_rhp not available'


      if(eos_author==3) then !333333>
       do k=1,kmax
       do j=jstr_,jend_
       do i=istr_,iend_

       tem_=anyv3d(i,j,k,id_tem)
       tem2_= tem_*tem_
       tem3_=tem2_*tem_
       tem4_=tem3_*tem_

       sal_=max(anyv3d(i,j,k,id_sal),zero)
       sal2_=sal_*sal_

       anyv3d(i,j,k,id_rhp)=                                  &
          (c1_jmfwg+c2_jmfwg *tem_+c3_jmfwg*tem2_             &
                   +c4_jmfwg *tem3_                           &
                   +c5_jmfwg *sal_                            &
                   +c6_jmfwg *sal_*tem_                       &
                   +c7_jmfwg *sal2_ )                         &
           /(1.    +c14_jmfwg*tem_ +c15_jmfwg*tem2_           &
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_           &
                   +c18_jmfwg*sal_                            &
                   +c19_jmfwg*tem_*sal_                       &
                   +c20_jmfwg*sal_*tem3_                      &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_))  &
           -rho

! Cas test exact lineaire:
!      anyv3d(i,j,k,id_rhp)=0.001*anyv3d(i,j,k,id_tem)

       enddo
       enddo
       enddo

      return
      endif                  !333333>

      if(eos_author==0) then !000000>
       cst1_=-rho*alp_t
       cst2_= rho*alp_s
       cst3_=-cst1_*t0-cst2_*s0
       do k=1,kmax
       do j=jstr_,jend_
       do i=istr_,iend_
        anyv3d(i,j,k,id_rhp)=cst1_*anyv3d(i,j,k,id_tem) &
                            +cst2_*anyv3d(i,j,k,id_sal) &
                            +cst3_
       enddo
       enddo
       enddo
      return
      endif                  !000000>

! Marque anyv3d(:,:,:,id_rhp) indisponible !30-09-15
      anyv3d(-1,-1,0,id_rhp)=-9999.
      
!#endif
      end subroutine equation_of_state_pot_anyv

!................................................................................

      subroutine equation_of_state_prs_anyv(istr_,iend_,jstr_,jend_)
      use module_principal
      use module_parallele
      implicit none
      double precision tem_,tem2_,tem3_,tem4_        &
                      ,sal_,sal2_,b0_,a0_,prs_,prs2_ &
                      ,prs3_
      integer istr_,iend_,jstr_,jend_

! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3


!     do k=1,kmax
      do j=jstr_,jend_
      do i=istr_,iend_

      x0=0.
      do k=kmax,1,-1
      x0=x0+h_w(i,j)*dsig_t(i,j,k)
!     prs_=-depth_t(i,j,k)
      prs_=x0-0.5*dsig_t(i,j,k)*h_w(i,j)
      prs2_= prs_*prs_
      prs3_=prs2_*prs_

      tem_=anyv3d(i,j,k,id_tem)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

      sal_=max(anyv3d(i,j,k,id_sal),zero)
      sal2_=sal_*sal_


      b0_=                                                   &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                           &
                   +c5_jmfwg*sal_                            &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.       +c14_jmfwg*tem_ +c15_jmfwg*tem2_          & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                           &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


      anyv3d(i,j,k,id_rhc)=                                   &

      ( b0_                                                   &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prs_          &
      +(c11_jmfwg+c12_jmfwg*tem2_)*                prs2_  )   &

      /( a0_                                                  &
        +c23_jmfwg*         prs_                              &
        +c24_jmfwg*tem3_*prs2_                                &
        +c25_jmfwg*tem_ *prs3_   )                            &

       -b0_/a0_

      enddo
      enddo
      enddo


      end subroutine equation_of_state_prs_anyv

!................................................................................

      subroutine equation_of_state_full_a1_a7(time_) !14-09-15
      use module_principal ; use module_parallele
      implicit none
      double precision tem_,sal_,prs_
      integer time_

! Ces lignes servent a verifier la disponibilite des tableaux generiques anyv3d  !30-09-15
      if(anyv3d(-1,-1,0,1)==-9999.)stop 'Err anyv3d1 not available'
      if(anyv3d(-1,-1,0,2)==-9999.)stop 'Err anyv3d2 not available'
      if(anyv3d(-1,-1,0,3)==-9999.)stop 'Err anyv3d3 not available'
      if(anyv3d(-1,-1,0,4)==-9999.)stop 'Err anyv3d4 not available'
      if(anyv3d(-1,-1,0,6)==-9999.)stop 'Err anyv3d5 not available'
      if(anyv3d(-1,-1,0,7)==-9999.)stop 'Err anyv3d6 not available'

      do j=0,jmax+1 ; do i=0,imax+1 !26-01-23

      x0=0.
      do k=kmax,1,-1

      x0=x0+h_w(i,j)*dsig_t(i,j,k)
      prs_=max(x0-0.5*dsig_t(i,j,k)*h_w(i,j),0.) !26-01-23

!     tem_=tem_t(i,j,k,time_)
!     sal_=max(sal_t(i,j,k,time_),zero)
! Appliquer A l'EOS des bornes de validites sur T et S            !07-02-20
      tem_=min(max(tem_t(i,j,k,time_),-3.),50.) !07-02-20
      sal_=min(max(sal_t(i,j,k,time_), 0.),50.) !07-02-20

! Compute density rho_    from (sal_   ,tem_   ,prs_   ) according
! to Jackett et al 2006 (JAOT):
!     tem2_= tem_*tem_
!     tem3_=tem2_*tem_
!     tem4_=tem3_*tem_
!     sal2_= sal_*sal_

      anyv3d(i,j,k,1)= c1_jmfwg                 &
                      +c2_jmfwg*tem_            &
                      +c3_jmfwg*tem_**2         &
                      +c4_jmfwg*tem_**3         &
                      +c5_jmfwg*sal_            &
                      +c6_jmfwg*sal_*tem_       &
                      +c7_jmfwg*sal_**2

      anyv3d(i,j,k,2)= c8_jmfwg                 &
                      +c9_jmfwg*tem_**2         &
                     +c10_jmfwg*sal_

      anyv3d(i,j,k,3)=c11_jmfwg                 &
                     +c12_jmfwg*tem_**2

      anyv3d(i,j,k,4)=1.                        &
                     +c14_jmfwg*tem_            &
                     +c15_jmfwg*tem_**2         &
                     +c16_jmfwg*tem_**3         &
                     +c17_jmfwg*tem_**4         &
                     +c18_jmfwg*sal_            &
                     +c19_jmfwg*tem_*sal_       &
                     +c20_jmfwg*sal_*tem_**3    &
                     +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem_**2)

!     anyv3d(i,j,k,5)=c23_jmfwg

      anyv3d(i,j,k,6)=c24_jmfwg*tem_**3

      anyv3d(i,j,k,7)=c25_jmfwg*tem_

! Potential density anomaly:
          rhp_t(i,j,k)=( anyv3d(i,j,k,1)             &
                        +anyv3d(i,j,k,2)*prs_        &
                        +anyv3d(i,j,k,3)*prs_**2 )   &
                      /( anyv3d(i,j,k,4)             &
                              +c23_jmfwg*prs_        &
                        +anyv3d(i,j,k,6)*prs_**2     &
                        +anyv3d(i,j,k,7)*prs_**3 )   &
                -rho

      enddo
      enddo
      enddo

! Marquer "indisponibles" les tableaux generiques !30-09-5
!     anyv3d(-1,-1,0,1:4)=-9999.  ; anyv3d(-1,-1,0,6:7)=-9999.

      i=imax/2 ; j=jmax/2 ; k=kmax/2 !16-02-17
      checkanyv3d(1)=anyv3d(i,j,k,1)
      checkanyv3d(2)=anyv3d(i,j,k,2)
      checkanyv3d(3)=anyv3d(i,j,k,3)
      checkanyv3d(4)=anyv3d(i,j,k,4)
      checkanyv3d(6)=anyv3d(i,j,k,6)
      checkanyv3d(7)=anyv3d(i,j,k,7)

      end subroutine equation_of_state_full_a1_a7

!................................................................................

      subroutine equation_of_state_full_anyv(istr_,iend_,jstr_,jend_)
      use module_principal
      use module_parallele
      implicit none
      double precision tem_,tem2_,tem3_,tem4_        &
                      ,sal_,sal2_,b0_,a0_,prs_,prs2_ &
                      ,prs3_
      integer istr_,iend_,jstr_,jend_

! JACKETT MCDOUGALL FEISTEL WRIGHT GRIFFIES, 2006
! Algorithms for Density, Potential Temperature, Conservative Temperature, and the
! Freezing Temperature of Seawater. JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY.
! A check value for this density equation is rho(35,25,2000)=1031.65056056576 kg m-3
! i.e. S=35 psu, Tpot=25°C, p=2000dbar. Other check values are
! rho(20,20,1000)=1017.72886801964 kg m-3 and rho(40,12,8000)=1062.95279820631 kg m-3


!     do k=1,kmax
      do j=jstr_,jend_
      do i=istr_,iend_

      x0=0.
      do k=kmax,1,-1
      x0=x0+h_w(i,j)*dsig_t(i,j,k)
!     prs_=-depth_t(i,j,k)
      prs_=x0-0.5*dsig_t(i,j,k)*h_w(i,j)
      prs2_= prs_*prs_
      prs3_=prs2_*prs_

      tem_=anyv3d(i,j,k,id_tem)
      tem2_= tem_*tem_
      tem3_=tem2_*tem_
      tem4_=tem3_*tem_

      sal_=max(anyv3d(i,j,k,id_sal),zero)
      sal2_=sal_*sal_


      b0_=                                                   &
           c1_jmfwg+c2_jmfwg*tem_+c3_jmfwg*tem2_             & ! b0
                   +c4_jmfwg*tem3_                           &
                   +c5_jmfwg*sal_                            &
                   +c6_jmfwg*sal_*tem_                       &
                   +c7_jmfwg*sal2_

      a0_=1.       +c14_jmfwg*tem_ +c15_jmfwg*tem2_          & ! /a0
                   +c16_jmfwg*tem3_+c17_jmfwg*tem4_          &
                   +c18_jmfwg*sal_                           &
                   +c19_jmfwg*tem_*sal_                      &
                   +c20_jmfwg*sal_*tem3_                     &
                   +(sal_**1.5)*(c21_jmfwg+c22_jmfwg*tem2_)


      anyv3d(i,j,k,id_rhf)=                                   &

      ( b0_                                                   &
      +(c8_jmfwg+c9_jmfwg*tem2_+c10_jmfwg*sal_)*prs_          &
      +(c11_jmfwg+c12_jmfwg*tem2_)*                prs2_  )   &

      /( a0_                                                  &
        +c23_jmfwg*         prs_                              &
        +c24_jmfwg*tem3_*prs2_                                &
        +c25_jmfwg*tem_ *prs3_   )                            &

       -rho                                                   &
       -rhcref_t(i,j,k)

      enddo
      enddo
      enddo


      end subroutine equation_of_state_full_anyv

