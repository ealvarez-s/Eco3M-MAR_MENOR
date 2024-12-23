      subroutine obc_river(case_,time1_,time2_)                                 !29/05/04
!______________________________________________________________________
! SYMPHONIE ocean model
! release 368 - last update: 07-04-23
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer case_,time1_,time2_
#ifdef synopsis
       subroutinetitle='obc_river'
       subroutinedescription= &
       'Update the fields at the river grid points'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         21/12/06: mis en service
!         21/04/07: passage � la coordonn�e curviligne
!         05-06-09  Parallelisation
! 2009.3  30-09-09  utilisation des nouveaux facteurs d'echelle verticale
!         11-10-09  tem_c et sal_c remplacent thz_c et shz_c
! 2010.7  22-02-10  velbar devient la variable d'etat
! 2010.8  24-03-10  Pour une conservation du flux de masse plus precise la
!                   CL considere flux_sum_t plutot que le debit de la riviere
!         03-05-10  nouveau schema forward
!         05-05-10  Profil vertical de vitesse donn� par C.L. gradient
! 2010.9  12-06-10  Flux advectif upwind � l'embouchure des fleuves
! 2010.13 01-11-10  river_inout remplace river_dom
! 2010.14 09-11-10  Possibilit� d'un contre-courant rentrant dans l'estuaire
!                   (voir obc 1) sous contraintes que le bilan de masse et
!                   de sel soient tous les 2 conserv�s
! 2010.20 19-04-11  Pas de temps modifiable
! S.26    19-04-14  river_inout change de nom
!         08-05-14  * C.L. T,S case1 rationnalis�.
!                   * C.L. u,v methode obc2 re-activee car obc3 conduit a trop de
!                     dillution
!         06-10-14  Modifs pour augmenter la robustesse de la conservation mpi
!                   des fleuves
!         09-01-15  nouvel algo obc 2 et suppression obc 1 et obc 3
!         22-01-15  C.L. vel_u vel_v harmonisee avec C.L. barotrope (mais en   
!                   encore trop simple)
!         20-07-15  adaptation de obc_river case3 a river_s non nul
!         14-07-16  seuil mini modifiE pour plus de precision sur le debit obtenu
!         16-03-17  river_obctype choix de la C.L. au point d'entree de la riviere
!         18-04-17  continuite mpi renforcee
!         13-06-17  modifier le volume du reservoir de la condition type 1 avec x10
! v251    09-04-19  subroutines obc_river_groundwater_w et obc_river_groundwater_s
! v252    17-04-19  ajout de la temperature dans la subroutine obc_river_groundwater_ts
! v296    18-02-21  ajout subroutine obc_river_surfriver_w !18-02-21
!                   ajout subroutine obc_river_surfriver_tem !18-02-21
! v297    05-03-21  rivieres surface et sous marine integrees dans le schema d'advection
! v301    05-04-21  Ajout d'une condition imposant la forme du profil vertical des 
!                   vitesses (pour la condition Bosphore de GLOBMED)
! v310    09-11-21  Avec approche VQS par limiteur de flux le courant au
!                   point source est nul pour k<kmerged (solution alternative plus simple
!                   que d'appliquer le limiteur sur la condition limite)
! v368    07-04-23  ajout subroutine obc_river_debug pour detecter contre courant 
!                   aux points des rivieres
!...............................................................................
!    _________                    .__                  .__             ! (�o�) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


!______________________________________________________________________________
! CONDITIONS LIMITES POUR LES RIVIERES:
!______________________________________________________________________________


! Si 0 riviere on resort
      if(nriver==0)return


!***************************************************************************
! CONDITIONS SUR VITESSES BAROCLINES:
! DEBUT
!***************************************************************************

      if(case_==2) then !-case2-case2-case2-case2->

!     const2=0.5/real(k2dfin)
      const2=1./(iteration2d_max_now+iteration2d_max_bef)   !19-04-11

      do 250 kr=1,nriver
      if(riverdir(kr)>0) then !m�v�m> !09-04-19

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
!     if(river_dom(kr) == par%rank) then !rrrrrrrrrrrrrrrrrr>
      if(rivervel_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>          !01-11-10
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Condition radiative sur les vitesses. La vitesse est ensuite ajust�e pour
! garantir la conservation du volume

      i=iriver(kr,2)
      j=jriver(kr,2)

! AIGUILLAGE TYPE DE C.L.
! OBCTYPE=0 : Pas de contre-courant
! OBCTYPE=1 : contre-courant possible (Sr tel que somme(U.Sr)=0)
! OBCTYPE=2 : Un profil dont la forme est imposEe

      if(river_obctype(kr)==0) then !000000000000> AIGUILLAGE TYPE DE C.L. Pas de contre-courant !16-03-17

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->

      j1=1
      if(riverdir(kr)==4)j1=-1

! Pratiquee au temps 2 pour conservation mpi: first guess !06-10-14
!  gradient nul sur vel(2) si sens riviere, 0 si contre-sens
      if(time2_==after) then !222222>
       stop 'remove call obc_river(2,2,2)in internal_mode.F90' !18-04-17
!      do k=kmin_v(i,j),kmax 
!       vel_v(i,j,k,after)=                              &
!       0.25*(       (vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now)) &  ! obc 2
!             +j1*abs(vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now))  ) !09-01-15
!      enddo             
      endif                  !222222>

! Equivalent du couple mode pour les fleuves - Correction facteur multiplicatif !06-10-14
! pour conserver le profil vertical du first guess - Appele pres couple_mode
      if(time2_==now) then !111111>

       sum1=0.
!      do k=kmin_v(i,j),kmax
       anyv1d(:,1)=0. !profil nul sous kmerged !09-11-21
       do k=kmerged_v(i,j),kmax !09-11-21
! Comme couple_mode peut changer le signe de vel on retablit ce dernier
        anyv1d(k,1)=                                             &
        0.25*(       (vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now)) &  ! obc 2
              +j1*abs(vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now))  ) !09-01-15

        anyv1d(k,1)=j1*max(0.001,j1*anyv1d(k,1))  !14-07-16
        sum1=sum1+anyv1d(k,1)*dz_v(i,j,k,now)*dx_v(i,j)
       enddo
! obc 2 ou 3:
       sum1=sign(dmax1(1.d-10,abs(sum1)),sum1)
       x1=const2*(fluxbar_sumt_v(i,j,0)+fluxbar_sumt_v(i,j,1))/sum1  !24-03-10
       do k=kmin_v(i,j),kmax
        veldxdz_v(i,j,k,now)=x1*anyv1d(k,1)*dx_v(i,j)*dz_v(i,j,k,now)
!       vel_v(i,j,k,now)=0. ! Pour le momnt on s'harmonise avec la C.L. barotrope !22-01-15
       enddo 

      endif              !111111>


      else !----------------------------------------------->

      i1=1
      if(riverdir(kr)==3)i1=-1

! Pratiquee au temps 2 pour conservation mpi: first guess !06-10-14
!  gradient nul sur vel(2) si sens riviere, 0 si contre-sens
      if(time2_==after) then !222222>
       stop 'remove call obc_river(2,2,2)in internal_mode.F90' !18-04-17
!      do k=kmin_u(i,j),kmax
!       vel_u(i,j,k,after)=       &
!       0.25*(       (vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now)) &  !  obc 2
!             +i1*abs(vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now))  ) !09-01-15
!      enddo 
      endif                 !222222>
    
! Equivalent du couple mode pour les fleuves - Correction facteur multiplicatif !06-10-14
! pour conserver le profil vertical du first guess - Appele pres couple_mode
      if(time2_==now) then !111111>

       sum1=0.
!      do k=kmin_u(i,j),kmax
       anyv1d(:,1)=0. !profil nul sous kmerged !09-11-21
       do k=kmerged_u(i,j),kmax !09-11-21
! Comme couple_mode peut changer le signe de vel on retablit ce dernier
        anyv1d(k,1)=       &
        0.25*(       (vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now)) &  !  obc 2
              +i1*abs(vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now))  ) !09-01-15
        anyv1d(k,1)=i1*max(0.001,i1*anyv1d(k,1))  !14-07-16
        sum1=sum1+anyv1d(k,1)*dz_u(i,j,k,now)*dy_u(i,j)
       enddo 
! Obc 2 ou 3:
      sum1=sign(dmax1(1.d-10,abs(sum1)),sum1)
      x1=const2*(fluxbar_sumt_u(i,j,0)+fluxbar_sumt_u(i,j,1))/sum1 !24-03-10
       do k=kmin_u(i,j),kmax
        veldydz_u(i,j,k,now)=x1*anyv1d(k,1)*dy_u(i,j)*dz_u(i,j,k,now)
!       vel_u(i,j,k,now)=0. ! Pour le moment on s'harmonise avec la C.L. barotrope !22-01-15
       enddo 

      endif              !111111>

      endif !---------------------------------------------->

      endif                         !000000000000> ! Aiguillage type de C.L.


      if(river_obctype(kr)==1) then !111111111111> AIGUILLAGE TYPE DE C.L. Possible contre-courant !16-03-17

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->

      j1=1
      if(riverdir(kr)==4)j1=-1

! Pratiquee au temps 2 pour conservation mpi: first guess !06-10-14
!  gradient nul sur vel(2) si sens riviere, 0 si contre-sens
      if(time2_==after) then !222222>
       stop 'remove call obc_river(2,2,2)in internal_mode.F90' !18-04-17
!      do k=kmin_v(i,j),kmax 
!       vel_v(i,j,k,after)=                                            &
!                 0.5*(vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now))      &
!        +sign(1.d0,   vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now) )*0.001
!      enddo             
      endif                  !222222>

! Equivalent du couple mode pour les fleuves - Correction facteur multiplicatif !06-10-14
! pour conserver le profil vertical du first guess - Appele pres couple_mode
      if(time2_==now) then !111111>

! v'=v+x1*|v| avec somme(v')=debit ==> ! x1=(debit-somme(v))/somme(|v| )
       sum1=0.
       sum2=small2
!      do k=kmin_v(i,j),kmax
       anyv1d(:,1)=0. !profil nul sous kmerged !09-11-21
       do k=kmerged_v(i,j),kmax !09-11-21

        anyv1d(k,1)=                                                   &
                  0.5*(vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now))      &
         +sign(1.d0,   vel_v(i,j+j1,k,after)+vel_v(i,j+j1,k,now) )*0.001


        sum1=sum1+    anyv1d(k,1) *dz_v(i,j,k,now)*dx_v(i,j) !somme(v)
        sum2=sum2+abs(anyv1d(k,1))*dz_v(i,j,k,now)*dx_v(i,j) !somme(|v|)
       enddo

! x1=(debit-somme(v))/somme(sqrt(v))
       x1=(const2*(fluxbar_sumt_v(i,j,0)+fluxbar_sumt_v(i,j,1))-sum1)/sum2

       do k=kmin_v(i,j),kmax
        veldxdz_v(i,j,k,now)=                          &
           (anyv1d(k,1)+x1*abs(anyv1d(k,1))) & !v'=v+x1*|v|
            *dz_v(i,j,k,now)                           &
            *dx_v(i,j)
!           vel_v(i,j,k,now)=0. ! Pour le moment on s'harmonise avec la C.L. barotrope !22-01-15
       enddo 

      endif              !111111>


      else !----------------------------------------------->

      i1=1
      if(riverdir(kr)==3)i1=-1

! Pratiquee au temps 2 pour conservation mpi: first guess !06-10-14
!  gradient nul sur vel(2) si sens riviere, 0 si contre-sens
      if(time2_==after) then !222222>
       stop 'remove call obc_river(2,2,2)in internal_mode.F90' !18-04-17
!      do k=kmin_u(i,j),kmax
!       vel_u(i,j,k,after)=                                    &
!            0.5*(vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now))   &
!      +sign(1.d0,vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now))*0.001
!      enddo 
      endif                 !222222>
    
! Equivalent du couple mode pour les fleuves - Correction facteur multiplicatif !06-10-14
! pour conserver le profil vertical du first guess - Appele pres couple_mode
      if(time2_==now) then !111111>

! v'=v+x1*|v| avec somme(v')=debit ==> ! x1=(debit-somme(v))/somme(|v| )
       sum1=0.
       sum2=small2
!      do k=kmin_u(i,j),kmax
       anyv1d(:,1)=0. !profil nul sous kmerged !09-11-21
       do k=kmerged_u(i,j),kmax !09-11-21

        anyv1d(k,1)=                                           &
             0.5*(vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now))   &
       +sign(1.d0,vel_u(i+i1,j,k,after)+vel_u(i+i1,j,k,now))*0.001

        sum1=sum1+    anyv1d(k,1) *dz_u(i,j,k,now)*dy_u(i,j)
        sum2=sum2+abs(anyv1d(k,1))*dz_u(i,j,k,now)*dy_u(i,j)
       enddo 

      x1=(const2*(fluxbar_sumt_u(i,j,0)+fluxbar_sumt_u(i,j,1))-sum1)/sum2

       do k=kmin_u(i,j),kmax
        veldydz_u(i,j,k,now)=                          &
           (anyv1d(k,1)+x1*abs(anyv1d(k,1))) & !v'=v+x1*|v|
            *dz_u(i,j,k,now)                           &
            *dy_u(i,j)
!           vel_u(i,j,k,now)=0. ! Pour le moment on s'harmonise avec la C.L. barotrope !22-01-15
       enddo 

      endif              !111111>

      endif !---------------------------------------------->

      endif                         !111111111111> ! Aiguillage type de C.L.


      if(river_obctype(kr)==2) then !2222222> AIGUILLAGE TYPE DE C.L.: PROFIL VERTICAL IMPOSE !05-04-21

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->

! Pour placer le point vitesse aval:
!     j1=1 ; if(riverdir(kr)==4)j1=-1
! Pour placer le point traceur aval:
      j1=0 ; if(riverdir(kr)==4)j1=-1

      if(time2_==now) then !111111>

       sum1=0.
!      do k=kmin_v(i,j),kmax
       anyv1d(:,1)=0. !profil nul sous kmerged !09-11-21
       do k=kmerged_v(i,j),kmax !09-11-21
        anyv1d(k,1)= &
        sigma_w(i,j+j1,k)+sigma_w(i,j+j1,k+1)-0.82 ! Un profil A 2 sens avec transport de fond = (-1)*transport net
        sum1=sum1+anyv1d(k,1)*dz_v(i,j,k,now)*dx_v(i,j)
       enddo
       x1=const2*(fluxbar_sumt_v(i,j,0)+fluxbar_sumt_v(i,j,1))/sum1  !24-03-10
       do k=kmin_v(i,j),kmax
        veldxdz_v(i,j,k,now)=x1*anyv1d(k,1)*dx_v(i,j)*dz_v(i,j,k,now)
       enddo 

      endif              !111111>


      else !----------------------------------------------->

! Pour placer le point vitesse aval:
!     i1=1 ; if(riverdir(kr)==3)i1=-1
! Pour placer le point traceur aval:
      i1=0 ; if(riverdir(kr)==3)i1=-1

      if(time2_==now) then !111111>

       sum1=0.
!      do k=kmin_u(i,j),kmax
       anyv1d(:,1)=0. !profil nul sous kmerged !09-11-21
       do k=kmerged_u(i,j),kmax !09-11-21
        anyv1d(k,1)=       &
        sigma_w(i+i1,j,k)+sigma_w(i+i1,j,k+1)-0.82 ! Un profil A 2 sens avec transport de fond = (-1)*transport net
        sum1=sum1+anyv1d(k,1)*dz_u(i,j,k,now)*dy_u(i,j)
       enddo 
       x1=const2*(fluxbar_sumt_u(i,j,0)+fluxbar_sumt_u(i,j,1))/sum1 !24-03-10
       do k=kmin_u(i,j),kmax
        veldydz_u(i,j,k,now)=x1*anyv1d(k,1)*dy_u(i,j)*dz_u(i,j,k,now)
       enddo 

      endif              !111111>

      endif !---------------------------------------------->

      endif                         !2222222> AIGUILLAGE TYPE DE C.L.: PROFIL VERTICAL IMPOSE 

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
      endif            !rrrrrrrrrrrrrrrrrr>
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif    ! test sur riverdir(kr)>0                 !m�v�m> !09-04-19
  250 continue ! fin de boucle sur KR

      endif             !-case2-case2-case2-case2->


!***************************************************************************
! CONDITIONS SUR VITESSES BAROCLINES:
! FIN
!***************************************************************************



!***************************************************************************
! CONDITIONS SUR T ET S:
! DEBUT
!***************************************************************************


      if(case_==1) then !-case1-case1-case1->

      do kr=1,nriver ! boucle sur KR d�but.
      if(riverdir(kr)>0) then !m�v�m> !09-04-19

      if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>          !01-11-10

      i=iriver(kr,1)
      j=jriver(kr,1)

       do k=1,kmax    ! boucle sur K d�but.
        tem_t(i,j,k,time1_:time2_)=river_t(kr,1)                       !11-10-09
        sal_t(i,j,k,time1_:time2_)=river_s(kr)
       enddo          ! boucle sur K fin.

      endif                              !rrrrrrrrrrrrrrrrrr>


! Si condition permettant contre-courant, alors river_s ne peut plus etre
! simplement 0 mais doit evoluer pour conserver le bilan de sel.
! Plus de details dans: !16-03-17
! https://docs.google.com/document/d/1dy5epHqM-nSF3xQsGFX2eS439aqu3Gv6t-Buw5YLT88/edit

#ifdef bidon
      if(river_obctype(kr)==1) then !obctype-obctype-> !16-03-17


      if(rivervel_inout(kr)==1) then !rivervel_inout>
       i=iriver(kr,2)
       j=jriver(kr,2)
       sum1=0.

       x10=1./(200.*dxdy_t(i,j)*5.) ! 1./(Volume of the reservoir) !13-06-17

       if(riverdir(kr)==1) then !1111>
!       sum0=0.
        do k=1,kmax
         sum1=sum1+0.5*( &
          (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*river_s(kr)    &
         +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*sal_t(i,j,k,0)) 
!        sum0=sum0+veldydz_u(i,j,k,1)
        enddo
!       river_s(kr)=river_s(kr)-dti_fw*sum1/(dxdy_t(i-1,j)*10.*100.) ! volume boite: dxdy*10m*N
        river_s(kr)=river_s(kr)-dti_fw*sum1*x10
!       if(i-1/=iriver(kr,1))stop 'PAS COMPRIS 1'
!       write(6,*)'sum0 sum1=',real(sum0),real(sum1),river_s(kr)
       endif                    !1111>
       if(riverdir(kr)==2) then !2222>
!       sum0=0.
        do k=1,kmax
         sum1=sum1+0.5*( &
          (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*river_s(kr)    &
         +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*sal_t(i,j,k,0)) 
!       sum0=sum0+veldxdz_v(i,j,k,1)
        enddo
!       river_s(kr)=river_s(kr)-dti_fw*sum1/(dxdy_t(i,j-1)*10.) ! volume boite: dxdy*10m
        river_s(kr)=river_s(kr)-dti_fw*sum1*x10
!       if(j-1/=jriver(kr,1))stop 'PAS COMPRIS 2'
!       write(6,*)'sum0 sum1=',real(sum0),real(sum1),river_s(kr)
       endif                    !2222> 
       if(riverdir(kr)==3) then !3333>
!       sum0=0.
        do k=1,kmax
         sum1=sum1+0.5*( &
          (veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*river_s(kr)    &
         +(veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*sal_t(i-1,j,k,0)) 
        sum0=sum0+veldydz_u(i,j,k,1)
        enddo
        river_s(kr)=river_s(kr)+dti_fw*sum1*x10
!       if(i/=iriver(kr,1))stop 'PAS COMPRIS 3'
!       write(6,*)'sum0 sum1=',real(sum0),real(sum1),river_s(kr)
       endif                    !3333>
       if(riverdir(kr)==4) then !4444>
!       sum0=0.
        do k=1,kmax
         sum1=sum1+0.5*( &
          (veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*river_s(kr)    &
         +(veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*sal_t(i,j-1,k,0)) 
!       sum0=sum0+veldxdz_v(i,j,k,1)
        enddo
        river_s(kr)=river_s(kr)+dti_fw*sum1*x10
!       write(6,*)'sum0 sum1=',real(sum0),real(sum1),river_s(kr)
!       if(j/=jriver(kr,1))stop 'PAS COMPRIS 4'
       endif                    !4444>


!      write(10+par%rank,*)elapsedtime_now/86400.,river_s(kr)

      endif                          !rivervel_inout>

      endif                         !obctype-obctype->
#endif

      endif                   !m�v�m> test riverdir(kr)>0 !09-04-19 
      enddo          ! boucle sur KR fin.
      return

      endif             !-case1-case1-case1->

!***************************************************************************
! CONDITIONS SUR T ET S:
! FIN
!***************************************************************************


!....................................                                 !12-06-10
! BEGIN T & S fluxes at river mounth
      if(case_==3) then !-case3-case3-case3->
      stop 'Passe pas par case_==3'
#ifdef bidon
!....................................

      do kr=1,nriver

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
!     if(river_dom(kr) == par%rank) then !rrrrrrrrrrrrrrrrrr>
!     if(rivertrc_inout(kr)==1)        then !rrrrrrrrrrrrrrrrrr>          !01-11-10
      if(rivervel_inout(kr)==1)        then !rrrrrrrrrrrrrrrrrr>          !08-05-14
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Note: dans le cas d'une condition radiative de courant permettant
! un contre-courant entrant dans l'estuaire, la condition aux limites
! sur S permet de conserver le bilan de sel, � savoir Somme(vel.S)=0.
! Si aucun contre-courant ne penetre dans l'estuaire, cette condition
! conduit � S=0
! Le principe general est le suivant: Si un contre-courant entrant existe, la conservation
! du volume conduit a ce que le debit sortant est plus grand que le debit du fleuve.
! Debit(entrant)+Debit(sortant)=Debit(Fleuve) dans obc case2 de cette meme routine.
! Imposer S=0 menerait a une exageration de la dessalure. La valeur de S est donc
! augmentee de telle sorte que le volume d'eau douce sortant soit bien celui du fleuve.
! Somme( Usortant.Ssortant ) + Somme (Uentrant.Socean) = 0
! => Ssortant=-Somme(Uentrant.Socean)/Somme( Usortant )
! Note1: hypothese d'homogeneisation qui conduit a considerer Ssortant constant sur la verticale
! Note2: Si aucun courant entrant alors Ssortant=0
! Note3: le principe est de corriger une valeur a priori de Ssortant. Jusqu'� maintenant l'appel
!        au case1 imposait comme valeur de depart 0. Mais on peut reflechir � partir d'un
!        gradient nul dont l'objectif serait d'ameliorer le comportement de la partie diffusive
!        du QUICK

      i=iriver(kr,2)
      j=jriver(kr,2)

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->

       j0=1
       if(riverdir(kr)==4)j0=-1

       sum1=0.
       sum2=small2
       do k=1,kmax
        anyv1d(k,1)=max(zero,j0*sign(un,veldxdz_v(i,j,k,1)))  ! =1 si courant sortant, =0 contre-courant ocean
        j1=j-( nint(sign(un,veldxdz_v(i,j,k,1))) + 1 )/2
        sum1=sum1+veldxdz_v(i,j,k,1)*(sal_t(i,j1,k,0)-river_s(kr)) !20-07-15
        sum2=sum2+veldxdz_v(i,j,k,1)*anyv1d(k,1)
       enddo
!      write(*,*)'---------------------------------'
!      write(*,*)'Verif river avant',sum1
        x1=-sum1/sum2
!      write(*,*)'DS=',x1
!      sum1=0.
       do k=1,kmax ! attention k demarre � 1 et non pas kmin_v
                   ! car les tableaux sont integralement consid�r�s dans
                   ! l'etape d'advection � suivre...

        j1=j-( nint(sign(un,veldxdz_v(i,j,k,1))) + 1 )/2
        sal_t(i,j1,k,0)=sal_t(i,j1,k,0)+x1*anyv1d(k,1) ! anyv1d ne modifie que Ssortant et pas Socean !
! sal(1)=sal(0) si point source, sal(1) inchange si point ocean. Pourquoi? Parce qu'
!       sal_t(i,j1,k,1)=anyv1d(k,1) *sal_t(i,j1,k,0)  &
!                  +(1.-anyv1d(k,1))*sal_t(i,j1,k,1)
! a l'iteration suivante, sal(1) sera sal(0), qui permettra au schema d'advection d'avoir une valeur
!  dans l'embouchure pour le calcul de l'hyperdiffusivite. Note qu'a priori le flux advectif sur 2 points
! sal(1) n'a pas besoin du point source.

        anyv3d(i,j,k,2)=-dti_lp*veldxdz_v(i,j,k,1)*tem_t(i,j1,k,0)
        anyv3d(i,j,k,5)=-dti_lp*veldxdz_v(i,j,k,1)*sal_t(i,j1,k,0)
!       sum1=sum1+anyv3d(i,j,k,5)
       enddo
!      write(*,*)'Verif river apres',-sum1/dti_lp

      else !----------------------------------------------->

       i0=1
       if(riverdir(kr)==3)i0=-1

       sum1=0.
       sum2=small2
       do k=1,kmax
        anyv1d(k,1)=max(zero,i0*sign(un,veldydz_u(i,j,k,1)))
        i1=i-( nint(sign(un,veldydz_u(i,j,k,1))) + 1 )/2
        sum1=sum1+veldydz_u(i,j,k,1)*(sal_t(i1,j,k,0)-river_s(kr)) !20-07-15
        sum2=sum2+veldydz_u(i,j,k,1)*anyv1d(k,1)
       enddo

        x1=-sum1/sum2
!      write(*,*)'---------------------------------'
!      write(*,*)'Verif river avant',sum1
!      write(*,*)'DS=',x1
!      sum1=0.

       do k=1,kmax ! attention k demarre � 1 et pas kmin_u
        i1=i-( nint(sign(un,veldydz_u(i,j,k,1))) + 1 )/2
        sal_t(i1,j,k,0)=sal_t(i1,j,k,0)+x1*anyv1d(k,1)
! sal(1)=sal(0) si point source, sal(1) inchange si point ocean. Pourquoi? Parce qu'
!       sal_t(i1,j,k,1)=anyv1d(k,1) *sal_t(i1,j,k,0)  &
!                  +(1.-anyv1d(k,1))*sal_t(i1,j,k,1)
! a l'iteration suivante, sal(1) sera sal(0), qui permettra au schema d'advection d'avoir une valeur
!  dans l'embouchure pour le calcul de l'hyperdiffusivite. Note qu'a priori le flux advectif sur 2 points
! sal(1) n'a pas besoin du point source.
        anyv3d(i,j,k,1)=-dti_lp*veldydz_u(i,j,k,1)*tem_t(i1,j,k,0)
        anyv3d(i,j,k,4)=-dti_lp*veldydz_u(i,j,k,1)*sal_t(i1,j,k,0)
!       sum1=sum1+anyv3d(i,j,k,4)
       enddo
!      write(*,*)'Verif river apres',-sum1/dti_lp

      endif !---------------------------------------------->

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
      endif            !rrrrrrrrrrrrrrrrrr>
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo ! fin de boucle sur KR

!....................................
! END of T & S fluxes at river mounth
      return
#endif
      endif           !-case3-case3-case3->
!....................................

      end subroutine obc_river

!........................................................................

      subroutine obc_river_groundwater_w !09-04-19
      use module_principal ; use module_parallele 
      implicit none

! Underwater sources. 
! Determination of vertical velocity from a "vertical" flux of fresh water
      do kr=1,nriver
       if(riverdir(kr)==-1) then !m�v�m> !09-04-19
        if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>  

         i=iriver(kr,1)
         j=jriver(kr,1)
         omega_w(i,j,1,1)=riverflux(kr,1)/dxdy_t(i,j)

        endif                              !rrrrrrrrrrrrrrrrrr>
       endif                     !m�v�m> !09-04-19
      enddo ! kr loop

      end subroutine obc_river_groundwater_w

!........................................................................
! supprime le 05-03-21 car rivieres surface et sous marine integrees dans le schema d'advection
#ifdef bidon
      subroutine obc_river_groundwater_ts !09-04-19
      use module_principal ; use module_parallele
      implicit none

! Underwater sources. 
! Update salinity (consistency with vertical water flux at bottom)
! Principe: l'etape precedente a fait l'hypothese d'un flux advectif vertical basE sur une
! hypothese de gradient nul dS/dz=0 (contingence temporaire liee A la matrice "figee" du systeme implicite).
! Le vrai flux etant nul (puisque salinite de la source=0) on retablit maintenant le vrai flux en enlevant
! le flux indument appliquE A l'etape precedente.
! Note: dans le cas des coord. fusionnees, l'etape de melange suivante redistribura la mise A jour
! sur les couches de k=1 A k=kmerged_t
! Cas de la temperature: quasi identique A celui de la salinite sauf que la temperature
! de la source sous marine n'est pas nulle et par consequent en plus d'enlever le flux de gradient
! nul on reinjecte un flux de temperature donnEe par river_t(kr,1)
      do kr=1,nriver
       if(riverdir(kr)==-1) then !m�v�m> !09-04-19
        if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>  

         i=iriver(kr,1)
         j=jriver(kr,1)
         k=kmerged_t(i,j)
! note: riverflux(kr,1)/dxdy_t(i,j) est omega_w(i,j,1,1) dans obc_river_groundwater_w

                 sal_t(i,j,k,2)=  &
           ( & !ooo>
                 sal_t(i,j,k,2)   &
         -dti_lp*sal_t(i,j,1,2)*riverflux(kr,1)/(dxdy_t(i,j)*dz_t(i,j,k,2)) &
           ) & !ooo>
               *mask_t(i,j,kmax)     &       
            *wetmask_t(i,j)          &   
        +(1.-wetmask_t(i,j))         &
                *sal_t(i,j,k,1) 

                 tem_t(i,j,k,2)=  & !17-04-19
           ( & !ooo>
                 tem_t(i,j,k,2)   &
        -dti_lp*(tem_t(i,j,1,2)-river_t(kr,1))*riverflux(kr,1)/(dxdy_t(i,j)*dz_t(i,j,k,2)) &
           ) & !ooo>
               *mask_t(i,j,kmax)     &       
            *wetmask_t(i,j)          &   
        +(1.-wetmask_t(i,j))         &
                *tem_t(i,j,k,1) 

        endif                              !rrrrrrrrrrrrrrrrrr>
       endif                     !m�v�m> !09-04-19
      enddo ! kr loop

      end subroutine obc_river_groundwater_ts
#endif
!........................................................................

      subroutine obc_river_surfriver_w !18-02-21
      use module_principal ; use module_parallele 
      implicit none

! ATTENTION ON N'A PAS CHANGE LE SIGNE DU DEBIT LU DANS LES FICHIERS (qui est donc positif) 
! ON CONSIDERE DONC -ABS(riverflux(kr,1)) POUR QUE LE FLUX SOIT DIRIGE VERS LE BAS

! Surface rivers
! Determination of vertical velocity from a "vertical" flux of fresh water
      do kr=1,nriver
       if(riverdir(kr)==0) then !m�v�m> 
        if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>  

         i=iriver(kr,1)
         j=jriver(kr,1)
         omega_w(i,j,kmax+1,1)=omega_w(i,j,kmax+1,1)*flag_omega_cumul   &
                              -abs(riverflux(kr,1))/dxdy_t(i,j)

         omega_evaprec_w(i,j,1)=omega_w(i,j,kmaxp1,1)

        endif                              !rrrrrrrrrrrrrrrrrr>
       endif                    !m�v�m> 
      enddo ! kr loop

! Ce flag indique que omega_w(i,j,kmaxp1,1) est desormais non nul et que 
! les autres operations sur omega de surface doivent incrementer omega_w 
! sans perdre la valeur acquise dans la presente subroutine
      flag_omega_cumul=1

      end subroutine obc_river_surfriver_w

!........................................................................
! supprime le 05-03-21 car rivieres surface et sous marine integrees dans le schema d'advection
#ifdef bidon
      subroutine obc_river_surfriver_tem !18-02-21
      use module_principal ; use module_parallele
      implicit none

! ATTENTION ON N'A PAS CHANGE LE SIGNE DU DEBIT LU DANS LES FICHIERS (qui est donc positif) 
! ON CONSIDERE DONC -ABS(riverflux(kr,1)) POUR QUE LE FLUX SOIT DIRIGE VERS LE BAS

! surface river sources. 
! L'etape precedente suppose gradient nul pour T, comme pour S. Cette condition
! indument appliquee est corrigee pour S (puisque la vraie condition est flux nulle)
! dans vertmix, mais pas pour T qui ne peut se satifaire d'une condition de gradient
! nul si on a specifie une temperature de la riviere). On corrige donc en tenant
! compte de la vraie temperature de la riviere.
      do kr=1,nriver
       if(riverdir(kr)==0) then !m�v�m> 
        if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>  

         i=iriver(kr,1)
         j=jriver(kr,1)
         k=kmax

! note: riverflux(kr,1)/dxdy_t(i,j) est omega_w(i,j,kmax+1,1) dans obc_river_surfriver_w
                 tem_t(i,j,k,2)=  & 
           ( & !ooo>
                 tem_t(i,j,k,2)   &
        -dti_lp*(tem_t(i,j,1,2)-river_t(kr,1))*abs(riverflux(kr,1))/(dxdy_t(i,j)*dz_t(i,j,k,2)) &
           ) & !ooo>
               *mask_t(i,j,kmax)     &       
            *wetmask_t(i,j)          &   
        +(1.-wetmask_t(i,j))         &
                *tem_t(i,j,k,1) 

        endif                              !rrrrrrrrrrrrrrrrrr>
       endif                   !m�v�m> 
      enddo ! kr loop

      end subroutine obc_river_surfriver_tem
#endif

!........................................................................

      subroutine obc_river_debug(jalon_) !07-04-23
      use module_principal ; use module_parallele 
      implicit none
      integer jalon_

! ICI ON VERIFIE QUE LA CONDITION 0 EST BIEN RESPECTEE (c.a.d. pas de
! contre courant)

      do kr=1,nriver
      if(riverdir(kr)>0) then !m�v�m> !09-04-19

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
!     if(river_dom(kr) == par%rank) then !rrrrrrrrrrrrrrrrrr>
      if(rivervel_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>          !01-11-10
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Condition radiative sur les vitesses. La vitesse est ensuite ajust�e pour
! garantir la conservation du volume

      i=iriver(kr,2)
      j=jriver(kr,2)

! OBCTYPE=0 : Pas de contre-courant
! OBCTYPE=1 : contre-courant possible (Sr tel que somme(U.Sr)=0)
! OBCTYPE=2 : Un profil dont la forme est imposEe

      if(river_obctype(kr)==0) then !000000000000> AIGUILLAGE TYPE DE C.L. Pas de contre-courant !16-03-17

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->

      j1=1
      if(riverdir(kr)==4)j1=-1

! Equivalent du couple mode pour les fleuves - Correction facteur multiplicatif !06-10-14
! pour conserver le profil vertical du first guess - Appele pres couple_mode

       if(j1==1.and.minval(veldxdz_v(i,j,:,:))<0.)  then !oooo>
        write(10+par%rank,*)'contrecourant NEGATIF dans riviere n�',kr
        write(10+par%rank,*)'i,j,loc _v',i,j
        write(10+par%rank,*)'i,j,glb _v',i+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'direction ',riverdir(kr)
        write(10+par%rank,*)'jalon_    ',jalon_
        do k=kmax,1,-1
         write(10+par%rank,*)'k,veldxdz_v',k,veldxdz_v(i,j,k,:)
        enddo
        flag_stop=1 
       endif                                             !oooo>
       if(j1==-1.and.maxval(veldxdz_v(i,j,:,:))>0.) then !oooo>
        write(10+par%rank,*)'contrecourant POSITIF dans riviere n�',kr
        write(10+par%rank,*)'i,j,loc _v',i,j
        write(10+par%rank,*)'i,j,glb _v',i+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'direction ',riverdir(kr)
        write(10+par%rank,*)'jalon_    ',jalon_
        do k=kmax,1,-1
         write(10+par%rank,*)'k,veldxdz_v',k,veldxdz_v(i,j,k,:)
        enddo
        flag_stop=1 
       endif                                            !oooo>

      else !----------------------------------------------->

      i1=1
      if(riverdir(kr)==3)i1=-1

       if(i1==1.and.minval(veldydz_u(i,j,:,:))<0.)  then !oooo>
        write(10+par%rank,*)'contrecourant NEGATIF dans riviere n�',kr
        write(10+par%rank,*)'i,j,loc _u',i,j
        write(10+par%rank,*)'i,j,glb _u',i+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'direction ',riverdir(kr)
        write(10+par%rank,*)'jalon_    ',jalon_
        do k=kmax,1,-1
         write(10+par%rank,*)'k,veldydz_u',k,veldydz_u(i,j,k,:)
        enddo
        flag_stop=1 
       endif                                             !oooo>
       if(i1==-1.and.maxval(veldydz_u(i,j,:,:))>0.) then !oooo>
        write(10+par%rank,*)'contrecourant POSITIF dans riviere n�',kr
        write(10+par%rank,*)'i,j,loc _u',i,j
        write(10+par%rank,*)'i,j,glb _u',i+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'direction ',riverdir(kr)
        write(10+par%rank,*)'jalon_    ',jalon_
        do k=kmax,1,-1
         write(10+par%rank,*)'k,veldxdz_v',k,veldydz_u(i,j,k,:)
        enddo
        flag_stop=1 
       endif                                            !oooo>

      endif !---------------------------------------------->

      endif                         !000000000000> ! Aiguillage type de C.L.


      if(river_obctype(kr)==1) then !111111111111> AIGUILLAGE TYPE DE C.L. Possible contre-courant !16-03-17
      endif                         !111111111111> ! Aiguillage type de C.L.


      if(river_obctype(kr)==2) then !2222222> AIGUILLAGE TYPE DE C.L.: PROFIL VERTICAL IMPOSE !05-04-21
      endif                         !2222222> AIGUILLAGE TYPE DE C.L.: PROFIL VERTICAL IMPOSE 

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
      endif            !rrrrrrrrrrrrrrrrrr>
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif    ! test sur riverdir(kr)>0                 !m�v�m> !09-04-19
      enddo ! fin de boucle sur KR

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop 'Err counter current in river. See fort.xxx files'

      end subroutine obc_river_debug !07-04-23
