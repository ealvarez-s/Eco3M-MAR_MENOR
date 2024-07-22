










      subroutine cellbox_thickness(case_,tstr_,tend_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 285 - last update: 06-06-20
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
      integer case_,tstr_,tend_,t_

!$ Computes the thickness of the cell boxes according to dz=(h+ssh)*dsig

!......................................................................
! Version date      Description des modifs
! 2009.3  30-09-09  Mise en service
!         08-10-09  - Passage au lagrangien
!                   - Les tableaux dz_ doient être déclarés de 0 à nr
!                     et initialisés en consequence
!         19-10-09  tstr_ & tend_ passés en argument
!         20-10-09  Mieux calculer l'etat initial pour ameliorer la conservation
!                   de l'etat initiale
!         15-12-09  Formalisation de l'approche hybride avec 2 choix possibles
!                   (nouvelle approche lagrangienne ou ancienne approche sigma)
!                   selon la valeur de la constante de temps de rappel vers
!                   le maillage de reference
!         16-12-09  suite du point precedent
! 2010.6  13-02-10  Filtre type asselin modifié sur dz
! 2010.8  10-03-10  beforebefore remplacé par before2
!         03-05-10  filtre d'asselin complété par un "moveforward" sauvant dz(t-2)
! 2010.9  28-05-10  debug test sur mask
!         13-06-10  ajout de ale_selected
! 2010.12 26-09-10  parametrage affiné
! 2010.13 11-10-10  ajout cellboxfactor1 & cellboxfactor2 définis dans set_parameters.F90
! 2010.15 27-12-10  dz en k=0 et kmax+1 = small3
! 2010.16 12-01-11  Filtre temporel ordre élevé
! 2010.18 26-03-11  small1 remplace small3
! 2010.23 24-05-11  initialiser dz_t(i,j,k,-1)
! S.26    24-03-14  ajout flag_nemoffline
!         02-09-14  prevoir le cas de dimensions non connues a l'avance pour
!                   dz_u et dz_v. Boucles separees pour dz_u et dz_v
!         28-01-16  filiere ale obsolete
! v285    05-06-20  suppression dsig_u et v
!         06-06-20  sigma_w (double precision) remplace dsig_t qui devient real*4
!...............................................................................
!    _________                    .__                  .__                     ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

      if(flag_nemoffline==1)return !24-03-14

!...........
! GEneric Modelling shortname conventions (m for "minus" & p for "plus")
! defined in set_parameters.F routine
      ip=ipt ; im=imt ; jp=jpt ; jm=jmt ; kp=kpt ; km=kmt
!...........

       if(case_==1) then !11111111111111>

! Here we have the choice between the "Arbitrary Lagrangian Eulerian" (ALE)
! approach provided that ale_selected=1 and the conventional sigma approach
! where dz depends on the SLA variations only provided that ale_selected=0
        if(ale_selected==1) then !--->
          call cellbox_thickness_ale
        else                     !--->
          call cellbox_thickness_sigma(2)
        endif                    !--->

       return
       endif              !11111111111111>


!              /    /   /


! Initial values:
       if(case_==0) then !>>>>>>>>>>>>>>>>

         dz_t(:,:,kmax+1,:)=small1
         dz_t(:,:,0     ,:)=small1
         dz_u(:,:,kmax+1,:)=small1
         dz_u(:,:,0     ,:)=small1
         dz_v(:,:,kmax+1,:)=small1
         dz_v(:,:,0     ,:)=small1

         do t_=tstr_,tend_
          call cellbox_thickness_sigma(t_)
         enddo

!        dz_t(:,:,kmax+1,2)=small1
!        dz_t(:,:,0     ,2)=small1
!        lb4=lbound(dz_t) 
!        do k=lb4(4),1
!         dz_t(:,:,:,k)=dz_t(:,:,:,2)
!        enddo

!        dz_u(:,:,kmax+1,2)=small1
!        dz_u(:,:,0     ,2)=small1
!        lb4=lbound(dz_u) 
!        do k=lb4(4),1
!         dz_u(:,:,:,k)=dz_u(:,:,:,2)
!        enddo

!        dz_v(:,:,kmax+1,2)=small1
!        dz_v(:,:,0     ,2)=small1
!        lb4=lbound(dz_v) 
!        do k=lb4(4),1
!         dz_v(:,:,:,k)=dz_v(:,:,:,2)
!        enddo

       return
       endif              !>>>>>>>>>>>>>>>>

!              /    /   /

! Updates arrays at the end of each iterative round:
       if(case_==2) then !22222222222222>
       stop 'ne pas passer par cellbox(2) carupdate dz dans moveforward'!01-06-10
       return
       endif              !2222222222222>

      end subroutine cellbox_thickness

!.....................................................................

      subroutine cellbox_thickness_ale
      stop 'cellbox_thickness_ale n est pas A jour' !28-01-16
      end subroutine cellbox_thickness_ale

!.......................................................................

      subroutine cellbox_thickness_sigma(t_)
      use module_principal
      use module_parallele
      implicit none
      integer t_,thz_

       lb3=lbound(hz_w) ; thz_=max(t_,lb3(3))
       do k=1,kmax
       do j=0,jmax+1
       do i=0,imax+1
!       dz_t(i,j,k,t_)=hz_w(i,j,thz_)*dsig_t(i,j,k)
! on remplace dsig_t par un delta de sigma_w car dsig_t passe en real*4
        dz_t(i,j,k,t_)=hz_w(i,j,thz_)*(sigma_w(i,j,k+1)-sigma_w(i,j,k)) !06-06-20
       enddo
       enddo
       enddo

       do k=1,kmax
       do j=0,jmax+1
       do i=1,imax+1
        dz_u(i,j,k,t_)=0.5*(dz_t(i,j,k,t_)+dz_t(i-1,j,k,t_)) !05-06-20
       enddo
       enddo
       enddo

       do k=1,kmax
       do j=1,jmax+1
       do i=0,imax+1
        dz_v(i,j,k,t_)=0.5*(dz_t(i,j,k,t_)+dz_t(i,j-1,k,t_))
       enddo
       enddo
       enddo

      return

      stop 'cellbox_thickness_sigma grid_or_wgrid not recognized'
      end subroutine cellbox_thickness_sigma
