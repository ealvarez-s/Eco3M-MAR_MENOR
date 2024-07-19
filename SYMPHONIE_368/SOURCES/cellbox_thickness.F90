      subroutine cellbox_thickness(case_,tstr_,tend_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 285 - last update: 06-06-20
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
      integer case_,tstr_,tend_,t_
#ifdef synopsis
       subroutinetitle='cellbox_thickness'
       subroutinedescription= &
      'Udapes cellbox thickness and stores the result in dz_t dz_u dz_v'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
#ifdef bidon
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='cellbox_thickness_ale'
       subroutinedescription= &
      'Udapes cellbox thickness dz_t using ale method'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ time scale of the restoring term:
       const1=dti_lp/(10.*86400.) ! new lagragian approach !26-09-10

!$ horizontal diffusion coef:
       const2=0.05                                        !26-09-10

       do j=1,jmax
        do i=1,imax
         if(mask_t(i,j,kmax)==1)then !mskmskmskmskmskmsk>                  !16-12-09

         sum0=0.
         sum1=small1
         sum2=small1
         do k=kmin_w(i,j),kmax                                          !08-10-09

!$ vertical scale factor of reference (initial grid):
         x0=hz_w(i,j,after)*dsig_t(i,j,k)

!$ First guess = adiabatic case + constraint 0.8*dz(t0) < dz(after) < 1.2*dz(t0)
!$                              + nudging toward the initial grid
!$                              + time filter
!$                              + spatial filter
                       dz_t(i   ,j   ,k   ,after)=               &

! min and max values:
       min(cellboxfactor2*x0,max(cellboxfactor1*x0,              &    !11-10-10

! Time filter: (Note que si before3=2 cela oblige à calculer le filre maintenant
! et non pas après l'étape d'ajustement à la contrainte barotrope)
        (tfd0*dz_t(i   ,j   ,k   ,now)                           &
        +tfd1*dz_t(i   ,j   ,k   ,before)                        &
        +tfd2*dz_t(i   ,j   ,k   ,before2)                       &
        +tfd3*dz_t(i   ,j   ,k   ,before3)                       &

! horizontal velocity divergence:
        -dti_lp*( veldydz_u(i+ip,j   ,k   ,now)                  &
                 -veldydz_u(i-im,j   ,k   ,now)                  &
                 +veldxdz_v(i   ,j+jp,k   ,now)                  &
                 -veldxdz_v(i   ,j-jm,k   ,now)    )/dxdy_t(i,j) &


! Restoring term:
        +const1*( x0 - dz_t(i   ,j   ,k   ,before))              &


! Spatial filter:
        -const2*(4.*dz_t(i  ,j  ,k,before)                       &   !13-02-10
                   -dz_t(i+1,j  ,k,before)                       &
                   -dz_t(i-1,j  ,k,before)                       &
                   -dz_t(i  ,j+1,k,before)                       &
                   -dz_t(i  ,j-1,k,before)                       &
                -4.*hz_w(i  ,j    ,before)*dsig_t(i  ,j  ,k)     &
                   +hz_w(i+1,j    ,before)*dsig_t(i+1,j  ,k)     &
                   +hz_w(i-1,j    ,before)*dsig_t(i-1,j  ,k)     &
                   +hz_w(i  ,j+1  ,before)*dsig_t(i  ,j+1,k)     &
                   +hz_w(i  ,j-1  ,before)*dsig_t(i  ,j-1,k))    &

                     )*tfd4                                      &

          ))

!$ adjustment function if sum(dz) > h+ssh
          anyv1d(k,1)=max(zero,dz_t(i,j,k,after)-x0)
!$ adjustment function if sum(dz) < h+ssh
          anyv1d(k,2)=max(zero,x0-dz_t(i,j,k,after))

!$ Check that sum(dz)=h+ssh:
          sum0=sum0+dz_t(i,j,k,after)
          sum1=sum1+anyv1d(k,1)
          sum2=sum2+anyv1d(k,2)

         enddo

!$ computes x1, the Lagrange coef ensuring that sum(dz)=h+ssh:
         if(sum0>hz_w(i,j,after)) then  !------>
          k1=1
          x1=(hz_w(i,j,after)-sum0)/sum1
         else                           !------>
          k1=2
          x1=(hz_w(i,j,after)-sum0)/sum2
         endif                          !------>

!$ adjusts dz so that sum(dz)=h+ssh:
!        sum0=0.
         do k=kmin_w(i,j),kmax
          dz_t(i,j,k,after)=dz_t(i,j,k,after)+x1*anyv1d(k,k1)
!         sum0=sum0+dz_c(i,j,k,after)
         enddo

         endif                    !mskmskmskmskmskmsk>
        enddo
       enddo

!$ Open boundaries:
      call obc_dz

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
       dz_u(i,j,k,after)=0.5*(dz_t(i,j,k,after)+dz_t(i-1,j,k,after))
       dz_v(i,j,k,after)=0.5*(dz_t(i,j,k,after)+dz_t(i,j-1,k,after))
      enddo
      enddo
      enddo

      if(fgrid_or_wgrid==fgrid_case) &
      stop 'cellbox_thickness_ale grid_or_wgrid=fgrid_case'

#endif
      end subroutine cellbox_thickness_ale

!.......................................................................

      subroutine cellbox_thickness_sigma(t_)
      use module_principal
      use module_parallele
      implicit none
      integer t_,thz_
#ifdef synopsis
       subroutinetitle='cellbox_thickness_ale'
       subroutinedescription= &
      'Udapes cellbox thickness dz using sigma method'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
