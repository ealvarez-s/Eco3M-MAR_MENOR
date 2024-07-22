










      subroutine obc_bio(ichoix)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 334 - last update: 04-03-22
!______________________________________________________________________
      use module_principal
      use module_parallele   !#mpi
      implicit none
      integer ichoix

! OUPUT: BIO_Z

!...............................................................................
! Version date      Description des modifications
! modifs: 25/11/01: prise en compte des petits fleuves
!         19/03/03: date de debut et de fin de calcul pour traceur
!                   traitement des sources déplacé dans un sous programme
!                   specifique: source_tracer.F
!         04/04/06: les dimensions de RIVER_BIO ont augmente
!         14/10/08: Modif sur cas2 (forcage externe): gradient nul quand le
!                   courant est sortant
!         16/02/09: Ce jour là je recupere pour symphonie2008 la routine
!                   utilisee par Pierre Amael dans symphonie2007 depuis le 14/10/08
! 2009.3  01-10-09: utilisation des nouveaux facteurs d'echelle verticale
! 2010.8  03-05-10  Terminologie & suppression biohz & nouveau schema forward
! 2010.10 18-06-10  Terminologie et parallelisation
! 2010.13 01-11-10  rivertrc_inout remplace river_dom
! 2010.25 26-02-12  suppression mbio nbio
!         28-03-12  echange compatible avec pgf Stelios
! S26.1   07-11-12  Possibilite de C.L. river par le haut
!         28-03-13  rap_biobc remplace rap_obc
!         06-02-14  ajout routine obc_bio_mpi
!         26-03-14  regroupement echanges mpi
!         19-04-14  river_inout renomme rivertrc_inout
!         14-07-14  nouveaux echanges
!         17-07-14  suite point precedent. ne pas passer par la routine si tableaux non alloues
!         20-07-14  ajout routine obc_bio_botsurf
!         25-07-14  Les obc particulieres sont deplacees dans subroutine usersbiobc.F90
!         25-12-14  Condition limite superieure pour fleuves dependante de flag_nemoffline
!         21-03-15  echange mpi fluxbio_w
!         08-12-15  prise en compte obcstatus
!         25-04-16  Cas du modele 1DV: !25-04-16
!         14-12-16  Un commentaire
! v334    04-03-22  fluxbio_w(i,j,vb,2)=abs(riverflux(kr,1))*inv_dxdy_t(i,j) & ! W (m/s)
!...............................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


!                     /     /     /

! Cas du modele 1DV: !25-04-16
      if(flag_1dv==1) then !1D1D1D>
        ub4=ubound(bio_t) ; lb4=lbound(bio_t)
        do vb=lb4(4),ub4(4) ; do k=lb4(3),ub4(3) ; do j=lb4(2),ub4(2) ; do i=lb4(1),ub4(1)
         bio_t(i,j,k,vb)=bio_t(2,2,k,vb)
        enddo ; enddo ; enddo ; enddo
       return
      endif                !1D1D1D>


!***************************************************************************
! CONDITIONS LIMITES OUVERTES
! DEBUT:
      if(ichoix==2) then                                             !07/04/06
!***************************************************************************

!........................................................
! Leap frog ou forward?
!     if(itimebio.eq.0) then ! CAS LEAPFROG                            !14/10/08
!          k1=1              ! dernier indice VHZ_X VHZ_Y
!     else                   ! CAS FORWARD
!          k1=-1
!     endif

! L'aiguillage est fonction du signe du courant sur la frontiere
      do k=1,kmax
       do i=1,imax

        anyv3d(i,0     ,k,1)=                                           &
                   max(sign(un, veldxdz_v(i,1     ,k,1)) ,zero)
        anyv3d(i,jmax+1,k,1)=                                           &
                   max(sign(un,-veldxdz_v(i,jmax+1,k,1)) ,zero)

        anyv3d(i,0     ,k,2)=1.-anyv3d(i,0     ,k,1)
        anyv3d(i,jmax+1,k,2)=1.-anyv3d(i,jmax+1,k,1)

       enddo

       do j=1,jmax

        anyv3d(0     ,j,k,1)=                                           &
                    max(sign(un, veldydz_u(1     ,j,k,1)) ,zero)
        anyv3d(imax+1,j,k,1)=                                           &
                    max(sign(un,-veldydz_u(imax+1,j,k,1)) ,zero)

        anyv3d(0     ,j,k,2)=1.-anyv3d(0     ,j,k,1)
        anyv3d(imax+1,j,k,2)=1.-anyv3d(imax+1,j,k,1)

       enddo

      enddo
!........................................................

       x2=    rap_biobc !28-03-13

!      write(6,*)'rap_biobc=',rap_biobc

       x1=(1.-rap_biobc)
!ccccccX3=1.    ! extrapolation par dérivée premiere nulle
       x3=2.    ! extrapolation par dérivée seconde  nulle
       x4=1.-x3

       i2=min0(2,imax+1)
       j2=min0(2,jmax+1)
       i3=max0(0,imax-1)
       j3=max0(0,jmax-1)

      do vb=1,vbmax   ! debut de boucle sur vb variable bio


! La condition "OBC=etat initial" (cas 0) equivaut a ne rien faire

! Forcage externe (cas2):
      if(biobc_type(vb).eq.2) then !222222222222222>

        if(obcstatus(jeq1)==1) then !----> !08-12-15
         do k=1,kmax
          do i=1,imax
!                                                    nouveau schema le !14/10/08

           bio_t(i,0     ,k,vb)=                                        &
          anyv3d(i,0     ,k,1)*(                         & ! 1 si courant entrant sinon 0
                                x2*biobc_i_t(i,k,vb,1,2)                &
                               +x1*biobc_i_t(i,k,vb,1,1) )              &
         +anyv3d(i,0     ,k,2)*(                         & ! 1 si courant sortant sinon 0
        x3*bio_t(i,1     ,k,vb)                                         &
       +x4*bio_t(i,j2    ,k,vb) )

          enddo
         enddo
        endif                      !---->

        if(obcstatus(jeqjmax)==1) then !----> !08-12-15
         do k=1,kmax
          do i=1,imax

           bio_t(i,jmax+1,k,vb)=                                        &
          anyv3d(i,jmax+1,k,1)*(                         & ! 1 si courant entrant sinon 0
                                x2*biobc_i_t(i,k,vb,2,2)                &
                               +x1*biobc_i_t(i,k,vb,2,1) )              &
         +anyv3d(i,jmax+1,k,2)*(                         & ! 1 si courant sortant sinon 0
        x3*bio_t(i,jmax  ,k,vb)                                         &
       +x4*bio_t(i,j3    ,k,vb) )

          enddo
         enddo
        endif                      !---->

        if(obcstatus(ieq1)==1) then !----> !08-12-15
         do k=1,kmax
          do j=1,jmax
!                                                    nouveau schema le !14/10/08
           bio_t(0     ,j,k,vb)=                                        &
          anyv3d(0     ,j,k,1)*(                                        &
                                x2*biobc_j_t(j,k,vb,1,2)                &
                               +x1*biobc_j_t(j,k,vb,1,1) )              &
         +anyv3d(0     ,j,k,2)*(                                        &
        x3*bio_t(1     ,j,k,vb)                                         &
       +x4*bio_t(i2    ,j,k,vb) )

          enddo
         enddo
        endif                      !---->

        if(obcstatus(ieqimax)==1) then !----> !08-12-15
         do k=1,kmax
          do j=1,jmax

           bio_t(imax+1,j,k,vb)=                                        &
          anyv3d(imax+1,j,k,1)*(                                        &
                                x2*biobc_j_t(j,k,vb,2,2)                &
                               +x1*biobc_j_t(j,k,vb,2,1) )              &
         +anyv3d(imax+1,j,k,2)*(                                        &
        x3*bio_t(imax  ,j,k,vb)                                         &
       +x4*bio_t(i3    ,j,k,vb) )

          enddo
         enddo
        endif                      !---->

      endif                        !222222222222222>

! Condition de gradient nul (cas1):
      if(biobc_type(vb).eq.1) then !111111111111111>
         do k=1,kmax
          do i=1,imax
           bio_t(i,0     ,k,vb)=bio_t(i,1   ,k,vb)
           bio_t(i,jmax+1,k,vb)=bio_t(i,jmax,k,vb)
          enddo
          do j=1,jmax
           bio_t(0     ,j,k,vb)=bio_t(1   ,j,k,vb)
           bio_t(imax+1,j,k,vb)=bio_t(imax,j,k,vb)
          enddo
         enddo
      endif                        !111111111111111>

      if(biobc_type(vb)>2)call usersbiobc !25-07-14

! Rangs supplementaires pour schemas d'advection "3 points" !14-07-14
         do k=1,kmax
          do i=1,imax
           bio_t(i,-1    ,k,vb)=bio_t(i,1   ,k,vb)
           bio_t(i,jmax+2,k,vb)=bio_t(i,jmax,k,vb)
          enddo
          do j=1,jmax
           bio_t(-1    ,j,k,vb)=bio_t(1   ,j,k,vb)
           bio_t(imax+2,j,k,vb)=bio_t(imax,j,k,vb)
          enddo
         enddo

      enddo           ! fin de boucle sur vb variable bio

!     call obc_bio_mpi('z1') !06-02-14
! Echange 'zb' pour schemas d'advection "3 points" !14-07-14
      call obc_bio_mpi('zb') !07-02-14

!***************************************************************************
! CONDITIONS LIMITES OUVERTES
! FIN.
      return
      endif
!***************************************************************************


!                     /     /     /


!***************************************************************************
! CONDITIONS LIMITES: FLEUVES
! DEBUT
      if(ichoix==1) then
!***************************************************************************

      do kr=1,nriver

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
!     if(river_dom(kr) == par%rank) then !rrrrrrrrrrrrrrrrrr>
      if(rivertrc_inout(kr)==1)        then !rrrrrrrrrrrrrrrrrr>          !01-11-10
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      i=iriver(kr,1)
      j=jriver(kr,1)
      if(riverdir(kr)>=1.and.riverdir(kr)<=4) then !-laterale->

       if(flag_nemoffline/=1) then !ooo>
! cas general:
        do vb=1,vbmax
         do k=kmin_w(i,j),kmax 
          bio_t(i,j,k,vb)=river_bio(vb,kr,1)                          !04/04/06
         enddo  ! k
        enddo  ! vb
       else                        !ooo>
! cas nemo:
        bio_t(i,j,:,:)=0. 
       endif                       !ooo>

      endif                                        !-laterale->

      if(riverdir(kr)==0) then                     !-surface-> !04-03-22
        do vb=1,vbmax
         fluxbio_w(i,j,vb,2)=abs(riverflux(kr,1))*invdxdy_t(i,j) & ! W (m/s)
                             *river_bio(vb,kr,1)                   ! BIO
        enddo  ! vb
      endif                                        !-surface->

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
      endif                              !rrrrrrrrrrrrrrrrrr>
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       enddo

!***************************************************************************
! CONDITIONS LIMITES: FLEUVES
! FIN
      return
      endif
!***************************************************************************


!                     /     /     /


!***************************************************************************
! CONDITIONS LIMITES: INITIALISATION
! DEBUT
      if(ichoix==0) then                                             !07/04/06
!***************************************************************************

! On s'assure de l'initialisation de la couronne pour que la condition
! au limite ouverte "persistence de l'etat initial" puisse fonctionner
! correctement

         do vb=1,vbmax
         do k=1,kmax
          do i=1,imax
           bio_t(i,0     ,k,vb)=bio_t(i,1   ,k,vb)
           bio_t(i,jmax+1,k,vb)=bio_t(i,jmax,k,vb)
          enddo
          do j=0,jmax+1
           bio_t(0     ,j,k,vb)=bio_t(1   ,j,k,vb)
           bio_t(imax+1,j,k,vb)=bio_t(imax,j,k,vb)
          enddo
         enddo
         enddo


!***************************************************************************
! CONDITIONS LIMITES: INITIALISATION
! FIN
      return
      endif
!***************************************************************************


      end subroutine obc_bio

!...............................................................

      subroutine obc_bio_mpi(txt_)
      use module_principal
      use module_parallele
      implicit none
      integer loop_
      character*2 txt_

      if(.not.allocated(bio_t))return !17-07-14

!!$! Nouvelle methode avec choix des voisins
                    if(txt_=='z1')           &
       call get_type_echange('z1','bio_z1'   &
                                  ,bio_t     &
                           ,lbound(bio_t)    &
                           ,ubound(bio_t),k0)

                    if(txt_=='z0')           &
       call get_type_echange('z0','bio_z0'   &
                                  ,bio_t     &
                           ,lbound(bio_t)    &
                           ,ubound(bio_t),k0)

                    if(txt_=='zb')           & !z1z2
       call get_type_echange('zb','bio_zb'   &
                                  ,bio_t     &
                           ,lbound(bio_t)    &
                           ,ubound(bio_t),k0)

                    if(txt_=='zc')           & !z0z1z2
       call get_type_echange('zc','bio_zc'   &
                                  ,bio_t     &
                           ,lbound(bio_t)    &
                           ,ubound(bio_t),k0)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(bio_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_bio_mpi

!.....................................................................................

      subroutine obc_bio_botsurf !20-07-14
      use module_principal
      implicit none

! Notons qu'avec mpi, la boucle sur vb est comparable en taille a la boucle
! sur i. Autant la mettre en boucle interieure.

! C.L. sous le fond pour le cas (toujours possible) de la grille sigma step
      do j=1,jmax ; do i=1,imax
       do k=1,kmin_w(i,j)-1
          do vb=1,vbmax
           bio_t(i,j,k,vb)=bio_t(i,j,kmin_w(i,j),vb)
          enddo
       enddo
      enddo ; enddo

! Surface kmax+1
! Note qu'il s'agit d'une C.L. par default qui peut etre ecrasee par
! une C.L. pour un flux de surface
      do j=1,jmax ; do i=1,imax
      if(omega_w(i,j,kmax+1,1)>=0.) then !------>
          do vb=1,vbmax
           bio_t(i,j,kmax+1,vb)=bio_t(i,j,kmax,vb) ! Gradient Nul
          enddo
      else                               !------>
          do vb=1,vbmax
           bio_t(i,j,kmax+1,vb)=0.                 ! L'eau si pure...
          enddo
      endif                              !------>
      enddo ; enddo

! Quelle consequence d'imposer 0 en kmax+1 en ce qui concerne le
! flux a omega(kmax) avec omega<0?
! La valeur du traceur donnee par schema up3 est
! (T(kmax-1)+T(kmax))/2+(1/6)*(-0+2T(kmax)-T(kmax-1))


      end subroutine obc_bio_botsurf

!.....................................................................................

      subroutine obc_bio_mpi_fluxbio(surf_or_bot_,exchangetype_) !21-03-15
      use module_principal ; use module_parallele ; use module_biology
      implicit none
      character(len=*),intent(in) :: surf_or_bot_  ! 'surface' or 'bottom'
      character(len=2),intent(in) :: exchangetype_ ! 'z0' etc..
      integer loop_
      flag=0

      if(surf_or_bot_=='surface') then !sssssssss>
       flag=1 ; k=2
      endif                            !sssssssss> 
      if(surf_or_bot_=='bottom') then  !bbbbbbbbb>
       flag=1 ; k=1
      endif                            !bbbbbbbbb> 
      write(texte30,'(a8,a2,i0)')'fluxbio_',exchangetype_,k

!!$! Nouvelle methode avec choix des voisins
      call get_type_echange(exchangetype_       &
                           ,trim(texte30)       &
                                  ,fluxbio_w    &
                           ,lbound(fluxbio_w)   &
                           ,ubound(fluxbio_w)   &
                           ,k                   &
                           ,k0)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(fluxbio_w,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      if(flag==0)stop 'obc_bio_mpi_fluxbio txt_ not recognized'
      end subroutine obc_bio_mpi_fluxbio
