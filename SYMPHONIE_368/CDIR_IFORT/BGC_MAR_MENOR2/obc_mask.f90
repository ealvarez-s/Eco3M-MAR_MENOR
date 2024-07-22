










      subroutine obc_mask(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 299 - last update: 18-03-21
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      implicit none
      integer case_

!...............................................................................
!         16/10/08: mise en service
!         06/03/09: regroupement des obc pour faciliter
!                   la parallelisation
!         09/03/09: tridimensionnaliser le masque
!         02/04/09: ajout d'une 2eme couronne extraperipherique
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S.26    16-07-14  echange mpi separee des c.l.
!         12-04-15  echange z3 requis pour modele IQS
!         19-06-17  ajout subroutine obc_mask2d_mpi
! v299    18-03-21  utiliser obcstatus
!...............................................................................
!  _________                    .__                  .__              ! m[°v°]m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____  
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!...............................................................................



      if(case_==1) then !11111111111111111>

! C.L. sur rang 1

      if(obcstatus(ieq1)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do j=1,jmax
          mask_t(0     ,j,k)=  mask_t(1   ,j,k)
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(ieqimax)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do j=1,jmax
          mask_t(imax+1,j,k)=  mask_t(imax,j,k)
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(jeq1)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do i=0,imax+1
          mask_t(i,0     ,k)=  mask_t(i,1   ,k)
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(jeqjmax)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do i=0,imax+1
          mask_t(i,jmax+1,k)=  mask_t(i,jmax,k)
       enddo
       enddo
      endif                                      !----------->

      return
      endif                !11111111111111111>


      if(case_==3) then !333333333333333333333>

! C.L. Sur rang 2

      if(obcstatus(ieq1)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do j=1,jmax
          mask_t(-1    ,j,k)=  mask_t(1   ,j,k)                        !02/04/09
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(ieqimax)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do j=1,jmax
          mask_t(imax+2,j,k)=  mask_t(imax,j,k)
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(jeq1)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do i=0,imax+1
          mask_t(i,-1    ,k)=  mask_t(i,1   ,k)                        !02/04/09
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(jeqjmax)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do i=0,imax+1
          mask_t(i,jmax+2,k)=  mask_t(i,jmax,k)
       enddo
       enddo
      endif                                      !----------->

      return
      endif                !3333333333333333333333>

      if(case_.eq.2) then !22222222222222222>

! Sécurité aux frontières ouvertes:
! les points sur la frontiere sont masqués si ils sont immediatement
! voisin d'un point en terre:
      if(obcstatus(jeq1)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do i=1,imax                                                   !30/10/01
        if( mask_t(i,2     ,k).eq.0)  mask_t(i,1   ,k)=0
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(jeqjmax)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do i=1,imax                                                   !30/10/01
        if( mask_t(i,jmax-1,k).eq.0)  mask_t(i,jmax,k)=0
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(ieq1)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do j=1,jmax                                                   !30/10/01
       if( mask_t(2     ,j,k).eq.0)  mask_t(1   ,j,k)=0
       enddo
       enddo
      endif                                      !----------->

      if(obcstatus(ieqimax)==1) then !----------->
       do k=1,kmax+1                                                        !09/03/09
       do j=1,jmax                                                   !30/10/01
        if( mask_t(imax-1,j,k).eq.0)  mask_t(imax,j,k)=0
       enddo
       enddo
      endif                                      !----------->

      return
      endif                !22222222222222222>

!______________________________________________________________________________


      if(case_.eq.0) then !00000000000000000>
      stop 'obc_morpho cas 0 pas prevu pour parallelisation'
      return
      endif                !00000000000000000>


      end subroutine obc_mask

!.............................................................................

      subroutine obc_mask_mpi
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer loop_,id_mask_t_zc_,id_mask_t_z3_

! Echange zc requis par schema d'advection_bioup2up3
! Echange z3 requis par modele IQS


!!$! Nouvelle methode avec choix des voisins
! zb = z1 & z2
       call get_type_echange('zc','mask_t_zc'  &
                                  ,mask_t      &
                           ,lbound(mask_t)     &
                           ,ubound(mask_t),id_mask_t_zc_)
       call get_type_echange('z3','mask_t_z3'  &
                                  ,mask_t      &
                           ,lbound(mask_t)     &
                           ,ubound(mask_t),id_mask_t_z3_)

      ! Echanges
      do loop_=1, subcycle_exchange
       call echange_voisin(mask_t,id_mask_t_zc_,mpi_neighbor_list(loop_))
       call echange_voisin(mask_t,id_mask_t_z3_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_mask_mpi

!.............................................................................

      subroutine obc_mask2d_mpi(txt_,k_)
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer k_,loop_
      character*2 txt_


!!$! Nouvelle methode avec choix des voisins
! zb = z1 & z2
       write(texte30,'(a,a,a,i0)')'mask_',txt_,'_',k_
       call get_type_echange(txt_,trim(texte30) &
                                  ,mask_t       &
                           ,lbound(mask_t)      &
                           ,ubound(mask_t)      &
                           ,k_                  &
                           ,k0)

      ! Echanges
      do loop_=1, subcycle_exchange
       call echange_voisin(mask_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_mask2d_mpi
