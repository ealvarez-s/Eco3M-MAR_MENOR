










      subroutine obc_h(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 338 - last update: 20-03-22
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer case_

!...............................................................................
! Version date      Description des modifications
!         09/03/09: mis en service
!         24/03/09: parallelisation
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S.26    02-07-14  nouveaux echanges
!         26-10-14  ajout subroutine obc_h_xyt1
!         12-11-14  ajout routine obc_h_xy_t_z0
!         09-01-15  OBC et mpiBC separees dans 2 subroutines differentes
! v338    20-03-22  elargir la condition de gradient nul sur H pour un
!                   meilleur fonctionnement des conditions radiatives
!...............................................................................
!    _________                    .__                  .__                     !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

! Condition de gradient nul sur les bords:
      do j=1,jmax
!        h_w(0,j)     =h_w(1,j)
!        h_w(imax+1,j)=h_w(imax,j)
         h_w(0:1,j)        =h_w(2,j)      !20-03-22
         h_w(imax:imax+1,j)=h_w(imax-1,j) !20-03-22
      enddo
      do i=0,imax+1
!        h_w(i,0)     =h_w(i,1)
!        h_w(i,jmax+1)=h_w(i,jmax)
         h_w(i,0:1)        =h_w(i,2)      !20-03-22
         h_w(i,jmax:jmax+1)=h_w(i,jmax-1) !20-03-22
      enddo

! mpi continuity:
      call obc_h_mpi_za !09-01-15

      end subroutine obc_h
!..........................................................................
      subroutine obc_h_mpi_za !09-01-15
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer loop_,var_id_

!!$! Nouvelle methode avec choix des voisins !02-07-14
! za = z0 & z1
       call get_type_echange('za','h_w_za',h_w     &
                                   ,lbound(h_w)    &
                                   ,ubound(h_w)    &
                                   ,var_id_)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(h_w,var_id_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_h_mpi_za
!..........................................................................
      subroutine obc_h_f
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer loop_

      do j=1,jmax+1
         h_f(0     ,j)=h_f(1     ,j)
         h_f(imax+2,j)=h_f(imax+1,j)
      enddo
      do i=0,imax+2
         h_f(i,0)     =h_f(i,1)
         h_f(i,jmax+2)=h_f(i,jmax+1)
      enddo
       call get_type_echange('r1','h_f_r1',h_f     &
                                   ,lbound(h_f)    &
                                   ,ubound(h_f)    &
                                   ,k0)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(h_f,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_h_f
!..........................................................................

!..........................................................................

      subroutine obc_h_xy_t_z0(k_) !12-11-14
      use module_principal
      use module_parallele      
      implicit none
      integer loop_,k_
      write(texte30,'(a10,i0)')'xy_t_z0_k=',k_
      call get_type_echange('z0',trim(texte30),xy_t        &
                                       ,lbound(xy_t)       &
                                       ,ubound(xy_t),k_,k0)
      ! Echanges
      do loop_=1,subcycle_exchange
       call echange_voisin(xy_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
      end subroutine obc_h_xy_t_z0
!..........................................................................

