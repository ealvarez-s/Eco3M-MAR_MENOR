      subroutine obc_depth(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release S.26 - last update: 17-11-18
!______________________________________________________________________
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! (°v°)
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! (°O°)
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      ! (°L°)
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     ! m°v°m
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal
      use module_parallele !#MPI
      implicit none
      integer case_,id_depth_,loop_
#ifdef synopsis
       subroutinetitle='obc_depth'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         09/03/09: mis en service
!         06/04/09: echange  depth_z
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
! 2010.8  09-05-10  schema d'advection "Modified up3" implique une condition
!                   "z2" sur depth_t
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S26     07-10-18  obc compatible avec diffusivite QUICK dans le plan horizontal
!         17-11-18  parce que le cas particulier des coins n'est pas traite et que 
!                   dans le cas des proc inutiles supprimes j'ai constate que cela 
!                   entraine une perte de continuite mpi, les conditions sont calculees 
!                   sans condition sur par%tvoisin
!...............................................................................

!................................................!
! Conditions aux limites laterales pour depth_t
      if(case_==2) then !2222222222222222222> !09-05-10


! West open boundary (if any) !26-11-14 
! test sur par%tvoin commentE le 17-11-18
!     if(par%tvoisin(ouest)==mpi_proc_null) then !----------->
       do k=1,kmax
       do j=0,jmax+1
         depth_t(-1,j,k)=depth_t(0,j,k)
       enddo
       enddo
!     endif                                      !----------->

! East open boundary (if any)
! test sur par%tvoin commentE le 17-11-18
!     if(par%tvoisin(est  )==mpi_proc_null) then !-----------> 
       do k=1,kmax
       do j=0,jmax+1
         depth_t(imax+2,j,k)=depth_t(imax+1,j,k)
       enddo
       enddo
!     endif                                      !----------->

! South open boundary (if any)
! test sur par%tvoin commentE le 17-11-18
!     if(par%tvoisin(sud  )==mpi_proc_null) then !----------->
       do k=1,kmax
       do i=-1,imax+2
         depth_t(i,-1,k)=depth_t(i,0,k)
       enddo
       enddo
!     endif                                      !----------->

! North open boundary (if any)
! test sur par%tvoin commentE le 17-11-18
!     if(par%tvoisin(nord )==mpi_proc_null) then !----------->
       do k=1,kmax
       do i=-1,imax+2
         depth_t(i,jmax+2,k)=depth_t(i,jmax+1,k)
       enddo
       enddo
!     endif                                      !----------->


      call get_type_echange('z2','depth_t_z2',depth_t     &
                                      ,lbound(depth_t)    &
                                      ,ubound(depth_t),id_depth_)
      do loop_=1, subcycle_exchange
         call echange_voisin(depth_t,id_depth_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() 


      return
      endif              !22222222222222222222>
!................................................!

!................................................!
! Conditions aux limites laterales pour PROF3D_Z
      if(case_==1) then !11111111111111111111>
      do k=0,kmax+1
       do i=1,imax
         depth_t(i,0     ,k)= depth_t(i,1   ,k)
         depth_t(i,jmax+1,k)= depth_t(i,jmax,k)
       enddo
       do j=0,jmax+1
         depth_t(0     ,j,k)= depth_t(1   ,j,k)
         depth_t(imax+1,j,k)= depth_t(imax,j,k)
       enddo
      enddo
      return
      endif              !11111111111111111111>
!................................................!

!................................................!
! Conditions aux limites laterales pour PROFWK_Z
      if(case_==0) then !00000000000000000000>
      do k=1,kmax+1
       do i=1,imax
         depth_w(i,0     ,k)= depth_w(i,1   ,k)
         depth_w(i,jmax+1,k)= depth_w(i,jmax,k)
       enddo
       do j=0,jmax+1
         depth_w(0     ,j,k)= depth_w(1   ,j,k)
         depth_w(imax+1,j,k)= depth_w(imax,j,k)
       enddo
      enddo

!#ifdef parallele
!       ub3=ubound(depth_w) ; lb3=lbound(depth_w)
!       call echange('z0',depth_w,lb3,ub3) !C.L. i=imax   i=1 j=jmax   j=1
!       call echange('z1',depth_w,lb3,ub3) !C.L. i=imax+1 i=0 j=jmax+1 j=0
!#endif
#ifdef parallele
      call get_type_echange('za','depth_w_za',depth_w,lbound(depth_w),ubound(depth_w),i1)
      do loop1=1,subcycle_exchange
        call echange_voisin(depth_w,i1,mpi_neighbor_list(loop1)) !31-07-14
      enddo
      call loc_wait()
#endif

      return
      endif              !00000000000000000000>
!................................................!

      end subroutine obc_depth

!......................................................................

      subroutine obc_depth_f(case_)
      use module_principal ; use module_parallele
      implicit none
      integer case_,loop_
#ifdef synopsis
       subroutinetitle='obc_depth_f'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do k=1,kmax+1
       do i=1,imax+1
         depth_f(i,0     ,k)= depth_f(i,1     ,k)
         depth_f(i,jmax+2,k)= depth_f(i,jmax+1,k)
       enddo
       do j=0,jmax+2
         depth_f(0     ,j,k)= depth_f(1     ,j,k)
         depth_f(imax+2,j,k)= depth_f(imax+1,j,k)
       enddo
      enddo
#ifdef parallele
      call get_type_echange('r1','depth_f_r1',depth_f,lbound(depth_f),ubound(depth_f),i1)
      do loop_=1,subcycle_exchange
        call echange_voisin(depth_f,i1,mpi_neighbor_list(loop_)) !31-07-14
      enddo
      call loc_wait()
#endif


      end subroutine obc_depth_f
