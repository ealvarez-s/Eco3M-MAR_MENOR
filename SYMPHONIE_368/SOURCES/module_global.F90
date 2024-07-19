      module module_global
      implicit none
!______________________________________________________________________
! SYMPHONIE ocean model
! release S26 - last update: 08-05-18
!______________________________________________________________________
!...............................................................................
!    _________                    .__                  .__             ! m[°v°]m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
! version date      Description des modifications
!         06/04/09: mise en service
! 2009.3  16-10-09: ajout glob_lon et glob_lat
! S26     08-05-18  declaration glob_mask kind=1
!                   declaration glob_h real
!______________________________________________________________________

!     double precision,dimension(:,:),allocatable ::                    &
      real,dimension(:,:),allocatable ::                                &
        glob_h                                                          &
       ,glob_h_u                                                        &
       ,glob_h_v                                                        &
       ,glob_u_x                                                        &
       ,glob_v_y                                                        &
       ,glob_sf0_r                                                      &
       ,glob_sf1_r                                                      &
       ,glob_v_i                                                        &
       ,glob_u_j                                                        &
       ,glob_stf_i                                                      &
       ,glob_stf_j                                                      &
       ,glob_lon                                                        &
       ,glob_lat

      real*4,dimension(:,:),allocatable ::                              &
        glob_anyvar2d

!     integer,dimension(:,:),allocatable ::                             &
      integer(kind=1),dimension(:,:),allocatable ::                     &
        glob_mask  ,                                                    &
        glob_kmin_w

contains

      subroutine allocate_global(txt_l1,txt_l2,dim_l1,dim_l2,dim_l3)
      use module_parameter
      implicit none
      integer ki1,dim_l1,dim_l2,dim_l3
      character txt_l1*1,txt_l2*20
#ifdef synopsis
       if(txt_l1=='a') then
        call main_synopsis('allocate_global','allocate global arrays')
       else
        call main_synopsis('allocate_global','Deallocate global arrays')
       endif
#endif

      if(txt_l2(1:8) =='glob_v_y') then  !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_v_y (1:dim_l1+2,1:dim_l2+2))
      glob_v_y=0.                                     ! 16-12-09
       if(size(glob_v_y(:,1))/=iglb+2)stop 'stop bad size glob_v_y'
       if(size(glob_v_y(1,:))/=jglb+2)stop 'stop bad size glob_v_y'
      endif

      if(txt_l1=='d')deallocate(glob_v_y)

      return
      endif                               !>>>>>>

      if(txt_l2(1:8) =='glob_u_x') then  !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_u_x (1:dim_l1+2,1:dim_l2+2))
      glob_u_x=0.                                                    ! 16-12-09
       if(size(glob_u_x(:,1))/=iglb+2)stop 'stop bad size glob_u_x'
       if(size(glob_u_x(1,:))/=jglb+2)stop 'stop bad size glob_u_x'
      endif

      if(txt_l1=='d')deallocate(glob_u_x)

      return
      endif                               !>>>>>>

      if(txt_l2(1:8) =='glob_h_v') then  !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_h_v (1:dim_l1+2,1:dim_l2+2))
      glob_h_v=0.                                                     ! 16-12-09
       if(size(glob_h_v(:,1))/=iglb+2)stop 'stop bad size glob_h_v'
       if(size(glob_h_v(1,:))/=jglb+2)stop 'stop bad size glob_h_v'
      endif

      if(txt_l1=='d')deallocate(glob_h_v)

      return
      endif                               !>>>>>>

      if(txt_l2(1:8) =='glob_h_u') then  !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_h_u (1:dim_l1+2,1:dim_l2+2))
      glob_h_u=0.                                                     ! 16-12-09
       if(size(glob_h_u(:,1))/=iglb+2)stop 'stop bad size glob_h_u'
       if(size(glob_h_u(1,:))/=jglb+2)stop 'stop bad size glob_h_u'
      endif

      if(txt_l1=='d')deallocate(glob_h_u)

      return
      endif                               !>>>>>>

      if(txt_l2(1:10)=='glob_sf0_r') then !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_sf0_r(1:dim_l1+2,1:dim_l2+2))
      glob_sf0_r=0.                                                   ! 16-12-09
       if(size(glob_sf0_r(:,1))/=iglb+2)stop 'stop bad size glob_sf0_r'
       if(size(glob_sf0_r(1,:))/=jglb+2)stop 'stop bad size glob_sf0_r'
      endif

      if(txt_l1=='d')deallocate(glob_sf0_r)

      return
      endif                               !>>>>>>

      if(txt_l2(1:10)=='glob_sf1_r') then !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_sf1_r(1:dim_l1+2,1:dim_l2+2))
      glob_sf1_r=0.                                                   ! 16-12-09
       if(size(glob_sf1_r(:,1))/=iglb+2)stop 'stop bad size glob_sf1_r'
       if(size(glob_sf1_r(1,:))/=jglb+2)stop 'stop bad size glob_sf1_r'
      endif

      if(txt_l1=='d')deallocate(glob_sf1_r)

      return
      endif                               !>>>>>>

      if(txt_l2(1:11)=='glob_mask  ') then !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_mask   (0:dim_l1+1,0:dim_l2+1))
      glob_mask=0                                                    !16-12-2009
       if(size(glob_mask  (:,1))/=iglb+2)stop 'stop bad size glob_mask  '
       if(size(glob_mask  (1,:))/=jglb+2)stop 'stop bad size glob_mask  '
      endif

      if(txt_l1=='d')deallocate(glob_mask  )

      return
      endif                                 !>>>>>>

      if(txt_l2(1:8)=='glob_h  ') then     !----->

      if(txt_l1=='a') then
      allocate(glob_h   (0:dim_l1+1,0:dim_l2+1))
      glob_h=0.                                                        ! 16-12-09
       if(size(glob_h  (:,1))/=iglb+2)stop 'stop bad size glob_h  '
       if(size(glob_h  (1,:))/=jglb+2)stop 'stop bad size glob_h  '
      endif

      if(txt_l1=='d')deallocate(glob_h  )

      return
      endif                                 !----->


      if(txt_l2(1:13)=='glob_anyvar2d') then !>>>>>>

      if(txt_l1=='a')  then
      allocate(glob_anyvar2d (0:dim_l1+1,0:dim_l2+1))
      glob_anyvar2d=0.                                                ! 16-12-09
       if(size(glob_anyvar2d(:,1))/=iglb+2)stop 'stop bad glob_anyvar2d'
       if(size(glob_anyvar2d(1,:))/=jglb+2)stop 'stop bad glob_anyvar2d'
      endif

      if(txt_l1=='d')deallocate(glob_anyvar2d)

      return
      endif                                   !>>>>>>

      if(txt_l2(1:11)=='glob_kmin_w') then   !----->

      if(txt_l1=='a') then
      allocate(glob_kmin_w (0:dim_l1+1,0:dim_l2+1))
      glob_kmin_w=0                                                   ! 16-12-09
       if(size(glob_kmin_w(:,1))/=iglb+2)stop 'stop bad size glob_kmin_w'
       if(size(glob_kmin_w(1,:))/=jglb+2)stop 'stop bad size glob_kmin_w'
      endif

      if(txt_l1=='d')deallocate(glob_kmin_w)

      return
      endif                                   !----->

      if(txt_l2(1:8)=='glob_lon') then   !----->

      if(txt_l1=='a') then
      allocate(glob_lon (0:dim_l1+1,0:dim_l2+1))
      glob_lon=0.                                                      ! 16-12-09
       if(size(glob_lon(:,1))/=iglb+2)stop 'stop bad size glob_lon'
       if(size(glob_lon(1,:))/=jglb+2)stop 'stop bad size glob_lon'
      endif

      if(txt_l1=='d')deallocate(glob_lon)

      return
      endif                                   !----->

      if(txt_l2(1:8)=='glob_lat') then   !----->

      if(txt_l1=='a') then
      allocate(glob_lat (0:dim_l1+1,0:dim_l2+1))
      glob_lat=0.                                                      ! 16-12-09
       if(size(glob_lat(:,1))/=iglb+2)stop 'stop bad size glob_lat'
       if(size(glob_lat(1,:))/=jglb+2)stop 'stop bad size glob_lat'
      endif

      if(txt_l1=='d')deallocate(glob_lat)

      return
      endif                                   !----->


      write(*,*)'Erreur sur arguments passes dans allocate_global:'
      write(*,'(a,a)')'txt_l1=',txt_l1
      write(*,'(a,a)')'txt_l2=',txt_l2
      write(*,*)'dim_l1=',dim_l1
      write(*,*)'dim_l2=',dim_l2
      write(*,*)'dim_l3=',dim_l3
      stop ' STOP dans allocate_global.F90'    !05-07-10

      end subroutine allocate_global

      end module module_global
