










      subroutine obc_turbulence
      implicit none
!______________________________________________________________________
! SYMPHONIE ocean model
! release 299 - last update: 18-03-21
!______________________________________________________________________
!.......................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!.......................................................................
! Version date      Description des modifications
! v299    18-03-21  utiliser obcstatus
!.......................................................................
      end subroutine obc_turbulence

!.......................................................................

      subroutine obc_turbulence_tken
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer loop_,id_tken_

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
! i=1
       do k=1,kmax+1
        do j=2,jmax-1
         tken_w(1   ,j,k)=tken_w(2     ,j,k)
        enddo
       enddo

      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
! i=imax
       do k=1,kmax+1
        do j=2,jmax-1
         tken_w(imax,j,k)=tken_w(imax-1,j,k)
        enddo
       enddo

      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
! j=1
       do k=1,kmax+1
        do i=1,imax ! note: etendre la boucle donne une condition naturelle pour deux coins
          tken_w(i,1   ,k)=tken_w(i,2     ,k)
        enddo
       enddo

      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
! j=jmax
       do k=1,kmax+1
        do i=1,imax ! note: etendre la boucle donne une condition naturelle pour deux coins
          tken_w(i,jmax,k)=tken_w(i,jmax-1,k)
        enddo
       enddo

      endif                                      !----------->

! mpi boundaries:
      call get_type_echange('z0','tken_z0',tken_w,lbound(tken_w),ubound(tken_w),id_tken_)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(tken_w,id_tken_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_turbulence_tken

!.......................................................................

      subroutine obc_turbulence_tkea
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer loop_,id_tkea_

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
! i=1
       do k=1,kmax+1
        do j=2,jmax-1
         tkea_w(1   ,j,k)=tkea_w(2     ,j,k)
        enddo
       enddo

      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
! i=imax
       do k=1,kmax+1
        do j=2,jmax-1
         tkea_w(imax,j,k)=tkea_w(imax-1,j,k)
        enddo
       enddo

      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
! j=1
       do k=1,kmax+1
        do i=1,imax ! note: etendre la boucle donne une condition naturelle pour deux coins
          tkea_w(i,1   ,k)=tkea_w(i,2     ,k)
        enddo
       enddo

      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
! j=jmax
       do k=1,kmax+1
        do i=1,imax ! note: etendre la boucle donne une condition naturelle pour deux coins
          tkea_w(i,jmax,k)=tkea_w(i,jmax-1,k)
        enddo
       enddo

      endif                                      !----------->

! mpi boundaries:
      call get_type_echange('z0','tkea_z0',tkea_w,lbound(tkea_w),ubound(tkea_w),id_tkea_)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(tkea_w,id_tkea_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges


      end subroutine obc_turbulence_tkea

!.................................................................

      subroutine obc_turbulence_epsa
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer loop_,id_epsa_

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
! i=1
       do k=1,kmax+1
        do j=2,jmax-1
         epsa_w(1   ,j,k)=epsa_w(2     ,j,k)
        enddo
       enddo

      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
! i=imax
       do k=1,kmax+1
        do j=2,jmax-1
         epsa_w(imax,j,k)=epsa_w(imax-1,j,k)
        enddo
       enddo

      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
! j=1
       do k=1,kmax+1
        do i=1,imax ! note: etendre la boucle donne une condition naturelle pour deux coins
          epsa_w(i,1   ,k)=epsa_w(i,2     ,k)
        enddo
       enddo

      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
! j=jmax
       do k=1,kmax+1
        do i=1,imax ! note: etendre la boucle donne une condition naturelle pour deux coins
          epsa_w(i,jmax,k)=epsa_w(i,jmax-1,k)
        enddo
       enddo

      endif                                      !----------->

! mpi boundaries:
      call get_type_echange('z0','epsa_z0',epsa_w,lbound(epsa_w),ubound(epsa_w),id_epsa_)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(epsa_w,id_epsa_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_turbulence_epsa
