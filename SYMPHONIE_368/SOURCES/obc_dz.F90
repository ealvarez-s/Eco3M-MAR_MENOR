      subroutine obc_dz
!______________________________________________________________________
! SYMPHONIE ocean model
! release 299 - last update: 18-03-21
!______________________________________________________________________
      use module_principal
      use module_parallele
      implicit none
      integer loop_,idi_dz_
#ifdef synopsis
       subroutinetitle='obc_dz'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ Description:
!$ dz boundary conditions for far external grid nodes. This is not a dynamic
!$ condition but a simple procedure required to compute
!$ half/half average values on _u nodes and _v nodes.


!...............................................................................
! Version date      Description des modifications
! 2009.3  05-12-09: mis en service
! S26     11-11-13: modif pour subcycling
!         24-03-14  Ne pas passer dans la C.L. si nemo offline
!         12-05-14  nouveaux echanges
! v299    18-03-21  utiliser obcstatus
!...............................................................................

      if(flag_nemoffline==0) then !>>>>>>>>>>>>>>>> !24-03-14

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       do k=1,kmax
        do j=1,jmax
         dz_t(0     ,j,k,2)=dz_t(1   ,j,k,2)
        enddo
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do k=1,kmax
        do i=0,imax+1
         dz_t(i,0     ,k,2)=dz_t(i,1   ,k,2)
        enddo
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       do k=1,kmax
        do j=1,jmax
         dz_t(imax+1,j,k,2)=dz_t(imax,j,k,2)
        enddo
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do k=1,kmax
        do i=0,imax+1
         dz_t(i,jmax+1,k,2)=dz_t(i,jmax,k,2)
        enddo
       enddo
      endif                                      !----------->

      endif                       !>>>>>>>>>>>>>>>>

! mpi boundaries:
! NOUVEAUX ECHANGES
#ifdef parallele
      call get_type_echange('z1','dz_z1_t2',dz_t,lbound(dz_t) &
                                                ,ubound(dz_t) &
                                                ,2,idi_dz_)    ! 'za' = 'z0' et 'z1'
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(dz_t,idi_dz_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif
! ANCIENS ECHANGES
!#ifdef parallele
!      ub4=ubound(dz_t) ; lb4=lbound(dz_t)
!      call echange('z1',dz_t,lb4,ub4,2) ! 2 pour echange 3eme arg = 2 !C.L. i=imax+1 i=0 j=jmax+1 j=0
!#endif

      end subroutine obc_dz
