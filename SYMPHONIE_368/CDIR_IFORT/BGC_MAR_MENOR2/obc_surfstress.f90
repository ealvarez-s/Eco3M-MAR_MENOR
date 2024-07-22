










      subroutine obc_surfstress
!______________________________________________________________________
! S model
! release S.26 - last update: 31-07-12
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none

!...............................................................................
! Version date      Description des modifications
! 2010.12 25-09-10  Mis en service
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S.26    31-07-14  Nouveaux echanges
!...............................................................................


      do j=2,jmax-1
       wstress_u(imax+1,j,1)=wstress_u(imax,j,1)
       wstress_u(1     ,j,1)=wstress_u(2   ,j,1)
      enddo
      do i=1,imax+1
       wstress_u(i,jmax,1)=wstress_u(i,jmax-1,1)
       wstress_u(i,1   ,1)=wstress_u(i,2     ,1)
      enddo
      do i=2,imax-1
       wstress_v(i,jmax+1,1)=wstress_v(i,jmax,1)
       wstress_v(i,1     ,1)=wstress_v(i,2   ,1)
      enddo
      do j=2,jmax+1
       wstress_v(imax,j,1)=wstress_v(imax-1,j,1)
       wstress_v(1   ,j,1)=wstress_v(2     ,j,1)
      enddo

!#ifdef 1
!      lb3=lbound(wstress_u) ; ub3=ubound(wstress_u)
!      call echange('x ',wstress_u,lb3,ub3,1) ! 2 pour echange 4eme arg = 2 ! C.L. i=imax+1 i=1 j=jmax j=1
!      lb3=lbound(wstress_v) ; ub3=ubound(wstress_v)
!      call echange('y ',wstress_v,lb3,ub3,1) ! 2 pour echange 4eme arg = 2 ! C.L. i=imax+1 i=1 j=jmax j=1
!#endif
      call get_type_echange('u1','wstress_u_u1_1',wstress_u,lbound(wstress_u),ubound(wstress_u),1,i1)
      call get_type_echange('v1','wstress_v_v1_1',wstress_v,lbound(wstress_v),ubound(wstress_v),1,i2)
      do loop1=1,subcycle_exchange
        call echange_voisin(wstress_u,i1,mpi_neighbor_list(loop1)) !31-07-14
        call echange_voisin(wstress_v,i2,mpi_neighbor_list(loop1)) !31-07-14
      enddo
      call loc_wait()

      end subroutine obc_surfstress
