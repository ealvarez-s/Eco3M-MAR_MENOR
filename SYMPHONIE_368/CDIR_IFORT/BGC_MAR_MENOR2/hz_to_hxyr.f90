










      subroutine hz_to_hxyr
!______________________________________________________________________
! S model
! release S.26 - last update: 05-08-15
!______________________________________________________________________

      use module_principal
      implicit none

!...............................................................................
! modifs: 24/03/04: mis en service
!         09/03/09: La condition aux limites pour H est passée dans obc_h.F
! S.26    25-11-14: cas h issue des points f
!         09-04-15  etendre la couverture de h_v et h_u
!         05-08-15  ajout commentaire
!...............................................................................

       do j=1,jmax+1 ; do i=1,imax+1
         h_f(i,j)=0.25*(h_w(i,j-1)+h_w(i,j)+h_w(i-1,j-1)+h_w(i-1,j))
        enddo ; enddo

       do j=0,jmax+1 ; do i=1,imax+1
         h_u(i,j)= 0.5*(h_w(i-1,j)+h_w(i,j))
        enddo ; enddo

       do j=1,jmax+1 ; do i=0,imax+1
         h_v(i,j)= 0.5*(h_w(i,j-1)+h_w(i,j))
        enddo ; enddo

!      if(fgrid_or_wgrid==fgrid_case)call hf_to_huvw ! commente le 05-08-15

      end subroutine hz_to_hxyr

!........................................................

      subroutine hf_to_huvw !25-11-14
      use module_principal
      use module_parallele
      implicit none

! Pour eviter division par zero dans routine sigma_levels:
      do j=1,jmax+1 ; do i=1,imax+1
       if(h_f(i,j)>=0.) then
        h_f(i,j)=max(h_f(i,j), 1.d-3)
       else
        h_f(i,j)=min(h_f(i,j),-1.d-3)
       endif
      enddo ; enddo
! obc and mpi continuity at f points:
      call obc_h_f

       do j=1,jmax   ; do i=1,imax+1
         h_u(i,j)= 0.5*(h_f(i,j  )+h_f(i  ,j+1))
       enddo ; enddo


       do j=1,jmax+1 ; do i=1,imax
         h_v(i,j)= 0.5*(h_f(i,j  )+h_f(i+1,j  ))
       enddo ; enddo

       do j=1,jmax ; do i=1,imax
         h_w(i,j)=0.25*(h_f(i,j  )+h_f(i+1,j  )   &
                       +h_f(i,j+1)+h_f(i+1,j+1))
       enddo ; enddo

       call obc_h(0)

! Blindage ssh pour eviter division par zero dans module_ogcm
!     do j=0,jmax+1                                                   !06-06-09
!     do i=0,imax+1
!      ssh_int_w(i,j,:)=max(ssh_int_w(i,j,:),0.001*wetdry_cst3-h_w(i,j))
!     enddo
!     enddo

!     do j=0,jmax+1                                                    !30-11-14
!     do i=1,imax+1
!      ssh_int_u(i,j,:)=max(ssh_int_u(i,j,:),0.001*wetdry_cst3-h_u(i,j))
!     enddo
!     enddo

!     do j=1,jmax+1                                                    !30-11-14
!     do i=0,imax+1
!      ssh_int_v(i,j,:)=max(ssh_int_v(i,j,:),0.001*wetdry_cst3-h_v(i,j))
!     enddo
!     enddo

      end subroutine hf_to_huvw
