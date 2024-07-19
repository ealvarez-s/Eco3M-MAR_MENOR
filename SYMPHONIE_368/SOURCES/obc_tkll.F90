      subroutine obc_tkll(ichoix)
!______________________________________________________________________
!
! S model
! release 2010.10  - last update: 24-06-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='obc_tkll'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! OUTPUT: H_X
! OUTPUT: H_Y
! OUTPUT: H_R

!...............................................................................
! modifs: 09/03/09: mis en service
!...............................................................................
!............................................!
! CONDITIONS LIMITES:
      do k=1,kmax+1
         do i=1,imax
          tkll_w(i,jmax,k)=tkll_w(i,jmax-1,k)
          tkle_w(i,jmax,k)=tkle_w(i,jmax-1,k)
          tkll_w(i,1   ,k)=tkll_w(i,2     ,k)
          tkle_w(i,1   ,k)=tkle_w(i,2     ,k)
         enddo
         do j=1,jmax
          tkll_w(imax,j,k)=tkll_w(imax-1,j,k)
          tkle_w(imax,j,k)=tkle_w(imax-1,j,k)
          tkll_w(1   ,j,k)=tkll_w(2     ,j,k)
          tkle_w(1   ,j,k)=tkle_w(2     ,j,k)
         enddo
      enddo

! cas du model_ 1D (I1D=0)
      if(flag3d.eq.0) then
      do k=1,kmax+1
      do j=1,jmax+1
      do i=1,imax+1
      tkll_w(i,j,k)=tkll_w(2,2,k)
      tkle_w(i,j,k)=tkle_w(2,2,k)
      enddo
      enddo
      enddo
! cas du model_ 1D: fin.
      end if

!............................................!


      return
      end
