      subroutine obc_tken(ichoix)
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
      use module_parallele !#MPI
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='obc_tken'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         09/03/09: mis en service
!         05-06-09: Parallelisation
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
!...............................................................................


!............................................!
! CONDITIONS LIMITES:
      do k=1,kmax+1
         do i=1,imax
          tken_w(i,jmax,k)=tken_w(i,jmax-1,k)
          tken_w(i,1   ,k)=tken_w(i,2     ,k)
         enddo
         do j=1,jmax
          tken_w(imax,j,k)=tken_w(imax-1,j,k)
          tken_w(1   ,j,k)=tken_w(2     ,j,k)
         enddo
      enddo

! cas du model_ 1D (I1D=0)
      if(flag3d.eq.0) then !---------1D
       do k=1,kmax+1
        do j=1,jmax+1
         do i=1,imax+1
          tken_w(i,j,k)=tken_w(2,2,k)
         enddo
        enddo
       enddo
! cas du model_ 1D: fin.
      end if            !---------1D
!............................................!

#ifdef parallele
      call echange('z0',tken_w,lbound(tken_w),ubound(tken_w)) !#MPI i=imax i=1 j=jmax  j=1
#endif

      end subroutine obc_tken
