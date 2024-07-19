      subroutine biology
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
#ifdef synopsis
       subroutinetitle='biology'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!......................................................................
! Version date     Description des modifications
! 2010.8  03-05-10 terminologie
!......................................................................


      do vb=1,vbmax
      do  k=1,kmax
      do  j=1,jmax
      do  i=1,imax
      tendancebio_t(i,j,k,vb)=0.
      enddo
      enddo
      enddo
      enddo

      end subroutine biology
