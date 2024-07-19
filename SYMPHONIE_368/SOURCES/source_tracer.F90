      subroutine source_tracer
!______________________________________________________________________
!
! S model
! release 2010.25  - last update: 27-02-12
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='source_tracer'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version Date      Description des modifications
!         19/03/03: mise en service
! 2010.8  03-05-10  Terminologie
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes
! 2010.25 27-02-12  Pas de calcul si hors domaine
!...............................................................................


!***************************************************************************
! Forçage des sources des tests de dispersion:
! DEBUT
!***************************************************************************

      do 4 kso=1,ksomax

       i=nint(socard(1,kso)) ; j=nint(socard(2,kso))
       if(i>=1.and.i<=imax.and.j>=1.and.j<=jmax) then!2222222>

        k1=nint(socard(3,kso))
        k2=nint(socard(4,kso))

        x1=0.
        if(elapsedtime_now>=socard(8,kso).and.   &
           elapsedtime_now<=socard(9,kso))x1=1.

         do vb=nint(socard(5,kso)),nint(socard(6,kso))
         do k=k1,k2
          tendancebio_t(i,j,k,vb)=socard(7,kso)*x1
         enddo
         enddo

       endif                                         !2222222>

    4 continue

!***************************************************************************
! Forçage des sources des tests de dispersion:
! FIN
!***************************************************************************

      end subroutine source_tracer
