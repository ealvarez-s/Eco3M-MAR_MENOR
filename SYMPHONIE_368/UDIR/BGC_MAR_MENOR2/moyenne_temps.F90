      subroutine moyenne_temps(ichoix)
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
       subroutinetitle='moyenne_temps'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!..............................................................................
! modifs 29/07/03: mis en service
!        25/04/07: passage coordonnee curviligne
!..............................................................................
      if(mi_onoff*grh_out_mi.eq.0)stop'erreur 1 dans moyenne_temps.f'

!******************************************************************************
! SOMMATION DES TABLEAUX DE VITESSE
! DEBUT:
      if(ichoix.eq.1)then
!******************************************************************************
      do k=1,kmax
       do j=1,jmax
        do i=1,imax
         vmi_u(i,j,k)=vmi_u(i,j,k)+veldydz_u(i,j,k,1)/dy_u(i,j)
         vmi_v(i,j,k)=vmi_v(i,j,k)+veldxdz_v(i,j,k,1)/dx_v(i,j)
        enddo
       enddo
      enddo
      sum_mi=sum_mi+1.
!******************************************************************************
! SOMMATION DES TABLEAUX DE VITESSE
! FIN.
      return
      endif
!******************************************************************************

!                      /     /    /

!******************************************************************************
! RESET DES TABLEAUX
! DEBUT:
      if(ichoix.eq.0)then
!******************************************************************************
      do k=1,kmax
       do j=1,jmax
        do i=1,imax
         vmi_u(i,j,k)=0.
         vmi_v(i,j,k)=0.
        enddo
       enddo
      enddo
      sum_mi=0.
!******************************************************************************
! RESET DES TABLEAUX
! FIN.
      return
      endif
!******************************************************************************


      return
      end
