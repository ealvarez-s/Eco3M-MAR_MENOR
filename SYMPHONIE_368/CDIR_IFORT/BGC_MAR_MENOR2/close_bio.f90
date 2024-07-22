










      subroutine close_bio
!______________________________________________________________________
!
! S model
! release 2010.22  - last update: 27-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none


!...............................................................................
! Version Date      Description des modifications
!         25/12/02: amenagement pour filiere restart; CALL HOT_RESTART(22)
! 2010.10 23-06-10  hot_restart renommé dyn_restart et bio_restart
! 2010.22 27-04-11  ichanel9 remplacé par restartfileperiod
!...............................................................................

!*******************************************************************************
! INITIALISATION DES TABLEAUX:
! DEBUT:
!*******************************************************************************

! Archivage de l'état final de la simulation:

      if(restartfileperiod>0.)call bio_restart('w')               !27-04-11

!*******************************************************************************
! INITIALISATION DES TABLEAUX:
! FIN.
!*******************************************************************************

      end subroutine close_bio
