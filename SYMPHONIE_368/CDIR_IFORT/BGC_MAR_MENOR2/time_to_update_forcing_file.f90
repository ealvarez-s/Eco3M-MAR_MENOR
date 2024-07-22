










      subroutine time_to_update_forcing_file(t0_loc,period_loc)
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 14-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      implicit none
      double precision t0_loc,period_loc

!...............................................................................
! Version date      Description des modifications
! 2010.19 15-04-11  Date de mise en service
!...............................................................................

! decision=1 if it is time to update the considered forcing file:

      if(  int((elapsedtime_now-t0_loc)/period_loc)           &
         /=int((elapsedtime_bef-t0_loc)/period_loc) ) then

       decision=1

      else

       decision=0

      endif

      end subroutine time_to_update_forcing_file
