










      subroutine botprescont(txt_loc)
!______________________________________________________________________
! S model
! release S.26 - last update: 10-09-16 
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
      double precision mean_ssh_loc,lagrange_coef_loc,ddz_loc    &
       ,rhp_avr_loc,restore_loc
      character*10 txt_loc


      stop 'botprescont: your choice is not valid'

      end subroutine botprescont
