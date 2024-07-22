










      subroutine latlon_to_ij(txt_)
!______________________________________________________________________
! S model
! release S.26  - last update: 10-09-13
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      use module_grid
      implicit none
      character*3 txt_

!__________________________________________________________________________
! Version date      Description des modifications
!         18/04/06: fonctions compatibles avec double precision
!         07/05/09: cas grille dlat dlon constant
! 2010.3  08-01-10: le pole de la grille n'est pas forcement le pole nord
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.13 21-10-10  choix entre indices globaux ou locaux
! S.26    02-06-13  northpole_lon northpole_lat renommes
!         05-06-13  cas de la grille bipole
!         10-09-13  reorientation vers nouvelle subroutine latlontoij
!__________________________________________________________________________

      call latlontoij(longi1,latit1,txt_) !10-09-13

      end subroutine latlon_to_ij
