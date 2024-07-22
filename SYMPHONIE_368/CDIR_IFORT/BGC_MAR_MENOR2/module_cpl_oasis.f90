










   module module_cpl_oasis
!______________________________________________________________________
! SYMPHONIE ocean model
! release 290  - last update: 02-10-20
!______________________________________________________________________
!...............................................................................
!    _________                    .__                  .__                     !m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................
!------------------------------------------------------------------------------!
! version   date    description                                                !
! S26       25-06-15  mise en service (Léo Seyfried)                           !
!           06-03-17  mise à jour pour la version 6.0 de WW3 + mise en forme   !
!                      (Fabien Rétif)                                          !
! v290      02-10-20 cas sym_sym
!------------------------------------------------------------------------------!

   implicit none

   integer :: istr_oasis,iend_oasis,jstr_oasis,jend_oasis  
   integer :: il_part_id  &   ! PartitionID
             ,var_nodims(2),var_actual_shape(1),var_type,ierror &
             ,oasis_get_count=0
   integer :: var_oasis_id(2)
   integer :: group_send_dim=0 , group_recv_dim=0 , oasismisspointmax=0 , kmax_recv 
   double precision, dimension(:,:,:), allocatable :: group_send_ocean, group_recv_ocean
   integer         , dimension(:,:)  , allocatable :: grid_msk_ocean &
                                                     ,vert_extrap_i  &
                                                     ,vert_extrap_j  &
                                                     ,oasisplugs
   double precision      :: oasis_time_next
   integer               :: il_compid                 ! Component model ID returned by oasis_init_comp
   character(len=6)      :: cl_model_name             ! Model name (same as in namcouple)   
   character(len=60)     :: oasis_dir
   integer               :: cl_model_num              ! Numero dans la liste des modeles couples
   integer               :: il_err                    ! Return error code


!#ifdef key_oasis
!  !! * Accessibility
!  PUBLIC cpl_oasis_init
!  PUBLIC cpl_oasis_grid
!  PUBLIC cpl_oasis_define 
!  PUBLIC cpl_oasis_snd
!  PUBLIC cpl_oasis_rcv
!  PUBLIC cpl_oasis_finalize   
!#endif

!pmx routines:
!  PUBLIC cpl_oasis_def_var_pmx
!  PUBLIC cpl_oasis_put_pmx
!  PUBLIC cpl_oasis_get_pmx
!  PUBLIC cpl_oasis_def_part_pmx
!  PUBLIC cpl_oasis_grids_pmx
!  PUBLIC cpl_oasis_debug_pmx

contains


!!======================================================================


!!======================================================================


!!======================================================================


!!======================================================================


!!======================================================================


!...............................................................
  subroutine cpl_oasis_def_part_pmx
  end subroutine cpl_oasis_def_part_pmx
!...............................................................
  subroutine cpl_oasis_grids_pmx
  end subroutine cpl_oasis_grids_pmx

!...............................................................

  subroutine cpl_oasis_def_var_pmx
  end subroutine cpl_oasis_def_var_pmx

!...............................................................

  subroutine cpl_oasis_get_pmx
  end subroutine cpl_oasis_get_pmx

!...............................................................

  subroutine cpl_oasis_put_pmx
  end subroutine cpl_oasis_put_pmx

!------------------------------------------------------------------------------

      subroutine cpl_oasis_def_oasisplugs

      end subroutine cpl_oasis_def_oasisplugs
!..................................................................
      subroutine cpl_oasis_debug_pmx !17-08-20
      end subroutine cpl_oasis_debug_pmx
!..................................................................
      subroutine cpl_oasis_reset_pmx !17-08-20
      end subroutine cpl_oasis_reset_pmx !17-08-20

!..................................................................

      subroutine cpl_oasis_mpi_delta
      end subroutine cpl_oasis_mpi_delta

!..................................................................

      subroutine cpl_oasis_tem_delta
      end subroutine cpl_oasis_tem_delta

!..................................................................

      subroutine cpl_oasis_sal_delta
      end subroutine cpl_oasis_sal_delta

!..................................................................

  subroutine cpl_oasis_init_pmx(id_lcomm)
  use module_principal ! , only : tmpdirname,texte4,k ! ,istr_oasis,iend_oasis,jstr_oasis,jend_oasis
  implicit none
! integer, INTENT(OUT) :: id_lcomm                   ! Model local communicator
  integer              :: id_lcomm                   ! Model local communicator
  end subroutine cpl_oasis_init_pmx

!..................................................................

  subroutine cpl_oasis_finalize
  implicit none


  end subroutine cpl_oasis_finalize

!..................................................................

end module module_cpl_oasis
