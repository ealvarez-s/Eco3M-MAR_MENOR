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

#ifdef key_oasis_sym_sym
   use mod_oasis                                      ! OASIS3-MCT module 
#endif
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

#ifdef key_oasis
   integer, public       :: il_nb_rcv, il_nb_snd      ! Number of coupling fields    
   integer, parameter    :: ip_maxfld=50              ! Maximum number of coupling fields
   
   type, public          :: CPL_FIELD                 ! Type for coupling field information
      character(len = 8) :: cl_field_name             ! Name of the coupling field   
      integer            :: il_field_id               ! Field ID
   end type CPL_FIELD
   type(CPL_FIELD), dimension(ip_maxfld), public :: rcv_fld, snd_fld   ! Coupling fields
#endif

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

#ifdef key_oasis
  subroutine cpl_oasis_init(id_lcomm)
  !==
  ! Initialisation de OASIS
  !==
   integer, INTENT(OUT) :: id_lcomm                   ! Model local communicator

!   PRINT_DBG*, 'enter CPL_OASIS_INIT routine'
   
   !! Initialise the coupling
   CALL oasis_init_comp (il_compid, cl_model_name, il_err)
   IF (il_err /= 0) THEN
      CALL oasis_abort(il_compid, 'cpl_oasis_init', 'Problem during oasis_init_comp')
   ENDIF
   
   !! Get the value of a local MPI communicator to be used by SYMPHONIE for its internal parallelisation
   CALL oasis_get_localcomm (id_lcomm, il_err)
   IF (il_err /= 0) THEN
      CALL oasis_abort(il_compid, 'cpl_oasis_init', 'Problem during oasis_get_localcomm')
   ENDIF

!   PRINT_DBG*, 'exit CPL_OASIS_INIT routine'

  end subroutine cpl_oasis_init
#endif

!!======================================================================

#ifdef key_oasis
  subroutine cpl_oasis_grid(ld_master,id_lcomm)
  !==
  ! Définition de la grille étendue avec les conditions aux limites
  ! Note: routine non testée jusqu'à présent
  !==

        USE module_principal
        USE module_parameter
        USE module_parallele

        LOGICAL, INTENT(IN) :: ld_master                  ! MASTER process or not
        INTEGER, INTENT(IN) :: id_lcomm                   ! Model local communicator

        if(.not.allocated(globlon_t))allocate(globlon_t(iglb  ,jglb  ))
        if(.not.allocated(globlat_t))allocate(globlat_t(iglb  ,jglb  ))

        do j=1,jmax ; do i=1,imax
        globlon_t(i+par%timax(1),j+par%tjmax(1))=lon_t(i,j)
        globlat_t(i+par%timax(1),j+par%tjmax(1))=lat_t(i,j)
        enddo           ; enddo
        ! Assemblage global:
        ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
        ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
        call par_gatherall_2d(globlon_t,lbound(globlon_t),ubound(globlon_t),ind,par%nbdom)
        call par_gatherall_2d(globlat_t,lbound(globlat_t),ubound(globlat_t),ind,par%nbdom)

        !oasis routine
        CALL oasis_write_grid ('symt', iglb, jglb,                &
            &                globlon_t(:,:),  &
            &                globlat_t(:,:)   &
            &                )

        deallocate(globlon_t)
        deallocate(globlat_t)

  end subroutine cpl_oasis_grid
#endif

!!======================================================================

#ifdef key_oasis
  subroutine cpl_oasis_define
  !==
  ! Définition du partitionnement
  !== 
   use module_parameter, only :  imax,jmax &
                                 ,iglb,jglb
   use module_parallele                             
   implicit none
   
   integer                 :: ib_i, rank
   integer                 :: il_part_id      ! PartitionID
   integer, dimension(5)   :: ila_paral       ! Description of the local partition in the global index space
   integer, dimension(2,2) :: ila_shape       ! Vector giving the min & max index for each dim of the fields
   integer, dimension(2)   :: ila_var_nodims  ! rank of fields & number of bundles (1 with OASIS3-MCT)
   logical                 :: ll_mpi_file     ! to check if there an mpi.txt file for domain decompasition

!   PRINT_DBG*, 'enter CPL_OASIS_DEFINE routine'

   !! Définition des partitions selon la méthode OASIS BOX. 
   !! OASIS utilise une définition de sa "grille" selon un vecteur 1D sans
   !! recoupement des processeurs. Ici chaque processeur déclare sa zone
   !! effective allant de 2 à imax-1 
   !! Extrait de la doc OASIS :
   !! ig_paral(1) = 2 (indicates a Box partition)
   !! ig_paral(2) = the upper left corner global offset
   !! ig_paral(3) = the local extent in x
   !! ig_paral(4) = the local extent in y
   !! ig_paral(5)= the global extent in x.

   !! Define the partition : OASIS BOX partition
   ila_paral(1) = 2
   ila_paral(2) = (iglb+2)*(par%tjmax(1)+2)+(par%timax(1)+2)
   ila_paral(3) = imax-2
   ila_paral(4) = jmax-2
   ila_paral(5) = iglb+2

   call oasis_def_partition (il_part_id, ila_paral, il_err, (iglb+2)*(jglb+2))
 
   if (il_err /= 0) then
      call oasis_abort(il_compid, 'cpl_oasis_define', 'Problem during oasis_def_partition')
   endif

   ila_shape(:,1) = (/1, imax-2 /)
   ila_shape(:,2) = (/1, jmax-2 /)

   !! Coupling fields declaration
   ila_var_nodims(1) = 2    ! rank of fields array
   ila_var_nodims(2) = 1    ! always 1 with OASIS3-MCT 2.0
      
   call get_list_exch_field (rcv_fld, snd_fld, il_nb_rcv, il_nb_snd)

   !! Send coupling fields  
   do ib_i = 1, il_nb_snd
      call oasis_def_var (snd_fld(ib_i)%il_field_id     &
           &            , snd_fld(ib_i)%cl_field_name   &
           &            , il_part_id                    & 
           &            , ila_var_nodims                &
           &            , OASIS_Out                     &
           &            , ila_shape                     &
           &            , OASIS_Real                    &
           &            , il_err )
      
      if (il_err /= 0) then
         call oasis_abort(il_compid, 'cpl_oasis_define', 'Problem during oasis_def_var')
      endif
   enddo

   !! Received coupling fields  
   do ib_i = 1, il_nb_rcv
      call oasis_def_var (rcv_fld(ib_i)%il_field_id    &
           &            , rcv_fld(ib_i)%cl_field_name  &
           &            , il_part_id                   & 
           &            , ila_var_nodims               &
           &            , OASIS_In                     &
           &            , ila_shape                    &
           &            , OASIS_Real                   &
           &            , il_err )
      
      if (il_err /= 0) then
         CALL oasis_abort(il_compid, 'cpl_oasis_define', 'Problem during oasis_def_var')
      endif         
   enddo
   
   !! End of definition phase
   call oasis_enddef(il_err)
   if (il_err /= 0) then
      CALL oasis_abort(il_compid, 'cpl_oasis_define', 'Problem during oasis_enddef')
   endif

!   PRINT_DBG*, 'exit CPL_OASIS_DEFINE routine'

  end subroutine cpl_oasis_define
#endif

!!======================================================================

#ifdef key_oasis
  subroutine cpl_oasis_snd(id_nb, id_time, rda_field, ld_action)
   implicit none
   !! Arguments
   INTEGER, INTENT(IN)   :: id_nb                         ! Number of the field to be send
   INTEGER, INTENT(IN)   :: id_time                       ! Atmosphere time-step in seconds
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rda_field  ! Coupling field array to be send
   LOGICAL, INTENT(OUT)  :: ld_action                     ! Action performed
   
   !! * Local declarations
   INTEGER :: il_info                                     ! OASIS3-MCT info argument

!   PRINT_DBG*, 'enter CPL_OASIS_SND routine'

   CALL oasis_put ( snd_fld(id_nb)%il_field_id &
   &              , id_time                    &
   &              , rda_field                  &
   &              , il_info                    &
   &                )

   ld_action = il_info == OASIS_Sent     .OR. il_info == OASIS_ToRest .OR.   &
   &           il_info == OASIS_SentOut  .OR. il_info == OASIS_ToRestOut

!   PRINT_DBG*, 'exit CPL_OASIS_SND routine'
  end subroutine cpl_oasis_snd    
#endif

!!======================================================================

#ifdef key_oasis
  subroutine cpl_oasis_rcv(id_nb, id_time, rda_field, ld_action)
   use module_parameter, only :  imax,jmax
   implicit none
   !! * Argument
   INTEGER, INTENT(IN)   :: id_nb                          ! Number of the field to be received
   INTEGER, INTENT(IN)   :: id_time                        ! Ocean time-step in seconds
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: rda_field   ! Coupling field array to be received
   LOGICAL, INTENT(OUT)  :: ld_action                      ! Action performed
!   REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE :: field
   !! * Local declarations
   INTEGER :: il_info,i,j                                      ! OASIS3-MCT info argument
  
!   PRINT_DBG*, 'enter CPL_OASIS_RCV routine'

!   allocate(field(imax-2,jmax-2))
   call oasis_get ( rcv_fld(id_nb)%il_field_id &
   &              , id_time                    &
   &              , rda_field                  &
   &              , il_info                    &
   &                )
   
   ld_action = il_info == OASIS_Recvd   .OR. il_info == OASIS_FromRest .OR.   &
   &           il_info == OASIS_RecvOut .OR. il_info == OASIS_FromRestOut!
!   do j=2,jmax-1
!   do i=2,imax-1
!   rda_field(i,j)=field(i-1,j-1)
!   enddo 
!   enddo
!   deallocate(field)
!   PRINT_DBG*, 'exit CPL_OASIS_RCV routine'

  end subroutine cpl_oasis_rcv
#endif

!!======================================================================

#ifdef key_oasis
  subroutine get_list_exch_field(rcv, snd, id_nb_rcv, id_nb_snd)
   !!==
   !! Définition des champs à échanger
   !! 1. Selon les clès de compilation, on lit les notebooks (32 pour un couplage avec ATM ou
   !! 33 pour WAVE)
   !! 2. On itère sur les champs et on les déclare à OASIS
   !! Note: Les noms des champs dans OASIS sont limités à 8 caractères
   use module_principal, only : nomfichier
   implicit none
   TYPE(CPL_FIELD), DIMENSION(ip_maxfld), INTENT (INOUT) :: rcv, snd
   INTEGER, INTENT(INOUT) :: id_nb_rcv, id_nb_snd 
   logical :: l_cpl_ocn_sst,l_cpl_ocn_maskt,l_cpl_ocn_masku,l_cpl_ocn_maskv &
                        ,l_cpl_ocn_ssh,l_cpl_ocn_cur,l_cpl_ocn_ucur,l_cpl_ocn_vcur &
                        ,l_cpl_atm_rain,l_cpl_atm_evap,l_cpl_atm_wat,l_cpl_atm_rad,l_cpl_atm_heat &
                        ,l_cpl_atm_ir,l_cpl_atm_lat,l_cpl_atm_sens,l_cpl_atm_pa &
                        ,l_cpl_atm_fwsu,l_cpl_atm_fwsv,l_cpl_ww3_dir     &
                        ,l_cpl_ww3_hs,l_cpl_ww3_hsw,l_cpl_ww3_foc,l_cpl_ww3_tw &
                        ,l_cpl_ww3_taw,l_cpl_ww3_two&
                        ,l_cpl_ww3_uss,l_cpl_ww3_mask
     
!   PRINT_DBG*, 'enter GET_LIST_EXCH_FIELD routine'
   ! Tous les champs sont initialisés à Faux
   l_cpl_ocn_sst=.FALSE. ; l_cpl_ocn_maskt=.FALSE. ; l_cpl_ocn_masku=.FALSE.
   l_cpl_ocn_maskv=.FALSE. ; l_cpl_ocn_ssh=.FALSE. ; l_cpl_ocn_cur=.FALSE. 
   l_cpl_atm_rain=.FALSE. ; l_cpl_atm_evap=.FALSE.
   l_cpl_atm_wat=.FALSE. ; l_cpl_atm_rad=.FALSE. ; l_cpl_atm_heat=.FALSE.
   l_cpl_atm_ir=.FALSE. ; l_cpl_atm_lat=.FALSE. ; l_cpl_atm_sens=.FALSE.
   l_cpl_atm_pa=.FALSE. ; l_cpl_atm_fwsu=.FALSE. ; l_cpl_atm_fwsv=.FALSE.
   l_cpl_ww3_dir=.FALSE.;l_cpl_ww3_hs=.FALSE.;l_cpl_ww3_hsw=.FALSE.
   l_cpl_ww3_foc=.FALSE.;l_cpl_ww3_tw=.FALSE.;l_cpl_ww3_taw=.FALSE.;
   l_cpl_ww3_two=.FALSE.; l_cpl_ww3_uss=.FALSE.; l_cpl_ww3_mask=.FALSE.

#ifdef key_oasis_symp_surfex
   namelist/notebook_cpl_surfex/l_cpl_ocn_sst,l_cpl_ocn_maskt,l_cpl_ocn_masku,l_cpl_ocn_maskv &
                        ,l_cpl_ocn_ssh,l_cpl_ocn_ucur,l_cpl_ocn_vcur &
                        ,l_cpl_atm_rain,l_cpl_atm_evap,l_cpl_atm_wat &
                        ,l_cpl_atm_rad,l_cpl_atm_heat &
                        ,l_cpl_atm_ir,l_cpl_atm_lat,l_cpl_atm_sens,l_cpl_atm_pa &
                        ,l_cpl_atm_fwsu,l_cpl_atm_fwsv
       open(100,file=nomfichier(32))
       read(100,nml=notebook_cpl_surfex)
       close(100)
#else 
#endif

#ifdef key_oasis_symp_ww3
   namelist/notebook_cpl_ww3/l_cpl_ocn_sst,l_cpl_ocn_ssh,l_cpl_ocn_cur &
                        ,l_cpl_ww3_dir,l_cpl_ww3_hs,l_cpl_ww3_hsw &
                        ,l_cpl_ww3_foc,l_cpl_ww3_tw &
                        ,l_cpl_ww3_taw,l_cpl_ww3_two &
                        ,l_cpl_ww3_uss,l_cpl_ww3_mask
       open(100,file=nomfichier(33))
       read(100,nml=notebook_cpl_ww3)
       close(100)
#else 
#endif

   !! Coupling fields send by S
   id_nb_snd = 0

   ! MARS Mask
   IF (l_cpl_ocn_maskt) THEN
     id_nb_snd=id_nb_snd+1
     snd(id_nb_snd)%cl_field_name='SYMP_MSK'
   END IF
   
   ! temp : sea surface temperature (Degrees C)
   IF (l_cpl_ocn_sst) THEN
     id_nb_snd=id_nb_snd+1
     snd(id_nb_snd)%cl_field_name='SYMP_SST'
   END IF
 
   ! ssh : sea surface height (m)
   IF (l_cpl_ocn_ssh) THEN
     id_nb_snd=id_nb_snd+1
     snd(id_nb_snd)%cl_field_name='SYMP_SSH'
   END IF
 
   ! sea surface currents (m.s-1)
   IF (l_cpl_ocn_cur) THEN
     ! uz : sea surface zonal currents (m.s-1)
     id_nb_snd=id_nb_snd+1
     snd(id_nb_snd)%cl_field_name='SYMP__UZ'
     
     ! vz : sea surface meridional currents (m.s-1)
     id_nb_snd=id_nb_snd+1
     snd(id_nb_snd)%cl_field_name='SYMP__VZ'
   END IF

   ! Atmospheric variables
   ! ---------------------

   !! Coupling fields received by S
   id_nb_rcv = 0   

   ! RAIN (mm.s-1)
   IF (l_cpl_atm_rain) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPRAIN'
   END IF

   ! Evaporation (mm.s-1)
   IF (l_cpl_atm_evap) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPEVAP'
   END IF
   
   ! Water flux
   IF (l_cpl_atm_wat) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_WAT'
   END IF

   ! Solar Flux (W.m-2)
   IF (l_cpl_atm_rad) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_QSR'
   END IF

   ! Total heat flux (LW + SW + Qlat + Qsens) (W.m-2)
   IF (l_cpl_atm_heat) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPHEAT'
   END IF

   ! Long-wave thermal flux (W.m-2)
   IF (l_cpl_atm_ir) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP__IR'
   END IF

   ! Latent heat flux (W.m-2)
   IF (l_cpl_atm_lat) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_LAT'
   END IF

   ! Sensible heat flux (W.m-2)
   IF (l_cpl_atm_sens) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_SEN'
   END IF
 
   ! Surface Pressure (Pa)
   IF (l_cpl_atm_pa) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP__PA'
   END IF

   ! Zonal Wind Stress (N.m-2)
   IF (l_cpl_atm_fwsu) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPTAUX'
   END IF

   ! Meridional Wind Stress (N.m-2)
   IF (l_cpl_atm_fwsv) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPTAUY' 
   END IF
   
   ! Wave variables
   ! ---------------------
   IF (l_cpl_ww3_dir) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPSDIR'
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPCDIR'
   END IF

   IF (l_cpl_ww3_hs) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP__HS' 
   END IF

   IF (l_cpl_ww3_hsw) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_HSW' 
   END IF

   IF (l_cpl_ww3_foc) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_FOC' 
   END IF

   IF (l_cpl_ww3_tw) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP__TW' 
   END IF

   IF (l_cpl_ww3_taw) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPTAWX' 
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPTAWY' 
   END IF

   IF (l_cpl_ww3_two) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPTWOX' 
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPTWOY' 
   END IF

   IF (l_cpl_ww3_uss) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPUSSX' 
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMPUSSY' 
   END IF

   IF (l_cpl_ww3_mask) THEN
     id_nb_rcv = id_nb_rcv +1
     rcv(id_nb_rcv)%cl_field_name='SYMP_MSK' 
   END IF

  end subroutine get_list_exch_field
#endif

!...............................................................
  subroutine cpl_oasis_def_part_pmx
#ifdef key_oasis_sym_sym
  use module_parallele ; use module_parameter ! ; use module_s
  use module_principal ! only : tmpdirname,texte4,k ! ,istr_oasis,iend_oasis,jstr_oasis,jend_oasis
  use module_s
  implicit none


! INTEGER               :: var_id_oasis(2)
! INTEGER               :: var_nodims(2)
! INTEGER               :: var_actual_shape(1)
! INTEGER               :: var_type
  INTEGER               :: ierror
! INTEGER               :: il_part_id
  INTEGER               :: ig_paral_size
  INTEGER, DIMENSION(:), ALLOCATABLE :: ig_paral,il_size_onerank,il_size_allrank 
!                                      ,kmax_reduce_loc,kmax_reduce_sum

  integer :: il_extentx, il_extenty, il_size, il_offsetx, il_offsety, il_offset 

!...................
! VERIFIER kmax_recv
! Partager kmax avec les autres modeles en employant une fonction bcast sur l'univers vu par les 2 codes (MPI_COMM_WORLD). 
! Le but est de verifier la valeur kmax_recv indiquee dans notebook_oasis
      call MPI_COMM_SIZE(MPI_COMM_WORLD,i1,ierror)   
      k1=kmax ; k2=kmax
      call mpi_bcast(k1,1,mpi_integer,   0,MPI_COMM_WORLD,ierror)
      call mpi_bcast(k2,1,mpi_integer,i1-1,MPI_COMM_WORLD,ierror)
! Dans MPI_COMM_WORLD le rank 0    est le rank 0       du modele 1
! Dans MPI_COMM_WORLD le rank i1-1 est le dernier rank du modele 2
! k1 doit donc etre egale A kmax_recv du modele 2
! k2 doit donc etre egale A kmax_recv du modele 1
! Note: cet algo ne marche pas si il y a plus que 2 modeles.
      if(cl_model_num==1.and.kmax_recv/=k2) then !m°v°m>
       write(10+par%rank,*)'Dans notebook_oasis du modele 1' &
       ,' kmax_recv ne correspond pas A kmax du modele 2'
       stop 'Err 605 kmax_recv see fort.xxx files'
      endif                                      !m°v°m>
      if(cl_model_num==2.and.kmax_recv/=k1) then !m°v°m>
       write(10+par%rank+i1,*)'Dans notebook_oasis du modele 2' &
       ,' kmax_recv ne correspond pas A kmax du modele 1'
       stop 'Err 612 kmax_recv see fort.xxx files'
      endif                                      !m°v°m>
!...................

#ifdef bidon
  istr_oasis=2
  iend_oasis=imax-1
  jstr_oasis=2
  jend_oasis=jmax-1
  if(1   +par%timax(1)==1)   istr_oasis=1
  if(imax+par%timax(1)==iglb)iend_oasis=imax
  if(1   +par%tjmax(1)==1)   jstr_oasis=1
  if(jmax+par%tjmax(1)==jglb)jend_oasis=jmax

  il_extentx=iend_oasis-istr_oasis+1
  il_extenty=jend_oasis-jstr_oasis+1
  il_size=il_extentx*il_extenty
  il_offsetx=istr_oasis-1+par%timax(1)
  il_offsety=jstr_oasis-1+par%tjmax(1)
! il_offset=(jstr_oasis-1+par%tjmax(1))*il_extentx ! Cas particulier des tranches de saucissons
  il_offset=(jstr_oasis-1+par%tjmax(1))*iglb &
           +(istr_oasis-1+par%timax(1))*il_extenty

!.........................
! il_offset obtenu par le cumul de il_size des rank precedents:
! allocate(il_size_onerank(0:nbdom-1)) ; il_size_onerank=0 
! allocate(il_size_allrank(0:nbdom-1)) ; il_size_allrank=0
! il_size_onerank(par%rank)=il_size
! call mpi_allreduce(il_size_onerank,il_size_allrank,nbdom,mpi_integer,mpi_sum,par%comm2d,ierr)
! il_offset=0
! do k=0,par%rank-1
!  il_offset=il_offset+il_size_allrank(k)
! enddo
! deallocate(il_size_onerank)
! deallocate(il_size_allrank)
!.........................
  allocate(ig_paral(3))
  ig_paral(1)=1
  ig_paral(2)=il_offset
  ig_paral(3)=il_size
  CALL oasis_def_partition (il_part_id, ig_paral, ierror)
#endif
  
  
!! Define the partition : OASIS BOX partition
  istr_oasis=2
  iend_oasis=imax-1
  jstr_oasis=2
  jend_oasis=jmax-1
  if(1   +par%timax(1)==1)   istr_oasis=0
  if(imax+par%timax(1)==iglb)iend_oasis=imax+1
  if(1   +par%tjmax(1)==1)   jstr_oasis=0
  if(jmax+par%tjmax(1)==jglb)jend_oasis=jmax+1

  allocate(ig_paral(5))
  ig_paral(1) = 2
  ig_paral(3) = iend_oasis-istr_oasis+1
  ig_paral(4) = jend_oasis-jstr_oasis+1

! Cas grille definie de 0 a iglb+1 et de 0 a jglb+1
  if(namsrc_nx(cl_model_num)/=iglb+2.or.namsrc_ny(cl_model_num)/=jglb+2) &
  stop 'err dimensions iglb+2,jglb+2 ne correspondent pas a la namcouple'
! ig_paral(2) = (iglb+2)*(par%tjmax(1)+2         )+(par%timax(1)+2)
  ig_paral(2) = (iglb+2)*(par%tjmax(1)+jstr_oasis)+(par%timax(1)+istr_oasis)
  ig_paral(5) = iglb+2
  call oasis_def_partition (il_part_id, ig_paral, il_err, (iglb+2)*(jglb+2))

! Cas grille definie de 1 a iglb   et de 1 a jglb  
! Verifier les dimensions donnees dans la namcouple
! if(namsrc_nx(cl_model_num)/=iglb.or.namsrc_ny(cl_model_num)/=jglb) &
! stop 'err dimensions iglb,jglb ne correspondent pas a la namcouple'
! ig_paral(2) = iglb*(par%tjmax(1)+jstr_oasis-1)+(par%timax(1)+istr_oasis-1)
! ig_paral(5) = iglb
! call oasis_def_partition (il_part_id, ig_paral, ierror , iglb*jglb)

  k=s_unit(7)
  write(texte4,'(i0)')par%rank
  open(unit=k,file=trim(tmpdirname)//trim(cl_model_name)//'_oasis_'//trim(texte4),position='append')
  write(k,'(a)')'........'
  write(k,*)'subroutine cpl_oasis_def_part_pmx'
  write(k,'(a,a)')'cl_model_name: ',trim(cl_model_name)
  write(k,*)'par%rank',par%rank
  write(k,*)'il_part_id',il_part_id
  write(k,*)'istr_oasis,iend_oasis',istr_oasis,iend_oasis
  write(k,*)'jstr_oasis,jend_oasis',jstr_oasis,jend_oasis
  write(k,*)'il_extentx',il_extentx
  write(k,*)'il_extenty',il_extenty
  write(k,*)'il_size   ',il_size
  write(k,*)'il_offsetx',il_offsetx
  write(k,*)'il_offsety',il_offsety
  write(k,*)'il_offset ',il_offset
  write(k,*)'ig_paral  ',ig_paral
  write(k,*)'Coupling periods (s) ',namflddti
  write(k,*)'Coupling lags    (s) ',namfldlag
  write(k,*)'oasis_def_partition => ierror ',ierror
  write(k,'(a)')'........'
  close(k)
  
  deallocate(ig_paral)

! call mpi_barrier(par%comm2d,k)      ! synchro processes

#endif
  end subroutine cpl_oasis_def_part_pmx
!...............................................................
  subroutine cpl_oasis_grids_pmx
#ifdef key_oasis_sym_sym
  use module_parameter
  use module_principal ! only : flag,rad2deg,lon_t,lat_t,i,j,lon_f,lat_f,mask_t
  use module_parallele 
  use module_s
  implicit none
  double precision , dimension(:,:)   , allocatable :: grid_lon_ocean, grid_lat_ocean
  double precision , dimension(:,:,:) , allocatable :: grid_clo_ocean, grid_cla_ocean
! integer          , dimension(:,:)   , allocatable :: grid_msk_ocean
  integer ncid_
  logical :: exists               ! file exist flag
  include 'netcdf.inc'

  allocate(grid_lon_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis))
  allocate(grid_lat_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis))
  allocate(grid_clo_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis,4))
  allocate(grid_cla_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis,4))

  allocate(grid_msk_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis)) ; grid_msk_ocean=1
  allocate(vert_extrap_i(istr_oasis:iend_oasis,jstr_oasis:jend_oasis)) ; vert_extrap_i=-999
  allocate(vert_extrap_j(istr_oasis:iend_oasis,jstr_oasis:jend_oasis)) ; vert_extrap_j=-999

  texte90=trim(oasis_dir)//'oasis_grid_msk.nc'
  inquire(file=trim(texte90),exist=exists)
  if(exists) then !pmx>
! Soit le fichier de masque existe (et on le lit)
       status=nf_open(trim(texte90),nf_nowrite,ncid_)     ; if(status/=0) stop 'err 657 nf_open'
       status=nf_inq_dimid(ncid_,'ni_t',k0)               ; if(status/=0) stop 'err 658 nf_inq_dimid'
       status=nf_inq_dimlen(ncid_,k0,i1)                  ; if(status/=0) stop 'err 659 nf_inq_dimlen'
       status=nf_inq_dimid(ncid_,'nj_t',k0)               ; if(status/=0) stop 'err 660 nf_inq_dimid'
       status=nf_inq_dimlen(ncid_,k0,j1)                  ; if(status/=0) stop 'err 661 nf_inq_dimlen'
       if(i1/=iglb+2.or.j1/=jglb+2)stop 'Err 657 i1/=iglb+2.or.j1/=jglb+2'
       varstart(1)=par%timax(1)+istr_oasis+1 ; varcount(1)=iend_oasis-istr_oasis+1
       varstart(2)=par%tjmax(1)+jstr_oasis+1 ; varcount(2)=jend_oasis-jstr_oasis+1
       status=nf_inq_varid(ncid_,'oasis_grid_msk',var_id) ; if(status/=0) stop 'err 662 nf_inq_varid'
       status=nf_get_vara_int(ncid_,var_id,varstart(1:2),varcount(1:2),grid_msk_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis))
                                                            if(status/=0) stop 'err 663 nf_get_vara_int'
!.... lignes suivantes pour verification....
! Au besoin on peut lire les lon,lat du fichier pour verifier egalitE avec lon_t,lat_t (ne pas lire sinon)
!      status=nf_inq_varid(ncid_,'lon',var_id)            ; if(status/=0) stop 'err 664 nf_inq_varid'
!      status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),grid_lon_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis))
!                                                           if(status/=0) stop 'err 665 nf_get_vara_double'
!      status=nf_inq_varid(ncid_,'lat',var_id)            ; if(status/=0) stop 'err 666 nf_inq_varid'
!      status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),grid_lat_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis))
!                                                           if(status/=0) stop 'err 667 nf_get_vara_double'
!...
       status=nf_close(ncid_)                             ; if(status/=0) stop 'err 668 nf_close'
  else            !pmx>
! Soit le fichier de masque n'existe pas et celui-ci est deduit de mask_t
       do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
            grid_msk_ocean(i,j)=1-mask_t(i,j,kmax) ! conventions oasis: 0=sea, 1=land
       enddo                      ; enddo
  endif           !pmx>

  do j=jstr_oasis,jend_oasis
  do i=istr_oasis,iend_oasis

      grid_lon_ocean(i,j)=rad2deg*lon_t(i,j)  
      grid_lat_ocean(i,j)=rad2deg*lat_t(i,j)  

! patrick: grid_lon_ocean pourrait correspondre A lon_t
!       et grid_clo_ocean serait lon_w quelque soit la valeur de nc_ocean allant
!       de 1 A 4. La difference viendrait dans le decallage comme dans le schema
!       suivant:
! clo(i,j,2),cla(i,j,2)                         clo(i,j,1),cla(i,j,1)
!                          lon(i,j),lat(i,j) 
! clo(i,j,3),cla(i,j,3)                         clo(i,j,4),cla(i,j,4)

! sachant que le decallage d'indices entre lon_t et lon_f est le suivant
! lon_f(i,j+1)                       lon_f(i+1,j+1)
!                    lon_t(i,j) 
! lon_f(i,j)                         lon_f(i+1,j  )

! on a donc clo(i,j,1)=lon_f(i+1,j+1) ; clo(i,j,2)=lon_f(i  ,j+1) ; clo(i,j,3)=lon_f(i  ,j  ) ; clo(i,j,4)=lon_f(i+1,j  )

      grid_clo_ocean(i,j,1)=rad2deg*lon_f(i+1,j+1)  
      grid_clo_ocean(i,j,2)=rad2deg*lon_f(i  ,j+1)  
      grid_clo_ocean(i,j,3)=rad2deg*lon_f(i  ,j  )  
      grid_clo_ocean(i,j,4)=rad2deg*lon_f(i+1,j  )  

  enddo
  enddo

!    k=s_unit(7)
!    open(unit=k,file='oasis_input_files/grid_msk_ocean1')
!680 read(k,*,end=681)deci,decj,x1,x2 
!     do j0=0,1 ; do i0=0,1
!      i=int(deci)+i0 ; j=int(decj)+j0
!      if(i>=istr_oasis.and.i<=iend_oasis.and. &
!         j>=jstr_oasis.and.j<=jend_oasis)grid_msk_ocean(i,j)=0
!     enddo     ; enddo
!    goto 680
!681 close(k)

  !!!!!!!!!!!!!!!!! OASIS_WRITE_GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL oasis_start_grids_writing(flag)
  CALL oasis_write_grid(trim(namsrcgrd(cl_model_num))  ,iglb+2,jglb+2, grid_lon_ocean, grid_lat_ocean, il_part_id)
  CALL oasis_write_corner(trim(namsrcgrd(cl_model_num)),iglb+2,jglb+2, 4, grid_clo_ocean, grid_cla_ocean, il_part_id)
  CALL oasis_write_mask(trim(namsrcgrd(cl_model_num))  ,iglb+2,jglb+2, grid_msk_ocean(:,:), il_part_id)
  CALL oasis_terminate_grids_writing()

  deallocate(grid_lon_ocean)
  deallocate(grid_lat_ocean)
  deallocate(grid_clo_ocean)
  deallocate(grid_cla_ocean)


      call cpl_oasis_def_oasisplugs ! a besoin de grid_msk_ocean

  deallocate(grid_msk_ocean)

! Afin d'approfondir localement les profils verticaux de T etS transmis au
! modele destinataire (afin de diminuer les risques que ce dernier se retrouve
! sans information lorsque sa bathy est plus profonde que celle du modele
! emetteur) on prevoit d'ajouter un niveau k supplementaire determinE A partir
! du voisin "en mer" le plus profond dans le voisinage immediat. Les coordonnes i,j de ce
! voisin sont sauvegardees dans vert_extrap_i et vert_extrap_j. On note que
! la routine presente est appelee seulement A l'etat initiale et apres que la
! bathy du modele ait ete initialisee (en clair nous connaissons h_w)
  do j=jstr_oasis,jend_oasis 
  do i=istr_oasis,iend_oasis
! first guess:
     vert_extrap_i(i,j)=i ; vert_extrap_j(i,j)=j ; x1=h_w(i,j)
     do j1=-1,1
     do i1=-1,1
     if(i+i1>=1.and.i+i1<=imax.and.j+j1>=1.and.j+j1<=jmax) then !-bounds->
         if(h_w(i+i1,j+j1)>x1.and.mask_t(i+i1,j+j1,kmax)==1) then !ooo>
           vert_extrap_i(i,j)=i+i1 ; vert_extrap_j(i,j)=j+j1 ; x1=h_w(i+i1,j+j1)
         endif                                                    !ooo>
     endif                                                      !-bounds->
     enddo
     enddo
!    if(cl_model_num==1.and.mask_t(i,j,kmax)==1) &
!    write(10+par%rank,'(2(1x,i6,i6,1x,f10.3),a)') &
!    i,j,h_w(i,j)        &
!    ,vert_extrap_i(i,j) &
!    ,vert_extrap_j(i,j) &
!    ,x1,'      trouvemoi'
  enddo        
  enddo

#endif
  end subroutine cpl_oasis_grids_pmx

!...............................................................

  subroutine cpl_oasis_def_var_pmx
#ifdef key_oasis_sym_sym
  use module_parameter
  use module_principal ! only : k
  use module_parallele 
  use module_s
  implicit none

  allocate(group_send_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis,group_send_dim)) ; group_send_ocean=0.
  allocate(group_recv_ocean(istr_oasis:iend_oasis,jstr_oasis:jend_oasis,group_recv_dim)) ; group_recv_ocean=0.

  !!!!!!!!!!!!!!!!!! OASIS_DEF_VAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  var_nodims(1) = 2       ! Rank of the field array ; not used anymore in OASIS3-MCT
! var_nodims(2) = 1       ! Number of bundle fields
  var_actual_shape(1) = 1 ! Not used anymore in OASIS3-MCT
  var_type = OASIS_Real
  !
  ! Declaration of the coupling fields
  write(texte30,'(a,i0)')'R',cl_model_num ! noms dans namcouple, R1,R2 etc..., imperativement courts
  var_nodims(2)=group_recv_dim
  CALL oasis_def_var (var_oasis_id(1),trim(texte30), il_part_id, var_nodims, OASIS_In, var_actual_shape, var_type, ierror)
! CALL oasis_def_var (var_oasis_id(1),'GROUP_RECV_OCEAN1', il_part_id, var_nodims, OASIS_In, var_actual_shape, var_type, ierror)
! CALL oasis_def_var (var_oasis_id(1),'GROUP_RECV_'//trim(cl_model_name) &
!                    ,il_part_id,var_nodims,OASIS_In,var_actual_shape,var_type,ierror)

  write(texte30,'(a,i0)')'S',cl_model_num ! noms dans namcouple, S1,S2 etc..., imperativement courts
  var_nodims(2)=group_send_dim
  CALL oasis_def_var (var_oasis_id(2),trim(texte30),il_part_id,var_nodims,OASIS_Out,var_actual_shape,var_type,ierror)
! CALL oasis_def_var (var_oasis_id(2),'GROUP_SEND_'//trim(cl_model_name)   &
!                    ,il_part_id,var_nodims,OASIS_Out,var_actual_shape,var_type,ierror)

  !!!!!!!!!!!!!!!!!! OASIS_ENDDEF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL oasis_enddef (ierror)


  k=s_unit(7)
  write(texte4,'(i0)')par%rank
  open(unit=k,file=trim(tmpdirname)//trim(cl_model_name)//'_oasis_'//trim(texte4),position='append')
  write(k,'(a)')'........'
  write(k,*)'subroutine cpl_oasis_def_var_pmx'
  write(k,'(a,a)')'cl_model_name: ',trim(cl_model_name)
  write(k,'(a,a)')'oasis_dir:     ',trim(oasis_dir)
  write(k,*)'cl_model_num:',cl_model_num
  write(k,*)'par%rank',par%rank
  write(k,*)'group_recv_dim=',group_recv_dim
  write(k,*)'group_send_dim=',group_send_dim
  write(k,*)'il_part_id',il_part_id
  write(k,*)'var_oasis_id GROUP_RECV_OCEAN ',var_oasis_id(1)
  write(k,*)'var_oasis_id GROUP_SEND_OCEAN ',var_oasis_id(2)
  write(k,*)'ierror(oasis_enddef)=',ierror
  write(k,'(a)')'........'
  close(k)

#endif
  end subroutine cpl_oasis_def_var_pmx

!...............................................................

  subroutine cpl_oasis_get_pmx
#ifdef key_oasis_sym_sym
  use module_principal ! only : itime,elapsedtime_now,k
  use module_parallele 
  use module_s
  implicit none
  integer :: t_=2
! integer cpl_freqs_pmx(1)

    flag=0
    itime=nint(elapsedtime_now)
    if(mod(itime,namflddti(cl_model_num))/=0) goto 926 ! Si pas de reception OASIS aller direct A l'interpolation temporelle

    !!!!!!!!!!!!!!!!!!!!!!!! OASIS_GET !!!!!!!!!!!!!!!!!!!!!!
    CALL oasis_get(var_oasis_id(1),itime, group_recv_ocean, flag)

  if(flag/=0) then !>>> FLAG/=0 >>>

  oasis_get_count=oasis_get_count+1

! call oasis_get_freqs(var_oasis_id(1),OASIS_In,1,cpl_freqs_pmx,k0)
  k=s_unit(7)
  write(texte4,'(i0)')par%rank
  open(unit=k,file=trim(tmpdirname)//trim(cl_model_name)//'_oasis_'//trim(texte4),position='append')
  write(k,'(a)')'........'
  write(k,*)'subroutine cpl_oasis_get_pmx'
  write(k,'(a,a)')'cl_model_name: ',trim(cl_model_name)
  write(k,*)'par%rank         ',par%rank
  write(k,*)'var_oasis_id     ',var_oasis_id(1)
  write(k,*)'flag             ',flag
  write(k,*)'oasis_get_count  ',oasis_get_count
  write(k,*)'elapsedtime_now  ',elapsedtime_now
  write(k,*)'itime            ',itime
  write(k,*)'group_recv_ocean ',real(group_recv_ocean(imax/2,jmax/2,1:group_recv_dim))
! write(k,*)'modulo',mod(itime,namflddti(cl_model_num)),itime,namflddti(cl_model_num) 
  write(k,*)'lon,lat(imax/2,jmax/2)',real(lon_t(imax/2,jmax/2)*rad2deg) &
                                    ,real(lat_t(imax/2,jmax/2)*rad2deg)  
  close(k)

! Modeles enfants seulement
    if(cl_model_num/=1) then !-enfants->

!..........
! Bascule:
         sshobc_w(:,:,0)=   sshobc_w(:,:,2)
      velbarobc_u(:,:,0)=velbarobc_u(:,:,2)
      velbarobc_v(:,:,0)=velbarobc_v(:,:,2)
         temobc_t(:,:,:,0)= temobc_t(:,:,:,2)
         salobc_t(:,:,:,0)= salobc_t(:,:,:,2)
         velobc_u(:,:,:,0)= velobc_u(:,:,:,2)
         velobc_v(:,:,:,0)= velobc_v(:,:,:,2)
! Note 1: un temps envisagE, le basculement de l'indice temps plutot que le basculement du tableau,
! est repoussee A plus tard car il faut modifier les OBC oU la derivee temporelle de velobc a besoin
! d'etre adaptee.
! Note 2: la bascule ne concerne que le modele enfant car le modele 1 (parent) continue sa logique d'imbrication dans les fichiers MERCATOR-COPERNICUS
!..........

! Pour le modele 2 (modele enfant) qui est en retard d'une periode de couplage par rapport au
! modele 1 (modele parent) le temps correspondant au champ recu d'OASIS c'est le temps present + la periode de couplage
     oasis_time_next=elapsedtime_now+namflddti(cl_model_num)

      k1=   kmax_recv+1  ! bundle start temperature
      k2=2*(kmax_recv+1) ! bundle start salinite
      k3=3*(kmax_recv+1) ! bundle start courant OE
      k4=4*(kmax_recv+1) ! bundle start courant SN
      k5=5*(kmax_recv+1) ! bundle start courant SSH

! SSH:
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       sshobc_w(i,j,t_)=group_recv_ocean(i,j,k5+1)
      enddo                      ; enddo
! Interpoler verticalement champs 3D sur la grille verticale du modele recepteur:
        do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
         do k=1,kmax

          do k0=1,kmax_recv
           if(depth_t(i,j,k)>=group_recv_ocean(i,j,k0).and.depth_t(i,j,k)<=group_recv_ocean(i,j,k0+1)) then !pmx>
            rap=(depth_t(i,j,k)-group_recv_ocean(i,j,k0))/(group_recv_ocean(i,j,k0+1)-group_recv_ocean(i,j,k0))
            temobc_t(i,j,k,t_)=rap*group_recv_ocean(i,j,k0+k1+1)+(1.-rap)*group_recv_ocean(i,j,k0+k1)
            salobc_t(i,j,k,t_)=rap*group_recv_ocean(i,j,k0+k2+1)+(1.-rap)*group_recv_ocean(i,j,k0+k2)
              anyv3d(i,j,k,1)= rap*group_recv_ocean(i,j,k0+k3+1)+(1.-rap)*group_recv_ocean(i,j,k0+k3)
              anyv3d(i,j,k,2)= rap*group_recv_ocean(i,j,k0+k4+1)+(1.-rap)*group_recv_ocean(i,j,k0+k4)
           endif                                                                                            !pmx>
          enddo !k0

! cas particulier surface recepteur au dessus surface emetteur
          if(depth_t(i,j,k)  >=group_recv_ocean(i,j,kmax_recv+1)) then !aaa>
            temobc_t(i,j,k,t_)=group_recv_ocean(i,j,kmax_recv+k1+1)
            salobc_t(i,j,k,t_)=group_recv_ocean(i,j,kmax_recv+k2+1)
              anyv3d(i,j,k,1)= group_recv_ocean(i,j,kmax_recv+k3+1)
              anyv3d(i,j,k,2)= group_recv_ocean(i,j,kmax_recv+k4+1)
          endif                                                        !aaa>
! cas particulier fond recepteur au dessous fond emetteur
          if(depth_t(i,j,k)  <=group_recv_ocean(i,j,1)) then !bbb>
            temobc_t(i,j,k,t_)=group_recv_ocean(i,j,1+k1)
            salobc_t(i,j,k,t_)=group_recv_ocean(i,j,1+k2)
              anyv3d(i,j,k,1)= group_recv_ocean(i,j,1+k3)
              anyv3d(i,j,k,2)= group_recv_ocean(i,j,1+k4)
          endif                                              !bbb>

         enddo  !k

!--- test de verification ---
!         if(i==imax/2.and.j==jmax/2.and.mask_t(i,j,kmax)==1.and.itime==600) then
!          do k0=1,kmax_recv+1
!           write(10+par%rank,*)real(group_recv_ocean(i,j,k0)),real(group_recv_ocean(i,j,k0+k1)),real(group_recv_ocean(i,j,k0+k2))
!          enddo
!          do k=1,kmax
!           write(100+par%rank,*)real(depth_t(i,j,k)),temobc_t(i,j,k,t_),salobc_t(i,j,k,t_),real(tem_t(i,j,k,1)),real(sal_t(i,j,k,1))
!          enddo
!         endif
!--- test de verification ---

        enddo ; enddo   ! i,j 

!--- Continuite mpi ----
!mpi: echange "z0" sur sshobc_w(:,:,t_),temobc_t(:,:,:,t_) et salobc_t(:,:,:,t_)
      call obc_sshobc_mpi(t_)
      call obc_scal_temobc(t_)
      call obc_scal_salobc(t_)
      call obc_mpi_anyv3d(0,1,'za') ! zb sur anyv3d(:,:,:,1) sans option "canal" (soit 1er arg nul)
      call obc_mpi_anyv3d(0,2,'za') ! zb sur anyv3d(:,:,:,2) sans option "canal" (soit 1er arg nul)

!--- Passer des composantes OE,SN aux composantes "alongaxis" du courant Et interpolation sur points u et v
! velobc_u
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        velobc_u(i,j,k,t_)=0.25*( &
              (gridrotcos_t(i,j)+gridrotcos_t(i-1,j))*(anyv3d(i,j,k,1)+anyv3d(i-1,j,k,1)) &
             -(gridrotsin_t(i,j)+gridrotsin_t(i-1,j))*(anyv3d(i,j,k,2)+anyv3d(i-1,j,k,2)) &
                                )*mask_u(i,j,kmax) 
       enddo       ; enddo         ; enddo
! velobc_v
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        velobc_v(i,j,k,t_)=0.25*( &
              (gridrotsin_t(i,j)+gridrotsin_t(i,j-1))*(anyv3d(i,j,k,1)+anyv3d(i,j-1,k,1)) &
             +(gridrotcos_t(i,j)+gridrotcos_t(i,j-1))*(anyv3d(i,j,k,2)+anyv3d(i,j-1,k,2)) &
                                )*mask_v(i,j,kmax) 
       enddo       ; enddo         ; enddo

! Composantes moyennees sur la verticale
! velbarobc_u
      k=1
      do j=1,jmax ; do i=1,imax+1
       velbarobc_u(i,j,t_)=velobc_u(i,j,k,t_)*dz_u(i,j,k,1)
      enddo ; enddo
      do k=2,kmax ; do j=1,jmax ; do i=1,imax+1
       velbarobc_u(i,j,t_)= &
       velbarobc_u(i,j,t_)+velobc_u(i,j,k,t_)*dz_u(i,j,k,1)
      enddo ; enddo ; enddo
      do j=1,jmax ; do i=1,imax+1
       velbarobc_u(i,j,t_)= &
       velbarobc_u(i,j,t_)/hz_u(i,j,1)
      enddo ; enddo
! velbarobc_v
      k=1
      do j=1,jmax+1 ; do i=1,imax
       velbarobc_v(i,j,t_)=velobc_v(i,j,k,t_)*dz_v(i,j,k,1)
      enddo ; enddo
      do k=2,kmax ; do j=1,jmax+1 ; do i=1,imax
       velbarobc_v(i,j,t_)= &
       velbarobc_v(i,j,t_)+velobc_v(i,j,k,t_)*dz_v(i,j,k,1)
      enddo ; enddo ; enddo
      do j=1,jmax+1 ; do i=1,imax
       velbarobc_v(i,j,t_)= &
       velbarobc_v(i,j,t_)/hz_v(i,j,1)
      enddo ; enddo

      if(oasis_get_count==1) then !-cas-particulier-1->

! Cas particulier du premier GET (expliquE plus loin) oU, pour les besoins
! d'une initialisation intermediare, OBC(0)=OBC(2)=OBC(1) puisque, suite
! au passage par quickinitial, seul OBC(1) est coherent (les autres temps sont nuls)
         sshobc_w(:,:,0)=   sshobc_w(:,:,1)
      velbarobc_u(:,:,0)=velbarobc_u(:,:,1)
      velbarobc_v(:,:,0)=velbarobc_v(:,:,1)
         temobc_t(:,:,:,0)= temobc_t(:,:,:,1)
         salobc_t(:,:,:,0)= salobc_t(:,:,:,1)
         velobc_u(:,:,:,0)= velobc_u(:,:,:,1)
         velobc_v(:,:,:,0)= velobc_v(:,:,:,1)

         sshobc_w(:,:,2)=   sshobc_w(:,:,1)
      velbarobc_u(:,:,2)=velbarobc_u(:,:,1)
      velbarobc_v(:,:,2)=velbarobc_v(:,:,1)
         temobc_t(:,:,:,2)= temobc_t(:,:,:,1)
         salobc_t(:,:,:,2)= salobc_t(:,:,:,1)
         velobc_u(:,:,:,2)= velobc_u(:,:,:,1)
         velobc_v(:,:,:,2)= velobc_v(:,:,:,1)
      endif                       !-cas-particulier-1->

      if(oasis_get_count==2) then !-cas-particulier-2->
! Cas particulier du deuxieme GET (expliquE plus loin) oU, pour les besoins
! d'une initialisation intermediare, OBC(0)=OBC(1)=OBC(2) car seul OBC(2) est coherent
         sshobc_w(:,:,0)=   sshobc_w(:,:,2)
      velbarobc_u(:,:,0)=velbarobc_u(:,:,2)
      velbarobc_v(:,:,0)=velbarobc_v(:,:,2)
         temobc_t(:,:,:,0)= temobc_t(:,:,:,2)
         salobc_t(:,:,:,0)= salobc_t(:,:,:,2)
         velobc_u(:,:,:,0)= velobc_u(:,:,:,2)
         velobc_v(:,:,:,0)= velobc_v(:,:,:,2)
      endif                       !-cas-particulier-2->

    endif                    !-enfants->

    if(cl_model_num==1.and. &
       oasis_symsym_retrots==1) then !-parent->

! Modele parent seulement

     oasis_time_next=elapsedtime_now+namflddti(cl_model_num)

      k1=   kmax_recv+1  ! bundle start temperature
      k2=2*(kmax_recv+1) ! bundle start salinite
      k3=3*(kmax_recv+1) ! bundle start courant OE
      k4=4*(kmax_recv+1) ! bundle start courant SN
      k5=5*(kmax_recv+1) ! bundle start courant SSH

! Interpoler verticalement champs 3D sur la grille verticale du modele recepteur:
        do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
         do k=1,kmax

          do k0=1,kmax_recv
           if(depth_t(i,j,k)>=group_recv_ocean(i,j,k0).and.depth_t(i,j,k)<=group_recv_ocean(i,j,k0+1)) then !pmx>
            rap=(depth_t(i,j,k)-group_recv_ocean(i,j,k0))/(group_recv_ocean(i,j,k0+1)-group_recv_ocean(i,j,k0))
            tem_delta_t(i,j,k)=rap*group_recv_ocean(i,j,k0+k1+1)+(1.-rap)*group_recv_ocean(i,j,k0+k1)
            sal_delta_t(i,j,k)=rap*group_recv_ocean(i,j,k0+k2+1)+(1.-rap)*group_recv_ocean(i,j,k0+k2)
           endif                                                                                            !pmx>
          enddo !k0

! cas particulier surface recepteur au dessus surface emetteur
          if(depth_t(i,j,k)  >=group_recv_ocean(i,j,kmax_recv+1)) then !aaa>
            tem_delta_t(i,j,k)=group_recv_ocean(i,j,kmax_recv+k1+1)
            sal_delta_t(i,j,k)=group_recv_ocean(i,j,kmax_recv+k2+1)
          endif                                                        !aaa>
! cas particulier fond recepteur au dessous fond emetteur
          if(depth_t(i,j,k)  <=group_recv_ocean(i,j,1)) then !bbb>
            tem_delta_t(i,j,k)=group_recv_ocean(i,j,1+k1)
            sal_delta_t(i,j,k)=group_recv_ocean(i,j,1+k2)
          endif                                              !bbb>

         enddo  !k

        enddo ; enddo   ! i,j 

!--- Continuite mpi ----
!mpi: echange "z0" sur tem_delta_t et sal_delta_t
      call cpl_oasis_mpi_delta

    endif                    !-parent->

  endif            !>>> FLAG/=0 >>>

#ifdef bidon
! Ces lignes servent à verifier que les modeles avancent de maniere homogene en 
! verifiant l'ecriture de la date des 2 codes dans un meme fichier:
      if(par%rank==0) then !>>>>
       open(unit=3,file='mouchard.out',position='append')
        write(3,'(a,1x,i2,1x,a,1x,i4,1x,i2,a1,i2,a1,i2,a1)')  &
        trim(cl_model_name), & 
         day_now,trim(month(month_now)),year_now                 &     !20-11-10
        ,hour_now,'h'                                      &
        ,minute_now,'m',second_now,'s'
       close(3)
      endif                !>>>>
#endif

#ifdef bidon
! Cas test echangeant les latitudes/longitudes pour verifier les fonctions d'interpolation
      if(itime==600) then !OOO>
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       if(grid_msk_ocean(i,j)==0) then !ooo>
         if(cl_model_num==1) then !>>>>
          write(10+par%rank,*)real(lon_t(i,j)*rad2deg) &
                             ,real(lat_t(i,j)*rad2deg) &
                             ,real(group_recv_ocean(i,j,1:2))
         else                     !>>>>
          write(100+par%rank,*)real(lon_t(i,j)*rad2deg) &
                              ,real(lat_t(i,j)*rad2deg) &
                              ,real(group_recv_ocean(i,j,1:2))
         endif                    !>>>>
       endif                           !ooo>
      enddo                      ; enddo
      endif               !OOO>
#endif

  926 continue

!..............................
! Mise A jour tableaux OBC "now"
! Note: ne concerne que le modele enfant (modele 2) car le modele parent (modele 1)
! conserve son protocole d'imbrication dans les fichiers MERCATOR
! Modeles enfants seulement:
      if(cl_model_num/=1) then !m°v°m>

        timeweightobc(:)=1.-(oasis_time_next-elapsedtime_now)/namflddti(cl_model_num)

!.................................
! INITIALISATION DU MODELE ENFANT
! A quel moment initialiser le modele enfant avec les champs du modele parent?
! oasis_get_count compte les appels A la fonction oasis GET
! Au premier appel   (oasis_get_count=1) elapsedtime_now=0, la fonction GET recoit des champs nuls (inutilisables donc)
! Au deuxieme appel  (oasis_get_count=2) les tabeaux "obc" sont renseignes pour l'echeance "next" (indice temporel=2)
!                                        mais pas pour l'echeance "now" (indice temporel=0), donc initialisation impossible.
! Au troisieme appel (oasis_get_count=3) les tabeaux "obc" sont totalement renseignEs (now et next)
! On initialise donc quand oasis_get_count=3 et aussi quand flag/=0 qui isole l'iteration correspond pil poil A l'appel A 
! la fonction GET (car oasis_get_count reste =3 pendant plusieurs iterations successives, or c'est la premiere seulement qui 
! doit declencher l'initialisation. On initialise avec le temps now, soit le temps 0 en argument de call initial_with_obc

! Une initialisation temporaire est tout de meme appliquee pour eviter les instabilites aux frontieres. Ainsi, 
! quand oasis_get_count=2, afin d'eviter les incoherences decoulant de que OBC(0)=0 alors que OBC(2) est "normal"
! on force OBC(0)=OBC(2) et on initialise temporairement le modele en attendant l'initialisation definitive en oasis_get_count=3

      if(flag/=0) then !-flag-non-nul-> !la fonction GET vient de fonctionner
       if(oasis_get_count==2.or. & ! initialisation intermediaire
          oasis_get_count==3     & ! initialisation veritable
                               ) then      !-count-bon-pour-init->
         call initial_with_obc(0)
         call graph_out
       endif                               !-count-bon-pour-init->
      endif            !-flag-non-nul->

!.................................

      endif                    !m°v°m>

#endif
  end subroutine cpl_oasis_get_pmx

!...............................................................

  subroutine cpl_oasis_put_pmx
#ifdef key_oasis_sym_sym
  use module_principal ! only : itime,elapsedtime_now,k,i,j
  use module_parallele 
  use module_s
  implicit none

    itime=nint(elapsedtime_now)
    if(mod(itime+dti_fw_i4,namflddti(cl_model_num))/=0)return

! NOTE: pour le moment modeles parent et enfant echangent les memes champs mais
! je prevois par ce codage la possibilite de differentier les 2
    if(cl_model_num==1) then !m°v°m>

     k1=0

! Modele parent:
!... depth_t ....
! niveau additionnel "k=0" pour approfondir le profil en utilisant le point le
! plus profond du voisinage immediat dont les coordonnees sont i1=vert_extrap_i,j1=vert_extrap_j,k=1
     k1=k1+1
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      i1=vert_extrap_i(i,j) ;  j1=vert_extrap_j(i,j)
      group_send_ocean(i,j,k1)=min(depth_t(i1,j1,1),depth_t(i,j,1))
     enddo                      ; enddo
     do k=1,kmax
     k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=depth_t(i,j,k)
      enddo                      ; enddo
     enddo

!... temperature ....
! niveau additionnel "k=0" pour approfondir le profil
     k1=k1+1
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      i1=vert_extrap_i(i,j) ;  j1=vert_extrap_j(i,j)
      group_send_ocean(i,j,k1)=tem_t(i1,j1,1,1)
     enddo                      ; enddo
     do k=1,kmax
     k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=tem_t(i,j,k,1)
      enddo                      ; enddo
     enddo

!... salinite ....
! niveau additionnel "k=0" pour approfondir le profil
     k1=k1+1
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      i1=vert_extrap_i(i,j) ;  j1=vert_extrap_j(i,j)
      group_send_ocean(i,j,k1)=sal_t(i1,j1,1,1)
     enddo                      ; enddo
     do k=1,kmax 
     k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=sal_t(i,j,k,1)
      enddo                      ; enddo
     enddo

!... Courant Ouest-Est au point t ...
! niveau additionnel "k=0" pour approfondir le profil
     k1=k1+1
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      i1=vert_extrap_i(i,j) ;  j1=vert_extrap_j(i,j)
      group_send_ocean(i,j,k1)=0.
     enddo                      ; enddo
     do k=1,kmax 
     k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)= &
       0.5*(vel_u(i,j,k,1)+vel_u(i+1,j,k,1))*gridrotcos_t(i,j) &
      +0.5*(vel_v(i,j,k,1)+vel_v(i,j+1,k,1))*gridrotsin_t(i,j)
      enddo                      ; enddo
     enddo

!... Courant Sud-Nord au point t ...
! niveau additionnel "k=0" pour approfondir le profil
     k1=k1+1
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      i1=vert_extrap_i(i,j) ;  j1=vert_extrap_j(i,j)
      group_send_ocean(i,j,k1)=0.
     enddo                      ; enddo
     do k=1,kmax ! Courant Sud-Nord au point t
     k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)= &
      -0.5*(vel_u(i,j,k,1)+vel_u(i+1,j,k,1))*gridrotsin_t(i,j) &
      +0.5*(vel_v(i,j,k,1)+vel_v(i,j+1,k,1))*gridrotcos_t(i,j)
      enddo                      ; enddo
     enddo

!... SSH ...
     k1=k1+1 
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=ssh_int_w(i,j,1)
     enddo                      ; enddo

! Specifique (pour le moment) au modele parent (donc modele 1)
! Les points emetteurs non definis (car mask_t=0) mais necessaires au modele destinataire
! prennent la valeur du point en mer le plus proche. 
! Sont concernes: depth_t, T, S, SSH (mais pas u,v dont la valeur nulle quand mask=0 convient)
! et par consequent la boucle verticale va de k1=1 A k1=3*(kmax+1)+1
      do k10=1,oasismisspointmax
        i=oasisplugs(k10,1) ; j=oasisplugs(k10,2) ; i4=oasisplugs(k10,3) ; j4=oasisplugs(k10,4)  
! boucher depth_t
          do k1=1,kmax 
!           group_send_ocean(i,j,k1)=group_send_ocean(i4,j4,k1)
            group_send_ocean(i,j,k1+1)=depth_t(i4,j4,k1)
          enddo
          group_send_ocean(i,j,1)=group_send_ocean(i,j,2) ! niv supl sous le fond
! boucher tem_t
          do k1=1,kmax
!           group_send_ocean(i,j,k1)=group_send_ocean(i4,j4,k1)
            group_send_ocean(i,j,kmax+1+k1+1)=tem_t(i4,j4,k1,1)
          enddo
          group_send_ocean(i,j,kmax+1+1)=group_send_ocean(i,j,kmax+1+2) ! niv supl sous le fond
! boucher sal_t
          do k1=1,kmax
!           group_send_ocean(i,j,k1)=group_send_ocean(i4,j4,k1)
            group_send_ocean(i,j,2*(kmax+1)+k1+1)=sal_t(i4,j4,k1,1)
          enddo
          group_send_ocean(i,j,2*(kmax+1)+1)=group_send_ocean(i,j,2*(kmax+1)+2) ! niv supl sous le fond
! boucher ssh_w
          k1=5*(kmax+1)+1    ! SSH
!         group_send_ocean(i,j,k1)=group_send_ocean(i4,j4,k1)
          group_send_ocean(i,j,k1)=ssh_int_w(i4,j4,1)
      enddo

    else                     !m°v°m>

! Modele enfant: envoie les anomalies (champs - champs OBC)

! Note 1: 0 est une valeur de masque continental coherente pour une anomalie.
! l'extrapolation verticale (sous la bathy) applique aussi cette coherence
! ainsi que la moyenne sur chaque pavE oU chaque "zero" compte dans la moyenne 
! Note 2: a l'approche des frontieres ouvertes l'anomalie tend vers 0 avec la
! multiplication par 1-sponge_t

     k1=0

!... depth_t ....
! niveau additionnel "k=0" pour approfondir le profil en utilisant le point le
! plus profond du voisinage immediat dont les coordonnees sont i1=vert_extrap_i,j1=vert_extrap_j,k=1
     k1=k1+1
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      i1=vert_extrap_i(i,j) ;  j1=vert_extrap_j(i,j)
      group_send_ocean(i,j,k1)=min(depth_t(i1,j1,1),depth_t(i,j,1))
     enddo                      ; enddo
     do k=1,kmax
     k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=depth_t(i,j,k)
      enddo                      ; enddo
     enddo

!... delta temperature ....
! niveau additionnel "k=0" pour approfondir le profil
      k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=0. ! zero est une extrapolation coherente pour une anomalie
      enddo                      ; enddo
      x0_r4=1./9. ! moyenne sur pavE de 9 points
      do k=1,kmax
      k1=k1+1
!     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       do j2=2,jmax-1             ; do i2=2,imax-1
! Moyenne sur un pavE de 9 points coherente avec saut de resolution de 3
!       x1_r4=small1 ; x2_r4=0.
        x2_r4=0.
        do j=j2-1,j2+1 ; do i=i2-1,i2+1
!        x1_r4=x1_r4+mask_t(i,j,kmax)
         x2_r4=x2_r4+mask_t(i,j,kmax)* &
                                (tem_t(i,j,k,1) &
                             -temobc_t(i,j,k,1))
        enddo          ; enddo            ! j2,i2 loops
        group_send_ocean(i2,j2,k1)=x2_r4*x0_r4*mask_t(i2,j2,kmax) &
                                        *(1.-sponge_t(i2,j2,1))
!       group_send_ocean(i2,j2,k1)=x2_r4/x1_r4
!       group_send_ocean(i,j,k1)=tem_t(i,j,k,1) &
!                            -temobc_t(i,j,k,1)
       enddo                      ; enddo ! j,i loops
      enddo ! k loop

!... delta salinite ....
! niveau additionnel "k=0" pour approfondir le profil
      k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=0. ! zero est une extrapolation coherente pour une anomalie
      enddo                      ; enddo
      x0_r4=1./9. ! moyenne sur pavE de 9 points
      do k=1,kmax 
      k1=k1+1
!      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       do j2=2,jmax-1             ; do i2=2,imax-1
! Moyenne sur un pavE de 9 points coherente avec saut de resolution de 3
!       x1_r4=small1 ; x2_r4=0.
        x2_r4=0.
        do j=j2-1,j2+1 ; do i=i2-1,i2+1
!        x1_r4=x1_r4+mask_t(i,j,kmax)
         x2_r4=x2_r4+mask_t(i,j,kmax)* &
                                (sal_t(i,j,k,1) &
                             -salobc_t(i,j,k,1))
        enddo          ; enddo            ! j2,i2 loops
        group_send_ocean(i2,j2,k1)=x2_r4*x0_r4*mask_t(i2,j2,kmax) &
                                        *(1.-sponge_t(i2,j2,1))
!       group_send_ocean(i2,j2,k1)=x2_r4/x1_r4
!       group_send_ocean(i,j,k1)=sal_t(i,j,k,1) &
!                            -salobc_t(i,j,k,1)
       enddo                      ; enddo
      enddo

!... delta Courant Ouest-Est au point t ...
! niveau additionnel "k=0" pour approfondir le profil
      k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=0.
      enddo                      ; enddo
      x0_r4=1./9. ! moyenne sur pavE de 9 points
      do k=1,kmax 
       k1=k1+1
!      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       do j2=2,jmax-1             ; do i2=2,imax-1
!       x1_r4=small1 ; x2_r4=0.
        x2_r4=0.
        do j=j2-1,j2+1 ; do i=i2-1,i2+1
!        x1_r4=x1_r4+mask_t(i,j,kmax)
         x2_r4=x2_r4+mask_t(i,j,kmax)*( &
              0.5*(    vel_u(i  ,j,k,1) &
                   -velobc_u(i  ,j,k,1) &
                      +vel_u(i+1,j,k,1) &
                   -velobc_u(i+1,j,k,1) &
                                 )*gridrotcos_t(i,j) &
             +0.5*(    vel_v(i,j  ,k,1) &
                   -velobc_v(i,j  ,k,1) &
                      +vel_v(i,j+1,k,1) &
                   -velobc_v(i,j+1,k,1) &
                                 )*gridrotsin_t(i,j) &
                                      )
        enddo          ; enddo            ! j2,i2 loops
        group_send_ocean(i2,j2,k1)=x2_r4*x0_r4*mask_t(i2,j2,kmax) &
                                        *(1.-sponge_t(i2,j2,1))
!       group_send_ocean(i2,j2,k1)=x2_r4/x1_r4
!      group_send_ocean(i,j,k1)= &
!      0.5*(    vel_u(i  ,j,k,1) &
!           -velobc_u(i  ,j,k,1) &
!              +vel_u(i+1,j,k,1) &
!           -velobc_u(i+1,j,k,1) &
!                               )*gridrotcos_t(i,j) &
!     +0.5*(    vel_v(i,j  ,k,1) &
!           -velobc_v(i,j  ,k,1) &
!              +vel_v(i,j+1,k,1) &
!           -velobc_v(i,j+1,k,1) &
!                               )*gridrotsin_t(i,j)
       enddo                      ; enddo
      enddo

!... delta Courant Sud-Nord au point t ...
! niveau additionnel "k=0" pour approfondir le profil
      k1=k1+1
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=0.
      enddo                      ; enddo
      x0_r4=1./9. ! moyenne sur pavE de 9 points
      do k=1,kmax ! Courant Sud-Nord au point t
      k1=k1+1
!      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       do j2=2,jmax-1             ; do i2=2,imax-1
!       x1_r4=small1 ; x2_r4=0.
        x2_r4=0.
        do j=j2-1,j2+1 ; do i=i2-1,i2+1
!        x1_r4=x1_r4+mask_t(i,j,kmax)
         x2_r4=x2_r4+mask_t(i,j,kmax)*( &
             -0.5*(    vel_u(i  ,j,k,1) &
                   -velobc_u(i  ,j,k,1) &
                      +vel_u(i+1,j,k,1) &
                   -velobc_u(i+1,j,k,1) &
                                 )*gridrotsin_t(i,j) &
             +0.5*(    vel_v(i,j  ,k,1) &
                   -velobc_v(i,j  ,k,1) &
                      +vel_v(i,j+1,k,1) &
                   -velobc_v(i,j+1,k,1) &
                                 )*gridrotcos_t(i,j) &
                                      )
        enddo          ; enddo            ! j2,i2 loops
        group_send_ocean(i2,j2,k1)=x2_r4*x0_r4*mask_t(i2,j2,kmax) &
                                        *(1.-sponge_t(i2,j2,1))
!       group_send_ocean(i2,j2,k1)=x2_r4/x1_r4
!      group_send_ocean(i,j,k1)= &
!     -0.5*(    vel_u(i  ,j,k,1) &
!           -velobc_u(i  ,j,k,1) &
!              +vel_u(i+1,j,k,1) &
!           -velobc_u(i+1,j,k,1) &
!                               )*gridrotsin_t(i,j) &
!     +0.5*(    vel_v(i,j  ,k,1) &
!           -velobc_v(i,j  ,k,1) &
!              +vel_v(i,j+1,k,1) &
!           -velobc_v(i,j+1,k,1) &
!                               )*gridrotcos_t(i,j)
       enddo                      ; enddo
      enddo

!... delta SSH ...
      k1=k1+1 
!     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
      x0_r4=1./9. ! moyenne sur pavE de 9 points
      do j2=2,jmax-1             ; do i2=2,imax-1
!       x1_r4=small1 ; x2_r4=0.
        x2_r4=0.
        do j=j2-1,j2+1 ; do i=i2-1,i2+1
         x1_r4=x1_r4+mask_t(i,j,kmax)
         x2_r4=x2_r4+mask_t(i,j,kmax)*( &
                                 ssh_int_w(i,j,1) &
                                 -sshobc_w(i,j,1) &
                                      )
        enddo          ; enddo            ! j2,i2 loops
       group_send_ocean(i2,j2,k1)=x2_r4*x0_r4*mask_t(i2,j2,kmax) &
                                       *(1.-sponge_t(i2,j2,1))
!      group_send_ocean(i2,j2,k1)=x2_r4/x1_r4
!      group_send_ocean(i,j,k1)=ssh_int_w(i,j,1) &
!                               -sshobc_w(i,j,1)
      enddo                      ; enddo


    endif                    !m°v°m>

! temporaire (supprimer plus tard) ajout de 3 variables de verifications lon,lat,temps
!... longitude ...
     k1=k1+1 
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=lon_t(i,j)*rad2deg
     enddo                      ; enddo
!... latitude ...
     k1=k1+1 
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=lat_t(i,j)*rad2deg
     enddo                      ; enddo
!... temps ...
     k1=k1+1 
     do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       group_send_ocean(i,j,k1)=elapsedtime_now
     enddo                      ; enddo
! temporaire (supprimer plus tard) ajout de 3 variables de verifications lon,lat,temps

    !!!!!!!!!!!!!!!!!!!!!!!! OASIS_PUT !!!!!!!!!!!!!!!!!!!!!!
    CALL oasis_put(var_oasis_id(2),itime, group_send_ocean,flag)

  if(flag/=0) then !>>>
  k=s_unit(7)
  write(texte4,'(i0)')par%rank
  open(unit=k,file=trim(tmpdirname)//trim(cl_model_name)//'_oasis_'//trim(texte4),position='append')
  write(k,'(a)')'........'
  write(k,*)'subroutine cpl_oasis_put_pmx'
  write(k,'(a,a)')'cl_model_name: ',trim(cl_model_name)
  write(k,*)'par%rank',par%rank
  write(k,*)'var_oasis_id ',var_oasis_id(2)
  write(k,*)'flag ',flag
  write(k,*)'elapsedtime_now ',elapsedtime_now
  write(k,*)'itime           ',itime
  write(k,*)'modulo',mod(itime+dti_fw_i4,namflddti(cl_model_num)),itime+dti_fw_i4,namflddti(cl_model_num) 
  close(k)
  endif            !>>>

#endif
  end subroutine cpl_oasis_put_pmx

!------------------------------------------------------------------------------

      subroutine cpl_oasis_def_oasisplugs
#ifdef key_oasis_sym_sym
      use module_principal ; use module_parallele
      implicit none

! Pour les points en terre qui ont besoin d'avoir une valeur de T et S,
! chercher le point en mer le plus proche

      oasismisspointmax=0
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
! formulation du dilem: un point du modele emetteur necesaire au modele destinataire n'est
! pas defini car situE dans le masque continental (du modele emetteur):
       if(grid_msk_ocean(i,j)==0.and. & ! convention oasis: point necessaire au modele destinataire
          mask_t(i,j,kmax)==0)        & ! convention symphonie: point en terre
          oasismisspointmax=oasismisspointmax+1
      enddo                      ; enddo
      allocate(oasisplugs(oasismisspointmax,4))

      oasismisspointmax=0
      do j=jstr_oasis,jend_oasis ; do i=istr_oasis,iend_oasis
       if(grid_msk_ocean(i,j)==0.and.mask_t(i,j,kmax)==0) then !m°v°m>

          oasismisspointmax=oasismisspointmax+1

! chercher le point le plus proche en mer
! Notes: on cherche de 1 a imax et pas plus car les champs peuvent ne pas etre definis en 0 ou imax+1
          i4=i ; j4=j ; dist1=1e10 ! reset
          do j1=1,jmax   ; do i1=1,imax ! i1,j1 loop
             if(        mask_t(i1,j1,kmax)==1      &
                                     ) then !>>>

                 dist2=real(i-i1)**2+real(j-j1)**2
                 if(dist2<dist1) then         !---->
                   i4=i1 ; j4=j1 ; dist1=dist2
                 endif                        !---->

             endif                          !>>>
          enddo                       ; enddo                       ! i1,j1 loop

! Si aucun boucheur n'est trouvE le modele s'arrete:
         if(i==i4.and.j==j4) then !-debug->
          flag_stop=1
          write(10+par%rank,*)'No oasisplugs for point:',i+par%timax(1),j+par%tjmax(1)
         endif                    !-debug->

         oasisplugs(oasismisspointmax,1)=i 
         oasisplugs(oasismisspointmax,2)=j
         oasisplugs(oasismisspointmax,3)=i4
         oasisplugs(oasismisspointmax,4)=j4

!       if(cl_model_num==1) &
!       write(10+par%rank,'(4i6,a)')oasisplugs(oasismisspointmax,1:4),' trouvemoi'
          
       endif                                                   !m°v°m>
      enddo                      ; enddo

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ;
      if(k0/=0)stop 'Err cpl_oasis_def_oasisplugs see fort.xxx files'

#endif

      end subroutine cpl_oasis_def_oasisplugs
!..................................................................
      subroutine cpl_oasis_debug_pmx !17-08-20
#ifdef key_oasis_sym_sym
      use module_principal ; use module_parallele
      implicit none
! write(k,*)'Coupling periods (s) ',namflddti
! write(k,*)'Coupling lags    (s) ',namfldlag

! Verifier la conformite d'un certains nombres de parametres vis A vis
! de certaines contraintes liees A l'utilisation d'OASIS

! Verifier que le pas de temps correspond A un nombre entier de secondes
      if(dti_fw-dti_fw_i4/=0) then !>>>
       write(10+par%rank,*)'dti_fw doit etre un nombre entier de secondes'
       write(10+par%rank,*)'dti_fw   ',dti_fw
       write(10+par%rank,*)'dti_fw_i4',dti_fw_i4
       write(10+par%rank,*)'modele numero ',cl_model_num
       stop 'Err 1349 cpl_oasis_debug_pmx see fort.xxx files'
      endif                        !>>>

! Verifier que la periode de couplage correspond A un nombre entier d'iterations
      if(namflddti(cl_model_num)/dti_fw_i4/= &
         namflddti(cl_model_num)/dti_fw) then !>>>>
       write(10+par%rank,*)'La periode de couplage oasis doit' &
       ,' correspondre à un nombre entier d''iterations'
       write(10+par%rank,*)'dti_fw_i4        ',dti_fw_i4
       write(10+par%rank,*)'Periode couplage ',namflddti(cl_model_num)
       write(10+par%rank,*)'modele numero    ',cl_model_num
       stop 'Err 1359 cpl_oasis_debug_pmx see fort.xxx files'
      else                                    !>>>>
! Si c'est ok alors dtobc (requis par les OBC) du modele enfant est donnE par la 
! periode de couplage d'OASIS (ne concerne pas le modele parent forcE par MERCATOR)
       if(cl_model_num==2)dtobc(:)=namflddti(cl_model_num)
      endif                                   !>>>>

! verifer que les lags sont egaux aux pas de temps.
      if(namfldlag(cl_model_num)/=dti_fw_i4) then !>>>
       write(10+par%rank,*)'Le "LAG" est different du pas de temps'
       write(10+par%rank,*)'dti_fw_i4 ',dti_fw_i4
       write(10+par%rank,*)'LAG       ',namfldlag(cl_model_num) 
       stop 'Err 1367 cpl_oasis_debug_pmx see fort.xxx files'
      endif                                       !>>>

! Forcages aux OBC:
! Maree:
      if(cl_model_num/=1.and.kmaxtide>0.and.tideforces<4) then !ooo>
       write(10+par%rank,*)'Pas de maree aux OBC pour modeles enfants'
       stop 'Err 1550 cpl_oasis_debug_pmx see fort.xxx files'
      endif                                                    !ooo>
! Copernicus:
      if(cl_model_num/=1.and.iobc_ogcm==1) then !ooo>
       write(10+par%rank,*)'Pas de fichiers OGCM aux OBC des modeles enfants'
       stop 'Err 1557 cpl_oasis_debug_pmx see fort.xxx files'
      endif                                                    !ooo>

#endif
      end subroutine cpl_oasis_debug_pmx
!..................................................................
      subroutine cpl_oasis_reset_pmx !17-08-20
#ifdef key_oasis_sym_sym
      use module_principal ; use module_parallele
      implicit none

! Verifier la conformite de certains parametres vis A vis des contraintes
! d'OASIS (si utilisE)
      call cpl_oasis_debug_pmx !17-08-20

      if(cl_model_num/=1)bi_onoff=0 ! modele enfant recoit du modele parent une SSH "brute" contenant l'effet atmospherique (donc on ne l'ajoute pas et bi_onoff=0)

      if(cl_model_num==1.and. &
         oasis_symsym_retrots==1) then !-allocate-delta->
       allocate(tem_delta_t(0:imax+1,0:jmax+1,kmax)) ; tem_delta_t=0.
       allocate(sal_delta_t(0:imax+1,0:jmax+1,kmax)) ; sal_delta_t=0.
      endif                            !-allocate-delta->

#endif
      end subroutine cpl_oasis_reset_pmx !17-08-20

!..................................................................

      subroutine cpl_oasis_mpi_delta
#ifdef key_oasis_sym_sym
      use module_principal ; use module_parallele
      implicit none
      integer id_tem_delta_z0_,loop_,id_sal_delta_z0_

#ifdef parallele
       write(texte30,'(a)')'tem_delta_z0_'
       call get_type_echange('z0',trim(texte30)   & !15-07-14
                                  ,tem_delta_t       &
                           ,lbound(tem_delta_t)      &
                           ,ubound(tem_delta_t)      &
                           ,id_tem_delta_z0_)

       write(texte30,'(a)')'sal_delta_z0_'
       call get_type_echange('z0',trim(texte30)   & !15-07-14
                                  ,sal_delta_t       &
                           ,lbound(sal_delta_t)      &
                           ,ubound(sal_delta_t)      &
                           ,id_sal_delta_z0_)

      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(tem_delta_t,id_tem_delta_z0_,mpi_neighbor_list(loop_))
        call echange_voisin(sal_delta_t,id_sal_delta_z0_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif
#endif
      end subroutine cpl_oasis_mpi_delta

!..................................................................

      subroutine cpl_oasis_tem_delta
#ifdef key_oasis_sym_sym
      use module_principal ; use module_parallele
      implicit none

! On ressort aussitot si 2 conditions:
       if(    oasis_symsym_retrots==0    & ! condition 1: on n'a pas demandE la retroaction sur T et S
          .or.cl_model_num/=1            & ! condition 2: le present modele n'est pas le modele parent
         )return

       x1_r4=0.5*dti_fw/namflddti(cl_model_num)
! Note: la multiplication par 0.5 pour ne corriger que la moitiE de l'erreur sur
! le lapse de temps correspondant A la periode de couplage

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        tridia_in(i,j,k,4)=tridia_in(i,j,k,4)+dz_t(i,j,k,1)*tem_delta_t(i,j,k)*x1_r4
       enddo       ; enddo       ; enddo

#endif
      end subroutine cpl_oasis_tem_delta

!..................................................................

      subroutine cpl_oasis_sal_delta
#ifdef key_oasis_sym_sym
      use module_principal ; use module_parallele
      implicit none

! On ressort aussitot sous certaines conditions:
       if(    oasis_symsym_retrots==0    & ! condition 1: on n'a pas demandE la retroaction sur T et S
          .or.cl_model_num/=1            & ! condition 2: le present modele n'est pas le modele parent
          .or.oasis_get_count<3          & ! condition 3: le spin-up des "get/put" d'oasis n'est pas terminE
         )return

       x1_r4=0.5*dti_fw/namflddti(cl_model_num)
! Note: la multiplication par 0.5 pour ne corriger que la moitiE de l'erreur sur
! le lapse de temps correspondant A la periode de couplage

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        tridia_in(i,j,k,4)=tridia_in(i,j,k,4)+dz_t(i,j,k,1)*sal_delta_t(i,j,k)*x1_r4
       enddo       ; enddo       ; enddo

#endif
      end subroutine cpl_oasis_sal_delta

!..................................................................

  subroutine cpl_oasis_init_pmx(id_lcomm)
  use module_principal ! , only : tmpdirname,texte4,k ! ,istr_oasis,iend_oasis,jstr_oasis,jend_oasis
  implicit none
! integer, INTENT(OUT) :: id_lcomm                   ! Model local communicator
  integer              :: id_lcomm                   ! Model local communicator
#ifdef key_oasis_sym_sym

! notebook_oasis
  namelist/notebook_oasis_sym_sym/cl_model_name          & 
                                 ,oasis_dir              &
                                 ,cl_model_num           &
                                 ,kmax_recv              &
                                 ,oasis_symsym_onoff     & 
                                 ,oasis_symsym_retrots

  !==
  ! Initialisation de OASIS
  !==

! notebook_oasis
  open(100,file=nomfichier(32)) ! 'notebook_oasis
  read(100,nml=notebook_oasis_sym_sym)
  close(100)

!   PRINT_DBG*, 'enter CPL_OASIS_INIT routine'
   
   !! Initialise the coupling
   CALL oasis_init_comp (il_compid, cl_model_name, il_err)

   IF (il_err /= 0) THEN
      CALL oasis_abort(il_compid, 'cpl_oasis_init_pmx', 'Problem during oasis_init_comp')
   ENDIF
   
   !! Get the value of a local MPI communicator to be used by SYMPHONIE for its internal parallelisation
   CALL oasis_get_localcomm (id_lcomm, il_err)
   IF (il_err /= 0) THEN
      CALL oasis_abort(il_compid, 'cpl_oasis_init_pmx', 'Problem during oasis_get_localcomm')
   ENDIF

!   PRINT_DBG*, 'exit CPL_OASIS_INIT routine'
   if(cl_model_num==1) then !m°v°m>
! Modele No1 = modele parent
    group_send_dim=1+5*(kmax+1)      ! kmax local 
    group_recv_dim=1+5*(kmax_recv+1) ! kmax de l'autre modele
   else                     !m°v°m>
! Modele No2 = modele enfant
    group_send_dim=1+5*(kmax+1)      ! kmax local 
    group_recv_dim=1+5*(kmax_recv+1) ! kmax de l'autre modele
   endif                    !m°v°m>

! temporairement (a supprimer plus tard) on ajoute 3 variables de verification: lon, lat, temps
    group_send_dim=group_send_dim+3
    group_recv_dim=group_recv_dim+3

#endif
  end subroutine cpl_oasis_init_pmx

!..................................................................

  subroutine cpl_oasis_finalize
  implicit none

#ifdef key_oasis
   call oasis_terminate (il_err)
   if (il_err /= 0) then
      CALL oasis_abort(il_compid, 'cpl_prism_finalize', 'Problem during oasis_terminate')
   endif  
#endif
#ifdef key_oasis_sym_sym
   call oasis_terminate (il_err)
   if (il_err /= 0) then
      CALL oasis_abort(il_compid, 'cpl_prism_finalize', 'Problem during oasis_terminate')
   endif  
#endif

  end subroutine cpl_oasis_finalize

!..................................................................

end module module_cpl_oasis
