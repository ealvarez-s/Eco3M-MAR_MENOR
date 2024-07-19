      subroutine reset
!______________________________________________________________________
! SYMPHONIE ocean model
! release 361 - last update: 01-01-23
!______________________________________________________________________

!...............................................................................
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................



      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='reset'
       subroutinedescription= &
       'Initial value of some constants and arrays '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!......................................................................
! Version date      Description des modifications
!         27/06/01: mise à zero des conditions aux limites:
!         17/12/01: bienvenue à la base orthonormee NOVECTOR(4,4)
!         28/02/02: reset de tridia_in à 1 pour eviter les divisions
!                   par zero sur les bords de domaine
!         01/03/02: tridia_in(2) seulement doit etre = 1
!         26/08/02: VMEAOBC remplacé par VMEAOBC
!                   & VHZOBC remplacé par VELOBC
!         23/09/02: debug 3eme argument flux atmos
!         26/12/02: prise en compte des nouvelles dimensions de OMEGA_Z
!         28/01/03: bienvenue à KLOW_X, KLOW_Y, KLOW_Z
!         19/03/03: bienvenue à PRECIPI_Z
!         02/11/03: plus de memoire plus NRJ3D
!         06/05/04: bienvenue NCMIN_AIRSEA NCMIN_RIVER NCMIN_AF
!         27/07/04: reset special de RIVER_T et RIVER_S pour debugeur
!                   dans set_rivers(3)
!         22/05/06: La dimension des tableaux OBC dépend de NOBCMIN & NOBCMAX
!         26/03/07: reset de tableaux meteo
!         09/05/07: reset des tableaux ADVMIX_X & ADVMIX_Y
!         19/01/08: reset de KSL à NR-1. Inportant au cas où l'appel à
!                   la routine KSL est deactive dans scalars.F
! 2009.3  02-10-09  Reset albebo
!         08-10-09  albedo devient albedo_z en 2d
!         09-10-09  reset de lagrange_zta
! 2010.7  22-02-10  reset timefilter_u et timefilter_v
!         26-02-10  reset vitesses de stokes
!         03-03-10  reset stokes vortex forces
!         04-03-10  reset sshstokes
! 2010.8  11-03-10  reset variables d'etat
!         18-03-10  reset wetmask=1.
!         19-03-10  reset heatrelax_w
!         03-05-10  nouvelles dim pour omega, veldxdz_v veldydz_u
! 2010.10 13-06-10  suppression ssh_ext_w
! 2010.11 16-07-10  temobc & salobc renommés temobc & salobc
! 2010.12 03-10-10  suppression obcmin et obcmax
! 2010.15 28-12-10  reset de upwindriver_t à 1 (sinon reste à zéro si
!                   aucun fleuve)
! 2010.18 23-02-11  Reset kmin en array syntax
! 2010.20 14-04-11  Reset elapsedtime
!         16-04-11  Reset elapsedtime_aft
! 2010.23 21-05-11  Reset tab1_code tab2_code tab_code_glb
!         02-06-11  la variable model_name porte le nom du modele
! 2010.24 03-12-11  Reset model name
!         15-12-11  reset array syntax pour pss_w evite erreur de mpi
!                   quand ifb=0
! 2010.25 04-01-12  reset model_name
!         03-04-12  Ajout du nom de la config dans les fichiers netcdf
!         31-05-12  model_name='S25.3'
!         31-05-12  model_name='S25.4'
! S26     05-11-12  ajout vert_axis_conv_direc=0 hori_axis_conv_start=0
!         09-03-13  ajout subcycle_echange et synchro
!         04-05-13  ajout type_structured type_unstructured
!         03-06-13  reset angle0
!         05-06-13  flag3d=1
!         30-06-13  subcycle_exchange
!         21-09-13  reset obc2dtype
!         17-11-13  obc2dtype devient real
!         16-01-14  fbtfiltercoef=0.01
!         14-02-14  reset de mpi_neighbor_list
!         16-07-14  offline_init_status=0
!         23-01-15  obcstatus
!         11-03-15  si wgrid alors allouer sigma_w
!         06-09-15  obcstatus depend de la periodicite en i ou j
!         16-02-16  norme gfortran
! v285    04-06-20  sigma_w alloue par module_principal
! v361    01-01-23  reset looplimit_hor=100 
!...............................................................................

      looplimit_hor=100 ! reset !01-01-23

      call reset_obcstatus !23-01-15

!     model_name='s26.bobshelf.20141113'
      open(unit=3,file='../../SOURCES/model_name',status='old')
       read(3,'(a)')model_name
      close(3)

      wgrid_case=0 !18-03-14
      fgrid_case=1
      fgrid_or_wgrid=wgrid_case
!     fgrid_or_wgrid=fgrid_case


!     if(fgrid_or_wgrid==wgrid_case)allocate(sigma_w(0:imax+1,0:jmax+1,0:kmax+1)) ; sigma_w=0 !11-03-15

      offline_init_status=0 !16-07-14

      mpi_neighbor_list=(/ ouest, est, nord, nordouest, nordest , sud, sudouest, sudest, ouest2, est2/) !14-02-14
      subcycle_exchange=8 !30-06-13
      subcycle_synchro=1
      fbtfiltercoef=0.01 !16-01-14

      obc2dtype=1.        !17-11-13
      type_structured=0   !04-05-13
      type_unstructured=1
      ww3_type_grid=0
      angle0=0.  !03-06-13
      flag3d=1      !05-06-13

      title_for_netcdf_files_txt='no title'
      open(unit=3,file='title_for_netcdf_files')            !03-04-12
       read(3,'(a)',end=68)title_for_netcdf_files_txt
   68 close(3)

!     vert_axis_conv_direc='up'                                   !05-11-12
!     vert_axis_conv_start='w'  ! fichier demarre par point w       !05-11-12
!     vert_axis_conv_end='w'    ! fichier finit   par point w       !05-11-12
!     hori_axis_conv_start='t'  ! fichier demarre par point t

      elapsedtime_now=0.                                              !14-04-11
      elapsedtime_bef=0.                                              !14-04-11
      elapsedtime_aft=0.                                              !14-04-11

      upwindriver_t(:,:)=1.                                             !28-12-10

      lagrange_ssh=0.                                                   !09-10-09

! Par defaut KSL_Z(I,J)=NR-1 sinon catastrophe si l'appel à la routine
! KSL est deactivé dans scalars.F
      do j=1,jmax                                                      !19/01/08
      do i=1,imax                                                      !19/01/08
       ksl_t(i,j)=kmax
      enddo
      enddo

! Seuil mini pour numero de record des fichiers de forcage:
      ncmin_airsea=-999
      ncmin_river=-999

!......................................................................!27/07/04
! reset special de RIVER_T et RIVER_S pour debugeur dans set_rivers(3)
      do k1=0,2
      do k2=1,dim_river
       river_s(k2)=-1.e10
       river_t(k2,k1)=-1.e10
      enddo
      enddo
!......................................................................

!***********************************************************************
! MISE A ZERO DE CERTAINS TABLEAUX:
!***********************************************************************
! Remarque: la mise à zero des tableaux peut paraitre superflu.
! Ce n'est pas exact dans le cas d'un ensemble de simulations genere
! au cours d'un seul run.
! reset des bornes inferieures des boucles verticales

      kmin_u(:,:)=1                                                    !23-02-11
      kmin_v(:,:)=1                                                    !23-02-11
      kmin_w(:,:)=1                                                    !23-02-11


! modif 27/06/01:
! TABLEAUX OBC:
      sumarchive=0.

      wetmask_u(:,:)=1.                                                !18-03-10
      wetmask_v(:,:)=1.                                                !18-03-10
      wetmask_t(:,:)=1.                                                !18-03-10


      pss_mean(:)=101300.       !15-12-11
      pss_w(:,:,:)=pss_mean(1)

      do k1=1,4
      do  k=1,kmax+1
      do  j=0,jmax+1
      do  i=0,imax+1
      tridia_in(i,j,k,k1)=0.                                           !28/02/02
      if(k1.eq.2)tridia_in(i,j,k,k1)=1.                                !01/03/02
      enddo
      enddo
      enddo
      enddo

!________________________________________________________________
! INITIALISATION DES MASQUES DES NIVEAUX VERTICAUX:
!     passehaut(1)=0.
!     passebas(1)=1.
!     do 40 k=2,kmax-1
!     passehaut(k)=1.
!     passebas(k)=1.
!  40 continue
!     passehaut(kmax)=1.
!     passehaut(kmax+1)=1.
!     passebas(kmax)=0.
!     passebas(kmax+1)=0.

!_______________________________________________________________       !17/12/01
! Construction d'une base orthonormée:
! Début:
!_______________________________________________________________       !17/12/01
      novector(1,1)=0.5
      novector(2,1)=0.5
      novector(3,1)=0.5
      novector(4,1)=0.5

      novector(1,2)=3.
      novector(2,2)=1.
      novector(3,2)=-1.
      novector(4,2)=-3.

      novector(1,3)=0.5
      novector(2,3)=-0.5
      novector(3,3)=-0.5
      novector(4,3)=0.5

      x1=1.
      novector(1,4)= 1./3.
      novector(2,4)=-1.
      novector(3,4)= 1.
      novector(4,4)=-1./3.


! normalisation des vecteurs:
      do j=1,4
       sum1=0.
       do i=1,4
       sum1=sum1+novector(i,j)**2
       enddo
       do i=1,4
       novector(i,j)=novector(i,j)/sqrt(sum1)
       enddo
      enddo

!c     DO J1=1,4
!c     WRITE(6,*)'___________________'
!c     DO J2=1,4
!c     SUM1=0.
!c     DO I=1,4
!c     SUM1=SUM1+NOVECTOR(I,J1)*NOVECTOR(I,J2)
!c     ENDDO
!c     WRITE(6,*)'verification',SUM1
!c     ENDDO
!c     ENDDO
!c     stop 'dans reset'
!_______________________________________________________________       !17/12/01
! Construction d'une base orthonormée:
! Fin.
!_______________________________________________________________       !17/12/01

!     CALL GRAPH_OUT
!     WRITE(6,*)'DANS SBR RESET'
!     WRITE(6,*)'-----------------------'
!     stop 'reset'

      end subroutine reset

!...................................................................

      subroutine reset_obcstatus
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='reset_obcstatus'
       subroutinedescription= &
       'lateral boundary status (obcstatus array)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! If obcstatus=1 apply radiative or sponge boundary condition !09-12-14
!    obcstatus=0 otherwise
      ieq1=ouest  !16-12-14
      ieqimax=est
      jeq1=sud
      jeqjmax=nord
      ieq1_jeq1=sudouest
      ieq1_jeqjmax=nordouest
      ieqimax_jeq1=sudest
      ieqimax_jeqjmax=nordest
      if(nordest  <0.or.nordest  >10)stop 'Error 59 main.F90'
      if(sudest   <0.or.sudest   >10)stop 'Error 59 main.F90'
      if(nordouest<0.or.nordouest>10)stop 'Error 59 main.F90'
      if(sudouest <0.or.sudouest >10)stop 'Error 59 main.F90'
      if(nord     <0.or.nord     >10)stop 'Error 59 main.F90'
      if(sud      <0.or.sud      >10)stop 'Error 59 main.F90'
      if(est      <0.or.est      >10)stop 'Error 59 main.F90'
      if(ouest    <0.or.ouest    >10)stop 'Error 59 main.F90'
      obcstatus=0

      if(.not.(iperiodicboundary)) then !>>>>>           !16-02-16
       if(     par%timax(1)==0   )obcstatus(ouest)=1     !16-12-14
       if(imax+par%timax(1)==iglb)obcstatus(est  )=1     !16-12-14
      endif                           !>>>>>

      if(.not.(jperiodicboundary)) then !>>>>>           !16-02-16
       if(     par%tjmax(1)==0   )obcstatus(sud  )=1
       if(jmax+par%tjmax(1)==jglb)obcstatus(nord )=1
      endif                           !>>>>>

      obcstatus(nordest  )=obcstatus(nord)*obcstatus(est)
      obcstatus(sudest   ) =obcstatus(sud)*obcstatus(est)
      obcstatus(nordouest)=obcstatus(nord)*obcstatus(ouest)
      obcstatus(sudouest) =obcstatus(sud) *obcstatus(ouest)

      end subroutine reset_obcstatus
!...................................................................
