      module module_my_outputs
!...............................................................................
! SYMPHONIE ocean model
! release 297 - last update: 06-03-21
!...............................................................................
! Version date      Description des modifications
!         30-07-15: ecrire la date sous forme calendaire
!         31-07-15: my_outputs_read_point_vs_time permet de relire les fichiers
!                   les stations julio et cie... 
!         19-10-15  debug ecriture netcdf
!         27-11-15  ecriture pour cas 1DV
!         30-11-15  ajout variables de la couche limite atmospherique de S26
!         14-02-16  Ajout des flux dans les sorties du modele 1DV
!         15-02-16  Ajout densite potentielle
!         30-03-16  - Debug signe N2
!                   - Ajout temlwf sallwf
!         02-04-16  attribut concernant temlwg et sallwf
!         05-04-16  ajout teta0, teta2 en Celsius
!         19-04-16  modif C.L. visu N2
!         22-04-16  Ajout d'un lien google
!         26-04-16  iturbulence peut etre egal A 2
!         23-08-16  debug boucle cas iturbulence=2
!         09-09-16  sorties "maregraphe"
!         10-10-16  info i,j dans attribut
!         28-10-16  Ajout de bio_t(1) dans sortie 1DV
!         25-01-17  ajout subroutine my_outputs_bio_sum
!         04-02-17  ajout subroutine my_outputs_tem_sum
!                   ajout subroutine my_outputs_sal_sum
!         16-03-17  ajout des "one cell river boxes" dans le bilan de sel de
!                   la subroutine my_outputs_sal_sum
!         19-03-17  modif nom du fichier salmean
!         02-07-17  ajout my_outputs_straitflux_i
!         20-05-18  ajout omega_w dans fichier 1DV
!         13-12-18  ajout along axis current dans fichier 1DV
!         14-01-19  ajout des flux barotropes dans les fichers 1DV
! v245    29-01-19  ajout subroutine my_outputs_obcsaltflux
! v249    10-03-19  ajout ssh_w_int_w(i,j,1) dans fichier 1DV
!         11-03-19  - call datetokount(year_,1,1,0,0,0) 
!                   - debug division par zero
! v250    26-03-19  arreter la simulation si le point 1DV est dans le land mask
!         27-03-19  Les stations PEACETIME pretes A l'emploi et des seuils min
!                   ajoutes A km_w et kh_w
!         31-03-19  ajout cas if(trim(locationconv_)=='latlon')
! v256    05-06-19  ajout subroutine my_outputs_glider
!         11-06-19  ajout my_outputs_zonesalttempflux
! v257    05-07-19  suite point precedent
! v258    17-07-19  zone_mask doit etre declare pour module_offline
! v259    05-09-19  ajout variable bio_t dans fichier 1DV
! v261    23-10-19  my_outputs_ssh_int_sum 
! v277    13-04-20  version 3D de my_outputs_zonesalttempflux
! v289    24-09-20  mise A jours my_outputs_zone1salttempflux
! v293    15-12-20  changement dans l'ecriture de la date dans "x1" pour
! v294    18-12-20  une nouvelle version du calcul des flux adaptee de la 
!                   version basique par Marine Herrmann 18-12-20     
! v295    29-01-21  depth_t remplace depth_w
! v296    02-03-21  subroutine my_outputs_timeabovewater !02-03-21
! v297    06-03-21  mise a jour du bilan de T et S suite aux modifs de l'advecion
!...............................................................................
!    _________                    .__                  .__             !(°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
      use module_principal
      use module_parallele
      use module_global
     use ModuleDeclaration
      implicit none
!     double precision                       &
!      cumulveltime
!     real*4,dimension(:,:,:),allocatable :: &
!      cumulvel_t
!     real*4,dimension(:,:,:,:),allocatable :: &
!      cumultime_t


#ifdef bidon
! declaration pour subroutine my_outputs_timeabovewater
!     integer , dimension(:,:),allocatable :: timeabovewater_w
      real*4  , dimension(:,:),allocatable :: timeabovewater_w
      integer :: timeabovewater_count=0
#endif

      integer :: ncid_,time_counter_=0 &
                ,i_,j_,flag_header_=0  &
                ,year_
      real*4,dimension(:,:),allocatable :: &
       angle_wave_beam_w
      double precision,dimension(:),allocatable :: &
       seconds_since_1jan

      double precision, dimension(:),allocatable :: zprofile_depth !,zprofile_tem,zprofile_sal

      character*14 :: date_
#ifdef bidon
      double precision :: flux_lowfreq=0.
#endif
#ifdef bidon
! Declaration pour la subroutine my_outputs_obcsaltflux
! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
      integer :: ifsbeg=-999,ifsend=-999,jfsbeg=-999,jfsend=-999,sponge_l_int=0 ! "free sponge (fs)" indice min max values
      integer :: nlayer_=2
      double precision , dimension(:) , allocatable ::      &
                          obcsaltfluxi1_,obcsaltfluxi1_glb_ &
                         ,obcsaltfluxi2_,obcsaltfluxi2_glb_ &
                         ,obcsaltfluxj1_,obcsaltfluxj1_glb_ &
                         ,obcsaltfluxj2_,obcsaltfluxj2_glb_ 
      real*4 , dimension(:) , allocatable ::  layerdepth_
#endif
#ifdef bidon
! Declaration pour la subroutine my_outputs_zonebioflux
      double precision , dimension(:,:)   , allocatable :: zone1biomasst0
      double precision , dimension(:,:,:) , allocatable ::      & !nb de limites*nb de couches,vb
                          zone1bioflux_glb,zone1bioflux_u,zone1bioflux_v
      double precision , dimension(:,:) , allocatable ::  & !vb
                          zone1biosurfflux_w,zone1biobotflux_w,zone1biorivflux_w &
                         ,zone1tendancebio_glb     &
                         ,zone1botsurfbio_glb
      double precision , dimension(:) , allocatable ::  & !vb
                          zone1biocumul_glb             &
                         ,zone1biocumul_loc              
#endif

contains

      subroutine my_outputs_driver
      implicit none
#ifdef synopsis
       subroutinetitle='my_outputs_driver'
       subroutinedescription= &
       'Driving routine of the "home made" users ' &
       //' output/diag/post-treatment etc... subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      flag_stop=0 !26-03-19

      if(initial_main_status==1)then !sssssss>


!        call my_outputs_timeabovewater
      
      
! Le tableau seconds_since_1jan va servir a construire une date en annee decimale
         if(.not.allocated(seconds_since_1jan)) then !>>>>>
           allocate(seconds_since_1jan(datesim(1,1):datesim(1,2)+1)) ! 1ere annee -> derniere annee+1
           seconds_since_1jan=0. !11-03-19
           do year_=datesim(1,1),datesim(1,2)+1
              call datetokount(year_,1,1,0,0,0) !11-03-19
            seconds_since_1jan(year_)=elapsedtime_out
           enddo
         endif                                       !>>>>>
         if(flag_1dv==1)call my_outputs_point_vs_time('1dv'  ,'ijglob',2.d0,2.d0)

!        texte250='julio_ref.nc'
!        texte250='julio_optics.nc'
!        texte250='julio_prodi_part2.nc'
!        texte250='julio_prodi_part3.nc'
!        texte250='julio_molecular.nc'
!        texte250='julio_clas.nc'
!        call my_outputs_read_point_vs_time(trim(texte250),5.255d0,43.135d0,2013,10,01,09,0,0)
!        call my_outputs_read_point_vs_time(trim(texte250),5.255d0,43.135d0,2013,12,10,11,0,0)
!        call my_outputs_read_point_vs_time(trim(texte250),5.255d0,43.135d0,2013,05,11,8,0,0)
!        call my_outputs_read_point_vs_time(trim(texte250),5.255d0,43.135d0,2014,03,13,11,0,0)
!        call my_outputs_read_point_vs_time(trim(texte250),5.255d0,43.135d0,2013,05,15,12,0,0)
!        call my_outputs_read_point_vs_time(trim(texte250),5.255d0,43.135d0,2011,10,18,8,0,0)

! EXTRA 1D vertical time serie outputs at specific locations:
! File is created in the RDIR/config/tmp directory
!        call my_outputs_point_vs_time('m°v°m','ijglob',757d0,124d0)
!        call my_outputs_point_vs_time('tonkin','lonlat',106.947173d0,19.219477d0)
!        call my_outputs_point_vs_time('station1','ijglob',77d0,525d0) !  attention A restpecter le format d0
!        call my_outputs_point_vs_time('station2','ijglob',94d0,331d0) !  attention A restpecter le format d0

! Simulations de Simon:
!         call my_outputs_point_vs_time('plateau','lonlat',-5.179655d0,47.138156d0) !  attention A restpecter le format d0
!         call my_outputs_point_vs_time('talus','lonlat',-5.643606d0,46.757172d0) !  attention A restpecter le format d0
!         call my_outputs_point_vs_time('plaine','lonlat',-5.709039d0,46.067091d0) !  attention A restpecter le format d0

#ifdef bidon
          call my_outputs_point_vs_time('st1','ijglob', 50d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st2','ijglob',100d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st3','ijglob',200d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st4','ijglob',300d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st5','ijglob',400d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st6','ijglob',500d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st7','ijglob',600d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st8','ijglob',700d0,586d0) !  attention A restpecter le format d0
          call my_outputs_point_vs_time('st9','ijglob',800d0,586d0) !  attention A restpecter le format d0
#endif

#ifdef bidon
         call my_outputs_point_vs_time('station1_vanuc','ijglob',70d0,460d0) !attention A restpecter le format d0
         call my_outputs_point_vs_time('station2_vanuc','ijglob',50d0,460d0) !attention A restpecter le format d0
         call my_outputs_point_vs_time('station3_vanuc','ijglob',37d0,460d0) !attention A restpecter le format d0
         call my_outputs_point_vs_time('stationTT_vanuc','ijglob',35d0,460d0)!attention A restpecter le format d0
         call my_outputs_point_vs_time('station_luoc_vanuc','ijglob',35d0,465d0)!attention A restpecter le format d0
         call my_outputs_point_vs_time('station_thaibinh','ijglob',35d0,470d0)!attention A restpecter le format d0
         call my_outputs_point_vs_time('station_luoc_upstream','ijglob',30d0,470d0)!attention A restpecter le format d0
#endif


#ifdef bidon
! 3 choix: ijglob, lonlat, latlon
        call my_outputs_point_vs_time('julio' ,'lonlat',5.255d0,43.135d0)
        call my_outputs_point_vs_time('solemio','lonlat',5.289d0,43.237d0)
        call my_outputs_point_vs_time('solemio2','lonlat',5.29167d0,43.2417d0)
        call my_outputs_point_vs_time('lion'   ,'lonlat',4.64d0 ,42.06d0)
        call my_outputs_point_vs_time('azur'   ,'lonlat',7.80d0,43.40d0)
        call my_outputs_point_vs_time('cassis' ,'lonlat',5.541596d0,43.146033d0)
        call my_outputs_point_vs_time('planier','lonlat',5.239798d0,43.202179d0)
        call my_outputs_point_vs_time('agde' ,'lonlat',3.537795d0,43.253423d0)
        call my_outputs_point_vs_time('riou' ,'lonlat',5.391274d0,43.17199d0)
        call my_outputs_point_vs_time('cros' ,'lonlat',6.378d0,43.0164d0)
        call my_outputs_point_vs_time('mejean' ,'lonlat',5.202154d0,43.31936d0)
        call my_outputs_point_vs_time('mesurho','lonlat',4.8662d0,43.3189d0)
        call my_outputs_point_vs_time('espiguette','lonlat',4.1333d0,43.425d0)
!        call my_outputs_point_vs_time('sete','lonlat',3.7797d0,43.3715d0)
        call my_outputs_point_vs_time('leucate','lonlat',3.1245d0,42.9162d0)
        call my_outputs_point_vs_time('sete','lonlat',3.7797d0,43.3715d0)
        call my_outputs_point_vs_time('sete2','lonlat',3.66167d0,43.3267d0)
        call my_outputs_point_vs_time('banyuls','lonlat',3.3278d0,42.5667d0)
        call my_outputs_point_vs_time('banyuls2','lonlat',3.145d0,42.4883d0)
        call my_outputs_point_vs_time('villefranche','lonlat',7.31667d0,43.6833d0)
!       call my_outputs_point_vs_time('lacazeduthiers','lonlat',3.3278d0,42.5667d0)
        call my_outputs_point_vs_time('lacazeduthiers','lonlat',3.55d0,42.4333d0)
        call my_outputs_point_vs_time('antares','lonlat',6.1667d0,42.8d0)
        call my_outputs_point_vs_time('mola','lonlat',3.533d0,42.45d0)
!!     call my_outputs_point_vs_time('planier','lonlat',5.239798d0,43.202179d0)
       call my_outputs_point_vs_time('capcreus300','lonlat',3.3278d0,42.5667d0)
       call my_outputs_point_vs_time('capcreus3002','lonlat',3.3217d0,42.39d0)
       call my_outputs_point_vs_time('CCC1000','lonlat',3.577d0,42.306d0)
       call my_outputs_point_vs_time('poem','lonlat',3.0669d0,42.7035d0)

#endif

#ifdef bidon
! SERIE PEACETIME
! 3 choix: ijglob, lonlat, latlon
      call my_outputs_point_vs_time('TMC_01_PEACETIME_ST01','latlon',41.8918d0, 6.3333d0)  !   May  12  2017    11:21:39  2500.337
      call my_outputs_point_vs_time('TMC_02_PEACETIME_ST02','latlon',40.5062d0, 6.7298d0)  !   May  13  2017    05:40:49  2699.025
      call my_outputs_point_vs_time('TMC_03_PEACETIME_ST03','latlon',39.1333d0, 7.6835d0)  !   May  14  2017    05:27:48  1249.154
      call my_outputs_point_vs_time('TMC_04_PEACETIME_ST03','latlon',39.1333d0, 7.6835d0)  !   May  14  2017    10:05:03  1247.304
      call my_outputs_point_vs_time('TMC_05_PEACETIME_ST04','latlon',37.9832d0, 7.9768d0)  !   May  15  2017    05:09:44  2664.976
      call my_outputs_point_vs_time('TMC_06_PEACETIME_ST05','latlon',38.9532d0, 11.0233d0) !    May  16  2017    02:04:26  2274.955
      call my_outputs_point_vs_time('TMC_07_PEACETIME_TYR','latlon',39.34d0, 12.5928d0)   ! May  17  2017    04:09:28  202.945
      call my_outputs_point_vs_time('TMC_08_PEACETIME_TYR','latlon',39.3403d0, 12.5943d0) !   May  17  2017    07:36:00  3299.755
      call my_outputs_point_vs_time('TMC_09_PEACETIME_TYR','latlon',39.3397d0, 12.5927d0) !   May  19  2017    03:53:36  20.121
      call my_outputs_point_vs_time('TMC_10_PEACETIME_ST06','latlon',38.8077d0, 14.4997d0) !   May  22  2017    03:48:48  2098.994
      call my_outputs_point_vs_time('TMC_11_PEACETIME_ST07','latlon',36.6582d0, 18.1548d0) !   May  23  2017    19:11:03  3506.335
      call my_outputs_point_vs_time('TMC_12_PEACETIME_ION','latlon',35.4892d0, 19.7762d0) !   May  25  2017    03:28:50  3057.979
      call my_outputs_point_vs_time('TMC_13_PEACETIME_ION','latlon',35.4892d0, 19.7762d0) !   May  25  2017    13:43:28  203.44
      call my_outputs_point_vs_time('TMC_14_PEACETIME_ION','latlon',35.4892d0, 19.7765d0) !   May  27  2017    06:24:17  205.405
      call my_outputs_point_vs_time('TMC_15_PEACETIME_ION','latlon',35.3607d0, 20.0088d0) !   May  29  2017    07:19:26  200.193
      call my_outputs_point_vs_time('TMC_16_PEACETIME_ST08','latlon',36.2103d0, 16.631d0)  !  May  30  2017    02:48:43  3000.114
      call my_outputs_point_vs_time('TMC_17_PEACETIME_ST09','latlon',38.1347d0, 5.8408d0)  !  Jun  01  2017    19:15:56  2758.184
      call my_outputs_point_vs_time('TMC_18_PEACETIME_FAST','latlon',37.946d0, 2.902d0)   ! Jun  02  2017    17:16:27  197.17
      call my_outputs_point_vs_time('TMC_19_PEACETIME_FAST','latlon',37.947d0, 2.9153d0)  !  Jun  03  2017    03:33:42  404.418
      call my_outputs_point_vs_time('TMC_20_PEACETIME_FAST','latlon',37.947d0, 2.9153d0)  !  Jun  03  2017    13:30:22  2700.635
      call my_outputs_point_vs_time('TMC_21_PEACETIME_FAST','latlon',37.947d0, 2.9153d0)  !  Jun  04  2017    18:57:18  197.476
      call my_outputs_point_vs_time('TMC_22_PEACETIME_FAST','latlon',37.9468d0, 2.9165d0) !   Jun  05  2017    06:54:33  200.916
      call my_outputs_point_vs_time('TMC_23_PEACETIME_FAST','latlon',37.9465d0, 2.9168d0) !   Jun  05  2017    13:50:30  404.095
      call my_outputs_point_vs_time('TMC_24_PEACETIME_FAST','latlon',37.9467d0, 2.9168d0) !   Jun  06  2017    02:23:26  400.022
      call my_outputs_point_vs_time('TMC_25_PEACETIME_FAST','latlon',37.9465d0, 2.9167d0) !   Jun  07  2017    03:38:54  403.521
      call my_outputs_point_vs_time('TMC_26_PEACETIME_ST10','latlon',37.457d0, 1.567d0)   ! Jun  08  2017    04:55:39  2000.673
      call my_outputs_point_vs_time('TMC_27_PEACETIME_FAST','latlon',37.9468d0, 2.916d0)  !  Jun  08  2017    20:09:30  401.58
#endif

#ifdef bidon
! Caroline Comby
       call my_outputs_point_vs_time('point01' ,'lonlat',5.2550d0,43.135d0)
       call my_outputs_point_vs_time('point02' ,'lonlat',7.9240d0,43.627d0)
       call my_outputs_point_vs_time('point03' ,'lonlat',8.2500d0,43.6600d0)
       call my_outputs_point_vs_time('point04' ,'lonlat',8.2500d0,43.5216d0)
       call my_outputs_point_vs_time('point05' ,'lonlat',8.3363d0,43.8718d0)
       call my_outputs_point_vs_time('point06' ,'lonlat',8.3984d0,43.7524d0)
       call my_outputs_point_vs_time('point07' ,'lonlat',8.4682d0,43.5216d0)
       call my_outputs_point_vs_time('point08' ,'lonlat',8.6500d0,43.9915d0)
       call my_outputs_point_vs_time('point09' ,'lonlat',8.6500d0,43.7524d0)
       call my_outputs_point_vs_time('point10' ,'lonlat',8.6500d0,43.5216d0)
       call my_outputs_point_vs_time('point11' ,'lonlat',8.6500d0,43.0000d0)
       call my_outputs_point_vs_time('point12' ,'lonlat',8.8805d0,43.9915d0)
       call my_outputs_point_vs_time('point13' ,'lonlat',8.8805d0,43.7432d0)
       call my_outputs_point_vs_time('point14' ,'lonlat',9.2102d0,43.7432d0)
#endif


! Mar Menor
      call my_outputs_point_vs_time('pointA' ,'lonlat',-0.783317d0,37.792949d0)
      call my_outputs_point_vs_time('pointB' ,'lonlat',-0.783191d0,37.698050d0)
      call my_outputs_point_vs_time('pointC' ,'lonlat',-0.752359d0,37.667696d0)
      call my_outputs_point_vs_time('merMed1' ,'lonlat',-0.659274d0,37.72370d0)
      call my_outputs_point_vs_time('merMed2' ,'lonlat',-0.466717d0,37.73270d0)


#ifdef bidon
      call my_outputs_glider
#endif


! TEST:
!      call latlontoij(lon_t(10,10),lat_t(10,10),'loc') ! returns deci, decj
!      write(10+par%rank,*)deci,decj

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!      stop 'cocote'
!#endif
        endif                          !sssssss>

        call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
        if(k0/=0)stop 'flag_stop=1 in my_outputs_driver' !26-03-19


      end subroutine my_outputs_driver

!.............................................................
#ifdef bidon
      subroutine my_outputs_read_point_vs_time(name_,lon_,lat_    &
      ,ye_,mo_,da_,ho_,mi_,se_)
      implicit none
      character(len=*) name_
      double precision deci_,decj_,deltlon_,deltlat_,dlon_di_,dlon_dj_ &
                      ,dlat_di_,dlat_dj_,lon_,lat_
      integer :: ncid_,dimtime_,year_,month_,day_,hour_,minute_,second_ &
             ,loop_,dimx_,dimy_,time_id_,ye_,mo_,da_,ho_,mi_,se_        &
             ,count_=0,var_dims,tabdim(4),var_type
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='my_outputs_read_point_vs_time'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      write(6,*)'*************************************'
      status=nf_open(trim(name_),nf_write,ncid_)
      if(status/=0)stop 'Err 74 nf_open'

! Dimensions horizontale
      status=nf_inq_dimid(ncid_,'ni_t',dim_x_id)
      if(status==0) then !>>>>
        status=nf_inq_dimlen(ncid_,dim_x_id,dimx_)
        if(status/=0)stop 'Err 75 nf_inq_dimlen x'
        write(6,*)'dimx_=',dimx_
      endif              !>>>>

      status=nf_inq_dimid(ncid_,'nj_t',dim_y_id)
      if(status==0) then !>>>>
        status=nf_inq_dimlen(ncid_,dim_y_id,dimy_)
        if(status/=0)stop 'Err 75 nf_inq_dimlen y'
        write(6,*)'dimy_=',dimy_
      endif              !>>>>

! Dimension verticale
      status=nf_inq_dimid(ncid_,'nk_t',dim_z_id)
      if(status/=0)stop 'Err 75 nf_inq_dimid z'
      status=nf_inq_dimlen(ncid_,dim_z_id,k0)
      if(status/=0)stop 'Err 75 nf_inq_dimlen z'
      if(k0/=kmax) then
       write(6,*)'k0 kmax',k0,kmax
       stop 'Err 82 k0/=kmax'
      endif

! Charger depth_t
      status=nf_inq_varid(ncid_,'depth_t',var_id)
      if(status/=0)stop 'Err 103 nf_inq_varid depth_t'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_dims==3) then !>>>>>
       varstart(1)=1 ; varstart(2)=1 ; varstart(3)=1
       varcount(1)=2 ; varcount(2)=2 ; varcount(3)=kmax
       i1=1 ; i2=2 ; j1=1 ; j2=2
      else                 !>>>>>
          if(var_dims==1) then !--->
            varstart(1)=1 ; varcount(1)=kmax ; i1=1 ; i2=1 ; j1=1 ; j2=1
          else                 !--->
           stop 'Err 121 var_dims'
          endif                !--->
      endif

      if(var_type==nf_double) then !r8r8r8>
           status=nf_get_vara_double(ncid_,var_id     &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                              ,anyv3d(i1:i2,j1:j2,1:kmax,1))
                             anyvar3d(i1:i2,j1:j2,1:kmax) &
                              =anyv3d(i1:i2,j1:j2,1:kmax,1)
      else                         !r8r8r8>
             status=nf_get_vara_real(ncid_,var_id     &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                              ,anyvar3d(i1:i2,j1:j2,1:kmax))
      endif                        !r8r8r8>
      if(status/=0)stop 'Err 105 get depth_t'
      do k=1,kmax ; do j=j1,j2 ; do i=i1,i2
       depth_t(i,j,k)=anyvar3d(i,j,k)
      enddo ; enddo ; enddo

! Charger depth_w
      status=nf_inq_varid(ncid_,'depth_w',var_id)
      if(status/=0)stop 'Err 103 nf_inq_varid depth_w'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_dims==3) then !>>>>>
       varstart(1)=1 ; varstart(2)=1 ; varstart(3)=1
       varcount(1)=2 ; varcount(2)=2 ; varcount(3)=kmax+1
       i1=1 ; i2=2 ; j1=1 ; j2=2
      else                 !>>>>>
          if(var_dims==1) then !--->
            varstart(1)=1 ; varcount(1)=kmax+1 ; i1=1 ; i2=1 ; j1=1 ; j2=1
          else                 !--->
           stop 'Err 121 var_dims'
          endif                !--->
      endif

      if(var_type==nf_double) then !r8r8r8>
           status=nf_get_vara_double(ncid_,var_id     &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                              ,anyv3d(i1:i2,j1:j2,1:kmax+1,1))
                             anyvar3d(i1:i2,j1:j2,1:kmax+1) &
                              =anyv3d(i1:i2,j1:j2,1:kmax+1,1)
      else                         !r8r8r8>
             status=nf_get_vara_real(ncid_,var_id     &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                              ,anyvar3d(i1:i2,j1:j2,1:kmax+1))
      endif                        !r8r8r8>
      if(status/=0)stop 'Err 105 get depth_w'
      do k=1,kmax+1 ; do j=j1,j2 ; do i=i1,i2
       depth_w(i,j,k)=anyvar3d(i,j,k)
      enddo ; enddo ; enddo

      if(var_dims/=1) then !oooooooo>

! Charger longitude
      varstart(1)=1 ; varstart(2)=1 
      varcount(1)=2 ; varcount(2)=2 
      status=nf_inq_varid(ncid_,'longitude_t',var_id)
      if(status/=0)stop 'Err 103 nf_inq_varid longitude_t'
      status=nf_get_vara_double(ncid_,var_id   &
                              ,varstart(1:2)   &
                              ,varcount(1:2)   &
                              ,lon_t(1:2,1:2))
      if(status/=0)stop 'Err 105 get lon_t'

! Charger latitude
      varstart(1)=1 ; varstart(2)=1 
      varcount(1)=2 ; varcount(2)=2 
      status=nf_inq_varid(ncid_,'latitude_t',var_id)
      if(status/=0)stop 'Err 103 nf_inq_varid latitude_t'
      status=nf_get_vara_double(ncid_,var_id   &
                              ,varstart(1:2)   &
                              ,varcount(1:2)   &
                              ,lat_t(1:2,1:2))
      if(status/=0)stop 'Err 105 get lat_t'

! deci_ et decj_ sont les indices decimaux de la station encadree par
! les 4 points d'indices (1:2,1;2)
      dlon_di_=0.5*(lon_t(2,1)-lon_t(1,1)  &    ! dlon/di
                   +lon_t(2,2)-lon_t(1,2))
      dlon_dj_=0.5*(lon_t(1,2)-lon_t(1,1)  &    ! dlon/dj
                   +lon_t(2,2)-lon_t(2,1))
      deltlon_=lon_-lon_t(1,1)                  ! longitude - longitude(1,1)
!     deltlon_=lon_t(1,1)-lon_t(1,1)  ! test benchmark
      dlat_di_=0.5*(lat_t(2,1)-lat_t(1,1)  &    ! dlat/di 
                   +lat_t(2,2)-lat_t(1,2))
      dlat_dj_=0.5*(lat_t(1,2)-lat_t(1,1)  &    ! dlat/dj
                   +lat_t(2,2)-lat_t(2,1))
      deltlat_=lat_-lat_t(1,1)                  ! latitude - latitude(1,1)
!     deltlat_=lat_t(1,1)-lat_t(1,1)   ! test benchmark
      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(deltlon_<-180.)deltlon_=deltlon_+360.
      if(deltlon_> 180.)deltlon_=deltlon_-360.
      deci_=1.                                      &
           +(deltlon_*dlat_dj_-deltlat_*dlon_dj_)   &
           /(dlon_di_*dlat_dj_-dlon_dj_*dlat_di_)
      decj_=1.                                      &
           +(dlon_di_*deltlat_-dlat_di_*deltlon_)   &
           /(dlon_di_*dlat_dj_-dlon_dj_*dlat_di_)

      deci_=min(max(deci_,1.),2.)
      decj_=min(max(decj_,1.),2.)

      else                 !oooooooo>

       deci_=1. ; decj_=1.

      endif                !oooooooo>

! Dimension temps 
      status=nf_inq_dimid(ncid_,'time',dim_t_id)
      if(status/=0)stop 'Err 75 nf_inq_dimid'

      status=nf_inq_dimlen(ncid_,dim_t_id,dimtime_)
      if(status/=0)stop 'Err 75 nf_inq_dimlen'

      status=nf_inq_varid(ncid_,'time',time_id_)
      if(status/=0)stop 'Err 76 nf_inq_varid'

      status=nf_get_att_text(ncid_,time_id_,'units',texte60)
      if(status/=0)stop 'Err 77 nf_get_att_text'


         read(texte60(14:17),*)year_
         read(texte60(19:20),*)month_
         read(texte60(22:23),*)day_
         read(texte60(25:26),*)hour_
         read(texte60(28:29),*)minute_
         read(texte60(31:32),*)second_

      write(6,*)'dimtime_=',dimtime_
      write(6,'(a)')trim(texte60)

      write(6,*)'year_',year_
      write(6,*)'month_',month_
      write(6,*)'day_',day_
      write(6,*)'hour_',hour_
      write(6,*)'minute_',minute_
      write(6,*)'second_',second_
      write(6,*)'datesim',datesim(1:6,1) 
      if(datesim(1,1)/=year_.or. &
         datesim(2,1)/=month_.or. &
         datesim(3,1)/=day_.or. &
         datesim(4,1)/=hour_.or. &
         datesim(5,1)/=minute_.or. &
         datesim(6,1)/=second_) then
      stop 'Err date initiale my_outputs_read_point_vs_time'
      endif

! Interpoler les 4 colonnes de profondeur a la longitude et
! latitude de la station
      if(var_dims/=1) then !---->
       i=int(deci_) ; j=int(decj_) ; rapi=deci_-i ; rapj=decj_-j
      else                 !---->
       i=1 ; j=1 ; rapi=0. ; rapj=0.
      endif
      do k=kmax+1,1,-1
                            depth_w(1  ,3  ,k)=  &
        (1.-rapi)*(1.-rapj)*depth_w(i  ,j  ,k)   &
       +(1.-rapi)*(   rapj)*depth_w(i  ,j+1,k)   &
       +(   rapi)*(1.-rapj)*depth_w(i+1,j  ,k)   &
       +(   rapi)*(   rapj)*depth_w(i+1,j+1,k)   
                            depth_t(1  ,3  ,k)=  &
        (1.-rapi)*(1.-rapj)*depth_t(i  ,j  ,k)   &
       +(1.-rapi)*(   rapj)*depth_t(i  ,j+1,k)   &
       +(   rapi)*(1.-rapj)*depth_t(i+1,j  ,k)   &
       +(   rapi)*(   rapj)*depth_t(i+1,j+1,k)   
      enddo

! A PARTIR D'ICI COMMENCE LA BOUCLE SUR LE TEMPS
! reset pour les calculs de moyennes
      km_w=0. ; count_=0 ; epsn_w=0.
      do loop_=1,dimtime_

         status=nf_get_vara_double(ncid_,time_id_,loop_,1,x3)
         if(status/=0)stop 'Err 119 nf_get_vara_double time'
         call elapsedtime2date(x3,year_,month_,day_,hour_,minute_,second_)

! Si le temps correspond a la date demandee faire l'analyse et cie...
!     if(year_==ye_.and.month_==mo_.and.day_==da_.and.hour_==ho_.and.minute_==mi_) then !>>>>>
      if(year_==ye_.and.month_==mo_.and.day_==da_.and.hour_==ho_) then !>>>>>
!     if(year_==ye_.and.month_==mo_.and.day_==da_               ) then !>>>>>

      elapsedtime_now=x3

      write(65,*)'---------------------------------------'
      write(65,*)trim(texte60),elapsedtime_now
      write(65,*)'date:',year_,month_,day_,hour_,minute_,second_

! Lire tem_t
      status=nf_inq_varid(ncid_,'tem_t',var_id)
      if(status/=0)stop 'Err 165 get tem_t'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(var_dims==2) then !vvvvvvv>
        varstart(1)=1      ; varstart(2)=loop_
        varcount(1)=kmax ; varcount(2)=1
        i1=1 ; i2=1 ; j1=1 ; j2=1
      else                 !vvvvvvv>
       if(var_dims==4) then !...>
        varstart(1)=1 ; varstart(2)=1 ; varstart(3)=1 ; varstart(4)=loop_
        varcount(1)=2 ; varcount(2)=2 ; varcount(3)=kmax ; varcount(4)=1
        i1=1 ; i2=2 ; j1=1 ; j2=2
       else                 !...>
        stop 'Err 288 var_dims?'
       endif                !...>
      endif                !vvvvvvv>

      if(var_type==nf_double) then !r8r8r8>
         status=nf_get_vara_double(ncid_,var_id       &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                    ,anyv3d(i1:i2,j1:j2,1:kmax,1))
                   anyvar3d(i1:i2,j1:j2,1:kmax)=    &
                     anyv3d(i1:i2,j1:j2,1:kmax,1) 
      else                         !r8r8r8>
        if(var_type==nf_real) then !r4r4>
         status=nf_get_vara_real(ncid_,var_id         &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                    ,anyvar3d(i1:i2,j1:j2,1:kmax))
        else                       !r4r4>
         stop 'Err 310 var_type?'
        endif                      !r4r4>
      endif                        !r8r8r8>
      do k=kmax,1,-1 ; do j=j1,j2 ; do i=i1,i2
       if(anyvar3d(i,j,k)/=-9999.) then ; tem_t(i,j,k,1)=anyvar3d(i,j,k)
                                   else ; tem_t(i,j,k,1)=   tem_t(i,j,k+1,1) ; endif
      enddo ; enddo ; enddo

! Lire km_w
      status=nf_inq_varid(ncid_,'km_w',var_id)
      if(status/=0)stop 'Err 165 get km_w'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(var_dims==2) then !vvvvvvv>
        varstart(1)=1      ; varstart(2)=loop_
        varcount(1)=kmax+1 ; varcount(2)=1
        i1=1 ; i2=1 ; j1=1 ; j2=1
      else                 !vvvvvvv>
       if(var_dims==4) then !...>
        varstart(1)=1 ; varstart(2)=1 ; varstart(3)=1 ; varstart(4)=loop_
        varcount(1)=2 ; varcount(2)=2 ; varcount(3)=kmax+1 ; varcount(4)=1
        i1=1 ; i2=2 ; j1=1 ; j2=2
       else                 !...>
        stop 'Err 288 var_dims?'
       endif                !...>
      endif                !vvvvvvv>

      if(var_type==nf_double) then !r8r8r8>
         status=nf_get_vara_double(ncid_,var_id       &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                    ,anyv3d(i1:i2,j1:j2,1:kmax+1,1))
                   anyvar3d(i1:i2,j1:j2,1:kmax+1)=    &
                     anyv3d(i1:i2,j1:j2,1:kmax+1,1) 
      else                         !r8r8r8>
        if(var_type==nf_real) then !r4r4>
         status=nf_get_vara_real(ncid_,var_id         &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                    ,anyvar3d(i1:i2,j1:j2,1:kmax+1))
        else                       !r4r4>
         stop 'Err 310 var_type?'
        endif                      !r4r4>
      endif                        !r8r8r8>
      do k=kmax+1,1,-1 ; do j=j1,j2 ; do i=i1,i2
       if(anyvar3d(i,j,k)/=-9999.) then ; km_w(i,j,k)=anyvar3d(i,j,k)
                                   else ; km_w(i,j,k)=    km_w(i,j,k+1) ; endif
      enddo ; enddo ; enddo

! Lire epsn_w
      if(.not.allocated(epsn_w))stop 'Err 252 activer scheme k-eps'
      status=nf_inq_varid(ncid_,'epsn_w',var_id)
      if(status/=0)stop 'Err 165 get epsn_w'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(var_dims==2) then !vvvvvvv>
        varstart(1)=1      ; varstart(2)=loop_
        varcount(1)=kmax+1 ; varcount(2)=1
        i1=1 ; i2=1 ; j1=1 ; j2=1
      else                 !vvvvvvv>
       if(var_dims==4) then !...>
        varstart(1)=1 ; varstart(2)=1 ; varstart(3)=1 ; varstart(4)=loop_
        varcount(1)=2 ; varcount(2)=2 ; varcount(3)=kmax+1 ; varcount(4)=1
        i1=1 ; i2=2 ; j1=1 ; j2=2
       else                 !...>
        stop 'Err 288 var_dims?'
       endif                !...>
      endif                !vvvvvvv>

      if(var_type==nf_double) then !r8r8r8>
         status=nf_get_vara_double(ncid_,var_id       &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                    ,anyv3d(i1:i2,j1:j2,1:kmax+1,1))
                   anyvar3d(i1:i2,j1:j2,1:kmax+1)=    &
                     anyv3d(i1:i2,j1:j2,1:kmax+1,1) 
      else                         !r8r8r8>
        if(var_type==nf_real) then !r4r4>
         status=nf_get_vara_real(ncid_,var_id         &
                              ,varstart(1:var_dims)   &
                              ,varcount(1:var_dims)   &
                    ,anyvar3d(i1:i2,j1:j2,1:kmax+1))
        else                         !r4r4>
         stop 'Err 310 var_type?'
        endif                        !r4r4>
      endif                        !r8r8r8>
      do k=kmax+1,1,-1 ; do j=j1,j2 ; do i=i1,i2
       if(anyvar3d(i,j,k)/=-9999.) then ; epsn_w(i,j,k)=anyvar3d(i,j,k)
                                   else ; epsn_w(i,j,k)=    epsn_w(i,j,k+1) ; endif
      enddo ; enddo ; enddo


! Legende: km_w(1,3,:)=moyenne des 4 points km_w(1:2,1:2,:)
!          km_w(3,1,:) moyenne temporelle (horaire? jour?) de km_w(1,3,:) 
!          km_w(3,3,:) moyenne temporelle (horaire? jour?) de km_w(1,3,:)**2
       count_=count_+1
       if(var_dims==4) then !>>>>
        i=int(deci_) ; j=int(decj_) ; rapi=deci_-i ; rapj=decj_-j
       else                !>>>>
        i=1 ; j=1 ; rapi=0. ; rapj=0.
       endif               !>>>>
       do k=kmax+1,1,-1
                                km_w(1  ,3  ,k)= &
         (1.-rapi)*(1.-rapj)*   km_w(i  ,j  ,k)  &
        +(1.-rapi)*(   rapj)*   km_w(i  ,j+1,k)  &
        +(   rapi)*(1.-rapj)*   km_w(i+1,j  ,k)  &
        +(   rapi)*(   rapj)*   km_w(i+1,j+1,k)  

        km_w(3,1,k)=km_w(3,1,k)+km_w(1,3,k)
        km_w(3,3,k)=km_w(3,3,k)+km_w(1,3,k)**2

                                epsn_w(1  ,3  ,k)= &
         (1.-rapi)*(1.-rapj)*   epsn_w(i  ,j  ,k)  &
        +(1.-rapi)*(   rapj)*   epsn_w(i  ,j+1,k)  &
        +(   rapi)*(1.-rapj)*   epsn_w(i+1,j  ,k)  &
        +(   rapi)*(   rapj)*   epsn_w(i+1,j+1,k)  

        epsn_w(3,1,k)=epsn_w(3,1,k)+epsn_w(1,3,k)
        epsn_w(3,3,k)=epsn_w(3,3,k)+epsn_w(1,3,k)**2

                                tem_t(1  ,3  ,k,1)= &
         (1.-rapi)*(1.-rapj)*   tem_t(i  ,j  ,k,1)  &
        +(1.-rapi)*(   rapj)*   tem_t(i  ,j+1,k,1)  &
        +(   rapi)*(1.-rapj)*   tem_t(i+1,j  ,k,1)  &
        +(   rapi)*(   rapj)*   tem_t(i+1,j+1,k,1)  

        tem_t(3,1,k,1)=tem_t(3,1,k,1)+tem_t(1,3,k,1)
        tem_t(3,3,k,1)=tem_t(3,3,k,1)+tem_t(1,3,k,1)**2

       enddo

      endif                                            !>>>>>

      enddo ! loop_

      if(count_/=0) then !oooooo>

      do k=kmax+1,1,-1
          km_w(3,1,k) =  km_w(3,1,k)  /real(count_)
          km_w(3,3,k) =  km_w(3,3,k)  /real(count_)
        epsn_w(3,1,k) =epsn_w(3,1,k)  /real(count_)
        epsn_w(3,3,k) =epsn_w(3,3,k)  /real(count_)
         tem_t(3,1,k,1)=tem_t(3,1,k,1)/real(count_)
         tem_t(3,3,k,1)=tem_t(3,3,k,1)/real(count_)
      enddo

      status=nf_close(ncid_)

      k0=index(name_,'.nc')
!     write(texte30,'(a1,i0)')'_',ye_*10000  +mo_*100  +da_
      if(ho_<10) then
       write(texte30,'(a1,i8,a2,i5)')'_',ye_*10000  +mo_*100  +da_ &
                                   ,'_0',ho_*10000  +mi_*100  +se_
      else
       write(texte30,'(a1,i8,a1,i6)')'_',ye_*10000  +mo_*100  +da_ &
                                   ,'_' ,ho_*10000  +mi_*100  +se_
      endif
      write(6,'(a)')trim(texte30)
      open(unit=66,file=name_(1:k0-1)//trim(texte30)//'.txt')
      write(66,*)'temps ecoule en secondes',elapsedtime_now
      if(var_dims==4) then !4444444>
       write(66,'(4a)')'z / km / standard_deviation_km /'  &
                ,' epsilon / standard_deviation_epsilon'   &
                ,' / km at 4 nearest points'               &
                ,' / epsilon at 4 nearest points'
       else                 !4444444>
       write(66,'(a,a)')'z / km / standard_deviation_km /' &
                ,' epsilon / standard_deviation_epsilon'
       endif                !4444444>
      
      do k=kmax+1,1,-1
       if(km_w(3,1,k)/=0.) then !>>>>
       if(var_dims==4) then !4444444>
        write(66,'(13(1x,e14.7))')                        &
                   depth_w(1,3,k)                         & ! z
                  ,km_w(3,1,k)                            & ! <km_w> (horaire)
                  ,sqrt( km_w(3,3,k)-km_w(3,1,k)**2)      & ! sqrt( <km_w**2>-<km_w>**2 )
                  ,epsn_w(3,1,k)                          & ! <epsn_w> (horaire)
                  ,sqrt( epsn_w(3,3,k)-epsn_w(3,1,k)**2)  & ! sqrt( <epsn_w**2>-<epsn_w>**2 )
                  ,km_w(1:2,1:2,k)                        & ! km_w (4 points instannes)
                  ,epsn_w(1:2,1:2,k)                        ! epsn (4 points instannes)
       else                 !4444444>
        write(66,'(5(1x,e14.7))')                         &
                   depth_w(1,3,k)                         & ! z
                  ,km_w(3,1,k)                            & ! <km_w> (horaire)
                  ,sqrt( km_w(3,3,k)-km_w(3,1,k)**2)      & ! sqrt( <km_w**2>-<km_w>**2 )
                  ,epsn_w(3,1,k)                          & ! <epsn_w> (horaire)
                  ,sqrt( epsn_w(3,3,k)-epsn_w(3,1,k)**2)    ! sqrt( <epsn_w**2>-<epsn_w>**2 )
       endif                 !444444>
    
       endif                    !>>>>
      enddo
      close(66)
! Parametres du niveau depth_t
      open(unit=66,file=name_(1:k0-1)//'_t'//trim(texte30)//'.txt')
      write(66,*)'temps ecoule en secondes',elapsedtime_now
      write(66,'(a)')'z(m) / T(°C)'
      do k=kmax,1,-1
       if(tem_t(3,1,k,1)/=0.) then !>>>>
        write(66,'(5(1x,e14.7))')                         &
                   depth_t(1,3,k)                         & ! z
                    ,tem_t(3,1,k,1)                         ! <tem_t> (horaire)
       endif                       !>>>>
      enddo
      close(66)

      endif              !oooooo>
      
      end subroutine my_outputs_read_point_vs_time
#endif
!.............................................................
!#ifdef bidon
      subroutine my_outputs_point_vs_time(name_,locationconv_,posx_,posy_)
      implicit none
      character(len=*) name_,locationconv_
      double precision :: lon_,lat_,posx_,posy_ &
                         ,filval_=-9999.d0
      real :: filvalr4_=-9999.
      integer ii_
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='my_outputs_point_vs_time'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
! https://docs.google.com/presentation/d/1fPdwjP79p_XLN9Dn__oRP1JzQ69ptKkZZQzPBHJ7Vb4/edit#slide=id.g44be21352a_0_9

! Eventuellement si on ne veut pas archiver à chaque iteration pour 
! limiter la taille du fichier, sortir si modulo(iteration3d,N)/=0
   
      if(modulo(iteration3d,2)/=0)return

! calculer les indices du point selon convention:

      if(trim(locationconv_)=='lonlat'.or. &
         trim(locationconv_)=='latlon') then  !-----> ! convention longitude latitude
         if(trim(locationconv_)=='lonlat') then ; lon_=posx_ ; lat_=posy_ ; endif
         if(trim(locationconv_)=='latlon') then ; lon_=posy_ ; lat_=posx_ ; endif !31-03-19
       ii_=1
       call latlontoij(lon_*deg2rad,lat_*deg2rad,'loc') ! returns deci, decj
       ii_=0 ; i_=nint(deci) ; j_=nint(decj) ! strategie un point (le + proche)
!              i_= int(deci) ; j_= int(decj) ! strategie quatre points encadrants
      else                                    !----->
       if(trim(locationconv_)=='ijglob') then !/////> ! convention longitude latitude
        ii_=0 ; i_=nint(posx_)-par%timax(1) ; j_=nint(posy_)-par%tjmax(1) ! strategie un point (le plus proche)
!               i_= int(posx_)-par%timax(1) ; j_= int(posy_)-par%tjmax(1) ! strategie quatre points encadrants
       else                                   !/////>
        stop 'STOP module_my_output err on locationconv_'
       endif                                  !/////>
      endif                                   !----->

! Si les indices (locaux) sont hors proc sortir:
!     write(10+par%rank,*)deci,decj
!     i1=211-par%timax(1)
!     j1=101-par%tjmax(1)
!     if(i1>2.and.i1<imax-1.and.j1>2.and.j1<jmax-1) then
!       write(10+par%rank,*)'i1,j1=',i1,j1
!       write(10+par%rank,*)'lon lat',lon_t(i1,j1)*rad2deg,lat_t(i1,j1)*rad2deg
!       call latlontoij(lon_t(i1,j1),lat_t(i1,j1),'loc') 
!       write(10+par%rank,*)'deci decj',deci,decj
!       call latlontoij(lon_t(10,10),lat_t(10,10),'loc') 
!       write(10+par%rank,*)'test 10 10=',deci,decj
!     endif

      if(i_<2.or.i_>imax-1.or.j_<2.or.j_>jmax-1)return

      flag_stop=0
      if(mask_t(i_,j_,kmax)==0) then !m°v°m> !26-03-19
       flag_stop=1
       write(6,'(a,a,a)')'Error: station ',trim(name_) &
      ,' is located in land mask. Details in fort.xxx files'
       write(10+par%rank,'(a,a,a)')'Error: station ',trim(name_) &
                        ,' is located in land mask'
       write(10+par%rank,*)'i,j glob:',i_+par%timax(1),j_+par%tjmax(1)
      endif                          !m°v°w>


! Verification: afficher lon lat au point i_,j_ et comparer a lon,lat bouee lion:
!     write(6,*)'lon lat S26',lon_t(i_,j_)*rad2deg,lat_t(i_,j_)*rad2deg
!     write(6,*)'lon lat BL ',4.64,42.06
!     stop ' test verification'

      write(texte30,'(a,a,a,i0,a)') &
      'tmp/',trim(name_),'_',par%rank,'.nc'

! Tout d'abord tester l'existence du fichier
      status=nf_open(trim(texte30),nf_write,ncid_)
      if(status==0) then !xxx>
       flag_header_=1 ! flag_header_=1 signifie fichier deja existant
       status=nf_close(ncid_) 
      else               !xxx>
       flag_header_=0 ! flag_header_=0 signifie fichier inexistant
      endif              !xxx>
      
      if(flag_header_==0) then !000000000000000>! flag_header_=0 signifie fichier inexistant

!     write(6,*)trim(name_),h_w(i_,j_)
!     do j1=int(decj),int(decj)+1
!     do i1=int(deci),int(deci)+1
!     write(6,*)trim(name_),i1+par%timax(1),j1+par%tjmax(1),h_w(i1,j1)
!     enddo
!     enddo
! Si le fichier n'existe pas commencer par faire l'entete du fichier
      status=nf_create(trim(texte30),nf_clobber,ncid_)
!     status=nf_create('toto.nc',nf_clobber,ncid_)
      if(status/=0) then
       write(6,'(a,a)')'Err nf_create file:',trim(texte30)
       stop ' STOP in module_my_output'
      endif

!     status=nf_def_dim(ncid_,'ni_t',2,i_t_dim)
      status=nf_def_dim(ncid_,'ni_t',1+ii_,i_t_dim)
      if(status/=0)stop 'erreur nf_def_dim ni_t'

!     status=nf_def_dim(ncid_,'nj_t',2,j_t_dim)
      status=nf_def_dim(ncid_,'nj_t',1+ii_,j_t_dim)
      if(status/=0)stop 'erreur nf_def_dim nj_t'

      status=nf_def_dim(ncid_,'nk_t',kmax,k_t_dim)
      if(status/=0)stop 'erreur nf_def_dim nk_t'

      status=nf_def_dim(ncid_,'nk_w',kmax+1,k_w_dim)
      if(status/=0)stop 'erreur nf_def_dim nk_t'

      status=nf_def_dim(ncid_,'time',nf_unlimited,time_dim)
      if(status/=0)stop ' erreur nf_def_dim time_dim'

      vardim(4)=time_dim
      vardim(3)= k_t_dim 
      vardim(2)= j_t_dim 
      vardim(1)= i_t_dim 


! Dans le cas de l'etude avec Andrea on sort les 4 points encadrants le
! point de station ce qui explique que l'on ajoute 2 dimensions
! horizontales supplementaires de dimensions reduites a 2.
! Certaines variables sont a 4 dimensions

      call elapsedtime2date(0.d0,i1,i2,i3,i4,i5,i6) ! enters 0. (i.e. elapsedtime_now=0)

! Time en secondes (egalement variable dimension)
      write(texte80(2),'(a14,i4,5(a1,i2))')                    & !units
      'seconds since ',i1,'-',i2,'-',i3,' ',i4,':',i5,':',i6
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0' 
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0' 
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'
      if(texte80(2)(32:32)==' ')texte80(2)(32:32)='0'

      status=nf_def_var(ncid_,'time',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time'
      status=nf_put_att_text(ncid_,var_id,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)stop ' erreur nf_put_att_text'
      status=nf_put_att_text(ncid_,var_id,'calendar',9,'gregorian')
      if(status/=0)stop ' erreur nf_put_att_text calendar'

! Time en jours
      write(texte80(2),'(a11,i4,5(a1,i2))')                    & !units
      'days since ',i1,'-',i2,'-',i3,' ',i4,':',i5,':',i6
      if(texte80(2)(17:17)==' ')texte80(2)(17:17)='0' 
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0' 
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0' 
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'

      status=nf_def_var(ncid_,'time_days',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_days'
      status=nf_put_att_text(ncid_,var_id,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)stop ' erreur nf_put_att_text'
      status=nf_put_att_text(ncid_,var_id,'calendar',9,'gregorian')
      if(status/=0)stop ' erreur nf_put_att_text calendar'
                                                    ! returns i1,i2,...=y,m,d,h,m,s
! Time en annee decimale:
      status=nf_def_var(ncid_,'time_years',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_years'
      status=nf_put_att_text(ncid_,var_id,'long_name',13 &
      ,'Decimal_years')

! Time date yyyymmddhhmmss:
      status=nf_def_var(ncid_,'time_calendar',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_calendar'
      status=nf_put_att_text(ncid_,var_id,'long_name',19 &
      ,'time_yyyymmddhhmmss')

      status=nf_def_var(ncid_,'longitude_t',nf_double,2,vardim(1:2),var_id)
      if(status/=0)stop ' erreur nf_def_var longitude_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',9 &
      ,'longitude')
      status=nf_put_att_text(ncid_,var_id,'units',6,'degree')

      status=nf_def_var(ncid_,'latitude_t',nf_double,2,vardim(1:2),var_id)
      if(status/=0)stop ' erreur nf_def_var latitude_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',8 &
      ,'latitude')
      status=nf_put_att_text(ncid_,var_id,'units',6,'degree')

      status=nf_def_var(ncid_,'h_w',nf_real,2,vardim(1:2),var_id)
      if(status/=0)stop ' erreur nf_def_var h_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',10 &
      ,'Bathymetry')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'depth_t',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var depth_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',14 &
      ,'levels_depth_t')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'depth_w',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var depth_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',14 &
      ,'levels_depth_w')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

!     write(6,*)'vardim(1:4)',vardim(1:4)
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'tem_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var tem_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',34 &
       ,'potential_temperature_of_sea_water')
      status=nf_put_att_text(ncid_,var_id,'units',1,'C')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

!      vardim(3)=k_t_dim 
!      status=nf_def_var(ncid_,'temobc_t',nf_real,4,vardim(1:4),var_id)
!      if(status/=0)stop ' erreur nf_def_var temobc_t'
!      status=nf_put_att_text(ncid_,var_id,'long_name',30 &
!       ,'Initial_and_border_temperature')
!      status=nf_put_att_text(ncid_,var_id,'units',1,'C')

!      if(allocated(temlwf_t)) then !>>>>>
!       vardim(3)=k_t_dim 
!       status=nf_def_var(ncid_,'temlwf_t',nf_real,4,vardim(1:4),var_id)
!       if(status/=0)stop ' erreur nf_def_var temlwf_t'
!       status=nf_put_att_text(ncid_,var_id,'long_name',33 &
!        ,'low_frequency_temperature_anomaly') !02-04-16
!       status=nf_put_att_text(ncid_,var_id,'units',1,'C')
!      endif                        !>>>>>

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'sal_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var sal_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',21 &
       ,'salinity_of_sea_water')
      status=nf_put_att_text(ncid_,var_id,'units',3,'PSU')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

!      vardim(3)=k_t_dim 
!      status=nf_def_var(ncid_,'salobc_t',nf_real,4,vardim(1:4),var_id)
!      if(status/=0)stop ' erreur nf_def_var tem_t'
!      status=nf_put_att_text(ncid_,var_id,'long_name',27 &
!       ,'Initial_and_border_salinity')
!      status=nf_put_att_text(ncid_,var_id,'units',3,'PSU')

!      if(allocated(sallwf_t)) then !>>>>>
!       vardim(3)=k_t_dim 
!       status=nf_def_var(ncid_,'sallwf_t',nf_real,4,vardim(1:4),var_id)
!       if(status/=0)stop ' erreur nf_def_var sallwf_t'
!       status=nf_put_att_text(ncid_,var_id,'long_name',30 &
!        ,'low_frequency_salinity_anomaly')
!       status=nf_put_att_text(ncid_,var_id,'units',3,'PSU')
!      endif                        !>>>>>

      vardim(3)=k_t_dim !15-02-16
      status=nf_def_var(ncid_,'rhp_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var rhp_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',17 &
       ,'Potential_density')
      status=nf_put_att_text(ncid_,var_id,'units',7,'kg/m**3')

      
!     if(allocated(bio_t)) then !>>>>>
!      vardim(3)=k_t_dim 
!      status=nf_def_var(ncid_,'tracer1',nf_real,4,vardim(1:4),var_id)
!      if(status/=0)stop ' erreur nf_def_var tracer1'
!      status=nf_put_att_text(ncid_,var_id,'long_name',16 &
!       ,'passive_tracer_1')
!      status=nf_put_att_text(ncid_,var_id,'units',9,'undefined')
!     endif                        !>>>>>

!      if(allocated(bio_t)) then ! bio_t !05-09-19
!        vardim(3)=k_t_dim
!        texte80(10)='none'
!        texte80(2)='none'         ! variable ; units
!        do vb=1,vbmax
!          write(texte80(1),'(a,i0)')'bio',vb
!          write(texte80(3),'(a,i0)')'biogeochemical_variable_number_',vb
!          status=nf_def_var(ncid_,texte80(1),nf_real,4,vardim(1:4),var_id)
!          if(status/=0)stop ' erreur nf_def_var BIO_T'
!          status=nf_put_att_text(ncid_,var_id,'long_name',16 &
!                                                ,texte80(3))
!          status=nf_put_att_text(ncid_,var_id,'units',9,'undefined')
!        enddo
!      endif                        ! bio_t

! BIO
     if(imodelbio==1) then !bbbbbbbbbbbbbb>

! INITRATE
     vardim(3)=k_t_dim 
     status=nf_def_var(ncid_,'INITRATE',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var INITRATE'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_nitrate')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! ICHLPICO
     vardim(3)=k_t_dim 
     status=nf_def_var(ncid_,'ICHLPPICO',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var ICHLPPICO'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_chlphypico')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')

! ICHLNANO
     vardim(3)=k_t_dim 
     status=nf_def_var(ncid_,'ICHLPNANO',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var ICHLPNANO'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_chlphynano')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')

! ICHLMICRO
     vardim(3)=k_t_dim
     status=nf_def_var(ncid_,'ICHLPMICRO',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var ICHLPMICRO'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_chlphymicro')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')


! ICPICO
     vardim(3)=k_t_dim 
     status=nf_def_var(ncid_,'ICPPICO',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var ICPPICO'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_cphypico')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! ICNANO
     vardim(3)=k_t_dim 
     status=nf_def_var(ncid_,'ICPNANO',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var ICPNANO'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_cphynano')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! ICMICRO
     vardim(3)=k_t_dim
     status=nf_def_var(ncid_,'ICPMICRO',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var ICPMICRO'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_cphymicro')
     status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! IDO
     vardim(3)=k_t_dim
     status=nf_def_var(ncid_,'IOXYGEN',nf_real,4,vardim(1:4),var_id)
     if(status/=0)stop ' erreur nf_def_var IOXYGEN'
     status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'concentration_oxygen')
     status=nf_put_att_text(ncid_,var_id,'units',8,'µmol.kg-1')


     endif                 !bbbbbbbbbbbbbb>

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'Veastward',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var Veastward'
      status=nf_put_att_text(ncid_,var_id,'long_name',17 &
       ,'Eastward_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'Vnorthward',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var Vnorthward'
      status=nf_put_att_text(ncid_,var_id,'long_name',18 &
       ,'Northward_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'vel_u',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var vel_u'
      status=nf_put_att_text(ncid_,var_id,'long_name',15 & !13-12-18
       ,'i_axis_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'vel_v',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var vel_v'
      status=nf_put_att_text(ncid_,var_id,'long_name',15 & !13-12-18
       ,'j_axis_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

!      vardim(3)=k_t_dim 
!      status=nf_def_var(ncid_,'V_OGCM_eastward',nf_real,4,vardim(1:4),var_id)
!      if(status/=0)stop ' erreur nf_def_var VOGCM_eastward'
!      status=nf_put_att_text(ncid_,var_id,'long_name',17 &
!       ,'OGCM_eastward_velocity')
!      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
!      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

!      vardim(3)=k_t_dim 
!      status=nf_def_var(ncid_,'V_OGCM_northward',nf_real,4,vardim(1:4),var_id)
!      if(status/=0)stop ' erreur nf_def_var VOGCM_northward'
!      status=nf_put_att_text(ncid_,var_id,'long_name',18 &
!       ,'OGCM_northward_velocity')
!      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
!      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'N2',nf_real,4,vardim(1:4),var_id)
! Definition de N2: !22-04-16
! https://docs.google.com/document/d/1BKPwtP5EYh7q_iz5xHUPfai06GCnGtBKuqdl07rakig/edit
      if(status/=0)stop ' erreur nf_def_var N2'
      status=nf_put_att_text(ncid_,var_id,'long_name',30 &
      ,'brunt_vaisala_frequency_square')
      status=nf_put_att_text(ncid_,var_id,'units',5,'s**-2')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'km_w',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var km_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',59 &
      ,'molecular+turbulent_vertical_mixing_coef_momentum_equations')
      status=nf_put_att_text(ncid_,var_id,'units',4,'m2/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)
      status=nf_put_att_real(ncid_,var_id,'Molecular viscosity',nf_real,1,kmol_m)

      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'kh_w',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var kh_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',38 &
      ,'turbulent_vertical_mixing_coef_tracers')
      status=nf_put_att_text(ncid_,var_id,'units',4,'m2/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)
      status=nf_put_att_real(ncid_,var_id,'Heat molecular diffusivity',nf_real,1,kmol_h)
      status=nf_put_att_real(ncid_,var_id,'Salt molecular diffusivity',nf_real,1,kmol_s)

      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'tken_w',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var tken_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',24 &
      ,'turbulent_kinetic_energy')
      status=nf_put_att_text(ncid_,var_id,'units',5,'m2/s2')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

!ctke2*sqrt(tken_w(i,j,k))/(tkle_w(i,j,k)+small2))
      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'epsn_w',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var eps'
      status=nf_put_att_text(ncid_,var_id,'long_name',24 &
      ,'tke_dissipation(epsilon)')
      status=nf_put_att_text(ncid_,var_id,'units',5,'m2/s3')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

! omega_w
      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'omega_w',nf_real,4,vardim(1:4),var_id) !20-05-18
      if(status/=0)stop ' erreur omega_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',26 & 
      ,'relative_vertical_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'Weastward',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var Weastward'
      status=nf_put_att_text(ncid_,var_id,'long_name',13 &
      ,'Eastward_Wind')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'Wnorthward',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var Wnorthward'
      status=nf_put_att_text(ncid_,var_id,'long_name',14 &
      ,'Northward_Wind')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')


      if(flag_abl==1) then !AAAAAAAAAAA>

      status=nf_def_var(ncid_,'teta2_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',31 &
        ,'potential_air_temperature_at_2m')
      status=nf_put_att_text(ncid_,var_id,'units',1,'C')

      status=nf_def_var(ncid_,'teta0_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',41 &
       ,'temperature_at_0m_in_meteorological_model')
      status=nf_put_att_text(ncid_,var_id,'units',1,'C')

      status=nf_def_var(ncid_,'teta2delta_t',nf_double,3,vardim(1:3),var_id) !30-11-15
      status=nf_put_att_text(ncid_,var_id,'long_name',39 &
        ,'potential_air_temperature_anomaly_at_2m')
      status=nf_put_att_text(ncid_,var_id,'units',1,'K')

      status=nf_def_var(ncid_,'q2_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',23 &
        ,'specific_humidity_at_2m')

      status=nf_def_var(ncid_,'q2delta_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',31 &
        ,'specific_humidity_anomaly_at_2m')

      if(flag_abl2==1) then !222222222>

      status=nf_def_var(ncid_,'uwind_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',9 &
        ,'10m_uwind')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')

      status=nf_def_var(ncid_,'uwinddelta_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',18 &
        ,'10m_uwind_anomalie')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')

      status=nf_def_var(ncid_,'vwind_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',9 &
        ,'10m_vwind')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')

      status=nf_def_var(ncid_,'vwinddelta_t',nf_double,3,vardim(1:3),var_id)
      status=nf_put_att_text(ncid_,var_id,'long_name',18 &
        ,'10m_vwind_anomalie')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')

      endif                 !222222222>

      endif                !AAAAAAAAAAA>

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'wstress_u',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var wstress_u'
      status=nf_put_att_text(ncid_,var_id,'long_name',20 &
      ,'Eastward_Wind_Stress')
      status=nf_put_att_text(ncid_,var_id,'units',6,'N/m**2')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'wstress_v',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var wstress_v'
      status=nf_put_att_text(ncid_,var_id,'long_name',21 &
      ,'Northward_Wind_Stress')
      status=nf_put_att_text(ncid_,var_id,'units',6,'N/m**2')

!      vardim(3)=time_dim
!      status=nf_def_var(ncid_,'Oi2Dflux',nf_real,3,vardim(1:3),var_id)
!      if(status/=0)stop ' erreur nf_def_var Oi2Dflux'
!      status=nf_put_att_text(ncid_,var_id,'long_name',29 &
!      ,'Along_Oi_axis_z_averaged_flux')
!      status=nf_put_att_text(ncid_,var_id,'units',6,'m**3/s')

!      vardim(3)=time_dim
!      status=nf_def_var(ncid_,'Oj2Dflux',nf_real,3,vardim(1:3),var_id)
!      if(status/=0)stop ' erreur nf_def_var Oj2Dflux'
!      status=nf_put_att_text(ncid_,var_id,'long_name',29 &
!      ,'Along_Oj_axis_z_averaged_flux')
!      status=nf_put_att_text(ncid_,var_id,'units',6,'m**3/s')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'ssr',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var ssr'
      status=nf_put_att_text(ncid_,var_id,'long_name',23 &
      ,'surface_solar_radiation')
      status=nf_put_att_text(ncid_,var_id,'units',6,'W/m**2')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'snsf',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var snsf'
      status=nf_put_att_text(ncid_,var_id,'long_name',22 &
      ,'surface_non_solar_flux')
      status=nf_put_att_text(ncid_,var_id,'units',6,'W/m**2')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'slhf',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var slhf'
      status=nf_put_att_text(ncid_,var_id,'long_name',24 &
      ,'surface_latent_heat_flux')
      status=nf_put_att_text(ncid_,var_id,'units',6,'W/m**2')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'sshf',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var sshf'
      status=nf_put_att_text(ncid_,var_id,'long_name',26 &
      ,'surface_sensible_heat_flux')
      status=nf_put_att_text(ncid_,var_id,'units',6,'W/m**2')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'precipitation',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var precipitation'
      status=nf_put_att_text(ncid_,var_id,'long_name',13 &
      ,'precipitation')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s') !27-03-19

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'ssh',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var ssh'
      status=nf_put_att_text(ncid_,var_id,'long_name',18 &
      ,'sea_surface_height')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

! Global attributs:
      k10=len(trim(name_))
      status=nf_put_att_text(ncid_,nf_global,'Point',k10,trim(name_))
      k10=len(trim(model_name))
      status=nf_put_att_text(ncid_,nf_global,'Model',k10,trim(model_name))
      write(texte30,'(a,i0,a,i0)')'global indexes: ',i_+par%timax(1),' ',j_+par%tjmax(1)
      k10=len(trim(texte30))
      status=nf_put_att_text(ncid_,nf_global,'Position',k10,trim(texte30)) !10-10-16

!     status=nf_put_att_double(ncid_,nf_global,'bathymetry',nf_double,1,h_w(i_,j_))
!     x1=lon_t(i_,j_)*rad2deg ; x2=lat_t(i_,j_)*rad2deg 
!     status=nf_put_att_double(ncid_,nf_global,'longitude',nf_double,1,x1)
!     status=nf_put_att_double(ncid_,nf_global,'latitude' ,nf_double,1,x2)

! 194 continue
      status=nf_enddef(ncid_)

      if(status/=0)stop 'Err module_my_outputs nf_enddef'

! Ecrire longitude_t
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyv3d(i,j,1,1)=lon_t(i,j)*rad2deg
      enddo        ; enddo
      status=nf_inq_varid(ncid_,'longitude_t',var_id)
      if(status/=0)stop 'erreur nf_inq_varid longitude_t'
      status=nf_put_vara_double(ncid_,var_id        &
                                     ,varstart(1:2) &
                                     ,varcount(1:2) &
                       ,anyv3d(i_:i_+ii_,j_:j_+ii_,1,1))
      if(status/=0)stop 'erreur nf_put_var_double longitude_t'

! Ecrire latitude_t
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyv3d(i,j,1,1)=lat_t(i,j)*rad2deg
      enddo        ; enddo
      status=nf_inq_varid(ncid_,'latitude_t',var_id)
      if(status/=0)stop 'erreur nf_inq_varid latitude_t'
      status=nf_put_vara_double(ncid_,var_id        &
                                     ,varstart(1:2) &
                                     ,varcount(1:2) &
                       ,anyv3d(i_:i_+ii_,j_:j_+ii_,1,1))
      if(status/=0)stop 'erreur nf_put_var_double latitude_t'

! Ecrire h_w
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar2d(i,j)=h_w(i,j)
      enddo        ; enddo
      status=nf_inq_varid(ncid_,'h_w',var_id)
      if(status/=0)stop 'erreur nf_inq_varid h_w'
      status=nf_put_vara_real(ncid_,var_id        &
                                   ,varstart(1:2) &
                                   ,varcount(1:2) &
                       ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'erreur nf_put_var_double h_w'

! Ecrire depth_t
      varstart(3)=1    ; varcount(3)=kmax
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,k)=depth_t(i,j,k)
      enddo       ; enddo        ; enddo
      status=nf_inq_varid(ncid_,'depth_t',var_id)
      if(status/=0)stop 'erreur nf_inq_varid depth_t'
      status=nf_put_vara_real(ncid_,var_id        &
                                   ,varstart(1:3) &
                                   ,varcount(1:3) &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'erreur nf_put_var_double depth_t'

! Ecrire depth_w
      varstart(3)=1    ; varcount(3)=kmax+1
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,k)=depth_w(i,j,k)
      enddo       ; enddo        ; enddo
      status=nf_inq_varid(ncid_,'depth_w',var_id)
      if(status/=0)stop 'erreur nf_inq_varid depth_w'
      status=nf_put_vara_real(ncid_,var_id        &
                                   ,varstart(1:3) &
                                   ,varcount(1:3) &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'erreur nf_put_var_double depth_w'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err module_my_outputs nf_close'

      endif                    !000000000000000>! flag_header_=0 signifie fichier inexistant

! Ecrire dans le fichier netcdf (en ajoutant au fichier existent au fur et a mesure que la
! simu avance)

      write(texte30,'(a,a,a,i0,a)') &
      'tmp/',trim(name_),'_',par%rank,'.nc'

      status=nf_open(trim(texte30),nf_write,ncid_) ! noter l'option "nf_write"
!     status=nf_open('toto.nc',nf_write,ncid_) ! noter l'option "nf_write"
      if(status/=0) then
       write(6,'(a,a)')'Err nf_open file:',trim(texte30)
       stop ' STOP in module_my_output'
      endif

      status=nf_inq_dimid(ncid_,'time',dim_t_id)
      if(status==0)status=nf_inq_dimlen(ncid_,dim_t_id,time_counter_)

      time_counter_=time_counter_+1

!...................
! Ecrire le temps en secondes (variable dimension):
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time'
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,elapsedtime_now)
      if(status/=0)stop 'Erreur nf_put_vara_double time'

!...................
! Ecrire le temps en jours
      x1=elapsedtime_now/86400.
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time_days',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time_days'
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,x1)
      if(status/=0)stop 'Erreur nf_put_vara_double time_days'

!...................
! Ecrire le temps en annees decimales: 
      x1=year_now                                                     &
       +(elapsedtime_now               -seconds_since_1jan(year_now)) &
       /(seconds_since_1jan(year_now+1)-seconds_since_1jan(year_now))
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time_years',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time_years'
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,x1)
      if(status/=0)stop 'Erreur nf_put_vara_double time_years'

!...................
! Ecrire le temps sous forme de date yyyymmddhhmmss !30-07-15
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time_calendar',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time_calendar'
      write(date_,'(I4.4,5(I2.2))')year_now,month_now,day_now,hour_now,minute_now,second_now !      15-12-20
      read(date_,'(f14.0)') x1                                                           !      15-12-20
         !  x1=second_now             &
         !    +minute_now*10         &
         !    +hour_now  *1000       &
         !    +day_now   *100000     &
         !    +month_now *10000000   &
         !    +year_now  *1000000000
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,x1)
      if(status/=0)stop 'Erreur nf_put_vara_double time_calendar'

!...................
! Ecrire le profil vertical de temperature:
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=tem_t(i,j,k,now)*mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'tem_t',var_id)
      if(status/=0)stop ' erreur nf_inq_varid depth_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real tem_t'

!!...................
!! temperature de mercator
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!      anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))                      &
!         +(     timeweightobc(trc_id) *temobc_t(i,j,k,2)               &
!           +(1.-timeweightobc(trc_id))*temobc_t(i,j,k,0))*mask_t(i,j,k)
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'temobc_t',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid temobc_t'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!
!!...................
!! temperature basse frequence
!      if(allocated(temlwf_t)) then !pmxpmxpmx>
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!       anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))+temlwf_t(i,j,k)*mask_t(i,j,k)
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'temlwf_t',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid temlwf_t'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!      endif                        !pmxpmxpmx>
!
!...................
! Ecrire le profil vertical de salinite:
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=sal_t(i,j,k,now)*mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'sal_t',var_id)
      if(status/=0)stop ' erreur nf_inq_varid sal_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real sal_t'

!!...................
!! salinite de mercator
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!      anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))                      &
!         +(     timeweightobc(trc_id) *salobc_t(i,j,k,2)               &
!           +(1.-timeweightobc(trc_id))*salobc_t(i,j,k,0))*mask_t(i,j,k)
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'salobc_t',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid salobc_t'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!
!!...................
!! salinity basse frequence
!      if(allocated(sallwf_t)) then !pmxpmxpmx> !30-03-16
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!       anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))+sallwf_t(i,j,k)*mask_t(i,j,k)
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'sallwf_t',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid sallwf_t'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!      endif                        !pmxpmxpmx>
!
!...................
! Ecrire le profil de densite potentielle
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=(rhp_t(i,j,k)+rho)*mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'rhp_t',var_id)
      if(status/=0)stop ' erreur nf_inq_varid rhp_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real rhp_t'

!...................
! traceur passif No 1
!     if(allocated(bio_t)) then !pmxpmxpmx> !27-10-16
!     varstart(4)=time_counter_  ; varcount(4)=1
!     varstart(3)=1              ; varcount(3)=kmax
!     varstart(1:2)=1            ; varcount(1:2)=1+ii_
!     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!      anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))+bio_t(i,j,k,1)*mask_t(i,j,k)
!     enddo ; enddo ; enddo
!     status=nf_inq_varid(ncid_,'tracer1',var_id)
!     if(status/=0)stop ' erreur nf_inq_varid tracer1'
!     status=nf_put_vara_real(ncid_,var_id                &
!                                  ,varstart(1:4)         &
!                                  ,varcount(1:4)         &
!                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!     endif                        !pmxpmxpmx>


!!...................
!! traceur passifs
!      if(allocated(bio_t)) then !pmxpmxpmx> !05-09-19
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do vb=1,vbmax
!       write(texte80(1),'(a,i0)')'bio',vb
!       do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
!        anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))+bio_t(i,j,k,vb)*mask_t(i,j,k)
!       enddo ; enddo ; enddo
!       status=nf_inq_varid(ncid_,texte80(1),var_id)
!       if(status/=0)stop ' erreur nf_inq_varid tracer bio'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!      enddo
!      endif                        !pmxpmxpmx>
!


!...................
     if(imodelbio==1) then !bbbbbbbbbbbbbb>
! Nitrate
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
     anyvar3d(i,j,k)=bio_t(i,j,k,INITRATE)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'INITRATE',var_id)
     if(status/=0)stop ' erreur nf_inq_varid INITRATE'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Chl PhyPico
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
     anyvar3d(i,j,k)=bio_t(i,j,k,Isynechl)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'ICHLPPICO',var_id)
     if(status/=0)stop ' erreur nf_inq_varid ICHLPPICO'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))


! Chl PhyNano
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
     anyvar3d(i,j,k)=bio_t(i,j,k,Inanochl)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'ICHLPNANO',var_id)
     if(status/=0)stop ' erreur nf_inq_varid ICHLPNANO'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Chl PhyMicro
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
     anyvar3d(i,j,k)=bio_t(i,j,k,Idiachl)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'ICHLPMICRO',var_id)
     if(status/=0)stop ' erreur nf_inq_varid ICHLPMICRO'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))


! C PhyPico
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
     anyvar3d(i,j,k)=bio_t(i,j,k,Isynec)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'ICPPICO',var_id)
     if(status/=0)stop ' erreur nf_inq_varid ICPPICO'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))


! C PhyNano
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
     anyvar3d(i,j,k)=bio_t(i,j,k,Inanoc)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'ICPNANO',var_id)
     if(status/=0)stop ' erreur nf_inq_varid ICPNANO'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! C PhyMicro
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
     anyvar3d(i,j,k)=bio_t(i,j,k,Idiac)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'ICPMICRO',var_id)
     if(status/=0)stop ' erreur nf_inq_varid ICPMICRO'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Oxygen
     varstart(4)=time_counter_  ; varcount(4)=1
     varstart(3)=1              ; varcount(3)=kmax
     varstart(1:2)=1            ; varcount(1:2)=1+ii_
     do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
     anyvar3d(i,j,k)=bio_t(i,j,k,Ioxygen)*mask_t(i,j,k) &
                             +filvalr4_*(1-mask_t(i,j,k))
     enddo ; enddo ; enddo
     status=nf_inq_varid(ncid_,'IOXYGEN',var_id)
     if(status/=0)stop ' erreur nf_inq_varid IOXYGEN'
     status=nf_put_vara_real(ncid_,var_id                &
                                  ,varstart(1:4)         &
                                  ,varcount(1:4)         &
                 ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))


   endif                 !bbbbbbbbbbbbbb>
      

!...................
! Ecrire le profil de courant OE
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
          anyvar3d(i,j,k)=0.5*(                        &
            ( vel_u(i  ,j  ,k,1)                       &
             +vel_u(i+1,j  ,k,1))*gridrotcos_t(i,j)    &
           +( vel_v(i  ,j  ,k,1)                       &
             +vel_v(i  ,j+1,k,1))*gridrotsin_t(i,j))*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'Veastward',var_id)
      if(status/=0)stop ' erreur nf_inq_varid Veastward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real Veastward'

!...................
! Ecrire le profil de courant SN
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
          anyvar3d(i,j,k)=0.5*(                        &
           -( vel_u(i  ,j  ,k,1)                       &
             +vel_u(i+1,j  ,k,1))*gridrotsin_t(i,j)    &
           +( vel_v(i  ,j  ,k,1)                       &
             +vel_v(i  ,j+1,k,1))*gridrotcos_t(i,j) )*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'Vnorthward',var_id)
      if(status/=0)stop ' erreur nf_inq_varid Vnorthward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real Vnorthward'

!...................
! Ecrire le profil de vel_u !13-12-18
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
          anyvar3d(i,j,k)=0.5*( vel_u(i  ,j  ,k,1)           &
                               +vel_u(i+1,j  ,k,1))*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'vel_u',var_id)
      if(status/=0)stop ' erreur nf_inq_varid vel_u'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real vel_u'

!...................
! Ecrire le profil de vel_v !13-12-18
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
          anyvar3d(i,j,k)=0.5*( vel_v(i  ,j  ,k,1)           &
                               +vel_v(i  ,j+1,k,1))*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'vel_v',var_id)
      if(status/=0)stop ' erreur nf_inq_varid vel_v'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real vel_v'

!!...................
!! Ecrire le profil de courant OGCM OE
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      x2=timeweightobc(vel_id) ; x0=1.-x2
!      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!          anyvar3d(i,j,k)=0.5*( &
!       x0*( ( velobc_u(i  ,j  ,k,0)                       &
!             +velobc_u(i+1,j  ,k,0))*gridrotcos_t(i,j)    &
!           +( velobc_v(i  ,j  ,k,0)                       &
!             +velobc_v(i  ,j+1,k,0))*gridrotsin_t(i,j))   &
!      +x2*( ( velobc_u(i  ,j  ,k,2)                       &
!             +velobc_u(i+1,j  ,k,2))*gridrotcos_t(i,j)    &
!           +( velobc_v(i  ,j  ,k,2)                       &
!             +velobc_v(i  ,j+1,k,2))*gridrotsin_t(i,j))   &
!                               )*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'V_OGCM_eastward',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid V_OGCM_eastward'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!      if(status/=0)stop 'Erreur nf_put_vara_real V_OGCM_eastward'
!
!!...................
!! Ecrire le profil de courant OGCM SN
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      x2=timeweightobc(vel_id) ; x0=1.-x2
!      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!          anyvar3d(i,j,k)=0.5*( &
!       x0*(-( velobc_u(i  ,j  ,k,0)                       &
!             +velobc_u(i+1,j  ,k,0))*gridrotsin_t(i,j)    &
!           +( velobc_v(i  ,j  ,k,0)                       &
!             +velobc_v(i  ,j+1,k,0))*gridrotcos_t(i,j) )  &
!      +x2*(-( velobc_u(i  ,j  ,k,2)                       &
!             +velobc_u(i+1,j  ,k,2))*gridrotsin_t(i,j)    &
!           +( velobc_v(i  ,j  ,k,2)                       &
!             +velobc_v(i  ,j+1,k,2))*gridrotcos_t(i,j) )  &
!                              )*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'V_OGCM_northward',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid V_OGCM_northward'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
!      if(status/=0)stop 'Erreur nf_put_vara_real V_OGCM_northward'

!...................

!...................
! Ecrire le profil vertical de N**2  brunt_vaisala_frequency_square  
! Definition de N2:
! https://docs.google.com/document/d/1BKPwtP5EYh7q_iz5xHUPfai06GCnGtBKuqdl07rakig/edit
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!     anyvar3d(i,j,1)=0. ; anyvar3d(i,j,kmax+1)=0.
      do k=2,kmax 
      anyvar3d(i,j,k)=-mask_t(i,j,k)*grav/rho*(rhp_t(i,j,k)  -rhp_t(i,j,k-1)) & !30-03-16
                                           /(depth_t(i,j,k)-depth_t(i,j,k-1)) & !29-01-21
         +filvalr4_*(1-mask_t(i,j,k))
      enddo
      anyvar3d(i,j,1)     =anyvar3d(i,j,2)   !19-04-16
      anyvar3d(i,j,kmax+1)=anyvar3d(i,j,kmax)
      enddo ; enddo
      status=nf_inq_varid(ncid_,'N2',var_id)
      if(status/=0)stop ' erreur nf_inq_varid N2'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real N2'

!...................
! Ecrire le profil vertical de km_w
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax+1
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!!      anyvar3d(i,j,k)=km_w(i,j,k)*mask_t(i,j,k) &
!       anyvar3d(i,j,k)=max(km_w(i,j,k),1.049e-6)*mask_t(i,j,k) & !27-03-19
!                                   +filvalr4_*(1-mask_t(i,j,k))
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'km_w',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid km_w'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
!      if(status/=0)stop 'Erreur nf_put_vara_real km_w'
!
!...................
! Ecrire le profil vertical de kh_w
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!      anyvar3d(i,j,k)=(kh_w(i,j,k)+kmol_h)*mask_t(i,j,k) &
       anyvar3d(i,j,k)=max(kh_w(i,j,k),1.49e-7)*mask_t(i,j,k) & !27-03-19
                                  +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'kh_w',var_id)
      if(status/=0)stop ' erreur nf_inq_varid kh_w'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real kh_w'

!...................
! Ecrire le profil vertical de tken_w
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax+1
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!       anyvar3d(i,j,k)=tken_w(i,j,k)*mask_t(i,j,k) &
!                       +filvalr4_*(1-mask_t(i,j,k))
!      enddo ; enddo ; enddo
 !     status=nf_inq_varid(ncid_,'tken_w',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid tken_w'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
!      if(status/=0)stop 'Erreur nf_put_vara_real tken_w'
!
!...................
! Ecrire le profil vertical de epsn_w
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax+1
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      if(iturbulence==1) then  !ooooo>
!      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!        anyvar3d(i,j,k)=epsn_w(i,j,k)*mask_t(i,j,k) &
!                        +filvalr4_*(1-mask_t(i,j,k))
!      enddo ; enddo ; enddo
!      else                     !ooooo>
!       if(iturbulence==0) then !mmmmm>
!       do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!!epsilon dans le schema de Gaspar=ctke2*(tken_w(i,j,k)**1.5)/(tkle_w(i,j,k)+small2)
!               anyvar3d(i,j,k)=        &
!        ctke2*(tken_w(i,j,k)**1.5)     &
!             /(tkle_w(i,j,k)+small2)   &
!              *mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
!       enddo ; enddo ; enddo
!       else                    !mmmmm>
!        do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_  !23-08-16
!         anyvar3d(i,j,k)=0. !26-04-16
!        enddo ; enddo ; enddo
!       endif                   !mmmmm>
!      endif                    !ooooo>
!      status=nf_inq_varid(ncid_,'epsn_w',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid epsn_w'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
!      if(status/=0)stop 'Erreur nf_put_vara_real epsn_w'
!
!...................
! Ecrire le profil vertical de omega_w !20-05-18
!      varstart(4)=time_counter_  ; varcount(4)=1
!      varstart(3)=1              ; varcount(3)=kmax+1
!      varstart(1:2)=1            ; varcount(1:2)=1+ii_
!      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!       anyvar3d(i,j,k)=omega_w(i,j,k,1)*mask_t(i,j,k) &
!                          +filvalr4_*(1-mask_t(i,j,k))
!      enddo ; enddo ; enddo
!      status=nf_inq_varid(ncid_,'omega_w',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid omega_w'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:4)         &
!                                   ,varcount(1:4)         &
!              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
!      if(status/=0)stop 'Erreur nf_put_vara_real omega_w'

! Ecrire le vent WE, SN
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= gridrotcos_t(i,j)*uwind_t(i,j,1)       &
                       +gridrotsin_t(i,j)*vwind_t(i,j,1)
       anyvar3d(i,j,2)=-gridrotsin_t(i,j)*uwind_t(i,j,1)       &
                       +gridrotcos_t(i,j)*vwind_t(i,j,1)
      enddo ; enddo

      status=nf_inq_varid(ncid_,'Weastward',var_id)
      if(status/=0)stop ' erreur nf_inq_varid Weastward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real Weastward'

      status=nf_inq_varid(ncid_,'Wnorthward',var_id)
      if(status/=0)stop ' erreur nf_inq_varid Wnorthward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,2))
      if(status/=0)stop 'Erreur nf_put_vara_real Wnorthward'

!...................
      if(flag_abl==1) then !AAAAAAAAAAA>

      x2=(elapsedtime_now           -airseafile_prvtime(t2m_id))      &
        /(airseafile_nextime(t2m_id)-airseafile_prvtime(t2m_id))
      x2=min(max(x2,0.d00),1.d00) ; x1=1.-x2

       do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
        xy_t(i,j,2)=teta2_t(i,j,now)-273.15                    !05-04-16
        xy_t(i,j,0)=x1*teta0_t(i,j,1)+x2*teta0_t(i,j,2)-273.15 !05-04-16
       enddo ; enddo

! Ecrire la temperature potentielle de l'air a 2m:
      status=nf_inq_varid(ncid_,'teta2_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid teta2'
      status=nf_put_vara_double(ncid_,var_id      &
                                   ,varstart(1:3) &
                                   ,varcount(1:3) &
                   ,xy_t(i_:i_+ii_,j_:j_+ii_,2))    !05-04-16
!               ,teta2_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double teta2'

! Ecrire la temperature de surface dans le modele meteo
      status=nf_inq_varid(ncid_,'teta0_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid teta0'
      status=nf_put_vara_double(ncid_,var_id      &
                                   ,varstart(1:3) &
                                   ,varcount(1:3) &
                   ,xy_t(i_:i_+ii_,j_:j_+ii_,0))    !05-04-16
!               ,teta2_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double teta0'

      status=nf_inq_varid(ncid_,'teta2delta_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid teta2delta'
      status=nf_put_vara_double(ncid_,var_id      &
                                   ,varstart(1:3) &
                                   ,varcount(1:3) &
           ,teta2delta_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_real teta2delta'

! Ecrire l'humidite specifique a 2m:
      status=nf_inq_varid(ncid_,'q2_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid q2'
      status=nf_put_vara_double(ncid_,var_id          &
                                   ,varstart(1:3)     &
                                   ,varcount(1:3)     &
                   ,q2_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double q2'

      status=nf_inq_varid(ncid_,'q2delta_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid q2delta'
      status=nf_put_vara_double(ncid_,var_id          &
                                   ,varstart(1:3)     &
                                   ,varcount(1:3)     &
              ,q2delta_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double q2delta'

      if(flag_abl2==1) then !222222222>

      status=nf_inq_varid(ncid_,'uwind_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid uwind'
      status=nf_put_vara_double(ncid_,var_id          &
                                   ,varstart(1:3)     &
                                   ,varcount(1:3)     &
                ,uwind_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double uwind'

      status=nf_inq_varid(ncid_,'uwinddelta_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid uwinddelta'
      status=nf_put_vara_double(ncid_,var_id          &
                                   ,varstart(1:3)     &
                                   ,varcount(1:3)     &
           ,uwinddelta_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double uwinddelta'

      status=nf_inq_varid(ncid_,'vwind_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid vwind'
      status=nf_put_vara_double(ncid_,var_id          &
                                   ,varstart(1:3)     &
                                   ,varcount(1:3)     &
                ,vwind_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double uwind'

      status=nf_inq_varid(ncid_,'vwinddelta_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid vwinddelta'
      status=nf_put_vara_double(ncid_,var_id          &
                                   ,varstart(1:3)     &
                                   ,varcount(1:3)     &
           ,vwinddelta_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_double uwinddelta'

      endif                 !222222222>

      endif                !AAAAAAAAAAA>

! wstress_u wstress_v
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= gridrotcos_t(i,j)*wstress_u(i,j,1)       &
                       +gridrotsin_t(i,j)*wstress_v(i,j,1)
       anyvar3d(i,j,2)=-gridrotsin_t(i,j)*wstress_u(i,j,1)       &
                       +gridrotcos_t(i,j)*wstress_v(i,j,1)
      enddo ; enddo

      status=nf_inq_varid(ncid_,'wstress_u',var_id)
      if(status/=0)stop ' erreur nf_inq_varid wstress_u'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real wstress_u'

      status=nf_inq_varid(ncid_,'wstress_v',var_id)
      if(status/=0)stop ' erreur nf_inq_varid wstress_v'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,2))
      if(status/=0)stop 'Erreur nf_put_vara_real wstress_v'

! Along axis barotropic fluxes !14-01-19
!      const2=1.0/real(max(iteration2d_max_now+iteration2d_max_bef,1)) !11-03-19
!      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!       anyvar3d(i,j,1)=const2*(fluxbar_sumt_u(i,j,0)+fluxbar_sumt_u(i,j,1))
!       anyvar3d(i,j,2)=const2*(fluxbar_sumt_v(i,j,0)+fluxbar_sumt_v(i,j,1))
!      enddo ; enddo
!
!      status=nf_inq_varid(ncid_,'Oi2Dflux',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid Oi2Dflux'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:3)         &
!                                   ,varcount(1:3)         &
!                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
!      if(status/=0)stop 'Erreur nf_put_vara_real Oi2Dflux'
!
!      status=nf_inq_varid(ncid_,'Oj2Dflux',var_id)
!      if(status/=0)stop ' erreur nf_inq_varid Oj2Dflux'
!      status=nf_put_vara_real(ncid_,var_id                &
!                                   ,varstart(1:3)         &
!                                   ,varcount(1:3)         &
!                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,2))
!      if(status/=0)stop 'Erreur nf_put_vara_real Oj2Dflux'

! ssr:
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= ssr_w(i,j,1) 
      enddo ; enddo
      status=nf_inq_varid(ncid_,'ssr',var_id)
      if(status/=0)stop ' erreur nf_inq_varid ssr'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real ssr'

! snsf:
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= snsf_w(i,j,1) 
      enddo ; enddo
      status=nf_inq_varid(ncid_,'snsf',var_id)
      if(status/=0)stop ' erreur nf_inq_varid snsf'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real snsf'


! slhf:
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= slhf_w(i,j,1) 
      enddo ; enddo
      status=nf_inq_varid(ncid_,'slhf',var_id)
      if(status/=0)stop ' erreur nf_inq_varid slhf'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real slhf'

! sshf:
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= sshf_w(i,j,1) 
      enddo ; enddo
      status=nf_inq_varid(ncid_,'sshf',var_id)
      if(status/=0)stop ' erreur nf_inq_varid sshf'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real sshf'

! precipitation:
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= precipi_w(i,j,1) 
      enddo ; enddo
      status=nf_inq_varid(ncid_,'precipitation',var_id)
      if(status/=0)stop ' erreur nf_inq_varid precipitation'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real precipitation'

! ssh: !23-10-18
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)=ssh_int_w(i,j,1) !10-03-19
      enddo ; enddo
      status=nf_inq_varid(ncid_,'ssh',var_id)
      if(status/=0)stop ' erreur nf_inq_varid ssh'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real ssh'
!...................

      status=nf_close(ncid_)

      end subroutine my_outputs_point_vs_time
!#endif
!.............................................................
#ifdef bidon
      subroutine my_outputs_profiles
      implicit none
      double precision lon_,lat_
#ifdef synopsis
       subroutinetitle='my_outputs_profiles'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Centroid des stations S:
      lon_=51.897823 ; lat_=-46.843384

      call latlontoij(lon_*deg2rad,lat_*deg2rad,'loc')
      if(deci>1.and.deci<imax.and.decj>1.and.decj<jmax)call my_outputs_write('S')

! Centroid des stations N:
      lon_=51.788300 ; lat_=-46.167126

      call latlontoij(lon_*deg2rad,lat_*deg2rad,'loc')
      if(deci>1.and.deci<imax.and.decj>1.and.decj<jmax)call my_outputs_write('N')

! stationN         290         268
! stationS          23         208

      end subroutine my_outputs_profiles
#endif
!.............................................................
#ifdef bidon
      subroutine my_outputs_cumulvel(case_)
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='my_outputs_cumulvel'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(cumulvel_t)) then !- reset ->
       allocate(cumulvel_t(1:imax,1:jmax,kmax))
                cumulvel_t=0.
       return ! appelé depuis initial_main on ne fait rien d'autre!
      endif                               !- reset ->

      if(case_==1) then !------------->

! Chemin positivé cumulé
       cumulveltime=cumulveltime+dti_fw
       do k=1,kmax
       do j=1,jmax
       do i=1,imax
         cumulvel_t(i,j,k)=                                          &
         cumulvel_t(i,j,k)                                           &
        +dti_fw*sqrt( (0.5*(vel_u(i,j,k,1)+vel_u(i+1,j,k,1)))**2     &
                     +(0.5*(vel_v(i,j,k,1)+vel_v(i,j+1,k,1)))**2   )

       enddo
       enddo
       enddo

      endif             !------------->


      if(case_==0) then !\\\\\\\\\\\\\>

       cumulveltime=0.
       cumulvel_t=0.

      endif             !\\\\\\\\\\\\\>

      end subroutine my_outputs_cumulvel
#endif
!.............................................................
#ifdef bidon
      subroutine my_outputs_write(txt_)
      implicit none
      character*1 txt_
#ifdef synopsis
       subroutinetitle='my_outputs_write'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       i=nint(deci) ; j=nint(decj)
       open(unit=3,file='tmp/station_'//txt_,position='append')
        write(3,'(e13.6,100(1x,e12.5))')elapsedtime_now/86400.   &
                                       ,vel_u(i,j,1:kmax,1)      &
                                       ,vel_v(i,j,1:kmax,1)
       close(3)

       if(iteration3d==0) then  !----->
       open(unit=3,file='tmp/station_depth_'//txt_,position='append')
        write(3,'(e13.6,100(1x,e12.5))')elapsedtime_now/86400.   &
                                       ,depth_u(i,j,1:kmax)      &
                                       ,depth_v(i,j,1:kmax)
       close(3)

       endif                    !----->

      end subroutine my_outputs_write
#endif
!.............................................................
#ifdef bidon
      subroutine my_outputs_cumultime
      implicit none
#ifdef synopsis
       subroutinetitle='my_outputs_cumultime'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(cumultime_t)) then !------->
               allocate(cumultime_t(1:imax,1:jmax,kmax,3))
                        cumultime_t=0.
       return ! appelé depuis initial_main on ne fait rien d'autre!
      endif                                !------->

       do k=1,kmax
       do j=1,jmax
       do i=1,imax

         x1=sqrt( (0.5*(vel_u(i,j,k,1)+vel_u(i+1,j,k,1)))**2     &
                 +(0.5*(vel_v(i,j,k,1)+vel_v(i,j+1,k,1)))**2   )

         if(x1>0.2)cumultime_t(i,j,k,1)=cumultime_t(i,j,k,1)+dti_fw
         if(x1>0.4)cumultime_t(i,j,k,2)=cumultime_t(i,j,k,2)+dti_fw
         if(x1>0.6)cumultime_t(i,j,k,3)=cumultime_t(i,j,k,3)+dti_fw

       enddo
       enddo
       enddo

      end subroutine my_outputs_cumultime
#endif
!.............................................................
#ifdef bidon
      subroutine my_outputs_wave_beam_angle
      implicit none
      real*4 nz2_,pulsation2_
#ifdef synopsis
       subroutinetitle='my_outputs_wave_beam_angle'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Wave beam angle (pairaud et al, 2010)
! angle=SQRT( (w**2-f**2)/(N(z)**2-w**2)
! w: pulsation = 2*pi/period
! f: coriolis parameter
! N(z): Brunt-Vaisala freq: N(z)=sqrt(-g/rho*drho/dz)

      allocate(angle_wave_beam_w(imax,jmax))

      pulsation2_=(2.*pi/43200.)**2 ! w**2 (pulsation au carré à la frequence demi-diurne)
!     do k=2,kmax
      k=2
      do j=1,jmax
      do i=1,imax

       nz2_=(grav/rho)*abs(rhp_t(i,j,k)-  rhp_t(i,j,k-1))       &
                       /(depth_t(i,j,k)-depth_t(i,j,k-1))

       angle_wave_beam_w(i,j)=                                  &
                    sqrt( abs(pulsation2_-coriolis_t(i,j)**2)   &
                         /abs(nz2_-pulsation2_) )

      enddo
      enddo
!     enddo


      end subroutine my_outputs_wave_beam_angle
#endif
!.............................................................

!#ifdef bidon
      subroutine my_outputs_globtke_vs_time
      implicit none
! Energie cinetique globale en fonction du temps
#ifdef synopsis
       subroutinetitle='my_outputs_globtke_vs_time'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      sum1=0.
      sum2=0.
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
       x1=mask_u(i,j,k)*dz_u(i,j,k,1)*dxdy_u(i,j)*mask_i_u(i)*mask_j_u(j)
       sum1=sum1+x1
       sum2=sum2+x1*vel_u(i,j,k,1)**2
      enddo ; enddo ; enddo
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
       x1=mask_v(i,j,k)*dz_v(i,j,k,1)*dxdy_v(i,j)*mask_i_v(i)*mask_j_v(j)
       sum1=sum1+x1
       sum2=sum2+x1*vel_v(i,j,k,1)**2
      enddo ; enddo ; enddo

#ifdef parallele
      call mpi_allreduce(sum1                            &
                        ,sum1glb,1,mpi_double_precision, & !05-10-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2                            & !26-04-11
                        ,sum2glb,1,mpi_double_precision, &
                         mpi_sum,par%comm2d,ierr)
#else
      sum2glb=sum2
      sum1glb=sum1
#endif

      if(par%rank==0) then !------------->
       open(unit=3,file='tmp/globtke_vs_time',position='append')
        write(3,'(2(1x,e14.7))')elapsedtime_now/86400.,sum2glb/sum1glb
       close(3)
      endif                !------------->

      end subroutine my_outputs_globtke_vs_time
!#endif

!.............................................................
!#ifdef bidon
      subroutine my_outputs_point_vs_time_2d(name_,locationconv_,posx_,posy_) !09-09-16
      implicit none
      character(len=*) name_,locationconv_
      double precision :: lon_,lat_,posx_,posy_  
!                        ,filval_=-9999.d0
!     real :: filvalr4_=-9999.
      integer ii_
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='my_outputs_point_vs_time'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Eventuellement si on ne veut pas archiver à chaque iteration pour 
! limiter la taille du fichier, sortir si modulo(iteration3d,N)/=0
   
      if(modulo(iteration3d,50)/=0)return

! calculer les indices du point selon convention:

      if(trim(locationconv_)=='lonlat') then  !-----> ! convention longitude latitude
       lon_=posx_ ; lat_=posy_ ; ii_=1
       call latlontoij(lon_*deg2rad,lat_*deg2rad,'loc') ! returns deci, decj
       ii_=0 ; i_=nint(deci) ; j_=nint(decj) ! strategie un point (le + proche)
!              i_= int(deci) ; j_= int(decj) ! strategie quatre points encadrants
      else                                    !----->
       if(trim(locationconv_)=='ijglob') then !/////> ! convention longitude latitude
        ii_=0 ; i_=nint(posx_)-par%timax(1) ; j_=nint(posy_)-par%tjmax(1) ! strategie un point (le plus proche)
!               i_= int(posx_)-par%timax(1) ; j_= int(posy_)-par%tjmax(1) ! strategie quatre points encadrants
       else                                   !/////>
        stop 'STOP module_my_output err on locationconv_'
       endif                                  !/////>
      endif                                   !----->


! Si les indices (locaux) sont hors proc sortir:
!     write(10+par%rank,*)deci,decj
!     i1=211-par%timax(1)
!     j1=101-par%tjmax(1)
!     if(i1>2.and.i1<imax-1.and.j1>2.and.j1<jmax-1) then
!       write(10+par%rank,*)'i1,j1=',i1,j1
!       write(10+par%rank,*)'lon lat',lon_t(i1,j1)*rad2deg,lat_t(i1,j1)*rad2deg
!       call latlontoij(lon_t(i1,j1),lat_t(i1,j1),'loc') 
!       write(10+par%rank,*)'deci decj',deci,decj
!       call latlontoij(lon_t(10,10),lat_t(10,10),'loc') 
!       write(10+par%rank,*)'test 10 10=',deci,decj
!     endif

      if(i_<2.or.i_>imax-1.or.j_<2.or.j_>jmax-1)return

! Verification: afficher lon lat au point i_,j_ et comparer a lon,lat bouee lion:
!     write(6,*)'lon lat S26',lon_t(i_,j_)*rad2deg,lat_t(i_,j_)*rad2deg
!     write(6,*)'lon lat BL ',4.64,42.06
!     stop ' test verification'

      write(texte30,'(a,a,a,i0,a)') &
      'tmp/',trim(name_),'_',par%rank,'_2d.nc'

! Tout d'abord tester l'existence du fichier
      status=nf_open(trim(texte30),nf_write,ncid_)
      if(status==0) then !xxx>
       flag_header_=1 ! flag_header_=1 signifie fichier deja existant
       status=nf_close(ncid_) 
      else               !xxx>
       flag_header_=0 ! flag_header_=0 signifie fichier inexistant
      endif              !xxx>
      
      if(flag_header_==0) then !000000000000000>! flag_header_=0 signifie fichier inexistant

!     write(6,*)trim(name_),h_w(i_,j_)
!     do j1=int(decj),int(decj)+1
!     do i1=int(deci),int(deci)+1
!     write(6,*)trim(name_),i1+par%timax(1),j1+par%tjmax(1),h_w(i1,j1)
!     enddo
!     enddo
! Si le fichier n'existe pas commencer par faire l'entete du fichier
      status=nf_create(trim(texte30),nf_clobber,ncid_)
!     status=nf_create('toto.nc',nf_clobber,ncid_)
      if(status/=0) then
       write(6,'(a,a)')'Err nf_create file:',trim(texte30)
       stop ' STOP in module_my_output'
      endif

!     status=nf_def_dim(ncid_,'ni_t',2,i_t_dim)
      status=nf_def_dim(ncid_,'ni_t',1+ii_,i_t_dim)
      if(status/=0)stop 'erreur nf_def_dim ni_t'

!     status=nf_def_dim(ncid_,'nj_t',2,j_t_dim)
      status=nf_def_dim(ncid_,'nj_t',1+ii_,j_t_dim)
      if(status/=0)stop 'erreur nf_def_dim nj_t'

      status=nf_def_dim(ncid_,'nk_t',kmax,k_t_dim)
      if(status/=0)stop 'erreur nf_def_dim nk_t'

      status=nf_def_dim(ncid_,'nk_w',kmax+1,k_w_dim)
      if(status/=0)stop 'erreur nf_def_dim nk_t'

      status=nf_def_dim(ncid_,'time',nf_unlimited,time_dim)
      if(status/=0)stop ' erreur nf_def_dim time_dim'

      vardim(4)=time_dim
      vardim(3)= k_t_dim 
      vardim(2)= j_t_dim 
      vardim(1)= i_t_dim 


! Dans le cas de l'etude avec Andrea on sort les 4 points encadrants le
! point de station ce qui explique que l'on ajoute 2 dimensions
! horizontales supplementaires de dimensions reduites a 2.
! Certaines variables sont a 4 dimensions

      call elapsedtime2date(0.d0,i1,i2,i3,i4,i5,i6) ! enters 0. (i.e. elapsedtime_now=0)
                                                    ! returns i1,i2,...=y,m,d,h,m,s
! Time en secondes (egalement variable dimension)
      write(texte80(2),'(a14,i4,5(a1,i2))')                    & !units
      'seconds since ',i1,'-',i2,'-',i3,' ',i4,':',i5,':',i6
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0' 
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0' 
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'
      if(texte80(2)(32:32)==' ')texte80(2)(32:32)='0'


      status=nf_def_var(ncid_,'time',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time'
      status=nf_put_att_text(ncid_,var_id,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)stop ' erreur nf_put_att_text'
      status=nf_put_att_text(ncid_,var_id,'calendar',9,'gregorian')
      if(status/=0)stop ' erreur nf_put_att_text calendar'

! Time en jours
      write(texte80(2),'(a11,i4,5(a1,i2))')                    & !units
      'days since ',i1,'-',i2,'-',i3,' ',i4,':',i5,':',i6
      if(texte80(2)(17:17)==' ')texte80(2)(17:17)='0' 
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0' 
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0' 
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'

      status=nf_def_var(ncid_,'time_days',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_days'
      status=nf_put_att_text(ncid_,var_id,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)stop ' erreur nf_put_att_text'
      status=nf_put_att_text(ncid_,var_id,'calendar',9,'gregorian')
      if(status/=0)stop ' erreur nf_put_att_text calendar'

! Time en annee decimale:
      status=nf_def_var(ncid_,'time_years',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_years'
      status=nf_put_att_text(ncid_,var_id,'long_name',13 &
      ,'Decimal_years')

! Time date yyyymmddhhmmss:
      status=nf_def_var(ncid_,'time_calendar',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_calendar'
      status=nf_put_att_text(ncid_,var_id,'long_name',19 &
      ,'time_yyyymmddhhmmss')

      status=nf_def_var(ncid_,'longitude_t',nf_double,2,vardim(1:2),var_id)
      if(status/=0)stop ' erreur nf_def_var longitude_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',9 &
      ,'longitude')
      status=nf_put_att_text(ncid_,var_id,'units',6,'degree')

      status=nf_def_var(ncid_,'latitude_t',nf_double,2,vardim(1:2),var_id)
      if(status/=0)stop ' erreur nf_def_var latitude_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',8 &
      ,'latitude')
      status=nf_put_att_text(ncid_,var_id,'units',6,'degree')

      status=nf_def_var(ncid_,'h_w',nf_real,2,vardim(1:2),var_id)
      if(status/=0)stop ' erreur nf_def_var h_w'
      status=nf_put_att_text(ncid_,var_id,'long_name',10 &
      ,'Bathymetry')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'ssh',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var ssh'
      status=nf_put_att_text(ncid_,var_id,'long_name',3 &
      ,'ssh')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'IB',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var IB'
      status=nf_put_att_text(ncid_,var_id,'long_name',17 &
      ,'inverse_barometer')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

      vardim(3)=time_dim
      status=nf_def_var(ncid_,'sshobc',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var sshobc'
      status=nf_put_att_text(ncid_,var_id,'long_name',6,'sshobc')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')

      if(kmaxtide>0) then !pmxpmx>
      vardim(3)=time_dim
      status=nf_def_var(ncid_,'sshtotide',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var sshtotide'
      status=nf_put_att_text(ncid_,var_id,'long_name',9,'sshtotide')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')
      endif               !pmxpmx>

#ifdef stokes
      vardim(3)=time_dim
      status=nf_def_var(ncid_,'sshstokes',nf_real,3,vardim(1:3),var_id)
      if(status/=0)stop ' erreur nf_def_var sshobc'
      status=nf_put_att_text(ncid_,var_id,'long_name',9,'sshstokes')
      status=nf_put_att_text(ncid_,var_id,'units',1,'m')
#endif

! Global attributs:
      k10=len(trim(name_))
      status=nf_put_att_text(ncid_,nf_global,'Point',k10,trim(name_))
      k10=len(trim(model_name))
      status=nf_put_att_text(ncid_,nf_global,'Model',k10,trim(model_name))
      write(texte30,'(a,i0,a,i0)')'global indexes: ',i_+par%timax(1),' ',j_+par%tjmax(1)
      k10=len(trim(texte30))
      status=nf_put_att_text(ncid_,nf_global,'Position',k10,trim(texte30))

!     status=nf_put_att_double(ncid_,nf_global,'bathymetry',nf_double,1,h_w(i_,j_))
!     x1=lon_t(i_,j_)*rad2deg ; x2=lat_t(i_,j_)*rad2deg 
!     status=nf_put_att_double(ncid_,nf_global,'longitude',nf_double,1,x1)
!     status=nf_put_att_double(ncid_,nf_global,'latitude' ,nf_double,1,x2)

! 194 continue
      status=nf_enddef(ncid_)

      if(status/=0)stop 'Err module_my_outputs nf_enddef'

! Ecrire longitude_t
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyv3d(i,j,1,1)=lon_t(i,j)*rad2deg
      enddo        ; enddo
      status=nf_inq_varid(ncid_,'longitude_t',var_id)
      if(status/=0)stop 'erreur nf_inq_varid longitude_t'
      status=nf_put_vara_double(ncid_,var_id        &
                                     ,varstart(1:2) &
                                     ,varcount(1:2) &
                       ,anyv3d(i_:i_+ii_,j_:j_+ii_,1,1))
      if(status/=0)stop 'erreur nf_put_var_double longitude_t'

! Ecrire latitude_t
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyv3d(i,j,1,1)=lat_t(i,j)*rad2deg
      enddo        ; enddo
      status=nf_inq_varid(ncid_,'latitude_t',var_id)
      if(status/=0)stop 'erreur nf_inq_varid latitude_t'
      status=nf_put_vara_double(ncid_,var_id        &
                                     ,varstart(1:2) &
                                     ,varcount(1:2) &
                       ,anyv3d(i_:i_+ii_,j_:j_+ii_,1,1))
      if(status/=0)stop 'erreur nf_put_var_double latitude_t'

! Ecrire h_w
      varstart(1:2)=1  ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar2d(i,j)=h_w(i,j)
      enddo        ; enddo
      status=nf_inq_varid(ncid_,'h_w',var_id)
      if(status/=0)stop 'erreur nf_inq_varid h_w'
      status=nf_put_vara_real(ncid_,var_id        &
                                   ,varstart(1:2) &
                                   ,varcount(1:2) &
                       ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'erreur nf_put_var_double h_w'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err module_my_outputs nf_close'

      endif                    !000000000000000>! flag_header_=0 signifie fichier inexistant

! Ecrire dans le fichier netcdf (en ajoutant au fichier existent au fur et a mesure que la
! simu avance)

      write(texte30,'(a,a,a,i0,a)') &
      'tmp/',trim(name_),'_',par%rank,'_2d.nc'

      status=nf_open(trim(texte30),nf_write,ncid_) ! noter l'option "nf_write"
!     status=nf_open('toto.nc',nf_write,ncid_) ! noter l'option "nf_write"
      if(status/=0) then
       write(6,'(a,a)')'Err nf_open file:',trim(texte30)
       stop ' STOP in module_my_output'
      endif

      status=nf_inq_dimid(ncid_,'time',dim_t_id)
      if(status==0)status=nf_inq_dimlen(ncid_,dim_t_id,time_counter_)

      time_counter_=time_counter_+1

!...................
! Ecrire le temps en secondes (variable dimension):
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time'
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,elapsedtime_now)
      if(status/=0)stop 'Erreur nf_put_vara_double time'

!...................
! Ecrire le temps en jours
      x1=elapsedtime_now/86400.
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time_days',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time_days'
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,x1)
      if(status/=0)stop 'Erreur nf_put_vara_double time_days'

!...................
! Ecrire le temps en annees decimales: 
      x1=year_now                                                     &
       +(elapsedtime_now               -seconds_since_1jan(year_now)) &
       /(seconds_since_1jan(year_now+1)-seconds_since_1jan(year_now))
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time_years',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time_years'
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,x1)
      if(status/=0)stop 'Erreur nf_put_vara_double time_years'

!...................
! Ecrire le temps sous forme de date yyyymmddhhmmss !30-07-15
      varstart(4)=time_counter_  ; varcount(4)=1
      status=nf_inq_varid(ncid_,'time_calendar',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid time_calendar'
      write(date_,'(I4.4,5(I2.2))')year_now,month_now,day_now,hour_now,minute_now,second_now !      15-12-20
      read(date_,'(f14.0)') x1                                                           !      15-12-20
        !  x1=second_now             &
        !    +minute_now*10         &
        !    +hour_now  *1000       &
        !    +day_now   *100000     &
        !    +month_now *10000000   &
        !    +year_now  *1000000000
      status=nf_put_vara_double(ncid_,var_id       &
                                     ,varstart(4)  &
                                     ,varcount(4)  &
                                     ,x1)
      if(status/=0)stop 'Erreur nf_put_vara_double time_calendar'

!...................
! Ecrire SSH
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!      anyvar2d(i,j)=ssh_int_w(i,j,1)
       anyvar2d(i,j)=ssh_w(i,j,1)
      enddo ; enddo
      status=nf_inq_varid(ncid_,'ssh',var_id)
      if(status/=0)stop ' erreur nf_inq_varid ssh'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                  ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'Erreur nf_put_vara_real ssh'

!...................
! Ecrire IB
      const1=-1./grav/rho     !17-05-11
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar2d(i,j)=const1*(pss_w(i,j,1)-pss_mean(1))
      enddo ; enddo
      status=nf_inq_varid(ncid_,'IB',var_id)
      if(status/=0)stop ' erreur nf_inq_varid IB'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                  ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'Erreur nf_put_vara_real IB'


!...................
! Ecrire sshobc
      x2=timeweightobc(ssh_id) ; x0=1.-x2
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar2d(i,j)=(x2*sshobc_w(i,j,2)+x0*sshobc_w(i,j,0))
      enddo ; enddo
      status=nf_inq_varid(ncid_,'sshobc',var_id)
      if(status/=0)stop ' erreur nf_inq_varid sshobc'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                  ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'Erreur nf_put_vara_real sshobc'

!...................
! Ecrire sshtotide_w
      if(kmaxtide>0) then !pmxpmx>
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       do ktide=1,kmaxtide
! Temps t+1:  (pour C.L. barotrope)
      time1=frqtide(ktide)*((elapsedtime_aft-ti0tide(ktide))+0.5*dte_lp) & !
            +v0tide(ktide)+utide(ktide,1)                   
      const3=cos(time1)*rampe                                         & !23/12/08
                       *ftide(ktide,1)                                  !21-12-09
      const4=sin(time1)*rampe                                         &
                       *ftide(ktide,1)                                  !21-12-09
        anyvar2d(i,j)=                                                &
        anyvar2d(i,j)*passetide(ktide,1)+                             &
       (sshtidecos_w(i,j,ktide)*const3+                               &
        sshtidesin_w(i,j,ktide)*const4)
       enddo ! ktide
      enddo ; enddo
      status=nf_inq_varid(ncid_,'sshtotide',var_id)
      if(status/=0)stop ' erreur nf_inq_varid sshtotide'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                  ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'Erreur nf_put_vara_real sshtotide'
      endif               !pmxpmx>

#ifdef stokes
!...................
! Ecrire sshstokes
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar2d(i,j)=sshstokes_w(i,j)
      enddo ; enddo
      status=nf_inq_varid(ncid_,'sshstokes',var_id)
      if(status/=0)stop ' erreur nf_inq_varid sshstokes'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                  ,anyvar2d(i_:i_+ii_,j_:j_+ii_))
      if(status/=0)stop 'Erreur nf_put_vara_real sshstokes'
#endif

      status=nf_close(ncid_)

      end subroutine my_outputs_point_vs_time_2d

!............................................................................
#ifdef bidon
      subroutine my_outputs_bio_sum !25-01-17
      implicit none

! Moyenne sur le volume (de tout le domaine) des variables bio_t
! Ecriture dans fichier tmp/biomeanXX (XX=numero "vb" de la variable)

      do vb=1,vbmax

       sum1=0.
       sum2=0.
       do k=1,kmax
       do j=1,jmax
       do i=1,imax
        x1=0.5*(dz_t(i,j,k,1)+dz_t(i,j,k,2))*dxdy_t(i,j) &
              *mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)      
        sum1=sum1+x1
        sum2=sum2+x1*bio_t(i,j,k,vb)
       enddo
       enddo
       enddo

       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       if(par%rank==0) then !>>>
        write(texte30,'(a,i0)')'tmp/biomean',vb
        open(unit=3,file=texte30,position='append')
         write(3,*)elapsedtime_now/86400.,sum2glb/sum1glb
        close(3)
       endif                !>>>

      enddo

      end subroutine my_outputs_bio_sum
#endif
!............................................................................
#ifdef bidon
      subroutine my_outputs_tem_sum !25-01-17
      implicit none

! Moyenne sur le volume (de tout le domaine) des variables bio_t
! Ecriture dans fichier tmp/biomeanXX (XX=numero "vb" de la variable)

       sum1=0.
       sum2=0.
       do k=1,kmax
       do j=1,jmax
       do i=1,imax
        x1=0.5*(dz_t(i,j,k,1)+dz_t(i,j,k,2))*dxdy_t(i,j) &
              *mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)      
        sum1=sum1+x1
        sum2=sum2+x1*tem_t(i,j,k,2)
       enddo
       enddo
       enddo

       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       if(par%rank==0) then !>>>
        write(texte30,'(a,i0)')'tmp/temmean',vb
        open(unit=3,file=texte30,position='append')
         write(3,*)elapsedtime_now/86400.,sum2glb/sum1glb
        close(3)
       endif                !>>>

      end subroutine my_outputs_tem_sum
#endif
!............................................................................
#ifdef bidon
      subroutine my_outputs_ssh_int_sum !23-10-19
      implicit none

! Moyenne sur le volume (de tout le domaine) des variables bio_t
! Ecriture dans fichier tmp/biomeanXX (XX=numero "vb" de la variable)

       sum1=0.
       sum2=0.
       do j=1,jmax
       do i=1,imax
        x1=dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)      
        sum1=sum1+x1
        sum2=sum2+x1*ssh_int_w(i,j,2)
       enddo
       enddo

       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       if(par%rank==0) then !>>>
        write(texte30,'(a,i0)')'tmp/ssh_intmean',vb
        open(unit=3,file=texte30,position='append')
         write(3,*)elapsedtime_now/86400.,sum2glb/sum1glb
        close(3)
       endif                !>>>

      end subroutine my_outputs_ssh_int_sum
#endif
!............................................................................
#ifdef bidon
      subroutine my_outputs_straitflux_i !02-07-17
      implicit none
! WARNING: DO NOT USE flux_lowfreq IF MORE THAN ONE SECTION !!!!

! Flux across an along Oj section
       i=447-par%timax(1) 
! Point 1 in global indexes:
       j1=319-par%tjmax(1) 
! Point 2 in global indexes:
       j2=330-par%tjmax(1)

! Give a name to the output file:
        texte30='tmp/flux_section_1'

       sum1=0.

       if(i>=1.and.i<=imax+1) then !m0v0m>

        do j=j1,j2
         if(j>=1.and.j<=jmax) then !pmxpmx>
          do k=1,kmax
           sum1=sum1+veldydz_u(i,j,k,1)*mask_u(i,j,k)*mask_i_u(i)*mask_j_u(j)
          enddo
         endif                     !pmxpmx>
        enddo

       endif                       !m0v0m>

       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)

! WARNING: DO NOT USE flux_lowfreq IF MORE THAN ONE SECTION !!!!
      x1=dti_fw/86400.
      flux_lowfreq=(1-x1)*flux_lowfreq+x1*sum1glb

       if(par%rank==0) then !>>>
        open(unit=3,file=texte30,position='append')
         write(3,*)elapsedtime_now/86400.,sum1glb,flux_lowfreq
        close(3)
       endif                !>>>

      end subroutine my_outputs_straitflux_i
#endif
!............................................................................

      subroutine my_outputs_mean_z_profile
      implicit none

      if(.not.allocated(zprofile_depth)) then !>>>
               allocate(zprofile_depth(kmax)) ; zprofile_depth=0.
!              allocate(zprofile_tem  (kmax)) ; zprofile_tem=0.
!              allocate(zprofile_sal  (kmax)) ; zprofile_sal=0.
      endif                                   !>>>

      if(par%rank==0) then
       open(unit=3,file='tmp/averaged_zprofile.txt',position='append')
       write(3,*)'Elapsed time in days',elapsedtime_now/86400.
      endif

!     x1=hmax/real(kmax)
!     zprofile_depth(kmax)=-0.5*x1
!     do k=kmax-1,1,-1
!      zprofile_depth(k)=zprofile_depth(k+1)-x1
!     enddo
       do k=1,kmax
       sum1=0.
       sum2=0.
       do j=1,jmax
       do i=1,imax
         sum1=sum1+h_w(i,j)*mask_t(i,j,kmax)
         sum2=sum2+h_w(i,j)*mask_t(i,j,kmax)*depth_t(i,j,k)
       enddo
       enddo
       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       zprofile_depth(k)=sum2glb/sum1glb
      enddo

!...............................................

      do k1=kmax,1,-1

       sum1=0.
       sum2=0.
       sum3=0.
       sum4=0.
       sum5=0.
       do k=1,kmax-1
       do j=1,jmax
       do i=1,imax

        if(depth_t(i,j,k  )<=zprofile_depth(k1).and. &
           depth_t(i,j,k+1)> zprofile_depth(k1)) then !ooo>

         rap=(zprofile_depth(k1)-depth_t(i,j,k  )) &
            /(  depth_t(i,j,k+1)-depth_t(i,j,k  ))

         sum1=sum1+mask_t(i,j,kmax)
         sum2=sum2+mask_t(i,j,kmax)*(rap*tem_t(i,j,k+1,1)+(1.-rap)*tem_t(i,j,k,1))

                       x1=timeweightobc(trc_id) *temobc_t(i,j,k  ,2) &
                     +(1.-timeweightobc(trc_id))*temobc_t(i,j,k  ,0)
                       x2=timeweightobc(trc_id) *temobc_t(i,j,k+1,2) &
                     +(1.-timeweightobc(trc_id))*temobc_t(i,j,k+1,0)
         sum3=sum3+mask_t(i,j,kmax)*(rap*x2+(1.-rap)*x1)


         sum4=sum4+mask_t(i,j,kmax)*(rap*sal_t(i,j,k+1,1)+(1.-rap)*sal_t(i,j,k,1))

                       x1=timeweightobc(trc_id) *salobc_t(i,j,k  ,2) &
                     +(1.-timeweightobc(trc_id))*salobc_t(i,j,k  ,0)
                       x2=timeweightobc(trc_id) *salobc_t(i,j,k+1,2) &
                     +(1.-timeweightobc(trc_id))*salobc_t(i,j,k+1,0)
         sum5=sum5+mask_t(i,j,kmax)*(rap*x2+(1.-rap)*x1)
        endif                                         !ooo>


       enddo
       enddo
       enddo

       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum4,sum4glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum5,sum5glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!      zprofile_tem(k1)=sum2glb/sum1glb
       if(par%rank==0)write(3,'(i3,5(1x,e14.7),a50)')k1              &
                                                 ,zprofile_depth(k1) &
                                                 ,sum2glb/sum1glb    &
                                                 ,sum3glb/sum1glb    &
                                                 ,sum4glb/sum1glb    &
                                                 ,sum5glb/sum1glb    &
                                                 ,' z tem temobc'

      enddo

      if(par%rank==0)close(3)

      end subroutine my_outputs_mean_z_profile

!............................................................................
#ifdef bidon
! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
      subroutine my_outputs_sal_sum !25-01-17
      implicit none
      double precision :: sumsal0_=0.
      integer(kind=1) :: loop_=1,loopmin_=1

! Moyenne sur le volume (de tout le domaine) des variables bio_t
! Ecriture dans fichier tmp/biomeanXX (XX=numero "vb" de la variable)

! si ifsbeg n'ont pas ete initialisee par la subroutine my_outputs_obcsaltflux mettre les valeurs par defaut:
       if(ifsbeg==-999) then !secu
        ifsbeg=1 ; ifsend=imax ; jfsbeg=1 ; jfsend=jmax !"free sponge (fs)" indice min max values
       endif                 !secu 

       loopmin_=1 
       if(iteration3d==0)loopmin_=0

       do loop_=loopmin_,1

        if(loop_==0) then !>>>
         k1=-1 ; k2=0
        else              !>>>
         k1=1  ; k2=2
        endif             !>>>

        sum1=0.  ; sum2=0.
        do k=1,kmax ; do j=jfsbeg,jfsend ; do i=ifsbeg,ifsend
         sum1=sum1+0.5*(dz_t(i,j,k,k1)+dz_t(i,j,k,k2))*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
         sum2=sum2+0.5*(dz_t(i,j,k,k1)+dz_t(i,j,k,k2))*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*sal_t(i,j,k,k2)
        enddo ; enddo ; enddo
        call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        if(loop_==0)sumsal0_=sum2glb

       enddo ! loop_

! Bilan des river boxes: !16-03-17
      sum3=0.
      do kr=1,nriver ! boucle sur KR début.
       if(rivertrc_inout(kr)==1)sum3=sum3+dxdy_t(iriver(kr,1),jriver(kr,1))*10.*river_s(kr)
      enddo

       call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)

       if(par%rank==0) then !>>>
        write(texte30,'(a)')'tmp/salmean' !19-03-17
        open(unit=3,file=texte30,position='append')
!        write(3,*)elapsedtime_now/86400.,sum2glb/sum1glb
!        write(3,*)real(elapsedtime_now/86400.),sum2glb,real(sum3glb),real(sum2glb/sum1glb)
!        write(3,*)real(elapsedtime_now/86400.),sum2glb,real(sum3glb),real(sum2glb/sum1glb)
         write(3,*)real(elapsedtime_now/86400.),sum2glb-sumsal0_
        close(3)
       endif                !>>>


      end subroutine my_outputs_sal_sum
#endif

!............................................................................
#ifdef bidon

! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
      subroutine my_outputs_obcsaltflux(txt_,id_flux_) !29-01-19
      implicit none
      integer id_flux_
      character(len=*)txt_

      flag=0
      if(txt_=='i')   flag=1
      if(txt_=='j')   flag=2
      if(txt_=='mpi') flag=3
      if(flag==0) stop 'Err 13 my_outputs_obcsaltflux'


! Initialisation
!     if(iteration3d==iteration3d_restart) then !initial>
       if(.not.allocated(obcsaltfluxi1_)) then !>>>>>>>>>

        allocate(obcsaltfluxi1_(nlayer_))        ; obcsaltfluxi1_=0.
        allocate(obcsaltfluxi2_     (nlayer_))   ; obcsaltfluxi2_=0.
        allocate(obcsaltfluxj1_     (nlayer_))   ; obcsaltfluxj1_=0.
        allocate(obcsaltfluxj2_     (nlayer_))   ; obcsaltfluxj2_=0. 
        allocate(obcsaltfluxi1_glb_ (nlayer_))   ; obcsaltfluxi1_glb_=0.
        allocate(obcsaltfluxi2_glb_ (nlayer_))   ; obcsaltfluxi2_glb_=0.
        allocate(obcsaltfluxj1_glb_ (nlayer_))   ; obcsaltfluxj1_glb_=0.
        allocate(obcsaltfluxj2_glb_ (nlayer_))   ; obcsaltfluxj2_glb_=0.
        allocate(layerdepth_        (nlayer_+1)) ; layerdepth_=0. ! from bottom to surface

        layerdepth_(1)=-hmax
        layerdepth_(2)=-50.
        layerdepth_(nlayer_+1)=999.  ! Grande valeur obligatoire


        sponge_l_int=ceiling(sponge_l)
!      sponge_l_int=10
        ifsbeg=min(max(1   +sponge_l_int-par%timax(1),1),imax+1)
        jfsbeg=min(max(1   +sponge_l_int-par%tjmax(1),1),jmax+1)
        ifsend=min(max(iglb-sponge_l_int-par%timax(1),0),imax)
        jfsend=min(max(jglb-sponge_l_int-par%tjmax(1),0),jmax)
! Le decallage de 1 des bornes pour que la boucle do i=imax+1,imax ne s'excute pas ....


!      write(10+par%rank,*)'ifsbeg,jfsbeg,ifsend,jfsend',ifsbeg,jfsbeg,ifsend,jfsend
!      write(10+par%rank,*)'imax,jmax',imax,jmax
!      write(10+par%rank,*)'sponge_l_int=',sponge_l_int
!      write(10+par%rank,*)'layerdepth_',layerdepth_
!      stop 'coco'
       endif                                  !>>>>>>>>>
!     endif                                     !initial>

      if(flag==1) then       !iii>

        i=1+sponge_l_int-par%timax(1)     ! i=1 si sans eponge et 1 proc
        if(i>=1.and.i<=imax+1) then !11111>
         do k1=1,nlayer_
         do k=1,kmax ; do j=jfsbeg,jfsend
          if(layerdepth_(k1)<=depth_u(i,j,k).and.depth_u(i,j,k)<layerdepth_(k1+1)) then !ooo>
           obcsaltfluxi1_(k1)=obcsaltfluxi1_(k1)-0.5*dti_lpsub*anyv3d(i,j,k,id_flux_)*mask_i_u(i)*mask_j_u(j)
          endif                                                                         !ooo>
         enddo       ; enddo 
         enddo
        endif                       !11111>

        i=iglb+1-sponge_l_int-par%timax(1) ! i=imax+1 si sans eponge et 1 proc
        if(i>=1.and.i<=imax+1) then !22222>
         do k1=1,nlayer_
         do k=1,kmax ; do j=jfsbeg,jfsend
          if(layerdepth_(k1)<=depth_u(i,j,k).and.depth_u(i,j,k)<layerdepth_(k1+1)) then !ooo>
           obcsaltfluxi2_(k1)=obcsaltfluxi2_(k1)+0.5*dti_lpsub*anyv3d(i,j,k,id_flux_)*mask_i_u(i)*mask_j_u(j)
          endif                                                                         !ooo>
         enddo       ; enddo 
         enddo
        endif                       !22222>

      endif                  !iii>

      if(flag==2) then       !jjj>

        j=1+sponge_l_int-par%tjmax(1) ! j=1 si sans eponge et 1 proc
        if(j>=1.and.j<=jmax+1) then !33333>
         do k1=1,nlayer_
         do k=1,kmax ; do i=ifsbeg,ifsend
          if(layerdepth_(k1)<=depth_v(i,j,k).and.depth_v(i,j,k)<layerdepth_(k1+1)) then !ooo>
           obcsaltfluxj1_(k1)=obcsaltfluxj1_(k1)-0.5*dti_lpsub*anyv3d(i,j,k,id_flux_)*mask_i_v(i)*mask_j_v(j)
          endif                                                                         !ooo>
         enddo       ; enddo 
         enddo
        endif                       !33333>

        j=jglb+1-sponge_l_int-par%tjmax(1) ! j=jmax+1 si sans eponge et 1 proc²
        if(j>=1.and.j<=jmax+1) then !44444>
         do k1=1,nlayer_
         do k=1,kmax ; do i=ifsbeg,ifsend
          if(layerdepth_(k1)<=depth_v(i,j,k).and.depth_v(i,j,k)<layerdepth_(k1+1)) then !ooo>
           obcsaltfluxj2_(k1)=obcsaltfluxj2_(k1)+0.5*dti_lpsub*anyv3d(i,j,k,id_flux_)*mask_i_v(i)*mask_j_v(j)
          endif                                                                         !ooo>
         enddo       ; enddo 
         enddo
        endif                       !44444>

      endif                  !jjj>

      if(flag==3) then !ooo>

         call mpi_allreduce(obcsaltfluxi1_,obcsaltfluxi1_glb_,nlayer_,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         call mpi_allreduce(obcsaltfluxi2_,obcsaltfluxi2_glb_,nlayer_,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         call mpi_allreduce(obcsaltfluxj1_,obcsaltfluxj1_glb_,nlayer_,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         call mpi_allreduce(obcsaltfluxj2_,obcsaltfluxj2_glb_,nlayer_,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(par%rank==0)  then !--->

! Toutes les couches verticales cumulees:
           write(texte30,'(a)')'tmp/obcflux'
           open(unit=3,file=texte30,position='append')
                      write(3,*)                    &
                      real(elapsedtime_now/86400.)  &
                     ,real(sum(obcsaltfluxi1_glb_)) & ! i=ifsbeg boundary
                     ,real(sum(obcsaltfluxi2_glb_)) & ! i=ifsend boundary
                     ,real(sum(obcsaltfluxj1_glb_)) & ! j=jfsend boundary
                     ,real(sum(obcsaltfluxj2_glb_))   ! j=jfsend boundary
           close(3)

! Couche par couche:
           do k1=1,nlayer_
           write(texte30,'(a,i0)')'tmp/obcflux',k1
           open(unit=3,file=texte30,position='append')
                      write(3,*)                   &
                      real(elapsedtime_now/86400.) &
                     ,real(obcsaltfluxi1_glb_(k1)) & ! i=ifsbeg boundary
                     ,real(obcsaltfluxi2_glb_(k1)) & ! i=ifsend boundary
                     ,real(obcsaltfluxj1_glb_(k1)) & ! j=jfsend boundary
                     ,real(obcsaltfluxj2_glb_(k1))   ! j=jfsend boundary
           close(3)

           enddo


         endif                 !--->
      endif            !ooo>

      end subroutine my_outputs_obcsaltflux
#endif

!............................................................................
!#ifdef bidon
      subroutine my_outputs_glider !05-06-19
      implicit none
      double precision :: lon_,lat_,z_
      REAL*8 lon_glider,lat_glider,z_glider,sum0,sum1,sum2,loop1_
      integer :: year_glider,month_glider,day_glider,hour_glider,minute_glider,second_glider,flag_glider=0, &
                 kv
      integer :: elapsedtime_glider_now_
 3135 continue

!      write(10+par%rank,*)'.................................'
!      write(10+par%rank,*)'elapsedtime_glider_now_',elapsedtime_glider_now_
!      write(10+par%rank,*)'elapsedtime_now        ',elapsedtime_now
!      write(10+par%rank,*)'DATE MODELE',month_now,day_now,hour_now,minute_now,second_now

!     if(flag_glider==0) open(unit=543,file='glidersort.txt')
      if(flag_glider==0) open(unit=543,file='gliderLCsort.txt')
      flag_glider=1

!      do while (elapsedtime_glider_now_<elapsedtime_now) !PPP>

!      read(543,*)lon_,lat_,z_,year_,month_,day_,hour_,minute_,second_
       read(543,*,end=3136)lon_glider,lat_glider,z_glider,year_glider,month_glider,day_glider,hour_glider,minute_glider,second_glider
!      write(10+par%rank,*)'DATE GLIDER',year_glider,month_glider,day_glider,hour_glider,minute_glider,second_glider
!      write(10+par%rank,*)'DANS GLIDER lon_glider',lon_glider
!      write(10+par%rank,*)'DANS GLIDER lat_glider',lat_glider
!      write(10+par%rank,*)'DANS GLIDER z_glider  ',z_glider

       call datetokount(year_glider,month_glider,day_glider,hour_glider,minute_glider,second_glider)  ! calcule elapsedtime_out:
                                                                ! temps en secondes ecoulees
                                                                ! depuis le debut du run
!      elapsedtime_glider_bef_=elapsedtime_glider_now_
       elapsedtime_glider_now_=elapsedtime_out
!      write(10+par%rank,*)'DOWHILE elapsedtime_glider_now_',elapsedtime_glider_now_

!      enddo                                              !PPP>

!    if(elapsedtime_glider_now_<elapsedtime_now)goto 2658




!     if( elapsedtime_now> elapsedtime_glider_bef_.and.  
!         elapsedtime_now<=elapsedtime_glider_now_) then !ooooo>

      lon_=lon_glider*deg2rad
      lat_=lat_glider*deg2rad
      z_=z_glider
!     write(10+par%rank,*)'------DANS GLIDER lon_ B',lon_*rad2deg
!     write(10+par%rank,*)'------DANS GLIDER lat_ B',lat_*rad2deg
!     write(10+par%rank,*)'------DANS GLIDER z_   B',z_
     call latlonztoijk(lon_,lat_,z_,'loc') ! calcule deci,decj,deck les indices decimaux

      i=nint(deci) ; j=nint(decj) ; k=nint(deck)

!     write(10+par%rank,*)'------DANS GLIDER i,j,k',i,j,k

! Ici on fait un test pour savoir si le glider est dans le sous_domaine
! ou dehors:
       if(i>1.and.i<imax.and.j>1.and.j<jmax) then ! test reussi

!     write(10+par%rank,*)'------DANS GLIDER je suis dans le dom'

       if(elapsedtime_now       < elapsedtime_glider_now_.and.  &
          elapsedtime_now+dti_fw>=elapsedtime_glider_now_) then !ooo>
         x0=deck-nint(deck)
         x1=(1.-x0)


! vitesse verticale (copie du graph_out.F90)
!        do kv=1,kmax+1  !---kv-->
!         if(kv==1)then
!          sum0=0.
!          sum1=0.
!         else
!          sum0=sum0+dz_t(i,j,kv-1,0)
!          sum1=sum1+dz_t(i,j,kv-1,1)
!         endif
!
!         sum2=0.
!         do loop1_=max0(kv-1,1),min0(kv,kmax)  !--loop1_-->
!         do i1=i,i+1
!           sum2=sum2+                                                   &
!           vel_u(i1,j,loop1_,1)*(depth_t(i1,j,loop1_)-depth_t(i1-1,j,loop1_)) &
!           /dx_u(i1,j)
!         enddo
!         do j1=j,j+1
!           sum2=sum2+                                                   &
!           vel_v(i,j1,loop1_,1)*(depth_t(i,j1,loop1_)-depth_t(i,j1-1,loop1_)) &
!           /dy_v(i,j1)
!         enddo
!         enddo                           !--loop1_--<


! w=omega+dz/dx*u+dz/dy*v+dz/dt:
!         anyvar3d(i,j,kv)=                                               &
!           mask_t(i,j,kv)*(                                              &
!          omega_w(i,j,kv,1)+0.25*sum2+(sum1-sum0)/dti_fw)

!        enddo      !---kv--<


!        open(unit=4,file='tmp/glider2_out.txt',position='append')
         open(unit=4,file='tmp/gliderLC_out.txt',position='append')
         if(k.lt.kmax) then

!         write(4,'(3(1x,e14.7),12(1x,i4),15(1x,e14.7),3(1x,i4))') &
          write(4,'(3(1x,e14.7),12(1x,i4),7(1x,e14.7))') &
          lon_*rad2deg,lat_*rad2deg,z_,year_glider,month_glider,day_glider,hour_glider,minute_glider,second_glider &
                                      ,year_now   ,month_now   ,day_now ,hour_now   ,minute_now   ,second_now    &
!        ,tem_t(i,j,k,1),sal_t(i,j,k,1) &
         ,tem_t(i,j,k+1,1)*x0+tem_t(i,j,k,1)*x1,sal_t(i,j,k+1,1)*x0+sal_t(i,j,k,1)*x1 &
!        ,bio_t(i,j,k,isynechl) &
         ,bio_t(i,j,k+1,isynechl)*x0+bio_t(i,j,k,isynechl)*x1 &
!         ,bio_t(i,j,k,inanochl) &
         ,bio_t(i,j,k+1,inanochl)*x0+bio_t(i,j,k,inanochl)*x1 &
!         ,bio_t(i,j,k,idiachl) &  
         ,bio_t(i,j,k+1,idiachl)*x0+bio_t(i,j,k,idiachl)*x1  &
!         ,bio_t(i,j,k,initrate) &
         ,bio_t(i,j,k+1,initrate)*x0+bio_t(i,j,k,initrate)*x1 &
!         ,bio_t(i,j,k,ioxygen) &
         ,bio_t(i,j,k+1,ioxygen)*x0+bio_t(i,j,k,ioxygen)*x1  !&
!!         ,0.5*(                        &
!!            ( vel_u(i  ,j  ,k,1)                       &
!!             +vel_u(i+1,j  ,k,1))*gridrotcos_t(i,j)    &
!!           +( vel_v(i  ,j  ,k,1)                       &
!!             +vel_v(i  ,j+1,k,1))*gridrotsin_t(i,j))*mask_t(i,j,k) &
!         ,(0.5*(                        &
!            ( vel_u(i  ,j  ,k+1,1)                       &
!             +vel_u(i+1,j  ,k+1,1))*gridrotcos_t(i,j)    &
!           +( vel_v(i  ,j  ,k+1,1)                       &
!             +vel_v(i  ,j+1,k+1,1))*gridrotsin_t(i,j))*mask_t(i,j,k+1))*x0 + &
!          (0.5*( &
!            ( vel_u(i  ,j  ,k,1)                       &
!             +vel_u(i+1,j  ,k,1))*gridrotcos_t(i,j)    &
!           +( vel_v(i  ,j  ,k,1)                       &
!             +vel_v(i  ,j+1,k,1))*gridrotsin_t(i,j))*mask_t(i,j,k))*x1   &
!!         ,0.5*(                        &
!!           -( vel_u(i  ,j  ,k,1)                       &
!!             +vel_u(i+1,j  ,k,1))*gridrotsin_t(i,j)    &
!!           +( vel_v(i  ,j  ,k,1)                       &
!!             +vel_v(i  ,j+1,k,1))*gridrotcos_t(i,j))*mask_t(i,j,k)
!         ,(0.5*(                        &
!           -( vel_u(i  ,j  ,k+1,1)                       &
!             +vel_u(i+1,j  ,k+1,1))*gridrotsin_t(i,j)    &
!           +( vel_v(i  ,j  ,k+1,1)                       &
!             +vel_v(i  ,j+1,k+1,1))*gridrotcos_t(i,j))*mask_t(i,j,k+1))*x0 + &
!           (0.5*(         &
!           -( vel_u(i  ,j  ,k,1)                       &
!             +vel_u(i+1,j  ,k,1))*gridrotsin_t(i,j)    &
!           +( vel_v(i  ,j  ,k,1)                       &
!             +vel_v(i  ,j+1,k,1))*gridrotcos_t(i,j))*mask_t(i,j,k))*x1  &
!         ,anyvar3d(i,j,k+1)*x0+anyvar3d(i,j,k)*x1 &
!         ,tem_t(i,j,k+1,1),x0,tem_t(i,j,k,1),x1,deck,int(deck),nint(deck),k 

          else

!         write(4,'(3(1x,e14.7),12(1x,i4),15(1x,e14.7),3(1x,i4))') &
          write(4,'(3(1x,e14.7),12(1x,i4),7(1x,e14.7))') &
          lon_*rad2deg,lat_*rad2deg,z_,year_glider,month_glider,day_glider,hour_glider,minute_glider,second_glider &
                                      ,year_now   ,month_now   ,day_now ,hour_now   ,minute_now   ,second_now    &
        ,tem_t(i,j,k,1),sal_t(i,j,k,1) &
        ,bio_t(i,j,k,isynechl) &
        ,bio_t(i,j,k,inanochl) &
         ,bio_t(i,j,k,idiachl) &  
         ,bio_t(i,j,k,initrate) &
         ,bio_t(i,j,k,ioxygen) !&
!         ,0.5*(                        &
!            ( vel_u(i  ,j  ,k,1)                       &
!             +vel_u(i+1,j  ,k,1))*gridrotcos_t(i,j)    &
!          +( vel_v(i  ,j  ,k,1)                       &
!             +vel_v(i  ,j+1,k,1))*gridrotsin_t(i,j))*mask_t(i,j,k) &
!         ,0.5*(                        &
!           -( vel_u(i  ,j  ,k,1)                       &
!             +vel_u(i+1,j  ,k,1))*gridrotsin_t(i,j)    &
!           +( vel_v(i  ,j  ,k,1)                       &
!             +vel_v(i  ,j+1,k,1))*gridrotcos_t(i,j))*mask_t(i,j,k)
!         ,anyvar3d(i,j,k) &
!         ,tem_t(i,j,k,1),x0,tem_t(i,j,k,1),x1,deck,int(deck),nint(deck),k

         endif


         close(4)
       endif                                                    !ooo>

       endif



       if(elapsedtime_glider_now_>elapsedtime_now+dti_fw) then
        backspace 543
        goto 3136     ! SORTIR
       else
        goto 3135     ! REVENIR AU DEBUT POUR LIRE LA DONNEE GLIDER SUIVANTE
       endif




 3136 continue
      end subroutine my_outputs_glider
!#endif

!............................................................................
#ifdef bidon
! Version de base ne distinguant pas le signe des flux
      subroutine my_outputs_zone1salttempflux(txt_,id_flux_) !13-04-20
      implicit none
      integer id_flux_
      character(len=*)txt_

! Details dans: https://docs.google.com/document/d/1bI1DeIjmaf8DNYUME7K8KsVX_lh4bnu5_ChKXJHXXlM/edit
! Details dans: https://docs.google.com/presentation/d/17KQ5fzQmqYz7Vp7yzZ6er4elShwLni25GOGIgFpBS3g/edit#slide=id.g516dd41227_0_6


! Initialisation
!      if(iteration3d==0) then !-initial->
!      if(.not.allocated(zone1_mask)) then !-initial->
       if(txt_=='init') then !-initial-> !05-07-19

! Nombre de couches
        zone1_nlayer=3

! Inverse ("1 sur") de l'epaisseur des couches du bilan
       zone1_inv_dz=1./50.
! Cas particulier zone1_nlayer=1: zone1_inv_dz=1./hmax

!      zone1_stretch_dz=1.  ! epaisseur constante = 1/zone1_inv_dz
       zone1_stretch_dz=0.5 ! epaisseur augmentant avec la profondeur (en surface 1/zone1_inv_dz)

! ICI CHACUN SE DEBROUILLE POUR LE LIRE SON FICHIER DE MASQUE AVEC LE
! BON FORMAT.....
!        allocate(glob_mask(0:iglb+1,0:jglb+1)) ; glob_mask=0
!        zone1_max=maxval(glob_mask)
!        allocate(zone1saltflux_glb(0:zone1_max)) ; zone1saltflux_glb=0.
!        allocate(zone1_mask(0:imax+1,0:jmax+1)) ; zone1_mask=0
!        do j=0,jmax+1 ; do i=0,imax+1
!          zone1_mask(i,j)=glob_mask(i+par%timax(1),j+par%tjmax(1))*mask_t(i,j,kmax)
!        enddo       ; enddo
!      deallocate(glob_mask)
! BIDOUILLE PATRICK
         allocate(zone1_mask(0:imax+1,0:jmax+1)) ; zone1_mask=0

!#ifdef bidon

! ATTENTION DE DEFINIR DE i=0 A i=iglb+1 et de j=0 A j=jglb+1
! sinon (par ex de i=1 A iglb) les flux aux OBC ne seront pas calculEs
! et le bilan ne sera pas équilibrer

! Zone -1 : 
         do j1=0,jglb+1 ; do i1=0,iglb+1
          if(j1>777.and.i1<615) then
            i=i1-par%timax(1)
            j=j1-par%tjmax(1)
            if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then
             zone1_mask(i,j)= (-1)*mask_t(i,j,kmax)
            endif
           endif
         enddo            ; enddo

! Zone  1 : 
         do j1=0,jglb+1 ; do i1=0,iglb+1
          if(i1>=615) then
            i=i1-par%timax(1)
            j=j1-par%tjmax(1)
            if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then
             zone1_mask(i,j)=    1*mask_t(i,j,kmax)
            endif
           endif
         enddo            ; enddo

! Zone 2 : 
         do j1=0,jglb+1 ; do i1=0,iglb+1
          if(j1<=777) then
            i=i1-par%timax(1)
            j=j1-par%tjmax(1)
            if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then
             zone1_mask(i,j)=    2*mask_t(i,j,kmax)
            endif
           endif
         enddo            ; enddo


!#endif

         k0=maxval(zone1_mask)
         call mpi_allreduce(k0,zone1_max,1,mpi_integer,mpi_max,par%comm2d ,ierr)
         allocate(zone1saltflux_glb (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1saltflux_glb=0.
         allocate(zone1saltflux_u   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_u=0.
         allocate(zone1saltflux_v   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_v=0.
         allocate(zone1tempflux_glb (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1tempflux_glb=0.
         allocate(zone1tempflux_u   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_u=0.
         allocate(zone1tempflux_v   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_v=0.
         allocate(zone1waterflux_glb(0:zone1_max+1,0:zone1_nlayer-1))   ; zone1waterflux_glb=0. !08-07-19
         allocate(zone1waterflux_u  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_u=0.
         allocate(zone1waterflux_v  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_v=0.

         allocate(zone1tempcumul_glb (0:zone1_nlayer-1))   ; zone1tempcumul_glb=0. !22-09-20 v289
         allocate(zone1tempcumul_loc (0:zone1_nlayer-1))   ; zone1tempcumul_loc=0.
         allocate(zone1tempmasst0    (0:zone1_nlayer-1))   ; zone1tempmasst0=0.

         allocate(zone1saltcumul_glb (0:zone1_nlayer-1))   ; zone1saltcumul_glb=0.
         allocate(zone1saltcumul_loc (0:zone1_nlayer-1))   ; zone1saltcumul_loc=0.
         allocate(zone1saltmasst0    (0:zone1_nlayer-1))   ; zone1saltmasst0=0.

         allocate(zone1watercumul_glb (0:zone1_nlayer-1))   ; zone1watercumul_glb=0.
         allocate(zone1watercumul_loc (0:zone1_nlayer-1))   ; zone1watercumul_loc=0.
         allocate(zone1watermasst0    (0:zone1_nlayer-1))   ; zone1watermasst0=0.
         
        
! La zone1 interieure est -1, les autres >=0.
! Identifier points de flux u
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax   ; do i=1,imax+1  
          if( (zone1_mask(i,j)==-1.and.zone1_mask(i-1,j)/=-1) &
          .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone1_u_max=k
         if(zone1_u_max>0) then !ppp>
          allocate(zone1_flux_u_node(zone1_u_max,3)) ; zone1_flux_u_node=0
          k=0 ! passage 2 pour renseigner tableau zone1_flux_u_node=0
          do j=1,jmax   ; do i=1,imax+1  
           if( (zone1_mask(i,j)==-1.and.zone1_mask(i-1,j)/=-1) &
           .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1) then !pmx>
            k=k+1
            zone1_flux_u_node(k,1)=i 
            zone1_flux_u_node(k,2)=j
            zone1_flux_u_node(k,3)=max(zone1_mask(i-1,j),zone1_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

! Identifier points de flux v
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax+1   ; do i=1,imax
          if( (zone1_mask(i,j)==-1.and.zone1_mask(i,j-1)/=-1) &
          .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone1_v_max=k
         if(zone1_v_max>0) then !ppp>
          allocate(zone1_flux_v_node(zone1_v_max,3)) ; zone1_flux_v_node=0
          k=0 ! passage 2 pour renseigner tableau zone1_flux_v_node=0
          do j=1,jmax+1   ; do i=1,imax
           if( (zone1_mask(i,j)==-1.and.zone1_mask(i,j-1)/=-1) &
           .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1) then !pmx>
            k=k+1
            zone1_flux_v_node(k,1)=i 
            zone1_flux_v_node(k,2)=j
            zone1_flux_v_node(k,3)=max(zone1_mask(i,j-1),zone1_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

      if(par%rank==0) then !#mpi-->>-->               
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine my_outputs_zone1salttempflux'
      write(3,*)'zone1_nlayer=',zone1_nlayer
!     write(3,*)'epaisseur des classes',1./zone1_inv_dz
      write(3,*)'No de classe ,  Zsup(m)  ,  Zinf(m)'
      k=0
      write(3,*)k,' surface      ',-real(k+1)**(1./zone1_stretch_dz)/zone1_inv_dz
      do k=1,zone1_nlayer-2
      write(3,*)k,-real(k  )**(1./zone1_stretch_dz)/zone1_inv_dz  &
                 ,-real(k+1)**(1./zone1_stretch_dz)/zone1_inv_dz
      enddo
      k=zone1_nlayer-1
      write(3,*)k,-real(k)**(1./zone1_stretch_dz)/zone1_inv_dz,-hmax
      close(3)
      endif                !#mpi-->>-->               


       return
       endif                   !-initial->


       if(txt_=='is') then !-flux-u-salt->
          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                           abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20 v289
                          ),0),zone1_nlayer-1)
            zone1saltflux_u(k2,k3)=zone1saltflux_u(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
           enddo
          enddo ! k1
       return
       endif            !-flux-u-salt->

       if(txt_=='it') then !-flux-u-temp->
          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1tempflux_u(k2,k3)=zone1tempflux_u(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
           enddo
          enddo ! k1
       return
       endif            !-flux-u-temp->

       if(txt_=='js') then !-flux-v-salt->
          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_v(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1saltflux_v(k2,k3)=zone1saltflux_v(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
           enddo
          enddo ! k1
       return
       endif            !-flux-v-salt->

       if(txt_=='jt') then !-flux-v-temp->
          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_v(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1tempflux_v(k2,k3)=zone1tempflux_v(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
           enddo
          enddo ! k1
       return
       endif            !-flux-v-temp->

       if(txt_=='k') then !-flux-w->
         do j=1,jmax ; do i=1,imax 
! Note: en surface omega se partage entre une partie implicite liee au wet/dry utilisant une condition !06-03-21
!       de gradient nul T ou S(:,:kmax,2)
!       et une partie explicite omega_evaprec_w pour laquelle Flux=0 pour S et flux T utilisant T(:,:,kmax+1,2)
!       et au fond le flux de sel est nul quelque soit omega et flux=omega(k=1)*T(:,:,0,2) si source sous marine

! salt
           zone1saltflux_w=zone1saltflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *sal_t(i,j,kmax,2)*(omega_w(i,j,kmax+1,1)-0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))) !isole w_surf wetdry

! temp
           zone1tempflux_w=zone1tempflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)*( & !m°v°m>
           ( & !ooo>
            tem_t(i,j,kmax  ,2)*(omega_w(i,j,kmax+1,1)-0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))) & !isole w_surf wetdry
           +tem_t(i,j,kmax+1,2)*(                      0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))) & !isole w_surf river
           -tem_t(i,j,0,2)*omega_w(i,j,1,1) &
           ) & !ooo>
                             -1./( & !pmx>                   ! 2.Dt / !24-04-16
                                   (rhp_t(i,j,kmax)+rho)    &! Surface density
                        *(4190.-5.7*sal_t(i,j,kmax,before)) &! Cp (heat *capacity)
                                 ) & !pmx>
               *(                  snsf_w(i,j,1)   &
                                  +slhf_w(i,j,1)   &
                                  +sshf_w(i,j,1)   &
                +(1.-albedo_w(i,j))*ssr_w(i,j,1)   &
                             +heatrelax_w(i,j,1))  &
                                                       )   !m°v°m>
         enddo       ; enddo
       return
       endif            !-flux-w->

       if(txt_=='fluxbar') then !-water-fluxes->

! note que pour etre coherent avec l'algo pour les flux de T et S (ou
! les tableaux des flux prennent la vitesse avec un signe -, on somme
! ici egalement avec un signe - devant fluxbar (dans le calcul de x0 en
! suivant)

!        x0=-dti_fw/real(iteration2d_max_now)
         x0=-dti_fw

          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
!           zone1waterflux_u(k2)=zone1waterflux_u(k2)+k0*x0*fluxbar_sumt_u(i,j,1)
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1waterflux_u(k2,k3)=zone1waterflux_u(k2,k3)+k0*x0*veldydz_u(i,j,k,1)
           enddo ! k
          enddo ! k1

          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
!           zone1waterflux_v(k2)=zone1waterflux_v(k2)+k0*x0*fluxbar_sumt_v(i,j,1)
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_v(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1waterflux_v(k2,k3)=zone1waterflux_v(k2,k3)+k0*x0*veldxdz_v(i,j,k,1)
           enddo ! k
          enddo ! k1

         do j=1,jmax ; do i=1,imax
           zone1waterflux_w=zone1waterflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *(omega_w(i,j,kmax+1,1)-omega_w(i,j,1,1))
         enddo       ; enddo

       return
       endif                    !-water-fluxes->

! mpi sum
       if(txt_=='mpi') then !--mpi-->

        do k3=0,zone1_nlayer-1

         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x1=   zone1saltflux_u(k2,k3)
           x3=   zone1tempflux_u(k2,k3)
           x5=  zone1waterflux_u(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x1=x1+ zone1saltflux_v(k2,k3)
           x3=x3+ zone1tempflux_v(k2,k3)
           x5=x5+zone1waterflux_v(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1saltflux_glb(k2,k3)=zone1saltflux_glb(k2,k3)+x2
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1tempflux_glb(k2,k3)=zone1tempflux_glb(k2,k3)+x4
          call mpi_allreduce(x5,x6,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb(k2,k3)=zone1waterflux_glb(k2,k3)+x6

         enddo ! k2

        enddo ! k3


         k3=0          ! k3 est la classe de profondeur de la couche de surface
         k2=zone1_max+1 ! echange vertical A travers la surface (justifiE si zone1s assechees)
         x1=zone1saltflux_w
         call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         zone1saltflux_glb(k2,k3)=zone1saltflux_glb(k2,k3)+x2
         x1=zone1tempflux_w
         call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         zone1tempflux_glb(k2,k3)=zone1tempflux_glb(k2,k3)+x2
         x1=zone1waterflux_w
         call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         zone1waterflux_glb(k2,k3)=zone1waterflux_glb(k2,k3)+x2

! reset des tableaux locaux
          zone1saltflux_u=0. ;  zone1saltflux_v=0. ;  zone1saltflux_w=0.
          zone1tempflux_u=0. ;  zone1tempflux_v=0. ;  zone1tempflux_w=0.
         zone1waterflux_u=0. ; zone1waterflux_v=0. ; zone1waterflux_w=0.

! integrale ssh
         zone1watercumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax

! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1watercumul_loc(k3)= &
         zone1watercumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)               +dz_t(i,j,k,1)               )

         enddo       ; enddo       ; enddo

         call mpi_allreduce(zone1watercumul_loc,zone1watercumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1watermasst0(k)=-zone1watercumul_glb(k)+sum(zone1waterflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcwaterflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1watermasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1watercumul_glb)+sum(zone1watermasst0) & ! colonne 2
                     ,sum(zone1waterflux_glb)                        & ! colonne 3
                     ,zone1waterflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1watercumul_glb(0:zone1_nlayer-1)+zone1watermasst0(0:zone1_nlayer-1)
           close(3)
         endif                !-rank0->

! integrale salinite
         zone1saltcumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax

! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1saltcumul_loc(k3)= &
         zone1saltcumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)*sal_t(i,j,k,2)+dz_t(i,j,k,1)*sal_t(i,j,k,1))

         enddo       ; enddo       ; enddo

         call mpi_allreduce(zone1saltcumul_loc,zone1saltcumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1saltmasst0(k)=-zone1saltcumul_glb(k)+sum(zone1saltflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20

         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcsaltflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1saltmasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                 &
                      real(elapsedtime_now/86400.)               & ! colonne 1
                     ,sum(zone1saltcumul_glb)+sum(zone1saltmasst0) & ! colonne 2
                     ,sum(zone1saltflux_glb)                      & ! colonne 3
                     ,zone1saltflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1saltcumul_glb(0:zone1_nlayer-1)+zone1saltmasst0(0:zone1_nlayer-1)
           close(3)
         endif                !-rank0->

! integrale temperature
         zone1tempcumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1tempcumul_loc(k3)= &
         zone1tempcumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)*tem_t(i,j,k,2)+dz_t(i,j,k,1)*tem_t(i,j,k,1))

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone1tempcumul_loc,zone1tempcumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1tempmasst0(k)=-zone1tempcumul_glb(k)+sum(zone1tempflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obctempflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1saltmasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                  &
                      real(elapsedtime_now/86400.)               & ! colonne 1
                     ,sum(zone1tempcumul_glb)+sum(zone1tempmasst0) & ! colonne 2
                     ,sum(zone1tempflux_glb)                      & ! colonne 3
                     ,zone1tempflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1tempcumul_glb(0:zone1_nlayer-1)+zone1tempmasst0(0:zone1_nlayer-1)
           close(3)
         endif                !-rank0->

       return
       endif               !--mpi-->

       stop 'Err undefined txt_ in my_outputs_zone1saltflux'

      end subroutine my_outputs_zone1salttempflux
#endif
!.........................................................................
#ifdef bidon
! VERSION FAISANT LA DISTINCTION SUR LE SIGNE (ENTRANT/SORTANT) DES FLUX
! adaptee de la version basique par Marine Herrmann 18-12-20     
      subroutine my_outputs_zone1salttempflux(txt_,id_flux_) !13-04-20
      implicit none
      integer id_flux_
      character(len=*)txt_

! Details dans: https://docs.google.com/document/d/1bI1DeIjmaf8DNYUME7K8KsVX_lh4bnu5_ChKXJHXXlM/edit
! Details dans: https://docs.google.com/presentation/d/17KQ5fzQmqYz7Vp7yzZ6er4elShwLni25GOGIgFpBS3g/edit#slide=id.g516dd41227_0_6


! Initialisation
!      if(iteration3d==0) then !-initial->
!      if(.not.allocated(zone1_mask)) then !-initial->
       if(txt_=='init') then !-initial-> !05-07-19

! Nombre de couches
        zone1_nlayer=3

! Inverse ("1 sur") de l'epaisseur des couches du bilan
       zone1_inv_dz=1./50.
! Cas particulier zone1_nlayer=1: zone1_inv_dz=1./hmax

!      zone1_stretch_dz=1.  ! epaisseur constante = 1/zone1_inv_dz
       zone1_stretch_dz=0.5 ! epaisseur augmentant avec la profondeur (en surface 1/zone1_inv_dz)


!#ifdef bidon


! lecture des fichiers de zone mask dans des fichiers ijh, meme
! structure que bathycote
! ATTENTION DE DEFINIR DE i=0 A i=iglb+1 et de j=0 A j=jglb+1
! sinon (par ex de i=1 A iglb) les flux aux OBC ne seront pas calculEs
! et le bilan ne sera pas équilibrer

       allocate(glob_mask(0:iglb+1,0:jglb+1)) ; glob_mask=0
       allocate(zone1_mask(0:imax+1,0:jmax+1)) ; zone1_mask=0


       texte250='zmask_SCS.ijh'

       write(6,'(a,a)')'read boxes for flux in ',trim(texte250)
       open(unit=3,file=texte250,recl=10000)
       do i=0,iglb+1
         read(3,'(100000i1)')(glob_mask(i,j),j=0,jglb+1)
       enddo
       close(3)
          
       do j1=0,jglb+1
       do i1=0,iglb+1
            i=i1-par%timax(1)
            j=j1-par%tjmax(1)
            if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then
             zone1_mask(i,j)= glob_mask(i1,j1)-1
           endif
       enddo
       enddo
 
       deallocate(glob_mask)

!#endif

         k0=maxval(zone1_mask)
         call mpi_allreduce(k0,zone1_max,1,mpi_integer,mpi_max,par%comm2d ,ierr)
         allocate(zone1saltflux_glb (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1saltflux_glb=0.
         allocate(zone1saltflux_u   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_u=0.
         allocate(zone1saltflux_v   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_v=0.
         allocate(zone1tempflux_glb (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1tempflux_glb=0.
         allocate(zone1tempflux_u   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_u=0.
         allocate(zone1tempflux_v   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_v=0.
         allocate(zone1waterflux_glb(0:zone1_max+1,0:zone1_nlayer-1))   ; zone1waterflux_glb=0. !08-07-19
         allocate(zone1waterflux_u  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_u=0.
         allocate(zone1waterflux_v  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_v=0.

         allocate(zone1tempcumul_glb (0:zone1_nlayer-1))   ; zone1tempcumul_glb=0. !22-09-20 v289
         allocate(zone1tempcumul_loc (0:zone1_nlayer-1))   ; zone1tempcumul_loc=0.
         allocate(zone1tempmasst0    (0:zone1_nlayer-1))   ; zone1tempmasst0=0.

         allocate(zone1saltcumul_glb (0:zone1_nlayer-1))   ; zone1saltcumul_glb=0.
         allocate(zone1saltcumul_loc (0:zone1_nlayer-1))   ; zone1saltcumul_loc=0.
         allocate(zone1saltmasst0    (0:zone1_nlayer-1))   ; zone1saltmasst0=0.

         allocate(zone1watercumul_glb (0:zone1_nlayer-1))   ; zone1watercumul_glb=0.
         allocate(zone1watercumul_loc (0:zone1_nlayer-1))   ; zone1watercumul_loc=0.
         allocate(zone1watermasst0    (0:zone1_nlayer-1))   ; zone1watermasst0=0.

         allocate(zone1saltflux_glb_in (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1saltflux_glb_in=0.
         allocate(zone1saltflux_u_in   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_u_in=0.
         allocate(zone1saltflux_v_in   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_v_in=0.
         allocate(zone1tempflux_glb_in (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1tempflux_glb_in=0.
         allocate(zone1tempflux_u_in   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_u_in=0.
         allocate(zone1tempflux_v_in   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_v_in=0.
         allocate(zone1waterflux_glb_in(0:zone1_max+1,0:zone1_nlayer-1))   ; zone1waterflux_glb_in=0. 
         allocate(zone1waterflux_u_in  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_u_in=0.
         allocate(zone1waterflux_v_in  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_v_in=0.

         allocate(zone1saltflux_glb_out (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1saltflux_glb_out=0.
         allocate(zone1saltflux_u_out   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_u_out=0.
         allocate(zone1saltflux_v_out   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1saltflux_v_out=0.
         allocate(zone1tempflux_glb_out (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1tempflux_glb_out=0.
         allocate(zone1tempflux_u_out   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_u_out=0.
         allocate(zone1tempflux_v_out   (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1tempflux_v_out=0.
         allocate(zone1waterflux_glb_out(0:zone1_max+1,0:zone1_nlayer-1))   ; zone1waterflux_glb_out=0. 
         allocate(zone1waterflux_u_out  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_u_out=0.
         allocate(zone1waterflux_v_out  (0:zone1_max  ,0:zone1_nlayer-1))   ; zone1waterflux_v_out=0.

         
        
! La zone1 interieure est -1, les autres >=0.
! Identifier points de flux u
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax   ; do i=1,imax+1  
          if( (zone1_mask(i,j)==-1.and.zone1_mask(i-1,j)/=-1) &
          .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone1_u_max=k
         if(zone1_u_max>0) then !ppp>
          allocate(zone1_flux_u_node(zone1_u_max,3)) ; zone1_flux_u_node=0
          k=0 ! passage 2 pour renseigner tableau zone1_flux_u_node=0
          do j=1,jmax   ; do i=1,imax+1  
           if( (zone1_mask(i,j)==-1.and.zone1_mask(i-1,j)/=-1) &
           .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1) then !pmx>
            k=k+1
            zone1_flux_u_node(k,1)=i 
            zone1_flux_u_node(k,2)=j
            zone1_flux_u_node(k,3)=max(zone1_mask(i-1,j),zone1_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

! Identifier points de flux v
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax+1   ; do i=1,imax
          if( (zone1_mask(i,j)==-1.and.zone1_mask(i,j-1)/=-1) &
          .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone1_v_max=k
         if(zone1_v_max>0) then !ppp>
          allocate(zone1_flux_v_node(zone1_v_max,3)) ; zone1_flux_v_node=0
          k=0 ! passage 2 pour renseigner tableau zone1_flux_v_node=0
          do j=1,jmax+1   ; do i=1,imax
           if( (zone1_mask(i,j)==-1.and.zone1_mask(i,j-1)/=-1) &
           .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1) then !pmx>
            k=k+1
            zone1_flux_v_node(k,1)=i 
            zone1_flux_v_node(k,2)=j
            zone1_flux_v_node(k,3)=max(zone1_mask(i,j-1),zone1_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

      if(par%rank==0) then !#mpi-->>-->               
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine my_outputs_zone1salttempflux'
      write(3,*)'zone1_nlayer=',zone1_nlayer
!     write(3,*)'epaisseur des classes',1./zone1_inv_dz
      write(3,*)'No de classe ,  Zsup(m)  ,  Zinf(m)'
      k=0
      write(3,*)k,' surface      ',-real(k+1)**(1./zone1_stretch_dz)/zone1_inv_dz
      do k=1,zone1_nlayer-2
      write(3,*)k,-real(k  )**(1./zone1_stretch_dz)/zone1_inv_dz  &
                 ,-real(k+1)**(1./zone1_stretch_dz)/zone1_inv_dz
      enddo
      k=zone1_nlayer-1
      write(3,*)k,-real(k)**(1./zone1_stretch_dz)/zone1_inv_dz,-hmax
      close(3)
      endif                !#mpi-->>-->               


       return
       endif                   !-initial->


       if(txt_=='is') then !-flux-u-salt->
          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                           abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20 v289
                          ),0),zone1_nlayer-1)
            zone1saltflux_u(k2,k3)=zone1saltflux_u(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)

            if (dti_lpsub*k0*anyv3d(i,j,k,id_flux_) > 0) then
            zone1saltflux_u_in(k2,k3) =zone1saltflux_u_in(k2,k3) +0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            else
            zone1saltflux_u_out(k2,k3)=zone1saltflux_u_out(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            endif

           enddo
          enddo ! k1
       return
       endif            !-flux-u-salt->

       if(txt_=='it') then !-flux-u-temp->
          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1tempflux_u(k2,k3)=zone1tempflux_u(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)

            if (dti_lpsub*k0*anyv3d(i,j,k,id_flux_) > 0) then
            zone1tempflux_u_in(k2,k3) =zone1tempflux_u_in(k2,k3) +0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            else
            zone1tempflux_u_out(k2,k3)=zone1tempflux_u_out(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            endif

           enddo
          enddo ! k1
       return
       endif            !-flux-u-temp->

       if(txt_=='js') then !-flux-v-salt->
          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_v(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1saltflux_v(k2,k3)=zone1saltflux_v(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)

            if (dti_lpsub*k0*anyv3d(i,j,k,id_flux_) > 0) then
            zone1saltflux_v_in(k2,k3) =zone1saltflux_v_in(k2,k3) +0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            else
            zone1saltflux_v_out(k2,k3)=zone1saltflux_v_out(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            endif

           enddo
          enddo ! k1
       return
       endif            !-flux-v-salt->

       if(txt_=='jt') then !-flux-v-temp->
          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_v(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1tempflux_v(k2,k3)=zone1tempflux_v(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)

            if (dti_lpsub*k0*anyv3d(i,j,k,id_flux_) > 0) then
            zone1tempflux_v_in(k2,k3) =zone1tempflux_v_in(k2,k3) +0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            else
            zone1tempflux_v_out(k2,k3)=zone1tempflux_v_out(k2,k3)+0.5*dti_lpsub*k0*anyv3d(i,j,k,id_flux_)
            endif

           enddo
          enddo ! k1
       return
       endif            !-flux-v-temp->

       if(txt_=='k') then !-flux-w->
         do j=1,jmax ; do i=1,imax 
! Note: en surface omega se partage entre une partie implicite liee au wet/dry utilisant une condition !06-03-21
!       de gradient nul T ou S(:,:kmax,2)
!       et une partie explicite omega_evaprec_w pour laquelle Flux=0 pour S et flux T utilisant T(:,:,kmax+1,2)
!       et au fond le flux de sel est nul quelque soit omega et flux=omega(k=1)*T(:,:,0,2) si source sous marine

! salt
           zone1saltflux_w=zone1saltflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *sal_t(i,j,kmax,2)*(omega_w(i,j,kmax+1,1)-0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))) !isole w_surf wetdry

! temp
           zone1tempflux_w=zone1tempflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)*( & !m°v°m>
           ( & !ooo>
            tem_t(i,j,kmax  ,2)*(omega_w(i,j,kmax+1,1)-0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))) & !isole w_surf wetdry
           +tem_t(i,j,kmax+1,2)*(                      0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))) & !isole w_surf river
           -tem_t(i,j,0,2)*omega_w(i,j,1,1) &
           ) & !ooo>
                             -1./( & !pmx>                   ! 2.Dt / !24-04-16
                                   (rhp_t(i,j,kmax)+rho)    &! Surface density
                        *(4190.-5.7*sal_t(i,j,kmax,before)) &! Cp (heat *capacity)
                                 ) & !pmx>
               *(                  snsf_w(i,j,1)   &
                                  +slhf_w(i,j,1)   &
                                  +sshf_w(i,j,1)   &
                +(1.-albedo_w(i,j))*ssr_w(i,j,1)   &
                             +heatrelax_w(i,j,1))  &
                                                       )   !m°v°m>
         enddo       ; enddo
       return
       endif            !-flux-w->

       if(txt_=='fluxbar') then !-water-fluxes->

! note que pour etre coherent avec l'algo pour les flux de T et S (ou
! les tableaux des flux prennent la vitesse avec un signe -, on somme
! ici egalement avec un signe - devant fluxbar (dans le calcul de x0 en
! suivant)

!        x0=-dti_fw/real(iteration2d_max_now)
         x0=-dti_fw

          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1waterflux_u(k2,k3)=zone1waterflux_u(k2,k3)+k0*x0*veldydz_u(i,j,k,1)


            if (k0*x0*veldydz_u(i,j,k,1) > 0) then
            zone1waterflux_u_in(k2,k3) = zone1waterflux_u_in(k2,k3)+k0*x0*veldydz_u(i,j,k,1)
            else
            zone1waterflux_u_out(k2,k3)= zone1waterflux_u_out(k2,k3)+k0*x0*veldydz_u(i,j,k,1)
            endif

           enddo ! k
          enddo ! k1

          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_v(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1waterflux_v(k2,k3)=zone1waterflux_v(k2,k3)+k0*x0*veldxdz_v(i,j,k,1)


            if (k0*x0*veldxdz_v(i,j,k,1) > 0) then
            zone1waterflux_v_in(k2,k3) = zone1waterflux_v_in(k2,k3)+k0*x0*veldxdz_v(i,j,k,1)
            else
            zone1waterflux_v_out(k2,k3)= zone1waterflux_v_out(k2,k3)+k0*x0*veldxdz_v(i,j,k,1)
            endif

           enddo ! k
          enddo ! k1

         do j=1,jmax ; do i=1,imax
           zone1waterflux_w=zone1waterflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *(omega_w(i,j,kmax+1,1)-omega_w(i,j,1,1))
         enddo       ; enddo

       return
       endif                    !-water-fluxes->

! mpi sum
       if(txt_=='mpi') then !--mpi-->

        do k3=0,zone1_nlayer-1

         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x1=   zone1saltflux_u(k2,k3)
           x3=   zone1tempflux_u(k2,k3)
           x5=  zone1waterflux_u(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x1=x1+ zone1saltflux_v(k2,k3)
           x3=x3+ zone1tempflux_v(k2,k3)
           x5=x5+zone1waterflux_v(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1saltflux_glb(k2,k3)=zone1saltflux_glb(k2,k3)+x2
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1tempflux_glb(k2,k3)=zone1tempflux_glb(k2,k3)+x4
          call mpi_allreduce(x5,x6,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb(k2,k3)=zone1waterflux_glb(k2,k3)+x6



          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x1=   zone1saltflux_u_in(k2,k3)
           x3=   zone1tempflux_u_in(k2,k3)
           x5=  zone1waterflux_u_in(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x1=x1+ zone1saltflux_v_in(k2,k3)
           x3=x3+ zone1tempflux_v_in(k2,k3)
           x5=x5+zone1waterflux_v_in(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1saltflux_glb_in(k2,k3)=zone1saltflux_glb_in(k2,k3)+x2
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1tempflux_glb_in(k2,k3)=zone1tempflux_glb_in(k2,k3)+x4
          call mpi_allreduce(x5,x6,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb_in(k2,k3)=zone1waterflux_glb_in(k2,k3)+x6



          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x1=   zone1saltflux_u_out(k2,k3)
           x3=   zone1tempflux_u_out(k2,k3)
           x5=  zone1waterflux_u_out(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x1=x1+ zone1saltflux_v_out(k2,k3)
           x3=x3+ zone1tempflux_v_out(k2,k3)
           x5=x5+zone1waterflux_v_out(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1saltflux_glb_out(k2,k3)=zone1saltflux_glb_out(k2,k3)+x2
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1tempflux_glb_out(k2,k3)=zone1tempflux_glb_out(k2,k3)+x4
          call mpi_allreduce(x5,x6,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb_out(k2,k3)=zone1waterflux_glb_out(k2,k3)+x6

         enddo ! k2

        enddo ! k3


         k3=0          ! k3 est la classe de profondeur de la couche de surface
         k2=zone1_max+1 ! echange vertical A travers la surface (justifiE si zone1s assechees)
         x1=zone1saltflux_w
         call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         zone1saltflux_glb(k2,k3)=zone1saltflux_glb(k2,k3)+x2
         x1=zone1tempflux_w
         call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         zone1tempflux_glb(k2,k3)=zone1tempflux_glb(k2,k3)+x2
         x1=zone1waterflux_w
         call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         zone1waterflux_glb(k2,k3)=zone1waterflux_glb(k2,k3)+x2

! reset des tableaux locaux
         zone1saltflux_u=0. ;  zone1saltflux_v=0. ;  zone1saltflux_w=0.
         zone1tempflux_u=0. ;  zone1tempflux_v=0. ;  zone1tempflux_w=0.
         zone1waterflux_u=0. ; zone1waterflux_v=0. ; zone1waterflux_w=0.

         zone1saltflux_u_in=0. ;  zone1saltflux_v_in=0. ; 
         zone1tempflux_u_in=0. ;  zone1tempflux_v_in=0. ;
         zone1waterflux_u_in=0. ; zone1waterflux_v_in=0. ; 

         zone1saltflux_u_out=0. ;  zone1saltflux_v_out=0. ;  
         zone1tempflux_u_out=0. ;  zone1tempflux_v_out=0. ; 
         zone1waterflux_u_out=0. ; zone1waterflux_v_out=0. ; 

! integrale ssh
         zone1watercumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax

! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1watercumul_loc(k3)= &
         zone1watercumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)               +dz_t(i,j,k,1)               )

         enddo       ; enddo       ; enddo

         call mpi_allreduce(zone1watercumul_loc,zone1watercumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1watermasst0(k)=-zone1watercumul_glb(k)+sum(zone1waterflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcwaterflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1watermasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1watercumul_glb)+sum(zone1watermasst0) & ! colonne 2
                     ,sum(zone1waterflux_glb)                        & ! colonne 3
                     ,zone1waterflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1watercumul_glb(0:zone1_nlayer-1)+zone1watermasst0(0:zone1_nlayer-1)
           close(3)

           write(texte30,'(a)')'tmp/obcwaterflux_zone1_in'
           open(unit=3,file=texte30,position='append')
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1waterflux_glb_in)                        & ! colonne 2
                     ,zone1waterflux_glb_in(0:zone1_max+1,0:zone1_nlayer-1) 
           close(3)

           write(texte30,'(a)')'tmp/obcwaterflux_zone1_out'
           open(unit=3,file=texte30,position='append')
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1waterflux_glb_out)                        & ! colonne 2
                     ,zone1waterflux_glb_out(0:zone1_max+1,0:zone1_nlayer-1) 
           close(3)


         endif                !-rank0->

! integrale salinite
         zone1saltcumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax

! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1saltcumul_loc(k3)= &
         zone1saltcumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)*sal_t(i,j,k,2)+dz_t(i,j,k,1)*sal_t(i,j,k,1))

         enddo       ; enddo       ; enddo

         call mpi_allreduce(zone1saltcumul_loc,zone1saltcumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1saltmasst0(k)=-zone1saltcumul_glb(k)+sum(zone1saltflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20

         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcsaltflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1saltmasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                 &
                      real(elapsedtime_now/86400.)               & ! colonne 1
                     ,sum(zone1saltcumul_glb)+sum(zone1saltmasst0) & ! colonne 2
                     ,sum(zone1saltflux_glb)                      & ! colonne 3
                     ,zone1saltflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1saltcumul_glb(0:zone1_nlayer-1)+zone1saltmasst0(0:zone1_nlayer-1)
           close(3)


           write(texte30,'(a)')'tmp/obcsaltflux_zone1_in'
           open(unit=3,file=texte30,position='append')
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1saltflux_glb_in)                        & ! colonne 2
                     ,zone1saltflux_glb_in(0:zone1_max+1,0:zone1_nlayer-1) 
           close(3)

           write(texte30,'(a)')'tmp/obcsaltflux_zone1_out'
           open(unit=3,file=texte30,position='append')
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1saltflux_glb_out)                        & ! colonne 2
                     ,zone1saltflux_glb_out(0:zone1_max+1,0:zone1_nlayer-1) 
           close(3)

         endif                !-rank0->

! integrale temperature
         zone1tempcumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1tempcumul_loc(k3)= &
         zone1tempcumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)*tem_t(i,j,k,2)+dz_t(i,j,k,1)*tem_t(i,j,k,1))

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone1tempcumul_loc,zone1tempcumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1tempmasst0(k)=-zone1tempcumul_glb(k)+sum(zone1tempflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obctempflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1saltmasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                  &
                      real(elapsedtime_now/86400.)               & ! colonne 1
                     ,sum(zone1tempcumul_glb)+sum(zone1tempmasst0) & ! colonne 2
                     ,sum(zone1tempflux_glb)                      & ! colonne 3
                     ,zone1tempflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1tempcumul_glb(0:zone1_nlayer-1)+zone1tempmasst0(0:zone1_nlayer-1)
           close(3)



           write(texte30,'(a)')'tmp/obctempflux_zone1_in'
           open(unit=3,file=texte30,position='append')
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1tempflux_glb_in)                        & ! colonne 2
                     ,zone1tempflux_glb_in(0:zone1_max+1,0:zone1_nlayer-1) 
           close(3)

           write(texte30,'(a)')'tmp/obctempflux_zone1_out'
           open(unit=3,file=texte30,position='append')
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1tempflux_glb_out)                        & ! colonne 2
                     ,zone1tempflux_glb_out(0:zone1_max+1,0:zone1_nlayer-1) 
           close(3)

         endif                !-rank0->

       return
       endif               !--mpi-->

       stop 'Err undefined txt_ in my_outputs_zone1saltflux'

      end subroutine my_outputs_zone1salttempflux
#endif
!.....................................................................
#ifdef bidon
      subroutine my_outputs_timeabovewater !02-03-21
      implicit none

      if(.not.allocated(timeabovewater_w)) then
       allocate(timeabovewater_w(0:imax+1,0:jmax+1)) ; timeabovewater_w=0.
      endif

! On cumule le temps des points asseches (c.a.d hz_w<0 )
      x1=dti_fw/86400.
!     timeabovewater_count=timeabovewater_count+1

! Si timeabovewater_w
      do j=0,jmax+1 ; do i=0,imax+1
      timeabovewater_w(i,j)= &
      timeabovewater_w(i,j)+x1*max(0.,sign(1.,-(ssh_w(i,j,1)+h_w(i,j)-0.02)))

!     timeabovewater_w(i,j)= &
!     timeabovewater_w(i,j)+max(0.,sign(1.,-(ssh_w(i,j,1)+h_w(i,j)-0.02)))

!      if(i+par%timax(1)==403.and.j+par%tjmax(1)==88) &
!      write(10+par%rank,*)real(elapsedtime_now/86400.) &
!        ,real(timeabovewater_w(i,j)) &
!        ,real(ssh_w(i,j,1))      &
!        ,real(h_w(i,j))

!      if(i+par%timax(1)==31.and.j+par%tjmax(1)==211) &
!      write(10+par%rank,*)real(elapsedtime_now/86400.) &
!        ,real(timeabovewater_w(i,j)) &
!        ,real(ssh_w(i,j,1))      &
!        ,real(h_w(i,j))

      enddo       ; enddo

      end subroutine my_outputs_timeabovewater
#endif
!.....................................................................
#ifdef bilanbio
! Version de base ne distinguant pas le signe des flux
      subroutine my_outputs_zone1bioflux(txt_,id_bio_) !13-04-20
      implicit none
      integer id_bio_
      character(len=*)txt_
      integer,dimension(:,:),allocatable :: glob_mask
      logical :: ldinmesh
      integer :: npoints_poly,i3,j3
      real ,dimension(:)  ,ALLOCATABLE :: xp,yp
      double precision :: x,y



! Details dans: https://docs.google.com/document/d/1bI1DeIjmaf8DNYUME7K8KsVX_lh4bnu5_ChKXJHXXlM/edit
! Details dans: https://docs.google.com/presentation/d/17KQ5fzQmqYz7Vp7yzZ6er4elShwLni25GOGIgFpBS3g/edit#slide=id.g516dd41227_0_6

! Initialisation
!      if(iteration3d==0) then !-initial->
!      if(.not.allocated(zone1_mask)) then !-initial->
       if(txt_=='init') then !-initial-> !05-07-19

! Nombre de couches
        zone1_nlayer=3
!       zone1_nlayer=2
!       zone1_nlayer=1

! Inverse ("1 sur") de l'epaisseur des couches du bilan
       zone1_inv_dz=1./150.
! Cas particulier zone1_nlayer=1: zone1_inv_dz=1./hmax

!      zone1_stretch_dz=1.  ! epaisseur constante = 1/zone1_inv_dz
       zone1_stretch_dz=0.7065 !0.414 ! epaisseur augmentant avec la profondeur (en surface 1/zone1_inv_dz)

! ICI CHACUN SE DEBROUILLE POUR LE LIRE SON FICHIER DE MASQUE AVEC LE
! BON FORMAT.....
         allocate(glob_mask(0:iglb+1,0:jglb+1)) ; glob_mask=0
         allocate(zone1_mask(0:imax+1,0:jmax+1)) ; zone1_mask=0

#ifdef in_out_flux
      open(unit=12,file='/tmpdir/culses/bassin2/GRAPHIQUES/simu2/mask_mldDCAf1000_budget2t.dat')
        do j=0,jglb+1
        read(12,*)(glob_mask(i,j),i=0,iglb+1)
        enddo
      close(12)

! Zone -1
        do j=0,jmax+1 ; do i=0,imax+1
        if(i+par%timax(1)<920) then
        zone1_mask(i,j)=(-1)*glob_mask(i+par%timax(1),j+par%tjmax(1))*mask_t(i,j,kmax)
        endif 
        enddo         ; enddo
! ATTENTION DE DEFINIR DE i=0 A i=iglb+1 et de j=0 A j=jglb+1
! sinon (par ex de i=1 A iglb) les flux aux OBC ne seront pas calculEs
! et le bilan ne sera pas équilibrer

! Zone  1 : 
        do j=0,jmax+1 ; do i=0,imax+1
            if(zone1_mask(i,j).ne.-1) then
             zone1_mask(i,j)=    1*mask_t(i,j,kmax)
            endif
         enddo            ; enddo
#endif

      open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_levantin.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))

        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone1_mask(i1,j1)=-1*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)


       open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_egee.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))


        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone1_mask(i1,j1)=2*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)

         
!        do j=0,jmax+1 ; do i=0,imax+1
!         if(i+par%timax(1)>478) then
!         if(i+par%timax(1)>565) then
!          zone1_mask(i,j)=mask_t(i,j,kmax)   ! zone 1
!         else
!          zone1_mask(i,j)=2*mask_t(i,j,kmax) ! zone 2
!         endif
!        enddo         ; enddo
!        do j=0,jmax+1 ; do i=0,imax+1
!         if(i+par%timax(1)>468.and.i+par%timax(1)<488.and. &
!            j+par%tjmax(1)>586.and.j+par%tjmax(1)<606) then
!         if(i+par%timax(1)>565.and.j+par%tjmax(1)>817) then
!          zone1_mask(i,j)=-1*mask_t(i,j,kmax) ! zone -1
!         endif
!        enddo         ; enddo

! Patrick
!         do j=0,jmax+1 ; do i=0,imax+1
!          if(j+par%tjmax(1)>800) then
!           zone1_mask(i,j)=-1*mask_t(i,j,kmax) ! zone -1

!          do j=0,jmax+1 ; do i=0,imax+1
!          if(i+par%timax(1).eq.953.and.j+par%tjmax(1).eq.125) then
!          print*,'passe ici5',zone1_mask(i,j)
!          endif
!          enddo ; enddo

          do j=0,jmax+1 ; do i=0,imax+1
          if(zone1_mask(i,j).ne.-1*mask_t(i,j,kmax).and.zone1_mask(i,j).ne.2*mask_t(i,j,kmax)) then
           zone1_mask(i,j)=mask_t(i,j,kmax)   ! zone 1
          endif
         enddo         ; enddo

         k0=maxval(zone1_mask)
         call mpi_allreduce(k0,zone1_max,1,mpi_integer,mpi_max,par%comm2d ,ierr)
         allocate(zone1bioflux_glb  (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_glb=0.
         allocate(zone1bioflux_u    (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_u=0.
         allocate(zone1bioflux_v    (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_v=0.

         allocate(zone1biocumul_glb   (0:zone1_nlayer-1))       ; zone1biocumul_glb=0. !22-09-20 v289
         allocate(zone1biocumul_loc   (0:zone1_nlayer-1))       ; zone1biocumul_loc=0.
         allocate(zone1biomasst0      (0:zone1_nlayer-1,vbmax)) ; zone1biomasst0=0.
         allocate(zone1tendancebio_glb(0:zone1_nlayer-1,vbmax)) ; zone1tendancebio_glb=0.
         allocate(zone1botsurfbio_glb (0:zone1_nlayer-1,vbmax)) ; zone1botsurfbio_glb=0.

#ifdef in_out_flux
         allocate(zone1bioflux_glb_in (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_glb_in=0.
         allocate(zone1bioflux_u_in (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_u_in=0.
         allocate(zone1bioflux_v_in (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_v_in=0.

         allocate(zone1bioflux_glb_out (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_glb_out=0.
         allocate(zone1bioflux_u_out (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_u_out=0.
         allocate(zone1bioflux_v_out (0:zone1_max,0:zone1_nlayer-1,vbmax))   ; zone1bioflux_v_out=0.
#endif

         allocate(zone1waterflux_glb (0:zone1_max+1,0:zone1_nlayer-1))   ; zone1waterflux_glb=0.
         allocate(zone1waterflux_u (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_u=0.
         allocate(zone1waterflux_v (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_v=0.

#ifdef in_out_flux
         allocate(zone1waterflux_glb_in (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_glb_in=0.
         allocate(zone1waterflux_u_in (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_u_in=0.
         allocate(zone1waterflux_v_in (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_v_in=0.

         allocate(zone1waterflux_glb_out (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_glb_out=0.
         allocate(zone1waterflux_u_out (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_u_out=0.
         allocate(zone1waterflux_v_out (0:zone1_max,0:zone1_nlayer-1))   ; zone1waterflux_v_out=0.
#endif

         allocate(zone1watercumul_glb   (0:zone1_nlayer-1)) ; zone1watercumul_glb=0. !22-09-20 v289
         allocate(zone1watercumul_loc   (0:zone1_nlayer-1)) ; zone1watercumul_loc=0.
         allocate(zone1watermasst0      (0:zone1_nlayer-1)) ; zone1watermasst0=0.

! La zone1 interieure est -1, les autres >=0.
! Identifier points de flux u
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax   ; do i=1,imax+1  
          if( (zone1_mask(i,j)==-1.and.zone1_mask(i-1,j)/=-1) &
          .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone1_u_max=k
         if(zone1_u_max>0) then !ppp>
          allocate(zone1_flux_u_node(zone1_u_max,3)) ; zone1_flux_u_node=0
          k=0 ! passage 2 pour renseigner tableau zone1_flux_u_node=0
          do j=1,jmax   ; do i=1,imax+1  
           if( (zone1_mask(i,j)==-1.and.zone1_mask(i-1,j)/=-1) &
           .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1) then !pmx>
            k=k+1
            zone1_flux_u_node(k,1)=i 
            zone1_flux_u_node(k,2)=j
            zone1_flux_u_node(k,3)=max(zone1_mask(i-1,j),zone1_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

! Identifier points de flux v
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax+1   ; do i=1,imax
          if( (zone1_mask(i,j)==-1.and.zone1_mask(i,j-1)/=-1) &
          .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone1_v_max=k
         if(zone1_v_max>0) then !ppp>
          allocate(zone1_flux_v_node(zone1_v_max,3)) ; zone1_flux_v_node=0
          k=0 ! passage 2 pour renseigner tableau zone1_flux_v_node=0
          do j=1,jmax+1   ; do i=1,imax
           if( (zone1_mask(i,j)==-1.and.zone1_mask(i,j-1)/=-1) &
           .or.(zone1_mask(i,j)/=-1.and.zone1_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1) then !pmx>
            k=k+1
            zone1_flux_v_node(k,1)=i 
            zone1_flux_v_node(k,2)=j
            zone1_flux_v_node(k,3)=max(zone1_mask(i,j-1),zone1_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

      if(par%rank==0) then !#mpi-->>-->               
      print*,'dir',trim(tmpdirname)//'messages'
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine my_outputs_zone1bioflux'
      write(3,*)'zone1_nlayer=',zone1_nlayer
!     write(3,*)'epaisseur des classes',1./zone1_inv_dz
      write(3,*)'No de classe ,  Zsup(m)  ,  Zinf(m)'
      k=0
      write(3,*)k,' surface      ',-real(k+1)**(1./zone1_stretch_dz)/zone1_inv_dz
      do k=1,zone1_nlayer-2
      write(3,*)k,-real(k  )**(1./zone1_stretch_dz)/zone1_inv_dz  &
                 ,-real(k+1)**(1./zone1_stretch_dz)/zone1_inv_dz
      enddo
      k=zone1_nlayer-1
      write(3,*)k,-real(k)**(1./zone1_stretch_dz)/zone1_inv_dz,-hmax
      close(3)
      endif                !#mpi-->>-->               



       return
       endif                   !-initial->


       if(txt_=='i_bio') then !-flux-u-bio->
          do k1=1,zone1_u_max
            i=zone1_flux_u_node(k1,1)
            j=zone1_flux_u_node(k1,2)
           k2=zone1_flux_u_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
!                           abs(depth_u(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
     abs(0.5*(depth_t(i-1,j,k)+depth_t(i,j,k))*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)
            zone1bioflux_u(k2,k3,vb)= &
            zone1bioflux_u(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                                )   !ooo>
           if(vb.eq.1) then !vb=1>
            zone1waterflux_u(k2,k3)= &
            zone1waterflux_u(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                   ) !ooo>
           endif            !vb=1>

#ifdef in_out_flux
            if((anyv3d(i,j,k,1)+anyv3d(i,j,k,2)+anyv3d(i,j,k,3)+anyv3d(i,j,k,4)).gt.0.D0) then
            zone1bioflux_u_in(k2,k3,vb)= &
            zone1bioflux_u_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                       )   !ooo>
                if(vb.eq.1) &
            zone1waterflux_u_in(k2,k3)= &
            zone1waterflux_u_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                      )  !ooo> 
             else
            zone1bioflux_u_out(k2,k3,vb)= &
            zone1bioflux_u_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                        )   !ooo>
                if(vb.eq.1) &
            zone1waterflux_u_out(k2,k3)= &
            zone1waterflux_u_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                       )  !ooo>

             endif !positif/negatif
#endif

           enddo !k
          enddo  ! k1
       return
       endif            !-flux-u-bio->

       if(txt_=='j_bio') then !-flux-v-bio->
          do k1=1,zone1_v_max
            i=zone1_flux_v_node(k1,1)
            j=zone1_flux_v_node(k1,2)
           k2=zone1_flux_v_node(k1,3) ! numero de la zone1 adjacente
           k0=sign(1,zone1_mask(i,j)-zone1_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
     abs(0.5*(depth_t(i,j-1,k)+depth_t(i,j,k))*zone1_inv_dz)**zone1_stretch_dz & !22-09-20   
                       ),0),zone1_nlayer-1)
            zone1bioflux_v(k2,k3,vb)= &
            zone1bioflux_v(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                    )   !ooo>
              if(vb.eq.1) then !AAA>
            zone1waterflux_v(k2,k3)= &
            zone1waterflux_v(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                   ) !ooo>

               endif           !AAA>

#ifdef in_out_flux
            if((anyv3d(i,j,k,5)+anyv3d(i,j,k,6)+anyv3d(i,j,k,7)+anyv3d(i,j,k,8)).gt.0.D0) then
            zone1bioflux_v_in(k2,k3,vb)= &
            zone1bioflux_v_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                       )   !ooo>
              if(vb.eq.1) &
            zone1waterflux_v_in(k2,k3)= &
            zone1waterflux_v_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                      ) !ooo>

             else
            zone1bioflux_v_out(k2,k3,vb)= &
            zone1bioflux_v_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                        )       !ooo>
              if(vb.eq.1) &
            zone1waterflux_v_out(k2,k3)= &
            zone1waterflux_v_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                       ) !ooo>
             endif
#endif
    
           enddo
          enddo ! k1
       return
       endif            !-flux-v-bio->

       if(txt_=='botsurf') then !-flux-w->
       do vb=1,vbmax
! integrale +dti_fw*fluxbio_w(i,j,vb,2)/dz_t(i,j,k,1)
! Note: attention operation repetee A chaque iteration principale
         zone1biocumul_loc=0.
! SURFACE:
         k=kmax
         do j=1,jmax ; do i=1,imax 
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1biocumul_loc(k3)= &
         zone1biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*fluxbio_w(i,j,vb,2) &
                                                          /dz_t(i,j,k,1)
         enddo       ; enddo
! FOND: (note: si plusieurs layer alors flux de surface et de fond sont distinguEs
         do j=1,jmax ; do i=1,imax 
         k=kmin_w(i,j)
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1biocumul_loc(k3)= &
         zone1biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*(fluxbio_w(i,j,vb,1)+wsed(1,vb)*bio_t(i,j,1,vb)) &
                                                           /dz_t(i,j,k,1)
         enddo       ; enddo
         call mpi_allreduce(zone1biocumul_loc,zone1biocumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone1_nlayer-1
         zone1botsurfbio_glb(k3,vb)= &
         zone1botsurfbio_glb(k3,vb)  &
           +zone1biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb

       do j=1,jmax ; do i=1,imax
           zone1waterflux_w=zone1waterflux_w  &
                         -max(-zone1_mask(i,j),0) & ! selectionne zone1_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *(omega_w(i,j,kmax+1,1)-omega_w(i,j,1,1))
       enddo       ; enddo

       return
       endif                   !-flux-w->

! tendancebio
       if(txt_=='tendancebio') then !--tendancebio-->
       do vb=1,vbmax

! integrale +dti_fw*tendancebio_t(i,j,k,vb)
! Note: attention operation repetee A chaque iteration principale
         zone1biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1biocumul_loc(k3)= &
         zone1biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*tendancebio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone1biocumul_loc,zone1biocumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone1_nlayer-1
         zone1tendancebio_glb(k3,vb)= &
         zone1tendancebio_glb(k3,vb)  &
           +zone1biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb
       return
       endif                        !--tendancebio-->

! mpi sum
       if(txt_=='mpi') then !--mpi-->
       do vb=1,vbmax


        do k3=0,zone1_nlayer-1

         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x3=   zone1bioflux_u(k2,k3,vb)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x3=x3+ zone1bioflux_v(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1bioflux_glb(k2,k3,vb)=zone1bioflux_glb(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3

#ifdef in_out_flux
        do k3=0,zone1_nlayer-1

         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x3=   zone1bioflux_u_in(k2,k3,vb)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x3=x3+ zone1bioflux_v_in(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1bioflux_glb_in(k2,k3,vb)=zone1bioflux_glb_in(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3


        do k3=0,zone1_nlayer-1

         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x3=   zone1bioflux_u_out(k2,k3,vb)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x3=x3+ zone1bioflux_v_out(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1bioflux_glb_out(k2,k3,vb)=zone1bioflux_glb_out(k2,k3,vb)+x4
                            
         enddo ! k2       

        enddo ! k3
#endif

! water
        if(vb.eq.1) then !vb=1>

        do k3=0,zone1_nlayer-1
         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x3=   zone1waterflux_u(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x3=x3+ zone1waterflux_v(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb(k2,k3)=zone1waterflux_glb(k2,k3)+x4

         enddo ! k2
        enddo ! k3

        k3=0           ! k3 est la classe de profondeur de la couche de surface
        k2=zone1_max+1 ! echange vertical A travers la surface
        x1=zone1waterflux_w
        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        zone1waterflux_glb(k2,k3)=zone1waterflux_glb(k2,k3)+x2

#ifdef in_out_flux
! water in
        do k3=0,zone1_nlayer-1
         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x3=   zone1waterflux_u_in(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x3=x3+ zone1waterflux_v_in(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb_in(k2,k3)=zone1waterflux_glb_in(k2,k3)+x4

         enddo ! k2
        enddo ! k3

! water out
        do k3=0,zone1_nlayer-1
         do k2=0,zone1_max

          x1=0. ; x3=0. ; x5=0.
          if(zone1_u_max>0) then !>>>
           x3=   zone1waterflux_u_out(k2,k3)
          endif                 !>>>
          if(zone1_v_max>0) then !>>>
           x3=x3+ zone1waterflux_v_out(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone1waterflux_glb_out(k2,k3)=zone1waterflux_glb_out(k2,k3)+x4

         enddo ! k2
        enddo ! k3
#endif

        endif !vb=1>     end water


!        k3=0           ! k3 est la classe de profondeur de la couche de surface
!        k2=zone1_max+1 ! echange vertical A travers la surface (justifiE si zone1s assechees)
!        x1=zone1bioflux_w
!        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!        zone1bioflux_glb(k2,k3,vb)=zone1bioflux_glb(k2,k3,vb)+x2

! integrale bio
         zone1biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1 pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1biocumul_loc(k3)= &
         zone1biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*bio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone1biocumul_loc,zone1biocumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1biomasst0(k,vb)=-zone1biocumul_glb(k)         &
                              +zone1tendancebio_glb(k,vb)      &
                               +zone1botsurfbio_glb(k,vb)      &
                            +sum(zone1bioflux_glb(:,k,vb))
           enddo
         endif                   !m°v°m> !22-09-20

! Integrale water
         zone1watercumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone1_nlayer-1
! pour z=zmax_zone1_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone1_inv_dz)**zone1_stretch_dz & !22-09-20
                          ),0),zone1_nlayer-1)

         zone1watercumul_loc(k3)= &
         zone1watercumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone1_mask(i,j)) & !=1 si zone1_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone1watercumul_loc,zone1watercumul_glb,zone1_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone1_nlayer-1
            zone1watermasst0(k)=-zone1watercumul_glb(k)         &
                            +sum(zone1waterflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcwaterflux_zone1'
           open(unit=3,file=texte30,position='append')
! Notes: - zone1watermasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone1watercumul_glb)+sum(zone1watermasst0) & ! colonne 2
                     ,sum(zone1waterflux_glb)                        & ! colonne 3
                     ,zone1waterflux_glb(0:zone1_max+1,0:zone1_nlayer-1) &
                     ,zone1watercumul_glb(0:zone1_nlayer-1)+zone1watermasst0(0:zone1_nlayer-1)
           close(3)
         endif                !-rank0->


         if(par%rank==0) then !-rank0->
           sum0=0. ; sum1=0.
           do k3=0,zone1_nlayer-1
             sum0=sum0+zone1biocumul_glb(k3)+zone1biomasst0(k3,vb)
             do k2=0,zone1_max
               sum1=sum1+zone1bioflux_glb(k2,k3,vb)
!              write(666,*)iteration3d,k2,k3,vb,zone1bioflux_glb(k2,k3,vb)
             enddo
             sum1=sum1+zone1tendancebio_glb(k3,vb) &
                       +zone1botsurfbio_glb(k3,vb)  
           enddo
           
!......................
! Ecrire un fichier pour verifier l'equilibrage
           write(texte60,'(a,i0)')'tmp/zone1_bilan_total_bio',vb
           open(unit=3,file=texte60,position='append')
! Notes: - zone1biomasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')    &
                      real(elapsedtime_now/86400.)  & ! colonne 1
                     ,sum0                          & ! colonne 2
                     ,sum1                            ! colonne 3
           close(3)
! Note: la colonne 2 et la colonne 3 doivent etre egales
! colonne 2: variation du contenu du traceur cumulee sur toutes les
!            couches si nlayer>1
! colonne 3: contribution cumulee de tous les termes du bilan (flux
!            lateraux, surface, tendancebio etc....)

!......................
! Ecrire un fichier par couche verticale (numero 0 pour couche de surface)
           do k=0,zone1_nlayer-1
           write(texte60,'(a,i0,a,i0)')'tmp/zone1_bilan_detail_bio',vb &
                                      ,'_layer',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & ! colonne 1
                     ,zone1biocumul_glb(k)+zone1biomasst0(k,vb) & ! colonne 2
                     ,zone1tendancebio_glb(k,vb)                & ! colonne 3
                      ,zone1botsurfbio_glb(k,vb)                & ! colonne 4
                     ,zone1bioflux_glb(0:zone1_max,k,vb)          ! colonne 5=zone0, colonne 6=zone1, etc...
           close(3)


#ifdef in_out_flux
! in/out flux
           write(texte60,'(a,i0,a,i0)')'tmp/zone1_bilan_detail_bio',vb &
                                      ,'_layer_inout',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone1bioflux_glb_in(0:zone1_max,k,vb)     & 
                     ,zone1bioflux_glb_out(0:zone1_max,k,vb)            
           close(3)
! water
        if(vb.eq.1) then
! Ecrire fichier detail
           write(texte60,'(a,i0)')'tmp/zone1_bilan_detail_water_layer_inout',k 
           open(unit=3,file=texte60,position='append')
                      write(3,'(120(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone1watercumul_glb(k)+zone1watermasst0(k) & !colonne 2
                     ,zone1waterflux_glb(0:zone1_max,k)        &
                     ,zone1waterflux_glb_in(0:zone1_max,k)     &
                     ,zone1waterflux_glb_out(0:zone1_max,k)    
           close(3)
         endif ! vb water
#endif

! Note:
! Si nlayer=1 alors la colonne 2 doit etre egale a la somme des colonnes
! suivantes (c.a.d. de la colonne 3 a la derniere colonne).
! Si nlayer>1 alors l'ecart entre la colonne 2 et la somme des autres
! colonnes s'explique par le flux vertical (advectif et diffusif) entre
! les couches (c'est donc comme cela que se deduit le flux vertical
! entre les couches).
! colonne 2: variation du contenu du traceur dans la couche concernee
! colonne 3: contribution de tendancebio_t
! colonne 4: contribution des flux de fond et de surface. Si nlayer=1
!            les 2 sont confondus. Si nlayer>1 alors la colonne 4 ne
!            contient que le flux de surface dans la couche 0 et que le
!            flux de fond dans la couche nlayer-1.
! colonne 5: contribution du flux lateral entre zone -1 et zone 0
!            la zone 0 est le masque continental donc il s'agit des
!            rivieres
! colonne 6: contribution du flux lateral entre zone -1 et zone 1
! colonne 7: etc...
           enddo ! fin de boucle sur k
         endif                !-rank0->

       enddo ! fin de boucle sur vb

! reset des tableaux locaux pour un nouveau cycle de cumul
       zone1bioflux_u=0. 
       zone1bioflux_v=0. 
       zone1bioflux_w=0.
       zone1waterflux_u=0.
       zone1waterflux_v=0.
       zone1waterflux_w=0.
#ifdef in_out_flux
       zone1bioflux_u_in=0.    
       zone1bioflux_v_in=0.  
       zone1bioflux_u_out=0.   
       zone1bioflux_v_out=0.  
       zone1waterflux_u_in=0.  
       zone1waterflux_v_in=0.  
       zone1waterflux_u_out=0. 
       zone1waterflux_v_out=0.   
       zone1waterflux_u_in=0.
       zone1waterflux_v_in=0.
       zone1waterflux_u_out=0.
       zone1waterflux_v_out=0.
#endif
       return
       endif               !--mpi-->

       stop 'Err undefined txt_ in my_outputs_zone1bioflux'

      end subroutine my_outputs_zone1bioflux
#endif


!.....................................................................
#ifdef bilanbio
! Version de base ne distinguant pas le signe des flux
      subroutine my_outputs_zone2bioflux(txt_,id_bio_) !13-04-20
      implicit none
      integer id_bio_
      character(len=*)txt_
      integer,dimension(:,:),allocatable :: glob_mask
      logical :: ldinmesh
      integer :: npoints_poly,i3,j3
      real ,dimension(:)  ,ALLOCATABLE :: xp,yp
      double precision :: x,y



! Details dans: https://docs.google.com/document/d/1bI1DeIjmaf8DNYUME7K8KsVX_lh4bnu5_ChKXJHXXlM/edit
! Details dans: https://docs.google.com/presentation/d/17KQ5fzQmqYz7Vp7yzZ6er4elShwLni25GOGIgFpBS3g/edit#slide=id.g516dd41227_0_6

! Initialisation
!      if(iteration3d==0) then !-initial->
!      if(.not.allocated(zone2_mask)) then !-initial->
       if(txt_=='init') then !-initial-> !05-07-19

! Nombre de couches
        zone2_nlayer=3
!       zone2_nlayer=2
!       zone2_nlayer=1

! Inverse ("1 sur") de l'epaisseur des couches du bilan
       zone2_inv_dz=1./150.
! Cas particulier zone2_nlayer=1: zone2_inv_dz=1./hmax

!      zone2_stretch_dz=1.  ! epaisseur constante = 1/zone2_inv_dz
       zone2_stretch_dz=0.7065 !0.414 ! epaisseur augmentant avec la profondeur (en surface 1/zone2_inv_dz)

! ICI CHACUN SE DEBROUILLE POUR LE LIRE SON FICHIER DE MASQUE AVEC LE
! BON FORMAT.....
         allocate(glob_mask(0:iglb+1,0:jglb+1)) ; glob_mask=0
         allocate(zone2_mask(0:imax+1,0:jmax+1)) ; zone2_mask=0

#ifdef in_out_flux
      open(unit=12,file='/tmpdir/culses/bassin2/GRAPHIQUES/simu2/mask_mldDCAf1000_budget2t.dat')
        do j=0,jglb+1
        read(12,*)(glob_mask(i,j),i=0,iglb+1)
        enddo
      close(12)

! Zone -1
        do j=0,jmax+1 ; do i=0,imax+1
        if(i+par%timax(1)<920) then
        zone2_mask(i,j)=(-1)*glob_mask(i+par%timax(1),j+par%tjmax(1))*mask_t(i,j,kmax)
        endif 
        enddo         ; enddo
! ATTENTION DE DEFINIR DE i=0 A i=iglb+1 et de j=0 A j=jglb+1
! sinon (par ex de i=1 A iglb) les flux aux OBC ne seront pas calculEs
! et le bilan ne sera pas équilibrer

! Zone  1 : 
        do j=0,jmax+1 ; do i=0,imax+1
            if(zone2_mask(i,j).ne.-1) then
             zone2_mask(i,j)=    1*mask_t(i,j,kmax)
            endif
         enddo            ; enddo
#endif

      open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_levantin.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))

        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone2_mask(i1,j1)=2*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)


       open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_egee.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))


        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone2_mask(i1,j1)=-1*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)

         
!        do j=0,jmax+1 ; do i=0,imax+1
!         if(i+par%timax(1)>478) then
!         if(i+par%timax(1)>565) then
!          zone2_mask(i,j)=mask_t(i,j,kmax)   ! zone 1
!         else
!          zone2_mask(i,j)=2*mask_t(i,j,kmax) ! zone 2
!         endif
!        enddo         ; enddo
!        do j=0,jmax+1 ; do i=0,imax+1
!         if(i+par%timax(1)>468.and.i+par%timax(1)<488.and. &
!            j+par%tjmax(1)>586.and.j+par%tjmax(1)<606) then
!         if(i+par%timax(1)>565.and.j+par%tjmax(1)>817) then
!          zone2_mask(i,j)=-1*mask_t(i,j,kmax) ! zone -1
!         endif
!        enddo         ; enddo

! Patrick
!         do j=0,jmax+1 ; do i=0,imax+1
!          if(j+par%tjmax(1)>800) then
!           zone2_mask(i,j)=-1*mask_t(i,j,kmax) ! zone -1

!          do j=0,jmax+1 ; do i=0,imax+1
!          if(i+par%timax(1).eq.953.and.j+par%tjmax(1).eq.125) then
!          print*,'passe ici5',zone2_mask(i,j)
!          endif
!          enddo ; enddo

          do j=0,jmax+1 ; do i=0,imax+1
          if(zone2_mask(i,j).ne.-1*mask_t(i,j,kmax).and.zone2_mask(i,j).ne.2*mask_t(i,j,kmax)) then
           zone2_mask(i,j)=mask_t(i,j,kmax)   ! zone 1
          endif
         enddo         ; enddo

         k0=maxval(zone2_mask)
         call mpi_allreduce(k0,zone2_max,1,mpi_integer,mpi_max,par%comm2d ,ierr)
         allocate(zone2bioflux_glb  (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_glb=0.
         allocate(zone2bioflux_u    (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_u=0.
         allocate(zone2bioflux_v    (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_v=0.

         allocate(zone2biocumul_glb   (0:zone2_nlayer-1))       ; zone2biocumul_glb=0. !22-09-20 v289
         allocate(zone2biocumul_loc   (0:zone2_nlayer-1))       ; zone2biocumul_loc=0.
         allocate(zone2biomasst0      (0:zone2_nlayer-1,vbmax)) ; zone2biomasst0=0.
         allocate(zone2tendancebio_glb(0:zone2_nlayer-1,vbmax)) ; zone2tendancebio_glb=0.
         allocate(zone2botsurfbio_glb (0:zone2_nlayer-1,vbmax)) ; zone2botsurfbio_glb=0.

#ifdef in_out_flux
         allocate(zone2bioflux_glb_in (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_glb_in=0.
         allocate(zone2bioflux_u_in (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_u_in=0.
         allocate(zone2bioflux_v_in (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_v_in=0.

         allocate(zone2bioflux_glb_out (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_glb_out=0.
         allocate(zone2bioflux_u_out (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_u_out=0.
         allocate(zone2bioflux_v_out (0:zone2_max,0:zone2_nlayer-1,vbmax))   ; zone2bioflux_v_out=0.
#endif

         allocate(zone2waterflux_glb (0:zone2_max+1,0:zone2_nlayer-1))   ; zone2waterflux_glb=0.
         allocate(zone2waterflux_u (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_u=0.
         allocate(zone2waterflux_v (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_v=0.

#ifdef in_out_flux
         allocate(zone2waterflux_glb_in (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_glb_in=0.
         allocate(zone2waterflux_u_in (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_u_in=0.
         allocate(zone2waterflux_v_in (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_v_in=0.

         allocate(zone2waterflux_glb_out (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_glb_out=0.
         allocate(zone2waterflux_u_out (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_u_out=0.
         allocate(zone2waterflux_v_out (0:zone2_max,0:zone2_nlayer-1))   ; zone2waterflux_v_out=0.
#endif

         allocate(zone2watercumul_glb   (0:zone2_nlayer-1)) ; zone2watercumul_glb=0. !22-09-20 v289
         allocate(zone2watercumul_loc   (0:zone2_nlayer-1)) ; zone2watercumul_loc=0.
         allocate(zone2watermasst0      (0:zone2_nlayer-1)) ; zone2watermasst0=0.

! La zone2 interieure est -1, les autres >=0.
! Identifier points de flux u
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax   ; do i=1,imax+1  
          if( (zone2_mask(i,j)==-1.and.zone2_mask(i-1,j)/=-1) &
          .or.(zone2_mask(i,j)/=-1.and.zone2_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone2_u_max=k
         if(zone2_u_max>0) then !ppp>
          allocate(zone2_flux_u_node(zone2_u_max,3)) ; zone2_flux_u_node=0
          k=0 ! passage 2 pour renseigner tableau zone2_flux_u_node=0
          do j=1,jmax   ; do i=1,imax+1  
           if( (zone2_mask(i,j)==-1.and.zone2_mask(i-1,j)/=-1) &
           .or.(zone2_mask(i,j)/=-1.and.zone2_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1) then !pmx>
            k=k+1
            zone2_flux_u_node(k,1)=i 
            zone2_flux_u_node(k,2)=j
            zone2_flux_u_node(k,3)=max(zone2_mask(i-1,j),zone2_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

! Identifier points de flux v
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax+1   ; do i=1,imax
          if( (zone2_mask(i,j)==-1.and.zone2_mask(i,j-1)/=-1) &
          .or.(zone2_mask(i,j)/=-1.and.zone2_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone2_v_max=k
         if(zone2_v_max>0) then !ppp>
          allocate(zone2_flux_v_node(zone2_v_max,3)) ; zone2_flux_v_node=0
          k=0 ! passage 2 pour renseigner tableau zone2_flux_v_node=0
          do j=1,jmax+1   ; do i=1,imax
           if( (zone2_mask(i,j)==-1.and.zone2_mask(i,j-1)/=-1) &
           .or.(zone2_mask(i,j)/=-1.and.zone2_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1) then !pmx>
            k=k+1
            zone2_flux_v_node(k,1)=i 
            zone2_flux_v_node(k,2)=j
            zone2_flux_v_node(k,3)=max(zone2_mask(i,j-1),zone2_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

      if(par%rank==0) then !#mpi-->>-->               
      print*,'dir',trim(tmpdirname)//'messages'
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine my_outputs_zone2bioflux'
      write(3,*)'zone2_nlayer=',zone2_nlayer
!     write(3,*)'epaisseur des classes',1./zone2_inv_dz
      write(3,*)'No de classe ,  Zsup(m)  ,  Zinf(m)'
      k=0
      write(3,*)k,' surface      ',-real(k+1)**(1./zone2_stretch_dz)/zone2_inv_dz
      do k=1,zone2_nlayer-2
      write(3,*)k,-real(k  )**(1./zone2_stretch_dz)/zone2_inv_dz  &
                 ,-real(k+1)**(1./zone2_stretch_dz)/zone2_inv_dz
      enddo
      k=zone2_nlayer-1
      write(3,*)k,-real(k)**(1./zone2_stretch_dz)/zone2_inv_dz,-hmax
      close(3)
      endif                !#mpi-->>-->               



       return
       endif                   !-initial->


       if(txt_=='i_bio') then !-flux-u-bio->
          do k1=1,zone2_u_max
            i=zone2_flux_u_node(k1,1)
            j=zone2_flux_u_node(k1,2)
           k2=zone2_flux_u_node(k1,3) ! numero de la zone2 adjacente
           k0=sign(1,zone2_mask(i,j)-zone2_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1 pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
!                           abs(depth_u(i,j,k)*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
     abs(0.5*(depth_t(i-1,j,k)+depth_t(i,j,k))*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
                          ),0),zone2_nlayer-1)
            zone2bioflux_u(k2,k3,vb)= &
            zone2bioflux_u(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                                )   !ooo>
           if(vb.eq.1) then !vb=1>
            zone2waterflux_u(k2,k3)= &
            zone2waterflux_u(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                   ) !ooo>
           endif            !vb=1>

#ifdef in_out_flux
            if((anyv3d(i,j,k,1)+anyv3d(i,j,k,2)+anyv3d(i,j,k,3)+anyv3d(i,j,k,4)).gt.0.D0) then
            zone2bioflux_u_in(k2,k3,vb)= &
            zone2bioflux_u_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                       )   !ooo>
                if(vb.eq.1) &
            zone2waterflux_u_in(k2,k3)= &
            zone2waterflux_u_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                      )  !ooo> 
             else
            zone2bioflux_u_out(k2,k3,vb)= &
            zone2bioflux_u_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                        )   !ooo>
                if(vb.eq.1) &
            zone2waterflux_u_out(k2,k3)= &
            zone2waterflux_u_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                       )  !ooo>

             endif !positif/negatif
#endif

           enddo !k
          enddo  ! k1
       return
       endif            !-flux-u-bio->

       if(txt_=='j_bio') then !-flux-v-bio->
          do k1=1,zone2_v_max
            i=zone2_flux_v_node(k1,1)
            j=zone2_flux_v_node(k1,2)
           k2=zone2_flux_v_node(k1,3) ! numero de la zone2 adjacente
           k0=sign(1,zone2_mask(i,j)-zone2_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1 pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
     abs(0.5*(depth_t(i,j-1,k)+depth_t(i,j,k))*zone2_inv_dz)**zone2_stretch_dz & !22-09-20   
                       ),0),zone2_nlayer-1)
            zone2bioflux_v(k2,k3,vb)= &
            zone2bioflux_v(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                    )   !ooo>
              if(vb.eq.1) then !AAA>
            zone2waterflux_v(k2,k3)= &
            zone2waterflux_v(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                   ) !ooo>

               endif           !AAA>

#ifdef in_out_flux
            if((anyv3d(i,j,k,5)+anyv3d(i,j,k,6)+anyv3d(i,j,k,7)+anyv3d(i,j,k,8)).gt.0.D0) then
            zone2bioflux_v_in(k2,k3,vb)= &
            zone2bioflux_v_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                       )   !ooo>
              if(vb.eq.1) &
            zone2waterflux_v_in(k2,k3)= &
            zone2waterflux_v_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                      ) !ooo>

             else
            zone2bioflux_v_out(k2,k3,vb)= &
            zone2bioflux_v_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                        )       !ooo>
              if(vb.eq.1) &
            zone2waterflux_v_out(k2,k3)= &
            zone2waterflux_v_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                       ) !ooo>
             endif
#endif
    
           enddo
          enddo ! k1
       return
       endif            !-flux-v-bio->

       if(txt_=='botsurf') then !-flux-w->
       do vb=1,vbmax
! integrale +dti_fw*fluxbio_w(i,j,vb,2)/dz_t(i,j,k,1)
! Note: attention operation repetee A chaque iteration principale
         zone2biocumul_loc=0.
! SURFACE:
         k=kmax
         do j=1,jmax ; do i=1,imax 
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1 pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
                          ),0),zone2_nlayer-1)

         zone2biocumul_loc(k3)= &
         zone2biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone2_mask(i,j)) & !=1 si zone2_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*fluxbio_w(i,j,vb,2) &
                                                          /dz_t(i,j,k,1)
         enddo       ; enddo
! FOND: (note: si plusieurs layer alors flux de surface et de fond sont distinguEs
         do j=1,jmax ; do i=1,imax 
         k=kmin_w(i,j)
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1 pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
                          ),0),zone2_nlayer-1)

         zone2biocumul_loc(k3)= &
         zone2biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone2_mask(i,j)) & !=1 si zone2_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*(fluxbio_w(i,j,vb,1)+wsed(1,vb)*bio_t(i,j,1,vb)) &
                                                           /dz_t(i,j,k,1)
         enddo       ; enddo
         call mpi_allreduce(zone2biocumul_loc,zone2biocumul_glb,zone2_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone2_nlayer-1
         zone2botsurfbio_glb(k3,vb)= &
         zone2botsurfbio_glb(k3,vb)  &
           +zone2biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb

       do j=1,jmax ; do i=1,imax
           zone2waterflux_w=zone2waterflux_w  &
                         -max(-zone2_mask(i,j),0) & ! selectionne zone2_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *(omega_w(i,j,kmax+1,1)-omega_w(i,j,1,1))
       enddo       ; enddo

       return
       endif                   !-flux-w->

! tendancebio
       if(txt_=='tendancebio') then !--tendancebio-->
       do vb=1,vbmax

! integrale +dti_fw*tendancebio_t(i,j,k,vb)
! Note: attention operation repetee A chaque iteration principale
         zone2biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1 pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
                          ),0),zone2_nlayer-1)

         zone2biocumul_loc(k3)= &
         zone2biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone2_mask(i,j)) & !=1 si zone2_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*tendancebio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone2biocumul_loc,zone2biocumul_glb,zone2_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone2_nlayer-1
         zone2tendancebio_glb(k3,vb)= &
         zone2tendancebio_glb(k3,vb)  &
           +zone2biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb
       return
       endif                        !--tendancebio-->

! mpi sum
       if(txt_=='mpi') then !--mpi-->
       do vb=1,vbmax


        do k3=0,zone2_nlayer-1

         do k2=0,zone2_max

          x1=0. ; x3=0. ; x5=0.
          if(zone2_u_max>0) then !>>>
           x3=   zone2bioflux_u(k2,k3,vb)
          endif                 !>>>
          if(zone2_v_max>0) then !>>>
           x3=x3+ zone2bioflux_v(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone2bioflux_glb(k2,k3,vb)=zone2bioflux_glb(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3

#ifdef in_out_flux
        do k3=0,zone2_nlayer-1

         do k2=0,zone2_max

          x1=0. ; x3=0. ; x5=0.
          if(zone2_u_max>0) then !>>>
           x3=   zone2bioflux_u_in(k2,k3,vb)
          endif                 !>>>
          if(zone2_v_max>0) then !>>>
           x3=x3+ zone2bioflux_v_in(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone2bioflux_glb_in(k2,k3,vb)=zone2bioflux_glb_in(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3


        do k3=0,zone2_nlayer-1

         do k2=0,zone2_max

          x1=0. ; x3=0. ; x5=0.
          if(zone2_u_max>0) then !>>>
           x3=   zone2bioflux_u_out(k2,k3,vb)
          endif                 !>>>
          if(zone2_v_max>0) then !>>>
           x3=x3+ zone2bioflux_v_out(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone2bioflux_glb_out(k2,k3,vb)=zone2bioflux_glb_out(k2,k3,vb)+x4
                            
         enddo ! k2       

        enddo ! k3
#endif

! water
        if(vb.eq.1) then !vb=1>

        do k3=0,zone2_nlayer-1
         do k2=0,zone2_max

          x1=0. ; x3=0. ; x5=0.
          if(zone2_u_max>0) then !>>>
           x3=   zone2waterflux_u(k2,k3)
          endif                 !>>>
          if(zone2_v_max>0) then !>>>
           x3=x3+ zone2waterflux_v(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone2waterflux_glb(k2,k3)=zone2waterflux_glb(k2,k3)+x4

         enddo ! k2
        enddo ! k3

        k3=0           ! k3 est la classe de profondeur de la couche de surface
        k2=zone2_max+1 ! echange vertical A travers la surface
        x1=zone2waterflux_w
        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        zone2waterflux_glb(k2,k3)=zone2waterflux_glb(k2,k3)+x2

#ifdef in_out_flux
! water in
        do k3=0,zone2_nlayer-1
         do k2=0,zone2_max

          x1=0. ; x3=0. ; x5=0.
          if(zone2_u_max>0) then !>>>
           x3=   zone2waterflux_u_in(k2,k3)
          endif                 !>>>
          if(zone2_v_max>0) then !>>>
           x3=x3+ zone2waterflux_v_in(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone2waterflux_glb_in(k2,k3)=zone2waterflux_glb_in(k2,k3)+x4

         enddo ! k2
        enddo ! k3

! water out
        do k3=0,zone2_nlayer-1
         do k2=0,zone2_max

          x1=0. ; x3=0. ; x5=0.
          if(zone2_u_max>0) then !>>>
           x3=   zone2waterflux_u_out(k2,k3)
          endif                 !>>>
          if(zone2_v_max>0) then !>>>
           x3=x3+ zone2waterflux_v_out(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone2waterflux_glb_out(k2,k3)=zone2waterflux_glb_out(k2,k3)+x4

         enddo ! k2
        enddo ! k3
#endif

        endif !vb=1>     end water


!        k3=0           ! k3 est la classe de profondeur de la couche de surface
!        k2=zone2_max+1 ! echange vertical A travers la surface (justifiE si zone2s assechees)
!        x1=zone2bioflux_w
!        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!        zone2bioflux_glb(k2,k3,vb)=zone2bioflux_glb(k2,k3,vb)+x2

! integrale bio
         zone2biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1 pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
                          ),0),zone2_nlayer-1)

         zone2biocumul_loc(k3)= &
         zone2biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone2_mask(i,j)) & !=1 si zone2_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*bio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone2biocumul_loc,zone2biocumul_glb,zone2_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone2_nlayer-1
            zone2biomasst0(k,vb)=-zone2biocumul_glb(k)         &
                              +zone2tendancebio_glb(k,vb)      &
                               +zone2botsurfbio_glb(k,vb)      &
                            +sum(zone2bioflux_glb(:,k,vb))
           enddo
         endif                   !m°v°m> !22-09-20

! Integrale water
         zone2watercumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone2_nlayer-1
! pour z=zmax_zone2_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone2_inv_dz)**zone2_stretch_dz & !22-09-20
                          ),0),zone2_nlayer-1)

         zone2watercumul_loc(k3)= &
         zone2watercumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone2_mask(i,j)) & !=1 si zone2_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone2watercumul_loc,zone2watercumul_glb,zone2_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone2_nlayer-1
            zone2watermasst0(k)=-zone2watercumul_glb(k)         &
                            +sum(zone2waterflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcwaterflux_zone2'
           open(unit=3,file=texte30,position='append')
! Notes: - zone2watermasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone2watercumul_glb)+sum(zone2watermasst0) & ! colonne 2
                     ,sum(zone2waterflux_glb)                        & ! colonne 3
                     ,zone2waterflux_glb(0:zone2_max+1,0:zone2_nlayer-1) &
                     ,zone2watercumul_glb(0:zone2_nlayer-1)+zone2watermasst0(0:zone2_nlayer-1)
           close(3)
         endif                !-rank0->


         if(par%rank==0) then !-rank0->
           sum0=0. ; sum1=0.
           do k3=0,zone2_nlayer-1
             sum0=sum0+zone2biocumul_glb(k3)+zone2biomasst0(k3,vb)
             do k2=0,zone2_max
               sum1=sum1+zone2bioflux_glb(k2,k3,vb)
!              write(666,*)iteration3d,k2,k3,vb,zone2bioflux_glb(k2,k3,vb)
             enddo
             sum1=sum1+zone2tendancebio_glb(k3,vb) &
                       +zone2botsurfbio_glb(k3,vb)  
           enddo
           
!......................
! Ecrire un fichier pour verifier l'equilibrage
           write(texte60,'(a,i0)')'tmp/zone2_bilan_total_bio',vb
           open(unit=3,file=texte60,position='append')
! Notes: - zone2biomasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')    &
                      real(elapsedtime_now/86400.)  & ! colonne 1
                     ,sum0                          & ! colonne 2
                     ,sum1                            ! colonne 3
           close(3)
! Note: la colonne 2 et la colonne 3 doivent etre egales
! colonne 2: variation du contenu du traceur cumulee sur toutes les
!            couches si nlayer>1
! colonne 3: contribution cumulee de tous les termes du bilan (flux
!            lateraux, surface, tendancebio etc....)

!......................
! Ecrire un fichier par couche verticale (numero 0 pour couche de surface)
           do k=0,zone2_nlayer-1
           write(texte60,'(a,i0,a,i0)')'tmp/zone2_bilan_detail_bio',vb &
                                      ,'_layer',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & ! colonne 1
                     ,zone2biocumul_glb(k)+zone2biomasst0(k,vb) & ! colonne 2
                     ,zone2tendancebio_glb(k,vb)                & ! colonne 3
                      ,zone2botsurfbio_glb(k,vb)                & ! colonne 4
                     ,zone2bioflux_glb(0:zone2_max,k,vb)          ! colonne 5=zone0, colonne 6=zone2, etc...
           close(3)


#ifdef in_out_flux
! in/out flux
           write(texte60,'(a,i0,a,i0)')'tmp/zone2_bilan_detail_bio',vb &
                                      ,'_layer_inout',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone2bioflux_glb_in(0:zone2_max,k,vb)     & 
                     ,zone2bioflux_glb_out(0:zone2_max,k,vb)            
           close(3)
! water
        if(vb.eq.1) then
! Ecrire fichier detail
           write(texte60,'(a,i0)')'tmp/zone2_bilan_detail_water_layer_inout',k 
           open(unit=3,file=texte60,position='append')
                      write(3,'(120(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone2watercumul_glb(k)+zone2watermasst0(k) & !colonne 2
                     ,zone2waterflux_glb(0:zone2_max,k)        &
                     ,zone2waterflux_glb_in(0:zone2_max,k)     &
                     ,zone2waterflux_glb_out(0:zone2_max,k)    
           close(3)
         endif ! vb water
#endif

! Note:
! Si nlayer=1 alors la colonne 2 doit etre egale a la somme des colonnes
! suivantes (c.a.d. de la colonne 3 a la derniere colonne).
! Si nlayer>1 alors l'ecart entre la colonne 2 et la somme des autres
! colonnes s'explique par le flux vertical (advectif et diffusif) entre
! les couches (c'est donc comme cela que se deduit le flux vertical
! entre les couches).
! colonne 2: variation du contenu du traceur dans la couche concernee
! colonne 3: contribution de tendancebio_t
! colonne 4: contribution des flux de fond et de surface. Si nlayer=1
!            les 2 sont confondus. Si nlayer>1 alors la colonne 4 ne
!            contient que le flux de surface dans la couche 0 et que le
!            flux de fond dans la couche nlayer-1.
! colonne 5: contribution du flux lateral entre zone -1 et zone 0
!            la zone 0 est le masque continental donc il s'agit des
!            rivieres
! colonne 6: contribution du flux lateral entre zone -1 et zone 1
! colonne 7: etc...
           enddo ! fin de boucle sur k
         endif                !-rank0->

       enddo ! fin de boucle sur vb

! reset des tableaux locaux pour un nouveau cycle de cumul
       zone2bioflux_u=0. 
       zone2bioflux_v=0. 
       zone2bioflux_w=0.
       zone2waterflux_u=0.
       zone2waterflux_v=0.
       zone2waterflux_w=0.
#ifdef in_out_flux
       zone2bioflux_u_in=0.    
       zone2bioflux_v_in=0.  
       zone2bioflux_u_out=0.   
       zone2bioflux_v_out=0.  
       zone2waterflux_u_in=0.  
       zone2waterflux_v_in=0.  
       zone2waterflux_u_out=0. 
       zone2waterflux_v_out=0.   
       zone2waterflux_u_in=0.
       zone2waterflux_v_in=0.
       zone2waterflux_u_out=0.
       zone2waterflux_v_out=0.
#endif
       return
       endif               !--mpi-->

       stop 'Err undefined txt_ in my_outputs_zone2bioflux'

      end subroutine my_outputs_zone2bioflux
#endif


!.....................................................................
#ifdef bilanbio
! Version de base ne distinguant pas le signe des flux
      subroutine my_outputs_zone3bioflux(txt_,id_bio_) !13-04-20
      implicit none
      integer id_bio_
      character(len=*)txt_
      integer,dimension(:,:),allocatable :: glob_mask
      logical :: ldinmesh
      integer :: npoints_poly,i3,j3
      real ,dimension(:)  ,ALLOCATABLE :: xp,yp
      double precision :: x,y



! Details dans: https://docs.google.com/document/d/1bI1DeIjmaf8DNYUME7K8KsVX_lh4bnu5_ChKXJHXXlM/edit
! Details dans: https://docs.google.com/presentation/d/17KQ5fzQmqYz7Vp7yzZ6er4elShwLni25GOGIgFpBS3g/edit#slide=id.g516dd41227_0_6

! Initialisation
!      if(iteration3d==0) then !-initial->
!      if(.not.allocated(zone3_mask)) then !-initial->
       if(txt_=='init') then !-initial-> !05-07-19

! Nombre de couches
        zone3_nlayer=3
!       zone3_nlayer=2
!       zone3_nlayer=1

! Inverse ("1 sur") de l'epaisseur des couches du bilan
       zone3_inv_dz=1./150.
! Cas particulier zone3_nlayer=1: zone3_inv_dz=1./hmax

!      zone3_stretch_dz=1.  ! epaisseur constante = 1/zone3_inv_dz
       zone3_stretch_dz=0.7065 !0.414 ! epaisseur augmentant avec la profondeur (en surface 1/zone3_inv_dz)

! ICI CHACUN SE DEBROUILLE POUR LE LIRE SON FICHIER DE MASQUE AVEC LE
! BON FORMAT.....
         allocate(glob_mask(0:iglb+1,0:jglb+1)) ; glob_mask=0
         allocate(zone3_mask(0:imax+1,0:jmax+1)) ; zone3_mask=0

#ifdef in_out_flux
      open(unit=12,file='/tmpdir/culses/bassin2/GRAPHIQUES/simu2/mask_mldDCAf1000_budget2t.dat')
        do j=0,jglb+1
        read(12,*)(glob_mask(i,j),i=0,iglb+1)
        enddo
      close(12)

! Zone -1
        do j=0,jmax+1 ; do i=0,imax+1
        if(i+par%timax(1)<920) then
        zone3_mask(i,j)=(-1)*glob_mask(i+par%timax(1),j+par%tjmax(1))*mask_t(i,j,kmax)
        endif 
        enddo         ; enddo
! ATTENTION DE DEFINIR DE i=0 A i=iglb+1 et de j=0 A j=jglb+1
! sinon (par ex de i=1 A iglb) les flux aux OBC ne seront pas calculEs
! et le bilan ne sera pas équilibrer

! Zone  1 : 
        do j=0,jmax+1 ; do i=0,imax+1
            if(zone3_mask(i,j).ne.-1) then
             zone3_mask(i,j)=    1*mask_t(i,j,kmax)
            endif
         enddo            ; enddo
#endif

      open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_levantin.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))

        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone3_mask(i1,j1)=1*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)


       open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_ionien.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))


        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone3_mask(i1,j1)=-1*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)


       open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_adria.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))


        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone3_mask(i1,j1)=3*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)



       open(unit=12,file='/tmpdir/jhabib/simulations/joelle_profil/POLY_BILAN/poly_egee.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))


        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone3_mask(i1,j1)=2*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)
         
!        do j=0,jmax+1 ; do i=0,imax+1
!         if(i+par%timax(1)>478) then
!         if(i+par%timax(1)>565) then
!          zone3_mask(i,j)=mask_t(i,j,kmax)   ! zone 1
!         else
!          zone3_mask(i,j)=2*mask_t(i,j,kmax) ! zone 2
!         endif
!        enddo         ; enddo
!        do j=0,jmax+1 ; do i=0,imax+1
!         if(i+par%timax(1)>468.and.i+par%timax(1)<488.and. &
!            j+par%tjmax(1)>586.and.j+par%tjmax(1)<606) then
!         if(i+par%timax(1)>565.and.j+par%tjmax(1)>817) then
!          zone3_mask(i,j)=-1*mask_t(i,j,kmax) ! zone -1
!         endif
!        enddo         ; enddo

! Patrick
!         do j=0,jmax+1 ; do i=0,imax+1
!          if(j+par%tjmax(1)>800) then
!           zone3_mask(i,j)=-1*mask_t(i,j,kmax) ! zone -1

!          do j=0,jmax+1 ; do i=0,imax+1
!          if(i+par%timax(1).eq.953.and.j+par%tjmax(1).eq.125) then
!          print*,'passe ici5',zone3_mask(i,j)
!          endif
!          enddo ; enddo

          do j=0,jmax+1 ; do i=0,imax+1
          if(zone3_mask(i,j).ne.-1*mask_t(i,j,kmax).and.zone3_mask(i,j).ne.2*mask_t(i,j,kmax).and. &
             zone3_mask(i,j).ne.3*mask_t(i,j,kmax).and.zone3_mask(i,j).ne.1*mask_t(i,j,kmax)) then
           zone3_mask(i,j)=4*mask_t(i,j,kmax)   ! zone 1
          endif
         enddo         ; enddo

         k0=maxval(zone3_mask)
         call mpi_allreduce(k0,zone3_max,1,mpi_integer,mpi_max,par%comm2d ,ierr)
         allocate(zone3bioflux_glb  (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_glb=0.
         allocate(zone3bioflux_u    (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_u=0.
         allocate(zone3bioflux_v    (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_v=0.

         allocate(zone3biocumul_glb   (0:zone3_nlayer-1))       ; zone3biocumul_glb=0. !22-09-20 v289
         allocate(zone3biocumul_loc   (0:zone3_nlayer-1))       ; zone3biocumul_loc=0.
         allocate(zone3biomasst0      (0:zone3_nlayer-1,vbmax)) ; zone3biomasst0=0.
         allocate(zone3tendancebio_glb(0:zone3_nlayer-1,vbmax)) ; zone3tendancebio_glb=0.
         allocate(zone3botsurfbio_glb (0:zone3_nlayer-1,vbmax)) ; zone3botsurfbio_glb=0.

#ifdef in_out_flux
         allocate(zone3bioflux_glb_in (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_glb_in=0.
         allocate(zone3bioflux_u_in (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_u_in=0.
         allocate(zone3bioflux_v_in (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_v_in=0.

         allocate(zone3bioflux_glb_out (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_glb_out=0.
         allocate(zone3bioflux_u_out (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_u_out=0.
         allocate(zone3bioflux_v_out (0:zone3_max,0:zone3_nlayer-1,vbmax))   ; zone3bioflux_v_out=0.
#endif

         allocate(zone3waterflux_glb (0:zone3_max+1,0:zone3_nlayer-1))   ; zone3waterflux_glb=0.
         allocate(zone3waterflux_u (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_u=0.
         allocate(zone3waterflux_v (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_v=0.

#ifdef in_out_flux
         allocate(zone3waterflux_glb_in (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_glb_in=0.
         allocate(zone3waterflux_u_in (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_u_in=0.
         allocate(zone3waterflux_v_in (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_v_in=0.

         allocate(zone3waterflux_glb_out (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_glb_out=0.
         allocate(zone3waterflux_u_out (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_u_out=0.
         allocate(zone3waterflux_v_out (0:zone3_max,0:zone3_nlayer-1))   ; zone3waterflux_v_out=0.
#endif

         allocate(zone3watercumul_glb   (0:zone3_nlayer-1)) ; zone3watercumul_glb=0. !22-09-20 v289
         allocate(zone3watercumul_loc   (0:zone3_nlayer-1)) ; zone3watercumul_loc=0.
         allocate(zone3watermasst0      (0:zone3_nlayer-1)) ; zone3watermasst0=0.

! La zone3 interieure est -1, les autres >=0.
! Identifier points de flux u
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax   ; do i=1,imax+1  
          if( (zone3_mask(i,j)==-1.and.zone3_mask(i-1,j)/=-1) &
          .or.(zone3_mask(i,j)/=-1.and.zone3_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone3_u_max=k
         if(zone3_u_max>0) then !ppp>
          allocate(zone3_flux_u_node(zone3_u_max,3)) ; zone3_flux_u_node=0
          k=0 ! passage 2 pour renseigner tableau zone3_flux_u_node=0
          do j=1,jmax   ; do i=1,imax+1  
           if( (zone3_mask(i,j)==-1.and.zone3_mask(i-1,j)/=-1) &
           .or.(zone3_mask(i,j)/=-1.and.zone3_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1) then !pmx>
            k=k+1
            zone3_flux_u_node(k,1)=i 
            zone3_flux_u_node(k,2)=j
            zone3_flux_u_node(k,3)=max(zone3_mask(i-1,j),zone3_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

! Identifier points de flux v
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax+1   ; do i=1,imax
          if( (zone3_mask(i,j)==-1.and.zone3_mask(i,j-1)/=-1) &
          .or.(zone3_mask(i,j)/=-1.and.zone3_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone3_v_max=k
         if(zone3_v_max>0) then !ppp>
          allocate(zone3_flux_v_node(zone3_v_max,3)) ; zone3_flux_v_node=0
          k=0 ! passage 2 pour renseigner tableau zone3_flux_v_node=0
          do j=1,jmax+1   ; do i=1,imax
           if( (zone3_mask(i,j)==-1.and.zone3_mask(i,j-1)/=-1) &
           .or.(zone3_mask(i,j)/=-1.and.zone3_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1) then !pmx>
            k=k+1
            zone3_flux_v_node(k,1)=i 
            zone3_flux_v_node(k,2)=j
            zone3_flux_v_node(k,3)=max(zone3_mask(i,j-1),zone3_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

      if(par%rank==0) then !#mpi-->>-->               
      print*,'dir',trim(tmpdirname)//'messages'
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine my_outputs_zone3bioflux'
      write(3,*)'zone3_nlayer=',zone3_nlayer
!     write(3,*)'epaisseur des classes',1./zone3_inv_dz
      write(3,*)'No de classe ,  Zsup(m)  ,  Zinf(m)'
      k=0
      write(3,*)k,' surface      ',-real(k+1)**(1./zone3_stretch_dz)/zone3_inv_dz
      do k=1,zone3_nlayer-2
      write(3,*)k,-real(k  )**(1./zone3_stretch_dz)/zone3_inv_dz  &
                 ,-real(k+1)**(1./zone3_stretch_dz)/zone3_inv_dz
      enddo
      k=zone3_nlayer-1
      write(3,*)k,-real(k)**(1./zone3_stretch_dz)/zone3_inv_dz,-hmax
      close(3)
      endif                !#mpi-->>-->               



       return
       endif                   !-initial->


       if(txt_=='i_bio') then !-flux-u-bio->
          do k1=1,zone3_u_max
            i=zone3_flux_u_node(k1,1)
            j=zone3_flux_u_node(k1,2)
           k2=zone3_flux_u_node(k1,3) ! numero de la zone3 adjacente
           k0=sign(1,zone3_mask(i,j)-zone3_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1 pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
!                           abs(depth_u(i,j,k)*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
     abs(0.5*(depth_t(i-1,j,k)+depth_t(i,j,k))*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
                          ),0),zone3_nlayer-1)
            zone3bioflux_u(k2,k3,vb)= &
            zone3bioflux_u(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                                )   !ooo>
           if(vb.eq.1) then !vb=1>
            zone3waterflux_u(k2,k3)= &
            zone3waterflux_u(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                   ) !ooo>
           endif            !vb=1>

#ifdef in_out_flux
            if((anyv3d(i,j,k,1)+anyv3d(i,j,k,2)+anyv3d(i,j,k,3)+anyv3d(i,j,k,4)).gt.0.D0) then
            zone3bioflux_u_in(k2,k3,vb)= &
            zone3bioflux_u_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                       )   !ooo>
                if(vb.eq.1) &
            zone3waterflux_u_in(k2,k3)= &
            zone3waterflux_u_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                      )  !ooo> 
             else
            zone3bioflux_u_out(k2,k3,vb)= &
            zone3bioflux_u_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                        )   !ooo>
                if(vb.eq.1) &
            zone3waterflux_u_out(k2,k3)= &
            zone3waterflux_u_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                       )  !ooo>

             endif !positif/negatif
#endif

           enddo !k
          enddo  ! k1
       return
       endif            !-flux-u-bio->

       if(txt_=='j_bio') then !-flux-v-bio->
          do k1=1,zone3_v_max
            i=zone3_flux_v_node(k1,1)
            j=zone3_flux_v_node(k1,2)
           k2=zone3_flux_v_node(k1,3) ! numero de la zone3 adjacente
           k0=sign(1,zone3_mask(i,j)-zone3_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1 pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
     abs(0.5*(depth_t(i,j-1,k)+depth_t(i,j,k))*zone3_inv_dz)**zone3_stretch_dz & !22-09-20   
                       ),0),zone3_nlayer-1)
            zone3bioflux_v(k2,k3,vb)= &
            zone3bioflux_v(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                    )   !ooo>
              if(vb.eq.1) then !AAA>
            zone3waterflux_v(k2,k3)= &
            zone3waterflux_v(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                   ) !ooo>

               endif           !AAA>

#ifdef in_out_flux
            if((anyv3d(i,j,k,5)+anyv3d(i,j,k,6)+anyv3d(i,j,k,7)+anyv3d(i,j,k,8)).gt.0.D0) then
            zone3bioflux_v_in(k2,k3,vb)= &
            zone3bioflux_v_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                       )   !ooo>
              if(vb.eq.1) &
            zone3waterflux_v_in(k2,k3)= &
            zone3waterflux_v_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                      ) !ooo>

             else
            zone3bioflux_v_out(k2,k3,vb)= &
            zone3bioflux_v_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                        )       !ooo>
              if(vb.eq.1) &
            zone3waterflux_v_out(k2,k3)= &
            zone3waterflux_v_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                       ) !ooo>
             endif
#endif
    
           enddo
          enddo ! k1
       return
       endif            !-flux-v-bio->

       if(txt_=='botsurf') then !-flux-w->
       do vb=1,vbmax
! integrale +dti_fw*fluxbio_w(i,j,vb,2)/dz_t(i,j,k,1)
! Note: attention operation repetee A chaque iteration principale
         zone3biocumul_loc=0.
! SURFACE:
         k=kmax
         do j=1,jmax ; do i=1,imax 
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1 pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
                          ),0),zone3_nlayer-1)

         zone3biocumul_loc(k3)= &
         zone3biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone3_mask(i,j)) & !=1 si zone3_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*fluxbio_w(i,j,vb,2) &
                                                          /dz_t(i,j,k,1)
         enddo       ; enddo
! FOND: (note: si plusieurs layer alors flux de surface et de fond sont distinguEs
         do j=1,jmax ; do i=1,imax 
         k=kmin_w(i,j)
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1 pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
                          ),0),zone3_nlayer-1)

         zone3biocumul_loc(k3)= &
         zone3biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone3_mask(i,j)) & !=1 si zone3_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*(fluxbio_w(i,j,vb,1)+wsed(1,vb)*bio_t(i,j,1,vb)) &
                                                           /dz_t(i,j,k,1)
         enddo       ; enddo
         call mpi_allreduce(zone3biocumul_loc,zone3biocumul_glb,zone3_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone3_nlayer-1
         zone3botsurfbio_glb(k3,vb)= &
         zone3botsurfbio_glb(k3,vb)  &
           +zone3biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb

       do j=1,jmax ; do i=1,imax
           zone3waterflux_w=zone3waterflux_w  &
                         -max(-zone3_mask(i,j),0) & ! selectionne zone3_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *(omega_w(i,j,kmax+1,1)-omega_w(i,j,1,1))
       enddo       ; enddo

       return
       endif                   !-flux-w->

! tendancebio
       if(txt_=='tendancebio') then !--tendancebio-->
       do vb=1,vbmax

! integrale +dti_fw*tendancebio_t(i,j,k,vb)
! Note: attention operation repetee A chaque iteration principale
         zone3biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1 pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
                          ),0),zone3_nlayer-1)

         zone3biocumul_loc(k3)= &
         zone3biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone3_mask(i,j)) & !=1 si zone3_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*tendancebio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone3biocumul_loc,zone3biocumul_glb,zone3_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone3_nlayer-1
         zone3tendancebio_glb(k3,vb)= &
         zone3tendancebio_glb(k3,vb)  &
           +zone3biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb
       return
       endif                        !--tendancebio-->

! mpi sum
       if(txt_=='mpi') then !--mpi-->
       do vb=1,vbmax


        do k3=0,zone3_nlayer-1

         do k2=0,zone3_max

          x1=0. ; x3=0. ; x5=0.
          if(zone3_u_max>0) then !>>>
           x3=   zone3bioflux_u(k2,k3,vb)
          endif                 !>>>
          if(zone3_v_max>0) then !>>>
           x3=x3+ zone3bioflux_v(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone3bioflux_glb(k2,k3,vb)=zone3bioflux_glb(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3

#ifdef in_out_flux
        do k3=0,zone3_nlayer-1

         do k2=0,zone3_max

          x1=0. ; x3=0. ; x5=0.
          if(zone3_u_max>0) then !>>>
           x3=   zone3bioflux_u_in(k2,k3,vb)
          endif                 !>>>
          if(zone3_v_max>0) then !>>>
           x3=x3+ zone3bioflux_v_in(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone3bioflux_glb_in(k2,k3,vb)=zone3bioflux_glb_in(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3


        do k3=0,zone3_nlayer-1

         do k2=0,zone3_max

          x1=0. ; x3=0. ; x5=0.
          if(zone3_u_max>0) then !>>>
           x3=   zone3bioflux_u_out(k2,k3,vb)
          endif                 !>>>
          if(zone3_v_max>0) then !>>>
           x3=x3+ zone3bioflux_v_out(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone3bioflux_glb_out(k2,k3,vb)=zone3bioflux_glb_out(k2,k3,vb)+x4
                            
         enddo ! k2       

        enddo ! k3
#endif

! water
        if(vb.eq.1) then !vb=1>

        do k3=0,zone3_nlayer-1
         do k2=0,zone3_max

          x1=0. ; x3=0. ; x5=0.
          if(zone3_u_max>0) then !>>>
           x3=   zone3waterflux_u(k2,k3)
          endif                 !>>>
          if(zone3_v_max>0) then !>>>
           x3=x3+ zone3waterflux_v(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone3waterflux_glb(k2,k3)=zone3waterflux_glb(k2,k3)+x4

         enddo ! k2
        enddo ! k3

        k3=0           ! k3 est la classe de profondeur de la couche de surface
        k2=zone3_max+1 ! echange vertical A travers la surface
        x1=zone3waterflux_w
        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        zone3waterflux_glb(k2,k3)=zone3waterflux_glb(k2,k3)+x2

#ifdef in_out_flux
! water in
        do k3=0,zone3_nlayer-1
         do k2=0,zone3_max

          x1=0. ; x3=0. ; x5=0.
          if(zone3_u_max>0) then !>>>
           x3=   zone3waterflux_u_in(k2,k3)
          endif                 !>>>
          if(zone3_v_max>0) then !>>>
           x3=x3+ zone3waterflux_v_in(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone3waterflux_glb_in(k2,k3)=zone3waterflux_glb_in(k2,k3)+x4

         enddo ! k2
        enddo ! k3

! water out
        do k3=0,zone3_nlayer-1
         do k2=0,zone3_max

          x1=0. ; x3=0. ; x5=0.
          if(zone3_u_max>0) then !>>>
           x3=   zone3waterflux_u_out(k2,k3)
          endif                 !>>>
          if(zone3_v_max>0) then !>>>
           x3=x3+ zone3waterflux_v_out(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone3waterflux_glb_out(k2,k3)=zone3waterflux_glb_out(k2,k3)+x4

         enddo ! k2
        enddo ! k3
#endif

        endif !vb=1>     end water


!        k3=0           ! k3 est la classe de profondeur de la couche de surface
!        k2=zone3_max+1 ! echange vertical A travers la surface (justifiE si zone3s assechees)
!        x1=zone3bioflux_w
!        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!        zone3bioflux_glb(k2,k3,vb)=zone3bioflux_glb(k2,k3,vb)+x2

! integrale bio
         zone3biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1 pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
                          ),0),zone3_nlayer-1)

         zone3biocumul_loc(k3)= &
         zone3biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone3_mask(i,j)) & !=1 si zone3_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*bio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone3biocumul_loc,zone3biocumul_glb,zone3_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone3_nlayer-1
            zone3biomasst0(k,vb)=-zone3biocumul_glb(k)         &
                              +zone3tendancebio_glb(k,vb)      &
                               +zone3botsurfbio_glb(k,vb)      &
                            +sum(zone3bioflux_glb(:,k,vb))
           enddo
         endif                   !m°v°m> !22-09-20

! Integrale water
         zone3watercumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone3_nlayer-1
! pour z=zmax_zone3_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone3_inv_dz)**zone3_stretch_dz & !22-09-20
                          ),0),zone3_nlayer-1)

         zone3watercumul_loc(k3)= &
         zone3watercumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone3_mask(i,j)) & !=1 si zone3_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone3watercumul_loc,zone3watercumul_glb,zone3_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone3_nlayer-1
            zone3watermasst0(k)=-zone3watercumul_glb(k)         &
                            +sum(zone3waterflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20
         if(par%rank==0) then !-rank0->
           write(texte30,'(a)')'tmp/obcwaterflux_zone3'
           open(unit=3,file=texte30,position='append')
! Notes: - zone3watermasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')                     &
                      real(elapsedtime_now/86400.)                   & ! colonne 1
                     ,sum(zone3watercumul_glb)+sum(zone3watermasst0) & ! colonne 2
                     ,sum(zone3waterflux_glb)                        & ! colonne 3
                     ,zone3waterflux_glb(0:zone3_max+1,0:zone3_nlayer-1) &
                     ,zone3watercumul_glb(0:zone3_nlayer-1)+zone3watermasst0(0:zone3_nlayer-1)
           close(3)
         endif                !-rank0->


         if(par%rank==0) then !-rank0->
           sum0=0. ; sum1=0.
           do k3=0,zone3_nlayer-1
             sum0=sum0+zone3biocumul_glb(k3)+zone3biomasst0(k3,vb)
             do k2=0,zone3_max
               sum1=sum1+zone3bioflux_glb(k2,k3,vb)
!              write(666,*)iteration3d,k2,k3,vb,zone3bioflux_glb(k2,k3,vb)
             enddo
             sum1=sum1+zone3tendancebio_glb(k3,vb) &
                       +zone3botsurfbio_glb(k3,vb)  
           enddo
           
!......................
! Ecrire un fichier pour verifier l'equilibrage
           write(texte60,'(a,i0)')'tmp/zone3_bilan_total_bio',vb
           open(unit=3,file=texte60,position='append')
! Notes: - zone3biomasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')    &
                      real(elapsedtime_now/86400.)  & ! colonne 1
                     ,sum0                          & ! colonne 2
                     ,sum1                            ! colonne 3
           close(3)
! Note: la colonne 2 et la colonne 3 doivent etre egales
! colonne 2: variation du contenu du traceur cumulee sur toutes les
!            couches si nlayer>1
! colonne 3: contribution cumulee de tous les termes du bilan (flux
!            lateraux, surface, tendancebio etc....)

!......................
! Ecrire un fichier par couche verticale (numero 0 pour couche de surface)
           do k=0,zone3_nlayer-1
           write(texte60,'(a,i0,a,i0)')'tmp/zone3_bilan_detail_bio',vb &
                                      ,'_layer',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & ! colonne 1
                     ,zone3biocumul_glb(k)+zone3biomasst0(k,vb) & ! colonne 2
                     ,zone3tendancebio_glb(k,vb)                & ! colonne 3
                      ,zone3botsurfbio_glb(k,vb)                & ! colonne 4
                     ,zone3bioflux_glb(0:zone3_max,k,vb)          ! colonne 5=zone0, colonne 6=zone3, etc...
           close(3)


#ifdef in_out_flux
! in/out flux
           write(texte60,'(a,i0,a,i0)')'tmp/zone3_bilan_detail_bio',vb &
                                      ,'_layer_inout',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone3bioflux_glb_in(0:zone3_max,k,vb)     & 
                     ,zone3bioflux_glb_out(0:zone3_max,k,vb)            
           close(3)
! water
        if(vb.eq.1) then
! Ecrire fichier detail
           write(texte60,'(a,i0)')'tmp/zone3_bilan_detail_water_layer_inout',k 
           open(unit=3,file=texte60,position='append')
                      write(3,'(120(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone3watercumul_glb(k)+zone3watermasst0(k) & !colonne 2
                     ,zone3waterflux_glb(0:zone3_max,k)        &
                     ,zone3waterflux_glb_in(0:zone3_max,k)     &
                     ,zone3waterflux_glb_out(0:zone3_max,k)    
           close(3)
         endif ! vb water
#endif

! Note:
! Si nlayer=1 alors la colonne 2 doit etre egale a la somme des colonnes
! suivantes (c.a.d. de la colonne 3 a la derniere colonne).
! Si nlayer>1 alors l'ecart entre la colonne 2 et la somme des autres
! colonnes s'explique par le flux vertical (advectif et diffusif) entre
! les couches (c'est donc comme cela que se deduit le flux vertical
! entre les couches).
! colonne 2: variation du contenu du traceur dans la couche concernee
! colonne 3: contribution de tendancebio_t
! colonne 4: contribution des flux de fond et de surface. Si nlayer=1
!            les 2 sont confondus. Si nlayer>1 alors la colonne 4 ne
!            contient que le flux de surface dans la couche 0 et que le
!            flux de fond dans la couche nlayer-1.
! colonne 5: contribution du flux lateral entre zone -1 et zone 0
!            la zone 0 est le masque continental donc il s'agit des
!            rivieres
! colonne 6: contribution du flux lateral entre zone -1 et zone 1
! colonne 7: etc...
           enddo ! fin de boucle sur k
         endif                !-rank0->

       enddo ! fin de boucle sur vb

! reset des tableaux locaux pour un nouveau cycle de cumul
       zone3bioflux_u=0. 
       zone3bioflux_v=0. 
       zone3bioflux_w=0.
       zone3waterflux_u=0.
       zone3waterflux_v=0.
       zone3waterflux_w=0.
#ifdef in_out_flux
       zone3bioflux_u_in=0.    
       zone3bioflux_v_in=0.  
       zone3bioflux_u_out=0.   
       zone3bioflux_v_out=0.  
       zone3waterflux_u_in=0.  
       zone3waterflux_v_in=0.  
       zone3waterflux_u_out=0. 
       zone3waterflux_v_out=0.   
       zone3waterflux_u_in=0.
       zone3waterflux_v_in=0.
       zone3waterflux_u_out=0.
       zone3waterflux_v_out=0.
#endif
       return
       endif               !--mpi-->

       stop 'Err undefined txt_ in my_outputs_zone3bioflux'

      end subroutine my_outputs_zone3bioflux

#endif


!.....................................................................
#ifdef bilanbio
! Version de base ne distinguant pas le signe des flux
      subroutine my_outputs_zone4bioflux(txt_,id_bio_) !13-04-20
      implicit none
      integer id_bio_
      character(len=*)txt_
      integer,dimension(:,:),allocatable :: glob_mask
      logical :: ldinmesh
      integer :: npoints_poly,i3,j3
      real ,dimension(:)  ,ALLOCATABLE :: xp,yp
      double precision :: x,y



! Details dans: https://docs.google.com/document/d/1bI1DeIjmaf8DNYUME7K8KsVX_lh4bnu5_ChKXJHXXlM/edit
! Details dans: https://docs.google.com/presentation/d/17KQ5fzQmqYz7Vp7yzZ6er4elShwLni25GOGIgFpBS3g/edit#slide=id.g516dd41227_0_6

! Initialisation
!      if(iteration3d==0) then !-initial->
!      if(.not.allocated(zone4_mask)) then !-initial->
       if(txt_=='init') then !-initial-> !05-07-19

! Nombre de couches
        zone4_nlayer=3
!       zone4_nlayer=2
!       zone4_nlayer=1

! Inverse ("1 sur") de l'epaisseur des couches du bilan
       zone4_inv_dz=1./150.
! Cas particulier zone4_nlayer=1: zone4_inv_dz=1./hmax

!      zone4_stretch_dz=1.  ! epaisseur constante = 1/zone4_inv_dz
       zone4_stretch_dz=0.7065 !0.414 ! epaisseur augmentant avec la profondeur (en surface 1/zone4_inv_dz)

! ICI CHACUN SE DEBROUILLE POUR LE LIRE SON FICHIER DE MASQUE AVEC LE
! BON FORMAT.....
         allocate(glob_mask(0:iglb+1,0:jglb+1)) ; glob_mask=0
         allocate(zone4_mask(0:imax+1,0:jmax+1)) ; zone4_mask=0

#ifdef in_out_flux
      open(unit=12,file='/tmpdir/culses/bassin2/GRAPHIQUES/simu2/mask_mldDCAf1000_budget2t.dat')
        do j=0,jglb+1
        read(12,*)(glob_mask(i,j),i=0,iglb+1)
        enddo
      close(12)

! Zone -1
        do j=0,jmax+1 ; do i=0,imax+1
        if(i+par%timax(1)<920) then
        zone4_mask(i,j)=(-1)*glob_mask(i+par%timax(1),j+par%tjmax(1))*mask_t(i,j,kmax)
        endif 
        enddo         ; enddo
! ATTENTION DE DEFINIR DE i=0 A i=iglb+1 et de j=0 A j=jglb+1
! sinon (par ex de i=1 A iglb) les flux aux OBC ne seront pas calculEs
! et le bilan ne sera pas équilibrer

! Zone  1 : 
        do j=0,jmax+1 ; do i=0,imax+1
            if(zone4_mask(i,j).ne.-1) then
             zone4_mask(i,j)=    1*mask_t(i,j,kmax)
            endif
         enddo            ; enddo
#endif

      open(unit=12,file='zone1.dat')
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))

        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=0,imax+1
        do j1=0,jmax+1
          if (mask_t(i1,j1,kmax)==1) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              zone4_mask(i1,j1)=-1*mask_t(i1,j1,kmax)
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)



          do j=0,jmax+1 ; do i=0,imax+1
          if(zone4_mask(i,j).ne.-1*mask_t(i,j,kmax)) then
           zone4_mask(i,j)=mask_t(i,j,kmax)   ! zone 1
          endif
         enddo         ; enddo

         k0=maxval(zone4_mask)
         call mpi_allreduce(k0,zone4_max,1,mpi_integer,mpi_max,par%comm2d ,ierr)
         allocate(zone4bioflux_glb  (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_glb=0.
         allocate(zone4bioflux_u    (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_u=0.
         allocate(zone4bioflux_v    (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_v=0.

         allocate(zone4biocumul_glb   (0:zone4_nlayer-1))       ; zone4biocumul_glb=0. !22-09-20 v289
         allocate(zone4biocumul_loc   (0:zone4_nlayer-1))       ; zone4biocumul_loc=0.
         allocate(zone4biomasst0      (0:zone4_nlayer-1,vbmax)) ; zone4biomasst0=0.
         allocate(zone4tendancebio_glb(0:zone4_nlayer-1,vbmax)) ; zone4tendancebio_glb=0.
         allocate(zone4botsurfbio_glb (0:zone4_nlayer-1,vbmax)) ; zone4botsurfbio_glb=0.

#ifdef in_out_flux
         allocate(zone4bioflux_glb_in (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_glb_in=0.
         allocate(zone4bioflux_u_in (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_u_in=0.
         allocate(zone4bioflux_v_in (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_v_in=0.

         allocate(zone4bioflux_glb_out (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_glb_out=0.
         allocate(zone4bioflux_u_out (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_u_out=0.
         allocate(zone4bioflux_v_out (0:zone4_max,0:zone4_nlayer-1,vbmax))   ; zone4bioflux_v_out=0.
#endif

         allocate(zone4waterflux_glb (0:zone4_max+1,0:zone4_nlayer-1))   ; zone4waterflux_glb=0.
         allocate(zone4waterflux_u (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_u=0.
         allocate(zone4waterflux_v (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_v=0.

#ifdef in_out_flux
         allocate(zone4waterflux_glb_in (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_glb_in=0.
         allocate(zone4waterflux_u_in (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_u_in=0.
         allocate(zone4waterflux_v_in (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_v_in=0.

         allocate(zone4waterflux_glb_out (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_glb_out=0.
         allocate(zone4waterflux_u_out (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_u_out=0.
         allocate(zone4waterflux_v_out (0:zone4_max,0:zone4_nlayer-1))   ; zone4waterflux_v_out=0.
#endif

         allocate(zone4watercumul_glb   (0:zone4_nlayer-1)) ; zone4watercumul_glb=0. !22-09-20 v289
         allocate(zone4watercumul_loc   (0:zone4_nlayer-1)) ; zone4watercumul_loc=0.
         allocate(zone4watermasst0      (0:zone4_nlayer-1)) ; zone4watermasst0=0.

! La zone4 interieure est -1, les autres >=0.
! Identifier points de flux u
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax   ; do i=1,imax+1  
          if( (zone4_mask(i,j)==-1.and.zone4_mask(i-1,j)/=-1) &
          .or.(zone4_mask(i,j)/=-1.and.zone4_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone4_u_max=k
         if(zone4_u_max>0) then !ppp>
          allocate(zone4_flux_u_node(zone4_u_max,3)) ; zone4_flux_u_node=0
          k=0 ! passage 2 pour renseigner tableau zone4_flux_u_node=0
          do j=1,jmax   ; do i=1,imax+1  
           if( (zone4_mask(i,j)==-1.and.zone4_mask(i-1,j)/=-1) &
           .or.(zone4_mask(i,j)/=-1.and.zone4_mask(i-1,j)==-1)) then !ooo>
           if(mask_i_u(i)*mask_j_u(j)==1) then !pmx>
            k=k+1
            zone4_flux_u_node(k,1)=i 
            zone4_flux_u_node(k,2)=j
            zone4_flux_u_node(k,3)=max(zone4_mask(i-1,j),zone4_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

! Identifier points de flux v
         k=0 ! passage 1 pour en determiner le nombre pour allocation A suivre
         do j=1,jmax+1   ; do i=1,imax
          if( (zone4_mask(i,j)==-1.and.zone4_mask(i,j-1)/=-1) &
          .or.(zone4_mask(i,j)/=-1.and.zone4_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1)k=k+1
          endif                                                   !ooo>
         enddo       ; enddo
         zone4_v_max=k
         if(zone4_v_max>0) then !ppp>
          allocate(zone4_flux_v_node(zone4_v_max,3)) ; zone4_flux_v_node=0
          k=0 ! passage 2 pour renseigner tableau zone4_flux_v_node=0
          do j=1,jmax+1   ; do i=1,imax
           if( (zone4_mask(i,j)==-1.and.zone4_mask(i,j-1)/=-1) &
           .or.(zone4_mask(i,j)/=-1.and.zone4_mask(i,j-1)==-1)) then !ooo>
           if(mask_i_v(i)*mask_j_v(j)==1) then !pmx>
            k=k+1
            zone4_flux_v_node(k,1)=i 
            zone4_flux_v_node(k,2)=j
            zone4_flux_v_node(k,3)=max(zone4_mask(i,j-1),zone4_mask(i,j))
           endif                               !pmx>
           endif                                                   !ooo>
          enddo       ; enddo
         endif                 !ppp>

      if(par%rank==0) then !#mpi-->>-->               
      print*,'dir',trim(tmpdirname)//'messages'
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine my_outputs_zone4bioflux'
      write(3,*)'zone4_nlayer=',zone4_nlayer
!     write(3,*)'epaisseur des classes',1./zone4_inv_dz
      write(3,*)'No de classe ,  Zsup(m)  ,  Zinf(m)'
      k=0
      write(3,*)k,' surface      ',-real(k+1)**(1./zone4_stretch_dz)/zone4_inv_dz
      do k=1,zone4_nlayer-2
      write(3,*)k,-real(k  )**(1./zone4_stretch_dz)/zone4_inv_dz  &
                 ,-real(k+1)**(1./zone4_stretch_dz)/zone4_inv_dz
      enddo
      k=zone4_nlayer-1
      write(3,*)k,-real(k)**(1./zone4_stretch_dz)/zone4_inv_dz,-hmax
      close(3)
      endif                !#mpi-->>-->               



       return
       endif                   !-initial->


       if(txt_=='i_bio') then !-flux-u-bio->
          do k1=1,zone4_u_max
            i=zone4_flux_u_node(k1,1)
            j=zone4_flux_u_node(k1,2)
           k2=zone4_flux_u_node(k1,3) ! numero de la zone4 adjacente
           k0=sign(1,zone4_mask(i,j)-zone4_mask(i-1,j))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1 pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
!                           abs(depth_u(i,j,k)*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
     abs(0.5*(depth_t(i-1,j,k)+depth_t(i,j,k))*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
                          ),0),zone4_nlayer-1)
            zone4bioflux_u(k2,k3,vb)= &
            zone4bioflux_u(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                                )   !ooo>
           if(vb.eq.1) then !vb=1>
            zone4waterflux_u(k2,k3)= &
            zone4waterflux_u(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                   ) !ooo>
           endif            !vb=1>

#ifdef in_out_flux
            if((anyv3d(i,j,k,1)+anyv3d(i,j,k,2)+anyv3d(i,j,k,3)+anyv3d(i,j,k,4)).gt.0.D0) then
            zone4bioflux_u_in(k2,k3,vb)= &
            zone4bioflux_u_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                       )   !ooo>
                if(vb.eq.1) &
            zone4waterflux_u_in(k2,k3)= &
            zone4waterflux_u_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                      )  !ooo> 
             else
            zone4bioflux_u_out(k2,k3,vb)= &
            zone4bioflux_u_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
                                           +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
                                                        )   !ooo>
                if(vb.eq.1) &
            zone4waterflux_u_out(k2,k3)= &
            zone4waterflux_u_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i  ,j,k,3) &
                                           +anyv3d(i  ,j,k,2) &
                                           +anyv3d(i  ,j,k,1) &
                                           +anyv3d(i  ,j,k,4) &
                                                       )  !ooo>

             endif !positif/negatif
#endif

           enddo !k
          enddo  ! k1
       return
       endif            !-flux-u-bio->

       if(txt_=='j_bio') then !-flux-v-bio->
          do k1=1,zone4_v_max
            i=zone4_flux_v_node(k1,1)
            j=zone4_flux_v_node(k1,2)
           k2=zone4_flux_v_node(k1,3) ! numero de la zone4 adjacente
           k0=sign(1,zone4_mask(i,j)-zone4_mask(i,j-1))
           do k=1,kmax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1 pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
     abs(0.5*(depth_t(i,j-1,k)+depth_t(i,j,k))*zone4_inv_dz)**zone4_stretch_dz & !22-09-20   
                       ),0),zone4_nlayer-1)
            zone4bioflux_v(k2,k3,vb)= &
            zone4bioflux_v(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                    )   !ooo>
              if(vb.eq.1) then !AAA>
            zone4waterflux_v(k2,k3)= &
            zone4waterflux_v(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                   ) !ooo>

               endif           !AAA>

#ifdef in_out_flux
            if((anyv3d(i,j,k,5)+anyv3d(i,j,k,6)+anyv3d(i,j,k,7)+anyv3d(i,j,k,8)).gt.0.D0) then
            zone4bioflux_v_in(k2,k3,vb)= &
            zone4bioflux_v_in(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                       )   !ooo>
              if(vb.eq.1) &
            zone4waterflux_v_in(k2,k3)= &
            zone4waterflux_v_in(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                      ) !ooo>

             else
            zone4bioflux_v_out(k2,k3,vb)= &
            zone4bioflux_v_out(k2,k3,vb)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
                                           +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
                                                        )       !ooo>
              if(vb.eq.1) &
            zone4waterflux_v_out(k2,k3)= &
            zone4waterflux_v_out(k2,k3)-dti_fwsubio*k0*( & !ooo>  
                                            anyv3d(i,j  ,k,7) &
                                           +anyv3d(i,j  ,k,6) &
                                           +anyv3d(i,j  ,k,5) &
                                           +anyv3d(i,j  ,k,8) &
                                                       ) !ooo>
             endif
#endif
    
           enddo
          enddo ! k1
       return
       endif            !-flux-v-bio->

       if(txt_=='botsurf') then !-flux-w->
       do vb=1,vbmax
! integrale +dti_fw*fluxbio_w(i,j,vb,2)/dz_t(i,j,k,1)
! Note: attention operation repetee A chaque iteration principale
         zone4biocumul_loc=0.
! SURFACE:
         k=kmax
         do j=1,jmax ; do i=1,imax 
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1 pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
                          ),0),zone4_nlayer-1)

         zone4biocumul_loc(k3)= &
         zone4biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone4_mask(i,j)) & !=1 si zone4_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*fluxbio_w(i,j,vb,2) &
                                                          /dz_t(i,j,k,1)
         enddo       ; enddo
! FOND: (note: si plusieurs layer alors flux de surface et de fond sont distinguEs
         do j=1,jmax ; do i=1,imax 
         k=kmin_w(i,j)
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1 pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
                          ),0),zone4_nlayer-1)

         zone4biocumul_loc(k3)= &
         zone4biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone4_mask(i,j)) & !=1 si zone4_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*(fluxbio_w(i,j,vb,1)+wsed(1,vb)*bio_t(i,j,1,vb)) &
                                                           /dz_t(i,j,k,1)
         enddo       ; enddo
         call mpi_allreduce(zone4biocumul_loc,zone4biocumul_glb,zone4_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone4_nlayer-1
         zone4botsurfbio_glb(k3,vb)= &
         zone4botsurfbio_glb(k3,vb)  &
           +zone4biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb

       do j=1,jmax ; do i=1,imax
           zone4waterflux_w=zone4waterflux_w  &
                         -max(-zone4_mask(i,j),0) & ! selectionne zone4_mask=-1
           *dti_fw*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) &
           *(omega_w(i,j,kmax+1,1)-omega_w(i,j,1,1))
       enddo       ; enddo

       return
       endif                   !-flux-w->

! tendancebio
       if(txt_=='tendancebio') then !--tendancebio-->
       do vb=1,vbmax

! integrale +dti_fw*tendancebio_t(i,j,k,vb)
! Note: attention operation repetee A chaque iteration principale
         zone4biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1 pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
                          ),0),zone4_nlayer-1)

         zone4biocumul_loc(k3)= &
         zone4biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone4_mask(i,j)) & !=1 si zone4_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dti_fw*tendancebio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone4biocumul_loc,zone4biocumul_glb,zone4_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)

         do k3=0,zone4_nlayer-1
         zone4tendancebio_glb(k3,vb)= &
         zone4tendancebio_glb(k3,vb)  &
           +zone4biocumul_glb(k3)
         enddo ! k3

       enddo ! fin de boucle vb
       return
       endif                        !--tendancebio-->

! mpi sum
       if(txt_=='mpi') then !--mpi-->
       do vb=1,vbmax


        do k3=0,zone4_nlayer-1

         do k2=0,zone4_max

          x1=0. ; x3=0. ; x5=0.
          if(zone4_u_max>0) then !>>>
           x3=   zone4bioflux_u(k2,k3,vb)
          endif                 !>>>
          if(zone4_v_max>0) then !>>>
           x3=x3+ zone4bioflux_v(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone4bioflux_glb(k2,k3,vb)=zone4bioflux_glb(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3

#ifdef in_out_flux
        do k3=0,zone4_nlayer-1

         do k2=0,zone4_max

          x1=0. ; x3=0. ; x5=0.
          if(zone4_u_max>0) then !>>>
           x3=   zone4bioflux_u_in(k2,k3,vb)
          endif                 !>>>
          if(zone4_v_max>0) then !>>>
           x3=x3+ zone4bioflux_v_in(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone4bioflux_glb_in(k2,k3,vb)=zone4bioflux_glb_in(k2,k3,vb)+x4

         enddo ! k2

        enddo ! k3


        do k3=0,zone4_nlayer-1

         do k2=0,zone4_max

          x1=0. ; x3=0. ; x5=0.
          if(zone4_u_max>0) then !>>>
           x3=   zone4bioflux_u_out(k2,k3,vb)
          endif                 !>>>
          if(zone4_v_max>0) then !>>>
           x3=x3+ zone4bioflux_v_out(k2,k3,vb)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone4bioflux_glb_out(k2,k3,vb)=zone4bioflux_glb_out(k2,k3,vb)+x4
                            
         enddo ! k2       

        enddo ! k3
#endif

! water
        if(vb.eq.1) then !vb=1>

        do k3=0,zone4_nlayer-1
         do k2=0,zone4_max

          x1=0. ; x3=0. ; x5=0.
          if(zone4_u_max>0) then !>>>
           x3=   zone4waterflux_u(k2,k3)
          endif                 !>>>
          if(zone4_v_max>0) then !>>>
           x3=x3+ zone4waterflux_v(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone4waterflux_glb(k2,k3)=zone4waterflux_glb(k2,k3)+x4

         enddo ! k2
        enddo ! k3


#ifdef in_out_flux
! water in
        do k3=0,zone4_nlayer-1
         do k2=0,zone4_max

          x1=0. ; x3=0. ; x5=0.
          if(zone4_u_max>0) then !>>>
           x3=   zone4waterflux_u_in(k2,k3)
          endif                 !>>>
          if(zone4_v_max>0) then !>>>
           x3=x3+ zone4waterflux_v_in(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone4waterflux_glb_in(k2,k3)=zone4waterflux_glb_in(k2,k3)+x4

         enddo ! k2
        enddo ! k3

! water out
        do k3=0,zone4_nlayer-1
         do k2=0,zone4_max

          x1=0. ; x3=0. ; x5=0.
          if(zone4_u_max>0) then !>>>
           x3=   zone4waterflux_u_out(k2,k3)
          endif                 !>>>
          if(zone4_v_max>0) then !>>>
           x3=x3+ zone4waterflux_v_out(k2,k3)
          endif                 !>>>
          call mpi_allreduce(x3,x4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
          zone4waterflux_glb_out(k2,k3)=zone4waterflux_glb_out(k2,k3)+x4

         enddo ! k2
        enddo ! k3
#endif

        endif !vb=1>     end water


!        k3=0           ! k3 est la classe de profondeur de la couche de surface
!        k2=zone4_max+1 ! echange vertical A travers la surface (justifiE si zone4s assechees)
!        x1=zone4bioflux_w
!        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!        zone4bioflux_glb(k2,k3,vb)=zone4bioflux_glb(k2,k3,vb)+x2

! integrale bio
         zone4biocumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1 pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
                          ),0),zone4_nlayer-1)

         zone4biocumul_loc(k3)= &
         zone4biocumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone4_mask(i,j)) & !=1 si zone4_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))*bio_t(i,j,k,vb)

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone4biocumul_loc,zone4biocumul_glb,zone4_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone4_nlayer-1
            zone4biomasst0(k,vb)=-zone4biocumul_glb(k)         &
                              +zone4tendancebio_glb(k,vb)      &
                               +zone4botsurfbio_glb(k,vb)      &
                            +sum(zone4bioflux_glb(:,k,vb))
           enddo
         endif                   !m°v°m> !22-09-20

! Integrale water
         zone4watercumul_loc=0.
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
! k3 est la classe de profondeur  k3=0 en surface, k3=zone4_nlayer-1
! pour z=zmax_zone4_nlayer.
            k3=min(max(int(   &
                            abs(depth_t(i,j,k)*zone4_inv_dz)**zone4_stretch_dz & !22-09-20
                          ),0),zone4_nlayer-1)

         zone4watercumul_loc(k3)= &
         zone4watercumul_loc(k3)  &
                   +mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
                                              *max(0,-zone4_mask(i,j)) & !=1 si zone4_mask==-1 , 0 sinon
            *0.5*(dz_t(i,j,k,2)+dz_t(i,j,k,1))

         enddo       ; enddo       ; enddo
         call mpi_allreduce(zone4watercumul_loc,zone4watercumul_glb,zone4_nlayer,mpi_double_precision,mpi_sum,par%comm2d,ierr)
         if(iteration3d==0) then !m°v°m> !22-09-20
           do k=0,zone4_nlayer-1
            zone4watermasst0(k)=-zone4watercumul_glb(k)         &
                            +sum(zone4waterflux_glb(:,k))
           enddo
         endif                   !m°v°m> !22-09-20

!         if(par%rank==0) then !-rank0->
!           write(texte30,'(a)')'tmp/obcwaterflux_zone4'
!           open(unit=3,file=texte30,position='append')
!! Notes: - zone4watermasst0 sert A transformer sum2 en une anomalie
!!        - colonnes 2 et 3 doivent etre identiques
!                      write(3,'(100(1x,e14.7))')                     &
!                      real(elapsedtime_now/86400.)                   & ! colonne 1
!                     ,sum(zone4watercumul_glb)+sum(zone4watermasst0) & ! colonne 2
!                     ,sum(zone4waterflux_glb)                        & ! colonne 3
!                     ,zone4waterflux_glb(0:zone4_max+1,0:zone4_nlayer-1) &
!                     ,zone4watercumul_glb(0:zone4_nlayer-1)+zone4watermasst0(0:zone4_nlayer-1)
!           close(3)
!         endif                !-rank0->


         if(par%rank==0) then !-rank0->
           sum0=0. ; sum1=0.
           do k3=0,zone4_nlayer-1
             sum0=sum0+zone4biocumul_glb(k3)+zone4biomasst0(k3,vb)
             do k2=0,zone4_max
               sum1=sum1+zone4bioflux_glb(k2,k3,vb)
!              write(666,*)iteration3d,k2,k3,vb,zone4bioflux_glb(k2,k3,vb)
             enddo
             sum1=sum1+zone4tendancebio_glb(k3,vb) &
                       +zone4botsurfbio_glb(k3,vb)  
           enddo
           
!......................
! Ecrire un fichier pour verifier l'equilibrage
           write(texte60,'(a,i0)')'tmp/zone4_bilan_total_bio',vb
           open(unit=3,file=texte60,position='append')
! Notes: - zone4biomasst0 sert A transformer sum2 en une anomalie
!        - colonnes 2 et 3 doivent etre identiques
                      write(3,'(100(1x,e14.7))')    &
                      real(elapsedtime_now/86400.)  & ! colonne 1
                     ,sum0                          & ! colonne 2
                     ,sum1                            ! colonne 3
           close(3)
! Note: la colonne 2 et la colonne 3 doivent etre egales
! colonne 2: variation du contenu du traceur cumulee sur toutes les
!            couches si nlayer>1
! colonne 3: contribution cumulee de tous les termes du bilan (flux
!            lateraux, surface, tendancebio etc....)

!......................
! Ecrire un fichier par couche verticale (numero 0 pour couche de surface)
           do k=0,zone4_nlayer-1
           write(texte60,'(a,i0,a,i0)')'tmp/zone4_bilan_detail_bio',vb &
                                      ,'_layer',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & ! colonne 1
                     ,zone4biocumul_glb(k)+zone4biomasst0(k,vb) & ! colonne 2
                     ,zone4tendancebio_glb(k,vb)                & ! colonne 3
                      ,zone4botsurfbio_glb(k,vb)                & ! colonne 4
                     ,zone4bioflux_glb(0:zone4_max,k,vb)          ! colonne 5=zone0, colonne 6=zone4, etc...
           close(3)

! water
        if(vb.eq.1) then
! Ecrire fichier detail
           write(texte60,'(a,i0)')'tmp/zone4_bilan_water_layer',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(120(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone4watercumul_glb(k)+zone4watermasst0(k) & !colonne 2
                     ,zone4waterflux_glb(0:zone4_max,k)        
           close(3)
         endif ! vb water



#ifdef in_out_flux
! in/out flux
           write(texte60,'(a,i0,a,i0)')'tmp/zone4_bilan_detail_bio',vb &
                                      ,'_layer_inout',k
           open(unit=3,file=texte60,position='append')
                      write(3,'(100(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone4bioflux_glb_in(0:zone4_max,k,vb)     & 
                     ,zone4bioflux_glb_out(0:zone4_max,k,vb)            
           close(3)
! water
        if(vb.eq.1) then
! Ecrire fichier detail
           write(texte60,'(a,i0)')'tmp/zone4_bilan_detail_water_layer_inout',k 
           open(unit=3,file=texte60,position='append')
                      write(3,'(120(1x,e14.7))')                &
                      real(elapsedtime_now/86400.)              & !colonne 1
                     ,zone4watercumul_glb(k)+zone4watermasst0(k) & !colonne 2
                     ,zone4waterflux_glb(0:zone4_max,k)        &
                     ,zone4waterflux_glb_in(0:zone4_max,k)     &
                     ,zone4waterflux_glb_out(0:zone4_max,k)    
           close(3)
         endif ! vb water
#endif

! Note:
! Si nlayer=1 alors la colonne 2 doit etre egale a la somme des colonnes
! suivantes (c.a.d. de la colonne 3 a la derniere colonne).
! Si nlayer>1 alors l'ecart entre la colonne 2 et la somme des autres
! colonnes s'explique par le flux vertical (advectif et diffusif) entre
! les couches (c'est donc comme cela que se deduit le flux vertical
! entre les couches).
! colonne 2: variation du contenu du traceur dans la couche concernee
! colonne 3: contribution de tendancebio_t
! colonne 4: contribution des flux de fond et de surface. Si nlayer=1
!            les 2 sont confondus. Si nlayer>1 alors la colonne 4 ne
!            contient que le flux de surface dans la couche 0 et que le
!            flux de fond dans la couche nlayer-1.
! colonne 5: contribution du flux lateral entre zone -1 et zone 0
!            la zone 0 est le masque continental donc il s'agit des
!            rivieres
! colonne 6: contribution du flux lateral entre zone -1 et zone 1
! colonne 7: etc...
           enddo ! fin de boucle sur k
         endif                !-rank0->

       enddo ! fin de boucle sur vb

! reset des tableaux locaux pour un nouveau cycle de cumul
       zone4bioflux_u=0. 
       zone4bioflux_v=0. 
       zone4bioflux_w=0.
       zone4waterflux_u=0.
       zone4waterflux_v=0.
       zone4waterflux_w=0.
#ifdef in_out_flux
       zone4bioflux_u_in=0.    
       zone4bioflux_v_in=0.  
       zone4bioflux_u_out=0.   
       zone4bioflux_v_out=0.  
       zone4waterflux_u_in=0.  
       zone4waterflux_v_in=0.  
       zone4waterflux_u_out=0. 
       zone4waterflux_v_out=0.   
       zone4waterflux_u_in=0.
       zone4waterflux_v_in=0.
       zone4waterflux_u_out=0.
       zone4waterflux_v_out=0.
#endif
       return
       endif               !--mpi-->

       stop 'Err undefined txt_ in my_outputs_zone4bioflux'

      end subroutine my_outputs_zone4bioflux
#endif

      end module module_my_outputs
