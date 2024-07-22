










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



      integer :: ncid_,time_counter_=0 &
                ,i_,j_,flag_header_=0  &
                ,year_
      real*4,dimension(:,:),allocatable :: &
       angle_wave_beam_w
      double precision,dimension(:),allocatable :: &
       seconds_since_1jan

      double precision, dimension(:),allocatable :: zprofile_depth !,zprofile_tem,zprofile_sal

      character*14 :: date_

contains

      subroutine my_outputs_driver
      implicit none

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








! Mar Menor
      call my_outputs_point_vs_time('pointA' ,'lonlat',-0.783317d0,37.792949d0)
      call my_outputs_point_vs_time('pointB' ,'lonlat',-0.783191d0,37.698050d0)
      call my_outputs_point_vs_time('pointC' ,'lonlat',-0.752359d0,37.667696d0)
      call my_outputs_point_vs_time('merMed1' ,'lonlat',-0.659274d0,37.72370d0)
      call my_outputs_point_vs_time('merMed2' ,'lonlat',-0.466717d0,37.73270d0)




! TEST:
!      call latlontoij(lon_t(10,10),lat_t(10,10),'loc') ! returns deci, decj
!      write(10+par%rank,*)deci,decj

!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!      stop 'cocote'
!#endif
        endif                          !sssssss>

        call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
        if(k0/=0)stop 'flag_stop=1 in my_outputs_driver' !26-03-19


      end subroutine my_outputs_driver

!.............................................................
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
!.............................................................
!.............................................................
!.............................................................
!.............................................................
!.............................................................

!#ifdef bidon
      subroutine my_outputs_globtke_vs_time
      implicit none
! Energie cinetique globale en fonction du temps

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

      call mpi_allreduce(sum1                            &
                        ,sum1glb,1,mpi_double_precision, & !05-10-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2                            & !26-04-11
                        ,sum2glb,1,mpi_double_precision, &
                         mpi_sum,par%comm2d,ierr)

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


      status=nf_close(ncid_)

      end subroutine my_outputs_point_vs_time_2d

!............................................................................
!............................................................................
!............................................................................
!............................................................................
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

!............................................................................

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
!.........................................................................
!.....................................................................
!.....................................................................


!.....................................................................


!.....................................................................


!.....................................................................

      end module module_my_outputs
