      module module_my_outputs
!...............................................................................
! S model
! release S26 - last update: 30-11-15
! Version date      Description des modifications
!         30-07-15: ecrire la date sous forme calendaire
!         31-07-15: my_outputs_read_point_vs_time permet de relire les fichiers
!                   les stations julio et cie... 
!         19-10-15  debug ecriture netcdf
!         27-11-15  ecriture pour cas 1DV
!         30-11-15  ajout variables de la couche limite atmospherique de S26
!         22-02-16  ajout variables bio (Alex)
!...............................................................................
      use module_principal
      use module_parallele
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

!     if(initial_main_status==1)call my_outputs_profiles !16-03-13
!     call my_outputs_cumulvel(1)
!     call my_outputs_cumultime
!     if(initial_main_status==1)call my_outputs_globtke_vs_time

      if(initial_main_status==1)then !sssssss>

      
! Le tableau seconds_since_1jan va servir a construire une date en annee decimale
         if(.not.allocated(seconds_since_1jan)) then !>>>>>
           allocate(seconds_since_1jan(datesim(1,1):datesim(1,2)+1)) ! 1ere annee -> derniere annee+1
           do year_=datesim(1,1),datesim(1,2)+1
              call datetokount(year_,0,0,0,0,0) 
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

!        call my_outputs_point_vs_time('julio'  ,'lonlat',5.255d0,43.135d0)
!        call my_outputs_point_vs_time('solemio','lonlat',5.289d0,43.237d0)
!        call my_outputs_point_vs_time('lion'   ,'lonlat',4.64d0 ,42.06d0)

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!      stop 'cocote'
!#endif
        endif                          !sssssss>

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

! Eventuellement si on ne veut pas archiver à chaque iteration pour 
! limiter la taille du fichier, sortir si modulo(iteration3d,N)/=0
      if(modulo(iteration3d,20)/=0)return

! calculer les indices du point selon convention:

      if(trim(locationconv_)=='lonlat') then  !-----> ! convention longitude latitude
       lon_=posx_ ; lat_=posy_ ; ii_=1
       call latlontoij(lon_*deg2rad,lat_*deg2rad,'loc') ! returns deci, decj
!      i_=nint(deci) ; j_=nint(decj) ! strategie un point (le + proche)
       i_= int(deci) ; j_= int(decj) ! strategie quatre points encadrants
      else                                    !----->
       if(trim(locationconv_)=='ijglob') then !/////> ! convention longitude latitude
        ii_=0
!       i_=nint(posx_)-par%timax(1) ; j_=nint(posy_)-par%tjmax(1) ! strategie un point (le plus proche)
        i_= int(posx_)-par%timax(1) ; j_= int(posy_)-par%tjmax(1) ! strategie quatre points encadrants
       else                                   !/////>
        stop 'STOP module_my_output err on locationconv_'
       endif                                  !/////>
      endif                                   !----->


! Si les indices (locaux) sont hors proc sortir:
      if(i_<2.or.i_>imax-1.or.j_<2.or.j_>jmax-1)return

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
                                                    ! returns i1,i2,...=y,m,d,h,m,s
! Time en secondes (egalement variable dimension)
      write(texte80(2),'(a13,i4,5(a1,i2))')                    & !units
      'seconds from ',i1,'-',i2,'-',i3,' ',i4,':',i5,':',i6

      status=nf_def_var(ncid_,'time',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time'
      status=nf_put_att_text(ncid_,var_id,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)stop ' erreur nf_put_att_text'

! Time en jours
      write(texte80(2),'(a10,i4,5(a1,i2))')                    & !units
      'days from ',i1,'-',i2,'-',i3,' ',i4,':',i5,':',i6

      status=nf_def_var(ncid_,'time_days',nf_double,1,vardim(4),var_id)
      if(status/=0)stop ' erreur nf_def_var time_days'
      status=nf_put_att_text(ncid_,var_id,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)stop ' erreur nf_put_att_text'

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

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'temobc_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var tem_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',30 &
       ,'Initial_and_border_temperature')
      status=nf_put_att_text(ncid_,var_id,'units',1,'C')

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'sal_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var sal_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',21 &
       ,'salinity_of_sea_water')
      status=nf_put_att_text(ncid_,var_id,'units',3,'PSU')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'salobc_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var tem_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',27 &
       ,'Initial_and_border_salinity')
      status=nf_put_att_text(ncid_,var_id,'units',1,'C')

      vardim(3)=k_t_dim !15-02-16
      status=nf_def_var(ncid_,'rhp_t',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var rhp_t'
      status=nf_put_att_text(ncid_,var_id,'long_name',17 &
       ,'Potential_density')
      status=nf_put_att_text(ncid_,var_id,'units',7,'kg/m**3')


! AJOUT VARIABLES BIO (Alex 22/02/2016)
      if(imodelbio==1) then !bbbbbbbbbbbbbb>

!!! Nutriments
! Nitrates
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'INITRATE',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var INITRATE'
      status=nf_put_att_text(ncid_,var_id,'long_name',21 &
       ,'concentration_nitrate')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! Ammonium
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'IAMMONIUM',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var IAMMONIUM'
      status=nf_put_att_text(ncid_,var_id,'long_name',22 &
       ,'concentration_ammonium')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! Phosphates
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'IPHOSPHATE',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var IPHOSPHATE'
      status=nf_put_att_text(ncid_,var_id,'long_name',23 &
       ,'concentration_phosphate')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! Silice
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'ISILICE',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var ISILICE'
      status=nf_put_att_text(ncid_,var_id,'long_name',20 &
       ,'concentration_silice')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

!!! Phytoplancton
!! Chlorophylle
! Picophytoplancton en chlorophylle
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'ISYNECHL',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var ISYNECHL'
      status=nf_put_att_text(ncid_,var_id,'long_name',36 &
       ,'concentration_chlorophylle_picophyto')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')

! Nanophytoplancton en chlorophylle
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'INANOCHL',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var INANOCHL'
      status=nf_put_att_text(ncid_,var_id,'long_name',36 &
       ,'concentration_chlorophylle_nanophyto')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')

! Microphytoplancton en chlorophylle
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'IDIACHL',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var IDIACHL'
      status=nf_put_att_text(ncid_,var_id,'long_name',37 &
       ,'concentration_chlorophylle_microphyto')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')

! Chlorophylle totale
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'ITOTALCHL',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var ITOTALCHL'
      status=nf_put_att_text(ncid_,var_id,'long_name',33 &
       ,'concentration_chlorophylle_totale')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mg.m-3')

!! Carbone
! Picophytoplancton en C
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'ISYNEC',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var ISYNEC'
      status=nf_put_att_text(ncid_,var_id,'long_name',25 &
       ,'concentration_C_picophyto')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! Nanophytoplancton en C
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'INANOC',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var INANOC'
      status=nf_put_att_text(ncid_,var_id,'long_name',25 &
       ,'concentration_C_nanophyto')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')

! Microphytoplancton en C
      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'IDIAC',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var IDIAC'
      status=nf_put_att_text(ncid_,var_id,'long_name',26 &
       ,'concentration_C_microphyto')
      status=nf_put_att_text(ncid_,var_id,'units',8,'mmol.m-3')



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
      status=nf_def_var(ncid_,'V_OGCM_eastward',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var VOGCM_eastward'
      status=nf_put_att_text(ncid_,var_id,'long_name',17 &
       ,'OGCM_eastward_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_t_dim 
      status=nf_def_var(ncid_,'V_OGCM_northward',nf_real,4,vardim(1:4),var_id)
      if(status/=0)stop ' erreur nf_def_var VOGCM_northward'
      status=nf_put_att_text(ncid_,var_id,'long_name',18 &
       ,'OGCM_northward_velocity')
      status=nf_put_att_text(ncid_,var_id,'units',3,'m/s')
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,filvalr4_)

      vardim(3)=k_w_dim 
      status=nf_def_var(ncid_,'N2',nf_real,4,vardim(1:4),var_id)
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
      status=nf_put_att_text(ncid_,var_id,'units',1,'K')

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
      if(status/=0)stop ' erreur nf_def_var sshf'
      status=nf_put_att_text(ncid_,var_id,'long_name',13 &
      ,'precipitation')
      status=nf_put_att_text(ncid_,var_id,'units',6,'m/s')

! Global attributs:
      k10=len(trim(name_))
      status=nf_put_att_text(ncid_,nf_global,'Point',k10,trim(name_))
      k10=len(trim(model_name))
      status=nf_put_att_text(ncid_,nf_global,'Model',k10,trim(model_name))

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
          x1=second_now             &
            +minute_now*100         &
            +hour_now  *10000       &
            +day_now   *1000000     &
            +month_now *100000000   &
            +year_now  *10000000000
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
      if(status/=0)stop' erreur nf_inq_varid depth_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real tem_t'

!...................
! temperature de mercator
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))                      &
         +(     timeweightobc(trc_id) *temobc_t(i,j,k,2)               &
           +(1.-timeweightobc(trc_id))*temobc_t(i,j,k,0))*mask_t(i,j,k)
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'temobc_t',var_id)
      if(status/=0)stop' erreur nf_inq_varid temobc_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

!...................
! Ecrire le profil vertical de salinite:
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=sal_t(i,j,k,now)*mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'sal_t',var_id)
      if(status/=0)stop' erreur nf_inq_varid sal_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real sal_t'

!...................
! salinite de mercator
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=filvalr4_*(1-mask_t(i,j,k))                      &
         +(     timeweightobc(trc_id) *salobc_t(i,j,k,2)               &
           +(1.-timeweightobc(trc_id))*salobc_t(i,j,k,0))*mask_t(i,j,k)
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'salobc_t',var_id)
      if(status/=0)stop' erreur nf_inq_varid salobc_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

!...................
! Ecrire le profil de densite potentielle
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_
      anyvar3d(i,j,k)=(rhp_t(i,j,k)+rho)*mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'rhp_t',var_id)
      if(status/=0)stop' erreur nf_inq_varid rhp_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real rhp_t'


!...................
      write(6,*)'imodelbio=',imodelbio
      if(imodelbio==1) then !bbbbbbbbbbbbbb>

! Nitrates
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,INITRATE)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'INITRATE',var_id)
      if(status/=0)stop' erreur nf_inq_varid INITRATE'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Ammonium
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,IAMMONIUM)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'IAMMONIUM',var_id)
      if(status/=0)stop' erreur nf_inq_varid IAMMONIUM'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Phosphates
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,IPHOSPHATE)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'IPHOSPHATE',var_id)
      if(status/=0)stop' erreur nf_inq_varid IPHOSPHATE'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Silice
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,ISILICE)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'ISILICE',var_id)
      if(status/=0)stop' erreur nf_inq_varid ISILICE'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Syne Chl
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,ISYNECHL)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'ISYNECHL',var_id)
      if(status/=0)stop' erreur nf_inq_varid ISYNECHL'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Nano Chl
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,INANOCHL)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'INANOCHL',var_id)
      if(status/=0)stop' erreur nf_inq_varid INANOCHL'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Micro Chl
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,IDIACHL)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'IDIACHL',var_id)
      if(status/=0)stop' erreur nf_inq_varid IDIACHL'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Total Chl
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=(bio_t(i,j,k,ISYNECHL)+bio_t(i,j,k,INANOCHL)+bio_t(i,j,k,IDIACHL))*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'ITOTALCHL',var_id)
      if(status/=0)stop' erreur nf_inq_varid ITOTALCHL'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Syne C
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,ISYNEC)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'ISYNEC',var_id)
      if(status/=0)stop' erreur nf_inq_varid ISYNEC'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Nano C
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,INANOC)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'INANOC',var_id)
      if(status/=0)stop' erreur nf_inq_varid INANOC'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                  ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))

! Micro C
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,k)=bio_t(i,j,k,IDIAC)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'IDIAC',var_id)
      if(status/=0)stop' erreur nf_inq_varid IDIAC'
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
      if(status/=0)stop' erreur nf_inq_varid Veastward'
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
      if(status/=0)stop' erreur nf_inq_varid Vnorthward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real Vnorthward'

!...................
! Ecrire le profil de courant OGCM OE
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      x2=timeweightobc(vel_id) ; x0=1.-x2
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
          anyvar3d(i,j,k)=0.5*( &
       x0*( ( velobc_u(i  ,j  ,k,0)                       &
             +velobc_u(i+1,j  ,k,0))*gridrotcos_t(i,j)    &
           +( velobc_v(i  ,j  ,k,0)                       &
             +velobc_v(i  ,j+1,k,0))*gridrotsin_t(i,j))   &
      +x2*( ( velobc_u(i  ,j  ,k,2)                       &
             +velobc_u(i+1,j  ,k,2))*gridrotcos_t(i,j)    &
           +( velobc_v(i  ,j  ,k,2)                       &
             +velobc_v(i  ,j+1,k,2))*gridrotsin_t(i,j))   &
                               )*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'V_OGCM_eastward',var_id)
      if(status/=0)stop' erreur nf_inq_varid V_OGCM_eastward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real V_OGCM_eastward'

!...................
! Ecrire le profil de courant OGCM SN
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      x2=timeweightobc(vel_id) ; x0=1.-x2
      do k=1,kmax ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
          anyvar3d(i,j,k)=0.5*( &
       x0*(-( velobc_u(i  ,j  ,k,0)                       &
             +velobc_u(i+1,j  ,k,0))*gridrotsin_t(i,j)    &
           +( velobc_v(i  ,j  ,k,0)                       &
             +velobc_v(i  ,j+1,k,0))*gridrotcos_t(i,j) )  &
      +x2*(-( velobc_u(i  ,j  ,k,2)                       &
             +velobc_u(i+1,j  ,k,2))*gridrotsin_t(i,j)    &
           +( velobc_v(i  ,j  ,k,2)                       &
             +velobc_v(i  ,j+1,k,2))*gridrotcos_t(i,j) )  &
                              )*mask_t(i,j,k)+filvalr4_*(1.-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'V_OGCM_northward',var_id)
      if(status/=0)stop' erreur nf_inq_varid V_OGCM_northward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
                ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax))
      if(status/=0)stop 'Erreur nf_put_vara_real V_OGCM_northward'

!...................

!...................
! Ecrire le profil vertical de N**2
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
      anyvar3d(i,j,1)=0. ; anyvar3d(i,j,kmax+1)=0.
      do k=2,kmax 
       anyvar3d(i,j,k)=mask_t(i,j,k)*grav/rho*(rhp_t(i,j,k)  -rhp_t(i,j,k-1)) &
                                           /(depth_w(i,j,k)-depth_w(i,j,k-1)) &
         +filvalr4_*(1-mask_t(i,j,k))
      enddo
      enddo ; enddo
      status=nf_inq_varid(ncid_,'N2',var_id)
      if(status/=0)stop' erreur nf_inq_varid N2'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real N2'

!...................
! Ecrire le profil vertical de km_w
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,k)=km_w(i,j,k)*mask_t(i,j,k) &
                     +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'km_w',var_id)
      if(status/=0)stop' erreur nf_inq_varid km_w'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real km_w'

!...................
! Ecrire le profil vertical de kh_w
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,k)=(kh_w(i,j,k)+kmol_h)*mask_t(i,j,k) &
                              +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'kh_w',var_id)
      if(status/=0)stop' erreur nf_inq_varid kh_w'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real kh_w'

!...................
! Ecrire le profil vertical de tken_w
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,k)=tken_w(i,j,k)*mask_t(i,j,k) &
                       +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      status=nf_inq_varid(ncid_,'tken_w',var_id)
      if(status/=0)stop' erreur nf_inq_varid tken_t'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real tken_w'

!...................
! Ecrire le profil vertical de epsn_w
      varstart(4)=time_counter_  ; varcount(4)=1
      varstart(3)=1              ; varcount(3)=kmax+1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      if(iturbulence==1) then  !ooooo>
      do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
        anyvar3d(i,j,k)=epsn_w(i,j,k)*mask_t(i,j,k) &
                        +filvalr4_*(1-mask_t(i,j,k))
      enddo ; enddo ; enddo
      else                     !ooooo>
       if(iturbulence==0) then !mmmmm>
       do k=1,kmax+1 ; do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
!epsilon dans le schema de Gaspar=ctke2*(tken_w(i,j,k)**1.5)/(tkle_w(i,j,k)+small2)
               anyvar3d(i,j,k)=        &
        ctke2*(tken_w(i,j,k)**1.5)     &
             /(tkle_w(i,j,k)+small2)   &
              *mask_t(i,j,k)+filvalr4_*(1-mask_t(i,j,k))
       enddo ; enddo ; enddo
       else                    !mmmmm>
        stop 'STOP Err my_output 296' 
       endif                   !mmmmm>
      endif                    !ooooo>
      status=nf_inq_varid(ncid_,'epsn_w',var_id)
      if(status/=0)stop' erreur nf_inq_varid epsn_w'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:4)         &
                                   ,varcount(1:4)         &
              ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1:kmax+1))
      if(status/=0)stop 'Erreur nf_put_vara_real epsn_w'

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
      if(status/=0)stop' erreur nf_inq_varid Weastward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real Weastward'

      status=nf_inq_varid(ncid_,'Wnorthward',var_id)
      if(status/=0)stop' erreur nf_inq_varid Wnorthward'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,2))
      if(status/=0)stop 'Erreur nf_put_vara_real Wnorthward'

!...................
      if(flag_abl==1) then !AAAAAAAAAAA>

! Ecrire la temperature potentielle de l'air a 2m:
      status=nf_inq_varid(ncid_,'teta2_t',var_id)
      if(status/=0)stop 'Erreur nf_inq_varid teta2'
      status=nf_put_vara_double(ncid_,var_id      &
                                   ,varstart(1:3) &
                                   ,varcount(1:3) &
                ,teta2_t(i_:i_+ii_,j_:j_+ii_,now))
      if(status/=0)stop 'Erreur nf_put_vara_real teta2'

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
      if(status/=0)stop' erreur nf_inq_varid wstress_u'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real wstress_u'

      status=nf_inq_varid(ncid_,'wstress_v',var_id)
      if(status/=0)stop' erreur nf_inq_varid wstress_v'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,2))
      if(status/=0)stop 'Erreur nf_put_vara_real wstress_v'

! ssr:
      varstart(3)=time_counter_  ; varcount(3)=1
      varstart(1:2)=1            ; varcount(1:2)=1+ii_
      do j=j_,j_+ii_ ; do i=i_,i_+ii_ 
       anyvar3d(i,j,1)= ssr_w(i,j,1) 
      enddo ; enddo
      status=nf_inq_varid(ncid_,'ssr',var_id)
      if(status/=0)stop' erreur nf_inq_varid ssr'
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
      if(status/=0)stop' erreur nf_inq_varid snsf'
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
      if(status/=0)stop' erreur nf_inq_varid slhf'
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
      if(status/=0)stop' erreur nf_inq_varid sshf'
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
      if(status/=0)stop' erreur nf_inq_varid precipitation'
      status=nf_put_vara_real(ncid_,var_id                &
                                   ,varstart(1:3)         &
                                   ,varcount(1:3)         &
                     ,anyvar3d(i_:i_+ii_,j_:j_+ii_,1))
      if(status/=0)stop 'Erreur nf_put_vara_real precipitation'

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

      end module module_my_outputs
