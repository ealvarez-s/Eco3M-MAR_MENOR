










      module module_depammo
!______________________________________________________________________
! S model
! release S26.1 - last update: 27-10-14
!______________________________________________________________________
!...............................................................................
      use module_principal
      use module_parallele
      use module_forcages
      implicit none
      include 'netcdf.inc'

      real*4,dimension(:,:),allocatable :: depammo_lonlat_r4 &
                      ,ij2depammo_i,ij2depammo_j

      integer flag_depammotime,flag_cumul,var_dims,var_type,lon_id,lat_id &
             ,tabdim(4)

      integer,dimension(:),allocatable :: varcode

      double precision lat_1_1_gridt,lat_1_1_gridu,lat_1_1_gridv &
                      ,lon_1_1_gridt,lon_1_1_gridu,lon_1_1_gridv &
                      ,lon_2_1_gridt,lat_1_2_gridt,scalefct

      real c_grid_shift

      character*10 :: depammovarname(4)



contains


!..............................................................................

      subroutine depammo_driver(case_)
      implicit none
      integer case_,loop_,t_,tstr_

      if(.not.allocated(varcode))allocate(varcode(ndepammo)) !08-03-13

      if(case_==1) then
       call depammo_initial(case_)
       call depammoflux_interp(0,case_)   ! read netcdf & grid interpolation
       call depammoflux_interp(2,case_)   ! read netcdf & grid interpolation
      endif            !- Initial phase ->

      if(case_==2)call depammoflux_interp(2,case_)   ! read netcdf & grid interpolation

      call depammoflux_flux_moveforward ! time linear interpolation
       

      end subroutine depammo_driver

!..............................................................................

      subroutine depammo_initial(case_)
      implicit none
      integer case_

! Initialiser les identifiants (type integer) des variables depammo
      call depammo_varidentifier

! A l'etat initial, faire la liste binaire des fichiers depammo
      call depammo_binary_file

! Reduire l'extraction des donnees a une zone englobant notre domaine
! numerique
! et definir les bouchetrous
      call depammoflux_extract_zone

      end subroutine depammo_initial

!..........................................................

      subroutine depammo_varidentifier
      implicit none


       do k=1,ndepammo

        k1=len_trim(depammo_file(k))
!        print*,'id',k1,k,depammo_file(k)
        if(depammo_file(k)(k1-2:k1)=='_dd')    dd_id=k
        if(depammo_file(k)(k1-3:k1)=='_wdr')   wdr_id=k
        if(depammo_file(k)(k1-3:k1)=='_wdw')   wdw_id=k

       enddo ! k

      end subroutine depammo_varidentifier

!..........................................................
      subroutine depammo_binary_file
      implicit none
      double precision time_,timeforecast_,time1_,time2_
      integer ncid_,flag_                                 &
             ,year_,month_,day_,hour_,minute_,second_     &
             ,loop1_,linemax_,loop2_,loop3_
      integer , dimension(:) , allocatable :: tabi4_
      double precision , dimension(:) , allocatable :: tabr8_

      cpu_seconds=MPI_Wtime ( )

! A l'etat initial, faire la liste binaire des fichiers meteo
! Reset zero (important car mpi_sum plus tard)
      depammofile_nextrec(:)=0
      depammofile_nextime(:)=0.

! Boucle sur les variables meteo: une variable = un fichier "liste" a traiter
      do loop1_=1,ndepammo

! Get meteovarname, the variable to be found in the meteo netcdf file
      call depammoflux_varname(loop1_) ! donne depammovarname le nom de la variable du fichier netcdf

! name of the record direct access binary file list:
        k2=index(depammo_file(loop1_),'/',back=.true.)  
        k3=index(depammo_file(loop1_),' ')
!      print*,k2,k3,depammo_file(loop1_)(k2+1:k3-1),'tmp/'//depammo_file(loop1_)(k2+1:k3-1)//'.binrec'
       depammobinreclist(loop1_)='tmp/'//depammo_file(loop1_)(k2+1:k3-1)//'.binrec'


        open(unit=3,file=depammo_file(loop1_))
! Combien y a t'il de lignes dans le fichier?
         linemax_=0
 1834    read(3,'(a)',end=1832)texte250 ; linemax_=linemax_+1 ; goto 1834
 1832    rewind(3)
         if(linemax_==0) then         !>>>>> 
          write(6,'(a,a,a)')'File ',trim(depammo_file(loop1_)),' is empty'
          stop ' stop depammo_binary_file_list 3284'
         endif                        !>>>>>

! Chaque proc va traiter une fraction du fichier. On commence par derouler les lignes !10-08-14
! qui ne concernent pas le proc
         do loop3_=1,int(real(par%rank  )/real(nbdom)*linemax_)
          read(3,*) ! lire ces lignes pour rien
         enddo

         write(texte60,'(a,i0)')'tmp/tmpdepammofile',par%rank
!        print*,'tmp/tmpdepammofile',texte60, &
!           int(real(par%rank)/real(nbdom)*linemax_)+1,    &
!           int(real(par%rank+1)/real(nbdom)*linemax_)
!        print*,'nbdom=',nbdom,'linemax=',linemax_ 

         open(unit=4,file=texte60,status='REPLACE')
          do loop3_=int(real(par%rank  )/real(nbdom)*linemax_)+1    &
                   ,int(real(par%rank+1)/real(nbdom)*linemax_)

! depammo file name given by the list:
          read(3,'(a)',end=332)texte80(1)
! open depammo file:
          status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
          if(status/=0) &
          stop 'error nf_open in subroutine depammo_binary_file'

! Get the scalefactor and the flag_cumul of the variable:
          call depammoflux_inquire_var(loop1_,ncid_) ! returns flag_cumul scalefct

! Get max_depammo_time_counter:
          max_depammo_time_counter=1 ; flag_depammotime=1
                       status=nf_inq_dimid(ncid_,'time_counter',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'TIME_COUNTER',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'time',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'TIME',dim_t_id)
          if(status/=0)flag_depammotime=0
          if(status==0)status=nf_inq_dimlen(ncid_,dim_t_id,max_depammo_time_counter)
! Get time units:
                      status=nf_inq_varid(ncid_,'time',var_id)
         if(status/=0)status=nf_inq_varid(ncid_,'time_counter',var_id)
         if(status/=0)stop 'erreur nf_inq_varid time depammo'
         texte60=''
         status=nf_get_att_text(ncid_,var_id,'units',texte60)
         if(status/=0)stop 'erreur nf_get_att_text depammo'
         k=index(texte60,' since ')
         read(texte60(k+7 :k+10),*)year_
         read(texte60(k+12:k+13),*)month_
         read(texte60(k+15:k+16),*)day_
         read(texte60(k+18:k+19),*)hour_
         read(texte60(k+21:k+22),*)minute_ ! note: les secondes ne sont pas toujours presentes dans fichiers ecmwf
!         print*,'verif lecture date',texte60,year_,month_,day_,hour_,minute_
!         stop'dans depammo binary file'
! Elapsed time (seconds) since the reference date of the units:
         call datetokount(year_,month_,day_,hour_,minute_,0) ! returns elapsedtime_out

! Get time of each field contained in the opened depammo netcdf file:
         do loop2_=1,max_depammo_time_counter
          status=nf_get_vara_double(ncid_,var_id,loop2_,1,x3)
          if(status/=0)stop 'erreur nf_get_vara_double depammo'
          x2=-999.
          if(index(texte60,'days')/=0)x2=86400.
          if(index(texte60,'hours')/=0)x2=3600.
          if(index(texte60,'seconds')/=0)x2=1.
          if(x2==-999.)stop 'depammo unites de temps incomprises'
          timeforecast_=x3*x2
          time_=elapsedtime_out+timeforecast_
!          print*,'unites de temps',texte60,x2
!          print*,'unites de temps 2 ',timeforecast_,time_
! write in the ascii temporary file:
          write(4,'(a)')trim(texte80(1))  ! Nom fichier netcdf
          write(4,*)loop2_                ! Numero d'echeance dans le fichier netdf
          write(4,*)time_                 ! temps en seconde du champ
          write(4,*)timeforecast_         ! temps en seconde depuis le debut du run ecmwf
          write(4,*)flag_cumul            ! Cumul ou pas cumul ?
          write(4,*)scalefct              ! facteur d'echelle aditionnel

!          print*,'write file',loop2_,trim(texte80(1)),time_, &
!                  timeforecast_,flag_cumul,scalefct

         enddo ! loop2_

! Close meteo netcdf file:
          status=nf_close(ncid_)

          enddo ! loop3_
 332     close(4)
         close(3)

! C'est le proc zero qui a la mission d'assembler tous les fichiers tmp en un seul fichier liste binaire
! La barriere suivante permet de s'assurer que tous les fichiers individuels ont bien ete faits au moment
! d'attaquer la concatenation
      call mpi_barrier(par%comm2d,k_out)  
       if(par%rank==0) then !00000000000>
        time2_=-999.
!        print*,'depammobinreclist(loop1_):',depammobinreclist(loop1_),loop1_
        open(unit=3,file=depammobinreclist(loop1_) &
            ,access='direct'                      &
            ,recl=232                             &
            ,form='unformatted')

        nc=1
        do loop2_=0,nbdom-1
         write(texte60,'(a,i0)')'tmp/tmpdepammofile',loop2_
         open(unit=4,file=texte60)

 3366     time1_=time2_
          read(4,'(a)',end=3367)texte80(1)  ! Nom fichier netcdf
          read(4,*)k10                      ! Numero d'echeance dans le fichier netdf
          read(4,*)time2_                   ! temps en seconde du champ
          read(4,*)timeforecast_            ! temps en seconde depuis le debut du run ecmwf
          read(4,*)flag_cumul               ! Cumul ou pas cumul ?
          read(4,*)scalefct                 ! facteur d'echelle aditionnel
!          print*,'ec',texte80(1)
!          print*,'ecbis',k10,time2_,timeforecast_,flag_cumul,scalefct
          if(time1_==-999.)time1_=time2_
          if(flag_cumul==1.or.flag_depammo_average==1) then !>>>
           time_=0.5*(time1_+time2_)
          else                                            !>>>
           time_=time2_
          endif                                           !>>>

!          print*,'time_',time_,loop2_,elapsedtime_now,loop1_

          if(time_<=elapsedtime_now) then !>>>> ce test garantit entre autres que model3d ne relancera pas une
           depammofile_nextrec(loop1_)=nc  !     procedure de remise a jour des champs a la premiere iteration
           depammofile_nextime(loop1_)=time_
          endif                           !>>>>

!          print*,'ec2',depammofile_nextrec(loop1_),depammofile_nextime(loop1_),loop1_,loop2_
        
          write(3,rec=nc)texte80(1)       & ! Nom fichier netcdf
                        ,k10              & ! Numero d'echeance dans le fichier netdf
                        ,time_            & ! temps en seconde dans le repere de S26
                        ,timeforecast_    & ! temps en seconde depuis le debut du run ecmwf
                        ,flag_cumul       & ! Cumul ou pas cumul ?
                        ,scalefct           ! facteur d'echelle aditionnel

          nc=nc+1


          goto 3366
 3367     close(4)

        enddo ! loop2_
        close(3)
       endif                !00000000000>
      call mpi_barrier(par%comm2d,k_out)  

!      print*,'loop1_,depammofile',depammofile_nextrec(loop1_),depammofile_nextime(loop1_),par%rank
      enddo ! loop1_
        
! Share with all other ranks:
        allocate(tabr8_(ndepammo))
        call mpi_allreduce(depammofile_nextime(1:ndepammo)  &
                          ,tabr8_(1:ndepammo)              &
                          ,ndepammo                        &
                          ,mpi_double_precision           &
                          ,mpi_sum                        &
                          ,par%comm2d,ierr)
!        print*,'tabr8_=',tabr8_(1:ndepammo)
        depammofile_nextime(1:ndepammo)=tabr8_(1:ndepammo)
        deallocate(tabr8_)
        allocate(tabi4_(ndepammo))
        call mpi_allreduce(depammofile_nextrec(1:ndepammo)  &
                          ,tabi4_(1:ndepammo)              &
                          ,ndepammo                        &
                          ,mpi_integer                    &
                          ,mpi_sum                        &
                          ,par%comm2d,ierr)
       depammofile_nextrec(1:ndepammo)=tabi4_(1:ndepammo)
       deallocate(tabi4_)

      call mpi_barrier(par%comm2d,k_out)  !13-10-10
      cpu_seconds=MPI_Wtime ( ) - cpu_seconds
      if(par%rank==0) &
      write(6,*)'depammoflux_binary_file_listcpu_seconds=',cpu_seconds

      end subroutine depammo_binary_file

!..........................................................

      subroutine depammoflux_varname(loop_)
      implicit none
      integer loop_

          if(loop_==dd_id )then;depammovarname(1)='ddflx';return;endif
          if(loop_==wdr_id)then;depammovarname(1)='wdrflx';return;endif
          if(loop_==wdw_id)then;depammovarname(1)='wdwflx';return;endif
      write(6,*)'loop',loop_,depammovarname(1),dd_id,wdr_id,wdw_id 
      stop ' STOP depammoflux_varname varname non trouve cas glorys'

      end subroutine depammoflux_varname
!..........................................................


      subroutine depammoflux_inquire_var(loop_,ncid_)
      implicit none
      integer loop_,ncid_,varid_

      flag_cumul=0 ! pas de cumul
      scalefct=1.  ! pas de chgt d'unites
      texte90=''   ! reset

      if( loop_==dd_id.or.   &
          loop_==wdr_id.or.  &
          loop_==wdw_id) then                          !*****>


                   status=nf_inq_varid(ncid_,depammovarname(1),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,depammovarname(2),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,depammovarname(3),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,depammovarname(4),varid_)
      if(status/=0) then !----->
       write(6,*)'dd_id wdr_id wdw_id loop_=' &
                 ,dd_id,wdr_id,wdw_id,loop_
       write(6,'(6(a,1x))')trim(depammovarname(1)) &
                          ,trim(depammovarname(2)) &
                          ,trim(depammovarname(3)) &
                          ,trim(depammovarname(4)) &
                          ,'not found in the netcdf file' & 
                          ,trim(texte80(1))
       stop 'Erreur nf_inq_varid depammoflux_inquire_var 1'
      endif              !----->

      status=nf_get_att_text(ncid_,varid_,'units',texte90);if(status/=0)stop 'echec attribut units pour dd'
      if(   texte90=='W m**-2 s'            &
        .or.texte90=='J m**-2')flag_cumul=1 ! cumul


      endif                                                           !*****>


      end subroutine depammoflux_inquire_var

!..........................................................

      subroutine depammoflux_extract_zone
      implicit none

!      if(    flag_meteodata=='ecmwf'                           &
!         .or.flag_meteodata=='glorys') then !ffffffffffffffff>
!
        texte30=depammobinreclist(1)
        nc=1
!
!      else                                  !ffffffffffffffff>
!
!      x2=(elapsedtime_now-airseadt(1,2))/airseadt(1,1)
!      nc=1+int(x2)
!
!        if(nc.le.0)then !---- debug ----->                           !06/05/04
!        if(par%rank==0)write(6,*)'Erreur l 3150'
!        if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 2'
!        if(par%rank==0)write(6,*)'car nc=',nc,' est <= 0 qui signifie que je suis'
!        if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
!        if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc+1 e 1'
!        ncmin_airsea=1
!        nc=ncmin_airsea
!        pause
!        endif            !---- debug ----->
!
!        k1=1 ; write(texte30,'(a14,i0)')'tmp/listemeteo',k1
!
!      endif                                 !ffffffffffffffff>

        open(unit=4,file=texte30                                        &
                   ,access='direct'          &
                   ,recl=232                 &
                   ,form='unformatted')
        read(4,rec=nc)texte80(1),i
        close(4)

        if(par%rank==0)write(6,'(a,a)')                    &
                'Fichier lu par depammoflux_extract_zone: ' &
               ,trim(texte80(1))

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec fichier 1'

                   status=nf_inq_dimid(ncid1,'lon',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'g0_lon_1',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jx',dim_x_id)
      if(status/=0)stop 'Erreur dim_x_id variable depammo longitude'

                   status=nf_inq_dimid(ncid1,'g0_lat_0',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'iy',dim_y_id)
      if(status/=0)stop 'Erreur dim_y_id variable depammo latitude'

      status=nf_inq_dimlen(ncid1,dim_x_id,depammo_imax)
      status=nf_inq_dimlen(ncid1,dim_y_id,depammo_jmax)

! on part sur l'hypothese que l'on connait deje depammo_kmax
      depammo_kmax=2

! une premiere allocation pour pouvoir charger les tableaux lon lat
      call allocate_forcages(1,13,depammo_imax,depammo_jmax,depammo_kmax) ! arg1=allouer arg2=meteo arg3,4,5=dimensions

! lire les longitudes:
                   status=nf_inq_varid(ncid1,'lon',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'g0_lon_1',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'nav_lon',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'xlon',lon_id)
      if(status/=0)stop 'Erreur lon_id variable meteo longitude'

      status=nf_inq_var(ncid1,lon_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'Erreur nf_inq_var meteo'

      if(var_type/=nf_float.and.var_type/=nf_double) &
      stop 'Erreur type sur longitude meteo'

      if(var_dims==1) then
       i1=depammo_imax ; j1=1 ; varstart(1)=1 ; varcount(1)=depammo_imax
      endif
      if(var_dims==2) then
       i1=depammo_imax ; j1=depammo_jmax
       varstart(1)=1 ; varcount(1)=depammo_imax ; varstart(2)=1 ; varcount(2)=depammo_jmax
      endif

      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lon_id,varstart(1:var_dims)   &
                                              ,varcount(1:var_dims)   &
                                             ,depammo_lon(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(depammo_lonlat_r4(depammo_imax,depammo_jmax))
       status=nf_get_vara_real(ncid1,lon_id,varstart(1:var_dims)      &
                                            ,varcount(1:var_dims)      &
                                           ,depammo_lonlat_r4(1:i1,1:j1))
       depammo_lon(1:i1,1:j1)=depammo_lonlat_r4(1:i1,1:j1)
       deallocate(depammo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture longitude meteo'

! lire les latitudes:
                   status=nf_inq_varid(ncid1,'lat',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'g0_lat_0',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'nav_lat',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'xlat',lat_id)
      if(status/=0)stop 'Erreur lat_id variable meteo latitude'

      status=nf_inq_var(ncid1,lat_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type/=nf_float.and.var_type/=nf_double) &
      stop 'Erreur type sur latitude meteo'

      if(var_dims==1) then
       i1=1 ; j1=depammo_jmax ; varstart(1)=1 ; varcount(1)=depammo_jmax
      endif
      if(var_dims==2) then
       i1=depammo_imax ; j1=depammo_jmax
       varstart(1)=1 ; varcount(1)=depammo_imax ; varstart(2)=1 ; varcount(2)=depammo_jmax
      endif

      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lat_id,varstart(1:var_dims)   &
                                              ,varcount(1:var_dims)   &
                                             ,depammo_lat(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(depammo_lonlat_r4(depammo_imax,depammo_jmax))
       status=nf_get_vara_real(ncid1,lat_id,varstart(1:var_dims)        &
                                            ,varcount(1:var_dims)        &
                                           ,depammo_lonlat_r4(1:i1,1:j1))
       depammo_lat(1:i1,1:j1)=depammo_lonlat_r4(1:i1,1:j1)
       deallocate(depammo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture latitude meteo'


      if(var_dims==1) then !1111111111111111111111>
       do j=1,depammo_jmax
       do i=1,depammo_imax
        depammo_lon(i,j)=depammo_lon(i,1)
        depammo_lat(i,j)=depammo_lat(1,j)
       enddo
       enddo
       depammo_lonstr=depammo_lon(1,1)
       depammo_latstr=depammo_lat(1,1)
       depammo_londlt=(depammo_lon(depammo_imax,1)-depammo_lon(1,1))/depammo_imax
       depammo_latdlt=(depammo_lat(1,depammo_jmax)-depammo_lat(1,1))/depammo_jmax
      endif                !1111111111111111111111>

! A partir des min & max des (lon,lat) des 2 grilles, reduire la zone d'extraction:
       depammozoom_istr=999999 ; depammozoom_iend=-999999 ! first guess
       depammozoom_jstr=999999 ; depammozoom_jend=-999999 ! first guess

       ksecu=0
       do j=1,depammo_jmax
       do i=1,depammo_imax

       if(depammo_lon(i,j)>=lonmin.and.depammo_lon(i,j)<=lonmax.and.   &
          depammo_lat(i,j)>=latmin.and.depammo_lat(i,j)<=latmax)then

         depammozoom_istr=min(depammozoom_istr,i)
         depammozoom_jstr=min(depammozoom_jstr,j)
         depammozoom_iend=max(depammozoom_iend,i)
         depammozoom_jend=max(depammozoom_jend,j)
         ksecu=1

       endif

       enddo
       enddo

! Si ksecu=0 c'est que le proc est si petit, qu'aucun point depammo ne s'est trouve e l'interieur
! On refait le test differement. On cherche le point le plus proche.
       if(ksecu==0)then !000000000000000>
        dist1=1.d20 ; i1=imax/2 ; j1=jmax/2
        do j=1,depammo_jmax
        do i=1,depammo_imax

         dist2=rayonterre*                                         &
         acos( sin(depammo_lat(i,j)*deg2rad)*sin(lat_t(i1,j1))     &
              +cos(depammo_lat(i,j)*deg2rad)*cos(lat_t(i1,j1))     &
              *cos(lon_t(i1,j1)-depammo_lon(i,j)*deg2rad))

         if(dist2<dist1)then !--->
          dist1=dist2 ; i2=i ; j2=j
         endif               !--->

        enddo
        enddo
        depammozoom_istr=i2 ; depammozoom_iend=i2
        depammozoom_jstr=j2 ; depammozoom_jend=j2
       endif            !000000000000000>


! on elargit un peu plus pour boucher les trous sans etre restreint par la taille
! reduite de la zone d'extraction, et puis aussi pour rattraper l'erreur liee au fait
! que plusieurs grilles depammo (point u, v, t) peuvent etre presentes.
       i0=10 ! elargissement e 10 lignes / 10 colonnes supplementaires - repere 1233
       depammozoom_istr=max0(depammozoom_istr-i0,1)
       depammozoom_iend=min0(depammozoom_iend+i0,depammo_imax)
       depammozoom_jstr=max0(depammozoom_jstr-i0,1)
       depammozoom_jend=min0(depammozoom_jend+i0,depammo_jmax)

       depammo_imax=depammozoom_iend-depammozoom_istr+1
       depammo_jmax=depammozoom_jend-depammozoom_jstr+1

!      write(*,*)'depammozoom_istr=',depammozoom_istr
!      write(*,*)'depammozoom_jstr=',depammozoom_jstr
!      write(*,*)'depammozoom_iend=',depammozoom_iend
!      write(*,*)'depammozoom_jend=',depammozoom_jend
!      write(*,*)'depammo_imax=',depammo_imax
!      write(*,*)'depammo_jmax=',depammo_jmax

! Refaire l'allocation dynamique en fonction des dimensions reduites:
      call allocate_forcages(2,13,0,0,0) ! arg1=desallouer arg2=depammo
      call allocate_forcages(1,13,depammo_imax,depammo_jmax,depammo_kmax) ! arg1=allouer arg2=depammo arg3,4,5=dimensions

! Fermer le fichier netcdf:
      status=nf_close(ncid1)

!      if(flag_depammodata=='glorys') then !---------->
!
       call depammoflux_sgrid_to_atmgrid_driver ! 

!        call airseaflux_bouchetrou_driver          ! glorys case
!
!      else                              !---------->
!
!! Identifier les points de la grille depammo e la fois en terre (sur grille depammo)
!! et en mer (sur grille symphonie):
!        call fichier_bouchetrou_depammo
!
!      endif                             !---------->

      end subroutine depammoflux_extract_zone

!..........................................................


      subroutine depammoflux_sgrid_to_atmgrid_driver
      implicit none
      integer loop_

        open(unit=4,file=depammobinreclist(1) &
                   ,access='direct'              &
                   ,recl=232                     &
                   ,form='unformatted')
        read(4,rec=1)texte80(1),i
        close(4)

       call depammoflux_sgrid_to_atmgrid_perform('t')



      end subroutine depammoflux_sgrid_to_atmgrid_driver

!----------------------------------------------------------------------------

      subroutine depammoflux_sgrid_to_atmgrid_perform(text_)
      implicit none
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)    &
                        ,dlon_di_,dlon_dj_,dlon_dm_
      character*1 text_

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status/=0)stop &
      'echec fichier depammoflux_sgrid_to_atmgrid_perform'

! Recharger la longitude et la latitude depuis la zone reduite:
! Lire Longitude:
      if(var_dims==1) then
       i1=depammo_imax ; j1=1 ; varstart(1)=depammozoom_istr ; varcount(1)=depammo_imax
      endif
      if(var_dims==2) then
       i1=depammo_imax ; j1=depammo_jmax
       varstart(1)=depammozoom_istr ; varcount(1)=depammo_imax
       varstart(2)=depammozoom_jstr ; varcount(2)=depammo_jmax
      endif
      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lon_id,varstart(1:var_dims)   &
                                             ,varcount(1:var_dims)   &
                                             ,depammo_lon(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(depammo_lonlat_r4(depammo_imax,depammo_jmax))
       status=nf_get_vara_real(ncid1,lon_id,varstart(1:var_dims)      &
                                           ,varcount(1:var_dims)      &
                                           ,depammo_lonlat_r4(1:i1,1:j1))
       depammo_lon(1:i1,1:j1)=depammo_lonlat_r4(1:i1,1:j1)
       deallocate(depammo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture longitude depammo'

! Lire Latitude:
      if(var_dims==1) then
       i1=1 ; j1=depammo_jmax ; varstart(1)=depammozoom_jstr ; varcount(1)=depammo_jmax
      endif
      if(var_dims==2) then
       i1=depammo_imax ; j1=depammo_jmax
       varstart(1)=depammozoom_istr ; varcount(1)=depammo_imax
       varstart(2)=depammozoom_jstr ; varcount(2)=depammo_jmax
      endif

      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lat_id,varstart(1:var_dims)   &
                                             ,varcount(1:var_dims)   &
                                             ,depammo_lat(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(depammo_lonlat_r4(depammo_imax,depammo_jmax))
       status=nf_get_vara_real(ncid1,lat_id,varstart(1:var_dims)        &
                                           ,varcount(1:var_dims)        &
                                           ,depammo_lonlat_r4(1:i1,1:j1))
       depammo_lat(1:i1,1:j1)=depammo_lonlat_r4(1:i1,1:j1)
       deallocate(depammo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture latitude depammo'

! Fermer le fichier netcdf:
      status=nf_close(ncid1)

       lon_1_1_gridt=depammo_lon(1,1) ; lat_1_1_gridt=depammo_lat(1,1)
       lon_2_1_gridt=depammo_lon(2,1) ; lat_1_2_gridt=depammo_lat(1,2)

!      write(10+par%rank,*),'passe ici'
      if(.not.allocated(ij2depammo_i))  allocate(ij2depammo_i   (0:imax+1,0:jmax+1))
      if(.not.allocated(ij2depammo_j))  allocate(ij2depammo_j   (0:imax+1,0:jmax+1))


      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

      x2=real(depammo_imax/2)
      x3=real(depammo_jmax/2)
      deci=x2
      decj=x3

! ETAPE 1: trouver les coordonnees dans la grille ORCA:

! First guess: centre du domaine:
      k10=0
 1456 continue

! Principe suppose une relation lineaire entre lat lon et indice de grille.
!      dlon/di*Di+dlon/dj*Dj=Dlon
!      dlat/di*Di+dlat/dj*Dj=Dlat
! On cherche Di et Dj correspondant e Dlon=lon(i,j)-londepammo(i0,j0)
!                                et e Dlat=lat(i,j)-latdepammo(i0,j0)

      i1=int(deci)
      j1=int(decj)

      do j2=0,1
      do i2=0,1
      i0=i1+i2
      j0=j1+j2

      dlon_di_=depammo_lon(i0+1,j0  )-depammo_lon(i0-1,j0  )
      dlon_dj_=depammo_lon(i0  ,j0+1)-depammo_lon(i0  ,j0-1)
      dlon_dm_=rad2deg*lon_t(i  ,j)-depammo_lon(i0  ,j0)

      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(dlon_dm_<-180.)dlon_dm_=dlon_dm_+360.
      if(dlon_dm_> 180.)dlon_dm_=dlon_dm_-360.

! Determinant principal:
      x1=( dlon_di_                                             &
          *(depammo_lat(i0  ,j0+1)-depammo_lat(i0  ,j0-1))          &
          -(depammo_lat(i0+1,j0  )-depammo_lat(i0-1,j0  ))          &
          *dlon_dj_ )*0.25

      deci_(i2,j2)=min(max(                                  &
       i0+( dlon_dm_                                            &
          *(depammo_lat(i0  ,j0+1)-depammo_lat(i0,j0-1))            &
          -(rad2deg*lat_t(i,j)  -depammo_lat(i0,j0))              &
          *dlon_dj_)/x1*0.5    &
                   ,2.0001d0),depammo_imax-1.0001d0)

      decj_(i2,j2)=min(max(                                  &
       j0+( dlon_di_                                            &
          *(rad2deg*lat_t(i,j)  -depammo_lat(i0  ,j0))            &
          -(depammo_lat(i0+1,j0  )-depammo_lat(i0-1,j0))            &
          *dlon_dm_   )/x1*0.5    &
                   ,2.0001d0),depammo_jmax-1.0001d0)

      enddo
      enddo

      rapi=deci-i1
      rapj=decj-j1

      deci=(1.-rapi)*(1.-rapj)*deci_(0,0)   &
          +(1.-rapi)*    rapj *deci_(0,1)   &
          +    rapi *    rapj *deci_(1,1)   &
          +    rapi *(1.-rapj)*deci_(1,0)
      decj=(1.-rapi)*(1.-rapj)*decj_(0,0)   &
          +(1.-rapi)*    rapj *decj_(0,1)   &
          +    rapi *    rapj *decj_(1,1)   &
          +    rapi *(1.-rapj)*decj_(1,0)

! Si le point vise est different du first guess refaire le calcul
! avec un first guess donne par le dernier point vise:
      if(sqrt( (deci-x2)**2+(decj-x3)**2 ).gt.0.001)then
       x2=deci
       x3=decj
       k10=k10+1
       if(k10>20)then !!!!!!>
         if(par%rank==0)write(6,*)'par%rank=',par%rank
         if(par%rank==0)write(6,*)'(i,j)=   ',i,j
         if(par%rank==0)write(6,*)'deci decj',deci,decj
         stop 'hr_to_lr ne converge pas dans le forfait'
       endif          !!!!!!>
       goto 1456
      endif

! ETAPE 2: CALCULER L'ANGLE D'ORIENTATION LOCALE DE LA GRILLE ORCA

      i1=int(deci) ; j1=int(decj)
      rapi=deci-i1 ; rapj=decj-j1
      do j2=0,1
      do i2=0,1
       i0=i1+i2
       j0=j1+j2

       dlon_di_=depammo_lon(i0+1,j0  )-depammo_lon(i0-1,j0  )
       if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
       if(dlon_di_> 180.)dlon_di_=dlon_di_-360.

       dy_(i2,j2)= depammo_lat(i0+1,j0)-depammo_lat(i0-1,j0)
       dx_(i2,j2)=dlon_di_           &
                            *cos(depammo_lat(i0,j0)*deg2rad)
      enddo
      enddo
      x1=(1.-rapi)*(1.-rapj)*dx_(0,0)   &
        +(1.-rapi)*    rapj *dx_(0,1)   &
        +    rapi *    rapj *dx_(1,1)   &
        +    rapi *(1.-rapj)*dx_(1,0)
      y1=(1.-rapi)*(1.-rapj)*dy_(0,0)   &
        +(1.-rapi)*    rapj *dy_(0,1)   &
        +    rapi *    rapj *dy_(1,1)   &
        +    rapi *(1.-rapj)*dy_(1,0)

      ij2depammo_i(i,j)=deci
      ij2depammo_j(i,j)=decj

      enddo
      enddo

      end subroutine depammoflux_sgrid_to_atmgrid_perform

!..........................................................

      subroutine depammoflux_interp(t_,case_)
      implicit none
      integer loop_,t_,vstart1_,vstart2_,flag_rotation_,filval_,case_
      double precision time1_,time2_
      real i_shift_,j_shift_,filvalr4_
      character*1 text_


      do loop_=1,ndepammo

      call depammoflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14

      if(decision==1) then                            !*******************>

      call depammoflux_get_time_from_binrecfile(loop_,t_,vstart1_,vstart2_,time1_,time2_)

! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      if(status.ne.0)stop 'echec fichier 2'

      if(loop_==dd_id) then
        text_='t' ; varcode(loop_)=1 ; 
                   status=nf_inq_varid(ncid1,'ddflx',var_id)
      endif
      if(loop_==wdr_id) then
        text_='t' ; varcode(loop_)=2
                   status=nf_inq_varid(ncid1,'wdrflx',var_id)
      endif
      if(loop_==wdw_id) then
        text_='t' ; varcode(loop_)=3
                   status=nf_inq_varid(ncid1,'wdwflx',var_id)
      endif

      if(status/=0) then
!      write(6,'(a,a)')'Pour variable:',trim(depammoflux_varname(var_id_))
       stop 'erreur varname depammoflux_interp'
      endif

      status=nf_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
!     if(var_type/=nf_short) then !------->
!      write(6,*)'loop_=',loop_
!      stop &
!     'cas var_type pas prevu airseaflux_glorys_interp'
!     endif                       !------->

      varstart(1)=depammozoom_istr ; varcount(1)=depammo_imax
      varstart(2)=depammozoom_jstr ; varcount(2)=depammo_jmax
      varstart(3)=vstart1_       ; varcount(3)=1
      varstart(4)=1              ; varcount(4)=1

      ksecu=0
      if(var_type==nf_real) then  !rrrrrrrr>
      ksecu=1
      status=nf_get_vara_real(ncid1,var_id,varstart(1:var_dims)        &
                                          ,varcount(1:var_dims)        &
                          ,depammo_var(1:depammo_imax,1:depammo_jmax))
! A priori pas de var_scalefactor var_addoffset dans archivage real....:
      var_scalefactor=1. ; var_addoffset=0.
      endif                       !rrrrrrrr>
      if(ksecu==0)stop' Erreur type 3905 module_depammo'

      texte90='' ; status=nf_get_att_text(ncid1,var_id,'units',texte90)
      if(status/=0)stop 'echec att units depammoflux_interp'

! Fermer le fichier netcdf
      status=nf_close(ncid1)


      if(varcode(loop_)==1) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2depammo_i(i,j) ; i1=int(deci) ; rapi=deci-i1
        decj=ij2depammo_j(i,j) ; j1=int(decj) ; rapj=decj-j1
            ddammo_w(i,j,0)=ddammo_w(i,j,2)
            ddammo_w(i,j,2)=((1.-rapi)*(1.-rapj)*depammo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *depammo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*depammo_var(i1+1,j1  )  &
                         +    rapi *    rapj *depammo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset

!        if(par%rank==41) then
!        print*,'sgrid',i,j,ddammo_w(i,j,0),ddammo_w(i,j,2)
!        print*,'sgrid',i,j,lat_t(i,j)*rad2deg,lon_t(i,j)*rad2deg
!        print*,'regcmgrid',i1,j1,depammo_var(i1  ,j1  )
!        print*,'regcmgrid',i1,j1,depammo_lat(i1  ,j1  ),depammo_lon(i1  ,j1  )
!        endif

       enddo ; enddo
      call get_type_echange('z0','ddammo_w_z0_2',ddammo_w,lbound(ddammo_w),ubound(ddammo_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(ddammo_w,i2,mpi_neighbor_list(loop3)) !31-07-14
      enddo
      call loc_wait()
      endif             !----->

      if(varcode(loop_)==2) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2depammo_i(i,j) ; i1=int(deci) ; rapi=deci-i1
        decj=ij2depammo_j(i,j) ; j1=int(decj) ; rapj=decj-j1
           wdrammo_w(i,j,0)=wdrammo_w(i,j,2)
           wdrammo_w(i,j,2)=((1.-rapi)*(1.-rapj)*depammo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *depammo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*depammo_var(i1+1,j1  )  &
                         +    rapi *    rapj *depammo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
      call get_type_echange('z0','wdrammo_w_z0_2',wdrammo_w,lbound(wdrammo_w),ubound(wdrammo_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(wdrammo_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      endif             !----->

      if(varcode(loop_)==3) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2depammo_i(i,j) ; i1=int(deci) ; rapi=deci-i1
        decj=ij2depammo_j(i,j) ; j1=int(decj) ; rapj=decj-j1
           wdwammo_w(i,j,0)=wdwammo_w(i,j,2)
           wdwammo_w(i,j,2)=((1.-rapi)*(1.-rapj)*depammo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *depammo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*depammo_var(i1+1,j1  )  &
                         +    rapi *    rapj *depammo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
      call get_type_echange('z0','wdwammo_w_z0_2',wdwammo_w,lbound(wdwammo_w),ubound(wdwammo_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(wdwammo_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      endif             !----->

      endif                                           !*******************>    !

      enddo ! find de boucle sur loop_


      end subroutine depammoflux_interp

!..........................................................


      subroutine depammoflux_flux_moveforward ! time linear interpolation
      implicit none
      integer loop_

      do loop_=1,ndepammo

      x2=(elapsedtime_now        -depammofile_prvtime(loop_))      &
        /(depammofile_nextime(loop_)-depammofile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2

!     print*,'x0',x0,x2

      if(varcode(loop_)==1) then
       do j=1,jmax ; do i=1,imax
        ddammo_w(i,j,1)=x0*ddammo_w(i,j,0)+x2*ddammo_w(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==2) then
       do j=1,jmax ; do i=1,imax
        wdrammo_w(i,j,1)=x0*wdrammo_w(i,j,0)+x2*wdrammo_w(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==3) then
       do j=1,jmax ; do i=1,imax
        wdwammo_w(i,j,1)=x0*wdwammo_w(i,j,0)+x2*wdwammo_w(i,j,2)
       enddo ; enddo
      endif

      enddo ! fin de boucle sur loop_


      end subroutine depammoflux_flux_moveforward ! time linear interpolation

!.............................................................

      subroutine depammoflux_decision(loop_,case_)
      implicit none
      integer loop_,case_

! decision=1 indique qu'il est temps de lire un nouveau fichier. decision=0 sinon.

      decision=0

! Si les proc ne sont pas synchro en raison du subcycling on s'abstient afin
! de garantir la continuite horizontale des champs interpoles:
      if(subcycle_synchro==0)return

! Jalon temps depasse = lire un nouveau fichier
      if(elapsedtime_now>depammofile_nextime(loop_))decision=1

! A l'etat initial lire imperativement
      if(case_==1)decision=1 ! Etat initial

      end subroutine depammoflux_decision

!..........................................................

      subroutine depammoflux_interp_driver(var_id_,t_,var_id_file_)
      implicit none
      integer var_id_,t_,vstart1_,vstart2_,ncid_,var_id_file_
      double precision time1_,time2_

! Obtenir le nom du fichier netcdf, texte80, et les numeros d'echeances time_ et time2_
      call depammoflux_get_time_from_binrecfile(var_id_file_,t_,vstart1_,vstart2_,time1_,time2_)

      call depammoflux_varname(var_id_) ! donne meteovarname le nom de la variable du fichier netcdf

      call depammoflux_read_only(vstart1_,1)

! Cas de variables cumules:
!      if(flag_cumul==1) then !1111111111111111>

!      k0=0
!      if(depammo_t0=='file'.and.vstart1_>vstart2_)k0=1
!      if(depammo_t0=='time'.and.time1_>time2_)    k0=1
     

!        if(k0==1) then !----------------> 

!            do j=1,depammo_jmax ; do i=1,depammo_imax
!             depammo_cum(i,j)=depammo_var(i,j)
!            enddo ; enddo
!            call depammoflux_read_only(vstart2_,2)

!            x1=scalefct/(time1_-time2_)
!            do j=1,depammo_jmax ; do i=1,depammo_imax
!             depammo_var(i,j)=(depammo_cum(i,j)-depammo_var(i,j))*x1
!            enddo ; enddo
!
!        else           !---------------->

!            if(time1_==0.)stop ' Stop depammoflux_interp_driver time1==0'
!            x1=scalefct/time1_
!            do j=1,depammo_jmax ; do i=1,depammo_imax
!             depammo_var(i,j)=depammo_var(i,j)*x1
!            enddo ; enddo
!
!        endif          !---------------->

!      endif                  !1111111111111111>

      call depammoflux_interpolation(var_id_) ! calcule xy_t(i,j,1), le champs interpole

! La ligne suivante est deplacee e la fin de depammoflux_get_time_from_binrecfile....
!     depammofile_nextrec(var_id_)=file_nextrec(var_id_)+1

      end subroutine depammoflux_interp_driver

!..........................................................


      subroutine depammoflux_get_time_from_binrecfile(loop_,t_,vstart1_,vstart2_,time1_,time2_)
      implicit none
      integer loop_,t_,vstart1_,vstart2_
      double precision time1_,time2_

      nc=depammofile_nextrec(loop_)
!     print*,'nc',nc

      open(unit=4,file=depammobinreclist(loop_)                     &
                 ,access='direct',recl=232,form='unformatted')
      read(4,rec=max(nc-1,1))texte80(2),vstart2_                   &
                                       ,depammofile_prvtime(loop_)  & ! Echeance avant
                                       ,time2_
      read(4,rec=nc         )texte80(1),vstart1_                   &
                                       ,depammofile_nextime(loop_)  & ! Echeance apres
                                       ,time1_                     &
                                       ,flag_cumul                 &
                                       ,scalefct
      close(4)

! Verifications de l'etat initial:
      if(iteration3d==0) then !0000000000000000000>
        if(t_==0.and.depammofile_nextime(loop_)>elapsedtime_now) then  !------->
         stop 'erreur1 depammoflux_get_time_from_binrecfile'
        endif                                                         !------->
        if(t_==2.and.depammofile_nextime(loop_)<elapsedtime_now) then  !------->
         stop 'erreur2 depammoflux_get_time_from_binrecfile'
        endif                                                         !------->
        if(t_==2.and.depammofile_nextime(loop_)==elapsedtime_now) then !------->
         stop 'erreur3 depammoflux_get_time_from_binrecfile'
        endif                                                         !------->
      endif                   !0000000000000000000>

! Move forward du prochain numero de record:
      depammofile_nextrec(loop_)=depammofile_nextrec(loop_)+1

      end subroutine depammoflux_get_time_from_binrecfile

!..............................................................................


      subroutine depammoflux_read_only(vstart_,ttx_)
      implicit none
      integer vstart_,ttx_,ncid_

      status=nf_open(trim(texte80(ttx_)),nf_nowrite,ncid_);if(status/=0)stop 'erreur nf_open depammoflux_read_only'
                   status=nf_inq_varid(ncid_,depammovarname(1),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,depammovarname(2),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,depammovarname(3),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,depammovarname(4),var_id)
      if(status/=0) then !------>
       write(6,'(6(a,1x))')trim(depammovarname(1)) &
                          ,trim(depammovarname(2)) &
                          ,trim(depammovarname(3)) &
                          ,trim(depammovarname(4)) &
                          ,'not found in the netcdf file' & 
                          ,trim(texte80(ttx_))
       stop 'Erreur nf_inq_varid depammovarname 807'
      endif

      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status.ne.0)stop 'echec nf_inq_var 822'

      varstart(3)=vstart_        ; varcount(3)=1
      varstart(2)=depammozoom_jstr ; varcount(2)=depammo_jmax ! lat
      varstart(1)=depammozoom_istr ; varcount(1)=depammo_imax ! lon

      if(var_type==nf_real) then  !rrrrr>
       status=nf_get_vara_real(ncid_,var_id        &
                                    ,varstart(1:var_dims) &
                                    ,varcount(1:var_dims) &
             ,depammo_var(1:depammo_imax,1:depammo_jmax))
       if(status/=0)stop 'Erreur nf_get_vara_real depammoflux_read_only'
      endif                       !rrrrr>
      if(var_type==nf_short) then !iiiii>
       status=nf_get_vara_int(ncid_,var_id        &
                                   ,varstart(1:var_dims) &
                                   ,varcount(1:var_dims) &
          ,depammo_short(1:depammo_imax,1:depammo_jmax))
       if(status/=0) then !----->
        write(6,'(a,a)')'netcdf file=',trim(texte80(ttx_))
        write(6,'(5a)')'depammovarname ',depammovarname(:)
        write(6,*)'depammo_imax,depammo_jmax ',depammo_imax,depammo_jmax
        write(6,*)'ubound(depammo_short) ',ubound(depammo_short)
        write(6,*)'vardims=',var_dims
        write(6,*)'vstart_ ',vstart_
        write(6,*)'varstart ',varstart(1:var_dims)
        write(6,*)'varcount ',varcount(1:var_dims)
        stop 'Erreur nf_get_vara_int depammoflux_read_only'
       endif              !----->
       status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
       if(status/=0) &
       stop 'error get scale_factor 843'
       status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
       if(status/=0)stop 'error get add_offset 843'
       do j=1,depammo_jmax ; do i=1,depammo_imax
        depammo_var(i,j)=depammo_short(i,j)*var_scalefactor+var_addoffset
       enddo ; enddo
      endif                      !iiiii>


      status=nf_close(ncid_)                           ;if(status/=0)stop 'erreur nf_close depammoflux_read_only'


      if(par%rank==0)write(6,'(a,a5,a,a,i4)')'depammoflux_read_only ' &
                                            ,depammovarname(1),' '     &
                                            ,trim(texte80(ttx_)),vstart_

      end subroutine depammoflux_read_only

!..............................................................................


      subroutine depammoflux_interpolation(var_id_)
      implicit none
      integer var_id_

!! Bouchage des trous pour certaines variables (vent exclus):
!      k0=1
!      if(var_id_==u10m_id) k0=0
!      if(var_id_==v10m_id) k0=0
!      if(var_id_==u100m_id)k0=0
!      if(var_id_==v100m_id)k0=0
!      if(var_id_==abl_id)  k0=0
!                       if(k0==1)call appli_bouchetrou_depammo

! Interpoler:
      const1=1./depammo_londlt
      const2=1./depammo_latdlt
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
       deci=1.+( lon_t(i,j)*rad2deg-depammo_lonstr)*const1
       i1=int(deci)
       rapi=deci-i1
       i1=i1-depammozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
       decj=1.+( lat_t(i,j)*rad2deg-depammo_latstr)*const2
       j1=int(decj)
       rapj=decj-j1
       j1=j1-depammozoom_jstr+1                              ! conservation la parallelisation

       xy_t(i,j,1)= (1.-rapi)*(1.-rapj)*depammo_var(i1  ,j1  )         &
                   +(1.-rapi)*    rapj *depammo_var(i1  ,j1+1)         &
                   +    rapi *(1.-rapj)*depammo_var(i1+1,j1  )         &
                   +    rapi *    rapj *depammo_var(i1+1,j1+1)

      enddo
      enddo

      end subroutine depammoflux_interpolation

!..............................................................................


!..............................................................................

      end module module_depammo
