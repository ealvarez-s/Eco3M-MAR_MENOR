     program lecture_notebook_rivers
      implicit none
      include 'netcdf.inc'
      character*100 file_grid_input,texte
      integer nriver,ncid2,imax,jmax,dim_x_id,dim_y_id,k,i,j,var_id,status
      real    runoff,freq,lect,annee,mois,jour,heure,minute,sec,ifile

      double precision ,dimension(:,:)  ,allocatable :: lon_t
      double precision ,dimension(:,:)  ,allocatable :: lat_t

! on lit lat,lon de globmed
      file_grid_input='/tmpdir/cestour/GLOBMED2/OFFLINE_SOL_ET_TIDE22/grid.nc'
      status=nf_open(file_grid_input,nf_nowrite,ncid2)
!  read idents of dimensions
       if(status.ne.0)stop'status open the file'
       status=nf_inq_dimid(ncid2,'ni_t',dim_x_id)
       if(status/=0)write(6,*)'status ni_t',status
       status=nf_inq_dimid(ncid2,'nj_t',dim_y_id)
       if(status/=0)write(6,*)'status nj_t',status
!  read dimensions
       status=nf_inq_dimlen(ncid2,dim_x_id,imax)
       if(status/=0)write(6,*)'status imax',status
       status=nf_inq_dimlen(ncid2,dim_y_id,jmax)
       if(status/=0)write(6,*)'status jmax',status
       write(6,*)imax,jmax


       if(.not.allocated(lat_t)) allocate(lat_t(0:imax-1,0:jmax-1))
       if(.not.allocated(lon_t)) allocate(lon_t(0:imax-1,0:jmax-1))
      status=nf_inq_varid(ncid2,'latitude_t',var_id)
       if(status/=0)write(6,*)'status lat_t',status
      status=nf_get_var_double(ncid2,var_id,lat_t)
       if(status/=0)write(6,*)'status lat_t 2 ',status
      status=nf_inq_varid(ncid2,'longitude_t',var_id)
       if(status/=0)write(6,*)'status lon_t',status
      status=nf_get_var_double(ncid2,var_id,lon_t)
       if(status/=0)write(6,*)'status lon_t 2 ',status
       status=nf_close(ncid2)



       open(unit=2,file='notebook_rivers7')
       open(unit=3,file='discharge_rivers7')

       read(2,*)nriver

       do k=1,nriver
       read(2,*)
       read(2,*)
       read(2,*)i
       read(2,*)j
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)runoff
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)
       read(2,*)ifile
       read(2,*)freq
       read(2,*)
       read(2,*)annee,mois,jour,heure,minute,sec

       
       print*,'i,j',i,j,lat_t(i,j),lon_t(i,j)
       write(3,'(2(i4,1x),11(f7.2,1x))')i,j,lat_t(i,j),lon_t(i,j),ifile,runoff,freq,annee,mois,jour,heure,minute,sec  


       enddo

       close(3)
       close(2)

       end
        

