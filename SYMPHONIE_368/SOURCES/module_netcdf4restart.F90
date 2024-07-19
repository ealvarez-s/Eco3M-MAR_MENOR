      module module_netcdf4restart
!______________________________________________________________________
! S model
! release S26  - last update: 14-12-15
!______________________________________________________________________
! version   date    description                                                !
! S.26    14-12-15 mise en service                                            !
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'
!     integer :: flag_netcdf_def=0
!     integer start(5),edge(5),shift_i_,shift_j_

      interface netcdf4restart_wrt
      module procedure netcdf4restart_wrt_5dr8
      module procedure netcdf4restart_wrt_5dr4
      module procedure netcdf4restart_wrt_4dr8
      module procedure netcdf4restart_wrt_4dr4
      module procedure netcdf4restart_wrt_3dr8
      module procedure netcdf4restart_wrt_3dr4
      module procedure netcdf4restart_wrt_3di4
      module procedure netcdf4restart_wrt_2dr8
      module procedure netcdf4restart_wrt_2dr4
      module procedure netcdf4restart_wrt_2di4
      module procedure netcdf4restart_wrt_1dr8
      module procedure netcdf4restart_wrt_1dr4
      module procedure netcdf4restart_wrt_1di4
      end interface 

      interface netcdf4restart_read
      module procedure netcdf4restart_read_5dr8
      module procedure netcdf4restart_read_5dr4
      module procedure netcdf4restart_read_4dr8
      module procedure netcdf4restart_read_4dr4
      module procedure netcdf4restart_read_3dr8
      module procedure netcdf4restart_read_3dr4
      module procedure netcdf4restart_read_3di4
      module procedure netcdf4restart_read_2dr8
      module procedure netcdf4restart_read_2dr4
      module procedure netcdf4restart_read_2di4
      module procedure netcdf4restart_read_1dr8
      module procedure netcdf4restart_read_1dr4
      module procedure netcdf4restart_read_1di4
      end interface 

      contains

!.......................................................................

      subroutine netcdf4restart_wrt_4dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_4dr8'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_4dr8'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_4dr8'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_4dr8'
      idim_=ub(4)-lb(4)+1 ; status=nf_def_dim(ncid_,'nt',idim_,dim_t_id)
      if(status/=0) stop 'Err 5 netcdf4restart_wrt_4dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id ; vardim(4)=dim_t_id
      status=nf_def_var(ncid_,nom,nf_double,4,vardim(1:4),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_4dr8'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_4dr8'

      status=nf_put_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_4dr8'


      end subroutine netcdf4restart_wrt_4dr8

!.......................................................................

      subroutine netcdf4restart_read_4dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      integer start(5),edge(5),shift_i_,shift_j_
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_4dr8'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_4dr8'

      status=nf_get_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_4dr8'

#ifdef bidon
      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

      shift_i_=1-lb(1)
      shift_j_=1-lb(2)
      istr=lb(1) ; iend=ub(1)
      jstr=lb(2) ; jend=ub(2)

      start(1)=istr+par%timax(1)+shift_i_
      start(2)=jstr+par%tjmax(1)+shift_j_
      edge(1)=iend-istr+1      
      edge(2)=jend-jstr+1

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_4dr8'


      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_4dr8'

      status=nf_get_vara_double(ncid_,var_id      &
                                     ,start(1:4)  &
                                      ,edge(1:4)  &
          ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_4dr8'
#endif



      end subroutine netcdf4restart_read_4dr8

!.......................................................................

      subroutine netcdf4restart_wrt_4dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      integer ncid_
!     integer(kind=MPI_OFFSET_KIND) idim_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_4dr4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_4dr4'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_4dr4'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_4dr4'
      idim_=ub(4)-lb(4)+1 ; status=nf_def_dim(ncid_,'nt',idim_,dim_t_id)
      if(status/=0) stop 'Err 5 netcdf4restart_wrt_4dr4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id ; vardim(4)=dim_t_id
      status=nf_def_var(ncid_,nom,nf_real,4,vardim(1:4),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_4dr4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_4dr4'

      status=nf_put_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_4dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_4dr4'


      end subroutine netcdf4restart_wrt_4dr4

!.......................................................................

      subroutine netcdf4restart_read_4dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      integer ncid_
!     integer(kind=MPI_OFFSET_KIND) idim_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_4dr4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_4dr4'

      status=nf_get_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_4dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_4dr4'

      end subroutine netcdf4restart_read_4dr4

!.......................................................................

      subroutine netcdf4restart_wrt_3dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      integer ncid_
!     integer(kind=MPI_OFFSET_KIND) idim_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_3dr8'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_3dr8'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_3dr8'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_3dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id 
      status=nf_def_var(ncid_,nom,nf_double,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_3dr8'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_3dr8'

      status=nf_put_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_3dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_3dr8'


      end subroutine netcdf4restart_wrt_3dr8

!.......................................................................

      subroutine netcdf4restart_read_3dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      integer ncid_
!     integer(kind=MPI_OFFSET_KIND) idim_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_3dr8'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_3dr8'

      status=nf_get_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_3dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_3dr8'

      end subroutine netcdf4restart_read_3dr8

!.......................................................................

      subroutine netcdf4restart_wrt_3dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      integer ncid_
      integer idim_


        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_3dr4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_3dr4'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_3dr4'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_3dr4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id 
      status=nf_def_var(ncid_,nom,nf_real,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_3dr4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_3dr4'

      status=nf_put_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_3dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_3dr4'


      end subroutine netcdf4restart_wrt_3dr4

!.......................................................................

      subroutine netcdf4restart_read_3dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_3dr4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_3dr4'

      status=nf_get_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_3dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_3dr4'

      end subroutine netcdf4restart_read_3dr4

!.......................................................................

      subroutine netcdf4restart_wrt_3di4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_3di4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_3di4'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_3di4'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_3di4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id 
      status=nf_def_var(ncid_,nom,nf_int,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_3di4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_3di4'

      status=nf_put_var_int(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_3di4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_3di4'


      end subroutine netcdf4restart_wrt_3di4

!.......................................................................

      subroutine netcdf4restart_read_3di4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_3di4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_3di4'

      status=nf_get_var_int(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_3di4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_3di4'

      end subroutine netcdf4restart_read_3di4

!.......................................................................

      subroutine netcdf4restart_wrt_2dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_2dr8'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_2dr8'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_2dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nf_def_var(ncid_,nom,nf_double,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_2dr8'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_2dr8'

      status=nf_put_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_2dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_2dr8'


      end subroutine netcdf4restart_wrt_2dr8

!.......................................................................

      subroutine netcdf4restart_read_2dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_2dr8'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_2dr8'

      status=nf_get_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_2dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_2dr8'

      end subroutine netcdf4restart_read_2dr8

!.......................................................................

      subroutine netcdf4restart_wrt_2dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_2dr4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_2dr4'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_2dr4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nf_def_var(ncid_,nom,nf_real,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_2dr4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_2dr4'

      status=nf_put_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_2dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_2dr4'


      end subroutine netcdf4restart_wrt_2dr4

!.......................................................................

      subroutine netcdf4restart_read_2dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_2dr4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_2dr4'

      status=nf_get_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_2dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_2dr4'

      end subroutine netcdf4restart_read_2dr4

!.......................................................................

      subroutine netcdf4restart_wrt_2di4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      integer ncid_
      integer idim_


        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_2di4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_2di4'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_2di4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nf_def_var(ncid_,nom,nf_int,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_2di4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_2di4'

      status=nf_put_var_int(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_2di4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_2di4'


      end subroutine netcdf4restart_wrt_2di4

!.......................................................................

      subroutine netcdf4restart_read_2di4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_2di4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_2di4'

      status=nf_get_var_int(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_2di4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_2di4'

      end subroutine netcdf4restart_read_2di4

!.......................................................................

      subroutine netcdf4restart_wrt_1dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(1) :: lb,ub
      double precision,dimension(lb(1):ub(1)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_1dr8'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_1dr8'

      vardim(1)=dim_x_id 
      status=nf_def_var(ncid_,nom,nf_double,1,vardim(1:1),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_1dr8'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_1dr8'

      status=nf_put_var_double(ncid_,var_id,tab(lb(1):ub(1)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_1dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_1dr8'


      end subroutine netcdf4restart_wrt_1dr8

!.......................................................................

      subroutine netcdf4restart_read_1dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(1) :: lb,ub
      double precision,dimension(lb(1):ub(1)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_1dr8'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_1dr8'

      status=nf_get_var_double(ncid_,var_id,tab(lb(1):ub(1)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_1dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_1dr8'

      end subroutine netcdf4restart_read_1dr8

!.......................................................................

      subroutine netcdf4restart_wrt_1dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(1) :: lb,ub
      real,dimension(lb(1):ub(1)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_1dr4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_1dr4'

      vardim(1)=dim_x_id 
      status=nf_def_var(ncid_,nom,nf_real,1,vardim(1:1),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_1dr4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_1dr4'

      status=nf_put_var_real(ncid_,var_id,tab(lb(1):ub(1)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_1dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_1dr4'


      end subroutine netcdf4restart_wrt_1dr4

!.......................................................................

      subroutine netcdf4restart_read_1dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(1) :: lb,ub
      real,dimension(lb(1):ub(1)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_1dr4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_1dr4'

      status=nf_get_var_real(ncid_,var_id,tab(lb(1):ub(1)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_1dr4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_1dr4'

      end subroutine netcdf4restart_read_1dr4

!.......................................................................

      subroutine netcdf4restart_wrt_1di4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(1) :: lb,ub
      integer,dimension(lb(1):ub(1)) :: tab
      integer ncid_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_1di4'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_1di4'

      vardim(1)=dim_x_id 
      status=nf_def_var(ncid_,nom,nf_int,1,vardim(1:1),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_1di4'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_1di4'

      status=nf_put_var_int(ncid_,var_id,tab(lb(1):ub(1)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_1di4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_1di4'


      end subroutine netcdf4restart_wrt_1di4

!.......................................................................

      subroutine netcdf4restart_read_1di4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(1) :: lb,ub
      integer,dimension(lb(1):ub(1)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_1di4'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_1di4'

      status=nf_get_var_int(ncid_,var_id,tab(lb(1):ub(1)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_1di4'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_1di4'

      end subroutine netcdf4restart_read_1di4

!.......................................................................

      subroutine netcdf4restart_wrt_5dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: vardim_,lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      integer ncid_,dim_5_id_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_4dr8'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_4dr8'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_4dr8'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_4dr8'
      idim_=ub(4)-lb(4)+1 ; status=nf_def_dim(ncid_,'nt',idim_,dim_t_id)
      if(status/=0) stop 'Err 5 netcdf4restart_wrt_4dr8'
      idim_=ub(5)-lb(5)+1 ; status=nf_def_dim(ncid_,'n5',idim_,dim_5_id_)
      if(status/=0) stop 'Err 5 netcdf4restart_wrt_4dr8'

      vardim_(1)=dim_x_id ; vardim_(2)=dim_y_id ; vardim_(3)=dim_z_id 
      vardim_(4)=dim_t_id ; vardim_(5)=dim_5_id_
      status=nf_def_var(ncid_,nom,nf_double,4,vardim_(1:5),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_4dr8'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_4dr8'

      status=nf_put_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_4dr8'


      end subroutine netcdf4restart_wrt_5dr8

!.......................................................................

      subroutine netcdf4restart_read_5dr8(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_4dr8'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_4dr8'

      status=nf_get_var_double(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_4dr8'

      end subroutine netcdf4restart_read_5dr8

!.......................................................................

      subroutine netcdf4restart_wrt_5dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: vardim_,lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      integer ncid_,dim_5_id_
      integer idim_

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         write(texte60,'(a,a,i0,a)')restartdir_out1//nom,'_',par%rank,'.nc'
        else
         write(texte60,'(a,a,i0,a)')restartdir_out2//nom,'_',par%rank,'.nc'
        endif

      status=nf_create(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_wrt_4dr8'

      idim_=ub(1)-lb(1)+1 ; status=nf_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdf4restart_wrt_4dr8'
      idim_=ub(2)-lb(2)+1 ; status=nf_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdf4restart_wrt_4dr8'
      idim_=ub(3)-lb(3)+1 ; status=nf_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdf4restart_wrt_4dr8'
      idim_=ub(4)-lb(4)+1 ; status=nf_def_dim(ncid_,'nt',idim_,dim_t_id)
      if(status/=0) stop 'Err 5 netcdf4restart_wrt_4dr8'
      idim_=ub(5)-lb(5)+1 ; status=nf_def_dim(ncid_,'n5',idim_,dim_5_id_)
      if(status/=0) stop 'Err 5 netcdf4restart_wrt_4dr8'

      vardim_(1)=dim_x_id ; vardim_(2)=dim_y_id ; vardim_(3)=dim_z_id 
      vardim_(4)=dim_t_id ; vardim_(5)=dim_5_id_
      status=nf_def_var(ncid_,nom,nf_real,4,vardim_(1:5),var_id)
      if(status/=0) stop 'Err 6 netcdf4restart_wrt_4dr8'

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_wrt_4dr8'

      status=nf_put_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
      if(status/=0)stop 'Err 8 netcdf4restart_wrt_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_wrt_4dr8'


      end subroutine netcdf4restart_wrt_5dr4

!.......................................................................

      subroutine netcdf4restart_read_5dr4(tab,nom,lb,ub)
      
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      integer ncid_
      integer idim_

      write(texte60,'(a,a,i0,a)') &
      'restart_input/'//nom,'_',par%rank,'.nc'

      status=nf_open(trim(texte60),nf_clobber,ncid_)
      if(status/=0) stop 'Err 1 netcdf4restart_read_4dr8'

      status=nf_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdf4restart_read_4dr8'

      status=nf_get_var_real(ncid_,var_id,tab(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
      if(status/=0)stop 'Err 8 netcdf4restart_read_4dr8'

      status=nf_close(ncid_)
      if(status/=0)stop 'Err 9 netcdf4restart_read_4dr8'

      end subroutine netcdf4restart_read_5dr4

!.......................................................................

      end module module_netcdf4restart
