      module module_netcdfrestart
!______________________________________________________________________
! SYMPHONIE Ocean Model
! release 366 - last update: 14-02-23
!______________________________________________________________________
! version   date    description                                               !
! S.26    14-12-15 mise en service                                            !
!         30-01-16 cas des grilles periodiques
!         04-02-16 permettre de repartir avec Gaspar (depuis K-eps)
!         02-04-16 Affichage ecran progression lecture/ecriture
!         13-09-16 flag_0status_option=1 ! 0=stop si variable non presente,  sinon ignore le fichier 
!         22-12-16 flag_0status_option modifiable depuis notebook_time
!         03-03-17 l'etiquette de compilation hotbin permet d'ecrire des fichiers
!                  restart "binaires" ancienne methode. L'interet est de pouvoir
!                  mieux traiter le debugage "checkmpi" car les fichieris restart 
!                  netcdf peuvent artificiellement retablir la continuite mpi 
!         14-03-17 L'etiquette de compilation hotbin est remplacee par flag_binary
!                  defini dans notebook_time.f
!         09-09-17 mask_u et cie integer kind=1
!         22-05-18 message pour aide au debugage
!         27-05-18 nfmpi_unlimited
! v260    05-10-19 ajout de nouvelles procedures
! v287    18-07-20 utiliser restartdir_out1 et cie...
! v295    29-12-20 utilisation de trim pour norme gfortran
! v366    14-02-23  flag_w_binary,flag_r_binary  
!______________________________________________________________________
!  _________                    .__                  .__             ! m°v°m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!        \/\/          \/|__|        \/            \/        \/      !
!______________________________________________________________________

      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
!     include 'netcdf.inc'
!     integer :: flag_netcdf_def=0
      integer(kind=MPI_OFFSET_KIND) start(5)
      integer(kind=MPI_OFFSET_KIND) edge(5)
      integer ncid_,shift_i_,shift_j_
      integer(kind=MPI_OFFSET_KIND)idim_,imax_,jmax_,kmax_
!     integer dynrstdim
      integer flag_header_   ! flag_header_=1 signifie fichier deja existant, 0 sinon.

! Variable desormais dans module_principal.F90 et modifiable depuis notebook_time !22-12-16
!     integer :: flag_0status_option=1 ! 0=stop si variable non presente,  sinon ignore le fichier !13-09-16

      interface netcdfrestart_wrt
      module procedure netcdfrestart_wrt_5dr8
      module procedure netcdfrestart_wrt_5dr4
      module procedure netcdfrestart_wrt_4dr8
      module procedure netcdfrestart_wrt_4dr4
      module procedure netcdfrestart_wrt_4di4
      module procedure netcdfrestart_wrt_3dr8
      module procedure netcdfrestart_wrt_3dr4
      module procedure netcdfrestart_wrt_3di4
      module procedure netcdfrestart_wrt_3di2
      module procedure netcdfrestart_wrt_2dr8
      module procedure netcdfrestart_wrt_2dr4
      module procedure netcdfrestart_wrt_2di4
      module procedure netcdfrestart_wrt_2di1
      end interface 

      interface netcdfrestart_read
      module procedure netcdfrestart_read_5dr8
      module procedure netcdfrestart_read_5dr4
      module procedure netcdfrestart_read_4dr8
      module procedure netcdfrestart_read_4dr4
      module procedure netcdfrestart_read_4di4
      module procedure netcdfrestart_read_3dr8
      module procedure netcdfrestart_read_3dr4
      module procedure netcdfrestart_read_3di4
      module procedure netcdfrestart_read_3di2
      module procedure netcdfrestart_read_2dr8
      module procedure netcdfrestart_read_2dr4
      module procedure netcdfrestart_read_2di4
      module procedure netcdfrestart_read_2di1
      end interface 

      contains

!.......................................................................

      subroutine netcdfrestart_bounds_wrt(lb1_,ub1_,lb2_,ub2_)
      implicit none
      integer lb1_,ub1_,lb2_,ub2_

       istr=2 ; iend=imax-1 ; jstr=2 ; jend=jmax-1

       if(par%tvoisin(ieq1   )==mpi_proc_null)istr=lb1_
       if(par%tvoisin(ieqimax)==mpi_proc_null)iend=ub1_
       if(par%tvoisin(jeq1   )==mpi_proc_null)jstr=lb2_
       if(par%tvoisin(jeqjmax)==mpi_proc_null)jend=ub2_

! Je garde la condition precedente dont je ne sais pas si elle ne serait pas
! indispensabe en cas de proc inutiles supprimEs... mais cette condition
! est prise en defaut en cas de grille periodique donc j'ajoute la condition
! suivante pour blinder le cas des grilles periodiques:
      if(1   +par%timax(1)==1)   istr=lb1_ !30-01-16
      if(imax+par%timax(1)==iglb)iend=ub1_
      if(1   +par%tjmax(1)==1)   jstr=lb2_ 
      if(jmax+par%tjmax(1)==jglb)jend=ub2_

! ATTENTION LES BORDS DU DOMAINE COMMENCENT A 1
! Par consequent on ajoute shift_i_ et shift_j_  A start(1) et start(2)
      shift_i_=1-lb1_
      shift_j_=1-lb2_

      start(1)=istr+par%timax(1)+shift_i_
      start(2)=jstr+par%tjmax(1)+shift_j_
       edge(1)=iend-istr+1      
       edge(2)=jend-jstr+1

      end subroutine netcdfrestart_bounds_wrt

!.......................................................................

      subroutine netcdfrestart_bounds_read(lb1_,ub1_,lb2_,ub2_)
      implicit none
      integer lb1_,ub1_,lb2_,ub2_

! ATTENTION LES BORDS DU DOMAINE COMMENCENT A 1
! Par consequent on ajoute shift_i_ et shift_j_  A start(1) et start(2)
      shift_i_=1-lb1_
      shift_j_=1-lb2_

      istr=lb1_ ; iend=ub1_
      jstr=lb2_ ; jend=ub2_

      start(1)=istr+par%timax(1)+shift_i_
      start(2)=jstr+par%tjmax(1)+shift_j_
       edge(1)=iend-istr+1
       edge(2)=jend-jstr+1

      end subroutine netcdfrestart_bounds_read

!.......................................................................

      subroutine netcdfrestart_inq_headfile
      implicit none

! Tester l'existence du fichier
      status=nfmpi_open(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status==0) then !xxx>
       flag_header_=1 ! flag_header_=1 signifie fichier deja existant
       status=nfmpi_close(ncid_)
      else               !xxx>
       flag_header_=0 ! flag_header_=0 signifie fichier inexistant
      endif              !xxx>

      end subroutine netcdfrestart_inq_headfile

!.......................................................................

      subroutine netcdfrestart_wrt_5dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab

      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      !dynrstdim=5
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1
      bstr=lb(5) ; bend=ub(5)
      start(5)=1 ; edge(5)=bend-bstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)

      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_5dr8'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_5dr8'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_5dr8'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_5dr8'

      idim_=ub(4)-lb(4)+1 
      status=nfmpi_def_dim(ncid_,'nt',idim_,dim_t_id)
      if(status/=0) stop 'Err 5 netcdfrestart_wrt_5dr8'

      idim_=ub(5)-lb(5)+1 
!     status=nfmpi_def_dim(ncid_,'nb',idim_,dim_b_id)
      status=nfmpi_def_dim(ncid_,'nb',nfmpi_unlimited,dim_b_id) !27-05-18
      if(status/=0) stop 'Err 5b netcdfrestart_wrt_5dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id 
      vardim(4)=dim_t_id ; vardim(5)=dim_b_id

      status=nfmpi_def_var(ncid_,nom,nf_double,5,vardim(1:5),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_5dr8'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_5dr8'

      status=nfmpi_put_vara_double_all(ncid_,var_id           &
                                            ,start(1:5)       &
                                             ,edge(1:5)       &
      ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend,bstr:bend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_5dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_5dr8'


      end subroutine netcdfrestart_wrt_5dr8

!.......................................................................

      subroutine netcdfrestart_read_5dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1
      bstr=lb(5) ; bend=ub(5)
      start(5)=1 ; edge(5)=bend-bstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_5dr8'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_5dr8'

      status=nfmpi_get_vara_double_all(ncid_,var_id      &
                                            ,start(1:5)       &
                                             ,edge(1:5)       &
      ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend,bstr:bend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_5dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_5dr8'

      end subroutine netcdfrestart_read_5dr8

!.......................................................................

      subroutine netcdfrestart_wrt_5dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=5
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1
      bstr=lb(5) ; bend=ub(5)
      start(5)=1 ; edge(5)=bend-bstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_5dr4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_5dr4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_5dr4'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_5dr4'

      idim_=ub(4)-lb(4)+1 
      status=nfmpi_def_dim(ncid_,'nt',idim_,dim_t_id)
      if(status/=0) stop 'Err 5 netcdfrestart_wrt_5dr4'

      idim_=ub(5)-lb(5)+1 
!     status=nfmpi_def_dim(ncid_,'nb',idim_,dim_b_id)
      status=nfmpi_def_dim(ncid_,'nb',nfmpi_unlimited,dim_b_id) !27-05-18
      if(status/=0) stop 'Err 5b netcdfrestart_wrt_5dr4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id 
      vardim(4)=dim_t_id ; vardim(5)=dim_b_id

      status=nfmpi_def_var(ncid_,nom,nf_real,5,vardim(1:5),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_5dr4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'
      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_5dr4'

      status=nfmpi_put_vara_real_all(ncid_,var_id           &
                                            ,start(1:5)       &
                                             ,edge(1:5)       &
      ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend,bstr:bend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_5dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_5dr4'


      end subroutine netcdfrestart_wrt_5dr4

!.......................................................................

      subroutine netcdfrestart_read_5dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(5) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1
      bstr=lb(5) ; bend=ub(5)
      start(5)=1 ; edge(5)=bend-bstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_5dr4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_5dr4'

      status=nfmpi_get_vara_real_all(ncid_,var_id      &
                                            ,start(1:5)       &
                                             ,edge(1:5)       &
      ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend,bstr:bend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_5dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_5dr4'

      end subroutine netcdfrestart_read_5dr4
!.......................................................................

      subroutine netcdfrestart_wrt_4dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=4
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_4dr8'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_4dr8'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_4dr8'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4dr8'

      idim_=ub(4)-lb(4)+1 
!!!!!!status=nfmpi_def_dim(ncid_,'nt',idim_,dim_t_id)
      status=nfmpi_def_dim(ncid_,'nt',nfmpi_unlimited,dim_t_id) !27-05-18
      if(status/=0) stop 'Err 5 netcdfrestart_wrt_4dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id ; vardim(4)=dim_t_id
      status=nfmpi_def_var(ncid_,nom,nf_double,4,vardim(1:4),var_id)
      call netcdfrestart_error('netcdfrestart_wrt_4dr8','nfmpi_def_var',nom,4)

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_4dr8'

       status=nfmpi_put_vara_double_all(ncid_,var_id      &
                                             ,start(1:4)  &
                                              ,edge(1:4)  &
           ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
       if(status/=0)stop 'Err 8 netcdfrestart_wrt_4dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_4dr8'


      end subroutine netcdfrestart_wrt_4dr8

!.......................................................................

      subroutine netcdfrestart_read_4dr8(tab,nom,lb,ub)
!      use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_4dr8'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_4dr8'

      status=nfmpi_get_vara_double_all(ncid_,var_id      &
                                            ,start(1:4)  &
                                             ,edge(1:4)  &
          ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_4dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_4dr8'

      end subroutine netcdfrestart_read_4dr8

!.......................................................................

      subroutine netcdfrestart_wrt_4dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=4
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_4dr4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_4dr4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_4dr4'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4dr4'

      idim_=ub(4)-lb(4)+1 
!     status=nfmpi_def_dim(ncid_,'nt',idim_,dim_t_id)
      status=nfmpi_def_dim(ncid_,'nt',nfmpi_unlimited,dim_t_id) !27-05-18
      if(status/=0) stop 'Err 5 netcdfrestart_wrt_4dr4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id ; vardim(4)=dim_t_id
      status=nfmpi_def_var(ncid_,nom,nf_real,4,vardim(1:4),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_4dr4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_4dr4'

      status=nfmpi_put_vara_real_all(ncid_,var_id      &
                                            ,start(1:4)  &
                                             ,edge(1:4)  &
          ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_4dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_4dr4'


      end subroutine netcdfrestart_wrt_4dr4

!.......................................................................

      subroutine netcdfrestart_read_4dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_4dr4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_4dr4'

      status=nfmpi_get_vara_real_all(ncid_,var_id      &
                                            ,start(1:4)  &
                                             ,edge(1:4)  &
          ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_4dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_4dr4'

      end subroutine netcdfrestart_read_4dr4

!.......................................................................

      subroutine netcdfrestart_wrt_4di4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=4
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_4di4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_4di4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_4di4'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4di4'

      idim_=ub(4)-lb(4)+1 
!     status=nfmpi_def_dim(ncid_,'nt',idim_,dim_t_id)
      status=nfmpi_def_dim(ncid_,'nt',nfmpi_unlimited,dim_t_id) !27-05-18
      if(status/=0) stop 'Err 5 netcdfrestart_wrt_4di4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id ; vardim(4)=dim_t_id
      status=nfmpi_def_var(ncid_,nom,nf_real,4,vardim(1:4),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_4di4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_4di4'

      status=nfmpi_put_vara_int_all(ncid_,var_id      &
                                            ,start(1:4)  &
                                             ,edge(1:4)  &
          ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_4di4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_4di4'


      end subroutine netcdfrestart_wrt_4di4

!.......................................................................

      subroutine netcdfrestart_read_4di4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(4) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1
      tstr=lb(4) ; tend=ub(4)
      start(4)=1 ; edge(4)=tend-tstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_4di4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_4di4'

      status=nfmpi_get_vara_int_all(ncid_,var_id      &
                                            ,start(1:4)  &
                                             ,edge(1:4)  &
          ,tab(istr:iend,jstr:jend,kstr:kend,tstr:tend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_4di4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_4di4'

      end subroutine netcdfrestart_read_4di4

!.......................................................................

      subroutine netcdfrestart_wrt_3dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=3
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_3dr8'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_3dr8'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_3dr8'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id
      status=nfmpi_def_var(ncid_,nom,nf_double,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_3dr8'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_3dr8'

      status=nfmpi_put_vara_double_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))

      call netcdfrestart_error('netcdfrestart_wrt_3dr8','nfmpi_put',nom,3)

      status=nfmpi_close(ncid_)

      end subroutine netcdfrestart_wrt_3dr8

!.......................................................................

      subroutine netcdfrestart_read_3dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:
      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)


      if(status/=0)then
         if(flag_0status_option==1) return !13-09-16
         if(nom=='tkle_w')return !04-02-16
         if(nom=='tkll_w')return !04-02-16
         write(6,'(a,a)')'Err reading restart variable ',nom
         stop 'Err 1 netcdfrestart_read_3dr8'
      endif

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_3dr8'

      status=nfmpi_get_vara_double_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_3dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_3dr8'

      end subroutine netcdfrestart_read_3dr8

!.......................................................................

      subroutine netcdfrestart_wrt_3dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=3
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_3dr4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_3dr4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_3dr4'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id
      status=nfmpi_def_var(ncid_,nom,nf_real,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_3dr4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_3dr4'

      status=nfmpi_put_vara_real_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_3dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_3dr4'


      end subroutine netcdfrestart_wrt_3dr4

!.......................................................................

      subroutine netcdfrestart_read_3dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_3dr4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_3dr4'

      status=nfmpi_get_vara_real_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_3dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_3dr4'

      end subroutine netcdfrestart_read_3dr4

!.......................................................................

      subroutine netcdfrestart_wrt_3di4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=3
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_3di4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_3di4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_3di4'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id
      status=nfmpi_def_var(ncid_,nom,nf_int,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_3di4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_3di4'

      status=nfmpi_put_vara_int_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_3di4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_3di4'


      end subroutine netcdfrestart_wrt_3di4

!.......................................................................

      subroutine netcdfrestart_read_3di4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)



      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_3di4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_3di4'

      status=nfmpi_get_vara_int_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_3di4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_3di4'

      end subroutine netcdfrestart_read_3di4

!.......................................................................

      subroutine netcdfrestart_wrt_3di2(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      integer(kind=1),dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=3
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_3di2'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_3di2'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_3di2'

      idim_=ub(3)-lb(3)+1 
      status=nfmpi_def_dim(ncid_,'nk',idim_,dim_z_id)
      if(status/=0) stop 'Err 4 netcdfrestart_wrt_4dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id ; vardim(3)=dim_z_id
      status=nfmpi_def_var(ncid_,nom,nf_byte,3,vardim(1:3),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_3di2'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_3di2'

      status=nfmpi_put_vara_int1_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_3di2'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_3di2'


      end subroutine netcdfrestart_wrt_3di2

!.......................................................................

      subroutine netcdfrestart_read_3di2(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(3) :: lb,ub
      integer(kind=1),dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      kstr=lb(3) ; kend=ub(3)
      start(3)=1 ; edge(3)=kend-kstr+1

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)



      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_3di2'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_3di2'

      status=nfmpi_get_vara_int1_all(ncid_,var_id      &
                                            ,start(1:3)  &
                                             ,edge(1:3)  &
                    ,tab(istr:iend,jstr:jend,kstr:kend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_3di2'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_3di2'

      end subroutine netcdfrestart_read_3di2

!.......................................................................

      subroutine netcdfrestart_wrt_2dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=2
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)

      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)

      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_2dr8'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_2dr8'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_2dr8'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nfmpi_def_var(ncid_,nom,nf_double,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_2dr8'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_2dr8'

      status=nfmpi_put_vara_double_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_2dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_2dr8'


      end subroutine netcdfrestart_wrt_2dr8

!.......................................................................

      subroutine netcdfrestart_read_2dr8(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

      !dynrstdim=2
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_2dr8'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_2dr8'

      status=nfmpi_get_vara_double_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_2dr8'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_2dr8'

      end subroutine netcdfrestart_read_2dr8

!.......................................................................

      subroutine netcdfrestart_wrt_2dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=2
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)

      if(flag_header_==0) then !00000>
      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_2dr4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_2dr4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_2dr4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nfmpi_def_var(ncid_,nom,nf_real,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_2dr4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_2dr4'

      status=nfmpi_put_vara_real_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_2dr4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_2dr4'


      end subroutine netcdfrestart_wrt_2dr4

!.......................................................................

      subroutine netcdfrestart_read_2dr4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      real,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0)stop 'Err 1 netcdfrestart_read_2dr4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_2dr4'

      status=nfmpi_get_vara_real_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_2dr4'


      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_2dr4'

      end subroutine netcdfrestart_read_2dr4

!.......................................................................

      subroutine netcdfrestart_wrt_2di4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=2
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_2di4'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_2di4'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_2di4'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nfmpi_def_var(ncid_,nom,nf_int,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_2di4'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_2di4'

      status=nfmpi_put_vara_int_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_2di4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_2di4'


      end subroutine netcdfrestart_wrt_2di4

!.......................................................................

      subroutine netcdfrestart_read_2di4(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_2di4'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_2di4'

      status=nfmpi_get_vara_int_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_2di4'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_2di4'

      end subroutine netcdfrestart_read_2di4

!.......................................................................

      subroutine netcdfrestart_wrt_2di1(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      integer(kind=1),dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_w_binary==1) then ; write(9)tab ; return ; endif !04-03-17

      !dynrstdim=2
! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
          write(9)tab !call netcdf4restart_wrt(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

      call netcdfrestart_bounds_wrt(lb(1),ub(1),lb(2),ub(2)) ! computes istr iend jstr jend start(1) start(2) edge(1) edge(2)

        if(index(dynrestartfilename,restartdir_out1)/=0) then
         !write(texte60,'(a,a)')restartdir_out1//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out1),trim(nom),'.nc' !29-12-20
        else
         !write(texte60,'(a,a)')restartdir_out2//nom,'.nc'
         write(texte60,'(a,a,a)')trim(restartdir_out2),trim(nom),'.nc' !29-12-20
        endif
        texte60=trim(texte60)
      call netcdfrestart_inq_headfile ! computes flag_header_ (does file already exist?)
      if(flag_header_==0) then !00000>

      status=nfmpi_create(par%comm2d,texte60,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0) stop 'Err 1 netcdfrestart_wrt_2di1'

      idim_= (iglb+(ub(1)-imax)) -(1+(lb(1)-1))+1
      status=nfmpi_def_dim(ncid_,'ni',idim_,dim_x_id)
      if(status/=0) stop 'Err 2 netcdfrestart_wrt_2di1'

      idim_= (jglb+(ub(2)-jmax)) -(1+(lb(2)-1))+1
      status=nfmpi_def_dim(ncid_,'nj',idim_,dim_y_id)
      if(status/=0) stop 'Err 3 netcdfrestart_wrt_2di1'

      vardim(1)=dim_x_id ; vardim(2)=dim_y_id 
      status=nfmpi_def_var(ncid_,nom,nf_int,2,vardim(1:2),var_id)
      if(status/=0) stop 'Err 6 netcdfrestart_wrt_2di1'

      status=nfmpi_enddef(ncid_)
      status=nfmpi_close(ncid_)

      endif                    !00000>

      status=nfmpi_open(par%comm2d,texte60,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)
      if(status/=0)STOP 'Err 29 nfmpi_open'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_wrt_2di1'

      status=nfmpi_put_vara_int1_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_wrt_2di1'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_wrt_2di1'


      end subroutine netcdfrestart_wrt_2di1

!.......................................................................

      subroutine netcdfrestart_read_2di1(tab,nom,lb,ub)
      !use module_netcdf4restart
      implicit none
      character(len=*) :: nom
      integer,dimension(2) :: lb,ub
      integer(kind=1),dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
      if(par%rank==0)write(6,'(a)')trim(nom) !01-04-16
      if(flag_r_binary==1) then ; read(9)tab ; return ; endif !04-03-17

! Si il ne s'agit pas d'un tableau dont les 2 premieres dimensions sont X,Y alors
! utiliser netcdf 4 ou chanel9 au choix:
      if(lb(1)>2.or.ub(1)<imax-1.or.lb(2)>2.or.ub(2)<jmax-1) then
         read(9)tab !call netcdf4restart_read(tab,nom,lb,ub) 
         return
      endif 
! Sinon utiliser netcdf parallele:

! computes istr iend jstr jend start(1) start(2) edge(1) edge(2) for "read" case:
      call netcdfrestart_bounds_read(lb(1),ub(1),lb(2),ub(2)) 

      write(texte60,'(a,a)')'restart_input/'//nom,'.nc'

      status=nfmpi_open(par%comm2d,texte60,nf_nowrite+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      if(status/=0.and.flag_0status_option==1) return !13-09-16
      if(status/=0) stop 'Err 1 netcdfrestart_read_2di1'

      status=nfmpi_inq_varid(ncid_,nom,var_id)
      if(status/=0)stop 'Err 7 netcdfrestart_read_2di1'

      status=nfmpi_get_vara_int1_all(ncid_,var_id      &
                                            ,start(1:2)  &
                                             ,edge(1:2)  &
                          ,tab(istr:iend,jstr:jend))
      if(status/=0)stop 'Err 8 netcdfrestart_read_2di1'

      status=nfmpi_close(ncid_)
      if(status/=0)stop 'Err 9 netcdfrestart_read_2di1'

      end subroutine netcdfrestart_read_2di1

!.......................................................................

      subroutine netcdfrestart_error(routinename_,functioname_,nom,dimcase_)
      implicit none
      character(len=*) routinename_,functioname_,nom
      integer dimcase_

      flag_stop=0
      if(status/=0) then !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS>
      flag_stop=1

        write(10+par%rank,'(a)') 'Error in module_netcdfrestart.F90'
        write(10+par%rank,*)'status= ',status
        write(10+par%rank,'(2a)') 'Subroutine name: ',routinename_
        write(10+par%rank,'(2a)') 'Function       : ',functioname_ 

      if(index(functioname_,'nfmpi_put')/=0) then !ooooo>

        status=nfmpi_inq_dimid(ncid_,'ni',dim_x_id)
        status=nfmpi_inq_dimlen(ncid_,dim_x_id,imax_)
        status=nfmpi_inq_dimid(ncid_,'nj',dim_y_id)
        status=nfmpi_inq_dimlen(ncid_,dim_y_id,jmax_)
        write(10+par%rank,'(2a)') 'Field name     : ',nom
        write(10+par%rank,*)'start(1:3)',start(1:3)
        write(10+par%rank,*)'edge(1:3)',edge(1:3)
        write(10+par%rank,*)'istr,iend',istr,iend
        write(10+par%rank,*)'jstr,jend',jstr,jend
        write(10+par%rank,*)'ni=',imax_
        write(10+par%rank,*)'nj=',jmax_

      if(dimcase_==3) then !33333333333>
        status=nfmpi_inq_dimid(ncid_,'nk',dim_z_id)
        status=nfmpi_inq_dimlen(ncid_,dim_z_id,kmax_)
        write(10+par%rank,*)'kstr,kend',kstr,kend
        write(10+par%rank,*)'nk=',kmax_
      endif                !33333333333>

      endif                                       !ooooo>

      if(index(functioname_,'nfmpi_def_var')/=0) then !xxxxx>
        write(10+par%rank,'(2a)') 'Field name     : ',nom
        write(10+par%rank,*)'vardim ',vardim(1:dimcase_)
        write(10+par%rank,*)'var_id',var_id
      endif                                           !xxxxx>

      endif              !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS>
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop 'Err netcdfrestart_error see fort... files'

      end subroutine netcdfrestart_error

!.......................................................................

      end module module_netcdfrestart
