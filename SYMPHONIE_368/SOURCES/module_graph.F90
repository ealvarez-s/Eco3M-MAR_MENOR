! --------------------------------------------------------
! Mdule de gestion du fichier texte decrivant les variables
! et routine d'ecriture
!--------------------------------------------------------

      module module_graph
      implicit none
       integer,parameter :: lunit=200
       integer,parameter :: maxvar=200
       integer,parameter :: nbgrid=4
       character(len=20) :: str20
       character(len=45) :: str45
                                !
!      -------------------------------------------------------------
!      Definit les numeros de type en fonction des compilateurs
!      Utile lors de la compilation avec -r8
!      TYPE real simple precision
       integer,parameter :: realsp=kind(1.0d00)/2
!      TYPE real double precision
       integer,parameter :: realdb=kind(1.0d00)
!      TYPE PAR DEFAUT POUR LES FOCNTIONS
       integer,parameter :: realdft=kind(1.0)
!      AMOD differ between g95 and f95
!      G95 AMOD(x,p) p is r4
!       integer,parameter :: AMODP=realSP
!      f95 AMOD(x,p) p is r4
       integer,parameter :: amodp=realdb
!      -------------------------------------------------------------

       type variable
       character(len=20) :: name
       character(len=20) :: unit
       character(len=45) :: stdname
       character(len=45) :: longname
       character(len=45) :: shortname
       integer           :: grid
       end type variable

       type(variable),dimension(maxvar) :: t_vars
       integer,dimension(maxvar) :: t_gridid
       integer,dimension(maxvar) :: t_varid
       character(len=60),dimension(maxvar) :: t_var_name
       character(len=60),dimension(nbgrid) :: t_grid_name
       integer :: nbvars

       contains

       subroutine readvarfile(file)
       implicit none
       character(len=*), intent(in) :: file
       integer ::        bcl,eof,grd
       bcl=0
       eof=0
       grd = -10
       !write(0,*) "FILE=",file
       open(lunit,file=file)
       do while(eof == 0)
          bcl=bcl+1
          !write(0,*) "bcl=",bcl
          read(lunit,*,iostat=eof) str20
          t_vars(bcl)%name = trim(str20) // char(0)
          !write(0,*) STR20
          read(lunit,*,iostat=eof) str20
          !write(0,*) STR20
          t_vars(bcl)%unit =  trim(str20) // char(0)
          read(lunit,*,iostat=eof) str45
          !write(0,*) STR45
          t_vars(bcl)%stdname = trim(str45) // char(0)
          read(lunit,*,iostat=eof) str45
          !write(0,*) STR45
          t_vars(bcl)%longname = trim(str45) // char(0)
          read(lunit,*,iostat=eof) str45
          !write(0,*) STR45
          t_vars(bcl)%shortname = trim(str45) // char(0)
          read(lunit,*,iostat=eof) grd
          !write(0,*) "grd=",grd
          t_vars(bcl)%grid = grd
          read(lunit,*,iostat=eof)
          !write(0,*) ""
          !read *
       enddo
       nbvars=bcl-1
       !write(0,*) "NBVARS=",NBVARS
       !write(0,*) T_VARS(1:NBVARS)%NAME
       !write(0,*) T_VARS(1:NBVARS)%UNIT
       close(lunit)
       end subroutine readvarfile
#ifdef synopsis
       subroutinetitle='readvarfile'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       integer function get_var_nb(name)
       implicit none
       character(len=*),intent(in) :: name
       integer bcl
       get_var_nb=0
       do bcl=1,nbvars
          if ( trim(t_vars(bcl)%name) == trim(name)) then
             get_var_nb=bcl
             return
          endif
       enddo
       if ((get_var_nb == nbvars) .or. (get_var_nb == 0)) then
          write(6,*) name," n a pas ete trouvee"
          stop
       endif
       end function get_var_nb

!     ---------------------------------------------------------
       subroutine poc_std_varxyt(file,filv,name,                        &
            scale,offset,v2d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,intent(in) :: filv,scale,offset
       real,dimension(:,:),intent(in) :: v2d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
#ifdef synopsis
       subroutinetitle='poc_std_varxyt'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       !write(0,*) NAME2," position=",NUM
       call pocstandardvariablexyt(                                     &
           file                                                         &
          ,real(filv,realsp)                                            &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,t_vars(num)%name                                             &
          ,t_vars(num)%unit                                             &
          ,real(scale,realsp)                                           &
          ,real(offset,realsp)                                          &
          ,t_vars(num)%stdname                                          &
          ,t_vars(num)%longname                                         &
          ,t_vars(num)%shortname                                        &
          ,status)

      call pocwritevar2d(                                               &
           file                                                         &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,ndim                                                         &
          ,tstart                                                       &
          ,real(v2d,realsp)                                             &
          ,status)

       end subroutine poc_std_varxyt
!     ---------------------------------------------------------

!     ---------------------------------------------------------
       subroutine poc_std_varxywt(file,filv,name,                       &
            scale,offset,v2d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,intent(in) :: filv,scale,offset
       real,dimension(:,:),intent(in) :: v2d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
#ifdef synopsis
       subroutinetitle='poc_std_varxywt'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       !write(0,*) NAME2," position=",NUM
       call pocstandardvariablexywt(                                    &
           file                                                         &
          ,real(filv,realsp)                                            &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,t_vars(num)%name                                             &
          ,t_vars(num)%unit                                             &
          ,real(scale,realsp)                                           &
          ,real(offset,realsp)                                          &
          ,t_vars(num)%stdname                                          &
          ,t_vars(num)%longname                                         &
          ,t_vars(num)%shortname                                        &
          ,status)

      call pocwritevar2d(                                               &
           file                                                         &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,ndim                                                         &
          ,tstart                                                       &
          ,real(v2d,realsp)                                             &
          ,status)

       end subroutine poc_std_varxywt
!     ---------------------------------------------------------

!     ---------------------------------------------------------
       subroutine poc_std_varxyzt(file,filv,name,                       &
            scale,offset,v3d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,intent(in) :: filv,scale,offset
       real,dimension(:,:,:),intent(in) :: v3d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
#ifdef synopsis
       subroutinetitle='poc_std_varxyzt'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       !write(0,*) NAME2," position=",NUM
       call pocstandardvariablexyzt(                                    &
           file                                                         &
          ,real(filv,realsp)                                            &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,t_vars(num)%name                                             &
          ,t_vars(num)%unit                                             &
          ,real(scale,realsp)                                           &
          ,real(offset,realsp)                                          &
          ,t_vars(num)%stdname                                          &
          ,t_vars(num)%longname                                         &
          ,t_vars(num)%shortname                                        &
          ,status)

      call pocwritevar3d(                                               &
           file                                                         &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,ndim                                                         &
          ,tstart                                                       &
          ,real(v3d,realsp)                                             &
          ,status)

       end subroutine poc_std_varxyzt
!     ---------------------------------------------------------

!     ---------------------------------------------------------
       subroutine poc_std_varxyzwt(file,filv,name,                      &
            scale,offset,v3d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,intent(in) :: filv,scale,offset
       real,dimension(:,:,:),intent(in) :: v3d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
#ifdef synopsis
       subroutinetitle='poc_std_varxyzwt'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       !write(0,*) NAME2," position=",NUM
       call pocstandardvariablexyzwt(                                   &
           file                                                         &
          ,real(filv,realsp)                                            &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,t_vars(num)%name                                             &
          ,t_vars(num)%unit                                             &
          ,real(scale,realsp)                                           &
          ,real(offset,realsp)                                          &
          ,t_vars(num)%stdname                                          &
          ,t_vars(num)%longname                                         &
          ,t_vars(num)%shortname                                        &
          ,status)

      call pocwritevar3d(                                               &
           file                                                         &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,ndim                                                         &
          ,tstart                                                       &
          ,real(v3d,realsp)                                             &
          ,status)

       end subroutine poc_std_varxyzwt
!     ---------------------------------------------------------

!     ---------------------------------------------------------
!         POC_STD_VAR???_C
!         POC_STD_VAR???_W
!     ---------------------------------------------------------
       subroutine poc_std_varxywt_c(file,filv,name,                     &
            scale,offset,v2d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,intent(in) :: filv,scale,offset
       real,dimension(:,:),intent(in) :: v2d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
#ifdef synopsis
       subroutinetitle='poc_std_varxywt_c'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       !write(0,*) NAME2," position=",NUM
       call pocstandardvariablexywt(                                    &
           file                                                         &
          ,real(filv,realsp)                                            &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,t_vars(num)%name                                             &
          ,t_vars(num)%unit                                             &
          ,real(scale,realsp)                                           &
          ,real(offset,realsp)                                          &
          ,t_vars(num)%stdname                                          &
          ,t_vars(num)%longname                                         &
          ,t_vars(num)%shortname                                        &
          ,status)

       end subroutine poc_std_varxywt_c

      subroutine  poc_std_2d_w(file,name,                               &
            v2d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,dimension(:,:),intent(in) :: v2d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
       !
       name2= trim(name) // char(0)
       num = get_var_nb(name2)
#ifdef synopsis
       subroutinetitle='poc_std_2d_w'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       call pocwritevar2d(                                              &
           file                                                         &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,ndim                                                         &
          ,tstart                                                       &
          ,real(v2d,realsp)                                             &
          ,status)

       end subroutine poc_std_2d_w
!     ---------------------------------------------------------
!     ---------------------------------------------------------
       subroutine poc_std_varxyzwt_c(file,filv,name,                    &
            scale,offset,v3d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,intent(in) :: filv,scale,offset
       real,dimension(:,:,:),intent(in) :: v3d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
#ifdef synopsis
       subroutinetitle='poc_std_varxyzwt_c'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       !write(0,*) NAME2," position=",NUM
       call pocstandardvariablexyzwt(                                   &
           file                                                         &
          ,real(filv,realsp)                                            &
          ,t_vars(num)%grid                                             &
          ,num                                                          &
          ,t_vars(num)%name                                             &
          ,t_vars(num)%unit                                             &
          ,real(scale,realsp)                                           &
          ,real(offset,realsp)                                          &
          ,t_vars(num)%stdname                                          &
          ,t_vars(num)%longname                                         &
          ,t_vars(num)%shortname                                        &
          ,status)

       end subroutine poc_std_varxyzwt_c

      subroutine  poc_std_3d_w(file,name,                               &
            v3d,ndim,tstart,status)
       implicit none
       character(len=*),intent(in) :: file,name
       real,dimension(:,:,:),intent(in) :: v3d
       integer,intent(in) :: ndim
       integer,dimension(:) :: tstart
       integer,intent(out) :: status
       !
       integer :: num
       character(len=20) :: name2
       !
       name2= trim(name) // char(0)
       num = get_var_nb(name2)
       call pocwritevar2d(                                              &
            file                                                        &
            ,t_vars(num)%grid                                           &
            ,num                                                        &
            ,ndim                                                       &
            ,tstart                                                     &
            ,real(v3d,realsp)                                           &
            ,status)
#ifdef synopsis
       subroutinetitle='poc_std_3d_w'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       end subroutine poc_std_3d_w
!     ---------------------------------------------------------

       end module module_graph

