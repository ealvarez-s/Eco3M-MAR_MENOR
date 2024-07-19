      module module_s
!______________________________________________________________________
! SYMPHONIE ocean model
! release 302  - last update: 20-05-21
!______________________________________________________________________
! version   date    description                                                !
! S26     05-08-14  mise en service                                            !
!         25-05-15  ajout routine cpu_s
!         03-07-15  ajout de s_matrick_reverser (inversion de matrice)
!         13-10-15  ajout subroutine s_smooth_xy_t
!         20-02-16  Afficher les cpu A l'iteration 0
!         29-05-17  Ajout subroutine s_stop_by_interrupteur_file
!         08-06-17  ajout message debug
! v250    18-03-19  ajout s_char2int_convert
! v272    14-01-20  interface renvoyant le type (r4 ou r8) d'une variable
! v274    10-02-20  subroutine donne_le_type_0di1(var_) !10-02-20
! v296    22-02-21  subroutine s_polygone !22-02-21
! v302    20-05-21  adaptation a la norme gfortran qui refuse .true. et .false. 
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal ; use module_parallele ; use module_parameter
      implicit none
      double precision scpu1,scpu2
      integer :: scpuk=0
      integer :: matdim=1+dim_ref/2
      integer(kind=1) :: flag_s_stop=0 &
                        ,s_typevar=0

! useful "s" functions 

      interface donne_le_type !14-01-20
       module procedure donne_le_type_0dr4
       module procedure donne_le_type_0dr8
       module procedure donne_le_type_0di1 !10-02-20
      end interface

contains

!......................................................................

      integer function s_unit(unit_)
!     use module_principal ; use module_parallele
      implicit none
      integer :: unit_       
      logical :: openedunit_

!      subroutinedescription=                                     &
!         'Inquires if submitted file unit is free and, if not,'  &
!      //' returns a new free file unit. Note: the returned unit' &
!      //' is >=7 what ever the received unit'


! Les lignes suivantes recherchent par ordre croissant une valeur de "unit" libre pour
! la fontion open(unit=  Pour eviter les unites speciales 5 et 6, on demarre a 7
      s_unit=max(unit_,7)            
      inquire(unit=s_unit,opened=openedunit_)
      do while (openedunit_)
       s_unit=s_unit+1
       inquire(unit=s_unit,opened=openedunit_)
      enddo

      end function s_unit

!......................................................................

      subroutine s_cpu(txt_,case_)
!     use module_principal ; use module_parallele
      implicit none
      character (len=*) txt_
      integer case_
#ifdef parallele

      if(case_==0.and.mod(iteration3d,100)/=0)return !20-02-16

      scpu1=scpu2
      scpu2=MPI_Wtime ( ) ! reset seconds_cpu_

      if(scpuk>0) then !------->
       if(par%rank==0)write(6,'(a,a20,f10.5)')'cpu seconds step ',txt_,scpu2-scpu1
      endif            !------->

      scpuk=scpuk+1

#endif
      end subroutine s_cpu

!......................................................................
!     subroutine s_matrix_reverser(matrix,inverse,errorflag) !03-06-15
      subroutine s_matrix_reverser(matrix,inverse,errorflag,lb_,ub_)
      implicit none
      integer,dimension(2) :: lb_,ub_
      integer :: errorflag,i_,j_,k_
!     real*8, dimension(matdim,matdim)   :: matrix,inverse
!     real*8, dimension(matdim,matdim*2) :: augmat
      real*8, dimension(lb_(1):ub_(1),lb_(2):ub_(2))   :: matrix,inverse
      real*8, dimension(lb_(1):ub_(1),lb_(2):ub_(2)*2) :: augmat
      real*8 x0_
      logical :: flag = .true.

      matdim=ub_(1)

! En entree: matrix 
! En sortie: "inverse" est la matrice inverse de "matrix"

! Augmentation de la matrice:
      do i_=1,matdim
      do j_=1,2*matdim
      if (j_<=matdim) then
      augmat(i_,j_) = matrix(i_,j_)
      else if ((i_+matdim) == j_) then
      augmat(i_,j_) = 1
      else
      augmat(i_,j_) = 0
      endif
      end do
      end do

! Triangularisation:
      do k_=1,matdim-1
      if (augmat(k_,k_) == 0) then
      flag = .false.
      do i_ = k_+1, matdim
      if (augmat(i_,k_) /= 0) then
      do j_ = 1,2*matdim
      augmat(k_,j_) = augmat(k_,j_)+augmat(i_,j_)
      end do
      flag = .true.
      exit
      endif

      if (flag .eqv. .false.) then
      write(6,*)'Err2 matrice ne peut etre inversee'
      inverse=0
      errorflag=-1
      return
      endif
      end do
      endif

      do j_ = k_+1, matdim
      x0_ = augmat(j_,k_)/augmat(k_,k_)
      do i_ = k_, 2*matdim
      augmat(j_,i_)=augmat(j_,i_)-x0_*augmat(k_,i_)
      end do
      end do
      end do

! Inversion reussie?
      do i_=1,matdim
      if (augmat(i_,i_) == 0) then
       inverse=0
       errorflag=-1
       write(6,*)'i,j,k',i,j,k
       stop 'Err1 matrice ne peut etre inversee'
      endif
      end do

! Normalisation unitaire de la diagonale principale:
      do i_=1,matdim
      x0_=augmat(i_,i_)
      do j_=i_,(2*matdim)
      augmat(i_,j_)=(augmat(i_,j_)/x0_)
      end do
      end do

! Partie droite de la matrice augmentee reduite a la matrice identite:
      do k_ =matdim-1,1,-1
      do i_=1,k_
      x0_=augmat(i_,k_+1)
      do j_=k_,(2*matdim)
      augmat(i_,j_)=augmat(i_,j_)-augmat(k_+1,j_)*x0_
      end do
      end do
      end do

! Finalisation du calcul:
      do i_ =1, matdim
      do j_ = 1, matdim
      inverse(i_,j_) = augmat(i_,j_+matdim)
      end do
      end do
      errorflag = 0

      end subroutine s_matrix_reverser

!......................................................................
#ifdef bidon
      subroutine s_smooth_xy_t(loopmax_) !13-10-15
      integer istr_,iend_,jstr_,jend_,loop_,loopmax_

! 5 points filter with impermeable solid boundaries

      do loop_=1,loopmax_

       do j=2,jmax-1 ; do i=2,imax-1
            xy_t(i  ,j  ,id_aft)=xy_t(i,j,id_now)                       &
         +((xy_t(i+1,j  ,id_now)-xy_t(i,j,id_now))*mask_u(i+1,j  ,kmax) &
          +(xy_t(i-1,j  ,id_now)-xy_t(i,j,id_now))*mask_u(i  ,j  ,kmax) &
          +(xy_t(i  ,j+1,id_now)-xy_t(i,j,id_now))*mask_v(i  ,j+1,kmax) &
          +(xy_t(i  ,j-1,id_now)-xy_t(i,j,id_now))*mask_v(i  ,j  ,kmax) &
          )*0.25
       enddo            ; enddo   

       do j=2,jmax-1 ; do i=2,imax-1
        xy_t(i,j,id_now)=xy_t(i,j,id_aft)
       enddo            ; enddo   

! 'z0' Open boundary conditions:

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then   !-----------> 
       do j=2,jmax-1
        xy_t(imax,j,id_now)=xy_t(imax-1,j,id_now)
       enddo
      endif                            !----------->

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then      !----------->
       do j=2,jmax-1
        xy_t(1,j,id_now)=xy_t(2,j,id_now)
       enddo
      endif                            !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then   !----------->
       do i=1,imax
        xy_t(i,jmax,id_now)=xy_t(i,jmax-1,id_now)
       enddo
      endif                            !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then      !----------->
       do i=1,imax
        xy_t(i,1,id_now)=xy_t(i,2,id_now)
       enddo
      endif                            !----------->

! 'z1' mpi continuity
       call obc_ext_xy_t('z0',id_now)      


      enddo ! loop_
      

      end subroutine s_smooth_xy_t
#endif
!......................................................................
      subroutine s_stop_by_interrupteur_file(message_txt_)
      implicit none
      character(len=*) message_txt_

! Programmer un arret de la simulation A la fin du cycle iteratif
! en passant par le fichier "interrupteur"
! Le motif de l'arret est ecrit dans message_txt_ et il est
! archivE dans un fichier fort***

       open(unit=3,file=trim(tmpdirname)//'interrupteur')
        write(3,'(i1,9x,i1,9x,i1)')1,1,0
       close(3)

       write(10+par%rank,'(a)') &
        'stop request in "interrupteur" file" with following message:'
       write(10+par%rank,'(a)')message_txt_
       write(10+par%rank,'(a,i0)')'iteration3d ',iteration3d !08-06-17

      end subroutine s_stop_by_interrupteur_file

!......................................................................

      subroutine s_stop(txt_,case_)
      implicit none
      character(len=*) :: txt_
      integer :: case_,flag_stop_glb_=0

      if(case_==1) then !case1>
       flag_stop=1
       write(10+par%rank,'(a)')'.......................'
       write(10+par%rank,'(a)')trim(txt_)
      endif             !case1>

      if(case_==2) then !case2>
       call mpi_allreduce(flag_stop,flag_stop_glb_,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       if(flag_stop_glb_/=0) stop 'Error see fort.xx files'
      endif             !case2>


      end subroutine s_stop

!...................................................................

      subroutine s_char2int_convert(txt_,year_,month_,day_,hour_,minute_,second_) !18-03-19
      implicit none
      character(len=*) :: txt_
      integer year_,month_,day_,hour_,minute_,second_  &
             ,intconverted_,factor_,ichr_,k1_,k2_,k_
      integer(kind=1) loop_

! En entree une chaine de characteres au format yyyymmdd-hhmmss (exemple txt_='20091122-035012')
! En sortie year_,month_,day_,hour_,minute_,second_ 

      do loop_=1,6

      if(loop_==1) then ; k1_=1  ; k2_=4  ; endif
      if(loop_==2) then ; k1_=5  ; k2_=6  ; endif
      if(loop_==3) then ; k1_=7  ; k2_=8  ; endif
      if(loop_==4) then ; k1_=10 ; k2_=11 ; endif
      if(loop_==5) then ; k1_=12 ; k2_=13 ; endif
      if(loop_==6) then ; k1_=14 ; k2_=15 ; endif

      factor_=1
      intconverted_=0
      do k_=k2_,k1_,-1
       if(txt_(k_:k_)=='0')ichr_=0
       if(txt_(k_:k_)=='1')ichr_=1
       if(txt_(k_:k_)=='2')ichr_=2
       if(txt_(k_:k_)=='3')ichr_=3
       if(txt_(k_:k_)=='4')ichr_=4
       if(txt_(k_:k_)=='5')ichr_=5
       if(txt_(k_:k_)=='6')ichr_=6
       if(txt_(k_:k_)=='7')ichr_=7
       if(txt_(k_:k_)=='8')ichr_=8
       if(txt_(k_:k_)=='9')ichr_=9
       intconverted_=intconverted_+factor_*ichr_
       factor_=factor_*10
      enddo

      if(loop_==1)year_  =intconverted_ 
      if(loop_==2)month_ =intconverted_ 
      if(loop_==3)day_   =intconverted_ 
      if(loop_==4)hour_  =intconverted_ 
      if(loop_==5)minute_=intconverted_ 
      if(loop_==6)second_=intconverted_ 

      enddo ! loop_

      end subroutine s_char2int_convert

!...................................................................

       subroutine donne_le_type_0dr4(var_) !14-01-20
       real var_
       s_typevar=4
       end subroutine donne_le_type_0dr4

       subroutine donne_le_type_0dr8(var_)
       double precision var_
       s_typevar=8
       end subroutine donne_le_type_0dr8

       subroutine donne_le_type_0di1(var_) !10-02-20
       integer(kind=1) var_
       s_typevar=1
       end subroutine donne_le_type_0di1

!...................................................................
      subroutine s_polygone !22-02-21
      implicit none
      real ,dimension(:)  ,allocatable ::   lonpol,latpol
      logical ldinmesh

! Les points de la grille A l'interieur du polygone sont masqueS

      open(unit=3,file='../../../GLOBMED/liste_des_polygones')
  370 read(3,'(a)',end=365)texte90  

       open(unit=4,file=trim(texte90))
        read(4,*)i0
        allocate(lonpol(i0))
        allocate(latpol(i0))
        do i=1,i0
         read(4,*)lonpol(i),latpol(i)
         lonpol(i)=lonpol(i)*deg2rad
         latpol(i)=latpol(i)*deg2rad
        enddo
       close(4)

       do j1=0,jmax+1 ; do i1=0,imax+1

       ldinmesh=.FALSE.
       j=i0
       do i=1,i0
        if ((((latpol(i)<=lat_t(i1,j1)).and.(lat_t(i1,j1)<latpol(j))).or.((latpol(j)<=lat_t(i1,j1)).and.(lat_t(i1,j1)<latpol(i)))) &
         .and.(lon_t(i1,j1)<(lonpol(j)-lonpol(i))*(lat_t(i1,j1)-latpol(i))/(latpol(j)-latpol(i))+lonpol(i))) then
        if(ldinmesh)then
         ldinmesh=.false.
        else
         ldinmesh=.true.
        endif
        endif
        j=i
       enddo
       if(ldinmesh) then !m°v°m>
         mask_t(i1,j1,:)=0
       endif                     !m°v°m>

       enddo ; enddo

       deallocate (lonpol)
       deallocate (latpol)
       goto 370


 365  close(3)
      end subroutine s_polygone
!...................................................................

      end module module_s
     
