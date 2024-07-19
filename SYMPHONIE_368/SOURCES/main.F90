      subroutine main
!______________________________________________________________________
! SYMPHONIE ocean model
! release 343 - last update: 03-04-22
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      implicit none

!.......................................................................
! Version  date     Description
!         23/08/01: amenagements sur OBD_DEMO
!         25/09/01: introduction du cas ISTREAMF>1
!         27/01/03: on appelle directement model_3D et non pas model_serie
!                  + bienvenue filiere offline pour model_ bio
!         27/04/05: Ajout d'un cas demonstration des conditions aux limites
!                   quand on utilise la filiere nest_inout
!         19/01/07: Mise en service de model_1D
!         26/03/07: appel e un model_2D "forcages et imbrications realistes"
!         03/04/07: filiere obc_demo supprimee
!         24/12/08: IAIRSEADEMO remplace par AIRSEAOPTION
!         24/03/09: parallelisation
! 2009.2  03-09-09  plus d'info dans les fichiers infoparallele
! 2009.3  05-10-09  ajout d'un "ifdef parallele"
!         06-11-09  ecrire les infos temporaires dans repertoire tmp
!         13-11-09  ajout d'une securite par rapport au decoupage en
!                   sous-domaines pour la parallelisation
!         14-11-09  "unique" devient "single"
!         15-11-09  "call initier" devient "call initial_main"
! 2010.6  04-02-10  airseademo supprime
! 2010.8  19-05-10  test nr=3 remplace par nr<=3
! 2010.13 15-10-10  conditions limites periodiques avec mpi
! 2010.25 02-02-12  commentaire devant appel a routine initialise_frontiere
!         25-02-12  allocation dynamique
!         31-05-12  allocation dynamique suite
! S.26    28-06-13  mpi subcycling
!         20-08-13  mpi subcycling
!         03-01-14  Une clef de compilation subcycling pour passer du cas subcycling au
!                   cas courant pendant la phase de developpement
!         07-01-14  suite subcycling
!         01-05-14  devient la version qui permettra l'elimination des
!                   procs inutiles et le mpi subcycling
!         16-05-14  logical fplan_grid
!         29-06-14  nbdom=par%nbdom
!         06-07-14  lonlatfile='none' deplace avant lecture notebookgrid
!         23-07-14  notebook_list au format namelist
!         02-08-14  par%comm2d (non connu a ce stade) remplace par mpi_comm_world
!                   Ecrire le synopsis de la simulation
!                   open(unit=3 remplace par une recherche d'unite libre
!         22-08-14  notebook_list renomme notebook_list.f
!         27-10-14  correction automatique du nom lonlatfile si erreur du user
!                   'none' --> 'nofile'
!         09-12-14  ajout obcstatus
!                   deduire dim_river de notebook_river
!         16-12-14  debug obcstatus et (en consequence) creation des variables ieq1, ieqimax etc....
!         23-01-15  c'est main_stop qui arrete la simulation
!         02-04-15  main.F90 devient la subroutine de s_model_main.F90 pour faciliter
!                   l'interfacage avec SDAP (Sequoia Data Assimilation Plateform)
!         25-06-15  oasis
!         12-07-15  ecrire kount dans repertoire tmp
!         16-05-18  modif format d'ecriture
!         03-10-18  ajout fplan1_grid  fplan2_grid discard_lonlat_periodicity
!         16-10-18  reset fplan2_grid
! v252    14-04-19  ajout webcanals_list dans notebook_grid
! v253    26-04-19  webcanals_list='none'              
! v262    23-10-19  Message ecran en fin de simulation
! v286    17-06-20  lecture notebook en serie
! v287    18-07-20  tmpdirname lu dans notebook_grid
! v343    03-04-22  procedure d'arret compatible avec occigen
!...............................................................................
!  _________                    .__                  .__              ! m[°v°]m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____  
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!...............................................................................

#ifdef parallele
         call mpi_init(ierr)                                       !31-05-12
#endif


#ifdef bidon
! IMPORTANT: A CAUSE d'OASIS CES LIGNES SONT FINALEMENT NEUTRALISEES
! En effet A ce stade nbdom sera le nombre total de rank, cumulEs sur
! tous les modeles, y compris les modeles qui ne sont pas symphonie.
! En consequence de quoi la barriere mpi_barrier est impossible car
! cela produit une attente infinie des modeles autres que symphonie qui
! n'ont pas ces lignes...

! Connaitre par%rank et nbdom pour lecture "en serie" des notebook dont
! la lecture est a suivre                                     !17-06-20
      call mpi_comm_rank(mpi_comm_world, par%rank , ierr)     !17-06-20
      call mpi_comm_size(mpi_comm_world, nbdom    , ierr)     !17-06-20
! La lecture en serie se fera par groupe, un groupe contient nbdom/kreadgroupmax ranks 
      kreadgroupmax=min(10,nbdom)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Read notebook_list notebook_grid (size of the grid, mpi dimensions,...)
! dim_river from notebook_river:
! Lecture en serie de groupes !17-06-20
       kread1=0
       do kreadgroup=1,kreadgroupmax
        kread2=nint( real(kreadgroup)/real(kreadgroupmax)*real(nbdom) ) -1
           if(par%rank>=kread1.and.par%rank<=kread2) then !pmx>
            call main_notebook_list_grid !23-01-15
           endif                                          !pmx>
           call mpi_barrier(mpi_comm_world,ierr)
         kread1=kread2+1
       enddo ! kreadgroup
       if(par%rank==0)write(6,*)'Lecture notebook depuis main.F90 OK!'
#endif

       call main_notebook_list_grid !23-01-15

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! mpi subdomains processing:
      call main_mpi_prep
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Initial oasis
      call initial_oasis !25-06-15
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Initial state of the ocean model
      call initial_main                                             !15-11-09
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Biogeochemical model 'offline mode' simulation                      !27/01/03
      if(ioffline==2)then
           call model_offline
           call main_stop
      endif
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! 1D vertical ocean model
      if(flag3d==0) then                                                !19/01/07
               call model_1d
      endif
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! 2D horizontal model based on the external mode equations
      if(i2dh==0) then
               call model_2d
      endif
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      if(i2dh==1) then
       if(kmax<=2) then   !----------->                                  !19-05-10
! 2D horizontal model based on the external mode equations
               call model_2d                                           !26/03/07
       else             !----------->
! 3D ocean model
               call model_3d                                           !27/01/03
       endif            !----------->
      endif
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Stop at the end of the simulation
      call main_stop                        !23-01-15

      end subroutine main

!.......................................................................

      subroutine main_default_mpi_map(rank_)
      use module_principal ; use module_parallele
      implicit none
      integer rank_
      integer,dimension(:,:),allocatable :: &
       rank_ij,imax_ij,jmax_ij
#ifdef synopsis
       subroutinetitle='main_default_mpi_map'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(rank_==0) then !000000000000>

      allocate(rank_ij        (-1:nbdom_imax  ,-1:nbdom_jmax ))
      allocate(imax_ij        ( 0:nbdom_imax-1,0:nbdom_jmax-1))
      allocate(jmax_ij        ( 0:nbdom_imax-1,0:nbdom_jmax-1))

      rank_ij=-2

      x1=real(iglb-2)/real(nbdom_imax) + 2.
      y1=real(jglb-2)/real(nbdom_jmax) + 2.
      write(6,*)'jmax(r8)=',y1
      write(6,*)'imax(r8)=',x1

      k0=0
      do i=0,nbdom_imax-1
      do j=0,nbdom_jmax-1

       imax_ij(i,j)=int(x1)
       jmax_ij(i,j)=int(y1)

       rank_ij(i,j)=k0
       k0=k0+1

      enddo
      enddo

! A ce stade tous les procs ont les memes dimensions mais
! a cause de l'arrondi à l'entier inferieur la somme sur i des imax_ij(i,j)-2*(nbdom_imax-1) n'est pas iglb
! Le deficit: iglb - [ nbdom_imax*(imax_ij(0,0)-2)+2 ] est reparti, "un a un" sur les premiers domaines:

      i0 = iglb - ( nbdom_imax*(   int(x1)   -2)+2 )
      j0 = jglb - ( nbdom_jmax*(   int(y1)   -2)+2 )

      write(6,*)'ecart sum i sur iglb ',i0
      write(6,*)'ecart sum j sur iglb ',j0

      do i=0,i0-1
      do j=0,nbdom_jmax-1
       imax_ij(i,j)=imax_ij(i,j)+1
      enddo
      enddo

      do i=0,nbdom_imax-1
      do j=0,j0-1
       jmax_ij(i,j)=jmax_ij(i,j)+1
      enddo
      enddo

! Verification
      do j=0,nbdom_jmax-1
       k0=2
       do i=0,nbdom_imax-1
        k0=k0+imax_ij(i,j)-2
       enddo
       if(k0/=iglb)stop 'erreur k0'
      enddo
      write(6,*)'k0 iglb = ',k0,iglb

      do i=0,nbdom_imax-1
       k0=2
       do j=0,nbdom_jmax-1
        k0=k0+jmax_ij(i,j)-2
       enddo
       if(k0/=jglb)stop 'erreur k0'
      enddo
      write(6,*)'k0 jglb = ',k0,jglb

! conditions cycliques direction Oi
      if(iperiodicboundary) then !>>>>>>>>
       do j=0,nbdom_jmax-1
! OBC: nbdom_imax
        rank_ij(nbdom_imax,j-1)=rank_ij(0,j-1)
        rank_ij(nbdom_imax,j  )=rank_ij(0,j  )
        rank_ij(nbdom_imax,j+1)=rank_ij(0,j+1)
! OBC: -1
        rank_ij(-1,j-1)=rank_ij(nbdom_imax-1,j-1)
        rank_ij(-1,j  )=rank_ij(nbdom_imax-1,j  )
        rank_ij(-1,j+1)=rank_ij(nbdom_imax-1,j+1)
       enddo
      endif                 !>>>>>>>>
      if(jperiodicboundary) then !>>>>>>>>
       do i=0,nbdom_imax-1
! OBC: nbdom_jmax
        rank_ij(i-1,nbdom_jmax)=rank_ij(i-1,0)
        rank_ij(i  ,nbdom_jmax)=rank_ij(i  ,0)
        rank_ij(i+1,nbdom_jmax)=rank_ij(i+1,0)
! OBC: -1
        rank_ij(i-1,-1)=rank_ij(i-1,nbdom_jmax-1)
        rank_ij(i  ,-1)=rank_ij(i  ,nbdom_jmax-1)
        rank_ij(i+1,-1)=rank_ij(i+1,nbdom_jmax-1)
       enddo
      endif                 !>>>>>>>>

      write(6,*)'Production fichier description_domaine.txt'
      open(unit= 3,file=mpi_map_file_name)

      write( 3,'(3i6,a60)')nbdom_imax,nbdom_jmax  &
                          ,nbdom_imax*nbdom_jmax, &
       ' ! Number of sub-domains in each direction & nbdom'
      write( 3,*)iglb,jglb,' ! iglb jglb'
      do i=0,nbdom_imax-1
      do j=0,nbdom_jmax-1

       i2=imax_ij(0,j)
       do i1=1,i
       i2=i2+imax_ij(i1,j)-2
       enddo

       j2=jmax_ij(i,0)
       do j1=1,j
       j2=j2+jmax_ij(i,j1)-2
       enddo

       if(j==nbdom_jmax-1.and.j2/=jglb)stop 'main.F90 erreur 165'
       if(i==nbdom_imax-1.and.i2/=iglb)stop 'main.F90 erreur 1 3'

      write( 3,*)'!---------------------------'
      write( 3,*)rank_ij(i,j),'             ! sub-domain order number'
      write( 3,*)i,j,' ! sub-domain (i,j) indexes in the mpi space'
      write( 3,*)1 &
      ,'             ! number of cycles in one principal time step'

      write( 3,'(i6,1x,i6,1x,i6,a)')     & !16-05-18
       i2-imax_ij(i,j),i2,imax_ij(i,j)   &
      ,'         ! i start i end imax'

      write( 3,'(i6,1x,i6,1x,i6,a)')     & !16-06-18
       j2-jmax_ij(i,j),j2,jmax_ij(i,j)   &
      ,'         ! j start j end jmax'

      write( 3,'(8(i6,1x),a40)')        & !16-05-18
                rank_ij(i-1,j  )    &
               ,rank_ij(i+1,j  )    &
               ,rank_ij(i  ,j+1)    &
               ,rank_ij(i  ,j-1)    &
               ,rank_ij(i-1,j-1)    &
               ,rank_ij(i+1,j-1)    &
               ,rank_ij(i-1,j+1)    &
               ,rank_ij(i+1,j+1)    &
               ,' ! Neighbors: w e n s ws es wn en'
      i1=i2-2
      j1=j2-2
      enddo
      enddo

      close( 3)

      deallocate(rank_ij)
      deallocate(imax_ij)
      deallocate(jmax_ij)

      endif             !000000000000>
#ifdef parallele
! Il faut que tous les rank attendent que le fichier description soit
! produit avant d'aller plus loin. Par consequent: barriere.
      call mpi_barrier(mpi_comm_symp,k_out)      !02-08-14
#endif

      end subroutine main_default_mpi_map

!............................................................................

      subroutine main_synopsis(title_,description_)
      use module_principal ; use module_parallele ; use module_s
      implicit none
      character(len=*),intent(in) :: title_,description_
      integer :: len0_,lenwrt_=80,loop_,fileunit_
      logical :: openedunit_
! Write the synopsis of the simulation !02-08-14

! Decommenter si on ne veut pas encombrer le sysnopsis avec les routines d'ecriture netcdf
      if(len(title_)>=6) then
       if(title_(1:6)=='netcdf')return
      endif

      if(par%rank==0) then !000000>

       fileunit_=s_unit(7) ! s_unit renvoie un file unit libre, si possible=argument (module_s)

       open(unit=fileunit_,file=trim(tmpdirname)//'synopsis'   &
                                          ,position='append') !05-08-14
        if(title_=='model_3d') then
           write(fileunit_,'(a)')'-------------------------------------'
           write(fileunit_,'(a,i0)')'BEGIN ITERATIVE MODE'
           write(fileunit_,'(a)')'-------------------------------------'
        endif
        if(title_=='model_2d') then
           write(fileunit_,'(a)')'-------------------------------------'
           write(fileunit_,'(a,i0)')'BEGIN ITERATIVE MODE'
           write(fileunit_,'(a)')'-------------------------------------'
        endif
        write(fileunit_,'(5x,a)')trim(title_)
!       write(fileunit_,*)'fileunit_=',fileunit_

        len0_=len(trim(description_))
        if(len0_/=0) then !kkkkkkkk>
           do loop_=1,len0_/lenwrt_+1
            write(fileunit_,'(10x,a)')  &
            description_(1+(loop_-1)*lenwrt_:min(loop_*lenwrt_,len0_))
           enddo
!           write(fileunit_,'(10x,a)')trim(description_)
        endif             !kkkkkkkk>

       close(fileunit_)
       endif                !00000>

      end subroutine main_synopsis

!.......................................................................

      subroutine main_stop
      use module_principal ; use module_parallele ; use module_cpl_oasis
      implicit none
#ifdef synopsis
       subroutinetitle='main_stop '
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

        if(par%rank==0) then !---------->
         open(unit=3,file=trim(tmpdirname)//'kount',position='append')  !12-07-15
          write(3,*)'subroutine main_stop:' &
                   ,' congratulations for this successful run !'     !26-02-13
          write(6,*)
          write(6,*)'subroutine main_stop:' &
                   ,' CONGRATULATIONS FOR THIS SUCCESSFUL RUN!'     !26-02-13
          write(6,*)
       write(6,*)'  ____  __ __  ____        ___   __  _ '
       write(6,*)' |    \|  |  ||    \      /   \ |  |/ ]'
       write(6,*)' |  D  )  |  ||  _  |    |     ||  | / '
       write(6,*)' |    /|  |  ||  |  |    |  O  ||    \ '
       write(6,*)' |    \|  :  ||  |  |    |     ||     |'
       write(6,*)' |  .  \     ||  |  |    |     ||  .  |'
       write(6,*)' |__|\_|\__,_||__|__|     \___/ |__|\_|'
         close(3)

        
        write(6,'(a)')
        write(6,'(a)')
        write(6,'(a,a)')'Open ../../SOURCES/model_name to see what''s ' &
        ,'new in this version of the model' !23-10-19
        

        endif               !---------->
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif

      call cpl_oasis_finalize
      stop

      end subroutine main_stop

!.......................................................................

      subroutine main_notebook_list_grid !23-01-15
      use module_principal
      implicit none
      namelist/notebook_list/directory,nomfichier  !27-07-14
      namelist /notebook_grid/ iglb,jglb,kmax,nbdom_imax,nbdom_jmax  &
      ,iperiodicboundary,jperiodicboundary,dxb,dyb,phi0,longi0       &
      ,latit0,grid_i0,grid_j0,rayonterre,northpole_lon               &
      ,northpole_lat,initialgridfile_txt,vert_axis_conv_direc        &
      ,vert_axis_conv_start,vert_axis_conv_end,hori_axis_conv_start  &
      ,southpole_lon,southpole_lat,lonlatfile,hori_axis_conv_end     &
      ,mpi_map_file_name,mpi_hole_plugging                           &
      ,fplan1_grid,fplan2_grid,discard_lonlat_periodicity            &   !03-10-18
      ,webcanals_list & !14-04-19
      ,tmpdirname !18-07-20

! notebook_list
      allocate(nomfichier(40)) ; nomfichier='s'
      open(100,file='notebook_list.f')             !22-08-14
      read(100,nml=notebook_list)
      close(100)
      do k=1,40
      nomfichier(k)=trim(directory)//'/'//trim(nomfichier(k))
      enddo

! notebook_grid
      lonlatfile='nofile'                                            !06-07-14
      fplan1_grid=.false.                                            !03-10-18
      fplan2_grid=0                                                  !16-10-18
      discard_lonlat_periodicity=0                                   !03-10-18
      mpi_map_file_name='default'                                    !06-07-14
      mpi_hole_plugging='none'                                       !06-07-14
      webcanals_list='none'                                          !26-04-19
      tmpdirname='tmp/' !18-07-20

      open(100,file=trim(nomfichier(2))) ! Lecture du namelist "notebook_grid" !27-07-14
      read(100,nml=notebook_grid)
      close(100)

      kmaxp1=kmax+1
      pi=dacos(-1.d0) ; deg2rad=pi/180.d0 ; rad2deg=180.d0/pi   !03-06-13
      longi0=longi0*deg2rad
      latit0=latit0*deg2rad
      phi0=phi0*deg2rad
      if(lonlatfile=='none')lonlatfile='nofile'  !27-10-14

! dim_river
      open(unit=3,file=nomfichier(4)) ! notebook_river
       read(3,*)dim_river ; dim_river=max(dim_river,1) !09-12-14
      close(3)

#ifndef parallele
      dom_c='single'                                                  !14-11-09
      if(nbdom_imax/=1.or.nbdom_jmax/=1) then !---------------------> !13-11-09
      write(6,*)'Si flag de compilation "parallele" non active dans'
      write(6,*)'makefile alors dans parameter fixer nbdom_imax=1'
      write(6,*)'                                  & nbdom_jmax=1'
      stop
      endif                                   !---------------------<
#endif

      end subroutine main_notebook_list_grid

!..............................................................................

      subroutine main_mpi_prep
      use module_principal ; use module_parallele ; use  module_cpl_oasis
      implicit none
      integer :: rank
#ifdef synopsis
       subroutinetitle='main_mpi_prep'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      mpi_comm_symp = mpi_comm_world

! Creation d'un communicateur local pour le couplage avec oasis3-mct
#ifdef key_oasis
     CALL cpl_oasis_init(mpi_comm_symp)
#endif
#ifdef key_oasis_sym_sym
     CALL cpl_oasis_init_pmx(mpi_comm_symp)
#endif


#ifdef parallele
!     call mpi_comm_rank(mpi_comm_world, rank, ierr)
      call mpi_comm_rank(mpi_comm_symp, rank, ierr)
#endif

! Default number of sub-domains along Oi and Oj axis:
      nbdom=nbdom_imax*nbdom_jmax
!     call mpi_comm_rank(mpi_comm_world, rank, ierr) !28-06-13
      call mpi_comm_rank(mpi_comm_symp, rank, ierr) !28-06-13

! Default regular full mpi map ("full" = including masked rank) archived in 'description_domaine.txt'
      if(mpi_map_file_name=='default') then !>>>>
         mpi_map_file_name=trim(tmpdirname)//'description_domaine.txt'
            call main_default_mpi_map(rank)
      endif                                 !>>>>

!....................
! Get imax, jmax the "local" domain size
      call initialise_parallel_manu(nbdom_imax,nbdom_jmax,iglb,jglb,  & !28-06-13
            imax,jmax,kmax+1,trim(mpi_map_file_name))
      nbdom=par%nbdom                                                   !29-06-14

      call mpi_allreduce(par%subcycle,i0,1,mpi_integer          &
                         ,mpi_max,par%comm2d ,ierr)
      if(i0==1) then
       subcycle_onoff=0 !07-01-14
      else
       subcycle_onoff=1 !07-01-14
      endif
      if(i0<1)stop ' STOP main.F90 i0<1'
      if(subcycle_onoff/=0)stop 'erreur main.F90 subcycle_onoff/=0'


! Verifier que les choix de grille periodique dans notebook_grid sont bien compatible
! avec le fichier de description du decoupage mpi: !20-08-13
       flag_stop=0
       if(par%rank==0) then !000000000000>
       if(iperiodicboundary.and.par%tvoisin(ouest)==mpi_proc_null) then !-verif->
        write(10+par%rank,*) &
        'inconsistent iperiodicboundary in notebook_grid'
        write(10+par%rank,*)' Explication possible:'
        write(10+par%rank,*) &
         ' Le domaine est circulaire periodique et dans le mEme temps' &
        ,' un fichier descrition mpi supprimant les domaines inutiles' &
        ,' est specifiE dans notebook_grid, ce qui n''est pas possible' &
        ,' avec la version actuelle. Changer pour une distribtion mpi' &
        ,' couvrant l''integralitE du domaine'
        flag_stop=1
       endif                                                            !-verif->
       if(jperiodicboundary.and.par%tvoisin(sud)==mpi_proc_null) then   !-verif->
        write(10+par%rank,*) &
        'inconsistent jperiodicboundary in notebook_grid'
        write(10+par%rank,*)' Explication possible:'
        write(10+par%rank,*) &
         ' Le domaine est circulaire periodique et dans le mEme temps' &
        ,' un fichier descrition mpi supprimant les domaines inutiles' &
        ,' est specifiE dans notebook_grid, ce qui n''est pas possible' &
        ,' avec la version actuelle. Changer pour une distribtion mpi' &
        ,' couvrant l''integralitE du domaine'
        flag_stop=1
       endif                                                            !-verif->
       endif                !000000000000>

       call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ;
       if(k0/=0) stop ' STOP dans main.F90 see fort.xxx files' !03-04-22

      end subroutine main_mpi_prep

!..........................................................................

      subroutine initial_oasis
      use  module_cpl_oasis
      implicit none
#ifdef key_oasis
     CALL cpl_oasis_grid(MASTER, mpi_comm_symp)
     CALL cpl_oasis_define
#endif
#ifdef key_oasis_sym_sym
     CALL cpl_oasis_def_part_pmx
#endif
        end subroutine initial_oasis

!..........................................................................


