        subroutine graph_out_bio
!______________________________________________________________________
! SYMPHONIE ocean model
! release 303 - last update: 21-07-21
!______________________________________________________________________
!...............................................................................
! Version date      Description des modifications:
! S.26    02-11-14  mise en service
!         08-11-14  filval deplace
!         20-12-14  masquage avec masque de surface
!         02-04-15  modif sur nfmpi a la demande de malek
! v261    22-10-19  modif de l'extension du nom
! v303    21-07-21  ajouter extension de nom: 
!...............................................................................
!    _________                    .__                  .__             !    (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
      use module_principal ; use module_parallele 
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_bio'
       subroutinedescription='Produces a netcdf output file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! si variable tracer non selectionnee dans notebook_rivers sortir:
       if(grh_out_var(54)/=1)return

! returns texte250 the bio netcdf file:
       call graph_out_bio_file_name 

! writes bio_t in the bionetcdf file:
       call graph_out_bio_write_var 

      end  subroutine graph_out_bio

!.......................................................................

      subroutine graph_out_bio_file_name
      use module_principal ; use module_parallele ; use module_s
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_bio_file_name'
       subroutinedescription='returns texte250 the bio netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!--------------------------------------------
! DEFINIR LE NOM DU FICHIER:
!.....calcul de la date                                            !27/01/03
      call elapsedtime2date(elapsedtime_now,i5,i6,i7,i3,i2,i1)     !01-04-13
      write(texte80(1),                            &
       '(1x,i2,1x,a9,1x,i4,1x,a,i2,a1,i2,a1,i2)')  &
       i7,month(i6),i5,'h:m:s ',i3,':',i2,':',i1

!.....écrire année année mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

!.....écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de séparation, "_"
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04


! Nom des fichiers de variables
      if(flag_nemoffline/=1) then !m°v°m> 
       texte250=dirgraph(1:lname4)//texte90(1:8)              &
                             //'_'//texte90(10:15)            &
                                  //'.ECO3M-S' !22-10-19  
      else                        !m°v°m> 
       texte250=dirgraph(1:lname4)//texte90(1:8)              &
                             //'_'//texte90(10:15)            &
                                  //'.ECO3M-S.BLOOM' !22-10-19  
      endif                       !m°v°m> 

! ajouter extension de nom: !21-07-21
      if(par%rank==0) then !m[°u°]m> 
       k=s_unit(7)
       open(unit=k,file='output_file_extension')
        texte30='' 
        read(k,*,end=471)texte30
  471  close(k)
      endif                !m[°u°]m> 
      call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) 
      texte250=trim(texte250)//trim(texte30)//'.nc' 

!.....puis afficher la date à l'écran:                                     !27/01/03
      if(par%rank==0) then !00000>
       write(6,*)'--------------------------------'
       write(6,*)'Biogeochemical netcdf file date='
       write(6,'(a33)')texte80(1)(1:33)
       write(6,'(a,a)')'File name: ',trim(texte250)
      endif                !00000>

      end subroutine graph_out_bio_file_name

!.......................................................................

      subroutine graph_out_bio_write_var
      use module_principal
      use module_parallele !#MPI
      use module_modeanalysis
      use module_my_outputs
      use pnetcdf
      implicit none
      integer ichoix,looproc_,looprocmax_,loopm_,loop1_,loop2_
      character*6 type_
      integer(kind=MPI_OFFSET_KIND) len_
      double precision :: scalarvalr8(1)
      integer :: scalarvalint(1)

!      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='graph_out_bio_write_var'
       subroutinedescription='writes bio_t in the bionetcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(bio_t))return

      filval=-9999.

      do loop_netcdf=0,1

      count_netcdfvar=0

      if(loop_netcdf==0) then  !§§§§§§§>
       status=nfmpi_create(par%comm2d,texte250 ,nf_clobber, MPI_INFO_NULL,ncid)
      else                     !§§§§§§§>
       status=nfmpi_open(par%comm2d,texte250 ,nf_write, MPI_INFO_NULL,ncid)
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop module_offline erreur 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! time:
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='time'                                                 ! variable
      call kount_to_date(0)                          !16-11-09  time origin corresponds to kount=0

      write(texte80(2),'(a14,i4,a1,i2,4(a1,i2))')                    & !units
      'seconds since ',i5,'-',i6,'-',i7,' ',i3,':',i2,':',i1
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0'
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0'
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'
      if(texte80(2)(32:32)==' ')texte80(2)(32:32)='0'
      texte80(8)=texte80(2)(14:33)                                      ! time_origin : kount=0
      texte80(9)='gregorian'                                            ! calendar
      texte80(3:4)='time'                                               ! long_name
      texte80(5)='T'  ; texte80(6)='time'                               ! axis ; associate
      texte80(7)='double'
      call netcdf_main('_t')

! write variables dimensions:
      call graph_out_trueaxis_nfmpi

      texte80(10)='none'

!.....
! write bio_t
      do vb=1,vbmax !VBVBVBVBVBVB>
      filval=-9999. !08-11-14

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax

         anyvar3d(i,j,k)=bio_t(i,j,k,vb)*mask_t(i,j,k)          & !20-12-14
                                     +(1-mask_t(i,j,k))*filval
        enddo ; enddo ; enddo

      endif                   !=======>
      write(texte80(1),'(a,i0)')'bio',vb
      texte80(2)='none'         ! variable ; units
      write(texte80(3),'(a,i0)')'biogeochemical_variable_number_',vb
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')


      enddo         !VBVBVBVBVBVB>


      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

!     status=nfmpi_put_att_text(ncid,nf_global,'Production',27        &
      status=nfmpi_put_att_text(ncid,nf_global,'Contact',int(22,8) & !02-04-15
      ,'LEGOS-SIROCCO Toulouse')
!     if(status/=0)stop ' Err1 graph_out_bio nfmpi_put_att_text'

!      scalarvalint(1)=modulo_biotimestep ; len_=1  !22-05-17
!      status=nfmpi_put_att_int(ncid_,nf_global                     &
!                              ,'modulo_biotimestep',nf_int &
!                              ,len_,scalarvalint)

       status=nfmpi_put_att_text(ncid,nf_global,'Biogeochemical_model',int(7,8) & 
       ,'ECO3M-S')

      scalarvalr8(1)=dti_fw*modulo_biotimestep ; len_=1
      status=nfmpi_put_att_double(ncid_,nf_global & !29-04-29
       ,'Biogeochemical_time_step',nf_double,len_,scalarvalr8)

      if(ioffline==2) then !ooo>
       if(flag_nemoffline==1) then !>>>
       status=nfmpi_put_att_text(ncid,nf_global,'Physical_model',int(52,8) & 
       ,'NEMO files processed by the SYMPHONIE-BLOOM platform')
       else                        !>>>
       status=nfmpi_put_att_text(ncid,nf_global,'Physical_model',int(24,8) & 
       ,'SYMPHONIE files off-line')
       endif                       !>>>
      else                 !ooo>
       status=nfmpi_put_att_text(ncid,nf_global,'Physical_model',int(17,8) & 
       ,'SYMPHONIE on line')
      endif                !ooo>

      scalarvalr8(1)=dti_fw ; len_=1
      status=nfmpi_put_att_double(ncid_,nf_global & 
       ,'Advection_time_step',nf_double,len_,scalarvalr8)

      scalarvalr8(1)=ratio_bionegdif ; len_=1
      status=nfmpi_put_att_double(ncid_,nf_global & 
       ,'Advection_parameter_ratio_bionegdif',nf_double,len_,scalarvalr8)

       scalarvalint(1)=flag_rmnegval ; len_=1  !02-08-21
       status=nfmpi_put_att_int(ncid_,nf_global                                     &
                               ,'Negative_value_remover_using_flag_rmnegval',nf_int &
                               ,len_,scalarvalint)

! Definition of variables: done.
      status=nfmpi_enddef(ncid)

      endif                  !>>>>>>>>>>>>>>>>>>>


      status=nfmpi_close(ncid)

      enddo  ! fin de boucle sur loop_netcdf

      end  subroutine graph_out_bio_write_var
