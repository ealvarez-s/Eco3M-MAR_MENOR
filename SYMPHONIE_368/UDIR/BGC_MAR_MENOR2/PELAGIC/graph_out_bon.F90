        subroutine graph_out
!______________________________________________________________________
! S model
! release S26 - last update: 10-09-14
!______________________________________________________________________

!......................................................................
! Version date      Description des modifications:
!         27/07/01: passage à la coordonnée sigma généralisée.
!         25/08/01: visualisation des variables "bio"
!         06/09/01: usage étendu de ANYVAR2D pour plus de compatibilité
!                   avec une éventuelle version "double precision"
!         29/11/01: ajout de SNSF_Z, SSR_Z, SLHF_Z, SSHF_Z
!         28/03/02: archivage vitesse verticale
!         17/04/02: debug numerotation 3Dgrph...out
!         20/04/02: idem
!         03/10/02: affichage rhp_c+rho-1000.
!         27/01/03: amenagements pour nouvelle facon de faire les sorties
!                   graphique. La construction du nom du fichier graphique
!                   est changée
!         10/07/03: suite, le nom porte la date.
!         11/07/03: suite, cas particulier de l'heure minute secondes nulles
!                   envisagé
!         28/07/03: NOUVELLE VERSION COUPLEE A initial_graph.F et notebook_graph
!         29/07/03: Bienvenue à CALL MOYENNE_TEMPS et à VMI_X et VMI_Y
!         01/08/03: bienvenue à KLOW_Z et à SPONGE_Z
!         02/08/03: bienvenue à VELOBC, ZTAOBC, VMEAOBC, TOBC, SOBC
!         07/08/03: debug dérouleur de boucle
!         09/08/03: impression relax_x relax_y
!         14/08/03: ZTA_EXT_Z(2) au lieu de ZTA_INT_Z(1)
!         21/08/03: debug division par zero
!         03/09/03: modification des conventions pour AMPTIDE, APOTIDE,
!                   PHITIDE, PPOTIDE
!         07/09/03: si affichage du forcage CALL DTVAR_OBC(0)
!         02/10/03: visu harmoniques 3D
!         10/10/03: visu STREAMF_R(0) et STREAMF_R(1)
!                  +debug i/o form='unformatted'
!         21/10/03: visu RELAX_R
!         16/03/04: visu niveaux verticaux grille stepped sigma
!         28/07/04: MTDE et NTDE remplacent MECO et NECO dans les tableaux
!                   de marée
!                   Caractere de separation ":" remplacé par "_" pour
!                   compatibilité avec Windows
!         12/01/01: impression variables calculées - variables forcantes
!         09/02/07: Module graph_out + fichier variable.txt dans notebook_graph
!         26/03/09: TEXTE250 remplace TEXTE60
!         15-05-09  Debug kmin dans vlxobc_x et vlyobc_y
!         01-06-09  Parallelisation: nom de fichier précédé du numero de sous-domaine
! 2009.2  07-09-09  bornes sur i & j pour wstress_x et wstress_y
! 2009.3  10-11-09  bornes sur i & j pour wstress_x et wstress_y
!         12-11-09  verification que la grille des variables de maree = 4
!         27-11-09  ajout sortie vitesse verticale
! 2010.2  21-12-09  modif dimensions parametres nodaux
!         22-12-09  - graph_out_mod renommé module_graph
!                   - ecrire dans tmp
!         24-12-09  - les parametres nodaux n'etantplus integres dans les
!                     amplitudes & phases initiales, ils n'ont plus à etre
!                     enlevés avant ecriture dans fichier
! 2010.3  07-01-10  modif seuil bancs decouvrants
!         13-01-10  modif seuil bancs decouvrants
! 2010.6  02-02-10  renomme lon_t lat_t
! 2010.7  16-02-10  vitesse forward = veldxdz(2)
!         25-02-10  Prise en compte des vitesses de Stokes
!         04-03-10  - Archivage parametres des vagues
!                   - relax_u relax_v relax_p supprimes
! 2010.8  24-03-10  hssh_w+h_w remplace ssh_ext_w
!         03-05-10  bio_c devient bio_t
!         07-05-10  t et s deviennent tem et sal
!         09-05-10  seul le proc 0 ecrit interrupteur
! 2010.9  29-05-10  suppression dsig_u dsig_v
!         07-06-10  ajout dir_wave
! 2010.10 23-06-10  operations sur interrupteur regroupees dans io_switch
! 2010.11 08-08-10  attention aux bornes des boucles....
! 2010.12 27-09-10  attention à ne pas ecrire les parametres des vagues
!                   sil ils ne sont pas dimensionnes
! 2010.13 02-11-10  nouvelle variable, changement ordre dom_c dans nom de fichier
! 2010.14 28-11-10  courant geostrophique
! 2010.22 30-04-11  Ecriture grille globale
! 2010.24 14-12-11  Ecriture traceurs passifs bio_t
!         15-12-11  Ecriture pricipitations
!         30-12-11  correction lien internet (fig5 et pas fig4)
! 2010.25 06-01-12  C.L. pour la visu du wstress
!         20-02-12  Ne pas utiliser de variables non allouees
!         24-02-12  Visu drifters
!         09-03-12  Stop du run_option==-1 introduit dans graph_out
!         20-03-12  debug visu wstress
!         21-06-12  Uilisation de pnetcdf au lieu de netcdf classique,
!                    uniquement pour graph_out_trueaxis_nfmpi
!         27-06-12  attribut axis "T" pour le temps
! S26     06-02-13  ecrire lon lat des points "f"
!         10-04-13  ecrire la vitesse orbitale des vagues
!         14-05-13 impression angle_wave_beam
!         09-06-13 affichage sponge_t "brut"
!         29-06-13 affichage CFL2D
!         05-09-13 ajout parametre vagues
!         02-10-13 le calcul du courant geostrophic dans une routine separee
!                  qui au passage prend mieux en compte les bords (y compris mpi)
!                  du domaine
!         06-05-14 ecrire le fichier de grille que si iteration3d=0
!         19-06-14 ecrire presgrad_u presgrad_v
!         21-06-14 seul le rank 0 ecrit a l'ecran
!         01-07-14 ajout tableau de nudging de la bio
!         05-07-14 k1 k2 remplaces par loop1_ loop2_
!         11-07-14 rap_obc devient timeweightobc(:)
!         14-07-14 elargir boucle sur bio_t pour eviter effets de bords sur graphhiques
!         30-07-14 Idem que le point precedent mais pour T et S
!         20-08-14 possibilite de conversion du courant en nbre de courant si lignes
!                  reperee par 20-08-14 decommentees
!         23-08-14 masquer les valeurs de T et S en zones decouvertes
!         10-09-14 possibilite de conversion du courant moyen en nbre de courant si lignes
!                  reperee par 20-08-14 decommentees
!......................................................................

      use module_principal
      use module_parallele !#MPI
      use module_modeanalysis
      use module_my_outputs
      use pnetcdf
      use ModuleDeclaration    !bio
      implicit none
      integer ichoix,looproc_,looprocmax_,loopm_,loop1_,loop2_
      character*6 type_
!      include 'netcdf.inc'
         integer ,parameter:: ichltotal =35
         integer indice_bio(35),kk      !bio
!        integer indice_bio(24),kk      !bio
#ifdef synopsis
       subroutinetitle='graph_out'
       subroutinedescription='Produces a netcdf output file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Archiver la latitude et la longitude en:
!     type_='real'
      type_='double'
      filval=-9999.

      if(iteration3d/=0) goto 116 ! N'ecrire le fichier de grille qu'à l'initialisation !06-05-14

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Initialisation: ECRIRE LE FICHIER DE GRILLE:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do loop_netcdf=0,1

!     if(loop_netcdf==0)looprocmax_=0
!     if(loop_netcdf==1)looprocmax_=par%nbdom
!     do looproc_   =0,looprocmax_
!     if(par%rank==looproc_   )then !procprocprocproc>

      count_netcdfvar=0

      texte80(3)='none'
      texte80(4)='none'

      texte250=dirgraph(1:lname4)//'grille.nc' !02-11-10
      texte60='grille.nc'                      !02-11-10

      if(loop_netcdf==0) then  !§§§§§§§>
!       status=nfmpi_create(par%comm2d,texte250, nf_clobber, MPI_INFO_NULL, ncid)
        status=nfmpi_create(par%comm2d,texte250, nf_clobber + NF_64BIT_OFFSET , MPI_INFO_NULL, ncid) !10-09-12
      else                     !§§§§§§§>
!       status=nfmpi_open(par%comm2d,texte250, nf_write, MPI_INFO_NULL, ncid)
        status=nfmpi_open(par%comm2d,texte250, nf_write + NF_64BIT_OFFSET , MPI_INFO_NULL, ncid) !10-09-12
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop graph_out erreur 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! write variables dimensions:
      call graph_out_trueaxis

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_t(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_t'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)   ! 'real'

      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_u(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_u'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_v(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_v'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_t(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_t'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_t')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_u(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_u'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_u')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_v(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_v'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
!      if(mask_t(i,j,k)==1) then
          anyvar3d(i,j,k)=depth_t(i,j,k)
!      else
!         anyvar3d(i,j,k)=-9999.
!      endif
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_t'
      texte80(2)='m'                        ! units
      texte80(3)='levels depths at rest'    ! long_name
      texte80(4)='levels_depths_at_rest'    ! standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
!      if(mask_u(i,j,k)==1) then
          anyvar3d(i,j,k)=depth_u(i,j,k)
!      else
!         anyvar3d(i,j,k)=-9999.
!      endif
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_u'
      texte80(2)='m'                        ! units
      texte80(3)='levels depths at rest'    ! long_name
      texte80(4)='levels_depths_at_rest'    ! standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
!      if(mask_v(i,j,k)==1) then
          anyvar3d(i,j,k)=depth_v(i,j,k)
!      else
!         anyvar3d(i,j,k)=-9999.
!      endif
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_v'
      texte80(2)='m'                        ! units
      texte80(3)='levels depths at rest'    ! long_name
      texte80(4)='levels_depths_at_rest'    ! standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_v')


      if(loop_netcdf==1) then !--------->
      do k=1,kmax+1
      do j=0,jmax+1
      do i=0,imax+1
          anyvar3d(i,j,k)=depth_w(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_w'
      texte80(2)='m'                        ! units
      texte80(3)='levels depths at rest'    ! long_name
      texte80(4)='levels_depths_at_rest'    ! standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=h_w(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_w' ; texte80(2)='m'                     ! variable ; units
      texte80(3)='undisturbed water depth'                  ! long_name
      texte80(4)='undisturbed_water_depth'                  !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_u(i,j,kmax  )==1) then
          anyvar2d(i,j)=dy_u(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dy_u'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oj axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oj_axis'   !  standard_name
      texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_u')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_v(i,j,kmax  )==1) then
          anyvar2d(i,j)=dx_v(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dx_v'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oi axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oi_axis'   !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_v')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=dxdy_t(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dxdy_t'   ; texte80(2)='m2'                 ! variable ; units
      texte80(3)='cell box horizontal area'                   ! long_name
      texte80(4)='cell_box_horizontal_area'                                !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')


      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
!         anyvar3d(i,j,k)=dz_t(i,j,max0(k,kmin_w(i,j)),1)
          anyvar3d(i,j,k)=dz_t(i,j,k,1)
       else
          anyvar3d(i,j,k)=-9999.
       endif
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dz_t'   ; texte80(2)='m'                            ! variable ; units
      texte80(3)='undisturbed vertical resolution'                ! long_name
      texte80(4)='undisturbed_vertical_resolution'                ! standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_t')


      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
          anyvar3d(i,j,k)=mask_t(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_t'   ; texte80(2)='none'                   ! variable ; units
      texte80(3)='sea land mask'                                  ! long_name
      texte80(4)='sea_land_mask'                           ! standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

      call netcdf_general_attributes(ncid)

! Definition of variables: done.
      status=nfmpi_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      status=nfmpi_close(ncid)


!     endif                         !procprocprocproc>
#ifdef parallele
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!     call barriere(looproc_   ,3)
#endif
!      enddo ! fin de boucle sur looproc_

      enddo  ! fin de boucle sur loop_netcdf

  116 continue !06-05-14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ECRIRE LE FICHIER DE VARIABLES:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!.....calcul de la date                                                    !27/01/03
      call elapsedtime2date(elapsedtime_now,i5,i6,i7,i3,i2,i1)     !01-04-13

      write(texte80(1),                                                 &
       '(1x,i2,1x,a9,1x,i4,1x,a,i2,a1,i2,a1,i2)')                       &
       i7,month(i6),i5,'h:m:s ',i3,':',i2,':',i1

!.....puis afficher la date à l'écran:                                     !27/01/03
      if(par%rank==0) then !>>>>>>>>
      write(6,*)'--------------------------------'
      write(6,*)'fichier graphique date:'
      write(6,'(a33)')texte80(1)(1:33)
      endif                !>>>>>>>>

!.....écrire année année mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

!.....écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de séparation, "_"
!     écrire ":" dans TEXTE90:
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04


! CONSTRUIRE LE NOM DU FICHIER DE SORTIE GRAPHIQUE                     !27/01/03
      k=int( elapsedtime_now / graphperiod )
      texte30='.nc'    !02-11-10

      texte250=dirgraph(1:lname4)//texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //texte30      !02-11-10

      texte60=                     texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //texte30      !02-11-10

      do loop_netcdf=0,1

!     if(loop_netcdf==0)looprocmax_=0
!     if(loop_netcdf==1)looprocmax_=par%nbdom
!     do looproc_   =0,looprocmax_
!     if(par%rank==looproc_   )then !procprocprocproc>

      count_netcdfvar=0

      if(loop_netcdf==0) then  !§§§§§§§>
!      status=nfmpi_create(par%comm2d,texte250,nf_clobber, MPI_INFO_NULL,ncid)
       status=nfmpi_create(par%comm2d,texte250,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
      else                     !§§§§§§§>
!      status=nfmpi_open(par%comm2d,texte250,nf_write, MPI_INFO_NULL,ncid)
       status=nfmpi_open(par%comm2d,texte250,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop graph_out erreur 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! Define time:
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='time'                                                 ! variable
      call kount_to_date(0)                          !16-11-09  time origin corresponds to kount=0
      write(texte80(2),'(a13,i4,a1,a3,4(a1,i2))')                    & !units
      'seconds from ',i5,'-',month(i6)(1:3),'-',i7,' ',i3,':',i2,':',i1
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0'
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'
      if(texte80(2)(32:32)==' ')texte80(2)(32:32)='0'
      texte80(8)=texte80(2)(14:33)                                      ! time_origin : kount=0
      texte80(9)='gregorian'                                            ! calendar
      texte80(3)='time'                                                 ! long_name
!     texte80(5)='time'  ; texte80(6)='time'                            ! axis ; associate
      texte80(5)='T'  ; texte80(6)='time'                            ! axis ; associate
      texte80(7)='double'
      call netcdf_main('_t')

! write variables dimensions:
!     call graph_out_variables_dimensions
      call graph_out_trueaxis

      texte80(10)='none'

      loop1_=0
      do 201 loop2_=1,grh_nb(1)
          loop1_=loop1_+1
      if(grh_out_var(loop1_).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! Dans cette section: des variables 2d "real" horizontales dependantes du temps:
      texte80(5)='TYX' ; texte80(7)='real'
      var_scalefactor=1. ; var_addoffset=0.


         if(loop2_.eq.1)then !------->

!******************************************************************************
! impression des Flux radiatifs point Z grille 4
!******************************************************************************

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=snsf_w(i,j,1)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='snsf' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=ssr_w(i,j,1)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='ssr' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=slhf_w(i,j,1)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='slhf' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=sshf_w(i,j,1)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='sshf' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !---------> !15-12-11
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=precipi_w(i,j,1)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='precipations' ; texte80(2)='m/s'
      call netcdf_main('_w')

      endif                     !------->

!******************************************************************************


!******************************************************************************
!     Impression maree
      if (kmaxtide>0) then
!******************************************************************************

      if(loop2_.eq.2)then           !2222222>
! ha
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax
         do i=1,imax
          if(mask_t(i,j,kmax  )==1) then
           anyvar2d(i,j)=sqrt(sshtidecosout_w(i,j,ktide)**2+             &
                              sshtidesinout_w(i,j,ktide)**2)
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Ha_'//texte3 ; texte80(2)='m'
        call netcdf_main('_w')
      enddo ! end of loop on ktide

      endif                     !2222222>

      if(loop2_.eq.3)then           !3333333>
! hg
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax
         do i=1,imax
          if(mask_t(i,j,kmax  )==1) then
           x1=-atan2(-sshtidesinout_w(i,j,ktide),sshtidecosout_w(i,j,ktide))
           anyvar2d(i,j)=mod( x1*180./pi,360.*un)
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Hg_'//texte3 ; texte80(2)='degrees'
        call netcdf_main('_w')
      enddo ! end of loop on ktide

      endif                     !3333333>

      if(loop2_.eq.4)then           !4444444>
! ha_pot
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax
         do i=1,imax
          if(mask_t(i,j,kmax  )==1) then
           anyvar2d(i,j)=sqrt(potidecos_w(i,j,ktide)**2+           &
                              potidesin_w(i,j,ktide)**2)
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Ha_pot'//texte3 ; texte80(2)='m'
        call netcdf_main('_w')
      enddo ! end of loop on ktide

      endif                     !4444444>

      if(loop2_.eq.5)then           !5555555>
! hg_pot
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax
         do i=1,imax
          if(mask_t(i,j,kmax  )==1) then
           x1=-atan2(-potidesin_w(i,j,ktide),potidecos_w(i,j,ktide))
           anyvar2d(i,j)=mod( x1*180./pi,360.*un )
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Hg_pot'//texte3 ; texte80(2)='degrees'
        call netcdf_main('_w')
      enddo ! end of loop on ktide

      endif                     !5555555>

      endif

!******************************************************************************
! impression elevation surface
!******************************************************************************
      if(loop2_.eq.6)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=ssh_int_w(i,j,2)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>

      texte80(1)='ssh'       ; texte80(2)='m'               ! variable ; units
      texte80(3)='sea surface height above geoid'       ! long_name
      texte80(4)='sea_surface_height_above_geoid'       ! standard_name
      call netcdf_main('_w')

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression du niveau de fond de la grille
!******************************************************************************

      if(loop2_.eq.7)then !------->
      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=kmin_w(i,j)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='kmin_w' ; texte80(2)='index'
      call netcdf_main('_w')

      endif            !------->
!******************************************************************************


!******************************************************************************
! impression de la zone eponge du mode interne
!******************************************************************************
      if(loop2_.eq.8)then !------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          if(sponge_t(i,j,1)>1.e-20) then
!            anyvar2d(i,j)=1./sponge_t(i,j,1)/3600./24.
             anyvar2d(i,j)=sponge_t(i,j,1)
          else
             anyvar2d(i,j)=0.
          endif
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='sponge_t' ; texte80(2)='none'
      call netcdf_main('_t')

       endif           !------->
!******************************************************************************


!******************************************************************************
! impression de l'elevation de la surface forcante
!******************************************************************************
      if(loop2_.eq.9)then !------->

      if(loop_netcdf==1) then !=======>
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax  ).eq.1) then

       anyvar2d(i,j)=timeweightobc(ssh_id) *sshobc_w(i,j,2)  &  !11-07-14
                +(1.-timeweightobc(ssh_id))*sshobc_w(i,j,0)

      else

       anyvar2d(i,j)=filval

      endif
      if(nest_onoff_in/=0.and.sponge_t(i,j,1)<small1)anyvar2d(i,j)=filval

      enddo
      enddo
      endif                  !=======>
      texte80(1)='sshobc' ; texte80(2)='m'
      call netcdf_main('_w')

      endif           !------->


!******************************************************************************
! impression de la zone de relaxation champs grANDe echelle methode MPV
!******************************************************************************
      if(loop2_.eq.10)then!------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=hssh_w(i,j,2)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>

      texte80(1)='hssh'       ; texte80(2)='m'               ! variable ; units
      texte80(3)='Layer_thickness'       ! long_name
      texte80(4)='Layer_thickness'       ! standard_name
      call netcdf_main('_w')
      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de la fonction de courant
!******************************************************************************
      if(loop2_.eq.11)then!------->
       write(6,*)'grille non prevue pour streamf stop '
       stop

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression elevation surface - elevation surface forcante
!******************************************************************************
      if(loop2_.eq.12)then!------->
      if(loop_netcdf==1) then !=======>
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax  ).eq.1) then

       anyvar2d(i,j)=ssh_int_w(i,j,2)                            &
           -(timeweightobc(ssh_id) *sshobc_w(i,j,2)              &
        +(1.-timeweightobc(ssh_id))*sshobc_w(i,j,0))

      else

       anyvar2d(i,j)=filval

      endif
      if(nest_onoff_in/=0.and.sponge_t(i,j,1)<small1)anyvar2d(i,j)=filval

      enddo
      enddo
      endif                  !=======>
      texte80(1)='ssh_sshobc' ; texte80(2)='m'
      call netcdf_main('_t')

      endif           !------->

!******************************************************************************
! impression des bancs decouvrants:
!******************************************************************************
      if(loop2_==13)then!------->
      if(loop_netcdf==1) then !=======>
       write(*,*)'impression bancs decouvrants'
       do j=0,jmax+1
       do i=0,imax+1
        if(mask_t(i,j,kmax  ).eq.1) then
         anyvar2d(i,j)=1.
         if(hssh_w(i,j,2)<0.5*wetdry_cst2)anyvar2d(i,j)=0.  !13-01-10
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='wetdry_area' ; texte80(2)='m'
      call netcdf_main('_t')
      endif         !------->

!******************************************************************************
! impression du barometre inverse:
!******************************************************************************
      if(loop2_==14)then!------->
      if(loop_netcdf==1) then !=======>
       write(*,*)'impression barometre inverse'
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=-1./grav/rho*(pss_w(i,j,1)-pss_mean(1))
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='barometre_inverse' ; texte80(2)='m'
      call netcdf_main('_w')
       endif         !------->

!******************************************************************************
! impression des parametres des vagues                                !04-03-10
!******************************************************************************
      if(loop2_==15)then!------->
#ifdef stokes

       write(*,*)'impression des parametres des vagues'

      if(loop_netcdf==1.and.allocated(hs_wave_t)) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=hs_wave_t(i,j,1)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='hs_wave_t' ; texte80(2)='m'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(ubw)) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=ubw(i,j)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='orbital_uv' ; texte80(2)='m'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(t_wave_t)) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=t_wave_t(i,j,1)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='t_wave_t' ; texte80(2)='s'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(dir_wave_t)) then !=======>    !07-06-10
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=mod(270.-dir_wave_t(i,j,1)*rad2deg             &
        -atan2(  lat_t(i+1,j)-lat_t(i-1,j)                            & ! rotation symphonie
               ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j)) )*rad2deg &
                          ,un*360.d0)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='dir_wave_t' ; texte80(2)='deg'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(kx_wave_t)) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=kx_wave_t(i,j,1)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='kx_wave_t' ; texte80(2)='1/m'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(ky_wave_t)) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=ky_wave_t(i,j,1)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='ky_wave_t' ; texte80(2)='1/m'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(foc_wave_t)) then !=======> !05-09-13
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=foc_wave_t(i,j,1)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='foc_wave_t' ; texte80(2)='W/m2'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=sshstokes_w(i,j)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='sshstokes_w' ; texte80(2)='m'
      call netcdf_main('_w')

      if(loop_netcdf==1.and.allocated(hsw_wave_t)) then !=======> !05-09-13
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=hsw_wave_t(i,j,1)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
      endif                  !=======>
      texte80(1)='hsw' ; texte80(2)='m'
      call netcdf_main('_w')
#endif
       endif         !------->

!******************************************************************************
! impression drifters
!******************************************************************************
      if(loop2_==16)then!------->
      if(loop_netcdf==1) then !=======>
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,kmax  )==1) then
         anyvar2d(i,j)=-hmax
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       do kbu=1,kbumax
        i=nint(drifter_l(kbu,1)-par%timax(1))
        j=nint(drifter_l(kbu,2)-par%tjmax(1))
        anyvar2d(i,j)=drifter_l(kbu,3)
       enddo
      endif                  !=======>
      texte80(1)='drifters' ; texte80(2)='m'
      call netcdf_main('_w')
      endif         !------->

!******************************************************************************
! impression Bottom Drag Coefficient !14-11-12
!******************************************************************************
      if(loop2_==17)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=cdb_t(i,j)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='cdb_t'       ; texte80(2)='none'  ! variable ; units
      texte80(3)='Bottom drag coefficient'          ! long_name
      texte80(4)='Bottom_drag_coefficient'          ! standard_name
      call netcdf_main('_t')

      endif                     !------->

!******************************************************************************
! impression coriolis
!******************************************************************************
      if(loop2_==18)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=coriolis_t(i,j)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='coriolis_t'       ; texte80(2)='s-1'  ! variable ; units
      texte80(3)='coriolis parameter at t location'          ! long_name
      texte80(4)='coriolis_parameter_at_t_location'          ! standard name
      call netcdf_main('_t')

      endif                     !------->


      if(loop2_==19)then        !------->

      do loopm_=1,countmodemax

      if(loop_netcdf==1) then !=======>
      if(allocated(uv_wmode_t)) then !)))))>
        do j=1,jmax
        do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=c_wave_mode(i,j,loopm_)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0)')'c_wave_mode',loopm_
      texte80(2)='m s**-1'  ! variable ; units
      call netcdf_main('_t')

      do ktide=1,kmaxtide

      if(loop_netcdf==1) then !=======>
      if(allocated(u3dmode_cos_t)) then !)))))>
        do j=1,jmax
        do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
           anyvar2d(i,j)=sqrt(u3dmode_cos_t(i,j,loopm_,ktide)**2 &
                             +u3dmode_sin_t(i,j,loopm_,ktide)**2)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'ua_m',loopm_,'tide',ktide
      texte80(2)='m s**-1'  ! variable ; units
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(u3dmode_cos_t)) then !)))))>
        do j=1,jmax
        do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
           x1=-atan2(-u3dmode_sin_t(i,j,loopm_,ktide)    &
                     ,u3dmode_cos_t(i,j,loopm_,ktide) )
           anyvar2d(i,j)=mod( x1*180./pi,360.*un)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'uh_m',loopm_,'tide',ktide
      texte80(2)='degrees'  ! variable ; units
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(v3dmode_cos_t)) then !)))))>
        do j=1,jmax
        do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
           anyvar2d(i,j)=sqrt(v3dmode_cos_t(i,j,loopm_,ktide)**2 &
                             +v3dmode_sin_t(i,j,loopm_,ktide)**2)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'va_m',loopm_,'tide',ktide
      texte80(2)='m s**-1'  ! variable ; units
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(v3dmode_cos_t)) then !)))))>
        do j=1,jmax
        do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
           x1=-atan2(-v3dmode_sin_t(i,j,loopm_,ktide)    &
                     ,v3dmode_cos_t(i,j,loopm_,ktide) )
           anyvar2d(i,j)=mod( x1*180./pi,360.*un)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'vh_m',loopm_,'tide',ktide
      texte80(2)='degrees'  ! variable ; units
      call netcdf_main('_t')


      enddo ! fin de boucle sur ktide
      enddo ! fin de boucle sur loopm_

      endif                     !------->


! Wave Beam Angle at first level above sea bottom
      if(loop2_==20)then        !------->

      if(loop_netcdf==1                    &
         .and.allocated(angle_wave_beam_w) &
                     ) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=angle_wave_beam_w(i,j)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='wave_beam_angle_w'   ; texte80(2)='none'   ! variable ; units
      texte80(3)='Wave_beam_angle_above_sea_bottom'          ! long_name
      texte80(4)='Wave_beam_angle_above_sea_bottom'          ! standard name
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          ip1=min(i+1,imax+1)
          jp1=min(j+1,jmax+1)
          x1=(h_u(ip1,j)-h_u(i,j))/dx_t(i,j)
          x2=(h_v(i,jp1)-h_v(i,j))/dy_t(i,j)
          anyvar2d(i,j)=sqrt( x1**2 + x2**2 )
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='Hslope'   ; texte80(2)='none'   ! variable ; units
      texte80(3)='bottom_slope'                  ! long_name
      texte80(4)='bottom_slope'                  ! standard name
      call netcdf_main('_t')

      if(loop_netcdf==1                    &
         .and.allocated(angle_wave_beam_w) &
                     ) then !=======>
        do j=0,jmax+1
        do i=0,imax+1

         if(mask_t(i,j,kmax  )==1) then
          ip1=min(i+1,imax+1)
          jp1=min(j+1,jmax+1)
          x1=(h_u(ip1,j)-h_u(i,j))/dx_t(i,j)
          x2=(h_v(i,jp1)-h_v(i,j))/dy_t(i,j)
          anyvar2d(i,j)=angle_wave_beam_w(i,j)-sqrt( x1**2 + x2**2 )
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='wavegeneration'   ; texte80(2)='none'   ! variable ; units
      texte80(3)='Wave_generation_area'          ! long_name
      texte80(4)='Wave_generation_area'          ! standard name
      call netcdf_main('_t')

      endif                     !------->

      if(loop2_==21)then !21212121>                            !29-06-13

       if(loop_netcdf==1) then !=======>
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax  ).eq.1) then !..>

          anyvar2d(i,j)=      &
           2./sqrt((1./dx_t(i,j))**2+(1./dy_t(i,j))**2) &
           /(2.*sqrt( grav*(h_w(i,j)+cfl_sshmax) ) + cfl_umax )

       else                             !..>

        anyvar2d(i,j)=filval

       endif                            !..>
       enddo
       enddo
       endif                  !=======>
       texte80(1)='CFL2D' ; texte80(2)='s'
       call netcdf_main('_w')



      endif          !21212121>


      if(loop2_==25)then !25252525>

          if(tps_ppb_2d.eq.0)tps_ppb_2d=1
           if(tps_strada_2d.eq.0)tps_strada_2d=1


!*****************************************************************
! production primaire
! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax

!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=ppb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='ppb' ; texte80(2)='mgC/m2'
        call netcdf_main('_t')
!        print*,'ppb tracee'
!*****************************************************************
!! production nouvelle
!! mise  zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax




!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=npb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
            enddo
            enddo
       endif
        texte80(1)='npb' ; texte80(2)='gC/m2'
       call netcdf_main('_t')
!        print*,'npb tracee'
!!*****************************************************************
!! production regeneree
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax



!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=rpb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif

        texte80(1)='rpb' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'rpb tracee'
!!***************************************************************
!! nitrification
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax



!            do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!            do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=nitrif2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='nitrif' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'nitrif tracee'
!!*****************************************************************
!! respiration
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax


!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!!! reduction des sorties a la zone de calcul
!!             do j=nbio1,nbio2
!!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=resp2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='resp' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!       print*,'resp tracee'
!!!*****************************************************************
!!! export sous 200 m de MOP
!!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax


!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,1)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='export_mopc' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'export_mopc tracee'
!!*****************************************************************
!! export sous 200 m de MOP
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax


!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,2)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='export_mopc_sed' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'export_mopc_sed tracee'
!!*****************************************************************
!! export sous 200 m de MOP
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax


!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,3)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='export_mopc_turb' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'export_mopc_turb tracee'
!!*****************************************************************
!! export sous 200 m de MOP
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax




!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!            enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,4)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='export_mopc_adv' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'export_mopc_adv tracee'
!*****************************************************************
!!*****************************************************************
!! export sous 200 m de MOD
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax


!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,5)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
             endif
        texte80(1)='export_modc' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'export_modc tracee'
!!*****************************************************************
!! export sous 200 m de MOD
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax



!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,6)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
           endif
        texte80(1)='export_modc_adv' ; texte80(2)='gC/m2'
        call netcdf_main('_t')
!        print*,'export_modc_adv tracee'
!!*****************************************************************
!!*****************************************************************
!! export sous 200 m de MOD
!! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax


!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
!              anyvar2d(i,j)=-9999.
!             enddo
!             enddo
!! reduction des sorties a la zone de calcul
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=exp2d(i,j,7)/tps_strada_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='export_modc_turb' ; texte80(2)='gC/m2'
        call netcdf_main('_t')

      endif          !25252525>


      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
  201 continue


! variables 2D vecteurs:

      do 202 loop2_=1,grh_nb(2)
      loop1_=loop1_+1

      if(grh_out_var(loop1_).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! Dans cette section: des variables 2d "real" horizontales dependantes du temps:
      texte80(5)='TYX' ; texte80(7)='real'

!******************************************************************************
! impression tension du vent
!******************************************************************************
      if(loop2_.eq.1)then !------->

      if(loop_netcdf==1) then !//////////////>
      if(allocated(twox_wave_t)) then !=======>

      do j=1,jmax                                                    !07-09-09
      do i=1,imax+1
         if(mask_u(i,j,kmax  )==1) then
! Vent seulement: penser à oter la partie "vagues"
           anyvar2d(i,j)=                                              &
#ifdef stokes
! Vent seulement: penser à oter la partie "vagues"
        -0.5*(twox_wave_t(i  ,j,1)-tawx_wave_t(i  ,j,1)                 &
             +twox_wave_t(i-1,j,1)-tawx_wave_t(i-1,j,1))+               &
#endif
          wstress_u(i,j,1)

         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo

      else                   !=======>                               !20-03-12

      do j=1,jmax                                                    !07-09-09
      do i=1,imax+1
         if(mask_u(i,j,kmax  )==1) then
! Vent seulement: penser à oter la partie "vagues"
           anyvar2d(i,j)=wstress_u(i,j,1)
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo

      endif                  !=======>
      endif                  !/////////////////////>

! Ajout de C.L. pour visu exploitable en mpi: !06-01-12
      anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
      anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12

      texte80(1)='wstress_u'  ; texte80(2)='N/m2'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !////////////>
      if(allocated(twoy_wave_t)) then !=======>

      do j=1,jmax+1
      do i=1,imax
         if(mask_v(i,j,kmax  )==1) then
           anyvar2d(i,j)=                                              &
#ifdef stokes
! Vent seulement: oter la partie "vagues"
       -0.5*(twoy_wave_t(i,j  ,1)-tawy_wave_t(i,j  ,1)         &
            +twoy_wave_t(i,j-1,1)-tawy_wave_t(i,j-1,1))+       &
#endif
          wstress_v(i,j,1)
!         if(i==imax/2.and.j==jmax/2)write(6,*)'V in graph', wstress_v(i,j,1)

         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo

      else                   !=======>

      do j=1,jmax+1
      do i=1,imax
         if(mask_v(i,j,kmax  )==1) then
           anyvar2d(i,j)=wstress_v(i,j,1)
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo

      endif                  !=======>
      endif                  !/////////////>

! Ajout de C.L. pour visu exploitable en mpi: !06-01-12
      anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
      anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12

      texte80(1)='wstress_v'  ; texte80(2)='N/m2'
      call netcdf_main('_v')

      if(iwve==1) then !ww3ww3ww3ww3ww3>
! Contribution des vagues à la tension de surface:

      if(loop_netcdf==1) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax+1
         if(mask_u(i,j,kmax  )==1) then
           anyvar2d(i,j)=wstress_u(i,j,1)
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
! Ajout de C.L. pour visu exploitable en mpi: !06-01-12
      anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
      anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12
      endif                  !=======>
      texte80(1)='twindx+twox-tawx'  ; texte80(2)='N/m2'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
      do j=1,jmax+1
      do i=1,imax
         if(mask_v(i,j,kmax  )==1) then
           anyvar2d(i,j)=wstress_v(i,j,1)
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
! Ajout de C.L. pour visu exploitable en mpi: !06-01-12
      anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
      anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12
      endif                  !=======>
      texte80(1)='twindy+twoy-tawy'  ; texte80(2)='N/m2'
      call netcdf_main('_v')



      if(loop_netcdf==1.and.allocated(twox_wave_t)) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
#ifdef stokes
           anyvar2d(i,j)=twox_wave_t(i  ,j,1)
#endif
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
      endif                  !=======>
      texte80(1)='twox'  ; texte80(2)='N/m2'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(tawx_wave_t)) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
#ifdef stokes
           anyvar2d(i,j)=tawx_wave_t(i  ,j,1)
#endif
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
      endif                  !=======>
      texte80(1)='tawx'  ; texte80(2)='N/m2'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(twox_wave_t)) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
#ifdef stokes
           anyvar2d(i,j)=twox_wave_t(i  ,j,1)-tawx_wave_t(i  ,j,1)
#endif
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
      endif                  !=======>
      texte80(1)='twox-tawx'  ; texte80(2)='N/m2'
      call netcdf_main('_t')



      if(loop_netcdf==1.and.allocated(twoy_wave_t)) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
#ifdef stokes
           anyvar2d(i,j)=twoy_wave_t(i  ,j,1)
#endif
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
      endif                  !=======>
      texte80(1)='twoy'  ; texte80(2)='N/m2'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(tawy_wave_t)) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
#ifdef stokes
           anyvar2d(i,j)=tawy_wave_t(i  ,j,1)
#endif
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
      endif                  !=======>
      texte80(1)='tawy'  ; texte80(2)='N/m2'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(twoy_wave_t)) then !=======>
      do j=1,jmax                                                    !07-09-09
      do i=1,imax
         if(mask_t(i,j,kmax  )==1) then
#ifdef stokes
           anyvar2d(i,j)=twoy_wave_t(i  ,j,1)-tawy_wave_t(i  ,j,1)
#endif
         else
           anyvar2d(i,j)=filval
         endif
      enddo
      enddo
      endif                  !=======>
      texte80(1)='twoy-tawy'  ; texte80(2)='N/m2'
      call netcdf_main('_t')



      endif            !ww3ww3ww3ww3ww3>


      endif                     !------->

!******************************************************************************


!******************************************************************************
! impression du courant moyen sur la colonne d'eau
!******************************************************************************
      if(loop2_.eq.2)then !2222222>
      if(loop_netcdf==1) then !=======>
          do j=1,jmax
          do i=1,imax+1
           if(mask_u(i,j,kmax  )==1) then
            anyvar2d(i,j)=velbar_u(i,j,1)+velbarstokes_u(i,j,1)   !25-02-10
! converti en nombre de courant: !20-08-14
!           anyvar2d(i,j)=anyvar2d(i,j)*dti_lp/dx_u(i,j)
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>

      texte80(1)='velbar_u' 
      texte80(2)='m/s'
!     texte80(2)='current number' !20-08-14
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
          do j=1,jmax+1
          do i=1,imax
           if(mask_v(i,j,kmax  )==1) then
            anyvar2d(i,j)=velbar_v(i,j,1)+velbarstokes_v(i,j,1)   !25-02-10
! converti en nombre de courant: !20-08-14
!           anyvar2d(i,j)=anyvar2d(i,j)*dti_lp/dy_v(i,j)
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>

      texte80(1)='velbar_v' 
      texte80(2)='m/s'
!     texte80(2)='current number' !20-08-14
      call netcdf_main('_v')

       endif          !2222222>


!******************************************************************************
! impression amplitude du courant de maree
      if(loop2_.eq.3)then !------->
      if (kmaxtide>0) then
!******************************************************************************

!Ua
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax
         do i=1,imax+1
          if(mask_u(i,j,kmax  )==1) then
           anyvar2d(i,j)=sqrt(veltidecosout_u(i,j,ktide)**2+             &
                              veltidesinout_u(i,j,ktide)**2)
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Ua_'//texte3 ; texte80(2)='m'
        call netcdf_main('_u')
      enddo ! end of loop on ktide

!Va
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax+1
         do i=1,imax
          if(mask_v(i,j,kmax  )==1) then
           anyvar2d(i,j)=sqrt(veltidecosout_v(i,j,ktide)**2+             &
                              veltidesinout_v(i,j,ktide)**2)
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Va_'//texte3 ; texte80(2)='m'
        call netcdf_main('_v')
      enddo ! end of loop on ktide


!******************************************************************************
! impression amplitude du courant de maree
      endif
      endif           !------->
!******************************************************************************


!******************************************************************************
! impression retard de phase du courant de maree
      if(loop2_.eq.4)then !------->
      if (kmaxtide>0) then
!******************************************************************************

!Ug
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax
         do i=1,imax+1
          if(mask_u(i,j,kmax  )==1) then
           x1=-atan2(-veltidesinout_u(i,j,ktide),veltidecosout_u(i,j,ktide))
           anyvar2d(i,j)=mod( x1*180./pi,360.*un )
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Ug_'//texte3 ; texte80(2)='degrees'
        call netcdf_main('_u')
      enddo ! end of loop on ktide

!Vg
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=1,jmax+1
         do i=1,imax
          if(mask_v(i,j,kmax  )==1) then
           x1=-atan2(-veltidesinout_v(i,j,ktide),veltidecosout_v(i,j,ktide))
           anyvar2d(i,j)=mod( x1*180./pi,360.*un )
          else
           anyvar2d(i,j)=-9999.
          endif
         enddo
         enddo
        endif                   !-------->
        if(ktide<100)write(texte3,'(i2)')ktide
        if(ktide<10) write(texte3,'(i1)')ktide
        texte80(1)='Vg_'//texte3 ; texte80(2)='degrees'
        call netcdf_main('_v')
      enddo ! end of loop on ktide

!******************************************************************************
! impression retard de phase du courant de maree
      endif
      endif           !------->
!******************************************************************************


!******************************************************************************
! impression courant moyen forcant
!******************************************************************************
      if(loop2_.eq.5)then        !------->

      if(loop_netcdf==1) then !=======>
          do j=1,jmax
          do i=1,imax
           if(mask_u(i,j,kmax  )==1) then
            anyvar2d(i,j)=timeweightobc(vel_id) *velbarobc_u(i,j,2)   &
                     +(1.-timeweightobc(vel_id))*velbarobc_u(i,j,0)
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>
      texte80(1)='velbarobc_u'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
          do j=1,jmax
          do i=1,imax
           if(mask_v(i,j,kmax  )==1) then
            anyvar2d(i,j)=timeweightobc(vel_id) *velbarobc_v(i,j,2)   &
                     +(1.-timeweightobc(vel_id))*velbarobc_v(i,j,0)
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>
      texte80(1)='velbarobc_v' ; texte80(2)='m/s'
      call netcdf_main('_u')

      endif                  !------->

!******************************************************************************


!******************************************************************************
! impression courant moyen - courant moyen forcant
!******************************************************************************
      if(loop2_.eq.6)then        !------->

      if(loop_netcdf==1) then !=======>
          do j=1,jmax
          do i=1,imax
           if(mask_u(i,j,kmax  )==1) then
            anyvar2d(i,j)=velbar_u(i,j,1)+velbarstokes_u(i,j,1) &       !25-02-10
               -(     timeweightobc(vel_id) *velbarobc_u(i,j,2) &
                 +(1.-timeweightobc(vel_id))*velbarobc_u(i,j,0))
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>
      texte80(1)='velbar_u-velbarobc_u'  ; texte80(2)='m/s'           ! variable ; units
      call netcdf_main('_u')


      if(loop_netcdf==1) then !=======>
          do j=1,jmax
          do i=1,imax
           if(mask_v(i,j,kmax  )==1) then
            anyvar2d(i,j)=velbar_v(i,j,1)+velbarstokes_v(i,j,1) &       !25-02-10
               -(     timeweightobc(vel_id) *velbarobc_v(i,j,2) &
                 +(1.-timeweightobc(vel_id))*velbarobc_v(i,j,0))
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>
      texte80(1)='velbar_v_velbarobc_v'  ; texte80(2)='m/s'           ! variable ; units
      call netcdf_main('_v')


      endif           !------->
!******************************************************************************


!******************************************************************************
! impression courant de Stokes moyen                                  !25-02-10
!******************************************************************************
      if(loop2_.eq.7)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=1,jmax
        do i=1,imax
         if(mask_u(i,j,kmax  )==1) then
          anyvar2d(i,j)=velbarstokes_u(i,j,1)
         else
          anyvar2d(i,j)=-9999.
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='velbarstokes_u' ; texte80(2)='m/s'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do j=1,jmax
        do i=1,imax
         if(mask_v(i,j,kmax  )==1) then
          anyvar2d(i,j)=velbarstokes_v(i,j,1)
         else
          anyvar2d(i,j)=-9999.
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='velbarstokes_v' ; texte80(2)='m/s'
      call netcdf_main('_v')

      endif                    !------->

!******************************************************************************

!******************************************************************************
! impression courant moyen - courant de Stokes moyen                  !25-02-10
!******************************************************************************
      if(loop2_.eq.8)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=1,jmax
        do i=1,imax
         if(mask_u(i,j,kmax  )==1) then
          anyvar2d(i,j)=velbar_u(i,j,1)                               !25-02-10
         else
          anyvar2d(i,j)=-9999.
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='ubar_hat' ; texte80(2)='m/s'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do j=1,jmax
        do i=1,imax
         if(mask_v(i,j,kmax  )==1) then
          anyvar2d(i,j)=velbar_v(i,j,1)                              !25-02-10
         else
          anyvar2d(i,j)=-9999.
         endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='vbar_hat' ; texte80(2)='m/s'
      call netcdf_main('_v')

      endif                    !------->

!******************************************************************************


      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
!******************************************************************************

  202 continue


! variables 3D 1/2 niveaux scalaires:

      do 203 loop2_=1,grh_nb(3)
       loop1_=loop1_+1

       if(grh_out_var(loop1_).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! Dans cette section: des variables "real" 3d dependantes du temps
      texte80(5)='TZYX'  ; texte80(7)='real'
      var_scalefactor=1. ; var_addoffset=0.

!******************************************************************************
! impression de la temperature
!******************************************************************************
      if(loop2_.eq.1)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=0,jmax+1 !30-07-14
        do i=0,imax+1
         if(mask_t(i,j,kmax)==1.and.wetmask_t(i,j)/=0.) then !23-08-14
          anyvar3d(i,j,k)=tem_t(i,j,max0(k,kmin_w(i,j)),1)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='tem'   ; texte80(2)='degrees_Celsius'
      texte80(3)='sea_water_potential_temperature'
      texte80(4)='sea_water_potential_temperature'
      call netcdf_main('_t')

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de la salinite
!******************************************************************************
      if(loop2_.eq.2)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=0,jmax+1 !30-07-14
        do i=0,imax+1
         if(mask_t(i,j,kmax)==1.and.wetmask_t(i,j)/=0.) then !23-08-14
           anyvar3d(i,j,k)=sal_t(i,j,max0(k,kmin_w(i,j)),1)
         else
           anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>

      texte80(1)='sal'       ; texte80(2)='1e-3'                  ! variable ; units
      texte80(3)='sea water salinity'                             ! long_name
      texte80(4)='sea_water_salinity'                             ! standard_name
      call netcdf_main('_t')

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de la densite
!******************************************************************************
      if(loop2_.eq.3)then !------->

      if(loop_netcdf==1) then !=======>
!     if(ioffline.eq.2)call density(0)
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
         if (mask_t(i,j,k     ).ne.0) then
!         anyvar3d(i,j,k)=rhp_t(i,j,max0(k,kmin_w(i,j)))+rho-1000.
          anyvar3d(i,j,k)=rhp_t(i,j,k)+rho-1000.
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='rhp'       ; texte80(2)='kg/m3'                     ! variable ; units
      call netcdf_main('_t')

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression des traceurs passifs bios
!******************************************************************************
      if(loop2_.eq.4)then !------->


         data indice_bio /iZooNanoC,iZooMicroC,iZooMesoC           &
                      ,iSyneC,iNanoC,iDiaC                         &
                      ,iSyneN,iNanoN,iDiaN                         &
                      ,iSyneP,iNanoP,iDiaP,iDiaSi                  &
                      ,iSyneChl,iNanoChl,iDiaChl,iBactC            &
                      ,iSMOPC,iSMOPN,iSMOPP,iSMOPSI,ISMOPCHL                &
                      ,iLMOPC,iLMOPN,iLMOPP,iLMOPSI                &
                      ,iMODC,iMODN,iMODP                           &
                      ,iNitrate,iAmmonium,iPhosphate,iSilice       &
                      ,iOxygen,ichltotal/   ! MODIF LE 26 03 2013 FAYCAL

       if (vbmax.gt.0) then !<><><><><>->

!     do vb=1,vbmax !- vb loop ->
!      do kk=1,23 !- vb loop ->
      do kk=1,35 !- vb loop ->

          vb=indice_bio(kk)

             if(vb.eq.1)texte80(1)='zoonanoc'
             if(vb.eq.2)texte80(1)='zoomicroc'
             if(vb.eq.3)texte80(1)='zoomesoc'
             if(vb.eq.4)texte80(1)='synec'
             if(vb.eq.5)texte80(1)='synen'
             if(vb.eq.6)texte80(1)='synep'
             if(vb.eq.7)texte80(1)='synechl'
             if(vb.eq.8)texte80(1)='nanoc'
             if(vb.eq.9)texte80(1)='nanon'
             if(vb.eq.10)texte80(1)='nanop'
             if(vb.eq.11)texte80(1)='nanochl'
             if(vb.eq.12)texte80(1)='diac'
             if(vb.eq.13)texte80(1)='dian'
             if(vb.eq.14)texte80(1)='diap'
             if(vb.eq.15)texte80(1)='diachl'
             if(vb.eq.16)texte80(1)='diasi'
             if(vb.eq.17)texte80(1)='bactc'
             if(vb.eq.18)texte80(1)='smopc'
             if(vb.eq.19)texte80(1)='smopn'
             if(vb.eq.20)texte80(1)='smopp'
             if(vb.eq.21)texte80(1)='smopchl'
             if(vb.eq.22)texte80(1)='smopsi'
             if(vb.eq.23)texte80(1)='lmopc'
             if(vb.eq.24)texte80(1)='lmopn'
             if(vb.eq.25)texte80(1)='lmopp'
             if(vb.eq.26)texte80(1)='lmopsi'
             if(vb.eq.27)texte80(1)='modc'
             if(vb.eq.28)texte80(1)='modn'
             if(vb.eq.29)texte80(1)='modp'
             if(vb.eq.30)texte80(1)='nitrate'
             if(vb.eq.31)texte80(1)='ammonium'
             if(vb.eq.32)texte80(1)='phosphate'
             if(vb.eq.33)texte80(1)='silice'
             if(vb.eq.34)texte80(1)='oxygen'
             if(vb.eq.35)texte80(1)='chltot'



      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
         if (mask_t(i,j,kmax+1)==1) then
          if(vb<=34) &
          anyvar3d(i,j,k)=bio_t(i,j,max0(k,kmin_w(i,j)),vb)   !14-12-11
          if(vb==35) &
               anyvar3d(i,j,k)=bio_t(i,j,max0(k,kmin_w(i,j)),isynechl) + &
                               bio_t(i,j,max0(k,kmin_w(i,j)),inanochl) + &
                               bio_t(i,j,max0(k,kmin_w(i,j)),idiachl)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>

      if(vb<100)write(texte3,'(i2)')vb
      if(vb<10) write(texte3,'(i1)')vb
!     texte80(1)='tracer_'//texte3 ; texte80(2)='unity not defined'
      texte80(2)='unity not defined'

      call netcdf_main('_t')

      enddo         !- vb loop ->

      endif                    !<><><><><>->
      endif           !------->



!******************************************************************************


!******************************************************************************
! impression de la temperature de forcage
!******************************************************************************
      if(loop2_.eq.5)then !------->

      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=0,jmax+1 !30-07-14
         do i=0,imax+1
          if (mask_t(i,j,kmax  ).ne.0) then
          anyvar3d(i,j,k)=timeweightobc(trc_id) *temobc_t(i,j,max0(k,kmin_w(i,j)),2) &
                     +(1.-timeweightobc(trc_id))*temobc_t(i,j,max0(k,kmin_w(i,j)),0)
          else
           anyvar3d(i,j,k)=-9999.
          endif
         enddo
         enddo
         enddo
      endif                  !=======>
      texte80(1)='temobc_t'     ; texte80(2)='degc'                       ! variable ; units
      call netcdf_main('_t')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de la salinite de forcage
!******************************************************************************
      if(loop2_.eq.6)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=0,jmax+1 !30-07-14
        do i=0,imax+1
         if (mask_t(i,j,kmax  ).ne.0) then
          anyvar3d(i,j,k)=timeweightobc(trc_id) *salobc_t(i,j,max0(k,kmin_w(i,j)),2) &
                     +(1.-timeweightobc(trc_id))*salobc_t(i,j,max0(k,kmin_w(i,j)),0)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='salobc_t'     ; texte80(2)='psu'                      ! variable ; units
      call netcdf_main('_t')

      endif           !------->


!******************************************************************************
! Bio: solution imposee dans couche lateralle de nudging
!******************************************************************************
      if(loop2_.eq.7)then !------->                                     !01-07-14

      do vb=1,vbmax           !-<>->
       if(loop_netcdf==1) then !=======>

       anyvar3d(:,:,:)=-9999.
       x1=(1.-rap_biobc)
       x2=    rap_biobc

! Si frontière ouest (i=1) ouverte:
       if(par%tvoisin(ouest) == mpi_proc_null) then !--------->

        if(biobc_type(vb)==2) then !*************>
         do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
         i=i0 !01-05-13
         do k=1,kmax
          do j=1,jmax
            anyvar3d(i,j,k)=             &
                            x1*bio_relax_west(j,i0,k,vb,1) &
                           +x2*bio_relax_west(j,i0,k,vb,2)
          enddo
         enddo
         enddo
        endif                      !*************>

       endif                                        !--------->

! Si frontière est (i=imax) ouverte:
       if(par%tvoisin(est  ) == mpi_proc_null) then !--------->

        if(biobc_type(vb)==2) then !*************>
         do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
         i=imax+1-i0 !01-05-13
         do k=1,kmax
          do j=1,jmax
            anyvar3d(i,j,k)=             &
                            x1*bio_relax_east(j,i0,k,vb,1) &
                           +x2*bio_relax_east(j,i0,k,vb,2)
          enddo
         enddo
         enddo
        endif                      !*************>

       endif                                        !--------->

! Si frontière sud (j=1) ouverte:
       if(par%tvoisin(sud  ) == mpi_proc_null) then !--------->
        i1=1 ; i2=imax
        if(par%tvoisin(ouest)==mpi_proc_null)i1=bio_relax_size+1
        if(par%tvoisin(est  )==mpi_proc_null)i2=imax-bio_relax_size

        if(biobc_type(vb)==2) then !*************>
         do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
         j=j0
         do k=1,kmax
          do i=i1,i2                                   !Ne pas appliquer 2 fois le rappel!01-05-13
            anyvar3d(i,j,k)=              &
                            x1*bio_relax_south(i,j0,k,vb,1) &
                           +x2*bio_relax_south(i,j0,k,vb,2)

          enddo
         enddo
         enddo

        endif                      !*************>

       endif                                        !--------->

! Si frontière nord (j=jmax) ouverte:
       if(par%tvoisin(nord ) == mpi_proc_null) then !--------->
        i1=1 ; i2=imax
        if(par%tvoisin(ouest)==mpi_proc_null)i1=bio_relax_size+1
        if(par%tvoisin(est  )==mpi_proc_null)i2=imax-bio_relax_size

        if(biobc_type(vb)==2) then !*************>
         do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
         j=jmax+1-j0
         do k=1,kmax
          do i=i1,i2                                   !Ne pas appliquer 2 fois le rappel!01-05-13
            anyvar3d(i,j,k)=              &
                            x1*bio_relax_north(i,j0,k,vb,1) &
                           +x2*bio_relax_north(i,j0,k,vb,2)
          enddo
         enddo
         enddo
        endif                      !*************>

       endif                                        !--------->

       endif                   !=======>

      write(texte30,'(a12,i0)')'tracer_relax',vb
      texte80(1)=trim(texte30) ; texte80(2)='unity not defined'

      call netcdf_main('_t')
      enddo                   !-<>->

      endif           !------->


!******************************************************************************
! impression des niveaux verticaux de grille
!******************************************************************************
      if(loop2_.eq.8)then !------->                                     !16/03/04

      if(loop_netcdf==1) then !=======>
          do k=1,kmax+1
          do j=1,jmax
          do i=1,imax
             if (mask_t(i,j,kmax  ).ne.0) then
                anyvar3d(i,j,k)=float(max0(k,kmin_w(i,j)))
             else
                anyvar3d(i,j,k)=-9999.
             endif
          enddo
          enddo
          enddo
      endif                  !=======>
      texte80(1)='vertical_level'  ; texte80(2)='index'               ! variable ; units
      call netcdf_main('_t')


      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de temperature - forcage temperature
!******************************************************************************
      if(loop2_.eq.9)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=0,jmax+1 !30-07-14
        do i=0,imax+1
         if(mask_t(i,j,kmax)==1.and.wetmask_t(i,j)/=0.) then !23-08-14
          anyvar3d(i,j,k)=tem_t(i,j,max0(k,kmin_w(i,j)),2)            &
                      - (timeweightobc(trc_id) *temobc_t(i,j,max0(k,kmin_w(i,j)),2) &
                    +(1.-timeweightobc(trc_id))*temobc_t(i,j,max0(k,kmin_w(i,j)),0))
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
!     texte80(1)='tem-temobc'  ; texte80(2)='degc'              ! variable ; units
      texte80(1)='tem_temobc'  ; texte80(2)='degc'              ! variable ; units
      call netcdf_main('_t')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de salinite - forcage salinite
!******************************************************************************
      if(loop2_.eq.10)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=0,jmax+1 !30-07-14
        do i=0,imax+1
         if(mask_t(i,j,kmax)==1.and.wetmask_t(i,j)/=0.) then !23-08-14
          anyvar3d(i,j,k)=sal_t(i,j,max0(k,kmin_w(i,j)),1)            &
                      - (timeweightobc(trc_id) *salobc_t(i,j,max0(k,kmin_w(i,j)),2) &
                    +(1.-timeweightobc(trc_id))*salobc_t(i,j,max0(k,kmin_w(i,j)),0))
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='sal_salobc'  ; texte80(2)='psu'               ! variable ; units
      call netcdf_main('_t')

      endif           !------->

      if(loop2_.eq.11)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
         do i=1,imax
         if (mask_t(i,j,kmax  ).ne.0) then
          anyvar3d(i,j,k)=dz_t(i,j,k,1)-hz_w(i,j,1)*dsig_t(i,j,k)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='dz-dzref'  ; texte80(2)='m'               ! variable ; units
      call netcdf_main('_t')

      endif           !------->

      if(loop2_.eq.12)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
         do i=1,imax
         if (mask_t(i,j,k).ne.0) then
! Cas où l'on a oté la partie lagrangienne barotrope:
!         anyvar3d(i,j,k)=       &
!              (omega_w(i,j,k+1,1)-omega_w(i,j,k,1)              )**2  &
!         /(   (omega_w(i,j,k+1,1)-omega_w(i,j,k,1)              )**2  &
!             +(( dz_t(i,j,k,2)-hz_w(i,j,2)*dsig_t(i,j,k)              &
!                -dz_t(i,j,k,0)+hz_w(i,j,0)*dsig_t(i,j,k))/dti_lp)**2  &
!             +small1    )
! Cas où l'on en prend en compte toute la partie lagrangienne:
          anyvar3d(i,j,k)=       &
               (omega_w(i,j,k+1,1)-omega_w(i,j,k,1)              )**2  &
          /(   (omega_w(i,j,k+1,1)-omega_w(i,j,k,1)              )**2  &
              +(( dz_t(i,j,k,2)                                        &
                 -dz_t(i,j,k,0)                          )/dti_lp)**2  &
              +small1    )
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='eulerian_ratio'  ; texte80(2)='none'               ! variable ; units
      call netcdf_main('_t')
      endif           !------->

      if(loop2_.eq.13)then !------->
      if(loop_netcdf==1) then !=======>
        do j=1,jmax
        do i=1,imax
        if (mask_t(i,j,kmax  ).ne.0) then
          sum1=ssh_int_w(i,j,1)
          sum2=sum1
          do k=kmax,1,-1
           sum1=sum1-dsig_t(i,j,k)*hz_w(i,j,1)
           anyvar3d(i,j,k)=depth_t(i,j,k)-0.5*(sum1+sum2)
           sum2=sum1
          enddo
        else
          do k=kmax,1,-1
           anyvar3d(i,j,k)=-9999.
          enddo
        endif
        enddo
        enddo
      endif                  !=======>
      texte80(1)='z-zref'  ; texte80(2)='m'               ! variable ; units
      call netcdf_main('_t')
      endif           !------->

      if(loop2_.eq.14)then !------->

      call modeanalysis_modeprojection

      do loopm_=1,countmodemax

      if(loop_netcdf==1) then !=======>
      if(allocated(uv_wmode_t)) then !)))))>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
        if (mask_t(i,j,kmax  ).ne.0) then
           anyvar3d(i,j,k)=uv_wmode_t(i,j,k,loopm_)
        else
           anyvar3d(i,j,k)=-9999.
        endif
        enddo
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0)')'uv_wmode_t',loopm_
      texte80(2)='none'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(uv_wmode_t)) then !)))))>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
        if (mask_t(i,j,kmax  ).ne.0) then
           anyvar3d(i,j,k)=uv_wmode_t(i,j,k,loopm_)*ucoefmode_t(i,j,loopm_)
        else
           anyvar3d(i,j,k)=-9999.
        endif
        enddo
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0)')'umode',loopm_
      texte80(2)='m.s-1'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(uv_wmode_t)) then !)))))>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
        if (mask_t(i,j,kmax  ).ne.0) then
           anyvar3d(i,j,k)=uv_wmode_t(i,j,k,loopm_)*vcoefmode_t(i,j,loopm_)
        else
           anyvar3d(i,j,k)=-9999.
        endif
        enddo
        enddo
        enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0)')'vmode',loopm_
      texte80(2)='m.s-1'
      call netcdf_main('_t')


      enddo ! fin de boucle sur loopm_

      endif           !------->


!******************************************************************************
! Orbital wave current amplitude !10-04-13
!******************************************************************************
      if(loop2_.eq.15)then !------->
      if(loop_netcdf==1) then !=======>
       if(allocated(t_wave_t).and.  &
          allocated(hs_wave_t).and. &
          allocated(k_wave_t)) then !"""""""">
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
         if (mask_t(i,j,k)==1) then

          anyvar3d(i,j,k)=                                             &
                  pi/t_wave_t(i,j,1)*hs_wave_t(i,j,1)/sqrt(2.)         &
           *cosh(    k_wave_t(i,j,1)*( h_w(i,j)+depth_t(i,j,k))       )&
           /sinh(    k_wave_t(i,j,1)* hz_w(i,j,1)                     )

         else

          anyvar3d(i,j,k)=-9999.

         endif
        enddo
        enddo
        enddo
       endif                        !"""""""">
      endif                  !=======>
      texte80(1)='UVorbital'  ; texte80(2)='meterpersecond'         ! variable ; units
      call netcdf_main('_t')

      endif           !------->
!******************************************************************************
!******************************************************************************


!******************************************************************************
       endif               !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  203 continue


!******************************************************************************
! variables 3D 1/2 niveaux vecteurs:
!******************************************************************************
      do 204 loop2_=1,grh_nb(4)
      loop1_=loop1_+1

      if(grh_out_var(loop1_).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! Dans cette section: des variables "real" 3d dependantes du temps
      texte80(5)='TZYX' ; texte80(7)='real'

!******************************************************************************
! impression du courant 3D
!******************************************************************************
      if(loop2_.eq.1)then !------->
      if(loop_netcdf==1) then !=======>

        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1                          !08-08-10
         if (mask_u(i,j,kmax  ).ne.0) then
          anyvar3d(i,j,k)=vel_u(i,j,k,1)  &
                   +velstokes_u(i,j,k,1)       !25-02-10
! converti en nombre de courant: !20-08-14
!         anyvar3d(i,j,k)=anyvar3d(i,j,k)*dti_lp/dx_u(i,j)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='vel_u'  
      texte80(2)='m s-1'    
!     texte80(2)='current number' !20-08-14
      texte80(3)='sea_water_x_velocity_at_u_location'
      texte80(4)='sea_water_x_velocity_at_u_location'
      call netcdf_main('_u')


      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1                          !08-08-10
        do i=1,imax
         if (mask_v(i,j,kmax  ).ne.0) then
          anyvar3d(i,j,k)=vel_v(i,j,k,1)  &
                   +velstokes_v(i,j,k,1)       !25-02-10
! converti en nombre de courant: !20-08-14
!         anyvar3d(i,j,k)=anyvar3d(i,j,k)*dti_lp/dy_v(i,j)
          else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='vel_v'  
      texte80(2)='m s-1' 
!     texte80(2)='current number' !20-08-14
      texte80(3)='sea_water_y_velocity_at_v_location'
      texte80(4)='sea_water_y_velocity_at_v_location'
      call netcdf_main('_v')

      endif           !------->

!******************************************************************************


!******************************************************************************
! impression du courant geostrophique
!******************************************************************************
      if(loop2_.eq.2)then !------->
      if(loop_netcdf==1) then !=======>
      call graph_out_geostrophic_current
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         if (mask_u(i,j,k)/=0) then
          anyvar3d(i,j,k)=anyv3d(i,j,k,3)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='u_geos'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
      call graph_out_geostrophic_current
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         if (mask_v(i,j,k)/=0) then
          anyvar3d(i,j,k)=anyv3d(i,j,k,4)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='v_geos'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_v')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression du courant forcant:
!******************************************************************************
      if(loop2_.eq.3)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1
           if (mask_u(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=timeweightobc(vel_id) *velobc_u(i,j,k,2)         &
                         +(1.-timeweightobc(vel_id))*velobc_u(i,j,k,0)
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='velobc_u'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_u')


      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1
        do i=1,imax
           if (mask_v(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=timeweightobc(vel_id) *velobc_v(i,j,k,2)         &
                         +(1.-timeweightobc(vel_id))*velobc_v(i,j,k,0)
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='velobc_v'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_v')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression courant - courant forcant
!******************************************************************************
      if(loop2_.eq.4)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1
           if (mask_u(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=vel_u(i,j,k,1)     &
                       +velstokes_u(i,j,k,1)     & !25-02-10
                -(timeweightobc(vel_id) *velobc_u(i,j,k,2)         &
             +(1.-timeweightobc(vel_id))*velobc_u(i,j,k,0))
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='vel_u_velobc_u'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_u')


      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1
        do i=1,imax
           if (mask_v(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=vel_v(i,j,k,1)     &
                       +velstokes_v(i,j,k,1)     & !25-02-10
                -(timeweightobc(vel_id) *velobc_v(i,j,k,2)         &
             +(1.-timeweightobc(vel_id))*velobc_v(i,j,k,0))
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='vel_v_velobc_v'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_v')

      endif           !------->
!******************************************************************************

!******************************************************************************
! impression courant de Stokes
!******************************************************************************
      if(loop2_.eq.5)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1
           if (mask_u(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=velstokes_u(i,j,k,1)   !25-02-10
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='velstokes_u'  ; texte80(2)='m/s'                    ! variable ; units
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1
        do i=1,imax
           if (mask_v(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=velstokes_v(i,j,k,1)   !25-02-10
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='velstokes_v'  ; texte80(2)='m/s'               ! variable ; units
      call netcdf_main('_v')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression courant - courant de Stokes
!******************************************************************************
      if(loop2_.eq.6)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1
           if (mask_u(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=vel_u(i,j,k,1)         !25-02-10
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='u_hat'  ; texte80(2)='m/s'                    ! variable ; units
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1
        do i=1,imax
           if (mask_v(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=vel_v(i,j,k,1)         !25-02-10
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='v_hat'  ; texte80(2)='m/s'                          ! variable ; units
      call netcdf_main('_v')

      endif           !------->

!******************************************************************************
! Baroclinic PGF !19-06-14
!******************************************************************************
      if(loop2_.eq.7)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1
           if (mask_u(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=presgrad_u(i,j,k)/dx_u(i,j)
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='presgrad_u'  ; texte80(2)='m/s'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1
        do i=1,imax
           if (mask_v(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=presgrad_v(i,j,k)/dy_v(i,j)
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='presgrad_v'  ; texte80(2)='m/s'
      call netcdf_main('_v')

      endif           !------->
!******************************************************************************


!******************************************************************************

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
!******************************************************************************

  204 continue


!******************************************************************************
! variables 3D niveaux entiers scalaires:
!******************************************************************************

      do 205 loop2_=1,grh_nb(5)
      loop1_=loop1_+1

      if(grh_out_var(loop1_).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! Dans cette section: des variables "real" 3d dependantes du temps
      texte80(5)='TZYX' ; texte80(7)='real'


!******************************************************************************
! impression de la diffusivite verticale:
!******************************************************************************
      if(loop2_.eq.1)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax+1
        do j=1,jmax
        do i=1,imax
           if (mask_t(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=kh_w(i,j,max0(k,kmin_w(i,j)))
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='kh'  ; texte80(2)='m2/s'                          ! variable ; units
      call netcdf_main('_w')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de l'energie cinetique turbulente
!******************************************************************************
      if(loop2_.eq.2)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax+1
        do j=1,jmax
        do i=1,imax
           if (mask_t(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=tken_w(i,j,max0(k,kmin_w(i,j)))
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='tken'  ; texte80(2)='(m/s)2'                     ! variable ; units
      call netcdf_main('_w')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de la longueur de melange
!******************************************************************************
      if(loop2_.eq.3)then !------->

      if(loop_netcdf==1) then !=======>

      if(iturbulence==0)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkll_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==1)anyvar3d(1:imax,1:jmax,1:kmax+1)=epsn_w(1:imax,1:jmax,1:kmax+1)

      do k=1,kmax+1
      do j=1,jmax
      do i=1,imax
       if(mask_t(i,j,kmax  )==0)anyvar3d(i,j,k)=-9999.
      enddo
      enddo
      enddo

      endif                  !=======>
      if(iturbulence==0)texte80(1)='tkll'
      if(iturbulence==1)texte80(1)='eps'
      if(iturbulence==0)texte80(2)='m'
      if(iturbulence==1)texte80(2)='m2/s3'
      call netcdf_main('_w')

      endif           !------->

!******************************************************************************


!******************************************************************************
! impression de la longueur de dissipation
!******************************************************************************
      if(loop2_.eq.4)then !------->

      if(loop_netcdf==1) then !=======>

      if(iturbulence==0)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkle_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==1)anyvar3d(1:imax,1:jmax,1:kmax+1)=  km_w(1:imax,1:jmax,1:kmax+1)

      do k=1,kmax+1
      do j=1,jmax
      do i=1,imax
       if(mask_t(i,j,kmax  )==0)anyvar3d(i,j,k)=-9999.
      enddo
      enddo
      enddo

      endif                  !=======>
      if(iturbulence==0)texte80(1)='tkle'
      if(iturbulence==1)texte80(1)='km'
      if(iturbulence==0)texte80(2)='m'
      if(iturbulence==1)texte80(2)='m2/s'
      call netcdf_main('_w')

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de omega
!******************************************************************************
      if(loop2_.eq.5)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax+1
        do j=1,jmax
        do i=1,imax
           if (mask_t(i,j,kmax  ).ne.0) then
              anyvar3d(i,j,k)=omega_w(i,j,max0(k,kmin_w(i,j)),1)
! converti en nombre de courant !20-08-14
!             anyvar3d(i,j,k)=anyvar3d(i,j,k)*dti_lp &
!            /(depth_t(i,j,k)-depth_t(i,j,k-1))
           else
              anyvar3d(i,j,k)=-9999.
           endif
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='omega'  
      texte80(2)='m s-1'           
!     texte80(2)='current number'   !20-08-14        
      texte80(3)='relative_vertical_sea_water_velocity_at_w_location'
      texte80(4)='relative_vertical_sea_water_velocity_at_w_location'
      call netcdf_main('_w')

      endif           !------->

!******************************************************************************
! impression de la vitesse verticale                                  !27-11-09
!******************************************************************************
      if(loop2_.eq.6)then !------->

      if(loop_netcdf==1) then !=======>
        do j=1,jmax
        do i=1,imax
        if (mask_t(i,j,kmax  ).ne.0) then !>>>>>

        do k=1,kmax+1  !---k-->

         if(k==1)then
          sum0=0.
          sum1=0.
         else
          sum0=sum0+dz_t(i,j,k-1,before)
          sum1=sum1+dz_t(i,j,k-1,now)
         endif

         sum2=0.
         do loop1_=max0(k-1,1),min0(k,kmax)  !--loop1_-->
         do i1=i,i+1
           sum2=sum2+                                                   &
           vel_u(i1,j,loop1_,1)*(depth_t(i1,j,loop1_)-depth_t(i1-1,j,loop1_))       &
           /dx_u(i1,j)
         enddo
         do j1=j,j+1
           sum2=sum2+                                                   &
           vel_v(i,j1,loop1_,1)*(depth_t(i,j1,loop1_)-depth_t(i,j1-1,loop1_))       &
           /dy_v(i,j1)
         enddo
         enddo                           !--loop1_--<

! w=omega+dz/dx*u+dz/dy*v+dz/dt:
         anyvar3d(i,j,k)=                                               &
           mask_t(i,j,k)*(                                              &
          omega_w(i,j,k,1)+0.25*sum2+(sum1-sum0)/dti_fw)


        enddo      !---k--<

        else                         !>>>>>

         do k=1,kmax+1
          anyvar3d(i,j,k)=-9999.
         enddo

        endif                        !>>>>>
        enddo
        enddo
      endif                  !=======>
      texte80(1)='w'  ; texte80(2)='m s-1'                     ! variable ; units
      texte80(3)='vertical_sea_water_velocity_at_w_location'
      texte80(4)='vertical_sea_water_velocity_at_w_location'
      call netcdf_main('_w')

      endif           !------->


      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
!******************************************************************************

  205 continue


!******************************************************************************
! variables 3D niveaux entiers vecteurs:
!******************************************************************************
      do 206 loop2_=1,grh_nb(6)
      loop1_=loop1_+1
      if(grh_out_var(loop1_).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! Dans cette section: des variables "real" 3d dependantes du temps
      texte80(5)='TZYX' ; texte80(7)='real'

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
  206 continue

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

      call netcdf_general_attributes(ncid)

! Definition of variables: done.
      status=nfmpi_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      status=nfmpi_close(ncid)


!     endif                         !procprocprocproc>
#ifdef parallele
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!     call barriere(looproc_   ,2)
#endif
!     enddo ! fin de boucle sur looproc_

      enddo  ! fin de boucle sur loop_netcdf

!******************************************************************************
! CREATION DU FICHIER GRAPHIQUE
! FIN.
!******************************************************************************

!...............................................................................
! Reset des tableaux de moyenne temporelle
! Début:
!...............................................................................
      if(grh_out_mi==1)call moyenne_temps(0)                         !29/07/03
!...............................................................................
! Reset des tableaux de moyenne temporelle
! Fin.
!...............................................................................

!_______________________________________________________!
      kpvwave=0                                         !
      call io_switch('w')                                             !23-06-10
      if(par%rank==0)then !>>>>>
       write(6,*)'-----------------------'
       write(6,*)'fichier graphique ecrit'
       write(6,'(a)')' nom:'
       write(6,'(a)')texte250
      endif               !>>>>>
!_______________________________________________________!

      if(run_option==-1) then  !------> !09-03-12
      if(par%rank==0)write(6,*)'RUN stopped after initial' & !21-06-14
                              ,' state as requested in notebooktime'
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      stop
      endif                    !------>


!  **  ECRITURE D'UN FICHIER PVWAVE   FIN
      end subroutine graph_out

!...................................................................................................

      subroutine graph_out_variables_dimensions
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_variables_dimensions'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=1+i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='imax_t' ; texte80(5)='X'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=1+j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='jmax_t' ; texte80(5)='Y'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='kmax_t' ; texte80(5)='Z'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=1+i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='imax_w' ; texte80(5)='X'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=1+j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='jmax_w' ; texte80(5)='Y'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='kmax_w' ; texte80(5)='Z'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='imax_u' ; texte80(5)='X'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=1+j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='jmax_u' ; texte80(5)='Y'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='kmax_u' ; texte80(5)='Z'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=1+i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='imax_v' ; texte80(5)='X'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='jmax_v' ; texte80(5)='Y'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='kmax_v' ; texte80(5)='Z'
      texte80(7)='integer'; texte80(2:4)='undefined'
      call netcdf_main('_v')


      end subroutine graph_out_variables_dimensions

!...................................................................................................

      subroutine graph_out_trueaxis
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_trueaxis'
       subroutinedescription='Writes the grid indexes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      texte80(2)='undefined'

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_t' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_t' ; texte80(5)='Y'
!     texte80(1)='jmax_t' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_t' ; texte80(5)='Z'
!     texte80(1)='kmax_t' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_w' ; texte80(5)='X'
!     texte80(1)='imax_w' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_w' ; texte80(5)='Y'
!     texte80(1)='jmax_w' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k-0.5 ; enddo
      endif                   !--------->
      texte80(1)='nk_w' ; texte80(5)='Z'
!     texte80(1)='kmax_w' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=-0.5+i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_u' ; texte80(5)='X'
!     texte80(1)='imax_u' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_u' ; texte80(5)='Y'
!     texte80(1)='jmax_u' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_u' ; texte80(5)='Z'
!     texte80(1)='kmax_u' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_v' ; texte80(5)='X'
!     texte80(1)='imax_v' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_v_location'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=-0.5+j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_v' ; texte80(5)='Y'
!     texte80(1)='jmax_v' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_v_location'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_v' ; texte80(5)='Z'
!     texte80(1)='kmax_v' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_v_location'
      call netcdf_main('_v')

      end subroutine graph_out_trueaxis

!...................................................................................................

      subroutine graph_out_trueaxis_nfmpi
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_trueaxis_nfmpi'
       subroutinedescription='Writes the grid indexes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      texte80(2)='undefined'

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_t' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_t' ; texte80(5)='Y'
!     texte80(1)='jmax_t' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_t' ; texte80(5)='Z'
!     texte80(1)='kmax_t' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_w' ; texte80(5)='X'
!     texte80(1)='imax_w' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_w' ; texte80(5)='Y'
!     texte80(1)='jmax_w' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k-0.5 ; enddo
      endif                   !--------->
      texte80(1)='nk_w' ; texte80(5)='Z'
!     texte80(1)='kmax_w' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=-0.5+i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_u' ; texte80(5)='X'
!     texte80(1)='imax_u' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_u' ; texte80(5)='Y'
!     texte80(1)='jmax_u' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_u' ; texte80(5)='Z'
!     texte80(1)='kmax_u' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_v' ; texte80(5)='X'
!     texte80(1)='imax_v' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_v_location'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=-0.5+j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_v' ; texte80(5)='Y'
!     texte80(1)='jmax_v' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_v_location'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_v' ; texte80(5)='Z'
!     texte80(1)='kmax_v' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_v_location'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
       do i=0,imax+1 ; anyv3d(i,1,1,1)=-0.5+i+par%timax(1) ; enddo
      endif                   !--------->
      texte80(1)='ni_f' ; texte80(5)='X'
      texte80(7)='real'; texte80(3:4)='x_grid_index_at_f_location'
      call netcdf_main('_f')

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; anyv3d(1,j,1,2)=-0.5+j+par%tjmax(1) ; enddo
      endif                   !--------->
      texte80(1)='nj_f' ; texte80(5)='Y'
      texte80(7)='real'; texte80(3:4)='y_grid_index_at_f_location'
      call netcdf_main('_f')

      if(loop_netcdf==1) then !--------->
       do k=0,kmax+1 ; anyv3d(1,1,k,3)=k ; enddo
      endif                   !--------->
      texte80(1)='nk_f' ; texte80(5)='Z'
      texte80(7)='real'; texte80(3:4)='z_grid_index_at_f_location'
      call netcdf_main('_f')

      end subroutine graph_out_trueaxis_nfmpi

!........................................................................................

      subroutine graph_out_geostrophic_current
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_geostrophic_current'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Tres important si frontiere ouvertes: reset a zero:
      anyv3d(:,:,:,1:4)=0.
      do k=1,kmax
      do j1=2,jmax-1
      do i1=2,imax
         anyv3d(i1,j1,k,1)=                                            &
                           (+presgrad_u(i1,j1,k)                       &
                      +grav*(ssh_int_w(i1,j1,1)-ssh_int_w(i1-1,j1,1))) &
                        /(0.5*(coriolis_t(i1,j1)+coriolis_t(i1-1,j1))) &
                        *mask_u(i1,j1,k)/dx_u(i1,j1)
      enddo
      enddo
      enddo
      do k=1,kmax
      do j1=2,jmax
      do i1=2,imax-1
         anyv3d(i1,j1,k,2)=                                            &
                           (-presgrad_v(i1,j1,k)                       &
                      -grav*(ssh_int_w(i1,j1,1)-ssh_int_w(i1,j1-1,1))) &
                        /(0.5*(coriolis_t(i1,j1)+coriolis_t(i1,j1-1))) &
                        *mask_v(i1,j1,k)/dy_v(i1,j1)
      enddo
      enddo
      enddo
!#ifdef parallele
!      ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
!      call echange('x ',anyv3d,lb4,ub4,1) ! 2 pour echange 4eme arg = 2
!      ! C.L. i=imax+1 i=1 j=jmax j=1
!      ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
!      call echange('y ',anyv3d,lb4,ub4,2) ! 2 pour echange 4eme arg = 2
!      ! C.L. i=imax+1 i=1 j=jmax j=1
!#endif
      call obc_int_anyv3d(1,'u1') ! Echange 'u1' sur anyv3d(:,:,:,1)
      call obc_int_anyv3d(2,'v1') ! Echange 'v1' sur anyv3d(:,:,:,2)

! Step 3: assemblage des elements du courant geostrophiques
! Composante u:
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
         anyv3d(i,j,k,3)=0.25*(                     &
                               anyv3d(i-1,j  ,k,2)  &
                              +anyv3d(i-1,j+1,k,2)  &
                              +anyv3d(i  ,j  ,k,2)  &
                              +anyv3d(i  ,j+1,k,2)  &
                              )
      enddo
      enddo
      enddo
! Composante v:
      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
         anyv3d(i,j,k,4)=0.25*(                     &
                               anyv3d(i+1,j  ,k,1)  &
                              +anyv3d(i+1,j-1,k,1)  &
                              +anyv3d(i  ,j  ,k,1)  &
                              +anyv3d(i  ,j-1,k,1)  &
                              )
      enddo
      enddo
      enddo
!#ifdef parallele
!      ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
!      call echange('x ',anyv3d,lb4,ub4,3) ! 2 pour echange 4eme arg = 2
!      ! C.L. i=imax+1 i=1 j=jmax j=1
!      ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
!      call echange('y ',anyv3d,lb4,ub4,4) ! 2 pour echange 4eme arg = 2
!      ! C.L. i=imax+1 i=1 j=jmax j=1
!#endif
      call obc_int_anyv3d(3,'u1') ! Echange 'u1' sur anyv3d(:,:,:,1)
      call obc_int_anyv3d(4,'v1') ! Echange 'v1' sur anyv3d(:,:,:,2)

      end subroutine graph_out_geostrophic_current
