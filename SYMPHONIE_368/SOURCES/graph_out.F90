        subroutine graph_out
!______________________________________________________________________
! SYMPHONIE ocean model
! release 368 - last update: 16-04-23
!______________________________________________________________________
!    _________                    .__                  .__     m[°v°]m !
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
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
!         20-10-14 ecrire densite potentielle zref=0 meme si eos_tkezref<0
!                  ecriture de uwind_t vwind_t
!         03-11-14 choix d'afficher les valeurs des traceurs dans les zones
!                  seches....
!         18-11-14 modif archivage rhp
!         25-11-14 suppression de lignes inutiles et amelioration d'un message d'erreur 
!                  a l'ecran
!         14-12-14 modif du temps pour rhpobc
!         20-12-14 bio_t "bidonne"
!         20-01-15 boucle sans test sur mask
!         24-01-15 texte80(11)='mask_obc_z1' pour masquer c.l. z1
!         29-01-15 courant geostrophique sans interpolation: courant u
!                  point v et courant sur point u....
!         20-02-15 echange mpi sur anyvar2d pour visu
!         16-03-15 suite point precedent
!         18-04-15 qq mise a jour pour nemo offline
!         24-05-15 ecriture z0_w
!         04-06-15 extension des boucles Ha Hg
!         04-07-15 aiguillage ssh_w ou ssh_int_w
!         21-07-15 la prise en compte de la subdivision du temps dans le mode externe
!                  dans le nom du fichier graphique permet de nommer distinctement
!                  les fichiers produits a l'interieur d'une meme sequence du mode
!                  externe
!         15-10-15 Ajout tableau temref et cie
!         12-11-15 ajout temref_z et salref_z
!         22-12-15 reduction de boucles
!         14-01-16 amenagement calcul courant geostrophique pour cas test baroclinic_jet
!         18-01-16 archivage courant de surface, vorticite de surface etc...
!         28-01-16 presgrad_u presgrad_v passent A 4 dimensions
!         09-02-16 ajout sponge_u sponge_v
!         11-02-16 ecriture time avec convention CF
!         17-02-16 extension des boucles Hapot Hgpot
!         20-02-16 C.L. parametre analyse harmonique barocline
!         17-03-16 C.L. uwind vwind
!         10-04-16 Ecriture de rhpref et cie en mEme temps que rhp
!         14-04-16 Ajout analyse 3d de la pression
!         19-04-16 Blinder les divisions par coriolis_t
!         08-05-16 masque de bord de domaine pour tableaux de flux
!         12-05-16 debug de la date dans le fichier netcdf
!         16-07-16 ajout sortie wetmask_u wetmask_v
!         08-09-16 modification de la convention ecriture date de reference dans 
!                  fichiers netcdf
!         13-09-16 modif de taille de boucle pour eviter depassement de memoire    
!         28-09-16 ajout wstresb
!         28-11-16 ajout Bulk Ridchardson Number
!         17-12-16 C.L. sortie eps
!         30-12-16 omega en m/s et omega en nombre de courant
!         12-01-17 u et v exprimes en nbre de courant
!         14-01-17 masquer bord z1 de rhpref_t
!         20-01-17 T et S, masque 3D
!         14-02-17 mask sur kh_w
!         29-03-17 cas iturbulence=2
!         25-04-17 dir_wave avec 3 conventions d'angle
!         03-05-17 ajout kmerged-kmin
!         25-05-17 ajout de ssh_tide, velbar_tide
!         23-10-17 blindage iteration2d_max_now=0
!         04-12-17 debug nom de fichier
!         17-04-18 tidecos, tidesin A la place de tidecosout, tidesinout
!         19-05-18 filval kh_w
!         05-07-18 nf_clobber+ NF_64BIT_DATA permet de creer des fichiers pour de grandes grilles
!         21-07-18 retour A NF_64BIT_OFFSET
! v259    17-08-19 Ajout d'une extension de nom de fichier
! v261    22-10-19 amelioration du point precedent
! v269    06-12-19 ajout velbot_u et v
! v285    07-06-20 call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20
! v287    14-08-20 obc_int_anyv3d renomme obc_mpi_anyv3d
! v290    30-10-20 drifters 2D et 3D
! v292    19-11-20 ajout total surface current, courant de maree, courant de stokes de surface
! v295    29-01-21 - retablir dz_t dans les sorties
!                  - (Brunt-Vaisala frequency)**2: !29-01-21
! v296    02-03-21  timeabovewater_w !02-03-21
! v299    18-03-21  ajout temf et salf
!         19-03-21  Brunt-Vaisala frequency suite....
! v300    24-03-21  rhp_1000
! v309    30-09-21  ajout temf/ratio_negdif_ver et salf/ratio_negdif_ver
!         02-10-21  dT/dxs , dT/dys , dS/dxs , dS/dys
! v315    10-12-21  call omega_vertical_velocity(2) !10-12-21
! v316    15-12-21  suite point precedent
! v328    23-02-22  appels aux EOS comme dans presgrad
! v332    28-02-22  ajout upwindwetdry_t(i,j)
! v335    04-03-22  mise A jour attribut rhp
! v340    26-03-22  sponge_u(:,:,0) sponge_v(:,:,0) sans dimensions
! v347    12-05-22  if(flag_bathy_update==1) ecrire bathy variable en
!                   meme temps que ssh
! v348    27-05-22  ajout sigma_fric_w
!         20-06-22  des visu pour le cas ifdef bilan_s_3d
!         21-06-22  des visu pour le cas ifdef bilan_s_3d
! v351    16-08-22  Au choix mask_t ou wetmask...
! v357    21-10-22  /max(relax_int,small2) 
! v359    31-10-22  ecrire temlwd, sallwd
! v363    13-01-23  ajout current_number
! v366    14-02-23  reecriture sans if
! v367    19-02-23  ecriture 2pi/k_wave_t
! v368    16-04-23  reecriture sans if 
!......................................................................

      use module_principal
      use module_parallele !#MPI
      use module_modeanalysis
      use module_my_outputs
      use module_s
      use module_parameter_sedim , only : l_sedim_s
      use pnetcdf
      implicit none
      integer ichoix,looproc_,looprocmax_,loopm_,loop1_,loop2_
      character*6 type_
!      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='graph_out'
       subroutinedescription='Produces a netcdf output file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(     imodeltrc==1  &
          .or.imodelbio==1  &
          .or.l_sedim_s)call graph_out_bio   !08-11-14
           if(l_sedim_s)call graph_out_sedim !14-12-16

      

! Archiver la latitude et la longitude en:
!     type_='real'
      type_='double'
      filval=-9999.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ECRIRE LE FICHIER DE VARIABLES:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!.....calcul de la date                                                    !27/01/03
      call elapsedtime2date(elapsedtime_now                         &
           +mod(iteration2d-1,max(iteration2d_max_now,1))*dte_fw    & !23-10-17
                           ,i5,i6,i7,i3,i2,i1)                 !01-04-13

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

      if(par%rank==0) then !m[°u°]m> !17-08-19
       k=s_unit(7)
       open(unit=k,file='output_file_extension')
        texte30='' !22-10-19
        read(k,*,end=471)texte30
  471  close(k)
      endif                !m[°u°]m> 
      call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20
      texte30=trim(texte30)//'_symphonie.nc' !22-10-19

! CONSTRUIRE LE NOM DU FICHIER DE SORTIE GRAPHIQUE                     !27/01/03
      k=int( elapsedtime_now / graphperiod )
!     texte30='.nc'    !02-11-10

      texte250=dirgraph(1:lname4)//texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //texte30      !02-11-10

      texte60=                     texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //texte30      !02-11-10

!     if(flag_nh2d==1) then
!      k10=1000000+iteration3d*iteration2d_max_now+iteration2d
!      write(texte30,'(i0,a)')k10,'.nc'
!      texte250=dirgraph(1:lname4)//trim(texte30)
!     endif
      if(flag_nh3d/=0) then
       k10=1000000+iteration3d
!      write(texte30,'(i0,a)')k10,'.nc'
!      texte250=dirgraph(1:lname4)//trim(texte30)
       write(texte250,'(i0)')k10
       texte250=dirgraph(1:lname4)//trim(texte250)//trim(texte30) !17-08-19
      endif

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
      if(status/=0) then !-----> !25-11-14
         write(6,'(a,a)')'File or directory not found:',trim(texte250)
         stop ' stop in subroutine graph_out'
      endif              !-----> 

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! Define time:
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='time'                                                 ! variable
      call kount_to_date(0)                          !16-11-09  time origin corresponds to kount=0

      write(texte80(2),'(a14,i4,5(a1,i2))')                    & !08-09-16
      'seconds since ',i5,'-',i6,'-',i7,' ',i3,':',i2,':',i1
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0'                  !12-05-16
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0'                  !12-05-16
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


!******************************************************************************
! impression des Flux radiatifs point Z grille 4
!******************************************************************************
      if(loop2_==1)then !pmxpmx>

      if(loop_netcdf==1) then !--------->
       do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)=snsf_w(i,j,1)*mask_t(i,j,kmax) &
                           +filval*(1-mask_t(i,j,kmax))
       enddo ; enddo
       anyvar2d(0:imax+1:imax+1,:)=filval !30-04-16
       anyvar2d(:,0:jmax+1:jmax+1)=filval !30-04-16
      endif                  !--------->
      texte80(1)='snsf' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
       do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)= ssr_w(i,j,1)*mask_t(i,j,kmax) &
                           +filval*(1-mask_t(i,j,kmax))
       enddo ; enddo
       anyvar2d(0:imax+1:imax+1,:)=filval !30-04-16
       anyvar2d(:,0:jmax+1:jmax+1)=filval !30-04-16
      endif                  !--------->
      texte80(1)='ssr' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
       do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)=slhf_w(i,j,1)*mask_t(i,j,kmax) &
                           +filval*(1-mask_t(i,j,kmax))
       enddo ; enddo
       anyvar2d(0:imax+1:imax+1,:)=filval !30-04-16
       anyvar2d(:,0:jmax+1:jmax+1)=filval !30-04-16
      endif                  !--------->
      texte80(1)='slhf' ; texte80(2)='w/m2'
      call netcdf_main('_w')


      if(loop_netcdf==1) then !--------->
       do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)=sshf_w(i,j,1)*mask_t(i,j,kmax) &
                           +filval*(1-mask_t(i,j,kmax))
       enddo ; enddo
       anyvar2d(0:imax+1:imax+1,:)=filval !30-04-16
       anyvar2d(:,0:jmax+1:jmax+1)=filval !30-04-16
      endif                  !--------->
      texte80(1)='sshf' ; texte80(2)='w/m2'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !---------> !15-12-11
       do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)=precipi_w(i,j,1)*mask_t(i,j,kmax) &
                              +filval*(1-mask_t(i,j,kmax))
       enddo ; enddo
       anyvar2d(0:imax+1:imax+1,:)=filval !30-04-16
       anyvar2d(:,0:jmax+1:jmax+1)=filval !30-04-16
      endif                  !--------->
      texte80(1)='precipitations' ; texte80(2)='m/s'
      call netcdf_main('_w')

      if(allocated(pss_w)) then !*****>
       if(loop_netcdf==1) then !---------> !15-12-11
        do j=1,jmax  ; do i=1,imax  
          anyvar2d(i,j)=                                             &
          -(1./grav/rho)*(pss_w(i,j,1)-pss_mean(1))*mask_t(i,j,kmax) & ! pss_w convertie en BI
!                         pss_w(i,j,1)*mask_t(i,j,kmax)              & ! pss_w brute
                                   -(1-mask_t(i,j,kmax))*9999.
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval !30-04-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval !30-04-16
       endif                  !--------->
       texte80(1)='inverse_barometer' ; texte80(2)='m'
!      texte80(1)='pss_w' ; texte80(2)='Pascal'
       call netcdf_main('_t')
      endif                       !*****>

      if(allocated(teta2_t)) then !*****>
       if(loop_netcdf==1) then !---------> !15-12-11
        do j=0,jmax+1
        do i=0,imax+1
          anyvar2d(i,j)=(teta2_t(i,j,1)-273.15)*mask_t(i,j,kmax) &
                                     -9999.*(1.-mask_t(i,j,kmax))
        enddo
        enddo
       endif                  !--------->
       texte80(1)='teta2_t' ; texte80(2)='degrees_Celsius'
       call netcdf_main('_t')
      endif                       !*****>

      if(allocated(teta0_t)) then !*****>
      if(allocated(airseafile_prvtime)) then !*****>

       if(loop_netcdf==1) then !---------> !15-12-11
        x2=(elapsedtime_now           -airseafile_prvtime(t2m_id))      &
          /(airseafile_nextime(t2m_id)-airseafile_prvtime(t2m_id))
        x2=min(max(x2,0.d00),1.d00) 
        x1=1.-x2
        do j=0,jmax+1
        do i=0,imax+1
          anyvar2d(i,j)= &
!         mask_t(i,j,kmax)*                             &
                           (x1*teta0_t(i,j,1)           &
!                       /((1.e5/(pss_w(i,j,0)))**.286)  & ! pour revenir a T
                           +x2*teta0_t(i,j,2)           &
!                       /((1.e5/(pss_w(i,j,2)))**.286)  & ! pour revenir a T
                                             -273.15)    
!         -(1-mask_t(i,j,kmax))*9999.
        enddo
        enddo
       endif                  !--------->
       texte80(1)='teta0_t' ; texte80(2)='degrees_Celsius'
       call netcdf_main('_t')
      endif                                  !*****>
      endif                       !*****>

      if(allocated(teta2delta_t)) then !*****>
       if(loop_netcdf==1) then !---------> !15-12-11
        do j=1,jmax  
        do i=1,imax  
          anyvar2d(i,j)=teta2delta_t(i,j,1)
        enddo
        enddo
       endif                  !--------->
       texte80(1)='teta2delta_t' ; texte80(2)='degree'
       call netcdf_main('_t')
      endif                           !*****>

      if(allocated(q2_t)) then    !*****>
       if(loop_netcdf==1) then !---------> !15-12-11
        do j=0,jmax+1
        do i=0,imax+1
          anyvar2d(i,j)=q2_t(i,j,1)*mask_t(i,j,kmax)  &
                         -9999.*(1.-mask_t(i,j,kmax))
        enddo
        enddo
       endif                  !--------->
       texte80(1)='q2_t' ; texte80(2)='kg/kg'
       call netcdf_main('_t')
      endif                       !*****>

      if(allocated(q2delta_t)) then    !*****>
       if(loop_netcdf==1) then !---------> !15-12-11
        do j=1,jmax  
        do i=1,imax  
          anyvar2d(i,j)=q2delta_t(i,j,1)
        enddo
        enddo
       endif                  !--------->
       texte80(1)='q2delta_t' ; texte80(2)='kg/kg'
       call netcdf_main('_t')
      endif                           !*****>

      if(allocated(ablheight_t)) then    !*****>
       if(loop_netcdf==1) then !---------> !15-12-11
        do j=0,jmax+1
        do i=0,imax+1
          anyvar2d(i,j)=ablheight_t(i,j,1)
        enddo
        enddo
       endif                  !--------->
       texte80(1)='ablheight_t' ; texte80(2)='m'
       call netcdf_main('_t')
      endif                           !*****>

      endif                     !pmxpmx>

!******************************************************************************


!******************************************************************************
!     Impression maree
      if (kmaxtide>0) then
!******************************************************************************

      if(loop2_.eq.2)then           !2222222>
! ha
      do ktide=1,kmaxtide
        if(loop_netcdf==1) then !-------->
         do j=0,jmax+1 !04-06-15
         do i=0,imax+1
           anyvar2d(i,j)=                                       &
           sqrt(sshtidecos_w(i,j,ktide)**2+                  &
                sshtidesin_w(i,j,ktide)**2)*mask_t(i,j,kmax) &
                                     -9999.*(1-mask_t(i,j,kmax))
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
         do j=0,jmax+1 !04-06-15
         do i=0,imax+1
           x1=-atan2(-sshtidesin_w(i,j,ktide),sshtidecos_w(i,j,ktide))
           anyvar2d(i,j)=mod( x1*180./pi,360.*un)*mask_t(i,j,kmax) & 
                                        -9999.*(1-mask_t(i,j,kmax))
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
         do j=0,jmax+1 !17-02-16
         do i=0,imax+1
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
         do j=0,jmax+1 !17-02-16
         do i=0,imax+1
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
        k0=0 ; if(ioffline==2)k0=1 !04-07-15
        if(iteration2d_max_now==0)k0=1
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=(                    &
                        ssh_w(i,j,2)*(1.-k0) &
                   +ssh_int_w(i,j,1)*    k0  &

!                       )*nint(wetmask_t(i,j)) & 
                        )*mask_t(i,j,kmax)     & !16-08-22
!              +filval*(1-nint(wetmask_t(i,j))) 
               +filval*(1-mask_t(i,j,kmax))      !16-08-22
        enddo ; enddo
      endif                  !=======>
      texte80(1)='ssh_w'       ; texte80(2)='m'               ! variable ; units
      texte80(3:4)='sea_surface_height_above_geoid' 
      call netcdf_main('_w')

      if(flag_bathy_update==1) then !ooo> !12-05-22
      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=h_w(i,j)*mask_t(i,j,kmax)     & !20-01-15
                      +filval*(1-mask_t(i,j,kmax)) 
        enddo ; enddo
      endif                  !=======>
      texte80(1)='h_w'       ; texte80(2)='m'               ! variable ; units
      texte80(3:4)='bathymetry' 
      call netcdf_main('_w')
      endif                         !ooo> !12-05-22

#ifdef bidon
! timeabovewater_w !02-03-21
      if(allocated(timeabovewater_w)) then !>>>>>>>
      if(loop_netcdf==1) then !=======>
        x1=86400./(elapsedtime_now-elapsedtime_rst+dti_fw)
!       x1=1./real(timeabovewater_count)
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=x1*real(timeabovewater_w(i,j))*mask_t(i,j,kmax)     & !20-01-15
                                            +filval*(1-mask_t(i,j,kmax)) 
        enddo ; enddo
      endif                  !=======>
      texte80(1)='timeabovewater_w'   ; texte80(2)='days'               ! variable ; units
      texte80(3:4)='timeabovewater_w' 
      call netcdf_main('_w')
      endif                                !>>>>>>>
#endif

!     if(flag_nh2d==1) then !pmx>
!     if(allocated(sshr4_w)) then !>>>
!     if(loop_netcdf==1) then !=======>
!     anyvar2d(0     ,:)=filval
!     anyvar2d(imax+1,:)=filval
!     anyvar2d(:,0     )=filval
!     anyvar2d(:,jmax+1)=filval
!       do j=1,jmax ; do i=1,imax
!         anyvar2d(i,j)=sshr4_w(i,j)*mask_t(i,j,kmax) &
!                         +filval*(1-mask_t(i,j,kmax))
!       enddo ; enddo
!     endif                  !=======>
!     texte80(1)='sshr4'       ; texte80(2)='m'               ! variable ; units
!     texte80(3:4)='sshr4'
!     call netcdf_main('_w')
!     endif                       !>>>
!     endif                 !pmx>

#ifdef bidon
      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=(                            &
                max(0.,(flag_nh2d *ssh_w(i,j,2)      &
                    +(1-flag_nh2d)*ssh_int_w(i,j,2)) &
                  /( abs(h_w(i,j))+0.001 ))          &
                        )*mask_t(i,j,kmax)     & !20-01-15
               +filval*(1-mask_t(i,j,kmax)) 
        enddo ; enddo
      endif                  !=======>
      texte80(1)='sshoverh'       ; texte80(2)='none'               ! variable ; units
      texte80(3:4)='sshoverh' 
      call netcdf_main('_w')
#endif

#ifdef bidon
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval
      anyvar2d(imax+1,:)=filval
      anyvar2d(:,0     )=filval
      anyvar2d(:,jmax+1)=filval
!     if(flag_nh2d==1) then !>>>
!       do j=1,jmax ; do i=1,imax
!         anyvar2d(i,j)=sqrt(                                  &
!           (0.5*(ssh_w(i+1,j,2)-ssh_w(i-1,j,2))/dx_t(i,j))**2 &
!          +(0.5*(ssh_w(i,j+1,2)-ssh_w(i,j-1,2))/dy_t(i,j))**2 &
!                       )*mask_t(i,j,kmax)     & !20-01-15
!              +filval*(1-mask_t(i,j,kmax)) 
!       enddo ; enddo
!     else                  !>>>
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=sqrt(                                          &
            (0.5*(ssh_int_w(i+1,j,2)-ssh_int_w(i-1,j,2))/dx_t(i,j))**2 &
           +(0.5*(ssh_int_w(i,j+1,2)-ssh_int_w(i,j-1,2))/dy_t(i,j))**2 &
                        )*mask_t(i,j,kmax)     & !20-01-15
               +filval*(1-mask_t(i,j,kmax)) 
        enddo ; enddo
!     endif                 !>>>
      endif                  !=======>
      texte80(1)='sshslope'       ; texte80(2)='none'               ! variable ; units
      texte80(3:4)='sshslope' 
      call netcdf_main('_w')
#endif

!     if(allocated(ssh_timeaveraged).and.iteration2d==iteration2d_max_now+1) then !pmx>
      if(allocated(ssh_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval
      anyvar2d(imax+1,:)=filval
      anyvar2d(:,0     )=filval
      anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)=ssh_timeaveraged(i,j)*mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='ssh_timeaveraged' ; texte80(2)='m'    ! variable ; units
      texte80(3:4)='ssh_timeaveraged' 
      call netcdf_main('_w')
      endif                                !pmx>

!     if(allocated(variance_timeaveraged)) then !pmx>
!     if(loop_netcdf==1) then !=======>
!     anyvar2d(0     ,:)=filval
!     anyvar2d(imax+1,:)=filval
!     anyvar2d(:,0     )=filval
!     anyvar2d(:,jmax+1)=filval
!       do j=1,jmax   ; do i=1,imax  
!         anyvar2d(i,j)=variance_timeaveraged(i,j)*mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
!       enddo ; enddo
!     endif                  !=======>
!     texte80(1)='ecarttypefois4' ; texte80(2)='m'    ! variable ; units
!     texte80(3:4)='significantheight' 
!     call netcdf_main('_w')
!     endif                                !pmx>

      if(allocated(sshmax_w)) then !pmx>
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0     )=filval ; anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax
          anyvar2d(i,j)=sshmax_w(i,j)*mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='sshmax_w' ; texte80(2)='m' ; texte80(3:4)='sshmax_w'
      call netcdf_main('_w')
      endif                       !pmx>

      if(allocated(sshmin_w)) then !pmx>
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0     )=filval ; anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax
          anyvar2d(i,j)=sshmin_w(i,j)*mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='sshmin_w' ; texte80(2)='m' ; texte80(3:4)='sshmin_w'
      call netcdf_main('_w')
      endif                        !pmx>

      if(allocated(sshlwf_w)) then !pmx>
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0     )=filval ; anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax
          anyvar2d(i,j)=sshlwf_w(i,j,1)*mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='sshlwf_w' ; texte80(2)='m' ; texte80(3:4)='sshlwf_w'
      call netcdf_main('_w')
      endif                        !pmx>

      if(allocated(variance_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0     )=filval ; anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax
          anyvar2d(i,j)=4.*sqrt(variance_timeaveraged(i,j)) &
           *mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='HS' ; texte80(2)='m' 
      texte80(3:4)='hauteur_significative'
      call netcdf_main('_w')
      endif                        !pmx>

      if(allocated(breaker2d_t)) then !pmx>
      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0     )=filval ; anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax
          anyvar2d(i,j)=breaker2d_t(i,j)*mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='breaker2d' ; texte80(2)='?' ; texte80(3:4)=texte80(1)
      call netcdf_main('_w')
      endif                           !pmx>

      if(loop_netcdf==1) then !=======>
      anyvar2d(0     ,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0     )=filval ; anyvar2d(:,jmax+1)=filval
        do j=1,jmax   ; do i=1,imax
          anyvar2d(i,j)=                                                       & 
         sqrt( (0.5*(ssh_int_w(i+1,j,1)-ssh_int_w(i-1,j,1))*invdx_t(i,j))**2   &
              +(0.5*(ssh_int_w(i,j+1,1)-ssh_int_w(i,j-1,1))*invdy_t(i,j))**2 ) &
                        *mask_t(i,j,kmax)+filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
      endif                  !=======>
      texte80(1)='sshslope' ; texte80(2)='?' ; texte80(3:4)=texte80(1)
      call netcdf_main('_w')

#ifdef bidon
      if(loop_netcdf==1) then !=======>
        if(.not.allocated(q_t))call q_allocate
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=q2davr_w(i,j)
!         anyvar2d(i,j)=(                    &
!                      q2davr_w(i,j)         &
!                       )*mask_t(i,j,kmax)   & !20-01-15
!              +filval*(1-mask_t(i,j,kmax)) 
        enddo ; enddo
      endif                  !=======>
      texte80(1)='q2davr'       ; texte80(2)='m'               ! variable ; units
      texte80(3:4)='q2davr' 
      call netcdf_main('_w')
#endif

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression du niveau de fond de la grille
!******************************************************************************

      if(loop2_==7)then !------->

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

      if(allocated(kmerged_t)) then !pmx> !03-05-17
       if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
          anyvar2d(i,j)=kmerged_t(i,j)-kmin_w(i,j)
         else
          anyvar2d(i,j)=filval
         endif
        enddo
        enddo
       endif                  !=======>
       texte80(1)='kmerged_minus_kmin_w' ; texte80(2)='index'
       call netcdf_main('_w')
      endif                         !pmx> !03-05-17

      endif            !------->
!******************************************************************************


!******************************************************************************
! impression de la zone eponge du mode interne
!******************************************************************************
      if(loop2_.eq.8)then !------->

       if(loop_netcdf==1) then !=======>
         do j=1,jmax  ; do i=1,imax  
          anyvar2d(i,j)=                       &
              mask_t(i,j,kmax)*sponge_t(i,j,1) &
          +(1-mask_t(i,j,kmax))*filval
         enddo ; enddo
         anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
         anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
       endif                  !=======>
       texte80(1)='sponge_t' ; texte80(2)='none'
       texte80(3:4)='dimensionless_sponge_t' 
       call netcdf_main('_t')

       if(loop_netcdf==1) then !=======>
         do j=1,jmax  ; do i=1,imax+1
          anyvar2d(i,j)=                         &
              mask_u(i,j,kmax)*sponge_u(i,j,1)/max(relax_int,small2) & !21-10-22
!             mask_u(i,j,kmax)/max(sponge_u(i,j,1),small2)/86400. &
          +(1-mask_u(i,j,kmax))*filval
         enddo ; enddo
       endif                  !=======>
       texte80(1)='sponge_u' ; texte80(2)='days'
       texte80(3:4)='sponge_u' 
       call netcdf_main('_u')

       if(loop_netcdf==1) then !=======>
         do j=1,jmax+1  ; do i=1,imax
          anyvar2d(i,j)=                       &
              mask_v(i,j,kmax)*sponge_v(i,j,1)/max(relax_int,small2) & !21-10-22
!             mask_v(i,j,kmax)/max(sponge_v(i,j,1),small2)/86400. &
          +(1-mask_v(i,j,kmax))*filval
         enddo ; enddo
       endif                  !=======>
       texte80(1)='sponge_v' ; texte80(2)='days'
       texte80(3:4)='sponge_v' 
       call netcdf_main('_v')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de l'elevation de la surface forcante
!******************************************************************************
      if(loop2_==9)then !------->

      if(loop_netcdf==1) then !=======>
      anyvar2d=filval ! reset pour les 3 champs A suivre
      do j=1,jmax  
      do i=1,imax  
       anyvar2d(i,j)=                                                &
       (    timeweightobc(ssh_id) *sshobc_w(i,j,2)                   &  !11-07-14
       +(1.-timeweightobc(ssh_id))*sshobc_w(i,j,0))*mask_t(i,j,kmax) &
                                                +(1-mask_t(i,j,kmax))*filval
      enddo
      enddo
      endif                  !=======>
      texte80(1)='ssh_ogcm' ; texte80(2)='m'
      call netcdf_main('_w')

      if(kmaxtide>0) then !--seulement-si-maree->

      if(loop_netcdf==1) then !=======>
      anyvar2d=filval
      do ktide=1,kmaxtide
       time1=frqtide(ktide)*((elapsedtime_aft-ti0tide(ktide))+0.5*dte_lp)+v0tide(ktide)+utide(ktide,1)
       const3=cos(time1)*rampe*ftide(ktide,1)  
       const4=sin(time1)*rampe*ftide(ktide,1)    
       do j=1,jmax ; do i=1,imax
        anyvar2d(i,j)=anyvar2d(i,j)*passetide(ktide,1)+(sshtidecos_w(i,j,ktide)*const3+sshtidesin_w(i,j,ktide)*const4) 
       enddo  ; enddo
      enddo
      do j=1,jmax ; do i=1,imax
       anyvar2d(i,j)=anyvar2d(i,j)*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
      enddo  ; enddo
      endif                  !=======>
      texte80(1)='ssh_tide' ; texte80(2)='m'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !=======>
      do j=1,jmax  
      do i=1,imax  
       anyvar2d(i,j)=(                                         &
       anyvar2d(i,j)+ & ! ssh tide de l'archivage precedent
              timeweightobc(ssh_id) *sshobc_w(i,j,2)           &
         +(1.-timeweightobc(ssh_id))*sshobc_w(i,j,0))*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
      enddo
      enddo
      endif                   !=======>
      texte80(1)='ssh_ogcm+tide' ; texte80(2)='m'
      call netcdf_main('_w')

      endif               !--seulement-si-maree->

      endif           !------->


!******************************************************************************
! impression de la zone de relaxation champs grANDe echelle methode MPV
!******************************************************************************
      if(loop2_.eq.10)then!------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1
        do i=0,imax+1
         if(mask_t(i,j,kmax  )==1) then
!         anyvar2d(i,j)=ssh_w(i,j,2)+h_w(i,j)
          anyvar2d(i,j)=ssh_int_w(i,j,2)+h_w(i,j)
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
! PAR2 and albedo !24-05-15
!******************************************************************************
      if(loop2_==11)then!------->

      if(loop_netcdf==1) then !=======>
       anyvar2d(0,:)=filval ; anyvar2d(:,0)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,jmax+1)=filval
       do j=1,jmax ; do i=1,imax
         anyvar2d(i,j)=mask_t(i,j,kmax)*(1./light_kpar2_w(i,j)) &
                   +(1-mask_t(i,j,kmax))*filval
       enddo ; enddo
      endif                  !=======>
      texte80(1)='PAR2' ; texte80(2)='m'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
       anyvar2d(0,:)=filval ; anyvar2d(:,0)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,jmax+1)=filval
       do j=1,jmax ; do i=1,imax
         anyvar2d(i,j)=mask_t(i,j,kmax)*albedo_w(i,j) &
                   +(1-mask_t(i,j,kmax))*filval
       enddo ; enddo
      endif                  !=======>
      texte80(1)='albedo_w' ; texte80(2)='ratio'
      call netcdf_main('_t')

      endif         !------->
!******************************************************************************


!******************************************************************************
! impression elevation surface - elevation surface forcante
!******************************************************************************
      if(loop2_.eq.12)then!------->
      if(loop_netcdf==1) then !=======>
      do j=1,jmax
      do i=1,imax
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
! C.L. pour visu cas baroclinic_jet
      anyvar2d(0       ,1:jmax)=anyvar2d(1       ,1:jmax)
      anyvar2d(imax+1  ,1:jmax)=anyvar2d(imax    ,1:jmax)
      anyvar2d(0:imax+1,0     )=anyvar2d(0:imax+1,1     )      
      anyvar2d(0:imax+1,jmax+1)=anyvar2d(0:imax+1,jmax  )
      endif                  !=======>
      texte80(1)='ssh_sshobc' ; texte80(2)='m'
      call netcdf_main('_t')

      endif           !------->

!******************************************************************************
! impression des bancs decouvrants:
!******************************************************************************
      if(loop2_==13)then!------->

       if(loop_netcdf==1) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax)*upwindwetdry_t(i,j) & !28-02-22
                   +(1-mask_t(i,j,kmax))*filval
       enddo ; enddo
       endif                  !=======>
       texte80(1)='upwindwetdry_t' ; texte80(2)='0=upwind'
       call netcdf_main('_t')

       if(loop_netcdf==1) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax)*wetmask_t(i,j) &
                   +(1-mask_t(i,j,kmax))*filval
       enddo ; enddo
       endif                  !=======>
       texte80(1)='wetmask_t' ; texte80(2)='0=dry_1=wet'
       call netcdf_main('_t')

       if(loop_netcdf==1) then !=======>
       do j=1,jmax   ; do i=1,imax  
         anyvar2d(i,j)=mask_t(i,j,kmax)*wetmask_wi_t(i,j) &
                   +(1-mask_t(i,j,kmax))*filval
       enddo ; enddo
       endif                  !=======>
       texte80(1)='wetmask_wi_t' ; texte80(2)='0=dry_1=wet'
       call netcdf_main('_t')

       if(loop_netcdf==1) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_u(i,j,kmax)*wetmask_u(i,j) &
                   +(1-mask_u(i,j,kmax))*filval
       enddo ; enddo
       endif                  !=======>
       texte80(1)='wetmask_u' ; texte80(2)='0=dry_1=wet' !16-07-16
       call netcdf_main('_u')

       if(loop_netcdf==1) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_v(i,j,kmax)*wetmask_v(i,j) &
                   +(1-mask_v(i,j,kmax))*filval
       enddo ; enddo
       endif                  !=======>
       texte80(1)='wetmask_v' ; texte80(2)='0=dry_1=wet'
       call netcdf_main('_v')

      endif         !------->

!******************************************************************************
! impression du barometre inverse:
!******************************************************************************
      if(loop2_==14)then!------->
      if(loop_netcdf==1) then !=======>
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

      if(par%rank==0)write(6,*)'impression des parametres des vagues'

      if(loop_netcdf==1.and.allocated(hs_wave_t)) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=hs_wave_t(i,j,1) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='hs_wave_t' ; texte80(2)='m'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(ubw)) then !=======>
       lb2=lbound(ubw) ; ub2=ubound(ubw) ; i1=lb2(1) ; i2=ub2(1) ; j1=lb2(2) ; j2=ub2(2)
       do j=j1,j2 ; do i=i1,i2
         anyvar2d(i,j)=ubw(i,j) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='orbital_uv' ; texte80(2)='m/s'
      texte80(11)='mask_obc_z1'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(ubw).and.allocated(t_wave_t))then !=======>
       lb2=lbound(ubw) ; ub2=ubound(ubw) ; i1=lb2(1) ; i2=ub2(1) ; j1=lb2(2) ; j2=ub2(2)
       do j=j1,j2 ; do i=i1,i2
         anyvar2d(i,j)=ubw(i,j)*t_wave_t(i,j,1)/(2.*pi) &
                   *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='demiexcursion' ; texte80(2)='m'
      texte80(11)='mask_obc_z1'
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(t_wave_t)) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=t_wave_t(i,j,1) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='t_wave_t' ; texte80(2)='s'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(dir_wave_t)) then !=======>    !07-06-10
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mod(270.-dir_wave_t(i,j,1)*rad2deg             &
        +atan2(gridrotsin_t(i,j),gridrotcos_t(i,j))*rad2deg           &
!       -atan2(  lat_t(i+1,j)-lat_t(i-1,j)                            & ! rotation symphonie
!              ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j)) )*rad2deg &
                          ,un*360.d0)                                 &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='dir_wave_conv1' ; texte80(2)='deg' !25-04-17
      texte80(3:4)='ww3_convention_0fromN_90fromE' 
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(dir_wave_t)) then !=======>    !07-06-10
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mod(     dir_wave_t(i,j,1)*rad2deg             &
        -atan2(gridrotsin_t(i,j),gridrotcos_t(i,j))*rad2deg           &
!       +atan2(  lat_t(i+1,j)-lat_t(i-1,j)                            & ! rotation symphonie
!              ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j)) )*rad2deg &
                          ,un*360.d0)                                 &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='dir_wave_conv2' ; texte80(2)='deg' !25-04-17
      texte80(3:4)='S_convention_0=O>E_90=S>N' 
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(dir_wave_t)) then !=======>    !07-06-10
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mod(     dir_wave_t(i,j,1)*rad2deg             &
                          ,un*360.d0)                                 &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='dir_wave_conv3' ; texte80(2)='deg' !25-04-17
      texte80(3:4)='0=along0i_90=alongOj' 
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

#ifdef bidon
      if(loop_netcdf==1.and.allocated(kx_wave_t)) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=kx_wave_t(i,j,1) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='kx_wave_t' ; texte80(2)='1/m'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(ky_wave_t)) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=ky_wave_t(i,j,1) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='ky_wave_t' ; texte80(2)='1/m'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')
#endif
      if(loop_netcdf==1.and.allocated(k_wave_t)) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=2*pi/max(k_wave_t(i,j,1),small1) & !19-02-23
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='wavelength' ; texte80(2)='m'
      texte80(11)='mask_obc_z1'
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1.and.allocated(foc_wave_t)) then !=======> !05-09-13
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=foc_wave_t(i,j,1) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='foc_wave_t' ; texte80(2)='W/m2'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=sshstokes_w(i,j) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='sshstokes_w' ; texte80(2)='m'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_w')

      if(loop_netcdf==1.and.allocated(hsw_wave_t)) then !=======> !05-09-13
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=hsw_wave_t(i,j,1) &
          *mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
       enddo
       enddo
      endif                  !=======>
      texte80(1)='hsw' ; texte80(2)='m'
      texte80(11)='mask_obc_z1' !24-01-15
      call netcdf_main('_w')
#endif
       endif         !------->

!******************************************************************************
! impression drifters 2D
!******************************************************************************
      if(loop2_==16)then!------->

      if(loop_netcdf==1) then !=======>
       do j=1,jmax ; do i=1,imax
        anyvar2d(i,j)=(1.-mask_t(i,j,kmax))*filval
       enddo       ; enddo
       do kbu=1,kbumax
        i=nint(drifter_l(kbu,1)-par%timax(1))
        j=nint(drifter_l(kbu,2)-par%tjmax(1))
        if(i>=2.and.i<=imax-1.and.j>=2.and.j<=jmax-1) then !>>>
         if(mask_t(i,j,kmax)==1)anyvar2d(i,j)=anyvar2d(i,j)+1.
        endif                                              !>>>
       
       enddo
       call obc_int_anyvar2d('za')

      endif                  !=======>
      texte80(1)='drifters2D' ; texte80(2)='nb/maille2d'
      texte80(3:4)='nbre_drifter_par_maille_horizontale'
      call netcdf_main('_w')


      endif         !------->

!******************************************************************************
! Bottom Drag Coefficient & Roughness
!******************************************************************************
      if(loop2_==17)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=mask_t(i,j,kmax) *cdb_t(i,j) &
                   +(1.-mask_t(i,j,kmax))*filval
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 
      endif                  !=======>
      texte80(1)='cdb_t'       ; texte80(2)='none'  ! variable ; units
      texte80(3)='Bottom drag coefficient'          ! long_name
      texte80(4)='Bottom_drag_coefficient'          ! standard_name
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=mask_t(i,j,kmax) *z0_w(i,j) & !24-05-15
                   +(1.-mask_t(i,j,kmax))*filval
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 
      endif                  !=======>
      texte80(1)='z0_w'       ; texte80(2)='none'  ! variable ; units
      texte80(3)='Bottom roughness'          ! long_name
      texte80(4)='Bottom_roughness'          ! standard_name
      call netcdf_main('_t')

      endif                     !------->

!******************************************************************************
! impression coriolis
!******************************************************************************
      if(loop2_==18)then        !------->

      if(loop_netcdf==1) then !=======>
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=coriolis_t(i,j)*mask_t(i,j,kmax)         &
                                    +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo
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
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
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
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
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
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
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
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
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
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'vh_m',loopm_,'tide',ktide
      texte80(2)='degrees'  ! variable ; units
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======> !14-04-16
      if(allocated(p3dmode_cos_t)) then !)))))>
        do j=1,jmax ; do i=1,imax
           anyvar2d(i,j)=                                              &
             mask_t(i,j,kmax)*sqrt(p3dmode_cos_t(i,j,loopm_,ktide)**2  &
                                  +p3dmode_sin_t(i,j,loopm_,ktide)**2) &
         +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'pa_m',loopm_,'tide',ktide
      texte80(2)='kg.m**-1.s**-2'  ! variable ; units
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======> !14-04-16
      if(allocated(p3dmode_cos_t)) then !)))))>
        do j=1,jmax ;  do i=1,imax
          anyvar2d(i,j)=filval*(1-mask_t(i,j,kmax))+mask_t(i,j,kmax)*( &!pmx>
               mod(-atan2(-p3dmode_sin_t(i,j,loopm_,ktide)             &
                          ,p3dmode_cos_t(i,j,loopm_,ktide) )*rad2deg,360.d0)   &
                                                                     )  !pmx>
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !20-02-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !20-02-16
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'ph_m',loopm_,'tide',ktide
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
        do j=1,jmax !13-09-16
        do i=1,imax
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

      if(loop2_==22)then !22-22-22-22> !18-01-16

       if(loop_netcdf==1) then !=======>
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=rhpzavr_w(i,j)*mask_t(i,j,kmax) &
                            +filval*(1-mask_t(i,j,kmax))
        enddo ; enddo
       endif                  !=======>
       texte80(1)='rhpzavr_w' ; texte80(2)='kg/m3'
       texte80(3)='z_averaged_sea_water_potential_density'
       texte80(4)='z_averaged_sea_water_potential_density'
       call netcdf_main('_w')

      endif              !22-22-22-22>

      if(loop2_==23)then !23-23-23-23> !18-01-16

       if(loop_netcdf==1) then !=======>
        k=kmax
        do j=1,jmax+1 ; do i=1,imax+1
          anyvar2d(i,j)=                                                 &
         mask_f(i,j,kmax)*( (vel_v(i,j,k,0)-vel_v(i-1,j,k,0))/dx_f(i,j)  &
                           -(vel_u(i,j,k,0)-vel_u(i,j-1,k,0))/dy_f(i,j)) &
         +(1-mask_f(i,j,kmax))*filval
        enddo ; enddo
       endif                  !=======>
       texte80(1)='surf_rel_vorticity' ; texte80(2)='s**-1'
       texte80(3)='Near_surface_relative_vorticity'
       texte80(4)='Near_surface_relative_vorticity'
       call netcdf_main('_f')

      endif              !23-23-23-23>

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

       if(allocated(wstresb_u)) then !ooooo> !28-09-16
         if(loop_netcdf==1) then !pmxpmx>
          do j=1,jmax ; do i=1,imax+1
           anyvar2d(i,j)=mask_u(i,j,kmax)*wstresb_u(i,j)+(1.-mask_u(i,j,kmax))*filval
          enddo       ; enddo
          call obc_int_anyvar2d('u1') 
         endif                   !pmxpmx>
         texte80(1)='wstresb_u'  ; texte80(2)='N/m2'
         call netcdf_main('_u')
       endif                         !ooooo>

       if(allocated(wstresb_v)) then !ooooo>
         if(loop_netcdf==1) then !pmxpmx>
          do j=1,jmax+1 ; do i=1,imax
           anyvar2d(i,j)=mask_v(i,j,kmax)*wstresb_v(i,j)+(1.-mask_v(i,j,kmax))*filval
          enddo       ; enddo
          call obc_int_anyvar2d('v1') 
         endif                   !pmxpmx>
         texte80(1)='wstresb_v'  ; texte80(2)='N/m2'
         call netcdf_main('_v')
       endif                         !ooooo>


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
!     anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
!     anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12
      call obc_int_anyvar2d('u1') !20-02-15 !16-03-15

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
!     anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
!     anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12
      call obc_int_anyvar2d('v1') !20-02-15!16-03-15

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
!     anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
!     anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12
      call obc_int_anyvar2d('u1') !20-02-15!16-03-15

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
!     anyvar2d(1,:)=anyvar2d(2,:)             !06-01-12
!     anyvar2d(:,1)=anyvar2d(:,2)             !06-01-12
      call obc_int_anyvar2d('v1') !20-02-15!16-03-15
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

      if(allocated(uwind_t)) then !=======> !20-10-14

       if(loop_netcdf==1) then !//////////////>
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=uwind_t(i,j,1)
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !17-03-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !17-03-16
       endif                   !//////////////>
      texte80(1)='uwind_t'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      endif                       !=======>

      if(allocated(uwinddelta_t)) then !=======> !20-10-14

       if(loop_netcdf==1) then !//////////////>
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=uwinddelta_t(i,j,1)
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !17-03-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !17-03-16
       endif                   !//////////////>
      texte80(1)='uwinddelta_t'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      endif                       !=======>

      if(allocated(uwind100_t)) then !=======> !20-10-14
      if(allocated(airseafile_prvtime)) then !*****>
       if(loop_netcdf==1) then !//////////////>

        x2=(elapsedtime_now           -airseafile_prvtime(t2m_id))      &
          /(airseafile_nextime(t2m_id)-airseafile_prvtime(t2m_id))
        x2=min(max(x2,0.d00),1.d00) 
        x1=1.-x2
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=x1*uwind100_t(i,j,1)+x2*uwind100_t(i,j,2)
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !17-03-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !17-03-16

       endif                   !//////////////>
      texte80(1)='uwind100_t'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      endif                                  !*****>
      endif                       !=======>


      if(allocated(vwind_t)) then !=======> !20-10-14

       if(loop_netcdf==1) then !//////////////>
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=vwind_t(i,j,1)
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !17-03-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !17-03-16
       endif                   !//////////////>
      texte80(1)='vwind_t'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      endif                       !=======>

      if(allocated(vwinddelta_t)) then !=======> !20-10-14

       if(loop_netcdf==1) then !//////////////>
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=vwinddelta_t(i,j,1)
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !17-03-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !17-03-16
       endif                   !//////////////>
      texte80(1)='vwinddelta_t'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      endif                       !=======>

      if(allocated(vwind100_t)) then !=======> !20-10-14
      if(allocated(airseafile_prvtime)) then !*****>
       if(loop_netcdf==1) then !//////////////>

        x2=(elapsedtime_now           -airseafile_prvtime(t2m_id))      &
          /(airseafile_nextime(t2m_id)-airseafile_prvtime(t2m_id))
        x2=min(max(x2,0.d00),1.d00) 
        x1=1.-x2
        do j=1,jmax ; do i=1,imax
          anyvar2d(i,j)=x1*vwind100_t(i,j,1)+x2*vwind100_t(i,j,2)
        enddo ; enddo
        anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 !17-03-16
        anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 !17-03-16

       endif                   !//////////////>
      texte80(1)='vwind100_t'  ; texte80(2)='m/s'
      call netcdf_main('_t')

      endif                                  !*****>
      endif                       !=======>

      endif                     !------->

!******************************************************************************


!******************************************************************************
! impression du courant moyen sur la colonne d'eau
!******************************************************************************
      if(loop2_.eq.2)then !2222222>

      if(loop_netcdf==1) then !=======>
          do j=1,jmax ; do i=1,imax+1
           if(mask_u(i,j,kmax  )==1) then
            anyvar2d(i,j)=velbar_u(i,j,1)+velbarstokes_u(i,j,1)   !25-02-10
! converti en nombre de courant: !20-08-14
!           anyvar2d(i,j)=anyvar2d(i,j)*dte_lp/dx_u(i,j)
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo ; enddo
      endif                  !=======>
      texte80(1)='velbar_u' ; texte80(2)='m/s'
      texte80(3:4)='along Oi axis depth-averaged total current'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
          do j=1,jmax+1
          do i=1,imax
           if(mask_v(i,j,kmax  )==1) then
            anyvar2d(i,j)=velbar_v(i,j,1)+velbarstokes_v(i,j,1)   !25-02-10
! converti en nombre de courant: !20-08-14
!           anyvar2d(i,j)=anyvar2d(i,j)*dte_lp/dy_v(i,j)
           else
            anyvar2d(i,j)=-9999.
           endif
          enddo
          enddo
      endif                  !=======>
      texte80(1)='velbar_v' ; texte80(2)='m/s'
      texte80(3:4)='along Oj axis depth-averaged total current'
      call netcdf_main('_v')

      if(allocated(flx2d_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
       do j=1,jmax ; do i=1,imax+1
        anyvar2d(i,j)=flx2d_timeaveraged(i,j)*invdy_u(i,j)*mask_u(i,j,kmax) &
                                                       +(1-mask_u(i,j,kmax))*filval
       enddo ; enddo
      endif                  !=======>
      texte80(1)='flx2d_timeaveraged' 
      texte80(2)='m**2/s'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
       do j=1,jmax+1 ; do i=1,imax
        anyvar2d(i,j)=fly2d_timeaveraged(i,j)*invdx_v(i,j)*mask_v(i,j,kmax) &
                                                       +(1-mask_v(i,j,kmax))*filval
       enddo ; enddo
      endif                  !=======>
      texte80(1)='fly2d_timeaveraged' 
      texte80(2)='m**2/s'
      call netcdf_main('_v')

      endif                                                                       !pmx>

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
           anyvar2d(i,j)=sqrt(veltidecos_u(i,j,ktide)**2+             &
                              veltidesin_u(i,j,ktide)**2)
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
           anyvar2d(i,j)=sqrt(veltidecos_v(i,j,ktide)**2+             &
                              veltidesin_v(i,j,ktide)**2)
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
           x1=-atan2(-veltidesin_u(i,j,ktide),veltidecos_u(i,j,ktide))
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
           x1=-atan2(-veltidesin_v(i,j,ktide),veltidecos_v(i,j,ktide))
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
      anyvar2d=filval ! reset pour les 3 champs A suivre
      do j=1,jmax ; do i=1,imax+1
       anyvar2d(i,j)=                                                 &
       (     timeweightobc(vel_id) *velbarobc_u(i,j,2)                &
        +(1.-timeweightobc(vel_id))*velbarobc_u(i,j,0))*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
      enddo       ; enddo
      endif                  !=======>
      texte80(1)='velbar_ogcm_u'  ; texte80(2)='m/s'
      call netcdf_main('_u')

      if(kmaxtide>0) then !--seulement-si-maree->

      if(loop_netcdf==1) then !=======>
      do ktide=1,kmaxtide
       time1=frqtide(ktide)*((elapsedtime_aft-ti0tide(ktide))+0.5*dte_lp)+v0tide(ktide)+utide(ktide,1)
       const3=cos(time1)*rampe*ftide(ktide,1)  
       const4=sin(time1)*rampe*ftide(ktide,1)    
       do j=1,jmax ; do i=1,imax+1
        anyvar2d(i,j)=anyvar2d(i,j)*passetide(ktide,1)+(veltidecos_u(i,j,ktide)*const3+veltidesin_u(i,j,ktide)*const4) 
       enddo  ; enddo
      enddo
      do j=1,jmax ; do i=1,imax+1
       anyvar2d(i,j)=anyvar2d(i,j)*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
      enddo  ; enddo
      endif                  !=======>
      texte80(1)='velbar_tide_u'  ; texte80(2)='m/s'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
      do j=1,jmax ; do i=1,imax+1
       anyvar2d(i,j)=(                                                &
       anyvar2d(i,j)+ & !velbar tide de l'archivage precedent
             timeweightobc(vel_id) *velbarobc_u(i,j,2)                &
        +(1.-timeweightobc(vel_id))*velbarobc_u(i,j,0))*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
      enddo       ; enddo
      endif                  !=======>
      texte80(1)='velbar_ogcm+tide_u'  ; texte80(2)='m/s'
      call netcdf_main('_u')

      endif               !--seulement-si-maree->

      if(loop_netcdf==1) then !=======>
      anyvar2d=filval ! reset pour les 3 champs A suivre
      do j=1,jmax+1 ; do i=1,imax
      anyvar2d(i,j)=    &
      (     timeweightobc(vel_id) *velbarobc_v(i,j,2)   &
       +(1.-timeweightobc(vel_id))*velbarobc_v(i,j,0))*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
      enddo ; enddo
      endif                  !=======>
      texte80(1)='velbar_ogcm_v' ; texte80(2)='m/s'
      call netcdf_main('_v')

      if(kmaxtide>0) then !--seulement-si-maree->

      if(loop_netcdf==1) then !=======>
      anyvar2d=filval
      do ktide=1,kmaxtide
       time1=frqtide(ktide)*((elapsedtime_aft-ti0tide(ktide))+0.5*dte_lp)+v0tide(ktide)+utide(ktide,1)
       const3=cos(time1)*rampe*ftide(ktide,1)  
       const4=sin(time1)*rampe*ftide(ktide,1)    
       do j=1,jmax+1 ; do i=1,imax
        anyvar2d(i,j)=anyvar2d(i,j)*passetide(ktide,1)+(veltidecos_v(i,j,ktide)*const3+veltidesin_v(i,j,ktide)*const4) 
       enddo  ; enddo
      enddo
      do j=1,jmax+1 ; do i=1,imax
       anyvar2d(i,j)=anyvar2d(i,j)*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
      enddo  ; enddo
      endif                  !=======>
      texte80(1)='velbar_tide_v' ; texte80(2)='m/s'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !=======>
      do j=1,jmax+1 ; do i=1,imax
      anyvar2d(i,j)=(   &
      anyvar2d(i,j)+    & ! velbar tide de l'archivage precedent
            timeweightobc(vel_id) *velbarobc_v(i,j,2)   &
       +(1.-timeweightobc(vel_id))*velbarobc_v(i,j,0))*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
      enddo ; enddo
      endif                  !=======>
      texte80(1)='velbar_ogcm+tide_v' ; texte80(2)='m/s'
      call netcdf_main('_v')

      endif                  !------->

      endif               !--seulement-si-maree->

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

      if(loop2_==9)then        !------->

      if(loop_netcdf==1) then !=======>
        anyvar2d=-9999.
        do j=2,jmax-1 ; do i=2,imax
          anyvar2d(i,j)=mask_u(i,j,kmax)*pres3d2d_u(i,j)/dx_u(i,j) &
                    -(1-mask_u(i,j,kmax))*9999.
        enddo ; enddo
      endif                  !=======>
      texte80(1)='pres3d2duoverdx' ; texte80(2)='??'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        anyvar2d=-9999.
        do j=2,jmax ; do i=2,imax-1
          anyvar2d(i,j)=mask_v(i,j,kmax)*pres3d2d_v(i,j)/dy_v(i,j) &
                    -(1-mask_v(i,j,kmax))*9999.
        enddo ; enddo
      endif                  !=======>
      texte80(1)='pres3d2dvoverdy' ; texte80(2)='??'
      call netcdf_main('_v')

      endif                    !------->

!******************************************************************************
! Surface total current                                               !19-11-20
!******************************************************************************
      if(loop2_==10)then        !------->

      if(loop_netcdf==1) then !=======>
       k=kmax
       do j=1,jmax ; do i=1,imax+1
        anyvar2d(i,j)=(       vel_u(i,j,kmax,1)                   &
                       +velstokes_u(i,j,kmax,1))*mask_u(i,j,kmax) & 
                                       -9999.*(1-mask_u(i,j,kmax))
       enddo ; enddo
      endif                  !=======>
      texte80(1)='surfvel_u' ; texte80(2)='m/s'
      texte80(3:4)='along Oi axis surface total current'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
       k=kmax
       do j=1,jmax+1 ; do i=1,imax
        anyvar2d(i,j)=(       vel_v(i,j,kmax,1)                   &
                       +velstokes_v(i,j,kmax,1))*mask_v(i,j,kmax) & 
                                       -9999.*(1-mask_v(i,j,kmax))
       enddo ; enddo
      endif                  !=======>
      texte80(1)='surfvel_v' ; texte80(2)='m/s'
      texte80(3:4)='along Oj axis surface total current'
      call netcdf_main('_v')

      endif                    !------->

!******************************************************************************
! courant de fond pour la coordonnee VQS                              !19-11-20
!******************************************************************************
      if(loop2_==11)then !11-11-11>

      if(loop_netcdf==1) then !=======>
       do j=1,jmax ; do i=1,imax+1
         anyvar2d(i,j)=velbot_u(i,j)*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
       enddo ; enddo
      endif                  !=======>
      texte80(1)='botvel_u' ; texte80(2)='m/s'
      texte80(3:4)='along Oi axis bottom current'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
       do j=1,jmax+1 ; do i=1,imax
         anyvar2d(i,j)=velbot_v(i,j)*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
       enddo ; enddo
      endif                  !=======>
      texte80(1)='botvel_v' ; texte80(2)='m/s'
      texte80(3:4)='along Oj axis bottom current'
      call netcdf_main('_v')

      endif              !11-11-11>

!******************************************************************************
! courant de maree                                                    !19-11-20
!******************************************************************************
      if(loop2_==12)then        !12-12-12>

!... utide
      if(loop_netcdf==1) then !=======>
      do ktide=1,kmaxtide 
       time1=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide)) & 
             +v0tide(ktide)+utide(ktide,1)                
       const3=cos(time1)*ftide(ktide,1)               
       const4=sin(time1)*ftide(ktide,1)            
 
       do j=1,jmax ; do i=1,imax+1
 
        anyvar2d(i,j)=                                          &
       (anyvar2d(i,j)*passetide(ktide,1)+                       &
          (veltidecosout_u(i,j,ktide)*const3+                   & 
           veltidesinout_u(i,j,ktide)*const4))*mask_u(i,j,kmax) &
                                    +filval*(1-mask_u(i,j,kmax))
       enddo ; enddo
      enddo ! ktide
      endif                  !=======>
      texte80(1)='tide_u'       ; texte80(2)='m/s'   ! variable ; units
      texte80(3:4)='along Oi axis tidal current'
      call netcdf_main('_u')
 
!... vtide
      if(loop_netcdf==1) then !=======>
      do ktide=1,kmaxtide 
       time1=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide)) & 
             +v0tide(ktide)+utide(ktide,1)                
       const3=cos(time1)*ftide(ktide,1)               
       const4=sin(time1)*ftide(ktide,1)            

       do j=1,jmax+1 ; do i=1,imax
 
         anyvar2d(i,j)=                                          &
        (anyvar2d(i,j)*passetide(ktide,1)+                       & 
           (veltidecosout_v(i,j,ktide)*const3+                   & 
            veltidesinout_v(i,j,ktide)*const4))*mask_v(i,j,kmax) &
                                     +filval*(1-mask_v(i,j,kmax))
 
       enddo ; enddo
      enddo ! ktide
      endif                  !=======>
      texte80(1)='tide_v'       ; texte80(2)='m/s'   ! variable ; units
      texte80(3:4)='along Oj axis tidal current'
      call netcdf_main('_v')

      endif                     !12-12-12>

!******************************************************************************
! Surface Stokes drift                                               !19-11-20
!******************************************************************************
      if(loop2_==13)then        !13-13-13>

      if(loop_netcdf==1) then !=======>
       k=kmax
       do j=1,jmax ; do i=1,imax+1
        anyvar2d(i,j)=(                                           &
                       +velstokes_u(i,j,kmax,1))*mask_u(i,j,kmax) & 
                                       -9999.*(1-mask_u(i,j,kmax))
       enddo ; enddo
      endif                  !=======>
      texte80(1)='surfstokes_u' ; texte80(2)='m/s'
      texte80(3:4)='along Oi axis surface Stokes drift'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
       k=kmax
       do j=1,jmax+1 ; do i=1,imax
        anyvar2d(i,j)=(                                           &
                       +velstokes_v(i,j,kmax,1))*mask_v(i,j,kmax) & 
                                       -9999.*(1-mask_v(i,j,kmax))
       enddo ; enddo
      endif                  !=======>
      texte80(1)='surfstokes_v' ; texte80(2)='m/s'
      texte80(3:4)='along Oj axis surface Stokes drift'
      call netcdf_main('_v')

      endif                     !13-13-13>

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
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=tem_t(i,j,k,1)*mask_t(i,j,k) &  !20-01-17
                               -9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='tem'   ; texte80(2)='degrees_Celsius'
      texte80(3)='sea_water_potential_temperature'
      texte80(4)='sea_water_potential_temperature'
      call netcdf_main('_t')

!     if(allocated(temf_t).and.ratio_negdif_ver>0.) then !>>> !30-09-21
!     if(loop_netcdf==1) then !=======>
!       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
!         anyvar3d(i,j,k)=temf_t(i,j,k)/ratio_negdif_ver*mask_t(i,j,k)-9999.*(1-mask_t(i,j,k))
!       enddo ; enddo ; enddo
!     endif                  !=======>
!     texte80(1)='temf_vert' ; texte80(2)='degrees_Celsius' ; texte80(3:4)=texte80(1)
!     call netcdf_main('_t')
!     endif                      !>>> !19-03-21

      if(allocated(temlwf_t)) then !>>> !31-10-22
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyvar3d(i,j,k)=temlwf_t(i,j,k)*mask_t(i,j,k) &  !20-01-17
                                -9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='temlwf_t' ; texte80(2)='degrees_Celsius' ; texte80(3:4)=texte80(1)
      call netcdf_main('_t')
      endif                        !>>>

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de la salinite
!******************************************************************************
      if(loop2_.eq.2)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=sal_t(i,j,k,2)*mask_t(i,j,k) & !20-01-17
                               -9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>

      texte80(1)='sal'       ; texte80(2)='1e-3'                  ! variable ; units
      texte80(3)='sea water salinity'                             ! long_name
      texte80(4)='sea_water_salinity'                             ! standard_name
      call netcdf_main('_t')

!     if(allocated(salf_t).and.ratio_negdif_ver>0.) then !>>> !30-09-21
!     if(loop_netcdf==1) then !=======>
!       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
!         anyvar3d(i,j,k)=salf_t(i,j,k)/ratio_negdif_ver*mask_t(i,j,k)-9999.*(1-mask_t(i,j,k))
!       enddo ; enddo ; enddo
!     endif                  !=======>
!     texte80(1)='salf_vert' ; texte80(2)='1e-3' ; texte80(3:4)=texte80(1)
!     call netcdf_main('_t')
!     endif                      !>>> !19-03-21

      if(allocated(temlwf_t)) then !>>> !31-10-22
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyvar3d(i,j,k)=sallwf_t(i,j,k)*mask_t(i,j,k) &  !20-01-17
                                -9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='salwf_t' ; texte80(2)='psu' ; texte80(3:4)=texte80(1)
      call netcdf_main('_t')
      endif                        !>>>

#ifdef bilan_s_3d
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
          sal_t(i,j,k,2)-dxdydzs0(i,j,k)/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dsal'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (flus(i+1,j,k)-flus(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dflusdi'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (flvs(i,j+1,k)-flvs(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dflvsdj'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (flws(i,j,k+1)-flws(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dflwsdk'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (fimplis(i,j,k+1)-fimplis(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dfimplisdk'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (assels(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='assels'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (fadvus(i+1,j,k)-fadvus(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dfadvusdi'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           (fadvvs(i,j+1,k)-fadvvs(i,j,k))/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dfadvvsdj'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           ( &
            (  flus(i+1,j,k)-  flus(i,j,k)) &
           -(fadvus(i+1,j,k)-fadvus(i,j,k)) &
           )/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dfdifusdi'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax   ; do i=1,imax  
          anyvar3d(i,j,k)=mask_t(i,j,k)*( &
           ( &
            (  flvs(i,j+1,k)-  flvs(i,j,k)) &
           -(fadvvs(i,j+1,k)-fadvvs(i,j,k)) &
           )/dz_t(i,j,k,2)/dxdy_t(i,j) &
                                        )-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dfdifvsdj'       ; texte80(2)='1e-3'       ! variable ; units
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')
#endif

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de la densite
!******************************************************************************
      if(loop2_.eq.3)then !------->

      if(loop_netcdf==1) then !=======>
        if(eos_author==0) then !eeee> !18-11-14
         call equation_of_state_linear(1) !20-01-16 
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=rhp_t(i,j,k)
         enddo       ; enddo       ; enddo
        else                   !eeee>
         if(eos_tkezref>=0.) then !-------->
          call equation_of_state_potzref_jmfwg('grp',1) ! returns potential density in anyv3d(:,:,:,1) !20-10-14
         else                     !-------->
           call equation_of_state_potloc_jmfwg('grp',1) !23-02-22
         endif                    !-------->
        endif                  !eeee>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyvar3d(i,j,k)=mask_t(i,j,k)*(anyv3d(i,j,k,1)+rho-1000.) &
                     -(1.-mask_t(i,j,k))*9999.
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='rhp_1000'       ; texte80(2)='kg/m3' ! variable ; units !24-03-21
      if(eos_tkezref>=0.) then !-------->
       write(texte30,'(i0)')nint(max(eos_tkezref,zero))                !20-10-14
       texte80(3)='sea_water_potential_density_anomaly_reference_level_' &
                //trim(texte30)//'_meters'
      else                     !-------->
       texte80(3)='sea_water_potential_density_anomaly_local_reference' !04-03-22
      endif                    !-------->
      texte80(4)=texte80(3)
      call netcdf_main('_t')

#ifdef bidon
      if(loop_netcdf==1) then !=======>
        if(eos_author==0) then !eeee> !18-11-14
         call equation_of_state_linear(1) !20-01-16 
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=rhp_t(i,j,k)
         enddo       ; enddo       ; enddo
        else                   !eeee>
         if(eos_tkezref>=0.) then !-------->
          call equation_of_state_potzref_jmfwg('grp',1) ! returns potential density in anyv3d(:,:,:,1) !20-10-14
         else                     !-------->
           call equation_of_state_potloc_jmfwg('grp',1) !23-02-22
         endif                    !-------->
        endif                  !eeee>
      do k=2,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1
! JOHN R. TAYLOR AND RAFFAELE FERRARI, 2010 !28-11-16
! Buoyancy and Wind-Driven Convection at Mixed Layer Density Fronts
! JOURNAL OF PHYSICAL OCEANOGRAPHY
! DOI: 10.1175/2010JPO4365.1
! Bulk Ridchardson Number: f**2*(rho0/grav)*(drhp/dz) / ( (drhp/dx)**2 + (drhp/dy)**2 )
! Or RIB=N**2/(dVG/dz)**2=(db/dz) / ( (db/dx)/f)**2=f**2 * db/dz / (db/dx)**2 
! b=-grav*rhp/rho0

      x1=-( anyv3d(i,j,k+1,1)-anyv3d(i,j,k-1,1)) &
         /(depth_t(i,j,k+1) -depth_t(i,j,k-1)) 

      anyvar3d(i,j,k)=mask_t(i,j,k)*(coriolis_t(i,j)**2)*(rho/grav)*x1        & !f**2*(rho/g)*(drho/dz)
       /( ((   anyv3d(i+1,j  ,k,1)-anyv3d(i-1,j  ,k,1)                        & 
            +(depth_t(i+1,j  ,k) -depth_t(i-1,j  ,k))*x1 )/(2*dx_t(i,j)))**2  &
         +((   anyv3d(i  ,j+1,k,1)-anyv3d(i  ,j-1,k,1)                        &
            +(depth_t(i  ,j+1,k) -depth_t(i  ,j-1,k))*x1 )/(2*dy_t(i,j)))**2  &
         +small2) &

                     -(1.-mask_t(i,j,k))*9999.
          
      enddo       ; enddo       ; enddo
      anyvar3d(:,:,kmax)=anyvar3d(:,:,max(kmax-1,1))
      anyvar3d(:,:,1   )=anyvar3d(:,:,2)
      anyvar3d(0     ,:,:)=anyvar3d(2     ,:,:)
      anyvar3d(1     ,:,:)=anyvar3d(2     ,:,:)
      anyvar3d(imax+1,:,:)=anyvar3d(imax-1,:,:)
      anyvar3d(imax  ,:,:)=anyvar3d(imax-1,:,:)
      anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax-1,:)
      anyvar3d(:,jmax  ,:)=anyvar3d(:,jmax-1,:)
      anyvar3d(:,0     ,:)=anyvar3d(:,2     ,:)
      anyvar3d(:,1     ,:)=anyvar3d(:,2     ,:)
      endif                  !=======>
      texte80(1)='Ridchardson'     ; texte80(2)='?'    
      texte80(3)='Bulk_Ridchardson_Number'
      texte80(4)=texte80(3)
      call netcdf_main('_t')
#endif

      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression des traceurs passifs bios
!******************************************************************************
#ifdef bidon 
! 20-12-14
      if(imodelbio/=1) then !bbbbb>
      if(loop2_.eq.4)then !------->
       if (vbmax.gt.0) then !<><><><><>->

      do vb=1,vbmax !- vb loop ->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=0,jmax+1 !14-07-14
        do i=0,imax+1
         if (mask_t(i,j,kmax  )==1) then
          anyvar3d(i,j,k)=bio_t(i,j,max0(k,kmin_w(i,j)),vb)   !14-12-11
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo
        enddo
        enddo
      endif                  !=======>

      write(texte30,'(a7,i0)')'tracer_',vb
      texte80(1)=trim(texte30) ; texte80(2)='unity not defined'

      call netcdf_main('_t')

      enddo         !- vb loop ->

      endif                    !<><><><><>->
      endif           !------->
      endif                 !bbbbb>
#endif



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

#ifdef bidonref
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax   !30-07-14
         do i=1,imax+1
          anyvar3d(i,j,k)=mask_u(i,j,k)*salref_u(i,j,k) &
                     -(1.-mask_u(i,j,k))*9999.
         enddo
         enddo
         enddo
      endif                  !=======>
      texte80(1)='salref_u'     ; texte80(2)='degc'                       ! variable ; units
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax+1 !30-07-14
         do i=1,imax  
          anyvar3d(i,j,k)=mask_v(i,j,k)*salref_v(i,j,k) &
                     -(1.-mask_v(i,j,k))*9999.
         enddo
         enddo
         enddo
      endif                  !=======>
      texte80(1)='salref_v'     ; texte80(2)='degc'                       ! variable ; units
      call netcdf_main('_v')
#endif

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
       if(obcstatus(ieq1)==1) then !--------->

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
       if(obcstatus(ieqimax)==1) then !--------->

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
       if(obcstatus(jeq1)==1) then !--------->
        i1=1 ; i2=imax
        if(obcstatus(ieq1)==1)i1=bio_relax_size+1
        if(obcstatus(ieqimax)==1)i2=imax-bio_relax_size

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
       if(obcstatus(jeqjmax)==1) then !--------->
        i1=1 ; i2=imax
        if(obcstatus(ieq1)==1)i1=bio_relax_size+1
        if(obcstatus(ieqimax)==1)i2=imax-bio_relax_size

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
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=(tem_t(i,j,k,2)                              &
           -(timeweightobc(trc_id) *temobc_t(i,j,k,2)                  &
        +(1.-timeweightobc(trc_id))*temobc_t(i,j,k,0)))*mask_t(i,j,k)  &
                                             -9999.*(1.-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='tem_temobc'  ; texte80(2)='degc'              ! variable ; units
      call netcdf_main('_t')

#ifdef bidonref
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=(temref_t(i,j,k)                             &
           -(timeweightobc(trc_id) *temobc_t(i,j,k,2)                  &
        +(1.-timeweightobc(trc_id))*temobc_t(i,j,k,0)))*mask_t(i,j,k)  &
                                             -9999.*(1.-mask_t(i,j,k))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='temref_temobc'  ; texte80(2)='degc'              ! variable ; units
      call netcdf_main('_t')
#endif

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de salinite - forcage salinite
!******************************************************************************
      if(loop2_.eq.10)then !------->

       if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=(sal_t(i,j,k,2)                              &
           -(timeweightobc(trc_id) *salobc_t(i,j,k,2)                  &
        +(1.-timeweightobc(trc_id))*salobc_t(i,j,k,0)))*mask_t(i,j,k)  &
                                             -9999.*(1.-mask_t(i,j,k))
        enddo ; enddo ; enddo
       endif                  !=======>
       texte80(1)='sal_salobc'  ; texte80(2)='psu'               ! variable ; units
       call netcdf_main('_t')

#ifdef bidonref
       if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=(salref_t(i,j,k)                             &
           -(timeweightobc(trc_id) *salobc_t(i,j,k,2)                  &
        +(1.-timeweightobc(trc_id))*salobc_t(i,j,k,0)))*mask_t(i,j,k)  &
                                             -9999.*(1.-mask_t(i,j,k))
        enddo ; enddo ; enddo
       endif                  !=======>
       texte80(1)='salref_salobc'  ; texte80(2)='psu'               ! variable ; units
       call netcdf_main('_t')
#endif

      endif           !------->

      if(loop2_.eq.11)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyvar3d(i,j,k)=dz_t(i,j,k,1)*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval !29-01-21
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='dz'  ; texte80(2)='m' ! variable ; units
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

!     call modeanalysis_modeprojection

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

      do ktide=1,kmaxtide

! Produire au temps elapsedtime_now le champ mode par mode A partir des composantes
! harmoniques aux frequences de la maree

      if(loop_netcdf==1) then !=======>
      if(allocated(uv_wmode_t)) then !)))))>
        time1=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide))+v0tide(ktide)+utide(ktide,1)       
        const3=cos(time1)*ftide(ktide,1) ; const4=sin(time1)*ftide(ktide,1)  
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
!       if (mask_t(i,j,kmax  ).ne.0) then
!          anyvar3d(i,j,k)=uv_wmode_t(i,j,k,loopm_)*ucoefmode_t(i,j,loopm_)
!       else
!          anyvar3d(i,j,k)=-9999.
!       endif
           anyvar3d(i,j,k)=                              &
           filval*(1-mask_t(i,j,kmax))+mask_t(i,j,kmax)* &
            uv_wmode_t(i,j,k,loopm_)*( & !pmx>        ! mode(z)* 
           u3dmode_cos_t(i,j,loopm_,ktide)*const3   & ![AmpCos*cosinus(wt)
          +u3dmode_sin_t(i,j,loopm_,ktide)*const4   & !+AmpSin*  sinus(wt)]
                                     )   !pmx>        !
                                                     
        enddo ; enddo ; enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'umode',loopm_,'tide',ktide
      texte80(2)='m.s-1'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(uv_wmode_t)) then !)))))>
        time1=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide))+v0tide(ktide)+utide(ktide,1)       
        const3=cos(time1)*ftide(ktide,1) ; const4=sin(time1)*ftide(ktide,1)  
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
!       if (mask_t(i,j,kmax  ).ne.0) then
!          anyvar3d(i,j,k)=uv_wmode_t(i,j,k,loopm_)*vcoefmode_t(i,j,loopm_)
!       else
!          anyvar3d(i,j,k)=-9999.
!       endif
           anyvar3d(i,j,k)=                              &
           filval*(1-mask_t(i,j,kmax))+mask_t(i,j,kmax)* &
            uv_wmode_t(i,j,k,loopm_)*( & !pmx>        ! mode(z)* 
           v3dmode_cos_t(i,j,loopm_,ktide)*const3   & ![AmpCos*cosinus(wt)
          +v3dmode_sin_t(i,j,loopm_,ktide)*const4   & !+AmpSin*  sinus(wt)]
                                     )   !pmx>        !
        enddo ; enddo ; enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'vmode',loopm_,'tide',ktide
      texte80(2)='m.s-1'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !=======>
      if(allocated(pcoefmode_t)) then !)))))> !14-04-16
        time1=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide))+v0tide(ktide)+utide(ktide,1)       
        const3=cos(time1)*ftide(ktide,1) ; const4=sin(time1)*ftide(ktide,1)  
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
           anyvar3d(i,j,k)=                              &
           filval*(1-mask_t(i,j,kmax))+mask_t(i,j,kmax)* &
            uv_wmode_t(i,j,k,loopm_)*( & !pmx>        ! mode(z)* 
           p3dmode_cos_t(i,j,loopm_,ktide)*const3   & ![AmpCos*cosinus(wt)
          +p3dmode_sin_t(i,j,loopm_,ktide)*const4   & !+AmpSin*  sinus(wt)]
                                     )   !pmx>        !
       enddo ; enddo ; enddo
      endif                          !)))))>
      endif                  !=======>
      write(texte80(1),'(a,i0,a,i0)')'pmode',loopm_,'tide',ktide
      texte80(2)='kg.m**-1.s**-2'
      call netcdf_main('_t')


      enddo ! fin de boucle sur ktide
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
! Potential density of the forcing field (time=2)
!******************************************************************************
      if(loop2_==16)then !------->

      if(loop_netcdf==1) then !=======>
        if(eos_tkezref>=0.) then !-------->
         call equation_of_state_potzref_jmfwg('obc',1) ! returns potential density in anyv3d(:,:,:,1) !20-10-14
        else                     !-------->
          call equation_of_state_potloc_jmfwg('obc',1) !23-02-22
        endif                    !-------->
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
         anyvar3d(i,j,k)= &
          (anyv3d(i,j,max(k,kmin_w(i,j)),1)+rho-1000.)*mask_t(i,j,kmax) &
                                            -9999.*(1.-mask_t(i,j,kmax))
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='rhpobc' ; texte80(2)='kg/m3'                        ! variable ; units
      write(texte30,'(i0)')nint(max(eos_tkezref,zero))                !20-10-14
      texte80(3)='forcing_sea_water_potential_density_reference_level_' &
               //trim(texte30)//'_meters'
      texte80(4)=texte80(3)
      call netcdf_main('_t')


      endif              !------->

      if(loop2_==17)then !17-17-17> 
      stop 'Cet espace est disponible'
      endif              !17-17-17>


! Vertical velocity !12-06-18
! https://docs.google.com/document/d/1P9YmC88Un2LPVI36eUNS0vY9uRph3bJc9JQ6-gy-8xo/edit
      if(loop2_==18)then !------->
       if(loop_netcdf==1) then !=======>
         anyvar3d=filval !15-12-21
!        call omega_vertical_velocity(1) !12-06-18
         call omega_vertical_velocity(2,2,imax-1,2,jmax-1) !10-12-21!15-12-21
!        anyvar3d(0     ,:     ,:)=-9999.
!        anyvar3d(imax+1,:     ,:)=-9999.
!        anyvar3d(:     ,0     ,:)=-9999.
!        anyvar3d(:     ,jmax+1,:)=-9999.
       endif                   !=======>
       texte80(1)='W' ; texte80(2)='m/s'  ! variable ; units
       texte80(3)='true_vertical_velocity'
       texte80(4)=texte80(3)
       call netcdf_main('_t')
      endif              !------->


!******************************************************************************
! impression drifters 3D
!******************************************************************************
      if(loop2_==19)then!------->

      if(loop_netcdf==1) then !=======>

       anyv3d(:,:,:,1)=0.
       do kbu=1,kbumax
        i=nint(drifter_l(kbu,1)-par%timax(1))
        j=nint(drifter_l(kbu,2)-par%tjmax(1))
        k=nint(drifter_l(kbu,3))
        if(i>=2.and.i<=imax-1.and.j>=2.and.j<=jmax-1) then !>>>
         anyv3d(i,j,k,1)=anyv3d(i,j,k,1)+1
        endif                                              !>>>
       enddo
       call obc_mpi_anyv3d(0,1,'z0')

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyvar3d(i,j,k)=mask_t(i,j,kmax)*anyv3d(i,j,k,1) &
                    +(1-mask_t(i,j,kmax))*filval
       enddo       ; enddo       ; enddo

      endif                  !=======>
      texte80(1)='drifters3D' ; texte80(2)='nb/maille'
      texte80(3:4)='nbre_drifter_par_maille'
      call netcdf_main('_t')

      endif         !------->


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

       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1        
        anyvar3d(i,j,k)=(       vel_u(i,j,k,1)                & !16-02-23
                         +velstokes_u(i,j,k,1))*mask_u(i,j,k) & 
                                     +filval*(1-mask_u(i,j,k))
       enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='vel_u'  
      texte80(2)='m s-1'    
      texte80(3)='sea_water_x_velocity_at_u_location'
      texte80(4)='sea_water_x_velocity_at_u_location'
      call netcdf_main('_u')

      if(allocated(flx3d_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1  
          anyvar3d(i,j,k)=flx3d_timeaveraged(i,j,k)*mask_u(i,j,kmax)   &
                                                +(1-mask_u(i,j,kmax))*filval
        enddo       ; enddo       ; enddo
      endif                  !=======>
        texte80(1)='flx3d_timeaveraged'  
        texte80(2)='m**2 s-1'    
        texte80(3:4)='flx3d_timeaveraged'
        call netcdf_main('_u')
      endif                                  !pmx>

      if(allocated(u_euler_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1  
          anyvar3d(i,j,k)=u_euler_timeaveraged(i,j,k)*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
        enddo       ; enddo       ; enddo
      endif                  !=======>
        texte80(1)='u_euler_timeaveraged'
        texte80(2)='m*s-1'    
        texte80(3:4)='u_euler_timeaveraged'
        call netcdf_main('_u')
      endif                                  !pmx>

#ifdef bilan_s_3d
      if(allocated(fadvus)) then !ooo> !20-06-22
        if(loop_netcdf==1) then !=======>
          do k=1,kmax ; do j=1,jmax ; do i=1,imax+1  
! attention au signe -
            anyvar3d(i,j,k)=-fadvus(i,j,k)*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
          enddo       ; enddo       ; enddo
        endif                   !=======>
        texte80(1)='flux_oi_S_adv'
        texte80(2)='psu*m**3'    
        texte80(3:4)=texte80(1)
        call netcdf_main('_u')

        if(loop_netcdf==1) then !=======>
          do k=1,kmax ; do j=1,jmax ; do i=1,imax+1  
! attention au signe -
            anyvar3d(i,j,k)=-(flus(i,j,k)-fadvus(i,j,k))*mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
          enddo       ; enddo       ; enddo
        endif                   !=======>
        texte80(1)='flux_oi_S_dif'
        texte80(2)='psu*m**3'    
        texte80(3:4)=texte80(1)
        call netcdf_main('_u')
      endif                    !ooo>
#endif

#ifdef bidon
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1                     
          anyvar3d(i,j,k)=veldydz_u(i,j,k,1)
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='veldydz_u' !12-01-17 
      texte80(2)='m**3/s'    
      texte80(3:4)='Oi_waterflux'
      call netcdf_main('_u')
#endif
!#ifdef bidon
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax+1                     
          anyvar3d(i,j,k)=veldydz_u(i,j,k,1) &
                              /dx_u(i,j)     &
                              /dy_u(i,j)     &
                              /dz_u(i,j,k,1) &
                               *dti_lp       &
                            *mask_u(i,j,kmax)+(1-mask_u(i,j,kmax))*filval
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='current_number_u' !13-01-23
      texte80(2)='?'    
      texte80(3:4)=texte80(1)
      call netcdf_main('_u')
!#endif

      if(loop_netcdf==1) then !=======>
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        anyvar3d(i,j,k)=(       vel_v(i,j,k,1)                & !16-04-23
                         +velstokes_v(i,j,k,1))*mask_v(i,j,k) & 
                                     +filval*(1-mask_v(i,j,k))
       enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='vel_v'  
      texte80(2)='m s-1' 
      texte80(3)='sea_water_y_velocity_at_v_location'
      texte80(4)='sea_water_y_velocity_at_v_location'
      call netcdf_main('_v')

      if(allocated(fly3d_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
          anyvar3d(i,j,k)=fly3d_timeaveraged(i,j,k)*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
        enddo       ; enddo         ; enddo
      endif                  !=======>
        texte80(1)='fly3d_timeaveraged'
        texte80(2)='m**2 s-1' 
        texte80(3:4)='fly3d_timeaveraged'
        call netcdf_main('_v')
      endif                                  !pmx>

      if(allocated(v_euler_timeaveraged)) then !pmx>
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
          anyvar3d(i,j,k)=v_euler_timeaveraged(i,j,k)*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
        enddo       ; enddo         ; enddo
      endif                  !=======>
        texte80(1)='v_euler_timeaveraged'
        texte80(2)='m*s-1' 
        texte80(3:4)='v_euler_timeaveraged'
        call netcdf_main('_v')
      endif                                  !pmx>

#ifdef bilan_s_3d
      if(allocated(fadvvs)) then !ooo> !20-06-22
        if(loop_netcdf==1) then !=======>
          do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax  
! attention au signe -
            anyvar3d(i,j,k)=-fadvvs(i,j,k)*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
          enddo       ; enddo       ; enddo
        endif                   !=======>
        texte80(1)='flux_oj_S_adv'
        texte80(2)='psu*m**3'    
        texte80(3:4)=texte80(1)
        call netcdf_main('_v')

        if(loop_netcdf==1) then !=======>
          do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax  
! attention au signe -
            anyvar3d(i,j,k)=-(flvs(i,j,k)-fadvvs(i,j,k))*mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
          enddo       ; enddo       ; enddo
        endif                   !=======>
        texte80(1)='flux_oj_S_dif'
        texte80(2)='psu*m**3'    
        texte80(3:4)=texte80(1)
        call netcdf_main('_v')
      endif                    !ooo>
#endif

#ifdef bidon
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1                          !08-08-10
        do i=1,imax
          anyvar3d(i,j,k)=veldxdz_v(i,j,k,1)
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='veldxdz_v' !12-01-17 
      texte80(2)='m**3/s'    
      texte80(3:4)='Oj_waterflux'
      call netcdf_main('_v')
#endif
!#ifdef bidon
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax+1
        do i=1,imax                     
          anyvar3d(i,j,k)=veldxdz_v(i,j,k,1) &
                              /dx_v(i,j)     &
                              /dy_v(i,j)     &
                              /dz_v(i,j,k,1) &
                               *dti_lp       &
                            *mask_v(i,j,kmax)+(1-mask_v(i,j,kmax))*filval
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='current_number_v' !13-01-23
      texte80(2)='?'    
      texte80(3:4)=texte80(1)
      call netcdf_main('_v')
!#endif

      endif           !------->

!******************************************************************************


!******************************************************************************
! impression du courant geostrophique
!******************************************************************************
      if(loop2_==2)then !------->

      if(loop_netcdf==1) then !=======>
      call graph_out_geostrophic_current('3d')
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         if (mask_u(i,j,k)/=0) then
          anyvar3d(i,j,k)=anyv3d(i,j,k,3)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='v_geos'  ; texte80(2)='m/s'               ! variable ; units
      texte80(3)='geostrophic_v_current_at_u_location' !29-01-15
      texte80(4)='geostrophic_v_current_at_u_location'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
      call graph_out_geostrophic_current('3d')
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         if (mask_v(i,j,k)/=0) then
          anyvar3d(i,j,k)=anyv3d(i,j,k,4)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='u_geos'  ; texte80(2)='m/s'               ! variable ; units
      texte80(3)='geostrophic_u_current_at_v_location' !29-01-15
      texte80(4)='geostrophic_u_current_at_v_location'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !=======>
      call graph_out_geostrophic_current('3d-2d')
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         if (mask_u(i,j,k)/=0) then
          anyvar3d(i,j,k)=anyv3d(i,j,k,3)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='v_geos0'  ; texte80(2)='m/s'               ! variable ; units
      texte80(3)='geostrophic_v_current_at_u_location_ssh_removed' !29-01-15
      texte80(4)=texte80(3)
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
      call graph_out_geostrophic_current('3d-2d')
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         if (mask_v(i,j,k)/=0) then
          anyvar3d(i,j,k)=anyv3d(i,j,k,4)
         else
          anyvar3d(i,j,k)=-9999.
         endif
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='u_geos0'  ; texte80(2)='m/s'               ! variable ; units
      texte80(3)='geostrophic_u_current_at_v_location_ssh_removed' !29-01-15
      texte80(4)=texte80(3)
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
              anyvar3d(i,j,k)=vel_u(i,j,k,2)     &
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
              anyvar3d(i,j,k)=vel_v(i,j,k,2)     &
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
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        anyvar3d(i,j,k)=velstokes_u(i,j,k,1)*mask_u(i,j,k) & 
                                  +filval*(1-mask_u(i,j,k)) 
       enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='velstokes_u'  ; texte80(2)='m/s'                    ! variable ; units
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        anyvar3d(i,j,k)=velstokes_v(i,j,k,1)*mask_v(i,j,k) & 
                                  +filval*(1-mask_v(i,j,k))
       enddo ; enddo ; enddo
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
              anyvar3d(i,j,k)=vel_u(i,j,k,2)         !25-02-10
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
              anyvar3d(i,j,k)=vel_v(i,j,k,2)         !25-02-10
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
        anyvar3d=-9999.
        do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
         anyvar3d(i,j,k)=presgrad_u(i,j,k,1)/dx_u(i,j)*mask_u(i,j,kmax) &
                        /max(coriolis_t(i,j),small1)                    &
                                           -9999.*(1-mask_u(i,j,kmax))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='presgrad_u'  ; texte80(2)='m/s'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        anyvar3d=-9999.
        do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
         anyvar3d(i,j,k)=presgrad_v(i,j,k,1)/dy_v(i,j)*mask_v(i,j,kmax) &
                        /max(coriolis_t(i,j),small1)                    &
                                           -9999.*(1-mask_v(i,j,kmax))
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='presgrad_v'  ; texte80(2)='m/s'
      call netcdf_main('_v')

      endif           !------->

!******************************************************************************
! dT/dxs , dT/dys , dS/dxs , dS/dys !02-10-21
!******************************************************************************
      if(loop2_==8)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         anyvar3d(i,j,k)=mask_u(i,j,kmax)*invdx_u(i,j)*(tem_t(i,j,k,1)-tem_t(i-1,j,k,1)) &
                    +(1.-mask_u(i,j,kmax))*filval
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='TgradX'  ; texte80(2)='d/m'
      texte80(3:4)='isosigma_X_gradient_of_temperature'

      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         anyvar3d(i,j,k)=mask_v(i,j,kmax)*invdy_v(i,j)*(tem_t(i,j,k,1)-tem_t(i,j-1,k,1)) &
                    +(1.-mask_v(i,j,kmax))*filval
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='TgradY'  ; texte80(2)='d/m'
      texte80(3:4)='isosigma_Y_gradient_of_temperature'
      call netcdf_main('_v')


      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         anyvar3d(i,j,k)=mask_u(i,j,kmax)*invdx_u(i,j)*(sal_t(i,j,k,1)-sal_t(i-1,j,k,1)) & 
                    +(1.-mask_u(i,j,kmax))*filval
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='SgradX'  ; texte80(2)='./m'
      texte80(3:4)='isosigma_X_gradient_of_salinity'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         anyvar3d(i,j,k)=mask_v(i,j,kmax)*invdy_v(i,j)*(sal_t(i,j,k,1)-sal_t(i,j-1,k,1)) & 
                    +(1.-mask_v(i,j,kmax))*filval
        enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='SgradY'  ; texte80(2)='./m'
      texte80(3:4)='isosigma_Y_gradient_of_salinity'
      call netcdf_main('_v')

      endif           !------->
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
        do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
         anyvar3d(i,j,k)=kh_w(i,j,k)*mask_t(i,j,k)-9999.*(1-mask_t(i,j,k)) !14-02-17
        enddo ; enddo ; enddo
        anyvar3d(0:imax+1:imax+1,:,:)=filval ! filval en i=0 et i=imax+1  !19-05-18
        anyvar3d(:,0:jmax+1:jmax+1,:)=filval ! filval en j=0 et j=jmax+1 
      endif                    !=======>
      texte80(1)='kh'  ; texte80(2)='m2/s'                          ! variable ; units
      call netcdf_main('_w')

#ifdef bidon
! Diffusite verticale numerique du schema d'advection upwind:
       if(loop_netcdf==1) then !=======>
        do k=2,kmax ; do j=1,jmax ; do i=1,imax
         anyvar3d(i,j,k)=0.5*abs(omega_w(i,j,k,1))*(depth_t(i,j,k)-depth_t(i,j,k-1))*mask_t(i,j,k)-9999.*(1-mask_t(i,j,k))
        enddo ; enddo ; enddo
        anyvar3d(:,:,1)=0.
        anyvar3d(:,:,kmax+1)=0.
        anyvar3d(0:imax+1:imax+1,:,:)=filval ! filval en i=0 et i=imax+1  !19-05-18
        anyvar3d(:,0:jmax+1:jmax+1,:)=filval ! filval en j=0 et j=jmax+1 
      endif                    !=======>
      texte80(1)='kh_upwind_advection'  ; texte80(2)='m2/s'                          ! variable ; units
      call netcdf_main('_w')
#endif

!      if(loop_netcdf==1) then !=======>
!       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
!        anyvar3d(i,j,k)=km_w(i,j,k)*mask_t(i,j,k)-9999.*(1-mask_t(i,j,k)) !14-02-17
!       enddo ; enddo ; enddo
!       anyvar2d(0:imax+1:imax+1,:)=filval ! filval en i=0 et i=imax+1 
!       anyvar2d(:,0:jmax+1:jmax+1)=filval ! filval en j=0 et j=jmax+1 
!     endif                  !=======>
!     texte80(1)='km'  ; texte80(2)='m2/s'                          ! variable ; units
!     call netcdf_main('_w')

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de l'energie cinetique turbulente
!******************************************************************************
      if(loop2_.eq.2)then !------->

      if(allocated(tken_w)) then
      if(loop_netcdf==1) then !=======>
       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
        anyvar3d(i,j,k)=tken_w(i,j,k)*mask_t(i,j,kmax)-9999.*(1-mask_t(i,j,kmax))
       enddo ; enddo ; enddo
       anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
       anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
       anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
       anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='tken'  ; texte80(2)='(m/s)2'                     ! variable ; units
      call netcdf_main('_w')
      endif

      endif           !------->
!******************************************************************************


!******************************************************************************
! impression de la longueur de melange
!******************************************************************************
      if(loop2_.eq.3)then !------->
      if(    iturbulence==0   &
         .or.iturbulence==1   &
         .or.iturbulence==4   & !09-02-19
         .or.iturbulence==5   & !09-02-19
                           ) then !pmxpmx> !29-03-17

      if(loop_netcdf==1) then !=======>

      if(iturbulence==0)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkll_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==1)anyvar3d(1:imax,1:jmax,1:kmax+1)=epsn_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==4)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkll_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==5)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkll_w(1:imax,1:jmax,1:kmax+1)

      do k=1,kmax+1
      do j=1,jmax
      do i=1,imax
       if(mask_t(i,j,kmax  )==0)anyvar3d(i,j,k)=-9999.
      enddo
      enddo
      enddo
      anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:) !17-12-16
      anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
      anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
      anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)

      endif                  !=======>
      texte80(1)='tkll'
      texte80(2)='m'
      if(iturbulence==1)texte80(1)='eps'
      if(iturbulence==1)texte80(2)='m2/s3'
      call netcdf_main('_w')

      endif                                     !pmxpmx> !29-03-17
      endif           !------->

!******************************************************************************


!******************************************************************************
! impression de la longueur de dissipation
!******************************************************************************
      if(loop2_.eq.4)then !------->
      if    (iturbulence==0    &
         .or.iturbulence==1    &
         .or.iturbulence==3    &
         .or.iturbulence==4    & !09-02-19
         .or.iturbulence==5    & !09-02-19
                           ) then !pmxpmx> !29-03-17

      if(loop_netcdf==1) then !=======>

      if(iturbulence==0)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkle_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==1)anyvar3d(1:imax,1:jmax,1:kmax+1)=  km_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==3)anyvar3d(1:imax,1:jmax,1:kmax+1)=  km_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==4)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkle_w(1:imax,1:jmax,1:kmax+1)
      if(iturbulence==5)anyvar3d(1:imax,1:jmax,1:kmax+1)=tkle_w(1:imax,1:jmax,1:kmax+1)

      do k=1,kmax+1
      do j=1,jmax
      do i=1,imax
       if(mask_t(i,j,kmax  )==0)anyvar3d(i,j,k)=-9999.
      enddo
      enddo
      enddo
      anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:) !17-12-16
      anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
      anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
      anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)

      endif                  !=======>
      texte80(1)='tkle'
      if(iturbulence==1)texte80(1)='km'
      if(iturbulence==3)texte80(1)='km'
      texte80(2)='m'
      if(iturbulence==1)texte80(2)='m2/s'
      if(iturbulence==3)texte80(2)='m2/s'
      call netcdf_main('_w')

      endif                                     !pmxpmx>
      endif                     !------->
!******************************************************************************


!******************************************************************************
! impression de omega
!******************************************************************************

      if(loop2_==5)then !------->

      if(loop_netcdf==1) then !=======>
        do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax

              anyvar3d(i,j,k)=omega_w(i,j,k,1)*mask_t(i,j,kmax)  &
                                     -9999.*(1-mask_t(i,j,kmax))
        enddo ; enddo ; enddo
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='omega'  ; texte80(2)='m s-1'           
      texte80(3:4)='relative_vertical_sea_water_velocity_at_w_location'
      call netcdf_main('_w')

      if(loop_netcdf==1) then !=======> !30-12-16
        do k=2,kmax   ; do j=1,jmax ; do i=1,imax

              anyvar3d(i,j,k)=omega_w(i,j,k,1)*mask_t(i,j,kmax)    &
                             *dti_lp                               &
                            /min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))    &
                                     -9999.*(1-mask_t(i,j,kmax))
        enddo ; enddo ; enddo
        anyvar3d(:,:,1)=0.
        anyvar3d(:,:,kmax+1)=0.
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='omega_dtodz'  ; texte80(2)='without_units'
      texte80(3:4)='vertical_current_number'
      call netcdf_main('_w')

      if(allocated(q_t)) then !>>>
      if(loop_netcdf==1) then !=======>
        do k=1,kmax   ; do j=1,jmax ; do i=1,imax

              anyvar3d(i,j,k)=q_t(i,j,k,1) !*mask_t(i,j,kmax)  &
                               !  -9999.*(1-mask_t(i,j,kmax))
        enddo ; enddo ; enddo
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='q_t'  ; texte80(2)='?'           
      texte80(3:4)='q_t'
      call netcdf_main('_t')

#ifdef bidon
      if(loop_netcdf==1) then !=======>
      k0=1 ; if(flag_nh3d==1)k0=0
        do k=1,kmax   ; do j=1,jmax ; do i=1,imax
              anyvar3d(i,j,k)=(q_t(i,j,k,1)         &
                         +grav*ssh_w(i,j,1)*k0      &
                     +grav*ssh_int_w(i,j,1)*(1-k0) ) !*mask_t(i,j,kmax)  &
        enddo ; enddo ; enddo
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='q+grav*ssh'  ; texte80(2)='?'           
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')
#endif
      endif                  !>>>

#ifdef bidon
      if(loop_netcdf==1) then !=======>
        do k=1,kmax   ; do j=1,jmax ; do i=1,imax

              anyvar3d(i,j,k)=qavr_t(i,j,k) !*mask_t(i,j,kmax)  &
                               !  -9999.*(1-mask_t(i,j,kmax))
        enddo ; enddo ; enddo
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='qavr_t'  ; texte80(2)='?'           
      texte80(3:4)='qavr_t'
      call netcdf_main('_t')
#endif

#ifdef bidon
!POUR UTILISER CES LIGNES IL FAUT CONNAITRE LA SOLUTION ANALYTIQUE....
! Il faut notamment connaitre la longueur d'ondes:
      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax

          x4=grav*ssh_w(i,j,1)*(cosh(kvector_*(h_w(i,j)+depth_t(i,j,k))) &
                               /cosh(kvector_*(h_w(i,j)               )) &
                                -1.)

              anyvar3d(i,j,k)= x4*mask_t(i,j,kmax)  &
                        -9999.*(1-mask_t(i,j,kmax))
        enddo ; enddo ; enddo
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='q_theory'  ; texte80(2)='?'           
      texte80(3:4)=texte80(1)
      call netcdf_main('_t')
#endif

      endif           !------->

!******************************************************************************
! (Brunt-Vaisala frequency)**2  !19-03-21
!******************************************************************************
      if(loop2_==6)then !------->
      if(loop_netcdf==1) then !=======>
        if(eos_author==0) then !eeee> !18-11-14
         call equation_of_state_linear(1) !20-01-16 
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=rhp_t(i,j,k)
         enddo       ; enddo       ; enddo
        else                   !eeee>
         if(eos_tkezref>=0.) then !-------->
          call equation_of_state_potzref_jmfwg('grp',1) ! returns potential density in anyv3d(:,:,:,1) !20-10-14
         else                     !-------->
           call equation_of_state_potloc_jmfwg('grp',1) !23-02-22
         endif                    !-------->
        endif                  !eeee>
        do k=2,kmax ; do j=1,jmax ; do i=1,imax
! (Brunt-Vaisala frequency)**2:                                !19-03-21
          anyvar3d(i,j,k)=mask_t(i,j,k)                      &
               *grav/rho*(anyv3d(i,j,k,1)-anyv3d(i,j,k-1,1)) &
                       /(depth_t(i,j,k) -depth_t(i,j,k-1))   & !19-03-21
                     -(1.-mask_t(i,j,k))*9999.
        enddo ; enddo ; enddo
        anyvar3d(:,:,1)     =anyvar3d(:,:,2)
        anyvar3d(:,:,kmax+1)=anyvar3d(:,:,kmax)
      endif                  !=======>
      texte80(1)='N2'       ; texte80(2)='s**-2'                ! variable ; units
      texte80(3:4)='(Brunt-Vaisala-frequency)**2'
      call netcdf_main('_w')
      endif             !------->
!******************************************************************************
! (Brunt-Vaisala frequency)**2  !19-03-21
!******************************************************************************

!#ifdef bidon
!******************************************************************************
! sigma_fric_w: fraction of the critical height for macroscale friction !27-05-22
!******************************************************************************
      if(allocated(sigma_fric_wu).and.allocated(sigma_fric_wv)) then !pmx>
      if(loop2_==7)then !------->
      if(loop_netcdf==1) then !=======>
        do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
              anyvar3d(i,j,k)=0.25*( sigma_fric_wu(i  ,j  ,k)  &
                                    +sigma_fric_wu(i+1,j  ,k)  &
                                    +sigma_fric_wv(i  ,j  ,k)  &
                                    +sigma_fric_wv(i  ,j+1,k))*mask_t(i,j,kmax)  &
                                                     -9999.*(1-mask_t(i,j,kmax))
        enddo ; enddo ; enddo
        anyvar3d(0     ,:,:)=anyvar3d(1   ,:,:)
        anyvar3d(imax+1,:,:)=anyvar3d(imax,:,:)
        anyvar3d(:,jmax+1,:)=anyvar3d(:,jmax,:)
        anyvar3d(:,0     ,:)=anyvar3d(:,1   ,:)
      endif                  !=======>
      texte80(1)='sigma_fric_w'  ; texte80(2)='none'           
      texte80(3:4)='sigma_fric_w'
      call netcdf_main('_w')
      endif             !------->
      endif                               !pmx>
!******************************************************************************
! sigma_fric_w: fraction of the critical height for macroscale friction !27-05-22
!******************************************************************************
!#endif

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
      restart_file_y_or_n=1
      if(restartfileperiod>0.)call dyn_restart('w')
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

#ifdef bidon
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

#endif
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

#ifdef bidon
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

#endif
      end subroutine graph_out_trueaxis_nfmpi
!........................................................................................

      subroutine graph_out_geostrophic_current(txt_)
      use module_principal ; use module_parallele
      implicit none
      integer :: flag_=-999
      character(len=*) :: txt_
#ifdef synopsis
       subroutinetitle='graph_out_geostrophic_current'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      if(txt_=='3d-2d')   flag_=0
      if(txt_=='3d')      flag_=1
      if(txt_=='standard')flag_=1 !14-01-16 ! pour cas baroclinic_jet

      if(flag_==-999) &
      stop 'graph_out_geostrophic_current txt_ not recognized'

      if(abs(lat_t(imax/2,jmax/2)*rad2deg)<5.) &
      stop 'Two small latitude for geostrophic current computing'

! Tres important si frontiere ouvertes: reset a zero:
      anyv3d(:,:,:,1:4)=0.
      do j1=2,jmax-1
      do i1=2,imax

       xy_u(i1,j1,1)=(1./rho)*( & !ppp>
        ssh_int_w(i1  ,j1,1)*(rho+rhpzavr_w(i1  ,j1)) &!grav/rho*d(ssh.rhp2d)/di
       -ssh_int_w(i1-1,j1,1)*(rho+rhpzavr_w(i1-1,j1)) &!
                              )   !ppp>

        do k=1,kmax
         anyv3d(i1,j1,k,3)=                                            &!29-01-15
                           (+presgrad_u(i1,j1,k,1)                     &
!               +flag_*grav*(ssh_int_w(i1,j1,1)-ssh_int_w(i1-1,j1,1))) &
                +flag_*grav*xy_u(i1,j1,1))                             &
                        /(0.5*(coriolis_t(i1,j1)+coriolis_t(i1-1,j1))) &
                        *mask_u(i1,j1,k)/dx_u(i1,j1)

        enddo
      enddo
      enddo
      do j1=2,jmax
      do i1=2,imax-1

       xy_v(i1,j1,1)=(1./rho)*( & !ppp>
        ssh_int_w(i1,j1  ,1)*(rho+rhpzavr_w(i1,j1  )) &!grav/rho*d(ssh.rhp2d)/dj
       -ssh_int_w(i1,j1-1,1)*(rho+rhpzavr_w(i1,j1-1)) &!
                                    )   !ppp>


        do k=1,kmax
         anyv3d(i1,j1,k,4)=                                            &!29-01-15
                           (-presgrad_v(i1,j1,k,1)                     &
!               -flag_*grav*(ssh_int_w(i1,j1,1)-ssh_int_w(i1,j1-1,1))) &
                -flag_*grav*xy_v(i1,j1,1))                             &
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
!     call obc_mpi_anyv3d(0,1,'u1') ! Echange 'u1' sur anyv3d(:,:,:,1)
!     call obc_mpi_anyv3d(0,2,'v1') ! Echange 'v1' sur anyv3d(:,:,:,2)
      call obc_mpi_anyv3d(0,3,'u1') ! Echange 'u1' sur anyv3d(:,:,:,1) !29-01-15
      call obc_mpi_anyv3d(0,4,'v1') ! Echange 'v1' sur anyv3d(:,:,:,2) !29-01-15

      if(txt_=='standard')then !sssssssss>  !14-01-16 ! pour cas baroclinic_jet
!_____________________________________________________________________________

! Step 3: assemblage des elements du courant geostrophiques
! Composante u:
      anyv3d(:,:,:,2)=anyv3d(:,:,:,4) !14-01-16
      anyv3d(:,:,:,1)=anyv3d(:,:,:,3) !14-01-16
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
      call obc_mpi_anyv3d(1,3,'u1') ! Echange 'u1' sur anyv3d(:,:,:,1)
      call obc_mpi_anyv3d(1,4,'v1') ! Echange 'v1' sur anyv3d(:,:,:,2)

      endif                    !sssssssss>  !14-01-16 ! pour cas baroclinic_jet

      end subroutine graph_out_geostrophic_current
