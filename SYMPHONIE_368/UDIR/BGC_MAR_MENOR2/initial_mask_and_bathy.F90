      subroutine initial_mask_and_bathy
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 - last update: 09-03-23
!______________________________________________________________________

      use module_principal ; use module_parallele ; use module_global ; use module_mangrove
      use module_ogcm  ; use module_grid ; use module_webcanals
      implicit none
      real rmax_
      double precision hmaxglob,hminglob
      integer :: flag_,len_,case_=4
      real h1_cl(825,450)

#ifdef synopsis
       subroutinetitle='initial_mask_and_bathy'
       subroutinedescription= &
       'Reads notebook_bathy, sea-land mask and bathymetry'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
!  _________                    .__                  .__              ! m[°v°]m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____  
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!...............................................................................
! Version  date     Description des modifications
!          03/09/01: construction des fleuves deplacee apres options
!                    "traitement" de la bathymetrie
!          18/09/01: HMAX ne concerne plus que les points de mer
!          05/10/01: le remplissage de mask aux points _x _y _r est fait
!                    par un appel au sousprogramme z_to_xyr.F
!                    De plus la chronologie des actions est modifiee afin
!                    que la bathy des canaux n'intervienne pas dans les
!                    procedures de lissage
!          30/10/01: ajout d'une securite sur le masque: un point de
!                    frontiere contigu e un point de terre est forcement
!                    en terre (sinon les obc ne marchent pas) et est donc
!                    automatiquement masque.
!          01/04/03: lecture de bathycote.dat en format libre pour avoir
!                    exactement la meme bathy dans la simu imbriquee
!                    "notebook_streamf". Sinon l'ecriture du fichier bathycote
!                    imbriquee en format 10.3 tronque la bathy et apres on a
!                    pas exactement la meme bathy ce qui est prejudiciable au
!                    bon fonctionnement de la grille hybride.
!          05/06/03: bienvenue filiere imbrication bathy "parent"
!          24/03/04: appel e HZ_TO_HXYR
!          15/12/04: limites extremes mask special schema d'advection
!          07/01/05: retour en arriere par rapport au point precedent
!          19/06/05: message e l'ecran:
!          18/04/06: fonctions compatibles avec double precision
!          19/01/07: compatible avec model_ 1D...
!          03/04/07: filiere model_in model_out definitivement supprimee
!          04/04/07: on ne lisse pas si RMAX sup ou egal e 1
!          07/06/07: texte90 remplace texte60 dans lecture du notebook_bathy
!          21/12/07: changement de l'ordre des options de lissage de la bathy
!                    On lisse d'abord et ensuite on applique le surplus de
!                    lissage necessaire pour respecter le critere concernant
!                    l'erreur de troncature.
!          03/01/08: correction format de lecture
!          16/10/08: ajout d'une routine de C.L. sur   mask_x &   mask_y
!          28/11/08: Ajout d'une securite: reperage des piscines
!          06/03/09: regroupement des conditions aux limites pour faciliter
!                    la parallelisation
!          09/03/09: Conditions aux limites pour h
!                    revoir l'ordre des choses: tridimensionnaliser le masque
!                    apres lecture fichier masque+bathy
!          06/04/09: Calculer HMAXGLOB & HMINGLOB
!          04-06-09  Calculer l'aire de la grille
!          10-06-09  Ameliorer le STOP des piscines
! 2009.2   04-09-09  format fichier sortie
! 2009.3   05-10-09  - ajout d'un "ifdef parallele"
!                    - array syntax remplacee par boucles
!          08-10-09  compatibilite avec f95: adequation des longueurs de
!                    chaines de charactere passees en argument de allocate_global
!          21-10-09  routine x_to_xyr renommee maskz_to_maskxyr et obc sur maskz
!                    appelees depuis maskz_to_maskxyr
!          04-11-09  ecrire les fichiers temporaire dans repertoire tmp
! 2010.2   20-12-09  subroutine bathycote renommee initial_mask_and_bathy
! 2010.3   12-01-10  lecture wetdry_cst2 wetdry_cst3
! 2010.5   01-02-10  interp_gebco est un interpolateur en ligne de la bathymetry
! 2010.6   02-02-10  renomme lon_t lat_t
! 2010.7   01-03-10  Ajout d'une routine produisant un fichier netcdf
!          05-03-10  finalisation du point precedent
! 2010.8   07-05-10  attribut axis remplace par content
!          09-05-10  seul le proc 0 ecrit interrupteur
! 2010.11  06-07-10  Produire un fichier de grille 2D COMODO
!          13-07-10  fichier topobug ecrit dans tmp et correction message d'erreur
!          15-07-10  "status" declare dans modume_principal
! 2010.13  14-10-10  Si masque rempli de zero alors stop
!          02-11-10  - plus d'info dans le fichier topo_bug
!                    - taille de boucles
! 2010.18  26-03-11  Debug rapport resolution modele/donnees
!                    une nouvelle subroutine (interp_gebco_japon)
! 2010.21  20-04-11  On n'arrete pas le run si le proc est vide
! 2010.22  26-04-11  Volume global de la grille au repos
!          30-04-11  Supprimer les bassins secondaires avec subroutine
!                    remove_secondary_bassins
!          10-05-11  used_unused_dom identifie les domaines mpi entierement masques
! 2010.25  05-03-12  reset biggest_bassin_
! S26      02-11-12  initialisation de la grille par fichier
!          23-06-13  zone mangrove
!          17-07-13  modif format fichier emodnet
!          19-07-13  h_inf_obc fait son entree dans notebook_bathy, h_sup sort
!          05-09-13  modif dans subroutine interp_emodnet
!          12-09-13  modifs notebook_tide
!          05-12-13  merge facile de la bathy avec celle de l'ogcm
!          06-12-13  seuil sur h_w pour eviter division par zero dans routine sigma_levels
!          19-03-14  ajout pour nemo offline
!          14-04-14  subroutine lissebathy cas 1, en argument, passer les indices globaux
!          21-04-14  amelioration filiere mask et bathy melanges avec ogcm
!          02-05-14  En cas de procedure offline, mask et h sont lu dans le fichier de
!                    de grille offline
!          26-06-14  Test sur sous-domaines actifs deplace
!          27-06-14  modif test aiguillage ioffline=2
!          05-07-14  ecriture d'un fichier listant les procs inutiles et designant
!                    les procs actifs charges de combler les trous des fichiers netcdf
!          09-07-14  suite point precedent
!          12-07-14  suite du point precedent
!          14-07-14  goto 126 devient goto 127
!          30-07-14  Meilleur partage par les procs actifs de la tache de comblement
!                    des zones occupees par les procs passifs
!          06-08-14  prevoir cas npasdom_=0
!          21-10-14  deplacement de l'appel a la routine de merge des bathy S et ogcm
!          13-12-14  melange avec la bathy de l'ogcm place apres lissebathy
!          01-02-15  ajout call obc_h(0)
!          19-04-15  Pourvoir lire glob_mask & glob_h depuis un fichier grid.nc
!          16-05-15  verification de la coherence de notebook_vertcoord et notebook_bathy
!          19-05-15  format lecture bathycote_in.dat
!          25-05-15  h_w=max(h_w,h0_w)
!          02-06-15  Si bathycote_in.dat='none' initiale h et mask avec h1d et 1
!          13-06-15  amelioration du cas bathycote au format netcdf
!          14-06-15  ajout d'un format case 3
!          07-07-15  Case 5 (bassin simplifie) enrichi
!          14-07-15  detection parametrage de notebook incoherent
!          21-07-15  Si cas VST il faut lire la bathy avant le calcul du pas de temps
!                    et autres variables dependant de h_w
!          07-11-15  fichier contenant la difference entre h lisse et h brute
!          09-01-16  amenagements sur case academique 5
!          12-02-16  lire hrmax dans notebook_bahy
!          16-02-16  format int autorise
!          03-04-16  flag_smooth_h_mask lu dans notebook_bathy
!          09-04-16  possibilite (si on enleve les commemtaires) de prendre en
!                    le rapport de frequence (inertie/onde) dans les vitesses
!                    de phase de la condition radiative barotrope (pour cas
!                    academiques par exemple)
!          29-06-16  possibilite d'appeler lissage cas 3 sur une fenetre i,j reduite
!          05-09-16  Filiere bouchage des trous, cas plus de trous que de boucheurs
!          10-09-16  Possibilite de lire le fichier grille de roms
!          16-01-17  upw_hrange1 et upw_hrange2 pour definir une zone d'advection upwind
!          16-02-17  Dans le cas d'une modelisation 1DV forcee par NEMO la bathy est 
!                    donnee par l'OGCM
!          07-03-17  Mise A jour canaux
!          02-04-17  Reset used_unused_dom
!          17-05-17  Pour contingence mpi la construction du reservoir utilise le tableau
!                    global glob_h. Celui-ci doit donc contenir la bathy post-traitee
!          28-05-17  flag_remove_secondary_bassins permet de choisir d'enlever les bassins
!                    secondaire ou non (si deja fait dans le bathy_maker)
!          01-06-17  ajout d'un correcteur de bathymetrie apres lissage
!          11-06-17  message ecran
!          19-06-17  routine remove_secondary_basin parallelisee
!          09-09-17  mask est kind=1
!          26-10-17  nouvelles dimensions pour sqr_hover
!          22-01-18  fichier d'entree format ijhlonlatmask
!          15-04-18  Obc treatment: if mask(i=2)=0 then mask(i=1)=0 etc...
!          05-05-18  - procedure arret d'urgence occigen
!                    - nouvelle version du bouchage des piscines
!                    - ajout flag1_smooth_h_mask flag3_smooth_h_mask
!          08-05-18 rendue compatibe avec glob_mask(kind=1)
!          08-05-18 subroutine load_mask_h_from_a_netcdffile rendue compatibe avec 
!                   glob_mask(kind=1) et glob_h real
!          09-05-18  - procedure arret d'urgence occigen
!          14-05-18 ajout call atlas
!          16-05-18  modif format d'ecriture
!          12-07-18  ajout mangrove_scheme
! v245     01-02-19  sqr_goverh_v(i,2)=sqrt(grav/max(h_v(i,jmax),wetdry_cst2)) 
! v246     14-02-19  ecriture du nom du fichier bathy landsea mask dans slurm
! v250     17-03-19  ajout sqr_goverh_u
! v252     14-04-19  call grid_webcanals_readfiles('hm')
! v257     18-06-19  Verifier la coherence du masque des canaux: call grid_webcanals_consistency
! v265     30-10-19  algo correction de la bathy avec meilleure prise en compte
!                    des chevauchements de correction
! v269     03-12-19  suppression sqr_goverh_v
! v270     12-12-19  modif debug
! v274     10-02-20  seul rank0 lit le fichier bathy mask puis bcast
! v276     12-03-20  plus de noms reconnus dans les fichiers netcdf
! v278     17-04-20  rankhmin rankhmax
! v287     17-07-20  utiliser tmpdirname
! v342     01-04-22  if(h_inf_obc==-9999.)h_inf_obc=h_inf !01-04-22
!          02-04-22  on etend de 1 point supplementaire la condition
!                    limite sur le land sea masque
! v347     28-04-22  Application des bornes hinf et hsup aprEs le lissage
!          17-05-22  flag_bathy_update et subroutine update_bathy pour faire varier la bathymetry
! v362     04-01-23  correction d'un bug dans la detection des piscines
! v364     16-01-23  Les perturbations de bathymetrie sont ajoutees via une somme ponderee exponentielle
! v365     31-01-23  retour A la version 362 mais avec un calcul plus precis de la distance
! v367     09-03-23  gerer le cas x10=0 dans initial_mask_and_bathy_corrector
!...............................................................................

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Land sea mask initialization
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      namelist/notebook_bathy/texte250,ioption,h_inf,h_inf_obc,rmax_ & !19-07-13
                             ,nsmooth,wetdry_cst1,wetdry_cst2        &
                             ,wetdry_cst3,mangrove_file_name         &
                             ,coef_diss_mangrove                     &
                             ,linear_coef_mangrove                   &
                             ,mergebathy_filename                    & !05-12-13
                             ,mergebathy_sponge,h1d,hrmax            & !12-02-16
                             ,flag1_smooth_h_mask                    & !05-05-18
                             ,flag3_smooth_h_mask                    & !05-05-18
                             ,flag_remove_secondary_bassins          & !28-05-17
                             ,mangrove_scheme                        & !12-07-18
                             ,flag_bathy_update !17-05-22

      namelist/notebook_bathy2/texte250 !01-06-17

      call allocate_global('a','glob_mask           ',iglb,jglb,0) !08-10-09
      call allocate_global('a','glob_h              ',iglb,jglb,0) !08-10-09

!.................
! lecture du namelist "notebook_bathy" 
      h_sup=99999. ; h1d=200. ;  mangrove_file_name='none'
      mergebathy_filename='nomerge' ; mergebathy_sponge=0 ; hrmax=20. !12-02-16
      flag_remove_secondary_bassins=1
      texte250='none'
      h_inf_obc=-9999. !01-04-22
      open(100,file=nomfichier(3)) ! lecture du namelist "notebook_bathy" !26-06-13
      read(100,nml=notebook_bathy)
      close(100)
      if(h_inf_obc==-9999.)h_inf_obc=h_inf !01-04-22
      if(mangrove_file_name=='none')coef_diss_mangrove=0.
!.................


!................................................................
! Particular case where the grid is loaded from an external file:
      if(ioffline==2) then !grdgrdgrdgrdgrdgrdgrdgrd> !27-06-14

! Case where the grid is loaded from a NEMO grid file
       if(initialgridfile_txt/='none') then ! nemo case nemo case >
        call initial_mask_h_from_file       ! nemo case
        goto 127
       endif                                ! nemo case nemo case >

! Case where the grid is loaded from a S grid file
       call grid_mask_h_offlinecase
       goto 126

      endif                !grdgrdgrdgrdgrdgrdgrdgrd>

!................................................................
! Standard case where the grid is computed from the notebook_grid
! polar/bipolar grid parameters:
! Cas du model_1D:
      if(flag3d.eq.0)then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        do i=0,imax+1
        do j=0,jmax+1
           mask_t(i,j,kmax+1)=1
              h_w(i,j)=h1d
        enddo
        enddo

      else             !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>             !19/01/07
! Cas des model_s 2D ou 3D:

      if(par%rank==0)write(6,'(a,a)')'sur le point de lire: ',texte250 !30/11/07

! Land/Sea Mask + Bathymetry File Loading:

! Before choosing default usual file format Case4, file name is
! analysed and special cases (if any) are applied instead.

      flag_=0 ; len_=len(trim(texte250))

      case_=4 ! default case

! Priority cases first:
! Case 6: sigstepfile
      if(flag_==0.and.ihybsig==2) then
        case_=6 ; flag_=1 
        if(ioption==1)stop 'ihybsig=2 & ioption=1 inconsistent choice'
      endif

!Case1: online interpolation of a GEBCO file:
      if(flag_==0.and.index(texte250,'gebco.txt')/=0) then
          case_=1 ; flag_=1
      endif          

!Case2: netcdf file:
      if(flag_==0.and.texte250(len_-2:len_)=='.nc')  then
         case_=2  ; flag_=1
      endif

!Case3: txt format, LSM followed by i,j,h:   !14-06-15
!     if(flag_==0.and.texte250(len_-3:len_)=='.ijh') then
      if(flag_==0.and.index(texte250,'.ijh')/=0) then !22-01-18
        case_=3  ; flag_=1
      endif

!Case4: usual default txt file (bathycote_in.dat):
      if(flag_==0.and.texte250/='none') then
        case_=4 ; flag_=1
      endif

!Case5: homogeneous case:
      if(flag_==0.and.texte250=='none') then
        case_=5 ; flag_=1
      endif

!.............................................................



!Case1: online interpolation of a GEBCO file:
      if(case_==1) then
         call interp_emodnet 
      endif

!Case2: netcdf file:
      if(case_==2) then
       call load_mask_h_from_a_netcdffile 
      endif

!Case3: txt format, LSM followed by i,j,h:   !14-06-15
      if(case_==3) then !--->
       if(par%rank==0) then !ooo> !10-02-20
       write(6,'(a,a)')'read bathy and land sea mask in ' ,trim(texte250) !14-02-19
       open(unit=3,file=texte250,recl=10000)
        do i=1,iglb
        read(3,'(100000i1)')(glob_mask(i,j),j=1,jglb)
        enddo
        do j=1,jglb
        do i=1,iglb
!        read(3,'(2i6,1x,f10.4)')i1,j1,glob_h(i,j)
                        read(3,*)i1,j1,glob_h(i,j)
!        if(i/=i1)stop 'Err Gus'
!        if(j/=j1)stop 'Err Tav'
        enddo
        enddo
       glob_mask(951:952,183)=0   ! Kasos Karpathos
! Dardanelles et alentours
       j=296
       do i=1028,1050
        glob_h(i,j)=50.
       enddo
       do i=1048,1054
       do j=292,296
        glob_h(i,j)=max(glob_h(i,j),50.)
       enddo
       enddo
       do i=1026,1029
       do j=297,299
        glob_h(i,j)=max(glob_h(i,j),50.)
       enddo
       enddo
! sortie Bosphore (traité comme fleuve sur 2 points)
       glob_h(1095,269:270)=50.

       close(3)


! Pour test de sensibilite: test sur Gibraltar origine
!       open(unit=3,file='espartel_inorigin')
!       if(par%rank==0)then
!       open(unit=3,file='espartel_19042019')
!        open(unit=3,file='espartel_20042019_etroit')
         open(unit=3,file='espartel_30042019')    ! utile pour TIDE8 et autres avant?
!        open(unit=3,file='espartel_04052019')    ! utile pour TIDE3 --> TIDE7
         do j=575,555,-1
!!                       write(3,'(i3,40(1x,f5.0))')j,(glob_h(i,j),i=68,98)
                         read(3,*)j1,(glob_h(i,j),i=68,98)
              if(j1/=j)stop 'concord pas'
         enddo
         close(3)
!       stop
!       endif




! ecriture Sicile
!      if(par%rank==0)then
!       open(unit=3,file='Sicily')
!       do j=486,464,-1
!                       write(3,'(i3,40(1x,f5.0))')j,(glob_h(i,j),i=642,661)
!                       read(3,*)j1,(glob_h(i,j),i=642,661)
!            if(j1/=j)stop 'concord pas'
!       enddo
!       close(3)
!      if(par%rank==0)then
        open(unit=3,file='Sicily_modif2')
        do j=490,455,-1
!                       write(3,'(i3,40(1x,f5.0))')j,(glob_h(i,j),i=642,661)
                        read(3,*)j1,(glob_h(i,j),i=642,661)
             if(j1/=j)stop 'concord pas'
        enddo
        close(3)
!      if(par%rank==0)then
        open(unit=3,file='ionian_modif4')
        do j=362,342,-1
!                       write(3,'(i3,40(1x,f5.0))')j,(glob_h(i,j),i=695,720)
                        read(3,*)j1,(glob_h(i,j),i=695,720)
             if(j1/=j)stop 'concord pas'
        enddo
        close(3)
!     if(par%rank==0)then
        open(unit=3,file='ionian_modif5')
        do j=350,320,-1
!                       write(3,'(i3,40(1x,f5.0))')j,(glob_h(i,j),i=660,690)
                        read(3,*)j1,(glob_h(i,j),i=660,690)
             if(j1/=j)stop 'concord pas'
        enddo
       close(3)
!     if(par%rank==0)then
!       open(unit=3,file='SudSardaigne')
        open(unit=3,file='SudSardaigne_step2')
        do j=562,517,-1
!                      write(3,'(i3,70(1x,f5.0))')j,(glob_h(i,j),i=600,650)
                        read(3,*)j1,(glob_h(i,j),i=600,650)
             if(j1/=j)stop 'concord pas'
        enddo
        close(3)
!     if(par%rank==0)then
        open(unit=3,file='SudSardaigne_plusouest')
        do j=570,540,-1
!                      write(3,'(i3,70(1x,f5.0))')j,(glob_h(i,j),i=575,600)
                        read(3,*)j1,(glob_h(i,j),i=575,600)
             if(j1/=j)stop 'concord pas'
        enddo
        close(3)

! 27/06/2019 corrections diverses
! entrée canal de Sicile
       glob_h(715:718,342)=470. ; glob_h(716,343)=470. 
       glob_h(716:719,341)=470. ; glob_h(715,341)=470. 
       glob_h(714:717,340)=470. ; glob_h(716:717,339)=470.
! Messine
       glob_h(775,430:431)=90.
! Peloponnese Cythere
       glob_h(895,263)=200. ; glob_h(896,264)=200. ;glob_h(897,265)=200.
! Cythere Anticythere
       do i=890,891
       do j=249,257
       glob_h(i,j)=min(glob_h(i,j),190.)
       enddo
       enddo
! Kasos Karpathos
!      glob_h(951:952,183)=50.   !  masqué maintenant
! Rhodes Karpathos
       glob_h(967,188)=450.
! Rhodes Turquie
       glob_h(988,190:191)=350. ; glob_h(987,191)=350. 
       glob_h(986,191:192)=350. ; glob_h(984:985,192)=350.
! Otrante
       do i=899,906
       do j=444,451
       glob_h(i,j)=max(glob_h(i,j),730.)
       enddo
       enddo
       do j=452,453
       do i=900,907
       glob_h(i,j)=max(glob_h(i,j),730.)
       enddo
       enddo
       j=454
       do i=899,908
       glob_h(i,j)=max(glob_h(i,j),730.)
       enddo
       glob_h(899,452:453)=730.
       do j=447,453
       glob_h(899,j)=0.5*(glob_h(899,j)+glob_h(898,j))
       enddo
       do j=444,451
       glob_h(906,j)=0.5*(glob_h(906,j)+glob_h(907,j))
       enddo


!      stop
!     endif

! ecriture Gibraltar Espartel Camarinal
!      if(par%rank==0)then
!       open(unit=3,file='espartel_incorrig')
!       do j=575,555,-1
!                       write(3,'(i3,40(1x,f5.0))')j,(glob_h(i,j),i=68,98)
!                       read(3,*)j1,(glob_h(i,j),i=68,98)
!            if(j1/=j)stop 'concord pas'
!       enddo
!       close(3)
!      endif

! ecriture temporaire bathycorse
!      open(unit=3,file='../../../GLOBMED/BATHYMASK/bathy_corse_new')
!       do j=725,699,-1
!                       read(3,'(i3,40(1x,f5.0))')j1,(glob_h(i,j),i=726,758)
!       enddo
!      close(3)
!      open(unit=3,file='../../../GLOBMED/BATHYMASK/bathy_corse_new2')
!       do j=735,705,-1
!                       read(3,'(i3,40(1x,f5.0))')j1,(glob_h(i,j),i=726,758)
!       enddo
!      close(3)

! Eciture nouveau fichier
!      if(par%rank==0)then
!      open(unit=3,file='../../../GLOBMED/BATHYMASK/bathy_GLOBMED_curvgib_dard_reparCors2fois_Sicil_1120_865.ijh',recl=10000)
!       do i=1,iglb
!       write(3,'(100000i1)')(glob_mask(i,j),j=1,jglb)
!       enddo
!       do j=1,jglb
!       do i=1,iglb
!                       write(3,*)i,j,glob_h(i,j)
!       enddo
!       enddo
!      close(3)
!      write(6,*)'fin ecriture nouveau fichier'
!      stop
!      endif

      endif                                               !--->

!Case4: usual default txt file (bathycote_in.dat):
      if(case_==4) then !---->
       if(par%rank==0) then !ooo> !10-02-20
       write(6,'(a,a)')'read bathy and land sea mask in ' ,trim(texte250) !14-02-19
       open(unit=3,file=texte250,recl=10000)                            !07/07/07
        do i=1,iglb
        read(3,'(100000i1)')(glob_mask  (i,j),j=1,jglb)     !19-05-15
        enddo
        do i=1,iglb
        read(3,*)(glob_h  (i,j),j=1,jglb)                               !01/04/03
        enddo
       close(3)
       endif                !ooo> !10-02-20
       s_typevar=0
       call donne_le_type(glob_mask(1,1))
       if(s_typevar==0)stop 'Err 328 s_typevar non prevu'
       if(s_typevar/=1)stop 'Err 328 glob_mask ne serait pas int*1 ?'
       if(s_typevar==1)call mpi_bcast(glob_mask,size(glob_mask),mpi_integer1,0,par%comm2d,ierr) !10-02-20
       s_typevar=0
       call donne_le_type(glob_h(1,1))
       if(s_typevar==0)stop 'Err 329 s_typevar non prevu'
       if(s_typevar/=4)stop 'Err 329 glob_h ne serait pas real*4 ?'
       if(s_typevar==4)call mpi_bcast(glob_h,size(glob_h),mpi_real,0,par%comm2d,ierr) !10-02-20
      endif             !---->

!Case5: homogeneous case:
      if(case_==5) then !---->
        glob_mask=1 ; glob_h=h1d ; flag_=1   !02-06-15
        glob_mask(1,:   )=0 ; glob_mask(iglb,:)=0
        glob_mask(:,jglb)=0 ; glob_mask(:   ,1)=0
        if(iperiodicboundary) then !09-01-16
         glob_mask(1   ,2:jglb-1)=1 
         glob_mask(iglb,2:jglb-1)=1 
        endif
        if(jperiodicboundary) then !09-01-16
         glob_mask(2:iglb-1,1   )=1 
         glob_mask(2:iglb-1,jglb)=1 
        endif
      endif             !---->

      if(case_==6) then
       texte250=trim(sigstepgridfile) 
       call load_mask_h_from_a_netcdffile 
      endif

! Obc treatment: if mask(i=2)=0 then mask(i=1)=0 etc...!15-04-18!02-04-22
      do j=1,jglb
       if(glob_mask(3     ,j)==0)glob_mask(2     ,j)=0 !02-04-22
       if(glob_mask(2     ,j)==0)glob_mask(1     ,j)=0
       if(glob_mask(iglb-2,j)==0)glob_mask(iglb-1,j)=0 !02-04-22
       if(glob_mask(iglb-1,j)==0)glob_mask(iglb  ,j)=0
      enddo
      do i=1,iglb
       if(glob_mask(i,3     )==0)glob_mask(i,2     )=0 !02-04-22
       if(glob_mask(i,2     )==0)glob_mask(i,1     )=0
       if(glob_mask(i,jglb-2)==0)glob_mask(i,jglb-1)=0 !02-04-22
       if(glob_mask(i,jglb-1)==0)glob_mask(i,jglb  )=0
      enddo

! Remove secondary bassins:
      if(flag_remove_secondary_bassins==1)call remove_secondary_bassins     !28-05-17

! Add Water ways:
!     call set_rivers(2) !07-03-17

!......................
! Passer de glob_h a h_w
! Passer de glob_mask a mask_t
! Reporter les valeurs de surface de mask sur tous les niveaux verticaux
      do j=0,jmax+1                                                     !05-10-09
      do i=0,imax+1
            mask_t(i             ,j             ,:)=    &
       glob_mask  (i+par%timax(1),j+par%tjmax(1))
               h_w(i             ,j             )=      &
          glob_h  (i+par%timax(1),j+par%tjmax(1))
      enddo
      enddo

      call obc_h(0) !01-02-15

! Melange bathy locale et bathy ogcm (attention appeler apres remove_secondary_bassins)!21-10-14
!     if(     mergebathy_filename/='nomerge' &
!        .and.mergebathy_sponge/=0            )call ogcm_get_bathy_ogcm !05-12-13

! Sea Mangrove Mask (if any):
       if(coef_diss_mangrove>0.)call mangrove_init_mask !23-06-13

      endif            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>              !19/01/07

!                           / / /


      if(     mergebathy_filename/='nomerge' &
         .and.mergebathy_sponge/=0            )call ogcm_get_bathy_ogcm !13-12-14

! Dans le cas d'une modelisation 1DV forcee par NEMO la bathy est donnee par l'OGCM !16-02-17
         if(flag_1dv==1.and.iobc_ogcm==1)call ogcm_get_bathy_ogcm


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Smoothing & threshold options on h:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! rdv les options sur H:
      if(ihybsig==1.and.ioption==0) &
      stop 'ihybsig==1.and.ioption==0 inconsistent choice' !14-07-15
      if(ioption.eq.1.and.flag3d.ne.0) then                               !19/01/07


!                       / / /

!...................................................
! LISSAGE. DEBUT:
!!!!     if(nsmooth.ne.0)call lissebathy(1,0.,1,imax,1,jmax,nsmooth)
!      if(nsmooth.ne.0)call lissebathy(1,0.,1,iglb,1,jglb,nsmooth)      !14-04-14
!      if(nsmooth.ne.0)call lissebathy(1,0.,575,750,550,740,nsmooth)      !14-04-14
       if(nsmooth.ne.0)call lissebathy(1,0.,560,690,565,730,nsmooth)      !14-04-14
! LISSAGE. FIN.
!...................................................

!                       / / /                                          !21/12/07

!.................................................................
! TRAITEMENT PREVENTIF DES ERREURS DE TRONCATURE LIEES AUX FORTES
! PENTES BATHYMETRIQUES. DEBUT:

! APPLICATION DU CRITERE dH/2H < RMAX
! (voir article Beckman and Haidvogel JPO 1993 pp 1736-1753)
! (voir article Mellor & cie JAOT 1994 pp 1126-1134)
! RMAX est compris entre 0 et 1
! Valeurs conseillees entre 0.1 (severe) et 0.6 (bien cool)
! Par defaut prendre RMAX=0.2
      rmax_=2.                                    ! claude
      if(rmax_.lt.1.)then
!       call lissebathy(3,rmax_,0,0   ,0,0   ,0)
        call lissebathy(3,rmax_,1  ,iglb,1  ,jglb,0)      !29-06-16
!       call lissebathy(3,0.1  ,230,342 ,374,452 ,0)      
      endif
      if(ihybsig==1) then !>>>
        call initial_main_inconsistencies4(rmax_)  !16-05-15
        call lissebathy(1,0.,1,iglb,1,jglb,nsmooth)
        h_w=max(h_w,h0_w) !25-05-15
      endif               !>>>

!................................................................
! Bornes mini et maxi (depend de la distance aux OBC):
! Appel deplacE le !28-04-22 apres le lissage pour que ce dernier 
! ne modifie pas le bornage
      call initial_sponge_hmin 

! TRAITEMENT PREVENTIF DES ERREURS DE TRONCATURE LIEES AUX FORTES
! PENTES BATHYMETRIQUES. FIN.
!.................................................................

!                       / / /

! rdv les options sur H:
      endif

!...............
! Smoothed bathymetry corrector
      texte250='none'
      open(100,file=nomfichier(3)) !01-06-17
      read(100,nml=notebook_bathy2)
      close(100)
      if(texte250/='none')call initial_mask_and_bathy_corrector

      endif    ! fin option deplacée ici
!     if(par%rank==0)write(6,*)'85 565 sortie de endif option ',glob_h(85,565)

! ecriture Gibraltar Espartel Camarinal
!       open(unit=3,file='espartel_out')
!       i1=78-par%timax(1) 
!       do j=559,567
!        j1=j-par%tjmax(1)
!       if(i1>0.and.i1<imax+1.and.j1>0.and.j1<jmax+1) &
!           write(3,*)i1+par%timax(1),j1+par%tjmax(1),h_w(i1,j1),glob_h(i1+par%timax(1),j1+par%tjmax(1))
!        enddo
!       i1=85-par%timax(1)
!       do j=561,569
!       j1=j-par%tjmax(1)
!       if(i1>0.and.i1<imax+1.and.j1>0.and.j1<jmax+1) &
!           write(3,*)i1+par%timax(1),j1+par%tjmax(1),h_w(i1,j1)
!       enddo
!       close(3)



!...............
! Archiver l'ecart A la bathy brute dans un fichier:
      write(texte30,'(a,i0)')trim(tmpdirname)//'err_h',par%rank !17-07-20
      open(unit=3,file=trim(texte30)) !07-11-15
       do j=0,jmax+1 ; do i=0,imax+1
        write(3,*)h_w(i             ,j             )   &
            -glob_h  (i+par%timax(1),j+par%tjmax(1))
       enddo         ; enddo
      close(3)

!................
! Add Water ways:
! Pour contingence mpi la construction du reservoir utilise le tableau
! global glob_h. Celui-ci doit donc contenir la bathy post-traitee, ce
! qui entraine la fonction par_gatherall_2d en suivant
      do j=1,jmax ; do i=1,imax
       glob_h(i+par%timax(1),j+par%tjmax(1))=h_w(i,j)
      enddo ; enddo
      ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
      ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
      lb2=lbound(glob_mask) ;  ub2=ubound(glob_mask)
      call  par_gatherall_2d(glob_h,lb2,ub2,ind,par%nbdom) !17-05-17
      call set_rivers(2) !07-03-17
!................
! Web of canals : lire la bathy des canaux
      if(nbcanal>0)call grid_webcanals_readfiles('hm') 

! Eviter division par zero dans routine sigma_levels: !06-12-13
      do j=1,jmax
      do i=1,imax
       if(h_w(i,j)>=0.) then
        h_w(i,j)=max(h_w(i,j), 1.d-3)
       else
        h_w(i,j)=min(h_w(i,j),-1.d-3)
       endif
      enddo
      enddo


! Conditions aux limites pour la bathy:                                !09/03/09
! en i=imax+1, i=0, j=0, j=jmax+1:
      call obc_h(0)

 126  continue

! Check the compatibility of the grid points river sources with the land-sea mask:
      call set_rivers(1)

 127  continue !15-07-14

! Remplissage des points de grilles peripheriques de H_Z
! Initialisation de la bathymetrie des vitesses et du point
!          "rotationnel" e partir de H_Z
      call hz_to_hxyr                                                  !24/03/04

! sqrt( H / grav) pour obc_ext: !17-07-13
! frequence de l'onde:
!     x0=2.*pi/43200. ! omega !09-04-16
      do i=1,imax ! 2,imax-1 !26-10-17
       sqr_hoverg_v(i,2)=sqrt(max(h_v(i,jmax),wetdry_cst2)/grav) &
!               *sqrt( 1.-(coriolis_t(i,jmax)/x0)**2 )           & !09-04-16
                        *obc2dtype ! 1 if standard 0 if clamped !12-09-13
       sqr_hoverg_v(i,1)=sqrt(max(h_v(i,2   ),wetdry_cst2)/grav) &
!               *sqrt( 1.-(coriolis_t(i,2   )/x0)**2 )           & !09-04-16
                        *obc2dtype ! 1 if standard 0 if clamped !12-09-13
      enddo
      do j=1,jmax ! 2,jmax-1 !26-10-17
       sqr_hoverg_u(j,2)=sqrt(max(h_u(imax,j),wetdry_cst2)/grav) &
!               *sqrt( 1.-(coriolis_t(imax,j)/x0)**2 )           & !09-04-16
                        *obc2dtype ! 1 if standard 0 if clamped !12-09-13
       sqr_hoverg_u(j,1)=sqrt(max(h_u(2   ,j),wetdry_cst2)/grav) &
!               *sqrt( 1.-(coriolis_t(2   ,j)/x0)**2 )           & !09-04-16
                        *obc2dtype ! 1 if standard 0 if clamped !12-09-13
      enddo

!...............................................................
! DEBUGAGE & CALCUL DE HMAX. DEBUT:
      x1=1.e10
      hmax=0.
      hmin=1.e10
      do 245 i=0,imax+1                                           !02-11-10
      do 245 j=0,jmax+1
      x1=min(x1,h_w(i,j))
         if(h_w(i,j).gt.hmax.and.  mask_t(i,j,kmax+1).eq.1) then  !18/09/01
           hmax=h_w(i,j)
           ihmax=i
           jhmax=j
         endif
         if(h_w(i,j).lt.hmin.and.  mask_t(i,j,kmax+1).eq.1) then  !03/12/01
           hmin=h_w(i,j)
           ihmin=i
           jhmin=j
         endif
  245 continue

      if(x1.lt.0.001) then !---------------------------------->
         open(unit=3,file=trim(tmpdirname)//'topo_bug'//dom_c//'.txt')              !13-07-10
         write(3,*)'legende: i j (loc) i j (glob) mask_t h_w'
         do 592 i=0,imax+1
         do 592 j=0,jmax+1
         if(h_w(i,j).lt.0.001)write(3,'(5(i4,1x),f12.5)')              &
                                            i,j,i+par%timax(1)         &
                                           ,j+par%tjmax(1)             &
                                           ,mask_t(i,j,kmax+1),h_w(i,j) !02-11-10
!        if(h_w(i,j).lt.0.001)write(3,*)i,j,  mask_t(i,j,kmax+1),h_w(i,j)
  592    continue
         close(3)

      endif               !----------------------------------->

      call mpi_allreduce(hmin,hminglob,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     hmin=hminglob
      call mpi_allreduce(hmax,hmaxglob,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     hmax=hmaxglob

! Trouver le rank du hmin: !17-04-20
      k0=-999
      if(hmin==hminglob)k0=par%rank
! Au cas oU ils seraient plusieurs, choisir le plus grand des ranks (mpi_max):
      call mpi_allreduce(k0,rankhmin,1,mpi_integer,mpi_max,par%comm2d ,ierr)

! Trouver le rank du hmax: !17-04-20
      k0=-999
      if(hmax==hmaxglob)k0=par%rank
      call mpi_allreduce(k0,rankhmax,1,mpi_integer,mpi_max,par%comm2d ,ierr)

      hmin=hminglob
      hmax=hmaxglob
      call mpi_bcast(ihmin,1,mpi_integer,rankhmin,par%comm2d,ierr)
      call mpi_bcast(jhmin,1,mpi_integer,rankhmin,par%comm2d,ierr)
      call mpi_bcast(ihmax,1,mpi_integer,rankhmax,par%comm2d,ierr)
      call mpi_bcast(jhmax,1,mpi_integer,rankhmax,par%comm2d,ierr)

! DEBUGAGE & CALCUL DE HMAX. FIN.
!...............................................................

!                       / / /

!...............................................................
! Dans le cas d'une imbrication type IARCHIVE=1
! on peut verifier que la zone s'imbrique bien en comparant
! le fichier toto avec le fichier bathycote_in.dat du petit
! domaine: il n'y doit pas y avoir de differences!
      if(iarchive.eq.1)stop 'filiere iarchive desactivee'

!                       / / /

!_____________________________________________________________________!
! ETAPE FINALE SUR LE MASQUAGE                                         !05/10/01
! Debut:
! remplissage des points de grille _x _y _z
      call maskt_to_maskuvp ! obc for mask_t & computes mask_u _y _r   !21-10-09

! Fin.
! ETAPE FINALE SUR LE MASQUAGE                                         !05/10/01
!_____________________________________________________________________!

!                       / / /

      if(par%rank==rankhmax) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine bathycote:'
      write(3,*)
      write(3,*)'maximum depth hmax=',hmax,' meters'
      write(3,*)'corresponding rank:',rankhmax
      write(3,*)'corresponding grid point loc:',ihmax,jhmax
      write(3,*)'corresponding grid point glb:',ihmax+par%timax(1),jhmax+par%tjmax(1)
      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif
      if(par%rank==rankhmin) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)
      write(3,*)'minimum depth hmin=',hmin,' meters'
      write(3,*)'corresponding rank:',rankhmin
      write(3,*)'corresponding grid point loc:',ihmin,jhmin
      write(3,*)'corresponding grid point glb:',ihmin+par%timax(1),jhmin+par%tjmax(1)
      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

!                       / / /

! Verifier la coherence du masque des canaux !18-06-19
      call grid_webcanals_consistency

!..............................................................................
! Debug: Verifier qu'il n'y a pas de piscine cachee dans le domaine:  !28/11/08
      flag_stop=0
      do j=1,jmax
      do i=1,imax
       ip1=min0(i+1,imax)
       jp1=min0(j+1,jmax)
       im1=max0(i-1,1)
       jm1=max0(j-1,1)
        if(  mask_t(i,j,kmax+1).eq.1) then
         ksecu=0
         if(  mask_t(ip1,j  ,kmax+1).eq.1)ksecu=1
         if(  mask_t(im1,j  ,kmax+1).eq.1)ksecu=1
         if(  mask_t(i  ,jp1,kmax+1).eq.1)ksecu=1
         if(  mask_t(i  ,jm1,kmax+1).eq.1)ksecu=1
         if(ksecu.eq.0) then
          flag_stop=1 !12-12-19
          write(10+par%rank,*) &
          'piscine par%rank,i,j,i+par%timax(1),j+par%tjmax(1)'    &
                  ,par%rank,i,j,i+par%timax(1),j+par%tjmax(1) !11-06-17
         endif
        endif
      enddo
      enddo
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) !05-05-18
      if(k0/=0) &
      stop 'Err: secondary bassin detected. See errors fort.xx files' !09-05-18

! Compute grid area & volume:           !26-04-11
      sum1=0.
      sum2=0.                           !26-04-11
      do j=1,jmax
      do i=1,imax
       sum1=sum1+mask_t(i,j,kmax+1)*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j) !ne pas sommer 2 fois zone de recouvrement !04-06-09
       sum2=sum2+mask_t(i,j,kmax+1)*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)*h_w(i,j)
      enddo
      enddo
        grid_area=sum1
#ifdef parallele
      call mpi_allreduce(grid_area                            &
                        ,grid_areaglb,1,mpi_double_precision, & !05-10-09
                         mpi_sum,par%comm2d,ierr)
      grid_area=grid_areaglb
      call mpi_allreduce(sum2                                   & !26-04-11
                        ,grid_volumeglb,1,mpi_double_precision, &
                         mpi_sum,par%comm2d,ierr)
#else
      grid_volumeglb=sum2
#endif

! Tant que glob_mask est toujours allouE appeler subroutine atlas
      call atlas !14-05-18

      call allocate_global('d','glob_mask           ',0,0,0)       !08-11-09
      call allocate_global('d','glob_h              ',0,0,0)       !08-11-09

! Analyse the mpi distribution and detect entirely masked subdomains.
! Propose a new distribution based on 100% active subdomains
! This new distribution is not currently applied but only stored in an
! description file that can be used for the next runs
      if(nbdom==nbdom_imax*nbdom_jmax) then ! m[0_0]m >
         call discard_unused_mpidom ! this routine is in the present file
      else                                  ! m[0_0]m >
         used_unused_dom=1 !02-04-17
      endif                                 ! m[0_0]m >


! Ajouter des criteres de bathymetrie dans le calcul de upwindriver_t
! Ici faire en sorte que les zones tres peu profondes suceptibles d'etre
! des zones intertidales soient "zone upwind"
      do j=0,jmax+1 ; do i=0,imax+1 

       upwindriver_t(i,j)=min(upwindriver_t(i,j)          & ! Prise de compte du critere precedent defini dans set_river(case=0)
       ,(h_w(i,j)-upw_hrange1)/(upw_hrange2-upw_hrange1))   ! Critere bathy 100% UPW si H<upw_hrange1 !16-01-17 

       upwindriver_t(i,j)=min(max(upwindriver_t(i,j),0.),1.) !  0<upwindriver_t<1

      enddo         ; enddo
      

      end subroutine initial_mask_and_bathy !21-06-14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine discard_unused_mpidom !21-06-14
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer,dimension(:),allocatable :: glb_used_unused_in  &
                                         ,glb_used_unused_out &
                                         ,new_rank
!                                        ,hole_per_rank
!     integer,dimension(:,:),allocatable :: hole_index
      integer nactdom_,npasdom_,n_act_ove_pas_,nrepldom1_,nrepldom2_ &
             ,istr_,iend_,jstr_,jend_
#ifdef synopsis
       subroutinetitle='discard_unused_mpidom'
       subroutinedescription= &
         ' Identifies the 100% inland subdomains and suggests an'       &
       //' alternative mpi map excluding unemployed subdomains'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      allocate(glb_used_unused_in (0:nbdom-1))      ; glb_used_unused_in=0
      allocate(glb_used_unused_out(0:nbdom-1))      ; glb_used_unused_out=0
      allocate(new_rank           (0:nbdom-1))      ; new_rank=-2
!     allocate(hole_per_rank      (0:nbdom-1))      ; hole_per_rank=0
!     allocate(hole_index         (0:nbdom-1,nbdom)); hole_index=0
! deallocate a la fin de cette subroutine...


      k0=0
      do j=0,jmax+1
      do i=0,imax+1
      k0=k0+mask_t(i,j,kmaxp1)
      enddo
      enddo

!>>> Check the number of sea node in the domain:
!     write(6,*)'Par%rank & Number of sea nodes=',par%rank,k0
      if(k0==0) then !wwww>                                          !14-10-10
       write(6,*)'WARNING: no sea in domain with rank=',par%rank
       used_unused_dom=0                                                !10-05-11
      else              !wwww>
       if(k0<0) stop ' k0<0 in discard_unused_mpidom'
       write(6,*)'Par%rank & Number of sea grid points=',par%rank,k0
       used_unused_dom=1                                                !10-05-11
      endif             !wwww>
!......................

      glb_used_unused_in(par%rank)=used_unused_dom

      call mpi_allreduce(glb_used_unused_in     &
                        ,glb_used_unused_out    &
                        ,nbdom,mpi_integer,mpi_sum,par%comm2d,ierr)

! glb_used_unused_out est global. Il est permet a chaque rank de savoir le
! status actif ou non actif de tous les autres rank
! new_rank est le futur nouveau numero de rank du rank actuel
       k0=-1
       do loop1=0,nbdom-1
        if(glb_used_unused_out(loop1)==1) then
          k0=k0+1
          new_rank(loop1)=k0
        endif
       enddo
! Au sortir de cette boucle k0+1 est le futur nombre total de sous domaines (actifs)

      nactdom_=k0+1          ! Number of ACTive DOMains !30-07-14
      npasdom_=nbdom-nactdom_! Number of PASsive DOMains
! Les doms actifs vont se partager les noeuds passifs dont ils rempliront l'espace fichier netcdf
! Cependant le ratio actifs/passifs ne tombe pas exactement sur une valeur entiere. Le partage
! ne sera donc pas tout a fait homogene. Certains procs passifs se verront remplaces par N
! proc actifs et d'autre par N+1 procs actifs, la somme donnant le nombre total de proc actif
!     N=nactdom_/npasdom_   ! Arrondi a la valeur inferieure du nbre d'actif remplacant un proc passif
! La difference entre npasdom_*N et nactdom_ donne le nombre de proc passif qui seront remplaces
! par N+1 procs actifs
      if(npasdom_==0) goto 580 !06-08-14

      n_act_ove_pas_=nactdom_/npasdom_ ! N
      nrepldom2_=nactdom_-npasdom_*n_act_ove_pas_
      nrepldom1_=npasdom_-nrepldom2_

! Seul le rank 0 cree le nouveau fichier description
      if(par%rank==0) then !000000000000>

       open(unit=3,file=trim(tmpdirname)//'description_domaine.txt') !09-07-14
       open(unit=4,file='description_domaine.next')
         read(3,*)i,j
        write(4,'(3i6,a60)')i,j,k0+1, &
       ' ! Number of sub-domains in each direction & nbdom'
         read(3,*)i,j
        write(4,*)i,j,' ! iglb jglb'

        do loop1=0,nbdom-1

        if(glb_used_unused_out(loop1)==0) then !>>>>>>>
         do k1=1,7 ; read(3,*) ; enddo
        else                                   !>>>>>>>
          read(3,*)
         write(4,'(a)')'------------------------'
          read(3,*)
         write(4,*)new_rank(loop1) &
                 ,'             ! sub-domain order number'

          read(3,'(a)')texte90
         write(4,'(a)')trim(texte90)
          read(3,'(a)')texte90
         write(4,'(a)')trim(texte90)
          read(3,'(a)')texte90
         write(4,'(a)')trim(texte90)
          read(3,'(a)')texte90
         write(4,'(a)')trim(texte90)

          read(3,*)i1,i2,i3,i4,i5,i6,i7,i8

         if(i1/=-2)i1=new_rank(i1)
         if(i2/=-2)i2=new_rank(i2)
         if(i3/=-2)i3=new_rank(i3)
         if(i4/=-2)i4=new_rank(i4)
         if(i5/=-2)i5=new_rank(i5)
         if(i6/=-2)i6=new_rank(i6)
         if(i7/=-2)i7=new_rank(i7)
         if(i8/=-2)i8=new_rank(i8)

         write(4,'(8(i6,1x),a40)')             & !16-05-18
                   i1,i2,i3,i4,i5,i6,i7,i8     &
               ,' ! Neighbors: w e n s ws es wn en'

        endif                                  !>>>>>>>


        enddo

       close(3)
       close(4)

       write(6,*)'Number of active  doms=',nactdom_
       write(6,*)'Number of passive doms=',npasdom_
       write(6,*)'int(actives/passives) =',n_act_ove_pas_
       write(6,*)'and concerned doms    =',nrepldom1_
       write(6,*)'actives/passives+1    =',n_act_ove_pas_+1
       write(6,*)'and concerned doms    =',nrepldom2_

       open(unit=3,file='description_trous.txt')
       if(npasdom_<nactdom_) then !case1>

       k0=0
       k2=nrepldom2_*(n_act_ove_pas_+1)
       k1=nrepldom1_* n_act_ove_pas_     ! reset compteurs des boucheurs
       do k=0,nbdom-1  ! > k >
        if(glb_used_unused_out(k)==0) then ! unused >

        write(6,*)'----------trou dom no ',k,par%gtjmax(k,:)

         if(k2/=0) then ! k2 k2 >
          j1=par%gtjmax(k,1)
          do k3=1,n_act_ove_pas_+1  !k3>
           j2=par%gtjmax(k,1)+real(k3)/real(n_act_ove_pas_+1)*(par%gtjmax(k,2)-par%gtjmax(k,1))
           write(6,*)'Boucheur no ',k0
           write(6,*)'Fraction    ',n_act_ove_pas_+1
           write(6,*)'Bornes      ',j1,j2
           write(3,*)'...........................'
           write(3,'(i6,a)')k0 &
        ,'      1  ! subdomain rank and number of related lost areas'
           write(3,'(2(i6,1x),a)')                    &
                  par%gtimax(k,1)                  &
                 ,par%gtimax(k,2)                  &
                 ,' ! i start i end (par%gtimax(1:2)) of the lost area'
           write(3,'(2(i6,1x),a)')  j1,j2                  &
                 ,' ! j start j end (par%gtjmax(1:2)) of the lost area'
           if(j2-j1<=2)flag=1

           k2=k2-1 ; k0=k0+1 ; j1=j2-2
          enddo                     !k3>
          goto 666
         endif          ! k2 k2 >

         if(k1/=0) then ! k2 k2 >
          j1=par%gtjmax(k,1)
          do k3=1,n_act_ove_pas_    !k3>
           j2=par%gtjmax(k,1)+real(k3)/real(n_act_ove_pas_)*(par%gtjmax(k,2)-par%gtjmax(k,1))
           write(6,*)'Boucheur no ',k0
           write(6,*)'Fraction    ',n_act_ove_pas_
           write(6,*)'Bornes      ',j1,j2
           write(3,*)'...........................'
           write(3,'(i6,a)')k0 &
        ,'      1  ! subdomain rank and number of related lost areas'
           write(3,'(2(i6,1x),a)')                    &
                  par%gtimax(k,1)                  &
                 ,par%gtimax(k,2)                  &
                 ,' ! i start i end (par%gtimax(1:2)) of the lost area'
           write(3,'(2(i6,1x),a)')  j1,j2                  &
                 ,' ! j start j end (par%gtjmax(1:2)) of the lost area'
           if(j2-j1<=2)flag=1
           k1=k1-1 ; k0=k0+1 ; j1=j2-2
          enddo                     !k3>
          goto 666
         endif          ! k2 k2 >

 666    continue
        endif                              ! unused >
       enddo           ! > k >

!      close(3)

! Renoncer au fichier des boucheurs de trous si les trous sont trop peu nombreux
! 645   if(flag==1) then !.........>
!       open(unit=3,file='description_trous.txt')
!       write(3,*)' Not enough holes. Files cancelled.'
!       close(3)
!      endif            !.........>

       else                       !case1> case2> !05-09-16

!     residu = nombre de trou - nombre de boucheur
!     residu decroit de un A chaque fois qu'un boucheur prend
!     en charge un trou en plus (en plus de 1)
!     residu_=npasdom_-nactdom_
!     write(6,*)'npasdom_',npasdom_

      k11=npasdom_
  925 continue
      k10=-1
      do k=0,nbdom-1 

! Identifier boucheurs:
       if(glb_used_unused_out(k)==1) then !boucheur>

! Chercher les trous
        do i1=0,nbdom_imax-1
        do j1=0,nbdom_jmax-1


! rank correspondant:
          k1=j1+i1*nbdom_jmax
         
          if(glb_used_unused_out(k1)==0) then !trous>
          k10=k10+1

! imin et imax du trou A boucher
                  istr_=par%gtimax(k1,1)
                  iend_=par%gtimax(k1,2)      
                  jstr_=par%gtjmax(k1,1)
                  jend_=par%gtjmax(k1,2)      

              glb_used_unused_out(k1)=2 ! Marquer ce trou comme Etant bouche
              k11=k11-1 ! compteur decroissant de trous restant A boucher
              goto 49 

          endif                               !trous>


        enddo !j1
        enddo  !i1

   49 continue

           write(3,*)'...........................'
           write(3,'(i6,a)')k10  &
        ,'         ! subdomain rank '
!       ,'      1  ! subdomain rank and number of related lost areas'
           write(3,'(2(i6,1x),a)')                    &
                  istr_                  &
                 ,iend_                  &
                 ,' ! i start i end (par%gtimax(1:2)) of the lost area'
           write(3,'(2(i6,1x),a)')  jstr_,jend_                  &
                 ,' ! j start j end (par%gtjmax(1:2)) of the lost area'

        if(k11==0)goto 926 ! sortir si tous les trous sont bouchEs

       endif                              !boucheur>

      enddo ! boucle sur les procs

!     k11 compteur decroissant de trous restant A boucher
!     s'il reste des trous repartir sur un tour supplementaire de
!     bouchage
       if(k11>0) then
          write(3,'(a)')'PERFORM'
          goto 925
       endif


       endif                             !case2>

 926   write(3,'(a)')'PERFORM'
       close(3)

!      write(6,*)'k11=',k11
!     stop 'jjtb'
      endif                !000000000000>

 580  continue

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif
      deallocate(glb_used_unused_in)
      deallocate(glb_used_unused_out)
      deallocate(new_rank)
!     deallocate(hole_per_rank)
!     deallocate(hole_index)

      end subroutine discard_unused_mpidom !21-06-14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef bidon
      subroutine interp_gebco
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      use module_global
      implicit none
      real lonmin_loc,lonmax_loc,latmin_loc,latmax_loc,res_loc
#ifdef synopsis
       subroutinetitle='interp_gebco'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      open(unit=3,file=texte250,recl=10000)                            !07/07/07

      read(3,*)
      read(3,*)
      read(3,*)max_x,max_y,   &
      lonmin_loc,latmin_loc,lonmax_loc,latmax_loc,res_loc
      res_loc=res_loc*deg2rad

      call allocate_forcages(1,7,max_x,max_y,1)

      write(*,*)max_x,max_y,  &
      lonmin_loc,latmin_loc,lonmax_loc,latmax_loc,res_loc

      do j=max_y,1,-1
      do i=1,max_x
       read(3,*)multi2d_lon(i,j),multi2d_lat(i,j),multi2d_var(i,j,1)
!      if(multi2d_var(i,j,1)<0)write(68,*)multi2d_lon(i,j),multi2d_lat(i,j)
      enddo
      enddo

      close(3)

      do j=1,jglb
      do i=1,iglb

! rapport d'echelle symphonie / gebco
      k=max0(1 , nint( sqrt(dxdy_t(i,j)) &
      /( rayonterre*res_loc*cos(lat_t(i,j)) ))) !26-03-11
      dist0=real(k)

! Indice decimale i (longitude) dans grille gebco:
              deci=1.+( lon_t(i,j)*rad2deg-lonmin_loc)              &
                              /(lonmax_loc-lonmin_loc)              &
                            *(max_x-1.)
! Indice decimale j (latitude) dans grille gebco:
              decj=1.+( lat_t(i,j)*rad2deg-latmin_loc)              &
                              /(latmax_loc-latmin_loc)              &
                            *(max_y-1.)

! On moyenne (k x k) points gebco autour de (i0,j0) point central
         i0=nint(deci)
         j0=nint(decj)
         sum1=small1
         sum2=0.
         do j1=max0(j0-k,1),min0(j0+k,max_y)
         do i1=max0(i0-k,1),min0(i0+k,max_x)

           x1=max(1.-sqrt(real((deci-i1)**2+(decj-j1)**2))/dist0,zero)
           sum1=sum1+x1
           sum2=sum2+x1*multi2d_var(i1  ,j1  ,1)

         enddo
         enddo

         glob_h(i,j)=sum2/sum1

         if(glob_h(i,j)>1.) then
            glob_mask(i,j)=1
         else
            glob_mask(i,j)=0
         endif

      enddo
      enddo

      call allocate_forcages(2,7,0,0,0)

! Supprimer les bassins secondaires:

      open(unit=3,file=trim(tmpdirname)//'piscines.txt')
      k=1
  618 k10=0
      sum1=1.
      k=k+1
      write(*,*)'Bassin No',k
  608 continue
      do j=1,jglb
      jp1=min0(j+1,jglb)
      jm1=max0(j-1,1)
      do i=1,iglb
       if(glob_mask(i,j)==1) then !>>>>>>
         if(k10==0)glob_mask(i,j)=k
         k10=1
         ksecu=0
         ip1=min0(i+1,iglb)
         im1=max0(i-1,1)
         if(glob_mask(ip1,j  )==k)ksecu=1
         if(glob_mask(im1,j  )==k)ksecu=1
         if(glob_mask(i  ,jm1)==k)ksecu=1
         if(glob_mask(i  ,jp1)==k)ksecu=1
         if(ksecu==1) then
             glob_mask(i,j)=k
             sum1=sum1+1
             goto 608
         endif
       endif                      !>>>>>>
      enddo
      enddo
      write(3,*)k,nint(sum1)
      if(k10==1)goto 618
      close(3)

     open(unit=3,file=trim(tmpdirname)//'piscines.txt')
       k3=0
       k4=0
  635  read(3,*,end=629)k1,k2
        if(k2>k3) then !---->
         k3=k2
         k4=k1
         goto 635
        endif          !---->
  629 close(3)

      write(*,*)'le plus grand bassin est le No',k4
      do j=1,jglb
      do i=1,iglb
        if(glob_mask(i,j)==k4) then
           glob_mask(i,j)=1
        else
           glob_mask(i,j)=0
        endif
      enddo
      enddo

! Archiver le resultat dans un fichier bathycote_in.dat
      write(*,*)'Archive dans tmp/bathycote_in.dat'
      open(unit=3,file=trim(tmpdirname)//'bathycote_in.dat',recl=10000)
        do i=1,iglb
        write(3,'(100000i1)')(glob_mask  (i,j),j=1,jglb)
        enddo
        do i=1,iglb
        write(3,*)(glob_h  (i,j),j=1,jglb)
        enddo
      close(3)

!     stop ' INTERPOLATION GEBCO TERMINEE'

      end subroutine interp_gebco
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef bidon
      subroutine save_grid_netcdf
!______________________________________________________________________
!
! S model
! release 2010.10  - last update: 24-06-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

!......................................................................
! Version date      Description des modifications:
! 2010.11 06-07-10  Application des recommandations COMODO
!......................................................................

      use module_principal
      use module_parallele !#MPI
      implicit none
      integer ichoix
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='save_grid_netcdf'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Storage of the 2D horizontal grid in a netcdf file


      filval=-9999.
      texte80(3)='none'
      texte80(4)='none'

! ECRIRE LE FICHIER DE GRILLE:

!......nom du fichier de grille
      texte250=trim(tmpdirname)//''//dom_c//'_'//'galeo_grid.nc'
      texte60='galeo_grid.nc'

      status=nf_create(texte250,nf_clobber,ncid)
      if(status/=0)stop ' stop save_grid_netcdf erreur 1'

      call netcdf_dim

! Reset VARDIM

      do loop_netcdf=0,1
      count_netcdfvar=0

      if(loop_netcdf==1) then !>>>>>>>>>>>>>>>>>>>>

      status=nf_put_att_text(ncid,nf_global,'production',48    &
       ,'Generalized ALE coordinate Ocean model')

      status=nf_put_att_text(ncid,nf_global,'Conventions',6   &
       ,'C0 grid')

!     status=nf_put_att_text(ncid,nf_global,'grid_type',55             &
!      ,'http://sirocco.omp.obs-mip.fr/images/c_grid_shift_4.png')

      status=nf_put_att_text(ncid,nf_global,'file_name',60,texte60)

! Definition of variables: done.
      status=nf_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lon_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_t'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lon_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_u'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lon_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_v'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lon_f(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_f'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_f')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lat_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_t'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_t')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lat_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_u'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_u')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lat_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_v'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
        anyv3d(i,j,1,1)=lat_f(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_f'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)='double'
      call netcdf_main('_f')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar2d(i,j)=h_w(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_t' ; texte80(2)='m'                     ! variable ; units
      texte80(3)='undisturbed water depth'                  ! long_name
      texte80(4)='undisturbed_water_depth'                  !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_u(i,j,kmax+1)==1) then
          anyvar2d(i,j)=h_u(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_u' ; texte80(2)='m'                     ! variable ; units
      texte80(3)='undisturbed water depth'                  ! long_name
      texte80(4)='undisturbed_water_depth'                  !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_v(i,j,kmax+1)==1) then
          anyvar2d(i,j)=h_v(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_v' ; texte80(2)='m'                     ! variable ; units
      texte80(3)='undisturbed water depth'                  ! long_name
      texte80(4)='undisturbed_water_depth'                  !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_f(i,j,kmax+1)==1) then
          anyvar2d(i,j)=h_f(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_f' ; texte80(2)='m'                     ! variable ; units
      texte80(3)='undisturbed water depth'                  ! long_name
      texte80(4)='undisturbed_water_depth'                  !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_f')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dy_t(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dy_t'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oj axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oj_axis'   !  standard_name
      texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_u(i,j,kmax+1)==1) then
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
       if(mask_v(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dy_v(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dy_v'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oj axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oj_axis'   !  standard_name
      texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_f(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dy_f(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dy_f'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oj axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oj_axis'   !  standard_name
      texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_f')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dx_t(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dx_t'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oi axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oi_axis'   !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_u(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dx_u(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dx_u'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oi axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oi_axis'   !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_v(i,j,kmax+1)==1) then
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
       if(mask_f(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dx_f(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dx_f'   ; texte80(2)='m'                    ! variable ; units
      texte80(3)='Horizontal resolution along Oi axis'   ! long_name
      texte80(4)='Horizontal_resolution_along_Oi_axis'   !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_f')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=mask_t(i,j,kmax+1)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_t'   ; texte80(2)='none'                   ! variable ; units
      texte80(3)='sea land mask'                                  ! long_name
      texte80(4)='sea_land_mask'                           ! standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=mask_u(i,j,kmax+1)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_u'   ; texte80(2)='none'                   ! variable ; units
      texte80(3)='sea land mask'                                  ! long_name
      texte80(4)='sea_land_mask'                           ! standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=mask_v(i,j,kmax+1)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_v'   ; texte80(2)='none'                   ! variable ; units
      texte80(3)='sea land mask'                                  ! long_name
      texte80(4)='sea_land_mask'                           ! standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=mask_f(i,j,kmax+1)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_f'   ; texte80(2)='none'                   ! variable ; units
      texte80(3)='sea land mask'                                  ! long_name
      texte80(4)='sea_land_mask'                           ! standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_f')

      enddo  ! fin de boucle sur loop_netcdf
      status=nf_close(ncid)

      end subroutine save_grid_netcdf
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef bidon
      subroutine interp_gebco_japon
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      use module_global
      implicit none
      double precision res_loc,h_loc,weight_loc
#ifdef synopsis
       subroutinetitle='interp_gebco_japon'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Reset des compteurs
      do j=1,jmax
      do i=1,imax
       xy_t(i,j,1)=0.
       xy_t(i,j,2)=0.
      enddo
      enddo

!     open(unit=3,file='../../../data/bathymetry/japon_gebco.txt',recl=10000)
      open(unit=3,file= &
      '../../../data/bathymetry/deux_bathy_du_japon.txt',recl=10000)
   99 read(3,*,end=100)x1,x2,h_loc,res_loc,weight_loc

      call latlontoij(x1*deg2rad,x2*deg2rad,'loc')
      i=nint(deci)          ; j=nint(decj)

      if(i<1   ) goto 99
      if(i>imax) goto 99
      if(j<1   ) goto 99
      if(j>jmax) goto 99

      i2=min(max(i,1),imax) ; j2=min(max(j,1),jmax)
! k = rapport de resolution "donnee/modele"
      rap=res_loc/sqrt(dxdy_t(i2,j2))
      k=nint(rap)
! Si k=0 la donnee ne se projete que sur une seule cellule (la plus proche)
! Ce cas ne se produit que si la grille des donnees est plus fine que celle du modele

!     write(67,*)k,i,j,sqrt( dxdy_t(i,j)),rayonterre*res_loc*deg2rad*cos(lat_t(i2,j2))


!     const1=(1./res_loc)**2
      do j1=j-k,j+k
       if(j1>=1.and.j1<=jmax) then !jjjjjjjjjj>
         do i1=i-k,i+k
          if(i1>=1.and.i1<=imax) then !iiiiiiiii>
            xy_t(i1,j1,1)=xy_t(i1,j1,1)+weight_loc
            xy_t(i1,j1,2)=xy_t(i1,j1,2)+weight_loc*h_loc
          endif                       !iiiiiiiii>
         enddo
       endif                      !jjjjjjjjjjj>
      enddo

      goto 99
  100 close(3)

      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmax+1)==1) then !1111111>
      if(xy_t(i,j,1)==0.)then
            write(*,*)i,j,xy_t(i,j,1),weight_loc
           stop 'Trou bathy'
      endif
!     write(66,*)i,j
!     if(xy_t(i,j,1)==0.)write(67,*)i,j
       x0=xy_t(i,j,2)/xy_t(i,j,1)
       if(x0<0)x0=h_w(i,j)
       h_w(i,j)=x0
      endif                        !1111111>
      enddo
      enddo

      write(*,*)'poids',xy_t(1,1,1)
      write(*,*)'poids',xy_t(imax/2,jmax/2,1)

! Archiver le resultat dans un fichier bathycote_in.dat
      write(*,*)'Archive dans tmp/bathycote_in.dat'
      open(unit=3,file=trim(tmpdirname)//'bathycote_in.dat',recl=10000)
        do i=1,imax
        write(3,'(100000i1)')(mask_t(i,j,kmax+1),j=1,jmax)
        enddo
        do i=1,imax
        write(3,*)(h_w(i,j),j=1,jmax)
        enddo
      close(3)

!     stop ' INTERPOLATION GEBCO TERMINEE'

      end subroutine interp_gebco_japon
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine initial_mask_h_from_file
      use module_principal
      use module_parallele !#MPI
      use pnetcdf
      implicit none
      integer ncid_
      integer(kind=MPI_OFFSET_KIND) start(4)
      integer(kind=MPI_OFFSET_KIND) edge(4)
      integer tabdim_(4),ndim_
#ifdef synopsis
       subroutinetitle='initial_mask_h_from_file'
       subroutinedescription=                                         &
          'Sea-land mask and bathymetry provided by reading an input' &
       //' netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      status=nfmpi_open(par%comm2d,initialgridfile_txt,nf_nowrite,MPI_INFO_NULL,ncid_)
      if(status/=0)stop 'erreur ouverture 4 fichier initialgridfile_txt'

! Lire mask_t:

       varstart(4)=1              ; varcount(4)=1        ! time
       varstart(3)=1              ; varcount(3)=kmax     ! k
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2   ! j
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2   ! i
       start(1:4) = varstart(1:4)
       edge(1:4)  = varcount(1:4)

                   status=nfmpi_inq_varid(ncid_,'mask_t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'tmask',var_id)
      if(status/=0)stop ' stop initial_mask var_id tmask'

      status=nfmpi_inq_vartype(ncid_,var_id,k0)

      if(vert_axis_conv_direc=='up') then
       k1=1 ; k2=kmax ; k3=1
      endif
      if(vert_axis_conv_direc=='dw') then
       k2=1 ; k1=kmax ; k3=-1
      endif

      if(k0==1) then  !----------------------------->
      status=nfmpi_get_vara_int1_all(ncid_,var_id,    & !09-09-17
                                 start(1:4),edge(1:4) &
                  ,mask_t(0:imax+1,0:jmax+1,k1:k2:k3))
      if(status/=0)stop 'initial_mask_h_from_file erreur lecture mask_t'
      else            !----------------------------->
        write(6,*)'k0=',k0
        stop 'initial_mask: cas mask non int pas prevu'
      endif           !----------------------------->

      do j=0,jmax+1
      do i=0,imax+1
       mask_t(i,j,kmax+1)=mask_t(i,j,kmax)
      enddo
      enddo

! Lire dz_t:
                   status=nfmpi_inq_varid(ncid_,'dz_t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'e3t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'e3t_0',var_id)
      if(status/=0)stop ' stop initial_mask_bathy var_id e3t'

      status=nfmpi_inq_var(ncid_,var_id,texte30,k0,ndim_,tabdim_,i0)
      if(status/=0)stop 'Stop nfmpi_inq_var initial_mask_h_from_file'

       if(ndim_==2) then !--------->
        start(2)=1              ; edge(2)=1        ! time
        start(1)=1              ; edge(1)=kmax     ! k
        i1=1 ; i2=1 ; j1=1 ; j2=1
       else              !--------->
        varstart(4)=1              ; varcount(4)=1        ! time
        varstart(3)=1              ; varcount(3)=kmax     ! k
        varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2   ! j
        varstart(1)=1+par%timax(1) ; varcount(1)=imax+2   ! i
        start(1:4)=varstart(1:4)
         edge(1:4)=varcount(1:4)
         i1=0 ; i2=imax+1 ; j1=0 ; j2=jmax+1
       endif             !--------->

      if(vert_axis_conv_direc=='up') then
       k1=1 ; k2=kmax ; k3=1
      endif
      if(vert_axis_conv_direc=='dw') then
       k2=1 ; k1=kmax ; k3=-1
      endif

      if(k0==6) then  !----------------------------->
      status=nfmpi_get_vara_double_all(ncid_,var_id,          &
                                 start(1:ndim_),edge(1:ndim_) &
                                ,dz_t(i1:i2,j1:j2,k1:k2:k3,1))
      if(status/=0)stop 'initial_mask_h_from_file erreur lecture dz_t'
       if(ndim_==2) then !>>>>>>
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         dz_t(i,j,k,1)=dz_t(i1,j1,k,1)
        enddo ; enddo ; enddo
       endif             !>>>>>
      else            !----------------------------->
        write(6,*)'k0=',k0
        stop 'initial_mask: cas e3t non double pas prevu'
      endif           !----------------------------->

! CAS DE LA CELLULE PARTIAL STEP
! Si le tableau e3t depend de x,y,z la cellule partial step est prise en compte
! et on passe (goto 1361)
       if(ndim_>2)goto 1361
! Si le tableau e3t du fichier netcdf ne depend pas de x et y, il faut trouver
! cette cellule dans les autres variables du fichier:

!#ifdef bidon
! Prendre en compte la cellule particuliere du partial step de nemo !19-03-14
      status=nfmpi_inq_varid(ncid_,'e3t_ps',var_id)
      if(status==0) then !pspspspspspspsps>
       status=nfmpi_inq_var(ncid_,var_id,texte30,k0,ndim_,tabdim_,i0)
       if(status/=0)stop 'Stop nfmpi_inq_var initial_mask_h_from_file'
        if(ndim_==3) then !33333333>
         start(3)=1              ; edge(3)=1        ! time
         start(2)=1+par%tjmax(1) ; edge(2)=jmax+2   ! j
         start(1)=1+par%timax(1) ; edge(1)=imax+2   ! i
         status=nfmpi_get_vara_double_all(ncid_,var_id,          &
                                    start(1:ndim_),edge(1:ndim_) &
                                     ,xy_t(0:imax+1,0:jmax+1,1))
! Pour savoir quel niveau vertical est concerne on analyse le masque nemo
! Ou alors on peut regarder la variable mbathy
        do j=0,jmax+1 ; do i=0,imax+1
         do k=2,kmax
          if(mask_t(i,j,k)==1.and.mask_t(i,j,k-1)==0) then
!          write(20+par%rank,*)i,j,dz_t(i,j,k,1),xy_t(i,j,1)
           dz_t(i,j,k,1)=xy_t(i,j,1)
          endif
         enddo
        enddo         ; enddo
        endif             !33333333>

      endif              !pspspspspspspsps>
!#endif

! Lire la variable depw (bathymetry?)
      status=nfmpi_inq_varid(ncid_,'hdepw',var_id)
      if(status==0) then !---------->
       status=nfmpi_inq_var(ncid_,var_id,texte30,k0,ndim_,tabdim_,i0)
        if(ndim_==3) then !33333333>
         start(3)=1              ; edge(3)=1        ! time
         start(2)=1+par%tjmax(1) ; edge(2)=jmax+2   ! j
         start(1)=1+par%timax(1) ; edge(1)=imax+2   ! i
         status=nfmpi_get_vara_double_all(ncid_,var_id,          &
                                    start(1:ndim_),edge(1:ndim_) &
                                     ,xy_t(0:imax+1,0:jmax+1,2))
!                                    ,h_w(0:imax+1,0:jmax+1))
        endif             !33333333>
!       do j=0,jmax+1 ; do i=0,imax+1
!        if(mask_t(i,j,kmaxp1)==1)write(10+par%rank,*)i,j,h_w(i,j),xy_t(i,j,2)
!       enddo         ; enddo

!       if(par%rank==0) then
!        i=6 ; j=70
!        write(6,*)'hihi ',xy_t(i,j,2),h_w(i,j)
!        do k=kmax,1,-1
!         write(6,*)'hoho ',k,mask_t(i,j,k),dz_t(i,j,k,1)
!        enddo
!        stop 'haha'
!       endif

#ifdef bidon
        do j=0,jmax+1
        do i=0,imax+1

         k=1
         sum1=dz_t(i,j,k,1)*mask_t(i,j,k)
         if(mask_t(i,j,k)==1)k0=1
         do k=2,kmax
          sum1=sum1+dz_t(i,j,k,1)*mask_t(i,j,k)
          if(mask_t(i,j,k)-mask_t(i,j,k-1)==1)k0=k
         enddo
         x1=sum1-h_w(i,j)
!        if(par%rank==0.and.i==6.and.j==70)write(6,*)'x1=',x1,kmax-k0+1
!        if(par%rank==0.and.i==17.and.j==85)write(6,*)'x1=',x1,kmax-k0+1
         if(mask_t(i,j,k)==1)write(30+par%rank,*)i,j,x1
        enddo
        enddo
        if(par%rank==0)stop 'bm'
#endif

      endif              !---------->

 1361 continue

      do j=0,jmax+1
      do i=0,imax+1
       h_w(i,j)=0.
      enddo
      enddo
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       h_w(i,j)=h_w(i,j)+mask_t(i,j,k)*dz_t(i,j,k,1)
      enddo
      enddo
      enddo
      do j=0,jmax+1
      do i=0,imax+1
       h_w(i,j)=max(h_w(i,j),1.d00)
      enddo
      enddo

!       if(par%rank==0) then
!        i=6 ; j=70
!        write(6,*)'hoho ',xy_t(i,j,2),h_w(i,j)
!        stop 'haha'
!       endif

      status=nfmpi_close(ncid_)

!     do j=0,jmax+1
!     do i=0,imax+1
!     if(mask_t(i,j,kmax)==1)write(10+par%rank,*)h_w(i,j),xy_t(i,j,2)
!     enddo
!     enddo
!     stop 'cocote'

      end subroutine initial_mask_h_from_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#ifdef bidon
      subroutine interp_emodnet
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      use module_global
      implicit none
      double precision resol_,h_,weight_,resolution_deg_
#ifdef synopsis
       subroutinetitle='interp_emodnet'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      weight_=1.

! Reset des compteurs
      do j=1,jmax
      do i=1,imax
       xy_t(i,j,1)=0.
       xy_t(i,j,2)=0.
       xy_t(i,j,3)=h_w(i,j)
      enddo
      enddo

!     open(unit=3,file='../../../data/bathymetry/japon_gebco.txt',recl=10000)
      open(unit=3,file= &
      '/home/marp/data/bathy_emodnet/bathy_emodnet.ascii')

      read(3,*)resol_
      resol_=resol_*deg2rad*rayonterre
      write(6,*)'resolution en metres=',resol_

   99 read(3,*,end=100)x1,x2,h_

      call latlontoij(x1*deg2rad,x2*deg2rad,'loc')
      i=nint(deci)          ; j=nint(decj)


      if(i<1   ) goto 99
      if(i>imax) goto 99
      if(j<1   ) goto 99
      if(j>jmax) goto 99

      i2=min(max(i,1),imax) ; j2=min(max(j,1),jmax)
! k = rapport de resolution "donnee/modele"
      rap=resol_/sqrt(dxdy_t(i2,j2))
      k=nint(rap)
! Si k=0 la donnee ne se projete que sur une seule cellule (la plus proche)
! Ce cas ne se produit que si la grille des donnees est plus fine que celle du modele

!     write(67,*)k,i,j,sqrt( dxdy_t(i,j)),rayonterre*resol_*deg2rad*cos(lat_t(i2,j2))


!     const1=(1./resol_)**2
      do j1=j-k,j+k
       if(j1>=1.and.j1<=jmax) then !jjjjjjjjjj>
         do i1=i-k,i+k
          if(i1>=1.and.i1<=imax) then !iiiiiiiii>
            xy_t(i1,j1,1)=xy_t(i1,j1,1)+weight_
            xy_t(i1,j1,2)=xy_t(i1,j1,2)+weight_*h_
          endif                       !iiiiiiiii>
         enddo
       endif                      !jjjjjjjjjjj>
      enddo

      goto 99
  100 close(3)

      mask_t=0
      h_w=0.1
      do j=1,jmax
      do i=1,imax
!     if(mask_t(i,j,kmax+1)==1) then !1111111>
      if(xy_t(i,j,1)==0.)then        !---------->
!           write(66,*)i,j,xy_t(i,j,1),weight_
!          stop 'Trou bathy'
      else                           !---------->
       x0=xy_t(i,j,2)/xy_t(i,j,1)
       if(x0>h_inf)mask_t(i,j,kmax+1)=1
       h_w(i,j)=max(x0,h_inf)
      endif                          !---------->
!     endif                        !1111111>
      enddo
      enddo

! Archiver le resultat dans un fichier bathycote_in.dat
      write(*,*)'Archive dans tmp/bathycote_in.dat'
      open(unit=3,file=trim(tmpdirname)//'bathycote_in.dat',recl=10000)
        do i=1,imax
        write(3,'(100000i1)')(mask_t(i,j,kmax+1),j=1,jmax) !16-05-15
        enddo
        do i=1,imax
        write(3,*)(h_w(i,j),j=1,jmax)
        enddo
      close(3)

! Decommenter pour faire graphique des differences
!     do j=1,jmax
!     do i=1,imax
!      h_w(i,j)=h_w(i,j)-xy_t(i,j,3)
!     enddo
!     enddo

!     stop ' INTERPOLATION GEBCO TERMINEE'

      end subroutine interp_emodnet
!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine load_mask_h_from_a_netcdffile !19-04-15
      use module_principal ; use module_parallele ; use module_global
      implicit none
      integer ncid_,istr_,istp_,jstr_,jstp_,var_dims_
!     integer*1,dimension(:,:),allocatable :: glob_mask_i1
!     real*8   ,dimension(:,:),allocatable :: glob_mask_r8
      real*4   ,dimension(:,:),allocatable :: glob_h_r4
      include 'netcdf.inc'

      status=nf_open(trim(texte250),nf_nowrite,ncid_)

      if(status==0)then    !sss>
       if(par%rank==0)write(6,'(a,a,a)')'nf_open ',trim(texte250),' OK!'
      else                 !sss>
       write(6,'(a,a)')trim(texte250),' not found' !06-07-14
       stop ' Err 1666 nf_open LSM+H netcdf file'
      endif                !sss>

! definir/verifier varstart:
                    status=nf_inq_dimid(ncid_,'ni_t',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'nx',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'nx_t',k0)   !12-03-20
       if(status/=0)status=nf_inq_dimid(ncid_,'xi_rho',k0) !10-09-16
       if(status/=0)stop 'Err1890'
       status=nf_inq_dimlen(ncid_,k0,i1) ; if(status/=0)stop 'Err1891'
       varstart(1)=1 ; varcount(1)=i1

                    status=nf_inq_dimid(ncid_,'nj_t',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'ny',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'ny_t',k0)    !12-03-20
       if(status/=0)status=nf_inq_dimid(ncid_,'eta_rho',k0) !10-09-16
       if(status/=0)stop 'Err1892'
       status=nf_inq_dimlen(ncid_,k0,j1)    ;if(status/=0)stop 'Err1893'
       varstart(2)=1 ; varcount(2)=j1

                    status=nf_inq_dimid(ncid_,'nk_t',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'nz',k0)
       if(status==0)then !>
         status=nf_inq_dimlen(ncid_,k0,k1) ;if(status/=0)stop 'Err1895'
         varstart(3)=k1 ; varcount(3)=1
       endif             !>

! Debug dims:
       istr_=-9999
       if(i1==iglb  .and.j1==jglb)   then !13-06-15
        istr_=1 ; istp_=iglb ; jstr_=1 ; jstp_=jglb
       endif
       if(i1==iglb+2.and.j1==jglb+2) then !13-06-15
        istr_=0 ; istp_=iglb+1 ; jstr_=0 ; jstp_=jglb+1
       endif
       if(istr_==-9999) then
         write(6,*)'iglb jglb i1,j1=',iglb,jglb,i1,j1
         stop 'Err dim2 netcdf LSM H File'
       endif

! Lire glob_mask
                    status=nf_inq_varid(ncid_,'landmask',var_id)  
       if(status/=0)status=nf_inq_varid(ncid_,'landmask_t',var_id)  
       if(status/=0)status=nf_inq_varid(ncid_,'mask2d_t',var_id)  
       if(status/=0)status=nf_inq_varid(ncid_,'mask_rho',var_id)  !10-09-16
       if(status/=0)status=nf_inq_varid(ncid_,'mask_t',var_id)  !10-09-16
       if(status/=0)stop 'Err 1667 nf_inq_varid mask_t'
       status=nf_inq_varndims(ncid_,var_id,var_dims_) ; if(status/=0)stop 'Err 1937 nf_inq_varndims'
       if(var_dims_/=2.and.var_dims_/=3)stop 'Err 1938 var_dims_'

       status=nf_inq_vartype(ncid_,var_id,k0)      
       if(k0/=nf_short.and.k0/=nf_byte.and.k0/=nf_int  &
                                      .and.k0/=nf_double) then !10-09-16
         write(6,*)'k0=',k0
         write(6,*)'nf_int=',nf_int
         write(6,*)'nf_short=',nf_short
         write(6,*)'nf_byte=' ,nf_byte
         write(6,*)'nf_double=' ,nf_double
         stop 'Err 1668 type mask non prevu'
       endif
!      if(k0==nf_short)then !>>>>> 
       if(k0==nf_short.or.k0==nf_byte)then !>>>>> 

       status=nf_get_vara_int1(ncid_,var_id,varstart(1:var_dims_)   &
                                           ,varcount(1:var_dims_)   &
                             ,glob_mask(istr_:istp_,jstr_:jstp_))  !13-06-15
       if(status/=0)stop 'Err 2177 nf_get_vara_int1 glob_mask'

       else                 !>>>>>

       if(k0==nf_double)stop 'Err 2126 k0==nf_double'
!      if(k0==nf_byte)  stop 'Err 2126 k0==nf_byte'

       if(k0==nf_int) then !integer> !08-05-18

         varstart(1)=1+par%timax(1) 
         varstart(2)=1+par%tjmax(1) 
         if(istp_==iglb.and.jstp_==jglb) then     !>>>
          varcount(2)=jmax ; varcount(1)=imax
          status=nf_get_vara_int(ncid_,var_id,varstart(1:2),varcount(1:2),anyv3dint(1:imax,1:jmax,1))
         endif                                    !>>>
         if(istp_==iglb+2.and.jstp_==jglb+2) then !>>>
          varcount(2)=jmax+2 ; varcount(1)=imax+2
          status=nf_get_vara_int(ncid_,var_id,varstart(1:2),varcount(1:2),anyv3dint(0:imax+1,0:jmax+1,1))
         endif                                    !>>>
         if(status/=0)stop 'Err 2142 nf_get_vara_int anyv3dint'
         do j=1,jmax   ; do i=1,imax  
          glob_mask(i+par%timax(1),j+par%tjmax(1))=anyv3dint(i,j,1)
         enddo ; enddo
#ifdef parallele
         ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax                 !23-09-09
         ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax  
         lb2=lbound(glob_mask) ; ub2=ubound(glob_mask)
         call  par_gatherall_2d(glob_mask,lb2,ub2,ind,par%nbdom)
#endif
       endif               !integer>

       endif                !>>>>>
       if(status/=0)stop ' Err1896 nf_get_vara_int glob_mask'

! Attention que les tools de Florent ont une convention -1/1 pour distinguer le masque terre/mer
! Retablir par consequent la convention standard 0/1:
        glob_mask=max(glob_mask,0)

! Lire glob_h
       varstart(1)=1 ; varstart(2)=1 
       if(istp_==iglb.and.jstp_==jglb) then     !>>>
         varcount(2)=jglb ; varcount(1)=iglb
       else                                     !>>>
         if(   istp_==iglb+1.and.   jstp_==jglb+1) then !ooo>
          varcount(2)=jglb+2 ; varcount(1)=iglb+2
         else                                           !ooo>
          write(6,*)'istp_,iglb+1,jstp_,jglb+1',istp_,iglb+1,jstp_,jglb+1
          stop 'Err 2166 in the definition of varcount'
         endif                                          !ooo>
       endif                                    !>>>
                    status=nf_inq_varid(ncid_,'bathymetry_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'h_w',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'h',var_id) !10-09-16
       if(status/=0)stop 'Err1896 nf_inq_varid h_w'
       status=nf_inq_vartype(ncid_,var_id,k0)

       if(k0==nf_real)then !>>>>>
        status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                            ,varcount(1:2)   &
                            ,glob_h(istr_:istp_,jstr_:jstp_))
       else                  !>>>>>
        if(k0==nf_double) then !double> !08-05-18

         varstart(1)=1+par%timax(1) 
         varstart(2)=1+par%tjmax(1) 
         ksecu=0
         if(istp_==iglb.and.jstp_==jglb) then     !>>>
          varcount(2)=jmax ; varcount(1)=imax ; ksecu=1
          status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),anyv3d(1:imax,1:jmax,1,1))
         endif                                    !>>>
         if(istp_==iglb+1.and.jstp_==jglb+1) then !>>>
          varcount(2)=jmax+2 ; varcount(1)=imax+2 ; ksecu=1
          status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),anyv3d(0:imax+1,0:jmax+1,1,1))
         endif                                    !>>>
         if(ksecu==0)stop 'Err 2196 in the definition of varcount'
         if(status/=0)stop 'Err 2142 nf_get_vara_double anyv3d'
         do j=1,jmax   ; do i=1,imax  
          glob_h(i+par%timax(1),j+par%tjmax(1))=anyv3d(i,j,1,1)
         enddo ; enddo

#ifdef parallele
         ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax                 !23-09-09
         ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax  
         lb2=lbound(glob_h) ; ub2=ubound(glob_h)
         call  par_gatherall_2d(glob_h,lb2,ub2,ind,par%nbdom)
#endif
        else                   !double>
         stop 'Err 2170 type bathy not available' !08-05-18
        endif                  !double>
       endif                 !>>>>>

       if(status/=0)stop 'Err1898 nf_get_vara_double glob_h'

      status=nf_close(ncid_) ; if(status/=0)stop 'Err1899 nf_close'


      if(par%rank==0) then !------------>
       open(unit=3,file=trim(tmpdirname)//'bathycote_in.dat',recl=10000) !13-06-15
        do i=1,iglb
        write(3,'(100000i1)')(glob_mask(i,j),j=1,jglb)
        enddo
        do i=1,iglb
        write(3,'(100000(1x,f10.4))')(glob_h(i,j),j=1,jglb)
        enddo
       close(3)
       write(6,'(a)')'file bathycote_in.dat created in tmp directory'
      endif                !------------>

      end subroutine load_mask_h_from_a_netcdffile

!.................................................................

      subroutine remove_secondary_bassins !08-05-18
      use module_principal
      use module_parallele !#MPI
      use module_global
      implicit none
      integer biggest_bassin_,loop_

! version rendue compatibe avec glob_mask(kind=1) le 08-05-18 
! l'ancienne version ne permettait pas de depasser un numero d'ile>127

      biggest_bassin_=0 ! 05-03-12
      sum0=0.

      loop_=1
      do j1=1,jglb
      do i1=1,iglb

      if(glob_mask(i1,j1)==1)then !11111111111111111111>

        loop_=loop_+1
        glob_mask(i1,j1)=2 ! flag d'un nouveau bassin A resencer

        sum1=0.
        sum2=0.
 1234   continue

        do j=1,jglb ; do i=1,iglb
         if(glob_mask(i,j)==2)then !---->
              if(glob_mask(min(i+1,iglb),j            )/=0)glob_mask(min(i+1,iglb),j            )=2
              if(glob_mask(max(i-1,1   ),j            )/=0)glob_mask(max(i-1,1   ),j            )=2
              if(glob_mask(i            ,min(j+1,jglb))/=0)glob_mask(i            ,min(j+1,jglb))=2
              if(glob_mask(i            ,max(j-1,1   ))/=0)glob_mask(i            ,max(j-1,1   ))=2
         endif                     !---->
        enddo       ; enddo

        sum2=sum1
        sum1=0.
        do j=1,jglb ; do i=1,iglb
         if(glob_mask(i,j)==2)sum1=sum1+1.
        enddo       ; enddo
        if(sum1/=sum2)goto 1234

        if(sum1>sum0) then !m°v°m>

          sum0=sum1 ; biggest_bassin_=loop_
! Si le bassin est le plus grand bassin de tous jusqu'A present identifiEs
! on lui applique la valeur 3 mais avant on deflague le bassin
! qui jusqu'A present detenait le record:
         do j=1,jglb ; do i=1,iglb
          if(glob_mask(i,j)==3)glob_mask(i,j)=0 ! deflaguer l'ancien record
          if(glob_mask(i,j)==2)glob_mask(i,j)=3 ! flaguer le nouveau record
         enddo       ; enddo

        else               !m°v°m> !04-01-23

         do j=1,jglb ; do i=1,iglb
          if(glob_mask(i,j)==2)glob_mask(i,j)=0 !04-01-23
         enddo       ; enddo
         
        endif              !m°v°m>


      if(par%rank==0)write(6,*)'Basin, size, main basin' &
       ,loop_,nint(sum1),biggest_bassin_
     

      endif                     !11111111111111111111>

      enddo
      enddo

      do j=1,jglb ; do i=1,iglb
       if(glob_mask(i,j)==3) then
          glob_mask(i,j)=1
       else
          glob_mask(i,j)=0
       endif
      enddo ; enddo

      if(par%rank==0) then !000000000000000>
       open(unit=3,file=trim(tmpdirname)//'mask_bassin1.out')
        do i=1,iglb
         write(3,'(1000i1)')(glob_mask(i,j),j=1,jglb)
        enddo
       close(3)
      endif                !000000000000000>
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

      end subroutine remove_secondary_bassins

!........................................................................

      subroutine initial_mask_and_bathy_corrector !01-06-17
      use module_principal ; use module_parallele
      implicit none
      real*4 :: delta_h_=0. , dist_h_
      double precision dist_ 

! https://docs.google.com/document/d/16A8Pb0JjOLYM1cI3CUiwJgDqkAfJhJbr7fkasAPnBoc/edit
! Corriger des biais introduits par le lissage de la bathymetrie avec
! des corrections de grande longueur d'onde


      ksecu=0

      if(index(texte250,'.ijDeltaH')/=0.or. & 
         index(texte250,'.lonlatDeltaH')/=0) then ! m[°v°]m >

      ksecu=1
! Le fichier de correction contient, pour chaque ligne: i,j,delta_h, D0 echelle de distance en indices
! La correction est dh=delta_h*exp( -(distance/d0)**2 )

#ifdef bidon
       ub2=ubound(h_w) ; lb2=lbound(h_w)
       if(par%rank==0)write(6,'(a,a)')'opening ',trim(texte250)
       open(unit=3,file=trim(texte250))
 2265   read(3,*,end=2263)i1,j1,delta_h_,i10
        i1=i1-par%timax(1) ; j1=j1-par%tjmax(1)
        x0=-(1./real(i10))**2
        do j=lb2(2),ub2(2) ; do i=lb2(1),ub2(1)
         h_w(i,j)=h_w(i,j)+delta_h_*exp(x0*( (i-i1)**2 + (j-j1)**2 ))
!        h_w(i,j)=1000.   +delta_h_*exp(x0*( (i-i1)**2 + (j-j1)**2 )) ! Pour test
        enddo              ; enddo
        goto 2265
 2263  close(3)
#endif
! Nouvel algo: si 2 corrections distinctes ont une zone de chevauchement !30-10-19
! on n'additionne pas les 2 corrections si amplification. Addition seulement si compensation.
! Pour cela on passe passe par xy_t(:,:,0) la prise en compte de toutes les delta
! avec la contrainte ci-dessus puis on ajoute A h_w
       xy_t(:,:,0)=0. ! reset de la correction
       ub2=ubound(h_w)  ; lb2=lbound(h_w)
       ub3=ubound(xy_t) ; lb3=lbound(xy_t)
       if(lb2(1)<lb3(1).or.lb2(2)<lb3(1).or. &
          ub2(1)>ub3(1).or.ub2(2)>ub3(2))    &
       stop 'Err 2300 bound in initial_mask_and_bathy_corrector'
       if(par%rank==0)write(6,'(a,a)')'opening ',trim(texte250)

       open(unit=3,file=trim(texte250))
 2265  continue

       if(index(texte250,'.ijDeltaH')/=0) then !-cas-indices->
        read(3,*,end=2263)i1,j1,delta_h_,i10

        x0=-(1./real(i10))**2
        do j=lb2(2),ub2(2) ; do i=lb2(1),ub2(1)
!22 454 10 3
         x10=delta_h_*exp(x0*( (i+par%timax(1)-i1)**2 + (j+par%tjmax(1)-j1)**2 ))

! Si la perturbation est negligeable ne rien faire !09-03-23
         if(abs(x10)>small1) then !pmx> !09-03-23

         if(xy_t(i,j,0)*x10<0.) then !oooo> 
! Si les perturbations se compensent on les additionne:
            xy_t(i,j,0)=x10+xy_t(i,j,0)
         else                        !oooo>
! Si les perturbations s'amplifient, on prend la plus grande des 2
          if(x10>0.) then !>>>>
            ! sous-cas perturbation positive
            xy_t(i,j,0)=max(x10,xy_t(i,j,0))
          endif           !>>>>
          if(x10<0.) then !>>>> !09-03-23 (le cas x10=0 ne doit pas jouer)
            ! sous-cas perturbation negative
            xy_t(i,j,0)=min(x10,xy_t(i,j,0))
          endif           !>>>>   
         endif                       !oooo>

         endif                    !pmx> !09-03-23
        enddo              ; enddo
       endif                                   !-cas-indices->

       if(index(texte250,'.lonlatDeltaH')/=0) then !-cas-lonlat->

        read(3,*,end=2263)longi1,latit1,delta_h_,dist_h_
        longi1=longi1*deg2rad
        latit1=latit1*deg2rad

!       call latlontoij(longi1,latit1,'glb')
!       i1=nint(deci)
!       j1=nint(decj)
!       x0=-(1./real(i10))**2
!       x1=1./dist_h_**2

        do j=lb2(2),ub2(2) ; do i=lb2(1),ub2(1)

        call lonlat2distance(longi1,latit1,lon_t(i,j),lat_t(i,j),dist_) !31-01-23

        x10=delta_h_*exp(-abs(dist_/dist_h_))

! Si la perturbation est negligeable ne rien faire !09-03-23
         if(abs(x10)>small1) then !pmx> !09-03-23

!       x0=-x1*dxdy_t(i,j) ! x0=-(dxdy/dist**2)
!       x10=delta_h_*exp(x0*( (i+par%timax(1)-i1)**2 + (j+par%tjmax(1)-j1)**2 ))

         if(xy_t(i,j,0)*x10<0.) then !oooo>
! Si les perturbations se compensent on les additionne:
            xy_t(i,j,0)=x10+xy_t(i,j,0)
         else                        !oooo>
! Si les perturbations s'amplifient, on prend la plus grande des 2
          if(x10>0.) then !>>>>
            ! sous-cas perturbation positive
            xy_t(i,j,0)=max(x10,xy_t(i,j,0))
          endif           !>>>>
          if(x10<0.) then !>>>> !09-03-23 (le cas x10=0 ne doit pas jouer)
            ! sous-cas perturbation negative
            xy_t(i,j,0)=min(x10,xy_t(i,j,0))
          endif           !>>>>   

         endif                       !oooo>

         endif                    !pmx> !09-03-23
        enddo              ; enddo
       endif                                       !-cas-lonlat->

       goto 2265

 2263  close(3)
! Une fois toutes les corrections chargees, les ajouter A la bathy
        do j=lb2(2),ub2(2) ; do i=lb2(1),ub2(1)
         h_w(i,j)=h_w(i,j)+xy_t(i,j,0)
        enddo              ; enddo

      endif                                   ! m[°v°]m >


      if(ksecu==0) then !pmxpmx>
       stop 'bathy corrector file name extension not recognised'
       call mpi_finalize
      endif             !pmxpmx>


      end subroutine initial_mask_and_bathy_corrector

!.................................................................

      subroutine update_bathy !17-05-22
      use module_principal ; use module_parallele
      implicit none

! Modification de la bathy en cours de simulation et ajustement de la SSH 
! pour conservation de la temperature et de la sainite

! Charger dans le tableau anyvar2d le delta de bathymetrie
      if(iteration3d>0) then
       do j=2,jmax-1
       do i=2,imax-1
! Pour faire un test academique:
        xy_t(i,j,0)=0.001*exp(-( (real(i+par%timax(1))-0.5*real(iglb))**2       &
                                +(real(j+par%tjmax(1))-0.5*real(jglb))**2)/100.)
       enddo
       enddo
      else
       xy_t(:,:,0)=0.
      endif

! Conditions aux limites laterales:
      if(obcstatus(ieq1)==1)   xy_t(0:1        ,:,0)=0.
      if(obcstatus(ieqimax)==1)xy_t(imax:imax+1,:,0)=0.

      if(obcstatus(jeq1)==1)   xy_t(:,0:1        ,0)=0.
      if(obcstatus(jeqjmax)==1)xy_t(:,jmax:jmax+1,0)=0.

! Assurer la continuite mpi
      call obc_ext_xy_t('za',0) 


! Mettre A jour la bathy:
      do j=0,jmax+1
      do i=0,imax+1

       h_w(i,j)=h_w(i,j)+xy_t(i,j,0)

      enddo
      enddo
      call hz_to_hxyr

! Ajuster la SSH pour maintenir l'epaisseur de la colonne d'eau inchangee: (noter le signe -)
! Note: Corriger le temps 1 revient A corriger le temps 2, corriger ssh_w entrainera la correction 
!       de ssh_int_w, hz_w, dz_t etc etc.....
      do j=0,jmax+1
      do i=0,imax+1

       ssh_w(i,j,1)=ssh_w(i,j,1)-xy_t(i,j,0)

      enddo
      enddo

!     call graph_out

      end subroutine update_bathy
