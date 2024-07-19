      module module_wave
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 - last update: 19-02-23
!______________________________________________________________________
!------------------------------------------------------------------------------!
! version   date    description                                                !
! 2008.14 24-06-09  mise en service                                            !
!                   note: pour s'adapter aux sorties de ww3, l'amplitude       !
!                   est remplacee par la hauteur significative                 !
!         25-06-09  attention on n'interpole pas un angle mais ses composantes !
!                   trigonometriques                                           !
!         26-06-09  cette remarque vaut egalement pour l'interpollation dans   !
!                   le temps                                                   !
!         27-06-09  - moyenne sur points r: prise en compte du masque          !
!                   - ne pas multiplier x10 par mask_u ou mask_v sinon     !
!                   l'integrale des termes 3d n'est plus coherente avec        !
!                   l'approche 2d (en plus ca engendre des discontinuités      !
!                   sur le champ 2d)                                           !
!                   - puisque la boucle verticale commence depuis 1 (et non    !
!                   pas depuis kmin, il est important que sigma des niveaux    !
!                   souterrains soient masqués, d'où la multiplication de      !
!                   sigma_z par mask_t                                       !
!         28-06-09  nous marquons une pause dans le developpement. un "return" !
!                   sort de la routine apres interpollation des champs de ww3  !
!                   de sorte que les parametres requis par le modele de        !
!                   transport sedimentaire a ses infos mais le calcul des      !
!                   tensions de radiations, (pas encore validé), n'est pas     !
!                   activé                                                     !
!         29-06-09  possibilité (à condition de decommenter les bonnes lignes) !
!                   d'interpolation bi-polynome ordre 3 à la place de bi-      !
!                   lineaire                                                   !
! 2008.15 26-08-09  - ajout cas swan                                           !
!                   - debug boucle sur i                                       !
!         28-08-09  - assistance aux utilisateurs quand la date en cours n'est !
!                   pas compatible avec le contenu des fichiers de vagues      !
! 2010.7  27-02-10  passage à la methode de Bennis et ardhuin 2010             !
!         04-03-10  sshstokes
! 2010.8  11-03-10  ecrire dans tmp
!         18-03-10  appliquer le wetmask sur les vitesses de stokes
!         21-05-10  - velbarstokes calculé dans z_averaged pour conservation
!                   - addition de 2 stress -> allocation de 2 tableaux
! 2010.12 19-09-10  ajout cas fichiers ww3 multifrequentiel
!         25-09-10  obc sur surfstress (parallelisation)
! 2010.13 20-10-10  - mise à jour nom de variable ww3
!                   - vus possibles aleas de la grille ww3 la compatibilité
!                   de la grille du modele oceanique avec celle de la grille ww3
!                   est verifiee pour chaque grille ww3
!         03-11-10  des arguments passés dans date_to_kount
! 2010.14 26-11-10  ajout de commentaires explicatifs pour la loi d'interpolation
!                   temporelle du courant de stokes
!         28-11-10  suite du point precedent: correction bug sur k2
! 2010.18 02-02-11  nouvel algo (+rapide+blindage) de l'inversion de la
!                   relation de dispersion
! 2010.20 15-04-11  renouvellement des fichiers sur la base d'un temps calculé
!                   en secondes
! 2010.22 06-05-11  ajout foc et hsw
!         07-05-11  Blindage calculs
! 2010.23 18-05-11  preciser dans notebook_wave le type de schema de conditions
!                   aux limites ouvertes lateralles
! 2010.25 23-02-12  waveforcing devient module_wave
!         29-01-13  MAJ routine wave_read_interp_netcdf et on enleve les potentiels sources
!                   d'instabilites
! S26     14-02-13  debug lignes d'interpolation hsw et J dans filiere ww3
!                   blindage calcul rotation symphonie
! S.26    11-04-13  Filiere IOWAGA
!         14-05-13  deplacement de boucles
!         07-09-13  filiere iowaga s26
!         25-05-14  mises a jour diverses, notamment sur la fabrication des listes
!                   binaires
!         20-10-14  tmpfile ouvert avec open status='replace'
!         25-01-15  ajout cas ww3_on_sgrid
!         03-02-15  Pour eviter division par zero dans momentum_equation seuillage
!                   sur hsw
!         17-02-15  Seuillage sur t_wave_t(i,j,2)
!         18-02-15  bancs decouvrants: hauteur des vagues bornee par h+ssh
!         27-02-15  Detection de listes de fichiers desordonnees
!         03-03-15  Detection de listes de fichiers desordonnees SUITE
!         06-03-15  Calcul Bernoulli Head zone decouverte
!         08-03-15  suppresion test redondant source de bug
!         21-03-15  argument 'za' dans call obc_int_anyvar2d
!         02-09-16  nouvel algo pourla lecture de txt_units
!         11-09-16  rendre compatible iwve=0 et option de compilation stokes
!         13-09-16  rendre possible un restart avec activation des vagues
!         23-09-16  un seuil plus robuste pour hsw_wave_t
!         24-09-16  Ajout d'un cas supplementaire (lon,lat real, structurE, grille differente)
!         30-09-16  version couplage ww3 par OASIS
!         30-10-16  Borner two
!         14-04-17  Seuil sur hw_z à 1m pour les bancs découvrants
!         12-09-18  3eme indice uss passe A 2 pour compatibilite avec cas oasis
! v250    25-03-19  ajout subroutines pour cas SETE 
!         27-03-19  x0=sqrt(2.*ww3efth(dirw,freq,1,1)*ww3deltafreq(freq)*ww3deltadir) 
! v259    01-10-19  Mises A jour interface WW3 Symphonie Phase Resolue
! v261    16-10-19  Appliquer une condition intertidale sur la SSH reconstituee:
! v267    17-11-19  do k=1,kmin_w(i,j)-1 
!                   cas iwve==2 (McWilliams99)
! v268    21-11-19  - liste unique dans notebook_wave (puisque WW3 regroupe toutes les variables dans un meme fichier)
!                   - lire seulement periode et courant de stokes de surface dans le cas iwve=2 (param McWilliams1999)
! v269    02-12-19  velwave_u et v de nouveau sur point u et v
! v280    02-05-20  subroutine wave_hamilton_ebersole !02-05-20
! v285    05-06-20  suppression zeroleveldepth_u et v
! v290    23-10-20  Si on n'appelle pas module_airseaflux remise A zero de wstress_u et wstress_v !23-10-20
!         30-10-20  etendre les dimensions de velstokes pour drifter
!         04-11-20  ajout subroutine wave_init_ww3plugs !04-11-20
! v292    20-11-20  iwve==3
! v311    11-11-21  adaptation aux fichiers de la simu WW3 SEA
! v350    25-07-22  ww3_iowaga remplacE par ww3_native_grid (pour plus de clartE)
! v366    14-02-23  champs copernicus
! v367    19-02-23  copernicus suite:  ww3varname(2)='VHM0' 
!------------------------------------------------------------------------------!
      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'

      character*60 txt_units
      character*10 ww3varname(5),varname_selected
      character*5 ww3flag1
      real*4 :: ww3mskval,ww3filval
      double precision,dimension(:),allocatable ::   &
        dist_4cadrans
      double precision,dimension(:,:),allocatable :: &
        ww3cg_lon                                    &
       ,ww3cg_lat
      real*4,dimension(:,:),allocatable :: &
        ww3_2dr4_lon                                 &
       ,ww3_2dr4_lat
      real*4,dimension(:,:,:),allocatable ::         &
       anyvwave
      real*4,dimension(:),allocatable ::             &
        ww3_lat                                      &
       ,ww3_lon                                      &
       ,ww3_res
      integer,dimension(:,:),allocatable ::          &
        ww3_short
      real*4,dimension(:,:),allocatable ::          &
        ww3_var
!     integer,dimension(:,:,:),allocatable ::        &
!       ww3_var2
      real,dimension(:),allocatable :: ar_send_sud,ar_recv_nord
      real,dimension(:),allocatable :: ar_send_est,ar_recv_ouest
      integer ww3_element                 &
             ,var_dims,var_type,tabdim(5) ! Que represente tabdim?
                                          ! Soit 1,2,3, l'ordre dans lequel
                                          ! sont listees les dimensions du fichiers
                                          ! netcdf. Tabdim indique (via le numero
                                          ! d'ordre) les dimensions dont dépend la
                                          ! variable
      integer :: reclen=270      &
                ,ww3_flagcumul   & ! ww3 variables are 0=instantaneous 1=time-averaged
                ,flag_ww3_dpt=0  &
                ,wavemodulo=10   &
!               ,wavemodulo=1    &
!               ,timewavenext=2  &
!               ,timewaveprev=0  &
                ,timewavenext=4  &
                ,timewaveprev=3  &
                ,timewaveprev2=0  &
                ,timewaveprev3=-1    

      integer(kind=1) :: flag_allocate_=0      &
                        ,flag_wavelength       & !11-11-21
                        ,ww3_on_sgrid_case     & !11-11-21
                        ,ww3_on_sgrid_casez0=0 & ! grille 1:iglb,   1:jglb   !11-11-21
                        ,ww3_on_sgrid_casez1=1   ! grille 0:iglb+1 ,0:jglb+1 !11-11-21

!     double precision , dimension(:)  , allocatable :: ww3puls,ww3dir
      real :: ww3deltadir
      real , dimension(:)  , allocatable :: ww3puls,ww3dir,ww3deltafreq
      real , dimension(:,:), allocatable :: ww3kvector
      real , dimension(:,:,:), allocatable :: ww3kvec_oj
      real , dimension(:,:,:,:), allocatable :: ww3phase
      real , dimension(:,:,:,:), allocatable :: ww3efth
!     real , dimension(:,:) , allocatable :: cumul1_mpi,cumul2_mpi

      integer , dimension(:,:) , allocatable :: ww3plugs               
      integer(kind=1) :: flag_ww3plugs=1
      integer :: ww3misspointmax=0

contains

!.................................................

      subroutine wave_driver
      implicit none
#ifdef synopsis
       subroutinetitle='wave_driver'
       subroutinedescription= &
       'Driver of the subroutines computing the wave-current effect'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


#ifdef key_oasis_symp_ww3
! Case symphonie-ww3-oasis
      call wave_driver_oasis_symp_ww3
#else

      if(initial_main_status==0) then !oooo>
! Initial state:
!      call wave_ww3_spec_sete_allocate

       call wave_initial 
       call wave_driver_s26

      else                            !oooo>
! Iterative phase:
       call wave_driver_s26

      endif                           !oooo>

#endif

      end subroutine wave_driver

!.................................................

      subroutine wave_driver_s26
      implicit none
#ifdef synopsis
       subroutinetitle='wave_driver_s26'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      ! TODO-CHECK THOM
      if (iwve==0) return
      if(initial_main_status==0) then !>>>>>>>>> !04-11-20

! Initial state - Interpolation of the IOWAGA data base:

       call wave_read_interp_netcdf_s26(0)

       call wave_linear_in_time_s26('initial1')
       wavefile_nextrec=wavefile_nextrec+1


       wavefile_prvtime=wavefile_nextime
       call wave_read_interp_netcdf_s26(2)

       call wave_linear_in_time_s26('initial2')
       wavefile_nextrec=wavefile_nextrec+1


      else                    !>>>>>>>>>

! Iterative phase - Interpolation of the IOWAGA data base provided that elapsedtime_now > wavefile_nextime
! et aussi iteration/=0 puisque desormais on fait la procedure d'initialisation (plus haut) depuis initial_main
! et ce afin que le boucheur de trou soit toujours fait avec le meme fichier ww3 (puisque le LSM de WW3 peut
! varier dans le temps). Cette maniere de faire consolide la restartabilite du code. A propos du LSM variable
! dans le temps, le modele prevoit le cas ou certains trous ne sont plus bouchEs quand LSM change.
       if(elapsedtime_now>=wavefile_nextime.and.iteration3d/=0) then !>>>>>>> !04-11-20
!      if(elapsedtime_now>=wavefile_nextime) then  !>>>>>>>>>

        call wave_read_interp_netcdf_s26(2)
        wavefile_nextrec=wavefile_nextrec+1

        if(wavefile_nextime<=elapsedtime_now)   then ! check the time of the next ww3 file
         stop 'error on time of the next iowaga file'
        endif                                        ! check the time of the next ww3 file

       endif                                       !>>>>>>>>>

      endif                   !>>>>>>>>>

! Wave parameters at intermediate time steps : linear in time interpolation:
      call wave_linear_in_time_s26('standard')


! depth-averaged Stokes current:
      call wave_zaveraged_velstokes

      if(iwve==1) then !m°v°m> !17-11-19

! Add the wave dissipation by bottom friction (uchiyama et al. 2010 - ums10)
      call wave_bottom_momentum_production

! Add the wave stress balance (tao-taw) to the wind stress at "velocity grid nodes": !21-05-10
      call wave_surface_momentum_production

! Compute the "Sj & Sshear" terms:
      call wave_sj_sshear

      endif            !m°v°m> !17-11-19
      if(iwve==3) then   ! pour calculer en plus du cas 2, tow pour MUSTANG et le mettre
                         ! dans les fichiers offline
! Compute ubw and fw pour le calcul de tow:
      call wave_ubw_fw
      endif

      end subroutine wave_driver_s26

!.......................................................

      subroutine wave_read_interp_netcdf_s26(t_)
      implicit none
      double precision dlon_di_
      real*8 valmin_
      integer t_,ncid_,time_,loop_,iter_loop_,varnum_,loopvar_
      double precision grav_over_2pi_,beta_,delta_c_,alpha_,c_,kn_ &
                      ,weight_,delta_ij_
#ifdef synopsis
       subroutinetitle='wave_read_interp_netcdf_s26'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      ! TODO_CHECK THOM
      if (iwve==0) return

      allocate(ww3_var(ww3_imax,ww3_jmax))

      ww3flag1='first'
      do loopvar_=1,ww3_varmax
       open(unit=10,file=trim(tmpdirname)//'wavelist')
       read(10,'(a)')texte30 !21-11-19
       close(10)
       ww3varname(:)='reset'
       if(iwve==1) then !-----Michaud2012----> !21-11-19
        if(loopvar_==1)call wave_get_dir(t_)
        if(loopvar_==2)call wave_get_hs(t_)
        if(loopvar_==3)call wave_get_foc(t_)
        if(loopvar_==4)call wave_get_t(t_)
        if(loopvar_==5)call wave_get_taw(t_)
        if(loopvar_==6)call wave_get_two(t_)
        if(loopvar_==7)call wave_get_uss(t_)
        if(loopvar_==8)call wave_get_hsw(t_)
       endif            !-----Michaud2012---->
       if(iwve==2.or.iwve==3) then !-----McWilliams1999----> !21-11-19
        if(loopvar_==1)call wave_get_t(t_)
        if(loopvar_==2)call wave_get_uss(t_)
       endif            !-----McWilliams1999----> !21-11-19

      if(iwve==2) then
         if(loopvar_==3.and.flag_wavelength==1)call wave_get_l(t_) !10-11-21
      endif

      if(iwve==3) then
        if(loopvar_==3)call wave_get_hs(t_)
      endif


       ww3flag1='other'
      enddo !loopvar_

! Par defaut on considere que la longueur d'onde est presente dans le
! fichier (c.a.d. flag_wavelength=1 et ww3_varmax=3) mais si le champ
! n'a pas EtE trouvE alors flag_wavelength=0 et ici ww3_varmax devient 2:
      if(iwve==2.and.flag_wavelength==0)ww3_varmax=2 !21-11-21

      deallocate(ww3_var)
      if(allocated(ww3_short))deallocate(ww3_short)

! Compute the wave vector norm from the wave period
      if(flag_wavelength==0)call wave_kn_from_period !12-11-21

! Deduire le courant de Stokes 3D du courant 2D et du vecteur d'ondes:
      call wave_us3d_from_uss

! Deduire la fonction J barotrope du vecteur d'onde et de la hauteur significative:
      if(iwve==1)call wave_j_from_kn_hs !12-11-21

      end subroutine wave_read_interp_netcdf_s26

!.................................................

      subroutine wave_read_only(time_)
      implicit none
      integer time_,ncid_,filvalshort_
#ifdef synopsis
       subroutinetitle='wave_read_only'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0) then
       write(6,'(a,a,a,i0,a,a)')'File:',trim(texte80(1)),' time=',time_ &
       ,' ww3varname(1)= ',trim(ww3varname(1))
      endif

      if(ww3_type_grid==type_structured) then   !ssssss>
       varstart(1)=ww3zoom_istr ; varcount(1)=ww3_imax
       varstart(2)=ww3zoom_jstr ; varcount(2)=ww3_jmax
       varstart(3)=time_        ; varcount(3)=1
      endif                                     !ssssss>
      if(ww3_type_grid==type_unstructured) then !uuuuuu>
       varstart(1)=1     ; varcount(1)=ww3_imax
       varstart(2)=time_ ; varcount(2)=1 ! note ww3_jmax=1 si type unstructured
      endif                                     !uuuuuu>

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_) ; if(status/=0)stop 'erreur nf_open wave_read_interp_netcdf_s26'
                   status=nf_inq_varid(ncid_,trim(ww3varname(1)),var_id)
               if(status==0)varname_selected=trim(ww3varname(1)) !16-02-23 
      if(status/=0)status=nf_inq_varid(ncid_,trim(ww3varname(2)),var_id)
               if(status==0)varname_selected=trim(ww3varname(2)) !16-02-23 
      if(status/=0)status=nf_inq_varid(ncid_,trim(ww3varname(3)),var_id)
               if(status==0)varname_selected=trim(ww3varname(3)) !16-02-23 
      if(status/=0)status=nf_inq_varid(ncid_,trim(ww3varname(4)),var_id)
               if(status==0)varname_selected=trim(ww3varname(4)) !16-02-23 
      if(status/=0)status=nf_inq_varid(ncid_,trim(ww3varname(5)),var_id)
               if(status==0)varname_selected=trim(ww3varname(5)) !16-02-23 
      if(status/=0) then

       if(forcedstatus==0) then !-------->

       write(6,'(8(a,1x))')'Failed finding ww3varname:'    &
                                         ,trim(ww3varname(1)) &
                                         ,trim(ww3varname(2)) &
                                         ,trim(ww3varname(3)) &
                                         ,trim(ww3varname(4)) &
                                         ,trim(ww3varname(5)) &
                                         ,'in file',trim(texte80(1))
       stop 'erreur nf_inq_varid 253'

       else                     !-------->
        forcedstatus=2 ; status=nf_close(ncid_) ; return
       endif                    !-------->

      endif

      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'erreur nf_inq_var netcdf_s26'

      if(var_type==nf_real) then !rrrrrrrrrrrrrrrrr>
       status=nf_get_vara_real(ncid_,var_id                          &
                                    ,varstart(1:var_dims)            &
                                    ,varcount(1:var_dims)            &
                                    ,ww3_var(1:ww3_imax,1:ww3_jmax))

       if(status/=0) then
        write(6,*)'var_dims varcount ww3_imax ww3_jmax ',var_dims &
        ,varcount(1:var_dims),ww3_imax,ww3_jmax
        stop 'erreur nf_get_vara_real netcdf_s26'
       endif
       status=nf_get_att_real(ncid_,var_id,'_FillValue',ww3filval) ; if(status/=0)stop 'erreur _FillValue netcdf_s26'
       status=nf_get_att_real(ncid_,var_id,'scale_factor',var_scalefactor)
       if(status/=0)stop 'erreur scalefactor netcdf_s26'

! Bidouille pour graphique necessitant ww3_lon et ww3_lat...
!      if(ww3varname(1)(1:2)=='hs')then
!       call wave_read_lonlat(ncid_)
!       do i=1,ww3_imax
!        if(ww3_var(i,1)/=ww3filval)write(20,*)i,ww3_lon(i),ww3_lat(i),ww3_var(i,1)
!       enddo
!       stop 'mimi'
!      endif

       status=nf_close(ncid_)
       return
      endif                      !rrrrrrrrrrrrrrrrr>

      if(var_type==nf_short) then !ssssssssssssssss>

       if(.not.allocated(ww3_short))allocate(ww3_short(ww3_imax,ww3_jmax))
       status=nf_get_vara_int(ncid_,var_id                          &
                                   ,varstart(1:var_dims)            &
                                   ,varcount(1:var_dims)            &
                                   ,ww3_short(1:ww3_imax,1:ww3_jmax))
       if(status/=0)stop 'erreur nf_get_vara_int netcdf_s26'
       status=nf_get_att_int(ncid_,var_id,'_FillValue',filvalshort_)
       if(status/=0)stop 'erreur _FillValue netcdf_s26'
       status=nf_get_att_real(ncid_,var_id,'scale_factor',var_scalefactor)

       ww3filval=filvalshort_
       do j=1,ww3_jmax
       do i=1,ww3_imax
        ww3_var(i,j)=ww3_short(i,j)
       enddo
       enddo

       status=nf_close(ncid_)
       return
      endif                       !ssssssssssssssss>

      stop 'Erreur var_type subroutine wave_read_only'
      end subroutine wave_read_only

!.........................................................

      subroutine wave_interp_only(flag2_)
      implicit none
      double precision weight_,delta_ij_,dlon_di_
      character*5 flag2_
#ifdef synopsis
       subroutinetitle='wave_interp_only'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(ww3_type_grid==type_structured) then   !ssssss>

! BOUCHETROU: sous quelles conditions?
!     1. flag_ww3plugs=1
!     2. grille structuree
!     3. on n'est pas dans le cas ww3_on_sgrid
! Note : si le LSM de WW3 change en cours de simu (car il n'est pas fixe dans le temps) alors dans la
!        partie 'first' en suivant on s'adapte au fait que des trous n'ont pas ete bouchEs par wave_appli_bouchetrou
       if(flag_ww3plugs==1.and.txt_wavemodeltype/='ww3_on_sgrid') then !--plugs-bouchetrou-->
        if(.not.allocated(ww3plugs))call wave_init_ww3plugs    ! A l'initialisation identifier les boucheurs
                                    call wave_appli_bouchetrou ! Appliquer le bouchage
       endif                                                           !--plugs-bouchetrou-->

       if(ww3flag1=='first') then !-------->

! Faire les poids pour toutes les variables si premiere variable
       if(txt_wavemodeltype/='ww3_on_sgrid') then !>>>>>>>>>>

        do j=0,jmax+1
        do i=0,imax+1
         deci=ij2ww3_i(i,j) ; i1=int(deci) ; rapi=deci-i1 ; i1=i1-ww3zoom_istr+1
         decj=ij2ww3_j(i,j) ; j1=int(decj) ; rapj=decj-j1 ; j1=j1-ww3zoom_jstr+1
         x1=(1.-rapi)*(1.-rapj)
         x2=(1.-rapi)*    rapj
         x3=    rapi *(1.-rapj)
         x4=    rapi *    rapj
         x5=small1
         if(ww3_var(i1  ,j1  )==ww3filval)x1=0.
         if(ww3_var(i1  ,j1+1)==ww3filval)x2=0.
         if(ww3_var(i1+1,j1  )==ww3filval)x3=0.
         if(ww3_var(i1+1,j1+1)==ww3filval)x4=0.
         x0=1./(x1+x2+x3+x4+x5)
         xy_t(i,j,1)=x1*x0
         xy_t(i,j,2)=x2*x0
         xy_t(i,j,3)=x3*x0
         xy_t(i,j,4)=x4*x0
         xy_t(i,j,5)=x5*x0
        enddo
        enddo

       endif                                      !>>>>>>>>>>

       endif                     !-------->

       if(txt_wavemodeltype=='ww3_on_sgrid') then !oooooooo>
!      if(iteration3d==0)call wave_ww3onsgrid_checkmask
         do j=1,ww3_jmax ; do i=1,ww3_imax
          if(ww3_var(i,j)==ww3filval)ww3_var(i,j)=ww3mskval
         enddo           ; enddo
       endif                                      !oooooooo>

       if(flag2_=='basic') then !nnnnnnnnn>
! Cas particulier des grilles identiques:
        if(txt_wavemodeltype=='ww3_on_sgrid') then !oooooooo>
         if(ww3_on_sgrid_case==ww3_on_sgrid_casez1) then !z1z1z1>
          i0=0 ; if(iperiodicboundary)i0=2
          do j=0,jmax+1 ; do i=i0,imax+1-i0
! Note: dans le cas wavemodeltype=='ww3_on_sgrid' (grille identique) et compte tenu
! de ce que l'indexation de ww3_var debute a (1,1) il faut (quelque soit par%rank) faire
! correspondre (i,j)=(0,0) a (i1,j1)=(1,1):
           anyvar2d(i,j)=var_scalefactor*ww3_var(i+1-i0,j+1)
          enddo         ; enddo
         else                                            !z1z1z1> !z0z0z0>
          i1=1 ; j1=1
          if(par%timax(1)==0)i1=0
          if(par%tjmax(1)==0)j1=0
          do j=0,jmax+1 ; do i=0,imax+1
           anyvar2d(i,j)=var_scalefactor*ww3_var(min(max(i+i1,1),ww3_imax)  &
                                                ,min(max(j+j1,1),ww3_jmax))  
!          if(i==imax/2.and.j==jmax/2.and.mask_t(i,j,kmax)==1)write(10+par%rank,*)anyvar2d(i,j)
          enddo         ; enddo
         endif                                                    !z0z0z0>
        else                                       !oooooooo>
! Cas general des grilles differentes
         do j=0,jmax+1
         do i=0,imax+1

          i1=int( ij2ww3_i(i,j) )-ww3zoom_istr+1 ! conserve la parallelisation !11-04-13
          j1=int( ij2ww3_j(i,j) )-ww3zoom_jstr+1 ! conserve la parallelisation !11-04-13

            anyvar2d(i,j)=var_scalefactor*(                          &
                               xy_t(i,j,1)*ww3_var(i1  ,j1  )        &
                              +xy_t(i,j,2)*ww3_var(i1  ,j1+1)        &
                              +xy_t(i,j,3)*ww3_var(i1+1,j1  )        &
                              +xy_t(i,j,4)*ww3_var(i1+1,j1+1)        &
                              +xy_t(i,j,5)*ww3mskval )
         enddo
         enddo
        endif                                      !oooooooo>
       endif                    !nnnnnnnnn>

       if(flag2_=='angle') then !aaaaaaaaa>

! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p
        x1=       270.*deg2rad ! convention meteo
        x2=var_scalefactor*deg2rad
! Cas particulier des grilles identiques:
        if(txt_wavemodeltype=='ww3_on_sgrid') then !oooooooo>

        if(ww3_on_sgrid_case==ww3_on_sgrid_casez0) & !12-11-21
        stop 'Cas angle reste A faire pour ww3_on_sgrid_casez0'

         i0=0 ; if(iperiodicboundary)i0=2
         do j=0,jmax+1 ; do i=i0,imax+1-i0
! Note: dans le cas wavemodeltype=='ww3_on_sgrid' (grille identique) et compte tenu
! de ce que l'indexation de ww3_var debute a (1,1) il faut (quelque soit par%rank) faire
! correspondre (i,j)=(0,0) a (i1,j1)=(1,1):
          anyvar2d(i,j)=atan2(sin(x1-x2*ww3_var(i+1-i0,j+1))  &
                             ,cos(x1-x2*ww3_var(i+1-i0,j+1)))
! enfin prendre en compte l'orientation locale de la grille symphonie:
          dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
          if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
          if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
          anyvar2d(i,j)=anyvar2d(i,j)                 &
          -atan2(  lat_t(i+1,j)-lat_t(i-1,j)          & ! S grid rotation
                 ,dlon_di_*cos(lat_t(i,j)) )
         enddo         ; enddo

        else                                       !oooooooo>

         do j=0,jmax+1 ; do i=0,imax+1

          i1=int( ij2ww3_i(i,j) )-ww3zoom_istr+1 ! conserve la parallelisation !11-04-13
          j1=int( ij2ww3_j(i,j) )-ww3zoom_jstr+1 ! conserve la parallelisation !11-04-13

         anyvar2d(i,j)=atan2(  & ! retrouver l'angle à partir de sinus et cosinus
! interpolation de la composante sinus:
          xy_t(i,j,1)*sin(x1-x2*ww3_var(i1  ,j1  ))  &
         +xy_t(i,j,2)*sin(x1-x2*ww3_var(i1  ,j1+1))  &
         +xy_t(i,j,3)*sin(x1-x2*ww3_var(i1+1,j1  ))  &
         +xy_t(i,j,4)*sin(x1-x2*ww3_var(i1+1,j1+1))  &
        ,                                            &
! interpolation de la composante cosinus:
          xy_t(i,j,1)*cos(x1-x2*ww3_var(i1  ,j1  ))  &
         +xy_t(i,j,2)*cos(x1-x2*ww3_var(i1  ,j1+1))  &
         +xy_t(i,j,3)*cos(x1-x2*ww3_var(i1+1,j1  ))  &
         +xy_t(i,j,4)*cos(x1-x2*ww3_var(i1+1,j1+1)) )

! enfin prendre en compte l'orientation locale de la grille symphonie:
         dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
         if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
         if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
         anyvar2d(i,j)=anyvar2d(i,j)                 &
         -atan2(  lat_t(i+1,j)-lat_t(i-1,j)          & ! S grid rotation
                ,dlon_di_*cos(lat_t(i,j)) )

         enddo         ; enddo
        endif                                      !oooooooo>

       endif                    !aaaaaaaaa>

       if(txt_wavemodeltype=='ww3_on_sgrid') then !>>>>>>
                             call wave_obc_anyvar2d_jeq1 ! repare artefact champs ww3...
        if(iperiodicboundary)call obc_int_anyvar2d('za') ! wave_mpi_anyvar2d !21-03-15
       endif                                      !>>>>>>

      endif                                     !ssssss>

      if(ww3_type_grid==type_unstructured) then !uuuuuu>

       if(flag2_=='basic') then !bbbbbbbbbbbbbbb>

       open(unit=3,file=trim(tmpdirname)//'ww3grid_to_sgrid_'//dom_c//'.out')
        xy_t(:,:,1)=0. ; xy_t(:,:,2)=0.
!       xy_t(:,:,3)=0. ! decommenter pour calculer le nombre de points ww3 impliqués

        do k0=1,ww3_imax

         read(3,*,end=276)k,deci,decj,delta_ij_               ! node,i,j,delta_ij
! Note: deci decj sont "globaux" pour conservation mpi
         if(ww3_var(k,1)==ww3filval)ww3_var(k,1)=ww3mskval ! verifier l'arrondi à l'entier

! Projection de la valeur dans la bulle d'influence.
! Ponderation en fonction de la distance au centre de la bulle
!           i1=int( deci-delta_ij_     )-par%timax(1)
!           i2=int( deci+delta_ij_+0.5 )-par%timax(1)
!           j1=int( decj-delta_ij_     )-par%tjmax(1)
!           j2=int( decj+delta_ij_+0.5 )-par%tjmax(1)
            i1=nint( deci-delta_ij_)-par%timax(1)
            i2=nint( deci+delta_ij_)-par%timax(1)
            j1=nint( decj-delta_ij_)-par%tjmax(1)
            j2=nint( decj+delta_ij_)-par%tjmax(1)

            i1=max( i1 , -1     ) ! Attention ce bornage
            i2=min( i2 , imax+2 ) ! a pour consequence que les
            j1=max( j1 , -1     ) ! extremes -1,imax+2,-1,jmax+2 sont "des poubelles"
            j2=min( j2 , jmax+2 ) ! et donc ne doivent pas entrer ensuite dans les calculs

            x0=1./delta_ij_
            do j=j1,j2
            do i=i1,i2
             weight_=max(1.-x0*sqrt((i+par%timax(1)-deci)**2     &
                                   +(j+par%tjmax(1)-decj)**2),zero)
             xy_t(i,j,1)=xy_t(i,j,1)+weight_
             xy_t(i,j,2)=xy_t(i,j,2)+weight_*ww3_var(k,1)
!            if(weight_>0.)xy_t(i,j,3)=xy_t(i,j,3)+1. ! decommenter pour calculer le nombre de points ww3 impliqués

            enddo
            enddo

!        endif                          !..............>
        enddo ! fin de boucle sur k0
  276   close(3)
!       sum1=0. ;  sum2=0.
        x1=ww3mskval*var_scalefactor
        do j=0,jmax+1
        do i=0,imax+1
         if(xy_t(i,j,1)/=0.) then !-------->
          anyvar2d(i,j)=var_scalefactor*xy_t(i,j,2)/xy_t(i,j,1)
         else                     !-------->
          anyvar2d(i,j)=x1
         endif                    !-------->
!        if(mask_t(i,j,kmaxp1)==1) then
!        sum1=sum1+1. ; sum2=sum2+xy_t(i,j,3)
!        endif
        enddo
        enddo
!       write(6,*)'nb moyen de points ww3 impliqués=',sum2/sum1
       endif                    !bbbbbbbbbbbbbbb>

       if(flag2_=='angle') then !aaaaaaaaa>

       open(unit=3,file=trim(tmpdirname)//'ww3grid_to_sgrid_'//dom_c//'.out')
        xy_t(:,:,1)=0. ; xy_t(:,:,2)=0. ; xy_t(:,:,3)=0.

! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p
        x1=       270.*deg2rad ! convention meteo
        x2=var_scalefactor*deg2rad

        do k0=1,ww3_imax

         read(3,*,end=277)k,deci,decj,delta_ij_               ! node,i,j,delta_ij
! Note: deci decj sont "globaux" pour conservation mpi
         if(ww3_var(k,1)==ww3filval)ww3_var(k,1)=ww3mskval ! verifier l'arrondi à l'entier

! Projection de la valeur dans la bulle d'influence.
! Ponderation en fonction de la distance au centre de la bulle
!           i1=int( deci-delta_ij_     )-par%timax(1)
!           i2=int( deci+delta_ij_+0.5 )-par%timax(1)
!           j1=int( decj-delta_ij_     )-par%tjmax(1)
!           j2=int( decj+delta_ij_+0.5 )-par%tjmax(1)
            i1=nint( deci-delta_ij_)-par%timax(1)
            i2=nint( deci+delta_ij_)-par%timax(1)
            j1=nint( decj-delta_ij_)-par%tjmax(1)
            j2=nint( decj+delta_ij_)-par%tjmax(1)

            i1=max( i1 , -1     ) ! Attention ce bornage
            i2=min( i2 , imax+2 ) ! a pour consequence que les
            j1=max( j1 , -1     ) ! extremes -1,imax+2,-1,jmax+2 sont "des poubelles"
            j2=min( j2 , jmax+2 ) ! et donc ne doivent pas entrer ensuite dans les calculs

            x0=1./delta_ij_
            x3=sin(x1-x2*ww3_var(k,1))
            x4=cos(x1-x2*ww3_var(k,1))
            do j=j1,j2
            do i=i1,i2
             weight_=max(1.-x0*sqrt((i+par%timax(1)-deci)**2     &
                                   +(j+par%tjmax(1)-decj)**2),zero)
             xy_t(i,j,1)=xy_t(i,j,1)+weight_
             xy_t(i,j,2)=xy_t(i,j,2)+weight_*x3
             xy_t(i,j,3)=xy_t(i,j,3)+weight_*x4

            enddo
            enddo

        enddo ! fin de boucle sur k0
  277   close(3)

        x1=ww3mskval*var_scalefactor
        do j=0,jmax+1
        do i=0,imax+1
         if(xy_t(i,j,1)/=0.) then !-------->
          anyvar2d(i,j)=atan2( xy_t(i,j,2)/xy_t(i,j,1) &
                              ,xy_t(i,j,3)/xy_t(i,j,1) )
! enfin prendre en compte l'orientation locale de la grille symphonie:
         dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
         if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
         if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
         anyvar2d(i,j)=anyvar2d(i,j)         &
         -atan2(  lat_t(i+1,j)-lat_t(i-1,j)          & ! S grid rotation
                ,dlon_di_*cos(lat_t(i,j)) )
         else                     !-------->
          anyvar2d(i,j)=x1
         endif                    !-------->
        enddo
        enddo


       endif                    !aaaaaaaaa>

      endif                                     !uuuuuu>

      end subroutine wave_interp_only

!.......................................................

      subroutine wave_get_time_from_binrecfile(t_,time_)
      implicit none
      integer t_,time_,varnum_
#ifdef synopsis
       subroutinetitle='wave_get_time_from_binrecfile'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      wavefile_nextrec=max(wavefile_nextrec,1)

  687 open(unit=4,file=trim(texte30),access='direct',recl=reclen       &
                                    ,form='unformatted')
      read(4,rec=wavefile_nextrec)texte250,time_,wavefile_nextime ;
      texte80(1)=trim(texte250)
      close(4)

      if(t_==2.and.wavefile_nextime<=elapsedtime_now)then !pmxpmx>
! Si le temps du fichier proposE est avant le temps present, essayer
! le fichier suivant:
        wavefile_nextrec=wavefile_nextrec+1
        goto 687 
      endif                                               !pmxpmx>

!     if(t_==0.and.x1>elapsedtime_now) then !------->
      if(t_==0.and.wavefile_nextime>elapsedtime_now) then !------->
         write(6,'(a,a)')'texte30=',trim(texte30)
         write(6,'(a,a)')'texte80=',trim(texte80(1))
         write(6,*)'wavefile_nextrec=',wavefile_nextrec
         write(6,*)'time_=           ',time_
         write(6,*)'wavefile_nextime=',wavefile_nextime
         write(6,*)'elapsedtime_now= ',elapsedtime_now
         stop 'erreur1 wave_get_time_from_binrecfile'
      endif                                               !------->

!     if(t_==2)wavefile_nextime=x1

      end subroutine wave_get_time_from_binrecfile

!.......................................................
      subroutine wave_get_dpt
      implicit none
#ifdef synopsis
       subroutinetitle='wave_get_dpt'
       subroutinedescription='Get native bathy of ww3 model'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(ww3_var))                  &
               allocate(ww3_var(ww3_imax,ww3_jmax))

      ww3varname(:)='dpt'
      forcedstatus=1 ! continue even if the variable is not found

      call wave_read_only(1)

      if(forcedstatus==2) then !>>> variable not found
         forcedstatus=0 ; deallocate(ww3_var) ; flag_ww3_dpt=0 ; return
      endif                    !>>>

      flag_ww3_dpt=1
      ww3mskval=-9999./var_scalefactor
      call wave_interp_only('basic')

      if(.not.allocated(dpt_wave_t))  &
               allocate(dpt_wave_t(0:imax+1,0:jmax+1))

      do j=0,jmax+1 ; do i=0,imax+1
       if(anyvar2d(i,j)/=-9999.) then !>>>
        dpt_wave_t(i,j)=anyvar2d(i,j)
       else                           !>>>
        dpt_wave_t(i,j)=h_w(i,j)
       endif                          !>>>
       if(isnan(dpt_wave_t(i,j)))dpt_wave_t(i,j)=h_w(i,j)
      enddo ; enddo

      deallocate(ww3_var)
      end subroutine wave_get_dpt
!.......................................................
      subroutine wave_get_dir(t_)
      implicit none
      integer t_,time_
#ifdef synopsis
       subroutinetitle='wave_get_dir'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      ww3varname(1)='dir'
      ww3mskval=0.0 ! la valeur de masque n'a pas de sens pour dir mais on
                    ! ne laisse pas pour autant ww3mskval non definie
      call wave_read_only(time_)

      call wave_interp_only('angle')

      do j=0,jmax+1 ; do i=0,imax+1
        dir_wave_t(i,j,1)=anyvar2d(i,j)
      enddo ; enddo

      end subroutine wave_get_dir
!.......................................................
      subroutine wave_get_hs(t_)
      implicit none
      integer t_,time_
#ifdef synopsis
       subroutinetitle='wave_get_hs'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)
      ww3varname(1)='hs'
      ww3varname(2)='VHM0' !19-02-23

      call wave_read_only(time_)

      ww3mskval=0.0/var_scalefactor
      call wave_interp_only('basic')

      do j=0,jmax+1 ; do i=0,imax+1
        hs_wave_t(i,j,2)=anyvar2d(i,j)
      enddo ; enddo

      end subroutine wave_get_hs
!.......................................................
      subroutine wave_get_hsw(t_)
      implicit none
      integer t_,time_
#ifdef synopsis
       subroutinetitle='wave_get_hsw'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)
      ww3varname(1)='hsw'
      ww3varname(2)='wch'

      call wave_read_only(time_)

      ww3mskval=0.0/var_scalefactor
      call wave_interp_only('basic')

      do j=0,jmax+1 ; do i=0,imax+1
        hsw_wave_t(i,j,2)=max(anyvar2d(i,j),1. ) !03-02-15
      enddo ; enddo

      end subroutine wave_get_hsw
!.......................................................
      subroutine wave_get_foc(t_)
      implicit none
      integer t_,time_
#ifdef synopsis
       subroutinetitle='wave_get_foc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      ww3varname(1)='foc'
      call wave_read_only(time_)

      ww3mskval=0.0/var_scalefactor
      call wave_interp_only('basic')

      do j=0,jmax+1 ; do i=0,imax+1
        foc_wave_t(i,j,2)=anyvar2d(i,j)
      enddo ; enddo

      end subroutine wave_get_foc
!.......................................................
      subroutine wave_get_t(t_)
      implicit none
      integer t_,time_
      real*8 valmin_
#ifdef synopsis
       subroutinetitle='wave_get_t'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      ww3varname(1)='fp' !21-11-19
      ww3varname(2)='t0m1' 
      ww3varname(3)='t01'
      ww3varname(4)='VTPK' !14-02-23
      call wave_read_only(time_)

      ww3mskval=1./var_scalefactor  ; valmin_=1.  !17-02-15
      call wave_interp_only('basic')

!     if(ww3varname(1)=='fp') then !m°v°m>
      if(varname_selected=='fp') then !m°v°m> !16-02-23

! Cas lecture frequence
       do j=0,jmax+1 ; do i=0,imax+1
        t_wave_t(i,j,2)=max(1./anyvar2d(i,j),valmin_) !17-02-15
       enddo ; enddo

      else                         !m°v°m>

! Cas lecture periode
       do j=0,jmax+1 ; do i=0,imax+1
        t_wave_t(i,j,2)=max(anyvar2d(i,j),valmin_) !17-02-15
       enddo ; enddo
! Si on se trouve ici c'est que la variable fp n'a pas ete trouvee dans le
! fichier ww3. Comme la chose est a priori etrange on prefere stopper le
! modele. Supprimer le stop en suivant si on est d'accord avec le fait
! de lire la periode plutot que la frequence.
      endif                        !m°v°m>

      end subroutine wave_get_t
!.......................................................
      subroutine wave_get_taw(t_)
      implicit none
      integer t_,time_,loop_
#ifdef synopsis
       subroutinetitle='wave_get_taw'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      do loop_=1,2 ! boucle sur composantes u et v
       if(loop_==1)ww3varname(1)='utaw'
       if(loop_==2)ww3varname(1)='vtaw'

       call wave_read_only(time_)

       var_scalefactor=rho*var_scalefactor
       ww3mskval=0./var_scalefactor

       call wave_interp_only('basic')

       do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,loop_)=anyvar2d(i,j)
       enddo ; enddo

      enddo ! fin de boucle sur loop_

! S grid rotation:
      do j=0,jmax+1
      do i=0,imax+1
        tawx_wave_t(i,j,2)=gridrotcos_t(i,j)*anyvar3d(i,j,1)-gridrotsin_t(i,j)*anyvar3d(i,j,2)
        tawy_wave_t(i,j,2)=gridrotsin_t(i,j)*anyvar3d(i,j,1)+gridrotcos_t(i,j)*anyvar3d(i,j,2)
      enddo
      enddo

      end subroutine wave_get_taw
!.......................................................
      subroutine wave_get_two(t_)
      implicit none
      integer t_,time_,loop_
#ifdef synopsis
       subroutinetitle='wave_get_two'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      do loop_=1,2 ! boucle sur composantes u et v
       if(loop_==1)ww3varname(1)='utwo'
       if(loop_==2)ww3varname(1)='vtwo'

       call wave_read_only(time_)

       var_scalefactor=rho*var_scalefactor
       ww3mskval=0./var_scalefactor

       call wave_interp_only('basic')

       do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,loop_)=anyvar2d(i,j)
       enddo ; enddo

      enddo ! fin de boucle sur loop_

! S grid rotation:
      do j=0,jmax+1
      do i=0,imax+1
        twox_wave_t(i,j,2)=gridrotcos_t(i,j)*anyvar3d(i,j,1)-gridrotsin_t(i,j)*anyvar3d(i,j,2)
        twoy_wave_t(i,j,2)=gridrotsin_t(i,j)*anyvar3d(i,j,1)+gridrotcos_t(i,j)*anyvar3d(i,j,2)
      enddo
      enddo

! Borner two: !30-10-16
      do j=0,jmax+1 ; do i=0,imax+1
       x0=max(sqrt( twox_wave_t(i,j,2)**2+twoy_wave_t(i,j,2)**2),small1)
       x1=min(x0,30.)/x0
       twox_wave_t(i,j,2)=x1*twox_wave_t(i,j,2)
       twoy_wave_t(i,j,2)=x1*twoy_wave_t(i,j,2)
!      if(x0>30)write(10+par%rank,*)sqrt(twox_wave_t(i,j,2)**2+twoy_wave_t(i,j,2)**2),twox_wave_t(i,j,2),twoy_wave_t(i,j,2)
      enddo         ; enddo

      end subroutine wave_get_two
!.......................................................
      subroutine wave_get_uss(t_)
      implicit none
      integer t_,time_,loop_
#ifdef synopsis
       subroutinetitle='wave_get_uss'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      do loop_=1,2 ! boucle sur composantes u et v
       if(loop_==1) then
         ww3varname(1)='uuss'
         ww3varname(2)='VSDX' !16-02-23
       endif
       if(loop_==2) then
         ww3varname(1)='vuss'
         ww3varname(2)='VSDY' !16-02-23
       endif

       call wave_read_only(time_)

!     if(par%rank==0) then
!      i=1 ; j=36 ; k=1
!      write(60,*)'wave_get_uss1'
!      write(60,*)'i global',i+par%timax(1)
!      write(60,*)'iglb    ',iglb
!      write(60,*)'ww3_var(i+1,j+1)',ww3_var(i+1,j+1)
!     endif
!     if(par%rank==6) then
!      i=206   ; j=36 ; k=1
!      write(66,*)'wave_get_uss1'
!      write(66,*)'i global',i+par%timax(1)
!      write(66,*)'iglb    ',iglb
!      write(66,*)'ww3_var(i+1,j+1)',ww3_var(i+1,j+1)
!     endif

       ww3mskval=0./var_scalefactor
       call wave_interp_only('basic')


       do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,loop_)=anyvar2d(i,j)
       enddo ; enddo

      enddo ! fin de boucle sur loop_

! S grid rotation:
      do j=0,jmax+1
      do i=0,imax+1
        uss_wave_t(i,j,2)=gridrotcos_t(i,j)*anyvar3d(i,j,1)-gridrotsin_t(i,j)*anyvar3d(i,j,2) !12-09-18
        vss_wave_t(i,j,2)=gridrotsin_t(i,j)*anyvar3d(i,j,1)+gridrotcos_t(i,j)*anyvar3d(i,j,2)
      enddo
      enddo

      end subroutine wave_get_uss
!.......................................................

      subroutine wave_driver_s25
      use module_principal ; use module_parallele ; use module_sedw
      implicit none
      include 'netcdf.inc'
      real lw0_
#ifdef synopsis
       subroutinetitle='wave_driver_s25'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!------------------------------------------------------------------------------!
!$ Begin the interpolation of the wave model fields
!------------------------------------------------------------------------------!

! est on de part et d'autre d'une nouvelle echeance?
      call time_to_update_forcing_file(wavedt(2),wavedt(1))        !15-04-11
      x2=(elapsedtime_now-wavedt(2))/wavedt(1)                     !15-04-11
      rap_wave=x2-int(x2)

! si la reponse est oui alors il faut lire un nouveau fichier
! a l'etat initial il faut egalement lire des fichiers

      if(decision == 1            & ! de part & d'autre de l'echeance?    !15-04-11
        .or.iteration3d== kount0  & ! a l'etat initial?
                                      ) then !**************************>


      if(iteration3d==kount0) then  !>>>>>>>>>

        if(txt_wavemodeltype=='ww3')  call wave_read_interp_netcdf(0)
        if(txt_wavemodeltype=='ww3cg')call wave_read_interp_netcdf_ww3cg(0)

        do k=1,kmax
        do i=1,imax+1
        do j=1,jmax+1
         velstokes_u(i,j,k,1)=velstokes_u(i,j,k,2)/dz_u(i,j,k,1) ! cas velstokes(2) transport sur dz
         velstokes_v(i,j,k,1)=velstokes_v(i,j,k,2)/dz_v(i,j,k,1) ! cas velstokes(2) transport sur dz
        enddo
        enddo
        enddo

        if(txt_wavemodeltype=='ww3')  call wave_read_interp_netcdf(2)
        if(txt_wavemodeltype=='ww3cg')call wave_read_interp_netcdf_ww3cg(2)

      else                    !>>>>>>>>>

! sinon on ne lit que le fichier "devant" apres avoir recopié le champ precedent:
       call wave_moveforward
       if(txt_wavemodeltype=='ww3')  call wave_read_interp_netcdf(2)
       if(txt_wavemodeltype=='ww3cg')call wave_read_interp_netcdf_ww3cg(2)

      endif                   !>>>>>>>>>

      endif                                  !**************************>

! entre 2 echeances: la valeur instantannée est donnée par l'interpolation
! lineaire dans le temps des 2 echeances qui encadrent le temps present
      call wave_linear_in_time_02
      call wave_linear_in_time_01

! depth-averaged Stokes current:
      call wave_zaveraged_velstokes !09-04-13

! add the wave dissipation by bottom friction (uchiyama et al. 2010 - ums10)
      call wave_bottom_momentum_production

! add the wave stress balance (tao-taw) to the wind stress at "velocity grid nodes": !21-05-10
      call wave_surface_momentum_production

! Compute the "Sj & Sshear" terms:
      call wave_sj_sshear

!--------------------------------------------------------------------!
!$ end of the interpolation of the wave model fields
!--------------------------------------------------------------------!

      end subroutine wave_driver_s25

!......................................................................

      subroutine wave_find_misval(misval_)
      use module_parallele                                             !11-03-10
      use module_principal
      implicit none
      double precision misval_
#ifdef synopsis
       subroutinetitle='wave_find_misval'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! reperer les trous: valeurs=misval et en même temps impliquées
! dans l'interpolation bilineaires
! un trou est "marqué" en lui attribuant l'opposé de misval (sign -)
! ceci suppose que misval soit different de zero. la ligne suivante
! est une securite par rapport à cela:
      if(misval_==0)stop 'stop dans wave_find_misval'

      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille ww3:
      deci=1.+(lon_t(i,j)*180./pi-var_lonmin)/(var_lonmax-var_lonmin) &
                           *(max_x-1.)
      i1=int(deci)

      decj=1.+(lat_t(i,j)*180./pi-var_latmin)/(var_latmax-var_latmin) &
                           *(max_y-1.)
      j1=int(decj)

! interpolation sur 16 points:
      do j2=j1-1,j1+2
      do i2=i1-1,i1+2

! verifier que ww3 englobe la grille symphonie (uniquement à l'initialisation)
      if(i2<1    )stop 'stop wave_read_interp car hors dommaine 1'
      if(i2>max_x)stop 'stop wave_read_interp car hors dommaine 2'
      if(j2<1    )stop 'stop wave_read_interp car hors dommaine 3'
      if(j2>max_y)stop 'stop wave_read_interp car hors dommaine 4'
      if(anyvwave(i2,j2,1)==misval_)anyvwave(i2,j2,1)=-misval_

      enddo ! fin de boucle sur i2
      enddo ! fin de boucle sur j2


      endif                        !>------------->
      enddo
      enddo

! trouver la valeur la plus proche de chaque trou pour combler le
! trou:
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_wave_trou.ascii')              !11-03-10
      do j=1,max_y
      do i=1,max_x                                                     !26-08-09

! reperage trou à boucher:
       if(anyvwave(i,j,1)==-misval_)then !>\\\\\\\\\\\\>

         k0=1
  693    continue
! chercher aux alentours une valeur pour combler:
         do j1=max0(j-k0,1),min0(j+k0,max_y)
         do i1=max0(i-k0,1),min0(i+k0,max_x)
             if(anyvwave(i1,j1,1)/=-misval_  &
           .and.anyvwave(i1,j1,1)/= misval_)then !-.-.-.-.-.->

! une valeur de bouchage est trouvée. on note les indices et on passe au suivant
              write(3,*)i,j,i1,j1
              goto 692

             endif                                    !-.-.-.-.-.->
         enddo
         enddo
! si on passe par là, c'est que l'on n'a rien trouvé...
! on augmente la distance de prospection et on recommence
         k0=k0+1
         goto 693

       endif                                  !>\\\\\\\\\\\\>

  692 continue
      enddo
      enddo
      close(3)

      end subroutine wave_find_misval

!............................................................................

      subroutine wave_fill_misval(k_)
      use module_parallele                                             !11-03-10
      use module_principal
      implicit none
      integer k_
#ifdef synopsis
       subroutinetitle='wave_fill_misval'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_wave_trou.ascii')              !11-03-10
       do k10=1,999999
        read(3,*,end=725)i,j,i1,j1
        anyvwave(i,j,k_)=anyvwave(i1,j1,k_)
       enddo
       if(k10==999998)stop 'stop wave_fill_misval boucle trop petite'
  725 close(3)

      end subroutine wave_fill_misval

!............................................................................

      subroutine wave_poly3(ichoix)
      use module_principal
      implicit none
      double precision,dimension(4,4):: y_
      double precision,dimension(4)::   x_
      double precision,dimension(4)::   val_
      integer ichoix
#ifdef synopsis
       subroutinetitle='wave_poly3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolation de la variable brute:
! debut:
        if(ichoix.eq.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        x_(1)=i1-1.
        x_(2)=i1
        x_(3)=i1+1.
        x_(4)=i1+2.

        do j2=1,4
        j0=j1-2+j2


        y_(1,1)=anyvwave(i1-1,j0,1)
        y_(2,1)=anyvwave(i1  ,j0,1)
        y_(3,1)=anyvwave(i1+1,j0,1)
        y_(4,1)=anyvwave(i1+2,j0,1)

          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1)) &
                        /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo

! val_ est l'interpolation polynome p3 direction x
       val_(j2)=y_(1,4)
       do i3=2,4
        val_(j2)=val_(j2)*(deci-x_(i3))+y_(i3,5-i3)
       enddo

!       if(rapi<0.01.and.rapj<0.01)
!    & write(6,*)val_(j2)*var_scalefactor

!      write(68,*)deci,val_(j2)
!      do loop2=0,100
!      x0=x_(1)+0.01*loop2*3.
!      x1=y_(1,4)
!      do i3=2,4
!       x1=x1*(x0-x_(i3))+y_(i3,5-i3)
!      enddo
!      write(66,*)x0,x1
!      enddo ! fin de boucle loop2
!      write(67,*)x_(1),y_(1,1)
!      write(67,*)x_(2),y_(2,1)
!      write(67,*)x_(3),y_(3,1)
!      write(67,*)x_(4),y_(4,1)
!      stop 'albert'

        enddo ! fin de boucle sur j2

! interpolation direction j:
        x_(1)=j1-1.
        x_(2)=j1
        x_(3)=j1+1.
        x_(4)=j1+2.
        y_(1,1)=val_(1)
        y_(2,1)=val_(2)
        y_(3,1)=val_(3)
        y_(4,1)=val_(4)
          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1)) &
                        /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo
        x1=y_(1,4)
        do i3=2,4
         x1=x1*(decj-x_(i3))+y_(i3,5-i3)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolation de la variable brute:
! fin.
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolation du cosinus de la variable:
! debut:
        if(ichoix.eq.1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        x_(1)=i1-1.
        x_(2)=i1
        x_(3)=i1+1.
        x_(4)=i1+2.

        do j2=1,4
        j0=j1-2+j2


        y_(1,1)=cos(anyvwave(i1-1,j0,1))
        y_(2,1)=cos(anyvwave(i1  ,j0,1))
        y_(3,1)=cos(anyvwave(i1+1,j0,1))
        y_(4,1)=cos(anyvwave(i1+2,j0,1))

          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1)) &
                        /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo

! val_ est l'interpolation polynome p3 direction x
        val_(j2)=y_(1,4)
        do i3=2,4
         val_(j2)=val_(j2)*(deci-x_(i3))+y_(i3,5-i3)
        enddo

        enddo ! fin de boucle sur j2

! interpolation direction j:
        x_(1)=j1-1.
        x_(2)=j1
        x_(3)=j1+1.
        x_(4)=j1+2.
        y_(1,1)=val_(1)
        y_(2,1)=val_(2)
        y_(3,1)=val_(3)
        y_(4,1)=val_(4)
          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1)) &
                        /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo
        x1=y_(1,4)
        do i3=2,4
         x1=x1*(decj-x_(i3))+y_(i3,5-i3)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolation du cosinus de la variable:
! fin.
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolation du sinus de la variable:
! debut:
        if(ichoix.eq.2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        x_(1)=i1-1.
        x_(2)=i1
        x_(3)=i1+1.
        x_(4)=i1+2.

        do j2=1,4
        j0=j1-2+j2


        y_(1,1)=sin(anyvwave(i1-1,j0,1))
        y_(2,1)=sin(anyvwave(i1  ,j0,1))
        y_(3,1)=sin(anyvwave(i1+1,j0,1))
        y_(4,1)=sin(anyvwave(i1+2,j0,1))

          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1)) &
                        /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo

! val_ est l'interpolation polynome p3 direction x
        val_(j2)=y_(1,4)
        do i3=2,4
         val_(j2)=val_(j2)*(deci-x_(i3))+y_(i3,5-i3)
        enddo

        enddo ! fin de boucle sur j2

! interpolation direction j:
        x_(1)=j1-1.
        x_(2)=j1
        x_(3)=j1+1.
        x_(4)=j1+2.
        y_(1,1)=val_(1)
        y_(2,1)=val_(2)
        y_(3,1)=val_(3)
        y_(4,1)=val_(4)
          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1)) &
                        /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo
        x1=y_(1,4)
        do i3=2,4
         x1=x1*(decj-x_(i3))+y_(i3,5-i3)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolation du cosinus de la variable:
! fin.
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end subroutine wave_poly3


!/////////////////////////////////////////////////////////////////////!


      subroutine wave_read_interp_netcdf(t_)
      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'
      double precision dlon_di_
      real*8 valmin_            !29-01-13
      integer t_
#ifdef synopsis
       subroutinetitle='wave_read_interp_netcdf'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      stop 'Ne pas utiliser filiere ww3 mais ww3cg'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cas du modele ww3
! debut:
!     if(txt_wavemodeltype=='ww3') then
      if(txt_wavemodeltype=='ww3'.or.txt_wavemodeltype=='swn') then     !26-08-09
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! trouver le fichier et le numero du champ à l'interieur du fichier
!!    x2=(real(kount  )-wavedt(2))/wavedt(1)+0.5*t_
      x2=(elapsedtime_now-wavedt(2))/wavedt(1)+0.5*t_                !15-04-11
      nc=1+int(x2)/dataperwavefile           ! numero du fichier dans la liste binaire
      nc1=1+int(x2)-(nc-1)*dataperwavefile   ! numero du champs dans le fichier

!.........
! verifier que le numero de record existe bien dans le fichier:        !28-08-09
      texte30=trim(tmpdirname)//''//dom_c//'_listwave_max.txt'
      open(unit=3,file=texte30)
      read(3,*)k
      close(3)
      if(nc>k.or.nc<1)then
       write(6,*)
       write(6,*)'probleme avec le numero de record ',nc
       write(6,*)'il n''y a pas de champs de vagues à la date en cours'
       write(6,*)'conseils:'
       write(6,*)'1-verifier les dates de notebook_time & notebook_wave'
       write(6,*)'2-verifier qu''il y a assez de champs de vagues'
       stop ' en attendant stop dans modele_wave.f'
      endif
!.........

!     call allocate_forcages(1,8,ww3_imax,ww3_jmax,1)
!     call allocate_forcages(1,8,ww3_imax,ww3_jmax,30)
      call wave_allocate_ww3('a',ww3_imax,ww3_jmax,30)
      do loop2=1,13 !06-05-11
      loop1=loop2
      if(loop2==1)loop1=01
      if(loop2==2)loop1=02
      if(loop2==3)loop1=03
      if(loop2==4)loop1=04                               !02-09-10
      if(loop2==5)loop1=04
      if(loop2==6)loop1=05
      if(loop2==7)loop1=05
      if(loop2==8)loop1=06                               !02-09-10
      if(loop2==9)loop1=06
      if(loop2==10)loop1=07      !06-05-11
      if(loop2==11)loop1=08      !06-05-11
      if(loop2==12)loop1=09
      if(loop2==13)loop1=10

      write(texte30,'(i2,a7)')loop1,'.binrec' ! nom de la liste binaire
      texte30=trim(tmpdirname)//''//dom_c//'_listewave'//texte30

        open(unit=3,file=texte30              &             ! ouvre la liste binaire
                   ,access='direct'           &
                   ,recl=80*8                 &
                   ,form='unformatted')
        read(3,rec=nc)texte80(1) ! lit le nom du fichier du modele de vague à interpoler
        if(par%rank==0)write(6,'(a)')trim(texte80(1))
        close(3)

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec ouverture fichier ww3'

      texte4='none'
      if(loop2==1)texte4='t'
      if(loop2==2)texte4='hs'
      if(loop2==3)texte4='dir'
      if(loop2==4)texte4='utwo'
      if(loop2==5)texte4='vtwo'
      if(loop2==6)texte4='utaw'
      if(loop2==7)texte4='vtaw'
!      if(loop2==8)texte4='uusf'                     !02-09-10
!      if(loop2==9)texte4='vusf'
      if(loop2==8)texte4='uuss'                     !02-09-10
      if(loop2==9)texte4='vuss'
      if(loop2==10)texte4='l'                       !20-10-10
      if(loop2==11)texte4='foc'   !06-05-11
      if(loop2==12)texte4='wch'   !06-05-11
      if(loop2==13)texte4='j'


      if(texte4=='none')stop 'Variable ww3 non identifiee'

      status=nf_inq_varid(ncid1,texte4,var_id)
!! noter l'erreur dans les fichiers ww3 qui oblige eventuellement à fixer texte4='hs' !06-05-11
!      if(status/=0.and.trim(texte4)=='l') then
!       texte4='hs'
!       status=nf_inq_varid(ncid1,texte4,var_id)
!      endif
      if(status.ne.0)then
       write(*,'(a,a)')'echec id variable ww3 ',trim(texte4)
       stop 'stop in subroutine wave'
      endif

      status=nf_get_att_real(ncid1,var_id,'_fillValue',filval)
      if(status/=0)stop 'echec get _fillValue fichier ww3'

      status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
      if(status/=0)stop 'echec get scale_factor fichier ww3'

      status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
      if(status/=0)stop 'echec get add_offset fichier ww3'


!......................................
! lire la variable dans le fichier ww3:
! les fichiers .usf ont une variable f en plus        !02-09-10
!      if(loop2==8.or.loop2==9) then    !*-*-*-*-*-*-*-*->
!
!                   status=nf_inq_dimid(ncid1,'f',dim_z_id)
!      if(status/=0)stop 'erreur lecture dimid freq ww3'
!      status=nf_inq_dimlen(ncid1,dim_z_id,ww3_fmax)
!      if(status/=0)stop 'erreur lecture dimlen z ww3'
!
!      if(loop2==8) then    !--------------->
!      status=nf_inq_varid(ncid1,'f',var_id)
!      if(status.ne.0)stop 'echec id variable ww3 frequence'
!      status=nf_get_var_real(ncid1,var_id,freq_wave(1:ww3_fmax))
!      if(status.ne.0)stop 'echec lecture variable netcdf ww3 frequence'
!      endif                !--------------->
!
!      varstart(4)=nc1          ; varcount(4)=1        ! time
!      varstart(3)=1            ; varcount(3)=ww3_fmax ! frequency
!      varstart(2)=ww3zoom_jstr ; varcount(2)=ww3_jmax ! lat
!      varstart(1)=ww3zoom_istr ; varcount(1)=ww3_imax ! lon
!      status=nf_inq_varid(ncid1,texte4,var_id)
!      if(status.ne.0)stop 'echec id variable ww3 3'
!      status=nf_get_vara_int(ncid1,var_id,varstart(1:4),varcount(1:4),ww3_var2)
!      if(status.ne.0)stop 'echec lecture variable netcdf ww3 2'
!
!      else                               !*-*-*-*-*-*-*-*->

      varstart(3)=nc1          ; varcount(3)=1        ! time
      varstart(2)=ww3zoom_jstr ; varcount(2)=ww3_jmax ! lat
      varstart(1)=ww3zoom_istr ; varcount(1)=ww3_imax ! lon
      status=nf_inq_varid(ncid1,texte4,var_id)
      if(status.ne.0)stop 'echec id variable ww3 2'
      status=nf_get_vara_int(ncid1,var_id,varstart(1:3),varcount(1:3),ww3_var)
      if(status.ne.0)stop 'echec lecture variable netcdf ww3 1'

!     endif                              !*-*-*-*-*-*-*-*->

!..........................
! interpoler la variable 1:
      if(loop2==1) then !1-1-1-1-1-1-1->
      if(texte4/='t')stop 'erreur wave_read_interp t'

       valmin_=0.001 ; ww3mskval=1.                   !29-01-13

       do j=0,jmax+1
       do i=0,imax+1

       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1                              ! conservation de la parallelisation

! indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1                              ! conservation de la parallelisation


! a la premiere variable on calcule (et on stocke pour les variables suivantes)
! les coefficient de la ponderation:
      x1=(1.-rapi)*(1.-rapj)
      x2=(1.-rapi)*    rapj
      x3=    rapi *(1.-rapj)
      x4=    rapi *    rapj
!      if(ww3_var(i1  ,j1  )==nint(filval))x1=0.
!      if(ww3_var(i1  ,j1+1)==nint(filval))x2=0.
!      if(ww3_var(i1+1,j1  )==nint(filval))x3=0.
!      if(ww3_var(i1+1,j1+1)==nint(filval))x4=0.
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor               !29-01-13
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor              !29-01-13
      x0=x1+x2+x3+x4+small1
      xy_t(i,j,1)=x1/x0
      xy_t(i,j,2)=x2/x0
      xy_t(i,j,3)=x3/x0
      xy_t(i,j,4)=x4/x0

! interpolation bilineaire:
        t_wave_t(i,j,t_)=max(var_scalefactor*                     &
                              (xy_t(i,j,1)*ww3_var(i1  ,j1  )        &
                              +xy_t(i,j,2)*ww3_var(i1  ,j1+1)        &
                              +xy_t(i,j,3)*ww3_var(i1+1,j1  )        &
                              +xy_t(i,j,4)*ww3_var(i1+1,j1+1)),valmin_)
! note: pour eviter les divisions par zero dans les zones masquées on a appliqué
! un seuil sur t
       else                         !>------------->
        t_wave_t(i,j,t_)=ww3mskval                          !29-01-13

       endif                        !>------------->
       enddo
       enddo

! a l'initialisation du modele (c.a.d. t_==0) on verifie la coherence
! des grilles du modele de circulation et du modele de vague. si le
! recoupement est mauvais, on propose un reajustement de la grille soumis
! à l'approbation du l'utilisateur:
      if(t_==0)call wave_inq_grid_consistency
!     call wave_inq_grid_consistency    !20-10-10

      endif             !1-1-1-1-1-1-1->


!..........................
! interpoler la variable 2:
      if(loop2==2) then !2-2-2-2-2-2-2->
      if(texte4/='hs')stop 'erreur wave_read_interp hs'

       ww3mskval=0.01                             !29-01-13
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1                              ! conservation de la parallelisation

! indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1                              ! conservation de la parallelisation

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor


! pour hs on utilise rapi et rapj plutôt que xy_t afin d'avoir une C.l. en terre
! de type: variable=0
! interpolation bilineaire:
        hs_wave_t(i,j,t_)=var_scalefactor*                     &
                          ((1.-rapi)*(1.-rapj)*ww3_var(i1  ,j1  ) &
                          +(1.-rapi)*    rapj *ww3_var(i1  ,j1+1) &
                          +    rapi *(1.-rapj)*ww3_var(i1+1,j1  ) &
                          +    rapi *    rapj *ww3_var(i1+1,j1+1))
! note: pour hs on a appliqué une condition limite hs=0 aux frontieres solides
! en incluant des valeurs nulles (masque) dans l'interpolation

       else                         !>------------->
        hs_wave_t(i,j,t_)=ww3mskval                        !29-01-13
       endif                        !>------------->
       enddo
       enddo
      if(t_==0)  call wave_inq_grid_consistency
      endif             !2-2-2-2-2-2-2->

!..........................
! interpoler la variable 3:
      if(loop2==3) then !3-3-3-3-3-3-3->
      if(texte4/='dir')stop 'erreur wave_read_interp 3'
! attention en ce qui concerne les conventions d'angle en meteorologie. celles ci
! sont appliquées par ww3. en meteo, un vent venant du nord a une direction de 0°,
! un vent venant de l'est a une direction de 90°. par consequent il faut corriger
! cette convention en faisant 270°-angle. ensuite on convertit en radians.
! d'autre part on n'interpole pas un angle mais ses composantes trigonometriques.
! puis on retrouve un angle avec la fonction atan2
! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p

! application du facteur d'echelle, prise en compte des conventions meteo ,
! conversion degres->radians, prise en compte de la rotation de la grille
! symphonie:
! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p
       x1=       270.*deg2rad ! convention meteo
       x2=var_scalefactor*deg2rad
!      do j=1,ww3_jmax
!      do i=1,ww3_imax
!       ww3_var(i,j)=x1-x2*ww3_var(i,j)
!      enddo
!      enddo

       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        dir_wave_t(i,j,1)=atan2(  & ! retrouver l'angle à partir de sinus et cosinus
! interpolation de la composante sinus:
         xy_t(i,j,1)*sin(x1-x2*ww3_var(i1  ,j1  ))  &
        +xy_t(i,j,2)*sin(x1-x2*ww3_var(i1  ,j1+1))  &
        +xy_t(i,j,3)*sin(x1-x2*ww3_var(i1+1,j1  ))  &
        +xy_t(i,j,4)*sin(x1-x2*ww3_var(i1+1,j1+1))  &
       ,                                            &
! interpolation de la composante cosinus:
         xy_t(i,j,1)*cos(x1-x2*ww3_var(i1  ,j1  ))  &
        +xy_t(i,j,2)*cos(x1-x2*ww3_var(i1  ,j1+1))  &
        +xy_t(i,j,3)*cos(x1-x2*ww3_var(i1+1,j1  ))  &
        +xy_t(i,j,4)*cos(x1-x2*ww3_var(i1+1,j1+1)) )

! interpolation bi-polynome ordre 3:                                   !29-06-09
!       call wave_poly3(1) ! interpoler cosinus
!       x10=x1
!       call wave_poly3(2) ! interpoler sinus
!       dir_wave_t(i,j,1)=atan2(x1,x10) ! atan2(sin,cos) retrouve l'angle à partir de sin et cos

! enfin prendre en compte l'orientation locale de la grille symphonie:
       dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)       !14-02-13
       if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
       if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
        dir_wave_t(i,j,1)=dir_wave_t(i,j,1)                   &
        -atan2(  lat_t(i+1,j)-lat_t(i-1,j)                    & ! rotation symphonie
               ,dlon_di_*cos(lat_t(i,j)) )          !14-02-13
!              ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j)) )
       else                         !>------------->
        dir_wave_t(i,j,1)=0.                              !29-01-13
       endif                        !>------------->
       enddo
       enddo
           if(t_==0)call wave_inq_grid_consistency
      endif             !3-3-3-3-3-3-3->


!..........................
! interpoler la variable 4:
      if(loop2==4) then !4-4-4-4-4-4-4->
      if(texte4/='utwo')stop 'erreur wave_read_interp 4'

! interpolation on ocean model "t" grid nodes:
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        twox_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                        !>------------->
       enddo
       enddo
       if(t_==0) call wave_inq_grid_consistency
      endif             !4-4-4-4-4-4-4->

!..........................
! interpoler la variable 5:
      if(loop2==5) then !5-5-5-5-5-5-5->
      if(texte4/='vtwo')stop 'erreur wave_read_interp 5'

! interpolation on ocean model "t" grid nodes:
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        twoy_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                      !>------------->
       enddo
       enddo
      if(t_==0) call wave_inq_grid_consistency
      endif             !5-5-5-5-5-5-5->

!..........................
! interpoler la variable 6:
      if(loop2==6) then !6-6-6-6-6-6-6->
      if(texte4/='utaw')stop 'erreur wave_read_interp utaw'

       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        tawx_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))



       endif                      !>------------->
       enddo
       enddo

       if(t_==0)  call wave_inq_grid_consistency
      endif             !6-6-6-6-6-6-6->

!..........................
! interpoler la variable 7:
      if(loop2==7) then !7-7-7-7-7-7-7->
      if(texte4/='vtaw')stop 'erreur wave_read_interp vtaw'

! interpolation on ocean model "t" grid nodes:
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        tawy_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                      !>------------->
       enddo
       enddo
         if(t_==0) call wave_inq_grid_consistency
      endif             !7-7-7-7-7-7-7->

!      if(loop2==8) then !8-8-8-8-8-8-8-8-8-8->
!      if(texte4/='uusf')stop 'erreur wave_read_interp usf'
!       do freq=1,ww3_fmax
!       do j=0,jmax+1
!       do i=0,imax+1
!       if(mask_t(i,j,kmax+1)==1) then !>------------->
!
!! indice decimale i (longitude) dans grille aladin:
!      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1
!
!! indice decimale j (latitude) dans grille aladin:
!      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
!      if(ww3_var(i1  ,j1  ,freq)==nint(filval))ww3_var2(i1  ,j1  ,freq)=small1
!      if(ww3_var(i1  ,j1+1,freq)==nint(filval))ww3_var2(i1  ,j1+1,freq)=small1
!      if(ww3_var(i1+1,j1  ,freq)==nint(filval))ww3_var2(i1+1,j1  ,freq)=small1
!      if(ww3_var(i1+1,j1+1,freq)==nint(filval))ww3_var2(i1+1,j1+1,freq)=small1
!
!              usf_wave_t(i,j,freq)=var_scalefactor*                    &
!                                 (xy_t(i,j,1)*ww3_var2(i1  ,j1  ,freq) &
!                                 +xy_t(i,j,2)*ww3_var2(i1  ,j1+1,freq) &
!                                 +xy_t(i,j,3)*ww3_var2(i1+1,j1  ,freq) &
!                                 +xy_t(i,j,4)*ww3_var2(i1+1,j1+1,freq))
!
!       endif                      !>------------->
!       enddo
!       enddo
!       enddo
!     if(t_==0)  call wave_inq_grid_consistency
!      endif             !8-8-8-8-8-8-8-8-8->


!..........................
! interpoler la variable 8:
      if(loop2==8) then !8-8-8-8-8-8-8-8-8-8->
      if(texte4/='uuss')stop 'erreur wave_read_interp uss'
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->
! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=small1
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=small1
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=small1
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=small1

          uss_wave_t(i,j,1)=var_scalefactor*                    &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))
       endif                      !>------------->
       enddo
       enddo
           if(t_==0)call wave_inq_grid_consistency
      endif             !8-8-8-8-8-8-8-8-8->

!!..........................
!! interpoler la variable 9:
!      if(loop2==9) then !9-9-9-9-9-9-9-9-9->
!      if(texte4/='vusf')stop 'erreur wave_read_interp vsf'
!
!! interpolation on ocean model "t" grid nodes:
!       do freq=1,ww3_fmax
!       do j=0,jmax+1
!       do i=0,imax+1
!       if(mask_t(i,j,kmax+1)==1) then !>------------->
!
!! indice decimale i (longitude) dans grille aladin:
!      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1
!
!! indice decimale j (latitude) dans grille aladin:
!      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
!
!      if(ww3_var(i1  ,j1  ,freq)==nint(filval))ww3_var2(i1  ,j1  ,freq)=small1
!      if(ww3_var(i1  ,j1+1,freq)==nint(filval))ww3_var2(i1  ,j1+1,freq)=small1
!      if(ww3_var(i1+1,j1  ,freq)==nint(filval))ww3_var2(i1+1,j1  ,freq)=small1
!      if(ww3_var(i1+1,j1+1,freq)==nint(filval))ww3_var2(i1+1,j1+1,freq)=small1
!
!              vsf_wave_t(i,j,freq)=var_scalefactor*                    &
!                                 (xy_t(i,j,1)*ww3_var2(i1  ,j1  ,freq) &
!                                 +xy_t(i,j,2)*ww3_var2(i1  ,j1+1,freq) &
!                                 +xy_t(i,j,3)*ww3_var2(i1+1,j1  ,freq) &
!                                 +xy_t(i,j,4)*ww3_var2(i1+1,j1+1,freq))
!
!       endif                      !>------------->
!       enddo
!       enddo
!       enddo
!      if(t_==0) call wave_inq_grid_consistency
!      endif             !9-9-9-9-9-9-9-9-9-9->

!..........................
! interpoler la variable 9:
      if(loop2==9) then !9-9-9-9-9-9-9-9-9->
      if(texte4/='vuss')stop 'erreur wave_read_interp vss'

! interpolation on ocean model "t" grid nodes:
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=small1
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=small1
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=small1
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=small1

              vss_wave_t(i,j,1)=var_scalefactor*                    &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                      !>------------->
       enddo
       enddo
        if(t_==0) call wave_inq_grid_consistency
      endif             !9-9-9-9-9-9-9-9-9-9->

!..........................
! interpoler la variable 10:
      if(loop2==10) then !10-10-10-10-10-10-10-10-10-10-10->
! noter prise en compte erreur attribut dans fichier ww3: !06-05-11
      if(texte4/='l'.and.texte4/='hs') then
        write(*,*)'erreur nom variable ww3: ',trim(texte4)
        stop 'erreur in subroutine wave'
      endif

       const2=2.*pi
       ww3mskval=10.

       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->


! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1


      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor


                                x1=var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

        k_wave_t(i,j,t_)=const2/x1
        kx_wave_t(i,j,t_)=k_wave_t(i,j,t_)*cos(dir_wave_t(i,j,1))
        ky_wave_t(i,j,t_)=k_wave_t(i,j,t_)*sin(dir_wave_t(i,j,1))
       else                         !>------------->
        k_wave_t(i,j,t_)=const2/ww3mskval
        kx_wave_t(i,j,t_)=k_wave_t(i,j,t_)*cos(dir_wave_t(i,j,1))
        ky_wave_t(i,j,t_)=k_wave_t(i,j,t_)*sin(dir_wave_t(i,j,1))

       endif                        !>------------->
       enddo
       enddo
         if(t_==0)    call wave_inq_grid_consistency
      endif            !10-10-10-10-10-10-10-10-->

!..........................
! interpoler la variable 11:
      if(loop2==11) then !11-11-11-11-11-11-11->
      if(texte4/='foc')stop 'erreur wave_read_interp foc'

! interpolation on ocean model "t" grid nodes:
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        foc_wave_t(i,j,t_)=var_scalefactor*                   &     !07-05-11
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                        !>------------->
       enddo
       enddo
         if(t_==0) call wave_inq_grid_consistency
      endif             !11-11-11-11-11-11-11->

!..........................
! interpoler la variable 12:
      if(loop2==12) then !12-12-12-12-12-12-12-->
      if(texte4/='wch') stop 'erreur wave_read_interp wch'

       ww3mskval=0.01                                              !29-01-13
       do j=0,jmax+1
       do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !>------------->


! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1


      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor          !29-01-13
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor


! interpolation bilineaire:
        hsw_wave_t(i,j,t_)=var_scalefactor*                      &     !07-05-11
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))
! note: pour hs on a appliqué une condition limite hs=0 aux frontieres solides
! en incluant des valeurs nulles (masque) dans l'interpolation

       else                         !>------------->                        !29-01-13
        hsw_wave_t(i,j,t_)=ww3mskval                                        !29-01-13
       endif                        !>------------->
       enddo
       enddo
         if(t_==0)  call wave_inq_grid_consistency
      endif             !12-12-12-12-12-12-12->

! interpoler la variable 13:
      if(loop2==13) then !13-13-13-13-13-13-13-->
      if(texte4/='j')stop 'erreur wave_read_interp j'

       do j=0,jmax+1
       do i=0,imax+1

       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1

! indice decimale j (latitude) dans grille aladin:
      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

! interpolation bilineaire:
        j_wave_t(i,j,t_)=var_scalefactor*                        &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                        !>------------->
       enddo
       enddo
      if(t_==0) call wave_inq_grid_consistency
      endif             !13-13-13-31-13-13-13->

! fermer le fichier netcdf ww3
      status=nf_close(ncid1)
      if(status/=0)stop 'erreur fermeture fichier ww3'
      enddo ! fin de boucle du loop1
!     call allocate_forcages(2,8,0,0,0) ! desalouer memoire
      call wave_allocate_ww3('d',0,0,0) ! desalouer memoire

!.........................................................


! grid rotation:
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !////////>

        x1=tawx_wave_t(i,j,t_)
        x2=tawy_wave_t(i,j,t_)
        tawx_wave_t(i,j,t_)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
        tawy_wave_t(i,j,t_)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2

        x1=twox_wave_t(i,j,t_)
        x2=twoy_wave_t(i,j,t_)
        twox_wave_t(i,j,t_)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
        twoy_wave_t(i,j,t_)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2

      endif                        !////////>
      enddo
      enddo
      do j=0,jmax+1
      do i=0,imax+1

       x1=uss_wave_t(i,j,1)
       x2=vss_wave_t(i,j,1)
       uss_wave_t(i,j,1)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
       vss_wave_t(i,j,1)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2

      enddo
      enddo

!      do freq=1,ww3_fmax
!      do j=0,jmax+1
!      do i=0,imax+1!!!

!       x1=usf_wave_t(i,j,freq)
!       x2=vsf_wave_t(i,j,freq)
!       usf_wave_t(i,j,freq)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
!       vsf_wave_t(i,j,freq)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2

!      enddo
!      enddo
!      enddo
!.........................................................


!! Compute the total 3d stokes velocity using the frequential components
!! and their respective z profil
!      do k=1,kmax
!      do j=0,jmax+1
!      do i=0,imax+1
!       anyv3d(i,j,k,1)=0.
!       anyv3d(i,j,k,2)=0.
!      enddo
!      enddo
!      enddo

      const1=0.5*grav/pi
      const2=2.*pi
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !mmmmmmmmmmmmmm>

!      do freq=1,ww3_fmax

!! x1 is the wave vector kn(freq)
!! x0 is the wave celerity at the frequency freq_wave(freq)
!! on veut résoudre x=f(x)=alpha*tanh(beta/x)
!! on dit que f(x+dx)=x+dx
!!           =f(x)+dx.f'(x)=x+dx
!!            dx=(x-f(x))/(f'(x)-1)
!! on estime une première valeur de x=x0
!      x3=const2*freq_wave(freq)*max(hz_w(i,j,1),0.01d0)         !02-02-11 beta
!      x0=const1/freq_wave(freq) ! c=g*t/2pi first guess for large depth
!
!      if(x3/x0>200.) then !fgfgfgfgfgfgfg>      !07-05-11
! Cas c=first guess valide:
!       x1=x0
!
!      else                !fgfgfgfgfgfgfg>
!
!! Cas calcul itératif:
!       x12=0.01                  ! dx, on estime une première valeur de dx
!       x1=x0+x12                 ! x0+dx
!         do k0=1,200
!          x12=(x1-x0*tanh(x3/x1))/(-x0*x3/(x1*cosh(x3/x1))**2-1.)
!          x1=max(x1+x12,1.d-4)
!          if(abs(x12)<0.01)goto 25
!         enddo
!         write(*,*)i,j,x1,x12
!       stop 'wave: le calcul de kn ne converge pas'
!   25  continue
!
!      endif               !fgfgfgfgfgfgfg>


! alternative si on ne dispose pas des vitesses de stokes sur les 30 fréquences:
       x1=sqrt(kx_wave_t(i,j,t_)**2+ky_wave_t(i,j,t_)**2+small1)

       if((h_w(i,j)+depth_w(i,j,kmax+1))*x1<6.)  then !hhhhhhhhhhhhhhhhh>

! Cas "petites profondeurs":
       if(h_w(i,j)+depth_w(i,j,kmax+1)>0.01) then !ssssssssssss> ! si dz=0 us=0
       x4=cosh(2.*x1*(h_w(i,j)+depth_w(i,j,kmax+1))) ! x4=cosh(2kn(h+ssh))
       do k=kmin_w(i,j),kmax

! méthode 1: valeur ponctuelle du profil:
!       x3=cosh(2.*x1*(depth_t(i,j,k)+h_w(i,j)))/x4
! méthode 2: profil moyénné sur dz
        x3=0.5/x1/x4/(depth_w(i,j,k+1)-depth_w(i,j,k))           &
           *( sinh(2.*x1*(depth_w(i,j,k+1)+h_w(i,j)))            &
             -sinh(2.*x1*(depth_w(i,j,k  )+h_w(i,j))))
!       x3=1. ! cas academique profil constant
!        anyv3d(i,j,k,1)=anyv3d(i,j,k,1)+usf_wave_t(i,j,freq)*x3
!        anyv3d(i,j,k,2)=anyv3d(i,j,k,2)+vsf_wave_t(i,j,freq)*x3
        anyv3d(i,j,k,1)=uss_wave_t(i,j,1)*x3
        anyv3d(i,j,k,2)=vss_wave_t(i,j,1)*x3

       enddo ! fin de boucle sur k
       endif                                      !ssssssssssss>

       else                                          !hhhhhhhhhhhhhhhhh>

! Cas "grandes profondeurs":
       do k=kmin_w(i,j),kmax

! méthode 1: valeur ponctuelle du profil exponentiel
!       x3=exp(2.*x1*(depth_t(i,j,k)-ssh_int_w(i,j,1)))
! méthode 2: profil exponentiel moyenné sur dz:
        x3=0.5/x1/(depth_w(i,j,k+1)-depth_w(i,j,k))                    &
                    *( exp(2.*x1*(depth_w(i,j,k+1)-ssh_int_w(i,j,1)))  &
                      -exp(2.*x1*(depth_w(i,j,k  )-ssh_int_w(i,j,1))))
!       x3=1. ! cas academique profil constant
!        anyv3d(i,j,k,1)=anyv3d(i,j,k,1)+usf_wave_t(i,j,freq)*x3
!        anyv3d(i,j,k,2)=anyv3d(i,j,k,2)+vsf_wave_t(i,j,freq)*x3
        anyv3d(i,j,k,1)=uss_wave_t(i,j,1)*x3
        anyv3d(i,j,k,2)=vss_wave_t(i,j,1)*x3
       enddo ! fin de boucle sur k

       endif                                         !hhhhhhhhhhhhhhhhh>

!      enddo ! fin de boucle sur freq

      endif                          !mmmmmmmmmmmmmm>
      enddo
      enddo


      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1

       velstokes_u(i,j,k,2)=                           &
              dz_u(i,j,k,1)*                           &
            mask_u(i,j,k)*0.5*( anyv3d(i  ,j  ,k,1)    &
                               +anyv3d(i-1,j  ,k,1))
       velstokes_v(i,j,k,2)=                           &
              dz_v(i,j,k,1)*                           &
            mask_v(i,j,k)*0.5*( anyv3d(i  ,j  ,k,2)    &
                               +anyv3d(i  ,j-1,k,2))

      enddo
      enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cas du modele ww3
! fin.
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end subroutine wave_read_interp_netcdf

!.....................................................................

      subroutine wave_inq_grid_consistency
      use module_principal
      use module_parallele
      use module_global
      implicit none
#ifdef synopsis
       subroutinetitle='wave_inq_grid_consistency'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! a l'initialisation du modele (c.a.d. t_==0) on verifie la coherence
! des grilles du modele de circulation et du modele de vague. si le
! recoupement est mauvais, on propose un reajustement de la grille soumis
! à l'approbation du l'utilisateur:

      ksecu=0
       do j=1,jmax
       do i=1,imax
       if(mask_t(i,j,kmax+1)==1) then !>------------->
        if(xy_t(i,j,1)+xy_t(i,j,2)+xy_t(i,j,3)+xy_t(i,j,4)==0.) then
         mask_t(i,j,kmax+1)=0
         write(*,*)'point ',i,j,' est masqué'
        write(*,*)'test pt ',i+par%timax(1),j+par%tjmax(1),' est masqué'
         ksecu=1
        endif
       endif                          !>------------->
       enddo
       enddo

       if(ksecu==1) then !stop-stop-stop>
       call allocate_global('a','glob_mask           ',iglb,jglb,0)
        do j=1,jmax
        do i=1,imax
         glob_mask(i+par%timax(1),j+par%tjmax(1))=mask_t(i,j,kmax+1)
        enddo
        enddo
#ifdef parallele
      ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
      ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
      lb2=lbound(glob_mask) ;  ub2=ubound(glob_mask)
      call  par_gatherall_2d(glob_mask,lb2,ub2,ind,par%nbdom)
#endif
       if(par%rank==0) then !r0-r0-r0-r0>
        write(*,*)'il existe des points de mer non encadrés par la'
        write(*,*)'grille du modele de vague. on masque ces points.'
        write(*,*)'un masque terre/mer compatible avec la grille'
        write(*,*)'du modele de vague est archivé dans:'
        texte30=trim(tmpdirname)//'mask_ww3compatible.txt'
        open(unit=3,file=texte30)
        do i=1,iglb
        write(3,'(1000i1)')(glob_mask(i,j),j=1,jglb)
        enddo
        close(3)
        write(*,'(a)')texte30
       endif                !r0-r0-r0-r0>
      call allocate_global('d','glob_mask           ',0,0,0)       !08-11-09
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
        stop ' stop dans wave'
        endif             !stop-stop-stop>


      end subroutine wave_inq_grid_consistency

!........................................................................

      subroutine  wave_grid_extraction_ww3
      use module_principal
      implicit none
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='wave_grid_extraction_ww3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec fichier'

                   status=nf_inq_dimid(ncid1,'nav_lon', dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'longitude',dim_x_id)
      if(status/=0)stop 'erreur lecture dimid longitude ww3'
                   status=nf_inq_dimid(ncid1,'nav_lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'latitude',dim_y_id)
      if(status/=0)stop 'erreur lecture dimid latitude ww3'

      status=nf_inq_dimlen(ncid1,dim_x_id,ww3_imax)
      if(status/=0)stop 'erreur lecture dimlen x ww3'
      status=nf_inq_dimlen(ncid1,dim_y_id,ww3_jmax)
      if(status/=0)stop 'erreur lecture dimlen y ww3'

!     call allocate_forcages(1,8,ww3_imax,ww3_jmax,1)
      call wave_allocate_ww3('a',ww3_imax,ww3_jmax,1)

! lire les longitudes et latitudes:
      varstart(1)=1 ; varcount(1)=ww3_imax
                   status=nf_inq_varid(ncid1,'nav_lon',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      if(status/=0)stop 'erreur var_id longitude ww3'
      status=nf_get_vara_real(ncid1,var_id,varstart(1),varcount(1),ww3_lon)
      if(status/=0)stop 'erreur lecture longitude ww3'

      varstart(1)=1 ; varcount(1)=ww3_jmax
                   status=nf_inq_varid(ncid1,'nav_lat',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      if(status/=0)stop 'erreur var_id latitude ww3'
      status=nf_get_vara_real(ncid1,var_id,varstart(1),varcount(1),ww3_lat)
      if(status/=0)stop 'erreur lecture latitude ww3'

! avant de redecouper la grille ww3 en sous domaines on memorise les dimensions natives:
      ww3full_imax=ww3_imax ; ww3full_jmax=ww3_jmax      ! conserve la parallelisation
      ww3_lonmin=ww3_lon(1)
      ww3_lonmax=ww3_lon(ww3full_imax)
      ww3_latmin=ww3_lat(1)
      ww3_latmax=ww3_lat(ww3full_jmax)
      ww3_dlon=(ww3_lonmax-ww3_lonmin)/(ww3full_imax-1) ! conserve la parallisation
      ww3_dlat=(ww3_latmax-ww3_latmin)/(ww3full_jmax-1) ! conserve la parallisation

! a partir des min & max des (lon,lat) des 2 grilles, reduire la zone d'extraction:
       ww3zoom_istr=1 ; ww3zoom_iend=ww3_imax ! first guess
       do i=1,ww3_imax-1
        if(ww3_lon(i)<lonmin.and.ww3_lon(i+1)>=lonmin)ww3zoom_istr=i
        if(ww3_lon(i)<lonmax.and.ww3_lon(i+1)>=lonmax)ww3zoom_iend=i+1
       enddo
       ww3zoom_jstr=1 ; ww3zoom_jend=ww3_jmax ! first guess
       do j=1,ww3_jmax-1
        if(ww3_lat(j)<latmin.and.ww3_lat(j+1)>=latmin)ww3zoom_jstr=j
        if(ww3_lat(j)<latmax.and.ww3_lat(j+1)>=latmax)ww3zoom_jend=j+1
       enddo
! on elargit un peu plus pour boucher les trous sans être restreint par la zone
! d'extraction (pour conserver la parallelisation)
       i0=5 ! elargissement à 1 ligne / 1 colonne supplementaire - repere 1498
       ww3zoom_istr=max0(ww3zoom_istr-i0,1)
       ww3zoom_iend=min0(ww3zoom_iend+i0,ww3full_imax)
       ww3zoom_jstr=max0(ww3zoom_jstr-i0,1)
       ww3zoom_jend=min0(ww3zoom_jend+i0,ww3full_jmax)

       ww3_imax=ww3zoom_iend-ww3zoom_istr+1
       ww3_jmax=ww3zoom_jend-ww3zoom_jstr+1

!     call allocate_forcages(2,8,0,0,0)
      call wave_allocate_ww3('d',0,0,0)


      status=nf_close(ncid1)


      end subroutine wave_grid_extraction_ww3

!...............................................................................

      subroutine wave_grid_extraction_ww3cg
      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'
      double precision deci_min_,deci_max_       &
                      ,decj_min_,decj_max_
      integer ncid_
#ifdef synopsis
       subroutinetitle='wave_grid_extraction_ww3cg'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
      if(status.ne.0)stop 'echec fichier'

                   status=nf_inq_dimid(ncid_,'nav_lon', dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'x',dim_x_id)
      if(status/=0) then !---------------->
          status=nf_inq_dimid(ncid_,'node',dim_x_id)
          if(status==0) then !----------------->
           ww3_type_grid=type_unstructured
          endif              !----------------->
      endif              !---------------->
      if(status/=0)stop 'erreur lecture dimid x ww3cg'


      if(ww3_type_grid==type_structured) then !------------------->

                   status=nf_inq_dimid(ncid_,'nav_lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'y',dim_y_id)
      if(status/=0)stop 'erreur lecture dimid y ww3cg'

      status=nf_inq_dimlen(ncid_,dim_x_id,ww3_imax)
      if(status/=0)stop 'erreur lecture dimlen x ww3cg'
      status=nf_inq_dimlen(ncid_,dim_y_id,ww3_jmax)
      if(status/=0)stop 'erreur lecture dimlen y ww3cg'

      endif                                   !------------------->

      if(ww3_type_grid==type_unstructured) then !------------------->
       status=nf_inq_dimlen(ncid_,dim_x_id,ww3_imax)
       if(status/=0)stop 'erreur lecture dimlen x ww3cg'
       ww3_jmax=1
!      status=nf_inq_dimid(ncid_,'element',dim_y_id)
!      if(status/=0)stop 'error nf_inq_dimid element'
!      status=nf_inq_dimlen(ncid_,dim_y_id,ww3_element)
!      if(status/=0)stop 'error nf_inq_dimlen element'
      endif                                     !------------------->

! lire les longitudes et latitudes:
      ww3zoom_istr=1 ; ww3zoom_jstr=1 ; call wave_read_lonlat(ncid_)

!     if(.not.allocated(ww3_short))allocate(ww3_short(ww3_element,3))
!     status=nf_inq_varid(ncid_,'tri',var_id)
!     if(status/=0)stop 'error nf_inq_varid tri'
!     varstart(1:2)=1 ; varcount(2)=ww3_element ; varcount(1)=3
!     status=nf_get_vara_int(ncid_,var_id         &
!                                 ,varstart(1:2)  &
!                                 ,varcount(1:2)  &
!                                 ,ww3_short(:,:))
!     if(status/=0)stop 'error nf_get_vara_int tri'
!     deallocate(ww3_short)

! Tableaux de passage de la grille S vers la grille WW3
! Placé ici les tableaux de passage donne la position sur la grille complete ww3
! ce qui peut sembler illogique (pourquoi ne pas faire la même chose mais sur
! la zone d'extraction reduite?) mais qui garantit la conservation de la parallelisation
! Attention qu'au moment de l'interpolation, ceci suppose d'ajouter un shift à deci decj
      if(ww3_type_grid==type_unstructured) then !-unstructured-unstructured-unstructured->

        call wave_ww3_grid_resolution ; ww3zoom_istr=1 ; ww3zoom_jstr=1

      endif                                     !-unstructured-unstructured-unstructured->

      if(ww3_type_grid==type_structured)   then !- structured - structured - structured ->

       call wave_sgrid_to_ww3grid !11-04-13


! A partir des min & max des (lon,lat) des 2 grilles, reduire la zone d'extraction:
       ww3zoom_istr=999999 ; ww3zoom_iend=-999999 ! first guess
       ww3zoom_jstr=999999 ; ww3zoom_jend=-999999 ! first guess

       ksecu=0
       do j=1,ww3_jmax
       do i=1,ww3_imax

       if(ww3cg_lon(i,j)>=lonmin.and.ww3cg_lon(i,j)<=lonmax.and.   &
          ww3cg_lat(i,j)>=latmin.and.ww3cg_lat(i,j)<=latmax)then

         ww3zoom_istr=min(ww3zoom_istr,i)
         ww3zoom_jstr=min(ww3zoom_jstr,j)
         ww3zoom_iend=max(ww3zoom_iend,i)
         ww3zoom_jend=max(ww3zoom_jend,j)
         ksecu=1

       endif

       enddo
       enddo

! Si ksecu=0 c'est que le proc est si petit, qu'aucun point ww3 ne s'est trouvé à l'interieur
! On refait le test differement. On cherche le point le plus proche.
       if(ksecu==0)then !000000000000000>
        dist1=1.d20 ; i1=imax/2 ; j1=jmax/2
        do j=1,ww3_jmax
        do i=1,ww3_imax

         dist2=rayonterre*                                         &
         acos( sin(ww3cg_lat(i,j)*deg2rad)*sin(lat_t(i1,j1))     &
              +cos(ww3cg_lat(i,j)*deg2rad)*cos(lat_t(i1,j1))     &
              *cos(lon_t(i1,j1)-ww3cg_lon(i,j)*deg2rad))

         if(dist2<dist1)then !--->
          dist1=dist2 ; i2=i ; j2=j
         endif               !--->

        enddo
        enddo
        ww3zoom_istr=i2 ; ww3zoom_iend=i2
        ww3zoom_jstr=j2 ; ww3zoom_jend=j2
       endif            !000000000000000>


! on elargit un peu plus pour boucher les trous sans être restreint par la taille
! reduite de la zone d'extraction, et puis aussi pour rattraper l'erreur liee au fait
! que plusieurs grilles ww3 (point u, v, t) peuvent être presentes.
       i0=10 ! elargissement à 10 lignes / 10 colonnes supplementaires - repere 1233
       ww3zoom_istr=max0(ww3zoom_istr-i0,1)
       ww3zoom_iend=min0(ww3zoom_iend+i0,ww3_imax)
       ww3zoom_jstr=max0(ww3zoom_jstr-i0,1)
       ww3zoom_jend=min0(ww3zoom_jend+i0,ww3_jmax)

       ww3_imax=ww3zoom_iend-ww3zoom_istr+1
       ww3_jmax=ww3zoom_jend-ww3zoom_jstr+1

      endif                                     !- structured - structured - structured ->

! Desallouer les tableaux
       deallocate(ww3_lon)
       deallocate(ww3_lat)
       deallocate(ww3cg_lon)
       deallocate(ww3cg_lat)

      status=nf_close(ncid_)

      end subroutine wave_grid_extraction_ww3cg

!------------------------------------------------------------------------------

      subroutine wave_ij_to_ijww3
      use module_principal
      use module_parallele !#mpi
      implicit none
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)
      integer loop1_,loop2_
#ifdef synopsis
       subroutinetitle='wave_ij_to_ijww3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')

      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

      x2=real(ww3_imax/2)
      x3=real(ww3_jmax/2)
      deci=x2
      decj=x3

! etape 1: trouver les coordonnées dans la grille orCa:

! first guess: centre du domaine:
      k10=0
 1456 continue

! principe suppose une relation lineaire entre lat lon et indice de grille.
!      dlon/di*di+dlon/dj*dj=dlon
!      dlat/di*di+dlat/dj*dj=dlat
! on cherche di et dj correspondant à dlon=lon(i,j)-lonogcm(i0,j0)
!                                et à dlat=lat(i,j)-latogcm(i0,j0)

      i1=int(deci)
      j1=int(decj)

      do j2=0,1
      do i2=0,1
      i0=i1+i2
      j0=j1+j2

! determinant principal:
      x1=( (ww3cg_lon(i0+1,j0  )-ww3cg_lon(i0-1,j0  ))          &
          *(ww3cg_lat(i0  ,j0+1)-ww3cg_lat(i0  ,j0-1))          &
          -(ww3cg_lat(i0+1,j0  )-ww3cg_lat(i0-1,j0  ))          &
          *(ww3cg_lon(i0  ,j0+1)-ww3cg_lon(i0  ,j0-1)) )*0.25

      deci_(i2,j2)=min(max(                                 &
       i0+((rad2deg*lon_t(i,j)     -ww3cg_lon(i0,j0  ))            &
          *(ww3cg_lat(i0  ,j0+1)-ww3cg_lat(i0,j0-1))           &
          -(rad2deg* lat_t(i,j)     -ww3cg_lat(i0,j0))              &
          *(ww3cg_lon(i0  ,j0+1)-ww3cg_lon(i0,j0-1)))/x1*0.5   &
                   ,2.0001d0),ww3_imax-1.0001d0)

      decj_(i2,j2)=min(max(                                &
       j0+((ww3cg_lon(i0+1,j0  )-ww3cg_lon(i0-1,j0))          &
          *(rad2deg* lat_t(i,j)     -ww3cg_lat(i0  ,j0))           &
          -(ww3cg_lat(i0+1,j0  )-ww3cg_lat(i0-1,j0))          &
          *(rad2deg* lon_t(i,j)     -ww3cg_lon(i0  ,j0)))/x1*0.5   &
                   ,2.0001d0),ww3_jmax-1.0001d0)

      enddo
      enddo

      rapi=deci-i1
      rapj=decj-j1

      deci=(1.-rapi)*(1.-rapj)*deci_(0,0)   &
          +(1.-rapi)*    rapj *deci_(0,1)   &
          +    rapi *    rapj *deci_(1,1)   &
          +    rapi *(1.-rapj)*deci_(1,0)
      decj=(1.-rapi)*(1.-rapj)*decj_(0,0)   &
          +(1.-rapi)*    rapj *decj_(0,1)   &
          +    rapi *    rapj *decj_(1,1)   &
          +    rapi *(1.-rapj)*decj_(1,0)

! si le point visé est different du first guess refaire le calcul
! avec un first guess donné par le dernier point visé:
      if(sqrt( (deci-x2)**2+(decj-x3)**2 ).gt.0.001)then
       x2=deci
       x3=decj
       k10=k10+1
       goto 1456
       if(k10>20)stop 'hr_to_lr ne converge pas dans le forfait'
      endif

! etape 2: CalCuler l'angle d'orientation loCale de la grille ww3

      i1=int(deci) ; j1=int(decj)
      rapi=deci-i1 ; rapj=decj-j1
      do j2=0,1
      do i2=0,1
       i0=i1+i2
       j0=j1+j2
       dy_(i2,j2)= ww3cg_lat(i0+1,j0)-ww3cg_lat(i0-1,j0)
       dx_(i2,j2)=(ww3cg_lon(i0+1,j0)-ww3cg_lon(i0-1,j0))           &
                            *cos(ww3cg_lat(i0,j0)*deg2rad)
      enddo
      enddo
      x1=(1.-rapi)*(1.-rapj)*dx_(0,0)   &
        +(1.-rapi)*    rapj *dx_(0,1)   &
        +    rapi *    rapj *dx_(1,1)   &
        +    rapi *(1.-rapj)*dx_(1,0)
      y1=(1.-rapi)*(1.-rapj)*dy_(0,0)   &
        +(1.-rapi)*    rapj *dy_(0,1)   &
        +    rapi *    rapj *dy_(1,1)   &
        +    rapi *(1.-rapj)*dy_(1,0)

      write(3,'(2(i4,1x),3(f9.4,1x))')i,j,deci,decj,atan2(y1,x1)*rad2deg

      enddo
      enddo

      close(3)

      stop 'NEW HELLO 3 !'
      end subroutine wave_ij_to_ijww3

!------------------------------------------------------------------------------

      subroutine wave_read_interp_netcdf_ww3cg(t_)
      use module_principal ; use module_parallele
      implicit none
      include 'netcdf.inc'
      double precision dlon_di_
      real*8 valmin_
      integer t_
#ifdef synopsis
       subroutinetitle='wave_read_interp_netcdf_ww3cg'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cas du modele ww3
! debut:
!     if(txt_wavemodeltype=='ww3') then
      if(txt_wavemodeltype=='ww3cg') then     !26-08-09         !29-01-13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! trouver le fichier et le numero du champ à l'interieur du fichier
!!    x2=(real(kount  )-wavedt(2))/wavedt(1)+0.5*t_
      x2=(elapsedtime_now-wavedt(2))/wavedt(1)+0.5*t_                !15-04-11
      nc=1+int(x2)/dataperwavefile           ! numero du fichier dans la liste binaire
      nc1=1+int(x2)-(nc-1)*dataperwavefile   ! numero du champs dans le fichier

!.........
! verifier que le numero de record existe bien dans le fichier:        !28-08-09
      texte30=trim(tmpdirname)//''//dom_c//'_listwave_max.txt'
      open(unit=3,file=texte30)
      read(3,*)k
      close(3)
      if(nc>k.or.nc<1)then
       write(*,*)
       write(*,*)'probleme avec le numero de record ',nc
       write(*,*)'il n''y a pas de champs de vagues à la date en cours'
       write(*,*)'conseils:'
       write(*,*)'1-verifier les dates de notebook_time & notebook_wave'
       write(*,*)'2-verifier qu''il y a assez de champs de vagues'
       stop ' en attendant stop dans modele_wave.f'
      endif
!.........

!      call allocate_forcages_ww3cg(1,10,ww3_imax,ww3_jmax,30)
!     call allocate_forcages(1,10,ww3_imax,ww3_jmax,30)
      call wave_allocate_ww3('a',ww3_imax,ww3_jmax,30)

      do loop2=1,13 !06-05-11
      loop1=loop2
      if(loop2==1)loop1=01
      if(loop2==2)loop1=02
      if(loop2==3)loop1=03
      if(loop2==4)loop1=04                               !02-09-10
      if(loop2==5)loop1=04
      if(loop2==6)loop1=05
      if(loop2==7)loop1=05
      if(loop2==8)loop1=06
      if(loop2==9)loop1=06      !06-05-11
      if(loop2==10)loop1=07      !06-05-11
      if(loop2==11)loop1=08      !06-05-11
      if(loop2==12)loop1=09
      if(loop2==13)loop1=10

      write(texte30,'(i2,a7)')loop1,'.binrec' ! nom de la liste binaire
      texte30=trim(tmpdirname)//''//dom_c//'_listewave'//texte30

        open(unit=3,file=texte30              &             ! ouvre la liste binaire
                   ,access='direct'           &
                   ,recl=80*8                 &
                   ,form='unformatted')
        read(3,rec=nc)texte80(1) ! lit le nom du fichier du modele de vague à interpoler
        if(par%rank==0)write(6,'(a)')trim(texte80(1))
        close(3)

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec ouverture fichier ww3'

      texte4='none'
      if(loop2==1)texte4='t'
      if(loop2==2)texte4='hs'
      if(loop2==3)texte4='dir'
      if(loop2==4)texte4='utwo'
      if(loop2==5)texte4='vtwo'
      if(loop2==6)texte4='utaw'
      if(loop2==7)texte4='vtaw'
!      if(loop2==8)texte4='uusf'                     !02-09-10
!      if(loop2==9)texte4='vusf'
      if(loop2==8)texte4='uuss'                     !02-09-10
      if(loop2==9)texte4='vuss'
      if(loop2==10)texte4='l'                       !20-10-10
      if(loop2==11)texte4='foc'   !06-05-11
      if(loop2==12)texte4='wch'   !06-05-11
      if(loop2==13)texte4='J'

      if(texte4=='none')stop 'Variable ww3 non identifiee'
      status=nf_inq_varid(ncid1,texte4,var_id)
!! noter l'erreur dans les fichiers ww3 qui oblige eventuellement à fixer texte4='hs' !06-05-11
!      if(status/=0.and.trim(texte4)=='l') then
!       texte4='hs'
!       status=nf_inq_varid(ncid1,texte4,var_id)
!      endif
      if(status.ne.0)then
       write(*,'(a,a)')'echec id variable ww3 ',trim(texte4)
       stop 'stop in subroutine wave'
      endif

      if(par%rank==0)write(6,'(a)')trim(texte80(1))
      status=nf_get_att_real(ncid1,var_id,'_fillValue',filval)
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'_FillValue',filval)
      if(status/=0)stop 'echec get _fillValue fichier ww3'

      status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
      if(status/=0)stop 'echec get scale_factor fichier ww3'

      status=nf_get_att_real(ncid1,var_id,'offset',var_addoffset)
      if(status/=0)stop 'echec get add_offset fichier ww3'

      varstart(3)=nc1          ; varcount(3)=1        ! time
      varstart(2)=ww3zoom_jstr ; varcount(2)=ww3_jmax ! lat
      varstart(1)=ww3zoom_istr ; varcount(1)=ww3_imax ! lon

      status=nf_inq_varid(ncid1,texte4,var_id)
      if(status.ne.0)stop 'echec id variable ww3 2'
     status=nf_get_vara_int(ncid1,var_id,varstart(1:3),varcount(1:3),ww3_var)
      if(status.ne.0)stop 'echec lecture variable netcdf ww3 1'
!      endif                              !*-*-*-*-*-*-*-*->


!..........................
! interpoler la variable 1:
      if(loop2==1) then !1-1-1-1-1-1-1->
      if(texte4/='t')stop 'erreur wave_read_interp t'

      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')

       valmin_=0.001 ; ww3mskval=1.
       do j=0,jmax+1
       do i=0,imax+1

       read(3,*)i3,j3,deci,decj,x0

       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1                              ! conservation de la parallelisation

! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1                              ! conservation de la parallelisation


! a la premiere variable on calcule (et on stocke pour les variables suivantes)
! les coefficient de la ponderation:
      x1=(1.-rapi)*(1.-rapj)
      x2=(1.-rapi)*    rapj
      x3=    rapi *(1.-rapj)
      x4=    rapi *    rapj
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor
      x0=x1+x2+x3+x4+small1
      xy_t(i,j,1)=x1/x0
      xy_t(i,j,2)=x2/x0
      xy_t(i,j,3)=x3/x0
      xy_t(i,j,4)=x4/x0
! interpolation bilineaire:
        t_wave_t(i,j,t_)=max(var_scalefactor*                          &
                              (xy_t(i,j,1)*ww3_var(i1  ,j1  )          &
                              +xy_t(i,j,2)*ww3_var(i1  ,j1+1)          &
                              +xy_t(i,j,3)*ww3_var(i1+1,j1  )          &
                              +xy_t(i,j,4)*ww3_var(i1+1,j1+1)),valmin_)
! note: pour eviter les divisions par zero dans les zones masquées on a appliqué
! un seuil sur t
       else                         !>------------->
        t_wave_t(i,j,t_)=ww3mskval
       endif                        !>------------->
       enddo
       enddo
       close(3)

! a l'initialisation du modele (c.a.d. t_==0) on verifie la coherence
! des grilles du modele de circulation et du modele de vague. si le
! recoupement est mauvais, on propose un reajustement de la grille soumis
! à l'approbation du l'utilisateur:
      if(t_==0)call wave_inq_grid_consistency    !20-10-10
      endif             !1-1-1-1-1-1-1->


!..........................
! interpoler la variable 2:
      if(loop2==2) then !2-2-2-2-2-2-2->
      if(texte4/='hs')stop 'erreur wave_read_interp hs'
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')

       ww3mskval=0.01
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1                              ! conservation de la parallelisation

! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1                              ! conservation de la parallelisation

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor

! pour hs on utilise rapi et rapj plutôt que xy_t afin d'avoir une C.l. en terre
! de type: variable=0
! interpolation bilineaire:
        hs_wave_t(i,j,t_)=var_scalefactor*                     &
                          ((1.-rapi)*(1.-rapj)*ww3_var(i1  ,j1  ) &
                          +(1.-rapi)*    rapj *ww3_var(i1  ,j1+1) &
                          +    rapi *(1.-rapj)*ww3_var(i1+1,j1  ) &
                          +    rapi *    rapj *ww3_var(i1+1,j1+1))
! note: pour hs on a appliqué une condition limite hs=0 aux frontieres solides
! en incluant des valeurs nulles (masque) dans l'interpolation

       else                         !>------------->
        hs_wave_t(i,j,t_)=ww3mskval
       endif                        !>------------->
       enddo
       enddo
       close(3)
      if(t_==0) call wave_inq_grid_consistency
      endif             !2-2-2-2-2-2-2->

!..........................
! interpoler la variable 3:
      if(loop2==3) then !3-3-3-3-3-3-3->
      if(texte4/='dir')stop 'erreur wave_read_interp 3'
! attention en ce qui concerne les conventions d'angle en meteorologie. celles ci
! sont appliquées par ww3. en meteo, un vent venant du nord a une direction de 0°,
! un vent venant de l'est a une direction de 90°. par consequent il faut corriger
! cette convention en faisant 270°-angle. ensuite on convertit en radians.
! d'autre part on n'interpole pas un angle mais ses composantes trigonometriques.
! puis on retrouve un angle avec la fonction atan2
! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p

! application du facteur d'echelle, prise en compte des conventions meteo ,
! conversion degres->radians, prise en compte de la rotation de la grille
! symphonie:
! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p
       x1=       270.*deg2rad ! convention meteo
       x2=var_scalefactor*deg2rad
!      do j=1,ww3_jmax
!      do i=1,ww3_imax
!       ww3_var(i,j)=x1-x2*ww3_var(i,j)
!      enddo
!      enddo
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
!      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
!      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.
        dir_wave_t(i,j,1)=atan2(  & ! retrouver l'angle à partir de sinus et cosinus
! interpolation de la composante sinus:
         xy_t(i,j,1)*sin(x1-x2*ww3_var(i1  ,j1  ))  &
        +xy_t(i,j,2)*sin(x1-x2*ww3_var(i1  ,j1+1))  &
        +xy_t(i,j,3)*sin(x1-x2*ww3_var(i1+1,j1  ))  &
        +xy_t(i,j,4)*sin(x1-x2*ww3_var(i1+1,j1+1))  &
       ,                                            &
! interpolation de la composante cosinus:
         xy_t(i,j,1)*cos(x1-x2*ww3_var(i1  ,j1  ))  &
        +xy_t(i,j,2)*cos(x1-x2*ww3_var(i1  ,j1+1))  &
        +xy_t(i,j,3)*cos(x1-x2*ww3_var(i1+1,j1  ))  &
        +xy_t(i,j,4)*cos(x1-x2*ww3_var(i1+1,j1+1)) )

! interpolation bi-polynome ordre 3:                                   !29-06-09
!       call wave_poly3(1) ! interpoler cosinus
!       x10=x1
!       call wave_poly3(2) ! interpoler sinus
!       dir_wave_t(i,j,1)=atan2(x1,x10) ! atan2(sin,cos) retrouve l'angle à partir de sin et cos

! enfin prendre en compte l'orientation locale de la grille symphonie:
        dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)       !14-02-13
        if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
        if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
        dir_wave_t(i,j,1)=dir_wave_t(i,j,1)                   &
        -atan2(  lat_t(i+1,j)-lat_t(i-1,j)                    & ! rotation symphonie
               ,dlon_di_*cos(lat_t(i,j)) )                      !14-02-13
!              ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j)) )

       else                         !>------------->
        dir_wave_t(i,j,1)=0.
       endif                        !>------------->
       enddo
       enddo
       close(3)
       if(t_==0)call wave_inq_grid_consistency
      endif             !3-3-3-3-3-3-3->


!..........................
! interpoler la variable 4:
      if(loop2==4) then !4-4-4-4-4-4-4->
      if(texte4/='utwo')stop 'erreur wave_read_interp 4'

! interpolation on ocean model "t" grid nodes:
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1


      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        twox_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                        !>------------->
       enddo
       enddo
       close(3)
       if(t_==0)call wave_inq_grid_consistency
      endif             !4-4-4-4-4-4-4->

!..........................
! interpoler la variable 5:
      if(loop2==5) then !5-5-5-5-5-5-5->
      if(texte4/='vtwo')stop 'erreur wave_read_interp 5'

! interpolation on ocean model "t" grid nodes:
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        twoy_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                      !>------------->
       enddo
       enddo
       close(3)
       if(t_==0)call wave_inq_grid_consistency
      endif             !5-5-5-5-5-5-5->

!..........................
! interpoler la variable 6:
      if(loop2==6) then !6-6-6-6-6-6-6->
      if(texte4/='utaw')stop 'erreur wave_read_interp utaw'

      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        tawx_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))



       endif                      !>------------->
       enddo
       enddo
       close(3)
       if(t_==0)call wave_inq_grid_consistency
      endif             !6-6-6-6-6-6-6->

!..........................
! interpoler la variable 7:
      if(loop2==7) then !7-7-7-7-7-7-7->
      if(texte4/='vtaw')stop 'erreur wave_read_interp vtaw'

! interpolation on ocean model "t" grid nodes:
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        tawy_wave_t(i,j,t_)=rho*var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                      !>------------->
       enddo
       enddo
       close(3)
       if(t_==0)call wave_inq_grid_consistency
      endif             !7-7-7-7-7-7-7->


!      if(loop2==8) then !8-8-8-8-8-8-8-8-8-8->
!      if(texte4/='uusf')stop 'erreur wave_read_interp usf'
!       do freq=1,ww3_fmax
!       do j=0,jmax+1
!       do i=0,imax+1
!       if(mask_t(i,j,kmax+1)==1) then !>------------->
!
!! indice decimale i (longitude) dans grille aladin:
!      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1
!
!! indice decimale j (latitude) dans grille aladin:
!      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
!
!      if(ww3_var2(i1  ,j1  ,freq)==nint(filval))ww3_var2(i1  ,j1  ,freq)=small1
!      if(ww3_var2(i1  ,j1+1,freq)==nint(filval))ww3_var2(i1  ,j1+1,freq)=small1
!      if(ww3_var2(i1+1,j1  ,freq)==nint(filval))ww3_var2(i1+1,j1  ,freq)=small1
!      if(ww3_var2(i1+1,j1+1,freq)==nint(filval))ww3_var2(i1+1,j1+1,freq)=small1
!              usf_wave_t(i,j,freq)=var_scalefactor*                    &
!                                 (xy_t(i,j,1)*ww3_var2(i1  ,j1  ,freq) &
!                                 +xy_t(i,j,2)*ww3_var2(i1  ,j1+1,freq) &
!                                 +xy_t(i,j,3)*ww3_var2(i1+1,j1  ,freq) &
!                                 +xy_t(i,j,4)*ww3_var2(i1+1,j1+1,freq))
!
!
!
!       endif                      !>------------->
!       enddo
!       enddo
!       enddo
!       if(t_==0) call wave_inq_grid_consistency
!      endif             !8-8-8-8-8-8-8-8-8->


!..........................
! interpoler la variable 8:
      if(loop2==8) then !8-8-8-8-8-8-8-8-8-8->
      if(texte4/='uuss')stop 'erreur wave_read_interp uss'
        open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
        do j=0,jmax+1
        do i=0,imax+1
        read(3,*)i3,j3,deci,decj,x0
        if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=small1
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=small1
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=small1
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=small1
         uss_wave_t(i,j,1)=var_scalefactor*                    &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))



        endif                      !>------------->
        enddo
        enddo
        close(3)
       if(t_==0)call wave_inq_grid_consistency
      endif             !8-8-8-8-8-8-8-8-8->


!      if(loop2==9) then !9-9-9-9-9-9-9-9-9->
!      if(texte4/='vusf')stop 'erreur wave_read_interp vsf'
!
!! interpolation on ocean model "t" grid nodes:
!       do freq=1,ww3_fmax
!       do j=0,jmax+1
!       do i=0,imax+1
!       if(mask_t(i,j,kmax+1)==1) then !>------------->
!
!! indice decimale i (longitude) dans grille aladin:
!      i1=int(1.+(lon_t(i,j)*rad2deg-ww3_lonmin)/ww3_dlon)-ww3zoom_istr+1
!
!! indice decimale j (latitude) dans grille aladin:
!      j1=int(1.+(lat_t(i,j)*rad2deg-ww3_latmin)/ww3_dlat)-ww3zoom_jstr+1
!      if(ww3_var2(i1  ,j1  ,freq)==nint(filval))ww3_var2(i1  ,j1  ,freq)=small1
!      if(ww3_var2(i1  ,j1+1,freq)==nint(filval))ww3_var2(i1  ,j1+1,freq)=small1
!      if(ww3_var2(i1+1,j1  ,freq)==nint(filval))ww3_var2(i1+1,j1  ,freq)=small1
!      if(ww3_var2(i1+1,j1+1,freq)==nint(filval))ww3_var2(i1+1,j1+1,freq)=small1
!              vsf_wave_t(i,j,freq)=var_scalefactor*                    &
!                                 (xy_t(i,j,1)*ww3_var2(i1  ,j1  ,freq) &
!                                 +xy_t(i,j,2)*ww3_var2(i1  ,j1+1,freq) &
!                                 +xy_t(i,j,3)*ww3_var2(i1+1,j1  ,freq) &
!                                 +xy_t(i,j,4)*ww3_var2(i1+1,j1+1,freq))
!
!       endif                      !>------------->
!       enddo
!       enddo
!       enddo
!        if(t_==0     call wave_inq_grid_consistency
!     endif             !9-9-9-9-9-9-9-9-9-9->

!..........................
! interpoler la variable 9:
      if(loop2==9) then !9-9-9-9-9-9-9-9-9->
      if(texte4/='vuss')stop 'erreur wave_read_interp vss'

! interpolation on ocean model "t" grid nodes:
        open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
        do j=0,jmax+1
        do i=0,imax+1
        read(3,*)i3,j3,deci,decj,x0
        if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1
      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=small1
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=small1
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=small1
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=small1
              vss_wave_t(i,j,1)=var_scalefactor*                    &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

        endif                      !>------------->
        enddo
        enddo
        close(3)
      if(t_==0) call wave_inq_grid_consistency
      endif             !9-9-9-9-9-9-9-9-9-9->


!..........................
! interpoler la variable 10:
      if(loop2==10) then !10-10-10-10-10->
!     if(texte4/='l')stop 'erreur wave_read_interp l'      !20-10-10
! noter prise en compte erreur attribut dans fichier ww3: !06-05-11
      if(texte4/='l'.and.texte4/='hs') then
        write(*,*)'erreur nom variable ww3: ',trim(texte4)
        stop 'erreur in subroutine wave'
      endif

       const2=2.*pi
       ww3mskval=10.
       open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor

                                x1=var_scalefactor*              &
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

        k_wave_t(i,j,t_)=const2/x1
        kx_wave_t(i,j,t_)=k_wave_t(i,j,t_)*cos(dir_wave_t(i,j,1))
        ky_wave_t(i,j,t_)=k_wave_t(i,j,t_)*sin(dir_wave_t(i,j,1))
      else                         !>------------->
        k_wave_t(i,j,t_)=const2/ww3mskval
        kx_wave_t(i,j,t_)=k_wave_t(i,j,t_)*cos(dir_wave_t(i,j,1))
        ky_wave_t(i,j,t_)=k_wave_t(i,j,t_)*sin(dir_wave_t(i,j,1))
      endif                        !>------------->
      enddo
      enddo
      close(3)
      if(t_==0)call wave_inq_grid_consistency
      endif      !10-10-10-10-10->

!..........................
! interpoler la variable 11:
      if(loop2==11) then !11-11-11-11-11-11->
      if(texte4/='foc')stop 'erreur wave_read_interp foc'

! interpolation on ocean model "t" grid nodes:
       open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1
! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

        foc_wave_t(i,j,t_)=var_scalefactor*                   &     !07-05-11
                                 (xy_t(i,j,1)*ww3_var(i1  ,j1  ) &
                                 +xy_t(i,j,2)*ww3_var(i1  ,j1+1) &
                                 +xy_t(i,j,3)*ww3_var(i1+1,j1  ) &
                                 +xy_t(i,j,4)*ww3_var(i1+1,j1+1))

       endif                        !>------------->
       enddo
       enddo
       close(3)
       if(t_==0)   call wave_inq_grid_consistency
      endif             !11-111-11-11-11-11-11-->

!..........................
! interpoler la variable 12:
      if(loop2==12) then !12-12-12-12-12-12->
      if(texte4/='wch')stop 'erreur wave_read_interp wch'

      ww3mskval=0.01
      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1                              ! conservation de la parallelisation

! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1                              ! conservation de la parallelisation

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=ww3mskval/var_scalefactor
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=ww3mskval/var_scalefactor

! pour hsw on utilise rapi et rapj plutôt que xy_t afin d'avoir une C.l. en terre
! de type: variable=0
! interpolation bilineaire:
        hsw_wave_t(i,j,t_)=var_scalefactor*                    &
                          ((1.-rapi)*(1.-rapj)*ww3_var(i1  ,j1  ) &
                          +(1.-rapi)*    rapj *ww3_var(i1  ,j1+1) &
                          +    rapi *(1.-rapj)*ww3_var(i1+1,j1  ) &
                          +    rapi *    rapj *ww3_var(i1+1,j1+1))
! note: pour hs on a appliqué une condition limite hs=0 aux frontieres solides
! en incluant des valeurs nulles (masque) dans l'interpolation

       else                         !>------------->
        hsw_wave_t(i,j,t_)=ww3mskval
       endif                        !>------------->
       enddo
       enddo
       close(3)
      if(t_==0) call wave_inq_grid_consistency
      endif             !12-12-12-12-12-12--->

! interpoler la variable 13:
      if(loop2==13) then !13-13-13-13-13-->
      if(texte4/='J')stop 'erreur wave_read_interp J'

      open(unit=3,file=trim(tmpdirname)//''//dom_c//'_ij2ijww3.txt')
       do j=0,jmax+1
       do i=0,imax+1
       read(3,*)i3,j3,deci,decj,x0
       if(mask_t(i,j,kmax+1)==1) then !>------------->

! indice decimale i (longitude) dans grille aladin:
      i1=int(deci)
      rapi=deci-i1
      i1=i1-ww3zoom_istr+1                              ! conservation de la parallelisation

! indice decimale j (latitude) dans grille aladin:
      j1=int(decj)
      rapj=decj-j1
      j1=j1-ww3zoom_jstr+1                              ! conservation de la parallelisation

      if(ww3_var(i1  ,j1  )==nint(filval))ww3_var(i1  ,j1  )=0.
      if(ww3_var(i1  ,j1+1)==nint(filval))ww3_var(i1  ,j1+1)=0.
      if(ww3_var(i1+1,j1  )==nint(filval))ww3_var(i1+1,j1  )=0.
      if(ww3_var(i1+1,j1+1)==nint(filval))ww3_var(i1+1,j1+1)=0.

! interpolation bilineaire:
        j_wave_t(i,j,t_)=var_scalefactor*                    &
                          ((1.-rapi)*(1.-rapj)*ww3_var(i1  ,j1  ) &
                          +(1.-rapi)*    rapj *ww3_var(i1  ,j1+1) &
                          +    rapi *(1.-rapj)*ww3_var(i1+1,j1  ) &
                          +    rapi *    rapj *ww3_var(i1+1,j1+1))

       endif                        !>------------->
       enddo
       enddo
       close(3)
      if(t_==0) call wave_inq_grid_consistency
      endif             !13-13-13-13-13-13->


! fermer le fichier netcdf ww3
      status=nf_close(ncid1)
      if(status/=0)stop 'erreur fermeture fichier ww3'
      enddo ! fin de boucle du loop1
!     call allocate_forcages(2,10,0,0,0)
      call wave_allocate_ww3('d',0,0,0)

!.........................................................

!      stop 'along axis or along ns-we?'

! grid rotation:
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !////////>

        x1=tawx_wave_t(i,j,t_)
        x2=tawy_wave_t(i,j,t_)
        tawx_wave_t(i,j,t_)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
        tawy_wave_t(i,j,t_)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2

        x1=twox_wave_t(i,j,t_)
        x2=twoy_wave_t(i,j,t_)
        twox_wave_t(i,j,t_)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
        twoy_wave_t(i,j,t_)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2
      endif                        !////////>
      enddo
      enddo

      do j=0,jmax+1
      do i=0,imax+1

       x1=uss_wave_t(i,j,1)
       x2=vss_wave_t(i,j,1)
       uss_wave_t(i,j,1)=gridrotcos_t(i,j)*x1-gridrotsin_t(i,j)*x2
       vss_wave_t(i,j,1)=gridrotsin_t(i,j)*x1+gridrotcos_t(i,j)*x2
      enddo
      enddo
!.........................................................


!! Compute the total 3d stokes velocity using the frequential components
!! and their respective z profil
!      do k=1,kmax
!      do j=0,jmax+1
!      do i=0,imax+1
!       anyv3d(i,j,k,1)=0.
!       anyv3d(i,j,k,2)=0.
!      enddo
!      enddo
!      enddo

      const1=0.5*grav/pi
      const2=2.*pi
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !mmmmmmmmmmmmmm>

!      do freq=1,ww3_fmax

!! x1 is the wave vector kn(freq)
!! x0 is the wave celerity at the frequency freq_wave(freq)
!! on veut résoudre x=f(x)=alpha*tanh(beta/x)
!! on dit que f(x+dx)=x+dx
!!            f(x)+dx f'(x)=x+dx
!!            dx=(x-f(x))/(f'(x)-1)
!! on estime une première valeur de x=x0
!      x3=const2*freq_wave(freq)*max(hz_w(i,j,1),0.01d0)         !02-02-11 beta
!      x0=const1/freq_wave(freq) ! c=g*t/2pi first guess for large depth
!
!      if(x3/x0>200.) then !fgfgfgfgfgfgfg>      !07-05-11
!! Cas c=first guess valide:
!       x1=x0
!
!      else                !fgfgfgfgfgfgfg>
!
!! Cas calcul itératif:
!       x12=0.01                  ! dx, on estime une première valeur de dx
!       x1=x0+x12                 ! x0+dx
!         do k0=1,200
!          x12=(x1-x0*tanh(x3/x1))/(-x0*x3/(x1*cosh(x3/x1))**2-1.)
!          x1=max(x1+x12,1.d-4)
!          if(abs(x12)<0.01)goto 25
!         enddo
!         write(*,*)i,j,x1,x12
!       stop 'wave: le calcul de kn ne converge pas'
!   25  continue
!
!      endif               !fgfgfgfgfgfgfg>
! alternative si on ne dispose pas des vitesses de stokes sur les 30 fréquences:
       x1=sqrt(kx_wave_t(i,j,t_)**2+ky_wave_t(i,j,t_)**2+small1)

       if((h_w(i,j)+depth_w(i,j,kmax+1))*x1<6.)  then !hhhhhhhhhhhhhhhhh>

! Cas "petites profondeurs":
       if(h_w(i,j)+depth_w(i,j,kmax+1)>0.01) then !ssssssssssss> ! si dz=0 us=0
       x4=cosh(2.*x1*(h_w(i,j)+depth_w(i,j,kmax+1))) ! x4=cosh(2k(h+ssh))
       do k=kmin_w(i,j),kmax

! méthode 1: valeur ponctuelle du profil:
!       x3=cosh(2.*x1*(depth_t(i,j,k)+h_w(i,j)))/x4
! méthode 2: profil moyénné sur dz
        x3=0.5/x1/x4/(depth_w(i,j,k+1)-depth_w(i,j,k))           &
           *( sinh(2.*x1*(depth_w(i,j,k+1)+h_w(i,j)))            &
             -sinh(2.*x1*(depth_w(i,j,k  )+h_w(i,j))))
!       x3=1. ! cas academique profil constant
!        anyv3d(i,j,k,1)=anyv3d(i,j,k,1)+usf_wave_t(i,j,freq)*x3
!        anyv3d(i,j,k,2)=anyv3d(i,j,k,2)+vsf_wave_t(i,j,freq)*x3
        anyv3d(i,j,k,1)=uss_wave_t(i,j,1)*x3
        anyv3d(i,j,k,2)=vss_wave_t(i,j,1)*x3

       enddo ! fin de boucle sur k
       endif                                      !ssssssssssss>

       else                                          !hhhhhhhhhhhhhhhhh>

! Cas "grandes profondeurs":
       do k=kmin_w(i,j),kmax

! méthode 1: valeur ponctuelle du profil exponentiel
!       x3=exp(2.*x1*(depth_t(i,j,k)-ssh_int_w(i,j,1)))
! méthode 2: profil exponentiel moyenné sur dz:
        x3=0.5/x1/(depth_w(i,j,k+1)-depth_w(i,j,k))                    &
                    *( exp(2.*x1*(depth_w(i,j,k+1)-ssh_int_w(i,j,1)))  &
                      -exp(2.*x1*(depth_w(i,j,k  )-ssh_int_w(i,j,1))))
!       x3=1. ! cas academique profil constant
!        anyv3d(i,j,k,1)=anyv3d(i,j,k,1)+usf_wave_t(i,j,freq)*x3
!        anyv3d(i,j,k,2)=anyv3d(i,j,k,2)+vsf_wave_t(i,j,freq)*x3
        anyv3d(i,j,k,1)=uss_wave_t(i,j,1)*x3
        anyv3d(i,j,k,2)=vss_wave_t(i,j,1)*x3

       enddo ! fin de boucle sur k

       endif                                         !hhhhhhhhhhhhhhhhh>

!      enddo ! fin de boucle sur freq

      endif                          !mmmmmmmmmmmmmm>
      enddo
      enddo


      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
       velstokes_u(i,j,k,2)=                           &
              dz_u(i,j,k,1)*                           &
            mask_u(i,j,k)*0.5*( anyv3d(i  ,j  ,k,1)    &
                               +anyv3d(i-1,j  ,k,1))
       velstokes_v(i,j,k,2)=                           &
              dz_v(i,j,k,1)*                           &
            mask_v(i,j,k)*0.5*( anyv3d(i  ,j  ,k,2)    &
                               +anyv3d(i  ,j-1,k,2))
      enddo
      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cas du modele ww3
! fin.
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end subroutine wave_read_interp_netcdf_ww3cg
!..................................................
      subroutine wave_initial
      implicit none
      integer year_,month_,day_,hour_,minute_,second_,loop_,ncid_   &
             ,debugi4_,linemax_,loop2_,loop1_
      double precision :: debugr8_,tm1_=0.,tm2_=0.
      namelist/notebook_wave2/texte80,txt_wavemodeltype,wave_obc_type &
                             ,ww3_flagcumul &
                             ,flag_wavelength !16-02-23
#ifdef synopsis
       subroutinetitle='wave_initial'
       subroutinedescription='Reads notebook_wave'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(iwve==0)return

! notebook_wave (Part II)
      txt_wavemodeltype='none' ; wave_obc_type=1 ; ww3_flagcumul=0 ; texte80='none'
      open(100,file=nomfichier(22)) ! notebook_wave
      read(100,nml=notebook_wave2)
      close(100)
      if(iwve==1)ww3_varmax=8  !-----Michaud2012----> !21-11-19
      if(iwve==2)ww3_varmax=3  !-----McWilliams1999----> !11-11-21
      if(iwve==3)ww3_varmax=3  !-----McWilliams1999----> !21-11-19

!     if(iwve==2)flag_wavelength=1 !Plus tard etendre flag_wavelength=1 aux autres cas que iwve=2
      if(iwve==1)flag_wavelength=0 !Plus tard etendre flag_wavelength=1 au cas iwve=1 !16-02-23


! tmp/wavelist est un ficher txt contenant les noms des listes acces direct
      if(par%rank==0)open(unit=10,file=trim(tmpdirname)//'wavelist')

      do loop1_=1,1 !ww3_varmax! depuis 21-11-19 liste unique

        open(unit=3,file=texte80(loop1_)) ! liste texte

!#ifdef bidon

! Combien y a t'il de lignes dans le fichier?
         linemax_=0
 1834    read(3,'(a)',end=1832)texte250 ; linemax_=linemax_+1 ; goto 1834
 1832    rewind(3)
         if(linemax_==0) then !>>>>> !23-04-14
          write(6,'(a,a,a)')'File ',trim(texte80(loop1_)),' is empty'
          stop ' stop module_wave erreur 1928'
         endif

! Chaque proc va traiter une fraction du fichier. On commence par derouler les lignes
! qui ne concernent pas le proc
         do loop2_=1,int(real(par%rank  )/real(nbdom)*linemax_)
          read(3,*) ! lire ces lignes pour rien
         enddo

         write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',par%rank
         open(unit=4,file=texte60,status='REPLACE') !20-10-14

         do loop2_=int(real(par%rank  )/real(nbdom)*linemax_)+1    &
                  ,int(real(par%rank+1)/real(nbdom)*linemax_)

         read(3,'(a)',end=3519)texte250 ! nom du fichier netcdf
         status=nf_open(trim(texte250),nf_nowrite,ncid_)
         if(status/=0) then
          write(6,'(a,a)')'Echec open fichier: ',trim(texte250)
          stop 'erreur module_wave 3521 nf_open'
         endif
         status=nf_inq_dimid(ncid_,'time',var_id)            ;if(status/=0)stop 'erreur nf_inq_dimid time wave_initial'
         status=nf_inq_dimlen(ncid_,var_id,max_time_counter) ;if(status/=0)stop 'erreur nf_inq_dimlen wave_initial'

         status=nf_inq_varid(ncid_,'time',var_id)            ;if(status/=0)stop 'erreur nf_inq_varid time wave_initial'
         txt_units=''
         status=nf_get_att_text(ncid_,var_id,'units',txt_units);if(status/=0)stop 'erreur nf_get_att_text wave_initial'

         k=1
         do while(txt_units(k:k)/='1'.and.txt_units(k:k)/='2') !02-09-16
          k=k+1
         enddo
         read(txt_units(k:k+3),*)year_
         read(txt_units(k+5:k+6),*)month_
         read(txt_units(k+8:k+9),*)day_
         read(txt_units(k+11:k+12),*)hour_
         read(txt_units(k+14:k+15),*)minute_
         read(txt_units(k+17:k+18),*)second_

         call datetokount(year_,month_,day_,hour_,minute_,second_) ! donne elapsedtime_out
         do loop_=1,max_time_counter
           status=nf_get_vara_double(ncid_,var_id,loop_,1,x1);if(status/=0)stop 'erreur nf_get_vara_double wave_initial'
           x2=-999.
           if(index(txt_units,'days')/=0)x2=86400.
           if(index(txt_units,'hours')/=0)x2=3600.
           if(index(txt_units,'seconds')/=0)x2=1.
           if(x2==-999.)stop 'wave unites de temps incomprises'
           x1=elapsedtime_out+x1*x2
!             if(x1<=elapsedtime_now) then !>>>>
!              wavefile_nextrec=nc ; wavefile_nextime=x1
!             endif                        !>>>>
           write(4,'(a)')trim(texte250)
           write(4,*)loop_,x1
         enddo ! fin de boucle sur loop_

         status=nf_close(ncid_)

         enddo ! loop2_
         close(4)

!#endif
 3519 close(3)

! C'est le proc zero qui a la mission d'assembler tous les fichiers tmp en un seul fichier liste binaire
! La barriere suivante permet de s'assurer que tous les fichiers individuels ont bien ete fait au moment
! d'attaquer la concatenation
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
      if(par%rank==0) then !00000000000>

      k1=index(texte80(loop1_),'/',back=.true.)
       k=index(texte80(loop1_),' ')

      texte90=trim(tmpdirname)//''//texte80(loop1_)(k1+1:k-1)//'.binrec'
         write(10,'(a)')trim(texte90)
       open(unit=3,file=trim(texte90) &
                  ,access='direct'    &
                  ,recl=reclen        &
                  ,form='unformatted')
        nc=1
        do loop2_=0,nbdom-1
         write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',loop2_
         open(unit=4,file=texte60)
 2016          tm1_=tm2_
               read(4,'(a)',end=2002)texte250
               read(4,*,end=2002)loop_,tm2_
               if(nc==1)tm1_=tm2_ !03-03-15

! Detection de listes desordonnees: !27-02-15 !03-03-15
               if(tm2_<tm1_) then   !dbdbdbdb>
                write(6,'(a,a)')'List: ',trim(texte90)
                write(6,'(a,a,a)')'File ',trim(texte250) &
                ,' does not respect the chronological order'
                write(6,*)'tm1_,tm2_ ',tm1_,tm2_
                stop ' Stop module_wave'
               endif                !dbdbdbdb>

                if(ww3_flagcumul==1)wavefile_nextime=0.5*(tm1_+tm2_)
                if(ww3_flagcumul==0)wavefile_nextime=tm2_

                write(3,rec=nc)texte250,loop_,wavefile_nextime

               if(ww3_flagcumul==1)then !>>>>
                if(nc==2.and.1.5*tm1_-0.5*tm2_<=elapsedtime_now)wavefile_nextrec=nc-1
               endif                    !>>>>
               if(wavefile_nextime<=elapsedtime_now) then !----->
                wavefile_nextrec=nc
               endif                                        !----->

               if(ww3_flagcumul==1)then !>>>>
               if(nc==2) then !222>
                 read(3,rec=1)texte250,i0,x0
                write(3,rec=1)texte250,i0,1.5*tm1_-0.5*tm2_
               endif          !222>
               endif                    !>>>>

               nc=nc+1
               goto 2016
 2002    close(4)
         enddo ! loop2_
        close(3)

       if(loop1_==1) then !----->
        debugr8_=wavefile_nextime
        debugi4_=wavefile_nextrec
       endif              !----->
       if(wavefile_nextime/=debugr8_.or.wavefile_nextrec/=debugi4_) then
         write(6,*)'loop1_=',loop1_
         write(6,*)'wavefile_nextime debugr8_',wavefile_nextime,debugr8_
         write(6,*)'wavefile_nextrec debugi4_',wavefile_nextrec,debugi4_
         stop 'routine wave_initial files do not start at the same time'
       endif

      endif                !00000000000>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      enddo ! loop1_

      if(par%rank==0)close(10)

#ifdef parallele
      call mpi_bcast(wavefile_nextime,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(wavefile_nextrec,1,mpi_integer         ,0,par%comm2d,ierr)
#endif

!. . . . . . . . .
! début  - definir la zone d'extraction des données -
!       loop1_=1
!       nc=1
!       write(texte30,'(i0,a7)')loop1_,'.binrec' ! nom de la liste binaire
!       texte30=trim(tmpdirname)//'listewave'//texte30
      open(unit=3,file=trim(tmpdirname)//'wavelist') ! liste des listes binaires
      read(3,'(a)')texte30
      close(3)

!     write(6,*)'wavefile_nextrec=',wavefile_nextrec

       nc=1
       open(unit=4,file=texte30                                    &
                  ,access='direct'                                 &
                  ,recl=reclen                                     &
                  ,form='unformatted')
       read(4,rec=nc)texte80(1)
       close(4)
!      write(6,'(a)')trim(texte80(1))

       if(par%rank==0) then !pmxpmx>
        write(6,*)'wavefile_nextime=',wavefile_nextime
        write(6,*)'wavefile_nextrec=',wavefile_nextrec
       endif                !pmxpmx>

! find the min and max (i,j) indexes in the ww3 grid:
!     if(txt_wavemodeltype=='ww3')         then;call wave_grid_extraction_ww3  ;return;endif
!     if(txt_wavemodeltype=='ww3cg')       then;call wave_grid_extraction_ww3cg;return;endif
      if(txt_wavemodeltype=='ww3_native_grid') then;call wave_grid_extraction_ww3cg;return;endif !25-07-22
      if(txt_wavemodeltype=='ww3_on_sgrid')    then;call wave_grid_extraction_wgrid;return;endif

      stop 'txt_wavemodeltype not recognized in module_wave'

      end subroutine wave_initial
!..................................................
      subroutine wave_allocate_ww3cg(txt_,dim1_,dim2_,dim3_)
      use module_principal
      implicit none
      character*1 txt_
      integer dim1_,dim2_,dim3_
#ifdef synopsis
       subroutinetitle='wave_allocate_ww3cg'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       if(txt_=='a') then !aaaaaaaaaaaa>
        allocate(ww3cg_lon(dim1_,dim2_))
        allocate(ww3cg_lat(dim1_,dim2_))
        allocate(ww3_var  (dim1_,dim2_))
!       allocate(ww3_var2 (dim1_,dim2_,dim3_))
       endif              !aaaaaaaaaaaa>
       if(txt_=='d') then !dddddddddddd>
        deallocate(ww3cg_lon)
        deallocate(ww3cg_lat)
        deallocate(ww3_var)
!       deallocate(ww3_var2)
       endif              !dddddddddddd>
      end subroutine wave_allocate_ww3cg
!..................................................
      subroutine wave_allocate_ww3(txt_,dim1_,dim2_,dim3_)
      use module_principal
      implicit none
      character*1 txt_
      integer dim1_,dim2_,dim3_
#ifdef synopsis
       subroutinetitle='wave_allocate_ww3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       if(txt_=='a') then !aaaaaaaaaaaa>
        allocate(ww3_lon (dim1_))
        allocate(ww3_lat (dim2_))
        allocate(ww3_var (dim1_,dim2_))
!       allocate(ww3_var2(dim1_,dim2_,dim3_))
       endif              !aaaaaaaaaaaa>
       if(txt_=='d') then !dddddddddddd>
        deallocate(ww3_lon)
        deallocate(ww3_lat)
        deallocate(ww3_var)
!       deallocate(ww3_var2)
       endif              !dddddddddddd>
      end subroutine wave_allocate_ww3

!..................................................

      subroutine wave_reset
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='wave_reset'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

        call wave_allocate_s

      end subroutine wave_reset

!..................................................

      subroutine wave_moveforward
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='wave_moveforward'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do j=0,jmax+1
      do i=0,imax+1
        t_wave_t(i,j,0)=     t_wave_t(i,j,2)
         j_wave_t(i,j,0)=    j_wave_t(i,j,2)
         hs_wave_t(i,j,0)=  hs_wave_t(i,j,2)
        hsw_wave_t(i,j,0)= hsw_wave_t(i,j,2)
        foc_wave_t(i,j,0)= foc_wave_t(i,j,2)
         k_wave_t(i,j,0) =   k_wave_t(i,j,2)
         kx_wave_t(i,j,0)=  kx_wave_t(i,j,2)
         ky_wave_t(i,j,0)=  ky_wave_t(i,j,2)
       twox_wave_t(i,j,0)=twox_wave_t(i,j,2)
       twoy_wave_t(i,j,0)=twoy_wave_t(i,j,2)
       tawx_wave_t(i,j,0)=tawx_wave_t(i,j,2)
       tawy_wave_t(i,j,0)=tawy_wave_t(i,j,2)
      enddo
      enddo

      end subroutine wave_moveforward

!..................................................

      subroutine wave_linear_in_time_02
      use module_principal
      implicit none
      real*4 c1_,c2_
#ifdef synopsis
       subroutinetitle='wave_linear_in_time_02'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! entre 2 echeances: la valeur instantannée est donnée par l'interpolation
! lineaire dans le temps des 2 echeances qui encadrent le temps present
      c1_=1.-rap_wave ; c2_=rap_wave

      do j=0,jmax+1
      do i=0,imax+1
          t_wave_t(i,j,1)=c1_*   t_wave_t(i,j,0)+c2_*   t_wave_t(i,j,2)
         hs_wave_t(i,j,1)=c1_*  hs_wave_t(i,j,0)+c2_*  hs_wave_t(i,j,2)
        hsw_wave_t(i,j,1)=c1_* hsw_wave_t(i,j,0)+c2_* hsw_wave_t(i,j,2)
        foc_wave_t(i,j,1)=c1_* foc_wave_t(i,j,0)+c2_* foc_wave_t(i,j,2)
         kx_wave_t(i,j,1)=c1_*  kx_wave_t(i,j,0)+c2_*  kx_wave_t(i,j,2)
         ky_wave_t(i,j,1)=c1_*  ky_wave_t(i,j,0)+c2_*  ky_wave_t(i,j,2)
          k_wave_t(i,j,1)=c1_*   k_wave_t(i,j,0)+c2_*   k_wave_t(i,j,2)
       twox_wave_t(i,j,1)=c1_*twox_wave_t(i,j,0)+c2_*twox_wave_t(i,j,2)
       twoy_wave_t(i,j,1)=c1_*twoy_wave_t(i,j,0)+c2_*twoy_wave_t(i,j,2)
       tawx_wave_t(i,j,1)=c1_*tawx_wave_t(i,j,0)+c2_*tawx_wave_t(i,j,2)
       tawy_wave_t(i,j,1)=c1_*tawy_wave_t(i,j,0)+c2_*tawy_wave_t(i,j,2)
          j_wave_t(i,j,1)=c1_*   j_wave_t(i,j,0)+c2_*   j_wave_t(i,j,2)
      enddo
      enddo

      do j=0,jmax+1
      do i=0,imax+1
        dir_wave_t(i,j,1)=atan2(ky_wave_t(i,j,1),kx_wave_t(i,j,1))
      enddo
      enddo

      end subroutine wave_linear_in_time_02

!..................................................

      subroutine wave_linear_in_time_01
      use module_principal
      implicit none
      double precision time_
#ifdef synopsis
       subroutinetitle='wave_linear_in_time_01'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! interpolation temporelle evitant le stokage de 2 echeances:      !26-11-10
! soit v(2) le courant de l'echeance 2 (echeance en avant).
! soit k2 la valeur de l'iteration correspondant à la date de l'echeance 2.
! soit v(kount) le courant interpolé temporellement pour l'iteration en cours.
! entre 2 echeances le courant varie lineairement et par consequent on a:

! ( v(2)-v(kount) )/( k2 - kount ) = ( v(2) - v(kount-1) )/( k2 - (kount-1) )

! on en deduit v(kount) en fonction de v(2) et de la valeur precedente v(kount-1):

! v(kount)=rap*v(2)+(1.-rap)*v(kount-1) avec rap=1/(k2-k+1)

! k2 est donné par (int((real(kount)-wavedt(2))/wavedt(1) )+1.)*wavedt(1)+wavedt(2)

!     k2=(int((real(kount)-wavedt(2))/wavedt(1))+1.)*wavedt(1)+wavedt(2)
!     rap=1./max(k2-kount+1.d0,un)                                      !28-11-10
      time_=(int( (elapsedtime_now-wavedt(2))              &    !17-04-11
                    /wavedt(1))+1.)*wavedt(1)+wavedt(2)
      rap=(elapsedtime_now-elapsedtime_bef)                   &    !17-04-11
             /max(time_-elapsedtime_bef,un)

      x2=rap ; x1=1.-rap
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
       velstokes_u(i,j,k,1)=(x1*velstokes_u(i,j,k,1) &
                            +x2*velstokes_u(i,j,k,2) &
                                      /dz_u(i,j,k,1) & ! /dz seulement si velstokes(2) homogene à transport
                                )*wetmask_u(i,j)

       velstokes_v(i,j,k,1)=(x1*velstokes_v(i,j,k,1) &
                            +x2*velstokes_v(i,j,k,2) &
                                      /dz_v(i,j,k,1) & ! /dz seulement si velstokes(2) homogene à transport
                                )*wetmask_v(i,j)
      enddo
      enddo
      enddo

      end subroutine wave_linear_in_time_01

!..................................................

      subroutine wave_ubw_fw
      use module_principal ; use module_sedw
      implicit none
      real lw0_,abwoverz0_,epswoversigma_,abw_
#ifdef synopsis
       subroutinetitle='wave_ubw_fw'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Vitesse orbitale près du fond.
! Loin du fond uorb(z) = sigma * Hs * cosh(k (z+h)) /((sqrt(8)*sinh(kD))

      do j=1,jmax
      do i=1,imax

        lw0_=grav*t_wave_t(i,j,1)**2/(2.*pi)

        if(h_w(i,j)/lw0_>2.) then !------------>

          ubw(i,j)=0.
          abw_=small1

        else                      !------------>

         ubw(i,j)=pi/t_wave_t(i,j,1)*hs_wave_t(i,j,1)/sqrt(2.)      &
                  /sinh(min(max(k_wave_t(i,j,1)*hz_w(i,j,2),1.0),50.d0))    ! on recalcule la vitesse orbitale
        ! On seuille k_wave_t * hz_w à 1 mètre

         abw_=small1+ubw(i,j)*t_wave_t(i,j,1)/(2.*pi)     ! et la demi excursion près du fond

        endif                        !------------>

!      if(ubw(i,j)< 0)stop 'ubw< 0'
!      if(abw_<=0)stop 'abw<=0'
!      if(t_wave_t(i,j,1)<=0)stop 't_wave_t(i,j,1)<=0'
!      if(z0_w(i,j)<=0.)stop 'z0_w(i,j)<=0'

       !abwoverz0_=abw_/z0_w(i,j)
       abwoverz0_=abw_/z0b
!       abwoverz0_=abw_/0.001

! Alternative 1
       if(abwoverz0_>200.and.abwoverz0_<11000) then  !.............>

              fw(i,j)=1.39*( abwoverz0_ )**(-0.52)

          else if (abwoverz0_<=200.) then

              fw(i,j)=18.*( abwoverz0_ )**(-1)

          else if (abwoverz0_>=11000) then

              fw(i,j)=0.112*( abwoverz0_ )**(-0.25)

          endif                                     !.............>

!! Alternative 2  pour info : à coder correctement
!        ks=z0b*30.
!        if(abw_/ks<=1.57)then
!         fw(i,j)=0.3
!        else
!         fw(i,j)=0.00251*exp(5.21*(abw_/ks)**(-0.19))
!        endif
       
!       if(i+par%timax(1)==141.and.j+par%tjmax(1)==140) then
!        write(69,*)'hz_w(i,j,2)=',hz_w(i,j,2)
!        write(69,*)'k_wave_t(i,j,1)=',k_wave_t(i,j,1)
!        write(69,*)'abwoverz0_=',abwoverz0_
!        write(69,*)'hs_wave_t(i,j,1)=',hs_wave_t(i,j,1)
!        write(69,*)'abw_=',abw_
!        write(69,*)'z0_w(i,j)=',z0_w(i,j)
!      !  write(69,*)'------------------'
!      endif
!       if(i+par%timax(1)==141.and.j+par%tjmax(1)==141) then
!        write(70,*)'hz_w(i,j,2)=',hz_w(i,j,2)
!        write(70,*)'k_wave_t(i,j,1)=',k_wave_t(i,j,1)
!        write(70,*)'abwoverz0_=',abwoverz0_
!        write(70,*)'hs_wave_t(i,j,1)=',hs_wave_t(i,j,1)
!        write(70,*)'abw_=',abw_
!        write(70,*)'z0_w(i,j)=',z0_w(i,j)
!      !  write(70,*)'------------------'
!      endif

      enddo
      enddo

      end subroutine wave_ubw_fw

!..................................................

      subroutine wave_ubw_fws
      use module_principal ; use module_sedw
      implicit none
      real lw0_,abwoverz0_,epswoversigma_,abw_
#ifdef synopsis
       subroutinetitle='wave_ubw_fws'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! as wave_ubw_fw but with z0_sf replacing z0_w

      do j=1,jmax
      do i=1,imax

        lw0_=grav*t_wave_t(i,j,1)**2/(2.*pi)

        if(h_w(i,j)/lw0_>2.) then !------------>

          ubw(i,j)=0.
          abw_=small1

        else                      !------------>

         ubw(i,j)=pi/t_wave_t(i,j,1)*hs_wave_t(i,j,1)/sqrt(2.)      &
                  /sinh(min(max(k_wave_t(i,j,1)*hz_w(i,j,2),1.0),50.d0))    ! on recalcule la vitesse orbitale

         abw_=small1+ubw(i,j)*t_wave_t(i,j,1)/(2.*pi)     ! et la demi excursion près du fond

        endif                        !------------>

!      if(ubw(i,j)< 0)stop 'ubw< 0'
!      if(abw_<=0)stop 'abw<=0'
!      if(t_wave_t(i,j,1)<=0)stop 't_wave_t(i,j,1)<=0'

!       abwoverz0_=abw_/z0_sf(i,j)
       abwoverz0_=abw_/z0b

       if(abwoverz0_>200.and.abwoverz0_<11000) then  !.............>

              fw(i,j)=1.39*( abwoverz0_ )**(-0.52)

          else if (abwoverz0_<=200.) then

              fw(i,j)=18.*( abwoverz0_ )**(-1)

          else if (abwoverz0_>=11000) then

              fw(i,j)=0.112*( abwoverz0_ )**(-0.25)

          endif                                     !.............>

      enddo
      enddo

      end subroutine wave_ubw_fws

!..................................................

      subroutine wave_bottom_momentum_production
      use module_principal ; use module_sedw
      implicit none
      real lw0_,abwoverz0_,epswoversigma_
#ifdef synopsis
       subroutinetitle='wave_bottom_momentum_production'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Compute ubw and fw:
      call wave_ubw_fw

! Compute the wave dissipation by bottom friction (uchiyama et al. 2010 - ums10)
      const1=rho/(2*sqrt(pi))/(2*pi)
      do j=1,jmax
      do i=1,imax

!      epswoversigma_=rho*fw(i,j)*(ubw(i,j))**3/(2*sqrt(pi)) & !epsilonw=rho* fw*(|uorb|)^3/(2*sqrt(pi)) eq (36) ums10
!                    /(2*pi/t_wave_t(i,j,1))                   ! / sigma
       epswoversigma_=const1*t_wave_t(i,j,1)*fw(i,j)*ubw(i,j)**3 ! equivalent au calcul ci dessus

       xy_t(i,j,1)=kx_wave_t(i,j,1)*epswoversigma_      !(epsilonw*kx/sigma) eq (57) ums10
       xy_t(i,j,2)=ky_wave_t(i,j,1)*epswoversigma_      !(epsilonw*ky/sigma)
       
!       if(i+par%timax(1)==141.and.j+par%tjmax(1)==140) then
!        write(69,*)'xy_t(i,j,1)=',xy_t(i,j,1)
!        write(69,*)'kx_wave_t(i,j,1)=',kx_wave_t(i,j,1)
!        write(69,*)'ky_wave_t(i,j,1)=',ky_wave_t(i,j,1)
!        write(69,*)'epswoversigma_=',epswoversigma_
!        write(69,*)'t_wave_t(i,j,1)=',t_wave_t(i,j,1)
!        write(69,*)'fw(i,j)=',fw(i,j)
!        write(69,*)'ubw(i,j)=',ubw(i,j)
!        write(69,*)'------------------'
!      endif
!       if(i+par%timax(1)==141.and.j+par%tjmax(1)==141) then
!        write(70,*)'xy_t(i,j,1)=',xy_t(i,j,1)
!        write(70,*)'kx_wave_t(i,j,1)=',kx_wave_t(i,j,1)
!        write(70,*)'ky_wave_t(i,j,1)=',ky_wave_t(i,j,1)
!        write(70,*)'epswoversigma_=',epswoversigma_
!        write(70,*)'t_wave_t(i,j,1)=',t_wave_t(i,j,1)
!        write(70,*)'fw(i,j)=',fw(i,j)
!        write(70,*)'ubw(i,j)=',ubw(i,j)
!        write(70,*)'------------------'
!      endif

      enddo
      enddo

! Interpolate on u nodes:
      do j=2,jmax-1
      do i=2,imax
       wstresb_u(i,j)=0.5*(xy_t(i,j,1)+xy_t(i-1,j,1))
      enddo
      enddo

! Interpolate on v nodes:
      do j=2,jmax
      do i=2,imax-1
       wstresb_v(i,j)=0.5*(xy_t(i,j,2)+xy_t(i,j-1,2))
      enddo
      enddo

      end subroutine wave_bottom_momentum_production

!..................................................

      subroutine wave_surface_momentum_production
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='wave_surface_momentum_production'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(iairsea==0) then !m°v°m> !23-10-20
! Si on n'appelle pas module_airseaflux remise A zero de wstress_u et wstress_v !23-10-20
       wstress_u(:,:,1)=0.  ; wstress_v(:,:,1)=0.      
      endif               !m°v°m> !23-10-20

      do j=2,jmax-1
      do i=2,imax
! totalstress=twind+tao-taw (u comp.):
       wstress_u(i,j,1)=                                               &
       wstress_u(i,j,1)+0.5*(twox_wave_t(i  ,j,1)-tawx_wave_t(i  ,j,1) &
                            +twox_wave_t(i-1,j,1)-tawx_wave_t(i-1,j,1))
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
       wstress_v(i,j,1)=                                               &
       wstress_v(i,j,1)+0.5*(twoy_wave_t(i,j  ,1)-tawy_wave_t(i,j  ,1) &
                            +twoy_wave_t(i,j-1,1)-tawy_wave_t(i,j-1,1))
      enddo
      enddo

! Boundary condition on wstress_u & wstress_v in order to compute wstress_w
      call obc_surfstress                      !25-09-10

! rebuild the total stress magnitude for the tKe equation !21-05-10
      do j=1,jmax
      do i=1,imax
       wstress_w(i,j)=sqrt( (0.5*( wstress_u(i  ,j  ,1)          &
                                  +wstress_u(i+1,j  ,1)))**2     &
                           +(0.5*( wstress_v(i  ,j  ,1)          &
                                  +wstress_v(i  ,j+1,1)))**2 )
      enddo
      enddo

      end subroutine wave_surface_momentum_production

!..................................................

      subroutine wave_sj_sshear
      use module_principal
      implicit none
      double precision ushear_,vshear_,sshearsurf_
#ifdef synopsis
       subroutinetitle='wave_sj_sshear'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! note: ipt=1 jpt=1 imt=0 jmt=0
       do j=1,jmax
       do i=1,imax
!! Les lignes ci-dessous sont commentees car entrainant trop d'instabilites.         !29-01-13
!       ushear_=0.25*( (vel_u(i+ipt,j,kmax,1)-vel_u(i+ipt,j,kmax-1,1)) &
!                  /(depth_u(i+ipt,j,kmax)-depth_u(i+ipt,j,kmax-1))   &
!                  +  (vel_u(i-imt,j,kmax,1)-vel_u(i-imt,j,kmax-1,1)) &
!                  /(depth_u(i-imt,j,kmax)-depth_u(i-imt,j,kmax-1)))  &
!              + 0.25*( (vel_u(i+ipt,j,kmax,0)-vel_u(i+ipt,j,kmax-1,0)) &
!                  /(depth_u(i+ipt,j,kmax)-depth_u(i+ipt,j,kmax-1))   &
!                  +  (vel_u(i-imt,j,kmax,0)-vel_u(i-imt,j,kmax-1,0)) &
!                  /(depth_u(i-imt,j,kmax)-depth_u(i-imt,j,kmax-1)))


!       vshear_=0.25*( (vel_v(i,j+jpt,kmax,1)-vel_v(i,j+jpt,kmax-1,1)) &
!                  /(depth_v(i,j+jpt,kmax)-depth_v(i,j+jpt,kmax-1))   &
!                  +  (vel_v(i,j-jmt,kmax,1)-vel_v(i,j-jmt,kmax-1,1)) &
!                  /(depth_v(i,j-jmt,kmax)-depth_v(i,j-jmt,kmax-1)))  &
!             +0.25*( (vel_v(i,j+jpt,kmax,0)-vel_v(i,j+jpt,kmax-1,0)) &
!                  /(depth_v(i,j+jpt,kmax)-depth_v(i,j+jpt,kmax-1))   &
!                  +  (vel_v(i,j-jmt,kmax,0)-vel_v(i,j-jmt,kmax-1,0)) &
!                  /(depth_v(i,j-jmt,kmax)-depth_v(i,j-jmt,kmax-1)))
!       sshearsurf_=                                                   & ! From Ardhuin et al 2008 Eq40 surface terms
!       -(0.0625*(hs_wave_t(i,j,1)**2))*(                              & ! -E=-(hs**2)/16*(
!                2.*pi/t_wave_t(i,j,1)                                 & ! sigma
!                     /k_wave_t(i,j,1)                                 & !      /k
!                                     *( kx_wave_t(i,j,1)*ushear_      & !        *(kx.du/dz
!                                       +ky_wave_t(i,j,1)*vshear_ )    & !                  +ky.dv/dz)
!                                   *tanh(k_wave_t(i,j,1)*hz_w(i,j,1)) & !                            *tanh(kD)
!                +0.5*(ushear_**2+vshear_**2)     )                      ! +m/2*(du/dz**2+du/dz**2)
! Limitations pour diminuer les instabilites
!       if (sshearsurf_.gt.0.) sshearsurf_=min(0.5,sshearsurf_)
!       if (sshearsurf_.lt.0.) sshearsurf_=max(-0.5,sshearsurf_)


!        sshstokes_w(i,j)=-(j_wave_t(i,j,1)+sshearsurf_)/grav
        sshstokes_w(i,j)=-(j_wave_t(i,j,1))/grav                   !29-01-13
       enddo
       enddo

      end subroutine wave_sj_sshear

!..................................................

      subroutine wave_allocate_s
!     use module_principal ; use module_sedw
      use module_principal
      implicit none
      integer lb_
#ifdef synopsis
       subroutinetitle='wave_allocate_s'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(initial_main_status==1)then !::::::::::>
       write(6,*)'If Restart files have been loaded it'
       write(6,*)'is too late for dynamic allocation in routine'
       write(6,*)'wave_allocate_s'
       stop 'stop in wave_allocate_s'
      endif                          !::::::::::>

      if(txt_wavemodeltype/='ww3_native_grid') then !25-07-22
       lb_=0
       lb_=1
!      stop 'Verifier que lb_ ne devrait pas etre 1'
      else
       lb_=1
       stop 'Verifier que lb_ ne devrait pas etre 0'
      endif

#ifdef key_oasis_symp_ww3
!Case symphonie-ww3-oasis :
        if(allocated(t_wave_t))deallocate(   t_wave_t)
        if(allocated(hs_wave_t))deallocate(  hs_wave_t)
        if(allocated(j_wave_t))deallocate(   j_wave_t)
        if(allocated(hsw_wave_t))deallocate( hsw_wave_t)
        if(allocated(foc_wave_t))deallocate( foc_wave_t)
        if(allocated(k_wave_t))deallocate(   k_wave_t)
        if(allocated(kx_wave_t))deallocate(  kx_wave_t)
        if(allocated(ky_wave_t))deallocate(  ky_wave_t)
        if(allocated(twox_wave_t))deallocate(twox_wave_t)
        if(allocated(twoy_wave_t))deallocate(twoy_wave_t)
        if(allocated(tawx_wave_t))deallocate(tawx_wave_t)
        if(allocated(tawy_wave_t))deallocate(tawy_wave_t)
        if(allocated(dir_wave_t))deallocate( dir_wave_t)
        if(allocated(uss_wave_t))deallocate( uss_wave_t)
        if(allocated(vss_wave_t))deallocate( vss_wave_t)
        if(allocated(mask_wave_t))deallocate( mask_wave_t)
        if(allocated(ubw))deallocate(ubw)
        if(allocated(fw))deallocate( fw)
        if(allocated(wstresb_u))deallocate(wstresb_u)
        if(allocated(wstresb_v))deallocate(wstresb_v)
!......
        if(.not.allocated(t_wave_t))allocate(   t_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; t_wave_t=20.0
        if(.not.allocated(hs_wave_t))allocate(  hs_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; hs_wave_t=0
        if(.not.allocated(j_wave_t))allocate(   j_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; j_wave_t=0
        if(.not.allocated(hsw_wave_t))allocate( hsw_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; hsw_wave_t=z0s/0.64 !17-03-29
        if(.not.allocated(foc_wave_t))allocate( foc_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; foc_wave_t=0
        if(.not.allocated(k_wave_t))allocate(   k_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; k_wave_t=0
        if(.not.allocated(kx_wave_t))allocate(  kx_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; kx_wave_t=0
        if(.not.allocated(ky_wave_t))allocate(  ky_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; ky_wave_t=0
        if(.not.allocated(twox_wave_t))allocate(twox_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; twox_wave_t=0
        if(.not.allocated(twoy_wave_t))allocate(twoy_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; twoy_wave_t=0
        if(.not.allocated(tawx_wave_t))allocate(tawx_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; tawx_wave_t=0
        if(.not.allocated(tawy_wave_t))allocate(tawy_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; tawy_wave_t=0
        if(.not.allocated(dir_wave_t))allocate( dir_wave_t(0:imax+1,0:jmax+1,1 )) ; dir_wave_t=0
        if(.not.allocated(uss_wave_t))allocate( uss_wave_t(0:imax+1,0:jmax+1,2  )) ; uss_wave_t=0
        if(.not.allocated(vss_wave_t))allocate( vss_wave_t(0:imax+1,0:jmax+1,2 )) ; vss_wave_t=0
        if(.not.allocated(mask_wave_t))allocate( mask_wave_t(0:imax+1,0:jmax+1)) ; mask_wave_t=0
        if(.not.allocated(ubw))allocate(ubw(imax,jmax)) ; ubw=0
        if(.not.allocated(fw))allocate( fw(imax,jmax)) ; fw=0
        if(.not.allocated(wstresb_u))allocate(wstresb_u(0:imax+1,0:jmax+1)) ; wstresb_u=0
        if(.not.allocated(wstresb_v))allocate(wstresb_v(0:imax+1,0:jmax+1)) ; wstresb_v=0
!......
#else
!Case symphonie-ww3 by files
        allocate(   t_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; t_wave_t=1.  !13-09-16
        allocate(  hs_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; hs_wave_t=0
        allocate(   j_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; j_wave_t=0
        allocate( hsw_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; hsw_wave_t=z0s/0.64 !11-09-16
        allocate( foc_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; foc_wave_t=0
        allocate(   k_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; k_wave_t=0
        allocate(  kx_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; kx_wave_t=0
        allocate(  ky_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; ky_wave_t=0
        allocate(twox_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; twox_wave_t=0
        allocate(twoy_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; twoy_wave_t=0
        allocate(tawx_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; tawx_wave_t=0
        allocate(tawy_wave_t(0:imax+1,0:jmax+1,lb_:2)) ; tawy_wave_t=0

        allocate( dir_wave_t(0:imax+1,0:jmax+1,1    )) ; dir_wave_t=0 
        allocate( uss_wave_t(0:imax+1,0:jmax+1,2:2  )) ; uss_wave_t=0 !12-09-18
        allocate( vss_wave_t(0:imax+1,0:jmax+1,2:2  )) ; vss_wave_t=0 !12-09-18
        allocate( mask_wave_t(0:imax+1,0:jmax+1)) ; mask_wave_t=0

        allocate(ubw(imax,jmax)) ; ubw=0
        allocate( fw(imax,jmax)) ; fw=0

        allocate(wstresb_u(0:imax+1,0:jmax+1)) ; wstresb_u=0
        allocate(wstresb_v(0:imax+1,0:jmax+1)) ; wstresb_v=0
#endif
      end subroutine wave_allocate_s

!.....................................................................

      subroutine wave_read_lonlat(ncid_)
      implicit none
      integer ncid_
#ifdef synopsis
       subroutinetitle='wave_read_lonlat'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(ww3_lon))allocate(ww3_lon  (ww3_imax))
      if(.not.allocated(ww3_lat)) then !---->
      if(ww3_type_grid==type_structured)  allocate(ww3_lat  (ww3_jmax))
      if(ww3_type_grid==type_unstructured)allocate(ww3_lat  (ww3_imax))
      endif                            !---->
      if(.not.allocated(ww3cg_lon))allocate(ww3cg_lon(ww3_imax,ww3_jmax))
      if(.not.allocated(ww3cg_lat))allocate(ww3cg_lat(ww3_imax,ww3_jmax))

      varstart(1)=ww3zoom_istr ; varcount(1)=ww3_imax
      varstart(2)=ww3zoom_jstr ; varcount(2)=ww3_jmax

! Longitude:
                   status=nf_inq_varid(ncid_,'nav_lon',var_id)
      if(status/=0)status=nf_inq_varid(ncid_,'longitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid_,'lon',var_id)
      if(status/=0)stop 'erreur var_id longitude ww3cg'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'erreur nf_inq_var wave_grid_extraction_ww3cg'

      ksecu=0

      if(ww3_type_grid==type_structured) then !ssssssssssssssssssssssssssssssssss>

      if(var_type==nf_real.and.var_dims==1) then   !---->
       ksecu=1
       status=nf_get_vara_real(ncid_,var_id,varstart(1:1),varcount(1:1),ww3_lon(1:ww3_imax))
       if(status/=0)stop 'Err4739 nf_get_vara_real ww3_lon'
!      if(txt_wavemodeltype=='ww3_on_sgrid') then !-debug->
!       do i=1,ww3_imax 
!        write(10+par%rank,*)i+ww3zoom_istr-1,i,ww3_lon(i)
!       enddo
!      endif                                      !-debug->
      endif                                        !---->
      if(var_type==nf_double.and.var_dims==2) then !......>
       ksecu=1
       status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),ww3cg_lon)
       if(status/=0)stop 'erreur lecture ww3cg_lon'
      endif                                        !......>
      if(var_type==nf_real.and.var_dims==2) then   !......>
       ksecu=1
       allocate(ww3_2dr4_lon(ww3_imax,ww3_jmax))
       status=nf_get_vara_real(ncid_,var_id,varstart(1:2),varcount(1:2),ww3_2dr4_lon)
       if(status/=0)call wave_error('nf_get_vara_real','longitude',4319)

       if(txt_wavemodeltype=='ww3_on_sgrid') then !jesuischarlie>
! ici on verifie que les longitude de la grille ww3 sont bien les memes que celles de S:
       i0=0 ; if(iperiodicboundary)i0=2
       do j=0,jmax+1 ; do i=i0,imax+1-i0
       i1=i+1-i0 ; j1=j+1

       if(abs(ww3_2dr4_lon(i1,j1)-lon_t(i,j)*rad2deg)>0.001) then
        write(10+par%rank,*)'par%rank=',par%rank
        write(10+par%rank,*)'i0=',i0
        write(10+par%rank,*)'i,j',i,j
        write(10+par%rank,*)'i1,j1',i1,j1
        write(10+par%rank,*)'ww3_2dr4_lon(i1,j1),lon_t(i,j)*rad2deg' &
                            ,ww3_2dr4_lon(i1,j1),lon_t(i,j)*rad2deg
        stop 'ww3 grid does not match S "w "grid see error files fort.n'
       endif
       enddo ; enddo

       endif                                   !jesuischarlie>
      endif                                        !......>

      endif                                   !ssssssssssssssssssssssssssssssssss>

      if(ww3_type_grid==type_unstructured) then !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>
      if(var_type==nf_real.and.var_dims==1) then   !---->
       ksecu=1
       status=nf_get_vara_real(ncid_,var_id,varstart(1:1),varcount(1:1),ww3_lon(1:ww3_imax))
       if(status/=0)stop 'Err4777 nf_get_vara_real ww3_lon'
      endif                                        !---->
      endif                                     !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>

      if(ksecu==0)stop 'wave_grid_extraction_ww3cg cas non prevu'

! Latitude:
                   status=nf_inq_varid(ncid_,'nav_lat',var_id)
      if(status/=0)status=nf_inq_varid(ncid_,'latitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid_,'lat',var_id)
      if(status/=0)stop 'erreur var_id latitude ww3cg'
      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'erreur nf_inq_var wave_grid_extraction_ww3cg'

      if(ww3_type_grid==type_structured) then !ssssssssssssssssssssssssssssssssss>
      if(var_type==nf_real.and.var_dims==1) then   !---->
       status=nf_get_vara_real(ncid_,var_id,varstart(2:2),varcount(2:2),ww3_lat(1:ww3_jmax))
       if(status/=0)stop 'erreur nf_get_vara_real ww3_lat'
!      if(txt_wavemodeltype=='ww3_on_sgrid') then !-debug->
!       do j=1,ww3_jmax 
!        write(10+par%rank,*)j+ww3zoom_jstr-1,j,ww3_lat(j)
!       enddo
!      endif                                      !-debug->
       do j=1,ww3_jmax ! 14-05-13
       do i=1,ww3_imax
        ww3cg_lon(i,j)=ww3_lon(i)
        ww3cg_lat(i,j)=ww3_lat(j)
       enddo
       enddo
      endif                                        !---->
      if(var_type==nf_double.and.var_dims==2) then !......>
       ksecu=1
       status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),ww3cg_lat)
       if(status/=0)stop 'erreur lecture ww3cg_lon'
      endif                                        !......>
      if(var_type==nf_real.and.var_dims==2) then   !......>
       ksecu=1
       allocate(ww3_2dr4_lat(ww3_imax,ww3_jmax))
       status=nf_get_vara_real(ncid_,var_id,varstart(1:2),varcount(1:2),ww3_2dr4_lat)
       if(status/=0)call wave_error('nf_get_vara_real','latgitude',4319)

       if(txt_wavemodeltype=='ww3_on_sgrid') then !jesuischarlie>
! ici on verifie que les longitude de la grille ww3 sont bien les memes que celles de S:
       i0=0 ; if(iperiodicboundary)i0=2
       do j=0,jmax+1 ; do i=i0,imax+1-i0
       i1=i+1-i0 ; j1=j+1
       if(abs(ww3_2dr4_lat(i1,j1)-lat_t(i,j)*rad2deg)>0.001) then
        write(10+par%rank,*)'par%rank=',par%rank
        write(10+par%rank,*)'i0=',i0
        write(10+par%rank,*)'i,j',i,j
        write(10+par%rank,*)'i1,j1',i1,j1
        write(10+par%rank,*)'ww3_2dr4_lat(i1,j1),lat_t(i,j)*rad2deg' &
                            ,ww3_2dr4_lat(i1,j1),lat_t(i,j)*rad2deg
        stop 'ww3 grid does not match S "w "grid see error files fort.n'
       endif
       enddo ; enddo
       else                                    !jesuischarlie>
           do j=1,ww3_jmax ; do i=1,ww3_imax !24-09-16
            ww3cg_lon(i,j)=ww3_2dr4_lon(i,j)
            ww3cg_lat(i,j)=ww3_2dr4_lat(i,j)
           enddo           ; enddo
           deallocate(ww3_2dr4_lon) ; deallocate(ww3_2dr4_lat)
       endif                                   !jesuischarlie>
      endif                                        !......>
      endif                                   !ssssssssssssssssssssssssssssssssss>

      if(ww3_type_grid==type_unstructured) then !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>
       status=nf_get_vara_real(ncid_,var_id,varstart(1:1),varcount(1:1),ww3_lat(1:ww3_imax))
       if(status/=0)stop 'erreur nf_get_vara_real ww3_lat'
      endif                                     !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>

      end subroutine wave_read_lonlat

!..............................................................

      subroutine wave_sgrid_to_ww3grid
      implicit none
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)    &
                        ,dlon_di_,dlon_dj_,dlon_dm_
#ifdef synopsis
       subroutinetitle='wave_sgrid_to_ww3grid'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      allocate(ij2ww3_i   (0:imax+1,0:jmax+1))
      allocate(ij2ww3_j   (0:imax+1,0:jmax+1))
      allocate(ij2ww3_teta(0:imax+1,0:jmax+1))

      if(txt_wavemodeltype=='ww3_on_sgrid') then !>>>>>>>>>>>>>>>>>>

! Cas trivial ou les grilles ww3 et S sont strictement les memes:
! ATTENTION, comme nous commencons les boucles a zero mais que le
! tableau ww3_var commence sa numerotation a 1, on applique un shift de 1.
! D'autre part nous raisonnons sur la grille globale et donc en plus nous
! appliquons les shifts par%timax(1) et par%tjmax(1)
       do j=0,jmax+1 ; do i=0,imax+1
        ij2ww3_i(i,j)=i+1+par%timax(1)
        ij2ww3_j(i,j)=j+1+par%tjmax(1)
       enddo         ; enddo

      else                                    !>>>>>>>>>>>>>>>>>>
! Cas general ou les grille ww3 et S sont differentes:

      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

      x2=real(ww3_imax/2)
      x3=real(ww3_jmax/2)
      deci=x2
      decj=x3

! ETAPE 1: trouver les coordonnées dans la grille ORCA:

! First guess: centre du domaine:
      k10=0
 1456 continue

! Principe suppose une relation lineaire entre lat lon et indice de grille.
!      dlon/di*Di+dlon/dj*Dj=Dlon
!      dlat/di*Di+dlat/dj*Dj=Dlat
! On cherche Di et Dj correspondant à Dlon=lon(i,j)-lonmeteo(i0,j0)
!                                et à Dlat=lat(i,j)-latmeteo(i0,j0)

      i1=int(deci)
      j1=int(decj)

      do j2=0,1
      do i2=0,1
      i0=i1+i2
      j0=j1+j2

      dlon_di_=ww3cg_lon(i0+1,j0  )-ww3cg_lon(i0-1,j0  )
      dlon_dj_=ww3cg_lon(i0  ,j0+1)-ww3cg_lon(i0  ,j0-1)
      dlon_dm_=rad2deg*lon_t(i  ,j)-ww3cg_lon(i0  ,j0)

      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(dlon_dm_<-180.)dlon_dm_=dlon_dm_+360.
      if(dlon_dm_> 180.)dlon_dm_=dlon_dm_-360.

! Determinant principal:
      x1=( dlon_di_                                             &
          *(ww3cg_lat(i0  ,j0+1)-ww3cg_lat(i0  ,j0-1))          &
          -(ww3cg_lat(i0+1,j0  )-ww3cg_lat(i0-1,j0  ))          &
          *dlon_dj_ )*0.25

      deci_(i2,j2)=min(max(                                  &
       i0+( dlon_dm_                                            &
          *(ww3cg_lat(i0  ,j0+1)-ww3cg_lat(i0,j0-1))            &
          -(rad2deg*lat_t(i,j)  -ww3cg_lat(i0,j0))              &
          *dlon_dj_)/x1*0.5    &
                   ,2.0001d0),ww3_imax-1.0001d0)

      decj_(i2,j2)=min(max(                                  &
       j0+( dlon_di_                                            &
          *(rad2deg*lat_t(i,j)  -ww3cg_lat(i0  ,j0))            &
          -(ww3cg_lat(i0+1,j0  )-ww3cg_lat(i0-1,j0))            &
          *dlon_dm_   )/x1*0.5    &
                   ,2.0001d0),ww3_jmax-1.0001d0)

!      if(par%rank==0)write(6,*)deci_(i2,j2),decj_(i2,j2)
      enddo
      enddo

      rapi=deci-i1
      rapj=decj-j1

      deci=(1.-rapi)*(1.-rapj)*deci_(0,0)   &
          +(1.-rapi)*    rapj *deci_(0,1)   &
          +    rapi *    rapj *deci_(1,1)   &
          +    rapi *(1.-rapj)*deci_(1,0)
      decj=(1.-rapi)*(1.-rapj)*decj_(0,0)   &
          +(1.-rapi)*    rapj *decj_(0,1)   &
          +    rapi *    rapj *decj_(1,1)   &
          +    rapi *(1.-rapj)*decj_(1,0)

! Si le point visé est different du first guess refaire le calcul
! avec un first guess donné par le dernier point visé:
      if(sqrt( (deci-x2)**2+(decj-x3)**2 ).gt.0.001)then
       x2=deci
       x3=decj
       k10=k10+1
       if(k10>20)then !!!!!!>
         if(par%rank==0)write(6,*)'par%rank=',par%rank
         if(par%rank==0)write(6,*)'(i,j)=   ',i,j
         if(par%rank==0)write(6,*)'deci decj',deci,decj
         stop 'hr_to_lr ne converge pas dans le forfait'
       endif          !!!!!!>
       goto 1456
      endif

! ETAPE 2: CALCULER L'ANGLE D'ORIENTATION LOCALE DE LA GRILLE ORCA

      i1=int(deci) ; j1=int(decj)
      rapi=deci-i1 ; rapj=decj-j1
      do j2=0,1
      do i2=0,1
       i0=i1+i2
       j0=j1+j2

       dlon_di_=ww3cg_lon(i0+1,j0  )-ww3cg_lon(i0-1,j0  )
       if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
       if(dlon_di_> 180.)dlon_di_=dlon_di_-360.

       dy_(i2,j2)= ww3cg_lat(i0+1,j0)-ww3cg_lat(i0-1,j0)
       dx_(i2,j2)=dlon_di_           &
                            *cos(ww3cg_lat(i0,j0)*deg2rad)
      enddo
      enddo
      x1=(1.-rapi)*(1.-rapj)*dx_(0,0)   &
        +(1.-rapi)*    rapj *dx_(0,1)   &
        +    rapi *    rapj *dx_(1,1)   &
        +    rapi *(1.-rapj)*dx_(1,0)
      y1=(1.-rapi)*(1.-rapj)*dy_(0,0)   &
        +(1.-rapi)*    rapj *dy_(0,1)   &
        +    rapi *    rapj *dy_(1,1)   &
        +    rapi *(1.-rapj)*dy_(1,0)

      ij2ww3_i(i,j)=deci
      ij2ww3_j(i,j)=decj
      ij2ww3_teta(i,j)=atan2(y1,x1)

      enddo
      enddo

      endif                                    !>>>>>>>>>>>>>>>>>>

      end subroutine wave_sgrid_to_ww3grid

!.........................................................

      subroutine wave_us3d_from_uss
      implicit none
      integer loop_ 
#ifdef synopsis
       subroutinetitle='wave_us3d_from_uss'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!      do loop_=0,nbdom-1

!      if(par%rank==loop_) then !>>>>>>
!      write(6,*)'loop_=',loop_

      const1=0.5*grav/pi
      const2=2.*pi
      do j=0,jmax+1
      do i=0,imax+1

! Cas "petites profondeurs":
       if((h_w(i,j)+depth_w(i,j,kmax+1))*k_wave_t(i,j,2)<6.)  then !hhhhhhhhhhhhhhhhh>

       x4=cosh(2.*k_wave_t(i,j,2)*(h_w(i,j)+depth_w(i,j,kmax+1))) ! x4=cosh(2kn(h+ssh))
       do k=kmin_w(i,j),kmax

! méthode 1: valeur ponctuelle du profil:
!       x3=cosh(2.*k_wave_t(i,j,2)*(depth_t(i,j,k)+h_w(i,j)))/x4
! méthode 2: profil moyénné sur dz
        x3=0.5/k_wave_t(i,j,2)/x4/(depth_w(i,j,k+1)-depth_w(i,j,k))           &
           *( sinh(2.*k_wave_t(i,j,2)*(depth_w(i,j,k+1)+h_w(i,j)))            &
             -sinh(2.*k_wave_t(i,j,2)*(depth_w(i,j,k  )+h_w(i,j))))
!       x3=1. ! cas academique profil constant
        anyv3d(i,j,k,1)=uss_wave_t(i,j,2)*x3
        anyv3d(i,j,k,2)=vss_wave_t(i,j,2)*x3

       enddo ! fin de boucle sur k

       else                                          !hhhhhhhhhhhhhhhhh>

! Cas "grandes profondeurs":
       do k=kmin_w(i,j),kmax

! méthode 1: valeur ponctuelle du profil exponentiel
!       x3=exp(2.*k_wave_t(i,j,2)*(depth_t(i,j,k)-ssh_int_w(i,j,1)))
! méthode 2: profil exponentiel moyenné sur dz:
       !if(par%rank==3) then  
       !   write(6,*)'k_wave_t(i,j,2)=',k_wave_t(i,j,2)
       !   write(6,*)'depth_w(i,j,k+1) en',i,j,k+1,depth_w(i,j,k+1)
       !   write(6,*)'depth_w(i,j,k) en',i,j,k,depth_w(i,j,k)
       !   write(6,*)'ssh_int_w(i,j,1) en',i,j,ssh_int_w(i,j,1)
       !endif

        x3=0.5/k_wave_t(i,j,2)/(depth_w(i,j,k+1)-depth_w(i,j,k))                    &
                    *( exp(2.*k_wave_t(i,j,2)*(depth_w(i,j,k+1)-ssh_int_w(i,j,1)))  &
                      -exp(2.*k_wave_t(i,j,2)*(depth_w(i,j,k  )-ssh_int_w(i,j,1))))
!       x3=1. ! cas academique profil constant
        anyv3d(i,j,k,1)=uss_wave_t(i,j,2)*x3
        anyv3d(i,j,k,2)=vss_wave_t(i,j,2)*x3

       enddo ! fin de boucle sur k
       endif                                         !hhhhhhhhhhhhhhhhh>

! Couches fusionnees:
       do k=1,kmin_w(i,j)-1 !17-11-19
        anyv3d(i,j,k,1)=anyv3d(i,j,kmin_w(i,j),1)
        anyv3d(i,j,k,2)=anyv3d(i,j,kmin_w(i,j),2)
       enddo

      enddo
      enddo

!     endif                    !>>>>>>
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      enddo ! boucle loop_

      do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1 !30-10-20

       velstokes_u(i,j,k,2)=                           &
            mask_u(i,j,k)*0.5*( anyv3d(i  ,j  ,k,1)    &
                               +anyv3d(i-1,j  ,k,1))

      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1 !30-10-20

       velstokes_v(i,j,k,2)=                           &
            mask_v(i,j,k)*0.5*( anyv3d(i  ,j  ,k,2)    &
                               +anyv3d(i  ,j-1,k,2))

      enddo ; enddo ; enddo

   
      end subroutine wave_us3d_from_uss

!.........................................................

      subroutine wave_kn_from_period
      implicit none
      integer iter_loop_
      double precision grav_over_2pi_,beta_,delta_c_,alpha_,c_
#ifdef synopsis
       subroutinetitle='wave_kn_from_period'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Compute the wave vector norm from the wave period


! Rustine cpl
!     do i=0,imax+1
!     do j=0,jmax+1
!     if (mask_t(i,j,kmax+1)==0) t_wave_t(i,j,1:2) = 20.0
!     enddo
!     enddo


#ifdef bidon
!**** RELATION DE DISPERSION CLASSIQUE - RESOLUTION ITERATIVE *****!

! S bathy or WW3 bathy ?
      if(flag_ww3_dpt==1) then !--ww3 native h available-->
       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,1)=2.*pi*max(dpt_wave_t(i,j),1.0)/t_wave_t(i,j,2)  ! beta ! on seuille à 1 mi pour les bancs découvrants
       enddo         ; enddo
      else                     !--ww3 native h not available-->
       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,1)=2.*pi*max(hz_w(i,j,1),1.0)/t_wave_t(i,j,2)      ! beta ! on seuille à 1m pour les bancs découvrants
       enddo         ; enddo
      endif                    !--ww3 native h not available-->


      grav_over_2pi_=grav/(2.*pi)
      do j=0,jmax+1
      do i=0,imax+1

! c_ is the phase speed
! On veut résoudre x=f(x)=alpha*tanh(beta/x) la vitesse de phase c_  etant l'inconnue x
! Le raisonnement est le suivant: f(x+dx)=x+dx
!                           soit: f(x)+dx f'(x)=x+dx
!                           soit: dx=(x-f(x))/(f'(x)-1)
! On estime une première valeur de x=alpha_ l'approximation pour c_ aux grandes profondeurs
!     beta_=2.*pi*max(hz_w(i,j,1),0.01d0)/t_wave_t(i,j,2)  ! beta
      beta_=xy_t(i,j,1)   ! beta
      alpha_=grav_over_2pi_*t_wave_t(i,j,2)         ! c=g*t/2pi first guess for large depth

      if(beta_/alpha_>200.) then !fgfgfgfgfgfgfg>      !07-05-11

! Cas c=first guess valide:
       c_=alpha_ ! phase speed for large depth

      else                       !fgfgfgfgfgfgfg>

! Cas calcul itératif:
       delta_c_=0.01                  ! dx, on estime une première valeur de dx
       c_=alpha_+delta_c_             ! alpha_+dx
         do iter_loop_=1,200

          delta_c_=(c_-alpha_*tanh(beta_/c_))                & !  (x-f(x))
                  /(-alpha_*beta_/(c_*cosh(beta_/c_))**2-1.)   !  /(f'(x)-1)
          c_=max(c_+delta_c_,1.d-4)
!         if(i==189.and.j==413)write(6,*)'c=',c_,alpha_,sqrt(grav*hz_w(i,j,1))
          if(abs(delta_c_)<0.01)goto 25
         enddo
         write(6,*)c_,delta_c_
         stop 'waveforcing: le calcul de kn ne converge pas'
   25  continue

      endif                      !fgfgfgfgfgfgfg>

!       k_wave_t(i,j,0)= k_wave_t(i,j,2)
!      kx_wave_t(i,j,0)=kx_wave_t(i,j,2)
!      ky_wave_t(i,j,0)=ky_wave_t(i,j,2)

!       k_wave_t(i,j,2)=2.*pi/(c_*t_wave_t(i,j,2)) ! Wave vector norm
        k_wave_t(i,j,2)=min(2.*pi/(c_*t_wave_t(i,j,2)),600.) ! 600=2pi/L avec L=1cm

!       if(i==imax/2.and.j==jmax/2.and.mask_t(i,j,kmax)==1) &
!       write(10+par%rank,*)'ANCIEN',k_wave_t(i,j,2),2*pi/k_wave_t(i,j,2)
        

      enddo
      enddo
#endif

!**** RELATION DE DISPERSION GUO 2002 - RESOLUTION DIRECTE *****!
      do j=0,jmax+1 ; do i=0,imax+1
! La formule de Guo 2002:
!      x=H*freq2pi/sqrt(grav*H)
!      kvectorG=(1./H)*( (x**2)*(1-exp(-x**2.4901))**(-1./2.4901) )
       x0=max(h_w(i,j),1.)
       x1=x0*(2.*pi/max(t_wave_t(i,j,2),0.001))/sqrt(grav*x0)
       k_wave_t(i,j,2)=min((1./x0)*( (x1**2)*(1-exp(-x1**2.4901))**(-1./2.4901) ),600.) ! 600=2pi/L avec L=1cm
!       if(i==imax/2.and.j==jmax/2.and.mask_t(i,j,kmax)==1) &
!       write(10+par%rank,*)'GUO   ',k_wave_t(i,j,2),2*pi/k_wave_t(i,j,2)
      enddo         ; enddo

      if(iwve==1) then !-----Michaud2012----> !21-11-19
       do j=0,jmax+1 ; do i=0,imax+1
        kx_wave_t(i,j,2)=k_wave_t(i,j,2)*cos(dir_wave_t(i,j,1))
        ky_wave_t(i,j,2)=k_wave_t(i,j,2)*sin(dir_wave_t(i,j,1))
       enddo         ; enddo
      endif            !-----Michaud2012----> !21-11-19

! Commente car seul k_wave_t (et non pas kx_wawe_t et ky_wave_t) compte
! pour le calcul du profil vertical du courant de stokes
!     if(iwve==2) then !-----McWilliams1999----> !21-11-19
!      do j=0,jmax+1 ; do i=0,imax+1
!       kx_wave_t(i,j,2)=k_wave_t(i,j,2)*uss_wave_t(i,j,2) &
!                                      /(uss_wave_t(i,j,2)**2+vss_wave_t(i,j,2)**2)
!       ky_wave_t(i,j,2)=k_wave_t(i,j,2)*vss_wave_t(i,j,2) &
!                 /(uss_wave_t(i,j,2)**2+vss_wave_t(i,j,2)**2)
!      enddo         ; enddo
!     endif            !-----McWilliams1999----> !21-11-19

      end subroutine wave_kn_from_period

!.........................................................

      subroutine wave_j_from_kn_hs
      implicit none
#ifdef synopsis
       subroutinetitle='wave_j_from_kn_hs'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       do j=0,jmax+1
       do i=0,imax+1

        x1=sqrt(kx_wave_t(i,j,2)**2+ky_wave_t(i,j,2)**2+small1) ! k

        j_wave_t(i,j,2)=grav*x1                  &  ! -J/g = -kE / sinh( 2kD ) ! Bennis et al, 2011
        *(0.25*                                  &
!              hs_wave_t(i,j,2)                  &
           min(hs_wave_t(i,j,2),0.5*hz_w(i,j,1)) &  ! 06-03-15
                               )**2              &  ! E=(hs/4)**2
          /sinh( min(2.*x1*max(hz_w(i,j,1),1.0),500.d0) )    ! sinh(2*k*D)

!        hs_wave_t(i,j,1)=min( hs_wave_t(i,j,1),0.5*hz_w(i,j,2)) !18-02-15

       enddo
       enddo

      end subroutine wave_j_from_kn_hs

!.........................................................

      subroutine wave_linear_in_time_s26(txt_)
      implicit none
      double precision interpcoef_
      real*4 c1_,c2_
      character*8 txt_
#ifdef synopsis
       subroutinetitle='wave_linear_in_time_s26'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(txt_=='standard') then !**************> ! Iterative phase
       ! Wave parameters at intermediate time steps : linear in time interpolation:
       ! between time=now and time=wavefile_nextime
       if(wavefile_nextime-elapsedtime_bef>0.) then !.....>

        interpcoef_=(wavefile_nextime-elapsedtime_now)    &
                   /(wavefile_nextime-elapsedtime_bef)

       else                                         !.....>
        stop 'wavefile_nextime-elapsedtime_bef<=0'
       endif                                        !.....>

       c1_=min(max(interpcoef_,0.d00),1.d00) ; c2_=1.-c1_

      endif                     !**************> ! Iterative phase

      if(txt_=='initial1') then !-- initialisation step1 -->
       c2_=1. ; c1_=0.
      endif                     !-- initialisation step1 -->
      if(txt_=='initial2') then !-- initialisation step2 -->
       c1_=(wavefile_nextime-elapsedtime_now)    &
          /(wavefile_nextime-wavefile_prvtime)
       c2_=1.-c1_
      endif                     !-- initialisation step2 -->

      if(iwve==1) then !-----Michaud2012----> !21-11-19
       do j=0,jmax+1 ; do i=0,imax+1
           t_wave_t(i,j,1)=c1_*   t_wave_t(i,j,1)+c2_*   t_wave_t(i,j,2)
          hs_wave_t(i,j,1)=c1_*  hs_wave_t(i,j,1)+c2_*  hs_wave_t(i,j,2)
         hsw_wave_t(i,j,1)=c1_* hsw_wave_t(i,j,1)+c2_* hsw_wave_t(i,j,2)
         foc_wave_t(i,j,1)=c1_* foc_wave_t(i,j,1)+c2_* foc_wave_t(i,j,2)
          kx_wave_t(i,j,1)=c1_*  kx_wave_t(i,j,1)+c2_*  kx_wave_t(i,j,2)
          ky_wave_t(i,j,1)=c1_*  ky_wave_t(i,j,1)+c2_*  ky_wave_t(i,j,2)
           k_wave_t(i,j,1)=c1_*   k_wave_t(i,j,1)+c2_*   k_wave_t(i,j,2)
        twox_wave_t(i,j,1)=c1_*twox_wave_t(i,j,1)+c2_*twox_wave_t(i,j,2)
        twoy_wave_t(i,j,1)=c1_*twoy_wave_t(i,j,1)+c2_*twoy_wave_t(i,j,2)
        tawx_wave_t(i,j,1)=c1_*tawx_wave_t(i,j,1)+c2_*tawx_wave_t(i,j,2)
        tawy_wave_t(i,j,1)=c1_*tawy_wave_t(i,j,1)+c2_*tawy_wave_t(i,j,2)
           j_wave_t(i,j,1)=c1_*   j_wave_t(i,j,1)+c2_*   j_wave_t(i,j,2)
       enddo ; enddo

       do j=0,jmax+1 ; do i=0,imax+1
         dir_wave_t(i,j,1)=atan2(ky_wave_t(i,j,1),kx_wave_t(i,j,1))
       enddo ; enddo
! Commentaire sur l'aspect discutable de l'usage du masque wetmask
! et du seuillage de hs_wave et hsw par hz dans
! ces circonstances rendues particulières par l'algo d'interpolation
! temporelle un peu special. En effet, on interpole entre "1" et "2" ce
! qui suppose qqpart que "1" ne doit pas être modifie par autre chose
! que l'interpolation. Or la on modifie par le wetmask et par le
! seuillage. Cela signifie qu'un decouvrement "laisse une memoire" dans le processus
! d'interpolation ce qui n'est pas terrible. A revoir plus tard par consequent.

! Cas des zones de decouvrements (borner la hauteur des vagues):
       do j=0,jmax+1 ; do i=0,imax+1
          hs_wave_t(i,j,1)=min( hs_wave_t(i,j,1),0.5*hz_w(i,j,2)) !18-02-15
       enddo ; enddo
      endif            !-----Michaud2012----> !21-11-19

      do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1 !30-10-20
       velstokes_u(i,j,k,1)=(c1_*velstokes_u(i,j,k,1) &
                            +c2_*velstokes_u(i,j,k,2) &
                                 )*wetmask_u(i,j)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1 !30-10-20
       velstokes_v(i,j,k,1)=(c1_*velstokes_v(i,j,k,1) &
                            +c2_*velstokes_v(i,j,k,2) &
                                 )*wetmask_v(i,j)
      enddo ; enddo ; enddo

      if(iwve==3) then !>>
       do j=0,jmax+1 ; do i=0,imax+1
          hs_wave_t(i,j,1)=c1_*  hs_wave_t(i,j,1)+c2_*  hs_wave_t(i,j,2)
          hs_wave_t(i,j,1)=min( hs_wave_t(i,j,1),0.5*hz_w(i,j,2)) !18-02-15
          t_wave_t(i,j,1)=c1_*   t_wave_t(i,j,1)+c2_*   t_wave_t(i,j,2)
          k_wave_t(i,j,1)=c1_*   k_wave_t(i,j,1)+c2_*   k_wave_t(i,j,2)
       enddo ; enddo
       do j=1,jmax ; do i=1,imax
          dir_wave_t(i,j,1)=atan2(0.5*(velstokes_v(i,j,kmax,1)+velstokes_v(i,j+1,kmax,1)) &
                    ,0.5*(velstokes_u(i,j,kmax,1)+velstokes_u(i+1,j,kmax,1)))
       enddo ; enddo
      endif            !>>

      end subroutine wave_linear_in_time_s26

!..................................................................

      subroutine wave_zaveraged_velstokes
      implicit none
#ifdef synopsis
       subroutinetitle='wave_zaveraged_velstokes'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Depth-averaged Stokes velocity:

! reset en K=1:
      do j=0,jmax+1 !30-10-20
      do i=1,imax+1
        velbarstokes_u(i,j,1)=velstokes_u(i,j,1,1)*dz_u(i,j,1,1) !26-05-10
      enddo
      enddo
      do j=1,jmax+1
      do i=0,imax+1 !30-10-20
        velbarstokes_v(i,j,1)=velstokes_v(i,j,1,1)*dz_v(i,j,1,1)
      enddo
      enddo

! sommation sur K suivants:
      do k=2,kmax
       do j=0,jmax+1 !30-10-20
       do i=1,imax+1
         velbarstokes_u(i,j,1)=velbarstokes_u(i,j,1)      &
                                 +velstokes_u(i,j,k,1)    &
                                        *dz_u(i,j,k,1)

       enddo
       enddo
       do j=1,jmax+1
       do i=0,imax+1 !30-10-20
         velbarstokes_v(i,j,1)=velbarstokes_v(i,j,1)+     &
                                         dz_v(i,j,k,1)    &
                                 *velstokes_v(i,j,k,1)
       enddo
       enddo
      enddo

! Normalisation:
      do j=0,jmax+1 !30-10-20
      do i=1,imax+1
       velbarstokes_u(i,j,1)=velbarstokes_u(i,j,1)     &
                                      /hz_u(i,j,1) !05-06-20

      enddo
      enddo
      do j=1,jmax+1
      do i=0,imax+1 !30-10-20
       velbarstokes_v(i,j,1)=velbarstokes_v(i,j,1)     &
                                      /hz_v(i,j,1)
      enddo
      enddo

      end subroutine wave_zaveraged_velstokes

!..................................................................

      subroutine wave_ww3_grid_resolution
      implicit none
      real*4 :: res_glb_          &
               ,ww3i_,ww3j_       &
               ,factor_=1.   ! =3.
!              ,ww3filval=-999.
      integer ncid_,loop_
#ifdef synopsis
       subroutinetitle='wave_ww3_grid_resolution'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      ww3filval=-999.
      allocate(ww3_res(ww3_imax)) ;   ww3_res=ww3filval

      open(unit=3,file=trim(tmpdirname)//'wavelist')
      read(3,'(a)')texte30
      close(3)
      open(unit=3,file=trim(texte30),access='direct',recl=reclen       &
                                    ,form='unformatted')
      read(3,rec=1)texte80(2),i,x1
      close(3)
      k=index(texte80(2),'/',back=.true.)
      texte80(2)(k+1: )='ww3_resolution.nc'

! critere 2: ne depend pas de la grille s et donc le fichier resolution est dans le repertoire io_waga
!     k=index(texte80(1),'/',back=.true.) ! chemin d'acces au repertoire des fichier oiwaga
!     status=nf_open(texte80(1)(1:k)//'ww3_resolution.nc',nf_nowrite,ncid_)
! critere 1: depend la grille s et donc le fichier resolution est dans le repertoire des notebook
!     k=index(nomfichier(1),'/',back=.true.) ! Critere 1: chemin d'acces au repertoire des notebook
!     status=nf_open(nomfichier(1)(1:k)//'ww3_resolution.nc',nf_nowrite,ncid_)
      status=nf_open(trim(texte80(2)),nf_nowrite,ncid_)


      if(status==0) then !********************************************************************>

       status=nf_inq_varid(ncid_,'resolution',var_id) ;  if(status/=0)stop 'erreur nf_inq_varid resolution'
       varstart(1)=1 ; varcount(1)=ww3_imax
       status=nf_get_vara_real(ncid_,var_id,varstart(1:1)    &
                                           ,varcount(1:1)    &
                                           ,ww3_res(1:ww3_imax)) ;  if(status/=0)stop 'erreur nf_get_vara_real ww3_res'
       status=nf_get_att_real(ncid_,var_id,'_FillValue',ww3filval)
       if(status/=0)stop 'erreur nf_get_att_real _FillValue'
       status=nf_close(ncid_)

! Maintenant on va convertir la distance d'influence de metres en points de grille
! Tout d'abord conversion de la resolution en metre en degres:

      do i=1,ww3_imax
       if(ww3_res(i)/=ww3filval)ww3_res(i)=(ww3_res(i)/rayonterre)*rad2deg
      enddo

! Ensuite conversion de la distance de degres en points de grille
! Stockage des points ww3 qui impactent le proc dans un fichier tmp
      open(unit=3,file=trim(tmpdirname)//'ww3grid_to_sgrid_'//dom_c//'.out')
      do i=1,ww3_imax

! Position du point ww3 sur la grille s:
       call latlontoij(ww3_lon(i)*deg2rad,ww3_lat(i)*deg2rad,'glb')
       ww3i_=deci ; ww3j_=decj

       if(ww3_res(i)/=ww3filval) then !--------->

! Position du point ww3 + delta de latitude correspondant à la distance d'influence du point
! Cette info et l'info precedente va servir à exprimer la distance d'influence en nombre de points
! sur la grille s
         call latlontoij(  ww3_lon(i)*deg2rad  &
                         ,(ww3_lat(i)+ww3_res(i))*deg2rad,'glb')

         ww3_res(i)=sqrt( (deci-ww3i_)**2 + (decj-ww3j_)**2 )*factor_  ! distance d'influence en nbre de points de grille s
         ww3_res(i)=max(ww3_res(i),real(un))


       else                         !--------->

         ww3_res(i)=1.

       endif                        !--------->

       if(     ww3i_+ww3_res(i)-par%timax(1)>=0            &
          .and.ww3i_-ww3_res(i)-par%timax(1)<=imax+1       &
          .and.ww3j_+ww3_res(i)-par%tjmax(1)>=0            &
          .and.ww3j_-ww3_res(i)-par%tjmax(1)<=jmax+1) then !>>>>>
           write(3,*)i,ww3i_,ww3j_,ww3_res(i) ! indices globaux pour conservation mpi
       endif                                  !>>>>>


      enddo
      close(3)

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif


      else               !******************************************************************************>

! Trouver la resolution en prenant la valeur maximum de deux criteres:
! Critere 1: En chaque point de grille S on cherche le point ww3 le plus proche.
! Un point ww3 se voit ainsi désigné comme le point ww3 le plus proche d'un lot de points S.
! A chaque fois, une distance est associé. La résolution de ce point ww3 est la plus
! grande de ces distances.
! Critere 2: En chaque point ww3, on cherche le point ww3 le plus proche et la distance donne la resolution.

! Chercher le point ww3 le plus proche dans 4 cadrans
      allocate(dist_4cadrans(4))
      do i=1,ww3_imax
      if(mod(i,100)==0)write(6,*)i,ww3_imax

       i0=1 ; j0=1
       do loop_=1,4 ! boucle sur 4 cadrans
       i0=-i0 ; j0=j0*i0

       dist_4cadrans(loop_)=1.d20
       do i1=1,ww3_imax
        if(i1/=i) then !XXXXX>
         if((ww3_lat(i1)-ww3_lat(i))*i0>=0.) then !----->
          if((ww3_lon(i1)-ww3_lon(i))*j0>=0.) then !----->

           dist_4cadrans(loop_)=min(dist_4cadrans(loop_),rayonterre    &
           *acos( sin(ww3_lat(i1)*deg2rad)*sin(ww3_lat(i)*deg2rad)     &
                 +cos(ww3_lat(i1)*deg2rad)*cos(ww3_lat(i)*deg2rad)     &
                 *cos(ww3_lon(i1)*deg2rad-ww3_lon(i)*deg2rad)))


          endif                                    !----->
         endif                                    !----->
        endif          !XXXXX>
       enddo ! i1

      enddo  ! loop_

! La plus petite des 4 distances
      dist1=min(min(min( dist_4cadrans(1)   &
                        ,dist_4cadrans(2))  &
                        ,dist_4cadrans(3))  &
                        ,dist_4cadrans(4))
! Retenir la plus grande des 4 distances en ecartant
! les points anormalement loins (liés aux bordures de domaines)
      ww3_res(i)=dist1
      do loop_=1,4
       if(dist_4cadrans(loop_)>ww3_res(i).and.  &
          dist_4cadrans(loop_)<10.*dist1)       &
          ww3_res(i)=dist_4cadrans(loop_)
      enddo

      enddo  ! i
      deallocate(dist_4cadrans)

      if(par%rank==0) then !0000000000000000000>
!     k=index(texte80(1),'/',back=.true.)    ! Critere 2: chemin d'acces au repertoire des fichier oiwaga
!     status=nf_create(texte80(1)(1:k)//'ww3_resolution.nc',nf_clobber,ncid_)
!     k=index(nomfichier(1),'/',back=.true.) ! Critere 1: chemin d'acces au repertoire des notebook
!     status=nf_create(nomfichier(1)(1:k)//'ww3_resolution.nc',nf_clobber,ncid_)
      status=nf_create(trim(texte80(2)),nf_clobber,ncid_)
      status=nf_def_dim(ncid_,'node',ww3_imax,dim_x_id)
      status=nf_def_var(ncid_,'resolution',nf_real,1,dim_x_id,var_id)
      status=nf_put_att_real(ncid_,var_id,'_FillValue',nf_real,1,ww3filval)
      status=nf_enddef(ncid_)
!..............
      status=nf_inq_varid(ncid_,'resolution',var_id)
      status=nf_put_var_real(ncid_,var_id,ww3_res(1:ww3_imax))
      status=nf_close(ncid_)
      endif                !0000000000000000000>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif

      stop 'STOP fichier ww3_resolution.nc produit'

      endif              !**************************************************************************>

      deallocate(ww3_res)
      end subroutine wave_ww3_grid_resolution

!...............................................................................
! WW3 SUR GRILLE SYMPHONIE:
      subroutine wave_grid_extraction_wgrid
      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'
      double precision deci_min_,deci_max_       &
                      ,decj_min_,decj_max_
      integer ncid_
#ifdef synopsis
       subroutinetitle='wave_grid_extraction_wgrid'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      ww3_type_grid=type_structured

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
      if(status.ne.0)stop 'echec fichier'

                   status=nf_inq_dimid(ncid_,'nav_lon',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'x',dim_x_id)
      if(status/=0)stop 'error module_wave nf_inq_dimid dim_x_id'

                   status=nf_inq_dimid(ncid_,'nav_lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'y',dim_y_id)
      if(status/=0)stop 'error module_wave nf_inq_dimid dim_y_id'

      status=nf_inq_dimlen(ncid_,dim_x_id,ww3_imax)
      if(status/=0)stop 'error module_wave nf_inq_dimlen ww3_imax'
      status=nf_inq_dimlen(ncid_,dim_y_id,ww3_jmax)
      if(status/=0)stop 'error module_wave nf_inq_dimlen ww3_jmax'

      if(iperiodicboundary) then !!! grille periodique Oi >>
       if(ww3_imax/=iglb-2)stop 'Err5610 module_wave ww3_imax/=iglb+2 &
        not consitent with ww3_on_sgrid'
      else                       !!! grille periodique Oi >>

       if(ww3_imax/=iglb+2.and.&
          ww3_imax/=iglb)stop 'Err5615 module_wave ww3_imax/=iglb+2 &
        or /=iglb not consitent with ww3_on_sgrid'

      endif                      !!! grille periodique Oi >>

      if(ww3_jmax/=jglb+2.and. &
         ww3_jmax/=jglb)stop 'Err5621 module_wave ww3_jmax/=jglb+2 &
       or /=jglb not consitent with ww3_on_sgrid'

      if(ww3_jmax==jglb+2)ww3_on_sgrid_case=ww3_on_sgrid_casez1 !11-11-21
      if(ww3_jmax==jglb  )ww3_on_sgrid_case=ww3_on_sgrid_casez0 !11-11-21

! lire les longitudes et latitudes:
! Dans le cas où les 2 grilles (ww3 et s) sont identiques on considere imax
! et jmax comme largeur d'intervale de lecture:
      if(ww3_on_sgrid_case==ww3_on_sgrid_casez1) then !/// iglb+2 /// jglb+2 ///>
       if(iperiodicboundary) then !!! grille periodique Oi >>
        ww3_imax=imax-2
       else                       !!! grille periodique Oi >>
        ww3_imax=imax+2
       endif                      !!! grille periodique Oi >>
       ww3_jmax=jmax+2
       ww3zoom_istr=1+par%timax(1)
       ww3zoom_jstr=1+par%tjmax(1)
       ww3zoom_iend=ww3zoom_istr+ww3_imax-1
       ww3zoom_jend=ww3zoom_jstr+ww3_jmax-1
      endif                                           !/// iglb+2 /// jglb+2 ///>
      if(ww3_on_sgrid_case==ww3_on_sgrid_casez0) then !/// iglb   /// jglb   ///> !11-11-21
       ww3_imax=imax+2
       ww3_jmax=jmax+2
! Cas particulier bord i=1, on lit de 1 A imax+1:
       if(par%timax(1)==0)ww3_imax=imax+1
! Cas particulier bord i=iglb, on lit de 0 A imax:
       if(imax+par%timax(1)==iglb)ww3_imax=imax+1
! Cas particulier bord j=1, on lit de 1 A jmax+1:
       if(par%tjmax(1)==0)ww3_jmax=jmax+1
! Cas particulier bord j=jglb, on lit de 0 A jmax:
       if(jmax+par%tjmax(1)==jglb)ww3_jmax=jmax+1

       ww3zoom_istr=0+par%timax(1)
       ww3zoom_jstr=0+par%tjmax(1)
! Cas particulier bord i=1, on lit de 1 A imax+1:
       if(par%timax(1)==0)ww3zoom_istr=1+par%timax(1)
! Cas particulier bord j=1, on lit de 1 A jmax+1:
       if(par%tjmax(1)==0)ww3zoom_jstr=1+par%tjmax(1)

       ww3zoom_iend=ww3zoom_istr+ww3_imax-1
       ww3zoom_jend=ww3zoom_jstr+ww3_jmax-1
      endif                                           !/// iglb   /// jglb   ///>
! On appelle wave_read_lonlat simplement pour verifier que les lon lat des 2
! grilles sont identiques (sinon le modele s'arrete)
      call wave_read_lonlat(ncid_)

! Tableaux de passage d'une grille a l'autre. Dans le cas wgrid il s'agit
! seulement d'une simple egalite...
      call wave_sgrid_to_ww3grid

! On s'assure de la bonne correspondance des grilles des subroutine wave_get_dir, wave_get_gs etc....
! en utilisant les memes relations de passage, appliquees aux longitudes latitudes qui doivent
! etre egales anx lon_t lat_t de S:
      if(ww3_on_sgrid_case==ww3_on_sgrid_casez0) then !z0z0z0>
! LE CAS z0 RESTE A FAIRE !!
      endif                                           !z0z0z0>
      if(ww3_on_sgrid_case==ww3_on_sgrid_casez1) then !z1z1z1>
       i0=0 ; if(iperiodicboundary)i0=2
       do j=0,jmax+1 ; do i=i0,imax+1-i0
       i1=i+1-i0 ; j1=j+1
         if(abs(ww3_2dr4_lon(i1,j1)-lon_t(i,j)*rad2deg)>0.001) then
          write(6,*)i,j,i1,j1,ww3_2dr4_lon(i1,j1),lon_t(i,j)*rad2deg
          stop ' Warning: ww3_2dr4_lon(i1,j1) /= lon_t(i,j)'
         endif
         if(abs(ww3_2dr4_lat(i1,j1)-lat_t(i,j)*rad2deg)>0.001) then
          write(6,*)i,j,i1,j1,ww3_2dr4_lat(i1,j1),lat_t(i,j)*rad2deg
          stop ' Warning: ww3_2dr4_lat(i1,j1) /= lat_t(i,j)'
         endif
       enddo ; enddo
      endif                                           !z1z1z1>

! Desallouer les tableaux
       if(allocated(ww3_lon))deallocate(ww3_lon)
       if(allocated(ww3_lat))deallocate(ww3_lat)
       if(allocated(ww3cg_lon))deallocate(ww3cg_lon)
       if(allocated(ww3cg_lat))deallocate(ww3cg_lat)
       if(allocated(ww3_2dr4_lon))deallocate(ww3_2dr4_lon)
       if(allocated(ww3_2dr4_lat))deallocate(ww3_2dr4_lat)

      status=nf_close(ncid_)

! Si disponible charger la bathy native de WW3
      call wave_get_dpt

      end subroutine wave_grid_extraction_wgrid

!------------------------------------------------------------------------------

      subroutine wave_error(txt_func_,txt_var_,call_code_)
      implicit none
      character(len=*),intent(in) :: txt_func_,txt_var_
      integer call_code_

      write(6,'(a)')'netcdf error status:'
      write(6,'(a,a)')'netcdf function:',txt_func_
      write(6,'(a,a)')'netcdf variable:',txt_var_
      write(6,'(a,i0)')'call code:',call_code_
      write(6,'(a,i0)')'par%rank',par%rank

      stop 'Stop in subroutine wave_error'
      end subroutine wave_error

!------------------------------------------------------------------------------
      subroutine wave_obc_anyvar2d_jeq1
      implicit none
#ifdef synopsis
       subroutinetitle='wave_obc_anyvar2d_jeq1'
       subroutinedescription= &
       'ww3 lateral open boundary conditions'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Apres constat que les champs ww3 n'etaient pas bons de j=0 à j=1 on applique une
! c.l. de gradient nul
      if(obcstatus(jeq1)==1) then !>>>>>>
       do i=0,imax+1
        anyvar2d(i,1)=anyvar2d(i,2)*mask_t(i,2,kmax)  &
                               +(1.-mask_t(i,2,kmax))*anyvar2d(i,1)
       enddo
      endif                       !>>>>>>

      end subroutine wave_obc_anyvar2d_jeq1

!------------------------------------------------------------------------------
#ifdef bidon
      subroutine wave_mpi_anyvar2d
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='wave_mpi_anyvar2d'
       subroutinedescription= &
       'Lateral boundary conditions ensuring mpi continuity'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
      call get_type_echange('za','anyvar2d_za' &
                                 ,anyvar2d     &
                          ,lbound(anyvar2d)    &
                          ,ubound(anyvar2d)    &
                                       ,k0)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(anyvar2d,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine wave_mpi_anyvar2d
#endif
!------------------------------------------------------------------------------

       subroutine wave_ww3onsgrid_checkmask
       implicit none
#ifdef synopsis
       subroutinetitle='wave_ww3onsgrid_checkmask'
       subroutinedescription= &
       'check that the ww3 mask is the same than in S26'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Dans le cas ww3_on_sgrid les masques des 2 grilles doivent et identique
         i0=0 ; if(iperiodicboundary)i0=2
         do j=0,jmax+1 ; do i=i0,imax+1-i0
          if(ww3_var(i+1-i0,j+1)==ww3filval.and.mask_t(i,j,kmax)==1)then!>>>>
           write(6,*)'ww3 mask and S26 mask are inconsistent'
           stop 'Stop wave_ww3onsgrid_checkmask'
          endif                                                         !>>>>
         enddo         ; enddo

       end subroutine wave_ww3onsgrid_checkmask

!------------------------------------------------------------------------------

      subroutine wave_notime_interp_cpl
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='wave_notime_interp_cpl'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do j=0,jmax+1
      do i=0,imax+1
        t_wave_t(i,j,1)=     t_wave_t(i,j,2)
        j_wave_t(i,j,1)=    j_wave_t(i,j,2)
         hs_wave_t(i,j,1)=  hs_wave_t(i,j,2)
        hsw_wave_t(i,j,1)= hsw_wave_t(i,j,2)
        foc_wave_t(i,j,1)= foc_wave_t(i,j,2)
        k_wave_t(i,j,1) =   k_wave_t(i,j,2)
        kx_wave_t(i,j,1)=  kx_wave_t(i,j,2)
        ky_wave_t(i,j,1)=  ky_wave_t(i,j,2)
       twox_wave_t(i,j,1)=twox_wave_t(i,j,2)
       twoy_wave_t(i,j,1)=twoy_wave_t(i,j,2)
       tawx_wave_t(i,j,1)=tawx_wave_t(i,j,2)
       tawx_wave_t(i,j,1)=tawx_wave_t(i,j,2)
       uss_wave_t(i,j,1)=uss_wave_t(i,j,2)
       vss_wave_t(i,j,1)=vss_wave_t(i,j,2)
      enddo
      enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
       velstokes_u(i,j,k,1)=velstokes_u(i,j,k,2)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
       velstokes_v(i,j,k,1)=velstokes_v(i,j,k,2)
      enddo ; enddo ; enddo

      do j=0,jmax+1
      do i=0,imax+1
        dir_wave_t(i,j,1)=atan2(ky_wave_t(i,j,1),kx_wave_t(i,j,1))
      
      enddo
      enddo

      end subroutine wave_notime_interp_cpl

!------------------------------------------------------------------------------

      subroutine wave_linear_interp_cpl
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='wave_linear_interp_cpl'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      alpha=1./(1.+max(0,wave_cpl_period_iter-wave_cpl_nextrec))
!      alpha=1.
      do j=0,jmax+1
      do i=0,imax+1
         t_wave_t(i,j,1)=alpha*t_wave_t(i,j,2)+(1.-alpha)*t_wave_t(i,j,1)
         j_wave_t(i,j,1)=alpha*j_wave_t(i,j,2)+(1.-alpha)*j_wave_t(i,j,1)
         hs_wave_t(i,j,1)=alpha*hs_wave_t(i,j,2)+(1.-alpha)*hs_wave_t(i,j,1)
         hsw_wave_t(i,j,1)=alpha*hsw_wave_t(i,j,2)+(1.-alpha)*hsw_wave_t(i,j,1)
         foc_wave_t(i,j,1)=alpha*foc_wave_t(i,j,2)+(1.-alpha)*foc_wave_t(i,j,1)
         k_wave_t(i,j,1)=alpha*k_wave_t(i,j,2)+(1.-alpha)*k_wave_t(i,j,1)
         kx_wave_t(i,j,1)=alpha*kx_wave_t(i,j,2)+(1.-alpha)*kx_wave_t(i,j,1)
         ky_wave_t(i,j,1)=alpha*ky_wave_t(i,j,2)+(1.-alpha)*ky_wave_t(i,j,1)
         twox_wave_t(i,j,1)=alpha*twox_wave_t(i,j,2)+(1.-alpha)*twox_wave_t(i,j,1)
         twoy_wave_t(i,j,1)=alpha*twoy_wave_t(i,j,2)+(1.-alpha)*twoy_wave_t(i,j,1)
         tawx_wave_t(i,j,1)=alpha*tawx_wave_t(i,j,2)+(1.-alpha)*tawx_wave_t(i,j,1)
         tawy_wave_t(i,j,1)=alpha*tawy_wave_t(i,j,2)+(1.-alpha)*tawy_wave_t(i,j,1)
         uss_wave_t(i,j,1)=alpha*uss_wave_t(i,j,2)+(1.-alpha)*uss_wave_t(i,j,1)
         vss_wave_t(i,j,1)=alpha*vss_wave_t(i,j,2)+(1.-alpha)*vss_wave_t(i,j,1)

!         if(i+par%timax(1)==136.and.j+par%tjmax(1)==119) then
!          write(69,*)'wave_cpl_nextrec=',wave_cpl_nextrec
!          write(69,*)'wave_cpl_period_iter=',wave_cpl_period_iter
!          write(69,*)'alpha=',alpha
!          write(69,*)'uss_wave_t(i,j,1)=',uss_wave_t(i,j,1)
!          write(69,*)'uss_wave_t(i,j,2)=',uss_wave_t(i,j,2)
!          write(69,*)'------------------'
!          write(66,*)uss_wave_t(i,j,1)
!         endif
      enddo
      enddo

      do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1 !30-10-20
       velstokes_u(i,j,k,1)=alpha*velstokes_u(i,j,k,2)+(1.-alpha)*velstokes_u(i,j,k,1)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1 !30-10-20
       velstokes_v(i,j,k,1)=alpha*velstokes_v(i,j,k,2)+(1.-alpha)*velstokes_v(i,j,k,1)
      enddo ; enddo ; enddo

      do j=0,jmax+1
      do i=0,imax+1
        dir_wave_t(i,j,1)=atan2(ky_wave_t(i,j,1),kx_wave_t(i,j,1))
      enddo
      enddo

      end subroutine wave_linear_interp_cpl

!------------------------------------------------------------------------------

      subroutine wave_post_reception_traitement
      use module_principal 
      implicit none
#ifdef synopsis
       subroutinetitle='wave_post_reception_traitement'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      ! 1. On applique le wetmask, les seuils et les rotations
      x1=3. ! Seuil max (m/s) appliquE au courant renvoyE au modele de vagues
      do j=0,jmax+1
      do i=0,imax+1
         if (mask_wave_t(i,j) == 0 .or. mask_wave_t(i,j) == -1) then
           t_wave_t(i,j,2)=20.0
           j_wave_t(i,j,2)=0.
           hs_wave_t(i,j,2)=0.
           hsw_wave_t(i,j,2)=0.
           foc_wave_t(i,j,2)=0.
           twox_wave_t(i,j,2)=0.
           twoy_wave_t(i,j,2)=0.
           tawx_wave_t(i,j,2)=0.
           tawy_wave_t(i,j,2)=0.
           uss_wave_t(i,j,2)=0.
           vss_wave_t(i,j,2)=0.
           dir_wave_t(i,j,1)=0.
         endif
     
         ! on recalcule la direction
         !dir = atan2(sin,cos)
         dir_wave_t(i,j,1)=atan2(xy_t(i,j,6), xy_t(i,j,5))
        
         ! S grid rotation
         dir_wave_t(i,j,1)=dir_wave_t(i,j,1)+atan2(gridrotsin_t(i,j),gridrotcos_t(i,j))

         xy_t(i,j,5)=gridrotcos_t(i,j)*tawx_wave_t(i,j,2)-gridrotsin_t(i,j)*tawy_wave_t(i,j,2)
         xy_t(i,j,6)=gridrotsin_t(i,j)*tawx_wave_t(i,j,2)+gridrotcos_t(i,j)*tawy_wave_t(i,j,2)
         tawx_wave_t(i,j,2) = xy_t(i,j,5)*rho
         tawy_wave_t(i,j,2) = xy_t(i,j,6)*rho

         xy_t(i,j,5)=gridrotcos_t(i,j)*twox_wave_t(i,j,2)-gridrotsin_t(i,j)*twoy_wave_t(i,j,2)
         xy_t(i,j,6)=gridrotsin_t(i,j)*twox_wave_t(i,j,2)+gridrotcos_t(i,j)*twoy_wave_t(i,j,2)
         twox_wave_t(i,j,2) = xy_t(i,j,5)*rho
         twoy_wave_t(i,j,2) = xy_t(i,j,6)*rho

         xy_t(i,j,5)=gridrotcos_t(i,j)*uss_wave_t(i,j,2)-gridrotsin_t(i,j)*vss_wave_t(i,j,2)
         xy_t(i,j,6)=gridrotsin_t(i,j)*uss_wave_t(i,j,2)+gridrotcos_t(i,j)*vss_wave_t(i,j,2)
         uss_wave_t(i,j,2) = xy_t(i,j,5)
         vss_wave_t(i,j,2) = xy_t(i,j,6)
     
!       uss_wave_t(i,j,1)=max(-0.9998,min(0.9998,uss_wave_t(i,j,1)))
!       vss_wave_t(i,j,1)=max(-0.9998,min(0.9998,vss_wave_t(i,j,1)))
        ! Seuiller à 3m/s:
        x2=min(1.,x1/sqrt(small1+uss_wave_t(i,j,2)**2+vss_wave_t(i,j,2)**2))
        uss_wave_t(i,j,2)=x2*uss_wave_t(i,j,2)
        vss_wave_t(i,j,2)=x2*vss_wave_t(i,j,2)

        ! Seuil HSW
        hsw_wave_t(i,j,2)=max(hsw_wave_t(i,j,2),1e-2) !03-02-15
      enddo
      enddo

      ! 2. On corrige les frontières ouvertes 
      !!==
      !! Les champs sont réellement calculés entre 2 et imax-1
      !! Il a été choisi de dubliquer les derniers points de grille sur les
      !! conditions aux limites. On applique par exemple sur le bord Est :
      !! imax = imax-1
      !! imax+1 = imax-1
      !! et pareil pour les autres bords sans tenir compte des coins
      !!==
      
      if(obcstatus(jeqjmax)==1) then
      ! j=jmax du processeur courant est la frontière Nord:
         do i=0,imax+1
           hs_wave_t(i,jmax,2)=hs_wave_t(i,jmax-1,2)
           hs_wave_t(i,jmax+1,2)=hs_wave_t(i,jmax-1,2)
     
           hsw_wave_t(i,jmax,2)=hsw_wave_t(i,jmax-1,2)
           hsw_wave_t(i,jmax+1,2)=hsw_wave_t(i,jmax-1,2)
     
           foc_wave_t(i,jmax,2)=foc_wave_t(i,jmax-1,2)
           foc_wave_t(i,jmax+1,2)=foc_wave_t(i,jmax-1,2)
     
           t_wave_t(i,jmax,2)=t_wave_t(i,jmax-1,2)
           t_wave_t(i,jmax+1,2)=t_wave_t(i,jmax-1,2)
     
           tawx_wave_t(i,jmax,2)=tawx_wave_t(i,jmax-1,2)
           tawx_wave_t(i,jmax+1,2)=tawx_wave_t(i,jmax-1,2)
           tawy_wave_t(i,jmax,2)=tawy_wave_t(i,jmax-1,2)
           tawy_wave_t(i,jmax+1,2)=tawy_wave_t(i,jmax-1,2)
     
           twox_wave_t(i,jmax,2)=twox_wave_t(i,jmax-1,2)
           twox_wave_t(i,jmax+1,2)=twox_wave_t(i,jmax-1,2)
           twoy_wave_t(i,jmax,2)=twoy_wave_t(i,jmax-1,2)
           twoy_wave_t(i,jmax+1,2)=twoy_wave_t(i,jmax-1,2)
     
           uss_wave_t(i,jmax,2)=uss_wave_t(i,jmax-1,2)
           uss_wave_t(i,jmax+1,2)=uss_wave_t(i,jmax-1,2)
           vss_wave_t(i,jmax,2)=vss_wave_t(i,jmax-1,2)
           vss_wave_t(i,jmax+1,2)=vss_wave_t(i,jmax-1,2)
     
           dir_wave_t(i,jmax,1)=dir_wave_t(i,jmax-1,1)
           dir_wave_t(i,jmax+1,1)=dir_wave_t(i,jmax-1,1)
     
           mask_wave_t(i,jmax)=mask_wave_t(i,jmax-1)
           mask_wave_t(i,jmax+1)=mask_wave_t(i,jmax-1)
         enddo
      endif
      if(obcstatus(jeq1)==1) then
      ! j=1 du processeur courant est la frontière Sud:
         do i=0,imax+1
           hs_wave_t(i,1,2)=hs_wave_t(i,2,2)
           hs_wave_t(i,0,2)=hs_wave_t(i,2,2)
     
           hsw_wave_t(i,1,2)=hsw_wave_t(i,2,2)
           hsw_wave_t(i,0,2)=hsw_wave_t(i,2,2)
     
           foc_wave_t(i,1,2)=foc_wave_t(i,2,2)
           foc_wave_t(i,0,2)=foc_wave_t(i,2,2)
     
           t_wave_t(i,1,2)=t_wave_t(i,2,2)
           t_wave_t(i,0,2)=t_wave_t(i,2,2)
     
           tawx_wave_t(i,1,2)=tawx_wave_t(i,2,2)
           tawx_wave_t(i,0,2)=tawx_wave_t(i,2,2)
           tawy_wave_t(i,1,2)=tawy_wave_t(i,2,2)
           tawy_wave_t(i,0,2)=tawy_wave_t(i,2,2)
     
           twox_wave_t(i,1,2)=twox_wave_t(i,2,2)
           twox_wave_t(i,0,2)=twox_wave_t(i,2,2)
           twoy_wave_t(i,1,2)=twoy_wave_t(i,2,2)
           twoy_wave_t(i,0,2)=twoy_wave_t(i,2,2)
     
           uss_wave_t(i,1,2)=uss_wave_t(i,2,2)
           uss_wave_t(i,0,2)=uss_wave_t(i,2,2)
           vss_wave_t(i,1,2)=vss_wave_t(i,2,2)
           vss_wave_t(i,0,2)=vss_wave_t(i,2,2)
     
           dir_wave_t(i,1,1)=dir_wave_t(i,2,1)
           dir_wave_t(i,0,1)=dir_wave_t(i,2,1)
     
           mask_wave_t(i,1)=mask_wave_t(i,2)
           mask_wave_t(i,0)=mask_wave_t(i,2)
         enddo
      endif
      if(obcstatus(ieqimax)==1) then
      ! i=imax du processeur courant est la frontière Est:
         do j=0,jmax+1
           hs_wave_t(imax,j,2)=hs_wave_t(imax-1,j,2)
           hs_wave_t(imax+1,j,2)=hs_wave_t(imax-1,j,2)
     
           hsw_wave_t(imax,j,2)=hsw_wave_t(imax-1,j,2)
           hsw_wave_t(imax+1,j,2)=hsw_wave_t(imax-1,j,2)
     
           foc_wave_t(imax,j,2)=foc_wave_t(imax-1,j,2)
           foc_wave_t(imax+1,j,2)=foc_wave_t(imax-1,j,2)
     
           t_wave_t(imax,j,2)=t_wave_t(imax-1,j,2)
           t_wave_t(imax+1,j,2)=t_wave_t(imax-1,j,2)
     
           tawx_wave_t(imax,j,2)=tawx_wave_t(imax-1,j,2)
           tawx_wave_t(imax+1,j,2)=tawx_wave_t(imax-1,j,2)
           tawy_wave_t(imax,j,2)=tawy_wave_t(imax-1,j,2)
           tawy_wave_t(imax+1,j,2)=tawy_wave_t(imax-1,j,2)
     
           twox_wave_t(imax,j,2)=twox_wave_t(imax-1,j,2)
           twox_wave_t(imax+1,j,2)=twox_wave_t(imax-1,j,2)
           twoy_wave_t(imax,j,2)=twoy_wave_t(imax-1,j,2)
           twoy_wave_t(imax+1,j,2)=twoy_wave_t(imax-1,j,2)
     
           uss_wave_t(imax,j,2)=uss_wave_t(imax-1,j,2)
           uss_wave_t(imax+1,j,2)=uss_wave_t(imax-1,j,2)
           vss_wave_t(imax,j,2)=vss_wave_t(imax-1,j,2)
           vss_wave_t(imax+1,j,2)=vss_wave_t(imax-1,j,2)
     
           dir_wave_t(imax,j,1)=dir_wave_t(imax-1,j,1)
           dir_wave_t(imax+1,j,1)=dir_wave_t(imax-1,j,1)
     
           mask_wave_t(imax,j)=mask_wave_t(imax-1,j)
           mask_wave_t(imax+1,j)=mask_wave_t(imax-1,j)
         enddo
      endif
      if (obcstatus(ieq1)==1) then
      ! i=1 du processeur courant est la frontière West:
         do j=0,jmax+1
           hs_wave_t(1,j,2)=hs_wave_t(2,j,2)
           hs_wave_t(0,j,2)=hs_wave_t(2,j,2)
     
           hsw_wave_t(1,j,2)=hsw_wave_t(2,j,2)
           hsw_wave_t(0,j,2)=hsw_wave_t(2,j,2)
     
           foc_wave_t(1,j,2)=foc_wave_t(2,j,2)
           foc_wave_t(0,j,2)=foc_wave_t(2,j,2)
     
           t_wave_t(1,j,2)=t_wave_t(2,j,2)
           t_wave_t(0,j,2)=t_wave_t(2,j,2)
     
           tawx_wave_t(1,j,2)=tawx_wave_t(2,j,2)
           tawx_wave_t(0,j,2)=tawx_wave_t(2,j,2)
           tawy_wave_t(1,j,2)=tawy_wave_t(2,j,2)
           tawy_wave_t(0,j,2)=tawy_wave_t(2,j,2)
     
           twox_wave_t(1,j,2)=twox_wave_t(2,j,2)
           twox_wave_t(0,j,2)=twox_wave_t(2,j,2)
           twoy_wave_t(1,j,2)=twoy_wave_t(2,j,2)
           twoy_wave_t(0,j,2)=twoy_wave_t(2,j,2)
     
           uss_wave_t(1,j,2)=uss_wave_t(2,j,2)
           uss_wave_t(0,j,2)=uss_wave_t(2,j,2)
           vss_wave_t(1,j,2)=vss_wave_t(2,j,2)
           vss_wave_t(0,j,2)=vss_wave_t(2,j,2)
     
           dir_wave_t(1,j,1)=dir_wave_t(2,j,1)
           dir_wave_t(0,j,1)=dir_wave_t(2,j,1)
     
           mask_wave_t(1,j)=mask_wave_t(2,j)
           mask_wave_t(0,j)=mask_wave_t(2,j)
         enddo
      endif

      ! 3. On fait les échanges MPI
      call get_type_echange('za','hs_wave_t_za_2',hs_wave_t,lbound(hs_wave_t),ubound(hs_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(hs_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','hsw_wave_t_za_2',hsw_wave_t,lbound(hsw_wave_t),ubound(hsw_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(hsw_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','dir_wave_t_za_1',dir_wave_t,lbound(dir_wave_t),ubound(dir_wave_t),1,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(dir_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','foc_wave_t_za_2',foc_wave_t,lbound(foc_wave_t),ubound(foc_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(foc_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','t_wave_t_za_2',t_wave_t,lbound(t_wave_t),ubound(t_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(t_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','tawx_wave_t_za_2',tawx_wave_t,lbound(tawx_wave_t),ubound(tawx_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(tawx_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      
      call get_type_echange('za','tawy_wave_t_za_2',tawy_wave_t,lbound(tawy_wave_t),ubound(tawy_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(tawy_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','twox_wave_t_za_2',twox_wave_t,lbound(twox_wave_t),ubound(twox_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(twox_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      
      call get_type_echange('za','twoy_wave_t_za_2',twoy_wave_t,lbound(twoy_wave_t),ubound(twoy_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(twoy_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','uss_wave_t_za_1',uss_wave_t,lbound(uss_wave_t),ubound(uss_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(uss_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','vss_wave_t_za_1',vss_wave_t,lbound(vss_wave_t),ubound(vss_wave_t),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(vss_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      call get_type_echange('za','mask_wave_t_za_1',mask_wave_t,lbound(mask_wave_t),ubound(mask_wave_t),i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(mask_wave_t,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()

      end subroutine wave_post_reception_traitement

!--------------------------------------------------------------------------------
#ifdef bidon
      subroutine wave_ww3_spec_sete_read !25-03-19
      use module_parallele
      implicit none

! Ce programme ne sert que pour les frontieres. Ne pas y entrer quand ce n'est pas utile
! permet d'eviter d'allouer de la memoire et surcharger la memoire du noeud (les
! coeurs de ce mEme noeud non concernEs n'allouent pas la memoire)
!     if(obcstatus(ieqimax)         &       
!       +obcstatus(ieq1)            &
!       +obcstatus(jeqjmax)         &
!       +obcstatus(jeq1)            &
!       +obcstatus(ieq1_jeq1)       &
!       +obcstatus(ieqimax_jeq1)    &
!       +obcstatus(ieqimax_jeqjmax) &
!       +obcstatus(ieq1_jeqjmax)/=0)flag_allocate_=1

!#ifdef bidon
      if( &
         obcstatus(ieqimax)         &       
        +obcstatus(ieq1)            &
        +obcstatus(jeqjmax)         &
        +obcstatus(jeq1)            &
!       +obcstatus(ieq1_jeq1)       &
!       +obcstatus(ieqimax_jeq1)    &
!       +obcstatus(ieqimax_jeqjmax) &
!       +obcstatus(ieq1_jeqjmax)    &
                                /=0)flag_allocate_=1
!#endif

!     write(6,*)'BIDOUILLE FLAG_ALLOCATE'
!     flag_allocate_=1

!      write(10+par%rank,*)'coco5 ',flag_allocate_


! Spectre complet pour simulations realistes
!     dirw_beg=4  ! 10 ! 4
!     dirw_end=15 ! 10 ! 15
!     freq_beg=1  ! 6  ! 1  
!     freq_end=13 ! 6  ! 13


! Spectre simplifie pour tests academiques
!     dirw_beg=13  ! Direction de i=imax vers i=1
!     dirw_end=13     
!     dirw_beg=10
!     dirw_end=10     

!     dirw_beg=16 ! avant le 01/07/2019
!     dirw_end=16     
      dirw_beg=17 ! avant le 01/07/2019
      dirw_end=17     
!     dirw_beg=18 ! avant le 01/07/2019
!     dirw_end=18     
!     dirw_beg=19 ! avant le 01/07/2019
!     dirw_end=19     
!     dirw_beg=21 ! avant le 01/07/2019
!     dirw_end=21     
      freq_beg=8  
      freq_end=8  

!     allocate(cumul1_mpi       (freq_beg:freq_end,dirw_beg:dirw_end)) ; cumul1_mpi=0. 
!     allocate(cumul2_mpi       (freq_beg:freq_end,dirw_beg:dirw_end)) ; cumul2_mpi=0. 

      texte90= &
      '../../../SETE/WW3_SPECTRAL_SETE_HELOISE/ww3.point36_2009_spec.nc'
!     '/home/marp/HELOISE/WW3_SPECTRAL_SETE_HELOISE/ww3.point36_2009_spec.nc'

      status=nf_open(trim(texte90),nf_nowrite,ncid)
      if(status/=0)stop 'Err nf_open wave_ww3_spec_sete_read'

! Station (on sait qu'il n'y a qu'une station...)
!     status=nf_inq_dimid(ncid,'station',dim_x_id)
!     if(status/=0)stop 'Err nf_inq_dimid station'
!     status=nf_inq_dimlen(ncid,dim_x_id,stationmax)
!     if(status/=0)stop 'nf_inq_dimlen stationmax'
!     write(6,*)'stationmax',stationmax

! dim time:
      status=nf_inq_dimid(ncid,'time',dim_z_id)
      if(status/=0)stop 'Err nf_inq_dimid time'
      status=nf_inq_dimlen(ncid,dim_z_id,timemax)
      if(status/=0)stop 'nf_inq_dimlen timemax'

! dim Frequences:
      status=nf_inq_dimid(ncid,'frequency',dim_y_id)
      if(status/=0)stop 'Err nf_inq_dimid frequency'
      status=nf_inq_dimlen(ncid,dim_y_id,freqmax)
      if(status/=0)stop 'nf_inq_dimlen'

! dim Directions:
      status=nf_inq_dimid(ncid,'direction',dim_x_id)
      if(status/=0)stop 'Err nf_inq_dimid direction'
      status=nf_inq_dimlen(ncid,dim_x_id,dirmax)
      if(status/=0)stop 'nf_inq_dimlen'

! Lire les frequences:
      allocate(ww3puls     (freqmax)) ; ww3puls=0.
      allocate(ww3deltafreq(freqmax)) ; ww3deltafreq=0.

      status=nf_inq_varid(ncid,'frequency',var_id)
      if(status/=0)stop 'Err nf_inq_varid'

      status=nf_inq_var(ncid,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'Err nf_inq_var'

!     write(6,*)'var_id=',var_id
!     write(6,*)'var_dims=',var_dims
!     write(6,*)'tabdim',tabdim
!     write(6,*)'i4',i4
!     write(6,*)'var_type=',var_type
!     if(var_type/=nf_real) stop 'var_type/=nf_real'
!     write(6,'(a,a)')'texte30=',trim(texte30)

      varstart(1:var_dims)=1
      varcount(1:var_dims)=freqmax
      status=nf_get_vara_real(ncid,var_id               &
                                  ,varstart(1:var_dims) &
                                  ,varcount(1:var_dims) &
                                   ,ww3puls(1:freqmax))
      if(status/=0)stop 'Err nf_get_vara_real'

! ww3deltafreq (si calcul placE avant la multiplication de ww3puls par 2pi alors ww3deltafreq est en s**-1) !27-03-19
      do freq=2,freqmax-1
       ww3deltafreq(freq)=0.5*(ww3puls(freq+1)-ww3puls(freq-1))
      enddo
      ww3deltafreq(1)      =ww3deltafreq(2)
      ww3deltafreq(freqmax)=ww3deltafreq(freqmax-1)

! on multiplie par 2pi pour passer de la frequence A la pulsation
      ww3puls(1:freqmax)=2.*pi*ww3puls(1:freqmax)

! Lire les directions: (attention aux conventions de direction, details dans:
! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p
      allocate(ww3dir(dirmax)) ; ww3dir=0.

      status=nf_inq_varid(ncid,'direction',var_id)
      if(status/=0)stop 'Err nf_inq_varid'

      status=nf_inq_var(ncid,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'Err nf_inq_var'

      varstart(1:var_dims)=1
      varcount(1:var_dims)=dirmax
      status=nf_get_vara_real(ncid,var_id               &
                                  ,varstart(1:var_dims) &
                                  ,varcount(1:var_dims) &
                                    ,ww3dir(1:dirmax))
      if(status/=0)stop 'Err nf_get_vara_real'

! L'angle archivE dans le fichier WW3 correspond au code de couleur magenta dans le doc suivant:
! https://docs.google.com/presentation/d/1LR_x8zJ6BpWutOMB634IthGejyfwlUZWegepkrGw7qI/edit#slide=id.p
! Il indique la destination dans la convention meteo. On le transforme en une destination dans la
! convention trigonometrique
      ww3dir=modulo(270-ww3dir+180.,360.)
! Passer des degres aux radians
      ww3dir=ww3dir*deg2rad


      ww3deltadir=ww3dir(2)-ww3dir(1)
      if(ww3deltadir<-pi)ww3deltadir=ww3deltadir+2.*pi
      if(ww3deltadir> pi)ww3deltadir=ww3deltadir-2.*pi

!     write(6,*)'BIDOUILLE PATRICK POUR EXPERIENCES SIMPLES'
!     ww3dir=-grid_angle_t(imax/2,jmax/2)-0.45*pi

!     write(6,*)'direction',direction


      if(par%rank==0) then
         write(6,*)'dirmax,freqmax,timemax',dirmax,freqmax,timemax
         write(6,*)'ww3deltadir*rad2deg=',ww3deltadir*rad2deg
         do freq=1,freqmax
          write(6,*)'w,T',freq,real(ww3puls(freq)/(2.*pi)) &
!                            ,real(0.0625*(1.1)**(freq-1)) &
!                            ,real(ww3deltafreq(freq))     &
                             ,real(2.*pi/ww3puls(freq))
         enddo
         do i=1,dirmax
          write(6,*)'i,ww3dir(i)',i,ww3dir(i)*rad2deg
         enddo
      endif

! Lire efth : sea surface wave directional variance spectral density
      allocate(ww3efth(dirmax,freqmax,1,1)) ; ww3efth=0.
      status=nf_inq_varid(ncid,'efth',var_id)
      if(status/=0)stop 'Err nf_inq_varid efth'
      varstart(1)=1 ; varcount(1)=dirmax
      varstart(2)=1 ; varcount(2)=freqmax
      varstart(3)=1 ; varcount(3)=1
      varstart(4)=91 ! une echeance au pif....
      varcount(4)=1
      status=nf_get_vara_real(ncid,var_id               &
                                  ,varstart(1:var_dims) &
                                  ,varcount(1:var_dims) &
                       ,ww3efth(1:dirmax,1:freqmax,1,1))
      if(status/=0)stop 'Err nf_get_vara_real efth'


! Module du vecteur d'onde par frequence le long de la frontiere j=jglb...
!     allocate(ww3kvector(0:imax+1,freqmax))   ; ww3kvector=0.
      allocate(ww3kvector(int(hmax)+1,freqmax))   ; ww3kvector=0. !03-04-19

      
      do freq=1,freqmax
!     do i=0,imax+1
      do i=1,int(hmax)+1 !03-04-19

       j=jmax
       x1=real(i) ! x1 est une bathy
!      x1=h_w(i,j)
!      x1=hmax
       x2=2*pi/ww3puls(freq) ! periode = 2pi/pulsation

!      call wave_freq2kvector(x1,x2,x3) ! x3 kvector
!      x10=x3
       call wave_freq2kvector_acm(x1,x2,x3) ! x3 kvector

       ww3kvector(i,freq)=x3  !&
                   !      /sqrt(1.-((dxb*x3)**2)/12.) ! Terme correction erreur time stepping Eq31 Marsaleix et al 2019 

!      if(par%rank==0)write(666,*)real(x1),real(x2),real(2*pi/x10),real(2*pi/x3)
!      if(par%rank==0)write(666,*)real(x1),real(x2),real(2*pi/x3),real(sqrt(1.-((dxb*x3)**2)/12.))

!      if(i+par%timax(1)==iglb/2.and.j+par%tjmax(1)==jglb) then !>>>
!       write(10+par%rank,*)'h periode longueur',real(x1),real(x2),real(2*pi/x3)
!       write(10+par%rank,*)'h p k L',real(x1),real(x2),ww3kvector(i,freq),real(2*pi/x3)
!      endif                                                    !>>>

      enddo
      enddo

      status=nf_close(ncid)

    
      if(flag_allocate_==1) then !ooooo>
        allocate(ww3phase(0:imax+1,0:jmax+1,freq_beg:freq_end,dirw_beg:dirw_end)) ; ww3phase=0
      endif                      !ooooo>

!     allocate(ww3kvec_oj(0:imax+1,freqmax,dirmax)) ; ww3kvec_oj=0

! Dephaser en i=0 avec une fonction aleatoire
!     x1_r4=0.
!     do dirw=dirw_beg,dirw_end
!     do freq=freq_beg,freq_end
!      call random_number(x1_r4)
!      cumul2_mpi(freq,dirw)=x1_r4*2.*pi
!     enddo
!     enddo


      call wave_phase_jeqjmax_increasing_i
      call wave_phase_decreasing_j
      deallocate(ar_send_sud)
      deallocate(ar_recv_nord)
      deallocate(ar_send_est)
      deallocate(ar_recv_ouest)

!     sshobc_w=0.
!     if(allocated(ww3phase)) then
!      do j=0,jmax+1
!      do i=0,imax+1
!        sshobc_w(i,j,1)=ww3phase(i,j,freq_beg,dirw_beg)*rad2deg
!        if(j==jmax/2)write(10+par%rank,*)i+par%timax(1),ww3phase(i,j,freq_beg,dirw_beg)
!      enddo
!      enddo
!     endif
!     call graph_out

!     stop 'coco3'

!     deallocate(cumul1_mpi)
!     deallocate(cumul2_mpi)

#ifdef bidon
      j=jmax
      do dirw=1,dirmax 
      do freq=1,freqmax
      do i=1,imax
! Koj=gridrotsin*Keastward+gridrotcos*Knorthward
!    =gridrotsin*K*cos(dir)+gridrotcos*K*sin(dir)
       ww3kvec_oj(i,freq,dirw)=                                 &
       ww3kvector(i,freq)*( gridrotsin_t(i,j)*cos(ww3dir(dirw)) &
                           +gridrotcos_t(i,j)*sin(ww3dir(dirw)))

      enddo
      enddo
      enddo
#endif
#ifdef bidon
      if(par%rank==0) then !>>>
       sum1=0. ; sum2=0
       do freq=1,freqmax
        do dirw=1,dirmax
          sum1=sum1+ww3deltafreq(freq)*ww3deltadir*sqrt(2.*ww3efth(dirw,freq,1,1))
          write(666,*)real(ww3dir(dirw)*rad2deg),real(2.*pi/ww3puls(freq)),ww3efth(dirw,freq,1,1)
!         if(dirw>=4.and.dirw<=15.and.freq<=13) then 
          if(freq<=13) then 
           sum2=sum2+ww3deltafreq(freq)*ww3deltadir*sqrt(2.*ww3efth(dirw,freq,1,1))
          endif
        enddo
       enddo
       write(666,*)'sum1 sum2=',sum1,sum2
      endif                !>>>
#endif

! IMPORTANT POUR LA SUITE
! Periode principale: definir cwavepeak, kvectorpeak, periodpeak pour H=hmax
       x1=0. ; x2=7.
!      do freq=1,freqmax ; do dirw=1,dirmax
       do freq=freq_beg,freq_end ; do dirw=dirw_beg,dirw_end
        if(ww3deltafreq(freq)*ww3deltadir*ww3efth(dirw,freq,1,1)>x1) then
         x1=ww3deltafreq(freq)*ww3deltadir*ww3efth(dirw,freq,1,1) 
         x2=2.*pi/ww3puls(freq)
        endif
       enddo ; enddo
       call mpi_allreduce(x2,x3,1,mpi_double_precision,mpi_max,par%comm2d ,ierr)
       periodpeak=x3
       freq2pipeak=2.*pi/x3

! Pour continuite mpi kvectorpeak est calcule avec la meme valeur de h
! (hmax) pour tous les rank
       x1=hmax ; x2=periodpeak ! periode = 2pi/pulsation
       call wave_freq2kvector(x1,x2,x3) ! x3 kvector
       kvectorpeak=x3

! x0 est mu=sqrt( 1-(c/alpha)**2 ) equation 21 dans Marsaleix et al, 2019
       x0=sqrt(1-(cwavepeak/acm_speed)**2 )
! Formule theorique (eq 28, Marsaleix et al, 2019) pour le parametre "r" en facteur du PGF non hydrostatique 
!      nhpgf_reduce=(x0**3)*(tanh(    kvectorpeak*hmax ) -    kvectorpeak*hmax) &
!                          /(tanh( x0*kvectorpeak*hmax ) - x0*kvectorpeak*hmax)



       cwavepeak=2.*pi/(kvectorpeak*periodpeak)


      k0=9999 ! sert A detecter le plus petit rank concernE par la frontiere
      if(flag_allocate_==1)k0=par%rank ! sert A detecter le plus petit rank concernE par les OBC pour diffusion d'info...
!     endif                               !avril> !03-04-19

! k1 est le plus petit rank concernE par la frontiere
      call mpi_allreduce(k0,k1,1,mpi_integer,mpi_min,par%comm2d,ierr)

      if(par%rank==k1) then !OOO> k1 est le plus petit rank concernE par la frontiere
       open(unit=3,file='tmp/messages',position='append')
        write(3,*)'..............................'
        write(3,*)'subroutine wave_ww3_spec_sete'
        write(3,*)'periodpeak      ',periodpeak
        write(3,*)'peak freq x 2pi ',freq2pipeak
        write(3,*)'Pour h=         ',hmax
        write(3,*)'kvectorpeak et L',kvectorpeak,2.*pi/kvectorpeak
        write(3,*)'cwavepeak       ',cwavepeak
        write(3,*)'acm_speed       ',acm_speed
        write(3,*)'recommended "r" parameter   ' &
                   ,(x0**3)*(tanh(    kvectorpeak*hmax ) -    kvectorpeak*hmax) &
                           /(tanh( x0*kvectorpeak*hmax ) - x0*kvectorpeak*hmax)
        write(3,*)'direction spectra from:',real(ww3dir(dirw_beg)*rad2deg) &
                                     ,' to',real(ww3dir(dirw_end)*rad2deg)
        write(3,*)'period spectra from:   ',real(2.*pi/ww3puls(freq_end)) &
                                  ,' to',real(2.*pi/ww3puls(freq_beg))
       close(3)
      endif                !OOO>

! Bidouille pour le cas du modele 2DV periodique:
!     if(iglb==4) then
!        ww3phase_i=0 
!        do freq=1,freqmax ; do i=0,imax+1
!          ww3kvector(i,freq)=ww3kvector(imax/2,freq)
!        enddo             ; enddo
!     endif

!#ifdef bidon
! Definir kvectorpeak_j valeur locale de kvectorpeak
! conditions aux limites
      if(flag_allocate_==1) then  !m°v°m> !24-05-19
       allocate(kvectorpeak_j(jmax,1:2)) ; kvectorpeak_j=0.
       x2=periodpeak ! periode = 2pi/pulsation
       do k=1,2
        if(k==1)i=1
        if(k==2)i=imax
        do j=1,jmax
          x1=h_w(i,j)
          call wave_freq2kvector(x1,x2,x3) ! x3 kvector
          kvectorpeak_j(j,k)=x3
!         if(k==1)write(10+par%rank,*)j,2.*pi/kvectorpeak_j(j,k),2.*pi/kvectorpeak
        enddo
       enddo

       allocate(kvectorpeak_i(imax,1:2))
       x2=periodpeak ! periode = 2pi/pulsation
       do k=1,2
        if(k==1)j=1
        if(k==2)j=jmax
        do i=1,imax
          x1=h_w(i,j)
          call wave_freq2kvector(x1,x2,x3) ! x3 kvector
          kvectorpeak_i(i,k)=x3
!         if(k==1)write(10+par%rank,*)j,2.*pi/kvectorpeak_j(j,k),2.*pi/kvectorpeak
        enddo
       enddo


      endif                       !m°v°m>
!#endif

      end subroutine wave_ww3_spec_sete_read

!..............................................................
      subroutine wave_ww3_spec_sete_allocate
      implicit none
       allocate(sshwave_j_w(0:jmax+1     ,4,-1:4,2)) ; sshwave_j_w=0.
       allocate(sshwave_i_w(0:imax+1     ,4,-1:4,2)) ; sshwave_i_w=0.
!      allocate(  sshlf_i_w(1:imax              ,2)) ; sshlf_i_w=0.
!      allocate(  sshlf_j_w(1:jmax              ,2)) ; sshlf_j_w=0.
       allocate(qwave_j_w   (0:jmax+1,kmax,4,-1:4,2)) ; qwave_j_w=0.
       allocate(qwave_i_w   (0:imax+1,kmax,4,-1:4,2)) ; qwave_i_w=0.
       allocate(velwave_i_v (0:imax+1,kmax,3,-1:4,2)) ; velwave_i_v=0.
       allocate(velwave_j_u (0:jmax+1,kmax,3,-1:4,2)) ; velwave_j_u=0.
       allocate(sshrelax_j_w(1:jmax,2)) ; sshrelax_j_w=1.
       allocate(sshrelax_i_w(1:imax,2)) ; sshrelax_i_w=1.
      end subroutine wave_ww3_spec_sete_allocate
!..............................................................

      subroutine wave_ww3_spec_sete_obc(istr_,iend_,jstr_,jend_,id_obc_)
      use module_parallele
      implicit none
!     integer :: v0or1_=-1,istr_,iend_,jstr_,jend_,id_obc_
      integer    v0or1_   ,istr_,iend_,jstr_,jend_,id_obc_

! id_obc_=1 obc en i=1
! id_obc_=2 obc en i=imax
! id_obc_=3 obc en j=1
! id_obc_=4 obc en j=jmax

!       if(id_obc_==1) then
!        istr_=1  
!        iend_=4 
!       endif
!       if(id_obc_==2) then
!        istr_=imax-3 ;  iend_=imax
!       endif
!       if(id_obc_==3) then
!        jstr_=1 ;  jend_=4 
!       endif
!       if(id_obc_==4) then
!        jstr_=jmax-3 ;  jend_=jmax
!       endif

! POUR FAIRE UN TEST VERIFIANT que les points de grille necessaires A obc_int ont bien ete defini
! en ce qui concerne le forcage (Commenter sinon)
!     if(iteration3d==0) then ! test A commenter
!       sshobc_w=1e10
!       velobc_v=1e10
!       velobc_u=1e10
!     endif                   ! test A commenter

      if(flag_allocate_/=1)return



      if(modulo(iteration3d,wavemodulo)==0) then !modulomodulo>


      v0or1_=-1
      do dirw=dirw_beg,dirw_end ! 1,dirmax  ! 13,13 ! 4,15 (directions entrantes)

! Annulation des amplitudes si ondes sortantes:
       i0=30 ! si distance A frontiere < i0 alors ssh=0
       j0=500
       if(par%coords(1)==0)j0=1200

       x1_r4=50. ! distance de transition

#ifdef bidon
       do j=jstr_,jend_ ! jmax,jmax !0,jmax+1
        x3= min(max(real(j+par%tjmax(1)-real(j0))/x1_r4,0.),1.) ! La cote est en j=1 on annule le forcage en s'approchant de j=1
        do i=istr_,iend_ ! 0,imax+1
#ifdef bidon
         anyvar2d(i,j)=x3                                              &
           *( min(max(real(iglb-i-par%timax(1)-i0)/x1_r4,0.),1.)       &
             *max(sign(1.,cos(ww3dir(dirw)+grid_angle_t(i,j))),0.)     &
          +1.-max(sign(1.,cos(ww3dir(dirw)+grid_angle_t(i,j))),0.) )   &

           *( min(max(real(i+par%timax(1)-1-i0)/x1_r4,0.),1.)          &
             *max(sign(1.,-cos(ww3dir(dirw)+grid_angle_t(i,j))),0.)    &
          +1.-max(sign(1.,-cos(ww3dir(dirw)+grid_angle_t(i,j))),0.) )  &

           *( min(max(real(jglb-j-par%tjmax(1)-i0)/x1_r4,0.),1.)       &
             *max(sign(1.,sin(ww3dir(dirw)+grid_angle_t(i,j))),0.)     &
          +1.-max(sign(1.,sin(ww3dir(dirw)+grid_angle_t(i,j))),0.) )   
#endif
! pour experiences simples:
!!       anyvar2d(i,j)= &
!!          ( min(max(real(jglb-j-par%tjmax(1)-i0)/x1_r4,0.),1.)       &
!!           *max(sign(1.,sin(ww3dir(dirw)+grid_angle_t(i,j))),0.)     &
!!        +1.-max(sign(1.,sin(ww3dir(dirw)+grid_angle_t(i,j))),0.) )   
!        anyvar2d(i,j)= &
!         (1.-max(sign(1., sin(ww3dir(dirw)+grid_angle_t(i,j))),0.) ) &
!        *(1.-max(sign(1.,-cos(ww3dir(dirw)+grid_angle_t(i,j))),0.))  &
!        *(1.-max(sign(1.,sin(ww3dir(dirw)+grid_angle_t(i,j))),0.) )   

!         if(j==jmax.and.i==imax/2)write(6,*)'coucou jmax',anyvar2d(i,j)
!         if(j==jmax/2.and.i==imax)write(6,*)'coucou imax',anyvar2d(i,j)
!         if(j==jmax/2.and.i==1   )write(6,*)'coucou i 1 ',anyvar2d(i,j)

!         if(i==imax/2.and.par%rank==3) then
!              write(6,*)'A',min(max(real(jglb-j-par%tjmax(1)-i0)/x1_r4,0.),1.)
!              write(6,*)'sinus',sin(ww3dir(dirw)+grid_angle_t(i,j))
!              write(6,*)'anyvar2d',anyvar2d(i,j)
!         endif

!     if(par%rank==3.and.j==jmax)write(10+par%rank,*)i,anyvar2d(i,j),ww3dir(dirw),grid_angle_t(i,j)*rad2deg

         enddo ! i
       enddo   ! j
#endif
       anyvar2d=0.
       if(id_obc_==4) then !444>
        do j=jstr_,jend_ ; do i=istr_,iend_ 
!        anyvar2d(i,j)=1.-max(sign(1., sin(ww3dir(dirw)+grid_angle_t(i,j))),0.) 
         anyvar2d(i,j)=1.-max(sign(1., sin(ww3dir(dirw)+grid_angle_t(i,j))+0.1 ),0.) 
! Note: +0.1 pour disqualifier les angles trop proches de la tangente

! Bidouille Patrick
!        anyvar2d(i,j)=anyvar2d(i,j)*0.5*(1+tanh(0.15*real(i+par%timax(1)-25)))

        enddo            ; enddo
       endif               !444>
       if(id_obc_==2) then !222>
        do j=jstr_,jend_ ; do i=istr_,iend_ 
!        anyvar2d(i,j)=1.-max(sign(1.,cos(ww3dir(dirw)+grid_angle_t(i,j))),0.) 
         anyvar2d(i,j)=1.-max(sign(1.,cos(ww3dir(dirw)+grid_angle_t(i,j))+0.1 ),0.) 
        enddo            ; enddo
       endif               !222>
       if(id_obc_==1) then !111>
        do j=jstr_,jend_ ; do i=istr_,iend_ 
!        anyvar2d(i,j)=1.-max(sign(1.,-cos(ww3dir(dirw)+grid_angle_t(i,j))),0.) 
         anyvar2d(i,j)=1.-max(sign(1.,-cos(ww3dir(dirw)+grid_angle_t(i,j))+0.1),0.) 
        enddo            ; enddo
       endif               !111>

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!     stop 'potiron'

!     if(iperiodicboundary) then 
! Supplement pour supprimer le forcage sur et A l'approche des CL en i=1 et i=iglb
!       do j=jmax,jmax ; do i=0,imax+1
!        anyvar2d(i,j)=anyvar2d(i,j)*0.25*        &
!            (1+tanh((i+par%timax(1)-100. )/20.)) &
!           *(1+tanh((3900.-i-par%timax(1))/20.))  
!       enddo          ; enddo
!       anyvar2d=1.
!      endif


      do freq=freq_beg,freq_end ! 1,freqmax ! 6,6   ! 1,13
      v0or1_=min(v0or1_+1,1)    ! v0or1_=0 puis v0or1_=1

      x0=sqrt(2.*ww3efth(dirw,freq,1,1)*ww3deltafreq(freq)*ww3deltadir) !27-03-19

      do j=jstr_,jend_   ; do i=istr_,iend_

! Decalage d'indices dans le temps et l'espace expliquE dans
! https://docs.google.com/document/d/1d7zaixIIH62c3ss0sD7ui_L3BV9MmesNJl9zfh28yUc/edit
! Elevation de la surface au meme temps que la ssh au t+1
       xy_t(i,j,1)=anyvar2d(i,j)*                               &
                   x0*sin( ww3puls(freq)*(elapsedtime_now+dti_fw) &
                          -ww3phase(i,j,freq,dirw) )
!                          -2.*pi*real(0*modulo(i+par%timax(1),iglb-2))/real(iglb-2))
! Elevation de la surface au meme temps que la ssh au t (pour les vitesses)
       xy_t(i,j,0)=anyvar2d(i,j)*                               &
                   x0*sin( ww3puls(freq)* elapsedtime_now       &
                          -ww3phase(i,j,freq,dirw) )

! note: parce que dans le sens perpendiculaire A la frontiere le vecteur d'onde est pris constant 
! ainsi que grid_angle et dy_t (pris en jmax constant expres) alors le calcul de
! la phase n'a pas besoin du cumul appliquE precedement dans la direction tangente 
      enddo  ! j
      enddo  ! i

!#ifdef bidon
! Profil du courant au point t
      do j=jstr_,jend_   ; do i=istr_,iend_
      k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
      kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

! mukvector THEORIQUE
!     mukvector=kvector ! mu*k avec mu**2=Eq21 MarsaleixetalOM2019
! mukvector ACM
      mukvector=sqrt(1-(ww3puls(freq)/(kvector*acm_speed))**2)*kvector ! mu*k avec mu**2=Eq21 MarsaleixetalOM2019

!     if(i+par%timax(1)==iglb.and.j+par%tjmax(1)==jglb) &
!     write(666,*)'mu',sqrt(1-(ww3puls(freq)/(kvector*acm_speed))**2)

      if(mukvector*h_w(i,j)<7.) then !pmx>

! cas kvector*H petit:
       do k=1,kmax
! Diviser par cosh ou sinh? voir discusion dans:
! https://docs.google.com/document/d/1d7zaixIIH62c3ss0sD7ui_L3BV9MmesNJl9zfh28yUc/edit
! Methode 1 vel=omega*cosh(k(z+h))/sinh(kh)*ssh
!       anyvar3d(i,j,k)=cosh(mukvector*(depth_t(i,j,k)+h_w(i,j))) &
!                      /sinh(mukvector*max(1.         ,h_w(i,j)))
! Methode 2 vel=(grav*k/omega)*cosh(k(z+h))/cosh(kh)*ssh=(grav*k/omega)*f(z)*ssh
        anyvar3d(i,j,k)=cosh(mukvector*(depth_t(i,j,k)+h_w(i,j))) &
                       /cosh(mukvector*max(1.         ,h_w(i,j)))
! Profil vertical apparaissant dans Eq22 et Eq23
       enddo

      else                                     !pmx>

! cas kvector*H grand:
       do k=1,kmax
! Diviser par cosh ou sinh? voir discusion dans:
! https://docs.google.com/document/d/1d7zaixIIH62c3ss0sD7ui_L3BV9MmesNJl9zfh28yUc/edit
! Methode 1 vel=omega*exp(kz)*ssh
!       anyvar3d(i,j,k)=exp(mukvector*depth_t(i,j,k))
! Methode 2 vel=omega*exp(kz)*ssh=(grav*k/omega)*f(z)*ssh avec f(z)=exp(kz)*(omega**2)/(grav*k)
        anyvar3d(i,j,k)=exp(mukvector*depth_t(i,j,k))*(ww3puls(freq)**2)/(grav*KVECTOR)
       enddo

      endif                                    !pmx>

      enddo         ; enddo ! i,j

!#ifdef bidon
      if(id_obc_==1.or.id_obc_==2) then !>>>>

! courant Oi 3D au point u au temps 2 du courant
       do j=jstr_,jend_   
       do i=istr_,iend_
! Composante Oi du courant de surface points t (voir eq26 papier nh)
! Methode 1 u=omega*f(z)*ssh
!       xy_t(i,j,2)=cos( ww3dir(dirw)+grid_angle_t(i,j) )*ww3puls(freq)*xy_t(i,j,1)
! Methode 2 u=cos(dir)*(grav*k/omega)*f(z)*ssh
        k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 
        xy_t(i,j,2)=cos( ww3dir(dirw)+grid_angle_t(i,j) ) &
                   *(grav*kvector/ww3puls(freq))*xy_t(i,j,0) ! g*k/w*ssh debut Eq23 Marsaleix et al OM2019
       enddo           
       enddo
       do i1=1,3
! Decalage d'indices dans le temps et l'espace expliquE dans
! https://docs.google.com/document/d/1d7zaixIIH62c3ss0sD7ui_L3BV9MmesNJl9zfh28yUc/edit
        if(id_obc_==1)i=i1+1
        if(id_obc_==2)i=imax-(i1-1) !02-12-19
        do j=jstr_,jend_   
        k0=max( int(h_u(i,j)),1) ; rap=max((h_u(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

! x1 THEORIQUE:
!       x1=1. ! x1=r/mu**2 dans Eq21
! x1 ACM:
        x1=nhpgf_reduce/(1.-(ww3puls(freq)/(kvector*acm_speed))**2) ! x1=r/mu**2 dans Eq21

        do k=1,kmax 
          velwave_j_u(j,k,i1,timewavenext,id_obc_)=          &
          velwave_j_u(j,k,i1,timewavenext,id_obc_)*v0or1_    &
           +0.5*(xy_t(i-1,j,2)+xy_t(i,j,2))*(1.+x1*(0.5*(anyvar3d(i-1,j,k)+anyvar3d(i,j,k))-1.)) !02-12-19
! Eq23: u=g*k/w*ssh*(1+r/mu**2*(cosh(mu*k*(z+h))/cosh(mu*k*h)-1))
         enddo

        enddo
       enddo
       do i1=1,4

        if(id_obc_==1)i=i1
        if(id_obc_==2)i=imax-(i1-1)

        do j=jstr_,jend_   
         sshwave_j_w(j,i1,timewavenext,id_obc_)= &
         sshwave_j_w(j,i1,timewavenext,id_obc_)*v0or1_+xy_t(i,j,1)
        enddo
        do j=jstr_,jend_   
         k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
         kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq)
! x1 THEORIQUE:
!       x1=grav*xy_t(i,j,1) 
! x1 ACM:
         x1=grav*xy_t(i,j,1)/(1-(ww3puls(freq)/(kvector*acm_speed))**2) ! g*ssh/mu**2 (debut Eq22)
         do k=1,kmax
! Methode 1: le profil f(z)=cosh(k(z+h))/sinh(kh)
!        qwave_j_w(j,k,i1,timewavenext,id_obc_)=         &
!        qwave_j_w(j,k,i1,timewavenext,id_obc_)*v0or1_   &
!        +grav*xy_t(i,j,1)*(anyvar3d(i,j,k)*(ww3puls(freq)**2)/(grav*KVECTOR)-1.)
! Methode 2: le profil f(z)=cosh(k(z+h))/cosh(kh)
          qwave_j_w(j,k,i1,timewavenext,id_obc_)=         &
          qwave_j_w(j,k,i1,timewavenext,id_obc_)*v0or1_   &
          +x1*(anyvar3d(i,j,k)-1.) ! Eq22
!         +grav*xy_t(i,j,1)*(anyvar3d(i,j,k)-1.)
         enddo
        enddo
       enddo ! i1

      endif                             !>>>>
      if(id_obc_==3.or.id_obc_==4) then !>>>>

! courant Oj 3D
       do j=jstr_,jend_   
       do i=istr_,iend_
! Composante Oj du courant 
! Methode 1 v=omega*f(z)*ssh
!       xy_t(i,j,2)=sin( ww3dir(dirw)+grid_angle_t(i,j) )*ww3puls(freq)*xy_t(i,j,1)
! Methode 2 v=sin(dir)*(grav*k/omega)*f(z)*ssh
        k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 
        xy_t(i,j,2)=sin( ww3dir(dirw)+grid_angle_t(i,j) ) &
                               *(grav*kvector/ww3puls(freq))*xy_t(i,j,0)
       enddo           
       enddo
! Nouvelles lignes
       do j1=1,3
! Decalage d'indices dans le temps et l'espace expliquE dans
! https://docs.google.com/document/d/1d7zaixIIH62c3ss0sD7ui_L3BV9MmesNJl9zfh28yUc/edit
        if(id_obc_==3)j=j1+1
        if(id_obc_==4)j=jmax-(j1-1) !02-12-19
        do i=istr_,iend_   
        k0=max( int(h_v(i,j)),1) ; rap=max((h_v(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

! x1 THEORIQUE:
!       x1=1.
! x1 ACM:
        x1=nhpgf_reduce/(1.-(ww3puls(freq)/(kvector*acm_speed))**2) ! x1=r/mu**2 dans Eq21

        do k=1,kmax 
          velwave_i_v(i,k,j1,timewavenext,id_obc_-2)=          &
          velwave_i_v(i,k,j1,timewavenext,id_obc_-2)*v0or1_    &
           +0.5*(xy_t(i,j-1,2)+xy_t(i,j,2))*(1.+x1*(0.5*(anyvar3d(i,j-1,k)+anyvar3d(i,j,k))-1.))  !02-12-19
! Eq23: v=g*k/w*ssh*(1+r/mu**2*(cosh(mu*k*(z+h))/cosh(mu*k*h)-1))
         enddo
        enddo
       enddo
       do j1=1,4
        if(id_obc_==3)j=j1
        if(id_obc_==4)j=jmax-(j1-1)
        do i=istr_,iend_   
         sshwave_i_w(i,j1,timewavenext,id_obc_-2)= &
         sshwave_i_w(i,j1,timewavenext,id_obc_-2)*v0or1_+xy_t(i,j,1)
        enddo
        do i=istr_,iend_   
         k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
         kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq)

! x1 THEORIQUE
!        x1=grav*xy_t(i,j,1)
! x1 ACM
         x1=grav*xy_t(i,j,1)/(1-(ww3puls(freq)/(kvector*acm_speed))**2) ! g*ssh/mu**2 (debut Eq22)

         do k=1,kmax
! Methode 1: le profil f(z)=cosh(k(z+h))/sinh(kh)
!        qwave_i_w(i,k,j1,timewavenext,id_obc_-2)=         &
!        qwave_i_w(i,k,j1,timewavenext,id_obc_-2)*v0or1_   &
!        +grav*xy_t(i,j,1)*(anyvar3d(i,j,k)*(ww3puls(freq)**2)/(grav*KVECTOR)-1.)
! Methode 2: le profil f(z)=cosh(k(z+h))/cosh(kh)
          qwave_i_w(i,k,j1,timewavenext,id_obc_-2)=         &
          qwave_i_w(i,k,j1,timewavenext,id_obc_-2)*v0or1_   &
          +x1*(anyvar3d(i,j,k)-1.) ! Eq22
!         +grav*xy_t(i,j,1)*(anyvar3d(i,j,k)-1.)
         enddo
        enddo
       enddo ! j1

      endif                             !>>>>


!#endif


#ifdef bidon
! Ici on verifie que la divergence du transport deduite de velobc 
! donne bien la variation de ssh escomptee
      x10=0.1 ! Pour verif, joue le role du pas de temps
      fluxbar_u=0. ; fluxbar_v=0. 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
       fluxbar_u(i,j,1)=fluxbar_u(i,j,1)+velobc_u(i,j,k,1)*dy_u(i,j)*dz_u(i,j,k,1)
      enddo ; enddo         ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
       fluxbar_v(i,j,1)=fluxbar_v(i,j,1)+velobc_v(i,j,k,1)*dx_v(i,j)*dz_v(i,j,k,1)
      enddo ; enddo         ; enddo
      do j=jmax/2,jmax/2   ; do i=2,imax-1
       write(10+par%rank,*) &
         (ssh_int_w(i,j,2)-ssh_int_w(i,j,0))/x10              &
                    ,-(                                       &
                        ( fluxbar_u(i+1,j  ,1)                &
                         -fluxbar_u(i  ,j  ,1)                &
                         +fluxbar_v(i  ,j+1,1)                &
                         -fluxbar_v(i  ,j  ,1))/dxdy_t(i,j)   &
                      )
      enddo ; enddo
#endif

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
      enddo ! boucle frequence
      enddo ! boucle direction

! 16-10-19 ! Appliquer une condition intertidale sur la SSH reconstituee:
      if(id_obc_==1.or.id_obc_==2) then !>>>>
       do i1=1,4
        if(id_obc_==1)i=i1
        if(id_obc_==2)i=imax-(i1-1)
        do j=jstr_,jend_
                sshwave_j_w(j,i1,timewavenext,id_obc_)= &
            max(sshwave_j_w(j,i1,timewavenext,id_obc_),-h_w(i,j))
        enddo
       enddo ! i1
      endif                             !>>>>
      if(id_obc_==3.or.id_obc_==4) then !>>>>
       do j1=1,4
        if(id_obc_==3)j=j1
        if(id_obc_==4)j=jmax-(j1-1)
        do i=istr_,iend_
             sshwave_i_w(i,j1,timewavenext,id_obc_-2)= &
         max(sshwave_i_w(i,j1,timewavenext,id_obc_-2),-h_w(i,j))
        enddo
       enddo ! j1
      endif                             !>>>>

! Suite: pour faire une reconstitution du signal
!      sum1=sum1+1.
!      sum2=sum2+sshobc_w(imax/2,jmax,1)**2
    
!      write(64,*)elapsedtime_now,sshobc_w(imax/2,jmax,1)
!     enddo ! loop1

!     write(6,*)'sum1',sum1
!     write(6,*)'sum2',sum2
!     write(6,*)'Variance SSH',sum2/sum1
!     write(6,*)'Ecart type SSH',sqrt(sum2/sum1)
!     write(6,*)'HS (4 fois ecart type)',4.*sqrt(sum2/sum1)

!     write(6,*)'Grid angle',grid_angle_t(imax/2,jmax/2)*rad2deg
!     i=imax/2
!     do j=1,jmax
!      write(65,*)lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
!     enddo

!     dirw=10
!     write(6,*)'cos sin',cos( ww3dir(dirw)+grid_angle_t(imax/2,jmax/2) ) &
!                        ,sin( ww3dir(dirw)+grid_angle_t(imax/2,jmax/2) )
!     write(6,*)'dir+gangle',(ww3dir(dirw)+grid_angle_t(imax/2,jmax/2))*rad2deg

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'STOP wave_ww3_spec_sete'

      endif                                      !modulomodulo>

! Interpolation temporelle A partir de XXXobc(0) et XXXobc(2):

!     if(par%rank==3)write(6,*)'Rampe',rampe,elapsedtime_now
!     rampe=1
!     if(par%rank==3)write(6,*)'BIDOUILLE Rampe',rampe,elapsedtime_now

! Les echeances de mise à jour sont timewaveprev et timewavenext (valeurs 0 ou 3 alternees)
! Le temps 2 est l'interpolation temporelle au temps 2 du modele (par ex vel_u(:,:,:,2))
! Le temps 1 correspond au temps 1 du modele (par ex vel_u(:,:,:,1))
!     x2_r4=real(modulo(iteration3d,wavemodulo))/real(wavemodulo)
!     x0_r4=1.-x2_r4
!     x2_r4=x2_r4*rampe
!     x0_r4=x0_r4*rampe

! Coefs d'une interpolation sur 4 temps:
! https://fr.wikipedia.org/wiki/Interpolation_lagrangienne
! voir aussi programme ../../tools/interpolation_polynomes_lagrange.F90
! Principe: le temps time_r4 est un temps adimensionnE allant de 0 A 1 entre
! les echeances timewaveprev2 et timewaveprev. Les echeances timewaveprev3 et timewavenext
! encadrant ces dernieres. AU final timewaveprev3<timewaveprev2<time_r4<timewaveprev<timewavenext
! Comme en vrai timewaveprev3 et cie ne sont que des entiers pour l'indexation des tableaux
! fortran; les temps associes pris en compte par l'algo de Lagragrange sont respectivement -1,0,1,2
      time_r4=real(modulo(iteration3d,wavemodulo))/real(wavemodulo)
      x0_r4=- time_r4   *(time_r4-1)*(time_r4-2)*0.166666666666667 &
            *rampe
      x1_r4= (time_r4+1)*(time_r4-1)*(time_r4-2)*0.5               &
            *rampe
      x2_r4=-(time_r4+1)* time_r4   *(time_r4-2)*0.5               &
            *rampe
      x3_r4= (time_r4+1)* time_r4   *(time_r4-1)*0.166666666666667 &
            *rampe

! Revenir A un algo sur 2 points lineaires entre timewaveprev2 et timewaveprev
!     x2_r4=real(modulo(iteration3d,wavemodulo))/real(wavemodulo)
!     x1_r4=1.-x2_r4
!     x3_r4=0.
!     x0_r4=0.
!     x2_r4=x2_r4*rampe
!     x1_r4=x1_r4*rampe

!      if(id_obc_==4)write(666,*)'time_r4',time_r4

      if(id_obc_==1.or.id_obc_==2) then !>>>>
       do i1=1,4
        do j=jstr_,jend_
         sshwave_j_w(j,i1,1,id_obc_)=sshwave_j_w(j,i1,2,id_obc_)
!        sshwave_j_w(j,i1,2,id_obc_)=x0_r4*sshwave_j_w(j,i1,timewaveprev,id_obc_) &
!                                   +x2_r4*sshwave_j_w(j,i1,timewavenext,id_obc_)  
         sshwave_j_w(j,i1,2,id_obc_)=x0_r4*sshwave_j_w(j,i1,timewaveprev3,id_obc_) &
                                    +x1_r4*sshwave_j_w(j,i1,timewaveprev2,id_obc_) & 
                                    +x2_r4*sshwave_j_w(j,i1,timewaveprev ,id_obc_) & 
                                    +x3_r4*sshwave_j_w(j,i1,timewavenext ,id_obc_)  
! bidouille pour graphique:
        if(id_obc_==1)i=i1
        if(id_obc_==2)i=imax-(i1-1)
         sshobc_w(i,j,:)=sshwave_j_w(j,i1,2,id_obc_)

        enddo
        do k=1,kmax
        do j=jstr_,jend_
         qwave_j_w(j,k,i1,1,id_obc_)=qwave_j_w(j,k,i1,2,id_obc_)
!        qwave_j_w(j,k,i1,2,id_obc_)=x0_r4*qwave_j_w(j,k,i1,timewaveprev,id_obc_) &
!                                   +x2_r4*qwave_j_w(j,k,i1,timewavenext,id_obc_)  
         qwave_j_w(j,k,i1,2,id_obc_)=x0_r4*qwave_j_w(j,k,i1,timewaveprev3,id_obc_) &
                                    +x1_r4*qwave_j_w(j,k,i1,timewaveprev2,id_obc_) & 
                                    +x2_r4*qwave_j_w(j,k,i1,timewaveprev ,id_obc_) & 
                                    +x3_r4*qwave_j_w(j,k,i1,timewavenext ,id_obc_)  
        enddo
        enddo
       enddo ! i1
       do i1=1,3
        do k=1,kmax
         do j=jstr_,jend_
          velwave_j_u(j,k,i1,1,id_obc_)=velwave_j_u(j,k,i1,2,id_obc_)
!         velwave_j_u(j,k,i1,2,id_obc_)=x0_r4*velwave_j_u(j,k,i1,timewaveprev,id_obc_) &
!                                      +x2_r4*velwave_j_u(j,k,i1,timewavenext,id_obc_)  
         velwave_j_u(j,k,i1,2,id_obc_)=x0_r4*velwave_j_u(j,k,i1,timewaveprev3,id_obc_) &
                                      +x1_r4*velwave_j_u(j,k,i1,timewaveprev2,id_obc_) & 
                                      +x2_r4*velwave_j_u(j,k,i1,timewaveprev ,id_obc_) & 
                                      +x3_r4*velwave_j_u(j,k,i1,timewavenext ,id_obc_)  
         enddo
        enddo
       enddo
      endif                             !>>>>
      if(id_obc_==3.or.id_obc_==4) then !>>>>
       do j1=1,4
        do i=istr_,iend_
         sshwave_i_w(i,j1,1,id_obc_-2)=sshwave_i_w(i,j1,2,id_obc_-2)
!        sshwave_i_w(i,j1,2,id_obc_-2)=x0_r4*sshwave_i_w(i,j1,timewaveprev,id_obc_-2) &
!                                     +x2_r4*sshwave_i_w(i,j1,timewavenext,id_obc_-2)  
         sshwave_i_w(i,j1,2,id_obc_-2)=x0_r4*sshwave_i_w(i,j1,timewaveprev3,id_obc_-2) &
                                      +x1_r4*sshwave_i_w(i,j1,timewaveprev2,id_obc_-2) & 
                                      +x2_r4*sshwave_i_w(i,j1,timewaveprev ,id_obc_-2) & 
                                      +x3_r4*sshwave_i_w(i,j1,timewavenext ,id_obc_-2)  
!        if(j1==1)write(10+id_obc_,*)i,sshwave_i_w(i,j1,0,id_obc_-2),sshwave_i_w(i,j1,3,id_obc_-2)
!        if(j1==4)write(20+id_obc_,*)i,sshwave_i_w(i,j1,0,id_obc_-2),sshwave_i_w(i,j1,3,id_obc_-2)

! bidouille pour graphique:
        if(id_obc_==3)j=j1
        if(id_obc_==4)j=jmax-(j1-1)
         sshobc_w(i,j,:)=sshwave_i_w(i,j1,2,id_obc_-2)
        enddo
        do k=1,kmax
        do i=istr_,iend_
         qwave_i_w(i,k,j1,1,id_obc_-2)=qwave_i_w(i,k,j1,2,id_obc_-2)
!        qwave_i_w(i,k,j1,2,id_obc_-2)=x0_r4*qwave_i_w(i,k,j1,timewaveprev,id_obc_-2) &
!                                     +x2_r4*qwave_i_w(i,k,j1,timewavenext,id_obc_-2)  
         qwave_i_w(i,k,j1,2,id_obc_-2)=x0_r4*qwave_i_w(i,k,j1,timewaveprev3,id_obc_-2) &
                                      +x1_r4*qwave_i_w(i,k,j1,timewaveprev2,id_obc_-2) & 
                                      +x2_r4*qwave_i_w(i,k,j1,timewaveprev ,id_obc_-2) & 
                                      +x3_r4*qwave_i_w(i,k,j1,timewavenext ,id_obc_-2)  
! bidouille pour graphique:
!        if(id_obc_==3)j=j1
!        if(id_obc_==4)j=jmax-(j1-1)
!        q_t(i,j,k,:)=qwave_i_w(i,k,j1,2,id_obc_-2)
        enddo
        enddo
       enddo ! j1
       do j1=1,3
        if(id_obc_==3)j=j1+1
        if(id_obc_==4)j=jmax-(j1-1)
        do k=1,kmax
         do i=istr_,iend_
          velwave_i_v(i,k,j1,1,id_obc_-2)=velwave_i_v(i,k,j1,2,id_obc_-2)
!         velwave_i_v(i,k,j1,2,id_obc_-2)=x0_r4*velwave_i_v(i,k,j1,timewaveprev,id_obc_-2) &
!                                        +x2_r4*velwave_i_v(i,k,j1,timewavenext,id_obc_-2)  
          velwave_i_v(i,k,j1,2,id_obc_-2)=x0_r4*velwave_i_v(i,k,j1,timewaveprev3,id_obc_-2) &
                                         +x1_r4*velwave_i_v(i,k,j1,timewaveprev2,id_obc_-2) & 
                                         +x2_r4*velwave_i_v(i,k,j1,timewaveprev ,id_obc_-2) & 
                                         +x3_r4*velwave_i_v(i,k,j1,timewavenext ,id_obc_-2)  
! bidouille pour graphique:
!        if(id_obc_==3)j=j1+1
!        if(id_obc_==4)j=jmax-(j1-1)
!        vel_v(i,j,k,:)=velwave_i_v(i,k,j1,2,id_obc_-2)
         enddo
        enddo
       enddo
      endif                             !>>>>

! Anciennes lignes:
! sshobc_w
!     do j=jstr_,jend_   ; do i=istr_,iend_
!      sshobc_w(i,j,1)=x0_r4*sshobc_w(i,j,timewaveprev)  &
!                     +x2_r4*sshobc_w(i,j,timewavenext)
!     enddo            ; enddo

!#ifdef bidon
! velobc_u
!      do k=1,kmax ; do j=jstr_,jend_   ; do i=istr_+1,iend_ !     do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1
!       velobc_u(i,j,k,1)=x0_r4*velobc_u(i,j,k,timewaveprev) &
!                        +x2_r4*velobc_u(i,j,k,timewavenext)
!      enddo       ; enddo         ; enddo
!#endif
! velobc_v
!      do k=1,kmax ; do j=jstr_+1,jend_ ; do i=istr_,iend_ !     do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1
!       velobc_v(i,j,k,1)=x0_r4*velobc_v(i,j,k,timewaveprev) &
!                        +x2_r4*velobc_v(i,j,k,timewavenext)
!      enddo       ; enddo        ; enddo

!     write(6,*)'cou',modulo(iteration3d,wavemodulo),iteration3d,wavemodulo,x2_r4
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      if(modulo(iteration3d,wavemodulo)==0)stop 'verifier x2 facteur sshobc non nul au depart'

#ifdef bidon
! Sur les bords Est et Ouest du domaine de Sete, on enleve du forcage la
! composante grande longueur d'onde du modele (obtenue par moyenne sur i ! de vel_u vel_v
! pour laisser se developper des courants sortants basse frequence
      if(obcstatus(ieqimax)==1.or.obcstatus(ieq1)==1) then !pmx>
! moyenne vel_u
       anyvar3d(1,1:jmax,1:kmax)=0.
        do k=1,kmax ; do j=1,jmax ; do i=2,imax  
         anyvar3d(1,j,k)=anyvar3d(1,j,k)+vel_u(i,j,k,0)
        enddo       ; enddo         ; enddo
        x0_r4=1./real(imax-1)
        do j=1,jmax
         x1_r4=x0_r4*min(max( real(1000-j-par%tjmax(1))/100.,0.),1.) !x1_r4 pour appliquer pres de la cote et passer progressivement A zero au large
         do k=1,kmax
          anyvar3d(1,j,k)=anyvar3d(1,j,k)*x1_r4
         enddo       
        enddo
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         velobc_u(i,j,k,1)=velobc_u(i,j,k,1)+anyvar3d(1,j,k)
        enddo       ; enddo       ; enddo

! moyenne vel_v
       anyvar3d(2,1:jmax+1,1:kmax)=0.
        do k=1,kmax ; do j=1,jmax+1 ; do i=2,imax-1  
         anyvar3d(2,j,k)=anyvar3d(2,j,k)+vel_v(i,j,k,0)
        enddo       ; enddo       ; enddo
        x0_r4=1./real(imax-2)
        do j=1,jmax+1 
         x1_r4=x0_r4*min(max( real(1000-j-par%tjmax(1))/100.,0.),1.) !x1_r4 pour appliquer pres de la cote et passer progressivement A zero au large
         do k=1,kmax
          anyvar3d(2,j,k)=anyvar3d(2,j,k)*x1_r4
         enddo       
        enddo
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         velobc_v(i,j,k,1)=velobc_v(i,j,k,1)+anyvar3d(2,j,k)
        enddo       ; enddo       ; enddo
      endif                                                !pmx>
#endif

      end subroutine wave_ww3_spec_sete_obc

!--------------------------------------------------------------------------------

      subroutine wave_phase_jeqjmax_increasing_i
      implicit none

! La phase en (i,j)=(0,0) est donnEe par une fonction aleatoire (voir cours F.Ardhuin)
      if(allocated(ww3phase)) then !pmx>
       do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end
        call random_number(x1_r4)
        ww3phase(0,0,freq,dirw)=x1_r4*2.*pi
       enddo ; enddo
      endif                        !pmx>

      if(nint(sponge_l)-jmax>0) &
      stop 'Err 24 wave_phase_jeqjmax_increasing_i sponge_l>jmax'

      do i10=0,nbdom_imax-1 !>>>>>i10>>>>>>

! ETAPE 1: Petit bord en i=0, de j=0 A jmax
      if(par%coords(2)==nbdom_jmax-1) then ! Etre sur le bord jglb
      if(par%coords(1)==0) then !ooo>

       i=0
       do j=1,jmax+1
        k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

        x1=lon_t(i,j)-lon_t(i,j-1)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi

        do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end

         ww3phase(i,j  ,freq,dirw)=                &
         ww3phase(i,j-1,freq,dirw)                 &
        +rayonterre*KVECTOR*(             &
             cos(ww3dir(dirw))*( x1 )*cos(0.5*(lat_t(i,j)+lat_t(i,j-1))) &
            +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i,j-1)) )

!        write(10+par%rank,*)i10,i+par%timax(1),j+par%tjmax(1),ww3phase(i,j,freq,dirw)

        enddo ; enddo ! dirw,freq

       enddo  ! j       

      endif                     !ooo>
      endif                       !Etre sur le bord jglb

! ETAPE 2: Petit bord en j=0, de i=0 A imax+1
      if(par%coords(2)==nbdom_jmax-1) then ! Etre sur le bord jglb

       j=0
       do i=1,imax+1
        k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

        x1=lon_t(i,j)-lon_t(i-1,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi

        do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end

         ww3phase(i  ,j,freq,dirw)=                &
         ww3phase(i-1,j,freq,dirw)                 &
        +rayonterre*KVECTOR*(             &
             cos(ww3dir(dirw))*( x1 )*cos(0.5*(lat_t(i,j)+lat_t(i-1,j))) &
            +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i-1,j)) )

        enddo ; enddo ! dirw,freq

       enddo  ! i       

      endif                       !Etre sur le bord jglb

! ETAPE 3: faire des grands bords "moyennes des 2 routes"

      if(par%coords(2)==nbdom_jmax-1) then ! Etre sur le bord jglb

        do j=1,jmax+1
        do i=1,imax+1
         k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
         kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

! route 1: integrer sur i
        x1=lon_t(i,j)-lon_t(i-1,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
! route 2: integrer sur j
        x2=lon_t(i,j)-lon_t(i,j-1)
        if(x2<-pi)x2=x2+2.*pi
        if(x2> pi)x2=x2-2.*pi

        do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end

         ww3phase(i  ,j,freq,dirw)= &

         0.5*( & !route 1
         ww3phase(i-1,j,freq,dirw)                      &
        +rayonterre*KVECTOR*(                  &
         cos(ww3dir(dirw))*( x1 )*cos(0.5*(lat_t(i,j)+lat_t(i-1,j))) &
        +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i-1,j)) ) &
             ) & !route 1

        +0.5*( & !route 2
         ww3phase(i,j-1,freq,dirw)                     &
        +rayonterre*KVECTOR*(                 &
             cos(ww3dir(dirw))*( x2 )*cos(0.5*(lat_t(i,j)+lat_t(i,j-1))) &
            +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i,j-1))) &
             )   !route 2

!       if(j+par%tjmax(1)==jglb)write(10+par%rank,*)i+par%timax(1),ww3phase(i  ,j,freq,dirw)
!       if(j+par%tjmax(1)==jglb)write(10+par%rank,*)i+par%timax(1) &
!           ,ww3phase(i-1,j,freq,dirw),ww3phase(i,j-1,freq,dirw)
        enddo ; enddo ! dirw,freq


        enddo ! i
        enddo ! j

       endif                       ! Etre sur le bord jglb

       call wave_mpi_eastward

      enddo                 !>>>>>i10>>>>>>

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'coco4'
      end subroutine wave_phase_jeqjmax_increasing_i

!............................................................................

      subroutine wave_phase_decreasing_j
      implicit none
      integer :: jend_


      do j10=nbdom_jmax-1,0,-1 !>>>>>j10>>>>>>
      if(allocated(ww3phase)) then !pmx>

! ETAPE 1: Bord i=0 
!          de j=0 A jmax-3 
!          Rien faire si par%coord(2)=nbdom_jmax-1 car deja fait

      jend_=jmax-3
      if(par%coords(2)==nbdom_jmax-1)jend_=-999

       i=0
       do j=jend_,0,-1
        k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

        x1=lon_t(i,j)-lon_t(i,j+1)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi

        do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end

         ww3phase(i,j  ,freq,dirw)=                &
         ww3phase(i,j+1,freq,dirw)                 &
        +rayonterre*KVECTOR*(             &
             cos(ww3dir(dirw))*( x1 )*cos(0.5*(lat_t(i,j)+lat_t(i,j+1))) &
            +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i,j+1)) )

        enddo ; enddo ! dirw,freq

       enddo  ! j 
        

! ETAPE 2: "2 routes"
!          de i=1,imax+1
!          de j=0 A jmax+1 CAS GENERAL
!          de j=0 A jmax-3 pour le dom juste en dessous de par%coord(2)=nbdom_jmax-1 soit nbdom_jmax-2
!          Rien faire si par%coord(2)=nbdom_jmax-1 car deja fait

       do j=jend_,0,-1
       do i=1,imax+1
        k0=max( int(h_w(i,j)),1) ; rap=max((h_w(i,j)),1.)-real(k0) 
        kvector=(1.-rap)*ww3kvector(k0,freq)+rap*ww3kvector(k0+1,freq) 

! route 1: integrer sur i
        x1=lon_t(i,j)-lon_t(i-1,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
! route 2: integrer sur j
        x2=lon_t(i,j)-lon_t(i,j+1)
        if(x2<-pi)x2=x2+2.*pi
        if(x2> pi)x2=x2-2.*pi

        do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end

         ww3phase(i  ,j,freq,dirw)= &

         0.5*( & !route 1
         ww3phase(i-1,j,freq,dirw)                      &
        +rayonterre*KVECTOR*(                  &
         cos(ww3dir(dirw))*( x1 )*cos(0.5*(lat_t(i,j)+lat_t(i-1,j))) &
        +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i-1,j)) ) &
             ) & !route 1

        +0.5*( & !route 2
         ww3phase(i,j+1,freq,dirw)                     &
        +rayonterre*KVECTOR*(                 &
             cos(ww3dir(dirw))*( x2 )*cos(0.5*(lat_t(i,j)+lat_t(i,j+1))) &
            +sin(ww3dir(dirw))*(lat_t(i,j)-lat_t(i,j+1))) &
             )   !route 2

        enddo ; enddo ! dirw,freq

       enddo   ! j
       enddo   ! i
      endif                        !pmx> allocated ww3phase

      call wave_mpi_southward(3,jmax+1)
      call wave_mpi_southward(2,jmax  )
      call wave_mpi_southward(1,jmax-1)
      call wave_mpi_southward(0,jmax-2)

      enddo                    !>>>>>j10>>>>>>

!     if(allocated(ww3phase)) then !ooo>
!      i=imax/2
!      do j=jmax+1,0,-1
!       write(10+par%rank,*)j+par%tjmax(1),ww3phase(i,j,freq_beg,dirw_beg)
!      enddo
!     endif                        !ooo>

      end subroutine wave_phase_decreasing_j

!..............................................................................

      subroutine wave_mpi_eastward
      implicit none
      integer,parameter :: nexchgmax_mpi=2 &
                          ,tagest_mpi      =1617  
      integer :: nexchg_=0,nb_send_est,nb_recv_ouest
      integer,dimension(nexchgmax_mpi   ) :: tabreq_est_mpi
      integer,dimension(mpi_status_size,nexchgmax_mpi   ) :: status_est_mpi
!     real,dimension(:),allocatable :: ar_send_est,ar_recv_ouest

      nb_send_est=(jmax+2)*(dirw_end-dirw_beg+1)*(freq_end-freq_beg+1)
      nb_recv_ouest=nb_send_est

      if(.not.allocated(ar_send_est))  allocate(ar_send_est  (nb_send_est))
      if(.not.allocated(ar_recv_ouest))allocate(ar_recv_ouest(nb_recv_ouest))

      ar_recv_ouest=0.
      ar_send_est=0.
      if(allocated(ww3phase)) then !m°v°m>
       k=0 ; i=imax-2
       do j=0,jmax+1 ; do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end
        k=k+1
        ar_send_est(k)=ww3phase(i,j,freq,dirw)
       enddo ; enddo ; enddo
      endif                        !m°v°m>
      

!.............................................................................................
! Les domaines indiquent à leurs voisins le nombre d'element s'appretant à traverser les
! frontieres. 
       nexchg_=0

! Je suis A l'ouest , j'enverrai A destination de l'Est
      if (par%tvoisin(est) /= mpi_proc_null) then !ooo>

! Comme ce test n'est pas suffisant dans le cas academique periodique j'ajoute comme
! condition supplementaire qu'on n'envoie pas A l'est si on est le dernier dom A l'est
       if(par%coords(1)/=nbdom_imax-1) then ! suplement cas periodique >

       nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
       call mpi_issend(nb_send_est,1,mpi_integer,par%tvoisin(est)  &
              ,tagest_mpi   ,par%comm2d,tabreq_est_mpi   (nexchg_   ),ierr)

      endif                                 ! suplement cas periodique >
      endif                                       !ooo>

! Je suis A l'Est, je recois en provenance de l'Ouest
      if (par%tvoisin(ouest) /= mpi_proc_null) then !ooo>
       if(par%coords(1)/=0           ) then ! suplement cas periodique >
      nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
       call mpi_irecv(nb_recv_ouest,1,mpi_integer,par%tvoisin(ouest)  &
              ,tagest_mpi   ,par%comm2d,tabreq_est_mpi   (nexchg_   ),ierr)
      endif                                 ! suplement cas periodique >
      endif                                         !ooo>

      ! La salle d'attente
      call mpi_waitall(nexchg_                  &
                      ,tabreq_est_mpi   (1:nexchg_   ) &
                      ,status_est_mpi  (:,1:nexchg_ ) &
                      ,ierr)


! Echange
      nexchg_=0
      if (par%tvoisin(est) /= mpi_proc_null) then !ooo>
      if(par%coords(1)/=nbdom_imax-1) then ! suplement cas periodique >
! Je suis A l'ouest , j'envoie A destination de l'Est
      nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
      call mpi_issend(ar_send_est(1)              &
                ,size(ar_send_est(1:nb_send_est)) &
                ,mpi_real                         &
                ,par%tvoisin(est)                 &
                ,tagest_mpi                       &
                ,par%comm2d                       &
                ,tabreq_est_mpi   (nexchg_   )    &
                ,ierr)
      endif                                ! suplement cas periodique >
      endif                                       !ooo>

      if (par%tvoisin(ouest) /= mpi_proc_null) then !ooo>
      if(par%coords(1)/=0           ) then ! suplement cas periodique >
! Je suis A l'Est , je recois en provenance de l'Ouest
      nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
!     call mpi_irecv(ar_recv_nord(1:k)  &
      call mpi_irecv(ar_recv_ouest(1)                &
               ,size(ar_recv_ouest(1:nb_recv_ouest)) &
               ,mpi_real                             &
               ,par%tvoisin(ouest)                   &
               ,tagest_mpi                           &
               ,par%comm2d                           &
               ,tabreq_est_mpi   (nexchg_   )        &
               ,ierr)
      endif                                ! suplement cas periodique >
      endif                                         !ooo>

      if(nexchg_   >nexchgmax_mpi   ) &
      stop 'ar_obc nexchg_   >nexchgmax_mpi   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                  &
                      ,tabreq_est_mpi   (1:nexchg_   ) &
                      ,status_est_mpi  (:,1:nexchg_ ) &
                      ,ierr)

      if (par%tvoisin(ouest) /= mpi_proc_null) then !pmx>
      if(par%coords(1)/=0           ) then ! suplement cas periodique >
       if(allocated(ww3phase)) then !m°v°m>
        k=0 ; i=0
        do j=0,jmax+1 ; do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end
         k=k+1
         ww3phase(i,j,freq,dirw)=ar_recv_ouest(k)
        enddo ; enddo ; enddo
       endif                        !m°v°m>
      endif                                ! suplement cas periodique >
      endif                                         !pmx>
     

!     deallocate(ar_send_est)
!     deallocate(ar_recv_ouest)
      end subroutine wave_mpi_eastward

!..............................................................................

      subroutine wave_mpi_southward(j1_,j2_)
      implicit none
      integer,parameter :: nexchgmax_mpi=2 &
                          ,tagsud_mpi      =1616  
      integer :: nexchg_=0,nb_send_sud,nb_recv_nord,j1_,j2_
      integer,dimension(nexchgmax_mpi   ) :: tabreq_sud_mpi
      integer,dimension(mpi_status_size,nexchgmax_mpi   ) :: status_sud_mpi
!     real,dimension(:),allocatable :: ar_send_sud,ar_recv_nord

      nb_send_sud=(imax+2)*(dirw_end-dirw_beg+1)*(freq_end-freq_beg+1)
      nb_recv_nord=nb_send_sud

      if(.not.allocated(ar_send_sud)) allocate(ar_send_sud (nb_send_sud))
      if(.not.allocated(ar_recv_nord))allocate(ar_recv_nord(nb_recv_nord))

      ar_recv_nord=0.
      ar_send_sud=0.
      if(allocated(ww3phase)) then !m°v°m>
       k=0 ; j=j1_ 
       do i=0,imax+1 ; do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end
        k=k+1
        ar_send_sud(k)=ww3phase(i,j,freq,dirw)
       enddo ; enddo ; enddo
      endif                        !m°v°m>
      

!.............................................................................................
! Les domaines indiquent à leurs voisins le nombre d'element s'appretant à traverser les
! frontieres. 
       nexchg_=0

      if (par%tvoisin(sud) /= mpi_proc_null) then
! Je suis au nord , j'enverrai A destination du sud
       nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
       call mpi_issend(nb_send_sud,1,mpi_integer,par%tvoisin(sud)  &
              ,tagsud_mpi   ,par%comm2d,tabreq_sud_mpi   (nexchg_   ),ierr)
      endif

      if (par%tvoisin(nord) /= mpi_proc_null) then
! Je suis au sud, je recois en provenance du nord
      nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
       call mpi_irecv(nb_recv_nord,1,mpi_integer,par%tvoisin(nord)  &
              ,tagsud_mpi   ,par%comm2d,tabreq_sud_mpi   (nexchg_   ),ierr)
      endif

      ! La salle d'attente
      call mpi_waitall(nexchg_                  &
                      ,tabreq_sud_mpi   (1:nexchg_   ) &
                      ,status_sud_mpi  (:,1:nexchg_ ) &
                      ,ierr)


! Echange
      nexchg_=0
      if (par%tvoisin(sud) /= mpi_proc_null) then
! Je suis nord , j'envoie A sud
      nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
      call mpi_issend(ar_send_sud(1)              &
                ,size(ar_send_sud(1:nb_send_sud)) &
                ,mpi_real                         &
                ,par%tvoisin(sud)                 &
                ,tagsud_mpi                       &
                ,par%comm2d                       &
                ,tabreq_sud_mpi   (nexchg_   )        &
                ,ierr)
      endif

      if (par%tvoisin(nord) /= mpi_proc_null) then
! Je suis sud , je recois nord
      nexchg_   =nexchg_   +1
      if(nexchg_ >nexchgmax_mpi)stop 'nexchg_>nexchgmax_mpi'
!     call mpi_irecv(ar_recv_nord(1:k)  &
      call mpi_irecv(ar_recv_nord(1)                &
               ,size(ar_recv_nord(1:nb_recv_nord)) &
               ,mpi_real                             &
               ,par%tvoisin(nord)                   &
               ,tagsud_mpi                           &
               ,par%comm2d                           &
               ,tabreq_sud_mpi   (nexchg_   )            &
               ,ierr)
      endif

      if(nexchg_   >nexchgmax_mpi   ) &
      stop 'ar_obc nexchg_   >nexchgmax_mpi   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                  &
                      ,tabreq_sud_mpi   (1:nexchg_   ) &
                      ,status_sud_mpi  (:,1:nexchg_ ) &
                      ,ierr)

      if (par%tvoisin(nord) /= mpi_proc_null) then !pmx>
       if(allocated(ww3phase)) then !m°v°m>
        k=0 ; j=j2_ ! jmax+1
        do i=0,imax+1 ; do dirw=dirw_beg,dirw_end ; do freq=freq_beg,freq_end
         k=k+1
         ww3phase(i,j,freq,dirw)=ar_recv_nord(k)
        enddo ; enddo ; enddo
       endif                        !m°v°m>
      endif                                        !pmx>

!     deallocate(ar_send_sud)
!     deallocate(ar_recv_nord)
      end subroutine wave_mpi_southward

!..............................................................................
      
      subroutine wave_ww3_spec_sete_driver
      implicit none

! Le probleme: Comme les tableaux sshobc etc sont une somme de plusieurs composantes spectrales
! il ne faut pas passer 2 fois sur le meme point de grille (i,j) d'ou le capilotractage
! sur jstr_

      if(modulo(iteration3d,wavemodulo)==0) then !modulomodulo>

!      if(timewavenext==3) then !>>>
!        timewavenext=0 ; timewaveprev=3
!      else                     !>>>
!        timewavenext=3 ; timewaveprev=0
!      endif 
!      if(iteration3d==0)timewaveprev=timewavenext

! sachant les valeurs initiales:
!timewavenext=4  
!timewaveprev=3  
!timewaveprev2=0  
!timewaveprev3=-1   
                   k0=timewaveprev3
        timewaveprev3=timewaveprev2
        timewaveprev2=timewaveprev 
        timewaveprev =timewavenext
        timewavenext =k0

!       write(666,*)'4 temps',timewaveprev3,timewaveprev2,timewaveprev,timewavenext 

      endif                                      !modulomodulo>
      

!     if(imax+par%timax(1)==iglb) then
!     call wave_ww3_spec_sete_obc(imax-3,imax,0,jmax+1,2)
!     endif
!     if(1   +par%timax(1)==1   ) then
!     call wave_ww3_spec_sete_obc(1,4,0,jmax+1,1)
!     endif
      if(jmax+par%tjmax(1)==jglb) then
         call wave_ww3_spec_sete_obc(0,imax+1,jmax-3,jmax,4) !jmax-2 si vel_v(j6)...
      endif

!     call graph_out
!     stop 'coco6'

!     if(iteration3d==1) then
!     call graph_out
!     stop 'ALBERT'
!     endif

      end subroutine wave_ww3_spec_sete_driver

!--------------------------------------------------------------------------------

#ifdef key_oasis_symp_ww3
      subroutine wave_driver_oasis_symp_ww3
      use module_cpl_ww3
      use module_parallele
      implicit none
   
      il_oasis_time=int((elapsedtime_now-elapsedtime_rst)*10) ! To deciseconds
      
      wave_cpl_ww3_sdir=0 ; wave_cpl_ww3_cdir=0 ; wave_cpl_ww3_hs=0
      wave_cpl_ww3_hsw=0 ; wave_cpl_ww3_foc=0 ; wave_cpl_ww3_tw=0
      wave_cpl_ww3_tawx=0 ; wave_cpl_ww3_tawy=0 ; wave_cpl_ww3_twox=0
      wave_cpl_ww3_twoy=0 ; wave_cpl_ww3_uss=0 ; wave_cpl_ww3_vss=0
      wave_cpl_ww3_msk=0 
      
      call rcv_fields_from_ww3(il_oasis_time) !03-04-15
      
      if(wave_cpl_ww3_sdir*wave_cpl_ww3_cdir*wave_cpl_ww3_hs*wave_cpl_ww3_hsw &
         *wave_cpl_ww3_foc*wave_cpl_ww3_tw*wave_cpl_ww3_tawx*wave_cpl_ww3_tawy &
         *wave_cpl_ww3_twox*wave_cpl_ww3_twoy*wave_cpl_ww3_uss*wave_cpl_ww3_vss &
         *wave_cpl_ww3_msk == 1) then
         !** We received fields **!  
         
         wave_cpl_period_iter=wave_cpl_nextrec 
         wave_cpl_nextrec=1

         ! Apply the post reception traitement (wetmask from ww3, threshold, rotation, MPI)
         call wave_post_reception_traitement

         ! Compute the wave vector norm from the wave period
         call wave_kn_from_period

         ! Deduire le courant de Stokes 3D du courant 2D et du vecteur d'ondes:
         call wave_us3d_from_uss

         ! Deduire la fonction J barotrope du vecteur d'onde et de la hauteur significative:
         call wave_j_from_kn_hs

      else
         !** We didn't receive fields **!
         wave_cpl_nextrec=wave_cpl_nextrec+1
      endif

! Wave parameters at intermediate time steps : linear in time interpolation:
!      call wave_notime_interp_cpl
      call wave_linear_interp_cpl

! depth-averaged Stokes current:
      call wave_zaveraged_velstokes

      if(iwve==1) then !m°v°m> !17-11-19
! Note: on ne passe par ces etapes dans le cas iwve==2 (c.a.d. McWilliams99)

! Add the wave dissipation by bottom friction (uchiyama et al. 2010 - ums10)
      call wave_bottom_momentum_production

! Add the wave stress balance (tao-taw) to the wind stress at "velocity grid nodes": !21-05-10
      call wave_surface_momentum_production

! Compute the "Sj & Sshear" terms:
      call wave_sj_sshear

      endif            !m°v°m>     !17-11-19

      end subroutine wave_driver_oasis_symp_ww3
#endif

!--------------------------------------------------------------------------------

      subroutine wave_hamilton_ebersole !02-05-20
      use module_parallele
      implicit none


       if(iteration3d==0) then !000000000>
          allocate(sshwave_j_w(0:jmax+1     ,1,2:2,2)) ; sshwave_j_w=0.
          allocate(qwave_j_w  (0:jmax+1,kmax,1,2:2,2)) ;   qwave_j_w=0.
         allocate(velwave_j_u (0:jmax+1,kmax,1,2:2,2)) ; velwave_j_u=0.
        allocate(kvectorpeak_j(jmax,1:2))              ; kvectorpeak_j=kvectorpeak

!       i=imax/2
!       j=jmax/2
!       do k=kmax,1,-1
!        write(10+par%rank,*)depth_t(i,j,k) &
!               ,h_w(i,j)*sig1dpom_t(i,j,k)+ssh_int_w(i,j,1)*(1.+sig1dpom_t(i,j,k))
!       enddo


       endif                   !000000000>

      if(obcstatus(ieq1)==1) then                !----------->

      loop1=1

      x0=min(1.,elapsedtime_now/(periodpeak*10.))       ! RAMPE SUR 1 PERIODE
      x1=kvectorpeak*cos(pi/180.*10.)             ! kvector "Oi" (angle de 10°)
      x2=x0*0.08                                  ! Amplitude SSH (avec rampe)
      x3=grav*x1/freq2pipeak                      ! grav*Kvector_oi/omega

! Elevation de la surface  
      i=1
       do j=1,jmax

        sshwave_j_w(j,1,2,loop1)= &

             x2*sin( & !Sinus>
                     freq2pipeak*iteration3d*dti_fw                         &
                    -2.*pi*real(modulo(j+par%tjmax(1),jglb-2))/real(jglb-2) &
                   )   & !Sinus
            -x0*0.0175

!       if(par%rank==0.and.j==jmax/2)write(10+par%rank,*)elapsedtime_now &
!           ,sshwave_j_w(j,1,2,loop1),ssh_int_w(i,j,1)

       enddo

       i=1
       do j=1,jmax

! X5 est l'elevation de la surface de forcage "wave" (sans la moyenne)
          x5=x2*sin( & !Sinus>
                     freq2pipeak*iteration3d*dti_fw                         &
                    -2.*pi*real(modulo(j+par%tjmax(1),jglb-2))/real(jglb-2) &
                   )   !Sinus

       do k=1,kmax

! x4 est depth_t si ssh=sshwave
!       x4=h_w(i,j)*sig1dpom_t(i,j,k)+(1.+sig1dpom_t(i,j,k))*sshwave_j_w(j,1,2,loop1)
        x4=h_w(i,j)*sig1dpom_t(i,j,k)+(1.+sig1dpom_t(i,j,k))*x5

        qwave_j_w(j,k,1,2,loop1)= &

        grav*x2*sin( & !Sinus
                     freq2pipeak*iteration3d*dti_fw                         &
                    -2.*pi*real(modulo(j+par%tjmax(1),jglb-2))/real(jglb-2) &
                   ) & !Sinus

               *( & !ooo>  
!                   cosh(kvectorpeak*(depth_t(i,j,k)+h_w(i,j))) &
                    cosh(kvectorpeak*(     x4+h_w(i,j))) &
                   /cosh(kvectorpeak*max(0.01,h_w(i,j))) &
                    -1. &
                )   !ooo>

!       if(par%rank==0.and.j==jmax/2.and.k==kmax)write(10+par%rank,*)real(elapsedtime_now) &
!       ,real(qwave_j_w(j,1,1,2,loop1)/grav),real(qwave_j_w(j,k,1,2,loop1)/grav)  &
!       ,real(q_t(i,j,1,1)/grav),real(q_t(i,j,kmax,1)/grav)


       enddo
       enddo

       i=2
       do j=1,jmax

! X5 est l'elevation de la surface de forcage au point de vitesse u et au temps de
! la vitesse u dont on va se servir pour calculer x4 = "depth_u"
! correspondant
          x5=x2*sin( & !Sinus
                     freq2pipeak*(iteration3d*dti_fw-0.5*dti_fw)              &
                    -2.*pi*real(modulo(j+par%tjmax(1),jglb-2))/real(jglb-2)   &
                    -x1*0.5*dxb                                               &
                   )   !Sinus

       do k=1,kmax

! x4 est depth_u si ssh=forcing ssh
        x4=h_w(i,j)*sig1dpom_t(i,j,k)+(1.+sig1dpom_t(i,j,k))*x5

        velwave_j_u(j,k,1,2,loop1)= & 

             x3*x5                                   &  ! grav*kvector_oi/omega*ssh
               *cosh(kvectorpeak*(     x4+h_u(i,j))) &
               /cosh(kvectorpeak*max(0.01,h_u(i,j)))

       enddo
       enddo


      endif

      end subroutine wave_hamilton_ebersole


#endif
!--------------------------------------------------------------------------------

      subroutine wave_freq2kvector_acm(h_,period_,kvector_) !26-09-19
      implicit none
      double precision :: h_,period_,kvector_,delta_pulsation_  &
                         ,cmin_,cmax_,c_,mu_,pulsation_,pulsation2_,kvector_acm_
     

! Formule 24 Marsaleix et al OM 2019

      pulsation_=2.*pi/period_

! trouver kvector_ACM
      k10=100
      delta_pulsation_=1.e10
      cmin_=0.0001 ; cmax_=sqrt(grav*hmax)
      do loop2=1,4
      flag=0
      do loop1=0,k10

       c_=cmin_+real(loop1)/real(k10)*(cmax_-cmin_)
       mu_=sqrt( 1-(c_/acm_speed)**2 )
       kvector_=pulsation_/c_

! Formule 24 Marsaleix et al OM 2019:
        pulsation2_=kvector_**2*grav*(1.-(nhpgf_reduce/mu_**2))*h_  &
             +(nhpgf_reduce/mu_**3)*grav*kvector_*tanh(mu_*kvector_*h_)

       if(pulsation2_>pulsation_**2)cmin_=c_
       if(pulsation2_<pulsation_**2.and.flag==0) then
        cmax_=c_ ; flag=1
       endif

       if( abs(pulsation2_-pulsation_**2)<delta_pulsation_) then !>>>
        delta_pulsation_=abs(pulsation2_-pulsation_**2)
        kvector_acm_=kvector_
       endif                                                     !>>>

      enddo ! loop1
!     write(6,*)'ACM cmin_,cmax_-cmin_',cmin_,cmax_-cmin_
      enddo ! loop2


      kvector_=kvector_acm_


      end subroutine wave_freq2kvector_acm

!--------------------------------------------------------------------------------

      subroutine wave_freq2kvector(h_,period_,kvector_) !25-03-19
      implicit none
      double precision beta_,delta_c_,alpha_,c_,h_,kvector_,period_


      beta_=2.*pi*max(h_,0.01d0)/period_
      alpha_=grav/(2.*pi)*period_                ! c=g*t/2pi first guess for large depth

      if(beta_/alpha_>200.) then !fgfgfgfgfgfgfg>      !07-05-11

! Cas c=first guess valide:
       c_=alpha_ ! phase speed for large depth

      else                       !fgfgfgfgfgfgfg>

! Cas calcul itératif:
       delta_c_=0.01                  ! dx, on estime une première valeur de dx
       c_=alpha_+delta_c_             ! alpha_+dx
         do loop1=1,200

          delta_c_=(c_-alpha_*tanh(beta_/c_))                & !  (x-f(x))
                  /(-alpha_*beta_/(c_*cosh(beta_/c_))**2-1.)   !  /(f'(x)-1)
          c_=max(c_+delta_c_,1.d-4)
          if(abs(delta_c_)<0.01)goto 25
         enddo
         write(6,*)c_,delta_c_
         stop 'waveforcing: le calcul de kn ne converge pas'
   25  continue
      endif                      !fgfgfgfgfgfgfg>

      kvector_=2.*pi/(c_*period_)

      end subroutine wave_freq2kvector

!--------------------------------------------------------------------------------
!.....................................................................

      subroutine wave_init_ww3plugs !04-11-20
      implicit none
      integer(kind=1) :: loop_,flag_restore_zero_=0

! Cette routine calcule ww3plugs, c.a.d. les indices des trous et de leur boucheur

!.....................
! Identifier les trous
!......
! etape 1: attention si ww3filval=0 (valeur incompatible avec l'algo):
      flag_restore_zero_=0 
      if(ww3filval==0.) then !-changer-temporairement-ww3filval->
       ww3filval=-9999. ; flag_restore_zero_=1 
        do j=1,ww3_jmax ; do i=1,ww3_imax
         if(ww3_var(i,j)==0.)ww3_var(i,j)=ww3filval
        enddo           ; enddo
      endif                  !-changer-temporairement-ww3filval->
!......
! etape 2: 
      do j=0,jmax+1 ; do i=0,imax+1
       if(mask_t(i,j,kmax)==1) then !lsm>
         deci=ij2ww3_i(i,j) ; i1=int(deci) ; rapi=deci-i1 ; i1=i1-ww3zoom_istr+1
         decj=ij2ww3_j(i,j) ; j1=int(decj) ; rapj=decj-j1 ; j1=j1-ww3zoom_jstr+1
         if(ww3_var(i1  ,j1  )==ww3filval)ww3_var(i1  ,j1  )=-ww3filval
         if(ww3_var(i1  ,j1+1)==ww3filval)ww3_var(i1  ,j1+1)=-ww3filval
         if(ww3_var(i1+1,j1  )==ww3filval)ww3_var(i1+1,j1  )=-ww3filval
         if(ww3_var(i1+1,j1+1)==ww3filval)ww3_var(i1+1,j1+1)=-ww3filval
       endif                        !lsm>
      enddo ; enddo

!.....................
! Identifier les boucheurs
      do 3645 loop_=0,1 
      if(loop_==1)allocate(ww3plugs(ww3misspointmax,4)) 
      ww3misspointmax=0    
! Premier passage (loop_=0) evaluation d'une nbre de trous
! Allocation du tableau A la dimension ww3misspointmax
! Second passage pour renseigner le tableau de relation (i,j) trou vers (i,j) boucheur

      flag_stop=0
      do 1862 j=1,ww3_jmax
      do 1862 i=1,ww3_imax

       if(ww3_var(i,j)==-ww3filval)then !eeeeeeeeeeeeeee>
        ksecu=0
        i10=1
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=j-i10
         j2=j+i10
         j3=k1*(j2-j0)+(1-k1)

         i0=i-i10+k1
         i2=i+i10-k1
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3

         flag_stop=0
         if(i1<1       ) flag_stop=1 
         if(i1>ww3_imax) flag_stop=1 
         if(j1<1       ) flag_stop=1 
         if(j1>ww3_jmax) flag_stop=1 
         if(flag_stop==1) then !pmx> 
          write(10+par%rank,*)'WW3 grid overflow in finding a land replacement value'
          write(10+par%rank,*)'par%rank=',par%rank
          write(10+par%rank,*)'WW3 grid relative to extraction zone i,j:',i,j
          write(10+par%rank,*)'is looking for fill value in        i1,j1:',i1,j1
          write(10+par%rank,*)'WW3 grid global i ,j :                   ', i+ww3zoom_istr-1,j +ww3zoom_jstr-1
          write(10+par%rank,*)'WW3 grid global i1,j1:                   ',i1+ww3zoom_istr-1,j1+ww3zoom_jstr-1
          write(10+par%rank,*)'ww3 local imin imax:',1,ww3_imax
          write(10+par%rank,*)'ww3 local jmin jmax:',1,ww3_jmax
          write(10+par%rank,*)'Extraction zoom coordinates in the ww3 global file:'
          write(10+par%rank,*)'Min Max i:',ww3zoom_istr,ww3zoom_iend
          write(10+par%rank,*)'Min Max j:',ww3zoom_jstr,ww3zoom_jend
          write(10+par%rank,*)
          if(i1<1)       write(10+par%rank,*)'PROBLEM: local i1<1'
          if(i1>ww3_imax)write(10+par%rank,*)'PROBLEM: local i1>ww3_imax'
          if(j1<1)       write(10+par%rank,*)'PROBLEM: local j1<1'
          if(j1>ww3_jmax)write(10+par%rank,*)'PROBLEM: local j1>ww3_jmax'
          write(10+par%rank,*)'SOLUTION: increase the value of ww3_enlarged_bounds' &
                             ,' in module_wave.F90'
          write(10+par%rank,*)
          goto 3678
         endif                 !pmx>
! si flag_stop=1 c'est qu'il n'y a pas de bouchage
! suffisament proche du trou et que l'algo cherche au dela de la zone d'extraction,
! ce qu'on ne permet pas pour avoir une parallelisation parfaite.
! Ce qu'on peu faire: on peut augmenter la taille de la
! zone d'extraction en augmentant ww3_enlarged_bounds. Mais cela revele surtout
! que les masques oceano et ww3 sont trop differents et qu'il faut sans doute
! ameliorer le masque oceano

          if(ww3_var(i1,j1)/=-ww3filval.and.  &
             ww3_var(i1,j1)/= ww3filval     )then !%%%%%%%%%%%%%>

           dist2=sqrt(real(i-i1)**2+real(j-j1)**2)
           ksecu=1
             if(dist2.lt.dist1) then              !>>>>>>>>>>>>>
              i4=i1
              j4=j1
              dist1=dist2
             endif                                !>>>>>>>>>>>>>

          endif                                   !%%%%%%%%%%%%%>
 1861    continue

 1863   continue
        i10=i10+1
        if(ksecu.eq.0)goto 1864
        ww3misspointmax=ww3misspointmax+1
        if(loop_==1) then !m°v°m>   
           ww3plugs(ww3misspointmax,1)=i
           ww3plugs(ww3misspointmax,2)=j
           ww3plugs(ww3misspointmax,3)=i4
           ww3plugs(ww3misspointmax,4)=j4
        endif             !m°v°m>

       endif                           !eeeeeeeeeeeeeee>
 1862 continue

 3645 continue ! boucle loop_

! si ww3filval a ete modifie (car nul initialement) alors retablir 0:
      if(flag_restore_zero_==1) then !--retablir-ww3filval-->
        do j=1,ww3_jmax ; do i=1,ww3_imax
         if(ww3_var(i,j)==ww3filval)ww3_var(i,j)=0.
        enddo           ; enddo
        ww3filval=0.
      endif                          !--retablir-ww3filval-->

  3678   call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ; if(k0/=0)stop 'Err 8125 grid overflow in module_wave. See fort.xxx error files' 

      end subroutine wave_init_ww3plugs

!.....................................................................

      subroutine wave_appli_bouchetrou !04-11-20
      implicit none

      do k=1,ww3misspointmax 
       ww3_var(ww3plugs(k,1),ww3plugs(k,2))=  & 
       ww3_var(ww3plugs(k,3),ww3plugs(k,4))
      enddo

      end subroutine wave_appli_bouchetrou

!.......................................................

      subroutine wave_get_l(t_) !10-11-21
      implicit none
      integer t_,time_
      real*8 valmin_
#ifdef synopsis
       subroutinetitle='wave_get_t'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call wave_get_time_from_binrecfile(t_,time_)

      ww3varname(1)='lm' 
      ww3varname(2)='lm' 
      ww3varname(3)='lm'
      call wave_read_only(time_)

      ww3mskval=1./var_scalefactor  ; valmin_=1. 
      call wave_interp_only('basic')

      do j=0,jmax+1 ; do i=0,imax+1
        k_wave_t(i,j,2)=2.*pi/max(anyvar2d(i,j),valmin_)
      enddo ; enddo

      end subroutine wave_get_l

!.......................................................

      end module module_wave
