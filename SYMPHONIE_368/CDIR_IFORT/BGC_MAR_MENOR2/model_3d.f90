










      subroutine model_3d
!______________________________________________________________________
! SYMHPHONE ocean model
! release 296 - last update: 27-02-21
!______________________________________________________________________
      use module_principal ; use module_parallele ; use module_wave
      use module_drifter ; use module_airseaflux ; use module_offline
      use module_modeanalysis ; use module_external_mode ; use module_s
      use module_my_outputs ; use module_cpl_oasis
! THOM SEDIM
      use module_parameter_sedim, only : l_sedim_s
      implicit none
      character*60 txt_spy_                     !27-05-11

!______________________________________________________________________
! Version date    Description des modifications
!         09/08/01: prise en compte des amenagements sur river_upd.F
!         12/08/01: amenagements sur le protocole de l'imbrication
!         13/08/01: amenagements sur le protocole de l'imbrication
!         21/08/01: amenagements sur le protocole de l'imbrication
!         25/08/01: appels à strada
!         06/12/01: appels à tkenergy_xing (en preparation)
!         15/07/02: amenagemenst pour introduction de la filiere
!                   "formules bulk" pour calculer les flux atmospheriques
!         23/08/02: mise à jour à partir de la version 2003. L'appel
!                   à la marée est réintroduit pour être promu 2003
!         10/09/02: amenagement d'un cas ichoix=1 pour update_obcforcingterms
!         07/10/02: appel à CALL BOUEES(1)
!         27/12/02: bienvenue à externe_fast
!         27/01/03: appel à la routine offline CALL OFFLINE_INOUT(1)
!                   appel à DATE_OUTPUT
!                   appel au schema de turbulence de Xing implémenté pour
!                   modele par Guillaume Reffray du LSEET.
!         25/03/03: Suppression de CALL ANIMATION
!         11/07/03: ajout d'un sous programme pour faire des sorties
!                   personnelles: CALL MY_OUTPUTS,
!                   suppressions appel DYNAMO et condition sur KSERIE
!         28/07/03: ajout d'un commentaire
!         29/07/03: bienvenue à CALL MOYENNE_TEMPS
!         05/08/03: il y avait un commentaire devant CALL THE_END(1)
!                   et il est maintenant enlevé.
!         28/08/03: Appel à NEST_INOUT(1)
!         31/08/03: Appel à NEST_INOUT(2)
!         29/09/03: appel à ANALYSEHARMONIQUE(1)
!         03/11/03: bilan d'energie pour turbulence
!         04/11/03: passage d'un argument dans CALL DENSITY
!                   + reorganisation: appel de tkenergy apres
!                   scalars pour connaitre T et S au temps t+1
!                   afin de calculer correctement le terme de
!                   flottabilité de tke pour être conservatif
!                   avec le melange 1DV de T et S. Note que ce
!                   déplacement ne pose plus de pb avec la filière
!                   forward car OMEGA depuis quelques mois OMEGA
!                   possede un 4eme argument qui permet de stocker
!                   les vitesses verticale "forward"
!         09/01/04: cas IAIRSEA=3
!         03/06/04: passage d'un argument dans my_outputs
!         08/07/04: routines non compatibles commentées
!                   & introduction routines "ondelettes" de Francis
!         14/11/04: possibilite equation d'etat non lineaire:
!         13/04/05: faire un fichier restart "chanel9" sans arreter
!                   la simulation
!         25/04/05: Possibilié de sortir des chanel9 intermediaires
!                   avec des noms datés.
!                   Ajout à la fin du fichier d'un sous programme
!                   pour un cas "année perpetuelle"
!         26/04/05: demander un fichier restart sans arreter le model_
!         17/01/06: ajout des forces de la houle (these Cléa Denamiel)
!                   l'appel au dragcoef est déplacé apres call vhz_to_vel
!                   pour eventuellement pouvoir appliquer les effets de
!                   la houle (et pour cela il nous faut vel(t) )
!         14/02/06: un argument apparait dans l'appel à dragcoef
!         02/07/06: La diffusion korizontale est appliquee sur les vitesses
!                   prises au temps t-1. Ceci entraine l'appel a la routine
!                   vishor apres la routine asselin_upd
!         15/09/06: Suite du point precedent. Separation des mises à jours
!                   avec d'un coté les vitesses et de l'autre les autres
!                   variables. Coef de diffusion pour les vitesses non
!                   constant
!         27/03/07: L'appel à dtvarobc est déplacé de sorte que le barometre
!                   inverse soit activé même si il n'y a que la meteo comme
!                   forcage externe
!         03/04/07: filiere model_in model_out definitivement supprimee
!                   Appel à chronos placé apres nest_inout et offline_inout
!         07/06/07: bornes min et max pour i et j passees en argument de la
!                   subroutine z_averaged
!         26/09/07: regroupement des ecritures des "sorties" en fin de cycle
!         17/01/08: IWAVE renommé IWVE
!         23/01/08: suppression routine periodic year
!         20/03/08: ajout 3eme arg dans vhz_to_vel pour ecologie offline
!                   ajout de la subroutine pour assimilation wp3
!         29/12/08: appel (par defaut commenté) à sse_correction.F pour
!                   projet ECOOP
! 2009.2  14-06-09  securite pour empecher d'appeler hot_restart quand ce
!                   n'est pas le moment
!         01-09-09  externe_fast renomme externe
!         16-09-09  Synchro des proc à la fin de chaque iteration
! 2009.3  30-09-09  appel à cellbox_thickness
!         01-10-09  l'appel à l'albedo est placé dans model_3D
!         02-10-09  vhz_to_vel renommé veldxdz_veldydz_to_vel
!         03-10-09  - internal_mode.F90 devient internal_mode.F
!                   - externe.F devient external_mode.F
!                   - tkenergy.F devient turbulent_kinetic_energy.F
!         05-10-09  ajout d'un "ifdef parallele"
!         08-10-09  advection verticale lagrangienne: deplacement de
!                   l'appel à cellbox_thickness
!         11-10-09  suppression tshz_to_ts
!         19-10-09  des arguments passés dans cellbox_thickness
! 2010.2  27-12-09  obc_af renommé read_ogcm_fields
! 2010.3  06-01-10  appel à wetdry_mask_airseafluxes
!         16-01-10  - suppression routine mixingcoef
!                   - routine difver_upd renommée vertmix_coef
! 2010.7  15-02-10  momentum_equation compte maintenant 2 parties:
!                   avant et apres le mode externe
!         16-02-10 move_ts.F90 renommé moveforward.F90
!         25-02-10 l'appel à model_wave se fait avec les autres forcages
!         28-02-10 supression du sous-programme de passage de vel"tilde"
!                  à vel
! 2010.8  12-03-10 mixer le time filter 3L et 4L avec tfc1 et tfc2
!         03-05-10 ajout routine stokesvortexforce & supression appel
!                  à momentum_equations
!         09-05-10 seul le proc 0 ecrit interrupteur
! 2010.9  01-06-10 - faire les archives graphiques avant moveforward
!                  - l'actualisation de dz passe dans moveforward
!                  autrement dit appel à cellbox(2) supprimé
!         06-06-10 model_wave renomme waveforcing
! 2010.10 21-06-10  operations sur interrupteur regroupees dans io_switch
!         23-06-10  hot_restart renommé dyn_restart et bio_restart
! 2010.11 23-07-10  - appel à equation_of_state (ex density) dans presgrad
!                   - les longueurs turbulentes sont deplacées apres
!                     les traceurs car rhp doit contenir la densité
!                     potentielle
!         05-08-10 presgrad renommé pressure_gradient
! 2010.13 31-10-10 amelioration des tests de conservation de la parallelisation
!         03-11-10 memoriser la date de l'iteration en cours
! 2010.14 05-11-10 suppression des anciens tfc1 et tfc2
!         14-12-10 cas ichanel9=2 obsolete
!         16-12-10 Routines turbulence appelée depuis internal mode
! 2010.20 16-04-11 Calculs sur la base d'un temps en secondes
!         19-04-11 appel à time_step et à time_step_next
! 2010.22 27-04-11 ichanel9 remplacé par restartfileperiod
!         28-04-11 Ne pas ecrire le fichier restart si Nan ou INF dans la
!                  solution
! 2010.23 12-05-11 Debug ecriture fichier restart cas "single"
!         15-05-11 routine model3d_produce_restart n'a pas besoin
!                  d'ecrire dans interrupteur
!         27-05-11 - subroutine model3d_produce_restart renommée dyn_restart_ask
!                  cette derniere étant incluse dans fichier dyn_restart.F90
!                  - variable subroutine_spy permet de savoir qui appelle
!                  la barriere
! 2010.25 01-02-12 routine drifter remplace routine bouee
!         09-02-12 verification modele 1DV
!         21_02-12 - modif argument dans dragcoef
!                  - appel a waveforcing devient appel a wave
!         09-03-12 Afin de pouvoir regarder l'initialisation des vagues
!                  l'appel a graph_out quand iteration3d=0 est place dans model3d
! 2010.25 08-06-12  use module_drifter
! S25.4   30-06-12  use module_airseaflux
! S26.1   04-09-12  Si iteration2d_max_now=1 pas de mode externe
!         30-01-13  module_albedo
!         10-02-13  appel a airseaflux_driver
!         26-02-13  appel modeanalysis
!         03-03-13  module_external_mode et cie
!         08-03-13  ajout my_outputs_driver
!         01-04-13  arg passes dans call elapsedtimetodate
!         09-04-13  Mise a jours des profondeurs avant l'appel aux forcages
!         02-10-13  On ecrit le fichier graphique de l'etat initial apres que
!                   le gradient de pression ait ete calcule afin de pouvoir visualiser
!                   le courant geostrophique initial
!         12-11-13  suppression de la barriere (subcycling)
!         11-07-14  suprime nest_inout
!         09-09-14  call momentum_before_external !09-09-14
!         02-11-14  appel graph_out_bio
!         08-11-14  appel graph_out_bio depuis graph_out
!         11-11-14  appel a graph_out dependant de iteration3d_restart
!         17-11-14  suite du point precedent
!         23-01-15  the_end(1) passe dans main.F90 
!         25-05-15  calcul cpu_seconds
!         29-06-15  clefs oasis
!         01-10-15  calcul cpu_seconds detaille dans internal_mode
!         27-11-15  check mip 1dv
!         20-02-16  calcul cpu harmonics 3D
!         16-05-16  ajout de l'appel à mang3dto2d pour une mangrove quadratique
!         29-09-16  clef oasis
!         21-12-16  subroutine model_3d_anyv_debug permet de verifier que les tableaux
!                   generiques ne laissent pas de traces sur la solution du modele
!         23-10-17  cas sans time splitting
!         18-08-18  passage au "no-splitting demain"
!         23-08-18  if(iteration2d_max_now>0)call z_averaged(1,imax+1,1,jmax+1,1)
!         06-09-18  condition d'appel au pgf
! v251    10-04-19  call obc_river_groundwater
! v256    06-06-19  if(ioffline.eq.1)call offline_inout(1) est deplacE
!                   dans internal_mode.F90 (voir explications dans cette routine)
! v267    16-11-19  cas iwve==2
! v269    06-12-19  call momentum_bottom_current
! v287    27-07-20  call cpl_oasis_get_pmx !27-07-20
! v291    07-11-20  call update_obcforcingterms(1)     !07-11-20
! v292    19-11-20  iwve==3
! v296    19-02-21  if(flag_surfriver==1)call obc_river_surfriver_w !19-02-21
!         27-02-21  flag_omega_cumul=0
!...............................................................................
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

       txt_spy_   ='model_3d'

!******************************************************
! Begin the iterative process
!******************************************************
 9999 continue
      scpuk=0 ; call s_cpu('init',0)

      restart_file_y_or_n=0 ! You can not write a restart file !14-06-09
      flag_omega_cumul=0                                       !27-02-21


      call time_step ! Gives the internal time step

!     call elapsedtimetodate(elapsedtime_now)                           !19-04-11
! changer pour elapsedtimetodate!
!     year_now=i5 ; month_now=i6 ; day_now=i7
!     hour_now=i3 ; minute_now=i2 ; second_now=i1                       !03-11-10
      call elapsedtime2date(elapsedtime_now,                &          !01-04-13
                            year_now,month_now,day_now,     &
                            hour_now,minute_now,second_now)

! level depth:
      call z_levels(0)            ! deplace le 09-04-13

!________________________________________________________________!
! Update forcing terms:
!________________________________________________________________!

! air/sea fluxes:
      call airseaflux_driver(2)            !10-02-13
!$ River water discharge:
      call river_upd(2)                                                !09/08/01
      if(flag_groundwater==1)call obc_river_groundwater_w !10-04-19
      if(flag_surfriver==1)call obc_river_surfriver_w !19-02-21

!$ Tidal harmonic components:
      if(kmaxtide>0) then !---> !26-02-13
         call update_tide ; call s_cpu('update_tide',0)
         if(flag_3dwaves_harmonics==1) then !hhh>
           call modeanalysis_harmonics(1)   !26-02-13
           call s_cpu('3Dharmonics',0)      !20-02-16
         endif                              !hhh>
      endif               !--->

!  Outer ocean fields:                                                 !17/06/01
      if(iobc_wv.eq.1.or.iobc_ogcm.eq.1   &
                     .or.nest_onoff_in.eq.1) then !$$$$$$$$$$$$$$$$$>  !31/08/03
!  General case:
      if(iobc_ogcm==1.and.iobc_lr==0)call read_ogcm_fields(2)
!  Particular case of self nesting:
!     if(nest_onoff_in.eq.1.and.kount.ne.kount0)call nest_inout(2)     !31/08/03
!     if(nest_onoff_in.eq.1.and.iteration3d.ne.kount0)call nest_inout(2)     !31/08/03

      endif                                       !$$$$$$$$$$$$$$$$$>

!  Wave/current forcing terms:
     if(iwve==1.or. & !16-11-19
        iwve==2.or. &
        iwve==3)call wave_driver ! obligatoirement apres airseaflux_driver

!  Forcing terms in open boundary conditions:
      if(oasis_symsym_onoff==1) then  !-avec-oasis-symsym->
       call cpl_oasis_get_pmx !27-07-20
       call cpl_oasis_put_pmx
       call update_obcforcingterms(iobc_ogcm) !si oasis alors iobc_ogcm=0 et obc(1) donnE sur toute la grille 27-07-20
      else                            !-sans-oasis-symsym-> !07-11-20
       call update_obcforcingterms(1)                       !07-11-20
      endif                           !-sans-oasis-symsym->
      call s_cpu('forcingterms',0)

!________________________________________________________________!

!  Depth average of the internal_mode velocities: (inutile si sans mode splitting)
!     if(iteration2d_max_now/=0)call z_averaged(1,imax+1,1,jmax+1,1) !07/06/07 !23-08-18
      if(flag_timesplitting==1) call z_averaged(1,imax+1,1,jmax+1,1) !07/06/07 !23-08-18

!  Bottom drag coefficient:
      call dragcoef

! Bottom current (including merged layers case) velbot_u, velbot_v
      call momentum_bottom_current !06-12-19

!call model_3d_anyv_debug !21-12-16
                                                                       !14/02/06
!  pressure gradient force related to the density anomaly:
!  Conditions pour ne pas entrer dans le calcul du PGF hydrostatique:
!    1: mode non-hydrostatique NE calculant PAS T et S 
!  Conditions pour entrer dans le calcul du PGF hydrostatique:
!    1: Mode hydrostatique (flag_nh3d=0)
!    2: Mode non-hydrostatique avec prise en compte de T et S (flag_nh3d=2)
      if(flag_nh3d/=flag_nh3d_nosplit_uv) then  !m[0v0]m> !06-09-18
         call pressure_gradient ; call s_cpu('PGFhydrostatic',0)
      endif                                     !m[0v0]m>
! Non-hydrostatic component 
! Conditions pour entrer dans la routine d'ajout de la composante NH
!    1: Mode non-hydrostatique flag_nh3d/=flag_nh3d_none
! note: cette etape inclut eventuellement l'ajout du gradient de ssh si SANS mode splitting
      if(flag_nh3d/=flag_nh3d_none)call pressure_gradient_add_nhpgf

      if(iwve==1) then
         call stokesforces ; call s_cpu('stokesforces',0)
      endif

      if(iteration3d==iteration3d_restart)call graph_out !17-11-14

! Cas avec time-splitting:
!     if(iteration2d_max_now/=0) then  !-time-splitting-case-> !18-08-18
      if(flag_timesplitting==1)  then  !-time-splitting-case-> !18-08-18

         call adve3dto2d ; call s_cpu('adve3dto2d',0)
         if(coef_diss_mangrove>0.)call mang3dto2d ! frozen terms from 3d to 2d for quadratic mangrove friction terms
!  Surface elevation & barotropic velocities:
         call external_mode_driver ; call s_cpu('external',0)

! CommentE le 18-08-18 avec le passage au "no-splitting demain"
!     else                            !m[0v0]m>
! Compute transport and sea surface elevation in case of no time splitting
!        call internal_mode_nosplitting !23-10-17

!  Thickness of the whole ocean layer:
      call z_thickness(2,2)

!  computes the wetdry mask and cancels airsea fluxes in dried aeras:   !06-01-10
      call wetdry_mask_airseafluxes

      call s_cpu('z_thickness',0)

!  Baroclinic velocities, internal/external modes coupling:
      call internal_mode

      endif                           !-time-splitting-case->

!     if(iteration2d_max_now==0) then !-NO-time-splitting-case-> !18-08-18
      if(flag_timesplitting==0)  then !-NO-time-splitting-case-> !18-08-18
! Compute transport and sea surface elevation in case of no time splitting
         call internal_mode_nosplitting  
      endif                           !-NO-time-splitting-case->


!$ update location of lagrangian drifters                              !01-02-12
      if(drifter_onoff==1)call drifter_update

!     call s_cpu('internal',0)

!  Temperature & salinity:
! Note: on ne calcule pas T et S dans le cas NH basE sur u,v uniquement
        if(flag_nh3d/=flag_nh3d_nosplit_uv)call scalars !23-08-18

!$ Passive tracers:
      if(imodeltrc==1.or.  &
         imodelbio==1.or.  &
         l_sedim_s) call strada(1)               !25/08/01

!$ Post-processing: time averaged 3d fields for biogeochemical "off line" simulations:
!     if(ioffline.eq.1)call offline_inout(1) commetE le !06-06-19

!$ Post-processing: output for graph                                    !01-06-10
                          call date_output(idate_output)
      if(kpvwave.eq.1)    call date_output(3)

      call my_outputs_driver !08-03-13

! Couplage surfex
! Couplage ww3

      call s_cpu('scalars',0)

!$ update time:
      call time_step_next ! update elapsedtime
      call chronos

!$ move forward the 3d variables:
      call moveforward ; call s_cpu('moveforward',0)

!***************************************************
!  Fin d'un cycle de calcul. Test de Fin.
!***************************************************
      restart_file_y_or_n=1 ! You can write a restart file             !14-06-09

!.......................................................!
! SORTIE FICHIER RESTART                                !
!$ Produce a restart file:

!ATTENTION DOUBLE PRECISION OBLIGATOIRE
      if(restartfileperiod>0.)                                 & !15-05-11
      call dyn_restart_ask(restartfileperiod,'modulo')  !ATTENTION DOUBLE PRECISION OBLIGATOIRE !27-05-11
!ATTENTION DOUBLE PRECISION OBLIGATOIRE

      if(give_chanel9==1) then    !%%%%%%%%%%%%%%%%%%%%%>
                                       call dyn_restart('w')
       if(imodeltrc==1.or.imodelbio==1)call bio_restart('w')
       give_chanel9=0
       call io_switch('w') ! (interrupteur)                                 !21-06-10
      endif                       !%%%%%%%%%%%%%%%%%%%%%>

      if(elapsedtime_now>=elapsedtime_end.or.kstop==1) then   !ksksksksksksks> !16-04-11
       if(restartfileperiod>0.) then      !>>>>>>>>>>>>>>>                     !27-04-11
                                        call dyn_restart('w')
        if(imodeltrc==1.or.imodelbio==1)call bio_restart('w')
       endif                              !>>>>>>>>>>>>>>>
      endif                                                   !ksksksksksksks>
!.......................................................!

!call model_3d_anyv_debug !21-12-16

      if(elapsedtime_now<elapsedtime_end.and.kstop.ne.1)goto 9999      !16-04-11

      call the_end(0) ! Post-Processing (tides) + netcdf file production

      end subroutine model_3d

!-----------------------------------------------------------------
!----------------------------------------------------------------------
   subroutine rotation_wind_surfex
    use module_principal
    use module_parameter
    use module_parallele
! Rotation point "t" :
!
       do j=0,jmax+1
       do i=0,imax+1
        xy_t(i,j,3)=gridrotcos_t(i,j)*taux_w(i,j,1)   &
                   -gridrotsin_t(i,j)*tauy_w(i,j,1)
        xy_t(i,j,4)=gridrotsin_t(i,j)*taux_w(i,j,1)   &
                   +gridrotcos_t(i,j)*tauy_w(i,j,1)
       enddo
       enddo
!     
      call get_type_echange('za','xy_t_za_3',xy_t,lbound(xy_t),ubound(xy_t),3,i3) !04-10-14
      call get_type_echange('za','xy_t_za_4',xy_t,lbound(xy_t),ubound(xy_t),4,i4)
      do loop3=1,subcycle_exchange
        call echange_voisin(xy_t,i3,mpi_neighbor_list(loop3)) !31-07-14
        call echange_voisin(xy_t,i4,mpi_neighbor_list(loop3))
      enddo
     call loc_wait()
! Interpolation point "u":
       do j=1,jmax
       do i=1,imax+1
        wstress_u(i,j,1)=0.5*(xy_t(i,j,3)+xy_t(i-1,j,3))
       enddo
       enddo
! Interpolation point "v":
       do j=1,jmax+1
       do i=1,imax
        wstress_v(i,j,1)=0.5*(xy_t(i,j,4)+xy_t(i,j-1,4))
       enddo
       enddo

! Module de la tension de vent:!
      do j=0,jmax
      do i=0,imax

      wstress_w(i,j)=sqrt(                                              &
       ((wstress_u(i,j,1)+wstress_u(i+1,j,1))/2.)**2+                   &
       ((wstress_v(i,j,1)+wstress_v(i,j+1,1))/2.)**2)

      enddo
      enddo


     end subroutine rotation_wind_surfex
!----------------------------------------------------------------------...............................................................
     subroutine model_3d_anyv_debug !21-12-16
     use module_principal ; use module_parallele
     implicit none

      if(par%rank==0) then
        write(6,*)'PASSE PAR model_3d_anyv_debug'
        write(6,*)'huge(anyv3d)',huge(anyv3d)
        write(6,*)'huge(anyvar3d)',huge(anyvar3d)
      endif

! Modifier les valeurs des tableaux generiques ne doit pas changer la trajectoire du modele....

! Cette routine s'utilise en 2 temps. Une premiere simulation sans l'appel A cette routine.
! Une seconde avec l'appel. Les resultats des 2 simulations doivent etre identiques

     anyv3d=huge(anyv3d)
     anyvar3d=huge(anyvar3d)
     anyvar2d=huge(anyvar2d)
     anyv1d=huge(anyv1d)
     xy_t=huge(xy_t)
     xy_u=huge(xy_u)
     xy_v=huge(xy_v)
     xy_f=huge(xy_f)

!#ifdef bidon
! Verifier l'etancheIte a la solution masquee pour les traceurs
! attention suppose pas de rivieres...
      if(nriver==0) then !>>>>>>>
       do k=1,kmax
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,k)==0) then
         tem_t(i,j,k,:)=1.
         sal_t(i,j,k,:)=20.
        endif
       enddo
       enddo
       enddo

      if(allocated(tken_w)) then !pmxpmx>
       do k=1,kmax+1
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,k)==0) then
         tken_w(i,j,k)=0.01
         km_w(i,j,k)=0.01
         kh_w(i,j,k)=0.01
        endif
       enddo
       enddo
       enddo
      endif                      !pmxpmx>

      if(allocated(epsn_w)) then !pmxpmx>
       do k=1,kmax+1
       do j=1,jmax
       do i=1,imax
        if(mask_t(i,j,k)==0) then
         epsn_w(i,j,k)=1.d-7
        endif
       enddo
       enddo
       enddo
      endif                      !pmxpmx>

      endif             !>>>>>>
!#endif


     end subroutine model_3d_anyv_debug
