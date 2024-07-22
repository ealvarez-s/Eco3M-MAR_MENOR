










      subroutine model_offline
!______________________________________________________________________
! SYMPHONIE ocean model
! release 303 - last update: 20-11-20
!______________________________________________________________________

      use module_principal  ; use module_parallele ; use module_drifter
      use module_airseaflux ; use module_offline   ; use module_s
      use module_wave
! THOM SEDIM
      use module_parameter_sedim, only : l_sedim_s
      implicit none

!...............................................................................
!Version  date      Description des modifications:
!         27/01/03: mise en service
!                   et appel à CALL DATE_OUTPUT (nouveauté par rapport
!                                                à SYMPHONIE2002)
!         09/01/04: cas IAIRSEA=3
!         16/02/09: debug: si model_offline physique et nesting bio
!                   il ne faut pas oublier le calcul de rap_obc afin
!                   que l'evolution des champs externes dans obc_bio
!                   puisse etre calculee
! 2009.2  16-09-09: Synchro des proc à la fin de chaque iteration
! 2009.3  01-10-09: Introduction des modifs de Caro et Claude
!                   standard netcdf compatible avec Ariane
!                   ajout de la securite RESTART_FILE_Y_OR_N pour controler
!                   le moment d'appel de hot_restart.F
!                   L'appel au calcul de l'albedo est placé dans model_3D
!                   et dans model_offline
!         05-10-09: ajout d'un "ifdef parallele"
! 2010.4  22-01-10: - prise en compte des bancs decouvrants
!                   - prise en compte de la coordonnee verticale generalisee
! 2010.8  09-05-10  seul le proc 0 ecrit interrupteur
! 2010.10 16-06-10  La mise à jour des facteurs d'echelle est directement
!                   faite dans offline_inout et non plus calculée par
!                   appel aux subroutine z_thickness et cellbox_thickness
!         21-06-10  operations sur interrupteur regroupees dans io_switch
!         23-06-10  hot_restart renommé dyn_restart et bio_restart
! 2010.18 27-03-11  Calcul date
! 2010.20 16-04-11  Calculs sur la base d'un temps en secondes
! 2010.21 20-04-11  suite
! 2010.22 27-04-11  ichanel9 remplacé par restartfileperiod
! 2010.25 01-02-12  appel a routine drifter
!         08-06-12  use module_drifter
!         30-06-12  use module_airseaflux
! S26.1   24-11-12  module_offline
!         30-01-13  module_albedo
!         10-02-13  appel a airseaflux_driver
!         01-04-13  appel a nouvelle routine elapsetime2date
!         28-03-14  possibilite de graphique a l'etat initial
!         02-11-14  appel graph_out_bio
!         08-11-14  appel graph_out_bio depuis graph_out
!         07-07-15  ajout call dyn_restart_ask
!         31-03-16  suppression call mpi_barrier(par%comm2d,k_out) 
!         12-02-17  jalons cpu
!         14-06-17  ajout modelisation sedimentaire offline
!         20-06-17  ajout module de vague et drag coef
!         10-07-17  ajout module_wave
! v267    19-11-19  appel aux vagues
! v269    06-12-19  Reste A coder (si besoin) en rapport avec velbot_u et _v
! v292    20-11-20  iwve==3
! v303    05-08-21  suppression d'un stop
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


!******************************************************
! Debut de la boucle temporelle.
!******************************************************

 9999 continue

      scpuk=0 ; call s_cpu('RESET',0) !11-02-17

      restart_file_y_or_n=0 ! you can not write a restart file !01-10-09
      call time_step        ! Gives the internal time step
      call elapsedtime2date(elapsedtime_now,                &          !01-04-13
                            year_now,month_now,day_now,     &
                            hour_now,minute_now,second_now)

!________________________________________________________________!
! MISE A JOUR DES FORCAGES:                                      !
! DEBUT:                                                         !
!________________________________________________________________!
! Calcul du forçage atmospherique à partir d'un fichier météo:
      call airseaflux_driver(2)                       !10-02-13
      call s_cpu('flux atmospherique',0) !11-02-17

!  Wave/current forcing terms:
     if(iwve==1.or. & !19-11-19
        iwve==2.or. & !20-11-20
        iwve==3)call wave_driver ! obligatoirement apres airseaflux_driver

! lecture et interpolation temporelle des fichiers physique offline
      call offline_inout(2)
      call s_cpu('physique offline',0) !11-02-17

! Bottom current (including merged layers case) velbot_u, velbot_v
      call momentum_bottom_current !06-12-19
! Ce stop en suivant pour verifier que velbot n'est pas necessaire A la
! filiere offline (par ex pour la remise en suspension des sediments...'
!     write(6,*)'flag_merged_levels',flag_merged_levels
!     i=imax/2 ; j=jmax/2 ; k=kmin_u(i,j)
!     if(mask_v(i,j,kmax)==1)write(6,*)'k',k,real(vel_v(i,j,k,0:1)),real(velbot_v(i,j))
!     stop 'Mettre A jour le codage de velbot offline si besoin' !06-12-19 !05-08-21

      if(iteration3d==0) then !>>>
                       call graph_out     !28-03-14
      endif                   !>>>

! Si nesting bio:                                                      !16/02/09
!     if(nest_onoff_in.eq.1)                               &
!     rap_obc=   (real(iteration3d  )-nest_dt_in(2))/nest_dt_in(1)    &
!          -int( (real(iteration3d  )-nest_dt_in(2))/nest_dt_in(1) )

!________________________________________________________________!
! MISE A JOUR DES FORCAGES:                                      !
! FIN.                                                           !
!________________________________________________________________!

!...............................................................................
! Mise à jour de la grille:
      call z_levels(1)
      call wetdry_mask_airseafluxes                                     !22-01-10
      call s_cpu('update grille',0) !11-02-17

! Wave/current forcing terms:
      if(iwve==1) then    ! m(0v0)m > !20-06-17
        call wave_driver  ! obligatoirement apres airseaflux_driver
        call dragcoef     ! Bottom drag coefficient
      endif               ! m(0v0)m >
!...............................................................................

      if(imodeltrc.eq.1.or.imodelbio.eq.1.or.l_sedim_s)call strada(1)  !14-06-17

!$ update location of lagrangian drifters                              !01-02-12
      if(drifter_onoff==1)call drifter_update

! Incrementer le compteur d'iteration, date, etc...
      call time_step_next ! update elapsedtime                         !21-04-11
      call chronos
      restart_file_y_or_n=1 ! Now you can write a restart file         !01-10-09

!........................................................!
! SORTIES GRAPHIQUE:                                     !
      if(flag3dbio/=1)        call date_output(idate_output)!             !01-10-09    !16-06-10 claude
      if(kpvwave==1.and.flag3dbio/=1) call date_output(3)   !             !01-10-09
!........................................................!


!............................................................!
! Restart file production 
      if(restartfileperiod>0.) &
      call dyn_restart_ask(restartfileperiod,'modulo')  !07_07-15

      if(give_chanel9.eq.1) then  !%%%%%%%%%%%%%%%%%%%%%>
                                          call dyn_restart('w')
      if(imodeltrc.eq.1.or.imodelbio.eq.1)call bio_restart('w')
      give_chanel9=0
      call io_switch('w') ! (interrupteur)                                 !21-06-10
      endif                       !%%%%%%%%%%%%%%%%%%%%%>

      if(  elapsedtime_now>=elapsedtime_end                     &      !16-04-11
       .or.kstop.eq.1) then !ksksksksksksksksksksksks>
                if(restartfileperiod>0.) then !>>>>>>>>>>>>>>>  !27-04-11
                                          call dyn_restart('w')
      if(imodeltrc.eq.1.or.imodelbio.eq.1)call bio_restart('w')
                endif                         !>>>>>>>>>>>>>>>
      endif                 !ksksksksksksksksksksksks>

!............................................................!


!***************************************************
!  Fin d'un cycle de calcul. Test de Fin.
!***************************************************


!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
      call s_cpu('outputs',0) !11-02-17
      if(elapsedtime_now<elapsedtime_end.and.kstop.ne.1)goto 9999



      call the_end(0)

!      **********   ***    ******    ****
!     *****               *******   ***
!    ****         ***    **** ***  ***
!   ********     ***    ***   **   **
!  ***          ****   ***    **  **
! **           ***    ***     *****
!****         ***    ***      ****


      return
      end subroutine model_offline
