










      subroutine model_2d
!______________________________________________________________________
! SYMPHONIE ocean model
! release 296 - last update: 18-02-21
!______________________________________________________________________

      use module_principal ; use module_parallele ; use module_airseaflux
      use module_external_mode ; use module_my_outputs ; use module_wave
      use module_modeanalysis
      use module_offline !cNA !17-12-20
! THOM SEDIM
      use module_parameter_sedim, only : l_sedim_s
      implicit none
      character*60 txt_spy_loc                                    !06-06-11

!______________________________________________________________________
! Version   date   Description
!         26/03/07: Mise en service
!                   Le model_2D.F est completement revu. Subsiste son
!                   ancienne version si I2DH=0. La nouveauté est un
!                   pseudo model_ 3D (pour pouvoir acceder aux filieres
!                   de forcage habituelles) detecté comme tel si NR=3.
!                   Le model_ est en realité 2D mais l'incrementation
!                   des iterations, l'appel aux forcages, aux emboitements
!                   se font comme dans le model_ 3D (sauf qu'il n'y a que
!                   NR=3 niveaux
!         27/03/07: L'appel à dtvarobc est déplacé de sorte que le barometre
!                   inverse soit activé même si il n'y a que la meteo comme
!                   forcage externe
!         03/04/07: filiere model_in model_out deefinitivement suprimee
!                   Appel à chronos placé apres nest_inout et offline_inout
!         24/08/07: On empeche l'acces à l'ancienne filiere 2D
! 2009.2  01-09-09: externe_fast renomme externe
!         16-09-09: Synchro des proc à la fin de chaque iteration
! 2009.3  03-10-09: externe.F devient external_mode.F
!         13-10-09: securite pour empecher d'appeler hot_restart quand ce
!                   n'est pas le moment
! 2010.2  27-12-09: iobc_ogcm renommé iobc_ogcm
! 2010.8  09-05-10  seul le proc 0 ecrit interrupteur
!         19-05-10  test nr=3 remplacé par nr<=3
! 2010.10 24-06-10  hot_restart renomme dyn_restart et bio_restart
! 2010.18 26-01-11  ajout calcul date
! 2010.20 16-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 27-04-11  ichanel9 remplacé par restartfileperiod
! 2010.23 06-06-11  Appels aux routines de time_step.F90
! 2010.24 15-12-11  Ajout potentiel generateur de maree
! 2010.25 12-04-12  appel graph_out a la fin de l'initialisation
! S25.4   30-06-12  use module_airseaflux
! S26     10-02-13  appel a airseaflux_driver
!         03-03-13  module_external_mode et cie...
!         01-04-13  appel nouvelle routine elapsedtime2date
!         23-07-13  appel my_outputs_driver
!         03-09-13  appel a wave_driver
!         11-07-14  supprime nest_inout
!         11-11-14  creation reguliere chanel9
!         28-11-14  suppression lignes inutiles
!         21-10-15  Le gradient de sshstokes est maintenant pris en compte
!                   dans la routine stokesforces
!         28-01-16  dif3d2d supprimE
!         08-06-17  ajout clef oasis
! v267    16-11-19  cas iwve==2
! v292    11-11-20  particules lagrangiennes
! v294    17-12-20  !cNA ajout des sorties offline if(ioffline==1)call offline_inout(1)
! v296    18-02-21  flag_omega_cumul=0 !18-02-21
!______________________________________________________________________

      txt_spy_loc='model_2d'              !06-06-11

! VERSION BASIQUE (ANCIENNE VERSION):
      if(i2dh.eq.0)then !0000000000000000000000>
      write(6,*)'le protocole 2d passe maintenant par'
      write(6,*)'un model_ "3d" avec nr=3'
      stop 'dans model_2d'
      endif             !0000000000000000000000>


!_____________________________________________________________________
! VERSION ELABOREE POUR FORCAGES REALISTES ET IMBRICATIONS:
! DEBUT
      if(kmax<=2) then                                                   !19-05-10
!_____________________________________________________________________


!******************************************************
! Initialisations speciales simulations 2D
! Debut.
!******************************************************
       do j=1,jmax+1
       do i=1,imax+1

! Assure que DIFVEL_Y et DIFVEL_X sont bien nuls:
             vel_v(i,j,kmin_v(i,j),1)=0.
           velavr_v(i,j,1)=0.
             vel_u(i,j,kmin_u(i,j),1)=0.
           velavr_u(i,j,1)=0.

! Pas de couplage interne/externe:
        pres3d2d_u(i,j)=0.
        adve3d2d_u(i,j)=0.
!        restoring3d2d_u(i,j)=0.
        pres3d2d_v(i,j)=0.
        adve3d2d_v(i,j)=0.
!        restoring3d2d_v(i,j)=0.

       enddo
       enddo
!******************************************************
! Initialisations speciales simulations 2D
! Fin.
!******************************************************


!******************************************************
! Debut de la boucle temporelle.
!******************************************************
 9999 continue
       restart_file_y_or_n=0 ! You can not write a restart file !13-10-09
       flag_omega_cumul=0 !18-02-21

      call time_step ! Gives the internal time step       !06-06-11

      call elapsedtime2date(elapsedtime_now,                &          !01-04-13
                            year_now,month_now,day_now,     &
                            hour_now,minute_now,second_now)

!________________________________________________________________!
! Update forcing terms:
!________________________________________________________________!

! air/sea fluxes:
      call airseaflux_driver(2)            !10-02-13
!$ River water discharge:
      call river_upd(2)                                                !09/08/01

!$ Tidal harmonic components:
      if(kmaxtide>0) then !---> !26-02-13
         call update_tide ; call s_cpu('update_tide',0)
         if(flag_3dwaves_harmonics==1) then !hhh>
           call modeanalysis_harmonics(1)   !26-02-13
           call s_cpu('3Dharmonics',0)      !20-02-16
         endif                              !hhh>
      endif               !--->

!$ Outer ocean fields:                                                 !17/06/01
      if(iobc_wv.eq.1.or.iobc_ogcm.eq.1   &
                     .or.nest_onoff_in.eq.1) then !$$$$$$$$$$$$$$$$$>  !31/08/03
!$ General case:
      if(iobc_ogcm==1.and.iobc_lr==0)call read_ogcm_fields(2)
!$ Particular case of self nesting:
!     if(nest_onoff_in.eq.1.and.kount.ne.kount0)call nest_inout(2)     !31/08/03
!     if(nest_onoff_in.eq.1.and.iteration3d.ne.kount0)call nest_inout(2)     !31/08/03

      endif                                       !$$$$$$$$$$$$$$$$$>

!$ Wave/current forcing terms:
     if(iwve==1.or. & !16-11-19
        iwve==2)call wave_driver ! obligatoirement apres airseaflux_driver 

!$ Forcing terms in open boundary conditions:
      call update_obcforcingterms(1)                                    !26/03/07
      call s_cpu('forcingterms',0)
!________________________________________________________________!

! Terms of the PGF that are frozen during the external mode sequence:   !16-12-11
      do 8 j=2,jmax-1
      do 8 i=2,imax

      pres3d2d_u(i,j)=                                                &
        (tidepotential_w(i-imu,j,1)-tidepotential_w(i+ipu,j,1)        & ! Tidal  equilibrium SSH gradient
!           +sshstokes_w(i-imu,j  )    -sshstokes_w(i+ipu,j  )        &
                                                              )*grav  & ! 21-10-15
       +(          pss_w(i+ipu,j,1)          -pss_w(i-imu,j,1))/rho     ! Atmos. pressure force

    8 continue

      do 7 j=2,jmax
      do 7 i=2,imax-1

      pres3d2d_v(i,j)=                                                &
        (tidepotential_w(i,j-jmv,1)-tidepotential_w(i,j+jpv,1)        & ! Tidal  equilibrium SSH gradient
!           +sshstokes_w(i,j-jmv  )    -sshstokes_w(i,j+jpv  )        &
                                                              )*grav  & ! 21-10-15
       +(          pss_w(i,j+jpv,1)          -pss_w(i,j-jmv,1))/rho     ! Atmos. pressure force

    7 continue

! Add vortex forces to pres3d2d_u & pres3d2d_v:
      if(iwve==1)call stokesforces_vortex_2dmode            !30-05-17

      if(iteration3d==iteration3d_restart)call graph_out    !11-11-14

! calcul du mode externe:
      call external_mode_driver !03-03-13

      if(drifter_onoff==1)call model_2d_drifter !11-11-20

! SORTIES GRAPHIQUE:
                          call date_output(idate_output)
      if(kpvwave.eq.1)    call date_output(3)

      call my_outputs_driver !23-07-13

!cNA ajout des sorties offline
       if(ioffline==1)call offline_inout(1) !17-12-20

! Couplage surfex !30-05-17
! Couplage ww3


! Incrementer le compteur d'iteration, date, etc...
      call time_step_next
      call chronos

!.......................................................!
! SORTIE FICHIER RESTART                                !

      restart_file_y_or_n=1 ! You can write a restart file              !13-10-09

!     if(ichanel9.eq.2)then !--------------------------->
!                                         call dyn_restart('w') !23-06-10
!     if(imodeltrc.eq.1.or.imodelbio.eq.1)call bio_restart('w')
!     endif                 !--------------------------->

      if(restartfileperiod>0.)                                 & !11-11-14
      call dyn_restart_ask(restartfileperiod,'modulo')

      if(give_chanel9.eq.1) then  !%%%%%%%%%%%%%%%%%%%%%>
                                          call dyn_restart('w')
      if(imodeltrc.eq.1.or.imodelbio.eq.1)call bio_restart('w')
      give_chanel9=0
      call io_switch('w') ! (interrupteur)                                 !21-06-10
      endif                       !%%%%%%%%%%%%%%%%%%%%%>

      if(  elapsedtime_now>=elapsedtime_end                &  !16-04-11
       .or.kstop.eq.1) then !ksksksksksksksksksksksks>
                if(restartfileperiod>0.) then !>>>>>>>>>>>>>>>  !27-04-11
                                          call dyn_restart('w')
      if(imodeltrc.eq.1.or.imodelbio.eq.1)call bio_restart('w')
                endif                  !>>>>>>>>>>>>>>>
      endif                 !ksksksksksksksksksksksks>
!.......................................................!

!***************************************************
!  Fin d'un cycle de calcul. Test de Fin.
!***************************************************


      if(elapsedtime_now<elapsedtime_end.and.kstop.ne.1)goto 9999  !16-04-11


      call the_end(0)
      call the_end(1)                                                  !05/08/03

!      **********   ***    ******    ****
!     *****               *******   ***
!    ****         ***    **** ***  ***
!   ********     ***    ***   **   **
!  ***          ****   ***    **  **
! **           ***    ***     *****
!****         ***    ***      ****



!_________________________________________________________
! VERSION ELABOREE POUR FORCAGES REALISTES ET IMBRICATIONS:
      return
      endif
!_________________________________________________________

      return
      end subroutine model_2d


!..................................................................

      subroutine model_2d_drifter !11-11-20
      use module_principal ; use module_drifter
      implicit none


! Charger les tableaux vel_u et vel_v utilises par module_drifter:

      call z_thickness(1,1) ! donne hz_u(i,j,1), hz_v(i,j,1)
      const2=1.0/real(iteration2d_max_now+iteration2d_max_bef)

      do j=1,jmax ; do i=1,imax+1
      vel_u(i,j,kmax,1)=                                &
       const2*(fluxbar_sumt_u(i,j,0)                    &
              +fluxbar_sumt_u(i,j,1))                   &
                        /dy_u(i,j)                      &
                        /hz_u(i,j,1)                    & 
              -velbarstokes_u(i,j,1)             
      enddo ; enddo
      do j=1,jmax+1 ; do i=1,imax               
      vel_v(i,j,kmax,1)=                                &
       const2*(fluxbar_sumt_v(i,j,0)                    &
              +fluxbar_sumt_v(i,j,1))                   &
                        /dx_v(i,j)                      &
                        /hz_v(i,j,1)                    &            
              -velbarstokes_v(i,j,1)                
      enddo ; enddo

!$ update location of lagrangian drifters
      call drifter_update

      end subroutine model_2d_drifter

!..................................................................
