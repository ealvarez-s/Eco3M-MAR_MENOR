      subroutine time_step
!______________________________________________________________________
! SYMPHONIE ocean model
! release 366 - last update: 14-02-23
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
#ifdef synopsis
       subroutinetitle='time_step'
       subroutinedescription= &
       'Updates dti_lp (usefull only if variable time step)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!_________________________________________________________________
! Version date      Description des modifications
!         12/08/01: amenagements sur calcul du plus grand multiple
!                   commun.
!         24/08/01: ICHANEL9 est lu dans notebook_time
!         10/09/01: si le numero de serie est > KMODEMAX on s'arrete
!                   proprement
!         12/11/01: ajout de relaxmin_ext et de relaxmin_int
!         27/02/02: debug relax_int et relax_ext
!         06/05/02: debug K5=1 si IARCHIVE=1
!         10/07/02: amenagements pour filiere "formules bulk" pour
!                   le calcul des flux à l'interface air/mer
!         23/08/02: mise à jour à partir de version 2003. Prise en
!                   compte de la marée et suppression RELAXMIN_EXT
!                   et RELAXMIN_INT pour promotion 2003
!         18/10/02: lorsque RELAX_INT ou RELAX_EXT=0 ps d'eponge
!                   respective.
!         27/01/03: la variable qui permet de selectionner le mode
!                   de sortie "graphique" du model_ s'appelle desormais
!                   IDATE_OUTPUT
!                   KOUNTGRAPH devient real
!         28/01/03: nouveaux schemas d'advection autorisés
!         25/03/03: lecture de notebook_technum supprimée
!         07/02/04: correction "bug de l'an 2000"
!         25/05/04: ajout CALL DATE_OUTPUT(4)
!                   & suppression dernier bouton de notebok_time
!         12/07/04: Si relax_es < 0 pas de correction de la moyenne
!                   de l'elevation de la surface
!         02/12/04: reset de kount_bio à kount0 et ks
!         16/02/06: RELAX_TS permet de differencier le rappel sur T & S
!                  de celui sur u et v
!         20/02/06: debug des ecritures dans fichier "messages"
!         18/04/06: fonctions compatibles avec double precision
!         08/04/07: Critere de stabilite adapté à la grille curviligne
!         17/04/07: Passage à coordonnees curvilignes (voir ajout dx_y et cie...)
!         21/01/08: si RELAX_BPC positif alors il doit être au moins egal au
!                   pas de temps barocline
!         22/01/08: renouvellement des forcages par modulo d'entier supprimé
!         23/01/08: renouvellement des forcages par modulo d'entier supprimé
!          09/04/08: test sur relax_bpc revu
!         04/12/08: possibilite de stop apres la phase d'initialisation
!         03/04/09: DTE parallelisation
!         24-05-09  augmenter la portee de DAYINDEX
! 2009.3  02-10-09  suppression dti et dte
!         05-10-09  ajout d'un "ifdef parallele"
! 2010.2  19-12-09  - time_step renommé initial_time_step
!                   - dayindex_size est la taille effective (nb de valeurs
!                   non nulles) du tableau day_index
!         22-12-09  prise en compte de tide_year_min et tide_year_max
! 2010.3  14-01-10  le calcul du pas de temps barotrope prend en compte
!                   un niveau moyen non nul, une vitesse de background max
!                   un facteur de reduction du critere cfl. Ces 3 nouveaux
!                   parametres sont lus dans notebook_time
! 2010.6  10-02-10  Ajout d'une verification de la logique des dates
! 2010.7  13-02-10  Lecture d'un notebook_time simplifié
! 2010.8  18-03-10  suppression obcdt
!         09-05-10  seul le proc 0 ecrit messages
! 2010.12 15-09-10  prendre en compte wave_year_min
!         20-09-10  Possibilité de calcul en simple precision
! 2010.20 14-04-11  dti_now = dti_fw
!         15-04-11  kountgraph renommé graphperiod
!         16-04-11  ajout elapsedtime_aft
! 2010.22 27-04-11  ichanel9 remplacé par restartfileperiod
! 2010.23 18-05-11  Suppression bornage sur relax_int et relax_ext
!         21-05-11  Partager elapsedtime_now avec tous les proc
! 2010.24 15-12-11  K2DFIN doit etre un nombre pair
! S26.1   04-09-12  autoriser le cas particulier iteration2d_max_now=1
!         01-02-13  ajout dte_fw
!         28-06-13  subcycling
!         01-11-13  petits amenagements pour dayindex
!         24-12-13  Le calcul du pas de temps ne tient plus compte des points en terre
!         29-12-13  Grille horizontale generalisee: la cfl se base sur
!                   le minmum entre dx et dy
!         17-02-14  stop si pas de temps nul
!         17-03-14  date rivieres
!         09-05-14  blindage fonction sqrt
!         06-09-14  ajout dti_lpbef et dtiratio
!         24-11-14  notebook_time passe au format namelist et ajout
!                   flag_sequoia
!         13-02-15  flag_dt_adjust et dti_fw dans notebook_time
!         17-02-15  departuredate permet de verifier qu'on n'a pas
!                   par erreur modifie la date de depart dans le cas d'un
!                   redemarrage d'un fichier restart
!         22-05-15  graphperiod defini dans initial_graph
!         13-06-15  message a l'ecran si dt=0
!         11-11-15  Aider A debusquer les erreurs commises sur les dates 
!                   de notebook_time
!                   subroutine time_step_dtimax projet de CFL pour l'advection
!         28-01-16  nouveau filtre 2dt-Ndt
!         31-01-16  - Procedure d'arrEt si nombre de courant > 5
!                   - call time_step_coeffilter
!         03-02-16  Reset dti_lpsub=dti_lp dti_fwsub=dti_fw 
!         09-02-16  mise A jour de la calibration du filtre FB
!                   et des noms de parametres associEs
!         14-02-16  Pas de cfl advection si modele 1DV 
!         01-04-16  modif du seuil d'arrEt sur dtimax
!         12-04-16  messages a l'ecran
!         03-05-16  coef du filtre laplacien pour time-stepping FB revus pour se 
!                   rapprocher d'un effet de filtrage equivalent avec le schema LF
!         14-07-16  Ne pas prendre en compte les points fleuves (masquEs) dans le calcul de dtimax
!         16-10-16  eviter cfl_reduce>1.9
!         22-12-16  flag_0status_option permet de choisir entre une interruption de la simulation
!                   ou une initialisation par defaut si le champ n'est pas present dans le fichier
!                   restart netcdf. Voir notebook_time.
!         01-03-17  ajout d'infos dans le fichier message
!         25-03-17  ajout flag_binary dans namelist
!         10-04-17  ajout restartfileunits
!         29-05-17  ajout iteration3d_max lu depuis notebook_time
!         23-10-17  Enlever time-splitting si iteration2d_max=0
!         10-11-17  modifs pour cas sans time-splittting
!         06-12-17  reset iteration2d_max_now
!         20-05-18  modif message
!         09-06-18  inv_dte_lp inv_dti_lp
!         08-07-18  cas NH avec dti_fw prescrit
!         17-08-18  cas flag_nh3d_uv > 1 (pour T et S NH) possible
!         23-08-18  supression dztimea dztimeb
!         19-10-18  ajout de spinup_forcing
!         07-12-18  la maille transverse d'un canal 2D ne limite pas le
!                   pas de temps
! v261    20-10-19  ajout simple_restart_file_txt 
!         23-10-19  Pas de temps variable en cours de simulation
! v262    24-10-19  suite point precedent
!                   subroutine time_step_dayindex 
! v264    27-10-19  suite du pas de temps variable en cours de simulation
! v268    24-11-19  initialiser date calendaire
! v286    17-06-20  lecture en serie de groupes de notebook_time
! v287    18-07-20  restartdir_out1,restartdir_out2,restartdir_in lus dans notebook_time
!         17-08-20  Des quantites de type integer pour les utilisations du coupleur OASIS
! v290    02-10-20  inv_dti_fw_p2=(1./dti_fw)**2 !02-10-20
! v292    22-11-20  relax_lwf=1./max(relax_lwf*86400., dti_fw) 
! v296    27-02-21  suppression d'un stop
! v303    19-07-21  modulo_biotimestep=20       
! v345    06-04-22  ajout simple_restart_biofile_txt  
! v361    31-12-22  inv_dti_fw
! v363    13-01-23  cfl_hsshmax pour prevoir les cas oU h_w au dessus du niveau 0
! v366    14-02-23  flag_w_binary,flag_r_binary  
!...............................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!......................................................................! m[°o°]m

      if(iteration2d_max_now/=0) then !>>> !23-01-17

! Time step of the internal mode dti_bef=t-(t-1) dti_now=(t+1)-t
      dti_bef=0.5*dte_lp*iteration2d_max_bef
      dti_now=0.5*dte_lp*iteration2d_max_now
      dti_fw=dti_now
! dti_lp=(t+1)-(t-1) dti_lpbef=t-(t-2)
      dti_lpbef=dti_lp          !06-09-14
      dti_lp=dti_bef+dti_now
      dtiratio=dti_lpbef/dti_lp !06-09-14

      endif                           !>>> !23-01-17

      rampe=min(max(elapsedtime_now/spinup_forcing,zero),un)    !03-04-11!19-10-18

      end subroutine time_step

!-----------------------------------------------------------------

      subroutine time_step_next
      use module_principal
      use module_parallele !#MPI
      implicit none
#ifdef synopsis
       subroutinetitle='time_step_next'
       subroutinedescription= &
       'Increments the elapsed time of the simulation ' &
       //'(elapsedtime_now etc...)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Number of iterations within an external mode sequence:
      iteration2d_max_bef=iteration2d_max_now

! Variable time step: !23-10-19
      if(variable_time_step_txt/='none') then !m°v°m> !23-10-19
! Update iteration2d_max_now:
       if(elapsedtime_now>dt_nextime) then !>>>
         if(dt_nextrec<=dt_nextrecmax)call time_step_variable ! compute iteration2d_max_now
       endif                               !>>>
      endif                                   !m°v°m> !23-10-19

! internal time step depending on the stability index:
      if(flag_dt_adjust==1) then !>>>>>>>> !13-02-15
!      x1=1.-stability_index
       x1=0.8-stability_index
       iteration2d_max_r8=min(max(                          &
       iteration2d_max_r8+max(x1,0.)**2-10.*max(-x1,0.)**2  &
        ,2.),iteration2d_upbound)
! Send iteration2d_max_r8 to the other mpi domains
       call mpi_bcast(iteration2d_max_r8,1,mpi_double_precision,0,par%comm2d,ierr)
       iteration2d_max_now=2*int(0.5*iteration2d_max_r8)
      endif                      !>>>>>>>> !13-02-15

! Time elapsed since the beginning of the run
      elapsedtime_bef=elapsedtime_now                                   !03-04-11
      elapsedtime_now=elapsedtime_aft
      elapsedtime_aft=elapsedtime_now+dti_now                           !16-04-11

! Des quantites de type integer pour les utilisations du coupleur OASIS ! !17-08-20
      dti_fw_i4=nint(dti_fw)                              !17-08-20
      elapsedtime_now_i4= elapsedtime_now_i4 +dti_fw_i4   !17-08-20
      elapsedtime_now_i8= elapsedtime_now_i8 +dti_fw_i4   !17-08-20
! Temps ecoule en version real*16 en prevision d'une possibilite de precision accrue sur cette information
      elapsedtime_now_r16=elapsedtime_now_r16+dti_fw      !17-08-20

#ifdef parallele
      if(subcycle_synchro==1)                                 &
      call mpi_bcast(elapsedtime_now,1,mpi_double_precision   &
                                    ,0,par%comm2d,ierr)                 !#mpi !21-05-11
#endif

      end subroutine time_step_next

!-----------------------------------------------------------------

      subroutine time_step_initial
      use module_principal ; use module_parallele
      use module_subcycle  ; use module_q
      use module_cpl_oasis
      implicit none
      double precision dteglob_
#ifdef synopsis
       subroutinetitle='time_step_initial'
       subroutinedescription=                                &
          'Computes time steps of external and'              &
       //' internal mode equations. Departure and end times.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Reads notebook_time:
      call time_step_read_notebook_time !24-11-14
      call time_step_coeffilter         !31-01-16
      call time_step_dayindex           !24-10-19

      dt_nextime=-999. 
      if(variable_time_step_txt/='none') then !m°v°m> !23-10-19
       if(par%rank==0) then !rank0>
       nc=0
       open(unit=9,file='tmp/messages',position='append')
       write(9,'(a)')'-----------------------------'
       write(9,'(a)')'subroutine time_step_initial'
       write(9,'(a,a)')'Liste des pas de temps dans ',trim(variable_time_step_txt)
       open(unit=4,file='tmp/variable_time_step.binrec' &
                        ,access='direct'                &
                        ,recl=16                        &
                        ,form='unformatted'             &
                        ,status='new'                   &
                        ,iostat=k0)
        open(unit=3,file=variable_time_step_txt) 
 235      read(3,*,end=236)x1                ! dti_fw
          read(3,*,end=236)i1,i2,i3,i4,i5,i6 ! y,m,d,h,m,s
          call datetokount(i1,i2,i3,i4,i5,i6)  ! donne elapsedtime_out
          if(elapsedtime_out>0.) then !--->
           nc=nc+1
           write(4,rec=nc)x1,elapsedtime_out
           dt_nextrecmax=nc
           if(nc==1) then !ooo>
                          dti_fw=x1 
                          dt_nextime=elapsedtime_out 
                          dt_nextrec=1 
           endif          !ooo>
           write(9,*)'x1=',x1,' jusque t=',elapsedtime_out
           write(9,*)'soit',i1,i2,i3,i4,i5,i6
          endif                       !--->
      
          goto 235
 236    close(3)
        close(4)
        close(9)
       endif                !rank0>
! rank0 envoie dt_nextrec,x1,x2 aux autres rank !24-10-19
       call mpi_bcast(dt_nextrec   ,1,mpi_integer         ,0,par%comm2d,ierr)
       call mpi_bcast(dt_nextrecmax,1,mpi_integer         ,0,par%comm2d,ierr)
       call mpi_bcast(dti_fw       ,1,mpi_double_precision,0,par%comm2d,ierr)
       call mpi_bcast(dt_nextime   ,1,mpi_double_precision,0,par%comm2d,ierr)
      endif                                   !m°v°m> !23-10-19

!$______________________________________________________________________________________
!$ computes the time step of the external mode according to Blumberg et Mellor (1987):
!$ dte_lp=2*dt where dt=1/(srqt( 1/dx**2 + 1/dy**2 )/( 2*sqrt(g*H)+Umax ) and H=h+SSHmax
      dte_lp=1.d10
      x2=sqrt(2.) !25-12-13
      do j=1,jmax
      do i=1,imax

       if(mask_t(i,j,kmaxp1)==1) then !2424242424>                   !24-12-13

! Pour ne pas etre restreint par la maille "transverse" des canaux (le   !07-12-18
! modele etant localement 2D, il n'y a pas de raisons pour que la maille
! perpendiculaire au mouvement limite le pas de temps), on ne considere
! pas dx si mask(i+1)=0 et mask(i-1) et pas dy si mask(j+1)=0 et mask(j-1)=0
         x0=1.d10
         if(mask_t(i+1,j  ,kmax)==1)x0=min(x0,dx_t(i,j)) !07-12-18
         if(mask_t(i-1,j  ,kmax)==1)x0=min(x0,dx_t(i,j))
         if(mask_t(i  ,j+1,kmax)==1)x0=min(x0,dy_t(i,j))
         if(mask_t(i  ,j-1,kmax)==1)x0=min(x0,dy_t(i,j))

!        dte_lp=min(dte_lp,                                &
!          x2*x0                                           & !29-12-13!07-12-18
!          /(2.*sqrt( max( grav*(h_w(i,j)+cfl_sshmax)      &
!                         ,small1 ) ) + cfl_umax ) )                 !09-05-14


! note: cfl_hsshmax pour prevoir les cas oU h_w au dessus du niveau 0 (cas des rivieres tres en amont de la cote) !13-01-23
         dte_lp=min(dte_lp,                                &
           x2*x0                                           & !29-12-13!07-12-18
           /(2.*sqrt( &
                      grav*max(h_w(i,j)+cfl_sshmax,cfl_hsshmax) & !13-01-23
                    ) + cfl_umax ) )                 !09-05-14

         if(dte_lp==0) then !>>>>>>>>
          write(6,*)'Err dte=0 point proc i j=',par%rank,i,j   &
         ,' dx,dy=',dx_t(i,j),dy_t(i,j) !13-06-15
          stop ' Stop dans subroutine time_step car dte_lp=0'
         endif              !>>>>>>>>

       endif                          !2424242424>

      enddo
      enddo
      dte_lp=cfl_reduce*dte_lp

#ifdef parallele
      call mpi_allreduce(dte_lp,dteglob_,1,mpi_double_precision,  & !03/04/09
                         mpi_min,par%comm2d ,ierr)
      call subcycle_dte(dteglob_) ! Subcycle case
#endif

! Cas reserves au modele NH:
! Assurer un nombre entier d'iterations en un temps donnE (dtmultiple):
!     if(flag_nh2d==1.or.flag_nh3d_uv>=1) then !m0v0m !06-12-17 !17-08-18
      if(flag_nh3d>=1) then !m0v0m !06-12-17 !17-08-18
       if(dtmultiple>0.) then    !pmx>
        k10=int( dtmultiple / dte_lp )+1
        dte_lp=dtmultiple/k10
        nh_period_modulo=2*k10 ! le 2* c'est pour aller avec dte_fw
! On aura la propriete nh_period_modulo*dte_fw=dtmultiple
       endif                     !pmx>
      endif                 !m0v0m               
      
      if(dte_lp==0)stop ' Stop dans subroutine time_step car dte_lp=0' !17-02-14
      dte_fw=0.5*dte_lp                       !01-02-13

! Internal time step:

      if(dti_fw<=0.) then !----------->

      if(iteration2d_max_now/=0) then !>>> !23-10-17
! Option 1: dti_fw is deduced from the ratio dt_internal/dt_external:
       dti_lp=dte_lp*iteration2d_max_now       !19-04-11
      else                            !>>>
! Si iteration2d_max_now=0 time-splitting annulE:
       dti_lp=dte_lp                       !23-10-17
      endif                           !>>>

       dti_fw=0.5*dti_lp


      else                !------------>

!     if(flag_nh2d==1.or.flag_nh3d_uv>=1) then !m0v0m !08-07-18 !17-08-18
      if(flag_nh3d>=1) then !m0v0m !08-07-18 !17-08-18
       if(dtmultiple>0.) then    !pmx>
        k10=int( dtmultiple /dti_fw )+1
        dti_fw=dtmultiple/k10
        nh_period_modulo=k10 ! le 2* c'est pour aller avec dte_fw
! On aura la propriete nh_period_modulo*dte_fw=dtmultiple
       endif                     !pmx>
      endif                 !m0v0m

! Option 2: the ratio dt_internal/dt_external is deduced from dti_fw
       dti_lp=2.*dti_fw

       if(iteration2d_max_now/=0) then !pmx> ! Cas standard avec time splitting !23-10-17
! Cas standard avec time splitting:
        iteration2d_max_now=int(dti_lp/dte_lp)+1
        if(mod(iteration2d_max_now,2)/=0)iteration2d_max_now=iteration2d_max_now+1
        dte_lp=dti_lp/iteration2d_max_now       !19-04-11
        dte_fw=0.5*dte_lp                       !01-02-13
       else                            !pmx>
! Cas particulier sans time-splitting:
        dte_fw=dti_fw ; dte_lp=dti_lp
       endif                           !pmx>

      endif               !------------>

      dti_lpbef=dti_lp                        !06-09-14
      dtiratio=dti_lpbef/dti_lp               !06-09-14
      dti_now=dti_fw                          !14-04-11
      dti_lpsub=dti_lp                        !03-02-16
      dti_fwsub=dti_fw                        !03-02-16
      if(dti_lp==0)stop ' Stop dans subroutine time_step car dti_lp=0' !17-02-14
      if(i2dh.eq.0)dti_fw =0.5*dte_lp
      if(i2dh.eq.0)dti_lp=     dte_lp

      inv_dte_lp=1./dte_lp !09-06-18
      inv_dti_lp=1./dti_lp !09-06-18
      inv_dti_fw=1./dti_fw !31-12-22

      inv_dti_fw_p2=(1./dti_fw)**2 !02-10-20


      iteration2d_max_r8=real(iteration2d_max_now)
      iteration2d_upbound=iteration2d_max_r8
      iteration2d_max_bef=iteration2d_max_now

       if(mod(iteration2d_max_now,2)==1) then                          !15-12-11
       if(iteration2d_max_now/=1) then              !04-09-12
       write(*,*)'Choose an even value for K2DFIN in:'
       write(*,'(A)')trim(nomfichier(1))
       stop 'STOP dans time_step.F90'
       endif
       endif

      elapsedtime_now=0.
      elapsedtime_bef=elapsedtime_now-dti_now
      elapsedtime_aft=elapsedtime_now+dti_now

! Des quantites de type integer pour les utilisations du coupleur OASIS !17-08-20
      dti_fw_i4=nint(dti_fw)                    !17-08-20
      elapsedtime_now_i8=nint(elapsedtime_now)  !17-08-20
      elapsedtime_now_i4=nint(elapsedtime_now)  !17-08-20
! Temps ecoule au format real*16 en prevision d'une necessite de precision accrue sur cette quantite
      elapsedtime_now_r16=elapsedtime_now       !17-08-20

!$______________________________________________________________________________________




!*******************************************************************************
! CALCUL DES DIFFERENTES ECHEANCES (En terme de numéro d'itération) A PARTIR
! D'UN CALENDRIER
!                                  DEBUT
!*******************************************************************************


      kount0=nint( (daysim(1) - daysim(0))*3600.*24./dti_fw )
      iteration3d= kount0
      if(iteration3d_max<1) then !---> !24-05-17
! Cas standrad:
       elapsedtime_end=(daysim(2)-daysim(0))*3600.*24.                 !16-04-11
      else                       !---> !24-05-17
! Cas ou la fin de la simulation est donnee par un nombre d'iterations !24-05-17
       elapsedtime_end=iteration3d_max*dti_fw
! entrainant la mise A jour de datesim(:,2)
       call elapsedtime2date(elapsedtime_end &
                            ,datesim(1,2)    &
                            ,datesim(2,2)    &
                            ,datesim(3,2)    &
                            ,datesim(4,2)    &
                            ,datesim(5,2)    &
                            ,datesim(6,2))
! Avertissement dans le cas du pas de temps variable car elapsedtime_end 
! ne peut pas etre donnE oar iteration3d_max*dti_fw
! commente le 27-02-21
!     if(variable_time_step_txt/='none') & !24-10-19
!     stop 'Ne pas utiliser iteration3d_max si dt variable'

      endif                      !---> !24-05-17

      nc=10                                                            !10/07/02
      do 350 k=1,nairsea
      nc=nc+1                                                          !10/07/02
  350 continue

      nc=nc+1                                                          !10/07/02
      nc=nc+1                                                          !10/07/02
      nc=nc+1                                                          !10/07/02

      do 354 k=1,nriver
      nc=nc+1                                                          !10/07/02
  354 continue

!     if(idate_output.eq.2)call date_output(4)                         !25/05/04

!*******************************************************************************
! CALCUL DES DIFFERENTES ECHEANCES (En terme de numéro d'itération) A PARTIR
! D'UN CALENDRIER
!                                  FIN
!*******************************************************************************


!.......................................................................!
! Normalisation des temps de relaxation avec le pas de temps:
      if(relax_ext.gt.small1)                                           &
!        relax_ext=1./max(relax_ext*86400.,0.5*dte_lp)                 !18/10/02
         relax_ext=1./(relax_ext*86400.)                               !18-05-11

!     if(relax_int.gt.small1)relax_int=1./max(relax_int*86400., dti_fw)!18/10/02
       if(relax_int.gt.small1)relax_int=1./(relax_int*86400.)          !18-05-11

       if(relax_ts.gt.small1) relax_ts=1./max(relax_ts *86400., dti_fw)!16/02/06

       relax_lwf=1./max(relax_lwf*86400., dti_fw) !(dti_fw pour const1=0.5) !22-11-20

      if(relax_es.ge.0.) then                                          !12/07/04
         relax_es= 0.5*dte_lp/ max( relax_es*86400. , 0.5*dte_lp)
      endif                                                            !12/07/04
      if(relax_bpc.gt.small1)relax_bpc=1./max(relax_bpc*86400.,dti_fw) !09/04/08
!.......................................................................!

      if(i2dh.eq.0) then
       stop 'Cas 2dh supprimé'              !16-04-11
!     k2dini=kount0
!     k2dfin=ks ! si mode externe seulement
      endif

! Reset des variables year_now,month_now,day_now,hour_now,minute_now,second_now qui peuvent
! etre requises en phase d'initialisation !24-11-19
     year_now  =datesim(1,1)
     month_now =datesim(2,1)
     day_now   =datesim(3,1)
     hour_now  =datesim(4,1)
     minute_now=datesim(5,1)
     second_now=datesim(6,1)

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine time_step:'
      write(3,*)
      write(3,*)'external time step dte_lp/2  =',dte_lp/2.,' seconds'

      if(i2dh.eq.1) then

      write(3,*)'internal time step dti_fw=',dti_fw ,' seconds'
      write(3,*)'iteration ration external/internal iteration2d_max',iteration2d_max_now
      if(itimets.eq.0)write(3,*)'leapfrog scheme for t & s'
      if(itimets.eq.1)write(3,*)'forward scheme for t & s'
!     if(iadvec_ts.eq.0)write(3,*)'advection scheme for t & s: beckers'
!     if(iadvec_ts.eq.1)write(3,*)'advection scheme for t & s: quick'
!     write(3,*)'internal iteration number at the end: ',ks

      endif

      write(3,*)'asselin coef. assel0=',assel0
      write(3,*)'relax conser toit (jours) ',0.5*dte_lp/relax_es/86400.
      if(relax_ext.gt.small1)                                          & !18/10/02
       write(3,*)'relax mode extern (jours) ',1./relax_ext/86400.      !20/02/06
      if(relax_int.gt.small1)                                          & !18/10/02
       write(3,*)'relax mode intern (jours) ',1./relax_int/86400.      !20/02/06
      if(relax_ts.gt.small1)                                           & !16/02/06
       write(3,*)'relax mode intern (jours) ',1./relax_ts/86400.       !20/02/06
      write(3,*)'largeur de l''eponge (en indices):',sponge_l

      iteration3d=kount0

      write(3,*)'----------------------'
      write(3,*)'Duree de la simulation' !01-03-17
      write(3,*)'  en secondes: ',elapsedtime_end
      write(3,*)'  en jours: ',elapsedtime_end/86400.
      write(3,*)'  en iterations du mode interne: ',nint(elapsedtime_end/dti_fw)
      write(3,*)'----------------------'

      if(iteration2d_max_now==0)write(3,*)'MODE3D SANS TIME-SPLITTING'
      if(dtmultiple>0.)then
       write(3,*)'dtmultiple=',dtmultiple,'secondes'
       write(3,*)'Nombre d iterations correspondanres=',nh_period_modulo
       write(3,*)'nh_period_modulo*dti_fw=',nh_period_modulo*dti_fw
      endif

      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

      end subroutine time_step_initial

!-----------------------------------------------------------------

      subroutine time_step_read_notebook_time
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='time_step_read_notebook_time'
       subroutinedescription= &
       'Reads notebook_time'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! notebook_time
      namelist/notebook_time/datesim,initial,restartfileperiod        &
       ,iteration2d_max_now,cfl_sshmax,cfl_umax,cfl_reduce,run_option &
       ,flag_sequoia,dti_fw,flag_dt_adjust,timestep_type              &
       ,flag_0status_option,restartfileunits                          & !10-04-17
       ,iteration3d_max,spinup_forcing                                & !24-05-17 !19-10-18
       ,simple_restart_file_txt                                       & !20-10-19
       ,simple_restart_biofile_txt                                    & !06-04-22
       ,variable_time_step_txt                                        & !23-10-19
       ,restartdir_out1,restartdir_out2,restartdir_in                 & !18-07-20
       ,modulo_biotimestep                                            & !19-07-21                                    
       ,cfl_hsshmax                                                   & !13-01-23 
       ,flag_w_binary,flag_r_binary                                     !14-02-23

        simple_restart_file_txt='none' !20-10-19
        spinup_forcing=86400.                                           !19-10-18
        dti_fw=-999. ; flag_dt_adjust=0 ; timestep_type=1 ; cfl_reduce=1.8
        flag_0status_option=0 ; iteration3d_max=-999 !22-12-16
        restartfileunits='d'  !10-04-17
        iteration2d_max_now=-999 !06-12-17
        restartdir_out1='restart_output/' !18-07-20
        restartdir_out2='restart_outbis/' !18-07-20
        restartdir_in='restart_input/'    !18-07-20
        modulo_biotimestep=20             !19-07-21


!.........
! Lecture en serie de groupes !17-06-20
       kread1=0
       do kreadgroup=1,kreadgroupmax
        kread2=nint(real(kreadgroup)/real(kreadgroupmax)*real(nbdom))-1
         if(par%rank>=kread1.and.par%rank<=kread2) then !pmx>
           open(100,file=nomfichier(1)) ! 'notebook_time
           read(100,nml=notebook_time)
           close(100)
         endif                                          !pmx>
         call mpi_barrier(par%comm2d,k_out)
         kread1=kread2+1
       enddo ! kreadgroup
       if(par%rank==0)write(6,*)'Lecture notebook_time OK!'
!.........


        if(cfl_reduce>1.901) & !16-10-16
        stop 'notebook_time  warning: avoid cfl_reduce>1.9'
        if(timestep_type==0)stop 'Err 685 timestep_type = leapfrog'

        write(dynrestartfilename,'(a,i0)') &
        trim(restartdir_out1)//'chanel9_',par%rank !18-07-20

! Afin d'eviter les fautes de parametrage de notebook_time dans le cas
! d'un redemarrage d'un fichier restart qui impose que seule la date de
! fin soit modifiee et surtout pas la date de depart, on archive cette derniere
! dans un tableau special qui sera archive dans le fichier reistart (alors que le
! tableau datesim ne l'est pas, ce qui pose pb). Au redemarrage il sera verifie
! que departuredate est bien identique a datesim(:,1)...
       departuredate(1:6)=datesim(1:6,1) !17-02-15

! Aider A la detection d'erreurs de date: !11-12-15
      do j=1,2
      k=1 ! k=1=OK k=0=Err
! annee 
      if(datesim(1,j)<=0)k=0 
! mois
      if(datesim(2,j)<=0)k=0
      if(datesim(2,j)>12)k=0
! jours
      if(datesim(3,j)<=0)k=0
      if(datesim(3,j)>31)k=0
! heures
      if(datesim(4,j)<0) k=0
      if(datesim(4,j)>23)k=0
! minutes
      if(datesim(5,j)<0) k=0
      if(datesim(5,j)>59) k=0
! secondes
      if(datesim(6,j)<0) k=0
      if(datesim(6,j)>59) k=0

      if(k==0) then
       write(6,*)'notebook_time Err on datesim(:,',j,') values=',datesim(:,j)
       stop 'STOP time_step_read_notebook_time'
      endif

      enddo

      end subroutine time_step_read_notebook_time
!-----------------------------------------------------------------
      subroutine time_step_dtimax !11-12-15
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='time_step_dtimax'
       subroutinedescription= &
       'Computes advection CFL'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      stop 'Je ne devrais pas pas passer par time_step_dtimax'

      if(flag_1dv==1)return ! Pas de cfl advection si modele 1DV !14-02-16

      x1=99999.        ! First guess dti_lpsub 
      x0=2.          & !
        *(1.-assel0)   !  Abaisser dti_lpmax de l'effet du filtre temporel 

      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       x2=-(x0*dxdy_t(i,j)*dz_t(i,j,k,0))       &
       /( &  !ooooo>
         -small2  & ! cste artificielle pour eviter singularite division par 0

       +(-veldydz_u(i+1,j  ,k,1)-abs(veldydz_u(i+1,j  ,k,1))  &
         +veldydz_u(i  ,j  ,k,1)-abs(veldydz_u(i  ,j  ,k,1))  &
         -veldxdz_v(i  ,j+1,k,1)-abs(veldxdz_v(i  ,j+1,k,1))  &
         +veldxdz_v(i  ,j  ,k,1)-abs(veldxdz_v(i  ,j  ,k,1))  &
        )*wetmask_t(i,j)   &
            *mask_t(i,j,k) & ! pour eviter prise en compte des points fleuves !14-07-16
        )    !ooooo>

      if(ceiling(dti_lp/x2)>=5) then !ooo>
! Si le nbre de sous-iteration advective donnE par ceiling(dti_lp/x2) devient enorme 
! c'est que la detection des zones sEches a echouE. On corrige wetmask_t
! en consequence:
       if(hz_w(i,j,2)<2.*wetdry_cst3) then !pmx>
        wetmask_t(i,j)=0.
       else                                !pmx>
        x1=min(x1,x2)

!....................................
        if(ceiling(dti_lp/x2)>7) then
        write(par%rank+10,*)'----------------------------' !12-04-16
        write(par%rank+10,*)'dti_lp/x2      ',dti_lp/x2
        write(par%rank+10,*)'i j glob k     ',i+par%timax(1),j+par%tjmax(1),k
        write(par%rank+10,*)'hz_w(i,j,0)    ',hz_w(i,j,0)
        write(par%rank+10,*)'hz_w(i,j,1)    ',hz_w(i,j,1)
        write(par%rank+10,*)'hz_w(i,j,2)    ',hz_w(i,j,2)
        write(par%rank+10,*)'2*wetdry_cst3  ',2.*wetdry_cst3
        write(par%rank+10,*)'dz_t(i,j,k,0)  ',dz_t(i,j,k,0)
        write(par%rank+10,*)'dz_t(i,j,k,1)  ',dz_t(i,j,k,1)
        write(par%rank+10,*)'dz_t(i,j,k,2)  ',dz_t(i,j,k,2)
!       write(par%rank+10,*)'vel_u(i,i+1)   ',vel_u(i:i+1,j,k,1)
!       write(par%rank+10,*)'vel_v(j:j+1)   ',vel_v(i,j:j+1,k,1)
        write(par%rank+10,*)'h_u(i:i+1)     ',h_u(i:i+1,j)
        write(par%rank+10,*)'h_v(j:j+1)     ',h_v(i,j:j+1)
        write(par%rank+10,*)'mask_u(i:i+1)  ',mask_u(i:i+1,j,k)
        write(par%rank+10,*)'mask_v(j:j+1)  ',mask_v(i,j:j+1,k)
        write(par%rank+10,*)'vel_u(i)       ',veldydz_u(i,j,k,1)/dy_u(i,j)/dz_u(i,j,k,1)
        write(par%rank+10,*)'vel_u(i+1)     ',veldydz_u(i+1,j,k,1)/dy_u(i+1,j)/dz_u(i+1,j,k,1)
        write(par%rank+10,*)'vel_v(j)       ',veldxdz_v(i,j,k,1)/dx_v(i,j)/dz_v(i,j,k,1)
        write(par%rank+10,*)'vel_v(j+1)     ',veldxdz_v(i,j+1,k,1)/dx_v(i,j+1)/dz_v(i,j+1,k,1)
        write(par%rank+10,*)'velavr_u(i:i+1)',velavr_u(i:i+1,j,1)
        write(par%rank+10,*)'velavr_v(j:j+1)',velavr_v(i,j:j+1,1)
        endif
!....................................

       endif                               !pmx>
      else                           !ooo>
! Sinon (cas normal) la cfl locale x2 est prise en compte dans la
! recherche de la cfl finale:
       x1=min(x1,x2)
      endif                          !ooo>

!     if(ceiling(dti_lp/x2)==3)write(30+par%rank,*)hz_w(i,j,2)
!     if(ceiling(dti_lp/x2)==4)write(40+par%rank,*)hz_w(i,j,2)
!     if(ceiling(dti_lp/x2)==5)write(50+par%rank,*)hz_w(i,j,2)
!     if(ceiling(dti_lp/x2)==5)write(50+par%rank,*)hz_w(i,j,2)
!     if(ceiling(dti_lp/x2)>7)write(10+par%rank,*)iteration3d,i+par%timax(1),j+par%tjmax(1),hz_w(i,j,2)

      enddo ; enddo ; enddo

      call mpi_allreduce(x1,dti_lpmax,1,mpi_double_precision,mpi_min,par%comm2d,ierr)


! dti_lpsub bornE par dti_lp:
      dti_lpsub=min(dti_lpmax,dti_lp)
      dti_fwsub=min(dti_lpmax,dti_fw) ! Oui c'est bien dti_lpmax.....
! dti_lpsub fraction entiere de dti_lp:
      dti_lpsub=dti_lp/ceiling(dti_lp/dti_lpsub)
      dti_fwsub=dti_fw/ceiling(dti_fw/dti_fwsub)

      if(par%rank==0) then !>>>>
       open(unit=3,file='tmp/dti',position='append')
!      write(3,*)real(elapsedtime_now/86400.),real(dti_lpmax),real(dti_lp),real(dti_lpsub),real(dti_fwsub)
       write(3,'(5(1x,e14.7),2(1x,i3))') &
        real(elapsedtime_now/86400.),real(dti_lpmax),real(dti_lp) &
       ,real(dti_lpsub),real(dti_fwsub)                           &
       ,nint(dti_lp/dti_lpsub)                                    &
       ,nint(dti_fw/dti_fwsub)                    
       close(3)
      endif                !>>>>

      if(nint(dti_lp/dti_lpsub)>10) then !01-04-16
      if(par%rank==0)write(6,*)'Run is stopped because of a big current' &
      ,' horizontal number is detected'
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif
      stop 'STOP subroutine time_step_dtimax'
      endif


      end subroutine time_step_dtimax
!............................................................................
      subroutine time_step_coeffilter !28-01-16
      use module_principal
      implicit none

! coefficients de l'ex filtre d'Asselin. Le coef est toujours utilisE pour
! son equivalence avec le filtre d'Asselin bien que le nouveau filtre soit
! un filtre de type laplacien:
      assel1=1.-assel0
      assel2=0.5*assel0
      assel3=1.-assel2

      if(timestep_type==timestep_leapfrog) then !>>>
! Cas du schema Leap-Frig: le coef est donnE dans notebook_visco:
       tfilterfb=0.
       tfilterlf=assel2
      endif                                     !>>>

      if(timestep_type==timestep_forwbckw) then !>>>
! Cas du schema FB le coef est donnE dans notebook_visco:
!      tfilterfb=0.5*assel2 ! moyenne demi-somme du filtre Laplacien !09-02-16
!      tfilterlf=0.01                                                !09-02-16
       tfilterfb=assel2/3. ! moyenne demi-somme du filtre Laplacien  !03-05-16
       tfilterlf=tfilterfb                                           !03-05-16
      endif                                     !>>>

#ifdef bidon
!....
! FD time filter parameters: !05-11-10
! NOTE: on prend assel2 qui est 2 fois plus petit que assel0 car le
! le filtre FD est 2 fois plus efficace que le filtre Laplacien. Cela
! permet de laisser le notebook_visco inchangé, que l'on utilise le
! filtre Laplacien ou le filtre FD, tout en gardant le même pouvoir
! filtrant du mode numérique
      tfalpha=0.545
!     tfalpha=0.5
      tfb0=1.+tfalpha*assel2/2. ! Attention tfb0 n'est pas concernee par la 1/2 somme des h
! external mode:
      tfa0=     -tfalpha*assel2/2./2. ! La 2eme division par zero est pour la 1/2 somme des h
      tfa1= (1.+tfalpha)*assel2/2./2.
      tfa2=-(2.-tfalpha)*assel2/2./2.
      tfa3= (1.-tfalpha)*assel2/2./2.
! internal mode, T & S:
      tfc0=     -tfalpha*assel2/2.
      tfc1= (1.+tfalpha)*assel2/2.
      tfc2=-(2.-tfalpha)*assel2/2.
      tfc3= (1.-tfalpha)*assel2/2.
! dz:
      tfd0= (1.+tfalpha*2.)*assel2/2.
      tfd1=           1.-3.*assel2/2.
      tfd2= (3.-tfalpha*2.)*assel2/2. !25-05-11
      tfd3= (   tfalpha-1.)*assel2/2.
      tfd4=1./tfb0
!....
! Laplacian filter:
      flag_lptimefilter=1                         !28-01-11
      if(flag_lptimefilter==1)then!lplplplplp>
      tfb0=1. ; tfd4=1.
      tfa1=assel2/2. ; tfa2=-tfa1 ; tfa0=0. ; tfa3=0.
      tfc1=assel2    ; tfc2=-tfc1 ; tfc0=0. ; tfc3=0.
      tfd0=assel2    ; tfd1=(1.-2.*assel2)  ; tfd2=assel2 ; tfd3=0.
      endif                       !lplplplplp>

! Filtre d'Asselin:                          !28-01-11
      flag_asselin=0                         !28-01-11
      if(flag_asselin==1)then!aaaaaaaa>
      tfb0=1.
      tfa1=0. ; tfa2=0. ; tfa0=0. ; tfa3=0.
      tfc1=0. ; tfc2=0. ; tfc0=0. ; tfc3=0.
      endif                  !aaaaaaaa>
#endif

      end subroutine time_step_coeffilter
!................................................................
      subroutine time_step_variable
      use module_principal ; use module_parallele 
      implicit none

       if(par%rank==0) then !rank0>
        open(unit=4,file='tmp/variable_time_step.binrec' &
                         ,access='direct'                &
                         ,recl=16                        &
                         ,form='unformatted'             &
                         ,iostat=k0)
 993    read(4,rec=dt_nextrec)x1,x2 !dti_fw,elapsedtime
        dt_nextrec=dt_nextrec+1
        if(elapsedtime_now>x2                      &
           .and.dt_nextrec<=dt_nextrecmax)goto 993  !27-10-19
        close(4)
       endif                !rank0>
! rank0 envoie dt_nextrec,x1,x2 aux autres rank
       call mpi_bcast(dt_nextrec,1,mpi_integer         ,0,par%comm2d,ierr)
       call mpi_bcast(x1        ,1,mpi_double_precision,0,par%comm2d,ierr)
       call mpi_bcast(x2        ,1,mpi_double_precision,0,par%comm2d,ierr)

       dt_nextime=x2
! x1 est dti_fw convertit en iteration2d:
       iteration2d_max_now=max(nint(2*x1/dte_lp),2)
       if(mod(iteration2d_max_now,2)/=0) then
        if(iteration2d_max_now*0.5*dte_lp<x1) then
          iteration2d_max_now=iteration2d_max_now+1
        else
          iteration2d_max_now=iteration2d_max_now-1
        endif
       endif

       if(par%rank==0) then !rank0>
         open(unit=3,file='tmp/messages',position='append')
         write(3,*) '------- new internal time step ------------'
         write(3,*) year_now,month_now,day_now,hour_now,minute_now,second_now
         write(3,*) 'Le pas de temps demandE est',x1
         write(3,*) 'Le nbre      d''iteration 2D le plus proche est:',nint(2*x1/dte_lp)
         write(3,*) 'Le nbre pair d''iteration 2D le plus proche est:',iteration2d_max_now
         write(3,*) 'Soit dti_fw=',iteration2d_max_now*0.5*dte_lp
         write(3,*) 'Prochain record',dt_nextrec
         write(3,*) 'temps actuel           ',elapsedtime_now
         write(3,*) 'Temps de remise A jour:',dt_nextime
         close(3)
       endif                !rank0>
      
      end subroutine time_step_variable
!................................................................
      subroutine time_step_dayindex !24-10-19
      use module_principal
      implicit none

! other resets & some warnings:
      iteration3d=0
      kount0=0
      do k1=1,6
       datesim(k1,0)=datesim(k1,1)  !13-02-10
      enddo
      do k1=1,6
      do k2=3,10
       datesim(k1,k2)=datesim(k1,1) !13-02-10
      enddo
      enddo

! REGLE DES ANNEES BISSEXTILES
! Une année est bissextile quand:
! Règle de base: elle est divisible par 4
! Exception1 à la règle de base: elle n'est pas bissextile
! si elle est divisible par 100.
! Exception2 à l'exception1: elle est tout de même bissextile
! si elle divisible par 400.

! JPM (jours par mois) est un tableau contenant le nb de jours
! pour chaque mois de l'année
      jpm(1) =31
      jpm(2) =28
      jpm(3) =31
      jpm(4) =30
      jpm(5) =31
      jpm(6) =30
      jpm(7) =31
      jpm(8) =31
      jpm(9) =30
      jpm(10)=31
      jpm(11)=30
      jpm(12)=31

!-----------------------------------------------------------------
! ici on passe en revue les années couvertes par la simulation.
! On en tire une valeur max et un valeur min qui servent à
! construire le tableau DAYINDEX.
!     k1=+999999
      k2=-999999

! En principe le parametre numbyears est de l'ordre de 50 ans
! le tableau dayindex couvrira le nbre d'annees correspondantes.
! Compte tenu de ce que les fichiers de forçages netcdf peuvent
! etre bien plus anciens que la date de depart du run, on applique
! un 'handicap' de quelques annees pour que le tableau dayindex
! puisse etre utilise pour ces fichiers
      k1=datesim(1,1)-5 !01-12-13

      do k=0,2                                                          !13-02-10
      k1=min0(k1,datesim(1,k))
      k2=max0(k2,datesim(1,k))
      enddo

      if(iairsea.ge.1) then ! --------------->                         !10/07/02
        k1=min(k1,airsea_year_min)                                     !17-09-13
        k2=max(k2,airsea_year_max)
      endif                 ! --------------->

      do k=1,nriver
        if(realriver(k).eq.1) then
          k1=min0(k1,dateriver(1,k))
          k2=max0(k2,dateriver(1,k))
        endif
      enddo

      if(kmaxtide.ne.0) then !**************>                          !23/08/02
          k1=min0(k1,tide_year_min)                                    !22-12-09
          k2=max0(k2,tide_year_max)
      endif                  !**************>                          !23/08/02

! REGLE DES ANNEES BISSEXTILES
! Une année est bissextile quand:
! Règle de base: elle est divisible par 4
! Exception1 à la règle de base: elle n'est pas bissextile
! si elle est divisible par 100.
! Exception2 à l'exception1: elle est tout de même bissextile
! si elle divisible par 400.
!.........................................................!
! DAYINDEX est un tableau contenant, pour tous les jours  !
! des années concernées par la simulation, un nbre entier !
! à 8 chiffres surlesquels est codée la date du jour en   !
! question. Ex: 20001130 est le code du 30 Novembre 2000. !
! Plus tard dans la simulation ce tableau permet un accés !
! rapide à la date au temps t à partir du nombre de jours !
! écoulés à depuis le début de la simulation.             !
      dayindex_size=0                                     !             !19-12-09
      if((k2-k1+1).gt.numbyears) then                     !
      write(6,*)                                        & !
       'memory size available for dayindex exceeded !'    !
      write(6,*)                                        & !
       'increase numbyears in file "parameter" up to'   & !
       ,k2-k1+1                                           !
      write(6,*)'min et max des annees:',k1,k2            !
      stop ' time_step.f'                                  !
      endif                                               !
      k2=k1-1+numbyears                                   !            !24-05-09
      j0=0                                                !
      do 151 k=0,k2-k1 ! boucle sur les annees            !
         if(mod(k+k1,4).eq.0) then                        !
               jpm(2)=29                                  !
         else                                             !
               jpm(2)=28                                  !
         endif                                            !
         if(mod(k+k1,100).eq.0)jpm(2)=28 ! Exception1     !            !07/02/04
         if(mod(k+k1,400).eq.0)jpm(2)=29 ! Exception2     !            !07/02/04
         do 157 i=1,12     ! boucle sur les mois          !
         do 157 j=1,jpm(i) ! boucle sur les jours du mois !
         j0=j0+1                                          !
           dayindex(j0)=(k+k1)*10000                    & !
                        +i*100                          & !
                        +j                                !
           dayindex_size=j0                               !            !19-12-09
  157 continue                                            !
  151 continue                                            !
!.........................................................!

! Le tableau DAYINDEX va maintenant servir à calculer le numéro
! d'itération de chaque échéance figurant dans le notebook_time.
      k6=10+nairsea                                                    !10/07/02
      k7=k6+3                                                          !10/07/02
      do 297 k3=0,2

              i1=datesim(1,k3)
              i2=datesim(2,k3)
              i3=datesim(3,k3)
              i4=datesim(4,k3)
              i5=datesim(5,k3)
              i6=datesim(6,k3)
      k1=0
      k0=i1*10000+i2*100+i3
      do 180 k=1,366*numbyears
      if(k0.eq.dayindex(k))then
      k1=k
      goto 181
      endif
  180 continue
  181 continue

      if(k1==0) then !----------->                                   !10-02-10
      write(6,*)'Date tres bizarre: ',i1,i2,i3,i4,i5,i6
      write(6,*)'k3=',k3
      write(6,*)'datesim=',datesim(1:6,k3)
      stop ' Arret dans routine initial_time_step'
      endif          !----------->

!---> DAYSIM est le numéro d'index (décimal) du jour référencé dans DAYINDEX
        daysim(k3)                                                      &
       =real(k1)+( real(i4)*3600.+real(i5)*60.+real(i6) )/86400.

  297 continue

      end subroutine time_step_dayindex
