      subroutine scalars
!______________________________________________________________________
! SYMPHONIE ocean model 
! release 365 - last update: 26-01-23
!______________________________________________________________________
      use module_principal                                        !27/01/03
      use module_parallele                                        !10-05-14
      use module_my_outputs
      use module_webcanals                                        !18-04-19
      implicit none
#ifdef synopsis
       subroutinetitle='scalars'
       subroutinedescription= &
      'Driving routine for temperature and salinity subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!...............................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!...............................................................................
! Version date      Description des modifications
!         26/12/02: on ne rappelle pas couple_modes et omega_upd. Les vitesses
!                   horizontales et verticale forward et leapfrog sont calculees
!                   une fois pour toute au premier appel de ces 2 routines
!         12/01/03: choix de plusieurs schemas d'advection dont le schema
!                   "QUICK_CA" implémenté pour modele par Guillaume Reffray
!                   du LSEET.
!         27/01/03: passé dans  version 2003
!         03/11/03: bilan d'energie sur terme de melange vertical....
!         29/05/04: possibilité de faire CALL OBC_SCAL avant la turbulence
!                   & nouveaux arguments dans OBC_SCAL
!         08/07/04: routines non compatibles commentées
!         16/06/06: ajout possibilite de convection automatique sur suggestion
!                   de Claude.
!         17/10/06: Nouvelle advection pour T et S
!         20/12/06: schema d'advection choisi par le passage en argument de
!                   IADVEC_TS
!         21/12/06: Les conditions aux limites sur les rivieres sont calculées
!                   avant l'appel à l'advection
!         10/12/07: ajout d'une routine calculant le niveau vertical correspondant
!                   à la base de la couche de mélange de surface
!         14/01/08: Ajout controle pression de fond dans la couche eponge
!         21/01/08: des parametres ajoutés dans notebook_spongelayer permette
!                   de choisir les conditions aux limites et les eponges
!         31/01/08: eponges calculees dans routine specifique -> appel à
!                   APPLY_SPONGE_LAYER. Deux cas possibles: 1=eponge classique
!                   2=eponge à densité inchangée
!         26/02/08: Dans cas 2. Avant on modifiait en dur les tableaux TOBC et SOBC
!                   (à densite inchangee) puis on relaxait vers ces derniers.
!                   Mais cela occasionnait une interference avec obc_scal qui utilise
!                   ces mêmes tableaux. Maintenant le calcul est le même mais temobc
!                   et salobc sont laissés intacts.
!         07/04/08: ajout RELAXTYPE_TS pour distinguer la relaxation des
!                   traceur en mode classique du mode "densité inchangée"
!         13/05/08: Deplacement de la condition au limite sur TEM et SAL juste
!                   avant l'advection. Consequence: interdit les obc sur THZ(2)...
!         05/11/08: Appel à K_SURFLAYER avant l'advection pour eventuellement
!                   utiliser ksl_z dans le schema d'advection.
!         01-06-09: Appel à obc_scal: passage en arg de obctype_ts pour faire
!                   un appel avec un arg different dans initial_with_obc
! 2009.3  02-10-09: dti_lp remplace dti
!                   utilisation des nouveaux facteurs d'echelle verticale
!         10-10-09: tem_c et sal_c remplacent thz_c et shz_c
!         15-10-09: manip sur i1d
! 2010.3  15-01-10  - subroutine advection_ts renommée advection_scal
!                   - subroutine temperature renommée vertmix_tem
!                   - subroutine salinity renommée vertmix_sal
! 2010.8  06-05-10  Pas d'argument passé dans advection_scal.F90
! 2010.11 16-07-10  temobc & salobc renommés temobc & salobc
!         23-07-10  equation d'etat donne densité potentielle
! 2010.14 16-12-10  Equation d'etat appelé au cas par cas depuis calcul PGF
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 22-04-11  Argument passé dans botprescont
! 2010.23 26-05-11  Possibilite de rappel de T et S vers OGCM sur toute la grille
! 2010.24 01-09-11  routine check mpi sur T et S
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S26     17-10-13  ajout d'une moyenne zonale de tem_t (choix 4) pour le cas test
!                   du jet barocline comodo
!         26-11-13  interdire l'option checkmpi quand moyenne zonale "desordonnee"
!         09-02-14  checkmpi seulement si subcycling desactive
!         01-04-14  si appel a deep_convection ne pas passer par k_surf
!         05-04-14  modif dans test advection et dans check mpi
!         08-05-14  arguments modifies dans call obc_river(1,0,1)
!         12-05-14  call obc_scal_mpi_now
!         02-07-14  nouveaux echanges pour check mpi: call obc_scal_mpi_z0_now
!         11-07-14  rap_obc devient timeweightobc(:)
!         14-07-14  ajout d'une routine de test de schema d'advection pour la
!                   bio; test_advection_bio.
!         31-07-14  Nouveaux echanges
!         15-08-14  reorganisation du calcul des coef de la matrice de melange
!                   prevoyant la possibilite d'une advection vertical implicite
!         26-02-15  plus d'affichage dans fichiers check_mpi
!         21-05-15  suppression d'une ecriture dans un fichier fort....
!         27-11-15  check mip 1DV case
!         11-01-16  Ajout subroutine scalars_ts_minmax
!         15-01-16  Modification de l'algo de rappel du cas baroclinic_jet
!         20-01-16  baroclinic_jet: moyenne zonale calculEe sur temobc-tem
!                   et ajout du cas type_ts=40 qui est comme type_ts=4 + rappel courant
!         28-01-16  nouveau schema de time-stepping hybride LF-FB
!         27-03-16  subroutine scalars_global_conservation permet de verifier les
!                   proprietes de conservation du schema d'advection et des flux de
!                   surface
!         30-03-16  L'algo de spectral nudging est modifiE. Les tableaux
!                   sallwf et temlwf representent desormais la Basse
!                   Frequence de sal_t-salobc_t et de tem_t-temobc_t
!         02-04-16  Initialisation de temlwf et sallwf dans initial_with_obc
!         06-08-16  ajout de test d'advection pour vel_u, tkea, epsa
!         13-11-16  subroutine pour tester advection tem sal
!         31-01-17  routine test_advection_vel_u updated
!         08-04-17  s-z coordinate
!         21-04-17  s-z coordinate suite
!         17-10-17  suite du point precedent
!         18-05-18  nouvel algo rappel T,S obc
!         25-05-18  suite point precedent
!         24-08-18  if(iteration2d_max_now/=0)call couple_modes_scalars !24-08-18
!         03-09-18  call couple_modes_scalars est maintenant appelE en dehors de scalars
!         21-10-18  ajout d'un cas n° 5 pour le rappel vers T,S
! v245    29-01-19  call my_outputs_sal_sum !29-01-19
! v247    01-03-19  call vertmix_merged_levels_ts(2) !01-03-19
! v248    04-03-19  if(tem_validmax<999.)call scalars_validbounds !04-03-19
! v252    18-04-19  La mer envoie tem_t(2),sal_t(2) aux extremitEs des canaux
! v253    05-05-19  sortir de scalars_k_surflayer si convect_yn=1
!                   if(convect_yn==3)call deep_convection 
! v256    11-06-19  call my_outputs_zone1salttempflux
! v257    07-05-19  call my_outputs_zone1salttempflux('fluxbar',0) 
! v261    20-10-19  - ajout rappel (optionnel) mer Marmara
!                   - routine min max temobc salobc
! v265    01-11-19  correction signe 
! v266    02-11-19  cas relaxtype_ts==6 rappel T,S upwind basE sur vel_u et vel_v
! v275    03-03-20  ajout subroutine scalars_eps_tken_minmax
! v292    25-11-20  amelioration du nudging "selectif" (voir aussi notebook_spongelayer)
! v301    21-03-21  - schema 5: avec le calcul mpi cela n'a pas de sens d'utiliser 
!                   spo_i et spo_j qui ont ete supprimES
!                   - flag_sponge_txt=='dxdy' n'est pas compatible avec certains schema de rappel
! v302    20-05-21  adaptation a la norme gfortran
! v309    19-08-21  adaptation mpi subroutine test_advection_temsal
!         23-08-21  mise A jour pour tester le schema vertical Lax-Wendroff avec 
!                   limiteur d'oscillation et diffusion negative
!         11-10-21  Relaxation dependante des vitesses: decalage d'indices entre traceurs 
!                   et vitesses conduisant A utiliser velobc_u(iglb+1,:) et velobc_v(:,jglb+1)
! v310    02-11-21  Melange des couches fusionnEes: supprimE le 02-11-21
! v325    12-02-22  ajout d'un projet de traitement de l'instabilite convective (non activE)
! v327    19-02-22  ajout subroutine scalars_instabilite_max
! v329    23-02-22  call vertmix_merged_levels_ts(2) !01-03-19!23-02-22
! v337    18-03-22  - suite du point precedent
!                   - ecrire iteration3d dans diag T,S min,max
! v345    06-04-22  Plus de diag dans routine minmax T et S
! v357    21-10-22  *sponge_t(i,j,1) !  note: sponge_t(:,:,1)=1 si rappel 'full'
! v359    27-10-22  Limiteur d'oscillation "28" dans nudging selectif
! v361    02-01-23  mises A jours des routines test_advection
! v365    26-01-23  call obc_scal_mpi_after      !26-01-23
!...............................................................................

!     call test_advection_temsal !13-11-16
!     call test_advection_bio !14-07-14
!     call test_advection_vel_u
!     stop 'iii'

      if(itimets.eq.1)stop 'plus de schema forward pour t & s svp!!'    !26/12/02
      itime=itimets

      if(timestep_type/=timestep_forwbckw) & !24-08-18
      stop 'Err 130 timestep_type/=timestep_forwbckw'

!     if(iteration2d_max_now/=0)call couple_modes_scalars !03-09-18

! Niveau vertical correspondant à la base de la couche de mélange      !05/11/08
      if(convect_yn==0)call scalars_k_surflayer  !05-05-19

      if(flag3d.eq.1) then ! §§§§§§§§§§§§§§§§§§§>

            call obc_river(1,0,1) ! C.L. T & S t-dt, t !08-05-14

      endif             ! §§§§§§§§§§§§§§§§§§§>                         !15-10-09

! Lateral Open Boundary Conditions:
!           call obc_scal(obctype_ts)                                  !01-06-09
! mpi exchange with mpi neighbours:
            if(iteration3d==0)call obc_scal_mpi_now

#ifdef checkmpi
      if(subcycle_onoff==0)call scalar_check_mpi_conservation          !09-02-14
#endif

! computes the coefficients of the "vertical mixing" matrix
! (the former are eventually completed by advection_scal) 
         call vertmix_matrix !15-08-14

! Computes advection diffusion T,S fluxes:
         call advection_scal                                        !06-05-10

! Updates T and S fields:
         call vertmix_tem(1) ! arg=flagsolver_ 1 = computes matrix pretreatment
         call vertmix_sal(0) ! arg=flagsolver_ 0 = does not compute matrix pretreatment

! Convection automatique si CONVECT_YN=1 (voir notebook_visco)
      if(convect_yn==1)call convect                                  !16/06/06
      if(convect_yn==3)call deep_convection !05-05-19

! Melange des couches fusionnEes: supprimE le 02-11-21, retabli le 23-02-22
!        call vertmix_merged_levels_ts(2) !01-03-19!23-02-22!18-03-22

! Calculer le rappel dans la zone eponge:                              !31/01/08
      if(relax_ts>0.)call scalars_spongelayer                            !07/04/08
! BUFFER ZONE MER MARMARA
!     call scalars_spongelayer_marmara

! Appliquer les bornes des intervalles de validite:   !04-03-19
      if(tem_validmax<999.)call scalars_validbounds   !04-03-19

! La mer envoie tem_t(2),sal_t(2) aux extremitEs des canaux !18-04-19
      if(nbcanal>0) then !m°v°m>
          call webcanals_gt1_to_gt2_temsal
          call webcanals_gt2_to_gt1_temsal
      endif              !m°v°m>


!     call scalars_eps_tken_minmax
!     call scalars_ts_minmax(2)                                        !11-01-16
!     call scalars_tsobc_minmax(1) 
!     call scalars_global_conservation                                  !27-03-16
!     call my_outputs_sal_sum     !29-01-19
!     call my_outputs_tem_sum     !29-01-19
!     call my_outputs_ssh_int_sum !29-01-19

#ifdef bilants6
      call my_outputs_zone1salttempflux('k',0) ; call my_outputs_zone2salttempflux('k',0) ; call my_outputs_zone3salttempflux('k',0) ; call my_outputs_zone4salttempflux('k',0) ; call my_outputs_zone5salttempflux('k',0) ; call my_outputs_zone6salttempflux('k',0) !11-06-19
      call my_outputs_zone1salttempflux('fluxbar',0) ; call my_outputs_zone2salttempflux('fluxbar',0) ; call my_outputs_zone3salttempflux('fluxbar',0) ; call my_outputs_zone4salttempflux('fluxbar',0) ; call my_outputs_zone5salttempflux('fluxbar',0) ; call my_outputs_zone6salttempflux('fluxbar',0) !08-07-19
      if(mod(iteration3d,10)==0) then
      call my_outputs_zone1salttempflux('mpi',0) ; call my_outputs_zone2salttempflux('mpi',0) ; call my_outputs_zone3salttempflux('mpi',0) ; call my_outputs_zone4salttempflux('mpi',0) ; call my_outputs_zone5salttempflux('mpi',0) ; call my_outputs_zone6salttempflux('mpi',0) !11-06-19
      endif 
#endif

!     sum1=0.
!     sum2=0.
!     sum3=0.
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)
!      sum2=sum2+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*tem_t(i,j,k,1)
!      sum3=sum3+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*sal_t(i,j,k,1)
!     enddo
!     enddo
!     enddo
!     write(6,*)'moyenne T et S',sum2/sum1,sum3/sum1

#ifdef checkmpi
      if(flag_1dv==1)call check_mpi_1dv_ts('A') !27-11-15
#endif


! Lateral Open Boundary Conditions: ! deplacee ici le 26-01-23
       call obc_scal(obctype_ts)                                  
! mpi exchange with mpi neighbours: ! deplacee ici le 26-01-23
       call obc_scal_mpi_after      !26-01-23

      end subroutine scalars

!-----------------------------------------------------------------------------
      subroutine scalars_crushed_cells
      use module_principal
      use module_parallele
      implicit none

      stop 'scalars_crushed_cells'

! Les couches ecrasees (epaisseur nulle) prennent la valeur de la premiere couche non nulle
!      if(par%rank==0)write(6,*)'bidouille scalars_crushed_cells'
       do j=1,jmax
       do i=1,imax
        do k=1,kmin_w(i,j)-1
           tem_t(i,j,k,2)=tem_t(i,j,kmin_w(i,j),2)
           sal_t(i,j,k,2)=sal_t(i,j,kmin_w(i,j),2)
!          tem_t(i,j,k,2)=tem_t(i,j,kmerged_t(i,j),2)
!          sal_t(i,j,k,2)=sal_t(i,j,kmerged_t(i,j),2)
        enddo
       enddo
       enddo

      end subroutine scalars_crushed_cells
!-----------------------------------------------------------------------------
      subroutine scalars_k_surflayer                                            !10/12/07
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='scalars_k_surflayer'
       subroutinedescription= &
       'Finds the nearest k-level to the base of the surface mixed' &
       //' layer.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! ksl_t(i,j) est le niveau vertical correspondant à la base de la couche
! de mélange:
! https://docs.google.com/presentation/d/19f3SSwAQ4flJZG0lJAP-axdwOeSqtoRbucJojyTE_Ac/edit#slide=id.g57193dca64_0_11

      do j=1,jmax
      do i=1,imax
      if(  mask_t(i,j,kmax).eq.1)then !>>>>>>>>>>>>>>>>>

       ksl_t(i,j)=kmax
       do k=kmax-1,kmin_w(i,j),-1
          if(rhp_t(i,j,k).le.rhp_t(i,j,kmax))then !*********>
            ksl_t(i,j)=k
          else                                    !*********>
            goto 95
          endif                                   !*********>
       enddo
   95 continue

      endif                           !>>>>>>>>>>>>>>>>>
      enddo
      enddo

#ifdef bidon
! APPELER CETTE ROUTINE JUSTE AVANT VERTMIX et APRES LA TURBULENCE
! projet nouveau algo: !12-02-22
! Cherche A partir de quelle profondeur le melange depuis la surface
! produit une densite potentielle plus faible que la densite du niveau
! vertical en dessous et augmente KH tant que cette profondeur n'est pas
! trouvee
      do j=1,jmax ; do i=1,imax
      if(rhp_t(i,j,kmax)>rhp_t(i,j,kmax-1)) then !m°v°m>

       sum1=rhp_t(i,j,kmax)*dz_t(i,j,kmax,1)

       do k=kmax-1,kmin_w(i,j)+1,-1

        kh_w(i,j,k+1)=max(10.,kh_w(i,j,k+1)) ! (2)
        sum1=sum1+rhp_t(i,j,k)*dz_t(i,j,k,1) ! (1)

! la densite moyenne de la couche consideree est sum1/(depth_w(i,j,kmax+1)-depth_w(i,j,k))

! Si ce test est reussi c'est que la stratif est stable et on sort:
!       if(rhp_t(i,j,k-1)>sum1/(depth_w(i,j,kmax+1)-depth_w(i,j,k)))goto 314
! aussi equivalent A:
        if(rhp_t(i,j,k-1)*(depth_w(i,j,kmax+1)-depth_w(i,j,k))>sum1)goto 314

       enddo
  314  continue
      endif                                      !m°v°m>
      enddo       ; enddo
#endif

      end subroutine scalars_k_surflayer                                            !10/12/07
!-----------------------------------------------------------------------------
      subroutine scalars_spongelayer
      use module_principal
      use module_parallele
      implicit none
!     real*4 factor_relax_,factor_lwf_
      integer :: loop_        &
                ,moduliter_=1
!     double precision,dimension(:,:),allocatable ::  iaveraged_in,iaveraged_out
#ifdef synopsis
       subroutinetitle='scalars_spongelayer'
       subroutinedescription= &
          'Apply lateral nudging layer within which the T and S'     &
       //' fields are guided by the outer forcing fields (generally' &
       //' MERCATOR)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...................................................................
! Full grid relaxation of T and S toward OGCM field - algo "Low Frequency"
      if(relaxtype_ts==2) then
!...................................................................

! https://docs.google.com/document/d/1bt_mnhZCXMuNnw3JqymtFjInqG6nNHQT0vuXX0HtZps/edit?usp=sharing
! https://docs.google.com/document/d/16WuI4w66CLc5pZbKrD2GOkwETiBmNRBD-ghwyMWMqQs/edit?usp=sharing

! Relax:
      x1=dti_lp*relax_ts        !25-11-20 
      x2=timeweightobc(trc_id) ; x0=1.-timeweightobc(trc_id) !30-03-16
! Low freq:
      const2=dti_lp*relax_lwf     !25-11-20
      const1=1.-const2

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Low frequency temperature:
       temlwf_t(i,j,k)=const1*temlwf_t(i,j,k)  &

                        +const2*( & !ooo>

          (x2*temobc_t(i,j,k,2)+x0*temobc_t(i,j,k,0))-tem_t(i,j,k,0) &

                                )   !ooo>

! Nudging process based on the difference between the Low Frequcy field and the reference field:
       tem_t(i,j,k,2)=tem_t(i,j,k,2)               &
                  +temlwf_t(i,j,k)*x1              &
                                  *sponge_t(i,j,1) & !21-10-22 !  note: sponge_t(:,:,1)=1 si rappel 'full'
      *max(0.d0,sign(1.d0,temlwf_t(i,j,k)*((x2*temobc_t(i,j,k,2)+x0*temobc_t(i,j,k,0))-tem_t(i,j,k,0)))) ! Limiteur "28"  !25-11-20
!     *0.5*(1.+tanh( temlwf_t(i,j,k)*((x2*temobc_t(i,j,k,2)+x0*temobc_t(i,j,k,0))-tem_t(i,j,k,0)) ))     ! Limiteur "29"  !27-10-22
! Note l'echelle de variance pour T est 1

!Voir limiteurs 28 et 29 dans: https://docs.google.com/document/d/16WuI4w66CLc5pZbKrD2GOkwETiBmNRBD-ghwyMWMqQs/edit?usp=sharing

! Low frequency salinity:
       sallwf_t(i,j,k)=const1*sallwf_t(i,j,k)  &

                        +const2*( & !ooo>

          (x2*salobc_t(i,j,k,2)+x0*salobc_t(i,j,k,0) )-sal_t(i,j,k,0) &

                                )   !ooo>

! Nudging process based on the difference between the Low Frequcy field and the reference field:
       sal_t(i,j,k,2)=sal_t(i,j,k,2)               &
                  +sallwf_t(i,j,k)*x1              &
                                  *sponge_t(i,j,1) & !21-10-22 !  note: sponge_t(:,:,1)=1 si rappel 'full'
      *max(0.d0,sign(1.d0,sallwf_t(i,j,k)*((x2*salobc_t(i,j,k,2)+x0*salobc_t(i,j,k,0))-sal_t(i,j,k,0))))      ! Limiteur "28"  
!     *0.5*(1.+tanh( 100.*sallwf_t(i,j,k)*((x2*salobc_t(i,j,k,2)+x0*salobc_t(i,j,k,0))-sal_t(i,j,k,0)) ))     ! Limiteur "29" !27-10-22
! Note l'echelle de variance pour S est 0.1, soit 1/0.1*1/0.1=100, d'oU le facteur 100.* 

      enddo
      enddo
      enddo

#ifdef bidon
! Ce petit programme est un mini-exemple du rappel du modele
! vers la basse frequence de l'ogcm. Attention pour simplifier
! on utilise un schema FB d'oU le 0.5*dti_lp et le 0.5*x1
      program toto
      implicit none
      double precision :: temobc,tem,temlwf=0.,pi=acos(-1.),dti_lp=500. &
      ,factor_relax_=0.6 , factor_lwf_=1.37 , relax_ts=1./(30.*86400.)  &
       ,time=0.,const1,const2,x1
       const2=dti_lp*relax_ts*factor_lwf_ 
       const1=1.-const2
       x1=dti_lp*relax_ts*factor_relax_*factor_lwf_
       temlwf=0.
       tem=10.
       do while(time<200.*86400)
        time=time+0.5*dti_lp
        temobc=20.+0.5*sin(2.*pi*time/(2.*86400.)) ! OGCM 
        temlwf=const1*temlwf+const2*(tem-temobc)   ! BF de S26
        tem=tem &                                                ! FB S26
               +0.5*dti_lp*sin(2.*pi*time/(2.*86400.))/100000. & !HF S26
               -0.5*x1*temlwf                                    !Restore OGCM BF
        write(66,'(4(1x,e14.7))')time,temobc,temlwf,tem
       enddo
      end
#endif
!...................................................................
! Full grid relaxation of T and S toward OGCM field - algo "Low Frequency"
      return
      endif
!...................................................................

!...................................................................
! Full grid relaxation of T and S toward OGCM field ! popular scheme    !26-05-11
      if(relaxtype_ts==3) then
!...................................................................

      const1= 1.-dti_lp*relax_ts
      const2=(1.-const1)*    timeweightobc(trc_id)
      const3=(1.-const1)*(1.-timeweightobc(trc_id))

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

       tem_t(i,j,k,2)=tem_t(i,j,k,2)*const1           &
                  +temobc_t(i,j,k,2)*const2           &
                  +temobc_t(i,j,k,0)*const3

       sal_t(i,j,k,2)=sal_t(i,j,k,2)*const1           &
                  +salobc_t(i,j,k,2)*const2           &
                  +salobc_t(i,j,k,0)*const3

      enddo
      enddo
      enddo

!...................................................................
! Full grid relaxation of T and S toward OGCM field ! popular scheme    !26-05-11
      return
      endif
!...................................................................

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! UPWIND basE sur velobc
      if(relaxtype_ts==1) then                                       !07/04/08
!...................................................................

      if(flag_sponge_txt=='dxdy') then
       stop 'Err 446 flag_sponge_txt==dxdy' !21-03-21
! Avec l'option 'dxdy' la notion de distance A la frontiere n'a plus de
! sens et avec le choix entre velobc_u ou velobc_v n'en a plus non plus. 
! Choisir un rappel classique A la place.
      endif

! Rappel dependant du signe et intensite du courant de reference
! https://docs.google.com/document/d/1bt_mnhZCXMuNnw3JqymtFjInqG6nNHQT0vuXX0HtZps/edit
        loop1=1 ! Frontiere i=1 !18-05-18

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dx_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_u(i  ,j,k,2)    & !11-10-21
                              +(1.-timeweightobc(vel_id))*velobc_u(i  ,j,k,0),0.)

           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_u(i  ,j,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_u(i  ,j,k,0),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=2 ! Frontiere i=imax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dx_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(    -timeweightobc(vel_id) *velobc_u(i+1,j,k,2)    & !01-11-19
                              -(1.-timeweightobc(vel_id))*velobc_u(i+1,j,k,0),0.)  !01-11-19
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(    -timeweightobc(vel_id) *velobc_u(i+1,j,k,2)    & !01-11-19
                              -(1.-timeweightobc(vel_id))*velobc_u(i+1,j,k,0),0.)  !01-11-19

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=3 ! Frontiere j=1

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dy_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_v(i,j  ,k,2)    &!11-10-21
                              +(1.-timeweightobc(vel_id))*velobc_v(i,j  ,k,0),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_v(i,j  ,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_v(i,j  ,k,0),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=4 ! Frontiere j=jmax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dy_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(    -timeweightobc(vel_id) *velobc_v(i,j+1,k,2)    &!11-10-21
                              -(1.-timeweightobc(vel_id))*velobc_v(i,j+1,k,0),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(    -timeweightobc(vel_id) *velobc_v(i,j+1,k,2)    &
                              -(1.-timeweightobc(vel_id))*velobc_v(i,j+1,k,0),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j


!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! UPWIND basE sur velobc
      return 
      endif ! cas relaxtype_ts=1
!...................................................................

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! UPWIND basE sur vel_u & vel_v (et non pas velobc comme cas relaxtype_ts=1)
      if(relaxtype_ts==6) then                                       !02-11-19
!...................................................................

      if(flag_sponge_txt=='dxdy') then
       stop 'Err 547 flag_sponge_txt==dxdy' !21-03-21
! Avec l'option 'dxdy' la notion de distance A la frontiere n'a plus de
! sens et avec le choix entre vel_u ou vel_v n'en a plus non plus. Choisir un rappel classique A la place 
      endif

! Rappel dependant du signe et intensite du courant de reference
! https://docs.google.com/document/d/1bt_mnhZCXMuNnw3JqymtFjInqG6nNHQT0vuXX0HtZps/edit
        loop1=1 ! Frontiere i=1 !18-05-18

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dx_t(i,j)*sponge_t(i,j,1)*0.5 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max( vel_u(i  ,j,k,2)    & !25-05-18!11-10-21
                              +vel_u(i  ,j,k,1),0.)

           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max( vel_u(i  ,j,k,2)    &
                              +vel_u(i  ,j,k,1),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=2 ! Frontiere i=imax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dx_t(i,j)*sponge_t(i,j,1)*0.5 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(-vel_u(i+1,j,k,2)    & !01-11-19
                              -vel_u(i+1,j,k,1),0.)  !01-11-19
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(-vel_u(i+1,j,k,2)    & !01-11-19
                              -vel_u(i+1,j,k,1),0.)  !01-11-19

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=3 ! Frontiere j=1

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dy_t(i,j)*sponge_t(i,j,1)*0.5 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max( vel_v(i,j  ,k,2)    & !11-10-21
                              +vel_v(i,j  ,k,1),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max( vel_v(i,j  ,k,2)    &
                              +vel_v(i,j  ,k,1),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=4 ! Frontiere j=jmax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dy_t(i,j)*sponge_t(i,j,1)*0.5 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(-vel_v(i,j+1,k,2)    & !11-10-21
                              -vel_v(i,j+1,k,1),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(-vel_v(i,j+1,k,2)    &
                              -vel_v(i,j+1,k,1),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! UPWIND basE sur vel_u & vel_v (et non pas velobc comme cas relaxtype_ts=1)
      return
      endif ! cas relaxtype_ts=6
!...................................................................

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! at unchanged potential density (provided that a linear EOS is used)
      if(relaxtype_ts==0) then                                       !07/04/08
!...................................................................

      const1=-rho*alp_t
      const2= rho*alp_s
      const4=dti_lp*relax_ts        ! 22-04-11
      do loop1=1,4                          ! boucle sur 4 frontières
       do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
       do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i
       if(  mask_t(i,j,kmax+1).eq.1) then !>>>>>>>>>>>>>>>>

       x5=const4*sponge_t(i,j,1)           !22-04-11
       x6=1.-x5                            !22-04-11

       do k=kmin_w(i,j),kmax                ! debut boucle k

        x1=abs(temobc_t(i,j,k,2)-temobc_t(i,j,k,0))+small1
        x2=abs(salobc_t(i,j,k,2)-salobc_t(i,j,k,0))+small1

        x3=(const1*(tem_t(i,j,k,2)-temobc_t(i,j,k,1))            &
           +const2*(sal_t(i,j,k,2)-salobc_t(i,j,k,1)))           &
         /(-const1*x1+const2*x2)

        tem_t(i,j,k,2)=tem_t(i,j,k,2)*x6+x5*(temobc_t(i,j,k,1)-x1*x3)  !26/02/08
        sal_t(i,j,k,2)=sal_t(i,j,k,2)*x6+x5*(salobc_t(i,j,k,1)+x2*x3)  !26/02/08

       enddo                                ! fin boucle k
       endif                              !>>>>>>>>>>>>>>>>
       enddo                                ! fin boucle i
       enddo                                ! fin boucle j
      enddo                                ! fin boucle loop1

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! at unchanged potential density (provided that a linear EOS is used)
      return
      endif
!...................................................................

!...................................................................
! Rappel moyenne zonale cas testjet barocline comodo !17-10-13!15-01-16
      if(relaxtype_ts==4.or.relaxtype_ts==40) then !20-01-16 Le cas 40 correspond a T + courant
      if(mod(iteration3d,moduliter_)/=0) return ! moduliter_ initialisE dans sa declaration
!...................................................................
      if(.not.allocated(iaveraged_in ))allocate(iaveraged_in (jglb,kmax))
      if(.not.allocated(iaveraged_out))allocate(iaveraged_out(jglb,kmax))

!     stop 'scalars'

!     write(6,*)'bidouille tem_t'
!     do k=1,kmax
!      do j=1,jmax
!       j1=j+par%tjmax(1)
!       do i=1,imax
!        tem_t(i,j,k,1)=real(j1)
!       enddo
!      enddo
!     enddo


      x2=timeweightobc(trc_id) ; x0=1.-x2
      iaveraged_in=0.
      do k=1,kmax
       do j=1,jmax
        j1=j+par%tjmax(1)
        do i=1,imax
        iaveraged_in(j1,k)=iaveraged_in(j1,k)+mask_i_w(i)*mask_j_w(j)*(&!ooo>!20-01-16
                           x2*temobc_t(i,j,k,2)                        &
                          +x0*temobc_t(i,j,k,0)                        &
                                -tem_t(i,j,k,0)                        &
                                                                      ) !ooo>
        enddo
       enddo
      enddo
      call mpi_allreduce( iaveraged_in         &
                         ,iaveraged_out        &
                         ,jglb*kmax            &
                         ,mpi_double_precision & 
                         ,mpi_sum              &
                         ,par%comm2d ,ierr)


      x0=dti_lp*moduliter_*relax_ts/real(iglb-2)
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       tem_t(i,j,k,2)=tem_t(i,j,k,2)+iaveraged_out(j+par%tjmax(1),k)*x0 
      enddo ; enddo ; enddo

!...................................................................
! Rappel moyenne zonale cas testjet barocline comodo !17-10-13
      return
      endif
!...................................................................

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! BASIC
      if(relaxtype_ts==5) then         !07/04/08 !21-10-18
!...................................................................

! Avec la parallelisation, cela n'a pas de sens d'utiliser spo_i et ! !21-03-21
! spo_j (sauf dans le cas du schema oU le rappel depend de velobc_u ou velobc_v)
        const4=dti_lp*relax_ts
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
         tem_t(i,j,k,2)=(tem_t(i,j,k,2)+const4*sponge_t(i,j,1)*temobc_t(i,j,k,1)) &
                       /(            1.+const4*sponge_t(i,j,1)                  )
         sal_t(i,j,k,2)=(sal_t(i,j,k,2)+const4*sponge_t(i,j,1)*salobc_t(i,j,k,1)) &
                       /(            1.+const4*sponge_t(i,j,1)                  )
        enddo ; enddo ; enddo

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! BASIC
      return
      endif
!...................................................................
      write(*,*)'Irrelevant value for relaxtype_ts given' !26-05-11
      write(*,*)'in notebook_spongelayer!'
      stop 'STOP in subroutine scalars_spongelayer'

      end subroutine scalars_spongelayer

!-----------------------------------------------------------------------------

      subroutine scalars_validbounds !04-03-19
      use module_principal ; use module_parallele
      implicit none

!     tem_t(:,:,:,2)=min(max(tem_t(:,:,:,2),-3.),40.) !04-03-19
!     sal_t(:,:,:,2)=min(max(sal_t(:,:,:,2), 0.),40.) !04-03-19
      tem_t(:,:,:,2)=min(max(tem_t(:,:,:,2),tem_validmin),tem_validmax) !04-03-19
      sal_t(:,:,:,2)=min(max(sal_t(:,:,:,2),sal_validmin),sal_validmax) !04-03-19

      end subroutine scalars_validbounds

!-----------------------------------------------------------------------------

      subroutine scalar_check_mpi_conservation
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
#ifdef synopsis
       subroutinetitle='scalar_check_mpi_conservation'
       subroutinedescription= &
       'Checks the mpi continuity of the scalar fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour verifier la conservation de la parallelisation:
      do k=1,kmax
        j1=1 ; j2=jmax
        do i=1,imax
         anyv3d(i,j1,k,1)=tem_t(i,j1,k,1)
         anyv3d(i,j2,k,1)=tem_t(i,j2,k,1)
         anyv3d(i,j1,k,2)=sal_t(i,j1,k,1)
         anyv3d(i,j2,k,2)=sal_t(i,j2,k,1)
        enddo
        i1=1 ; i2=imax
        do j=1,jmax
         anyv3d(i1,j,k,1)=tem_t(i1,j,k,1)
         anyv3d(i2,j,k,1)=tem_t(i2,j,k,1)
         anyv3d(i1,j,k,2)=sal_t(i1,j,k,1)
         anyv3d(i2,j,k,2)=sal_t(i2,j,k,1)
        enddo
      enddo

      call obc_scal_checkmpi !02-07-14

      ksecu=0
      do k=1,kmax

      j1=1 ; j2=jmax
      do i=1,imax

      if(mask_t(i,j1,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i,j1,k,1)/=tem_t(i,j1,k,1))then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i,j1,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i+par%timax(1),j1+par%tjmax(1),k
        write(10+par%rank,*)'tem_t ',anyv3d(i,j1,k,1),tem_t(i,j1,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i,j1,k,1)-tem_t(i,j1,k,1)
        write(10+par%rank,*)'temobc',j1,j1+1,temobc_t(i,j1:j1+1,k,1)
        write(10+par%rank,*)'mask_t',j1,j1+1,  mask_t(i,j1:j1+1,k)
        write(10+par%rank,*)'veldydz_u',veldydz_u(i:i+1,j1,k,1)
        write(10+par%rank,*)'veldydz_v',veldxdz_v(i,j1:j1+1,k,1)
        write(10+par%rank,*)'velobc_u ',velobc_u(i:i+1,j1,k,1)
        write(10+par%rank,*)'velobc_v ',velobc_v(i,j1:j1+1,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      if(mask_t(i,j2,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i,j2,k,1)/=tem_t(i,j2,k,1))then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i,j2,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(nord)
        write(10+par%rank,*)'coordonnées ',i,2
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i+par%timax(1),j2+par%tjmax(1),k
        write(10+par%rank,*)'tem_t ',anyv3d(i,j2,k,1),tem_t(i,j2,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i,j2,k,1)-tem_t(i,j2,k,1)
        write(10+par%rank,*)'temobc',j2-1,j2,temobc_t(i,j2-1:j2,k,1)
        write(10+par%rank,*)'mask_t',j2-1,j2,  mask_t(i,j2-1:j2,k)
        write(10+par%rank,*)'veldydz_u',veldydz_u(i:i+1,j2,k,1)
        write(10+par%rank,*)'veldxdz_v',veldxdz_v(i,j2:j2+1,k,1)
        write(10+par%rank,*)'velobc_u ',velobc_u(i:i+1,j2,k,1)
        write(10+par%rank,*)'velobc_v ',velobc_v(i,j2:j2+1,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! i

      i1=1 ; i2=imax
      do j=1,jmax

      if(mask_t(i1,j,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i1,j,k,1)/=tem_t(i1,j,k,1)) then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i1,j,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(ouest)
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i1+par%timax(1),j+par%tjmax(1),k
        write(10+par%rank,*)'tem_t ',anyv3d(i1,j,k,1),tem_t(i1,j,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i1,j,k,1),tem_t(i1,j,k,1)
        write(10+par%rank,*)'temobc',i1,i1+1,temobc_t(i1:i1+1,j,k,1)
        write(10+par%rank,*)'mask_t',i1,i1+1,  mask_t(i1:i1+1,j,k)
        write(10+par%rank,*)'veldydz_u(i,i+1)',veldydz_u(i1:i1+1,j,k,1)
        write(10+par%rank,*)'veldxdz_v(j,j+1)',veldxdz_v(i1,j:j+1,k,1)
        write(10+par%rank,*)'velobc_u(i,i+1) ',velobc_u(i1:i1+1,j,k,1)
        write(10+par%rank,*)'velobc_v(j,j+1) ',velobc_v(i1,j:j+1,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      if(mask_t(i2,j,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i2,j,k,1)/=tem_t(i2,j,k,1)) then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i2,j,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(est)
        write(10+par%rank,*)'coordonnées',2,j
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i2+par%timax(1),j+par%tjmax(1),k
        write(10+par%rank,*)'tem_t ',anyv3d(i2,j,k,1),tem_t(i2,j,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i2,j,k,1)-tem_t(i2,j,k,1)
        write(10+par%rank,*)'temobc',i2-1,i2,temobc_t(i2-1:i2,j,k,1)
        write(10+par%rank,*)'mask_t',i2-1,i2,  mask_t(i2-1:i2,j,k)
        write(10+par%rank,*)'veldydz_u(i,i+1)',veldydz_u(i2:i2+1,j,k,1)
        write(10+par%rank,*)'veldxdz_v(j,j+1)',veldxdz_v(i2,j:j+1,k,1)
        write(10+par%rank,*)'velobc_u(i,i+1) ',velobc_u(i2:i2+1,j,k,1)
        write(10+par%rank,*)'velobc_v(j,j+1) ',velobc_v(i2,j:j+1,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! j

      enddo ! kmax

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(ksecu==1) then
        write(6,'(a,i0,a)')   &
        ' mpi conservation error for tem_t Details in fort.' &
         ,10+par%rank,' available in RDIR/CONFIG directory'
        stop ' STOP in routine scalar_check_mpi_conservation!'
      endif

      ksecu=0
      do k=1,kmax

      j1=1 ; j2=jmax
      do i=1,imax

      if(mask_t(i,j1,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i,j1,k,2)/=sal_t(i,j1,k,1))then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i,j1,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*)'sal_t ',anyv3d(i,j1,k,2),sal_t(i,j1,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i,j1,k,2)-sal_t(i,j1,k,1)
        write(10+par%rank,*)'salobc local ',salobc_t(i,j1,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      if(mask_t(i,j2,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i,j2,k,2)/=sal_t(i,j2,k,1))then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i,j2,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(nord)
        write(10+par%rank,*)'coordonnées ',i,2
        write(10+par%rank,*)'sal_t ',anyv3d(i,j2,k,2),sal_t(i,j2,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i,j2,k,2)-sal_t(i,j2,k,1)
        write(10+par%rank,*)'salobc local ',salobc_t(i,j2,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! i

      i1=1 ; i2=imax
      do j=1,jmax

      if(mask_t(i1,j,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i1,j,k,2)/=sal_t(i1,j,k,1)) then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i1,j,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(ouest)
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*)'sal_t ',anyv3d(i1,j,k,2),sal_t(i1,j,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i1,j,k,2),sal_t(i1,j,k,1)
        write(10+par%rank,*)'salobc local ',salobc_t(i1,j,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      if(mask_t(i2,j,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3d(i2,j,k,2)/=sal_t(i2,j,k,1)) then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i2,j,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(est)
        write(10+par%rank,*)'coordonnées',2,j
        write(10+par%rank,*)'sal_t ',anyv3d(i2,j,k,2),sal_t(i2,j,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i2,j,k,2)-sal_t(i2,j,k,1)
        write(10+par%rank,*)'salobc local ',salobc_t(i2,j,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! j

      enddo ! kmax

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(ksecu==1) then
        write(6,'(a,i0,a)')   &
        ' mpi conservation error for sal_t Details in fort.' &
         ,10+par%rank,' available in RDIR/CONFIG directory'
        stop ' STOP in routine scalar_check_mpi_conservation!'
      endif

      end subroutine scalar_check_mpi_conservation

!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine scalars_global_conservation !27-03-16
      use module_principal ; use module_parallele
      implicit none
      double precision :: sum1_cumul=0.  & ! pour cumuler l'effet de omega sur bilan ssh
                         ,sum2_cumul=0.  & ! pour cumuler l'effet de omega sur bilan T
                         ,sum3_cumul=0.  & ! pour cumuler l'effet de omega sur bilan S
                         ,sum_ssh0=0.    & ! cumul ssh t0
                         ,sum_tem0=0.    & ! cumul tem t0
                         ,sum_sal0=0.      ! cumul sal t0

! Cette routine permet de verifier les proprietes de conservation du schema d'advection
! et des flux en surface
! Elle suppose que les flux exterieurs (obc, riviere) soient nuls

      if(iteration3d==0) then !000000>
       sum1=0. ; k=kmax
       do j=1,jmax ; do i=1,imax
        sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*ssh_int_w(i,j,1)
       enddo       ; enddo
       sum2=0. ; sum3=0.
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        sum2=sum2+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,1)*tem_t(i,j,k,1)
        sum3=sum3+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,1)*sal_t(i,j,k,1)
       enddo       ; enddo       ; enddo
       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       sum_ssh0=sum1glb ; sum_tem0=sum2glb ; sum_sal0=sum3glb
      endif                   !000000>

!.......................
! Conservation de la ssh
      sum1=0. ; sum2=0. ; sum3=0. ; k=kmax
      do j=1,jmax ; do i=1,imax
       sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)
       sum2=sum2+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*ssh_int_w(i,j,2)
       sum3=sum3  &
         -dti_fw*mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*omega_w(i,j,kmax+1,1)
      enddo       ; enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      sum1_cumul=sum1_cumul+sum3glb
      if(par%rank==0) then !ooooo>
       open(unit=3,file='tmp/ssh_mean',position='append')
       write(3,'(5(1x,e14.7))')elapsedtime_now/86400.      &
                              ,(sum2glb         )/sum1glb  &
                              ,(        sum_ssh0)/sum1glb  &
                              ,sum1_cumul/sum1glb
       close(3)
      endif                !ooooo>
!.......................


!.......................
! T et S moyen:
      sum1=0. ; sum2=0. ; sum3=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)
       sum2=sum2+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*tem_t(i,j,k,2)
       sum3=sum3+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*sal_t(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)

! Effet cumul omega, cumul flux radiatif
      sum5=0.  ; sum6=0.  ; k=kmax
      do j=1,jmax ; do i=1,imax
       sum5=sum5-dti_fw*mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*omega_w(i,j,kmax+1,1)*tem_t(i,j,kmax,1)
       sum6=sum6+dti_fw*mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j) &
               /(rho*cp)*(                         &
           (1.-albedo_w(i,j))*ssr_w(i,j,1)         &
                +snsf_w(i,j,1)                     &
                +slhf_w(i,j,1)                     &
                +sshf_w(i,j,1)                     &
           +heatrelax_w(i,j,1))                      
      enddo       ; enddo
      call mpi_allreduce(sum5,sum5glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum6,sum6glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      sum2_cumul=sum2_cumul+sum5glb
      sum3_cumul=sum3_cumul+sum6glb

! Note sur le cas du bilan de sel: il n'y a pas de flux de sel A travers la surface
! de l'ocean car d'une part la pluie n'est pas salee et d'autre l'evaporation n'enleve
! que de l'eau douce. Autrement dit la quantite de sel ne change pas du fait des
! echange A travers la surface. Ce qui change c'est le volume d'eau autrement dit la
! concentration de sel (mais pas la quantite de sel totale). Dans notre discretisation
! ceci se traduit par le fait que le flux de sel par omega est exactement opposE au flux par
! slhf+precipitation. Si le domaine est ferme la quantite de sel doit restee inchangee
! meme si precipi+slhf est non nul. Par contre le volume change (si ces derniers sont non nuls)
! et par consequent la concentration de sel change. Il n'y a donc pas lieu de comptabiliser
! l'echange de surface pour le sel.... MEme genres de commentaires en ce qui concerne l'apport
! par les fleuves dans la mesure oU la valeur amont source de sel est nulle, ce qui revient
! A ce que les fleuves n'apportent pas de sel. 

! Le raisonnement est peu different pour la temperature: le volume d'eau
! retirE ou ajoutE A la surface a une temperature (la pluie est chaude ou froide)
! et donc il y a potentiellement un effet de modification de la temperature de
! surface par la pluie par exemple. A defaut de connaitre la temperature de la pluie
! on suppose que celle ci a la meme temperature que la surface de l'ocean. Idem quand
! omega est dirigEe vers le haut. Par contre l'evaporation consomme de la chaleur et
! donc, en plus de l'apport/retrait de chaleur par omega, on enleve de la chaleur par
! le flux slhf....
! 

      if(par%rank==0) then !pppp>
       write(6,*)'Tmean Smean=',sum2glb/sum1glb,sum3glb/sum1glb
       open(unit=3,file='tmp/ts_mean',position='append')
       write(3,'(10(1x,e14.7))')elapsedtime_now/86400.           &
! Temperature:
                              ,sum2glb/sum1glb                   & ! T moyen
                              ,sum_tem0/sum1glb                  & ! T moyen t0
                              ,sum2_cumul/sum1glb                & ! Effet omega sur T moyen
                              ,sum3_cumul/sum1glb                & ! Effet flux radiatif sur T moyen
! salinite:
                              ,sum3glb/sum1glb                   & ! S moyen
                              ,sum3glb                           & ! S cumul 
                              ,sum_sal0                            ! S cumul t0
       close(3)
      endif                !pppp>
      
      end subroutine scalars_global_conservation !27-03-16
#endif
!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine scalars_ts_minmax(t_) !11-01-16
      use module_principal ; use module_parallele
      implicit none
      integer t_

       x1=9999. ; x2=-x1 ; x3=9999. ; x4=-x1 
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
       if(mask_t(i,j,k)==1) then

        if(tem_t(i,j,k,t_)<x1) then
         x1=tem_t(i,j,k,t_) ; i1=i ; j1=j ; k1=k
        endif
        if(tem_t(i,j,k,t_)>x2) then
         x2=tem_t(i,j,k,t_) ; i2=i ; j2=j ; k2=k
        endif
        if(sal_t(i,j,k,t_)<x3) then
         x3=sal_t(i,j,k,t_) ; i3=i ; j3=j ; k3=k
        endif
        if(sal_t(i,j,k,t_)>x4) then
         x4=sal_t(i,j,k,t_) ; i4=i ; j4=j ; k4=k
        endif
       endif

       enddo ; enddo ; enddo

! Min et Max de T et S par proc ecrits dans des fichiers:
!      write(10+par%rank,*)iteration3d
!      write(10+par%rank,*)i1,j1,k1,x1
!      write(10+par%rank,*)i2,j2,k2,x2
!      write(10+par%rank,*)i3,j3,k3,x3
!      write(10+par%rank,*)i4,j4,k4,x4

      call mpi_allreduce(x1,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN T=',x0
      if(x1==x0) then !>>>
           write(6,*)'MIN T val',iteration3d,tem_t(i1,j1,k1,t_) !18-03-22
           write(6,*)'MIN T rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN T glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN T dz_t',dz_t(i1,j1,k1,t_)
           write(6,*)'MIN T h_w',h_w(i1,j1)
           write(6,*)'MIN T wetmask_t',wetmask_t(i1,j1) !06-04-22
      endif           !>>>

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX T=',x0
      if(x2==x0) then !>>>
           write(6,*)'MAX T val',iteration3d,tem_t(i2,j2,k2,t_)
           write(6,*)'MAX T rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX T glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX T dz_t',dz_t(i2,j2,k2,t_)
           write(6,*)'MAX T h_w',h_w(i2,j2)
           write(6,*)'MAX T wetmask_t',wetmask_t(i2,j2)
      endif           !>>>

      call mpi_allreduce(x3,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN S=',x0
      if(x3==x0) then !>>>
           write(6,*)'MIN S val',iteration3d,sal_t(i3,j3,k3,t_)
           write(6,*)'MIN S rank i j k',par%rank,i3,j3,k3
           write(6,*)'MIN S glob i j k',i3+par%timax(1),j3+par%tjmax(1),k3
           write(6,*)'MIN S dz_t',dz_t(i3,j3,k3,t_)
           write(6,*)'MIN S h_w',h_w(i3,j3)
           write(6,*)'MIN S wetmask_t',wetmask_t(i3,j3)
      endif           !>>>

      call mpi_allreduce(x4,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX S=',x0
      if(x4==x0) then !>>>
           write(6,*)'MAX S val',iteration3d,sal_t(i4,j4,k4,t_)
           write(6,*)'MAX S rank i j k',par%rank,i4,j4,k4
           write(6,*)'MAX S glob i j k',i4+par%timax(1),j4+par%tjmax(1),k4
           write(6,*)'MAX S dz_t',dz_t(i4,j4,k4,t_)
           write(6,*)'MAX S h_w',h_w(i4,j4)
           write(6,*)'MAX S wetmask_t',wetmask_t(i4,j4)
      endif           !>>>

      end subroutine scalars_ts_minmax
#endif
!-----------------------------------------------------------------------------

#ifdef bidon
      subroutine scalars_tsobc_minmax(t_) !11-01-16
      use module_principal ; use module_parallele
      implicit none
      integer t_

       x1_r4=9999. ; x2_r4=-x1_r4 ; x3_r4=9999. ; x4_r4=-x1_r4 
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
       if(mask_t(i,j,k)==1) then

        if(temobc_t(i,j,k,t_)<x1_r4) then
         x1_r4=temobc_t(i,j,k,t_) ; i1=i ; j1=j ; k1=k
        endif
        if(temobc_t(i,j,k,t_)>x2_r4) then
         x2_r4=temobc_t(i,j,k,t_) ; i2=i ; j2=j ; k2=k
        endif
        if(salobc_t(i,j,k,t_)<x3_r4) then
         x3_r4=salobc_t(i,j,k,t_) ; i3=i ; j3=j ; k3=k
        endif
        if(salobc_t(i,j,k,t_)>x4_r4) then
         x4_r4=salobc_t(i,j,k,t_) ; i4=i ; j4=j ; k4=k
        endif
       endif

       enddo ; enddo ; enddo

! Min et Max de T et S par proc ecrits dans des fichiers:
!      write(10+par%rank,*)iteration3d
!      write(10+par%rank,*)i1,j1,k1,x1_r4
!      write(10+par%rank,*)i2,j2,k2,x2_r4
!      write(10+par%rank,*)i3,j3,k3,x3_r4
!      write(10+par%rank,*)i4,j4,k4,x4_r4

      call mpi_allreduce(x1_r4,x0_r4,1,mpi_real,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN OBC T=',x0_r4
      if(x1_r4==x0_r4) then !>>>
           write(6,*)'MIN OBC T val',x0_r4,temobc_t(i1,j1,k1,t_)
           write(6,*)'MIN OBC T rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN OBC T glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN OBC T dz_t',dz_t(i1,j1,k1,t_)
           write(6,*)'MIN OBC T h_w',h_w(i1,j1)
      endif           !>>>

      call mpi_allreduce(x2_r4,x0_r4,1,mpi_real,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX OBC T=',x0_r4
      if(x2_r4==x0_r4) then !>>>
           write(6,*)'MAX OBC T val',x0_r4,temobc_t(i2,j2,k2,t_)
           write(6,*)'MAX OBC T rank i j k mask',par%rank,i2,j2,k2,mask_t(i2,j2,k2)
           write(6,*)'MAX OBC T glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX OBC T dz_t',dz_t(i2,j2,k2,t_)
           write(6,*)'MAX OBC T h_w',h_w(i2,j2)
      endif           !>>>

      call mpi_allreduce(x3_r4,x0_r4,1,mpi_real,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN OBC S=',x0_r4
      if(x3_r4==x0_r4) then !>>>
           write(6,*)'MIN OBC S val',x0_r4,salobc_t(i3,j3,k3,t_)
           write(6,*)'MIN OBC S rank i j k',par%rank,i3,j3,k3
           write(6,*)'MIN OBC S glob i j k',i3+par%timax(1),j3+par%tjmax(1),k3
           write(6,*)'MIN OBC S dz_t',dz_t(i3,j3,k3,t_)
           write(6,*)'MIN OBC S h_w',h_w(i3,j3)
      endif           !>>>

      call mpi_allreduce(x4_r4,x0_r4,1,mpi_real,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX OBC S=',x0_r4
      if(x4_r4==x0_r4) then !>>>
           write(6,*)'MAX OBC S val',x0_r4,salobc_t(i4,j4,k4,t_)
           write(6,*)'MAX OBC S rank i j k',par%rank,i4,j4,k4
           write(6,*)'MAX OBC S glob i j k',i4+par%timax(1),j4+par%tjmax(1),k4
           write(6,*)'MAX OBC S dz_t',dz_t(i4,j4,k4,t_)
           write(6,*)'MAX OBC S h_w',h_w(i4,j4)
      endif           !>>>

      end subroutine scalars_tsobc_minmax
#endif

!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine test_advection_bio !14-07-14
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='test_advection_bio'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!     flag=1 ! Cas 2D horizontal
      flag=0 ! Cas 2D vertical OXZ (eventuellement reduire jglb pour gagner cpu)
! ATTENTION CAS 2D OXZ Visualiser avec ncview en dehors des C.L. qui restent figees.....

      ub4=ubound(bio_t)
      lb4=lbound(bio_t)
      if(ub4(2)<jmax+2)stop 'Erreur max dimension 2 bio_t'
      if(lb4(2)>-1)    stop 'Erreur min dimension 2 bio_t'
      if(ub4(1)<imax+2)stop 'Erreur max dimension 1 bio_t'
      if(lb4(1)>-1)    stop 'Erreur min dimension 1 bio_t'

! initialisation bio_t(1)
      vb=1

      x1=real(iglb)/4.
      x2=real(jglb)/2.
      x3=0. !x3=0.2 ! Amplitude du bruit

      write(6,*)'Position initiale=',x1,x2

      do k=1,kmax
      do i1=-1,iglb+2
      do j1=-1,jglb+2
       i=i1-par%timax(1)
       j=j1-par%tjmax(1)
       if(i>=-1.and.i<=imax+2.and. &
          j>=-1.and.j<=jmax+2) then !--------------->

      if(flag==1) then ! 2D horizontal

! Initialiser un tourbillon chapeau pointu
         dist0=sqrt( (real(i1)-x1)**2 + (real(j1)-x2)**2 )
         dist1=max(5.-dist0,zero)/5.
         bio_t(i,j,k,:)=dist1 ! +x3*cos(pi*i1)*cos(pi*j1)*dist1 ! + du bruit

! Initialiser un front/tourbillon gaussien
!        dist0=sqrt( (real(i1)-x1)**2 + (real(j1)-x2)**2 )*dxb
!        dist0=(real(i1)-x1)*dxb
!        bio_t(i,j,k,:)=exp(-(dist0**2)/(6.*50.*280.45))

      endif            ! 2D horizontal

      if(flag==0) then ! 2D vertical OXZ

         dist0=sqrt( (real(i1)-x1)**2 + (real(k)-0.5*real(kmax))**2 )
!        dist1=max(20.-dist0,zero)/20.
         dist1=    min(max( real(kmax-10-k) , 0.),20.)/20.
         bio_t(i,j,k,:)=sin(pi*dist1)**2

!        if(vb==1.and.i==imax/2.and.j==jmax/2)write(6,*)k,dist1     ,sin(pi*dist1)**2

!        dist0=sqrt( (real(i1)-x1-5)**2 + (real(k)-0.5*real(kmax))**2 )
!        dist1=max(5.-dist0,zero)/5.
!        bio_t(i,j,k,vbmax)=dist1

      endif            ! 2D vertical OXZ


! Initialiser un front linéaire
!        x1=real(i1)/real(iglb)+real(j1)/real(jglb)
!        bio_t(i,j,k,:)=x1
!        x1=  real(j1-i1+iglb) / real(iglb+jglb) - 0.5
!        bio_t(i,j,k,:)=tanh(20.*x1)

       endif                        !--------------->

      enddo
      enddo
      enddo

      call graph_out

      if(flag==1) then !------->
      j=jmax/2 ; k=1
      do i=1,imax
       write(65,*)i,bio_t(i,j,k,1)
      enddo
      i=imax/4 ; k=1
      do j=1,jmax
       write(68,*)j,bio_t(i,j,k,1)
      enddo
      endif            !------->

      if(flag==0) then !------->
       j=jmax/2 ; i=imax/4
       do k=1,kmax
        write(65,*)k,bio_t(i,j,k,1)
       enddo
      endif            !------->


! ANNULER LE MELANGE:
      kh_w=0.

! allouer correctement anyv3d:
           call advection_quickest_checkdim
           call advection_checkdim_anyv1d

! un tour en key
      key=1000
!     key=5000
!     key=300
!     key=4*250
!     key=4*150
! ITERATION3D doit commencer A zero pour l'allocation dynamique correcte de anyv3d
!     do 2000 iteration3d=0,key/2      ! quart de un tour
!     do 2000 iteration3d=0,0           ! quart de un tour
      do 2000 iteration3d=1,500         ! quart de un tour


!     do 2000 iteration3d=0,1*key/4     ! 3/4 de tour
!     do 2000 iteration3d=0,key         ! un tour

      elapsedtime_now=elapsedtime_now+dti_fw

!     sum1=0.
!     sum2=0.
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      sum1=sum1+dxdy_t(i,j)*0.5*( dz_t(i,j,k,0) &
!                                 +dz_t(i,j,k,0)) 
!      sum2=sum2+dxdy_t(i,j)*0.5*( dz_t(i,j,k,0) &
!                                 +dz_t(i,j,k,0))&
!                                *bio_t(i,j,k,1)
!     enddo ; enddo       ; enddo
!     write(6,*)sum2/sum1
      

      if(par%rank==0)then
          if(mod(iteration3d,10)==0)write(6,*)iteration3d,key,elapsedtime_now
!         write(6,*)iteration3d,key
      endif

!..............................................................
!___________________________________________________
! Module du courant horizontal
!     x3=0.1*dxb/dti_fw ! Nb de courant * DX / DT
!     x3=0.99999*dxb/dti_fw ! Nb de courant * DX / DT
!     x3=1.99999*dxb/dti_fw ! Nb de courant * DX / DT
!     x3=0.5*dxb/dti_fw ! Nb de courant * DX / DT
      x3=0.
      if(iteration3d>20)x3=-x3

! Module du courant vertical (formule valide si DZ constant!)
      i=imax/2 ; j=jmax/2 ; k=kmax/2
!     x4=0.9999*h_w(i,j)*dsig_t(i,j,k)/dti_fw
!     x4=-1.99999    *h_w(i,j)*dsig_t(i,j,k)/dti_fw
      wsed(:)=-0.1*h_w(i,j)*dsig_t(i,j,k)/dti_fw
!     x4=-0.1   *h_w(i,j)*dsig_t(i,j,k)/dti_fw
      x4=0.
!     write(6,*)'DZ=',h_w(i,j)*dsig_t(i,j,k)
!     write(6,*)'DT=',dti_fw

! direction du courant:
      x2=2.*pi*real(iteration3d)/real(key)  ! rotation
!     x2=0.25*pi
!     x2=0.5*pi
!     x2=0.

! Marche pour 2DH et 2DV d'ou les boucles surdimensionnees dans 3 axes
      do k=1,kmax+1
       do i=0,imax+1
       do j=0,jmax+1

       if(flag==1) then ! 2D horizontal
! Test dans le plan horizontal OXY
!       vel_u(i,j,k,1)=  x3*sin(x2)
!       vel_v(i,j,k,1)=  x3*cos(x2)
!       if(iteration3d<=20) then
!        vel_u(i,j,k,1)= x3
!       else
!        vel_u(i,j,k,1)=-x3
!       endif
        vel_u(i,j,k,1)=x3
        vel_v(i,j,k,1)=x3
!       vel_u(i,j,k,1)=x3 
!       vel_v(i,j,k,1)=0.
        omega_w(i,j,k,1)=0.
       endif            ! 2D horizontal

       if(flag==0) then ! 2D Yvertical OXZ
! Test dans le plan vertical OXZ
!       vel_u(i,j,k,1)=  x3*sin(x2)
        vel_u(i,j,k,1)=x3
        omega_w(i,j,k,1)=x4 !*cos(x2)
        vel_v(i,j,k,1)=0.
       endif            ! 2D vertical OXZ

!       veldydz_u(i,j,k,1)=vel_u(i,j,k,1)*dy_u(i,j)*dz_u(i,j,k,1)       !30-09-09
!       veldxdz_v(i,j,k,1)=vel_v(i,j,k,1)*dx_v(i,j)*dz_v(i,j,k,1)       !30-09-09

       enddo
       enddo
      enddo

      do k=1,kmax+1 ; do i=1,imax+1 ; do j=0,jmax+1
        veldydz_u(i,j,k,1)=vel_u(i,j,k,1)*dy_u(i,j)*dz_u(i,j,k,1)      
      enddo         ; enddo         ; enddo
      do k=1,kmax+1 ; do i=0,imax+1 ; do j=1,jmax+1
        veldxdz_v(i,j,k,1)=vel_v(i,j,k,1)*dx_v(i,j)*dz_v(i,j,k,1)       
      enddo         ; enddo         ; enddo

!___________________________________________________
!..............................................................

      call advection_bio
      call mixsed_bio


 2000 continue

      if(flag==1) then !------->
       j=jmax/2 ; k=1
       do i=1,imax
        write(66,*)i,bio_t(i,j,k,1)
       enddo
       i=imax/4 ; k=1
       do j=1,jmax
        write(67,*)j,bio_t(i,j,k,1)
       enddo
      endif            !------->
      if(flag==0) then !------->
       j=jmax/2 ; i=imax/4
       do k=1,kmax
        write(66,*)k,bio_t(i,j,k,1)
       enddo
      endif            !------->


      call graph_out
      write(6,*)'dans sbr test_advection_bio'
      write(6,*)'-----------------------'
      stop 'test_advection_bio'


      end subroutine test_advection_bio
#endif
!---------------------------------------------------------------------------------------------------------
#ifdef bidon
      subroutine scalars_spongelayer_marmara
      use module_principal ; use module_parallele
      implicit none
      real*4 :: c1_,c2_,c3_
      integer :: istr_,jstr_,iend_,jend_

! Ancien domaine
!     if(iglb/=1120.or.jglb/=865) & 
!     stop 'scalars_spongelayer_marmara: iglb/=1120.or.jglb/=865'
!     istr_=max(1056-par%timax(1),1)
!     iend_=min(1099-par%timax(1),imax)
!     jstr_=max( 247-par%tjmax(1),1)
!     jend_=min( 299-par%tjmax(1),jmax)
! Nouveau domaine
      if(iglb/=1670.or.jglb/=1200) & 
      stop 'scalars_spongelayer_marmara: iglb/=1670.or.jglb/=1200'
      istr_=max(1311-par%timax(1),1)
      iend_=min(1473-par%timax(1),imax)
      jstr_=max(1062-par%tjmax(1),1)
      jend_=min(1164-par%tjmax(1),jmax)


      if(elapsedtime_now<864000.)then  ! temps relax court les 10 premiers jours
      c1_= 1.-dti_lp/(86400.*3)
      else
      c1_= 1.-dti_lp/(86400.*20.)     ! ensuite 20 jours
      endif
      c2_=(1.-c1_)*    timeweightobc(trc_id)
      c3_=(1.-c1_)*(1.-timeweightobc(trc_id))

!     ksecu=0
      do k=1,kmax
      do j=jstr_,jend_
      do i=istr_,iend_

       if (depth_t(i,j,k)>-15.)then
       tem_t(i,j,k,2)=tem_t(i,j,k,2)*c1_           &
                  +temobc_t(i,j,k,2)*c2_           &
                  +temobc_t(i,j,k,0)*c3_
       else
       tem_t(i,j,k,2)=tem_t(i,j,k,2)*c1_           &
                  +14.5*(1.-c1_)           
       endif

       if (depth_t(i,j,k)>-15.)then
       sal_t(i,j,k,2)=sal_t(i,j,k,2)*c1_           &
                  +25.*(1.-c1_)
       else
       sal_t(i,j,k,2)=sal_t(i,j,k,2)*c1_           &
                  +38.6*(1.-c1_)
       endif

!     ksecu=1
      enddo
      enddo
      enddo
!     if(ksecu==1)write(10+par%rank,*)'Je sponge marmara'

      end subroutine scalars_spongelayer_marmara
#endif
!.......................................................................
#ifdef bidon
      subroutine scalars_eps_tken_minmax  !03-03-20
      use module_principal ; use module_parallele
      implicit none

       x1=9999. ; x2=-x1 ; x3=9999. ; x4=-x1 
       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       if(mask_t(i,j,k)==1) then

        if(epsa_w(i,j,k)<x1) then
         x1=epsa_w(i,j,k) ; i1=i ; j1=j ; k1=k
        endif
        if(epsa_w(i,j,k)>x2) then
         x2=epsa_w(i,j,k) ; i2=i ; j2=j ; k2=k
        endif
        if(tkea_w(i,j,k)<x3) then
         x3=tkea_w(i,j,k) ; i3=i ; j3=j ; k3=k
        endif
        if(tkea_w(i,j,k)>x4) then
         x4=tkea_w(i,j,k) ; i4=i ; j4=j ; k4=k
        endif
       endif

       enddo ; enddo ; enddo

! Min et Max de T et S par proc ecrits dans des fichiers:
!      write(10+par%rank,*)iteration3d
!      write(10+par%rank,*)i1,j1,k1,x1
!      write(10+par%rank,*)i2,j2,k2,x2
!      write(10+par%rank,*)i3,j3,k3,x3
!      write(10+par%rank,*)i4,j4,k4,x4

      call mpi_allreduce(x1,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN EPS=',x0
      if(x1==x0) then !>>>
           write(6,*)'MIN EPS val',x0,epsa_w(i1,j1,k1)
           write(6,*)'MIN EPS rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN EPS glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN EPS h_w',h_w(i1,j1)
      endif           !>>>

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX EPS=',x0
      if(x2==x0) then !>>>
           write(6,*)'MAX EPS val',x0,epsa_w(i2,j2,k2)
           write(6,*)'MAX EPS rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX EPS glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX EPS h_w',h_w(i2,j2)
      endif           !>>>

      call mpi_allreduce(x3,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN TKE=',x0
      if(x3==x0) then !>>>
           write(6,*)'MIN TKE val',x0,tkea_w(i3,j3,k3)
           write(6,*)'MIN TKE rank i j k',par%rank,i3,j3,k3
           write(6,*)'MIN TKE glob i j k',i3+par%timax(1),j3+par%tjmax(1),k3
           write(6,*)'MIN TKE h_w',h_w(i3,j3)
      endif           !>>>

      call mpi_allreduce(x4,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX TKE=',x0
      if(x4==x0) then !>>>
           write(6,*)'MAX TKE val',x0,tkea_w(i4,j4,k4)
           write(6,*)'MAX TKE rank i j k',par%rank,i4,j4,k4
           write(6,*)'MAX TKE glob i j k',i4+par%timax(1),j4+par%tjmax(1),k4
           write(6,*)'MAX TKE h_w',h_w(i4,j4)
      endif           !>>>

      end subroutine scalars_eps_tken_minmax
#endif
!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine scalars_instabilite_max(t_) !19-02-22
      use module_principal ; use module_parallele
      implicit none
      integer t_

! Trouve le point d'instabilite statique maximale

       x2=-999.
       do j=1,jmax ; do i=1,imax
       if(mask_t(i,j,kmax)==1) then

        do k=kmin_w(i,j)+1,kmax

         if(   (rhp_t(i,j,k)  -rhp_t(i,j,k-1))  &
            /(depth_t(i,j,k)-depth_t(i,j,k-1))  >x2) then
 
            x2=(rhp_t(i,j,k)  -rhp_t(i,j,k-1))  &
            /(depth_t(i,j,k)-depth_t(i,j,k-1)) 
 
           i2=i ; j2=j ; k2=k
         endif

        enddo
       endif
       enddo ; enddo

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
      if(x2==x0) then !>>>
           x2=(rhp_t(i2,j2,k2)  -rhp_t(i2,j2,k2-1))  &
           /(depth_t(i2,j2,k2)-depth_t(i2,j2,k2-1)) 
           write(6,*)'MAX INSTAB drhp/dz',x0,x2
           write(6,*)'MAX INSTAB drhp/dk',rhp_t(i2,j2,k2)-rhp_t(i2,j2,k2-1)
           write(6,*)'MAX INSTAB rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX INSTAB glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX INSTAB 2 rhp',rhp_t(i2,j2,k2-1)  ,rhp_t(i2,j2,k2)
           write(6,*)'MAX INSTAB 2 tem',tem_t(i2,j2,k2-1,1),tem_t(i2,j2,k2,1)
           write(6,*)'MAX INSTAB 2 sal',sal_t(i2,j2,k2-1,1),sal_t(i2,j2,k2,1)
      endif           !>>>

      end subroutine scalars_instabilite_max
#endif
!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine test_advection_vel_u !02-01-23
      use module_principal
      use module_parallele
      implicit none
      real :: current_number_=0.1
      integer loop_,looprank_
#ifdef synopsis
       subroutinetitle='test_advection_bio'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      do loop_=3,3 ! 25,25
!     current_number_=real(loop_)/20.
      if(loop_==1)current_number_=1.999999
      if(loop_==2)current_number_=0.999999
      if(loop_==3)current_number_=0.499999
      if(loop_==4)current_number_=0.2499999

      flag=1 ! Cas 2D horizontal
!     flag=0 ! Cas 2D vertical OXZ (eventuellement reduire jglb pour gagner cpu)
! ATTENTION CAS 2D OXZ Visualiser avec ncview en dehors des C.L. qui restent figees.....

      ub4=ubound(vel_u)
      lb4=lbound(vel_u)
      if(ub4(2)<jmax+1)stop 'Erreur max dimension 2 vel_u'
      if(lb4(2)>0)     stop 'Erreur min dimension 2 vel_u'
      if(ub4(1)<imax+2)stop 'Erreur max dimension 1 vel_u'
      if(lb4(1)>0)     stop 'Erreur min dimension 1 vel_u'

      x1=real(iglb)/4.
      x2=real(jglb)/2.
      x3=0. !x3=0.2 ! Amplitude du bruit

      write(6,*)'Position initiale=',x1,x2

      vel_u=0.
      vel_v=0.
      do k=1,kmax
      do i1=0,iglb+2
      do j1=0,jglb+1
       i=i1-par%timax(1)
       j=j1-par%tjmax(1)
       if(i>=0.and.i<=imax+2.and. &
          j>=0.and.j<=jmax+1) then !--------------->

      if(flag==1) then ! 2D horizontal
! Initialiser un tourbillon:
!        dist0=sqrt( (real(i1)-x1)**2 + (real(j1)-x2)**2 )
!        dist1=max(10.-dist0,zero)/10.
!        vel_u(i,j,k,:)=dist1+x3*cos(pi*i1)*cos(pi*j1)*dist1
! Cas unidirectionnel Oi
!        vel_u(i,j,k,:)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
! Cas unidirectionnel Oj
!        vel_u(i,j,k,:)=exp(-( (real(j1)-real(jglb)/2.0)/5.)**2)
         vel_u(i,j,k,:)=0.

      endif

      if(flag==0) then ! 2D vertical OXZ
         dist0=sqrt( (real(i1)-x1)**2 + (real(k)-0.25*real(kmax))**2 )
         dist1=max(10.-dist0,zero)/10.
         vel_u(i,j,k,:)=dist1
      endif

       endif                        !--------------->

       if(i>=0.and.i<=imax+1.and. &
          j>=0.and.j<=jmax+2) then !--------------->


!        vel_v(i,j,k,:)=0.
! Cas unidirectionnel Oi
         vel_v(i,j,k,:)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
! Cas unidirectionnel Oj
!        vel_v(i,j,k,:)=exp(-( (real(j1)-real(jglb)/2.0)/5.)**2)


       endif                        !--------------->

      enddo
      enddo
      enddo

! POUR VERIFIER QU'UN CHAMPS CONSTANt RESTE INCHANGE
!     vel_u=1.
!     vel_v=1.

       if(flag==1) then !111>
        k=kmax
        i=imax/2 
        j=jmax/2 
        do looprank_=0,nbdom-1
         if(par%rank==looprank_) then !RANK>
          if(par%rank==0) then
           open(unit=66,file='vel_init')
          else
           open(unit=66,file='vel_init',position='append')
          endif
!         do i=2,imax
!          write(66,*)i+par%timax(1),vel_u(i,j,k,2)
!         enddo
!         do j=2,jmax-1
!          write(66,*)j+par%tjmax(1),vel_u(i,j,k,2)
!         enddo
          do i=2,imax-1
           write(66,*)i+par%timax(1),vel_v(i,j,k,2)
          enddo
!         do j=2,jmax
!          write(66,*)j+par%tjmax(1),vel_v(i,j,k,2)
!         enddo
          close(66)
         endif                        !RANK>
#ifdef parallele
        call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
         enddo ! looprank_
       endif            !111> 


      if(flag==0) then !000> vertical
       i=nint(real(iglb)/4.) ; j=jmax/2
       open(unit=66,file='vel_init')
       do k=1,kmax
        write(66,*)k,vel_u(i,j,k,2),vel_v(i,j,k,2)
       enddo
       close(66)
      endif            !000> vertical

      sum1=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
       sum1=sum1+mask_i_u(i)*mask_j_u(j)*mask_u(i,j,k)*dxdy_u(i,j)*dz_u(i,j,k,2)*vel_u(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme vel_u initiale',sum2
      sum1=0.
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
       sum1=sum1+mask_i_v(i)*mask_j_v(j)*mask_v(i,j,k)*dxdy_v(i,j)*dz_v(i,j,k,2)*vel_v(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme vel_v initiale',sum2

! un tour en key
!     key=1000
      key=5000
!     key=300
!     key=4*250
!     key=4*150
!     do 2000 iteration3d=0,key/2      ! quart de un tour
!     do 2000 iteration3d=0,0           ! quart de un tour
!     do 2000 iteration3d=0,3*key/4     ! 3/4 de tour
!     do 2000 iteration3d=0,key         ! un tour
      elapsedtime_now=0.
      call graph_out
!     do 2000 iteration3d=1,nint(20./current_number_)
      do 2000 iteration3d=1,iteration3d_max
      elapsedtime_now=elapsedtime_now+dti_fw

      if(par%rank==0)then
          if(mod(iteration3d,1 )==0)write(6,*)iteration3d,key,elapsedtime_now
!         write(6,*)iteration3d,key
      endif

!..............................................................
!___________________________________________________
! Module du courant horizontal
      if(flag==1) then
       x3=current_number_*dxb/dti_lp ! Nb de courant * DX / DT
       x4=0.
      endif
      if(flag==0) then
       x3=0.
       x4=current_number_*h_w(i,j)*dsig_t(i,j,k)/dti_lp
      endif

! direction du courant:
!     x2=2.*pi*real(iteration3d)/real(key)  ! rotation
!     x2=0.25*pi
!     x2=0.5*pi
      x2=0.

! Vitesses advectantes:
      if(flag==1) then !000>
       omega_w=0.
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
!        veldydz_u(i,j,k,1)=x3*sin(x2)*dy_u(i,j)*dz_u(i,j,k,1)      
!        veldydz_u(i,j,k,1)=x3        *dy_u(i,j)*dz_u(i,j,k,1)      
         veldydz_u(i,j,k,1)=+1.9999*dx_u(i,j)*dy_u(i,j)*dz_u(i,j,k,1)/dti_lp
!        veldydz_u(i,j,k,1)=0.
       enddo       ; enddo       ; enddo
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!        veldxdz_v(i,j,k,1)=x3*cos(x2)*dx_v(i,j)*dz_v(i,j,k,1)   
!        veldxdz_v(i,j,k,1)=x3        *dx_v(i,j)*dz_v(i,j,k,1)   
!        veldxdz_v(i,j,k,1)=+1.9999*dx_v(i,j)*dy_v(i,j)*dz_v(i,j,k,1)/dti_lp
         veldxdz_v(i,j,k,1)=0. 
       enddo       ; enddo       ; enddo
      endif            !000>

      if(flag==0) then !111>
         veldydz_u=0.
         veldxdz_v=0.
         do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
          omega_w(i,j,k,1)=x4
         enddo       ; enddo       ; enddo
      endif            !111>

#ifdef bidon
! Dans l'hypothese ou veldtodx veldtody requis:
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
! u*dt/dx:
        veldtodx_u(i,j,k,1)=veldydz_u(i,j,k,1)*dti_lp  &
                            /( dxdy_u(i,j)             &
                                *dz_u(i,j,k,1) )
      enddo
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
! v*dt/dy
        veldtody_v(i,j,k,1)=veldxdz_v(i,j,k,1)*dti_lp  &
                            /( dxdy_v(i,j)             &
                                *dz_v(i,j,k,1) )
      enddo
      enddo
      enddo
#endif


!     write(6,*)veldtodx_u(imax/2,jmax/2,kmax/2,1)
!     write(6,*)veldtody_v(imax/2,jmax/2,kmax/2,1)
!     stop 'koko'

! Termes annulEs dans momentum_equation:
           cdb_t=0.
!          vel_v=0.
      coriolis_t=0.
            km_w=0.
       wetmask_u=1.
       wetmask_v=1.
       wetmask_t=1.
       gradssh_u=0.
       gradssh_v=0.
      presgrad_u=0.
      presgrad_v=0.
      velobc_u=0.
      velobc_v=0.
      sponge_u=0.
      sponge_v=0.
      tfilterfb=0.
      tfilterlf=0.
 

!___________________________________________________
!..............................................................

      call momentum_equations('after external mode ')
      if(timestep_type==timestep_leapfrog)call obc_int_mpi(2,1) !28-03-16
      if(timestep_type==timestep_forwbckw)call obc_int_mpi(2,12)

      vel_u(:,:,:,-1)=vel_u(:,:,:,0)
      vel_u(:,:,:,0)= vel_u(:,:,:,1)
      vel_u(:,:,:,1)= vel_u(:,:,:,2)

      vel_v(:,:,:,-1)=vel_v(:,:,:,0)
      vel_v(:,:,:,0)= vel_v(:,:,:,1)
      vel_v(:,:,:,1)= vel_v(:,:,:,2)

 2000 continue

       if(flag==1) then !111>
        k=kmax
        i=imax/2 
        j=jmax/2 
        do looprank_=0,nbdom-1
         if(par%rank==looprank_) then !RANK>
          if(par%rank==0) then
           open(unit=66,file='vel_fin')
          else
           open(unit=66,file='vel_fin',position='append')
          endif
!         do i=2,imax
!          write(66,*)i+par%timax(1),vel_u(i,j,k,2)
!         enddo
!         do j=2,jmax-1
!          write(66,*)j+par%tjmax(1),vel_u(i,j,k,2)
!         enddo
          do i=2,imax-1
           write(66,*)i+par%timax(1),vel_v(i,j,k,2)
          enddo
!         do j=2,jmax
!          write(66,*)j+par%tjmax(1),vel_v(i,j,k,2)
!         enddo
          close(66)
         endif                        !RANK>
#ifdef parallele
        call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
         enddo ! looprank_
       endif            !111> 

      if(flag==0) then !000> vertical
       i=nint(real(iglb)/4.) ; j=jmax/2
       write(texte90,'(a,i0)')'vel_',loop_
       open(unit=66,file=texte90)
       do k=1,kmax
        write(66,*)k,vel_u(i,j,k,2),vel_v(i,j,k,2)
       enddo
       close(66)
      endif            !000> vertical

!     write(6,*)'poid schema centrE:',2*anyv3d(imax/2,jmax/2,kmax,2)

      enddo ! fin de boucle sur loop_

      sum1=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
       sum1=sum1+mask_i_u(i)*mask_j_u(j)*mask_u(i,j,k)*dxdy_u(i,j)*dz_u(i,j,k,2)*vel_u(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme vel_u finale  ',sum2
      sum1=0.
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
       sum1=sum1+mask_i_v(i)*mask_j_v(j)*mask_v(i,j,k)*dxdy_v(i,j)*dz_v(i,j,k,2)*vel_v(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme vel_v finale  ',sum2

      call graph_out
      write(6,*)'dans sbr test_advection_bio'
      write(6,*)'-----------------------'
      stop 'test_advection_vel_u'


      end subroutine test_advection_vel_u

#endif
!---------------------------------------------------------------------------------------------------------
#ifdef bidon
      subroutine test_advection_temsal !02-01-23
      use module_principal
      use module_parallele
      implicit none
      real :: current_number_=0.1
      integer loop_,loopmax_,looprank_

! Derniere mise A jour le 23-08-21 pour tester le schema vertical Lax-Wendroff avec limiteur d'oscillation et diffusion negative

! UNE CHOSE IMPORTANTE POUR COMPRENDRE CERTAINS TEST.
! Nous sommes en Leap-Frog. Le nombre de courant s'exprime dx (ou dz) divise par dti_lp soit 2 fois dt
! Autrement dit si le nombre de courant est 1 et que nous faisons 10 iterations, la forme se deplace
! de 10*dt*dx/(2dt) soit 5dx, donc 5 mailles (et non pas 10 comme intuitivement)
! Attention aux initialisation un peu simple (T(0)=T(1)) et au leap-frog qui supervise 2 trajectoires,
! une paire, une impaire. iParfois il faut comparer les paires entre elles, les impaires entre elles.

      do loop_=1,1 ! 25,25


      flag=1 ! Cas 2D horizontal
!     flag=0 ! Cas 2D vertical OXZ (eventuellement reduire jglb pour gagner cpu)
! ATTENTION CAS 2D OXZ Visualiser avec ncview en dehors des C.L. qui restent figees.....

      x1=real(iglb)/4.
      x2=real(jglb)/2.
      x3=0. !x3=0.2 ! Amplitude du bruit

!     write(6,*)'Position initiale=',x1,x2

      do k=1,kmax

       do i1=-1,iglb+2
       do j1=-1,jglb+2
        i=i1-par%timax(1)
        j=j1-par%tjmax(1)

         if(i>=-1.and.i<=imax+2.and.&
            j>=-1.and.j<=jmax+2) then !pmx> !19-08-21>

         if(flag==1) then !>>>
! Cas unidirectionnel Oi
!        if(i1>iglb/2-10.and.i1<iglb/2+10) then
!         tem_t(i,j,k,:)=1.
!         sal_t(i,j,k,:)=1.
!        else
!         tem_t(i,j,k,:)=0.
!         sal_t(i,j,k,:)=0.
!        endif
          tem_t(i,j,k,:)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
          sal_t(i,j,k,:)=35.
!         tem_t(i,j,k,:)=15.
!         sal_t(i,j,k,:)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
! Cas unidirectionnel Oj
!         tem_t(i,j,k,:)=exp(-( (real(j1)-real(jglb)/2.0)/5.)**2)
!         sal_t(i,j,k,:)=35.
! Cas unidirectionnel 2D
!         tem_t(i,j,k,:)=exp(-( (real(i1)-real(iglb)/4.0)**2 &
!                              +(real(j1)-real(jglb)/2.0)**2 &
!                             )/5.) 
!         sal_t(i,j,k,:)=0.
!        if(j1>jglb/2-10.and.j1<jglb/2+10) then
!         tem_t(i,j,k,:)=1.
!         sal_t(i,j,k,:)=1.
!        else
!         tem_t(i,j,k,:)=0.
!         sal_t(i,j,k,:)=0.
!        endif
! Pour verifier conservation:
         endif            !>>>


         if(flag==0) then !>>>
! Tester l'advection verticale:
! Une gaussienne:
!        dist1=    min(max( real(k-10) , 0.),20.)/20.
!        tem_t(i,j,k,:)=sin(pi*dist1)**2
!        sal_t(i,j,k,:)=sin(pi*dist1)**2
! Un front:
!        if(k>kmax/2) then
!        if(depth_t(i,j,k)>-hmax/2) then  
         if(depth_t(i,j,k)>-hmax/2.and. &
            depth_t(i,j,k)<-hmax/4) then
          tem_t(i,j,k,:)=1.
          sal_t(i,j,k,:)=1.
         else
          tem_t(i,j,k,:)=0.
          sal_t(i,j,k,:)=0.
         endif

         endif            !>>>

         endif                        !pmx> !19-08-21>
       enddo
       enddo

       do i1=0,iglb+1
       do j1=0,jglb+1
        i=i1-par%timax(1)
        j=j1-par%tjmax(1)
        if(i>=0.and.i<=imax+1.and.&
           j>=0.and.j<=jmax+1) then !pmx> !19-08-21>
         temobc_t(i,j,k,:)=tem_t(i,j,k,1)
         salobc_t(i,j,k,:)=sal_t(i,j,k,1)
        endif                       !pmx> !19-08-21>
       enddo
       enddo

      enddo ! boucle k

      if(flag==1) then !>>>

      sum1=0.
      sum3=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
         sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*tem_t(i,j,k,2)
         sum3=sum3+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*sal_t(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme temperature initiale',sum2
      if(par%rank==0)write(6,*)'Somme salinite    initiale',sum4

      k=kmax
      i=imax/2 
      j=jmax/2 
      do looprank_=0,nbdom-1
       if(par%rank==looprank_) then !RANK>
        if(par%rank==0) then
         open(unit=66,file='tem_init')
        else
         open(unit=66,file='tem_init',position='append')
        endif
        do i=2,imax-1
         write(66,*)i+par%timax(1),real(tem_t(i,j,k,1)),real(sal_t(i,j,k,1))
        enddo
!       do j=2,jmax-1
!        write(66,*)j+par%tjmax(1),real(tem_t(i,j,k,1)),real(sal_t(i,j,k,1))
!       enddo
        close(66)
       endif                        !RANK>
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      enddo ! looprank_
      endif            !>>>

      if(flag==0) then !>>>
      i=nint(real(iglb)/2.)-par%timax(1)
      j=nint(real(jglb)/2.)-par%tjmax(1)
      if(i>1.and.i<imax.and.&
         j>1.and.j<jmax) then !pmx> !19-08-21>
      write(texte30,'(a,i0)')'tem_init_rank',par%rank
      open(unit=66,file=trim(texte30))
      do k=1,kmax
       write(66,*)depth_t(i,j,k),tem_t(i,j,k,1),sal_t(i,j,k,1)
      enddo
      close(66)
      endif                   !pmx> !19-08-21>
      endif            !>>>

! Pour que les C.L. fonctionnent:
!      salobc_t(0:imax+1,0:jmax+1,1:kmax+1,0:2) &
!        =sal_t(0:imax+1,0:jmax+1,1:kmax+1,0:2)  
!      temobc_t(0:imax+1,0:jmax+1,1:kmax+1,0:2) &
!        =tem_t(0:imax+1,0:jmax+1,1:kmax+1,0:2)  


      call graph_out
! un tour en key
      elapsedtime_now=0
      key=1000
!     key=5000
!     key=300
!     key=4*250
!     key=4*150
!     do 2000 iteration3d=0,key/2      ! quart de un tour
!     do 2000 iteration3d=1-1,key-1           ! quart de un tour
!     do 2000 iteration3d=0,3*key/4     ! 3/4 de tour
!     do 2000 iteration3d=0,key         ! un tour
!     do 2000 iteration3d=1,nint(20./current_number_)

!     do 2000 iteration3d=1,nint(100./dti_fw/0.001)

      if(flag==0) then !00000000000>
       veldydz_u=0.
       veldxdz_v=0.
       do k=1,kmax+1   ; do j=1,jmax ; do i=1,imax
        omega_w(i,j,k,1)=-0.0011
       enddo       ; enddo       ; enddo
!      omega_w(:,:,1     ,1)=0.
!      omega_w(:,:,kmax+1,1)=0.
      endif            !00000000000>
      if(flag==1) then !11111111111>
! Cas unidirectionnel Oi
      omega_w=0.
      veldxdz_v=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        veldydz_u(i,j,k,1)=+(0.9999*dx_u(i,j)/dti_lp)*dy_u(i,j)*dz_u(i,j,k,1)      
!       veldydz_u(i,j,k,1)=-(1.  *dx_u(i,j)/dti_lp)*dy_u(i,j)*dz_u(i,j,k,1)      
      enddo       ; enddo       ; enddo
! Cas unidirectionnel Oj
!     omega_w=0.
!     veldydz_u=0.
!     do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!       veldxdz_v(i,j,k,1)=-(1.0*dy_v(i,j)/dti_lp)*dx_v(i,j)*dz_v(i,j,k,1)      
!     enddo       ; enddo       ; enddo
      endif            !11111111111>

!     loopmax_=100000
!     loopmax_=300000

      do 2000 iteration3d=0,iteration3d_max

      elapsedtime_now=elapsedtime_now+dti_fw

      if(par%rank==0)then
!         if(mod(iteration3d,10)==0)write(6,*)iteration3d,key
          if(mod(iteration3d,10)==0)write(6,*)iteration3d
!         write(6,*)iteration3d
      endif


! Ces lignes servent A Activer l'aiguillage "anti-boules-d'eau-dense"
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
         rhp_t(i,j,k)=tem_t(i,j,k,1)
        enddo ; enddo ; enddo

#ifdef bidon
!..............................................................
!___________________________________________________
! Module du courant horizontal
      if(flag==1) then !1111111111111111111111111>
       x3=current_number_*dxb/dti_lp ! Nb de courant * DX / DT
!      write(6,*)'x3=',x3
! direction du courant:
      x2=2.*pi*real(iteration3d)/real(key)  ! rotation
      if(mod(iteration3d,50)==0)write(6,*)'sin,cos',sin(x2),cos(x2)
!     x2=0.25*pi
!     x2=0.5*pi
!     x2=0.
!  wsed(:)=-0.1*h_w(i,j)*dsig_t(i,j,k)/dti_fw
! Vitesses advectantes:
      omega_w=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        veldydz_u(i,j,k,1)=x3*sin(x2)*dy_u(i,j)*dz_u(i,j,k,1)      
!       veldydz_u(i,j,k,1)=x3        *dy_u(i,j)*dz_u(i,j,k,1)      
        vel_u(i,j,k,:)=veldydz_u(i,j,k,1)/dy_u(i,j)/dz_u(i,j,k,1)
      enddo       ; enddo       ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        veldxdz_v(i,j,k,1)=x3*cos(x2)*dx_v(i,j)*dz_v(i,j,k,1)   
!       veldxdz_v(i,j,k,1)=x3        *dx_v(i,j)*dz_v(i,j,k,1)   
        vel_v(i,j,k,:)=veldxdz_v(i,j,k,1)/dx_v(i,j)/dz_v(i,j,k,1)
      enddo       ; enddo       ; enddo
!     dz_t(:,:,:,2)=dz_t(:,:,:,0)
!     call omega
!     write(6,*)'w',omega_w(30,30,kmax+1,1)
! cas unidirectionnel Oi
!     if(modulo(iteration3d+1,1000)==0)veldydz_u=-veldydz_u
! cas unidirectionnel Oj
!     if(modulo(iteration3d+1,1000)==0)veldxdz_v=-veldxdz_v
!     do k=1,kmax 
!     do j=1,jmax
!     do i=1,imax
!      dz_t(i,j,k,2)=dz_t(i,j,k,0)+dti_lp*omega_w(i,j,kmax+1,1)/real(kmax)
!     enddo
!     enddo
!     enddo
!     call omega
!     write(6,*)'w',omega_w(30,30,kmax+1,1),dz_t(30,30,kmax,2)
!      i=imax/2 ; j=jmax/2 ; k=kmax/2
!      write(6,*)(veldydz_u(i,j,k,1)/dy_u(i,j)/dz_u(i,j,k,1))*dti_fw/dx_u(i,j)
!      stop 'coco second'
      endif            !1111111111111111111111111>
#endif

!#ifdef bidon
      if(flag==0) then !00000000000>
!      veldydz_u=0.
!      veldxdz_v=0.
!      do k=2,kmax   ; do j=1,jmax ; do i=1,imax
!       omega_w(i,j,k,1)= current_number_*h_w(i,j)*dsig_t(i,j,k)/dti_lp
!       omega_w(i,j,k,1)=-0.0011
!      enddo       ; enddo       ; enddo
!      omega_w(:,:,1     ,1)=0.
!      omega_w(:,:,kmax+1,1)=0.
       if(modulo(iteration3d+1,500)==0)omega_w=-omega_w       
      endif            !00000000000>
!#endif

! Termes annulEs dans l'equation de tke de gaspar
            kh_w=0.
       wetmask_t=1.
       upwindriver_t=1.
       omega_evaprec_w=0.
#ifdef bidonref
       temref_t=0.
       temref_u=0.
       temref_v=0.
       salref_t=0.
       salref_u=0.
       salref_v=0.
#endif


!___________________________________________________
!..............................................................

!     write(66,*)'-----------'
!     write(67,*)'-----------'
         call vertmix_matrix
         call obc_scal_mpi_now
         call advection_scal                                        !06-05-10
         call vertmix_tem(1) ! arg=flagsolver_ 1 = computes matrix pretreatment
         call vertmix_sal(0)

      dz_t(:,:,:,0)=dz_t(:,:,:,1)
      dz_t(:,:,:,1)=dz_t(:,:,:,2)
      tem_t(:,:,:,-1)=tem_t(:,:,:,0)
      tem_t(:,:,:, 0)=tem_t(:,:,:,1)
      tem_t(:,:,:, 1)=tem_t(:,:,:,2)
      sal_t(:,:,:,-1)=sal_t(:,:,:,0)
      sal_t(:,:,:, 0)=sal_t(:,:,:,1)
      sal_t(:,:,:, 1)=sal_t(:,:,:,2)

! Verif de conservation:
!      sum1=0.
!      sum2=0.
!      do k=1,kmax
!      do j=1,jmax ; do i=1,imax
!      sum1=sum1+tem_t(i,j,k,1)*dxdy_t(i,j)*dz_t(i,j,k,1)
!      sum2=sum2+sal_t(i,j,k,1)*dxdy_t(i,j)*dz_t(i,j,k,1)
!        if(tem_t(i,j,k,2)/=sal_t(i,j,k,2)) then
!         write(6,*)i,j,k,tem_t(i,j,k,2),sal_t(i,j,k,2)
!         stop 'tem/=sal'
!        endif
!      enddo ; enddo
!      enddo

 
!      if(par%rank==0)write(68,*)iteration3d,tem_t(imax/2,jmax/2,kmax,2) 

!     if(flag==1) then !>>>
!      k=kmax
!      i=nint(real(iglb)/4.)
!      write(texte90,'(a,i0)')'tem_',loop_
!      open(unit=66,file=texte90)
!      x1=0.
!      x2=0.
!      do j=1,jmax
!       write(66,*)j,tem_t(i,j,k,1),sal_t(i,j,k,1)
!       x1=max(x1,tem_t(i,j,k,1))
!       x2=max(x2,sal_t(i,j,k,1))
!      enddo
!      close(66)
!     endif

      if(mod(iteration3d,10000)==0.or. &
                    iteration3d==iteration3d_max)then !ecrire fichier>

       if(flag==0) then !000>
        i=nint(real(iglb)/2.)-par%timax(1)
        j=nint(real(jglb)/2.)-par%tjmax(1)
         if(i>1.and.i<imax.and.&
            j>1.and.j<jmax) then !pmx> !19-08-21>
         write(texte90,'(a,i0)')'tem_',loop_
         open(unit=66,file=texte90)
!        x1=0.
!        x2=0.
         do k=1,kmax
!         write(66,*)depth_t(i,j,k),tem_t(i,j,k,1),sal_t(i,j,k,1)
          write(66,*)depth_t(i,j,k),tem_t(i,j,k,1),temf_t(i,j,k)/ratio_negdif_ver
!         x1=max(x1,tem_t(i,j,k,1))
!         x2=max(x2,sal_t(i,j,k,1))
         enddo
         close(66)
         endif                   !pmx> !19-08-21>
       endif            !000> 

       if(flag==1) then !111>
        k=kmax
        i=imax/2 
        j=jmax/2 
        do looprank_=0,nbdom-1
         if(par%rank==looprank_) then !RANK>
          if(par%rank==0) then
           open(unit=66,file='tem_fin')
          else
           open(unit=66,file='tem_fin',position='append')
          endif
          do i=2,imax-1
           write(66,*)i+par%timax(1),real(tem_t(i,j,k,1)),real(sal_t(i,j,k,1)),real(rhp_t(i,j,k))
          enddo
!         do j=2,jmax-1
!          write(66,*)j+par%tjmax(1),real(tem_t(i,j,k,1)),real(sal_t(i,j,k,1))
!         enddo
          close(66)
         endif                        !RANK>
#ifdef parallele
        call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
         enddo ! looprank_
       endif            !111> 

      endif                                    !ecrire fichier>

!      if(par%rank==0)write(68,*)iteration3d*dti_fw/86400 &
!       ,tem_t(imax/2,jmax/2,9,1) &
!      ,temf_t(imax/2,jmax/2,9)/ratio_negdif_ver

 2000 continue

!     write(67,*)current_number_,x1,x2

      enddo ! fin de boucle sur loop_

      sum1=0.
      sum3=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
         sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*tem_t(i,j,k,2)
         sum3=sum3+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*sal_t(i,j,k,2)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme temperature finale  ',sum2
      if(par%rank==0)write(6,*)'Somme salinite    finale  ',sum4

      k=kmax

      call graph_out
      write(6,*)'dans sbr test_advection_temsal'
      write(6,*)'-----------------------'
      stop 'test_advection_temsal'


      end subroutine test_advection_temsal
#endif
!.......................................................................
#ifdef bidon
      subroutine test_advection_keps !02-01-23
      use module_principal
      use module_parallele
      implicit none
      real :: current_number_=0.1
      integer loop_,loopmax_,looprank_

! Derniere mise A jour le 23-08-21 pour tester le schema vertical Lax-Wendroff avec limiteur d'oscillation et diffusion negative

! UNE CHOSE IMPORTANTE POUR COMPRENDRE CERTAINS TEST.
! Nous sommes en Leap-Frog. Le nombre de courant s'exprime dx (ou dz) divise par dti_lp soit 2 fois dt
! Autrement dit si le nombre de courant est 1 et que nous faisons 10 iterations, la forme se deplace
! de 10*dt*dx/(2dt) soit 5dx, donc 5 mailles (et non pas 10 comme intuitivement)
! Attention aux initialisation un peu simple (T(0)=T(1)) et au leap-frog qui supervise 2 trajectoires,
! une paire, une impaire. iParfois il faut comparer les paires entre elles, les impaires entre elles.

!     looplimit_hor=1
        emin=0.
        epsmin=0.

      do loop_=1,1 ! 25,25


      flag=1 ! Cas 2D horizontal
!     flag=0 ! Cas 2D vertical OXZ (eventuellement reduire jglb pour gagner cpu)
! ATTENTION CAS 2D OXZ Visualiser avec ncview en dehors des C.L. qui restent figees.....

      x1=real(iglb)/4.
      x2=real(jglb)/2.
      x3=0. !x3=0.2 ! Amplitude du bruit

!     write(6,*)'Position initiale=',x1,x2

      do k=2,kmax

       do i1=-1,iglb+2
       do j1=-1,jglb+2
        i=i1-par%timax(1)
        j=j1-par%tjmax(1)

         if(i>= 1.and.i<=imax  .and.&
            j>= 1.and.j<=jmax  ) then !pmx> !19-08-21>

         if(flag==1) then !>>>
! Cas unidirectionnel Oi
!         tkea_w(i,j,k)=15.
!         epsa_w(i,j,k)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
!         tkea_w(i,j,k)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
!         tken_w(i,j,k)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
!         epsa_w(i,j,k)=10.
          tken_w(i,j,k)=exp(-( (real(i1)-real(iglb)/2.0)/5.)**2)
! Cas unidirectionnel Oj
!         tkea_w(i,j,k)=15.
!         epsa_w(i,j,k)=exp(-( (real(j1)-real(jglb)/2.0)/5.)**2)
!         tkea_w(i,j,k)=exp(-( (real(j1)-real(jglb)/2.0)/5.)**2)
!         tken_w(i,j,k)=exp(-( (real(j1)-real(jglb)/2.0)/5.)**2)
!         epsa_w(i,j,k)=10.
         endif            !>>>

         endif                        !pmx> !19-08-21>
       enddo
       enddo

      enddo ! boucle k

      if(flag==1) then !>>>

      sum1=0.
      sum3=0.
      sum5=0.
      do k=2,kmax ; do j=1,jmax ; do i=1,imax
         sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*tkea_w(i,j,k)
         sum3=sum3+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*epsa_w(i,j,k)
         sum5=sum5+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*dz_t(i,j,k,2)*tken_w(i,j,k)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum5,sum6,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme tkea initiale',sum2
      if(par%rank==0)write(6,*)'Somme epsa initiale',sum4
      if(par%rank==0)write(6,*)'Somme tken initiale',sum6

      k=kmax
      i=imax/2 
      j=jmax/2 
      do looprank_=0,nbdom-1
       if(par%rank==looprank_) then !RANK>
        if(par%rank==0) then
         open(unit=66,file='tem_init')
        else
         open(unit=66,file='tem_init',position='append')
        endif
        do i=2,imax-1
         write(66,*)i+par%timax(1),real(tkea_w(i,j,k)),real(epsa_w(i,j,k)),real(tken_w(i,j,k))
        enddo
!       do j=2,jmax-1
!        write(66,*)j+par%tjmax(1),real(tkea_w(i,j,k)),real(epsa_w(i,j,k)),real(tken_w(i,j,k))
!       enddo
        close(66)
       endif                        !RANK>
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      enddo ! looprank_
      endif            !>>>

      if(flag==0) then !>>>
      i=nint(real(iglb)/2.)-par%timax(1)
      j=nint(real(jglb)/2.)-par%tjmax(1)
      if(i>1.and.i<imax.and.&
         j>1.and.j<jmax) then !pmx> !19-08-21>
      write(texte30,'(a,i0)')'tem_init_rank',par%rank
      open(unit=66,file=trim(texte30))
      do k=2,kmax
       write(66,*)depth_t(i,j,k),tkea_w(i,j,k),epsa_w(i,j,k)
      enddo
      close(66)
      endif                   !pmx> !19-08-21>
      endif            !>>>

      call graph_out
! un tour en key
      elapsedtime_now=0
      key=1000
!     key=5000
!     key=300
!     key=4*250
!     key=4*150
!     do 2000 iteration3d=0,key/2      ! quart de un tour
!     do 2000 iteration3d=1-1,key-1           ! quart de un tour
!     do 2000 iteration3d=0,3*key/4     ! 3/4 de tour
!     do 2000 iteration3d=0,key         ! un tour
!     do 2000 iteration3d=1,nint(20./current_number_)

!     do 2000 iteration3d=1,nint(100./dti_fw/0.001)

      if(flag==1) then !11111111111>
! Cas unidirectionnel Oi
      omega_w=0.
      veldxdz_v=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
!       veldydz_u(i,j,k,1)=+(2.*dx_u(i,j)/dti_fw)*dy_u(i,j)*dz_u(i,j,k,1)      
        veldydz_u(i,j,k,1)=-(2.  *dx_u(i,j)/dti_fw)*dy_u(i,j)*dz_u(i,j,k,1)      
      enddo       ; enddo       ; enddo
! Cas unidirectionnel Oj
!     omega_w=0.
!     veldydz_u=0.
!     do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!       veldxdz_v(i,j,k,1)=+(2.0*dy_v(i,j)/dti_fw)*dx_v(i,j)*dz_v(i,j,k,1)      
!       veldxdz_v(i,j,k,1)=-(2.0*dy_v(i,j)/dti_fw)*dx_v(i,j)*dz_v(i,j,k,1)      
!     enddo       ; enddo       ; enddo
      endif            !11111111111>

!     loopmax_=100000
!     loopmax_=300000

      do 2000 iteration3d=0,iteration3d_max

      elapsedtime_now=elapsedtime_now+dti_fw

      if(par%rank==0)then
!         if(mod(iteration3d,10)==0)write(6,*)iteration3d,key
          if(mod(iteration3d,10)==0)write(6,*)iteration3d
!         write(6,*)iteration3d
      endif


! Termes annulEs dans l'equation de tke de gaspar
            kh_w=0.
       wetmask_t=1.
!      tfilterfb=0.
!      tfilterlf=0.
       upwindriver_t=1.
       omega_evaprec_w=0.


!___________________________________________________
!..............................................................

! POUR schema k-epsilon
!     tken_w=tkea_w
!     call turbulence_adv_tkea
!     call obc_turbulence_tkea

!     epsn_w=epsa_w
!     call turbulence_adv_epsa
!     call obc_turbulence_epsa

! POUR schema tke gaspar
      call turbulence_adv_tken
      call obc_turbulence_tken

      if(mod(iteration3d,10000)==0.or. &
                    iteration3d==iteration3d_max)then !ecrire fichier>

       if(flag==1) then !111>
        k=kmax
        i=imax/2 
        j=jmax/2 
        do looprank_=0,nbdom-1
         if(par%rank==looprank_) then !RANK>
          if(par%rank==0) then
           open(unit=66,file='tem_fin')
          else
           open(unit=66,file='tem_fin',position='append')
          endif
          do i=2,imax-1
           write(66,*)i+par%timax(1),real(tkea_w(i,j,k)),real(epsa_w(i,j,k)),real(tken_w(i,j,k))
          enddo
!         do j=2,jmax-1
!          write(66,*)j+par%tjmax(1),real(tkea_w(i,j,k)),real(epsa_w(i,j,k)),real(tken_w(i,j,k))
!         enddo
          close(66)
         endif                        !RANK>
#ifdef parallele
        call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
         enddo ! looprank_
       endif            !111> 

      endif                                    !ecrire fichier>

!      if(par%rank==0)write(68,*)iteration3d*dti_fw/86400 &
!       ,tem_t(imax/2,jmax/2,9,1) &
!      ,temf_t(imax/2,jmax/2,9)/ratio_negdif_ver

 2000 continue

!     write(67,*)current_number_,x1,x2

      enddo ! fin de boucle sur loop_

!0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1)) 

      sum1=0.
      sum3=0.
      sum5=0.
      do k=2,kmax ; do j=1,jmax ; do i=1,imax
         sum1=sum1+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1))*tkea_w(i,j,k)
         sum3=sum3+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1))*epsa_w(i,j,k)
         sum5=sum5+mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1))*tken_w(i,j,k)
      enddo       ; enddo       ; enddo
      call mpi_allreduce(sum1,sum2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum4,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum5,sum6,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)'Somme tkea finale  ',sum2
      if(par%rank==0)write(6,*)'Somme epsa finale  ',sum4
      if(par%rank==0)write(6,*)'Somme tken finale  ',sum6

      k=kmax

      call graph_out
      write(6,*)'dans sbr test_advection_keps'
      write(6,*)'-----------------------'
      stop 'test_advection_keps'


      end subroutine test_advection_keps
#endif
