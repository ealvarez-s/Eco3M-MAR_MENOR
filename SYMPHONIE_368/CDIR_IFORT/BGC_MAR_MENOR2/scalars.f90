










      subroutine scalars
!______________________________________________________________________
! SYMPHONIE ocean model 
! release m�v�m - last update: 21-10-18
!______________________________________________________________________
      use module_principal                                        !27/01/03
      use module_parallele                                        !10-05-14
      use module_my_outputs
      implicit none


!...............................................................................
!    _________                    .__                  .__             ! m�v�m
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
!                   "QUICK_CA" impl�ment� pour modele par Guillaume Reffray
!                   du LSEET.
!         27/01/03: pass� dans  version 2003
!         03/11/03: bilan d'energie sur terme de melange vertical....
!         29/05/04: possibilit� de faire CALL OBC_SCAL avant la turbulence
!                   & nouveaux arguments dans OBC_SCAL
!         08/07/04: routines non compatibles comment�es
!         16/06/06: ajout possibilite de convection automatique sur suggestion
!                   de Claude.
!         17/10/06: Nouvelle advection pour T et S
!         20/12/06: schema d'advection choisi par le passage en argument de
!                   IADVEC_TS
!         21/12/06: Les conditions aux limites sur les rivieres sont calcul�es
!                   avant l'appel � l'advection
!         10/12/07: ajout d'une routine calculant le niveau vertical correspondant
!                   � la base de la couche de m�lange de surface
!         14/01/08: Ajout controle pression de fond dans la couche eponge
!         21/01/08: des parametres ajout�s dans notebook_spongelayer permette
!                   de choisir les conditions aux limites et les eponges
!         31/01/08: eponges calculees dans routine specifique -> appel �
!                   APPLY_SPONGE_LAYER. Deux cas possibles: 1=eponge classique
!                   2=eponge � densit� inchang�e
!         26/02/08: Dans cas 2. Avant on modifiait en dur les tableaux TOBC et SOBC
!                   (� densite inchangee) puis on relaxait vers ces derniers.
!                   Mais cela occasionnait une interference avec obc_scal qui utilise
!                   ces m�mes tableaux. Maintenant le calcul est le m�me mais temobc
!                   et salobc sont laiss�s intacts.
!         07/04/08: ajout RELAXTYPE_TS pour distinguer la relaxation des
!                   traceur en mode classique du mode "densit� inchang�e"
!         13/05/08: Deplacement de la condition au limite sur TEM et SAL juste
!                   avant l'advection. Consequence: interdit les obc sur THZ(2)...
!         05/11/08: Appel � K_SURFLAYER avant l'advection pour eventuellement
!                   utiliser ksl_z dans le schema d'advection.
!         01-06-09: Appel � obc_scal: passage en arg de obctype_ts pour faire
!                   un appel avec un arg different dans initial_with_obc
! 2009.3  02-10-09: dti_lp remplace dti
!                   utilisation des nouveaux facteurs d'echelle verticale
!         10-10-09: tem_c et sal_c remplacent thz_c et shz_c
!         15-10-09: manip sur i1d
! 2010.3  15-01-10  - subroutine advection_ts renomm�e advection_scal
!                   - subroutine temperature renomm�e vertmix_tem
!                   - subroutine salinity renomm�e vertmix_sal
! 2010.8  06-05-10  Pas d'argument pass� dans advection_scal.F90
! 2010.11 16-07-10  temobc & salobc renomm�s temobc & salobc
!         23-07-10  equation d'etat donne densit� potentielle
! 2010.14 16-12-10  Equation d'etat appel� au cas par cas depuis calcul PGF
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 22-04-11  Argument pass� dans botprescont
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
!         21-10-18  ajout d'un cas n� 5 pour le rappel vers T,S
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

! Niveau vertical correspondant � la base de la couche de m�lange      !05/11/08
      call scalars_k_surflayer

      if(flag3d.eq.1) then ! �������������������>

            call obc_river(1,0,1) ! C.L. T & S t-dt, t !08-05-14

      endif             ! �������������������>                         !15-10-09

! Lateral Open Boundary Conditions:
            call obc_scal(obctype_ts)                                  !01-06-09
! mpi exchange with mpi neighbours:
            call obc_scal_mpi_now


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

! Controle de la pression de fond dans la couche eponge:
      if(relax_bpc>0.)call botprescont('eulerian  ')                !22-04-11

! Calculer le rappel dans la zone eponge:                              !31/01/08
      if(relax_ts>0.)call scalars_spongelayer                            !07/04/08
! BUFFER ZONE MER MARMARA
      call scalars_spongelayer_marmara


!     if(flag_merged_levels==1)call scalars_crushed_cells     !17-10-17

!     call scalars_ts_minmax(2)                                        !11-01-16
!     call scalars_tsobc_minmax(2)
!     call scalars_global_conservation                                  !27-03-16
!     call my_outputs_tem_sum

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

      if(convect_yn==3)return !01-04-14

! KSL_Z(I,J) est le niveau vertical correspondant � la base de la couche
! de m�lange:

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

      end subroutine scalars_k_surflayer                                            !10/12/07
!-----------------------------------------------------------------------------
      subroutine scalars_spongelayer
      use module_principal
      use module_parallele
      implicit none
      real*4 factor_relax_,factor_lwf_
      integer :: loop_        &
                ,moduliter_=1
!     double precision,dimension(:,:),allocatable ::  iaveraged_in,iaveraged_out

!...................................................................
! Full grid relaxation of T and S toward OGCM field - algo "Low Frequency"
      if(relaxtype_ts==2) then
!...................................................................

! Les valeurs factor_relax_=0.6 et factor_lwf_1.37 permettent au schema
! de reduire de 90% l'ecart � la reference en un temps equivalent �
! celui que prendrait l'algo classique relaxtype_ts==2
      factor_relax_=0.6 ; factor_lwf_=1.37
!     factor_relax_=1.  ; factor_lwf_=1.

! Relax:
      x1=dti_lp*relax_ts*factor_relax_*factor_lwf_
!     x2=x1*    timeweightobc(trc_id) ; x0=x1*(1.-timeweightobc(trc_id))
      x2=timeweightobc(trc_id) ; x0=1.-timeweightobc(trc_id) !30-03-16
! Low freq:
      const2=dti_lp*relax_ts*factor_lwf_
      const1=1.-const2

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Low frequency temperature:
!      temlwf_t(i,j,k)=const1*temlwf_t(i,j,k)+const2*tem_t(i,j,k,0)
       temlwf_t(i,j,k)=const1*temlwf_t(i,j,k)  &

                        +const2*( & !ooo>

          tem_t(i,j,k,0)-(x2*temobc_t(i,j,k,2)+x0*temobc_t(i,j,k,0) ) &

                                )   !ooo>

! Nudging process based on the difference between the Low Frequcy field and the reference field:
       tem_t(i,j,k,2)=tem_t(i,j,k,2)              &
!                 +temobc_t(i,j,k,2)*x2           & !30-03-16
!                 +temobc_t(i,j,k,0)*x0           &
                  -temlwf_t(i,j,k)*x1

! Low frequency salinity:
!      sallwf_t(i,j,k)=const1*sallwf_t(i,j,k)+const2*sal_t(i,j,k,0)
       sallwf_t(i,j,k)=const1*sallwf_t(i,j,k)  &

                        +const2*( & !ooo>

          sal_t(i,j,k,0)-(x2*salobc_t(i,j,k,2)+x0*salobc_t(i,j,k,0) ) &

                                )   !ooo>

! Nudging process based on the difference between the Low Frequcy field and the reference field:
       sal_t(i,j,k,2)=sal_t(i,j,k,2)              &
!                 +salobc_t(i,j,k,2)*x2           &
!                 +salobc_t(i,j,k,0)*x0           &
                  -sallwf_t(i,j,k)*x1


      enddo
      enddo
      enddo

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
! UPWIND
      if(relaxtype_ts==1) then                                       !07/04/08
!...................................................................


        loop1=1 ! Frontiere i=1 !18-05-18

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dx_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_u(i+1,j,k,2)    & !25-05-18
                              +(1.-timeweightobc(vel_id))*velobc_u(i+1,j,k,0),0.)

           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_u(i+1,j,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_u(i+1,j,k,0),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=2 ! Frontiere i=imax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dx_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_u(i+1,j,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_u(i+1,j,k,0),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_u(i+1,j,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_u(i+1,j,k,0),0.) 

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=3 ! Frontiere j=1

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dy_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_v(i,j+1,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_v(i,j+1,k,0),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(     timeweightobc(vel_id) *velobc_v(i,j+1,k,2)    &
                              +(1.-timeweightobc(vel_id))*velobc_v(i,j+1,k,0),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=4 ! Frontiere j=jmax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=dti_lp/dy_t(i,j)*sponge_t(i,j,1) 

          do k=kmin_w(i,j),kmax                ! debut boucle k

           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   &
                         *max(    -timeweightobc(vel_id) *velobc_v(i,j  ,k,2)    &
                              -(1.-timeweightobc(vel_id))*velobc_v(i,j  ,k,0),0.)
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))   &
                         *max(    -timeweightobc(vel_id) *velobc_v(i,j  ,k,2)    &
                              -(1.-timeweightobc(vel_id))*velobc_v(i,j  ,k,0),0.)

          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j


!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! UPWIND
      return
      endif
!...................................................................

!...................................................................
! Relaxation of T and S toward OGCM field within a lateral nudging layer
! at unchanged potential density (provided that a linear EOS is used)
      if(relaxtype_ts==0) then                                       !07/04/08
!...................................................................

      const1=-rho*alp_t
      const2= rho*alp_s
      const4=dti_lp*relax_ts        ! 22-04-11
      do loop1=1,4                          ! boucle sur 4 fronti�res
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

        const4=dti_lp*relax_ts       ! 22-04-11

        loop1=1 ! Frontiere i=1 !18-05-18

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=const4*sponge_t(i,j,1)           !22-04-11
          do k=kmin_w(i,j),kmax                ! debut boucle k
           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))    
          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=2 ! Frontiere i=imax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=const4*sponge_t(i,j,1)           !22-04-11
          do k=kmin_w(i,j),kmax                ! debut boucle k
           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))    
          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=3 ! Frontiere j=1

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=const4*sponge_t(i,j,1)           !22-04-11
          do k=kmin_w(i,j),kmax                ! debut boucle k
           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))    
          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j

        loop1=4 ! Frontiere j=jmax

        do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
        do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

          x0=const4*sponge_t(i,j,1)           !22-04-11
          do k=kmin_w(i,j),kmax                ! debut boucle k
           tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0*(temobc_t(i,j,k,1)-tem_t(i,j,k,2))   
           sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0*(salobc_t(i,j,k,1)-sal_t(i,j,k,2))    
          enddo                                ! fin boucle k

        enddo                                ! fin boucle i
        enddo                                ! fin boucle j


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

      subroutine scalar_check_mpi_conservation
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none

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
        write(10+par%rank,*)'coordonn�es ',i,j0-1   ! jmax-1
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
        write(10+par%rank,*)'coordonn�es ',i,2
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
        write(10+par%rank,*)'coordonn�es ',i0-1,j ! imax-1,j
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
        write(10+par%rank,*)'coordonn�es',2,j
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

      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
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
        write(10+par%rank,*)'coordonn�es ',i,j0-1   ! jmax-1
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
        write(10+par%rank,*)'coordonn�es ',i,2
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
        write(10+par%rank,*)'coordonn�es ',i0-1,j ! imax-1,j
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
        write(10+par%rank,*)'coordonn�es',2,j
        write(10+par%rank,*)'sal_t ',anyv3d(i2,j,k,2),sal_t(i2,j,k,1)
        write(10+par%rank,*)'delta ',anyv3d(i2,j,k,2)-sal_t(i2,j,k,1)
        write(10+par%rank,*)'salobc local ',salobc_t(i2,j,k,1)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! j

      enddo ! kmax

      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
      if(ksecu==1) then
        write(6,'(a,i0,a)')   &
        ' mpi conservation error for sal_t Details in fort.' &
         ,10+par%rank,' available in RDIR/CONFIG directory'
        stop ' STOP in routine scalar_check_mpi_conservation!'
      endif

      end subroutine scalar_check_mpi_conservation

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

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
           write(6,*)'MIN T val',x0,tem_t(i1,j1,k1,t_)
           write(6,*)'MIN T rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN T mask i j k',par%rank,mask_t(i1,j1,k1)
           write(6,*)'MIN T glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN T dz_t',dz_t(i1,j1,k1,t_)
           write(6,*)'MIN T h_w',h_w(i1,j1)
      endif           !>>>

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX T=',x0
      if(x2==x0) then !>>>
           write(6,*)'MAX T val',x0,tem_t(i2,j2,k2,t_)
           write(6,*)'MAX T rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX T mask i j k',par%rank,mask_t(i2,j2,k2)
           write(6,*)'MAX T glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX T dz_t',dz_t(i2,j2,k2,t_)
           write(6,*)'MAX T h_w',h_w(i2,j2)
      endif           !>>>

      call mpi_allreduce(x3,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN S=',x0
      if(x3==x0) then !>>>
           write(6,*)'MIN S val',x0,sal_t(i3,j3,k3,t_)
           write(6,*)'MIN S rank i j k',par%rank,i3,j3,k3
           write(6,*)'MIN S mask i j k',par%rank,mask_t(i3,j3,k3)
           write(6,*)'MIN S glob i j k',i3+par%timax(1),j3+par%tjmax(1),k3
           write(6,*)'MIN S dz_t',dz_t(i3,j3,k3,t_)
           write(6,*)'MIN S h_w',h_w(i3,j3)
      endif           !>>>

      call mpi_allreduce(x4,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX S=',x0
      if(x4==x0) then !>>>
           write(6,*)'MAX S val',x0,sal_t(i4,j4,k4,t_)
           write(6,*)'MAX S rank i j k',par%rank,i4,j4,k4
           write(6,*)'MAX S mask i j k',par%rank,mask_t(i4,j4,k4)
           write(6,*)'MAX S glob i j k',i4+par%timax(1),j4+par%tjmax(1),k4
           write(6,*)'MAX S dz_t',dz_t(i4,j4,k4,t_)
           write(6,*)'MAX S h_w',h_w(i4,j4)
      endif           !>>>

      end subroutine scalars_ts_minmax
!-----------------------------------------------------------------------------

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

!-----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!#ifdef bidon
      subroutine test_advection_vel_u
      use module_principal
      use module_parallele
      implicit none
      real :: current_number_=0.1
      integer loop_


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
         dist0=sqrt( (real(i1)-x1)**2 + (real(j1)-x2)**2 )
         dist1=max(10.-dist0,zero)/10.
         vel_u(i,j,k,:)=dist1+x3*cos(pi*i1)*cos(pi*j1)*dist1
      endif

      if(flag==0) then ! 2D vertical OXZ
         dist0=sqrt( (real(i1)-x1)**2 + (real(k)-0.25*real(kmax))**2 )
         dist1=max(10.-dist0,zero)/10.
         vel_u(i,j,k,:)=dist1
      endif

        if(i<iglb+2)vel_v(i,j,k,:)=vel_u(i,j,k,:)


! Initialiser un front lin�aire
!        x1=real(i1)/real(iglb)+real(j1)/real(jglb)
!        vel_u(i,j,k,:)=x1
!        x1=  real(j1-i1+iglb) / real(iglb+jglb) - 0.5
!        vel_u(i,j,k,:)=tanh(20.*x1)


       endif                        !--------------->

      enddo
      enddo
      enddo


! POUR VERIFIER QU'UN CHAMPS CONSTANt RESTE INCHANGE
!     vel_u=1.
!     vel_v=1.


      if(flag==1) then !hor
      k=kmax
      i=nint(real(iglb)/4.)
      open(unit=66,file='vel_init')
      do j=1,jmax
       write(66,*)j,vel_u(i,j,k,2),vel_v(i,j,k,2)
      enddo
      close(66)
      endif            !hor
      if(flag==0) then !000> vertical
       i=nint(real(iglb)/4.) ; j=jmax/2
       open(unit=66,file='vel_init')
       do k=1,kmax
        write(66,*)k,vel_u(i,j,k,2),vel_v(i,j,k,2)
       enddo
       close(66)
      endif            !000> vertical


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
      do 2000 iteration3d=1,nint(20./current_number_)
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
         veldydz_u(i,j,k,1)=x3        *dy_u(i,j)*dz_u(i,j,k,1)      
!        veldydz_u(i,j,k,1)=0.
       enddo       ; enddo       ; enddo
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!        veldxdz_v(i,j,k,1)=x3*cos(x2)*dx_v(i,j)*dz_v(i,j,k,1)   
         veldxdz_v(i,j,k,1)=x3        *dx_v(i,j)*dz_v(i,j,k,1)   
!        veldxdz_v(i,j,k,1)=0. 
       enddo       ; enddo       ; enddo
      endif            !000>

      if(flag==0) then !111>
         veldydz_u=0.
         veldxdz_v=0.
         do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
          omega_w(i,j,k,1)=x4
         enddo       ; enddo       ; enddo
      endif            !111>



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

      vel_u(:,:,:,-1)=vel_u(:,:,:,0)
      vel_u(:,:,:,0)= vel_u(:,:,:,1)
      vel_u(:,:,:,1)= vel_u(:,:,:,2)

      vel_v(:,:,:,-1)=vel_v(:,:,:,0)
      vel_v(:,:,:,0)= vel_v(:,:,:,1)
      vel_v(:,:,:,1)= vel_v(:,:,:,2)

 2000 continue

      if(flag==1) then !111> horizontal
       k=kmax
       i=nint(real(iglb)/4.)
       write(texte90,'(a,i0)')'vel_',loop_
       open(unit=66,file=texte90)
       x1=0.
       do j=1,jmax
        write(66,*)j,vel_u(i,j,k,2),vel_v(i,j,k,2)
        x1=max(x1,vel_u(i,j,k,2))
       enddo
       close(66)
       write(67,*)current_number_,x1,  anyv3d(imax/2,jmax/2,kmax,2)
      endif             !111>

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

      call graph_out
      write(6,*)'dans sbr test_advection_bio'
      write(6,*)'-----------------------'
      stop 'test_advection_vel_u'


      end subroutine test_advection_vel_u

!#endif
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!#ifdef bidon
      subroutine test_advection_temsal
      use module_principal
      use module_parallele
      implicit none
      real :: current_number_=0.1
      integer loop_

! UNE CHOSE IMPORTANTE POUR COMPRENDRE CERTAINS TEST.
! Nous sommes en Leap-Frog. Le nombre de courant s'exprime dx (ou dz) divise par dti_lp soit 2 fois dt
! Autrement dit si le nombre de courant est 1 et que nous faisons 10 iterations, la forme se deplace
! de 10*dt*dx/(2dt) soit 5dx, donc 5 mailles (et non pas 10 comme intuitivement)
! Attention aux initialisation un peu simple (T(0)=T(1)) et au leap-frog qui supervise 2 trajectoires,
! une paire, une impaire. iParfois il faut comparer les paires entre elles, les impaires entre elles.

      do loop_=1,3 ! 25,25
       if(loop_==3) then
         current_number_=0.999999*1
         key=30    
       endif
       if(loop_==2) then
         current_number_=0.999999*0.5
         key=30/0.5
       endif
       if(loop_==1) then
         current_number_=0.999999*0.1
         key=30/0.1
       endif

!     current_number_=0.5

!     flag=1 ! Cas 2D horizontal
      flag=0 ! Cas 2D vertical OXZ (eventuellement reduire jglb pour gagner cpu)
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

         if(flag==1) then !>>>
! Tester l'advection horizontale:
! Initialiser un tourbillon:
         dist0=sqrt( (real(i1)-x1)**2 + (real(j1)-x2)**2 )
         dist1=max(15.-dist0,zero)/15.
         tem_t(i,j,k,:)=1.+dist1+x3*cos(pi*i1)*cos(pi*j1)*dist1
! Initialiser un tourbillon:
         dist0=sqrt( (real(i1)-x1)**2 + (real(j1)-x2-10.)**2 )
         dist1=max(15.-dist0,zero)/15.
         sal_t(i,j,k,:)=1.+dist1+x3*cos(pi*i1)*cos(pi*j1)*dist1
         endif            !>>>


         if(flag==0) then !>>>
! Tester l'advection verticale:
! Une gaussienne:
         dist1=    min(max( real(k-10) , 0.),20.)/20.
         tem_t(i,j,k,:)=sin(pi*dist1)**2
         sal_t(i,j,k,:)=sin(pi*dist1)**2
! Un front:
         if(k>kmax/2) then
          tem_t(i,j,k,:)=0.
          sal_t(i,j,k,:)=0.
         else
          tem_t(i,j,k,:)=1.
          sal_t(i,j,k,:)=1.
         endif

         endif            !>>>

       enddo
       enddo

       do i1=0,iglb+1
       do j1=0,jglb+1
        i=i1-par%timax(1)
        j=j1-par%tjmax(1)
         temobc_t(i,j,k,:)=tem_t(i,j,k,1)
         salobc_t(i,j,k,:)=sal_t(i,j,k,1)
       enddo
       enddo

      enddo

      if(flag==1) then !>>>
      k=kmax
      i=nint(real(iglb)/4.)
      open(unit=66,file='tem_init')
      do j=1,jmax
       write(66,*)j,tem_t(i,j,k,1),sal_t(i,j,k,1)
      enddo
      close(66)
      endif            !>>>
      if(flag==0) then !>>>
      i=nint(real(iglb)/4.)
      j=nint(real(jglb)/2.)
      open(unit=66,file='tem_init')
      do k=1,kmax
       write(66,*)k,tem_t(i,j,k,1),sal_t(i,j,k,1)
      enddo
      close(66)
      endif            !>>>

      call graph_out

! un tour en key
      elapsedtime_now=0
!     key=1000
!     key=5000
!     key=300
!     key=4*250
!     key=4*150
!     do 2000 iteration3d=0,key/2      ! quart de un tour
      do 2000 iteration3d=1-1,key-1           ! quart de un tour
!     do 2000 iteration3d=0,3*key/4     ! 3/4 de tour
!     do 2000 iteration3d=0,key         ! un tour
!     do 2000 iteration3d=1,nint(20./current_number_)

      elapsedtime_now=elapsedtime_now+dti_fw

      if(par%rank==0)then
!         if(mod(iteration3d,10)==0)write(6,*)iteration3d,key
          write(6,*)iteration3d
      endif


! Ces lignes servent A Activer l'aiguillage "anti-boules-d'eau-dense"
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
         rhp_t(i,j,k)=tem_t(i,j,k,1)
        enddo ; enddo ; enddo

!..............................................................
!___________________________________________________
! Module du courant horizontal
      if(flag==1) then !1111111111111111111111111>
       x3=current_number_*dxb/dti_lp ! Nb de courant * DX / DT
!      write(6,*)'x3=',x3
! direction du courant:
!     x2=2.*pi*real(iteration3d)/real(key)  ! rotation
!     x2=0.25*pi
!     x2=0.5*pi
      x2=0.
!  wsed(:)=-0.1*h_w(i,j)*dsig_t(i,j,k)/dti_fw
! Vitesses advectantes:
!     omega_w=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        veldydz_u(i,j,k,1)=x3*sin(x2)*dy_u(i,j)*dz_u(i,j,k,1)      
!       veldydz_u(i,j,k,1)=x3        *dy_u(i,j)*dz_u(i,j,k,1)      
      enddo       ; enddo       ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        veldxdz_v(i,j,k,1)=x3*cos(x2)*dx_v(i,j)*dz_v(i,j,k,1)   
!       veldxdz_v(i,j,k,1)=x3        *dx_v(i,j)*dz_v(i,j,k,1)   
      enddo       ; enddo       ; enddo
      dz_t(:,:,:,2)=dz_t(:,:,:,0)
      call omega
!     write(6,*)'w',omega_w(30,30,kmax+1,1)
      do k=1,kmax 
      do j=1,jmax
      do i=1,imax
       dz_t(i,j,k,2)=dz_t(i,j,k,0)+dti_lp*omega_w(i,j,kmax+1,1)/real(kmax)
      enddo
      enddo
      enddo
      call omega
!     write(6,*)'w',omega_w(30,30,kmax+1,1),dz_t(30,30,kmax,2)
      endif            !1111111111111111111111111>

      if(flag==0) then !00000000000>
       veldydz_u=0.
       veldxdz_v=0.
       do k=2,kmax   ; do j=1,jmax ; do i=1,imax
        omega_w(i,j,k,1)= current_number_*h_w(i,j)*dsig_t(i,j,k)/dti_lp
       enddo       ; enddo       ; enddo
       omega_w(:,:,1     ,1)=0.
       omega_w(:,:,kmax+1,1)=0.
      endif            !00000000000>

! Termes annulEs dans l'equation de tke de gaspar
            kh_w=0.
       wetmask_t=1.
       tfilterfb=0.
       tfilterlf=0.
       upwindriver_t=1.


!___________________________________________________
!..............................................................

!     write(66,*)'-----------'
!     write(67,*)'-----------'
         call vertmix_matrix
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

! Verifi de conservation:
       sum1=0.
       sum2=0.
       do k=1,kmax
       do j=1,jmax ; do i=1,imax
       sum1=sum1+tem_t(i,j,k,1)*dxdy_t(i,j)*dz_t(i,j,k,1)
       sum2=sum2+sal_t(i,j,k,1)*dxdy_t(i,j)*dz_t(i,j,k,1)
!        if(tem_t(i,j,k,2)/=sal_t(i,j,k,2)) then
!         write(6,*)i,j,k,tem_t(i,j,k,2),sal_t(i,j,k,2)
!         stop 'tem/=sal'
!        endif
       enddo ; enddo
       enddo
!      write(6,*)'sums=',sum1,sum2

 2000 continue

      if(flag==1) then !>>>
       k=kmax
       i=nint(real(iglb)/4.)
       write(texte90,'(a,i0)')'tem_',loop_
       open(unit=66,file=texte90)
       x1=0.
       x2=0.
       do j=1,jmax
        write(66,*)j,tem_t(i,j,k,1),sal_t(i,j,k,1)
        x1=max(x1,tem_t(i,j,k,1))
        x2=max(x2,sal_t(i,j,k,1))
       enddo
       close(66)
      endif

      if(flag==0) then !>>>
       i=nint(real(iglb)/4.)
       j=nint(real(jglb)/2.)
       write(texte90,'(a,i0)')'tem_',loop_
       open(unit=66,file=texte90)
       x1=0.
       x2=0.
       do k=1,kmax
        write(66,*)k,tem_t(i,j,k,1),sal_t(i,j,k,1)
        x1=max(x1,tem_t(i,j,k,1))
        x2=max(x2,sal_t(i,j,k,1))
       enddo
       close(66)
      endif

!     write(67,*)current_number_,x1,x2

      enddo ! fin de boucle sur loop_

      call graph_out
      write(6,*)'dans sbr test_advection_temsal'
      write(6,*)'-----------------------'
      stop 'test_advection_temsal'


      end subroutine test_advection_temsal
!#endif


      subroutine scalars_spongelayer_marmara
      use module_principal ; use module_parallele
      implicit none
      real*4 :: c1_,c2_,c3_
      integer :: istr_,jstr_,iend_,jend_

      if(iglb/=1120.or.jglb/=865) &
      stop 'scalars_spongelayer_marmara: iglb/=1120.or.jglb/=865'

      istr_=max(1056-par%timax(1),1)
      iend_=min(1099-par%timax(1),imax)
      jstr_=max( 247-par%tjmax(1),1)
      jend_=min( 299-par%tjmax(1),jmax)

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

!-----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
