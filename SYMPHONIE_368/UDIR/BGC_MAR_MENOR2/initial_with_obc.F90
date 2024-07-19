      subroutine initial_with_obc(ki1)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 310 - last update: 26-10-21
!______________________________________________________________________

      use module_principal ; use module_parallele
      implicit none
      integer ki1,loop_
#ifdef synopsis
       subroutinetitle='initial_with_obc'
       subroutinedescription= &
       'Initial state provided from interpolated OGCM (often' &
       //' MERCATOR) fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         22/07/01: passage à la coordonnée sigma généralisée
!         10/08/01: ajout initialisation de la densité
!         27/10/01: KI1 passé en argument
!         30/10/01: Initialisation de T & S aux sources des fleuves
!         31/10/01: prise en compte des masques:
!         14/03/02: amenagement pour model_ inverse. Initialisation plus complete
!         26/08/02: VBAROBC remplacé par VMEAOBC
!                   + inversion des boucles 55 & 57
!                   + modifs mineures
!         04/11/03: passage d'un argument dans CALL DENSITY
!         23/03/04: amelioration de la conservation de la masse à l'instant
!                   initial
!         14/11/04: l'initialisation de la densite est maintenant faite dans
!                   initial_state_eq.F
!         27/07/05: remaniement des boucles sur T et S pour s'assurer que les
!                   conditions aux limites uptream sont bien initialisées...
!         15/09/05: initialisation tableaux ztarefobc et cie...
!         11/08/06: Appel à vmto_vmhz supprimé
!         02/10/06: initialisation des tableaux de passage des modes externes
!                   jumaux
!         17/04/07: passage coordonnees curvilignes
!         07/06/07: initialisation de vmean(2)
!         24/02/08: modif reset suite à redimensionnement de VMEAN_
!         13/05/08: sauvegarde de tem_c et sal_c dans shz_c et dans thz_c
!                   pour advection
!         31-05-09  Parallelisation
!         01-06-09  Parallelisation: appel à obc_scal
! 2009.3  30-09-09  utilisation des nouveaux facteurs d'echelle verticale
!         11-10-09  tem_c et sal_x remplacent thz_z et shz_z
!         14-10-09  Pour le schema d'advection de t et s il faut initialiser
!                   tem_c(before) sur une grille aggrandie
!         19-10-09  mettre toute la grille à jour
! 2010.7  16-02-10  modif dimension veldxdz
! 2010.8  10-03-10  ajout echeance "-1"
!         22-03-10  ajout hssh_w
!         03-05-10  nouvelles dimensions pour veldxdz et veldydz
! 2010.10 13-06-10  suppression ssh_ext_w
! 2010.11 16-07-10  temobc & salobc renommés temobc & salobc
! 2010.12 03-09-10  call z_levels(0) remplace call z_levels(1)
! 2010.23 24-05-11  initialiser dz_t(i,j,k,-1)
! S26     27-01-13  advection qdm o4 entraine des echanges supplementaires
!         17-07-13  Etat initial hssh positive
!         09-10-13  Possibilite d'un nouvel algo de nudging
!         12-05-14  call obc_scal_mpi_now
!         02-07-14  call obc_int_mpi4d !02-07-14
!                   call obc_ext_mpi3d !02-07-14
!         30-07-14  call obc_scal avant call des echanges
!         11-11-14  adapter l'initialisation a la 3eme dim de hssh_w
!         28-11-14  suppression ssh_int_u ssh_int_v
!         11-04-15  boucles ssh_w
!         17-12-15  boucles f90
!         03-02-16  call couple_modes_currentnumber pour reset veldtodx...
!         02-04-16  - allocation de temlwf et sallwf conditionnEe A la valeur
!                   de relax_ts
!                   - Desormais anomalies temlwf et sallwf sont initialisEs A zEro
!         04-09-16  initialisation ssh en harmonie avec z_thickness.F90
! v310    26-10-21  Possibilite d'appel A vertmix_merged_levels_obc (actuellement 
!                   sans objet si T et S ne sont pas melangEs dans la couche merged)
!...............................................................................
!    _________                    .__                  .__             ! m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

! Ce programme initialise le modèle avec la famille de
! tableaux OBC dernier indice 1.

!     call vertmix_merged_levels_obc(ki1) !26-10-21

!*****************************************************************
! Initialisation de l'élévation de surface:
!*****************************************************************
      do j=1,jmax
      do i=1,imax
!      hssh_w(i,j,2)=                                                &
!      max(sshobc_w(i,j,ki1)*mask_t(i,j,kmax+1)+h_w(i,j),wetdry_cst2) ! 17-07-13
               ssh_w(i,j,2)                                          & !27-11-14
!      =max(sshobc_w(i,j,ki1)*mask_t(i,j,kmax+1),wetdry_cst2-h_w(i,j))
       =max(sshobc_w(i,j,ki1)*mask_t(i,j,kmax+1),0.001*wetdry_cst3-h_w(i,j)) !04-09-16
! MOD Thom...
       !=max(sshobc_w(i,j,ki1)*mask_t(i,j,kmax+1),wetdry_cst3-h_w(i,j)) !02-03-24
      enddo
      enddo

      call obc_ssh ! c.l. en i=0, =imax+1, j=0, =jmax+1 pour 1/2 somme en x et y

      lb3=lbound(ssh_w) ; ub3=ubound(ssh_w) !11-04-15
      do k=lb3(3),1
      do j=lb3(2),ub3(2)
      do i=lb3(1),ub3(1)
        ssh_w(i,j,k)=ssh_w(i,j,2)
      enddo ; enddo ; enddo

      lb3=lbound(ssh_int_w) 
      do k=lb3(3),2
       do j=0,jmax+1 ; do i=0,imax+1
        ssh_int_w(i,j,k)=ssh_w(i,j,2)
       enddo ; enddo
      enddo

!     lb3=lbound(ssh_int_u)
!     do k=lb3(3),2
!      do j=0,jmax+1 ; do i=1,imax+1
!       ssh_int_u(i,j,k)=max(0.5*(ssh_int_w(i,j,2)+ssh_int_w(i-1,j,2)) &
!                            ,0.001*wetdry_cst3-h_u(i,j) )
!      enddo  ; enddo
!     enddo

!     lb3=lbound(ssh_int_v)
!     do k=lb3(3),2
!      do j=1,jmax+1 ; do i=0,imax+1
!       ssh_int_v(i,j,k)=max(0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j-1,2)) &
!                            ,0.001*wetdry_cst3-h_v(i,j) )
!      enddo ; enddo
!     enddo

!......Mise à jour de la grille:
      call z_thickness(0,2)                                             !31-05_09
      call cellbox_thickness(0,-1,2)                                    !24-05-11
      call z_levels(0)                                                  !03-09-10

!*****************************************************************
! Initialisation de T, S et V:
!*****************************************************************

       lb4=lbound(tem_t) ; ub4=ubound(tem_t) !17-12-15
       do j=lb4(2),ub4(2)
       do i=lb4(1),ub4(1)
        i1=min(max(i,1),imax)
        j1=min(max(j,1),jmax)
        do k=1,kmax
         tem_t(i,j,k,:)=temobc_t(i1,j1,k,ki1)
         sal_t(i,j,k,:)=salobc_t(i1,j1,k,ki1)
        enddo
       enddo
       enddo

       if(     relaxtype_ts==2     &
          .and.relax_ts>0.    ) then !22222> !02-04-16
        if(.not.allocated(temlwf_t))allocate(temlwf_t(imax,jmax,kmax))
        if(.not.allocated(sallwf_t))allocate(sallwf_t(imax,jmax,kmax))
        temlwf_t=0.                          !02-04-16
        sallwf_t=0.
!       do k=1,kmax ; do j=1,jmax ; do i=1,imax
!        temlwf_t(i,j,k)=temobc_t(i,j,k,ki1)
!        sallwf_t(i,j,k)=salobc_t(i,j,k,ki1)
!       enddo       ; enddo       ; enddo
       endif                         !22222>

       lb4=lbound(vel_u) ; ub4=ubound(vel_u) !17-12-15
       do j=lb4(2),ub4(2)
       do i=lb4(1),ub4(1)
        i1=min(max(i,1),imax+1)
        j1=min(max(j,1),jmax)
        do k=1,kmax
         vel_u(i,j,k,:)=velobc_u(i1,j1,k,ki1)
        enddo
       enddo
       enddo

       lb4=lbound(vel_v) ; ub4=ubound(vel_v) !17-12-15
       do j=lb4(2),ub4(2)
       do i=lb4(1),ub4(1)
        i1=min(max(i,1),imax)
        j1=min(max(j,1),jmax+1)
        do k=1,kmax
         vel_v(i,j,k,:)=velobc_v(i1,j1,k,ki1)
        enddo
       enddo
       enddo

!#endif
! Echanges mpi u1 u2 v1 v2 sur vel_u(:,:,:,:) et vel_v(:,:,:,:):
      call obc_int_mpi4d !02-07-14

      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
        veldydz_u(i,j,k,1)=vel_u(i,j,k,1)                   &
                           *dy_u(i,j)                       & !30-09-09
                           *dz_u(i,j,k,1)
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
        veldxdz_v(i,j,k,1)=vel_v(i,j,k,1)                   &
                           *dx_v(i,j)                       & !30-09-09
                           *dz_v(i,j,k,1)
      enddo
      enddo
      enddo

!..............................................................
! Computes dimensionless velocity and advection stability index:
!     call couple_modes_currentnumber !03-02-16

! Le tableau thz_c(2) contient temporairement tem_c pour le schema
! d'advection. L'ordre du schema impose de connaitre des conditions
! au limites eloignées:
!     call obc_scal(1)! C.L. sur tem_t en i=-1,0,imax+1,imax+2 !01-06-09
      call obc_scal(obctype_ts) !30-07-14
      call obc_scal_mpi_now     !12-05-14
      do k=1,kmax
      do j=-1,jmax+2
      do i=-1,imax+2
        tem_t(i,j,k,0)=tem_t(i,j,k,1)                                  !14-10-09
        sal_t(i,j,k,0)=sal_t(i,j,k,1)
      enddo
      enddo
      enddo

!*****************************************************************
! Initialisation du transport et de la vitesse moyenne:
!*****************************************************************
      do k1=0,2                                                     !11-03-10
      do j=1,jmax+1
      do i=1,imax+1
          velbar_u(i,j,k1)=velbarobc_u(i,j,ki1)*mask_u(i,j,kmax+1)       !26/08/02
          velbar_v(i,j,k1)=velbarobc_v(i,j,ki1)*mask_v(i,j,kmax+1)       !26/08/02
      enddo
      enddo
      enddo
!#ifdef parallele
!      ub3=ubound(velbar_u) ; lb3=lbound(velbar_u)
!      call echange('x ',velbar_u,lb3,ub3) ! 2 pour echange 3eme arg = 2 ! C.L. i=imax+1 i=1 j=jmax j=1
!      call echange('u1',velbar_u,lb3,ub3) ! Pour schema adv/dif ordre 4
!      call echange('u2',velbar_u,lb3,ub3) ! Pour schema adv/dif ordre 4

!      ub3=ubound(velbar_v) ; lb3=lbound(velbar_v)
!      call echange('y ',velbar_v,lb3,ub3) ! 2 pour echange 3eme arg = 2 !C.L. i=imax+1 i=1 j=jmax j=1
!      call echange('v1',velbar_v,lb3,ub3) ! Pour schema adv/dif ordre 4
!      call echange('v2',velbar_v,lb3,ub3) ! Pour schema adv/dif ordre 4
!#endif
! Echanges mpi u1 u2 v1 v2 sur velbar_u(:,:,:) et velbar_v(:,:,:):
      call obc_ext_mpi3d !02-07-14

      do j=1,jmax ; do i=1,imax+1
             velavr_u(i,j,1)=velbar_u(i,j,1)
            fluxbar_u(i,j,:)=velbar_u(i,j,1)*hz_u(i,j,1)*dy_u(i,j) !11-11-14
       fluxbar_sumt_u(i,j,:)=0.
      enddo ; enddo

      do j=1,jmax+1 ; do i=1,imax
             velavr_v(i,j,1)=velbar_v(i,j,1)
            fluxbar_v(i,j,:)=velbar_v(i,j,1)*hz_v(i,j,1)*dx_v(i,j) !11-11-14
       fluxbar_sumt_v(i,j,:)=0.
      enddo ; enddo



! attention que les fichiers de grande echelle peuvent contenir
! des valeurs aberantes dans le masque. Attention donc aux points
! sources. Donc apres lecture des fichiers on reinitialise T & S
! aux sources pour eviter les catastrophes:
      call set_rivers(3)                                               !30/10/01

      do k=0,1
      do j=1,jmax+1
      do i=1,imax+1

!      fluxbarsave_u(i,j,0,k)=fluxbar_u(i,j,0)
!      fluxbarsave_u(i,j,1,k)=fluxbar_u(i,j,1)
!      fluxbarsave_v(i,j,0,k)=fluxbar_v(i,j,0)
!      fluxbarsave_v(i,j,1,k)=fluxbar_v(i,j,1)

!      sshsave_w(i,j,0,k)=ssh_ext_w(i,j,0)
!      sshsave_w(i,j,1,k)=ssh_ext_w(i,j,1)

      enddo
      enddo
      enddo

! Initialisation des tableaux de conditions aux limites:               !15/09/05
      do i=1,imax
       sshrefobc_i(i,1)=ssh_int_w(i,1   ,2)
       sshrefobc_i(i,2)=ssh_int_w(i,jmax,2)
       vbrrefobc_i(i,1)=   velbar_v(i,2   ,2)
       vbrrefobc_i(i,2)=   velbar_v(i,jmax,2)
      enddo
      do j=1,jmax
       sshrefobc_j(j,1)=ssh_int_w(1   ,j,2)
       sshrefobc_j(j,2)=ssh_int_w(imax,j,2)
       vbrrefobc_j(j,1)=   velbar_u(2   ,j,2)
       vbrrefobc_j(j,2)=   velbar_u(imax,j,2)
      enddo

      end subroutine initial_with_obc
