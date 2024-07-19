      subroutine couple_modes
!______________________________________________________________________
! SYMPHONIE ocean model
! release 338 - last update: 19-03-22
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='couple_modes'
       subroutinedescription= &
       'Replace the depth-averaged part of the internal-mode'      &
       //'-computed velocity with the depth-averaged current computed'&
       //' by the external mode equations'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
!    _________                    .__                  .__                     !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         22/11/02: annulation du schema forward
!         26/12/02: amenagements pour remettre en service le schema forward
!         10/02/03: les vitesses doivent être multipliées par leur mask
!                   pour les besoins de la grille hybride
!         10/07/03: debug du point precedent: multiplier par mask annule
!                   la vitesse au point source du fleuve, donc on change la
!                   règle habituelle de l'imbrication de boucle 3D
!                   (do k, do j, do i) et à la place on fait
!                   un calcul s'appuyant sur KLOW plutot que sur mask
!         14/06/04: pour garantir rigoureusement la conservation de la masse
!                   à la premiere iteration qui est un cas particulier quelques
!                   hypotheses ont été ajoutées sur ZTA_INT et VBAR_SUMT
!         19/09/05: suite du point precedent. Pour garantir rigoureusement la
!                   conservation de la chaleur et de la salinite, l'ajustement
!                   initial du point precedent est etendue à T et à S
!         22/07/06: ajout du masque, mais dans le doute...
!         07/08/06: Modes externes parallelles
!         11/08/06: Appel à vm_to_vmhz supprimé
!         20/04/07: Passage à coordonnée curviligne
!         07/06/07: vmean(0) represente maintenant le courant moyen de l'iteration
!                   precedente
!         24/02/08: tableau xy_x(i,j,2) remplacé par hzdy_x(i,j,1)...
! 2009.3  30-09-09  utilisation des nouveaux facteurs d'echelle verticale
!         11-10-09  suppression thz et shz
!         20-10-09  dans le cas de la nouvelle grille hybride où la premiere cellule
!                   de courant sous le fond, n'est pas forcement nulle, le calcul des
!                   moyennes impose de diviser par zerolevel (l'epaisseur où le
!                   courant est non nul, et non pas l'epaisseur totale hz.
! 2010.7  16-02-10  Le flux forward est maintenant stocké dans veldxdz(2)
!         26-02-10  ajout velstokes
!         28-02-10  ajout velstokes suite
! 2010.8  03-05-10  Avec nouveau schema forward plus besoin de veldxdz(2)
! 2010.14 12-11-10  modif sur boucle pour conservation mpi
! 2010.20 17-04-11  Calculs sur la base d'un temps en secondes
!         19-04-11  - Pas de temps modifiable
!                   - Suppression lignes mode parallele "double"
! S26.1   27-01-13  Advection qdm ordre 4 entraine des echanges supplementaires
!         14-05-13  ajout routine couple_modes_instant
!         18-01-14  ajout routine couple_modes_3v
!         02-07-14  call obc_ext_mpi_xy(0)
!         21-08-14  ajout subroutine couple_modes_3dback2d
!                         subroutine couple_modes_currentnumber
!                         subroutine couple_modes_bounds
!         06-09-14  ajout d'un diag sur nombre de courant
!         29-10-14  reshape loop
!         31-10-14  index de stabilite pondere par wetwask
!         03-11-14  attention a multiplication par mask!
!         17-11-14  alternative (non activee actuellement) au calcul du
!                   coef de stabilite
!         28-11-14  suppression ssh_int_u ssh_int_v
!         13-02-15  stability_index est calcule seulement si flag_dt_adjust=1
!         15-07-15  couple_modes_constant_corr  et couple_modes_zprofil_corr  
!                   permettent de choisir la forme du profil de courant qui 
!                   corrige la moyenne verticale du courant 3d
!         05-08-15  methode vst continue vs discontinue
!         16-10-15  correction d'un bug introduit au point precedent 
!         28-11-15  ajout flag_1dv
!         21-12-15  boucles fortran
!         28-01-16  Amenagements pour nouveau time stepping hybride LF-FB
!         02-02-16  le calcul du nombre de courant doit etre remis a jour pour
!                   l'advection de T,S si shema LFFB
!         14-02-16  Amenagement modele 1DV
!         28-03-16  Leap-Frog correction de vel(t) avec moyenne temporelle du 
!                   mode externe - implication mpi
!         07-08-16  appel A couple_modes_currentnumber est commentE car transfere
!                   dans le schema d'advection
!         04-09-16  time-splitting LF corriger veldydz seulement pour eviter singularite du wetdry
!         29-09-16  correction du courant moyen instantannE de vel(2) en
!                   prenant compte d'un profil dans la couche de fond
!         02-10-16  call wetdry_mask_airseafluxes 
!         03-10-16  modif par rapport aux flux entrant des fleuves
!         15-12-16  dans couple_modes_scalars calculer w en passant par le hub omega.F90
!         05-02-17  ligne inutile commentee
!         21-04-17  s-z coordinate
! v247    25-02-19  wetdry schema 6
! v285    05-06-20  suppression zeroleveldepth
! v310    02-11-21  if(flag_merged_levels==1) then !- REDUCTEUR DE FLUX - COUCHE KMERGED-1 ---> !02-11-21
!         04-11-21  if(flag_merged_levels==1) then !- REDUCTEUR DE FLUX - COUCHE KMERGED-1 ---> !04-11-21
! v311    14-11-21  min(hz_u(i,j,1)*dsigmerged_u(i,j),dz_u(i,j,k,1)) et idem pour _v
! v312    17-11-21  nouvelle formule pour dsigmerged_u
! v326    12-12-21  pourquoi ne pas filtrer le mode numerique dans couple_instant?
!                   On le fait donc desormais (sans plus d'argument que l'intuition)
! v338    19-03-22  ameliorer le test pour garantir que les flux de volume sont bien nuls
!                   en dessous de kmin
!...............................................................................

! Corrige la valeur du courant moyen (sur z) issue des calculs
! du mode interne par celle calculée par le mode externe.
! Ref: Blumberg et Mellor 1987.

!......................................
! Particular case of the initial state:
      if(iteration3d==kount0) then !00000000000000>

!      do j=0,jmax+1                                                    !12-11-10
!      do i=0,imax+1
             ssh_int_w(:,:,0)=  2.*ssh_int_w(:,:,1)-ssh_int_w(:,:,2) !21-12-15
        fluxbar_sumt_u(:,:,0)=fluxbar_sumt_u(:,:,1)
        fluxbar_sumt_v(:,:,0)=fluxbar_sumt_v(:,:,1)
!      enddo
!      enddo

       call z_thickness(0,0)                                            !19/09/05
       call cellbox_thickness(0,0,0)                                    !19/09/05
       call wetdry_mask_airseafluxes !02-10-16

      endif                        !00000000000000>


!......................................................
! Save the "old" z-averaged current in temporary array:
      do j=1,jmax+1
      do i=1,imax+1
! boucle stricte: do i=1,imax+1 | j=1,jmax
       xy_u(i,j,1)=velavr_u(i,j,1)
! boucle stricte: do i=1,imax | j=1,jmax+1
       xy_v(i,j,1)=velavr_v(i,j,1)
      enddo
      enddo


!..........................................................
! Save the "new" z-averaged current in velavr_u & velavr_v:
      const2=1.0/real(iteration2d_max_now+iteration2d_max_bef)    !19-04-11

      do j=1,jmax                                                  !29-10-14
      do i=1,imax+1

! vitesse moyenne "leapfrog ou centrée":
      velavr_u(i,j,1)=                                  &
       const2*(fluxbar_sumt_u(i,j,0)                    &
              +fluxbar_sumt_u(i,j,1))                   &
                        /dy_u(i,j)                      &
                        /hz_u(i,j,1)                    & !20-10-09 !05-06-20
              -velbarstokes_u(i,j,1)                      !28-02-10

      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax                                               !29-10-14

      velavr_v(i,j,1)=                                  &
       const2*(fluxbar_sumt_v(i,j,0)                    &
              +fluxbar_sumt_v(i,j,1))                   &
                        /dx_v(i,j)                      &
                        /hz_v(i,j,1)                    &               !20-10-09
              -velbarstokes_v(i,j,1)                                    !28-02-10

      enddo
      enddo

!..........................................................................
! Replace the "old" z-averaged component of the 3d current by the "new" one
! in vel_u and vel_v arrays:
!     if(flag_1dv==0) then ! 3D case > !28-11-15! A commenter si modele 1dv pilotE par velbar
       if(zprofile2d3d==0) then !>>>
        call couple_modes_constant_corr ! with a constant (z independent) value
       else                     !>>>
        call couple_modes_zprofil_corr  ! with a z-dependent Ekman bottom profile
       endif                    !>>>
!     endif                ! 3D case >

      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1

         veldydz_u(i,j,k,1)=(              &
            anyv3d(i,j,k,id_u_now)         &
      +velstokes_u(i,j,k,1) )              & !26-02-10
             *dy_u(i,j)                    &
             *dz_u(i,j,k,1)                & !30-09-09
           *mask_u(i,j,k ) 

       enddo ; enddo ; enddo

       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax

         veldxdz_v(i,j,k,1)=(              &
            anyv3d(i,j,k,id_v_now)         &
      +velstokes_v(i,j,k,1) )              & !26-02-10
             *dx_v(i,j)                    & !16-10-15
             *dz_v(i,j,k,1)                & !30-09-09
           *mask_v(i,j,k )

       enddo ; enddo ; enddo 

      if(flag_merged_levels==1) then !- REDUCTEUR DE FLUX - COUCHE KMERGED-1 ---> !04-11-21

      do j=1,jmax ; do i=1,imax+1
!     if(kmin_u(i,j)/=kmerged_u(i,j)) then !m°v°m>

       sum1=0. ; sum2=0.
       do k=1,kmerged_u(i,j)-2 
        sum1=sum1+veldydz_u(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite
        veldydz_u(i,j,k,1)=0.        ! nouveau flux
       enddo

       k=kmerged_u(i,j)-1
       sum1=sum1+veldydz_u(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite

      if(kmin_u(i,j)/=kmerged_u(i,j)) then !m°v°m> deplacE le !19-03-22
! Formule standard mais dz est remplacE par h*dsigmerged
         veldydz_u(i,j,k,1)=(               &
            anyv3d(i,j,k,id_u_now)          &
      +velstokes_u(i,j,k,1) )               & !26-02-10
             *dy_u(i,j)                     &
             *dz_u(i,j,k,1)                 & !30-09-09
     *dsigmerged_u(i,j)                     & !17-11-21
!        *min(hz_u(i,j,1)*dsigmerged_u(i,j),dz_u(i,j,k,1)) & !04-11-21!14-11-21
           *mask_u(i,j,k ) 
      else                                 !m°v°m>
         veldydz_u(i,j,k,1)=0.  !19-03-22           
      endif                                !m°v°m>
      sum2=sum2+veldydz_u(i,j,k,1) ! Erreur de bilan de tranport du nouveau flux A corriger ensuite

! La couche kmerged_u recupere le defaut de bilan de transport sum2-sum1
       k=kmerged_u(i,j)
         veldydz_u(i,j,k,1)= &
         veldydz_u(i,j,k,1)+sum1-sum2

!     endif                                !m°v°m>
      enddo ; enddo

      do j=1,jmax+1 ; do i=1,imax
!     if(kmin_v(i,j)/=kmerged_v(i,j)) then !m°v°m>

       sum1=0. ; sum2=0.
       do k=1,kmerged_v(i,j)-2 
        sum1=sum1+veldxdz_v(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite
        veldxdz_v(i,j,k,1)=0.        ! nouveau flux
       enddo

       k=kmerged_v(i,j)-1
       sum1=sum1+veldxdz_v(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite
      if(kmin_v(i,j)/=kmerged_v(i,j)) then !m°v°m> !19-03-22
! Formule standard mais dz est remplacE par h*dsigmerged
         veldxdz_v(i,j,k,1)=(               &
            anyv3d(i,j,k,id_v_now)          &
      +velstokes_v(i,j,k,1) )               & !26-02-10
             *dx_v(i,j)                     & !16-10-15
             *dz_v(i,j,k,1)                 & !30-09-09
     *dsigmerged_v(i,j)                     & !17-11-21
!        *min(hz_v(i,j,1)*dsigmerged_v(i,j),dz_v(i,j,k,1)) & !14-11-21
           *mask_v(i,j,k )
      else                                 !m°v°m>
         veldxdz_v(i,j,k,1)=0. !19-03-22
      endif                                !m°v°m>
      sum2=sum2+veldxdz_v(i,j,k,1) ! Erreur de bilan de tranport du nouveau flux A corriger ensuite

! La couche kmerged_v recupere le defaut de bilan de transport sum2-sum1
       k=kmerged_v(i,j)
         veldxdz_v(i,j,k,1)= &
         veldxdz_v(i,j,k,1)+sum1-sum2

!     endif                                !m°v°m>
      enddo ; enddo

      endif                          !- REDUCTEUR DE FLUX - COUCHE KMERGED-1 --->

! L'echange est inutile si on ne corrige pas vel(1)
! Echanges au temps t_=1 de type=2  (u2 et v2)
      if(timestep_type==timestep_leapfrog)call obc_int_mpi(1,2)  ! (t_,obctype_) avec: !28-03-16
                                                                 ! t_=1 ou 2
                                                                 ! obctype=1  pour (u1,v1) 
                                                                 ! obctype=2  pour (u2,v2) 
                                                                 ! obctype=12 pour (u1,v1) & (u2,v2)

!..........
! Limiteur:
!     call couple_modes_bounds

!..............................................................
! Computes dimensionless velocity and advection stability index:
!     call couple_modes_currentnumber ! commentee le 07-08-16 

      end subroutine couple_modes

!......................................................................
#ifdef bidon
      subroutine couple_modes_currentnumber !21-08-14
      use module_principal
      use module_parallele
      use module_s
      implicit none
      integer t1_,t2_
#ifdef synopsis
       subroutinetitle='couple_modes_currentnumber'
       subroutinedescription=                                         &
          'Computes the dimensionless velocity field u*dt/dx v*dt/dy' &
       //' derived from the "advecting" velocity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      stop 'je ne dois plus passer par couple_modes_currentnumber'
#ifdef bidon
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

! u*dt/dx:
!     do j=1,jmax
!     do i=1,imax+1
!      x1=dti_lp/dx_u(i,j)
!      do k=1,kmax
!         veldtodx_u(i,j,k,1)=               &
!             (vel_u(i,j,k,1)                &
!       +velstokes_u(i,j,k,1) )*x1
!      enddo
!     enddo
!     enddo

! v*dt/dy
!     do j=1,jmax+1
!     do i=1,imax
!      x1=dti_lp/dy_v(i,j)
!      do k=1,kmax
!         veldtody_v(i,j,k,1)=               &
!             (vel_v(i,j,k,1)                &
!       +velstokes_v(i,j,k,1) )*x1
!      enddo
!     enddo
!     enddo

      if(flag_dt_adjust==0)return !13-02-15
      stop 'OPTION REMPLACEE PAR ADVECTION SEPAREE'

#ifdef bidon
! Calcul alternatif au coef de stabilite !17-11-14
! Stability index:
      x1=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       x1=max(x1,                                                      &!10-11-14

                 wetmask_t(i,j)*dti_lp*(                               &

        ((veldydz_u(i+1,j,k,1)+abs(veldydz_u(i+1,j,k,1)))              &
        -(veldydz_u(i  ,j,k,1)-abs(veldydz_u(i  ,j,k,1)))              &
        +(veldxdz_v(i,j+1,k,1)+abs(veldxdz_v(i,j+1,k,1)))              &
        -(veldxdz_v(i,j  ,k,1)-abs(veldxdz_v(i,j  ,k,1))))/dxdy_t(i,j) &
!         +(omega_w(i,j,k+1,1)+abs(  omega_w(i,j,k+1,1)))              &
!         -(omega_w(i,j,k  ,1)-abs(  omega_w(i,j,k  ,1)))              &

                                       )/(2.*dz_t(i,j,k,2))            &
              )

      enddo ; enddo ; enddo
      call mpi_allreduce(x1,stability_index,1,mpi_double_precision,mpi_max,par%comm2d ,ierr)
#endif

! Ecrire le nombre de courant max et moyenne dans un fichier !06-09-14
      sum1=small1 ; sum2=0. ; x1=0.
      do k=1,kmax ; do j=1,jmax   ; do i=1,imax+1
       sum1=sum1+mask_u(i,j,k)
       sum2=sum2+mask_u(i,j,k)*abs(veldtodx_u(i,j,k,1))
                     x1=max(x1,abs(veldtodx_u(i,j,k,1)*wetmask_u(i,j))) !31-10-14
      enddo       ; enddo         ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
       sum1=sum1+mask_v(i,j,k)
       sum2=sum2+mask_v(i,j,k)*abs(veldtody_v(i,j,k,1))
                     x1=max(x1,abs(veldtody_v(i,j,k,1)*wetmask_v(i,j))) !31-10-14
      enddo       ; enddo         ; enddo
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(x1,stability_index,1,mpi_double_precision,mpi_max,par%comm2d ,ierr)
      if(par%rank==0) then !-------->
       k0=s_unit(7) ! s_unit renvoie un file unit libre, si possible=argument (module_s)
       open(unit=k0,file=trim(tmpdirname)//'current_number_mean_max',position='append')
       write(k0,'(5(e14.7,1x),i3)')elapsedtime_now/86400.   &
                                   ,sum2glb/sum1glb          &
                                   ,stability_index          &
                                   ,dti_fw                   &
                                   ,real(iteration2d_max_r8) &
                                   ,iteration2d_max_now
       close(k0)
      endif                !------->

#endif
      end subroutine couple_modes_currentnumber
#endif
!......................................................................
!#ifdef bidon
      subroutine couple_modes_instant(t1_,t2_)
      use module_principal
      use module_parallele
      implicit none
      integer t1_,t2_
#ifdef synopsis
       subroutinetitle='couple_modes_instant'
       subroutinedescription= &
       'Replace the depth-averaged part of the 3D velocity at'       &
       //' time=after with the latest (i.e. not averaged over time)' &
       //' depth-averaged current computed by the external mode'     &
       //' equations.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
  
     if(ihybsig/=0)stop 'couple_modes_instant ihybsig/=0'

! reset en K=1:
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=vel_u(i,j,1,t1_)*dz_u(i,j,1,t1_)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=vel_v(i,j,1,t1_)*dz_v(i,j,1,t1_)
      enddo
      enddo

! sommation sur K suivants:
      do k=2,kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=xy_u(i,j,1)          &
                 +vel_u(i,j,k,t1_)      &
                  *dz_u(i,j,k,t1_)
      enddo
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax+1
      do i=1,imax
         xy_v(i,j,1)=xy_v(i,j,1)          &
                   +vel_v(i,j,k,t1_)      &
                    *dz_v(i,j,k,t1_)
      enddo
      enddo
      enddo

! Delta de courant moyen:
!     sum1=0.
!     sum2=0.
      do j=1,jmax
      do i=1,imax+1
!      if(iteration3d==3.and.mask_u(i,j,kmax)==1) then
!         write(10+par%rank,*)i,j,velbar_u(i,j,t2_),xy_u(i,j,1)/hz_u(i,j,t1_)
!         sum1=sum1+abs( velbar_u(i,j,t2_)-xy_u(i,j,1)/hz_u(i,j,t1_) )
!      endif

!      xy_u(i,j,2)=velbar_u(i,j,t2_)-xy_u(i,j,1)/hz_u(i,j,t1_)
       xy_u(i,j,2)=0.5*(velbar_u(i,j,2)+velbar_u(i,j,1))-xy_u(i,j,1)/hz_u(i,j,t1_) !12-12-21

!      sum1=sum1+1.
!      sum2=sum2+xy_u(i,j,2)**2
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
!      if(iteration3d==3.and.mask_v(i,j,kmax)==1) then
!          write(20+par%rank,*)i,j,velbar_v(i,j,t2_),xy_v(i,j,1)/hz_v(i,j,t1_)
!         sum2=sum2+abs( velbar_v(i,j,t2_)-xy_v(i,j,1)/hz_v(i,j,t1_) )
!      endif

!      xy_v(i,j,2)=velbar_v(i,j,t2_)-xy_v(i,j,1)/hz_v(i,j,t1_)
       xy_v(i,j,2)=0.5*(velbar_v(i,j,2)+velbar_v(i,j,1))-xy_v(i,j,1)/hz_v(i,j,t1_) !12-12-21

!      sum1=sum1+1.
!      sum2=sum2+xy_v(i,j,2)**2
      enddo
      enddo
!     write(40+par%rank,*)elapsedtime_now/86400.,sqrt(sum2/sum1)

!     if(iteration3d==3) then
!      write(6,*)'sum1 sum2 ',par%rank,sum1,sum2
!     endif
!     if(iteration3d==4)stop 'bac'

! Correction du courant 3D:
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
       vel_u(i,j,k,t1_)=vel_u(i,j,k,t1_)+xy_u(i,j,2)
      enddo
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
       vel_v(i,j,k,t1_)=vel_v(i,j,k,t1_)+xy_v(i,j,2)
      enddo
      enddo
      enddo

! Verification
!     i=imax/2 ; j=jmax/2
!     sum1=0. ; sum2=0.
!     do k=1,kmax
!      sum1=sum1+dz_u(i,j,k,t1_)*vel_u(i,j,k,t1_)
!      sum2=sum2+dz_v(i,j,k,t1_)*vel_v(i,j,k,t1_)
!     enddo
!     write(6,*)'verif u',par%rank,sum1/hz_u(i,j,t1_),velbar_u(i,j,t2_)
!     write(6,*)'verif v',par%rank,sum2/hz_v(i,j,t1_),velbar_v(i,j,t2_)

      end subroutine couple_modes_instant
!#endif
!............................................................................
#ifdef bidon
      subroutine couple_modes_bounds !21-08-14
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='couple_modes_bound'
       subroutinedescription= &
          'Applies ranging bounds on the "advecting" velocity field' &
       //' veldydz_u and veldxdz_v'

       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Limites:
      x1=1. ! upper/lower bound for the dimensionless velocity (u*dt/dx,v*dt/dy)
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
          veldydz_u(i,j,k,1)=max(min(veldydz_u(i,j,k,1)          &
                         ,x1*dxdy_u(i,j)*dz_u(i,j,k,2)/dti_lp)  &
                        ,-x1*dxdy_u(i,j)*dz_u(i,j,k,2)/dti_lp)
!      if(veldydz_u(i,j,k,1)>+x1*dxdy_u(i,j)*dz_u(i,j,k,2)/dti_lp) then
!         veldydz_u(i,j,k,1)=+x1*dxdy_u(i,j)*dz_u(i,j,k,2)/dti_lp
!         write(6,*)'A ',par%rank,i,j,k
!         stop
!      endif
!      if(veldydz_u(i,j,k,1)<-x1*dxdy_u(i,j)*dz_u(i,j,k,2)/dti_lp) then
!         veldydz_u(i,j,k,1)=-x1*dxdy_u(i,j)*dz_u(i,j,k,2)/dti_lp
!         write(6,*)'B ',par%rank,i,j,k
!         stop
!      endif

      enddo
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
          veldxdz_v(i,j,k,1)=max(min(veldxdz_v(i,j,k,1)          &
                         ,x1*dxdy_v(i,j)*dz_v(i,j,k,2)/dti_lp)  &
                        ,-x1*dxdy_v(i,j)*dz_v(i,j,k,2)/dti_lp)
!      if(veldxdz_v(i,j,k,1)>+x1*dxdy_v(i,j)*dz_v(i,j,k,2)/dti_lp) then
!         veldxdz_v(i,j,k,1)=+x1*dxdy_v(i,j)*dz_v(i,j,k,2)/dti_lp
!         write(6,*)'C ',par%rank,i,j,k
!         stop
!      endif
!      if(veldxdz_v(i,j,k,1)<-x1*dxdy_v(i,j)*dz_v(i,j,k,2)/dti_lp) then
!         veldxdz_v(i,j,k,1)=-x1*dxdy_v(i,j)*dz_v(i,j,k,2)/dti_lp
!         write(6,*)'D ',par%rank,i,j,k
!         stop
!      endif

      enddo
      enddo
      enddo

      end subroutine couple_modes_bounds
#endif
!...................................................................................
#ifdef bidon
      subroutine couple_modes_3dback2d !21-08-14
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='couple_modes_3dback2d'
       subroutinedescription=                                        &
          'Adjusts the velocity field computed by the external mode' &
       //' with the z-averaged internal mode velocity field at'      &
       //' time=after'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! reset en K=1:
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=vel_u(i,j,1,2)*dz_u(i,j,1,2)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=vel_v(i,j,1,2)*dz_v(i,j,1,2)
      enddo
      enddo

! sommation sur K suivants:
      do k=2,kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=xy_u(i,j,1)        &
                 +vel_u(i,j,k,2)      &
                  *dz_u(i,j,k,2)
      enddo
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax+1
      do i=1,imax
         xy_v(i,j,1)=xy_v(i,j,1)        &
                   +vel_v(i,j,k,2)      &
                    *dz_v(i,j,k,2)
      enddo
      enddo
      enddo

! Delta de courant moyen:
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,2)=-velbar_u(i,j,2)+xy_u(i,j,1)/hz_u(i,j,2)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,2)=-velbar_v(i,j,2)+xy_v(i,j,1)/hz_v(i,j,2)
      enddo
      enddo

! Correction du courant 2D:
      do loop1=-1,2
      do j=1,jmax
      do i=1,imax+1
       velbar_u(i,j,loop1)=velbar_u(i,j,loop1)+xy_u(i,j,2)
      enddo
      enddo
      enddo
      do loop1=-1,2
      do j=1,jmax+1
      do i=1,imax
       velbar_v(i,j,loop1)=velbar_v(i,j,loop1)+xy_v(i,j,2)
      enddo
      enddo
      enddo

      end subroutine couple_modes_3dback2d
#endif
!............................................................................
#ifdef bidon
      subroutine couple_modes_3v !18-01-14
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='couple_modes_3v'
       subroutinedescription= &
       'Builds the tracer advection velocity with a 3-times (before-' &
       //'now-after) average of the current field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Advection de T, S, Tke... avec la moyenne des trois vitesses vel(0) vel(1) vel(2)
! On ajuste la moyenne verticale pour que le courant barotrope equilibre la variation
! de ssh entre t=0 et t=2 en se servant de velavr
! Le but est de produire un nouveau veldxdz_v & veldydz_u
! On calcule ensuite omega (on peut verifier que omega(kmax+1) est bien nul

! reset en K=1:
      k=1
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=(0.25*( vel_u(i,j,k,0)     &
                          +vel_u(i,j,k,2))    &
                    +0.5*  vel_u(i,j,k,1)     &
                          )*dz_u(i,j,k,1)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=(0.25*( vel_v(i,j,k,0)     &
                          +vel_v(i,j,k,2))    &
                    +0.5*  vel_v(i,j,k,1)     &
                          )*dz_v(i,j,k,1)
      enddo
      enddo

! sommation sur K suivants:
      do k=2,kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=                           &
       xy_u(i,j,1)+(0.25*( vel_u(i,j,k,0)     &
                          +vel_u(i,j,k,2))    &
                    +0.5*  vel_u(i,j,k,1)     &
                          )*dz_u(i,j,k,1)
      enddo
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=                           &
       xy_v(i,j,1)+(0.25*( vel_v(i,j,k,0)     &
                          +vel_v(i,j,k,2))    &
                    +0.5*  vel_v(i,j,k,1)     &
                          )*dz_v(i,j,k,1)
      enddo
      enddo
      enddo

! Delta de courant moyen: (voir calcul de velavr_ dans couple_modes)
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,0)=velavr_u(i,j,1)-xy_u(i,j,1)/hz_u(i,j,1)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,0)=velavr_v(i,j,1)-xy_v(i,j,1)/hz_v(i,j,1)
      enddo
      enddo
! Echange mpi u1 u2 v1 v2 sur xy_u(:,:,0) et xy_v(:,:,0)
      call obc_ext_mpi_xy(0) !02-07-14
!#ifdef parallele
!      ub3=ubound(xy_u) ; lb3=lbound(xy_u)
!      call echange('u1',xy_u,lb3,ub3,0) ! Pour schema adv/dif ordre 4 !27-01-13
!      call echange('u2',xy_u,lb3,ub3,0) ! Pour schema adv/dif ordre 4

!      ub3=ubound(xy_v) ; lb3=lbound(xy_v)
!      call echange('v1',xy_v,lb3,ub3,0) ! Pour schema adv/dif ordre 4
!      call echange('v2',xy_v,lb3,ub3,0) ! Pour schema adv/dif ordre 4
!#endif

      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
        veldydz_u(i,j,k,1)=(                         &

                   (0.25*( vel_u(i,j,k,0)            &
                          +vel_u(i,j,k,2))           &
                    +0.5*  vel_u(i,j,k,1)            &
                           +xy_u(i,j,0))             &
                         *mask_u(i,j,k)              &

                     +velstokes_u(i,j,k,1)           & !26-02-10
                           )*dy_u(i,j)*dz_u(i,j,k,1)   !30-09-09
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
        veldxdz_v(i,j,k,1)=(                         &

                   (0.25*( vel_v(i,j,k,0)            &
                          +vel_v(i,j,k,2))           &
                    +0.5*  vel_v(i,j,k,1)            &
                           +xy_v(i,j,0))             &
                         *mask_v(i,j,k)              &

                     +velstokes_v(i,j,k,1)           & !26-02-10
                           )*dx_v(i,j)*dz_v(i,j,k,1)   !30-09-09
      enddo
      enddo
      enddo

      call omega
      end subroutine couple_modes_3v
#endif
!.......................................................................
      subroutine couple_modes_constant_corr
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='couple_modes_constant_corr'
       subroutinedescription=' '
#endif

! Ici j'ai imaginE corriger la valeur instantanne barotrope du courant 3D au temps 1 (pour economiser de la cpu)
! mais je crois que corriger vel(1) et ne pas corriger vel(2) risque d'aboutir A une moyenne bancale vel(1)+vel(2)
! qui sera utilisEe pour l'advection FB de T et S. Je crois donc plus propre de corriger vel(2) juste apres le
! calcul des moments afin que vel(1) et vel(2) soit plus homogene pour l'advection de T,S

      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,0)=-xy_u(i,j,1)+velavr_u(i,j,1) ! pour corriger veldydz avec moyenne temporelle de velbar
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,0)=-xy_v(i,j,1)+velavr_v(i,j,1)
      enddo
      enddo

        sum0=0.
        sum2=0.

      id_u_now=6
      do k=1,kmax
      do j=1,jmax  
      do i=1,imax+1
! attention a ne pas tout multiplier par mask_t sinon pas de river flux !03-11-14
!       vel_u(i,j,k,1)=        vel_u(i,j,k,1)+xy_u(i,j,0) ! correction moyennee sur cycle barotrope
        anyv3d(i,j,k,id_u_now)=(vel_u(i,j,k,1)+xy_u(i,j,0))*mask_u(i,j,k)

!       if(i==25.and.j==jmax/2) then
!           sum0=sum0+velavr_u(i,j,1)*dz_u(i,j,k,1)*mask_u(i,j,k)
!           sum2=sum2+                dz_u(i,j,k,1)*mask_u(i,j,k)
!           write(6,*)k,mask_u(i,j,k),dz_u(i,j,k,1),sum2
!       endif

      enddo
      enddo
      enddo
      checkanyv3d(id_u_now)=anyv3d(imax/2,jmax/2,kmax,id_u_now)
! checkanyv3d sert A verifier plus tard que anyv3d(:,:,:,,id_u_now) n'a pas EtE corrompu
! par des archivages intermediairess

      id_v_now=7
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax  
! attention a ne pas tout multiplier par mask_t sinon pas de river flux !03-11-14
!       vel_v(i,j,k,1)=        vel_v(i,j,k,1)+xy_v(i,j,0) ! correction moyennee sur cycle barotrope
        anyv3d(i,j,k,id_v_now)=(vel_v(i,j,k,1)+xy_v(i,j,0))*mask_v(i,j,k)
      enddo
      enddo
      enddo
      checkanyv3d(id_v_now)=anyv3d(imax/2,jmax/2,kmax,id_v_now)
! checkanyv3d sert A verifier plus tard que anyv3d(:,:,:,,id_u_now) n'a pas EtE corrompu
! par des archivages intermediairess

! Le Leap-Frog "POM like" suppose de corriger egalement vel(t) option
! eventuellement discutable. Si on ne corrige pas vel(t) alors penser
! A corriger vel(t+1) avec la valeur instantannEe de velbar. Operation
! A effectuer depuis la fin de internal_mode.
! Ces lignes sont finalement commentEes pour eviter la singularitE de la ! !04-09-16
! zone wetdry qui, pour un flux donnE, peut amener A des courants ! infinis
! si hz_u(:,:,1)=0. Ceci implique de corriger vel(t+1) avec la valeur
! instantannEe de velbar
!     if(timestep_type==timestep_leapfrog) then !>>> !28-03-16
!      id_u_now=6
!      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
!       vel_u(i,j,k,1)=anyv3d(i,j,k,id_u_now)
!      enddo ; enddo ; enddo
!      id_v_now=7
!      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax  
!       vel_v(i,j,k,1)=anyv3d(i,j,k,id_v_now)
!      enddo ; enddo ; enddo
!     endif                                     !>>>

      end subroutine couple_modes_constant_corr

!.......................................................................

      subroutine couple_modes_zprofil_corr !14-07-15
      use module_principal
      use module_parallele
      implicit none
!     real*4 inv_ekman_depth
#ifdef synopsis
       subroutinetitle='couple_modes_zprofil_corr'
       subroutinedescription=' '
#endif
! Mes notes:
! https://docs.google.com/document/d/1kXMG-3WaF1Miif5-gj1kzN8dSpD5Mhhowg0IzoM1Afk/edit

!     inv_ekman_depth=1./10. ! d**-1=sqrt(f/(2.Kz)) ! sqrt(0.5*10e-4/5e-2)
!     inv_ekman_depth=1./10. ! 1/(Ekman depth)

      stop 'couple_modes_zprofil_corr n est pas A jour!'

!.............
! U component:
      k=kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,3)=                                             &
                      dz_u(i,j,k,1)                             &
!       *(1.-exp(-(depth_u(i,j,k)-depth_u(i,j,kmin_u(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_u(i,j,k)+h_u(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_u(i,j,k)                               &
                +(1-mask_u(i,j,k))*small1

      enddo
      enddo
      do k=kmax-1,1,-1
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,3)=                                             &
       xy_u(i,j,3)+                                             &
                      dz_u(i,j,k,1)                             &
!       *(1.-exp(-(depth_u(i,j,k)-depth_u(i,j,kmin_u(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_u(i,j,k)+h_u(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_u(i,j,k)
      enddo
      enddo
      enddo

      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,4)=    (-xy_u(i,j,1)     &
                    +velavr_u(i,j,1))    &
                        *hz_u(i,j,1)     &
                        /xy_u(i,j,3)
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1

                   vel_u(i,j,k,1)=                            &
                   vel_u(i,j,k,1)                             &
                   +xy_u(i,j,4)                               &
!     *(1.-exp(-(depth_u(i,j,k)-depth_u(i,j,kmin_u(i,j)-1))*inv_ekman_depth )) & ! methode 1
      *(1.-exp(-(depth_u(i,j,k)+h_u(i,j))*inv_ekman_depth )) & ! methode 2
                 *mask_u(i,j,k)

      enddo
      enddo
      enddo

!.............
! V component:
      k=kmax
      do j=1,jmax+1
      do i=1,imax
! Note on suppose que le courant est nul en z=-h_v ce qui
! n'est pas toujours vrai en coordonnee hybride mais bon....
       xy_v(i,j,3)=                                             &
                      dz_v(i,j,k,1)                             &
!       *(1.-exp(-(depth_v(i,j,k)-depth_v(i,j,kmin_v(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_v(i,j,k)+h_v(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_v(i,j,k)                               &
                +(1-mask_v(i,j,k))*small1
      enddo
      enddo
      do k=kmax-1,1,-1
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,3)=                                             &
       xy_v(i,j,3)+                                             &
                      dz_v(i,j,k,1)                             &
!       *(1.-exp(-(depth_v(i,j,k)-depth_v(i,j,kmin_v(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_v(i,j,k)+h_v(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_v(i,j,k)
      enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,4)=(    -xy_v(i,j,1)    &
                    +velavr_v(i,j,1) )  &
                        *hz_v(i,j,1)    &
                        /xy_v(i,j,3)
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
                   vel_v(i,j,k,1)=                            &
                   vel_v(i,j,k,1)                             &
                   +xy_v(i,j,4)                               &
!     *(1.-exp(-(depth_v(i,j,k)-depth_v(i,j,kmin_v(i,j)-1))*inv_ekman_depth )) & ! methode 1
      *(1.-exp(-(depth_v(i,j,k)+h_v(i,j))*inv_ekman_depth )) & ! methode 2
                 *mask_v(i,j,k)

      enddo
      enddo
      enddo

      end subroutine couple_modes_zprofil_corr

!.............................................................................................
#ifdef bidon
      subroutine couple_modes_3ddrives2d
      use module_principal
      use module_parallele
      implicit none

! Produit restoring3d2d_u et restoring3d2d_v les termes de rappel du mode externe vers
! le mode interne quand celui est en mode maItre.

      k=1
      do j=2,jmax-1 ; do i=2,imax
       xy_u(i,j,1)=dz_u(i,j,k,2)*vel_u(i,j,k,2)*mask_u(i,j,k)
      enddo       ; enddo
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax
       xy_u(i,j,1)=xy_u(i,j,1)+dz_u(i,j,k,2)*vel_u(i,j,k,2)*mask_u(i,j,k)
      enddo       ; enddo       ; enddo
      do j=2,jmax-1 ; do i=2,imax
       xy_u(i,j,1)=xy_u(i,j,1)/hz_u(i,j,2)
      enddo       ; enddo
! note: en toute rigueur il faudrait diviser par hz_u mais ce tableau n'est pas
! disponible au temps 2. Ce n'est A priori pas grave car il ne s'agit que de produire qu'un
! terme de rappel autrement dit une correction de toute faCon approximative

      k=1
      do j=2,jmax ; do i=2,imax-1
       xy_v(i,j,1)=dz_v(i,j,k,2)*vel_v(i,j,k,2)*mask_v(i,j,k)
      enddo       ; enddo
      do k=2,kmax ; do j=2,jmax ; do i=2,imax-1
       xy_v(i,j,1)=xy_v(i,j,1)+dz_v(i,j,k,2)*vel_v(i,j,k,2)*mask_v(i,j,k)
      enddo       ; enddo       ; enddo
      do j=2,jmax ; do i=2,imax-1
       xy_v(i,j,1)=xy_v(i,j,1)/hz_v(i,j,2)
      enddo       ; enddo

       do j=2,jmax-1 ; do i=2,imax
        restoring3d2d_u(i,j)=0.05*(xy_u(i,j,1)-0.5*(velbar_u(i,j,0)+velbar_u(i,j,1)))/dti_fw
       enddo       ; enddo

       do j=2,jmax ; do i=2,imax-1
        restoring3d2d_v(i,j)=0.05*(xy_v(i,j,1)-0.5*(velbar_v(i,j,0)+velbar_v(i,j,1)))/dti_fw
       enddo       ; enddo

       
      end subroutine couple_modes_3ddrives2d
#endif
!.............................................................................................

      subroutine couple_modes_scalars
      use module_principal
      use module_parallele
      implicit none

! composante U:

      k0=0
      do k=1,kmax
       do j=1,jmax
       do i=1,imax+1
        xy_u(i,j,1)=xy_u(i,j,1)*k0                &
       +dz_u(i,j,k,1)*(  0.5*( vel_u(i,j,k,1)     &
                              +vel_u(i,j,k,2))    &
                        +velstokes_u(i,j,k,1)  )  &
                             *mask_u(i,j,k)
       enddo
       enddo
      k0=1
      enddo

      const2=1.0/real(iteration2d_max_now+iteration2d_max_bef)    !19-04-11
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,2)=( & !gggg>
                      const2*(fluxbar_sumt_u(i,j,0)   &
                             +fluxbar_sumt_u(i,j,1))  &
                                       /dy_u(i,j)     &

                                      -xy_u(i,j,1)    &
                   ) & !gggg>
                          /hz_u(i,j,1)    
      enddo
      enddo

      do k=1,kmax
! si method vst discontinue k1=k sinon (sigma ou vst continue) k1=kmax
!     k1=vststep*k+(1-vststep)*kmax !05-08-15
       do j=1,jmax
       do i=1,imax+1

         veldydz_u(i,j,k,1)=( & !rrrr>
       0.5*( vel_u(i,j,k,1)           &
            +vel_u(i,j,k,2))          &
             +xy_u(i,j,2)             &
      +velstokes_u(i,j,k,1) ) & !rrrr>
             *dy_u(i,j)       &
             *dz_u(i,j,k,1)   &    
        *wetmask_u(i,j)       & !25-02-19 wetdry schema 6
! La ligne suivante (capilotractee) sert A ne pas perdre le flux aux points fleuves
! calculE par obc_river appelE depuis internal_mode !03-10-16
!          *mask_u(i,j,k1)+(1-mask_u(i,j,k1))*veldydz_u(i,j,k,1)
           *mask_u(i,j,k )+(1-mask_u(i,j,k ))*veldydz_u(i,j,k,1)
        
       enddo
       enddo
      enddo


! composante V:
      k0=0
      do k=1,kmax
       do j=1,jmax+1
       do i=1,imax
        xy_v(i,j,1)=xy_v(i,j,1)*k0                &
       +dz_v(i,j,k,1)*(  0.5*( vel_v(i,j,k,1)     &
                              +vel_v(i,j,k,2))    &
                        +velstokes_v(i,j,k,1)  )  &
                             *mask_v(i,j,k)
       enddo
       enddo
      k0=1
      enddo

      const2=1.0/real(iteration2d_max_now+iteration2d_max_bef)    !19-04-11
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,2)=( & !gggg>
                      const2*(fluxbar_sumt_v(i,j,0)   &
                             +fluxbar_sumt_v(i,j,1))  &
                                       /dx_v(i,j)     &

                                      -xy_v(i,j,1)    &
                   ) & !gggg>
                          /hz_v(i,j,1)    
      enddo
      enddo

      do k=1,kmax
! si method vst discontinue k1=k sinon (sigma ou vst continue) k1=kmax
!     k1=vststep*k+(1-vststep)*kmax !05-08-15
       do j=1,jmax+1
       do i=1,imax

         veldxdz_v(i,j,k,1)=( & !rrrr>
       0.5*( vel_v(i,j,k,1)           &
            +vel_v(i,j,k,2))          &
             +xy_v(i,j,2)             &
      +velstokes_v(i,j,k,1) ) & !rrrr>
             *dx_v(i,j)       &
             *dz_v(i,j,k,1)   &
        *wetmask_v(i,j)       & !25-02-19 wetdry schema 6
!          *mask_v(i,j,k1)+(1-mask_v(i,j,k1))*veldxdz_v(i,j,k,1)
           *mask_v(i,j,k )+(1-mask_v(i,j,k ))*veldxdz_v(i,j,k,1)

       enddo
       enddo
      enddo


      if(flag_merged_levels==1) then !- REDUCTEUR DE FLUX - COUCHE KMERGED-1 ---> !02-11-21

      do j=1,jmax ; do i=1,imax+1
!     if(kmin_u(i,j)/=kmerged_u(i,j)) then !m°v°m>

       sum1=0. ; sum2=0.
       do k=1,kmerged_u(i,j)-2 
        sum1=sum1+veldydz_u(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite
        veldydz_u(i,j,k,1)=0.        ! nouveau flux
       enddo

       k=kmerged_u(i,j)-1
       sum1=sum1+veldydz_u(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite

      if(kmin_u(i,j)/=kmerged_u(i,j)) then !m°v°m> !19-03-22
! Formule standard mais dz est remplacE par h*dsigmerged
         veldydz_u(i,j,k,1)=( & !rrrr> ! nouveau flux
       0.5*( vel_u(i,j,k,1)           &
            +vel_u(i,j,k,2))          &
             +xy_u(i,j,2)             &
      +velstokes_u(i,j,k,1) ) & !rrrr>
             *dy_u(i,j)       &
             *dz_u(i,j,k,1)   &    
     *dsigmerged_u(i,j)       & !17-11-21
!        *min(hz_u(i,j,1)*dsigmerged_u(i,j),dz_u(i,j,k,1))  &!14-11-21
        *wetmask_u(i,j)       & !25-02-19 wetdry schema 6
           *mask_u(i,j,k )+(1-mask_u(i,j,k ))*veldydz_u(i,j,k,1)
      else                                 !m°v°m>
         veldydz_u(i,j,k,1)=0.   !19-03-22
      endif                                !m°v°m>
      sum2=sum2+veldydz_u(i,j,k,1) ! Erreur de bilan de tranport du nouveau flux A corriger ensuite

! La couche kmerged_u recupere le defaut de bilan de transport sum2-sum1
       k=kmerged_u(i,j)
         veldydz_u(i,j,k,1)= &
         veldydz_u(i,j,k,1)+sum1-sum2

!     endif                                !m°v°m>
      enddo ; enddo

      do j=1,jmax+1 ; do i=1,imax
!     if(kmin_v(i,j)/=kmerged_v(i,j)) then !m°v°m>

       sum1=0. ; sum2=0.
       do k=1,kmerged_v(i,j)-2 
        sum1=sum1+veldxdz_v(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite
        veldxdz_v(i,j,k,1)=0.        ! nouveau flux
       enddo

       k=kmerged_v(i,j)-1
       sum1=sum1+veldxdz_v(i,j,k,1) ! Bilan de transport de l'ancien flux A recuperer ensuite

      if(kmin_v(i,j)/=kmerged_v(i,j)) then !m°v°m> !19-03-22
! Formule standard mais dz est remplacE par h*dsigmerged
         veldxdz_v(i,j,k,1)=( & !rrrr> ! nouveau flux
       0.5*( vel_v(i,j,k,1)           &
            +vel_v(i,j,k,2))          &
             +xy_v(i,j,2)             &
      +velstokes_v(i,j,k,1) ) & !rrrr>
             *dy_v(i,j)       &
             *dz_v(i,j,k,1)   &    
     *dsigmerged_v(i,j)       & !17-11-21
!        *min(hz_v(i,j,1)*dsigmerged_v(i,j),dz_v(i,j,k,1))  & !14-11-21
        *wetmask_v(i,j)       & !25-02-19 wetdry schema 6
           *mask_v(i,j,k )+(1-mask_v(i,j,k ))*veldxdz_v(i,j,k,1)
      else                                 !m°v°m>
         veldxdz_v(i,j,k,1)=0.    !19-03-22
      endif                                !m°v°m>
      sum2=sum2+veldxdz_v(i,j,k,1) ! Erreur de bilan de tranport du nouveau flux A corriger ensuite

! La couche kmerged_v recupere le defaut de bilan de transport sum2-sum1
       k=kmerged_v(i,j)
         veldxdz_v(i,j,k,1)= &
         veldxdz_v(i,j,k,1)+sum1-sum2

!     endif                                !m°v°m>
      enddo ; enddo

      endif                          !- REDUCTEUR DE FLUX - COUCHE KMERGED-1 --->

      call omega !15-12-16

      end subroutine couple_modes_scalars

!.......................................................................

      subroutine couple_modes_zprofil_instant !29-09-16
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='couple_modes_zprofil_instant'
       subroutinedescription=' '
#endif
! Mes notes:
! https://docs.google.com/document/d/1kXMG-3WaF1Miif5-gj1kzN8dSpD5Mhhowg0IzoM1Afk/edit

!     inv_ekman_depth=1./10. ! d**-1=sqrt(f/(2.Kz)) ! sqrt(0.5*10e-4/5e-2)
!     inv_ekman_depth=1./10. ! 1/(Ekman depth)

      if(ihybsig/=0)stop 'couple_modes_instant ihybsig/=0'

!.....................................................
! Step 1: depth-averaged current before correction
! reset en K=1:
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=vel_u(i,j,1,2)*dz_u(i,j,1,2)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=vel_v(i,j,1,2)*dz_v(i,j,1,2)
      enddo
      enddo

! sommation sur K suivants:
      do k=2,kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=xy_u(i,j,1)          &
                 +vel_u(i,j,k,2)      &
                  *dz_u(i,j,k,2)
      enddo
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax+1
      do i=1,imax
         xy_v(i,j,1)=xy_v(i,j,1)          &
                   +vel_v(i,j,k,2)      &
                    *dz_v(i,j,k,2)
      enddo
      enddo
      enddo

      do j=1,jmax !07-04-14
      do i=1,imax+1
         xy_u(i,j,1)=xy_u(i,j,1)  & 
                    /hz_u(i,j,2)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
         xy_v(i,j,1)=xy_v(i,j,1)  & 
                    /hz_v(i,j,2)
      enddo
      enddo


!............................................................
! Step2 : Vertical integral of the z profile of the corretion:
!.............
! U component:
      k=kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,3)=                                             &
                      dz_u(i,j,k,2)                             &
!       *(1.-exp(-(depth_u(i,j,k)-depth_u(i,j,kmin_u(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_u(i,j,k)+h_u(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_u(i,j,k)                               &
                +(1-mask_u(i,j,k))*small1

      enddo
      enddo
      do k=kmax-1,1,-1
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,3)=                                             &
       xy_u(i,j,3)+                                             &
                      dz_u(i,j,k,2)                             &
!       *(1.-exp(-(depth_u(i,j,k)-depth_u(i,j,kmin_u(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_u(i,j,k)+h_u(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_u(i,j,k)
      enddo
      enddo
      enddo

      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,4)=    (-xy_u(i,j,1)     &
                    +velbar_u(i,j,2))    &
                        *hz_u(i,j,2)     &
                        /xy_u(i,j,3)
      enddo
      enddo

!............................................................
! Step3 : Update U
!     sum1=0. 
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1

!      x1=vel_u(i,j,k,2) !05-02-17 

                   vel_u(i,j,k,2)=                            &
                   vel_u(i,j,k,2)                             &
                   +xy_u(i,j,4)                               &
!     *(1.-exp(-(depth_u(i,j,k)-depth_u(i,j,kmin_u(i,j)-1))*inv_ekman_depth )) & ! methode 1
      *(1.-exp(-(depth_u(i,j,k)+h_u(i,j))*inv_ekman_depth )) & ! methode 2
                 *mask_u(i,j,k)

      enddo
      enddo
      enddo

!............................................................
! Step2 : Vertical integral of the z profile of the corretion:
!.............
! V component:
      k=kmax
      do j=1,jmax+1
      do i=1,imax
! Note on suppose que le courant est nul en z=-h_v ce qui
! n'est pas toujours vrai en coordonnee hybride mais bon....
       xy_v(i,j,3)=                                             &
                      dz_v(i,j,k,2)                             &
!       *(1.-exp(-(depth_v(i,j,k)-depth_v(i,j,kmin_v(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_v(i,j,k)+h_v(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_v(i,j,k)                               &
                +(1-mask_v(i,j,k))*small1
      enddo
      enddo
      do k=kmax-1,1,-1
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,3)=                                             &
       xy_v(i,j,3)+                                             &
                      dz_v(i,j,k,2)                             &
!       *(1.-exp(-(depth_v(i,j,k)-depth_v(i,j,kmin_v(i,j)-1))*inv_ekman_depth )) & ! methode 1
        *(1.-exp(-(depth_v(i,j,k)+h_v(i,j))*inv_ekman_depth )) & ! methode 2
                   *mask_v(i,j,k)
      enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,4)=(    -xy_v(i,j,1)    &
                    +velbar_v(i,j,2) )  &
                        *hz_v(i,j,2)    &
                        /xy_v(i,j,3)
      enddo
      enddo

!............................................................
! Step3 : Update V
!     sum1=0.
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
                   vel_v(i,j,k,2)=                            &
                   vel_v(i,j,k,2)                             &
                   +xy_v(i,j,4)                               &
!     *(1.-exp(-(depth_v(i,j,k)-depth_v(i,j,kmin_v(i,j)-1))*inv_ekman_depth )) & ! methode 1
      *(1.-exp(-(depth_v(i,j,k)+h_v(i,j))*inv_ekman_depth )) & ! methode 2
                 *mask_v(i,j,k)

!      if(i==imax/2.and.j==jmax/2)sum1=sum1+vel_v(i,j,k,2)*dz_v(i,j,k,2)
!      if(i==imax/2.and.j==jmax/2.and.k==kmax.and.mask_v(i,j,kmax)==1)write(10+par%rank,*)sum1/hz_v(i,j,2),velbar_v(i,j,2)
      enddo
      enddo
      enddo

      end subroutine couple_modes_zprofil_instant
