      subroutine pressure_gradient
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 - last update: 07-03-23
!______________________________________________________________________

      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient'
       subroutinedescription= &
       'Driver of the pressure gradient subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         06/03/03: passage à la coordonnee verticale hybride
!                   modifs sur calcul intermediare: stockage dans X2 plutot que
!                   dans XY_Z(I,J,2)
!         09/05/07: Choix entre un schema d'ordre 1 et un schema d'ordre 3
!         27/12/07: Schema ordre 1 retenu. Le terme de correction lie à la
!                   pente des sigma a une nouvelle ponderation deduite du
!                   nouveau concept d'initialisation plus proche de la definition
!                   d'une variable en volume fini: rho(i,j,k)=<rho>
!         24/02/08: Ajout de conditions aux limites pour la pression barocline
! 2009.3  30-09-09: On ne retient qu'un seul schema
!                   anyv3d(i,j,k,2) remplace pressure_z(i,j,k)
!                   coupure du pgf en cas de model_ 1D (*flag3d=0) (déplacé de uv_upd.F)
!         02-10-09: utilisation des nouveaux facteurs d'echelle verticale
!         06-12-09: numerotation generique
! 2010.11 25-07-10  anyv3d(i,j,k,0) est la densite integrant l'effet de compression
!         08-08-10  3 schémas de PGF selon le type d'equation d'état
!         09-08-10  - time_   =now
!                   - alpha=0 si k=kmax
!         10-08-10  schema pom
! 2010.14 21-11-10  Orientation vers le pgf_eos1 si effet de compression dans
!                   l'equation d'état
!         16-12-10  - Ajout cas EOS Wright (1997) et EOS McDougall et al. (2003)
!                   - Introduction de l'appel à l'EOS au cas par cas
! 2010.15 09-01-10  - Suite du point precedent: correction d'un bug dans
!                   le calcul de z*
!                   - Pour chaque pgf on a codé une version hybride où
!                   la partie compressible du pgf s'appuie sur la
!                   formulation géopotentielle du schema primitif direct
! 2010.25 22-02-12  Le gradient de sshstokes est maintenant pris en compte
!                   dans la routine stokesforces
!         03-08-12  L'ajout des termes 2d meteo & maree ainsi que la
!                   moyenne verticale pour le mode externe sont calcules dans la
!                   subroutine pressure_gradient_2d_terms
! S.26    16-06-14  PGF marsaleix et al 2001: pression de surface identique a droite
!                   et a gauche
!         19-06-14  Version compactee du PGF marsaleix et al 2001. Possibilite
!                   d'appel d'une densite potentielle a z reference locale pour le
!                   schema de tke
!         20-06-14  call equation_of_state_potzref2
!         20-10-14  legeres modifs ordre d'appel du synopsis
!         10-07-15  rhpzavr est une moyenne verticale ponderee de rgp dependant de x,y
!         14-09-15  pressure_gradient_initref permet de prendre en compte un etat de
!                   reference
!         16-09-15  Adaptation de la discretisation de rho moyen dependant de (x,y) au
!                   nouvel algo de l'etat de reference
!         15-10-15  PGF Marsaleix et al 2009 + etat de reference
!         31-10-15  suite du point precedent. Plus stricte dans les zones de bancs decouvrants
!         14-11-15  une nouvelle version de la routine "homogeneous"
!         18-11-15  Modif pour que rhcref soit calcule meme si flag_refstate==0 
!         29-11-15  Cas 1DV
!         30-11-15  Cas 1DV interpolation temporelle de velobc(0) et velobc(2)
!         28-01-16  nouveau time-stepping LF-FB
!         16-02-16  aiguillage 3D 1D dans pressure_gradient_eos0
!         14-04-16  call modeanalysis_pmodeprojection 
!         24-04-16  Systematiquement calculer la densite potentielle maintenant que celle-ci 
!                   est egalement utilisee dans le calcul de omega_w(:,:,kmax+1) et le calcul 
!                   de rho.cp dans vertmix.F90
!         25-04-16  Suite du point precedent
!         15-09-16  methode equivalente mais + rapide pour calcul rhpzavr_w
!         08-04-17  s-z coordinate
!         21-04-17  s-z coordinate suite
!         17-10-17  suite du poin precedent
!         20-08-18  gradient de pression au temps now
!         23-08-18  if(iteration2d_max_now>1)call pressure_gradient_2d_terms 
!         06-09-18  condition d'ajout de nhpgf
! v285    09-06-20  - if(flag_merged_levels==1)call pressure_gradient_vqs
!                   - Profil lineaire reconstituE dans la couche fusionnEe
!                     par turbulence_rhp_linear_profil
! v293    08-10-20  ajout pss_w dans "cas sans split sans T,S" de pressure_gradient_add_nhpgf
! v310    29-10-21  - reconstruction ordre 2 du profil de densite/pression A la
!                     profondeur intermediaire virtuelle du PGF VQS
!         02-11-21  - par contre on ne reconstruit pas un profil lineaire sur les niveaux 
!                     kmin_w et kmerged_t si ces niveaux ne sont plus melangEs
!         03-11-21  - retenir 100% du calcul PGF(kmerged-1) si dz(kmerged-1)>=dz(kmerged)
!                     et prendre la valeur de PGF(kmerged) si dz(kmerged-1)=0
!                   - dz methode 2
! v339    24-03-22  l'analyse harmonique porte sur rho_t la densitE utilisEe par le gradient 
!                   de pression
! v365    26-01-23  modif rho_t
! v367    07-03-23  activer calcul rho_t
!.......................................................................
!    _________                    .__                  .__             ! 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !

! Computes the Pressure Gradient Force related to the sea water density.
! The scheme depends on the type of equation of state.

       if(timestep_type==timestep_forwbckw) then !fbfbfb> !28-01-16
        presgrad_u(:,:,:,0)=presgrad_u(:,:,:,1)
        presgrad_v(:,:,:,0)=presgrad_v(:,:,:,1)
       endif                                     !fbfbfb>

       if(flag_merged_levels==0)call pressure_gradient_eos0
!      if(flag_merged_levels==1)call pressure_gradient_eos1
       if(flag_merged_levels==1)call pressure_gradient_vqs !09-06-20

       if(iteration3d==0.or.timestep_type==timestep_leapfrog) then !lplplp>
! A lEtat initial pour les 2 methodes de time-stepping ou bien pour le
! schema Leap-frog (uniquement):
        presgrad_u(:,:,:,0)=presgrad_u(:,:,:,1)
        presgrad_v(:,:,:,0)=presgrad_v(:,:,:,1)
       endif                                                       !lplplp>

! Tide potential + atmospheric pressure gradient, depth-averaged pgf !23-08-18
                                 call pressure_gradient_add_2d_terms  ! Tide potential + atmospheric pressure graident
!      if(iteration2d_max_now/=0)call pressure_gradient_zaveragedpgf  ! z-averaged pressure gradient
       if(flag_timesplitting==1) call pressure_gradient_zaveragedpgf  ! z-averaged pressure gradient

! In case of no compression effects in the equation of state (eos_comprs=0):
#ifdef bidon
      if(eos_comprs==0) then
       call pressure_gradient_eos0 ! No compression effects in the equation of state
       call pressure_gradient_2d_terms
       return
      else
        if(eos_author==0) then
           call pressure_gradient_eos0       ! Linear EOS
           call pressure_gradient_2d_terms
           return
        endif
        if(eos_author==1) then
!          call pressure_gradient_m91                     ! Mellor (1991) EOS
           call pressure_gradient_m91_hyb !09-01-11       ! Mellor (1991) EOS
           call pressure_gradient_2d_terms
           return
        endif
        if(eos_author==2) then
!          call pressure_gradient_w97                 ! Wright (1997) EOS
           call pressure_gradient_w97_hyb  !09-01-11  ! Wright (1997) EOS
           call pressure_gradient_2d_terms
           return
        endif
        if(eos_author==3) then
!          call pressure_gradient_jmfwg06                 ! Jackett et al (2006) EOS
           call pressure_gradient_jmfwg06_hyb !09-01-11   ! Jackett et al (2006) EOS
           call pressure_gradient_2d_terms
           return
        endif
      endif
      stop'No PGF selected'
#endif

      end subroutine pressure_gradient

!.................................................................................

      subroutine pressure_gradient_add_2d_terms !03-08-12
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_2d_terms'
       subroutinedescription= &
       'Depth-averaged residual of the sea-water-dependent 3D pressure'&
       //' gradient that will be used as a frozen forcing term in the' &
       //' barotropic momentum equations.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! This subroutine adds the atmospheric and tidal 2d Pressure Gradient
! Forces (PGF) to the "sea-water-density-PGF" and then computes the
! depth-averaged PGF residual term used by the external mode momentum
! equations.

! U component:
! GEneric MOdelling shortname conventions (ipu=0 imu=1)

      do j=2,jmax-1
      do i=2,imax

! 2d atmospheric and tidal terms:
            xy_u(i,j,1)=                                              &
        (tidepotential_w(i-imu,j,1)-tidepotential_w(i+ipu,j,1))*grav  & ! Tidal  equilibrium SSH gradient
       +(          pss_w(i+ipu,j,1)          -pss_w(i-imu,j,1))/rho     ! Atmos. pressure force

      enddo
      enddo

      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax

! Add the 2d atmospheric and tidal terms:
      presgrad_u(i,j,k,1)=presgrad_u(i,j,k,1)+xy_u(i,j,1)

      enddo
      enddo
      enddo

! V component:
! GEneric MOdelling shortname conventions (jpv=0 jmv=1)

      do j=2,jmax
      do i=2,imax-1

! 2d atmospheric and tidal terms:
            xy_v(i,j,1)=                                              &
        (tidepotential_w(i,j-jmv,1)-tidepotential_w(i,j+jpv,1))*grav  & ! Tidal  equilibrium SSH gradient
       +(          pss_w(i,j+jpv,1)          -pss_w(i,j-jmv,1))/rho     ! Atmos. pressure force

      enddo
      enddo

      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1

! Add the 2d atmospheric and tidal terms:
      presgrad_v(i,j,k,1)=presgrad_v(i,j,k,1)+xy_v(i,j,1)

      enddo
      enddo
      enddo

      end subroutine pressure_gradient_add_2d_terms

!.................................................................................

      subroutine pressure_gradient_zaveragedpgf
      use module_principal
      use module_parallele
      implicit none

! This subroutine computes the
! depth-averaged PGF residual term used by the external mode momentum
! equations.

! U component:
! GEneric MOdelling shortname conventions (ipu=0 imu=1)
      do j=2,jmax-1 ; do i=2,imax
! Reset the Depth-averaged PGF frozen term for external mode U-momentum equation:
       pres3d2d_u(i,j)=0.
      enddo         ; enddo

      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
! Compute the depth-averaged PGF frozen term for external mode U-momentum equation:
      pres3d2d_u(i,j)=pres3d2d_u(i,j)+dz_u(i,j,k,1)*presgrad_u(i,j,k,1)  &
                                                       *mask_u(i,j,k)  !06/03/03
      enddo ; enddo ; enddo

      do j=2,jmax-1 ; do i=2,imax
        pres3d2d_u(i,j)                                                 &
       =pres3d2d_u(i,j)                                                 &
             /hz_u(i,j,1)
      enddo ; enddo

! V component:
! GEneric MOdelling shortname conventions (jpv=0 jmv=1)
      do j=2,jmax ; do i=2,imax-1
! Reset the depth-averaged PGF frozen term for external mode V-momentum equation:
       pres3d2d_v(i,j)=0.
      enddo ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
! Compute the depth-averaged PGF frozen term for external mode V-momentum equation:
      pres3d2d_v(i,j)=pres3d2d_v(i,j)+dz_v(i,j,k,1)*presgrad_v(i,j,k,1)  &
                                                       *mask_v(i,j,k)  !06/03/03
      enddo ; enddo ; enddo

      do j=2,jmax ; do i=2,imax-1
        pres3d2d_v(i,j)                                                 &
       =pres3d2d_v(i,j)                                                 &
             /hz_v(i,j,1)
      enddo ; enddo

      end subroutine pressure_gradient_zaveragedpgf

!...............................................................................

      subroutine pressure_gradient_eos2
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_eos2'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!$ This routine computes the pressure gradient force related to the potential density
!$ according to the method described in:
!$ Marsaleix P., Auclair F., Estournel C., 2009. Low-order pressure gradient
!$ schemes in sigma coordinate models: The seamount test revisited.
!$ Ocean Modelling, 30, 169-177. http://dx.doi.org/10.1016/j.ocemod.2009.06.011
!$ The pressure gradient force related to the compression terms of the state equation
!$ is computed with a trapezoidal method

! Compute anyv3d(i,j,k,0) the part of the density related to pressure:
      call equation_of_state('compression terms',now)

!  Computes the Hydrostatic pressure:
!  - using the rectangular method as far as the potential part of the density is concerned
!  - using the trapezoidal method as far as the "pressured" part of the density is concerned

! Note: rhp_t, the potential density - rho, has been computed at the previous step in scalars.F

! Reset of the summation at sea surface:
! At sea surface z=depth_w(kmax+1) and compression vanishes (i.e. anyv3d(i,j,k,0)=0)
      k=kmax
      const3=0.5*grav
      do j=1,jmax
      do i=1,imax

       xy_t(i,j,1)=grav*rhp_t(i,j,k)*dz_t(i,j,k,1)             & ! potential part (rectangular)
                  +const3*                   anyv3d(i,j,k,0)   & ! compression part (trapezoidal)
                        *(depth_w(i,j,k+1) -depth_t(i,j,k))

! Hydrostatic pressure at "t" nodes & at k=kmax
       anyv3d(i,j,k,2)=xy_t(i,j,1)-const3*rhp_t(i,j,k)*dz_t(i,j,k,1)

      enddo
      enddo

! Summation over the kmax-1 remaining levels:
      const3=0.5*grav
      do 10 k=kmax-1,1,-1
       do 10 j=1,jmax
       do 10 i=1,imax

       xy_t(i,j,1)=xy_t(i,j,1)                                 &
                  +grav*rhp_t(i,j,k)*dz_t(i,j,k,1)             & ! potential part (rectangular)
                  +const3*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0))  & ! compression part (trapezoidal)
                        *(depth_t(i,j,k+1) -depth_t(i,j,k))

!$ Hydrostatic pressure at "t" nodes:
       anyv3d(i,j,k,2)=xy_t(i,j,1)-const3*rhp_t(i,j,k)*dz_t(i,j,k,1)

   10 continue

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

      do 14 k=kmax,1,-1
      do 14 j=2,jmax-1
      do 14 i=2,imax

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                              &
                   +const1*( anyv3d(i+ip,j,k,2)-anyv3d(i-im,j,k,2))   &
                   +const2*(depth_t(i+ip,j,k) -depth_t(i-im,j,k))*(   &
                           (  rhp_t(i+ip,j,k)  *  dz_t(i+ip,j,k,1)    & ! The weighted average for potential density
                             +rhp_t(i-im,j,k)  *  dz_t(i-im,j,k,1))   & ! is consistent with the "rectangular" pressure
                           /(  dz_t(i+ip,j,k,1) + dz_t(i-im,j,k,1))   &
                       +0.5*(anyv3d(i+ip,j,k,0)+anyv3d(i-im,j,k,0)))    ! The half half average for "pressured" density
                                                                        ! is consistent with the "trapezoidal" pressure
   14 continue

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09
      const2= grav/(rho)*flag3d                                            !30-09-09

      do 19 k=kmax,1,-1
      do 19 j=2,jmax
      do 19 i=2,imax-1

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                              &
                   +const1*( anyv3d(i,j+jp,k,2)-anyv3d(i,j-jm,k,2))   &
                   +const2*(depth_t(i,j+jp,k) -depth_t(i,j-jm,k))*(   &
                           (  rhp_t(i,j+jp,k)  *  dz_t(i,j+jp,k,1)    & !27/12/07
                             +rhp_t(i,j-jm,k)  *  dz_t(i,j-jm,k,1))   &
                           /(  dz_t(i,j+jp,k,1)+  dz_t(i,j-jm,k,1))   &
                       +0.5*(anyv3d(i,j+jp,k,0)+anyv3d(i,j-jm,k,0)))

   19 continue

      end subroutine pressure_gradient_eos2

!..............................................................................

      subroutine pressure_gradient_m91
      use module_principal
      implicit none
      integer time_
#ifdef synopsis
       subroutinetitle='pressure_gradient_m91'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      time_   =now
      call equation_of_state('potential density',now)

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=0.
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference:
!     anyv3d(i,j,k,1)=depth_t(i-im,j,k)                    &
!                -2.*(depth_t(i-im,j,k)-depth_u(i,j,k))    &
!                       *dz_t(i-im,j,k,1)                  &
!                      /(dz_t(i-im,j,k,1)                  &
!                       +dz_t(i+ip,j,k,1))
! same as:
!     anyv3d(i,j,k,1)=depth_t(i+ip,j,k)                    &
!                +2.*(depth_t(i-im,j,k)-depth_u(i,j,k))    &
!                       *dz_t(i+ip,j,k,1)                  &
!                      /(dz_t(i-im,j,k,1)                  &
!                       +dz_t(i+ip,j,k,1))
! Same as:
      anyv3d(i,j,k,1)=(depth_t(i-im,j,k)*dz_t(i+ip,j,k,1)  &
                      +depth_t(i+ip,j,k)*dz_t(i-im,j,k,1)) &
                     /(                  dz_t(i+ip,j,k,1)  &
                                        +dz_t(i-im,j,k,1))

! Hydrostatic pressure at i=i-im:

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.25*depth_w(i-im,j,k+1)    &
        +0.25*depth_w(i+ip,j,k+1)    &
         +0.5*anyv3d(i,j,k,1)

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i-im,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i-im,j,1)=                                                &
                     grav*(depth_w(i-im,j,k+1)                       &
                          -anyv3d(i,j,k  ,1))                        &
                                *(     alpha *rhp_t(i-im,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i-im,j,k)        &
                                     +2000.*x0*(5.-x0) )

! Hydrostatic pressure at i=i+ip:
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i+ip,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i+ip,j,2)=                                                &
                     grav*(depth_w(i+ip,j,k+1)                       &
                          -anyv3d(i,j,k  ,1))                        &
                                *(     alpha *rhp_t(i+ip,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i+ip,j,k)        &
                                     +2000.*x0*(5.-x0) )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &
                +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1

      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference:
!     anyv3d(i,j,k,1)=depth_t(i-im,j,k)                    &
!                -2.*(depth_t(i-im,j,k)-depth_u(i,j,k))    &
!                       *dz_t(i-im,j,k,1)                  &
!                      /(dz_t(i-im,j,k,1)                  &
!                       +dz_t(i+ip,j,k,1))
! same as:
!     anyv3d(i,j,k,1)=depth_t(i+ip,j,k)                    &
!                +2.*(depth_t(i-im,j,k)-depth_u(i,j,k))    &
!                       *dz_t(i+ip,j,k,1)                  &
!                      /(dz_t(i-im,j,k,1)                  &
!                       +dz_t(i+ip,j,k,1))
! Same as:
      anyv3d(i,j,k,1)=(depth_t(i-im,j,k)*dz_t(i+ip,j,k,1)  &
                      +depth_t(i+ip,j,k)*dz_t(i-im,j,k,1)) &
                     /(                  dz_t(i+ip,j,k,1)  &
                                        +dz_t(i-im,j,k,1))

! Hydrostatic pressure at i=i-im:

! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=(anyv3d(i,j,k+1,1)-depth_w(i-im,j,k+1))    &  !09-08-10
           /(anyv3d(i,j,k+1,1)- anyv3d(i   ,j,k,1))

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.5*(anyv3d(i,j,k+1,1)+anyv3d(i,j,k,1))
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i-im,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i-im,j,1)=                                                &
      xy_t(i-im,j,1)+grav*(anyv3d(i,j,k+1,1)                         &
                          -anyv3d(i,j,k  ,1))                        &
                                *(     alpha *rhp_t(i-im,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i-im,j,k)        &
                                     +2000.*x0*(5.-x0) )

! Hydrostatic pressure at i=i+ip:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,1)-depth_w(i+ip,j,k+1))          &
           /(anyv3d(i,j,k+1,1)- anyv3d(i   ,j,k,1))

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i+ip,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i+ip,j,2)=                                                &
      xy_t(i+ip,j,2)+grav*(anyv3d(i,j,k+1,1)                         &
                          -anyv3d(i,j,k  ,1))                        &
                                *(     alpha *rhp_t(i+ip,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i+ip,j,k)        &
                                     +2000.*x0*(5.-x0) )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &
                +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1

      enddo
      enddo
      enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging:
      alpha=0.
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference:
!     anyv3d(i,j,k,1)=depth_t(i,j-jm,k)                    &
!                -2.*(depth_t(i,j-jm,k)-depth_v(i,j,k))    &
!                       *dz_t(i,j-jm,k,1)                  &
!                      /(dz_t(i,j-jm,k,1)                  &
!                       +dz_t(i,j+jp,k,1))
! same as
      anyv3d(i,j,k,1)=(depth_t(i,j-jm,k)*dz_t(i,j+jp,k,1)   &
                      +depth_t(i,j+jp,k)*dz_t(i,j-jm,k,1))  &
                     /(                  dz_t(i,j+jp,k,1)   &
                                        +dz_t(i,j-jm,k,1))

! Hydrostatic pressure at j=j-jm:

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.25*depth_w(i,j-jm,k+1)    &
        +0.25*depth_w(i,j+jp,k+1)    &
          +0.5*anyv3d(i,j,k,1)

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j-jm,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j-jm,1)=                                                &
                     grav*(depth_w(i,j-jm,k+1)                       &
                           -anyv3d(i,j   ,k  ,1))                    &
                                *(     alpha *rhp_t(i,j-jm,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j-jm,k)        &
                                     +2000.*x0*(5.-x0) )

! Hydrostatic pressure at j=j+jp:
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j+jp,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j+jp,2)=                                                &
                     grav*(depth_w(i,j+jp,k+1)                       &
                           -anyv3d(i,j   ,k  ,1))                    &
                                *(     alpha *rhp_t(i,j+jp,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j+jp,k)        &
                                     +2000.*x0*(5.-x0) )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                      &
                +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1


      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference:
!     anyv3d(i,j,k,1)=depth_t(i,j-jm,k)                    &
!                -2.*(depth_t(i,j-jm,k)-depth_v(i,j,k))    &
!                       *dz_t(i,j-jm,k,1)                  &
!                      /(dz_t(i,j-jm,k,1)                  &
!                       +dz_t(i,j+jp,k,1))
! same as:
!     anyv3d(i,j,k,1)=depth_t(i,j+jp,k)                    &
!                +2.*(depth_t(i,j-jm,k)-depth_v(i,j,k))    &
!                       *dz_t(i,j+jp,k,1)                  &
!                      /(dz_t(i,j-jm,k,1)                  &
!                       +dz_t(i,j+jp,k,1))
! same as
      anyv3d(i,j,k,1)=(depth_t(i,j-jm,k)*dz_t(i,j+jp,k,1)   &
                      +depth_t(i,j+jp,k)*dz_t(i,j-jm,k,1))  &
                     /(                  dz_t(i,j+jp,k,1)   &
                                        +dz_t(i,j-jm,k,1))

! Hydrostatic pressure at j=j-jm:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,1)-depth_w(i,j-jm,k+1))          &
           /(anyv3d(i,j,k+1,1)- anyv3d(i,j,k,1))

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.5*(anyv3d(i,j,k+1,1)+anyv3d(i,j,k,1))
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j-jm,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j-jm,1)=                                                &
      xy_t(i,j-jm,1)+grav*(anyv3d(i,j,k+1,1)                         &
                          -anyv3d(i,j,k  ,1))                        &
                                *(     alpha *rhp_t(i,j-jm,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j-jm,k)        &
                                     +2000.*x0*(5.-x0) )

! Hydrostatic pressure at j=j+jp:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,1)-depth_w(i,j+jp,k+1))          &
           /(anyv3d(i,j,k+1,1)- anyv3d(i,j,k,1))

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j+jp,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j+jp,2)=                                                &
      xy_t(i,j+jp,2)+grav*(anyv3d(i,j,k+1,1)                         &
                          -anyv3d(i,j,k  ,1))                        &
                                *(     alpha *rhp_t(i,j+jp,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j+jp,k)        &
                                     +2000.*x0*(5.-x0) )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                      &
                +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1


      enddo
      enddo
      enddo

      end subroutine pressure_gradient_m91

!..............................................................................

      subroutine pressure_gradient_eos1
      use module_principal ; use module_parallele ; use module_modeanalysis
      implicit none
      integer :: id_prsref_=5

! ATTENTION le rho1 dependant de (x,y) n'est pas utilisE comme dans:
!https://docs.google.com/document/d/1_SdjMHFzI2EVxk17kXwjsXi4ra1nHiU_Wyg2LyL9KZ8/edit?usp=sharing
! AUTREMENT DIT IL N'EST PAS COMPATIBLE AVEC L'ETAT DE REFERENCE 3D (amenager ce dernier ou
! alors utiliser pressure_gradient_eos0)

      if(flag_refstate==1) &
      stop 'Err flag_refstate=1 dans pressure_gradient_eos1'


! Compute rhp_t(:,:,:)
      if(eos_comprs==0) then !000000>
! Compute rhp the potential density at z=0m and the a1, a2, ....,a7 coefficients of the Jackett et al 2006 (JAOT) EOS
       if(eos_author==0)call equation_of_state_linear(now)
       if(eos_author==3)call equation_of_state_rhp_a1_a7(now) ! rhp=densite potentielle (zref=0)
       if(eos_author==2)stop 'Err 697 pressure_gradient'
       if(eos_author==1)stop 'Err 698 pressure_gradient'
       if(flag_tide3d_analysis==1 &
             .or.drifter_onoff==1) then !---> !07-03-23

        rho_t(0:imax+1,0:jmax+1,1:kmax)=rhp_t(0:imax+1,0:jmax+1,1:kmax)+rho !26-01-23
       endif                            !--->
      endif                  !000000>

      if(eos_comprs==1) then !111111>
       if(eos_author==0)stop 'Err 699 pressure_gradient'
       if(eos_author==3)call equation_of_state_full_a1_a7(now) ! rhp=densite potentielle (zref=0)
       if(eos_author==2)stop 'Err 680 pressure_gradient'
       if(eos_author==1)stop 'Err 681 pressure_gradient'
       if(flag_tide3d_analysis==1 &
             .or.drifter_onoff==1) then !---> !07-03-23

        rho_t(0:imax+1,0:jmax+1,1:kmax)=rhp_t(0:imax+1,0:jmax+1,1:kmax)+rho !26-01-23
       endif                            !--->
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhp_t(i,j,k)=rhp_t(i,j,k)-rhcref_t(i,j,k) ! Ote l'etat de reference de compression
       enddo       ; enddo       ; enddo
      endif                  !111111>

! En cas d'analyse harmonique de la pression A partir de rhp avant que rhmean ne soit otE:
      if(flag_3dwaves_pharmonics==1)call modeanalysis_pmodeprojection !14-04-16

! remove a depth-weighted average of rhp
      if(rhp_zavr_xy==1)call pressure_gradient_rhpzavr

!$ Computes the Hydrostatic pressure using the rectangular method:
      id_prs=0 ! attention anyv3d(:,:,:,1 a 4 et de 6 a 7) est pris par les coef de l'EOS

! verifier la disponibilite de anyv3d(:,:,:,id_prs)
      if(anyv3d(-1,-1,0,id_prs)==-9999.) &
      stop 'Err anyv3d id_prs not available'

      if(flag_1dv==0) then !3D3D3D3D3D> !16-02-16

      do j=1,jmax ; do i=1,imax
       anyv3d(i,j,kmax+1,id_prs)=0.
       anyv3d(i,j,kmax+1,id_prsref_)=0.
      enddo       ; enddo
      do k=kmax,1,-1 ; do j=1,jmax ; do i=1,imax
             anyv3d(i,j,k  ,id_prs)= &
             anyv3d(i,j,k+1,id_prs)  &
        +grav*rhp_t(i,j,k)           &
              *dz_t(i,j,k,1)

             anyv3d(i,j,k  ,id_prsref_)= &
             anyv3d(i,j,k+1,id_prsref_)  &
       +grav*(rhp_t(i,j,k)               &
#ifdef bidonref
          -rhpref_t(i,j,k)               &
#endif
                          )              &
              *dz_t(i,j,k,1)

!      if(i==21.and.j==jmax/2)write(6,*)k,kmin_w(i,j),real(anyv3d(i,j,k  ,id_prs)),real(depth_w(i,j,k))

      enddo ; enddo ; enddo

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const2=1./rho
      do k=kmax,1,-1 ; do j=2,jmax-1 ; do i=2,imax

!dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=  &
          mask_u(i,j,k)*const2*( & !PMXPMX>

          +0.5*( anyv3d(i+ipu,j,k+1,id_prs)-anyv3d(i-imu,j,k+1,id_prs)  &
                +anyv3d(i+ipu,j,k  ,id_prs)-anyv3d(i-imu,j,k  ,id_prs)) &

      +(depth_t(i+ipu,j,k)-depth_t(i-imu,j,k))*( & !ooo>

       (  anyv3d(i+ipu,j,k,id_prsref_)-anyv3d(i+ipu,j,k+1,id_prsref_)  &
         +anyv3d(i-imu,j,k,id_prsref_)-anyv3d(i-imu,j,k+1,id_prsref_)  &
       )/(dz_t(i+ipu,j,k,1)+dz_t(i-imu,j,k,1))                         &

#ifdef bidonref
         +grav*rhpref_u(i,j,k)                 &
#endif
                                             ) & !ooo>

! rhzavr (if any) contribution:
           -grav*depth_u(i,j,k)*( rhpzavr_w(i+ipu,j)                       &
                                 -rhpzavr_w(i-imu,j))  ) !PMXPMX>

!      if(j==jmax/2.and.i==24) then
!          write(67,*)k,presgrad_u(i,j,k,1)
!      endif

      enddo ; enddo ; enddo

!     j=jmax/2
!     do i=1,imax
!       k=kmin_u(i,j)   
!       write(50,*)i,presgrad_u(i,j,k,1)
!     enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const2=1./rho
      do k=kmax,1,-1 ; do j=2,jmax ; do i=2,imax-1

! dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=  &
          mask_v(i,j,k)*const2*( & !PMXPMX>

          +0.5*( anyv3d(i,j+jpv,k+1,id_prs)-anyv3d(i,j-jmv,k+1,id_prs)  &
                +anyv3d(i,j+jpv,k  ,id_prs)-anyv3d(i,j-jmv,k  ,id_prs)) &

      +(depth_t(i,j+jpv,k)-depth_t(i,j-jmv,k))*( & !ooo>
           
        (   anyv3d(i,j+jpv,k,id_prsref_)-anyv3d(i,j+jpv,k+1,id_prsref_) &
           +anyv3d(i,j-jmv,k,id_prsref_)-anyv3d(i,j-jmv,k+1,id_prsref_) &
           )/(dz_t(i,j+jpv,k,1)+dz_t(i,j-jmv,k,1))                      &

#ifdef bidonref
         +grav*rhpref_v(i,j,k)                   &
#endif
                                               ) & !ooo>

! rhzavr (if any) contribution:
           -grav* depth_v(i,j,k)*( rhpzavr_w(i,j+jpv)                       &
                                  -rhpzavr_w(i,j-jmv))  ) !PMXPMX>

      enddo ; enddo ; enddo

      else                 !3D3D3D3D3D> or !1D1D1D1D1D> !16-02-16

       call pressure_gradient_1dvcase !29-11-15

      endif                                !1D1D1D1D1D>


! Calculer la pression au niveau "z1"
!     do j=1,jmax ; do i=1,imax

! Profondeur "z(fond+1 continu"
!      depth_w(i,j,0)=depth_w(i,j,kmin_w(i,j)) &
!                  +( depth_w(i,j,kmin_w(i,j)+2)-depth_w(i,j,kmin_w(i,j)))*botlevmerged_w(i,j)  

!      rap=min(max( (depth_w(i,j,0)            -depth_w(i,j,kmin_w(i,j)+1) ) &
!                  /(depth_w(i,j,kmin_w(i,j)+2)-depth_w(i,j,kmin_w(i,j)+1) ),0.),1.)

!      anyv3d(i,j,0,id_prs)=(1.-rap)*anyv3d(i,j,kmin_w(i,j)+1,id_prs) &
!                              +rap *anyv3d(i,j,kmin_w(i,j)+2,id_prs)  

!     enddo        ; enddo

!#ifdef bidon
! COUCHES MERGED:
      k1=1 ; k2=0 ; k3=1 ; k4=0 ; rap=1.
      do j=2,jmax-1
      do i=2,imax

!     do k=1,kmerged_u(i,j)-1
!      presgrad_u(i,j,k,1)=presgrad_u(i,j,kmerged_u(i,j),1)
!     enddo

!#ifdef bidon

#ifdef bidon
!     rap=min(dz_u(i,j,kmin_u(i,j),1)/dz_u(i,j,kmin_u(i,j)+1,1),1.)
      rap=1.
! Cas standard: k et k+1 en i-imu (k1,k2) et en i+ipu (k3,k4)
      if(kmin_w(i-1,j)==kmin_w(i,j)) then     !pmx>
       k1=1 ; k2=0 ; k3=1 ; k4=0
      else                                    !pmx>
       if(kmin_w(i-1,j)<kmin_w(i,j)) then !oo>
        k1=kmin_w(i  ,j) ; k2=k1+1 ; k3=1 ; k4=0
       else                               !oo>
        k1=1 ; k2=0 ; k3=kmin_w(i-1,j) ; k4=k3+1
       endif                              !oo>
      endif                                   !pmx>
#endif

!     k2=kmerged_t(i-imu,j)+1
!     k4=kmerged_t(i+ipu,j)+1

      k=kmerged_u(i,j)-1
      k4=k+1
      k2=k+1
      if(dsig_t(i,j,k)<dsig_t(i-1,j,k)) then !>>>
       k1=0 ; k3=k  
       depth_w(i-1,j,k1)=depth_w(i-1,j,k+1)-dsig_t(i,j,k)*hz_w(i-1,j,1)
       anyv3d(i-imu,j,k1,id_prs)=( &
            (depth_w(i-1,j,k1 )-depth_w(i-1,j,k ))*anyv3d(i-imu,j,k+1,id_prs) &
           +(depth_w(i-1,j,k+1)-depth_w(i-1,j,k1))*anyv3d(i-imu,j,k  ,id_prs) &
          )/(depth_w(i-1,j,k+1)-depth_w(i-1,j,k  ))
    
!      if(i==24.and.j==jmax/2) then
!       write(6,*)'kmin_u(i,j)',kmin_u(i,j)
!       write(6,*)'kmerged_u(i,j)',kmerged_u(i,j)
!       write(6,*)'2dz i',depth_w(i,j,k+1)-depth_w(i,j,k),hz_w(i,j,1)*dsig_t(i,j,k)
!       write(6,*)'2dz i-1',depth_w(i-1,j,k+1)-depth_w(i-1,j,k),hz_w(i-1,j,1)*dsig_t(i-1,j,k)
!       write(6,*)'2z i  ',depth_w(i  ,j,k:k+1)
!       write(6,*)'2z i-1',depth_w(i-1,j,k:k+1)
!       write(6,*)'zk1 i-1',depth_w(i-1,j,k1)
!       write(6,*)'2 dsig',dsig_t(i-1:i,j,k)
!       write(6,*)'2P i-1',anyv3d(i-1,j,k:k+1,id_prs)
!       write(6,*)'Pk1 i-1',anyv3d(i-imu,j,k1,id_prs)
!      endif

      else                                   !>>>
       k3=0 ; k1=k  
       depth_w(i,j,k3)=depth_w(i,j,k+1)-dsig_t(i-1,j,k)*hz_w(i,j,1)
       anyv3d(i+ipu,j,k3,id_prs)=( &
            (depth_w(i,j,k3 )-depth_w(i,j,k ))*anyv3d(i+ipu,j,k+1,id_prs) &
           +(depth_w(i,j,k+1)-depth_w(i,j,k3))*anyv3d(i+ipu,j,k  ,id_prs) &
          )/(depth_w(i,j,k+1)-depth_w(i,j,k))
!      if(i==37.and.j==jmax/2) then
!       write(6,*)'----------------------------'
!       write(6,*)'kmin_u(i,j)',kmin_u(i,j)
!       write(6,*)'kmerged_u(i,j)',kmerged_u(i,j)
!       write(6,*)'2dz i',depth_w(i,j,k+1)-depth_w(i,j,k),hz_w(i,j,1)*dsig_t(i,j,k)
!       write(6,*)'2dz i-1',depth_w(i-1,j,k+1)-depth_w(i-1,j,k),hz_w(i-1,j,1)*dsig_t(i-1,j,k)
!       write(6,*)'2z i  ',depth_w(i  ,j,k:k+1)
!       write(6,*)'2z i-1',depth_w(i-1,j,k:k+1)
!       write(6,*)'zk3 i',depth_w(i,j,k3)
!       write(6,*)'2P i',anyv3d(i,j,k:k+1,id_prs)
!       write(6,*)'Pk3 i',anyv3d(i+ipu,j,k3,id_prs)
!       stop 'PAS OK'
!      endif
      endif                                  !>>>


!dx_u(i,j) * dp/dx:
!     presgrad_u(i,j,1:kmin_u(i,j),1)=const2*( & !PMXPMX>
      presgrad_u(i,j,1:kmerged_u(i,j)-1,1)=const2*( & !PMXPMX>

      +0.5*( anyv3d(i+ipu,j,k4,id_prs)-anyv3d(i-imu,j,k2,id_prs)   & ! delta P(bottom+1)
            +anyv3d(i+ipu,j,k3,id_prs)-anyv3d(i-imu,j,k1,id_prs) ) & ! delta P(bottom)

         +0.5*(depth_w(i+ipu,j,k4)-depth_w(i-imu,j,k2)    &          ! Z(bottom+1)
              +depth_w(i+ipu,j,k3)-depth_w(i-imu,j,k1))*( & !ooo>    ! Z(bottom)

        (  anyv3d(i+ipu,j,k3,id_prs)-anyv3d(i+ipu,j,k4,id_prs) &
          +anyv3d(i-imu,j,k1,id_prs)-anyv3d(i-imu,j,k2,id_prs) &
          )/(depth_w(i+ipu,j,k4)-depth_w(i+ipu,j,k3)           &
            +depth_w(i-imu,j,k2)-depth_w(i-imu,j,k1))          &

                                                      ) & !ooo>


           -grav*0.25*(depth_w(i+ipu,j,k4)     & ! < depth_u(i,j,k) >
                      +depth_w(i-imu,j,k2)     &
                      +depth_w(i+ipu,j,k3)     &
                      +depth_w(i-imu,j,k1))    &
                               *( rhpzavr_w(i+ipu,j)     &
                                 -rhpzavr_w(i-imu,j))  )   !PMXPMX>
                                                       !*rap+(1.-rap)*presgrad_u(i,j,kmin_u(i,j)+1,1)


      enddo
      enddo

      k1=1 ; k2=0 ; k3=1 ; k4=0 ; rap=1.
      do j=2,jmax
      do i=2,imax-1

      k=kmerged_v(i,j)-1
      k4=k+1
      k2=k+1
      if(dsig_t(i,j,k)<dsig_t(i,j-1,k)) then !>>>
       k1=0 ; k3=k  
       depth_w(i,j-1,k1)=depth_w(i,j-1,k+1)-dsig_t(i,j,k)*hz_w(i,j-1,1)
       anyv3d(i,j-jmv,k1,id_prs)=( &
            (depth_w(i,j-1,k1 )-depth_w(i,j-1,k ))*anyv3d(i,j-jmv,k+1,id_prs) &
           +(depth_w(i,j-1,k+1)-depth_w(i,j-1,k1))*anyv3d(i,j-jmv,k  ,id_prs) &
          )/(depth_w(i,j-1,k+1)-depth_w(i,j-1,k  ))
    
      else                                   !>>>
       k3=0 ; k1=k  
       depth_w(i,j,k3)=depth_w(i,j,k+1)-dsig_t(i,j-1,k)*hz_w(i,j,1)
       anyv3d(i,j+jpv,k3,id_prs)=( &
            (depth_w(i,j,k3 )-depth_w(i,j,k ))*anyv3d(i,j+jpv,k+1,id_prs) &
           +(depth_w(i,j,k+1)-depth_w(i,j,k3))*anyv3d(i,j+jpv,k  ,id_prs) &
          )/(depth_w(i,j,k+1)-depth_w(i,j,k))
      endif                                  !>>>


      presgrad_v(i,j,1:kmerged_v(i,j)-1,1)=const2*( & !PMXPMX>

      +0.5*( anyv3d(i,j+jpv,k4,id_prs)-anyv3d(i,j-jmv,k2,id_prs)   & ! delta P(bottom+1)
            +anyv3d(i,j+jpv,k3,id_prs)-anyv3d(i,j-jmv,k1,id_prs) ) & ! delta P(bottom)

         +0.5*(depth_w(i,j+jpv,k4)-depth_w(i,j-jmv,k2)    &          ! Z(bottom+1)
              +depth_w(i,j+jpv,k3)-depth_w(i,j-jmv,k1))*( & !ooo>    ! Z(bottom)

        (  anyv3d(i,j+jpv,k3,id_prs)-anyv3d(i,j+jpv,k4,id_prs) &
          +anyv3d(i,j-jmv,k1,id_prs)-anyv3d(i,j-jmv,k2,id_prs) &
          )/(depth_w(i,j+jpv,k4)-depth_w(i,j+jpv,k3)           &
            +depth_w(i,j-jmv,k2)-depth_w(i,j-jmv,k1))          &

                                                      ) & !ooo>


           -grav*0.25*(depth_w(i,j+jpv,k4)     & ! < depth_v(i,j,k) >
                      +depth_w(i,j-jmv,k2)     &
                      +depth_w(i,j+jpv,k3)     &
                      +depth_w(i,j-jmv,k1))    &
                               *( rhpzavr_w(i,j+jpv)     &
                                 -rhpzavr_w(i,j-jmv))  )   !PMXPMX>
                                                       !*rap+(1.-rap)*presgrad_v(i,j,kmin_v(i,j)+1,1)


      enddo
      enddo

!....................................
! The TKE scheme is coming next
! rhp_t possibly needs to be redefined if the level of reference of the potential
! density in the TKE scheme is different from that used in the PGF:

!       if(eos_author/=0) then
! Systematiquement calculer la densite potentielle maintenant que celle-ci est egalement !24-04-16
! utilisee dans le calcul de omega_w(:,:,kmax+1) et le calcul de rho.cp dans vertmix.F90
       if(eos_author==0) then !>>>>> !25-04-16
           call equation_of_state_linear(now)
       else                   !>>>>> !25-04-16
        if(eos_tkezref>=0.) then !-------->
          call equation_of_state_potzref2                 ! z reference = eos_tkezref !20-06-14
        else                     !-------->
          call equation_of_state_potloc2                  ! reference is current depth
        endif                    !-------->
       endif                  !>>>>> !25-04-16

! Marquer "disponibles" les tableaux generiques 
      anyv3d(-1,-1,0,1:4)=0. ; anyv3d(-1,-1,0,6:7)=0.
      anyv3d(-1,-1,0,id_prs)=0.

      end subroutine pressure_gradient_eos1

!.....................................................................................

!.....................................................................................

      subroutine pressure_gradient_pom                                 !10-08_10
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_pom'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!$ This routine computes the pressure gradient force related to the potential density
!$ according to the method described in:
!$ Marsaleix P., Auclair F., Estournel C., 2009. Low-order pressure gradient
!$ schemes in sigma coordinate models: The seamount test revisited.
!$ Ocean Modelling, 30, 169-177. http://dx.doi.org/10.1016/j.ocemod.2009.06.011
!$ The pressure gradient force related to the compression terms of the state equation
!$ is computed with a trapezoidal method

! Compute anyv3d(i,j,k,0) the part of the density related to pressure:
      call equation_of_state('compression terms',now)
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
        anyv3d(i,j,k,0)    &    ! full density
       = rhp_t(i,j,k)      &    !=potential density
       +anyv3d(i,j,k,0)         !+compression terms
      enddo
      enddo
      enddo

!  Computes the Hydrostatic pressure:
!  - using the rectangular method as far as the potential part of the density is concerned
!  - using the trapezoidal method as far as the "pressured" part of the density is concerned

! Note: rhp_t, the potential density - rho, has been computed at the previous step in scalars.F

! Reset of the summation at sea surface:
! At sea surface z=depth_w(kmax+1) and compression vanishes (i.e. anyv3d(i,j,k,0)=0)
      k=kmax
      const3=0.5*grav
      do j=1,jmax
      do i=1,imax

       anyv3d(i,j,k,2)=                                        &
                  +const3*                   anyv3d(i,j,k,0)   & ! compression part (trapezoidal)
                        *(depth_w(i,j,k+1) -depth_t(i,j,k))

      enddo
      enddo

! Summation over the kmax-1 remaining levels:
      const3=0.5*grav
      do 10 k=kmax-1,1,-1
       do 10 j=1,jmax
       do 10 i=1,imax

       anyv3d(i,j,k  ,2)=                                      &
       anyv3d(i,j,k+1,2)                                       &
                  +const3*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0))  & ! compression part (trapezoidal)
                        *(depth_t(i,j,k+1) -depth_t(i,j,k))

   10 continue

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

      do 14 k=kmax,1,-1
      do 14 j=2,jmax-1
      do 14 i=2,imax

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                              &
                   +const1*( anyv3d(i+ip,j,k,2)-anyv3d(i-im,j,k,2))   &
                   +const2*(depth_t(i+ip,j,k) -depth_t(i-im,j,k))*    &
                        0.5*(anyv3d(i+ip,j,k,0)+anyv3d(i-im,j,k,0))     ! The half half average for "pressured" density
                                                                        ! is consistent with the "trapezoidal" pressure
   14 continue

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09
      const2= grav/(rho)*flag3d                                            !30-09-09

      do 19 k=kmax,1,-1
      do 19 j=2,jmax
      do 19 i=2,imax-1

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                              &
                   +const1*( anyv3d(i,j+jp,k,2)-anyv3d(i,j-jm,k,2))   &
                   +const2*(depth_t(i,j+jp,k) -depth_t(i,j-jm,k))*    &
                        0.5*(anyv3d(i,j+jp,k,0)+anyv3d(i,j-jm,k,0))

   19 continue

      end subroutine pressure_gradient_pom

!------------------------------------------------------------------------

      subroutine pressure_gradient_jmfwg06
      use module_principal
      implicit none
      double precision a2_   ,a3_   ,a4_   ,rhc_   ,tem_   ,sal_     &
       ,tem2_   ,tem3_   ,prs_   ,a0_   ,a1_   ,a5_   ,a6_   ,a7_    &
       ,tem4_   ,sal2_   ,prs2_   ,prs3_
      integer time_
#ifdef synopsis
       subroutinetitle='pressure_gradient_jmfwg06'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      time_   =now

      a5_   =c23_jmfwg

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_   =tem_t(i,j,k,time_   )
      sal_   =max(sal_t(i,j,k,time_   ),zero)

! Compute density rho_    from (sal_   ,tem_   ,prs_   ) according
! to Jackett et al 2006 (JAOT):
      tem2_   = tem_   *tem_
      tem3_   =tem2_   *tem_
      tem4_   =tem3_   *tem_
      sal2_   =sal_   *sal_

      anyv3d(i,j,k,1)= c1_jmfwg                     &
                      +c2_jmfwg*tem_                &
                      +c3_jmfwg*tem2_               &
                      +c4_jmfwg*tem3_               &
                      +c5_jmfwg*sal_                &
                      +c6_jmfwg*sal_   *tem_        &
                      +c7_jmfwg*sal2_

      anyv3d(i,j,k,2)= c8_jmfwg                     &
                      +c9_jmfwg*tem2_               &
                     +c10_jmfwg*sal_

      anyv3d(i,j,k,3)=c11_jmfwg                    &
                     +c12_jmfwg*tem2_

      anyv3d(i,j,k,4)=1.                           &
                     +c14_jmfwg*tem_               &
                     +c15_jmfwg*tem2_              &
                     +c16_jmfwg*tem3_              &
                     +c17_jmfwg*tem4_              &
                     +c18_jmfwg*sal_               &
                     +c19_jmfwg*tem_   *sal_       &
                     +c20_jmfwg*sal_   *tem3_      &
                     +(sal_   **1.5)*(c21_jmfwg+c22_jmfwg*tem2_   )

!     anyv3d(i,j,k,5)=c23_jmfwg

      anyv3d(i,j,k,6)=c24_jmfwg*tem3_

      anyv3d(i,j,k,7)=c25_jmfwg*tem_

! Potential density anomaly:
      rhp_t(i,j,k)=anyv3d(i,j,k,1)/anyv3d(i,j,k,4) - rho

      enddo
      enddo
      enddo

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=0.
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i-im,j,k)*dz_t(i+ip,j,k,1)  & !09-01-11
                      +depth_t(i+ip,j,k)*dz_t(i-im,j,k,1)) &
                     /(                  dz_t(i+ip,j,k,1)  &
                                        +dz_t(i-im,j,k,1))
! Hydrostatic pressure at i=i-im:

! Compression terms:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =-0.25*depth_w(i-im,j,k+1)    &
              -0.25*depth_w(i+ip,j,k+1)    &
               -0.5*anyv3d(i,j,k,0)
      prs2_   = prs_   *prs_
      prs3_   =prs2_   *prs_

      a1_   =alpha*anyv3d(i-im,j,k+1,1)+(1.-alpha)*anyv3d(i-im,j,k,1)
      a2_   =alpha*anyv3d(i-im,j,k+1,2)+(1.-alpha)*anyv3d(i-im,j,k,2)
      a3_   =alpha*anyv3d(i-im,j,k+1,3)+(1.-alpha)*anyv3d(i-im,j,k,3)
      a4_   =alpha*anyv3d(i-im,j,k+1,4)+(1.-alpha)*anyv3d(i-im,j,k,4)
!     a5_   =alpha*anyv3d(i-im,j,k+1,5)+(1.-alpha)*anyv3d(i-im,j,k,5)
      a6_   =alpha*anyv3d(i-im,j,k+1,6)+(1.-alpha)*anyv3d(i-im,j,k,6)
      a7_   =alpha*anyv3d(i-im,j,k+1,7)+(1.-alpha)*anyv3d(i-im,j,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_


      xy_t(i-im,j,1)=                                                &
                     grav*(depth_w(i-im,j,k+1)                       &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i-im,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i-im,j,k)        &
                                     +rhc_    )

! Hydrostatic pressure at i=i+ip:
      a1_   =alpha*anyv3d(i+ip,j,k+1,1)+(1.-alpha)*anyv3d(i+ip,j,k,1)
      a2_   =alpha*anyv3d(i+ip,j,k+1,2)+(1.-alpha)*anyv3d(i+ip,j,k,2)
      a3_   =alpha*anyv3d(i+ip,j,k+1,3)+(1.-alpha)*anyv3d(i+ip,j,k,3)
      a4_   =alpha*anyv3d(i+ip,j,k+1,4)+(1.-alpha)*anyv3d(i+ip,j,k,4)
      a6_   =alpha*anyv3d(i+ip,j,k+1,6)+(1.-alpha)*anyv3d(i+ip,j,k,6)
      a7_   =alpha*anyv3d(i+ip,j,k+1,7)+(1.-alpha)*anyv3d(i+ip,j,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_

      xy_t(i+ip,j,2)=                                                &
                     grav*(depth_w(i+ip,j,k+1)                       &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i+ip,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i+ip,j,k)        &
                                     +rhc_    )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &
                +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1

      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i-im,j,k)*dz_t(i+ip,j,k,1)  & !09-01-11
                      +depth_t(i+ip,j,k)*dz_t(i-im,j,k,1)) &
                     /(                  dz_t(i+ip,j,k,1)  &
                                        +dz_t(i-im,j,k,1))

! Hydrostatic pressure at i=i-im:

! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i-im,j,k+1))    &  !09-08-10
           /(anyv3d(i,j,k+1,0)- anyv3d(i   ,j,k,0))

! Compression terms:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =-0.5*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0))
      prs2_   = prs_   *prs_
      prs3_   =prs2_   *prs_

      a1_   =alpha*anyv3d(i-im,j,k+1,1)+(1.-alpha)*anyv3d(i-im,j,k,1)
      a2_   =alpha*anyv3d(i-im,j,k+1,2)+(1.-alpha)*anyv3d(i-im,j,k,2)
      a3_   =alpha*anyv3d(i-im,j,k+1,3)+(1.-alpha)*anyv3d(i-im,j,k,3)
      a4_   =alpha*anyv3d(i-im,j,k+1,4)+(1.-alpha)*anyv3d(i-im,j,k,4)
      a6_   =alpha*anyv3d(i-im,j,k+1,6)+(1.-alpha)*anyv3d(i-im,j,k,6)
      a7_   =alpha*anyv3d(i-im,j,k+1,7)+(1.-alpha)*anyv3d(i-im,j,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_

      xy_t(i-im,j,1)=                                                &
      xy_t(i-im,j,1)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i-im,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i-im,j,k)        &
                                     +rhc_    )

! Hydrostatic pressure at i=i+ip:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i+ip,j,k+1))          &
           /(anyv3d(i,j,k+1,0)- anyv3d(i   ,j,k,0))

      a1_   =alpha*anyv3d(i+ip,j,k+1,1)+(1.-alpha)*anyv3d(i+ip,j,k,1)
      a2_   =alpha*anyv3d(i+ip,j,k+1,2)+(1.-alpha)*anyv3d(i+ip,j,k,2)
      a3_   =alpha*anyv3d(i+ip,j,k+1,3)+(1.-alpha)*anyv3d(i+ip,j,k,3)
      a4_   =alpha*anyv3d(i+ip,j,k+1,4)+(1.-alpha)*anyv3d(i+ip,j,k,4)
      a6_   =alpha*anyv3d(i+ip,j,k+1,6)+(1.-alpha)*anyv3d(i+ip,j,k,6)
      a7_   =alpha*anyv3d(i+ip,j,k+1,7)+(1.-alpha)*anyv3d(i+ip,j,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_

      xy_t(i+ip,j,2)=                                                &
      xy_t(i+ip,j,2)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i+ip,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i+ip,j,k)        &
                                     +rhc_    )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &
                +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1

      enddo
      enddo
      enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging:
      alpha=0.
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i,j-jm,k)*dz_t(i,j+jp,k,1)   & !09-01-11
                      +depth_t(i,j+jp,k)*dz_t(i,j-jm,k,1))  &
                     /(                  dz_t(i,j+jp,k,1)   &
                                        +dz_t(i,j-jm,k,1))

! Hydrostatic pressure at j=j-jm:

! Compression terms:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =-0.25*depth_w(i,j-jm,k+1)    &
              -0.25*depth_w(i,j+jp,k+1)    &
               -0.5*anyv3d(i,j,k,0)
      prs2_   = prs_   *prs_
      prs3_   =prs2_   *prs_

      a1_   =alpha*anyv3d(i,j-jm,k+1,1)+(1.-alpha)*anyv3d(i,j-jm,k,1)
      a2_   =alpha*anyv3d(i,j-jm,k+1,2)+(1.-alpha)*anyv3d(i,j-jm,k,2)
      a3_   =alpha*anyv3d(i,j-jm,k+1,3)+(1.-alpha)*anyv3d(i,j-jm,k,3)
      a4_   =alpha*anyv3d(i,j-jm,k+1,4)+(1.-alpha)*anyv3d(i,j-jm,k,4)
      a6_   =alpha*anyv3d(i,j-jm,k+1,6)+(1.-alpha)*anyv3d(i,j-jm,k,6)
      a7_   =alpha*anyv3d(i,j-jm,k+1,7)+(1.-alpha)*anyv3d(i,j-jm,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_

      xy_t(i,j-jm,1)=                                                &
                     grav*(depth_w(i,j-jm,k+1)                       &
                           -anyv3d(i,j   ,k  ,0))                    &
                                *(     alpha *rhp_t(i,j-jm,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j-jm,k)        &
                                     +rhc_    )

! Hydrostatic pressure at j=j+jp:
      a1_   =alpha*anyv3d(i,j+jp,k+1,1)+(1.-alpha)*anyv3d(i,j+jp,k,1)
      a2_   =alpha*anyv3d(i,j+jp,k+1,2)+(1.-alpha)*anyv3d(i,j+jp,k,2)
      a3_   =alpha*anyv3d(i,j+jp,k+1,3)+(1.-alpha)*anyv3d(i,j+jp,k,3)
      a4_   =alpha*anyv3d(i,j+jp,k+1,4)+(1.-alpha)*anyv3d(i,j+jp,k,4)
      a6_   =alpha*anyv3d(i,j+jp,k+1,6)+(1.-alpha)*anyv3d(i,j+jp,k,6)
      a7_   =alpha*anyv3d(i,j+jp,k+1,7)+(1.-alpha)*anyv3d(i,j+jp,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_

      xy_t(i,j+jp,2)=                                                &
                     grav*(depth_w(i,j+jp,k+1)                       &
                           -anyv3d(i,j   ,k  ,0))                    &
                                *(     alpha *rhp_t(i,j+jp,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j+jp,k)        &
                                     +rhc_    )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                      &
                +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1


      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i,j-jm,k)*dz_t(i,j+jp,k,1)   & !09-01-11
                      +depth_t(i,j+jp,k)*dz_t(i,j-jm,k,1))  &
                     /(                  dz_t(i,j+jp,k,1)   &
                                        +dz_t(i,j-jm,k,1))

! Hydrostatic pressure at j=j-jm:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i,j-jm,k+1))          &
           /(anyv3d(i,j,k+1,0)- anyv3d(i,j,k,0))

! Compression terms:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =-0.5*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0))
      prs2_   = prs_   *prs_
      prs3_   =prs2_   *prs_

      a1_   =alpha*anyv3d(i,j-jm,k+1,1)+(1.-alpha)*anyv3d(i,j-jm,k,1)
      a2_   =alpha*anyv3d(i,j-jm,k+1,2)+(1.-alpha)*anyv3d(i,j-jm,k,2)
      a3_   =alpha*anyv3d(i,j-jm,k+1,3)+(1.-alpha)*anyv3d(i,j-jm,k,3)
      a4_   =alpha*anyv3d(i,j-jm,k+1,4)+(1.-alpha)*anyv3d(i,j-jm,k,4)
      a6_   =alpha*anyv3d(i,j-jm,k+1,6)+(1.-alpha)*anyv3d(i,j-jm,k,6)
      a7_   =alpha*anyv3d(i,j-jm,k+1,7)+(1.-alpha)*anyv3d(i,j-jm,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_


      xy_t(i,j-jm,1)=                                                &
      xy_t(i,j-jm,1)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i,j-jm,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j-jm,k)        &
                                     +rhc_    )

! Hydrostatic pressure at j=j+jp:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i,j+jp,k+1))          &
           /(anyv3d(i,j,k+1,0)- anyv3d(i,j,k,0))

      a1_   =alpha*anyv3d(i,j+jp,k+1,1)+(1.-alpha)*anyv3d(i,j+jp,k,1)
      a2_   =alpha*anyv3d(i,j+jp,k+1,2)+(1.-alpha)*anyv3d(i,j+jp,k,2)
      a3_   =alpha*anyv3d(i,j+jp,k+1,3)+(1.-alpha)*anyv3d(i,j+jp,k,3)
      a4_   =alpha*anyv3d(i,j+jp,k+1,4)+(1.-alpha)*anyv3d(i,j+jp,k,4)
      a6_   =alpha*anyv3d(i,j+jp,k+1,6)+(1.-alpha)*anyv3d(i,j+jp,k,6)
      a7_   =alpha*anyv3d(i,j+jp,k+1,7)+(1.-alpha)*anyv3d(i,j+jp,k,7)
      rhc_   =(a1_   +a2_   *prs_   +a3_   *prs2_   )                  &
             /(a4_   +a5_   *prs_   +a6_   *prs2_   +a7_   *prs3_   )  &
              -a1_   /a4_

      xy_t(i,j+jp,2)=                                                &
      xy_t(i,j+jp,2)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i,j+jp,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j+jp,k)        &
                                     +rhc_    )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                      &
                +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1


      enddo
      enddo
      enddo

      end subroutine pressure_gradient_jmfwg06

!..............................................................................

      subroutine pressure_gradient_w97
      use module_principal
      implicit none
      double precision a2_   ,a3_   ,a4_   ,rhc_   ,tem_   ,sal_    &
       ,tem2_   ,tem3_   ,prs_   ,a0_   ,l_   ,p0_
      integer time_
#ifdef synopsis
       subroutinetitle='pressure_gradient_w97'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      time_   =now

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      tem_   =tem_t(i,j,k,time_   )
      tem2_   = tem_   *tem_
      tem3_   =tem2_   *tem_

      sal_   =sal_t(i,j,k,time_   )

! Pression en Pa:
!     prs_   =-rho*grav*depth_t(i,j,k)
!     prs_   =  -10000.*depth_t(i,j,k)

! Coef alpha0:
      a0_   =a0_w97+a1_w97*tem_   +a2_w97*sal_
! Coef p0:
      p0_   =b0_w97+b1_w97*tem_   +b2_w97*tem2_   +b3_w97*tem3_      &
                   +b4_w97*sal_   +b5_w97*tem_   *sal_
! Coef lambda:
      l_   = c0_w97+c1_w97*tem_   +c2_w97*tem2_   +c3_w97*tem3_      &
                   +c4_w97*sal_   +c5_w97*tem_   *sal_

      x0=l_   +a0_   *p0_
      anyv3d(i,j,k,2)=l_
      anyv3d(i,j,k,3)=x0**2
      anyv3d(i,j,k,4)=a0_   *x0

! Potential density anomaly:
      rhp_t(i,j,k)=p0_   /(l_   +a0_   *p0_   ) - rho

!     rho_   =anyv3d(i,j,k,2)*prs_        &
!           /(anyv3d(i,j,k,3)+prs_   *anyv3d(i,j,k,4))
      enddo
      enddo
      enddo

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=0.
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i-im,j,k)*dz_t(i+ip,j,k,1)  &
                      +depth_t(i+ip,j,k)*dz_t(i-im,j,k,1)) &
                     /(                  dz_t(i+ip,j,k,1)  &
                                        +dz_t(i-im,j,k,1))

! Hydrostatic pressure at i=i-im:

! Compression terms of the wright97 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.25*depth_w(i-im,j,k+1)    &
        +0.25*depth_w(i+ip,j,k+1)    &
         +0.5*anyv3d(i,j,k,0)

      a2_   =alpha*anyv3d(i-im,j,k+1,2)+(1.-alpha)*anyv3d(i-im,j,k,2)
      a3_   =alpha*anyv3d(i-im,j,k+1,3)+(1.-alpha)*anyv3d(i-im,j,k,3)
      a4_   =alpha*anyv3d(i-im,j,k+1,4)+(1.-alpha)*anyv3d(i-im,j,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i-im,j,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i-im,j,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i-im,j,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i-im,j,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i-im,j,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i-im,j,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2



      xy_t(i-im,j,1)=                                                &
                     grav*(depth_w(i-im,j,k+1)                       &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i-im,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i-im,j,k)        &
                                     +rhc_    )

! Hydrostatic pressure at i=i+ip:
      a2_   =alpha*anyv3d(i+ip,j,k+1,2)+(1.-alpha)*anyv3d(i+ip,j,k,2)
      a3_   =alpha*anyv3d(i+ip,j,k+1,3)+(1.-alpha)*anyv3d(i+ip,j,k,3)
      a4_   =alpha*anyv3d(i+ip,j,k+1,4)+(1.-alpha)*anyv3d(i+ip,j,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i+ip,j,2)=                                                &
                     grav*(depth_w(i+ip,j,k+1)                       &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i+ip,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i+ip,j,k)        &
                                     +rhc_    )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &
                +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1

      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax-1
      do i=2,imax

! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i-im,j,k)*dz_t(i+ip,j,k,1)  &
                      +depth_t(i+ip,j,k)*dz_t(i-im,j,k,1)) &
                     /(                  dz_t(i+ip,j,k,1)  &
                                        +dz_t(i-im,j,k,1))

! Hydrostatic pressure at i=i-im:

! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i-im,j,k+1))    &  !09-08-10
           /(anyv3d(i,j,k+1,0)- anyv3d(i   ,j,k,0))

! Compression terms of the w97 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.5*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0))
!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i-im,j,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i-im,j,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i-im,j,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i-im,j,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i-im,j,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i-im,j,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2
      a2_   =alpha*anyv3d(i-im,j,k+1,2)+(1.-alpha)*anyv3d(i-im,j,k,2)
      a3_   =alpha*anyv3d(i-im,j,k+1,3)+(1.-alpha)*anyv3d(i-im,j,k,3)
      a4_   =alpha*anyv3d(i-im,j,k+1,4)+(1.-alpha)*anyv3d(i-im,j,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i-im,j,1)=                                                &
      xy_t(i-im,j,1)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i-im,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i-im,j,k)        &
                                     +rhc_    )

! Hydrostatic pressure at i=i+ip:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i+ip,j,k+1))          &
           /(anyv3d(i,j,k+1,0)- anyv3d(i   ,j,k,0))

!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i+ip,j,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i+ip,j,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i+ip,j,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i+ip,j,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i+ip,j,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i+ip,j,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2
      a2_   =alpha*anyv3d(i+ip,j,k+1,2)+(1.-alpha)*anyv3d(i+ip,j,k,2)
      a3_   =alpha*anyv3d(i+ip,j,k+1,3)+(1.-alpha)*anyv3d(i+ip,j,k,3)
      a4_   =alpha*anyv3d(i+ip,j,k+1,4)+(1.-alpha)*anyv3d(i+ip,j,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i+ip,j,2)=                                                &
      xy_t(i+ip,j,2)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i+ip,j,k+1)      &
                                  +(1.-alpha)*rhp_t(i+ip,j,k)        &
                                     +rhc_    )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &
                +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1

      enddo
      enddo
      enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging:
      alpha=0.
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i,j-jm,k)*dz_t(i,j+jp,k,1)   &
                      +depth_t(i,j+jp,k)*dz_t(i,j-jm,k,1))  &
                     /(                  dz_t(i,j+jp,k,1)   &
                                        +dz_t(i,j-jm,k,1))

! Hydrostatic pressure at j=j-jm:

! Compression terms of the w97 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.25*depth_w(i,j-jm,k+1)    &
        +0.25*depth_w(i,j+jp,k+1)    &
          +0.5*anyv3d(i,j,k,0)

!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i,j-jm,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i,j-jm,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i,j-jm,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j-jm,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i,j-jm,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j-jm,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2
      a2_   =alpha*anyv3d(i,j-jm,k+1,2)+(1.-alpha)*anyv3d(i,j-jm,k,2)
      a3_   =alpha*anyv3d(i,j-jm,k+1,3)+(1.-alpha)*anyv3d(i,j-jm,k,3)
      a4_   =alpha*anyv3d(i,j-jm,k+1,4)+(1.-alpha)*anyv3d(i,j-jm,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j-jm,1)=                                                &
                     grav*(depth_w(i,j-jm,k+1)                       &
                           -anyv3d(i,j   ,k  ,0))                    &
                                *(     alpha *rhp_t(i,j-jm,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j-jm,k)        &
                                     +rhc_    )

! Hydrostatic pressure at j=j+jp:
!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i,j+jp,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i,j+jp,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i,j+jp,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j+jp,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i,j+jp,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j+jp,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2
      a2_   =alpha*anyv3d(i,j+jp,k+1,2)+(1.-alpha)*anyv3d(i,j+jp,k,2)
      a3_   =alpha*anyv3d(i,j+jp,k+1,3)+(1.-alpha)*anyv3d(i,j+jp,k,3)
      a4_   =alpha*anyv3d(i,j+jp,k+1,4)+(1.-alpha)*anyv3d(i,j+jp,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j+jp,2)=                                                &
                     grav*(depth_w(i,j+jp,k+1)                       &
                           -anyv3d(i,j   ,k  ,0))                    &
                                *(     alpha *rhp_t(i,j+jp,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j+jp,k)        &
                                     +rhc_    )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                      &
                +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1


      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference:
      anyv3d(i,j,k,0)=(depth_t(i,j-jm,k)*dz_t(i,j+jp,k,1)   &
                      +depth_t(i,j+jp,k)*dz_t(i,j-jm,k,1))  &
                     /(                  dz_t(i,j+jp,k,1)   &
                                        +dz_t(i,j-jm,k,1))

! Hydrostatic pressure at j=j-jm:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i,j-jm,k+1))          &
           /(anyv3d(i,j,k+1,0)- anyv3d(i,j,k,0))

! Compression terms of the w97 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.5*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0))
!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i,j-jm,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i,j-jm,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i,j-jm,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j-jm,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i,j-jm,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j-jm,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2
      a2_   =alpha*anyv3d(i,j-jm,k+1,2)+(1.-alpha)*anyv3d(i,j-jm,k,2)
      a3_   =alpha*anyv3d(i,j-jm,k+1,3)+(1.-alpha)*anyv3d(i,j-jm,k,3)
      a4_   =alpha*anyv3d(i,j-jm,k+1,4)+(1.-alpha)*anyv3d(i,j-jm,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j-jm,1)=                                                &
      xy_t(i,j-jm,1)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i,j-jm,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j-jm,k)        &
                                     +rhc_    )

! Hydrostatic pressure at j=j+jp:

! Coefficient for rho, T, S averaging:
      alpha=(anyv3d(i,j,k+1,0)-depth_w(i,j+jp,k+1))          &
           /(anyv3d(i,j,k+1,0)- anyv3d(i,j,k,0))

!     x0=-x1/(1402.3                                                &
!               +1.34*(     alpha *sal_t(i,j+jp,k+1,time_   )       &
!                      +(1.-alpha)*sal_t(i,j+jp,k  ,time_   ) )     &
!                    +(     alpha *tem_t(i,j+jp,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j+jp,k  ,time_   ) )     &
!        *(4.55-0.045*(     alpha *tem_t(i,j+jp,k+1,time_   )       &
!                      +(1.-alpha)*tem_t(i,j+jp,k  ,time_   )  ) )  &
!      +x1*(15.e-9*x1-0.00821) )**2
      a2_   =alpha*anyv3d(i,j+jp,k+1,2)+(1.-alpha)*anyv3d(i,j+jp,k,2)
      a3_   =alpha*anyv3d(i,j+jp,k+1,3)+(1.-alpha)*anyv3d(i,j+jp,k,3)
      a4_   =alpha*anyv3d(i,j+jp,k+1,4)+(1.-alpha)*anyv3d(i,j+jp,k,4)
      prs_   =-x1*10000.
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j+jp,2)=                                                &
      xy_t(i,j+jp,2)+grav*(anyv3d(i,j,k+1,0)                         &
                          -anyv3d(i,j,k  ,0))                        &
                                *(     alpha *rhp_t(i,j+jp,k+1)      &
                                  +(1.-alpha)*rhp_t(i,j+jp,k)        &
                                     +rhc_    )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                      &
                +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1


      enddo
      enddo
      enddo

      end subroutine pressure_gradient_w97

!------------------------------------------------------------------------

      subroutine pressure_gradient_jmfwg06_hyb
      use module_principal
      use module_parallele
      implicit none
      double precision w1_,w2_
#ifdef synopsis
       subroutinetitle='pressure_gradient_jmfwg06_hyb'             !20-10-14
       subroutinedescription= &
      'Computes the 3D pressure gradient according to Marsaleix' &
      //' et al, Ocean Modelling, 2011'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!---------------------------------------------------------
! Pour arbitrairement annuler la densite potentielle et
! calculer le PGF OM2011 avec uniquement la densite de compression:
!     write(6,*)'BIDOUILLE dans pressure_gradient_jmfwg06_hyb'
!     write(6,*)'BIDOUILLE dans pressure_gradient_jmfwg06_hyb'
!     rhp_t(:,:,:)=0.
!---------------------------------------------------------

! Compute rhp the potential density at z=0m and the a1, a2, ....,a7 coefficients of the
! Jackett et al 2006 (JAOT) EOS
      call equation_of_state_rhp_a1_a7(now)

#ifdef bidon
      write(6,*)'BIDOUILLE'
!     x2=1./real(iglb-1)/dxb
      x2=1./real(jglb-1)/dyb
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      rhp_t(i,j,k)=x2*real(i+par%timax(1)-1)*dxb
       rhp_t(i,j,k)=x2*real(j+par%tjmax(1)-1)*dyb
      enddo ; enddo ; enddo
#endif

      if(rhp_zavr_xy==1)call pressure_gradient_rhpzavr

! Compute anyv3d(i,j,k,0) the hydrostatic pressure associated to
! rhp the potential density (at z=0m)
      x0=0.5*grav
      do j=1,jmax
      do i=1,imax
       anyv3d(i,j,kmax,0)=x0*rhp_t(i,j,kmax)*dz_t(i,j,kmax,1)
      enddo
      enddo
      do k=kmax-1,1,-1
      do j=1,jmax
      do i=1,imax
       anyv3d(i,j,k  ,0)=                                     &
       anyv3d(i,j,k+1,0)+x0*( rhp_t(i,j,k  )*dz_t(i,j,k  ,1)  &
                             +rhp_t(i,j,k+1)*dz_t(i,j,k+1,1))
      enddo
      enddo
      enddo

!     i=3 ; j=3
!     do k=kmax,1,-1
!      write(6,*)depth_t(i,j,k),anyv3d(i,j,k  ,0),anyv3d(i,j,k  ,5)
!     enddo
!     stop'hanhan'


!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

      do j=2,jmax-1
      do i=1,imax
       xy_t(i-imu,j,1)=0.
       xy_t(i+ipu,j,2)=0.
      enddo
      enddo
      do k=kmax,1,-1
      k0=1 ; if(k==kmax)k0=0
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference: depth_u(i,j,k)

! Hydrostatic pressure at i=i-imu:

! Coefficient for rho, T, S averaging (w1_=0 if k=kmax):
      w1_=k0*(depth_u(i,j,k+1)-depth_w(i-imu,j,k+1))    &  !09-08-10
            /(depth_u(i,j,k+1)-depth_u(i   ,j,k))
      w2_=1.-w1_

! Compression terms:
! x1, the pressure (decibar), is the same for i=i-imu and i+ipu in order to
! remove "sigma type" troncation errors
!     prs_   =-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k))
!     prs2_   = (-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))   *(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))
!     prs3_   =prs2_   *(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))

      xy_t(i-imu,j,1)=                                      &
      xy_t(i-imu,j,1)+grav*(depth_u(i,j,k+1)                &
                           -depth_u(i,j,k  ))               &

      *(( (w1_*anyv3d(i-imu,j,k+1,1)+w2_*anyv3d(i-imu,j,k,1))          & !rhc_
        + (w1_*anyv3d(i-imu,j,k+1,2)+w2_*anyv3d(i-imu,j,k,2))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))     &
        + (w1_*anyv3d(i-imu,j,k+1,3)+w2_*anyv3d(i-imu,j,k,3))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))**2) &
       /( (w1_*anyv3d(i-imu,j,k+1,4)+w2_*anyv3d(i-imu,j,k,4))          &
                                                  + c23_jmfwg*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))     &
        + (w1_*anyv3d(i-imu,j,k+1,6)+w2_*anyv3d(i-imu,j,k,6))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))**2  &
        + (w1_*anyv3d(i-imu,j,k+1,7)+w2_*anyv3d(i-imu,j,k,7))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))**3) &

        - (w1_*anyv3d(i-imu,j,k+1,1)+w2_*anyv3d(i-imu,j,k,1))        &
        / (w1_*anyv3d(i-imu,j,k+1,4)+w2_*anyv3d(i-imu,j,k,4)))

! Hydrostatic pressure at i=i+ipu:

! Coefficient for rho, T, S averaging:
      w1_=k0*(depth_u(i,j,k+1)-depth_w(i+ipu,j,k+1))          &
            /(depth_u(i,j,k+1)-depth_u(i   ,j,k))
      w2_=1.-w1_

      xy_t(i+ipu,j,2)=                                    &
      xy_t(i+ipu,j,2)+grav*(depth_u(i,j,k+1)              &
                           -depth_u(i,j,k  ))             &

      *(( (w1_*anyv3d(i+ipu,j,k+1,1)+w2_*anyv3d(i+ipu,j,k,1))          & !rhc_
        + (w1_*anyv3d(i+ipu,j,k+1,2)+w2_*anyv3d(i+ipu,j,k,2))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))     &
        + (w1_*anyv3d(i+ipu,j,k+1,3)+w2_*anyv3d(i+ipu,j,k,3))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))**2) &
       /( (w1_*anyv3d(i+ipu,j,k+1,4)+w2_*anyv3d(i+ipu,j,k,4))          &
                                                   +c23_jmfwg*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))     &
        + (w1_*anyv3d(i+ipu,j,k+1,6)+w2_*anyv3d(i+ipu,j,k,6))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))**2  &
        + (w1_*anyv3d(i+ipu,j,k+1,7)+w2_*anyv3d(i+ipu,j,k,7))*(-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k)))**3) &

        - (w1_*anyv3d(i+ipu,j,k+1,1)+w2_*anyv3d(i+ipu,j,k,1))          &
        / (w1_*anyv3d(i+ipu,j,k+1,4)+w2_*anyv3d(i+ipu,j,k,4)))

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                    &

! pgf compression terms:
           +(xy_t(i+ipu,j,2)-xy_t(i-imu,j,1))*const1          &

! pgf potential density terms:
           +const1*( anyv3d(i+ipu,j,k,0)-anyv3d(i-imu,j,k,0)) &
           +const2*(depth_t(i+ipu,j,k) -depth_t(i-imu,j,k))   &
                  *(  rhp_t(i+ipu,j,k)  *  dz_t(i+ipu,j,k,1)  &
                     +rhp_t(i-imu,j,k)  *  dz_t(i-imu,j,k,1)) &
                   /(  dz_t(i+ipu,j,k,1) + dz_t(i-imu,j,k,1)) &

! rhzavr (if any) contribution:
           -const2*depth_u(i,j,k)*( rhpzavr_w(i+ipu,j)         &
                                   -rhpzavr_w(i-imu,j))

      enddo
      enddo
      enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

      do j=2,jmax
      do i=2,imax-1
       xy_t(i,j-jmv,1)=0.
       xy_t(i,j+jpv,2)=0.
      enddo
      enddo
      do k=kmax,1,-1
      k0=1  ; if(k==kmax)k0=0
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference: depth_v

! Hydrostatic pressure at j=j-jmv:

! Coefficient for rho, T, S averaging:
      w1_=k0*(depth_v(i,j,k+1)-depth_w(i,j-jmv,k+1))          &
            /(depth_v(i,j,k+1)-depth_v(i,j   ,k))
      w2_=1.-w1_

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
!     prs_   =-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k))
!     prs2_   = (-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))   *(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))
!     prs3_   =prs2_   *(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))

      xy_t(i,j-jmv,1)=                                                &
      xy_t(i,j-jmv,1)+grav*(depth_v(i,j,k+1)                         &
                           -depth_v(i,j,k  ))                        &

      *(( (w1_*anyv3d(i,j-jmv,k+1,1)+w2_*anyv3d(i,j-jmv,k,1))          & !rhc_
        + (w1_*anyv3d(i,j-jmv,k+1,2)+w2_*anyv3d(i,j-jmv,k,2))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))     &
        + (w1_*anyv3d(i,j-jmv,k+1,3)+w2_*anyv3d(i,j-jmv,k,3))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))**2) &
       /( (w1_*anyv3d(i,j-jmv,k+1,4)+w2_*anyv3d(i,j-jmv,k,4))          &
                                                  + c23_jmfwg*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))     &
        + (w1_*anyv3d(i,j-jmv,k+1,6)+w2_*anyv3d(i,j-jmv,k,6))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))**2  &
        + (w1_*anyv3d(i,j-jmv,k+1,7)+w2_*anyv3d(i,j-jmv,k,7))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))**3) &

        - (w1_*anyv3d(i,j-jmv,k+1,1)+w2_*anyv3d(i,j-jmv,k,1))        &
        / (w1_*anyv3d(i,j-jmv,k+1,4)+w2_*anyv3d(i,j-jmv,k,4)))

! Hydrostatic pressure at j=j+jpv:

! Coefficient for rho, T, S averaging:
      w1_=k0*(depth_v(i,j,k+1)-depth_w(i,j+jpv,k+1))          &
            /(depth_v(i,j,k+1)-depth_v(i,j   ,k))
      w2_=1.-w1_

      xy_t(i,j+jpv,2)=                                               &
      xy_t(i,j+jpv,2)+grav*(depth_v(i,j,k+1)                         &
                           -depth_v(i,j,k  ))                        &

      *(( (w1_*anyv3d(i,j+jpv,k+1,1)+w2_*anyv3d(i,j+jpv,k,1))          & !rhc_
        + (w1_*anyv3d(i,j+jpv,k+1,2)+w2_*anyv3d(i,j+jpv,k,2))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))     &
        + (w1_*anyv3d(i,j+jpv,k+1,3)+w2_*anyv3d(i,j+jpv,k,3))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))**2) &
       /( (w1_*anyv3d(i,j+jpv,k+1,4)+w2_*anyv3d(i,j+jpv,k,4))          &
                                                  + c23_jmfwg*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))     &
        + (w1_*anyv3d(i,j+jpv,k+1,6)+w2_*anyv3d(i,j+jpv,k,6))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))**2  &
        + (w1_*anyv3d(i,j+jpv,k+1,7)+w2_*anyv3d(i,j+jpv,k,7))*(-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k)))**3) &

        - (w1_*anyv3d(i,j+jpv,k+1,1)+w2_*anyv3d(i,j+jpv,k,1))        &
        / (w1_*anyv3d(i,j+jpv,k+1,4)+w2_*anyv3d(i,j+jpv,k,4))  )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                  &

! pgf compression terms:
          +(xy_t(i,j+jpv,2)-xy_t(i,j-jmv,1))*const1         &

! pgf potential density terms:
         +const1*( anyv3d(i,j+jpv,k,0)-anyv3d(i,j-jmv,k,0)) &
         +const2*(depth_t(i,j+jpv,k) -depth_t(i,j-jmv,k))   &
                *(  rhp_t(i,j+jpv,k)  *  dz_t(i,j+jpv,k,1)  &
                   +rhp_t(i,j-jmv,k)  *  dz_t(i,j-jmv,k,1)) &
                 /(  dz_t(i,j+jpv,k,1)+  dz_t(i,j-jmv,k,1)) & 

! rhzavr (if any) contribution:
         -const2*depth_v(i,j,k)*( rhpzavr_w(i,j+jpv)      &
                                 -rhpzavr_w(i,j-jmv))

      enddo
      enddo
      enddo

#ifdef bidon
      write(6,*)'BIDOUILLE pressure_gradient_2d_terms'
      call pressure_gradient_2d_terms
      j=jmax/2 ; k=kmax/2 ; i=imax/2
      x1=0.5*grav/rho
!     x2=1./real(iglb-1)/dxb     ! pente rhpzavr
!     x3=0.001/real(iglb-1)/dxb ! pente ssh
      x2=1./real(jglb-1)/dyb     ! pente rhpzavr
      x3=0.01/real(jglb-1)/dyb ! pente ssh
!     do i=2,imax
      do j=2,jmax
   
!       x0=real(i+par%timax(1)-1)*dxb
        x0=real(j+par%tjmax(1)-1)*dyb

!.......................
        write(66,*)                                            &
!       (presgrad_u(i,j,k,1)                                     &
!                      +x1*(                                   &
!       (ssh_w(i  ,j,0)+ssh_w(i  ,j,1))*(rho+rhpzavr_w(i  ,j)) &!grav/rho*d(ssh.rhp2d)/di
!      -(ssh_w(i-1,j,0)+ssh_w(i-1,j,1))*(rho+rhpzavr_w(i-1,j)) &!
!                          ) &
!        )/dxb         & 

        (presgrad_v(i,j,k,1)                                     &
                       +x1*(                                   &
        (ssh_w(i  ,j,0)+ssh_w(i  ,j,1))*(rho+rhpzavr_w(i  ,j)) &!grav/rho*d(ssh.rhp2d)/di
       -(ssh_w(i,j-1,0)+ssh_w(i,j-1,1))*(rho+rhpzavr_w(i,j-1)) &!
                           ) &
         )/dyb         & 

       ,grav*(                   &
           x3*rho                &
        +2*x2*x3*x0              &
          -depth_u(i,j,k)*x2     &
!         -depth_v(i,j,k)*x2     &
        )/rho


!.......................
       x5=(0.25*(ssh_w(i,j,0)+ssh_w(i-1,j,0)            &
                +ssh_w(i,j,1)+ssh_w(i-1,j,1))+h_u(i,j)) 
       write(67,*) &
                     x5*(x1*( & !ppp>

            (ssh_w(i  ,j,0)+ssh_w(i  ,j,1))*(rho+rhpzavr_w(i  ,j))  & !grav/rho*d(ssh.rhp2d)/di
           -(ssh_w(i-1,j,0)+ssh_w(i-1,j,1))*(rho+rhpzavr_w(i-1,j))  &

                           ) & !ppp>
           +pres3d2d_u(i,j)                                   &
         ) & !mmmm>
         /dx_u(i,j)        &
       ,(x5*grav*(x3*rho+2.*x2*x3*x0)   &
          +grav*0.5*x2*( h_u(i,j)**2   &
                        -(0.5*(ssh_w(i  ,j,0)+ssh_w(i-1,j,0)))**2 ))/rho

!.......................
        x5=(0.25*(ssh_w(i,j,0)+ssh_w(i,j-1,0)            &
                 +ssh_w(i,j,1)+ssh_w(i,j-1,1))+h_v(i,j))
       write(68,*) &
                     x5*(x1*( & !ppp>

            (ssh_w(i  ,j,0)+ssh_w(i  ,j,1))*(rho+rhpzavr_w(i  ,j))  & !grav/rho*d(ssh.rhp2d)/dj
           -(ssh_w(i,j-1,0)+ssh_w(i,j-1,1))*(rho+rhpzavr_w(i,j-1))  &

                           ) & !ppp>

           +pres3d2d_v(i,j)                                   &

         ) & !mmmm>
         /dy_v(i,j) &
       ,(x5*grav*(x3*rho+2.*x2*x3*x0)   &
           +grav*0.5*x2*( h_v(i,j)**2   &
                -(0.5*(ssh_w(i,j,0)+ssh_w(i,j-1,0)))**2 ))/rho

      enddo


#endif

!     cpu_seconds=MPI_Wtime ( ) - cpu_seconds
!     if(par%rank==0)write(20,*)cpu_seconds

!-----------------------------
! The TKE scheme is coming next
! The eos potential density in the PGF is referenced at z=0m
! If the level of reference of the potential density in the TKE scheme
! is different then recalculate rhp_t
      if(eos_tkezref/=0.) then !eos-tke-eos-tke->
            if(eos_tkezref<0.) then !>>>>>
               call equation_of_state_potloc2  ! reference is current depth
            else                    !>>>>>
               call equation_of_state_potzref2 ! z reference = eos_tkezref !20-06-14
            endif                   !>>>>w
      endif                    !eos-tke-eos-tke->
!-----------------------------

      end subroutine pressure_gradient_jmfwg06_hyb

!..............................................................................

      subroutine pressure_gradient_w97_hyb
      use module_principal
      implicit none
      double precision a2_   ,a3_   ,a4_   ,rhc_   ,tem_   ,sal_    &
       ,tem2_   ,tem3_   ,prs_   ,a0_   ,l_   ,p0_
      integer time_
#ifdef synopsis
       subroutinetitle='pressure_gradient_w97_hyb'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      time_   =now

      do j=1,jmax
      do i=1,imax
      xy_t(i,j,1)=0.
      enddo
      enddo

      do k=kmax,1,-1
      do j=1,jmax
      do i=1,imax

      tem_   =tem_t(i,j,k,time_   )
      tem2_   = tem_   *tem_
      tem3_   =tem2_   *tem_

      sal_   =sal_t(i,j,k,time_   )

! Pression en Pa:
!     prs_   =-rho*grav*depth_t(i,j,k)
!     prs_   =  -10000.*depth_t(i,j,k)

! Coef alpha0:
      a0_   =a0_w97+a1_w97*tem_   +a2_w97*sal_
! Coef p0:
      p0_   =b0_w97+b1_w97*tem_   +b2_w97*tem2_   +b3_w97*tem3_      &
                   +b4_w97*sal_   +b5_w97*tem_   *sal_
! Coef lambda:
      l_   = c0_w97+c1_w97*tem_   +c2_w97*tem2_   +c3_w97*tem3_      &
                   +c4_w97*sal_   +c5_w97*tem_   *sal_

      x0=l_   +a0_   *p0_
      anyv3d(i,j,k,2)=l_
      anyv3d(i,j,k,3)=x0**2
      anyv3d(i,j,k,4)=a0_   *x0

! Potential density anomaly:
      rhp_t(i,j,k)=p0_   /(l_   +a0_   *p0_   ) - rho

!     rho_   =anyv3d(i,j,k,2)*prs_        &
!           /(anyv3d(i,j,k,3)+prs_   *anyv3d(i,j,k,4))
!$ Hydrostatic pressure related to the potential density:
       x2=grav*rhp_t(i,j,k)*dz_t(i,j,k,1)                              !02-10-09
       xy_t(i,j,1)=xy_t(i,j,1)+x2                                      !06/03/03
       anyv3d(i,j,k,5)=xy_t(i,j,1)-0.5*x2         ! p(i,j,k)

      enddo
      enddo
      enddo

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=0.
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference: depth_u(i,j,k)

! Hydrostatic pressure at i=i-im:

! Compression termscheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =(-0.25*depth_w(i-im,j,k+1)    &
               -0.25*depth_w(i+ip,j,k+1)    &
                -0.5*depth_u(i   ,j,k  ))*10000.

      a2_   =alpha*anyv3d(i-im,j,k+1,2)+(1.-alpha)*anyv3d(i-im,j,k,2)
      a3_   =alpha*anyv3d(i-im,j,k+1,3)+(1.-alpha)*anyv3d(i-im,j,k,3)
      a4_   =alpha*anyv3d(i-im,j,k+1,4)+(1.-alpha)*anyv3d(i-im,j,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i-im,j,1)=                                                &
                     grav*(depth_w(i-im,j,k+1)                       &
                          -depth_u(i   ,j,k  ))                      &
                                     *rhc_

! Hydrostatic pressure at i=i+ip:
      a2_   =alpha*anyv3d(i+ip,j,k+1,2)+(1.-alpha)*anyv3d(i+ip,j,k,2)
      a3_   =alpha*anyv3d(i+ip,j,k+1,3)+(1.-alpha)*anyv3d(i+ip,j,k,3)
      a4_   =alpha*anyv3d(i+ip,j,k+1,4)+(1.-alpha)*anyv3d(i+ip,j,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i+ip,j,2)=                                                &
                     grav*(depth_w(i+ip,j,k+1)                       &
                          -depth_u(i   ,j,k  ))                      &
                                     *rhc_

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                      &

! pgf compression terms:
             +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1          &

! pgf potential density terms:
             +const1*( anyv3d(i+ip,j,k,5)-anyv3d(i-im,j,k,5)) &
             +const2*(depth_t(i+ip,j,k) -depth_t(i-im,j,k))   &
                    *(  rhp_t(i+ip,j,k)  *  dz_t(i+ip,j,k,1)  &!27/12/07
                       +rhp_t(i-im,j,k)  *  dz_t(i-im,j,k,1)) &
                     /(  dz_t(i+ip,j,k,1) + dz_t(i-im,j,k,1))

      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference: depth_u(i,j,k)

! Hydrostatic pressure at i=i-im:

! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=(depth_u(i,j,k+1)-depth_w(i-im,j,k+1))    &  !09-08-10
           /(depth_u(i,j,k+1)-depth_u(i   ,j,k))

! Compression terms:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =-0.5*(depth_u(i,j,k+1)+depth_u(i,j,k))*10000.

      a2_   =alpha*anyv3d(i-im,j,k+1,2)+(1.-alpha)*anyv3d(i-im,j,k,2)
      a3_   =alpha*anyv3d(i-im,j,k+1,3)+(1.-alpha)*anyv3d(i-im,j,k,3)
      a4_   =alpha*anyv3d(i-im,j,k+1,4)+(1.-alpha)*anyv3d(i-im,j,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i-im,j,1)=                                      &
      xy_t(i-im,j,1)+grav*(depth_u(i,j,k+1)                &
                          -depth_u(i,j,k  ))               &
                                     *rhc_

! Hydrostatic pressure at i=i+ip:

! Coefficient for rho, T, S averaging:
      alpha=(depth_u(i,j,k+1)-depth_w(i+ip,j,k+1))          &
           /(depth_u(i,j,k+1)-depth_u(i   ,j,k))

      a2_   =alpha*anyv3d(i+ip,j,k+1,2)+(1.-alpha)*anyv3d(i+ip,j,k,2)
      a3_   =alpha*anyv3d(i+ip,j,k+1,3)+(1.-alpha)*anyv3d(i+ip,j,k,3)
      a4_   =alpha*anyv3d(i+ip,j,k+1,4)+(1.-alpha)*anyv3d(i+ip,j,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i+ip,j,2)=                                    &
      xy_t(i+ip,j,2)+grav*(depth_u(i,j,k+1)              &
                          -depth_u(i,j,k  ))             &
                                     *rhc_

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                    &

! pgf compression terms:
           +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1          &

! pgf potential density terms:
           +const1*( anyv3d(i+ip,j,k,5)-anyv3d(i-im,j,k,5)) &
           +const2*(depth_t(i+ip,j,k) -depth_t(i-im,j,k))   &
                  *(  rhp_t(i+ip,j,k)  *  dz_t(i+ip,j,k,1)  &
                     +rhp_t(i-im,j,k)  *  dz_t(i-im,j,k,1)) &
                   /(  dz_t(i+ip,j,k,1) + dz_t(i-im,j,k,1))

      enddo
      enddo
      enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging:
      alpha=0.
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference: depth_v

! Hydrostatic pressure at j=j-jm:

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =(-0.25*depth_w(i,j-jm,k+1)    &
               -0.25*depth_w(i,j+jp,k+1)    &
                -0.5*depth_v(i,j   ,k))*10000.

      a2_   =alpha*anyv3d(i,j-jm,k+1,2)+(1.-alpha)*anyv3d(i,j-jm,k,2)
      a3_   =alpha*anyv3d(i,j-jm,k+1,3)+(1.-alpha)*anyv3d(i,j-jm,k,3)
      a4_   =alpha*anyv3d(i,j-jm,k+1,4)+(1.-alpha)*anyv3d(i,j-jm,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j-jm,1)=                                                &
                     grav*(depth_w(i,j-jm,k+1)                       &
                          -depth_v(i,j   ,k  ))                      &
                                     *rhc_

! Hydrostatic pressure at j=j+jp:
      a2_   =alpha*anyv3d(i,j+jp,k+1,2)+(1.-alpha)*anyv3d(i,j+jp,k,2)
      a3_   =alpha*anyv3d(i,j+jp,k+1,3)+(1.-alpha)*anyv3d(i,j+jp,k,3)
      a4_   =alpha*anyv3d(i,j+jp,k+1,4)+(1.-alpha)*anyv3d(i,j+jp,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j+jp,2)=                                              &
                     grav*(depth_w(i,j+jp,k+1)                     &
                          -depth_v(i,j   ,k  ))                    &
                                     *rhc_

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                  &

! pgf compression terms:
          +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1         &

! pgf potential density terms:
         +const1*( anyv3d(i,j+jp,k,5)-anyv3d(i,j-jm,k,5)) &
         +const2*(depth_t(i,j+jp,k) -depth_t(i,j-jm,k))   &
                *(  rhp_t(i,j+jp,k)  *  dz_t(i,j+jp,k,1)  &
                   +rhp_t(i,j-jm,k)  *  dz_t(i,j-jm,k,1)) &
                 /(  dz_t(i,j+jp,k,1)+  dz_t(i,j-jm,k,1))


      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference: depth_v

! Hydrostatic pressure at j=j-jm:

! Coefficient for rho, T, S averaging:
      alpha=(depth_v(i,j,k+1)-depth_w(i,j-jm,k+1))          &
           /(depth_v(i,j,k+1)-depth_v(i,j   ,k))

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      prs_   =-0.5*(depth_v(i,j,k+1)+depth_v(i,j,k))*10000

      a2_   =alpha*anyv3d(i,j-jm,k+1,2)+(1.-alpha)*anyv3d(i,j-jm,k,2)
      a3_   =alpha*anyv3d(i,j-jm,k+1,3)+(1.-alpha)*anyv3d(i,j-jm,k,3)
      a4_   =alpha*anyv3d(i,j-jm,k+1,4)+(1.-alpha)*anyv3d(i,j-jm,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j-jm,1)=                                                &
      xy_t(i,j-jm,1)+grav*(depth_v(i,j,k+1)                         &
                          -depth_v(i,j,k  ))                        &
                                     *rhc_

! Hydrostatic pressure at j=j+jp:

! Coefficient for rho, T, S averaging:
      alpha=(depth_v(i,j,k+1)-depth_w(i,j+jp,k+1))          &
           /(depth_v(i,j,k+1)-depth_v(i,j   ,k))

      a2_   =alpha*anyv3d(i,j+jp,k+1,2)+(1.-alpha)*anyv3d(i,j+jp,k,2)
      a3_   =alpha*anyv3d(i,j+jp,k+1,3)+(1.-alpha)*anyv3d(i,j+jp,k,3)
      a4_   =alpha*anyv3d(i,j+jp,k+1,4)+(1.-alpha)*anyv3d(i,j+jp,k,4)
      rhc_   =a2_   *prs_   /(a3_   +prs_   *a4_   )

      xy_t(i,j+jp,2)=                                               &
      xy_t(i,j+jp,2)+grav*(depth_v(i,j,k+1)                         &
                          -depth_v(i,j,k  ))                        &
                                     *rhc_

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                  &

! pgf compression terms:
          +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1         &

! pgf potential density terms:
         +const1*( anyv3d(i,j+jp,k,5)-anyv3d(i,j-jm,k,5)) &
         +const2*(depth_t(i,j+jp,k) -depth_t(i,j-jm,k))   &
                *(  rhp_t(i,j+jp,k)  *  dz_t(i,j+jp,k,1)  &
                   +rhp_t(i,j-jm,k)  *  dz_t(i,j-jm,k,1)) &
                 /(  dz_t(i,j+jp,k,1)+  dz_t(i,j-jm,k,1))


      enddo
      enddo
      enddo

      end subroutine pressure_gradient_w97_hyb

!-------------------------------------------------------------------------

      subroutine pressure_gradient_m91_hyb
      use module_principal
      implicit none
      integer time_
#ifdef synopsis
       subroutinetitle='pressure_gradient_m91_hyb'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      time_   =now
      call equation_of_state('potential density',now)

      do j=1,jmax
      do i=1,imax
      xy_t(i,j,1)=0.
      enddo
      enddo

      do k=kmax,1,-1
      do j=1,jmax
      do i=1,imax

!$ Hydrostatic pressure related to the potential density:
       x2=grav*rhp_t(i,j,k)*dz_t(i,j,k,1)                              !02-10-09
       xy_t(i,j,1)=xy_t(i,j,1)+x2                                      !06/03/03
       anyv3d(i,j,k,5)=xy_t(i,j,1)-0.5*x2         ! p(i,j,k)

      enddo
      enddo
      enddo


!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=0.
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference: depth_u

! Hydrostatic pressure at i=i-im:

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.25*depth_w(i-im,j,k+1)    &
        +0.25*depth_w(i+ip,j,k+1)    &
         +0.5*depth_u(i   ,j,k)

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i-im,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i-im,j,1)=                                                &
                     grav*(depth_w(i-im,j,k+1)                       &
                          -depth_u(i   ,j,k  ))                      &
                                *(2000.*x0*(5.-x0) )

! Hydrostatic pressure at i=i+ip:
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i+ip,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i+ip,j,2)=                                                &
                     grav*(depth_w(i+ip,j,k+1)                       &
                          -depth_u(i   ,j,k  ))                      &
                                *(2000.*x0*(5.-x0) )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                    &

! pgf compression terms:
           +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1          &

! pgf potential density terms:
           +const1*( anyv3d(i+ip,j,k,5)-anyv3d(i-im,j,k,5)) &
           +const2*(depth_t(i+ip,j,k) -depth_t(i-im,j,k))   &
                  *(  rhp_t(i+ip,j,k)  *  dz_t(i+ip,j,k,1)  &
                     +rhp_t(i-im,j,k)  *  dz_t(i-im,j,k,1)) &
                   /(  dz_t(i+ip,j,k,1) + dz_t(i-im,j,k,1))

      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax-1
      do i=2,imax


! Depth of the horizontal plan of reference: depth_u

! Hydrostatic pressure at i=i-im:

! Coefficient for rho, T, S averaging (alpha=0 if k=kmax):
      alpha=(depth_u(i,j,k+1)-depth_w(i-im,j,k+1))    &  !09-08-10
           /(depth_u(i,j,k+1)-depth_u(i   ,j,k))

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.5*(depth_u(i,j,k+1)+depth_u(i,j,k))
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i-im,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i-im,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i-im,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i-im,j,1)=                                               &
      xy_t(i-im,j,1)+grav*(depth_u(i,j,k+1)                         &
                          -depth_u(i,j,k  ))                        &
                                *(2000.*x0*(5.-x0) )

! Hydrostatic pressure at i=i+ip:

! Coefficient for rho, T, S averaging:
      alpha=(depth_u(i,j,k+1)-depth_w(i+ip,j,k+1))          &
           /(depth_u(i,j,k+1)-depth_u(i   ,j,k))

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i+ip,j,k  ,time_   ) )     &
                     +(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i+ip,j,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i+ip,j,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i+ip,j,2)=                                                &
      xy_t(i+ip,j,2)+grav*(depth_u(i,j,k+1)                        &
                          -depth_u(i,j,k  ))                       &
                                *(2000.*x0*(5.-x0) )

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                    &

! pgf compression terms:
           +(xy_t(i+ip,j,2)-xy_t(i-im,j,1))*const1          &

! pgf potential density terms:
           +const1*( anyv3d(i+ip,j,k,5)-anyv3d(i-im,j,k,5)) &
           +const2*(depth_t(i+ip,j,k) -depth_t(i-im,j,k))   &
                  *(  rhp_t(i+ip,j,k)  *  dz_t(i+ip,j,k,1)  &
                     +rhp_t(i-im,j,k)  *  dz_t(i-im,j,k,1)) &
                   /(  dz_t(i+ip,j,k,1) + dz_t(i-im,j,k,1))

      enddo
      enddo
      enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

! k=kmax, the first level under the surface, is a particular case:
      k=kmax
! Coefficient for rho, T, S averaging:
      alpha=0.
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference: depth_v

! Hydrostatic pressure at j=j-jm:

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.25*depth_w(i,j-jm,k+1)    &
        +0.25*depth_w(i,j+jp,k+1)    &
         +0.5*depth_v(i,j   ,k)

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j-jm,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j-jm,1)=                                              &
                     grav*(depth_w(i,j-jm,k+1)                     &
                          -depth_v(i,j   ,k  ))                    &
                                *(2000.*x0*(5.-x0) )

! Hydrostatic pressure at j=j+jp:
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j+jp,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j+jp,2)=                                              &
                     grav*(depth_w(i,j+jp,k+1)                     &
                          -depth_v(i,j   ,k  ))                    &
                                *(2000.*x0*(5.-x0) )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                  &

! pgf compression terms:
          +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1         &

! pgf potential density terms:
         +const1*( anyv3d(i,j+jp,k,5)-anyv3d(i,j-jm,k,5)) &
         +const2*(depth_t(i,j+jp,k) -depth_t(i,j-jm,k))   &
                *(  rhp_t(i,j+jp,k)  *  dz_t(i,j+jp,k,1)  &
                   +rhp_t(i,j-jm,k)  *  dz_t(i,j-jm,k,1)) &
                 /(  dz_t(i,j+jp,k,1)+  dz_t(i,j-jm,k,1))


      enddo
      enddo

! Now the other vertical levels (from k=kmax-1 to k=1):
      do k=kmax-1,1,-1
      do j=2,jmax
      do i=2,imax-1

! Depth of the horizontal plan of reference: depth_v

! Hydrostatic pressure at j=j-jm:

! Coefficient for rho, T, S averaging:
      alpha=(depth_v(i,j,k+1)-depth_w(i,j-jm,k+1))          &
           /(depth_v(i,j,k+1)-depth_v(i,j,k))

! Compression terms of the Mellor 91 scheme:
! x1, the pressure (decibar), is the same for i=i-im and i+ip in order to
! remove "sigma type" troncation errors
      x1=0.5*(depth_v(i,j,k+1)+depth_v(i,j,k))
      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j-jm,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j-jm,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j-jm,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j-jm,1)=                                               &
      xy_t(i,j-jm,1)+grav*(depth_v(i,j,k+1)                         &
                          -depth_v(i,j,k  ))                        &
                                *(2000.*x0*(5.-x0) )

! Hydrostatic pressure at j=j+jp:

! Coefficient for rho, T, S averaging:
      alpha=(depth_v(i,j,k+1)-depth_w(i,j+jp,k+1))          &
           /(depth_v(i,j,k+1)-depth_v(i,j   ,k))

      x0=-x1/(1402.3                                                &
                +1.34*(     alpha *sal_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*sal_t(i,j+jp,k  ,time_   ) )     &
                     +(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   ) )     &
         *(4.55-0.045*(     alpha *tem_t(i,j+jp,k+1,time_   )       &
                       +(1.-alpha)*tem_t(i,j+jp,k  ,time_   )  ) )  &
       +x1*(15.e-9*x1-0.00821) )**2

      xy_t(i,j+jp,2)=                                               &
      xy_t(i,j+jp,2)+grav*(depth_v(i,j,k+1)                       &
                          -depth_v(i,j,k  ))                      &
                                *(2000.*x0*(5.-x0) )

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                  &

! pgf compression terms:
          +(xy_t(i,j+jp,2)-xy_t(i,j-jm,1))*const1         &

! pgf potential density terms:
         +const1*( anyv3d(i,j+jp,k,5)-anyv3d(i,j-jm,k,5)) &
         +const2*(depth_t(i,j+jp,k) -depth_t(i,j-jm,k))   &
                *(  rhp_t(i,j+jp,k)  *  dz_t(i,j+jp,k,1)  &
                   +rhp_t(i,j-jm,k)  *  dz_t(i,j-jm,k,1)) &
                 /(  dz_t(i,j+jp,k,1)+  dz_t(i,j-jm,k,1))


      enddo
      enddo
      enddo

      end subroutine pressure_gradient_m91_hyb

!..........................................................................

      subroutine pressure_gradient_rhpzavr
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_rhpzavr'
       subroutinedescription= &
       'Computes a sigma weighted vertical average of rhp'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Explications pour rhpzavr dans:
!https://docs.google.com/document/d/1_SdjMHFzI2EVxk17kXwjsXi4ra1nHiU_Wyg2LyL9KZ8/edit?usp=sharing

! Pourquoi la moyenne de la densite potentielle et non pas la densite totale?
! Parce que le terme de compression est calcule en considerant que la pression
! est donnee par la profondeur (depth_u ou depth_v) ce qui a pour effet de
! negliger le rôle de la ssh dans l'effet de compression.... Inutile donc
! de vouloir enlever ce dernier dans rhpzavr

! Si on prend le parti de ce que rhpzavr represente la moyenne de 
! l'anomalie de densite (densite potentielle -rho) cela
! devrait permettre de revenir automatiquement a l'ancienne methode
! en specifiant simplement rhzavr=0. Dans le cas general la densité potentielle
! totale est donc: rhp+rho+rhzavr


! Etape 1: calculer rhpzavr_w
! Methode 1 !15-09-16
      do j=1,jmax ; do i=1,imax
       rhpzavr_w(i,j)=0.
      enddo       ; enddo
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhpzavr_w(i,j)=rhpzavr_w(i,j)                   &
                       +(depth_t(i,j,k)-depth_w(i,j,1)) &
                          *rhp_t(i,j,k)                 &
                           *dz_t(i,j,k,1)
      enddo       ; enddo      ; enddo
      do j=1,jmax ; do i=1,imax
       rhpzavr_w(i,j)=rhpzavr_w(i,j)*2./hz_w(i,j,1) & ! OUI OUI IL FAUT LA MULTIPLICATION PAR 2
                      /(depth_w(i,j,kmaxp1)-depth_w(i,j,1))
      enddo       ; enddo

! Methode 2 (equivalente A methode 1)
!     do j=1,jmax
!     do i=1,imax
!      sum1=0.
!      do k=kmin_w(i,j),kmax
!       sum1=sum1+(depth_t(i,j,k)     -depth_w(i,j,kmin_w(i,j))) &
!                /(depth_w(i,j,kmaxp1)-depth_w(i,j,kmin_w(i,j))) &
!                  * rhp_t(i,j,k)                                &
!                    *dz_t(i,j,k,1)                               
!      enddo
!      rhpzavr_w(i,j)=2.*sum1/hz_w(i,j,1) ! OUI OUI IL FAUT LA MULTIPLICATION PAR 2
!      if(par%rank==114.and.mask_t(i,j,kmax)==1)write(65,*)rhpzavr_w(i,j),2.*sum1/hz_w(i,j,1)
!     enddo
!     enddo

! Etape 2: soustraire A rhp_t:
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhp_t(i,j,k)=rhp_t(i,j,k)-rhpzavr_w(i,j)
      enddo       ; enddo      ; enddo

      end subroutine pressure_gradient_rhpzavr

!................................................................
#ifdef bidon
      subroutine pressure_gradient_rhpzavrstep2
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_rhpzavrstep2'
       subroutinedescription=''
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! rhzavr (if any) contribution:
      const2=grav/(rho)*flag3d 
      do k=kmax,1,-1
      do j=2,jmax-1
      do i=2,imax

       presgrad_u(i,j,k,1)=                                       &
       presgrad_u(i,j,k,1)                                        &
            -const2*depth_u(i,j,k)*( rhpzavr_w(i+ipu,j)         &
                                    -rhpzavr_w(i-imu,j))

      enddo
      enddo
      enddo
      do k=kmax,1,-1
      do j=2,jmax
      do i=2,imax-1
       presgrad_v(i,j,k,1)=                                 &
       presgrad_v(i,j,k,1)                                  &
         -const2*depth_v(i,j,k)*( rhpzavr_w(i,j+jpv)      &
                                 -rhpzavr_w(i,j-jmv))
      enddo
      enddo
      enddo
      end subroutine pressure_gradient_rhpzavrstep2
#endif
!...............................................................
      subroutine pressure_gradient_initref !14-09-15
      use module_principal
      use module_parallele
      implicit none
      integer loop_
      real dtem_,dsal_

!     allocate(trefco_t(0:imax+1,0:jmax+1,0:dim_ref)) ; trefco_t=0.
!     allocate(srefco_t(0:imax+1,0:jmax+1,0:dim_ref)) ; srefco_t=0.

! Un etat de reference 1DV z (independant de x et y) de densite de compression
      if(eos_comprs==1) then !---->
       id_tem=0 ; id_sal=1
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_tem)=t0
        anyv3d(i,j,k,id_sal)=s0
       enddo ; enddo ; enddo
      ! Calculer rhcref
       id_rhc=3
       call equation_of_state_prs_anyv(1,imax,1,jmax)
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhcref_t(i,j,k)=anyv3d(i,j,k,id_rhc)
        !write(75,*)depth_t(i,j,k),rhcref_t(i,j,k)
       enddo       ; enddo       ; enddo
      else                   !---->
       rhcref_t=0.
      endif                  !---->
      if(flag_refstate==0)return !18-11-15 rhcref_t est calcule meme si flag_refstate=0

#ifdef bidonref

! continuite mpi + obc pour temref et salref: (ce serait mieux de le faire dans module_ogcm dans la foulee
!                                              de l'interpolation de temref et salref)
      call pressure_gradient_ts_ref

! DENSITE POTENTIELLE DE REFERENCE ASSOCIEE
       id_rhp=2
       call equation_of_state_pot_anyv(1,imax,1,jmax)
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhpref_t(i,j,k)=anyv3d(i,j,k,id_rhp)
       enddo       ; enddo       ; enddo

! CALCULER rhpref_u et rhpref_v annulant le PGF pour la densite potentielle de reference
       id_prs=4
       call pressure_gradient_cancelpgf(0) ! 0=pas d'effet de compression

! Obtenir temref_u salref_u temref_v salref_v correspondant aux potentielles rhpref_u rhpref_v
! en inversant l'EOS(rhpref_u,rhpref_v) autour de temref_t et de salref_t:
       call pressure_gradient_eosinverse

! CALCULER rhpref_u et rhpref_v annulant le PGF pour la densite de reference complete
      if(eos_comprs==1) then !>>>>>>>>
       id_tem=0 ; id_sal=1 ; id_rhf=2
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_tem)=temref_t(i,j,k)
        anyv3d(i,j,k,id_sal)=salref_t(i,j,k)
       enddo ; enddo ; enddo
       call equation_of_state_full_anyv(1,imax,1,jmax) ! (inclut la soustraction de l'etat 1Dz ref de comprs)
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
!       anyv3d(i,j,k,id_rhf)=anyv3d(i,j,k,id_rhf)-rhcref_t(i,j,k) ! ote l'etat 1Dz de compression
             rhpref_t(i,j,k)=anyv3d(i,j,k,id_rhf)
       enddo ; enddo ; enddo
! CALCULER rhpref_u et rhpref_v annulant le PGF pour la densite de reference complete
       id_prs=4
       call pressure_gradient_cancelpgf(1)  ! 1=effet de compression active
      endif                  !>>>>>>>>

!     deallocate(trefco_t)
!     deallocate(srefco_t)

#endif
      end subroutine pressure_gradient_initref
!...................................................................................
#ifdef bidonref
      subroutine pressure_gradient_cancelpgf(flag_)
      use module_principal
      use module_parallele
      implicit none
      integer :: id_rhorefuv_=0 &
                ,id_prsref_=5   &
                ,flag_
 
! Cette routine calcule rhpref_u et rhpref_v pour annuler le PGF
! Le tableau anyv3d(i,j,k,id_rho_) envisage 2 types de densite de
! reference: une densite potentielle de reference et une densite
! de compression de reference


! Un champ "2 tableaux" de densite donnant un champ de densite localement
! identique au champs de reference mais sans le gradient horizontal
      id_tem=0 ; id_sal=1 ; id_rhom=2 ; id_rhop=3
      call pressure_gradient_homogeneousref('u',flag_)

! CALCULER rhpref_u
      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=const1
      x0=0.01 ! est un delta z critique qui prend le relais du gradient de depth si celui ci est trop petit
!     do k=1,kmax
      do j=2,jmax-1
      do i=2,imax

      anyv3d(i+ipu,j,kmax+1,id_prs)=0.
      anyv3d(i-imu,j,kmax+1,id_prs)=0.
      anyv3d(i+ipu,j,kmax+1,id_prsref_)=0.
      anyv3d(i-imu,j,kmax+1,id_prsref_)=0.

      do k=kmax,1,-1

! On calcule 2 fois la pression (en i+ipu et en i-imu) pour pouvoir envisager le dedoublement
! d'une meme colonne de densite de reference lorsque l'etat de reference a une dependance vis
! a vis de x ou y

           anyv3d(i+ipu,j,k  ,id_prs)  &
          =anyv3d(i+ipu,j,k+1,id_prs)  &
       +grav*dz_t(i+ipu,j,k,1)         &
          *anyv3d(i+ipu,j,k,id_rhop)

           anyv3d(i+ipu,j,k  ,id_prsref_)  &
          =anyv3d(i+ipu,j,k+1,id_prsref_)  &
       +grav*dz_t(i+ipu,j,k,1)             &
      *(   anyv3d(i+ipu,j,k,id_rhop)       &
        -rhpref_t(i+ipu,j,k)        )

           anyv3d(i-imu,j,k  ,id_prs)  &
          =anyv3d(i-imu,j,k+1,id_prs)  &
       +grav*dz_t(i-imu,j,k,1)         &
          *anyv3d(i-imu,j,k,id_rhom)

           anyv3d(i-imu,j,k  ,id_prsref_)  &
          =anyv3d(i-imu,j,k+1,id_prsref_)  &
       +grav*dz_t(i-imu,j,k,1)             &
      *(   anyv3d(i-imu,j,k,id_rhom)       &
        -rhpref_t(i-imu,j,k)        )


      enddo


! k0 est le Niveau de la couche merged:
      k0=max(kmin_w(i,j),kmin_w(i-1,j))
      do k=kmax,1,-1

! x0 est une echelle  en dessous de laquelle on considere que dz/di est trop petit
! pour qu'il y ait de l'erreur a corriger. Un "petit" ecart peut etre defini comme
! une fraction de la resolution verticale par exemple x0=0.1*dz/dk  
!     x0= 0.01*dz_u(i,j,k,1)
      x0=0.001*dz_u(i,j,k,1)

      if(k>k0) then !>>>
       kp1=k+1  ; km1=k  ! Cas normal
      else          !>>>
       kp1=k0+1 ; km1=1  ! Cas couche merged entre 1 et k0+1
      endif         !>>>

      x1=0.5*(depth_w(i+ipu,j,kp1)-depth_w(i-imu,j,kp1)  & ! <dz/di>
             +depth_w(i+ipu,j,km1)-depth_w(i-imu,j,km1) )

      x2= depth_w(i+ipu,j,kp1)-depth_w(i+ipu,j,km1)   &      ! <dz/dk>
         +depth_w(i-imu,j,kp1)-depth_w(i-imu,j,km1)    

! Voici la formule qu'on utilise pour eviter la singularite:
! rhp=(dz**2 * formule_de_base + x0**2 * rhp_brute ) / (dz**2+x0**2)
                  anyv3d(i,j,k,id_rhorefuv_)=(                        &

      -x1*( &

              0.5*( anyv3d(i+ipu,j,kp1,id_prs)-anyv3d(i-imu,j,kp1,id_prs)  &
                   +anyv3d(i+ipu,j,km1,id_prs)-anyv3d(i-imu,j,km1,id_prs)) &


             +x1*(         &

         ( anyv3d(i+ipu,j,km1,id_prsref_)-anyv3d(i+ipu,j,kp1,id_prsref_)   &
          +anyv3d(i-imu,j,km1,id_prsref_)-anyv3d(i-imu,j,kp1,id_prsref_) ) &
                 /x2   &

           )                    &

                          )/grav                  &

      +(x0**2)*(anyv3d(i+ipu,j,km1,id_prs)-anyv3d(i+ipu,j,kp1,id_prs)  &
               +anyv3d(i-imu,j,km1,id_prs)-anyv3d(i-imu,j,kp1,id_prs))/grav &
             /x2   &
                                      )                             &
            /(x0**2+x1**2)

        x3=min(anyv3d(i-imu,j,k,id_rhom),anyv3d(i+ipu,j,k,id_rhop))
        x4=max(anyv3d(i-imu,j,k,id_rhom),anyv3d(i+ipu,j,k,id_rhop))
! cas particulier de la couche merged, la borne max A considerer est celle des
! niveaux les plus profonds 'donc k=1)
!       if(k<=k0)x4=max(anyv3d(i-imu,j,1,id_rhom),anyv3d(i+ipu,j,1,id_rhop))
!       if(k<=k0) then !zzz>
! couche merged = first gusse sans condition:
!         anyv3d(i,j,k,id_rhorefuv_)=                            & 
!        (anyv3d(i+ipu,j,km1,id_prs)-anyv3d(i+ipu,j,kp1,id_prs)  &
!        +anyv3d(i-imu,j,km1,id_prs)-anyv3d(i-imu,j,kp1,id_prs))/grav/x2
!       endif          !zzz>

! Par rapport  au pb de la singularite on peut envisager 2 approches
! Approche 1: borner aux valeurs extremes autorisees:
        anyv3d(i,j,k,id_rhorefuv_)=min(max(anyv3d(i,j,k,id_rhorefuv_),x3),x4) !Approche 1
! Approche 2: considerer que revenir au first guess est la solution la plus sure
! voir parfois plus precise
!      if(anyv3d(i,j,k,id_rhorefuv_)>x4.or. &
!         anyv3d(i,j,k,id_rhorefuv_)<x3) then !pmxpxm>

!         anyv3d(i,j,k,id_rhorefuv_)=                            & !Approche 2
!        (anyv3d(i+ipu,j,km1,id_prs)-anyv3d(i+ipu,j,kp1,id_prs)  &
!        +anyv3d(i-imu,j,km1,id_prs)-anyv3d(i-imu,j,kp1,id_prs))/grav/x2

!      endif                                  !pmxpmx> 

       enddo
       enddo
       enddo

!...........
! Sur les bords:
      if(obcstatus(ieq1)==1) then ; do k=1,kmax ; do j=2,jmax-1
         anyv3d(1,j,k,0)=anyv3d(2,j,k,0)
      enddo ; enddo ; endif
      if(obcstatus(ieqimax)==1) then ; do k=1,kmax ; do j=2,jmax-1
         anyv3d(imax+1,j,k,0)=anyv3d(imax,j,k,0)
      enddo ; enddo ; endif
      if(obcstatus(jeq1)==1) then ; do k=1,kmax ; do i=1,imax+1
         anyv3d(i,1,k,0)=anyv3d(i,2,k,0)
      enddo ; enddo ; endif
      if(obcstatus(jeqjmax)==1) then ; do k=1,kmax ; do i=1,imax+1
         anyv3d(i,jmax,k,0)=anyv3d(i,jmax-1,k,0)
      enddo ; enddo ; endif
!...........
! Conditions mpi:
       call obc_int_anyv3d(id_rhorefuv_,'u1')
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        rhpref_u(i,j,k)=anyv3d(i,j,k,id_rhorefuv_)
       enddo ; enddo ; enddo

!.................................................................
! CALCULER rhpref_v
      id_tem=0 ; id_sal=1 ; id_rhom=2 ; id_rhop=3 ; id_tem2=5 ; id_sal2=6
      call pressure_gradient_homogeneousref('v',flag_)

! CALCULER rhpref_v
      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=const1
      x0=0.01 ! est un delta z critique qui prend le relais du gradient de depth si celui ci est trop petit
      do j=2,jmax
      do i=2,imax-1

      anyv3d(i,j+jpv,kmax+1,id_prs)=0.
      anyv3d(i,j-jmv,kmax+1,id_prs)=0.
      anyv3d(i,j+jpv,kmax+1,id_prsref_)=0.
      anyv3d(i,j-jmv,kmax+1,id_prsref_)=0.

      do k=kmax,1,-1

! On calcule 2 fois la pression (en i+ipu et en i-imu) pour pouvoir envisager le dedoublement
! d'une meme colonne de densite de reference lorsque l'etat de reference a une dependance vis
! a vis de x ou y

           anyv3d(i,j+jpv,k  ,id_prs)  &
          =anyv3d(i,j+jpv,k+1,id_prs)  &
       +grav*dz_t(i,j+jpv,k,1)         &
          *anyv3d(i,j+jpv,k,id_rhop)

           anyv3d(i,j+jpv,k  ,id_prsref_)  &
          =anyv3d(i,j+jpv,k+1,id_prsref_)  &
       +grav*dz_t(i,j+jpv,k,1)             &
      *(   anyv3d(i,j+jpv,k,id_rhop)       &
        -rhpref_t(i,j+jpv,k)        )

           anyv3d(i,j-jmv,k  ,id_prs)  &
          =anyv3d(i,j-jmv,k+1,id_prs)  &
       +grav*dz_t(i,j-jmv,k,1)         &
          *anyv3d(i,j-jmv,k,id_rhom)

           anyv3d(i,j-jmv,k  ,id_prsref_)  &
          =anyv3d(i,j-jmv,k+1,id_prsref_)  &
       +grav*dz_t(i,j-jmv,k,1)             &
      *(   anyv3d(i,j-jmv,k,id_rhom)       &
        -rhpref_t(i,j-jmv,k)        )

      enddo


! k0 est le Niveau de la couche merged:
      k0=max(kmin_w(i,j),kmin_w(i,j-1))
      do k=kmax,1,-1

! x0 est une echelle  en dessous de laquelle on considere que dz/di est trop petit
! pour qu'il y ait de l'erreur a corriger. Un "petit" ecart peut etre defini comme
! une fraction de la resolution verticale par exemple x0=0.1*dz/dk  
!     x0= 0.01*dz_v(i,j,k,1)
      x0=0.001*dz_v(i,j,k,1)

      if(k>k0) then !>>>
       kp1=k+1  ; km1=k  ! Cas normal
      else          !>>>
       kp1=k0+1 ; km1=1  ! Cas couche merged entre 1 et k0+1
      endif         !>>>

      x1=0.5*(depth_w(i,j+jpv,kp1)-depth_w(i,j-jmv,kp1)  & ! <dz/di>
             +depth_w(i,j+jpv,km1)-depth_w(i,j-jmv,km1) )

      x2= depth_w(i,j+jpv,kp1)-depth_w(i,j+jpv,km1)   &      ! <dz/dk>
         +depth_w(i,j-jmv,kp1)-depth_w(i,j-jmv,km1)    

! Voici la formule qu'on utilise pour eviter la singularite:
! rhp=(dz**2 * formule_de_base + x0**2 * rhp_brute ) / (dz**2+x0**2)
                  anyv3d(i,j,k,id_rhorefuv_)=(                        &

      -x1*( &

              0.5*( anyv3d(i,j+jpv,kp1,id_prs)-anyv3d(i,j-jmv,kp1,id_prs)  &
                   +anyv3d(i,j+jpv,km1,id_prs)-anyv3d(i,j-jmv,km1,id_prs)) &


             +x1*(         &

         ( anyv3d(i,j+jpv,km1,id_prsref_)-anyv3d(i,j+jpv,kp1,id_prsref_)   &
          +anyv3d(i,j-jmv,km1,id_prsref_)-anyv3d(i,j-jmv,kp1,id_prsref_) ) &
                 /x2   &

           )                    &

                          )/grav                  &

      +(x0**2)*(anyv3d(i,j+jpv,km1,id_prs)-anyv3d(i,j+jpv,kp1,id_prs)  &
               +anyv3d(i,j-jmv,km1,id_prs)-anyv3d(i,j-jmv,kp1,id_prs))/grav &
             /x2   &
                                      )                             &
            /(x0**2+x1**2)

        x3=min(anyv3d(i,j-jmv,k,id_rhom),anyv3d(i,j+jpv,k,id_rhop))
        x4=max(anyv3d(i,j-jmv,k,id_rhom),anyv3d(i,j+jpv,k,id_rhop))
! cas particulier de la couche merged, la borne max A considerer est celle des
! niveaux les plus profonds 'donc k=1)
!       if(k<=k0)x4=max(anyv3d(i,j-jmv,1,id_rhom),anyv3d(i,j+jpv,1,id_rhop))
!       if(k<=k0) then !zzz>
! couche merged = first guess sans condition:
!         anyv3d(i,j,k,id_rhorefuv_)=                            & 
!        (anyv3d(i,j+jpv,km1,id_prs)-anyv3d(i,j+jpv,kp1,id_prs)  &
!        +anyv3d(i,j-jmv,km1,id_prs)-anyv3d(i,j-jmv,kp1,id_prs))/grav/x2
!       endif          !zzz>

! Par rapport  au pb de la singularite on peut envisager 2 approches
! Approche 1: borner aux valeurs extremes autorisees:
        anyv3d(i,j,k,id_rhorefuv_)=min(max(anyv3d(i,j,k,id_rhorefuv_),x3),x4) !Approche 1
! Approche 2: considerer que revenir au first guess est la solution la plus sure
! voir parfois plus precise
!      if(anyv3d(i,j,k,id_rhorefuv_)>x4.or. &
!         anyv3d(i,j,k,id_rhorefuv_)<x3) then !pmxpxm>

!         anyv3d(i,j,k,id_rhorefuv_)=                            & !Approche 2
!        (anyv3d(i,j+jpv,km1,id_prs)-anyv3d(i,j+jpv,kp1,id_prs)  &
!        +anyv3d(i,j-jmv,km1,id_prs)-anyv3d(i,j-jmv,kp1,id_prs))/grav/x2

!      endif                                  !pmxpmx> 

      enddo
      enddo
      enddo

!...........
! Sur les bords:
      if(obcstatus(jeq1)==1) then ; do k=1,kmax ; do i=2,imax-1
         anyv3d(i,1,k,0)=anyv3d(i,2,k,0)
      enddo ; enddo ; endif
      if(obcstatus(jeqjmax)==1) then ; do k=1,kmax ; do i=2,imax-1
         anyv3d(i,jmax+1,k,0)=anyv3d(i,jmax,k,0)
      enddo ; enddo ; endif
      if(obcstatus(ieq1)==1) then ; do k=1,kmax ; do j=1,jmax+1
         anyv3d(1,j,k,0)=anyv3d(2,j,k,0)
      enddo ; enddo ; endif
      if(obcstatus(ieqimax)==1) then ; do k=1,kmax ; do j=1,jmax+1
         anyv3d(imax,j,k,0)=anyv3d(imax-1,j,k,0)
      enddo ; enddo ; endif
!...........
! Continuite mpi:
       call obc_int_anyv3d(id_rhorefuv_,'v1')
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        rhpref_v(i,j,k)=anyv3d(i,j,k,id_rhorefuv_)
       enddo ; enddo ; enddo

      end subroutine pressure_gradient_cancelpgf
#endif
!...................................................................
#ifdef bidonref
      subroutine pressure_gradient_eosinverse
      use module_principal
      use module_parallele
      implicit none
      integer loop_
      real dtem_,dsal_

!...................................................................
! Trouver temref_u et salref_u autour de tem et sal first guess local
       id_tem=0 ; id_sal=1 ; id_rhp=2

       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        anyv3d(i,j,k,id_tem)=( dz_t(i  ,j,k,1)*(temref_t(i  ,j,k))  &
                              +dz_t(i-1,j,k,1)*(temref_t(i-1,j,k))) &
                            /( dz_t(i  ,j,k,1)                      &
                              +dz_t(i-1,j,k,1))
        anyv3d(i,j,k,id_sal)=( dz_t(i  ,j,k,1)*(salref_t(i  ,j,k))  &
                              +dz_t(i-1,j,k,1)*(salref_t(i-1,j,k))) &
                            /( dz_t(i  ,j,k,1)                      &
                              +dz_t(i-1,j,k,1))
!     if(i==238.and.j==1.and.k==1)write(6,*)'T S REF',anyv3d(i,j,k,id_tem),anyv3d(i,j,k,id_sal)
!     if(i==238.and.j==1.and.k==1)write(6,*)'rhpref_u',rhpref_u(i,j:j+1,k)
       enddo ; enddo ; enddo

       dtem_=0.1
       dsal_=0.03
       do loop_=1,3  !llll>

       id_rhp=2
       call equation_of_state_pot_anyv(1,imax+1,1,jmax)

!...........
! Pour le moment seule une condition mpi a ete appliquee sur rhpref_u et rhpref_v
! on applique une condition de bords de domaine sans consequence pour le PGF mais
! qui a de l'importance pour l'advection de temref_u salref_u. On choisit arbitrairement
! rhpref_u = < rhpref_t > qui conduira a temref_u = < temref_t > etc...
      if(obcstatus(ieq1)==1) then ; do k=1,kmax ; do j=2,jmax-1
         rhpref_u(1,j,k)=anyv3d(1,j,k,id_rhp)
      enddo ; enddo ; endif
      if(obcstatus(ieqimax)==1) then ; do k=1,kmax ; do j=2,jmax-1
         rhpref_u(imax+1,j,k)=anyv3d(imax+1,j,k,id_rhp)
      enddo ; enddo ; endif
      if(obcstatus(jeq1)==1) then ; do k=1,kmax ; do i=1,imax+1
         rhpref_u(i,1,k)=anyv3d(i,1,k,id_rhp)
      enddo ; enddo ; endif
      if(obcstatus(jeqjmax)==1) then ; do k=1,kmax ; do i=1,imax+1
         rhpref_u(i,jmax,k)=anyv3d(i,jmax,k,id_rhp)
      enddo ; enddo ; endif
!     i=258-par%timax(1)
!     if(i>1.and.i<imax) then
!     write(6,*)'riri1',par%rank,real(rhpref_u(i,1:3,kmax)+rho-1000.)
!     endif
!...........

! Ajouter une perturbation sur tem et recalculer l'eos:
       id_rhp=3
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        anyv3d(i,j,k,id_tem)=anyv3d(i,j,k,id_tem)+dtem_
       enddo ; enddo ; enddo
       call equation_of_state_pot_anyv(1,imax+1,1,jmax)

! Oter la perturbation de tem et ajouter une perturbation sur sal
! et recalculer l'eos
       id_rhp=4
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        anyv3d(i,j,k,id_tem)=anyv3d(i,j,k,id_tem)-dtem_
        anyv3d(i,j,k,id_sal)=anyv3d(i,j,k,id_sal)+dsal_
       enddo ; enddo ; enddo
       call equation_of_state_pot_anyv(1,imax+1,1,jmax)

! Oter la perturbation de sal
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        anyv3d(i,j,k,id_sal)=anyv3d(i,j,k,id_sal)-dsal_
       enddo ; enddo ; enddo


       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1

! dT=x1.dV dS=x2.dV dV=drhp/(x1**2+x2**2)
         x1=(anyv3d(i,j,k,3)-anyv3d(i,j,k,2))/dtem_ ! coef contraction tem
         x2=(anyv3d(i,j,k,4)-anyv3d(i,j,k,2))/dsal_ ! coef contraction sal

! Deduire temref_u et salref_u

! Methode des "coef de contraction" dT=x1.dV dS=x2.dV dV=drhp/(x1**2+x2**2)
!     temref_u(i,j,k)=anyv3d(i,j,k,id_tem)+x1/(x1**2+x2**2+small1)*(rhpref_u(i,j,k)-anyv3d(i,j,k,2))
!     salref_u(i,j,k)=anyv3d(i,j,k,id_sal)+x2/(x1**2+x2**2+small1)*(rhpref_u(i,j,k)-anyv3d(i,j,k,2))
! Methode des "variations locales"
! dS/dT=-abs(VariationS)/abs(VariationT)
!      x4=abs(   sal_t(i,j,k,1) -sal_t(i-1,j,k,1))
       x4=abs(salref_t(i,j,k)-salref_t(i-1,j,k))
!      x3=abs(   tem_t(i,j,k,1) -tem_t(i-1,j,k,1))
       x3=abs(temref_t(i,j,k)-temref_t(i-1,j,k))

      temref_u(i,j,k)=anyv3d(i,j,k,id_tem)-x3/(-x1*x3+x2*x4+small1)*(rhpref_u(i,j,k)-anyv3d(i,j,k,2))
      salref_u(i,j,k)=anyv3d(i,j,k,id_sal)+x4/(-x1*x3+x2*x4+small1)*(rhpref_u(i,j,k)-anyv3d(i,j,k,2))


!     if(i==238.and.j==1.and.k==1)write(6,*)'delta',rhpref_u(i,j,k),anyv3d(i,j,k,2)
!     if(i==238.and.j==1.and.k==1)write(6,*)'temref_u(i,j,k)',temref_u(i,j,k),anyv3d(i,j,k,id_tem)
!     if(i==238.and.j==1.and.k==1)write(6,*)'salref_u(i,j,k)',salref_u(i,j,k),anyv3d(i,j,k,id_sal)

! Pour une autre iteration:
!        write(65,*)anyv3d(i,j,k,id_tem)-temref_u(i,j,k),anyv3d(i,j,k,id_sal)-salref_u(i,j,k)
         anyv3d(i,j,k,id_tem)=temref_u(i,j,k)
         anyv3d(i,j,k,id_sal)=salref_u(i,j,k)

       enddo ; enddo ; enddo

       enddo ! loop_

! retrouver la coherence entre tem et rhp:
       id_rhp=4
       call equation_of_state_pot_anyv(1,imax+1,1,jmax)

!      sum1=0.
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
!      sum1=sum1+mask_u(i,j,k)*abs( rhpref_u(i,j,k) - anyv3d(i,j,k,4) )
!       write(64,*)rhpref_u(i,j,k),anyv3d(i,j,k,2),anyv3d(i,j,k,4)
        rhpref_u(i,j,k)=anyv3d(i,j,k,id_rhp)
!       if(i==imax.and.j==jmax/2)   &
!       write(6,*)depth_u(i,j,k),rhpref_u(i,j,k),temref_u(i,j,k)
       enddo ; enddo ; enddo

!...................................................................
! Trouver temref_v et salref_v autour de tem et sal first guess local
!      id_tem=4 
       id_sal=1 ; id_rhp=2
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        anyv3d(i,j,k,id_tem)=( dz_t(i,j  ,k,1)*(temref_t(i,j  ,k))  &
                              +dz_t(i,j-1,k,1)*(temref_t(i,j-1,k))) &
                            /( dz_t(i,j  ,k,1)                      &
                              +dz_t(i,j-1,k,1))
        anyv3d(i,j,k,id_sal)=( dz_t(i,j  ,k,1)*(salref_t(i,j  ,k))  &
                              +dz_t(i,j-1,k,1)*(salref_t(i,j-1,k))) &
                            /( dz_t(i,j  ,k,1)                      &
                              +dz_t(i,j-1,k,1))
       enddo ; enddo ; enddo

       dtem_=0.1
       dsal_=0.03
       do loop_=1,3  !llll>

       id_rhp=2
       call equation_of_state_pot_anyv(1,imax,1,jmax+1)

!...........
! Pour le moment seule une condition mpi a ete appliquee sur rhpref_v et rhpref_v
! on applique une condition de bords de domaine sans consequence pour le PGF mais
! qui a de l'importance pour l'advection de temref_v salref_v. On choisit arbitrairement
! rhpref_v = < rhpref_t > qui conduira a temref_v = < temref_t > etc...
      if(obcstatus(jeq1)==1) then ; do k=1,kmax ; do i=2,imax-1
         rhpref_v(i,1,k)=anyv3d(i,1,k,id_rhp)
      enddo ; enddo ; endif
      if(obcstatus(jeqjmax)==1) then ; do k=1,kmax ; do i=2,imax-1
         rhpref_v(i,jmax+1,k)=anyv3d(i,jmax+1,k,id_rhp)
      enddo ; enddo ; endif
      if(obcstatus(ieq1)==1) then ; do k=1,kmax ; do j=1,jmax+1
         rhpref_v(1,j,k)=anyv3d(1,j,k,id_rhp)
      enddo ; enddo ; endif
      if(obcstatus(ieqimax)==1) then ; do k=1,kmax ; do j=1,jmax+1
         rhpref_v(imax,j,k)=anyv3d(imax,j,k,id_rhp)
      enddo ; enddo ; endif
!...........

! Ajouter une perturbation sur tem et recalculer l'eos:
       id_rhp=3
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        anyv3d(i,j,k,id_tem)=anyv3d(i,j,k,id_tem)+dtem_
       enddo ; enddo ; enddo
       call equation_of_state_pot_anyv(1,imax,1,jmax+1)

! Oter la perturbation de tem et ajouter une perturbation sur sal
! et recalculer l'eos
       id_rhp=4
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        anyv3d(i,j,k,id_tem)=anyv3d(i,j,k,id_tem)-dtem_
        anyv3d(i,j,k,id_sal)=anyv3d(i,j,k,id_sal)+dsal_
       enddo ; enddo ; enddo
       call equation_of_state_pot_anyv(1,imax,1,jmax+1)

! Oter la perturbation de sal
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
        anyv3d(i,j,k,id_sal)=anyv3d(i,j,k,id_sal)-dsal_
       enddo ; enddo ; enddo


       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax

! dT=x1.dV dS=x2.dV dV=drhp/(x1**2+x2**2)
         x1=(anyv3d(i,j,k,3)-anyv3d(i,j,k,2))/dtem_ ! coef contraction tem
         x2=(anyv3d(i,j,k,4)-anyv3d(i,j,k,2))/dsal_ ! coef contraction sal

! Deduire temref_v et salref_v

! Methode des "coef de contraction" dT=x1.dV dS=x2.dV dV=drhp/(x1**2+x2**2)
!     temref_v(i,j,k)=anyv3d(i,j,k,id_tem)+x1/(x1**2+x2**2+small1)*(rhpref_v(i,j,k)-anyv3d(i,j,k,2))
!     salref_v(i,j,k)=anyv3d(i,j,k,id_sal)+x2/(x1**2+x2**2+small1)*(rhpref_v(i,j,k)-anyv3d(i,j,k,2))
! Methode des "variations locales"
! dS/dT=-abs(VariationS)/abs(VariationT)
!      x4=abs(   sal_t(i,j,k,1) -sal_t(i,j-1,k,1))
       x4=abs(salref_t(i,j,k)-salref_t(i,j-1,k))
!      x3=abs(   tem_t(i,j,k,1) -tem_t(i,j-1,k,1))
       x3=abs(temref_t(i,j,k)-temref_t(i,j-1,k))

      temref_v(i,j,k)=anyv3d(i,j,k,id_tem)-x3/(-x1*x3+x2*x4+small1)*(rhpref_v(i,j,k)-anyv3d(i,j,k,2))
      salref_v(i,j,k)=anyv3d(i,j,k,id_sal)+x4/(-x1*x3+x2*x4+small1)*(rhpref_v(i,j,k)-anyv3d(i,j,k,2))

! Pour une autre iteration:
!        write(65,*)anyv3d(i,j,k,id_tem)-temref_v(i,j,k),anyv3d(i,j,k,id_sal)-salref_v(i,j,k)
         anyv3d(i,j,k,id_tem)=temref_v(i,j,k)
         anyv3d(i,j,k,id_sal)=salref_v(i,j,k)

       enddo ; enddo ; enddo

       enddo ! loop_

! retrouver la coherence entre tem et rhp:
       id_rhp=4
       call equation_of_state_pot_anyv(1,imax,1,jmax+1)

!      sum1=0.
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!      sum1=sum1+mask_v(i,j,k)*abs( rhpref_v(i,j,k) - anyv3d(i,j,k,4) )
!       write(64,*)rhpref_v(i,j,k),anyv3d(i,j,k,2),anyv3d(i,j,k,4)
        rhpref_v(i,j,k)=anyv3d(i,j,k,id_rhp)
       enddo ; enddo ; enddo

      end subroutine pressure_gradient_eosinverse
#endif
!.....................................................................
#ifdef bidon
      subroutine pressure_gradient_coref
      use module_principal ; use module_parallele ; use module_s
      implicit none
      double precision, dimension(:,:), allocatable :: matrix_,invrsmat_,matrix2_
      double precision, dimension(:)  , allocatable :: mat_rhs_,mat_sol_,mat_rhs2_
      integer loop_,ligne_,colonne_,loop2_
      real scz2_,scz4_,scz6_,scz8_
      real delta_

! First guess
        trefco_t(:,:,0)=1.
        trefco_t(:,:,1)=1.
        trefco_t(:,:,2)=1./740. ! 1./400.
        trefco_t(:,:,3)=1.
        trefco_t(:,:,4)=1./355. ! 1./200.
        trefco_t(:,:,5)=1.
        trefco_t(:,:,6)=1./55.  ! 1./100.
        trefco_t(:,:,7)=1.
        trefco_t(:,:,8)=1./40   ! 1./50.

        srefco_t(:,:,0)=1.
        srefco_t(:,:,1)=1.
        srefco_t(:,:,2)=1./620. ! 1./400.
        srefco_t(:,:,3)=1.
        srefco_t(:,:,4)=1./375. ! 1./200.
        srefco_t(:,:,5)=1.
        srefco_t(:,:,6)=1./55.  ! 1./100
        srefco_t(:,:,7)=1.
        srefco_t(:,:,8)=1./40.


      allocate(matrix_  (1+dim_ref/2     ,1+dim_ref/2))
      allocate(invrsmat_(1+dim_ref/2     ,1+dim_ref/2))
      allocate(matrix2_ (1+dim_ref/2+kmax,1+dim_ref/2))
      allocate(mat_rhs2_(1+dim_ref/2+kmax))
      allocate(mat_rhs_ (1+dim_ref/2     ))
        if(par%rank==0) then
         write(6,*)'-------------------------------------'
         write(6,*)'CALCUL REFERENCE TEMPERATURE'
         write(6,*)'-------------------------------------'
        endif

!      do scz8_=25.,95.,10.
!      do scz6_=40.,180.,20.
!      do scz4_=155.,395.,20.
!      do scz2_=300.,800.,40.
!      trefco_t(:,:,8)=1./scz8_
!      trefco_t(:,:,6)=1./scz6_
!      trefco_t(:,:,4)=1./scz4_
!      trefco_t(:,:,2)=1./scz2_

! TEMPERATURE:
        do j=1,jmax ; do i=1,imax
        if(mask_t(i,j,kmax)==1) then !>>>>>>>>

! Equations supplementaires de stabilite:
          do ligne_=1,1+dim_ref/2
            do colonne_=1,1+dim_ref/2
             matrix2_(kmax+ligne_,colonne_)=0.
            enddo
            matrix2_(kmax+ligne_,ligne_)=1.*0.0001 ! 1 x poids
           mat_rhs2_(kmax+ligne_)=0.
          enddo

! Equations du systeme x0+x1*exp1+x2*exp2....=T
         do ligne_=1,kmax
         k=ligne_

          mat_rhs2_(ligne_)=                &
!                           tem_t(i,j,k,1)  &
                         temref_t(i,j,k)    &
                          *dsig_t(i,j,k)   ! poids

         colonne_=1 ! la colonne 1 est la colonne de la constante
         matrix2_(ligne_,colonne_)=dsig_t(i,j,k)   ! 1 x poids

         do colonne_=2,1+dim_ref/2 ! colonnes 2,3,...

            k1=colonne_*2-2 ! k1=2,4,6,etc.. les numeros des echelles de prof

            matrix2_(ligne_,colonne_)=   dsig_t(i,j,k)    & ! poids
                                      /trefco_t(i,j,k1)   &
               *( exp(depth_w(i,j,k+1)*trefco_t(i,j,k1))  &
                 -exp(depth_w(i,j,k  )*trefco_t(i,j,k1))) &
               /(     depth_w(i,j,k+1)                    &
                     -depth_w(i,j,k  )       )

          enddo ! colonne_

!         write(6,*)real(matrix2_(ligne_,:))


         enddo  ! ligne_

         matrix_ =matmul(transpose(matrix2_),matrix2_)
         mat_rhs_=matmul(transpose(matrix2_),mat_rhs2_)

         lb2=lbound(matrix_) ; ub2=ubound(matrix_)
         call s_matrix_reverser(matrix_,invrsmat_,k0,lb2,ub2)

! On met temporairement la solution dans le tableaux mat_rhs_
         mat_rhs_=matmul(invrsmat_,mat_rhs_)
         trefco_t(i,j,0)=mat_rhs_(1)
         do colonne_=2,1+dim_ref/2
          k1=colonne_*2-3 ! k1=1,3,5,etc.. les numeros des amplitudes
          trefco_t(i,j,k1)=mat_rhs_(colonne_)
         enddo

        endif                        !>>>>>>>>
        enddo         ; enddo

       call pressure_gradient_tref_errors
!      enddo !scz2_
!      enddo !scz4_
!      enddo !scz6_
!      enddo !scz8_

        if(par%rank==0) then
         write(6,*)'-------------------------------------'
         write(6,*)'CALCUL REFERENCE SALINITE'
         write(6,*)'-------------------------------------'
        endif

!      do scz8_=25.,95.,10.
!      do scz6_=40.,180.,20.
!      do scz4_=155.,395.,20.
!      do scz2_=300.,800.,40.
!      srefco_t(:,:,8)=1./scz8_
!      srefco_t(:,:,6)=1./scz6_
!      srefco_t(:,:,4)=1./scz4_
!      srefco_t(:,:,2)=1./scz2_

! SALINITY:
        do j=1,jmax ; do i=1,imax
        if(mask_t(i,j,kmax)==1) then !>>>>>>>>

! Equations supplementaires de stabilite:
          do ligne_=1,1+dim_ref/2
            do colonne_=1,1+dim_ref/2
             matrix2_(kmax+ligne_,colonne_)=0.
            enddo
            matrix2_(kmax+ligne_,ligne_)=1.*0.0001 ! 1 x poids
           mat_rhs2_(kmax+ligne_)=0.
          enddo

! Equations du syssale x0+x1*exp1+x2*exp2....=T
         do ligne_=1,kmax
         k=ligne_

          mat_rhs2_(ligne_)=                  &
!                             sal_t(i,j,k,1)  &
                           salref_t(i,j,k)    &
                            *dsig_t(i,j,k)   ! poids

         colonne_=1 ! la colonne 1 est la colonne de la constante
         matrix2_(ligne_,colonne_)=dsig_t(i,j,k)   ! 1 x poids

         do colonne_=2,1+dim_ref/2 ! colonnes 2,3,...

            k1=colonne_*2-2 ! k1=2,4,6,etc.. les numeros des echelles de prof

            matrix2_(ligne_,colonne_)=   dsig_t(i,j,k)    & ! poids
                                      /srefco_t(i,j,k1)   &
               *( exp(depth_w(i,j,k+1)*srefco_t(i,j,k1))  &
                 -exp(depth_w(i,j,k  )*srefco_t(i,j,k1))) &
               /(     depth_w(i,j,k+1)                    &
                     -depth_w(i,j,k  )       )

          enddo ! colonne_

!         write(6,*)real(matrix2_(ligne_,:))


         enddo  ! ligne_

         matrix_ =matmul(transpose(matrix2_),matrix2_)
         mat_rhs_=matmul(transpose(matrix2_),mat_rhs2_)

         if(matdim/=1+dim_ref/2) then !>>>>
           write(6,*)'matdim=               ',matdim
           write(6,*)'1+dim_ref/2, dim_ref/2',1+dim_ref/2,dim_ref
           stop 'Err matdim=/dim_ref'
         endif                        !>>>>

         lb2=lbound(matrix_) ; ub2=ubound(matrix_)
         call s_matrix_reverser(matrix_,invrsmat_,k0,lb2,ub2)

! On met temporairement la solution dans le tableaux mat_rhs_
         mat_rhs_=matmul(invrsmat_,mat_rhs_)
         srefco_t(i,j,0)=mat_rhs_(1)
         do colonne_=2,1+dim_ref/2
          k1=colonne_*2-3 ! k1=1,3,5,etc.. les numeros des amplitudes
          srefco_t(i,j,k1)=mat_rhs_(colonne_)
         enddo

        endif                        !>>>>>>>>
        enddo         ; enddo

       call pressure_gradient_sref_errors
!      enddo !scz2_
!      enddo !scz4_
!      enddo !scz6_
!      enddo !scz8_

      deallocate(matrix_)
      deallocate(matrix2_)
      deallocate(invrsmat_)
      deallocate(mat_rhs2_)
      deallocate(mat_rhs_)

      call pressure_gradient_ts_ref

      end subroutine pressure_gradient_coref
#endif
!.....................................................................
#ifdef bidon
      subroutine pressure_gradient_tref_errors
      use module_principal ; use module_parallele 
      implicit none
      integer loop_

      id_tem2=2
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_tem2)=trefco_t(i,j,0)
      enddo ; enddo ; enddo
      do loop_=2,dim_ref,2
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

          anyv3d(i,j,k,id_tem2)=                                    &
          anyv3d(i,j,k,id_tem2)+                                    &
       +trefco_t(i,j,loop_-1)/trefco_t(i,j,loop_)                   &
                      *( exp(depth_w(i,j,k+1)*trefco_t(i,j,loop_))  &
                        -exp(depth_w(i,j,k  )*trefco_t(i,j,loop_))) &
                      /(     depth_w(i,j,k+1)                       &
                            -depth_w(i,j,k  )       )
       enddo ; enddo ; enddo
      enddo 

! Ecart a first guess temref:
       sum0=small1  ; sum1=0.
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        sum0=sum0+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)
        sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*abs(anyv3d(i,j,k,id_tem2)-temref_t(i,j,k))
!       if(i==167-par%timax(1).and.j==110-par%tjmax(1))write(67,*)depth_t(i,j,k),anyv3d(i,j,k,id_tem2),temref_t(i,j,k)
       enddo ; enddo ; enddo
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,  & !#MPI 01-06-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)

      if(par%rank==0) then !>>>>
        write(6,*)'Erreur sur temref first guess=',sum1glb/sum0glb
        open(unit=3,file='tmp/temref_errors',access='append')
        write(3,'(10(e14.7,1x))')(1./trefco_t(1,1,k),k=2,dim_ref,2),sum1glb/sum0glb
        close(3)
      endif                !>>>>


      end subroutine pressure_gradient_tref_errors
#endif
!.....................................................................
#ifdef bidon
      subroutine pressure_gradient_sref_errors
      use module_principal ; use module_parallele 
      implicit none
      integer loop_

      id_sal2=2
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_sal2)=srefco_t(i,j,0)
      enddo ; enddo ; enddo
      do loop_=2,dim_ref,2
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

          anyv3d(i,j,k,id_sal2)=                                    &
          anyv3d(i,j,k,id_sal2)+                                    &
       +srefco_t(i,j,loop_-1)/srefco_t(i,j,loop_)                   &
                      *( exp(depth_w(i,j,k+1)*srefco_t(i,j,loop_))  &
                        -exp(depth_w(i,j,k  )*srefco_t(i,j,loop_))) &
                      /(     depth_w(i,j,k+1)                       &
                            -depth_w(i,j,k  )       )
       enddo ; enddo ; enddo
      enddo 

! Ecart a first guess salref:
       sum0=small1  ; sum1=0.
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        sum0=sum0+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)
        sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*abs(anyv3d(i,j,k,id_sal2)-salref_t(i,j,k))
!       if(i==167-par%timax(1).and.j==110-par%tjmax(1))write(68,*)depth_t(i,j,k),anyv3d(i,j,k,id_sal2),salref_t(i,j,k)
       enddo ; enddo ; enddo
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,  & !#MPI 01-06-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)

      if(par%rank==0) then !>>>>
        write(6,*)'Erreur sur salref first guess=',sum1glb/sum0glb
        open(unit=3,file='tmp/salref_errors',access='append')
        write(3,'(10(e14.7,1x))')(1./srefco_t(1,1,k),k=2,dim_ref,2),sum1glb/sum0glb
        close(3)
      endif                !>>>>


      end subroutine pressure_gradient_sref_errors
#endif
!.....................................................................
#ifdef bidonref
      subroutine pressure_gradient_ts_ref
      use module_principal ; use module_parallele 
      implicit none
      integer loop_

      id_tem=0 ; id_sal=1
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_tem)=temref_t(i,j,k)
        anyv3d(i,j,k,id_sal)=salref_t(i,j,k)
       enddo ; enddo ; enddo

! OBC
      do loop_=0,1  ! O=tem ; 1=sal
       do k=1,kmax
       do i=1,imax
        anyv3d(i,jmax+2,k,loop_)=anyv3d(i,jmax,k,loop_)
        anyv3d(i,jmax+1,k,loop_)=anyv3d(i,jmax,k,loop_)
        anyv3d(i,    -1,k,loop_)=anyv3d(i,1   ,k,loop_)
        anyv3d(i,     0,k,loop_)=anyv3d(i,1   ,k,loop_)
       enddo
       enddo
       do k=1,kmax
       do j=-1,jmax+2
        anyv3d(imax+2,j,k,loop_)=anyv3d(imax,j,k,loop_)
        anyv3d(imax+1,j,k,loop_)=anyv3d(imax,j,k,loop_)
        anyv3d(    -1,j,k,loop_)=anyv3d(1   ,j,k,loop_)
        anyv3d(     0,j,k,loop_)=anyv3d(1   ,j,k,loop_)
       enddo
       enddo
       call obc_int_anyv3d(loop_,'zc')
      enddo ! loop_ 
! Pourquoi il faut etendre les boucles de -1 a imax+2 (avec obc mpi qui va avec) ?
! A cause de l'advection up3 et des routines de soustraction/addition de temref_t salref_t
      id_tem=0 ; id_sal=1
      do k=1,kmax
      do j=-1,jmax+2
      do i=-1,imax+2
        temref_t(i,j,k)=anyv3d(i,j,k,id_tem)
        salref_t(i,j,k)=anyv3d(i,j,k,id_sal)
      enddo
      enddo
      enddo

      end subroutine pressure_gradient_ts_ref
#endif
!.........................................................................
#ifdef bidon
      subroutine pressure_gradient_rhpref_errors
      use module_principal
      use module_parallele
      implicit none

! A partir de temref et salref on calcule la densite potentielle que l'on compare
! rhpref qui contient encore a ce stade le first guess de la densite potentielle
! de reference
      id_tem=0 ; id_sal=1
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_tem)=temref_t(i,j,k)
        anyv3d(i,j,k,id_sal)=salref_t(i,j,k)
       enddo ; enddo ; enddo
! DENSITE POTENTIELLE DE REFERENCE ASSOCIEE
       id_rhp=2
       call equation_of_state_pot_anyv(1,imax,1,jmax)

       sum0=small1  ; sum1=0.
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        sum0=sum0+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)
        sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*abs(anyv3d(i,j,k,id_rhp)-rhpref_t(i,j,k))
!       if(i==167-par%timax(1).and.j==110-par%tjmax(1))write(69,*)depth_t(i,j,k),anyv3d(i,j,k,id_rhp),rhpref_t(i,j,k)
       enddo ; enddo ; enddo
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,  & !#MPI 01-06-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)

      if(par%rank==0)write(6,*)'Erreur sur rhpref first guess=',sum1glb/sum0glb

      end subroutine pressure_gradient_rhpref_errors
#endif
!.........................................................................
#ifdef bidonref
      subroutine pressure_gradient_homogeneousref(txt_,flag_) !14-11-15
      use module_principal
      use module_parallele
      implicit none
      integer flag_,loop_
      character*1 txt_

      if(txt_=='u') then !uuuuuuuuuuuu>


       do loop_=-1,0

       if(loop_==-1)call quick_tsref_interp(1) ! arg 1 = interpolation de la moyenne i-1 , i sur i1=i-1
       if(loop_== 0)call quick_tsref_interp(2) ! arg 2 = interpolation de la moyenne i-1 , i sur i1=i

       if(loop_==-1) then !mmmmmmm>
         if(flag_==0) then !>>>>
          id_rhp=id_rhom
          call equation_of_state_pot_anyv (1,imax-1,2,jmax-1) ! Arg tiennent compte de i-1
         endif             !>>>>
         if(flag_==1) then !>>>>
          id_rhf=id_rhom
          call equation_of_state_full_anyv(1,imax-1,2,jmax-1) ! Arg tiennent compte de i-1
         endif             !>>>>
       endif              !mmmmmmm>

       if(loop_==0)  then !ppppppp>
         if(flag_==0) then !>>>>
          id_rhp=id_rhop
          call equation_of_state_pot_anyv (2,imax,2,jmax-1)
         endif             !>>>>
         if(flag_==1) then !>>>>
          id_rhf=id_rhop
          call equation_of_state_full_anyv(2,imax,2,jmax-1)
         endif             !>>>>
       endif              !ppppppp>


!      if(loop_==-1) then
!       i=22 ; j=22
!       do k=kmax,1,-1
!       write(6,*)k,anyv3d(i,j,k,id_tem),temref_t(i,j,k)
!       write(6,*)k,anyv3d(i,j,k,id_sal),salref_t(i,j,k)
!       enddo
!      stop 'cocote_u'
!      endif

       enddo  ! loop_

!...................
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!       stop 'mimi'


      endif              !uuuuuuuuuuuu>

      if(txt_=='v') then !vvvvvvvvvvvv>

       do loop_=-1,0

       if(loop_==-1)call quick_tsref_interp(3) ! arg 3 = interpolation de la moyenne j-1 , j sur j1=j-1
       if(loop_== 0)call quick_tsref_interp(4) ! arg 4 = interpolation de la moyenne j-1 , j sur j1=j

       if(loop_==-1) then !mmmmmmm>
         if(flag_==0) then !>>>>
          id_rhp=id_rhom
          call equation_of_state_pot_anyv (2,imax-1,1,jmax-1) ! Arg tiennent compte de j-1
         endif             !>>>>
         if(flag_==1) then !>>>>
          id_rhf=id_rhom
          call equation_of_state_full_anyv(2,imax-1,1,jmax-1) ! Arg tiennent compte de j-1
         endif             !>>>>
       endif              !mmmmmmm>

       if(loop_==0)  then !ppppppp>
         if(flag_==0) then !>>>>
          id_rhp=id_rhop
          call equation_of_state_pot_anyv (2,imax-1,2,jmax)
         endif             !>>>>
         if(flag_==1) then !>>>>
          id_rhf=id_rhop
          call equation_of_state_full_anyv(2,imax-1,2,jmax)
         endif             !>>>>
       endif              !ppppppp>

!      if(loop_==0) then
!       i=22 ; j=22
!       do k=kmax,1,-1
!       write(6,*)k,anyv3d(i,j,k,id_tem),temref_t(i,j,k)
!       write(6,*)k,anyv3d(i,j,k,id_sal),salref_t(i,j,k)
!       enddo
!      stop 'cocote_v'
!      endif

       enddo  ! loop_

!...................
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!       stop 'mimi v'


      endif              !vvvvvvvvvvvv>

      end subroutine pressure_gradient_homogeneousref
#endif
!....................................................................
      subroutine pressure_gradient_1dvcase
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_1dvcase'
       subroutinedescription=''
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       if(iobc_ogcm==0)return

! 1DV case:
! Pressure gradient consistent with geostrophic current
        do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
         presgrad_u(i,j,k,1)= coriolis_t(i,j)*dx_u(i,j)*(                 &
                               timeweightobc(vel_id) *velobc_v(i,j,k,2) &!30-11-15
                          +(1.-timeweightobc(vel_id))*velobc_v(i,j,k,0))

        enddo       ; enddo         ; enddo
        do k=1,kmax ; do j=2,jmax   ; do i=2,imax-1
         presgrad_v(i,j,k,1)=-coriolis_t(i,j)*dy_v(i,j)*(                 &
                               timeweightobc(vel_id) *velobc_u(i,j,k,2) &
                          +(1.-timeweightobc(vel_id))*velobc_u(i,j,k,0))
        enddo       ; enddo         ; enddo

      end subroutine pressure_gradient_1dvcase

!..............................................................................

      subroutine pressure_gradient_eos0
      use module_principal ; use module_parallele ; use module_modeanalysis
      implicit none
#ifdef synopsis
       subroutinetitle='pressure_gradient_eos0'
       subroutinedescription= &
       'Computes the pressure gradient force according to the method'  &
       //' described in: Marsaleix P., Auclair F., Estournel C., 2009' &
       //' Low-order pressure gradient schemes in sigma coordinate'    &
       //' models: The seamount test revisited. Ocean Modelling, 30,'  &
       //' 169-177. http://dx.doi.org/10.1016/j.ocemod.2009.06.011'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!-----------------------------------
! Pour arbirtrairement eliminer la densite potentielle
! et calculer le PGF OM2009 avec la densite de compression:
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     write(6,*)'BIDOUILLE'
!     call equation_of_state_pressure_jmfwg(now)
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      rhp_t(i,j,k)=anyv3d(i,j,k,0)
!     enddo ; enddo ; enddo
!-----------------------------------

! Details sur rho qui depend de (x,y) dans:
!https://docs.google.com/document/d/1_SdjMHFzI2EVxk17kXwjsXi4ra1nHiU_Wyg2LyL9KZ8/edit?usp=sharing

! Compute the potential density:
!     if(iteration3d==0) then !>>>
!      if(flag_refstate==1)call pressure_gradient_initref !15-10-15
!     endif                   !>>>

! Compute rhp_t(:,:,:)
      if(eos_comprs==0) then !000000>
! Compute rhp the potential density at z=0m and the a1, a2, ....,a7 coefficients of the Jackett et al 2006 (JAOT) EOS
       if(eos_author==0)call equation_of_state_linear(now)
       if(eos_author==3)call equation_of_state_rhp_a1_a7(now) ! rhp=densite potentielle (zref=0)
       if(eos_author==2)stop 'Err 697 pressure_gradient'
       if(eos_author==1)stop 'Err 698 pressure_gradient'
      endif                  !000000>

      if(eos_comprs==1) then !111111>
       if(eos_author==0)stop 'Err 699 pressure_gradient'
       if(eos_author==3)call equation_of_state_full_a1_a7(now) ! rhp=densite potentielle (zref=0)
       if(eos_author==2)stop 'Err 680 pressure_gradient'
       if(eos_author==1)stop 'Err 681 pressure_gradient'
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhp_t(i,j,k)=rhp_t(i,j,k)-rhcref_t(i,j,k) ! Ote l'etat de reference de compression
       enddo       ; enddo       ; enddo
      endif                  !111111>

! En cas d'analyse harmonique de la pression A partir de rhp avant que rhmean ne soit otE:
      if(flag_3dwaves_pharmonics==1)call modeanalysis_pmodeprojection !14-04-16

!     write(6,*)'BIDOUILLE PGF'
!     x1=1./real(imax+jmax)
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      rhp_t(i,j,k)=real(i+j)*dxb*x1
!     enddo
!     enddo
!     enddo

! remove a depth-weighted average of rhp
      if(rhp_zavr_xy==1)call pressure_gradient_rhpzavr

!$ Computes the Hydrostatic pressure using the rectangular method:
      id_prs=0 ! attention anyv3d(:,:,:,1 a 4 et de 6 a 7) est pris par les coef de l'EOS

! verifier la disponibilite de anyv3d(:,:,:,id_prs)
      if(anyv3d(-1,-1,0,id_prs)==-9999.) &
      stop 'Err anyv3d id_prs not available'

      if(flag_1dv==0) then !3D3D3D3D3D> !16-02-16

      do 9 j=1,jmax
      do 9 i=1,imax
      xy_t(i,j,1)=0.
    9 continue

      do 10 k=kmax,1,-1
       do 10 j=1,jmax
       do 10 i=1,imax

!$ Hydrostatic pressure:
       x2=grav*rhp_t(i,j,k)*dz_t(i,j,k,1)              !02-10-09
       xy_t(i,j,1)=xy_t(i,j,1)+x2                      !06/03/03
       anyv3d(i,j,k,id_prs)=xy_t(i,j,1)-0.5*x2         ! p(i,j,k)

   10 continue

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      const1=1./(  rho)*flag3d                                           !30-09-09
      const2=grav/(rho)*flag3d                                           !30-09-09

      do 14 k=kmax,1,-1
      do 14 j=2,jmax-1
      do 14 i=2,imax

!$ dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=                                               &
            +const1*( anyv3d(i+ip,j,k,id_prs)-anyv3d(i-im,j,k,id_prs)) &
            +const2*(depth_t(i+ip,j,k)      -depth_t(i-im,j,k))*(      &

          (  dz_t(i+ip,j,k,1)*(rhp_t(i+ip,j,k)       &
#ifdef bidonref
                           -rhpref_t(i+ip,j,k)       &
#endif
                          +rhpzavr_w(i+ip,j)   )     &!16-09-15
            +dz_t(i-im,j,k,1)*(rhp_t(i-im,j,k)       &
#ifdef bidonref
                           -rhpref_t(i-im,j,k)       &
#endif
                          +rhpzavr_w(i-im,j))        &
!         (  dz_t(i+ip,j,k,1)*(rhp_t(i+ip,j,k)-rhpref_t(i+ip,j,k))     &
!           +dz_t(i-im,j,k,1)*(rhp_t(i-im,j,k)-rhpref_t(i-im,j,k))     &
          )/(dz_t(i+ip,j,k,1)+dz_t(i-im,j,k,1))                        &

#ifdef bidonref
          +rhpref_u(i,j,k)            &
#endif
                             )        &

! rhzavr (if any) contribution:
!          -const2*depth_u(i,j,k)*( rhpzavr_w(i+ipu,j)                 &
!                                  -rhpzavr_w(i-imu,j))                &


           -const2*( rhpzavr_w(i+ip,j)*depth_t(i+ip,j,k)              & !16-09-15
                    -rhpzavr_w(i-im,j)*depth_t(i-im,j,k))   

!       if(i+par%timax(1)==260.and.j+par%tjmax(1)==40) then
!        write(52,*)depth_u(i,j,k),presgrad_u(i,j,k,1)
!       endif

   14 continue


!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      const1=1./(   rho)*flag3d                                            !30-09-09
      const2= grav/(rho)*flag3d                                            !30-09-09

      do 19 k=kmax,1,-1
      do 19 j=2,jmax
      do 19 i=2,imax-1

!$ dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=                                               &

            +const1*( anyv3d(i,j+jp,k,id_prs)-anyv3d(i,j-jm,k,id_prs)) &
            +const2*(depth_t(i,j+jp,k)      -depth_t(i,j-jm,k))*(      &

          (  dz_t(i,j+jp,k,1)*(rhp_t(i,j+jp,k)     &
#ifdef bidonref
                           -rhpref_t(i,j+jp,k)     &
#endif
                          +rhpzavr_w(i,j+jp)     ) & !16-09-15
            +dz_t(i,j-jm,k,1)*(rhp_t(i,j-jm,k)     &
#ifdef bidonref
                           -rhpref_t(i,j-jm,k)     &
#endif
                          +rhpzavr_w(i,j-jm)     ) &
          )/(dz_t(i,j+jp,k,1)+dz_t(i,j-jm,k,1))                        &

#ifdef bidonref
          +rhpref_v(i,j,k)       &
#endif
                             )                                         &

! rhzavr (if any) contribution:
         -const2*( rhpzavr_w(i,j+jpv)*depth_t(i,j+jpv,k)              & !16-09-15
                  -rhpzavr_w(i,j-jmv)*depth_t(i,j-jmv,k))

   19 continue

      else                 !3D3D3D3D3D> or !1D1D1D1D1D> !16-02-16

       call pressure_gradient_1dvcase !29-11-15

      endif                                !1D1D1D1D1D>

!....................................
! The TKE scheme is coming next
! rhp_t possibly needs to be redefined if the level of reference of the potential
! density in the TKE scheme is different from that used in the PGF:
!     if(eos_pgfzref/=eos_tkezref) then !rhprhprhp> !20-06-14

!       if(eos_author/=0) then
! Systematiquement calculer la densite potentielle maintenant que celle-ci est egalement !24-04-16
! utilisee dans le calcul de omega_w(:,:,kmax+1) et le calcul de rho.cp dans vertmix.F90
       if(eos_author==0) then !>>>>> !25-04-16
           call equation_of_state_linear(now)
       else                   !>>>>> !25-04-16
        if(eos_tkezref>=0.) then !-------->
          call equation_of_state_potzref2                 ! z reference = eos_tkezref !20-06-14
        else                     !-------->
          call equation_of_state_potloc2                  ! reference is current depth
        endif                    !-------->
       endif                  !>>>>> !25-04-16

! Marquer "disponibles" les tableaux generiques 
      anyv3d(-1,-1,0,1:4)=0. ; anyv3d(-1,-1,0,6:7)=0.
      anyv3d(-1,-1,0,id_prs)=0.


!     endif                             !rhprhprhp>
!....................................
      end subroutine pressure_gradient_eos0

!.....................................................................................

      subroutine pressure_gradient_add_nhpgf
      use module_principal ; use module_parallele ; use module_q
      implicit none

! compute nhpgf_u nhpgf_v
      call q_nhpgf3d

! Cas du NH sans mode splitting et sans T et S
      if(flag_nh3d==flag_nh3d_nosplit_uv) then !m°v°m> !23-08-18

       if(timestep_type==timestep_forwbckw) then !fbfbfb> 
! Bascule puisque le pgf n'a jamais ete calculE si flag_nh3d=1
        presgrad_u(:,:,:,0)=presgrad_u(:,:,:,1)
        presgrad_v(:,:,:,0)=presgrad_v(:,:,:,1)
       endif                                     !fbfbfb>


       inv_rho=1./rho
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        presgrad_u(i,j,k,1)=grav*(ssh_int_w(i,j,1)-ssh_int_w(i-1,j,1))  & !20-08-18
                        +inv_rho*(    pss_w(i,j,1)    -pss_w(i-1,j,1))  & !08-10-20
          +nhpgf_u(i,j,k)
       enddo ; enddo ; enddo

       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        presgrad_v(i,j,k,1)=grav*(ssh_int_w(i,j,1)-ssh_int_w(i,j-1,1)) & !20-08-18
                        +inv_rho*(    pss_w(i,j,1)    -pss_w(i,j-1,1)) & !08-10-20
          +nhpgf_v(i,j,k)
       enddo ; enddo ; enddo

      endif                                    !m°v°m>

! Cas du NH sans mode splitting AVEC T et S
      if(flag_nh3d==flag_nh3d_nosplit_tsuv) then !m°-°m> !23-08-18

! Pas de bascule pgf(0)=pgf(1) puisque deja calculE
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        presgrad_u(i,j,k,1)= &
        presgrad_u(i,j,k,1)+grav*(ssh_int_w(i,j,1)-ssh_int_w(i-1,j,1)) & !20-08-18
          +nhpgf_u(i,j,k)
       enddo ; enddo ; enddo

       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        presgrad_v(i,j,k,1)= &
        presgrad_v(i,j,k,1)+grav*(ssh_int_w(i,j,1)-ssh_int_w(i,j-1,1)) & !20-08-18
          +nhpgf_v(i,j,k)
       enddo ; enddo ; enddo

      endif                                      !m°-°m> !23-08-18

! Cas du NH AVEC time splitting AVEC T et S
      if(flag_nh3d==flag_nh3d_timesplit_tsuv) then !w°o°w> !06-09-18

! Pas de bascule pgf(0)=pgf(1) puisque deja calculE
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        presgrad_u(i,j,k,1)= &
        presgrad_u(i,j,k,1)  & !20-08-18
          +nhpgf_u(i,j,k)
       enddo ; enddo ; enddo

       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        presgrad_v(i,j,k,1)= &
        presgrad_v(i,j,k,1)  & !20-08-18
          +nhpgf_v(i,j,k)
       enddo ; enddo ; enddo

      endif                                        !w°o°w> !06-09-18

      end subroutine pressure_gradient_add_nhpgf

!.....................................................................................

      subroutine pressure_gradient_vqs !09-06-20
      use module_principal ; use module_parallele ; use module_modeanalysis
      implicit none

! Explications concernant la prise en compte de la densite moyenne 2D dans le calcul du PGF:
! https://docs.google.com/document/d/1_SdjMHFzI2EVxk17kXwjsXi4ra1nHiU_Wyg2LyL9KZ8/edit?usp=sharing

! Explications concernant le cas particulier du calcul de la couche kmin
! https://docs.google.com/presentation/d/1PohrXJbiEMyCNGmGO4lh9n0OQyRR8apDQm-A9o5sLZk/edit#slide=id.g898a38edfc_0_25

      if(flag_refstate==1) &
      stop 'Err flag_refstate=1 dans pressure_gradient_vqs'


! Compute rhp_t(:,:,:)
      if(eos_comprs==0) then !000000>
! Compute rhp the potential density at z=0m and the a1, a2, ....,a7 coefficients of the Jackett et al 2006 (JAOT) EOS
       if(eos_author==0)call equation_of_state_linear(now)
       if(eos_author==3)call equation_of_state_rhp_a1_a7(now) ! rhp=densite potentielle (zref=0)
       if(eos_author==2)stop 'Err 697 pressure_gradient'
       if(eos_author==1)stop 'Err 698 pressure_gradient'
       if(flag_tide3d_analysis==1 &
             .or.drifter_onoff==1) then !---> !07-03-23

        rho_t(0:imax+1,0:jmax+1,1:kmax)=rhp_t(0:imax+1,0:jmax+1,1:kmax)+rho !26-01-23
       endif                            !--->
      endif                  !000000>

      if(eos_comprs==1) then !111111>
       if(eos_author==0)stop 'Err 699 pressure_gradient'
       if(eos_author==3)call equation_of_state_full_a1_a7(now) ! rhp=densite potentielle (zref=0)
       if(eos_author==2)stop 'Err 680 pressure_gradient'
       if(eos_author==1)stop 'Err 681 pressure_gradient'
       if(flag_tide3d_analysis==1 &
             .or.drifter_onoff==1) then !---> !07-03-23
        rho_t(0:imax+1,0:jmax+1,1:kmax)=rhp_t(0:imax+1,0:jmax+1,1:kmax)+rho !26-01-23
       endif                            !--->
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        rhp_t(i,j,k)=rhp_t(i,j,k)-rhcref_t(i,j,k) ! Ote l'etat de reference de compression
       enddo       ; enddo       ; enddo
      endif                  !111111>

!**********************************
! Profil academique
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      rhp_t(i,j,k)=-0.001*depth_t(i,j,k)
!     enddo
!     enddo
!     enddo
! melange couche fusionnee
!     do j=1,jmax ; do i=1,imax
!        sum1=0. ; sum2=0. 
!        do k=1,kmerged_t(i,j)
!         sum1=sum1+dz_t(i,j,k,1)
!         sum2=sum2+dz_t(i,j,k,1)*rhp_t(i,j,k)
!        enddo
!        x1=sum2/max(sum1,small1)*mask_t(i,j,kmax)
!        do k=1,kmerged_t(i,j)
!         rhp_t(i,j,k)=x1
!        enddo
!     enddo ; enddo
!**********************************

! Profil lineaire reconstituE dans la couche fusionnE !09-06-20 (supprimeE le 02-11-21)
!     call turbulence_rhp_linear_profil               !09-06-20!02-11-21

! En cas d'analyse harmonique de la pression A partir de rhp avant que rhmean ne soit otE:
      if(flag_3dwaves_pharmonics==1)call modeanalysis_pmodeprojection !14-04-16

! remove a depth-weighted average of rhp
      if(rhp_zavr_xy==1)call pressure_gradient_rhpzavr

!$ Computes the Hydrostatic pressure using the rectangular method:
      id_prs=0 ! attention anyv3d(:,:,:,1 a 4 et de 6 a 7) est pris par les coef de l'EOS

! verifier la disponibilite de anyv3d(:,:,:,id_prs)
      if(anyv3d(-1,-1,0,id_prs)==-9999.) &
      stop 'Err anyv3d id_prs not available'

      if(flag_1dv==0) then !3D3D3D3D3D> !16-02-16

      do j=1,jmax ; do i=1,imax
       anyv3d(i,j,kmax+1,id_prs)=0.
      enddo       ; enddo
      do k=kmax,1,-1 ; do j=1,jmax ; do i=1,imax

             anyv3d(i,j,k  ,id_prs)= &
             anyv3d(i,j,k+1,id_prs)  &
        +grav*rhp_t(i,j,k)           &
              *dz_t(i,j,k,1)

      enddo ; enddo ; enddo
! Refait le profil de pression dans la couche fusionnEe et en deduit rhp_t en k=kmin_w et kmerged_t: 
!     call pressuregradient_remap_prs !29-10-21

!$ Computes the Ox component of the pressure gradient force multiplied by dx_u:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      ip=ipu ; im=imu ! ip=0 im=1                                     !04-03-10

      inv_rho=1./rho
      do k=kmax,1,-1 ; do j=2,jmax-1 ; do i=2,imax

!dx_u(i,j) * dp/dx:
      presgrad_u(i,j,k,1)=  &
          mask_u(i,j,k)*inv_rho*( & !PMXPMX>

          +0.5*( anyv3d(i+ipu,j,k+1,id_prs)-anyv3d(i-imu,j,k+1,id_prs)  &
                +anyv3d(i+ipu,j,k  ,id_prs)-anyv3d(i-imu,j,k  ,id_prs)) &

      +(depth_t(i+ipu,j,k)-depth_t(i-imu,j,k))*( & !ooo>

       (  anyv3d(i+ipu,j,k,id_prs)-anyv3d(i+ipu,j,k+1,id_prs)  &
         +anyv3d(i-imu,j,k,id_prs)-anyv3d(i-imu,j,k+1,id_prs)  &
       )/(dz_t(i+ipu,j,k,1)+dz_t(i-imu,j,k,1))                         &

                                             ) & !ooo>

! rhzavr (if any) contribution:
           -grav*depth_u(i,j,k)*( rhpzavr_w(i+ipu,j)                       &
                                 -rhpzavr_w(i-imu,j))  ) !PMXPMX>

      enddo ; enddo ; enddo

!$ Computes the Oy component of the pressure gradient force multiplied by dy_v:

!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
      jp=jpv ; jm=jmv  ! jp=0 jm=1                                      !04-03-10

      inv_rho=1./rho
      do k=kmax,1,-1 ; do j=2,jmax ; do i=2,imax-1

! dy_v(i,j) * dp/dy:
      presgrad_v(i,j,k,1)=  &
          mask_v(i,j,k)*inv_rho*( & !PMXPMX>

          +0.5*( anyv3d(i,j+jpv,k+1,id_prs)-anyv3d(i,j-jmv,k+1,id_prs)  &
                +anyv3d(i,j+jpv,k  ,id_prs)-anyv3d(i,j-jmv,k  ,id_prs)) &

      +(depth_t(i,j+jpv,k)-depth_t(i,j-jmv,k))*( & !ooo>
           
        (   anyv3d(i,j+jpv,k,id_prs)-anyv3d(i,j+jpv,k+1,id_prs) &
           +anyv3d(i,j-jmv,k,id_prs)-anyv3d(i,j-jmv,k+1,id_prs) &
           )/(dz_t(i,j+jpv,k,1)+dz_t(i,j-jmv,k,1))                      &

                                               ) & !ooo>

! rhzavr (if any) contribution:
           -grav* depth_v(i,j,k)*( rhpzavr_w(i,j+jpv)                       &
                                  -rhpzavr_w(i,j-jmv))  ) !PMXPMX>

      enddo ; enddo ; enddo

      else                 !3D3D3D3D3D> or !1D1D1D1D1D> !16-02-16

       call pressure_gradient_1dvcase !29-11-15

      endif                                !1D1D1D1D1D>


! COUCHES MERGED:
      k1=1 ; k2=0 ; k3=1 ; k4=0 ; rap=1.
      do j=2,jmax-1 ; do i=2,imax

!     if(kmerged_u(i,j)/=kmin_u(i,j)) then !m°v°m> !option deux couches regroupees  !29-10-21
! note parce que l'option 2 couches regroupees peut conduire A construire le PGF sur !29-20-21
! kmin_u-1:kmin_u si jamais kmerged_u=kmin_u on conditionne le calcul A la condition que kmerged_u/=kmin_u

      k=kmerged_u(i,j)-1
      k4=k+1 ! option deux couches distinctes
      k2=k+1
!     k4=k+2 ! option deux couches regroupees (kmerged-1 et kmerged) !29-10-21
!     k2=k+2

      if(dsig_t(i,j,k)<dsig_t(i-1,j,k)) then !>>> ! La plus grande couche est diminuee CSP1
!     if(dsig_t(i,j,k)>dsig_t(i-1,j,k)) then !>>> ! La plus petite couche est agrandie CSP2

       k1=0 ; k3=k  

       depth_w(i-1,j,k1)=depth_w(i-1,j,k+1) &
!                          -dz_t(i,j,k,1)             ! dz methode 1
                         -dsig_t(i,j,k)*hz_w(i-1,j,1) ! dz methode 2 !03-11-21

! interpolation ordre 1
!        rhp_t(i-1,j,k1)=  &
!        rhp_t(i-1,j,k+1)-(rhp_t(i-1,j,k+1)  -rhp_t(i-1,j,k)) &
!                      /(depth_t(i-1,j,k+1)-depth_t(i-1,j,k)) &
!                      *(depth_t(i-1,j,k+1)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1)))
! interpolation ordre 2 !29-10-21
      rhp_t(i-1,j,k1)=( &
         (  (  rhp_t(i-1,j,k+1)                                           ) &
           /(depth_t(i-1,j,k+1)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))) &
           -(                                            -rhp_t(i-1,j,k))   &
        /max(0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))-depth_t(i-1,j,k),small1)   &
         ) /(depth_t(i-1,j,k+1)-depth_t(i-1,j,k))                           &
        -(  (  rhp_t(i-1,j,k+2)  -rhp_t(i-1,j,k+1)) &
           /(depth_t(i-1,j,k+2)-depth_t(i-1,j,k+1)) &
           -(  rhp_t(i-1,j,k+1)                   ) &
           /(depth_t(i-1,j,k+1)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))) &
         ) /(depth_t(i-1,j,k+2)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))) &
                      ) &
                     /( &
          1./(depth_t(i-1,j,k+1)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))) &
            /(depth_t(i-1,j,k+2)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))) &
         +1./(depth_t(i-1,j,k+1)-0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))) &
            /(depth_t(i-1,j,k+1)-depth_t(i-1,j,k)) &
         +1./(depth_t(i-1,j,k+1)-depth_t(i-1,j,k)) &
       /max(0.5*(depth_w(i-1,j,k+1)+depth_w(i-1,j,k1))-depth_t(i-1,j,k),small1)  &
                      )

             anyv3d(i-1,j,k1 ,id_prs)=    &
             anyv3d(i-1,j,k+1,id_prs)     &
        +grav*rhp_t(i-1,j,k1)             &
!             *dz_t(i,j,k,1)             ! dz methode 1
            *dsig_t(i,j,k)*hz_w(i-1,j,1) ! dz methode 2 !03-11-21

      else                                   !>>>

       k3=0 ; k1=k  

       depth_w(i,j,k3)=depth_w(i,j,k+1)   &
!                        -dz_t(i-1,j,k,1)           ! dz methode 1
                       -dsig_t(i-1,j,k)*hz_w(i,j,1) ! dz methode 2 !03-11-21

! interpolation ordre 1
!        rhp_t(i,j,k3)=  &
!        rhp_t(i,j,k+1)-(rhp_t(i,j,k+1)  -rhp_t(i,j,k)) &
!                    /(depth_t(i,j,k+1)-depth_t(i,j,k)) &
!                    *(depth_t(i,j,k+1)-0.5*(depth_w(i,j,k+1)+depth_w(i,j,k3)))
! interpolation ordre 2 !29-10-21
      rhp_t(i  ,j,k3)=( &
         (  (  rhp_t(i  ,j,k+1)                                           ) &
           /(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
           -(                                            -rhp_t(i  ,j,k))   &
        /max(0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))-depth_t(i  ,j,k),small1)   &
         ) /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k))                           &
        -(  (  rhp_t(i  ,j,k+2)  -rhp_t(i  ,j,k+1)) &
           /(depth_t(i  ,j,k+2)-depth_t(i  ,j,k+1)) &
           -(  rhp_t(i  ,j,k+1)                   ) &
           /(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
         ) /(depth_t(i  ,j,k+2)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
                      ) &
                     /( &
          1./(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
            /(depth_t(i  ,j,k+2)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
         +1./(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
            /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k)) &
         +1./(depth_t(i  ,j,k+1)-depth_t(i  ,j,k)) &
        /max(0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))-depth_t(i  ,j,k),small1)  &
                      )

             anyv3d(i,j,k3 ,id_prs)=    &
             anyv3d(i,j,k+1,id_prs)     &
        +grav*rhp_t(i,j,k3)             &
!             *dz_t(i-1,j,k,1)           ! dz methode 1
            *dsig_t(i-1,j,k)*hz_w(i,j,1) ! dz methode 2 !03-11-21

      endif                                  !>>>

!dx_u(i,j) * dp/dx:
                               xy_u(i,j,1)=inv_rho*( & !PMXPMX> option deux couches distinctes
!     presgrad_u(i,j,1:kmerged_u(i,j)-1,1)=inv_rho*( & !PMXPMX> option deux couches distinctes
!     presgrad_u(i,j,1:kmerged_u(i,j)  ,1)=inv_rho*( & !PMXPMX> option deux couches regroupees !29-10-21

      +0.5*( anyv3d(i+ipu,j,k4,id_prs)-anyv3d(i-imu,j,k2,id_prs)   & ! delta P(bottom+1)
            +anyv3d(i+ipu,j,k3,id_prs)-anyv3d(i-imu,j,k1,id_prs) ) & ! delta P(bottom)

         +0.5*(depth_w(i+ipu,j,k4)-depth_w(i-imu,j,k2)    &          ! Z(bottom+1)
              +depth_w(i+ipu,j,k3)-depth_w(i-imu,j,k1))*( & !ooo>    ! Z(bottom)

        (  anyv3d(i+ipu,j,k3,id_prs)-anyv3d(i+ipu,j,k4,id_prs) &
          +anyv3d(i-imu,j,k1,id_prs)-anyv3d(i-imu,j,k2,id_prs) &
          )/(depth_w(i+ipu,j,k4)-depth_w(i+ipu,j,k3)           &
            +depth_w(i-imu,j,k2)-depth_w(i-imu,j,k1))          &

                                                      ) & !ooo>


           -grav*0.25*(depth_w(i+ipu,j,k4)     & ! < depth_u(i,j,k) >
                      +depth_w(i-imu,j,k2)     &
                      +depth_w(i+ipu,j,k3)     &
                      +depth_w(i-imu,j,k1))    &
                               *( rhpzavr_w(i+ipu,j)     &
                                 -rhpzavr_w(i-imu,j))  ) & !PMXPMX>
! lignes suivantes pour retenir 100% du calcul si dz(kmerged-1)>=dz(kmerged)
! et prendre la valeur de PGF(kmerged) si dz(kmerged-1)=0: !03-11-21
           *pgfratio_u(i,j)+(1.-pgfratio_u(i,j))*presgrad_u(i,j,kmerged_u(i,j),1)
!          *min(1.,dz_u(i,j,kmerged_u(i,j)-1,1)/dz_u(i,j,kmerged_u(i,j),1))  & 
!      +(1.-min(1.,dz_u(i,j,kmerged_u(i,j)-1,1)/dz_u(i,j,kmerged_u(i,j),1)))*presgrad_u(i,j,kmerged_u(i,j),1)

!     endif                                !m°v°m> !option deux couches regroupees  !29-10-21
      enddo ; enddo

      do j=2,jmax-1 ; do i=2,imax
       do k=1,kmerged_u(i,j)-1 
        presgrad_u(i,j,k,1)=xy_u(i,j,1)
       enddo
      enddo ; enddo

      k1=1 ; k2=0 ; k3=1 ; k4=0 ; rap=1.
      do j=2,jmax ; do i=2,imax-1

!     if(kmerged_v(i,j)/=kmin_v(i,j)) then !m°v°m> !option deux couches regroupees  !29-10-21
! note parce que l'option 2 couches regroupees peut conduire A construire le PGF sur !29-20-21
! kmin_v-1:kmin_v si jamais kmerged_v=kmin_v on conditionne le calcul A la condition que kmerged_v/=kmin_v

      k=kmerged_v(i,j)-1
      k4=k+1 ! option deux couches distinctes
      k2=k+1
!     k4=k+2 ! option deux couches regroupees !29-10-21
!     k2=k+2

      if(dsig_t(i,j,k)<dsig_t(i,j-1,k)) then !>>> La plus grande couche est diminuee CSP1
!     if(dsig_t(i,j,k)>dsig_t(i,j-1,k)) then !>>> La plus petite couche est agrandie CSP2

       k1=0 ; k3=k  

       depth_w(i,j-1,k1)=depth_w(i,j-1,k+1) &
!                          -dz_t(i,j,k,1)             ! dz methode 1
                         -dsig_t(i,j,k)*hz_w(i,j-1,1) ! dz methode 2

! interpolation ordre 1
!        rhp_t(i,j-1,k1)=  &
!        rhp_t(i,j-1,k+1)-(rhp_t(i,j-1,k+1)  -rhp_t(i,j-1,k)) &
!                      /(depth_t(i,j-1,k+1)-depth_t(i,j-1,k)) &
!                      *(depth_t(i,j-1,k+1)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1)))
! interpolation ordre 2 !29-10-21
      rhp_t(i,j-1,k1)=( &
         (  (  rhp_t(i,j-1,k+1)                                           ) &
           /(depth_t(i,j-1,k+1)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))) &
           -(                                            -rhp_t(i,j-1,k))   &
           /max(0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))-depth_t(i,j-1,k),small1)   &
         ) /(depth_t(i,j-1,k+1)-depth_t(i,j-1,k))                           &
        -(  (  rhp_t(i,j-1,k+2)  -rhp_t(i,j-1,k+1)) &
           /(depth_t(i,j-1,k+2)-depth_t(i,j-1,k+1)) &
           -(  rhp_t(i,j-1,k+1)                   ) &
           /(depth_t(i,j-1,k+1)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))) &
         ) /(depth_t(i,j-1,k+2)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))) &
                      ) &
                     /( &
          1./(depth_t(i,j-1,k+1)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))) &
            /(depth_t(i,j-1,k+2)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))) &
         +1./(depth_t(i,j-1,k+1)-0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))) &
            /(depth_t(i,j-1,k+1)-depth_t(i,j-1,k)) &
         +1./(depth_t(i,j-1,k+1)-depth_t(i,j-1,k)) &
        /max(0.5*(depth_w(i,j-1,k+1)+depth_w(i,j-1,k1))-depth_t(i,j-1,k),small1)  &
                      )

             anyv3d(i,j-1,k1 ,id_prs)=    &
             anyv3d(i,j-1,k+1,id_prs)     &
        +grav*rhp_t(i,j-1,k1)             &
!             *dz_t(i,j,k,1)             ! dz methode 1
            *dsig_t(i,j,k)*hz_w(i,j-1,1) ! dz methode 2
    
      else                                   !>>>

       k3=0 ; k1=k  

       depth_w(i,j,k3)=depth_w(i,j,k+1) &
!                        -dz_t(i,j-1,k,1)           ! dz methode 1
                       -dsig_t(i,j-1,k)*hz_w(i,j,1) ! dz methode 2

! interpolation ordre 1
!        rhp_t(i,j,k3)=  &
!        rhp_t(i,j,k+1)-(rhp_t(i,j,k+1)  -rhp_t(i,j,k)) &
!                    /(depth_t(i,j,k+1)-depth_t(i,j,k)) &
!                    *(depth_t(i,j,k+1)-0.5*(depth_w(i,j,k+1)+depth_w(i,j,k3)))
! interpolation ordre 2 !29-10-21
      rhp_t(i  ,j,k3)=( &
         (  (  rhp_t(i  ,j,k+1)                                           ) &
           /(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
           -(                                            -rhp_t(i  ,j,k))   &
        /max(0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))-depth_t(i  ,j,k),small1)   &
         ) /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k))                           &
        -(  (  rhp_t(i  ,j,k+2)  -rhp_t(i  ,j,k+1)) &
           /(depth_t(i  ,j,k+2)-depth_t(i  ,j,k+1)) &
           -(  rhp_t(i  ,j,k+1)                   ) &
           /(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
         ) /(depth_t(i  ,j,k+2)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
                      ) &
                     /( &
          1./(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
            /(depth_t(i  ,j,k+2)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
         +1./(depth_t(i  ,j,k+1)-0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))) &
            /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k)) &
         +1./(depth_t(i  ,j,k+1)-depth_t(i  ,j,k)) &
        /max(0.5*(depth_w(i  ,j,k+1)+depth_w(i  ,j,k3))-depth_t(i  ,j,k),small1)  &
                      )



             anyv3d(i,j,k3 ,id_prs)=    &
             anyv3d(i,j,k+1,id_prs)     &
        +grav*rhp_t(i,j,k3)             &
!             *dz_t(i,j-1,k,1)           ! dz methode 1
            *dsig_t(i,j-1,k)*hz_w(i,j,1) ! dz methode 2

      endif                                  !>>>


                               xy_v(i,j,1)=inv_rho*( & !PMXPMX> option deux couches distinctes
!     presgrad_v(i,j,1:kmerged_v(i,j)-1,1)=inv_rho*( & !PMXPMX> option deux couches distinctes
!     presgrad_v(i,j,1:kmerged_v(i,j)  ,1)=inv_rho*( & !PMXPMX> option deux couches regroupees !29-10-21

      +0.5*( anyv3d(i,j+jpv,k4,id_prs)-anyv3d(i,j-jmv,k2,id_prs)   & ! delta P(bottom+1)
            +anyv3d(i,j+jpv,k3,id_prs)-anyv3d(i,j-jmv,k1,id_prs) ) & ! delta P(bottom)

         +0.5*(depth_w(i,j+jpv,k4)-depth_w(i,j-jmv,k2)    &          ! Z(bottom+1)
              +depth_w(i,j+jpv,k3)-depth_w(i,j-jmv,k1))*( & !ooo>    ! Z(bottom)

        (  anyv3d(i,j+jpv,k3,id_prs)-anyv3d(i,j+jpv,k4,id_prs) &
          +anyv3d(i,j-jmv,k1,id_prs)-anyv3d(i,j-jmv,k2,id_prs) &
          )/(depth_w(i,j+jpv,k4)-depth_w(i,j+jpv,k3)           &
            +depth_w(i,j-jmv,k2)-depth_w(i,j-jmv,k1))          &

                                                      ) & !ooo>


           -grav*0.25*(depth_w(i,j+jpv,k4)     & ! < depth_v(i,j,k) >
                      +depth_w(i,j-jmv,k2)     &
                      +depth_w(i,j+jpv,k3)     &
                      +depth_w(i,j-jmv,k1))    &
                               *( rhpzavr_w(i,j+jpv)     &
                                 -rhpzavr_w(i,j-jmv))  ) & !PMXPMX>
! lignes suivantes pour retenir 100% du calcul si dz(kmerged-1)>=dz(kmerged)
! et prendre la valeur de PGF(kmerged) si dz(kmerged-1)=0: !03-11-21
           *pgfratio_v(i,j)+(1.-pgfratio_v(i,j))*presgrad_v(i,j,kmerged_v(i,j),1)
!          *min(1.,dz_v(i,j,kmerged_v(i,j)-1,1)/dz_v(i,j,kmerged_v(i,j),1))  &
!      +(1.-min(1.,dz_v(i,j,kmerged_v(i,j)-1,1)/dz_v(i,j,kmerged_v(i,j),1)))*presgrad_v(i,j,kmerged_v(i,j),1)

!     endif                                !m°v°m> !option deux couches regroupees  !29-10-21
      enddo ; enddo
      do j=2,jmax ; do i=2,imax-1
       do k=1,kmerged_v(i,j)-1 
        presgrad_v(i,j,k,1)=xy_v(i,j,1)
       enddo
      enddo ; enddo

!....................................
! The TKE scheme is coming next
! rhp_t possibly needs to be redefined if the level of reference of the potential
! density in the TKE scheme is different from that used in the PGF:

!       if(eos_author/=0) then
! Systematiquement calculer la densite potentielle maintenant que celle-ci est egalement !24-04-16
! utilisee dans le calcul de omega_w(:,:,kmax+1) et le calcul de rho.cp dans vertmix.F90
       if(eos_author==0) then !>>>>> !25-04-16
           call equation_of_state_linear(now)
       else                   !>>>>> !25-04-16
        if(eos_tkezref>=0.) then !-------->
          call equation_of_state_potzref2                 ! z reference = eos_tkezref !20-06-14
        else                     !-------->
          call equation_of_state_potloc2                  ! reference is current depth
        endif                    !-------->
       endif                  !>>>>> !25-04-16

! Marquer "disponibles" les tableaux generiques 
      anyv3d(-1,-1,0,1:4)=0. ; anyv3d(-1,-1,0,6:7)=0.
      anyv3d(-1,-1,0,id_prs)=0.

      end subroutine pressure_gradient_vqs

!.....................................................................................

      subroutine pressuregradient_remap_prs !29-10-21
      use module_principal ; use module_parallele
      implicit none

! Deduire la pression au niveau k-1 A partir des 4 pressions des niveaux k-2,k,k+1,k+2

! Methode 1: bloc 1 utilise k-1,k,k+1,k+2 
!            bloc 2 utilise k-2,k-1,k,k+1 
        do j=1,jmax ; do i=1,imax

        k=kmerged_t(i,j)+1

!       anyv3d(i,j,k-1,id_prs)=                                          &
!       anyv3d(i,j,k-1,id_prs)*(1-max(kmerged_t(i,j)-kmin_w(i,j),0))     &
!       +max(kmerged_t(i,j)-kmin_w(i,j),0)*0.25*( & !m°v°m>


        anyv3d(i,j,k-1,id_prs)=     &
        anyv3d(i,j,k-1,id_prs)*(1-max(kmerged_t(i,j)-kmin_w(i,j),0))     &
                                 +max(kmerged_t(i,j)-kmin_w(i,j),0)*( & !oooo>

           2./( depth_w(i,j,k+2)+depth_w(i,j,k  ) &
               -depth_w(i,j,k+1)-depth_w(i,j,k-1))*( & !AAA>

           2./(depth_w(i,j,k+2)-depth_w(i,j,k))*( & !BBB>
         ( anyv3d(i,j,k+2,id_prs)-anyv3d(i,j,k+1,id_prs)) &
        /(depth_w(i,j,k+2)      -depth_w(i,j,k+1)       ) &
        -( anyv3d(i,j,k+1,id_prs)-anyv3d(i,j,k  ,id_prs)) &
        /(depth_w(i,j,k+1)      -depth_w(i,j,k  )       ) &
                                                ) & !BBB>

          -2./(depth_w(i,j,k+1)-depth_w(i,j,k-1))*( & !CCC>
         ( anyv3d(i,j,k+1,id_prs)-anyv3d(i,j,k  ,id_prs)) &
        /(depth_w(i,j,k+1)      -depth_w(i,j,k  )       ) &
        -( anyv3d(i,j,k  ,id_prs)                       ) &
        /(depth_w(i,j,k  )      -depth_w(i,j,k-1)       ) &
                                                  ) & !CCC>

                                                   ) & !AAA>

          -2./( depth_w(i,j,k+1)+depth_w(i,j,k-1) &
               -depth_w(i,j,k  )-depth_w(i,j,k-2))*( & !AAA>

           2./(depth_w(i,j,k+1)-depth_w(i,j,k-1))*( & !BBB>
         ( anyv3d(i,j,k+1,id_prs)-anyv3d(i,j,k  ,id_prs)) &
        /(depth_w(i,j,k+1)      -depth_w(i,j,k  )       ) &
        -( anyv3d(i,j,k  ,id_prs)                       ) &
        /(depth_w(i,j,k  )      -depth_w(i,j,k-1)       ) &
                                                ) & !BBB>

          -2./(depth_w(i,j,k  )-depth_w(i,j,k-2))*( & !CCC>
         ( anyv3d(i,j,k  ,id_prs)                       ) &
        /(depth_w(i,j,k  )      -depth_w(i,j,k-1)       ) &
        -(                       -anyv3d(i,j,k-2,id_prs)) &
        /(depth_w(i,j,k-1)      -depth_w(i,j,k-2)       ) &
                                                  ) & !CCC>

                                                   ) & !AAA>

                                                                    ) & !oooo>
            /( & !>>>>

           4./(depth_w(i,j,k  )-depth_w(i,j,k-1)) &
             /(depth_w(i,j,k+1)-depth_w(i,j,k-1)) &
            /( depth_w(i,j,k+2)+depth_w(i,j,k  )  &
              -depth_w(i,j,k+1)-depth_w(i,j,k-1)) &         

          +4./( depth_w(i,j,k+1)+depth_w(i,j,k-1)  &
               -depth_w(i,j,k  )-depth_w(i,j,k-2))*( & !pmx>

           1./(depth_w(i,j,k  )-depth_w(i,j,k-1))*( & !xmp>
                            1./(depth_w(i,j,k+1)-depth_w(i,j,k-1)) &
                           +1./(depth_w(i,j,k  )-depth_w(i,j,k-2)) &
                                                  ) & !xmp

          +1./(depth_w(i,j,k-1)-depth_w(i,j,k-2)) &
             /(depth_w(i,j,k  )-depth_w(i,j,k-2)) &

                                                   ) & !pmx>
             )   !>>>>

              rhp_t(i,j,k-1)=invgrav*(anyv3d(i,j,k-1,id_prs)-anyv3d(i,j,k  ,id_prs))/dz_t(i,j,k-1,1)
              rhp_t(i,j,k-2)=invgrav*(anyv3d(i,j,k-2,id_prs)-anyv3d(i,j,k-1,id_prs))/dz_t(i,j,k-2,1)

        enddo ; enddo

      end subroutine pressuregradient_remap_prs
