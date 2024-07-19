      subroutine z_levels(ichoix)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 296 - last update: 28-02-21
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='z_levels'
       subroutinedescription= &
      'Updates the depth of grid points and stores the result in depth_'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ Computes the depth of cell boxes

!______________________________________________________________________
!Version date      Description des modifs
!        01/11/01: inversion de l'ordre des boucles: J passe avant I
!        09/03/09: Condition à la limite regroupee dans obc_z pour clarifier
!                  le travail de parallelisation
!        06/04/09: calcul basé sur hz_z, bornes boucle Z, parallelisation
! 2009.3 30-09-09: suppression sigma_c et amenagement des boucles en
!                  consequence
!        07-10-09: la generalisation du maillage vertical conduit à
!                  calculer depth_z à partir de dz_c
! S.26   09-04-13: c'est desormais ici que l'on calcule zeroleveldepth
!                  pour pouvoir utiliser rapidement cette info dans de
!                  nouvelles fonctions comme par ex les vagues
!        02-11-14: blindage wetdry
!        05-08-15  - menage
!                  - methode vst continue 
!        21-04-17  coordonnee s-z
!        17-10-17  coordonnee s-z
!        07-10-18  call obc_depth(2)
!        16-11-18  suite et complement du point precedent
! v285   05-06-20  depth_u depth_v calcul compact
!                  suppression zeroleveldepth_u et _v
!        06-06-20  depth_w calculE avec sigma_w
! v296   28-02-21  nouveau calcul empechant z(k)-z(k-1)<=0
!...............................................................................
!    _________                    .__                  .__                     ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

! depth_w at box interfaces starting from bottom:
#ifdef bidon
! commencer depuis la surface sinon en coordonnee sigma-step ca merde
       do j=0,jmax+1                                               !06/04/09
        do i=0,imax+1
          depth_w(i,j,kmax+1)=hz_w(i,j,1)-h_w(i,j)
        enddo
       enddo

       do k=kmax,1,-1
        do j=0,jmax+1                                                   !06/04/09
         do i=0,imax+1
          depth_w(i,j,k)=depth_w(i,j,k+1)-dz_t(i,j,k,1)
         enddo
        enddo
       enddo
#endif

! depth_w !06-06-20
      do k=1,kmax+1 ; do j=0,jmax+1 ; do i=0,imax+1
! Attention avec la formule du 06-06-20 z(k)-z(k-1) pouvait etre nul ou
! negatif car ssh_int_w n'est pas seuillE et peut passer sous h_w.
!      depth_w(i,j,k)=ssh_int_w(i,j,1)*sigma_w(i,j,k)+h_w(i,j)*(sigma_w(i,j,k)-1.) !06-06-20 
! La nouvelle formule utilise hz qui est sujet A seuillage donc evite ce defaut:
       depth_w(i,j,k)=hz_w(i,j,1)*sigma_w(i,j,k)-h_w(i,j) !28-02-21
      enddo         ; enddo         ; enddo

! depth at mid-box level:
      do k=1,kmax
       do j=0,jmax+1                                                    !06/04/09
        do i=0,imax+1
          depth_t(i,j,k)=0.5*(depth_w(i,j,k)+depth_w(i,j,k+1))          !30-09-09
        enddo
       enddo
      enddo

! Vertical extrapolation required for the computation of turbulent length scales:
! Ces lignes sont deplacees avant call obc_depth(2) le !16-11-18
       do j=0,jmax+1 ; do i=0,imax+1
          depth_t(i,j,kmax+1)=depth_w(i,j,kmax+1) 
          depth_t(i,j,0     )=depth_w(i,j,1)            
       enddo         ; enddo

! continuitE mpi "z2" (i=-1 etc...) sur depth_t pour algo de correction de pente
! du terme de diffusion du quickest
      call obc_depth(2)  !07-10-18

! u grid point:
       do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1
          depth_u(i,j,k)=0.5*(depth_t(i,j,k)+depth_t(i-1,j,k)) !05-06-20
       enddo       ; enddo         ; enddo

! v grid point:
       do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1
          depth_v(i,j,k)=0.5*(depth_t(i,j,k)+depth_t(i,j-1,k)) !05-06-20
       enddo       ; enddo         ; enddo

      end subroutine z_levels
