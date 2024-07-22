










      subroutine omega
!______________________________________________________________________
! SYMPHONIE ocean model
! release 315 - last update: 10-12-21
!______________________________________________________________________
      use module_principal ; use module_parallele
      implicit none


!...............................................................................
! Version date      Description des modifications
!         14/07/01: passage à la coordonnée sigma généralisée
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         26/12/02: amenagements remise en service du schema forward
!         07/08/06: Modes externes parallelles
!         17/04/07: Passage aux coordonnees curvilignes
!         28/05/07: ajout de la multiplication par   mask_c(k-1)
! 2009.3  30-09-09  omega_x et omega_y sont supprimés
!         16-11-09  calculer la vitesse verticale "schema temporel forward"
!                   que si l'advection des traceurs passifs est activée
! 2010.7  16-02-10  vitesse forward = veldxdz(2)
!         28-02-10  - Le schema forward de la TKE impose de calculer omega(2)
!                   systematiquement
!                   - essai de suppression de *mask_t
! 2010.8  03-05-10  Avec nouveau schema forward pas besoin de omega(2)
!         18-05-10  Une astuce pour compenser les pertes de precision machine
!                   dans l'equation de continuité
!         21-05-10  En coordonnee hydride, dz sous le fond n'est pas exactement
!                   nul pour eviter les divisions par zero. La precision sur omega
!                   est donc renforcée par une multiplication par mask_t
! 2010.9  31-05-10  imc et cie remplacé par imt et cie....
! 2010.12 27-09-10  argument key supprimé
! S26     09-09-14  ajout use module_parallele
!         31-10-15  calcul de omega(kmax+1) pour bancs decouvrants
!         17-10-17  modif couches fusionnees
!         12-06-18  ajout subroutine omega_vertical_velocity
! v315    10-12-21  if(case_==2) ajoutE dans omega_vertical_velocity
!...............................................................................
!    _________                    .__                  .__             !   (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !   m°v°m
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................




!$ GEneric MOdelling shortname conventions (m for "minus" & p for "plus")
!     ipt=1 ; imt=0 ; jpt=1 ; jmt=0 ; kpt=1 ; kmt=0
!$:...........

! On calcule en kmax+1 car dans les bancs decouvrants ou la ssh est
! seuillee le bilan de masse peut parfois ne pas etre satisfait ce qui resulte
! en une valeur de omega(kmax+1) differente de ce qu'on stipule dans
! module_airseaflux. De part la condition de flux nul (gradient nul)
! associee a omega(kmax+1) ceci est en mesure d'eviter des valeurs
! etranges dans les zones decouvertes.
      do k=1,kmax ! kmax-1 ! 31-10-15
      do j=1,jmax
      do i=1,imax

             omega_w(i    ,j    ,k+kpt,1  )=                    &
             omega_w(i    ,j    ,k-kmt,1  )                     &
        -( veldydz_u(i+ipt,j    ,k    ,now)                     &
          -veldydz_u(i-imt,j    ,k    ,now)                     &
          +veldxdz_v(i    ,j+jpt,k    ,now)                     &
          -veldxdz_v(i    ,j-jmt,k    ,now)    )/dxdy_t(i,j)    &
        -(      dz_t(i    ,j    ,k    ,after )                  &
               -dz_t(i    ,j    ,k    ,before) )/dti_lp ! 17-10-17 *mask_t(i,j,k)   !21-05-10

      enddo
      enddo
      enddo

! La suite n'est pas utilisee

      end subroutine omega

!..............................................................

      subroutine omega_vertical_velocity(case_,istart_,iend_,jstart_,jend_)
      use module_principal ; use module_parallele
      implicit none
      integer case_,istart_,iend_,jstart_,jend_
      ksecu=0

! Cet algo est par Knut Klingbeil, Hans Burchard, Ocean Modelling 65 (2013) 64¿77 https://doi.org/10.1016/j.ocemod.2013.02.002
! Details dans
! https://docs.google.com/document/d/1P9YmC88Un2LPVI36eUNS0vY9uRph3bJc9JQ6-gy-8xo/edit
 
      if(case_==1) then !pmx>
       ksecu=1
! calculer depth_w en differents temps
       do j=0,jmax+1 ; do i=0,imax+1
        anyv3d(i,j,1,0:2)=-h_w(i,j)
       enddo         ; enddo
       do j=0,jmax+1 ; do i=0,imax+1
          do k=2,kmax+1
            anyv3d(i,j,k,0)=anyv3d(i,j,k-1,0)+dsig_t(i,j,k-1)*(h_w(i,j)+ssh_int_w(i,j,0))
            anyv3d(i,j,k,1)=anyv3d(i,j,k-1,0)+dsig_t(i,j,k-1)*(h_w(i,j)+ssh_int_w(i,j,1))
            anyv3d(i,j,k,2)=anyv3d(i,j,k-1,0)+dsig_t(i,j,k-1)*(h_w(i,j)+ssh_int_w(i,j,2))
!             if(i==imax/2.and.j==jmax/2)write(6,*)k,real(depth_w(i,j,k)),real(anyv3d(i,j,k,0)),real(ssh_int_w(i,j,0))
          enddo
       enddo         ; enddo
      
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

        anyvar3d(i,j,k)=( & !PMX>

                (1./(anyv3d(i,j,k+1,1)-anyv3d(i,j,k,1)))*( & ! 1./Dz

         (1./dxdy_t(i,j))*( & !ooo>
           veldydz_u(i+1,j,k,1)*0.25*(anyv3d(i,j,k,1)+anyv3d(i+1,j,k,1)+anyv3d(i,j,k+1,1)+anyv3d(i+1,j,k+1,1)) & !dy*dz*u*z  i+1
          -veldydz_u(i  ,j,k,1)*0.25*(anyv3d(i,j,k,1)+anyv3d(i-1,j,k,1)+anyv3d(i,j,k+1,1)+anyv3d(i-1,j,k+1,1)) & !dy*dz*u*z  i
          +veldxdz_v(i,j+1,k,1)*0.25*(anyv3d(i,j,k,1)+anyv3d(i,j+1,k,1)+anyv3d(i,j,k+1,1)+anyv3d(i,j+1,k+1,1)) & !dx*dz*v*z  j+1
          -veldxdz_v(i,j  ,k,1)*0.25*(anyv3d(i,j,k,1)+anyv3d(i,j-1,k,1)+anyv3d(i,j,k+1,1)+anyv3d(i,j-1,k+1,1)) & !dx*dz*v*z  j
                          ) & !ooo>

        +omega_w(i,j,k+1,1)*anyv3d(i,j,k+1,1) & !omega*z k+1
        -omega_w(i,j,k  ,1)*anyv3d(i,j,k  ,1) & !omega*z k

        +(          & !ttt>
            (anyv3d(i,j,k+1,2)-anyv3d(i,j,k,2))*0.5*(anyv3d(i,j,k+1,2)+anyv3d(i,j,k,2)) & ! dz*z t+1
           -(anyv3d(i,j,k+1,0)-anyv3d(i,j,k,0))*0.5*(anyv3d(i,j,k+1,0)+anyv3d(i,j,k,0)) & ! dz*z t+1
         )/dti_lp   & !ttt>

                                                         ) &
                         ) & !PMX>
                          *mask_t(i,j,kmax)+(1.-mask_t(i,j,kmax))*(-9999.)
       enddo       ; enddo       ; enddo
      endif             !pmx>

      if(case_==2) then !m°v°m> !10-12-21
       ksecu=1
       do k=1,kmax ; do j=jstart_,jend_ ; do i=istart_,iend_
        anyvar3d(i,j,k)=( & !PMX>
                                      (1./dz_t(i,j,k,1))*( & ! 1./Dz
         (1./dxdy_t(i,j))*( & !ooo>
           veldydz_u(i+1,j,k,1)*0.5*(depth_u(i,j,k)+depth_u(i+1,j,k)) & !dy*dz*u*z  i+1
          -veldydz_u(i  ,j,k,1)*0.5*(depth_u(i,j,k)+depth_u(i-1,j,k)) & !dy*dz*u*z  i
          +veldxdz_v(i,j+1,k,1)*0.5*(depth_v(i,j,k)+depth_v(i,j+1,k)) & !dx*dz*v*z  j+1
          -veldxdz_v(i,j  ,k,1)*0.5*(depth_v(i,j,k)+depth_v(i,j-1,k)) & !dx*dz*v*z  j
                          ) & !ooo>
        +omega_w(i,j,k+1,1)*depth_w(i,j,k+1) & !omega*z k+1
        -omega_w(i,j,k  ,1)*depth_w(i,j,k  ) & !omega*z k
        +(          & !ttt>
            dz_t(i,j,k,after) *(hz_w(i,j,after) *0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))-h_w(i,j)) & !dz(t+1)*z(t+1)
           -dz_t(i,j,k,before)*(hz_w(i,j,before)*0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))-h_w(i,j)) & !dz(t-1)*z(t-1)
         )/dti_lp   & !ttt>
                                                         ) & ! 1./Dz
                         ) & !PMX>
                          *mask_t(i,j,kmax)+(1.-mask_t(i,j,kmax))*(-9999.)
       enddo       ; enddo       ; enddo
      endif             !m°v°m> !10-12-21

      if(ksecu==0) &
      stop ' Error omega_vertical_velocity case not recognized'
      end subroutine omega_vertical_velocity

!..............................................................
