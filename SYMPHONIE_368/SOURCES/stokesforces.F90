      subroutine stokesforces
      use module_principal
      implicit none
!______________________________________________________________________
! SYMPHONIE OCEAN MODEL
! release S.26  - last update: 11-09-18
!______________________________________________________________________
!...............................................................................
! Version date      Description des modifs
! 2010.8  03-05-10  Mise en service
!         21-05-10  Ajuster rigoureusement les boucles pour eviter division par
!                   zero
! 2010.25 23-02-12  Ajout Sshear 3D
! S.26    09-04-13  vitesse prises au temps t_=now
!                   ajout subroutine stokesforces_vortex_external
!         30-05-17  ajout subroutine stokesforces_vortex_2dmode pour les cas kmax=1
!         11-09-18  modif C.L. fond et surface + reorganisation boucle vericale
!                   pour continuite mpi
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
#ifdef synopsis
       subroutinetitle='stokesforces'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call stokesforces_vortex
      call stokesforces_sshear

! Compute the depth averaged forcing terms for the external mode computation:
      call stokesforces_depth_averaged

      end subroutine stokesforces

!...........................................................................

      subroutine stokesforces_vortex
      use module_principal
      use module_parallele
      implicit none
      integer t_
#ifdef synopsis
       subroutinetitle='stokesforces_vortex'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      t_=now ! 09-04-13

!$ u component of the velocity:

      do k=1,kmax   ! debut boucle 1 K

      kp1=min0(k+1,kmax)
      km1=max0(k-1,1   )

!$ Vortex forces:

      do j=2,jmax      !21-05-10
      do i=2,imax

! Stokes Vortex force 1: (2*Dx)*( dv/dx*Vs ) at p point:
      ustokesvortex_f(i,j,1)=                                         &

          (vel_v(i,j,k,t_)-  vel_v(i-1,j,k,t_)                        &

       -(    depth_v(i  ,j  ,k)   -depth_v(i-1,j  ,k)  )              &
      *0.5*(  (vel_v(i  ,j  ,kp1,t_)-vel_v(i  ,j  ,km1,t_))           &
           /(depth_v(i  ,j  ,kp1) -depth_v(i  ,j  ,km1)  )            &
             +(vel_v(i-1,j  ,kp1,t_)-vel_v(i-1,j  ,km1,t_))           &
           /(depth_v(i-1,j  ,kp1) -depth_v(i-1,j  ,km1)  ) ))         &

                    *(   velstokes_v(i,j,k,1)+ velstokes_v(i-1,j,k,1))

      enddo
      enddo

      do j=2,jmax-1     !21-05-10
      do i=1,imax

! Stokes Vortex force 2: (2*Dx)*( du/dx*Us ) at t point:
      ustokesvortex_t(i,j,1)=                                         &

          (vel_u(i+1,j,k,t_)-vel_u(i,j,k,t_)                          &

       -(    depth_u(i+1,j  ,k)   -depth_u(i  ,j  ,k)  )              &
      *0.5*(  (vel_u(i  ,j  ,kp1,t_)-vel_u(i  ,j  ,km1,t_))           &
           /(depth_u(i  ,j  ,kp1) -depth_u(i  ,j  ,km1)  )            &
             +(vel_u(i+1,j  ,kp1,t_)-vel_u(i+1,j  ,km1,t_))           &
           /(depth_u(i+1,j  ,kp1) -depth_u(i+1,j  ,km1)  ) ))         &

                    *(   velstokes_u(i+1,j,k,1)+ velstokes_u(i,j,k,1))

      enddo
      enddo



      do j=2,jmax-1
      do i=2,imax

! Vortex force: Dx*(du/dx*Us+dv/dx*Vs)
      stokesforces_u(i,j,k)=                                     &
           0.25*( ustokesvortex_t(i  ,j  ,1)                     &
                 +ustokesvortex_t(i-1,j  ,1)                     &
                 +ustokesvortex_f(i  ,j  ,1)                     &
                 +ustokesvortex_f(i  ,j+1,1))

      enddo
      enddo
      enddo ! fin boucle 1 K

!$ v component of the velocity:

      do k=1,kmax   ! debut boucle 2 K

      kp1=min0(k+1,kmax)
      km1=max0(k-1,1   )

!$ Compute the full velocity diffusive fluxes & rotation terms:

!............
! v equation:

!$ Vortex forces: 3d terms

      do j=2,jmax    !21-05-10
      do i=2,imax

!............
! v equation:
! Stokes Vortex force 3:(2*Dy)*( du/dy*Us )at p point:
      vstokesvortex_f(i,j,1)=                                         &

          (vel_u(i,j,k,t_)-  vel_u(i,j-1,k,t_)                        &

       -(    depth_u(i  ,j  ,k)   -depth_u(i  ,j-1,k)  )              &
      *0.5*(  (vel_u(i  ,j  ,kp1,t_)-vel_u(i  ,j  ,km1,t_))           &
           /(depth_u(i  ,j  ,kp1) -depth_u(i  ,j  ,km1)  )            &
             +(vel_u(i  ,j-1,kp1,t_)-vel_u(i  ,j-1,km1,t_))           &
           /(depth_u(i  ,j-1,kp1) -depth_u(i  ,j-1,km1)  ) ))         &

                    *(   velstokes_u(i,j,k,1)+ velstokes_u(i,j-1,k,1))

      enddo
      enddo

      do j=1,jmax       !21-05-10
      do i=2,imax-1

! Stokes Vortex force 4: (2*Dy)*( dv/dy*Vs ) at p point:
      vstokesvortex_t(i,j,1)=                                         &

          (vel_v(i,j+1,k,t_)-  vel_v(i,j,k,t_)                        &

       -(    depth_v(i  ,j+1,k)   -depth_v(i  ,j  ,k)  )              &
      *0.5*(  (vel_v(i  ,j  ,kp1,t_)-vel_v(i  ,j  ,km1,t_))           &
           /(depth_v(i  ,j  ,kp1) -depth_v(i  ,j  ,km1)  )            &
             +(vel_v(i  ,j+1,kp1,t_)-vel_v(i  ,j+1,km1,t_))           &
           /(depth_v(i  ,j+1,kp1) -depth_v(i  ,j+1,km1)  ) ))         &

                    *(   velstokes_v(i,j+1,k,1) +velstokes_v(i,j,k,1))

      enddo
      enddo


      do j=2,jmax
      do i=2,imax-1


! Vortex Force: Dy*( du/dy*Us+dv/dy*Vs )
      stokesforces_v(i,j,k)=                  &
           0.25*( vstokesvortex_t(i  ,j  ,1)  &
                 +vstokesvortex_t(i  ,j-1,1)  &
                 +vstokesvortex_f(i  ,j  ,1)  &
                 +vstokesvortex_f(i+1,j  ,1))

      enddo
      enddo
      enddo ! fin boucle 2 K


      end subroutine stokesforces_vortex

!...................................................................................

      subroutine stokesforces_sshear
      use module_principal
      implicit none
      double precision vsdvdz_,dsheardi_,dsheardj_
      integer t_
#ifdef synopsis
       subroutinetitle='stokesforces_sshear'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      t_=now ! before

! Sshear Force:(Ardhuin et al, 2008, Eq 40)
! Note that we assumed that the vertical velocity term could be neglected


! commentE le 11-09-18 pour continuitE mpi !11-09-18
! Extrapolate the velocity at kmax+1 and 0 levels:
!     do j=1,jmax
!     do i=1,imax+1
!            vel_u(i,j,kmin_u(i,j)-1,t_)=      vel_u(i,j,kmin_u(i,j),t_)
!            vel_u(i,j,kmax+1       ,t_)=      vel_u(i,j,kmax       ,t_)
!      velstokes_u(i,j,kmin_u(i,j)-1,1) =velstokes_u(i,j,kmin_u(i,j),1)
!      velstokes_u(i,j,kmax+1       ,1) =velstokes_u(i,j,kmax       ,1)
!     enddo
!     enddo
!     do j=1,jmax+1
!     do i=1,imax
!            vel_v(i,j,kmin_v(i,j)-1,t_)=      vel_v(i,j,kmin_v(i,j),t_)
!            vel_v(i,j,kmax+1       ,t_)=      vel_v(i,j,kmax       ,t_)
!      velstokes_v(i,j,kmin_v(i,j)-1,1) =velstokes_v(i,j,kmin_v(i,j),1)
!      velstokes_v(i,j,kmax+1       ,1) =velstokes_v(i,j,kmax       ,1)
!     enddo
!     enddo

! Compute Sshear at "t" nodes:
      do j=1,jmax
      do i=1,imax
! reset Sshear vertical summation:
      anyv3d(i,j,kmax+1,1)=0.   ! Sshear (neglecting ws term)
      enddo
      enddo

! Boucles reorganisEes le !11-09-18 
      do k=kmax,1,-1
      kp1=min(k+1,kmax)
      km1=max(k-1,1)
       do j=1,jmax
       do i=1,imax

!vsdvdz_=-ustokes*d(ulag)/dz-vstokes*d(vlag)/dz
!      vsdvdz_=                                                      &
       anyv3d(i,j,k,1)=anyv3d(i,j,k+1,1)+dz_t(i,j,k,now)*            &  ! Sshear=VertSum(-ustokes*d(ulag)/dz-vstokes*d(vlag)/dz)
      (-0.5*( velstokes_u(i+ipt,j,k,1)*(velstokes_u(i+ipt,j,kp1,1)   &
                                             +vel_u(i+ipt,j,kp1,t_)  &
                                       -velstokes_u(i+ipt,j,km1,1)   &
                                             -vel_u(i+ipt,j,km1,t_)) &
             +velstokes_u(i-imt,j,k,1)*(velstokes_u(i-imt,j,kp1,1)   &
                                             +vel_u(i-imt,j,kp1,t_)  &
                                       -velstokes_u(i-imt,j,km1,1)   &
                                             -vel_u(i-imt,j,km1,t_)) &
             +velstokes_v(i,j+jpt,k,1)*(velstokes_v(i,j+jpt,kp1,1)   &
                                             +vel_v(i,j+jpt,kp1,t_)  &
                                       -velstokes_v(i,j+jpt,km1,1)   &
                                             -vel_v(i,j+jpt,km1,t_)) &
             +velstokes_v(i,j-jmt,k,1)*(velstokes_v(i,j-jmt,kp1,1)   &
                                             +vel_v(i,j-jmt,kp1,t_)  &
                                       -velstokes_v(i,j-jmt,km1,1)   &
                                             -vel_v(i,j-jmt,km1,t_)))&
       /(depth_t(i,j,kp1)-depth_t(i,j,km1)))


!      anyv3d(i,j,k,1)=anyv3d(i,j,k+1,1)+vsdvdz_*dz_t(i,j,k,now)  ! Sshear=VertSum(-ustokes*d(ulag)/dz-vstokes*d(vlag)/dz)

       enddo
       enddo
      enddo

      do j=1,jmax
      do i=1,imax
      anyv3d(i,j,0,1)=anyv3d(i,j,1,1)
      enddo
      enddo

! Compute dx*d(Sshear/dx)=d(Sshear)/di-dz/di*d(Sshear)/dz at "u" nodes
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
       dsheardi_=                                                  &
           anyv3d(i+ipu,j,k,1)-anyv3d(i-imu,j,k,1)                 & !  d(Shear)/di
       -( depth_t(i+ipu,j,k)                                       & ! -dz/di*
         -depth_t(i-imu,j,k) )*0.5*(anyv3d(i+ipu,j,k+1,1)          & !  d(Shear)
                                   +anyv3d(i-imu,j,k+1,1)          &
                                   -anyv3d(i+ipu,j,k-1,1)          &
                                   -anyv3d(i-imu,j,k-1,1))         &
                                               /(depth_u(i,j,k+1)  & ! /dz
                                                -depth_u(i,j,k-1))

        stokesforces_u(i,j,k)=                           &
        stokesforces_u(i,j,k)                            &
       -dsheardi_                                        &
       +grav*(sshstokes_w(i+ipu,j)-sshstokes_w(i-imu,j))   ! -dJ/di -dSsurf/di

      enddo
      enddo
      enddo

! Compute dy*d(Sshear/dy)=d(Sshear)/dj-dz/dj*d(Sshear)/dz at "v" nodes
      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
       dsheardj_=                                                  &
           anyv3d(i,j+jpv,k,1)-anyv3d(i,j-jmv,k,1)                 & !  d(Shear)/dj
       -( depth_t(i,j+jpv,k)                                       & ! -dz/dj*
         -depth_t(i,j-jmv,k) )*0.5*(anyv3d(i,j+jpv,k+1,1)          & !  d(Shear)
                                   +anyv3d(i,j-jmv,k+1,1)          &
                                   -anyv3d(i,j+jpv,k-1,1)          &
                                   -anyv3d(i,j-jmv,k-1,1))         &
                                               /(depth_v(i,j,k+1)  & ! /dz
                                                -depth_v(i,j,k-1))

        stokesforces_v(i,j,k)=                          &
        stokesforces_v(i,j,k)                           &
       -dsheardj_                                       &
       +grav*(sshstokes_w(i,j+jpv)-sshstokes_w(i,j-jmv))   ! -dJ/dj -dSsurf/dj

      enddo
      enddo
      enddo

      end subroutine stokesforces_sshear

!.........................................................................

      subroutine stokesforces_depth_averaged
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='stokesforces_depth_averaged'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Note that the shear force does not need to be recomputed during the
! external mode cycle since d(velbar)/dz=0
! On the other hand the vortex force is possibly updated during the external
! mode cycle with the updated depth-averaged current:
! VortexForce(t+dt)=VortexForce(t)+Delta_VortexForce(delta_velbar)

      do j=2,jmax-1
      do i=2,imax
        stokesforces3d2d_u(i,j)=0.
      enddo
      enddo
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
       stokesforces3d2d_u(i,j)=                                  &
       stokesforces3d2d_u(i,j)+                                  &
                   mask_u(i,j,k)                                 &
                    *dz_u(i,j,k,now)                             &
          *stokesforces_u(i,j,k)
      enddo
      enddo
      enddo
      do j=2,jmax-1
      do i=2,imax
       stokesforces3d2d_u(i,j)=stokesforces3d2d_u(i,j)/hz_u(i,j,1)
      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1
        stokesforces3d2d_v(i,j)=0.
      enddo
      enddo
      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
       stokesforces3d2d_v(i,j)=                                  &
       stokesforces3d2d_v(i,j)+                                  &
                   mask_v(i,j,k)                                 &
                    *dz_v(i,j,k,now)                             &
          *stokesforces_v(i,j,k)
      enddo
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
       stokesforces3d2d_v(i,j)=stokesforces3d2d_v(i,j)/hz_v(i,j,1)
      enddo
      enddo

      end subroutine stokesforces_depth_averaged

!...........................................................................

      subroutine stokesforces_vortex_external
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='stokesforces_vortex_external'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!$ u component of the velocity:

!$ Vortex forces:

      do j=2,jmax      !21-05-10
      do i=2,imax

! Stokes Vortex force 1: (2*Dx)*( dv/dx*Vs ) at p point:
      ustokesvortex_f(i,j,1)=                                  &

               ( velbar_v(i,j,1)    -  velbar_v(i-1,j,1)      &
                -velbar_v(i,j,0)    +  velbar_v(i-1,j,0))     &
         *(velbarstokes_v(i,j,1)+velbarstokes_v(i-1,j,1))

      enddo
      enddo

      do j=2,jmax-1     !21-05-10
      do i=1,imax

! Stokes Vortex force 2: (2*Dx)*( du/dx*Us ) at t point:
      ustokesvortex_t(i,j,1)=                                 &

               ( velbar_u(i+1,j,1)     - velbar_u(i,j,1)     &
                -velbar_u(i+1,j,0)     + velbar_u(i,j,0))    &
         *(velbarstokes_u(i+1,j,1)+velbarstokes_u(i,j,1))

      enddo
      enddo



      do j=2,jmax-1
      do i=2,imax

! Vortex force: Dx*(du/dx*Us+dv/dx*Vs)
      frozenterm3d_u(i,j,1)=                                     &
      frozenterm3d_u(i,j,1)                                      &
          -0.25*( ustokesvortex_t(i  ,j  ,1)                     &
                 +ustokesvortex_t(i-1,j  ,1)                     &
                 +ustokesvortex_f(i  ,j  ,1)                     &
                 +ustokesvortex_f(i  ,j+1,1))*wetmask_u(i,j)
      enddo
      enddo

!$ v component of the velocity:

!$ Compute the full velocity diffusive fluxes & rotation terms:

!............
! v equation:

!$ Vortex forces: 3d terms

      do j=2,jmax    !21-05-10
      do i=2,imax

!............
! v equation:
! Stokes Vortex force 3:(2*Dy)*( du/dy*Us )at p point:
      vstokesvortex_f(i,j,1)=                                  &

               ( velbar_u(i,j,1)    -  velbar_u(i,j-1,1)       &
                -velbar_u(i,j,0)    +  velbar_u(i,j-1,0))      &
         *(velbarstokes_u(i,j,1)+velbarstokes_u(i,j-1,1))

      enddo
      enddo

      do j=1,jmax       !21-05-10
      do i=2,imax-1

! Stokes Vortex force 4: (2*Dy)*( dv/dy*Vs ) at p point:
      vstokesvortex_t(i,j,1)=                                         &

               ( velbar_v(i,j+1,1)    -  velbar_v(i,j,1)              &
                -velbar_v(i,j+1,0)    +  velbar_v(i,j,0))             &
         *(velbarstokes_v(i,j+1,1)+velbarstokes_v(i,j,1))

      enddo
      enddo


      do j=2,jmax
      do i=2,imax-1


! Vortex Force: Dy*( du/dy*Us+dv/dy*Vs )
      frozenterm3d_v(i,j,1)=                                  &
      frozenterm3d_v(i,j,1)                                   &
          -0.25*( vstokesvortex_t(i  ,j  ,1)                  &
                 +vstokesvortex_t(i  ,j-1,1)                  &
                 +vstokesvortex_f(i  ,j  ,1)                  &
                 +vstokesvortex_f(i+1,j  ,1))*wetmask_v(i,j)
      enddo
      enddo

      end subroutine stokesforces_vortex_external

!...................................................................................

      subroutine stokesforces_vortex_2dmode !30-05-17
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='stokesforces_vortex_2dmode'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! CAS DU MODELE 2D (kmax=1)
! Cette version remet A jour le forcage une fois par pas de temps
! principal. Une version remettant A jour au pas de temps du mode
! externe pourrait etre testee


!$ u component of the velocity:

!$ Vortex forces:

      do j=2,jmax      !21-05-10
      do i=2,imax

! Stokes Vortex force 1: (2*Dx)*( dv/dx*Vs ) at p point:
      ustokesvortex_f(i,j,1)=                                  &
               ( velbar_v(i,j,1)    -  velbar_v(i-1,j,1))     &
         *(velbarstokes_v(i,j,1)+velbarstokes_v(i-1,j,1))

      enddo
      enddo

      do j=2,jmax-1     !21-05-10
      do i=1,imax

! Stokes Vortex force 2: (2*Dx)*( du/dx*Us ) at t point:
      ustokesvortex_t(i,j,1)=                                 &
               ( velbar_u(i+1,j,1)     - velbar_u(i,j,1))    &
         *(velbarstokes_u(i+1,j,1)+velbarstokes_u(i,j,1))

      enddo
      enddo



      do j=2,jmax-1
      do i=2,imax

! Vortex force: Dx*(du/dx*Us+dv/dx*Vs)
      pres3d2d_u(i,j)=  &
      pres3d2d_u(i,j)   &
      +wetmask_u(i,j)*( & !pmxpmx>
          -0.25*( ustokesvortex_t(i  ,j  ,1)                       &
                 +ustokesvortex_t(i-1,j  ,1)                       &
                 +ustokesvortex_f(i  ,j  ,1)                       &
                 +ustokesvortex_f(i  ,j+1,1))                      &
                 +grav*(sshstokes_w(i-imu,j)-sshstokes_w(i+ipu,j)) &
                      )   !pmxpmx>
      enddo
      enddo

!$ v component of the velocity:

!$ Compute the full velocity diffusive fluxes & rotation terms:

!............
! v equation:

!$ Vortex forces: 3d terms

      do j=2,jmax    !21-05-10
      do i=2,imax

!............
! v equation:
! Stokes Vortex force 3:(2*Dy)*( du/dy*Us )at p point:
      vstokesvortex_f(i,j,1)=                                  &
               ( velbar_u(i,j,1)    -  velbar_u(i,j-1,1))      &
         *(velbarstokes_u(i,j,1)+velbarstokes_u(i,j-1,1))

      enddo
      enddo

      do j=1,jmax       !21-05-10
      do i=2,imax-1

! Stokes Vortex force 4: (2*Dy)*( dv/dy*Vs ) at p point:
      vstokesvortex_t(i,j,1)=                                         &
               ( velbar_v(i,j+1,1)    -  velbar_v(i,j,1))             &
         *(velbarstokes_v(i,j+1,1)+velbarstokes_v(i,j,1))

      enddo
      enddo


      do j=2,jmax
      do i=2,imax-1


! Vortex Force: Dy*( du/dy*Us+dv/dy*Vs )
      pres3d2d_v(i,j)=  &
      pres3d2d_v(i,j)   &
      +wetmask_v(i,j)*( & !pmxpmx>
          -0.25*( vstokesvortex_t(i  ,j  ,1)                       &
                 +vstokesvortex_t(i  ,j-1,1)                       &
                 +vstokesvortex_f(i  ,j  ,1)                       &
                 +vstokesvortex_f(i+1,j  ,1))                      &   
                 +grav*(sshstokes_w(i,j-jmv)-sshstokes_w(i,j+jpv)) &
                      )   !pmxpmx>
      enddo
      enddo

      end subroutine stokesforces_vortex_2dmode

!...................................................................................
