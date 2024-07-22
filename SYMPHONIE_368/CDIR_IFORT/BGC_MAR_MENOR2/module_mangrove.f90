










      module module_mangrove
!______________________________________________________________________
! SYMPHONIE ocean model 
! release S26.1 - last update: 12-07-18
!______________________________________________________________________
! version date      Description des modifications
!         01-05-13: mise en service
!         17-06-13: Debug nom de fichier mangrove
!         28-02-17  2D3D: possibilite de choix de la Methode Eq. 7
!         12-07-18  - mangrove_scheme permet de choisir entre Eq 7 et Eq 9
!                   - methode Eq 7 modifiee pour garantir propriete mang3dto2d=0
!                     si option "lineraire" choisie
!...............................................................................
!    _________                    .__                  .__             !(°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
      use module_principal
      use module_parallele
      implicit none

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine mangrove_init_mask
      implicit none

       allocate(r_mangrove  (1:imax  ,1:jmax,kmax)) ; r_mangrove=0.
       allocate(rmangrovebar(0:imax+1,0:jmax+1)   ) ; rmangrovebar=0. 
       allocate(mang3dto2d_u(0:imax+1,0:jmax+1)   ) ; mang3dto2d_u=0.
       allocate(mang3dto2d_v(0:imax+1,0:jmax+1)   ) ; mang3dto2d_v=0

! lecture du mask_mangrove global
      allocate(glob_mask_mangrove (0:iglb+1,0:jglb+1))
      glob_mask_mangrove(:,:)=0

      if(par%rank==0)write(6,'(a,a)')'sur le point de lire: ',mangrove_file_name
      open(unit=3,file=mangrove_file_name,recl=10000)
      do i=1,iglb
        read(3,'(1000i1)')(glob_mask_mangrove  (i,j),j=1,jglb)
      enddo
      close(3)

! allocation des mask_mangrove_t
      allocate(mask_mangrove_t (1:imax,1:jmax))

! ecriture des mask_mangrove_t

      do j=1,jmax
      do i=1,imax
          mask_mangrove_t(i,j)=                                        &
                      glob_mask_mangrove(i+par%timax(1),j+par%tjmax(1))
      enddo
      enddo

      deallocate(glob_mask_mangrove)

      end subroutine mangrove_init_mask

!......................................................................

      subroutine mangrove_2d_friction
      implicit none

! Action de la mangrove sur le mode externe
! Mangrove effect on momentum external mode as a linear friction term - implicit scheme -
! Mangrove effect on momentum external mode as a quadratic friction term - implicit scheme -

! the equation is now d(ubar)/dt = (ubar(t+1)-ubar(t-1))/dte = A - (rbar)*(ubar(t+1))
! with A = -1/H*SUM[ ((r-rbar)/(1+2*r*dti)*u(t-1)*dz ]

! We have to take into acount the first step
!      if (iteration3d==0) call mangrove_init_coef


! u component:
      do j=2,jmax-1
      do i=2,imax
          velbar_u(i,j,2)=(velbar_u(i,j,2)+mang3dto2d_u(i,j)*dte_lp)/  &
               (1.+0.5*(rmangrovebar(i,j)+rmangrovebar(i-1,j)*dte_lp))
      enddo
      enddo

! v component:
      do j=2,jmax
      do i=2,imax-1
          velbar_v(i,j,2)=(velbar_v(i,j,2)+mang3dto2d_v(i,j)*dte_lp)/  &
               (1.+0.5*(rmangrovebar(i,j)+rmangrovebar(i,j-1)*dte_lp))
      enddo
      enddo

      end subroutine mangrove_2d_friction

!----------------------------------------------------------------------

      subroutine mangrove_3d_friction
      implicit none

      double precision :: r_u, r_v
!      double precision,allocatable,dimension(:,:,:) :: r_mangrove

! Action de la mangrove sur le mode interne
! Mangrove effect on momentum internal mode as a linear friction term - implicit scheme -
! Mangrove effect on momentum internal mode as a quatratic friction term - implicit scheme -


! u component:
! calcul of the zonal velocity with mangrove dissipation correction 
! and verticaly averaging of the zonal dissipative coefficient for the external mode correction 

! rmangrovebar = 1/H*SUM [ r*dz ] = 1/H SUM[ cd(i,j,k)|u|(i,j,k)*dz ]
! mang3dto2d_u = -1/H*SUM[ (r-rbar_u)/(1+2*r*dt)*u*dz ]
! mang3dto2d_v = -1/H*SUM[ (r-rbar_v)/(1+2*r*dt)*v*dz ]
! where rbar_u=0.5*(rmangrovebar(i)+rmangrovebar(i-1))
! and   rbar_v=0.5*(rmangrovebar(j)+rmangrovebar(j-1))

! Evaluation de rmangrovebar = 1/H*SUM [ r*dz ] = 1/H SUM[ cd(i,j,k)|u|(i,j,k)*dz ]
!      allocate(r_mangrove(1:imax,1:jmax,1:kmax)) ; r_mangrove(:,:,:)=0.

! if linear friction linear_coef_mangrove=1
!      linear_coef_mangrove=0
      
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
!         r_u=0.5*(r_mangrove(i,j,k)+r_mangrove(i-1,j,k))
!         vel_u(i,j,k,2)=vel_u(i,j,k,2)/(1.+r_u*dti_lp )
          vel_u(i,j,k,2)=vel_u(i,j,k,2)/(1.+0.5*(r_mangrove(i,j,k)+r_mangrove(i-1,j,k))*dti_lp)
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
!         r_v=0.5*(r_mangrove(i,j,k)+r_mangrove(i,j-1,k))
!         vel_v(i,j,k,2)=vel_v(i,j,k,2)/(1.+r_v*dti_lp )
          vel_v(i,j,k,2)=vel_v(i,j,k,2)/(1.+0.5*(r_mangrove(i,j,k)+r_mangrove(i,j-1,k))*dti_lp)
      enddo
      enddo
      enddo


      end subroutine mangrove_3d_friction

!----------------------------------------------------------------------

      subroutine mang3dto2d
      implicit none

      double precision :: r_u, r_v

! initialisation of vertically averaged friction coefficient and frozen term for mangrove
! Mangrove effect  as a quatratic friction term


! verticaly averaging of the zonal dissipative coefficient for the external mode correction 
! rmangrovebar = 1/H*SUM [ r*dz ] = 1/H SUM[ cd(i,j,k)|u|(i,j,k)*dz ]
! mang3dto2d_u = -1/H*SUM[ (r-rbar_u)/(1+2*r*dt)*u*dz ]
! mang3dto2d_v = -1/H*SUM[ (r-rbar-v)/(1+2*r*dt)*v*dz ]
! where rbar_u=0.5*(rmangrovebar(i)+rmangrovebar(i-1))
! and   rbar_v=0.5*(rmangrovebar(j)+rmangrovebar(j-1))

! Evaluation de rmangrovebar = 1/H*SUM [ r*dz ] = 1/H SUM[ cd(i,j,k)|u|(i,j,k)*dz ]
!      allocate(r_mangrove(1:imax,1:jmax,1:kmax)) ; r_mangrove(:,:,:)=0.

!      print*,'linear_coef_mangrove :',linear_coef_mangrove

      rmangrovebar(:,:)=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
          r_mangrove(i,j,k)=coef_diss_mangrove*mask_mangrove_t(i,j)*(   &
                       linear_coef_mangrove + (1-linear_coef_mangrove)* &
                       sqrt(0.5*(vel_u(i,j,k,1)**2+vel_u(i+1,j,k,1)**2  &
                                +vel_v(i,j,k,1)**2+vel_v(i,j+1,k,1)**2))&
                                                                    )
          rmangrovebar(i,j)=   &
          rmangrovebar(i,j)    &
           +r_mangrove(i,j,k)  &
                 *dz_t(i,j,k,1)     !Sum( r.dz )
      enddo
      enddo
      enddo
      do j=1,jmax
      do i=1,imax
          rmangrovebar(i,j)=   &
          rmangrovebar(i,j)    &
                 /hz_w(i,j,1)       ! Sum( r.dz ) / H
      enddo
      enddo


      if(mangrove_scheme==1) then !-Methode Equation 7->

! Methode Equation 7 du document: !28-02-17
! https://docs.google.com/document/d/1V9cnPdDir-x3Ic3DjFHDpJ8osOkKobqs4pt5z4iHTPU/edit
! avec courant 3D pris A la moyenne de l'echeance before et now pour filtrer le mode numrique

! La modification apportee apres 12-07-18 a pour objectif de garantir la
! propriete mang3dto2d=0 si option lineaire choisie (c.a.d. si r=rbar)

! composante u:
      mang3dto2d_u(:,:)=0.
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax

! Avant 12-07-18
!         mang3dto2d_u(i,j)=                          & ! Sum( -dz*r.u )
!         mang3dto2d_u(i,j)                           &
!       -0.5*(r_mangrove(i,j,k)+r_mangrove(i-1,j,k))  &
!            *0.5*(vel_u(i,j,k,0)+vel_u(i,j,k,1))     & ! courant moyenne tems before et now
!                  *dz_u(i,j,k,1)                      

! Apres 12-07-18
          mang3dto2d_u(i,j)=                               & ! Sum( -dz*r.u )
          mang3dto2d_u(i,j)                                &
        -0.5*(   r_mangrove(i,j,k)+  r_mangrove(i-1,j,k)   &
              -rmangrovebar(i,j)  -rmangrovebar(i-1,j)  )  &
             *0.5*(vel_u(i,j,k,0)+vel_u(i,j,k,1))          & ! courant moyenne tems before et now
                   *dz_u(i,j,k,1)                      

      enddo
      enddo
      enddo
      do j=2,jmax-1
      do i=2,imax
! Avant 12-07-18
!        mang3dto2d_u(i,j)=   &  
!        mang3dto2d_u(i,j)    &
!               /hz_u(i,j,1)  &                                    ! Sum(-dz*r*u)/H
!     +0.5*(rmangrovebar(i,j)+rmangrovebar(i-1,j))*velbar_u(i,j,1) ! +rbar.ubar(t0)

! Apres 12-07-18
         mang3dto2d_u(i,j)=   &  
         mang3dto2d_u(i,j)    &
                /hz_u(i,j,1)  
      enddo
      enddo

! composante v:
      mang3dto2d_v(:,:)=0.
      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1

! Avant 12-07-18
!         mang3dto2d_v(i,j)=                          & ! Sum( -dz*r.v )
!         mang3dto2d_v(i,j)                           &
!       -0.5*(r_mangrove(i,j,k)+r_mangrove(i,j-1,k))  &
!            *0.5*(vel_v(i,j,k,0)+vel_v(i,j,k,1))     & ! courant moyenne tems before et now
!                  *dz_v(i,j,k,1)                      

! Apres 12-07-18
          mang3dto2d_v(i,j)=                               & ! Sum( -dz*r.v )
          mang3dto2d_v(i,j)                                &
        -0.5*(   r_mangrove(i,j,k)+  r_mangrove(i,j-1,k)   &
              -rmangrovebar(i,j)  -rmangrovebar(i,j-1  ))  &
             *0.5*(vel_v(i,j,k,0)+vel_v(i,j,k,1))     & ! courant moyenne tems before et now
                   *dz_v(i,j,k,1)                      

      enddo
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
! Avant 12-07-18
!        mang3dto2d_v(i,j)=   &  
!        mang3dto2d_v(i,j)    &
!               /hz_v(i,j,1)  &                                    ! Sum(-dz*r*v)/H
!     +0.5*(rmangrovebar(i,j)+rmangrovebar(i,j-1))*velbar_v(i,j,1) ! +rbar.vbar(t0)

! Apres 12-07-18
         mang3dto2d_v(i,j)=   &  
         mang3dto2d_v(i,j)    &
                /hz_v(i,j,1)  
      enddo
      enddo

      else                        !-Methode Equation 9->

! Methode Equation 9 du document:
! https://docs.google.com/document/d/1V9cnPdDir-x3Ic3DjFHDpJ8osOkKobqs4pt5z4iHTPU/edit
      mang3dto2d_u(:,:)=0.
      mang3dto2d_v(:,:)=0.
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
          r_u=0.5*(r_mangrove(i,j,k)+r_mangrove(i-1,j,k))
          mang3dto2d_u(i,j) = mang3dto2d_u(i,j)                       &
                             -(r_u-0.5*(rmangrovebar(i,j)             &
                                       +rmangrovebar(i-1,j)) )        &
                              *vel_u(i,j,k,1)*dz_u(i,j,k,2)           &
                              /(1.+r_u*dti_lp)/hz_u(i,j,2)

      enddo
      enddo
      enddo

      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
          r_v=0.5*(r_mangrove(i,j,k)+r_mangrove(i,j-1,k))
          mang3dto2d_v(i,j) = mang3dto2d_v(i,j)                       &
                             -(r_v-0.5*(rmangrovebar(i,j)             &
                                       +rmangrovebar(i,j-1)) )        &
                              *vel_v(i,j,k,1)*dz_v(i,j,k,2)           &
                              /(1.+r_v*dti_lp)/hz_v(i,j,2)

      enddo
      enddo
      enddo

      endif                       !-Methode Equation 9->

      end subroutine mang3dto2d

!----------------------------------------------------------------------

      end module module_mangrove
