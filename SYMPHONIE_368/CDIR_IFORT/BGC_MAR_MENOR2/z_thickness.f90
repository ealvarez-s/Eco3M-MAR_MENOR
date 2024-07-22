










      subroutine z_thickness(t1_,t2_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 296 - last update: 28-02-21
!______________________________________________________________________
      use module_principal
      use module_parallele
      implicit none
      integer t1_,t2_

!......................................................................
! modifs 01/11/01: inversion de l'ordre des boucles: J passe avant I
!        08/06/07: bienvenue à hzdy_x et hzdx_y
!        06/04/09: bornes, obc, parallelisation
!        06-06-09: boucle de calcul de hz_z etendue à i=0 & j=0 pour
!                  parallelisation. Suppression de call obc_hz devenu
!                  inutile
!        14-03-14  Forcer hz_w a rester strictement positif
!        28-11-14  suppression hzdy_u hzdx_v
!        10-01-16  amenagement commentaires
! v296   28-02-21  seuiller moins petit en zone intertidale
!...............................................................................
!    _________                    .__                  .__                     ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................
! Computes the thickness of the water column

      do k=t1_,t2_

       do j=0,jmax+1 ; do i=0,imax+1
        hz_w(i,j,k)=max(      wetdry_cst3,h_w(i,j)+ssh_int_w(i,j,k)) !28-02-21
!       hz_w(i,j,k)=max(0.001*wetdry_cst3,h_w(i,j)+ssh_int_w(i,j,k))
!       hz_w(i,j,k)=h_w(i,j)+ssh_int_w(i,j,k)
       enddo         ; enddo

       do j=0,jmax+1 ; do i=1,imax+1
        hz_u(i,j,k)=max(      wetdry_cst3,h_u(i,j)+0.5*(ssh_int_w(i,j,k)+ssh_int_w(i-1,j,k))) !28-02-21
!       hz_u(i,j,k)=max(0.001*wetdry_cst3,h_u(i,j)+0.5*(ssh_int_w(i,j,k)+ssh_int_w(i-1,j,k)))
!       hz_u(i,j,k)=h_u(i,j)+0.5*(ssh_int_w(i,j,k)+ssh_int_w(i-1,j,k))
       enddo         ; enddo

       do j=1,jmax+1 ; do i=0,imax+1
        hz_v(i,j,k)=max(      wetdry_cst3,h_v(i,j)+0.5*(ssh_int_w(i,j,k)+ssh_int_w(i,j-1,k)))
!       hz_v(i,j,k)=max(0.001*wetdry_cst3,h_v(i,j)+0.5*(ssh_int_w(i,j,k)+ssh_int_w(i,j-1,k)))
!       hz_v(i,j,k)=h_v(i,j)+0.5*(ssh_int_w(i,j,k)+ssh_int_w(i,j-1,k))
       enddo         ; enddo

      enddo ! fin de boucle sur K

      end subroutine z_thickness

!........................................................................

      subroutine z_thickness_changebottomdepth
      use module_principal
      use module_parallele
      implicit none

! https://docs.google.com/document/d/1qxWX2Uf1KKAflyokBXaOqacVvez086WwzLCwesaV45U/edit#

    


!************************************************************
! Ici introduire l'increment de bathymetry
! Ici donner un increment de bathymetry (tableau xy_t(:,:,1))
      i0=iglb/2 ; j0=jglb/2
      x0=(1./86400.)*dti_fw ! 1m par jour....
      do j=0,jmax+1 ; do i=0,imax+1
        i1=i+par%timax(1)
        j1=j+par%tjmax(1)
        dist1=real(i1-i0)**2+real(j1-j0)**2
        xy_t(i,j,1)=x0*exp(-0.001*dist1)
      enddo ; enddo 
! Ici introduire l'increment de bathymetry
!************************************************************

! Apporter l'increment A h_w, h_u, h_v
      do j=0,jmax+1 ; do i=0,imax+1
       h_w(i,j)=h_w(i,j)+xy_t(i,j,1)
       if(i+par%timax(1)==iglb/2.and.j+par%tjmax(1)==jglb/2)write(10+par%rank,*)elapsedtime_now,h_w(i,j)
      enddo ; enddo 
      do j=0,jmax+1 ; do i=1,imax+1
       h_u(i,j)=0.5*(h_w(i,j)+h_w(i-1,j))
      enddo ; enddo 
      do j=0,jmax+1 ; do i=1,imax+1
       h_v(i,j)=0.5*(h_w(i,j)+h_w(i,j-1))
      enddo ; enddo 

! Apporter le mEme increment A la SSH pour ne pas changer le volume d'eau
      if(flag_timesplitting==1)  then  !-time-splitting-case-> 
       do j=0,jmax+1 ; do i=0,imax+1
        ssh_w(i,j,1)=ssh_w(i,j,1)+xy_t(i,j,1)
        ssh_w(i,j,2)=ssh_w(i,j,2)+xy_t(i,j,1)
       enddo ; enddo 
      endif                           !-time-splitting-case->

      if(flag_timesplitting==0)  then !-NO-time-splitting-case-> 
       do j=0,jmax+1 ; do i=0,imax+1
        ssh_int_w(i,j,1)=ssh_int_w(i,j,1)+xy_t(i,j,1)
       enddo ; enddo 
      endif                           !-NO-time-splitting-case->

      end subroutine z_thickness_changebottomdepth
