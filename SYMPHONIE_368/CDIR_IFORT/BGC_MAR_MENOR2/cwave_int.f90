










      subroutine cwave_int
!______________________________________________________________________
! SYMPHONIE ocean model
! last update: 25-03-19
!______________________________________________________________________
!...............................................................................
! modifs 15/12/02: mise en service
!        10/05/06: fonctions compatibles avec double precision
!        22/01/09: critere CFL tient compte du nouveau systeme horizontal
! v250   19-03-19  ajout subroutine cwave_int_shortwaves 
!        21-03-19  suite ajout subroutine cwave_int_shortwaves c adimensionnee 
!                  par dt/dx
!        25-03-19  x1=max(h_w(i1,j),wetdry_cst1) 
!...............................................................................
!    _________                    .__                  .__             !
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !m°v°m
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


      use module_principal
      implicit none

!...............................................................................
! Description:
! Calcul de la vitesse des ondes internes
!...............................................................................

      do loop1=1,2 ! debut boucle loop1

       const3=1. ! vitesse de phase en m/s

       if(loop1.eq.1) then
        i1=1
        j1=1
        i2=1
        j2=1
       else
        i1=imax
        j1=jmax
        i2=imax+1
        j2=jmax+1
       endif

       do i=2,imax
        cwi_int_u(i,loop1)=min(const3,dy_u(i,j1)/dti_fw)
        cwi_int_v(i,loop1)=min(const3,dy_v(i,j2)/dti_fw)
       enddo
       do j=2,jmax
        cwj_int_v(j,loop1)=min(const3,dx_v(i1,j)/dti_fw)
        cwj_int_u(j,loop1)=min(const3,dx_u(i2,j)/dti_fw)
       enddo

      enddo ! fin boucle loop1

      end subroutine cwave_int

!...............................................................

      subroutine cwave_int_shortwaves !19-03-19
      use module_principal ; use module_wave ; use module_parallele
      implicit none
      integer loop_

! ATTENTION DANS CETTE VERSION LA VITESSE DE PHASE EST ADIMENSIONNEE PAR DT/DX

      do loop_=1,2 ! debut boucle loop_


       if(loop_.eq.1) then
        i1=1
        j1=1
        i2=1
        j2=1
       else
        i1=imax
        j1=jmax
        i2=imax+1
        j2=jmax+1
       endif

       do i=1,imax
        x1=max(h_w(i,j1),wetdry_cst1) !25-03-19
        x2=periodpeak
        call wave_freq2kvector(x1,x2,x3) ! x3 kvector
        xy_t(i,j1,1)=2.*pi/(periodpeak*x3) !c=2*pi/(period*kvector)
       enddo
       do i=2,imax
!       cwi_int_u(i,loop_)=min(0.5*(xy_t(i-1,j1,1)+xy_t(i,j1,1)),dy_u(i,j1)/dti_fw)
!       cwi_int_v(i,loop_)=min(                    xy_t(i,j1,1) ,dy_v(i,j2)/dti_fw)
        cwi_int_u(i,loop_)=0.5*(xy_t(i-1,j1,1)+xy_t(i,j1,1))*dti_fw/dy_u(i,j1) !21-03-19
        cwi_int_v(i,loop_)=                    xy_t(i,j1,1) *dti_fw/dy_v(i,j2)
       enddo

       do j=1,jmax
        x1=max(h_w(i1,j),wetdry_cst1) !25-03-19
        x2=periodpeak
        call wave_freq2kvector(x1,x2,x3) ! x3 kvector
        xy_t(i1,j,1)=2.*pi/(periodpeak*x3) !c=2*pi/(period*kvector)
       enddo
       do j=2,jmax
        cwj_int_v(j,loop_)=0.5*(xy_t(i1,j-1,1)+xy_t(i1,j,1))*dti_fw/dx_v(i1,j)
        cwj_int_u(j,loop_)=                    xy_t(i1,j,1) *dti_fw/dx_u(i2,j)
       enddo

      enddo ! fin boucle loop_

      end subroutine cwave_int_shortwaves

!...............................................................
