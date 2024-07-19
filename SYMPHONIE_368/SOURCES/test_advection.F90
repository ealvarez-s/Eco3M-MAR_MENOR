      subroutine test_advection
!______________________________________________________________________
!
! S model
! release 2010.10  - last update: 24-06-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
      integer key
#ifdef synopsis
       subroutinetitle='test_advection'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! initialisation temperature
      x1=real(imax)/4.
      x2=real(jmax)/2.
      do k=1,kmax
      do i=0,imax+1
      do j=0,jmax+1
      dist0=sqrt( (real(i)-x1)**2 + (real(j)-x2)**2 )
      dist1=max(20.-dist0,0.)/20.

      tem_t(i,j,k,-1)=dist1+0.2*cos(pi*i)*cos(pi*j)*dist1
      tem_t(i,j,k,0 )=dist1+0.2*cos(pi*i)*cos(pi*j)*dist1
      tem_t(i,j,k,1 )=dist1+0.2*cos(pi*i)*cos(pi*j)*dist1
      tem_t(i,j,k,2 )=dist1+0.2*cos(pi*i)*cos(pi*j)*dist1

      enddo
      enddo
      enddo



! un tour en KEY
!     KEY=4*250 ! nb de courant =0.1
      key=4*150 ! nb de courant =0.3
!     DO 2000 KOUNT=0,KEY/4 ! quart de tour
      do 2000 kount=0,key   ! un tour
!     DO 2000 KOUNT=0,0    ! un coup
!     WRITE(6,*)KOUNT,KEY

!..............................................................
!___________________________________________________
! module du courant
!     X4=0.1 ! nb de courant = 0.1
      x4=0.3 ! nb de courant = 0.3
      x3=x4*dxb/dti_lp
! direction du vent
      x2=2.*pi*real(kount)/real(key)

      do k=1,kmax
       do i=0,imax+1
       do j=0,jmax+1

        vel_u(i,j,k,1)=x3*sin(x2)
        vel_v(i,j,k,1)=x3*cos(x2)

        veldydz_u(i,j,k,1)=vel_u(i,j,k,1)*dy_u(i,j)*dz_u(i,j,k,1)       !30-09-09
        veldxdz_v(i,j,k,1)=vel_v(i,j,k,1)*dx_v(i,j)*dz_v(i,j,k,1)       !30-09-09


       enddo
       enddo
      enddo

      if(mod(kount,20).eq.0)                                            &
       write(6,*)'%tour=',100.*real(kount)/real(key)
      if(mod(kount,20).eq.0)                                            &
       write(6,*)'vel_u vel_v module',vel_u(1,1,1,1),vel_v(1,1,1,1)     &
       ,sqrt(vel_u(1,1,1,1)**2+vel_v(1,1,1,1)**2)
!___________________________________________________
!..............................................................


      call advection_scal
      call moveforward



      sum1=0.
      do j=1,jmax
      do i=1,imax
      sum1=sum1+cos(pi*i)*cos(pi*j)*tem_t(i,j,kmax,1)
      enddo
      enddo
      write(6,*)'bruit',kount,key,sum1/real(imax)/real(jmax)
!     write(66,*)sum1/real(imax)/real(jmax)

 2000 continue

! initialisation temperature
      x1=real(imax)/4.
      x2=real(jmax)/2.
      do k=1,kmax
      do i=0,imax+1
      do j=0,jmax+1
      dist0=sqrt( (real(i)-x1)**2 + (real(j)-x2)**2 )
      dist1=max(20.-dist0,0.)/20.
      tem_t(i,j,k,0 )=dist1+0.2*cos(pi*i)*cos(pi*j)*dist1
      enddo
      enddo
      enddo

! CRITERE DE QUALITE:
      sum1=0.
      sum2=0.
      k=kmax
      do i=1,imax
      do j=1,jmax
      sum1=sum1+ tem_t(i,j,k,1)**2
      sum2=sum2+(tem_t(i,j,k,1)-tem_t(i,j,k,0))**2
      enddo
      enddo
      write(6,*)'critere qualite: 0=parfait, 1=diffusif, >>1=instable'
      write(6,*)'---------------> ',sqrt(sum2)/sqrt(sum1)

      call graph_out
      write(6,*)'dans sbr test_advection'
      write(6,*)'-----------------------'
      stop 'test_advection'


      return
      end
