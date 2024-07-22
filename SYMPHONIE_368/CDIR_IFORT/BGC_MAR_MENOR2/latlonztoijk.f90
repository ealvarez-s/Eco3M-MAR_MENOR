










      subroutine latlonztoijk(lon_,lat_,depth_,txt_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 256 - last update: 05-06-19
!______________________________________________________________________
!__________________________________________________________________________
! Version date      Description des modifications
! 2010.13 03-11-10  Mise en service
! v256    05-06-19  forcer le passage en indices "locaux"
!__________________________________________________________________________
!...............................................................................
!    _________                    .__                  .__             ! 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal
      use module_parallele !#mpi
      implicit none
      double precision lon_,lat_,depth_
      character*3 txt_

! This routine returns deci decj decki, the (real) grid indexes corresponding to
! lon_, lat_ (longitude and latitude in radians) and depth_ (depth in m).


      flag_stop=0
      if(txt_/='loc') then !>>> !05-06-19
        flag_stop=1
!        write(10+par%rank,*) &
!        'Use latlonztoijk subroutine with "loc" mode'
      endif                !>>>
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop 'Err latlonztoijk see fort.xxx files'

      call latlontoij(lon_,lat_,txt_)

      i=int(deci) ; j=int(decj) ; rapi=deci-i  ; rapj=decj-j

      deck=1.

      if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !m°v°m> !05-06-19

      x1=(1.-rapi)*(1.-rapj)*depth_t(i  ,j  ,1)                &
        +(1.-rapi)*    rapj *depth_t(i  ,j+1,1)                &
        +    rapi *(1.-rapj)*depth_t(i+1,j  ,1)                &
        +    rapi *    rapj *depth_t(i+1,j+1,1)

      do k=1,kmax-1

       x2=(1.-rapi)*(1.-rapj)*depth_t(i  ,j  ,k+1)                &
         +(1.-rapi)*    rapj *depth_t(i  ,j+1,k+1)                &
         +    rapi *(1.-rapj)*depth_t(i+1,j  ,k+1)                &
         +    rapi *    rapj *depth_t(i+1,j+1,k+1)

       if(depth_>=x1.and.depth_<=x2) then  !>>>>>>>

        deck=k+(depth_-x1)/(x2-x1)

        goto 100

       endif                                     !>>>>>>>

       x1=x2

      enddo

      if(depth_>x2)deck=kmax

  100 continue


      endif                                              !m°v°m> !05-06-19

      end subroutine latlonztoijk
