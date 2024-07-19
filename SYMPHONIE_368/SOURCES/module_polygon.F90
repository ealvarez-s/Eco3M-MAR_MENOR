      module module_polygon
!______________________________________________________________________
! SYMHPHONE ocean model
! release 363 - last update: 15-01-23
!______________________________________________________________________
      use module_principal
      implicit none
      integer irep
!...............................................................................
! Version date    Description des modifications
! v363  15-01-23  Mise en service
!...............................................................................
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      contains 

!...................................
#ifdef bidon
      subroutine polygon_driver(name_poly,type_poly)
      implicit none
      integer type_poly,ichoix
      character(len=*) name_poly
      double precision xlon,xlat

       ichoix=0 ; xlon=0. ; xlat=0.

       write(6,'(a,a)')'Polygon file ',trim(name_poly)
       call polygon_in(ichoix,irep,xlon,xlat,name_poly)
       ichoix=1
       do i=1,imax
       do j=1,jmax
        call polygon_in(ichoix,irep,lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg,name_poly)
!       if (type_poly==0) then
!         if(irep==1)mask_t(i,j)=0
!       else
!         if(irep==1)mask_t(i,j)=1
!       endif
!       if(irep==1)mask_t(i,j)=type_poly ! 1 ou 0
       enddo
       enddo
      end subroutine polygon_driver
#endif
!...................................

      subroutine polygon_in(ichoix,irep,x,y,name_poly)
      implicit none
      integer npoints_poly,irep,ichoix,i,j
      logical ldinmesh
      double precision x,y
      character(len=*) name_poly
      real ,dimension(:)  ,allocatable ::   xp,yp
      save xp,yp,npoints_poly

      irep=0

!----------------------------------------
! Lecture, allocation
      if(ichoix==0)then !>>>
       open(unit=12,file=trim(name_poly))
       read(12,*)npoints_poly
       if(allocated(xp)) deallocate (xp)
       if(allocated(yp)) deallocate (yp)
       allocate(xp(npoints_poly))
       allocate(yp(npoints_poly))
       do i=1,npoints_poly
       read(12,*)xp(i),yp(i)
       enddo
       return
      endif             !>>>

!----------------------------------------
! Algo: dans ou hors du polygone
      if(ichoix==1)then !ooo>
        ldinmesh=.FALSE.
       j=npoints_poly
      do i=1,npoints_poly
       if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
        .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
       if(ldinmesh)then
        ldinmesh=.false.
       else
        ldinmesh=.true.
       endif
       endif
       j=i
      enddo
      endif              !oooo>
      if(ldinmesh)irep=1
      return
      end subroutine polygon_in

!...................................

      end module module_polygon
