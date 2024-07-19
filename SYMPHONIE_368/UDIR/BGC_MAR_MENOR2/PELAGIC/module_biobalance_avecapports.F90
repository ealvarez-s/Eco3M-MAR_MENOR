      module module_biobalance
      use module_principal
      use module_parallele
      use module_biology
      use ModuleDeclaration
      implicit none
      include 'netcdf.inc'

contains


!...........................................

      subroutine biobalance_gateA
      implicit none


      end subroutine biobalance_gateA

!...........................................

      subroutine biobalance_gateB
      implicit none


      end subroutine biobalance_gateB

!...........................................

      subroutine biobalance_gateC
      implicit none

      double precision fluxNemed,fluxPemed,fluxSiemed, &
                       fluxNwmed,fluxPwmed,fluxSiwmed,fluxAMwmed

#ifdef benthic
! budget at the sea/sediment interface 
     call budget_benthic


! flux at the vertical boundaries    
      call BCFlux
#endif

! flux of oxygen
      fluxbio_w(:,:,:,2)=0.   !25-12-14   Mise à 0 car on va sommer oxygene et apports riviere
      call oxygen_surfaceflux
! SEUIL MAX
!      fluxNwmed  = 72825.e9/365/24/3600/sumareawmed ! 10042  !  72825.e9/365/24/3600/7.2281e11 !mmolN m-2 s-1
!      fluxPwmed  = 697.e9/365/24/3600/sumareawmed !  357 ! 697.e9/365/24/3600/7.2281e11 !mmolP m-2 s-1 
!      fluxSiwmed = 2201.e9/365/24/3600/sumareawmed !  709 !2201.e9/365/24/3600/7.2281e11 !mmolSi m-2 s-1
!      fluxNemed  = 73621.e9/365/24/3600/sumareaemed ! 6064 !73621.e9/365/24/3600/1.2775e12 !mmolN m-2 s-1
!      fluxPemed  = 957.e9/365/24/3600/sumareaemed !  379  !957.e9/365/24/3600/1.2775e12 !mmolP m-2 s-1
!      fluxSiemed = 5153.e9/365/24/3600/sumareaemed !  700 ! 5153.e9/365/24/3600/1.2775e12 !mmolSi m-2 s-1
!
!      fluxAMwmed  = 19904.e9/365/24/3600/sumareawmed  ! MOOSE
! SEUIL MIN
      fluxNwmed  = 10042.e9/365.0/24.0/3600.0/sumareawmed ! RIBERA 2003
      fluxPwmed  = 357.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
      fluxSiwmed = 709.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
      fluxNemed  = 6064.e9/365.0/24.0/3600.0/sumareaemed ! RIBERA 2003
      fluxPemed  = 379.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
      fluxSiemed = 700.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003

      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE



!      do vb=1,vbmax
      do i=1,imax
      do j=1,jmax

! bassin ouest
!        if  (lon_t(i,j)*rad2deg<10.or.    &
!        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42.or. &
!        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25.or. &
!        lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25)THEN


      IF(mask_t(i,j,kmax+1).EQ.1.and.lon_t(i,j)*rad2deg>-5.6)then !debut test

      if ((lon_t(i,j)*rad2deg<10).or.    &
        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42).or. &
        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25).or. &
        (lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25))then


             fluxbio_w(i,j,iNitrate,2)  = fluxNwmed
             fluxbio_w(i,j,iPhosphate,2)= fluxPwmed
             fluxbio_w(i,j,iSilice,2)   = fluxSiwmed

             fluxbio_w(i,j,iAmmonium,2)   = fluxAMwmed

             fluxbio_w(i,j,ismopn,2)   = fluxNwmed*.3
!             fluxbio_w(i,j,ilmopn,2)   = fluxNwmed*.3
             fluxbio_w(i,j,imodn,2)   = fluxNwmed*.4

             fluxbio_w(i,j,ismopp,2)   = fluxPwmed!*3
!             fluxbio_w(i,j,ilmopp,2)   = fluxPwmed*3
             fluxbio_w(i,j,imodp,2)   = fluxPwmed

             fluxbio_w(i,j,ismopsi,2)   = fluxSiwmed*.3
!             fluxbio_w(i,j,ilmopsi,2)   = fluxSiwmed*.3

!bassin est
         else
             fluxbio_w(i,j,iNitrate,2)  = fluxNemed
             fluxbio_w(i,j,iPhosphate,2)= fluxPemed
             fluxbio_w(i,j,iSilice,2)   = fluxSiemed

             fluxbio_w(i,j,iAmmonium,2)   = fluxAMwmed

             fluxbio_w(i,j,ismopn,2)   = fluxNemed*.3
!             fluxbio_w(i,j,ilmopn,2)   = fluxNemed*.3
             fluxbio_w(i,j,imodn,2)   = fluxNemed*.4

             fluxbio_w(i,j,ismopp,2)   = fluxPemed!*3
!             fluxbio_w(i,j,ilmopp,2)   = fluxPemed*3
             fluxbio_w(i,j,imodp,2)   = fluxPemed

             fluxbio_w(i,j,ismopsi,2)   = fluxSiemed*.3
!             fluxbio_w(i,j,ilmopsi,2)   = fluxSiemed*.3

         endif
        ENDIF
      enddo
      enddo
!      enddo









! flux des rivieres en condition limite de surface turbulence dans le cas nemoffline
      if(flag_nemoffline==1) then !wwwwwwwwwwwwwwwww>       !25-12-14
         do kr=1,nriver
          if(rivertrc_inout(kr)==1) then
           i=iriver(kr,1)
           j=jriver(kr,1)
           do vb=1,vbmax
            fluxbio_w(i,j,vb,2)=fluxbio_w(i,j,vb,2)+ &
             runoff_w(kr,1)*river_bio(vb,kr,1)   
           enddo
          endif
         enddo
      endif                       !wwwwwwwwwwwwwwwww>

      call Budget_River

        if(iteration3d/=0)then
        x2=elapsedtime_now/86400.
        x1=elapsedtime_bef/86400.
        j2=int(x2)
        j1=int(x1)
        if(j1/=j2)call sortie_sedimentation
        endif

      end subroutine biobalance_gateC


!...........................................

      subroutine biobalance_gateD
      implicit none
      call budget_export_bottom

      end subroutine biobalance_gateD

!...........................................

      subroutine biobalance_gateE
      implicit none


      end subroutine biobalance_gateE

!...........................................

      subroutine biobalance_gateF
      implicit none


      end subroutine biobalance_gateF

!...........................................

      end module module_biobalance
