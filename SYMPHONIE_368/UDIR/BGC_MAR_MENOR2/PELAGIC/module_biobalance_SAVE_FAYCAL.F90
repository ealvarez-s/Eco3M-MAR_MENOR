      module module_biobalance
      use module_principal
      use module_parallele
      use ModuleComBlock
      use UserDeclaration
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
 !   call budget_benthic


! flux at the vertical boundaries    
      call BCFlux
#endif

! flux of oxygen
      fluxbio_w(:,:,:,2)=0.   !25-12-14   Mise à 0 car on va sommer oxygene et apports riviere
      call oxygen_surfaceflux

!flux atmospheriques
      fluxNwmed  = 10042.e9/365.0/24.0/3600.0/sumareawmed ! RIBERA 2003
      fluxPwmed  = 357.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
      fluxSiwmed = 709.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
      fluxNemed  = 6064.e9/365.0/24.0/3600.0/sumareaemed ! RIBERA 2003
      fluxPemed  = 379.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
      fluxSiemed = 700.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE

!!! Modif Alex 25/01/2018 On fait un test sur les dépôts atmosph
!      fluxNwmed  = 40000.e9/365.0/24.0/3600.0/sumareawmed ! MID
!      fluxPwmed  = 550.e9/365.0/24.0/3600.0/sumareawmed !  MID
!      fluxSiwmed = 709.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxNemed  = 23000.e9/365.0/24.0/3600.0/sumareaemed ! LOW MID
!      fluxPemed  = 670.e9/365.0/24.0/3600.0/sumareaemed !  MID
!      fluxSiemed = 700.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE

!! Modif Alex 01/12 On passe à l'échelle Mid de Ribera 2003
!      fluxNwmed  = 41434.e9/365.0/24.0/3600.0/sumareawmed ! RIBERA 2003
!      fluxPwmed  = 527.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxSiwmed = 1455.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxNemed  = 39843.e9/365.0/24.0/3600.0/sumareaemed ! RIBERA 2003
!      fluxPemed  = 668.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxSiemed = 2581.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE

!! Modif Alex 01/12 On passe à l'échelle High de Ribera 2003
!      fluxNwmed  = 72825.e9/365.0/24.0/3600.0/sumareawmed ! RIBERA 2003
!      fluxPwmed  = 697.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxSiwmed = 2201.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxNemed  = 73621.e9/365.0/24.0/3600.0/sumareaemed ! RIBERA 2003
!      fluxPemed  = 957.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxSiemed = 5153.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE

       do i=1,imax
      do j=1,jmax

      IF(mask_t(i,j,kmax+1).EQ.1.and.lon_t(i,j)*rad2deg>-5.6)then !debut test

! bassin ouest
        if  ((lon_t(i,j)*rad2deg<10.and.lon_t(i,j)*rad2deg>-5.6).or.    &
        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42.or. &
        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25.or. &
        lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25) then


             fluxbio_w(i,j,iNitrate,2)  = fluxNwmed
             fluxbio_w(i,j,iPhosphate,2)= fluxPwmed
             fluxbio_w(i,j,iSilice,2)   = fluxSiwmed

             fluxbio_w(i,j,iAmmonium,2)   = fluxAMwmed
             fluxbio_w(i,j,ismopn,2)   = fluxNwmed*.3
             fluxbio_w(i,j,imodn,2)   = fluxNwmed*.4
             fluxbio_w(i,j,ismopp,2)   = fluxPwmed!*3
             fluxbio_w(i,j,imodp,2)   = fluxPwmed
             fluxbio_w(i,j,ismopsi,2)   = fluxSiwmed*.3

!bassin est
         else
!           if(lon_t(i,j)*rad2deg>-5.6) then
             fluxbio_w(i,j,iNitrate,2)  = fluxNemed
             fluxbio_w(i,j,iPhosphate,2)= fluxPemed
             fluxbio_w(i,j,iSilice,2)   = fluxSiemed

             fluxbio_w(i,j,iAmmonium,2)   = fluxAMwmed
             fluxbio_w(i,j,ismopn,2)   = fluxNemed*.3
             fluxbio_w(i,j,imodn,2)   = fluxNemed*.4
             fluxbio_w(i,j,ismopp,2)   = fluxPemed!*3
             fluxbio_w(i,j,imodp,2)   = fluxPemed
             fluxbio_w(i,j,ismopsi,2)   = fluxSiemed*.3
!           endif

         endif

         endif
      enddo
      enddo

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
