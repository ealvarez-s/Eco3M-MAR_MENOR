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
                       fluxNwmed,fluxPwmed,fluxSiwmed,fluxAMwmed, &
                       fluxPwmedmin,fluxPemedmin,fluxPwmedp,fluxPemedp

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
      if(idust==1.or.idepnit==1.or.idepammo==1) call AtmDeposition


!      fluxNwmed  = 10042.e9/365.0/24.0/3600.0/sumareawmed ! RIBERA 2003
!      fluxPwmed  = 357.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxSiwmed = 709.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
!      fluxNemed  = 6064.e9/365.0/24.0/3600.0/sumareaemed ! RIBERA 2003
!      fluxPemed  = 379.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxSiemed = 700.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
!      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE

! Caroline: je prends pour OP 100% de DIP déduit des depots de dust
! d'ALADIN (% dans Powley et al. 2017, et données Frioul et Béar),
!  que je repartis de facon homogène. On pourrait aussi mettre
! une augmentation en fonction de la latitude.
      fluxPwmed  = 2.9e12/31/365.0/24.0/3600.0/sumareawmed 
      fluxPemed  = 8.4e12/31/365.0/24.0/3600.0/sumareaemed

!On pourrait aussi mettre une augmentation en fonction de la latitude.
      fluxPwmedmin  = 1.9e12/31/365.0/24.0/3600.0/sumareawmed 
      fluxPemedmin  = 5.4e12/31/365.0/24.0/3600.0/sumareaemed
      fluxPwmedp    = 1.e12/31/365.0/24.0/3600.0/sumareawmed
      fluxPemedp    = 1.e12/31/365.0/24.0/3600.0/sumareaemed


       do i=1,imax
      do j=1,jmax

      IF(mask_t(i,j,kmax+1).EQ.1.and.lon_t(i,j)*rad2deg>-5.6)then !debut test

! bassin ouest
        if  ((lon_t(i,j)*rad2deg<=10.and.lon_t(i,j)*rad2deg>-5.6).or.    &
        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<=15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<=42.or. &
        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25.or. &
        lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25) then


!             fluxbio_w(i,j,iNitrate,2)  = fluxNwmed
!             fluxbio_w(i,j,iPhosphate,2)= fluxPwmed
             fluxbio_w(i,j,iSilice,2)   = fluxSiwmed

!             fluxbio_w(i,j,iAmmonium,2)   = fluxAMwmed
!             fluxbio_w(i,j,ismopn,2)   = fluxNwmed*.3
!             fluxbio_w(i,j,imodn,2)   = fluxNwmed*.4
!             fluxbio_w(i,j,ismopp,2)   = fluxPwmed!*3
!             fluxbio_w(i,j,imodp,2)   = fluxPwmed
             fluxbio_w(i,j,ismopsi,2)   = fluxSiwmed*.3

! caroline 
! homogène :
!            fluxbio_w(i,j,ismopp,2)   = 0.
!            fluxbio_w(i,j,imodp,2)   = fluxPwmed
! variation en fonction de la latitude :
             fluxbio_w(i,j,ismopp,2)   = 0.
             fluxbio_w(i,j,imodp,2) = fluxPwmedmin+(lat_t(i,j)*rad2deg-34.5)/5*fluxPwmedp



!bassin est
         else
!           if(lon_t(i,j)*rad2deg>-5.6) then
!             fluxbio_w(i,j,iNitrate,2)  = fluxNemed
!            fluxbio_w(i,j,iPhosphate,2)= fluxPemed
             fluxbio_w(i,j,iSilice,2)   = fluxSiemed

!             fluxbio_w(i,j,iAmmonium,2)   = fluxAMwmed
!             fluxbio_w(i,j,ismopn,2)   = fluxNemed*.3
!             fluxbio_w(i,j,imodn,2)   = fluxNemed*.4
!             fluxbio_w(i,j,ismopp,2)   = fluxPemed!*3
!             fluxbio_w(i,j,imodp,2)   = fluxPemed
             fluxbio_w(i,j,ismopsi,2)   = fluxSiemed*.3
!           endif

! caroline 
! homogène :
!            fluxbio_w(i,j,ismopp,2)   = 0.
!            fluxbio_w(i,j,imodp,2)   = fluxPwmed
! variation en fonction de la latitude :
             fluxbio_w(i,j,ismopp,2)   = 0.
             fluxbio_w(i,j,imodp,2) = fluxPemedmin+(lat_t(i,j)*rad2deg-30)*5.5/10*fluxPemedp


         endif

!! MODIF ALEX: ON and IN from map in Kanakidou 2012
!         if(lat_t(i,j)*rad2deg>38) then
!
!             fluxbio_w(i,j,iNitrate,2)  = 0.00000030199801889299 !200/14/365/24/3600*2/3 ! 2/3 NO3-
!             fluxbio_w(i,j,iAmmonium,2)  = 0.00000015099900944649 !200/14/365/24/3600*1/3 ! 1/3 NH4+

! 23/02/2018 MODIF ALEX: IN from map in Kanakidou 2012 and Frioul + CapBear Obs
! for lat > 38 (see Aude Carreric report)
         if(lat_t(i,j)*rad2deg>38) then
             fluxbio_w(i,j,iNitrate,2)  = ((70*lat_t(i,j)*rad2deg-2460)/14/365/24/3600)*2/3 !2/3 NO3-
             fluxbio_w(i,j,iAmmonium,2)  = ((70*lat_t(i,j)*rad2deg-2460)/14/365/24/3600)*1/3 !1/3 NH4+

         else !  if(lat_t(i,j)*rad2deg<=38) then
             fluxbio_w(i,j,iNitrate,2)  = ((9.6277278562 * lat_t(i,j)*rad2deg-165.8536585366) &
                                          /14/365/24/3600) * 2/3 ! 2/3 NO3-
	     fluxbio_w(i,j,iAmmonium,2)  = ((9.6277278562 * lat_t(i,j)*rad2deg-165.8536585366) &
                                          /14/365/24/3600) * 1/3 ! 1/3 NH4+
         endif

! ON from map in Kanakidou 2012
             fluxbio_w(i,j,imodn,2)   = ((3.223726628 * lat_t(i,j)*rad2deg-82.3887814313) &
                                          /14/365/24/3600)

!! ALEX: DOC deposition from map in Kanakidou 2012
!         if(lon_t(i,j)*rad2deg<=22.or.(lon_t(i,j)*rad2deg>22.and. &
!            lat_t(i,j)*rad2deg<=37)) then
!             fluxbio_w(i,j,imodc,2)  = ((0.0773694391 * lat_t(i,j)*rad2deg-2.1873307544) &
!                                          *1000/12/365/24/3600) ! Atomic mass C=12
!         else
!             fluxbio_w(i,j,imodc,2)  = ((0.1833 * lat_t(i,j)*rad2deg-6.11667) &
!                                          *1000/12/365/24/3600) ! Atomic mass C=12
!         endif


         AtmDepDOP(i,j)=fluxbio_w(i,j,imodp,2)
         AtmDepDON(i,j)=fluxbio_w(i,j,imodn,2) 
!         AtmDepDOC(i,j)=fluxbio_w(i,j,imodc,2) ! Pas de DOC 
         AtmDepAmmo(i,j)=fluxbio_w(i,j,iammonium,2)
         AtmDepNit(i,j)=fluxbio_w(i,j,initrate,2)

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
