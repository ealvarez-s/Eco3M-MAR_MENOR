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

! flux de CO2
     call AirSeaCO2Flux

!flux atmospheriques
      if(idust==1.or.idepnit==1.or.idepammo==1) call AtmDeposition

      fluxNwmed  = 10042.e9/365.0/24.0/3600.0/sumareawmed ! RIBERA 2003
      fluxPwmed  = 357.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
      fluxSiwmed = 709.e9/365.0/24.0/3600.0/sumareawmed !  RIBERA 2003
      fluxNemed  = 6064.e9/365.0/24.0/3600.0/sumareaemed ! RIBERA 2003
      fluxPemed  = 379.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
      fluxSiemed = 700.e9/365.0/24.0/3600.0/sumareaemed !  RIBERA 2003
      fluxAMwmed  = 1300.e9/365./24./3600./sumareawmed  ! MOYENNE MOOSE

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
!      call budget_export_bottom

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
