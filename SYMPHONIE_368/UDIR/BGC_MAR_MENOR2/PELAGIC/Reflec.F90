      SUBROUTINE Reflec

!_____________________________________________________________________*
! 3d ecosystem model                                                  *
!                                                                     *
! LAST REVISION: 6 NOVEMBER 2009                                      *
!                                                                     *
! Implementation: Caroline Ulses                                      *
!                 LA                                                  *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Computes the benthic to water column fluxes                         *
! Here we assume all nitrogen, phosphorus and silicium is returned    *
! to the water column. There is no denitrification                    * 
!_____________________________________________________________________*

!_____________________________________________________________________*
!                                                                     *
! Modifications:                                                      *
!                                                                     *
!_____________________________________________________________________*



!---------------------------------------------------------------------*
! Declarations:


! Global variables

      USE ModuleComBlock
      USE UserDeclaration
      USE ModuleDeclaration
      USE MODULE_PRINCIPAL
      IMPLICIT NONE

! Local variables 
      DOUBLE PRECISION :: SEC2DAY
      REAL pNOs

! Conversion /d into /s
      SEC2DAY = 86400.

!---------------------------------------------------------------------*

       tps_benth_2d=tps_benth_2d+1


       DO J=1,jmax
       DO I=1,imax

! Nutrient fluxes is function of depostion.
! All deposition is returned to the water column.

! part of nitrogen returned as nitrate in the water column.
       pNOs=0.70       

       fluxbio_w(I,J,iNitrate  ,1) =    pNOs  * NDepo(I,J)
       fluxbio_w(I,J,iAmmonium ,1) = (1-pNOs) * NDepo(I,J)
       fluxbio_w(I,J,iPhosphate,1) =            PDepo(I,J)
       fluxbio_w(I,J,iSilice   ,1) =           SiDepo(I,J)

! 2024/12/16 EVA Paste from SimpleDiaMeta 23/03/2017 Ajout variables 2D depots benthiques Alex
! 2024/12/16 EVA Paste from SimpleDiaMeta 05/04/2017 Passage en mmol/m2/d => SEC2DAY Alex

       NO3efflux2d(i,j)=NO3efflux2d(i,j)+(fluxbio_w(i,j,iNitrate,1))*SEC2DAY
       NH4efflux2d(i,j)=NH4efflux2d(i,j)+(fluxbio_w(i,j,iAmmonium,1))*SEC2DAY
       Pefflux2d(i,j)=Pefflux2d(i,j)+fluxbio_w(i,j,iPhosphate,1)*SEC2DAY
       Siefflux2d(i,j)=Siefflux2d(i,j)+fluxbio_w(i,j,iSilice,1)*SEC2DAY
       O2influx2d(i,j)=O2influx2d(i,j)+fluxbio_w(i,j,iOxygen,1)*SEC2DAY
       DICefflux2d(i,j)=DICefflux2d(i,j)+fluxbio_w(i,j,iDIC,1)*SEC2DAY

       ENDDO
       ENDDO


       RETURN
       END
