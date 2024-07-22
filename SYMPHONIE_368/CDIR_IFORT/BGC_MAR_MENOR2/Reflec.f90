










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
! Computes the 1 to water column fluxes                         *
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

      REAL pNOs

!---------------------------------------------------------------------*



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

       ENDDO
       ENDDO


       RETURN
       END
