      SUBROUTINE SimpleDiaMeta

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
! For nitrogen :                                                      *
!                                                                     *
! Following Soetaert et al 2000                                       *
! FluxNit = Nmin *  pNit - Cmin *pdeNit * 0.8                         *
! FluxNit = Nmin *( pNit - C:N  *pdeNit * 0.8)                        *
! FluxAmmo = Nmin * (1- pNit)                                         *
!                                                                     *
! - pNit : following Soetaert et al 2000 0.4 < pNit < 1               *
!   ==> deduced after calibration 
! - Nmin= SedON * Kn
! - Cmin=Nmin * C:N
! - pDeNit : according to 1D modelling (Pastor et al, Denis) :
!   pdeNit = 0.05-0.07  A verifier
! - Kn : according to Pastor et al (Rhone prodelta) 2 classes ,
!  classe of fast detritus : 36.5/an ou 0.09 /j
!  according to Denis (SOFI) : 0.08 /j
!                                                                     *
! For phosphorus and silicium, we return a certain percent of         *
! the mineralisation                                                  *
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
      DOUBLE PRECISION :: SEC2DAY,FDetFlux,SDetFlux,NCDepo,AmmAdjust,  &
         NO3Flux,NH3Flux,PFlux,SiFlux,O2Flux  

      DOUBLE PRECISION :: CMinFast,CMinSlow,CMin,NMin,PMin,SiMin       &
                         ,Nitcons,O2cons

! Conversion /d into /s
      SEC2DAY = 86400.

!---------------------------------------------------------------------*



!---------------------------------------------------------------------*

       tps_benth_2d=tps_benth_2d+1

!       print*,'passe ici Simple'  


       DO J=NBIO1,NBIO2
       DO I=MBIO1,MBIO2
       IF(mask_t(i,j,kmax+1)==1) THEN !>>>>>>>>>>>>>>

      IF(NumBDet.EQ.2) THEN

!    Convert CDepo and NDepo into
! a deposition of fast and slow detritus: SdetFlux, FdetFlux
! Based on the NC ratio of Fast- and Slow detritus
! Units in mol C/m2/s
      FDetFlux = 0.
      SDetFlux = 0.
      AmmAdjust= 0.
      NCDepo   = 0.

      IF (CDepo(I,J) .GT. 0.) THEN                
        NCDepo = NDepo(I,J) / CDepo(I,J)  ! N/C ratio of deposited detritus
      ENDIF

      IF (NCDepo .GT. NCrFD) THEN            ! All deposition as fast detritus
        FDetFlux = CDepo(I,J)
        AmmAdjust = (NCDepo-NCrFD)*CDepo(I,J)    ! 'instantaneous' mineralisation
      ELSEIF (NCDepo .LE. NCrFD .AND. NCDepo .GT. NCrSD) THEN    ! Both SDET and FDET are deposited
        FDetFlux = (NCDepo - NCrSD)/(NCrFD-NCrSD)*CDepo(I,J)
        SDetFlux =  CDepo(I,J)  - FDetFlux
      ELSE       ! Only SDET
        SDetFlux = CDepo(I,J)
      ENDIF

      IF (FDetFlux .LT.0) FDetFlux = 0.
      IF (SDetFlux .LT.0) SDetFlux = 0.


! Degradation rate at the bottom water temperature
       IF(BTfunc.EQ.1) THEN
       kFast= cFdet20SED * exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
                         * LOG(BQ10)/10.)
       kSlow= cSdet20SED * exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
                         *LOG(BQ10)/10.)
!      kSiDiss=rSiMin    * exp((tem_t(i,j,kmin_w(i,j),1))-20.) &
!                        *LOG(BQ10)/10.)
       ENDIF

! Mineralisation rate of fast and slow detritus in mmol/m2/d
      CMinFast = kFast * CBFDet(I,J)
      CMinSlow = kSlow * CBSDet(I,J)
      CMin     = ( CMinFast + CMinSlow )
      NMin     = (CMinFast*NCrFD + CMinSlow*NCrSD)
      
        PMin    = DecayRate *  PBDet(I,J)
       SiMin    = DecayRate * SiBDet(I,J)

       ELSEIF(NumBDet.EQ.1) THEN


       IF(BTfunc.EQ.1)                                             &
       DecayRate = cDet20SED *exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
                             *LOG(BQ10)/10.)


        CMin    = DecayRate *  CBDet(I,J)
        NMin    = DecayRate *  NBDet(I,J)
        PMin    = DecayRate *  PBDet(I,J)
       SiMin    = DecayRate * SiBDet(I,J)

       ENDIF

       if (CMin.LE.0.OR.NMin.LE.0.) then
       NO3Flux = 0.
       NH3Flux = 0.
       GOTO 100
       endif

! Nitrogen fluxes in mmol/m2/d
      IF(NBDet(I,J).GT.0) THEN
      Nitcons = pdeNit * CMin * 0.8
      DenitrificationB(I,J) = CMin * pdeNit 
      NitrificationB(I,J) = NMin * pNit

! Nitrate and ammonium flux into the water column mmol/m2/d 

      NO3Flux = NitrificationB(I,J) - Nitcons

      NH3Flux = NMin - NitrificationB(I,J)

! Oxygen consumed in sediment
! assume that there is no anoxic remineralisation (no realistic on the Rhone prodelta for instance)

      O2cons =   CMin * ( 1-pdeNit ) *1. &  !(OCratio=1)       !12/06/13
               + NitrificationB(I,J) * 2.  !(O2Nratio=2)      !12/06/13

      O2Flux = - O2cons                                        !12/06/13
 
      ELSE

      NO3Flux = 0.
      NH3Flux = 0.

      ENDIF

  100 CONTINUE

! Phosphorus and silicium fluxes in mmol/m2/d

! Test modif Alex 22/03/2017 forcing P sed => water col according to N:P ratio
! EVA commented this part, PFlux gets pPMin from notebook_benthic
!! Bassin W 
!      IF(mask_t(i,j,kmax+1).EQ.1.and.lon_t(i,j)*rad2deg>-5.6)then !debut test
!
!      if ((lon_t(i,j)*rad2deg<10).or.    &
!        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42).or. &
!        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25).or. &
!        (lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25))then
!           PFlux =  (NO3Flux + NH3Flux)/22 
!! Bassin E
!       else 
!           PFlux =  (NO3Flux + NH3Flux)/27
!       endif
!
!       ENDIF 

!! EVA: PFlux gets pPMin from notebook_benthic 
      PFlux =  PBDet(I,J) * DecayRate *  pPMin
      SiFlux = SiBDet(I,J) * DecayRate * pSiMin


!      IF(I.EQ.142.AND.J.EQ.209) THEN
!         print*,'inside1',NitrificationB(I,J),pNit,pdeNit,CMin
!     &   ,CBDet(I,J)
!!        STOP'dans SimpleDiaMeta'
!         ENDIF


! Update the benthic compartment in mmol/m2/s

! Rate of change of the benthic pool
      IF(NumBDet.EQ.1)  &
       dCBDet(I,J) =  CDepo(I,J) -  CMin / SEC2DAY

      IF(NumBDet.EQ.2) THEN
       dCBFDet(I,J) =  FDetFlux -  CMinFast / SEC2DAY
       dCBSDet(I,J) =  SDetFlux -  CMinSlow / SEC2DAY
      ENDIF 

       dNBDet(I,J) =  NDepo(I,J) -  NMin / SEC2DAY
       dPBDet(I,J) =  PDepo(I,J) -  PMin / SEC2DAY
       dSiBDet(I,J)= SiDepo(I,J) - SiMin / SEC2DAY

! Benthic pool in mmol/m2
      IF(NumBDet.EQ.1)  &
       CBDet(I,J) =  CBDet(I,J) +  dCBDet(I,J) * dti_fw

      IF(NumBDet.EQ.2) THEN
      CBFDet(I,J) =  CBFDet(I,J) +  dCBFDet(I,J) * dti_fw
      CBSDet(I,J) =  CBSDet(I,J) +  dCBSDet(I,J) * dti_fw
      ENDIF

!       IF(I.EQ.142.AND.J.EQ.209) THEN
!         print*,'inside1',NitrificationB(I,J),pNit,pdeNit,CMin  &
!         ,CBDet(I,J)
!        STOP'dans SimpleDiaMeta'
!         ENDIF

      NBDet(I,J) =  NBDet(I,J) +  dNBDet(I,J) * dti_fw
      PBDet(I,J) =  PBDet(I,J) +  dPBDet(I,J) * dti_fw
      SiBDet(I,J)= SiBDet(I,J) + dSiBDet(I,J) * dti_fw

! Nutrient fluxes in mmol/m2/s

      fluxbio_w(i,j,iNitrate  ,1) = NO3Flux / SEC2DAY
      fluxbio_w(i,j,iAmmonium ,1) = NH3Flux / SEC2DAY
      fluxbio_w(i,j,iPhosphate,1) =   PFlux / SEC2DAY
      fluxbio_w(i,j,iSilice   ,1) =  SiFlux / SEC2DAY
      fluxbio_w(i,j,iOxygen   ,1) =  O2Flux / SEC2DAY             !12/06/13
      fluxbio_w(i,j,iDIC      ,1) =  CMin / SEC2DAY 


! 23/03/2017 Ajout variables 2D depots benthiques Alex 
! 05/04/2017 Passage en mmol/m2/d => SEC2DAY Alex  

      NO3efflux2d(i,j)=NO3efflux2d(i,j)+(fluxbio_w(i,j,iNitrate,1)) * SEC2DAY 
      NH4efflux2d(i,j)=NH4efflux2d(i,j)+(fluxbio_w(i,j,iAmmonium,1)) * SEC2DAY 
      Pefflux2d(i,j)=Pefflux2d(i,j)+fluxbio_w(i,j,iPhosphate,1)*SEC2DAY
      Siefflux2d(i,j)=Siefflux2d(i,j)+fluxbio_w(i,j,iSilice,1)*SEC2DAY
      O2influx2d(i,j)=O2influx2d(i,j)+fluxbio_w(i,j,iSilice,1)*SEC2DAY
      DICefflux2d(i,j)=DICefflux2d(i,j)+fluxbio_w(i,j,iDIC,1)*SEC2DAY

! Fin ajout variables 2D

!       IF(I.EQ.142.AND.J.EQ.209) THEN
!         print*,'inside2',NitrificationB(I,J),pNit,pdeNit,CMin   &
!         ,CBDet(I,J),CDepo(I,J),CMin / SEC2DAY
!        STOP'dans SimpleDiaMeta'
!         ENDIF


!       IF(I.EQ.142.AND.J.EQ.209) THEN
!       print*,'inside',PBDet(I,J),PMin,pPMin,PFlux
!       STOP'dans SimpleDiaMeta'
!       ENDIF

      CMin_out(I,J) = CMin
      NMin_out(I,J) = NMin

      ENDIF
      ENDDO
      ENDDO 

      RETURN
      END 
      
