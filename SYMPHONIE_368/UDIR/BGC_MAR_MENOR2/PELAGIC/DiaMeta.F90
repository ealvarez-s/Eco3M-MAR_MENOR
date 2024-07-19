      SUBROUTINE DiaMeta

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
! Meta-model made by Kert Buis : parametrisations with regression     *
! model against full diagenetic model                                 *
! Parametrisations of denitrification (organic matter mineralisation  *
! performed by consumation of nitrate ), nitrification (ammonium      *
! oxidised into nitrate)  and oxic mineralisation                     * 
! In pratice the diagenetic model has been run multiple times unsing  *
! random combinaisons of determining factors and the results were     *
! statically analysed. The rangeof each factor used in the MonteCarlo *
! simulations were determined from the benthic database constituted   *
! in the frame of the MetroMed project                                *  
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
         NO3Flux,NH3Flux,PFlux,SiFlux,ODURed,ODUsolid

      DOUBLE PRECISION :: CMinFast,CMinSlow,CMin,NMin,PMin,SiMin,pAmmon  &
         ,smallb,DepthB,AmmoB,NitrateB,Nitcons

! Conversion /d into /s
      SEC2DAY = 86400.

      smallb = 1.D-8 

!---------------------------------------------------------------------*

!---------------------------------------------------------------------*


      DO J=NBIO1,NBIO2 
      DO I=MBIO1,MBIO2          !09/09/08
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

       IF (BTfunc.EQ.1) THEN
       kFast= cFDet20SED * exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
                         *LOG(BQ10)/10.)
       kSlow= cSDet20SED * exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
                         *LOG(BQ10)/10.)
!      kSiDiss=rSiMin    * exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
!                        *LOG(BQ10)/10.)
       ENDIF       

! Mineralisation rate of fast and slow detritus in mmol/m2/d
       CMinFast = kFast * CBFDET(I,J)
       CMinSlow = kSlow * CBSDET(I,J)
       CMin     = ( CMinFast + CMinSlow )
       NMin     = (CMinFast*NCrFD + CMinSlow*NCrSD)


       ELSEIF(NumBDet.EQ.1) THEN

       IF (BTfunc.EQ.1)                                         &
       DecayRate=cDet20SED*exp((tem_t(i,j,kmin_w(i,j),1)-20.)  &
                          *LOG(BQ10)/10.)

       
 
        CMin    = DecayRate *  CBDet(I,J)
        NMin    = DecayRate *  NBDet(I,J)
       SiMin    = DecayRate * SiBDet(I,J)

       ENDIF

       if (CMin.LE.0.OR.NMin.LE.0.) then 
       NO3Flux = 0.
       NH3Flux = 0.
       GOTO 100
       endif  




! Nitrogen fluxes

! computation of pNit and pdeNit : Kert Buis, MetroMed project 

      DepthB = - depth_w(i,j,1)
      AmmoB  = BIO_t(I,J,kmin_w(I,J),iAmmonium)
      NitrateB = BIO_t(I,J,kmin_w(I,J),iNitrate)


!      if(AmmoB.LT.0)       AmmoB = 0.
!      if(NitrateB.LT.0) NitrateB = 0.
      if(AmmoB.LT.0.001)    AmmoB = 0.001
      if(AmmoB.GT.5.)       AmmoB = 5.
      if(NitrateB.LT.0.012) NitrateB = 0.012
      if(NitrateB.GT.40.) NitrateB = 40.


!      if(CMin.EQ.0.) THEN
!      print*,'CMin=0',I,J,CMin,DepthB,AmmoB,NitrateB
!      print*,'CMin=0 2',CMinFast, CMinSlow, CBFDET(I,J),CBSDET(I,J)
!      endif



! denitrification part (range = 0-1)
      pdeNit =                                               &
            EXP(-1.74053                                     &
                -1.15259*log(CMin+smallb)                    &
                 + 0.06542*log(CMin+smallb)*log(CMin+smallb) &
                 - 0.08106*log(DepthB+smallb)                &
                         *log(DepthB+smallb)                 &
                 + 0.36404*log(CMin+smallb)                  &
                         *log(DepthB+smallb)                 &
                 - 0.03926*log(NitrateB)                     &
                         *log(DecayRate+smallb)              &
                 - 0.08480*log(DecayRate+smallb)             &
                         *log(DepthB+smallb))                &
                       /CMin
! 
      IF( pdeNit.GT.1) pdeNit = 1.
     

! oxic mineralisation part : range 0-1
      pAmmon=                                                &
             EXP(3.17784                                     &
                 + .38381*log(CMin+smallb)*log(O2bw+smallb)  &
                 + .33529*log(CMin+smallb)*log(DecayRate+smallb) &
                 + .09857*log(DepthB+smallb)                 &
                         *log(DepthB+smallb)                 &
                 + .04058*log(DepthB+smallb)                 &
                         *log(DecayRate+smallb)              &
                 - .12045*log(CMin+smallb)*log(CMin+smallb)  &
                 -1.02150*log(DepthB+smallb))                &
                       /CMin



!         IF(I.EQ.142.AND.J.EQ.209) THEN
!         print*,'inside',PROFWK_Z(I,J,1),smallb
!         STOP'dans DiaMeta'
!         ENDIF 


      IF (pAmmon.gt.1) pAmmon=1

      IF (pAmmon+pdeNit.gt.1) then
        pAmmon=1-pdeNit
      ENDIF

! ammonium nitrified, pNit*NMin in mmolN/m2/d

      pNit =                                                 &
             EXP(-1.57401                                    &
                 + .13422*log(CMin+smallb)*log(O2bw+smallb)  &
                 - .03160*log(DepthB+smallb)                 &
                         *log(AmmoB)                         &
                 - .00756*log(DepthB+smallb)                 &
                         *log(DepthB+smallb)                 &
                 - .14603*log(CMin+smallb)                   &
                         *log(AmmoB)                         &
                 + .03514*log(AmmoB)                         &
                         *log(AmmoB)                         &
                 + .49373*log(AmmoB))                        &
                 /NMin

       if(pNit.GT.1) pNit = 1.

!       IF(I.EQ.142.AND.J.EQ.209) THEN
!         print*,'inside',NitrificationB(I,J),NitrificationB(I,J)/NMin   &
!         ,pAmmon,pdeNit,CMin                                            &
!         ,CBDet(I,J),AmmoB,NitrateB,DepthB                              &
!         ,PROFWK_Z(I,J,1), smallb, O2bw
!!        STOP'dans DiaMeta'
!         ENDIF


      IF(NBDet(I,J).GT.0) THEN

! Nitrogen fluxes in mmol/m2/d

! Nitrate consumed by denitrification mM/m2/d
      Nitcons = pdeNit * CMin * 0.8
      DenitrificationB(I,J) = pdeNit * CMin
      NitrificationB(I,J) = pNit * NMin

! Nitrate and ammonium flux into the water column mmol/m2/d
      NO3Flux =  NitrificationB(I,J) - Nitcons 
      NH3Flux =  NMin - NitrificationB(I,J) !+ AmmAdjust   

! oxygen consumed by mineralisation mM/m2/d
      O2Min(I,J) =  pAmmon * CMin
      ODURed = 1.0 - DenitrificationB(I,J) - pAmmon

! oxygen consumed by reduced compounds mM/m2/d
      ODUsolid = 0  !a verifier...
      O2ODU(I,J) = ODURed * CMin * (1-ODUsolid)

      ELSE

      NO3Flux = 0.
      NH3Flux = 0.

      ENDIF


  100 CONTINUE

! Phosphorus and silicium fluxes in mmol/m2/d

       PFlux =  PBDet(I,J) * DecayRate *  pPMin
      SiFlux = SiBDet(I,J) * DecayRate * pSiMin


! Update the benthic compartment in mmol/m2/s
 
! Rate of change of the benthic pool
      IF(NumBDet.EQ.1)                               &
       dCBDet(I,J) =  CDepo(I,J) -  CMin / SEC2DAY

      IF(NumBDet.EQ.2) THEN
       dCBFDet(I,J) =  FDetFlux -  CMinFast / SEC2DAY
       dCBSDet(I,J) =  SDetFlux -  CMinSlow / SEC2DAY
      ENDIF

      dNBDet(I,J) =  NDepo(I,J) -  NMin / SEC2DAY
      dPBDet(I,J) =  PDepo(I,J) -  PMin / SEC2DAY
      dSiBDet(I,J)= SiDepo(I,J) - SiMin / SEC2DAY

! Benthic pool in mmol/m2
      IF(NumBDet.EQ.1)                                   &
       CBDet(I,J) =  CBDet(I,J) +  dCBDet(I,J) * DTI_FW

      IF(NumBDet.EQ.2) THEN
      CBFDet(I,J) =  CBFDet(I,J) +  dCBFDet(I,J) * DTI_FW
      CBSDet(I,J) =  CBSDet(I,J) +  dCBSDet(I,J) * DTI_FW
      ENDIF

      NBDet(I,J) =  NBDet(I,J) +  dNBDet(I,J) * DTI_FW
      PBDet(I,J) =  PBDet(I,J) +  dPBDet(I,J) * DTI_FW
      SiBDet(I,J)= SiBDet(I,J) + dSiBDet(I,J) * DTI_FW

! Nutrient fluxes in mmol/m2/s
! Conversion /d into /s
!     SEC2DAY = 86400.

      fluxbio_w(I,J,iNitrate  ,1) = NO3Flux / SEC2DAY
      fluxbio_w(I,J,iAmmonium ,1) = NH3Flux / SEC2DAY
      fluxbio_w(I,J,iPhosphate,1) =   PFlux / SEC2DAY
      fluxbio_w(I,J,iSilice   ,1) =  SiFlux / SEC2DAY

!       IF(I.EQ.142.AND.J.EQ.209) THEN
!         print*,'inside2'                    &
!         ,fluxbio_w(I,J,iNitrate  ,1)        &
!         ,fluxbio_w(I,J,iAmmonium  ,1)       &
!         ,fluxbio_w(I,J,iPhosphate  ,1)      &
!         ,fluxbio_w(I,J,iSilice  ,1)     
!!        STOP'dans DiaMeta'
!         ENDIF

      CMin_out(I,J) = CMin
      NMin_out(I,J) = NMin


      ENDIF
      ENDDO
      ENDDO 

      RETURN
      END 
      
