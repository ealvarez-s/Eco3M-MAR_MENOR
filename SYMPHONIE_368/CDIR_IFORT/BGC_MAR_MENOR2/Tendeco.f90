










    SUBROUTINE Tendeco(fBIO,fTemperature,fSalinity,fDensity,fSolarradiation,fAlbedoBio,    &
                       fDepthInt,fDepth                                                    & !7_8
                      ,fPH,fmask,ipoint,jpoint,oTENDANCEBIO                                &  !7
                      ,oPrimProdTot,oPrimProdDia,oPrimProdNano,oPrimProdSyne               & !,fmask !11
                      ,oNewProdTot,oNewProdDia,oNewProdNano,oNewProdSyne                   & ! 15
                      ,oRecProdTot,oRecProdDia,oRecProdNano,oRecProdSyne,oRecProdBact      & ! 20
                      ,oRespTot,oExuCTot,oRespPhyto,oRespDia,oRespNano,oRespSyne           & !,oNitrif !26
                      ,oNetPrimProdTot,oNetPrimProdDia,oNetPrimProdNano,oNetPrimProdSyne   & !30
                      ,oTotalGrazC,oGbact,oPPBi,oRespi                                     & !33
                      ,oNitSurf,oRespZoo,oRemPOC,oRemPON,oLostPOC,oLostPON                 & !39
!                      ,oUptBactDON,oExcHeteroN,oExcHeteroNtop,oExcbactNH4top               &
!                      ,oExcBactNH4N,oExcZooAmmoN,oExczooNH4top,oNitrifCOLUMN               & 
                      ,oUptBactDON                                                         & !40
                      ,oTotalGrazCPur,oTotalGrazCHerbi,oNitrifi                            & !oNitrifi,oUptNiti !42
                      ,oExcBactNH4NTOPLAYER,oExcBactNH4NINT,oExcBactNH4NDEEP               & !45
                      ,oExcZooAmmoNTOPLAYER,oExcZooAmmoNINT,oExcZooAmmoNDEEP               & !48
                      ,oExcZooPO4PTOPLAYER,oExcZooPO4PINT,oExcZooPO4PDEEP                  & !51
                      ,oExcBactPO4PTOPLAYER,oExcBactPO4PINT,oExcBactPO4PDEEP               & !54
                      ,oExuSiTOPLAYER,oExuSiINT,oExuSiDEEP                                 &
                      ,oUptNitTOPLAYER,oUptNitINT,oUptNitDEEP                              &
                      ,oUptAmmoTOPLAYER,oUptAmmoINT,oUptAmmoDEEP                           &
                      ,oUptPTOPLAYER,oUptPINT,oUptPDEEP                                    &
                      ,oUptSiTOPLAYER,oUptSiINT,oUptSiDEEP                                 &
                      ,oNitrifTOPLAYER,oNitrifINT,oNitrifDEEP                              &
                      ,oRemSMOPSiTOPLAYER,oRemSMOPSiINT,oRemSMOPSiDEEP                     &
                      ,oRemLMOPSiTOPLAYER,oRemLMOPSiINT,oRemLMOPSiDEEP                     &
                      ,oUptBactPTOPLAYER,oUptBactPINT,oUptBactPDEEP                        &
                      ,oUptBactNTOPLAYER,oUptBactNINT,oUptBactNDEEP,oPH,oPCO2              &
                      ,oRespPi,oRespZi,oRespBi,oGbaci,oGrazHerbii,oGrazToti                &
                      ,oPPBpi,oPPBni,oPPBmi,oRespPpi,oRespPni,oRespPmi,oUptNiti,oUptAmmoi  &
                      ,oUptPi,oUptPBi)                     


!_____________________________________________________________________
! 3d ecosystem model                                                  *
!                                                                     *
! LAST REVISION: 24 MARCH 2008                                        *
!                                                                     *
! Implementation: Caroline Ulses                                      *
!                 NIOO-CEME                                           *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! 1. Transfers symphonie variables in eco3m variables                 *
! 2. Calculates PAR                                                   *
! 3. Calculates the rates of change                                   *
! 4. Transfers eco3m variables in symphonie variables                 *
!                                                                     *
!_____________________________________________________________________*

!_____________________________________________________________________*
!                                                                     *
! Modifications:                                                      *
!                                                                     *
! 21/09/07: Addition of outputs variables                             *
! 24/03/08: Replace 1.e-8 by 1.d-8
!                                                                     *
!_____________________________________________________________________*



!---------------------------------------------------------------------*
! Declarations:                                                       
                                                                     
                                                                     
! Global variables                                                    
                                                                      
      USE ModuleComBlock                                                   
      USE UserDeclaration                                                 
      USE ModuleDeclaration  
      
      use module_principal, only : ssr_w
     
      IMPLICIT NONE                                                        
                                                                      
! Local variables                                                     
                                                                      
      INTEGER :: I,KB,KOUNTBIO,fkount,oIEupho,ipoint,jpoint                                                  
      DOUBLE PRECISION :: fBIO(NumpelagicBoxes,1:xsvardeclare)             
      DOUBLE PRECISION :: fTemperature(NumpelagicBoxes)                   
      DOUBLE PRECISION, INTENT(IN) :: fSolarRadiation
      DOUBLE PRECISION :: fAlbedoBio                                 
      DOUBLE PRECISION :: fDepthInt(NumpelagicBoxesP1)                          
      DOUBLE PRECISION :: fmask(NumpelagicBoxes)
      DOUBLE PRECISION :: fSalinity(NumpelagicBoxes)
      DOUBLE PRECISION :: fDensity(NumpelagicBoxes)
      DOUBLE PRECISION :: fDepth(NumpelagicBoxes)
      DOUBLE PRECISION :: fPH(NumpelagicBoxes)
      DOUBLE PRECISION :: oTENDANCEBIO(NumpelagicBoxes,1:xsvardeclare)
      DOUBLE PRECISION :: oPrimProdTot
      DOUBLE PRECISION :: oPrimProdDia
      DOUBLE PRECISION :: oPrimProdNano
      DOUBLE PRECISION :: oPrimProdSyne
      DOUBLE PRECISION :: oNewProdTot
      DOUBLE PRECISION :: oNewProdDia
      DOUBLE PRECISION :: oNewProdNano
      DOUBLE PRECISION :: oNewProdSyne
      DOUBLE PRECISION :: oRecProdTot
      DOUBLE PRECISION :: oRecProdDia
      DOUBLE PRECISION :: oRecProdSyne
      DOUBLE PRECISION :: oRecProdNano
      DOUBLE PRECISION :: oRecProdBact
     ! DOUBLE PRECISION :: oNitrif
      DOUBLE PRECISION :: oRespTot
      DOUBLE PRECISION :: oExuCTot
      DOUBLE PRECISION :: oRespPhyto
      DOUBLE PRECISION :: oRespSyne
      DOUBLE PRECISION :: oRespNano
      DOUBLE PRECISION :: oRespDia
      DOUBLE PRECISION :: oNetPrimProdTot
      DOUBLE PRECISION :: oNetPrimProdSyne
      DOUBLE PRECISION :: oNetPrimProdNano
      DOUBLE PRECISION :: oNetPrimProdDia
      DOUBLE PRECISION :: oTotalGrazC
      DOUBLE PRECISION :: oGBact
      DOUBLE PRECISION :: oPPBi(NumpelagicBoxes)
      DOUBLE PRECISION :: oRespi(NumpelagicBoxes)
      DOUBLE PRECISION :: oNitSurf
      DOUBLE PRECISION :: oRespZoo
      DOUBLE PRECISION :: oRemPOC
      DOUBLE PRECISION :: oRemPON
      DOUBLE PRECISION :: oLostPOC
      DOUBLE PRECISION :: oLostPON
      DOUBLE PRECISION :: oUptBactDON
     ! DOUBLE PRECISION :: oExcHeteroN
     ! DOUBLE PRECISION :: oExcHeteroNtop
     ! DOUBLE PRECISION :: oExcbactNH4top
     ! DOUBLE PRECISION :: oExcBactNH4N
     ! DOUBLE PRECISION :: oExcZooAmmoN
     ! DOUBLE PRECISION :: oExczooNH4top
     ! DOUBLE PRECISION :: oNitrifCOLUMN
      DOUBLE PRECISION :: oTotalGrazCPur
      DOUBLE PRECISION :: oTotalGrazCHerbi
      DOUBLE PRECISION :: oNitrifi(NumpelagicBoxes)
     ! DOUBLE PRECISION :: oUptNiti(NumpelagicBoxes)
      DOUBLE PRECISION :: oExcBactNH4NTOPLAYER
      DOUBLE PRECISION :: oExcBactNH4NINT
      DOUBLE PRECISION :: oExcBactNH4NDEEP
      DOUBLE PRECISION :: oExcZooAmmoNTOPLAYER
      DOUBLE PRECISION :: oExcZooAmmoNINT
      DOUBLE PRECISION :: oExcZooAmmoNDEEP
      DOUBLE PRECISION :: oExcZooPO4PTOPLAYER
      DOUBLE PRECISION :: oExcZooPO4PINT
      DOUBLE PRECISION :: oExcZooPO4PDEEP
      DOUBLE PRECISION :: oExcBactPO4PTOPLAYER
      DOUBLE PRECISION :: oExcBactPO4PINT
      DOUBLE PRECISION :: oExcBactPO4PDEEP
      DOUBLE PRECISION :: oExuSiTOPLAYER
      DOUBLE PRECISION :: oExuSiINT
      DOUBLE PRECISION :: oExuSiDEEP
      DOUBLE PRECISION :: oUptNitTOPLAYER
      DOUBLE PRECISION :: oUptNitINT
      DOUBLE PRECISION :: oUptNitDEEP
      DOUBLE PRECISION :: oUptAmmoTOPLAYER
      DOUBLE PRECISION :: oUptAmmoINT
      DOUBLE PRECISION :: oUptAmmoDEEP
      DOUBLE PRECISION :: oUptPTOPLAYER
      DOUBLE PRECISION :: oUptPINT
      DOUBLE PRECISION :: oUptPDEEP
      DOUBLE PRECISION :: oUptSiTOPLAYER
      DOUBLE PRECISION :: oUptSiINT
      DOUBLE PRECISION :: oUptSiDEEP
      DOUBLE PRECISION :: oNitrifTOPLAYER
      DOUBLE PRECISION :: oNitrifINT
      DOUBLE PRECISION :: oNitrifDEEP
      DOUBLE PRECISION :: oRemSMOPSiTOPLAYER
      DOUBLE PRECISION :: oRemSMOPSiINT
      DOUBLE PRECISION :: oRemSMOPSiDEEP
      DOUBLE PRECISION :: oRemLMOPSiTOPLAYER
      DOUBLE PRECISION :: oRemLMOPSiINT
      DOUBLE PRECISION :: oRemLMOPSiDEEP
      DOUBLE PRECISION :: oUptBactPTOPLAYER
      DOUBLE PRECISION :: oUptBactPINT
      DOUBLE PRECISION :: oUptBactPDEEP
      DOUBLE PRECISION :: oUptBactNTOPLAYER
      DOUBLE PRECISION :: oUptBactNINT
      DOUBLE PRECISION :: oUptBactNDEEP
      DOUBLE PRECISION :: oPH(NumpelagicBoxes)
      DOUBLE PRECISION :: oPCO2
      DOUBLE PRECISION :: oRespPi(NumpelagicBoxes)
      DOUBLE PRECISION :: oRespZi(NumpelagicBoxes)
      DOUBLE PRECISION :: oRespBi(NumpelagicBoxes)
      DOUBLE PRECISION :: oGbaci(NumpelagicBoxes)
      DOUBLE PRECISION :: oGrazHerbii(NumpelagicBoxes)
      DOUBLE PRECISION :: oGrazToti(NumpelagicBoxes)
      DOUBLE PRECISION :: oRespPpi(NumpelagicBoxes)
      DOUBLE PRECISION :: oRespPni(NumpelagicBoxes)
      DOUBLE PRECISION :: oRespPmi(NumpelagicBoxes)
      DOUBLE PRECISION :: oPPBpi(NumpelagicBoxes)
      DOUBLE PRECISION :: oPPBni(NumpelagicBoxes)
      DOUBLE PRECISION :: oPPBmi(NumpelagicBoxes)
      DOUBLE PRECISION :: oUptNiti(NumpelagicBoxes)
      DOUBLE PRECISION :: oUptAmmoi(NumpelagicBoxes)
      DOUBLE PRECISION :: oUptPi(NumpelagicBoxes)
      DOUBLE PRECISION :: oUptPBi(NumpelagicBoxes)
!---------------------------------------------------------------------*

!---------------------------------------------------------------------*
!                                                                     *
! 1. Symphonie to Eco3m                                               *
!                                                                     *
!---------------------------------------------------------------------*
      

! States variables

      DO  I = 1 , NumpelagicBoxes    

        ZooNanoC(I)  = fBIO(NumpelagicBoxes-I+1,iZooNanoC)
        ZooMicroC(I) = fBIO(NumpelagicBoxes-I+1,iZooMicroC)
        ZooMesoC(I)  = fBIO(NumpelagicBoxes-I+1,iZooMesoC)
        SyneC(I)     = fBIO(NumpelagicBoxes-I+1,iSyneC)
        SyneN(I)     = fBIO(NumpelagicBoxes-I+1,iSyneN)
        SyneP(I)     = fBIO(NumpelagicBoxes-I+1,iSyneP)
        SyneChl(I)   = fBIO(NumpelagicBoxes-I+1,iSyneChl)
        NanoC(I)     = fBIO(NumpelagicBoxes-I+1,iNanoC)
        NanoN(I)     = fBIO(NumpelagicBoxes-I+1,iNanoN)
        NanoP(I)     = fBIO(NumpelagicBoxes-I+1,iNanoP)
        NanoChl(I)   = fBIO(NumpelagicBoxes-I+1,iNanoChl)
        DiaC(I)      = fBIO(NumpelagicBoxes-I+1,iDiaC)
        DiaN(I)      = fBIO(NumpelagicBoxes-I+1,iDiaN)
        DiaP(I)      = fBIO(NumpelagicBoxes-I+1,iDiaP)
        DiaChl(I)    = fBIO(NumpelagicBoxes-I+1,iDiaChl)
        DiaSi(I)     = fBIO(NumpelagicBoxes-I+1,iDiaSi)
        BactC(I)     = fBIO(NumpelagicBoxes-I+1,iBactC)
        S_MOPC(I)    = fBIO(NumpelagicBoxes-I+1,iSMOPC)
        S_MOPN(I)    = fBIO(NumpelagicBoxes-I+1,iSMOPN)
        S_MOPP(I)    = fBIO(NumpelagicBoxes-I+1,iSMOPP)
        S_MOPChl(I)  = fBIO(NumpelagicBoxes-I+1,iSMOPChl)
        S_MOPSi(I)   = fBIO(NumpelagicBoxes-I+1,iSMOPSi)
        L_MOPC(I)    = fBIO(NumpelagicBoxes-I+1,iLMOPC)
        L_MOPN(I)    = fBIO(NumpelagicBoxes-I+1,iLMOPN)
        L_MOPP(I)    = fBIO(NumpelagicBoxes-I+1,iLMOPP)
        L_MOPSi(I)   = fBIO(NumpelagicBoxes-I+1,iLMOPSi)
        MODC(I)      = fBIO(NumpelagicBoxes-I+1,iMODC)
        MODN(I)      = fBIO(NumpelagicBoxes-I+1,iMODN)
        MODP(I)      = fBIO(NumpelagicBoxes-I+1,iMODP)
        Nitrate(I)   = fBIO(NumpelagicBoxes-I+1,iNitrate)
        Ammonium(I)  = fBIO(NumpelagicBoxes-I+1,iAmmonium)
        Phosphate(I) = fBIO(NumpelagicBoxes-I+1,iPhosphate)
        Silice(I)    = fBIO(NumpelagicBoxes-I+1,iSilice)
        Oxygen(I)    = fBIO(NumpelagicBoxes-I+1,ioxygen)
        DIC(I)       = fBIO(NumpelagicBoxes-I+1,iDIC)
        Alkalinity(I)= fBIO(NumpelagicBoxes-I+1,iAlkalinity)

        maskbio(I)   =fmask(NumpelagicBoxes-I+1)

      ENDDO


! Solar radiation & Albedo

      SolarRadiation = fSolarRadiation
      AlbedoBio = fAlbedoBio


! Temperature and DepthInt


!!!! ????Convention pour les prof a verifier
      DO I = 1 , NumPelagicBoxesP1
 
        DepthInt(I)    = - fDepthInt(NumpelagicBoxesP1-I+1)

      ENDDO

      DO I = 1 , NumPelagicBoxes

        Temperature(I) = max(1.d-10,fTemperature(NumpelagicBoxes-I+1))         !24/03/08
        Salinity(I) = max(1.d-10,fSalinity(NumpelagicBoxes-I+1))
        DensEco(I) = max(1.d-10,fDensity(NumpelagicBoxes-I+1))

        Depth(I)= - fDepth(NumpelagicBoxes-I+1) 

        Thickness(I)   =   fDepthInt(NumpelagicBoxesP1-I+1)                  &
                         - fDepthInt(NumpelagicBoxesP1-I)

        PH(I)        = fPH(NumpelagicBoxes-I+1)

      ENDDO


!       WRITE(6,*)'fBIO'
!       WRITE(6,*)(fBIO(NumpelagicBoxes,KB),KB=1,28)
!       WRITE(6,*)'Temperature Tendeco'
!       WRITE(6,*)(Temperature(I),I=1,NumpelagicBoxes)
!       WRITE(6,*)'DepthInt'
!       WRITE(6,*)(DepthInt(I),I=1,NumpelagicBoxesP1)
!       WRITE(6,*)'Thickness'
!       WRITE(6,*)(Thickness(I),I=1,NumpelagicBoxes)
!       WRITE(6,*)'SolarRadiation',SolarRadiation
!       WRITE(6,*)'AlbedoBio',AlbedoBio


       if(ipoint.eq.195.and.jpoint.eq.202) then
!         print*,'Tendeco'
!         print*,(fBIO(I,30),I=1,NumpelagicBoxes)
!         print*,'Temperature Tendeco'
!         print*,(Temperature(I),I=1,NumpelagicBoxes)
!         print*,'Salinity Tendeco'
!         print*,(Salinity(I),I=1,NumpelagicBoxes)
!         print*,'DepthInt'
!         print*,(DepthInt(I),I=1,NumpelagicBoxesP1)
!         print*,'Thickness'
!         print*,(Thickness(I),I=1,NumpelagicBoxes)
         print*,'SolarRadiation',SolarRadiation
         print*,'AlbedoBio',AlbedoBio
!               STOP'dans Tendeco'
        endif 

!---------------------------------------------------------------------*
!                                                                     *
! 1. Symphonie to Eco3m                                               *
! END                                                                 *
!---------------------------------------------------------------------*


!---------------------------------------------------------------------*
!                                                                     *
! 2. Calculation of the Photosynthetically active radiation (PAR)     *
! BEGINNING                                                           *
!---------------------------------------------------------------------*


      CALL DynamicsLight(oIEupho,ipoint,jpoint)

!---------------------------------------------------------------------*
!                                                                     *
! 2. Calculation of the Photosynthetically active radiation (PAR)     *
! END                                                                 *
!---------------------------------------------------------------------*



!---------------------------------------------------------------------*
!                                                                     *
! 3. Calculation of the rates of change                               *
! BEGINNING                                                           *
!---------------------------------------------------------------------*

      CALL DynamicsEco3m

      CALL CO2_Pressure_PH

    if(ipoint.eq.195.and.jpoint.eq.202) then
!    print*,'Tendeco',pCO2(1),pH(1),CO2(1),K1CO2(1),K2CO2(1), &
!           DIC(1)*1000/DensEco(1), &
!           Alkalinity(1)*1000/DensEco(1)
    do I=1,NumpelagicBoxes
    print*,'Tendeco2',dNitrate(I)
    enddo
    !    print*,'Tendeco3',K0_CO2(1),CO2(1)*1000/DensEco(1)/K0_CO2(1),pCO2(1),DensEco(1)
    endif


!---------------------------------------------------------------------*
!                                                                     *
! 3. Calculation of the rates o change                                *
! END                                                                 *
!---------------------------------------------------------------------*



!---------------------------------------------------------------------*
!                                                                     *
! 4. Eco3m to Symphonie                                               *
! BEGINNING                                                           *
!---------------------------------------------------------------------*
  
       DO I = 1 , NumpelagicBoxes

         oTENDANCEBIO(NumpelagicBoxes-I+1,iZOONanoC)  = dZOONanoC(I)     
         oTENDANCEBIO(NumpelagicBoxes-I+1,iZOOMicroC) = dZOOMicroC(I)     
         oTENDANCEBIO(NumpelagicBoxes-I+1,iZOOMesoC)  = dZOOMesoC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSyneC)     = dSyneC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSyneN)     = dSyneN(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSyneP)     = dSyneP(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSyneChl)   = dSyneChl(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iNanoC)     = dNanoC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iNanoN)     = dNanoN(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iNanoP)     = dNanoP(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iNanoChl)   = dNanoChl(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iDIAC)      = dDIAC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iDIAN)      = dDIAN(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iDIAP)      = dDIAP(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iDIAChl)    = dDIAChl(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iDIASi)     = dDIASi(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iBactC)     = dBactC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSMOPC)     = dS_MOPC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSMOPN)     = dS_MOPN(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSMOPP)     = dS_MOPP(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSMOPChl)   = dS_MOPChl(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSMOPSi)    = dS_MOPSi(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iLMOPC)     = dL_MOPC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iLMOPN)     = dL_MOPN(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iLMOPP)     = dL_MOPP(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iLMOPSi)    = dL_MOPSi(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iMODC)      = dMODC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iMODN)      = dMODN(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iMODP)      = dMODP(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iNitrate)   = dNitrate(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iAmmonium)  = dAmmonium(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iPhosphate) = dPhosphate(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iSilice)    = dSilice(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,ioxygen)    = doxygen(I) 
         oTENDANCEBIO(NumpelagicBoxes-I+1,iDIC)       = dDIC(I)
         oTENDANCEBIO(NumpelagicBoxes-I+1,iAlkalinity)= dAlkalinity(I)
         oPH(NumpelagicBoxes-I+1)=PH(I)
!        oPCO2(NumpelagicBoxes-I+1)=pCO2(I)


! Flux from Eco3m to Diagnostics
         oPPBi(NumpelagicBoxes-I+1)       = PPBi(I)
         oRespi(NumpelagicBoxes-I+1)      = RespCi(I)       
         oRespPi(NumpelagicBoxes-I+1)     = RespCPi(I)
         oRespZi(NumpelagicBoxes-I+1)     = RespCZi(I)
         oRespBi(NumpelagicBoxes-I+1)     = RespCBi(I)
         oNitrifi(NumpelagicBoxes-I+1)    = Nitrifi(I)
         oGbaci(NumpelagicBoxes-I+1)      = Gbaci(I)
         oGrazHerbii(NumpelagicBoxes-I+1) = GrazHerbii(I)
         oGrazToti(NumpelagicBoxes-I+1)   = GrazCi(I)
         oUptNiti(NumpelagicBoxes-I+1)    = vUptNit(I)         
         oUptAmmoi(NumpelagicBoxes-I+1)   = vUptAmmo(I)
         oPPBpi(NumpelagicBoxes-I+1)      = PPBpi(I)
         oRespPpi(NumpelagicBoxes-I+1)    = RespCPpi(I)
         oPPBni(NumpelagicBoxes-I+1)      = PPBni(I)
         oRespPni(NumpelagicBoxes-I+1)    = RespCPni(I)
         oPPBmi(NumpelagicBoxes-I+1)      = PPBmi(I)
         oRespPmi(NumpelagicBoxes-I+1)    = RespCPmi(I)
         oUptPi(NumpelagicBoxes-I+1)      = vUptPhytoP(I)
         oUptPBi(NumpelagicBoxes-I+1)     = vUptBactP(I)

       ENDDO

! Conversion in s
       DO KB = 1 , xsvardeclare
       DO I  = 1 , NumpelagicBoxes

         oTENDANCEBIO(I,KB) = oTENDANCEBIO(I,KB) / 3600.

       ENDDO
       ENDDO


              if(ipoint.eq.195.and.jpoint.eq.202) then
                print*,'Tendeco'
                print*,(oTENDANCEBIO(I,30),I=1,NumpelagicBoxes)
       !         print*,'Temperature Tendeco'
       !         print*,(Temperature(I),I=1,NumpelagicBoxes)
       !         print*,'Salinity Tendeco'
       !         print*,(Salinity(I),I=1,NumpelagicBoxes)
       !         print*,'DepthInt'
       !         print*,(DepthInt(I),I=1,NumpelagicBoxesP1)
       !         print*,'Thickness'
       !         print*,(Thickness(I),I=1,NumpelagicBoxes)
       !         print*,'SolarRadiation',SolarRadiation

       !WRITE(6,*)'AlbedoBio',AlbedoBio
       !               STOP'dans Tendeco'
              endif


       oPCO2=pCO2(1)


! Euphotic Layer Depth
!      oEuphoticLayerDepth=EuphoticLayerDepth 

! Brut primary production
       oPrimProdTot  = TotalPPB
       oPrimProdSyne = TotalPPBSyne
       oPrimProdNano = TotalPPBNano
       oPrimProdDia  = TotalPPBDia

! net production
       oNetPrimProdTot  = TotalNetPPB
       oNetPrimProdSyne = TotalNetPPBSyne
       oNetPrimProdNano = TotalNetPPBNano
       oNetPrimProdDia  = TotalNetPPBDia

! new production: uptake of nitrate
       oNewProdTot  = TotalUptNit
       oNewProdSyne = TotalUptSyneNit
       oNewProdNano = TotalUptNanoNit
       oNewProdDia  = TotalUptDiaNit

! recycled production: ammonium uptake
       oRecProdTot  = TotalUptAmmo
       oRecProdSyne = TotalUptSyneAmmo
       oRecProdDia  = TotalUptDiaAmmo
       oRecProdNano = TotalUptNanoAmmo
       oRecProdBact = TotalUptBactAmmo

! Grazing ZooC                            !10/06/10
       oTotalGrazC = TotalGrazC           !10/06/10
! claude
       oTotalGrazCPur = TotalGrazCPur
       oTotalGrazCHerbi = TotalGrazCHerbi  

! nitrification COMMENTE ALEX 11/01/18
!      oNitrif  = TotalNitrif
!      oNitrifCOLUMN = TotalNitrifColumn

! respiration
       oRespTot   = TotalResp
       oRespPhyto = TotalRespPhyto
       oRespSyne  = TotalRespSyne
       oRespNano  = TotalRespNano
       oRespDia   = TotalRespDia
       oRespZoo   = TotalRespZoo

! exudation
       oExuCTot = TotalExuC
 
! croissance bacterienne 
       oGBact  = TotalGrowthBact

! remineralisation de POC et PON
       oRemPOC = TotalRemSMOPC + TotalRemLMOPC
       oRemPON = TotalRemSMOPN + TotalRemLMOPN

! lost of POC in DOC : messyfeed + mortbact
       oLostPOC = TotalMessyfeedC + TotalMortBactC

! lost of PON in DON : messyfeedN + mortbact + exsu phyto
       oLostPON = TotalMessyfeedN + TotalMortBactC * NCBact / 12. * 14. + TotalExuN

! Uptake of DON by bacteria (DON -> PON)
       oUptBactDON = TotalUptBactMODN

! COMMENTE PAR ALEX 11/01/18
!! Excretion of DIN by heterotrophs (PON -> DIN)
!      oExcHeteroN = TotalBactExcNH4 + TotalExcZooAmmo
!      oExcBactNH4N = TotalBactExcNH4
!      oExcZooAmmoN = TotalExcZooAmmo
!      oExcHeteroNtop = TOTALBACTEXCNH4TOPLAYER + TOTALEXCZOOAMMOTOPLAYER 
!      oExcbactNH4top = TOTALBACTEXCNH4TOPLAYER
!      oExczooNH4top = TOTALEXCZOOAMMOTOPLAYER
! FIN COMMENTE PAR ALEX 11/01/18

! BY LAYER (11/01/2018)

       oExcBactNH4NTOPLAYER = TotalBactExcNH4TOPLAYER
       oExcZooAmmoNTOPLAYER = TotalExcZooAmmoTOPLAYER
       oExcZooPO4PTOPLAYER = TotalExcZooPO4TOPLAYER
       oExcBactPO4PTOPLAYER = TotalBactExcPO4TOPLAYER
       oExuSiTOPLAYER = TotalExuSiTOPLAYER
! new production: uptake of nitrate
       oUptNitTOPLAYER  = TotalUptNitTOPLAYER
! recycled production: ammonium uptake
       oUptAmmoTOPLAYER  = TotalUptAmmoTOPLAYER
       oUptBactNTOPLAYER  = TotalBactUptAmmoTOPLAYER
! uptake of phosphate
       oUptPTOPLAYER  = TotalUptPTOPLAYER
       oUptBactPTOPLAYER = TotalBactUptPTOPLAYER
! uptake of silicate
       oUptSiTOPLAYER  = TotalUptSiTOPLAYER
! nitrification
       oNitrifTOPLAYER  = TotalNitrifTOPLAYER
! remineralisation Si
       oRemSMOPSiTOPLAYER = TotalRemSMOPSiTOPLAYER
       oRemLMOPSiTOPLAYER = TotalRemLMOPSiTOPLAYER

       oExcBactNH4NINT = TotalBactExcNH4INT
       oExcZooAmmoNINT = TotalExcZooAmmoINT
       oExcZooPO4PINT = TotalExcZooPO4INT
       oExcBactPO4PINT = TotalBactExcPO4INT
       oExuSiINT = TotalExuSiINT
! new production: uptake of nitrate
       oUptNitINT  = TotalUptNitINT
! recycled production: ammonium uptake
       oUptAmmoINT  = TotalUptAmmoINT
       oUptBactNINT  = TotalBactUptAmmoINT
! uptake of phosphate
       oUptPINT  = TotalUptPINT
       oUptBactPINT  = TotalBactUptPINT
! uptake of silicate
       oUptSiINT  = TotalUptSiINT
! nitrification
       oNitrifINT  = TotalNitrifINT
! remineralisation Si
       oRemSMOPSiINT = TotalRemSMOPSiINT
       oRemLMOPSiINT = TotalRemLMOPSiINT

       oExcBactNH4NDEEP = TotalBactExcNH4DEEP
       oExcZooAmmoNDEEP = TotalExcZooAmmoDEEP
       oExcZooPO4PDEEP = TotalExcZooPO4DEEP
       oExcBactPO4PDEEP = TotalBactExcPO4DEEP
       oExuSiDEEP = TotalExuSiDEEP
! new production: uptake of nitrate
       oUptNitDEEP  = TotalUptNitDEEP
! recycled production: ammonium uptake
       oUptAmmoDEEP  = TotalUptAmmoDEEP
       oUptBactNDEEP  = TotalBactUptAmmoDEEP
! uptake of phosphate
       oUptPDEEP  = TotalUptPDEEP
       oUptBactPDEEP  = TotalBactUptPDEEP
! uptake of silicate
       oUptSiDEEP  = TotalUptSiDEEP
! nitrification
       oNitrifDEEP  = TotalNitrifDEEP
! remineralisation Si
       oRemSMOPSiDEEP = TotalRemSMOPSiDEEP
       oRemLMOPSiDEEP = TotalRemLMOPSiDEEP

!      KOUNTBIO=KOUNTBIO+1
!       write(6,*)'oTendancebio'
!     OPEN(UNIT=3,FILE='Verif_Tendance',ACCESS='APPEND')
!       WRITE(6,*)(oTENDANCEBIO(NumpelagicBoxes,KB),KB=1,28)
!      WRITE(6,*)oEuphoticLayerDepth
!     CLOSE(3)
!      OPEN(UNIT=3,FILE='Verif_Eupho',ACCESS='APPEND')
!      WRITE(3,103)KOUNTBIO,PAR_Z(1),SolarRadiation,oEuphoticLayerDepth  &
!                ,(PAR_Z(I),I=1,NumpelagicBoxes)
!      CLOSE(3)
!      WRITE(6,*)oEuphoticLayerDepth  

 103  FORMAT(I3,43(1X,F8.4))

!      STOP'dans Tendeco'


!---------------------------------------------------------------------*
!                                                                     *
! 4. Eco3m to Symphonie                                               *
! END                                                                 *
!---------------------------------------------------------------------*


    
   END SUBROUTINE Tendeco
