










!######################################################################
!
!
! 1-D ecosystem model
!
!
! REFERENCE: none
!
! Implementation: Karline Soetaert - Caroline ULSES
!                 NIOO-CEME
!
!   modif: Fayçal Kessouri (Oxygen)
!
!######################################################################

!---------------------------------------------------------------------
! Modifications:
!
! 18/07/07: Bug dans la partie Messy Feeding du Microzooplancton
!           corrigé
! 30/07/07: - Fonctions gml de Caperon et Meyer pour
!             Nanophytoplancton et Diatomees
!           - Multiplication par une fonction temperature et PAR dans
!             le calcul de la nitrification
!           - Verification du signe de l'excretion du zooplancton et
!             des bacteries
!           - Test sur le signe de fZooMicro et fZooMeso
!           - Mortalite des phytoplanctons en fonction de la
!             temperature
!           - Ajout de ZoopropMesoMOP dans la mortalite du
!             Mesozooplancton
! 03/08/07: - Bug dans la somme de l'exudation en carbone corrigé
!           - Fontion gml: dans condition AND au lieu de OR
!           - excrétion: on enlève les cas redondants
! 21/08/07: - Toute la chlorophylle broutee part dans l'egestion
!           - Correction du calcul de la respiration basale des Zoo
! 22/09/08: - Je corrige le bug sur le ratio ChlN
!           - Tests sur fonction mortalite du mesozoo
! 27/11/12: - Rajout au programme les fonctions oxygen 
!             et Oxygen Demand Unit comme variables d'etat
!--------------------------------------------------------------------


             !<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>!
             !<<<                                            >>>!
             !<<<        MAIN SUBROUTINES AND MODULES        >>>!
             !<<<                                            >>>!
             !<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>!


!********************************************************************


!********************************************************************
!  INITIALISATION 
!********************************************************************

!   I remove this section for the 3d : the initialisation is in
!   InitPelagic.F

!********************************************************************
!  DYNAMICS
!********************************************************************

   SUBROUTINE DynamicsEco3m

   USE ModuleComBlock
   use module_parallele
   IMPLICIT NONE

   INTEGER :: I

   DOUBLE PRECISION :: fZooMicro,fZooMeso,fZooNano,TempFac, PARZ
   DOUBLE PRECISION :: ZooNanoGrazZooNanoC,ZooNanograzNanoC,ZooNanoGrazSyneC,ZooNanoGrazBactC,ZooNanoGrazC
   DOUBLE PRECISION :: ZooNanoGrazZooNanoN,ZooNanograzNanoN,ZooNanoGrazSyneN,ZooNanoGrazBactN,ZooNanoGrazN
   DOUBLE PRECISION :: ZooNanoGrazZooNanoP,ZooNanograzNanoP,ZooNanoGrazSyneP,ZooNanoGrazBactP,ZooNanoGrazP
   DOUBLE PRECISION :: ZooNanograzDiaChl,ZooNanograzNanoChl,ZooNanoGrazSyneChl,ZooNanoGrazSMOPChl,ZooNanoGrazChl
   DOUBLE PRECISION :: ZooNanograzDiaSi,ZooNanograzSMOPSi,ZooNanoGrazSi
   DOUBLE PRECISION :: ZooMicroGrazZooMicroC,ZooMicrograzDiaC,ZooMicrograzZooNanoC,ZooMicrograzNanoC,ZooMicroGrazSyneC
   DOUBLE PRECISION :: ZooMicroGrazBactC,ZooMicroGrazSMOPC,ZooMicroGrazC
   DOUBLE PRECISION :: ZooMicroGrazZooMicroN,ZooMicrograzDiaN,ZooMicrograzZooNanoN,ZooMicrograzNanoN,ZooMicroGrazSyneN
   DOUBLE PRECISION :: ZooMicroGrazBactN,ZooMicroGrazSMOPN,ZooMicroGrazN
   DOUBLE PRECISION :: ZooMicroGrazZooMicroP,ZooMicrograzDiaP,ZooMicrograzZooNanoP,ZooMicrograzNanoP,ZooMicroGrazSyneP
   DOUBLE PRECISION :: ZooMicroGrazBactP,ZooMicroGrazSMOPP,ZooMicroGrazP
   DOUBLE PRECISION :: ZooMicrograzDiaChl,ZooMicrograzNanoChl,ZooMicroGrazSyneChl,ZooMicroGrazSMOPChl,ZooMicroGrazChl
   DOUBLE PRECISION :: ZooMicrograzDiaSi,ZooMicrograzSMOPSi,ZooMicroGrazSi
   DOUBLE PRECISION :: ZooMesograzZooMicroC,ZooMesoGrazDiaC,ZooMesoGrazNanoC,ZooMesoGrazLMOPC,ZooMesoGrazSMOPC,ZooMesoGrazC
   DOUBLE PRECISION :: ZooMesograzZooMicroN,ZooMesoGrazDiaN,ZooMesoGrazNanoN,ZooMesoGrazLMOPN,ZooMesoGrazSMOPN,ZooMesoGrazN
   DOUBLE PRECISION :: ZooMesograzZooMicroP,ZooMesoGrazDiaP,ZooMesoGrazNanoP,ZooMesoGrazLMOPP,ZooMesoGrazSMOPP,ZooMesoGrazP
   DOUBLE PRECISION :: ZooMesoGrazDiaChl,ZooMesoGrazNanoChl,ZooMesoGrazSMOPChl,ZooMesoGrazChl
   DOUBLE PRECISION :: ZooMesograzDiaSi,ZooMesoGrazSMOPSi,ZooMesoGrazLMOPSi,ZooMesoGrazSi
   DOUBLE PRECISION :: ZoograzNanoC,ZooGrazZooNanoC,ZooGrazSyneC,ZooGrazSMOPC,ZooGrazZooMicroC,ZooGrazDiaC,ZooGrazLMOPC
   DOUBLE PRECISION :: ZooGrazBactC,ZooGrazC
   DOUBLE PRECISION :: ZoograzNanoN,ZooGrazZooNanoN,ZooGrazSyneN,ZooGrazSMOPN,ZooGrazZooMicroN,ZooGrazDiaN,ZooGrazLMOPN
   DOUBLE PRECISION :: ZooGrazBactN,ZooGrazN
   DOUBLE PRECISION :: ZoograzNanoP,ZooGrazZooNanoP,ZooGrazSyneP,ZooGrazSMOPP,ZooGrazZooMicroP,ZooGrazDiaP,ZooGrazLMOPP
   DOUBLE PRECISION :: ZooGrazBactP,ZooGrazP
   DOUBLE PRECISION :: ZoograzNanoChl,ZooGrazSyneChl,ZooGrazDiaChl,ZooGrazSMOPChl,ZooGrazChl
   DOUBLE PRECISION :: ZoograzSMOPSi,ZoograzLMOPSi,ZoograzDiaSi,ZooGrazSi

   DOUBLE PRECISION :: ZooNanoNanoMessyFeedMODC,NanoNanoMessyFeedMODC,SyneNanoMessyFeedMODC,BactNanoMessyFeedMODC
   DOUBLE PRECISION :: NanoMessyFeedCMOD
   DOUBLE PRECISION :: ZooNanoNanoMessyFeedMODN,NanoNanoMessyFeedMODN,SyneNanoMessyFeedMODN,BactNanoMessyFeedMODN
   DOUBLE PRECISION :: NanoMessyFeedNMOD
   DOUBLE PRECISION :: ZooNanoNanoMessyFeedMODP,NanoNanoMessyFeedMODP,SyneNanoMessyFeedMODP,BactNanoMessyFeedMODP
   DOUBLE PRECISION :: NanoMessyFeedPMOD
   DOUBLE PRECISION :: ZooMicroMicroMessyFeedMODC,DiaMicroMessyFeedMODC,ZooNanoMicroMessyFeedMODC,NanoMicroMessyFeedMODC
   DOUBLE PRECISION :: SyneMicroMessyFeedMODC,BactMicroMessyFeedMODC,SMOPMicroMessyFeedMODC
   DOUBLE PRECISION :: MicroMessyFeedCMOD
   DOUBLE PRECISION :: ZooMicroMicroMessyFeedMODN,DiaMicroMessyFeedMODN,ZooNanoMicroMessyFeedMODN,NanoMicroMessyFeedMODN
   DOUBLE PRECISION :: SyneMicroMessyFeedMODN,BactMicroMessyFeedMODN,SMOPMicroMessyFeedMODN
   DOUBLE PRECISION :: MicroMessyFeedNMOD
   DOUBLE PRECISION :: ZooMicroMicroMessyFeedMODP,DiaMicroMessyFeedMODP,ZooNanoMicroMessyFeedMODP,NanoMicroMessyFeedMODP
   DOUBLE PRECISION :: SyneMicroMessyFeedMODP,BactMicroMessyFeedMODP,SMOPMicroMessyFeedMODP
   DOUBLE PRECISION :: MicroMessyFeedPMOD
   DOUBLE PRECISION :: ZooMicroMesoMessyFeedMODC,DiaMesoMessyFeedMODC,NanoMesoMessyFeedMODC,LMOPMesoMessyFeedMODC &
                       ,SMOPMesoMessyFeedMODC
   DOUBLE PRECISION :: MesoMessyFeedCMOD
   DOUBLE PRECISION :: ZooMicroMesoMessyFeedMODN,DiaMesoMessyFeedMODN,NanoMesoMessyFeedMODN,LMOPMesoMessyFeedMODN  &
                      ,SMOPMesoMessyFeedMODN
   DOUBLE PRECISION :: MesoMessyFeedNMOD
   DOUBLE PRECISION :: ZooMicroMesoMessyFeedMODP,DiaMesoMessyFeedMODP,NanoMesoMessyFeedMODP,LMOPMesoMessyFeedMODP  &
                        ,SMOPMesoMessyFeedMODP
   DOUBLE PRECISION :: MesoMessyFeedPMOD
   DOUBLE PRECISION :: MessyFeedCMOD,MessyFeedNMOD,MessyFeedPMOD

   DOUBLE PRECISION :: ZooNanoNanoEgeSMOPC,NanoNanoEgeSMOPC,SyneNanoEgeSMOPC,BactNanoEgeSMOPC
   DOUBLE PRECISION :: NanoEgeCSMOP
   DOUBLE PRECISION :: ZooNanoNanoEgeSMOPN,NanoNanoEgeSMOPN,SyneNanoEgeSMOPN,BactNanoEgeSMOPN
   DOUBLE PRECISION :: NanoEgeNSMOP
   DOUBLE PRECISION :: ZooNanoNanoEgeSMOPP,NanoNanoEgeSMOPP,SyneNanoEgeSMOPP,BactNanoEgeSMOPP
   DOUBLE PRECISION :: NanoEgePSMOP
   DOUBLE PRECISION :: NanoNanoEgeSMOPChl,SyneNanoEgeSMOPChl
   DOUBLE PRECISION :: NanoEgeChlSMOP
   DOUBLE PRECISION :: ZooMicroMicroEgeSMOPC,DiaMicroEgeSMOPC,ZooNanoMicroEgeSMOPC,NanoMicroEgeSMOPC,SyneMicroEgeSMOPC
   DOUBLE PRECISION :: BactMicroEgeSMOPC,SMOPMicroEgeSMOPC
   DOUBLE PRECISION :: MicroEgeCSMOP
   DOUBLE PRECISION :: ZooMicroMicroEgeSMOPN,DiaMicroEgeSMOPN,ZooNanoMicroEgeSMOPN,NanoMicroEgeSMOPN,SyneMicroEgeSMOPN
   DOUBLE PRECISION :: BactMicroEgeSMOPN,SMOPMicroEgeSMOPN
   DOUBLE PRECISION :: MicroEgeNSMOP
   DOUBLE PRECISION :: ZooMicroMicroEgeSMOPP,DiaMicroEgeSMOPP,ZooNanoMicroEgeSMOPP,NanoMicroEgeSMOPP,SyneMicroEgeSMOPP
   DOUBLE PRECISION :: BactMicroEgeSMOPP,SMOPMicroEgeSMOPP
   DOUBLE PRECISION :: MicroEgePSMOP
   DOUBLE PRECISION :: DiaMicroEgeSMOPChl,NanoMicroEgeSMOPChl,SyneMicroEgeSMOPChl,SMOPMicroEgeSMOPChl
   DOUBLE PRECISION :: MicroEgeChlSMOP
   DOUBLE PRECISION :: DiaMicroEgeSMOPSi,DiaMicroEgeLMOPSi
   DOUBLE PRECISION :: MicroEgeSiSMOP,MicroEgeSiLMOP
   DOUBLE PRECISION :: ZooMicroMesoEgeSMOPC,DiaMesoEgeSMOPC,NanoMesoEgeSMOPC,LMOPMesoEgeSMOPC,SMOPMesoEgeSMOPC
   DOUBLE PRECISION :: MesoEgeCSMOP
   DOUBLE PRECISION :: ZooMicroMesoEgeSMOPN,DiaMesoEgeSMOPN,NanoMesoEgeSMOPN,LMOPMesoEgeSMOPN,SMOPMesoEgeSMOPN
   DOUBLE PRECISION :: MesoEgeNSMOP
   DOUBLE PRECISION :: ZooMicroMesoEgeSMOPP,DiaMesoEgeSMOPP,NanoMesoEgeSMOPP,LMOPMesoEgeSMOPP,SMOPMesoEgeSMOPP
   DOUBLE PRECISION :: MesoEgePSMOP
   DOUBLE PRECISION :: DiaMesoEgeSMOPChl,NanoMesoEgeSMOPChl,SMOPMesoEgeSMOPChl
   DOUBLE PRECISION :: MesoEgeChlSMOP
   DOUBLE PRECISION :: LMOPMesoEgeSMOPSi,DiaMesoEgeSMOPSi,DiaMesoEgeLMOPSi
   DOUBLE PRECISION :: MEsoEgeSiSMOP,MesoEgeSiLMOP
   DOUBLE PRECISION :: EgeCSMOP,EgeC,EgeNSMOP,EgeN,EgePSMOP,EgeP,EgeChlSMOP,EgeChl,EgeSiLMOP,EgeSiSMOP,EgeSi

   DOUBLE PRECISION :: NPZooNano,ZooNanoFoodN,ZooNanoFoodP,ZooNanoFoodC,NCNanoFood,PCNanoFood,NPNanoFood
   DOUBLE PRECISION :: ZooNanoExcNH4,ZooNanoExcPO4,ZooNanoExcCO2
   DOUBLE PRECISION :: NPZooMicro,ZooMicroFoodN,ZooMicroFoodP,ZooMicroFoodC,NCMicroFood,PCMicroFood,NPMicroFood
   DOUBLE PRECISION :: ZooMicroExcNH4,ZooMicroExcPO4,ZooMicroExcCO2
   DOUBLE PRECISION :: NPZooMeso,ZooMesoFoodN,ZooMesoFoodP,ZooMesoFoodC,NCMesoFood,PCMesoFood,NPMesoFood
   DOUBLE PRECISION :: ZooMesoExcNH4,ZooMesoExcPO4,ZooMesoExcCO2
   DOUBLE PRECISION :: ZooExcNH4,ZooExcPO4,ZooExcCO2
   DOUBLE PRECISION :: ZooNanoGrazNetGrowthEffC,ZooMicroGrazNetGrowthEffC,ZooMesoGrazNetGrowthEffC

   DOUBLE PRECISION :: ZooMesoPredC,ZooMesoPredSMOPC,ZooMesoPredLMOPC
   DOUBLE PRECISION :: ZooMesoPredN,ZooMesoPredSMOPN,ZooMesoPredLMOPN
   DOUBLE PRECISION :: ZooMesoPredP,ZooMesoPredSMOPP,ZooMesoPredLMOPP
   DOUBLE PRECISION :: ZooMicroMortC,ZooMicroMortSMOPC,ZooMicroMortLMOPC
   DOUBLE PRECISION :: ZooMicroMortN,ZooMicroMortSMOPN,ZooMicroMortLMOPN
   DOUBLE PRECISION :: ZooMicroMortP,ZooMicroMortSMOPP,ZooMicroMortLMOPP
   DOUBLE PRECISION :: ZooNanoMortC,ZooNanoMortSMOPC,ZooNanoMortLMOPC
   DOUBLE PRECISION :: ZooNanoMortN,ZooNanoMortSMOPN,ZooNanoMortLMOPN
   DOUBLE PRECISION :: ZooNanoMortP,ZooNanoMortSMOPP,ZooNanoMortLMOPP

   DOUBLE PRECISION :: LimitationBactMODC,GrowthBact,RespBact,NCMOD,PCMOD,UptBactMODC,UptBactMODN,UptBactMODP
   DOUBLE PRECISION :: NPBact,BactFoodN,BactFoodP,BactFoodC,NCBactFood,PCBactFood,NPBactFood
   DOUBLE PRECISION :: BactExcNH4,BactExcPO4,UptBactAmmoS,UptBactAmmo,UptBactPS,UptBactP,LimitationBactAmmo  &
                       ,LimitationBactP
   DOUBLE PRECISION :: TestBactP,TestBactN
   DOUBLE PRECISION :: MortBactC,MortBactN,MortBactP

   DOUBLE PRECISION :: TempfacNit, Nitrification
   DOUBLE PRECISION :: PPBSyne,PPBNano,PPBDia,MuphySyne,MuphyNano,MuphyDia,fppbSyne,fppbNano,fppbDia
   DOUBLE PRECISION :: TempfacppbSyne,TempfacppbNano,TempfacppbDia
   DOUBLE PRECISION :: ExuSyneC,ExuNanoC,ExuDiaC,ExuC
   DOUBLE PRECISION :: RespSyneC,RespNanoC,RespDiaC
   DOUBLE PRECISION :: ExuSyneNit,ExuNanoNit,ExuDiaNit,ExuSyneAmmo,ExuNanoAmmo,ExuDiaAmmo,ExuSyneN,ExuNanoN,ExuDiaN,ExuN 
   DOUBLE PRECISION :: ExuSyneP,ExuNanoP,ExuDiaP,ExuP
   DOUBLE PRECISION :: rhoChlSyne,UptNVelSyne,ChlprodSyne
   DOUBLE PRECISION :: rhoChlNano,UptNVelNano,ChlprodNano
   DOUBLE PRECISION :: rhoChlDia,UptNVelDia,ChlprodDia
   DOUBLE PRECISION :: LimitationSyneNitrate,LimitationSyneAmmo,LimitationSyneP,InhibitionAmmonium
   DOUBLE PRECISION :: LimitationNanoNitrate,LimitationNanoAmmo,LimitationNanoP
   DOUBLE PRECISION :: LimitationDiaNitrate,LimitationDiaAmmo,LimitationDiaP,LimitationDiaSi
   DOUBLE PRECISION :: UptmaxSyneNit,UptmaxSyneAmmo,UptmaxSyneP
   DOUBLE PRECISION :: UptmaxNanoNit,UptmaxNanoAmmo,UptmaxNanoP
   DOUBLE PRECISION :: UptmaxDiaNit,UptmaxDiaAmmo,UptmaxDiaP,UptmaxDiaSi
   DOUBLE PRECISION :: UptSyneNit,UptSyneAmmo,UptSyneP,UptSyneN
   DOUBLE PRECISION :: UptNanoNit,UptNanoAmmo,UptNanoP,UptNanoN
   DOUBLE PRECISION :: UptDiaNit,UptDiaAmmo,UptDiaP,UptDiaSi,UptDiaN
   DOUBLE PRECISION :: UptNit,UptAmmo,UptP
   DOUBLE PRECISION :: fUptSyne,fUptDia,fUptNano 
   DOUBLE PRECISION :: ExcSyneAmmo,ExcSynePO,ExcNanoAmmo,ExcNanoPO,ExcDiaAmmo,ExcDiaPO,ExcPO,ExcDiaSi,ExcAmmo
   
   DOUBLE PRECISION :: NCDia,NCNano,NCSyne,NCSMOP,NCLMOP,PCDia,PCNano,PCSyne,PCSMOP,PCLMOP,ChlCDia,ChlCNano  &
                        ,ChlCSyne,ChlCSMOP
   DOUBLE PRECISION :: SiCDia, NCDiaeff, PCDiaeff,SiCDiaeff,gmlDia,func
   DOUBLE PRECISION :: NCNanoeff, PCNanoeff,gmlNano
   DOUBLE PRECISION :: NCSyneeff, PCSyneeff,gmlSyne
   DOUBLE PRECISION :: ChlNDia,ChlNNano,ChlNSyne
   DOUBLE PRECISION :: NCLimSyne,PCLimSyne,NCLimNano,PCLimNano,NCLimDia,PCLimDia,SiCLimDia  
   DOUBLE PRECISION :: SiCSMOP,SiCLMOP,brol

   DOUBLE PRECISION :: MortDiaC,MortDiaN,MortDiaP,MortDiaChl,MortDiaSi
   DOUBLE PRECISION :: MortSyneC,MortSyneN,MortSyneP,MortSyneChl
   DOUBLE PRECISION :: MortNanoC,MortNanoN,MortNanoP,MortNanoChl
   DOUBLE PRECISION :: SMOPMortPhytoC,SMOPMortPhytoN,SMOPMortPhytoP,SMOPMortPhytoChl,SMOPMortPhytoSi
   DOUBLE PRECISION :: LMOPMortPhytoC,LMOPMortPhytoN,LMOPMortPhytoP,LMOPMortPhytoChl,LMOPMortPhytoSi

   DOUBLE PRECISION :: TempfacRem,RemSMOPC,RemSMOPN,RemSMOPP,RemSMOPSi,RemSMOPChl
   DOUBLE PRECISION :: RemLMOPC,RemLMOPN,RemLMOPP,RemLMOPSi,RemMODC,RemMODN,RemMODP
   DOUBLE PRECISION :: EPAR_SURF,MuphyNRSyne,MuphyNRNano,MuphyNRDia,UptmaxSyneN,UptmaxNanoN,UptmaxDiaN,ExuDiaSi,   &
                       RespSyneNit,RespSyneAmmo,RespSyneP,RespNanoNit,RespNanoAmmo,RespNanoP,RespDiaNit  &
                      ,RespDiaAmmo,RespDiaP,RespDiaSi
   DOUBLE PRECISION :: Test,Tempfaczoo, Tempfacbact 

   DOUBLE PRECISION :: IntNitrate,IntPhosphate

   DOUBLE PRECISION, EXTERNAL :: f_pproduction

   DOUBLE PRECISION :: Tempfacoxid

   DOUBLE PRECISION :: limitOxygen_oxicmin_satOxygen,limitOxygen_oduoxid_satOxygen,limitOxygen_Ammoniumoxid_satOxygen,  &
                         limitOxygen_anoxmin_inOxygen,limitOxygen_oduoxidNit_inOxygen,limitOxygen_oduoxidNit_satOxygen, &
                         limitNitrate_anoxmin_inOxygen,limitNitrate_oduoxidNit_satNit,limitIron_solform_satIron
   DOUBLE PRECISION :: respBACoxygen,anoxicrespBAC,ODUoxidOxygen,ODUoxidNitrate,AmmoniumoxidOxygen,solidODU          
   DOUBLE PRECISION :: ZooGrazHerbiC




!-----------------------------------------------------------------------------!

      UptBactP=0.
      ExcAmmo=0.

      IntNanoZoo   = 0. 
      IntMicroZoo  = 0.
      IntMesoZoo   = 0.
      IntNanoZoo   = 0. 
      IntMicroZoo  = 0.
      IntMesoZoo   = 0.
      IntNitrate   = 0.
      IntPhosphate = 0.   

!Phyto
      TotalPPBSyne         =  0.D0
      TotalPPBNano         =  0.D0
      TotalPPBDia          =  0.D0
      TotalPPB             =  0.D0

      TotalNetPPBSyne      =  0.D0
      TotalNetPPBNano      =  0.D0
      TotalNetPPBDia       =  0.D0
      TotalNetPPB          =  0.D0

      TotalUptN            =  0.D0

      TotalUptNit          =  0.D0
      TotalUptSyneNit      =  0.D0
      TotalUptDiaNit       =  0.D0
      TotalUptNanoNit      =  0.D0

      TotalUptAmmo         =  0.D0
      TotalUptSyneAmmo     =  0.D0
      TotalUptDiaAmmo      =  0.D0
      TotalUptNanoAmmo     =  0.D0
      TotalUptBactAmmo     =  0.D0

      TotalExcZooP         =  0.D0

      TotalUptBactMODC     =  0.D0
      TotalUptBactMODN     =  0.D0
      TotalUptBactMODP     =  0.D0
      TotalUptBactP        =  0.D0
      TotalBactExcNH4INT   =  0.D0
      TotalBactExcPO4INT   =  0.D0
      TotalRespBact        =  0.D0
      TotalGrowthBact      =  0.D0
      TotalMortBactC       =  0.D0


      TotalRespPhyto       =  0.D0
      TotalRespSyne        =  0.D0
      TotalRespNano        =  0.D0
      TotalRespDia         =  0.D0

      TotalNitrif          =  0.D0
      NewPPB               =  0.D0
      TotalResp            =  0.D0
      TotalRespZoo         =  0.D0 
      TotalExuC            =  0.D0
      TotalExuN            =  0.D0
      TotalExuP            =  0.D0
      
      TotalGrazC           =  0.D0
      TotalGrazCPur        =  0.D0
      TotalGrazCHerbi      =  0.D0
      TotalGrazN           =  0.D0
      TotalGrazP           =  0.D0
      
      TotalEgestionC       =  0.D0
      TotalEgestionN       =  0.D0
      TotalEgestionP       =  0.D0
      
      TotalMessyFeedC      =  0.D0
      TotalMessyFeedN      =  0.D0
      TotalMessyFeedP      =  0.D0
      
      TotalUptPhytoP       =  0.D0
      TotalUptSyneP        =  0.D0
      TotalUptNanoP        =  0.D0
      TotalUptDiaP         =  0.D0
      
      TotalRemSMOPC        =  0.D0
      TotalRemSMOPN        =  0.D0
      TotalRemSMOPP        =  0.D0
      
      TotalRemLMOPC        =  0.D0
      TotalRemLMOPN        =  0.D0
      TotalRemLMOPP        =  0.D0
      
      TotalZooGrazSMOPC    =  0.D0
      TotalZooGrazSMOPN    =  0.D0
      TotalZooGrazSMOPP    =  0.D0
      
      TotalZooGrazLMOPC    =  0.D0
      TotalZooGrazLMOPN    =  0.D0
      TotalZooGrazLMOPP    =  0.D0
      
      TotalZooGrazBact     =  0.D0
      
      TotalZooGrazPhytoC   =  0.D0
      TotalZooGrazPhytoN   =  0.D0
      TotalZooGrazPhytoP   =  0.D0
      
      
      TotalZooMesoPredSMOPC  =  0.D0
      TotalZooMesoPredSMOPN  =  0.D0
      TotalZooMesoPredSMOPP  =  0.D0
      
      TotalZooMesoPredLMOPC  =  0.D0
      TotalZooMesoPredLMOPN  =  0.D0
      TotalZooMesoPredLMOPP  =  0.D0
      
      TotalZooMicroMortSMOPC =  0.D0
      TotalZooMicroMortSMOPN =  0.D0
      TotalZooMicroMortSMOPP =  0.D0
      
      TotalZooNanoMortSMOPC  =  0.D0
      TotalZooNanoMortSMOPN  =  0.D0
      TotalZooNanoMortSMOPP  =  0.D0
      
      TotalSMOPMortPhytoC    =  0.D0
      TotalSMOPMortPhytoN    =  0.D0
      TotalSMOPMortPhytoP    =  0.D0
      
      TotalLMOPMortPhytoC    =  0.D0
      TotalLMOPMortPhytoN    =  0.D0
      TotalLMOPMortPhytoP    =  0.D0
      
      TotalMortPhytoC  =  0.D0
      TotalMortPhytoN  =  0.D0
      TotalMortPhytoP  =  0.D0
      
      TotalNitrateSurf =  0.D0
      
      IntSyne     =  0.
      IntNano     =  0.
      IntDia      =  0.
      PhytoChl    =  0.
      MaxDiaChl   =  0.
      MaxNanoChl  =  0.
      MaxSyneChl  =  0. 


      TotalBactExcNH4TOPLAYER  =  0.D0
      TotalBactExcNH4INT       =  0.D0
      TotalBactExcNH4DEEP      =  0.D0

      TotalExcZooAmmoTOPLAYER  =  0.D0
      TotalExcZooAmmoINT       =  0.D0
      TotalExcZooAmmoDEEP      =  0.D0

      TotalBactExcPO4TOPLAYER  =  0.D0
      TotalBactExcPO4INT       =  0.D0
      TotalBactExcPO4DEEP      =  0.D0

      TotalExcZooPO4TOPLAYER   =  0.D0
      TotalExcZooPO4INT        =  0.D0
      TotalExcZooPO4DEEP       =  0.D0
 
      TotalExuSiTOPLAYER  =  0.D0 
      TotalExuSiINT       =  0.D0
      TotalExuSiDEEP      =  0.D0

      TotalUptNitTOPLAYER  =  0.D0
      TotalUptNitINT       =  0.D0
      TotalUptNitDEEP      =  0.D0 
      TotalUptAmmoTOPLAYER =  0.D0
      TotalUptAmmoINT      =  0.D0
      TotalUptAmmoDEEP     =  0.D0
      TotalUptPTOPLAYER    =  0.D0
      TotalUptPINT         =  0.D0
      TotalUptPDEEP        =  0.D0
      TotalUptSiTOPLAYER   =  0.D0
      TotalUptSiINT        =  0.D0
      TotalUptSiDEEP       =  0.D0

      TotalBactUptPTOPLAYER     =  0.D0
      TotalBactUptPINT          =  0.D0
      TotalBactUptPDEEP         =  0.D0
      TotalBactUptAmmoTOPLAYER  =  0.D0
      TotalBactUptAmmoINT       =  0.D0
      TotalBactUptAmmoDEEP      =  0.D0


      TotalNitrifTOPLAYER =  0.D0
      TotalNitrifINT      =  0.D0
      TotalNitrifDEEP     =  0.D0    

      TotalRemSMOPSiTOPLAYER =  0.D0
      TotalRemSMOPSiINT      =  0.D0
      TotalRemSMOPSiDEEP     =  0.D0
      TotalRemLMOPSiTOPLAYER =  0.D0
      TotalRemLMOPSiINT      =  0.D0
      TotalRemLMOPSiDEEP     =  0.D0

      TotalNitrifCOLUMN      =  0.D0

      prefNanoBact   = 1.D0 - prefNanoZooNano   - prefNanoNano    - prefNanoSyne 
      prefMicroSMOP  = 1.D0 - prefMicroZooMicro - prefMicroDia    - prefMicroZooNano  &
                            - prefMicroNano   - prefMicroSyne  - prefMicroBact  ! Preference factor for small particulate organic matter
      prefMesoSMOP   = 1.D0 - prefMesoZooMicro  - prefMesoDia     - prefMesoLMOP   ! Preference factor for small particulate organic matter

     if(prefMicroSMOP < 0.D0) then
       write(6,*)'Error: prefMicroSMOP <0'
     endif
     if(prefMesoSMOP < 0.D0) then
       write(6,*)'Error: prefMesoSMOP <0'
     endif  

     DO I=1,NumpelagicBoxes  !->NumpelagicBoxes
       if(maskbio(I)==1) then
         PARZ = max(1d-20, PAR_Z(I)) 

         Chl(I)       = 0.D0
         NOP(I)       = 0.D0
         COP(I)       = 0.D0
         ZooNanoC(I)  =max(0.D0,ZooNanoC(I)) 
         ZooMicroC(I) =max(0.D0,ZooMicroC(I))
         ZooMesoC(I)  =max(0.D0,ZooMesoC(I)) 
         BactC(I)     =max(0.D0,BactC(I))
         DiaC(I)      =max(0.D0,DiaC(I))
         DiaN(I)      =max(0.D0,DiaN(I))
         DiaP(I)      =max(0.D0,DiaP(I))
         DiaChl(I)    =max(0.D0,DiaChl(I)) 
         DiaSi(I)     =max(0.D0,DiaSi(I))
         NanoC(I)     =max(0.D0,NanoC(I))
         NanoN(I)     =max(0.D0,NanoN(I))
         NanoP(I)     =max(0.D0,NanoP(I))
         NanoChl(I)   =max(0.D0,NanoChl(I)) 
         SyneC(I)     =max(0.D0,SyneC(I))
         SyneN(I)     =max(0.D0,SyneN(I))
         SyneP(I)     =max(0.D0,SyneP(I))
         SyneChl(I)   =max(0.D0,SyneChl(I))
         S_MOPC(I)    =max(0.D0,S_MOPC(I))
         S_MOPN(I)    =max(0.D0,S_MOPN(I))
         S_MOPP(I)    =max(0.D0,S_MOPP(I))
         S_MOPChl(I)  =max(0.D0,S_MOPChl(I)) 
         S_MOPSi(I)   =max(0.D0,S_MOPSi(I))
         L_MOPC(I)    =max(0.D0,L_MOPC(I))
         L_MOPN(I)    =max(0.D0,L_MOPN(I))
         L_MOPP(I)    =max(0.D0,L_MOPP(I))
         L_MOPSi(I)   =max(0.D0,L_MOPSi(I))
         MODC(I)      =max(0.D0,MODC(I))
         MODN(I)      =max(0.D0,MODN(I))
         MODP(I)      =max(0.D0,MODP(I))
         Nitrate(I)   =max(0.D0,Nitrate(I))
         Ammonium(I)  =max(0.D0,Ammonium(I))
         Phosphate(I) =max(0.D0,Phosphate(I))
         Silice(I)    =max(0.D0,Silice(I))
         Oxygen(I)    =max(0.D0,Oxygen(I))
         ODU(I)       =max(0.D0,ODU(I))

! The quota. 
         if (DIAC(I) > 0) then
           NCDia      = DiaN(I)   / DiaC(I) 
           PCDia      = DiaP(I)   / DiaC(I) 
           ChlCDia    = DiaChl(I)  / DiaC(I) 
           SiCDia     = DiaSi(I)   / DiaC(I) 
         else
           NcDia      = 0.D0
           PcDia      = 0.D0
           ChlCDia    = 0.D0
           SiCDia     = 0.D0
         endif

         if (NanoC(I) > 0) then
           NCNano     = NanoN(I)   /NanoC(I) 
           PCNano     = NanoP(I)   /NanoC(I) 
           ChlCNano   = NanoChl(I) /NanoC(I) 
         else
           NCNano     = 0.D0
           PCNano     = 0.D0
           ChlCNano   = 0.D0
         endif
     
         if (SYNEC(I) > 0) then
           NCSyne     = SyneN(I)   /SyneC(I) 
           PCSyne     = SyneP(I)   /SyneC(I) 
           ChlCSyne   = SyneChl(I) /SyneC(I) 
         else
           NCSyne     = 0.D0
           PCSyne     = 0.D0
           ChlCSyne   = 0.D0
         endif

         if (S_MOPC(I) > 0) then
           NCSmop     = S_MOPN(I)   / S_MOPC(I) 
           PCSmop     = S_MOPP(I)   / S_MOPC(I) 
           ChlCSmop   = S_MOPChl(I) / S_MOPC(I) 
           SiCSmop    = S_MOPSi(I)  / S_MOPC(I) 
         else
           NCSmop     = 0.D0
           PCSmop     = 0.D0
           ChlCSmop   = 0.D0
           SiCSmop    = 0.D0
         endif

         if (L_MOPC(I) > 0) then
           NCLmop     = L_MOPN(I)   / L_MOPC(I) 
           PCLmop     = L_MOPP(I)   / L_MOPC(I) 
           SiCLmop    = L_MOPSi(I)  / L_MOPC(I) 
         else
           NCLmop     = 0.D0
           PCLmop     = 0.D0
           SiCLmop    = 0.D0
         endif

         if (MODC(I) > 0) then
           NCMOD      = MODN(I)   / MODC(I) 
           PCMOD      = MODP(I)   / MODC(I) 
         else
           NCMOD      = 0.D0
           PCMOD      = 0.D0
         endif

         if (SYNEN(I) > 0) then
           ChlNSyne   = SyneChl(I)   /SyneN(I) 
         else
           ChlNSyne   = 0.D0
         endif

         if (NanoN(I) > 0) then
           ChlNNano   = NanoChl(I)   /NanoN(I) 
         else
           ChlNNano   = 0.D0
         endif


         if (DIAN(I) > 0) then
           ChlNDia    = DiaChl(I)  / DiaN(I) 
         else
           ChlNDia    = 0.D0
         endif


         Tempfaczoo   = Q10_Zoo**((Temperature(I)-TREF_Zoo)/10.)

!~~~~~~~~~~!
!  GRAZING !
!~~~~~~~~~~!


! Nano-zooplankton grazing

         Test = prefNanoBact*BactC(I)+prefNanoSyne*SyneC(I)+prefNanoZooNano*ZooNanoC(I)+prefNanoNano*NanoC(I)   

         if ( Test > 0) then 

           fZooNano = Tempfaczoo  &
              / (kZooNano*(prefNanoBact*BactC(I)+prefNanoSyne*SyneC(I)+prefNanoNano*NanoC(I) +prefNanoZooNano*ZooNanoC(I) ) +   &
                 prefNanoBact*BactC(I)**2 + prefNanoSyne*SyneC(I)**2 +prefNanoNano*NanoC(I)**2 +prefNanoZooNano*ZooNanoC(I)**2 )

         else
           fZooNano = 0.
         endif

         ZooNanoGrazZooNanoC = grazNano  * ZooNanoC(I) * ZooNanoC(I)**2.D0   * prefNanoZooNano  * fZooNano  
         ZooNanoGrazNanoC    = grazNano  * ZooNanoC(I) *    NanoC(I)**2.D0   * prefNanoNano     * fZooNano  
         ZooNanoGrazSyneC    = grazNano  * ZooNanoC(I) *    SyneC(I)**2.D0   * prefNanoSyne     * fZooNano  
         ZooNanoGrazBactC    = grazNano  * ZooNanoC(I) *    BactC(I)**2.D0   * prefNanoBact     * fZooNano  
         ZooNanoGrazC        = ZooNanoGrazZooNanoC + ZooNanoGrazNanoC + ZooNanoGrazSyneC + ZooNanoGrazBactC 

         ZooNanoGrazZooNanoN = NCZooNano * ZooNanoGrazZooNanoC
         ZooNanoGrazNanoN    = NCNano    * ZooNanoGrazNanoC
         ZooNanoGrazSyneN    = NCSyne    * ZooNanoGrazSyneC  
         ZooNanoGrazBactN    = NCBact    * ZooNanoGrazBactC
         ZooNanoGrazN        = ZooNanoGrazZooNanoN + ZooNanoGrazNanoN + ZooNanoGrazSyneN + ZooNanoGrazBactN  

         ZooNanoGrazZooNanoP = PCZooNano * ZooNanoGrazZooNanoC
         ZooNanoGrazNanoP    = PCNano    * ZooNanoGrazNanoC
         ZooNanoGrazSyneP    = PCSyne    * ZooNanoGrazSyneC 
         ZooNanoGrazBactP    = PCBact    * ZooNanoGrazBactC
         ZooNanoGrazP        = ZooNanoGrazZooNanoP + ZooNanoGrazNanoP + ZooNanoGrazSyneP + ZooNanoGrazBactP  

         ZooNanoGrazNanoChl  = ChlCNano  * ZooNanoGrazNanoC
         ZooNanoGrazSyneChl  = ChlCSyne  * ZooNanoGrazSyneC  
         ZooNanoGrazChl      = ZooNanoGrazNanoChl + ZooNanoGrazSyneChl 


! Micro-zooplankton grazing
        
         Test = prefMicroZooMicro*ZooMicroC(I) + prefMicroDia*DiaC(I) +prefMicroZooNano*ZooNanoC(I)  &
                   +prefMicroNano*NanoC(I) +prefMicroSyne*SyneC(I) + prefMicroBact*BactC(I) + prefMicroSMOP*S_MOPC(I)   

         if ( Test > 0) then 

!fZooMicro =  1.D0 / (kzooMicro*(prefMicroZooMicro*ZooMicroC(I)    + prefMicroDia*DiaC(I)    +prefMicroZooNano*ZooNanoC(I)     +prefMicroNano*NanoC(I)    +prefMicroSyne*SyneC(I)    + prefMicroBact*BactC(I)     + prefMicroSMOP*S_MOPC(I)) +   &
!                                       prefMicroZooMicro*ZooMicroC(I)**2 + prefMicroDia*DiaC(I)**2 +prefMicroZooNano*ZooNanoC(I)**2  +prefMicroNano*NanoC(I)**2 +prefMicroSyne*SyneC(I)**2 + prefMicroBact*BactC(I)**2. + prefMicroSMOP*S_MOPC(I)**2)

           fZooMicro =  Tempfaczoo & 
                   / (kzooMicro*(prefMicroZooMicro*ZooMicroC(I)    + prefMicroDia*DiaC(I)    +prefMicroZooNano*ZooNanoC(I)   &  
                    +prefMicroNano*NanoC(I)+prefMicroSyne*SyneC(I)    + prefMicroBact*BactC(I)     + prefMicroSMOP*S_MOPC(I)) +   &
                           prefMicroZooMicro*ZooMicroC(I)**2 + prefMicroDia*DiaC(I)**2 +prefMicroZooNano*ZooNanoC(I)**2  & 
                    +prefMicroNano*NanoC(I)**2 +prefMicroSyne*SyneC(I)**2 + prefMicroBact*BactC(I)**2. + prefMicroSMOP*S_MOPC(I)**2)


         else
           fZooMicro = 0.
         endif

         ZooMicroGrazZooMicroC = grazMicro * ZooMicroC(I) * ZooMicroC(I)**2.D0 * prefMicroZooMicro  * fZooMicro  
         ZooMicroGrazDiaC      = grazMicro * ZooMicroC(I) * DiaC(I)**2.D0      * prefMicroDia       * fZooMicro  
         ZooMicroGrazZooNanoC  = grazMicro * ZooMicroC(I) * ZooNanoC(I)**2.D0  * prefMicroZooNano   * fZooMicro
         ZooMicroGrazNanoC     = grazMicro * ZooMicroC(I) * NanoC(I)**2.D0     * prefMicroNano      * fZooMicro  
         ZooMicroGrazSyneC     = grazMicro * ZooMicroC(I) * SyneC(I)**2.D0     * prefMicroSyne      * fZooMicro  
         ZooMicroGrazBactC     = grazMicro * ZooMicroC(I) * BactC(I)**2.D0     * prefMicroBact      * fZooMicro  
         ZooMicroGrazSMOPC     = grazMicro * ZooMicroC(I) * S_MOPC(I)**2.D0    * prefMicroSMOP      * fZooMicro  
         ZooMicroGrazC         = ZooMicroGrazZooMicroC + ZooMicroGrazDiaC + ZooMicroGrazZooNanoC & 
                               + ZooMicroGrazNanoC + ZooMicroGrazSyneC + ZooMicroGrazBactC + ZooMicroGrazSMOPC 

         ZooMicroGrazZooMicroN = NCZooMicro * ZooMicroGrazZooMicroC
         ZooMicroGrazDiaN      = NCDia      * ZooMicroGrazDiaC
         ZooMicroGrazZooNanoN  = NCZooNano  * ZooMicroGrazZooNanoC 
         ZooMicroGrazNanoN     = NCNano     * ZooMicroGrazNanoC
         ZooMicroGrazSyneN     = NCSyne     * ZooMicroGrazSyneC  
         ZooMicroGrazBactN     = NCBact     * ZooMicroGrazBactC  
         ZooMicroGrazSMOPN     = NCSMOP     * ZooMicroGrazSMOPC
         ZooMicroGrazN         = ZooMicroGrazZooMicroN + ZooMicroGrazDiaN + ZooMicroGrazZooNanoN & 
                               + ZooMicroGrazNanoN + ZooMicroGrazSyneN + ZooMicrograzBactN + ZooMicroGrazSMOPN 

         ZooMicroGrazZooMicroP = PCZooMicro * ZooMicroGrazZooMicroC
         ZooMicroGrazDiaP      = PCDia      * ZooMicroGrazDiaC
         ZooMicroGrazZooNanoP  = PCZooNano  * ZooMicroGrazZooNanoC 
         ZooMicroGrazNanoP     = PCNano     * ZooMicroGrazNanoC
         ZooMicroGrazSyneP     = PCSyne     * ZooMicroGrazSyneC 
         ZooMicroGrazBActP     = PCBact     * ZooMicroGrazBactC
         ZooMicroGrazSMOPP     = PCSMOP     * ZooMicroGrazSMOPC
         ZooMicroGrazP         = ZooMicroGrazZooMicroP + ZooMicroGrazDiaP + ZooMicroGrazZooNanoP & 
                              + ZooMicroGrazNanoP + ZooMicroGrazSyneP + ZooMicroGrazBactP + ZooMicroGrazSMOPP 

         ZooMicroGrazDiaChl   = ChlCDia   * ZooMicroGrazDiaC
         ZooMicroGrazNanoChl  = ChlCNano  * ZooMicroGrazNanoC
         ZooMicroGrazSyneChl  = ChlCSyne  * ZooMicroGrazSyneC  
         ZooMicroGrazSMOPChl  = ChlCSmop  * ZooMicroGrazSMOPC
         ZooMicroGrazChl      = ZooMicroGrazDiaChl + ZooMicroGrazNanoChl + ZooMicroGrazSyneChl + ZooMicroGrazSMOPChl

         ZooMicroGrazDiaSi   = SiCDia  * ZooMicroGrazDiaC
         ZooMicroGrazSi      = ZooMicroGrazDiaSi 


! Meso-zooplankton grazing

         Test = prefMesoZooMicro*ZooMicroC(I)+prefMesoDia*DiaC(I)+ prefMesoSMOP*S_MOPC(I)    + prefMesoLMOP*L_MOPC(I)
!        Test = prefMesoZooMicro*ZooMicroC(I)    +prefMesoDia*DiaC(I)    + prefMesoNano*NanoC(I)     + prefMesoSMOP*S_MOPC(I)    + prefMesoLMOP*L_MOPC(I)


         if ( Test > 0) then 
   
           fZooMeso  =  Tempfaczoo & 
             / (kzooMeso*(prefMesoZooMicro*ZooMicroC(I)+prefMesoDia*DiaC(I)+prefMesoSMOP*S_MOPC(I)+prefMesoLMOP*L_MOPC(I))+  &
                prefMesoZooMicro*ZooMicroC(I)**2 +prefMesoDia*DiaC(I)**2 + prefMesoSMOP*S_MOPC(I)**2 + prefMesoLMOP*L_MOPC(I)**2)
!fZooMeso  =  Tempfaczoo / (kzooMeso*(prefMesoZooMicro*ZooMicroC(I)    +prefMesoDia*DiaC(I)    + prefMesoNanoC*NanoC(I)    + prefMesoSMOP*S_MOPC(I)    + prefMesoLMOP*L_MOPC(I)) +   &
!                                       prefMesoZooMicro*ZooMicroC(I)**2 +prefMesoDia*DiaC(I)**2 +prefMesoNano*NanoC(I)**2 + prefMesoSMOP*S_MOPC(I)**2 + prefMesoLMOP*L_MOPC(I)**2)


         else
           fZooMeso = 0.
         endif

         ZooMesoGrazZooMicroC = grazMeso  * ZooMesoC(I) * ZooMicroC(I)**2.D0 * prefMesoZooMicro  * fZooMeso  
         ZooMesoGrazDiaC      = grazMeso  * ZooMesoC(I) *      DiaC(I)**2.D0 * prefMesoDia       * fZooMeso  
!        ZooMesoGrazNanoC     = grazMeso  * ZooMesoC(I) *     NanoC(I)**2.D0 * prefMesoNano      * fZooMeso  
         ZooMesoGrazSMOPC     = grazMeso  * ZooMesoC(I) *    S_MOPC(I)**2.D0 * prefMesoSMOP      * fZooMeso  
         ZooMesoGrazLMOPC     = grazMeso  * ZooMesoC(I) *    L_MOPC(I)**2.D0 * prefMesoLMOP      * fZooMeso
         ZooMesoGrazC         = ZooMesoGrazZooMicroC + ZooMesoGrazDiaC & 
                              + ZooMesoGrazSMOPC + ZooMesoGrazLMOPC !+ ZooMesoGrazNanoC

         ZooMesoGrazZooMicroN = NCZooMicro * ZooMesoGrazZooMicroC
         ZooMesoGrazDiaN      = NCDia      * ZooMesoGrazDiaC
!        ZooMesoGrazNanoN     = NCNano     * ZooMesoGrazNanoC
         ZooMesoGrazSMOPN     = NCSMOP     * ZooMesoGrazSMOPC  
         ZooMesoGrazLMOPN     = NCLMOP     * ZooMesoGrazLMOPC
         ZooMesoGrazN         = ZooMesoGrazZooMicroN + ZooMesoGrazDiaN & 
                              + ZooMesoGrazSMOPN + ZooMesoGrazLMOPN !+ ZooMesoGrazNanoN

         ZooMesoGrazZooMicroP = PCZooMicro * ZooMesoGrazZooMicroC
         ZooMesoGrazDiaP      = PCDia      * ZooMesoGrazDiaC
!        ZooMesoGrazNanoP     = PCNano     * ZooMesoGrazNanoC
         ZooMesoGrazSMOPP     = PCSMOP     * ZooMesoGrazSMOPC  
         ZooMesoGrazLMOPP     = PCLMOP     * ZooMesoGrazLMOPC
         ZooMesoGrazP         = ZooMesoGrazZooMicroP + ZooMesoGrazDiaP & 
                              + ZooMesoGrazSMOPP + ZooMesoGrazLMOPP ! + ZooMesoGrazNanoP

         ZooMesoGrazDiaChl    = ChlCDia      * ZooMesoGrazDiaC
!        ZooMesoGrazNanoChl   = ChlCNano     * ZooMesoGrazNanoC
         ZooMesoGrazSMOPChl   = ChlCSMOP     * ZooMesoGrazSMOPC  
         ZooMesoGrazChl       = ZooMesoGrazDiaChl + ZooMesoGrazSMOPChl !+ ZooMesoGrazNanoChl

         ZooMesoGrazDiaSi  = SiCDia  * ZooMesoGrazDiaC
         ZooMesoGrazSi     = ZooMesoGrazDiaSi 


! Total grazing

         ZooGrazZooMicroC = ZooMesoGrazZooMicroC + ZooMicroGrazZooMicroC
         ZooGrazDiaC      = ZooMesoGrazDiaC      + ZooMicroGrazDiaC
         ZooGrazZooNanoC  =                      + ZooMicroGrazZooNanoC   + ZooNanoGrazZooNanoC 
         ZooGrazNanoC     =                        ZooMicroGrazNanoC      + ZooNanoGrazNanoC
         ZooGrazSyneC     =                        ZooMicroGrazSyneC      + ZooNanoGrazSyneC
         ZooGrazBactC     =                        ZooMicroGrazBactC      + ZooNanoGrazBactC 
         ZooGrazSMOPC     = ZooMesoGrazSMOPC     + ZooMicroGrazSMOPC     
         ZooGrazLMOPC     = ZooMesoGrazLMOPC
         ZooGrazC         = ZooGrazZooMicroC + ZooGrazDiaC + ZooGrazZooNanoC + ZooGrazNanoC & 
                           + ZooGrazSyneC + ZooGrazBactC + ZooGrazSMOPC + ZooGrazLMOPC 
! claude 
         ZooGrazHerbiC = ZooGrazDiaC + ZooGrazNanoC + ZooGrazSyneC


 
         ZooGrazZooMicroN = ZooMesoGrazZooMicroN + ZooMicroGrazZooMicroN
         ZooGrazDiaN      = ZooMesoGrazDiaN      + ZooMicroGrazDiaN 
         ZooGrazZooNanoN  =                      + ZooMicroGrazZooNanoN   + ZooNanoGrazZooNanoN 
         ZooGrazNanoN     =                        ZooMicroGrazNanoN      + ZooNanoGrazNanoN
         ZooGrazSyneN     =                        ZooMicroGrazSyneN      + ZooNanoGrazSyneN
         ZooGrazBactN     =                        ZooMicroGrazBactN      + ZooNanoGrazBactN 
         ZooGrazSMOPN     = ZooMesoGrazSMOPN     + ZooMicroGrazSMOPN  
         ZooGrazLMOPN     = ZooMesoGrazLMOPN
         ZooGrazN         = ZooGrazZooMicroN + ZooGrazDiaN + ZooGrazZooNanoN + ZooGrazNanoN & 
                           + ZooGrazSyneN + ZooGrazBactN + ZooGrazSMOPN + ZooGrazLMOPN 

         ZooGrazZooMicroP = ZooMesoGrazZooMicroP + ZooMicroGrazZooMicroP
         ZooGrazDiaP      = ZooMesoGrazDiaP      + ZooMicroGrazDiaP
         ZooGrazZooNanoP  =                      + ZooMicroGrazZooNanoP   + ZooNanoGrazZooNanoP 
         ZooGrazNanoP     =                        ZooMicroGrazNanoP      + ZooNanoGrazNanoP
         ZooGrazSyneP     =                        ZooMicroGrazSyneP      + ZooNanoGrazSyneP
         ZooGrazBactP     =                        ZooMicroGrazBactP      + ZooNanoGrazBactP
         ZooGrazSMOPP     = ZooMesoGrazSMOPP     + ZooMicroGrazSMOPP  
         ZooGrazLMOPP     = ZooMesoGrazLMOPP
         ZooGrazP         = ZooGrazZooMicroP + ZooGrazDiaP + ZooGrazZooNanoP + ZooGrazNanoP & 
                           + ZooGrazSyneP + ZooGrazSMOPP + ZooGrazLMOPP  + ZooGrazBactP 

         ZooGrazDiaChl    = ZooMesoGrazDiaChl      + ZooMicroGrazDiaChl
         ZooGrazNanoChl   =                          ZooMicroGrazNanoChl  + ZooNanoGrazNanoChl
         ZooGrazSyneChl   =                          ZooMicroGrazSyneChl  + ZooNanoGrazSyneChl
         ZooGrazSMOPChl   = ZooMesoGrazSMOPChl     + ZooMicroGrazSMOPChl  
         ZooGrazChl       = ZooGrazDiaChl + ZooGrazNanoChl + ZooGrazSyneChl + ZooGrazSMOPChl 

         ZooGrazDiaSi   = ZooMesoGrazDiaSi   +  ZooMicroGrazDiaSi
         ZooGrazSi      = ZooGrazDiaSi 


!~~~~~~~~~~~~~~~!
! MESSY FEEDING !
!~~~~~~~~~~~~~~~!

! Messy feeding of Nanozooplankton

        
         ZooNanoNanoMessyFeedMODC =  MessyfeedingfracC  *  ZooNanoGrazZooNanoC  
         NanoNanoMessyFeedMODC    =  MessyfeedingfracC  *  ZooNanoGrazNanoC  
         SyneNanoMessyFeedMODC    =  MessyfeedingfracC  *  ZooNanoGrazSyneC
         BactNanoMessyFeedMODC    =  MessyfeedingfracC  *  ZooNanoGrazBactC 
         NanoMessyFeedCMOD = ZooNanoNanoMessyFeedMODC  + NanoNanoMessyFeedMODC & 
                           + SyneNanoMessyFeedMODC  + BactNanoMessyFeedMODC

         ZooNanoNanoMessyFeedMODN = MessyfeedingfracN  *   ZooNanoGrazZooNanoN 
         NanoNanoMessyFeedMODN    = MessyfeedingfracN  *   ZooNanoGrazNanoN  
         SyneNanoMessyFeedMODN    = MessyfeedingfracN  *   ZooNanoGrazSyneN
         BactNanoMessyFeedMODN    = MessyfeedingfracN  *   ZooNanoGrazBactN 
         NanoMessyFeedNMOD = ZooNanoNanoMessyFeedMODN  + NanoNanoMessyFeedMODN & 
                           + SyneNanoMessyFeedMODN + BactNanoMessyFeedMODN

         ZooNanoNanoMessyFeedMODP = MessyfeedingfracP   * ZooNanoGrazZooNanoP 
         NanoNanoMessyFeedMODP    = MessyfeedingfracP   * ZooNanoGrazNanoP  
         SyneNanoMessyFeedMODP    = MessyfeedingfracP   * ZooNanoGrazSyneP
         BactNanoMessyFeedMODP    = MessyfeedingfracP   * ZooNanoGrazBactP 
         NanoMessyFeedPMOD = ZooNanoNanoMessyFeedMODP  + NanoNanoMessyFeedMODP &
                           + SyneNanoMessyFeedMODP + BactNanoMessyFeedMODP


! Messy feeding of Microzooplankton

        
         ZooMicroMicroMessyFeedMODC =  MessyfeedingfracC  *  ZooMicroGrazZooMicroC  
         DiaMicroMessyFeedMODC      =  MessyfeedingfracC  *  ZooMicroGrazDiaC  
         ZooNanoMicroMessyFeedMODC  =  MessyfeedingfracC  *  ZooMicroGrazZooNanoC  
         NanoMicroMessyFeedMODC     =  MessyfeedingfracC  *  ZooMicroGrazNanoC  
         SyneMicroMessyFeedMODC     =  MessyfeedingfracC  *  ZooMicroGrazSyneC
         BactMicroMessyFeedMODC     =  MessyfeedingfracC  *  ZooMicroGrazBactC 
         SMOPMicroMessyFeedMODC     =  MessyfeedingfracC  *  ZooMicroGrazSMOPC 
         MicroMessyFeedCMOD = ZooMicroMicroMessyFeedMODC  + DiaMicroMessyFeedMODC  + ZooNanoMicroMessyFeedMODC & 
                            + NanoMicroMessyFeedMODC  + SyneMicroMessyFeedMODC + BactMicroMessyFeedMODC  + SMOPMicroMessyFeedMODC

         ZooMicroMicroMessyFeedMODN = MessyfeedingfracN  *   ZooMicroGrazZooMicroN 
         DiaMicroMessyFeedMODN      = MessyfeedingfracN  *   ZooMicroGrazDiaN  
         ZooNanoMicroMessyFeedMODN  = MessyfeedingfracN  *   ZooMicroGrazZooNanoN  
         NanoMicroMessyFeedMODN     = MessyfeedingfracN  *   ZooMicroGrazNanoN  
         SyneMicroMessyFeedMODN     = MessyfeedingfracN  *   ZooMicroGrazSyneN
         BactMicroMessyFeedMODN     = MessyfeedingfracN  *   ZooMicroGrazBactN 
         SMOPMicroMessyFeedMODN     = MessyfeedingfracN  *   ZooMicroGrazSMOPN 
         MicroMessyFeedNMOD = ZooMicroMicroMessyFeedMODN + DiaMicroMessyFeedMODN + ZooNanoMicroMessyFeedMODN & 
                            + NanoMicroMessyFeedMODN + SyneMicroMessyFeedMODN + BactMicroMessyFeedMODN + SMOPMicroMessyFeedMODN
 

         ZooMicroMicroMessyFeedMODP = MessyfeedingfracP   * ZooMicroGrazZooMicroP 
         DiaMicroMessyFeedMODP      = MessyfeedingfracP   * ZooMicroGrazDiaP  
         ZooNanoMicroMessyFeedMODP  = MessyfeedingfracP   * ZooMicroGrazZooNanoP  
         NanoMicroMessyFeedMODP     = MessyfeedingfracP   * ZooMicroGrazNanoP  
         SyneMicroMessyFeedMODP     = MessyfeedingfracP   * ZooMicroGrazSyneP
         BactMicroMessyFeedMODP     = MessyfeedingfracP   * ZooMicroGrazBactP 
         SMOPMicroMessyFeedMODP     = MessyfeedingfracP   * ZooMicroGrazSMOPP 
         MicroMessyFeedPMOD = ZooMicroMicroMessyFeedMODP  + DiaMicroMessyFeedMODP + ZooNanoMicroMessyFeedMODP & 
                            + NanoMicroMessyFeedMODP + SyneMicroMessyFeedMODP + BactMicroMessyFeedMODP + SMOPMicroMessyFeedMODP





! Messy feeding of MesoZooplankton


         ZooMicroMesoMessyFeedMODC  = MessyfeedingfracC   * ZooMesoGrazZooMicroC  
         DiaMesoMessyFeedMODC       = MessyfeedingfracC   * ZooMesoGrazDiaC  
!        NanoMesoMessyFeedMODC      = MessyfeedingfracC   * ZooMesoGrazNanoC  
         SMOPMesoMessyFeedMODC      = MessyfeedingfracC   * ZooMesoGrazSMOPC 
         LMOPMesoMessyFeedMODC      = MessyfeedingfracC   * ZooMesoGrazLMOPC
         MesoMessyFeedCMOD  = ZooMicroMesoMessyFeedMODC  + DiaMesoMessyFeedMODC  & 
                            + LMOPMesoMessyFeedMODC  + SMOPMesoMessyFeedMODC !+ NanoMesoMessyFeedMODC


         ZooMicroMesoMessyFeedMODN = MessyfeedingfracN   * ZooMesoGrazZooMicroN  
         DiaMesoMessyFeedMODN      = MessyfeedingfracN   * ZooMesoGrazDiaN  
!        NanoMesoMessyFeedMODN     = MessyfeedingfracN   * ZooMesoGrazNanoN  
         SMOPMesoMessyFeedMODN     = MessyfeedingfracN   * ZooMesoGrazSMOPN 
         LMOPMesoMessyFeedMODN     = MessyfeedingfracN   * ZooMesoGrazLMOPN
         MesoMessyFeedNMOD  = ZooMicroMesoMessyFeedMODN  + DiaMesoMessyFeedMODN  &
                            + LMOPMesoMessyFeedMODN + SMOPMesoMessyFeedMODN !+ NanoMesoMessyFeedMODN
        
        
         ZooMicroMesoMessyFeedMODP  = MessyfeedingfracP    * ZooMesoGrazZooMicroP  
         DiaMesoMessyFeedMODP       = MessyfeedingfracP    * ZooMesoGrazDiaP  
!        NanoMesoMessyFeedMODP      = MessyfeedingfracP    * ZooMesoGrazNanoP  
         SMOPMesoMessyFeedMODP      = MessyfeedingfracP    * ZooMesoGrazSMOPP 
         LMOPMesoMessyFeedMODP      = MessyfeedingfracP    * ZooMesoGrazLMOPP
         MesoMessyFeedPMOD  = ZooMicroMesoMessyFeedMODP  + DiaMesoMessyFeedMODP  & 
                            + LMOPMesoMessyFeedMODP + SMOPMesoMessyFeedMODP !+ NanoMesoMessyFeedMODP
       

! Total Messy feeding

         MessyFeedCMOD        = NanoMessyFeedCMOD + MicroMessyFeedCMOD  + MesoMessyFeedCMOD 
         MessyFeedNMOD        = NanoMessyFeedNMOD + MicroMessyFeedNMOD  + MesoMessyFeedNMOD 
         MessyFeedPMOD        = NanoMessyFeedPMOD + MicroMessyFeedPMOD  + MesoMessyFeedPMOD 
        

!~~~~~~~~~~!
! EGESTION !
!~~~~~~~~~~!

! Egestion of NanoZooplankton
 
         ZooNanoNanoEgeSMOPC = (1. - NanoAssEffZooNano)   *  (1. - MessyFeedingfracC) * ZooNanoGrazZooNanoC  
         NanoNanoEgeSMOPC    = (1. - NanoAssEffNano)      *  (1. - MessyFeedingfracC) * ZooNanoGrazNanoC  
         SyneNanoEgeSMOPC    = (1. - NanoAssEffSyne)      *  (1. - MessyFeedingfracC) * ZooNanoGrazSyneC
         BactNanoEgeSMOPC    = (1. - NanoAssEffBact)      *  (1. - MessyFeedingfracC) * ZooNanoGrazBactC              
         NanoEgeCSMOP = ZooNanoNanoEgeSMOPC  + NanoNanoEgeSMOPC + SyneNanoEgeSMOPC + BactNanoEgeSMOPC


         ZooNanoNanoEgeSMOPN = (1. - NanoAssEffZooNano)   * (1. - MessyFeedingfracN)  * ZooNanoGrazZooNanoN  
         NanoNanoEgeSMOPN    = (1. - NanoAssEffNano)      * (1. - MessyFeedingfracN)  * ZooNanoGrazNanoN  
         SyneNanoEgeSMOPN    = (1. - NanoAssEffSyne)      * (1. - MessyFeedingfracN)  * ZooNanoGrazSyneN 
         BactNanoEgeSMOPN    = (1. - NanoAssEffBact)      * (1. - MessyFeedingfracN)  * ZooNanoGrazBactN  
         NanoEgeNSMOP = ZooNanoNanoEgeSMOPN  + NanoNanoEgeSMOPN + SyneNanoEgeSMOPN + BactNanoEgeSMOPN


         ZooNanoNanoEgeSMOPP  = (1. - NanoAssEffZooNano)   * (1. - MessyFeedingfracP)  * ZooNanoGrazZooNanoP  
         NanoNanoEgeSMOPP     = (1. - NanoAssEffNano)      * (1. - MessyFeedingfracP)  * ZooNanoGrazNanoP  
         SyneNanoEgeSMOPP     = (1. - NanoAssEffSyne)      * (1. - MessyFeedingfracP)  * ZooNanoGrazSyneP 
         BactNanoEgeSMOPP     = (1. - NanoAssEffBact)      * (1. - MessyFeedingfracP)  * ZooNanoGrazBactP 
         NanoEgePSMOP = ZooNanoNanoEgeSMOPP  + NanoNanoEgeSMOPP + SyneNanoEgeSMOPP + BactNanoEgeSMOPP


! pourquoi pas 1 au lieu de 1-AssEff
         NanoNanoEgeSMOPChl =  ZooNanoGrazNanoChl  
         SyneNanoEgeSMOPChl =  ZooNanoGrazSyneChl 

         NanoEgeChlSMOP     =  NanoNanoEgeSMOPChl + SyneNanoEgeSMOPChl 


! Egestion of MicroZooplankton
 
         ZooMicroMicroEgeSMOPC = (1. - MicroAssEffZooMicro)  *  (1. - MessyFeedingfracC) * ZooMicroGrazZooMicroC  
         DiaMicroEgeSMOPC      = (1. - MicroAssEffDia)       *  (1. - MessyFeedingfracC) * ZooMicroGrazDiaC  
         ZooNanoMicroEgeSMOPC  = (1. - MicroAssEffZooNano)   *  (1. - MessyFeedingfracC) * ZooMicroGrazZooNanoC  
         NanoMicroEgeSMOPC     = (1. - MicroAssEffNano)      *  (1. - MessyFeedingfracC) * ZooMicroGrazNanoC  
         SyneMicroEgeSMOPC     = (1. - MicroAssEffSyne)      *  (1. - MessyFeedingfracC) * ZooMicroGrazSyneC
         BactMicroEgeSMOPC     = (1. - MicroAssEffBact)      *  (1. - MessyFeedingfracC) * ZooMicroGrazBactC            
         SMOPMicroEgeSMOPC     = (1. - MicroAssEffSMOP)      *  (1. - MessyFeedingfracC) * ZooMicroGrazSMOPC            

         MicroEgeCSMOP       = ZooMicroMicroEgeSMOPC + DiaMicroEgeSMOPC + ZooNanoMicroEgeSMOPC & 
                             + NanoMicroEgeSMOPC + SyneMicroEgeSMOPC & 
                             + BactMicroEgeSMOPC + SMOPMicroEgeSMOPC


         ZooMicroMicroEgeSMOPN = (1. - MicroAssEffZooMicro)  * (1. - MessyFeedingfracN)  * ZooMicroGrazZooMicroN  
         DiaMicroEgeSMOPN      = (1. - MicroAssEffDia)       * (1. - MessyFeedingfracN)  * ZooMicroGrazDiaN  
         ZooNanoMicroEgeSMOPN  = (1. - MicroAssEffZooNano)   * (1. - MessyFeedingfracN)  * ZooMicroGrazZooNanoN  
         NanoMicroEgeSMOPN     = (1. - MicroAssEffNano)      * (1. - MessyFeedingfracN)  * ZooMicroGrazNanoN  
         SyneMicroEgeSMOPN     = (1. - MicroAssEffSyne)      * (1. - MessyFeedingfracN)  * ZooMicroGrazSyneN 
         BactMicroEgeSMOPN     = (1. - MicroAssEffBact)      * (1. - MessyFeedingfracN)  * ZooMicroGrazBactN            
         SMOPMicroEgeSMOPN     = (1. - MicroAssEffSMOP)      * (1. - MessyFeedingfracN)  * ZooMicroGrazSMOPN            

         MicroEgeNSMOP = ZooMicroMicroEgeSMOPN + DiaMicroEgeSMOPN + ZooNanoMicroEgeSMOPN &
                       + NanoMicroEgeSMOPN + SyneMicroEgeSMOPN & 
                       + BactMicroEgeSMOPN + SMOPMicroEgeSMOPN


         ZooMicroMicroEgeSMOPP = (1. - MicroAssEffZooMicro)  * (1. - MessyFeedingfracP)  * ZooMicroGrazZooMicroP  
         DiaMicroEgeSMOPP      = (1. - MicroAssEffDia)       * (1. - MessyFeedingfracP)  * ZooMicroGrazDiaP  
         ZooNanoMicroEgeSMOPP  = (1. - MicroAssEffZooNano)   * (1. - MessyFeedingfracP)  * ZooMicroGrazZooNanoP  
         NanoMicroEgeSMOPP     = (1. - MicroAssEffNano)      * (1. - MessyFeedingfracP)  * ZooMicroGrazNanoP  
         SyneMicroEgeSMOPP     = (1. - MicroAssEffSyne)      * (1. - MessyFeedingfracP)  * ZooMicroGrazSyneP 
         BactMicroEgeSMOPP     = (1. - MicroAssEffBact)      * (1. - MessyFeedingfracP)  * ZooMicroGrazBactP 
         SMOPMicroEgeSMOPP     = (1. - MicroAssEffSMOP)      * (1. - MessyFeedingfracP)  * ZooMicroGrazSMOPP            

         MicroEgePSMOP       = ZooMicroMicroEgeSMOPP + DiaMicroEgeSMOPP + ZooNanoMicroEgeSMOPP  + NanoMicroEgeSMOPP & 
                  + SyneMicroEgeSMOPP + BactMicroEgeSMOPP + SMOPMicroEgeSMOPP


! pourquoi pas 1 au lieu de 1-AssEff
         DiaMicroEgeSMOPChl  =     ZooMicroGrazDiaChl  
         NanoMicroEgeSMOPChl =     ZooMicroGrazNanoChl  
         SyneMicroEgeSMOPChl =     ZooMicroGrazSyneChl 
         SMOPMicroEgeSMOPChl =     ZooMicroGrazSMOPChl 

         MicroEgeChlSMOP     = DiaMicroEgeSMOPChl + NanoMicroEgeSMOPChl + SyneMicroEgeSMOPChl + SMOPMicroEgeSMOPChl
         
         DiaMicroEgeSMOPSi   =        DiaPropMOPSi  *   ZooMicroGrazDiaSi 
         DiaMicroEgeLMOPSi   =  (1. - DiaPropMOPSi) *   ZooMicroGrazDiaSi            

         MicroEgeSiSMOP      = DiaMicroEgeSMOPSi 
         MicroEgeSiLMOP      = DiaMicroEgeLMOPSi


! Egestion of MesoZooplankton

         ZooMicroMesoEgeSMOPC = (1. - MesoAssEffZooMicro)  *  (1. - MessyFeedingfracC)  * ZooMesoGrazZooMicroC  
         DiaMesoEgeSMOPC = (1. - MesoAssEffDia)       *  (1. - MessyFeedingfracC)  * ZooMesoGrazDiaC  
!        NanoMesoEgeSMOPC = (1. - MesoAssEffNano)      *  (1. - MessyFeedingfracC)  * ZooMesoGrazNanoC  
         SMOPMesoEgeSMOPC = (1. - MesoAssEffSMOP)      *  (1. - MessyFeedingfracC)  * ZooMesoGrazSMOPC
         LMOPMesoEgeSMOPC = (1. - MesoAssEffLMOP)      *  (1. - MessyFeedingfracC)  * ZooMesoGrazLMOPC
            
         MesoEgeCSMOP       = ZooMicroMesoEgeSMOPC + DiaMesoEgeSMOPC + LMOPMesoEgeSMOPC + SMOPMesoEgeSMOPC !+ NanoMesoEgeSMOPC

         ZooMicroMesoEgeSMOPN = (1. - MesoAssEffZooMicro)  *  (1. - MessyFeedingfracN) * ZooMesoGrazZooMicroN  
         DiaMesoEgeSMOPN = (1. - MesoAssEffDia)       *  (1. - MessyFeedingfracN) * ZooMesoGrazDiaN  
!        NanoMesoEgeSMOPN = (1. - MesoAssEffNano)      *  (1. - MessyFeedingfracN) * ZooMesoGrazNanoN  
         SMOPMesoEgeSMOPN = (1. - MesoAssEffSMOP)      *  (1. - MessyFeedingfracN) * ZooMesoGrazSMOPN
         LMOPMesoEgeSMOPN = (1. - MesoAssEffLMOP)      *  (1. - MessyFeedingfracN) * ZooMesoGrazLMOPN

         MesoEgeNSMOP      = ZooMicroMesoEgeSMOPN + DiaMesoEgeSMOPN + LMOPMesoEgeSMOPN + SMOPMesoEgeSMOPN !+ NanoMesoEgeSMOPN

         ZooMicroMesoEgeSMOPP = (1. - MesoAssEffZooMicro)  *  (1. - MessyFeedingfracP) * ZooMesoGrazZooMicroP  
         DiaMesoEgeSMOPP = (1. - MesoAssEffDia)       *  (1. - MessyFeedingfracP) * ZooMesoGrazDiaP  
!        NanoMesoEgeSMOPP = (1. - MesoAssEffNano)      *  (1. - MessyFeedingfracP) * ZooMesoGrazNanoP  
         SMOPMesoEgeSMOPP = (1. - MesoAssEffSMOP)      *  (1. - MessyFeedingfracP) * ZooMesoGrazSMOPP
         LMOPMesoEgeSMOPP = (1. - MesoAssEffLMOP)      *  (1. - MessyFeedingfracP) * ZooMesoGrazLMOPP

         MesoEgePSMOP      = ZooMicroMesoEgeSMOPP + DiaMesoEgeSMOPP + LMOPMesoEgeSMOPP + SMOPMesoEgeSMOPP !+ NanoMesoEgeSMOPP



! pourquoi pas 1 au lieu de 1-AssEff
         DiaMesoEgeSMOPChl  =  ZooMesoGrazDiaChl  
!        NanoMesoEgeSMOPChl =  ZooMesoGrazNanoChl  
         SMOPMesoEgeSMOPChl =  ZooMesoGrazSMOPChl  


         MesoEgeChlSMOP      =  DiaMesoEgeSMOPChl + SMOPMesoEgeSMOPChl !+ NanoMesoEgeSMOPChl
 

         DiaMesoEgeSMOPSi   =        DiaPropMOPSi  *   ZooMesoGrazDiaSi 
         DiaMesoEgeLMOPSi   =  (1. - DiaPropMOPSi) *   ZooMesoGrazDiaSi            

         MesoEgeSiSMOP      = DiaMesoEgeSMOPSi 
         MesoEgeSiLMOP      = DiaMesoEgeLMOPSi
 
 
! Total Egestion

         EgeCSMOP       = MicroEgeCSMOP   + MesoEgeCSMOP   + NanoEgeCSMOP
         EgeNSMOP       = MicroEgeNSMOP   + MesoEgeNSMOP   + NanoEgeNSMOP
         EgePSMOP       = MicroEgePSMOP   + MesoEgePSMOP   + NanoEgePSMOP 
         EgeChlSMOP     = MicroEgeChlSMOP + MesoEgeChlSMOP + NanoEgeChlSMOP
         EgeSiSMOP      = MicroEgeSiSMOP  + MesoEgeSiSMOP 
         EgeSiLMOP      = MicroEgeSiLMOP  + MesoEgeSiLMOP 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! RESPIRATION & EXCRETION OF DISSOLVED INORGANIC MATTER  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! Activity Respiration

         ZooNanoGrazNetGrowthEffC = kcNetGrowthEff * ( ZooNanoGrazC - NanoEgeCSMOP - NanoMessyFeedCMOD )    



! Stoichiometric fluxes

         if(abs(NCZooNano)<1.e-2.or.abs(NCZooNano)>10.or.abs(PCZooNano)<1.e-2.or.abs(PCZooNano)>10.) &
           write(30+par%rank,*)par%rank,NCZooNano,PCZooNano
         NPZooNano = NCZooNano / PCZooNano

         ZooNanoFoodC = ZooNanoGrazNetGrowthEffC  
!        ZooNanoFoodC = ZooNanoGrazNetGrowthEffC - NanoEgeCSMOP - NanoMessyFeedCMOD
         ZooNanoFoodN = ZooNanoGrazN             - NanoEgeNSMOP - NanoMessyFeedNMOD
         ZooNanoFoodP = ZooNanoGrazP             - NanoEgePSMOP - NanoMessyFeedPMOD
 

         if ( ZooNanoFoodC > 0.D0) then
           NCNanoFood =  ZooNanoFoodN  /  ZooNanoFoodC 
           PCNanoFood =  ZooNanoFoodP  /  ZooNanoFoodC 
         else 
           NCNanoFood = 0.
           PCNanoFood = 0.
         endif


         if ( ZooNanoFoodP > 0.D0) then
           NPNanoFood =  ZooNanoFoodN  /  ZooNanoFoodP
         else 
           NPNanoFood = 0.
         endif


         if (NCNanoFood <= NCZooNano.and.PCNanoFood <= PCZooNano) then

           if ( NPNanoFood > NPZooNano ) then
 
             ZooNanoExcNH4 =  ZooNanoFoodN - ZooNanoFoodP * NPZooNano               
             ZooNanoExcCO2 =  ZooNanoFoodC - ZooNanoFoodP / PCZooNano
             ZooNanoExcPO4 =  0.D0

           elseif ( NPNanoFood < NPZooNano ) then

             ZooNanoExcPO4 =  ZooNanoFoodP - ZooNanoFoodN / NPZooNano               
             ZooNanoExcCO2 =  ZooNanoFoodC - ZooNanoFoodN / NCZooNano
             ZooNanoExcNH4 =  0.D0

          endif

         elseif (NCNanoFood > NCZooNano.and.PCNanoFood > PCZooNano) then 

           ZooNanoExcCO2 =  0.D0
           ZooNanoExcNH4 =  ZooNanoFoodN - ZooNanoFoodC * NCZooNano
           ZooNanoExcPO4 =  ZooNanoFoodP - ZooNanoFoodC * PCZooNano
       
         elseif (NCNanoFood > NCZooNano.and.PCNanoFood <= PCZooNano) then

           ZooNanoExcNH4 =  ZooNanoFoodN - ZooNanoFoodP / PCZooNano * NCZooNano           
           ZooNanoExcPO4 =  0.D0
           ZooNanoExcCO2 =  ZooNanoFoodC - ZooNanoFoodP / PCZooNano

         else 

           ZooNanoExcNH4 =  0.D0           
           ZooNanoExcPO4 =  ZooNanoFoodP - ZooNanoFoodN / NCZooNano * PCZooNano
           ZooNanoExcCO2 =  ZooNanoFoodC - ZooNanoFoodN / NCZooNano

         endif


! Activity Respiration

        ZooMicroGrazNetGrowthEffC = kcNetGrowthEff * ( ZooMicroGrazC - MicroEgeCSMOP - MicroMessyFeedCMOD ) 
!       ZooMicroGrazNetGrowthEffC = kcNetGrowthEff *   ZooMicroGrazC 


! Stoichiometric fluxes 

         NPZooMicro = NCZooMicro / PCZooMicro

         ZooMicroFoodC = ZooMicroGrazNetGrowthEffC 
!        ZooMicroFoodC = ZooMicroGrazNetGrowthEffC - MicroEgeCSMOP - MicroMessyFeedCMOD
         ZooMicroFoodN = ZooMicroGrazN             - MicroEgeNSMOP - MicroMessyFeedNMOD
         ZooMicroFoodP = ZooMicroGrazP             - MicroEgePSMOP - MicroMessyFeedPMOD
 


         if ( ZooMicroFoodC > 0.D0) then
           NCMicroFood =  ZooMicroFoodN  /  ZooMicroFoodC 
           PCMicroFood =  ZooMicroFoodP  /  ZooMicroFoodC 
         else 
           NCMicroFood = 0.
           PCMicroFood = 0. 
         endif


         if ( ZooMicroFoodP > 0.D0) then
           NPMicroFood =  ZooMicroFoodN  /  ZooMicroFoodP
         else 
           NPMicroFood = 0.
         endif


         if (NCMicroFood <= NCZooMicro.and.PCMicroFood <= PCZooMicro) then

           if ( NPMicroFood > NPZooMicro ) then
 
              ZooMicroExcNH4 =  ZooMicroFoodN - ZooMicroFoodP * NPZooMicro               
              ZooMicroExcCO2 =  ZooMicroFoodC - ZooMicroFoodP / PCZooMicro
              ZooMicroExcPO4 =  0.D0

           elseif ( NPMicroFood < NPZooMicro ) then

              ZooMicroExcPO4 =  ZooMicroFoodP - ZooMicroFoodN / NPZooMicro               
              ZooMicroExcCO2 =  ZooMicroFoodC - ZooMicroFoodN / NCZooMicro
              ZooMicroExcNH4 =  0.D0

           endif

         elseif (NCMicroFood > NCZooMicro.and.PCMicroFood > PCZooMicro) then 

           ZooMicroExcCO2 =  0.D0
           ZooMicroExcNH4 =  ZooMicroFoodN - ZooMicroFoodC * NCZooMicro
           ZooMicroExcPO4 =  ZooMicroFoodP - ZooMicroFoodC * PCZooMicro
           
         elseif (NCMicroFood > NCZooMicro.and.PCMicroFood <= PCZooMicro) then

           ZooMicroExcNH4 =  ZooMicroFoodN - ZooMicroFoodP / PCZooMicro * NCZooMicro           
           ZooMicroExcPO4 =  0.D0
           ZooMicroExcCO2 =  ZooMicroFoodC - ZooMicroFoodP / PCZooMicro

         else 

           ZooMicroExcNH4 =  0.D0           
           ZooMicroExcPO4 =  ZooMicroFoodP - ZooMicroFoodN / NCZooMicro* PCZooMicro
           ZooMicroExcCO2 =  ZooMicroFoodC - ZooMicroFoodN / NCZooMicro

         endif

! Activity respiration

         ZooMesoGrazNetGrowthEffC = kcNetGrowthEff * ( ZooMesoGrazC - MesoEgeCSMOP - MesoMessyFeedCMOD ) 
!        ZooMesoGrazNetGrowthEffC = kcNetGrowthEff *   ZooMesoGrazC


! Stoichiometric fluxes
         NPZooMeso = NCZooMeso / PCZooMeso


         ZooMesoFoodC = ZooMesoGrazNetGrowthEffC 
!        ZooMesoFoodC = ZooMesoGrazNetGrowthEffC - MesoEgeCSMOP - MesoMessyFeedCMOD
         ZooMesoFoodN = ZooMesoGrazN             - MesoEgeNSMOP - MesoMessyFeedNMOD
         ZooMesoFoodP = ZooMesoGrazP             - MesoEgePSMOP - MesoMessyFeedPMOD
 
        

         if ( ZooMesoFoodC > 0.D0) then
           NCMesoFood =  ZooMesoFoodN  /  ZooMesoFoodC 
           PCMesoFood =  ZooMesoFoodP  /  ZooMesoFoodC 
         else 
           NCMesoFood = 0.
           PCMesoFood = 0.
         endif
 

         if ( ZooMesoFoodP > 0.D0) then
           NPMesoFood =  ZooMesoFoodN  /  ZooMesoFoodP
         else 
           NPMesoFood = 0.
         endif


         if (NCMesoFood <= NCZooMeso.and.PCMesoFood <= PCZooMeso) then

           if ( NPMesoFood > NPZooMeso) then
 
             ZooMesoExcNH4 =  ZooMesoFoodN - ZooMesoFoodP * NPZooMeso               
             ZooMesoExcCO2 =  ZooMesoFoodC - ZooMesoFoodP / PCZooMeso
             ZooMesoExcPO4 =  0.D0

           elseif ( NPMesoFood < NPZooMeso ) then

             ZooMesoExcPO4 =  ZooMesoFoodP - ZooMesoFoodN / NPZooMeso               
             ZooMesoExcCO2 =  ZooMesoFoodC - ZooMesoFoodN / NCZooMeso
             ZooMesoExcNH4 =  0.D0

           endif

         elseif (NCMesoFood > NCZooMeso.and.PCMesoFood > PCZooMeso) then 

           ZooMesoExcCO2 =  0.D0
           ZooMesoExcNH4 =  ZooMesoFoodN - ZooMesoFoodC * NCZooMeso
           ZooMesoExcPO4 =  ZooMesoFoodP - ZooMesoFoodC * PCZooMeso

         elseif (NCMesoFood > NCZooMeso.and.PCMesoFood <= PCZooMeso) then

           ZooMesoExcNH4 =  ZooMesoFoodN - ZooMesoFoodP / PCZooMeso * NCZooMeso           
           ZooMesoExcPO4 =  0.D0
           ZooMesoExcCO2 =  ZooMesoFoodC - ZooMesoFoodP / PCZooMeso

         else 

           ZooMesoExcNH4 =  0.D0           
           ZooMesoExcPO4 =  ZooMesoFoodP - ZooMesoFoodN / NCZooMeso* PCZooMeso
           ZooMesoExcCO2 =  ZooMesoFoodC - ZooMesoFoodN / NCZooMeso

         endif


! Total excretion

         ZooExcNH4 =  ZooNanoExcNH4 + ZooMicroExcNH4 + ZooMesoExcNH4           
         ZooExcPO4 =  ZooNanoExcPO4 + ZooMicroExcPO4 + ZooMesoExcPO4 
         ZooExcCO2 =  ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2 

         if(ZooExcNH4.LT.-1e-20.OR.ZooExcPO4.LT.-1e-20.OR.ZooExcCO2.LT.-1e-20) then
           write(6,*)'ZooExcNH4',ZooExcNH4,ZooExcPO4,ZooExcCO2
           write(6,*)'ZooExcCO2',ZooNanoExcCO2,ZooMicroExcCO2,ZooMesoExcCO2
           write(6,*)'ZooFood',ZooMesoFoodC,ZooMesoGrazNetGrowthEffC,MesoEgeCSMOP,MesoMessyFeedCMOD 
           write(6,*)'ZooGraz',ZooMesoGrazZooMicroC,ZooMesoGrazDiaC,ZooMesoGrazSMOPC,ZooMesoGrazLMOPC,ZooMesoGrazC
           write(6,*)'grazMeso',grazMeso,ZooMesoC(I),ZooMicroC(I),prefMesoZooMicro,fZooMeso
           write(6,*)'autres',DiaC(I),prefMesoDia,S_MOPC(I),prefMesoSMOP,L_MOPC(I),prefMesoLMOP
          stop 'Excretion Zoo negative'
         endif 


!~~~~~~~~~~~~!
! MORTALITY  !
!~~~~~~~~~~~~!

! Natural mortality for Nano Zooplankton

         ZooNanoMortC       =  NanoMortRate *     ZooNanoC(I) * Tempfaczoo !** 2.   
!         ZooNanoMortC       =  NanoMortRate *     ZooNanoC(I)**2. * Tempfaczoo !** 2.
         ZooNanoMortSMOPC   =  ZooNanoMortC *    ZooPropMOPC   
         ZooNanoMortLMOPC   =  ZooNanoMortC * (1-ZooPropMOPC)  

         ZooNanoMortN       =  ZooNanoMortC * NCZooNano
         ZooNanoMortSMOPN   =  ZooNanoMortN *    ZooPropMOPN   
         ZooNanoMortLMOPN   =  ZooNanoMortN * (1-ZooPropMOPN)  

         ZooNanoMortP       =  ZooNanoMortC * PCZooNano
         ZooNanoMortSMOPP   =  ZooNanoMortP *    ZooPropMOPP  
         ZooNanoMortLMOPP   =  ZooNanoMortP * (1-ZooPropMOPP) 


! Natural mortality for Micro Zooplankton

         ZooMicroMortC       =  MicroMortRate *     ZooMicroC(I) * Tempfaczoo !** 2. 
!         ZooMicroMortC       =  MicroMortRate *     ZooMicroC(I)**2. * Tempfaczoo !** 2.
         ZooMicroMortSMOPC   =  ZooMicroMortC *    ZooPropMOPC   
         ZooMicroMortLMOPC   =  ZooMicroMortC * (1-ZooPropMOPC)  

         ZooMicroMortN       =  ZooMicroMortC * NCZooMicro
         ZooMicroMortSMOPN   =  ZooMicroMortN *    ZooPropMOPN   
         ZooMicroMortLMOPN   =  ZooMicroMortN * (1-ZooPropMOPN)  

         ZooMicroMortP       =  ZooMicroMortC * PCZooMicro
         ZooMicroMortSMOPP   =  ZooMicroMortP *    ZooPropMOPP  
         ZooMicroMortLMOPP   =  ZooMicroMortP * (1-ZooPropMOPP) 


! Predation on Meso Zooplankton 
 
!        ZooMesoPredC       =  MesoPredRate *     ZooMesoC(I)** 2. * Tempfaczoo / (Kclos + ZooMesoC(I)) 
!        ZooMesoPredC       =  MesoPredRate *     ZooMesoC(I)** 2. * Tempfaczoo !/ (Kclos + ZooMesoC(I))      !22/09/08
         ZooMesoPredC       =  MesoPredRate *     ZooMesoC(I)** 2. * Tempfaczoo / (Kclos + ZooMesoC(I)) !2021-03-23
!        ZooMesoPredC       =  2.3e-7       *     ZooMesoC(I)      * Tempfaczoo !/ (Kclos + ZooMesoC(I))
         ZooMesoPredSMOPC   =  ZooMesoPredC *    ZooPropMesoMOPC   
         ZooMesoPredLMOPC   =  ZooMesoPredC * (1-ZooPropMesoMOPC)  

         ZooMesoPredN       =  ZooMesoPredC * NCZooMeso
         ZooMesoPredSMOPN   =  ZooMesoPredN *    ZooPropMesoMOPN   
         ZooMesoPredLMOPN   =  ZooMesoPredN * (1-ZooPropMesoMOPN)  

         ZooMesoPredP       =  ZooMesoPredC * PCZooMeso
         ZooMesoPredSMOPP   =  ZooMesoPredP *    ZooPropMesoMOPP  
         ZooMesoPredLMOPP   =  ZooMesoPredP * (1-ZooPropMesoMOPP) 



!~~~~~~~~~~~!
! BACTERIA  !
!~~~~~~~~~~~!

         Tempfacbact = Q10_B**((Temperature(I)-TREF_B)/10.)


! Uptake of Dissolved Organic Matter

         LimitationBactMODC = MODC(I)     / ( MODC(I)     + KMODCBact)

         UptBactMODC   =  UptmaxBact * BactC(I)  * LimitationBactMODC !* Tempfaczoo
         UptBactMODN   =  NCMOD * UptBactMODC
         UptBactMODP   =  PCMOD * UptBactMODC


! Uptake of Dissolved Inorganic Matter

         LimitationBactAmmo = Ammonium(I)     / ( Ammonium(I)      + KAmmoBact)
         LimitationBactP    = Phosphate(I)    / ( Phosphate(I)     + KPBact)

         UptBactAmmoS   =  UptmaxBact * BactC(I) * NCBact * LimitationBactAmmo !* Tempfacbact
         UptBactPS      =  UptmaxBact * BactC(I) * PCBact * LimitationBactP !* Tempfacbact



! Growth 

        GrowthBact =      EffGrowthBact   * UptBactMODC


! Excretion

         NPBact= NCBact / PCBact

         BactFoodC = EffGrowthBact   * UptBactMODC 
         BactFoodN = UptBactMODN
         BactFoodP = UptBactMODP

         TestBactN = BactFoodN - BactFoodC * NCBact  
         TestBactP = BactFoodP - BactFoodC * PCBact

         if ( BactFoodC > 0.D0) then
           NCBactFood =  BactFoodN  /  BactFoodC 
           PCBactFood =  BactFoodP  /  BactFoodC 
         else 
           NCBactFood = 0.
           PCBactFood = 0.
         endif


         if ( BactFoodP > 0.D0) then
           NPBactFood =  BactFoodN  /  BactFoodP
         else 
           NPBactFood = 0.
         endif


! Limitation en P et N: Uptake PO4 et Ammo
         if (NCBactFood <= NCBact.and.PCBactFood <= PCBact) then


           if ( NPBactFood > NPBact ) then

             if ( abs(TestBactP) .gt. UptBactPS) then
               UptBactP    =    UptBactPS
             else
               UptBactP    =  - TestBactP
             endif

             if((UptBactP  + UptBactMODP ) * NPBact .gt. UptBactMODN ) then   !24-08-10

                UptBactAmmo = (UptBactP  + UptBactMODP ) * NPBact - UptBactMODN

               if (UptBactAmmo.gt.UptBactAmmoS) then
                 UptBactAmmo  =   UptBactAmmoS
                 UptBactP     =  (UptBactAmmo + UptBactMODN ) / NPBact  - UptBactMODP

               endif

               GrowthBact  = (UptBactP  + UptBactMODP ) / PCBact
               BactExcNH4  = 0.
               BactExcPO4  = 0.

             else

               UptBactAmmo = 0.
               GrowthBact  = (UptBactP  + UptBactMODP ) / PCBact
               BactExcNH4  = UptBactMODN - (UptBactP  + UptBactMODP ) * NPBact
               BactExcPO4  = 0.

             endif

           elseif ( NPBactFood < NPBact ) then
             if ( abs(TestBactN) .gt. UptBactAmmoS) then
                UptBactAmmo    =    UptBactAmmoS
             else
                UptBactAmmo    =  - TestBactN
             endif

             if ((UptBactAmmo  + UptBactMODN ) /  NPBact .gt. UptBactMODP ) then  !24-08-10

               UptBactP = (UptBactAmmo  + UptBactMODN ) /  NPBact - UptBactMODP

               if (UptBactP.gt.UptBactPS) then
                 UptBactP    =   UptBactPS
                 UptBactAmmo =  (UptBactP + UptBactMODP ) * NPBact  - UptBactMODN
               endif

               GrowthBact  = (UptBactAmmo  + UptBactMODN ) / NCBact
               BactExcNH4  = 0.
               BactExcPO4  = 0.

             else

               UptBactP = 0.
               GrowthBact  = (UptBactAmmo  + UptBactMODN ) / NCBact
               BactExcNH4  = 0.
               BactExcPO4  = UptBactMODP - (UptBactAmmo  + UptBactMODN ) /  NPBact

             endif

             LimitationBact(I)=1.

           endif

! Limitation en C : on ne consomme pas d'Ammo et de PO4, on n'en excrete
         elseif (NCBactFood > NCBact.and.PCBactFood > PCBact) then 

           UptBactP    =  0.D0
           UptBactAmmo =  0.D0
           BactExcNH4  =  TestBactN           
           BactExcPO4  =  TestBactP
           
           LimitationBact(I)=2.

!write(6,*)'doit passer par la'

         elseif (NCBactFood > NCBact.and.PCBactFood <= PCBact) then

           UptBactAmmo =  0.D0

           if ( abs(TestBactP) .gt. UptBactPS) then
             UptBactP =   UptBactPS
           else
             UptBactP = - TestBactP
           endif

           GrowthBact  = (UptBactP  + UptBactMODP ) / PCBact
           BactExcPO4  = 0.
           BactExcNH4  = UptBactMODN - GrowthBact * NCBact

           LimitationBact(I)=3.

         else 

           UptBactP =  0.D0

           if ( abs(TestBactN) .gt. UptBactAmmoS) then
             UptBactAmmo    =    UptBactAmmoS
           else
             UptBactAmmo    =  - TestBactN
           endif

           GrowthBact  = (UptBactAmmo  + UptBactMODN ) / NCBact
           BactExcPO4  = UptBactMODP - GrowthBact * PCBact 
           BactExcNH4  = 0.

           LimitationBact(I)=4.

         endif

! Respiration 

         RespBact = ( 1. - EffGrowthBact ) / EffGrowthBact * GrowthBact 
 

         if(BactExcNH4.LT.0.OR.BactExcPO4.LT.0.) then
           write(6,*)'BactExcNH4',BactExcNH4,BactExcPO4
           stop 'Excretion negative'
         endif 

! Mortality (bacteria become DOC, DON and DOP)

         MortBactC = MortBactRate * BactC(I) * Tempfacbact

         MortBactN = MortBactC * NCBact
         MortBactP = MortBactC * PCBact

!STOP'Verif Bact'

!~~~~~~~~~~~~~~~~~~~~~~~~~~!
! MINERALISATION OF MOP    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~!

!  Reference temperature = 19 dg C 
         TempfacRem   = Q10_Rem**((Temperature(I)-TREF_REM)/10.)

! Fonction de Karline
!        TempfacRem    = Exp(log(Q10_Rem)*(Temperature(I)-20.)/10.)

         RemSMOPC     = RemRateSMOPC   * TempfacRem * S_MOPC(I)
         RemSMOPN     = RemRateSMOPN   * TempfacRem * S_MOPN(I)
         RemSMOPP     = RemRateSMOPP   * TempfacRem * S_MOPP(I)
         RemSMOPSi    = RemRateSMOPSi  * TempfacRem * S_MOPSi(I)
         RemSMOPChl   = RemRateSMOPchl * TempfacRem * S_MOPChl(I)

         RemLMOPC     = RemRateLMOPC   * TempfacRem * L_MOPC(I)
         RemLMOPN     = RemRateLMOPN   * TempfacRem * L_MOPN(I)
         RemLMOPP     = RemRateLMOPP   * TempfacRem * L_MOPP(I)
         RemLMOPSi    = RemRateLMOPSi  * TempfacRem * L_MOPSi(I)



!~~~~~~~~~~~~~~~~~!
! NITRIFICATION   !
!~~~~~~~~~~~~~~~~~!

! Fonction de Fred
         TempfacNit   = Q10_Nit**((Temperature(I)-TREF_NIT)/10.)
         EPAR_SURF     = max(1d-8,PAR_Z(1))
         Nitrification = NitriRate   * Ammonium(I) * TempfacNit * (1. - PARZ/EPAR_SURF)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! GROSS PRIMARY PRODUCTION    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! Fonctions from Droop 
         NCSyneeff = max(NCSyne,minNCSyne)
         PCSyneeff = max(PCSyne,minPCSyne)

!        if(NCSyneeff/maxNCSyne.GE.1.and.PCSyneeff/maxPCSyne.GE.1) then
!gmlSyne = 1. 
!        else
         if (NCSyneeff/maxNCSyne < PCSyneeff/maxPCSyne) then
           gmlSyne = 1. - minNCSyne / NCSyneeff
           gmlSyne = min(gmlSyne,1.D0)
         else
           gmlSyne = 1. - minPCSyne / PCSyneeff
           gmlSyne = min(gmlSyne,1.D0)
         endif
!endif


!        NCNanoeff = max(NCNano,minNCNano)
!        PCNanoeff = max(PCNano,minPCNano)


!!        if(NCNanoeff/maxNCNano.GE.1.and.PCNanoeff/maxPCNano.GE.1) then
!!gmlNano = 1. 
!!        else
!        if (NCNanoeff/maxNCNano < PCNanoeff/maxPCNano) then
!        gmlNano = 1. - minNCNano / NCNanoeff
!   gmlNano = min(gmlNano,1.)
!else
!gmlNano = 1. - minPCNano / PCNanoeff
!   gmlNano = min(gmlNano,1.)
!endif
!!    endif


!        NCDiaeff = max(NCDia,minNCDia)
!        PCDiaeff = max(PCDia,minPCDia)
!SiCDiaeff= max(SiCDia,minSiCDia)

!        func = min(NCDiaeff/maxNCDia,PCDiaeff/maxPCDia,SiCDiaeff/maxSiCDia)

!!        if(NCDiaeff/maxNCDia.GE.1.and.PCDiaeff/maxPCDia.GE.1.and.SiCDiaeff/maxSiCDia.GE.1) then
!!gmlDia = 1.
!        if (func == NCDiaeff/maxNCDia) then
!        gmlDia = 1. - minNCDia  / NCDiaeff
!   gmlDia = min(gmlDia,1.)
!elseif (func == PCDiaeff/maxPCDia) then
!gmlDia = 1. - minPCDia  / PCDiaeff
!gmlDia = min(gmlDia,1.)
!else
!        gmlDia = 1. - minSiCDia / SiCDiaeff
!gmlDia = min(gmlDia,1.)
!endif


! Fonctions from Caperon and Meyer

!        NCSyneeff = max(NCSyne,minNCSyne)
!        PCSyneeff = max(PCSyne,minPCSyne)


!        if(NCSyneeff/maxNCSyne.GE.1.or.PCSyneeff/maxPCSyne.GE.1) then
!gmlSyne = 1. 
!        else
!if (NCSyneeff/maxNCSyne < PCSyneeff/maxPCSyne) then
!        gmlSyne = (NCSyneeff - minNCSyne) / (NCSyneeff - minNCSyne + 0.0114)
!else
!gmlSyne = (PCSyneeff - minPCSyne) / (PCSyneeff - minPCSyne + 0.00034)
!    endif
!endif

         NCNanoeff = max(NCNano,minNCNano)
         PCNanoeff = max(PCNano,minPCNano)


!        if(NCNanoeff/maxNCNano.GE.1.and.PCNanoeff/maxPCNano.GE.1) then
!gmlNano = 1. 
!        else
         if (NCNanoeff/maxNCNano < PCNanoeff/maxPCNano) then
           gmlNano = (NCNanoeff - minNCNano) / (NCNanoeff - minNCNano + 0.0072)
           gmlNano = min(gmlNano,1.D0)
         else
           gmlNano = (PCNanoeff - minPCNano) / (PCNanoeff - minPCNano + 0.0002)
           gmlNano = min(gmlNano,1.D0)
         endif
!endif


         NCDiaeff = max(NCDia,minNCDia)
         PCDiaeff = max(PCDia,minPCDia)
         SiCDiaeff= max(SiCDia,minSiCDia)

         func = min(NCDiaeff/maxNCDia,PCDiaeff/maxPCDia,SiCDiaeff/maxSiCDia)
 
!        if(NCDiaeff/maxNCDia.GE.1.and.PCDiaeff/maxPCDia.GE.1.and.SiCDiaeff/maxSiCDia.GE.1) then
!gmlDia = 1. 
!        else
         if (func == NCDiaeff/maxNCDia) then
           gmlDia = (NCDiaeff - minNCDia) / (NCDiaeff - minNCDia + 0.002)
           gmlDia = min(gmlDia,1.D0)
         elseif(func == PCDiaeff/maxPCDia) then
           gmlDia = (PCDiaeff - minPCDia) / (PCDiaeff - minPCDia + 0.0005)
           gmlDia = min(gmlDia,1.D0)
         else
           gmlDia = (SiCDiaeff - minSiCDia) / (SiCDiaeff - minSiCDia + 0.004) &
                            * NCDiaeff**10. / (NCDiaeff**10. + 0.1**10.)
           gmlDia = min(gmlDia,1.D0)
         endif
!endif
                           

         TempfacppbSyne    = Exp(log(Q10_P)*(Temperature(I)-TREF_P)/10.)
         TempfacppbNano    = Exp(log(Q10_P)*(Temperature(I)-TREF_P)/10.)
         TempfacppbDia     = Exp(log(Q10_P)*(Temperature(I)-TREF_P)/10.)



!        endif

         ! mmolC/mgChl/hr
         ! Test Alex 14/04/2016 de suppression de Gamma_M : ici et dans la fonction f_pproduction
         fppbSyne   = f_pproduction(AbsChlSyne,PhymaxSyne,Sig_Ps2Syne,PARZ,kd,krep,Tau_renewSyne,Gamma_M,TempfacppbSyne) ! Gamma_M,TempfacppbSyne
         fppbNano   = f_pproduction(AbsChlNano,PhymaxNano,Sig_Ps2Nano,PARZ,kd,krep,Tau_renewNano,Gamma_M,TempfacppbNano) !
         fppbDia    = f_pproduction(AbsChlDia ,PhymaxDia ,Sig_Ps2Dia ,PARZ,kd,krep,Tau_renewDia,Gamma_M,TempfacppbDia) !     
      


         ! Internal quota limited photosynthesis
         PPBSyne   = SyneChl(I) * fppbSyne      ! mmolC/m3/hr
         PPBNano   = NanoChl(I) * fppbNano
         PPBDia    = DiaChl(I)  * fppbDia   

         ! Internal quota replete growth rate 
         MuphyNRSyne = fppbSyne * ChlCSyne   ! molC/molC/hr    
         MuphyNRNano = fppbNano * ChlCNano
         MuphyNRDia  = fppbDia  * ChlCDia  

         ! Internal quota limited growth rate
         MuphySyne   = MuphyNRSyne * gmlSyne  ! /hr        
         MuphyNano   = MuphyNRNano * gmlNano
         MuphyDia    = MuphyNRDia  * gmlDia



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! EXUDATION OF DISSOLVED ORGANIC CARBON !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

         ExuSyneC =  PPBSyne * (1. - gmlSyne) ! mmolC/m3/hr
         ExuNanoC =  PPBNano * (1. - gmlNano)
         ExuDiaC  =  PPBDia  * (1. - gmlDia)
         ExuC     = ExuSyneC + ExuNanoC + ExuDiaC 


!~~~~~~~~~~~~~!
! RESPIRATION !
!~~~~~~~~~~~~~! 

         RespSyneC = RespCostSyne * PPBSyne * gmlSyne
         RespNanoC = RespCostNano * PPBNano * gmlNano
         RespDiaC  = RespCostDia  * PPBDia  * gmlDia


!~~~~~~~~~~~~~~~~~~~~~!
! NUTRIMENT UPTAKE    !
!~~~~~~~~~~~~~~~~~~~~~!

         LimitationSyneAmmo    =                    Ammonium(I)  / (Ammonium(I)  + KAmmoSyne)
         LimitationSyneNitrate =                    Nitrate(I)   / (Nitrate(I)   + KNitrateSyne)
         LimitationSyneP       =                    Phosphate(I) / (Phosphate(I) + KPSyne)

         LimitationNanoAmmo    =                    Ammonium(I)  / (Ammonium(I)  + KAmmoNano)
         LimitationNanoNitrate =                    Nitrate(I)   / (Nitrate(I)   + KNitrateNano)
         LimitationNanoP       =                    Phosphate(I) / (Phosphate(I) + KPNano)

         LimitationDiaAmmo     =                    Ammonium(I)  / (Ammonium(I)  + KAmmoDia)
         LimitationDiaNitrate  =                    Nitrate(I)   / (Nitrate(I)   + KNitrateDia) 
         LimitationDiaP        =                    Phosphate(I) / (Phosphate(I) + KPDia)
         LimitationDiaSi       =                    Silice(I)    / (Silice(I)    + KSiDia)
         
         InhibitionAmmonium = (1. - InhibAmmo      * Ammonium(I) / (Ammonium(I)  + KInhibAmmo))

         NCSyneeff = min(maxNCSyne,NCSyne)
         NCSyneeff = max(minNCSyne,NCSyneeff)

         PCSyneeff = min(maxPCSyne,PCSyne)
         PCSyneeff = max(minPCSyne,PCSyneeff)

         NCNanoeff = min(maxNCNano,NCNano)
         NCNanoeff = max(minNCNano,NCNanoeff)

         PCNanoeff = min(maxPCNano,PCNano)
         PCNanoeff = max(minPCNano,PCNanoeff)


         NCDiaeff = min(maxNCDia,NCDia)
         NCDiaeff = max(minNCDia,NCDiaeff)

         PCDiaeff = min(maxPCDia,PCDia)
         PCDiaeff = max(minPCDia,PCDiaeff)

         SiCDiaeff = min(maxSiCDia,SiCDia)
         SiCDiaeff = max(minSiCDia,SiCDiaeff)
         NCLimSyne =  ( (maxNCSyne - NCSyneeff)  / (maxNCSyne - minNCSyne) ) ** NGeider !2. 
         PCLimSyne =  ( (maxPCSyne - PCSyneeff)  / (maxPCSyne - minPCSyne) ) ** NGeider !2.
         
         NCLimNano =  ( (maxNCNano - NCNanoeff)  / (maxNCNano - minNCNano) ) ** NGeider !2. 
         PCLimNano =  ( (maxPCNano - PCNanoeff)  / (maxPCNano - minPCNano) ) ** NGeider !2.
         
         NCLimDia  =  ( (maxNCDia  -  NCDiaeff)  / (maxNCDia  -  minNCDia) ) ** NGeider !2. 
         PCLimDia  =  ( (maxPCDia  -  PCDiaeff)  / (maxPCDia  -  minPCDia) ) ** NGeider !2.
         SiCLimDia =  ( (maxSiCDia -  SiCDiaeff) / (maxSiCDia -  minSiCDia)) ** NGeider !2.


! uptake of nutrients depends on light !
         UptmaxSyneN  = muphyNRSyne   * maxNCSyne !* (1- NCSyne / maxNCSyne)
         UptmaxSyneP  = muphyNRSyne   * maxPCSyne !* (1- PCSyne / maxPCSyne)
 
         UptmaxNanoN  = muphyNRNano   * maxNCNano !* (1- NCNano / maxNCNano)
         UptmaxNanoP  = muphyNRNano   * maxPCNano !* (1- PCNano / maxPCNano)

         UptmaxDiaN   = muphyNRDia   * maxNCDia   !* (1- NCDia / maxNCDia)
         UptmaxDiaP   = muphyNRDia   * maxPCDia   !* (1- PCDia / maxPCDia)
         UptmaxDiaSi  = muphyNRDia   * maxSiCDia  !* (1- SiCDia / maxSiCDia)


         UptmaxSyneNit  = UptmaxSyneN !* 1./3.
         UptmaxSyneAmmo = UptmaxSyneN !* 2./3.
         UptmaxNanoNit  = UptmaxNanoN !* 1./3.
         UptmaxNanoAmmo = UptmaxNanoN !* 2./3.
         UptmaxDiaNit   = UptmaxDiaN  !* 1./3.
         UptmaxDiaAmmo  = UptmaxDiaN  !* 2./3.
         
                
         UptSyneNit  = UptmaxSyneNit  *  SyneC(I)  * LimitationSyneNitrate * InhibitionAmmonium  !mmolN/m3/hr
         UptSyneAmmo = UptmaxSyneAmmo *  SyneC(I)  * LimitationSyneAmmo
         UptSyneP    = UptmaxSyneP    *  SyneC(I)  * LimitationSyneP
         UptSyneN    = UptSyneNit + UptSyneAmmo 
         
         UptNanoNit  = UptmaxNanoNit  *  NanoC(I)  * LimitationNanoNitrate  * InhibitionAmmonium  !mmolN/m3/hr
         UptNanoAmmo = UptmaxNanoAmmo *  NanoC(I)  * LimitationNanoAmmo
         UptNanoP    = UptmaxNanoP    *  NanoC(I)  * LimitationNanoP
         UptNanoN    = UptNanoNit + UptNanoAmmo 
         
         
         UptDiaNit   = UptmaxDiaNit  *  DiaC(I)   * LimitationDiaNitrate !* InhibitionAmmonium
         UptDiaAmmo  = UptmaxDiaAmmo *  DiaC(I)   * LimitationDiaAmmo
         UptDiaP     = UptmaxDiaP    *  DiaC(I)   * LimitationDiaP
         UptDiaSi    = UptmaxDiaSi   *  DiaC(I)   * LimitationDiaSi
         UptDiaN     = UptDiaNit  + UptDiaAmmo
         
          
         UptNit  = UptSyneNit  + UptNanoNit  + UptDiaNit
         UptAmmo = UptSyneAmmo + UptNanoAmmo + UptDiaAmmo
         UptP    = UptSyneP    + UptNanoP    + UptDiaP
         
         
         ExcAmmo = 0.
         ExcPO    = 0.
         ExcDiaSi = 0.



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! CHLOROPHYLLE PRODUCTION    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~!



! Case of nutrient limitation Geider et al 1998


! Division par gamma_M  il me semble??? ajouter GML ??? division par tfunc T ?????
         if(ChlNSyne < ChlNSyneMax ) then                                 ! 22/09/08
           if(ChlCSyne > 0.D0) then
             rhoChlSyne = ChlNSyneMax * MuphySyne  / & 
                         ( AbsChlSyne * PhyMaxSyne * PARZ * ChlCSyne * 3600.)
           else
             rhoChlSyne = 0.D0
           endif
           UptNVelSyne = ( UptSyneAmmo + UptSyneNit ) * ( 1. - ChlNSyne / ChlNSyneMax ) &
                                                  / ( 1.05 - ChlNSyne / ChlNSyneMax )

           ChlprodSyne = rhoChlSyne * UptNVelSyne

         else
           ChlprodSyne = 0.
         endif


! Division par gamma_M  il me semble??? ajouter GML ??? division par tfunc T ?????
         if(ChlNNano < ChlNNanoMax) then                                  ! 22/09/08
           if(ChlCNano > 0.D0) then
             rhoChlNano = ChlNNanoMax * MuphyNano  / & 
                         ( AbsChlNano * PhyMaxNano * PARZ * ChlCNano * 3600.)
           else
             rhoChlNano = 0.D0
           endif
           UptNVelNano = ( UptNanoAmmo + UptNanoNit ) * ( 1. - ChlNNano / ChlNNanoMax ) &
                                                  / ( 1.05 - ChlNNano / ChlNNanoMax )

           ChlprodNano = rhoChlNano * UptNVelNano
         else
            ChlprodNano = 0.
         endif 

! Division par gamma_M  il me semble??? ajouter GML ??? division par tfunc T ?????
         if(ChlNDia < ChlNDiaMax ) then                                   ! 22/09/08 
           if(ChlCDia > 0.D0) then
             rhoChlDia   = ChlNDiaMax * MuphyDia   / & 
                          ( AbsChlDia * PhyMaxDia * PARZ * ChlCDia * 3600.)
           else
             rhoChlDia = 0.D0
           endif

           UptNVelDia  = ( UptDiaAmmo + UptDiaNit ) * ( 1. - ChlNDia / ChlNDiaMax ) &
                                                    / ( 1.05 - ChlNDia / ChlNDiaMax )

           ChlprodDia = rhoChlDia * UptNVelDia
         else
           ChlprodDia = 0.
         endif



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! EXUDATION OF DISSOLVED ORGANIC MATTER   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         ExuSyneAmmo = UptSyneAmmo * (1. - NCLimSyne)
         ExuSyneP    = UptSyneP    * (1. - PCLimSyne)
         ExuSyneN    = ExuSyneNit + ExuSyneAmmo

         ExuNanoNit  = UptNanoNit  * (1. - NCLimNano)
         ExuNanoAmmo = UptNanoAmmo * (1. - NCLimNano)
         ExuNanoP    = UptNanoP    * (1. - PCLimNano)
         ExuNanoN    = ExuNanoNit + ExuNanoAmmo

         ExuDiaNit   = UptDiaNit   * (1. - NCLimDia )
         ExuDiaAmmo  = UptDiaAmmo  * (1. - NCLimDia ) 
         ExuDiaP     = UptDiaP     * (1. - PCLimDia )
         ExuDiaSi    = UptDiaSi    * (1. - SiCLimDia)
         ExuDiaN     = ExuDiaNit  + ExuDiaAmmo

         ExuN  =  ExuDiaN + ExuNanoN + ExuSyneN
         ExuP  =  ExuDiaP + ExuNanoP + ExuSyneP


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! RESPIRATION LIEE A L'UPTAKE !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         RespSyneNit  = UptSyneNit  * fRespNit
         RespSyneAmmo = UptSyneAmmo * fRespAmmo
         RespSyneP    = UptSyneP    * fRespP

         RespNanoNit  = UptNanoNit  * fRespNit
         RespNanoAmmo = UptNanoAmmo * fRespAmmo
         RespNanoP    = UptNanoP    * fRespP

         RespDiaNit   = UptDiaNit   * fRespNit    ! mmolC/m3/yr
         RespDiaAmmo  = UptDiaAmmo  * fRespAmmo
         RespDiaP     = UptDiaP     * fRespP
         RespDiaSi    = UptDiaSi    * fRespSi


!~~~~~~~~~~~~~~~~~~~!
! NATURAL MORTALITY !
!~~~~~~~~~~~~~~~~~~~!
         MortDiaC   = MortDia *  DiaC(I) * TempfacppbDia
         MortDiaN   = MortDia *  DiaN(I) * TempfacppbDia
         MortDiaP   = MortDia *  DiaP(I) * TempfacppbDia
         MortDiaChl = MortDia *  DiaChl(I) * TempfacppbDia
         MortDiaSi  = MortDia *  DiaSi(I) * TempfacppbDia

         MortNanoC   = MortNano  * NanoC(I) * TempfacppbNano
         MortNanoN   = MortNano  * NanoN(I) * TempfacppbNano
         MortNanoP   = MortNano  * NanoP(I) * TempfacppbNano
         MortNanoChl = MortNano  * NanoChl(I) * TempfacppbNano

         MortSyneC   = MortSyne  * SyneC(I) * TempfacppbSyne
         MortSyneN   = MortSyne  * SyneN(I) * TempfacppbSyne 
         MortSyneP   = MortSyne  * SyneP(I) * TempfacppbSyne
         MortSyneChl = MortSyne  * SyneChl(I) * TempfacppbSyne


         SMOPMortPhytoC   = PropDia * MortDiaC   + MortNanoC   + MortSyneC
         SMOPMortPhytoN   = PropDia * MortDiaN   + MortNanoN   + MortSyneN
         SMOPMortPhytoP   = PropDia * MortDiaP   + MortNanoP   + MortSyneP
         SMOPMortPhytoChl = PropDia * MortDiaChl + MortNanoChl + MortSyneChl
         SMOPMortPhytoSi  = PropDia * MortDiaSi  

         LMOPMortPhytoC   = (1. - PropDia) * MortDiaC   
         LMOPMortPhytoN   = (1. - PropDia) * MortDiaN   
         LMOPMortPhytoP   = (1. - PropDia) * MortDiaP   
         LMOPMortPhytoChl = (1. - PropDia) * MortDiaChl 
         LMOPMortPhytoSi  = (1. - PropDia) * MortDiaSi  


!~~~~~~~~~~~~~~~~~~~!
! OXYGEN            !
!~~~~~~~~~~~~~~~~~~~!

! fonction temperature factor for oxidation

         Tempfacoxid=Q10_P**((Temperature(I)-20.)/10.)

! limitation par des formules chimiques

         limitOxygen_oxicmin_satOxygen      =  oxygen(I) / (oxygen(I) + kOxicMinSatOX)   
         limitOxygen_oduoxid_satOxygen      =  oxygen(I) / (oxygen(I) + kODUoxidSatOX)
         limitOxygen_Ammoniumoxid_satOxygen =  oxygen(I) / (oxygen(I) + kNHsoxidSatOX)
         limitOxygen_anoxmin_inOxygen       =  oxygen(I) / (oxygen(I) + kAnoxMinInOX)
         limitOxygen_oduoxidNit_inOxygen    =  oxygen(I) / (oxygen(I) + kODUoxidNOsInOX)
         limitOxygen_oduoxidNit_satOxygen   =  oxygen(I) / (oxygen(I) + kODUoxidNOsSatNOs)
         
         limitNitrate_anoxmin_inOxygen      = nitrate(I) / (nitrate(I) + kAnoxMinInNOs)
         limitNitrate_oduoxidNit_satNit     = nitrate(I) / (nitrate(I) + kODUoxidNOsSatNOs)
         
         limitIron_solform_satIron    = iron    / (iron    + kSolFormSatIron)


! respiration bacterienne avec limitatin de l'O2


         respBACoxygen       = respBact * limitOxygen_oxicmin_satOxygen                                              !mmolC/(m3)/d
         anoxicrespBAC       = respBact * (1. - limitNitrate_anoxmin_inOxygen) * (1. - limitOxygen_anoxmin_inOxygen)     !mmolC/(m3)/d

! Oxidation

         ODUoxidOxygen       = Tempfacoxid * nuOXYGENODU  * limitOxygen_oduoxid_satOxygen                                           ! /d
         ODUoxidNitrate      = Tempfacoxid * nuNITODU * limitNitrate_oduoxidNit_satNit * (1. - limitOxygen_oduoxidNit_inOxygen)    ! /d
         AmmoniumoxidOxygen  = Tempfacoxid * nuOXYGENNHs  * limitOxygen_Ammoniumoxid_satOxygen                                            ! /d
         
         solidODU            = Rsolid      * limitIron_solform_satiron              * gammaPOC_ODU * anoxicrespBAC   ! mmolODU/(m3)/d



!!=====================!!
!!=====================!!
!! THE RATES OF CHANGE !!
!!=====================!!
!!=====================!!

! Zooplancton


         dZOONanoC(I)   =  ZooNanoGrazNetGrowthEffC   - ZooGrazZooNanoC  -ZooNanoMortC    -  ZooNanoExcCO2
         dZOOMicroC(I)  =  ZooMicroGrazNetGrowthEffC  - ZooGrazZooMicroC -ZooMicroMortC  -  ZooMicroExcCO2
         dZOOMesoC(I)   =  ZooMesoGrazNetGrowthEffC   - ZooMesopredC                      -  ZooMesoExcCO2


! Phytoplankton 

         dDIAC(I)   =  - ZooGrazDiaC    + PPBDia      - ExuDiaC  - MortDiaC &
                       - RespDiaC - RespDiaNit - RespDiaAmmo - RespDiaP - RespDiaSi  
         dDIAN(I)   =  - ZooGrazDiaN    + UptDiaN     - ExuDiaN  - MortDiaN
         dDIAP(I)   =  - ZooGrazDiaP    + UptDiaP     - ExuDiaP  - MortDiaP
         dDIASi(I)  =  - ZooGrazDiaSi   + UptDiaSi    - ExuDiaSi - MortDiaSi  
         dDIAChl(I) =  - ZooGrazDiaChl  + ChlprodDia             - MortDiaChl

         dNanoC(I)   =  - ZooGrazNanoC   + PPBNano     - ExuNanoC - MortNanoC &
                        - RespNanoC - RespNanoNit - RespNanoAmmo - RespNanoP  
         dNanoN(I)   =  - ZooGrazNanoN   + UptNanoN    - ExuNanoN - MortNanoN                                
         dNanoP(I)   =  - ZooGrazNanoP   + UptNanoP    - ExuNanoP - MortNanoP
         dNanoChl(I) =  - ZooGrazNanoChl + ChlprodNano            - MortNanoChl

         dSyneC(I)   = - ZooGrazSyneC   + PPBSyne     - ExuSyneC - MortSyneC & 
                       - RespSyneC - RespSyneNit - RespSyneAmmo - RespSyneP  
         dSyneN(I)   = - ZooGrazSyneN   + UptSyneN    - ExuSyneN - MortSyneN                                
         dSyneP(I)   = - ZooGrazSyneP   + UptSyneP    - ExuSyneP - MortSyneP
         dSyneChl(I) = - ZooGrazSyneChl + ChlprodSyne            - MortSyneChl

! Bacteria

         dBactC(I)   = - ZooGrazBactC + GrowthBact - MortBactC


! Particulate Organic Matter
 
         dS_MOPC(I)  = - ZooGrazSMOPC   + EgeCSMOP   + ZooMesoPredSMOPC + ZooMicroMortSMOPC &
                      + ZooNanoMortSMOPC - RemSMOPC + SMOPMortPhytoC
         dS_MOPN(I)  = - ZooGrazSMOPN   + EgeNSMOP   + ZooMesoPredSMOPN + ZooMicroMortSMOPN &
                       + ZooNanoMortSMOPN - RemSMOPN  + SMOPMortPhytoN
         dS_MOPP(I)  = - ZooGrazSMOPP   + EgePSMOP   + ZooMesoPredSMOPP + ZooMicroMortSMOPP &
                      + ZooNanoMortSMOPP - RemSMOPP  + SMOPMortPhytoP 
         dS_MOPSi(I) =   EgeSiSMOP       - RemSMOPSi + SMOPMortPhytoSi
         dS_MOPChl(I)= - ZooGrazSMOPChl + EgeChlSMOP     - RemSMOPChl + SMOPMortPhytoChl 

         dL_MOPC(I)  = - ZooGrazLMOPC + ZooMesoPredLMOPC  + ZooMicroMortLMOPC + ZooNanoMortLMOPC &
                       - RemLMOPC + LMOPMortPhytoC
         dL_MOPN(I)  = - ZooGrazLMOPN               + ZooMesoPredLMOPN  + ZooMicroMortLMOPN + ZooNanoMortLMOPN & 
                       - RemLMOPN + LMOPMortPhytoN
         dL_MOPP(I)  = - ZooGrazLMOPP               + ZooMesoPredLMOPP  + ZooMicroMortLMOPP + ZooNanoMortLMOPP & 
                        - RemLMOPP + LMOPMortPhytoP
         dL_MOPSi(I) =  EgeSiLMOP    - RemLMOPSi + LMOPMortPhytoSi

! Dissolved Oragnic Matter

         dMODC(I)    =  MessyFeedCMOD         + ExuC - GrowthBact - RespBact + MortBactC + RemSMOPC + RemLMOPC                          
         dMODN(I)    =  MessyFeedNMOD         + ExuN            - UptBactMODN + MortBactN + RemSMOPN + RemLMOPN                          
         dMODP(I)    =  MessyFeedPMOD         + ExuP            - UptBactMODP + MortBactP + RemSMOPP + RemLMOPP                          


! Dissolved Inorganic Matter

!        if(UptBactP.ne.0.or.ExcPO.ne.0) then
!        print*,'uptbactP,ExcPO',UptBactP,ExcPO,UptP
!        stop'dans Eco3M.F90'
!        endif    

         dNitrate(I)  =                          + Nitrification - UptNit 
         dAmmonium(I) =  ZooExcNH4  + BactExcNH4 - Nitrification - UptAmmo               + ExcAmmo - UptBactAmmo
         dPhosphate(I)=  ZooExcPO4  + BactExcPO4                 - UptP                  + ExcPO   - UptBactP
         dSilice(I)   =      ExuDiaSi                            - UptDiaSi  + RemSMOPSi + RemLMOPSi 

! Oxygen and ODU

    ! taux du ODU = "oxygen demand unit"


!       dODU(I)       = gammaPOC_ODU * anoxicrespBAC                                       &
!                     - ODUoxidNitrate * ODU(I)                                            & ! a verifier
!                     - ODUoxidOxygen  * ODU(I)                                            & ! a verifier
!                     - solidODU                                                                                        ! mmolODU/m3


    ! concentration en oxygen


!       dOxygen(I)    = (ChlprodNano+ChlprodSyne+ChlprodDia)                                   * gammaC_OX       & !modif caroline
        dOxygen(I)    = (PPBSyne+PPBNano+PPBDia)                                               * gammaC_OX       & ! modif caroline
                      - (RespDiaC      + RespNanoC      + RespSyneC)                           * gammaC_OX      &  
                      - (RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi                                      &  ! modif caroline
                        +RespNanoNit + RespNanoAmmo + RespNanoP                                                 &  ! modif caroline
                        +RespSyneNit + RespSyneAmmo + RespSyneP)                               * gammaC_OX      &  ! modif caroline
                      - (ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2                                         & 
                        +( ZooNanoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff))                      &
                        +( ZooMicroGrazNetGrowthEffC / kcNetGrowthEff *(1- kcNetGrowthEff))                      &
                        +( ZooMesoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)) )  * gammaC_OX       &
!                     - (respBACoxygen)                                                        * gammaC_OX       & ! modif caroline
                      - respBACt                                                               * gammaC_OX       & ! modif caroline
!                     - ((UptSyneNit+UptNanoNit+UptDiaNit) - (AmmoniumoxidOxygen*Ammonium(I)))    * gammaNHs_OX  & ! modif caroline
                      + ((UptSyneNit+UptNanoNit+UptDiaNit) - Nitrification) * gammaNHs_OX          ! modif caroline
!                      - Nitrification * gammaNHs_OX                              !modif caroline alex 10/04 (correction unit) 
!                    - (ODUoxidOxygen * ODU(I))                                               * gammaODU_OX          ! mmolO2/m3  modif caroline              

         dDIC(I) =    - ( PPBSyne + PPBNano + PPBDia) & 
                      + (RespBact &
                      + RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi      &
                      + RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP &
                      + RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP &
                      + ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2 &
                      + ZooNanoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)    &
                      + ZooMicroGrazNetGrowthEffC / kcNetGrowthEff *(1- kcNetGrowthEff)    &
                      + ZooMesoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff) )



! a voir : consommation de nutriments 
         PPBi(I)      = (PPBSyne + PPBNano + PPBDia)* 12. * 24. ! mgC/m3/d
         PPBpi(I)     = (PPBSyne                   )* 12. * 24. ! mgC/m3/d
         PPBni(I)     = (PPBNano                   )* 12. * 24. ! mgC/m3/d
         PPBmi(I)     = (PPBDia                    )* 12. * 24. ! mgC/m3/d
         GrazCi(I)    = (ZooNanoGrazNetGrowthEffC + ZooMicroGrazNetGrowthEffC + ZooMesoGrazNetGrowthEffC)* Thickness(I) * 12. * 24.
         RespCi(I)    = (RespBact &
                      + RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi      &
                      + RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP &
                      + RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP &
                      + ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2 &
                      + ZooNanoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)    &
                      + ZooMicroGrazNetGrowthEffC / kcNetGrowthEff *(1- kcNetGrowthEff)    &
                      + ZooMesoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff))*  12. * 24.
         RespCPi(I)   = (RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi      &
                       + RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP &
                       + RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP)* 12. * 24.
         RespCPpi(I)  = (RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP)* 12.* 24.
         RespCPni(I)  = (RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP)* 12.* 24.
         RespCPmi(I)  = (RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP+RespDiaSi)* 12.*24.
         RespCZi(I)   = (ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2 &
                        + ZooNanoGrazNetGrowthEffC  / kcNetGrowthEff *(1-kcNetGrowthEff)    &
                        + ZooMicroGrazNetGrowthEffC / kcNetGrowthEff *(1-kcNetGrowthEff)    &
                      + ZooMesoGrazNetGrowthEffC  / kcNetGrowthEff *(1-kcNetGrowthEff))*  12. * 24.
         RespCBi(I)   = (RespBact) * 12. * 24.
         ExuCi(I)     = (ExuC) * 12. * 24.
         MortCi(I)    = (MortSyneC + MortNanoC + MortDiaC) * 12. * 24. ! mgC/m2/d

         Nitrifi(I)   = Nitrification * 24.
         GrazHerbii(I)= ZooGrazHerbiC * Thickness(I) * 12. * 24.
         Gbaci(I)     = GrowthBact    * Thickness(I) * 12. * 24. 

!---------!
! BUDGET  !
!---------!



! Phyt
!       if(I.le.27) then
         TotalPPBSyne  =  TotalPPBSyne     + PPBSyne      * Thickness(I) * 12. * 24.! mgC/m2/d
         TotalPPBDia   =  TotalPPBDia      + PPBDia       * Thickness(I) * 12. * 24.! mgC/m2/d
         TotalPPBNano  =  TotalPPBNano     + PPBNano      * Thickness(I) * 12. * 24.! mgC/m2/d
         TotalPPB      = TotalPPB  + (PPBSyne + PPBNano + PPBDia ) * Thickness(I) * 12. * 24.! mgC/m2/d


         vPPBSyne(I) = PPBSyne* 24.
         vPPBNano(I) = PPBNano* 24.
         vPPBDia(I)  = PPBDia    * 24. 
         vPPBTotale(I)  = (PPBSyne + PPBNano + PPBDia) * 24.
         
         
         v2PPBSyne(I) = PPBSyne * Thickness(I) * 12. * 24.
         v2PPBNano(I) = PPBNano * Thickness(I) * 12. * 24.
         v2PPBDia(I)  = PPBDia  * Thickness(I) * 12. * 24.
         v2PPBTotale(I)  = (PPBSyne + PPBNano + PPBDia) * Thickness(I) * 12. * 24.
 

         TotalNetPPBSyne   =  TotalNetPPBSyne  &
                           + (PPBSyne - ExuSyneC - RespSyneC - RespSyneNit - RespSyneAmmo - RespSyneP) &
                           * Thickness(I) * 12. * 24.! mgC/m2/d
         TotalNetPPBDia    =  TotalNetPPBDia   & 
                           + (PPBDia  - ExuDiaC  - RespDiaC - RespDiaNit - RespDiaAmmo &
                                      - RespDiaP - RespDiaSi) &!
                           * Thickness(I) * 12. * 24.! mgC/m2/d
         TotalNetPPBNano   =  TotalNetPPBNano &
                           + (PPBNano - ExuNanoC - RespNanoC - RespNanoNit - RespNanoAmmo - RespNanoP) &
                           * Thickness(I) * 12. * 24.! mgC/m2/d
         TotalNetPPB       =  TotalNetPPBSyne + TotalNetPPBNano + TotalNetPPBDia! mgC/m2/d
!       endif


! New preoduction
!       if(I.LE.27) THEN
         TotalUptNit     = TotalUptNit     +UptNit     * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptSyneNit = TotalUptSyneNit +UptSyneNit * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptDiaNit  = TotalUptDiaNit  +UptDiaNit  * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptNanoNit = TotalUptNanoNit +UptNanoNit * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d

! Recycled production
         TotalUptSyneAmmo = TotalUptSyneAmmo +UptSyneAmmo * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptDiaAmmo  = TotalUptDiaAmmo  +UptDiaAmmo  * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptNanoAmmo = TotalUptNanoAmmo +UptNanoAmmo * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptBactAmmo = TotalUptBactAmmo +UptBactAmmo * Thickness(I) * 14. * 24. !*365./ 1000.! mgN/m2/d
         TotalUptAmmo = TotalUptSyneAmmo + TotalUptDiaAmmo + TotalUptNanoAmmo + TotalUptBactAmmo! mgN/m2/d

         TotalUptN = TotalUptN + (UptNit + UptAmmo)       * Thickness(I) * 14. * 24. !* 365. / 1000.! mgN/m2/d


! Nitrification
         TotalNitrif = TotalNitrif + Nitrification        * Thickness(I) * 14. * 24.! mgN/m2/d

         NewPPB = NewPPB + Max(0.D0, (UptNit - Nitrification)) * Thickness(I) * 14. * 24.

         if(TotalUptN.GT.0.) then
           fratio1 = TotalUptNit / TotalUptN
           fratio2 = NewPPB      / TotalUptN     
         endif
!       endif 
! Nitrification ALLCOLUMN
         if(I.le.NumPelagicBoxes-1) then
           TotalNitrifCOLUMN = TotalNitrifCOLUMN + Nitrification * Thickness(I) * 14. *24. ! mgN/m2/d ! Exudation 
         endif

         TotalExuC = TotalExuC + ExuC *12. *24. * Thickness(I) ! mgC/m2/d
         TotalExuN = TotalExuN + ExuN *14. *24. * Thickness(I) ! mgN/m2/d
         TotalExuP = TotalExuP + ExuP *31. *24. * Thickness(I) ! mgP/m2/d



         vUptNit(I)      = UptNit       * 24. ! mmolN/m3/d
         vUptSyneNit(I)  = UptSyneNit   * 24. ! mmolN/m3/d
         vUptDiaNit(I)   = UptDiaNit    * 24. ! mmolN/m3/d
         vUptNanoNit(I)  = UptNanoNit   * 24. ! mmolN/m3/d
         vUptSyneAmmo(I) = UptSyneAmmo  * 24. ! mmolN/m3/d
         vUptDiaAmmo(I)  = UptDiaAmmo   * 24. ! mmolN/m3/d
         vUptNanoAmmo(I) = UptNanoAmmo  * 24. ! mmolN/m3/d
         vUptBactAmmo(I) = UptBactAmmo  * 24. ! mmolN/m3/d
         vUptAmmo(I)     = UptAmmo      * 24. ! mmolN/m3/d
!         vNitrif(I)= Nitrification * 24.! mmolN/m3/d
         vZooExcNH4(I)   = ZooExcNH4    * 24. ! mmolN/m3/d
         vBactExcNH4(I)  = BactExcNH4   * 24. ! mmolN/m3/d
         vdAmmo(I) = ( ZooExcNH4  + BactExcNH4 - Nitrification - UptAmmo + ExcAmmo - UptBactAmmo)*24
         vBactExcPO4(I)  = BactExcPO4   * 24. ! mmolP/m3/d
         vUptBactP(I)    = UptBactP     * 24. ! mmolP/m3/d 
         vUptPhytoP(I)   = UptP         * 24. ! mmolP/m3/d 
!        vGrazHerbi(I)   = ZooGrazHerbiC  * 24. 


! Uptake phosphate
         TotalUptPhytoP= TotalUptPhytoP +UptP    * Thickness(I) * 31. * 24. !*365./ 1000.! mgP/m2/d
         TotalUptSyneP = TotalUptSyneP +UptSyneP * Thickness(I) * 31. * 24. !*365./ 1000.! mgP/m2/d
         TotalUptDiaP  = TotalUptDiaP  +UptDiaP  * Thickness(I) * 31. * 24. !*365./ 1000.! mgP/m2/d
         TotalUptNanoP = TotalUptNanoP +UptNanoP * Thickness(I) * 31. * 24. !*365./ 1000.! mgP/m2/d

! Grazing-Assimilation
         TotalGrazC  = TotalGrazC + (ZooMicroGrazNetGrowthEffC + ZooMicroGrazNetGrowthEffC + ZooMesoGrazNetGrowthEffC ) &
                      * Thickness(I) * 12. * 24. !* 365. / 1000. ! mgC/m2/d
! claude
         TotalGrazCPur = TotalGrazCPur + (ZooNanoGrazC + ZooMicroGrazC + ZooMesoGrazC) &
                      * Thickness(I) * 12. * 24. !* 365. / 1000.         ! mgC/m2/d
         TotalGrazCHerbi = TotalGrazCHerbi + ZooGrazHerbiC  &
                      * Thickness(I) * 12. * 24. !* 365. / 1000.         ! mgC/m2/d
! fin claude
         TotalGrazN  = TotalGrazN + (ZooGrazN - EgeNSMOP - MessyFeedNMOD) * Thickness(I) * 14. *24. ! mgN/m2/d
         TotalGrazP  = TotalGrazP + (ZooGrazP - EgePSMOP - MessyFeedPMOD) * Thickness(I) * 31. *24. ! mgP/m2/d


! Egestion
         TotalEgestionC   =TotalEgestionC + EgeCSMOP *12. *24. * Thickness(I) ! mgC/m2/d
         TotalEgestionN   =TotalEgestionN + EgeNSMOP *14. *24. * Thickness(I) ! mgN/m2/d
         TotalEgestionP   =TotalEgestionP + EgePSMOP *31. *24. * Thickness(I) ! mgP/m2/d


! Messyfeeding
         TotalMessyfeedC  =TotalMessyfeedC + MessyfeedCMOD *12. *24. * Thickness(I) ! mgC/m2/d
         TotalMessyfeedN  =TotalMessyfeedN + MessyfeedNMOD *14. *24. * Thickness(I) ! mgN/m2/d
         TotalMessyfeedP  =TotalMessyfeedP + MessyfeedPMOD *31. *24. * Thickness(I) ! mgP/m2/d


! Excretion Zoo
!        if(I.le.NumPelagicBoxes-1) then
!TotalExcZooAmmo  = TotalExcZooAmmo + ZooExcNH4 *14. *24. * Thickness(I) ! mgN/m2/d
!TotalExcZooP = TotalExcZooP    + ZooExcPO4 *31. *24. * Thickness(I) ! mgN/m2/d
!        endif

!       if(I.le.27) then
! Respiration
         TotalResp = TotalResp +  (                                               &
               RespBact                                                           &
             + RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi      &
             + RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP                 &
             + RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP                 &
             + ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2                     &
             + ZooNanoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)    &
             + ZooMicroGrazNetGrowthEffC / kcNetGrowthEff *(1- kcNetGrowthEff)    &
             + ZooMesoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)    &
                                  )*12. * 24. *Thickness(I)! mgC/m2/d
         TotalRespPhyto=TotalRespPhyto+ (                                         &
               RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi      &
             + RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP                 &
             + RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP                 &
                                  )*12. * 24. *Thickness(I)! mgC/m2/d

         TotalRespSyne=TotalRespSyne+ (                                           &
                RespSyneC + RespSyneNit + RespSyneAmmo + RespSyneP                &
                                   )*12. * 24. *Thickness(I)! mgC/m2/d

         TotalRespNano=TotalRespNano+ (                                           &
              + RespNanoC + RespNanoNit + RespNanoAmmo + RespNanoP                &
                                   )*12. * 24. *Thickness(I)! mgC/m2/d

         TotalRespDia=TotalRespDia+ (                                             &
                RespDiaC  + RespDiaNit  + RespDiaAmmo  + RespDiaP + RespDiaSi     &
                                   )*12. * 24. *Thickness(I)! mgC/m2/d

         TotalRespZoo = TotalRespZoo + (                                          &
               ZooNanoExcCO2 + ZooMicroExcCO2 + ZooMesoExcCO2                     &
             + ZooNanoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)    &
             + ZooMicroGrazNetGrowthEffC / kcNetGrowthEff *(1- kcNetGrowthEff)    &
             + ZooMesoGrazNetGrowthEffC  / kcNetGrowthEff *(1- kcNetGrowthEff)    &
                                  )*12. * 24. *Thickness(I)
!       endif

         if(I.le.NumPelagicBoxes-1) then
! Bacteries
           TotalRespBact    = TotalRespBact    + RespBact    *12. * 24. *Thickness(I)
           TotalGrowthBact  = TotalGrowthBact  + GrowthBact  *12. * 24. *Thickness(I)
           TotalUptBactMODC = TotalUptBactMODC + UptBactMODC *12. * 24. *Thickness(I)
           TotalUptBactMODN = TotalUptBactMODN + UptBactMODN *14. * 24. *Thickness(I)
           TotalUptBactMODP = TotalUptBactMODP + UptBactMODP *31. * 24. *Thickness(I)
           TotalUptBactP    = TotalUptBactP    + UptBactP    *31. * 24. *Thickness(I)
           TotalUptBactAmmo = TotalUptBactAmmo + UptBactAmmo *14. * 24. *Thickness(I)
!          TotalBactExcNH4  = TotalBactExcNH4  + BactExcNH4  *14. * 24. *Thickness(I)
!          TotalBactExcPO4  = TotalBactExcPO4  + BactExcPO4  *31. * 24. *Thickness(I)
           TotalMortBactC = TotalMortBactC   + MortBactC   *12. * 24. *Thickness(I)
         endif

! ALL PROCESSES BY LAYER
! EXCRETION NH4 and PO4 + EXUDATION SI04
         if(I.LE.14) THEN
           TotalBactExcNH4TOPLAYER  = TotalBactExcNH4TOPLAYER  + BactExcNH4  *14. * 24. *Thickness(I) ! mgN/m2/d
           TotalExcZooAmmoTOPLAYER  = TotalExcZooAmmoTOPLAYER + ZooExcNH4 *14. *24. * Thickness(I) ! mgN/m2/d
           TotalBactExcPO4TOPLAYER  = TotalBactExcPO4TOPLAYER  + BactExcPO4  *31. * 24. *Thickness(I) ! mgN/m2/d
           TotalExcZooPO4TOPLAYER  = TotalExcZooPO4TOPLAYER + ZooExcPO4 *31. *24. * Thickness(I) ! mgN/m2/d
           TotalExuSiTOPLAYER  = TotalExuSiTOPLAYER + ExuDiaSi *28. *24. * Thickness(I) ! mgN/m2/d

           TotalUptNitTOPLAYER  = TotalUptNitTOPLAYER + UptNit *14. *24. * Thickness(I) ! mgN/m2/d
           TotalUptAmmoTOPLAYER  = TotalUptAmmoTOPLAYER  + UptAmmo  *14. * 24. *Thickness(I) ! mgN/m2/
           TotalBactUptAmmoTOPLAYER  = TotalBactUptAmmoTOPLAYER  + UptBactAmmo  *14. * 24. *Thickness(I) ! mgN/m2/
           TotalUptPTOPLAYER  = TotalUptPTOPLAYER + UptP *31. *24. * Thickness(I) ! mgN/m2/d
           TotalBactUptPTOPLAYER  = TotalBactUptPTOPLAYER + UptBactP *31. *24. * Thickness(I) ! mgN/m2/d
           TotalUptSiTOPLAYER  = TotalUptSiTOPLAYER  + UptDiaSi  *28. * 24. *Thickness(I) ! mgN/m2/

           TotalNitrifTOPLAYER = TotalNitrifTOPLAYER + Nitrification *14. * 24. *Thickness(I) ! mgN/m2/d

           TotalRemSMOPSiTOPLAYER = TotalRemSMOPSiTOPLAYER + RemSMOPSi *28. * 24. *Thickness(I) ! mgN/m2/ 
           TotalRemLMOPSiTOPLAYER = TotalRemLMOPSiTOPLAYER + RemLMOPSi *28. * 24. *Thickness(I) ! mgN/m2/

         elseif(I.GT.14.AND.I.LE.21) THEN
           TotalBactExcNH4INT  = TotalBactExcNH4INT  + BactExcNH4  *14. * 24. *Thickness(I) ! mgN/m2/d
           TotalExcZooAmmoINT  = TotalExcZooAmmoINT + ZooExcNH4 *14. *24. * Thickness(I) ! mgN/m2/d
           TotalBactExcPO4INT  = TotalBactExcPO4INT  + BactExcPO4  *31. * 24. *Thickness(I) ! mgN/m2/d
           TotalExcZooPO4INT  = TotalExcZooPO4INT + ZooExcPO4 *31. *24. * Thickness(I) ! mgN/m2/d
           TotalExuSiINT  = TotalExuSiINT + ExuDiaSi *28. *24. * Thickness(I) ! mgN/m2/d

           TotalUptNitINT  = TotalUptNitINT + UptNit *14. *24. * Thickness(I) ! mgN/m2/d
           TotalUptAmmoINT  = TotalUptAmmoINT  + UptAmmo  *14. * 24. *Thickness(I) ! mgN/m2/d
           TotalBactUptAmmoINT  = TotalBactUptAmmoINT  + UptBactAmmo *14. * 24. *Thickness(I) ! mgN/m2/
           TotalUptPINT  = TotalUptPINT + UptP *31. *24. * Thickness(I) ! mgN/m2/d
           TotalBactUptPINT  = TotalBactUptPINT + UptBactP *31. *24. * Thickness(I) ! mgN/m2/d
           TotalUptSiINT  = TotalUptSiINT  + UptDiaSi  *28. * 24. *Thickness(I) ! mgN/m2/d

           TotalNitrifINT = TotalNitrifINT + Nitrification *14. * 24. *Thickness(I) ! mgN/m2/d

           TotalRemSMOPSiINT = TotalRemSMOPSiINT + RemSMOPSi *28. * 24. *Thickness(I) ! mgN/m2/ 
           TotalRemLMOPSiINT = TotalRemLMOPSiINT + RemLMOPSi *28. * 24. *Thickness(I) ! mgN/m2/

         else
           TotalBactExcNH4DEEP  = TotalBactExcNH4DEEP  + BactExcNH4  *14. * 24. *Thickness(I) ! mgN/m2/d
           TotalExcZooAmmoDEEP  = TotalExcZooAmmoDEEP + ZooExcNH4 *14. *24. * Thickness(I) ! mgN/m2/d
           TotalBactExcPO4DEEP  = TotalBactExcPO4DEEP  + BactExcPO4  *31. * 24. *Thickness(I) ! mgN/m2/d
           TotalExcZooPO4DEEP  = TotalExcZooPO4DEEP + ZooExcPO4 *31. *24. * Thickness(I) ! mgN/m2/d
           TotalExuSiDEEP  = TotalExuSiDEEP + ExuDiaSi *28. *24. * Thickness(I) ! mgN/m2/d

           TotalUptNitDEEP  = TotalUptNitDEEP + UptNit *14. *24. * Thickness(I) ! mgN/m2/d
           TotalUptAmmoDEEP  = TotalUptAmmoDEEP  + UptAmmo  *14. * 24. *Thickness(I) ! mgN/m2/
           TotalBactUptAmmoDEEP  = TotalBactUptAmmoDEEP  + UptBactAmmo  *14. * 24. *Thickness(I) ! mgN/m2/
           TotalUptPDEEP  = TotalUptPDEEP + UptP *31. *24. * Thickness(I) ! mgN/m2/d
           TotalBactUptPDEEP  = TotalBactUptPDEEP + UptBactP *31. *24. * Thickness(I) ! mgN/m2/d
           TotalUptSiDEEP  = TotalUptSiDEEP  + UptDiaSi  *28. * 24. *Thickness(I) ! mgN/m2/

           TotalNitrifDEEP = TotalNitrifDEEP + Nitrification *14. * 24. *Thickness(I) ! mgN/m2/d

           TotalRemSMOPSiDEEP = TotalRemSMOPSiDEEP + RemSMOPSi *28. * 24. *Thickness(I) ! mgN/m2/ 
           TotalRemLMOPSiDEEP = TotalRemLMOPSiDEEP + RemLMOPSi *28. * 24. *Thickness(I) ! mgN/m2/
         endif

! Particulate matter
        !if(I.LE.27) THEN
        TotalRemSMOPC= TotalRemSMOPC     + RemSMOPC *12. * 24. *Thickness(I)
        TotalRemSMOPN= TotalRemSMOPN     + RemSMOPN *14. * 24. *Thickness(I)
        TotalRemSMOPP= TotalRemSMOPP     + RemSMOPP *31. * 24. *Thickness(I)
        
        TotalRemLMOPC= TotalRemLMOPC     + RemLMOPC *12. * 24. *Thickness(I)
        TotalRemLMOPN= TotalRemLMOPN     + RemLMOPN *14. * 24. *Thickness(I)
        TotalRemLMOPP= TotalRemLMOPP     + RemLMOPP *31. * 24. *Thickness(I)
        
        TotalZooGrazSMOPC= TotalZooGrazSMOPC     + ZooGrazSMOPC *12. * 24. *Thickness(I)
        TotalZooGrazSMOPN= TotalZooGrazSMOPN     + ZooGrazSMOPN *14. * 24. *Thickness(I)
        TotalZooGrazSMOPP= TotalZooGrazSMOPP     + ZooGrazSMOPP *31. * 24. *Thickness(I)
        
        TotalZooGrazLMOPC= TotalZooGrazLMOPC     + ZooGrazLMOPC *12. * 24. *Thickness(I)
        TotalZooGrazLMOPN= TotalZooGrazLMOPN     + ZooGrazLMOPN *14. * 24. *Thickness(I)
        TotalZooGrazLMOPP= TotalZooGrazLMOPP     + ZooGrazLMOPP *31. * 24. *Thickness(I)
        
        TotalZooGrazBact= TotalZooGrazBact + ZooGrazBactC *12. * 24. *Thickness(I)
        
        TotalZooGrazPhytoC  = TotalZooGrazPhytoC + (ZooGrazDiaC + ZooGrazNanoC + ZooGrazSyneC ) * 12. * 24. *Thickness(I)
        TotalZooGrazPhytoN  = TotalZooGrazPhytoN + (ZooGrazDiaN + ZooGrazNanoN + ZooGrazSyneN ) * 14. * 24. *Thickness(I)
        TotalZooGrazPhytoP  = TotalZooGrazPhytoP + (ZooGrazDiaP + ZooGrazNanoP + ZooGrazSyneP ) * 31. * 24. *Thickness(I)
        
        
        TotalZooMesoPredSMOPC= TotalZooMesoPredSMOPC     + ZooMesoPredSMOPC *12. * 24. *Thickness(I)
        TotalZooMesoPredSMOPN= TotalZooMesoPredSMOPN     + ZooMesoPredSMOPN *14. * 24. *Thickness(I)
        TotalZooMesoPredSMOPP= TotalZooMesoPredSMOPP     + ZooMesoPredSMOPP *31. * 24. *Thickness(I)
        
        TotalZooMesoPredLMOPC= TotalZooMesoPredLMOPC     + ZooMesoPredLMOPC *12. * 24. *Thickness(I)
        TotalZooMesoPredLMOPN= TotalZooMesoPredLMOPN     + ZooMesoPredLMOPN *14. * 24. *Thickness(I)
        TotalZooMesoPredLMOPP= TotalZooMesoPredLMOPP     + ZooMesoPredLMOPP *31. * 24. *Thickness(I)
        
        TotalZooMicroMortSMOPC= TotalZooMicroMortSMOPC     + ZooMicroMortSMOPC *12. * 24. *Thickness(I)
        TotalZooMicroMortSMOPN= TotalZooMicroMortSMOPN     + ZooMicroMortSMOPN *14. * 24. *Thickness(I)
        TotalZooMicroMortSMOPP= TotalZooMicroMortSMOPP     + ZooMicroMortSMOPP *31. * 24. *Thickness(I)
        
        TotalZooNanoMortSMOPC= TotalZooNanoMortSMOPC     + ZooNanoMortSMOPC *12. * 24. *Thickness(I)
        TotalZooNanoMortSMOPN= TotalZooNanoMortSMOPN     + ZooNanoMortSMOPN *14. * 24. *Thickness(I)
        TotalZooNanoMortSMOPP= TotalZooNanoMortSMOPP     + ZooNanoMortSMOPP *31. * 24. *Thickness(I)
        
        
        TotalSMOPMortPhytoC = TotalSMOPMortPhytoC     + SMOPMortPhytoC *12. * 24. *Thickness(I)
        TotalSMOPMortPhytoN = TotalSMOPMortPhytoN     + SMOPMortPhytoN *14. * 24. *Thickness(I)
        TotalSMOPMortPhytoP = TotalSMOPMortPhytoP     + SMOPMortPhytoP *31. * 24. *Thickness(I)
        
        TotalLMOPMortPhytoC = TotalLMOPMortPhytoC     + LMOPMortPhytoC *12. * 24. *Thickness(I)
        TotalLMOPMortPhytoN = TotalLMOPMortPhytoN     + LMOPMortPhytoN *14. * 24. *Thickness(I)
        TotalLMOPMortPhytoP = TotalLMOPMortPhytoP     + LMOPMortPhytoP *31. * 24. *Thickness(I)
        
        TotalMortPhytoC = TotalMortPhytoC + ( SMOPMortPhytoC + LMOPMortPhytoC) *12. * 24. *Thickness(I)
        TotalMortPhytoN = TotalMortPhytoN + ( SMOPMortPhytoN + LMOPMortPhytoN) *14. * 24. *Thickness(I)
        TotalMortPhytoP = TotalMortPhytoP + ( SMOPMortPhytoP + LMOPMortPhytoP) *31. * 24. *Thickness(I)


!endif


! Output variables
         Chl(I)      = SyneChl(I) + NanoChl(I) + DiaChl(I) + S_MOPChl(I)
         NOP(I)      = S_MOPN(I) + L_MOPN(I) + ZooMicroC(I) * NCZooMicro + ZooMesoC(I) * NCZooMeso & 
                                 + ZooNanoC(I) * NCZooNano + DiaN(I) + NanoN(I) + SyneN(I)
         COP(I)      = S_MOPC(I) + L_MOPC(I) + ZooMicroC(I)              + ZooMesoC(I)             & 
                                 + ZooNanoC(I)             + DiaC(I) + NanoC(I) + SyneC(I)

! Maximum Chl and Maximum Chl Depth

         MaxDiaChl = MAX(MaxDiaChl,DiaChl(I))
         IF(MaxDiaChl .EQ. DiaChl(I)) DepthMaxDiaChl = Depth(I)  

         MaxNanoChl = MAX(MaxNanoChl,NanoChl(I))
         IF(MaxNanoChl .EQ. NanoChl(I)) DepthMaxNanoChl = Depth(I)  

         MaxSyneChl = MAX(MaxSyneChl,SyneChl(I))
         IF(MaxSyneChl .EQ. SyneChl(I)) DepthMaxSyneChl = Depth(I)  


        if (I.EQ.4) then
          SurfDiaChl = DiaChl(I)
          SurfNanoChl = NanoChl(I)
          SurfSyneChl = SyneChl(I)
        endif   

        NCDiaratio(I)    = NCDia
        NCSyneratio(I)   = NCSyne
        NCNanoratio(I)   = NCNano
        PCDiaratio(I)    = PCDia
        PCSyneratio(I)   = PCSyne
        PCNanoratio(I)   = PCNano
        ChlNDiaratio(I)  = ChlNDia
        ChlNSyneratio(I) = ChlNSyne
        ChlNNanoratio(I) = ChlNNano
        ChlCDiaratio(I)  = ChlCDia
        ChlCSyneratio(I) = ChlCSyne
        ChlCNanoratio(I) = ChlCNano
        CChlDiaratio(I)  = 1/ChlCDia  * 12. ! mgC/mgChl
        CChlSyneratio(I) = 1/ChlCSyne * 12. ! mgC/mgChl
        CChlNanoratio(I) = 1/ChlCNano * 12. ! mgC/mgChl
        NCMODratio(I)    = NCMOD
        if ((Nitrate(I)/Phosphate(I)).LE.50.)  then
          NPDI(I) =Nitrate(I) / Phosphate(I) 
        else
          NPDI(I) = 0.
        endif
        IF(MessyFeedCMOD > 0.) NCMessyFeedMOD(I)= MessyFeedNMOD / MessyFeedCMOD
        IF(ExuC > 0.)  NCExu(I) = ExuN / ExuC
        IF(ExuSyneC > 0.)  NCExuSyne(I) = ExuSyneN / ExuSyneC
        IF(ExuNanoC > 0.)  NCExuNano(I) = ExuNanoN / ExuNanoC
        IF(ExuDiaC > 0.)   NCExuDia(I)  = ExuDiaN  / ExuDiaC
        IF(RemSMOPC > 0) NCRemSMOP(I) = RemSMOPN /RemSMOPC
        IF(RemLMOPC > 0) NCRemLMOP(I) = RemLMOPN /RemLMOPC
        
        if(DiaN(I) > 0.) then
          CNDiaratio(I)    = DiaC(I) * 12. / (DiaN(I) * 14.)
        else
          CNDiaratio(I) = 0.
        endif
        IF(NanoN(I) > 0.) then
          CNNanoratio(I)    = NanoC(I) * 12. / (NanoN(I) * 14.)
        else
          CNNanoratio(I) = 0.
        endif
        IF(SyneN(I) > 0.) then
          CNSyneratio(I)    = SyneC(I) * 12. / (SyneN(I) * 14.)
        else
          CNSyneratio(I) = 0.
        endif

! Biomass Phyto and Zoo integrated 0-200m
!   if(I.LE.27) then
         IntSyne  = IntSyne   +  SyneChl(I)                           * Thickness(I)
         IntNano  = IntNano   +  NanoChl(I)                           * Thickness(I)
         IntDia   = IntDia    +   DiaChl(I)                           * Thickness(I)
         PhytoChl = PhytoChl  + (SyneChl(I) + NanoChl(I) + DiaChl(I)) * Thickness(I)
!endif

!   if(I.LE.27) then
         IntNanoZoo  = IntNanoZoo  +   ZooNanoC(I)       * Thickness(I)
         IntMicroZoo = IntMicroZoo +   ZooMicroC(I)      * Thickness(I)
         IntMesoZoo  = IntMesoZoo  +   ZooMesoC(I)       * Thickness(I)
!ZooC        = ZooC        + (zooNanoC(I) + ZooMicroC(I) + ZooMesoC(I) ) * Thickness(I)
!endif

!if(I.LE.27) then

         if(Phosphate(I).GT.0.AND.(Nitrate(I)/Phosphate(I)).GT.50.) then
           IntNitrate   = IntNitrate   +  Nitrate(I)     * Thickness(I)
           IntPhosphate = IntPhosphate +  Phosphate(I)   * Thickness(I)
         endif

!endif
     
       endif  ! mask debut de boucle sur box  
     enddo !->NumpelagicBoxes ! with boxes

     if(IntPhosphate.GT.0.) then 
       IntNPDI = IntNitrate / IntPhosphate
     else
       IntNPDI = 0.
     endif 

! Quantity of surface nutrient


! Attention ces lignes doivent etre deplacees dans la boucle DO I -->  ENDDO (au dessus...)
!      if(I.LE.27) then
!        TotalNitrateSurf   = TotalNitrateSurf   +  Nitrate(I)     * Thickness(I)
!        TotalPhosphateSurf = TotalPhosphateSurf +  Phosphate(I)   * Thickness(I)
!        TotalSiliceSurf    = TotalSiliceSurf    +  Silice(I)      * Thickness(I)
!      endif


    CALL Budget

   RETURN
   END SUBROUTINE DynamicsEco3m

!**********************************************************************
!  END OF SIMULATION
!**********************************************************************

   SUBROUTINE FinaliseEco3m
   END SUBROUTINE FinaliseEco3m



!**********************************************************************
!  FUNCTIONS
!**********************************************************************



   DOUBLE PRECISION FUNCTION  f_pproduction (AbsChl,PhyMax,Sig_Ps2,EPAR,kd,kr,Tau_renew,Gamma_M,Tempfacppb) ! Gamma_M,Tempfacppb

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT (IN) :: AbsChl,Phymax,Sig_Ps2,EPAR,kd,kr,Tau_renew,Gamma_M,Tempfacppb ! Gamma_M,Tempfacppb
   DOUBLE PRECISION :: P0
   ! return value f_pproduction s in mmolC/mg Chl
   ! Test Alex 20/04/2016 suppression de partie !!! 24/05/2016 * 1.55 sans arrondi pour Phymax ou AbsChl 
 

   p0 =  1.D0 / ( 1 + Sig_PS2 * EPAR * Tau_renew + kd/kr * (Sig_PS2 * EPAR)**2. * Tau_renew ) 
    
   ! Alex 22/09/2016 ajout à nouveau du boost + 30% (* 1.3) pour simu 3+10 ans 
   ! Alex 03/10/2016 retrait du boost + 30% (* 1.3) pour test fuite atlantique
   ! valeurs negatives 
   f_pproduction = Tempfacppb * Gamma_M * AbsChl * Phymax * p0 * EPAR * 3600. ! * 1.3                 
!  f_pproduction =              Gamma_M * AbsChl * Phymax * p0 * EPAR * 3600.

   RETURN
   END FUNCTION f_pproduction


   SUBROUTINE Budget

   USE ModuleComBlock

   INTEGER :: I

   TotalN   = 0.D0
   TotalP   = 0.D0
   TotalSi  = 0.D0

   do I = 1,NumPelagicBoxes

    TotalN   = TotalN  + (ZooNanoC(I) * NCZooNano +ZooMicroC(I) * NCZooMicro + ZooMesoC(I) * NCZooMeso & 
          +BactC(I) * NCBact + DiaN(I)+NanoN(I)+SyneN(I)+S_MOPN(I)+L_MOPN(I)+MODN(I)+Nitrate(I)+Ammonium(I))*Thickness(I)
    TotalP   = TotalP  + (ZooNanoC(I) * PCZooNano +ZooMicroC(I) * PCZooMicro + ZooMesoC(I) * PCZooMeso & 
          +BactC(I) * PCBact + DiaP(I)+NanoP(I)+SyneP(I)+S_MOPP(I)+L_MOPP(I)+MODP(I)+Phosphate(I))*Thickness(I)
    TotalSi  = TotalSi + (DiaSi(I)+S_MOPSi(I)+L_MOPSi(I)+Silice(I))*Thickness(I)

   enddo


   END SUBROUTINE Budget
