      SUBROUTINE Eco3mParameters
!______________________________________________________________________
! 3d ecosystem model
!
! LAST REVISION: 3 MAY 2007
!
! Implementation: Caroline Ulses
!                 NIOO-CEME
!
! 
!______________________________________________________________________


!______________________________________________________________________
!
! 1. Reads the pelagic parameters in notebooks
!   1.a Zooplankton
!   1.b Phytoplankton
!   1.c Bacteria
!   1.d Remineralisation
!   1.e Oxygen
!   1.f Bio (sinking velocity)
! 2. Converts in hours
! 3. Writes them in an output file 
!
!______________________________________________________________________

!______________________________________________________________________
!
! Modifications:
! 21/02/07: include symphonie2007.h
! 10/09/08: remplacer symphonie2007.h par MODULE PRINCIPAL
! 28/11/12: lecture du notebook_oxygen - Faycal KESSOURI
! 16/02/18: lecture du notebook_bio - Alex
!______________________________________________________________________
      USE ModuleComBlock
      USE MODULE_PRINCIPAL
      use module_parallele
      use ModuleDeclaration
      IMPLICIT NONE

      INTEGER LDIR

!======================================================================
! 1.a Reading of the zooplankton parameter notebook: notebook_zooplankton
! BEGINNING:
!======================================================================


      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read notebook_zooplankton'     & 
                        ,NOMFICHIER(24)
      
      OPEN(UNIT=3,FILE=NOMFICHIER(24)) ! open the notebook
      
      DO K=1,11
      READ(3,*)
      ENDDO

      READ(3,*)GrazNano      
      READ(3,*)PrefNanoSyne      
      READ(3,*)PrefNanoNano       
      READ(3,*)PrefNanoZooNano
      READ(3,*)kZooNano          

      READ(3,*)

      READ(3,*)GrazMicro
      READ(3,*)PrefMicroBact
      READ(3,*)PrefMicroSyne
      READ(3,*)PrefMicroNano
      READ(3,*)PrefMicroDia
      READ(3,*)PrefMicroZooMicro
      READ(3,*)PrefMicroZooNano
      READ(3,*)kZooMicro

      READ(3,*)

      READ(3,*)GrazMeso
      READ(3,*)PrefMesoDia
      READ(3,*)PrefMesoZooMicro
      READ(3,*)PrefMesoLMOP
      READ(3,*)kZooMeso

      READ(3,*)

      READ(3,*)kcNetGrowthEff

      READ(3,*)
      READ(3,*)

      READ(3,*)MessyFeedingfracC
      READ(3,*)MessyFeedingfracN
      READ(3,*)MessyFeedingfracP

      READ(3,*)
      READ(3,*)

      READ(3,*)NanoAssEffBact  
      READ(3,*)NanoAssEffSyne
      READ(3,*)NanoAssEffNano
      READ(3,*)NanoAssEffZooNano

      READ(3,*)

      READ(3,*)MicroAssEffBact
      READ(3,*)MicroAssEffSyne
      READ(3,*)MicroAssEffNano
      READ(3,*)MicroAssEffDia
      READ(3,*)MicroAssEffZooNano
      READ(3,*)MicroAssEffZooMicro
      READ(3,*)MicroAssEffSMOP

      READ(3,*)

      READ(3,*)MesoAssEffDia
      READ(3,*)MesoAssEffZooMicro
      READ(3,*)MesoAssEffSMOP
      READ(3,*)MesoAssEffLMOP
      
      READ(3,*)

      READ(3,*)DiaPropMOPSi

      READ(3,*)
      READ(3,*)

      READ(3,*)NCZooNano        
      READ(3,*)PCZooNano

      READ(3,*)

      READ(3,*)NCZooMicro
      READ(3,*)PCZooMicro

      READ(3,*)

      READ(3,*)NCZooMeso
      READ(3,*)PCZooMeso

      READ(3,*)
      READ(3,*)

      READ(3,*)NanoMortRate   
      READ(3,*)MicroMortRate
      READ(3,*)MesoPredRate

      READ(3,*)

      READ(3,*)ZooPropMOPC     
      READ(3,*)ZooPropMOPN     
      READ(3,*)ZooPropMOPP     

      READ(3,*)

      READ(3,*)ZooPropMesoMOPC
      READ(3,*)ZooPropMesoMOPN
      READ(3,*)ZooPropMesoMOPP
      READ(3,*)Kclos

      READ(3,*)

      READ(3,*)Q10_Zoo
      READ(3,*)TREF_Zoo

      CLOSE(3)

      prefNanoBact   = 1.D0 - prefNanoSyne  - prefNanoNano      &
                            - prefNanoZooNano

      prefMicroSMOP  = 1.D0 - prefMicroBact - prefMicroSyne     &
                            - prefMicroNano - prefMicroDia      &
                            - prefMicroZooNano - prefMicroZooMicro   

      prefMesoSMOP   = 1.D0 - prefMesoDia   - prefMesoZooMicro  &
                            - prefMesoLMOP

      if(prefNanoBact .lt. 0.D0) then
         write(6,*)'Error: prefNanoBact <0'
      endif
      if(prefMicroSMOP .lt. 0.D0) then
         write(6,*)'Error: prefMicroSMOP <0'
      endif
      if(prefMesoSMOP .lt. 0.D0) then
         write(6,*)'Error: prefMesoSMOP <0'
      endif


!======================================================================
! 1.a Reading of the zooplankton parameter notebook: notebook_zooplankton
! END:
!======================================================================      


!======================================================================
! 1.b Reading of the phytoplankton parameter notebook: notebook_phytoplankton
! BEGINNING:
!======================================================================



      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read notebook_phytoplankton'     & 
                        ,NOMFICHIER(25)

      OPEN(UNIT=3,FILE=NOMFICHIER(25)) ! open the notebook

      DO K=1,11
      READ(3,*)
      ENDDO

      READ(3,*)PhyMaxSyne
      READ(3,*)PhyMaxNano
      READ(3,*)PhyMaxDia
      READ(3,*)AbsChlSyne
      READ(3,*)AbsChlNano
      READ(3,*)AbsChlDia
      READ(3,*)Tau_renewSyne
      READ(3,*)Tau_renewNano
      READ(3,*)Tau_renewDia
      READ(3,*)Sig_PS2Syne
      READ(3,*)Sig_PS2Nano
      READ(3,*)Sig_PS2Dia
      READ(3,*)kd
      READ(3,*)krep
      READ(3,*)NGeider

      READ(3,*)

      READ(3,*)minNCSyne
      READ(3,*)maxNCSyne
      READ(3,*)minPCSyne
      READ(3,*)maxPCSyne

      READ(3,*)

      READ(3,*)minNCNano
      READ(3,*)maxNCNano
      READ(3,*)minPCNano
      READ(3,*)maxPCNano

      READ(3,*)

      READ(3,*)minNCDia
      READ(3,*)maxNCDia
      READ(3,*)minPCDia
      READ(3,*)maxPCDia
      READ(3,*)minSiCDia
      READ(3,*)maxSiCDia

      READ(3,*)
      READ(3,*)

      READ(3,*)RespCostSyne
      READ(3,*)RespCostNano
      READ(3,*)RespCostDia

      READ(3,*)
      READ(3,*)

      READ(3,*)ChlNSyneMax
      READ(3,*)ChlNNanoMax
      READ(3,*)ChlNDiaMax

      READ(3,*)
      READ(3,*)

      READ(3,*)KNitrateSyne
      READ(3,*)KNitrateNano
      READ(3,*)KNitrateDia
      READ(3,*)KAmmoSyne
      READ(3,*)KAmmoNano
      READ(3,*)KAmmoDia
      READ(3,*)KPSyne
      READ(3,*)KPNano
      READ(3,*)KPDia
      READ(3,*)KSiDia
      READ(3,*)KInhibAmmo
      READ(3,*)InhibAmmo
      READ(3,*)fRespNit
      READ(3,*)fRespAmmo
      READ(3,*)fRespP
      READ(3,*)fRespSi

      READ(3,*)
      READ(3,*)

      READ(3,*)MortSyne
      READ(3,*)MortNano
      READ(3,*)MortDia

      READ(3,*)

      READ(3,*)PropDia

      READ(3,*)
      READ(3,*)

      READ(3,*)Q10_P
      READ(3,*)TREF_P



      CLOSE(3)

!=======================================================================
! 1.b Reading of the phytoplankton parameter notebook: notebook_phytoplankton
! END:
!======================================================================

!======================================================================
! 1.c Reading of the bacteria parameter notebook: notebook_bacteria    
! BEGINNING:
!======================================================================


      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read notebook_bacteria'         &
                        ,NOMFICHIER(26)

      OPEN(UNIT=3,FILE=NOMFICHIER(26)) ! open the notebook

      DO K=1,10
      READ(3,*)
      ENDDO

      READ(3,*)UptmaxBact
      READ(3,*)KMODCBact
      READ(3,*)KAmmoBact
      READ(3,*)KPBact
      READ(3,*)NCBact
      READ(3,*)PCBact
      READ(3,*)EffGrowthBact
      READ(3,*)MortBactRate
      READ(3,*)Q10_B
      READ(3,*)TREF_B

      CLOSE(3) 

!======================================================================
! 1.c Reading of the bacteria parameter notebook: notebook_bacteria
! END:
!======================================================================

!=======================================================================
! 1.d Reading of the remineralisation parameter notebook: notebook_remineralisation
! BEGINNING:
!======================================================================



      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read noteook_remineralisation'    &
                        ,NOMFICHIER(27)

      OPEN(UNIT=3,FILE=NOMFICHIER(27)) ! open the notebook

      DO K=1,11
      READ(3,*)
      ENDDO

      READ(3,*)RemRateSMOPC
      READ(3,*)RemRateSMOPN
      READ(3,*)RemRateSMOPP
      READ(3,*)RemRateSMOPSi
      READ(3,*)RemRateSMOPChl

      READ(3,*)

      READ(3,*)RemRateLMOPC
      READ(3,*)RemRateLMOPN
      READ(3,*)RemRateLMOPP
      READ(3,*)RemRateLMOPSi

      READ(3,*)

      READ(3,*)Q10_Rem
      READ(3,*)TREF_Rem

      READ(3,*)
      READ(3,*)

      READ(3,*)Q10_Nit
      READ(3,*)NitriRate
      READ(3,*)Gamma_M
      READ(3,*)TREF_Nit

      CLOSE(3)

!=======================================================================
! 1.d Reading of the remineralisation parameter notebook: notebook_remineralisation
! END:
!======================================================================



!======================================================================
! 1.e Reading of the oxygen parameter notebook: notebook_oxygen    
! BEGINNING:
!======================================================================


      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read notebook_oxygen'         &
                        ,NOMFICHIER(31)

      OPEN(UNIT=3,FILE=NOMFICHIER(31)) ! open the notebook

      DO K=1,13
      READ(3,*)
      ENDDO

      READ(3,*) ICHOIXFS
      ichoixFSCO2=ICHOIXFS
      ichoixK1K2C=2
      READ(3,*)
      READ(3,*)

      READ(3,*)nuOXYGENNHs
      READ(3,*)nuNITNHs
      READ(3,*)nuOXYGENODU
      READ(3,*)nuNITODU
      READ(3,*)

      READ(3,*)kNHsoxidSatOX
      READ(3,*)kNHsoxidNOsInOX
      READ(3,*)kNHsoxidNOsSatNOs
      READ(3,*)

      READ(3,*)kODUoxidSatOX
      READ(3,*)kODUoxidNOsInOX
      READ(3,*)kODUoxidNOsSatNOs
      READ(3,*)

      READ(3,*)kOxicMinSatOX
      READ(3,*)

      READ(3,*)kDenitInOX
      READ(3,*)kDenitSatNOs
      READ(3,*)

      READ(3,*)kAnoxMinInOX
      READ(3,*)kAnoxMinInNOs
      READ(3,*)

      READ(3,*)kSolFormSatIron
      READ(3,*)

      READ(3,*)Rsolid
      READ(3,*)

      READ(3,*)gammaC_OX
      READ(3,*)gammaNHs_OX
      READ(3,*)gammaODU_OX
      READ(3,*)

      READ(3,*)gammaPOC_ODU
      READ(3,*)

      READ(3,*)gammaPOC_NOs
      READ(3,*)gammaNHs_NOs
      READ(3,*)gammaODU_NOs
      READ(3,*)

      READ(3,*)iron
      do i=1,6
      READ(3,*)
      enddo

      READ(3,*)ksatMortZoo
      READ(3,*)AnoxMort
      do i=1,3
      READ(3,*)
      enddo

      READ(3,*)ksatOxicHydrol


      CLOSE(3)

!=======================================================================
! 1.f Reading of the bio parameter notebook:
! notebook_bio
! BEGINNING:
!======================================================================



!      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read noteook_bio'&
!                        ,NOMFICHIER(12)
!
!      OPEN(UNIT=3,FILE=NOMFICHIER(12)) ! open the notebook
!
!      DO K=1,17 ! doesn't read anything on lines 1 to 17
!      READ(3,*)
!      ENDDO
!      
!      DO vb=18,32
!      READ(3,*)wsed(1,vb)
!      ENDDO
!     
! 
!      
!      CLOSE(3)

!=======================================================================
! 1.f Reading of the bio parameter notebook:
! notebook_bio
! END:
!======================================================================


!=======================================================================
! 2. Conversion in hours
! BEGINNING:
!======================================================================
     
! Zooplankton dynamics

      GrazNano      = GrazNano      * 3600.
      GrazMicro     = GrazMicro     * 3600.
      GrazMeso      = GrazMeso      * 3600. 
      NanoMortRate  = NanoMortRate  * 3600.
      MicroMortRate = MicroMortRate * 3600.
      MesoPredRate  = MesoPredRate  * 3600.


! Phytoplankton dynamics

      MortSyne = MortSyne * 3600.
      MortNano = MortNano * 3600.
      MortDia  = MortDia  * 3600.


! Bacteria dynamics

      UptmaxBact   = UptmaxBact   * 3600.
      MortBactRate = MortBactRate * 3600.


! Remineralisation

      RemRateSMOPC   = RemRateSMOPC   * 3600.
      RemRateSMOPN   = RemRateSMOPN   * 3600.
      RemRateSMOPP   = RemRateSMOPP   * 3600.
      RemRateSMOPSi  = RemRateSMOPSi  * 3600.
      RemRateSMOPChl = RemRateSMOPChl * 3600.
      RemRateLMOPC   = RemRateLMOPC   * 3600.
      RemRateLMOPN   = RemRateLMOPN   * 3600.
      RemRateLMOPP   = RemRateLMOPP   * 3600.
      RemRateLMOPSi  = RemRateLMOPSi  * 3600.
      NitriRate      = NitriRate      * 3600.

!=======================================================================
! 2. Conversion in hours
! END:
!======================================================================


!======================================================================
! 3. Writing of the pelagic parameters in an output file
! BEGINNING :
!======================================================================

      if(par%rank==0)print*,'Avant nomfichier Eco3m Parameters'
      TEXTE90=nomfichier(21)
      OPEN(UNIT=3,FILE=TEXTE90)
      READ(3,*)
      READ(3,*)
      READ(3,'(A)')DIRGRAPH                                            !13/04/06
      CLOSE(3)

      DO 20 K=1,90
        IF(DIRGRAPH(K:K).EQ.' ') THEN
        LDIR=K-1
        GOTO 21
        ENDIF
   20 CONTINUE
   21 CONTINUE

      TEXTE90=DIRGRAPH(1:LDIR)//'/pelagic_parameters'
      if(par%rank==0)WRITE(6,*)'Writing of the parameters in the directory'           &
                ,TEXTE90 

      if(par%rank==0)then
      OPEN(UNIT=2,FILE=TEXTE90)

      WRITE(2,*)'======================='
      WRITE(2,*)' Zooplankton dynamics  '
      WRITE(2,*)'======================='

      WRITE(2,*)

      WRITE(2,*)'Grazing'
      WRITE(2,*)'-------'
      WRITE(2,*)'GrazNano         =    ',GrazNano/3600.  ,'    /s'     &
                                        ,GrazNano        ,'    /hr'
      WRITE(2,*)'PrefNanoSyne     =    ',PrefNanoSyne	 ,'    - '
      WRITE(2,*)'PrefNanoNano     =    ',PrefNanoNano    ,'     - '
      WRITE(2,*)'PrefNanoZooNano  =    ',PrefNanoZooNano ,'     - '
      WRITE(2,*)'PrefNanoBact     =    ',PrefNanoBact    ,'     - '
      WRITE(2,*)'kZooNano         =    ',kZooNano        ,'   mmolC/m3'

      WRITE(2,*)

      WRITE(2,*)'GrazMicro         =    ',GrazMicro/3600.  ,'     /s'  &
                                         ,GrazMicro        ,'     /hr'
      WRITE(2,*)'PrefMicroBact     =    ',PrefMicroBact    ,'     - '
      WRITE(2,*)'PrefMicroSyne     =    ',PrefMicroSyne    ,'     - '
      WRITE(2,*)'PrefMicroNano     =    ',PrefMicroNano    ,'     - '
      WRITE(2,*)'PrefMicroDia      =    ',PrefMicroDia     ,'     - '
      WRITE(2,*)'PrefMicroZooMicro =    ',PrefMicroZooMicro,'     - '
      WRITE(2,*)'PrefMicroZooNano  =    ',PrefMicroZooNano ,'     - '
      WRITE(2,*)'PrefMicroSMOP     =    ',PrefMicroSMOP    ,'     - '
      WRITE(2,*)'kZooMicro         =    ',kZooMicro        ,'  mmolC/m3'

      WRITE(2,*)

      WRITE(2,*)'GrazMeso         =    ',GrazMeso/3600.    ,'     /s'  &
                                        ,GrazMeso          ,'     /hr'
      WRITE(2,*)'PrefMesoDia      =    ',PrefMesoDia       ,'     - '
      WRITE(2,*)'PrefMesoZooMicro =    ',PrefMesoZooMicro  ,'     - '
      WRITE(2,*)'PrefMesoSMOP     =    ',PrefMesoSMOP      ,'     - '
      WRITE(2,*)'PrefMesoLMOP     =    ',PrefMesoLMOP      ,'     - '
      WRITE(2,*)'kZooMeso         =    ',kZooMeso          ,'  mmolC/m3'

      WRITE(2,*)

      WRITE(2,*)'kcNetGrowthEff  =    ',kcNetGrowthEff   ,'     - '

      WRITE(2,*)

      WRITE(2,*)'Messy Feeding'
      WRITE(2,*)'-------------'
      WRITE(2,*)'MessyFeedingfracC =    ',MessyFeedingfracC,'     - '
      WRITE(2,*)'MessyFeedingfracN =    ',MessyFeedingfracN,'     - '
      WRITE(2,*)'MessyFeedingfracP =    ',MessyFeedingfracP,'     - '

      WRITE(2,*)

      WRITE(2,*)'Egestion'
      WRITE(2,*)'--------'
      WRITE(2,*)'NanoAssEffBact      =     ',NanoAssEffBact    ,'  -'
      WRITE(2,*)'NanoAssEffSyne      =     ',NanoAssEffSyne    ,'  -'
      WRITE(2,*)'NanoAssEffNano      =     ',NanoAssEffNano    ,'  -'
      WRITE(2,*)'NanoAssEffZooNano   =     ',NanoAssEffZooNano ,'  -'

      WRITE(2,*)

      WRITE(2,*)'MicroAssEffBact      =     ',MicroAssEffBact    ,' -'
      WRITE(2,*)'MicroAssEffSyne      =     ',MicroAssEffSyne    ,' -'
      WRITE(2,*)'MicroAssEffNano      =     ',MicroAssEffNano    ,' -'
      WRITE(2,*)'MicroAssEffDia       =     ',MicroAssEffDia     ,' -'
      WRITE(2,*)'MicroAssEffZooNano   =     ',MicroAssEffZooNano ,' -'
      WRITE(2,*)'MicroAssEffSMOP      =     ',MicroAssEffSMOP    ,' -'
      WRITE(2,*)'MicroAssEffZooMicro  =     ',MicroAssEffZooMicro,' -'

      WRITE(2,*)

      WRITE(2,*)'MesoAssEffDia      =     ',MesoAssEffDia     	,' -'
      WRITE(2,*)'MesoAssEffZooMicro =     ',MesoAssEffZooMicro  ,' -'
      WRITE(2,*)'MesoAssEffSMOP     =     ',MesoAssEffSMOP    	,' -'
      WRITE(2,*)'MesoAssEffLMOP     =     ',MesoAssEffLMOP 	,' -'

      WRITE(2,*)

      WRITE(2,*)'DiaPropMOPSi    =     ',DiaPropMOPSi  ,'	-'

      WRITE(2,*)

      WRITE(2,*)'Excretion'
      WRITE(2,*)'---------'
      WRITE(2,*)'NCZooNano     =       ',NCZooNano  ,'    mmolN/mmolC'  
      WRITE(2,*)'PCZooNano     =       ',PCZooNano  ,'    mmolP/mmolC'  

      WRITE(2,*)

      WRITE(2,*)'NCZooMicro     =       ',NCZooMicro  ,' mmolN/mmolC'
      WRITE(2,*)'PCZooMicro     =       ',PCZooMicro  ,' mmolP/mmolC'

      WRITE(2,*)

      WRITE(2,*)'NCZooMeso     =       ',NCZooMeso  ,'   mmolN/mmolC'
      WRITE(2,*)'PCZooMeso     =       ',PCZooMeso  ,'   mmolP/mmolC'

      WRITE(2,*)

      WRITE(2,*)'Mortality'
      WRITE(2,*)'---------'
      WRITE(2,*)'NanoMortRate     =',NanoMortRate/3600.,'m3/s/mmolC'   &
                                    ,NanoMortRate      ,'m3/hr/mmolC'
      WRITE(2,*)'MicroMortRate    =',MicroMortRate/3600.,'m3/s/mmolC'  &
                                    ,MicroMortRate      ,'m3/hr/mmolC'
      WRITE(2,*)'MesoPredRate     =',MesoPredRate/3600.,'m3/s/mmolC'   &
                                    ,MesoPredRate      ,'m3/hr/mmolC'
      WRITE(2,*)'ZooPropMOPC  =       ',ZooPropMOPC ,'	-'
      WRITE(2,*)'ZooPropMOPN  =       ',ZooPropMOPN ,'	-'
      WRITE(2,*)'ZooPropMOPP  =       ',ZooPropMOPP ,'	-'

      WRITE(2,*)'ZooPropMesoMOPC  =       ',ZooPropMesoMOPC ,'  -'
      WRITE(2,*)'ZooPropMesoMOPN  =       ',ZooPropMesoMOPN ,'  -'
      WRITE(2,*)'ZooPropMespMOPP  =       ',ZooPropMesoMOPP ,'  -'
      WRITE(2,*)'Kclos		=	  ',Kclos           ,'-'

      WRITE(2,*)

      WRITE(2,*)'Q10_Zoo   =       ',Q10_Zoo ,'  /C'
      WRITE(2,*)'TREF_Zoo  =       ',TREF_Zoo ,'  C'


      WRITE(2,*)
      WRITE(2,*)

      WRITE(2,*)'======================='
      WRITE(2,*)'Phytoplankton dynamics '
      WRITE(2,*)'======================='

      WRITE(2,*)'PhyMaxSyne	=',PhyMaxSyne   ,'	mmolC/J'
      WRITE(2,*)'PhyMaxNano     =',PhyMaxNano   ,'       mmolC/J'
      WRITE(2,*)'PhyMaxDia	=',PhyMaxDia    ,'	mmolC/J'
      WRITE(2,*)'AbsChlSyne	=',AbsChlSyne   ,'	m2/mgChl'
      WRITE(2,*)'AbsChlNano     =',AbsChlNano   ,'       m2/mgChl'
      WRITE(2,*)'AbsChlDia	=',AbsChlDia    ,'	m2/mgChl'
      WRITE(2,*)'Tau_renewSyne	=',Tau_renewSyne,'	sec'
      WRITE(2,*)'Tau_renewNano  =',Tau_renewNano,'       sec'
      WRITE(2,*)'Tau_renewDia	=',Tau_renewDia ,'	sec'
      WRITE(2,*)'Sig_PS2Syne	=',Sig_PS2Syne  ,'	m2/J'
      WRITE(2,*)'Sig_PS2Nano    =',Sig_PS2Nano  ,'       m2/J' 		
      WRITE(2,*)'Sig_PS2Dia	=',Sig_PS2Dia   ,'	m2/J'
      WRITE(2,*)'kd		=',kd	       ,'	-'	
      WRITE(2,*)'krep		=	',krep	       ,'	-'	
      WRITE(2,*)'NGeider        =       ',NGeider         ,'       -'

      WRITE(2,*)

      WRITE(2,*)'minNCSyne	=	',minNCSyne,'	mmolN/mmolC'
      WRITE(2,*)'maxNCSyne	=	',maxNCSyne,'	mmolN/mmolC'
      WRITE(2,*)'minPCSyne	=	',minPCSyne,'	mmolP/mmolC'
      WRITE(2,*)'maxPCSyne	=	',maxPCSyne,'	mmolP/mmolC'

      WRITE(2,*)

      WRITE(2,*)'minNCNano      =       ',minNCNano,'   mmolN/mmolC'
      WRITE(2,*)'maxNCNano      =       ',maxNCNano,'   mmolN/mmolC'
      WRITE(2,*)'minPCNano      =       ',minPCNano,'   mmolP/mmolC'
      WRITE(2,*)'maxPCNano      =       ',maxPCNano,'   mmolP/mmolC'

      WRITE(2,*)

      WRITE(2,*)'minNCDia	=	',minNCDia,'	mmolN/mmolC'
      WRITE(2,*)'maxNCDia	=	',maxNCDia,'	mmolN/mmolC'
      WRITE(2,*)'minPCDia	=	',minPCDia,'	mmolP/mmolC'
      WRITE(2,*)'maxPCDia	=	',maxPCDia,'	mmolP/mmolC'
      WRITE(2,*)'minSiCDia	=	',minSiCDia,'	mmolSi/mmolC'
      WRITE(2,*)'maxSiCDia	=	',maxSiCDia,'	mmolsi/mmolC'

      WRITE(2,*)

      WRITE(2,*)'RespCostSyne  	=	',RespCostSyne,'	mmolC'
      WRITE(2,*)'RespCostNano   =       ',RespCostNano,'        mmolC'
      WRITE(2,*)'RespCostDia	=	',RespCostDia ,'	mmolC'

      WRITE(2,*)

      WRITE(2,*)'ChlNSyneMax	=	',ChlNSyneMax,'	mgChl/mmolN'
      WRITE(2,*)'ChlNNanoMax    =       ',ChlNNanoMax,' mgChl/mmolN'
      WRITE(2,*)'ChlNDiaMax	=	',ChlNDiaMax ,'	mgChl/mmolN'

      WRITE(2,*)

      WRITE(2,*)'KNitrateSyne	=	',KNitrateSyne,'mmolN/m3'
      WRITE(2,*)'KNitrateNano   =       ',KNitrateNano,'mmolN/m3'
      WRITE(2,*)'KNitrateDia	=	',KNitrateDia ,'mmolN/m3'
      WRITE(2,*)'KAmmoSyne	=	',KAmmoSyne   ,'mmolN/m3'
      WRITE(2,*)'KAmmoNano      =       ',KAmmoNano   ,'mmolN/m3'
      WRITE(2,*)'KAmmoDia	=	',KAmmoDia    ,'mmolN/m3'
      WRITE(2,*)'KPSyne		=	',KPSyne      ,'mmolP/m3'
      WRITE(2,*)'KPNano         =       ',KPNano      ,'mmolP/m3'
      WRITE(2,*)'KPDia		=	',KPDia	      ,'mmolP/m3'
      WRITE(2,*)'KSiDia		=	',KSiDia      ,'mmolSi/m3'
      WRITE(2,*)'KInhibAmmo	=	',KInhibAmmo  ,'mmolN/m3'
      WRITE(2,*)'InhibAmmo	=	',InhibAmmo   ,'-'
      WRITE(2,*)'fRespNit	=	',fRespNit    ,'molC/molN'	
      WRITE(2,*)'fRespAmmo	=	',fRespAmmo   ,'molC/molN'
      WRITE(2,*)'fRespP		=	',fRespP      ,'molP/molN'
      WRITE(2,*)'fRespSi	=	',fRespSi     ,'molSi/molN'


      WRITE(2,*)'MortSyne      	=       ',MortSyne   ,'molC/molN'
      WRITE(2,*)'MortNano       =       ',MortNano      ,'molP/molN'
      WRITE(2,*)'MortDia        =       ',MortDia     ,'molSi/molN'

      WRITE(2,*)

      WRITE(2,*)'PropDia	=	',PropDia     ,'-'

      WRITE(2,*)

      WRITE(2,*)'Q10_P		=	',Q10_P       ,'/degC'
      WRITE(2,*)'TREF_P          =       ',TREF_P       ,'degC'


      WRITE(2,*)'=========='
      WRITE(2,*)' Bacteria '
      WRITE(2,*)'=========='

      WRITE(2,*)'UptmaxBact   =       ',UptmaxBact/3600.,'/s'          &
                                       ,UptmaxBact      ,'/hr'
      WRITE(2,*)'KMODCBact    =       ',KMODCBact       ,'mmolC/m3'
      WRITE(2,*)'KAmmoBact    =       ',KAmmoBact       ,'mmolN/m3'
      WRITE(2,*)'KPBact       =       ',KPBact       	,'mmolP/m3'
      WRITE(2,*)'NCBact       =       ',NCBact          ,'molN/molC'
      WRITE(2,*)'PCBact       =       ',PCBact          ,'molP/molC'
      WRITE(2,*)'EffGrowthBact=       ',EffGrowthBact   ,'-'
      WRITE(2,*)'MortBactRate =       ',MortBactRate/3600. ,'/s'       &
                                       ,MortBactRate      ,'/hr'
      WRITE(2,*)'Q10_B   	=       ',Q10_B ,'  /C'
      WRITE(2,*)'TREF_B  	=       ',TREF_B ,'  C'


      WRITE(2,*)

      WRITE(2,*)'======================='
      WRITE(2,*)' Remineralisation '
      WRITE(2,*)'======================='

      WRITE(2,*)'RemRateSMOPC	=	',RemRateSMOPC/3600.,'	/s'     &
                                         ,RemRateSMOPC,'	/hr'
      WRITE(2,*)'RemRateSMOPN	=	',RemRateSMOPN/3600.,'  /s'     &
                                         ,RemRateSMOPN,'	/hr'
      WRITE(2,*)'RemRateSMOPP	=	',RemRateSMOPP/3600.,'  /s'     &
                                         ,RemRateSMOPP,'	/hr'
      WRITE(2,*)'RemRateSMOPSi	=	',RemRateSMOPSi/3600.,' /s'     &
                                         ,RemRateSMOPSi,'	/hr'
      WRITE(2,*)'RemRateSMOPChl	=	',RemRateSMOPChl/3600.,'/s'     &
                                         ,RemRateSMOPChl,'	/hr'

      WRITE(2,*)

      WRITE(2,*)'RemRateLMOPC	=	',RemRateLMOPC/3600.,'  /s'     &
                                         ,RemRateLMOPC ,'	/hr'
      WRITE(2,*)'RemRateLMOPN	=	',RemRateLMOPN/3600.,'  /s'     &
                                         ,RemRateLMOPN ,'	/hr'
      WRITE(2,*)'RemRateLMOPP	=	',RemRateLMOPP/3600.,'  /s'     &
                                         ,RemRateLMOPP ,'	/hr'
      WRITE(2,*)'RemRateLMOPSi	=	',RemRateLMOPSi/3600.,' /s'     &
                                         ,RemRateLMOPSi,'	/hr'
      WRITE(2,*)

      WRITE(2,*)'Q10_Rem        =       ',Q10_Rem,'             /C'
      WRITE(2,*)'TREF_Rem        =       ',TREF_Rem,'             C'

      WRITE(2,*)
      WRITE(2,*)

      WRITE(2,*)'Q10_Nit	=	',Q10_Nit,'			/C'
      WRITE(2,*)'NitriRate	=	',NitriRate/3600.,'   	/s'     &
                                         ,NitriRate,'		/hr'
      WRITE(2,*)'Gamma_M	=	',Gamma_M,'		-'
      WRITE(2,*)'TREF_Nit	=	',TREF_Nit,'		-'

       
      WRITE(2,*)'======================='
      WRITE(2,*)'Oxygen '
      WRITE(2,*)'======================='

      WRITE(2,*)'nuOXYGENNHs     =       ',nuOXYGENNHs,'        /d '
      WRITE(2,*)'nuNITNHs        =       ',nuNITNHs,'           /d'
      WRITE(2,*)'nuOXYGENODU     =       ',nuOXYGENODU,'         /d'
      WRITE(2,*)'nuNITODU        =       ',nuNITODU,'            /d'
      WRITE(2,*)'kNHsoxidSatOX   =       ',kNHsoxidSatOX,'    mmolO2/m3'
      WRITE(2,*)'kNHsoxidNOsInOX =       ',kNHsoxidNOsInOX,'  mmolO2/m3'
      WRITE(2,*)'kNHsoxidNOsSatNOs =     ',kNHsoxidNOsSatNOs,' mmolN/m3'
      WRITE(2,*)'kODUoxidSatOX   =       ',kODUoxidSatOX,'    mmolO2/m3'
      WRITE(2,*)'kODUoxidNOsInOX =       ',kODUoxidNOsInOX,'  mmolO2/m3'
      WRITE(2,*)'kODUoxidNOsSatNOs =     ',kODUoxidNOsSatNOs,' mmolN/m3'
      WRITE(2,*)'kOxicMinSatOX   =       ',kOxicMinSatOX,'   mmolO2/m3'
      WRITE(2,*)'kDenitInOX      =       ',kDenitInOX,'      mmolO2/m3'
      WRITE(2,*)'kDenitSatNOs    =       ',kDenitSatNOs,'    mmolN/m3'
      WRITE(2,*)'kAnoxMinInOX    =       ',kAnoxMinInOX,'    mmolO2/m3'
      WRITE(2,*)'kAnoxMinInNOs   =       ',kAnoxMinInNOs,'   mmolN/m3'
      WRITE(2,*)'kSolFormSatIron =       ',kSolFormSatIron,' mmolFE/m3'
      WRITE(2,*)'Rsolid          =       ',Rsolid
      WRITE(2,*)'gammaC_OX       =       ',gammaC_OX,'       mmolC/m3'
      WRITE(2,*)'gammaNHs_OX     =       ',gammaNHs_OX,'     mmolNH/m3'
      WRITE(2,*)'gammaODU_OX     =       ',gammaODU_OX,'     mmolODU/m3'
      WRITE(2,*)'gammaPOC_ODU    =       ',gammaPOC_ODU,'    mmolC/m3'
      WRITE(2,*)'gammaPOC_NOs    =       ',gammaPOC_NOs,'    mmolC/m3'
      WRITE(2,*)'gammaNHs_NOs    =       ',gammaNHs_NOs,'    mmolNH/m3'
      WRITE(2,*)'gammaODU_NOs    =       ',gammaODU_NOs,'    mmolODU/m3'
      WRITE(2,*)'ksatMortZoo     =       ',ksatMortZoo,'     mmolc/m3'
      WRITE(2,*)'AnoxMort        =       ',AnoxMort,'             /d'
      WRITE(2,*)'ksatOxicHydrol  =       ',ksatOxicHydrol,'  mmolO2/m3'


      WRITE(2,*)'=========='
      WRITE(2,*)' Bio '
      WRITE(2,*)'=========='

      WRITE(2,*)'Vitesses de chutes en m/s (convention de signe: <0) '
      WRITE(2,*)' '
      WRITE(2,*)'diac     =       ',wsed(1,idiac)*86400    ,'m/j'
      WRITE(2,*)'dian     =       ',wsed(1,idian)*86400    ,'m/j'
      WRITE(2,*)'diap     =       ',wsed(1,idiap)*86400    ,'m/j'
      WRITE(2,*)'diachl   =       ',wsed(1,idiachl)*86400    ,'m/j'
      WRITE(2,*)'diasi    =       ',wsed(1,idiasi)*86400    ,'m/j'
      WRITE(2,*)
      WRITE(2,*)'bactc    =       ',wsed(1,ibactc)*86400    ,'m/j'
      WRITE(2,*)
      WRITE(2,*)'smopc    =       ',wsed(1,ismopc)*86400    ,'m/j'
      WRITE(2,*)'smopn    =       ',wsed(1,ismopn)*86400    ,'m/j'
      WRITE(2,*)'smopp    =       ',wsed(1,ismopp)*86400    ,'m/j'
      WRITE(2,*)'smopchl  =       ',wsed(1,ismopchl)*86400    ,'m/j'
      WRITE(2,*)'smopsi   =       ',wsed(1,ismopsi)*86400    ,'m/j'
      WRITE(2,*)
      WRITE(2,*)'lmopc    =       ',wsed(1,ilmopc)*86400    ,'m/j'
      WRITE(2,*)'lmopn    =       ',wsed(1,ilmopn)*86400    ,'m/j'
      WRITE(2,*)'lmopp    =       ',wsed(1,ilmopp)*86400    ,'m/j'
      WRITE(2,*)'lmopsi   =       ',wsed(1,ilmopsi)*86400    ,'m/j'

 
      CLOSE(2)
      endif
 
!======================================================================
! 3. Writing of the pelagic parameters in an output file
! END:
!======================================================================

      RETURN
      END !SUBROUTINE Eco3mParameters 
