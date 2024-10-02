!######################################################################
!
! 1-D ecosystem model
!
! Implementation: Karline Soetaert - Caroline Ulses
!                 NIOO-CEME
!
! Light profile in water column
!
! Calculates:
!    Photosynthetically active radiation in each box, W/m2     
!    HeatInput to each box, W/m2
!
! Calls:
!    CalculateChlorophyll(I)
! 
! Uses:
!    General, Pelagic, Radiation
!
! Functions/subroutines: none
!
! Revised: 15-10-2004
!
!######################################################################      


             !<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>!
             !<<<                                            >>>!
             !<<<        MAIN SUBROUTINES AND MODULES        >>>!
             !<<<                                            >>>!
             !<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>!


!********************************************************************
!  INITIALISATION
!********************************************************************

      SUBROUTINE InitialiseLight
      USE ModuleComBlock
      USE UserDeclaration
      IMPLICIT NONE

! Background attenuation coefficient for shorter wavelength PAR               
! a function of water depth

        HeatInput(:) = 0.D0
        PAR_Z(:)       = 0.D0
   
      END SUBROUTINE InitialiseLight

!********************************************************************

      SUBROUTINE InitialiseLight2
      END SUBROUTINE InitialiseLight2


!********************************************************************
!  DYNAMICS
!********************************************************************

   !SUBROUTINE DynamicsLight(IEupho)
      SUBROUTINE DynamicsLight(IEupho,Iptloc_,Jptloc_)
      use module_parallele
!----------------------------------------------------------------------*
! Calculates the depth gradient of photosynthetically active radiation * 
! (PAR) as a function of total radiation, the background extinction    *
! coefficient of long- and shortwave PAR and chlorophyll               *
! Also calculates the heat input to each box.                          *
!----------------------------------------------------------------------*

      USE MODULECOMBLOCK
!     USE MODULEGENERAL
      IMPLICIT NONE

      DOUBLE PRECISION ::                   &                               
            HeatRad,                        & ! W/m2   
            RADtop,                         & ! W/m2   Radiation entering water column
            ParTop,PARtopLong,ParTopShort,  & ! W/m2   PAR at top of each water column box
            ParBot,PARbotLong,ParBotShort,  & ! W/m2   PAR at base of each water column box
            ExtincLong,ExtincShort,         & ! /m     Light extinction coeffcient
            ExtincCDOM                        ! /m     Light extinction due to CDOM, function of MODC

      DOUBLE PRECISION, EXTERNAL :: sineWave
      DOUBLE PRECISION, EXTERNAL :: CalculateChlorophyll
      DOUBLE PRECISION           :: Parlong,ParShort
      INTEGER                    :: I,IEupho
      INTEGER                    :: Iptloc_,Jptloc_
      INTEGER                    :: EuphoCalculation

!-------------------------------------------------------------------
! PAR at the sea surface (in W/m2); some radiation is reflected, some is heat


      EuphoCalculation = 1 

      IF (SolarRadiation == 0.D0) THEN

        Par_z(:)       = 0.D0
        HeatInput(:) = 0.D0     

        RETURN

      ENDIF

      RADtop  = SolarRadiation  * (1.d0-AlbedoBio) 

      PARtop  = RADtop * (1.D0-pHeat)
      HeatRad = RADtop - ParTop

      ! Part of the PAR belongs to the long wavelength (ParTopLong) 
      ! and has extinction coeff. kdLong
      ! The remainder (ParTopshort) has extinction coeff. kdShort

      ParTopLong  = pLong*Partop        
      ParTopshort = Partop - ParTopLong

      ! Loop over all water column boxes
      DO I = 1, NumPelagicBoxes

        ! The extinction coefficient (1/m) in water column box I
        ! a background extinction and a chlorophyll-specific term

        !Ref: Para LMGEM : dans le modele on n'a que le labile: coherent?
        ExtincCDOM=max(0.D0,0.0068*MODC(I)-0.4579)

        DiaChl(I)        =max(0.D0,DiaChl(I))
        NanoChl(I)       =max(0.D0,NanoChl(I))
        SyneChl(I)       =max(0.D0,SyneChl(I))
        MIP(I)           =max(0.D0,MIP(I))

        ExtincLong  = kro + krp*(DiaChl(I)+NanoChl(I)+SyneChl(I))**lr
        ExtincShort = kbo + kbp*(DiaChl(I)+NanoChl(I)+SyneChl(I))**lb   
#ifdef mes
        ExtincShort = ExtincShort + kbn*(MIP(I)*1000.) &  !Babin et al 2003, la MIP est en g/L, dans formule g/m3 
                                  + ExtincCDOM 
#endif


       ! IF vertical: Par at bottom of box; 
        IF (PelagicVertical) THEN

          IF(ExtincLong  *Thickness(I)<50.and.ParTopLong>1.d-50) THEN
            ParBotLong = ParTopLong  * EXP(-ExtincLong  *Thickness(I))
          ELSE
            ParBotLong = 0.D0
          ENDIF

          IF(ExtincShort *Thickness(I)<50.and.ParTopShort>1.d-50) THEN
            ParBotShort= ParTopShort * EXP(-ExtincShort *Thickness(I))
          ELSE
            ParBotShort = 0.D0
          ENDIF

          ParBot = ParBotShort + ParBotLong

        ENDIF

           
       ! PAR in the middle of the box
        IF(ExtincLong  *Thickness(I)<50.and.ParTopLong>1.d-50) THEN
          ParLong = PARtopLong  * EXP(-ExtincLong  * Thickness(I)*0.5D0)  
        ELSE
          ParLong = 0.D0
        ENDIF
        IF(ExtincShort  *Thickness(I)<50.and.ParTopShort>1.d-50) THEN
          ParShort = PARtopShort * EXP(-ExtincShort * Thickness(I)*0.5D0)
        ELSE
          ParShort = 0.D0
        ENDIF

        PAR_Z(I) =  ParLong  + ParShort 

          ! Calculation of the Euphotic Layer Depth  
        !  IF (PelagicVertical) THEN
          IF(EuphoCalculation.EQ.1.AND.ParBot.LT.0.01*PAR_Z(1)) THEN

           ParTop = ParTopLong + ParTopShort

           EuphoticLayerDepth =  DepthInt(I) + Thickness(I)*log( 0.01*PAR_Z(1) / ParTop )                & 
                                                           /log( ParBot      / ParTop )
           !EuphoCalculation = 0 

           IEupho=I !Indice du bas de la couche euphotique
          ENDIF      
        !  ENDIF !PelagicVertical


         ! PAR also decreases due to calcite scattering
         ! Uses photon budget from COPT model of Tyrell et al (1999)
         ! An extra 1.6% of incident light is scattered for every 8mmol(CaCO3-C)/m3 
         !      PAR_Z(I)  = PAR_Z(I) - (PAR_Z(I) * (0.016/8) * CalculateCalcite(I))

        IF (I == 1) THEN
          ParSurfaceLong   = ParTopLong
          ParSurfaceShort  = ParTopShort
          kdSurfaceLong    = ExtincLong
          kdSurfaceShort   = ExtincShort
        ENDIF
     

        IF (PelagicVertical) THEN
         ! Absorbed radiation is heatinput
          heatInput(I) =   ParTopLong + ParTopShort                       &
                         - ParBotLong - ParBotShort
         ! Par at bottom of this box = par at top of next box
          ParTopLong  = ParBotLong
          ParTopShort = ParBotShort
        ELSE
          HeatInput(I) = HeatRad
        ENDIF 

      ENDDO

     ! All heat radiation is absorbed within the uppermost depth increment
     ! Ivanoff, 1977
      IF (PelagicVertical) THEN
 
        HeatInput(1) = HeatInput(1) + HeatRad

       ! All heat that remains in the bottom compartment is absorbed; 
       ! (to close the heat budget)
        HeatInput(NumPelagicBoxes) = HeatInput(NumPelagicBoxes)            &
                                   + ParbotLong +ParBotShort
 
        PARbottom = ParBotLong + ParBotShort
      ENDIF

      END SUBROUTINE DynamicsLight


!**********************************************************************
! END OF SIMULATION
!**********************************************************************
 
      SUBROUTINE FinaliseLight
      END SUBROUTINE FinaliseLight
