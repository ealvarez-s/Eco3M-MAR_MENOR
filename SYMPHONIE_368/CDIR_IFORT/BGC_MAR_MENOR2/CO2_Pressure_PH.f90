










!      subroutine CO2_Pressure_PH

! Written by Karline Soetaert - Jim Greenwood
!                 NIOO-CEME
! Implemented in S model by Caroline Ulses
!----------------------------------------------------------------!
! Based on alkalinity and DIC, calculates pH and pCO2            !
! using the equations in                                         ! 
! Millero, F.J. (1995), Geochim. Cosmochim. Acta 59(4): 661-677  !
! Archer, D. (1995).  Ecological applications 5(3): 724-739.     !
!                                                                ! 
! Alkalinity is defined as OH- + HCO3- + 2 CO3-- + BOH + si      !
!                        + HPO4-- + 2 PO4---                     !
!                        - H                                     !

! Based on this and the sumCO2 (DIC) the speciation is           !
! estimated based on the reactions:                              !
! CO2 + H2O <-> HCO3-+H+    dissociation ct K1_CO2               !
! HCO3-     <-> CO3-- +H+   dissociation ct K2_CO2               !

! Equilibrium expressed in terms of SumCO2 and H+:               !
! CO2   = H2    /(H2+K1_CO2H+K2_CO2K1_CO2)       *SumCO2         !
! HCO3- = K1_CO2H  /(H2+K1_CO2H+K2_CO2K1_CO2)    *SumCO2         !
! CO3-- = K2_CO2K1_CO2/(H2+K1_CO2H+K2_CO2K1_CO2) *SumCO2         !
! SumCO2 in mmol/m3                                              !

! Boric acid:
! H3BO3 + H2O <-> H+ + B(OH)4-           ct K_BOH3               !
! B(OH)4- = K_BOH3/(K_BOH3+H+) * SumBorate                       !

! silicic acid:                          ct K_Si                 !
! SiO2 + H2O <-> H+ + Si(OH)3O-                                  !
! Si(OH)3O-   = K_Si/(K_Si+H+) * SumSi

! Ammonia:                                                       !
! NH4+   <-> H+ + NH3                    ct K_NH4                !
! NH3     =  K_NH4 / (H+ + K_NH4) * SUMNHx                       !   
! phosphoric acid:                                               !
! H3PO4- <-> H+ + H2PO4-                 ct K1_H3PO4             !
! H2PO4- <-> H+ + HPO4--                 ct K2_H3PO4             !
! HPO4-- <-> H+ + PO4---                 ct K3_H3PO4             !
! H2PO4- = K1_H3PO4*H2                / Denom                    !
! HPO4-- = K2_H3PO4*K1_H3PO4*H        / Denom                    !
! PO4--- = K2_H3PO4*K1_H3PO4*K3_H3PO4 / Denom                    !
!         DENOM =(H3+H2*K1_H3PO4+H*K1_H3PO4K2_H3PO4+             !
!                               K1_H3PO4*K2_H3PO4*K3_H3PO4)*SUMP ! 
!----------------------------------------------------------------!
   MODULE ModuleEquilibrium

!----------------------------------------------------------------!
! The summed concentrations of species that affect pH            !
!----------------------------------------------------------------!
    DOUBLE PRECISION  :: SumCO2, SumP, SumSi, &
                         SumNHx, SumNit, SumBorate


    DOUBLE PRECISION  :: HPlus,SolvePHs,     &
                         TotalAlkalinity

   END MODULE ModuleEquilibrium

!************************************************************!

!   MODULE ModuleDissociation

!   DOUBLE PRECISION :: K_H2Os      ,& !  mol/l           Dissociation ct of water
!                       K1_CO2s     ,& !  mol/l           Dissociation ct of carbonic acid
!                       K2_CO2s     ,& !  mol/l           Dissociation ct of carbonic acid
!                       K_BOH3s     ,& !  mol/l           Dissociation ct of boric acid
!                       K1_H3PO4s   ,& !  mol/l           Dissociation ct of phosphorus
!                       K2_H3PO4s   ,& !  mol/l           
!                       K3_H3PO4s   ,& !  mol/l           
!                       K_HNO3s     ,& !  mol/l           
!                       K_NH4s      ,& !  mol/l           Dissociation ct of ammonium
!                       K_H2Ss      ,& ! mol/l           
!                       K_Aragonites,& !  mol/l           
!                       K_Calcites  ,& !  mol/l           
!                       K_Sis       ,& !  mol/l           Dissociation ct of silicic acid
!                       IonicStrength ! mol/kgLiq       

!   END MODULE ModuleDissociation

      subroutine CO2_Pressure_PH


!     use ModuleDeclaration
      USE ModuleComBlock
      use module_biology


      implicit none
      double precision :: HenryCst


      integer           :: I,II,NumIteration
      double precision  :: Resolution
      double precision :: Temp,   & ! Temperature, dgC
                                   S       ! Salinity, PSU

      double precision  :: T,P,R,LnK,LnKp,IMcal,DeltaV,DeltaK, &
                           lnK1fac,lnK2fac,lnKBfac,pK1,pK2

      double precision :: K_H2O      ,& !  mol/l           Dissociation ct of water
                       K1_CO2     ,& !  mol/l           Dissociation ct of carbonic acid
                       K2_CO2     ,& !  mol/l           Dissociation ct of carbonic acid
                       K_BOH3     ,& !  mol/l           Dissociation ct of boric acid
                       K1_H3PO4   ,& !  mol/l           Dissociation ct of phosphorus
                       K2_H3PO4   ,& !  mol/l           
                       K3_H3PO4   ,& !  mol/l           
                       K_HNO3     ,& !  mol/l           
                       K_NH4      ,& !  mol/l           Dissociation ct of ammonium
                       K_H2S      ,& ! mol/l           
                       K_Aragonite,& !  mol/l           
                       K_Calcite  ,& !  mol/l           
                       K_Si       ,& !  mol/l           Dissociation ct of silicic acid
                       IonicStrength, & ! mol/kgLiq
                       KS,KF,TS,TF,pKS,lnKF,HF ! mol/kg 

    double precision   :: SumCO2, SumP, SumSi, &
                         SumNHx, SumNit, SumBorate,HSO4


    double precision   :: HPlus,SolvePHs,     &
                         TotalAlkalinity

    integer            :: Iterations

    integer            :: MaxIterations 
    double precision   :: Residual, dResidual, dx,SWStoTOT

    double precision   :: Hini
    double precision   :: Tolerance 

    double precision   :: HOld, HNew
    double precision   :: dHplus 

    double precision   :: CoAlkalinity, NonCOAlkalinity, aa, bb, cc
    double precision   :: CalculateNonCOAlkalinity
    double precision   :: Denominator
    double precision   :: BorateAlk, Palk, SiAlk, WaterAlk,               &
                       NitrateAlk, AmmoniumAlk

    double precision   :: pHTol,pHx,pHGuess,ln10,deltapH,Denom,CAlk,Balk,OH, &
                          Phostop,PhosBot,Hfree,Slope, &
                          b,Delta,P1atm,FugFac,factor,FREEtoTOT  

    DOUBLE PRECISION, EXTERNAL    :: SolvePH

      meanIterations  = 0
      meanResolution  = 0.D0
      MaxIterations = 100
      Tolerance = 1.d-10 
      dHplus = 1d-12
      Hplus = 8.
      dx = 0.
      NumIteration = 1

      do I=1, NumPelagicBoxes

      ! estimating pH using Newton-Raphson method   
      ! Hplus is known from previous estimate, in mol/l 
      ! Remark that the other species are in µmol/l)

! The dissociation constants
!        CALL Thermo_2(Temperature(I),Salinity(I),Depth(I))

!-----------------------------------------------------------------
! Calculates the thermodynamic constants, 
! (dissociation cts, solubility products) in
! seawater at temperature (T, in deg C), salinity (S, in PSU), 
! and pressure (P, in atm) conditions.  
! The  calculated constants are on the SWS pH scale.
!
! K_H2O    = dissociation constant of H2O
! K1_CO2   = 1st dissociation constant of CO2(aq)
! K2_CO2   = 2nd dissociation constant of CO2(aq)
! K_BOH3   = dissociation constant of B(OH)3
! K1_H3PO4 = 1st dissociation constant of H2PO4
! K2_H3PO4 = 2nd dissociation constant of H3PO4
! K3_H3PO4 = 3nd dissociation constant of H3PO4
! K_HNO3   = dissociation constant of HNO3
! K_NH4    = dissociation constant of NH4
! K_H2S    = 1st dissociation constant of H2S
! 
! Also calculates ionic strength, of the solution (mol/kg-solution)
!
! SLIGHTLY DIFFERENT FROM Thermo 
!-----------------------------------------------------------------

      if(maskbio(I)==1) then

      ! Temperature in dg Kelvin
      T = Temperature(I) + 273.15
      temp = Temperature(I) 
      S = Salinity(I)

      ! P is applied pressure (in Bars) = (Total Pressure-1)
      P = 0.1*Depth(I)*1.01326

      R = 83.145 ! ml bar-1 K-1 mol-1

!==============================================!
! apparent thermodynamic EQUILIBRIUM CONSTANTS !
!==============================================!

! Dans Soetaert and Greenwood
!  Ionic strength (Lyman and Fleming, 1940)
!  (in units of mol kg-water).

!      IMcal            = 0.00147 + 0.019885*S + 0.000038*S*S
!      IonicStrength = IMcal*(1.0-S/1000.0)

! Dans CO2SYS
! This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
       IonicStrength  = 19.924 * S / (1000 - 1.005 * S)


!----------------
! Dissociation of carbonic acid constants 
       if(ichoixK1K2C.eq.1) then

! K1_CO2
! Roy et al : dans Raick
      lnK = 2.83655 - 2307.1266/T -1.5529413 * LOG(T)  &
           + (-0.20760841 - 4.0484/T) * SQRT(S)       &
          + 0.08468345 * S - 0.00654208 * (S**1.5)
      
      !Pressure Correction from Morse and Mackenzie 1990
      !DeltaV = 21.07 + (0.1850*temp-(0.002248*(temp**2)))
      !DeltaK = -0.00530 + (0.1310*temp- 1.160e-06*(temp**2))
     
      ! Pressure Correction from Millero 1995     
      ! Incl.typo corrections from Lewis & Wallace (CO2SYS)   
      
      DeltaV = -25.5 + (0.1271*temp)
      DeltaK = -3.08 + (0.877*temp)
      DeltaK = DeltaK/1000
      
      !lnKp = (-(DeltaV/(R*T))*P + (0.5*(DeltaK/(R*T)))*(P**2)) + lnK 
	  lnkp=lnk
      K1_CO2 = EXP(lnKp)

! K2_CO2
      lnK = -9.226508 - 3351.6106/T - 0.2005743* LOG(T)  &
            + (-0.106901773 - 23.9722/T) * SQRT(S)       &
            + 0.1130822 * S - 0.00846934 * (S**1.5)

      !Pressure Correction from Morse and Mackenzie 1990
      !DeltaV = -8.74 + 0.3*temp-0.004064*(temp**2)
      !DeltaK = -0.01104 +0.224*temp-3.056e-6*(temp**2) 

      ! Pressure Correction from Millero 1995 
      ! Incl.typo corrections from Lewis & Wallace (CO2SYS)     
      DeltaV = -15.82 - (0.0219*temp)
      DeltaK = 1.13 - (0.1475d-3*temp)
      DeltaK = DeltaK / 1000.

    !  lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
           lnkp=lnk

      K2_CO2 = EXP(lnKp)

      elseif(ichoixK1K2C.eq.2) then

! De CO2SYS.m
! MEHRBACH refit BY DICKSON AND MILLERO
! Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
! (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
! refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
! on the SWS pH scale in mol/kg-SW.
! Mehrbach et al gave results on the NBS scale.
! The 2s precision in pK1 is .011, or 2.6% in K1.
! The 2s precision in pK2 is .020, or 4.6% in K2.
! Valid for salinity 20-40.
! This is in Table 4 on p. 1739.

    pK1 = 3670.7/T - 62.008 + 9.7944 *log(T) &
             - 0.0118*S + 0.000116 * (S**2.)

    K1_CO2 = 10**(-pK1) ! this is on the SWS pH scale in mol/kg-SW


! This is in Table 4 on p. 1739.
    pK2 = 1394.7 /T + 4.777 - 0.0184 *S + 0.000118 * (S**2.)
    K2_CO2 = 10**(-pK2) ! this is on the SWS pH scale in mol/kg-SW

    

! PressureEffectsOnK1
    DeltaV  = -25.5 + 0.1271 * temp
    !                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
    DeltaK   = (-3.08 + 0.0877 * temp) /1000
    !                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; %
    !                 Millero,1979
    lnK1fac = (- DeltaV + 0.5 *DeltaK * P) * P / (R*T)

    K1_CO2  = K1_CO2 * exp(lnK1fac)

! PressureEffectsOnK2
    DeltaV  = -15.82 - 0.0219 * temp
    !                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
    DeltaK   = (1.13 - 0.1475 * temp) /1000
    !                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero,1979
    lnK2fac = (- DeltaV + 0.5 *DeltaK *P) *P / (R*T)

    K2_CO2  = K2_CO2 * exp(lnK2fac)

    endif 

!----------------
! Dissociation of water in seawater
! Idem dans Soetart et Greenwood et CO2SYS
! CO2SYS : "Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
! his check value of 1.6 umol/kg-SW should be 6.2"

      lnK = 148.9802 - 13847.26/T - 23.6521*LOG(T)      &
            +(-5.977+118.67/T+1.0495 * LOG(T)) * SQRT(S) &
            - 0.01615*S
      

      !Pressure Correction from Millero 1995
      !Typos from Lewis and Wallace (CO2SYS)

! Dans Soetart et Greenwood   pas de correction   
!      DeltaV = -25.6 + (0.2324*temp) + (-0.0036246*(temp**2))
!      DeltaK = -5.13 + (0.0794*temp)
!      DeltaK = DeltaK/1000
     ! lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
!	   lnkp=lnk

! Dans CO2SYS :
       DeltaV = -20.02 + 0.1119 *temp - 0.001409 * temp**2
       DeltaK = -5.13 + (0.0794*temp) /1000
       lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) 

      K_H2O = EXP(lnK) ! SWS pH scale

      K_H2O = K_H2O * exp(lnKp) ! SWS pH scale

!----------------
! Dissociation of ammonia in seawater
      lnK = -6285.33/T+0.0001635*T-0.25444      &
            +(0.46532-123.7184/T)* SQRT(S)      &
            +(-0.01992+3.17556/T)*S
      

      !Pressure Correction from Millero 1995
      !Typos corrections from Lewis and Wallace (CO2SYS)
      
      DeltaV = -26.43 + (0.0889*temp) + (-0.000905*(temp**2))
      DeltaK = -5.03 + (0.0814*temp)
      DeltaK = DeltaK/1000    
      
      !lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
	   lnkp=lnk

      K_NH4 = EXP(lnKp) 

! and of nitrate, a constant
      K_HNO3   = 23.44

!----------------
! Dissociation of boric acid

! Soetart and Greenwood
!      lnK = (-8966.9 - 2890.51*SQRT(S) - 77.942*S        &
!             +1.726*S**1.5 - 0.0993*S*S)/T               &
!             +(148.0248 + 137.194*SQRT(S) + 1.62247*S)   &
!             +(-24.4344 - 25.085 *SQRT(S) - 0.2474*S)*LOG(T) &
!             + 0.053105*SQRT(S)*T 
      !Pressure Correction
!      DeltaV = -29.48 + (-0.1622*temp) + (0.002608*(temp**2))
!      DeltaK = -2.84
!      DeltaK = DeltaK/1000
!!      lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
!      lnkp=lnk
!      K_BOH3 = EXP(lnKp)
! CO2SYS

! Riley, J. P., Deep-Sea Research 12:219-220, 1965:
! this is .000068.*Sali./35. = .00000195.*Sali
      TF = (0.000067 /18.998) *(S /1.80655) ! in mol/kg-SW

! Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
! this is .02824.*Sali./35. = .0008067.*Sali
      TS = (0.14 /96.062) *(S /1.80655) ! in mol/kg-SW

! Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
      pKS = 647.59 / T - 6.3451 + 0.019085 *T - 0.5208 *sqrt(IonicStrength)
      KS = 10**(-pKS) &          ! this is on the free pH scale in mol/kg-H2O
         * (1 - 0.001005*S)    ! convert to mol/kg-SW

! Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
      lnKF = 1590.2 /T - 12.641 + 1.525 *IonicStrength**0.5
      KF   = exp(lnKF) &                 ! this is on the free pH scale in mol/kg-H2O
             *(1 - 0.001005 *S)          ! convert to mol/kg-SW

      SWStoTOT  = (1 + TS /KS) / (1 + TS /KS + TF /KF)

      lnK = (-8966.9 - 2890.51*SQRT(S) - 77.942*S        &
             +1.728*S**1.5 - 0.0996*S*S)/T               &
             +(148.0248 + 137.1942*SQRT(S) + 1.62142*S)   &
             +(-24.4344 - 25.085 *SQRT(S) - 0.2474*S)*LOG(T) &
             + 0.053105*SQRT(S)*T

      K_BOH3 = EXP(lnK) &
               /SWStoTOT ! convert to SWS pH scale  

      !Pressure Correction
      DeltaV = -29.48 + (-0.1622*temp) + (0.002608*(temp**2))
      DeltaK = -2.84 
      DeltaK = DeltaK/1000    
      lnKBfac = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2))
      K_BOH3 = K_BOH3 * exp(lnKBfac) 

!----------------
! Dissociation of silicic acid

! Soetaert and Greenwood     	   
!	  K_Si = 4.E-10
!
!      !Pressure Correction - Assumed to be the same as Boric Acid (See Millero 1995)
!
!      
!      DeltaV = -29.48 + (-0.1622*temp) + (0.002608*(temp**2))
!      DeltaK = -2.84 
!      DeltaK = DeltaK/1000 
!      !Delta V and Delta K estimated from values for boric acid         
!      
!     ! lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + LOG(K_Si) 
!		lnkp=lnk
!      K_Si = EXP(lnKp) 

! Dans CO2SYS
       lnK = -8904.2 / T + 117.4 - 19.334 *LOG(T) +  & 
            (-458.79 / T + 3.5913) *sqrt(IonicStrength) + &
            (188.74 / T - 1.5998) *IonicStrength + &
            (-12.1652 / T + 0.07871) *IonicStrength**2

       K_Si = EXP(lnK) &
             *(1 - 0.001005 *S) ! convert to mol/kg-SW

! PressureEffectsOnKSi:
!  "The only mention of this is Millero, 1995 where it is stated that the
!   values have been estimated from the values of boric acid. HOWEVER,
!    there is no listing of the values in the table.
!    I used the values for boric acid from above."
       DeltaV = -29.48 + (-0.1622*temp) + (0.002608*(temp**2))
       DeltaK = -2.84 
       DeltaK = DeltaK/1000
       lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2))

       K_Si = K_Si * exp(lnKp)

!----------------
! Dissociation of phosphoric acid
! K1_H3PO4 is included for completeness
       
      lnK = 115.54-4576.752/T-18.453*LOG(T)+(0.69171-106.736/T)*(S**0.5)     &
              +(-0.01844-0.65643/T)*S

     !Pressure correction for K1_H3PO4      
 
      DeltaV = -14.51 + (0.1211*temp) + (-0.000321*(temp**2))
      DeltaK = -2.67 + (0.0427*temp)
      DeltaK = DeltaK/1000    
      
! Dans Soetart la correction de pression est commentée       
!	  !lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
!        lnkp=lnk
!      K1_H3PO4 = EXP(lnKp)

! Dans CO2SYS
      K1_H3PO4 = EXP(lnK)
      lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2))  
      K1_H3PO4 = K1_H3PO4 *  exp(lnKp)                         

! K2_H3PO4
      
      lnK = 172.1033-8814.715/T-27.927*LOG(T)+(1.3566-160.340/T)*(S**0.5)    &
	        +(-0.05778+0.37335/T)*S
	   
      !Pressure Correction 
      
      DeltaV = -23.12 + (0.1758*temp) + (-0.002647*(temp**2))
      DeltaK = -5.15 + (0.09*temp)
      DeltaK = DeltaK/1000    

! Dans Soetart la correction de pression est commentée      
!      !lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
!		lnkp=lnk
!      K2_H3PO4 = EXP(lnKp) 

! Dans CO2SYS
      K2_H3PO4 = EXP(lnK)
      lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2))
      K2_H3PO4 = K2_H3PO4 * exp(lnKp)

! K3_H3PO4

      lnK = -18.126-3070.75/T+(2.81197+17.27039/T)*(S**0.5)     &
	        +(-0.09984-44.99486/T)*S

     !Pressure Correction
      DeltaV = -26.57 + (0.2020*temp) + (-0.003042*(temp**2))
      DeltaK = -4.08 + (0.0714*temp)
      DeltaK = DeltaK/1000    

! Dans Soetart la correction de pression est commentée      
!     ! lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 
!        lnkp=lnk
!      K3_H3PO4 = EXP(lnKp) 

! Dans CO2SYS
       K3_H3PO4 = EXP(lnK)
       lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2))
       K3_H3PO4 = K3_H3PO4 * exp(lnKp)      

!Dissociation of hydrogen sulfide

      lnK = 225.838-13275.3/T-34.6435*LOG(T)+0.3449*(S**0.5)-0.0274*S 

      !Pressure correction

      DeltaV = -14.8 + (0.0020*temp) + (-0.0004*(temp**2))
      DeltaK = 2.89 + (0.054*temp)
      DeltaK = DeltaK/1000    
      
      lnKp = (-(DeltaV/(R*T))*P + (0.5*DeltaK/(R*T))*(P**2)) + lnK 

      K_H2S = EXP(lnKp)   


! End of calculation of dissociation constants

! Total concentration of species contributing to alkalinity 
! Pour CO2SYS je convertis en mol/kg
!    SumBorate       = 4.106e2 * Salinity(I) / 35.    ! µmol/l ou mmol/m3
! Uppstrom, L., Deep-Sea Research 21:161-162, 1974
      SumBorate       = 4.157e2 * Salinity(I) / 35.   *1000/DensEco(I)/1e6 ! µmol/l ou mmol/kg
      SumNHx          = Ammonium          (I) *1000/DensEco(I)/1e6
      SumNit          = Nitrate           (I) *1000/DensEco(I) /1e6
      SumCO2          = DIC               (I) *1000/DensEco(I)/1e6
      SumP            = Phosphate         (I) *1000/DensEco(I)/1e6
      SumSi           = Silice            (I) *1000/DensEco(I)/1e6
      TotalAlkalinity = Alkalinity        (I) *1000/DensEco(I)/1e6
!      Hplus           = 10**          (-pH(I))
!      Hplus           = 10**          (-8.) ! First guess pris dans CO2SYS


!        pH(I) = EstimatepH(NumIteration,Resolution)
! Estimatation of PH
!        Hini       = HPlus
         Resolution = 1.D0
!        numIteration = MaxIterations

!Newton-Raphson method
! CO2SYS mthod
        pHGuess     = 8
        pHTol       = 0.0001 ! tolerance for iterations end
        ln10        = log(10.)
        pHx         =pHGuess
        deltapH     = pHTol+1

        do while (abs(deltapH) > pHTol)
        HPlus        = 10**(-pHx)
        Denom     = (HPlus*HPlus + K1_CO2*HPlus + K1_CO2*K2_CO2)
        CAlk      = SumCO2 *K1_CO2*(HPlus + 2 *K2_CO2) /Denom
        BAlk      = SumBorate *K_BOH3 /(K_BOH3 + HPlus)
        OH        = K_H2O /HPlus
!       PhosTop   = K1_H3PO4 *K2_H3PO4 *HPlus + 2.*K1_H3PO4 *K2_H3PO4 * K3_H3PO4 - HPlus*HPlus*HPlus ! CO2Sys
        PhosTop   = K1_H3PO4 *HPlus *HPlus + 2* K1_H3PO4*K2_H3PO4*HPlus + 3 *K1_H3PO4 *K2_H3PO4 * K3_H3PO4 ! Karline
        PhosBot   = HPlus*HPlus*HPlus + K1_H3PO4*HPlus*HPlus + K1_H3PO4*K2_H3PO4*HPlus + K1_H3PO4*K2_H3PO4*K3_H3PO4
        PAlk      = SumP * PhosTop / PhosBot
        SiAlk     = SumSi *K_Si /(K_Si + HPlus)
        AmmoniumAlk = SumNHx - (K_NH4 /(K_NH4+Hplus)* SumNHx)
        ! test CO2sys -- ajouté par Kim !!!  AmmoniumAlk = AmmoniumAlk+3
        NitrateAlk  =                                 SumNit
        FREEtoTOT = (1 + TS /KS) ! pH scale conversion factor
       Hfree     = HPlus/FREEtoTOT ! for H on the total scale
       HSO4      = TS /(1 + KS /Hfree) ! since KS is on the free scale
       HF  =          TF /(1 + KF /Hfree) ! since KF is on the free scale
    Residual  = TotalAlkalinity -  CAlk - ( BAlk + OH + PAlk + SiAlk + NitrateAlk - AmmoniumAlk - HPlus )
    ! find Slope dTA/dpH
    ! (this is not exact, but keeps all important terms)
    Slope     = ln10 *(SumCO2 *K1_CO2*HPlus*(HPlus*HPlus + K1_CO2*K2_CO2 + 4*HPlus*K2_CO2) /Denom /Denom+ BAlk *HPlus/(K_BOH3 + HPlus) + OH + HPlus)
    deltapH   = Residual /Slope ! this is Newton's method
    ! to keep the jump from being too big;
    do while (abs(deltapH) > 1)
        deltapH=deltapH /2
    enddo
    pHx       = pHx + deltapH ! Is on the same scale as K1 and K2 were calculated...
enddo


! Utile????
      KNHs(I)   = K_NH4
      K1H3PO4(I)= K1_H3PO4
      K2H3PO4(I)= K2_H3PO4
      K3H3PO4(I)= K3_H3PO4
!      K1CO2(I) = K1_CO2
!      K2CO2(I) = K2_CO2

        meanIterations  = MeanIterations + NumIteration
        meanResolution  = MeanResolution + Resolution

!Calculate Carbonate speciation

! Solubility constant
      ! Temperature in dg Kelvin
      T = Temperature(I) + 273.15
      S = Salinity (I)

      HenryCst = exp( -60.2409 + 9345.17/T                 &
            +23.3585* DLOG(T/100.) +                     &
            S*( 0.023517-2.3656d-4*T+4.7036d-7*T*T  ))

! Calculate CO2 species

!      PAlkCO2sys= SumP-2.d0*(Hplus*Hplus*K1_H3PO4 + K1_H3PO4*K2_H3PO4*Hplus + K1_H3PO4*K2_H3PO4*K3_H3PO4)/PhosBot*SumP
      HPlus = 10**(-pHx)
      Denominator  = Hplus*(K1_CO2+Hplus) + K1_CO2*K2_CO2
! je convertis en mmol/m3
      CO2(I)          = Hplus *HPlus /denominator * SumCO2 / 1000 * DensEco(I)*1.e6
      HCO3(I)         = Hplus *K1_CO2/denominator * SumCO2 / 1000 * DensEco(I)*1.e6
      CO3(I)          = K2_CO2*K1_CO2/denominator * SumCO2 / 1000 * DensEco(I)*1.e6
      HPlus = 10**(-pHx)*1000 / DensEco(I) ! mmol/m3
      pH(I)           = -LOG10(10**(-pHx)) ! idem que K1 et K2

! directement dans graph_out.F90
     factor = -log(SWStoTOT) /log(0.1)
      pHT(I)          =  pH(I)    - factor ! total scale => directement dans
!     TA(I) =    (TotalAlkalinity - SumNit + SumNHx - SumP )*1.e6 ! µmol/kg 

! Vérif termes supplémentaires ou non compris dans TACO2SYS
!      if(iipoint.eq.557.and.jjpoint.eq.486) then
!      print*,'CO2_Pressure_PH Alk',TotalAlkalinity*1.e6, &
!           (TotalAlkalinity - SumNit + SumNHx - SumP )*1.e6
!    print*,'CO2_Pressure_PH Terme',2*TS,SumSi * K_Si / (K_Si+HPlus)*1.e6,HSO4*1.e6,HF*1.e6
!    print*,'CO2_Pressure_PH pH',pH(I),pHT(I),I,Depth(I)
!    endif
 


     ! partial CO2 pressure is calculated based on the solubility equation
     ! Convert from atm to µatm and from mmol/m3 to mol/kg
     ! K0, moles l/atm or µmoles/l/µatm
         !KO is in mol/kg/atm or µmol/kg/µatm

          pCO2(I) = CO2(I)*1000/DensEco(I)/HenryCst

          Delta = (57.7 - 0.118 *T)
          b = -1636.75 + 12.0408 *T - 0.0327957 *T**2 + 3.16528 *0.00001 *T**3
! For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
         P1atm = 1.01325 ! in bar
         FugFac = exp((b + 2 *Delta) *P1atm / (R*T))

         pCO2(I) = pCO2(I) / FugFac 

        K0_CO2(I)=HenryCst

!         print*,'co2_pressure',I,CO2(I),SumCO2,pCO2(I),pH(I)

   endif ! mask

    ENDDO
    meanIterations =      meanIterations /NumPelagicBoxes
    meanResolution =      meanResolution /NumPelagicBoxes


   RETURN
   END SUBROUTINE CO2_Pressure_PH

!*************************************************************!
! CALCULCATE THE PH BY SOLVING A QUADRATIC EQUATION
!*************************************************************!

   DOUBLE PRECISION FUNCTION SolvePH(K_H2Os,HPluss,K_NH4s,SumNHxs,SumNits,K_BOH3s, &
                                     SumBorates,K1_H3PO4s,K2_H3PO4s,K3_H3PO4s,SumPs, &
                                     K_Sis,SumSis,TotalAlkalinitys,SumCO2s,K1_CO2s,K2_CO2s)     

   USE MODULECOMBLOCK
   USE ModuleEquilibrium
!  USE ModuleDissociation
   IMPLICIT NONE
   DOUBLE PRECISION  :: CoAlkalinity, NonCOAlkalinity, aa, bb, cc
   DOUBLE PRECISION, EXTERNAL :: CalculateNonCOAlkalinity
   DOUBLE PRECISION :: WaterAlk,AmmoniumAlk,NitrateAlk,BorateAlk,PAlk,SiAlk
   DOUBLE PRECISION :: Denominator,Hpluss,TotalAlkalinitys
   DOUBLE PRECISION :: K_H2Os      ,& !  mol/l           Dissociation ct of water
                       K1_CO2s     ,& !  mol/l           Dissociation ct of carbonic acid
                       K2_CO2s     ,& !  mol/l           Dissociation ct of carbonic acid
                       K_BOH3s     ,& !  mol/l           Dissociation ct of boric acid
                       K1_H3PO4s   ,& !  mol/l           Dissociation ct of phosphorus
                       K2_H3PO4s   ,& !  mol/l
                       K3_H3PO4s   ,& !  mol/l
                       K_HNO3s     ,& !  mol/l
                       K_NH4s      ,& !  mol/l           Dissociation ct of ammonium
                       K_H2Ss      ,& ! mol/l
                       K_Aragonites,& !  mol/l
                       K_Calcites  ,& !  mol/l
                       K_Sis       ,& !  mol/l           Dissociation ct of silicic acid
                       SumNHxs     ,&
                       SumBorates  ,&
                       SumPs       ,&
                       SumSis      ,&
                       SumNits     ,&
                       SumCO2s


! 1. Estimate the contribution of the less important species directly

!        print*,'water',K_H2Os,HPluss
      ! Water (assume that water concentration is 1mol; convert to µmol)
        WaterAlk    = K_H2Os/HPluss*1.e6

      ! Ammonium
        AmmoniumAlk = SumNHxs - (K_NH4s /(K_NH4s+Hpluss)* SumNHxs)
        ! test CO2sys -- ajouté par Kim !!!  AmmoniumAlk = AmmoniumAlk+3
      ! Nitrate 
        NitrateAlk  =                                 SumNits
      ! It is currently assumed that HNO3 is fully dissociated 

      ! Borate
!        print*,'borate',K_BOH3s,Hpluss,SumBorates
        BorateAlk   = K_BOH3s / (K_BOH3s+Hpluss)       * SumBorates

      ! Phosphate
        Denominator = (Hpluss**3)+(HPluss*HPluss*K1_H3PO4s)+                    &
                              (HPluss*K1_H3PO4s*K2_H3PO4s)+(K1_H3PO4s*K2_H3PO4s*K3_H3PO4s)
      !                H2PO4-      + 2*HPO4--    + 3*PO4---
      PAlk        = (Hpluss*HPluss*K1_H3PO4s + 2*K1_H3PO4s*K2_H3PO4s*HPluss + 3*K1_H3PO4s*K2_H3PO4s*K3_H3PO4s)  &
                   /Denominator                  * SumPs
                
      ! Silicate
        SiAlk       = K_Sis / (K_Sis+Hpluss)           * SumSis


      ! Alkalinity component not due to CO species 
      !(Hplus is mol/l, others µmol/l)
         NonCOAlkalinity = SiAlk + PAlk                  &
                                  + BorateAlk                     &
                                  + WaterAlk  + NitrateAlk        &
                                  - HPluss*1e6 - AmmoniumAlk


      ! 2. Estimate H+ by solving the equation for the DIC species

      ! 3. carbonate alkalinity
        COAlkalinity = TotalAlkalinitys - NonCOAlkalinity

        !COalkalinity equals HCO3-  +2* CO3--

        ! HCO3- = K1_CO2H  /(H2+K1_CO2H+K2_CO2K1_CO2) *SumCO2
        ! CO3-- = K2_CO2K1_CO2/(H2+K1_CO2H+K2_CO2K1_CO2) *SumCO2
        ! This gives rise to a quadratic equation in H

        aa = COAlkalinity
        bb = (COAlkalinity-1.d0*SumCO2s)*K1_CO2s
        cc = (COAlkalinity - 2.d0*SumCO2s)*K1_CO2s*K2_CO2s

! new guess of H+
       SolvePH = (-bb + SQRT(BB*BB-4.D0*aa*cc) )/(2.D0*aa)

! ps: SolvePH and HPlus should converge


          RETURN
   END FUNCTION SOLVEPH
