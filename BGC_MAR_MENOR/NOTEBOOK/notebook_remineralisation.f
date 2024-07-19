&notebook_remineralisation
!
! 3-D ecosystem model - Eco3Med: Remineralisation characteristics
!
! Implementation: Karline Soetaert - Caroline Ulses
!                 NIOO-CEME
!
! PARAMETERS
!
!###################################################################
!
! Remineralisation
 RemRateSMOPC  = 2.e-6       ! (unit s)     : Remineralisation rate for Small Particulate Organic Carbon      ! 5.e-7   1.25e-5
 RemRateSMOPN  = 2.4e-6      ! (unit s)     : Remineralisation rate for Small Particulate Organic Nitrogen    ! 6.e-7   1.52e-5
 RemRateSMOPP  = 2.8e-6      ! (unit s)     : Remineralisation rate for Small Particulate Organic Phosphate   ! 7.e-7   1.81e-5
 RemRateSMOPSi = 2.312e-7    ! (unit s)     : Remineralisation rate for Small Particulate Organic Silicium     
 RemRateSMOPChl= 4.64e-6     ! (unit s)     : Remineralisation rate for Small Particulate Organic Chlorophylle      

 RemRateLMOPC  = 2.e-6       ! (unit s)     : Remineralisation rate for Large Particulate Organic Carbon      
 RemRateLMOPN  = 2.4e-6      ! (unit s)     : Remineralisation rate for Large Particulate Organic Nitrogen     
 RemRateLMOPP  = 2.8e-6      ! (unit s)     : Remineralisation rate for Large Particulate Organic Phosphate     
 RemRateLMOPSi = 2.312e-7    ! (unit s)     : Remineralisation rate for Large Particulate Organic Silicium      

 Q10_Rem       = 2.95        ! (unit degC)  : Temperature coefficient for remineralisation
 TREF_Rem      = 20.0        ! (unit C      : Reference temperature for remineralisation

! Nitrification
 Q10_Nit       = 2.37        ! (unit degC)  : Temperature coefficient for nitrification dans fonction Fred
 NitriRate     = 5.91e-7     ! (unit s)     : Maximum Nitrification Rate at 0 dg C 
 Gamma_M       = 1.0         ! (unit -)     : To account for the effect of mixing on primary production
 TREF_Nit      = 10.0        ! (unit degC)  : Reference temperature for nitrification

/
