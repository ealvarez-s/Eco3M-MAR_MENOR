&notebook_zooplk
!
! 3-D ecosystem model - Eco3Med: Zooplankton characteristics 
!
! Implementation: Karline Soetaert - Caroline Ulses
!                 NIOO-CEME
!
! PARAMETERS
!
!###################################################################
!
! Grazing
6.75e-5  /s              grazNano        Maximal grazing rate of NanoZooplancton at 15.5 dg C
0.65            -       prefNanoSyne   Preference factor of NanoZooplancton for Synechococcus
0.0             -       prefNanoNano    Preference factor of NanoZooplancton for Nanophytoplancton
0.0        -prefNanoZooNano Preference factor of NanoZooplancton for NanoZooplancton
5.            mmolC/m3 kZooNano Semi-saturation constant for grazing

4.5e-5  /s              grazMicro          Maximal grazing rate of MicroZooplancton at 15.5 dg C 
0.08   -      prefMicroBact    Preference factor of MicroZooplancton for Bacteria
0.06            -      prefMicroSyne    Preference factor of MicroZooplancton for Synechococcus
0.3     !.25    -      prefMicroNano    Preference factor of MicroZooplancton for NanoPhytoplancton
0.15      !.2  -      prefMicroDia     Preference factor of MicroZooplancton for Diatoms
0.12            -      prefMicroZooMicro Preference factor of MicroZooplancton for MicroZooplancton
0.25            -      prefMicroZooNano Preference factor of MicroZooplancton for NanoZooplancton
8.5            mmolC/m3         kZooMicro       Semi-saturation constant for grazing

2.25e-5  /s             grazMeso            Maximal grazing rate of MesoZooplancton at 15.5 dg C 
0.5           -      prefMesoDia    Preference factor of MesoZooplancton for Diatoms
0.45           -      pefMesoZooMicro    Preference factor of MesoZooplancton for MicroZooplancton
0.0            -      prefMesoLMOP Preference factor of MesoZooplancton for Large Particulate Organic Matter
20.            mmolC/m3  kZooMeso                Semi-saturation constant for grazing

0.8 -kcNetGrowthEff Net Growth Efficiency

! Messy feeding
0.23 - MessyFeedingfracC Fraction of messy feeding	
0.23            -       MessyFeedingfracN       Fraction of messy feeding
0.23            -       MessyFeedingfracP       Fraction of messy feeding

! Egestion
0.6             -       NanoAssEffBact  Assimilation efficiency of NanoZooplankton for Bacteria
0.6             -       NanoAssEffSyne   Assimilation efficiency of NanoZooplankton for Synechococcus 
0.6             -       NanoAssEffNano  Assimilation efficiency of NanoZooplankton for NanoPhytoplankton
0.6             -  NanoAssEffZooNano Assimilation efficiency of NanoZooplankton for NanoZooplankton   

0.6             -       MicroAssEffBact         Assimilation efficiency of MicroZooplankton for Bacteria
0.6             -       MicroAssEffSyne         Assimilation efficiency of MicroZooplankton for Synechococcus
0.6             -       MicroAssEffNano         Assimilation efficiency of MicroZooplankton for NanoPhytoplankton
0.6             -       MicroAssEffDia       Assimilation efficiency of MicroZooplankton for Diatoms
0.6             -       MicroAssEffZooNano      Assimilation efficiency of MicroZooplankton for NanoZooplankton
0.6             -       MicroAssEffZooMicro     Assimilation efficiency of MicroZooplankton for MicroZooplankton
0.6             -       MicroAssEffSMOP         Assimilation efficiency of MicroZooplankton for Small Particulate Organic Matter

0.6             -       MesoAssEffDia          Assimilation efficiency of MesoZooplankton for Diatoms
0.6            -       MesoAssEffZooMicro      Assimilation efficiency of MesoZooplankton for MicroZooplankton
0.6             -       MesoAssEffSMOP          Assimilation efficiency of MesoZooplankton for Small Particulate Organic Matter
0.6             -       MesoAssEffLMOP       Assimilation efficiency of MesoZooplankton for Large Particulate Organic Matter

0.8            - DiaPropMOPSi Ratio small/large organic matter in residu of egestion

! Excretion  
0.18            mmolN/mmolC NCZooNano Fixed N/C Ratio for NanoZooplankton 
0.013           mmolP/mmolC PCZooNano Fixed P/C Ratio for NanoZooplankton

0.18            mmolN/mmolC     NCZooMicro      Fixed N/C Ratio for MicroZooplankton
0.013           mmolP/mmolC     PCZooMicro      Fixed P/C Ratio for MicroZooplankton

0.18            mmolN/mmolC      NCZooMeso       Fixed N/C Ratio for MesoZooplankton
0.013           mmolP/mmolC      PCZooMeso       Fixed P/C Ratio for MesoZooplankton

! Mortality  
1.044e-6 !2.e-6    !1.044e-6     ! -45%   2.e-6 m3/s/mmolC      NanoMortRate   NanoZooplankton mortality rate
3.8e-7 !2.e-6   !3.8e-7     ! -80%  2.e-6    m3/s/mmolC      MicroMortRate   MicroZooplankton mortality rate
6e-6 !6.045e-7  !3.4722e-06 !6.045e-7	   ! -45%   9.e-7	m3/s/mmolC	MesoPredRate	Predation mortality rate for MesoZooplankton 

1.            	-		ZooPropMOPC	Ratio small/large organic matter in cadaver of zooplankton C
1.            	-		ZooPropMOPN	Ratio small/large organic matter in cadaver of zooplankton N
1.            	-		ZooPropMOPP	Ratio small/large organic matter in cadaver of zooplankton P  

0.90            	-               ZooPropMesoMOPC     Ratio small/large organic matter in cadaver of zooplankton C     ! .95
0.90            	-               ZooPropMesoMOPN     Ratio small/large organic matter in cadaver of zooplankton N     ! .95
0.90            	-               ZooPropMesoMOPP     Ratio small/large organic matter in cadaver of zooplankton P     ! .95
1.1		mmolC/m3	Kclos			Half-saturation for closure		

2.		/C		Q10_Z		Temperature coefficient for Zooplankton processes
18.	  !20	C		TREF_Z		Reference temperature for Zooplankton processes	
