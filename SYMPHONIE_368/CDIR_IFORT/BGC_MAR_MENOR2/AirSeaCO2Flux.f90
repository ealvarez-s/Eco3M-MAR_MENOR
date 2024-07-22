










      subroutine AirSeaCO2Flux
      use module_principal
      use module_biology
      use ModuleDeclaration
      use module_parallele
      implicit none
      double precision ::Temp,Wind,kbul,Sal,   &
                       Rho0,     &  ! density of water, kg/m3 
                       pCO2Air,  &  ! partial pressure of CO2 in air (atm)     
                       pCO2Wat,  &  ! partial pressure of CO2 in water (atm
                       piston       ! piston velocity (m/d))
      double precision :: T,TT,TTT
      double precision :: TransferRate, HenryCst,lnK, &
                                 Schmidt    &  ! The Schmidt number
           ,B_A1,B_A2,B_A3,B_B1,B_B2,B_B3 &
          ,Temperature,Salinity,Ts,Kw,Pression,PressionVap,density_l &
          ,TemperatureK,TemperatureAir,Xi,zbub,BetaX,DeltaE,Kbw,vt,dm 
      double precision :: Cdg,Uastar,Uwstar,Fs,Fp,Fc,Kpg,Kc,rho_a,XCO2, &
                         deltaP,betag,bunsen,lnbunsen,diff,Md,Mv,visc, &
                         Ft
 

!      double precision :: pi,Jour

!       MeanpCO2Air        = 391.5      !360  µatm     Atmospheric CO2 partial pressure
!       AmplitudepCO2Air   = 0.001     !0.01          -
!       PeriodpCO2Air      = 365    ! hour     Periodicity in pCO2Air variation
!       PhasepCO2Air       = 1            !-1.432        radans   Phase
      kbul = 1.d0 ! tuning coefficient of the aire sea exchanges    

!----------------------------------------------------------------!
! Calculates the exchange of CO2 across the air-sea interface    !
! assuming a certain piston velocity (m/d), the partial pressure !
! of CO2 in the air and surface water, and the solubility of     !
! CO2 in seawater                                                !
!----------------------------------------------------------------!

! Update of the pCO2Air
      call atm_bio_upd

      TPS_ASF_2D=TPS_ASF_2D+1

! Calculate the exchanges
      do i=1,imax
      do j=1,jmax

      if(mask_t(i,j,kmax)==1) then

      Temp=tem_t(i,j,kmax,1)
      Sal =sal_t(i,j,kmax,1) 
      Rho0=rhp_t(i,j,kmax)+rho !1000 a verifier
      pCO2Wat = pCO2W(i,j)
!      pCO2Air = atm_bio(1)  

      TT= Temp*Temp
      TTT = TT*Temp

      Temperature= max(tem_t(i,j,kmax,1),0.)
      TemperatureK=Temperature+273.15
      Salinity=max(sal_t(i,j,kmax,1),0.)
      Pression=pss_w(i,j,1) ! en Pa
      PressionVap=q2_t(i,j,1)*pss_w(i,j,1)/0.622 ! en Pa
      TemperatureAir=teta2_t(i,j,1)
      Density_l=rhp_t(i,j,kmax)+rho

      XCO2 = atm_bio(1)*1e-6
      pCO2Air = atm_bio(1)*1e-6 * (Pression-PressionVap) *9.869233e-6*1.e6 ! en µatm, conversion de la pression de Pa en atm

! 1. The Schmidt number (Jahne et al 1984 -> Essayer avec les valeurs de Wanninkhof 2014)
      Schmidt = 2073.1-125.62*Temp+3.6276*TT -0.043219*TTT


! K0 is in mol/kg/atm, or µmol/kg/µatm
! pCO2 is in µatm, rho is in kg/m3

      T = Temp + 273.15

      lnK = -60.2409 + 9345.17/T                 &
            +23.3585* LOG(T/100) +                     &
            Sal*( 0.023517-2.3656e-4*T+4.7036e-7*T*T  )
      HenryCst = EXP(lnK)


!-----------------------------------------------------------------------------
! 2. The piston velocty, in cm/hour according to Wanninkhof (1992)
! Journ. Geophys. Res. 97: 7327-7382.

! Either we use bulk formulae and we know the wind             
      if(iairsea/=1)then
        Wind=sqrt(uwind_t(i,j,1)**2+vwind_t(i,j,1)**2)
      else
! Or fluxes are prescribed and we compute the wind using the windstress 
        x1=wstress_w(i,j)
        Wind=-5.8546*x1**6+46.5303*x1**5-143.3061*x1**4+216.0191*x1**3-167.8414*x1**2+73.5452*x1+2.2032
      endif


       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(i2==969.and.j2==145) then
       print*,'dans AirSeaFlux,Wind,pCO2W(i,j),Temp,Schmidt,wstress,iairsea,mask'
       print*,'dans AirSeaFlux',Wind,pCO2W(i,j),Temp,Schmidt,wstress_w(i,j),iairsea
!       print*,'dans AirSeaFlux',i,j,mask_t(i,j,kmax)
!       print*,'dans AirSeaFlux',atm_bio(1),pCO2Air,PressionVap,Pression
!       print*,'dans AirSeaFlux',bio_t(i,j,kmax,idic)/(rhp_t(i,j,kmax)+rho)*1000., &
!       (bio_t(i,j,kmax,ialkalinity) &
!                            +bio_t(i,j,kmax,iAmmonium)  &
!                            -bio_t(i,j,kmax,iNitrate)   &
!                            -bio_t(i,j,kmax,iPhosphate)) &
!                            /(rhp_t(i,j,kmax)+rho)*1000.
       print*,'dans AirSeaFlux spH,pHCO2',spH(i,j,kmax),-LOG10(10**(-spH(i,j,kmax))*1e3/(rhp_t(i,j,max0(k,kmin_w(i,j)))+rho))
!       print*,'dans AirSeaFlux T,S,Pho,Sil',tem_t(i,j,kmax,1),sal_t(i,j,kmax,1), &
!                           bio_t(i,j,kmax,iPhosphate)/(rhp_t(i,j,max0(k,kmin_w(i,j)))+rho)*1000., &
!                           bio_t(i,j,kmax,isilice)/(rhp_t(i,j,max0(k,kmin_w(i,j)))+rho)*1000.

       endif

      if(ichoixFSCO2.eq.1) then
! Wanninkhof and McGillis 1999
! steady/short-term wind
      Piston = (0.0283 * (wind**3)) * (Schmidt / 660.)**(-0.5)     ! Wanninkhof and McGillis 1999, ce qui a été utilisé jusqu'à présent

      elseif(ichoixFSCO2.eq.2) then
! Wanninkhof 1992 used in most of models
!       Piston = 2.5*(0.5246+0.016256*temp+0.00049946*TT) +kbul*0.3*wind*wind &
!            *1/SQRT(Schmidt/660)
       Piston = kbul*0.31*wind*wind*1/SQRT(Schmidt/660)


      elseif(ichoixFSCO2.eq.3) then
! Wanninkhof et al. 2009  
!      Piston = kbul*0.31*wind*wind*1/SQRT(Schmidt/660)
       Piston = (3+0.1*Wind + 0.064*Wind**2. + 0.011 * (Wind**3.)) * (Schmidt /660.)**(-0.5)

      elseif(ichoixFSCO2.eq.4) then
! Liss and Merlivat 1986
      if(Wind.le.3.6) then
      Piston = 0.17*Wind * (Schmidt / 600.)**(-0.6667)
      elseif(Wind.gt.3.6.and.Wind.lt.13) then
      Piston = ((2.85*Wind) -9.65) * (Schmidt/600)**(-0.5)
      else
      Piston =  (5.9*Wind-49.3) * (Schmidt/600)**(-0.5)
      endif

      elseif(ichoixFSCO2.eq.5) then
! Nightingale et al. 2000
      Piston = (0.222*Wind**2+0.333*Wind)* (Schmidt / 600.)** (-0.5)

! Wanninkhof et al 2014
      elseif(ichoixFSCO2.eq.6) then
      Piston = (0.251 * (Wind**2.)) * (Schmidt / 660.)**(-0.5)  ! 0.31 short-term wind mise à jour de W92


      elseif(ichoixFSCO2.eq.7) then
!! Ho et al. 2006
      Piston = (0.266*Wind**2.)* (Schmidt / 600.)** (-0.5)

!----------------------------------------------------
      elseif(ichoixFSCO2.ge.8.and.ichoixFSCO2.le.10) then
! Liang et al 2013 avec paramétrisation de 2 types de bulles (généralement + faible que
! les autres)

      if(ichoixFSCO2.eq.8) betag=1 ! Liang et al 2013
      if(ichoixFSCO2.eq.9) betag=0.53 ! Yang et al 2017
      if(ichoixFSCO2.eq.10) betag=0.29 ! Bushinsky and Emerson 2018 

! Mass transfer
        if(Wind<11) then
        Cdg=0.0012
        elseif(Wind>=11.and.Wind<=20) then
        Cdg=(0.49+0.065*Wind)*1e-3
        else
        Cdg=0.0018
        endif

        Uastar=Cdg**0.5*Wind
        Kw= 1.3*1e-4*1.04*Uastar*(Schmidt / 660.)**(-0.5) ! m s-1
        Fs= Kw * HenryCst * (pCO2Air - pCO2Wat) *1e-3*Rho0 !  en mmol m-2 s-1 avec pCO2 en µatmHenryCst en µmol/kg/µatm, Kw en m s-1 et Rho0 en kg/m3

! Echanges dus à l'eclatement  des petites bulles (d'apres Atamanchuck
! et al 2020 : contribution faible)

!       rho_a=1.293
        Md=0.0289652 
        Mv=0.018016
        rho_a=(Pression*Md+PressionVap*Mv)/(8.31*TemperatureK)
        Uwstar= sqrt (rho_a /density_l)*Uastar
        Kc = 5.56  * Uwstar**3.86
        Fc =  betag * Kc * XCO2 *1e3 ! en mmol m-2 s-1


!       i2=i+par%timax(1)      ! ds le repere global
!       j2=j+par%tjmax(1)
!       if(i2==557.and.j2==486) then
!       print*,'dans AirSeaFlux'
!       print*,'dans AirSeaFlux',Wind,pCO2W(i,j),Temp,Schmidt,wstress_w(i,j),iairsea
!       print*,'dans AirSeaFlux',i,j,mask_t(i,j,kmax)
!       print*,'dans AirSeaFlux',atm_bio(1),pCO2Air,PressionVap,Pression,rho_a,density_l
!       endif


! Echanges dus aux grosses bulles poussées dans l'océan

        Kpg =  5.5  * Uwstar**2.76 * (Schmidt/660)**(-0.6667)
        deltaP = 1.5244*Uwstar**(1.06)
        Fp = betag * Kpg*((1+deltaP)*HenryCst * (pCO2Air - pCO2Wat) *1e-3*Rho0)  ! en mmol m-2 s-1


!----------------------------------------------------
      elseif(ichoixFSCO2.eq.11) then

! Stanley et al 2009
! Diffusive
       Kw= 0.97*8.6e-7*Wind**2*(Schmidt / 660.)**(-0.5) ! m s-1
       Fs= Kw * HenryCst * (pCO2Air - pCO2Wat) *1e-3*Rho0  !en mmol m-2 s-1 avec pCO2 en µatm, HenryCst en µmol/kg/µatm, Kw en m s-1 et Rho0 en kg/m3        
   

! Bubbles completly trapped  
      Xi=XCO2 !ppm !0.20946 pour 02 !faction molaire du gaz 
      Pression=pss_w(i,j,1)
      PressionVap=q2_t(i,j,1)*pss_w(i,j,1)/0.622

! R=8.31 en m3 Pa mol-1 K-1

        if(Wind.gt.2.27) then
!       Fc=9.1e-11*(Wind-2.27)**3*(Pression-PressionVap)*Xi/8.31/TemperatureAir*1e3 !       ! en mmol m-2 s-1
        Fc=9.1e-11*(Wind-2.27)**3.*(Pression-PressionVap)*Xi/8.31/TemperatureK*1e3 ! en mmol m-2 s-1
        else
        Fc=0.
        endif

! Bubbles partially dissolved
        if(Wind.gt.2.27) then

! Bunsen coefficient of solubility (Weiss, 1970)
!    Wanninkhof 1992 : coef bunsen=K0/(rho*volume molaire du gaz) et 
! " The Ostwald solubility coefficient a is related to Bunsen coef B by a = B/T (T in
! degrees Kelvin)." 

!        lnbunsen= B_A1+ B_A2*(100./TemperatureK) + B_A3 * log(TemperatureK/100.) +     &
!         (Salinity * (B_B1 + B_B2 * (TemperatureK/100.) + B_B3 * (TemperatureK/100.)**2.))
!        Bunsen = exp(lnbunsen)
!

         Bunsen = HenryCst * density_l * 1e-3 * 22.414*0.99498

!! size of bubbles
          if(Wind.ge.0.37) then
          zbub = 0.15*Wind-0.55
          else
          zbub = 0.
          endif

       visc=1e-4*(17.91-0.5381*Temperature+0.00694*Temperature**2+0.02305*Salinity)/density_l
       Diff=visc/Schmidt ! en m2/s


        Fp=2.3e-3*Bunsen*Diff**(0.6667)*(Wind-2.27)**3.*(Pression-PressionVap)*Xi/8.31/TemperatureK*1e3*(1+density_l*9.81*zbub/(Pression-PressionVap)-pCO2Wat/pCO2Air)
!!        ! en mmol m-2 s-1
!        Fp=2.3e-3*(Wind-2.27)**3*bunsen*(Diff*1e-4)**(0.6667)*Pression*Xi/8.31/TemperatureK*1e3*(1+density_l*9.81*zbub/Pression-Oxygen/OxygenSat)
!
!!         if(i+par%timax(1)==569.and.j+par%tjmax(1)==484)then


!       i2=i+par%timax(1)      ! ds le repere global
!       j2=j+par%tjmax(1)
!       if(i2==557.and.j2==486) then
!       print*,'dans AirSeaFlux'
!       print*,'dans AirSeaFlux',Wind,pCO2W(i,j),Temp,Schmidt,wstress_w(i,j),iairsea
!       print*,'dans AirSeaFlux',i,j,mask_t(i,j,kmax)
!       print*,'dans AirSeaFlux',atm_bio(1),pCO2Air,PressionVap,Pression,rho_a,density_l
!       print*,'dans AirSeaFlux',Fs,Fp,Fc
!       endif

        endif 
!----------------------------------------------------
       elseif(ichoixFSCO2.eq.12) then
! Woolf 1997
! Mass transfer
        if(Wind<11) then
        Cdg=0.0012;
        elseif(Wind>=11.and.Wind<=20) then
        Cdg=(0.49+0.065*Wind)*1e-3;
        else
        Cdg=0.0018;
        endif

        Uastar=Cdg**0.5*Wind;
        DeltaE=0.01*(Wind/49)**2 ! commentaire Dariia Atamanchuk : 
! "eq.A4 where U10i is the ten meter wind speed at which the value of De is
! 1%. In their model U10,N2 7.2 m s21,
! U10,O25 9.0 m s21 and U10,CO2 49 m s21." vient de Emmerson and
! Bushinsky (2016) A4 p 4372

! Bunsen coefficient of solubility (Weiss, 1970)
!    Wanninkhof 1992 : coef bunsen=K0/(rho*volume molaire du gaz) et 
! " The Ostwald solubility coefficient a is related to Bunsen coef B by
! a = B/T (T in degrees Kelvin)."
! Ostwald solubility
        BetaX = HenryCst * density_l * 1e-3 * 22.414*0.99498 / TemperatureK 
! Ostwald solubility, alpha * T/Tstp
! Dariia Atamanchuck prend 23.1 pour Vm et D. Nicholson 22.414*0.99498

        Kw= 1.57*1e-4*Uastar*(Schmidt / 600.)**(-0.5) ! m s-1
        Fs= Kw * HenryCst * (pCO2Air - pCO2Wat) *1e-3*Rho0  
!en mmol m-2 s-1 avec pCO2 en µatm, HenryCst en µmol/kg/µatm, Kw en m s-1 et Rho0 en kg/m3)

        Kbw = 9.4e-3*Wind**3.41/BetaX*(1+(14*BetaX*(Schmidt)**(-0.5))**(-0.8333))**(-1.2)/3600/100; !page 196 (w97), eq. a7 (e16) in m s-1

        Ft= (Kw + Kbw) *((1+DeltaE)*pCO2Air-pCO2Wat)* HenryCst *1e-3*Rho0 !  en mmol m-2 s-1

        Fp = Ft - Fs  !  en mmol m-2 s-1 avec oxygen en mmol m-3 et kw en m s-1

        Fc=0.

       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(i2==557.and.j2==486) then
       print*,'dans AirSeaFlux'
       print*,'dans AirSeaFlux',Wind,pCO2W(i,j),Temp,Schmidt,wstress_w(i,j),iairsea
       print*,'dans AirSeaFlux',i,j,mask_t(i,j,kmax)
       print*,'dans AirSeaFlux',atm_bio(1),pCO2Air,PressionVap,Pression,rho_a,density_l
       print*,'dans AirSeaFlux',Fs,Fp,Fc
       endif



       endif
 

!     Piston = (2.5*(0.5246+0.016256*temp+0.00049946*TT) +15*wind*wind) &
!      *1/SQRT(Schmidt/660)

!-----------------------------------------------------------------------------
! The piston velocity in cm/hour according to Wanninkhof and McGillis 1999
! long term wind (>1month)
        !Piston = (1.09*wind - 0.33*wind*wind + 0.078*wind*wind*wind) &
        !               *1/sqrt(Schmidt/660)

        !Piston = (4.09*wind - 0.01*wind*wind + 2.9*wind*wind*wind) &
        !               *1/sqrt(Schmidt/660)

!-----------------------------------------------------------------------------
! The piston velocity in cm/hour according to Nightingale 2000

        !Piston = (0.23*wind*wind + 0.1*wind) &
        !               *1/sqrt(Schmidt/660)
!-----------------------------------------------------------------------------
! The piston velocity in cm/hour according to McGillis 2001

        !Piston = (3.3 + 0.026*wind*wind*wind) &
        !               *1/sqrt(Schmidt/660)


! Faycal pour l'oxygene
!      Kw= (0.0283 * (Wind**3)) * (Sc / 660.)**(-0.5)     ! Wanninkhof and McGillis 1999
     ! Kw= (0.031 * (Wind**2)) * (Sc / 660.)**-0.5      ! Wanninkhof 1992 used in most of models

! 3. convert from cm/hour to m/sec
      Piston = Piston *0.01/3600.

!------------------------!
! Transfer in mol/kg/sec !
!------------------------!

! pCO2 Air in µatm
!      Jour=iteration3d*dti_fw/86400
!      pi= atan(1)*4
!      pCO2Air = MeanpCO2air* AmplitudepCO2air *sin(Jour * 2*pi/PeriodpCO2air + PhasepCO2air)
!      pCO2Air =391.


      if(ichoixFSCO2.lt.8) then
      TransferRate = Piston * HenryCst * (pCO2Air - pCO2Wat)
! Convert to mmol/m2/sec
       CO2_AirSeaExchange(i,j) = TransferRate /1e3*Rho0
       ASF2D(i,j) = ASF2D(i,j) + CO2_AirSeaExchange(i,j)
       fluxbio_w(i,j,iDIC,2) = TransferRate /1e3*Rho0

       else
       CO2_AirSeaExchange(i,j) = Fs + Fp + Fc
       ASF2D(i,j) = ASF2D(i,j) + CO2_AirSeaExchange(i,j)
       fluxbio_w(i,j,iDIC,2) = CO2_AirSeaExchange(i,j)

       endif




      endif !mask

      enddo
      enddo

      
  
      return
      end subroutine AirSeaCO2Flux







