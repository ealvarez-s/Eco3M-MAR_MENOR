      subroutine oxygen_surfaceflux
      use module_principal
      use module_biology
      use ModuleDeclaration
      use module_parallele
      implicit none
      real ox_A0,ox_A1,ox_A2,ox_A3,ox_A4,ox_A5,ox_B0,ox_B1,ox_B2,ox_B3,ox_A,ox_B,ox_C,ox_D    &
          ,ox_C0,B_A1,B_A2,B_A3,B_B1,B_B2,B_B3                                                &
          ,Temperature,Salinity,OxygenSat,Oxygen,Ts,Kw,Sc,Wind,Pression,PressionVap,density_l &
          ,TemperatureK,TemperatureAir,Xi,zbub,BetaX,DeltaE,Kbw,vt,dm
      real Cdg,Uastar,Uwstar,Fs,Fp,Fc,Kpg,Kc,rho_a,XO2,deltaP,betag,bunsen,lnbunsen,diff

          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! oxygene à saturation !!!!!!!!!!!!!!!!!!!!!!!!

      TPS_O2F_2D=TPS_O2F_2D+1



! ecriture de l'algorithme de l'oxygene à saturation
! equation of Garcia and Gordon 1992


! les constantes:

ox_A0= 2.00907
ox_A1= 3.22014
ox_A2= 4.05010
ox_A3= 4.94457
ox_A4= -0.256847
ox_A5= 3.88767

ox_B0= -6.24523e-3
ox_B1= -7.37614e-3
ox_B2= -1.03410e-2
ox_B3= -8.17083e-3

ox_C0= -4.88682e-7

! coeff pour le calcul du nombre de schmidt: from Keeling et al 1998
ox_A=1638
ox_B=81.83
ox_C=1.483
ox_D=0.008004

! Bunsen coeffcient (Weiss, 1970)
B_A1= -58.3877
B_A2= 85.8079
B_A3= 23.8439

B_B1= -0.034892
B_B2=  0.015568
B_B3= -0.0019387

!et d'après Winninkhof 1992
!ox_A=1953.4
!ox_B=128
!ox_C=3.9012
!ox_D=0.050091


      do i=1,imax
      do j=1,jmax

      if(mask_t(i,j,kmax)==1)then

       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)


      Temperature= max(tem_t(i,j,kmax,1),0.)
      TemperatureK=Temperature+273.15
      Salinity=max(sal_t(i,j,kmax,1),0.)
      Oxygen=bio_t(i,j,kmax,iOxygen)        
!      print*,'oxygensurf1',i2,j2,pss_w(i,j,1)
!      print*,'oxygensurf2',i2,j2,q2_t(i,j,1)
      Pression=pss_w(i,j,1)
      
      PressionVap=q2_t(i,j,1)*pss_w(i,j,1)/0.622 
      TemperatureAir=teta2_t(i,j,1)
      Density_l=rhp_t(i,j,kmax)+rho

! ici Temperature est la temperature de surface
! ici Oxygen et l'oxygen de surface

      Ts = log((298.5-Temperature)/(273.15+Temperature))

! x1 au lieu de l
      x1= (ox_A0 + ox_A1*Ts + ox_A2 * Ts**2 + ox_A3 *Ts**3 + ox_A4 *Ts**4 + ox_A5 * Ts**5) +    &
         (Salinity * (ox_B0 + ox_B1 * Ts + ox_B2 * Ts**2 + ox_B3 * Ts**3)) + (ox_C0 * Salinity**2)

       OxygenSat = (1000./22.3916) * exp(x1)

       OxySat(i,j) = OxygenSat

! 1000=rho0? unité?

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! piston velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  piston velocity
!  equation de Wanninkhof and McGillis 1999
      
      !Sc is the Schmidt number

      Sc = ox_A - (ox_B * Temperature) + (ox_C * (Temperature**2)) - (ox_D * (Temperature**3))

      if(iairsea/=1)then    ! soit on est en formule bulk et on connait le vent
      Wind=sqrt(uwind_t(i,j,1)**2+vwind_t(i,j,1)**2)
      else
! soit on est en flux prescrit et on calcule le vent a partir du windstress 
      x1=wstress_w(i,j)
      Wind=-5.8546*x1**6+46.5303*x1**5-143.3061*x1**4+216.0191*x1**3-167.8414*x1**2+73.5452*x1+2.2032
      endif
 
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)

!       if(i2.eq.549.and.j2.eq.504) &
!        print*,'wind',uwind_t(i,j,1),vwind_t(i,j,1),wstress_w(i,j)


! Il existe plusieurs estimations du coef d'echange.
! Nightingale 2000 se trouve entre W92 et LM86

      if(ichoixFS.eq.1) then
! Wanninkhof and McGillis 1999
      Kw = (0.0283 * (Wind**3)) * (Sc / 660.)**(-0.5)     ! Wanninkhof and McGillis 1999, ce qui a été utilisé jusqu'à présent

      elseif(ichoixFS.eq.2) then
! Wanninkhof 1992 used in most of models
      Kw = (0.31 * (Wind**2)) * (Sc / 660.)**(-0.5)  ! 0.31 short-term wind   

      elseif(ichoixFS.eq.3) then
! Wanninkhof et al 2009
      Kw = (3+0.1*Wind + 0.064*Wind**2 + 0.011 * (Wind**3)) * (Sc / 660.)**(-0.5)

      elseif(ichoixFS.eq.4) then
! Liss and Merlivat 1986
      if(Wind.le.3.6) then
      Kw = 0.17*Wind * (Sc / 600.)**(-0.6667)
      elseif(Wind.gt.3.6.and.Wind.lt.13) then
      Kw = ((2.85*Wind) -9.65) * (Sc/600)**(-0.5)
      else
      Kw =  (5.9*Wind-49.3) * (Sc/600)**(-0.5)
      endif

      elseif(ichoixFS.eq.5) then
! Nightingale et al. 2000
      Kw = (0.222*Wind**2+0.333*Wind)* (Sc / 600.)** (-0.5)

! Wanninkhof et al 2014
      elseif(ichoixFS.eq.6) then
      Kw = (0.251 * (Wind**2)) * (Sc / 660.)**(-0.5)  ! 0.31 short-term wind mise à jour de W92

      elseif(ichoixFSCO2.eq.7) then
!! Ho et al. 2006
      Kw = (0.266*Wind**2.)* (Sc / 600.)** (-0.5)


      elseif(ichoixFS.ge.8.and.ichoixFS.le.10) then
! Liang et al 2013 avec paramétrisation de 2 types de bulles (généralement + faible que
! les autres)

      if(ichoixFS.eq.8) betag=1 ! Liang et al 2013
      if(ichoixFS.eq.9) betag=0.53 ! Yang et al 2017
      if(ichoixFS.eq.10) betag=0.29 ! Bushinsky and Emerson 2018 

! Mass transfer
        if(Wind<11) then
        Cdg=0.0012
        elseif(Wind>=11.and.Wind<=20) then
        Cdg=(0.49+0.065*Wind)*1e-3
        else
        Cdg=0.0018
        endif

        Uastar=Cdg**0.5*Wind
        Kw= 1.3*1e-4*1.04*Uastar*(Sc / 660.)**(-0.5) ! m s-1
        Fs= Kw * (OxygenSat -  Oxygen) !  en mmol m-2 s-1 avec oxygen en mmol m-3 et kw en m s-1

! Echanges dus à l'eclatement  des petites bulles

        rho_a=1.293
        Uwstar= sqrt (rho_a /density_l)*Uastar
        XO2= 20.95/100;
        Kc = 5.56  * Uwstar**3.86
        Fc =  betag * Kc * XO2 *1e3 ! en mol m-2 s-1

! Echanges dus aux grosses bulles poussées dans l'océan

        Kpg =  5.5  * Uwstar**2.76 * (Sc/660)**(-0.6667)
        deltaP = 1.5244*Uwstar**(1.06)
        Fp = betag * Kpg*((1+deltaP)*OxygenSat-Oxygen) ! en mmol m-2 s-1


        elseif(ichoixFS.eq.11) then
! Stanley et al 2009
! Diffusive
        Kw= 0.97*8.6e-7*Wind**2*(Sc / 660.)**(-0.5) ! m s-1
        Fs= Kw * (OxygenSat -  Oxygen) !  en mmol m-2 s-1 avec oxygen en mmol m-3 et kw en m s-1

! Bubbles completly trapped  
        Xi=0.20946 !faction molaire du gaz 

        if(Wind.gt.2.27) then
!       Fc=9.1e-11*(Wind-2.27)**3*(Pression-PressionVap)*Xi/8.31/TemperatureAir*1e3 ! en mmol m-2 s-1
        Fc=9.1e-11*(Wind-2.27)**3*(Pression-PressionVap)*Xi/8.31/TemperatureK*1e3 ! en mmol m-2 s-1
        else
        Fc=0.
        endif 

! Bubbles partially dissolved
        if(Wind.gt.2.27) then

! Bunsen coefficient of solubility (Weiss, 1970)
        lnbunsen= B_A1+ B_A2*(100./TemperatureK) + B_A3 * log(TemperatureK/100.) +     &
         (Salinity * (B_B1 + B_B2 * (TemperatureK/100.) + B_B3 * (TemperatureK/100.)**2.)) 
        Bunsen = exp(lnbunsen)

! oxygen diffusion Wilke and Chang (1955) http://compost.css.cornell.edu/oxygen/oxygen.diff.water.html
        Diff =  7.4e-8*TemperatureK*(2.26*18.)**(0.5)/(1.*(25.6)**(0.6)) ! cm2/s

! size of bubbles
          if(Wind.ge.0.37) then
          zbub = 0.15*Wind-0.55
          else
          zbub = 0.
          endif

!        Fp=2.3e-3*(Wind-2.27)**3*bunsen*(Diff*1e-4)**(0.6667)*Pression*Xi/8.31/TemperatureAir*1e3*(1+density_l*9.81*zbub/Pression-Oxygen/OxygenSat) ! en mmol m-2 s-1
        Fp=2.3e-3*(Wind-2.27)**3*bunsen*(Diff*1e-4)**(0.6667)*Pression*Xi/8.31/TemperatureK*1e3*(1+density_l*9.81*zbub/Pression-Oxygen/OxygenSat)

!         if(i+par%timax(1)==569.and.j+par%tjmax(1)==484)then
!         write(6,*)'S09d',2.3e-3*(Wind-2.27)**3,bunsen*(Diff*1e-4)**(2/3),Pression*Xi/8.31/TemperatureAir*1e3,(1+density_l*9.81*zbub/Pression-Oxygen/OxygenSat)
!         write(6,*)'S09d2',bunsen,Diff,(Diff*1e-4)**(2/3),(Diff*1e-4)**(0.6667)
!         endif
        else
        Fp=0.
        endif

        elseif(ichoixFS.eq.12) then
! Woolf 1997 (aide avec les codes de Dariia Atamanchuk)
! Mass transfer
        if(Wind<11) then
        Cdg=0.0012;
        elseif(Wind>=11.and.Wind<=20) then
        Cdg=(0.49+0.065*Wind)*1e-3;
        else
        Cdg=0.0018;
        endif

        Uastar=Cdg**0.5*Wind;
        DeltaE=0.01*(Wind/9)**2 ! page 189 (w97), eq. a4 (e16)
        BetaX = 2.53e-2*TemperatureK/273.15 ! Ostwald solubility, alpha * T/Tstp

        Kw= 1.57*1e-4*Uastar*(Sc / 600.)**(-0.5) ! m s-1
        Fs= Kw * ((1+DeltaE)*(OxygenSat -  Oxygen)) !  en mmol m-2 s-1 avec oxygen en mmol m-3 et kw en m s-1

        Kbw = 9.4e-3*Wind**3.41/BetaX*(1+(14*BetaX*(Sc)**(-0.5))**(-0.8333))**(-1.2)/3600/100; !page 196 (w97), eq. a7 (e16)
        Fc =  Kbw * ((1+DeltaE)*(OxygenSat -  Oxygen)) !  en mmol m-2 s-1 avec oxygen en mmol m-3 et kw en m s-1

        Fp=0.



        elseif(ichoixFS.eq.13) then
! Vagle et al 2010 (codes de Dariia Atamanchuk)

! Mass transfer
              if(Wind.le.3.6) then
      Kw = 0.17*Wind * (Sc / 600.)**(-0.6667) ! Liss and Merlivat
      elseif(Wind.gt.3.6.and.Wind.le.12) then
      Kw = ((2.85*Wind) -9.65) * (Sc/600)**(-0.5) ! Liss and Merlivat
      else
      Kw =  0.078*(Wind**2.62)*(Sc/660)**(-0.5)
      endif

! bubble penetration depth as function of wind speed, given on page 6

        dm = -0.83+0.481*Wind
        if(dm.lt.0.7.or.Wind.lt.12) dm=0.
        vt = 0.377*(dm**2.33) !eq. 20, coefficients on page 12

        Xi=0.20946
!        R=.082057338
! volume of air from ideal gas law; used to convert air injection to mol/m^2 of O2
!        V = R*TemperatureK/Pression

        Fs = Kw/360000. *(OxygenSat - Oxygen) 

        Fc = vt*Pression/1.01325e5*Xi/8.2057338/TemperatureK/360000/24*1000 ! pas clair pour moi

        Fp = 0.

        endif !ichoixFS

        if(ichoixFS.lt.8) then
!    fluxbio_w(i,j,iOxygen,2) = Kw/360000. * (abs( Oxygen-OxygenSat ))    !    conversion de cm/hr en m/s 
        fluxbio_w(i,j,iOxygen,2) = Kw/360000. * (OxygenSat -  Oxygen)    ! conversion de cm/hr en m/s modif caroline a verifier convention 
        else
        fluxbio_w(i,j,iOxygen,2) = Fs + Fp + Fc
        endif !ichoixFS

     O2flux2D(i,j)=O2flux2D(i,j)+fluxbio_w(i,j,iOxygen,2)

     if(i+par%timax(1)==569.and.j+par%tjmax(1)==484)then
!     write(6,*)'L13',fluxbio_w(i,j,iOxygen,2),Fs,Fp,Fc,Wind,Sc,Uastar,Kw
!     stop
!    write(6,*)'S09',fluxbio_w(i,j,iOxygen,2),Fs,Fp,Fc,Wind,Pression,PressionVap
!     write(6,*)'S09',fluxbio_w(i,j,iOxygen,2),Fs,Fp,Fc,Wind
!     write(6,*)'S09',q2_t(i,j,1),teta2_t(i,j,1),TemperatureAir,OxygenSat,Oxygen
!     write(6,*)'S09',density_l,bunsen,Xi,TemperatureK,Diff,zbub,Sc,Kw
     
!     stop
!     write(6,*)'S09',fluxbio_w(i,j,iOxygen,2),Fs,Fp,Fc,O2flux2D(i,j),Wind,Sc
     endif


     endif

     enddo
     enddo



!-----------------------------------------------------------------------------



  
     return
     end subroutine oxygen_surfaceflux







