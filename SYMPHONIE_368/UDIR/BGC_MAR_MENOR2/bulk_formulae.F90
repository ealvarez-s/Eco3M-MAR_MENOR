      subroutine bulk_formulae(ichoix,tsst_,tvar_,relativewind_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 362 - last update: 06-01-23
!______________________________________________________________________
      use module_principal
      implicit none
      integer  ichoix,tsst_,tvar_
      real*4 relativewind_
#ifdef synopsis
       subroutinetitle='bulk_formulae'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version Date      Description des modifications
!         12/07/02  mise en service
!         06/01/04  Claude ajoute la parametrisation "MFSTEP" pour le calcul du
!                   flux net infrarouge (net longwave flux)
!         26/01/06  Un calcul different dans le cas de la houle
!                   + reorganisation de la routine en differentes subroutines
!         18/04/06  fonctions compatibles avec double precision
!         10/05/06  fonctions compatibles avec double precision
!         02/06/06  debug lie à faute sur equation 77 Geernaert page 109
!                   Au passage on note que diverses valeurs ont ete proposees
!                   pour beta. Ici beta=7 mais on trouve aussi beta=5. A tester.
!         08/06/06  rationnalisation du calcul (via approximations legitimes)
!                   proposée par Guillaume Reffray dans le cadre du projet
!                   "OPA Symphonique" conduit par MERCATOR & POC.
!         29/07/06  amelioration des commentaires
!         17/01/08  IWAVE renommé IWVE
! 2010.8  19-03-10  Tsurf prise au temps t-1. kmax remplace kmax
! 2010.9  06-06-10  model_wave renommé waveforcing
! 2010.12 09-09-10  on appelle iter_bulk_2006 quelque soit la valeur de iwve
! S26     13-10-14  choix du 4eme argument de tem_t pour passer la SST du modele
!                   meteo. Par defaut 4eme arg=0 sinon (=2) il s'agit de la SST
!                   du modele meteo
!         23-11-15  ote une partie du courant de surface au vent
!         27-11-15  cas 1Dv
!         01-03-16  passer relativewind_ en argument de bukl_relative_wind
!         24-04-16  la constante rhoair remplacee par airdensity=p/(R.T)
!         25-04-16  SUite du point precedent: passer par une autre variable
!                   que rhoair car celle ci existe en tant que constante dans d'autres
!                   subroutines ce qui crEE un probleme de conservation mpi. On crEe
!                   donc airdensity
!         02-05-16  Ajout des formules bulk de Moon
!         08-05-16  - Ajout des formules bulk COARE
!                   - Calcul si mask_t=1
!                   - Seuils vent min et ustar min revus
!         28-11-16  flag_net_ir permet de distinger cas IR atmos du cas IR atmos/ocean
!         25-01-18  si flag_wstressbulk=0 ne pas calculer le stress avec les formules bulk
! v309    18-09-21  formule CORE: 
!                   - Filtre EMA 
!                   - suppresion du cas neutre si forfait max atteint
!                   - forfait max grand A l'etat initial
! v362    06-01-23  Ajout formules bulk ecume
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! (°v°)
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! (°O°)
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      ! (°L°)
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

! A la toute premiere iteration (i.e. à l'initialisation) la premiere
! estimation fait l'hypothèse de la stabilité:

! remove a fraction (relativewind_) of the surface current to the wind velocity
      if(relativewind_/=0.)call bukl_relative_wind(tsst_,tvar_,relativewind_) !01-03-16

! Initial state, reset,...:
      if(ichoix==1) then !-reset-> ! Initial State>
       if(bulk_scheme==bulk_core) call initial_bulk_2006(tsst_,tvar_)
       if(bulk_scheme==bulk_moon) call initial_bulk_moon(tsst_,tvar_) !02-05-16
       if(bulk_scheme==bulk_coare)call initial_bulk_coare(tsst_,tvar_)!08-05-16
       if(bulk_scheme==bulk_ecume)call initial_bulk_ecume(tsst_,tvar_)!06-01-23
      endif              !-reset-> ! Initial State>

! Update wstress_w slhf_w sshf_w:
       if(bulk_scheme==bulk_core) call iter_bulk_2006(tsst_,tvar_) !09-09-10
       if(bulk_scheme==bulk_moon) call iter_bulk_moon(tsst_,tvar_) !02-05-16
       if(bulk_scheme==bulk_coare)call iter_bulk_coare(tsst_,tvar_)!08-05-16
       if(bulk_scheme==bulk_ecume)call iter_bulk_ecume(tsst_,tvar_)!06-01-23
! Update snsf_w:
      call snsf_bulk(tsst_,tvar_) !02-05-16

! Repartition de l'info depuis points _Z vers points _X et _Y
      if(flag_abl==0) then !zzzzzz>
        call ztoxy_bulk ! standard case
      else                 !zzzzzz>
        call bulk_stress_w2uv(tvar_) ! Couche Limite Atmospherique !24-11-15
      endif                !zzzzzz>

      if(flag_1dv==1)call bulk_1dv !27-11-15

      end subroutine bulk_formulae

!................................................................................

      subroutine snsf_bulk(tsst_,tvar_)
      use module_principal
      implicit none
      integer tsst_,tvar_

      if(iairsea==2) then !ooooo> !03-05-16
       if(flag_net_ir==0) then !pmxpmx> !28-11-16
        do j=1,jmax ; do i=1,imax
         snsf_w(i,j,1)=snsf_w(i,j,1)-stefan*(tem_t(i,j,kmax,tsst_)+celsius2kelvin)**4 
        enddo       ; enddo
       endif                   !pmxpmx> !28-11-16
!      call snsf_netparam_bulk(tsst_,tvar_)

      return
      endif               !ooooo>

      if(iairsea==3) then !ooooo>
       do j=1,jmax ; do i=1,imax

       x1=pss_w(i,j,tvar_)*1.608*q2_t(i,j,tvar_)/    & ! Pression de vapeur fonction de humidite specifique
        (1.+0.608*q2_t(i,j,tvar_))/                  & ! Queney page 115               060104
        1.e2                                           ! exprimee en millibars         060104

       snsf_w(i,j,1)=-(                                               & ! flux IR net                   060104
        (1.-0.75*(snsf_w(i,j,1)**3.4))*                               & ! SNSF : ici nebulosite         060104
        (stefan*teta2_t(i,j,tvar_)**4*                                & ! Temperature de l'air (on approxime temperature etat e
        (0.4-0.05*sqrt(x1))+                                          & !060104
        4.*stefan*(teta2_t(i,j,tvar_)**3)*                            & !060104
        (celsius2kelvin+tem_t(i,j,kmax,tsst_)-teta2_t(i,j,tvar_))))     !060104

       enddo       ; enddo
      return
      endif               !ooooo>

      stop ' Err 119 snsf_bulk iairsea unrecognized' 
      end subroutine snsf_bulk
!................................................................................
#ifdef bidon
      subroutine snsf_netparam_bulk(tsst_,tvar_)
      use module_principal ; use module_parallele
      implicit none
      integer tsst_,tvar_

! Cette routine a pour but de corriger le flux IR descendant en raison
! du changement de SST, de T2M et de Q2M. Il existe des formules donnant
! le flux IR net A la surface de l'ocean, voir la reference suivante:
! BIGNAMI et al, JOURNAL OF GEOPHYSICAL RESEARCH, VOL.  1OO, NO. C2, PAGES 2501-2514, FEBRUARY 15, 1995
! Table 3
! Clark formula: eps.sigma.SST**4*(0.39-0.05*sqrt(e))*(1-0.69*c**2)
!            +4.*eps*sigma.SST**3*(SST-TA)

! J'ai testE cette routine en omettant la couverture nuageuse que je
! n'avais pas A disposition et comparant au flux net calculE avec le
! snsf et le sima*SST**4 et je constate que ce n'est pas terrible (sans
! doute du fait de l'omission des nuages) A savoir qu'il y a de gros
! biais. En revanche c'est plutot bien correlE. C'est donc correct pour
! sans servir de cette formule pour construire une perturbation. On la
! calcule donc 2 fois (avec les champs UNCPL et les champs CPL) et la
! difference sert A incrementer le SNSF fourni par le modele meteo.
! Principe general: si la SST refroidit, l'air se refroidit aussi et
! par consequent IRdown (positif vers le bas) doit diminuer (increment
! negatif)


      x2=(elapsedtime_now           -airseafile_prvtime(t2m_id))      &
        /(airseafile_nextime(t2m_id)-airseafile_prvtime(t2m_id))
! Bornage car l'attente de la synchro subcycling peut faire sortir
! (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x1=1.-x2 ; x0=x1
      do j=1,jmax ; do i=1,imax
      if(mask_t(i,j,kmax)==1) then !>>>>>

! Clark formula: eps.sigma.SST**4*(0.39-0.05*sqrt(e))*(1-0.69*c**2)
!            +4.*eps*sigma.SST**3*(SST-TA)

!---------------
! Flux IRnet calculE A partir des champs du modele METEO non-couple
       x3=x1*teta0_t(i,j,1)+x2*teta0_t(i,j,2) ! SST modele meteo
       x4=x0*   q2_t(i,j,0)+x2*   q2_t(i,j,2) ! Q2 non couple
       x5=x0*teta2_t(i,j,0)+x2*teta2_t(i,j,2) ! T2 non couple
!      x4=q2_t(i,j,1)
!      x5=teta2_t(i,j,1)

! x6=pressure vapor in mbars
       x6=0.01*pss_w(i,j,tvar_)*x4 &
                  /(0.622+0.378*x4)

        x7=                  &
        -0.98*stefan*x3**3*( &
                     x3      &
                    *(0.39-0.05*sqrt(x6))                        &    
                  ! *(1.-0.69*c**2) ! if assuming no cloud cover
                 +4.*(x3-x5) )

!---------------
! Flux IRnet calculE A partir des champs S couplEs A SABL

! x6=pressure vapor in mbars
       x6=0.01*pss_w(i,j,tvar_)*q2_t(i,j,tvar_) &
                  /(0.622+0.378*q2_t(i,j,tvar_))

!                      xy_t(i,j,2)=                              &
        x8=                                                      &
        -0.98*stefan*(tem_t(i,j,kmax,tsst_)+celsius2kelvin)**3*( &
                     (tem_t(i,j,kmax,tsst_)+celsius2kelvin)      &
                    *(0.39-0.05*sqrt(x6))                        &    
                  ! *(1.-0.69*c**2) ! if assuming no cloud cover
             +4.*(    tem_t(i,j,kmax,tsst_)+celsius2kelvin       &
                   -teta2_t(i,j,tvar_))                        )

          if(tvar_==1)write(10+par%rank,*)snsf_w(i,j,1),x7,x8

! x7=IRnet UNCPL = IRdownUNCPL-sig.SST0**4
! x8=IRnet   CPL = IRdownCPL  -sig.SST **4
! IRdownCPL - IRdown = x8-x7+sigma*( SST**4 - SST0**4 )


! IRdownCPL - IRdown = est l'increment de flux IRdown lie aux
! modifications de la SST, de la temperature et l'humidite de l'air
      if(tvar_==1)write(300+par%rank,*)x8-x7+stefan*((tem_t(i,j,kmax,tsst_)+celsius2kelvin)**4-x3**4) &
                                      ,tem_t(i,j,kmax,tsst_)+celsius2kelvin-x3

      endif                        !>>>>>
      enddo       ; enddo


#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
      if(tvar_==1)stop 'vivi'
      end subroutine snsf_netparam_bulk
#endif
!................................................................................     

      subroutine initial_bulk_coare(tsst_,tvar_)
      use module_principal
      implicit none
      integer tsst_,tvar_
! Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and J. B. Edson, 2003 
! Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm. 
! J. Climate, 16, 571–591, 
! doi:10.1175/1520-0442 

! Constants involved in Coare Scheme
       z2m=2.
       z10m=10.
       pss0=1.e5
       karman=0.4
       boltz=1.380658e-23
       planck=6.6260755e-34
       avogadro= 6.0221367e+23
       xmd=28.9644e-3
       xmv=18.0153e-3
       xrd=avogadro*boltz/xmd
       xrv=avogadro*boltz/xmv
       xcpd=7.*xrd /2.
       xcpv=4.*xrv
       xcl=4.218e+3
       celsius2kelvin=273.16 ! Oui ca finit bien par un 6 dans SURFEX....
       xlvtt=2.5008e+6
       xestt=611.14
       xgamw=(xcl-xcpv)/xrv
       xbetaw=(xlvtt/xrv)+(xgamw*celsius2kelvin)
       xalpw=log(xestt)+(xbetaw/celsius2kelvin)+(xgamw *log(celsius2kelvin))

      end subroutine initial_bulk_coare

!................................................................................
      subroutine initial_bulk_ecume(tsst_,tvar_) !06-01-23
      use module_principal
      use module_ecume
      implicit none
      integer tsst_,tvar_
     z2m=2.
     z10m=10.
     pss0=1.e5
     XRD=6.0221367E+23*1.380658E-23/28.9644E-3
     XRV=6.0221367E+23*1.380658E-23/18.0153E-3
!    write(6,*)'XRD= ',XRD
!    write(6,*)'XRV= ',XRV
     XGAMW = (4.218E+3 - 4.* XRV) / XRV
     XBETAW = (2.5008E+6/XRV) + (XGAMW * 273.16)
     XALPW = LOG(611.14) + (XBETAW /273.16) + (XGAMW *LOG(273.16))
     XP00= pss0
     XCPD=7.* (6.0221367E+23*1.380658E-23 / 28.9644E-3)  /2.
!    write(6,*)'XCPD= ',XCPD
     XKARMAN =0.4
     PZREF=z2m
     PUREF=z10m
     XZ0=0.001
     XLVTT  = 2.5008E+6
     XCPV   = 4.* XRV
     XCL    = 4.218E+3
     XG = 9.80665
     XTT      = 273.16
     NITERMAX = 5
     NITERSUP = 5
     NITERFL=NITERMAX
     ZUTU = 40.0
     ZUTT = 14.4
     ZUTQ = 10.0
     ZETV    = XRV/XRD-1.0   !~0.61 (cf Liu et al. 1979)
     ZLMOMIN = -200.0
     ZLMOMAX = 0.25
     ZBTA    = 16.0
     ZGMA    = 7.0 
     ZP00    = 1013.25E+02
     XPI      = 2.*ASIN(1.)
     ZSQR3   = SQRT(3.0)
     ZCOEFU = (/ 1.00E-03, 3.66E-02, -1.92E-03, 2.32E-04, -7.02E-06,  6.40E-08 /)
     ZCOEFT = (/ 5.36E-03, 2.90E-02, -1.24E-03, 4.50E-04, -2.06E-05,       0.0 /)
     ZCOEFQ = (/ 1.00E-03, 3.59E-02, -2.87E-04,      0.0,       0.0,       0.0 /)
!

     ZDUSR0   = 1.E-06
     ZDTSR0   = 1.E-06
     ZDQSR0   = 1.E-09

     ZCDIRU = ZCOEFU(1) + 2.0*ZCOEFU(2)*ZUTU + 3.0*ZCOEFU(3)*ZUTU**2   &
                        + 4.0*ZCOEFU(4)*ZUTU**3 + 5.0*ZCOEFU(5)*ZUTU**4
     ZCDIRT = ZCOEFT(1) + 2.0*ZCOEFT(2)*ZUTT + 3.0*ZCOEFT(3)*ZUTT**2   &
                        + 4.0*ZCOEFT(4)*ZUTT**3
     ZCDIRQ = ZCOEFQ(1) + 2.0*ZCOEFQ(2)*ZUTQ

     ZORDOU = ZCOEFU(0) + ZCOEFU(1)*ZUTU + ZCOEFU(2)*ZUTU**2 + ZCOEFU(3)*ZUTU**3   &
                        + ZCOEFU(4)*ZUTU**4 + ZCOEFU(5)*ZUTU**5
     ZORDOT = ZCOEFT(0) + ZCOEFT(1)*ZUTT + ZCOEFT(2)*ZUTT**2 + ZCOEFT(3)*ZUTT**3   &
                        + ZCOEFT(4)*ZUTT**4
     ZORDOQ = ZCOEFQ(0) + ZCOEFQ(1)*ZUTQ + ZCOEFQ(2)*ZUTQ**2

     XUNDEF   = -999.

     PCD   = XUNDEF
     PCH   = XUNDEF
     PCE   = XUNDEF
     PCHN   = XUNDEF
     PCEN   = XUNDEF
     PCDN  = XUNDEF
     ZUSR  = XUNDEF
     ZTSR  = XUNDEF
     ZQSR  = XUNDEF
     ZTAU  = XUNDEF
     ZHF   = XUNDEF
     ZEF   = XUNDEF
     PSFTH   = XUNDEF
     PSFTQ   = XUNDEF
     PUSTAR  = XUNDEF
     PRESA   = XUNDEF
     PRI     = XUNDEF


     OPCVFLX  = .TRUE.
     LPRECIP = .FALSE.
     LPWEBB = .FALSE.
  
     ZTAUR=0.
     ZRF=0.
     ZEFWEBB=0.


           end subroutine initial_bulk_ecume
!................................................................................
      subroutine iter_bulk_ecume(tsst_,tvar_) !06-01-23
      use module_principal; use module_parallele
      use module_ecume
      implicit none
      integer tsst_,tvar_

! passage variable modele à ECUME
      do j=1,jmax
      do i=1,imax
      if(  mask_t(i,j,kmax+1).eq.1) then !%%%%%%%%%%%%%%%%%%%%%%%%%%%>

      sst_kelvin=tem_t(i,j,kmax,tsst_)+celsius2kelvin
      prs_atm_z= pss_w(i,j,tvar_)*exp(-grav/xrd/(0.5*(teta2_t(i,j,tvar_)+sst_kelvin)*( 1.+((xrv/xrd)-1.)*q2_t(i,j,tvar_) )) * (z2m))
      PSST=sst_kelvin    ! passage variable ECUME
      PPS=pss_w(i,j,tvar_) ! passage variable ECUME
      PPA=PPS
     delta_u=max(0.5,sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2))
     PVMOD=delta_u      ! passage variable ECUME
     ZVMOD=MAX(PVMOD , 0.1 * MIN(10.,PUREF) )
     PTA=teta2_t(i,j,tvar_)  ! passage variable ECUME
     PQA=q2_t(i,j,tvar_)     ! passage variable ECUME
     ZSSS=sal_t(i,j,kmax,1)  ! passage variable ECUME
     PRAIN=0.               !!!!!!!! A VOIR
     PRHOA=PPS / XRD / PTA / ( 1.+((XRV/XRD)-1.)*PQA )

!    write(6,*)'sst_kelvin,prs_atm_z,PPS,PTA,PQA ',sst_kelvin,prs_atm_z,PPS,PTA,PQA

!---------------------------
! CAS TEST
!     PSST=280.
!     PPS=100000. ; PPA=PPS
!     PQA=5.559e-3
!     delta_u=20. ; PVMOD=delta_u
!     ZVMOD=MAX(PVMOD , 0.1 * MIN(10.,PUREF) )
!     PTA=286.
!     PRHOA=PPS / XRD / PTA / ( 1.+((XRV/XRD)-1.)*PQA )
!     ZSSS=38.3
! la reponse
!!!! sshf_w(i,j,1),slhf_w(i,j,1),wstress_w(i,j)   287.1609 -44.80267 0.883524835109711

!---------------------------


      PEXNA=(prs_atm_z/pss0)**(xrd/xcpd)
      exner_sea_z=(pss_w(i,j,tvar_)/pss0)**(xrd/xcpd)
      PEXNX=exner_sea_z

!    write(6,*)'PEXNA,PEXNX,pss_w(i,j,tvar_),prs_atm_z,pss0,xrd,xcpd ', &
!     PEXNA,PEXNX,pss_w(i,j,tvar_),prs_atm_z,pss0,xrd,xcpd

     ZFOES  = 0.98*EXP( XALPW - XBETAW/PSST  - XGAMW*LOG(PSST ) )
     PQSAT  = XRD/XRV*ZFOES /PPS  / (1.+(XRD/XRV-1.)*ZFOES /PPS )
     ZPA     = XP00*(PEXNA **(XCPD/XRD))
     ZFOESA  = EXP( XALPW - XBETAW/PTA  - XGAMW*LOG(PTA ) )
     ZQSATA  = XRD/XRV*ZFOESA /ZPA  / (1.+(XRD/XRV-1.)*ZFOESA /ZPA )

!       2.2. Gradients at the air-sea interface
!
      ZDU  = ZVMOD                !one assumes u is measured / sea surface current
!     write(6,*)'PTA,PEXNA,PSST,PEXNX ',PTA,PEXNA,PSST,PEXNX
      ZDT  = PTA /PEXNA -PSST /PEXNX 
!     write(6,*)'ZDT ',ZDT
      ZDQ  = PQA -PQSAT 

!       2.3. Latent heat of vaporisation
!
      ZLVA  = XLVTT+(XCPV-XCL)*(PTA  -XTT)                !of pure water at atm level
      ZLVS  = XLVTT+(XCPV-XCL)*(PSST -XTT)                !of pure water at sea surface
      ZLVS  = ZLVS *(1.0-1.00472E-3*ZSSS)            !of seawater at sea surface

!       2.4. Specific heat of moist air (Businger 1982)
!
!ZCPA = XCPD*(1.0+(XCPV/XCPD-1.0)*PQA)
      ZCPA = XCPD

!       2.4b Kinematic viscosity of dry air (Andreas 1989, CRREL Rep. 89-11)
!
      ZVISA  = 1.326E-05*(1.0+6.542E-03*(PTA -XTT)+8.301E-06*(PTA -XTT)**2   &
           -4.84E-09*(PTA -XTT)**3)

!       2.5. Initial guess
!
      ZDDU  = ZDU 
      ZDDT  = ZDT 
      ZDDQ  = ZDQ 
      ZDDU  = SIGN(MAX(ABS(ZDDU ),10.0*ZDUSR0),ZDDU )
      ZDDT  = SIGN(MAX(ABS(ZDDT ),10.0*ZDTSR0),ZDDT )
      ZDDQ  = SIGN(MAX(ABS(ZDDQ ),10.0*ZDQSR0),ZDDQ )
!
      JCV   = -1
      ZUSR  = 0.04*ZDDU 
      ZTSR  = 0.04*ZDDT 
      ZQSR  = 0.04*ZDDQ 
      ZDELTAU10N  = ZDDU 
      ZDELTAT10N  = ZDDT 
      ZDELTAQ10N  = ZDDQ 
      JITER  = 99

!
! In the following, we suppose that Richardson number PRI < XRIMAX
! If not true, Monin-Obukhov theory can't (and therefore shouldn't) be applied !
!-------------------------------------------------------------------------------
!
!       3.   ITERATIVE LOOP TO COMPUTE U*, T*, Q*.
!       ------------------------------------------
!
     DO JJ=1,NITERFL
      IF (JCV  == -1) THEN
        ZUSR0 =ZUSR 
        ZTSR0 =ZTSR 
        ZQSR0 =ZQSR 
        IF (JJ == NITERMAX+1 .OR. JJ == NITERMAX+NITERSUP) THEN
          ZDELTAU10N  = 0.5*(ZDUSTO +ZDELTAU10N )    !forced convergence
          ZDELTAT10N  = 0.5*(ZDTSTO +ZDELTAT10N )
          ZDELTAQ10N  = 0.5*(ZDQSTO +ZDELTAQ10N )
          IF (JJ == NITERMAX+NITERSUP) JCV =3
        ENDIF
        ZDUSTO  = ZDELTAU10N 
        ZDTSTO  = ZDELTAT10N 
        ZDQSTO  = ZDELTAQ10N 
     
!
!       3.1. Neutral parameter for wind speed (ECUME_V8 formulation)
!
       IF (ZDELTAU10N  <= ZUTU) THEN
         ZPARUN  = ZCOEFU(0) + ZCOEFU(1)*ZDELTAU10N       &
                                  + ZCOEFU(2)*ZDELTAU10N **2   &
                                  + ZCOEFU(3)*ZDELTAU10N **3   &
                                  + ZCOEFU(4)*ZDELTAU10N **4   &
                                  + ZCOEFU(5)*ZDELTAU10N **5
       ELSE
         ZPARUN  = ZCDIRU*(ZDELTAU10N -ZUTU) + ZORDOU
       ENDIF
       PCDN  = (ZPARUN /ZDELTAU10N )**2
!
!       3.2. Neutral parameter for temperature (ECUME_V8 formulation)
!
       IF (ZDELTAU10N  <= ZUTT) THEN
         ZPARTN  = ZCOEFT(0) + ZCOEFT(1)*ZDELTAU10N       &
                                  + ZCOEFT(2)*ZDELTAU10N **2   &
                                  + ZCOEFT(3)*ZDELTAU10N **3   &
                                  + ZCOEFT(4)*ZDELTAU10N **4
       ELSE
         ZPARTN  = ZCDIRT*(ZDELTAU10N -ZUTT) + ZORDOT
       ENDIF
!
!       3.3. Neutral parameter for humidity (ECUME_V8 formulation)
!
       IF (ZDELTAU10N  <= ZUTQ) THEN
         ZPARQN  = ZCOEFQ(0) + ZCOEFQ(1)*ZDELTAU10N       &
                                  + ZCOEFQ(2)*ZDELTAU10N **2
       ELSE
         ZPARQN  = ZCDIRQ*(ZDELTAU10N -ZUTQ) + ZORDOQ
       ENDIF
!
!       3.4. Scaling parameters U*, T*, Q*
!
       ZUSR  = ZPARUN 
       ZTSR  = ZPARTN *ZDELTAT10N /ZDELTAU10N 
       ZQSR  = ZPARQN *ZDELTAQ10N /ZDELTAU10N 

!       3.5. Obukhovs stability param. z/l following Liu et al. (JAS, 1979)
!
! For U
       ZLMOU = PUREF *XG*XKARMAN*(ZTSR /PTA    &
            +ZETV*ZQSR /(1.0+ZETV*PQA ))/ZUSR **2
! For T/Q
       ZLMOT = ZLMOU*(PZREF /PUREF )
       ZLMOU = MAX(MIN(ZLMOU,ZLMOMAX),ZLMOMIN)
       ZLMOT = MAX(MIN(ZLMOT,ZLMOMAX),ZLMOMIN)

!
!       3.6. Stability function psi (see Liu et al, 1979 ; Dyer and Hicks, 1970)
!            Modified to include convective form following Fairall (unpublished)
!
! For U
      IF (ZLMOU == 0.0) THEN
        ZPSI_U = 0.0
      ELSEIF (ZLMOU > 0.0) THEN
        ZPSI_U = -ZGMA*ZLMOU
      ELSE
        ZCHIK  = (1.0-ZBTA*ZLMOU)**0.25
        ZPSIK  = 2.0*LOG((1.0+ZCHIK)/2.0)  &
                  +LOG((1.0+ZCHIK**2)/2.0) &
                  -2.0*ATAN(ZCHIK)+0.5*XPI
        ZCHIC  = (1.0-12.87*ZLMOU)**(1.0/3.0)       !for very unstable conditions
        ZPSIC  = 1.5*LOG((ZCHIC**2+ZCHIC+1.0)/3.0)   &
                  -ZSQR3*ATAN((2.0*ZCHIC+1.0)/ZSQR3) &
                  +XPI/ZSQR3
        ZPSI_U = ZPSIC+(ZPSIK-ZPSIC)/(1.0+ZLMOU**2) !match Kansas & free-conv. forms
      ENDIF
      ZPSIU = ZPSI_U
! For T/Q
      IF (ZLMOT == 0.0) THEN
        ZPSI_T = 0.0
      ELSEIF (ZLMOT > 0.0) THEN
        ZPSI_T = -ZGMA*ZLMOT
      ELSE
        ZCHIK  = (1.0-ZBTA*ZLMOT)**0.25
        ZPSIK  = 2.0*LOG((1.0+ZCHIK**2)/2.0)
        ZCHIC  = (1.0-12.87*ZLMOT)**(1.0/3.0)       !for very unstable conditions
        ZPSIC  = 1.5*LOG((ZCHIC**2+ZCHIC+1.0)/3.0)   &
                  -ZSQR3*ATAN((2.0*ZCHIC+1.0)/ZSQR3) &
                  +XPI/ZSQR3
        ZPSI_T = ZPSIC+(ZPSIK-ZPSIC)/(1.0+ZLMOT**2) !match Kansas & free-conv. forms
      ENDIF
      ZPSIT = ZPSI_T

!
!       3.7. Update air-sea gradients
!
      ZDDU  = ZDU 
      ZDDT  = ZDT 
      ZDDQ  = ZDQ 
      ZDDU  = SIGN(MAX(ABS(ZDDU ),10.0*ZDUSR0),ZDDU )
      ZDDT  = SIGN(MAX(ABS(ZDDT ),10.0*ZDTSR0),ZDDT )
      ZDDQ  = SIGN(MAX(ABS(ZDDQ ),10.0*ZDQSR0),ZDDQ )
      ZLOGUS10   = LOG(PUREF /10.0)
      ZLOGTS10   = LOG(PZREF /10.0)
      ZDELTAU10N  = ZDDU -ZUSR *(ZLOGUS10-ZPSI_U)/XKARMAN
      ZDELTAT10N  = ZDDT -ZTSR *(ZLOGTS10-ZPSI_T)/XKARMAN
      ZDELTAQ10N  = ZDDQ -ZQSR *(ZLOGTS10-ZPSI_T)/XKARMAN
      ZDELTAU10N  = SIGN(MAX(ABS(ZDELTAU10N ),10.0*ZDUSR0),   &
                              ZDELTAU10N )
      ZDELTAT10N  = SIGN(MAX(ABS(ZDELTAT10N ),10.0*ZDTSR0),   &
                              ZDELTAT10N )
      ZDELTAQ10N  = SIGN(MAX(ABS(ZDELTAQ10N ),10.0*ZDQSR0),   &
                            ZDELTAQ10N )
!
!       3.8. Test convergence for U*, T*, Q*
!
     IF (ABS(ZUSR -ZUSR0 ) < ZDUSR0 .AND.   &
         ABS(ZTSR -ZTSR0 ) < ZDTSR0 .AND.   &
         ABS(ZQSR -ZQSR0 ) < ZDQSR0) THEN
       JCV  = 1                                     !free convergence
       IF (JJ >= NITERMAX+1) JCV  = 2               !leaded convergence
     ENDIF
     JITER  = JJ
     ENDIF
!
     ENDDO
!
!-------------------------------------------------------------------------------
!
!       4.   COMPUTATION OF TURBULENT FLUXES AND EXCHANGE COEFFICIENTS.
!       ---------------------------------------------------------------
!
!
!       4.1. Surface turbulent fluxes
!            (ATM CONV.: ZTAU<<0 ; ZHF,ZEF<0 if atm looses heat)
!
       ZTAU  = -PRHOA *ZUSR **2
       ZHF   = -PRHOA *ZCPA *ZUSR *ZTSR 
       ZEF   = -PRHOA *ZLVS *ZUSR *ZQSR 
!!  print*, PRHOA ,ZPARUN 
!
!       4.2. Exchange coefficients PCD, PCH, PCE
!
       PCD  = (ZUSR /ZDDU )**2
       PCH  = (ZUSR *ZTSR )/(ZDDU *ZDDT )
       PCE  = (ZUSR *ZQSR )/(ZDDU *ZDDQ )
     !
!       4.3. Stochastic perturbation of turbulent fluxes
!
!  IF( OPERTFLUX )THEN
!    ZTAU  = ZTAU * ( 1. + PPERTFLUX  / 2. )
!    ZHF   = ZHF *  ( 1. + PPERTFLUX  / 2. )
!    ZEF   = ZEF *  ( 1. + PPERTFLUX  / 2. )
!  ENDIF
!
!
!-------------------------------------------------------------------------------
!       5.   COMPUTATION OF FLUX CORRECTIONS DUE TO RAINFALL.
!            (ATM conv: ZRF<0 if atm. looses heat, ZTAUR<<0)
!       -----------------------------------------------------
!
     IF (LPRECIP) THEN
!
!       5.1. Momentum flux due to rainfall (ZTAUR, N/m2)
!
! See pp3752 in FBR96.
        ZTAUR = -0.85*PRAIN*PVMOD
!
!       5.2. Sensible heat flux due to rainfall (ZRF, W/m2)
!
! See Eq.12 in GoF95 with ZCPWA as specific heat of water at atm level (J/kg/K),
! ZDQSDT from Clausius-Clapeyron relation, ZDWAT as water vapor diffusivity 
! (Eq.13-3 of Pruppacher and Klett, 1978), ZDTMP as heat diffusivity, and ZBULB
! as wet-bulb factor (Eq.11 in GoF95).
!
        ZTAC   = PTA -XTT
        ZCPWA  = 4217.51 -3.65566*ZTAC +0.1381*ZTAC**2       &
                  -2.8309E-03*ZTAC**3 +3.42061E-05*ZTAC**4   &
                  -2.18107E-07*ZTAC**5 +5.74535E-10*ZTAC**6
        ZDQSDT = (ZLVA *ZQSATA )/(XRV*PTA **2)
        ZDWAT  = 2.11E-05*(ZP00/PPA )*(PTA /XTT)**1.94
        ZDTMP  = (1.0+3.309E-03*ZTAC-1.44E-06*ZTAC**2)   &
                  *0.02411/(PRHOA *ZCPA )
        ZBULB  = 1.0/(1.0+ZDQSDT*(ZLVA *ZDWAT)/(ZCPA *ZDTMP))
        ZRF  = PRAIN *ZCPWA*ZBULB*((PSST -PTA )   &
                +(PQSAT -PQA )*(ZLVA *ZDWAT)/(ZCPA *ZDTMP))
!
     ENDIF
!
!-------------------------------------------------------------------------------
!
!       6.   COMPUTATION OF WEBB CORRECTION TO LATENT HEAT FLUX (ZEFWEBB, W/m2).
!       ------------------------------------------------------------------------
!
! See Eq.21 and Eq.22 in FBR96.
     IF (LPWEBB) THEN
        ZWW = (1.0+ZETV)*(ZUSR *ZQSR )   &
               +(1.0+(1.0+ZETV)*PQA )*(ZUSR *ZTSR )/PTA 
        ZEFWEBB  = -PRHOA *ZLVS *ZWW*PQA 
     ENDIF
!
!-------------------------------------------------------------------------------
!
!       7.   FINAL STEP : TOTAL SURFACE FLUXES AND DERIVED DIAGNOSTICS. 
!       ---------------------------------------------------------------
!
!       7.1. Richardson number
! je commente des choses pas utiles ici
!
!    CALL SURFACE_RI(PSST,PQSAT,PEXNX,PEXNA,PTA,PQA,   &
!               PZREF,PUREF,ZDIRCOSZW,PVMOD,PRI)
!
!
!       7.2. Friction velocity which contains correction due to rain
!
     ZUSTAR2 = -(ZTAU+ZTAUR)/PRHOA       !>>0 as ZTAU<<0 & ZTAUR<=0
!
     IF (LPRECIP) THEN
       PCD = ZUSTAR2/ZDDU**2
     ENDIF
     !
     PUSTAR = SQRT(ZUSTAR2)                    !>>0
!
!       7.3. Aerodynamical conductance and resistance
!
!    ZAC  (:) = PCH(:)*ZDDU(:)
!    PRESA(:) = 1.0/ZAC(:)
!
!       7.4. Total surface fluxes
!
! claude: je change les signes car on veut la convention inverse
     PSFTH =  -(ZHF+ZRF)   ! chaleur sensible
     PSFTQ = -(ZEF+ZEFWEBB) ! chaleur latente 
     PTAU =  PRHOA *PUSTAR * PUSTAR  ! Qté de mouvement
!!PSFTQ(:) = (ZEF(:)+ZEFWEBB(:))/ZLVS(:)
     sshf_w(i,j,1)=PSFTH
     slhf_w(i,j,1)=PSFTQ
     wstress_w(i,j)=PTAU

!     if(i+par%timax(1)==10.and.j+par%tjmax(1)==10)then
!      write(50,*)'sshf_w(i,j,1),slhf_w(i,j,1),wstress_w(i,j)', &
!       sshf_w(i,j,1),slhf_w(i,j,1),wstress_w(i,j)
!      write(60,*)iteration3d,PTA,PPS,PVMOD,PSST,PQA
!     endif

!
! claude: commente
!       7.5. Charnock number
!
!    ZCHARN = MIN(0.018,MAX(0.011,0.011+(0.007/8.0)*(ZDDU-10.0)))
!
!       7.6. Roughness lengths Z0 and Z0H over sea
!! Formulation si KZ0 == 2
!     PZ0SEA   = PUREF /EXP(XKARMAN*ZDDU /PUSTAR +ZPSIU )
!     Z0TSEA        = PZREF /EXP(XKARMAN*ZDDT /ZTSR   +ZPSIT )
!     Z0QSEA        = PZREF /EXP(XKARMAN*ZDDQ /ZQSR   +ZPSIT )
!     PZ0HSEA  = 0.5*(Z0TSEA+Z0QSEA)
!     PCHN  = (XKARMAN/LOG(PZREF /PZ0SEA ))*(XKARMAN/LOG(PZREF /Z0TSEA))
!     PCEN  = (XKARMAN/LOG(PZREF /PZ0SEA ))*(XKARMAN/LOG(PZREF /Z0QSEA))

     endif   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ENDDO
     ENDDO
     end subroutine iter_bulk_ecume
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     SUBROUTINE SURFACE_RI
!      ZTHVA =PTA /PEXNA *( 1.+(XRV/XRD-1.)*PQA  )
!      ZTHVS =PTG /PEXNX *( 1.+(XRV/XRD-1.)*PQS  )
!!
!      ZVMOD  = WIND_THRESHOLD(PVMOD ,PUREF )
!!
!! Richardson's number
!     PRI  = XG * PDIRCOSZW  * PUREF  * PUREF  &
!     * (ZTHVA -ZTHVS ) / (0.5 * (ZTHVA +ZTHVS ) ) &
!     / (ZVMOD *ZVMOD ) /PZREF 
!
!     END SUBROUTINE SURFACE_RI


      subroutine iter_bulk_coare(tsst_,tvar_)
      use module_principal ; use module_parallele
      implicit none
      integer tsst_,tvar_,loop_
      real x_claude
! Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and J. B. Edson, 2003 
! Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm. 
! J. Climate, 16, 571–591, 
! http://dx.doi.org/10.1175/1520-0442

!     do loop2_=0,nbdom-1
!     if(par%rank==loop2_) then
!     if(par%rank==0)write(6,*)'BIDOUILLE VENT'
!     uwind_t=0. ; vwind_t=0.

      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmax)==1) then !pmxpmx>

      z0_u=0.0001
      sst_kelvin=tem_t(i,j,kmax,tsst_)+celsius2kelvin

      x1=0.98*exp( xalpw - xbetaw/sst_kelvin - xgamw*log(sst_kelvin)  ) / pss_w(i,j,tvar_)
! vapor pressure reduction of 2% over saline seawater could have a significant 
! impact on the computation of surface latent heat flux under strong wind 
! conditions (Zeng et al, 1998). 
      qsat_sea_z= (xrd/xrv)*x1 / (1.+((xrd/xrv)-1.)*x1)
!
      prs_atm_z= pss_w(i,j,tvar_)*exp(-grav/xrd/(0.5*(teta2_t(i,j,tvar_)+sst_kelvin)*( 1.+((xrv/xrd)-1.)*q2_t(i,j,tvar_) )) * (z2m))
      airdensity=prs_atm_z / xrd /  teta2_t(i,j,tvar_) / ( 1.+((xrv/xrd)-1.)*q2_t(i,j,tvar_) )

      exner_atm_z=(prs_atm_z/pss0)**(xrd/xcpd)
      exner_sea_z=(pss_w(i,j,tvar_)/pss0)**(xrd/xcpd)

! set a non-zero value for the temperature gradient
      tem_atm_z=teta2_t(i,j,tvar_)
      if (  (teta2_t(i,j,tvar_)*exner_sea_z/exner_atm_z-sst_kelvin)  == 0. ) tem_atm_z=teta2_t(i,j,tvar_)-1e-3


      delta_u=max(0.5,sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2))
!
!
!      initial guess
!    

      delta_t= -(tem_atm_z/exner_atm_z) + (sst_kelvin/exner_sea_z) !potential temperature difference
      delta_q= qsat_sea_z - q2_t(i,j,tvar_)                        !specific humidity difference
!
      ustar_t(i,j,1)=0.035*delta_u

!     write(10+par%rank,*)0.
!     write(10+par%rank,*)ustar_t(i,j,1)

      visa= 1.326d-5*(1.+6.542d-3*(tem_atm_z-celsius2kelvin)+&
            8.301d-6*(tem_atm_z-celsius2kelvin)**2.-4.84d-9*(tem_atm_z-celsius2kelvin)**3.) !Andrea (1989) CRREL Rep. 89-11
! 
      z0_u= 0.011*ustar_t(i,j,1)*ustar_t(i,j,1)/grav+0.11*visa/ustar_t(i,j,1)

      x1  = (karman/log(z10m/Z0_u))**2
      x2  = 0.00115/sqrt(x1)
      z0_t= 10./exp(karman/x2)
!
      x1 =(karman/log(z10m/Z0_u))**2               !drag coefficient
      x2 = karman/log(z10m/z0_u)                   !temperature transfer coefficient
      x3 = karman*x2/x1                            !z/L vs Rib linear coef.
!
      x4 = -z10m/(600*0.004*1.2**3)                !saturation or plateau Rib
      x5  = -grav*z10m*(delta_t+zrvsrdm1*tem_atm_z*delta_q)/&
               (tem_atm_z*delta_u**2)
!
      if (x5<0.) then
         zl_10 = x3*x5/(1.+x5/x4)                   !Unstable G and F
      else
         zl_10 = x3*x5/(1.+27./9.*x5/x3)            !Stable 
      endif
!     
      zl_2=zl_10*z2m/z10m
!      
      if(zl_10>0.) then  !>>>
        psifunctu=-((1.+1.*zl_10)**1. + 0.6667*(zl_10-14.28)/exp(min(50.,0.35*zl_10)) + 8.525)
      else               !>>>
        x1=(1.-15.*zl_10)**0.25 
        x2= 2.*log((1.+x1)/2.)+log((1.+x1*x1)/2.)-2.*atan(x1)+2.*atan(1.)
        x3=(1.-10.15*zl_10)**0.3333
        x4=1.5*log((x3*x3+x3+1.)/3.)-(3.**0.5)*atan((2.*x3+1.)/(3.**0.5))+4.*atan(1.)/(3.**0.5)
        x5=zl_10*zl_10/(1.+zl_10*zl_10)
        psifunctu=(1.-x5)*x2+x5*x4
      endif 

      if(zl_2>0.) then  !>>>
        psifunctt=-((1.+2.*zl_2/3.)**1.5+0.6667*(zl_2-14.28)/exp(min(50.,0.35*zl_2)) + 8.525)
      else               !>>>
        x1=(1.-15.*zl_2)**0.5 
        x2= 2.*log((1.+x1)/2.)
        x3=(1.-34.15*zl_2)**0.3333
        x4=1.5*log((x3*x3+x3+1.)/3.)-(3.**0.5)*atan((2.*x3+1.)/(3.**0.5))+4.*atan(1.)/(3.**0.5)
        x5=zl_2*zl_2/(1.+zl_2*zl_2)
        psifunctt=(1.-x5)*x2+x5*x4
      endif    
!
!     ustar_t(i,j,1)= delta_u*karman/(log(z10m/z0_u)-psifunctu)
!     ustar_t(i,j,1)=max(delta_u*karman/(log(z10m/z0_u)-psifunctu),1.d-3)
      ustar_t(i,j,1)=    delta_u*karman/(log(z10m/z0_u)-psifunctu)       
      qstar_t(i,j,1)=   -delta_q*karman/(log( z2m/z0_t)-psifunctt)
      tetastar_t(i,j,1)=-delta_t*karman/(log( z2m/z0_t)-psifunctt)
      charnock = 0.011
      if (delta_u>10.) charnock = 0.011 + (0.018-0.011)*(delta_u-10.)/(18-10)
      if (delta_u>18.) charnock = 0.018

!     write(10+par%rank,*)ustar_t(i,j,1)

!      iterative loop
      loop1=3
      if(zl_10>50.) loop1=1
!     if(delta_u<=1.) loop1=10
!     if(delta_u<=0.5) then
!      loop1=0 ; ustar_t(i,j,1)=0.
!     endif

      do loop_=1,loop1 ! begin of iterative loop

      z0_u= charnock*ustar_t(i,j,1)*ustar_t(i,j,1)/grav+ 0.11*visa/ustar_t(i,j,1) !Smith 1988    
      z0_t= min(1.15d-4 , 5.5d-5/(z0_u*ustar_t(i,j,1)/visa)**0.6)
   
      zl_10= karman*grav*z10m*&
              (tetastar_t(i,j,1)*(1.+zrvsrdm1*q2_t(i,j,tvar_))+zrvsrdm1*tem_atm_z*qstar_t(i,j,1))/&
              (tem_atm_z*ustar_t(i,j,1)*ustar_t(i,j,1)*(1.+zrvsrdm1*q2_t(i,j,tvar_)))

      zl_2=zl_10*z2m/z10m
!      
      if(zl_10>0.) then  !>>>
        psifunctu=-((1.+1.*zl_10)**1. + 0.6667*(zl_10-14.28)/exp(min(50.,0.35*zl_10)) + 8.525)
      else               !>>>
        x1=(1.-15.*zl_10)**0.25 
        x2= 2.*log((1.+x1)/2.)+log((1.+x1*x1)/2.)-2.*atan(x1)+2.*atan(1.)
        x3=(1.-10.15*zl_10)**0.3333
        x4=1.5*log((x3*x3+x3+1.)/3.)-(3.**0.5)*atan((2.*x3+1.)/(3.**0.5))+4.*atan(1.)/(3.**0.5)
        x5=zl_10*zl_10/(1.+zl_10*zl_10)
        psifunctu=(1.-x5)*x2+x5*x4
      endif 

      if(zl_2>0.) then  !>>>
        psifunctt=-((1.+2.*zl_2/3.)**1.5+0.6667*(zl_2-14.28)/exp(min(50.,0.35*zl_2)) + 8.525)
      else               !>>>
        x1=(1.-15.*zl_2)**0.5 
        x2= 2.*log((1.+x1)/2.)
        x3=(1.-34.15*zl_2)**0.3333
        x4=1.5*log((x3*x3+x3+1.)/3.)-(3.**0.5)*atan((2.*x3+1.)/(3.**0.5))+4.*atan(1.)/(3.**0.5)
        x5=zl_2*zl_2/(1.+zl_2*zl_2)
        psifunctt=(1.-x5)*x2+x5*x4
      endif    
!
!        ustar_t(i,j,1)=max(delta_u*karman/(log(z10m/z0_u)-psifunctu),1.d-3)
         ustar_t(i,j,1)=    delta_u*karman/(log(z10m/z0_u)-psifunctu)
         qstar_t(i,j,1)=   -delta_q*karman/(log( z2m/z0_t)-psifunctt)
      tetastar_t(i,j,1)=   -delta_t*karman/(log( z2m/z0_t)-psifunctt)

!     write(10+par%rank,*)ustar_t(i,j,1)

      enddo ! fin de boucle iterative

        x_claude=0.125*max(delta_u,11.)-0.375   ! TIDE35_occ   ! 1.5 36_occ 37 38 40
!       x_claude=1.

        wstress_w(i,j)=airdensity*ustar_t(i,j,1)   *ustar_t(i,j,1) 
         sshf_w(i,j,1)=airdensity*ustar_t(i,j,1)*tetastar_t(i,j,1)*xcpd*x_claude
         slhf_w(i,j,1)=airdensity*ustar_t(i,j,1) *qstar_t(i,j,1)*xlvtt*x_claude


!-------------------------------------------------------------------------------

      endif                        !pmxpmx>
      enddo ! fin de boucle sur i
      enddo ! fin de boucle sur j

!      endif
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)
!#endif
!      enddo ! loop2_
!      stop 'koko'
      end subroutine iter_bulk_coare

!................................................................................

      subroutine initial_bulk_moon(tsst_,tvar_)
      use module_principal
      implicit none
      integer tsst_,tvar_
! Moon, I.-J., I. Ginis, T. Hara, and B. Thomas, 2007: A physics-based
! parameterization of air¿sea momentum flux at high wind speeds and its
! impact on hurricane intensity predictions. Mon. Wea. Rev.,135, 2869¿2878
! http://dx.doi.org/10.1175/MWR3432.1

! Constants involved in Moon Scheme
       z2m=2.
       z10m=10.
       pss0=1.e5
       karman=0.4
       boltz=1.380658e-23
       planck=6.6260755e-34
       avogadro= 6.0221367e+23
       xmd=28.9644e-3
       xmv=18.0153e-3
       xrd=avogadro*boltz/xmd
       xrv=avogadro*boltz/xmv
       xcpd=7.*xrd /2.
       xcpv=4.*xrv
       xcl=4.218e+3
       celsius2kelvin=273.16 ! Oui ca finit bien par un 6 dans SURFEX....
       xlvtt=2.5008e+6
       xestt=611.14
       xgamw=(xcl-xcpv)/xrv
       xbetaw=(xlvtt/xrv)+(xgamw*celsius2kelvin)
       xalpw=log(xestt)+(xbetaw/celsius2kelvin)+(xgamw *log(celsius2kelvin))
       zrvsrdm1=xrv/xrd-1. ! 0.607766

      end subroutine initial_bulk_moon

!................................................................................

      subroutine iter_bulk_moon(tsst_,tvar_)
      use module_principal ; use module_parallele
      implicit none
      integer tsst_,tvar_,loop_
! Moon, I.-J., I. Ginis, T. Hara, and B. Thomas, 2007: A physics-based
! parameterization of air¿sea momentum flux at high wind speeds and its
! impact on hurricane intensity predictions. Mon. Wea. Rev.,135, 2869¿2878
! http://dx.doi.org/10.1175/MWR3432.1

      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmax)==1) then !pmxpmx>

! Dans la cas oU la couche limite simplifiee SABL est utilisee on
! veille A ce que 0<q<qsat et on corrige delta_q en consequence
!     if(q2_t(i,j,tvar_)<0) then !sabl>
!      if(tvar_==1)q2delta_t(i,j,1)=q2delta_t(i,j,1)-q2_t(i,j,1) ! dq+q=0
!      q2_t(i,j,tvar_)=0.
!     endif                      !sabl>

      z0_u=0.0001
      sst_kelvin=tem_t(i,j,kmax,tsst_)+celsius2kelvin

      x1=0.98*exp( xalpw - xbetaw/sst_kelvin - xgamw*log(sst_kelvin)  ) / pss_w(i,j,tvar_)
! vapor pressure reduction of 2% over saline seawater could have a significant 
! impact on the computation of surface latent heat flux under strong wind 
! conditions (Zeng et al, 1998). 
      qsat_sea_z= (xrd/xrv)*x1 / (1.+((xrd/xrv)-1.)*x1)
      prs_atm_z= pss_w(i,j,tvar_)*exp(-grav/xrd/(0.5*(teta2_t(i,j,tvar_)+sst_kelvin)*( 1.+((xrv/xrd)-1.)*q2_t(i,j,tvar_) )) * (z2m))

      airdensity=prs_atm_z / xrd /  teta2_t(i,j,tvar_) / ( 1.+((xrv/xrd)-1.)*q2_t(i,j,tvar_) )

      exner_atm_z=(prs_atm_z/pss0)**(xrd/xcpd)
      exner_sea_z=(pss_w(i,j,tvar_)/pss0)**(xrd/xcpd)

! set a non-zero value for the temperature gradient
      tem_atm_z=teta2_t(i,j,tvar_)
      if ((teta2_t(i,j,tvar_)*exner_sea_z/exner_atm_z-sst_kelvin)==0.) tem_atm_z=teta2_t(i,j,tvar_)-1e-3


      delta_u=max(6.7e-5,sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2))
!
! specific humidity at saturation at the atm. level 
!     x1=exp( xalpw - xbetaw/tem_atm_z - xgamw*log(tem_atm_z)  ) / prs_atm_z
!     qsat_atm_z=(xrd/xrv)*x1/(1.+((xrd/xrv)-1.)*x1)

! Dans la cas oU la couche limite simplifiee SABL est utilisee on
! veille A ce que 0<q<qsat et on corrige delta_q en consequence
      

!
!      2.2       initial guess
!    

      delta_t= -(tem_atm_z/exner_atm_z) + (sst_kelvin/exner_sea_z) !potential temperature difference
      delta_q= qsat_sea_z - q2_t(i,j,tvar_)                        !specific humidity difference

!      2.4 initialization of scaling params (u*,t*,q*)
      ustar_t(i,j,1)=0.35*delta_u
      qstar_t(i,j,1)=    -delta_q*karman/log(z2m/z0_u)
      tetastar_t(i,j,1)=    -delta_t*karman/log(z2m/z0_u)

      delta_u_n=delta_u
      do loop_=1,5 ! begin of iterative loop

      zl_10=max(-200.,min(0.25,                                                                                 &
      karman*grav*z10m/tem_atm_z*(tetastar_t(i,j,1)*(1.+zrvsrdm1*q2_t(i,j,tvar_))+zrvsrdm1*tem_atm_z*qstar_t(i,j,1)) &
                                 /(ustar_t(i,j,1)*ustar_t(i,j,1))                                             &
                                 /(1.+zrvsrdm1*q2_t(i,j,tvar_))))


      zl_2=zl_10*z2m/z10m

      z0_u=0.0185*(ustar_t(i,j,1)*ustar_t(i,j,1))/grav

      if( delta_u_n>12.5)z0_u=0.001*(0.085*delta_u_n-0.58)

                  z0_q=z0_u
                  z0_t=z0_u

      if(zl_10>0.) then  !>>>
        psifunctu=-5.*zl_10
      else               !>>>
        x1=(1.-16.*zl_10)**.25
        psifunctu=2.*log((1.+x1)/2.)+log((1.+x1*x1)/2.)-2.*atan(x1)+2.*atan(1.)
      endif              !>>>

      if(zl_2>0.) then   !ooo>
        psifunctt=-5.*zl_2
      else               !ooo>
        x1=(1.-16.*zl_2)**.25
        psifunctt=2.*log((1.+x1)/2.)+log((1.+x1*x1)/2.)-2.*atan(x1)+2.*atan(1.)
      endif              !ooo>


      delta_u_n=max(delta_u*log(z10m/z0_u)/(log(z10m/z0_u)-psifunctu) , 0.1 )

         ustar_t(i,j,1)=max(delta_u*karman/(log(z10m/z0_u)-psifunctu),0.005)
      tetastar_t(i,j,1)=   -delta_t*karman/(log(z2m /z0_t)-psifunctt) 
         qstar_t(i,j,1)=   -delta_q*karman/(log(z2m /z0_q)-psifunctt) 

      enddo

      wstress_w(i,j)=  airdensity*ustar_t(i,j,1)   *ustar_t(i,j,1) 
         sshf_w(i,j,1)=airdensity*ustar_t(i,j,1)*tetastar_t(i,j,1)*xcpd
         slhf_w(i,j,1)=airdensity*ustar_t(i,j,1)   *qstar_t(i,j,1)*xlvtt

! cd=(ustar_t(i,j,1)/delta_u)**2.
! ch=ustar_t(i,j,1)*tetastar_t(i,j,1)/(delta_u*(tem_atm_z*exner_sea_z/exner_atm_z-sst_kelvin))
! ce=ustar_t(i,j,1)*qstar_t(i,j,1)/(delta_u*(qsat_sea_z-qsat_atm_z))

! cdn=(karman/log(z10m/z0_u))**2.
! chn=(karman/log(z10m/z0_u))*(karman/log(z10m/z0_t))
! cen=(karman/log(z10m/z0_u))*(karman/log(z10m/z0_q))

!-------------------------------------------------------------------------------

!     if(i+par%timax(1)==379.and.j+par%tjmax(1)==116) then
!       write(66,*)'------------------------------------------'
!       write(66,*)'tsst_ tvar_',tsst_,tvar_
!       write(66,*)'ustart     ',ustar_t(i,j,1)
!       if(zl_10>0.) then 
!         write(66,*)'1+5.*z50m/L',1.+5.*zl_10*50./10.
!         write(66,*)'K          ',karman*ustar_t(i,j,1)*50. &
!             /(1.+5.*(zl_10*50./10.))
!       else
!         write(66,*)'K instable ',karman*ustar_t(i,j,1)*50. &
!         *sqrt(1-9.*zl_10*50./10.)/0.74
!       endif
!     endif 

      xy_t(i,j,1)=zl_10

      endif                        !pmxpmx>

      enddo ! fin de boucle sur i
      enddo ! find de boucle sur j

      if(flag_abl==1) then !>>>>>>>>
       do j=1,jmax
       do i=1,imax

        if(xy_t(i,j,1)>0.) then !>>>>>>>>>
         kz_abl_w(i,j,1)=mask_t(i,j,kmax)*karman*zw_abl(2)*ustar_t(i,j,1) &
                              *(1.-zw_abl(2)/ablheight_t(i,j,tvar_))**1.5 &
               /(1.+5.*xy_t(i,j,1)*zw_abl(2)/10.)
        else                    !>>>>>>>>>
         kz_abl_w(i,j,1)=mask_t(i,j,kmax)*karman*zw_abl(2)*ustar_t(i,j,1) &
                              *(1.-zw_abl(2)/ablheight_t(i,j,tvar_))**1.5 &
           *sqrt(1.-9.*xy_t(i,j,1)*zw_abl(2)/10.)
        endif                   !>>>>>>>>>

! Cas 1ere couche = couche de surface = 10% couche limite
!       kz_abl_w(i,j,1)=mask_t(i,j,kmax)*karman*ustar_t(i,j,1)*0.1*ablheight_t(i,j,1) 

       enddo
       enddo
      endif                !>>>>>>>>
     

      end subroutine iter_bulk_moon
!................................................................................
#ifdef bidon
      subroutine initial_bulk(tsst_,tvar_)
      use module_principal
      implicit none
      integer tsst_,tvar_
#ifdef synopsis
       subroutinetitle='initial_bulk'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!*******************************************************************************
! INITIALISATION DE TABLEAUX ET DE CONSTANTES
! DEBUT:
!*******************************************************************************

!...............................................................................
! Constantes intervenant dans le calcul des flux air/mer
! calculés par les formules "bulk"
      karman=0.4            ! Constante de Von Karman
      stefan=5.67032e-8     ! constante de stephan (W/m**2/K**4)
      z_2=2.                ! niveau 2  mètres
      z_10=10.              ! niveau 10 mètres
      falpha=16.            ! coef fonction universelle instable (1-ALPHA Z/L)
      fbeta=7.              ! coef fonction universelle   stable (1+BETA  Z/L)
      cp_air=1002.93        ! capacité calorifique de l'air
      cen=1.2e-3            ! Smith 1980) Dalton number (à la neutralité)
!...............................................................................

! premiere estimation: on se place à la stabilité:
      chn=1.2e-3 ! Smith (1980) instable
      do j=1,jmax
      do i=1,imax

      uv_10=                                             & !vent à 10m instantanné
       max(un,                                           & !seuil pour eviter zero
       sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2)         & !
          )                                              !à partir compo. u & v

            teta10_t(i,j)=teta2_t(i,j,tvar_)
               q10_t(i,j)=   q2_t(i,j,tvar_)

! Pour USTAR:
! à la neutralité on remplace CD de la formule suivante:
! USTAR_Z(I,J,1)=(CD*UV_10**2)**0.5=SQRT(CD)*UV_10
! par CDN, sachant que CDN=MAX(0.93E-3,0.61E-3+6.3E-5*UV_10) ! Smith (1980)
! pour finalement obtenir:
      cdn=max(un*0.93e-3,0.61e-3+6.3e-5*uv_10)              ! Smith (1980)
      ustar_t(i,j,1)=sqrt(cdn)*uv_10

! Pour QSTAR:
! à la neutralité on remplace CE de la formule suivante (Geernaert p111):
! QSTAR_Z(I,J,1)=(CE*UV_10*(Q10_Z(I,J)-Q_0))/USTAR_Z(I,J,1)
! par CEN, sachant que Q_0 est donné par:
! trouver Q_0
!............> début:
       pvs_0=                                      & ! Pres Vap Saturante 0 metre
       610.78*exp(17.27*tem_t(i,j,kmax,tsst_)          & !
                      /(tem_t(i,j,kmax,tsst_)+237.29)) ! Nb: 237.29=273.15-35.86

      r_0=                   & ! rapport de melange saturant à 0 metre
        0.622*pvs_0          & ! 0.622 * Pression Vapeur Saturante 0 metre
       /(pss_w(i,j,tvar_)-pvs_0) ! diff de Pression Surface & Pr Vap sat 0 m


        q_0=               & ! humidite spec 0 metre (saturante)
        r_0/(1.+r_0)       & ! R_0 est le rapport de melange à 0 metre
       *0.98               ! (0.98 pour eau salée)

! trouver Q_0
!............> Fin.
! pour finalement obtenir:
      qstar_t(i,j,1)=cen*uv_10*( q10_t(i,j) - q_0)/ustar_t(i,j,1)

! Pour TETASTAR:
! même raisonnement en remplaçant CH par sa valeur à la neutralité c.a.d.
! CHN=1.2E-3 ! Smith (1980) instable
      tetastar_t(i,j,1)=                                                &
       chn*uv_10*( teta10_t(i,j)                                       & !Tem pot 10m
                 -(tem_t(i,j,kmax,tsst_)+273.15)*(1.e5/pss_w(i,j,tvar_))**.286 & !Tem pot  0m
                )/ustar_t(i,j,1)

!     IF(I.EQ.50.AND.J.EQ.50) THEN
!     WRITE(6,*)'Q_0=',Q_0
!     WRITE(6,*)'UV_10=',UV_10
!     WRITE(6,*)'CDN=',CDN
!     WRITE(6,*)'CEN=',CEN
!     WRITE(6,*)'CHN=',CHN
!     WRITE(6,*)'U*=',USTAR_Z(I,J,1)
!     WRITE(6,*)'Q*=',QSTAR_Z(I,J,1)
!     WRITE(6,*)'T*=',TETASTAR_Z(I,J,1)
!     STOP 'fifi'
!     WRITE(6,*)'T* m1=',TETASTAR_Z(I,J,1)
!     WRITE(6,*)'T* m2=',
!    1 CHN*UV_10*( TETA10_Z(I,J)                                    !Tem pot 10m
!    2           -(TEM_Z(I,J,NR-1)+273.15)                          !hyp : la Tem est pot à 0m
!    3           )/USTAR_Z(I,J,1)
!     WRITE(6,*)'SST=',TEM_Z(I,J,NR-1)+273.15
!     WRITE(6,*)'SSP=',(TEM_Z(I,J,NR-1)+273.15)
!    2           *(1.E5/PSS_Z(I,J,1))**.286
!     WRITE(6,*)'T10=',TETA10_Z(I,J)
!     WRITE(6,*)'T*=',TETASTAR_Z(I,J,1)
!     WRITE(6,*)'CHN*UV_10=',CHN*UV_10
!     WRITE(6,*)'DELTAT=',
!    1             TETA10_Z(I,J)                                    !Tem pot 10m
!    2           -(TEM_Z(I,J,NR-1)+273.15)*(1.E5/PSS_Z(I,J,1))**.286!Tem pot  0m
!     STOP 'fifi'
!     ENDIF

      enddo
      enddo

!*******************************************************************************
! INITIALISATION DE TABLEAUX ET DE CONSTANTES
! FIN.
!*******************************************************************************
      end subroutine initial_bulk

!.............................................................................

!______________________________________________________________________________
      subroutine iter_bulk(tsst_,tvar_)
!______________________________________________________________________________

      use module_principal
      implicit none
      integer tsst_,tvar_
#ifdef synopsis
       subroutinetitle='iter_bulk'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      const1=log(z_10/z_2)     ! log rapport des hauteurs 10m/2m
      const2=pi/2.             ! log rapport des hauteurs 10m/2m
      const3=rhoair*grav*z_2   ! dens air surface * G * 2m
      const4=10.               ! largeur de la jonction continue zl>0 zl<0 pour CHN
      loop2=6                  ! nbre max d'iterations permis

      do 1964 j=1,jmax
      do 1964 i=1,imax

      if(  mask_t(i,j,kmax+1).eq.1) then !%%%%%%%%%%%%%%%%%%%%%%%%%%%>

      uv_10=                                             & !vent à 10m instantanné
       max(un,                                           & !seuil pour eviter zero
       sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2) & !
          )                                              !à partir compo. u & v

      ro=                                      & ! densite de l'air à 2m
        (pss_w(i,j,tvar_)-const3)              & ! pres 2m (pres 0m + hypo hydros)
       *3.4897714e-3                           & ! *0.029/8.31
       /(teta2_t(i,j,tvar_)*(1.+0.608*q2_t(i,j,tvar_)))! /temp virt 2m

!                         .  .  .

! trouver Q_0
!............> début:

       pvs_0=                                      & ! Pres Vap Saturante 0 metre
       610.78*exp(17.27*tem_t(i,j,kmax,tsst_)          & !
                      /(tem_t(i,j,kmax,tsst_)+237.29)) ! Nb: 237.29=273.15-35.86

      r_0=                   & ! rapport de melange saturant à 0 metre
        0.622*pvs_0          & ! 0.622 * Pression Vapeur Saturante 0 metre
       /(pss_w(i,j,tvar_)-pvs_0) ! diff de Pression Surface & Pr Vap sat 0 m


        q_0=               & ! humidite spec 0 metre (saturante)
        r_0/(1.+r_0)       & ! R_0 est le rapport de melange à 0 metre
       *0.98               ! (0.98 pour eau salée)

! trouver Q_0
!............> Fin.

!                         .  .  .

! Début de la partie itérative:

      do 198 loop1=1,loop2 ! compteur d'iteration incrementé

!                         .  .  .

! Trouver le rapport z/L à 10m & 2m
!............> début:


! z/L= -g*Kappa*z*<w'Tu'>/(Tu*ustar**3) (eq 69)
! Tu = T*(1+0.608*Q) (eq17 page 96)
! Tu=bar(Tu)+ Tu'= (bar(T)+T')*(1+0.608*(bar(Q)+Q')
! <w'Tu'>=<w'*((bar(T)+T')*(1+0.608*(bar(Q)+Q'))-bar(Tu))>
! on néglige <w'T'Q'> et on obtient:
! Z/L= -g*Kappa*z*(tstar*(1+0.608*Q)                                   !29/07/06
!                 +qstar*0.608*T)/(T*(1+0.608*Q))/ustar**2

      x1=karman*grav*(                                                  &
                     tetastar_t(i,j,1)*(1.+0.608*   q10_t(i,j))         &
                    +   qstar_t(i,j,1)*    0.608*teta10_t(i,j)          &
                  )/(  teta10_t(i,j)*(1.+0.608*q10_t(i,j))              &
                       *ustar_t(i,j,1)**2          )
      zl_10=min(10.*x1,un*10.)
      zl_2= min( 2.*x1,un*10.)

! Trouver le rapport z/L à 10m & 2m
!............> Fin.

!                         .  .  .

! trouver les fonctions universelles psim, phim etc....
!............> début:

      if(zl_10.gt.0    ) then !§§§§§§§§§§§§§§§§§§§§§§§§§>
! demander à Claude pourquoi il y a 1.E-6 dans son programme

!     CHN=0.7E-3                                ! Smith (1980) stable

!........................................
!ccc   Plutôt que d'écrire:
!ccc   PHIM   =1.+FBETA*ZL_10
!ccc   PHIH   =PHIM
!ccc   PSIH_10=1.-PHIH
!ccc   on prefere l'expression equivalente et plus compacte:
       psih_10=-fbeta*zl_10

! mêmes remarques pour:
       psih_2= -fbeta*zl_2
       psim_10=-fbeta*zl_10
       psim_2= -fbeta*zl_2
!........................................
      else                    !§§§§§§§§§§§§§§§§§§§§§§§§§>

!     CHN=1.2E-3                                ! Smith (1980) instable

! calculs à 10m:
      phim=(1.-falpha*zl_10)**(-0.25) ! fonction universelle phim

      psim_10= 2.*log((1.+phim**(-1))/2.)                               &
                 +log((1.+phim**(-2))/2.)                               &
              -2.*atan(   phim**(-1)    )                               &
              +const2                        ! pi/2.

! Avant debug:
!     PSIH_10= 2.*LOG((1.+PHIM**(-0.5))/2.)
      psih_10= 2.*log((1.+phim**(-2.0))/2.)                            !02/06/06

!c est équivalent (mais plus rapide) à
!c    PHIQ=PHIM**0.5
!c    PSIH_10= 2.*LOG((1.+PHIQ**(-1))/2.)

! calculs à 2m:
      phim=(1.-falpha*zl_2)**(-0.25) ! fonction universelle phim

      psim_2=  2.*log((1.+phim**(-1))/2.)                               &
                 +log((1.+phim**(-2))/2.)                               &
              -2.*atan(   phim**(-1)    )                               &
              +const2                        ! pi/2.

! Avant debug:
!     PSIH_2= 2.*LOG((1.+PHIM**(-0.5))/2.)
      psih_2= 2.*log((1.+phim**(-2.0))/2.)                             !02/06/06

!c est équivalent (mais plus rapide) à
!c    PHIQ=PHIM**0.5
!c    PSIH_2= 2.*LOG((1.+PHIQ**(-1))/2.)


      endif                   !§§§§§§§§§§§§§§§§§§§§§§§§§>

! trouver les fonctions universelles psim, phim etc....
!............> Fin.

!                         .  .  .

! Calculs des Coef
! Début: ..........>

! l'operation de seuillage à 0.1E-3 n'est pas exactement equivalente
! à celle de Claude & Raoul. Leur demander pourquoi.
!c    CEN=1.2E-3    ! calcul placé en dehors des boucles ! Smith (1980)
!     CHN=0.7E-3  ! Smith (1980) stable   (ZL_10>0)
!     CHN=1.2E-3  ! Smith (1980) instable (ZL_10<0)
! Pour eviter les effets de seuil qui peuvent retarder la convergence
! on passe de l'un à l'autre cas via une transition douce (fonction TANH)
! dont la largeur est reglee par CONST4:
      chn=                                                              &
        0.7e-3*(1.+tanh(zl_10*const4))/2.                               &
       +1.2e-3*(1.-tanh(zl_10*const4))/2.

      cdn=max(un*0.93e-3,0.61e-3+6.3e-5*uv_10)              ! Smith (1980)

      ce=                                        & ! CE is the Dalton number.
         max(un*0.1e-3,                                                 &
         cen/( 1.- psim_10*cdn**0.5/karman       & ! (eq 106)
              -psih_10*cen*cdn**(-0.5)/karman    & ! en toute rigueur PSIQ_10
              +cen*psim_10*psih_10/karman**2   ) & ! devrait figurer à la place
            )                                    ! de PSIH_10 mais comme ils
                                                 ! sont egaux c'est equivalent

      ch=                                        & ! (eq 105)
         max(un*0.1e-3,                                                 &
         chn/( 1.- psim_10*cdn**0.5/karman                              &
              -psih_10*chn*cdn**(-0.5)/karman                           &
              +chn*psim_10*psih_10/karman**2)                           &
            )

      cd=                                                               &
         max(un*0.1e-3,                                                 &
         (cdn**(-0.5)-psim_10/karman)**(-2)       & ! (eq 103)
            )

! Calculs des Coef
! ..........> Fin.

!                         .  .  .

! mise à jour des variables "STAR" iteratives
! Début: ..........>

!        USTAR_Z(I,J,2)=(CD*UV_10**2)**0.5=SQRT(CD)*UV_10
         ustar_t(i,j,2)=sqrt(cd)*uv_10
         qstar_t(i,j,2)= ce*uv_10*( q10_t(i,j) - q_0)/ustar_t(i,j,2)
      tetastar_t(i,j,2)=                                                &
       ch*uv_10*(  teta10_t(i,j)                                       & !Tem pot 10m
                 -(tem_t(i,j,kmax,tsst_)+273.15)*(1.e5/pss_w(i,j,tvar_))**.286 & !Tem pot  0m
                )/ustar_t(i,j,2)

! suggestion tester la convergence si on divise par USTAR_Z(I,J,1 ou 2)
! mise à jour des variables "STAR" iteratives
! ..........> Fin.

!                         .  .  .

! trouver Q10_Z(I,J)
!............> début:

       q10_t(i,j)=                & ! humidite specifique à 10m
        q2_t(i,j,tvar_)               & ! humidite specifique à 2m
         +( qstar_t(i,j,2)/karman)& ! rapport echelle d'humidite / von karman
         *( const1                & ! log rapport des hauteurs 10m/2m
           -psih_10+psih_2 )      ! difference des fonctions psi 10m - 2m
                          ! NOTE: ! en toute rigueur c'est PSIQ_10-PSIQ_2
                                  ! qu'il faudrait prendre mais comme
                                  ! celles ci sont egales, on gagne du temps
                                  ! à ne pas les reecrire
! suggestion tester la convergence si on utilise QSTAR_Z(I,J,1 ou 2)

! trouver Q10_Z(I,J)
!............> Fin.

!                         .  .  .

! trouver TETA10_Z(I,J)
!............> début:

       teta10_t(i,j)=                   & ! temp potentielle 10m
        teta2_t(i,j,tvar_)                  & ! temp potentielle 2m
            +( tetastar_t(i,j,2)/karman)& ! rap echelle tem pot / karman
            *( const1                   & ! log rapport des hauteurs 10m/2m
              -psih_10+psih_2 )         ! difference des fonctions psi 10m - 2m


! trouver TETA10_Z(I,J)
!............> Fin.

!     IF(I.EQ.50.AND.J.EQ.50) THEN
!     WRITE(6,*)
!     WRITE(6,*)'LOOP1 LOOP2 = ',LOOP1,LOOP2
!     WRITE(6,*)'ZL_10  =',ZL_10
!     WRITE(6,*)'PSIH_10=',PSIH_10
!     WRITE(6,*)'PSIH_2 =',PSIH_2
!     WRITE(6,*)'PSIM_10=',PSIM_10
!     WRITE(6,*)'PSIM_2 =',PSIM_2
!     WRITE(6,*)'CD  =',CD
!     WRITE(6,*)'CE  =',CE
!     WRITE(6,*)'CH  =',CH
!     WRITE(6,*)'CDN =',CDN
!     WRITE(6,*)'CEN =',CEN
!     WRITE(6,*)'CHN =',CHN
!     WRITE(6,*)'poid1=',(1.+TANH(ZL_10*CONST4))*0.5
!     WRITE(6,*)'poid2=',(1.-TANH(ZL_10*CONST4))*0.5
!     WRITE(6,*)'U*=',USTAR_Z(I,J,2)
!     WRITE(6,*)'Q*=',QSTAR_Z(I,J,2)
!     WRITE(6,*)'T*=',TETASTAR_Z(I,J,2)
!     WRITE(6,*)'Q10=',Q10_Z(I,J)
!     WRITE(6,*)'T10=',TETA10_Z(I,J)
!     ENDIF

!                         .  .  .

! Test de convergence et
! mise à jour des variables itérées:
!............> début:
      if(    abs(   ustar_t(i,j,1)-   ustar_t(i,j,2)).gt.1.e-4          &
         .or.abs(   qstar_t(i,j,1)-   qstar_t(i,j,2)).gt.1.e-4          &
         .or.abs(tetastar_t(i,j,1)-tetastar_t(i,j,2)).gt.1.e-4          &
        ) then !---------------------------------------------------->

         ustar_t(i,j,1)=ustar_t(i,j,2)
         qstar_t(i,j,1)=qstar_t(i,j,2)
      tetastar_t(i,j,1)=tetastar_t(i,j,2)

       else    !---------------------------------------------------->

       goto 199 ! sortie de la boucle iterative si test OK

       endif   !---------------------------------------------------->


! Test de convergence et
! mise à jour des variables itérées:
!............> Fin.

  198 continue  ! fin du forfait max pour la boucle iterative
  199 continue  ! sortie de boucle iterative si test convergence OK

!c    X10=X10+1.
!c    X11=X11+REAL(LOOP1)
!c    XY_Z(I,J,1)=REAL(LOOP1)

! Nous sommes sortis de la boucle iterative et nous calculons les flux:

      slhf_w(i,j,1)=                          & ! flux de chaleur latente
        ro                                    & ! densité de l'air
       *(2500.9-2.365*tem_t(i,j,kmax,tsst_))*1.e3 & ! Chaleur latente de vaporisation
       *ce                                    & ! CE is the Dalton number
       *uv_10                                 & ! vitesse du vent à 10 mètres
       *(q10_t(i,j)   -q_0)                   ! Delta d'humidité spécifique

      sshf_w(i,j,1)=                          & ! flux de chaleur sensible
        ro                                    & ! densité de l'air
       *cp_air                                & ! capacité calorifique de l'air
       *ch                                    & ! CH is the Stanton number
       *uv_10                                 & ! vitesse du vent à 10 mètres
       *(  teta10_t(i,j)                      & ! temp potentielle à 10m
         -(tem_t(i,j,kmax,tsst_)+273.15)          & ! temp potentielle à 0m
         *(1.e5/pss_w(i,j,tvar_))**.286 )         ! temp potentielle à 0m suite.


!     if (iairsea.eq.2)                                              & !060104
!       snsf_w(i,j,1)=                       & !  flux IR net
!       snsf_w(i,j,1)                        & ! =flux IR descendant
!      -stefan*(273.15+tem_t(i,j,kmax,tsst_))**4 ! +flux IR montant (loi de stefan)

!     if(iairsea.eq.3) then                  ! calcul Firnet par formule May reprise Pinardi 060104
!      x1=pss_w(i,j,tvar_)*1.608*q2_t(i,j,tvar_)/    & ! Pression de vapeur fonction de humidite specifique
!       (1.+0.608*q2_t(i,j,tvar_))/              & ! Queney page 115               060104
!       1.e2                                 ! exprimee en millibars         060104
!      snsf_w(i,j,1)=-(                      & ! flux IR net                   060104
!       (1.-0.75*(snsf_w(i,j,1)**3.4))*      & ! SNSF : ici nebulosite         060104
!       (stefan*teta2_t(i,j,tvar_)**4*           & ! Temperature de l'air (on approxime temperature etat e
!       (0.4-0.05*sqrt(x1))+                                            & !060104
!       4.*stefan*(teta2_t(i,j,tvar_)**3)*                                  & !060104
!       (273.15+tem_t(i,j,kmax,tsst_)-teta2_t(i,j,tvar_))))                     !060104
!     endif                                                             !060104

      wstress_w(i,j)=                        & ! module de la tension du vent
        ro                                   & ! densité de l'air
       *ustar_t(i,j,2)**2                    ! echelle de vitesse


      if(loop1.eq.loop2+1) then !:::::::::::::::::::::::::::::::::::::::::::::>

! Si on se trouve ici cela signifie que l'algo n'a pas réussi à converger
! correctement dans le forfait imparti. On reset les variables iteratives
! aux valeurs prescrites à la neutralité (Gradients verticaux nuls +
! coef CD, CE & CH remplacés par CDN, CEN & CHN):

            teta10_t(i,j)=teta2_t(i,j,tvar_)
               q10_t(i,j)=   q2_t(i,j,tvar_)
           ustar_t(i,j,1)=sqrt(cdn)*uv_10
           qstar_t(i,j,1)=cen*uv_10*( q10_t(i,j) - q_0)/ustar_t(i,j,1)
        tetastar_t(i,j,1)=                                             &
       chn*uv_10*( teta10_t(i,j)                                       & !Tem pot 10m
                 -(tem_t(i,j,kmax,tsst_)+273.15)*(1.e5/pss_w(i,j,tvar_))**.286 & !Tem pot  0m
                )/ustar_t(i,j,1)

      endif                     !:::::::::::::::::::::::::::::::::::::::::::::>

!     IF(I.EQ.50.AND.J.EQ.50)THEN
!     WRITE(6,*)
!     WRITE(6,*)'SLHF=',SLHF_Z(I,J,1)
!     WRITE(6,*)'SSHF=',SSHF_Z(I,J,1)
!     WRITE(6,*)'SNSF=',SNSF_Z(I,J,1)
!     WRITE(6,*)'WSTR=',WSTRESS_Z(I,J)
!     STOP 'coco'
!     ENDIF
!     WRITE(67,*)SLHF_Z(I,J,1),SSHF_Z(I,J,1)
!    &          ,SNSF_Z(I,J,1),WSTRESS_Z(I,J)

      endif                          !%%%%%%%%%%%%%%%%%%%%%%%%%%%>
 1964 continue

!cccccSTOP 'groscoco'
      end subroutine iter_bulk
#endif
!______________________________________________________________________________
      subroutine ztoxy_bulk
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='ztoxy_bulk'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(flag_wstressbulk==1) then !pmx> 25-01-18

! Repartir l'info de Z vers X et Y:


      do 1968 j=2,jmax
      do 1968 i=2,imax

         wstress_u(i  ,j,1)=0.5*(                                     &

         wstress_w(i  ,j)*uwind_t(i  ,j,1)                            &
         /max(small1,sqrt(uwind_t(i  ,j,1)**2+vwind_t(i  ,j,1)**2))   &

        +wstress_w(i-1,j)*uwind_t(i-1,j,1)                            &
         /max(small1,sqrt(uwind_t(i-1,j,1)**2+vwind_t(i-1,j,1)**2))   &

                                )

         wstress_v(i,j  ,1)=0.5*(                                     &

         wstress_w(i,j  )*vwind_t(i,j  ,1)                            &
         /max(small1,sqrt(uwind_t(i,j  ,1)**2+vwind_t(i,j  ,1)**2))   &

        +wstress_w(i,j-1)*vwind_t(i,j-1,1)                            &
         /max(small1,sqrt(uwind_t(i,j-1,1)**2+vwind_t(i,j-1,1)**2))   &

                                )
! claude
          wstress_u(i  ,j,1)= wstress_u(i  ,j,1)*1.1
          wstress_v(i  ,j,1)= wstress_v(i  ,j,1)*1.1
 1968 continue

      endif                        !pmx> 25-01-18

      if(flag_wstressbulk==0) then !pmx> 25-01-18

! Si le stress est donnE par le modele meteo retrouver wstress_w (qui est modifiE par
! les routines bulk) A partir de wstress_u et wstress_v (donnEs par modele meteo)
       do j=1,jmax ; do i=1,imax
        wstress_w(i,j)=sqrt( &
                             (0.5*(wstress_u(i,j,1)+wstress_u(i+1,j,1)))**2 &
                            +(0.5*(wstress_v(i,j,1)+wstress_v(i,j+1,1)))**2 &
                           )
       enddo       ; enddo

      endif                        !pmx> 25-01-18

      end subroutine ztoxy_bulk

!.............................................................................

      subroutine bukl_relative_wind(tsst_,tvar_,relativewind_) !23-11-15
      use module_principal
      implicit none
      integer  ichoix,tsst_,tvar_
      real*4 relativewind_
#ifdef synopsis
       subroutinetitle='bukl_relative_wind'
       subroutinedescription= &
       'removes a fraction of the surface current to the wind velocity'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Le vent relatif n'est applique que dans le cas standard c.a.d. tvar_=1
! Le cas tvar_=2 correspond a l'estimation du flux de qdm calcule dans
! le modele meteo qui suppose un courant de surface nul, enfin on espere :)
! D'autre part ce serait un bug que de l'oter pour tvar_=2 car cela
! reviendrait a enlever le courant 2 fois (au coup suivant il serait de
! nouveau enlever a tvar_=1)

       if(tvar_==1) then !>>>>>>>>
        k=kmax
        do j=0,jmax+1 ; do i=0,imax+1
         uwind_t(i,j,tvar_)=uwind_t(i,j,tvar_) &
           -relativewind_*0.5*(vel_u(i,j,k,before)+vel_u(i+1,j,k,before))
         vwind_t(i,j,tvar_)=vwind_t(i,j,tvar_) &
           -relativewind_*0.5*(vel_v(i,j,k,before)+vel_v(i,j+1,k,before))
        enddo       ; enddo
       endif             !>>>>>>>>

      end subroutine bukl_relative_wind

!.............................................................................

      subroutine bulk_stress_w2uv(tvar_) ! 24-11-15
      use module_principal
      use module_parallele
      implicit none
      integer tvar_
#ifdef synopsis
       subroutinetitle='bulk_stress_w2uv'
       subroutinedescription= &
       'computes wstress_u wstress_v from wstress_w'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Dans la cas de la couche limite simplifiee le calcul des anomalies
! de vent aux points "t" sur 1:imax et 1:jmax necessite la connaissance
! de wstress_u et wstress_v sur les bords.

      do j=1,jmax
      do i=1,imax

       xy_t(i,j,0)=0.5*wstress_w(i,j) &
           /(small1+sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2))

      enddo
      enddo
! Conditions limites z1 & continuite mpi
      if(obcstatus(ieq1)==1)   xy_t(0     ,:     ,0)=xy_t(1   ,:   ,0)
      if(obcstatus(ieqimax)==1)xy_t(imax+1,:     ,0)=xy_t(imax,:   ,0)
      if(obcstatus(jeq1)==1)   xy_t(:     ,0     ,0)=xy_t(:   ,1   ,0)
      if(obcstatus(jeqjmax)==1)xy_t(:     ,jmax+1,0)=xy_t(:   ,jmax,0)
      call obc_ext_xy_t('z1',0)

      do j=1,jmax  
      do i=1,imax+1
         wstress_u(i,j,1)=xy_t(i  ,j,0)*uwind_t(i  ,j,tvar_) &
                         +xy_t(i-1,j,0)*uwind_t(i-1,j,tvar_)  
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
         wstress_v(i,j,1)=xy_t(i,j  ,0)*vwind_t(i,j  ,tvar_) &
                         +xy_t(i,j-1,0)*vwind_t(i,j-1,tvar_)  
      enddo
      enddo

      end subroutine bulk_stress_w2uv 
!......................................................................
      subroutine bulk_1dv ! 27-11-15
      use module_principal
      implicit none

          slhf_w(:,:,1)=   slhf_w(2,2,1) 
          snsf_w(:,:,1)=   snsf_w(2,2,1) 
          sshf_w(:,:,1)=   sshf_w(2,2,1) 
       wstress_w(:,:  )=wstress_w(2,2  )
       wstress_u(:,:,1)=wstress_u(2,2,1)
       wstress_v(:,:,1)=wstress_v(2,2,1)

      end subroutine bulk_1dv ! 27-11-15

!.........................................................................

      subroutine initial_bulk_2006(tsst_,tvar_)
      use module_principal
      implicit none
      integer tsst_,tvar_
#ifdef synopsis
       subroutinetitle='initial_bulk_2006'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Proposé par Guillaume Reffray: greffray@mercator-ocean.fr
!        08/06/06: rationnalisation du calcul (via approximations legitimes)
!                  proposée par Guillaume Reffray dans le cadre du projet
!                  "OPA Symphonique" conduit par MERCATOR & POC.
! Document de reference:
! Diurnal to decadal global forcing for oceran and sea-ice models:
! The data sets and flux climatologies. Auteurs: W.G. Large & S.G. Yeager
! http://www.cgd.ucar.edu/oce/pubs/04pubs_files/TN460.pdf



!...............................................................................
! Constantes intervenant dans le calcul des flux air/mer
! calculés par les formules "bulk"
      karman=0.4            ! Constante de Von Karman
      stefan=5.67032e-8     ! constante de stephan (W/m**2/K**4)
      z_2=2.                ! niveau 2  mètres
      z_10=10.              ! niveau 10 mètres
      falpha=16.            ! coef fonction universelle instable (1-ALPHA Z/L)
      fbeta=7.              ! coef fonction universelle   stable (1+BETA  Z/L)
      cp_air=1000.5         ! capacité calorifique de l'air
      lv=2.5e6              ! chaleur latente de vaporisation (J/kg)
!...............................................................................

! premiere estimation: on se place à la stabilité:

      do j=1,jmax
      do i=1,imax

      uv_10=                                             & !vent à 10m instantanné
       max(un,                                           & !seuil pour eviter zero
       sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2)         & !
          )                                              !à partir compo. u & v

      teta10_t(i,j)=teta2_t(i,j,tvar_)
         q10_t(i,j)=   q2_t(i,j,tvar_)

! Pour USTAR:
! à la neutralité on remplace CD de la formule suivante:
! USTAR_Z(I,J,1)=(CD*UV_10**2)**0.5=SQRT(CD)*UV_10
! par CDN, sachant que CDN=MAX(0.93E-3,0.61E-3+6.3E-5*UV_10) ! Smith (1980)
! pour finalement obtenir:

      cdn  = 1e-3 * ( 2.7/uv_10 + 0.142 + uv_10/13.09 )       !  \\ L & Y eq. (6a)
      cen  = 1e-3 * ( 34.6 * sqrt(cdn) )                      !  \\ L & Y eq. (6b)

                                                              !  \\ L & Y eq. (6c), (6d)
! X1=Delta de temperature potentielle entre les niveaux 10m et 0m:
      x1=teta10_t(i,j)                                          & ! T pot 10m
          -(tem_t(i,j,kmax,tsst_)+273.15)*(1.e5/pss_w(i,j,tvar_))**.286 ! T pot  0m

      if(x1.gt.0.0) then          !---->
! Cas stable:
       chn  = 18.0*1e-3*sqrt(cdn)
      else                        !---->
! Cas instable:
       chn =  32.7*1e-3*sqrt(cdn)
      endif                       !---->


      !! Initializing transfert coefficients with their first guess neutral equivalents :

      ustar_t(i,j,1)=sqrt(cdn)*uv_10

! Pour QSTAR:
! à la neutralité on remplace CE de la formule suivante (Geernaert p111):
! QSTAR_Z(I,J,1)=(CE*UV_10*(Q10_Z(I,J)-Q_0))/USTAR_Z(I,J,1)
! par CEN, sachant que Q_0 est donné par:
! trouver Q_0
!............> début:

      airdensity=pss_w(i,j,1)/(287.058*teta2_t(i,j,tvar_)) !25-04-16

      q_0 = 0.98*640380/airdensity*exp(-5107.4/(tem_t(i,j,kmax,tsst_)+273.15))

! trouver Q_0
!............> Fin.
! pour finalement obtenir:
      qstar_t(i,j,1)=cen*uv_10*( q10_t(i,j) - q_0)/ustar_t(i,j,1)

! Pour TETASTAR:
! même raisonnement en remplaçant CH par sa valeur à la neutralité c.a.d.
! CHN=1.2E-3 ! Smith (1980) instable
      tetastar_t(i,j,1)=chn*uv_10*x1/ustar_t(i,j,1) ! X1=delta de T pot entre 10m et 0m

!     IF(I.EQ.50.AND.J.EQ.50) THEN
!     WRITE(6,*)'Q_0=',Q_0
!     WRITE(6,*)'UV_10=',UV_10
!     WRITE(6,*)'CDN=',CDN
!     WRITE(6,*)'CEN=',CEN
!     WRITE(6,*)'CHN=',CHN
!     WRITE(6,*)'U*=',USTAR_Z(I,J,1)
!     WRITE(6,*)'Q*=',QSTAR_Z(I,J,1)
!     WRITE(6,*)'T*=',TETASTAR_Z(I,J,1)
!     WRITE(6,*)'T10M=',TETA10_Z(I,J)
!     WRITE(6,*)'SST=',TEM_Z(I,J,NR-1)+273.15
!     WRITE(6,*)'T*=',TETASTAR_Z(I,J,1)
!     WRITE(6,*)'CHN*UV_10=',CHN*UV_10
!     WRITE(6,*)'DELTAT=',X1
!     STOP 'fifi2006'
!     ENDIF

      enddo
      enddo

      end subroutine initial_bulk_2006

!............................................................................

      subroutine iter_bulk_2006(tsst_,tvar_)
      use module_principal
      use module_parallele
      implicit none
      integer tsst_,tvar_
#ifdef synopsis
       subroutinetitle='iter_bulk_2006'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Proposé par Guillaume Reffray: greffray@mercator-ocean.fr
!        08/06/06: rationnalisation du calcul (via approximations legitimes)
!                  proposée par Guillaume Reffray dans le cadre du projet
!                  "OPA Symphonique" conduit par MERCATOR & POC.
! Document de reference:
! Diurnal to decadal global forcing for oceran and sea-ice models:
! The data sets and flux climatologies. Auteurs: W.G. Large & S.G. Yeager
! http://www.cgd.ucar.edu/oce/pubs/04pubs_files/TN460.pdf

      const1=log(z_10/z_2)       ! log rapport des hauteurs 10m/2m
      const2=pi/2.               ! log rapport des hauteurs 10m/2m
      const4=10.                 ! largeur de la jonction continue zl>0 zl<0 pour CHN
      const6=0.1e-3              ! Seuil mini pour coef

! Nombre d'iterations maximum autorisE
      if(iteration3d==0) then !m°v°m> !18-09-21
       loop2=50 ! Un grand nombre A l'etat initial pour atteindre la convergence
      else                    !m°v°m>
       loop2=6  ! un petit nombre est suffisant pour converger ensuite
      endif                   !m°v°m> 

      do 1964 j=1,jmax
      do 1964 i=1,imax

      if(  mask_t(i,j,kmax+1).eq.1) then !%%%%%%%%%%%%%%%%%%%%%%%%%%%>

            ustar_bef=   ustar_t(i,j,1)
            qstar_bef=   qstar_t(i,j,1)
         tetastar_bef=tetastar_t(i,j,1)

      airdensity=pss_w(i,j,1)/(287.058*teta2_t(i,j,tvar_)) !25-04-16

! ajout Patrick 08/06/06                                 !08/06/06
      sst_kelvin=tem_t(i,j,kmax,tsst_)+273.15                        ! temp      0m
      sst1000hpa_kelvin=sst_kelvin*(1.e5/pss_w(i,j,tvar_))**.286    ! temp pot  0m

      uv_10=                                             & !vent à 10m instantanné
       max(un,                                           & !seuil pour eviter zero
       sqrt(uwind_t(i,j,tvar_)**2+vwind_t(i,j,tvar_)**2)         & !
          )                                              !à partir compo. u & v

! trouver Q_0
!............> début:

!cccccQ_0 = 0.98*640380/RHOAIR*EXP(-5107.4/(TEM_Z(I,J,NR-1)+273.15))
      q_0 =0.98*640380./airdensity*exp(-5107.4/sst_kelvin)  !08/06/06 !24-04-16 !25-04-16

! trouver Q_0
!............> Fin.

!                         .  .  .

! Début de la partie itérative:

      do 198 loop1=1,loop2 ! compteur d'iteration incrementé

! Trouver le rapport z/L à 10m & 2m
!............> début:

      x1=karman*grav*(                                                  &
                     tetastar_t(i,j,1)*(1.+0.608*   q10_t(i,j))         &
                    +   qstar_t(i,j,1)*    0.608*teta10_t(i,j)          &
                  )/(  teta10_t(i,j  )*(1.+0.608*   q10_t(i,j))         &
                       *ustar_t(i,j,1)**2          )
      zl_10=min(10.*x1,10.*un)
      zl_2= min( 2.*x1,10.*un)

! Trouver le rapport z/L à 10m & 2m
!............> Fin.

!                         .  .  .

! trouver les fonctions universelles psim, phim etc....
!............> début:

      if(zl_10.gt.0    ) then !§§§§§§§§§§§§§§§§§§§§§§§§§>

!........................................
!ccc   Plutôt que d'écrire:
!ccc   PHIM   =1.+FBETA*ZL_10
!ccc   PHIH   =PHIM
!ccc   PSIH_10=1.-PHIH
!ccc   on prefere l'expression equivalente et plus compacte:
       psih_10=-fbeta*zl_10

! mêmes remarques pour:
       psih_2= -fbeta*zl_2
       psim_10=-fbeta*zl_10
       psim_2= -fbeta*zl_2

!ccc   CHN  = 18.*1E-3*SQRT(CDN)               ! \\ L & Y eq. (6c), (6d)

!........................................
      else                    !§§§§§§§§§§§§§§§§§§§§§§§§§>

! calculs à 10m:
      phim=(1.-falpha*zl_10)**(-0.25) ! fonction universelle phim

      psim_10= 2.*log((1.+phim**(-1))/2.)                               &
                 +log((1.+phim**(-2))/2.)                               &
              -2.*atan(   phim**(-1)    )                               &
              +const2                        ! pi/2.

      psih_10= 2.*log((1.+phim**(-2))/2.)
!c est équivalent (mais plus rapide) à
!c    PHIQ=PHIM**0.5
!c    PSIH_10= 2.*LOG((1.+PHIQ**(-1))/2.)

! calculs à 2m:
      phim=(1.-falpha*zl_2)**(-0.25) ! fonction universelle phim

      psim_2=  2.*log((1.+phim**(-1))/2.)                               &
                 +log((1.+phim**(-2))/2.)                               &
              -2.*atan(   phim**(-1)    )                               &
              +const2                        ! pi/2.

      psih_2= 2.*log((1.+phim**(-2.))/2.)
!c est équivalent (mais plus rapide) à
!c    PHIQ=PHIM**0.5
!c    PSIH_2= 2.*LOG((1.+PHIQ**(-1))/2.)

!ccc  CHN  = 32.7*1E-3*SQRT(CDN)               ! \\ L & Y eq. (6c), (6d)

      endif                   !§§§§§§§§§§§§§§§§§§§§§§§§§>

! trouver les fonctions universelles psim, phim etc....
!............> Fin.
!                         .  .  .
! Calculs des Coef
! Début: ..........>

! l'operation de seuillage à 0.1E-3 n'est pas exactement equivalente
! à celle de Claude & Raoul. Leur demander pourquoi.
!c    CEN=1.2E-3    ! calcul placé en dehors des boucles ! Smith (1980)
!     CHN=0.7E-3  ! Smith (1980) stable   (ZL_10>0)
!     CHN=1.2E-3  ! Smith (1980) instable (ZL_10<0)
! Pour eviter les effets de seuil qui peuvent retarder la convergence
! on passe de l'un à l'autre cas via une transition douce (fonction TANH)
! dont la largeur est réglée par CONST4:

      cdn  = 1e-3 * ( 2.7/uv_10 + 0.142 + uv_10/13.09 )       !  \\ L & Y eq. (6a)
      cen  = 1e-3 * ( 34.6 * sqrt(cdn) )                      !  \\ L & Y eq. (6b)
      chn  = 0.5*1.e-3*sqrt(cdn)*(                            & !  \\ L & Y combine:
        18.0*(1.+tanh(zl_10*const4))                          & !   eq. 6c et
       +32.7*(1.-tanh(zl_10*const4))  )                       !   eq. 6d

!...........................................
! méthode Geernaert:
!     CD=
!    1   MAX(CONST6,
!    2   (CDN**(-0.5)-PSIM_10/KARMAN)**(-2)
!    3      )

!     CH=
!    1   MAX(CONST6,
!    2   CHN/( 1.- PSIM_10*CDN**0.5/KARMAN
!    3        -PSIH_10*CHN*CDN**(-0.5)/KARMAN
!    4        +CHN*PSIM_10*PSIH_10/KARMAN**2)
!    5      )

!     CE=                                        ! CE is the Dalton number.
!    1   MAX(CONST6,
!    2   CEN/( 1.- PSIM_10*CDN**0.5/KARMAN
!    3        -PSIH_10*CEN*CDN**(-0.5)/KARMAN    ! en toute rigueur PSIQ_10
!    4        +CEN*PSIM_10*PSIH_10/KARMAN**2   ) ! devrait figurer à la place
!    5      )                                    ! de PSIH_10 mais comme ils
!...........................................
! méthode L & Y (equivalente à la precedente)
      cd=max(const6,                                                    &
        cdn*(1-sqrt(cdn)*karman**(-1)*psim_10)**(-2) & ! L & Y (10a) pour zu=10m
            )

      ch=max(const6,                                                    &
        chn*sqrt(cd/cdn)*(1-chn*karman**(-1)*cdn**(-0.5)*psih_10)**(-1) &
            )

      ce=max(const6,                                                    &
        cen*sqrt(cd/cdn)*(1-cen*karman**(-1)*cdn**(-0.5)*psih_10)**(-1) &
            )
!...........................................

! Calculs des Coef
! ..........> Fin.

!                         .  .  .

! mise à jour des variables "STAR" iteratives
! Début: ..........>

! Appliquer un filtre EMA sur qstar, tetastar (mais pas sur ustar) q10_t et teta10_t
! accelere la convergence  pour les cas de vents faibles !18-09-21
         ustar_t(i,j,1)=sqrt(cd)*uv_10   
         qstar_t(i,j,1)= ce*uv_10*( q10_t(i,j) - q_0)/ustar_t(i,j,1) & 
                        *0.5+0.5*qstar_t(i,j,1)    ! filtre EMA !18-09-21
      tetastar_t(i,j,1)= ch*uv_10*( teta10_t(i,j)-sst1000hpa_kelvin)/ustar_t(i,j,1) & 
                        *0.5+0.5*tetastar_t(i,j,1) ! filtre EMA !18-09-21

! suggestion tester la convergence si on divise par USTAR_Z(I,J,1 ou 2)
! mise à jour des variables "STAR" iteratives
! ..........> Fin.

!                         .  .  .

! trouver Q10_Z(I,J)
!............> début:

       q10_t(i,j)=(                  & ! humidite specifique à 10m
        q2_t(i,j,tvar_)              & ! humidite specifique à 2m
         +( qstar_t(i,j,1)/karman)   & ! rapport echelle d'humidite / von karman
         *( const1                   & ! log rapport des hauteurs 10m/2m
           -psih_10+psih_2 )         & ! difference des fonctions psi 10m - 2m
                  )*0.5+0.5*q10_t(i,j) ! filtre EMA !18-09-21
                          ! NOTE: ! en toute rigueur c'est PSIQ_10-PSIQ_2
                                  ! qu'il faudrait prendre mais comme
                                  ! celles ci sont egales, on gagne du temps
                                  ! à ne pas les reecrire
! suggestion tester la convergence si on utilise QSTAR_Z(I,J,1 ou 2)

! trouver Q10_Z(I,J)
!............> Fin.

!                         .  .  .

! trouver TETA10_Z(I,J)
!............> début:

       teta10_t(i,j)=(                      & ! temp potentielle 10m
        teta2_t(i,j,tvar_)                  & ! temp potentielle 2m
            +( tetastar_t(i,j,1)/karman)    & ! rap echelle tem pot / karman
            *( const1                       & ! log rapport des hauteurs 10m/2m
              -psih_10+psih_2 )             & ! difference des fonctions psi 10m - 2m
                     )*0.5+0.5*teta10_t(i,j)  ! filtre EMA !18-09-21


! trouver TETA10_Z(I,J)
!............> Fin.

!     IF(I.EQ.50.AND.J.EQ.50) THEN
!     WRITE(6,*)
!     WRITE(6,*)'LOOP1 LOOP2 = ',LOOP1,LOOP2
!     WRITE(6,*)'ZL_10  =',ZL_10
!     WRITE(6,*)'PSIH_10=',PSIH_10
!     WRITE(6,*)'PSIH_2 =',PSIH_2
!     WRITE(6,*)'PSIM_10=',PSIM_10
!     WRITE(6,*)'PSIM_2 =',PSIM_2
!     WRITE(6,*)'CD  =',CD
!     WRITE(6,*)'CE  =',CE
!     WRITE(6,*)'CH  =',CH
!     WRITE(6,*)'CDN =',CDN
!     WRITE(6,*)'CEN =',CEN
!     WRITE(6,*)'CHN =',CHN
!     WRITE(6,*)'poid1=',(1.+TANH(ZL_10*CONST4))*0.5
!     WRITE(6,*)'poid2=',(1.-TANH(ZL_10*CONST4))*0.5
!     WRITE(6,*)'U*=',USTAR_Z(I,J,2)
!     WRITE(6,*)'Q*=',QSTAR_Z(I,J,2)
!     WRITE(6,*)'T*=',TETASTAR_Z(I,J,2)
!     WRITE(6,*)'Q10=',Q10_Z(I,J)
!     WRITE(6,*)'T10=',TETA10_Z(I,J)
!     STOP 'toto2006'
!     ENDIF

!                         .  .  .

! Test de convergence et
! mise à jour des variables itérées:
!............> début:
      if(    abs(   ustar_bef-   ustar_t(i,j,1)).gt.1.e-4          &
         .or.abs(   qstar_bef-   qstar_t(i,j,1)).gt.1.e-4          &
         .or.abs(tetastar_bef-tetastar_t(i,j,1)).gt.1.e-4          &
        ) then !---------------------------------------------------->

            ustar_bef=   ustar_t(i,j,1)
            qstar_bef=   qstar_t(i,j,1)
         tetastar_bef=tetastar_t(i,j,1)

       else    !---------------------------------------------------->

       goto 199 ! sortie de la boucle iterative si test OK

       endif   !---------------------------------------------------->


! Test de convergence et
! mise à jour des variables itérées:
!............> Fin.

  198 continue  ! fin du forfait max pour la boucle iterative
  199 continue  ! sortie de boucle iterative si test convergence OK


! Nous sommes sortis de la boucle iterative et nous calculons les flux:

      slhf_w(i,j,1)=                         & ! flux de chaleur latente
        airdensity                           & ! densité de l'air
       *lv                                   & ! Chaleur latente de vaporisation
       *ce                                   & ! CE is the Dalton number
       *uv_10                                & ! vitesse du vent à 10 mètres
       *(q10_t(i,j)   -q_0)                    ! Delta d'humidité spécifique

      sshf_w(i,j,1)=                         & ! flux de chaleur sensible
        airdensity                           & ! densité de l'air
       *cp_air                               & ! capacité calorifique de l'air
       *ch                                   & ! CH is the Stanton number
       *uv_10                                & ! vitesse du vent à 10 mètres
       *(teta10_t(i,j)-sst1000hpa_kelvin)      ! delta temp potentielle 10m-0m


!     if (iairsea.eq.2)                      & !060104
!       snsf_w(i,j,1)=                       & !  flux IR net
!       snsf_w(i,j,1)                        & ! =flux IR descendant
!      -stefan*sst_kelvin**4                         ! +flux IR montant (loi de stefan)

!     if(iairsea.eq.3) then                  ! calcul Firnet par formule May reprise Pinardi 060104
!      x1=pss_w(i,j,tvar_)*1.608*q2_t(i,j,tvar_)/    & ! Pression de vapeur fonction de humidite specifique
!       (1.+0.608*q2_t(i,j,tvar_))/                  & ! Queney page 115               060104
!       1.e2                                           ! exprimee en millibars         060104
!      snsf_w(i,j,1)=-(                              & ! flux IR net                   060104
!       (1.-0.75*(snsf_w(i,j,1)**3.4))*              & ! SNSF : ici nebulosite         060104
!       (stefan*teta2_t(i,j,tvar_)**4*               & ! Temperature de l'air (on approxime temperature etat e
!       (0.4-0.05*sqrt(x1))+                                          & !060104
!       4.*stefan*(teta2_t(i,j,tvar_)**3)*                            & !060104
!       (273.15+tem_t(i,j,kmax,tsst_)-teta2_t(i,j,tvar_))))             !060104
!     endif                                                             !060104

      wstress_w(i,j)=                        & ! module de la tension du vent
        airdensity                           & ! densité de l'air
       *ustar_t(i,j,1)**2                      ! echelle de vitesse


! ATTENTION CES LIGNES SONT IMPORTANTES: 
! IL A ETE VERIFIE QUE LE SCHEMA RISQUE D'ETRE INSTABLE SANS ELLES !18-09-21
      if(loop1.eq.loop2+1) then !:::::::::::::::::::::::::::::::::::::::::::::>

! Si on se trouve ici cela signifie que l'algo n'a pas réussi à converger
! correctement dans le forfait imparti. On reset les variables iteratives
! aux valeurs prescrites à la neutralité (Gradients verticaux nuls +
! coef CD, CE & CH remplacés par CDN, CEN & CHN):

            teta10_t(i,j)=teta2_t(i,j,tvar_)
               q10_t(i,j)=   q2_t(i,j,tvar_)
           ustar_t(i,j,1)=sqrt(cdn)*uv_10
           qstar_t(i,j,1)=cen*uv_10*( q10_t(i,j) - q_0)/ustar_t(i,j,1)
        tetastar_t(i,j,1)=chn*uv_10*(teta10_t(i,j)-sst1000hpa_kelvin )/ustar_t(i,j,1)

      endif                     !:::::::::::::::::::::::::::::::::::::::::::::>


!     IF(I.EQ.50.AND.J.EQ.50)THEN
!     WRITE(6,*)
!     WRITE(6,*)'SLHF=',SLHF_Z(I,J,1)
!     WRITE(6,*)'SSHF=',SSHF_Z(I,J,1)
!     WRITE(6,*)'SNSF=',SNSF_Z(I,J,1)
!     WRITE(6,*)'WSTR=',WSTRESS_Z(I,J)
!     STOP'coco2006'
!     ENDIF
!     WRITE(66,*)SLHF_Z(I,J,1),SSHF_Z(I,J,1)
!    &          ,SNSF_Z(I,J,1),WSTRESS_Z(I,J)

      endif                          !%%%%%%%%%%%%%%%%%%%%%%%%%%%>
     if(i+par%timax(1)==10.and.j+par%tjmax(1)==10)then
      write(50,*)'sshf_w(i,j,1),slhf_w(i,j,1),wstress_w(i,j)', &
       sshf_w(i,j,1),slhf_w(i,j,1),wstress_w(i,j)
!      write(60,*)PTA,PPS,PVMOD,PSST,PQA
     endif
 1964 continue

      end subroutine iter_bulk_2006
