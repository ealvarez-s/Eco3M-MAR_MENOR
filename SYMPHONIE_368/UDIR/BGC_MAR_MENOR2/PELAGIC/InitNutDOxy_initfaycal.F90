
      SUBROUTINE InitNutDOxy

!_____________________________________________________________________*
! 3d ecosystem model                                                  *
!                                                                     *
! LAST REVISION: 16 MAY 2007                                           *
!                                                                     *
! Implementation: Caroline Ulses                                      *
!                 NIOO-CEME                                           *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! 1. Initialises the general parameters                               *
! 2. Reads the light and pelagic parameters                           *
! 3. Initialises the variables                                        *
! 4. Writes the initial profiles for the Femme Environment            *
!_____________________________________________________________________*

!_____________________________________________________________________*
!                                                                     *
! Modifications:                                                      *
! 21/02/07: include symphonie2007.h dans la version 2004_v7           *
! 16/05/07: version a 3 phytos, 3 zoo et bacteries                    *
! 03/12/07: modifs Claude pour frontieres Est                         *
! 18/07/08: ajout d'une option nudging                                *
!  21/03/2013: modification dans les initialisations bio              *
!_____________________________________________________________________*



!---------------------------------------------------------------------*
! Declarations:


! Global variables

      use ModuleComBlock
      use UserDeclaration
      use ModuleDeclaration
      use module_principal
      use module_parallele !#MPI
      use module_global
      use module_offline
      implicit none
!      include 'netcdf.inc'

! Local variables
      integer ChoiceInitPelagic,NumVar(20),LDir,K_pr  &
             ,InitFemme,ifemme,JFemme,nmaxlpc,nmax,nmax1,nmax2,nmax3    &
             ,isc2160,jsc2160,isw2060,jsw2060,nmax4,nmax5,nmax6         &
             ,isc2240,jsc2240,ij,InitOC

      character DATE*6,DirInitFile*90,Fic(20)*80      &
               ,NameVar(vbmax)*9,FicInitVar*80        &
               ,FicLPC_Pho*80,FicLPC_Nit*80


      double precision InitValue(20,3000),Prof(20,3000),Prof_oxyg(20,3000),Prof_oxyg_WEST(20,3000)       &
                      ,Prof_oxyg_EST(20,3000),InitValue_oxyg_EST(20,3000) &
                      ,InitValue_GLT_DIC(20,3000)    &
                      ,InitValue_EST_WEST_DIC(20,3000)    &
                      ,InitValue_EST_EST_DIC(20,3000)    &
                      ,InitValue_EST_DIC(20,3000)    &
                      ,InitValue_EST_WEST_ALK(20,3000)    &
                      ,InitValue_EST_EST_ALK(20,3000)    &
                      ,InitValue_EST_ALK(20,3000)    &
                      ,InitValue_WEST_DIC(20,3000)    &
                      ,InitValue_EST(20,3000)   &
                      ,InitValue_WEST(20,3000)   &
                      ,InitValue_cyp(20,3000)   &
                      ,InitValue_ion(20,3000)   &
                      ,InitValue_lib(20,3000)   &
                      ,InitValue_tun(20,3000)   &
                      ,InitValue_sic(20,3000)   &
                      ,InitValue_nad(20,3000)   &
                      ,InitValue_sad(20,3000)   &
                      ,InitValue_mlt(20,3000)   &
                      ,InitValue_ege(20,3000)   &
                      ,InitValue_oxyg_cyp(20,3000)   &
                      ,InitValue_oxyg_ion(20,3000)   &
                      ,InitValue_oxyg_lib(20,3000)   &
                      ,InitValue_oxyg_tun(20,3000)   &
                      ,InitValue_oxyg_sic(20,3000)   &
                      ,InitValue_oxyg_nad(20,3000)   &
                      ,InitValue_oxyg_sad(20,3000)   &
                      ,InitValue_oxyg_mlt(20,3000)   &
                      ,InitValue_oxyg_ege(20,3000)   &
                      ,CoefDiaC,CoefZooc,CoefSyneC           &
                      ,MasseMolC,MasseMolN                   &
                      ,RedfieldNC,RedfieldPC                 &
                      ,PropMOP                               &
                      ,InitLPOSi,InitAmmonium                &
                      ,InitValueN(99),ProfN(99),CoefNanoC    &
                      ,TChl(0:imax+1,0:jmax+1,kmax)          &
                      ,CNSMOP,ChlNDia                        &
                      ,XP1,XP2,YP1,YP2                       &
                      ,InitValueLPC(20,3000),ProfLPC(20,3000),lon,lat &
                      ,salz,salmu,salsigma,p1,p2,p3,p4,p5&
                      ,InitPH(1,22),ProfPH(1,22)            &
                      ,InitTAlk(5,22),ProfTALk(5,22)            &
                      ,InitDIC(5,22),ProfDIC(5,22)

      real InitPH3D(iglb+2,jglb+2,kmax),InitPH2D((iglb+2)*(jglb+2),kmax), &
           InitPH3D2(iglb+2,jglb+2,kmax),InitMODC(iglb+2,jglb+2,kmax)            
 
       integer l,m,nmax10y(14),nmax_oxy_lev,nmax_oxy_ion,nmax_oxy_tyr, &
               nmax_oxy_alb,nmax_oxy_alg
       real Prof10y(14,2511),var10y(14,2511,10)
       character*58 filename


       InitOC=1



      nmax1=68
      nmax2=68
      nmax3=68
      nmax4=21
      nmax5=22
!     nmax6=183
!     nmax6=2311
      nmax6=2351

       

! to get initial temperature and salinity fields  
           if(ioffline==2) call offline_inout(2)   
!$ water density:                                                      !14/11/04
!      if(eqs_state1.eq.0)call density(0) !lineaire
!      if(eqs_state1.eq.1)call density(3) !Non line sans pr
       call equation_of_state('potential density',1)

! nutrients


!!!!!! !!!!!!!!!!!!!!!!!! LOADING PROFILES FROM OBSERVATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!load oxygene EST et OUEST

!!      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_oxyg_est.txt')
!!      read(3,*)nmax3
!     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_oxyg_EST2.txt')
!      do k=1,nmax4
!      read(3,*)InitValue_EST(1,k),Prof_oxyg(1,k)
!      enddo
!      close(3)


      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitRegCTDOxyMeteor_lev.txt')
      read(3,*)nmax_oxy_lev
      do k=1,nmax_oxy_lev
      read(3,*)InitValue_EST(6,k),Prof_oxyg_EST(6,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitRegCTDOxyMeteor_ion.txt')
      read(3,*)nmax_oxy_ion
      do k=1,nmax_oxy_ion
      read(3,*)InitValue_EST(5,k),Prof_oxyg_EST(5,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitRegCTDOxyMeteor_tyr.txt')
      read(3,*)nmax_oxy_tyr
      do k=1,nmax_oxy_tyr
      read(3,*)InitValue_EST(4,k),Prof_oxyg_EST(4,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitRegCTDOxyMeteor_alg.txt')
      read(3,*)nmax_oxy_alg
      do k=1,nmax_oxy_alg
      read(3,*)InitValue_EST(3,k),Prof_oxyg_EST(3,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitRegCTDOxyMeteor_alb.txt')
      read(3,*)nmax_oxy_alb
      do k=1,nmax_oxy_alb
      read(3,*)InitValue_EST(2,k),Prof_oxyg_EST(2,k)
      enddo
      close(3)

!     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/OXYGEN_WEST.txt')
      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitOxyDyf201108.txt')
      do k=1,nmax6
      read(3,*)InitValue_WEST(1,k),Prof_oxyg_WEST(1,k)
      enddo
      close(3)



! load nitrate

! cyp

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_cyp.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_cyp(1,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_cyp.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_cyp(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_cyp.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_cyp(3,k)
      enddo
      close(3)


! ion

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_ion.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_ion(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_ion.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_ion(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_ion.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_ion(3,k)
      enddo
      close(3)


! lib

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_lib.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_lib(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_lib.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_lib(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_lib.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_lib(3,k)
      enddo
      close(3)

! tun

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_tun.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_tun(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_tun.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_tun(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_tun.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_tun(3,k)
      enddo
      close(3)


! sic

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_sic.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_sic(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_sic.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_sic(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_sic.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_sic(3,k)
      enddo
      close(3)

! nad

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_nad.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_nad(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_nad.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_nad(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_nad.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_nad(3,k)
      enddo
      close(3)


! sad

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_sad.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_sad(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_sad.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_sad(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_sad.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_sad(3,k)
      enddo
      close(3)

! mlt

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_mlt.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_mlt(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_mlt.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_mlt(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_mlt.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_mlt(3,k)
      enddo
      close(3)


! ege

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_nit_ege.txt')
!      read(3,*)nmax1
      do k=1,nmax1
      read(3,*)Prof(1,k),InitValue_ege(1,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_phos_ege.txt')
!      read(3,*)nmax2
      do k=1,nmax2
      read(3,*)Prof(2,k),InitValue_ege(2,k)
      enddo
      close(3)


      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_sil_ege.txt')
!      read(3,*)nmax3
      do k=1,nmax3
      read(3,*)Prof(3,k),InitValue_ege(3,k)
      enddo
      close(3)



!!!!!!!!!!!!!!!!!!!! END PROFILE LOADING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! bioregion 11 tyrrhenien

! coeff pour la relation oxygène/salinité polynôme 4
      salmu = 38.49 ! moyenne des salinité des mesures ctd
      salsigma = 0.085254 ! ecart type des salinité
      p1 = -8.9604d-5
      p2 = -0.016762
      p3 = -0.9261
      p4 = -15.234
      p5 = 203.9


      do i=mbio1,mbio2
      do j=nbio1,nbio2
      do k=kmin_w(i,j),kmax


       if (mask_t(i,j,k)==1) then

!       print*,'Passe ici dans InitNut',lon_t(i,j)/deg2rad,lat_t(i,j)/deg2rad
!       stop'dans InitNut'

         if (                              &
         lon_t(i,j)>=9.5*deg2rad.and.                            & 
         lon_t(i,j)<12.5*deg2rad.and.                            &
         lat_t(i,j)>41.5*deg2rad.and.               &
         lat_t(i,j)<44.3*deg2rad)then  !.and.lat_t(i,j)<44.*deg2rad)then




         bio_t(i,j,k,iNITRATE)= (0.0400278*(sal_t(i,j,k,1)**2)-57.0955)**2
         if(bio_t(i,j,k,iNITRATE)>=7.and.depth_t(i,j,k)<-500)           &
         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513
         if(depth_t(i,j,k)<-800)                                        &
         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513

         if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
         if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08097*(rhp_t(i,j,k)+rho-1000.)-30.8847)

         bio_t(i,j,k,iPHOSPHATE)= (-9.93055 + (0.00700501*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iPHOSPHATE)>=7.and.depth_t(i,j,k)<-500)           &
         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
         if(depth_t(i,j,k)<-800)                                        &
         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974

         if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
         if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.937182*(rhp_t(i,j,k)+rho-1000.)-29.4776)

         bio_t(i,j,k,iSILICE)= (-62.7026 + (0.0437829*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iSILICE)>=7.and.depth_t(i,j,k)<-500)           &
         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
         if(depth_t(i,j,k)<-800)                                        &
         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212

         if(bio_t(i,j,k,iSILICE)>=7.5) bio_t(i,j,k,iSILICE)=7.5
         if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.08691*(rhp_t(i,j,k)+rho-1000.)-31.0381)


! oxygene
!        bio_t(i,j,k,34)=200.
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5


          if (bio_t(i,j,k,iOXYGEN)>165.and.depth_t(i,j,k)<-200) bio_t(i,j,k,iOXYGEN)=165 

          if (depth_t(i,j,k)<-800.and.depth_t(i,j,k)>-1000) bio_t(i,j,k,iOXYGEN)=185

         if (depth_t(i,j,k)<=-1000.and.depth_t(i,j,k)>=-1500) bio_t(i,j,k,iOXYGEN)=190

!        if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=195
         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond


         if(bio_t(i,j,k,iOXYGEN)>=250) bio_t(i,j,k,iOXYGEN)=250

         if(bio_t(i,j,k,iOXYGEN)<=165) bio_t(i,j,k,iOXYGEN)=165


!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

!         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond 
!end oxygen




        elseif (                               &
        lon_t(i,j)>=9.5*deg2rad.and.                                 &
        lon_t(i,j)<16.*deg2rad.and.                                 &
        lat_t(i,j)<=41.5*deg2rad.and.lat_t(i,j)>38.*deg2rad)then



         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513

         if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
         if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08097*(rhp_t(i,j,k)+rho-1000.)-30.8847)

         bio_t(i,j,k,iPHOSPHATE)= (-9.93055 + (0.00700501*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iPHOSPHATE)>=7.and.depth_t(i,j,k)<-500)           &
         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
         if(depth_t(i,j,k)<-800)                                        &
         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974

         if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
         if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.937182*(rhp_t(i,j,k)+rho-1000.)-29.4776)

         bio_t(i,j,k,iSILICE)= (-62.7026 + (0.0437829*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iSILICE)>=7.and.depth_t(i,j,k)<-500)           &
         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
         if(depth_t(i,j,k)<-800)                                        &
         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212

         if(bio_t(i,j,k,iSILICE)>=7.5) bio_t(i,j,k,iSILICE)=7.5
         if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.08691*(rhp_t(i,j,k)+rho-1000.)-31.0381)


! oxygene
!        bio_t(i,j,k,34)=200.
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5


         if (bio_t(i,j,k,iOXYGEN)>165.and.depth_t(i,j,k)<-200) bio_t(i,j,k,iOXYGEN)=165 
         if (depth_t(i,j,k)<-800.and.depth_t(i,j,k)>-1000) bio_t(i,j,k,iOXYGEN)=185
         if (depth_t(i,j,k)<=-1000.and.depth_t(i,j,k)>=-1500) bio_t(i,j,k,iOXYGEN)=190
!        if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=195
         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond


         if(bio_t(i,j,k,iOXYGEN)>=250) bio_t(i,j,k,iOXYGEN)=250

         if(bio_t(i,j,k,iOXYGEN)<=165) bio_t(i,j,k,iOXYGEN)=165


!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

!         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond 
!end oxygen




! bioregion 12  subbassin algero-provençal

         elseif (                               &
         lon_t(i,j)>-2.*deg2rad.and.                                 &
         lon_t(i,j)<9.5*deg2rad.and.                                 &
         lat_t(i,j)<39.5*deg2rad)then    ! 38.


         if(depth_t(i,j,k)<-400) bio_t(i,j,k,iNITRATE)=9.
!         if(depth_t(i,j,k)<-1800)                                        &
!         bio_t(i,j,k,iNITRATE)=1./(36.0319 - (1.23423*(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iNITRATE)>=10.) bio_t(i,j,k,iNITRATE)=10.
         if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(0.690112*(rhp_t(i,j,k)+rho-1000.)-19.7654)

         bio_t(i,j,k,iPHOSPHATE)= (-7.0582 + (0.00517724*sal_t(i,j,k,1)**2))**2
!         if(bio_t(i,j,k,iPHOSPHATE)>=9.and.depth_t(i,j,k)<-1000)           &
!         bio_t(i,j,k,iPHOSPHATE)=exp(164.933 - (4828.84/(rhp_t(i,j,k)+rho-1000.)))
         if(depth_t(i,j,k)<-1800)                                        &
         bio_t(i,j,k,iPHOSPHATE)=exp(164.933 - (4828.84/(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
         if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.560291*(rhp_t(i,j,k)+rho-1000.)-18.9134)
        
        bio_t(i,j,k,iSILICE)= (-33.0298 + (0.0242586*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iSILICE)<=8.5.and.depth_t(i,j,k)<-100)           &
         bio_t(i,j,k,iSILICE)=8.5
!         bio_t(i,j,k,iSILICE)= 1./(-37.9453 + (1108.41/(rhp_t(i,j,k)+rho-1000.)))
         if(depth_t(i,j,k)<-300)                                        &
         bio_t(i,j,k,iNITRATE)=9.
!         bio_t(i,j,k,iSILICE)= 1./(-37.9453 + (1108.41/(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iSILICE)>=8.5) bio_t(i,j,k,iSILICE)=8.5
         if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(0.769009*(rhp_t(i,j,k)+rho-1000.)-21.6826)


! oxygene
!        bio_t(i,j,k,34)=200.
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5

 
         if (bio_t(i,j,k,iOXYGEN)>165.and.depth_t(i,j,k)<-200) bio_t(i,j,k,iOXYGEN)=165 
         if (depth_t(i,j,k)<-800.and.depth_t(i,j,k)>-1000) bio_t(i,j,k,iOXYGEN)=185
         if (depth_t(i,j,k)<=-1000.and.depth_t(i,j,k)>=-1500) bio_t(i,j,k,iOXYGEN)=190
!        if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=195
         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond


         if(bio_t(i,j,k,iOXYGEN)>=250) bio_t(i,j,k,iOXYGEN)=250

         if(bio_t(i,j,k,iOXYGEN)<=165) bio_t(i,j,k,iOXYGEN)=165


!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

!         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond 
!end oxygen





!         if(depth_t(i,j,k)>-60.)then
!         bio_t(i,j,k,iNITRATE)=0.05
!         bio_t(i,j,k,iPHOSPHATE)=0.02
!         bio_t(i,j,k,iSILICE)=0.05
!         endif


! sicile
        elseif(                                  &
        lon_t(i,j)>=9.5*deg2rad.and.lon_t(i,j)<11.*deg2rad.and.         &
        lat_t(i,j)>36.5*deg2rad.and.lat_t(i,j)<=38.*deg2rad)then

         
!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_sic(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_sic(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_sic(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_sic(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif


!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_sic(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_sic(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_sic(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_sic(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_sic(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_sic(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_sic(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_sic(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice


!         if(depth_t(i,j,k)>-60.)then
!         bio_t(i,j,k,iNITRATE)=0.05
!         bio_t(i,j,k,iPHOSPHATE)=0.02
!         bio_t(i,j,k,iSILICE)=0.05
!         endif

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen




! malte
        elseif (                                 &
        lon_t(i,j)>10.*deg2rad.and.lon_t(i,j)<=15.*deg2rad.and.        &
        lat_t(i,j)>34.*deg2rad.and.lat_t(i,j)<=38.*deg2rad)then
 


      
!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_mlt(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_mlt(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_mlt(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_mlt(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_mlt(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_mlt(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_mlt(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_mlt(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_mlt(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_mlt(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_mlt(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_mlt(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen




! libye ouest
! libye est
        elseif(                                  &
        lon_t(i,j)>15.*deg2rad.and.lon_t(i,j)<=30.*deg2rad.and.        &
        lat_t(i,j)<=37.*deg2rad)then



!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_lib(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_lib(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_lib(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_lib(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_lib(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_lib(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_lib(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_lib(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_lib(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_lib(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_lib(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_lib(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



! tun
        elseif(                                  &
        lon_t(i,j)>9.5*deg2rad.and.lon_t(i,j)<=15.*deg2rad.and.        &
!        lat_t(i,j)>36.*deg2rad.and.lat_t(i,j)<38.*deg2rad)then
        lat_t(i,j)<=34.*deg2rad)then


       
!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_tun(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_tun(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_tun(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_tun(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_tun(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_tun(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_tun(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_tun(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_tun(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_tun(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_tun(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_tun(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



! ionnienne
        elseif(                                  &
        lon_t(i,j)>15.*deg2rad.and.lon_t(i,j)<22.*deg2rad.and.        &
        lat_t(i,j)>37.*deg2rad.and.lat_t(i,j)<39.*deg2rad)then
!        lat_t(i,j)<39.*deg2rad)then



!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_ion(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_ion(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_ion(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_ion(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_ion(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_ion(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_ion(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_ion(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_ion(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_ion(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_ion(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_ion(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



! cypre
        elseif(lon_t(i,j)>=30.*deg2rad)then!.and.   &
!        lat_t(i,j)<36.*deg2rad)then
!         lat_t(i,j)<34.*deg2rad)then



!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_cyp(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_cyp(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_cyp(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_cyp(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_cyp(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_cyp(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_cyp(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_cyp(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_cyp(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_cyp(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_cyp(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_cyp(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



! NAD
        elseif(                                  &
        lon_t(i,j)>11.*deg2rad.and.lon_t(i,j)<20.*deg2rad.and.        &
        lat_t(i,j)>41.*deg2rad )then


!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_nad(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_nad(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_nad(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_nad(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_nad(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_nad(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_nad(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_nad(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_nad(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_nad(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_nad(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_nad(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



! SAD
        elseif(                                  &
        lon_t(i,j)>16.*deg2rad.and.lon_t(i,j)<20.5*deg2rad.and.        &
        lat_t(i,j)>39.*deg2rad.and.lat_t(i,j)<41.*deg2rad)then



!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_sad(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_sad(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_sad(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_sad(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_sad(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_sad(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_sad(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_sad(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_sad(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_sad(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_sad(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_sad(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



! ege
        elseif(                                   &
        lon_t(i,j)>=21.*deg2rad.and.lon_t(i,j)<30.*deg2rad.and.      &
        lat_t(i,j)>=37.*deg2rad)then



!nitrate
             if    (-depth_t(i,j,k)<=Prof(1,1   ))then
               bio_t(i,j,k,iNitrate) = InitValue_ege(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof(1,nmax1))then
               bio_t(i,j,k,iNitrate) = InitValue_ege(1,nmax1)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iNitrate)=(                              &
              ( Prof(1,K_pr  )+depth_t(i,j,k))*InitValue_ege(1,K_pr-1) &
             +(-Prof(1,K_pr-1)-depth_t(i,j,k))*InitValue_ege(1,K_pr  ) &
             )/(Prof(1,K_pr  )-Prof(1,K_pr-1))
           endif
!end nitrate
! phosphate
             if    (-depth_t(i,j,k)<=Prof(2,1   ))then
               bio_t(i,j,k,iPhosphate) = InitValue_ege(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(2,nmax2))then
               bio_t(i,j,k,iPhosphate) = InitValue_ege(2,nmax2)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iPhosphate)=(                              &
              ( Prof(2,K_pr  )+depth_t(i,j,k))*InitValue_ege(2,K_pr-1) &
             +(-Prof(2,K_pr-1)-depth_t(i,j,k))*InitValue_ege(2,K_pr  ) &
             )/(Prof(2,K_pr  )-Prof(2,K_pr-1))
           endif
!end phosphate

!silice
             if    (-depth_t(i,j,k)<=Prof(3,1   ))then
               bio_t(i,j,k,iSilice) = InitValue_ege(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof(3,nmax3))then
               bio_t(i,j,k,iSilice) = InitValue_ege(3,nmax3)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iSilice)=(                              &
              ( Prof(3,K_pr  )+depth_t(i,j,k))*InitValue_ege(3,K_pr-1) &
             +(-Prof(3,K_pr-1)-depth_t(i,j,k))*InitValue_ege(3,K_pr  ) &
             )/(Prof(3,K_pr  )-Prof(3,K_pr-1))
           endif
!end silice

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg(1,nmax4))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(1,nmax4)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              ( Prof_oxyg(1,K_pr  )+depth_t(i,j,k))*InitValue_EST(1,K_pr-1) &
             +(-Prof_oxyg(1,K_pr-1)-depth_t(i,j,k))*InitValue_EST(1,K_pr  ) &
             )/(Prof_oxyg(1,K_pr  )-Prof_oxyg(1,K_pr-1))
           endif
!end oxygen



!!  changement dans laboran et gibraltar: tester a gibralat la relation dans alboran

! alboran
        elseif(                                  &
        lon_t(i,j)<=-2.*deg2rad.and.lon_t(i,j)>-5.*deg2rad)then



         bio_t(i,j,k,iNITRATE)= (0.0389*(sal_t(i,j,k,1)**2)-49.4811)**2
         if(bio_t(i,j,k,iNITRATE)>=7.and.depth_t(i,j,k)<-250)           &
        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924 
         if(depth_t(i,j,k)<-300)                                        &
        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924 
         if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(2.00125*(rhp_t(i,j,k)+rho-1000.)-53.8274) 

         bio_t(i,j,k,iPHOSPHATE)= (-9.09731 + (0.00656935*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iPHOSPHATE)>=7.and.depth_t(i,j,k)<-250)           &
        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102   
         if(depth_t(i,j,k)<-300)                                        &
        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102   
        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(1.86852*(rhp_t(i,j,k)+rho-1000.)-52.7586) 

         bio_t(i,j,k,iSILICE)= (-21.0933 + (0.0162298*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iSILICE)>=7.and.depth_t(i,j,k)<-250)           &
        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489  
         if(depth_t(i,j,k)<-300)                                        &
        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489  
        if(bio_t(i,j,k,iSILICE)>=7.5) bio_t(i,j,k,iSILICE)=7.5
        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.78793*(rhp_t(i,j,k)+rho-1000.)-48.4092)  


! oxygene
!        bio_t(i,j,k,34)=200.
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5



         if (bio_t(i,j,k,iOXYGEN)>165.and.depth_t(i,j,k)<-200) bio_t(i,j,k,iOXYGEN)=165 
         if (depth_t(i,j,k)<-800.and.depth_t(i,j,k)>-1000) bio_t(i,j,k,iOXYGEN)=185
         if (depth_t(i,j,k)<=-1000.and.depth_t(i,j,k)>=-1500) bio_t(i,j,k,iOXYGEN)=190
!        if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=195
         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond


         if(bio_t(i,j,k,iOXYGEN)>=250) bio_t(i,j,k,iOXYGEN)=250

         if(bio_t(i,j,k,iOXYGEN)<=165) bio_t(i,j,k,iOXYGEN)=165


!         if(depth_t(i,j,k)>-60.)then
!         bio_t(i,j,k,iNITRATE)=0.05
!         bio_t(i,j,k,iPHOSPHATE)=0.02
!         bio_t(i,j,k,iSILICE)=0.05
!         endif


!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

!         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond 
!end oxygen




!! gibraltar
        elseif(                                  &
        lon_t(i,j)<=-5.*deg2rad)then
!


        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489


!
!!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!!inflow
        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iNITRATE)=0.4
        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iPHOSPHATE)=0.12
        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iSILICE)=0.79
!!outflow
        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iNITRATE)=8.
        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iPHOSPHATE)=0.4
        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iSILICE)=6.47
!!!!!!!!!!!!!!!!!!fin gibralar

! oxygene
!        bio_t(i,j,k,34)=200.
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5

         if (bio_t(i,j,k,iOXYGEN)>165.and.depth_t(i,j,k)<-200) bio_t(i,j,k,iOXYGEN)=165 
    
         if (depth_t(i,j,k)<-800.and.depth_t(i,j,k)>-1000) bio_t(i,j,k,iOXYGEN)=185
        
         if (depth_t(i,j,k)<=-1000.and.depth_t(i,j,k)>=-1500) bio_t(i,j,k,iOXYGEN)=190
        
!        if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=195
         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond         


         if(bio_t(i,j,k,iOXYGEN)>=250) bio_t(i,j,k,iOXYGEN)=250
         
         if(bio_t(i,j,k,iOXYGEN)<=165) bio_t(i,j,k,iOXYGEN)=165

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

!         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond 
!end oxygen




!region nord occidental (provençale + ligure)

       else
         


         bio_t(i,j,k,iNITRATE)= exp(0.0508111*(sal_t(i,j,k,1))**2 - 73.2698)
!         if(bio_t(i,j,k,iNITRATE)>=9.and.depth_t(i,j,k)<-250)           &
!         bio_t(i,j,k,iNITRATE)=1./(36.0319 - (1.23423*(rhp_t(i,j,k)+rho-1000.))) 
         if(bio_t(i,j,k,iNITRATE)<=8.5.and.depth_t(i,j,k)<-200)           &
          bio_t(i,j,k,iNITRATE)=8.5
         if(depth_t(i,j,k)<-400)                                        &
          bio_t(i,j,k,iNITRATE)=9.25
!         bio_t(i,j,k,iNITRATE)=1./(36.0319 - (1.23423*(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iNITRATE)>=9.5) bio_t(i,j,k,iNITRATE)=9.5
      if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(0.690112*(rhp_t(i,j,k)+rho-1000.)-19.7654)

         bio_t(i,j,k,iPHOSPHATE)= (-16.5165 + (0.0115393*sal_t(i,j,k,1)**2))**2
         if(bio_t(i,j,k,iPHOSPHATE)<=0.42.and.depth_t(i,j,k)<-200)           &
         bio_t(i,j,k,iPHOSPHATE)=0.42
!         bio_t(i,j,k,iPHOSPHATE)=exp(164.933 - (4828.84/(rhp_t(i,j,k)+rho-1000.)))  
!         if(depth_t(i,j,k)<-300)                                        &
!         bio_t(i,j,k,iPHOSPHATE)=exp(164.933 - (4828.84/(rhp_t(i,j,k)+rho-1000.)))  
          if(bio_t(i,j,k,iPHOSPHATE)>=0.42) bio_t(i,j,k,iPHOSPHATE)=0.42
       if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.560291*(rhp_t(i,j,k)+rho-1000.)-18.9134)

         bio_t(i,j,k,iSILICE)= 1./(22.9274 - (0.0153768*(sal_t(i,j,k,1))**2))
         if(bio_t(i,j,k,iSILICE)<=8.and.depth_t(i,j,k)<-100)           &
         bio_t(i,j,k,iSILICE)=8.
!         if(bio_t(i,j,k,iSILICE)>=9.and.depth_t(i,j,k)<-250)           &
!         bio_t(i,j,k,iSILICE)= 1./(-37.9453 + (1108.41/(rhp_t(i,j,k)+rho-1000.)))
         if(depth_t(i,j,k)<-200)                                        &
          bio_t(i,j,k,iSILICE)= 8.
!         bio_t(i,j,k,iSILICE)= 1./(-37.9453 + (1108.41/(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
       if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(0.769009*(rhp_t(i,j,k)+rho-1000.)-21.6826)



!      endif

! oxygene
!        bio_t(i,j,k,34)=200.
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5

         if (bio_t(i,j,k,iOXYGEN)>165.and.depth_t(i,j,k)<-200) bio_t(i,j,k,iOXYGEN)=165 
       
         if (depth_t(i,j,k)<-800.and.depth_t(i,j,k)>-1000) bio_t(i,j,k,iOXYGEN)=185
        
         if (depth_t(i,j,k)<=-1000.and.depth_t(i,j,k)>=-1500) bio_t(i,j,k,iOXYGEN)=190
        
!        if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=195
         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond         


         if(bio_t(i,j,k,iOXYGEN)>=250) bio_t(i,j,k,iOXYGEN)=250
         
         if(bio_t(i,j,k,iOXYGEN)<=165) bio_t(i,j,k,iOXYGEN)=165

!oxygen
             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

!         if (depth_t(i,j,k)<=-1500) bio_t(i,j,k,iOXYGEN)=205 ! pour une initialisation en 2010 apres les evts de formation d'eau dense 2005 2006 2009 2010 on prefere partir avec une + forte valeur au fond 
!end oxygen

         


           endif


! Other regions for oxygen

! Alboran                      
        if(lon_t(i,j)<=-1.*deg2rad.and.lon_t(i,j)>-6.*deg2rad)then

       print*,'Passe ici dans InitNut',lon_t(i,j)/deg2rad,lat_t(i,j)/deg2rad
       print*,'-depth_t(i,j,k)',k,-depth_t(i,j,k),Prof_oxyg_EST(2,1   )
       print*,'2',InitValue_EST(2,1),InitValue_EST(2,nmax_oxy_alb)
!      stop'dans InitNut'

             if    (-depth_t(i,j,k)<=Prof_oxyg_EST(2,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EST(2,nmax_oxy_alb))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(2,nmax_oxy_alb)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EST(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EST(2,K_pr)+depth_t(i,j,k))*InitValue_EST(2,K_pr-1) &
             +(-Prof_oxyg_EST(2,K_pr-1)-depth_t(i,j,k))*InitValue_EST(2,K_pr) &
             )/(Prof_oxyg_EST(2,K_pr  )-Prof_oxyg_EST(2,K_pr-1))

       print*,'3',InitValue_EST(2,K_pr),InitValue_EST(2,K_pr-1),InitValue_EST(2,K_pr)

           endif 

        endif

! Algero-Provençal                      
        if(lon_t(i,j)<=10.*deg2rad.and.lon_t(i,j)>-1.*deg2rad)then

             if    (-depth_t(i,j,k)<=Prof_oxyg_EST(3,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EST(3,nmax_oxy_alg))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(3,nmax_oxy_alg)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EST(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EST(3,K_pr)+depth_t(i,j,k))*InitValue_EST(3,K_pr-1) &
             +(-Prof_oxyg_EST(3,K_pr-1)-depth_t(i,j,k))*InitValue_EST(3,K_pr) &
             )/(Prof_oxyg_EST(3,K_pr  )-Prof_oxyg_EST(3,K_pr-1))
           endif

        endif

! Tyrhennienne + sicile                      
        if(lon_t(i,j)<=17.*deg2rad.and.lon_t(i,j)>10.*deg2rad)then

             if    (-depth_t(i,j,k)<=Prof_oxyg_EST(4,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(4,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EST(4,nmax_oxy_tyr))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(4,nmax_oxy_tyr)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EST(4,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EST(4,K_pr)+depth_t(i,j,k))*InitValue_EST(4,K_pr-1) &
             +(-Prof_oxyg_EST(4,K_pr-1)-depth_t(i,j,k))*InitValue_EST(4,K_pr) &
             )/(Prof_oxyg_EST(4,K_pr  )-Prof_oxyg_EST(4,K_pr-1))
           endif

        endif

! Ionienne                       
        if(lon_t(i,j)<=27.*deg2rad.and.lon_t(i,j)>17.*deg2rad)then
               
             if    (-depth_t(i,j,k)<=Prof_oxyg_EST(5,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(5,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EST(5,nmax_oxy_ion))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(5,nmax_oxy_ion)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EST(5,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EST(5,K_pr)+depth_t(i,j,k))*InitValue_EST(5,K_pr-1) &
             +(-Prof_oxyg_EST(5,K_pr-1)-depth_t(i,j,k))*InitValue_EST(5,K_pr) &
             )/(Prof_oxyg_EST(5,K_pr  )-Prof_oxyg_EST(5,K_pr-1))
           endif

        endif

! Levantine                      
        if(lon_t(i,j)>27.*deg2rad)then
               
             if    (-depth_t(i,j,k)<=Prof_oxyg_EST(6,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(6,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EST(6,nmax_oxy_lev))then
               bio_t(i,j,k,iOxygen) = InitValue_EST(6,nmax_oxy_lev)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EST(6,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EST(6,K_pr)+depth_t(i,j,k))*InitValue_EST(6,K_pr-1) &
             +(-Prof_oxyg_EST(6,K_pr-1)-depth_t(i,j,k))*InitValue_EST(6,K_pr) &
             )/(Prof_oxyg_EST(6,K_pr  )-Prof_oxyg_EST(6,K_pr-1))
           endif

        endif

! deep NV                       
        if(lon_t(i,j)>2.*deg2rad.and.lon_t(i,j)<9.*deg2rad.and. &
           lat_t(i,j)>41.*deg2rad.and.h_w(i,j).gt.2000)then

             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax6))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax6)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))
           endif

        endif



          bio_t(i,j,k,iOxygen)       =bio_t(i,j,k,iOxygen)  &
                            *(rhp_t(i,j,k)+rho) /1000.



      endif ! mask 
      enddo
      enddo
      enddo


       if(par%rank==0) write(6,*)'Passe par InitNutDoxy'
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!09-05-10
#endif

      END SUBROUTINE InitNutDOxy
