
      SUBROUTINE InitDOxy

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
                      ,Prof_oxyg_EAST(20,3000),InitValue_oxyg_EAST(20,3000) &
                      ,InitValue_GLT_DIC(20,3000)    &
                      ,InitValue_EST_WEST_DIC(20,3000)    &
                      ,InitValue_EST_EST_DIC(20,3000)    &
                      ,InitValue_EST_DIC(20,3000)    &
                      ,InitValue_EST_WEST_ALK(20,3000)    &
                      ,InitValue_EST_EST_ALK(20,3000)    &
                      ,InitValue_EST_ALK(20,3000)    &
                      ,InitValue_WEST_DIC(20,3000)    &
                      ,InitValue_EAST(20,3000)   &
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
               nmax_oxy_alb,nmax_oxy_alg,nmax_oxy_cen,nmax_oxy_gol
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

!     open(unit=2,file='/tmpdir/culses/bassin2/MERCATOFF/DyfamedmeanpH.txt'))
      open(unit=2,file='/tmpdir/jhabib/simulations/initdyfamed/DyfamedmeanpH.txt')
      do k=1,nmax5
      read(2,*)ProfPH(1,k),InitPH(1,k)
      enddo
      close(2)


!      open(unit=2,file='../../../MERCATOFF/BassinInitco2syspH2.txt')
!      do k=1,kmax     
!      read(2,*)(InitPH2D(ij,k),ij=1,(iglb+2)*(jglb+2))
!      enddo
!      close(2)

!      do i=1,iglb+2
!      do j=1,jglb+2
!      do k=1,kmax 
!      ij=j+(jglb+2)*(i-1)
!      InitPH3D(i,j,k)=InitPH2D(ij,k)
!      enddo
!      enddo
!      enddo

!       do k=1,kmax
!       print*,'phinit',k,iglb+2,jglb+2,(iglb+2)*(jglb+2),InitPH3D(52,82,k) 
!       enddo
       

! to get initial temperature and salinity fields  
           if(ioffline==2) call offline_inout(2)   
!$ water density:                                                      !14/11/04
!      if(eqs_state1.eq.0)call density(0) !lineaire
!      if(eqs_state1.eq.1)call density(3) !Non line sans pr
       call equation_of_state('potential density',1)

! nutrients


!!!!!! !!!!!!!!!!!!!!!!!! LOADING PROFILES FROM OBSERVATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!load oxygene EST et OUEST

!!      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/interp_oxyg_est.txt')
!!      read(3,*)nmax3
!     open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/interp_oxyg_EST2.txt')
!      do k=1,nmax4
!      read(3,*)InitValue_EST(1,k),Prof_oxyg(1,k)
!      enddo
!      close(3)


      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyMeteor_lev2.txt')
      read(3,*)nmax_oxy_lev
      do k=1,nmax_oxy_lev
      read(3,*)InitValue_EAST(1,k),Prof_oxyg_EAST(1,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyBoum_ion.txt')
      read(3,*)nmax_oxy_ion
      do k=1,nmax_oxy_ion
      read(3,*)InitValue_EAST(2,k),Prof_oxyg_EAST(2,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyMeteor_tyr2.txt')
      read(3,*)nmax_oxy_tyr
      do k=1,nmax_oxy_tyr
      read(3,*)InitValue_WEST(4,k),Prof_oxyg_WEST(4,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyMeteor_alg.txt')
      read(3,*)nmax_oxy_alg
      do k=1,nmax_oxy_alg
      read(3,*)InitValue_WEST(2,k),Prof_oxyg_WEST(2,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyBoum_cen.txt')
      read(3,*)nmax_oxy_cen
      do k=1,nmax_oxy_cen
      read(3,*)InitValue_EAST(3,k),Prof_oxyg_EAST(3,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyMeteor_alb.txt')
      read(3,*)nmax_oxy_alb
      do k=1,nmax_oxy_alb
      read(3,*)InitValue_WEST(1,k),Prof_oxyg_WEST(1,k)
      enddo
      close(3)

!     open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitRegCTDOxyBoum_gol.txt')
      open(unit=3,file='/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitOxyDyf201108.txt') !11/12/17
      read(3,*)nmax_oxy_gol
      do k=1,nmax_oxy_gol
      read(3,*)InitValue_WEST(3,k),Prof_oxyg_WEST(3,k)
      enddo
      close(3)



!!!!!!!!!!!!!!!!!!!! END PROFILE LOADING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do i=mbio1,mbio2
      do j=nbio1,nbio2
      do k=kmin_w(i,j),kmax


       if (mask_t(i,j,k)==1) then

!       print*,'Passe ici dans
!       InitNut',lon_t(i,j)/deg2rad,lat_t(i,j)/deg2rad
!       stop'dans InitNut'


! bioregion 11 tyrrhenien



! Alboran                      
        if(lon_t(i,j)<=-1.*deg2rad.and.lon_t(i,j)>-6.*deg2rad)then

!      print*,'Passe ici dans InitNut',lon_t(i,j)/deg2rad,lat_t(i,j)/deg2rad
!      print*,'-depth_t(i,j,k)',k,-depth_t(i,j,k),Prof_oxyg_EAST(2,1   )
!      print*,'2',InitValue_EAST(2,1),InitValue_EAST(2,nmax_oxy_alb)
!      stop'dans InitNut'

             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(1,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(1,nmax_oxy_alb))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(1,nmax_oxy_alb)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(1,K_pr)+depth_t(i,j,k))*InitValue_WEST(1,K_pr-1) &
             +(-Prof_oxyg_WEST(1,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(1,K_pr) &
             )/(Prof_oxyg_WEST(1,K_pr  )-Prof_oxyg_WEST(1,K_pr-1))

!       print*,'3',InitValue_WEST(1,K_pr),InitValue_WEST(1,K_pr-1),InitValue_WEST(1,K_pr)

           endif 

        if (depth_t(i,j,k)<-300.and.depth_t(i,j,k)>-500)  &   !ajout le 11/12/17
            bio_t(i,j,k,iOXYGEN)=165


        endif

! Algero-Provençal                      
        if(lon_t(i,j)<=10.9*deg2rad.and.lon_t(i,j)>-1.*deg2rad)then

             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(2,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(2,nmax_oxy_alg))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(2,nmax_oxy_alg)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(2,K_pr)+depth_t(i,j,k))*InitValue_WEST(2,K_pr-1) &
             +(-Prof_oxyg_WEST(2,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(2,K_pr) &
             )/(Prof_oxyg_WEST(2,K_pr  )-Prof_oxyg_WEST(2,K_pr-1))
           endif


        if (depth_t(i,j,k).le.-300.and.depth_t(i,j,k).ge.-500)  &   !ajout le 11/12/17
            bio_t(i,j,k,iOXYGEN)=169.

        endif

! Tyrhennienne                       
        if(lon_t(i,j)<=16.4*deg2rad.and.lon_t(i,j)>9.05*deg2rad.and. &
           lat_t(i,j)<=42.83*deg2rad)then

!             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(4,1   ))then
!               bio_t(i,j,k,iOxygen) = InitValue_WEST(4,1)
!
!             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(4,nmax_oxy_tyr))then
!               bio_t(i,j,k,iOxygen) = InitValue_WEST(4,nmax_oxy_tyr)
!
!             ELSE
!
!              K_pr=1
!              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(4,K_pr))
!               K_pr=K_pr+1
!              enddo
!
!              bio_t(i,j,k,iOxygen)=(                              &
!              (Prof_oxyg_WEST(4,K_pr)+depth_t(i,j,k))*InitValue_WEST(4,K_pr-1) &
!             +(-Prof_oxyg_WEST(4,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(4,K_pr) &
!             )/(Prof_oxyg_WEST(4,K_pr  )-Prof_oxyg_WEST(4,K_pr-1))
!           endif


             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(2,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(2,nmax_oxy_alg))then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(2,nmax_oxy_alg)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(2,K_pr)+depth_t(i,j,k))*InitValue_WEST(2,K_pr-1) &
             +(-Prof_oxyg_WEST(2,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(2,K_pr) &
             )/(Prof_oxyg_WEST(2,K_pr  )-Prof_oxyg_WEST(2,K_pr-1))
           endif


        endif

! Centrale est                       
        if(lon_t(i,j)<=16.2*deg2rad.and.lon_t(i,j)>9.05*deg2rad.and. &
           lat_t(i,j)<=37.83*deg2rad)then

                     if    (-depth_t(i,j,k)<=Prof_oxyg_EAST(3,1   )) then
                            bio_t(i,j,k,iOxygen) =InitValue_EAST(3,1)
                     ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EAST(3,nmax_oxy_cen)) then
                            bio_t(i,j,k,iOxygen) =InitValue_EAST(3,nmax_oxy_cen)
                     ELSE
                     K_pr=1
                     do while(-depth_t(i,j,k)>Prof_oxyg_EAST(3,K_pr))
                        K_pr=K_pr+1
                     enddo
                     bio_t(i,j,k,iOxygen)=( &
                     (Prof_oxyg_EAST(3,K_pr)+depth_t(i,j,k))*InitValue_EAST(3,K_pr-1) &
                    +(-Prof_oxyg_EAST(3,K_pr-1)-depth_t(i,j,k))*InitValue_EAST(3,K_pr) &
                        )/(Prof_oxyg_EAST(3,K_pr)-Prof_oxyg_EAST(3,K_pr-1))
                     endif

        endif        


! Ionienne                       
        if(lon_t(i,j)<=27.*deg2rad.and.lon_t(i,j)>16.2*deg2rad)then
               
             if    (-depth_t(i,j,k)<=Prof_oxyg_EAST(2,1   ))then
               bio_t(i,j,k,iOxygen) = InitValue_EAST(2,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EAST(2,nmax_oxy_ion))then
               bio_t(i,j,k,iOxygen) = InitValue_EAST(2,nmax_oxy_ion)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EAST(2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EAST(2,K_pr)+depth_t(i,j,k))*InitValue_EAST(2,K_pr-1) &
             +(-Prof_oxyg_EAST(2,K_pr-1)-depth_t(i,j,k))*InitValue_EAST(2,K_pr) &
             )/(Prof_oxyg_EAST(2,K_pr  )-Prof_oxyg_EAST(2,K_pr-1))
           endif

        endif

! Levantine                      
        if(lon_t(i,j)>27.*deg2rad)then
               
             if    (-depth_t(i,j,k)<=Prof_oxyg_EAST(1,1   )) then
               bio_t(i,j,k,iOxygen) = InitValue_EAST(1,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_EAST(1,nmax_oxy_lev)) then
               bio_t(i,j,k,iOxygen) = InitValue_EAST(1,nmax_oxy_lev)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_EAST(1,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_EAST(1,K_pr)+depth_t(i,j,k))*InitValue_EAST(1,K_pr-1) &
             +(-Prof_oxyg_EAST(1,K_pr-1)-depth_t(i,j,k))*InitValue_EAST(1,K_pr) &
             )/(Prof_oxyg_EAST(1,K_pr  )-Prof_oxyg_EAST(1,K_pr-1))
           endif

        endif

! deep NV                       
        if(lon_t(i,j)<10.*deg2rad.and. &
           lat_t(i,j)>40.*deg2rad.and.h_w(i,j)>2000)then

             if    (-depth_t(i,j,k)<=Prof_oxyg_WEST(3,1   )) then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(3,1)

             ELSEif(-depth_t(i,j,k)>=Prof_oxyg_WEST(3,nmax_oxy_gol)) then
               bio_t(i,j,k,iOxygen) = InitValue_WEST(3,nmax_oxy_gol)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof_oxyg_WEST(3,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,iOxygen)=(                              &
              (Prof_oxyg_WEST(3,K_pr)+depth_t(i,j,k))*InitValue_WEST(3,K_pr-1) &
             +(-Prof_oxyg_WEST(3,K_pr-1)-depth_t(i,j,k))*InitValue_WEST(3,K_pr) &
             )/(Prof_oxyg_WEST(3,K_pr  )-Prof_oxyg_WEST(3,K_pr-1))
           endif

        if (depth_t(i,j,k).le.-300.and.depth_t(i,j,k).ge.-500)  & !ajout le 11/12/17
            bio_t(i,j,k,iOXYGEN)=173.

        endif



          bio_t(i,j,k,iOxygen)       =bio_t(i,j,k,iOxygen)  &
                            *(rhp_t(i,j,k)+rho) /1000.


      endif ! mask 
      enddo
      enddo
      enddo


       if(par%rank==0) write(6,*)'Passe par InitDOxy'
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!09-05-10
#endif

      END SUBROUTINE InitDOxy
