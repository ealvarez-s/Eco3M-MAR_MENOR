










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
! 21/03/2013: modification dans les initialisations bio               *
! 13/10/2017: ajout Caroline condition pour corriger valeurs négatives*      
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
      ! modif polygone et initialisation Alex, modifié pour Joelle
      integer,parameter :: nmaxprof=300, npoly=10
      integer :: npoints_poly, n, polygon_m(mbio1:mbio2,nbio1:nbio2),nmaxprofk(11,14)

      double precision :: x,y, r1, r2, deltap
      double precision,dimension(:),allocatable :: xp,yp

      logical :: ldinmesh
      character(100) :: name_poly

      integer ChoiceInitPelagic,NumVar(20),LDir,K_pr  &
             ,InitFemme,ifemme,JFemme,nmaxlpc,nmax,nmax1,nmax2,nmax3    &
             ,isc2160,jsc2160,isw2060,jsw2060,nmax4,nmax5,nmax6         &
             ,isc2240,jsc2240,ij,InitOC

      character DATE*6,DirInitFile*90,Fic(20)*80      &
               ,NameVar(vbmax)*9,FicInitVar*80        &
               ,FicLPC_Pho*80,FicLPC_Nit*80


      double precision InitValue(20,3000),Prof(20,3000),Prof_oxyg(20,3000),Prof_oxyg_WEST(20,3000)       &
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
                      ,InitValue_EST_DOC(20,3000)   &
                      ,InitValue_WEST_DOC(20,3000)   &
                      ,InitValue_EST_bact(20,3000)   &
                      ,InitValue_WEST_bact(20,3000)   &
!                      ,Prof_EST_doc(20,3000),Prof_WEST_doc(20,3000)       &
                      ,Init_DOC_Refr,Init_DON_Refr &
!                      ,Prof_EST_bact(20,3000),Prof_WEST_bact(20,3000)       &
                      ! New stuff 
                      ! valeur de 1 à 9 pour les sous-régions: ALB, ALG, GOL,
                      ! TYR, CMED, ION, ADR, AEG, LEV
                      ,InitVal_NIT(npoly,nmaxprof)   &
                      ,InitVal_PHO(npoly,nmaxprof)   &
                      ,InitVal_SIL(npoly,nmaxprof)   &
                      ,Init_PRF(npoly,nmaxprof)        &
                      ,InitVal_DIC(npoly,nmaxprof)   &
                      ,Init_PRF_DIC(npoly,nmaxprof)        &
                      ,InitVal_ALK(npoly,nmaxprof)   &
                      ,Init_PRF_ALK(npoly,nmaxprof)        &
                      ,InitVal_OXY(npoly,nmaxprof)        &
                      ,InitVal_bact(npoly,nmaxprof)   &
                      ,Init_PRF_bact(npoly,nmaxprof)        &
                      ,InitVal_DOC(npoly,nmaxprof)   &
                      ,Init_PRF_DOC(npoly,nmaxprof)        &
                      ,InitVal_DON(npoly,nmaxprof)   &
                      ,Init_PRF_DON(npoly,nmaxprof)        &
                      ,InitVal_DOP(npoly,nmaxprof)   &
                      ,Init_PRF_DOP(npoly,nmaxprof)        &
                      ,InitVal_CHL(npoly,nmaxprof)   &
                      ,InitVal_MICRO(npoly,nmaxprof)   &
                      ,InitVal_NANO(npoly,nmaxprof)   &
                      ,InitVal_PICO(npoly,nmaxprof)   &
                      ,Init_PRF_PHO(npoly,nmaxprof)        &
                      ,Init_PRF_SIL(npoly,nmaxprof)  &
                      ,Init_PRF_CHL(npoly,nmaxprof)        &
                      ,Init_PRF_OXY(npoly,nmaxprof)  &
                      ,Init_PRF_MICRO(npoly,nmaxprof)  & 
                      ! New stuff
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
           InitPH3D2(iglb+2,jglb+2,kmax),InitMODC(iglb+2,jglb+2,kmax), &
           perchl(iglb+2,jglb+2,kmax,3)            
 
       integer l,m,nmax10y(14)
       real Prof10y(14,2511),var10y(14,2511,10)
       character*60 filename


       InitOC=0
       Init_DOC_Refr=37.5 !34
       Init_DON_Refr=3.1 !2.5


!      open(unit=10,file='../../../MERCATOFF/InitMODCrestart')
!!      print*,'dim',kmax,jglb+2,iglb+2
!      do k=1,kmax
!      do j=1,jglb+2
!      do i=1,iglb+2
!      read(10,*)InitMODC(i,j,k)
!      enddo
!      enddo
!      enddo 


!---------------------------------------------------------------------*


!=====================================================================*
! 3. Initialisation of the variables                                  *
! BEGINNING                                                           *
!=====================================================================*

! State variables
!-----------------

! claude : Caroline a moins de fichiers et moins de NumVar


      nmax1=68
      nmax2=68
      nmax3=68
      nmax4=21
      nmax5=22
!     nmax6=183
      nmax6=2311

!     open(unit=2,file='/tmpdir/culses/bassin2/MERCATOFF/DyfamedmeanpH.txt')
      open(unit=2,file='/tmpdir/jhabib/simulations/initdyfamed/DyfamedmeanpH.txt')
      do k=1,nmax5
      read(2,*)ProfPH(1,k),InitPH(1,k)
      enddo
      close(2)
       

! to get initial temperature and salinity fields  
           if(ioffline==2) call offline_inout(2)   
!$ water density:                                                      !14/11/04
!      if(eqs_state1.eq.0)call density(0) !lineaire
!      if(eqs_state1.eq.1)call density(3) !Non line sans pr
       call equation_of_state('potential density',1)

! nutrients


!!!!!! !!!!!!!!!!!!!!!!!! LOADING PROFILES FROM OBSERVATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Lecture des profils de PROFILNUTS : Joelle ! 2020-04-07
! ALB Modif Alex 06/04/2018 *nit.txt remplace par *nit2.txt, pareil pour phos
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Albo_nit.dat') !max + faible que pour Alex
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALB_nit.dat') 
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALB_nit.dat')
      read(3,*)nmaxprofk(1,1)
      do k=1,nmaxprofk(1,1)
      read(3,*)Init_PRF(1,k),InitVal_NIT(1,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Albo_phos.dat') !max + faible que pour Alex
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALB_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALB_phos.dat')
      read(3,*)nmaxprofk(1,2)
      do k=1,nmaxprofk(1,2)
      read(3,*)Init_PRF_PHO(1,k),InitVal_PHO(1,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Albo_sil.dat') ! diminution vers le fond à voir si réaliste
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALB_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALB_sil.dat')
      read(3,*)nmaxprofk(1,3)
      do k=1,nmaxprofk(1,3)
      read(3,*)Init_PRF_SIL(1,k),InitVal_SIL(1,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/moy_alboran.dat') 
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_alboran.dat')
      read(3,*)nmaxprofk(1,4)
      do k=1,nmaxprofk(1,4)
      read(3,*)InitVal_DIC(1,k),Init_PRF_DIC(1,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_alboran.dat')
      read(3,*)nmaxprofk(1,5)
      do k=1,nmaxprofk(1,5)
      read(3,*)InitVal_ALK(1,k),Init_PRF_ALK(1,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALB_oxy.dat')
      read(3,*)nmaxprofk(1,6)
      do k=1,nmaxprofk(1,6)
      read(3,*)Init_PRF_OXY(1,k),InitVal_OXY(1,k)
      enddo
      close(3)

! ALG
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Alg_nit.dat') ! + faible que dans Alex
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALG_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALG_nit.dat')
      read(3,*)nmaxprofk(2,1)
      do k=1,nmaxprofk(2,1)
      read(3,*)Init_PRF(2,k),InitVal_NIT(2,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Alg_phos.dat') ! + faible que dans Alex
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALG_phos.dat')
    open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALG_phos.dat')
      read(3,*)nmaxprofk(2,2)
      do k=1,nmaxprofk(2,2)
      read(3,*)Init_PRF_PHO(2,k),InitVal_PHO(2,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Alg_sil.dat') ! pas de max à 9 à 1000m comme dans Alex, en meilleur accord avec Tanhua?
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALG_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALG_sil.dat')
      read(3,*)nmaxprofk(2,3)
      do k=1,nmaxprofk(2,3)
      read(3,*)Init_PRF_SIL(2,k),InitVal_SIL(2,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_algerian.dat') 
      read(3,*)nmaxprofk(2,4)
      do k=1,nmaxprofk(2,4)
      read(3,*)InitVal_DIC(2,k),Init_PRF_DIC(2,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_algerian.dat')
      read(3,*)nmaxprofk(2,5)
      do k=1,nmaxprofk(2,5)
      read(3,*)InitVal_ALK(2,k),Init_PRF_ALK(2,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALG_oxy.dat')
      read(3,*)nmaxprofk(2,6)
      do k=1,nmaxprofk(2,6)
      read(3,*)Init_PRF_OXY(2,k),InitVal_OXY(2,k)
      enddo
      close(3)

! GOL
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Gdl_nit.dat') ! pas de max dans les eaux intermédiaires
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/GDL_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/GDL_nit.dat')
      read(3,*)nmaxprofk(3,1)
      do k=1,nmaxprofk(3,1)
      read(3,*)Init_PRF(3,k),InitVal_NIT(3,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Gdl_phos.dat') ! pas de max dans les eaux intermédiaires (dans Alex le pic était à 0.5!)
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/GDL_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/GDL_phos.dat')
      read(3,*)nmaxprofk(3,2)
      do k=1,nmaxprofk(3,2)
      read(3,*)Init_PRF_PHO(3,k),InitVal_PHO(3,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Gdl_sil.dat') ! idem que nitrate
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/GDL_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/GDL_sil.dat')
      read(3,*)nmaxprofk(3,3)
      do k=1,nmaxprofk(3,3)
      read(3,*)Init_PRF_SIL(3,k),InitVal_SIL(3,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_gdl.dat')  
      read(3,*)nmaxprofk(3,4)
      do k=1,nmaxprofk(3,4)
      read(3,*)InitVal_DIC(3,k),Init_PRF_DIC(3,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_gdl.dat')
      read(3,*)nmaxprofk(3,5)
      do k=1,nmaxprofk(3,5)
      read(3,*)InitVal_ALK(3,k),Init_PRF_ALK(3,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/GDL_oxy.dat')
      read(3,*)nmaxprofk(3,6)
      do k=1,nmaxprofk(3,6)
      read(3,*)Init_PRF_OXY(3,k),InitVal_OXY(3,k)
      enddo
      close(3)
! TYR
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Tyr_nit.dat') ! un peu plus faible que dans Alex 
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/TYRR_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/TYRR_nit.dat')     
      read(3,*)nmaxprofk(4,1)
      do k=1,nmaxprofk(4,1)
      read(3,*)Init_PRF(4,k),InitVal_NIT(4,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Tyr_phos.dat') ! un peu plus faible que dans Alex
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/TYRR_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/TYRR_phos.dat')
      read(3,*)nmaxprofk(4,2)
      do k=1,nmaxprofk(4,2)
      read(3,*)Init_PRF_PHO(4,k),InitVal_PHO(4,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Tyr_sil.dat') ! un peu plus faible que dans Alex
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/TYRR_sil.dat')
     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/TYRR_sil.dat')
      read(3,*)nmaxprofk(4,3)
      do k=1,nmaxprofk(4,3)
      read(3,*)Init_PRF_SIL(4,k),InitVal_SIL(4,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_tyrr.dat')     
      read(3,*)nmaxprofk(4,4)
      do k=1,nmaxprofk(4,4)
      read(3,*)InitVal_DIC(4,k),Init_PRF_DIC(4,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_tyrr.dat')
      read(3,*)nmaxprofk(4,5)
      do k=1,nmaxprofk(4,5)
      read(3,*)InitVal_ALK(4,k),Init_PRF_ALK(4,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/TYRR_oxy.dat')
      read(3,*)nmaxprofk(4,6)
      do k=1,nmaxprofk(4,6)
      read(3,*)Init_PRF_OXY(4,k),InitVal_OXY(4,k)
      enddo
      close(3)

! ALGEROBAL
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Algerob_nit.dat') ! assez similaire à GdL
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALGEROBAL_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALGEROB_nit.dat')
      read(3,*)nmaxprofk(5,1)
      do k=1,nmaxprofk(5,1)
      read(3,*)Init_PRF(5,k),InitVal_NIT(5,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Algerob_phos.dat')
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALGEROBAL_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALGEROB_phos.dat')
      read(3,*)nmaxprofk(5,2)
      do k=1,nmaxprofk(5,2)
      read(3,*)Init_PRF_PHO(5,k),InitVal_PHO(5,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Algerob_sil.dat')
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALGEROBAL_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ALGEROB_sil.dat')  
      read(3,*)nmaxprofk(5,3)
      do k=1,nmaxprofk(5,3)
      read(3,*)Init_PRF_SIL(5,k),InitVal_SIL(5,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_algerobal.dat')
      read(3,*)nmaxprofk(5,4)
      do k=1,nmaxprofk(5,4)
      read(3,*)InitVal_DIC(5,k),Init_PRF_DIC(5,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_algerobal.dat')
      read(3,*)nmaxprofk(5,5)
      do k=1,nmaxprofk(5,5)
      read(3,*)InitVal_ALK(5,k),Init_PRF_ALK(5,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ALGEROBAL_oxy.dat')
      read(3,*)nmaxprofk(5,6)
      do k=1,nmaxprofk(5,6)
      read(3,*)Init_PRF_OXY(5,k),InitVal_OXY(5,k)
      enddo
      close(3)
! ION
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Ion_nit.dat') ! max plus profond que pour Alex et Faycal, ce qui est en accord avec Tanhua et al 2013
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ION_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ION_nit.dat')
      read(3,*)nmaxprofk(6,1)
      do k=1,nmaxprofk(6,1)
      read(3,*)Init_PRF(6,k),InitVal_NIT(6,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Ion_phos.dat') ! idem que pour nitrate, la valeur du max en meilleur accord avec Tanhua et al 2013
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ION_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ION_phos.dat')
      read(3,*)nmaxprofk(6,2)
      do k=1,nmaxprofk(6,2)
      read(3,*)Init_PRF_PHO(6,k),InitVal_PHO(6,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Ion_sil.dat') ! semble mieux(moins fort) dans la couche de surface qu'A et F, dans Faycal concentration vraiment tres faible au fond
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ION_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ION_sil.dat')
      read(3,*)nmaxprofk(6,3)
      do k=1,nmaxprofk(6,3)
      read(3,*)Init_PRF_SIL(6,k),InitVal_SIL(6,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_ionien.dat')
      read(3,*)nmaxprofk(6,4)
      do k=1,nmaxprofk(6,4)
      read(3,*)InitVal_DIC(6,k),Init_PRF_DIC(6,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_ionien.dat')
      read(3,*)nmaxprofk(6,5)
      do k=1,nmaxprofk(6,5)
      read(3,*)InitVal_ALK(6,k),Init_PRF_ALK(6,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ION_oxy.dat')
      read(3,*)nmaxprofk(6,6)
      do k=1,nmaxprofk(6,6)
      read(3,*)Init_PRF_OXY(6,k),InitVal_OXY(6,k)
      enddo
      close(3)
! ADR
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Adria_nit.dat') ! concentration plus faible que dans Alex et F., en coherence avec Tanhua
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ADRIA_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ADRIA_nit.dat')
      read(3,*)nmaxprofk(7,1)
      do k=1,nmaxprofk(7,1)
      read(3,*)Init_PRF(7,k),InitVal_NIT(7,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Adria_phos.dat') ! idem que nit
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ADRIA_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ADRIA_phos.dat')
      read(3,*)nmaxprofk(7,2)
      do k=1,nmaxprofk(7,2)
      read(3,*)Init_PRF_PHO(7,k),InitVal_PHO(7,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Adria_sil.dat') ! idem que nit
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ADRIA_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/ADRIA_sil.dat')
      read(3,*)nmaxprofk(7,3)
      do k=1,nmaxprofk(7,3)
      read(3,*)Init_PRF_SIL(7,k),InitVal_SIL(7,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_adria.dat')
      read(3,*)nmaxprofk(7,4)
      do k=1,nmaxprofk(7,4)
      read(3,*)InitVal_DIC(7,k),Init_PRF_DIC(7,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_adria.dat')
      read(3,*)nmaxprofk(7,5)
      do k=1,nmaxprofk(7,5)
      read(3,*)InitVal_ALK(7,k),Init_PRF_ALK(7,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/ADRIA_oxy.dat')
      read(3,*)nmaxprofk(7,6)
      do k=1,nmaxprofk(7,6)
      read(3,*)Init_PRF_OXY(7,k),InitVal_OXY(7,k)
      enddo
      close(3)

! AEG
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Egee_nit.dat') ! zone pas profonde, plus faible
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/EGEE_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/EGEE_nit.dat')
      read(3,*)nmaxprofk(8,1)
      do k=1,nmaxprofk(8,1)
      read(3,*)Init_PRF(8,k),InitVal_NIT(8,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Egee_phos.dat') ! zone pas profonde, plus faible
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/EGEE_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/EGEE_phos.dat')
      read(3,*)nmaxprofk(8,2)
      do k=1,nmaxprofk(8,2)
      read(3,*)Init_PRF_PHO(8,k),InitVal_PHO(8,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Egee_sil.dat') ! zone pas profonde, comme Faycal
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/EGEE_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/EGEE_sil.dat')      
      read(3,*)nmaxprofk(8,3)
      do k=1,nmaxprofk(8,3)
      read(3,*)Init_PRF_SIL(8,k),InitVal_SIL(8,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_egee.dat')
      read(3,*)nmaxprofk(8,4)
      do k=1,nmaxprofk(8,4)
      read(3,*)InitVal_DIC(8,k),Init_PRF_DIC(8,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_egee.dat')
      read(3,*)nmaxprofk(8,5)
      do k=1,nmaxprofk(8,5)
      read(3,*)InitVal_ALK(8,k),Init_PRF_ALK(8,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/EGEE_oxy.dat')
      read(3,*)nmaxprofk(8,6)
      do k=1,nmaxprofk(8,6)
      read(3,*)Init_PRF_OXY(8,k),InitVal_OXY(8,k)
      enddo
      close(3)
! LEV
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Lev_nit.dat') ! max + prononcé et semble mieux localisé que dans Alex
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LEV_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/LEV_nit.dat')
      read(3,*)nmaxprofk(9,1)
      do k=1,nmaxprofk(9,1)
      read(3,*)Init_PRF(9,k),InitVal_NIT(9,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Lev_phos.dat') ! idem que nitrate
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LEV_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/LEV_phos.dat')
      read(3,*)nmaxprofk(9,2)
      do k=1,nmaxprofk(9,2)
      read(3,*)Init_PRF_PHO(9,k),InitVal_PHO(9,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Lev_sil.dat') 
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LEV_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/LEV_sil.dat')
      read(3,*)nmaxprofk(9,3)
      do k=1,nmaxprofk(9,3)
      read(3,*)Init_PRF_SIL(9,k),InitVal_SIL(9,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_levantin.dat')
      read(3,*)nmaxprofk(9,4)
      do k=1,nmaxprofk(9,4)
      read(3,*)InitVal_DIC(9,k),Init_PRF_DIC(9,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_levantin.dat')
      read(3,*)nmaxprofk(9,5)
      do k=1,nmaxprofk(9,5)
      read(3,*)InitVal_ALK(9,k),Init_PRF_ALK(9,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LEV_oxy2.dat')
      read(3,*)nmaxprofk(9,6)
      do k=1,nmaxprofk(9,6)
      read(3,*)Init_PRF_OXY(9,k),InitVal_OXY(9,k)
      enddo
      close(3)
! CRETE 
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Crete_nit.dat') ! nouvelle région, plus élevée que Egée d'Alex et Faycal
!      read(3,*)nmaxprofk(10,1)
!      do k=1,nmaxprofk(10,1)
!      read(3,*)Init_PRF(10,k),InitVal_NIT(10,k)
!      enddo
!      close(3)
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Crete_phos.dat') ! 
!      read(3,*)nmaxprofk(10,2)
!      do k=1,nmaxprofk(10,2)
!      read(3,*)Init_PRF(10,k),InitVal_PHO(10,k)
!      enddo
!      close(3)
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Crete_sil.dat')
!      read(3,*)nmaxprofk(10,3)
!      do k=1,nmaxprofk(10,3)
!      read(3,*)Init_PRF(10,k),InitVal_SIL(10,k)
!      enddo
!      close(3)
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/moy_crete.dat')
!      read(3,*)nmaxprofk(10,4)
!      do k=1,nmaxprofk(10,4)
!      read(3,*)InitVal_DIC(10,k),Init_PRF_DIC(10,k)
!      enddo
!      close(3)
!      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/caro_moy_Talk_alboran.dat')
!      read(3,*)nmaxprofk(1,5)
!      do k=1,nmaxprofk(1,5)
!      read(3,*)InitVal_ALK(1,k),Init_PRF_ALK(1,k)
!      enddo
!      close(3)

! Ligurie
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Lig_nit.dat') ! nouvelle région, similaire à GdL ex et Faycal
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LIG_nit.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/LIG_nit.dat')
      read(3,*)nmaxprofk(10,1)
      do k=1,nmaxprofk(10,1)
      read(3,*)Init_PRF(10,k),InitVal_NIT(10,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Lig_phos.dat') !
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LIG_phos.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/LIG_phos.dat')
      read(3,*)nmaxprofk(10,2)
      do k=1,nmaxprofk(10,2)
      read(3,*)Init_PRF_PHO(10,k),InitVal_PHO(10,k)
      enddo
      close(3)
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/Lig_sil.dat')
!     open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LIG_sil.dat')
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_022021/LIG_sil.dat')
      read(3,*)nmaxprofk(10,3)
      do k=1,nmaxprofk(10,3)
      read(3,*)Init_PRF_SIL(10,k),InitVal_SIL(10,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_DIC_ligure.dat')
      read(3,*)nmaxprofk(10,4)
      do k=1,nmaxprofk(10,4)
      read(3,*)InitVal_DIC(10,k),Init_PRF_DIC(10,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS/caro_moy_Talk_ligure.dat')
      read(3,*)nmaxprofk(10,5)
      do k=1,nmaxprofk(10,5)
      read(3,*)InitVal_ALK(10,k),Init_PRF_ALK(10,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/jhabib/simulations/joelle_profil/PROFILNUTS_092020/LIG_oxy.dat')
      read(3,*)nmaxprofk(10,6)
      do k=1,nmaxprofk(10,6)
      read(3,*)Init_PRF_OXY(10,k),InitVal_OXY(10,k)
      enddo
      close(3)
!      Init_PRF(:,:)=-Init_PRF(:,:)

! Initialisation et lecture des polygones
! Attention à l'ordre des polygones dans le fichier
! ici l'ordre est :
! ALB, ALG, GOL, TYR, ALGEROB, ION, ADR, AEG, LEV, Crete, Lig

      polygon_m(:,:)=0
!     open(unit=11,file='/tmpdir/alex/BASSIN_3D/POLYGONE/liste_poly')
!     open(unit=11,file='/tmpdir/jhabib/simulations/joelle_profil/POLYGONE_20/liste_poly')
      open(unit=11,file='/tmpdir/jhabib/simulations/joelle_profil/POLYGONE_202008/liste_poly')
!  nouvel ordre : ALB, ALG, GOL, TYR, ALGEROB, ION, ADR, AEG, LEV,  Lig (pareil
!  qu'avant avec le Crete en moins)
      do n=1,npoly

        read(11,'(a)')name_poly
        open(unit=12,file=trim(name_poly))
        read(12,*)npoints_poly
        if(allocated(xp))then
         deallocate(xp)
         deallocate(yp)
        endif
        allocate(xp(npoints_poly))
        allocate(yp(npoints_poly))

        do i=1,npoints_poly
        read(12,*)xp(i),yp(i)
        enddo

      ! On repère les points correspondants aux polygones                                                                   
        do i1=mbio1,mbio2
        do j1=nbio1,nbio2
          if (mask_t(i1,j1,kmax)==1.and.polygon_m(i1,j1)==0) then
            x=lon_t(i1,j1)*rad2deg
            y=lat_t(i1,j1)*rad2deg
            ldinmesh=.false.
            j=npoints_poly
            do i=1,npoints_poly
              if ((((YP(I)<=Y).and.(Y<YP(J))).or.((YP(J)<=Y).and.(Y<YP(I)))) &
               .and.(X<(XP(J)-XP(I))*(Y-YP(I))/(YP(J)-YP(I))+XP(I))) then
                if(ldinmesh)then
                  ldinmesh=.false.
                else
                  ldinmesh=.true.
                endif
              endif
              j=i
            enddo


            if(ldinmesh)then
              polygon_m(i1,j1)=n
            endif
            ! polygon_m(i1,j1) contient alors le numero du polygone et donc du
            ! profil a interpoler

          endif
        enddo
        enddo
        close(12)
      enddo
      close(11)

!!!!!! !!!!!!!!!!!!!!!!!! LOADING PROFILES FROM 10-year simulations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!      do l=1,14
!!     nmax10y=2511
!
!      if(l<10) then
!!      write (filename,"(A50,I1,A4)")'/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitOC2_',l,'.txt'
!!      else
!!      write(filename,"(A50,I2,A4)")'/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitOC2_',l,'.txt'
!      write(filename,"(A52,I1,A4)")'/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitOCM8_',l,'.txt'
!      else
!      write(filename,"(A52,I2,A4)")'/tmpdir/jhabib/simulations/SIMED/PROFILNUTS/InitOCM8_',l,'.txt'
!      endif
!
!      open(unit=3,file=filename)
!      read(3,*)nmax10y(l)
!!      print*,'l',l,nmax10y(l)
!      do k=1,nmax10y(l)
!      read(3,*)Prof10y(l,k),(var10y(l,k,m),m=1,10)
!! 13/10 ajout Caroline correction valeurs négatives      
!      do m=1,10
!      if(var10y(l,k,m)<0) var10y(l,k,m)=0.
!      enddo
!      enddo
!      close(3)
!
!      enddo



!!!!!!!!!!!!!!!!!!!! END PROFILE LOADING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      do i=mbio1,mbio2
      do j=nbio1,nbio2
      do k=kmin_w(i,j),kmax


       if (mask_t(i,j,k)==1) then

       ! lecture et interpolation des profils pour nitrate, silice phosphate PUIS
       ! Initialisation de l'oxygene, zoo, ALK, etc

       ! interpolation des profils de Nitrate, Phosphate et Silice

         if(polygon_m(i,j).eq.0) &
         print*,'polygone',i,j,i+par%timax(1),j+par%tjmax(1)


         if(nmaxprofk(polygon_m(i,j),3).eq.0) &
         print*,'nmaxprofk',polygon_m(i,j)

         if (depth_t(i,j,k)<=Init_PRF(polygon_m(i,j),1)) then
           bio_t(i,j,k,iNitrate)   = InitVal_NIT(polygon_m(i,j),1)
           
         elseif(depth_t(i,j,k)>=Init_PRF(polygon_m(i,j),nmaxprofk(polygon_m(i,j),1)))then
           bio_t(i,j,k,iNitrate)   = InitVal_NIT(polygon_m(i,j),nmaxprofk(polygon_m(i,j),1))
         else
  
           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF(polygon_m(i,j),K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF(polygon_m(i,j),K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF(polygon_m(i,j),K_pr)
           deltap=Init_PRF(polygon_m(i,j),K_pr)-Init_PRF(polygon_m(i,j),K_pr-1) 

           bio_t(i,j,k,iNitrate)  =(r1*InitVal_NIT(polygon_m(i,j),K_pr  )        &
                                   +r2*InitVal_NIT(polygon_m(i,j),K_pr-1))/deltap
          endif 


! Pour le phosphate
         if (depth_t(i,j,k)<=Init_PRF_PHO(polygon_m(i,j),1)) then
           bio_t(i,j,k,iPhosphate)        = InitVal_PHO(polygon_m(i,j),1)
         elseif(depth_t(i,j,k)>=Init_PRF_PHO(polygon_m(i,j),nmaxprofk(polygon_m(i,j),2)))then
           bio_t(i,j,k,iPhosphate)    = InitVal_PHO(polygon_m(i,j),nmaxprofk(polygon_m(i,j),2))
         else

           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF_PHO(polygon_m(i,j),K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF_PHO(polygon_m(i,j),K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF_PHO(polygon_m(i,j),K_pr)
           deltap=Init_PRF_PHO(polygon_m(i,j),K_pr)-Init_PRF_PHO(polygon_m(i,j),K_pr-1)

           bio_t(i,j,k,iPhosphate)  =(r1*InitVal_PHO(polygon_m(i,j),K_pr  ) &
                                   +r2*InitVal_PHO(polygon_m(i,j),K_pr-1))/deltap

          endif

! Pour le silicate
         if (depth_t(i,j,k)<=Init_PRF_SIL(polygon_m(i,j),1)) then
           bio_t(i,j,k,iSilice)        = InitVal_SIL(polygon_m(i,j),1)
         elseif(depth_t(i,j,k)>=Init_PRF_SIL(polygon_m(i,j),nmaxprofk(polygon_m(i,j),3)))then
           bio_t(i,j,k,iSilice)    = InitVal_SIL(polygon_m(i,j),nmaxprofk(polygon_m(i,j),3))
         else

           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF_SIL(polygon_m(i,j),K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF_SIL(polygon_m(i,j),K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF_SIL(polygon_m(i,j),K_pr)
           deltap=Init_PRF_SIL(polygon_m(i,j),K_pr)-Init_PRF_SIL(polygon_m(i,j),K_pr-1)

           bio_t(i,j,k,iSilice)  =(r1*InitVal_SIL(polygon_m(i,j),K_pr  ) &
                                   +r2*InitVal_SIL(polygon_m(i,j),K_pr-1))/deltap

          endif

! Pour le oxygene
         if (depth_t(i,j,k)<=Init_PRF_OXY(polygon_m(i,j),1)) then
           bio_t(i,j,k,ioxygen)        = InitVal_OXY(polygon_m(i,j),1)
         elseif(depth_t(i,j,k)>=Init_PRF_OXY(polygon_m(i,j),nmaxprofk(polygon_m(i,j),6)))then
           bio_t(i,j,k,ioxygen)    = InitVal_oxy(polygon_m(i,j),nmaxprofk(polygon_m(i,j),6))
         else
         
           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF_OXY(polygon_m(i,j),K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF_OXY(polygon_m(i,j),K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF_OXY(polygon_m(i,j),K_pr)
           deltap=Init_PRF_OXY(polygon_m(i,j),K_pr)-Init_PRF_OXY(polygon_m(i,j),K_pr-1)

           bio_t(i,j,k,ioxygen)  =(r1*InitVal_OXY(polygon_m(i,j),K_pr  ) &
                                   +r2*InitVal_OXY(polygon_m(i,j),K_pr-1))/deltap

          endif

! Pour le DIC
         if (depth_t(i,j,k)<=Init_PRF_DIC(polygon_m(i,j),1)) then
           bio_t(i,j,k,iDIC)        = InitVal_DIC(polygon_m(i,j),1)
         elseif(depth_t(i,j,k)>=Init_PRF_DIC(polygon_m(i,j),nmaxprofk(polygon_m(i,j),4)))then
           bio_t(i,j,k,iDIC)    = InitVal_DIC(polygon_m(i,j),nmaxprofk(polygon_m(i,j),4))
         else

           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF_DIC(polygon_m(i,j),K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF_DIC(polygon_m(i,j),K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF_DIC(polygon_m(i,j),K_pr)
           deltap=Init_PRF_DIC(polygon_m(i,j),K_pr)-Init_PRF_DIC(polygon_m(i,j),K_pr-1)

           bio_t(i,j,k,iDIC)  =(r1*InitVal_DIC(polygon_m(i,j),K_pr  ) &
                                   +r2*InitVal_DIC(polygon_m(i,j),K_pr-1))/deltap



          endif

! Pour l'alcalinité
         if (depth_t(i,j,k)<=Init_PRF_ALK(polygon_m(i,j),1)) then
           bio_t(i,j,k,iALKALINITY)        = InitVal_ALK(polygon_m(i,j),1)
         elseif(depth_t(i,j,k)>=Init_PRF_ALK(polygon_m(i,j),nmaxprofk(polygon_m(i,j),5)))then
           bio_t(i,j,k,iALKALINITY)    = InitVal_ALK(polygon_m(i,j),nmaxprofk(polygon_m(i,j),5))
         else

           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF_ALK(polygon_m(i,j),K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF_ALK(polygon_m(i,j),K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF_ALK(polygon_m(i,j),K_pr)
           deltap=Init_PRF_ALK(polygon_m(i,j),K_pr)-Init_PRF_ALK(polygon_m(i,j),K_pr-1)

           bio_t(i,j,k,iALKALINITY)  =(r1*InitVal_ALK(polygon_m(i,j),K_pr  ) &
                                   +r2*InitVal_ALK(polygon_m(i,j),K_pr-1))/deltap

          endif


          bio_t(i,j,k,iOxygen)       =bio_t(i,j,k,iOxygen)  &
                            *(rhp_t(i,j,k)+rho) /1000.

          bio_t(i,j,k,iDIC)       =bio_t(i,j,k,iDIC)  &
                            *(rhp_t(i,j,k)+rho) /1000.

          bio_t(i,j,k,iAlkalinity)=bio_t(i,j,k,iAlkalinity) &
                            *(rhp_t(i,j,k)+rho) /1000.




! We deduce the excess of negative charge using Total Alkalinity from
! InitPelagic
! and Nitrate, Ammonium and Phosphate from obcdriver
         bio_t(i,j,k,iAlkalinity)= bio_t(i,j,k,iAlkalinity)  & ! total alkalinity
                                  -bio_t(i,j,k,iAmmonium)    &
                                  +bio_t(i,j,k,iPhosphate)   &
                                  +bio_t(i,j,k,iNitrate)


         pCO2W(i,j)=410.


! 1D
             if    (-depth_t(i,j,k)<=ProfPH(1,1   ))then
               spH(i,j,k) = InitPH(1,1)

             ELSEif(-depth_t(i,j,k)>=ProfPH(1,22))then
               spH(i,j,k) = InitPH(1,22)

             ELSE

              K_pr=1
             do while(-depth_t(i,j,k)>ProfPH(1,K_pr))
               K_pr=K_pr+1
              enddo

              spH(i,j,k)=(                              &
              ( ProfPH(1,K_pr  )+depth_t(i,j,k))*InitPH(1,K_pr-1) &
             +(-ProfPH(1,K_pr-1)-depth_t(i,j,k))*InitPH(1,K_pr  ) &
             )/(ProfPH(1,K_pr  )-ProfPH(1,K_pr-1))
           endif

      endif ! mask 
      enddo
      enddo
      enddo


!=====================================================================*
! 3. Initialisation of the variables                                  *
! END                                                                 *
!=====================================================================*


!=====================================================================*
! 5. Initialisation of output variables                               *
! BEGINNING                                                           *
!=====================================================================*
! mise a zero des diagnostics 0D

       do vb=1,vbmax

        EXPORT_200(vb)=0.
        SUM_EXPORT_200(vb)=0.
        EXPORT_200(vbmax+vb)    =0.
        SUM_EXPORT_200(vbmax+vb)=0.
        EXPORT_200(2*vbmax+vb)    =0.
        SUM_EXPORT_200(2*vbmax+vb)=0.
        EXPORT_plat(vb)     =0.
        SUM_EXPORT_plat(vb) =0.
        EXPORTGDL_200(vb)=0.
        SUM_EXPORTGDL_200(vb)=0.
        EXPORTGDL_200(vbmax+vb)    =0.
        SUM_EXPORTGDL_200(vbmax+vb)=0.
        EXPORTGDL_200(2*vbmax+vb)    =0.
        SUM_EXPORTGDL_200(2*vbmax+vb)=0.
        EXPORTGDL_plat(vb)     =0.
        SUM_EXPORTGDL_plat(vb) =0.
        EXPORTDYF_200(vb)=0.
        SUM_EXPORTDYF_200(vb)=0.
        EXPORTDYF_200(vbmax+vb)    =0.
        SUM_EXPORTDYF_200(vbmax+vb)=0.
        EXPORTDYF_200(2*vbmax+vb)    =0.
        SUM_EXPORTDYF_200(2*vbmax+vb)=0.
        EXPORT_200_1000m(vb)=0.
        SUM_EXPORT_200_1000m(vb)=0.
        EXPORTDYF_1000(vb)=0.
        SUM_EXPORTDYF_1000(vb)=0.
        EXPORTDYF_1000(vbmax+vb)    =0.
        SUM_EXPORTDYF_1000(vbmax+vb)=0.
        EXPORTDYF_1000(2*vbmax+vb)    =0.
        SUM_EXPORTDYF_1000(2*vbmax+vb)=0.
        EXPORTDYF_BOT(vb)=0.
        SUM_EXPORTDYF_BOT(vb)=0.


! Apport rivieres
        BIL_FLEUVES(vb)=0.
          do kr=1,nriver
           BIL_FLEUVE(kr,vb)=0.
          enddo

! Export frontieres
        SUM_OB_SOUTH(vb)=0.
        SUM_OB_WEST(vb)=0.
        SUM_OB_EAST(vb)=0.
        FLUX_OB_SOUTH(vb)=0.
        FLUX_OB_WEST(vb)=0.
        FLUX_OB_EAST(vb)=0.

          do k=1,kmax
           SUM_OB_SOUTH_VERT(k,vb)=0.
           SUM_OB_WEST_VERT(k,vb)=0.
           SUM_OB_EAST_VERT(k,vb)=0.
           FLUX_OB_SOUTH_VERT(k,vb)=0.
           FLUX_OB_WEST_VERT(k,vb)=0.
           FLUX_OB_EAST_VERT(k,vb)=0.
!           do J=nbio1+1,nbio2
!            SUM_OB_WEST_PT(J-nbio1+1,k,vb)=0.
!            FLUX_OB_WEST_PT(J-nbio1+1,k,vb)=0.
!           enddo
!           do I=mbio1+1,mbio2
!            SUM_OB_SOUTH_PT(I-mbio1+1,k,vb)=0.
!            FLUX_OB_SOUTH_PT(I-mbio1+1,k,vb)=0.
!           enddo
          enddo

      enddo !vb
      
      TPS_PPB=0
      TPS_STRADA=0
      do k=1,4
      SUMPPBTOTAL(k)   = 0.
      PPBTOTAL(k)      = 0.
      SUMNPBTOTAL(k)   = 0.
      NPBTOTAL(k)      = 0.
      SUMRPBTOTAL(k)   = 0.
      RPBTOTAL(k)      = 0.
      SUMEXPTOTAL(k)   = 0.
      SUMEXPTOTAL(k+4) = 0.
      EXPTOTAL(k)      = 0.
      EXPTOTAL(k+4)    = 0.
      SUMPPBLIG(k)     = 0.
      SUMPPBDYF(k)     = 0.
      SUMPPBGDL(k)     = 0.
      PPBLIG(k)        = 0.
      PPBDYF(k)        = 0.
      PPBGDL(k)        = 0.
      SUMNETPPBTOTAL(k)= 0.
      NETPPBTOTAL(k)   = 0.
      SUMRESPPHYTOTOTAL(k)= 0.
      RESPPHYTOTOTAL(k)= 0.
      SUMNETPPBLIG(k)  = 0.
      SUMNETPPBDYF(k)  = 0.
      SUMNETPPBGDL(k)  = 0.
      NETPPBLIG(k)     = 0.
      NETPPBDYF(k)     = 0.
      NETPPBGDL(k)     = 0.
      SUMPPBTOTAL_plat(k) = 0.
      PPBTOTAL_plat(k) = 0.
      enddo

      DO K=1,3
      SUMPPBGDL_platcl(K)   = 0.
      PPBGDL_platcl(K)      = 0.
      ENDDO

      SUMPPBGDL_plat   = 0.
      SUMRPBGDL_plat   = 0.
      SUMNPBGDL_plat   = 0.
      SUMNETPPBGDL_plat= 0.
      PPBGDL_plat      = 0.
      RPBGDL_plat      = 0.
      NPBGDL_plat      = 0.
      NETPPBGDL_plat   = 0.
      SUMGBACGDL_plat  = 0.
      SUMUPTBACTDONGDL_plat    = 0.
      SUMEXCHETERONGDL_plat    = 0.
      SUMREMPOCGDL_plat = 0.
      SUMREMPONGDL_plat = 0.
      SUMLOSTPOCGDL_plat = 0.
      SUMLOSTPONGDL_plat = 0.
      SUMRESPPHYTOGDL_plat = 0.
      SUMRESPTOTGDL_plat = 0.
      SUMRESPZOOGDL_plat = 0.
      SUMRPBTOTAL(5)   = 0.
      RPBTOTAL(5)      = 0.
      SUMNITRifTOTAL   = 0.
      NITRifTOTAL      = 0.
      SUMRESPTOTTOTAL  = 0.
      RESPTOTTOTAL     = 0.
      SUMEXUCTOTTOTAL  = 0.
      EXUCTOTTOTAL     = 0.
      SUMGBACTOTAL     = 0. 
      SUMGBACDYF       = 0.
      SUMGBACGDL       = 0.  
      SUMGRAZCTOTAL    = 0.
      GRAZCTOTAL       = 0.

! Exportation flux
      SUMEXPDYF   = 0.
      SUMEXPGDL   = 0.

      SUM_EXPORTC_BOT = 0.
      SUM_EXPORTN_BOT = 0.
      SUM_EXPORTP_BOT = 0.
      SUM_EXPORTSI_BOT = 0.
      EXPORTC_BOT = 0.
      EXPORTN_BOT = 0.
      EXPORTP_BOT = 0.
      EXPORTSI_BOT = 0.
      SUM_SEDC    = 0.
      SUM_SEDP    = 0.
      SUM_SEDN    = 0.
      SUM_SEDSI   = 0.
          SEDC    = 0.
          SEDP    = 0.
          SEDN    = 0.
          SEDSI   = 0.

! mise a zero des diagnostics 2D

      TPS_PPB_2D=0
      TPS_STRADA_2D=0
      do j1=1,jmax ! debut boucle sur j
      do i1=1,imax ! debut boucle sur i
      do k=1,4 !8 attention ici...
        PPB2D(i1,j1,k)   =0.
      enddo
      do k=1,4
        NPB2D(i1,j1,k)   =0.
        RPB2D(i1,j1,k)   =0.
        CHL2D(i1,j1,k)   =0.
      enddo
      do k=1,27 ! 40   mais de qui se moque-t-on... 40...
        EXP2D(i1,j1,k)   =0.
      enddo
!      do k=1,5
!        NITRif2D(i1,j1,k)  =0.
!        RESP2D(i1,j1,k)    =0.
!      enddo
        EXUCTOT2D(i1,j1) =0.
        RPB2D(i1,j1,5)   =0.
!claude verifier qu'il faut bien cidessous sum_ devant export caroline n'a pas
        SUM_EXPORTC_BOT_2D(i1,j1)    =0.
        SUM_EXPORTN_BOT_2D(i1,j1)    =0.
        SUM_EXPORTP_BOT_2D(i1,j1)    =0.
        SUM_EXPORTSI_BOT_2D(i1,j1)    =0.
      enddo    ! fin de boucle i1
      enddo    ! fin de boucle j1


!       
!      do vb=1,vbmax
!      sumtot_bio(1,vb)=0.
!      sumcoltot_bio(1,vb)=0.
!      sumshelf_bio(1,vb)=0.
!      sumcolshelf_bio(1,vb)=0.
!      sumdeep_bio(1,vb)=0.
!      sumcoldeep_bio(1,vb)=0.
!      enddo  

!! Integration area
!!! Dyfamed
!      IDYF = 276
!      JDYF = 222
!
!
!!! Gulf of Lion ds le repere global
!      i1gdl = 113
!      i2gdl = 214
!      j1gdl = 172
!      j2gdl = 232
!
!! bilan GDL inutile 
!! il faudrait prendre la zone : I=113 -> 214
!!                               J=157 -> 223 
!! Medoc
!      I1MEDOC = 152
!      I2MEDOC = 215
!      J1MEDOC = 122
!      J2MEDOC = 183
!
!      ILIG = 220 

      sumareagdl   = 0.
      sumareagdl_plat = 0.
      sumareatotal = 0.
      sumareatotal_plat = 0.
      sumarealig   = 0.
      sumarea1000m = 0.
      sumareamedoc = 0.
      sumareadyf = 0.
      sumareadiag = 0. 
      sumareawmed = 0.
      sumareaemed = 0.


      do i1=1,imax ! debut boucle sur j
      do j1=1,jmax ! debut boucle sur i
      i2=i1+par%timax(1)      ! ds le repere global
      j2=j1+par%tjmax(1)

!      if (j2>=j1domain.and.j2<=j2domain.and.i2>=i1domain.and.i2<=i2domain) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(mask_t(i1,j1,kmax+1)==1) then !>>>>>>>>>>>>>>
       x2=dxdy_t(i1,j1)*mask_i_w(i1)*mask_j_w(j1)    ! surface de la maille

        if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then
!       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then
         sumareagdl = sumareagdl + x2
          if(-depth_w(i1,j1,kmin_w(i1,j1))<=200.) then
          sumareagdl_plat = sumareagdl_plat + x2
          endif
       endif

       if (j2>=jdyf1.and.j2<=jdyf2.and.i2>=idyf1.and.i2<=idyf2) then
         sumareadyf = sumareadyf + x2
       endif

       if (j2>=j1medoc.and.j2<=j2medoc.and.i2>=i1medoc.and.i2<=i2medoc) then
         if(depth_t(i1,j1,1)<-2000.)                    &
         sumareamedoc = sumareamedoc + x2
       endif

       sumareatotal = sumareatotal + x2

! bassin ouest
            if(lon_t(i1,j1)<10.*deg2rad.and.lon_t(i1,j1)>-5.6*deg2rad) then
!             if(lon_t(i1,j1)<10.*deg2rad) then
             sumareawmed = sumareawmed + x2
         elseif(lon_t(i1,j1)>10.*deg2rad.and.lon_t(i1,j1)<15.*deg2rad.and. &
                lat_t(i1,j1)>37.*deg2rad.and.lat_t(i1,j1)<42.*deg2rad)then 
             sumareawmed = sumareawmed + x2
         elseif(lon_t(i1,j1)>10.*deg2rad.and.lon_t(i1,j1)<12.25*deg2rad.and. &
                lat_t(i1,j1)>42.*deg2rad.and.lat_t(i1,j1)<44.25*deg2rad)then
             sumareawmed = sumareawmed + x2
         elseif(lon_t(i1,j1)>15.*deg2rad.and.lon_t(i1,j1)<16.25*deg2rad.and. &
                lat_t(i1,j1)>38.*deg2rad.and.lat_t(i1,j1)<40.25*deg2rad)then
             sumareawmed = sumareawmed + x2
!bassin est
         else
             if(lon_t(i1,j1)>5.6*deg2rad) then
             sumareaemed = sumareaemed + x2
             endif
!          if(lon_t(i1,j1)<-5.6)'passe ici 5.6')
         endif


       if(-depth_w(i1,j1,kmin_w(i1,j1))<=200.) then
       sumareatotal_plat = sumareatotal_plat + x2
       endif

!      if (i2>=ILIG) sumarealig = sumarealig + x2      ! repere global
       if (j2>=j1ligure.and.j2<=j2ligure.and. &
           i2>=i1ligure.and.i2<=i2ligure) then
       if(depth_t(i1,j1,1)<-1000.)                    &
       sumarea1000m = sumarea1000m + x2
       if(depth_t(i1,j1,1)<-2000.)                    &
       sumarealig = sumarealig + x2
       endif

      if (j2>=j1diag.and.j2<=j2diag.and. &
          i2>=i1diag.and.i2<=i2diag) then
       if(-depth_w(i1,j1,kmin_w(i1,j1))>=1500.)                    &
       sumareadiag = sumareadiag + x2
       endif


      endif ! test masque
!      endif ! test domain
      enddo    ! fin de boucle i1
      enddo    ! fin de boucle j1


      SUM_BIOPOC_P1=0.
      SUM_BIODOC_P1=0.
      SUM_BIOPON_P1=0.
      SUM_BIODON_P1=0.
      SUM_BIONO3_P1=0.
      SUM_BIONH3_P1=0.

      SUM_HBIO_P1=0.
      SUM_HBION_P1=0.

! claude: a partir de ci-dessous les sumareatotal representent la somme sur tous les procs
      call mpi_allreduce(sumareatotal,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareatotal=sum1glb
      call mpi_allreduce(sumareatotal_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareatotal_plat=sum1glb
      call mpi_allreduce(sumarea1000m,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumarea1000m=sum1glb
      call mpi_allreduce(sumareagdl,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareagdl=sum1glb
      call mpi_allreduce(sumareagdl_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareagdl_plat=sum1glb
      call mpi_allreduce(sumarealig,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumarealig=sum1glb
      call mpi_allreduce(sumareamedoc,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareamedoc=sum1glb
      call mpi_allreduce(sumareadyf,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareadyf=sum1glb
      call mpi_allreduce(sumareadiag,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareadiag=sum1glb
      call mpi_allreduce(sumareawmed,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareawmed=sum1glb
      call mpi_allreduce(sumareaemed,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sumareaemed=sum1glb


      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=4,file='messages',position='append')
      write(4,*)'------------------------'
      write(4,*)'Aire des domaines de cumul de diag BIO (sumarea en m2)'
      write(4,*)'sumareatotal',sumareatotal
      write(4,*)'sumareagdl'  ,sumareagdl
      write(4,*)'sumareatotal_plat',sumareatotal_plat
      write(4,*)'sumareagdl_plat',sumareagdl_plat
      write(4,*)'sumarealig',sumarealig
      write(4,*)'sumarea1000m',sumarea1000m
      write(4,*)'sumareamedoc',sumareamedoc
      write(4,*)'sumareadyf',sumareadyf
      write(4,*)'sumareaemed',sumareaemed
      write(4,*)'sumareawmed',sumareawmed
      close(4)
      endif !#mpi-->>-->


!=====================================================================*
! 5. Initialisation of output variables                               *
! END                                                                 *
!=====================================================================*


       if(par%rank==0) write(6,*)'Passe par InitPelagic',sumareaemed,sumareawmed,deg2rad,rad2deg
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10

      END SUBROUTINE InitNutDOxy
