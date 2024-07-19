      SUBROUTINE InitPelagic

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
      ! modif polygone et initialisation Alex
      integer,parameter :: nmaxprof=43, npoly=9
      integer :: npoints_poly, n, polygon_m(mbio1:mbio2,nbio1:nbio2)

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
                      ! New stuff 
                      ! valeur de 1 à 9 pour les sous-régions: ALB, ALG, GOL,
                      ! TYR, CMED, ION, ADR, AEG, LEV
                      ,InitVal_NIT(npoly,nmaxprof)   &
                      ,InitVal_PHO(npoly,nmaxprof)   &
                      ,InitVal_SIL(npoly,nmaxprof)   &
                      ,Init_PRF(nmaxprof)        &
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
           InitPH3D2(iglb+2,jglb+2,kmax),InitMODC(iglb+2,jglb+2,kmax)            
 
       integer l,m,nmax10y(14)
       real Prof10y(14,2511),var10y(14,2511,10)
       character*60 filename


       InitOC=1

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
! 1. Initialisation of the general parameters                         *
! BEGINNING                                                           *
!=====================================================================*


      if(NumPelagicBoxes.NE.kmax) then
      write(6,*)'Parameter Error 1:'
      write(6,*)'NumpelagicBoxes should be equal to kmax'
      write(6,*)'The value of ' 
      write(6,*)'NumPelagicBoxes in xModuledeclaration.f90 is '      &
                ,NumPelagicBoxes
      write(6,*)'And of vbmax+1 in parameter is ',vbmax+1
      write(6,*)'Correct and compile again'
      STOP 'Parameter Error1'
      endif 


      if(xsvardeclare.NE.vbmax) then
      write(6,*)'Parameter Error 2:'
      write(6,*)'xsvardeclare should be equal to vbmax'
      write(6,*)'The value of'
      write(6,*)'xsvardeclare in xModuledeclaration.f90 is '         &
                ,xsvardeclare
      write(6,*)'And of vbmax (Nombre de variables) in notebook_bio is'&
                       ,vbmax
      write(6,*)'Correct and compile again'
      STOP 'Parameter Error2'
      endif

!=====================================================================*
! 1. Initialisation of the general parameters                         *
! END                                                                 *
!=====================================================================*


!=====================================================================*
! 2. Reading of the light and pelagic parameters                      *
! BEGINNING                                                           *
!=====================================================================*

      CALL LightParameters

      CALL Eco3mParameters

!=====================================================================*
! 2. Reading of the light and pelagic parameters                      *
! END                                                                 *
!=====================================================================*


!=====================================================================*
! 3. Initialisation of the variables                                  *
! BEGINNING                                                           *
!=====================================================================*

! State variables
!-----------------

      do  I = 1 , NumpelagicBoxes

      ZooNanoC(I)  = 0.
      ZooMicroC(I) = 0.
      ZooMesoC(I)  = 0.
      SyneC(I)     = 0.
      SyneN(I)     = 0.
      SyneChl(I)   = 0.
      SyneP(I)     = 0.
      NanoC(I)     = 0.
      NanoN(I)     = 0.
      NanoChl(I)   = 0.
      NanoP(I)     = 0.
      DiaC(I)      = 0.
      DiaN(I)      = 0.
      DiaChl(I)    = 0.
      DiaSi(I)     = 0.
      DiaP(I)      = 0.
      BactC(I)     = 0.
      S_MOPC(I)    = 0.
      S_MOPN(I)    = 0.
      S_MOPP(I)    = 0.
      S_MOPChl(I)  = 0.
      S_MOPSi(I)   = 0.
      L_MOPC(I)    = 0.
      L_MOPN(I)    = 0.
      L_MOPSi(I)   = 0.
      L_MOPP(I )   = 0.
      MODC(I)      = 0.
      MODN(I)      = 0.
      MODP(I)      = 0.
      Nitrate(I)   = 0.
      Ammonium(I)  = 0.
      Silice(I)    = 0.
      Phosphate(I) = 0.
      Oxygen(I)    = 0.
!      ODU(I)       = 0.
!      MIP(I)       = 0.
      enddo


! Reading of the notebook

      if(par%rank==0) write(6,'(A,A)')'Ready to read  &
                       notebook_initpelagic',nomfichier(28)

      open(unit=3,file=nomfichier(28)) ! open the notebook

      do k=1,10
      read(3,*)
      enddo 

      read(3,*)ChoiceInitPelagic

      if(ChoiceInitPelagic>=1) then    ! 17/02/09
      read(3,*)
      read(3,*)
      read(3,'(A)')DirInitFile 
      read(3,*)DATE
      read(3,*)CoefDiaC
      read(3,*)CoefZooC
      read(3,*)CoefSyneC
      read(3,*)CoefNanoC
      read(3,*)PropMOP
      read(3,*)InitLPOSi
      read(3,*)InitAmmonium
      ELSE
       do k=1,10
        read(3,*)
       enddo
      endif

      read(3,*)
      read(3,*)
      read(3,*)InitFemme

      if(InitFemme==1) then
      read(3,*)ifemme,JFemme
      endif  

    
      read(3,*)                                                         ! 18/07/08
      read(3,*)                                                         ! 18/07/08
      read(3,*)Inudging                                                 ! 18/07/08   
      if(Inudging == 1)then
      write(*,*)'Nudging commenté dans biology.F90'
      write(*,*)'tres gourmand en memoire:tableau bioinit'
      write(*,*)'pas encore mis en service'
      stop 'Inudging dans InitPelagic'
      endif
      read(3,*)I1,I2,I3                                                 ! 18/07/08
      I4=0.                                                             ! 18/07/08
      I5=0.                                                             ! 18/07/08
      I6=0.                                                             ! 18/07/08
      CALL DATETOKOUNT(I1,I2,I3,I4,I5,I6)                               ! 18/07/08
      knudgBegin=XDTK_OUT                                               ! 18/07/08
      read(3,*)I1,I2,I3                                                 ! 18/07/08
      CALL DATETOKOUNT(I1,I2,I3,I4,I5,I6)                               ! 18/07/08
      knudgEnd=XDTK_OUT                                                 ! 18/07/08
      Nudgperiod=(knudgEnd-knudgBegin)*DTI_FW/86400.                    ! 18/07/08

      close(3) 

      if(2/(CoefDiaC*(minNCDia+maxNCDia))>ChlNDiaMax) then
      write(6,*)'Error initialisation:'
      write(6,*)'ChlNDia > ChlNDiaMax'
      write(6,*)'ChlNDia=',2./(CoefDiaC*(minNCDia+maxNCDia))
      write(6,*)'ChlNDiaMax=',ChlNDiaMax
      write(6,*)'Change CoefDiaC in InitPelagic.F'
      STOP'donc dans InitPelagic'
      endif

      if(2/(CoefNanoC*(minNCNano+maxNCNano))>ChlNNanoMax) then
      write(6,*)'Error initialisation:'
      write(6,*)'ChlNNano > ChlNNanoMax'
      write(6,*)'ChlNNano=',2./(CoefNanoC*(minNCNano+maxNCNano))
      write(6,*)'ChlNNanoMax=',ChlNNanoMax
      write(6,*)'Change CoefNanoC in InitPelagic.F'
      STOP'donc dans InitPelagic'
      endif

      if(2/(CoefSyneC*(minNCSyne+maxNCSyne))>ChlNSyneMax) then
      write(6,*)'Error initialisation:'
      write(6,*)'ChlNSyne > ChlNSyneMax'
      write(6,*)'Change CoefSyneC in InitPelagic.F'
      STOP'donc dans InitPelagic'
      endif

! Reading of notebook_biobcforcing


      if(ChoiceInitPelagic==3) then

      write(6,'(A,A)')'Ready to read notebook_biobcforcing'   &
                        ,nomfichier(29)

      open(unit=3,file=nomfichier(29)) ! open the notebook      

      read(3,*) 
      read(3,*)
      read(3,*)
      read(3,*)

! section with forcing by 3D fields:
      read(3,*)IBIOBC_AF
      read(3,*)INTERP_LR

      if(INTERP_LR==1.and.IBIOBC_AF==1) then
      read(3,'(A)')biobcfile        
      write(6,*)'BIOBCfile',biobcfile

      read(3,*)    BIOBCINFO
      read(3,*) 
      read(3,*) 
      read(3,*) 
      read(3,*)DATEBIOBC(1)           &
              ,DATEBIOBC(2)           &
              ,DATEBIOBC(3)           &
              ,DATEBIOBC(4)           &
              ,DATEBIOBC(5)           &
              ,DATEBIOBC(6)


      endif ! INTERP_LR & IBIOBC_AF
 
      close(3)

      endif ! ChoiceInitPelagic


! Constant parameters
      MassemolC  = 12. !g/mol
      MassemolN  = 14. !g/mol

      RedfieldNC  = 1/6.625
      RedfieldPC  = 1/106.


! directory of the initial profil files
       do 20 K=1,90
         if(DirInitFile(K:K).eq.' ') then !-------->
          LDir=K-1
          GOTO 21
         endif                         !-------->
   20  CONTINUE
   21  CONTINUE
       if(DirInitFile(LDir:LDir).eq.'/')LDir=LDir-1
       LDir=LDir+1
       DirInitFile(LDir:LDir)='/'


! a faire:
! check date date doit correspondre a moins de 16 j 

      if(ChoiceInitPelagic>=1) then

      write(Fic(1),'(A,A6,A8 )')                                           &
                            DirInitFile(1:LDIR),DATE,'_DiaInit'
      write(Fic(2),'(A,A6,A9 )')                                           &
                            DirInitFile(1:LDIR),DATE,'_NanoInit'            
      write(Fic(3),'(A,A6,A9 )')                                           &
                            DirInitFile(1:LDIR),DATE,'_SyneInit'            
      write(Fic(4),'(A,A6,A17)')                                           &
                            DirInitFile(1:LDIR),DATE,'_NitrateInitaugpa'    
      write(Fic(5),'(A,A6,A19)')                                           &
                          DirInitFile(1:LDIR),DATE,'_PhosphateInitaugpa'    
      write(Fic(6),'(A,A6,A16)')                                           &
                            DirInitFile(1:LDIR),DATE,'_SiliceInitaugpa'     
!      write(Fic(4),'(A,A6,A12)')
!     &                      DirInitFile(1:LDIR),DATE,'_NitrateInit'
!      write(Fic(5),'(A,A6,A14)')
!     &                      DirInitFile(1:LDIR),DATE,'_PhosphateInit'
!      write(Fic(6),'(A,A6,A11)')
!     &                      DirInitFile(1:LDIR),DATE,'_SiliceInit'
      write(Fic(7),'(A,A6,A9 )')                                       &
                            DirInitFile(1:LDIR),DATE,'_MODCInit'
      write(Fic(8),'(A,A6,A10 )')                                      &
                            DirInitFile(1:LDIR),DATE,'_BactCInit'
      write(Fic(9),'(A,A6,A12 )')                                      &
                            DirInitFile(1:LDIR),DATE,'_ZooNanoInit'
      write(Fic(10),'(A,A6,A13 )')                                     &
                            DirInitFile(1:LDIR),DATE,'_ZooMicroInit'
      write(Fic(11),'(A,A6,A12 )')                                     &
                            DirInitFile(1:LDIR),DATE,'_ZooMesoInit'
      WRITE(Fic(12),'(A,A6,A8 )')                                      &
                            DirInitFile(1:LDIR),DATE,'_POCInit'
      WRITE(Fic(13),'(A,A6,A8 )')                                      &
                            DirInitFile(1:LDIR),DATE,'_PONInit'
!      write(Fic(12),'(A,A6,A9 )')                                      &
!                            DirInitFile(1:LDIR),DATE,'_SPOCInit'
!      write(Fic(13),'(A,A6,A9 )')                                      &
!                            DirInitFile(1:LDIR),DATE,'_SPONInit'
!      write(Fic(14),'(A,A6,A9)')                                       &
!                            DirInitFile(1:LDIR),DATE,'_SPOPInit'
!      write(Fic(15),'(A,A6,A11)')                                      &
!                            DirInitFile(1:LDIR),DATE,'_SPOChlInit'
!      write(Fic(16),'(A,A6,A10)')                                      &
!                            DirInitFile(1:LDIR),DATE,'_SPOSiInit'
!      write(Fic(17),'(A,A6,A9 )')                                      &
!                            DirInitFile(1:LDIR),DATE,'_LPOCInit'
!      write(Fic(18),'(A,A6,A9 )')                                      &
!                            DirInitFile(1:LDIR),DATE,'_LPONInit'
!      write(Fic(19),'(A,A6,A9 )')                                      &
!                            DirInitFile(1:LDIR),DATE,'_LPOPInit'
!      write(Fic(20),'(A,A6,A10 )')                                     &
!                            DirInitFile(1:LDIR),DATE,'_LPOSiInit'
      write(FicLPC_Nit,'(A,A17)') &
                           DirInitFile(1:LDIR),'Moogl3InitNitrate'
      write(FicLPC_Pho,'(A,A19)') &
                         DirInitFile(1:LDIR),'Moogl3InitPhosphate'
!     write(FicLPC_Nit,'(A,A17)')                                    &
!                            DirInitFile(1:LDIR),'MoogliInitNitrate'
!      write(FicLPC_Pho,'(A,A19)')                                    &
!                          DirInitFile(1:LDIR),'MoogliInitPhosphate'

! claude : Caroline a moins de fichiers et moins de NumVar

      NumVar(1) = iDiaChl
      NumVar(2) = iNanoChl
      NumVar(3) = iSyneChl
      NumVar(4) = iNitrate
      NumVar(5) = iPhosphate
      NumVar(6) = iSilice
      NumVar(7) = iMODC
      NumVar(8) = iBactC
      NumVar(9) = iZooNanoC
      NumVar(10) = iZooMicroC
      NumVar(11) = iZooMesoC
      NumVar(12) = iSMOPC
      NumVar(13) = iSMOPN
!      NumVar(14) = iSMOPP
!      NumVar(15) = iSMOPChl
!      NumVar(16) = iSMOPSi
!      NumVar(17) = iLMOPC
!      NumVar(18) = iLMOPN
!      NumVar(19) = iLMOPP
!      NumVar(20) = iLMOPSi



!! Initialisation of diatomees Chl, syneChl, POC, PON, MOD and nutrients
      do k2=1,13 !20

!     write(*,*)'Fic',K2,Fic(K2)
      open(unit=3,file=Fic(K2))
!     write(*,*)'Fic',K2,Fic(K2)
      read(3,*)x1     
      nmax=int(x1)
      do k=1,nmax 
      read(3,*)Prof(K2,K),InitValue(K2,K)
      enddo
      close(3)    

      if(ChoiceInitPelagic==2) then

       if(NumVar(K2)==iNitrate)then
      open(unit=3,file=FicLPC_Nit)
      read(3,*)nmaxlpc
      do k=1,nmaxlpc
      read(3,*)ProfLPC(K2,k),InitValueLPC(K2,k)
      enddo
      close(3)
      endif

       if(NumVar(K2)==iPhosphate)then
      open(unit=3,file=FicLPC_Pho)
      read(3,*)nmaxlpc
      do k=1,nmaxlpc
      read(3,*)ProfLPC(K2,k),InitValueLPC(K2,k)
      enddo
      close(3)
       endif


      endif

         do  k=1,kmax
         do  i=0,imax+1 ! surdimensionne pour englober obc
         do  j=0,jmax+1

           if(mask_t(i,j,k)==1)then

             if    (-depth_t(i,j,k)<=Prof(K2,1   ))then
               bio_t(i,j,k,NumVar(K2)) = InitValue(K2,1)

             ELSEif(-depth_t(i,j,k)>=Prof(K2,nmax))then
               bio_t(i,j,k,NumVar(K2)) = InitValue(K2,nmax)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof(K2,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,NumVar(K2))=(                              &
              ( Prof(K2,K_pr  )+depth_t(i,j,k))*InitValue(K2,K_pr-1) &
             +(-Prof(K2,K_pr-1)-depth_t(i,j,k))*InitValue(K2,K_pr  ) &
             )/(Prof(K2,K_pr  )-Prof(K2,K_pr-1))

            endif !prof3D

            if(ChoiceInitPelagic==2) then

            if(NumVar(K2)==iNitrate.OR.NumVar(K2)==iPhosphate) then
             if    (-depth_t(i,j,k)<=ProfLPC(K2,1   ))then
               anyvar3d(i,j,k) = InitValueLPC(K2,1)

             ELSEif(-depth_t(i,j,k)>=ProfLPC(K2,nmaxlpc))then
               anyvar3d(i,j,k) = InitValueLPC(K2,nmaxlpc)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>ProfLPC(K2,K_pr))
               K_pr=K_pr+1
              enddo

              anyvar3d(i,j,k)=(                         &
              ( ProfLPC(K2,K_pr  )+depth_t(i,j,k))      &
                    *InitValueLPC(K2,K_pr-1)            &
             +(-ProfLPC(K2,K_pr-1)-depth_t(i,j,k))      &
                    *InitValueLPC(K2,K_pr  )            &
             )/(ProfLPC(K2,K_pr  )-ProfLPC(K2,K_pr-1))


           endif !prof3D

           endif ! ChoiceInitPelagic 2 

            endif !numvar



           if(ChoiceInitPelagic==2) then
! Idem pour frontiere Est
      YP1=1.
      YP2=1.


      if(i>mbio2.or.(h_w(i,j)>0..and.h_w(i,j)<1800.))then
               Y1=MIN(1.D0,MAX(0.D0,(FLOAT(J)-YP1)/YP2))
 
              if(NumVar(K2)==iPhosphate.OR.NumVar(K2)==iNitrate)then
                 bio_t(i,j,k,NumVar(K2))=                            &
                 (1.-Y1)*bio_t(i,j,k,NumVar(K2))                     &
                 +Y1*anyvar3d(i,j,k)          
               endif !NumVar
              if(NumVar(K2)==iSilice)then
                 bio_t(i,j,k,NumVar(K2))=bio_t(i,j,k,iNitrate)/1.4
              endif


      endif !I

      endif ! ChoiceInitPelagic 2

           endif !  mask_t
         enddo
         enddo
         enddo


         enddo !K2




      do  k=1,kmax
      do  i=0,imax+1 ! surdimensionne pour englober obc
      do  j=0,jmax+1

           if(mask_t(i,j,k)==1)then


!! Deduction of other Dia and Zoo from DiaChl
         bio_t(i,j,k,iDiaC)  =                 CoefDiaC        &
       * bio_t(i,j,k,iDiaChl)
         bio_t(i,j,k,iDiaN)  = (  minNCDia +  maxNCDia ) / 2.  &
       * bio_t(i,j,k,iDiaC)
         bio_t(i,j,k,iDiaP)  = (  minPCDia +  maxPCDia ) / 2.  &
       * bio_t(i,j,k,iDiaC)
         bio_t(i,j,k,iDiaSi) = ( minSiCDia + maxSiCDia ) / 2.  &
       * bio_t(i,j,k,iDiaC)


!      if(I==144.and.J==31.) write(6,*)bio_t(i,j,k,iDiaC),'DiaC'
!     &     ,bio_t(i,j,k,iDiaChl),CoefDiaC


!! Deduction of other Nano from NanoChl
         bio_t(i,j,k,iNanoC) =                  CoefNanoC       &
       * bio_t(i,j,k,iNanoChl)
         bio_t(i,j,k,iNanoN) = (  minNCNano +  maxNCNano ) / 2. &
       * bio_t(i,j,k,iNanoC)
         bio_t(i,j,k,iNanoP) = (  minPCNano +  maxPCNano ) / 2. &
       * bio_t(i,j,k,iNanoC)


!! Deduction of other Syne from SyneChl
         bio_t(i,j,k,iSyneC) =                  CoefSyneC       &
       * bio_t(i,j,k,iSyneChl)
         bio_t(i,j,k,iSyneN) = (  minNCSyne +  maxNCSyne ) / 2. &
       * bio_t(i,j,k,iSyneC)
         bio_t(i,j,k,iSyneP) = (  minPCSyne +  maxPCSyne ) / 2. &
       * bio_t(i,j,k,iSyneC)


! claude : caroline a commenté ci-dessous la reduction des MOD
!! Comme j'initialise avec une moyenne de decembre, je reduis les MOD
       bio_t(i,j,k,iMODC) = bio_t(i,j,k,iMODC) / 4.

!! Deduction of MODP from MODC 
      bio_t(i,j,k,iMODN)   =  RedfieldNC * bio_t(i,j,k,iMODC)
      bio_t(i,j,k,iMODP)   =  RedfieldPC * bio_t(i,j,k,iMODC)


!! Initialisation of Ammonium
      bio_t(i,j,k,iAmmonium)   = InitAmmonium  
      
! claude : Caroline n'initalise plus les bacteries
!! Bacteria are initialized according to MODC
      bio_t(i,j,k,iBactC)   = bio_t(i,j,k,iMODC)   / 200. * 3.
 
!claude : caroline initialise les SMOP et les LMOP
      bio_t(i,j,k,iSMOPC) = max(ZERO,                     &
                              bio_t(i,j,k,iSMOPC)         &
                          -(  bio_t(i,j,k,iDiaC)          & 
                            + bio_t(i,j,k,iNanoC)         & 
                            + bio_t(i,j,k,iSyneC)         &
                            + bio_t(i,j,k,iZooMicroC)     &
                            + bio_t(i,j,k,iZooNanoC)      &
                            + bio_t(i,j,k,iZooMesoC)      &
                            + bio_t(i,j,k,iBactC)    ))



!      IF(I.EQ.276.AND.J.EQ.222) THEN
!      print*,K,PROF3D_Z(I,J,K),bio_t(i,j,k,iSMOPC),bio_t(i,j,k,iDiaC)
!     &                       + bio_t(i,j,k,iNanoC)
!     &                       + bio_t(i,j,k,iSyneC)
!     &                       + bio_t(i,j,k,iZooMicroC)
!     &                       + bio_t(i,j,k,iZooNanoC)
!     &                       + bio_t(i,j,k,iZooMesoC)
!     &                       + bio_t(i,j,k,iBactC)
!      ENDIF

      bio_t(i,j,k,iSMOPN) = max(ZERO,                              &
                              bio_t(i,j,k,iSMOPN)                  &
                          -(  bio_t(i,j,k,iDiaN)                   &
                            + bio_t(i,j,k,iNanoN)                  &
                            + bio_t(i,j,k,iSyneN)                  &
                            + bio_t(i,j,k,iZooMicroC)*NCZooMicro   &
                            + bio_t(i,j,k,iZooNanoC)*NCZooNano     &
                            + bio_t(i,j,k,iZooMesoC)*NCZooMeso     &
                            + bio_t(i,j,k,iBactC)*NCBact    ))  

!!! Repartition in small and large particules
      bio_t(i,j,k,iLMOPC) = (1 - PropMOP) * bio_t(i,j,k,iSMOPC)
      bio_t(i,j,k,iSMOPC) =      PropMOP  * bio_t(i,j,k,iSMOPC)
      bio_t(i,j,k,iLMOPN) = (1 - PropMOP) * bio_t(i,j,k,iSMOPN)
      bio_t(i,j,k,iSMOPN) =      PropMOP  * bio_t(i,j,k,iSMOPN)

!!! Deduction of other small and large particulate organic matter
      bio_t(i,j,k,iLMOPP) =    RedfieldPC * bio_t(i,j,k,iLMOPC)
      bio_t(i,j,k,iSMOPP) =    RedfieldPC * bio_t(i,j,k,iSMOPC)

!!! Initialisation of small and large particulate organic silice
      bio_t(i,j,k,iLMOPSi) =                  InitLPOSi
      bio_t(i,j,k,iSMOPSi) =   DiaPropMOPSi * InitLPOSi   ! DiaPropMOPSi Ratio small/large organic m


        endif !  mask

      enddo
      enddo 
      enddo 

! Interpolation of nutrient fields from a general bio model
         if(ChoiceInitPelagic==3) then
          write(*,*) 'a reprendre de SYMPHONIE2008 --> stop'
         stop
!        CALL OBC_AF_BIO(1)      !  Spatial interpolation
!        CALL DTVAR_OBC_BIO(0)   ! Temporal interpolation
         endif

! density-dependent initialisation ofnutrients
      if(ChoiceInitPelagic==4) then   

      nmax1=68
      nmax2=68
      nmax3=68
      nmax4=21
      nmax5=22
!     nmax6=183
      nmax6=2311

      open(unit=2,file='/tmpdir/culses/bassin2/MERCATOFF/DyfamedmeanpH.txt')
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


!load oxygene EST et OUEST

!      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_oxyg_est.txt')
!      read(3,*)nmax3
     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_oxyg_EST2.txt')
      do k=1,nmax4
      read(3,*)InitValue_EST(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

!     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/OXYGEN_WEST.txt')
      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitOxyDyf201008.txt')
      do k=1,nmax6
      read(3,*)InitValue_WEST(1,k),Prof_oxyg_WEST(1,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_DIC_EST.txt')
      do k=1,nmax4
      read(3,*)InitValue_EST_DIC(1,k),Prof_oxyg(1,k)
      enddo
      close(3)
! dic
      open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_DIC_WEST.txt')
      do k=1,nmax4
      read(3,*)InitDIC(1,k),ProfDIC(1,k)
      enddo
      close(3)

     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_DIC_EST_WEST.txt')
      do k=1,nmax4
      read(3,*)InitDIC(3,k),ProfDIC(3,k)
      enddo
      close(3)

     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_DIC_EST_EST.txt')
      do k=1,nmax4
      read(3,*)InitDIC(4,k),ProfDIC(4,k)
      enddo
      close(3)


     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_DIC_GLT.txt')
      do k=1,nmax4
      read(3,*)InitDIC(5,k),ProfDIC(5,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/Dewex/Data/DyfamedmeanDIC.txt')
      do k=1,nmax5
      read(3,*)ProfDIC(2,k),InitDIC(2,k)
      enddo
      close(3)


! alk
     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_ALK_WEST.txt')
      do k=1,nmax4
      read(3,*)InitTAlk(1,k),ProfTAlk(1,k)
      enddo
      close(3)

     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_ALK_EST_WEST.txt')
      do k=1,nmax4
      read(3,*)InitTAlk(3,k),ProfTAlk(3,k)
      enddo
      close(3)

     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_ALK_EST_EST.txt')
      do k=1,nmax4
      read(3,*)InitTAlk(4,k),ProfTAlk(4,k)
      enddo
      close(3)

     open(unit=3,file='/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/interp_ALK_GLT.txt')
      do k=1,nmax4
      read(3,*)InitTAlk(5,k),ProfTAlk(5,k)
      enddo
      close(3)

      open(unit=3,file='/tmpdir/culses/Dewex/Data/DyfamedmeanTAlk.txt')
      do k=1,nmax5
      read(3,*)ProfTAlk(2,k),InitTAlk(2,k)
      enddo
      close(3)

! load Oxygene Bassin Feyçal 

! cyp

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_cyp.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_cyp(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! ion

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_ion.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_ion(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! lib

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_lib.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_lib(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! tun

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_tun.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_tun(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! sic

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_sic.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_sic(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! nad

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_nad.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_nad(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! sad

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_sad.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_sad(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! mlt

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_mlt.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_mlt(1,k),Prof_oxyg(1,k)
      enddo
      close(3)


! ege

      open(unit=3,file='/tmpdir/culses/bassin/NEMOFF1/PROFILNUTS/interp_oxyg_ege.txt')
      do k=1,nmax4
      read(3,*)InitValue_oxyg_ege(1,k),Prof_oxyg(1,k)
      enddo
      close(3)

! Lecture des profils de PROFILNUTS : Alex
! ALB
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ALB_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(1,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ALB_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(1,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ALB_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(1,k)
      enddo
      close(3)
! ALG
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ALG_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(2,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ALG_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(2,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ALG_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(2,k)
      enddo
      close(3)
! GOL
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/GOL_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(3,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/GOL_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(3,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/GOL_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(3,k)
      enddo
      close(3)
! TYR
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/TYR_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(4,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/TYR_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(4,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/TYR_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(4,k)
      enddo
      close(3)
! CMED
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/CMED_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(5,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/CMED_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(5,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/CMED_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(5,k)
      enddo
      close(3)
! ION
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ION_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(6,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ION_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(6,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ION_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(6,k)
      enddo
      close(3)
! ADR
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ADR_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(7,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ADR_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(7,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/ADR_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(7,k)
      enddo
      close(3)
! AEG
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/AEG_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(8,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/AEG_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(8,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/AEG_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(8,k)
      enddo
      close(3)
! LEV
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/LEV_nit.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_NIT(9,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/LEV_phos.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_PHO(9,k)
      enddo
      close(3)
      open(unit=3,file='/tmpdir/alex/BASSIN_3D/PROFILNUTS/LEV_sil.txt')
      do k=1,nmaxprof
      read(3,*)Init_PRF(k),InitVal_SIL(9,k)
      enddo
      close(3)

! Initialisation et lecture des polygones
! Attention à l'ordre des polygones dans le fichier
! ici l'ordre est :
! ALB, ALG, GOL, TYR, CMED, ION, ADR, AEG et LEV

      polygon_m(:,:)=0
      open(unit=11,file='/tmpdir/alex/BASSIN_3D/POLYGONE/liste_poly')
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



      do l=1,14
!     nmax10y=2511

      if(l<10) then
!      write (filename,"(A50,I1,A4)")'/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitOC2_',l,'.txt'
!      else
!      write(filename,"(A50,I2,A4)")'/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitOC2_',l,'.txt'
      write(filename,"(A52,I1,A4)")'/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitOCM8_',l,'.txt'
      else
      write(filename,"(A52,I2,A4)")'/tmpdir/culses/bassin2/NEMOFF1/PROFILNUTS/InitOCM8_',l,'.txt'
      endif

      open(unit=3,file=filename)
      read(3,*)nmax10y(l)
!      print*,'l',l,nmax10y(l)
      do k=1,nmax10y(l)
      read(3,*)Prof10y(l,k),(var10y(l,k,m),m=1,10)
! 13/10 ajout Caroline correction valeurs négatives      
      do m=1,10
      if(var10y(l,k,m)<0) var10y(l,k,m)=0.
      enddo
      enddo
      close(3)

      enddo



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

       ! lecture et interpolation des profils pour nitrate, silice phosphate PUIS
       ! Initialisation de l'oxygene, zoo, ALK, etc

       ! interpolation des profils de Nitrate, Phosphate et Silice

         if (depth_t(i,j,k)<=Init_PRF(1)) then
           bio_t(i,j,k,iNitrate)   = InitVal_NIT(polygon_m(i,j),1)
           bio_t(i,j,k,iPhosphate) = InitVal_PHO(polygon_m(i,j),1)
           bio_t(i,j,k,iSilice)    = InitVal_SIL(polygon_m(i,j),1)
         elseif(depth_t(i,j,k)>=Init_PRF(nmaxprof))then
           bio_t(i,j,k,iNitrate)   = InitVal_NIT(polygon_m(i,j),nmaxprof)
           bio_t(i,j,k,iPhosphate) = InitVal_PHO(polygon_m(i,j),nmaxprof)
           bio_t(i,j,k,iSilice)    = InitVal_SIL(polygon_m(i,j),nmaxprof)
         else
  
           K_pr=1
           do while(depth_t(i,j,k)>Init_PRF(K_pr))
             K_pr=K_pr+1
           enddo

           r1=-Init_PRF(K_pr-1)+depth_t(i,j,k)
           r2=-depth_t(i,j,k)+Init_PRF(K_pr)
           deltap=Init_PRF(K_pr)-Init_PRF(K_pr-1) 

           bio_t(i,j,k,iNitrate)  =(r1*InitVal_NIT(polygon_m(i,j),K_pr  )        &
                                   +r2*InitVal_NIT(polygon_m(i,j),K_pr-1))/deltap
           bio_t(i,j,k,iPhosphate)=(r1*InitVal_PHO(polygon_m(i,j),K_pr  )        &
                                   +r2*InitVal_PHO(polygon_m(i,j),K_pr-1))/deltap
           bio_t(i,j,k,iSilice)   =(r1*InitVal_SIL(polygon_m(i,j),K_pr  )        &
                                   +r2*InitVal_SIL(polygon_m(i,j),K_pr-1))/deltap
         endif



         if (                              &
         lon_t(i,j)>=9.5*deg2rad.and.                            & 
         lon_t(i,j)<12.5*deg2rad.and.                            &
         lat_t(i,j)>41.5*deg2rad.and.               &
         lat_t(i,j)<44.3*deg2rad)then  !.and.lat_t(i,j)<44.*deg2rad)then


         if(InitOC==1) then
         l=1

!         print*,'depth init1',-depth_t(i,j,k),Prof10y(l,1),Prof10y(l,nmax10y(l)),l,nmax10y(l)

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
!         print*,'depth init',-depth_t(i,j,k),Prof10y(l,K_pr)
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


!              print*,'i,j',i,j,var10y(l,K_pr-1,1),bio_t(i,j,k,izoonanoc)

              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


     !    bio_t(i,j,k,iNITRATE)= (0.0400278*(sal_t(i,j,k,1)**2)-57.0955)**2
     !    if(bio_t(i,j,k,iNITRATE)>=7.and.depth_t(i,j,k)<-500)           &
     !    bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513
     !    if(depth_t(i,j,k)<-800)                                        &
     !    bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513

     !    if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
     !    if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08097*(rhp_t(i,j,k)+rho-1000.)-30.8847)

     !    bio_t(i,j,k,iPHOSPHATE)= (-9.93055 + (0.00700501*sal_t(i,j,k,1)**2))**2
     !    if(bio_t(i,j,k,iPHOSPHATE)>=7.and.depth_t(i,j,k)<-500)           &
     !    bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
     !    if(depth_t(i,j,k)<-800)                                        &
     !    bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974

     !    if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
     !    if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.937182*(rhp_t(i,j,k)+rho-1000.)-29.4776)

     !    bio_t(i,j,k,iSILICE)= (-62.7026 + (0.0437829*sal_t(i,j,k,1)**2))**2
     !    if(bio_t(i,j,k,iSILICE)>=7.and.depth_t(i,j,k)<-500)           &
     !    bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
     !    if(depth_t(i,j,k)<-800)                                        &
     !    bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212

     !    if(bio_t(i,j,k,iSILICE)>=7.5) bio_t(i,j,k,iSILICE)=7.5
     !    if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.08691*(rhp_t(i,j,k)+rho-1000.)-31.0381)


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



         if(InitOC==1) then
         l=2

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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


         if(InitOC==1) then
         l=3

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC




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






! sicile
        elseif(                                  &
        lon_t(i,j)>=9.5*deg2rad.and.lon_t(i,j)<11.*deg2rad.and.         &
        lat_t(i,j)>36.5*deg2rad.and.lat_t(i,j)<=38.*deg2rad)then

         
         if(InitOC==1) then
         l=4

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC



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
 

         if(InitOC==1) then
         l=5

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC

      

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


         if(InitOC==1) then
         l=6

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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


         if(InitOC==1) then
         l=7

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC
       

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


         if(InitOC==1) then
         l=8

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC

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


         if(InitOC==1) then
         l=9

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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


         if(InitOC==1) then
         l=10

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC

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


         if(InitOC==1) then
         l=11

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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


         if(InitOC==1) then
         l=12

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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


         if(InitOC==1) then
         l=13

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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




!! gibraltar
        elseif(                                  &
        lon_t(i,j)<=-5.*deg2rad)then
!

         if(InitOC==1) then
         l=13

         !print*,'passe l=13'

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


!
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
         

         if(InitOC==1) then
         l=14

         if(-depth_t(i,j,k)<=Prof10y(l,1))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,1,1)
               bio_t(i,j,k,izoomicroc) = var10y(l,1,2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,1,3)
               bio_t(i,j,k,isynec)     = var10y(l,1,4)
               bio_t(i,j,k,inanoc)     = var10y(l,1,5)
               bio_t(i,j,k,idiac)      = var10y(l,1,6)
               bio_t(i,j,k,ibactc)     = var10y(l,1,7)
               bio_t(i,j,k,ismopc)     = var10y(l,1,8)
               bio_t(i,j,k,ilmopc)     = var10y(l,1,9)
               bio_t(i,j,k,imodc)      = var10y(l,1,10)

             ELSEif(-depth_t(i,j,k)>=Prof10y(l,nmax10y(l)))then
               bio_t(i,j,k,izoonanoc)  = var10y(l,nmax10y(l),1)
               bio_t(i,j,k,izoomicroc) = var10y(l,nmax10y(l),2)
               bio_t(i,j,k,izoomesoc)  = var10y(l,nmax10y(l),3)
               bio_t(i,j,k,isynec)     = var10y(l,nmax10y(l),4)
               bio_t(i,j,k,inanoc)     = var10y(l,nmax10y(l),5)
               bio_t(i,j,k,idiac)      = var10y(l,nmax10y(l),6)
               bio_t(i,j,k,ibactc)     = var10y(l,nmax10y(l),7)
               bio_t(i,j,k,ismopc)     = var10y(l,nmax10y(l),8)
               bio_t(i,j,k,ilmopc)     = var10y(l,nmax10y(l),9)
               bio_t(i,j,k,imodc)      = var10y(l,nmax10y(l),10)

             ELSE

              K_pr=1
              do while(-depth_t(i,j,k)>Prof10y(l,K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,izoonanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,1) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,1) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomicroc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,2) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,2) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,izoomesoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,3) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,3) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,isynec)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,4) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,4) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,inanoc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,5) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,5) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,idiac)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,6) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,6) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ibactc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,7) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,7) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ismopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,8) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,8) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,ilmopc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,9) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,9) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))


              bio_t(i,j,k,imodc)=(                              &
              ( Prof10y(l,K_pr  )+depth_t(i,j,k))*var10y(l,K_pr-1,10) &
             +(-Prof10y(l,K_pr-1)-depth_t(i,j,k))*var10y(l,K_pr  ,10) &
             )/(Prof10y(l,K_pr  )-Prof10y(l,K_pr-1))

           endif


         endif !InitOC


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

          if (InitOC.eq.1) then

!! Deduction of other Dia and Zoo from DiaChl
         bio_t(i,j,k,iDiaChl)  = bio_t(i,j,k,iDiaC)/        CoefDiaC        
         bio_t(i,j,k,iDiaN)  = (  minNCDia +  maxNCDia ) / 2.  &
       * bio_t(i,j,k,iDiaC)
         bio_t(i,j,k,iDiaP)  = (  minPCDia +  maxPCDia ) / 2.  &
       * bio_t(i,j,k,iDiaC)
         bio_t(i,j,k,iDiaSi) = ( minSiCDia + maxSiCDia ) / 2.  &
       * bio_t(i,j,k,iDiaC)


!! Deduction of other Nano from NanoChl
         bio_t(i,j,k,iNanoChl) = bio_t(i,j,k,iNanoC) /    CoefNanoC       
         bio_t(i,j,k,iNanoN) = (  minNCNano +  maxNCNano ) / 2. &
       * bio_t(i,j,k,iNanoC)
         bio_t(i,j,k,iNanoP) = (  minPCNano +  maxPCNano ) / 2. &
       * bio_t(i,j,k,iNanoC)


!! Deduction of other Syne from SyneChl
         bio_t(i,j,k,iSyneChl) = bio_t(i,j,k,iSyneC) /  CoefSyneC       
         bio_t(i,j,k,iSyneN) = (  minNCSyne +  maxNCSyne ) / 2. &
       * bio_t(i,j,k,iSyneC)
         bio_t(i,j,k,iSyneP) = (  minPCSyne +  maxPCSyne ) / 2. &
       * bio_t(i,j,k,iSyneC)


!! Deduction of MODP from MODC 
      bio_t(i,j,k,iMODN)   =  RedfieldNC * bio_t(i,j,k,iMODC)
      bio_t(i,j,k,iMODP)   =  RedfieldPC * bio_t(i,j,k,iMODC)


!! Deduction of SMOP and LMOP from SMOPC 
      bio_t(i,j,k,iSMOPN) = RedfieldNC * bio_t(i,j,k,iSMOPC)
      bio_t(i,j,k,iSMOPP) = RedfieldPC * bio_t(i,j,k,iSMOPC)
      bio_t(i,j,k,iSMOPChl) = bio_t(i,j,k,iSMOPC) /CoefSyneC

      bio_t(i,j,k,iLMOPN) = RedfieldNC * bio_t(i,j,k,iLMOPC)
      bio_t(i,j,k,iLMOPP) = RedfieldPC * bio_t(i,j,k,iLMOPC)
      bio_t(i,j,k,iLMOPSi) =  bio_t(i,j,k,iLMOPN)                

          endif



          bio_t(i,j,k,iOxygen)       =bio_t(i,j,k,iOxygen)  &
                            *(rhp_t(i,j,k)+rho) /1000.




      endif ! mask 
      enddo
      enddo
      enddo


      endif     ! ChoiceInitPelagic==4


      endif ! InitPelagic




! Rates of change
!----------------

      do  I = 1 , NumpelagicBoxes

      dZooNanoC(I)  = 0.
      dZooMicroC(I) = 0.
      dZooMesoC(I)  = 0.
      dSyneC(I)     = 0.
      dSyneN(I)     = 0.
      dSyneChl(I)   = 0.
      dSyneP(I)     = 0.
      dNanoC(I)     = 0.
      dNanoN(I)     = 0.
      dNanoChl(I)   = 0.
      dNanoP(I)     = 0.
      dDiaC(I)      = 0.
      dDiaN(I)      = 0.
      dDiaChl(I)    = 0.
      dDiaSi(I)     = 0.
      dDiaP(I)      = 0.
      dBactC(I)     = 0.
      dS_MOPC(I)    = 0.
      dS_MOPN(I)    = 0.
      dS_MOPP(I)    = 0.
      dS_MOPChl(I)  = 0.
      dS_MOPSi(I)   = 0.
      dL_MOPC(I)    = 0.
      dL_MOPN(I)    = 0.
      dL_MOPSi(I)   = 0.
      dL_MOPP(I )   = 0.
      dMODC(I)      = 0.
      dMODN(I)      = 0.
      dMODP(I)      = 0.
      dNitrate(I)   = 0.
      dAmmonium(I)  = 0.
      dSilice(I)    = 0.
      dPhosphate(I) = 0.
      dOxygen(I)    = 0.
!     dODU(I)       = 0.

      enddo


! Flux from Eco3m to Diagnostics
!-------------------------------

      do  I = 1 , NumpelagicBoxes

      PPBi(I)   = 0.
      GrazCi(I) = 0.
      RespCi(I) = 0.
      ExuCi(I)  = 0.
      MortCi(I) = 0.

      enddo

!=====================================================================*
! 3. Initialisation of the variables                                  *
! END                                                                 *
!=====================================================================*


!=====================================================================*
! 3,5. Definition of different coordinates (i,j) for diagnostics      *
! BEGINNING                                                           *
!=====================================================================*

!Point Dyfamed
      longi1=7.87*deg2rad
      latit1 =43.42*deg2rad
!      call latlon_to_ij('glb')
      call latlontoij(longi1,latit1,'glb')
      idyf=NINT(deci)
      jdyf=NINT(decj)
!      print*,'passe ici',deci,decj,longi1,latit1

 
      if((idyf-par%timax(1)).ge.1.and.(idyf-par%timax(1)).le.imax.and.    &
         (jdyf-par%tjmax(1)).ge.1.and.(jdyf-par%tjmax(1)).le.jmax)   then
       print*,'dyfamed',idyf,jdyf,idyf-par%timax(1),jdyf-par%tjmax(1),                                    &
            depth_w(idyf-par%timax(1),jdyf-par%tjmax(1),1),                &
            mask_t(idyf-par%timax(1),jdyf-par%tjmax(1),kmin_w(idyf-par%timax(1),jdyf-par%tjmax(1)))

!      stop'dans InitPelagic'
      endif

      idyf3=NINT(deci)
      jdyf3=NINT(decj)-2


!      print*,'passe ici',deci,decj,longi1,latit1

!Point 1 d'une zone entourant Dyfamed
      longi1=7.7*deg2rad
      latit1 =43.3*deg2rad
      call latlon_to_ij('glb')
      idyf1=NINT(deci)
      jdyf1=NINT(decj)

!Point 2 d'une zone entourant Dyfamed
      longi1=8*deg2rad
      latit1 =43.5*deg2rad
      call latlon_to_ij('glb')
      idyf2=NINT(deci)
      jdyf2=NINT(decj)

!Point GDL 1
      longi1=3.63*deg2rad
      latit1 =43.1*deg2rad
      call latlon_to_ij('glb')
      igdl1=NINT(deci)
      jgdl1=NINT(decj)

!Point GDL 2
      longi1=3.4*deg2rad
      latit1 =42.8*deg2rad
      call latlon_to_ij('glb')
      igdl2=NINT(deci)
      jgdl2=NINT(decj)




!18 Points de la grille Symphonie les plus proches 
!des points de mesure BIO 
!par prelevements bouteilles CASCADE
! 1 à 9 radiale L
! 10 à 18 radiale M

!CL010
      longi1=5.501*deg2rad
      latit1 =42.0349*deg2rad
      call latlon_to_ij('glb')
      icas1=NINT(deci)
      jcas1=NINT(decj)
!      print*,'passe ici2',deci,decj,longi1,latit1   
      

!CL12_1
      longi1=5.9039*deg2rad
      latit1 =42.036*deg2rad
      call latlon_to_ij('glb')
      icas2=NINT(deci)
      jcas2=NINT(decj)
!CL12_2
      longi1=5.9039*deg2rad
      latit1 =42.036*deg2rad
      call latlon_to_ij('glb')
      icas3=NINT(deci)
      jcas3=NINT(decj)

!CL01
      longi1=3.4883*deg2rad
      latit1 =42.0371*deg2rad
      call latlon_to_ij('glb')
      icas4=NINT(deci)
      jcas4=NINT(decj)
!CL03
      longi1=3.883*deg2rad
      latit1 =42.0338*deg2rad
      call latlon_to_ij('glb')
      icas5=NINT(deci)
      jcas5=NINT(decj)

!CL05
      longi1=4.2892*deg2rad
      latit1 =42.0319*deg2rad
      call latlon_to_ij('glb')
      icas6=NINT(deci)
      jcas6=NINT(decj)

!CL08
      longi1=5.0982*deg2rad
      latit1 =42.0324*deg2rad
      call latlon_to_ij('glb')
      icas7=NINT(deci)
      jcas7=NINT(decj)

!CS2400_1
      longi1=4.6955*deg2rad
      latit1 =42.0286*deg2rad
      call latlon_to_ij('glb')
      icas8=NINT(deci)
      jcas8=NINT(decj)

!CS2400_2
      longi1=4.692*deg2rad
      latit1 =42.0372*deg2rad
      call latlon_to_ij('glb')
      icas9=NINT(deci)
      jcas9=NINT(decj)

!CM01_1
      longi1=4.6988*deg2rad
      latit1 =41.1371*deg2rad
      call latlon_to_ij('glb')
      icas10=NINT(deci)
      jcas10=NINT(decj)

!CM03
      longi1=4.6927*deg2rad
      latit1 =41.4363*deg2rad
      call latlon_to_ij('glb')
      icas11=NINT(deci)
      jcas11=NINT(decj)

!CM05
      longi1=4.6981*deg2rad
      latit1 =41.738*deg2rad
      call latlon_to_ij('glb')
      icas12=NINT(deci)
      jcas12=NINT(decj)

!CM08_01
      longi1=4.7009*deg2rad
      latit1 =42.3303*deg2rad
      call latlon_to_ij('glb')
      icas13=NINT(deci)
      jcas13=NINT(decj)

!CM08_03
      longi1=4.7009*deg2rad
      latit1 =42.3303*deg2rad
      call latlon_to_ij('glb')
      icas14=NINT(deci)
      jcas14=NINT(decj)

!CM10
      longi1=4.6947*deg2rad
      latit1 =42.6294*deg2rad
      call latlon_to_ij('glb')
      icas15=NINT(deci)
      jcas15=NINT(decj)

!CM12
      longi1=4.7002*deg2rad
      latit1 =42.9310*deg2rad
      call latlon_to_ij('glb')
      icas16=NINT(deci)
      jcas16=NINT(decj)

!CS2400_7
      longi1=4.692*deg2rad
      latit1 =42.0372*deg2rad
      call latlon_to_ij('glb')
      icas17=NINT(deci)
      jcas17=NINT(decj)

!CS2400_9
      longi1=4.6955*deg2rad
      latit1 =42.0286*deg2rad
      call latlon_to_ij('glb')
      icas18=NINT(deci)
      jcas18=NINT(decj)

! points de mesure MOSE 2011 sur 2radiales


      open(unit=1,file='LON_LAT_MOOSE2011.txt')
      
      do i=1,23
      read(1,*)lon,lat

      longi1=lon*deg2rad
      latit1 =lat*deg2rad
      call latlon_to_ij('glb')
      imoose(i)=NINT(deci)
      jmoose(i)=NINT(decj)
      enddo
      close(1)

!Point MEDOC
      longi1=5.*deg2rad
      latit1 =42.*deg2rad
!      longi1=4.6833*deg2rad  !diag17 point pour comparaison piege marion stabhloz-xavier
!      latit1 =42.41667*deg2rad !diag17 point pour comparaison piege marion stabhloz-xavier
      call latlon_to_ij('glb')
      imedoc=NINT(deci)
      jmedoc=NINT(decj)

      longi1=4.316667*deg2rad  !point pour comparaison piege marion stabhloz-xavier
      latit1 =42.116667*deg2rad !point pour comparaison piege marion stabhloz-xavier
      call latlon_to_ij('glb')
      isw2060=NINT(deci)
      jsw2060=NINT(decj)

      longi1=4.35*deg2rad  !point pour comparaison piege marion stabhloz-xavier
      latit1 =42.25*deg2rad !point pour comparaison piege marion stabhloz-xavier
      call latlon_to_ij('glb')
      isc2160=NINT(deci)
      jsc2160=NINT(decj)

      longi1=4.55*deg2rad  !point pour comparaison piege marion stabhloz-xavier
      latit1 =42.166667*deg2rad !point pour comparaison piege marion stabhloz-xavier
      call latlon_to_ij('glb')
      isc2240=NINT(deci)
      jsc2240=NINT(decj)

!Point Villefranche SOMLIT
      longi1=7.32*deg2rad
      latit1 =43.68*deg2rad
      call latlon_to_ij('glb')
      ivill=NINT(deci)
      jvill=NINT(decj)
      if((ivill-par%timax(1)).ge.1.and.(ivill-par%timax(1)).le.imax.and.    &
         (jvill-par%tjmax(1)).ge.1.and.(jvill-par%tjmax(1)).le.jmax)   &
       print*,'villefranche',ivill,jvill,                                    &
            depth_w(ivill-par%timax(1),jvill-par%tjmax(1),1),                &
            mask_t(ivill-par%timax(1),jvill-par%tjmax(1),kmin_w(ivill-par%timax(1),jvill-par%tjmax(1)))

!      ivill=269
!      jvill=319
!      if((ivill-par%timax(1)).ge.1.and.(ivill-par%timax(1)).le.imax.and.    &
!         (jvill-par%tjmax(1)).ge.1.and.(jvill-par%tjmax(1)).le.jmax)   &
!      print*,'villefranche',ivill,jvill,                                    &
!            depth_w(ivill-par%timax(1),jvill-par%tjmax(1),1),                &
!            mask_t(ivill-par%timax(1),jvill-par%tjmax(1),kmin_w(ivill-par%timax(1),jvill-par%tjmax(1)))

!      ivill=271
!      jvill=319
!      if((ivill-par%timax(1)).ge.1.and.(ivill-par%timax(1)).le.imax.and.    &
!         (jvill-par%tjmax(1)).ge.1.and.(jvill-par%tjmax(1)).le.jmax)   &
!      print*,'villefranche',ivill,jvill,                                    &
!            depth_w(ivill-par%timax(1),jvill-par%tjmax(1),1),                &
!            mask_t(ivill-par%timax(1),jvill-par%tjmax(1),kmin_w(ivill-par%timax(1),jvill-par%tjmax(1)))
!
!      ivill=269
!      jvill=318
!      if((ivill-par%timax(1)).ge.1.and.(ivill-par%timax(1)).le.imax.and.    &
!         (jvill-par%tjmax(1)).ge.1.and.(jvill-par%tjmax(1)).le.jmax)   &
!      print*,'villefranche',ivill,jvill,                                    &
!            depth_w(ivill-par%timax(1),jvill-par%tjmax(1),1),                &
!            mask_t(ivill-par%timax(1),jvill-par%tjmax(1),kmin_w(ivill-par%timax(1),jvill-par%tjmax(1)))


!Point Marseille SOMLIT
      longi1=5.28*deg2rad
      latit1 =43.23*deg2rad
      call latlon_to_ij('glb')
      imars=NINT(deci)
      jmars=NINT(decj)
      if((imars-par%timax(1)).ge.1.and.(imars-par%timax(1)).le.imax.and.    &
         (jmars-par%tjmax(1)).ge.1.and.(jmars-par%tjmax(1)).le.jmax)   &
      print*,'marseille',imars,jmars,                                    &
            depth_w(imars-par%timax(1),jmars-par%tjmax(1),1)
!     & ,mask_t(imars,jmars,kmax)

!Point Banyuls SOMLIT
      longi1=3.13*deg2rad
      latit1 =42.48*deg2rad
      call latlon_to_ij('glb')
      iban=NINT(deci)
      jban=NINT(decj)
      if((iban-par%timax(1)).ge.1.and.(iban-par%timax(1)).le.imax.and.    &
         (jban-par%tjmax(1)).ge.1.and.(jban-par%tjmax(1)).le.jmax)   &
      print*,'banyuls',iban,jban,                                    &
            depth_w(iban-par%timax(1),jban-par%tjmax(1),1),                &
            mask_t(iban-par%timax(1),jban-par%tjmax(1),kmin_w(iban-par%timax(1),jban-par%tjmax(1)))
!     & ,mask_t(iban,jban,kmax)
!      iban=131
!      jban=264
!      if((iban-par%timax(1)).ge.1.and.(iban-par%timax(1)).le.imax.and.    &
!         (jban-par%tjmax(1)).ge.1.and.(jban-par%tjmax(1)).le.jmax)   &
!      print*,'banyuls',iban,jban,                                    &
!            depth_w(iban-par%timax(1),jban-par%tjmax(1),1),                &
!            mask_t(iban-par%timax(1),jban-par%tjmax(1),kmin_w(iban-par%timax(1),jban-par%tjmax(1)))
!      iban=130
!      jban=263
!      if((iban-par%timax(1)).ge.1.and.(iban-par%timax(1)).le.imax.and.    &
!         (jban-par%tjmax(1)).ge.1.and.(jban-par%tjmax(1)).le.jmax)   &
!      print*,'banyuls',iban,jban,                                    &
!            depth_w(iban-par%timax(1),jban-par%tjmax(1),1)

!      stop'points somlit'

!Point in the Northern-Current 1 : Large of Villefranche
      longi1=7.32*deg2rad
      latit1 =43.51*deg2rad
      call latlon_to_ij('glb')
      icn1=NINT(deci)
      jcn1=NINT(decj)

!Point in the Northern-Current 2 : Upwelling
      longi1=5.94*deg2rad
      latit1 =42.92*deg2rad
      call latlon_to_ij('glb')
      icn2=NINT(deci)
      jcn2=NINT(decj)

!Entire domain
      longi1=0.*PI/180.
      latit1 =40.*PI/180.
      call latlon_to_ij('glb')
      i1domain=NINT(deci)
      j1domain=NINT(decj)
      longi1=9.18*PI/180.
      latit1 =44.50*PI/180.
      call latlon_to_ij('glb')
      i2domain=NINT(deci)
      j2domain=NINT(decj)
      i1domain=1
      i2domain=iglb
      j1domain=1
      j2domain=jglb



!GDL                  
      longi1=3.*deg2rad
      latit1 =42.30*deg2rad
      call latlon_to_ij('glb')
      i1gdl=NINT(deci)
      j1gdl=NINT(decj)
      longi1=6.*deg2rad
      latit1 =43.65*deg2rad
      call latlon_to_ij('glb')
      i2gdl=NINT(deci)
      j2gdl=NINT(decj)

!MEDOC box
      longi1=4.*deg2rad
      latit1 =41.*deg2rad
      call latlon_to_ij('glb')
      i1medoc=NINT(deci)
      j1medoc=NINT(decj)
      longi1=6.*deg2rad
      latit1 =42.50*deg2rad
      call latlon_to_ij('glb')
      i2medoc=NINT(deci)
      j2medoc=NINT(decj)

!LIGURE box
      longi1=7.*deg2rad
      latit1 =42.*deg2rad
      call latlon_to_ij('glb')
      i1ligure=NINT(deci)
      j1ligure=NINT(decj)
      longi1=9.*deg2rad
      latit1 =43.50*deg2rad
      call latlon_to_ij('glb')
      i2ligure=NINT(deci)
      j2ligure=NINT(decj)


! Diag box (same as MLD caculation) 
      longi1=2.*deg2rad
      latit1 =40.*deg2rad
      call latlon_to_ij('glb')
      i1diag=NINT(deci)
      j1diag=NINT(decj)
      longi1=9.*deg2rad
      latit1 =45.*deg2rad
      call latlon_to_ij('glb')
      i2diag=NINT(deci)
      j2diag=NINT(decj)

!Frontiere Ouest d'un domaine "LIGURE"
      longi1=6.18*deg2rad
      latit1 =42.*deg2rad
      call latlon_to_ij('glb')
      ilig=NINT(deci)

!Western boundary
      longi1=2.*deg2rad
      call latlon_to_ij('glb')
      iobwest=NINT(deci)

!Eastern boundary
      longi1=9.*deg2rad
      call latlon_to_ij('glb')
      iobeast=NINT(deci)

!Southern boundary
      latit1=40.*deg2rad
      call latlon_to_ij('glb')
      jobsouth=NINT(decj)

!Southern boundary
      latit1=45.*deg2rad
      call latlon_to_ij('glb')
      jobnorth=NINT(decj)

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=4,file='messages',position='append')
      write(4,*)'-----------------------------------------------------'
      write(4,*)'subroutine InitPelagic'
      write(4,*)'Points de diag bio (i,j) dans le repere global'
      write(4,*)'idyf,jdyf',idyf,jdyf
      write(4,*)'igdl1,jgdl1',igdl1,jgdl1
      write(4,*)'igdl2,jgdl2',igdl2,jgdl2
      write(4,*)'imedoc,jmedoc',imedoc,jmedoc
      write(4,*)'ivill,jvill',ivill,jvill
      write(4,*)'imars,jmars',imars,jmars
      write(4,*)'iban,jban',iban,jban
      write(4,*)'icn1,jcn1',icn1,jcn1
      write(4,*)'icn2,jcn2',icn2,jcn2
      write(4,*)'i1domain,j1domain',i1domain,j1domain
      write(4,*)'i2domain,j2domain',i2domain,j2domain
      write(4,*)'i1gdl,j1gdl',i1gdl,j1gdl
      write(4,*)'i2gdl,j2gdl',i2gdl,j2gdl
      write(4,*)'i1medoc,j1medoc',i1medoc,j1medoc
      write(4,*)'i2medoc,j2medoc',i2medoc,j2medoc
      write(4,*)'i1ligure,j1ligure',i1ligure,j1ligure
      write(4,*)'i2ligure,j2ligure',i2ligure,j2ligure
      write(4,*)'ilig',ilig
      write(4,*)'iobwest',iobwest
      write(4,*)'iobeast',iobeast
      write(4,*)'jobsouth',jobsouth
      write(4,*)'isw2060,jsw2060',isw2060,jsw2060
      write(4,*)'isc2160,jsc2160',isc2160,jsc2160
      write(4,*)'isc2240,jsc2240',isc2240,jsc2240
      close(4)
      endif !#mpi-->>-->


!=====================================================================*
! 3,5. Definition des coordonnees (i,j) des points de diagnostiques   *
! END                                                                 *
!=====================================================================*


!=====================================================================*
! 4. Writes the initial profiles for the Femme Environment            *
! BEGINNING                                                           *
!=====================================================================*

      if(InitFemme==1) then    

      NameVar(iZooNanoC)='ZooNanoC'
      NameVar(iZooMicroC)='ZooMicroC'
      NameVar(iZooMesoC)='ZooMesoC'
      NameVar(iSyneC)='SyneC'
      NameVar(iSyneN)='SyneN'
      NameVar(iSyneP)='SyneP'
      NameVar(iSyneChl)='SyneChl'
      NameVar(iNanoC)='NanoC'
      NameVar(iNanoN)='NanoN'
      NameVar(iNanoP)='NanoP'
      NameVar(iNanoChl)='NanoChl'
      NameVar(iDiaC)='DiaC'
      NameVar(iDiaN)='DiaN'
      NameVar(iDiaP)='DiaP'
      NameVar(iDiaChl)='DiaChl'
      NameVar(iDiaSi)='DiaSi'
      NameVar(iBactC)='BactC'
      NameVar(iSMOPC)='S_MOPC'
      NameVar(iSMOPN)='S_MOPN'
      NameVar(iSMOPP)='S_MOPP'
      NameVar(iSMOPCHL)='S_MOPChl'
      NameVar(iSMOPSi)='S_MOPSi'
      NameVar(iLMOPC)='L_MOPC'
      NameVar(iLMOPN)='L_MOPN'
      NameVar(iLMOPP)='L_MOPP'
      NameVar(iLMOPSi)='L_MOPSi'
      NameVar(iMODC)='MODC'
      NameVar(iMODN)='MODN'
      NameVar(iMODP)='MODP'
      NameVar(iNitrate)='Nitrate'
      NameVar(iAmmonium)='Ammonium'
      NameVar(iPhosphate)='Phosphate'
      NameVar(iSilice)='Silice'
# ifdef mes
      NameVar(iMIP)='MIP'
# endif 

      i1=ifemme-par%timax(1)
      j1=jfemme-par%tjmax(1)
      if(i1>=1.and.i1<=imax-2.and.j1>=1.and.j1<=jmax-2)then
      open(unit=3,file='Init_Dyfamed0304')

      write(3,*)
      write(3,*)
      write(3,*)'#####'
      write(3,'(A24)')'@name    Depth    @value'
      write(3,*)


      do vb=1,vbmax 
      do k=1,kmax
      write(3,*)NameVar(vb),-depth_t(i1,j1,k)       &
                                ,bio_t(i1,j1,k,vb)
      enddo
      enddo
      close(3)
      endif

! Verif profils regression
      i1=147-par%timax(1)
      j1=32-par%tjmax(1)
      if(i1>=1.and.i1<=imax-2.and.j1>= 1.and.j1<=jmax-2)then
       open(unit=3,file='VerifReg_Interieur')
      do k=1,kmax
      write(3,*)-depth_t(i1,j1,k)                           &
               ,bio_t(i1,j1,k,iNitrate),bio_t(i1,j1,k,iPhosphate)
      enddo
      close(3)
      endif

      i1=173-par%timax(1)
      j1=32-par%tjmax(1)
      if(i1>=1.and.i1<=imax-2.and.j1>= 1.and.j1<=jmax-2)then
       open(unit=3,file='VerifReg_FrontiereNR')
      do k=1,kmax
      write(3,*)-depth_t(i1,j1,k)                                 &
               ,bio_t(i1,j1,k,iNitrate),bio_t(i1,j1,k,iPhosphate) &
               ,(18.1+15*exp(depth_t(i1,j1,k)/300.))
      enddo
      close(3)
      endif

      i1=173-par%timax(1)
      j1=49-par%tjmax(1)
      if(i1>=1.and.i1<=imax-2.and.j1>= 1.and.j1<=jmax-2)then
       open(unit=3,file='VerifReg_FrontiereR')
      do k=1,kmax
      write(3,*)-depth_t(i1,j1,k)                                 &
               ,bio_t(i1,j1,k,iNitrate),bio_t(i1,j1,k,iPhosphate)
      enddo
      close(3)
      endif

!      STOP'dans Initpelagic'


! Outputs for Matlab
! on lit juste notebook_graph pour connaitre le nom du repertoire des graphiques
      if(par%rank==0) print*,'Avant nomfichier InitPelagic'
      TEXTE90=nomfichier(21)
      open(unit=3,file=TEXTE90)
      read(3,*)
      read(3,*)
      read(3,'(A)')dirgraph
      close(3)
      do k=1,90
      if(dirgraph(K:K)==' ') then
          lname4=K-1
          GOTO 23
      endif
      enddo
   23 CONTINUE
!      if(dirgraph(lname4:lname4)=='/')lname4=lname4-1
!      lname4=lname4+1
!      dirgraph(lname4:lname4)='/'
!
!      TEXTE90=dirgraph(1:lname4)//'Prof'
!
!      if(par%rank==0)then
!      open(unit=3,file=TEXTE90) 
!      do k=1,kmax
!      write(3,*)-depth_t(idyf,jdyf,K)
!      enddo
!      close(3)
!      endif
!      
!!Medoc
!      i1=181-par%timax(1)
!      j1=157-par%tjmax(1)
!      if(i1>=1.and.i1<=mbio-2.and.j1>= 1.and.j1<=nbio-2)then
!      do vb=1,vbmax
!      FicInitVar=dirgraph(1:lname4)//'Ini1'//NameVar(vb)
!!      write(6,*)'FicInitVar: ',FicInitVar     
!      open(unit=3,file=FicInitVar) 
!      do k=1,kmax
!        write(3,*)-depth_t(i1,j1,k),bio_t(i1,j1,K,vb)
!!          -depth_t(I2,J2,K),bio_t(I2,J2,K,vb)
!      enddo
!      close(3)
!      enddo
!      endif
!
!!Dyf
!      i1=276-par%timax(1)
!      j1=222-par%tjmax(1)
!      if(i1>=1.and.i1<=mbio-2.and.j1>= 1.and.j1<=nbio-2)then
!      do vb=1,vbmax
!      FicInitVar=dirgraph(1:lname4)//'Ini1'//NameVar(vb)
!!      write(6,*)'FicInitVar: ',FicInitVar
!      open(unit=3,file=FicInitVar)
!      do k=1,kmax
!        write(3,*)-depth_t(i1,j1,k),bio_t(i1,j1,K,vb)
!!          -depth_t(I2,J2,K),bio_t(I2,J2,K,vb)
!      enddo
!      close(3)
!      enddo
!      endif
!
!
!! NC
!      i1=212-par%timax(1)
!      j1=199-par%tjmax(1)
!      if(i1>=1.and.i1<=mbio-2.and.j1>= 1.and.j1<=nbio-2)then
!      do vb=1,vbmax
!      FicInitVar=dirgraph(1:lname4)//'Ini2'//NameVar(vb)
!!      write(6,*)'FicInitVar: ',FicInitVar
!      open(unit=3,file=FicInitVar)
!      do k=1,kmax
!        write(3,*)-depth_t(i1,j1,k),bio_t(i1,j1,K,vb)
!      enddo
!      close(3)
!      enddo
!      endif
!
!! NC
!      i1=250-par%timax(1)
!      j1=221-par%tjmax(1)
!      if(i1>=1.and.i1<=mbio-2.and.j1>= 1.and.j1<=nbio-2)then
!      do vb=1,vbmax
!      FicInitVar=dirgraph(1:lname4)//'Ini2'//NameVar(vb)
!!      write(6,*)'FicInitVar: ',FicInitVar
!      open(unit=3,file=FicInitVar)
!      do k=1,kmax
!        write(3,*)-depth_t(i1,j1,k),bio_t(i1,j1,K,vb)
!      enddo
!      close(3)
!      enddo
!      endif
!
!
!! GDL : 2 points
!      i1=135-par%timax(1)
!      j1=207-par%tjmax(1)
!      if(i1>=1.and.i1<=mbio-2.and.j1>= 1.and.j1<=nbio-2)then
!      do vb=1,vbmax
!      FicInitVar=dirgraph(1:lname4)//'Ini3'//NameVar(vb)
!!      write(6,*)'FicInitVar: ',FicInitVar
!      open(unit=3,file=FicInitVar)
!      do k=1,kmax
!        write(3,*)-depth_t(i1,j1,k),bio_t(i1,j1,K,vb)        
!      enddo
!      close(3)
!      enddo
!      endif
!
!      i1=130-par%timax(1)
!      j1=190-par%tjmax(1)
!      if(i1>=1.and.i1<=mbio-2.and.j1>= 1.and.j1<=nbio-2)then
!      do vb=1,vbmax
!      FicInitVar=dirgraph(1:lname4)//'Ini3'//NameVar(vb)
!      open(unit=3,file=FicInitVar,access='append')  ! pour le coller au precedent si celui
!! ci a été créé avec le meme proc
!      do k=1,kmax
!        write(3,*)-depth_t(i1,j1,k),bio_t(i1,j1,K,vb)
!      enddo
!      close(3)
!      enddo
!      endif

      endif

!=====================================================================*
! 4. Writes the initial profiles for the Femme Environment            *
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

#ifdef parallele
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
#endif


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
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

      END SUBROUTINE InitPelagic
