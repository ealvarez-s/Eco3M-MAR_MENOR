      SUBROUTINE InitBenthic

!_____________________________________________________________________*
! 3d ecosystem model                                                  *
!                                                                     *
! LAST REVISION: 5 NOVEMBER 2009                                      *
!                                                                     *
! Implementation: Caroline Ulses                                      *
!                 LA                                                  *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! 1. Initialises the benthic variables                                *
!_____________________________________________________________________*

!_____________________________________________________________________*
!                                                                     *
! Modifications:                                                      *
!                                                                     *
!_____________________________________________________________________*



!---------------------------------------------------------------------*
! Declarations:


! Global variables

      USE ModuleComBlock
      USE UserDeclaration
      USE ModuleDeclaration
      USE MODULE_PRINCIPAL
      use module_parallele !#mpi
      use module_global
      IMPLICIT NONE


! Local variables
      INTEGER NumVar(3),LDir,NRECBENTHstart,nrecBenthend,LRECB,NRECBENTH


       real*4,dimension(0:imax+1,0:jmax+1) ::                            &
           SedStart,SedEnd

      CHARACTER DirInitBenthicFile*90,FicBenthic(4)*82


!     REAL ANYVAR2DB(MECO,NECO,4)
       DOUBLE PRECISION NCBDet

!      DOUBLE PRECISION NCBDet,SedStart(mbio2-mbio1+1,nbio2-nbio1+1), &
!                              SedEnd(mbio2-mbio1+1,nbio2-nbio1+1)

!---------------------------------------------------------------------*

       LRECB=(MBIO2-MBIO1+1)*(NBIO2-NBIO1+1)*4

!=====================================================================*
! 1. Initialisation of the variables                                  *
! BEGINNING                                                           *
!=====================================================================*


! Reading of the notebook

      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read'       &
                        ,NOMFICHIER(30)

      OPEN(UNIT=3,FILE=NOMFICHIER(30)) ! open the notebook

      DO K=1,11
      READ(3,*)
      ENDDO 

      READ(3,*)IBenthic
      READ(3,*)

      READ(3,*)
      READ(3,'(A)')DirInitBenthicFile 
      READ(3,*)NRECBENTHStart
      READ(3,*)NRECBENTHEnd

      READ(3,*)
      READ(3,*)

      READ(3,*)NumBDET
      READ(3,*)BTfunc
      READ(3,*)
      READ(3,*)DecayRate
      READ(3,*)cdet20SED
      READ(3,*)
      READ(3,*)kFast
      READ(3,*)kSlow
      READ(3,*)pFast
      READ(3,*)cFdet20SED
      READ(3,*)cSdet20SED
      READ(3,*)BQ10
      READ(3,*)NCrFD
      READ(3,*)NCrSD

      READ(3,*)
      READ(3,*)
      READ(3,*)pPMin
      READ(3,*)pSiMin

      READ(3,*)
      READ(3,*)

      READ(3,*)pNit
      READ(3,*)pdeNit

      READ(3,*)
      READ(3,*)

      READ(3,*)O2bw

      CLOSE(3) 

      IF(IBenthic.EQ.2.OR.IBenthic.EQ.3) THEN  !================================>

!     print*,'passe IBenthic',MBIO1,NBIO1,MBIO2,NBIO2       
!     print*,'passe IBenthic',IBenthic                      


! directory of the initial element sediment content files
       DO 20 K=1,90
         IF(DirInitBenthicFile(K:K).EQ.' ') THEN !-------->
          LDir=K-1
          GOTO 21
         ENDIF                         !-------->
   20  CONTINUE
   21  CONTINUE
       IF(DirInitBenthicFile(LDir:LDir).EQ.'/')LDir=LDir-1
       LDir=LDir+1
       DirInitBenthicFile(LDir:LDir)='/'

     print*,'pres a lire FicBenthic(1)'           !!!!! fay

      WRITE(FicBenthic(1),'(A,A17)')         &                         
                 DirInitBenthicFile(1:LDIR),'sedC'//dom_c//'.binrec'
      WRITE(FicBenthic(2),'(A,A17)')         &
                 DirInitBenthicFile(1:LDIR),'sedN'//dom_c//'.binrec'
      WRITE(FicBenthic(3),'(A,A17)')         &
                 DirInitBenthicFile(1:LDIR),'sedP'//dom_c//'.binrec'
      WRITE(FicBenthic(4),'(A,A18)')         &
                 DirInitBenthicFile(1:LDIR),'sedSI'//dom_c//'.binrec'  


!! Initialisation of N, P, Si content in the sediment
      DO K2=1,4
      OPEN(UNIT=3,FILE=FicBenthic(K2)          &
                 ,ACCESS='DIRECT'              &
                 ,RECL=LRECB                   &
                 ,FORM='UNFORMATTED') 
      write(6,*)'Fic',K2,FicBenthic(K2),NRECBENTH,LRECB,    &
      	size(anyvar2d(MBIO1:MBIO2,NBIO1:NBIO2))

! On lit le depot annuel : donnees en mmol/m2 *86400
      READ(3,REC=nrecBenthStart)  &
       SedStart(MBIO1:MBIO2,NBIO1:NBIO2)
      READ(3,REC=nrecBenthEnd)  &
       SedEnd(MBIO1:MBIO2,NBIO1:NBIO2)
      
! On en deduit la teneur en matiere organique dans le sediment
      DO J=NBIO1,NBIO2
      DO I=MBIO1,MBIO2
      if(k2==1)  &
      CBDet(I,J)  = (SedEnd(i,j)-SedStart(i,j))  &
            /((nrecBenthEnd-nrecBenthStart)*nint(86400/dti_fw)*dti_fw/86400) &
            /DecayRate   ! en mmolC/m2 !11/7/13

      if(k2==2)  &
      NBDet(I,J)  = (SedEnd(i,j)-SedStart(i,j))  &
            /((nrecBenthEnd-nrecBenthStart)*nint(86400/dti_fw)*dti_fw/86400) &
            /DecayRate   ! en mmolN/m2 !11/7/13

      if(k2==3)  &
      PBDet(I,J)  = (SedEnd(i,j)-SedStart(i,j))  &
            /((nrecBenthEnd-nrecBenthStart)*nint(86400/dti_fw)*dti_fw/86400) &
            /DecayRate ! en mmolP/m2 !11/7/13

      if(k2==4)  &
      SiBDet(I,J) = (SedEnd(i,j)-SedStart(i,j))  &
            /((nrecBenthEnd-nrecBenthStart)*nint(86400/dti_fw)*dti_fw/86400) &
            /DecayRate ! en mmolSi/m2 !11/7/13
      ENDDO
      ENDDO
      ENDDO    !K2
      CLOSE(3)


      DO J=NBIO1,NBIO2
      DO I=MBIO1,MBIO2
      IF (CBDet(I,J) .GT. 0.) THEN
        NCBDet = NBDet(I,J) / CBDet(I,J)  ! N/C ratio of deposited detritus
      ENDIF

      IF (NCBDet .GT. NCrFD) THEN            ! All deposition as fast detritus
        CBFDet(I,J) = CBDet(I,J)
      ELSEIF (NCBDet .LE. NCrFD .AND. NCBDet .GT. NCrSD) THEN    ! Both SDETand FDET are deposited
        CBFDet(I,J) = (NCBDet - NCrSD)/(NCrFD-NCrSD)*CBDet(I,J)
        CBSDet(I,J) = CBDet(I,J)  - CBFDet(I,J)
      ELSE       ! Only SDET
        CBSDet(I,J) = CBDet(I,J)
      ENDIF

      IF (CBFDet(I,J) .LT.0) CBFDet(I,J) = 0.
      IF (CBSDet(I,J) .LT.0) CBSDet(I,J) = 0.



      ENDDO
      ENDDO

      ENDIF !========================================================>

!=====================================================================*
! 3. Initialisation of the variables                                  *
! END                                                                 *
!=====================================================================*

!=====================================================================*
! 5. Initialisation of output variables                               *
! BEGINNING                                                           *
!=====================================================================*
! mise a zero des diagnostics 0D

      SUM_NITRIF = 0.
      SUM_DENITRIF = 0.
      SUM_NO3Flux = 0.
      SUM_NH3Flux = 0.
      NitrifSMean   = 0.
      DenitrifSMean = 0.
      NO3FluxSMean  = 0.
      NH3FluxSMean  = 0.
      SUM_NITRIF_est = 0.
      SUM_DENITRIF_est = 0.
      SUM_NO3Flux_est = 0.
      SUM_NH3Flux_est = 0.
      NitrifSMean_est   = 0.
      DenitrifSMean_est = 0.
      NO3FluxSMean_est  = 0.
      NH3FluxSMean_est  = 0.

      SUMT_BENT = 0.
      TPS_BENT = 0.

      DO J=MAX0(NBIO1,1),NBIO2
      DO I=MAX0(MBIO1,1),MBIO2
      SUMT_NITRIF(I,J) = 0.
      SUMT_DENITRIF(I,J)= 0.
      SUMT_NH3Flux(I,J) = 0.
      SUMT_NO3Flux(I,J) = 0.
      ENDDO
      ENDDO 


!=====================================================================*
! 5. Initialisation of output variables                               *
! END                                                                 *
!=====================================================================*

!     print*,'SUMAREAGDL_plat',SUMAREAGDL_plat
!     print*,'SUMAREATOTAL',SUMAREATOTAL
       if(par%rank==0)write(6,*)'Passe par InitBenthic' 

!       STOP'dans InitBenthic'

      END SUBROUTINE InitBenthic
