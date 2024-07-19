      SUBROUTINE BUDGET_BENTHIC


!_____________________________________________________________________*
!                                                                     *
! LAST REVISION: 8 JANUARY 2010                                       *
!                                                                     *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Calculates the budget term for the benthic compartiment             *
!                                                                     *
!_____________________________________________________________________*

!_____________________________________________________________________*
!                                                                     *
! Modifications:                                                      *
!                                                                     *
!                                                                     *
!_____________________________________________________________________*



!---------------------------------------------------------------------*
! Declarations:


! Global variables

      USE MODULE_PRINCIPAL
      USE ModuleDeclaration
      use module_parallele !#MPI
      use module_global

      IMPLICIT NONE

! Local variables
      REAL SEC2DAY

      DOUBLE PRECISION   &
       MOYNITRIF,        &
       MOYDENITRIF,      &
       MOYNO3Flux,       &
       MOYNH3Flux,       &
       MOYNMin,          &
       MOYCMin,          &
       MOYNDepo,         &
       MOYCDepo,         &
       MOYAnoxMin,       &
       MOYNITRIF_plat,   &
       MOYDENITRIF_plat, &
       MOYNO3Flux_plat,  &
       MOYNH3Flux_plat,  &
       MOYNMin_plat,     &
       MOYCMin_plat,     &
       MOYNDepo_plat,    &
       MOYCDepo_plat,    &
       MOYAnoxMin_plat

      INTEGER NMOD
!---------------------------------------------------------------------*

      SEC2DAY = 86400.
      LREC=(MBIO2-MBIO1+1)*(NBIO2-NBIO1+1)*4

! calcul de l'export en dessous de 200 m
      IF(I1DBIO.NE.1)THEN
      TPS_BENT=TPS_BENT+1

      DO J=MAX0(NBIO1,1),NBIO2 ! debut boucle sur j
      DO I=MAX0(MBIO1,1),MBIO2 ! debut boucle sur i
!     print*,'i,j,k',i,j,kmax,nbio1,nbio2,mbio1,mbio2
      IF(mask_t(i,j,kmax+1).EQ.1.and.lon_t(i,j)*rad2deg>-5.6)THEN !debut test

      x2=dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)   !surface de la maille masquée des pts fantome de la parall

! Moyenne spatiale cumulée dans le temps sur toute la simulation
      SUM_NITRIF = SUM_NITRIF + NitrificationB(I,J) / SEC2DAY *x2
      SUM_DENITRIF = SUM_DENITRIF + DenitrificationB(I,J) / SEC2DAY *x2
      SUM_NO3Flux = SUM_NO3Flux + fluxbio_w(I,J,iNitrate,1) *x2
      SUM_NH3Flux = SUM_NH3Flux + fluxbio_w(I,J,iAmmonium,1) *x2
      SUM_NMin    = SUM_NMin    + NMin_out(I,J) / SEC2DAY *x2
      SUM_CMin    = SUM_CMin    + CMin_out(I,J) / SEC2DAY *x2
      SUM_NDepo    = SUM_NDepo   + NDepo(I,J) *x2
      SUM_CDepo    = SUM_CDepo   + CDepo(I,J) *x2
      SUM_AnoxMin  = SUM_AnoxMin + AnoxMin(I,J) *x2

! Moyenne spatiale, moyennee sur une periode de temps  
      NitrifSMean   = NitrifSMean + NitrificationB(I,J)/ SEC2DAY *x2
      DenitrifSMean = DenitrifSMean + DenitrificationB(I,J)/ SEC2DAY *x2 
      NO3FluxSMean  = NO3FluxSMean + fluxbio_w(I,J,iNitrate,1) *x2
      NH3FluxSMean  = NH3FluxSMean + fluxbio_w(I,J,iAmmonium,1) *x2
      NMinSMean    =  NMinSMean    + NMin_out(I,J) / SEC2DAY *x2
      CMinSMean    =  CMinSMean    + CMin_out(I,J) / SEC2DAY *x2
      NDepoSMean    = NDepoSMean   + NDepo(I,J) *x2
      CDepoSMean    = CDepoSMean   + CDepo(I,J) *x2
      AnoxMinSMean  = AnoxMinSMean + AnoxMin(I,J) *x2 

! Cas plateau GdL

       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)

      if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

       IF(-depth_w(I,J,kmin_w(I,J)).LE.200) THEN

! Moyenne spatiale cumulée dans le temps sur toute la simulation
      SUM_NITRIF_plat = SUM_NITRIF_plat + NitrificationB(I,J)/ SEC2DAY *x2
      SUM_DENITRIF_plat=SUM_DENITRIF_plat+DenitrificationB(I,J)/ SEC2DAY *x2
      SUM_NO3Flux_plat = SUM_NO3Flux_plat + fluxbio_w(I,J,iNitrate,1) *x2
      SUM_NH3Flux_plat = SUM_NH3Flux_plat + fluxbio_w(I,J,iAmmonium,1) *x2
      SUM_NMin_plat    = SUM_NMin_plat    + NMin_out(I,J) / SEC2DAY *x2
      SUM_CMin_plat    = SUM_CMin_plat    + CMin_out(I,J) / SEC2DAY *x2
      SUM_NDepo_plat    = SUM_NDepo_plat  + NDepo(I,J) *x2
      SUM_CDepo_plat    = SUM_CDepo_plat  + CDepo(I,J) *x2
      SUM_AnoxMin_plat  = SUM_AnoxMin_plat + AnoxMin(I,J) *x2

! Moyenne spatiale, moyennee sur une periode de temps
      NitrifSMean_plat=NitrifSMean_plat + NitrificationB(I,J)/ SEC2DAY *x2
      DenitrifSMean_plat=DenitrifSMean_plat       &
             +DenitrificationB(I,J)/SEC2DAY *x2
      NO3FluxSMean_plat  = NO3FluxSMean_plat + fluxbio_w(I,J,iNitrate,1) *x2
      NH3FluxSMean_plat  = NH3FluxSMean_plat+ fluxbio_w(I,J,iAmmonium,1) *x2
      NMinSMean_plat    =  NMinSMean_plat    + NMin_out(I,J) / SEC2DAY *x2
      CMinSMean_plat    =  CMinSMean_plat    + CMin_out(I,J) / SEC2DAY *x2
      NDepoSMean_plat    = NDepoSMean_plat   + NDepo(I,J) *x2
      CDepoSMean_plat    = CDepoSMean_plat   + CDepo(I,J) *x2
      AnoxMinSMean_plat  = AnoxMinSMean_plat + AnoxMin(I,J) *x2

      ENDIF     !depth_w 
      endif                                                             !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


! Carte flux cumules
      SUMT_NITRIF(I,J) = SUMT_NITRIF(I,J) + NitrificationB(I,J) *x2
      SUMT_DENITRIF(I,J)= SUMT_DENITRIF(I,J)+ DenitrificationB(I,J) *x2 
      SUMT_NH3Flux(I,J) = SUMT_NH3Flux(I,J) + fluxbio_w(I,J,iAmmonium,1) *x2
      SUMT_NO3Flux(I,J) = SUMT_NO3Flux(I,J) + fluxbio_w(I,J,iNitrate,1) *x2
      SUMT_NMin(I,J) = SUMT_NMin(I,J) + NMin_out(I,J) *x2
      SUMT_CMin(I,J) = SUMT_CMin(I,J) + CMin_out(I,J) *x2
      SUMT_NDepo(I,J) = SUMT_NDepo(I,J) + NDepo(I,J) *x2
      SUMT_CDepo(I,J) = SUMT_CDepo(I,J) + CDepo(I,J) *x2
      SUMT_AnoxMin(I,J) = SUMT_AnoxMin(I,J) + AnoxMin(I,J) *x2


      ENDIF !fin test mask
      ENDDO !fin boucle i
      ENDDO !fin boucle j
      ENDIF !I1DBIO


! Sortie
       NMOD=12.

      IF(MOD((iteration3d),NMOD).EQ.0.and.iteration3d/=0) THEN

      MOYNITRIF         =      SUM_NITRIF  / SUMAREATOTAL  
         NitrifSMean    =     NitrifSMean  / SUMAREATOTAL/ TPS_BENT 

      MOYDENITRIF       =      SUM_DENITRIF  / SUMAREATOTAL
       DeNitrifSMean    =     DenitrifSMean  / SUMAREATOTAL/ TPS_BENT

      MOYNO3Flux         =      SUM_NO3Flux  / SUMAREATOTAL
         NO3FluxSMean    =     NO3FluxSMean  / SUMAREATOTAL/ TPS_BENT

      MOYNH3Flux         =      SUM_NH3Flux  / SUMAREATOTAL
         NH3FluxSMean    =     NH3FluxSMean  / SUMAREATOTAL/ TPS_BENT

      MOYNMin           =      SUM_NMin  / SUMAREATOTAL
         NMinSMean    =     NMinSMean  / SUMAREATOTAL/ TPS_BENT

      MOYCMin           =      SUM_CMin  / SUMAREATOTAL
         CMinSMean    =     CMinSMean  / SUMAREATOTAL/ TPS_BENT

      MOYNDepo           =      SUM_NDepo  / SUMAREATOTAL
         NDepoSMean    =     NDepoSMean  / SUMAREATOTAL/ TPS_BENT

      MOYCDepo           =      SUM_CDepo  / SUMAREATOTAL
         CDepoSMean    =     CDepoSMean  / SUMAREATOTAL/ TPS_BENT

      MOYAnoxMin           =      SUM_AnoxMin  / SUMAREATOTAL
         AnoxMinSMean    =     AnoxMinSMean  / SUMAREATOTAL/ TPS_BENT

      MOYNITRIF_plat     =      SUM_NITRIF_plat  / SUMAREAGDL_plat
        NitrifSMean_plat=  NitrifSMean_plat  / SUMAREAGDL_plat/ TPS_BENT

      MOYDENITRIF_plat   =   SUM_DENITRIF_plat  / SUMAREAGDL_plat
      DeNitrifSMean_plat =DenitrifSMean_plat  /SUMAREAGDL_plat/ TPS_BENT

      MOYNO3Flux_plat     =      SUM_NO3Flux_plat  / SUMAREAGDL_plat
      NO3FluxSMean_plat  =NO3FluxSMean_plat  / SUMAREAGDL_plat/ TPS_BENT

      MOYNH3Flux_plat    =    SUM_NH3Flux_plat  / SUMAREAGDL_plat
      NH3FluxSMean_plat =NH3FluxSMean_plat  / SUMAREAGDL_plat/ TPS_BENT


      MOYNMin_plat           =      SUM_NMin_plat  / SUMAREAGDL_plat
         NMinSMean_plat =    NMinSMean_plat  / SUMAREAGDL_plat/ TPS_BENT

      MOYCMin_plat        =      SUM_CMin_plat  / SUMAREAGDL_plat
         CMinSMean_plat  =   CMinSMean_plat  / SUMAREAGDL_plat/ TPS_BENT

      MOYNDepo_plat           =      SUM_NDepo_plat  / SUMAREAGDL_plat
         NDepoSMean_plat = NDepoSMean_plat  / SUMAREAGDL_plat / TPS_BENT

      MOYCDepo_plat           =      SUM_CDepo_plat  / SUMAREAGDL_plat
         CDepoSMean_plat =  CDepoSMean_plat  / SUMAREAGDL_plat/ TPS_BENT

      MOYAnoxMin_plat         =      SUM_AnoxMin_plat  / SUMAREAGDL_plat
      AnoxMinSMean_plat = AnoxMinSMean_plat  / SUMAREAGDL_plat/ TPS_BENT

 
#ifdef parallele
      call mpi_allreduce(MOYNITRIF,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNITRIF=sum1glb
      call mpi_allreduce(NitrifSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NitrifSMean=sum1glb
      call mpi_allreduce(MOYDENITRIF,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYDENITRIF=sum1glb
      call mpi_allreduce(DeNitrifSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      DeNitrifSMean=sum1glb
      call mpi_allreduce(MOYNO3Flux,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNO3Flux=sum1glb
      call mpi_allreduce(NO3FluxSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NO3FluxSMean=sum1glb
      call mpi_allreduce(MOYNH3Flux,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNH3Flux=sum1glb
      call mpi_allreduce(NH3FluxSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NH3FluxSMean=sum1glb
      call mpi_allreduce(MOYNMin,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNMin=sum1glb
      call mpi_allreduce(NMinSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NMinSMean=sum1glb
      call mpi_allreduce(MOYCMin,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYCMin=sum1glb
      call mpi_allreduce(CMinSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      CMinSMean=sum1glb
      call mpi_allreduce(MOYNDepo,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNDepo=sum1glb
      call mpi_allreduce(NDepoSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NDepoSMean=sum1glb
      call mpi_allreduce(MOYCDepo,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYCDepo=sum1glb
      call mpi_allreduce(CDepoSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      CDepoSMean=sum1glb
      call mpi_allreduce(MOYAnoxMin,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYAnoxMin=sum1glb
      call mpi_allreduce(AnoxMinSMean,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      AnoxMinSMean=sum1glb
#endif


! on archive
        if(par%rank==0)then                           !nnnnnnnnnnnnnnnn
        OPEN(UNIT=34                                  &
            ,FILE=DIRGRAPH(1:LNAME4)//'Benth.out'     &
           ,POSITION='APPEND')
!        WRITE(34,103)FLOAT(KOUNT)*DTI_FW/86400.   &
        WRITE(34,103)elapsedtime_now/86400.   &
                 ,MOYNITRIF            &
                 ,NitrifSMean          &
                 ,MOYDENITRIF          &
                 ,DeNitrifSMean        &
                 ,MOYNO3Flux           &
                 ,NO3FluxSMean         &
                 ,MOYNH3Flux           &
                 ,NH3FluxSMean         &
                 ,MOYNMin              &
                 ,NMinSMean            &
                 ,MOYCMin              &
                 ,CMinSMean            &
                 ,MOYNDepo             &
                 ,NDepoSMean           &
                 ,MOYCDepo             &
                 ,CDepoSMean           &
                 ,MOYAnoxMin           &
                 ,AnoxMinSMean

        CLOSE(34)
        endif                                         !nnnnnnnnnnnnnnnnn

#ifdef parallele
      call mpi_allreduce(MOYNITRIF_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNITRIF_plat=sum1glb
      call mpi_allreduce(NitrifSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NitrifSMean_plat=sum1glb
      call mpi_allreduce(MOYDENITRIF_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYDENITRIF_plat=sum1glb
      call mpi_allreduce(DeNitrifSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      DeNitrifSMean_plat=sum1glb
      call mpi_allreduce(MOYNO3Flux_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNO3Flux_plat=sum1glb
      call mpi_allreduce(NO3FluxSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NO3FluxSMean_plat=sum1glb
      call mpi_allreduce(MOYNH3Flux_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNH3Flux_plat=sum1glb
      call mpi_allreduce(NH3FluxSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NH3FluxSMean_plat=sum1glb
      call mpi_allreduce(MOYNMin_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNMin_plat=sum1glb
      call mpi_allreduce(NMinSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NMinSMean_plat=sum1glb
      call mpi_allreduce(MOYCMin_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYCMin_plat=sum1glb
      call mpi_allreduce(CMinSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      CMinSMean_plat=sum1glb
      call mpi_allreduce(MOYNDepo_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNDepo_plat=sum1glb
      call mpi_allreduce(NDepoSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NDepoSMean_plat=sum1glb
      call mpi_allreduce(MOYCDepo_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYCDepo_plat=sum1glb
      call mpi_allreduce(CDepoSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      CDepoSMean_plat=sum1glb
      call mpi_allreduce(MOYAnoxMin_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYAnoxMin_plat=sum1glb
      call mpi_allreduce(AnoxMinSMean_plat,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      AnoxMinSMean_plat=sum1glb
#endif

!        if(par%rank==0)then                           !nnnnnnnnnnnnnnnn
!        OPEN(UNIT=34                    &
!            ,FILE=DIRGRAPH(1:LNAME4)//'Benth_plat.res'  &
!           ,POSITION='APPEND')
!        WRITE(34,103)elapsedtime_now/86400.  &
!                 ,MOYNITRIF_plat       &
!                 ,NitrifSMean_plat     &
!                 ,MOYDENITRIF_plat     &
!                 ,DeNitrifSMean_plat   &
!                 ,MOYNO3Flux_plat      &
!                 ,NO3FluxSMean_plat    &
!                 ,MOYNH3Flux_plat      &
!                 ,NH3FluxSMean_plat    &
!                 ,MOYNMin_plat         &
!                 ,NMinSMean_plat       &
!                 ,MOYCMin_plat         &
!                 ,CMinSMean_plat       &
!                 ,MOYNDepo_plat        &
!                 ,NDepoSMean_plat      &
!                 ,MOYCDepo_plat        &
!                 ,CDepoSMean_plat      &
!                 ,MOYAnoxMin_plat      &
!                 ,AnoxMinSMean_plat
!
!        CLOSE(34)
!        endif                                       !nnnnnnnnnnnnnnnnnnn


 103   FORMAT(19(E12.3,' '))

! on remet à zero
         TPS_BENT = 0
         NitrifSMean=0.
       DeNitrifSMean=0.
         NO3FluxSMean=0.
         NH3FluxSMean=0.
         NMinSMean=0.
         CMinSMean=0.
         NDepoSMean=0.
         CDepoSMean=0.
         AnoxMinSMean=0.
          NitrifSMean_plat=0.
       DeNitrifSMean_plat=0.
         NO3FluxSMean_plat=0.
         NH3FluxSMean_plat=0.
         NMinSMean_plat=0.
         CMinSMean_plat=0.
         NDepoSMean_plat=0.
         CDepoSMean_plat=0.
         AnoxMinSMean_plat=0.

      ENDIF ! mod


!      IF(MOD((KOUNT),360).EQ.0.) THEN
!   
!      X1=86400./DTI_FW
!      NC=MAX(1,INT(FLOAT(KOUNT)/X1))
!      CALL KOUNT_TO_DATE(KOUNT)
!      if(par%rank==0)then
!      OPEN(UNIT=3,FILE='date_Bent',POSITION='APPEND')
!      WRITE(3,'(a25,1x,I6,1X,I4,1X,I4,1X,I2,1X,I2,1X,I2,1X,I2)')  &
!       'ecriture Benth kount rec',KOUNT,NC,I5,I6,I7,I3,I2
!      CLOSE(3)
!      endif
!
!      OPEN(UNIT=3,FILE=                       &
!       DIRGRAPH(1:LNAME4)//'NitrifiB'//dom_c//'.binrec'  &
!       ,ACCESS='DIRECT',RECL=LREC             &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_NITRIF(I,J)
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=       &
!       DIRGRAPH(1:LNAME4)//'DenitB'//dom_c//'.binrec' &
!       ,ACCESS='DIRECT',RECL=LREC          &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_DENITRIF(I,J) 
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                     &
!       DIRGRAPH(1:LNAME4)//'NO3Flux'//dom_c//'.binrec' &
!       ,ACCESS='DIRECT',RECL=LREC           &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_NO3Flux(I,J) *SEC2DAY
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                      &
!       DIRGRAPH(1:LNAME4)//'NH3Flux'//dom_c//'.binrec'  &
!       ,ACCESS='DIRECT',RECL=LREC            &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_NH3Flux(I,J) *SEC2DAY
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                     &
!       DIRGRAPH(1:LNAME4)//'NMin'//dom_c//'.binrec'    &
!       ,ACCESS='DIRECT',RECL=LREC           &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_NMin(I,J) 
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                &
!       DIRGRAPH(1:LNAME4)//'CMin'//dom_c//'.binrec'  &
!       ,ACCESS='DIRECT',RECL=LREC         &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_CMin(I,J) 
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                      &
!       DIRGRAPH(1:LNAME4)//'NDepo'//dom_c//'.binrec'    &
!       ,ACCESS='DIRECT',RECL=LREC            &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_NDepo(I,J) 
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                      &
!       DIRGRAPH(1:LNAME4)//'CDepo'//dom_c//'.binrec'    &
!       ,ACCESS='DIRECT',RECL=LREC            &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=SUMT_CDepo(I,J) 
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!      OPEN(UNIT=3,FILE=                     &
!       DIRGRAPH(1:LNAME4)//'pAnOxic'//dom_c//'.binrec' &
!       ,ACCESS='DIRECT',RECL=LREC           &
!       ,FORM='UNFORMATTED')
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      IF(SUMT_CMin(I,J).GT.0) THEN
!      ANYVAR2D(I,J)=SUMT_AnoxMin(I,J) / SUMT_CMin(I,J)
!      ELSE
!      ANYVAR2D(I,J)=0.
!      ENDIF 
!      ENDDO
!      ENDDO
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!
!
!
!      ENDIF   !mod





 100  FORMAT(F7.2,1X,6(F11.4,1X))

      RETURN
      END
