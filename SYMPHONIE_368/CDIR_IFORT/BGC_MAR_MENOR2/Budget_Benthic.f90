










      SUBROUTINE BUDGET_BENTHIC


!_____________________________________________________________________*
!                                                                     *
! LAST REVISION: 8 JANUARY 2010                                       *
!                                                                     *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Calculates the budget term for the 1 compartiment             *
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
       MOYNITRIF_est,   &
       MOYDENITRIF_est, &
       MOYNO3Flux_est,  &
       MOYNH3Flux_est,  &
       MOYNMin_est,     &
       MOYCMin_est,     &
       MOYNDepo_est,    &
       MOYCDepo_est,    &
       MOYAnoxMin_est

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
      IF(mask_t(i,j,kmax+1).EQ.1.and.lon_t(i,j)*rad2deg>-5.6)then !debut test

!        if  (lon_t(i,j)*rad2deg<10.or.    &
!        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42.or. &
!        lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25.or. &
!        lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25)THEN
!

      if ((lon_t(i,j)*rad2deg<10).or.    &
        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42).or. &
        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25).or. &
        (lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25))then


!
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

! EST_MED
        else
      x2=dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)   !surface de la maille masquée des pts fantome de la parall

! Moyenne spatiale cumulée dans le temps sur toute la simulation
      SUM_NITRIF_est = SUM_NITRIF_est + NitrificationB(I,J) / SEC2DAY *x2
      SUM_DENITRIF_est = SUM_DENITRIF_est + DenitrificationB(I,J) / SEC2DAY *x2
      SUM_NO3Flux_est = SUM_NO3Flux_est + fluxbio_w(I,J,iNitrate,1) *x2
      SUM_NH3Flux_est = SUM_NH3Flux_est + fluxbio_w(I,J,iAmmonium,1) *x2
      SUM_NMin_est    = SUM_NMin_est    + NMin_out(I,J) / SEC2DAY *x2
      SUM_CMin_est    = SUM_CMin_est    + CMin_out(I,J) / SEC2DAY *x2
      SUM_NDepo_est    = SUM_NDepo_est   + NDepo(I,J) *x2
      SUM_CDepo_est    = SUM_CDepo_est   + CDepo(I,J) *x2
      SUM_AnoxMin_est  = SUM_AnoxMin_est + AnoxMin(I,J) *x2

! Moyenne spatiale, moyennee sur une periode de temps  
      NitrifSMean_est   = NitrifSMean_est + NitrificationB(I,J)/ SEC2DAY *x2
      DenitrifSMean_est = DenitrifSMean_est + DenitrificationB(I,J)/ SEC2DAY *x2 
      NO3FluxSMean_est  = NO3FluxSMean_est + fluxbio_w(I,J,iNitrate,1) *x2
      NH3FluxSMean_est  = NH3FluxSMean_est + fluxbio_w(I,J,iAmmonium,1) *x2
      NMinSMean_est    =  NMinSMean_est    + NMin_out(I,J) / SEC2DAY *x2
      CMinSMean_est    =  CMinSMean_est    + CMin_out(I,J) / SEC2DAY *x2
      NDepoSMean_est    = NDepoSMean_est   + NDepo(I,J) *x2
      CDepoSMean_est    = CDepoSMean_est   + CDepo(I,J) *x2
      AnoxMinSMean_est  = AnoxMinSMean_est + AnoxMin(I,J) *x2


! Carte flux cumules
      SUMT_NITRIF_est(I,J) = SUMT_NITRIF_est(I,J) + NitrificationB(I,J) *x2
      SUMT_DENITRIF_est(I,J)= SUMT_DENITRIF_est(I,J)+ DenitrificationB(I,J) *x2
      SUMT_NH3Flux_est(I,J) = SUMT_NH3Flux_est(I,J) + fluxbio_w(I,J,iAmmonium,1) *x2
      SUMT_NO3Flux_est(I,J) = SUMT_NO3Flux_est(I,J) + fluxbio_w(I,J,iNitrate,1) *x2
      SUMT_NMin_est(I,J) = SUMT_NMin_est(I,J) + NMin_out(I,J) *x2
      SUMT_CMin_est(I,J) = SUMT_CMin_est(I,J) + CMin_out(I,J) *x2
      SUMT_NDepo_est(I,J) = SUMT_NDepo_est(I,J) + NDepo(I,J) *x2
      SUMT_CDepo_est(I,J) = SUMT_CDepo_est(I,J) + CDepo(I,J) *x2
      SUMT_AnoxMin_est(I,J) = SUMT_AnoxMin_est(I,J) + AnoxMin(I,J) *x2




        endif ! delimite geographique
      ENDIF !fin test mask
      ENDDO !fin boucle i
      ENDDO !fin boucle j
      ENDIF !I1DBIO


! Sortie
       NMOD=12.

      IF(MOD((iteration3d),NMOD).EQ.0.and.iteration3d/=0) THEN

      MOYNITRIF         =      SUM_NITRIF  / sumareawmed  
         NitrifSMean    =     NitrifSMean  / sumareawmed/ TPS_BENT 

      MOYDENITRIF       =      SUM_DENITRIF  / sumareawmed
       DeNitrifSMean    =     DenitrifSMean  / sumareawmed/ TPS_BENT

      MOYNO3Flux         =      SUM_NO3Flux  / sumareawmed
         NO3FluxSMean    =     NO3FluxSMean  / sumareawmed/ TPS_BENT

      MOYNH3Flux         =      SUM_NH3Flux  / sumareawmed
         NH3FluxSMean    =     NH3FluxSMean  / sumareawmed/ TPS_BENT

      MOYNMin           =      SUM_NMin  / sumareawmed
         NMinSMean    =     NMinSMean  / sumareawmed/ TPS_BENT

      MOYCMin           =      SUM_CMin  / sumareawmed
         CMinSMean    =     CMinSMean  / sumareawmed/ TPS_BENT

      MOYNDepo           =      SUM_NDepo  / sumareawmed
         NDepoSMean    =     NDepoSMean  / sumareawmed/ TPS_BENT

      MOYCDepo           =      SUM_CDepo  / sumareawmed
         CDepoSMean    =     CDepoSMean  / sumareawmed/ TPS_BENT

      MOYAnoxMin           =      SUM_AnoxMin  / sumareawmed
         AnoxMinSMean    =     AnoxMinSMean  / sumareawmed/ TPS_BENT


! EST_MED

      MOYNITRIF_est         =      SUM_NITRIF_est  / sumareaemed
         NitrifSMean_est    =     NitrifSMean_est  / sumareaemed/ TPS_BENT

      MOYDENITRIF_est       =      SUM_DENITRIF_est  / sumareaemed
       DeNitrifSMean_est    =     DenitrifSMean_est  / sumareaemed/ TPS_BENT

      MOYNO3Flux_est         =      SUM_NO3Flux_est  / sumareaemed
         NO3FluxSMean_est    =     NO3FluxSMean_est  / sumareaemed/ TPS_BENT

      MOYNH3Flux_est         =      SUM_NH3Flux_est  / sumareaemed
         NH3FluxSMean_est    =     NH3FluxSMean_est  / sumareaemed/ TPS_BENT

      MOYNMin_est           =      SUM_NMin_est  / sumareaemed
         NMinSMean_est    =     NMinSMean_est  / sumareaemed/ TPS_BENT

      MOYCMin_est           =      SUM_CMin_est  / sumareaemed
         CMinSMean_est    =     CMinSMean_est  / sumareaemed/ TPS_BENT

      MOYNDepo_est           =      SUM_NDepo_est  / sumareaemed
         NDepoSMean_est    =     NDepoSMean_est  / sumareaemed/ TPS_BENT

      MOYCDepo_est           =      SUM_CDepo_est  / sumareaemed
         CDepoSMean_est    =     CDepoSMean_est  / sumareaemed/ TPS_BENT

      MOYAnoxMin_est           =      SUM_AnoxMin_est  / sumareaemed
         AnoxMinSMean_est    =     AnoxMinSMean_est  / sumareaemed/ TPS_BENT



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



! on archive
        if(par%rank==0)then                           !nnnnnnnnnnnnnnnn
        OPEN(UNIT=34                                  &
            ,FILE=DIRGRAPH(1:LNAME4)//'Benth_ouest.out'     &
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


! EST_MED
      call mpi_allreduce(MOYNITRIF_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNITRIF_est=sum1glb
      call mpi_allreduce(NitrifSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NitrifSMean_est=sum1glb
      call mpi_allreduce(MOYDENITRIF_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYDENITRIF_est=sum1glb
      call mpi_allreduce(DeNitrifSMean_est,sum1glb,1,mpi_double_precision, & !#MPI
                         mpi_sum,par%comm2d,ierr)
      DeNitrifSMean_est=sum1glb
      call mpi_allreduce(MOYNO3Flux_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNO3Flux_est=sum1glb
      call mpi_allreduce(NO3FluxSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NO3FluxSMean_est=sum1glb
      call mpi_allreduce(MOYNH3Flux_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNH3Flux_est=sum1glb
      call mpi_allreduce(NH3FluxSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NH3FluxSMean_est=sum1glb
      call mpi_allreduce(MOYNMin_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNMin_est=sum1glb
      call mpi_allreduce(NMinSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NMinSMean_est=sum1glb
      call mpi_allreduce(MOYCMin_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYCMin_est=sum1glb
      call mpi_allreduce(CMinSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      CMinSMean_est=sum1glb
      call mpi_allreduce(MOYNDepo_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYNDepo_est=sum1glb
      call mpi_allreduce(NDepoSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      NDepoSMean_est=sum1glb
      call mpi_allreduce(MOYCDepo_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYCDepo_est=sum1glb
      call mpi_allreduce(CDepoSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      CDepoSMean_est=sum1glb
      call mpi_allreduce(MOYAnoxMin_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      MOYAnoxMin_est=sum1glb
      call mpi_allreduce(AnoxMinSMean_est,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      AnoxMinSMean_est=sum1glb



! on archive
        if(par%rank==0)then                           !nnnnnnnnnnnnnnnn
        OPEN(UNIT=34                                  &
            ,FILE=DIRGRAPH(1:LNAME4)//'Benth_est.out'     &
           ,POSITION='APPEND')
!        WRITE(34,103)FLOAT(KOUNT)*DTI_FW/86400.   &
        WRITE(34,103)elapsedtime_now/86400.   &
                 ,MOYNITRIF_est            &
                 ,NitrifSMean_est          &
                 ,MOYDENITRIF_est          &
                 ,DeNitrifSMean_est        &
                 ,MOYNO3Flux_est           &
                 ,NO3FluxSMean_est         &
                 ,MOYNH3Flux_est           &
                 ,NH3FluxSMean_est         &
                 ,MOYNMin_est              &
                 ,NMinSMean_est            &
                 ,MOYCMin_est              &
                 ,CMinSMean_est            &
                 ,MOYNDepo_est             &
                 ,NDepoSMean_est           &
                 ,MOYCDepo_est             &
                 ,CDepoSMean_est           &
                 ,MOYAnoxMin_est           &
                 ,AnoxMinSMean_est

        CLOSE(34)
        endif                                         !nnnnnnnnnnnnnnnnn




 103   FORMAT(19(E12.3,' '))

! on remet à zero
         TPS_BENT = 0.
         NitrifSMean=0.
       DeNitrifSMean=0.
         NO3FluxSMean=0.
         NH3FluxSMean=0.
         NMinSMean=0.
         CMinSMean=0.
         NDepoSMean=0.
         CDepoSMean=0.
         AnoxMinSMean=0.

         NitrifSMean_est=0.
       DeNitrifSMean_est=0.
         NO3FluxSMean_est=0.
         NH3FluxSMean_est=0.
         NMinSMean_est=0.
         CMinSMean_est=0.
         NDepoSMean_est=0.
         CDepoSMean_est=0.
         AnoxMinSMean_est=0.


      ENDIF ! mod


 100  FORMAT(F7.2,1X,6(F11.4,1X))

      RETURN
      END
