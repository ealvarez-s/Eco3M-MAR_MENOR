










      SUBROUTINE SORTIE_SEDIMENTATION

      USE ModuleDeclaration
      USE MODULE_PRINCIPAL
      use module_parallele !#MPI
      use module_global
      IMPLICIT NONE

      IF(I1DBIO.EQ.0)THEN      !°°°°°°°°°°°°°°°°°°°°°°°  ECRITURE POUR 3D seulement
      
      LREC=(MBIO2-MBIO1+1)*(NBIO2-NBIO1+1)*4

      X1=86400./DTI_FW
      NC=MAX(1,INT(FLOAT(iteration3d-KOUNT0)/X1)) 
!     CALL KOUNT_TO_DATE(iteration3d)
      call elapsedtime2date(elapsedtime_now,year_now,month_now,day_now,hour_now,minute_now,second_now)
      if(par%rank==0)then
      OPEN(UNIT=3,FILE='date_Sed',POSITION='APPEND')
      WRITE(3,'(a25,1x,I6,1X,I4,1X,I4,1X,I2,1X,I2,1X,I2,1X,I2)')       &
       'ecriture Sed kount rec',iteration3d,NC,year_now,month_now,day_now,hour_now,minute_now
      CLOSE(3)
      endif


! A la fin ecriture de la sedimentation de MOP C
      DO I=MBIO1,MBIO2
      DO J=NBIO1,NBIO2
      ANYVAR2D(I,J)=SUM_EXPORTC_BOT_2D(I,J)           ! claude ici ajoute encore sum_
      ENDDO
      ENDDO
      OPEN(UNIT=3,FILE=                &
       DIRGRAPH(1:LNAME4)//'sedC'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,FORM=         &
       'UNFORMATTED')
      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
      CLOSE(3)
! Sedimentation de N
      DO I=MBIO1,MBIO2
      DO J=NBIO1,NBIO2
      ANYVAR2D(I,J)=SUM_EXPORTN_BOT_2D(I,J)
      ENDDO
      ENDDO
      OPEN(UNIT=3,FILE=    &
       DIRGRAPH(1:LNAME4)//'sedN'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,FORM=      &
       'UNFORMATTED')
      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
      CLOSE(3)
! Sedimentation de N
      DO I=MBIO1,MBIO2
      DO J=NBIO1,NBIO2
      ANYVAR2D(I,J)=SUM_EXPORTP_BOT_2D(I,J)
      ENDDO
      ENDDO
      OPEN(UNIT=3,FILE=      &
       DIRGRAPH(1:LNAME4)//'sedP'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,FORM=        &
       'UNFORMATTED')
      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
      CLOSE(3)
! Sedimentation de N
      DO I=MBIO1,MBIO2
      DO J=NBIO1,NBIO2
      ANYVAR2D(I,J)=SUM_EXPORTSI_BOT_2D(I,J)
      ENDDO
      ENDDO
      OPEN(UNIT=3,FILE=        &
      DIRGRAPH(1:LNAME4)//'sedSI'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,  &
       FORM='UNFORMATTED')
      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
      CLOSE(3)
!! Export à 200m de C
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=EXP2D(I,J,1)/TPS_STRADA_2D
!      ENDDO
!      ENDDO
!      OPEN(UNIT=3,FILE=        &
!       DIRGRAPH(1:LNAME4)//'ExportC'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,  &
!       FORM='UNFORMATTED')
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!! Export à 200m de C
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=EXP2D(I,J,2)/TPS_STRADA_2D
!      ENDDO
!      ENDDO
!      OPEN(UNIT=3,FILE=        &
!       DIRGRAPH(1:LNAME4)//'ExpsedC'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,  &
!       FORM='UNFORMATTED')
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!! Export à 200m de C
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=EXP2D(I,J,3)/TPS_STRADA_2D
!      ENDDO
!      ENDDO
!      OPEN(UNIT=3,FILE=        &
!       DIRGRAPH(1:LNAME4)//'ExpturC'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,  &
!       FORM='UNFORMATTED')
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!! Export à 200m de C
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=EXP2D(I,J,4)/TPS_STRADA_2D
!      ENDDO
!      ENDDO
!      OPEN(UNIT=3,FILE=        &
!       DIRGRAPH(1:LNAME4)//'ExpadvC'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,  &
!       FORM='UNFORMATTED')
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)
!! Export à 200m de C
!      DO I=MBIO1,MBIO2
!      DO J=NBIO1,NBIO2
!      ANYVAR2D(I,J)=EXP2D(I,J,5)/TPS_STRADA_2D
!      ENDDO
!      ENDDO
!      OPEN(UNIT=3,FILE=        &
!       DIRGRAPH(1:LNAME4)//'ExpmodC'//dom_c//'.binrec',ACCESS='DIRECT',RECL=LREC,  &
!       FORM='UNFORMATTED')
!      WRITE(3,REC=NC)ANYVAR2D(MBIO1:MBIO2,NBIO1:NBIO2)
!      CLOSE(3)

!      tps_strada_2d=0
!      do j=nbio1,nbio2 ! debut boucle sur j
!      do i=mbio1,mbio2 ! debut boucle sur i
!        exp2d(i,j,1)=0.
!        exp2d(i,j,2)=0.
!        exp2d(i,j,3)=0.
!        exp2d(i,j,4)=0.
!        exp2d(i,j,5)=0.
!        exp2d(i,j,6)=0.
!        exp2d(i,j,7)=0.
!      enddo    ! fin de boucle i1
!      enddo    ! fin de boucle j1
 
      ENDIF

 100  FORMAT(F7.2,1X,F11.4,1X,F13.5,1X,34(F11.4,1X))

      RETURN
      END

