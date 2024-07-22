










      SUBROUTINE BCFlux

!_____________________________________________________________________*
! 3d ecosystem model                                                  *
!                                                                     *
! LAST REVISION: 6 NOVEMBER 2009                                      *
!                                                                     *
! Implementation: Caroline Ulses                                      *
!                 LA                                                  *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Computes the flux of biogeochemical variables at                    *
! the Boundaries                                                      *
!                                                                     *
! 1. At the air-sea interface                                         * 
! 2. At the water-sediment interface                                  *
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
      use module_parallele !#MPI
      use module_global

      IMPLICIT NONE

! Local variables
      REAL LATOBS,LONOBS,SEC2DAY

      CHARACTER STATION(9)*9,NSTATION(14)*9

      DOUBLE PRECISION X15,X16,X17,X18,X19 

!---------------------------------------------------------------------*



!=====================================================================*
! 1. Flux at the air-sea interface                                    *
! BEGINNING                                                           *
!=====================================================================*

!      fluxbio_w(:,:,:,2)=0.       !25-12-2014

!=====================================================================*
! 1. Flux at the air-sea interface                                    *
! END                                                                 *
!=====================================================================*

!=====================================================================*
! 2. Flux at the water-sediment interface                             *
! BEGINNING                                                           *
!=====================================================================*

! 2.1 Compute the deposition rate at the bottom
! ----------------------------------------------

       DO J=1,jmax
       DO I=1,imax

       CDepo(I,J) = -(wsed(1,iLMOPC)*bio_t(I,J,kmin_w(i,j),iLMOPC)+    &
                      wsed(1,iSMOPC)*bio_t(I,J,kmin_w(i,j),iSMOPC)+     &
                      wsed(1,iDIAC) *bio_t(I,J,kmin_w(i,j),iDIAC )) 

       NDepo(I,J) = -(wsed(1,iLMOPN)*bio_t(I,J,kmin_w(i,j),iLMOPN)+    &
                      wsed(1,iSMOPN)*bio_t(I,J,kmin_w(i,j),iSMOPN)+     &
                      wsed(1,iDIAN) *bio_t(I,J,kmin_w(i,j),iDIAN ))

       PDepo(I,J) = -(wsed(1,iLMOPP)*bio_t(I,J,kmin_w(i,j),iLMOPP)+    &
                      wsed(1,iSMOPP)*bio_t(I,J,kmin_w(i,j),iSMOPP)+     &
                      wsed(1,iDIAP) *bio_t(I,J,kmin_w(i,j),iDIAP ))

       SiDepo(I,J) = -(wsed(1,iLMOPSI)*bio_t(I,J,kmin_w(i,j),iLMOPSI)+ &
                       wsed(1,iDIASI) *bio_t(I,J,kmin_w(i,j),iDIASI ))


       ENDDO
       ENDDO

! 2.2 Coupling with the 1 compartment
! -----------------------------------------

! No condition 
! Stock the deposition rate for initialisation of meta or full diagenetic model
       IF(IBenthic.EQ.0) THEN
!      CALL StockDepo

! Reflective condition 
! Stock the deposition rate for initialisation of meta or full diagenetic model
       ELSEIF(IBenthic.EQ.1) THEN
       CALL Reflec

! Simplified meta diagenetic model (pNit and pdeNit are constants)
       ELSEIF(IBenthic.EQ.2) THEN
!       print*,'sur le point d entrer dans SimpleDiaMeta'
       CALL SimpleDiaMeta

! Meta diagenetic model
       ELSEIF(IBenthic.EQ.3) THEN
       print*,'sur le point d entrer dans DiaMeta' 
       CALL DiaMeta

       ENDIF


!=====================================================================*
! 2. Flux at the water-sediment interface                             *
! END                                                                 *
!=====================================================================*


!       write(6,*)'Passe par InitBenthic' 

      IF(MOD((iteration3d-KOUNT0),72).EQ.0) THEN

      IF(IBENTHIC.GE.2) THEN

! Sorties aux stations CHACCRA

!      STATION(1)='Station A'
!      STATION(2)='Station B'
!      STATION(3)='Station N'
!      STATION(4)='Station C'
!      STATION(5)='Station F'
!      STATION(6)='Station J'
!      STATION(7)='Station I'
!      STATION(8)='Station K'
!      STATION(9)='Station L'

      NSTATION(1)='StatA.res'
      NSTATION(2)='StatC.res'
      NSTATION(3)='StatF.res'
      NSTATION(4)='StatJ.res'
      NSTATION(5)='Medoc.res'
      NSTATION(6)='Dyfam.res'
      NSTATION(7)='GdLPI.res'
      NSTATION(8)='GdLPM.res'
      NSTATION(9)='GdLP3.res'
      NSTATION(10)='StatB.res'
      NSTATION(11)='StatN.res'
      NSTATION(12)='StatK.res'
      NSTATION(13)='StatL.res'
      NSTATION(14)='StatI.res'


!     OPEN(UNIT=11,FILE='StationChaccra',ACCESS='APPEND')
!       write(11,*)'STATION(K1),DECI,DECJ,Prof,NitFlux,AmmoFlux,Nitrifica
!     & tion, denitrification'

      DO K1=1,14

      IF(K1.EQ.1) THEN
        LATOBS=43+18.5/60.
        LONOBS=4+51.1/60
      ELSEIF(K1.EQ.2) THEN
        LATOBS=43+16.433/60.
        LONOBS=4+46.544/60
      ELSEIF(K1.EQ.3) THEN
        LATOBS=43+9.925/60.
        LONOBS=4+39.038/60
      ELSEIF(K1.EQ.4) THEN
        LATOBS=43+15.993/60.
        LONOBS=4+58.05/60
      ELSEIF(K1.EQ.5) THEN
        LATOBS=42.
        LONOBS=5.
      ELSEIF(K1.EQ.6) THEN
        LATOBS=43.+25./60.
        LONOBS=7.+52/60.
      ELSEIF(K1.EQ.7) THEN
        LATOBS=43.076
        LONOBS=3.390448
      ELSEIF(K1.EQ.8) THEN
        LATOBS=42.9
        LONOBS=4.020855
      ELSEIF(K1.EQ.9) THEN
        LATOBS=43.36
        LONOBS=4.381089
      ELSEIF(K1.EQ.10) THEN
        LATOBS=43+18.2/60.
        LONOBS=4+49.765/60
      ELSEIF(K1.EQ.11) THEN
        LATOBS=43+17.566/60.
        LONOBS=4+47.726/60
      ELSEIF(K1.EQ.12) THEN
        LATOBS=43+17.91/60.
        LONOBS=4+51.345/60
      ELSEIF(K1.EQ.13) THEN
        LATOBS=43+18.191/60.
        LONOBS=4+53./60
      ELSEIF(K1.EQ.14) THEN
        LATOBS=43+15.81/60.
        LONOBS=4+52.916/60
      ENDIF

        LATIT1=LATOBS*PI/180.
        LONGI1=LONOBS*PI/180.
        CALL LATLON_TO_IJ('glb')     ! resultat dans le repere global ....

! position en indice decimal le long de l'axe:

       I1=INT(DECI)
       J1=INT(DECJ)
! passe dans le repere du proc
       i=i1-par%timax(1)
       j=j1-par%tjmax(1)
      if(i>=1.and.i<=imax-2.and.j>= 1.and.j<=jmax-2)then
      OPEN(UNIT=11,FILE=NSTATION(K1),POSITION='APPEND')

       RAPI=DECI    -FLOAT(I)

      X4=(1.-RAPI)*(1.-RAPJ)*depth_w(I  ,J  ,1) &
        +(1.-RAPI)*    RAPJ *depth_w(I  ,J+1,1) &
        +    RAPI *(1.-RAPJ)*depth_w(I+1,J  ,1) &
        +    RAPI *    RAPJ *depth_w(I+1,J+1,1)

! Conversion /d into /s
      SEC2DAY = 86400.

      X8=(1.-RAPI)*(1.-RAPJ)*fluxbio_w(I  ,J  ,iNitrate  ,1) &
        +(1.-RAPI)*    RAPJ *fluxbio_w(I  ,J+1,iNitrate  ,1) &
        +    RAPI *(1.-RAPJ)*fluxbio_w(I+1,J  ,iNitrate  ,1) &
        +    RAPI *    RAPJ *fluxbio_w(I+1,J+1,iNitrate  ,1)

      X8=X8*SEC2DAY

      X9=(1.-RAPI)*(1.-RAPJ)*fluxbio_w(I  ,J  ,iAmmonium ,1) &
        +(1.-RAPI)*    RAPJ *fluxbio_w(I  ,J+1,iAmmonium ,1) &
        +    RAPI *(1.-RAPJ)*fluxbio_w(I+1,J  ,iAmmonium ,1) &
        +    RAPI *    RAPJ *fluxbio_w(I+1,J+1,iAmmonium ,1)

      X9=X9*SEC2DAY

      X10=(1.-RAPI)*(1.-RAPJ)*NitrificationB(I  ,J ) &
        +(1.-RAPI)*    RAPJ *NitrificationB(I  ,J+1) &
        +    RAPI *(1.-RAPJ)*NitrificationB(I+1,J  ) &
        +    RAPI *    RAPJ *NitrificationB(I+1,J+1)

      X11=(1.-RAPI)*(1.-RAPJ)*DenitrificationB(I  ,J ) &
        +(1.-RAPI)*    RAPJ *DenitrificationB(I  ,J+1) &
        +    RAPI *(1.-RAPJ)*DenitrificationB(I+1,J  ) &
        +    RAPI *    RAPJ *DenitrificationB(I+1,J+1)


      X12 = (1.-RAPI)*(1.-RAPJ)*CDepo(I  ,J ) &
        +(1.-RAPI)*    RAPJ *CDepo(I  ,J+1) &
        +    RAPI *(1.-RAPJ)*CDepo(I+1,J  ) &
        +    RAPI *    RAPJ *CDepo(I+1,J+1)

      X12=X12*SEC2DAY

      X13 = (1.-RAPI)*(1.-RAPJ)*NDepo(I  ,J )  &
        +(1.-RAPI)*    RAPJ *NDepo(I  ,J+1)  &
        +    RAPI *(1.-RAPJ)*NDepo(I+1,J  )  &
        +    RAPI *    RAPJ *NDepo(I+1,J+1)

      X13=X13*SEC2DAY

      X14 = (1.-RAPI)*(1.-RAPJ)*CMin_out(I  ,J ) &
        +(1.-RAPI)*    RAPJ *CMin_out(I  ,J+1)  &
        +    RAPI *(1.-RAPJ)*CMin_out(I+1,J  )  &
        +    RAPI *    RAPJ *CMin_out(I+1,J+1)


      X15 = (1.-RAPI)*(1.-RAPJ)*NMin_out(I  ,J )  &
        +(1.-RAPI)*    RAPJ *NMin_out(I  ,J+1)    &
        +    RAPI *(1.-RAPJ)*NMin_out(I+1,J  )    &
        +    RAPI *    RAPJ *NMin_out(I+1,J+1) 


      X16 = (1.-RAPI)*(1.-RAPJ)*PPB2D(I  ,J ,1)  &
        +(1.-RAPI)*    RAPJ *PPB2D(I  ,J+1,1)    &
        +    RAPI *(1.-RAPJ)*PPB2D(I+1,J  ,1)    &
        +    RAPI *    RAPJ *PPB2D(I+1,J+1,1)

!	write(6,*)'TPS_PPB_2D =', TPS_PPB_2D ! Test Alex 01/09/17

!      X16=X16/TPS_PPB_2D
      
      if(TPS_PPB_2D>0) THEN
	X16=X16/TPS_PPB_2D
      ELSE
	X16=-9999.
      ENDIF

      X17 = (1.-RAPI)*(1.-RAPJ)*NPB2D(I  ,J ,1)  & 
        +(1.-RAPI)*    RAPJ *NPB2D(I  ,J+1,1)    &
        +    RAPI *(1.-RAPJ)*NPB2D(I+1,J  ,1)    &
        +    RAPI *    RAPJ *NPB2D(I+1,J+1,1)

      IF(TPS_PPB_2D.GT.0) X17=X17/TPS_PPB_2D

      X18 = (1.-RAPI)*(1.-RAPJ)*RPB2D(I  ,J ,1)  &
        +(1.-RAPI)*    RAPJ *RPB2D(I  ,J+1,1)    &
        +    RAPI *(1.-RAPJ)*RPB2D(I+1,J  ,1)    &
        +    RAPI *    RAPJ *RPB2D(I+1,J+1,1)

      IF(TPS_PPB_2D.GT.0) X18=X18/TPS_PPB_2D

      X19 = (1.-RAPI)*(1.-RAPJ)*AnoxMin(I  ,J )   &
        +(1.-RAPI)*    RAPJ *AnoxMin(I  ,J+1)     &
        +    RAPI *(1.-RAPJ)*AnoxMin(I+1,J  )     &
        +    RAPI *    RAPJ *AnoxMin(I+1,J+1)

!       IF(TPS_PPB_2D.GT.0) X19=X19/TPS_PPB_2D   

! On archive
      WRITE(11,100)FLOAT(iteration3d-KOUNT0)*DTI_FW/86400.,DECI,DECJ,  &
                 X4,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18    &
                ,X19 



  100 FORMAT (F7.2,1X,15(F11.4,1X))

      CLOSE(11)


         OPEN(UNIT=34,FILE=                                  &
              DIRGRAPH(1:LNAME4)//'bio3d'//NSTATION(K1)      &
              ,POSITION='APPEND')
         DO K=kmin_w(I,J),kmax
            WRITE(34,102)FLOAT(iteration3d-KOUNT0)*DTI_FW/86400.   &
                 ,depth_t(I,J,K)                            &
                 ,TEM_t(I,J,K,1)                             &
                 ,(BIO_t(I,J,K,vB),vB=1,vBMAX)                
         ENDDO
         CLOSE(34)

      endif

 102  FORMAT(F7.2,1X,F11.4,1X,F13.5,1X,33(F11.4,1X))

      ENDDO ! K1

      ENDIF ! ibenth

      ENDIF ! mod


! Verification
!     CALL GRAPH_OUT

!     STOP'dans BCFlux'

      END SUBROUTINE BCFlux
