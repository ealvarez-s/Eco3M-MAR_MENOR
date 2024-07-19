      SUBROUTINE BUDGET_EXPORT_BOTTOM


!_____________________________________________________________________*
!                                                                     *
! LAST REVISION: 27-12-2014                                           *
!                                                                     *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Calculates the export towards the bottom at 200m depth              *
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
      use ModuleDeclaration
      use module_principal
      use module_parallele !#MPI
      use module_global
      IMPLICIT NONE

! Local variables
      REAL VIT,EXP_SED,EXP_ADV,EXP_TURB,EXP_SED_BOT

!---------------------------------------------------------------------*

! calcul de l'export en dessous de 200 m
      if(I1DBIO/=1)then
      tps_strada=tps_strada+1
      tps_strada_2d=tps_strada_2d+1



!       do i=mbio1,mbio2
!       do j=nbio1,nbio2
       
!       Test Alex 12/10/2016
       do i=1,imax
       do j=1,jmax
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if (j2>=j1domain.and.j2<=j2domain.and.i2>=i1domain.and.i2<=i2domain) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! ATTENTION TEST PLUS BAS INCLU DANS CETTE ZONE
      if(mask_t(I,J,kmax+1)==1)then !debut test

! calcul du niveau juste au dessous de 200m, ou fond
      k=kmax
     do while(k>kmin_w(i,j).and.depth_t(i,j,k)>-200) ! pour diag12
        k=k-1
      enddo

! surface de la maille
      x2=dxdy_t(i,j)*mask_i_w(i)*mask_j_w(j)


!      k=30


      DO vb=1,vbmax      ! debut boucle vb


!      print*,'export_bottom',i,j,k,vb
      if(dabs(bio_t(i,j,k,vb))>1.d-30)then    !ce 09-12-2014
      EXP_SED=-wsed(1,vb)*bio_t(i,j,k,vb)
      else
      EXP_SED=0.
      endif


! export advectif
      EXP_ADV=0.
!     if(abs(bio_t(i,j,k,vb)+bio_t(i,j,k-1,vb))>1.e-10.and.abs(omega_w(i,j,k,1))>1.e-10)then
!      if(abs(bio_t(i,j,k,vb))<1.e4.and.abs(omega_w(i,j,k,1))<1.e4)then
     if(dabs(bio_t(i,j,k,vb)+bio_t(i,j,k-1,vb))>1.d-30)then
! Test Alex Floating overflow 11/01/18
!     write(6,*)'i =',i
!     write(6,*)'j =',j
!     write(6,*)'k =',k
!     write(6,*)'k-1 =',k-1
!     write(6,*)'vb =',vb
!     write(6,*)'bio_t k =', -bio_t(i,j,k,vb)
!     write(6,*)'bio_t k-1 =', bio_t(i,j,k-1,vb)
!     write(6,*)'omega_w =', omega_w(i,j,k,1)
      EXP_ADV=-0.5*(bio_t(i,j,k,vb)*mask_t(i,j,k)    &
                   +bio_t(i,j,k-1,vb)*mask_t(i,j,k-1)) &
                   *omega_w(i,j,k,1)
!     write(6,*)'exp_adv =', -0.5*(bio_t(i,j,k,vb)+bio_t(i,j,k-1,vb))*omega_w(i,j,k,1)
      else
!      write(6,*)'bio_t(i,j,k,vb) =', bio_t(i,j,k,vb)
!      write(6,*)'omega_w(i,j,k,1)', omega_w(i,j,k,1)
      EXP_ADV=0.
      endif


! export turbulent
      EXP_TURB=0.
!      if(k/=kmin_w(i,j).and.abs(bio_t(i,j,k,vb)-bio_t(i,j,k-1,vb))>1.e-30.and.mask_t(I,J,k-1)==1)then
      if(k/=kmin_w(i,j).and.dabs(bio_t(i,j,k,vb)-bio_t(i,j,k-1,vb))>1.d-30)then !ce 09-12-2014
      EXP_TURB=                                        &
             (bio_t(i,j,k  ,vb)                        &
             -bio_t(i,j,k-1,vb))                       &
            /(depth_t(i,j,k)-depth_t(i,j,k-1))         &
             *kh_w(i,j,k)                              
      else
      EXP_TURB=0.
      endif


! Modif by Alex 12/01/2018 => eliminating the details
! import de Nitrate de mod en mgC m-2 j-1
      if(vb==iNitrate) then
         EXP2D(I,J,1)= EXP2D(I,J,1)+(EXP_ADV+EXP_TURB)*86400
      endif

! export de lammonium en mmolN m-2 j-1
       if(vb==iAMMONIUM) then
         EXP2D(I,J,2)= EXP2D(I,J,2)+(EXP_ADV+EXP_TURB)*86400
       endif

! export de mod en mmolN m-2 j-1
       if(vb==iMODN) then
         EXP2D(I,J,3)= EXP2D(I,J,3)+(EXP_ADV+EXP_TURB)*86400
       endif

! export de PON en mmolN m-2 j-1
       if(vb==iDIAN.or.vb==iNANON.or.vb==iSYNEN.or. &
          vb==iLMOPN.or.vb==iSMOPN) then
         EXP2D(I,J,4)= EXP2D(I,J,4)+(EXP_SED+EXP_ADV+EXP_TURB)*86400
       endif

! export de ZOON en mmolN m-2 j-1
      if(vb==iZOONANOC.or.vb==iZOOMESOC.or.vb==iZOOMICROC) then
         EXP2D(I,J,4)= EXP2D(I,J,4)+(EXP_ADV+EXP_TURB)*0.18*86400 
       endif

! export de BACTN en mmolN m-2 j-1
       if(vb==iBACTC) then
         EXP2D(I,J,4)= EXP2D(I,J,4)+(EXP_ADV+EXP_TURB)*0.232*86400 ! correction de (EXP_ADV)*0.232*86400 23-3-2018
       endif

!!!!!!!!!!

!!! P BUDGET !!!
! import de Phosphate de mod en mgC m-2 j-1
      if(vb==iPhosphate) then
         EXP2D(I,J,5)= EXP2D(I,J,5)+(EXP_ADV+EXP_TURB)*86400
      endif

! export de mod en mmolP m-2 j-1
       if(vb==iMODP) then
         EXP2D(I,J,6)= EXP2D(I,J,6)+(EXP_ADV+EXP_TURB)*86400
       endif

! export de POP en mmolP m-2 j-1
       if(vb==iDIAP.or.vb==iNANOP.or.vb==iSYNEP.or. &
          vb==iLMOPP.or.vb==iSMOPP) then
         EXP2D(I,J,7)= EXP2D(I,J,7)+(EXP_SED+EXP_ADV+EXP_TURB)*86400
       endif

! export de ZOOP en mmolP m-2 j-1
      if(vb==iZOONANOC.or.vb==iZOOMESOC.or.vb==iZOOMICROC) then
         EXP2D(I,J,7)= EXP2D(I,J,7)+(EXP_ADV+EXP_TURB)*0.013*86400
       endif

! export de BACTP en mmolP m-2 j-1
       if(vb==iBACTC) then
         EXP2D(I,J,7)= EXP2D(I,J,7)+(EXP_ADV+EXP_TURB)*0.022*86400 ! correction de (EXP_ADV)*0.22*86400 23-3-2018
       endif

!!!!!!!!!!

!!! SI BUDGET !!!
! import de Silicate en mgC m-2 j-1
      if(vb==iSilice) then
         EXP2D(I,J,8)= EXP2D(I,J,8)+(EXP_ADV+EXP_TURB)*86400
      endif             
      
! export de POSi en mmolSi m-2 j-1
       if(vb==iDIASi.or. &
          vb==iLMOPSi.or.vb==iSMOPSi) then
         EXP2D(I,J,9)= EXP2D(I,J,9)+(EXP_SED+EXP_ADV+EXP_TURB)*86400
       endif

!!! C BUDGET
! export de POC
       if(vb==iBACTC.or.vb==iZOONANOC.or.vb==iZOOMESOC.or.vb==iZOOMICROC.or. &
          vb==iDIAC.or.vb==iNANOC.or.vb==iSYNEC.or. &
          vb==iLMOPC.or.vb==iSMOPC) then
         EXP2D(I,J,19)= EXP2D(I,J,19)+(EXP_ADV+EXP_TURB+EXP_SED)*86400
       endif

! export de DOC
       if(vb==iMODC) then
         EXP2D(I,J,20)= EXP2D(I,J,20)+(EXP_ADV+EXP_TURB)*86400
       endif


      enddo !fin boucle vb


! What does that mean...?
!    k=max(23,kmax)


      DO vb=1,vbmax      ! debut boucle vb


! export par sedimentation
      if(dabs(bio_t(i,j,k,vb))>1.d-30)then    !ce 09-12-2014
      EXP_SED=-wsed(1,vb)*bio_t(i,j,k,vb)
      else
      EXP_SED=0.
      endif

! export advectif
      EXP_ADV=0.
      if(dabs(bio_t(i,j,k,vb)+bio_t(i,j,k-1,vb))>1.d-30)then
      EXP_ADV=-0.5*(bio_t(i,j,k,vb)*mask_t(i,j,k)    &
                   +bio_t(i,j,k-1,vb)*mask_t(i,j,k-1)) &
                   *omega_w(i,j,k,1)
      else
      EXP_ADV=0.
      endif
         
      
! export turbulent
      EXP_TURB=0.
!      if(k/=kmin_w(i,j).and.abs(bio_t(i,j,k,vb)-bio_t(i,j,k-1,vb))>1.e-30.and.mask_t(I,J,k-1)==1)then
      if(k/=kmin_w(i,j).and.dabs(bio_t(i,j,k,vb)-bio_t(i,j,k-1,vb))>1.d-30)then !ce 09-12-2014
      EXP_TURB=                                        &
             (bio_t(i,j,k  ,vb)                        &
             -bio_t(i,j,k-1,vb))                       &
            /(depth_t(i,j,k)-depth_t(i,j,k-1))         &
             *kh_w(i,j,k)                              
      else
      EXP_TURB=0.
      endif

! Modif by Alex 12/01/2018 => eliminating the details

!!! N BUDGET !!!
! import de Nitrate de mod en mgC m-2 j-1
      if(vb==iNitrate) then
         EXP2D(I,J,10)= EXP2D(I,J,10)+(EXP_ADV+EXP_TURB)*86400
      endif

! export de lammonium en mmolN m-2 j-1
       if(vb==iAMMONIUM) then
         EXP2D(I,J,11)= EXP2D(I,J,11)+(EXP_ADV+EXP_TURB)*86400 
       endif

! export de mod en mmolN m-2 j-1
       if(vb==iMODN) then
         EXP2D(I,J,12)= EXP2D(I,J,12)+(EXP_ADV+EXP_TURB)*86400
       endif

! export de PON en mmolN m-2 j-1
       if(vb==iDIAN.or.vb==iNANON.or.vb==iSYNEN.or. &
          vb==iLMOPN.or.vb==iSMOPN) then
         EXP2D(I,J,13)= EXP2D(I,J,13)+(EXP_SED+EXP_ADV+EXP_TURB)*86400
       endif

! export de ZOON en mmolN m-2 j-1
      if(vb==iZOONANOC.or.vb==iZOOMESOC.or.vb==iZOOMICROC) then
         EXP2D(I,J,13)= EXP2D(I,J,13)+(EXP_ADV+EXP_TURB)*0.18*86400 
       endif

! export de BACTN en mmolN m-2 j-1
       if(vb==iBACTC) then
         EXP2D(I,J,13)= EXP2D(I,J,13)+(EXP_ADV+EXP_TURB)*0.232*86400 ! correction de (EXP_ADV)*0.232*86400 23-3-2018
       endif

!!!!!!!!!!

!!! P BUDGET !!!
! import de Phosphate de mod en mgC m-2 j-1
      if(vb==iPhosphate) then
         EXP2D(I,J,14)= EXP2D(I,J,14)+(EXP_ADV+EXP_TURB)*86400
      endif             
      
! export de mod en mmolP m-2 j-1
       if(vb==iMODP) then
         EXP2D(I,J,15)= EXP2D(I,J,15)+(EXP_ADV+EXP_TURB)*86400
       endif
      
! export de POP en mmolP m-2 j-1
       if(vb==iDIAP.or.vb==iNANOP.or.vb==iSYNEP.or. &
          vb==iLMOPP.or.vb==iSMOPP) then
         EXP2D(I,J,16)= EXP2D(I,J,16)+(EXP_SED+EXP_ADV+EXP_TURB)*86400
       endif

! export de ZOOP en mmolP m-2 j-1
      if(vb==iZOONANOC.or.vb==iZOOMESOC.or.vb==iZOOMICROC) then
         EXP2D(I,J,16)= EXP2D(I,J,16)+(EXP_ADV+EXP_TURB)*0.013*86400 
       endif

! export de BACTP en mmolP m-2 j-1
       if(vb==iBACTC) then
         EXP2D(I,J,16)= EXP2D(I,J,16)+(EXP_ADV+EXP_TURB)*0.022*86400 ! correction de (EXP_ADV)*0.22*86400 23-3-2018
       endif

!!!!!!!!!!

!!! SI BUDGET !!!
! import de Silicate en mgC m-2 j-1
      if(vb==iSilice) then
         EXP2D(I,J,17)= EXP2D(I,J,17)+(EXP_ADV+EXP_TURB)*86400
      endif

! export de POSi en mmolSi m-2 j-1
       if(vb==iDIASi.or. &
          vb==iLMOPSi.or.vb==iSMOPSi) then
         EXP2D(I,J,18)= EXP2D(I,J,18)+(EXP_SED+EXP_ADV+EXP_TURB)*86400
       endif

!!! C BUDGET
! export de POC
       if(vb==iBACTC.or.vb==iZOONANOC.or.vb==iZOOMESOC.or.vb==iZOOMICROC.or. &
          vb==iDIAC.or.vb==iNANOC.or.vb==iSYNEC.or. &
          vb==iLMOPC.or.vb==iSMOPC) then
         EXP2D(I,J,21)= EXP2D(I,J,21)+(EXP_ADV+EXP_TURB+EXP_SED)*86400
       endif

! export de DOC
       if(vb==iMODC) then
         EXP2D(I,J,22)= EXP2D(I,J,22)+(EXP_ADV+EXP_TURB)*86400
       endif


      enddo !fin boucle vb

      endif !fin test
      endif             !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      enddo !fin boucle i
      enddo !fin boucle j
      endif


 34      format(1(F13.5,1X),33(E,1X))


!     endif ! kount éééééééééééééééééééééééééééééééééééééééééééééééééééééééééééé

      RETURN
      END
