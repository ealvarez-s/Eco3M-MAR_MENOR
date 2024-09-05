      SUBROUTINE BIOLOGY
!______________________________________________________________________
! MODELE SYMPHONIE : VERSION 2007
! DERNIERE REVISION: 21 FEVRIER 2007
!      Equipe d'Océanographie Cotière du Laboratoire d'Aérologie       
!                Laboratoire d'Aérologie                               
!    Ckmax+1S - Université Paul Sabatier - Observatoire Midi Pyrénées      
!         14 Avenue Edouard Belin - 31400 Toulouse - FRANCE            
! Contact:     Patrick.Marsaleix@aero.obs-mip.fr                       
!               Claude.Estournel@aero.obs-mip.fr                       
!                Francis.Auclair@aero.obs-mip.fr
!______________________________________________________________________

      use module_principal
      use ModuleDeclaration
      use module_parallele !#MPI
      use module_global
      IMPLICIT NONE
      REAL*4 :: MLD
      INTEGER :: KEUPHO
      real sum_pt
      double precision :: FLUX_PPBi(kmax),FLUX_Respi(kmax),ppbint,nitrif_trans(kmax), &
!                          uptnit_trans(kmax)
                          ExcbactNH4NTOPLAYER,ExcbactNH4NINT,ExcbactNH4NDEEP,         &
                          ExcZooAmmoNTOPLAYER,ExcZooAmmoNINT,ExcZooAmmoNDEEP,           &
                          ExcZooPO4PTOPLAYER,ExcZooPO4PINT,ExcZooPO4PDEEP,           &
                          ExcbactPO4PTOPLAYER,ExcbactPO4PINT,ExcbactPO4PDEEP,         &
                          ExuSiTOPLAYER,ExuSiINT,ExuSiDEEP,           &
                          UptNitTOPLAYER,UptNitINT,UptNitDEEP,         &
                          UptAmmoTOPLAYER,UptAmmoINT,UptAmmoDEEP,         &
                          UptPTOPLAYER,UptPINT,UptPDEEP,         &
                          UptSiTOPLAYER,UptSiINT,UptSiDEEP,         &
                          NitrifTOPLAYER,NitrifINT,NitrifDEEP,              &
                          RemSMOPSiTOPLAYER,RemSMOPSiINT,RemSMOPSiDEEP,              &
                          RemLMOPSiTOPLAYER,RemLMOPSiINT,RemLMOPSiDEEP, &
                          UptbactPTOPLAYER,UptbactPINT,UptbactPDEEP, &
                          UptbactNTOPLAYER,UptbactNINT,UptbactNDEEP,bpCO2, &
                          zint,zdeep, &
                          FLUX_RespPi(kmax),FLUX_RespZi(kmax),FLUX_RespBi(kmax), &
                          gbac_trans(kmax),grazherbi_trans(kmax),graztot_trans(kmax),  &
                          FLUX_RespPpi(kmax),FLUX_RespPni(kmax),FLUX_RespPmi(kmax),  &
                          FLUX_PPBpi(kmax),FLUX_PPBni(kmax),FLUX_PPBmi(kmax), &
                          FLUX_UptNiti(kmax),FLUX_UptAmmoi(kmax),FLUX_UptPi(kmax),FLUX_UptPbi(kmax) , &
                          TOT_COL_N, TOT_COL_P, TOT_COL_Si , &
                          bEuphoticLayerDepth , &
                          bPAR(kmax)


      double precision :: mask_trans(kmax),bpH(kmax)
!     integer indice_bio(9),kk      !bio
!......................................................................
! modifs: 14/12/01: mise en service
!         06/04/06: modifs pour ECO3M
!         25/05/06: modif pour calcul epaisseur couche melange
!         21/02/07: include symphonie2007.h pour la version 2004
!         07/03/07: seuil sur variable BIO
!         04/05/07: seuil sur variable BIO 1D SMALL1
!         23/05/07: boucle sur MBIO1-1 
!         09/09/08: modif pour frontiere Est
!......................................................................



!******************************************************************************
! MISE A JOUR DES TENDANCES DU MODELE BIO AU PAS DE TEMPS BIO:
! DEBUT:
      if(  int(elapsedtime_now/dti_bio)       &
         /=int(elapsedtime_bef/dti_bio)) then !-------------->  !16-04-11

      if(par%rank==0)write(6,*)'passe dans biology'

!******************************************************************************


!******************************************************
! comptage du temps
!******************************************************
      TPS_PPB    = TPS_PPB    + 1
      TPS_PPB_2D = TPS_PPB_2D + 1
      TPS_PPB_3D = TPS_PPB_3D + 1
!=================================
! DEBUT 3D
!=================================


      if(I1DBIO.NE.1)then !->IF_I1DBIO
      do j1=nbio1,nbio2
      do i1=mbio1,mbio2

        if(mask_t(i1,j1,kmax+1).EQ.1) then !>>>>>>>>>>>>>>

!******************************************************
! APPEL A TENDECO
!******************************************************
      
!     if(i2.eq.899.and.j2.eq.121) then
!     do k=1,kmax
!     write(55,*)'pH SL avant tendeco',k,sPH(i1,j1,k),pCO2W(i1,j1)
!     enddo
!     endif


          do k=1,kmax
            tem_trans(k)=   tem_t(i1,j1,k,1)
            sal_trans(k)=   sal_t(i1,j1,k,1)
            rhp_trans(k)=   rhp_t(i1,j1,k)+rho
            ph_trans(k)=sPH(i1,j1,k)
            p3d_trans(k)=depth_t(i1,j1,k)
            pwk_trans(k)=depth_w(i1,j1,k)
            mask_trans(k)=mask_t(i1,j1,k)            
            do vb=1,vbmax
              bio_trans(k,vb)=MAX(bio_t(i1,j1,k,vb),SMALL1)               !07/03/07
            enddo
          enddo
          pwk_trans(kmax+1)=depth_w(i1,j1,kmax+1)
!           if (par%rank==0) write(6,*)'avant tendeco'
          CALL TENDECO( bio_trans             &
                       ,tem_trans             &
                       ,sal_trans             &
                       ,rhp_trans             &
                       ,dble(ssr_w(i1,j1,1))  &
                       ,albedo_w(i1,j1)       &
                       ,pwk_trans             &
                       ,p3d_trans             &
                       ,ph_trans              &
                       ,mask_trans            & !10
                       ,i1+par%timax(1)       &
                       ,j1+par%tjmax(1)       &
                       ,tdcbio_trans          &
                       ,PPB(1)                &
                       ,PPB(2)                &
                       ,PPB(3)                &
                       ,PPB(4)                &
                       ,NPB(1)                &
                       ,NPB(2)                &
                       ,NPB(3)                &
                       ,NPB(4)                &
                       ,RPB(1)                &
                       ,RPB(2)                &
                       ,RPB(3)                &
                       ,RPB(4)                &
                       ,RPB(5)                &
!                      ,NITRIF                &
                       ,RESPTOT               &
                       ,EXUCTOT               &
                       ,RESPPHYTO(1)          &
                       ,RESPPHYTO(2)          &
                       ,RESPPHYTO(3)          &
                       ,RESPPHYTO(4)          &
                       ,NETPPB(1)             &
                       ,NETPPB(2)             &
                       ,NETPPB(3)             &
                       ,NETPPB(4)             &
                       ,GRAZZOOC              &
                       ,GBAC                  &
                       ,FLUX_PPBi             &
                       ,FLUX_Respi            &
                       ,TOTALNITSURF          &
                       ,RESPZOO               &           
                       ,REMPOC                &            
                       ,REMPON                &          
                       ,LOSTPOC               &          
                       ,LOSTPON               &         
                       ,UPTBACTDON            &
!                      ,EXCHETERON            &
!                      ,ExcHeteroNtop         &
!                      ,ExcbactNH4top         &
!                      ,ExcBactNH4N           &
!                      ,ExcZooAmmoN           &
!                      ,ExczooNH4top          &
!                      ,NitrifCOLUMN          &
                       ,GRAZZOOCPUR           &
                       ,GRAZZOOCHERBI         &
                       ,Nitrif_trans          &
!                      ,UptNit_trans          &
                       ,ExcbactNH4NTOPLAYER   &
                       ,ExcbactNH4NINT        &
                       ,ExcbactNH4NDEEP       &
                       ,ExcZooAmmoNTOPLAYER   &
                       ,ExcZooAmmoNINT        &
                       ,ExcZooAmmoNDEEP       &
                       ,ExcZooPO4PTOPLAYER    &
                       ,ExcZooPO4PINT         &
                       ,ExcZooPO4PDEEP        &
                       ,ExcbactPO4PTOPLAYER   &
                       ,ExcbactPO4PINT        &
                       ,ExcbactPO4PDEEP       &
                       ,ExuSiTOPLAYER         &
                       ,ExuSiINT              &
                       ,ExuSiDEEP             &
                       ,UptNitTOPLAYER        &
                       ,UptNitINT             &
                       ,UptNitDEEP            &
                       ,UptAmmoTOPLAYER       &
                       ,UptAmmoINT            &
                       ,UptAmmoDEEP           &
                       ,UptPTOPLAYER          &
                       ,UptPINT               &
                       ,UptPDEEP              &
                       ,UptSiTOPLAYER         &
                       ,UptSiINT              &
                       ,UptSiDEEP             &
                       ,NitrifTOPLAYER        &
                       ,NitrifINT             &
                       ,NitrifDEEP            &
                       ,RemSMOPSiTOPLAYER     &
                       ,RemSMOPSiINT          &
                       ,RemSMOPSiDEEP         &
                       ,RemLMOPSiTOPLAYER     &
                       ,RemLMOPSiINT          &
                       ,RemLMOPSiDEEP         &
                       ,UptbactPTOPLAYER      &
                       ,UptbactPINT           &
                       ,UptbactPDEEP          &
                       ,UptbactNTOPLAYER      &
                       ,UptbactNINT           &
                       ,UptbactNDEEP          &
                       ,bpH                   &
                       ,bpCO2                 &
                       ,FLUX_RespPi           &
                       ,FLUX_RespZi           &
                       ,FLUX_RespBi           &
                       ,gbac_trans            &
                       ,grazherbi_trans       &
                       ,graztot_trans         &
                       ,FLUX_PPBpi            &
                       ,FLUX_PPBni            &
                       ,FLUX_PPBmi            &
                       ,FLUX_RespPpi          &
                       ,FLUX_RespPni          &
                       ,FLUX_RespPmi          &
                       ,FLUX_UptNiti          &
                       ,FLUX_UptAmmoi         &
                       ,FLUX_UptPi            &
                       ,FLUX_UptPbi           &
                       ,TOT_COL_N             &
                       ,TOT_COL_P             &
                       ,TOT_COL_Si            &
                       ,bEuphoticLayerDepth   &
                       ,bPAR                  )

          do vb=1,vbmax
          do k=1,kmax
            tendancebio_t(i1,j1,k,vb)=tdcbio_trans(k,vb)

!           if(lon_t(i1,j1)<-6.87*deg2rad) tendancebio_t(i1,j1,k,vb)=0. !2021/12/22  
!           if(lon_t(i,j)<-5.52*deg2rad) tendancebio_t(i1,j1,k,vb)=0. !2021/12/25

          enddo
          enddo

          do k=1,kmax
!           print*,i1,j1,k,oPH(k)
            sPH(i1,j1,k)=bpH(k)
            sPAR(i1,j1,k)=bPAR(k)
!           pCO2W(i1,j1,k)=oPCO2(k)
!           print*,'apres tendeco',k,sPH(i1,j1,k)
          enddo
          pCO2W(i1,j1)=bPCO2      

!       call equation_of_state('potential density',1)
!
!
!      data indice_bio /iSyneChl,iNanoChl,iDiaChl &
!                      ,iNitrate,iPhosphate,iOxygen &
!                      ,idic,iAlkalinity/
!
!      do kk=1,8 !- vb loop ->
!
!          vb=indice_bio(kk)
!
!      do k=1,kmax
!       if(vb==ioxygen.or.vb==idic) then
!       mean_vb(i1,j1,k,kk)=mean_vb(i1,j1,k,kk)+bio_t(i1,j1,k,vb) &
!                            /(rhp_t(i1,j1,k)+rho)*1000. !pour avoir des umol/kg
!       elseif(vb==ialkalinity) then
!       mean_vb(i1,j1,k,kk)=mean_vb(i1,j1,k,kk)+(bio_t(i1,j1,k,vb) &
!                            +bio_t(i1,j1,k,iAmmonium)  &
!                            -bio_t(i1,j1,k,iNitrate)   &
!                            -bio_t(i1,j1,k,iPhosphate)) &
!                            /(rhp_t(i1,j1,k)+rho)*1000. 
!       else
!       mean_vb(i1,j1,k,kk)=mean_vb(i1,j1,k,kk)+bio_t(i1,j1,k,vb)
!       endif
!
!      enddo !k
!
!      enddo !kk
!
!      do k=1,kmax 
!       mean_pH(i1,j1,k)=mean_pH(i1,j1,k)+sPH(i1,j1,k)
!      enddo
!
!
! 100   continue

!     print*,'tendecoen cours', i1,j1,par%rank,PPB(1),PPB(2),PPB(3),PPB(4)
!*****************************************************
!    CALCUL DES DIAGNOSTICS SUR MBIO1 MBIO2,NMIO1,NBIO2
!*****************************************************
          i2=i1+par%timax(1)      ! ds le repere global
          j2=j1+par%tjmax(1)

     if(i2.eq.195.and.j2.eq.202) then
!     do k=1,kmax
!     print*,'biolo apres',k,tendancebio_t(i1,j1,k,30)
!     enddo
     print*,'biolo apres',TOT_COL_N,TOT_COL_P,TOT_COL_Si
!     stop'biolo'
     endif

      TOT_COL_N_t(i1,j1)=TOT_COL_N_t(i1,j1)+TOT_COL_N
      TOT_COL_P_t(i1,j1)=TOT_COL_P_t(i1,j1)+TOT_COL_P
      TOT_COL_Si_t(i1,j1)=TOT_COL_Si_t(i1,j1)+TOT_COL_Si

      EuphoticLayerDepth_t(i1,j1)=EuphoticLayerDepth_t(i1,j1)+bEuphoticLayerDepth
    
!      if(i2.eq.335.and.j2.eq.32) then
!      write(55,*),UptPTOPLAYER,ExcbactPO4PTOPLAYER,ExcZooPO4PTOPLAYER,UptbactPTOPLAYER
!      endif

          do k=1,4
            PPB2D(i1,j1,k) =    PPB2D(i1,j1,k) +    PPB(K)      ! en mgC/m2/j
          enddo
!          if(par%rank==0)write(6,*)'apres PPB2D'

          do k=1,4
            NPB2D(i1,j1,k) = NPB2D(i1,j1,k) + NPB(K)
          enddo

! Net Primary Production
          do k=1,4
            NETPPB2D(i1,j1,k) = NETPPB2D(i1,j1,k) + NETPPB(K)
          enddo

! Recycled production
!      print*,'RPB'
          do k=1,5
            if(abs(rpb(k))>1.e-30)then
              RPB2D(i1,j1,k) = RPB2D(i1,j1,k) + RPB(K)
            endif
          enddo

          GBAC2D(i1,j1)= GBAC2D(i1,j1) + GBAC
!         if(par%rank==0)write(6,*)'apres GBACD'

      
          CHLsurf(i1,j1)= CHLsurf(i1,j1) + bio_t(i1,j1,kmax,isynechl) &
                                         + bio_t(i1,j1,kmax,inanochl) &
                                         + bio_t(i1,j1,kmax,idiachl) 

          CHLpsurf(i1,j1)= CHLpsurf(i1,j1) + bio_t(i1,j1,kmax,isynechl) 
          CHLnsurf(i1,j1)= CHLnsurf(i1,j1) + bio_t(i1,j1,kmax,inanochl) 
          CHLmsurf(i1,j1)= CHLmsurf(i1,j1) + bio_t(i1,j1,kmax,idiachl)

          NITsurf(i1,j1)= NITsurf(i1,j1) + bio_t(i1,j1,kmax,initrate) 

          PHOsurf(i1,j1)= PHOsurf(i1,j1) + bio_t(i1,j1,kmax,iphosphate)

! Integration sur les couches de surface et intermediaires
! Inverse ("1 sur") de l'epaisseur des couches du bilan
          zone4_inv_dz=1./150.
! Cas particulier zone4_nlayer=1: zone4_inv_dz=1./hmax

!         zone4_stretch_dz=1.  ! epaisseur constante = 1/zone4_inv_dz
          zone4_stretch_dz=0.706

          zint=-real(1)**(1./zone4_stretch_dz)/zone4_inv_dz
          zdeep=-real(2)**(1./zone4_stretch_dz)/zone4_inv_dz

!      print*,'verif biology',zint,zdeep

          do k=kmax+1,2,-1  !->Kloop_kmaxp1->2
            if(mask_t(i1,j1,k)==1) then !->If_mask_t

              if(depth_w(i1,j1,k)>zint) then
                CHLSW(i1,j1)= CHLSW(i1,j1) + (bio_t(i1,j1,k-1,isynechl) &
                                           + bio_t(i1,j1,k-1,inanochl)  &
                                           + bio_t(i1,j1,k-1,idiachl))  &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CHLpSW(i1,j1)= CHLpSW(i1,j1) + bio_t(i1,j1,k-1,isynechl) &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CHLnSW(i1,j1)= CHLnSW(i1,j1) + bio_t(i1,j1,k-1,inanochl) &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CHLmSW(i1,j1)= CHLmSW(i1,j1) + bio_t(i1,j1,k-1,idiachl) &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)



                NITSW(i1,j1)= NITSW(i1,j1) + bio_t(i1,j1,k-1,initrate) &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                PHOSW(i1,j1)= PHOSW(i1,j1) + bio_t(i1,j1,k-1,iphosphate) &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                GPPSW(i1,j1)= GPPSW(i1,j1) + FLUX_PPBi(k-1) &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CRPSW(i1,j1)= CRPSW(i1,j1) + FLUX_RespPi(k-1)          &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CRZSW(i1,j1)= CRZSW(i1,j1) + FLUX_RespZi(k-1)          &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CRBSW(i1,j1)= CRBSW(i1,j1) + FLUX_RespBi(k-1)          & 
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                NITRIFSW(i1,j1)= NITRIFSW(i1,j1) + nitrif_trans(k-1)         & 
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                GBSW(i1,j1)= GBSW(i1,j1) + gbac_trans(k-1)             &          
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                GrazHerbiSW(i1,j1)= GrazHerbiSW(i1,j1) + grazherbi_trans(k-1) &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                GrazTotSW(i1,j1)= GrazHerbiSW(i1,j1) + graztot_trans(k-1) &
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)


                GPPPSW(i1,j1)= GPPPSW(i1,j1) + FLUX_PPBpi(k-1) & 
                                   * (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CRPPSW(i1,j1)= CRPPSW(i1,j1) + FLUX_RespPpi(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                GPPNSW(i1,j1)= GPPNSW(i1,j1) + FLUX_PPBni(k-1) &                       
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CRPNSW(i1,j1)= CRPNSW(i1,j1) + FLUX_RespPni(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                GPPMSW(i1,j1)= GPPMSW(i1,j1) + FLUX_PPBmi(k-1) &                       
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                CRPMSW(i1,j1)= CRPMSW(i1,j1) + FLUX_RespPmi(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)


                UPTNITSW(i1,j1)= UPTNITSW(i1,j1) + FLUX_UptNiti(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                UPTAMMOSW(i1,j1)= UPTAMMOSW(i1,j1) + FLUX_UptAmmoi(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)


                UPTPSW(i1,j1)= UPTPSW(i1,j1) + FLUX_UptPi(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                UPTPBSW(i1,j1)= UPTPBSW(i1,j1) + FLUX_UptPbi(k-1)          &
                                   *(depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

                dzSW(i1,j1)= dzSW(i1,j1) + (depth_w(i1,j1,k)-max(depth_w(i1,j1,k-1),zint))*mask_t(i1,j1,k-1)

              endif


              if(depth_w(i1,j1,k)>zdeep.and.depth_w(i1,j1,k)<zint)then
                GPPIW(i1,j1)= GPPIW(i1,j1) + FLUX_PPBi(k-1) &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1) 

                CRPIW(i1,j1)= CRPIW(i1,j1) + FLUX_RespPi(k-1) &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                CRZIW(i1,j1)= CRZIW(i1,j1) + FLUX_RespZi(k-1)          &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                CRBIW(i1,j1)= CRBIW(i1,j1) + FLUX_RespBi(k-1)          &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                NITRIFIW(i1,j1)= NITRIFIW(i1,j1) + nitrif_trans(k-1)         &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                GBIW(i1,j1)= GBIW(i1,j1) + gbac_trans(k-1)             &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                GrazHerbiIW(i1,j1)= GrazHerbiIW(i1,j1) + grazherbi_trans(k-1) &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                GrazTotIW(i1,j1)= GraztotSW(i1,j1) + graztot_trans(k-1) &
                                   * (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                dzIW(i1,j1)= dzIW(i1,j1)+ (min(depth_w(i1,j1,k),zint)-max(depth_w(i1,j1,k-1),zdeep))*mask_t(i1,j1,k-1)

                if(depth_w(i1,j1,k+1)>zint) then
                  GPPIW(i1,j1)= GPPIW(i1,j1) + FLUX_PPBi(k) *(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  CRPIW(i1,j1)= CRPIW(i1,j1) + FLUX_RespPi(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  CRZIW(i1,j1)= CRZIW(i1,j1) + FLUX_RespZi(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  CRBIW(i1,j1)= CRBIW(i1,j1) + FLUX_RespBi(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  NITRIFIW(i1,j1)= NITRIFIW(i1,j1) + nitrif_trans(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  GBIW(i1,j1)= GBIW(i1,j1) + gbac_trans(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  GrazHerbiIW(i1,j1)= GrazHerbiIW(i1,j1) + grazherbi_trans(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  GrazTotIW(i1,j1)= GraztotSW(i1,j1) + graztot_trans(k)*(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                  dzIW(i1,j1)=dzIW(i1,j1)+(zint-max(depth_w(i1,j1,k),zdeep))*mask_t(i1,j1,k)
                endif


              endif

            endif !->if_mask_t
          enddo !->Kloop_kmaxp1->2

! COMMENTE ALEX 11/01/18
! Nitrification, Respiration, Exudation
!      NITRIF2D(i1,j1) = NITRIF2D(i1,j1)  + NITRIF
!      NITRIFCOLUMN2D(i1,j1) = NITRIFCOLUMN2D(i1,j1)  + NITRIFCOLUMN
!
          do k=1,kmax
            ppb3d(i1,j1,k)=ppb3d(i1,j1,k)+FLUX_PPBi(k)
            resp3d(i1,j1,k)=resp3d(i1,j1,k)+FLUX_Respi(k)
            nitrif3d(i1,j1,k)=nitrif3d(i1,j1,k)+nitrif_trans(k)/14 &
                                               /(depth_w(i1,j1,k+1)-depth_w(i1,j1,k))
!            uptnit3d(i1,j1,k)=uptnit3d(i1,j1,k)+uptnit_trans(k)
          enddo      


          RESP2D(i1,j1)   = RESP2D(i1,j1) + RESPTOT  



!         EXCHETERON2D(i1,j1)   = EXCHETERON2D(i1,j1) + EXCHETERON
!
!          ExcHeteroNtop2D(i1,j1)   = ExcHeteroNtop2D(i1,j1) + ExcHeteroNtop
!          ExcbactNH4top2D(i1,j1)   = ExcbactNH4top2D(i1,j1) + ExcbactNH4top
!
!          ExcBactNH4N2D(i1,j1)   = ExcBactNH4N2D(i1,j1) + ExcBactNH4N
!          ExcZooAmmoN2D(i1,j1)   = ExcZooAmmoN2D(i1,j1) + ExcZooAmmoN
!
!          ExczooNH4top2D(i1,j1)   = ExczooNH4top2D(i1,j1) + ExczooNH4top

          ExcbactNH4NTOPLAYER2D(i1,j1) = ExcbactNH4NTOPLAYER2D(i1,j1) + ExcbactNH4NTOPLAYER
          ExcZooAmmoNTOPLAYER2D(i1,j1) = ExcZooAmmoNTOPLAYER2D(i1,j1) + ExcZooAmmoNTOPLAYER
          ExcZooPO4PTOPLAYER2D(i1,j1)  = ExcZooPO4PTOPLAYER2D(i1,j1) + ExcZooPO4PTOPLAYER
          ExcbactPO4PTOPLAYER2D(i1,j1) = ExcbactPO4PTOPLAYER2D(i1,j1) + ExcbactPO4PTOPLAYER
          ExuSiTOPLAYER2D(i1,j1)       = ExuSiTOPLAYER2D(i1,j1) + ExuSiTOPLAYER
          UptNitTOPLAYER2D(i1,j1)      = UptNitTOPLAYER2D(i1,j1) + UptNitTOPLAYER
          UptAmmoTOPLAYER2D(i1,j1)     = UptAmmoTOPLAYER2D(i1,j1) + UptAmmoTOPLAYER 
          UptPTOPLAYER2D(i1,j1)        = UptPTOPLAYER2D(i1,j1) + UptPTOPLAYER
          UptSiTOPLAYER2D(i1,j1)       = UptSiTOPLAYER2D(i1,j1) + UptSiTOPLAYER
          NitrifTOPLAYER2D(i1,j1)      = NitrifTOPLAYER2D(i1,j1) + NitrifTOPLAYER
          RemSMOPSiTOPLAYER2D(i1,j1)   = RemSMOPSiTOPLAYER2D(i1,j1) + RemSMOPSiTOPLAYER
          RemLMOPSiTOPLAYER2D(i1,j1)   = RemLMOPSiTOPLAYER2D(i1,j1) + RemLMOPSiTOPLAYER

          UptBactPTOPLAYER2D(i1,j1)    = UptBactPTOPLAYER2D(i1,j1) + UptbactPTOPLAYER
          UptBactNTOPLAYER2D(i1,j1)    = UptBactNTOPLAYER2D(i1,j1) + UptbactNTOPLAYER

          ExcbactNH4NINT2D(i1,j1) = ExcbactNH4NINT2D(i1,j1) + ExcbactNH4NINT
          ExcZooAmmoNINT2D(i1,j1) = ExcZooAmmoNINT2D(i1,j1) + ExcZooAmmoNINT
          ExcZooPO4PINT2D(i1,j1)  = ExcZooPO4PINT2D(i1,j1) + ExcZooPO4PINT
          ExcbactPO4PINT2D(i1,j1) = ExcbactPO4PINT2D(i1,j1) + ExcbactPO4PINT
          ExuSiINT2D(i1,j1)       = ExuSiINT2D(i1,j1) + ExuSiINT
          UptNitINT2D(i1,j1)      = UptNitINT2D(i1,j1) + UptNitINT
          UptAmmoINT2D(i1,j1)     = UptAmmoINT2D(i1,j1) + UptAmmoINT
          UptPINT2D(i1,j1)        = UptPINT2D(i1,j1) + UptPINT
          UptSiINT2D(i1,j1)       = UptSiINT2D(i1,j1) + UptSiINT
          NitrifINT2D(i1,j1)      = NitrifINT2D(i1,j1) + NitrifINT
          RemSMOPSiINT2D(i1,j1)   = RemSMOPSiINT2D(i1,j1) + RemSMOPSiINT
          RemLMOPSiINT2D(i1,j1)   = RemLMOPSiINT2D(i1,j1) + RemLMOPSiINT

          UptBactPINT2D(i1,j1)   = UptBactPINT2D(i1,j1) + UptbactPINT
          UptBactNINT2D(i1,j1)   = UptBactNINT2D(i1,j1) + UptbactNINT

          ExcbactNH4NDEEP2D(i1,j1)  = ExcbactNH4NDEEP2D(i1,j1) + ExcbactNH4NDEEP
          ExcZooAmmoNDEEP2D(i1,j1)  = ExcZooAmmoNDEEP2D(i1,j1) + ExcZooAmmoNDEEP
          ExcZooPO4PDEEP2D(i1,j1)   = ExcZooPO4PDEEP2D(i1,j1) + ExcZooPO4PDEEP
          ExcbactPO4PDEEP2D(i1,j1)  = ExcbactPO4PDEEP2D(i1,j1) + ExcbactPO4PDEEP
          ExuSiDEEP2D(i1,j1)        = ExuSiDEEP2D(i1,j1) + ExuSiDEEP
          UptNitDEEP2D(i1,j1)       = UptNitDEEP2D(i1,j1) + UptNitDEEP 
          UptAmmoDEEP2D(i1,j1)      = UptAmmoDEEP2D(i1,j1) + UptAmmoDEEP
          UptPDEEP2D(i1,j1)         = UptPDEEP2D(i1,j1) + UptPDEEP 
          UptSiDEEP2D(i1,j1)        = UptSiDEEP2D(i1,j1) + UptSiDEEP
          NitrifDEEP2D(i1,j1)       = NitrifDEEP2D(i1,j1) + NitrifDEEP 
          RemSMOPSiDEEP2D(i1,j1)    = RemSMOPSiDEEP2D(i1,j1) + RemSMOPSiDEEP
          RemLMOPSiDEEP2D(i1,j1)    = RemLMOPSiDEEP2D(i1,j1) + RemLMOPSiDEEP

          UptBactPDEEP2D(i1,j1)  = UptBactPDEEP2D(i1,j1) + UptbactPDEEP
          UptBactNDEEP2D(i1,j1)  = UptBactNDEEP2D(i1,j1) + UptbactNDEEP

          GRAZ2D(i1,j1)     = GRAZ2D(i1,j1) + GRAZZOOC
          GRAZPUR2D(i1,j1)  = GRAZPUR2D(i1,j1) + GRAZZOOCPUR    !claude
          GRAZHERBI2D(i1,j1)= GRAZHERBI2D(i1,j1) + GRAZZOOCHERBI  !claude
!      if(par%rank==0)write(6,*)'apres GRAZHERBI'

!******************************************************
! fin diagnostics bio 
!******************************************************


        endif                          !>>>>>>>>>>>>>>
      enddo    ! fin de boucle i1
      enddo    ! fin de boucle j1
!       if (par%rank==0) write(6,*)'apres tendeco'


!******************************************************
! Nudging Begin                                                       !18/07/08
!******************************************************

!     if(Inudging.EQ.1) then
!     if(KOUNT.GT.knudgBegin.and.kount<=knudgEnd) then

! West boundary
!     DO I1=MAX(MBIO1-5,1),MBIO1-1 
!     DO J1=MAX(NBIO1-5,1),NBIO2

! All the domain
!      DO I1=MAX(MBIO1-5,1),MBIO2
!      DO J1=MAX(NBIO1-5,1),NBIO2


!     do k=1,kmax


!     tendancebio_t(i1,j1,k,iNitrate)=tendancebio_t(i1,j1,k,iNitrate)
!    1  +(BIOINIT_Z(i1,j1,k,iNitrate)-bio_t(i1,j1,k,iNitrate))
!    1   *Nudgperiod /24. /3600. 

!     tendancebio_t(i1,j1,k,iPhosphate)=
!    1    tendancebio_t(i1,j1,k,iPhosphate)
!    1  +(BIOINIT_Z(i1,j1,k,iPhosphate)-bio_t(i1,j1,k,iPhosphate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iSilice)=tendancebio_t(i1,j1,k,iSilice)
!    1  +(BIOINIT_Z(i1,j1,k,iSilice)-bio_t(i1,j1,k,iSilice))
!    1   *Nudgperiod /24. /3600.

 

!     enddo

!     enddo    ! fin de boucle i1
!     enddo    ! fin de boucle j1

! South boundary
!     DO I1=MBIO1,MBIO2
!     DO J1=MAX(NBIO1-5,1),NBIO1-1

!     do k=1,kmax
!     tendancebio_t(i1,j1,k,iNitrate)=tendancebio_t(i1,j1,k,iNitrate)
!    1  +(BIOINIT_Z(i1,j1,k,iNitrate)-bio_t(i1,j1,k,iNitrate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iPhosphate)
!    1  =tendancebio_t(i1,j1,k,iPhosphate)
!    1  +(BIOINIT_Z(i1,j1,k,iPhosphate)-bio_t(i1,j1,k,iPhosphate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iSilice)=tendancebio_t(i1,j1,k,iSilice)
!    1  +(BIOINIT_Z(i1,j1,k,iSilice)-bio_t(i1,j1,k,iSilice))
!    1   *Nudgperiod /24. /3600.
!     enddo

!     enddo    ! fin de boucle i1
!     enddo    ! fin de boucle j1

! East boundary
!     DO I1=MBIO2+1,MIN0(MBIO2+5,MBIO)
!     DO J1=NBIO1,NBIO2

!     do k=1,kmax
!     tendancebio_t(i1,j1,k,iNitrate)=tendancebio_t(i1,j1,k,iNitrate)
!    1  +(BIOINIT_Z(i1,j1,k,iNitrate)-bio_t(i1,j1,k,iNitrate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iPhosphate)
!    1  =tendancebio_t(i1,j1,k,iPhosphate)
!    1  +(BIOINIT_Z(i1,j1,k,iPhosphate)-bio_t(i1,j1,k,iPhosphate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iSilice)=tendancebio_t(i1,j1,k,iSilice)
!    1  +(BIOINIT_Z(i1,j1,k,iSilice)-bio_t(i1,j1,k,iSilice))
!    1   *Nudgperiod /24. /3600.
!     enddo

!     enddo    ! fin de boucle i1
!     enddo    ! fin de boucle j1



!     endif
!     endif


!******************************************************
! Nudging End                                                          ! 18/07/08     
!******************************************************
      endif !->IF_I1DBIO


!      print*,'apres flux'

 103  FORMAT(29(F18.4,1X)) 


!=================================
! FIN 3D
!=================================

!=================================
! DEBUT 1D
!=================================


      if(I1DBIO.EQ.1)then !->model_1DBIO
        I1=IPTBIO
        J1=JPTBIO

!      MLD=-depth_w(i1,j1,kmax)
!      do k=kmax,1,-1                                   ! 25/05/06
!      if(kz_w(i1,j1,k).LT.1.E-4)then
!       MLD=-depth_w(i1,j1,k)
!       GOTO 16
!      endif
!      enddo
! 16   continue


        do k=1,kmax
          tem_trans(k)=   tem_t(i1,j1,k,1)
          p3d_trans(k)=depth_t(i1,j1,k)
          pwk_trans(k)=depth_w(i1,j1,k)
          mask_trans(k)=mask_t(i1,j1,k)
          do vb=1,vbmax
            bio_trans(k,vb)=MAX(bio_t(i1,j1,k,vb),SMALL1)
          enddo
        enddo
        pwk_trans(kmax+1)=depth_w(i1,j1,kmax+1)


        CALL TENDECO( bio_trans             &
                     ,tem_trans             &
                     ,sal_trans             &
                     ,rhp_trans             &
                     ,dble(ssr_w(i1,j1,1))   &
                     ,albedo_w(i1,j1)       &
                     ,pwk_trans             &
                     ,p3d_trans             &
                     ,ph_trans              &
                     ,mask_trans            &
                     ,tdcbio_trans          &
                     ,PPB(1)                &
                     ,PPB(2)                &
                     ,PPB(3)                &
                     ,PPB(4)                &
                     ,NPB(1)                &
                     ,NPB(2)                &
                     ,NPB(3)                &
                     ,NPB(4)                &
                     ,RPB(1)                &
                     ,RPB(2)                &
                     ,RPB(3)                &
                     ,RPB(4)                &
                     ,RPB(5)                &
!                    ,NITRIF                &
                     ,RESPTOT               &
                     ,EXUCTOT               &
                     ,RESPPHYTO(1)          &
                     ,RESPPHYTO(2)          &
                     ,RESPPHYTO(3)          &
                     ,RESPPHYTO(4)          &
                     ,NETPPB(1)             &
                     ,NETPPB(2)             &
                     ,NETPPB(3)             &
                     ,NETPPB(4)             &
                     ,GRAZZOOC              &
                     ,GBAC                  &
                     ,FLUX_PPBi             &
                     ,FLUX_Respi            &
                     ,TOTALNITSURF          &
                     ,RESPZOO               &           
                     ,REMPOC                &            
                     ,REMPON                &          
                     ,LOSTPOC               &          
                     ,LOSTPON               &         
                     ,UPTBACTDON            &
!                    ,EXCHETERON            &
!                    ,ExcHeteroNtop         &
!                    ,ExcbactNH4top         &
!                    ,ExcBactNH4N           &
!                    ,ExcZooAmmoN           &
!                    ,ExczooNH4top          &
!                    ,NitrifCOLUMN          &
                     ,GRAZZOOCPUR           &
                     ,GRAZZOOCHERBI         &
                     ,Nitrif_trans          &
!                    ,UptNit_trans          &
                     ,ExcbactNH4NTOPLAYER   &
                     ,ExcbactNH4NINT        &
                     ,ExcbactNH4NDEEP       &
                     ,ExcZooAmmoNTOPLAYER   &
                     ,ExcZooAmmoNINT        &
                     ,ExcZooAmmoNDEEP       &
                     ,ExcZooPO4PTOPLAYER    &
                     ,ExcZooPO4PINT         &
                     ,ExcZooPO4PDEEP        &
                     ,ExcbactPO4PTOPLAYER   &
                     ,ExcbactPO4PINT        &
                     ,ExcbactPO4PDEEP       &
                     ,ExuSiTOPLAYER         &
                     ,ExuSiINT              &
                     ,ExuSiDEEP             &
                     ,UptNitTOPLAYER        &
                     ,UptNitINT             &
                     ,UptNitDEEP            &
                     ,UptAmmoTOPLAYER       &
                     ,UptAmmoINT            &
                     ,UptAmmoDEEP           &
                     ,UptPTOPLAYER          &
                     ,UptPINT               &
                     ,UptPDEEP              &
                     ,UptSiTOPLAYER         &
                     ,UptSiINT              &
                     ,UptSiDEEP             &
                     ,NitrifTOPLAYER        &
                     ,NitrifINT             &
                     ,NitrifDEEP            &
                     ,RemSMOPSiTOPLAYER     &
                     ,RemSMOPSiINT          &
                     ,RemSMOPSiDEEP         &
                     ,RemLMOPSiTOPLAYER     &
                     ,RemLMOPSiINT          &
                     ,RemLMOPSiDEEP         &
                     ,UptbactPTOPLAYER      &
                     ,UptbactPINT           &
                     ,UptbactPDEEP          &
                     ,UptbactNTOPLAYER      &
                     ,UptbactNINT           &
                     ,UptbactNDEEP          &
                     ,bpH                   &
                     ,bpCO2                 & 
                     ,FLUX_RespPi           &
                     ,FLUX_RespZi           &
                     ,FLUX_RespBi           &
                     ,gbac_trans            &  
                     ,grazherbi_trans       &  
                     ,graztot_trans         &
                     ,FLUX_PPBpi            &
                     ,FLUX_PPBni            &
                     ,FLUX_PPBmi            &
                     ,FLUX_RespPpi          &
                     ,FLUX_RespPni          &
                     ,FLUX_RespPmi          &
                     ,FLUX_UptNiti          &
                     ,FLUX_UptAmmoi         &
                     ,FLUX_UptPi            &
                     ,FLUX_UptPbi           &
                     ,TOT_COL_N             &
                     ,TOT_COL_P             &
                     ,TOT_COL_Si            &
                     ,bEuphoticLayerDepth   &
                     ,bPAR        )


        do vb=1,vbmax
        do k=1,kmax
          tendancebio_t(i1,j1,k,vb)=tdcbio_trans(k,vb)
        enddo
        enddo


! Nudging: beginning
!     do k=1,kmax
!     tendancebio_t(i1,j1,k,iNitrate)=tendancebio_t(i1,j1,k,iNitrate)
!    1  +(BIOINIT_Z(i1,j1,k,iNitrate)-bio_t(i1,j1,k,iNitrate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iPhosphate)=
!    1    tendancebio_t(i1,j1,k,iPhosphate)
!    1  +(BIOINIT_Z(i1,j1,k,iPhosphate)-bio_t(i1,j1,k,iPhosphate))
!    1   *Nudgperiod /24. /3600.

!     tendancebio_t(i1,j1,k,iSilice)=tendancebio_t(i1,j1,k,iSilice)
!    1  +(BIOINIT_Z(i1,j1,k,iSilice)-bio_t(i1,j1,k,iSilice))
!    1   *Nudgperiod /24. /3600.
!     enddo
! Nudging: end

      endif !->model_1DBIO

!=================================
! FIN 1D
!=================================


!      if(par%rank==0)write(6,*)'sors de biology'
!******************************************************************************
! MISE A JOUR DES TENDANCES DU MODELE BIO AU PAS DE TEMPS BIO:
! FIN.
      endif !-------------->
!******************************************************************************
      RETURN
    
      END
