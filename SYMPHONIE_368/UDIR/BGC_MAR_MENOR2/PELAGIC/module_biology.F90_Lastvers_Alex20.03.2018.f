      module module_biology
!.......................................................................
! ECO3M-S model
! release S26.1 - last update: 12-02-17
!.......................................................................
! S26  12-02-17  wsed dimension verticale
!.......................................................................
      use module_parameter

      double precision,allocatable,dimension(:,:,:,:) ::                &
       bio_t

      double precision,allocatable,dimension(:,:,:,:) ::                &
       biosum_t

      double precision,allocatable,dimension(:,:,:,:) ::                &
       tendancebio_t

      double precision,allocatable,dimension(:,:,:,:) ::                &
       fluxbio_w

      double precision,allocatable,dimension(:,:,:,:) ::                &
       nitobc_t,                                                        &
       ammobc_t,                                                        &
       phoobc_t,                                                        &
       silobc_t

      double precision,allocatable,dimension(:,:,:,:,:)  ::             &
       biobc_i_t

      double precision,allocatable,dimension(:,:,:,:,:) ::              &
       biobc_j_t

!      double precision,allocatable,dimension(:,:,:) ::                  &
!       nitrif3d                                                         &
!      ,uptnit3d

      double precision,allocatable,dimension(:,:) ::                    &
!       nitrif2d                                                         &
      resp2d                                                           &
      ,graz2d                                                           &
      ,gbac2d                                                           &
      ,exuctot2d                                                        &
      ,sumnitratesurf                                                   & 
      ,ppbsurf2d                                                        &
      ,ppbmax2d                                                         &
      ,profppbmax2d                                                     &
      ,sum_exportc_bot_2d                                               &
      ,sum_exportn_bot_2d                                               &
      ,sum_exportp_bot_2d                                               &
      ,sum_exportsi_bot_2d                                              &
!      ,NITRIFCOLUMN2D                                                   &
!      ,EXCHETERON2D                                                     &
!      ,ExcHeteroNtop2D                                                  &
!      ,ExcbactNH4top2D                                                  &
!      ,ExcBactNH4N2D                                                    &
!      ,ExcZooAmmoN2D                                                    &
!      ,ExczooNH4top2D                                                   &
      ,GRAZHERBI2D                                                      &
      ,GRAZPUR2D                                                        &
!      ,netppb2d                                                         &
      ,NO3efflux2d                                                      &
      ,NH4efflux2d                                                      &
      ,Pefflux2d                                                        &
      ,Siefflux2d                                                       &
      ,ExcbactNH4NTOPLAYER2D                                            &
      ,ExcZooAmmoNTOPLAYER2D                                            & 
      ,ExcZooPO4PTOPLAYER2D                                             &
      ,ExcbactPO4PTOPLAYER2D                                            &
      ,ExuSiTOPLAYER2D                                                  &       
      ,UptNitTOPLAYER2D                                                 &
      ,UptAmmoTOPLAYER2D                                                &
      ,UptPTOPLAYER2D                                                   &
      ,UptSiTOPLAYER2D                                                  &
      ,NitrifTOPLAYER2D                                                 &
      ,RemSMOPSiTOPLAYER2D                                              &
      ,RemLMOPSiTOPLAYER2D                                              &
      ,ExcbactNH4NINT2D                                            &
      ,ExcZooAmmoNINT2D                                            &
      ,ExcZooPO4PINT2D                                             &
      ,ExcbactPO4PINT2D                                            &
      ,ExuSiINT2D                                                  &      
      ,UptNitINT2D                                                 &
      ,UptAmmoINT2D                                                &
      ,UptPINT2D                                                   &
      ,UptSiINT2D                                                  &
      ,NitrifINT2D                                                 &
      ,RemSMOPSiINT2D                                              &
      ,RemLMOPSiINT2D                                              &
      ,ExcbactNH4NDEEP2D                                            &
      ,ExcZooAmmoNDEEP2D                                            &
      ,ExcZooPO4PDEEP2D                                             &
      ,ExcbactPO4PDEEP2D                                            &
      ,ExuSiDEEP2D                                                  &      
      ,UptNitDEEP2D                                                 &
      ,UptAmmoDEEP2D                                                &
      ,UptPDEEP2D                                                   &
      ,UptSiDEEP2D                                                  &
      ,NitrifDEEP2D                                                 &
      ,RemSMOPSiDEEP2D                                              &
      ,RemLMOPSiDEEP2D                                              
                                                     
      double precision,allocatable,dimension(:,:) ::                    &
       CBDet,                                                           &
       CBFDet,                                                          &
       CBSDet,                                                          &
       NBDet,                                                           &
       PBDet,                                                           &
       SiBDet,                                                          &
       dCBDet,                                                          &
       dCBFDet,                                                         &
       dCBSDet,                                                         &
       dNBDet,                                                          &
       dPBDet,                                                          &
       dSiBDet,                                                         &
       CDepo,                                                           &
       NDepo,                                                           &
       PDepo,                                                           &
       SiDepo,                                                          &
       NitrificationB,                                                  &
       DenitrificationB,                                                &
       O2Min,                                                           &
       AnoxMin,                                                         &
       O2ODU,                                                           &
       SUMT_NITRIF_est,                                                     &
       SUMT_DENITRIF_est,                                                   &
       SUMT_NO3Flux_est,                                                    &
       SUMT_NH3Flux_est,                                                    &
       SUMT_NMin_est,                                                       &
       SUMT_CMin_est,                                                       &
       SUMT_NDepo_est,                                                      &
       SUMT_CDepo_est,                                                      &
       SUMT_NITRIF,                                                     &
       SUMT_DENITRIF,                                                   &
       SUMT_NO3Flux,                                                    &
       SUMT_NH3Flux,                                                    &
       SUMT_NMin,                                                       &
       SUMT_CMin,                                                       &
       SUMT_NDepo,                                                      &
       SUMT_CDepo,                                                      &
       CMin_out,                                                        &
       NMin_out,                                                        &
       SUMT_O2Min,                                                      &
       SUMT_AnoxMin,                                                    &
       SUMT_AnoxMin_est


      double precision,allocatable,dimension(:,:,:) ::                  &
       ppb2d,                                                           &
       npb2d,                                                           &
       netppb2d,                                                           &
       chl2d

      double precision,allocatable,dimension(:,:,:) ::                  &
       rpb2d

      double precision,allocatable,dimension(:,:,:) ::                  &
       exp2d



      double precision,allocatable,dimension(:,:) ::                    &
       bio_trans,                                                       &
       wsed,                                                            &
       tdcbio_trans

      double precision,allocatable,dimension(:) ::                      &
       pwk_trans

      double precision,allocatable,dimension(:) ::                      &
       tem_trans,                                                       &
       p3d_trans                                                       

      integer,allocatable,dimension(:) ::                               &
       radionucleide


!      double precision,allocatable,dimension(:,:,:) ::                 &
!      par_bio

      double precision,allocatable,dimension(:) ::                      &
!      wsed,                                                            &
       TauR,                                                            &
       alpha_radio,                                                     &
       radio_coef,                                                      &
       export_plat,                                                     &
       sum_export_plat,                                                 &
       exportgdl_plat,                                                  &
       sum_exportgdl_plat,                                              &
       export_200_1000m,                                                &
       sum_export_200_1000m,                                            &
       bil_fleuves,                                                     &
       bil_fleuves_est,                                                 &
       sum_ob_west,                                                     &
       sum_ob_south,                                                    &
       flux_ob_south,                                                   &
       flux_ob_west,                                                    &
       exportdyf_bot,                                                   &
       sum_exportdyf_bot,                                               &
       ins_fleuves,                                                     &
       ins_fleuves_est,                                                     &
       sum_ob_east,                                                     &
       sum_ob_north,                                                    &
       flux_ob_east,                                                    &
       flux_ob_north,                                                   &
       sumbio_p,                                                        &
       sumbion_p,                                                       &
       par_sum_ob_west,                                                 &
       par_sum_ob_south,                                                &
       par_sum_ob_east,                                                 &
       par_sum_export_plat,                                             &
       par_sum_exportgdl_plat,                                          &
       par_sum_exportdyf_bot,                                           &
       par_sum_export_200_1000m,                                        &
       par_bil_fleuves,                                                 &
       exportdiag_bot,                                                  &
       sum_exportdiag_bot,                                              &
       par_sum_exportdiag_bot,                                          &
       bil_fleuves_area,                                                &
       bil_fleuves_area_est
      double precision,allocatable,dimension(:) ::                      &
       exportdyf_100,                                                   &
       sum_exportdyf_100,                                               &
       par_sum_exportdyf_100


      double precision,allocatable,dimension(:) ::                      &
       export_200,                                                      &
       sum_export_200,                                                  &
       exportgdl_200,                                                   &
       sum_exportgdl_200,                                               &
       exportdyf_200,                                                   &
       sum_exportdyf_200,                                               &
       exportdyf3_200,                                                  &
       sum_exportdyf3_200,                                              &
       exportdyf_1000,                                                  &
       sum_exportdyf_1000,                                              &
       par_sum_export_200,                                              &
       par_sum_exportgdl_200,                                           &
       par_sum_exportdyf_200,                                           &
       par_sum_exportdyf3_200,                                          &
       par_sum_exportdyf_1000,                                          &
       exportdiag_200,                                                  &
       sum_exportdiag_200,                                              &
       par_sum_exportdiag_200

      double precision,allocatable,dimension(:,:) ::                    &
       bil_fleuve,                                                      &
       ins_fleuve,                                                      &
       par_bil_fleuve

      double precision,allocatable,dimension(:) ::                      &
       bioriv_info

      double precision,allocatable,dimension(:,:) ::                    &
       bioriv_dt

      double precision,allocatable,dimension(:,:) ::                    &
       socard

      double precision,allocatable,dimension(:) ::                      &
       obc_bio_info

      double precision,allocatable,dimension(:,:) ::                    &
       obc_bio_dt

      double precision :: &
       bioatm_info     &
      ,dustpfr  &
      ,psolublefr

      double precision,allocatable,dimension(:) :: &
       bioatm_dt

      double precision,allocatable,dimension(:,:) ::                    &
       flux_ob_west_vert,                                               &
       flux_ob_south_vert,                                              &
       flux_ob_east_vert,                                               &
       sum_ob_west_vert,                                                &
       sum_ob_south_vert,                                               &
       sum_ob_east_vert,                                                &
       par_sum_ob_west_vert,                                            &
       par_sum_ob_south_vert,                                           &
       par_sum_ob_east_vert

      double precision,allocatable,dimension(:,:,:) ::                  &
       flux_ob_west_pt,                                                 &
       sum_ob_west_pt

      double precision,allocatable,dimension(:,:,:) ::                  &
       flux_ob_south_pt,                                                &
       sum_ob_south_pt

      character,dimension(:),allocatable :: &
       biobcfile*90                         &
      ,obc_bio_file*90                      &
      ,biovarname*10                        &
      ,dust_file*90                         &
      ,depnit_file*90                       &
      ,depammo_file*90                       &
      ,dustbinreclist*30                    &
      ,depnitbinreclist*30                  &
      ,depammobinreclist*30

      double precision,dimension(:),allocatable ::   &
      dustfile_nextime                                                 &
     ,dustfile_prvtime                                                 &
     ,depnitfile_nextime                                               &
     ,depnitfile_prvtime                                               &
     ,depammofile_nextime                                              &
     ,depammofile_prvtime


       double precision ::  &
       dust_lonmin                                                  &
      ,dust_latmin                                                  &
      ,dust_lonmax                                                  &
      ,dust_latmax                                                  &
      ,dust_resol                                                   &
      ,dust_resol_u                                                 &
      ,dust_resol_v                                                 &
      ,dust_lonstr                                                  &
      ,dust_lonend                                                  &
      ,dust_londlt                                                  &
      ,dust_latstr                                                  &
      ,dust_latend                                                  &
      ,dust_latdlt                                                  &
      ,depnit_lonmin                                                  &
      ,depnit_latmin                                                  &
      ,depnit_lonmax                                                  &
      ,depnit_latmax                                                  &
      ,depnit_resol                                                   &
      ,depnit_resol_u                                                 &
      ,depnit_resol_v                                                 &
      ,depnit_lonstr                                                  &
      ,depnit_lonend                                                  &
      ,depnit_londlt                                                  &
      ,depnit_latstr                                                  &
      ,depnit_latend                                                  &
      ,depnit_latdlt                                                &
      ,depammo_lonmin                                                  &
      ,depammo_latmin                                                  &
      ,depammo_lonmax                                                  &
      ,depammo_latmax                                                  &
      ,depammo_resol                                                   &
      ,depammo_resol_u                                                 &
      ,depammo_resol_v                                                 &
      ,depammo_lonstr                                                  &
      ,depammo_lonend                                                  &
      ,depammo_londlt                                                  &
      ,depammo_latstr                                                  &
      ,depammo_latend                                                  &
      ,depammo_latdlt

      double precision,dimension(:,:,:),allocatable ::              &
       dd_w                                                         &
      ,wdr_w                                                        &
      ,wdw_w                                                        &
      ,ddnit_w                                                         &
      ,wdrnit_w                                                        &
      ,wdwnit_w                                                        &
      ,ddammo_w &
      ,wdrammo_w &
      ,wdwammo_w

      double precision,dimension(:,:),allocatable ::              &
       AtmDepNit   &
      ,AtmDepAmmo   &
      ,AtmDepDOP   &
      ,AtmDepDON   &
      ,AtmDepDOC   &
      ,AtmDepPho

      character :: &
       bioatm_file*90


      double precision,dimension(2) ::                                  &
       nest_bio_out,                                                    &
       biobcafdt

      double precision,dimension(3) ::                                  &
       sumppbgdl_platcl,                                                &
       ppbgdl_platcl,                                                   &
       moyppbgdl_platcl,                                                &
       sumnpbgdl_platcl,                                                &
       npbgdl_platcl,                                                   &
       moynpbgdl_platcl,                                                &
       sumrpbgdl_platcl,                                                &
       rpbgdl_platcl,                                                   &
       moyrpbgdl_platcl

      double precision,dimension(4) ::                                  &
       respphyto,                                                       &
       netppb,                                                          &
       sumppbtotal,                                                     &
       sumppbtotal_plat,                                                &
       sumnpbtotal,                                                     &
       sumnetppbtotal,                                                  &
       sumrespphytototal,                                               &
       sumppbdyf,                                                       &
       sumppbgdl,                                                       &
       sumppbmedoc,                                                     &
       sumppblig,                                                       &
       sumppbdiag,                                                      &
       sumnpbdyf,                                                       &
       sumnpbmedoc,                                                     &
       sumnpblig,                                                       &
       sumnpbdiag,                                                      &
       sumnetppbdyf,                                                    &
       sumnetppbgdl,                                                    &
       sumnetppblig,                                                    &
       ppbtotal,                                                        &
       ppbtotal_plat,                                                   &
       npbtotal,                                                        &
       netppbtotal,                                                     &
       respphytototal,                                                  &
       ppbdyf,                                                          &
       ppbgdl,                                                          &
       ppbmedoc,                                                        &
       ppblig,                                                          &
       ppbdiag,                                                         &
       npbdyf,                                                          &
       npbmedoc,                                                        &
       npblig,                                                          &
       npbdiag,                                                         &
       netppbdyf,                                                       &
       netppbgdl,                                                       &
       netppblig,                                                       &
       moyppbtotal,                                                     &
       moyppbtotal_plat,                                                &
       moynpbtotal,                                                     &
       moyrespphytototal,                                               &
       moynetppblig,                                                    &
       moynetppbgdl,                                                    &
       moynetppbdyf,                                                    &
       moynetppbtotal,                                                  &
       moyppbdyf,                                                       &
       moyppbgdl,                                                       &
       moyppbmedoc,                                                     &
       moyppblig,                                                       &
       moyppbdiag,                                                      &
       moynpbdyf,                                                       &
       moynpbmedoc,                                                     &
       moynpblig,                                                       &
       moynpbdiag,                                                      &
       ppb,                                                             &
       npb

      double precision,dimension(5) ::                                  &
       sumrpbtotal,                                                     &
       rpbtotal,                                                        &
       moyrpbtotal,                                                     &
       sumrpblig,                                                       &
       rpblig,                                                          &
       moyrpblig,                                                       &
       sumrpbdyf,                                                       &
       rpbdyf,                                                          &
       moyrpbdyf,                                                       &
       sumrpbmedoc,                                                     &
       rpbmedoc,                                                        &
       moyrpbmedoc,                                                     &
       sumrpbdiag,                                                      &
       rpbdiag,                                                         &
       moyrpbdiag,                                                      &
       rpb


      double precision,dimension(8) ::                                  &
       exptotal,                                                        &
       sumexptotal

      double precision ::                                               &
       bionest_sum,                                                     &
       dti_bio,                                                         &
       ssr_trans,                                                       &
       eupho,                                                           &
       nitrif,                                                          &
       resptot,                                                         &
       exuctot,                                                         &
       sumnitriftotal,                                                  &
       sumresptottotal,                                                 &
       sumexuctottotal,                                                 &
       sumgrazctotal,                                                   &
       nitriftotal,                                                     &
       resptottotal,                                                    &
       exuctottotal,                                                    &
       moynitriftotal,                                                  &
       moyresptottotal,                                                 &
       moyexuctottotal,                                                 &
       moyexpdyf,                                                       &
       moyexpgdl,                                                       &
       expdyf,                                                          &
       expgdl,                                                          &
       sumexpdyf,                                                       &
       sumexpgdl,                                                       &
       sumareagdl,                                                      &
       sumareatotal,                                                    &
       sumarealig,                                                      &
       sumarea1000m,                                                    &
       sumareaemed,                                                     &
       sumareawmed,                                                     &
       expflux,                                                         &
       biobcinfo,                                                       &
       rap_biobc             =0.,                                       &
       dt_biobc,                                                        &
       datenudgbegin                                                    &
       ,datenudgend                                                     &
       ,nudgperiod                                                      &
       ,totalnitsurf                                                    &
       ,totalnitsurfdyf,                                                &
       grazctotal,                                                      &
       moygrazctotal,                                                   &
       sumareagdl_plat,                                                 &
       sumareamedoc,                                                    &
       sumareatotal_plat,                                               &
       sumareadyf,                                                      &
       sumareadiag,                                                     &
       grazzooc,                                                        &
       grazzoocpur,                                                     &
       grazzoocherbi,                                                   &
       gbac,                                                            &
       gbactotal,                                                       &
!       ExcHeteroNtop,                                                   &
!       ExcbactNH4top,                                                   &
!       ExcBactNH4N,                                                     &
!       ExcZooAmmoN,                                                     &
!       ExczooNH4top,                                                    &
!       NitrifCOLUMN,                                                    &
       sumgbactotal,                                                    &
       moygbactotal,                                                    &
       gbacdyf,                                                         &
       gbacgdl,                                                         &
       sumgbacdyf,                                                      &
       sumgbacgdl,                                                      &
       moygbacdyf,                                                      &
       moygbacgdl,                                                      &
       respzoototal,                                                    &
       sumrespzoototal,                                                 &
       moyrespzoototal,                                                 &
       respzoo,                                                         &
       rempoc,                                                          &
       rempon,                                                          &
       lostpoc,                                                         &
       lostpon,                                                         &
       uptbactdon,                                                      &
       excheteron,                                                      &
       sum_exportc_bot,                                                 &
       sum_exportn_bot,                                                 &
       sum_exportp_bot,                                                 &
       sum_exportsi_bot,                                                &
       exportc_bot,                                                     &
       exportn_bot,                                                     &
       exportp_bot,                                                     &
       exportsi_bot,                                                    &
       ppbgdl_plat,                                                     &
       rpbgdl_plat,                                                     &
       npbgdl_plat,                                                     &
       netppbgdl_plat,                                                  &
       moyppbgdl_plat,                                                  &
       moyrpbgdl_plat,                                                  &
       moynpbgdl_plat,                                                  &
       moynetppbgdl_plat,                                               &
       sumppbgdl_plat,                                                  &
       sumrpbgdl_plat,                                                  &
       sumnpbgdl_plat,                                                  &
       sumnetppbgdl_plat,                                               &
       SUMEXCHETERONGDL_plat,                                           &
       INSEXCHETERONGDL_plat,                                           &
       SUMUPTBACTDONGDL_plat,                                           &
       INSUPTBACTDONGDL_plat,                                           &
       SUMLOSTPOCGDL_plat,                                              &
       SUMLOSTPONGDL_plat,                                              &
       INSLOSTPOCGDL_plat,                                              &
       INSLOSTPONGDL_plat,                                              &
       SUMREMPOCGDL_plat,                                               &
       SUMREMPONGDL_plat,                                               &
       INSPPBGDL_plat,                                                  &
       gbacgdl_plat,                                                    &
       insgbacgdl_plat,                                                 &
       SUMRESPPHYTOGDL_plat,                                            &
       SUMRESPZOOGDL_plat,                                              &
       INSRESPPHYTOGDL_plat,                                            &
       INSRESPZOOGDL_PLAT,                                              &
       INSNETPPBGDL_plat,                                               &
       INSNPBGDL_plat,                                                  &
       INSRPBGDL_plat,                                                  &
       SUMGBACGDL_plat,                                                 &
       INSRESPBACTGDL_plat,                                             &
       SUMRESPTOTGDL_plat,                                              &
       INSRESPTOTGDL_plat,                                              &
       BIL_TOUTFLEUVEPOC,                                               &
       BIL_TOUTFLEUVEDOC,                                               &
       BIL_TOUTFLEUVEPON,                                               &
       BIL_TOUTFLEUVEDON,                                               &
       BIL_TOUTFLEUVENO3,                                               &
       BIL_TOUTFLEUVENH3,                                               &
       INS_TOUTFLEUVEPOC,                                               &
       INS_TOUTFLEUVEDOC,                                               &
       INS_TOUTFLEUVEPON,                                               &
       INS_TOUTFLEUVEDON,                                               &
       INS_TOUTFLEUVENO3,                                               &
       INS_TOUTFLEUVENH3,                                               &
       SUM_BIOPOC_P,                                                    &
       SUM_BIODOC_P,                                                    &
       SUM_BIOPON_P,                                                    &
       SUM_BIODON_P,                                                    &
       SUM_BIONO3_P,                                                    &
       SUM_BIONH3_P,                                                    &
       SUM_HBIO_P,                                                      &
       SUM_HBION_P,                                                     &
       SUM_BIOPOC_P1,                                                   &
       SUM_BIODOC_P1,                                                   &
       SUM_BIOPON_P1,                                                   &
       SUM_BIODON_P1,                                                   &
       SUM_BIONO3_P1,                                                   &
       SUM_BIONH3_P1,                                                   &
       SUM_HBIO_P1,                                                     &
       SUM_HBION_P1,                                                    &
       bil_ppb,                                                         &
       bil_npb,                                                         &
       bil_rpb,                                                         &
       bil_netppb,                                                      &
       bil_respphyto,                                                   &
       bil_rempoc,                                                      &
       bil_rempon,                                                      &
       bil_gbac,                                                        &
       bil_lostpoc,                                                     &
       bil_lostpon,                                                     &
       bil_uptbactdon,                                                  &
       bil_excheteron,                                                  &
       bil_respzoo,                                                     &
       bil_respbact,                                                    &
       ins_ppb,                                                         &
       ins_npb,                                                         &
       ins_rpb,                                                         &
       ins_respphyto,                                                   &
       ins_gbac,                                                        &
       ins_lostpoc,                                                     &
       ins_lostpon,                                                     &
       ins_uptbactdon,                                                  &
       ins_excheteron,                                                  &
       ins_respzoo,                                                     &
       ins_respbact,                                                    &
       ins_netppb,                                                      &
       SUM_SEDC,                                                        &
       SUM_SEDN,                                                        &
       SUM_SEDP,                                                        &
       SUM_SEDSI,                                                       &
       SEDC,                                                            &
       SEDN,                                                            &
       SEDP,                                                            &
       SEDSI,                                                           &
       par_sum_exportc_bot,                                             &
       par_sum_exportn_bot,                                             &
       par_sum_exportp_bot,                                             &
       par_sum_exportsi_bot

      double precision ::                                               &
       DecayRate,                                                       &
       cDet20SED,                                                       &
       kFast,                                                           &
       kSlow,                                                           &
       pFast,                                                           &
       cFDet20SED,                                                      &
       cSDet20SED,                                                      &
       BQ10,                                                            &
       NCrFD,                                                           &
       NCrSD,                                                           &
       pNit,                                                            &
       pdeNit,                                                          &
       pPMin,                                                           &
       pSiMin,                                                          &
       O2bw,                                                            &
       NitrifSMean,                                                     &
       DeNitrifSMean,                                                   &
       NO3FluxSMean,                                                    &
       NH3FluxSMean,                                                    &
       NMinSMean,                                                       &
       CMinSMean,                                                       &
       NDepoSMean,                                                      &
       CDepoSMean,                                                      &
       AnoxMinSMean,                                                    &
       NitrifSMean_est,                                                &
       DeNitrifSMean_est,                                              &
       NO3FluxSMean_est,                                               &
       NH3FluxSMean_est,                                               &
       NMinSMean_est,                                                  &
       CMinSMean_est,                                                  &
       NDepoSMean_est,                                                 &
       CDepoSMean_est,                                                 &
       AnoxMinSMean_est,                                               &
       SUM_NITRIF,                                                      &
       SUM_DENITRIF,                                                    &
       SUM_NO3Flux,                                                     &
       SUM_NH3Flux,                                                     &
       SUM_NMin,                                                        &
       SUM_CMin,                                                        &
       SUM_NDepo,                                                       &
       SUM_CDepo,                                                       &
       SUM_AnoxMin,                                                     &
       SUM_NITRIF_est,                                                 &
       SUM_DENITRIF_est,                                               &
       SUM_NO3Flux_est,                                                &
       SUM_NH3Flux_est,                                                &
       SUM_NMin_est,                                                   &
       SUM_CMin_est,                                                   &
       SUM_NDepo_est,                                                  &
       SUM_CDepo_est,                                                  &
       SUM_AnoxMin_est
!______________________________________________________________________
! LES NOMBRES FORCEMENT EN SIMPLE PRECISION:

      real*4,allocatable,dimension(:,:,:) ::                            &
       river_bio
      real*4,allocatable,dimension(:,:,:) ::                            &
       fluxmoyen_w
      real*4                                                            &
       count_output
      real*4,allocatable,dimension(:,:,:) ::                            &
       obc_bio_data
      real*4,allocatable,dimension(:) :: &
       atm_bio
!______________________________________________________________________
! LES NOMBRES ENTIERS:

      integer,allocatable,dimension(:,:) ::                             &
       bioriv_date

      integer,allocatable,dimension(:) ::                               &
       biobc_type

      integer,allocatable,dimension(:) ::                               &
       ild_to_sd

      integer,dimension(:),allocatable ::                              &  
      dustfile_nextrec                                                 &
     ,dustfile_prvtrec                                                 & 
     ,depnitfile_nextrec                                               &
     ,depnitfile_prvtrec                                               &
     ,depammofile_nextrec &
     ,depammofile_prvtrec

      integer,allocatable,dimension(:,:) ::                             &
       obc_bio_date

      integer ::                                                        &
        vb                                                              &
       ,vbmax                                                           &
       ,vbmax_eco3ms   =0                                               &
       ,kso                                                             &
       ,ksomax                                                          &
       ,iadvec_bio                                                      &
       ,i1dbio                                                          &
       ,iptbio                                                          &
       ,jptbio                                                          &
       ,mbio1                                                           &
       ,mbio2                                                           &
       ,nbio1                                                           &
       ,nbio2                                                           &
       ,mbio1_glob                                                      &
       ,mbio2_glob                                                      &
       ,nbio1_glob                                                      &
       ,nbio2_glob                                                      &
       ,tps_ppb                                                         &
       ,tps_ppb_2d                                                      &
       ,tps_strada                                                      &
       ,tps_strada_2d                                                   &
       ,tps_sed                                                         &
       ,tps_benth_2d                                                    &
       ,jgdl                                                            &
       ,idyf                                                            &
       ,jdyf                                                            &
       ,ilig                                                            &
       ,inudging                                                        &
       ,knudgbegin                                                      &
       ,knudgend                                                        &
       ,interp_lr                                                       &
       ,ibiobc_af                                                       &
       ,idyf1                                                           &
       ,jdyf1                                                           &
       ,idyf2                                                           &
       ,jdyf2                                                           &
       ,idyf3                                                           &
       ,jdyf3                                                           &
       ,igdl1                                                           &
       ,jgdl1                                                           &
       ,igdl2                                                           &
       ,jgdl2                                                           &
       ,imedoc                                                          &
       ,jmedoc                                                          &
       ,ivill                                                           &
       ,jvill                                                           &
       ,imars                                                           &
       ,jmars                                                           &
       ,iban                                                            &
       ,jban                                                            &
       ,icn1                                                            &
       ,jcn1                                                            &
       ,icn2                                                            &
       ,jcn2                                                            &
       ,i1domain                                                        &
       ,j1domain                                                        &
       ,i2domain                                                        &
       ,j2domain                                                        &
       ,i1gdl                                                           &
       ,j1gdl                                                           &
       ,i2gdl                                                           &
       ,j2gdl                                                           &
       ,i1medoc                                                         &
       ,j1medoc                                                         &
       ,i2medoc                                                         &
       ,j2medoc                                                         &
       ,i1ligure                                                        &
       ,j1ligure                                                        &
       ,i2ligure                                                        &
       ,j2ligure                                                        &
       ,i1diag                                                          &
       ,j1diag                                                          &
       ,i2diag                                                          &
       ,j2diag                                                          &
       ,iobwest                                                         &
       ,iobeast                                                         &
       ,jobsouth                                                        &
       ,jobnorth                                                        &
       ,gridbio                                                         &
       ,IBenthic                                                        &
       ,NumBDET                                                         &
       ,BTfunc                                                          &
       ,TPS_BENT                                                        &
       ,SUMT_BENT                                                       &
       ,icas1                                                           &
       ,jcas1                                                           &
       ,icas2                                                           &
       ,jcas2                                                           &
       ,icas3                                                           &
       ,jcas3                                                           &
       ,icas4                                                           &
       ,jcas4                                                           &
       ,icas5                                                           &
       ,jcas5                                                           &
       ,icas6                                                           &
       ,jcas6                                                           &
       ,icas7                                                           &
       ,jcas7                                                           &
       ,icas8                                                           &
       ,jcas8                                                           &
       ,icas9                                                           &
       ,jcas9                                                           &
       ,icas10                                                          &
       ,jcas10                                                          &
       ,icas11                                                          &
       ,jcas11                                                          &
       ,icas12                                                          &
       ,jcas12                                                          &
       ,icas13                                                          &
       ,jcas13                                                          &
       ,icas14                                                          &
       ,jcas14                                                          &
       ,icas15                                                          &
       ,jcas15                                                          &
       ,icas16                                                          &
       ,jcas16                                                          &
       ,icas17                                                          &
       ,jcas17                                                          &
       ,icas18                                                          &
       ,jcas18 &
       ,idust &
       ,ichoicedust  &
       ,idepnit                                                        &
       ,idepammo &
       ,ndust                                                       &
       ,ndepnit                                                       &
       ,ndepammo                                                       &
       ,flag_dust_average                                           &
       ,flag_depnit_average                                         &
       ,flag_depammo_average                                         &
          ,dust_imax                                               &
          ,dust_jmax                                               &
          ,dust_kmax                                               &
          ,dustzoom_istr                                           &
          ,dustzoom_iend                                           &
          ,dustzoom_jstr                                           &
          ,dustzoom_jend                                           &
          ,dustfull_imax                                           &
          ,dustfull_jmax                                           &
          ,depnit_imax                                               &
          ,depnit_jmax                                               &
          ,depnit_kmax                                               &
          ,depnitzoom_istr                                           &
          ,depnitzoom_iend                                           &
          ,depnitzoom_jstr                                           &
          ,depnitzoom_jend                                           &
          ,depnitfull_imax                                           &
          ,depnitfull_jmax                                           &
          ,depammo_imax                                               &
          ,depammo_jmax                                               &
          ,depammo_kmax                                               &
          ,depammozoom_istr                                           &
          ,depammozoom_iend                                           &
          ,depammozoom_jstr                                           &
          ,depammozoom_jend                                           &
          ,depammofull_imax                                           &
          ,depammofull_jmax                                                          

      integer ::          &
           dd_id    =0   & ! norestart
          ,wdr_id     =0   & ! norestart
          ,wdw_id   =0    ! norestart                                                          

      integer ::          &
          max_dust_time_counter      &
         ,max_depnit_time_counter   &
         ,max_depammo_time_counter


       integer,dimension(23) ::                                         &
        imoose                                                          &
       ,jmoose                                                          


      integer,dimension(6) ::                                           &
        datebiobc

contains

      subroutine biology_allocate
      implicit none
      
      allocate(bio_t        (-1:imax+2,-1:jmax+2,0:kmax+1,vbmax)) ;  bio_t=0.
      allocate(biosum_t     (0:imax+1,0:jmax+1,kmax,vbmax))       ;  biosum_t=0.
      allocate(tendancebio_t(0:imax+1,0:jmax+1,kmax,vbmax))       ;  tendancebio_t=0.
      allocate(fluxbio_w    (0:imax+1,0:jmax+1,vbmax,2))          ;   fluxbio_w=0.

      allocate(nitobc_t(0:imax+1,0:jmax+1,kmax,0:2),      &
               ammobc_t(0:imax+1,0:jmax+1,kmax,0:2),      &
               phoobc_t(0:imax+1,0:jmax+1,kmax,0:2),      &
               silobc_t(0:imax+1,0:jmax+1,kmax,0:2))
      nitobc_t=0.
      ammobc_t=0.
      phoobc_t=0.
      silobc_t=0.

      allocate(biobc_i_t(0:imax+1,kmax,vbmax,2,2)) ; biobc_i_t=0.
      allocate(biobc_j_t(0:jmax+1,kmax,vbmax,2,2)) ; biobc_j_t=0.

!      allocate(nitrif3d(0:imax+1,0:jmax+1,kmax)                    &
!              ,uptnit3d(0:imax+1,0:jmax+1,kmax))
!      nitrif3d=0
!      uptnit3d=0   

!      allocate(nitrif2d(0:imax+1,0:jmax+1)                    &
      allocate(resp2d(0:imax+1,0:jmax+1)                      &
                ,graz2d(0:imax+1,0:jmax+1)                    &
                ,gbac2d(0:imax+1,0:jmax+1)                    &
!                ,NITRIFCOLUMN2D(0:imax+1,0:jmax+1)            &
!                ,EXCHETERON2D(0:imax+1,0:jmax+1)              &
!                ,ExcHeteroNtop2D(0:imax+1,0:jmax+1)           &
!                ,ExcbactNH4top2D(0:imax+1,0:jmax+1)           &
!                ,ExcBactNH4N2D(0:imax+1,0:jmax+1)             &
!                ,ExcZooAmmoN2D(0:imax+1,0:jmax+1)             &
!                ,ExczooNH4top2D(0:imax+1,0:jmax+1)            &
             ,exuctot2d(0:imax+1,0:jmax+1)                    &
                ,grazpur2d(0:imax+1,0:jmax+1)                 &
                ,grazherbi2d(0:imax+1,0:jmax+1)               &
        ,sumnitratesurf(0:imax+1,0:jmax+1)                    &
             ,ppbsurf2d(0:imax+1,0:jmax+1)                    &
              ,ppbmax2d(0:imax+1,0:jmax+1)                    &
          ,profppbmax2d(0:imax+1,0:jmax+1)                    &
      ,sum_exportc_bot_2d(0:imax+1,0:jmax+1)                  &
      ,sum_exportn_bot_2d(0:imax+1,0:jmax+1)                  &
      ,sum_exportp_bot_2d(0:imax+1,0:jmax+1)                  &
      ,sum_exportsi_bot_2d(0:imax+1,0:jmax+1)                 &
      ,NO3efflux2d(0:imax+1,0:jmax+1)                           &
      ,NH4efflux2d(0:imax+1,0:jmax+1)                           &
      ,Pefflux2d(0:imax+1,0:jmax+1)                           &
      ,Siefflux2d(0:imax+1,0:jmax+1)                          &
      ,ExcbactNH4NTOPLAYER2D(0:imax+1,0:jmax+1)                                             &
      ,ExcZooAmmoNTOPLAYER2D(0:imax+1,0:jmax+1)                                             &
      ,ExcZooPO4PTOPLAYER2D(0:imax+1,0:jmax+1)                                              &
      ,ExcbactPO4PTOPLAYER2D(0:imax+1,0:jmax+1)                                             &
      ,ExuSiTOPLAYER2D(0:imax+1,0:jmax+1)                                                   &
      ,UptNitTOPLAYER2D(0:imax+1,0:jmax+1)                                                  &
      ,UptAmmoTOPLAYER2D(0:imax+1,0:jmax+1)                                                 &
      ,UptPTOPLAYER2D(0:imax+1,0:jmax+1)                                                    &
      ,UptSiTOPLAYER2D(0:imax+1,0:jmax+1)                                                   &
      ,NitrifTOPLAYER2D(0:imax+1,0:jmax+1)                                                  &
      ,RemSMOPSiTOPLAYER2D(0:imax+1,0:jmax+1)                                               &
      ,RemLMOPSiTOPLAYER2D(0:imax+1,0:jmax+1)                                               &
      ,ExcbactNH4NINT2D(0:imax+1,0:jmax+1)                                             &
      ,ExcZooAmmoNINT2D(0:imax+1,0:jmax+1)                                             &
      ,ExcZooPO4PINT2D(0:imax+1,0:jmax+1)                                              &
      ,ExcbactPO4PINT2D(0:imax+1,0:jmax+1)                                             &
      ,ExuSiINT2D(0:imax+1,0:jmax+1)                                                   &
      ,UptNitINT2D(0:imax+1,0:jmax+1)                                                  &
      ,UptAmmoINT2D(0:imax+1,0:jmax+1)                                                 &
      ,UptPINT2D(0:imax+1,0:jmax+1)                                                    &
      ,UptSiINT2D(0:imax+1,0:jmax+1)                                                   &
      ,NitrifINT2D(0:imax+1,0:jmax+1)                                                  &
      ,RemSMOPSiINT2D(0:imax+1,0:jmax+1)                                               &
      ,RemLMOPSiINT2D(0:imax+1,0:jmax+1)                                               &
      ,ExcbactNH4NDEEP2D(0:imax+1,0:jmax+1)                                             &
      ,ExcZooAmmoNDEEP2D(0:imax+1,0:jmax+1)                                             &
      ,ExcZooPO4PDEEP2D(0:imax+1,0:jmax+1)                                              &
      ,ExcbactPO4PDEEP2D(0:imax+1,0:jmax+1)                                             &
      ,ExuSiDEEP2D(0:imax+1,0:jmax+1)                                                   &
      ,UptNitDEEP2D(0:imax+1,0:jmax+1)                                                  &
      ,UptAmmoDEEP2D(0:imax+1,0:jmax+1)                                                 &
      ,UptPDEEP2D(0:imax+1,0:jmax+1)                                                    &
      ,UptSiDEEP2D(0:imax+1,0:jmax+1)                                                   &
      ,NitrifDEEP2D(0:imax+1,0:jmax+1)                                                  &
      ,RemSMOPSiDEEP2D(0:imax+1,0:jmax+1)                                               &
      ,RemLMOPSiDEEP2D(0:imax+1,0:jmax+1))
!      ,netppb2d(0:imax+1,0:jmax+1))

!      nitrif2d=0.
      resp2d=0.
      graz2d=0.
      grazpur2d=0.
      grazherbi2d=0.
      gbac2d=0.
!      NITRIFCOLUMN2D=0.    
!      EXCHETERON2D=0.   
!      ExcHeteroNtop2D=0.  
!      ExcbactNH4top2D=0.  
!      ExcBactNH4N2D=0.  
!      ExcZooAmmoN2D=0. 
!      ExczooNH4top2D=0.
      exuctot2d=0.
      sumnitratesurf=0.
      ppbsurf2d=0.
      ppbmax2d=0.
      profppbmax2d=0.
      sum_exportc_bot_2d=0.
      sum_exportn_bot_2d=0.
      sum_exportp_bot_2d=0.
      sum_exportsi_bot_2d=0.
      NO3efflux2d=0.
      NH4efflux2d=0.
      Pefflux2d=0.
      Siefflux2d=0.
!      netppb2d=0.
      ExcbactNH4NTOPLAYER2D=0.  
      ExcZooAmmoNTOPLAYER2D=0.
      ExcZooPO4PTOPLAYER2D=0.
      ExcbactPO4PTOPLAYER2D=0.  
      ExuSiTOPLAYER2D=0.        
      UptNitTOPLAYER2D=0.       
      UptAmmoTOPLAYER2D=0.      
      UptPTOPLAYER2D=0.         
      UptSiTOPLAYER2D=0.     
      NitrifTOPLAYER2D=0.     
      RemSMOPSiTOPLAYER2D=0.    
      RemLMOPSiTOPLAYER2D=0.    
      ExcbactNH4NINT2D=0.       
      ExcZooAmmoNINT2D=0.       
      ExcZooPO4PINT2D=0.        
      ExcbactPO4PINT2D=0.       
      ExuSiINT2D=0.             
      UptNitINT2D=0.            
      UptAmmoINT2D=0.           
      UptPINT2D=0.              
      UptSiINT2D=0.             
      NitrifINT2D=0.            
      RemSMOPSiINT2D=0.         
      RemLMOPSiINT2D=0.         
      ExcbactNH4NDEEP2D=0.      
      ExcZooAmmoNDEEP2D=0.      
      ExcZooPO4PDEEP2D=0.       
      ExcbactPO4PDEEP2D=0.      
      ExuSiDEEP2D=0.            
      UptNitDEEP2D=0.           
      UptAmmoDEEP2D=0.          
      UptPDEEP2D=0.             
      UptSiDEEP2D=0.            
      NitrifDEEP2D=0.           
      RemSMOPSiDEEP2D=0.        

      allocate(CBDet(imax,jmax),                                        &
       CBFDet(imax,jmax),                                               &
       CBSDet(imax,jmax),                                               &
       NBDet(imax,jmax),                                                &
       PBDet(imax,jmax),                                                &
       SiBDet(imax,jmax),                                               &
       dCBDet(imax,jmax),                                               &
       dCBFDet(imax,jmax),                                              &
       dCBSDet(imax,jmax),                                              &
       dNBDet(imax,jmax),                                               &
       dPBDet(imax,jmax),                                               &
       dSiBDet(imax,jmax),                                              &
       CDepo(imax,jmax),                                                &
       NDepo(imax,jmax),                                                &
       PDepo(imax,jmax),                                                &
       SiDepo(imax,jmax),                                               &
       NitrificationB(imax,jmax),                                       &
       DenitrificationB(imax,jmax),                                     &
       O2Min(imax,jmax),                                                &
       AnoxMin(imax,jmax),                                              &
       O2ODU(imax,jmax),                                                &
       SUMT_NITRIF_est(imax,jmax),                                          &
       SUMT_DENITRIF_est(imax,jmax),                                        &
       SUMT_NO3Flux_est(imax,jmax),                                         &
       SUMT_NH3Flux_est(imax,jmax),                                         &
       SUMT_NMin_est(imax,jmax),                                            &
       SUMT_CMin_est(imax,jmax),                                            &
       SUMT_NDepo_est(imax,jmax),                                           &
       SUMT_CDepo_est(imax,jmax),                                           &
       SUMT_NITRIF(imax,jmax),                                          &
       SUMT_DENITRIF(imax,jmax),                                        &
       SUMT_NO3Flux(imax,jmax),                                         &
       SUMT_NH3Flux(imax,jmax),                                         &
       SUMT_NMin(imax,jmax),                                            &
       SUMT_CMin(imax,jmax),                                            &
       SUMT_NDepo(imax,jmax),                                           &
       SUMT_CDepo(imax,jmax),                                           &
       CMin_out(imax,jmax),                                             &
       NMin_out(imax,jmax),                                             &
       SUMT_O2Min(imax,jmax),                                           &
       SUMT_AnoxMin(imax,jmax),                                          &
       SUMT_AnoxMin_est(imax,jmax))
       CBDet=0.
       CBFDet=0.
       CBSDet=0.
       NBDet=0.
       PBDet=0.
       SiBDet=0.
       dCBDet=0.
       dCBFDet=0.
       dCBSDet=0.
       dNBDet=0.
       dPBDet=0.
       dSiBDet=0.
       CDepo=0.
       NDepo=0.
       PDepo=0.
       SiDepo=0.
       NitrificationB=0.
       DenitrificationB=0.
       O2Min=0.
       AnoxMin=0.
       O2ODU=0.
       SUMT_NITRIF_est=0.
       SUMT_DENITRIF_est=0.
       SUMT_NO3Flux_est=0.
       SUMT_NH3Flux_est=0.
       SUMT_NMin_est=0.
       SUMT_CMin_est=0.
       SUMT_NDepo_est=0.
       SUMT_CDepo_est=0.
       SUMT_NITRIF=0.
       SUMT_DENITRIF=0.
       SUMT_NO3Flux=0.
       SUMT_NH3Flux=0.
       SUMT_NMin=0.
       SUMT_CMin=0.
       SUMT_NDepo=0.
       SUMT_CDepo=0.
       CMin_out=0.
       NMin_out=0.
       SUMT_O2Min=0.
       SUMT_AnoxMin=0.
       SUMT_AnoxMin_est=0.

      allocate(ppb2d(0:imax+1,0:jmax+1,4),             &
               netppb2d(0:imax+1,0:jmax+1,4),             &
               npb2d(0:imax+1,0:jmax+1,4),             &
               chl2d(0:imax+1,0:jmax+1,4))
      ppb2d=0.
      netppb2d=0.
      npb2d=0.
      chl2d=0.

      allocate(rpb2d(0:imax+1,0:jmax+1,5))
      rpb2d=0.
      allocate(exp2d(0:imax+1,0:jmax+1,40))
      exp2d=0.

      allocate(bio_trans(kmax,vbmax),         &
            tdcbio_trans(kmax,vbmax))
      bio_trans=0.
      tdcbio_trans=0.

      allocate(pwk_trans(kmax+1))
      pwk_trans=0.

      allocate(tem_trans(kmax),          &
               p3d_trans(kmax))
      tem_trans=0.
      p3d_trans=0.

      allocate(radionucleide(vbmax))
      radionucleide=0

!     allocate(par_bio(imax,jmax,kmax))
!     par_bio=0.

      allocate(wsed(1:kmax+1,vbmax)) ; wsed=0. !12-02-17

      allocate(                                                   &
                       tauR(vbmax),                               &
                alpha_radio(vbmax),                               &
                 radio_coef(vbmax),                               &
                export_plat(vbmax),                               &
            sum_export_plat(vbmax),                               &
             exportgdl_plat(vbmax),                               &
         sum_exportgdl_plat(vbmax),                               &
           export_200_1000m(vbmax),                               &
       sum_export_200_1000m(vbmax),                               &
                bil_fleuves(vbmax),                               &
                bil_fleuves_est(vbmax),                           &
                sum_ob_west(vbmax),                               &
               sum_ob_south(vbmax),                               &
              flux_ob_south(vbmax),                               &
               flux_ob_west(vbmax),                               &
       exportdyf_bot(vbmax),                                      &
       sum_exportdyf_bot(vbmax),                                  &
       ins_fleuves(vbmax),                                        &
       ins_fleuves_est(vbmax),                                        &
       sum_ob_east(vbmax),                                        &
       sum_ob_north(vbmax),                                       &
       flux_ob_east(vbmax),                                       &     
       flux_ob_north(vbmax),                                      &
       sumbio_p(vbmax),                                           &
       sumbion_p(vbmax),                                          &
       par_sum_ob_west(vbmax),                                    &
       par_sum_ob_south(vbmax),                                   &
       par_sum_ob_east(vbmax),                                    &
       par_sum_export_plat(vbmax),                                &
       par_sum_exportgdl_plat(vbmax),                             &
       par_sum_exportdyf_bot(vbmax),                              &
       par_sum_export_200_1000m(vbmax),                           &
       par_bil_fleuves(vbmax),                                    &
       exportdiag_bot(vbmax),                                     &
       sum_exportdiag_bot(vbmax),                                 &
       par_sum_exportdiag_bot(vbmax),                             &
       bil_fleuves_area(vbmax),                                  &
       bil_fleuves_area_est(vbmax))
      export_plat=0.
      sum_export_plat=0.
      exportgdl_plat=0.
      sum_exportgdl_plat=0.
      export_200_1000m=0.
      sum_export_200_1000m=0.
      bil_fleuves=0.
      bil_fleuves_est=0.
      sum_ob_west=0.
      sum_ob_south=0.
      flux_ob_west=0.
      exportdyf_bot=0.
      sum_exportdyf_bot=0.
      ins_fleuves=0.
      ins_fleuves_est=0.
       sum_ob_east=0.
       sum_ob_north=0.
       flux_ob_east=0.
       flux_ob_north=0.
       sumbio_p=0.
       sumbion_p=0.
       par_sum_ob_west=0.
       par_sum_ob_south=0.
       par_sum_ob_east=0.
       par_sum_export_plat=0.
       par_sum_exportgdl_plat=0.
       par_sum_exportdyf_bot=0.
       par_sum_export_200_1000m=0.
       par_bil_fleuves=0.
       exportdiag_bot=0.
       sum_exportdiag_bot=0.
       par_sum_exportdiag_bot=0.
       bil_fleuves_area=0.
       bil_fleuves_area_est=0.

       allocate(exportdyf_100(vbmax*2),                             &
       sum_exportdyf_100(vbmax*2),                                  &
       par_sum_exportdyf_100(vbmax*2))
       exportdyf_100=0.
       sum_exportdyf_100=0.
       par_sum_exportdyf_100=0.


      allocate(export_200(vbmax*3),                                  &
           sum_export_200(vbmax*3),                                  &
            exportgdl_200(vbmax*3),                                  &
        sum_exportgdl_200(vbmax*3),                                  &
            exportdyf_200(vbmax*3),                                  &
        sum_exportdyf_200(vbmax*3),                                  &
           exportdyf3_200(vbmax*3),                                  &
       sum_exportdyf3_200(vbmax*3),                                  &
       exportdyf_1000(vbmax*3),                                      &
       sum_exportdyf_1000(vbmax*3),                                  &
       par_sum_export_200(vbmax*3),                                  &
       par_sum_exportgdl_200(vbmax*3),                               &
       par_sum_exportdyf_200(vbmax*3),                               &
       par_sum_exportdyf3_200(vbmax*3),                              &
       par_sum_exportdyf_1000(vbmax*3),                              &
       exportdiag_200(vbmax*3),                                      &
       sum_exportdiag_200(vbmax*3),                                  &
       par_sum_exportdiag_200(vbmax*3))
        export_200=0.
        sum_export_200=0.
        exportgdl_200=0.
        sum_exportgdl_200=0.
        exportdyf_200=0.
        sum_exportdyf_200=0.
       exportdyf3_200=0.
       sum_exportdyf3_200=0.
       exportdyf_1000=0.
       sum_exportdyf_1000=0.
       par_sum_export_200=0.
       par_sum_exportgdl_200=0.
       par_sum_exportdyf_200=0.
       par_sum_exportdyf3_200=0.
       par_sum_exportdyf_1000=0.
       exportdiag_200=0.
       sum_exportdiag_200=0.
       par_sum_exportdiag_200=0.

      allocate(bil_fleuve(dim_river,vbmax),                           &
               ins_fleuve(dim_river,vbmax),                           &
               par_bil_fleuve(dim_river,vbmax))
      bil_fleuve=0.
      ins_fleuve=0.
      par_bil_fleuve=0.

      allocate(bioriv_info(dim_river))
      bioriv_info=0.

      allocate(bioriv_dt(dim_river,2))
      bioriv_dt=0.

      allocate(obc_bio_info(4))
      obc_bio_info=0.

      allocate(obc_bio_dt(6,2))
      obc_bio_dt=0.

      allocate(flux_ob_west_vert(kmax,vbmax),                     &
              flux_ob_south_vert(kmax,vbmax),                     &
              flux_ob_east_vert(kmax,vbmax),                      &
                sum_ob_west_vert(kmax,vbmax),                     &
               sum_ob_south_vert(kmax,vbmax),                     &    
               sum_ob_east_vert(kmax,vbmax),                     &
               par_sum_ob_west_vert(kmax,vbmax),                  &
              par_sum_ob_south_vert(kmax,vbmax),                  &
              par_sum_ob_east_vert(kmax,vbmax))
      flux_ob_west_vert=0.
      flux_ob_south_vert=0.
      flux_ob_east_vert=0.
      sum_ob_west_vert=0.
      sum_ob_south_vert=0.
      sum_ob_east_vert=0.
      par_sum_ob_west_vert=0.
      par_sum_ob_south_vert=0.
      par_sum_ob_east_vert=0.


      allocate(flux_ob_west_pt(0:jmax+1,kmax,vbmax),              &
                sum_ob_west_pt(0:jmax+1,kmax,vbmax))
      flux_ob_west_pt=0.
      sum_ob_west_pt=0.
      

      allocate(flux_ob_south_pt(0:imax+1,kmax,vbmax),             &
                sum_ob_south_pt(0:imax+1,kmax,vbmax))
      flux_ob_south_pt=0.
      sum_ob_south_pt=0.

      allocate(biobcfile(dim_river))           ; biobcfile='s'
      allocate(obc_bio_file(4))           ; biobcfile='s'
      allocate(biovarname(vbmax))              ; biovarname='s'
      bioatm_file='s'
      allocate(dust_file(3))              ;dust_file='s'
      allocate(depnit_file(3))              ;depnit_file='s'
      allocate(depammo_file(3))              ;depammo_file='s'
      allocate(dustbinreclist(3))           ;dustbinreclist='reset'
      allocate(depnitbinreclist(3))           ;depnitbinreclist='reset'
      allocate(depammobinreclist(3))     ;depammobinreclist='reset'
      allocate(dustfile_prvtime(3)) ;dustfile_prvtime=0.
      allocate(dustfile_nextime(3)) ;dustfile_nextime=0.
      allocate(dustfile_prvtrec(3)) ;dustfile_prvtrec=0
      allocate(dustfile_nextrec(3)) ;dustfile_nextrec=0
      allocate(depnitfile_prvtime(3)) ;depnitfile_prvtime=0.
      allocate(depnitfile_nextime(3)) ;depnitfile_nextime=0.
      allocate(depnitfile_prvtrec(3)) ;depnitfile_prvtrec=0
      allocate(depnitfile_nextrec(3)) ;depnitfile_nextrec=0
      allocate(depammofile_prvtime(3)) ;depammofile_prvtime=0.
      allocate(depammofile_nextime(3)) ;depammofile_nextime=0.
      allocate(depammofile_prvtrec(3)) ;depammofile_prvtrec=0
      allocate(depammofile_nextrec(3)) ;depammofile_nextrec=0


      allocate(dd_w               (0:imax+1,0:jmax+1,0:2)) ; dd_w=0
      allocate(wdr_w               (0:imax+1,0:jmax+1,0:2)) ; wdr_w=0
      allocate(wdw_w               (0:imax+1,0:jmax+1,0:2)) ; wdw_w=0
      allocate(ddnit_w               (0:imax+1,0:jmax+1,0:2)) ; ddnit_w=0
      allocate(wdrnit_w               (0:imax+1,0:jmax+1,0:2)) ; wdrnit_w=0
      allocate(wdwnit_w               (0:imax+1,0:jmax+1,0:2)) ; wdwnit_w=0
      allocate(ddammo_w               (0:imax+1,0:jmax+1,0:2)) ; ddammo_w=0
      allocate(wdrammo_w               (0:imax+1,0:jmax+1,0:2)) ; wdrammo_w=0
      allocate(wdwammo_w               (0:imax+1,0:jmax+1,0:2)) ; wdwammo_w=0

      allocate(AtmDepNit               (0:imax+1,0:jmax+1)) ; AtmDepNit=0
      allocate(AtmDepPho               (0:imax+1,0:jmax+1)) ; AtmDepPho=0
      allocate(AtmDepAmmo               (0:imax+1,0:jmax+1)) ; AtmDepAmmo=0
      allocate(AtmDepDOP               (0:imax+1,0:jmax+1)) ; AtmDepDOP=0
      allocate(AtmDepDON               (0:imax+1,0:jmax+1)) ; AtmDepDON=0
      allocate(AtmDepDOC               (0:imax+1,0:jmax+1)) ; AtmDepDOC=0
      allocate(river_bio(vbmax,dim_river,0:2))
      river_bio=0.
      allocate(atm_bio(0:2))
      atm_bio=0.
      allocate(fluxmoyen_w(imax,jmax,25))
      fluxmoyen_w=0.
      allocate(bioriv_date(6,dim_river))
      bioriv_date=0.
      allocate(biobc_type(vbmax))
      biobc_type=0.
      allocate(ild_to_sd(dim_river))
      ild_to_sd=0.
      allocate(obc_bio_data(kmax,4,0:2))
      obc_bio_data=0.
      allocate(obc_bio_date(6,4))
      obc_bio_date=0.

      end subroutine biology_allocate

!....................................................................

      subroutine biology_allocate_socard
      implicit none

       allocate(socard(9,ksomax))

      end subroutine biology_allocate_socard

!....................................................................

      end module module_biology
