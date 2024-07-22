










      subroutine dyn_restart(txt_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 - last update: 19-02-23
!______________________________________________________________________
      use module_principal ; use module_parallele ; use module_s
      use module_netcdfrestart
      implicit none
      character*1 txt_
      integer debug_
      integer,dimension(:,:),allocatable :: &
       datesimbis
      
!...............................................................................
! modifs: 24/12/02: derniere mise à jour
!                   + bienvenue à restart de la biologie: ichoix=11 & ichoix=22
!         25/12/02: include fichier common bio
!         25/01/03: bienvenue à  à DFVOFL_Z
!         19/03/03: suppression EVA_Z
!         21/08/03: suppression TIDESTEP1
!         11/06/04: changement de noms: chanel9_in  -> chanel9_input
!                                       chanel9_out -> chanel9_output
!                   afin d'éviter les accidents malheureux (rm *out)
!         25/04/05: attribuer un nom "daté" au fichier restart
!         26/04/05: suite du point precedent
!         02/07/06: schema tke forward permet d'eliminer le tableau TKENHZ_Z
!         04/07/06: eliminination serie tableaux ADHY
!         17/01/08: IWAVE renommé IWVE
!         01-06-09: Parallelisation
!         14-06-09: Pas d'ecriture si la sequence de l'iteration en cours
!                   n'est pas achevee
! 2010.23 12-05-11  Modification du nom d'ecriture du fichier restart
!         27-05-11  - Procedure de verification avant ecriture du fichier
!                   restart incluse dans dyn_restart.F90
!                   - la variable txt_spy_loc permet de savoir qui
!                   appelle la barriere
!         28-05-11  Fichier "check": verifier hssh-h
! S.26    01-04-13  appel a nouvelle routine elapsedtime2date
!         17-10-13  sauver alternativement dans restart_output ou restart_bis
!         17-01-15  dom_c remplace par par%rank
!         17-02-15  departuredate permet de verifier que la date de depart
!                   n'a pas ete accidentellement modifiee dans le cas d'un
!                   redemarrage d'un fichier restart
!         25-05-15  seconds_cpu_
!         13-10-15  detection des incoherences de date dans notebook_time
!         06-05-16  correction du nom du repertoire du fichier restart
!         10-04-17  ajout restartfileunits
!         08-06-17  produire un fichier graphique si NaN detectEs
!         04-06-18  ne pas ecriture de fichier restart si run_option=-1
!         14-11-18  Modif debug datesimbis
! v290    05-10-20  if(restartfileunits=='h')x1=period_*3600. !05-10-20
!...............................................................................
!    _________                    .__                  .__             ! m[0v0]m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!...............................................................................
                                                                                
!************************************************************
! Lecture d'un fichier chanel9
!************************************************************
                                                             
      debug_=0
              
      if (txt_=='r'              &
      .or.txt_=='R') then !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr>
      debug_=1
              
! Check dates consistenties using a datesimbis, a copy of datesim
      lb2=lbound(datesim) ; ub2=ubound(datesim)
      allocate(datesimbis(lb2(1):ub2(1),lb2(2):ub2(2))) ; datesimbis=datesim
                                                                            
      write(texte60,'(a,i0)')trim(restartdir_in)//'chanel9_',par%rank !06-05-16
      open(unit=9,file=texte60,access='sequential',form='unformatted')
                                                                      
      if(par%rank==0)write(6,*)'lecture chanel9 debut:'
                                                       
! insere partie_lit_1 debut: SURTOUT NE PAS EFFACER CETTE LIGNE
      read(9)ioffline_prv
      call dyn_restart_allo('a') ; call dyn_restart_bounds_r
      read(9)zone1bioflux_w
      if(allocated(zone1bioflux_glb))call netcdfrestart_read(zone1bioflux_glb,'zone1bioflux_glb',lbound(zone1bioflux_glb),ubound(zone1bioflux_glb))
      if(allocated(zone1bioflux_u))call netcdfrestart_read(zone1bioflux_u,'zone1bioflux_u',lbound(zone1bioflux_u),ubound(zone1bioflux_u))
      if(allocated(zone1bioflux_v))call netcdfrestart_read(zone1bioflux_v,'zone1bioflux_v',lbound(zone1bioflux_v),ubound(zone1bioflux_v))
      if(allocated(zone1tendancebio_glb))call netcdfrestart_read(zone1tendancebio_glb,'zone1tendancebio_glb',lbound(zone1tendancebio_glb),ubound(zone1tendancebio_glb))
      if(allocated(zone1botsurfbio_glb))call netcdfrestart_read(zone1botsurfbio_glb,'zone1botsurfbio_glb',lbound(zone1botsurfbio_glb),ubound(zone1botsurfbio_glb))
      if(allocated(zone1bioflux_glb_in))call netcdfrestart_read(zone1bioflux_glb_in,'zone1bioflux_glb_in',lbound(zone1bioflux_glb_in),ubound(zone1bioflux_glb_in))
      if(allocated(zone1bioflux_u_in))call netcdfrestart_read(zone1bioflux_u_in,'zone1bioflux_u_in',lbound(zone1bioflux_u_in),ubound(zone1bioflux_u_in))
      if(allocated(zone1bioflux_v_in))call netcdfrestart_read(zone1bioflux_v_in,'zone1bioflux_v_in',lbound(zone1bioflux_v_in),ubound(zone1bioflux_v_in))
      if(allocated(zone1bioflux_glb_out))call netcdfrestart_read(zone1bioflux_glb_out,'zone1bioflux_glb_out',lbound(zone1bioflux_glb_out),ubound(zone1bioflux_glb_out))
      if(allocated(zone1bioflux_u_out))call netcdfrestart_read(zone1bioflux_u_out,'zone1bioflux_u_out',lbound(zone1bioflux_u_out),ubound(zone1bioflux_u_out))
      if(allocated(zone1bioflux_v_out))call netcdfrestart_read(zone1bioflux_v_out,'zone1bioflux_v_out',lbound(zone1bioflux_v_out),ubound(zone1bioflux_v_out))
      if(allocated(zone1biomasst0))call netcdfrestart_read(zone1biomasst0,'zone1biomasst0',lbound(zone1biomasst0),ubound(zone1biomasst0))
      if(allocated(zone1biocumul_glb))read(9)zone1biocumul_glb
      if(allocated(zone1biocumul_loc))read(9)zone1biocumul_loc
      if(allocated(analysis3dmatrix))call netcdfrestart_read(analysis3dmatrix,'analysis3dmatrix',lbound(analysis3dmatrix),ubound(analysis3dmatrix))
      if(allocated(ema1_s_i))call netcdfrestart_read(ema1_s_i,'ema1_s_i',lbound(ema1_s_i),ubound(ema1_s_i))
      if(allocated(ema2_s_i))call netcdfrestart_read(ema2_s_i,'ema2_s_i',lbound(ema2_s_i),ubound(ema2_s_i))
      if(allocated(ema1_s_j))call netcdfrestart_read(ema1_s_j,'ema1_s_j',lbound(ema1_s_j),ubound(ema1_s_j))
      if(allocated(ema2_s_j))call netcdfrestart_read(ema2_s_j,'ema2_s_j',lbound(ema2_s_j),ubound(ema2_s_j))
      if(allocated(c_wave_mode))call netcdfrestart_read(c_wave_mode,'c_wave_mode',lbound(c_wave_mode),ubound(c_wave_mode))
      if(allocated(pcoefmode_t))call netcdfrestart_read(pcoefmode_t,'pcoefmode_t',lbound(pcoefmode_t),ubound(pcoefmode_t))
      if(allocated(ucoefmode_t))call netcdfrestart_read(ucoefmode_t,'ucoefmode_t',lbound(ucoefmode_t),ubound(ucoefmode_t))
      if(allocated(vcoefmode_t))call netcdfrestart_read(vcoefmode_t,'vcoefmode_t',lbound(vcoefmode_t),ubound(vcoefmode_t))
      if(allocated(rhrefmode_t))call netcdfrestart_read(rhrefmode_t,'rhrefmode_t',lbound(rhrefmode_t),ubound(rhrefmode_t))
      if(allocated(ema1_q_i))call netcdfrestart_read(ema1_q_i,'ema1_q_i',lbound(ema1_q_i),ubound(ema1_q_i))
      if(allocated(ema2_q_i))call netcdfrestart_read(ema2_q_i,'ema2_q_i',lbound(ema2_q_i),ubound(ema2_q_i))
      if(allocated(ema1_q_j))call netcdfrestart_read(ema1_q_j,'ema1_q_j',lbound(ema1_q_j),ubound(ema1_q_j))
      if(allocated(ema2_q_j))call netcdfrestart_read(ema2_q_j,'ema2_q_j',lbound(ema2_q_j),ubound(ema2_q_j))
      if(allocated(uv_wmode_t))call netcdfrestart_read(uv_wmode_t,'uv_wmode_t',lbound(uv_wmode_t),ubound(uv_wmode_t))
      if(allocated(analysis_p3d_t))call netcdfrestart_read(analysis_p3d_t,'analysis_p3d_t',lbound(analysis_p3d_t),ubound(analysis_p3d_t))
      if(allocated(analysis_u3d_t))call netcdfrestart_read(analysis_u3d_t,'analysis_u3d_t',lbound(analysis_u3d_t),ubound(analysis_u3d_t))
      if(allocated(analysis_v3d_t))call netcdfrestart_read(analysis_v3d_t,'analysis_v3d_t',lbound(analysis_v3d_t),ubound(analysis_v3d_t))
      if(allocated(p3dmode_cos_t))call netcdfrestart_read(p3dmode_cos_t,'p3dmode_cos_t',lbound(p3dmode_cos_t),ubound(p3dmode_cos_t))
      if(allocated(p3dmode_sin_t))call netcdfrestart_read(p3dmode_sin_t,'p3dmode_sin_t',lbound(p3dmode_sin_t),ubound(p3dmode_sin_t))
      if(allocated(u3dmode_cos_t))call netcdfrestart_read(u3dmode_cos_t,'u3dmode_cos_t',lbound(u3dmode_cos_t),ubound(u3dmode_cos_t))
      if(allocated(u3dmode_sin_t))call netcdfrestart_read(u3dmode_sin_t,'u3dmode_sin_t',lbound(u3dmode_sin_t),ubound(u3dmode_sin_t))
      if(allocated(v3dmode_cos_t))call netcdfrestart_read(v3dmode_cos_t,'v3dmode_cos_t',lbound(v3dmode_cos_t),ubound(v3dmode_cos_t))
      if(allocated(v3dmode_sin_t))call netcdfrestart_read(v3dmode_sin_t,'v3dmode_sin_t',lbound(v3dmode_sin_t),ubound(v3dmode_sin_t))
      if(allocated(analysetide3d_u))call netcdfrestart_read(analysetide3d_u,'analysetide3d_u',lbound(analysetide3d_u),ubound(analysetide3d_u))
      if(allocated(analysetide3d_v))call netcdfrestart_read(analysetide3d_v,'analysetide3d_v',lbound(analysetide3d_v),ubound(analysetide3d_v))
      if(allocated(analysetide3d_t))call netcdfrestart_read(analysetide3d_t,'analysetide3d_t',lbound(analysetide3d_t),ubound(analysetide3d_t))
      if(allocated(vel3dtidecosout_u))call netcdfrestart_read(vel3dtidecosout_u,'vel3dtidecosout_u',lbound(vel3dtidecosout_u),ubound(vel3dtidecosout_u))
      if(allocated(vel3dtidesinout_u))call netcdfrestart_read(vel3dtidesinout_u,'vel3dtidesinout_u',lbound(vel3dtidesinout_u),ubound(vel3dtidesinout_u))
      if(allocated(vel3dtidecosout_v))call netcdfrestart_read(vel3dtidecosout_v,'vel3dtidecosout_v',lbound(vel3dtidecosout_v),ubound(vel3dtidecosout_v))
      if(allocated(vel3dtidesinout_v))call netcdfrestart_read(vel3dtidesinout_v,'vel3dtidesinout_v',lbound(vel3dtidesinout_v),ubound(vel3dtidesinout_v))
      if(allocated(rhptidecosout_t))call netcdfrestart_read(rhptidecosout_t,'rhptidecosout_t',lbound(rhptidecosout_t),ubound(rhptidecosout_t))
      if(allocated(rhptidesinout_t))call netcdfrestart_read(rhptidesinout_t,'rhptidesinout_t',lbound(rhptidesinout_t),ubound(rhptidesinout_t))
      if(allocated(bio_relax_north))call netcdfrestart_read(bio_relax_north,'bio_relax_north',lbound(bio_relax_north),ubound(bio_relax_north))
      if(allocated(bio_relax_south))call netcdfrestart_read(bio_relax_south,'bio_relax_south',lbound(bio_relax_south),ubound(bio_relax_south))
      if(allocated(bio_relax_east))call netcdfrestart_read(bio_relax_east,'bio_relax_east',lbound(bio_relax_east),ubound(bio_relax_east))
      if(allocated(bio_relax_west))call netcdfrestart_read(bio_relax_west,'bio_relax_west',lbound(bio_relax_west),ubound(bio_relax_west))
      if(allocated(bio_relax_full))call netcdfrestart_read(bio_relax_full,'bio_relax_full',lbound(bio_relax_full),ubound(bio_relax_full))
      if(allocated(analysetide_u))call netcdfrestart_read(analysetide_u,'analysetide_u',lbound(analysetide_u),ubound(analysetide_u))
      if(allocated(analysetide_v))call netcdfrestart_read(analysetide_v,'analysetide_v',lbound(analysetide_v),ubound(analysetide_v))
      if(allocated(analysetide_w))call netcdfrestart_read(analysetide_w,'analysetide_w',lbound(analysetide_w),ubound(analysetide_w))
      if(allocated(analysetide2_w))call netcdfrestart_read(analysetide2_w,'analysetide2_w',lbound(analysetide2_w),ubound(analysetide2_w))
      if(allocated(cocosisi_f))call netcdfrestart_read(cocosisi_f,'cocosisi_f',lbound(cocosisi_f),ubound(cocosisi_f))
      if(allocated(gridrotcos_f))call netcdfrestart_read(gridrotcos_f,'gridrotcos_f',lbound(gridrotcos_f),ubound(gridrotcos_f))
      if(allocated(gridrotsin_f))call netcdfrestart_read(gridrotsin_f,'gridrotsin_f',lbound(gridrotsin_f),ubound(gridrotsin_f))
      if(allocated(tideanalysismatrix))call netcdfrestart_read(tideanalysismatrix,'tideanalysismatrix',lbound(tideanalysismatrix),ubound(tideanalysismatrix))
      if(allocated(iaveraged_in))call netcdfrestart_read(iaveraged_in,'iaveraged_in',lbound(iaveraged_in),ubound(iaveraged_in))
      if(allocated(iaveraged_out))call netcdfrestart_read(iaveraged_out,'iaveraged_out',lbound(iaveraged_out),ubound(iaveraged_out))
      if(allocated(iaverag1d_in))call netcdfrestart_read(iaverag1d_in,'iaverag1d_in',lbound(iaverag1d_in),ubound(iaverag1d_in))
      if(allocated(iaverag1d_out))call netcdfrestart_read(iaverag1d_out,'iaverag1d_out',lbound(iaverag1d_out),ubound(iaverag1d_out))
      if(allocated(light_kpar2_w))call netcdfrestart_read(light_kpar2_w,'light_kpar2_w',lbound(light_kpar2_w),ubound(light_kpar2_w))
      if(allocated(tidepotential_w))call netcdfrestart_read(tidepotential_w,'tidepotential_w',lbound(tidepotential_w),ubound(tidepotential_w))
      if(allocated(sshtidecos_w))call netcdfrestart_read(sshtidecos_w,'sshtidecos_w',lbound(sshtidecos_w),ubound(sshtidecos_w))
      if(allocated(sshtidesin_w))call netcdfrestart_read(sshtidesin_w,'sshtidesin_w',lbound(sshtidesin_w),ubound(sshtidesin_w))
      if(allocated(veltidecos_u))call netcdfrestart_read(veltidecos_u,'veltidecos_u',lbound(veltidecos_u),ubound(veltidecos_u))
      if(allocated(veltidesin_u))call netcdfrestart_read(veltidesin_u,'veltidesin_u',lbound(veltidesin_u),ubound(veltidesin_u))
      if(allocated(veltidecos_v))call netcdfrestart_read(veltidecos_v,'veltidecos_v',lbound(veltidecos_v),ubound(veltidecos_v))
      if(allocated(veltidesin_v))call netcdfrestart_read(veltidesin_v,'veltidesin_v',lbound(veltidesin_v),ubound(veltidesin_v))
      if(allocated(potidecos_w))call netcdfrestart_read(potidecos_w,'potidecos_w',lbound(potidecos_w),ubound(potidecos_w))
      if(allocated(potidesin_w))call netcdfrestart_read(potidesin_w,'potidesin_w',lbound(potidesin_w),ubound(potidesin_w))
      if(allocated(veltidecosout_u))call netcdfrestart_read(veltidecosout_u,'veltidecosout_u',lbound(veltidecosout_u),ubound(veltidecosout_u))
      if(allocated(veltidesinout_u))call netcdfrestart_read(veltidesinout_u,'veltidesinout_u',lbound(veltidesinout_u),ubound(veltidesinout_u))
      if(allocated(veltidecosout_v))call netcdfrestart_read(veltidecosout_v,'veltidecosout_v',lbound(veltidecosout_v),ubound(veltidecosout_v))
      if(allocated(veltidesinout_v))call netcdfrestart_read(veltidesinout_v,'veltidesinout_v',lbound(veltidesinout_v),ubound(veltidesinout_v))
      if(allocated(sshtidecosout_w))call netcdfrestart_read(sshtidecosout_w,'sshtidecosout_w',lbound(sshtidecosout_w),ubound(sshtidecosout_w))
      if(allocated(sshtidesinout_w))call netcdfrestart_read(sshtidesinout_w,'sshtidesinout_w',lbound(sshtidesinout_w),ubound(sshtidesinout_w))
      if(allocated(ssh3dtidecosout_w))call netcdfrestart_read(ssh3dtidecosout_w,'ssh3dtidecosout_w',lbound(ssh3dtidecosout_w),ubound(ssh3dtidecosout_w))
      if(allocated(ssh3dtidesinout_w))call netcdfrestart_read(ssh3dtidesinout_w,'ssh3dtidesinout_w',lbound(ssh3dtidesinout_w),ubound(ssh3dtidesinout_w))
      if(allocated(passetide))call netcdfrestart_read(passetide,'passetide',lbound(passetide),ubound(passetide))
      if(allocated(ftide))call netcdfrestart_read(ftide,'ftide',lbound(ftide),ubound(ftide))
      if(allocated(utide))call netcdfrestart_read(utide,'utide',lbound(utide),ubound(utide))
      if(allocated(northflux_sumt_v))call netcdfrestart_read(northflux_sumt_v,'northflux_sumt_v',lbound(northflux_sumt_v),ubound(northflux_sumt_v))
      if(allocated(southflux_sumt_v))call netcdfrestart_read(southflux_sumt_v,'southflux_sumt_v',lbound(southflux_sumt_v),ubound(southflux_sumt_v))
      if(allocated(drifter_l))call netcdfrestart_read(drifter_l,'drifter_l',lbound(drifter_l),ubound(drifter_l))
      if(allocated(velbot3d2d_u))call netcdfrestart_read(velbot3d2d_u,'velbot3d2d_u',lbound(velbot3d2d_u),ubound(velbot3d2d_u))
      if(allocated(velbot3d2d_v))call netcdfrestart_read(velbot3d2d_v,'velbot3d2d_v',lbound(velbot3d2d_v),ubound(velbot3d2d_v))
      if(allocated(velbot_u))call netcdfrestart_read(velbot_u,'velbot_u',lbound(velbot_u),ubound(velbot_u))
      if(allocated(velbot_v))call netcdfrestart_read(velbot_v,'velbot_v',lbound(velbot_v),ubound(velbot_v))
      if(allocated(kvectorpeak_i))call netcdfrestart_read(kvectorpeak_i,'kvectorpeak_i',lbound(kvectorpeak_i),ubound(kvectorpeak_i))
      if(allocated(kvectorpeak_j))call netcdfrestart_read(kvectorpeak_j,'kvectorpeak_j',lbound(kvectorpeak_j),ubound(kvectorpeak_j))
      if(allocated(drifter_send_nord))read(9)drifter_send_nord
      if(allocated(drifter_recv_nord))read(9)drifter_recv_nord
      if(allocated(drifter_send_sud))read(9)drifter_send_sud
      if(allocated(drifter_recv_sud))read(9)drifter_recv_sud
      if(allocated(drifter_send_est))read(9)drifter_send_est
      if(allocated(drifter_recv_est))read(9)drifter_recv_est
      if(allocated(drifter_send_ouest))read(9)drifter_send_ouest
      if(allocated(drifter_recv_ouest))read(9)drifter_recv_ouest
      if(allocated(drifter_send_nordest))read(9)drifter_send_nordest
      if(allocated(drifter_send_nordouest))read(9)drifter_send_nordouest
      if(allocated(drifter_send_sudest))read(9)drifter_send_sudest
      if(allocated(drifter_send_sudouest))read(9)drifter_send_sudouest
      if(allocated(drifter_recv_nordest))read(9)drifter_recv_nordest
      if(allocated(drifter_recv_nordouest))read(9)drifter_recv_nordouest
      if(allocated(drifter_recv_sudest))read(9)drifter_recv_sudest
      if(allocated(drifter_recv_sudouest))read(9)drifter_recv_sudouest
      if(allocated(drifter_send_canal))call netcdfrestart_read(drifter_send_canal,'drifter_send_canal',lbound(drifter_send_canal),ubound(drifter_send_canal))
      if(allocated(drifter_recv_canal))call netcdfrestart_read(drifter_recv_canal,'drifter_recv_canal',lbound(drifter_recv_canal),ubound(drifter_recv_canal))
      if(allocated(qwave_j_w))call netcdfrestart_read(qwave_j_w,'qwave_j_w',lbound(qwave_j_w),ubound(qwave_j_w))
      if(allocated(qwave_i_w))call netcdfrestart_read(qwave_i_w,'qwave_i_w',lbound(qwave_i_w),ubound(qwave_i_w))
      if(allocated(velwave_i_u))call netcdfrestart_read(velwave_i_u,'velwave_i_u',lbound(velwave_i_u),ubound(velwave_i_u))
      if(allocated(velwave_j_u))call netcdfrestart_read(velwave_j_u,'velwave_j_u',lbound(velwave_j_u),ubound(velwave_j_u))
      if(allocated(velwave_i_v))call netcdfrestart_read(velwave_i_v,'velwave_i_v',lbound(velwave_i_v),ubound(velwave_i_v))
      if(allocated(velwave_j_v))call netcdfrestart_read(velwave_j_v,'velwave_j_v',lbound(velwave_j_v),ubound(velwave_j_v))
      if(allocated(q_t))call netcdfrestart_read(q_t,'q_t',lbound(q_t),ubound(q_t))
      if(allocated(sshwave_j_w))call netcdfrestart_read(sshwave_j_w,'sshwave_j_w',lbound(sshwave_j_w),ubound(sshwave_j_w))
      if(allocated(sshwave_i_w))call netcdfrestart_read(sshwave_i_w,'sshwave_i_w',lbound(sshwave_i_w),ubound(sshwave_i_w))
      if(allocated(mc))call netcdfrestart_read(mc,'mc',lbound(mc),ubound(mc))
      if(allocated(anyv3dr4))call netcdfrestart_read(anyv3dr4,'anyv3dr4',lbound(anyv3dr4),ubound(anyv3dr4))
      if(allocated(mc_out))call netcdfrestart_read(mc_out,'mc_out',lbound(mc_out),ubound(mc_out))
      if(allocated(velbarwave_i_v))call netcdfrestart_read(velbarwave_i_v,'velbarwave_i_v',lbound(velbarwave_i_v),ubound(velbarwave_i_v))
      if(allocated(qavr_t))call netcdfrestart_read(qavr_t,'qavr_t',lbound(qavr_t),ubound(qavr_t))
      if(allocated(nhpgf_u))call netcdfrestart_read(nhpgf_u,'nhpgf_u',lbound(nhpgf_u),ubound(nhpgf_u))
      if(allocated(nhpgf_v))call netcdfrestart_read(nhpgf_v,'nhpgf_v',lbound(nhpgf_v),ubound(nhpgf_v))
      if(allocated(dsigr4_t))call netcdfrestart_read(dsigr4_t,'dsigr4_t',lbound(dsigr4_t),ubound(dsigr4_t))
      if(allocated(sig1dpom_t))call netcdfrestart_read(sig1dpom_t,'sig1dpom_t',lbound(sig1dpom_t),ubound(sig1dpom_t))
      if(allocated(obc_q_i))call netcdfrestart_read(obc_q_i,'obc_q_i',lbound(obc_q_i),ubound(obc_q_i))
      if(allocated(obc_ub_i))call netcdfrestart_read(obc_ub_i,'obc_ub_i',lbound(obc_ub_i),ubound(obc_ub_i))
      if(allocated(dsig_w))call netcdfrestart_read(dsig_w,'dsig_w',lbound(dsig_w),ubound(dsig_w))
      if(allocated(retarded_pss_w))call netcdfrestart_read(retarded_pss_w,'retarded_pss_w',lbound(retarded_pss_w),ubound(retarded_pss_w))
      if(allocated(velwf_u))call netcdfrestart_read(velwf_u,'velwf_u',lbound(velwf_u),ubound(velwf_u))
      if(allocated(velwf_v))call netcdfrestart_read(velwf_v,'velwf_v',lbound(velwf_v),ubound(velwf_v))
      if(allocated(kw_w))call netcdfrestart_read(kw_w,'kw_w',lbound(kw_w),ubound(kw_w))
      if(allocated(q2davr_w))call netcdfrestart_read(q2davr_w,'q2davr_w',lbound(q2davr_w),ubound(q2davr_w))
      if(allocated(breaker2d_t))call netcdfrestart_read(breaker2d_t,'breaker2d_t',lbound(breaker2d_t),ubound(breaker2d_t))
      if(allocated(invdx_t))call netcdfrestart_read(invdx_t,'invdx_t',lbound(invdx_t),ubound(invdx_t))
      if(allocated(invdx_f))call netcdfrestart_read(invdx_f,'invdx_f',lbound(invdx_f),ubound(invdx_f))
      if(allocated(invdx_u))call netcdfrestart_read(invdx_u,'invdx_u',lbound(invdx_u),ubound(invdx_u))
      if(allocated(invdy_u))call netcdfrestart_read(invdy_u,'invdy_u',lbound(invdy_u),ubound(invdy_u))
      if(allocated(invdy_t))call netcdfrestart_read(invdy_t,'invdy_t',lbound(invdy_t),ubound(invdy_t))
      if(allocated(invdy_f))call netcdfrestart_read(invdy_f,'invdy_f',lbound(invdy_f),ubound(invdy_f))
      if(allocated(invdy_v))call netcdfrestart_read(invdy_v,'invdy_v',lbound(invdy_v),ubound(invdy_v))
      if(allocated(invdx_v))call netcdfrestart_read(invdx_v,'invdx_v',lbound(invdx_v),ubound(invdx_v))
      if(allocated(invdxdy_t))call netcdfrestart_read(invdxdy_t,'invdxdy_t',lbound(invdxdy_t),ubound(invdxdy_t))
      if(allocated(invdxdy_u))call netcdfrestart_read(invdxdy_u,'invdxdy_u',lbound(invdxdy_u),ubound(invdxdy_u))
      if(allocated(invdxdy_v))call netcdfrestart_read(invdxdy_v,'invdxdy_v',lbound(invdxdy_v),ubound(invdxdy_v))
      if(allocated(nhpgf2d_u))call netcdfrestart_read(nhpgf2d_u,'nhpgf2d_u',lbound(nhpgf2d_u),ubound(nhpgf2d_u))
      if(allocated(nhpgf2d_v))call netcdfrestart_read(nhpgf2d_v,'nhpgf2d_v',lbound(nhpgf2d_v),ubound(nhpgf2d_v))
      if(allocated(sshmax_w))call netcdfrestart_read(sshmax_w,'sshmax_w',lbound(sshmax_w),ubound(sshmax_w))
      if(allocated(sshmin_w))call netcdfrestart_read(sshmin_w,'sshmin_w',lbound(sshmin_w),ubound(sshmin_w))
      if(allocated(sshmin_tmp_w))call netcdfrestart_read(sshmin_tmp_w,'sshmin_tmp_w',lbound(sshmin_tmp_w),ubound(sshmin_tmp_w))
      if(allocated(rwavebreak_t))call netcdfrestart_read(rwavebreak_t,'rwavebreak_t',lbound(rwavebreak_t),ubound(rwavebreak_t))
      if(allocated(breaker_t))call netcdfrestart_read(breaker_t,'breaker_t',lbound(breaker_t),ubound(breaker_t))
      if(allocated(hs_w))call netcdfrestart_read(hs_w,'hs_w',lbound(hs_w),ubound(hs_w))
      if(allocated(cx_upper_t))call netcdfrestart_read(cx_upper_t,'cx_upper_t',lbound(cx_upper_t),ubound(cx_upper_t))
      if(allocated(cy_upper_t))call netcdfrestart_read(cy_upper_t,'cy_upper_t',lbound(cy_upper_t),ubound(cy_upper_t))
      if(allocated(c_lower_t))call netcdfrestart_read(c_lower_t,'c_lower_t',lbound(c_lower_t),ubound(c_lower_t))
      if(allocated(slopemax_u))call netcdfrestart_read(slopemax_u,'slopemax_u',lbound(slopemax_u),ubound(slopemax_u))
      if(allocated(pgfwave_u))call netcdfrestart_read(pgfwave_u,'pgfwave_u',lbound(pgfwave_u),ubound(pgfwave_u))
      if(allocated(dt_over_retartedtime))read(9)dt_over_retartedtime
      if(allocated(sshrelax_j_w))call netcdfrestart_read(sshrelax_j_w,'sshrelax_j_w',lbound(sshrelax_j_w),ubound(sshrelax_j_w))
      if(allocated(sshrelax_i_w))call netcdfrestart_read(sshrelax_i_w,'sshrelax_i_w',lbound(sshrelax_i_w),ubound(sshrelax_i_w))
      read(9)asselin_nh
      read(9)wetdry_cstnh
      read(9)cfl_nh
      read(9)cwavepeak
      read(9)periodpeak
      read(9)freq2pipeak
      read(9)kvectorpeak
      read(9)kvector
      read(9)mukvector
      read(9)nhpgf_reduce
      read(9)acm_speed
      read(9)inv_brkh
      read(9)inv_brkslope2
      read(9)inv_brkvar
      read(9)brk_crit_slope
      read(9)brk_crit_var
      read(9)brk_crit_h
      read(9)brk_crit_r
      read(9)nh2d_graph_period
      read(9)nh2d_graph_spinup
      read(9)sum_avr_hs
      read(9)nh_frozensigma
      read(9)nh_wavebreakscheme
      read(9)flag_adve2d
      read(9)flag_adve3d
      read(9)flag_nh3d
      read(9)flag_nh3d_none
      read(9)flag_nh3d_nosplit_uv
      read(9)flag_nh3d_nosplit_tsuv
      read(9)flag_nh3d_timesplit_tsuv
      read(9)flag_timesplitting
      read(9)flag_timesplitting_adve_uv
      read(9)flagspo_i1
      read(9)flagspo_i2
      read(9)flagspo_j1
      read(9)flagspo_j2
      read(9)flag_ksloffline
      read(9)flag_groundwater
      read(9)flag_surfriver
      read(9)ofl_reversedtime
      read(9)obcnh_scheme
      read(9)obc_scheme
      read(9)drifter_random_w
      read(9)flag_buo_w
      read(9)flag_ogcmtidemixing
      read(9)flag_ogcm_instab
      read(9)flag_ofactor
      read(9)flag_omega_cumul
      read(9)flag_negdif_ver
      read(9)flag_negdif_hor
      read(9)flag_z0_macro
      read(9)oasis_symsym_onoff
      read(9)oasis_symsym_retrots
      if(allocated(variance_timeaveraged))call netcdfrestart_read(variance_timeaveraged,'variance_timeaveraged',lbound(variance_timeaveraged),ubound(variance_timeaveraged))
      if(allocated(ssh_timeaveraged))call netcdfrestart_read(ssh_timeaveraged,'ssh_timeaveraged',lbound(ssh_timeaveraged),ubound(ssh_timeaveraged))
      if(allocated(flx2d_timeaveraged))call netcdfrestart_read(flx2d_timeaveraged,'flx2d_timeaveraged',lbound(flx2d_timeaveraged),ubound(flx2d_timeaveraged))
      if(allocated(fly2d_timeaveraged))call netcdfrestart_read(fly2d_timeaveraged,'fly2d_timeaveraged',lbound(fly2d_timeaveraged),ubound(fly2d_timeaveraged))
      if(allocated(flx3d_timeaveraged))call netcdfrestart_read(flx3d_timeaveraged,'flx3d_timeaveraged',lbound(flx3d_timeaveraged),ubound(flx3d_timeaveraged))
      if(allocated(fly3d_timeaveraged))call netcdfrestart_read(fly3d_timeaveraged,'fly3d_timeaveraged',lbound(fly3d_timeaveraged),ubound(fly3d_timeaveraged))
      if(allocated(u_euler_timeaveraged))call netcdfrestart_read(u_euler_timeaveraged,'u_euler_timeaveraged',lbound(u_euler_timeaveraged),ubound(u_euler_timeaveraged))
      if(allocated(v_euler_timeaveraged))call netcdfrestart_read(v_euler_timeaveraged,'v_euler_timeaveraged',lbound(v_euler_timeaveraged),ubound(v_euler_timeaveraged))
      if(allocated(tkeb_w))call netcdfrestart_read(tkeb_w,'tkeb_w',lbound(tkeb_w),ubound(tkeb_w))
      if(allocated(tken_w))call netcdfrestart_read(tken_w,'tken_w',lbound(tken_w),ubound(tken_w))
      if(allocated(tkea_w))call netcdfrestart_read(tkea_w,'tkea_w',lbound(tkea_w),ubound(tkea_w))
      if(allocated(tkle_w))call netcdfrestart_read(tkle_w,'tkle_w',lbound(tkle_w),ubound(tkle_w))
      if(allocated(tkll_w))call netcdfrestart_read(tkll_w,'tkll_w',lbound(tkll_w),ubound(tkll_w))
      if(allocated(epsb_w))call netcdfrestart_read(epsb_w,'epsb_w',lbound(epsb_w),ubound(epsb_w))
      if(allocated(epsn_w))call netcdfrestart_read(epsn_w,'epsn_w',lbound(epsn_w),ubound(epsn_w))
      if(allocated(epsa_w))call netcdfrestart_read(epsa_w,'epsa_w',lbound(epsa_w),ubound(epsa_w))
      if(allocated(gradssh_u))call netcdfrestart_read(gradssh_u,'gradssh_u',lbound(gradssh_u),ubound(gradssh_u))
      if(allocated(gradssh_v))call netcdfrestart_read(gradssh_v,'gradssh_v',lbound(gradssh_v),ubound(gradssh_v))
      if(allocated(pgf_u))call netcdfrestart_read(pgf_u,'pgf_u',lbound(pgf_u),ubound(pgf_u))
      if(allocated(pgf_v))call netcdfrestart_read(pgf_v,'pgf_v',lbound(pgf_v),ubound(pgf_v))
      if(allocated(dz_u))call netcdfrestart_read(dz_u,'dz_u',lbound(dz_u),ubound(dz_u))
      if(allocated(dz_v))call netcdfrestart_read(dz_v,'dz_v',lbound(dz_v),ubound(dz_v))
      if(allocated(dz_t))call netcdfrestart_read(dz_t,'dz_t',lbound(dz_t),ubound(dz_t))
      if(allocated(depth_t))call netcdfrestart_read(depth_t,'depth_t',lbound(depth_t),ubound(depth_t))
      if(allocated(depth_u))call netcdfrestart_read(depth_u,'depth_u',lbound(depth_u),ubound(depth_u))
      if(allocated(depth_v))call netcdfrestart_read(depth_v,'depth_v',lbound(depth_v),ubound(depth_v))
      if(allocated(depth_w))call netcdfrestart_read(depth_w,'depth_w',lbound(depth_w),ubound(depth_w))
      if(allocated(depth_f))call netcdfrestart_read(depth_f,'depth_f',lbound(depth_f),ubound(depth_f))
      if(allocated(stokesforces_u))call netcdfrestart_read(stokesforces_u,'stokesforces_u',lbound(stokesforces_u),ubound(stokesforces_u))
      if(allocated(stokesforces_v))call netcdfrestart_read(stokesforces_v,'stokesforces_v',lbound(stokesforces_v),ubound(stokesforces_v))
      if(allocated(r_mangrove))call netcdfrestart_read(r_mangrove,'r_mangrove',lbound(r_mangrove),ubound(r_mangrove))
      if(allocated(sigma_w))call netcdfrestart_read(sigma_w,'sigma_w',lbound(sigma_w),ubound(sigma_w))
      if(allocated(dsig_f))call netcdfrestart_read(dsig_f,'dsig_f',lbound(dsig_f),ubound(dsig_f))
      if(allocated(rhp_t))call netcdfrestart_read(rhp_t,'rhp_t',lbound(rhp_t),ubound(rhp_t))
      if(allocated(rhcref_t))call netcdfrestart_read(rhcref_t,'rhcref_t',lbound(rhcref_t),ubound(rhcref_t))
      if(allocated(vel_u))call netcdfrestart_read(vel_u,'vel_u',lbound(vel_u),ubound(vel_u))
      if(allocated(vel_v))call netcdfrestart_read(vel_v,'vel_v',lbound(vel_v),ubound(vel_v))
      if(allocated(presgrad_u))call netcdfrestart_read(presgrad_u,'presgrad_u',lbound(presgrad_u),ubound(presgrad_u))
      if(allocated(presgrad_v))call netcdfrestart_read(presgrad_v,'presgrad_v',lbound(presgrad_v),ubound(presgrad_v))
      if(allocated(omega_w))call netcdfrestart_read(omega_w,'omega_w',lbound(omega_w),ubound(omega_w))
      if(allocated(veldydz_u))call netcdfrestart_read(veldydz_u,'veldydz_u',lbound(veldydz_u),ubound(veldydz_u))
      if(allocated(veldxdz_v))call netcdfrestart_read(veldxdz_v,'veldxdz_v',lbound(veldxdz_v),ubound(veldxdz_v))
      if(allocated(tem_t))call netcdfrestart_read(tem_t,'tem_t',lbound(tem_t),ubound(tem_t))
      if(allocated(sal_t))call netcdfrestart_read(sal_t,'sal_t',lbound(sal_t),ubound(sal_t))
      if(allocated(tridia_in))call netcdfrestart_read(tridia_in,'tridia_in',lbound(tridia_in),ubound(tridia_in))
      if(allocated(hr_z2lr_w))call netcdfrestart_read(hr_z2lr_w,'hr_z2lr_w',lbound(hr_z2lr_w),ubound(hr_z2lr_w))
      if(allocated(anyv3d))call netcdfrestart_read(anyv3d,'anyv3d',lbound(anyv3d),ubound(anyv3d))
      if(allocated(sponge_t))call netcdfrestart_read(sponge_t,'sponge_t',lbound(sponge_t),ubound(sponge_t))
      if(allocated(sponge_u))call netcdfrestart_read(sponge_u,'sponge_u',lbound(sponge_u),ubound(sponge_u))
      if(allocated(sponge_v))call netcdfrestart_read(sponge_v,'sponge_v',lbound(sponge_v),ubound(sponge_v))
      if(allocated(uwind_t))call netcdfrestart_read(uwind_t,'uwind_t',lbound(uwind_t),ubound(uwind_t))
      if(allocated(vwind_t))call netcdfrestart_read(vwind_t,'vwind_t',lbound(vwind_t),ubound(vwind_t))
      if(allocated(uwind100_t))call netcdfrestart_read(uwind100_t,'uwind100_t',lbound(uwind100_t),ubound(uwind100_t))
      if(allocated(vwind100_t))call netcdfrestart_read(vwind100_t,'vwind100_t',lbound(vwind100_t),ubound(vwind100_t))
      if(allocated(sshf_w))call netcdfrestart_read(sshf_w,'sshf_w',lbound(sshf_w),ubound(sshf_w))
      if(allocated(slhf_w))call netcdfrestart_read(slhf_w,'slhf_w',lbound(slhf_w),ubound(slhf_w))
      if(allocated(ssr_w))call netcdfrestart_read(ssr_w,'ssr_w',lbound(ssr_w),ubound(ssr_w))
      if(allocated(ssr24prv_w))call netcdfrestart_read(ssr24prv_w,'ssr24prv_w',lbound(ssr24prv_w),ubound(ssr24prv_w))
      if(allocated(snsf_w))call netcdfrestart_read(snsf_w,'snsf_w',lbound(snsf_w),ubound(snsf_w))
      if(allocated(precipi_w))call netcdfrestart_read(precipi_w,'precipi_w',lbound(precipi_w),ubound(precipi_w))
      if(allocated(taux_w))call netcdfrestart_read(taux_w,'taux_w',lbound(taux_w),ubound(taux_w))
      if(allocated(tauy_w))call netcdfrestart_read(tauy_w,'tauy_w',lbound(tauy_w),ubound(tauy_w))
      if(allocated(sigma_f))call netcdfrestart_read(sigma_f,'sigma_f',lbound(sigma_f),ubound(sigma_f))
      if(allocated(km_w))call netcdfrestart_read(km_w,'km_w',lbound(km_w),ubound(km_w))
      if(allocated(kh_w))call netcdfrestart_read(kh_w,'kh_w',lbound(kh_w),ubound(kh_w))
      if(allocated(rho_t))call netcdfrestart_read(rho_t,'rho_t',lbound(rho_t),ubound(rho_t))
      if(allocated(omega_evaprec_w))call netcdfrestart_read(omega_evaprec_w,'omega_evaprec_w',lbound(omega_evaprec_w),ubound(omega_evaprec_w))
      if(allocated(tridia_out))call netcdfrestart_read(tridia_out,'tridia_out',lbound(tridia_out),ubound(tridia_out))
      if(allocated(fluxbar_sumt_u))call netcdfrestart_read(fluxbar_sumt_u,'fluxbar_sumt_u',lbound(fluxbar_sumt_u),ubound(fluxbar_sumt_u))
      if(allocated(fluxbar_sumt_v))call netcdfrestart_read(fluxbar_sumt_v,'fluxbar_sumt_v',lbound(fluxbar_sumt_v),ubound(fluxbar_sumt_v))
      if(allocated(velbar_u))call netcdfrestart_read(velbar_u,'velbar_u',lbound(velbar_u),ubound(velbar_u))
      if(allocated(uflux2d_f))call netcdfrestart_read(uflux2d_f,'uflux2d_f',lbound(uflux2d_f),ubound(uflux2d_f))
      if(allocated(vflux2d_f))call netcdfrestart_read(vflux2d_f,'vflux2d_f',lbound(vflux2d_f),ubound(vflux2d_f))
      if(allocated(vlxbar_f))call netcdfrestart_read(vlxbar_f,'vlxbar_f',lbound(vlxbar_f),ubound(vlxbar_f))
      if(allocated(velbar_v))call netcdfrestart_read(velbar_v,'velbar_v',lbound(velbar_v),ubound(velbar_v))
      if(allocated(vlybar_f))call netcdfrestart_read(vlybar_f,'vlybar_f',lbound(vlybar_f),ubound(vlybar_f))
      if(allocated(flux2d_u))call netcdfrestart_read(flux2d_u,'flux2d_u',lbound(flux2d_u),ubound(flux2d_u))
      if(allocated(flux2d_v))call netcdfrestart_read(flux2d_v,'flux2d_v',lbound(flux2d_v),ubound(flux2d_v))
      if(allocated(ssh_w))call netcdfrestart_read(ssh_w,'ssh_w',lbound(ssh_w),ubound(ssh_w))
      if(allocated(hssh_f))call netcdfrestart_read(hssh_f,'hssh_f',lbound(hssh_f),ubound(hssh_f))
      if(allocated(hz_u))call netcdfrestart_read(hz_u,'hz_u',lbound(hz_u),ubound(hz_u))
      if(allocated(hz_v))call netcdfrestart_read(hz_v,'hz_v',lbound(hz_v),ubound(hz_v))
      if(allocated(hz_w))call netcdfrestart_read(hz_w,'hz_w',lbound(hz_w),ubound(hz_w))
      if(allocated(heatrelax_w))call netcdfrestart_read(heatrelax_w,'heatrelax_w',lbound(heatrelax_w),ubound(heatrelax_w))
      if(allocated(xy_u))call netcdfrestart_read(xy_u,'xy_u',lbound(xy_u),ubound(xy_u))
      if(allocated(xy_v))call netcdfrestart_read(xy_v,'xy_v',lbound(xy_v),ubound(xy_v))
      if(allocated(xy_t))call netcdfrestart_read(xy_t,'xy_t',lbound(xy_t),ubound(xy_t))
      if(allocated(xy_f))call netcdfrestart_read(xy_f,'xy_f',lbound(xy_f),ubound(xy_f))
      if(allocated(pss_w))call netcdfrestart_read(pss_w,'pss_w',lbound(pss_w),ubound(pss_w))
      if(allocated(wstress_u))call netcdfrestart_read(wstress_u,'wstress_u',lbound(wstress_u),ubound(wstress_u))
      if(allocated(wstress_v))call netcdfrestart_read(wstress_v,'wstress_v',lbound(wstress_v),ubound(wstress_v))
      if(allocated(ssh_int_w))call netcdfrestart_read(ssh_int_w,'ssh_int_w',lbound(ssh_int_w),ubound(ssh_int_w))
      if(allocated(ustokesvortex_t))call netcdfrestart_read(ustokesvortex_t,'ustokesvortex_t',lbound(ustokesvortex_t),ubound(ustokesvortex_t))
      if(allocated(ustokesvortex_f))call netcdfrestart_read(ustokesvortex_f,'ustokesvortex_f',lbound(ustokesvortex_f),ubound(ustokesvortex_f))
      if(allocated(vstokesvortex_t))call netcdfrestart_read(vstokesvortex_t,'vstokesvortex_t',lbound(vstokesvortex_t),ubound(vstokesvortex_t))
      if(allocated(vstokesvortex_f))call netcdfrestart_read(vstokesvortex_f,'vstokesvortex_f',lbound(vstokesvortex_f),ubound(vstokesvortex_f))
      if(allocated(velavr_u))call netcdfrestart_read(velavr_u,'velavr_u',lbound(velavr_u),ubound(velavr_u))
      if(allocated(velavr_v))call netcdfrestart_read(velavr_v,'velavr_v',lbound(velavr_v),ubound(velavr_v))
      if(allocated(timefilter_u))call netcdfrestart_read(timefilter_u,'timefilter_u',lbound(timefilter_u),ubound(timefilter_u))
      if(allocated(timefilter_v))call netcdfrestart_read(timefilter_v,'timefilter_v',lbound(timefilter_v),ubound(timefilter_v))
      if(allocated(teta0_t))call netcdfrestart_read(teta0_t,'teta0_t',lbound(teta0_t),ubound(teta0_t))
      if(allocated(teta2_t))call netcdfrestart_read(teta2_t,'teta2_t',lbound(teta2_t),ubound(teta2_t))
      if(allocated(q2_t))call netcdfrestart_read(q2_t,'q2_t',lbound(q2_t),ubound(q2_t))
      if(allocated(teta2delta_t))call netcdfrestart_read(teta2delta_t,'teta2delta_t',lbound(teta2delta_t),ubound(teta2delta_t))
      if(allocated(q2delta_t))call netcdfrestart_read(q2delta_t,'q2delta_t',lbound(q2delta_t),ubound(q2delta_t))
      if(allocated(uwinddelta_t))call netcdfrestart_read(uwinddelta_t,'uwinddelta_t',lbound(uwinddelta_t),ubound(uwinddelta_t))
      if(allocated(vwinddelta_t))call netcdfrestart_read(vwinddelta_t,'vwinddelta_t',lbound(vwinddelta_t),ubound(vwinddelta_t))
      if(allocated(tetastar_t))call netcdfrestart_read(tetastar_t,'tetastar_t',lbound(tetastar_t),ubound(tetastar_t))
      if(allocated(ustar_t))call netcdfrestart_read(ustar_t,'ustar_t',lbound(ustar_t),ubound(ustar_t))
      if(allocated(qstar_t))call netcdfrestart_read(qstar_t,'qstar_t',lbound(qstar_t),ubound(qstar_t))
      if(allocated(frozenterm2d_u))call netcdfrestart_read(frozenterm2d_u,'frozenterm2d_u',lbound(frozenterm2d_u),ubound(frozenterm2d_u))
      if(allocated(frozenterm2d_v))call netcdfrestart_read(frozenterm2d_v,'frozenterm2d_v',lbound(frozenterm2d_v),ubound(frozenterm2d_v))
      if(allocated(frozenterm3d_u))call netcdfrestart_read(frozenterm3d_u,'frozenterm3d_u',lbound(frozenterm3d_u),ubound(frozenterm3d_u))
      if(allocated(frozenterm3d_v))call netcdfrestart_read(frozenterm3d_v,'frozenterm3d_v',lbound(frozenterm3d_v),ubound(frozenterm3d_v))
      if(allocated(fluxbar_u))call netcdfrestart_read(fluxbar_u,'fluxbar_u',lbound(fluxbar_u),ubound(fluxbar_u))
      if(allocated(fluxbar_v))call netcdfrestart_read(fluxbar_v,'fluxbar_v',lbound(fluxbar_v),ubound(fluxbar_v))
      if(allocated(lon_t))call netcdfrestart_read(lon_t,'lon_t',lbound(lon_t),ubound(lon_t))
      if(allocated(lat_t))call netcdfrestart_read(lat_t,'lat_t',lbound(lat_t),ubound(lat_t))
      if(allocated(lon_u))call netcdfrestart_read(lon_u,'lon_u',lbound(lon_u),ubound(lon_u))
      if(allocated(lat_u))call netcdfrestart_read(lat_u,'lat_u',lbound(lat_u),ubound(lat_u))
      if(allocated(lon_v))call netcdfrestart_read(lon_v,'lon_v',lbound(lon_v),ubound(lon_v))
      if(allocated(lat_v))call netcdfrestart_read(lat_v,'lat_v',lbound(lat_v),ubound(lat_v))
      if(allocated(lon_f))call netcdfrestart_read(lon_f,'lon_f',lbound(lon_f),ubound(lon_f))
      if(allocated(lat_f))call netcdfrestart_read(lat_f,'lat_f',lbound(lat_f),ubound(lat_f))
      if(allocated(globlon_t))call netcdfrestart_read(globlon_t,'globlon_t',lbound(globlon_t),ubound(globlon_t))
      if(allocated(globlat_t))call netcdfrestart_read(globlat_t,'globlat_t',lbound(globlat_t),ubound(globlat_t))
      if(allocated(globlon_u))call netcdfrestart_read(globlon_u,'globlon_u',lbound(globlon_u),ubound(globlon_u))
      if(allocated(globlat_u))call netcdfrestart_read(globlat_u,'globlat_u',lbound(globlat_u),ubound(globlat_u))
      if(allocated(globlon_v))call netcdfrestart_read(globlon_v,'globlon_v',lbound(globlon_v),ubound(globlon_v))
      if(allocated(globlat_v))call netcdfrestart_read(globlat_v,'globlat_v',lbound(globlat_v),ubound(globlat_v))
      if(ioffline_prv>=1.and.allocated(fluxbarsum_ofl_u))call netcdfrestart_read(fluxbarsum_ofl_u,'fluxbarsum_ofl_u',lbound(fluxbarsum_ofl_u),ubound(fluxbarsum_ofl_u))
      if(ioffline_prv>=1.and.allocated(fluxbarsum_ofl_v))call netcdfrestart_read(fluxbarsum_ofl_v,'fluxbarsum_ofl_v',lbound(fluxbarsum_ofl_v),ubound(fluxbarsum_ofl_v))
      if(allocated(dxdy_u))call netcdfrestart_read(dxdy_u,'dxdy_u',lbound(dxdy_u),ubound(dxdy_u))
      if(allocated(dxdy_v))call netcdfrestart_read(dxdy_v,'dxdy_v',lbound(dxdy_v),ubound(dxdy_v))
      if(allocated(dxdy_t))call netcdfrestart_read(dxdy_t,'dxdy_t',lbound(dxdy_t),ubound(dxdy_t))
      if(allocated(dxdy_e_t))call netcdfrestart_read(dxdy_e_t,'dxdy_e_t',lbound(dxdy_e_t),ubound(dxdy_e_t))
      if(allocated(dx_e_f))call netcdfrestart_read(dx_e_f,'dx_e_f',lbound(dx_e_f),ubound(dx_e_f))
      if(allocated(dy_e_f))call netcdfrestart_read(dy_e_f,'dy_e_f',lbound(dy_e_f),ubound(dy_e_f))
      if(allocated(dx_u))call netcdfrestart_read(dx_u,'dx_u',lbound(dx_u),ubound(dx_u))
      if(allocated(dx_v))call netcdfrestart_read(dx_v,'dx_v',lbound(dx_v),ubound(dx_v))
      if(allocated(dx_t))call netcdfrestart_read(dx_t,'dx_t',lbound(dx_t),ubound(dx_t))
      if(allocated(dx_f))call netcdfrestart_read(dx_f,'dx_f',lbound(dx_f),ubound(dx_f))
      if(allocated(dy_u))call netcdfrestart_read(dy_u,'dy_u',lbound(dy_u),ubound(dy_u))
      if(allocated(dy_v))call netcdfrestart_read(dy_v,'dy_v',lbound(dy_v),ubound(dy_v))
      if(allocated(dy_t))call netcdfrestart_read(dy_t,'dy_t',lbound(dy_t),ubound(dy_t))
      if(allocated(dy_f))call netcdfrestart_read(dy_f,'dy_f',lbound(dy_f),ubound(dy_f))
      if(allocated(h_w))call netcdfrestart_read(h_w,'h_w',lbound(h_w),ubound(h_w))
      if(allocated(hcopy_w))call netcdfrestart_read(hcopy_w,'hcopy_w',lbound(hcopy_w),ubound(hcopy_w))
      if(allocated(h0_w))call netcdfrestart_read(h0_w,'h0_w',lbound(h0_w),ubound(h0_w))
      if(allocated(h_u))call netcdfrestart_read(h_u,'h_u',lbound(h_u),ubound(h_u))
      if(allocated(h_v))call netcdfrestart_read(h_v,'h_v',lbound(h_v),ubound(h_v))
      if(allocated(h_f))call netcdfrestart_read(h_f,'h_f',lbound(h_f),ubound(h_f))
      if(allocated(coriolis_t))call netcdfrestart_read(coriolis_t,'coriolis_t',lbound(coriolis_t),ubound(coriolis_t))
      if(allocated(coriolis_f))call netcdfrestart_read(coriolis_f,'coriolis_f',lbound(coriolis_f),ubound(coriolis_f))
      if(allocated(rhpzavr_w))call netcdfrestart_read(rhpzavr_w,'rhpzavr_w',lbound(rhpzavr_w),ubound(rhpzavr_w))
      if(allocated(q10_t))call netcdfrestart_read(q10_t,'q10_t',lbound(q10_t),ubound(q10_t))
      if(allocated(teta10_t))call netcdfrestart_read(teta10_t,'teta10_t',lbound(teta10_t),ubound(teta10_t))
      if(allocated(fric_u))call netcdfrestart_read(fric_u,'fric_u',lbound(fric_u),ubound(fric_u))
      if(allocated(fric_v))call netcdfrestart_read(fric_v,'fric_v',lbound(fric_v),ubound(fric_v))
      if(allocated(fric_t))call netcdfrestart_read(fric_t,'fric_t',lbound(fric_t),ubound(fric_t))
      if(allocated(cdb_t))call netcdfrestart_read(cdb_t,'cdb_t',lbound(cdb_t),ubound(cdb_t))
      if(allocated(cdb_f))call netcdfrestart_read(cdb_f,'cdb_f',lbound(cdb_f),ubound(cdb_f))
      if(allocated(xflux_t))call netcdfrestart_read(xflux_t,'xflux_t',lbound(xflux_t),ubound(xflux_t))
      if(allocated(xflux_f))call netcdfrestart_read(xflux_f,'xflux_f',lbound(xflux_f),ubound(xflux_f))
      if(allocated(yflux_t))call netcdfrestart_read(yflux_t,'yflux_t',lbound(yflux_t),ubound(yflux_t))
      if(allocated(yflux_f))call netcdfrestart_read(yflux_f,'yflux_f',lbound(yflux_f),ubound(yflux_f))
      if(allocated(pres3d2d_u))call netcdfrestart_read(pres3d2d_u,'pres3d2d_u',lbound(pres3d2d_u),ubound(pres3d2d_u))
      if(allocated(pres3d2d_v))call netcdfrestart_read(pres3d2d_v,'pres3d2d_v',lbound(pres3d2d_v),ubound(pres3d2d_v))
      if(allocated(adve3d2d_u))call netcdfrestart_read(adve3d2d_u,'adve3d2d_u',lbound(adve3d2d_u),ubound(adve3d2d_u))
      if(allocated(adve3d2d_v))call netcdfrestart_read(adve3d2d_v,'adve3d2d_v',lbound(adve3d2d_v),ubound(adve3d2d_v))
      if(allocated(rmangrovebar))call netcdfrestart_read(rmangrovebar,'rmangrovebar',lbound(rmangrovebar),ubound(rmangrovebar))
      if(allocated(mang3dto2d_u))call netcdfrestart_read(mang3dto2d_u,'mang3dto2d_u',lbound(mang3dto2d_u),ubound(mang3dto2d_u))
      if(allocated(mang3dto2d_v))call netcdfrestart_read(mang3dto2d_v,'mang3dto2d_v',lbound(mang3dto2d_v),ubound(mang3dto2d_v))
      if(allocated(restoring3d2d_u))call netcdfrestart_read(restoring3d2d_u,'restoring3d2d_u',lbound(restoring3d2d_u),ubound(restoring3d2d_u))
      if(allocated(restoring3d2d_v))call netcdfrestart_read(restoring3d2d_v,'restoring3d2d_v',lbound(restoring3d2d_v),ubound(restoring3d2d_v))
      if(allocated(stokesforces3d2d_u))call netcdfrestart_read(stokesforces3d2d_u,'stokesforces3d2d_u',lbound(stokesforces3d2d_u),ubound(stokesforces3d2d_u))
      if(allocated(stokesforces3d2d_v))call netcdfrestart_read(stokesforces3d2d_v,'stokesforces3d2d_v',lbound(stokesforces3d2d_v),ubound(stokesforces3d2d_v))
      if(allocated(wstress_w))call netcdfrestart_read(wstress_w,'wstress_w',lbound(wstress_w),ubound(wstress_w))
      if(allocated(z0_w))call netcdfrestart_read(z0_w,'z0_w',lbound(z0_w),ubound(z0_w))
      if(allocated(albedo_w))call netcdfrestart_read(albedo_w,'albedo_w',lbound(albedo_w),ubound(albedo_w))
      if(allocated(gridrotcos_t))call netcdfrestart_read(gridrotcos_t,'gridrotcos_t',lbound(gridrotcos_t),ubound(gridrotcos_t))
      if(allocated(gridrotsin_t))call netcdfrestart_read(gridrotsin_t,'gridrotsin_t',lbound(gridrotsin_t),ubound(gridrotsin_t))
      if(allocated(grid_angle_t))call netcdfrestart_read(grid_angle_t,'grid_angle_t',lbound(grid_angle_t),ubound(grid_angle_t))
      if(allocated(sshstokes_w))call netcdfrestart_read(sshstokes_w,'sshstokes_w',lbound(sshstokes_w),ubound(sshstokes_w))
      if(allocated(cwi_int_u))call netcdfrestart_read(cwi_int_u,'cwi_int_u',lbound(cwi_int_u),ubound(cwi_int_u))
      if(allocated(cwi_int_v))call netcdfrestart_read(cwi_int_v,'cwi_int_v',lbound(cwi_int_v),ubound(cwi_int_v))
      if(allocated(cwj_int_u))call netcdfrestart_read(cwj_int_u,'cwj_int_u',lbound(cwj_int_u),ubound(cwj_int_u))
      if(allocated(cwj_int_v))call netcdfrestart_read(cwj_int_v,'cwj_int_v',lbound(cwj_int_v),ubound(cwj_int_v))
      if(allocated(sshrefobc_i))call netcdfrestart_read(sshrefobc_i,'sshrefobc_i',lbound(sshrefobc_i),ubound(sshrefobc_i))
      if(allocated(vbrrefobc_i))call netcdfrestart_read(vbrrefobc_i,'vbrrefobc_i',lbound(vbrrefobc_i),ubound(vbrrefobc_i))
      if(allocated(sshrefobc_j))call netcdfrestart_read(sshrefobc_j,'sshrefobc_j',lbound(sshrefobc_j),ubound(sshrefobc_j))
      if(allocated(vbrrefobc_j))call netcdfrestart_read(vbrrefobc_j,'vbrrefobc_j',lbound(vbrrefobc_j),ubound(vbrrefobc_j))
      if(allocated(anyv1d))call netcdfrestart_read(anyv1d,'anyv1d',lbound(anyv1d),ubound(anyv1d))
      if(allocated(gridcard))call netcdfrestart_read(gridcard,'gridcard',lbound(gridcard),ubound(gridcard))
      if(allocated(novector))call netcdfrestart_read(novector,'novector',lbound(novector),ubound(novector))
      if(allocated(proficard))call netcdfrestart_read(proficard,'proficard',lbound(proficard),ubound(proficard))
      if(allocated(airseainfo))call netcdfrestart_read(airseainfo,'airseainfo',lbound(airseainfo),ubound(airseainfo))
      if(allocated(airseadt))call netcdfrestart_read(airseadt,'airseadt',lbound(airseadt),ubound(airseadt))
      if(allocated(riverdt))call netcdfrestart_read(riverdt,'riverdt',lbound(riverdt),ubound(riverdt))
      if(allocated(river_t))call netcdfrestart_read(river_t,'river_t',lbound(river_t),ubound(river_t))
      if(allocated(riverflux))call netcdfrestart_read(riverflux,'riverflux',lbound(riverflux),ubound(riverflux))
      if(allocated(mask_t))call netcdfrestart_read(mask_t,'mask_t',lbound(mask_t),ubound(mask_t))
      if(allocated(mask_f))call netcdfrestart_read(mask_f,'mask_f',lbound(mask_f),ubound(mask_f))
      if(allocated(mask_u))call netcdfrestart_read(mask_u,'mask_u',lbound(mask_u),ubound(mask_u))
      if(allocated(mask_v))call netcdfrestart_read(mask_v,'mask_v',lbound(mask_v),ubound(mask_v))
      if(allocated(mask_vqs_tke_w))call netcdfrestart_read(mask_vqs_tke_w,'mask_vqs_tke_w',lbound(mask_vqs_tke_w),ubound(mask_vqs_tke_w))
      if(allocated(canaldir))call netcdfrestart_read(canaldir,'canaldir',lbound(canaldir),ubound(canaldir))
      if(allocated(wetmask_wi_t))call netcdfrestart_read(wetmask_wi_t,'wetmask_wi_t',lbound(wetmask_wi_t),ubound(wetmask_wi_t))
      if(allocated(lonlat2ij_t))call netcdfrestart_read(lonlat2ij_t,'lonlat2ij_t',lbound(lonlat2ij_t),ubound(lonlat2ij_t))
      if(allocated(canalmpioverlap))call netcdfrestart_read(canalmpioverlap,'canalmpioverlap',lbound(canalmpioverlap),ubound(canalmpioverlap))
      if(allocated(sodate))call netcdfrestart_read(sodate,'sodate',lbound(sodate),ubound(sodate))
      if(allocated(canalrank))call netcdfrestart_read(canalrank,'canalrank',lbound(canalrank),ubound(canalrank))
      if(allocated(canalrankbis))call netcdfrestart_read(canalrankbis,'canalrankbis',lbound(canalrankbis),ubound(canalrankbis))
      if(allocated(i_canalcoord))call netcdfrestart_read(i_canalcoord,'i_canalcoord',lbound(i_canalcoord),ubound(i_canalcoord))
      if(allocated(j_canalcoord))call netcdfrestart_read(j_canalcoord,'j_canalcoord',lbound(j_canalcoord),ubound(j_canalcoord))
      if(allocated(kmin_u))call netcdfrestart_read(kmin_u,'kmin_u',lbound(kmin_u),ubound(kmin_u))
      if(allocated(kmin_v))call netcdfrestart_read(kmin_v,'kmin_v',lbound(kmin_v),ubound(kmin_v))
      if(allocated(kmin_w))call netcdfrestart_read(kmin_w,'kmin_w',lbound(kmin_w),ubound(kmin_w))
      if(allocated(kundermin_t))call netcdfrestart_read(kundermin_t,'kundermin_t',lbound(kundermin_t),ubound(kundermin_t))
      if(allocated(kundermin_u))call netcdfrestart_read(kundermin_u,'kundermin_u',lbound(kundermin_u),ubound(kundermin_u))
      if(allocated(kundermin_v))call netcdfrestart_read(kundermin_v,'kundermin_v',lbound(kundermin_v),ubound(kundermin_v))
      if(allocated(kmerged_t))call netcdfrestart_read(kmerged_t,'kmerged_t',lbound(kmerged_t),ubound(kmerged_t))
      if(allocated(kmerged_u))call netcdfrestart_read(kmerged_u,'kmerged_u',lbound(kmerged_u),ubound(kmerged_u))
      if(allocated(kmerged_v))call netcdfrestart_read(kmerged_v,'kmerged_v',lbound(kmerged_v),ubound(kmerged_v))
      if(allocated(ksl_t))call netcdfrestart_read(ksl_t,'ksl_t',lbound(ksl_t),ubound(ksl_t))
      if(allocated(glob_mask_mangrove))call netcdfrestart_read(glob_mask_mangrove,'glob_mask_mangrove',lbound(glob_mask_mangrove),ubound(glob_mask_mangrove))
      if(allocated(mask_mangrove_t))call netcdfrestart_read(mask_mangrove_t,'mask_mangrove_t',lbound(mask_mangrove_t),ubound(mask_mangrove_t))
      if(allocated(mask_wave_t))call netcdfrestart_read(mask_wave_t,'mask_wave_t',lbound(mask_wave_t),ubound(mask_wave_t))
      if(allocated(upwindriver_t))call netcdfrestart_read(upwindriver_t,'upwindriver_t',lbound(upwindriver_t),ubound(upwindriver_t))
      if(allocated(upwindwetdry_t))call netcdfrestart_read(upwindwetdry_t,'upwindwetdry_t',lbound(upwindwetdry_t),ubound(upwindwetdry_t))
      if(allocated(kmergedr4_u))call netcdfrestart_read(kmergedr4_u,'kmergedr4_u',lbound(kmergedr4_u),ubound(kmergedr4_u))
      if(allocated(kmergedr4_v))call netcdfrestart_read(kmergedr4_v,'kmergedr4_v',lbound(kmergedr4_v),ubound(kmergedr4_v))
      if(allocated(dsigmerged_u))call netcdfrestart_read(dsigmerged_u,'dsigmerged_u',lbound(dsigmerged_u),ubound(dsigmerged_u))
      if(allocated(dsigmerged_v))call netcdfrestart_read(dsigmerged_v,'dsigmerged_v',lbound(dsigmerged_v),ubound(dsigmerged_v))
      if(allocated(pgfratio_u))call netcdfrestart_read(pgfratio_u,'pgfratio_u',lbound(pgfratio_u),ubound(pgfratio_u))
      if(allocated(pgfratio_v))call netcdfrestart_read(pgfratio_v,'pgfratio_v',lbound(pgfratio_v),ubound(pgfratio_v))
      if(allocated(sshr4_w))call netcdfrestart_read(sshr4_w,'sshr4_w',lbound(sshr4_w),ubound(sshr4_w))
      if(allocated(wetmask_u))call netcdfrestart_read(wetmask_u,'wetmask_u',lbound(wetmask_u),ubound(wetmask_u))
      if(allocated(wetmask_v))call netcdfrestart_read(wetmask_v,'wetmask_v',lbound(wetmask_v),ubound(wetmask_v))
      if(allocated(wetmask_t))call netcdfrestart_read(wetmask_t,'wetmask_t',lbound(wetmask_t),ubound(wetmask_t))
      if(allocated(botlevmerged_w))call netcdfrestart_read(botlevmerged_w,'botlevmerged_w',lbound(botlevmerged_w),ubound(botlevmerged_w))
      if(allocated(maxbotstress_aft_w))call netcdfrestart_read(maxbotstress_aft_w,'maxbotstress_aft_w',lbound(maxbotstress_aft_w),ubound(maxbotstress_aft_w))
      if(allocated(maxbotstress_bef_w))call netcdfrestart_read(maxbotstress_bef_w,'maxbotstress_bef_w',lbound(maxbotstress_bef_w),ubound(maxbotstress_bef_w))
      if(allocated(maxbotstress_w))call netcdfrestart_read(maxbotstress_w,'maxbotstress_w',lbound(maxbotstress_w),ubound(maxbotstress_w))
      if(allocated(stresswave_w))call netcdfrestart_read(stresswave_w,'stresswave_w',lbound(stresswave_w),ubound(stresswave_w))
      if(allocated(stressc_w))call netcdfrestart_read(stressc_w,'stressc_w',lbound(stressc_w),ubound(stressc_w))
      if(allocated(sqr_hoverg_u))call netcdfrestart_read(sqr_hoverg_u,'sqr_hoverg_u',lbound(sqr_hoverg_u),ubound(sqr_hoverg_u))
      if(allocated(sqr_hoverg_v))call netcdfrestart_read(sqr_hoverg_v,'sqr_hoverg_v',lbound(sqr_hoverg_v),ubound(sqr_hoverg_v))
      if(allocated(temobc_t))call netcdfrestart_read(temobc_t,'temobc_t',lbound(temobc_t),ubound(temobc_t))
      if(allocated(salobc_t))call netcdfrestart_read(salobc_t,'salobc_t',lbound(salobc_t),ubound(salobc_t))
      if(allocated(velobc_u))call netcdfrestart_read(velobc_u,'velobc_u',lbound(velobc_u),ubound(velobc_u))
      if(allocated(velobc_v))call netcdfrestart_read(velobc_v,'velobc_v',lbound(velobc_v),ubound(velobc_v))
      if(ioffline_prv>=1.and.allocated(temofl_t))call netcdfrestart_read(temofl_t,'temofl_t',lbound(temofl_t),ubound(temofl_t))
      if(ioffline_prv>=1.and.allocated(salofl_t))call netcdfrestart_read(salofl_t,'salofl_t',lbound(salofl_t),ubound(salofl_t))
      if(ioffline_prv>=1.and.allocated(dzofl_t))call netcdfrestart_read(dzofl_t,'dzofl_t',lbound(dzofl_t),ubound(dzofl_t))
      if(ioffline_prv>=1.and.allocated(velofl_u))call netcdfrestart_read(velofl_u,'velofl_u',lbound(velofl_u),ubound(velofl_u))
      if(ioffline_prv>=1.and.allocated(velofl_v))call netcdfrestart_read(velofl_v,'velofl_v',lbound(velofl_v),ubound(velofl_v))
      if(ioffline_prv>=1.and.allocated(tkeofl_w))call netcdfrestart_read(tkeofl_w,'tkeofl_w',lbound(tkeofl_w),ubound(tkeofl_w))
      if(ioffline_prv>=1.and.allocated(bioofl_t))call netcdfrestart_read(bioofl_t,'bioofl_t',lbound(bioofl_t),ubound(bioofl_t))
      if(ioffline_prv>=1.and.allocated(dfvofl_w))call netcdfrestart_read(dfvofl_w,'dfvofl_w',lbound(dfvofl_w),ubound(dfvofl_w))
      if(allocated(uwindabl_t))call netcdfrestart_read(uwindabl_t,'uwindabl_t',lbound(uwindabl_t),ubound(uwindabl_t))
      if(allocated(vwindabl_t))call netcdfrestart_read(vwindabl_t,'vwindabl_t',lbound(vwindabl_t),ubound(vwindabl_t))
      if(ioffline_prv>=1.and.allocated(w0mofl_w))call netcdfrestart_read(w0mofl_w,'w0mofl_w',lbound(w0mofl_w),ubound(w0mofl_w))
      if(ioffline_prv>=1.and.allocated(w_keq1_ofl_w))call netcdfrestart_read(w_keq1_ofl_w,'w_keq1_ofl_w',lbound(w_keq1_ofl_w),ubound(w_keq1_ofl_w))
      if(ioffline_prv>=1.and.allocated(kslofl_t))call netcdfrestart_read(kslofl_t,'kslofl_t',lbound(kslofl_t),ubound(kslofl_t))
      if(allocated(ablheight_t))call netcdfrestart_read(ablheight_t,'ablheight_t',lbound(ablheight_t),ubound(ablheight_t))
      if(allocated(wwindabl_w))call netcdfrestart_read(wwindabl_w,'wwindabl_w',lbound(wwindabl_w),ubound(wwindabl_w))
      if(allocated(kz_abl_w))call netcdfrestart_read(kz_abl_w,'kz_abl_w',lbound(kz_abl_w),ubound(kz_abl_w))
      if(allocated(upwzone0_t))call netcdfrestart_read(upwzone0_t,'upwzone0_t',lbound(upwzone0_t),ubound(upwzone0_t))
      if(allocated(velstokes_u))call netcdfrestart_read(velstokes_u,'velstokes_u',lbound(velstokes_u),ubound(velstokes_u))
      if(allocated(velstokes_v))call netcdfrestart_read(velstokes_v,'velstokes_v',lbound(velstokes_v),ubound(velstokes_v))
      if(allocated(velbarstokes_u))call netcdfrestart_read(velbarstokes_u,'velbarstokes_u',lbound(velbarstokes_u),ubound(velbarstokes_u))
      if(allocated(velbarstokes_v))call netcdfrestart_read(velbarstokes_v,'velbarstokes_v',lbound(velbarstokes_v),ubound(velbarstokes_v))
      if(allocated(nhp1_t))call netcdfrestart_read(nhp1_t,'nhp1_t',lbound(nhp1_t),ubound(nhp1_t))
      if(allocated(nhp2_t))call netcdfrestart_read(nhp2_t,'nhp2_t',lbound(nhp2_t),ubound(nhp2_t))
      if(allocated(temf_t))call netcdfrestart_read(temf_t,'temf_t',lbound(temf_t),ubound(temf_t))
      if(allocated(salf_t))call netcdfrestart_read(salf_t,'salf_t',lbound(salf_t),ubound(salf_t))
      if(allocated(temlwf_t))call netcdfrestart_read(temlwf_t,'temlwf_t',lbound(temlwf_t),ubound(temlwf_t))
      if(allocated(sallwf_t))call netcdfrestart_read(sallwf_t,'sallwf_t',lbound(sallwf_t),ubound(sallwf_t))
      if(allocated(sshlwf_w))call netcdfrestart_read(sshlwf_w,'sshlwf_w',lbound(sshlwf_w),ubound(sshlwf_w))
      if(allocated(t_wave_t))call netcdfrestart_read(t_wave_t,'t_wave_t',lbound(t_wave_t),ubound(t_wave_t))
      if(allocated(hs_wave_t))call netcdfrestart_read(hs_wave_t,'hs_wave_t',lbound(hs_wave_t),ubound(hs_wave_t))
      if(allocated(hsw_wave_t))call netcdfrestart_read(hsw_wave_t,'hsw_wave_t',lbound(hsw_wave_t),ubound(hsw_wave_t))
      if(allocated(foc_wave_t))call netcdfrestart_read(foc_wave_t,'foc_wave_t',lbound(foc_wave_t),ubound(foc_wave_t))
      if(allocated(k_wave_t))call netcdfrestart_read(k_wave_t,'k_wave_t',lbound(k_wave_t),ubound(k_wave_t))
      if(allocated(kx_wave_t))call netcdfrestart_read(kx_wave_t,'kx_wave_t',lbound(kx_wave_t),ubound(kx_wave_t))
      if(allocated(ky_wave_t))call netcdfrestart_read(ky_wave_t,'ky_wave_t',lbound(ky_wave_t),ubound(ky_wave_t))
      if(allocated(twox_wave_t))call netcdfrestart_read(twox_wave_t,'twox_wave_t',lbound(twox_wave_t),ubound(twox_wave_t))
      if(allocated(twoy_wave_t))call netcdfrestart_read(twoy_wave_t,'twoy_wave_t',lbound(twoy_wave_t),ubound(twoy_wave_t))
      if(allocated(tawx_wave_t))call netcdfrestart_read(tawx_wave_t,'tawx_wave_t',lbound(tawx_wave_t),ubound(tawx_wave_t))
      if(allocated(tawy_wave_t))call netcdfrestart_read(tawy_wave_t,'tawy_wave_t',lbound(tawy_wave_t),ubound(tawy_wave_t))
      if(allocated(usf_wave_t))call netcdfrestart_read(usf_wave_t,'usf_wave_t',lbound(usf_wave_t),ubound(usf_wave_t))
      if(allocated(vsf_wave_t))call netcdfrestart_read(vsf_wave_t,'vsf_wave_t',lbound(vsf_wave_t),ubound(vsf_wave_t))
      if(allocated(dir_wave_t))call netcdfrestart_read(dir_wave_t,'dir_wave_t',lbound(dir_wave_t),ubound(dir_wave_t))
      if(allocated(uss_wave_t))call netcdfrestart_read(uss_wave_t,'uss_wave_t',lbound(uss_wave_t),ubound(uss_wave_t))
      if(allocated(j_wave_t))call netcdfrestart_read(j_wave_t,'j_wave_t',lbound(j_wave_t),ubound(j_wave_t))
      if(allocated(vss_wave_t))call netcdfrestart_read(vss_wave_t,'vss_wave_t',lbound(vss_wave_t),ubound(vss_wave_t))
      if(allocated(ubw))call netcdfrestart_read(ubw,'ubw',lbound(ubw),ubound(ubw))
      if(allocated(fw))call netcdfrestart_read(fw,'fw',lbound(fw),ubound(fw))
      if(allocated(dpt_wave_t))call netcdfrestart_read(dpt_wave_t,'dpt_wave_t',lbound(dpt_wave_t),ubound(dpt_wave_t))
      if(allocated(wstresb_u))call netcdfrestart_read(wstresb_u,'wstresb_u',lbound(wstresb_u),ubound(wstresb_u))
      if(allocated(wstresb_v))call netcdfrestart_read(wstresb_v,'wstresb_v',lbound(wstresb_v),ubound(wstresb_v))
      if(allocated(ij2ww3_i))call netcdfrestart_read(ij2ww3_i,'ij2ww3_i',lbound(ij2ww3_i),ubound(ij2ww3_i))
      if(allocated(ij2ww3_j))call netcdfrestart_read(ij2ww3_j,'ij2ww3_j',lbound(ij2ww3_j),ubound(ij2ww3_j))
      if(allocated(ij2ww3_teta))call netcdfrestart_read(ij2ww3_teta,'ij2ww3_teta',lbound(ij2ww3_teta),ubound(ij2ww3_teta))
      if(allocated(slhf_aver_w))call netcdfrestart_read(slhf_aver_w,'slhf_aver_w',lbound(slhf_aver_w),ubound(slhf_aver_w))
      if(allocated(sshf_aver_w))call netcdfrestart_read(sshf_aver_w,'sshf_aver_w',lbound(sshf_aver_w),ubound(sshf_aver_w))
      if(allocated(snsf_aver_w))call netcdfrestart_read(snsf_aver_w,'snsf_aver_w',lbound(snsf_aver_w),ubound(snsf_aver_w))
      if(allocated(ssr_aver_w))call netcdfrestart_read(ssr_aver_w,'ssr_aver_w',lbound(ssr_aver_w),ubound(ssr_aver_w))
      if(allocated(precipi_aver_w))call netcdfrestart_read(precipi_aver_w,'precipi_aver_w',lbound(precipi_aver_w),ubound(precipi_aver_w))
      if(allocated(wstress_aver_u))call netcdfrestart_read(wstress_aver_u,'wstress_aver_u',lbound(wstress_aver_u),ubound(wstress_aver_u))
      if(allocated(wstress_aver_v))call netcdfrestart_read(wstress_aver_v,'wstress_aver_v',lbound(wstress_aver_v),ubound(wstress_aver_v))
      if(allocated(hsedofl_t))call netcdfrestart_read(hsedofl_t,'hsedofl_t',lbound(hsedofl_t),ubound(hsedofl_t))
      read(9)zone1saltflux_w
      read(9)zone1tempflux_w
      read(9)zone1waterflux_w
      read(9)zone1_nlayer
      read(9)zone1_max
      read(9)zone1_u_max
      read(9)zone1_v_max
      read(9)zone1_inv_dz
      read(9)zone1_stretch_dz
      if(allocated(zone1saltflux_glb))call netcdfrestart_read(zone1saltflux_glb,'zone1saltflux_glb',lbound(zone1saltflux_glb),ubound(zone1saltflux_glb))
      if(allocated(zone1saltflux_u))call netcdfrestart_read(zone1saltflux_u,'zone1saltflux_u',lbound(zone1saltflux_u),ubound(zone1saltflux_u))
      if(allocated(zone1saltflux_v))call netcdfrestart_read(zone1saltflux_v,'zone1saltflux_v',lbound(zone1saltflux_v),ubound(zone1saltflux_v))
      if(allocated(zone1tempflux_glb))call netcdfrestart_read(zone1tempflux_glb,'zone1tempflux_glb',lbound(zone1tempflux_glb),ubound(zone1tempflux_glb))
      if(allocated(zone1tempflux_u))call netcdfrestart_read(zone1tempflux_u,'zone1tempflux_u',lbound(zone1tempflux_u),ubound(zone1tempflux_u))
      if(allocated(zone1tempflux_v))call netcdfrestart_read(zone1tempflux_v,'zone1tempflux_v',lbound(zone1tempflux_v),ubound(zone1tempflux_v))
      if(allocated(zone1waterflux_glb))call netcdfrestart_read(zone1waterflux_glb,'zone1waterflux_glb',lbound(zone1waterflux_glb),ubound(zone1waterflux_glb))
      if(allocated(zone1waterflux_u))call netcdfrestart_read(zone1waterflux_u,'zone1waterflux_u',lbound(zone1waterflux_u),ubound(zone1waterflux_u))
      if(allocated(zone1waterflux_v))call netcdfrestart_read(zone1waterflux_v,'zone1waterflux_v',lbound(zone1waterflux_v),ubound(zone1waterflux_v))
      if(allocated(zone1saltflux_glb_in))call netcdfrestart_read(zone1saltflux_glb_in,'zone1saltflux_glb_in',lbound(zone1saltflux_glb_in),ubound(zone1saltflux_glb_in))
      if(allocated(zone1saltflux_u_in))call netcdfrestart_read(zone1saltflux_u_in,'zone1saltflux_u_in',lbound(zone1saltflux_u_in),ubound(zone1saltflux_u_in))
      if(allocated(zone1saltflux_v_in))call netcdfrestart_read(zone1saltflux_v_in,'zone1saltflux_v_in',lbound(zone1saltflux_v_in),ubound(zone1saltflux_v_in))
      if(allocated(zone1tempflux_glb_in))call netcdfrestart_read(zone1tempflux_glb_in,'zone1tempflux_glb_in',lbound(zone1tempflux_glb_in),ubound(zone1tempflux_glb_in))
      if(allocated(zone1tempflux_u_in))call netcdfrestart_read(zone1tempflux_u_in,'zone1tempflux_u_in',lbound(zone1tempflux_u_in),ubound(zone1tempflux_u_in))
      if(allocated(zone1tempflux_v_in))call netcdfrestart_read(zone1tempflux_v_in,'zone1tempflux_v_in',lbound(zone1tempflux_v_in),ubound(zone1tempflux_v_in))
      if(allocated(zone1waterflux_glb_in))call netcdfrestart_read(zone1waterflux_glb_in,'zone1waterflux_glb_in',lbound(zone1waterflux_glb_in),ubound(zone1waterflux_glb_in))
      if(allocated(zone1waterflux_u_in))call netcdfrestart_read(zone1waterflux_u_in,'zone1waterflux_u_in',lbound(zone1waterflux_u_in),ubound(zone1waterflux_u_in))
      if(allocated(zone1waterflux_v_in))call netcdfrestart_read(zone1waterflux_v_in,'zone1waterflux_v_in',lbound(zone1waterflux_v_in),ubound(zone1waterflux_v_in))
      if(allocated(zone1saltflux_glb_out))call netcdfrestart_read(zone1saltflux_glb_out,'zone1saltflux_glb_out',lbound(zone1saltflux_glb_out),ubound(zone1saltflux_glb_out))
      if(allocated(zone1saltflux_u_out))call netcdfrestart_read(zone1saltflux_u_out,'zone1saltflux_u_out',lbound(zone1saltflux_u_out),ubound(zone1saltflux_u_out))
      if(allocated(zone1saltflux_v_out))call netcdfrestart_read(zone1saltflux_v_out,'zone1saltflux_v_out',lbound(zone1saltflux_v_out),ubound(zone1saltflux_v_out))
      if(allocated(zone1tempflux_glb_out))call netcdfrestart_read(zone1tempflux_glb_out,'zone1tempflux_glb_out',lbound(zone1tempflux_glb_out),ubound(zone1tempflux_glb_out))
      if(allocated(zone1tempflux_u_out))call netcdfrestart_read(zone1tempflux_u_out,'zone1tempflux_u_out',lbound(zone1tempflux_u_out),ubound(zone1tempflux_u_out))
      if(allocated(zone1tempflux_v_out))call netcdfrestart_read(zone1tempflux_v_out,'zone1tempflux_v_out',lbound(zone1tempflux_v_out),ubound(zone1tempflux_v_out))
      if(allocated(zone1waterflux_glb_out))call netcdfrestart_read(zone1waterflux_glb_out,'zone1waterflux_glb_out',lbound(zone1waterflux_glb_out),ubound(zone1waterflux_glb_out))
      if(allocated(zone1waterflux_u_out))call netcdfrestart_read(zone1waterflux_u_out,'zone1waterflux_u_out',lbound(zone1waterflux_u_out),ubound(zone1waterflux_u_out))
      if(allocated(zone1waterflux_v_out))call netcdfrestart_read(zone1waterflux_v_out,'zone1waterflux_v_out',lbound(zone1waterflux_v_out),ubound(zone1waterflux_v_out))
      if(allocated(zone1tempcumul_glb))read(9)zone1tempcumul_glb
      if(allocated(zone1tempcumul_loc))read(9)zone1tempcumul_loc
      if(allocated(zone1tempmasst0))read(9)zone1tempmasst0
      if(allocated(zone1saltcumul_glb))read(9)zone1saltcumul_glb
      if(allocated(zone1saltcumul_loc))read(9)zone1saltcumul_loc
      if(allocated(zone1saltmasst0))read(9)zone1saltmasst0
      if(allocated(zone1watercumul_glb))read(9)zone1watercumul_glb
      if(allocated(zone1watercumul_loc))read(9)zone1watercumul_loc
      if(allocated(zone1watermasst0))read(9)zone1watermasst0
      if(allocated(zone1_mask))call netcdfrestart_read(zone1_mask,'zone1_mask',lbound(zone1_mask),ubound(zone1_mask))
      if(allocated(zone1_flux_u_node))call netcdfrestart_read(zone1_flux_u_node,'zone1_flux_u_node',lbound(zone1_flux_u_node),ubound(zone1_flux_u_node))
      if(allocated(zone1_flux_v_node))call netcdfrestart_read(zone1_flux_v_node,'zone1_flux_v_node',lbound(zone1_flux_v_node),ubound(zone1_flux_v_node))
      read(9)zone2saltflux_w
      read(9)zone2tempflux_w
      read(9)zone2waterflux_w
      read(9)zone2_nlayer
      read(9)zone2_max
      read(9)zone2_u_max
      read(9)zone2_v_max
      read(9)zone2_inv_dz
      read(9)zone2_stretch_dz
      if(allocated(zone2saltflux_glb))call netcdfrestart_read(zone2saltflux_glb,'zone2saltflux_glb',lbound(zone2saltflux_glb),ubound(zone2saltflux_glb))
      if(allocated(zone2saltflux_u))call netcdfrestart_read(zone2saltflux_u,'zone2saltflux_u',lbound(zone2saltflux_u),ubound(zone2saltflux_u))
      if(allocated(zone2saltflux_v))call netcdfrestart_read(zone2saltflux_v,'zone2saltflux_v',lbound(zone2saltflux_v),ubound(zone2saltflux_v))
      if(allocated(zone2tempflux_glb))call netcdfrestart_read(zone2tempflux_glb,'zone2tempflux_glb',lbound(zone2tempflux_glb),ubound(zone2tempflux_glb))
      if(allocated(zone2tempflux_u))call netcdfrestart_read(zone2tempflux_u,'zone2tempflux_u',lbound(zone2tempflux_u),ubound(zone2tempflux_u))
      if(allocated(zone2tempflux_v))call netcdfrestart_read(zone2tempflux_v,'zone2tempflux_v',lbound(zone2tempflux_v),ubound(zone2tempflux_v))
      if(allocated(zone2waterflux_glb))call netcdfrestart_read(zone2waterflux_glb,'zone2waterflux_glb',lbound(zone2waterflux_glb),ubound(zone2waterflux_glb))
      if(allocated(zone2waterflux_u))call netcdfrestart_read(zone2waterflux_u,'zone2waterflux_u',lbound(zone2waterflux_u),ubound(zone2waterflux_u))
      if(allocated(zone2waterflux_v))call netcdfrestart_read(zone2waterflux_v,'zone2waterflux_v',lbound(zone2waterflux_v),ubound(zone2waterflux_v))
      if(allocated(zone2saltflux_glb_in))call netcdfrestart_read(zone2saltflux_glb_in,'zone2saltflux_glb_in',lbound(zone2saltflux_glb_in),ubound(zone2saltflux_glb_in))
      if(allocated(zone2saltflux_u_in))call netcdfrestart_read(zone2saltflux_u_in,'zone2saltflux_u_in',lbound(zone2saltflux_u_in),ubound(zone2saltflux_u_in))
      if(allocated(zone2saltflux_v_in))call netcdfrestart_read(zone2saltflux_v_in,'zone2saltflux_v_in',lbound(zone2saltflux_v_in),ubound(zone2saltflux_v_in))
      if(allocated(zone2tempflux_glb_in))call netcdfrestart_read(zone2tempflux_glb_in,'zone2tempflux_glb_in',lbound(zone2tempflux_glb_in),ubound(zone2tempflux_glb_in))
      if(allocated(zone2tempflux_u_in))call netcdfrestart_read(zone2tempflux_u_in,'zone2tempflux_u_in',lbound(zone2tempflux_u_in),ubound(zone2tempflux_u_in))
      if(allocated(zone2tempflux_v_in))call netcdfrestart_read(zone2tempflux_v_in,'zone2tempflux_v_in',lbound(zone2tempflux_v_in),ubound(zone2tempflux_v_in))
      if(allocated(zone2waterflux_glb_in))call netcdfrestart_read(zone2waterflux_glb_in,'zone2waterflux_glb_in',lbound(zone2waterflux_glb_in),ubound(zone2waterflux_glb_in))
      if(allocated(zone2waterflux_u_in))call netcdfrestart_read(zone2waterflux_u_in,'zone2waterflux_u_in',lbound(zone2waterflux_u_in),ubound(zone2waterflux_u_in))
      if(allocated(zone2waterflux_v_in))call netcdfrestart_read(zone2waterflux_v_in,'zone2waterflux_v_in',lbound(zone2waterflux_v_in),ubound(zone2waterflux_v_in))
      if(allocated(zone2saltflux_glb_out))call netcdfrestart_read(zone2saltflux_glb_out,'zone2saltflux_glb_out',lbound(zone2saltflux_glb_out),ubound(zone2saltflux_glb_out))
      if(allocated(zone2saltflux_u_out))call netcdfrestart_read(zone2saltflux_u_out,'zone2saltflux_u_out',lbound(zone2saltflux_u_out),ubound(zone2saltflux_u_out))
      if(allocated(zone2saltflux_v_out))call netcdfrestart_read(zone2saltflux_v_out,'zone2saltflux_v_out',lbound(zone2saltflux_v_out),ubound(zone2saltflux_v_out))
      if(allocated(zone2tempflux_glb_out))call netcdfrestart_read(zone2tempflux_glb_out,'zone2tempflux_glb_out',lbound(zone2tempflux_glb_out),ubound(zone2tempflux_glb_out))
      if(allocated(zone2tempflux_u_out))call netcdfrestart_read(zone2tempflux_u_out,'zone2tempflux_u_out',lbound(zone2tempflux_u_out),ubound(zone2tempflux_u_out))
      if(allocated(zone2tempflux_v_out))call netcdfrestart_read(zone2tempflux_v_out,'zone2tempflux_v_out',lbound(zone2tempflux_v_out),ubound(zone2tempflux_v_out))
      if(allocated(zone2waterflux_glb_out))call netcdfrestart_read(zone2waterflux_glb_out,'zone2waterflux_glb_out',lbound(zone2waterflux_glb_out),ubound(zone2waterflux_glb_out))
      if(allocated(zone2waterflux_u_out))call netcdfrestart_read(zone2waterflux_u_out,'zone2waterflux_u_out',lbound(zone2waterflux_u_out),ubound(zone2waterflux_u_out))
      if(allocated(zone2waterflux_v_out))call netcdfrestart_read(zone2waterflux_v_out,'zone2waterflux_v_out',lbound(zone2waterflux_v_out),ubound(zone2waterflux_v_out))
      if(allocated(zone2tempcumul_glb))read(9)zone2tempcumul_glb
      if(allocated(zone2tempcumul_loc))read(9)zone2tempcumul_loc
      if(allocated(zone2tempmasst0))read(9)zone2tempmasst0
      if(allocated(zone2saltcumul_glb))read(9)zone2saltcumul_glb
      if(allocated(zone2saltcumul_loc))read(9)zone2saltcumul_loc
      if(allocated(zone2saltmasst0))read(9)zone2saltmasst0
      if(allocated(zone2watercumul_glb))read(9)zone2watercumul_glb
      if(allocated(zone2watercumul_loc))read(9)zone2watercumul_loc
      if(allocated(zone2watermasst0))read(9)zone2watermasst0
      if(allocated(zone2_mask))call netcdfrestart_read(zone2_mask,'zone2_mask',lbound(zone2_mask),ubound(zone2_mask))
      if(allocated(zone2_flux_u_node))call netcdfrestart_read(zone2_flux_u_node,'zone2_flux_u_node',lbound(zone2_flux_u_node),ubound(zone2_flux_u_node))
      if(allocated(zone2_flux_v_node))call netcdfrestart_read(zone2_flux_v_node,'zone2_flux_v_node',lbound(zone2_flux_v_node),ubound(zone2_flux_v_node))
      read(9)zone3saltflux_w
      read(9)zone3tempflux_w
      read(9)zone3waterflux_w
      read(9)zone3_nlayer
      read(9)zone3_max
      read(9)zone3_u_max
      read(9)zone3_v_max
      read(9)zone3_inv_dz
      read(9)zone3_stretch_dz
      if(allocated(zone3saltflux_glb))call netcdfrestart_read(zone3saltflux_glb,'zone3saltflux_glb',lbound(zone3saltflux_glb),ubound(zone3saltflux_glb))
      if(allocated(zone3saltflux_u))call netcdfrestart_read(zone3saltflux_u,'zone3saltflux_u',lbound(zone3saltflux_u),ubound(zone3saltflux_u))
      if(allocated(zone3saltflux_v))call netcdfrestart_read(zone3saltflux_v,'zone3saltflux_v',lbound(zone3saltflux_v),ubound(zone3saltflux_v))
      if(allocated(zone3tempflux_glb))call netcdfrestart_read(zone3tempflux_glb,'zone3tempflux_glb',lbound(zone3tempflux_glb),ubound(zone3tempflux_glb))
      if(allocated(zone3tempflux_u))call netcdfrestart_read(zone3tempflux_u,'zone3tempflux_u',lbound(zone3tempflux_u),ubound(zone3tempflux_u))
      if(allocated(zone3tempflux_v))call netcdfrestart_read(zone3tempflux_v,'zone3tempflux_v',lbound(zone3tempflux_v),ubound(zone3tempflux_v))
      if(allocated(zone3waterflux_glb))call netcdfrestart_read(zone3waterflux_glb,'zone3waterflux_glb',lbound(zone3waterflux_glb),ubound(zone3waterflux_glb))
      if(allocated(zone3waterflux_u))call netcdfrestart_read(zone3waterflux_u,'zone3waterflux_u',lbound(zone3waterflux_u),ubound(zone3waterflux_u))
      if(allocated(zone3waterflux_v))call netcdfrestart_read(zone3waterflux_v,'zone3waterflux_v',lbound(zone3waterflux_v),ubound(zone3waterflux_v))
      if(allocated(zone3saltflux_glb_in))call netcdfrestart_read(zone3saltflux_glb_in,'zone3saltflux_glb_in',lbound(zone3saltflux_glb_in),ubound(zone3saltflux_glb_in))
      if(allocated(zone3saltflux_u_in))call netcdfrestart_read(zone3saltflux_u_in,'zone3saltflux_u_in',lbound(zone3saltflux_u_in),ubound(zone3saltflux_u_in))
      if(allocated(zone3saltflux_v_in))call netcdfrestart_read(zone3saltflux_v_in,'zone3saltflux_v_in',lbound(zone3saltflux_v_in),ubound(zone3saltflux_v_in))
      if(allocated(zone3tempflux_glb_in))call netcdfrestart_read(zone3tempflux_glb_in,'zone3tempflux_glb_in',lbound(zone3tempflux_glb_in),ubound(zone3tempflux_glb_in))
      if(allocated(zone3tempflux_u_in))call netcdfrestart_read(zone3tempflux_u_in,'zone3tempflux_u_in',lbound(zone3tempflux_u_in),ubound(zone3tempflux_u_in))
      if(allocated(zone3tempflux_v_in))call netcdfrestart_read(zone3tempflux_v_in,'zone3tempflux_v_in',lbound(zone3tempflux_v_in),ubound(zone3tempflux_v_in))
      if(allocated(zone3waterflux_glb_in))call netcdfrestart_read(zone3waterflux_glb_in,'zone3waterflux_glb_in',lbound(zone3waterflux_glb_in),ubound(zone3waterflux_glb_in))
      if(allocated(zone3waterflux_u_in))call netcdfrestart_read(zone3waterflux_u_in,'zone3waterflux_u_in',lbound(zone3waterflux_u_in),ubound(zone3waterflux_u_in))
      if(allocated(zone3waterflux_v_in))call netcdfrestart_read(zone3waterflux_v_in,'zone3waterflux_v_in',lbound(zone3waterflux_v_in),ubound(zone3waterflux_v_in))
      if(allocated(zone3saltflux_glb_out))call netcdfrestart_read(zone3saltflux_glb_out,'zone3saltflux_glb_out',lbound(zone3saltflux_glb_out),ubound(zone3saltflux_glb_out))
      if(allocated(zone3saltflux_u_out))call netcdfrestart_read(zone3saltflux_u_out,'zone3saltflux_u_out',lbound(zone3saltflux_u_out),ubound(zone3saltflux_u_out))
      if(allocated(zone3saltflux_v_out))call netcdfrestart_read(zone3saltflux_v_out,'zone3saltflux_v_out',lbound(zone3saltflux_v_out),ubound(zone3saltflux_v_out))
      if(allocated(zone3tempflux_glb_out))call netcdfrestart_read(zone3tempflux_glb_out,'zone3tempflux_glb_out',lbound(zone3tempflux_glb_out),ubound(zone3tempflux_glb_out))
      if(allocated(zone3tempflux_u_out))call netcdfrestart_read(zone3tempflux_u_out,'zone3tempflux_u_out',lbound(zone3tempflux_u_out),ubound(zone3tempflux_u_out))
      if(allocated(zone3tempflux_v_out))call netcdfrestart_read(zone3tempflux_v_out,'zone3tempflux_v_out',lbound(zone3tempflux_v_out),ubound(zone3tempflux_v_out))
      if(allocated(zone3waterflux_glb_out))call netcdfrestart_read(zone3waterflux_glb_out,'zone3waterflux_glb_out',lbound(zone3waterflux_glb_out),ubound(zone3waterflux_glb_out))
      if(allocated(zone3waterflux_u_out))call netcdfrestart_read(zone3waterflux_u_out,'zone3waterflux_u_out',lbound(zone3waterflux_u_out),ubound(zone3waterflux_u_out))
      if(allocated(zone3waterflux_v_out))call netcdfrestart_read(zone3waterflux_v_out,'zone3waterflux_v_out',lbound(zone3waterflux_v_out),ubound(zone3waterflux_v_out))
      if(allocated(zone3tempcumul_glb))read(9)zone3tempcumul_glb
      if(allocated(zone3tempcumul_loc))read(9)zone3tempcumul_loc
      if(allocated(zone3tempmasst0))read(9)zone3tempmasst0
      if(allocated(zone3saltcumul_glb))read(9)zone3saltcumul_glb
      if(allocated(zone3saltcumul_loc))read(9)zone3saltcumul_loc
      if(allocated(zone3saltmasst0))read(9)zone3saltmasst0
      if(allocated(zone3watercumul_glb))read(9)zone3watercumul_glb
      if(allocated(zone3watercumul_loc))read(9)zone3watercumul_loc
      if(allocated(zone3watermasst0))read(9)zone3watermasst0
      if(allocated(zone3_mask))call netcdfrestart_read(zone3_mask,'zone3_mask',lbound(zone3_mask),ubound(zone3_mask))
      if(allocated(zone3_flux_u_node))call netcdfrestart_read(zone3_flux_u_node,'zone3_flux_u_node',lbound(zone3_flux_u_node),ubound(zone3_flux_u_node))
      if(allocated(zone3_flux_v_node))call netcdfrestart_read(zone3_flux_v_node,'zone3_flux_v_node',lbound(zone3_flux_v_node),ubound(zone3_flux_v_node))
      read(9)zone4saltflux_w
      read(9)zone4tempflux_w
      read(9)zone4waterflux_w
      read(9)zone4_nlayer
      read(9)zone4_max
      read(9)zone4_u_max
      read(9)zone4_v_max
      read(9)zone4_inv_dz
      read(9)zone4_stretch_dz
      if(allocated(zone4saltflux_glb))call netcdfrestart_read(zone4saltflux_glb,'zone4saltflux_glb',lbound(zone4saltflux_glb),ubound(zone4saltflux_glb))
      if(allocated(zone4saltflux_u))call netcdfrestart_read(zone4saltflux_u,'zone4saltflux_u',lbound(zone4saltflux_u),ubound(zone4saltflux_u))
      if(allocated(zone4saltflux_v))call netcdfrestart_read(zone4saltflux_v,'zone4saltflux_v',lbound(zone4saltflux_v),ubound(zone4saltflux_v))
      if(allocated(zone4tempflux_glb))call netcdfrestart_read(zone4tempflux_glb,'zone4tempflux_glb',lbound(zone4tempflux_glb),ubound(zone4tempflux_glb))
      if(allocated(zone4tempflux_u))call netcdfrestart_read(zone4tempflux_u,'zone4tempflux_u',lbound(zone4tempflux_u),ubound(zone4tempflux_u))
      if(allocated(zone4tempflux_v))call netcdfrestart_read(zone4tempflux_v,'zone4tempflux_v',lbound(zone4tempflux_v),ubound(zone4tempflux_v))
      if(allocated(zone4waterflux_glb))call netcdfrestart_read(zone4waterflux_glb,'zone4waterflux_glb',lbound(zone4waterflux_glb),ubound(zone4waterflux_glb))
      if(allocated(zone4waterflux_u))call netcdfrestart_read(zone4waterflux_u,'zone4waterflux_u',lbound(zone4waterflux_u),ubound(zone4waterflux_u))
      if(allocated(zone4waterflux_v))call netcdfrestart_read(zone4waterflux_v,'zone4waterflux_v',lbound(zone4waterflux_v),ubound(zone4waterflux_v))
      if(allocated(zone4saltflux_glb_in))call netcdfrestart_read(zone4saltflux_glb_in,'zone4saltflux_glb_in',lbound(zone4saltflux_glb_in),ubound(zone4saltflux_glb_in))
      if(allocated(zone4saltflux_u_in))call netcdfrestart_read(zone4saltflux_u_in,'zone4saltflux_u_in',lbound(zone4saltflux_u_in),ubound(zone4saltflux_u_in))
      if(allocated(zone4saltflux_v_in))call netcdfrestart_read(zone4saltflux_v_in,'zone4saltflux_v_in',lbound(zone4saltflux_v_in),ubound(zone4saltflux_v_in))
      if(allocated(zone4tempflux_glb_in))call netcdfrestart_read(zone4tempflux_glb_in,'zone4tempflux_glb_in',lbound(zone4tempflux_glb_in),ubound(zone4tempflux_glb_in))
      if(allocated(zone4tempflux_u_in))call netcdfrestart_read(zone4tempflux_u_in,'zone4tempflux_u_in',lbound(zone4tempflux_u_in),ubound(zone4tempflux_u_in))
      if(allocated(zone4tempflux_v_in))call netcdfrestart_read(zone4tempflux_v_in,'zone4tempflux_v_in',lbound(zone4tempflux_v_in),ubound(zone4tempflux_v_in))
      if(allocated(zone4waterflux_glb_in))call netcdfrestart_read(zone4waterflux_glb_in,'zone4waterflux_glb_in',lbound(zone4waterflux_glb_in),ubound(zone4waterflux_glb_in))
      if(allocated(zone4waterflux_u_in))call netcdfrestart_read(zone4waterflux_u_in,'zone4waterflux_u_in',lbound(zone4waterflux_u_in),ubound(zone4waterflux_u_in))
      if(allocated(zone4waterflux_v_in))call netcdfrestart_read(zone4waterflux_v_in,'zone4waterflux_v_in',lbound(zone4waterflux_v_in),ubound(zone4waterflux_v_in))
      if(allocated(zone4saltflux_glb_out))call netcdfrestart_read(zone4saltflux_glb_out,'zone4saltflux_glb_out',lbound(zone4saltflux_glb_out),ubound(zone4saltflux_glb_out))
      if(allocated(zone4saltflux_u_out))call netcdfrestart_read(zone4saltflux_u_out,'zone4saltflux_u_out',lbound(zone4saltflux_u_out),ubound(zone4saltflux_u_out))
      if(allocated(zone4saltflux_v_out))call netcdfrestart_read(zone4saltflux_v_out,'zone4saltflux_v_out',lbound(zone4saltflux_v_out),ubound(zone4saltflux_v_out))
      if(allocated(zone4tempflux_glb_out))call netcdfrestart_read(zone4tempflux_glb_out,'zone4tempflux_glb_out',lbound(zone4tempflux_glb_out),ubound(zone4tempflux_glb_out))
      if(allocated(zone4tempflux_u_out))call netcdfrestart_read(zone4tempflux_u_out,'zone4tempflux_u_out',lbound(zone4tempflux_u_out),ubound(zone4tempflux_u_out))
      if(allocated(zone4tempflux_v_out))call netcdfrestart_read(zone4tempflux_v_out,'zone4tempflux_v_out',lbound(zone4tempflux_v_out),ubound(zone4tempflux_v_out))
      if(allocated(zone4waterflux_glb_out))call netcdfrestart_read(zone4waterflux_glb_out,'zone4waterflux_glb_out',lbound(zone4waterflux_glb_out),ubound(zone4waterflux_glb_out))
      if(allocated(zone4waterflux_u_out))call netcdfrestart_read(zone4waterflux_u_out,'zone4waterflux_u_out',lbound(zone4waterflux_u_out),ubound(zone4waterflux_u_out))
      if(allocated(zone4waterflux_v_out))call netcdfrestart_read(zone4waterflux_v_out,'zone4waterflux_v_out',lbound(zone4waterflux_v_out),ubound(zone4waterflux_v_out))
      if(allocated(zone4tempcumul_glb))read(9)zone4tempcumul_glb
      if(allocated(zone4tempcumul_loc))read(9)zone4tempcumul_loc
      if(allocated(zone4tempmasst0))read(9)zone4tempmasst0
      if(allocated(zone4saltcumul_glb))read(9)zone4saltcumul_glb
      if(allocated(zone4saltcumul_loc))read(9)zone4saltcumul_loc
      if(allocated(zone4saltmasst0))read(9)zone4saltmasst0
      if(allocated(zone4watercumul_glb))read(9)zone4watercumul_glb
      if(allocated(zone4watercumul_loc))read(9)zone4watercumul_loc
      if(allocated(zone4watermasst0))read(9)zone4watermasst0
      if(allocated(zone4_mask))call netcdfrestart_read(zone4_mask,'zone4_mask',lbound(zone4_mask),ubound(zone4_mask))
      if(allocated(zone4_flux_u_node))call netcdfrestart_read(zone4_flux_u_node,'zone4_flux_u_node',lbound(zone4_flux_u_node),ubound(zone4_flux_u_node))
      if(allocated(zone4_flux_v_node))call netcdfrestart_read(zone4_flux_v_node,'zone4_flux_v_node',lbound(zone4_flux_v_node),ubound(zone4_flux_v_node))
      read(9)zone5saltflux_w
      read(9)zone5tempflux_w
      read(9)zone5waterflux_w
      read(9)zone5_nlayer
      read(9)zone5_max
      read(9)zone5_u_max
      read(9)zone5_v_max
      read(9)zone5_inv_dz
      read(9)zone5_stretch_dz
      if(allocated(zone5saltflux_glb))call netcdfrestart_read(zone5saltflux_glb,'zone5saltflux_glb',lbound(zone5saltflux_glb),ubound(zone5saltflux_glb))
      if(allocated(zone5saltflux_u))call netcdfrestart_read(zone5saltflux_u,'zone5saltflux_u',lbound(zone5saltflux_u),ubound(zone5saltflux_u))
      if(allocated(zone5saltflux_v))call netcdfrestart_read(zone5saltflux_v,'zone5saltflux_v',lbound(zone5saltflux_v),ubound(zone5saltflux_v))
      if(allocated(zone5tempflux_glb))call netcdfrestart_read(zone5tempflux_glb,'zone5tempflux_glb',lbound(zone5tempflux_glb),ubound(zone5tempflux_glb))
      if(allocated(zone5tempflux_u))call netcdfrestart_read(zone5tempflux_u,'zone5tempflux_u',lbound(zone5tempflux_u),ubound(zone5tempflux_u))
      if(allocated(zone5tempflux_v))call netcdfrestart_read(zone5tempflux_v,'zone5tempflux_v',lbound(zone5tempflux_v),ubound(zone5tempflux_v))
      if(allocated(zone5waterflux_glb))call netcdfrestart_read(zone5waterflux_glb,'zone5waterflux_glb',lbound(zone5waterflux_glb),ubound(zone5waterflux_glb))
      if(allocated(zone5waterflux_u))call netcdfrestart_read(zone5waterflux_u,'zone5waterflux_u',lbound(zone5waterflux_u),ubound(zone5waterflux_u))
      if(allocated(zone5waterflux_v))call netcdfrestart_read(zone5waterflux_v,'zone5waterflux_v',lbound(zone5waterflux_v),ubound(zone5waterflux_v))
      if(allocated(zone5saltflux_glb_in))call netcdfrestart_read(zone5saltflux_glb_in,'zone5saltflux_glb_in',lbound(zone5saltflux_glb_in),ubound(zone5saltflux_glb_in))
      if(allocated(zone5saltflux_u_in))call netcdfrestart_read(zone5saltflux_u_in,'zone5saltflux_u_in',lbound(zone5saltflux_u_in),ubound(zone5saltflux_u_in))
      if(allocated(zone5saltflux_v_in))call netcdfrestart_read(zone5saltflux_v_in,'zone5saltflux_v_in',lbound(zone5saltflux_v_in),ubound(zone5saltflux_v_in))
      if(allocated(zone5tempflux_glb_in))call netcdfrestart_read(zone5tempflux_glb_in,'zone5tempflux_glb_in',lbound(zone5tempflux_glb_in),ubound(zone5tempflux_glb_in))
      if(allocated(zone5tempflux_u_in))call netcdfrestart_read(zone5tempflux_u_in,'zone5tempflux_u_in',lbound(zone5tempflux_u_in),ubound(zone5tempflux_u_in))
      if(allocated(zone5tempflux_v_in))call netcdfrestart_read(zone5tempflux_v_in,'zone5tempflux_v_in',lbound(zone5tempflux_v_in),ubound(zone5tempflux_v_in))
      if(allocated(zone5waterflux_glb_in))call netcdfrestart_read(zone5waterflux_glb_in,'zone5waterflux_glb_in',lbound(zone5waterflux_glb_in),ubound(zone5waterflux_glb_in))
      if(allocated(zone5waterflux_u_in))call netcdfrestart_read(zone5waterflux_u_in,'zone5waterflux_u_in',lbound(zone5waterflux_u_in),ubound(zone5waterflux_u_in))
      if(allocated(zone5waterflux_v_in))call netcdfrestart_read(zone5waterflux_v_in,'zone5waterflux_v_in',lbound(zone5waterflux_v_in),ubound(zone5waterflux_v_in))
      if(allocated(zone5saltflux_glb_out))call netcdfrestart_read(zone5saltflux_glb_out,'zone5saltflux_glb_out',lbound(zone5saltflux_glb_out),ubound(zone5saltflux_glb_out))
      if(allocated(zone5saltflux_u_out))call netcdfrestart_read(zone5saltflux_u_out,'zone5saltflux_u_out',lbound(zone5saltflux_u_out),ubound(zone5saltflux_u_out))
      if(allocated(zone5saltflux_v_out))call netcdfrestart_read(zone5saltflux_v_out,'zone5saltflux_v_out',lbound(zone5saltflux_v_out),ubound(zone5saltflux_v_out))
      if(allocated(zone5tempflux_glb_out))call netcdfrestart_read(zone5tempflux_glb_out,'zone5tempflux_glb_out',lbound(zone5tempflux_glb_out),ubound(zone5tempflux_glb_out))
      if(allocated(zone5tempflux_u_out))call netcdfrestart_read(zone5tempflux_u_out,'zone5tempflux_u_out',lbound(zone5tempflux_u_out),ubound(zone5tempflux_u_out))
      if(allocated(zone5tempflux_v_out))call netcdfrestart_read(zone5tempflux_v_out,'zone5tempflux_v_out',lbound(zone5tempflux_v_out),ubound(zone5tempflux_v_out))
      if(allocated(zone5waterflux_glb_out))call netcdfrestart_read(zone5waterflux_glb_out,'zone5waterflux_glb_out',lbound(zone5waterflux_glb_out),ubound(zone5waterflux_glb_out))
      if(allocated(zone5waterflux_u_out))call netcdfrestart_read(zone5waterflux_u_out,'zone5waterflux_u_out',lbound(zone5waterflux_u_out),ubound(zone5waterflux_u_out))
      if(allocated(zone5waterflux_v_out))call netcdfrestart_read(zone5waterflux_v_out,'zone5waterflux_v_out',lbound(zone5waterflux_v_out),ubound(zone5waterflux_v_out))
      if(allocated(zone5tempcumul_glb))read(9)zone5tempcumul_glb
      if(allocated(zone5tempcumul_loc))read(9)zone5tempcumul_loc
      if(allocated(zone5tempmasst0))read(9)zone5tempmasst0
      if(allocated(zone5saltcumul_glb))read(9)zone5saltcumul_glb
      if(allocated(zone5saltcumul_loc))read(9)zone5saltcumul_loc
      if(allocated(zone5saltmasst0))read(9)zone5saltmasst0
      if(allocated(zone5watercumul_glb))read(9)zone5watercumul_glb
      if(allocated(zone5watercumul_loc))read(9)zone5watercumul_loc
      if(allocated(zone5watermasst0))read(9)zone5watermasst0
      if(allocated(zone5_mask))call netcdfrestart_read(zone5_mask,'zone5_mask',lbound(zone5_mask),ubound(zone5_mask))
      if(allocated(zone5_flux_u_node))call netcdfrestart_read(zone5_flux_u_node,'zone5_flux_u_node',lbound(zone5_flux_u_node),ubound(zone5_flux_u_node))
      if(allocated(zone5_flux_v_node))call netcdfrestart_read(zone5_flux_v_node,'zone5_flux_v_node',lbound(zone5_flux_v_node),ubound(zone5_flux_v_node))
      read(9)zone6saltflux_w
      read(9)zone6tempflux_w
      read(9)zone6waterflux_w
      read(9)zone6_nlayer
      read(9)zone6_max
      read(9)zone6_u_max
      read(9)zone6_v_max
      read(9)zone6_inv_dz
      read(9)zone6_stretch_dz
      if(allocated(zone6saltflux_glb))call netcdfrestart_read(zone6saltflux_glb,'zone6saltflux_glb',lbound(zone6saltflux_glb),ubound(zone6saltflux_glb))
      if(allocated(zone6saltflux_u))call netcdfrestart_read(zone6saltflux_u,'zone6saltflux_u',lbound(zone6saltflux_u),ubound(zone6saltflux_u))
      if(allocated(zone6saltflux_v))call netcdfrestart_read(zone6saltflux_v,'zone6saltflux_v',lbound(zone6saltflux_v),ubound(zone6saltflux_v))
      if(allocated(zone6tempflux_glb))call netcdfrestart_read(zone6tempflux_glb,'zone6tempflux_glb',lbound(zone6tempflux_glb),ubound(zone6tempflux_glb))
      if(allocated(zone6tempflux_u))call netcdfrestart_read(zone6tempflux_u,'zone6tempflux_u',lbound(zone6tempflux_u),ubound(zone6tempflux_u))
      if(allocated(zone6tempflux_v))call netcdfrestart_read(zone6tempflux_v,'zone6tempflux_v',lbound(zone6tempflux_v),ubound(zone6tempflux_v))
      if(allocated(zone6waterflux_glb))call netcdfrestart_read(zone6waterflux_glb,'zone6waterflux_glb',lbound(zone6waterflux_glb),ubound(zone6waterflux_glb))
      if(allocated(zone6waterflux_u))call netcdfrestart_read(zone6waterflux_u,'zone6waterflux_u',lbound(zone6waterflux_u),ubound(zone6waterflux_u))
      if(allocated(zone6waterflux_v))call netcdfrestart_read(zone6waterflux_v,'zone6waterflux_v',lbound(zone6waterflux_v),ubound(zone6waterflux_v))
      if(allocated(zone6saltflux_glb_in))call netcdfrestart_read(zone6saltflux_glb_in,'zone6saltflux_glb_in',lbound(zone6saltflux_glb_in),ubound(zone6saltflux_glb_in))
      if(allocated(zone6saltflux_u_in))call netcdfrestart_read(zone6saltflux_u_in,'zone6saltflux_u_in',lbound(zone6saltflux_u_in),ubound(zone6saltflux_u_in))
      if(allocated(zone6saltflux_v_in))call netcdfrestart_read(zone6saltflux_v_in,'zone6saltflux_v_in',lbound(zone6saltflux_v_in),ubound(zone6saltflux_v_in))
      if(allocated(zone6tempflux_glb_in))call netcdfrestart_read(zone6tempflux_glb_in,'zone6tempflux_glb_in',lbound(zone6tempflux_glb_in),ubound(zone6tempflux_glb_in))
      if(allocated(zone6tempflux_u_in))call netcdfrestart_read(zone6tempflux_u_in,'zone6tempflux_u_in',lbound(zone6tempflux_u_in),ubound(zone6tempflux_u_in))
      if(allocated(zone6tempflux_v_in))call netcdfrestart_read(zone6tempflux_v_in,'zone6tempflux_v_in',lbound(zone6tempflux_v_in),ubound(zone6tempflux_v_in))
      if(allocated(zone6waterflux_glb_in))call netcdfrestart_read(zone6waterflux_glb_in,'zone6waterflux_glb_in',lbound(zone6waterflux_glb_in),ubound(zone6waterflux_glb_in))
      if(allocated(zone6waterflux_u_in))call netcdfrestart_read(zone6waterflux_u_in,'zone6waterflux_u_in',lbound(zone6waterflux_u_in),ubound(zone6waterflux_u_in))
      if(allocated(zone6waterflux_v_in))call netcdfrestart_read(zone6waterflux_v_in,'zone6waterflux_v_in',lbound(zone6waterflux_v_in),ubound(zone6waterflux_v_in))
      if(allocated(zone6saltflux_glb_out))call netcdfrestart_read(zone6saltflux_glb_out,'zone6saltflux_glb_out',lbound(zone6saltflux_glb_out),ubound(zone6saltflux_glb_out))
      if(allocated(zone6saltflux_u_out))call netcdfrestart_read(zone6saltflux_u_out,'zone6saltflux_u_out',lbound(zone6saltflux_u_out),ubound(zone6saltflux_u_out))
      if(allocated(zone6saltflux_v_out))call netcdfrestart_read(zone6saltflux_v_out,'zone6saltflux_v_out',lbound(zone6saltflux_v_out),ubound(zone6saltflux_v_out))
      if(allocated(zone6tempflux_glb_out))call netcdfrestart_read(zone6tempflux_glb_out,'zone6tempflux_glb_out',lbound(zone6tempflux_glb_out),ubound(zone6tempflux_glb_out))
      if(allocated(zone6tempflux_u_out))call netcdfrestart_read(zone6tempflux_u_out,'zone6tempflux_u_out',lbound(zone6tempflux_u_out),ubound(zone6tempflux_u_out))
      if(allocated(zone6tempflux_v_out))call netcdfrestart_read(zone6tempflux_v_out,'zone6tempflux_v_out',lbound(zone6tempflux_v_out),ubound(zone6tempflux_v_out))
      if(allocated(zone6waterflux_glb_out))call netcdfrestart_read(zone6waterflux_glb_out,'zone6waterflux_glb_out',lbound(zone6waterflux_glb_out),ubound(zone6waterflux_glb_out))
      if(allocated(zone6waterflux_u_out))call netcdfrestart_read(zone6waterflux_u_out,'zone6waterflux_u_out',lbound(zone6waterflux_u_out),ubound(zone6waterflux_u_out))
      if(allocated(zone6waterflux_v_out))call netcdfrestart_read(zone6waterflux_v_out,'zone6waterflux_v_out',lbound(zone6waterflux_v_out),ubound(zone6waterflux_v_out))
      if(allocated(zone6tempcumul_glb))read(9)zone6tempcumul_glb
      if(allocated(zone6tempcumul_loc))read(9)zone6tempcumul_loc
      if(allocated(zone6tempmasst0))read(9)zone6tempmasst0
      if(allocated(zone6saltcumul_glb))read(9)zone6saltcumul_glb
      if(allocated(zone6saltcumul_loc))read(9)zone6saltcumul_loc
      if(allocated(zone6saltmasst0))read(9)zone6saltmasst0
      if(allocated(zone6watercumul_glb))read(9)zone6watercumul_glb
      if(allocated(zone6watercumul_loc))read(9)zone6watercumul_loc
      if(allocated(zone6watermasst0))read(9)zone6watermasst0
      if(allocated(zone6_mask))call netcdfrestart_read(zone6_mask,'zone6_mask',lbound(zone6_mask),ubound(zone6_mask))
      if(allocated(zone6_flux_u_node))call netcdfrestart_read(zone6_flux_u_node,'zone6_flux_u_node',lbound(zone6_flux_u_node),ubound(zone6_flux_u_node))
      if(allocated(zone6_flux_v_node))call netcdfrestart_read(zone6_flux_v_node,'zone6_flux_v_node',lbound(zone6_flux_v_node),ubound(zone6_flux_v_node))
      if(allocated(anyvar3d))call netcdfrestart_read(anyvar3d,'anyvar3d',lbound(anyvar3d),ubound(anyvar3d))
      if(allocated(sigma_fric_wu))call netcdfrestart_read(sigma_fric_wu,'sigma_fric_wu',lbound(sigma_fric_wu),ubound(sigma_fric_wu))
      if(allocated(sigma_fric_wv))call netcdfrestart_read(sigma_fric_wv,'sigma_fric_wv',lbound(sigma_fric_wv),ubound(sigma_fric_wv))
      if(allocated(dsig_t))call netcdfrestart_read(dsig_t,'dsig_t',lbound(dsig_t),ubound(dsig_t))
      if(allocated(anyv3dint))call netcdfrestart_read(anyv3dint,'anyv3dint',lbound(anyv3dint),ubound(anyv3dint))
      if(allocated(sshobc_w))call netcdfrestart_read(sshobc_w,'sshobc_w',lbound(sshobc_w),ubound(sshobc_w))
      if(allocated(velbarobc_u))call netcdfrestart_read(velbarobc_u,'velbarobc_u',lbound(velbarobc_u),ubound(velbarobc_u))
      if(allocated(velbarobc_v))call netcdfrestart_read(velbarobc_v,'velbarobc_v',lbound(velbarobc_v),ubound(velbarobc_v))
      if(ioffline_prv>=1.and.allocated(sshofl_w))call netcdfrestart_read(sshofl_w,'sshofl_w',lbound(sshofl_w),ubound(sshofl_w))
      if(ioffline_prv>=1.and.allocated(velbarofl_u))call netcdfrestart_read(velbarofl_u,'velbarofl_u',lbound(velbarofl_u),ubound(velbarofl_u))
      if(ioffline_prv>=1.and.allocated(velbarofl_v))call netcdfrestart_read(velbarofl_v,'velbarofl_v',lbound(velbarofl_v),ubound(velbarofl_v))
      if(allocated(tem_delta_t))call netcdfrestart_read(tem_delta_t,'tem_delta_t',lbound(tem_delta_t),ubound(tem_delta_t))
      if(allocated(sal_delta_t))call netcdfrestart_read(sal_delta_t,'sal_delta_t',lbound(sal_delta_t),ubound(sal_delta_t))
      if(allocated(anyvar2d))call netcdfrestart_read(anyvar2d,'anyvar2d',lbound(anyvar2d),ubound(anyvar2d))
      if(allocated(riverdir))read(9)riverdir
      if(allocated(l_river))read(9)l_river
      if(allocated(river_no))read(9)river_no
      if(allocated(rivertrc_inout))read(9)rivertrc_inout
      if(allocated(rivervel_inout))read(9)rivervel_inout
      if(allocated(riverupwdist))read(9)riverupwdist
      if(allocated(nest_len_in))read(9)nest_len_in
      if(allocated(wsed_explicit))read(9)wsed_explicit
      if(allocated(jpm))read(9)jpm
      if(allocated(iriver))call netcdfrestart_read(iriver,'iriver',lbound(iriver),ubound(iriver))
      if(allocated(jriver))call netcdfrestart_read(jriver,'jriver',lbound(jriver),ubound(jriver))
      if(allocated(rankcoords))call netcdfrestart_read(rankcoords,'rankcoords',lbound(rankcoords),ubound(rankcoords))
      if(allocated(spo_i1_u))read(9)spo_i1_u
      if(allocated(spo_i2_u))read(9)spo_i2_u
      if(allocated(spo_j1_u))read(9)spo_j1_u
      if(allocated(spo_j2_u))read(9)spo_j2_u
      if(allocated(spo_i1_v))read(9)spo_i1_v
      if(allocated(spo_i2_v))read(9)spo_i2_v
      if(allocated(spo_j1_v))read(9)spo_j1_v
      if(allocated(spo_j2_v))read(9)spo_j2_v
      if(allocated(spo_i1_t))read(9)spo_i1_t
      if(allocated(spo_i2_t))read(9)spo_i2_t
      if(allocated(spo_j1_t))read(9)spo_j1_t
      if(allocated(spo_j2_t))read(9)spo_j2_t
      if(allocated(obcstreamf))read(9)obcstreamf
      if(allocated(mask_i_w))read(9)mask_i_w
      if(allocated(mask_i_v))read(9)mask_i_v
      if(allocated(mask_i_u))read(9)mask_i_u
      if(allocated(mask_j_w))read(9)mask_j_w
      if(allocated(mask_j_u))read(9)mask_j_u
      if(allocated(mask_j_v))read(9)mask_j_v
      if(allocated(l2ij_out_u))call netcdfrestart_read(l2ij_out_u,'l2ij_out_u',lbound(l2ij_out_u),ubound(l2ij_out_u))
      if(allocated(l2ij_out_v))call netcdfrestart_read(l2ij_out_v,'l2ij_out_v',lbound(l2ij_out_v),ubound(l2ij_out_v))
      if(allocated(l2ij_out_w))call netcdfrestart_read(l2ij_out_w,'l2ij_out_w',lbound(l2ij_out_w),ubound(l2ij_out_w))
      if(allocated(l2ij_in_u))call netcdfrestart_read(l2ij_in_u,'l2ij_in_u',lbound(l2ij_in_u),ubound(l2ij_in_u))
      if(allocated(l2ij_in_v))call netcdfrestart_read(l2ij_in_v,'l2ij_in_v',lbound(l2ij_in_v),ubound(l2ij_in_v))
      if(allocated(l2ij_in_w))call netcdfrestart_read(l2ij_in_w,'l2ij_in_w',lbound(l2ij_in_w),ubound(l2ij_in_w))
      if(allocated(fillmask_t))call netcdfrestart_read(fillmask_t,'fillmask_t',lbound(fillmask_t),ubound(fillmask_t))
      if(allocated(dayindex))read(9)dayindex
      if(allocated(varid))read(9)varid
      if(allocated(grh_nb))read(9)grh_nb
      if(allocated(vardim))read(9)vardim
      if(allocated(varstart))read(9)varstart
      if(allocated(varcount))read(9)varcount
      if(allocated(departuredate))read(9)departuredate
      if(allocated(datesim))call netcdfrestart_read(datesim,'datesim',lbound(datesim),ubound(datesim))
      if(allocated(dateobc))call netcdfrestart_read(dateobc,'dateobc',lbound(dateobc),ubound(dateobc))
      if(allocated(dateairsea))call netcdfrestart_read(dateairsea,'dateairsea',lbound(dateairsea),ubound(dateairsea))
      if(allocated(dateriver))call netcdfrestart_read(dateriver,'dateriver',lbound(dateriver),ubound(dateriver))
      if(allocated(cgridshift))read(9)cgridshift
      if(allocated(runoff_dt))read(9)runoff_dt
      if(allocated(runoff_w))call netcdfrestart_read(runoff_w,'runoff_w',lbound(runoff_w),ubound(runoff_w))
      if(allocated(frqtide))read(9)frqtide
      if(allocated(v0tide))read(9)v0tide
      if(allocated(ti0tide))read(9)ti0tide
      if(allocated(equitide))read(9)equitide
      if(allocated(tideanaweight))read(9)tideanaweight
      if(allocated(airseafile_prvtime))read(9)airseafile_prvtime
      if(allocated(airseafile_nextime))read(9)airseafile_nextime
      if(allocated(ogcm_readtime_next))read(9)ogcm_readtime_next
      if(allocated(ogcm_readtime_prev))read(9)ogcm_readtime_prev
      if(allocated(ogcm_period_prev))read(9)ogcm_period_prev
      if(allocated(ogcm_period_next))read(9)ogcm_period_next
      if(allocated(timeweightobc))read(9)timeweightobc
      if(allocated(dtobc))read(9)dtobc
      if(allocated(dt_abl_max))read(9)dt_abl_max
      if(allocated(zref_z))read(9)zref_z
      if(allocated(checkxyt))read(9)checkxyt
      if(allocated(checkxyf))read(9)checkxyf
      if(allocated(checkanyv3d))read(9)checkanyv3d
      if(allocated(tab1_code))read(9)tab1_code
      if(allocated(tab2_code))read(9)tab2_code
      if(allocated(tab_code_glb))read(9)tab_code_glb
      if(allocated(ub2))read(9)ub2
      if(allocated(lb2))read(9)lb2
      if(allocated(ub3))read(9)ub3
      if(allocated(lb3))read(9)lb3
      if(allocated(ub4))read(9)ub4
      if(allocated(lb4))read(9)lb4
      if(allocated(nutide))read(9)nutide
      if(allocated(nest_len_out))read(9)nest_len_out
      if(allocated(ind))read(9)ind
      if(allocated(airseafile_prvtrec))read(9)airseafile_prvtrec
      if(allocated(airseafile_nextrec))read(9)airseafile_nextrec
      if(allocated(mpi_neighbor_list))read(9)mpi_neighbor_list
      if(allocated(ogcm_rec_next))read(9)ogcm_rec_next
      if(allocated(ogcm_rec_prev))read(9)ogcm_rec_prev
      if(allocated(obcstatus))read(9)obcstatus
      read(9)drifter_out_sampling
      read(9)obc2dtype
      read(9)coef_diss_mangrove
      read(9)expnum
      read(9)cst_c0cub
      read(9)relativewind
      read(9)offset_sshobc
      read(9)biharm_2dfactor
      read(9)checkr0
      read(9)checkr1
      read(9)checkr2
      read(9)checkr3
      read(9)checkr4
      read(9)sponge_l
      read(9)sponge_dx_critic
      read(9)sponge_dx_width
      read(9)dzsurfmin
      if(allocated(drifter_send_order_nord))read(9)drifter_send_order_nord
      if(allocated(drifter_send_order_sud))read(9)drifter_send_order_sud
      if(allocated(drifter_send_order_est))read(9)drifter_send_order_est
      if(allocated(drifter_send_order_ouest))read(9)drifter_send_order_ouest
      if(allocated(drifter_send_order_out))read(9)drifter_send_order_out
      if(allocated(drifter_send_order_nordest))read(9)drifter_send_order_nordest
      if(allocated(drifter_send_order_nordouest))read(9)drifter_send_order_nordouest
      if(allocated(drifter_send_order_sudest))read(9)drifter_send_order_sudest
      if(allocated(drifter_send_order_sudouest))read(9)drifter_send_order_sudouest
      if(allocated(nd_send_canal))read(9)nd_send_canal
      if(allocated(nd_recv_canal))read(9)nd_recv_canal
      if(allocated(zw_abl))read(9)zw_abl
      if(allocated(drifter_send_order_canal))call netcdfrestart_read(drifter_send_order_canal,'drifter_send_order_canal',lbound(drifter_send_order_canal),ubound(drifter_send_order_canal))
      read(9)drifter_onoff
      if(allocated(river_obctype))read(9)river_obctype
      if(allocated(river_s))read(9)river_s
      if(allocated(river_tmin))read(9)river_tmin
      if(allocated(river_tmax))read(9)river_tmax
      if(allocated(hriver))read(9)hriver
      if(allocated(friver))read(9)friver
      if(allocated(realriver))read(9)realriver
      if(allocated(riverinfo))read(9)riverinfo
      if(allocated(river_timeref))read(9)river_timeref
      if(allocated(daysim))read(9)daysim
      if(allocated(tdate_output))read(9)tdate_output
      if(allocated(obcinfo))read(9)obcinfo
      if(allocated(offlinedt))read(9)offlinedt
      if(allocated(pss_mean))read(9)pss_mean
      if(allocated(nest_dt_in))read(9)nest_dt_in
      if(allocated(nest_dt_out))read(9)nest_dt_out
      if(allocated(wavedt))read(9)wavedt
      if(allocated(obc_hf_dt))read(9)obc_hf_dt
      if(allocated(glob_dte_lp))read(9)glob_dte_lp
      if(allocated(glob_dte_lp_tmp))read(9)glob_dte_lp_tmp
      if(allocated(alongresolriver))call netcdfrestart_read(alongresolriver,'alongresolriver',lbound(alongresolriver),ubound(alongresolriver))
      if(allocated(crossresolriver))call netcdfrestart_read(crossresolriver,'crossresolriver',lbound(crossresolriver),ubound(crossresolriver))
      read(9)ww3_varmax
      read(9)ww3_type_grid
      read(9)type_unstructured
      read(9)type_structured
      read(9)vststep
      read(9)zprofile2d3d
      read(9)timestep_type
      read(9)timestep_leapfrog
      read(9)timestep_forwbckw
      read(9)bulk_core
      read(9)bulk_moon
      read(9)bulk_coare
      read(9)bulk_ecume
      read(9)bulk_scheme
      if(allocated(freq_wave))read(9)freq_wave
      read(9)tideana_spinup
      read(9)tideana_delta
      read(9)ogcm_time_lag
      read(9)albedo_val
      read(9)tfilterfb
      read(9)tfilterlf
      read(9)ktide
      read(9)tideforces
      read(9)tideana_yesno
      read(9)ioffline
      read(9)tideanalysis_count
      read(9)ibl1_advbio
      read(9)ibl2_advbio
      read(9)jbl1_advbio
      read(9)jbl2_advbio
      read(9)ogcm_time_shift
      read(9)ieq1
      read(9)ieqimax
      read(9)jeq1
      read(9)jeqjmax
      read(9)ieq1_jeq1
      read(9)ieq1_jeqjmax
      read(9)ieqimax_jeq1
      read(9)ieqimax_jeqjmax
      read(9)dim_varid
      read(9)albedo_constant
      read(9)albedo_apel1987
      read(9)albedo_br1982
      read(9)albedo_br1986
      read(9)albedo_case
      read(9)il_oasis_time
      read(9)flag_steric_effect
      read(9)flag_remove_secondary_bassins
      read(9)one_kind1
      read(9)signe
      read(9)flag_maxbotstress
      read(9)flag_offline_binary
      read(9)flag_wstressbulk
      read(9)flag_write
      read(9)mangrove_scheme
      read(9)flag_meteo_land_plug
      read(9)flag_meteo_land_plug_wind
      read(9)flag_upwind_obc
      read(9)discard_lonlat_periodicity
      read(9)fplan2_grid
      read(9)flag_ts_effectivedensity
      read(9)flag_ts_quicklim
      read(9)flag_z2dv_outputs
      read(9)flag_ogcmname2date
      read(9)flag_meteoname2date
      read(9)ofl_bio
      read(9)drifter_output_files
      read(9)flag_tide3d_analysis
      read(9)flag_bathy_update
      read(9)loop1
      read(9)loop2
      read(9)loop3
      read(9)loopmaxtke
      read(9)loopmaxbio
      read(9)loopmaxts
      read(9)nairsea
      read(9)bi_onoff
      read(9)nc_or_bin_airsea
      read(9)loop_netcdf
      read(9)count_netcdfvar
      read(9)wavefile_prvtrec
      read(9)wavefile_nextrec
      read(9)wave_cpl_nextrec
      read(9)wave_cpl_period_iter
      read(9)wave_cpl_ww3_sdir
      read(9)wave_cpl_ww3_cdir
      read(9)wave_cpl_ww3_hs
      read(9)wave_cpl_ww3_hsw
      read(9)wave_cpl_ww3_foc
      read(9)wave_cpl_ww3_tw
      read(9)wave_cpl_ww3_tawx
      read(9)wave_cpl_ww3_tawy
      read(9)wave_cpl_ww3_twox
      read(9)wave_cpl_ww3_twoy
      read(9)wave_cpl_ww3_uss
      read(9)wave_cpl_ww3_vss
      read(9)wave_cpl_ww3_msk
      read(9)ofl_rec_now
      read(9)ofl_rec_max
      read(9)flag_kz_enhanced
      read(9)flag_net_ir
      read(9)flag_dt_adjust
      read(9)tide_interpolation
      read(9)tide_flagrotation
      read(9)flag_ssr24avr
      read(9)flag_abl
      read(9)flag_abl2
      read(9)flag_sequoia
      read(9)dimssr24prv
      read(9)flag_meteo_average
      read(9)flag_nemoffline
      read(9)trc_id
      read(9)vel_id
      read(9)ssh_id
      read(9)var_num
      read(9)flag_p0m_filter
      read(9)flag_refstate
      read(9)flag_linearfric
      read(9)linear_coef_mangrove
      read(9)flag_1dv
      read(9)flag_merged_levels
      read(9)flag1_smooth_h_mask
      read(9)flag3_smooth_h_mask
      read(9)flag_0status_option
      read(9)flag_rmnegval
      read(9)timemax
      read(9)dirmax
      read(9)freqmax
      read(9)var_misval
      read(9)un_r8
      read(9)heure
      read(9)suma
      read(9)sumb
      read(9)sum_mi
      read(9)grid_area
      read(9)grid_areaglb
      read(9)grid_volumeglb
      read(9)sumarchive
      read(9)lonmin
      read(9)lonmax
      read(9)latmin
      read(9)latmax
      read(9)wavefile_prvtime
      read(9)wavefile_nextime
      read(9)ofl_period_prev
      read(9)ofl_period_now
      read(9)ofl_period_next
      read(9)ofl_nextrec_time
      read(9)ofl_writime
      read(9)ofl_readtime_next
      read(9)ofl_readtime_prev
      read(9)biobc_nextfiletime
      read(9)biobc_prevfiletime
      read(9)stability_index
      read(9)iteration2d_max_r8
      read(9)iteration2d_upbound
      read(9)coef_linearfric
      read(9)gh
      read(9)gm
      read(9)nn
      read(9)tke2overeps
      read(9)tkeovereps
      read(9)gravoverrho
      read(9)sh
      read(9)sm
      read(9)heatfluxbias
      read(9)check0
      read(9)check1
      read(9)check2
      read(9)check3
      read(9)check4
      read(9)check5
      read(9)check6
      read(9)check7
      read(9)check8
      read(9)zero
      read(9)un
      read(9)cdseuil
      read(9)cdb_2dh
      read(9)small1
      read(9)deci
      read(9)decj
      read(9)deck
      read(9)rap1
      read(9)rap2
      read(9)rapi
      read(9)rapj
      read(9)rapk
      read(9)rap
      read(9)const0
      read(9)const1
      read(9)const2
      read(9)const3
      read(9)const4
      read(9)const5
      read(9)const6
      read(9)const7
      read(9)const8
      read(9)const9
      read(9)uv_10
      read(9)karman
      read(9)stefan
      read(9)z2m
      read(9)z10m
      read(9)pss0
      read(9)boltz
      read(9)planck
      read(9)avogadro
      read(9)ce
      read(9)cen
      read(9)ch
      read(9)chn
      read(9)cd
      read(9)cdn
      read(9)z_2
      read(9)z_10
      read(9)q_0
      read(9)r_0
      read(9)pvs_0
      read(9)psih_10
      read(9)psih_2
      read(9)phim
      read(9)zl_10
      read(9)zl_2
      read(9)falpha
      read(9)fbeta
      read(9)ro
      read(9)cp_air
      read(9)lv
      read(9)psim_10
      read(9)psim_2
      read(9)dte_lp
      read(9)inv_dte_lp
      read(9)inv_dti_lp
      read(9)inv_dti_fw
      read(9)inv_dti_fw_p2
      read(9)dte_fw
      read(9)dti_lpbef
      read(9)dti_lp
      read(9)dti_lpmax
      read(9)dti_lpsub
      read(9)dti_fwsub
      read(9)dti_fwsubio
      read(9)dti_fw
      read(9)dti_bef
      read(9)dti_now
      read(9)dtiratio
      read(9)dt_drf
      read(9)fbtfiltercoef
      read(9)assel0
      read(9)assel1
      read(9)assel2
      read(9)assel3
      read(9)wetdry_cst1
      read(9)wetdry_cst2
      read(9)wetdry_cst3
      read(9)h_inf
      read(9)h_inf_obc
      read(9)h_sup
      read(9)dist
      read(9)dist0
      read(9)dist1
      read(9)dist2
      read(9)dist3
      read(9)t0surf
      read(9)s0surf
      read(9)deg2rad
      read(9)rad2deg
      read(9)zmin
      read(9)zmax
      read(9)coa
      read(9)cob
      read(9)coc
      read(9)sol1
      read(9)sol2
      read(9)discri
      read(9)profm1
      read(9)profp1
      read(9)hmax
      read(9)hstepmin
      read(9)hstepmax
      read(9)grav
      read(9)cfl_sshmax
      read(9)cfl_hsshmax
      read(9)cfl_umax
      read(9)cfl_reduce
      read(9)relax_es
      read(9)relax_ext
      read(9)relax_int
      read(9)relax_ts
      read(9)relax_bpc
      read(9)momentum_input_depth
      read(9)z0s
      read(9)z1
      read(9)z2
      read(9)z3
      read(9)z4
      read(9)z5
      read(9)y0
      read(9)y2
      read(9)y3
      read(9)y4
      read(9)y5
      read(9)x0
      read(9)x1
      read(9)x2
      read(9)x3
      read(9)x4
      read(9)x5
      read(9)x6
      read(9)x7
      read(9)x8
      read(9)x9
      read(9)x10
      read(9)x11
      read(9)x12
      read(9)x13
      read(9)x14
      read(9)x20
      read(9)x21
      read(9)x22
      read(9)x33
      read(9)x44
      read(9)area
      read(9)tem_validmin
      read(9)tem_validmax
      read(9)sal_validmin
      read(9)sal_validmax
      read(9)xmd
      read(9)xmv
      read(9)xrd
      read(9)xrv
      read(9)xcpd
      read(9)xcpv
      read(9)xcl
      read(9)celsius2kelvin
      read(9)xlvtt
      read(9)xestt
      read(9)xgamw
      read(9)xbetaw
      read(9)xalpw
      read(9)zrvsrdm1
      read(9)z0_u
      read(9)qsat_sea_z
      read(9)airdensity
      read(9)sst_kelvin
      read(9)prs_atm_z
      read(9)tem_atm_z
      read(9)exner_atm_z
      read(9)delta_u
      read(9)delta_t
      read(9)delta_q
      read(9)delta_u_n
      read(9)exner_sea_z
      read(9)psifunctt
      read(9)psifunctu
      read(9)z0_q
      read(9)z0_t
      read(9)sst1000hpa_kelvin
      read(9)visa
      read(9)charnock
      read(9)qsat_atm_z
      read(9)ustar_bef
      read(9)qstar_bef
      read(9)tetastar_bef
      read(9)rayonterre
      read(9)northpole_lon
      read(9)northpole_lat
      read(9)southpole_lon
      read(9)southpole_lat
      read(9)phi0
      read(9)longi
      read(9)latit
      read(9)longi0
      read(9)latit0
      read(9)longi1
      read(9)latit1
      read(9)angle0
      read(9)alp_t
      read(9)alp_s
      read(9)pi
      read(9)rho
      read(9)inv_rho
      read(9)rhoair
      read(9)valmax
      read(9)vis
      read(9)tkee1
      read(9)tkee2
      read(9)tkee3
      read(9)tkeg2
      read(9)tkeg3
      read(9)tkeg4
      read(9)tkeg5
      read(9)tkeg6
      read(9)tkeb1
      read(9)tkeb2
      read(9)ctke1
      read(9)ctke2
      read(9)t0
      read(9)s0
      read(9)cp
      read(9)light_kpar1
      read(9)light_att1
      read(9)light_att2
      read(9)light_rat1
      read(9)light_rat2
      read(9)light_att2_val1
      read(9)light_att2_h1
      read(9)light_att2_val2
      read(9)light_att2_h2
      read(9)t0_base
      read(9)s0_base
      read(9)rho_base
      read(9)alp_t_base
      read(9)alp_s_base
      read(9)tfb0
      read(9)meteo_lonmin
      read(9)meteo_latmin
      read(9)meteo_lonmax
      read(9)meteo_latmax
      read(9)meteo_resol
      read(9)meteo_resol_u
      read(9)meteo_resol_v
      read(9)meteo_lonstr
      read(9)meteo_lonend
      read(9)meteo_londlt
      read(9)meteo_latstr
      read(9)meteo_latend
      read(9)meteo_latdlt
      read(9)var_lonmin
      read(9)var_latmin
      read(9)var_lonmax
      read(9)var_latmax
      read(9)ww3_lonmin
      read(9)ww3_latmin
      read(9)ww3_lonmax
      read(9)ww3_latmax
      read(9)ww3_dlon
      read(9)ww3_dlat
      read(9)tide_lonmin
      read(9)tide_latmin
      read(9)tide_lonmax
      read(9)tide_latmax
      read(9)tide_dlon
      read(9)tide_dlat
      read(9)dxb
      read(9)dyb
      read(9)dxa
      read(9)dya
      read(9)hmin
      read(9)h1d
      read(9)dlon
      read(9)dlat
      read(9)epsi
      read(9)lagrange_ssh
      read(9)diffu
      read(9)difnorm
      read(9)lup
      read(9)ldown
      read(9)zup
      read(9)zdown
      read(9)rbase
      read(9)small
      read(9)small3
      read(9)hgesig
      read(9)pgesig
      read(9)windfactor
      read(9)xdtk_out
      read(9)tfond
      read(9)sfond
      read(9)rfond
      read(9)c1streamf
      read(9)c2streamf
      read(9)rampe
      read(9)rampe_wind
      read(9)y1
      read(9)obc_hf_reset
      read(9)rho_0d
      read(9)tem_0d
      read(9)sal_0d
      read(9)rho_tmp
      read(9)cst_adv_hor
      read(9)cst_adv_ver
      read(9)cst_adv_vel
      read(9)ssh_avr_nest_out
      read(9)graph_nextime
      read(9)tidenodal_prev_rdv
      read(9)tidenodal_next_rdv
      read(9)tideana_modulo
      read(9)tideana_nextime
      read(9)cellboxfactor1
      read(9)cellboxfactor2
      read(9)rap_wave
      read(9)ihmax
      read(9)jhmax
      read(9)ihmin
      read(9)jhmin
      read(9)convect_yn
      read(9)nbvstepmin
      read(9)bio_relax_size
      read(9)unit_r4
      read(9)x0_r4
      read(9)x1_r4
      read(9)x2_r4
      read(9)x3_r4
      read(9)x4_r4
      read(9)time_r4
      read(9)discharge
      read(9)filval
      read(9)var_validmin
      read(9)var_validmax
      read(9)vdw_loc
      read(9)vup_loc
      read(9)zdw_loc
      read(9)var_scalefactor
      read(9)inv_scalefactor
      read(9)var_addoffset
      read(9)zup_loc
      read(9)hrmax
      read(9)relax_bio
      read(9)kmol_m
      read(9)kmol_h
      read(9)kmol_s
      read(9)upw_hrange1
      read(9)upw_hrange2
      read(9)inv_ekman_depth
      read(9)constant_km
      read(9)constant_kh
      read(9)z0b
      read(9)z0b_land
      read(9)z0b_rivers
      read(9)zlevel_land
      read(9)invloopmaxts
      read(9)invloopmaxu
      read(9)invloopmaxv
      read(9)sqrtgrav
      read(9)invgrav
      read(9)spinup_forcing
      read(9)relax_lwf
      read(9)dz_vertical_incr_fact
      read(9)ratio_negdif_ver
      read(9)ratio_negdif_hor
      read(9)ratio_bionegdif
      read(9)vqs_cst1
      read(9)vqs_cst2
      read(9)vqs_cst3
      read(9)ema_mu
      read(9)dz_over_z0_min
      read(9)coastal_viscosity
      read(9)quick_coef
      read(9)i
      read(9)j
      read(9)k
      read(9)flag3d
      read(9)lrec
      read(9)iteration3d
      read(9)compt1
      read(9)compt2
      read(9)compt3
      read(9)compt4
      read(9)kount0
      read(9)kount1
      read(9)kount2
      read(9)kount3
      read(9)kount4
      read(9)kount5
      read(9)kount6
      read(9)kount7
      read(9)kount8
      read(9)kount9
      read(9)kountmod
      read(9)kountrdv1
      read(9)kountrdv2
      read(9)kountrdv3
      read(9)kountrdv4
      read(9)substep_advbio
      read(9)subcycle_exchange
      read(9)subcycle_onoff
      read(9)subcycle_synchro
      read(9)subcycle_modulo
      read(9)dt_drf_over_dti_fw
      read(9)iteration3d_restart
      read(9)filvalshort
      read(9)quick_filter_points
      read(9)status
      read(9)forcedstatus
      read(9)decision
      read(9)ncid1
      read(9)ncid2
      read(9)dim_x_id
      read(9)dim_y_id
      read(9)dim_z_id
      read(9)dim_t_id
      read(9)dim_b_id
      read(9)max_x
      read(9)max_y
      read(9)max_z
      read(9)max_time_counter
      read(9)max_meteo_time_counter
      read(9)var_id
      read(9)var_nftype
      read(9)meteo_imax
      read(9)meteo_jmax
      read(9)meteo_kmax
      read(9)meteozoom_istr
      read(9)meteozoom_iend
      read(9)meteozoom_jstr
      read(9)meteozoom_jend
      read(9)meteofull_imax
      read(9)meteofull_jmax
      read(9)tide_imax
      read(9)tide_jmax
      read(9)tide_kmax
      read(9)tidezoom_istr
      read(9)tidezoom_iend
      read(9)tidezoom_jstr
      read(9)tidezoom_jend
      read(9)tidezoom_istr_t
      read(9)tidezoom_iend_t
      read(9)tidezoom_jstr_t
      read(9)tidezoom_jend_t
      read(9)tidezoom_istr_u
      read(9)tidezoom_iend_u
      read(9)tidezoom_jstr_u
      read(9)tidezoom_jend_u
      read(9)tidezoom_istr_v
      read(9)tidezoom_iend_v
      read(9)tidezoom_jstr_v
      read(9)tidezoom_jend_v
      read(9)tidefull_imax
      read(9)tidefull_jmax
      read(9)ww3_imax
      read(9)ww3_jmax
      read(9)ww3_kmax
      read(9)ww3_fmax
      read(9)ww3zoom_istr
      read(9)ww3zoom_iend
      read(9)ww3zoom_jstr
      read(9)ww3zoom_jend
      read(9)ww3full_imax
      read(9)ww3full_jmax
      read(9)jour
      read(9)kstop
      read(9)istr
      read(9)jstr
      read(9)kstr
      read(9)tstr
      read(9)bstr
      read(9)iend
      read(9)jend
      read(9)kend
      read(9)tend
      read(9)bend
      read(9)dimend
      read(9)kpvwave
      read(9)give_chanel9
      read(9)i0
      read(9)j0
      read(9)iteration2d
      read(9)iteration2d_begin
      read(9)iteration2d_max_now
      read(9)iteration2d_max_bef
      read(9)i2dh
      read(9)l1
      read(9)l1_sca
      read(9)l1_vec
      read(9)l2
      read(9)l2_sca
      read(9)l2_vec
      read(9)l3
      read(9)l3_sca
      read(9)l3_vec
      read(9)len1
      read(9)nc
      read(9)nc1
      read(9)itimets
      read(9)itimebio
      read(9)iadvec_ts_hor
      read(9)iadvec_ts_hor_upwind
      read(9)iadvec_ts_hor_quickest
      read(9)iadvec_ts_hor_quickest2
      read(9)iadvec_ts_ver
      read(9)iadvec_ts_ver_quickest
      read(9)iadvec_ts_ver_c2
      read(9)iadvec_ts_ver_quickest2
      read(9)iturbulence
      read(9)istreamf
      read(9)itime
      read(9)kts
      read(9)kuv
      read(9)ko
      read(9)jm1
      read(9)jm2
      read(9)jp1
      read(9)jp2
      read(9)im1
      read(9)im2
      read(9)ip1
      read(9)ip2
      read(9)iwind
      read(9)ip
      read(9)im
      read(9)jp
      read(9)jm
      read(9)kp
      read(9)km
      read(9)nbcanal
      read(9)gridtype1
      read(9)gridtype2
      read(9)point1
      read(9)point2
      read(9)sender
      read(9)receiver
      read(9)ipnoc
      read(9)i1
      read(9)i2
      read(9)i3
      read(9)i4
      read(9)i5
      read(9)i6
      read(9)i7
      read(9)i8
      read(9)i9
      read(9)i10
      read(9)i11
      read(9)j1
      read(9)j2
      read(9)j3
      read(9)j4
      read(9)j5
      read(9)j6
      read(9)j7
      read(9)j8
      read(9)j9
      read(9)j10
      read(9)j11
      read(9)k0
      read(9)k1
      read(9)k2
      read(9)k3
      read(9)k4
      read(9)k5
      read(9)k6
      read(9)k7
      read(9)k8
      read(9)k9
      read(9)kr
      read(9)kp1
      read(9)kp2
      read(9)km1
      read(9)km2
      read(9)key
      read(9)flag
      read(9)nbomax
      read(9)nbobuffermax
      read(9)nbobuffermax_c
      read(9)flag_stop
      read(9)itest
      read(9)itest1
      read(9)itest2
      read(9)itest3
      read(9)ioption
      read(9)nsmooth
      read(9)nriver
      read(9)ian0
      read(9)imois0
      read(9)ijour0
      read(9)iheure0
      read(9)iminute0
      read(9)iseconde0
      read(9)iref
      read(9)jref
      read(9)lname1
      read(9)lname2
      read(9)lname3
      read(9)lname4
      read(9)kmode
      read(9)kmodemax
      read(9)fgrid_or_wgrid
      read(9)fgrid_case
      read(9)wgrid_case
      read(9)typegrid
      read(9)typegrid_monopole
      read(9)typegrid_file
      read(9)typegrid_bipole
      read(9)istar
      read(9)jstar
      read(9)istop
      read(9)jstop
      read(9)isend
      read(9)jsend
      read(9)irecv
      read(9)jrecv
      read(9)izoomin
      read(9)izoomax
      read(9)jzoomin
      read(9)jzoomax
      read(9)iairsea
      read(9)ialbedo
      read(9)airseaoption
      read(9)iobc_f
      read(9)iobc_wv
      read(9)iobc_ogcm
      read(9)iobc_lr
      read(9)obc_option
      read(9)iobc_demo_wv
      read(9)iarchive
      read(9)imodeltrc
      read(9)imodelbio
      read(9)multiple
      read(9)kvarmax
      read(9)igesig
      read(9)isigfile
      read(9)ksecu
      read(9)kmaxtide
      read(9)kmaxtidep1
      read(9)kminserie
      read(9)nzctdmax
      read(9)nsctdmax
      read(9)ktctdmin
      read(9)ktctdmax
      read(9)kdtk_out
      read(9)k10
      read(9)k11
      read(9)nbinco
      read(9)nbequa
      read(9)nbsparse
      read(9)kland
      read(9)tidestep2
      read(9)tidestep3
      read(9)ncfilm_max
      read(9)in_out_tide
      read(9)kmax_dof
      read(9)k_in
      read(9)k_out
      read(9)used_unused_dom
      read(9)mergebathy_sponge
      read(9)rankhmax
      read(9)rankhmin
      read(9)id_tem
      read(9)id_tem2
      read(9)id_dtem
      read(9)id_sal
      read(9)id_sal2
      read(9)id_rhp
      read(9)id_rhop
      read(9)id_rhom
      read(9)id_rhf
      read(9)id_rhc
      read(9)id_rhpa
      read(9)id_rhpb
      read(9)id_z
      read(9)id_prs
      read(9)id_now
      read(9)id_aft
      read(9)id_eost
      read(9)id_eoss
      read(9)id_ssh
      read(9)id_breaker
      read(9)id_dz
      read(9)id_zt
      read(9)id_tdiv
      read(9)id_sdiv
      read(9)id_udiv
      read(9)id_vdiv
      read(9)id_bdiv
      read(9)id_bnegdif
      read(9)id_biobefo
      read(9)id_bioaftr
      read(9)id_u_now
      read(9)id_v_now
      read(9)id_u_rot
      read(9)id_v_rot
      read(9)id_ffreq
      read(9)id_coriolis1
      read(9)id_coriolis2
      read(9)id_flx
      read(9)id_fly
      read(9)id_uflx
      read(9)id_ufly
      read(9)id_vflx
      read(9)id_vfly
      read(9)id_tcn
      read(9)id_scn
      read(9)id_gradssh
      read(9)id_ofactort
      read(9)id_ofactors
      read(9)id_hybcoefu
      read(9)id_hybcoefv
      read(9)id_dxdydz
      read(9)id_wb
      read(9)id_prod
      read(9)id_buoy
      read(9)iDtOvRhCp
      read(9)id_ncu
      read(9)id_ncv
      read(9)id_webio
      read(9)id_varbef2
      read(9)id_kh_over_dz
      read(9)id_bihar_lim
      read(9)id_veltot
      read(9)id_velexp
      read(9)id_velimp
      read(9)id_wdrifter
      read(9)kreadgroup
      read(9)kreadgroupmax
      read(9)kread1
      read(9)kread2
      read(9)modulo_biotimestep
      read(9)obctime_bef2
      read(9)obctime_bef
      read(9)obctime_aft
      read(9)obctime_aft2
      read(9)obctime_order
      read(9)sch_imp_ts_u_loc
      read(9)sch_imp_ts_v_loc
      read(9)sch_imp_ts_u_glb
      read(9)sch_imp_ts_v_glb
      read(9)sch_imp_tke_u_loc
      read(9)sch_imp_tke_v_loc
      read(9)sch_imp_tke_u_glb
      read(9)sch_imp_tke_v_glb
      read(9)looplimit_hor
      read(9)ihybsig
      read(9)nhybsig
      read(9)tke_surf
      read(9)grh_out_mi
      read(9)nest_onoff_in
      read(9)nest_onoff_demo
      read(9)nest_full_in
      read(9)nest_full_out
      read(9)nest_onoff_out
      read(9)ncmin_airsea
      read(9)ncmin_river
      read(9)eos_author
      read(9)eos_comprs
      read(9)eos_linear
      read(9)obcfreeorfix
      read(9)iwve
      read(9)wave_obc_type
      read(9)dataperwavefile
      read(9)sp_or_db
      read(9)kbu
      read(9)kbumax
      read(9)kbumax_glb
      read(9)n_element
      read(9)relaxtype_ts
      read(9)obctype_ts
      read(9)obctype_p
      read(9)restart_file_y_or_n
      read(9)ncid
      read(9)x_xhl_dim
      read(9)y_xhl_dim
      read(9)z_xhl_dim
      read(9)x_yhl_dim
      read(9)y_yhl_dim
      read(9)z_yhl_dim
      read(9)x_zhl_dim
      read(9)y_zhl_dim
      read(9)z_zhl_dim
      read(9)x_zl_dim
      read(9)y_zl_dim
      read(9)z_zl_dim
      read(9)time_dim
      read(9)i_t_dim
      read(9)j_t_dim
      read(9)k_t_dim
      read(9)i_w_dim
      read(9)j_w_dim
      read(9)k_w_dim
      read(9)i_u_dim
      read(9)j_u_dim
      read(9)k_u_dim
      read(9)i_v_dim
      read(9)j_v_dim
      read(9)k_v_dim
      read(9)i_f_dim
      read(9)j_f_dim
      read(9)k_f_dim
      read(9)dayindex_size
      read(9)tide_year_min
      read(9)airsea_year_min
      read(9)tide_year_max
      read(9)airsea_year_max
      read(9)wave_year_min
      read(9)irelaxsst
      read(9)removetide
      read(9)ofl_rotation
      read(9)ofl_rhp
      read(9)ofl_tke
      read(9)ofl_surflux
      read(9)ale_selected
      read(9)flag_asselin
      read(9)freq
      read(9)dirw
      read(9)dirw_beg
      read(9)dirw_end
      read(9)freq_beg
      read(9)freq_end
      read(9)year_now
      read(9)month_now
      read(9)day_now
      read(9)hour_now
      read(9)minute_now
      read(9)second_now
      read(9)nd_send_est
      read(9)nd_send_ouest
      read(9)nd_send_nord
      read(9)nd_send_sud
      read(9)nd_send_out
      read(9)nd_recv_est
      read(9)nd_recv_ouest
      read(9)nd_recv_nord
      read(9)nd_recv_sud
      read(9)nd_send_sudouest
      read(9)nd_send_sudest
      read(9)nd_send_nordouest
      read(9)nd_send_nordest
      read(9)nd_recv_sudouest
      read(9)nd_recv_sudest
      read(9)nd_recv_nordouest
      read(9)nd_recv_nordest
      read(9)initial_main_status
      read(9)offline_init_status
      read(9)meteo_sealand_mask
      read(9)grid_i0
      read(9)grid_j0
      read(9)ifb
      read(9)dbefore
      read(9)dnow
      read(9)dafter
      read(9)dim_airsea
      read(9)rhp_zavr_xy
      read(9)ssr_id
      read(9)ir_id
      read(9)rain_id
      read(9)t2m_id
      read(9)t0m_id
      read(9)abl_id
      read(9)dp2m_id
      read(9)u10m_id
      read(9)v10m_id
      read(9)u100m_id
      read(9)v100m_id
      read(9)p0m_id
      read(9)ustrs_id
      read(9)vstrs_id
      read(9)slhf_id
      read(9)netir_id
      read(9)sshf_id
      read(9)t_flux_cumul
      read(9)s_flux_cumul
      read(9)cumuldeltaflux
      read(9)som0
      read(9)som2
      read(9)ssh_reservoir
      read(9)idealflux
      read(9)sum0
      read(9)sum1
      read(9)sum2
      read(9)sum3
      read(9)sum4
      read(9)sum5
      read(9)sum6
      read(9)sum7
      read(9)sum8
      read(9)sum9
      read(9)sum10
      read(9)sum11
      read(9)sum12
      read(9)sum13
      read(9)time0
      read(9)time1
      read(9)time2
      read(9)small2
      read(9)x1_r8
      read(9)x2_r8
      read(9)x3_r8
      read(9)x4_r8
      read(9)sum0glb
      read(9)sum1glb
      read(9)sum2glb
      read(9)sum3glb
      read(9)sum4glb
      read(9)sum5glb
      read(9)sum6glb
      read(9)emin
      read(9)epsmin
      read(9)elapsedtime_out
      read(9)elapsedtime_rst
      read(9)elapsedtime_now
      read(9)elapsedtime_last_writing
      read(9)elapsedtime_bef
      read(9)elapsedtime_aft
      read(9)cpu_seconds
      read(9)alpha
      read(9)eos_pgfzref
      read(9)eos_tkezref
      read(9)filvalr8
      read(9)dti_fw_i4
      read(9)elapsedtime_now_i4
      read(9)elapsedtime_now_i8
      read(9)elapsedtime_now_r16
! insere partie_lit_1 fin. SURTOUT NE PAS EFFACER CETTE LIGNE

      close(9)
      if(par%rank==0)write(6,*)'lecture chanel9 fin.'
      call dyn_restart_allo('d') ! desalloue des tableaux sous certaines conditions

      if(datesimbis(1,1)/=datesim(1,1).or. &
         datesimbis(2,1)/=datesim(2,1).or. &
         datesimbis(3,1)/=datesim(3,1).or. &
         datesimbis(4,1)/=datesim(4,1).or. &
         datesimbis(5,1)/=datesim(5,1).or. &
         datesimbis(6,1)/=datesim(6,1))then !>>>
        if(par%rank==0) then !www>
         write(6,'(a)')'Inconsistent departure date in notebook_time'
         write(6,'(a,6(i4,1x))')'Replace ',datesimbis(:,1)
         write(6,'(a,6(i4,1x))')'with    ',datesim(:,1) !14-11-18
        endif                !www>
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
      stop 'Err 1163 dyn_restart'
      endif                                 !>>>
      datesim=datesimbis

      deallocate(datesimbis)
      endif                  !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr>

      if (txt_=='w'          &
      .or.txt_=='W') then !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww>

      if(run_option==-1)return !04-06-18 ne pas ecriture de fichier restart si run_option=-1


      debug_=1

      call s_cpu('chanel9_writting_bef',1)

      if(restart_file_y_or_n==0)then
       write(*,*)
       write(*,*)'erreur dans dyn_restart.f: la routine est appellee'
       write(*,*)'alors que la sequence de calcul pour l''iteration'
       write(*,*)'en cours n''est pas terminee.'
       write(*,*)'remplacez votre appel par procedure similaire à'
       write(*,*)'celle utilisee dans date_output.f.'
       stop ' stop dans dyn_restart.f'
      endif

! Check the principal fields and write the restart file only if the check list is OK:
      call dyn_restart_check              !27-05-11

!     if(dynrestartfilename(1:14)==restartdir_out1) then !----> !17-10-13
!        dynrestartfilename(1:14)= restartdir_out2
!     else                                                !---->
!        dynrestartfilename(1:14)= restartdir_out1
!     endif                                               !---->

      k1=0 ; k2=0
      k1=index(dynrestartfilename,trim(restartdir_out1))
      k2=index(dynrestartfilename,trim(restartdir_out2))
      if(k1/=0)dynrestartfilename(k1:k1+len(trim(restartdir_out2))-1)=trim(restartdir_out2)
      if(k2/=0)dynrestartfilename(k2:k2+len(trim(restartdir_out1))-1)=trim(restartdir_out1)

      call dyn_restart_bounds_w

      open(unit=9,file=dynrestartfilename                      &     !12-05-11
       ,access='sequential',form='unformatted')

      if(par%rank==0) then !00000>
       write(6,'(a12,a50)')'ecriture de ',dynrestartfilename          !12-05-11
       write(6,*)'debut:'
      endif                !00000>

! insere partie_ecrit_1 debut: SURTOUT NE PAS EFFACER CETTE LIGNE
      write(9)ioffline_prv
      write(9)zone1bioflux_w
      if(allocated(zone1bioflux_glb))call netcdfrestart_wrt(zone1bioflux_glb,'zone1bioflux_glb',lbound(zone1bioflux_glb),ubound(zone1bioflux_glb))
      if(allocated(zone1bioflux_u))call netcdfrestart_wrt(zone1bioflux_u,'zone1bioflux_u',lbound(zone1bioflux_u),ubound(zone1bioflux_u))
      if(allocated(zone1bioflux_v))call netcdfrestart_wrt(zone1bioflux_v,'zone1bioflux_v',lbound(zone1bioflux_v),ubound(zone1bioflux_v))
      if(allocated(zone1tendancebio_glb))call netcdfrestart_wrt(zone1tendancebio_glb,'zone1tendancebio_glb',lbound(zone1tendancebio_glb),ubound(zone1tendancebio_glb))
      if(allocated(zone1botsurfbio_glb))call netcdfrestart_wrt(zone1botsurfbio_glb,'zone1botsurfbio_glb',lbound(zone1botsurfbio_glb),ubound(zone1botsurfbio_glb))
      if(allocated(zone1bioflux_glb_in))call netcdfrestart_wrt(zone1bioflux_glb_in,'zone1bioflux_glb_in',lbound(zone1bioflux_glb_in),ubound(zone1bioflux_glb_in))
      if(allocated(zone1bioflux_u_in))call netcdfrestart_wrt(zone1bioflux_u_in,'zone1bioflux_u_in',lbound(zone1bioflux_u_in),ubound(zone1bioflux_u_in))
      if(allocated(zone1bioflux_v_in))call netcdfrestart_wrt(zone1bioflux_v_in,'zone1bioflux_v_in',lbound(zone1bioflux_v_in),ubound(zone1bioflux_v_in))
      if(allocated(zone1bioflux_glb_out))call netcdfrestart_wrt(zone1bioflux_glb_out,'zone1bioflux_glb_out',lbound(zone1bioflux_glb_out),ubound(zone1bioflux_glb_out))
      if(allocated(zone1bioflux_u_out))call netcdfrestart_wrt(zone1bioflux_u_out,'zone1bioflux_u_out',lbound(zone1bioflux_u_out),ubound(zone1bioflux_u_out))
      if(allocated(zone1bioflux_v_out))call netcdfrestart_wrt(zone1bioflux_v_out,'zone1bioflux_v_out',lbound(zone1bioflux_v_out),ubound(zone1bioflux_v_out))
      if(allocated(zone1biomasst0))call netcdfrestart_wrt(zone1biomasst0,'zone1biomasst0',lbound(zone1biomasst0),ubound(zone1biomasst0))
      if(allocated(zone1biocumul_glb))write(9)zone1biocumul_glb
      if(allocated(zone1biocumul_loc))write(9)zone1biocumul_loc
      if(allocated(analysis3dmatrix))call netcdfrestart_wrt(analysis3dmatrix,'analysis3dmatrix',lbound(analysis3dmatrix),ubound(analysis3dmatrix))
      if(allocated(ema1_s_i))call netcdfrestart_wrt(ema1_s_i,'ema1_s_i',lbound(ema1_s_i),ubound(ema1_s_i))
      if(allocated(ema2_s_i))call netcdfrestart_wrt(ema2_s_i,'ema2_s_i',lbound(ema2_s_i),ubound(ema2_s_i))
      if(allocated(ema1_s_j))call netcdfrestart_wrt(ema1_s_j,'ema1_s_j',lbound(ema1_s_j),ubound(ema1_s_j))
      if(allocated(ema2_s_j))call netcdfrestart_wrt(ema2_s_j,'ema2_s_j',lbound(ema2_s_j),ubound(ema2_s_j))
      if(allocated(c_wave_mode))call netcdfrestart_wrt(c_wave_mode,'c_wave_mode',lbound(c_wave_mode),ubound(c_wave_mode))
      if(allocated(pcoefmode_t))call netcdfrestart_wrt(pcoefmode_t,'pcoefmode_t',lbound(pcoefmode_t),ubound(pcoefmode_t))
      if(allocated(ucoefmode_t))call netcdfrestart_wrt(ucoefmode_t,'ucoefmode_t',lbound(ucoefmode_t),ubound(ucoefmode_t))
      if(allocated(vcoefmode_t))call netcdfrestart_wrt(vcoefmode_t,'vcoefmode_t',lbound(vcoefmode_t),ubound(vcoefmode_t))
      if(allocated(rhrefmode_t))call netcdfrestart_wrt(rhrefmode_t,'rhrefmode_t',lbound(rhrefmode_t),ubound(rhrefmode_t))
      if(allocated(ema1_q_i))call netcdfrestart_wrt(ema1_q_i,'ema1_q_i',lbound(ema1_q_i),ubound(ema1_q_i))
      if(allocated(ema2_q_i))call netcdfrestart_wrt(ema2_q_i,'ema2_q_i',lbound(ema2_q_i),ubound(ema2_q_i))
      if(allocated(ema1_q_j))call netcdfrestart_wrt(ema1_q_j,'ema1_q_j',lbound(ema1_q_j),ubound(ema1_q_j))
      if(allocated(ema2_q_j))call netcdfrestart_wrt(ema2_q_j,'ema2_q_j',lbound(ema2_q_j),ubound(ema2_q_j))
      if(allocated(uv_wmode_t))call netcdfrestart_wrt(uv_wmode_t,'uv_wmode_t',lbound(uv_wmode_t),ubound(uv_wmode_t))
      if(allocated(analysis_p3d_t))call netcdfrestart_wrt(analysis_p3d_t,'analysis_p3d_t',lbound(analysis_p3d_t),ubound(analysis_p3d_t))
      if(allocated(analysis_u3d_t))call netcdfrestart_wrt(analysis_u3d_t,'analysis_u3d_t',lbound(analysis_u3d_t),ubound(analysis_u3d_t))
      if(allocated(analysis_v3d_t))call netcdfrestart_wrt(analysis_v3d_t,'analysis_v3d_t',lbound(analysis_v3d_t),ubound(analysis_v3d_t))
      if(allocated(p3dmode_cos_t))call netcdfrestart_wrt(p3dmode_cos_t,'p3dmode_cos_t',lbound(p3dmode_cos_t),ubound(p3dmode_cos_t))
      if(allocated(p3dmode_sin_t))call netcdfrestart_wrt(p3dmode_sin_t,'p3dmode_sin_t',lbound(p3dmode_sin_t),ubound(p3dmode_sin_t))
      if(allocated(u3dmode_cos_t))call netcdfrestart_wrt(u3dmode_cos_t,'u3dmode_cos_t',lbound(u3dmode_cos_t),ubound(u3dmode_cos_t))
      if(allocated(u3dmode_sin_t))call netcdfrestart_wrt(u3dmode_sin_t,'u3dmode_sin_t',lbound(u3dmode_sin_t),ubound(u3dmode_sin_t))
      if(allocated(v3dmode_cos_t))call netcdfrestart_wrt(v3dmode_cos_t,'v3dmode_cos_t',lbound(v3dmode_cos_t),ubound(v3dmode_cos_t))
      if(allocated(v3dmode_sin_t))call netcdfrestart_wrt(v3dmode_sin_t,'v3dmode_sin_t',lbound(v3dmode_sin_t),ubound(v3dmode_sin_t))
      if(allocated(analysetide3d_u))call netcdfrestart_wrt(analysetide3d_u,'analysetide3d_u',lbound(analysetide3d_u),ubound(analysetide3d_u))
      if(allocated(analysetide3d_v))call netcdfrestart_wrt(analysetide3d_v,'analysetide3d_v',lbound(analysetide3d_v),ubound(analysetide3d_v))
      if(allocated(analysetide3d_t))call netcdfrestart_wrt(analysetide3d_t,'analysetide3d_t',lbound(analysetide3d_t),ubound(analysetide3d_t))
      if(allocated(vel3dtidecosout_u))call netcdfrestart_wrt(vel3dtidecosout_u,'vel3dtidecosout_u',lbound(vel3dtidecosout_u),ubound(vel3dtidecosout_u))
      if(allocated(vel3dtidesinout_u))call netcdfrestart_wrt(vel3dtidesinout_u,'vel3dtidesinout_u',lbound(vel3dtidesinout_u),ubound(vel3dtidesinout_u))
      if(allocated(vel3dtidecosout_v))call netcdfrestart_wrt(vel3dtidecosout_v,'vel3dtidecosout_v',lbound(vel3dtidecosout_v),ubound(vel3dtidecosout_v))
      if(allocated(vel3dtidesinout_v))call netcdfrestart_wrt(vel3dtidesinout_v,'vel3dtidesinout_v',lbound(vel3dtidesinout_v),ubound(vel3dtidesinout_v))
      if(allocated(rhptidecosout_t))call netcdfrestart_wrt(rhptidecosout_t,'rhptidecosout_t',lbound(rhptidecosout_t),ubound(rhptidecosout_t))
      if(allocated(rhptidesinout_t))call netcdfrestart_wrt(rhptidesinout_t,'rhptidesinout_t',lbound(rhptidesinout_t),ubound(rhptidesinout_t))
      if(allocated(bio_relax_north))call netcdfrestart_wrt(bio_relax_north,'bio_relax_north',lbound(bio_relax_north),ubound(bio_relax_north))
      if(allocated(bio_relax_south))call netcdfrestart_wrt(bio_relax_south,'bio_relax_south',lbound(bio_relax_south),ubound(bio_relax_south))
      if(allocated(bio_relax_east))call netcdfrestart_wrt(bio_relax_east,'bio_relax_east',lbound(bio_relax_east),ubound(bio_relax_east))
      if(allocated(bio_relax_west))call netcdfrestart_wrt(bio_relax_west,'bio_relax_west',lbound(bio_relax_west),ubound(bio_relax_west))
      if(allocated(bio_relax_full))call netcdfrestart_wrt(bio_relax_full,'bio_relax_full',lbound(bio_relax_full),ubound(bio_relax_full))
      if(allocated(analysetide_u))call netcdfrestart_wrt(analysetide_u,'analysetide_u',lbound(analysetide_u),ubound(analysetide_u))
      if(allocated(analysetide_v))call netcdfrestart_wrt(analysetide_v,'analysetide_v',lbound(analysetide_v),ubound(analysetide_v))
      if(allocated(analysetide_w))call netcdfrestart_wrt(analysetide_w,'analysetide_w',lbound(analysetide_w),ubound(analysetide_w))
      if(allocated(analysetide2_w))call netcdfrestart_wrt(analysetide2_w,'analysetide2_w',lbound(analysetide2_w),ubound(analysetide2_w))
      if(allocated(cocosisi_f))call netcdfrestart_wrt(cocosisi_f,'cocosisi_f',lbound(cocosisi_f),ubound(cocosisi_f))
      if(allocated(gridrotcos_f))call netcdfrestart_wrt(gridrotcos_f,'gridrotcos_f',lbound(gridrotcos_f),ubound(gridrotcos_f))
      if(allocated(gridrotsin_f))call netcdfrestart_wrt(gridrotsin_f,'gridrotsin_f',lbound(gridrotsin_f),ubound(gridrotsin_f))
      if(allocated(tideanalysismatrix))call netcdfrestart_wrt(tideanalysismatrix,'tideanalysismatrix',lbound(tideanalysismatrix),ubound(tideanalysismatrix))
      if(allocated(iaveraged_in))call netcdfrestart_wrt(iaveraged_in,'iaveraged_in',lbound(iaveraged_in),ubound(iaveraged_in))
      if(allocated(iaveraged_out))call netcdfrestart_wrt(iaveraged_out,'iaveraged_out',lbound(iaveraged_out),ubound(iaveraged_out))
      if(allocated(iaverag1d_in))call netcdfrestart_wrt(iaverag1d_in,'iaverag1d_in',lbound(iaverag1d_in),ubound(iaverag1d_in))
      if(allocated(iaverag1d_out))call netcdfrestart_wrt(iaverag1d_out,'iaverag1d_out',lbound(iaverag1d_out),ubound(iaverag1d_out))
      if(allocated(light_kpar2_w))call netcdfrestart_wrt(light_kpar2_w,'light_kpar2_w',lbound(light_kpar2_w),ubound(light_kpar2_w))
      if(allocated(tidepotential_w))call netcdfrestart_wrt(tidepotential_w,'tidepotential_w',lbound(tidepotential_w),ubound(tidepotential_w))
      if(allocated(sshtidecos_w))call netcdfrestart_wrt(sshtidecos_w,'sshtidecos_w',lbound(sshtidecos_w),ubound(sshtidecos_w))
      if(allocated(sshtidesin_w))call netcdfrestart_wrt(sshtidesin_w,'sshtidesin_w',lbound(sshtidesin_w),ubound(sshtidesin_w))
      if(allocated(veltidecos_u))call netcdfrestart_wrt(veltidecos_u,'veltidecos_u',lbound(veltidecos_u),ubound(veltidecos_u))
      if(allocated(veltidesin_u))call netcdfrestart_wrt(veltidesin_u,'veltidesin_u',lbound(veltidesin_u),ubound(veltidesin_u))
      if(allocated(veltidecos_v))call netcdfrestart_wrt(veltidecos_v,'veltidecos_v',lbound(veltidecos_v),ubound(veltidecos_v))
      if(allocated(veltidesin_v))call netcdfrestart_wrt(veltidesin_v,'veltidesin_v',lbound(veltidesin_v),ubound(veltidesin_v))
      if(allocated(potidecos_w))call netcdfrestart_wrt(potidecos_w,'potidecos_w',lbound(potidecos_w),ubound(potidecos_w))
      if(allocated(potidesin_w))call netcdfrestart_wrt(potidesin_w,'potidesin_w',lbound(potidesin_w),ubound(potidesin_w))
      if(allocated(veltidecosout_u))call netcdfrestart_wrt(veltidecosout_u,'veltidecosout_u',lbound(veltidecosout_u),ubound(veltidecosout_u))
      if(allocated(veltidesinout_u))call netcdfrestart_wrt(veltidesinout_u,'veltidesinout_u',lbound(veltidesinout_u),ubound(veltidesinout_u))
      if(allocated(veltidecosout_v))call netcdfrestart_wrt(veltidecosout_v,'veltidecosout_v',lbound(veltidecosout_v),ubound(veltidecosout_v))
      if(allocated(veltidesinout_v))call netcdfrestart_wrt(veltidesinout_v,'veltidesinout_v',lbound(veltidesinout_v),ubound(veltidesinout_v))
      if(allocated(sshtidecosout_w))call netcdfrestart_wrt(sshtidecosout_w,'sshtidecosout_w',lbound(sshtidecosout_w),ubound(sshtidecosout_w))
      if(allocated(sshtidesinout_w))call netcdfrestart_wrt(sshtidesinout_w,'sshtidesinout_w',lbound(sshtidesinout_w),ubound(sshtidesinout_w))
      if(allocated(ssh3dtidecosout_w))call netcdfrestart_wrt(ssh3dtidecosout_w,'ssh3dtidecosout_w',lbound(ssh3dtidecosout_w),ubound(ssh3dtidecosout_w))
      if(allocated(ssh3dtidesinout_w))call netcdfrestart_wrt(ssh3dtidesinout_w,'ssh3dtidesinout_w',lbound(ssh3dtidesinout_w),ubound(ssh3dtidesinout_w))
      if(allocated(passetide))call netcdfrestart_wrt(passetide,'passetide',lbound(passetide),ubound(passetide))
      if(allocated(ftide))call netcdfrestart_wrt(ftide,'ftide',lbound(ftide),ubound(ftide))
      if(allocated(utide))call netcdfrestart_wrt(utide,'utide',lbound(utide),ubound(utide))
      if(allocated(northflux_sumt_v))call netcdfrestart_wrt(northflux_sumt_v,'northflux_sumt_v',lbound(northflux_sumt_v),ubound(northflux_sumt_v))
      if(allocated(southflux_sumt_v))call netcdfrestart_wrt(southflux_sumt_v,'southflux_sumt_v',lbound(southflux_sumt_v),ubound(southflux_sumt_v))
      if(allocated(drifter_l))call netcdfrestart_wrt(drifter_l,'drifter_l',lbound(drifter_l),ubound(drifter_l))
      if(allocated(velbot3d2d_u))call netcdfrestart_wrt(velbot3d2d_u,'velbot3d2d_u',lbound(velbot3d2d_u),ubound(velbot3d2d_u))
      if(allocated(velbot3d2d_v))call netcdfrestart_wrt(velbot3d2d_v,'velbot3d2d_v',lbound(velbot3d2d_v),ubound(velbot3d2d_v))
      if(allocated(velbot_u))call netcdfrestart_wrt(velbot_u,'velbot_u',lbound(velbot_u),ubound(velbot_u))
      if(allocated(velbot_v))call netcdfrestart_wrt(velbot_v,'velbot_v',lbound(velbot_v),ubound(velbot_v))
      if(allocated(kvectorpeak_i))call netcdfrestart_wrt(kvectorpeak_i,'kvectorpeak_i',lbound(kvectorpeak_i),ubound(kvectorpeak_i))
      if(allocated(kvectorpeak_j))call netcdfrestart_wrt(kvectorpeak_j,'kvectorpeak_j',lbound(kvectorpeak_j),ubound(kvectorpeak_j))
      if(allocated(drifter_send_nord))write(9)drifter_send_nord
      if(allocated(drifter_recv_nord))write(9)drifter_recv_nord
      if(allocated(drifter_send_sud))write(9)drifter_send_sud
      if(allocated(drifter_recv_sud))write(9)drifter_recv_sud
      if(allocated(drifter_send_est))write(9)drifter_send_est
      if(allocated(drifter_recv_est))write(9)drifter_recv_est
      if(allocated(drifter_send_ouest))write(9)drifter_send_ouest
      if(allocated(drifter_recv_ouest))write(9)drifter_recv_ouest
      if(allocated(drifter_send_nordest))write(9)drifter_send_nordest
      if(allocated(drifter_send_nordouest))write(9)drifter_send_nordouest
      if(allocated(drifter_send_sudest))write(9)drifter_send_sudest
      if(allocated(drifter_send_sudouest))write(9)drifter_send_sudouest
      if(allocated(drifter_recv_nordest))write(9)drifter_recv_nordest
      if(allocated(drifter_recv_nordouest))write(9)drifter_recv_nordouest
      if(allocated(drifter_recv_sudest))write(9)drifter_recv_sudest
      if(allocated(drifter_recv_sudouest))write(9)drifter_recv_sudouest
      if(allocated(drifter_send_canal))call netcdfrestart_wrt(drifter_send_canal,'drifter_send_canal',lbound(drifter_send_canal),ubound(drifter_send_canal))
      if(allocated(drifter_recv_canal))call netcdfrestart_wrt(drifter_recv_canal,'drifter_recv_canal',lbound(drifter_recv_canal),ubound(drifter_recv_canal))
      if(allocated(qwave_j_w))call netcdfrestart_wrt(qwave_j_w,'qwave_j_w',lbound(qwave_j_w),ubound(qwave_j_w))
      if(allocated(qwave_i_w))call netcdfrestart_wrt(qwave_i_w,'qwave_i_w',lbound(qwave_i_w),ubound(qwave_i_w))
      if(allocated(velwave_i_u))call netcdfrestart_wrt(velwave_i_u,'velwave_i_u',lbound(velwave_i_u),ubound(velwave_i_u))
      if(allocated(velwave_j_u))call netcdfrestart_wrt(velwave_j_u,'velwave_j_u',lbound(velwave_j_u),ubound(velwave_j_u))
      if(allocated(velwave_i_v))call netcdfrestart_wrt(velwave_i_v,'velwave_i_v',lbound(velwave_i_v),ubound(velwave_i_v))
      if(allocated(velwave_j_v))call netcdfrestart_wrt(velwave_j_v,'velwave_j_v',lbound(velwave_j_v),ubound(velwave_j_v))
      if(allocated(q_t))call netcdfrestart_wrt(q_t,'q_t',lbound(q_t),ubound(q_t))
      if(allocated(sshwave_j_w))call netcdfrestart_wrt(sshwave_j_w,'sshwave_j_w',lbound(sshwave_j_w),ubound(sshwave_j_w))
      if(allocated(sshwave_i_w))call netcdfrestart_wrt(sshwave_i_w,'sshwave_i_w',lbound(sshwave_i_w),ubound(sshwave_i_w))
      if(allocated(mc))call netcdfrestart_wrt(mc,'mc',lbound(mc),ubound(mc))
      if(allocated(anyv3dr4))call netcdfrestart_wrt(anyv3dr4,'anyv3dr4',lbound(anyv3dr4),ubound(anyv3dr4))
      if(allocated(mc_out))call netcdfrestart_wrt(mc_out,'mc_out',lbound(mc_out),ubound(mc_out))
      if(allocated(velbarwave_i_v))call netcdfrestart_wrt(velbarwave_i_v,'velbarwave_i_v',lbound(velbarwave_i_v),ubound(velbarwave_i_v))
      if(allocated(qavr_t))call netcdfrestart_wrt(qavr_t,'qavr_t',lbound(qavr_t),ubound(qavr_t))
      if(allocated(nhpgf_u))call netcdfrestart_wrt(nhpgf_u,'nhpgf_u',lbound(nhpgf_u),ubound(nhpgf_u))
      if(allocated(nhpgf_v))call netcdfrestart_wrt(nhpgf_v,'nhpgf_v',lbound(nhpgf_v),ubound(nhpgf_v))
      if(allocated(dsigr4_t))call netcdfrestart_wrt(dsigr4_t,'dsigr4_t',lbound(dsigr4_t),ubound(dsigr4_t))
      if(allocated(sig1dpom_t))call netcdfrestart_wrt(sig1dpom_t,'sig1dpom_t',lbound(sig1dpom_t),ubound(sig1dpom_t))
      if(allocated(obc_q_i))call netcdfrestart_wrt(obc_q_i,'obc_q_i',lbound(obc_q_i),ubound(obc_q_i))
      if(allocated(obc_ub_i))call netcdfrestart_wrt(obc_ub_i,'obc_ub_i',lbound(obc_ub_i),ubound(obc_ub_i))
      if(allocated(dsig_w))call netcdfrestart_wrt(dsig_w,'dsig_w',lbound(dsig_w),ubound(dsig_w))
      if(allocated(retarded_pss_w))call netcdfrestart_wrt(retarded_pss_w,'retarded_pss_w',lbound(retarded_pss_w),ubound(retarded_pss_w))
      if(allocated(velwf_u))call netcdfrestart_wrt(velwf_u,'velwf_u',lbound(velwf_u),ubound(velwf_u))
      if(allocated(velwf_v))call netcdfrestart_wrt(velwf_v,'velwf_v',lbound(velwf_v),ubound(velwf_v))
      if(allocated(kw_w))call netcdfrestart_wrt(kw_w,'kw_w',lbound(kw_w),ubound(kw_w))
      if(allocated(q2davr_w))call netcdfrestart_wrt(q2davr_w,'q2davr_w',lbound(q2davr_w),ubound(q2davr_w))
      if(allocated(breaker2d_t))call netcdfrestart_wrt(breaker2d_t,'breaker2d_t',lbound(breaker2d_t),ubound(breaker2d_t))
      if(allocated(invdx_t))call netcdfrestart_wrt(invdx_t,'invdx_t',lbound(invdx_t),ubound(invdx_t))
      if(allocated(invdx_f))call netcdfrestart_wrt(invdx_f,'invdx_f',lbound(invdx_f),ubound(invdx_f))
      if(allocated(invdx_u))call netcdfrestart_wrt(invdx_u,'invdx_u',lbound(invdx_u),ubound(invdx_u))
      if(allocated(invdy_u))call netcdfrestart_wrt(invdy_u,'invdy_u',lbound(invdy_u),ubound(invdy_u))
      if(allocated(invdy_t))call netcdfrestart_wrt(invdy_t,'invdy_t',lbound(invdy_t),ubound(invdy_t))
      if(allocated(invdy_f))call netcdfrestart_wrt(invdy_f,'invdy_f',lbound(invdy_f),ubound(invdy_f))
      if(allocated(invdy_v))call netcdfrestart_wrt(invdy_v,'invdy_v',lbound(invdy_v),ubound(invdy_v))
      if(allocated(invdx_v))call netcdfrestart_wrt(invdx_v,'invdx_v',lbound(invdx_v),ubound(invdx_v))
      if(allocated(invdxdy_t))call netcdfrestart_wrt(invdxdy_t,'invdxdy_t',lbound(invdxdy_t),ubound(invdxdy_t))
      if(allocated(invdxdy_u))call netcdfrestart_wrt(invdxdy_u,'invdxdy_u',lbound(invdxdy_u),ubound(invdxdy_u))
      if(allocated(invdxdy_v))call netcdfrestart_wrt(invdxdy_v,'invdxdy_v',lbound(invdxdy_v),ubound(invdxdy_v))
      if(allocated(nhpgf2d_u))call netcdfrestart_wrt(nhpgf2d_u,'nhpgf2d_u',lbound(nhpgf2d_u),ubound(nhpgf2d_u))
      if(allocated(nhpgf2d_v))call netcdfrestart_wrt(nhpgf2d_v,'nhpgf2d_v',lbound(nhpgf2d_v),ubound(nhpgf2d_v))
      if(allocated(sshmax_w))call netcdfrestart_wrt(sshmax_w,'sshmax_w',lbound(sshmax_w),ubound(sshmax_w))
      if(allocated(sshmin_w))call netcdfrestart_wrt(sshmin_w,'sshmin_w',lbound(sshmin_w),ubound(sshmin_w))
      if(allocated(sshmin_tmp_w))call netcdfrestart_wrt(sshmin_tmp_w,'sshmin_tmp_w',lbound(sshmin_tmp_w),ubound(sshmin_tmp_w))
      if(allocated(rwavebreak_t))call netcdfrestart_wrt(rwavebreak_t,'rwavebreak_t',lbound(rwavebreak_t),ubound(rwavebreak_t))
      if(allocated(breaker_t))call netcdfrestart_wrt(breaker_t,'breaker_t',lbound(breaker_t),ubound(breaker_t))
      if(allocated(hs_w))call netcdfrestart_wrt(hs_w,'hs_w',lbound(hs_w),ubound(hs_w))
      if(allocated(cx_upper_t))call netcdfrestart_wrt(cx_upper_t,'cx_upper_t',lbound(cx_upper_t),ubound(cx_upper_t))
      if(allocated(cy_upper_t))call netcdfrestart_wrt(cy_upper_t,'cy_upper_t',lbound(cy_upper_t),ubound(cy_upper_t))
      if(allocated(c_lower_t))call netcdfrestart_wrt(c_lower_t,'c_lower_t',lbound(c_lower_t),ubound(c_lower_t))
      if(allocated(slopemax_u))call netcdfrestart_wrt(slopemax_u,'slopemax_u',lbound(slopemax_u),ubound(slopemax_u))
      if(allocated(pgfwave_u))call netcdfrestart_wrt(pgfwave_u,'pgfwave_u',lbound(pgfwave_u),ubound(pgfwave_u))
      if(allocated(dt_over_retartedtime))write(9)dt_over_retartedtime
      if(allocated(sshrelax_j_w))call netcdfrestart_wrt(sshrelax_j_w,'sshrelax_j_w',lbound(sshrelax_j_w),ubound(sshrelax_j_w))
      if(allocated(sshrelax_i_w))call netcdfrestart_wrt(sshrelax_i_w,'sshrelax_i_w',lbound(sshrelax_i_w),ubound(sshrelax_i_w))
      write(9)asselin_nh
      write(9)wetdry_cstnh
      write(9)cfl_nh
      write(9)cwavepeak
      write(9)periodpeak
      write(9)freq2pipeak
      write(9)kvectorpeak
      write(9)kvector
      write(9)mukvector
      write(9)nhpgf_reduce
      write(9)acm_speed
      write(9)inv_brkh
      write(9)inv_brkslope2
      write(9)inv_brkvar
      write(9)brk_crit_slope
      write(9)brk_crit_var
      write(9)brk_crit_h
      write(9)brk_crit_r
      write(9)nh2d_graph_period
      write(9)nh2d_graph_spinup
      write(9)sum_avr_hs
      write(9)nh_frozensigma
      write(9)nh_wavebreakscheme
      write(9)flag_adve2d
      write(9)flag_adve3d
      write(9)flag_nh3d
      write(9)flag_nh3d_none
      write(9)flag_nh3d_nosplit_uv
      write(9)flag_nh3d_nosplit_tsuv
      write(9)flag_nh3d_timesplit_tsuv
      write(9)flag_timesplitting
      write(9)flag_timesplitting_adve_uv
      write(9)flagspo_i1
      write(9)flagspo_i2
      write(9)flagspo_j1
      write(9)flagspo_j2
      write(9)flag_ksloffline
      write(9)flag_groundwater
      write(9)flag_surfriver
      write(9)ofl_reversedtime
      write(9)obcnh_scheme
      write(9)obc_scheme
      write(9)drifter_random_w
      write(9)flag_buo_w
      write(9)flag_ogcmtidemixing
      write(9)flag_ogcm_instab
      write(9)flag_ofactor
      write(9)flag_omega_cumul
      write(9)flag_negdif_ver
      write(9)flag_negdif_hor
      write(9)flag_z0_macro
      write(9)oasis_symsym_onoff
      write(9)oasis_symsym_retrots
      if(allocated(variance_timeaveraged))call netcdfrestart_wrt(variance_timeaveraged,'variance_timeaveraged',lbound(variance_timeaveraged),ubound(variance_timeaveraged))
      if(allocated(ssh_timeaveraged))call netcdfrestart_wrt(ssh_timeaveraged,'ssh_timeaveraged',lbound(ssh_timeaveraged),ubound(ssh_timeaveraged))
      if(allocated(flx2d_timeaveraged))call netcdfrestart_wrt(flx2d_timeaveraged,'flx2d_timeaveraged',lbound(flx2d_timeaveraged),ubound(flx2d_timeaveraged))
      if(allocated(fly2d_timeaveraged))call netcdfrestart_wrt(fly2d_timeaveraged,'fly2d_timeaveraged',lbound(fly2d_timeaveraged),ubound(fly2d_timeaveraged))
      if(allocated(flx3d_timeaveraged))call netcdfrestart_wrt(flx3d_timeaveraged,'flx3d_timeaveraged',lbound(flx3d_timeaveraged),ubound(flx3d_timeaveraged))
      if(allocated(fly3d_timeaveraged))call netcdfrestart_wrt(fly3d_timeaveraged,'fly3d_timeaveraged',lbound(fly3d_timeaveraged),ubound(fly3d_timeaveraged))
      if(allocated(u_euler_timeaveraged))call netcdfrestart_wrt(u_euler_timeaveraged,'u_euler_timeaveraged',lbound(u_euler_timeaveraged),ubound(u_euler_timeaveraged))
      if(allocated(v_euler_timeaveraged))call netcdfrestart_wrt(v_euler_timeaveraged,'v_euler_timeaveraged',lbound(v_euler_timeaveraged),ubound(v_euler_timeaveraged))
      if(allocated(tkeb_w))call netcdfrestart_wrt(tkeb_w,'tkeb_w',lbound(tkeb_w),ubound(tkeb_w))
      if(allocated(tken_w))call netcdfrestart_wrt(tken_w,'tken_w',lbound(tken_w),ubound(tken_w))
      if(allocated(tkea_w))call netcdfrestart_wrt(tkea_w,'tkea_w',lbound(tkea_w),ubound(tkea_w))
      if(allocated(tkle_w))call netcdfrestart_wrt(tkle_w,'tkle_w',lbound(tkle_w),ubound(tkle_w))
      if(allocated(tkll_w))call netcdfrestart_wrt(tkll_w,'tkll_w',lbound(tkll_w),ubound(tkll_w))
      if(allocated(epsb_w))call netcdfrestart_wrt(epsb_w,'epsb_w',lbound(epsb_w),ubound(epsb_w))
      if(allocated(epsn_w))call netcdfrestart_wrt(epsn_w,'epsn_w',lbound(epsn_w),ubound(epsn_w))
      if(allocated(epsa_w))call netcdfrestart_wrt(epsa_w,'epsa_w',lbound(epsa_w),ubound(epsa_w))
      if(allocated(gradssh_u))call netcdfrestart_wrt(gradssh_u,'gradssh_u',lbound(gradssh_u),ubound(gradssh_u))
      if(allocated(gradssh_v))call netcdfrestart_wrt(gradssh_v,'gradssh_v',lbound(gradssh_v),ubound(gradssh_v))
      if(allocated(pgf_u))call netcdfrestart_wrt(pgf_u,'pgf_u',lbound(pgf_u),ubound(pgf_u))
      if(allocated(pgf_v))call netcdfrestart_wrt(pgf_v,'pgf_v',lbound(pgf_v),ubound(pgf_v))
      if(allocated(dz_u))call netcdfrestart_wrt(dz_u,'dz_u',lbound(dz_u),ubound(dz_u))
      if(allocated(dz_v))call netcdfrestart_wrt(dz_v,'dz_v',lbound(dz_v),ubound(dz_v))
      if(allocated(dz_t))call netcdfrestart_wrt(dz_t,'dz_t',lbound(dz_t),ubound(dz_t))
      if(allocated(depth_t))call netcdfrestart_wrt(depth_t,'depth_t',lbound(depth_t),ubound(depth_t))
      if(allocated(depth_u))call netcdfrestart_wrt(depth_u,'depth_u',lbound(depth_u),ubound(depth_u))
      if(allocated(depth_v))call netcdfrestart_wrt(depth_v,'depth_v',lbound(depth_v),ubound(depth_v))
      if(allocated(depth_w))call netcdfrestart_wrt(depth_w,'depth_w',lbound(depth_w),ubound(depth_w))
      if(allocated(depth_f))call netcdfrestart_wrt(depth_f,'depth_f',lbound(depth_f),ubound(depth_f))
      if(allocated(stokesforces_u))call netcdfrestart_wrt(stokesforces_u,'stokesforces_u',lbound(stokesforces_u),ubound(stokesforces_u))
      if(allocated(stokesforces_v))call netcdfrestart_wrt(stokesforces_v,'stokesforces_v',lbound(stokesforces_v),ubound(stokesforces_v))
      if(allocated(r_mangrove))call netcdfrestart_wrt(r_mangrove,'r_mangrove',lbound(r_mangrove),ubound(r_mangrove))
      if(allocated(sigma_w))call netcdfrestart_wrt(sigma_w,'sigma_w',lbound(sigma_w),ubound(sigma_w))
      if(allocated(dsig_f))call netcdfrestart_wrt(dsig_f,'dsig_f',lbound(dsig_f),ubound(dsig_f))
      if(allocated(rhp_t))call netcdfrestart_wrt(rhp_t,'rhp_t',lbound(rhp_t),ubound(rhp_t))
      if(allocated(rhcref_t))call netcdfrestart_wrt(rhcref_t,'rhcref_t',lbound(rhcref_t),ubound(rhcref_t))
      if(allocated(vel_u))call netcdfrestart_wrt(vel_u,'vel_u',lbound(vel_u),ubound(vel_u))
      if(allocated(vel_v))call netcdfrestart_wrt(vel_v,'vel_v',lbound(vel_v),ubound(vel_v))
      if(allocated(presgrad_u))call netcdfrestart_wrt(presgrad_u,'presgrad_u',lbound(presgrad_u),ubound(presgrad_u))
      if(allocated(presgrad_v))call netcdfrestart_wrt(presgrad_v,'presgrad_v',lbound(presgrad_v),ubound(presgrad_v))
      if(allocated(omega_w))call netcdfrestart_wrt(omega_w,'omega_w',lbound(omega_w),ubound(omega_w))
      if(allocated(veldydz_u))call netcdfrestart_wrt(veldydz_u,'veldydz_u',lbound(veldydz_u),ubound(veldydz_u))
      if(allocated(veldxdz_v))call netcdfrestart_wrt(veldxdz_v,'veldxdz_v',lbound(veldxdz_v),ubound(veldxdz_v))
      if(allocated(tem_t))call netcdfrestart_wrt(tem_t,'tem_t',lbound(tem_t),ubound(tem_t))
      if(allocated(sal_t))call netcdfrestart_wrt(sal_t,'sal_t',lbound(sal_t),ubound(sal_t))
      if(allocated(tridia_in))call netcdfrestart_wrt(tridia_in,'tridia_in',lbound(tridia_in),ubound(tridia_in))
      if(allocated(hr_z2lr_w))call netcdfrestart_wrt(hr_z2lr_w,'hr_z2lr_w',lbound(hr_z2lr_w),ubound(hr_z2lr_w))
      if(allocated(anyv3d))call netcdfrestart_wrt(anyv3d,'anyv3d',lbound(anyv3d),ubound(anyv3d))
      if(allocated(sponge_t))call netcdfrestart_wrt(sponge_t,'sponge_t',lbound(sponge_t),ubound(sponge_t))
      if(allocated(sponge_u))call netcdfrestart_wrt(sponge_u,'sponge_u',lbound(sponge_u),ubound(sponge_u))
      if(allocated(sponge_v))call netcdfrestart_wrt(sponge_v,'sponge_v',lbound(sponge_v),ubound(sponge_v))
      if(allocated(uwind_t))call netcdfrestart_wrt(uwind_t,'uwind_t',lbound(uwind_t),ubound(uwind_t))
      if(allocated(vwind_t))call netcdfrestart_wrt(vwind_t,'vwind_t',lbound(vwind_t),ubound(vwind_t))
      if(allocated(uwind100_t))call netcdfrestart_wrt(uwind100_t,'uwind100_t',lbound(uwind100_t),ubound(uwind100_t))
      if(allocated(vwind100_t))call netcdfrestart_wrt(vwind100_t,'vwind100_t',lbound(vwind100_t),ubound(vwind100_t))
      if(allocated(sshf_w))call netcdfrestart_wrt(sshf_w,'sshf_w',lbound(sshf_w),ubound(sshf_w))
      if(allocated(slhf_w))call netcdfrestart_wrt(slhf_w,'slhf_w',lbound(slhf_w),ubound(slhf_w))
      if(allocated(ssr_w))call netcdfrestart_wrt(ssr_w,'ssr_w',lbound(ssr_w),ubound(ssr_w))
      if(allocated(ssr24prv_w))call netcdfrestart_wrt(ssr24prv_w,'ssr24prv_w',lbound(ssr24prv_w),ubound(ssr24prv_w))
      if(allocated(snsf_w))call netcdfrestart_wrt(snsf_w,'snsf_w',lbound(snsf_w),ubound(snsf_w))
      if(allocated(precipi_w))call netcdfrestart_wrt(precipi_w,'precipi_w',lbound(precipi_w),ubound(precipi_w))
      if(allocated(taux_w))call netcdfrestart_wrt(taux_w,'taux_w',lbound(taux_w),ubound(taux_w))
      if(allocated(tauy_w))call netcdfrestart_wrt(tauy_w,'tauy_w',lbound(tauy_w),ubound(tauy_w))
      if(allocated(sigma_f))call netcdfrestart_wrt(sigma_f,'sigma_f',lbound(sigma_f),ubound(sigma_f))
      if(allocated(km_w))call netcdfrestart_wrt(km_w,'km_w',lbound(km_w),ubound(km_w))
      if(allocated(kh_w))call netcdfrestart_wrt(kh_w,'kh_w',lbound(kh_w),ubound(kh_w))
      if(allocated(rho_t))call netcdfrestart_wrt(rho_t,'rho_t',lbound(rho_t),ubound(rho_t))
      if(allocated(omega_evaprec_w))call netcdfrestart_wrt(omega_evaprec_w,'omega_evaprec_w',lbound(omega_evaprec_w),ubound(omega_evaprec_w))
      if(allocated(tridia_out))call netcdfrestart_wrt(tridia_out,'tridia_out',lbound(tridia_out),ubound(tridia_out))
      if(allocated(fluxbar_sumt_u))call netcdfrestart_wrt(fluxbar_sumt_u,'fluxbar_sumt_u',lbound(fluxbar_sumt_u),ubound(fluxbar_sumt_u))
      if(allocated(fluxbar_sumt_v))call netcdfrestart_wrt(fluxbar_sumt_v,'fluxbar_sumt_v',lbound(fluxbar_sumt_v),ubound(fluxbar_sumt_v))
      if(allocated(velbar_u))call netcdfrestart_wrt(velbar_u,'velbar_u',lbound(velbar_u),ubound(velbar_u))
      if(allocated(uflux2d_f))call netcdfrestart_wrt(uflux2d_f,'uflux2d_f',lbound(uflux2d_f),ubound(uflux2d_f))
      if(allocated(vflux2d_f))call netcdfrestart_wrt(vflux2d_f,'vflux2d_f',lbound(vflux2d_f),ubound(vflux2d_f))
      if(allocated(vlxbar_f))call netcdfrestart_wrt(vlxbar_f,'vlxbar_f',lbound(vlxbar_f),ubound(vlxbar_f))
      if(allocated(velbar_v))call netcdfrestart_wrt(velbar_v,'velbar_v',lbound(velbar_v),ubound(velbar_v))
      if(allocated(vlybar_f))call netcdfrestart_wrt(vlybar_f,'vlybar_f',lbound(vlybar_f),ubound(vlybar_f))
      if(allocated(flux2d_u))call netcdfrestart_wrt(flux2d_u,'flux2d_u',lbound(flux2d_u),ubound(flux2d_u))
      if(allocated(flux2d_v))call netcdfrestart_wrt(flux2d_v,'flux2d_v',lbound(flux2d_v),ubound(flux2d_v))
      if(allocated(ssh_w))call netcdfrestart_wrt(ssh_w,'ssh_w',lbound(ssh_w),ubound(ssh_w))
      if(allocated(hssh_f))call netcdfrestart_wrt(hssh_f,'hssh_f',lbound(hssh_f),ubound(hssh_f))
      if(allocated(hz_u))call netcdfrestart_wrt(hz_u,'hz_u',lbound(hz_u),ubound(hz_u))
      if(allocated(hz_v))call netcdfrestart_wrt(hz_v,'hz_v',lbound(hz_v),ubound(hz_v))
      if(allocated(hz_w))call netcdfrestart_wrt(hz_w,'hz_w',lbound(hz_w),ubound(hz_w))
      if(allocated(heatrelax_w))call netcdfrestart_wrt(heatrelax_w,'heatrelax_w',lbound(heatrelax_w),ubound(heatrelax_w))
      if(allocated(xy_u))call netcdfrestart_wrt(xy_u,'xy_u',lbound(xy_u),ubound(xy_u))
      if(allocated(xy_v))call netcdfrestart_wrt(xy_v,'xy_v',lbound(xy_v),ubound(xy_v))
      if(allocated(xy_t))call netcdfrestart_wrt(xy_t,'xy_t',lbound(xy_t),ubound(xy_t))
      if(allocated(xy_f))call netcdfrestart_wrt(xy_f,'xy_f',lbound(xy_f),ubound(xy_f))
      if(allocated(pss_w))call netcdfrestart_wrt(pss_w,'pss_w',lbound(pss_w),ubound(pss_w))
      if(allocated(wstress_u))call netcdfrestart_wrt(wstress_u,'wstress_u',lbound(wstress_u),ubound(wstress_u))
      if(allocated(wstress_v))call netcdfrestart_wrt(wstress_v,'wstress_v',lbound(wstress_v),ubound(wstress_v))
      if(allocated(ssh_int_w))call netcdfrestart_wrt(ssh_int_w,'ssh_int_w',lbound(ssh_int_w),ubound(ssh_int_w))
      if(allocated(ustokesvortex_t))call netcdfrestart_wrt(ustokesvortex_t,'ustokesvortex_t',lbound(ustokesvortex_t),ubound(ustokesvortex_t))
      if(allocated(ustokesvortex_f))call netcdfrestart_wrt(ustokesvortex_f,'ustokesvortex_f',lbound(ustokesvortex_f),ubound(ustokesvortex_f))
      if(allocated(vstokesvortex_t))call netcdfrestart_wrt(vstokesvortex_t,'vstokesvortex_t',lbound(vstokesvortex_t),ubound(vstokesvortex_t))
      if(allocated(vstokesvortex_f))call netcdfrestart_wrt(vstokesvortex_f,'vstokesvortex_f',lbound(vstokesvortex_f),ubound(vstokesvortex_f))
      if(allocated(velavr_u))call netcdfrestart_wrt(velavr_u,'velavr_u',lbound(velavr_u),ubound(velavr_u))
      if(allocated(velavr_v))call netcdfrestart_wrt(velavr_v,'velavr_v',lbound(velavr_v),ubound(velavr_v))
      if(allocated(timefilter_u))call netcdfrestart_wrt(timefilter_u,'timefilter_u',lbound(timefilter_u),ubound(timefilter_u))
      if(allocated(timefilter_v))call netcdfrestart_wrt(timefilter_v,'timefilter_v',lbound(timefilter_v),ubound(timefilter_v))
      if(allocated(teta0_t))call netcdfrestart_wrt(teta0_t,'teta0_t',lbound(teta0_t),ubound(teta0_t))
      if(allocated(teta2_t))call netcdfrestart_wrt(teta2_t,'teta2_t',lbound(teta2_t),ubound(teta2_t))
      if(allocated(q2_t))call netcdfrestart_wrt(q2_t,'q2_t',lbound(q2_t),ubound(q2_t))
      if(allocated(teta2delta_t))call netcdfrestart_wrt(teta2delta_t,'teta2delta_t',lbound(teta2delta_t),ubound(teta2delta_t))
      if(allocated(q2delta_t))call netcdfrestart_wrt(q2delta_t,'q2delta_t',lbound(q2delta_t),ubound(q2delta_t))
      if(allocated(uwinddelta_t))call netcdfrestart_wrt(uwinddelta_t,'uwinddelta_t',lbound(uwinddelta_t),ubound(uwinddelta_t))
      if(allocated(vwinddelta_t))call netcdfrestart_wrt(vwinddelta_t,'vwinddelta_t',lbound(vwinddelta_t),ubound(vwinddelta_t))
      if(allocated(tetastar_t))call netcdfrestart_wrt(tetastar_t,'tetastar_t',lbound(tetastar_t),ubound(tetastar_t))
      if(allocated(ustar_t))call netcdfrestart_wrt(ustar_t,'ustar_t',lbound(ustar_t),ubound(ustar_t))
      if(allocated(qstar_t))call netcdfrestart_wrt(qstar_t,'qstar_t',lbound(qstar_t),ubound(qstar_t))
      if(allocated(frozenterm2d_u))call netcdfrestart_wrt(frozenterm2d_u,'frozenterm2d_u',lbound(frozenterm2d_u),ubound(frozenterm2d_u))
      if(allocated(frozenterm2d_v))call netcdfrestart_wrt(frozenterm2d_v,'frozenterm2d_v',lbound(frozenterm2d_v),ubound(frozenterm2d_v))
      if(allocated(frozenterm3d_u))call netcdfrestart_wrt(frozenterm3d_u,'frozenterm3d_u',lbound(frozenterm3d_u),ubound(frozenterm3d_u))
      if(allocated(frozenterm3d_v))call netcdfrestart_wrt(frozenterm3d_v,'frozenterm3d_v',lbound(frozenterm3d_v),ubound(frozenterm3d_v))
      if(allocated(fluxbar_u))call netcdfrestart_wrt(fluxbar_u,'fluxbar_u',lbound(fluxbar_u),ubound(fluxbar_u))
      if(allocated(fluxbar_v))call netcdfrestart_wrt(fluxbar_v,'fluxbar_v',lbound(fluxbar_v),ubound(fluxbar_v))
      if(allocated(lon_t))call netcdfrestart_wrt(lon_t,'lon_t',lbound(lon_t),ubound(lon_t))
      if(allocated(lat_t))call netcdfrestart_wrt(lat_t,'lat_t',lbound(lat_t),ubound(lat_t))
      if(allocated(lon_u))call netcdfrestart_wrt(lon_u,'lon_u',lbound(lon_u),ubound(lon_u))
      if(allocated(lat_u))call netcdfrestart_wrt(lat_u,'lat_u',lbound(lat_u),ubound(lat_u))
      if(allocated(lon_v))call netcdfrestart_wrt(lon_v,'lon_v',lbound(lon_v),ubound(lon_v))
      if(allocated(lat_v))call netcdfrestart_wrt(lat_v,'lat_v',lbound(lat_v),ubound(lat_v))
      if(allocated(lon_f))call netcdfrestart_wrt(lon_f,'lon_f',lbound(lon_f),ubound(lon_f))
      if(allocated(lat_f))call netcdfrestart_wrt(lat_f,'lat_f',lbound(lat_f),ubound(lat_f))
      if(allocated(globlon_t))call netcdfrestart_wrt(globlon_t,'globlon_t',lbound(globlon_t),ubound(globlon_t))
      if(allocated(globlat_t))call netcdfrestart_wrt(globlat_t,'globlat_t',lbound(globlat_t),ubound(globlat_t))
      if(allocated(globlon_u))call netcdfrestart_wrt(globlon_u,'globlon_u',lbound(globlon_u),ubound(globlon_u))
      if(allocated(globlat_u))call netcdfrestart_wrt(globlat_u,'globlat_u',lbound(globlat_u),ubound(globlat_u))
      if(allocated(globlon_v))call netcdfrestart_wrt(globlon_v,'globlon_v',lbound(globlon_v),ubound(globlon_v))
      if(allocated(globlat_v))call netcdfrestart_wrt(globlat_v,'globlat_v',lbound(globlat_v),ubound(globlat_v))
      if(allocated(fluxbarsum_ofl_u))call netcdfrestart_wrt(fluxbarsum_ofl_u,'fluxbarsum_ofl_u',lbound(fluxbarsum_ofl_u),ubound(fluxbarsum_ofl_u))
      if(allocated(fluxbarsum_ofl_v))call netcdfrestart_wrt(fluxbarsum_ofl_v,'fluxbarsum_ofl_v',lbound(fluxbarsum_ofl_v),ubound(fluxbarsum_ofl_v))
      if(allocated(dxdy_u))call netcdfrestart_wrt(dxdy_u,'dxdy_u',lbound(dxdy_u),ubound(dxdy_u))
      if(allocated(dxdy_v))call netcdfrestart_wrt(dxdy_v,'dxdy_v',lbound(dxdy_v),ubound(dxdy_v))
      if(allocated(dxdy_t))call netcdfrestart_wrt(dxdy_t,'dxdy_t',lbound(dxdy_t),ubound(dxdy_t))
      if(allocated(dxdy_e_t))call netcdfrestart_wrt(dxdy_e_t,'dxdy_e_t',lbound(dxdy_e_t),ubound(dxdy_e_t))
      if(allocated(dx_e_f))call netcdfrestart_wrt(dx_e_f,'dx_e_f',lbound(dx_e_f),ubound(dx_e_f))
      if(allocated(dy_e_f))call netcdfrestart_wrt(dy_e_f,'dy_e_f',lbound(dy_e_f),ubound(dy_e_f))
      if(allocated(dx_u))call netcdfrestart_wrt(dx_u,'dx_u',lbound(dx_u),ubound(dx_u))
      if(allocated(dx_v))call netcdfrestart_wrt(dx_v,'dx_v',lbound(dx_v),ubound(dx_v))
      if(allocated(dx_t))call netcdfrestart_wrt(dx_t,'dx_t',lbound(dx_t),ubound(dx_t))
      if(allocated(dx_f))call netcdfrestart_wrt(dx_f,'dx_f',lbound(dx_f),ubound(dx_f))
      if(allocated(dy_u))call netcdfrestart_wrt(dy_u,'dy_u',lbound(dy_u),ubound(dy_u))
      if(allocated(dy_v))call netcdfrestart_wrt(dy_v,'dy_v',lbound(dy_v),ubound(dy_v))
      if(allocated(dy_t))call netcdfrestart_wrt(dy_t,'dy_t',lbound(dy_t),ubound(dy_t))
      if(allocated(dy_f))call netcdfrestart_wrt(dy_f,'dy_f',lbound(dy_f),ubound(dy_f))
      if(allocated(h_w))call netcdfrestart_wrt(h_w,'h_w',lbound(h_w),ubound(h_w))
      if(allocated(hcopy_w))call netcdfrestart_wrt(hcopy_w,'hcopy_w',lbound(hcopy_w),ubound(hcopy_w))
      if(allocated(h0_w))call netcdfrestart_wrt(h0_w,'h0_w',lbound(h0_w),ubound(h0_w))
      if(allocated(h_u))call netcdfrestart_wrt(h_u,'h_u',lbound(h_u),ubound(h_u))
      if(allocated(h_v))call netcdfrestart_wrt(h_v,'h_v',lbound(h_v),ubound(h_v))
      if(allocated(h_f))call netcdfrestart_wrt(h_f,'h_f',lbound(h_f),ubound(h_f))
      if(allocated(coriolis_t))call netcdfrestart_wrt(coriolis_t,'coriolis_t',lbound(coriolis_t),ubound(coriolis_t))
      if(allocated(coriolis_f))call netcdfrestart_wrt(coriolis_f,'coriolis_f',lbound(coriolis_f),ubound(coriolis_f))
      if(allocated(rhpzavr_w))call netcdfrestart_wrt(rhpzavr_w,'rhpzavr_w',lbound(rhpzavr_w),ubound(rhpzavr_w))
      if(allocated(q10_t))call netcdfrestart_wrt(q10_t,'q10_t',lbound(q10_t),ubound(q10_t))
      if(allocated(teta10_t))call netcdfrestart_wrt(teta10_t,'teta10_t',lbound(teta10_t),ubound(teta10_t))
      if(allocated(fric_u))call netcdfrestart_wrt(fric_u,'fric_u',lbound(fric_u),ubound(fric_u))
      if(allocated(fric_v))call netcdfrestart_wrt(fric_v,'fric_v',lbound(fric_v),ubound(fric_v))
      if(allocated(fric_t))call netcdfrestart_wrt(fric_t,'fric_t',lbound(fric_t),ubound(fric_t))
      if(allocated(cdb_t))call netcdfrestart_wrt(cdb_t,'cdb_t',lbound(cdb_t),ubound(cdb_t))
      if(allocated(cdb_f))call netcdfrestart_wrt(cdb_f,'cdb_f',lbound(cdb_f),ubound(cdb_f))
      if(allocated(xflux_t))call netcdfrestart_wrt(xflux_t,'xflux_t',lbound(xflux_t),ubound(xflux_t))
      if(allocated(xflux_f))call netcdfrestart_wrt(xflux_f,'xflux_f',lbound(xflux_f),ubound(xflux_f))
      if(allocated(yflux_t))call netcdfrestart_wrt(yflux_t,'yflux_t',lbound(yflux_t),ubound(yflux_t))
      if(allocated(yflux_f))call netcdfrestart_wrt(yflux_f,'yflux_f',lbound(yflux_f),ubound(yflux_f))
      if(allocated(pres3d2d_u))call netcdfrestart_wrt(pres3d2d_u,'pres3d2d_u',lbound(pres3d2d_u),ubound(pres3d2d_u))
      if(allocated(pres3d2d_v))call netcdfrestart_wrt(pres3d2d_v,'pres3d2d_v',lbound(pres3d2d_v),ubound(pres3d2d_v))
      if(allocated(adve3d2d_u))call netcdfrestart_wrt(adve3d2d_u,'adve3d2d_u',lbound(adve3d2d_u),ubound(adve3d2d_u))
      if(allocated(adve3d2d_v))call netcdfrestart_wrt(adve3d2d_v,'adve3d2d_v',lbound(adve3d2d_v),ubound(adve3d2d_v))
      if(allocated(rmangrovebar))call netcdfrestart_wrt(rmangrovebar,'rmangrovebar',lbound(rmangrovebar),ubound(rmangrovebar))
      if(allocated(mang3dto2d_u))call netcdfrestart_wrt(mang3dto2d_u,'mang3dto2d_u',lbound(mang3dto2d_u),ubound(mang3dto2d_u))
      if(allocated(mang3dto2d_v))call netcdfrestart_wrt(mang3dto2d_v,'mang3dto2d_v',lbound(mang3dto2d_v),ubound(mang3dto2d_v))
      if(allocated(restoring3d2d_u))call netcdfrestart_wrt(restoring3d2d_u,'restoring3d2d_u',lbound(restoring3d2d_u),ubound(restoring3d2d_u))
      if(allocated(restoring3d2d_v))call netcdfrestart_wrt(restoring3d2d_v,'restoring3d2d_v',lbound(restoring3d2d_v),ubound(restoring3d2d_v))
      if(allocated(stokesforces3d2d_u))call netcdfrestart_wrt(stokesforces3d2d_u,'stokesforces3d2d_u',lbound(stokesforces3d2d_u),ubound(stokesforces3d2d_u))
      if(allocated(stokesforces3d2d_v))call netcdfrestart_wrt(stokesforces3d2d_v,'stokesforces3d2d_v',lbound(stokesforces3d2d_v),ubound(stokesforces3d2d_v))
      if(allocated(wstress_w))call netcdfrestart_wrt(wstress_w,'wstress_w',lbound(wstress_w),ubound(wstress_w))
      if(allocated(z0_w))call netcdfrestart_wrt(z0_w,'z0_w',lbound(z0_w),ubound(z0_w))
      if(allocated(albedo_w))call netcdfrestart_wrt(albedo_w,'albedo_w',lbound(albedo_w),ubound(albedo_w))
      if(allocated(gridrotcos_t))call netcdfrestart_wrt(gridrotcos_t,'gridrotcos_t',lbound(gridrotcos_t),ubound(gridrotcos_t))
      if(allocated(gridrotsin_t))call netcdfrestart_wrt(gridrotsin_t,'gridrotsin_t',lbound(gridrotsin_t),ubound(gridrotsin_t))
      if(allocated(grid_angle_t))call netcdfrestart_wrt(grid_angle_t,'grid_angle_t',lbound(grid_angle_t),ubound(grid_angle_t))
      if(allocated(sshstokes_w))call netcdfrestart_wrt(sshstokes_w,'sshstokes_w',lbound(sshstokes_w),ubound(sshstokes_w))
      if(allocated(cwi_int_u))call netcdfrestart_wrt(cwi_int_u,'cwi_int_u',lbound(cwi_int_u),ubound(cwi_int_u))
      if(allocated(cwi_int_v))call netcdfrestart_wrt(cwi_int_v,'cwi_int_v',lbound(cwi_int_v),ubound(cwi_int_v))
      if(allocated(cwj_int_u))call netcdfrestart_wrt(cwj_int_u,'cwj_int_u',lbound(cwj_int_u),ubound(cwj_int_u))
      if(allocated(cwj_int_v))call netcdfrestart_wrt(cwj_int_v,'cwj_int_v',lbound(cwj_int_v),ubound(cwj_int_v))
      if(allocated(sshrefobc_i))call netcdfrestart_wrt(sshrefobc_i,'sshrefobc_i',lbound(sshrefobc_i),ubound(sshrefobc_i))
      if(allocated(vbrrefobc_i))call netcdfrestart_wrt(vbrrefobc_i,'vbrrefobc_i',lbound(vbrrefobc_i),ubound(vbrrefobc_i))
      if(allocated(sshrefobc_j))call netcdfrestart_wrt(sshrefobc_j,'sshrefobc_j',lbound(sshrefobc_j),ubound(sshrefobc_j))
      if(allocated(vbrrefobc_j))call netcdfrestart_wrt(vbrrefobc_j,'vbrrefobc_j',lbound(vbrrefobc_j),ubound(vbrrefobc_j))
      if(allocated(anyv1d))call netcdfrestart_wrt(anyv1d,'anyv1d',lbound(anyv1d),ubound(anyv1d))
      if(allocated(gridcard))call netcdfrestart_wrt(gridcard,'gridcard',lbound(gridcard),ubound(gridcard))
      if(allocated(novector))call netcdfrestart_wrt(novector,'novector',lbound(novector),ubound(novector))
      if(allocated(proficard))call netcdfrestart_wrt(proficard,'proficard',lbound(proficard),ubound(proficard))
      if(allocated(airseainfo))call netcdfrestart_wrt(airseainfo,'airseainfo',lbound(airseainfo),ubound(airseainfo))
      if(allocated(airseadt))call netcdfrestart_wrt(airseadt,'airseadt',lbound(airseadt),ubound(airseadt))
      if(allocated(riverdt))call netcdfrestart_wrt(riverdt,'riverdt',lbound(riverdt),ubound(riverdt))
      if(allocated(river_t))call netcdfrestart_wrt(river_t,'river_t',lbound(river_t),ubound(river_t))
      if(allocated(riverflux))call netcdfrestart_wrt(riverflux,'riverflux',lbound(riverflux),ubound(riverflux))
      if(allocated(mask_t))call netcdfrestart_wrt(mask_t,'mask_t',lbound(mask_t),ubound(mask_t))
      if(allocated(mask_f))call netcdfrestart_wrt(mask_f,'mask_f',lbound(mask_f),ubound(mask_f))
      if(allocated(mask_u))call netcdfrestart_wrt(mask_u,'mask_u',lbound(mask_u),ubound(mask_u))
      if(allocated(mask_v))call netcdfrestart_wrt(mask_v,'mask_v',lbound(mask_v),ubound(mask_v))
      if(allocated(mask_vqs_tke_w))call netcdfrestart_wrt(mask_vqs_tke_w,'mask_vqs_tke_w',lbound(mask_vqs_tke_w),ubound(mask_vqs_tke_w))
      if(allocated(canaldir))call netcdfrestart_wrt(canaldir,'canaldir',lbound(canaldir),ubound(canaldir))
      if(allocated(wetmask_wi_t))call netcdfrestart_wrt(wetmask_wi_t,'wetmask_wi_t',lbound(wetmask_wi_t),ubound(wetmask_wi_t))
      if(allocated(lonlat2ij_t))call netcdfrestart_wrt(lonlat2ij_t,'lonlat2ij_t',lbound(lonlat2ij_t),ubound(lonlat2ij_t))
      if(allocated(canalmpioverlap))call netcdfrestart_wrt(canalmpioverlap,'canalmpioverlap',lbound(canalmpioverlap),ubound(canalmpioverlap))
      if(allocated(sodate))call netcdfrestart_wrt(sodate,'sodate',lbound(sodate),ubound(sodate))
      if(allocated(canalrank))call netcdfrestart_wrt(canalrank,'canalrank',lbound(canalrank),ubound(canalrank))
      if(allocated(canalrankbis))call netcdfrestart_wrt(canalrankbis,'canalrankbis',lbound(canalrankbis),ubound(canalrankbis))
      if(allocated(i_canalcoord))call netcdfrestart_wrt(i_canalcoord,'i_canalcoord',lbound(i_canalcoord),ubound(i_canalcoord))
      if(allocated(j_canalcoord))call netcdfrestart_wrt(j_canalcoord,'j_canalcoord',lbound(j_canalcoord),ubound(j_canalcoord))
      if(allocated(kmin_u))call netcdfrestart_wrt(kmin_u,'kmin_u',lbound(kmin_u),ubound(kmin_u))
      if(allocated(kmin_v))call netcdfrestart_wrt(kmin_v,'kmin_v',lbound(kmin_v),ubound(kmin_v))
      if(allocated(kmin_w))call netcdfrestart_wrt(kmin_w,'kmin_w',lbound(kmin_w),ubound(kmin_w))
      if(allocated(kundermin_t))call netcdfrestart_wrt(kundermin_t,'kundermin_t',lbound(kundermin_t),ubound(kundermin_t))
      if(allocated(kundermin_u))call netcdfrestart_wrt(kundermin_u,'kundermin_u',lbound(kundermin_u),ubound(kundermin_u))
      if(allocated(kundermin_v))call netcdfrestart_wrt(kundermin_v,'kundermin_v',lbound(kundermin_v),ubound(kundermin_v))
      if(allocated(kmerged_t))call netcdfrestart_wrt(kmerged_t,'kmerged_t',lbound(kmerged_t),ubound(kmerged_t))
      if(allocated(kmerged_u))call netcdfrestart_wrt(kmerged_u,'kmerged_u',lbound(kmerged_u),ubound(kmerged_u))
      if(allocated(kmerged_v))call netcdfrestart_wrt(kmerged_v,'kmerged_v',lbound(kmerged_v),ubound(kmerged_v))
      if(allocated(ksl_t))call netcdfrestart_wrt(ksl_t,'ksl_t',lbound(ksl_t),ubound(ksl_t))
      if(allocated(glob_mask_mangrove))call netcdfrestart_wrt(glob_mask_mangrove,'glob_mask_mangrove',lbound(glob_mask_mangrove),ubound(glob_mask_mangrove))
      if(allocated(mask_mangrove_t))call netcdfrestart_wrt(mask_mangrove_t,'mask_mangrove_t',lbound(mask_mangrove_t),ubound(mask_mangrove_t))
      if(allocated(mask_wave_t))call netcdfrestart_wrt(mask_wave_t,'mask_wave_t',lbound(mask_wave_t),ubound(mask_wave_t))
      if(allocated(upwindriver_t))call netcdfrestart_wrt(upwindriver_t,'upwindriver_t',lbound(upwindriver_t),ubound(upwindriver_t))
      if(allocated(upwindwetdry_t))call netcdfrestart_wrt(upwindwetdry_t,'upwindwetdry_t',lbound(upwindwetdry_t),ubound(upwindwetdry_t))
      if(allocated(kmergedr4_u))call netcdfrestart_wrt(kmergedr4_u,'kmergedr4_u',lbound(kmergedr4_u),ubound(kmergedr4_u))
      if(allocated(kmergedr4_v))call netcdfrestart_wrt(kmergedr4_v,'kmergedr4_v',lbound(kmergedr4_v),ubound(kmergedr4_v))
      if(allocated(dsigmerged_u))call netcdfrestart_wrt(dsigmerged_u,'dsigmerged_u',lbound(dsigmerged_u),ubound(dsigmerged_u))
      if(allocated(dsigmerged_v))call netcdfrestart_wrt(dsigmerged_v,'dsigmerged_v',lbound(dsigmerged_v),ubound(dsigmerged_v))
      if(allocated(pgfratio_u))call netcdfrestart_wrt(pgfratio_u,'pgfratio_u',lbound(pgfratio_u),ubound(pgfratio_u))
      if(allocated(pgfratio_v))call netcdfrestart_wrt(pgfratio_v,'pgfratio_v',lbound(pgfratio_v),ubound(pgfratio_v))
      if(allocated(sshr4_w))call netcdfrestart_wrt(sshr4_w,'sshr4_w',lbound(sshr4_w),ubound(sshr4_w))
      if(allocated(wetmask_u))call netcdfrestart_wrt(wetmask_u,'wetmask_u',lbound(wetmask_u),ubound(wetmask_u))
      if(allocated(wetmask_v))call netcdfrestart_wrt(wetmask_v,'wetmask_v',lbound(wetmask_v),ubound(wetmask_v))
      if(allocated(wetmask_t))call netcdfrestart_wrt(wetmask_t,'wetmask_t',lbound(wetmask_t),ubound(wetmask_t))
      if(allocated(botlevmerged_w))call netcdfrestart_wrt(botlevmerged_w,'botlevmerged_w',lbound(botlevmerged_w),ubound(botlevmerged_w))
      if(allocated(maxbotstress_aft_w))call netcdfrestart_wrt(maxbotstress_aft_w,'maxbotstress_aft_w',lbound(maxbotstress_aft_w),ubound(maxbotstress_aft_w))
      if(allocated(maxbotstress_bef_w))call netcdfrestart_wrt(maxbotstress_bef_w,'maxbotstress_bef_w',lbound(maxbotstress_bef_w),ubound(maxbotstress_bef_w))
      if(allocated(maxbotstress_w))call netcdfrestart_wrt(maxbotstress_w,'maxbotstress_w',lbound(maxbotstress_w),ubound(maxbotstress_w))
      if(allocated(stresswave_w))call netcdfrestart_wrt(stresswave_w,'stresswave_w',lbound(stresswave_w),ubound(stresswave_w))
      if(allocated(stressc_w))call netcdfrestart_wrt(stressc_w,'stressc_w',lbound(stressc_w),ubound(stressc_w))
      if(allocated(sqr_hoverg_u))call netcdfrestart_wrt(sqr_hoverg_u,'sqr_hoverg_u',lbound(sqr_hoverg_u),ubound(sqr_hoverg_u))
      if(allocated(sqr_hoverg_v))call netcdfrestart_wrt(sqr_hoverg_v,'sqr_hoverg_v',lbound(sqr_hoverg_v),ubound(sqr_hoverg_v))
      if(allocated(temobc_t))call netcdfrestart_wrt(temobc_t,'temobc_t',lbound(temobc_t),ubound(temobc_t))
      if(allocated(salobc_t))call netcdfrestart_wrt(salobc_t,'salobc_t',lbound(salobc_t),ubound(salobc_t))
      if(allocated(velobc_u))call netcdfrestart_wrt(velobc_u,'velobc_u',lbound(velobc_u),ubound(velobc_u))
      if(allocated(velobc_v))call netcdfrestart_wrt(velobc_v,'velobc_v',lbound(velobc_v),ubound(velobc_v))
      if(allocated(temofl_t))call netcdfrestart_wrt(temofl_t,'temofl_t',lbound(temofl_t),ubound(temofl_t))
      if(allocated(salofl_t))call netcdfrestart_wrt(salofl_t,'salofl_t',lbound(salofl_t),ubound(salofl_t))
      if(allocated(dzofl_t))call netcdfrestart_wrt(dzofl_t,'dzofl_t',lbound(dzofl_t),ubound(dzofl_t))
      if(allocated(velofl_u))call netcdfrestart_wrt(velofl_u,'velofl_u',lbound(velofl_u),ubound(velofl_u))
      if(allocated(velofl_v))call netcdfrestart_wrt(velofl_v,'velofl_v',lbound(velofl_v),ubound(velofl_v))
      if(allocated(tkeofl_w))call netcdfrestart_wrt(tkeofl_w,'tkeofl_w',lbound(tkeofl_w),ubound(tkeofl_w))
      if(allocated(bioofl_t))call netcdfrestart_wrt(bioofl_t,'bioofl_t',lbound(bioofl_t),ubound(bioofl_t))
      if(allocated(dfvofl_w))call netcdfrestart_wrt(dfvofl_w,'dfvofl_w',lbound(dfvofl_w),ubound(dfvofl_w))
      if(allocated(uwindabl_t))call netcdfrestart_wrt(uwindabl_t,'uwindabl_t',lbound(uwindabl_t),ubound(uwindabl_t))
      if(allocated(vwindabl_t))call netcdfrestart_wrt(vwindabl_t,'vwindabl_t',lbound(vwindabl_t),ubound(vwindabl_t))
      if(allocated(w0mofl_w))call netcdfrestart_wrt(w0mofl_w,'w0mofl_w',lbound(w0mofl_w),ubound(w0mofl_w))
      if(allocated(w_keq1_ofl_w))call netcdfrestart_wrt(w_keq1_ofl_w,'w_keq1_ofl_w',lbound(w_keq1_ofl_w),ubound(w_keq1_ofl_w))
      if(allocated(kslofl_t))call netcdfrestart_wrt(kslofl_t,'kslofl_t',lbound(kslofl_t),ubound(kslofl_t))
      if(allocated(ablheight_t))call netcdfrestart_wrt(ablheight_t,'ablheight_t',lbound(ablheight_t),ubound(ablheight_t))
      if(allocated(wwindabl_w))call netcdfrestart_wrt(wwindabl_w,'wwindabl_w',lbound(wwindabl_w),ubound(wwindabl_w))
      if(allocated(kz_abl_w))call netcdfrestart_wrt(kz_abl_w,'kz_abl_w',lbound(kz_abl_w),ubound(kz_abl_w))
      if(allocated(upwzone0_t))call netcdfrestart_wrt(upwzone0_t,'upwzone0_t',lbound(upwzone0_t),ubound(upwzone0_t))
      if(allocated(velstokes_u))call netcdfrestart_wrt(velstokes_u,'velstokes_u',lbound(velstokes_u),ubound(velstokes_u))
      if(allocated(velstokes_v))call netcdfrestart_wrt(velstokes_v,'velstokes_v',lbound(velstokes_v),ubound(velstokes_v))
      if(allocated(velbarstokes_u))call netcdfrestart_wrt(velbarstokes_u,'velbarstokes_u',lbound(velbarstokes_u),ubound(velbarstokes_u))
      if(allocated(velbarstokes_v))call netcdfrestart_wrt(velbarstokes_v,'velbarstokes_v',lbound(velbarstokes_v),ubound(velbarstokes_v))
      if(allocated(nhp1_t))call netcdfrestart_wrt(nhp1_t,'nhp1_t',lbound(nhp1_t),ubound(nhp1_t))
      if(allocated(nhp2_t))call netcdfrestart_wrt(nhp2_t,'nhp2_t',lbound(nhp2_t),ubound(nhp2_t))
      if(allocated(temf_t))call netcdfrestart_wrt(temf_t,'temf_t',lbound(temf_t),ubound(temf_t))
      if(allocated(salf_t))call netcdfrestart_wrt(salf_t,'salf_t',lbound(salf_t),ubound(salf_t))
      if(allocated(temlwf_t))call netcdfrestart_wrt(temlwf_t,'temlwf_t',lbound(temlwf_t),ubound(temlwf_t))
      if(allocated(sallwf_t))call netcdfrestart_wrt(sallwf_t,'sallwf_t',lbound(sallwf_t),ubound(sallwf_t))
      if(allocated(sshlwf_w))call netcdfrestart_wrt(sshlwf_w,'sshlwf_w',lbound(sshlwf_w),ubound(sshlwf_w))
      if(allocated(t_wave_t))call netcdfrestart_wrt(t_wave_t,'t_wave_t',lbound(t_wave_t),ubound(t_wave_t))
      if(allocated(hs_wave_t))call netcdfrestart_wrt(hs_wave_t,'hs_wave_t',lbound(hs_wave_t),ubound(hs_wave_t))
      if(allocated(hsw_wave_t))call netcdfrestart_wrt(hsw_wave_t,'hsw_wave_t',lbound(hsw_wave_t),ubound(hsw_wave_t))
      if(allocated(foc_wave_t))call netcdfrestart_wrt(foc_wave_t,'foc_wave_t',lbound(foc_wave_t),ubound(foc_wave_t))
      if(allocated(k_wave_t))call netcdfrestart_wrt(k_wave_t,'k_wave_t',lbound(k_wave_t),ubound(k_wave_t))
      if(allocated(kx_wave_t))call netcdfrestart_wrt(kx_wave_t,'kx_wave_t',lbound(kx_wave_t),ubound(kx_wave_t))
      if(allocated(ky_wave_t))call netcdfrestart_wrt(ky_wave_t,'ky_wave_t',lbound(ky_wave_t),ubound(ky_wave_t))
      if(allocated(twox_wave_t))call netcdfrestart_wrt(twox_wave_t,'twox_wave_t',lbound(twox_wave_t),ubound(twox_wave_t))
      if(allocated(twoy_wave_t))call netcdfrestart_wrt(twoy_wave_t,'twoy_wave_t',lbound(twoy_wave_t),ubound(twoy_wave_t))
      if(allocated(tawx_wave_t))call netcdfrestart_wrt(tawx_wave_t,'tawx_wave_t',lbound(tawx_wave_t),ubound(tawx_wave_t))
      if(allocated(tawy_wave_t))call netcdfrestart_wrt(tawy_wave_t,'tawy_wave_t',lbound(tawy_wave_t),ubound(tawy_wave_t))
      if(allocated(usf_wave_t))call netcdfrestart_wrt(usf_wave_t,'usf_wave_t',lbound(usf_wave_t),ubound(usf_wave_t))
      if(allocated(vsf_wave_t))call netcdfrestart_wrt(vsf_wave_t,'vsf_wave_t',lbound(vsf_wave_t),ubound(vsf_wave_t))
      if(allocated(dir_wave_t))call netcdfrestart_wrt(dir_wave_t,'dir_wave_t',lbound(dir_wave_t),ubound(dir_wave_t))
      if(allocated(uss_wave_t))call netcdfrestart_wrt(uss_wave_t,'uss_wave_t',lbound(uss_wave_t),ubound(uss_wave_t))
      if(allocated(j_wave_t))call netcdfrestart_wrt(j_wave_t,'j_wave_t',lbound(j_wave_t),ubound(j_wave_t))
      if(allocated(vss_wave_t))call netcdfrestart_wrt(vss_wave_t,'vss_wave_t',lbound(vss_wave_t),ubound(vss_wave_t))
      if(allocated(ubw))call netcdfrestart_wrt(ubw,'ubw',lbound(ubw),ubound(ubw))
      if(allocated(fw))call netcdfrestart_wrt(fw,'fw',lbound(fw),ubound(fw))
      if(allocated(dpt_wave_t))call netcdfrestart_wrt(dpt_wave_t,'dpt_wave_t',lbound(dpt_wave_t),ubound(dpt_wave_t))
      if(allocated(wstresb_u))call netcdfrestart_wrt(wstresb_u,'wstresb_u',lbound(wstresb_u),ubound(wstresb_u))
      if(allocated(wstresb_v))call netcdfrestart_wrt(wstresb_v,'wstresb_v',lbound(wstresb_v),ubound(wstresb_v))
      if(allocated(ij2ww3_i))call netcdfrestart_wrt(ij2ww3_i,'ij2ww3_i',lbound(ij2ww3_i),ubound(ij2ww3_i))
      if(allocated(ij2ww3_j))call netcdfrestart_wrt(ij2ww3_j,'ij2ww3_j',lbound(ij2ww3_j),ubound(ij2ww3_j))
      if(allocated(ij2ww3_teta))call netcdfrestart_wrt(ij2ww3_teta,'ij2ww3_teta',lbound(ij2ww3_teta),ubound(ij2ww3_teta))
      if(allocated(slhf_aver_w))call netcdfrestart_wrt(slhf_aver_w,'slhf_aver_w',lbound(slhf_aver_w),ubound(slhf_aver_w))
      if(allocated(sshf_aver_w))call netcdfrestart_wrt(sshf_aver_w,'sshf_aver_w',lbound(sshf_aver_w),ubound(sshf_aver_w))
      if(allocated(snsf_aver_w))call netcdfrestart_wrt(snsf_aver_w,'snsf_aver_w',lbound(snsf_aver_w),ubound(snsf_aver_w))
      if(allocated(ssr_aver_w))call netcdfrestart_wrt(ssr_aver_w,'ssr_aver_w',lbound(ssr_aver_w),ubound(ssr_aver_w))
      if(allocated(precipi_aver_w))call netcdfrestart_wrt(precipi_aver_w,'precipi_aver_w',lbound(precipi_aver_w),ubound(precipi_aver_w))
      if(allocated(wstress_aver_u))call netcdfrestart_wrt(wstress_aver_u,'wstress_aver_u',lbound(wstress_aver_u),ubound(wstress_aver_u))
      if(allocated(wstress_aver_v))call netcdfrestart_wrt(wstress_aver_v,'wstress_aver_v',lbound(wstress_aver_v),ubound(wstress_aver_v))
      if(allocated(hsedofl_t))call netcdfrestart_wrt(hsedofl_t,'hsedofl_t',lbound(hsedofl_t),ubound(hsedofl_t))
      write(9)zone1saltflux_w
      write(9)zone1tempflux_w
      write(9)zone1waterflux_w
      write(9)zone1_nlayer
      write(9)zone1_max
      write(9)zone1_u_max
      write(9)zone1_v_max
      write(9)zone1_inv_dz
      write(9)zone1_stretch_dz
      if(allocated(zone1saltflux_glb))call netcdfrestart_wrt(zone1saltflux_glb,'zone1saltflux_glb',lbound(zone1saltflux_glb),ubound(zone1saltflux_glb))
      if(allocated(zone1saltflux_u))call netcdfrestart_wrt(zone1saltflux_u,'zone1saltflux_u',lbound(zone1saltflux_u),ubound(zone1saltflux_u))
      if(allocated(zone1saltflux_v))call netcdfrestart_wrt(zone1saltflux_v,'zone1saltflux_v',lbound(zone1saltflux_v),ubound(zone1saltflux_v))
      if(allocated(zone1tempflux_glb))call netcdfrestart_wrt(zone1tempflux_glb,'zone1tempflux_glb',lbound(zone1tempflux_glb),ubound(zone1tempflux_glb))
      if(allocated(zone1tempflux_u))call netcdfrestart_wrt(zone1tempflux_u,'zone1tempflux_u',lbound(zone1tempflux_u),ubound(zone1tempflux_u))
      if(allocated(zone1tempflux_v))call netcdfrestart_wrt(zone1tempflux_v,'zone1tempflux_v',lbound(zone1tempflux_v),ubound(zone1tempflux_v))
      if(allocated(zone1waterflux_glb))call netcdfrestart_wrt(zone1waterflux_glb,'zone1waterflux_glb',lbound(zone1waterflux_glb),ubound(zone1waterflux_glb))
      if(allocated(zone1waterflux_u))call netcdfrestart_wrt(zone1waterflux_u,'zone1waterflux_u',lbound(zone1waterflux_u),ubound(zone1waterflux_u))
      if(allocated(zone1waterflux_v))call netcdfrestart_wrt(zone1waterflux_v,'zone1waterflux_v',lbound(zone1waterflux_v),ubound(zone1waterflux_v))
      if(allocated(zone1saltflux_glb_in))call netcdfrestart_wrt(zone1saltflux_glb_in,'zone1saltflux_glb_in',lbound(zone1saltflux_glb_in),ubound(zone1saltflux_glb_in))
      if(allocated(zone1saltflux_u_in))call netcdfrestart_wrt(zone1saltflux_u_in,'zone1saltflux_u_in',lbound(zone1saltflux_u_in),ubound(zone1saltflux_u_in))
      if(allocated(zone1saltflux_v_in))call netcdfrestart_wrt(zone1saltflux_v_in,'zone1saltflux_v_in',lbound(zone1saltflux_v_in),ubound(zone1saltflux_v_in))
      if(allocated(zone1tempflux_glb_in))call netcdfrestart_wrt(zone1tempflux_glb_in,'zone1tempflux_glb_in',lbound(zone1tempflux_glb_in),ubound(zone1tempflux_glb_in))
      if(allocated(zone1tempflux_u_in))call netcdfrestart_wrt(zone1tempflux_u_in,'zone1tempflux_u_in',lbound(zone1tempflux_u_in),ubound(zone1tempflux_u_in))
      if(allocated(zone1tempflux_v_in))call netcdfrestart_wrt(zone1tempflux_v_in,'zone1tempflux_v_in',lbound(zone1tempflux_v_in),ubound(zone1tempflux_v_in))
      if(allocated(zone1waterflux_glb_in))call netcdfrestart_wrt(zone1waterflux_glb_in,'zone1waterflux_glb_in',lbound(zone1waterflux_glb_in),ubound(zone1waterflux_glb_in))
      if(allocated(zone1waterflux_u_in))call netcdfrestart_wrt(zone1waterflux_u_in,'zone1waterflux_u_in',lbound(zone1waterflux_u_in),ubound(zone1waterflux_u_in))
      if(allocated(zone1waterflux_v_in))call netcdfrestart_wrt(zone1waterflux_v_in,'zone1waterflux_v_in',lbound(zone1waterflux_v_in),ubound(zone1waterflux_v_in))
      if(allocated(zone1saltflux_glb_out))call netcdfrestart_wrt(zone1saltflux_glb_out,'zone1saltflux_glb_out',lbound(zone1saltflux_glb_out),ubound(zone1saltflux_glb_out))
      if(allocated(zone1saltflux_u_out))call netcdfrestart_wrt(zone1saltflux_u_out,'zone1saltflux_u_out',lbound(zone1saltflux_u_out),ubound(zone1saltflux_u_out))
      if(allocated(zone1saltflux_v_out))call netcdfrestart_wrt(zone1saltflux_v_out,'zone1saltflux_v_out',lbound(zone1saltflux_v_out),ubound(zone1saltflux_v_out))
      if(allocated(zone1tempflux_glb_out))call netcdfrestart_wrt(zone1tempflux_glb_out,'zone1tempflux_glb_out',lbound(zone1tempflux_glb_out),ubound(zone1tempflux_glb_out))
      if(allocated(zone1tempflux_u_out))call netcdfrestart_wrt(zone1tempflux_u_out,'zone1tempflux_u_out',lbound(zone1tempflux_u_out),ubound(zone1tempflux_u_out))
      if(allocated(zone1tempflux_v_out))call netcdfrestart_wrt(zone1tempflux_v_out,'zone1tempflux_v_out',lbound(zone1tempflux_v_out),ubound(zone1tempflux_v_out))
      if(allocated(zone1waterflux_glb_out))call netcdfrestart_wrt(zone1waterflux_glb_out,'zone1waterflux_glb_out',lbound(zone1waterflux_glb_out),ubound(zone1waterflux_glb_out))
      if(allocated(zone1waterflux_u_out))call netcdfrestart_wrt(zone1waterflux_u_out,'zone1waterflux_u_out',lbound(zone1waterflux_u_out),ubound(zone1waterflux_u_out))
      if(allocated(zone1waterflux_v_out))call netcdfrestart_wrt(zone1waterflux_v_out,'zone1waterflux_v_out',lbound(zone1waterflux_v_out),ubound(zone1waterflux_v_out))
      if(allocated(zone1tempcumul_glb))write(9)zone1tempcumul_glb
      if(allocated(zone1tempcumul_loc))write(9)zone1tempcumul_loc
      if(allocated(zone1tempmasst0))write(9)zone1tempmasst0
      if(allocated(zone1saltcumul_glb))write(9)zone1saltcumul_glb
      if(allocated(zone1saltcumul_loc))write(9)zone1saltcumul_loc
      if(allocated(zone1saltmasst0))write(9)zone1saltmasst0
      if(allocated(zone1watercumul_glb))write(9)zone1watercumul_glb
      if(allocated(zone1watercumul_loc))write(9)zone1watercumul_loc
      if(allocated(zone1watermasst0))write(9)zone1watermasst0
      if(allocated(zone1_mask))call netcdfrestart_wrt(zone1_mask,'zone1_mask',lbound(zone1_mask),ubound(zone1_mask))
      if(allocated(zone1_flux_u_node))call netcdfrestart_wrt(zone1_flux_u_node,'zone1_flux_u_node',lbound(zone1_flux_u_node),ubound(zone1_flux_u_node))
      if(allocated(zone1_flux_v_node))call netcdfrestart_wrt(zone1_flux_v_node,'zone1_flux_v_node',lbound(zone1_flux_v_node),ubound(zone1_flux_v_node))
      write(9)zone2saltflux_w
      write(9)zone2tempflux_w
      write(9)zone2waterflux_w
      write(9)zone2_nlayer
      write(9)zone2_max
      write(9)zone2_u_max
      write(9)zone2_v_max
      write(9)zone2_inv_dz
      write(9)zone2_stretch_dz
      if(allocated(zone2saltflux_glb))call netcdfrestart_wrt(zone2saltflux_glb,'zone2saltflux_glb',lbound(zone2saltflux_glb),ubound(zone2saltflux_glb))
      if(allocated(zone2saltflux_u))call netcdfrestart_wrt(zone2saltflux_u,'zone2saltflux_u',lbound(zone2saltflux_u),ubound(zone2saltflux_u))
      if(allocated(zone2saltflux_v))call netcdfrestart_wrt(zone2saltflux_v,'zone2saltflux_v',lbound(zone2saltflux_v),ubound(zone2saltflux_v))
      if(allocated(zone2tempflux_glb))call netcdfrestart_wrt(zone2tempflux_glb,'zone2tempflux_glb',lbound(zone2tempflux_glb),ubound(zone2tempflux_glb))
      if(allocated(zone2tempflux_u))call netcdfrestart_wrt(zone2tempflux_u,'zone2tempflux_u',lbound(zone2tempflux_u),ubound(zone2tempflux_u))
      if(allocated(zone2tempflux_v))call netcdfrestart_wrt(zone2tempflux_v,'zone2tempflux_v',lbound(zone2tempflux_v),ubound(zone2tempflux_v))
      if(allocated(zone2waterflux_glb))call netcdfrestart_wrt(zone2waterflux_glb,'zone2waterflux_glb',lbound(zone2waterflux_glb),ubound(zone2waterflux_glb))
      if(allocated(zone2waterflux_u))call netcdfrestart_wrt(zone2waterflux_u,'zone2waterflux_u',lbound(zone2waterflux_u),ubound(zone2waterflux_u))
      if(allocated(zone2waterflux_v))call netcdfrestart_wrt(zone2waterflux_v,'zone2waterflux_v',lbound(zone2waterflux_v),ubound(zone2waterflux_v))
      if(allocated(zone2saltflux_glb_in))call netcdfrestart_wrt(zone2saltflux_glb_in,'zone2saltflux_glb_in',lbound(zone2saltflux_glb_in),ubound(zone2saltflux_glb_in))
      if(allocated(zone2saltflux_u_in))call netcdfrestart_wrt(zone2saltflux_u_in,'zone2saltflux_u_in',lbound(zone2saltflux_u_in),ubound(zone2saltflux_u_in))
      if(allocated(zone2saltflux_v_in))call netcdfrestart_wrt(zone2saltflux_v_in,'zone2saltflux_v_in',lbound(zone2saltflux_v_in),ubound(zone2saltflux_v_in))
      if(allocated(zone2tempflux_glb_in))call netcdfrestart_wrt(zone2tempflux_glb_in,'zone2tempflux_glb_in',lbound(zone2tempflux_glb_in),ubound(zone2tempflux_glb_in))
      if(allocated(zone2tempflux_u_in))call netcdfrestart_wrt(zone2tempflux_u_in,'zone2tempflux_u_in',lbound(zone2tempflux_u_in),ubound(zone2tempflux_u_in))
      if(allocated(zone2tempflux_v_in))call netcdfrestart_wrt(zone2tempflux_v_in,'zone2tempflux_v_in',lbound(zone2tempflux_v_in),ubound(zone2tempflux_v_in))
      if(allocated(zone2waterflux_glb_in))call netcdfrestart_wrt(zone2waterflux_glb_in,'zone2waterflux_glb_in',lbound(zone2waterflux_glb_in),ubound(zone2waterflux_glb_in))
      if(allocated(zone2waterflux_u_in))call netcdfrestart_wrt(zone2waterflux_u_in,'zone2waterflux_u_in',lbound(zone2waterflux_u_in),ubound(zone2waterflux_u_in))
      if(allocated(zone2waterflux_v_in))call netcdfrestart_wrt(zone2waterflux_v_in,'zone2waterflux_v_in',lbound(zone2waterflux_v_in),ubound(zone2waterflux_v_in))
      if(allocated(zone2saltflux_glb_out))call netcdfrestart_wrt(zone2saltflux_glb_out,'zone2saltflux_glb_out',lbound(zone2saltflux_glb_out),ubound(zone2saltflux_glb_out))
      if(allocated(zone2saltflux_u_out))call netcdfrestart_wrt(zone2saltflux_u_out,'zone2saltflux_u_out',lbound(zone2saltflux_u_out),ubound(zone2saltflux_u_out))
      if(allocated(zone2saltflux_v_out))call netcdfrestart_wrt(zone2saltflux_v_out,'zone2saltflux_v_out',lbound(zone2saltflux_v_out),ubound(zone2saltflux_v_out))
      if(allocated(zone2tempflux_glb_out))call netcdfrestart_wrt(zone2tempflux_glb_out,'zone2tempflux_glb_out',lbound(zone2tempflux_glb_out),ubound(zone2tempflux_glb_out))
      if(allocated(zone2tempflux_u_out))call netcdfrestart_wrt(zone2tempflux_u_out,'zone2tempflux_u_out',lbound(zone2tempflux_u_out),ubound(zone2tempflux_u_out))
      if(allocated(zone2tempflux_v_out))call netcdfrestart_wrt(zone2tempflux_v_out,'zone2tempflux_v_out',lbound(zone2tempflux_v_out),ubound(zone2tempflux_v_out))
      if(allocated(zone2waterflux_glb_out))call netcdfrestart_wrt(zone2waterflux_glb_out,'zone2waterflux_glb_out',lbound(zone2waterflux_glb_out),ubound(zone2waterflux_glb_out))
      if(allocated(zone2waterflux_u_out))call netcdfrestart_wrt(zone2waterflux_u_out,'zone2waterflux_u_out',lbound(zone2waterflux_u_out),ubound(zone2waterflux_u_out))
      if(allocated(zone2waterflux_v_out))call netcdfrestart_wrt(zone2waterflux_v_out,'zone2waterflux_v_out',lbound(zone2waterflux_v_out),ubound(zone2waterflux_v_out))
      if(allocated(zone2tempcumul_glb))write(9)zone2tempcumul_glb
      if(allocated(zone2tempcumul_loc))write(9)zone2tempcumul_loc
      if(allocated(zone2tempmasst0))write(9)zone2tempmasst0
      if(allocated(zone2saltcumul_glb))write(9)zone2saltcumul_glb
      if(allocated(zone2saltcumul_loc))write(9)zone2saltcumul_loc
      if(allocated(zone2saltmasst0))write(9)zone2saltmasst0
      if(allocated(zone2watercumul_glb))write(9)zone2watercumul_glb
      if(allocated(zone2watercumul_loc))write(9)zone2watercumul_loc
      if(allocated(zone2watermasst0))write(9)zone2watermasst0
      if(allocated(zone2_mask))call netcdfrestart_wrt(zone2_mask,'zone2_mask',lbound(zone2_mask),ubound(zone2_mask))
      if(allocated(zone2_flux_u_node))call netcdfrestart_wrt(zone2_flux_u_node,'zone2_flux_u_node',lbound(zone2_flux_u_node),ubound(zone2_flux_u_node))
      if(allocated(zone2_flux_v_node))call netcdfrestart_wrt(zone2_flux_v_node,'zone2_flux_v_node',lbound(zone2_flux_v_node),ubound(zone2_flux_v_node))
      write(9)zone3saltflux_w
      write(9)zone3tempflux_w
      write(9)zone3waterflux_w
      write(9)zone3_nlayer
      write(9)zone3_max
      write(9)zone3_u_max
      write(9)zone3_v_max
      write(9)zone3_inv_dz
      write(9)zone3_stretch_dz
      if(allocated(zone3saltflux_glb))call netcdfrestart_wrt(zone3saltflux_glb,'zone3saltflux_glb',lbound(zone3saltflux_glb),ubound(zone3saltflux_glb))
      if(allocated(zone3saltflux_u))call netcdfrestart_wrt(zone3saltflux_u,'zone3saltflux_u',lbound(zone3saltflux_u),ubound(zone3saltflux_u))
      if(allocated(zone3saltflux_v))call netcdfrestart_wrt(zone3saltflux_v,'zone3saltflux_v',lbound(zone3saltflux_v),ubound(zone3saltflux_v))
      if(allocated(zone3tempflux_glb))call netcdfrestart_wrt(zone3tempflux_glb,'zone3tempflux_glb',lbound(zone3tempflux_glb),ubound(zone3tempflux_glb))
      if(allocated(zone3tempflux_u))call netcdfrestart_wrt(zone3tempflux_u,'zone3tempflux_u',lbound(zone3tempflux_u),ubound(zone3tempflux_u))
      if(allocated(zone3tempflux_v))call netcdfrestart_wrt(zone3tempflux_v,'zone3tempflux_v',lbound(zone3tempflux_v),ubound(zone3tempflux_v))
      if(allocated(zone3waterflux_glb))call netcdfrestart_wrt(zone3waterflux_glb,'zone3waterflux_glb',lbound(zone3waterflux_glb),ubound(zone3waterflux_glb))
      if(allocated(zone3waterflux_u))call netcdfrestart_wrt(zone3waterflux_u,'zone3waterflux_u',lbound(zone3waterflux_u),ubound(zone3waterflux_u))
      if(allocated(zone3waterflux_v))call netcdfrestart_wrt(zone3waterflux_v,'zone3waterflux_v',lbound(zone3waterflux_v),ubound(zone3waterflux_v))
      if(allocated(zone3saltflux_glb_in))call netcdfrestart_wrt(zone3saltflux_glb_in,'zone3saltflux_glb_in',lbound(zone3saltflux_glb_in),ubound(zone3saltflux_glb_in))
      if(allocated(zone3saltflux_u_in))call netcdfrestart_wrt(zone3saltflux_u_in,'zone3saltflux_u_in',lbound(zone3saltflux_u_in),ubound(zone3saltflux_u_in))
      if(allocated(zone3saltflux_v_in))call netcdfrestart_wrt(zone3saltflux_v_in,'zone3saltflux_v_in',lbound(zone3saltflux_v_in),ubound(zone3saltflux_v_in))
      if(allocated(zone3tempflux_glb_in))call netcdfrestart_wrt(zone3tempflux_glb_in,'zone3tempflux_glb_in',lbound(zone3tempflux_glb_in),ubound(zone3tempflux_glb_in))
      if(allocated(zone3tempflux_u_in))call netcdfrestart_wrt(zone3tempflux_u_in,'zone3tempflux_u_in',lbound(zone3tempflux_u_in),ubound(zone3tempflux_u_in))
      if(allocated(zone3tempflux_v_in))call netcdfrestart_wrt(zone3tempflux_v_in,'zone3tempflux_v_in',lbound(zone3tempflux_v_in),ubound(zone3tempflux_v_in))
      if(allocated(zone3waterflux_glb_in))call netcdfrestart_wrt(zone3waterflux_glb_in,'zone3waterflux_glb_in',lbound(zone3waterflux_glb_in),ubound(zone3waterflux_glb_in))
      if(allocated(zone3waterflux_u_in))call netcdfrestart_wrt(zone3waterflux_u_in,'zone3waterflux_u_in',lbound(zone3waterflux_u_in),ubound(zone3waterflux_u_in))
      if(allocated(zone3waterflux_v_in))call netcdfrestart_wrt(zone3waterflux_v_in,'zone3waterflux_v_in',lbound(zone3waterflux_v_in),ubound(zone3waterflux_v_in))
      if(allocated(zone3saltflux_glb_out))call netcdfrestart_wrt(zone3saltflux_glb_out,'zone3saltflux_glb_out',lbound(zone3saltflux_glb_out),ubound(zone3saltflux_glb_out))
      if(allocated(zone3saltflux_u_out))call netcdfrestart_wrt(zone3saltflux_u_out,'zone3saltflux_u_out',lbound(zone3saltflux_u_out),ubound(zone3saltflux_u_out))
      if(allocated(zone3saltflux_v_out))call netcdfrestart_wrt(zone3saltflux_v_out,'zone3saltflux_v_out',lbound(zone3saltflux_v_out),ubound(zone3saltflux_v_out))
      if(allocated(zone3tempflux_glb_out))call netcdfrestart_wrt(zone3tempflux_glb_out,'zone3tempflux_glb_out',lbound(zone3tempflux_glb_out),ubound(zone3tempflux_glb_out))
      if(allocated(zone3tempflux_u_out))call netcdfrestart_wrt(zone3tempflux_u_out,'zone3tempflux_u_out',lbound(zone3tempflux_u_out),ubound(zone3tempflux_u_out))
      if(allocated(zone3tempflux_v_out))call netcdfrestart_wrt(zone3tempflux_v_out,'zone3tempflux_v_out',lbound(zone3tempflux_v_out),ubound(zone3tempflux_v_out))
      if(allocated(zone3waterflux_glb_out))call netcdfrestart_wrt(zone3waterflux_glb_out,'zone3waterflux_glb_out',lbound(zone3waterflux_glb_out),ubound(zone3waterflux_glb_out))
      if(allocated(zone3waterflux_u_out))call netcdfrestart_wrt(zone3waterflux_u_out,'zone3waterflux_u_out',lbound(zone3waterflux_u_out),ubound(zone3waterflux_u_out))
      if(allocated(zone3waterflux_v_out))call netcdfrestart_wrt(zone3waterflux_v_out,'zone3waterflux_v_out',lbound(zone3waterflux_v_out),ubound(zone3waterflux_v_out))
      if(allocated(zone3tempcumul_glb))write(9)zone3tempcumul_glb
      if(allocated(zone3tempcumul_loc))write(9)zone3tempcumul_loc
      if(allocated(zone3tempmasst0))write(9)zone3tempmasst0
      if(allocated(zone3saltcumul_glb))write(9)zone3saltcumul_glb
      if(allocated(zone3saltcumul_loc))write(9)zone3saltcumul_loc
      if(allocated(zone3saltmasst0))write(9)zone3saltmasst0
      if(allocated(zone3watercumul_glb))write(9)zone3watercumul_glb
      if(allocated(zone3watercumul_loc))write(9)zone3watercumul_loc
      if(allocated(zone3watermasst0))write(9)zone3watermasst0
      if(allocated(zone3_mask))call netcdfrestart_wrt(zone3_mask,'zone3_mask',lbound(zone3_mask),ubound(zone3_mask))
      if(allocated(zone3_flux_u_node))call netcdfrestart_wrt(zone3_flux_u_node,'zone3_flux_u_node',lbound(zone3_flux_u_node),ubound(zone3_flux_u_node))
      if(allocated(zone3_flux_v_node))call netcdfrestart_wrt(zone3_flux_v_node,'zone3_flux_v_node',lbound(zone3_flux_v_node),ubound(zone3_flux_v_node))
      write(9)zone4saltflux_w
      write(9)zone4tempflux_w
      write(9)zone4waterflux_w
      write(9)zone4_nlayer
      write(9)zone4_max
      write(9)zone4_u_max
      write(9)zone4_v_max
      write(9)zone4_inv_dz
      write(9)zone4_stretch_dz
      if(allocated(zone4saltflux_glb))call netcdfrestart_wrt(zone4saltflux_glb,'zone4saltflux_glb',lbound(zone4saltflux_glb),ubound(zone4saltflux_glb))
      if(allocated(zone4saltflux_u))call netcdfrestart_wrt(zone4saltflux_u,'zone4saltflux_u',lbound(zone4saltflux_u),ubound(zone4saltflux_u))
      if(allocated(zone4saltflux_v))call netcdfrestart_wrt(zone4saltflux_v,'zone4saltflux_v',lbound(zone4saltflux_v),ubound(zone4saltflux_v))
      if(allocated(zone4tempflux_glb))call netcdfrestart_wrt(zone4tempflux_glb,'zone4tempflux_glb',lbound(zone4tempflux_glb),ubound(zone4tempflux_glb))
      if(allocated(zone4tempflux_u))call netcdfrestart_wrt(zone4tempflux_u,'zone4tempflux_u',lbound(zone4tempflux_u),ubound(zone4tempflux_u))
      if(allocated(zone4tempflux_v))call netcdfrestart_wrt(zone4tempflux_v,'zone4tempflux_v',lbound(zone4tempflux_v),ubound(zone4tempflux_v))
      if(allocated(zone4waterflux_glb))call netcdfrestart_wrt(zone4waterflux_glb,'zone4waterflux_glb',lbound(zone4waterflux_glb),ubound(zone4waterflux_glb))
      if(allocated(zone4waterflux_u))call netcdfrestart_wrt(zone4waterflux_u,'zone4waterflux_u',lbound(zone4waterflux_u),ubound(zone4waterflux_u))
      if(allocated(zone4waterflux_v))call netcdfrestart_wrt(zone4waterflux_v,'zone4waterflux_v',lbound(zone4waterflux_v),ubound(zone4waterflux_v))
      if(allocated(zone4saltflux_glb_in))call netcdfrestart_wrt(zone4saltflux_glb_in,'zone4saltflux_glb_in',lbound(zone4saltflux_glb_in),ubound(zone4saltflux_glb_in))
      if(allocated(zone4saltflux_u_in))call netcdfrestart_wrt(zone4saltflux_u_in,'zone4saltflux_u_in',lbound(zone4saltflux_u_in),ubound(zone4saltflux_u_in))
      if(allocated(zone4saltflux_v_in))call netcdfrestart_wrt(zone4saltflux_v_in,'zone4saltflux_v_in',lbound(zone4saltflux_v_in),ubound(zone4saltflux_v_in))
      if(allocated(zone4tempflux_glb_in))call netcdfrestart_wrt(zone4tempflux_glb_in,'zone4tempflux_glb_in',lbound(zone4tempflux_glb_in),ubound(zone4tempflux_glb_in))
      if(allocated(zone4tempflux_u_in))call netcdfrestart_wrt(zone4tempflux_u_in,'zone4tempflux_u_in',lbound(zone4tempflux_u_in),ubound(zone4tempflux_u_in))
      if(allocated(zone4tempflux_v_in))call netcdfrestart_wrt(zone4tempflux_v_in,'zone4tempflux_v_in',lbound(zone4tempflux_v_in),ubound(zone4tempflux_v_in))
      if(allocated(zone4waterflux_glb_in))call netcdfrestart_wrt(zone4waterflux_glb_in,'zone4waterflux_glb_in',lbound(zone4waterflux_glb_in),ubound(zone4waterflux_glb_in))
      if(allocated(zone4waterflux_u_in))call netcdfrestart_wrt(zone4waterflux_u_in,'zone4waterflux_u_in',lbound(zone4waterflux_u_in),ubound(zone4waterflux_u_in))
      if(allocated(zone4waterflux_v_in))call netcdfrestart_wrt(zone4waterflux_v_in,'zone4waterflux_v_in',lbound(zone4waterflux_v_in),ubound(zone4waterflux_v_in))
      if(allocated(zone4saltflux_glb_out))call netcdfrestart_wrt(zone4saltflux_glb_out,'zone4saltflux_glb_out',lbound(zone4saltflux_glb_out),ubound(zone4saltflux_glb_out))
      if(allocated(zone4saltflux_u_out))call netcdfrestart_wrt(zone4saltflux_u_out,'zone4saltflux_u_out',lbound(zone4saltflux_u_out),ubound(zone4saltflux_u_out))
      if(allocated(zone4saltflux_v_out))call netcdfrestart_wrt(zone4saltflux_v_out,'zone4saltflux_v_out',lbound(zone4saltflux_v_out),ubound(zone4saltflux_v_out))
      if(allocated(zone4tempflux_glb_out))call netcdfrestart_wrt(zone4tempflux_glb_out,'zone4tempflux_glb_out',lbound(zone4tempflux_glb_out),ubound(zone4tempflux_glb_out))
      if(allocated(zone4tempflux_u_out))call netcdfrestart_wrt(zone4tempflux_u_out,'zone4tempflux_u_out',lbound(zone4tempflux_u_out),ubound(zone4tempflux_u_out))
      if(allocated(zone4tempflux_v_out))call netcdfrestart_wrt(zone4tempflux_v_out,'zone4tempflux_v_out',lbound(zone4tempflux_v_out),ubound(zone4tempflux_v_out))
      if(allocated(zone4waterflux_glb_out))call netcdfrestart_wrt(zone4waterflux_glb_out,'zone4waterflux_glb_out',lbound(zone4waterflux_glb_out),ubound(zone4waterflux_glb_out))
      if(allocated(zone4waterflux_u_out))call netcdfrestart_wrt(zone4waterflux_u_out,'zone4waterflux_u_out',lbound(zone4waterflux_u_out),ubound(zone4waterflux_u_out))
      if(allocated(zone4waterflux_v_out))call netcdfrestart_wrt(zone4waterflux_v_out,'zone4waterflux_v_out',lbound(zone4waterflux_v_out),ubound(zone4waterflux_v_out))
      if(allocated(zone4tempcumul_glb))write(9)zone4tempcumul_glb
      if(allocated(zone4tempcumul_loc))write(9)zone4tempcumul_loc
      if(allocated(zone4tempmasst0))write(9)zone4tempmasst0
      if(allocated(zone4saltcumul_glb))write(9)zone4saltcumul_glb
      if(allocated(zone4saltcumul_loc))write(9)zone4saltcumul_loc
      if(allocated(zone4saltmasst0))write(9)zone4saltmasst0
      if(allocated(zone4watercumul_glb))write(9)zone4watercumul_glb
      if(allocated(zone4watercumul_loc))write(9)zone4watercumul_loc
      if(allocated(zone4watermasst0))write(9)zone4watermasst0
      if(allocated(zone4_mask))call netcdfrestart_wrt(zone4_mask,'zone4_mask',lbound(zone4_mask),ubound(zone4_mask))
      if(allocated(zone4_flux_u_node))call netcdfrestart_wrt(zone4_flux_u_node,'zone4_flux_u_node',lbound(zone4_flux_u_node),ubound(zone4_flux_u_node))
      if(allocated(zone4_flux_v_node))call netcdfrestart_wrt(zone4_flux_v_node,'zone4_flux_v_node',lbound(zone4_flux_v_node),ubound(zone4_flux_v_node))
      write(9)zone5saltflux_w
      write(9)zone5tempflux_w
      write(9)zone5waterflux_w
      write(9)zone5_nlayer
      write(9)zone5_max
      write(9)zone5_u_max
      write(9)zone5_v_max
      write(9)zone5_inv_dz
      write(9)zone5_stretch_dz
      if(allocated(zone5saltflux_glb))call netcdfrestart_wrt(zone5saltflux_glb,'zone5saltflux_glb',lbound(zone5saltflux_glb),ubound(zone5saltflux_glb))
      if(allocated(zone5saltflux_u))call netcdfrestart_wrt(zone5saltflux_u,'zone5saltflux_u',lbound(zone5saltflux_u),ubound(zone5saltflux_u))
      if(allocated(zone5saltflux_v))call netcdfrestart_wrt(zone5saltflux_v,'zone5saltflux_v',lbound(zone5saltflux_v),ubound(zone5saltflux_v))
      if(allocated(zone5tempflux_glb))call netcdfrestart_wrt(zone5tempflux_glb,'zone5tempflux_glb',lbound(zone5tempflux_glb),ubound(zone5tempflux_glb))
      if(allocated(zone5tempflux_u))call netcdfrestart_wrt(zone5tempflux_u,'zone5tempflux_u',lbound(zone5tempflux_u),ubound(zone5tempflux_u))
      if(allocated(zone5tempflux_v))call netcdfrestart_wrt(zone5tempflux_v,'zone5tempflux_v',lbound(zone5tempflux_v),ubound(zone5tempflux_v))
      if(allocated(zone5waterflux_glb))call netcdfrestart_wrt(zone5waterflux_glb,'zone5waterflux_glb',lbound(zone5waterflux_glb),ubound(zone5waterflux_glb))
      if(allocated(zone5waterflux_u))call netcdfrestart_wrt(zone5waterflux_u,'zone5waterflux_u',lbound(zone5waterflux_u),ubound(zone5waterflux_u))
      if(allocated(zone5waterflux_v))call netcdfrestart_wrt(zone5waterflux_v,'zone5waterflux_v',lbound(zone5waterflux_v),ubound(zone5waterflux_v))
      if(allocated(zone5saltflux_glb_in))call netcdfrestart_wrt(zone5saltflux_glb_in,'zone5saltflux_glb_in',lbound(zone5saltflux_glb_in),ubound(zone5saltflux_glb_in))
      if(allocated(zone5saltflux_u_in))call netcdfrestart_wrt(zone5saltflux_u_in,'zone5saltflux_u_in',lbound(zone5saltflux_u_in),ubound(zone5saltflux_u_in))
      if(allocated(zone5saltflux_v_in))call netcdfrestart_wrt(zone5saltflux_v_in,'zone5saltflux_v_in',lbound(zone5saltflux_v_in),ubound(zone5saltflux_v_in))
      if(allocated(zone5tempflux_glb_in))call netcdfrestart_wrt(zone5tempflux_glb_in,'zone5tempflux_glb_in',lbound(zone5tempflux_glb_in),ubound(zone5tempflux_glb_in))
      if(allocated(zone5tempflux_u_in))call netcdfrestart_wrt(zone5tempflux_u_in,'zone5tempflux_u_in',lbound(zone5tempflux_u_in),ubound(zone5tempflux_u_in))
      if(allocated(zone5tempflux_v_in))call netcdfrestart_wrt(zone5tempflux_v_in,'zone5tempflux_v_in',lbound(zone5tempflux_v_in),ubound(zone5tempflux_v_in))
      if(allocated(zone5waterflux_glb_in))call netcdfrestart_wrt(zone5waterflux_glb_in,'zone5waterflux_glb_in',lbound(zone5waterflux_glb_in),ubound(zone5waterflux_glb_in))
      if(allocated(zone5waterflux_u_in))call netcdfrestart_wrt(zone5waterflux_u_in,'zone5waterflux_u_in',lbound(zone5waterflux_u_in),ubound(zone5waterflux_u_in))
      if(allocated(zone5waterflux_v_in))call netcdfrestart_wrt(zone5waterflux_v_in,'zone5waterflux_v_in',lbound(zone5waterflux_v_in),ubound(zone5waterflux_v_in))
      if(allocated(zone5saltflux_glb_out))call netcdfrestart_wrt(zone5saltflux_glb_out,'zone5saltflux_glb_out',lbound(zone5saltflux_glb_out),ubound(zone5saltflux_glb_out))
      if(allocated(zone5saltflux_u_out))call netcdfrestart_wrt(zone5saltflux_u_out,'zone5saltflux_u_out',lbound(zone5saltflux_u_out),ubound(zone5saltflux_u_out))
      if(allocated(zone5saltflux_v_out))call netcdfrestart_wrt(zone5saltflux_v_out,'zone5saltflux_v_out',lbound(zone5saltflux_v_out),ubound(zone5saltflux_v_out))
      if(allocated(zone5tempflux_glb_out))call netcdfrestart_wrt(zone5tempflux_glb_out,'zone5tempflux_glb_out',lbound(zone5tempflux_glb_out),ubound(zone5tempflux_glb_out))
      if(allocated(zone5tempflux_u_out))call netcdfrestart_wrt(zone5tempflux_u_out,'zone5tempflux_u_out',lbound(zone5tempflux_u_out),ubound(zone5tempflux_u_out))
      if(allocated(zone5tempflux_v_out))call netcdfrestart_wrt(zone5tempflux_v_out,'zone5tempflux_v_out',lbound(zone5tempflux_v_out),ubound(zone5tempflux_v_out))
      if(allocated(zone5waterflux_glb_out))call netcdfrestart_wrt(zone5waterflux_glb_out,'zone5waterflux_glb_out',lbound(zone5waterflux_glb_out),ubound(zone5waterflux_glb_out))
      if(allocated(zone5waterflux_u_out))call netcdfrestart_wrt(zone5waterflux_u_out,'zone5waterflux_u_out',lbound(zone5waterflux_u_out),ubound(zone5waterflux_u_out))
      if(allocated(zone5waterflux_v_out))call netcdfrestart_wrt(zone5waterflux_v_out,'zone5waterflux_v_out',lbound(zone5waterflux_v_out),ubound(zone5waterflux_v_out))
      if(allocated(zone5tempcumul_glb))write(9)zone5tempcumul_glb
      if(allocated(zone5tempcumul_loc))write(9)zone5tempcumul_loc
      if(allocated(zone5tempmasst0))write(9)zone5tempmasst0
      if(allocated(zone5saltcumul_glb))write(9)zone5saltcumul_glb
      if(allocated(zone5saltcumul_loc))write(9)zone5saltcumul_loc
      if(allocated(zone5saltmasst0))write(9)zone5saltmasst0
      if(allocated(zone5watercumul_glb))write(9)zone5watercumul_glb
      if(allocated(zone5watercumul_loc))write(9)zone5watercumul_loc
      if(allocated(zone5watermasst0))write(9)zone5watermasst0
      if(allocated(zone5_mask))call netcdfrestart_wrt(zone5_mask,'zone5_mask',lbound(zone5_mask),ubound(zone5_mask))
      if(allocated(zone5_flux_u_node))call netcdfrestart_wrt(zone5_flux_u_node,'zone5_flux_u_node',lbound(zone5_flux_u_node),ubound(zone5_flux_u_node))
      if(allocated(zone5_flux_v_node))call netcdfrestart_wrt(zone5_flux_v_node,'zone5_flux_v_node',lbound(zone5_flux_v_node),ubound(zone5_flux_v_node))
      write(9)zone6saltflux_w
      write(9)zone6tempflux_w
      write(9)zone6waterflux_w
      write(9)zone6_nlayer
      write(9)zone6_max
      write(9)zone6_u_max
      write(9)zone6_v_max
      write(9)zone6_inv_dz
      write(9)zone6_stretch_dz
      if(allocated(zone6saltflux_glb))call netcdfrestart_wrt(zone6saltflux_glb,'zone6saltflux_glb',lbound(zone6saltflux_glb),ubound(zone6saltflux_glb))
      if(allocated(zone6saltflux_u))call netcdfrestart_wrt(zone6saltflux_u,'zone6saltflux_u',lbound(zone6saltflux_u),ubound(zone6saltflux_u))
      if(allocated(zone6saltflux_v))call netcdfrestart_wrt(zone6saltflux_v,'zone6saltflux_v',lbound(zone6saltflux_v),ubound(zone6saltflux_v))
      if(allocated(zone6tempflux_glb))call netcdfrestart_wrt(zone6tempflux_glb,'zone6tempflux_glb',lbound(zone6tempflux_glb),ubound(zone6tempflux_glb))
      if(allocated(zone6tempflux_u))call netcdfrestart_wrt(zone6tempflux_u,'zone6tempflux_u',lbound(zone6tempflux_u),ubound(zone6tempflux_u))
      if(allocated(zone6tempflux_v))call netcdfrestart_wrt(zone6tempflux_v,'zone6tempflux_v',lbound(zone6tempflux_v),ubound(zone6tempflux_v))
      if(allocated(zone6waterflux_glb))call netcdfrestart_wrt(zone6waterflux_glb,'zone6waterflux_glb',lbound(zone6waterflux_glb),ubound(zone6waterflux_glb))
      if(allocated(zone6waterflux_u))call netcdfrestart_wrt(zone6waterflux_u,'zone6waterflux_u',lbound(zone6waterflux_u),ubound(zone6waterflux_u))
      if(allocated(zone6waterflux_v))call netcdfrestart_wrt(zone6waterflux_v,'zone6waterflux_v',lbound(zone6waterflux_v),ubound(zone6waterflux_v))
      if(allocated(zone6saltflux_glb_in))call netcdfrestart_wrt(zone6saltflux_glb_in,'zone6saltflux_glb_in',lbound(zone6saltflux_glb_in),ubound(zone6saltflux_glb_in))
      if(allocated(zone6saltflux_u_in))call netcdfrestart_wrt(zone6saltflux_u_in,'zone6saltflux_u_in',lbound(zone6saltflux_u_in),ubound(zone6saltflux_u_in))
      if(allocated(zone6saltflux_v_in))call netcdfrestart_wrt(zone6saltflux_v_in,'zone6saltflux_v_in',lbound(zone6saltflux_v_in),ubound(zone6saltflux_v_in))
      if(allocated(zone6tempflux_glb_in))call netcdfrestart_wrt(zone6tempflux_glb_in,'zone6tempflux_glb_in',lbound(zone6tempflux_glb_in),ubound(zone6tempflux_glb_in))
      if(allocated(zone6tempflux_u_in))call netcdfrestart_wrt(zone6tempflux_u_in,'zone6tempflux_u_in',lbound(zone6tempflux_u_in),ubound(zone6tempflux_u_in))
      if(allocated(zone6tempflux_v_in))call netcdfrestart_wrt(zone6tempflux_v_in,'zone6tempflux_v_in',lbound(zone6tempflux_v_in),ubound(zone6tempflux_v_in))
      if(allocated(zone6waterflux_glb_in))call netcdfrestart_wrt(zone6waterflux_glb_in,'zone6waterflux_glb_in',lbound(zone6waterflux_glb_in),ubound(zone6waterflux_glb_in))
      if(allocated(zone6waterflux_u_in))call netcdfrestart_wrt(zone6waterflux_u_in,'zone6waterflux_u_in',lbound(zone6waterflux_u_in),ubound(zone6waterflux_u_in))
      if(allocated(zone6waterflux_v_in))call netcdfrestart_wrt(zone6waterflux_v_in,'zone6waterflux_v_in',lbound(zone6waterflux_v_in),ubound(zone6waterflux_v_in))
      if(allocated(zone6saltflux_glb_out))call netcdfrestart_wrt(zone6saltflux_glb_out,'zone6saltflux_glb_out',lbound(zone6saltflux_glb_out),ubound(zone6saltflux_glb_out))
      if(allocated(zone6saltflux_u_out))call netcdfrestart_wrt(zone6saltflux_u_out,'zone6saltflux_u_out',lbound(zone6saltflux_u_out),ubound(zone6saltflux_u_out))
      if(allocated(zone6saltflux_v_out))call netcdfrestart_wrt(zone6saltflux_v_out,'zone6saltflux_v_out',lbound(zone6saltflux_v_out),ubound(zone6saltflux_v_out))
      if(allocated(zone6tempflux_glb_out))call netcdfrestart_wrt(zone6tempflux_glb_out,'zone6tempflux_glb_out',lbound(zone6tempflux_glb_out),ubound(zone6tempflux_glb_out))
      if(allocated(zone6tempflux_u_out))call netcdfrestart_wrt(zone6tempflux_u_out,'zone6tempflux_u_out',lbound(zone6tempflux_u_out),ubound(zone6tempflux_u_out))
      if(allocated(zone6tempflux_v_out))call netcdfrestart_wrt(zone6tempflux_v_out,'zone6tempflux_v_out',lbound(zone6tempflux_v_out),ubound(zone6tempflux_v_out))
      if(allocated(zone6waterflux_glb_out))call netcdfrestart_wrt(zone6waterflux_glb_out,'zone6waterflux_glb_out',lbound(zone6waterflux_glb_out),ubound(zone6waterflux_glb_out))
      if(allocated(zone6waterflux_u_out))call netcdfrestart_wrt(zone6waterflux_u_out,'zone6waterflux_u_out',lbound(zone6waterflux_u_out),ubound(zone6waterflux_u_out))
      if(allocated(zone6waterflux_v_out))call netcdfrestart_wrt(zone6waterflux_v_out,'zone6waterflux_v_out',lbound(zone6waterflux_v_out),ubound(zone6waterflux_v_out))
      if(allocated(zone6tempcumul_glb))write(9)zone6tempcumul_glb
      if(allocated(zone6tempcumul_loc))write(9)zone6tempcumul_loc
      if(allocated(zone6tempmasst0))write(9)zone6tempmasst0
      if(allocated(zone6saltcumul_glb))write(9)zone6saltcumul_glb
      if(allocated(zone6saltcumul_loc))write(9)zone6saltcumul_loc
      if(allocated(zone6saltmasst0))write(9)zone6saltmasst0
      if(allocated(zone6watercumul_glb))write(9)zone6watercumul_glb
      if(allocated(zone6watercumul_loc))write(9)zone6watercumul_loc
      if(allocated(zone6watermasst0))write(9)zone6watermasst0
      if(allocated(zone6_mask))call netcdfrestart_wrt(zone6_mask,'zone6_mask',lbound(zone6_mask),ubound(zone6_mask))
      if(allocated(zone6_flux_u_node))call netcdfrestart_wrt(zone6_flux_u_node,'zone6_flux_u_node',lbound(zone6_flux_u_node),ubound(zone6_flux_u_node))
      if(allocated(zone6_flux_v_node))call netcdfrestart_wrt(zone6_flux_v_node,'zone6_flux_v_node',lbound(zone6_flux_v_node),ubound(zone6_flux_v_node))
      if(allocated(anyvar3d))call netcdfrestart_wrt(anyvar3d,'anyvar3d',lbound(anyvar3d),ubound(anyvar3d))
      if(allocated(sigma_fric_wu))call netcdfrestart_wrt(sigma_fric_wu,'sigma_fric_wu',lbound(sigma_fric_wu),ubound(sigma_fric_wu))
      if(allocated(sigma_fric_wv))call netcdfrestart_wrt(sigma_fric_wv,'sigma_fric_wv',lbound(sigma_fric_wv),ubound(sigma_fric_wv))
      if(allocated(dsig_t))call netcdfrestart_wrt(dsig_t,'dsig_t',lbound(dsig_t),ubound(dsig_t))
      if(allocated(anyv3dint))call netcdfrestart_wrt(anyv3dint,'anyv3dint',lbound(anyv3dint),ubound(anyv3dint))
      if(allocated(sshobc_w))call netcdfrestart_wrt(sshobc_w,'sshobc_w',lbound(sshobc_w),ubound(sshobc_w))
      if(allocated(velbarobc_u))call netcdfrestart_wrt(velbarobc_u,'velbarobc_u',lbound(velbarobc_u),ubound(velbarobc_u))
      if(allocated(velbarobc_v))call netcdfrestart_wrt(velbarobc_v,'velbarobc_v',lbound(velbarobc_v),ubound(velbarobc_v))
      if(allocated(sshofl_w))call netcdfrestart_wrt(sshofl_w,'sshofl_w',lbound(sshofl_w),ubound(sshofl_w))
      if(allocated(velbarofl_u))call netcdfrestart_wrt(velbarofl_u,'velbarofl_u',lbound(velbarofl_u),ubound(velbarofl_u))
      if(allocated(velbarofl_v))call netcdfrestart_wrt(velbarofl_v,'velbarofl_v',lbound(velbarofl_v),ubound(velbarofl_v))
      if(allocated(tem_delta_t))call netcdfrestart_wrt(tem_delta_t,'tem_delta_t',lbound(tem_delta_t),ubound(tem_delta_t))
      if(allocated(sal_delta_t))call netcdfrestart_wrt(sal_delta_t,'sal_delta_t',lbound(sal_delta_t),ubound(sal_delta_t))
      if(allocated(anyvar2d))call netcdfrestart_wrt(anyvar2d,'anyvar2d',lbound(anyvar2d),ubound(anyvar2d))
      if(allocated(riverdir))write(9)riverdir
      if(allocated(l_river))write(9)l_river
      if(allocated(river_no))write(9)river_no
      if(allocated(rivertrc_inout))write(9)rivertrc_inout
      if(allocated(rivervel_inout))write(9)rivervel_inout
      if(allocated(riverupwdist))write(9)riverupwdist
      if(allocated(nest_len_in))write(9)nest_len_in
      if(allocated(wsed_explicit))write(9)wsed_explicit
      if(allocated(jpm))write(9)jpm
      if(allocated(iriver))call netcdfrestart_wrt(iriver,'iriver',lbound(iriver),ubound(iriver))
      if(allocated(jriver))call netcdfrestart_wrt(jriver,'jriver',lbound(jriver),ubound(jriver))
      if(allocated(rankcoords))call netcdfrestart_wrt(rankcoords,'rankcoords',lbound(rankcoords),ubound(rankcoords))
      if(allocated(spo_i1_u))write(9)spo_i1_u
      if(allocated(spo_i2_u))write(9)spo_i2_u
      if(allocated(spo_j1_u))write(9)spo_j1_u
      if(allocated(spo_j2_u))write(9)spo_j2_u
      if(allocated(spo_i1_v))write(9)spo_i1_v
      if(allocated(spo_i2_v))write(9)spo_i2_v
      if(allocated(spo_j1_v))write(9)spo_j1_v
      if(allocated(spo_j2_v))write(9)spo_j2_v
      if(allocated(spo_i1_t))write(9)spo_i1_t
      if(allocated(spo_i2_t))write(9)spo_i2_t
      if(allocated(spo_j1_t))write(9)spo_j1_t
      if(allocated(spo_j2_t))write(9)spo_j2_t
      if(allocated(obcstreamf))write(9)obcstreamf
      if(allocated(mask_i_w))write(9)mask_i_w
      if(allocated(mask_i_v))write(9)mask_i_v
      if(allocated(mask_i_u))write(9)mask_i_u
      if(allocated(mask_j_w))write(9)mask_j_w
      if(allocated(mask_j_u))write(9)mask_j_u
      if(allocated(mask_j_v))write(9)mask_j_v
      if(allocated(l2ij_out_u))call netcdfrestart_wrt(l2ij_out_u,'l2ij_out_u',lbound(l2ij_out_u),ubound(l2ij_out_u))
      if(allocated(l2ij_out_v))call netcdfrestart_wrt(l2ij_out_v,'l2ij_out_v',lbound(l2ij_out_v),ubound(l2ij_out_v))
      if(allocated(l2ij_out_w))call netcdfrestart_wrt(l2ij_out_w,'l2ij_out_w',lbound(l2ij_out_w),ubound(l2ij_out_w))
      if(allocated(l2ij_in_u))call netcdfrestart_wrt(l2ij_in_u,'l2ij_in_u',lbound(l2ij_in_u),ubound(l2ij_in_u))
      if(allocated(l2ij_in_v))call netcdfrestart_wrt(l2ij_in_v,'l2ij_in_v',lbound(l2ij_in_v),ubound(l2ij_in_v))
      if(allocated(l2ij_in_w))call netcdfrestart_wrt(l2ij_in_w,'l2ij_in_w',lbound(l2ij_in_w),ubound(l2ij_in_w))
      if(allocated(fillmask_t))call netcdfrestart_wrt(fillmask_t,'fillmask_t',lbound(fillmask_t),ubound(fillmask_t))
      if(allocated(dayindex))write(9)dayindex
      if(allocated(varid))write(9)varid
      if(allocated(grh_nb))write(9)grh_nb
      if(allocated(vardim))write(9)vardim
      if(allocated(varstart))write(9)varstart
      if(allocated(varcount))write(9)varcount
      if(allocated(departuredate))write(9)departuredate
      if(allocated(datesim))call netcdfrestart_wrt(datesim,'datesim',lbound(datesim),ubound(datesim))
      if(allocated(dateobc))call netcdfrestart_wrt(dateobc,'dateobc',lbound(dateobc),ubound(dateobc))
      if(allocated(dateairsea))call netcdfrestart_wrt(dateairsea,'dateairsea',lbound(dateairsea),ubound(dateairsea))
      if(allocated(dateriver))call netcdfrestart_wrt(dateriver,'dateriver',lbound(dateriver),ubound(dateriver))
      if(allocated(cgridshift))write(9)cgridshift
      if(allocated(runoff_dt))write(9)runoff_dt
      if(allocated(runoff_w))call netcdfrestart_wrt(runoff_w,'runoff_w',lbound(runoff_w),ubound(runoff_w))
      if(allocated(frqtide))write(9)frqtide
      if(allocated(v0tide))write(9)v0tide
      if(allocated(ti0tide))write(9)ti0tide
      if(allocated(equitide))write(9)equitide
      if(allocated(tideanaweight))write(9)tideanaweight
      if(allocated(airseafile_prvtime))write(9)airseafile_prvtime
      if(allocated(airseafile_nextime))write(9)airseafile_nextime
      if(allocated(ogcm_readtime_next))write(9)ogcm_readtime_next
      if(allocated(ogcm_readtime_prev))write(9)ogcm_readtime_prev
      if(allocated(ogcm_period_prev))write(9)ogcm_period_prev
      if(allocated(ogcm_period_next))write(9)ogcm_period_next
      if(allocated(timeweightobc))write(9)timeweightobc
      if(allocated(dtobc))write(9)dtobc
      if(allocated(dt_abl_max))write(9)dt_abl_max
      if(allocated(zref_z))write(9)zref_z
      if(allocated(checkxyt))write(9)checkxyt
      if(allocated(checkxyf))write(9)checkxyf
      if(allocated(checkanyv3d))write(9)checkanyv3d
      if(allocated(tab1_code))write(9)tab1_code
      if(allocated(tab2_code))write(9)tab2_code
      if(allocated(tab_code_glb))write(9)tab_code_glb
      if(allocated(ub2))write(9)ub2
      if(allocated(lb2))write(9)lb2
      if(allocated(ub3))write(9)ub3
      if(allocated(lb3))write(9)lb3
      if(allocated(ub4))write(9)ub4
      if(allocated(lb4))write(9)lb4
      if(allocated(nutide))write(9)nutide
      if(allocated(nest_len_out))write(9)nest_len_out
      if(allocated(ind))write(9)ind
      if(allocated(airseafile_prvtrec))write(9)airseafile_prvtrec
      if(allocated(airseafile_nextrec))write(9)airseafile_nextrec
      if(allocated(mpi_neighbor_list))write(9)mpi_neighbor_list
      if(allocated(ogcm_rec_next))write(9)ogcm_rec_next
      if(allocated(ogcm_rec_prev))write(9)ogcm_rec_prev
      if(allocated(obcstatus))write(9)obcstatus
      write(9)drifter_out_sampling
      write(9)obc2dtype
      write(9)coef_diss_mangrove
      write(9)expnum
      write(9)cst_c0cub
      write(9)relativewind
      write(9)offset_sshobc
      write(9)biharm_2dfactor
      write(9)checkr0
      write(9)checkr1
      write(9)checkr2
      write(9)checkr3
      write(9)checkr4
      write(9)sponge_l
      write(9)sponge_dx_critic
      write(9)sponge_dx_width
      write(9)dzsurfmin
      if(allocated(drifter_send_order_nord))write(9)drifter_send_order_nord
      if(allocated(drifter_send_order_sud))write(9)drifter_send_order_sud
      if(allocated(drifter_send_order_est))write(9)drifter_send_order_est
      if(allocated(drifter_send_order_ouest))write(9)drifter_send_order_ouest
      if(allocated(drifter_send_order_out))write(9)drifter_send_order_out
      if(allocated(drifter_send_order_nordest))write(9)drifter_send_order_nordest
      if(allocated(drifter_send_order_nordouest))write(9)drifter_send_order_nordouest
      if(allocated(drifter_send_order_sudest))write(9)drifter_send_order_sudest
      if(allocated(drifter_send_order_sudouest))write(9)drifter_send_order_sudouest
      if(allocated(nd_send_canal))write(9)nd_send_canal
      if(allocated(nd_recv_canal))write(9)nd_recv_canal
      if(allocated(zw_abl))write(9)zw_abl
      if(allocated(drifter_send_order_canal))call netcdfrestart_wrt(drifter_send_order_canal,'drifter_send_order_canal',lbound(drifter_send_order_canal),ubound(drifter_send_order_canal))
      write(9)drifter_onoff
      if(allocated(river_obctype))write(9)river_obctype
      if(allocated(river_s))write(9)river_s
      if(allocated(river_tmin))write(9)river_tmin
      if(allocated(river_tmax))write(9)river_tmax
      if(allocated(hriver))write(9)hriver
      if(allocated(friver))write(9)friver
      if(allocated(realriver))write(9)realriver
      if(allocated(riverinfo))write(9)riverinfo
      if(allocated(river_timeref))write(9)river_timeref
      if(allocated(daysim))write(9)daysim
      if(allocated(tdate_output))write(9)tdate_output
      if(allocated(obcinfo))write(9)obcinfo
      if(allocated(offlinedt))write(9)offlinedt
      if(allocated(pss_mean))write(9)pss_mean
      if(allocated(nest_dt_in))write(9)nest_dt_in
      if(allocated(nest_dt_out))write(9)nest_dt_out
      if(allocated(wavedt))write(9)wavedt
      if(allocated(obc_hf_dt))write(9)obc_hf_dt
      if(allocated(glob_dte_lp))write(9)glob_dte_lp
      if(allocated(glob_dte_lp_tmp))write(9)glob_dte_lp_tmp
      if(allocated(alongresolriver))call netcdfrestart_wrt(alongresolriver,'alongresolriver',lbound(alongresolriver),ubound(alongresolriver))
      if(allocated(crossresolriver))call netcdfrestart_wrt(crossresolriver,'crossresolriver',lbound(crossresolriver),ubound(crossresolriver))
      write(9)ww3_varmax
      write(9)ww3_type_grid
      write(9)type_unstructured
      write(9)type_structured
      write(9)vststep
      write(9)zprofile2d3d
      write(9)timestep_type
      write(9)timestep_leapfrog
      write(9)timestep_forwbckw
      write(9)bulk_core
      write(9)bulk_moon
      write(9)bulk_coare
      write(9)bulk_ecume
      write(9)bulk_scheme
      if(allocated(freq_wave))write(9)freq_wave
      write(9)tideana_spinup
      write(9)tideana_delta
      write(9)ogcm_time_lag
      write(9)albedo_val
      write(9)tfilterfb
      write(9)tfilterlf
      write(9)ktide
      write(9)tideforces
      write(9)tideana_yesno
      write(9)ioffline
      write(9)tideanalysis_count
      write(9)ibl1_advbio
      write(9)ibl2_advbio
      write(9)jbl1_advbio
      write(9)jbl2_advbio
      write(9)ogcm_time_shift
      write(9)ieq1
      write(9)ieqimax
      write(9)jeq1
      write(9)jeqjmax
      write(9)ieq1_jeq1
      write(9)ieq1_jeqjmax
      write(9)ieqimax_jeq1
      write(9)ieqimax_jeqjmax
      write(9)dim_varid
      write(9)albedo_constant
      write(9)albedo_apel1987
      write(9)albedo_br1982
      write(9)albedo_br1986
      write(9)albedo_case
      write(9)il_oasis_time
      write(9)flag_steric_effect
      write(9)flag_remove_secondary_bassins
      write(9)one_kind1
      write(9)signe
      write(9)flag_maxbotstress
      write(9)flag_offline_binary
      write(9)flag_wstressbulk
      write(9)flag_write
      write(9)mangrove_scheme
      write(9)flag_meteo_land_plug
      write(9)flag_meteo_land_plug_wind
      write(9)flag_upwind_obc
      write(9)discard_lonlat_periodicity
      write(9)fplan2_grid
      write(9)flag_ts_effectivedensity
      write(9)flag_ts_quicklim
      write(9)flag_z2dv_outputs
      write(9)flag_ogcmname2date
      write(9)flag_meteoname2date
      write(9)ofl_bio
      write(9)drifter_output_files
      write(9)flag_tide3d_analysis
      write(9)flag_bathy_update
      write(9)loop1
      write(9)loop2
      write(9)loop3
      write(9)loopmaxtke
      write(9)loopmaxbio
      write(9)loopmaxts
      write(9)nairsea
      write(9)bi_onoff
      write(9)nc_or_bin_airsea
      write(9)loop_netcdf
      write(9)count_netcdfvar
      write(9)wavefile_prvtrec
      write(9)wavefile_nextrec
      write(9)wave_cpl_nextrec
      write(9)wave_cpl_period_iter
      write(9)wave_cpl_ww3_sdir
      write(9)wave_cpl_ww3_cdir
      write(9)wave_cpl_ww3_hs
      write(9)wave_cpl_ww3_hsw
      write(9)wave_cpl_ww3_foc
      write(9)wave_cpl_ww3_tw
      write(9)wave_cpl_ww3_tawx
      write(9)wave_cpl_ww3_tawy
      write(9)wave_cpl_ww3_twox
      write(9)wave_cpl_ww3_twoy
      write(9)wave_cpl_ww3_uss
      write(9)wave_cpl_ww3_vss
      write(9)wave_cpl_ww3_msk
      write(9)ofl_rec_now
      write(9)ofl_rec_max
      write(9)flag_kz_enhanced
      write(9)flag_net_ir
      write(9)flag_dt_adjust
      write(9)tide_interpolation
      write(9)tide_flagrotation
      write(9)flag_ssr24avr
      write(9)flag_abl
      write(9)flag_abl2
      write(9)flag_sequoia
      write(9)dimssr24prv
      write(9)flag_meteo_average
      write(9)flag_nemoffline
      write(9)trc_id
      write(9)vel_id
      write(9)ssh_id
      write(9)var_num
      write(9)flag_p0m_filter
      write(9)flag_refstate
      write(9)flag_linearfric
      write(9)linear_coef_mangrove
      write(9)flag_1dv
      write(9)flag_merged_levels
      write(9)flag1_smooth_h_mask
      write(9)flag3_smooth_h_mask
      write(9)flag_0status_option
      write(9)flag_rmnegval
      write(9)timemax
      write(9)dirmax
      write(9)freqmax
      write(9)var_misval
      write(9)un_r8
      write(9)heure
      write(9)suma
      write(9)sumb
      write(9)sum_mi
      write(9)grid_area
      write(9)grid_areaglb
      write(9)grid_volumeglb
      write(9)sumarchive
      write(9)lonmin
      write(9)lonmax
      write(9)latmin
      write(9)latmax
      write(9)wavefile_prvtime
      write(9)wavefile_nextime
      write(9)ofl_period_prev
      write(9)ofl_period_now
      write(9)ofl_period_next
      write(9)ofl_nextrec_time
      write(9)ofl_writime
      write(9)ofl_readtime_next
      write(9)ofl_readtime_prev
      write(9)biobc_nextfiletime
      write(9)biobc_prevfiletime
      write(9)stability_index
      write(9)iteration2d_max_r8
      write(9)iteration2d_upbound
      write(9)coef_linearfric
      write(9)gh
      write(9)gm
      write(9)nn
      write(9)tke2overeps
      write(9)tkeovereps
      write(9)gravoverrho
      write(9)sh
      write(9)sm
      write(9)heatfluxbias
      write(9)check0
      write(9)check1
      write(9)check2
      write(9)check3
      write(9)check4
      write(9)check5
      write(9)check6
      write(9)check7
      write(9)check8
      write(9)zero
      write(9)un
      write(9)cdseuil
      write(9)cdb_2dh
      write(9)small1
      write(9)deci
      write(9)decj
      write(9)deck
      write(9)rap1
      write(9)rap2
      write(9)rapi
      write(9)rapj
      write(9)rapk
      write(9)rap
      write(9)const0
      write(9)const1
      write(9)const2
      write(9)const3
      write(9)const4
      write(9)const5
      write(9)const6
      write(9)const7
      write(9)const8
      write(9)const9
      write(9)uv_10
      write(9)karman
      write(9)stefan
      write(9)z2m
      write(9)z10m
      write(9)pss0
      write(9)boltz
      write(9)planck
      write(9)avogadro
      write(9)ce
      write(9)cen
      write(9)ch
      write(9)chn
      write(9)cd
      write(9)cdn
      write(9)z_2
      write(9)z_10
      write(9)q_0
      write(9)r_0
      write(9)pvs_0
      write(9)psih_10
      write(9)psih_2
      write(9)phim
      write(9)zl_10
      write(9)zl_2
      write(9)falpha
      write(9)fbeta
      write(9)ro
      write(9)cp_air
      write(9)lv
      write(9)psim_10
      write(9)psim_2
      write(9)dte_lp
      write(9)inv_dte_lp
      write(9)inv_dti_lp
      write(9)inv_dti_fw
      write(9)inv_dti_fw_p2
      write(9)dte_fw
      write(9)dti_lpbef
      write(9)dti_lp
      write(9)dti_lpmax
      write(9)dti_lpsub
      write(9)dti_fwsub
      write(9)dti_fwsubio
      write(9)dti_fw
      write(9)dti_bef
      write(9)dti_now
      write(9)dtiratio
      write(9)dt_drf
      write(9)fbtfiltercoef
      write(9)assel0
      write(9)assel1
      write(9)assel2
      write(9)assel3
      write(9)wetdry_cst1
      write(9)wetdry_cst2
      write(9)wetdry_cst3
      write(9)h_inf
      write(9)h_inf_obc
      write(9)h_sup
      write(9)dist
      write(9)dist0
      write(9)dist1
      write(9)dist2
      write(9)dist3
      write(9)t0surf
      write(9)s0surf
      write(9)deg2rad
      write(9)rad2deg
      write(9)zmin
      write(9)zmax
      write(9)coa
      write(9)cob
      write(9)coc
      write(9)sol1
      write(9)sol2
      write(9)discri
      write(9)profm1
      write(9)profp1
      write(9)hmax
      write(9)hstepmin
      write(9)hstepmax
      write(9)grav
      write(9)cfl_sshmax
      write(9)cfl_hsshmax
      write(9)cfl_umax
      write(9)cfl_reduce
      write(9)relax_es
      write(9)relax_ext
      write(9)relax_int
      write(9)relax_ts
      write(9)relax_bpc
      write(9)momentum_input_depth
      write(9)z0s
      write(9)z1
      write(9)z2
      write(9)z3
      write(9)z4
      write(9)z5
      write(9)y0
      write(9)y2
      write(9)y3
      write(9)y4
      write(9)y5
      write(9)x0
      write(9)x1
      write(9)x2
      write(9)x3
      write(9)x4
      write(9)x5
      write(9)x6
      write(9)x7
      write(9)x8
      write(9)x9
      write(9)x10
      write(9)x11
      write(9)x12
      write(9)x13
      write(9)x14
      write(9)x20
      write(9)x21
      write(9)x22
      write(9)x33
      write(9)x44
      write(9)area
      write(9)tem_validmin
      write(9)tem_validmax
      write(9)sal_validmin
      write(9)sal_validmax
      write(9)xmd
      write(9)xmv
      write(9)xrd
      write(9)xrv
      write(9)xcpd
      write(9)xcpv
      write(9)xcl
      write(9)celsius2kelvin
      write(9)xlvtt
      write(9)xestt
      write(9)xgamw
      write(9)xbetaw
      write(9)xalpw
      write(9)zrvsrdm1
      write(9)z0_u
      write(9)qsat_sea_z
      write(9)airdensity
      write(9)sst_kelvin
      write(9)prs_atm_z
      write(9)tem_atm_z
      write(9)exner_atm_z
      write(9)delta_u
      write(9)delta_t
      write(9)delta_q
      write(9)delta_u_n
      write(9)exner_sea_z
      write(9)psifunctt
      write(9)psifunctu
      write(9)z0_q
      write(9)z0_t
      write(9)sst1000hpa_kelvin
      write(9)visa
      write(9)charnock
      write(9)qsat_atm_z
      write(9)ustar_bef
      write(9)qstar_bef
      write(9)tetastar_bef
      write(9)rayonterre
      write(9)northpole_lon
      write(9)northpole_lat
      write(9)southpole_lon
      write(9)southpole_lat
      write(9)phi0
      write(9)longi
      write(9)latit
      write(9)longi0
      write(9)latit0
      write(9)longi1
      write(9)latit1
      write(9)angle0
      write(9)alp_t
      write(9)alp_s
      write(9)pi
      write(9)rho
      write(9)inv_rho
      write(9)rhoair
      write(9)valmax
      write(9)vis
      write(9)tkee1
      write(9)tkee2
      write(9)tkee3
      write(9)tkeg2
      write(9)tkeg3
      write(9)tkeg4
      write(9)tkeg5
      write(9)tkeg6
      write(9)tkeb1
      write(9)tkeb2
      write(9)ctke1
      write(9)ctke2
      write(9)t0
      write(9)s0
      write(9)cp
      write(9)light_kpar1
      write(9)light_att1
      write(9)light_att2
      write(9)light_rat1
      write(9)light_rat2
      write(9)light_att2_val1
      write(9)light_att2_h1
      write(9)light_att2_val2
      write(9)light_att2_h2
      write(9)t0_base
      write(9)s0_base
      write(9)rho_base
      write(9)alp_t_base
      write(9)alp_s_base
      write(9)tfb0
      write(9)meteo_lonmin
      write(9)meteo_latmin
      write(9)meteo_lonmax
      write(9)meteo_latmax
      write(9)meteo_resol
      write(9)meteo_resol_u
      write(9)meteo_resol_v
      write(9)meteo_lonstr
      write(9)meteo_lonend
      write(9)meteo_londlt
      write(9)meteo_latstr
      write(9)meteo_latend
      write(9)meteo_latdlt
      write(9)var_lonmin
      write(9)var_latmin
      write(9)var_lonmax
      write(9)var_latmax
      write(9)ww3_lonmin
      write(9)ww3_latmin
      write(9)ww3_lonmax
      write(9)ww3_latmax
      write(9)ww3_dlon
      write(9)ww3_dlat
      write(9)tide_lonmin
      write(9)tide_latmin
      write(9)tide_lonmax
      write(9)tide_latmax
      write(9)tide_dlon
      write(9)tide_dlat
      write(9)dxb
      write(9)dyb
      write(9)dxa
      write(9)dya
      write(9)hmin
      write(9)h1d
      write(9)dlon
      write(9)dlat
      write(9)epsi
      write(9)lagrange_ssh
      write(9)diffu
      write(9)difnorm
      write(9)lup
      write(9)ldown
      write(9)zup
      write(9)zdown
      write(9)rbase
      write(9)small
      write(9)small3
      write(9)hgesig
      write(9)pgesig
      write(9)windfactor
      write(9)xdtk_out
      write(9)tfond
      write(9)sfond
      write(9)rfond
      write(9)c1streamf
      write(9)c2streamf
      write(9)rampe
      write(9)rampe_wind
      write(9)y1
      write(9)obc_hf_reset
      write(9)rho_0d
      write(9)tem_0d
      write(9)sal_0d
      write(9)rho_tmp
      write(9)cst_adv_hor
      write(9)cst_adv_ver
      write(9)cst_adv_vel
      write(9)ssh_avr_nest_out
      write(9)graph_nextime
      write(9)tidenodal_prev_rdv
      write(9)tidenodal_next_rdv
      write(9)tideana_modulo
      write(9)tideana_nextime
      write(9)cellboxfactor1
      write(9)cellboxfactor2
      write(9)rap_wave
      write(9)ihmax
      write(9)jhmax
      write(9)ihmin
      write(9)jhmin
      write(9)convect_yn
      write(9)nbvstepmin
      write(9)bio_relax_size
      write(9)unit_r4
      write(9)x0_r4
      write(9)x1_r4
      write(9)x2_r4
      write(9)x3_r4
      write(9)x4_r4
      write(9)time_r4
      write(9)discharge
      write(9)filval
      write(9)var_validmin
      write(9)var_validmax
      write(9)vdw_loc
      write(9)vup_loc
      write(9)zdw_loc
      write(9)var_scalefactor
      write(9)inv_scalefactor
      write(9)var_addoffset
      write(9)zup_loc
      write(9)hrmax
      write(9)relax_bio
      write(9)kmol_m
      write(9)kmol_h
      write(9)kmol_s
      write(9)upw_hrange1
      write(9)upw_hrange2
      write(9)inv_ekman_depth
      write(9)constant_km
      write(9)constant_kh
      write(9)z0b
      write(9)z0b_land
      write(9)z0b_rivers
      write(9)zlevel_land
      write(9)invloopmaxts
      write(9)invloopmaxu
      write(9)invloopmaxv
      write(9)sqrtgrav
      write(9)invgrav
      write(9)spinup_forcing
      write(9)relax_lwf
      write(9)dz_vertical_incr_fact
      write(9)ratio_negdif_ver
      write(9)ratio_negdif_hor
      write(9)ratio_bionegdif
      write(9)vqs_cst1
      write(9)vqs_cst2
      write(9)vqs_cst3
      write(9)ema_mu
      write(9)dz_over_z0_min
      write(9)coastal_viscosity
      write(9)quick_coef
      write(9)i
      write(9)j
      write(9)k
      write(9)flag3d
      write(9)lrec
      write(9)iteration3d
      write(9)compt1
      write(9)compt2
      write(9)compt3
      write(9)compt4
      write(9)kount0
      write(9)kount1
      write(9)kount2
      write(9)kount3
      write(9)kount4
      write(9)kount5
      write(9)kount6
      write(9)kount7
      write(9)kount8
      write(9)kount9
      write(9)kountmod
      write(9)kountrdv1
      write(9)kountrdv2
      write(9)kountrdv3
      write(9)kountrdv4
      write(9)substep_advbio
      write(9)subcycle_exchange
      write(9)subcycle_onoff
      write(9)subcycle_synchro
      write(9)subcycle_modulo
      write(9)dt_drf_over_dti_fw
      write(9)iteration3d_restart
      write(9)filvalshort
      write(9)quick_filter_points
      write(9)status
      write(9)forcedstatus
      write(9)decision
      write(9)ncid1
      write(9)ncid2
      write(9)dim_x_id
      write(9)dim_y_id
      write(9)dim_z_id
      write(9)dim_t_id
      write(9)dim_b_id
      write(9)max_x
      write(9)max_y
      write(9)max_z
      write(9)max_time_counter
      write(9)max_meteo_time_counter
      write(9)var_id
      write(9)var_nftype
      write(9)meteo_imax
      write(9)meteo_jmax
      write(9)meteo_kmax
      write(9)meteozoom_istr
      write(9)meteozoom_iend
      write(9)meteozoom_jstr
      write(9)meteozoom_jend
      write(9)meteofull_imax
      write(9)meteofull_jmax
      write(9)tide_imax
      write(9)tide_jmax
      write(9)tide_kmax
      write(9)tidezoom_istr
      write(9)tidezoom_iend
      write(9)tidezoom_jstr
      write(9)tidezoom_jend
      write(9)tidezoom_istr_t
      write(9)tidezoom_iend_t
      write(9)tidezoom_jstr_t
      write(9)tidezoom_jend_t
      write(9)tidezoom_istr_u
      write(9)tidezoom_iend_u
      write(9)tidezoom_jstr_u
      write(9)tidezoom_jend_u
      write(9)tidezoom_istr_v
      write(9)tidezoom_iend_v
      write(9)tidezoom_jstr_v
      write(9)tidezoom_jend_v
      write(9)tidefull_imax
      write(9)tidefull_jmax
      write(9)ww3_imax
      write(9)ww3_jmax
      write(9)ww3_kmax
      write(9)ww3_fmax
      write(9)ww3zoom_istr
      write(9)ww3zoom_iend
      write(9)ww3zoom_jstr
      write(9)ww3zoom_jend
      write(9)ww3full_imax
      write(9)ww3full_jmax
      write(9)jour
      write(9)kstop
      write(9)istr
      write(9)jstr
      write(9)kstr
      write(9)tstr
      write(9)bstr
      write(9)iend
      write(9)jend
      write(9)kend
      write(9)tend
      write(9)bend
      write(9)dimend
      write(9)kpvwave
      write(9)give_chanel9
      write(9)i0
      write(9)j0
      write(9)iteration2d
      write(9)iteration2d_begin
      write(9)iteration2d_max_now
      write(9)iteration2d_max_bef
      write(9)i2dh
      write(9)l1
      write(9)l1_sca
      write(9)l1_vec
      write(9)l2
      write(9)l2_sca
      write(9)l2_vec
      write(9)l3
      write(9)l3_sca
      write(9)l3_vec
      write(9)len1
      write(9)nc
      write(9)nc1
      write(9)itimets
      write(9)itimebio
      write(9)iadvec_ts_hor
      write(9)iadvec_ts_hor_upwind
      write(9)iadvec_ts_hor_quickest
      write(9)iadvec_ts_hor_quickest2
      write(9)iadvec_ts_ver
      write(9)iadvec_ts_ver_quickest
      write(9)iadvec_ts_ver_c2
      write(9)iadvec_ts_ver_quickest2
      write(9)iturbulence
      write(9)istreamf
      write(9)itime
      write(9)kts
      write(9)kuv
      write(9)ko
      write(9)jm1
      write(9)jm2
      write(9)jp1
      write(9)jp2
      write(9)im1
      write(9)im2
      write(9)ip1
      write(9)ip2
      write(9)iwind
      write(9)ip
      write(9)im
      write(9)jp
      write(9)jm
      write(9)kp
      write(9)km
      write(9)nbcanal
      write(9)gridtype1
      write(9)gridtype2
      write(9)point1
      write(9)point2
      write(9)sender
      write(9)receiver
      write(9)ipnoc
      write(9)i1
      write(9)i2
      write(9)i3
      write(9)i4
      write(9)i5
      write(9)i6
      write(9)i7
      write(9)i8
      write(9)i9
      write(9)i10
      write(9)i11
      write(9)j1
      write(9)j2
      write(9)j3
      write(9)j4
      write(9)j5
      write(9)j6
      write(9)j7
      write(9)j8
      write(9)j9
      write(9)j10
      write(9)j11
      write(9)k0
      write(9)k1
      write(9)k2
      write(9)k3
      write(9)k4
      write(9)k5
      write(9)k6
      write(9)k7
      write(9)k8
      write(9)k9
      write(9)kr
      write(9)kp1
      write(9)kp2
      write(9)km1
      write(9)km2
      write(9)key
      write(9)flag
      write(9)nbomax
      write(9)nbobuffermax
      write(9)nbobuffermax_c
      write(9)flag_stop
      write(9)itest
      write(9)itest1
      write(9)itest2
      write(9)itest3
      write(9)ioption
      write(9)nsmooth
      write(9)nriver
      write(9)ian0
      write(9)imois0
      write(9)ijour0
      write(9)iheure0
      write(9)iminute0
      write(9)iseconde0
      write(9)iref
      write(9)jref
      write(9)lname1
      write(9)lname2
      write(9)lname3
      write(9)lname4
      write(9)kmode
      write(9)kmodemax
      write(9)fgrid_or_wgrid
      write(9)fgrid_case
      write(9)wgrid_case
      write(9)typegrid
      write(9)typegrid_monopole
      write(9)typegrid_file
      write(9)typegrid_bipole
      write(9)istar
      write(9)jstar
      write(9)istop
      write(9)jstop
      write(9)isend
      write(9)jsend
      write(9)irecv
      write(9)jrecv
      write(9)izoomin
      write(9)izoomax
      write(9)jzoomin
      write(9)jzoomax
      write(9)iairsea
      write(9)ialbedo
      write(9)airseaoption
      write(9)iobc_f
      write(9)iobc_wv
      write(9)iobc_ogcm
      write(9)iobc_lr
      write(9)obc_option
      write(9)iobc_demo_wv
      write(9)iarchive
      write(9)imodeltrc
      write(9)imodelbio
      write(9)multiple
      write(9)kvarmax
      write(9)igesig
      write(9)isigfile
      write(9)ksecu
      write(9)kmaxtide
      write(9)kmaxtidep1
      write(9)kminserie
      write(9)nzctdmax
      write(9)nsctdmax
      write(9)ktctdmin
      write(9)ktctdmax
      write(9)kdtk_out
      write(9)k10
      write(9)k11
      write(9)nbinco
      write(9)nbequa
      write(9)nbsparse
      write(9)kland
      write(9)tidestep2
      write(9)tidestep3
      write(9)ncfilm_max
      write(9)in_out_tide
      write(9)kmax_dof
      write(9)k_in
      write(9)k_out
      write(9)used_unused_dom
      write(9)mergebathy_sponge
      write(9)rankhmax
      write(9)rankhmin
      write(9)id_tem
      write(9)id_tem2
      write(9)id_dtem
      write(9)id_sal
      write(9)id_sal2
      write(9)id_rhp
      write(9)id_rhop
      write(9)id_rhom
      write(9)id_rhf
      write(9)id_rhc
      write(9)id_rhpa
      write(9)id_rhpb
      write(9)id_z
      write(9)id_prs
      write(9)id_now
      write(9)id_aft
      write(9)id_eost
      write(9)id_eoss
      write(9)id_ssh
      write(9)id_breaker
      write(9)id_dz
      write(9)id_zt
      write(9)id_tdiv
      write(9)id_sdiv
      write(9)id_udiv
      write(9)id_vdiv
      write(9)id_bdiv
      write(9)id_bnegdif
      write(9)id_biobefo
      write(9)id_bioaftr
      write(9)id_u_now
      write(9)id_v_now
      write(9)id_u_rot
      write(9)id_v_rot
      write(9)id_ffreq
      write(9)id_coriolis1
      write(9)id_coriolis2
      write(9)id_flx
      write(9)id_fly
      write(9)id_uflx
      write(9)id_ufly
      write(9)id_vflx
      write(9)id_vfly
      write(9)id_tcn
      write(9)id_scn
      write(9)id_gradssh
      write(9)id_ofactort
      write(9)id_ofactors
      write(9)id_hybcoefu
      write(9)id_hybcoefv
      write(9)id_dxdydz
      write(9)id_wb
      write(9)id_prod
      write(9)id_buoy
      write(9)iDtOvRhCp
      write(9)id_ncu
      write(9)id_ncv
      write(9)id_webio
      write(9)id_varbef2
      write(9)id_kh_over_dz
      write(9)id_bihar_lim
      write(9)id_veltot
      write(9)id_velexp
      write(9)id_velimp
      write(9)id_wdrifter
      write(9)kreadgroup
      write(9)kreadgroupmax
      write(9)kread1
      write(9)kread2
      write(9)modulo_biotimestep
      write(9)obctime_bef2
      write(9)obctime_bef
      write(9)obctime_aft
      write(9)obctime_aft2
      write(9)obctime_order
      write(9)sch_imp_ts_u_loc
      write(9)sch_imp_ts_v_loc
      write(9)sch_imp_ts_u_glb
      write(9)sch_imp_ts_v_glb
      write(9)sch_imp_tke_u_loc
      write(9)sch_imp_tke_v_loc
      write(9)sch_imp_tke_u_glb
      write(9)sch_imp_tke_v_glb
      write(9)looplimit_hor
      write(9)ihybsig
      write(9)nhybsig
      write(9)tke_surf
      write(9)grh_out_mi
      write(9)nest_onoff_in
      write(9)nest_onoff_demo
      write(9)nest_full_in
      write(9)nest_full_out
      write(9)nest_onoff_out
      write(9)ncmin_airsea
      write(9)ncmin_river
      write(9)eos_author
      write(9)eos_comprs
      write(9)eos_linear
      write(9)obcfreeorfix
      write(9)iwve
      write(9)wave_obc_type
      write(9)dataperwavefile
      write(9)sp_or_db
      write(9)kbu
      write(9)kbumax
      write(9)kbumax_glb
      write(9)n_element
      write(9)relaxtype_ts
      write(9)obctype_ts
      write(9)obctype_p
      write(9)restart_file_y_or_n
      write(9)ncid
      write(9)x_xhl_dim
      write(9)y_xhl_dim
      write(9)z_xhl_dim
      write(9)x_yhl_dim
      write(9)y_yhl_dim
      write(9)z_yhl_dim
      write(9)x_zhl_dim
      write(9)y_zhl_dim
      write(9)z_zhl_dim
      write(9)x_zl_dim
      write(9)y_zl_dim
      write(9)z_zl_dim
      write(9)time_dim
      write(9)i_t_dim
      write(9)j_t_dim
      write(9)k_t_dim
      write(9)i_w_dim
      write(9)j_w_dim
      write(9)k_w_dim
      write(9)i_u_dim
      write(9)j_u_dim
      write(9)k_u_dim
      write(9)i_v_dim
      write(9)j_v_dim
      write(9)k_v_dim
      write(9)i_f_dim
      write(9)j_f_dim
      write(9)k_f_dim
      write(9)dayindex_size
      write(9)tide_year_min
      write(9)airsea_year_min
      write(9)tide_year_max
      write(9)airsea_year_max
      write(9)wave_year_min
      write(9)irelaxsst
      write(9)removetide
      write(9)ofl_rotation
      write(9)ofl_rhp
      write(9)ofl_tke
      write(9)ofl_surflux
      write(9)ale_selected
      write(9)flag_asselin
      write(9)freq
      write(9)dirw
      write(9)dirw_beg
      write(9)dirw_end
      write(9)freq_beg
      write(9)freq_end
      write(9)year_now
      write(9)month_now
      write(9)day_now
      write(9)hour_now
      write(9)minute_now
      write(9)second_now
      write(9)nd_send_est
      write(9)nd_send_ouest
      write(9)nd_send_nord
      write(9)nd_send_sud
      write(9)nd_send_out
      write(9)nd_recv_est
      write(9)nd_recv_ouest
      write(9)nd_recv_nord
      write(9)nd_recv_sud
      write(9)nd_send_sudouest
      write(9)nd_send_sudest
      write(9)nd_send_nordouest
      write(9)nd_send_nordest
      write(9)nd_recv_sudouest
      write(9)nd_recv_sudest
      write(9)nd_recv_nordouest
      write(9)nd_recv_nordest
      write(9)initial_main_status
      write(9)offline_init_status
      write(9)meteo_sealand_mask
      write(9)grid_i0
      write(9)grid_j0
      write(9)ifb
      write(9)dbefore
      write(9)dnow
      write(9)dafter
      write(9)dim_airsea
      write(9)rhp_zavr_xy
      write(9)ssr_id
      write(9)ir_id
      write(9)rain_id
      write(9)t2m_id
      write(9)t0m_id
      write(9)abl_id
      write(9)dp2m_id
      write(9)u10m_id
      write(9)v10m_id
      write(9)u100m_id
      write(9)v100m_id
      write(9)p0m_id
      write(9)ustrs_id
      write(9)vstrs_id
      write(9)slhf_id
      write(9)netir_id
      write(9)sshf_id
      write(9)t_flux_cumul
      write(9)s_flux_cumul
      write(9)cumuldeltaflux
      write(9)som0
      write(9)som2
      write(9)ssh_reservoir
      write(9)idealflux
      write(9)sum0
      write(9)sum1
      write(9)sum2
      write(9)sum3
      write(9)sum4
      write(9)sum5
      write(9)sum6
      write(9)sum7
      write(9)sum8
      write(9)sum9
      write(9)sum10
      write(9)sum11
      write(9)sum12
      write(9)sum13
      write(9)time0
      write(9)time1
      write(9)time2
      write(9)small2
      write(9)x1_r8
      write(9)x2_r8
      write(9)x3_r8
      write(9)x4_r8
      write(9)sum0glb
      write(9)sum1glb
      write(9)sum2glb
      write(9)sum3glb
      write(9)sum4glb
      write(9)sum5glb
      write(9)sum6glb
      write(9)emin
      write(9)epsmin
      write(9)elapsedtime_out
      write(9)elapsedtime_rst
      write(9)elapsedtime_now
      write(9)elapsedtime_last_writing
      write(9)elapsedtime_bef
      write(9)elapsedtime_aft
      write(9)cpu_seconds
      write(9)alpha
      write(9)eos_pgfzref
      write(9)eos_tkezref
      write(9)filvalr8
      write(9)dti_fw_i4
      write(9)elapsedtime_now_i4
      write(9)elapsedtime_now_i8
      write(9)elapsedtime_now_r16
! insere partie_ecrit_1 fin. SURTOUT NE PAS EFFACER CETTE LIGNE

      close(9)
      if(par%rank==0)write(6,*)'ecriture chanel9 fin.'
      call s_cpu('chanel9_writting_aft',1)

      return
      endif  !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww>

      if(debug_==0)stop 'mauvais argument passe dans dyn_restart'
      end subroutine dyn_restart

!.......................................................................

      subroutine dyn_restart_ask(period_,single_or_modulo_loc)    !27-05
      use module_principal
      use module_parallele                                        !#MPI
      implicit none
      double precision period_
      character txt6_loc*6,single_or_modulo_loc*6


      ksecu=0
                               x1=period_*86400.
      if(restartfileunits=='h')x1=period_*3600. !05-10-20
      if(restartfileunits=='m')x1=period_*60.
      if(restartfileunits=='s')x1=period_     !10-04-17
      if(single_or_modulo_loc=='modulo') then !mmmmmmmmmmmmm>
        if(int(elapsedtime_now/x1)/=int(elapsedtime_bef/x1))ksecu=1
      endif                                   !mmmmmmmmmmmmm>
      if(single_or_modulo_loc=='single') then !sssssssssssss>
        if(elapsedtime_now>x1.and.elapsedtime_bef<=x1)ksecu=1
      endif                                   !sssssssssssss>

      call mpi_bcast(ksecu,1,mpi_integer,0,par%comm2d,ierr) !#mpi !21-05
      if(ksecu==0) return

      give_chanel9=1

      return

      call elapsedtime2date(elapsedtime_now,i5,i6,i7,i3,i2,i1)     !01-0
!.....écrire année année mois dans TEXTE90:
      i0=i7+100*i6+10000*i5
      write(texte30(1:8),'(i8)')i0

!.....écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte30(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de séparation, "_"
!     écrire ":" dans TEXTE90:
      write(texte30(9:9),'(a1)')'_'

      if(period_==1.) txt6_loc='_j0p1_'
      if(period_==7.) txt6_loc='_j0p7_'
      if(period_==8.) txt6_loc='_j0p8_'
      if(period_==15.)txt6_loc='j0p15_'
!     texte90=txt6_loc//dom_c//'_'//texte30(1:15)
      write(texte90,'(a,i0,a,a)')txt6_loc,par%rank,'_'//texte30(1:15) !17-01-15
      dynrestartfilename=trim(restartdir_out1)//'chanel9'//trim(texte90) !18-07-20

      end subroutine dyn_restart_ask

!.......................................................................

      subroutine dyn_restart_check                                !27-05
      use module_principal
      use module_parallele                                        !#MPI
      implicit none
      character txt_spy_loc*60

      txt_spy_loc='dyn_restart_check'     !27-05-11

! Ne pas ecrire le fichier si solution bugguee:
      sum1=0. ;  sum0=small1
      do j=1,jmax
      do i=1,imax
      sum0=sum0+dxdy_t(i,j)
      sum1=sum1+dxdy_t(i,j)*ssh_w(i,j,2) !28-05-11
      enddo
      enddo
      sum1=sum1/sum0


      sum2=0. ;  sum3=0.; sum0=small1
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       sum0=sum0+mask_t(i,j,k)*dz_t(i,j,k,2)*dxdy_t(i,j)
       sum2=sum2+mask_t(i,j,k)*dz_t(i,j,k,2)*dxdy_t(i,j)*tem_t(i,j,k,2)
       sum3=sum3+mask_t(i,j,k)*dz_t(i,j,k,2)*dxdy_t(i,j)*sal_t(i,j,k,2)
      enddo
      enddo
      enddo
      sum2=sum2/sum0 ;  sum3=sum3/sum0

      sum4=0. ; sum0=small1
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
       sum0=sum0+mask_u(i,j,k)*dz_u(i,j,k,2)*dxdy_u(i,j)
       sum4=sum4+mask_u(i,j,k)*dz_u(i,j,k,2)*dxdy_u(i,j)*vel_u(i,j,k,2)
      enddo
      enddo
      enddo
      sum4=sum4/sum0

      sum5=0. ; sum0=small1
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
       sum0=sum0+mask_v(i,j,k)*dz_v(i,j,k,2)*dxdy_v(i,j)
       sum5=sum5+mask_v(i,j,k)*dz_v(i,j,k,2)*dxdy_v(i,j)*vel_v(i,j,k,2)
      enddo
      enddo
      enddo
      sum5=sum5/sum0

!     open(unit=3     &
!         ,file='restart_output/checkrestartfile'//dom_c//'.ascii')
      write(texte60,'(a,i0)')trim(restartdir_out1)//'checkrestartfile',par%rank !18-07-20


      ksecu=0
      open(unit=3,file=trim(texte60)) !17-01-15
       write(3,*)sum1
       write(3,*)sum2
       write(3,*)sum3
       write(3,*)sum4
       write(3,*)sum5
      close(3)
      open(unit=3,file=texte60) !17-01-15
      do k=1,5
       read(3,'(a)')texte90
       if(index(texte90,'N')/=0.or.index(texte90,'n')/=0) then !%%%%%>
        ksecu=1
       endif                                                   !%%%%%>
      enddo
      close(3)
      call mpi_allreduce(ksecu,k10,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k10/=0) then !  m[°v°]m  > !08-06-17
        call graph_out
        stop ' Run stopped by dyn_restart_check after NaN detected'
      endif           !  m[°v°]m  >

!     call barriere(iteration3d,1,txt_spy_loc)  !27-05-11
      call barriere(2216,20130309,txt_spy_loc)  !27-05-11

      end subroutine dyn_restart_check
!.......................................................................
      subroutine dyn_restart_bounds_w
      use module_principal
      use module_parallele !#MPI
      implicit none

      k=0   ; k=index(dynrestartfilename,trim(restartdir_out1)) ; k2=len(trim(restartdir_out1))
      if(k==0) then
          k=index(dynrestartfilename,trim(restartdir_out2))
          k2=len(trim(restartdir_out2))
      endif
      write(texte60,'(a,i0)') &
      dynrestartfilename(1:k+k2-1)//'dynrestartbounds',par%rank

      open(unit=3,file=trim(texte60))
! jalon1 dyn_restart_check ne pas effacer
      write(3,*)'ioffline_prv'
      write(3,*)'zone1bioflux_w'
      if(allocated(zone1bioflux_glb))write(3,*)'zone1bioflux_glb',1+ubound(zone1bioflux_glb)-lbound(zone1bioflux_glb)
      if(allocated(zone1bioflux_u))write(3,*)'zone1bioflux_u',1+ubound(zone1bioflux_u)-lbound(zone1bioflux_u)
      if(allocated(zone1bioflux_v))write(3,*)'zone1bioflux_v',1+ubound(zone1bioflux_v)-lbound(zone1bioflux_v)
      if(allocated(zone1tendancebio_glb))write(3,*)'zone1tendancebio_glb',1+ubound(zone1tendancebio_glb)-lbound(zone1tendancebio_glb)
      if(allocated(zone1botsurfbio_glb))write(3,*)'zone1botsurfbio_glb',1+ubound(zone1botsurfbio_glb)-lbound(zone1botsurfbio_glb)
      if(allocated(zone1bioflux_glb_in))write(3,*)'zone1bioflux_glb_in',1+ubound(zone1bioflux_glb_in)-lbound(zone1bioflux_glb_in)
      if(allocated(zone1bioflux_u_in))write(3,*)'zone1bioflux_u_in',1+ubound(zone1bioflux_u_in)-lbound(zone1bioflux_u_in)
      if(allocated(zone1bioflux_v_in))write(3,*)'zone1bioflux_v_in',1+ubound(zone1bioflux_v_in)-lbound(zone1bioflux_v_in)
      if(allocated(zone1bioflux_glb_out))write(3,*)'zone1bioflux_glb_out',1+ubound(zone1bioflux_glb_out)-lbound(zone1bioflux_glb_out)
      if(allocated(zone1bioflux_u_out))write(3,*)'zone1bioflux_u_out',1+ubound(zone1bioflux_u_out)-lbound(zone1bioflux_u_out)
      if(allocated(zone1bioflux_v_out))write(3,*)'zone1bioflux_v_out',1+ubound(zone1bioflux_v_out)-lbound(zone1bioflux_v_out)
      if(allocated(zone1biomasst0))write(3,*)'zone1biomasst0',1+ubound(zone1biomasst0)-lbound(zone1biomasst0)
      if(allocated(analysis3dmatrix))write(3,*)'analysis3dmatrix',1+ubound(analysis3dmatrix)-lbound(analysis3dmatrix)
      if(allocated(ema1_s_i))write(3,*)'ema1_s_i',1+ubound(ema1_s_i)-lbound(ema1_s_i)
      if(allocated(ema2_s_i))write(3,*)'ema2_s_i',1+ubound(ema2_s_i)-lbound(ema2_s_i)
      if(allocated(ema1_s_j))write(3,*)'ema1_s_j',1+ubound(ema1_s_j)-lbound(ema1_s_j)
      if(allocated(ema2_s_j))write(3,*)'ema2_s_j',1+ubound(ema2_s_j)-lbound(ema2_s_j)
      if(allocated(c_wave_mode))write(3,*)'c_wave_mode',1+ubound(c_wave_mode)-lbound(c_wave_mode)
      if(allocated(pcoefmode_t))write(3,*)'pcoefmode_t',1+ubound(pcoefmode_t)-lbound(pcoefmode_t)
      if(allocated(ucoefmode_t))write(3,*)'ucoefmode_t',1+ubound(ucoefmode_t)-lbound(ucoefmode_t)
      if(allocated(vcoefmode_t))write(3,*)'vcoefmode_t',1+ubound(vcoefmode_t)-lbound(vcoefmode_t)
      if(allocated(rhrefmode_t))write(3,*)'rhrefmode_t',1+ubound(rhrefmode_t)-lbound(rhrefmode_t)
      if(allocated(ema1_q_i))write(3,*)'ema1_q_i',1+ubound(ema1_q_i)-lbound(ema1_q_i)
      if(allocated(ema2_q_i))write(3,*)'ema2_q_i',1+ubound(ema2_q_i)-lbound(ema2_q_i)
      if(allocated(ema1_q_j))write(3,*)'ema1_q_j',1+ubound(ema1_q_j)-lbound(ema1_q_j)
      if(allocated(ema2_q_j))write(3,*)'ema2_q_j',1+ubound(ema2_q_j)-lbound(ema2_q_j)
      if(allocated(uv_wmode_t))write(3,*)'uv_wmode_t',1+ubound(uv_wmode_t)-lbound(uv_wmode_t)
      if(allocated(analysis_p3d_t))write(3,*)'analysis_p3d_t',1+ubound(analysis_p3d_t)-lbound(analysis_p3d_t)
      if(allocated(analysis_u3d_t))write(3,*)'analysis_u3d_t',1+ubound(analysis_u3d_t)-lbound(analysis_u3d_t)
      if(allocated(analysis_v3d_t))write(3,*)'analysis_v3d_t',1+ubound(analysis_v3d_t)-lbound(analysis_v3d_t)
      if(allocated(p3dmode_cos_t))write(3,*)'p3dmode_cos_t',1+ubound(p3dmode_cos_t)-lbound(p3dmode_cos_t)
      if(allocated(p3dmode_sin_t))write(3,*)'p3dmode_sin_t',1+ubound(p3dmode_sin_t)-lbound(p3dmode_sin_t)
      if(allocated(u3dmode_cos_t))write(3,*)'u3dmode_cos_t',1+ubound(u3dmode_cos_t)-lbound(u3dmode_cos_t)
      if(allocated(u3dmode_sin_t))write(3,*)'u3dmode_sin_t',1+ubound(u3dmode_sin_t)-lbound(u3dmode_sin_t)
      if(allocated(v3dmode_cos_t))write(3,*)'v3dmode_cos_t',1+ubound(v3dmode_cos_t)-lbound(v3dmode_cos_t)
      if(allocated(v3dmode_sin_t))write(3,*)'v3dmode_sin_t',1+ubound(v3dmode_sin_t)-lbound(v3dmode_sin_t)
      if(allocated(analysetide3d_u))write(3,*)'analysetide3d_u',1+ubound(analysetide3d_u)-lbound(analysetide3d_u)
      if(allocated(analysetide3d_v))write(3,*)'analysetide3d_v',1+ubound(analysetide3d_v)-lbound(analysetide3d_v)
      if(allocated(analysetide3d_t))write(3,*)'analysetide3d_t',1+ubound(analysetide3d_t)-lbound(analysetide3d_t)
      if(allocated(vel3dtidecosout_u))write(3,*)'vel3dtidecosout_u',1+ubound(vel3dtidecosout_u)-lbound(vel3dtidecosout_u)
      if(allocated(vel3dtidesinout_u))write(3,*)'vel3dtidesinout_u',1+ubound(vel3dtidesinout_u)-lbound(vel3dtidesinout_u)
      if(allocated(vel3dtidecosout_v))write(3,*)'vel3dtidecosout_v',1+ubound(vel3dtidecosout_v)-lbound(vel3dtidecosout_v)
      if(allocated(vel3dtidesinout_v))write(3,*)'vel3dtidesinout_v',1+ubound(vel3dtidesinout_v)-lbound(vel3dtidesinout_v)
      if(allocated(rhptidecosout_t))write(3,*)'rhptidecosout_t',1+ubound(rhptidecosout_t)-lbound(rhptidecosout_t)
      if(allocated(rhptidesinout_t))write(3,*)'rhptidesinout_t',1+ubound(rhptidesinout_t)-lbound(rhptidesinout_t)
      if(allocated(bio_relax_north))write(3,*)'bio_relax_north',1+ubound(bio_relax_north)-lbound(bio_relax_north)
      if(allocated(bio_relax_south))write(3,*)'bio_relax_south',1+ubound(bio_relax_south)-lbound(bio_relax_south)
      if(allocated(bio_relax_east))write(3,*)'bio_relax_east',1+ubound(bio_relax_east)-lbound(bio_relax_east)
      if(allocated(bio_relax_west))write(3,*)'bio_relax_west',1+ubound(bio_relax_west)-lbound(bio_relax_west)
      if(allocated(bio_relax_full))write(3,*)'bio_relax_full',1+ubound(bio_relax_full)-lbound(bio_relax_full)
      if(allocated(analysetide_u))write(3,*)'analysetide_u',1+ubound(analysetide_u)-lbound(analysetide_u)
      if(allocated(analysetide_v))write(3,*)'analysetide_v',1+ubound(analysetide_v)-lbound(analysetide_v)
      if(allocated(analysetide_w))write(3,*)'analysetide_w',1+ubound(analysetide_w)-lbound(analysetide_w)
      if(allocated(analysetide2_w))write(3,*)'analysetide2_w',1+ubound(analysetide2_w)-lbound(analysetide2_w)
      if(allocated(cocosisi_f))write(3,*)'cocosisi_f',1+ubound(cocosisi_f)-lbound(cocosisi_f)
      if(allocated(gridrotcos_f))write(3,*)'gridrotcos_f',1+ubound(gridrotcos_f)-lbound(gridrotcos_f)
      if(allocated(gridrotsin_f))write(3,*)'gridrotsin_f',1+ubound(gridrotsin_f)-lbound(gridrotsin_f)
      if(allocated(tideanalysismatrix))write(3,*)'tideanalysismatrix',1+ubound(tideanalysismatrix)-lbound(tideanalysismatrix)
      if(allocated(iaveraged_in))write(3,*)'iaveraged_in',1+ubound(iaveraged_in)-lbound(iaveraged_in)
      if(allocated(iaveraged_out))write(3,*)'iaveraged_out',1+ubound(iaveraged_out)-lbound(iaveraged_out)
      if(allocated(iaverag1d_in))write(3,*)'iaverag1d_in',1+ubound(iaverag1d_in)-lbound(iaverag1d_in)
      if(allocated(iaverag1d_out))write(3,*)'iaverag1d_out',1+ubound(iaverag1d_out)-lbound(iaverag1d_out)
      if(allocated(light_kpar2_w))write(3,*)'light_kpar2_w',1+ubound(light_kpar2_w)-lbound(light_kpar2_w)
      if(allocated(tidepotential_w))write(3,*)'tidepotential_w',1+ubound(tidepotential_w)-lbound(tidepotential_w)
      if(allocated(sshtidecos_w))write(3,*)'sshtidecos_w',1+ubound(sshtidecos_w)-lbound(sshtidecos_w)
      if(allocated(sshtidesin_w))write(3,*)'sshtidesin_w',1+ubound(sshtidesin_w)-lbound(sshtidesin_w)
      if(allocated(veltidecos_u))write(3,*)'veltidecos_u',1+ubound(veltidecos_u)-lbound(veltidecos_u)
      if(allocated(veltidesin_u))write(3,*)'veltidesin_u',1+ubound(veltidesin_u)-lbound(veltidesin_u)
      if(allocated(veltidecos_v))write(3,*)'veltidecos_v',1+ubound(veltidecos_v)-lbound(veltidecos_v)
      if(allocated(veltidesin_v))write(3,*)'veltidesin_v',1+ubound(veltidesin_v)-lbound(veltidesin_v)
      if(allocated(potidecos_w))write(3,*)'potidecos_w',1+ubound(potidecos_w)-lbound(potidecos_w)
      if(allocated(potidesin_w))write(3,*)'potidesin_w',1+ubound(potidesin_w)-lbound(potidesin_w)
      if(allocated(veltidecosout_u))write(3,*)'veltidecosout_u',1+ubound(veltidecosout_u)-lbound(veltidecosout_u)
      if(allocated(veltidesinout_u))write(3,*)'veltidesinout_u',1+ubound(veltidesinout_u)-lbound(veltidesinout_u)
      if(allocated(veltidecosout_v))write(3,*)'veltidecosout_v',1+ubound(veltidecosout_v)-lbound(veltidecosout_v)
      if(allocated(veltidesinout_v))write(3,*)'veltidesinout_v',1+ubound(veltidesinout_v)-lbound(veltidesinout_v)
      if(allocated(sshtidecosout_w))write(3,*)'sshtidecosout_w',1+ubound(sshtidecosout_w)-lbound(sshtidecosout_w)
      if(allocated(sshtidesinout_w))write(3,*)'sshtidesinout_w',1+ubound(sshtidesinout_w)-lbound(sshtidesinout_w)
      if(allocated(ssh3dtidecosout_w))write(3,*)'ssh3dtidecosout_w',1+ubound(ssh3dtidecosout_w)-lbound(ssh3dtidecosout_w)
      if(allocated(ssh3dtidesinout_w))write(3,*)'ssh3dtidesinout_w',1+ubound(ssh3dtidesinout_w)-lbound(ssh3dtidesinout_w)
      if(allocated(passetide))write(3,*)'passetide',1+ubound(passetide)-lbound(passetide)
      if(allocated(ftide))write(3,*)'ftide',1+ubound(ftide)-lbound(ftide)
      if(allocated(utide))write(3,*)'utide',1+ubound(utide)-lbound(utide)
      if(allocated(northflux_sumt_v))write(3,*)'northflux_sumt_v',1+ubound(northflux_sumt_v)-lbound(northflux_sumt_v)
      if(allocated(southflux_sumt_v))write(3,*)'southflux_sumt_v',1+ubound(southflux_sumt_v)-lbound(southflux_sumt_v)
      if(allocated(drifter_l))write(3,*)'drifter_l',1+ubound(drifter_l)-lbound(drifter_l)
      if(allocated(velbot3d2d_u))write(3,*)'velbot3d2d_u',1+ubound(velbot3d2d_u)-lbound(velbot3d2d_u)
      if(allocated(velbot3d2d_v))write(3,*)'velbot3d2d_v',1+ubound(velbot3d2d_v)-lbound(velbot3d2d_v)
      if(allocated(velbot_u))write(3,*)'velbot_u',1+ubound(velbot_u)-lbound(velbot_u)
      if(allocated(velbot_v))write(3,*)'velbot_v',1+ubound(velbot_v)-lbound(velbot_v)
      if(allocated(kvectorpeak_i))write(3,*)'kvectorpeak_i',1+ubound(kvectorpeak_i)-lbound(kvectorpeak_i)
      if(allocated(kvectorpeak_j))write(3,*)'kvectorpeak_j',1+ubound(kvectorpeak_j)-lbound(kvectorpeak_j)
      if(allocated(drifter_send_canal))write(3,*)'drifter_send_canal',1+ubound(drifter_send_canal)-lbound(drifter_send_canal)
      if(allocated(drifter_recv_canal))write(3,*)'drifter_recv_canal',1+ubound(drifter_recv_canal)-lbound(drifter_recv_canal)
      if(allocated(qwave_j_w))write(3,*)'qwave_j_w',1+ubound(qwave_j_w)-lbound(qwave_j_w)
      if(allocated(qwave_i_w))write(3,*)'qwave_i_w',1+ubound(qwave_i_w)-lbound(qwave_i_w)
      if(allocated(velwave_i_u))write(3,*)'velwave_i_u',1+ubound(velwave_i_u)-lbound(velwave_i_u)
      if(allocated(velwave_j_u))write(3,*)'velwave_j_u',1+ubound(velwave_j_u)-lbound(velwave_j_u)
      if(allocated(velwave_i_v))write(3,*)'velwave_i_v',1+ubound(velwave_i_v)-lbound(velwave_i_v)
      if(allocated(velwave_j_v))write(3,*)'velwave_j_v',1+ubound(velwave_j_v)-lbound(velwave_j_v)
      if(allocated(q_t))write(3,*)'q_t',1+ubound(q_t)-lbound(q_t)
      if(allocated(sshwave_j_w))write(3,*)'sshwave_j_w',1+ubound(sshwave_j_w)-lbound(sshwave_j_w)
      if(allocated(sshwave_i_w))write(3,*)'sshwave_i_w',1+ubound(sshwave_i_w)-lbound(sshwave_i_w)
      if(allocated(mc))write(3,*)'mc',1+ubound(mc)-lbound(mc)
      if(allocated(anyv3dr4))write(3,*)'anyv3dr4',1+ubound(anyv3dr4)-lbound(anyv3dr4)
      if(allocated(mc_out))write(3,*)'mc_out',1+ubound(mc_out)-lbound(mc_out)
      if(allocated(velbarwave_i_v))write(3,*)'velbarwave_i_v',1+ubound(velbarwave_i_v)-lbound(velbarwave_i_v)
      if(allocated(qavr_t))write(3,*)'qavr_t',1+ubound(qavr_t)-lbound(qavr_t)
      if(allocated(nhpgf_u))write(3,*)'nhpgf_u',1+ubound(nhpgf_u)-lbound(nhpgf_u)
      if(allocated(nhpgf_v))write(3,*)'nhpgf_v',1+ubound(nhpgf_v)-lbound(nhpgf_v)
      if(allocated(dsigr4_t))write(3,*)'dsigr4_t',1+ubound(dsigr4_t)-lbound(dsigr4_t)
      if(allocated(sig1dpom_t))write(3,*)'sig1dpom_t',1+ubound(sig1dpom_t)-lbound(sig1dpom_t)
      if(allocated(obc_q_i))write(3,*)'obc_q_i',1+ubound(obc_q_i)-lbound(obc_q_i)
      if(allocated(obc_ub_i))write(3,*)'obc_ub_i',1+ubound(obc_ub_i)-lbound(obc_ub_i)
      if(allocated(dsig_w))write(3,*)'dsig_w',1+ubound(dsig_w)-lbound(dsig_w)
      if(allocated(retarded_pss_w))write(3,*)'retarded_pss_w',1+ubound(retarded_pss_w)-lbound(retarded_pss_w)
      if(allocated(velwf_u))write(3,*)'velwf_u',1+ubound(velwf_u)-lbound(velwf_u)
      if(allocated(velwf_v))write(3,*)'velwf_v',1+ubound(velwf_v)-lbound(velwf_v)
      if(allocated(kw_w))write(3,*)'kw_w',1+ubound(kw_w)-lbound(kw_w)
      if(allocated(q2davr_w))write(3,*)'q2davr_w',1+ubound(q2davr_w)-lbound(q2davr_w)
      if(allocated(breaker2d_t))write(3,*)'breaker2d_t',1+ubound(breaker2d_t)-lbound(breaker2d_t)
      if(allocated(invdx_t))write(3,*)'invdx_t',1+ubound(invdx_t)-lbound(invdx_t)
      if(allocated(invdx_f))write(3,*)'invdx_f',1+ubound(invdx_f)-lbound(invdx_f)
      if(allocated(invdx_u))write(3,*)'invdx_u',1+ubound(invdx_u)-lbound(invdx_u)
      if(allocated(invdy_u))write(3,*)'invdy_u',1+ubound(invdy_u)-lbound(invdy_u)
      if(allocated(invdy_t))write(3,*)'invdy_t',1+ubound(invdy_t)-lbound(invdy_t)
      if(allocated(invdy_f))write(3,*)'invdy_f',1+ubound(invdy_f)-lbound(invdy_f)
      if(allocated(invdy_v))write(3,*)'invdy_v',1+ubound(invdy_v)-lbound(invdy_v)
      if(allocated(invdx_v))write(3,*)'invdx_v',1+ubound(invdx_v)-lbound(invdx_v)
      if(allocated(invdxdy_t))write(3,*)'invdxdy_t',1+ubound(invdxdy_t)-lbound(invdxdy_t)
      if(allocated(invdxdy_u))write(3,*)'invdxdy_u',1+ubound(invdxdy_u)-lbound(invdxdy_u)
      if(allocated(invdxdy_v))write(3,*)'invdxdy_v',1+ubound(invdxdy_v)-lbound(invdxdy_v)
      if(allocated(nhpgf2d_u))write(3,*)'nhpgf2d_u',1+ubound(nhpgf2d_u)-lbound(nhpgf2d_u)
      if(allocated(nhpgf2d_v))write(3,*)'nhpgf2d_v',1+ubound(nhpgf2d_v)-lbound(nhpgf2d_v)
      if(allocated(sshmax_w))write(3,*)'sshmax_w',1+ubound(sshmax_w)-lbound(sshmax_w)
      if(allocated(sshmin_w))write(3,*)'sshmin_w',1+ubound(sshmin_w)-lbound(sshmin_w)
      if(allocated(sshmin_tmp_w))write(3,*)'sshmin_tmp_w',1+ubound(sshmin_tmp_w)-lbound(sshmin_tmp_w)
      if(allocated(rwavebreak_t))write(3,*)'rwavebreak_t',1+ubound(rwavebreak_t)-lbound(rwavebreak_t)
      if(allocated(breaker_t))write(3,*)'breaker_t',1+ubound(breaker_t)-lbound(breaker_t)
      if(allocated(hs_w))write(3,*)'hs_w',1+ubound(hs_w)-lbound(hs_w)
      if(allocated(cx_upper_t))write(3,*)'cx_upper_t',1+ubound(cx_upper_t)-lbound(cx_upper_t)
      if(allocated(cy_upper_t))write(3,*)'cy_upper_t',1+ubound(cy_upper_t)-lbound(cy_upper_t)
      if(allocated(c_lower_t))write(3,*)'c_lower_t',1+ubound(c_lower_t)-lbound(c_lower_t)
      if(allocated(slopemax_u))write(3,*)'slopemax_u',1+ubound(slopemax_u)-lbound(slopemax_u)
      if(allocated(pgfwave_u))write(3,*)'pgfwave_u',1+ubound(pgfwave_u)-lbound(pgfwave_u)
      if(allocated(sshrelax_j_w))write(3,*)'sshrelax_j_w',1+ubound(sshrelax_j_w)-lbound(sshrelax_j_w)
      if(allocated(sshrelax_i_w))write(3,*)'sshrelax_i_w',1+ubound(sshrelax_i_w)-lbound(sshrelax_i_w)
      write(3,*)'asselin_nh'
      write(3,*)'wetdry_cstnh'
      write(3,*)'cfl_nh'
      write(3,*)'cwavepeak'
      write(3,*)'periodpeak'
      write(3,*)'freq2pipeak'
      write(3,*)'kvectorpeak'
      write(3,*)'kvector'
      write(3,*)'mukvector'
      write(3,*)'nhpgf_reduce'
      write(3,*)'acm_speed'
      write(3,*)'inv_brkh'
      write(3,*)'inv_brkslope2'
      write(3,*)'inv_brkvar'
      write(3,*)'brk_crit_slope'
      write(3,*)'brk_crit_var'
      write(3,*)'brk_crit_h'
      write(3,*)'brk_crit_r'
      write(3,*)'nh2d_graph_period'
      write(3,*)'nh2d_graph_spinup'
      write(3,*)'sum_avr_hs'
      write(3,*)'nh_frozensigma'
      write(3,*)'nh_wavebreakscheme'
      write(3,*)'flag_adve2d'
      write(3,*)'flag_adve3d'
      write(3,*)'flag_nh3d'
      write(3,*)'flag_nh3d_none'
      write(3,*)'flag_nh3d_nosplit_uv'
      write(3,*)'flag_nh3d_nosplit_tsuv'
      write(3,*)'flag_nh3d_timesplit_tsuv'
      write(3,*)'flag_timesplitting'
      write(3,*)'flag_timesplitting_adve_uv'
      write(3,*)'flagspo_i1'
      write(3,*)'flagspo_i2'
      write(3,*)'flagspo_j1'
      write(3,*)'flagspo_j2'
      write(3,*)'flag_ksloffline'
      write(3,*)'flag_groundwater'
      write(3,*)'flag_surfriver'
      write(3,*)'ofl_reversedtime'
      write(3,*)'obcnh_scheme'
      write(3,*)'obc_scheme'
      write(3,*)'drifter_random_w'
      write(3,*)'flag_buo_w'
      write(3,*)'flag_ogcmtidemixing'
      write(3,*)'flag_ogcm_instab'
      write(3,*)'flag_ofactor'
      write(3,*)'flag_omega_cumul'
      write(3,*)'flag_negdif_ver'
      write(3,*)'flag_negdif_hor'
      write(3,*)'flag_z0_macro'
      write(3,*)'oasis_symsym_onoff'
      write(3,*)'oasis_symsym_retrots'
      if(allocated(variance_timeaveraged))write(3,*)'variance_timeaveraged',1+ubound(variance_timeaveraged)-lbound(variance_timeaveraged)
      if(allocated(ssh_timeaveraged))write(3,*)'ssh_timeaveraged',1+ubound(ssh_timeaveraged)-lbound(ssh_timeaveraged)
      if(allocated(flx2d_timeaveraged))write(3,*)'flx2d_timeaveraged',1+ubound(flx2d_timeaveraged)-lbound(flx2d_timeaveraged)
      if(allocated(fly2d_timeaveraged))write(3,*)'fly2d_timeaveraged',1+ubound(fly2d_timeaveraged)-lbound(fly2d_timeaveraged)
      if(allocated(flx3d_timeaveraged))write(3,*)'flx3d_timeaveraged',1+ubound(flx3d_timeaveraged)-lbound(flx3d_timeaveraged)
      if(allocated(fly3d_timeaveraged))write(3,*)'fly3d_timeaveraged',1+ubound(fly3d_timeaveraged)-lbound(fly3d_timeaveraged)
      if(allocated(u_euler_timeaveraged))write(3,*)'u_euler_timeaveraged',1+ubound(u_euler_timeaveraged)-lbound(u_euler_timeaveraged)
      if(allocated(v_euler_timeaveraged))write(3,*)'v_euler_timeaveraged',1+ubound(v_euler_timeaveraged)-lbound(v_euler_timeaveraged)
      if(allocated(tkeb_w))write(3,*)'tkeb_w',1+ubound(tkeb_w)-lbound(tkeb_w)
      if(allocated(tken_w))write(3,*)'tken_w',1+ubound(tken_w)-lbound(tken_w)
      if(allocated(tkea_w))write(3,*)'tkea_w',1+ubound(tkea_w)-lbound(tkea_w)
      if(allocated(tkle_w))write(3,*)'tkle_w',1+ubound(tkle_w)-lbound(tkle_w)
      if(allocated(tkll_w))write(3,*)'tkll_w',1+ubound(tkll_w)-lbound(tkll_w)
      if(allocated(epsb_w))write(3,*)'epsb_w',1+ubound(epsb_w)-lbound(epsb_w)
      if(allocated(epsn_w))write(3,*)'epsn_w',1+ubound(epsn_w)-lbound(epsn_w)
      if(allocated(epsa_w))write(3,*)'epsa_w',1+ubound(epsa_w)-lbound(epsa_w)
      if(allocated(gradssh_u))write(3,*)'gradssh_u',1+ubound(gradssh_u)-lbound(gradssh_u)
      if(allocated(gradssh_v))write(3,*)'gradssh_v',1+ubound(gradssh_v)-lbound(gradssh_v)
      if(allocated(pgf_u))write(3,*)'pgf_u',1+ubound(pgf_u)-lbound(pgf_u)
      if(allocated(pgf_v))write(3,*)'pgf_v',1+ubound(pgf_v)-lbound(pgf_v)
      if(allocated(dz_u))write(3,*)'dz_u',1+ubound(dz_u)-lbound(dz_u)
      if(allocated(dz_v))write(3,*)'dz_v',1+ubound(dz_v)-lbound(dz_v)
      if(allocated(dz_t))write(3,*)'dz_t',1+ubound(dz_t)-lbound(dz_t)
      if(allocated(depth_t))write(3,*)'depth_t',1+ubound(depth_t)-lbound(depth_t)
      if(allocated(depth_u))write(3,*)'depth_u',1+ubound(depth_u)-lbound(depth_u)
      if(allocated(depth_v))write(3,*)'depth_v',1+ubound(depth_v)-lbound(depth_v)
      if(allocated(depth_w))write(3,*)'depth_w',1+ubound(depth_w)-lbound(depth_w)
      if(allocated(depth_f))write(3,*)'depth_f',1+ubound(depth_f)-lbound(depth_f)
      if(allocated(stokesforces_u))write(3,*)'stokesforces_u',1+ubound(stokesforces_u)-lbound(stokesforces_u)
      if(allocated(stokesforces_v))write(3,*)'stokesforces_v',1+ubound(stokesforces_v)-lbound(stokesforces_v)
      if(allocated(r_mangrove))write(3,*)'r_mangrove',1+ubound(r_mangrove)-lbound(r_mangrove)
      if(allocated(sigma_w))write(3,*)'sigma_w',1+ubound(sigma_w)-lbound(sigma_w)
      if(allocated(dsig_f))write(3,*)'dsig_f',1+ubound(dsig_f)-lbound(dsig_f)
      if(allocated(rhp_t))write(3,*)'rhp_t',1+ubound(rhp_t)-lbound(rhp_t)
      if(allocated(rhcref_t))write(3,*)'rhcref_t',1+ubound(rhcref_t)-lbound(rhcref_t)
      if(allocated(vel_u))write(3,*)'vel_u',1+ubound(vel_u)-lbound(vel_u)
      if(allocated(vel_v))write(3,*)'vel_v',1+ubound(vel_v)-lbound(vel_v)
      if(allocated(presgrad_u))write(3,*)'presgrad_u',1+ubound(presgrad_u)-lbound(presgrad_u)
      if(allocated(presgrad_v))write(3,*)'presgrad_v',1+ubound(presgrad_v)-lbound(presgrad_v)
      if(allocated(omega_w))write(3,*)'omega_w',1+ubound(omega_w)-lbound(omega_w)
      if(allocated(veldydz_u))write(3,*)'veldydz_u',1+ubound(veldydz_u)-lbound(veldydz_u)
      if(allocated(veldxdz_v))write(3,*)'veldxdz_v',1+ubound(veldxdz_v)-lbound(veldxdz_v)
      if(allocated(tem_t))write(3,*)'tem_t',1+ubound(tem_t)-lbound(tem_t)
      if(allocated(sal_t))write(3,*)'sal_t',1+ubound(sal_t)-lbound(sal_t)
      if(allocated(tridia_in))write(3,*)'tridia_in',1+ubound(tridia_in)-lbound(tridia_in)
      if(allocated(hr_z2lr_w))write(3,*)'hr_z2lr_w',1+ubound(hr_z2lr_w)-lbound(hr_z2lr_w)
      if(allocated(anyv3d))write(3,*)'anyv3d',1+ubound(anyv3d)-lbound(anyv3d)
      if(allocated(sponge_t))write(3,*)'sponge_t',1+ubound(sponge_t)-lbound(sponge_t)
      if(allocated(sponge_u))write(3,*)'sponge_u',1+ubound(sponge_u)-lbound(sponge_u)
      if(allocated(sponge_v))write(3,*)'sponge_v',1+ubound(sponge_v)-lbound(sponge_v)
      if(allocated(uwind_t))write(3,*)'uwind_t',1+ubound(uwind_t)-lbound(uwind_t)
      if(allocated(vwind_t))write(3,*)'vwind_t',1+ubound(vwind_t)-lbound(vwind_t)
      if(allocated(uwind100_t))write(3,*)'uwind100_t',1+ubound(uwind100_t)-lbound(uwind100_t)
      if(allocated(vwind100_t))write(3,*)'vwind100_t',1+ubound(vwind100_t)-lbound(vwind100_t)
      if(allocated(sshf_w))write(3,*)'sshf_w',1+ubound(sshf_w)-lbound(sshf_w)
      if(allocated(slhf_w))write(3,*)'slhf_w',1+ubound(slhf_w)-lbound(slhf_w)
      if(allocated(ssr_w))write(3,*)'ssr_w',1+ubound(ssr_w)-lbound(ssr_w)
      if(allocated(ssr24prv_w))write(3,*)'ssr24prv_w',1+ubound(ssr24prv_w)-lbound(ssr24prv_w)
      if(allocated(snsf_w))write(3,*)'snsf_w',1+ubound(snsf_w)-lbound(snsf_w)
      if(allocated(precipi_w))write(3,*)'precipi_w',1+ubound(precipi_w)-lbound(precipi_w)
      if(allocated(taux_w))write(3,*)'taux_w',1+ubound(taux_w)-lbound(taux_w)
      if(allocated(tauy_w))write(3,*)'tauy_w',1+ubound(tauy_w)-lbound(tauy_w)
      if(allocated(sigma_f))write(3,*)'sigma_f',1+ubound(sigma_f)-lbound(sigma_f)
      if(allocated(km_w))write(3,*)'km_w',1+ubound(km_w)-lbound(km_w)
      if(allocated(kh_w))write(3,*)'kh_w',1+ubound(kh_w)-lbound(kh_w)
      if(allocated(rho_t))write(3,*)'rho_t',1+ubound(rho_t)-lbound(rho_t)
      if(allocated(omega_evaprec_w))write(3,*)'omega_evaprec_w',1+ubound(omega_evaprec_w)-lbound(omega_evaprec_w)
      if(allocated(tridia_out))write(3,*)'tridia_out',1+ubound(tridia_out)-lbound(tridia_out)
      if(allocated(fluxbar_sumt_u))write(3,*)'fluxbar_sumt_u',1+ubound(fluxbar_sumt_u)-lbound(fluxbar_sumt_u)
      if(allocated(fluxbar_sumt_v))write(3,*)'fluxbar_sumt_v',1+ubound(fluxbar_sumt_v)-lbound(fluxbar_sumt_v)
      if(allocated(velbar_u))write(3,*)'velbar_u',1+ubound(velbar_u)-lbound(velbar_u)
      if(allocated(uflux2d_f))write(3,*)'uflux2d_f',1+ubound(uflux2d_f)-lbound(uflux2d_f)
      if(allocated(vflux2d_f))write(3,*)'vflux2d_f',1+ubound(vflux2d_f)-lbound(vflux2d_f)
      if(allocated(vlxbar_f))write(3,*)'vlxbar_f',1+ubound(vlxbar_f)-lbound(vlxbar_f)
      if(allocated(velbar_v))write(3,*)'velbar_v',1+ubound(velbar_v)-lbound(velbar_v)
      if(allocated(vlybar_f))write(3,*)'vlybar_f',1+ubound(vlybar_f)-lbound(vlybar_f)
      if(allocated(flux2d_u))write(3,*)'flux2d_u',1+ubound(flux2d_u)-lbound(flux2d_u)
      if(allocated(flux2d_v))write(3,*)'flux2d_v',1+ubound(flux2d_v)-lbound(flux2d_v)
      if(allocated(ssh_w))write(3,*)'ssh_w',1+ubound(ssh_w)-lbound(ssh_w)
      if(allocated(hssh_f))write(3,*)'hssh_f',1+ubound(hssh_f)-lbound(hssh_f)
      if(allocated(hz_u))write(3,*)'hz_u',1+ubound(hz_u)-lbound(hz_u)
      if(allocated(hz_v))write(3,*)'hz_v',1+ubound(hz_v)-lbound(hz_v)
      if(allocated(hz_w))write(3,*)'hz_w',1+ubound(hz_w)-lbound(hz_w)
      if(allocated(heatrelax_w))write(3,*)'heatrelax_w',1+ubound(heatrelax_w)-lbound(heatrelax_w)
      if(allocated(xy_u))write(3,*)'xy_u',1+ubound(xy_u)-lbound(xy_u)
      if(allocated(xy_v))write(3,*)'xy_v',1+ubound(xy_v)-lbound(xy_v)
      if(allocated(xy_t))write(3,*)'xy_t',1+ubound(xy_t)-lbound(xy_t)
      if(allocated(xy_f))write(3,*)'xy_f',1+ubound(xy_f)-lbound(xy_f)
      if(allocated(pss_w))write(3,*)'pss_w',1+ubound(pss_w)-lbound(pss_w)
      if(allocated(wstress_u))write(3,*)'wstress_u',1+ubound(wstress_u)-lbound(wstress_u)
      if(allocated(wstress_v))write(3,*)'wstress_v',1+ubound(wstress_v)-lbound(wstress_v)
      if(allocated(ssh_int_w))write(3,*)'ssh_int_w',1+ubound(ssh_int_w)-lbound(ssh_int_w)
      if(allocated(ustokesvortex_t))write(3,*)'ustokesvortex_t',1+ubound(ustokesvortex_t)-lbound(ustokesvortex_t)
      if(allocated(ustokesvortex_f))write(3,*)'ustokesvortex_f',1+ubound(ustokesvortex_f)-lbound(ustokesvortex_f)
      if(allocated(vstokesvortex_t))write(3,*)'vstokesvortex_t',1+ubound(vstokesvortex_t)-lbound(vstokesvortex_t)
      if(allocated(vstokesvortex_f))write(3,*)'vstokesvortex_f',1+ubound(vstokesvortex_f)-lbound(vstokesvortex_f)
      if(allocated(velavr_u))write(3,*)'velavr_u',1+ubound(velavr_u)-lbound(velavr_u)
      if(allocated(velavr_v))write(3,*)'velavr_v',1+ubound(velavr_v)-lbound(velavr_v)
      if(allocated(timefilter_u))write(3,*)'timefilter_u',1+ubound(timefilter_u)-lbound(timefilter_u)
      if(allocated(timefilter_v))write(3,*)'timefilter_v',1+ubound(timefilter_v)-lbound(timefilter_v)
      if(allocated(teta0_t))write(3,*)'teta0_t',1+ubound(teta0_t)-lbound(teta0_t)
      if(allocated(teta2_t))write(3,*)'teta2_t',1+ubound(teta2_t)-lbound(teta2_t)
      if(allocated(q2_t))write(3,*)'q2_t',1+ubound(q2_t)-lbound(q2_t)
      if(allocated(teta2delta_t))write(3,*)'teta2delta_t',1+ubound(teta2delta_t)-lbound(teta2delta_t)
      if(allocated(q2delta_t))write(3,*)'q2delta_t',1+ubound(q2delta_t)-lbound(q2delta_t)
      if(allocated(uwinddelta_t))write(3,*)'uwinddelta_t',1+ubound(uwinddelta_t)-lbound(uwinddelta_t)
      if(allocated(vwinddelta_t))write(3,*)'vwinddelta_t',1+ubound(vwinddelta_t)-lbound(vwinddelta_t)
      if(allocated(tetastar_t))write(3,*)'tetastar_t',1+ubound(tetastar_t)-lbound(tetastar_t)
      if(allocated(ustar_t))write(3,*)'ustar_t',1+ubound(ustar_t)-lbound(ustar_t)
      if(allocated(qstar_t))write(3,*)'qstar_t',1+ubound(qstar_t)-lbound(qstar_t)
      if(allocated(frozenterm2d_u))write(3,*)'frozenterm2d_u',1+ubound(frozenterm2d_u)-lbound(frozenterm2d_u)
      if(allocated(frozenterm2d_v))write(3,*)'frozenterm2d_v',1+ubound(frozenterm2d_v)-lbound(frozenterm2d_v)
      if(allocated(frozenterm3d_u))write(3,*)'frozenterm3d_u',1+ubound(frozenterm3d_u)-lbound(frozenterm3d_u)
      if(allocated(frozenterm3d_v))write(3,*)'frozenterm3d_v',1+ubound(frozenterm3d_v)-lbound(frozenterm3d_v)
      if(allocated(fluxbar_u))write(3,*)'fluxbar_u',1+ubound(fluxbar_u)-lbound(fluxbar_u)
      if(allocated(fluxbar_v))write(3,*)'fluxbar_v',1+ubound(fluxbar_v)-lbound(fluxbar_v)
      if(allocated(lon_t))write(3,*)'lon_t',1+ubound(lon_t)-lbound(lon_t)
      if(allocated(lat_t))write(3,*)'lat_t',1+ubound(lat_t)-lbound(lat_t)
      if(allocated(lon_u))write(3,*)'lon_u',1+ubound(lon_u)-lbound(lon_u)
      if(allocated(lat_u))write(3,*)'lat_u',1+ubound(lat_u)-lbound(lat_u)
      if(allocated(lon_v))write(3,*)'lon_v',1+ubound(lon_v)-lbound(lon_v)
      if(allocated(lat_v))write(3,*)'lat_v',1+ubound(lat_v)-lbound(lat_v)
      if(allocated(lon_f))write(3,*)'lon_f',1+ubound(lon_f)-lbound(lon_f)
      if(allocated(lat_f))write(3,*)'lat_f',1+ubound(lat_f)-lbound(lat_f)
      if(allocated(globlon_t))write(3,*)'globlon_t',1+ubound(globlon_t)-lbound(globlon_t)
      if(allocated(globlat_t))write(3,*)'globlat_t',1+ubound(globlat_t)-lbound(globlat_t)
      if(allocated(globlon_u))write(3,*)'globlon_u',1+ubound(globlon_u)-lbound(globlon_u)
      if(allocated(globlat_u))write(3,*)'globlat_u',1+ubound(globlat_u)-lbound(globlat_u)
      if(allocated(globlon_v))write(3,*)'globlon_v',1+ubound(globlon_v)-lbound(globlon_v)
      if(allocated(globlat_v))write(3,*)'globlat_v',1+ubound(globlat_v)-lbound(globlat_v)
      if(allocated(fluxbarsum_ofl_u))write(3,*)'fluxbarsum_ofl_u',1+ubound(fluxbarsum_ofl_u)-lbound(fluxbarsum_ofl_u)
      if(allocated(fluxbarsum_ofl_v))write(3,*)'fluxbarsum_ofl_v',1+ubound(fluxbarsum_ofl_v)-lbound(fluxbarsum_ofl_v)
      if(allocated(dxdy_u))write(3,*)'dxdy_u',1+ubound(dxdy_u)-lbound(dxdy_u)
      if(allocated(dxdy_v))write(3,*)'dxdy_v',1+ubound(dxdy_v)-lbound(dxdy_v)
      if(allocated(dxdy_t))write(3,*)'dxdy_t',1+ubound(dxdy_t)-lbound(dxdy_t)
      if(allocated(dxdy_e_t))write(3,*)'dxdy_e_t',1+ubound(dxdy_e_t)-lbound(dxdy_e_t)
      if(allocated(dx_e_f))write(3,*)'dx_e_f',1+ubound(dx_e_f)-lbound(dx_e_f)
      if(allocated(dy_e_f))write(3,*)'dy_e_f',1+ubound(dy_e_f)-lbound(dy_e_f)
      if(allocated(dx_u))write(3,*)'dx_u',1+ubound(dx_u)-lbound(dx_u)
      if(allocated(dx_v))write(3,*)'dx_v',1+ubound(dx_v)-lbound(dx_v)
      if(allocated(dx_t))write(3,*)'dx_t',1+ubound(dx_t)-lbound(dx_t)
      if(allocated(dx_f))write(3,*)'dx_f',1+ubound(dx_f)-lbound(dx_f)
      if(allocated(dy_u))write(3,*)'dy_u',1+ubound(dy_u)-lbound(dy_u)
      if(allocated(dy_v))write(3,*)'dy_v',1+ubound(dy_v)-lbound(dy_v)
      if(allocated(dy_t))write(3,*)'dy_t',1+ubound(dy_t)-lbound(dy_t)
      if(allocated(dy_f))write(3,*)'dy_f',1+ubound(dy_f)-lbound(dy_f)
      if(allocated(h_w))write(3,*)'h_w',1+ubound(h_w)-lbound(h_w)
      if(allocated(hcopy_w))write(3,*)'hcopy_w',1+ubound(hcopy_w)-lbound(hcopy_w)
      if(allocated(h0_w))write(3,*)'h0_w',1+ubound(h0_w)-lbound(h0_w)
      if(allocated(h_u))write(3,*)'h_u',1+ubound(h_u)-lbound(h_u)
      if(allocated(h_v))write(3,*)'h_v',1+ubound(h_v)-lbound(h_v)
      if(allocated(h_f))write(3,*)'h_f',1+ubound(h_f)-lbound(h_f)
      if(allocated(coriolis_t))write(3,*)'coriolis_t',1+ubound(coriolis_t)-lbound(coriolis_t)
      if(allocated(coriolis_f))write(3,*)'coriolis_f',1+ubound(coriolis_f)-lbound(coriolis_f)
      if(allocated(rhpzavr_w))write(3,*)'rhpzavr_w',1+ubound(rhpzavr_w)-lbound(rhpzavr_w)
      if(allocated(q10_t))write(3,*)'q10_t',1+ubound(q10_t)-lbound(q10_t)
      if(allocated(teta10_t))write(3,*)'teta10_t',1+ubound(teta10_t)-lbound(teta10_t)
      if(allocated(fric_u))write(3,*)'fric_u',1+ubound(fric_u)-lbound(fric_u)
      if(allocated(fric_v))write(3,*)'fric_v',1+ubound(fric_v)-lbound(fric_v)
      if(allocated(fric_t))write(3,*)'fric_t',1+ubound(fric_t)-lbound(fric_t)
      if(allocated(cdb_t))write(3,*)'cdb_t',1+ubound(cdb_t)-lbound(cdb_t)
      if(allocated(cdb_f))write(3,*)'cdb_f',1+ubound(cdb_f)-lbound(cdb_f)
      if(allocated(xflux_t))write(3,*)'xflux_t',1+ubound(xflux_t)-lbound(xflux_t)
      if(allocated(xflux_f))write(3,*)'xflux_f',1+ubound(xflux_f)-lbound(xflux_f)
      if(allocated(yflux_t))write(3,*)'yflux_t',1+ubound(yflux_t)-lbound(yflux_t)
      if(allocated(yflux_f))write(3,*)'yflux_f',1+ubound(yflux_f)-lbound(yflux_f)
      if(allocated(pres3d2d_u))write(3,*)'pres3d2d_u',1+ubound(pres3d2d_u)-lbound(pres3d2d_u)
      if(allocated(pres3d2d_v))write(3,*)'pres3d2d_v',1+ubound(pres3d2d_v)-lbound(pres3d2d_v)
      if(allocated(adve3d2d_u))write(3,*)'adve3d2d_u',1+ubound(adve3d2d_u)-lbound(adve3d2d_u)
      if(allocated(adve3d2d_v))write(3,*)'adve3d2d_v',1+ubound(adve3d2d_v)-lbound(adve3d2d_v)
      if(allocated(rmangrovebar))write(3,*)'rmangrovebar',1+ubound(rmangrovebar)-lbound(rmangrovebar)
      if(allocated(mang3dto2d_u))write(3,*)'mang3dto2d_u',1+ubound(mang3dto2d_u)-lbound(mang3dto2d_u)
      if(allocated(mang3dto2d_v))write(3,*)'mang3dto2d_v',1+ubound(mang3dto2d_v)-lbound(mang3dto2d_v)
      if(allocated(restoring3d2d_u))write(3,*)'restoring3d2d_u',1+ubound(restoring3d2d_u)-lbound(restoring3d2d_u)
      if(allocated(restoring3d2d_v))write(3,*)'restoring3d2d_v',1+ubound(restoring3d2d_v)-lbound(restoring3d2d_v)
      if(allocated(stokesforces3d2d_u))write(3,*)'stokesforces3d2d_u',1+ubound(stokesforces3d2d_u)-lbound(stokesforces3d2d_u)
      if(allocated(stokesforces3d2d_v))write(3,*)'stokesforces3d2d_v',1+ubound(stokesforces3d2d_v)-lbound(stokesforces3d2d_v)
      if(allocated(wstress_w))write(3,*)'wstress_w',1+ubound(wstress_w)-lbound(wstress_w)
      if(allocated(z0_w))write(3,*)'z0_w',1+ubound(z0_w)-lbound(z0_w)
      if(allocated(albedo_w))write(3,*)'albedo_w',1+ubound(albedo_w)-lbound(albedo_w)
      if(allocated(gridrotcos_t))write(3,*)'gridrotcos_t',1+ubound(gridrotcos_t)-lbound(gridrotcos_t)
      if(allocated(gridrotsin_t))write(3,*)'gridrotsin_t',1+ubound(gridrotsin_t)-lbound(gridrotsin_t)
      if(allocated(grid_angle_t))write(3,*)'grid_angle_t',1+ubound(grid_angle_t)-lbound(grid_angle_t)
      if(allocated(sshstokes_w))write(3,*)'sshstokes_w',1+ubound(sshstokes_w)-lbound(sshstokes_w)
      if(allocated(cwi_int_u))write(3,*)'cwi_int_u',1+ubound(cwi_int_u)-lbound(cwi_int_u)
      if(allocated(cwi_int_v))write(3,*)'cwi_int_v',1+ubound(cwi_int_v)-lbound(cwi_int_v)
      if(allocated(cwj_int_u))write(3,*)'cwj_int_u',1+ubound(cwj_int_u)-lbound(cwj_int_u)
      if(allocated(cwj_int_v))write(3,*)'cwj_int_v',1+ubound(cwj_int_v)-lbound(cwj_int_v)
      if(allocated(sshrefobc_i))write(3,*)'sshrefobc_i',1+ubound(sshrefobc_i)-lbound(sshrefobc_i)
      if(allocated(vbrrefobc_i))write(3,*)'vbrrefobc_i',1+ubound(vbrrefobc_i)-lbound(vbrrefobc_i)
      if(allocated(sshrefobc_j))write(3,*)'sshrefobc_j',1+ubound(sshrefobc_j)-lbound(sshrefobc_j)
      if(allocated(vbrrefobc_j))write(3,*)'vbrrefobc_j',1+ubound(vbrrefobc_j)-lbound(vbrrefobc_j)
      if(allocated(anyv1d))write(3,*)'anyv1d',1+ubound(anyv1d)-lbound(anyv1d)
      if(allocated(gridcard))write(3,*)'gridcard',1+ubound(gridcard)-lbound(gridcard)
      if(allocated(novector))write(3,*)'novector',1+ubound(novector)-lbound(novector)
      if(allocated(proficard))write(3,*)'proficard',1+ubound(proficard)-lbound(proficard)
      if(allocated(airseainfo))write(3,*)'airseainfo',1+ubound(airseainfo)-lbound(airseainfo)
      if(allocated(airseadt))write(3,*)'airseadt',1+ubound(airseadt)-lbound(airseadt)
      if(allocated(riverdt))write(3,*)'riverdt',1+ubound(riverdt)-lbound(riverdt)
      if(allocated(river_t))write(3,*)'river_t',1+ubound(river_t)-lbound(river_t)
      if(allocated(riverflux))write(3,*)'riverflux',1+ubound(riverflux)-lbound(riverflux)
      if(allocated(mask_t))write(3,*)'mask_t',1+ubound(mask_t)-lbound(mask_t)
      if(allocated(mask_f))write(3,*)'mask_f',1+ubound(mask_f)-lbound(mask_f)
      if(allocated(mask_u))write(3,*)'mask_u',1+ubound(mask_u)-lbound(mask_u)
      if(allocated(mask_v))write(3,*)'mask_v',1+ubound(mask_v)-lbound(mask_v)
      if(allocated(mask_vqs_tke_w))write(3,*)'mask_vqs_tke_w',1+ubound(mask_vqs_tke_w)-lbound(mask_vqs_tke_w)
      if(allocated(canaldir))write(3,*)'canaldir',1+ubound(canaldir)-lbound(canaldir)
      if(allocated(wetmask_wi_t))write(3,*)'wetmask_wi_t',1+ubound(wetmask_wi_t)-lbound(wetmask_wi_t)
      if(allocated(lonlat2ij_t))write(3,*)'lonlat2ij_t',1+ubound(lonlat2ij_t)-lbound(lonlat2ij_t)
      if(allocated(canalmpioverlap))write(3,*)'canalmpioverlap',1+ubound(canalmpioverlap)-lbound(canalmpioverlap)
      if(allocated(sodate))write(3,*)'sodate',1+ubound(sodate)-lbound(sodate)
      if(allocated(canalrank))write(3,*)'canalrank',1+ubound(canalrank)-lbound(canalrank)
      if(allocated(canalrankbis))write(3,*)'canalrankbis',1+ubound(canalrankbis)-lbound(canalrankbis)
      if(allocated(i_canalcoord))write(3,*)'i_canalcoord',1+ubound(i_canalcoord)-lbound(i_canalcoord)
      if(allocated(j_canalcoord))write(3,*)'j_canalcoord',1+ubound(j_canalcoord)-lbound(j_canalcoord)
      if(allocated(kmin_u))write(3,*)'kmin_u',1+ubound(kmin_u)-lbound(kmin_u)
      if(allocated(kmin_v))write(3,*)'kmin_v',1+ubound(kmin_v)-lbound(kmin_v)
      if(allocated(kmin_w))write(3,*)'kmin_w',1+ubound(kmin_w)-lbound(kmin_w)
      if(allocated(kundermin_t))write(3,*)'kundermin_t',1+ubound(kundermin_t)-lbound(kundermin_t)
      if(allocated(kundermin_u))write(3,*)'kundermin_u',1+ubound(kundermin_u)-lbound(kundermin_u)
      if(allocated(kundermin_v))write(3,*)'kundermin_v',1+ubound(kundermin_v)-lbound(kundermin_v)
      if(allocated(kmerged_t))write(3,*)'kmerged_t',1+ubound(kmerged_t)-lbound(kmerged_t)
      if(allocated(kmerged_u))write(3,*)'kmerged_u',1+ubound(kmerged_u)-lbound(kmerged_u)
      if(allocated(kmerged_v))write(3,*)'kmerged_v',1+ubound(kmerged_v)-lbound(kmerged_v)
      if(allocated(ksl_t))write(3,*)'ksl_t',1+ubound(ksl_t)-lbound(ksl_t)
      if(allocated(glob_mask_mangrove))write(3,*)'glob_mask_mangrove',1+ubound(glob_mask_mangrove)-lbound(glob_mask_mangrove)
      if(allocated(mask_mangrove_t))write(3,*)'mask_mangrove_t',1+ubound(mask_mangrove_t)-lbound(mask_mangrove_t)
      if(allocated(mask_wave_t))write(3,*)'mask_wave_t',1+ubound(mask_wave_t)-lbound(mask_wave_t)
      if(allocated(upwindriver_t))write(3,*)'upwindriver_t',1+ubound(upwindriver_t)-lbound(upwindriver_t)
      if(allocated(upwindwetdry_t))write(3,*)'upwindwetdry_t',1+ubound(upwindwetdry_t)-lbound(upwindwetdry_t)
      if(allocated(kmergedr4_u))write(3,*)'kmergedr4_u',1+ubound(kmergedr4_u)-lbound(kmergedr4_u)
      if(allocated(kmergedr4_v))write(3,*)'kmergedr4_v',1+ubound(kmergedr4_v)-lbound(kmergedr4_v)
      if(allocated(dsigmerged_u))write(3,*)'dsigmerged_u',1+ubound(dsigmerged_u)-lbound(dsigmerged_u)
      if(allocated(dsigmerged_v))write(3,*)'dsigmerged_v',1+ubound(dsigmerged_v)-lbound(dsigmerged_v)
      if(allocated(pgfratio_u))write(3,*)'pgfratio_u',1+ubound(pgfratio_u)-lbound(pgfratio_u)
      if(allocated(pgfratio_v))write(3,*)'pgfratio_v',1+ubound(pgfratio_v)-lbound(pgfratio_v)
      if(allocated(sshr4_w))write(3,*)'sshr4_w',1+ubound(sshr4_w)-lbound(sshr4_w)
      if(allocated(wetmask_u))write(3,*)'wetmask_u',1+ubound(wetmask_u)-lbound(wetmask_u)
      if(allocated(wetmask_v))write(3,*)'wetmask_v',1+ubound(wetmask_v)-lbound(wetmask_v)
      if(allocated(wetmask_t))write(3,*)'wetmask_t',1+ubound(wetmask_t)-lbound(wetmask_t)
      if(allocated(botlevmerged_w))write(3,*)'botlevmerged_w',1+ubound(botlevmerged_w)-lbound(botlevmerged_w)
      if(allocated(maxbotstress_aft_w))write(3,*)'maxbotstress_aft_w',1+ubound(maxbotstress_aft_w)-lbound(maxbotstress_aft_w)
      if(allocated(maxbotstress_bef_w))write(3,*)'maxbotstress_bef_w',1+ubound(maxbotstress_bef_w)-lbound(maxbotstress_bef_w)
      if(allocated(maxbotstress_w))write(3,*)'maxbotstress_w',1+ubound(maxbotstress_w)-lbound(maxbotstress_w)
      if(allocated(stresswave_w))write(3,*)'stresswave_w',1+ubound(stresswave_w)-lbound(stresswave_w)
      if(allocated(stressc_w))write(3,*)'stressc_w',1+ubound(stressc_w)-lbound(stressc_w)
      if(allocated(sqr_hoverg_u))write(3,*)'sqr_hoverg_u',1+ubound(sqr_hoverg_u)-lbound(sqr_hoverg_u)
      if(allocated(sqr_hoverg_v))write(3,*)'sqr_hoverg_v',1+ubound(sqr_hoverg_v)-lbound(sqr_hoverg_v)
      if(allocated(temobc_t))write(3,*)'temobc_t',1+ubound(temobc_t)-lbound(temobc_t)
      if(allocated(salobc_t))write(3,*)'salobc_t',1+ubound(salobc_t)-lbound(salobc_t)
      if(allocated(velobc_u))write(3,*)'velobc_u',1+ubound(velobc_u)-lbound(velobc_u)
      if(allocated(velobc_v))write(3,*)'velobc_v',1+ubound(velobc_v)-lbound(velobc_v)
      if(allocated(temofl_t))write(3,*)'temofl_t',1+ubound(temofl_t)-lbound(temofl_t)
      if(allocated(salofl_t))write(3,*)'salofl_t',1+ubound(salofl_t)-lbound(salofl_t)
      if(allocated(dzofl_t))write(3,*)'dzofl_t',1+ubound(dzofl_t)-lbound(dzofl_t)
      if(allocated(velofl_u))write(3,*)'velofl_u',1+ubound(velofl_u)-lbound(velofl_u)
      if(allocated(velofl_v))write(3,*)'velofl_v',1+ubound(velofl_v)-lbound(velofl_v)
      if(allocated(tkeofl_w))write(3,*)'tkeofl_w',1+ubound(tkeofl_w)-lbound(tkeofl_w)
      if(allocated(bioofl_t))write(3,*)'bioofl_t',1+ubound(bioofl_t)-lbound(bioofl_t)
      if(allocated(dfvofl_w))write(3,*)'dfvofl_w',1+ubound(dfvofl_w)-lbound(dfvofl_w)
      if(allocated(uwindabl_t))write(3,*)'uwindabl_t',1+ubound(uwindabl_t)-lbound(uwindabl_t)
      if(allocated(vwindabl_t))write(3,*)'vwindabl_t',1+ubound(vwindabl_t)-lbound(vwindabl_t)
      if(allocated(w0mofl_w))write(3,*)'w0mofl_w',1+ubound(w0mofl_w)-lbound(w0mofl_w)
      if(allocated(w_keq1_ofl_w))write(3,*)'w_keq1_ofl_w',1+ubound(w_keq1_ofl_w)-lbound(w_keq1_ofl_w)
      if(allocated(kslofl_t))write(3,*)'kslofl_t',1+ubound(kslofl_t)-lbound(kslofl_t)
      if(allocated(ablheight_t))write(3,*)'ablheight_t',1+ubound(ablheight_t)-lbound(ablheight_t)
      if(allocated(wwindabl_w))write(3,*)'wwindabl_w',1+ubound(wwindabl_w)-lbound(wwindabl_w)
      if(allocated(kz_abl_w))write(3,*)'kz_abl_w',1+ubound(kz_abl_w)-lbound(kz_abl_w)
      if(allocated(upwzone0_t))write(3,*)'upwzone0_t',1+ubound(upwzone0_t)-lbound(upwzone0_t)
      if(allocated(velstokes_u))write(3,*)'velstokes_u',1+ubound(velstokes_u)-lbound(velstokes_u)
      if(allocated(velstokes_v))write(3,*)'velstokes_v',1+ubound(velstokes_v)-lbound(velstokes_v)
      if(allocated(velbarstokes_u))write(3,*)'velbarstokes_u',1+ubound(velbarstokes_u)-lbound(velbarstokes_u)
      if(allocated(velbarstokes_v))write(3,*)'velbarstokes_v',1+ubound(velbarstokes_v)-lbound(velbarstokes_v)
      if(allocated(nhp1_t))write(3,*)'nhp1_t',1+ubound(nhp1_t)-lbound(nhp1_t)
      if(allocated(nhp2_t))write(3,*)'nhp2_t',1+ubound(nhp2_t)-lbound(nhp2_t)
      if(allocated(temf_t))write(3,*)'temf_t',1+ubound(temf_t)-lbound(temf_t)
      if(allocated(salf_t))write(3,*)'salf_t',1+ubound(salf_t)-lbound(salf_t)
      if(allocated(temlwf_t))write(3,*)'temlwf_t',1+ubound(temlwf_t)-lbound(temlwf_t)
      if(allocated(sallwf_t))write(3,*)'sallwf_t',1+ubound(sallwf_t)-lbound(sallwf_t)
      if(allocated(sshlwf_w))write(3,*)'sshlwf_w',1+ubound(sshlwf_w)-lbound(sshlwf_w)
      if(allocated(t_wave_t))write(3,*)'t_wave_t',1+ubound(t_wave_t)-lbound(t_wave_t)
      if(allocated(hs_wave_t))write(3,*)'hs_wave_t',1+ubound(hs_wave_t)-lbound(hs_wave_t)
      if(allocated(hsw_wave_t))write(3,*)'hsw_wave_t',1+ubound(hsw_wave_t)-lbound(hsw_wave_t)
      if(allocated(foc_wave_t))write(3,*)'foc_wave_t',1+ubound(foc_wave_t)-lbound(foc_wave_t)
      if(allocated(k_wave_t))write(3,*)'k_wave_t',1+ubound(k_wave_t)-lbound(k_wave_t)
      if(allocated(kx_wave_t))write(3,*)'kx_wave_t',1+ubound(kx_wave_t)-lbound(kx_wave_t)
      if(allocated(ky_wave_t))write(3,*)'ky_wave_t',1+ubound(ky_wave_t)-lbound(ky_wave_t)
      if(allocated(twox_wave_t))write(3,*)'twox_wave_t',1+ubound(twox_wave_t)-lbound(twox_wave_t)
      if(allocated(twoy_wave_t))write(3,*)'twoy_wave_t',1+ubound(twoy_wave_t)-lbound(twoy_wave_t)
      if(allocated(tawx_wave_t))write(3,*)'tawx_wave_t',1+ubound(tawx_wave_t)-lbound(tawx_wave_t)
      if(allocated(tawy_wave_t))write(3,*)'tawy_wave_t',1+ubound(tawy_wave_t)-lbound(tawy_wave_t)
      if(allocated(usf_wave_t))write(3,*)'usf_wave_t',1+ubound(usf_wave_t)-lbound(usf_wave_t)
      if(allocated(vsf_wave_t))write(3,*)'vsf_wave_t',1+ubound(vsf_wave_t)-lbound(vsf_wave_t)
      if(allocated(dir_wave_t))write(3,*)'dir_wave_t',1+ubound(dir_wave_t)-lbound(dir_wave_t)
      if(allocated(uss_wave_t))write(3,*)'uss_wave_t',1+ubound(uss_wave_t)-lbound(uss_wave_t)
      if(allocated(j_wave_t))write(3,*)'j_wave_t',1+ubound(j_wave_t)-lbound(j_wave_t)
      if(allocated(vss_wave_t))write(3,*)'vss_wave_t',1+ubound(vss_wave_t)-lbound(vss_wave_t)
      if(allocated(ubw))write(3,*)'ubw',1+ubound(ubw)-lbound(ubw)
      if(allocated(fw))write(3,*)'fw',1+ubound(fw)-lbound(fw)
      if(allocated(dpt_wave_t))write(3,*)'dpt_wave_t',1+ubound(dpt_wave_t)-lbound(dpt_wave_t)
      if(allocated(wstresb_u))write(3,*)'wstresb_u',1+ubound(wstresb_u)-lbound(wstresb_u)
      if(allocated(wstresb_v))write(3,*)'wstresb_v',1+ubound(wstresb_v)-lbound(wstresb_v)
      if(allocated(ij2ww3_i))write(3,*)'ij2ww3_i',1+ubound(ij2ww3_i)-lbound(ij2ww3_i)
      if(allocated(ij2ww3_j))write(3,*)'ij2ww3_j',1+ubound(ij2ww3_j)-lbound(ij2ww3_j)
      if(allocated(ij2ww3_teta))write(3,*)'ij2ww3_teta',1+ubound(ij2ww3_teta)-lbound(ij2ww3_teta)
      if(allocated(slhf_aver_w))write(3,*)'slhf_aver_w',1+ubound(slhf_aver_w)-lbound(slhf_aver_w)
      if(allocated(sshf_aver_w))write(3,*)'sshf_aver_w',1+ubound(sshf_aver_w)-lbound(sshf_aver_w)
      if(allocated(snsf_aver_w))write(3,*)'snsf_aver_w',1+ubound(snsf_aver_w)-lbound(snsf_aver_w)
      if(allocated(ssr_aver_w))write(3,*)'ssr_aver_w',1+ubound(ssr_aver_w)-lbound(ssr_aver_w)
      if(allocated(precipi_aver_w))write(3,*)'precipi_aver_w',1+ubound(precipi_aver_w)-lbound(precipi_aver_w)
      if(allocated(wstress_aver_u))write(3,*)'wstress_aver_u',1+ubound(wstress_aver_u)-lbound(wstress_aver_u)
      if(allocated(wstress_aver_v))write(3,*)'wstress_aver_v',1+ubound(wstress_aver_v)-lbound(wstress_aver_v)
      if(allocated(hsedofl_t))write(3,*)'hsedofl_t',1+ubound(hsedofl_t)-lbound(hsedofl_t)
      write(3,*)'zone1saltflux_w'
      write(3,*)'zone1tempflux_w'
      write(3,*)'zone1waterflux_w'
      write(3,*)'zone1_nlayer'
      write(3,*)'zone1_max'
      write(3,*)'zone1_u_max'
      write(3,*)'zone1_v_max'
      write(3,*)'zone1_inv_dz'
      write(3,*)'zone1_stretch_dz'
      if(allocated(zone1saltflux_glb))write(3,*)'zone1saltflux_glb',1+ubound(zone1saltflux_glb)-lbound(zone1saltflux_glb)
      if(allocated(zone1saltflux_u))write(3,*)'zone1saltflux_u',1+ubound(zone1saltflux_u)-lbound(zone1saltflux_u)
      if(allocated(zone1saltflux_v))write(3,*)'zone1saltflux_v',1+ubound(zone1saltflux_v)-lbound(zone1saltflux_v)
      if(allocated(zone1tempflux_glb))write(3,*)'zone1tempflux_glb',1+ubound(zone1tempflux_glb)-lbound(zone1tempflux_glb)
      if(allocated(zone1tempflux_u))write(3,*)'zone1tempflux_u',1+ubound(zone1tempflux_u)-lbound(zone1tempflux_u)
      if(allocated(zone1tempflux_v))write(3,*)'zone1tempflux_v',1+ubound(zone1tempflux_v)-lbound(zone1tempflux_v)
      if(allocated(zone1waterflux_glb))write(3,*)'zone1waterflux_glb',1+ubound(zone1waterflux_glb)-lbound(zone1waterflux_glb)
      if(allocated(zone1waterflux_u))write(3,*)'zone1waterflux_u',1+ubound(zone1waterflux_u)-lbound(zone1waterflux_u)
      if(allocated(zone1waterflux_v))write(3,*)'zone1waterflux_v',1+ubound(zone1waterflux_v)-lbound(zone1waterflux_v)
      if(allocated(zone1saltflux_glb_in))write(3,*)'zone1saltflux_glb_in',1+ubound(zone1saltflux_glb_in)-lbound(zone1saltflux_glb_in)
      if(allocated(zone1saltflux_u_in))write(3,*)'zone1saltflux_u_in',1+ubound(zone1saltflux_u_in)-lbound(zone1saltflux_u_in)
      if(allocated(zone1saltflux_v_in))write(3,*)'zone1saltflux_v_in',1+ubound(zone1saltflux_v_in)-lbound(zone1saltflux_v_in)
      if(allocated(zone1tempflux_glb_in))write(3,*)'zone1tempflux_glb_in',1+ubound(zone1tempflux_glb_in)-lbound(zone1tempflux_glb_in)
      if(allocated(zone1tempflux_u_in))write(3,*)'zone1tempflux_u_in',1+ubound(zone1tempflux_u_in)-lbound(zone1tempflux_u_in)
      if(allocated(zone1tempflux_v_in))write(3,*)'zone1tempflux_v_in',1+ubound(zone1tempflux_v_in)-lbound(zone1tempflux_v_in)
      if(allocated(zone1waterflux_glb_in))write(3,*)'zone1waterflux_glb_in',1+ubound(zone1waterflux_glb_in)-lbound(zone1waterflux_glb_in)
      if(allocated(zone1waterflux_u_in))write(3,*)'zone1waterflux_u_in',1+ubound(zone1waterflux_u_in)-lbound(zone1waterflux_u_in)
      if(allocated(zone1waterflux_v_in))write(3,*)'zone1waterflux_v_in',1+ubound(zone1waterflux_v_in)-lbound(zone1waterflux_v_in)
      if(allocated(zone1saltflux_glb_out))write(3,*)'zone1saltflux_glb_out',1+ubound(zone1saltflux_glb_out)-lbound(zone1saltflux_glb_out)
      if(allocated(zone1saltflux_u_out))write(3,*)'zone1saltflux_u_out',1+ubound(zone1saltflux_u_out)-lbound(zone1saltflux_u_out)
      if(allocated(zone1saltflux_v_out))write(3,*)'zone1saltflux_v_out',1+ubound(zone1saltflux_v_out)-lbound(zone1saltflux_v_out)
      if(allocated(zone1tempflux_glb_out))write(3,*)'zone1tempflux_glb_out',1+ubound(zone1tempflux_glb_out)-lbound(zone1tempflux_glb_out)
      if(allocated(zone1tempflux_u_out))write(3,*)'zone1tempflux_u_out',1+ubound(zone1tempflux_u_out)-lbound(zone1tempflux_u_out)
      if(allocated(zone1tempflux_v_out))write(3,*)'zone1tempflux_v_out',1+ubound(zone1tempflux_v_out)-lbound(zone1tempflux_v_out)
      if(allocated(zone1waterflux_glb_out))write(3,*)'zone1waterflux_glb_out',1+ubound(zone1waterflux_glb_out)-lbound(zone1waterflux_glb_out)
      if(allocated(zone1waterflux_u_out))write(3,*)'zone1waterflux_u_out',1+ubound(zone1waterflux_u_out)-lbound(zone1waterflux_u_out)
      if(allocated(zone1waterflux_v_out))write(3,*)'zone1waterflux_v_out',1+ubound(zone1waterflux_v_out)-lbound(zone1waterflux_v_out)
      if(allocated(zone1_mask))write(3,*)'zone1_mask',1+ubound(zone1_mask)-lbound(zone1_mask)
      if(allocated(zone1_flux_u_node))write(3,*)'zone1_flux_u_node',1+ubound(zone1_flux_u_node)-lbound(zone1_flux_u_node)
      if(allocated(zone1_flux_v_node))write(3,*)'zone1_flux_v_node',1+ubound(zone1_flux_v_node)-lbound(zone1_flux_v_node)
      write(3,*)'zone2saltflux_w'
      write(3,*)'zone2tempflux_w'
      write(3,*)'zone2waterflux_w'
      write(3,*)'zone2_nlayer'
      write(3,*)'zone2_max'
      write(3,*)'zone2_u_max'
      write(3,*)'zone2_v_max'
      write(3,*)'zone2_inv_dz'
      write(3,*)'zone2_stretch_dz'
      if(allocated(zone2saltflux_glb))write(3,*)'zone2saltflux_glb',1+ubound(zone2saltflux_glb)-lbound(zone2saltflux_glb)
      if(allocated(zone2saltflux_u))write(3,*)'zone2saltflux_u',1+ubound(zone2saltflux_u)-lbound(zone2saltflux_u)
      if(allocated(zone2saltflux_v))write(3,*)'zone2saltflux_v',1+ubound(zone2saltflux_v)-lbound(zone2saltflux_v)
      if(allocated(zone2tempflux_glb))write(3,*)'zone2tempflux_glb',1+ubound(zone2tempflux_glb)-lbound(zone2tempflux_glb)
      if(allocated(zone2tempflux_u))write(3,*)'zone2tempflux_u',1+ubound(zone2tempflux_u)-lbound(zone2tempflux_u)
      if(allocated(zone2tempflux_v))write(3,*)'zone2tempflux_v',1+ubound(zone2tempflux_v)-lbound(zone2tempflux_v)
      if(allocated(zone2waterflux_glb))write(3,*)'zone2waterflux_glb',1+ubound(zone2waterflux_glb)-lbound(zone2waterflux_glb)
      if(allocated(zone2waterflux_u))write(3,*)'zone2waterflux_u',1+ubound(zone2waterflux_u)-lbound(zone2waterflux_u)
      if(allocated(zone2waterflux_v))write(3,*)'zone2waterflux_v',1+ubound(zone2waterflux_v)-lbound(zone2waterflux_v)
      if(allocated(zone2saltflux_glb_in))write(3,*)'zone2saltflux_glb_in',1+ubound(zone2saltflux_glb_in)-lbound(zone2saltflux_glb_in)
      if(allocated(zone2saltflux_u_in))write(3,*)'zone2saltflux_u_in',1+ubound(zone2saltflux_u_in)-lbound(zone2saltflux_u_in)
      if(allocated(zone2saltflux_v_in))write(3,*)'zone2saltflux_v_in',1+ubound(zone2saltflux_v_in)-lbound(zone2saltflux_v_in)
      if(allocated(zone2tempflux_glb_in))write(3,*)'zone2tempflux_glb_in',1+ubound(zone2tempflux_glb_in)-lbound(zone2tempflux_glb_in)
      if(allocated(zone2tempflux_u_in))write(3,*)'zone2tempflux_u_in',1+ubound(zone2tempflux_u_in)-lbound(zone2tempflux_u_in)
      if(allocated(zone2tempflux_v_in))write(3,*)'zone2tempflux_v_in',1+ubound(zone2tempflux_v_in)-lbound(zone2tempflux_v_in)
      if(allocated(zone2waterflux_glb_in))write(3,*)'zone2waterflux_glb_in',1+ubound(zone2waterflux_glb_in)-lbound(zone2waterflux_glb_in)
      if(allocated(zone2waterflux_u_in))write(3,*)'zone2waterflux_u_in',1+ubound(zone2waterflux_u_in)-lbound(zone2waterflux_u_in)
      if(allocated(zone2waterflux_v_in))write(3,*)'zone2waterflux_v_in',1+ubound(zone2waterflux_v_in)-lbound(zone2waterflux_v_in)
      if(allocated(zone2saltflux_glb_out))write(3,*)'zone2saltflux_glb_out',1+ubound(zone2saltflux_glb_out)-lbound(zone2saltflux_glb_out)
      if(allocated(zone2saltflux_u_out))write(3,*)'zone2saltflux_u_out',1+ubound(zone2saltflux_u_out)-lbound(zone2saltflux_u_out)
      if(allocated(zone2saltflux_v_out))write(3,*)'zone2saltflux_v_out',1+ubound(zone2saltflux_v_out)-lbound(zone2saltflux_v_out)
      if(allocated(zone2tempflux_glb_out))write(3,*)'zone2tempflux_glb_out',1+ubound(zone2tempflux_glb_out)-lbound(zone2tempflux_glb_out)
      if(allocated(zone2tempflux_u_out))write(3,*)'zone2tempflux_u_out',1+ubound(zone2tempflux_u_out)-lbound(zone2tempflux_u_out)
      if(allocated(zone2tempflux_v_out))write(3,*)'zone2tempflux_v_out',1+ubound(zone2tempflux_v_out)-lbound(zone2tempflux_v_out)
      if(allocated(zone2waterflux_glb_out))write(3,*)'zone2waterflux_glb_out',1+ubound(zone2waterflux_glb_out)-lbound(zone2waterflux_glb_out)
      if(allocated(zone2waterflux_u_out))write(3,*)'zone2waterflux_u_out',1+ubound(zone2waterflux_u_out)-lbound(zone2waterflux_u_out)
      if(allocated(zone2waterflux_v_out))write(3,*)'zone2waterflux_v_out',1+ubound(zone2waterflux_v_out)-lbound(zone2waterflux_v_out)
      if(allocated(zone2_mask))write(3,*)'zone2_mask',1+ubound(zone2_mask)-lbound(zone2_mask)
      if(allocated(zone2_flux_u_node))write(3,*)'zone2_flux_u_node',1+ubound(zone2_flux_u_node)-lbound(zone2_flux_u_node)
      if(allocated(zone2_flux_v_node))write(3,*)'zone2_flux_v_node',1+ubound(zone2_flux_v_node)-lbound(zone2_flux_v_node)
      write(3,*)'zone3saltflux_w'
      write(3,*)'zone3tempflux_w'
      write(3,*)'zone3waterflux_w'
      write(3,*)'zone3_nlayer'
      write(3,*)'zone3_max'
      write(3,*)'zone3_u_max'
      write(3,*)'zone3_v_max'
      write(3,*)'zone3_inv_dz'
      write(3,*)'zone3_stretch_dz'
      if(allocated(zone3saltflux_glb))write(3,*)'zone3saltflux_glb',1+ubound(zone3saltflux_glb)-lbound(zone3saltflux_glb)
      if(allocated(zone3saltflux_u))write(3,*)'zone3saltflux_u',1+ubound(zone3saltflux_u)-lbound(zone3saltflux_u)
      if(allocated(zone3saltflux_v))write(3,*)'zone3saltflux_v',1+ubound(zone3saltflux_v)-lbound(zone3saltflux_v)
      if(allocated(zone3tempflux_glb))write(3,*)'zone3tempflux_glb',1+ubound(zone3tempflux_glb)-lbound(zone3tempflux_glb)
      if(allocated(zone3tempflux_u))write(3,*)'zone3tempflux_u',1+ubound(zone3tempflux_u)-lbound(zone3tempflux_u)
      if(allocated(zone3tempflux_v))write(3,*)'zone3tempflux_v',1+ubound(zone3tempflux_v)-lbound(zone3tempflux_v)
      if(allocated(zone3waterflux_glb))write(3,*)'zone3waterflux_glb',1+ubound(zone3waterflux_glb)-lbound(zone3waterflux_glb)
      if(allocated(zone3waterflux_u))write(3,*)'zone3waterflux_u',1+ubound(zone3waterflux_u)-lbound(zone3waterflux_u)
      if(allocated(zone3waterflux_v))write(3,*)'zone3waterflux_v',1+ubound(zone3waterflux_v)-lbound(zone3waterflux_v)
      if(allocated(zone3saltflux_glb_in))write(3,*)'zone3saltflux_glb_in',1+ubound(zone3saltflux_glb_in)-lbound(zone3saltflux_glb_in)
      if(allocated(zone3saltflux_u_in))write(3,*)'zone3saltflux_u_in',1+ubound(zone3saltflux_u_in)-lbound(zone3saltflux_u_in)
      if(allocated(zone3saltflux_v_in))write(3,*)'zone3saltflux_v_in',1+ubound(zone3saltflux_v_in)-lbound(zone3saltflux_v_in)
      if(allocated(zone3tempflux_glb_in))write(3,*)'zone3tempflux_glb_in',1+ubound(zone3tempflux_glb_in)-lbound(zone3tempflux_glb_in)
      if(allocated(zone3tempflux_u_in))write(3,*)'zone3tempflux_u_in',1+ubound(zone3tempflux_u_in)-lbound(zone3tempflux_u_in)
      if(allocated(zone3tempflux_v_in))write(3,*)'zone3tempflux_v_in',1+ubound(zone3tempflux_v_in)-lbound(zone3tempflux_v_in)
      if(allocated(zone3waterflux_glb_in))write(3,*)'zone3waterflux_glb_in',1+ubound(zone3waterflux_glb_in)-lbound(zone3waterflux_glb_in)
      if(allocated(zone3waterflux_u_in))write(3,*)'zone3waterflux_u_in',1+ubound(zone3waterflux_u_in)-lbound(zone3waterflux_u_in)
      if(allocated(zone3waterflux_v_in))write(3,*)'zone3waterflux_v_in',1+ubound(zone3waterflux_v_in)-lbound(zone3waterflux_v_in)
      if(allocated(zone3saltflux_glb_out))write(3,*)'zone3saltflux_glb_out',1+ubound(zone3saltflux_glb_out)-lbound(zone3saltflux_glb_out)
      if(allocated(zone3saltflux_u_out))write(3,*)'zone3saltflux_u_out',1+ubound(zone3saltflux_u_out)-lbound(zone3saltflux_u_out)
      if(allocated(zone3saltflux_v_out))write(3,*)'zone3saltflux_v_out',1+ubound(zone3saltflux_v_out)-lbound(zone3saltflux_v_out)
      if(allocated(zone3tempflux_glb_out))write(3,*)'zone3tempflux_glb_out',1+ubound(zone3tempflux_glb_out)-lbound(zone3tempflux_glb_out)
      if(allocated(zone3tempflux_u_out))write(3,*)'zone3tempflux_u_out',1+ubound(zone3tempflux_u_out)-lbound(zone3tempflux_u_out)
      if(allocated(zone3tempflux_v_out))write(3,*)'zone3tempflux_v_out',1+ubound(zone3tempflux_v_out)-lbound(zone3tempflux_v_out)
      if(allocated(zone3waterflux_glb_out))write(3,*)'zone3waterflux_glb_out',1+ubound(zone3waterflux_glb_out)-lbound(zone3waterflux_glb_out)
      if(allocated(zone3waterflux_u_out))write(3,*)'zone3waterflux_u_out',1+ubound(zone3waterflux_u_out)-lbound(zone3waterflux_u_out)
      if(allocated(zone3waterflux_v_out))write(3,*)'zone3waterflux_v_out',1+ubound(zone3waterflux_v_out)-lbound(zone3waterflux_v_out)
      if(allocated(zone3_mask))write(3,*)'zone3_mask',1+ubound(zone3_mask)-lbound(zone3_mask)
      if(allocated(zone3_flux_u_node))write(3,*)'zone3_flux_u_node',1+ubound(zone3_flux_u_node)-lbound(zone3_flux_u_node)
      if(allocated(zone3_flux_v_node))write(3,*)'zone3_flux_v_node',1+ubound(zone3_flux_v_node)-lbound(zone3_flux_v_node)
      write(3,*)'zone4saltflux_w'
      write(3,*)'zone4tempflux_w'
      write(3,*)'zone4waterflux_w'
      write(3,*)'zone4_nlayer'
      write(3,*)'zone4_max'
      write(3,*)'zone4_u_max'
      write(3,*)'zone4_v_max'
      write(3,*)'zone4_inv_dz'
      write(3,*)'zone4_stretch_dz'
      if(allocated(zone4saltflux_glb))write(3,*)'zone4saltflux_glb',1+ubound(zone4saltflux_glb)-lbound(zone4saltflux_glb)
      if(allocated(zone4saltflux_u))write(3,*)'zone4saltflux_u',1+ubound(zone4saltflux_u)-lbound(zone4saltflux_u)
      if(allocated(zone4saltflux_v))write(3,*)'zone4saltflux_v',1+ubound(zone4saltflux_v)-lbound(zone4saltflux_v)
      if(allocated(zone4tempflux_glb))write(3,*)'zone4tempflux_glb',1+ubound(zone4tempflux_glb)-lbound(zone4tempflux_glb)
      if(allocated(zone4tempflux_u))write(3,*)'zone4tempflux_u',1+ubound(zone4tempflux_u)-lbound(zone4tempflux_u)
      if(allocated(zone4tempflux_v))write(3,*)'zone4tempflux_v',1+ubound(zone4tempflux_v)-lbound(zone4tempflux_v)
      if(allocated(zone4waterflux_glb))write(3,*)'zone4waterflux_glb',1+ubound(zone4waterflux_glb)-lbound(zone4waterflux_glb)
      if(allocated(zone4waterflux_u))write(3,*)'zone4waterflux_u',1+ubound(zone4waterflux_u)-lbound(zone4waterflux_u)
      if(allocated(zone4waterflux_v))write(3,*)'zone4waterflux_v',1+ubound(zone4waterflux_v)-lbound(zone4waterflux_v)
      if(allocated(zone4saltflux_glb_in))write(3,*)'zone4saltflux_glb_in',1+ubound(zone4saltflux_glb_in)-lbound(zone4saltflux_glb_in)
      if(allocated(zone4saltflux_u_in))write(3,*)'zone4saltflux_u_in',1+ubound(zone4saltflux_u_in)-lbound(zone4saltflux_u_in)
      if(allocated(zone4saltflux_v_in))write(3,*)'zone4saltflux_v_in',1+ubound(zone4saltflux_v_in)-lbound(zone4saltflux_v_in)
      if(allocated(zone4tempflux_glb_in))write(3,*)'zone4tempflux_glb_in',1+ubound(zone4tempflux_glb_in)-lbound(zone4tempflux_glb_in)
      if(allocated(zone4tempflux_u_in))write(3,*)'zone4tempflux_u_in',1+ubound(zone4tempflux_u_in)-lbound(zone4tempflux_u_in)
      if(allocated(zone4tempflux_v_in))write(3,*)'zone4tempflux_v_in',1+ubound(zone4tempflux_v_in)-lbound(zone4tempflux_v_in)
      if(allocated(zone4waterflux_glb_in))write(3,*)'zone4waterflux_glb_in',1+ubound(zone4waterflux_glb_in)-lbound(zone4waterflux_glb_in)
      if(allocated(zone4waterflux_u_in))write(3,*)'zone4waterflux_u_in',1+ubound(zone4waterflux_u_in)-lbound(zone4waterflux_u_in)
      if(allocated(zone4waterflux_v_in))write(3,*)'zone4waterflux_v_in',1+ubound(zone4waterflux_v_in)-lbound(zone4waterflux_v_in)
      if(allocated(zone4saltflux_glb_out))write(3,*)'zone4saltflux_glb_out',1+ubound(zone4saltflux_glb_out)-lbound(zone4saltflux_glb_out)
      if(allocated(zone4saltflux_u_out))write(3,*)'zone4saltflux_u_out',1+ubound(zone4saltflux_u_out)-lbound(zone4saltflux_u_out)
      if(allocated(zone4saltflux_v_out))write(3,*)'zone4saltflux_v_out',1+ubound(zone4saltflux_v_out)-lbound(zone4saltflux_v_out)
      if(allocated(zone4tempflux_glb_out))write(3,*)'zone4tempflux_glb_out',1+ubound(zone4tempflux_glb_out)-lbound(zone4tempflux_glb_out)
      if(allocated(zone4tempflux_u_out))write(3,*)'zone4tempflux_u_out',1+ubound(zone4tempflux_u_out)-lbound(zone4tempflux_u_out)
      if(allocated(zone4tempflux_v_out))write(3,*)'zone4tempflux_v_out',1+ubound(zone4tempflux_v_out)-lbound(zone4tempflux_v_out)
      if(allocated(zone4waterflux_glb_out))write(3,*)'zone4waterflux_glb_out',1+ubound(zone4waterflux_glb_out)-lbound(zone4waterflux_glb_out)
      if(allocated(zone4waterflux_u_out))write(3,*)'zone4waterflux_u_out',1+ubound(zone4waterflux_u_out)-lbound(zone4waterflux_u_out)
      if(allocated(zone4waterflux_v_out))write(3,*)'zone4waterflux_v_out',1+ubound(zone4waterflux_v_out)-lbound(zone4waterflux_v_out)
      if(allocated(zone4_mask))write(3,*)'zone4_mask',1+ubound(zone4_mask)-lbound(zone4_mask)
      if(allocated(zone4_flux_u_node))write(3,*)'zone4_flux_u_node',1+ubound(zone4_flux_u_node)-lbound(zone4_flux_u_node)
      if(allocated(zone4_flux_v_node))write(3,*)'zone4_flux_v_node',1+ubound(zone4_flux_v_node)-lbound(zone4_flux_v_node)
      write(3,*)'zone5saltflux_w'
      write(3,*)'zone5tempflux_w'
      write(3,*)'zone5waterflux_w'
      write(3,*)'zone5_nlayer'
      write(3,*)'zone5_max'
      write(3,*)'zone5_u_max'
      write(3,*)'zone5_v_max'
      write(3,*)'zone5_inv_dz'
      write(3,*)'zone5_stretch_dz'
      if(allocated(zone5saltflux_glb))write(3,*)'zone5saltflux_glb',1+ubound(zone5saltflux_glb)-lbound(zone5saltflux_glb)
      if(allocated(zone5saltflux_u))write(3,*)'zone5saltflux_u',1+ubound(zone5saltflux_u)-lbound(zone5saltflux_u)
      if(allocated(zone5saltflux_v))write(3,*)'zone5saltflux_v',1+ubound(zone5saltflux_v)-lbound(zone5saltflux_v)
      if(allocated(zone5tempflux_glb))write(3,*)'zone5tempflux_glb',1+ubound(zone5tempflux_glb)-lbound(zone5tempflux_glb)
      if(allocated(zone5tempflux_u))write(3,*)'zone5tempflux_u',1+ubound(zone5tempflux_u)-lbound(zone5tempflux_u)
      if(allocated(zone5tempflux_v))write(3,*)'zone5tempflux_v',1+ubound(zone5tempflux_v)-lbound(zone5tempflux_v)
      if(allocated(zone5waterflux_glb))write(3,*)'zone5waterflux_glb',1+ubound(zone5waterflux_glb)-lbound(zone5waterflux_glb)
      if(allocated(zone5waterflux_u))write(3,*)'zone5waterflux_u',1+ubound(zone5waterflux_u)-lbound(zone5waterflux_u)
      if(allocated(zone5waterflux_v))write(3,*)'zone5waterflux_v',1+ubound(zone5waterflux_v)-lbound(zone5waterflux_v)
      if(allocated(zone5saltflux_glb_in))write(3,*)'zone5saltflux_glb_in',1+ubound(zone5saltflux_glb_in)-lbound(zone5saltflux_glb_in)
      if(allocated(zone5saltflux_u_in))write(3,*)'zone5saltflux_u_in',1+ubound(zone5saltflux_u_in)-lbound(zone5saltflux_u_in)
      if(allocated(zone5saltflux_v_in))write(3,*)'zone5saltflux_v_in',1+ubound(zone5saltflux_v_in)-lbound(zone5saltflux_v_in)
      if(allocated(zone5tempflux_glb_in))write(3,*)'zone5tempflux_glb_in',1+ubound(zone5tempflux_glb_in)-lbound(zone5tempflux_glb_in)
      if(allocated(zone5tempflux_u_in))write(3,*)'zone5tempflux_u_in',1+ubound(zone5tempflux_u_in)-lbound(zone5tempflux_u_in)
      if(allocated(zone5tempflux_v_in))write(3,*)'zone5tempflux_v_in',1+ubound(zone5tempflux_v_in)-lbound(zone5tempflux_v_in)
      if(allocated(zone5waterflux_glb_in))write(3,*)'zone5waterflux_glb_in',1+ubound(zone5waterflux_glb_in)-lbound(zone5waterflux_glb_in)
      if(allocated(zone5waterflux_u_in))write(3,*)'zone5waterflux_u_in',1+ubound(zone5waterflux_u_in)-lbound(zone5waterflux_u_in)
      if(allocated(zone5waterflux_v_in))write(3,*)'zone5waterflux_v_in',1+ubound(zone5waterflux_v_in)-lbound(zone5waterflux_v_in)
      if(allocated(zone5saltflux_glb_out))write(3,*)'zone5saltflux_glb_out',1+ubound(zone5saltflux_glb_out)-lbound(zone5saltflux_glb_out)
      if(allocated(zone5saltflux_u_out))write(3,*)'zone5saltflux_u_out',1+ubound(zone5saltflux_u_out)-lbound(zone5saltflux_u_out)
      if(allocated(zone5saltflux_v_out))write(3,*)'zone5saltflux_v_out',1+ubound(zone5saltflux_v_out)-lbound(zone5saltflux_v_out)
      if(allocated(zone5tempflux_glb_out))write(3,*)'zone5tempflux_glb_out',1+ubound(zone5tempflux_glb_out)-lbound(zone5tempflux_glb_out)
      if(allocated(zone5tempflux_u_out))write(3,*)'zone5tempflux_u_out',1+ubound(zone5tempflux_u_out)-lbound(zone5tempflux_u_out)
      if(allocated(zone5tempflux_v_out))write(3,*)'zone5tempflux_v_out',1+ubound(zone5tempflux_v_out)-lbound(zone5tempflux_v_out)
      if(allocated(zone5waterflux_glb_out))write(3,*)'zone5waterflux_glb_out',1+ubound(zone5waterflux_glb_out)-lbound(zone5waterflux_glb_out)
      if(allocated(zone5waterflux_u_out))write(3,*)'zone5waterflux_u_out',1+ubound(zone5waterflux_u_out)-lbound(zone5waterflux_u_out)
      if(allocated(zone5waterflux_v_out))write(3,*)'zone5waterflux_v_out',1+ubound(zone5waterflux_v_out)-lbound(zone5waterflux_v_out)
      if(allocated(zone5_mask))write(3,*)'zone5_mask',1+ubound(zone5_mask)-lbound(zone5_mask)
      if(allocated(zone5_flux_u_node))write(3,*)'zone5_flux_u_node',1+ubound(zone5_flux_u_node)-lbound(zone5_flux_u_node)
      if(allocated(zone5_flux_v_node))write(3,*)'zone5_flux_v_node',1+ubound(zone5_flux_v_node)-lbound(zone5_flux_v_node)
      write(3,*)'zone6saltflux_w'
      write(3,*)'zone6tempflux_w'
      write(3,*)'zone6waterflux_w'
      write(3,*)'zone6_nlayer'
      write(3,*)'zone6_max'
      write(3,*)'zone6_u_max'
      write(3,*)'zone6_v_max'
      write(3,*)'zone6_inv_dz'
      write(3,*)'zone6_stretch_dz'
      if(allocated(zone6saltflux_glb))write(3,*)'zone6saltflux_glb',1+ubound(zone6saltflux_glb)-lbound(zone6saltflux_glb)
      if(allocated(zone6saltflux_u))write(3,*)'zone6saltflux_u',1+ubound(zone6saltflux_u)-lbound(zone6saltflux_u)
      if(allocated(zone6saltflux_v))write(3,*)'zone6saltflux_v',1+ubound(zone6saltflux_v)-lbound(zone6saltflux_v)
      if(allocated(zone6tempflux_glb))write(3,*)'zone6tempflux_glb',1+ubound(zone6tempflux_glb)-lbound(zone6tempflux_glb)
      if(allocated(zone6tempflux_u))write(3,*)'zone6tempflux_u',1+ubound(zone6tempflux_u)-lbound(zone6tempflux_u)
      if(allocated(zone6tempflux_v))write(3,*)'zone6tempflux_v',1+ubound(zone6tempflux_v)-lbound(zone6tempflux_v)
      if(allocated(zone6waterflux_glb))write(3,*)'zone6waterflux_glb',1+ubound(zone6waterflux_glb)-lbound(zone6waterflux_glb)
      if(allocated(zone6waterflux_u))write(3,*)'zone6waterflux_u',1+ubound(zone6waterflux_u)-lbound(zone6waterflux_u)
      if(allocated(zone6waterflux_v))write(3,*)'zone6waterflux_v',1+ubound(zone6waterflux_v)-lbound(zone6waterflux_v)
      if(allocated(zone6saltflux_glb_in))write(3,*)'zone6saltflux_glb_in',1+ubound(zone6saltflux_glb_in)-lbound(zone6saltflux_glb_in)
      if(allocated(zone6saltflux_u_in))write(3,*)'zone6saltflux_u_in',1+ubound(zone6saltflux_u_in)-lbound(zone6saltflux_u_in)
      if(allocated(zone6saltflux_v_in))write(3,*)'zone6saltflux_v_in',1+ubound(zone6saltflux_v_in)-lbound(zone6saltflux_v_in)
      if(allocated(zone6tempflux_glb_in))write(3,*)'zone6tempflux_glb_in',1+ubound(zone6tempflux_glb_in)-lbound(zone6tempflux_glb_in)
      if(allocated(zone6tempflux_u_in))write(3,*)'zone6tempflux_u_in',1+ubound(zone6tempflux_u_in)-lbound(zone6tempflux_u_in)
      if(allocated(zone6tempflux_v_in))write(3,*)'zone6tempflux_v_in',1+ubound(zone6tempflux_v_in)-lbound(zone6tempflux_v_in)
      if(allocated(zone6waterflux_glb_in))write(3,*)'zone6waterflux_glb_in',1+ubound(zone6waterflux_glb_in)-lbound(zone6waterflux_glb_in)
      if(allocated(zone6waterflux_u_in))write(3,*)'zone6waterflux_u_in',1+ubound(zone6waterflux_u_in)-lbound(zone6waterflux_u_in)
      if(allocated(zone6waterflux_v_in))write(3,*)'zone6waterflux_v_in',1+ubound(zone6waterflux_v_in)-lbound(zone6waterflux_v_in)
      if(allocated(zone6saltflux_glb_out))write(3,*)'zone6saltflux_glb_out',1+ubound(zone6saltflux_glb_out)-lbound(zone6saltflux_glb_out)
      if(allocated(zone6saltflux_u_out))write(3,*)'zone6saltflux_u_out',1+ubound(zone6saltflux_u_out)-lbound(zone6saltflux_u_out)
      if(allocated(zone6saltflux_v_out))write(3,*)'zone6saltflux_v_out',1+ubound(zone6saltflux_v_out)-lbound(zone6saltflux_v_out)
      if(allocated(zone6tempflux_glb_out))write(3,*)'zone6tempflux_glb_out',1+ubound(zone6tempflux_glb_out)-lbound(zone6tempflux_glb_out)
      if(allocated(zone6tempflux_u_out))write(3,*)'zone6tempflux_u_out',1+ubound(zone6tempflux_u_out)-lbound(zone6tempflux_u_out)
      if(allocated(zone6tempflux_v_out))write(3,*)'zone6tempflux_v_out',1+ubound(zone6tempflux_v_out)-lbound(zone6tempflux_v_out)
      if(allocated(zone6waterflux_glb_out))write(3,*)'zone6waterflux_glb_out',1+ubound(zone6waterflux_glb_out)-lbound(zone6waterflux_glb_out)
      if(allocated(zone6waterflux_u_out))write(3,*)'zone6waterflux_u_out',1+ubound(zone6waterflux_u_out)-lbound(zone6waterflux_u_out)
      if(allocated(zone6waterflux_v_out))write(3,*)'zone6waterflux_v_out',1+ubound(zone6waterflux_v_out)-lbound(zone6waterflux_v_out)
      if(allocated(zone6_mask))write(3,*)'zone6_mask',1+ubound(zone6_mask)-lbound(zone6_mask)
      if(allocated(zone6_flux_u_node))write(3,*)'zone6_flux_u_node',1+ubound(zone6_flux_u_node)-lbound(zone6_flux_u_node)
      if(allocated(zone6_flux_v_node))write(3,*)'zone6_flux_v_node',1+ubound(zone6_flux_v_node)-lbound(zone6_flux_v_node)
      if(allocated(anyvar3d))write(3,*)'anyvar3d',1+ubound(anyvar3d)-lbound(anyvar3d)
      if(allocated(sigma_fric_wu))write(3,*)'sigma_fric_wu',1+ubound(sigma_fric_wu)-lbound(sigma_fric_wu)
      if(allocated(sigma_fric_wv))write(3,*)'sigma_fric_wv',1+ubound(sigma_fric_wv)-lbound(sigma_fric_wv)
      if(allocated(dsig_t))write(3,*)'dsig_t',1+ubound(dsig_t)-lbound(dsig_t)
      if(allocated(anyv3dint))write(3,*)'anyv3dint',1+ubound(anyv3dint)-lbound(anyv3dint)
      if(allocated(sshobc_w))write(3,*)'sshobc_w',1+ubound(sshobc_w)-lbound(sshobc_w)
      if(allocated(velbarobc_u))write(3,*)'velbarobc_u',1+ubound(velbarobc_u)-lbound(velbarobc_u)
      if(allocated(velbarobc_v))write(3,*)'velbarobc_v',1+ubound(velbarobc_v)-lbound(velbarobc_v)
      if(allocated(sshofl_w))write(3,*)'sshofl_w',1+ubound(sshofl_w)-lbound(sshofl_w)
      if(allocated(velbarofl_u))write(3,*)'velbarofl_u',1+ubound(velbarofl_u)-lbound(velbarofl_u)
      if(allocated(velbarofl_v))write(3,*)'velbarofl_v',1+ubound(velbarofl_v)-lbound(velbarofl_v)
      if(allocated(tem_delta_t))write(3,*)'tem_delta_t',1+ubound(tem_delta_t)-lbound(tem_delta_t)
      if(allocated(sal_delta_t))write(3,*)'sal_delta_t',1+ubound(sal_delta_t)-lbound(sal_delta_t)
      if(allocated(anyvar2d))write(3,*)'anyvar2d',1+ubound(anyvar2d)-lbound(anyvar2d)
      if(allocated(iriver))write(3,*)'iriver',1+ubound(iriver)-lbound(iriver)
      if(allocated(jriver))write(3,*)'jriver',1+ubound(jriver)-lbound(jriver)
      if(allocated(rankcoords))write(3,*)'rankcoords',1+ubound(rankcoords)-lbound(rankcoords)
      if(allocated(l2ij_out_u))write(3,*)'l2ij_out_u',1+ubound(l2ij_out_u)-lbound(l2ij_out_u)
      if(allocated(l2ij_out_v))write(3,*)'l2ij_out_v',1+ubound(l2ij_out_v)-lbound(l2ij_out_v)
      if(allocated(l2ij_out_w))write(3,*)'l2ij_out_w',1+ubound(l2ij_out_w)-lbound(l2ij_out_w)
      if(allocated(l2ij_in_u))write(3,*)'l2ij_in_u',1+ubound(l2ij_in_u)-lbound(l2ij_in_u)
      if(allocated(l2ij_in_v))write(3,*)'l2ij_in_v',1+ubound(l2ij_in_v)-lbound(l2ij_in_v)
      if(allocated(l2ij_in_w))write(3,*)'l2ij_in_w',1+ubound(l2ij_in_w)-lbound(l2ij_in_w)
      if(allocated(fillmask_t))write(3,*)'fillmask_t',1+ubound(fillmask_t)-lbound(fillmask_t)
      if(allocated(datesim))write(3,*)'datesim',1+ubound(datesim)-lbound(datesim)
      if(allocated(dateobc))write(3,*)'dateobc',1+ubound(dateobc)-lbound(dateobc)
      if(allocated(dateairsea))write(3,*)'dateairsea',1+ubound(dateairsea)-lbound(dateairsea)
      if(allocated(dateriver))write(3,*)'dateriver',1+ubound(dateriver)-lbound(dateriver)
      if(allocated(runoff_w))write(3,*)'runoff_w',1+ubound(runoff_w)-lbound(runoff_w)
      write(3,*)'drifter_out_sampling'
      write(3,*)'obc2dtype'
      write(3,*)'coef_diss_mangrove'
      write(3,*)'expnum'
      write(3,*)'cst_c0cub'
      write(3,*)'relativewind'
      write(3,*)'offset_sshobc'
      write(3,*)'biharm_2dfactor'
      write(3,*)'checkr0'
      write(3,*)'checkr1'
      write(3,*)'checkr2'
      write(3,*)'checkr3'
      write(3,*)'checkr4'
      write(3,*)'sponge_l'
      write(3,*)'sponge_dx_critic'
      write(3,*)'sponge_dx_width'
      write(3,*)'dzsurfmin'
      if(allocated(drifter_send_order_canal))write(3,*)'drifter_send_order_canal',1+ubound(drifter_send_order_canal)-lbound(drifter_send_order_canal)
      write(3,*)'drifter_onoff'
      if(allocated(alongresolriver))write(3,*)'alongresolriver',1+ubound(alongresolriver)-lbound(alongresolriver)
      if(allocated(crossresolriver))write(3,*)'crossresolriver',1+ubound(crossresolriver)-lbound(crossresolriver)
      write(3,*)'ww3_varmax'
      write(3,*)'ww3_type_grid'
      write(3,*)'type_unstructured'
      write(3,*)'type_structured'
      write(3,*)'vststep'
      write(3,*)'zprofile2d3d'
      write(3,*)'timestep_type'
      write(3,*)'timestep_leapfrog'
      write(3,*)'timestep_forwbckw'
      write(3,*)'bulk_core'
      write(3,*)'bulk_moon'
      write(3,*)'bulk_coare'
      write(3,*)'bulk_ecume'
      write(3,*)'bulk_scheme'
      write(3,*)'tideana_spinup'
      write(3,*)'tideana_delta'
      write(3,*)'ogcm_time_lag'
      write(3,*)'albedo_val'
      write(3,*)'tfilterfb'
      write(3,*)'tfilterlf'
      write(3,*)'ktide'
      write(3,*)'tideforces'
      write(3,*)'tideana_yesno'
      write(3,*)'ioffline'
      write(3,*)'tideanalysis_count'
      write(3,*)'ibl1_advbio'
      write(3,*)'ibl2_advbio'
      write(3,*)'jbl1_advbio'
      write(3,*)'jbl2_advbio'
      write(3,*)'ogcm_time_shift'
      write(3,*)'ieq1'
      write(3,*)'ieqimax'
      write(3,*)'jeq1'
      write(3,*)'jeqjmax'
      write(3,*)'ieq1_jeq1'
      write(3,*)'ieq1_jeqjmax'
      write(3,*)'ieqimax_jeq1'
      write(3,*)'ieqimax_jeqjmax'
      write(3,*)'dim_varid'
      write(3,*)'albedo_constant'
      write(3,*)'albedo_apel1987'
      write(3,*)'albedo_br1982'
      write(3,*)'albedo_br1986'
      write(3,*)'albedo_case'
      write(3,*)'il_oasis_time'
      write(3,*)'flag_steric_effect'
      write(3,*)'flag_remove_secondary_bassins'
      write(3,*)'one_kind1'
      write(3,*)'signe'
      write(3,*)'flag_maxbotstress'
      write(3,*)'flag_offline_binary'
      write(3,*)'flag_wstressbulk'
      write(3,*)'flag_write'
      write(3,*)'mangrove_scheme'
      write(3,*)'flag_meteo_land_plug'
      write(3,*)'flag_meteo_land_plug_wind'
      write(3,*)'flag_upwind_obc'
      write(3,*)'discard_lonlat_periodicity'
      write(3,*)'fplan2_grid'
      write(3,*)'flag_ts_effectivedensity'
      write(3,*)'flag_ts_quicklim'
      write(3,*)'flag_z2dv_outputs'
      write(3,*)'flag_ogcmname2date'
      write(3,*)'flag_meteoname2date'
      write(3,*)'ofl_bio'
      write(3,*)'drifter_output_files'
      write(3,*)'flag_tide3d_analysis'
      write(3,*)'flag_bathy_update'
      write(3,*)'loop1'
      write(3,*)'loop2'
      write(3,*)'loop3'
      write(3,*)'loopmaxtke'
      write(3,*)'loopmaxbio'
      write(3,*)'loopmaxts'
      write(3,*)'nairsea'
      write(3,*)'bi_onoff'
      write(3,*)'nc_or_bin_airsea'
      write(3,*)'loop_netcdf'
      write(3,*)'count_netcdfvar'
      write(3,*)'wavefile_prvtrec'
      write(3,*)'wavefile_nextrec'
      write(3,*)'wave_cpl_nextrec'
      write(3,*)'wave_cpl_period_iter'
      write(3,*)'wave_cpl_ww3_sdir'
      write(3,*)'wave_cpl_ww3_cdir'
      write(3,*)'wave_cpl_ww3_hs'
      write(3,*)'wave_cpl_ww3_hsw'
      write(3,*)'wave_cpl_ww3_foc'
      write(3,*)'wave_cpl_ww3_tw'
      write(3,*)'wave_cpl_ww3_tawx'
      write(3,*)'wave_cpl_ww3_tawy'
      write(3,*)'wave_cpl_ww3_twox'
      write(3,*)'wave_cpl_ww3_twoy'
      write(3,*)'wave_cpl_ww3_uss'
      write(3,*)'wave_cpl_ww3_vss'
      write(3,*)'wave_cpl_ww3_msk'
      write(3,*)'ofl_rec_now'
      write(3,*)'ofl_rec_max'
      write(3,*)'flag_kz_enhanced'
      write(3,*)'flag_net_ir'
      write(3,*)'flag_dt_adjust'
      write(3,*)'tide_interpolation'
      write(3,*)'tide_flagrotation'
      write(3,*)'flag_ssr24avr'
      write(3,*)'flag_abl'
      write(3,*)'flag_abl2'
      write(3,*)'flag_sequoia'
      write(3,*)'dimssr24prv'
      write(3,*)'flag_meteo_average'
      write(3,*)'flag_nemoffline'
      write(3,*)'trc_id'
      write(3,*)'vel_id'
      write(3,*)'ssh_id'
      write(3,*)'var_num'
      write(3,*)'flag_p0m_filter'
      write(3,*)'flag_refstate'
      write(3,*)'flag_linearfric'
      write(3,*)'linear_coef_mangrove'
      write(3,*)'flag_1dv'
      write(3,*)'flag_merged_levels'
      write(3,*)'flag1_smooth_h_mask'
      write(3,*)'flag3_smooth_h_mask'
      write(3,*)'flag_0status_option'
      write(3,*)'flag_rmnegval'
      write(3,*)'timemax'
      write(3,*)'dirmax'
      write(3,*)'freqmax'
      write(3,*)'var_misval'
      write(3,*)'un_r8'
      write(3,*)'heure'
      write(3,*)'suma'
      write(3,*)'sumb'
      write(3,*)'sum_mi'
      write(3,*)'grid_area'
      write(3,*)'grid_areaglb'
      write(3,*)'grid_volumeglb'
      write(3,*)'sumarchive'
      write(3,*)'lonmin'
      write(3,*)'lonmax'
      write(3,*)'latmin'
      write(3,*)'latmax'
      write(3,*)'wavefile_prvtime'
      write(3,*)'wavefile_nextime'
      write(3,*)'ofl_period_prev'
      write(3,*)'ofl_period_now'
      write(3,*)'ofl_period_next'
      write(3,*)'ofl_nextrec_time'
      write(3,*)'ofl_writime'
      write(3,*)'ofl_readtime_next'
      write(3,*)'ofl_readtime_prev'
      write(3,*)'biobc_nextfiletime'
      write(3,*)'biobc_prevfiletime'
      write(3,*)'stability_index'
      write(3,*)'iteration2d_max_r8'
      write(3,*)'iteration2d_upbound'
      write(3,*)'coef_linearfric'
      write(3,*)'gh'
      write(3,*)'gm'
      write(3,*)'nn'
      write(3,*)'tke2overeps'
      write(3,*)'tkeovereps'
      write(3,*)'gravoverrho'
      write(3,*)'sh'
      write(3,*)'sm'
      write(3,*)'heatfluxbias'
      write(3,*)'check0'
      write(3,*)'check1'
      write(3,*)'check2'
      write(3,*)'check3'
      write(3,*)'check4'
      write(3,*)'check5'
      write(3,*)'check6'
      write(3,*)'check7'
      write(3,*)'check8'
      write(3,*)'zero'
      write(3,*)'un'
      write(3,*)'cdseuil'
      write(3,*)'cdb_2dh'
      write(3,*)'small1'
      write(3,*)'deci'
      write(3,*)'decj'
      write(3,*)'deck'
      write(3,*)'rap1'
      write(3,*)'rap2'
      write(3,*)'rapi'
      write(3,*)'rapj'
      write(3,*)'rapk'
      write(3,*)'rap'
      write(3,*)'const0'
      write(3,*)'const1'
      write(3,*)'const2'
      write(3,*)'const3'
      write(3,*)'const4'
      write(3,*)'const5'
      write(3,*)'const6'
      write(3,*)'const7'
      write(3,*)'const8'
      write(3,*)'const9'
      write(3,*)'uv_10'
      write(3,*)'karman'
      write(3,*)'stefan'
      write(3,*)'z2m'
      write(3,*)'z10m'
      write(3,*)'pss0'
      write(3,*)'boltz'
      write(3,*)'planck'
      write(3,*)'avogadro'
      write(3,*)'ce'
      write(3,*)'cen'
      write(3,*)'ch'
      write(3,*)'chn'
      write(3,*)'cd'
      write(3,*)'cdn'
      write(3,*)'z_2'
      write(3,*)'z_10'
      write(3,*)'q_0'
      write(3,*)'r_0'
      write(3,*)'pvs_0'
      write(3,*)'psih_10'
      write(3,*)'psih_2'
      write(3,*)'phim'
      write(3,*)'zl_10'
      write(3,*)'zl_2'
      write(3,*)'falpha'
      write(3,*)'fbeta'
      write(3,*)'ro'
      write(3,*)'cp_air'
      write(3,*)'lv'
      write(3,*)'psim_10'
      write(3,*)'psim_2'
      write(3,*)'dte_lp'
      write(3,*)'inv_dte_lp'
      write(3,*)'inv_dti_lp'
      write(3,*)'inv_dti_fw'
      write(3,*)'inv_dti_fw_p2'
      write(3,*)'dte_fw'
      write(3,*)'dti_lpbef'
      write(3,*)'dti_lp'
      write(3,*)'dti_lpmax'
      write(3,*)'dti_lpsub'
      write(3,*)'dti_fwsub'
      write(3,*)'dti_fwsubio'
      write(3,*)'dti_fw'
      write(3,*)'dti_bef'
      write(3,*)'dti_now'
      write(3,*)'dtiratio'
      write(3,*)'dt_drf'
      write(3,*)'fbtfiltercoef'
      write(3,*)'assel0'
      write(3,*)'assel1'
      write(3,*)'assel2'
      write(3,*)'assel3'
      write(3,*)'wetdry_cst1'
      write(3,*)'wetdry_cst2'
      write(3,*)'wetdry_cst3'
      write(3,*)'h_inf'
      write(3,*)'h_inf_obc'
      write(3,*)'h_sup'
      write(3,*)'dist'
      write(3,*)'dist0'
      write(3,*)'dist1'
      write(3,*)'dist2'
      write(3,*)'dist3'
      write(3,*)'t0surf'
      write(3,*)'s0surf'
      write(3,*)'deg2rad'
      write(3,*)'rad2deg'
      write(3,*)'zmin'
      write(3,*)'zmax'
      write(3,*)'coa'
      write(3,*)'cob'
      write(3,*)'coc'
      write(3,*)'sol1'
      write(3,*)'sol2'
      write(3,*)'discri'
      write(3,*)'profm1'
      write(3,*)'profp1'
      write(3,*)'hmax'
      write(3,*)'hstepmin'
      write(3,*)'hstepmax'
      write(3,*)'grav'
      write(3,*)'cfl_sshmax'
      write(3,*)'cfl_hsshmax'
      write(3,*)'cfl_umax'
      write(3,*)'cfl_reduce'
      write(3,*)'relax_es'
      write(3,*)'relax_ext'
      write(3,*)'relax_int'
      write(3,*)'relax_ts'
      write(3,*)'relax_bpc'
      write(3,*)'momentum_input_depth'
      write(3,*)'z0s'
      write(3,*)'z1'
      write(3,*)'z2'
      write(3,*)'z3'
      write(3,*)'z4'
      write(3,*)'z5'
      write(3,*)'y0'
      write(3,*)'y2'
      write(3,*)'y3'
      write(3,*)'y4'
      write(3,*)'y5'
      write(3,*)'x0'
      write(3,*)'x1'
      write(3,*)'x2'
      write(3,*)'x3'
      write(3,*)'x4'
      write(3,*)'x5'
      write(3,*)'x6'
      write(3,*)'x7'
      write(3,*)'x8'
      write(3,*)'x9'
      write(3,*)'x10'
      write(3,*)'x11'
      write(3,*)'x12'
      write(3,*)'x13'
      write(3,*)'x14'
      write(3,*)'x20'
      write(3,*)'x21'
      write(3,*)'x22'
      write(3,*)'x33'
      write(3,*)'x44'
      write(3,*)'area'
      write(3,*)'tem_validmin'
      write(3,*)'tem_validmax'
      write(3,*)'sal_validmin'
      write(3,*)'sal_validmax'
      write(3,*)'xmd'
      write(3,*)'xmv'
      write(3,*)'xrd'
      write(3,*)'xrv'
      write(3,*)'xcpd'
      write(3,*)'xcpv'
      write(3,*)'xcl'
      write(3,*)'celsius2kelvin'
      write(3,*)'xlvtt'
      write(3,*)'xestt'
      write(3,*)'xgamw'
      write(3,*)'xbetaw'
      write(3,*)'xalpw'
      write(3,*)'zrvsrdm1'
      write(3,*)'z0_u'
      write(3,*)'qsat_sea_z'
      write(3,*)'airdensity'
      write(3,*)'sst_kelvin'
      write(3,*)'prs_atm_z'
      write(3,*)'tem_atm_z'
      write(3,*)'exner_atm_z'
      write(3,*)'delta_u'
      write(3,*)'delta_t'
      write(3,*)'delta_q'
      write(3,*)'delta_u_n'
      write(3,*)'exner_sea_z'
      write(3,*)'psifunctt'
      write(3,*)'psifunctu'
      write(3,*)'z0_q'
      write(3,*)'z0_t'
      write(3,*)'sst1000hpa_kelvin'
      write(3,*)'visa'
      write(3,*)'charnock'
      write(3,*)'qsat_atm_z'
      write(3,*)'ustar_bef'
      write(3,*)'qstar_bef'
      write(3,*)'tetastar_bef'
      write(3,*)'rayonterre'
      write(3,*)'northpole_lon'
      write(3,*)'northpole_lat'
      write(3,*)'southpole_lon'
      write(3,*)'southpole_lat'
      write(3,*)'phi0'
      write(3,*)'longi'
      write(3,*)'latit'
      write(3,*)'longi0'
      write(3,*)'latit0'
      write(3,*)'longi1'
      write(3,*)'latit1'
      write(3,*)'angle0'
      write(3,*)'alp_t'
      write(3,*)'alp_s'
      write(3,*)'pi'
      write(3,*)'rho'
      write(3,*)'inv_rho'
      write(3,*)'rhoair'
      write(3,*)'valmax'
      write(3,*)'vis'
      write(3,*)'tkee1'
      write(3,*)'tkee2'
      write(3,*)'tkee3'
      write(3,*)'tkeg2'
      write(3,*)'tkeg3'
      write(3,*)'tkeg4'
      write(3,*)'tkeg5'
      write(3,*)'tkeg6'
      write(3,*)'tkeb1'
      write(3,*)'tkeb2'
      write(3,*)'ctke1'
      write(3,*)'ctke2'
      write(3,*)'t0'
      write(3,*)'s0'
      write(3,*)'cp'
      write(3,*)'light_kpar1'
      write(3,*)'light_att1'
      write(3,*)'light_att2'
      write(3,*)'light_rat1'
      write(3,*)'light_rat2'
      write(3,*)'light_att2_val1'
      write(3,*)'light_att2_h1'
      write(3,*)'light_att2_val2'
      write(3,*)'light_att2_h2'
      write(3,*)'t0_base'
      write(3,*)'s0_base'
      write(3,*)'rho_base'
      write(3,*)'alp_t_base'
      write(3,*)'alp_s_base'
      write(3,*)'tfb0'
      write(3,*)'meteo_lonmin'
      write(3,*)'meteo_latmin'
      write(3,*)'meteo_lonmax'
      write(3,*)'meteo_latmax'
      write(3,*)'meteo_resol'
      write(3,*)'meteo_resol_u'
      write(3,*)'meteo_resol_v'
      write(3,*)'meteo_lonstr'
      write(3,*)'meteo_lonend'
      write(3,*)'meteo_londlt'
      write(3,*)'meteo_latstr'
      write(3,*)'meteo_latend'
      write(3,*)'meteo_latdlt'
      write(3,*)'var_lonmin'
      write(3,*)'var_latmin'
      write(3,*)'var_lonmax'
      write(3,*)'var_latmax'
      write(3,*)'ww3_lonmin'
      write(3,*)'ww3_latmin'
      write(3,*)'ww3_lonmax'
      write(3,*)'ww3_latmax'
      write(3,*)'ww3_dlon'
      write(3,*)'ww3_dlat'
      write(3,*)'tide_lonmin'
      write(3,*)'tide_latmin'
      write(3,*)'tide_lonmax'
      write(3,*)'tide_latmax'
      write(3,*)'tide_dlon'
      write(3,*)'tide_dlat'
      write(3,*)'dxb'
      write(3,*)'dyb'
      write(3,*)'dxa'
      write(3,*)'dya'
      write(3,*)'hmin'
      write(3,*)'h1d'
      write(3,*)'dlon'
      write(3,*)'dlat'
      write(3,*)'epsi'
      write(3,*)'lagrange_ssh'
      write(3,*)'diffu'
      write(3,*)'difnorm'
      write(3,*)'lup'
      write(3,*)'ldown'
      write(3,*)'zup'
      write(3,*)'zdown'
      write(3,*)'rbase'
      write(3,*)'small'
      write(3,*)'small3'
      write(3,*)'hgesig'
      write(3,*)'pgesig'
      write(3,*)'windfactor'
      write(3,*)'xdtk_out'
      write(3,*)'tfond'
      write(3,*)'sfond'
      write(3,*)'rfond'
      write(3,*)'c1streamf'
      write(3,*)'c2streamf'
      write(3,*)'rampe'
      write(3,*)'rampe_wind'
      write(3,*)'y1'
      write(3,*)'obc_hf_reset'
      write(3,*)'rho_0d'
      write(3,*)'tem_0d'
      write(3,*)'sal_0d'
      write(3,*)'rho_tmp'
      write(3,*)'cst_adv_hor'
      write(3,*)'cst_adv_ver'
      write(3,*)'cst_adv_vel'
      write(3,*)'ssh_avr_nest_out'
      write(3,*)'graph_nextime'
      write(3,*)'tidenodal_prev_rdv'
      write(3,*)'tidenodal_next_rdv'
      write(3,*)'tideana_modulo'
      write(3,*)'tideana_nextime'
      write(3,*)'cellboxfactor1'
      write(3,*)'cellboxfactor2'
      write(3,*)'rap_wave'
      write(3,*)'ihmax'
      write(3,*)'jhmax'
      write(3,*)'ihmin'
      write(3,*)'jhmin'
      write(3,*)'convect_yn'
      write(3,*)'nbvstepmin'
      write(3,*)'bio_relax_size'
      write(3,*)'unit_r4'
      write(3,*)'x0_r4'
      write(3,*)'x1_r4'
      write(3,*)'x2_r4'
      write(3,*)'x3_r4'
      write(3,*)'x4_r4'
      write(3,*)'time_r4'
      write(3,*)'discharge'
      write(3,*)'filval'
      write(3,*)'var_validmin'
      write(3,*)'var_validmax'
      write(3,*)'vdw_loc'
      write(3,*)'vup_loc'
      write(3,*)'zdw_loc'
      write(3,*)'var_scalefactor'
      write(3,*)'inv_scalefactor'
      write(3,*)'var_addoffset'
      write(3,*)'zup_loc'
      write(3,*)'hrmax'
      write(3,*)'relax_bio'
      write(3,*)'kmol_m'
      write(3,*)'kmol_h'
      write(3,*)'kmol_s'
      write(3,*)'upw_hrange1'
      write(3,*)'upw_hrange2'
      write(3,*)'inv_ekman_depth'
      write(3,*)'constant_km'
      write(3,*)'constant_kh'
      write(3,*)'z0b'
      write(3,*)'z0b_land'
      write(3,*)'z0b_rivers'
      write(3,*)'zlevel_land'
      write(3,*)'invloopmaxts'
      write(3,*)'invloopmaxu'
      write(3,*)'invloopmaxv'
      write(3,*)'sqrtgrav'
      write(3,*)'invgrav'
      write(3,*)'spinup_forcing'
      write(3,*)'relax_lwf'
      write(3,*)'dz_vertical_incr_fact'
      write(3,*)'ratio_negdif_ver'
      write(3,*)'ratio_negdif_hor'
      write(3,*)'ratio_bionegdif'
      write(3,*)'vqs_cst1'
      write(3,*)'vqs_cst2'
      write(3,*)'vqs_cst3'
      write(3,*)'ema_mu'
      write(3,*)'dz_over_z0_min'
      write(3,*)'coastal_viscosity'
      write(3,*)'quick_coef'
      write(3,*)'i'
      write(3,*)'j'
      write(3,*)'k'
      write(3,*)'flag3d'
      write(3,*)'lrec'
      write(3,*)'iteration3d'
      write(3,*)'compt1'
      write(3,*)'compt2'
      write(3,*)'compt3'
      write(3,*)'compt4'
      write(3,*)'kount0'
      write(3,*)'kount1'
      write(3,*)'kount2'
      write(3,*)'kount3'
      write(3,*)'kount4'
      write(3,*)'kount5'
      write(3,*)'kount6'
      write(3,*)'kount7'
      write(3,*)'kount8'
      write(3,*)'kount9'
      write(3,*)'kountmod'
      write(3,*)'kountrdv1'
      write(3,*)'kountrdv2'
      write(3,*)'kountrdv3'
      write(3,*)'kountrdv4'
      write(3,*)'substep_advbio'
      write(3,*)'subcycle_exchange'
      write(3,*)'subcycle_onoff'
      write(3,*)'subcycle_synchro'
      write(3,*)'subcycle_modulo'
      write(3,*)'dt_drf_over_dti_fw'
      write(3,*)'iteration3d_restart'
      write(3,*)'filvalshort'
      write(3,*)'quick_filter_points'
      write(3,*)'status'
      write(3,*)'forcedstatus'
      write(3,*)'decision'
      write(3,*)'ncid1'
      write(3,*)'ncid2'
      write(3,*)'dim_x_id'
      write(3,*)'dim_y_id'
      write(3,*)'dim_z_id'
      write(3,*)'dim_t_id'
      write(3,*)'dim_b_id'
      write(3,*)'max_x'
      write(3,*)'max_y'
      write(3,*)'max_z'
      write(3,*)'max_time_counter'
      write(3,*)'max_meteo_time_counter'
      write(3,*)'var_id'
      write(3,*)'var_nftype'
      write(3,*)'meteo_imax'
      write(3,*)'meteo_jmax'
      write(3,*)'meteo_kmax'
      write(3,*)'meteozoom_istr'
      write(3,*)'meteozoom_iend'
      write(3,*)'meteozoom_jstr'
      write(3,*)'meteozoom_jend'
      write(3,*)'meteofull_imax'
      write(3,*)'meteofull_jmax'
      write(3,*)'tide_imax'
      write(3,*)'tide_jmax'
      write(3,*)'tide_kmax'
      write(3,*)'tidezoom_istr'
      write(3,*)'tidezoom_iend'
      write(3,*)'tidezoom_jstr'
      write(3,*)'tidezoom_jend'
      write(3,*)'tidezoom_istr_t'
      write(3,*)'tidezoom_iend_t'
      write(3,*)'tidezoom_jstr_t'
      write(3,*)'tidezoom_jend_t'
      write(3,*)'tidezoom_istr_u'
      write(3,*)'tidezoom_iend_u'
      write(3,*)'tidezoom_jstr_u'
      write(3,*)'tidezoom_jend_u'
      write(3,*)'tidezoom_istr_v'
      write(3,*)'tidezoom_iend_v'
      write(3,*)'tidezoom_jstr_v'
      write(3,*)'tidezoom_jend_v'
      write(3,*)'tidefull_imax'
      write(3,*)'tidefull_jmax'
      write(3,*)'ww3_imax'
      write(3,*)'ww3_jmax'
      write(3,*)'ww3_kmax'
      write(3,*)'ww3_fmax'
      write(3,*)'ww3zoom_istr'
      write(3,*)'ww3zoom_iend'
      write(3,*)'ww3zoom_jstr'
      write(3,*)'ww3zoom_jend'
      write(3,*)'ww3full_imax'
      write(3,*)'ww3full_jmax'
      write(3,*)'jour'
      write(3,*)'kstop'
      write(3,*)'istr'
      write(3,*)'jstr'
      write(3,*)'kstr'
      write(3,*)'tstr'
      write(3,*)'bstr'
      write(3,*)'iend'
      write(3,*)'jend'
      write(3,*)'kend'
      write(3,*)'tend'
      write(3,*)'bend'
      write(3,*)'dimend'
      write(3,*)'kpvwave'
      write(3,*)'give_chanel9'
      write(3,*)'i0'
      write(3,*)'j0'
      write(3,*)'iteration2d'
      write(3,*)'iteration2d_begin'
      write(3,*)'iteration2d_max_now'
      write(3,*)'iteration2d_max_bef'
      write(3,*)'i2dh'
      write(3,*)'l1'
      write(3,*)'l1_sca'
      write(3,*)'l1_vec'
      write(3,*)'l2'
      write(3,*)'l2_sca'
      write(3,*)'l2_vec'
      write(3,*)'l3'
      write(3,*)'l3_sca'
      write(3,*)'l3_vec'
      write(3,*)'len1'
      write(3,*)'nc'
      write(3,*)'nc1'
      write(3,*)'itimets'
      write(3,*)'itimebio'
      write(3,*)'iadvec_ts_hor'
      write(3,*)'iadvec_ts_hor_upwind'
      write(3,*)'iadvec_ts_hor_quickest'
      write(3,*)'iadvec_ts_hor_quickest2'
      write(3,*)'iadvec_ts_ver'
      write(3,*)'iadvec_ts_ver_quickest'
      write(3,*)'iadvec_ts_ver_c2'
      write(3,*)'iadvec_ts_ver_quickest2'
      write(3,*)'iturbulence'
      write(3,*)'istreamf'
      write(3,*)'itime'
      write(3,*)'kts'
      write(3,*)'kuv'
      write(3,*)'ko'
      write(3,*)'jm1'
      write(3,*)'jm2'
      write(3,*)'jp1'
      write(3,*)'jp2'
      write(3,*)'im1'
      write(3,*)'im2'
      write(3,*)'ip1'
      write(3,*)'ip2'
      write(3,*)'iwind'
      write(3,*)'ip'
      write(3,*)'im'
      write(3,*)'jp'
      write(3,*)'jm'
      write(3,*)'kp'
      write(3,*)'km'
      write(3,*)'nbcanal'
      write(3,*)'gridtype1'
      write(3,*)'gridtype2'
      write(3,*)'point1'
      write(3,*)'point2'
      write(3,*)'sender'
      write(3,*)'receiver'
      write(3,*)'ipnoc'
      write(3,*)'i1'
      write(3,*)'i2'
      write(3,*)'i3'
      write(3,*)'i4'
      write(3,*)'i5'
      write(3,*)'i6'
      write(3,*)'i7'
      write(3,*)'i8'
      write(3,*)'i9'
      write(3,*)'i10'
      write(3,*)'i11'
      write(3,*)'j1'
      write(3,*)'j2'
      write(3,*)'j3'
      write(3,*)'j4'
      write(3,*)'j5'
      write(3,*)'j6'
      write(3,*)'j7'
      write(3,*)'j8'
      write(3,*)'j9'
      write(3,*)'j10'
      write(3,*)'j11'
      write(3,*)'k0'
      write(3,*)'k1'
      write(3,*)'k2'
      write(3,*)'k3'
      write(3,*)'k4'
      write(3,*)'k5'
      write(3,*)'k6'
      write(3,*)'k7'
      write(3,*)'k8'
      write(3,*)'k9'
      write(3,*)'kr'
      write(3,*)'kp1'
      write(3,*)'kp2'
      write(3,*)'km1'
      write(3,*)'km2'
      write(3,*)'key'
      write(3,*)'flag'
      write(3,*)'nbomax'
      write(3,*)'nbobuffermax'
      write(3,*)'nbobuffermax_c'
      write(3,*)'flag_stop'
      write(3,*)'itest'
      write(3,*)'itest1'
      write(3,*)'itest2'
      write(3,*)'itest3'
      write(3,*)'ioption'
      write(3,*)'nsmooth'
      write(3,*)'nriver'
      write(3,*)'ian0'
      write(3,*)'imois0'
      write(3,*)'ijour0'
      write(3,*)'iheure0'
      write(3,*)'iminute0'
      write(3,*)'iseconde0'
      write(3,*)'iref'
      write(3,*)'jref'
      write(3,*)'lname1'
      write(3,*)'lname2'
      write(3,*)'lname3'
      write(3,*)'lname4'
      write(3,*)'kmode'
      write(3,*)'kmodemax'
      write(3,*)'fgrid_or_wgrid'
      write(3,*)'fgrid_case'
      write(3,*)'wgrid_case'
      write(3,*)'typegrid'
      write(3,*)'typegrid_monopole'
      write(3,*)'typegrid_file'
      write(3,*)'typegrid_bipole'
      write(3,*)'istar'
      write(3,*)'jstar'
      write(3,*)'istop'
      write(3,*)'jstop'
      write(3,*)'isend'
      write(3,*)'jsend'
      write(3,*)'irecv'
      write(3,*)'jrecv'
      write(3,*)'izoomin'
      write(3,*)'izoomax'
      write(3,*)'jzoomin'
      write(3,*)'jzoomax'
      write(3,*)'iairsea'
      write(3,*)'ialbedo'
      write(3,*)'airseaoption'
      write(3,*)'iobc_f'
      write(3,*)'iobc_wv'
      write(3,*)'iobc_ogcm'
      write(3,*)'iobc_lr'
      write(3,*)'obc_option'
      write(3,*)'iobc_demo_wv'
      write(3,*)'iarchive'
      write(3,*)'imodeltrc'
      write(3,*)'imodelbio'
      write(3,*)'multiple'
      write(3,*)'kvarmax'
      write(3,*)'igesig'
      write(3,*)'isigfile'
      write(3,*)'ksecu'
      write(3,*)'kmaxtide'
      write(3,*)'kmaxtidep1'
      write(3,*)'kminserie'
      write(3,*)'nzctdmax'
      write(3,*)'nsctdmax'
      write(3,*)'ktctdmin'
      write(3,*)'ktctdmax'
      write(3,*)'kdtk_out'
      write(3,*)'k10'
      write(3,*)'k11'
      write(3,*)'nbinco'
      write(3,*)'nbequa'
      write(3,*)'nbsparse'
      write(3,*)'kland'
      write(3,*)'tidestep2'
      write(3,*)'tidestep3'
      write(3,*)'ncfilm_max'
      write(3,*)'in_out_tide'
      write(3,*)'kmax_dof'
      write(3,*)'k_in'
      write(3,*)'k_out'
      write(3,*)'used_unused_dom'
      write(3,*)'mergebathy_sponge'
      write(3,*)'rankhmax'
      write(3,*)'rankhmin'
      write(3,*)'id_tem'
      write(3,*)'id_tem2'
      write(3,*)'id_dtem'
      write(3,*)'id_sal'
      write(3,*)'id_sal2'
      write(3,*)'id_rhp'
      write(3,*)'id_rhop'
      write(3,*)'id_rhom'
      write(3,*)'id_rhf'
      write(3,*)'id_rhc'
      write(3,*)'id_rhpa'
      write(3,*)'id_rhpb'
      write(3,*)'id_z'
      write(3,*)'id_prs'
      write(3,*)'id_now'
      write(3,*)'id_aft'
      write(3,*)'id_eost'
      write(3,*)'id_eoss'
      write(3,*)'id_ssh'
      write(3,*)'id_breaker'
      write(3,*)'id_dz'
      write(3,*)'id_zt'
      write(3,*)'id_tdiv'
      write(3,*)'id_sdiv'
      write(3,*)'id_udiv'
      write(3,*)'id_vdiv'
      write(3,*)'id_bdiv'
      write(3,*)'id_bnegdif'
      write(3,*)'id_biobefo'
      write(3,*)'id_bioaftr'
      write(3,*)'id_u_now'
      write(3,*)'id_v_now'
      write(3,*)'id_u_rot'
      write(3,*)'id_v_rot'
      write(3,*)'id_ffreq'
      write(3,*)'id_coriolis1'
      write(3,*)'id_coriolis2'
      write(3,*)'id_flx'
      write(3,*)'id_fly'
      write(3,*)'id_uflx'
      write(3,*)'id_ufly'
      write(3,*)'id_vflx'
      write(3,*)'id_vfly'
      write(3,*)'id_tcn'
      write(3,*)'id_scn'
      write(3,*)'id_gradssh'
      write(3,*)'id_ofactort'
      write(3,*)'id_ofactors'
      write(3,*)'id_hybcoefu'
      write(3,*)'id_hybcoefv'
      write(3,*)'id_dxdydz'
      write(3,*)'id_wb'
      write(3,*)'id_prod'
      write(3,*)'id_buoy'
      write(3,*)'iDtOvRhCp'
      write(3,*)'id_ncu'
      write(3,*)'id_ncv'
      write(3,*)'id_webio'
      write(3,*)'id_varbef2'
      write(3,*)'id_kh_over_dz'
      write(3,*)'id_bihar_lim'
      write(3,*)'id_veltot'
      write(3,*)'id_velexp'
      write(3,*)'id_velimp'
      write(3,*)'id_wdrifter'
      write(3,*)'kreadgroup'
      write(3,*)'kreadgroupmax'
      write(3,*)'kread1'
      write(3,*)'kread2'
      write(3,*)'modulo_biotimestep'
      write(3,*)'obctime_bef2'
      write(3,*)'obctime_bef'
      write(3,*)'obctime_aft'
      write(3,*)'obctime_aft2'
      write(3,*)'obctime_order'
      write(3,*)'sch_imp_ts_u_loc'
      write(3,*)'sch_imp_ts_v_loc'
      write(3,*)'sch_imp_ts_u_glb'
      write(3,*)'sch_imp_ts_v_glb'
      write(3,*)'sch_imp_tke_u_loc'
      write(3,*)'sch_imp_tke_v_loc'
      write(3,*)'sch_imp_tke_u_glb'
      write(3,*)'sch_imp_tke_v_glb'
      write(3,*)'looplimit_hor'
      write(3,*)'ihybsig'
      write(3,*)'nhybsig'
      write(3,*)'tke_surf'
      write(3,*)'grh_out_mi'
      write(3,*)'nest_onoff_in'
      write(3,*)'nest_onoff_demo'
      write(3,*)'nest_full_in'
      write(3,*)'nest_full_out'
      write(3,*)'nest_onoff_out'
      write(3,*)'ncmin_airsea'
      write(3,*)'ncmin_river'
      write(3,*)'eos_author'
      write(3,*)'eos_comprs'
      write(3,*)'eos_linear'
      write(3,*)'obcfreeorfix'
      write(3,*)'iwve'
      write(3,*)'wave_obc_type'
      write(3,*)'dataperwavefile'
      write(3,*)'sp_or_db'
      write(3,*)'kbu'
      write(3,*)'kbumax'
      write(3,*)'kbumax_glb'
      write(3,*)'n_element'
      write(3,*)'relaxtype_ts'
      write(3,*)'obctype_ts'
      write(3,*)'obctype_p'
      write(3,*)'restart_file_y_or_n'
      write(3,*)'ncid'
      write(3,*)'x_xhl_dim'
      write(3,*)'y_xhl_dim'
      write(3,*)'z_xhl_dim'
      write(3,*)'x_yhl_dim'
      write(3,*)'y_yhl_dim'
      write(3,*)'z_yhl_dim'
      write(3,*)'x_zhl_dim'
      write(3,*)'y_zhl_dim'
      write(3,*)'z_zhl_dim'
      write(3,*)'x_zl_dim'
      write(3,*)'y_zl_dim'
      write(3,*)'z_zl_dim'
      write(3,*)'time_dim'
      write(3,*)'i_t_dim'
      write(3,*)'j_t_dim'
      write(3,*)'k_t_dim'
      write(3,*)'i_w_dim'
      write(3,*)'j_w_dim'
      write(3,*)'k_w_dim'
      write(3,*)'i_u_dim'
      write(3,*)'j_u_dim'
      write(3,*)'k_u_dim'
      write(3,*)'i_v_dim'
      write(3,*)'j_v_dim'
      write(3,*)'k_v_dim'
      write(3,*)'i_f_dim'
      write(3,*)'j_f_dim'
      write(3,*)'k_f_dim'
      write(3,*)'dayindex_size'
      write(3,*)'tide_year_min'
      write(3,*)'airsea_year_min'
      write(3,*)'tide_year_max'
      write(3,*)'airsea_year_max'
      write(3,*)'wave_year_min'
      write(3,*)'irelaxsst'
      write(3,*)'removetide'
      write(3,*)'ofl_rotation'
      write(3,*)'ofl_rhp'
      write(3,*)'ofl_tke'
      write(3,*)'ofl_surflux'
      write(3,*)'ale_selected'
      write(3,*)'flag_asselin'
      write(3,*)'freq'
      write(3,*)'dirw'
      write(3,*)'dirw_beg'
      write(3,*)'dirw_end'
      write(3,*)'freq_beg'
      write(3,*)'freq_end'
      write(3,*)'year_now'
      write(3,*)'month_now'
      write(3,*)'day_now'
      write(3,*)'hour_now'
      write(3,*)'minute_now'
      write(3,*)'second_now'
      write(3,*)'nd_send_est'
      write(3,*)'nd_send_ouest'
      write(3,*)'nd_send_nord'
      write(3,*)'nd_send_sud'
      write(3,*)'nd_send_out'
      write(3,*)'nd_recv_est'
      write(3,*)'nd_recv_ouest'
      write(3,*)'nd_recv_nord'
      write(3,*)'nd_recv_sud'
      write(3,*)'nd_send_sudouest'
      write(3,*)'nd_send_sudest'
      write(3,*)'nd_send_nordouest'
      write(3,*)'nd_send_nordest'
      write(3,*)'nd_recv_sudouest'
      write(3,*)'nd_recv_sudest'
      write(3,*)'nd_recv_nordouest'
      write(3,*)'nd_recv_nordest'
      write(3,*)'initial_main_status'
      write(3,*)'offline_init_status'
      write(3,*)'meteo_sealand_mask'
      write(3,*)'grid_i0'
      write(3,*)'grid_j0'
      write(3,*)'ifb'
      write(3,*)'dbefore'
      write(3,*)'dnow'
      write(3,*)'dafter'
      write(3,*)'dim_airsea'
      write(3,*)'rhp_zavr_xy'
      write(3,*)'ssr_id'
      write(3,*)'ir_id'
      write(3,*)'rain_id'
      write(3,*)'t2m_id'
      write(3,*)'t0m_id'
      write(3,*)'abl_id'
      write(3,*)'dp2m_id'
      write(3,*)'u10m_id'
      write(3,*)'v10m_id'
      write(3,*)'u100m_id'
      write(3,*)'v100m_id'
      write(3,*)'p0m_id'
      write(3,*)'ustrs_id'
      write(3,*)'vstrs_id'
      write(3,*)'slhf_id'
      write(3,*)'netir_id'
      write(3,*)'sshf_id'
      write(3,*)'t_flux_cumul'
      write(3,*)'s_flux_cumul'
      write(3,*)'cumuldeltaflux'
      write(3,*)'som0'
      write(3,*)'som2'
      write(3,*)'ssh_reservoir'
      write(3,*)'idealflux'
      write(3,*)'sum0'
      write(3,*)'sum1'
      write(3,*)'sum2'
      write(3,*)'sum3'
      write(3,*)'sum4'
      write(3,*)'sum5'
      write(3,*)'sum6'
      write(3,*)'sum7'
      write(3,*)'sum8'
      write(3,*)'sum9'
      write(3,*)'sum10'
      write(3,*)'sum11'
      write(3,*)'sum12'
      write(3,*)'sum13'
      write(3,*)'time0'
      write(3,*)'time1'
      write(3,*)'time2'
      write(3,*)'small2'
      write(3,*)'x1_r8'
      write(3,*)'x2_r8'
      write(3,*)'x3_r8'
      write(3,*)'x4_r8'
      write(3,*)'sum0glb'
      write(3,*)'sum1glb'
      write(3,*)'sum2glb'
      write(3,*)'sum3glb'
      write(3,*)'sum4glb'
      write(3,*)'sum5glb'
      write(3,*)'sum6glb'
      write(3,*)'emin'
      write(3,*)'epsmin'
      write(3,*)'elapsedtime_out'
      write(3,*)'elapsedtime_rst'
      write(3,*)'elapsedtime_now'
      write(3,*)'elapsedtime_last_writing'
      write(3,*)'elapsedtime_bef'
      write(3,*)'elapsedtime_aft'
      write(3,*)'cpu_seconds'
      write(3,*)'alpha'
      write(3,*)'eos_pgfzref'
      write(3,*)'eos_tkezref'
      write(3,*)'filvalr8'
      write(3,*)'dti_fw_i4'
      write(3,*)'elapsedtime_now_i4'
      write(3,*)'elapsedtime_now_i8'
      write(3,*)'elapsedtime_now_r16'
! jalon2 dyn_restart_check ne pas effacer
      close(3)
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !
      end subroutine dyn_restart_bounds_w
!.......................................................................
      subroutine dyn_restart_bounds_r
      use module_principal
      use module_parallele !#MPI
      implicit none
      write(texte60,'(a,i0)')trim(restartdir_in)//'dynrestartbounds',par%rank !18-07-20
      open(unit=3,file=trim(texte60))
! jalon3 dyn_restart_check ne pas effacer
      write(3,*)'ioffline_prv'
      write(3,*)'zone1bioflux_w'
      if(allocated(zone1bioflux_glb))write(3,*)'zone1bioflux_glb',1+ubound(zone1bioflux_glb)-lbound(zone1bioflux_glb)
      if(allocated(zone1bioflux_u))write(3,*)'zone1bioflux_u',1+ubound(zone1bioflux_u)-lbound(zone1bioflux_u)
      if(allocated(zone1bioflux_v))write(3,*)'zone1bioflux_v',1+ubound(zone1bioflux_v)-lbound(zone1bioflux_v)
      if(allocated(zone1tendancebio_glb))write(3,*)'zone1tendancebio_glb',1+ubound(zone1tendancebio_glb)-lbound(zone1tendancebio_glb)
      if(allocated(zone1botsurfbio_glb))write(3,*)'zone1botsurfbio_glb',1+ubound(zone1botsurfbio_glb)-lbound(zone1botsurfbio_glb)
      if(allocated(zone1bioflux_glb_in))write(3,*)'zone1bioflux_glb_in',1+ubound(zone1bioflux_glb_in)-lbound(zone1bioflux_glb_in)
      if(allocated(zone1bioflux_u_in))write(3,*)'zone1bioflux_u_in',1+ubound(zone1bioflux_u_in)-lbound(zone1bioflux_u_in)
      if(allocated(zone1bioflux_v_in))write(3,*)'zone1bioflux_v_in',1+ubound(zone1bioflux_v_in)-lbound(zone1bioflux_v_in)
      if(allocated(zone1bioflux_glb_out))write(3,*)'zone1bioflux_glb_out',1+ubound(zone1bioflux_glb_out)-lbound(zone1bioflux_glb_out)
      if(allocated(zone1bioflux_u_out))write(3,*)'zone1bioflux_u_out',1+ubound(zone1bioflux_u_out)-lbound(zone1bioflux_u_out)
      if(allocated(zone1bioflux_v_out))write(3,*)'zone1bioflux_v_out',1+ubound(zone1bioflux_v_out)-lbound(zone1bioflux_v_out)
      if(allocated(zone1biomasst0))write(3,*)'zone1biomasst0',1+ubound(zone1biomasst0)-lbound(zone1biomasst0)
      if(allocated(analysis3dmatrix))write(3,*)'analysis3dmatrix',1+ubound(analysis3dmatrix)-lbound(analysis3dmatrix)
      if(allocated(ema1_s_i))write(3,*)'ema1_s_i',1+ubound(ema1_s_i)-lbound(ema1_s_i)
      if(allocated(ema2_s_i))write(3,*)'ema2_s_i',1+ubound(ema2_s_i)-lbound(ema2_s_i)
      if(allocated(ema1_s_j))write(3,*)'ema1_s_j',1+ubound(ema1_s_j)-lbound(ema1_s_j)
      if(allocated(ema2_s_j))write(3,*)'ema2_s_j',1+ubound(ema2_s_j)-lbound(ema2_s_j)
      if(allocated(c_wave_mode))write(3,*)'c_wave_mode',1+ubound(c_wave_mode)-lbound(c_wave_mode)
      if(allocated(pcoefmode_t))write(3,*)'pcoefmode_t',1+ubound(pcoefmode_t)-lbound(pcoefmode_t)
      if(allocated(ucoefmode_t))write(3,*)'ucoefmode_t',1+ubound(ucoefmode_t)-lbound(ucoefmode_t)
      if(allocated(vcoefmode_t))write(3,*)'vcoefmode_t',1+ubound(vcoefmode_t)-lbound(vcoefmode_t)
      if(allocated(rhrefmode_t))write(3,*)'rhrefmode_t',1+ubound(rhrefmode_t)-lbound(rhrefmode_t)
      if(allocated(ema1_q_i))write(3,*)'ema1_q_i',1+ubound(ema1_q_i)-lbound(ema1_q_i)
      if(allocated(ema2_q_i))write(3,*)'ema2_q_i',1+ubound(ema2_q_i)-lbound(ema2_q_i)
      if(allocated(ema1_q_j))write(3,*)'ema1_q_j',1+ubound(ema1_q_j)-lbound(ema1_q_j)
      if(allocated(ema2_q_j))write(3,*)'ema2_q_j',1+ubound(ema2_q_j)-lbound(ema2_q_j)
      if(allocated(uv_wmode_t))write(3,*)'uv_wmode_t',1+ubound(uv_wmode_t)-lbound(uv_wmode_t)
      if(allocated(analysis_p3d_t))write(3,*)'analysis_p3d_t',1+ubound(analysis_p3d_t)-lbound(analysis_p3d_t)
      if(allocated(analysis_u3d_t))write(3,*)'analysis_u3d_t',1+ubound(analysis_u3d_t)-lbound(analysis_u3d_t)
      if(allocated(analysis_v3d_t))write(3,*)'analysis_v3d_t',1+ubound(analysis_v3d_t)-lbound(analysis_v3d_t)
      if(allocated(p3dmode_cos_t))write(3,*)'p3dmode_cos_t',1+ubound(p3dmode_cos_t)-lbound(p3dmode_cos_t)
      if(allocated(p3dmode_sin_t))write(3,*)'p3dmode_sin_t',1+ubound(p3dmode_sin_t)-lbound(p3dmode_sin_t)
      if(allocated(u3dmode_cos_t))write(3,*)'u3dmode_cos_t',1+ubound(u3dmode_cos_t)-lbound(u3dmode_cos_t)
      if(allocated(u3dmode_sin_t))write(3,*)'u3dmode_sin_t',1+ubound(u3dmode_sin_t)-lbound(u3dmode_sin_t)
      if(allocated(v3dmode_cos_t))write(3,*)'v3dmode_cos_t',1+ubound(v3dmode_cos_t)-lbound(v3dmode_cos_t)
      if(allocated(v3dmode_sin_t))write(3,*)'v3dmode_sin_t',1+ubound(v3dmode_sin_t)-lbound(v3dmode_sin_t)
      if(allocated(analysetide3d_u))write(3,*)'analysetide3d_u',1+ubound(analysetide3d_u)-lbound(analysetide3d_u)
      if(allocated(analysetide3d_v))write(3,*)'analysetide3d_v',1+ubound(analysetide3d_v)-lbound(analysetide3d_v)
      if(allocated(analysetide3d_t))write(3,*)'analysetide3d_t',1+ubound(analysetide3d_t)-lbound(analysetide3d_t)
      if(allocated(vel3dtidecosout_u))write(3,*)'vel3dtidecosout_u',1+ubound(vel3dtidecosout_u)-lbound(vel3dtidecosout_u)
      if(allocated(vel3dtidesinout_u))write(3,*)'vel3dtidesinout_u',1+ubound(vel3dtidesinout_u)-lbound(vel3dtidesinout_u)
      if(allocated(vel3dtidecosout_v))write(3,*)'vel3dtidecosout_v',1+ubound(vel3dtidecosout_v)-lbound(vel3dtidecosout_v)
      if(allocated(vel3dtidesinout_v))write(3,*)'vel3dtidesinout_v',1+ubound(vel3dtidesinout_v)-lbound(vel3dtidesinout_v)
      if(allocated(rhptidecosout_t))write(3,*)'rhptidecosout_t',1+ubound(rhptidecosout_t)-lbound(rhptidecosout_t)
      if(allocated(rhptidesinout_t))write(3,*)'rhptidesinout_t',1+ubound(rhptidesinout_t)-lbound(rhptidesinout_t)
      if(allocated(bio_relax_north))write(3,*)'bio_relax_north',1+ubound(bio_relax_north)-lbound(bio_relax_north)
      if(allocated(bio_relax_south))write(3,*)'bio_relax_south',1+ubound(bio_relax_south)-lbound(bio_relax_south)
      if(allocated(bio_relax_east))write(3,*)'bio_relax_east',1+ubound(bio_relax_east)-lbound(bio_relax_east)
      if(allocated(bio_relax_west))write(3,*)'bio_relax_west',1+ubound(bio_relax_west)-lbound(bio_relax_west)
      if(allocated(bio_relax_full))write(3,*)'bio_relax_full',1+ubound(bio_relax_full)-lbound(bio_relax_full)
      if(allocated(analysetide_u))write(3,*)'analysetide_u',1+ubound(analysetide_u)-lbound(analysetide_u)
      if(allocated(analysetide_v))write(3,*)'analysetide_v',1+ubound(analysetide_v)-lbound(analysetide_v)
      if(allocated(analysetide_w))write(3,*)'analysetide_w',1+ubound(analysetide_w)-lbound(analysetide_w)
      if(allocated(analysetide2_w))write(3,*)'analysetide2_w',1+ubound(analysetide2_w)-lbound(analysetide2_w)
      if(allocated(cocosisi_f))write(3,*)'cocosisi_f',1+ubound(cocosisi_f)-lbound(cocosisi_f)
      if(allocated(gridrotcos_f))write(3,*)'gridrotcos_f',1+ubound(gridrotcos_f)-lbound(gridrotcos_f)
      if(allocated(gridrotsin_f))write(3,*)'gridrotsin_f',1+ubound(gridrotsin_f)-lbound(gridrotsin_f)
      if(allocated(tideanalysismatrix))write(3,*)'tideanalysismatrix',1+ubound(tideanalysismatrix)-lbound(tideanalysismatrix)
      if(allocated(iaveraged_in))write(3,*)'iaveraged_in',1+ubound(iaveraged_in)-lbound(iaveraged_in)
      if(allocated(iaveraged_out))write(3,*)'iaveraged_out',1+ubound(iaveraged_out)-lbound(iaveraged_out)
      if(allocated(iaverag1d_in))write(3,*)'iaverag1d_in',1+ubound(iaverag1d_in)-lbound(iaverag1d_in)
      if(allocated(iaverag1d_out))write(3,*)'iaverag1d_out',1+ubound(iaverag1d_out)-lbound(iaverag1d_out)
      if(allocated(light_kpar2_w))write(3,*)'light_kpar2_w',1+ubound(light_kpar2_w)-lbound(light_kpar2_w)
      if(allocated(tidepotential_w))write(3,*)'tidepotential_w',1+ubound(tidepotential_w)-lbound(tidepotential_w)
      if(allocated(sshtidecos_w))write(3,*)'sshtidecos_w',1+ubound(sshtidecos_w)-lbound(sshtidecos_w)
      if(allocated(sshtidesin_w))write(3,*)'sshtidesin_w',1+ubound(sshtidesin_w)-lbound(sshtidesin_w)
      if(allocated(veltidecos_u))write(3,*)'veltidecos_u',1+ubound(veltidecos_u)-lbound(veltidecos_u)
      if(allocated(veltidesin_u))write(3,*)'veltidesin_u',1+ubound(veltidesin_u)-lbound(veltidesin_u)
      if(allocated(veltidecos_v))write(3,*)'veltidecos_v',1+ubound(veltidecos_v)-lbound(veltidecos_v)
      if(allocated(veltidesin_v))write(3,*)'veltidesin_v',1+ubound(veltidesin_v)-lbound(veltidesin_v)
      if(allocated(potidecos_w))write(3,*)'potidecos_w',1+ubound(potidecos_w)-lbound(potidecos_w)
      if(allocated(potidesin_w))write(3,*)'potidesin_w',1+ubound(potidesin_w)-lbound(potidesin_w)
      if(allocated(veltidecosout_u))write(3,*)'veltidecosout_u',1+ubound(veltidecosout_u)-lbound(veltidecosout_u)
      if(allocated(veltidesinout_u))write(3,*)'veltidesinout_u',1+ubound(veltidesinout_u)-lbound(veltidesinout_u)
      if(allocated(veltidecosout_v))write(3,*)'veltidecosout_v',1+ubound(veltidecosout_v)-lbound(veltidecosout_v)
      if(allocated(veltidesinout_v))write(3,*)'veltidesinout_v',1+ubound(veltidesinout_v)-lbound(veltidesinout_v)
      if(allocated(sshtidecosout_w))write(3,*)'sshtidecosout_w',1+ubound(sshtidecosout_w)-lbound(sshtidecosout_w)
      if(allocated(sshtidesinout_w))write(3,*)'sshtidesinout_w',1+ubound(sshtidesinout_w)-lbound(sshtidesinout_w)
      if(allocated(ssh3dtidecosout_w))write(3,*)'ssh3dtidecosout_w',1+ubound(ssh3dtidecosout_w)-lbound(ssh3dtidecosout_w)
      if(allocated(ssh3dtidesinout_w))write(3,*)'ssh3dtidesinout_w',1+ubound(ssh3dtidesinout_w)-lbound(ssh3dtidesinout_w)
      if(allocated(passetide))write(3,*)'passetide',1+ubound(passetide)-lbound(passetide)
      if(allocated(ftide))write(3,*)'ftide',1+ubound(ftide)-lbound(ftide)
      if(allocated(utide))write(3,*)'utide',1+ubound(utide)-lbound(utide)
      if(allocated(northflux_sumt_v))write(3,*)'northflux_sumt_v',1+ubound(northflux_sumt_v)-lbound(northflux_sumt_v)
      if(allocated(southflux_sumt_v))write(3,*)'southflux_sumt_v',1+ubound(southflux_sumt_v)-lbound(southflux_sumt_v)
      if(allocated(drifter_l))write(3,*)'drifter_l',1+ubound(drifter_l)-lbound(drifter_l)
      if(allocated(velbot3d2d_u))write(3,*)'velbot3d2d_u',1+ubound(velbot3d2d_u)-lbound(velbot3d2d_u)
      if(allocated(velbot3d2d_v))write(3,*)'velbot3d2d_v',1+ubound(velbot3d2d_v)-lbound(velbot3d2d_v)
      if(allocated(velbot_u))write(3,*)'velbot_u',1+ubound(velbot_u)-lbound(velbot_u)
      if(allocated(velbot_v))write(3,*)'velbot_v',1+ubound(velbot_v)-lbound(velbot_v)
      if(allocated(kvectorpeak_i))write(3,*)'kvectorpeak_i',1+ubound(kvectorpeak_i)-lbound(kvectorpeak_i)
      if(allocated(kvectorpeak_j))write(3,*)'kvectorpeak_j',1+ubound(kvectorpeak_j)-lbound(kvectorpeak_j)
      if(allocated(drifter_send_canal))write(3,*)'drifter_send_canal',1+ubound(drifter_send_canal)-lbound(drifter_send_canal)
      if(allocated(drifter_recv_canal))write(3,*)'drifter_recv_canal',1+ubound(drifter_recv_canal)-lbound(drifter_recv_canal)
      if(allocated(qwave_j_w))write(3,*)'qwave_j_w',1+ubound(qwave_j_w)-lbound(qwave_j_w)
      if(allocated(qwave_i_w))write(3,*)'qwave_i_w',1+ubound(qwave_i_w)-lbound(qwave_i_w)
      if(allocated(velwave_i_u))write(3,*)'velwave_i_u',1+ubound(velwave_i_u)-lbound(velwave_i_u)
      if(allocated(velwave_j_u))write(3,*)'velwave_j_u',1+ubound(velwave_j_u)-lbound(velwave_j_u)
      if(allocated(velwave_i_v))write(3,*)'velwave_i_v',1+ubound(velwave_i_v)-lbound(velwave_i_v)
      if(allocated(velwave_j_v))write(3,*)'velwave_j_v',1+ubound(velwave_j_v)-lbound(velwave_j_v)
      if(allocated(q_t))write(3,*)'q_t',1+ubound(q_t)-lbound(q_t)
      if(allocated(sshwave_j_w))write(3,*)'sshwave_j_w',1+ubound(sshwave_j_w)-lbound(sshwave_j_w)
      if(allocated(sshwave_i_w))write(3,*)'sshwave_i_w',1+ubound(sshwave_i_w)-lbound(sshwave_i_w)
      if(allocated(mc))write(3,*)'mc',1+ubound(mc)-lbound(mc)
      if(allocated(anyv3dr4))write(3,*)'anyv3dr4',1+ubound(anyv3dr4)-lbound(anyv3dr4)
      if(allocated(mc_out))write(3,*)'mc_out',1+ubound(mc_out)-lbound(mc_out)
      if(allocated(velbarwave_i_v))write(3,*)'velbarwave_i_v',1+ubound(velbarwave_i_v)-lbound(velbarwave_i_v)
      if(allocated(qavr_t))write(3,*)'qavr_t',1+ubound(qavr_t)-lbound(qavr_t)
      if(allocated(nhpgf_u))write(3,*)'nhpgf_u',1+ubound(nhpgf_u)-lbound(nhpgf_u)
      if(allocated(nhpgf_v))write(3,*)'nhpgf_v',1+ubound(nhpgf_v)-lbound(nhpgf_v)
      if(allocated(dsigr4_t))write(3,*)'dsigr4_t',1+ubound(dsigr4_t)-lbound(dsigr4_t)
      if(allocated(sig1dpom_t))write(3,*)'sig1dpom_t',1+ubound(sig1dpom_t)-lbound(sig1dpom_t)
      if(allocated(obc_q_i))write(3,*)'obc_q_i',1+ubound(obc_q_i)-lbound(obc_q_i)
      if(allocated(obc_ub_i))write(3,*)'obc_ub_i',1+ubound(obc_ub_i)-lbound(obc_ub_i)
      if(allocated(dsig_w))write(3,*)'dsig_w',1+ubound(dsig_w)-lbound(dsig_w)
      if(allocated(retarded_pss_w))write(3,*)'retarded_pss_w',1+ubound(retarded_pss_w)-lbound(retarded_pss_w)
      if(allocated(velwf_u))write(3,*)'velwf_u',1+ubound(velwf_u)-lbound(velwf_u)
      if(allocated(velwf_v))write(3,*)'velwf_v',1+ubound(velwf_v)-lbound(velwf_v)
      if(allocated(kw_w))write(3,*)'kw_w',1+ubound(kw_w)-lbound(kw_w)
      if(allocated(q2davr_w))write(3,*)'q2davr_w',1+ubound(q2davr_w)-lbound(q2davr_w)
      if(allocated(breaker2d_t))write(3,*)'breaker2d_t',1+ubound(breaker2d_t)-lbound(breaker2d_t)
      if(allocated(invdx_t))write(3,*)'invdx_t',1+ubound(invdx_t)-lbound(invdx_t)
      if(allocated(invdx_f))write(3,*)'invdx_f',1+ubound(invdx_f)-lbound(invdx_f)
      if(allocated(invdx_u))write(3,*)'invdx_u',1+ubound(invdx_u)-lbound(invdx_u)
      if(allocated(invdy_u))write(3,*)'invdy_u',1+ubound(invdy_u)-lbound(invdy_u)
      if(allocated(invdy_t))write(3,*)'invdy_t',1+ubound(invdy_t)-lbound(invdy_t)
      if(allocated(invdy_f))write(3,*)'invdy_f',1+ubound(invdy_f)-lbound(invdy_f)
      if(allocated(invdy_v))write(3,*)'invdy_v',1+ubound(invdy_v)-lbound(invdy_v)
      if(allocated(invdx_v))write(3,*)'invdx_v',1+ubound(invdx_v)-lbound(invdx_v)
      if(allocated(invdxdy_t))write(3,*)'invdxdy_t',1+ubound(invdxdy_t)-lbound(invdxdy_t)
      if(allocated(invdxdy_u))write(3,*)'invdxdy_u',1+ubound(invdxdy_u)-lbound(invdxdy_u)
      if(allocated(invdxdy_v))write(3,*)'invdxdy_v',1+ubound(invdxdy_v)-lbound(invdxdy_v)
      if(allocated(nhpgf2d_u))write(3,*)'nhpgf2d_u',1+ubound(nhpgf2d_u)-lbound(nhpgf2d_u)
      if(allocated(nhpgf2d_v))write(3,*)'nhpgf2d_v',1+ubound(nhpgf2d_v)-lbound(nhpgf2d_v)
      if(allocated(sshmax_w))write(3,*)'sshmax_w',1+ubound(sshmax_w)-lbound(sshmax_w)
      if(allocated(sshmin_w))write(3,*)'sshmin_w',1+ubound(sshmin_w)-lbound(sshmin_w)
      if(allocated(sshmin_tmp_w))write(3,*)'sshmin_tmp_w',1+ubound(sshmin_tmp_w)-lbound(sshmin_tmp_w)
      if(allocated(rwavebreak_t))write(3,*)'rwavebreak_t',1+ubound(rwavebreak_t)-lbound(rwavebreak_t)
      if(allocated(breaker_t))write(3,*)'breaker_t',1+ubound(breaker_t)-lbound(breaker_t)
      if(allocated(hs_w))write(3,*)'hs_w',1+ubound(hs_w)-lbound(hs_w)
      if(allocated(cx_upper_t))write(3,*)'cx_upper_t',1+ubound(cx_upper_t)-lbound(cx_upper_t)
      if(allocated(cy_upper_t))write(3,*)'cy_upper_t',1+ubound(cy_upper_t)-lbound(cy_upper_t)
      if(allocated(c_lower_t))write(3,*)'c_lower_t',1+ubound(c_lower_t)-lbound(c_lower_t)
      if(allocated(slopemax_u))write(3,*)'slopemax_u',1+ubound(slopemax_u)-lbound(slopemax_u)
      if(allocated(pgfwave_u))write(3,*)'pgfwave_u',1+ubound(pgfwave_u)-lbound(pgfwave_u)
      if(allocated(sshrelax_j_w))write(3,*)'sshrelax_j_w',1+ubound(sshrelax_j_w)-lbound(sshrelax_j_w)
      if(allocated(sshrelax_i_w))write(3,*)'sshrelax_i_w',1+ubound(sshrelax_i_w)-lbound(sshrelax_i_w)
      write(3,*)'asselin_nh'
      write(3,*)'wetdry_cstnh'
      write(3,*)'cfl_nh'
      write(3,*)'cwavepeak'
      write(3,*)'periodpeak'
      write(3,*)'freq2pipeak'
      write(3,*)'kvectorpeak'
      write(3,*)'kvector'
      write(3,*)'mukvector'
      write(3,*)'nhpgf_reduce'
      write(3,*)'acm_speed'
      write(3,*)'inv_brkh'
      write(3,*)'inv_brkslope2'
      write(3,*)'inv_brkvar'
      write(3,*)'brk_crit_slope'
      write(3,*)'brk_crit_var'
      write(3,*)'brk_crit_h'
      write(3,*)'brk_crit_r'
      write(3,*)'nh2d_graph_period'
      write(3,*)'nh2d_graph_spinup'
      write(3,*)'sum_avr_hs'
      write(3,*)'nh_frozensigma'
      write(3,*)'nh_wavebreakscheme'
      write(3,*)'flag_adve2d'
      write(3,*)'flag_adve3d'
      write(3,*)'flag_nh3d'
      write(3,*)'flag_nh3d_none'
      write(3,*)'flag_nh3d_nosplit_uv'
      write(3,*)'flag_nh3d_nosplit_tsuv'
      write(3,*)'flag_nh3d_timesplit_tsuv'
      write(3,*)'flag_timesplitting'
      write(3,*)'flag_timesplitting_adve_uv'
      write(3,*)'flagspo_i1'
      write(3,*)'flagspo_i2'
      write(3,*)'flagspo_j1'
      write(3,*)'flagspo_j2'
      write(3,*)'flag_ksloffline'
      write(3,*)'flag_groundwater'
      write(3,*)'flag_surfriver'
      write(3,*)'ofl_reversedtime'
      write(3,*)'obcnh_scheme'
      write(3,*)'obc_scheme'
      write(3,*)'drifter_random_w'
      write(3,*)'flag_buo_w'
      write(3,*)'flag_ogcmtidemixing'
      write(3,*)'flag_ogcm_instab'
      write(3,*)'flag_ofactor'
      write(3,*)'flag_omega_cumul'
      write(3,*)'flag_negdif_ver'
      write(3,*)'flag_negdif_hor'
      write(3,*)'flag_z0_macro'
      write(3,*)'oasis_symsym_onoff'
      write(3,*)'oasis_symsym_retrots'
      if(allocated(variance_timeaveraged))write(3,*)'variance_timeaveraged',1+ubound(variance_timeaveraged)-lbound(variance_timeaveraged)
      if(allocated(ssh_timeaveraged))write(3,*)'ssh_timeaveraged',1+ubound(ssh_timeaveraged)-lbound(ssh_timeaveraged)
      if(allocated(flx2d_timeaveraged))write(3,*)'flx2d_timeaveraged',1+ubound(flx2d_timeaveraged)-lbound(flx2d_timeaveraged)
      if(allocated(fly2d_timeaveraged))write(3,*)'fly2d_timeaveraged',1+ubound(fly2d_timeaveraged)-lbound(fly2d_timeaveraged)
      if(allocated(flx3d_timeaveraged))write(3,*)'flx3d_timeaveraged',1+ubound(flx3d_timeaveraged)-lbound(flx3d_timeaveraged)
      if(allocated(fly3d_timeaveraged))write(3,*)'fly3d_timeaveraged',1+ubound(fly3d_timeaveraged)-lbound(fly3d_timeaveraged)
      if(allocated(u_euler_timeaveraged))write(3,*)'u_euler_timeaveraged',1+ubound(u_euler_timeaveraged)-lbound(u_euler_timeaveraged)
      if(allocated(v_euler_timeaveraged))write(3,*)'v_euler_timeaveraged',1+ubound(v_euler_timeaveraged)-lbound(v_euler_timeaveraged)
      if(allocated(tkeb_w))write(3,*)'tkeb_w',1+ubound(tkeb_w)-lbound(tkeb_w)
      if(allocated(tken_w))write(3,*)'tken_w',1+ubound(tken_w)-lbound(tken_w)
      if(allocated(tkea_w))write(3,*)'tkea_w',1+ubound(tkea_w)-lbound(tkea_w)
      if(allocated(tkle_w))write(3,*)'tkle_w',1+ubound(tkle_w)-lbound(tkle_w)
      if(allocated(tkll_w))write(3,*)'tkll_w',1+ubound(tkll_w)-lbound(tkll_w)
      if(allocated(epsb_w))write(3,*)'epsb_w',1+ubound(epsb_w)-lbound(epsb_w)
      if(allocated(epsn_w))write(3,*)'epsn_w',1+ubound(epsn_w)-lbound(epsn_w)
      if(allocated(epsa_w))write(3,*)'epsa_w',1+ubound(epsa_w)-lbound(epsa_w)
      if(allocated(gradssh_u))write(3,*)'gradssh_u',1+ubound(gradssh_u)-lbound(gradssh_u)
      if(allocated(gradssh_v))write(3,*)'gradssh_v',1+ubound(gradssh_v)-lbound(gradssh_v)
      if(allocated(pgf_u))write(3,*)'pgf_u',1+ubound(pgf_u)-lbound(pgf_u)
      if(allocated(pgf_v))write(3,*)'pgf_v',1+ubound(pgf_v)-lbound(pgf_v)
      if(allocated(dz_u))write(3,*)'dz_u',1+ubound(dz_u)-lbound(dz_u)
      if(allocated(dz_v))write(3,*)'dz_v',1+ubound(dz_v)-lbound(dz_v)
      if(allocated(dz_t))write(3,*)'dz_t',1+ubound(dz_t)-lbound(dz_t)
      if(allocated(depth_t))write(3,*)'depth_t',1+ubound(depth_t)-lbound(depth_t)
      if(allocated(depth_u))write(3,*)'depth_u',1+ubound(depth_u)-lbound(depth_u)
      if(allocated(depth_v))write(3,*)'depth_v',1+ubound(depth_v)-lbound(depth_v)
      if(allocated(depth_w))write(3,*)'depth_w',1+ubound(depth_w)-lbound(depth_w)
      if(allocated(depth_f))write(3,*)'depth_f',1+ubound(depth_f)-lbound(depth_f)
      if(allocated(stokesforces_u))write(3,*)'stokesforces_u',1+ubound(stokesforces_u)-lbound(stokesforces_u)
      if(allocated(stokesforces_v))write(3,*)'stokesforces_v',1+ubound(stokesforces_v)-lbound(stokesforces_v)
      if(allocated(r_mangrove))write(3,*)'r_mangrove',1+ubound(r_mangrove)-lbound(r_mangrove)
      if(allocated(sigma_w))write(3,*)'sigma_w',1+ubound(sigma_w)-lbound(sigma_w)
      if(allocated(dsig_f))write(3,*)'dsig_f',1+ubound(dsig_f)-lbound(dsig_f)
      if(allocated(rhp_t))write(3,*)'rhp_t',1+ubound(rhp_t)-lbound(rhp_t)
      if(allocated(rhcref_t))write(3,*)'rhcref_t',1+ubound(rhcref_t)-lbound(rhcref_t)
      if(allocated(vel_u))write(3,*)'vel_u',1+ubound(vel_u)-lbound(vel_u)
      if(allocated(vel_v))write(3,*)'vel_v',1+ubound(vel_v)-lbound(vel_v)
      if(allocated(presgrad_u))write(3,*)'presgrad_u',1+ubound(presgrad_u)-lbound(presgrad_u)
      if(allocated(presgrad_v))write(3,*)'presgrad_v',1+ubound(presgrad_v)-lbound(presgrad_v)
      if(allocated(omega_w))write(3,*)'omega_w',1+ubound(omega_w)-lbound(omega_w)
      if(allocated(veldydz_u))write(3,*)'veldydz_u',1+ubound(veldydz_u)-lbound(veldydz_u)
      if(allocated(veldxdz_v))write(3,*)'veldxdz_v',1+ubound(veldxdz_v)-lbound(veldxdz_v)
      if(allocated(tem_t))write(3,*)'tem_t',1+ubound(tem_t)-lbound(tem_t)
      if(allocated(sal_t))write(3,*)'sal_t',1+ubound(sal_t)-lbound(sal_t)
      if(allocated(tridia_in))write(3,*)'tridia_in',1+ubound(tridia_in)-lbound(tridia_in)
      if(allocated(hr_z2lr_w))write(3,*)'hr_z2lr_w',1+ubound(hr_z2lr_w)-lbound(hr_z2lr_w)
      if(allocated(anyv3d))write(3,*)'anyv3d',1+ubound(anyv3d)-lbound(anyv3d)
      if(allocated(sponge_t))write(3,*)'sponge_t',1+ubound(sponge_t)-lbound(sponge_t)
      if(allocated(sponge_u))write(3,*)'sponge_u',1+ubound(sponge_u)-lbound(sponge_u)
      if(allocated(sponge_v))write(3,*)'sponge_v',1+ubound(sponge_v)-lbound(sponge_v)
      if(allocated(uwind_t))write(3,*)'uwind_t',1+ubound(uwind_t)-lbound(uwind_t)
      if(allocated(vwind_t))write(3,*)'vwind_t',1+ubound(vwind_t)-lbound(vwind_t)
      if(allocated(uwind100_t))write(3,*)'uwind100_t',1+ubound(uwind100_t)-lbound(uwind100_t)
      if(allocated(vwind100_t))write(3,*)'vwind100_t',1+ubound(vwind100_t)-lbound(vwind100_t)
      if(allocated(sshf_w))write(3,*)'sshf_w',1+ubound(sshf_w)-lbound(sshf_w)
      if(allocated(slhf_w))write(3,*)'slhf_w',1+ubound(slhf_w)-lbound(slhf_w)
      if(allocated(ssr_w))write(3,*)'ssr_w',1+ubound(ssr_w)-lbound(ssr_w)
      if(allocated(ssr24prv_w))write(3,*)'ssr24prv_w',1+ubound(ssr24prv_w)-lbound(ssr24prv_w)
      if(allocated(snsf_w))write(3,*)'snsf_w',1+ubound(snsf_w)-lbound(snsf_w)
      if(allocated(precipi_w))write(3,*)'precipi_w',1+ubound(precipi_w)-lbound(precipi_w)
      if(allocated(taux_w))write(3,*)'taux_w',1+ubound(taux_w)-lbound(taux_w)
      if(allocated(tauy_w))write(3,*)'tauy_w',1+ubound(tauy_w)-lbound(tauy_w)
      if(allocated(sigma_f))write(3,*)'sigma_f',1+ubound(sigma_f)-lbound(sigma_f)
      if(allocated(km_w))write(3,*)'km_w',1+ubound(km_w)-lbound(km_w)
      if(allocated(kh_w))write(3,*)'kh_w',1+ubound(kh_w)-lbound(kh_w)
      if(allocated(rho_t))write(3,*)'rho_t',1+ubound(rho_t)-lbound(rho_t)
      if(allocated(omega_evaprec_w))write(3,*)'omega_evaprec_w',1+ubound(omega_evaprec_w)-lbound(omega_evaprec_w)
      if(allocated(tridia_out))write(3,*)'tridia_out',1+ubound(tridia_out)-lbound(tridia_out)
      if(allocated(fluxbar_sumt_u))write(3,*)'fluxbar_sumt_u',1+ubound(fluxbar_sumt_u)-lbound(fluxbar_sumt_u)
      if(allocated(fluxbar_sumt_v))write(3,*)'fluxbar_sumt_v',1+ubound(fluxbar_sumt_v)-lbound(fluxbar_sumt_v)
      if(allocated(velbar_u))write(3,*)'velbar_u',1+ubound(velbar_u)-lbound(velbar_u)
      if(allocated(uflux2d_f))write(3,*)'uflux2d_f',1+ubound(uflux2d_f)-lbound(uflux2d_f)
      if(allocated(vflux2d_f))write(3,*)'vflux2d_f',1+ubound(vflux2d_f)-lbound(vflux2d_f)
      if(allocated(vlxbar_f))write(3,*)'vlxbar_f',1+ubound(vlxbar_f)-lbound(vlxbar_f)
      if(allocated(velbar_v))write(3,*)'velbar_v',1+ubound(velbar_v)-lbound(velbar_v)
      if(allocated(vlybar_f))write(3,*)'vlybar_f',1+ubound(vlybar_f)-lbound(vlybar_f)
      if(allocated(flux2d_u))write(3,*)'flux2d_u',1+ubound(flux2d_u)-lbound(flux2d_u)
      if(allocated(flux2d_v))write(3,*)'flux2d_v',1+ubound(flux2d_v)-lbound(flux2d_v)
      if(allocated(ssh_w))write(3,*)'ssh_w',1+ubound(ssh_w)-lbound(ssh_w)
      if(allocated(hssh_f))write(3,*)'hssh_f',1+ubound(hssh_f)-lbound(hssh_f)
      if(allocated(hz_u))write(3,*)'hz_u',1+ubound(hz_u)-lbound(hz_u)
      if(allocated(hz_v))write(3,*)'hz_v',1+ubound(hz_v)-lbound(hz_v)
      if(allocated(hz_w))write(3,*)'hz_w',1+ubound(hz_w)-lbound(hz_w)
      if(allocated(heatrelax_w))write(3,*)'heatrelax_w',1+ubound(heatrelax_w)-lbound(heatrelax_w)
      if(allocated(xy_u))write(3,*)'xy_u',1+ubound(xy_u)-lbound(xy_u)
      if(allocated(xy_v))write(3,*)'xy_v',1+ubound(xy_v)-lbound(xy_v)
      if(allocated(xy_t))write(3,*)'xy_t',1+ubound(xy_t)-lbound(xy_t)
      if(allocated(xy_f))write(3,*)'xy_f',1+ubound(xy_f)-lbound(xy_f)
      if(allocated(pss_w))write(3,*)'pss_w',1+ubound(pss_w)-lbound(pss_w)
      if(allocated(wstress_u))write(3,*)'wstress_u',1+ubound(wstress_u)-lbound(wstress_u)
      if(allocated(wstress_v))write(3,*)'wstress_v',1+ubound(wstress_v)-lbound(wstress_v)
      if(allocated(ssh_int_w))write(3,*)'ssh_int_w',1+ubound(ssh_int_w)-lbound(ssh_int_w)
      if(allocated(ustokesvortex_t))write(3,*)'ustokesvortex_t',1+ubound(ustokesvortex_t)-lbound(ustokesvortex_t)
      if(allocated(ustokesvortex_f))write(3,*)'ustokesvortex_f',1+ubound(ustokesvortex_f)-lbound(ustokesvortex_f)
      if(allocated(vstokesvortex_t))write(3,*)'vstokesvortex_t',1+ubound(vstokesvortex_t)-lbound(vstokesvortex_t)
      if(allocated(vstokesvortex_f))write(3,*)'vstokesvortex_f',1+ubound(vstokesvortex_f)-lbound(vstokesvortex_f)
      if(allocated(velavr_u))write(3,*)'velavr_u',1+ubound(velavr_u)-lbound(velavr_u)
      if(allocated(velavr_v))write(3,*)'velavr_v',1+ubound(velavr_v)-lbound(velavr_v)
      if(allocated(timefilter_u))write(3,*)'timefilter_u',1+ubound(timefilter_u)-lbound(timefilter_u)
      if(allocated(timefilter_v))write(3,*)'timefilter_v',1+ubound(timefilter_v)-lbound(timefilter_v)
      if(allocated(teta0_t))write(3,*)'teta0_t',1+ubound(teta0_t)-lbound(teta0_t)
      if(allocated(teta2_t))write(3,*)'teta2_t',1+ubound(teta2_t)-lbound(teta2_t)
      if(allocated(q2_t))write(3,*)'q2_t',1+ubound(q2_t)-lbound(q2_t)
      if(allocated(teta2delta_t))write(3,*)'teta2delta_t',1+ubound(teta2delta_t)-lbound(teta2delta_t)
      if(allocated(q2delta_t))write(3,*)'q2delta_t',1+ubound(q2delta_t)-lbound(q2delta_t)
      if(allocated(uwinddelta_t))write(3,*)'uwinddelta_t',1+ubound(uwinddelta_t)-lbound(uwinddelta_t)
      if(allocated(vwinddelta_t))write(3,*)'vwinddelta_t',1+ubound(vwinddelta_t)-lbound(vwinddelta_t)
      if(allocated(tetastar_t))write(3,*)'tetastar_t',1+ubound(tetastar_t)-lbound(tetastar_t)
      if(allocated(ustar_t))write(3,*)'ustar_t',1+ubound(ustar_t)-lbound(ustar_t)
      if(allocated(qstar_t))write(3,*)'qstar_t',1+ubound(qstar_t)-lbound(qstar_t)
      if(allocated(frozenterm2d_u))write(3,*)'frozenterm2d_u',1+ubound(frozenterm2d_u)-lbound(frozenterm2d_u)
      if(allocated(frozenterm2d_v))write(3,*)'frozenterm2d_v',1+ubound(frozenterm2d_v)-lbound(frozenterm2d_v)
      if(allocated(frozenterm3d_u))write(3,*)'frozenterm3d_u',1+ubound(frozenterm3d_u)-lbound(frozenterm3d_u)
      if(allocated(frozenterm3d_v))write(3,*)'frozenterm3d_v',1+ubound(frozenterm3d_v)-lbound(frozenterm3d_v)
      if(allocated(fluxbar_u))write(3,*)'fluxbar_u',1+ubound(fluxbar_u)-lbound(fluxbar_u)
      if(allocated(fluxbar_v))write(3,*)'fluxbar_v',1+ubound(fluxbar_v)-lbound(fluxbar_v)
      if(allocated(lon_t))write(3,*)'lon_t',1+ubound(lon_t)-lbound(lon_t)
      if(allocated(lat_t))write(3,*)'lat_t',1+ubound(lat_t)-lbound(lat_t)
      if(allocated(lon_u))write(3,*)'lon_u',1+ubound(lon_u)-lbound(lon_u)
      if(allocated(lat_u))write(3,*)'lat_u',1+ubound(lat_u)-lbound(lat_u)
      if(allocated(lon_v))write(3,*)'lon_v',1+ubound(lon_v)-lbound(lon_v)
      if(allocated(lat_v))write(3,*)'lat_v',1+ubound(lat_v)-lbound(lat_v)
      if(allocated(lon_f))write(3,*)'lon_f',1+ubound(lon_f)-lbound(lon_f)
      if(allocated(lat_f))write(3,*)'lat_f',1+ubound(lat_f)-lbound(lat_f)
      if(allocated(globlon_t))write(3,*)'globlon_t',1+ubound(globlon_t)-lbound(globlon_t)
      if(allocated(globlat_t))write(3,*)'globlat_t',1+ubound(globlat_t)-lbound(globlat_t)
      if(allocated(globlon_u))write(3,*)'globlon_u',1+ubound(globlon_u)-lbound(globlon_u)
      if(allocated(globlat_u))write(3,*)'globlat_u',1+ubound(globlat_u)-lbound(globlat_u)
      if(allocated(globlon_v))write(3,*)'globlon_v',1+ubound(globlon_v)-lbound(globlon_v)
      if(allocated(globlat_v))write(3,*)'globlat_v',1+ubound(globlat_v)-lbound(globlat_v)
      if(ioffline_prv>=1.and.allocated(fluxbarsum_ofl_u))write(3,*)'fluxbarsum_ofl_u',1+ubound(fluxbarsum_ofl_u)-lbound(fluxbarsum_ofl_u)
      if(ioffline_prv>=1.and.allocated(fluxbarsum_ofl_v))write(3,*)'fluxbarsum_ofl_v',1+ubound(fluxbarsum_ofl_v)-lbound(fluxbarsum_ofl_v)
      if(allocated(dxdy_u))write(3,*)'dxdy_u',1+ubound(dxdy_u)-lbound(dxdy_u)
      if(allocated(dxdy_v))write(3,*)'dxdy_v',1+ubound(dxdy_v)-lbound(dxdy_v)
      if(allocated(dxdy_t))write(3,*)'dxdy_t',1+ubound(dxdy_t)-lbound(dxdy_t)
      if(allocated(dxdy_e_t))write(3,*)'dxdy_e_t',1+ubound(dxdy_e_t)-lbound(dxdy_e_t)
      if(allocated(dx_e_f))write(3,*)'dx_e_f',1+ubound(dx_e_f)-lbound(dx_e_f)
      if(allocated(dy_e_f))write(3,*)'dy_e_f',1+ubound(dy_e_f)-lbound(dy_e_f)
      if(allocated(dx_u))write(3,*)'dx_u',1+ubound(dx_u)-lbound(dx_u)
      if(allocated(dx_v))write(3,*)'dx_v',1+ubound(dx_v)-lbound(dx_v)
      if(allocated(dx_t))write(3,*)'dx_t',1+ubound(dx_t)-lbound(dx_t)
      if(allocated(dx_f))write(3,*)'dx_f',1+ubound(dx_f)-lbound(dx_f)
      if(allocated(dy_u))write(3,*)'dy_u',1+ubound(dy_u)-lbound(dy_u)
      if(allocated(dy_v))write(3,*)'dy_v',1+ubound(dy_v)-lbound(dy_v)
      if(allocated(dy_t))write(3,*)'dy_t',1+ubound(dy_t)-lbound(dy_t)
      if(allocated(dy_f))write(3,*)'dy_f',1+ubound(dy_f)-lbound(dy_f)
      if(allocated(h_w))write(3,*)'h_w',1+ubound(h_w)-lbound(h_w)
      if(allocated(hcopy_w))write(3,*)'hcopy_w',1+ubound(hcopy_w)-lbound(hcopy_w)
      if(allocated(h0_w))write(3,*)'h0_w',1+ubound(h0_w)-lbound(h0_w)
      if(allocated(h_u))write(3,*)'h_u',1+ubound(h_u)-lbound(h_u)
      if(allocated(h_v))write(3,*)'h_v',1+ubound(h_v)-lbound(h_v)
      if(allocated(h_f))write(3,*)'h_f',1+ubound(h_f)-lbound(h_f)
      if(allocated(coriolis_t))write(3,*)'coriolis_t',1+ubound(coriolis_t)-lbound(coriolis_t)
      if(allocated(coriolis_f))write(3,*)'coriolis_f',1+ubound(coriolis_f)-lbound(coriolis_f)
      if(allocated(rhpzavr_w))write(3,*)'rhpzavr_w',1+ubound(rhpzavr_w)-lbound(rhpzavr_w)
      if(allocated(q10_t))write(3,*)'q10_t',1+ubound(q10_t)-lbound(q10_t)
      if(allocated(teta10_t))write(3,*)'teta10_t',1+ubound(teta10_t)-lbound(teta10_t)
      if(allocated(fric_u))write(3,*)'fric_u',1+ubound(fric_u)-lbound(fric_u)
      if(allocated(fric_v))write(3,*)'fric_v',1+ubound(fric_v)-lbound(fric_v)
      if(allocated(fric_t))write(3,*)'fric_t',1+ubound(fric_t)-lbound(fric_t)
      if(allocated(cdb_t))write(3,*)'cdb_t',1+ubound(cdb_t)-lbound(cdb_t)
      if(allocated(cdb_f))write(3,*)'cdb_f',1+ubound(cdb_f)-lbound(cdb_f)
      if(allocated(xflux_t))write(3,*)'xflux_t',1+ubound(xflux_t)-lbound(xflux_t)
      if(allocated(xflux_f))write(3,*)'xflux_f',1+ubound(xflux_f)-lbound(xflux_f)
      if(allocated(yflux_t))write(3,*)'yflux_t',1+ubound(yflux_t)-lbound(yflux_t)
      if(allocated(yflux_f))write(3,*)'yflux_f',1+ubound(yflux_f)-lbound(yflux_f)
      if(allocated(pres3d2d_u))write(3,*)'pres3d2d_u',1+ubound(pres3d2d_u)-lbound(pres3d2d_u)
      if(allocated(pres3d2d_v))write(3,*)'pres3d2d_v',1+ubound(pres3d2d_v)-lbound(pres3d2d_v)
      if(allocated(adve3d2d_u))write(3,*)'adve3d2d_u',1+ubound(adve3d2d_u)-lbound(adve3d2d_u)
      if(allocated(adve3d2d_v))write(3,*)'adve3d2d_v',1+ubound(adve3d2d_v)-lbound(adve3d2d_v)
      if(allocated(rmangrovebar))write(3,*)'rmangrovebar',1+ubound(rmangrovebar)-lbound(rmangrovebar)
      if(allocated(mang3dto2d_u))write(3,*)'mang3dto2d_u',1+ubound(mang3dto2d_u)-lbound(mang3dto2d_u)
      if(allocated(mang3dto2d_v))write(3,*)'mang3dto2d_v',1+ubound(mang3dto2d_v)-lbound(mang3dto2d_v)
      if(allocated(restoring3d2d_u))write(3,*)'restoring3d2d_u',1+ubound(restoring3d2d_u)-lbound(restoring3d2d_u)
      if(allocated(restoring3d2d_v))write(3,*)'restoring3d2d_v',1+ubound(restoring3d2d_v)-lbound(restoring3d2d_v)
      if(allocated(stokesforces3d2d_u))write(3,*)'stokesforces3d2d_u',1+ubound(stokesforces3d2d_u)-lbound(stokesforces3d2d_u)
      if(allocated(stokesforces3d2d_v))write(3,*)'stokesforces3d2d_v',1+ubound(stokesforces3d2d_v)-lbound(stokesforces3d2d_v)
      if(allocated(wstress_w))write(3,*)'wstress_w',1+ubound(wstress_w)-lbound(wstress_w)
      if(allocated(z0_w))write(3,*)'z0_w',1+ubound(z0_w)-lbound(z0_w)
      if(allocated(albedo_w))write(3,*)'albedo_w',1+ubound(albedo_w)-lbound(albedo_w)
      if(allocated(gridrotcos_t))write(3,*)'gridrotcos_t',1+ubound(gridrotcos_t)-lbound(gridrotcos_t)
      if(allocated(gridrotsin_t))write(3,*)'gridrotsin_t',1+ubound(gridrotsin_t)-lbound(gridrotsin_t)
      if(allocated(grid_angle_t))write(3,*)'grid_angle_t',1+ubound(grid_angle_t)-lbound(grid_angle_t)
      if(allocated(sshstokes_w))write(3,*)'sshstokes_w',1+ubound(sshstokes_w)-lbound(sshstokes_w)
      if(allocated(cwi_int_u))write(3,*)'cwi_int_u',1+ubound(cwi_int_u)-lbound(cwi_int_u)
      if(allocated(cwi_int_v))write(3,*)'cwi_int_v',1+ubound(cwi_int_v)-lbound(cwi_int_v)
      if(allocated(cwj_int_u))write(3,*)'cwj_int_u',1+ubound(cwj_int_u)-lbound(cwj_int_u)
      if(allocated(cwj_int_v))write(3,*)'cwj_int_v',1+ubound(cwj_int_v)-lbound(cwj_int_v)
      if(allocated(sshrefobc_i))write(3,*)'sshrefobc_i',1+ubound(sshrefobc_i)-lbound(sshrefobc_i)
      if(allocated(vbrrefobc_i))write(3,*)'vbrrefobc_i',1+ubound(vbrrefobc_i)-lbound(vbrrefobc_i)
      if(allocated(sshrefobc_j))write(3,*)'sshrefobc_j',1+ubound(sshrefobc_j)-lbound(sshrefobc_j)
      if(allocated(vbrrefobc_j))write(3,*)'vbrrefobc_j',1+ubound(vbrrefobc_j)-lbound(vbrrefobc_j)
      if(allocated(anyv1d))write(3,*)'anyv1d',1+ubound(anyv1d)-lbound(anyv1d)
      if(allocated(gridcard))write(3,*)'gridcard',1+ubound(gridcard)-lbound(gridcard)
      if(allocated(novector))write(3,*)'novector',1+ubound(novector)-lbound(novector)
      if(allocated(proficard))write(3,*)'proficard',1+ubound(proficard)-lbound(proficard)
      if(allocated(airseainfo))write(3,*)'airseainfo',1+ubound(airseainfo)-lbound(airseainfo)
      if(allocated(airseadt))write(3,*)'airseadt',1+ubound(airseadt)-lbound(airseadt)
      if(allocated(riverdt))write(3,*)'riverdt',1+ubound(riverdt)-lbound(riverdt)
      if(allocated(river_t))write(3,*)'river_t',1+ubound(river_t)-lbound(river_t)
      if(allocated(riverflux))write(3,*)'riverflux',1+ubound(riverflux)-lbound(riverflux)
      if(allocated(mask_t))write(3,*)'mask_t',1+ubound(mask_t)-lbound(mask_t)
      if(allocated(mask_f))write(3,*)'mask_f',1+ubound(mask_f)-lbound(mask_f)
      if(allocated(mask_u))write(3,*)'mask_u',1+ubound(mask_u)-lbound(mask_u)
      if(allocated(mask_v))write(3,*)'mask_v',1+ubound(mask_v)-lbound(mask_v)
      if(allocated(mask_vqs_tke_w))write(3,*)'mask_vqs_tke_w',1+ubound(mask_vqs_tke_w)-lbound(mask_vqs_tke_w)
      if(allocated(canaldir))write(3,*)'canaldir',1+ubound(canaldir)-lbound(canaldir)
      if(allocated(wetmask_wi_t))write(3,*)'wetmask_wi_t',1+ubound(wetmask_wi_t)-lbound(wetmask_wi_t)
      if(allocated(lonlat2ij_t))write(3,*)'lonlat2ij_t',1+ubound(lonlat2ij_t)-lbound(lonlat2ij_t)
      if(allocated(canalmpioverlap))write(3,*)'canalmpioverlap',1+ubound(canalmpioverlap)-lbound(canalmpioverlap)
      if(allocated(sodate))write(3,*)'sodate',1+ubound(sodate)-lbound(sodate)
      if(allocated(canalrank))write(3,*)'canalrank',1+ubound(canalrank)-lbound(canalrank)
      if(allocated(canalrankbis))write(3,*)'canalrankbis',1+ubound(canalrankbis)-lbound(canalrankbis)
      if(allocated(i_canalcoord))write(3,*)'i_canalcoord',1+ubound(i_canalcoord)-lbound(i_canalcoord)
      if(allocated(j_canalcoord))write(3,*)'j_canalcoord',1+ubound(j_canalcoord)-lbound(j_canalcoord)
      if(allocated(kmin_u))write(3,*)'kmin_u',1+ubound(kmin_u)-lbound(kmin_u)
      if(allocated(kmin_v))write(3,*)'kmin_v',1+ubound(kmin_v)-lbound(kmin_v)
      if(allocated(kmin_w))write(3,*)'kmin_w',1+ubound(kmin_w)-lbound(kmin_w)
      if(allocated(kundermin_t))write(3,*)'kundermin_t',1+ubound(kundermin_t)-lbound(kundermin_t)
      if(allocated(kundermin_u))write(3,*)'kundermin_u',1+ubound(kundermin_u)-lbound(kundermin_u)
      if(allocated(kundermin_v))write(3,*)'kundermin_v',1+ubound(kundermin_v)-lbound(kundermin_v)
      if(allocated(kmerged_t))write(3,*)'kmerged_t',1+ubound(kmerged_t)-lbound(kmerged_t)
      if(allocated(kmerged_u))write(3,*)'kmerged_u',1+ubound(kmerged_u)-lbound(kmerged_u)
      if(allocated(kmerged_v))write(3,*)'kmerged_v',1+ubound(kmerged_v)-lbound(kmerged_v)
      if(allocated(ksl_t))write(3,*)'ksl_t',1+ubound(ksl_t)-lbound(ksl_t)
      if(allocated(glob_mask_mangrove))write(3,*)'glob_mask_mangrove',1+ubound(glob_mask_mangrove)-lbound(glob_mask_mangrove)
      if(allocated(mask_mangrove_t))write(3,*)'mask_mangrove_t',1+ubound(mask_mangrove_t)-lbound(mask_mangrove_t)
      if(allocated(mask_wave_t))write(3,*)'mask_wave_t',1+ubound(mask_wave_t)-lbound(mask_wave_t)
      if(allocated(upwindriver_t))write(3,*)'upwindriver_t',1+ubound(upwindriver_t)-lbound(upwindriver_t)
      if(allocated(upwindwetdry_t))write(3,*)'upwindwetdry_t',1+ubound(upwindwetdry_t)-lbound(upwindwetdry_t)
      if(allocated(kmergedr4_u))write(3,*)'kmergedr4_u',1+ubound(kmergedr4_u)-lbound(kmergedr4_u)
      if(allocated(kmergedr4_v))write(3,*)'kmergedr4_v',1+ubound(kmergedr4_v)-lbound(kmergedr4_v)
      if(allocated(dsigmerged_u))write(3,*)'dsigmerged_u',1+ubound(dsigmerged_u)-lbound(dsigmerged_u)
      if(allocated(dsigmerged_v))write(3,*)'dsigmerged_v',1+ubound(dsigmerged_v)-lbound(dsigmerged_v)
      if(allocated(pgfratio_u))write(3,*)'pgfratio_u',1+ubound(pgfratio_u)-lbound(pgfratio_u)
      if(allocated(pgfratio_v))write(3,*)'pgfratio_v',1+ubound(pgfratio_v)-lbound(pgfratio_v)
      if(allocated(sshr4_w))write(3,*)'sshr4_w',1+ubound(sshr4_w)-lbound(sshr4_w)
      if(allocated(wetmask_u))write(3,*)'wetmask_u',1+ubound(wetmask_u)-lbound(wetmask_u)
      if(allocated(wetmask_v))write(3,*)'wetmask_v',1+ubound(wetmask_v)-lbound(wetmask_v)
      if(allocated(wetmask_t))write(3,*)'wetmask_t',1+ubound(wetmask_t)-lbound(wetmask_t)
      if(allocated(botlevmerged_w))write(3,*)'botlevmerged_w',1+ubound(botlevmerged_w)-lbound(botlevmerged_w)
      if(allocated(maxbotstress_aft_w))write(3,*)'maxbotstress_aft_w',1+ubound(maxbotstress_aft_w)-lbound(maxbotstress_aft_w)
      if(allocated(maxbotstress_bef_w))write(3,*)'maxbotstress_bef_w',1+ubound(maxbotstress_bef_w)-lbound(maxbotstress_bef_w)
      if(allocated(maxbotstress_w))write(3,*)'maxbotstress_w',1+ubound(maxbotstress_w)-lbound(maxbotstress_w)
      if(allocated(stresswave_w))write(3,*)'stresswave_w',1+ubound(stresswave_w)-lbound(stresswave_w)
      if(allocated(stressc_w))write(3,*)'stressc_w',1+ubound(stressc_w)-lbound(stressc_w)
      if(allocated(sqr_hoverg_u))write(3,*)'sqr_hoverg_u',1+ubound(sqr_hoverg_u)-lbound(sqr_hoverg_u)
      if(allocated(sqr_hoverg_v))write(3,*)'sqr_hoverg_v',1+ubound(sqr_hoverg_v)-lbound(sqr_hoverg_v)
      if(allocated(temobc_t))write(3,*)'temobc_t',1+ubound(temobc_t)-lbound(temobc_t)
      if(allocated(salobc_t))write(3,*)'salobc_t',1+ubound(salobc_t)-lbound(salobc_t)
      if(allocated(velobc_u))write(3,*)'velobc_u',1+ubound(velobc_u)-lbound(velobc_u)
      if(allocated(velobc_v))write(3,*)'velobc_v',1+ubound(velobc_v)-lbound(velobc_v)
      if(ioffline_prv>=1.and.allocated(temofl_t))write(3,*)'temofl_t',1+ubound(temofl_t)-lbound(temofl_t)
      if(ioffline_prv>=1.and.allocated(salofl_t))write(3,*)'salofl_t',1+ubound(salofl_t)-lbound(salofl_t)
      if(ioffline_prv>=1.and.allocated(dzofl_t))write(3,*)'dzofl_t',1+ubound(dzofl_t)-lbound(dzofl_t)
      if(ioffline_prv>=1.and.allocated(velofl_u))write(3,*)'velofl_u',1+ubound(velofl_u)-lbound(velofl_u)
      if(ioffline_prv>=1.and.allocated(velofl_v))write(3,*)'velofl_v',1+ubound(velofl_v)-lbound(velofl_v)
      if(ioffline_prv>=1.and.allocated(tkeofl_w))write(3,*)'tkeofl_w',1+ubound(tkeofl_w)-lbound(tkeofl_w)
      if(ioffline_prv>=1.and.allocated(bioofl_t))write(3,*)'bioofl_t',1+ubound(bioofl_t)-lbound(bioofl_t)
      if(ioffline_prv>=1.and.allocated(dfvofl_w))write(3,*)'dfvofl_w',1+ubound(dfvofl_w)-lbound(dfvofl_w)
      if(allocated(uwindabl_t))write(3,*)'uwindabl_t',1+ubound(uwindabl_t)-lbound(uwindabl_t)
      if(allocated(vwindabl_t))write(3,*)'vwindabl_t',1+ubound(vwindabl_t)-lbound(vwindabl_t)
      if(ioffline_prv>=1.and.allocated(w0mofl_w))write(3,*)'w0mofl_w',1+ubound(w0mofl_w)-lbound(w0mofl_w)
      if(ioffline_prv>=1.and.allocated(w_keq1_ofl_w))write(3,*)'w_keq1_ofl_w',1+ubound(w_keq1_ofl_w)-lbound(w_keq1_ofl_w)
      if(ioffline_prv>=1.and.allocated(kslofl_t))write(3,*)'kslofl_t',1+ubound(kslofl_t)-lbound(kslofl_t)
      if(allocated(ablheight_t))write(3,*)'ablheight_t',1+ubound(ablheight_t)-lbound(ablheight_t)
      if(allocated(wwindabl_w))write(3,*)'wwindabl_w',1+ubound(wwindabl_w)-lbound(wwindabl_w)
      if(allocated(kz_abl_w))write(3,*)'kz_abl_w',1+ubound(kz_abl_w)-lbound(kz_abl_w)
      if(allocated(upwzone0_t))write(3,*)'upwzone0_t',1+ubound(upwzone0_t)-lbound(upwzone0_t)
      if(allocated(velstokes_u))write(3,*)'velstokes_u',1+ubound(velstokes_u)-lbound(velstokes_u)
      if(allocated(velstokes_v))write(3,*)'velstokes_v',1+ubound(velstokes_v)-lbound(velstokes_v)
      if(allocated(velbarstokes_u))write(3,*)'velbarstokes_u',1+ubound(velbarstokes_u)-lbound(velbarstokes_u)
      if(allocated(velbarstokes_v))write(3,*)'velbarstokes_v',1+ubound(velbarstokes_v)-lbound(velbarstokes_v)
      if(allocated(nhp1_t))write(3,*)'nhp1_t',1+ubound(nhp1_t)-lbound(nhp1_t)
      if(allocated(nhp2_t))write(3,*)'nhp2_t',1+ubound(nhp2_t)-lbound(nhp2_t)
      if(allocated(temf_t))write(3,*)'temf_t',1+ubound(temf_t)-lbound(temf_t)
      if(allocated(salf_t))write(3,*)'salf_t',1+ubound(salf_t)-lbound(salf_t)
      if(allocated(temlwf_t))write(3,*)'temlwf_t',1+ubound(temlwf_t)-lbound(temlwf_t)
      if(allocated(sallwf_t))write(3,*)'sallwf_t',1+ubound(sallwf_t)-lbound(sallwf_t)
      if(allocated(sshlwf_w))write(3,*)'sshlwf_w',1+ubound(sshlwf_w)-lbound(sshlwf_w)
      if(allocated(t_wave_t))write(3,*)'t_wave_t',1+ubound(t_wave_t)-lbound(t_wave_t)
      if(allocated(hs_wave_t))write(3,*)'hs_wave_t',1+ubound(hs_wave_t)-lbound(hs_wave_t)
      if(allocated(hsw_wave_t))write(3,*)'hsw_wave_t',1+ubound(hsw_wave_t)-lbound(hsw_wave_t)
      if(allocated(foc_wave_t))write(3,*)'foc_wave_t',1+ubound(foc_wave_t)-lbound(foc_wave_t)
      if(allocated(k_wave_t))write(3,*)'k_wave_t',1+ubound(k_wave_t)-lbound(k_wave_t)
      if(allocated(kx_wave_t))write(3,*)'kx_wave_t',1+ubound(kx_wave_t)-lbound(kx_wave_t)
      if(allocated(ky_wave_t))write(3,*)'ky_wave_t',1+ubound(ky_wave_t)-lbound(ky_wave_t)
      if(allocated(twox_wave_t))write(3,*)'twox_wave_t',1+ubound(twox_wave_t)-lbound(twox_wave_t)
      if(allocated(twoy_wave_t))write(3,*)'twoy_wave_t',1+ubound(twoy_wave_t)-lbound(twoy_wave_t)
      if(allocated(tawx_wave_t))write(3,*)'tawx_wave_t',1+ubound(tawx_wave_t)-lbound(tawx_wave_t)
      if(allocated(tawy_wave_t))write(3,*)'tawy_wave_t',1+ubound(tawy_wave_t)-lbound(tawy_wave_t)
      if(allocated(usf_wave_t))write(3,*)'usf_wave_t',1+ubound(usf_wave_t)-lbound(usf_wave_t)
      if(allocated(vsf_wave_t))write(3,*)'vsf_wave_t',1+ubound(vsf_wave_t)-lbound(vsf_wave_t)
      if(allocated(dir_wave_t))write(3,*)'dir_wave_t',1+ubound(dir_wave_t)-lbound(dir_wave_t)
      if(allocated(uss_wave_t))write(3,*)'uss_wave_t',1+ubound(uss_wave_t)-lbound(uss_wave_t)
      if(allocated(j_wave_t))write(3,*)'j_wave_t',1+ubound(j_wave_t)-lbound(j_wave_t)
      if(allocated(vss_wave_t))write(3,*)'vss_wave_t',1+ubound(vss_wave_t)-lbound(vss_wave_t)
      if(allocated(ubw))write(3,*)'ubw',1+ubound(ubw)-lbound(ubw)
      if(allocated(fw))write(3,*)'fw',1+ubound(fw)-lbound(fw)
      if(allocated(dpt_wave_t))write(3,*)'dpt_wave_t',1+ubound(dpt_wave_t)-lbound(dpt_wave_t)
      if(allocated(wstresb_u))write(3,*)'wstresb_u',1+ubound(wstresb_u)-lbound(wstresb_u)
      if(allocated(wstresb_v))write(3,*)'wstresb_v',1+ubound(wstresb_v)-lbound(wstresb_v)
      if(allocated(ij2ww3_i))write(3,*)'ij2ww3_i',1+ubound(ij2ww3_i)-lbound(ij2ww3_i)
      if(allocated(ij2ww3_j))write(3,*)'ij2ww3_j',1+ubound(ij2ww3_j)-lbound(ij2ww3_j)
      if(allocated(ij2ww3_teta))write(3,*)'ij2ww3_teta',1+ubound(ij2ww3_teta)-lbound(ij2ww3_teta)
      if(allocated(slhf_aver_w))write(3,*)'slhf_aver_w',1+ubound(slhf_aver_w)-lbound(slhf_aver_w)
      if(allocated(sshf_aver_w))write(3,*)'sshf_aver_w',1+ubound(sshf_aver_w)-lbound(sshf_aver_w)
      if(allocated(snsf_aver_w))write(3,*)'snsf_aver_w',1+ubound(snsf_aver_w)-lbound(snsf_aver_w)
      if(allocated(ssr_aver_w))write(3,*)'ssr_aver_w',1+ubound(ssr_aver_w)-lbound(ssr_aver_w)
      if(allocated(precipi_aver_w))write(3,*)'precipi_aver_w',1+ubound(precipi_aver_w)-lbound(precipi_aver_w)
      if(allocated(wstress_aver_u))write(3,*)'wstress_aver_u',1+ubound(wstress_aver_u)-lbound(wstress_aver_u)
      if(allocated(wstress_aver_v))write(3,*)'wstress_aver_v',1+ubound(wstress_aver_v)-lbound(wstress_aver_v)
      if(allocated(hsedofl_t))write(3,*)'hsedofl_t',1+ubound(hsedofl_t)-lbound(hsedofl_t)
      write(3,*)'zone1saltflux_w'
      write(3,*)'zone1tempflux_w'
      write(3,*)'zone1waterflux_w'
      write(3,*)'zone1_nlayer'
      write(3,*)'zone1_max'
      write(3,*)'zone1_u_max'
      write(3,*)'zone1_v_max'
      write(3,*)'zone1_inv_dz'
      write(3,*)'zone1_stretch_dz'
      if(allocated(zone1saltflux_glb))write(3,*)'zone1saltflux_glb',1+ubound(zone1saltflux_glb)-lbound(zone1saltflux_glb)
      if(allocated(zone1saltflux_u))write(3,*)'zone1saltflux_u',1+ubound(zone1saltflux_u)-lbound(zone1saltflux_u)
      if(allocated(zone1saltflux_v))write(3,*)'zone1saltflux_v',1+ubound(zone1saltflux_v)-lbound(zone1saltflux_v)
      if(allocated(zone1tempflux_glb))write(3,*)'zone1tempflux_glb',1+ubound(zone1tempflux_glb)-lbound(zone1tempflux_glb)
      if(allocated(zone1tempflux_u))write(3,*)'zone1tempflux_u',1+ubound(zone1tempflux_u)-lbound(zone1tempflux_u)
      if(allocated(zone1tempflux_v))write(3,*)'zone1tempflux_v',1+ubound(zone1tempflux_v)-lbound(zone1tempflux_v)
      if(allocated(zone1waterflux_glb))write(3,*)'zone1waterflux_glb',1+ubound(zone1waterflux_glb)-lbound(zone1waterflux_glb)
      if(allocated(zone1waterflux_u))write(3,*)'zone1waterflux_u',1+ubound(zone1waterflux_u)-lbound(zone1waterflux_u)
      if(allocated(zone1waterflux_v))write(3,*)'zone1waterflux_v',1+ubound(zone1waterflux_v)-lbound(zone1waterflux_v)
      if(allocated(zone1saltflux_glb_in))write(3,*)'zone1saltflux_glb_in',1+ubound(zone1saltflux_glb_in)-lbound(zone1saltflux_glb_in)
      if(allocated(zone1saltflux_u_in))write(3,*)'zone1saltflux_u_in',1+ubound(zone1saltflux_u_in)-lbound(zone1saltflux_u_in)
      if(allocated(zone1saltflux_v_in))write(3,*)'zone1saltflux_v_in',1+ubound(zone1saltflux_v_in)-lbound(zone1saltflux_v_in)
      if(allocated(zone1tempflux_glb_in))write(3,*)'zone1tempflux_glb_in',1+ubound(zone1tempflux_glb_in)-lbound(zone1tempflux_glb_in)
      if(allocated(zone1tempflux_u_in))write(3,*)'zone1tempflux_u_in',1+ubound(zone1tempflux_u_in)-lbound(zone1tempflux_u_in)
      if(allocated(zone1tempflux_v_in))write(3,*)'zone1tempflux_v_in',1+ubound(zone1tempflux_v_in)-lbound(zone1tempflux_v_in)
      if(allocated(zone1waterflux_glb_in))write(3,*)'zone1waterflux_glb_in',1+ubound(zone1waterflux_glb_in)-lbound(zone1waterflux_glb_in)
      if(allocated(zone1waterflux_u_in))write(3,*)'zone1waterflux_u_in',1+ubound(zone1waterflux_u_in)-lbound(zone1waterflux_u_in)
      if(allocated(zone1waterflux_v_in))write(3,*)'zone1waterflux_v_in',1+ubound(zone1waterflux_v_in)-lbound(zone1waterflux_v_in)
      if(allocated(zone1saltflux_glb_out))write(3,*)'zone1saltflux_glb_out',1+ubound(zone1saltflux_glb_out)-lbound(zone1saltflux_glb_out)
      if(allocated(zone1saltflux_u_out))write(3,*)'zone1saltflux_u_out',1+ubound(zone1saltflux_u_out)-lbound(zone1saltflux_u_out)
      if(allocated(zone1saltflux_v_out))write(3,*)'zone1saltflux_v_out',1+ubound(zone1saltflux_v_out)-lbound(zone1saltflux_v_out)
      if(allocated(zone1tempflux_glb_out))write(3,*)'zone1tempflux_glb_out',1+ubound(zone1tempflux_glb_out)-lbound(zone1tempflux_glb_out)
      if(allocated(zone1tempflux_u_out))write(3,*)'zone1tempflux_u_out',1+ubound(zone1tempflux_u_out)-lbound(zone1tempflux_u_out)
      if(allocated(zone1tempflux_v_out))write(3,*)'zone1tempflux_v_out',1+ubound(zone1tempflux_v_out)-lbound(zone1tempflux_v_out)
      if(allocated(zone1waterflux_glb_out))write(3,*)'zone1waterflux_glb_out',1+ubound(zone1waterflux_glb_out)-lbound(zone1waterflux_glb_out)
      if(allocated(zone1waterflux_u_out))write(3,*)'zone1waterflux_u_out',1+ubound(zone1waterflux_u_out)-lbound(zone1waterflux_u_out)
      if(allocated(zone1waterflux_v_out))write(3,*)'zone1waterflux_v_out',1+ubound(zone1waterflux_v_out)-lbound(zone1waterflux_v_out)
      if(allocated(zone1_mask))write(3,*)'zone1_mask',1+ubound(zone1_mask)-lbound(zone1_mask)
      if(allocated(zone1_flux_u_node))write(3,*)'zone1_flux_u_node',1+ubound(zone1_flux_u_node)-lbound(zone1_flux_u_node)
      if(allocated(zone1_flux_v_node))write(3,*)'zone1_flux_v_node',1+ubound(zone1_flux_v_node)-lbound(zone1_flux_v_node)
      write(3,*)'zone2saltflux_w'
      write(3,*)'zone2tempflux_w'
      write(3,*)'zone2waterflux_w'
      write(3,*)'zone2_nlayer'
      write(3,*)'zone2_max'
      write(3,*)'zone2_u_max'
      write(3,*)'zone2_v_max'
      write(3,*)'zone2_inv_dz'
      write(3,*)'zone2_stretch_dz'
      if(allocated(zone2saltflux_glb))write(3,*)'zone2saltflux_glb',1+ubound(zone2saltflux_glb)-lbound(zone2saltflux_glb)
      if(allocated(zone2saltflux_u))write(3,*)'zone2saltflux_u',1+ubound(zone2saltflux_u)-lbound(zone2saltflux_u)
      if(allocated(zone2saltflux_v))write(3,*)'zone2saltflux_v',1+ubound(zone2saltflux_v)-lbound(zone2saltflux_v)
      if(allocated(zone2tempflux_glb))write(3,*)'zone2tempflux_glb',1+ubound(zone2tempflux_glb)-lbound(zone2tempflux_glb)
      if(allocated(zone2tempflux_u))write(3,*)'zone2tempflux_u',1+ubound(zone2tempflux_u)-lbound(zone2tempflux_u)
      if(allocated(zone2tempflux_v))write(3,*)'zone2tempflux_v',1+ubound(zone2tempflux_v)-lbound(zone2tempflux_v)
      if(allocated(zone2waterflux_glb))write(3,*)'zone2waterflux_glb',1+ubound(zone2waterflux_glb)-lbound(zone2waterflux_glb)
      if(allocated(zone2waterflux_u))write(3,*)'zone2waterflux_u',1+ubound(zone2waterflux_u)-lbound(zone2waterflux_u)
      if(allocated(zone2waterflux_v))write(3,*)'zone2waterflux_v',1+ubound(zone2waterflux_v)-lbound(zone2waterflux_v)
      if(allocated(zone2saltflux_glb_in))write(3,*)'zone2saltflux_glb_in',1+ubound(zone2saltflux_glb_in)-lbound(zone2saltflux_glb_in)
      if(allocated(zone2saltflux_u_in))write(3,*)'zone2saltflux_u_in',1+ubound(zone2saltflux_u_in)-lbound(zone2saltflux_u_in)
      if(allocated(zone2saltflux_v_in))write(3,*)'zone2saltflux_v_in',1+ubound(zone2saltflux_v_in)-lbound(zone2saltflux_v_in)
      if(allocated(zone2tempflux_glb_in))write(3,*)'zone2tempflux_glb_in',1+ubound(zone2tempflux_glb_in)-lbound(zone2tempflux_glb_in)
      if(allocated(zone2tempflux_u_in))write(3,*)'zone2tempflux_u_in',1+ubound(zone2tempflux_u_in)-lbound(zone2tempflux_u_in)
      if(allocated(zone2tempflux_v_in))write(3,*)'zone2tempflux_v_in',1+ubound(zone2tempflux_v_in)-lbound(zone2tempflux_v_in)
      if(allocated(zone2waterflux_glb_in))write(3,*)'zone2waterflux_glb_in',1+ubound(zone2waterflux_glb_in)-lbound(zone2waterflux_glb_in)
      if(allocated(zone2waterflux_u_in))write(3,*)'zone2waterflux_u_in',1+ubound(zone2waterflux_u_in)-lbound(zone2waterflux_u_in)
      if(allocated(zone2waterflux_v_in))write(3,*)'zone2waterflux_v_in',1+ubound(zone2waterflux_v_in)-lbound(zone2waterflux_v_in)
      if(allocated(zone2saltflux_glb_out))write(3,*)'zone2saltflux_glb_out',1+ubound(zone2saltflux_glb_out)-lbound(zone2saltflux_glb_out)
      if(allocated(zone2saltflux_u_out))write(3,*)'zone2saltflux_u_out',1+ubound(zone2saltflux_u_out)-lbound(zone2saltflux_u_out)
      if(allocated(zone2saltflux_v_out))write(3,*)'zone2saltflux_v_out',1+ubound(zone2saltflux_v_out)-lbound(zone2saltflux_v_out)
      if(allocated(zone2tempflux_glb_out))write(3,*)'zone2tempflux_glb_out',1+ubound(zone2tempflux_glb_out)-lbound(zone2tempflux_glb_out)
      if(allocated(zone2tempflux_u_out))write(3,*)'zone2tempflux_u_out',1+ubound(zone2tempflux_u_out)-lbound(zone2tempflux_u_out)
      if(allocated(zone2tempflux_v_out))write(3,*)'zone2tempflux_v_out',1+ubound(zone2tempflux_v_out)-lbound(zone2tempflux_v_out)
      if(allocated(zone2waterflux_glb_out))write(3,*)'zone2waterflux_glb_out',1+ubound(zone2waterflux_glb_out)-lbound(zone2waterflux_glb_out)
      if(allocated(zone2waterflux_u_out))write(3,*)'zone2waterflux_u_out',1+ubound(zone2waterflux_u_out)-lbound(zone2waterflux_u_out)
      if(allocated(zone2waterflux_v_out))write(3,*)'zone2waterflux_v_out',1+ubound(zone2waterflux_v_out)-lbound(zone2waterflux_v_out)
      if(allocated(zone2_mask))write(3,*)'zone2_mask',1+ubound(zone2_mask)-lbound(zone2_mask)
      if(allocated(zone2_flux_u_node))write(3,*)'zone2_flux_u_node',1+ubound(zone2_flux_u_node)-lbound(zone2_flux_u_node)
      if(allocated(zone2_flux_v_node))write(3,*)'zone2_flux_v_node',1+ubound(zone2_flux_v_node)-lbound(zone2_flux_v_node)
      write(3,*)'zone3saltflux_w'
      write(3,*)'zone3tempflux_w'
      write(3,*)'zone3waterflux_w'
      write(3,*)'zone3_nlayer'
      write(3,*)'zone3_max'
      write(3,*)'zone3_u_max'
      write(3,*)'zone3_v_max'
      write(3,*)'zone3_inv_dz'
      write(3,*)'zone3_stretch_dz'
      if(allocated(zone3saltflux_glb))write(3,*)'zone3saltflux_glb',1+ubound(zone3saltflux_glb)-lbound(zone3saltflux_glb)
      if(allocated(zone3saltflux_u))write(3,*)'zone3saltflux_u',1+ubound(zone3saltflux_u)-lbound(zone3saltflux_u)
      if(allocated(zone3saltflux_v))write(3,*)'zone3saltflux_v',1+ubound(zone3saltflux_v)-lbound(zone3saltflux_v)
      if(allocated(zone3tempflux_glb))write(3,*)'zone3tempflux_glb',1+ubound(zone3tempflux_glb)-lbound(zone3tempflux_glb)
      if(allocated(zone3tempflux_u))write(3,*)'zone3tempflux_u',1+ubound(zone3tempflux_u)-lbound(zone3tempflux_u)
      if(allocated(zone3tempflux_v))write(3,*)'zone3tempflux_v',1+ubound(zone3tempflux_v)-lbound(zone3tempflux_v)
      if(allocated(zone3waterflux_glb))write(3,*)'zone3waterflux_glb',1+ubound(zone3waterflux_glb)-lbound(zone3waterflux_glb)
      if(allocated(zone3waterflux_u))write(3,*)'zone3waterflux_u',1+ubound(zone3waterflux_u)-lbound(zone3waterflux_u)
      if(allocated(zone3waterflux_v))write(3,*)'zone3waterflux_v',1+ubound(zone3waterflux_v)-lbound(zone3waterflux_v)
      if(allocated(zone3saltflux_glb_in))write(3,*)'zone3saltflux_glb_in',1+ubound(zone3saltflux_glb_in)-lbound(zone3saltflux_glb_in)
      if(allocated(zone3saltflux_u_in))write(3,*)'zone3saltflux_u_in',1+ubound(zone3saltflux_u_in)-lbound(zone3saltflux_u_in)
      if(allocated(zone3saltflux_v_in))write(3,*)'zone3saltflux_v_in',1+ubound(zone3saltflux_v_in)-lbound(zone3saltflux_v_in)
      if(allocated(zone3tempflux_glb_in))write(3,*)'zone3tempflux_glb_in',1+ubound(zone3tempflux_glb_in)-lbound(zone3tempflux_glb_in)
      if(allocated(zone3tempflux_u_in))write(3,*)'zone3tempflux_u_in',1+ubound(zone3tempflux_u_in)-lbound(zone3tempflux_u_in)
      if(allocated(zone3tempflux_v_in))write(3,*)'zone3tempflux_v_in',1+ubound(zone3tempflux_v_in)-lbound(zone3tempflux_v_in)
      if(allocated(zone3waterflux_glb_in))write(3,*)'zone3waterflux_glb_in',1+ubound(zone3waterflux_glb_in)-lbound(zone3waterflux_glb_in)
      if(allocated(zone3waterflux_u_in))write(3,*)'zone3waterflux_u_in',1+ubound(zone3waterflux_u_in)-lbound(zone3waterflux_u_in)
      if(allocated(zone3waterflux_v_in))write(3,*)'zone3waterflux_v_in',1+ubound(zone3waterflux_v_in)-lbound(zone3waterflux_v_in)
      if(allocated(zone3saltflux_glb_out))write(3,*)'zone3saltflux_glb_out',1+ubound(zone3saltflux_glb_out)-lbound(zone3saltflux_glb_out)
      if(allocated(zone3saltflux_u_out))write(3,*)'zone3saltflux_u_out',1+ubound(zone3saltflux_u_out)-lbound(zone3saltflux_u_out)
      if(allocated(zone3saltflux_v_out))write(3,*)'zone3saltflux_v_out',1+ubound(zone3saltflux_v_out)-lbound(zone3saltflux_v_out)
      if(allocated(zone3tempflux_glb_out))write(3,*)'zone3tempflux_glb_out',1+ubound(zone3tempflux_glb_out)-lbound(zone3tempflux_glb_out)
      if(allocated(zone3tempflux_u_out))write(3,*)'zone3tempflux_u_out',1+ubound(zone3tempflux_u_out)-lbound(zone3tempflux_u_out)
      if(allocated(zone3tempflux_v_out))write(3,*)'zone3tempflux_v_out',1+ubound(zone3tempflux_v_out)-lbound(zone3tempflux_v_out)
      if(allocated(zone3waterflux_glb_out))write(3,*)'zone3waterflux_glb_out',1+ubound(zone3waterflux_glb_out)-lbound(zone3waterflux_glb_out)
      if(allocated(zone3waterflux_u_out))write(3,*)'zone3waterflux_u_out',1+ubound(zone3waterflux_u_out)-lbound(zone3waterflux_u_out)
      if(allocated(zone3waterflux_v_out))write(3,*)'zone3waterflux_v_out',1+ubound(zone3waterflux_v_out)-lbound(zone3waterflux_v_out)
      if(allocated(zone3_mask))write(3,*)'zone3_mask',1+ubound(zone3_mask)-lbound(zone3_mask)
      if(allocated(zone3_flux_u_node))write(3,*)'zone3_flux_u_node',1+ubound(zone3_flux_u_node)-lbound(zone3_flux_u_node)
      if(allocated(zone3_flux_v_node))write(3,*)'zone3_flux_v_node',1+ubound(zone3_flux_v_node)-lbound(zone3_flux_v_node)
      write(3,*)'zone4saltflux_w'
      write(3,*)'zone4tempflux_w'
      write(3,*)'zone4waterflux_w'
      write(3,*)'zone4_nlayer'
      write(3,*)'zone4_max'
      write(3,*)'zone4_u_max'
      write(3,*)'zone4_v_max'
      write(3,*)'zone4_inv_dz'
      write(3,*)'zone4_stretch_dz'
      if(allocated(zone4saltflux_glb))write(3,*)'zone4saltflux_glb',1+ubound(zone4saltflux_glb)-lbound(zone4saltflux_glb)
      if(allocated(zone4saltflux_u))write(3,*)'zone4saltflux_u',1+ubound(zone4saltflux_u)-lbound(zone4saltflux_u)
      if(allocated(zone4saltflux_v))write(3,*)'zone4saltflux_v',1+ubound(zone4saltflux_v)-lbound(zone4saltflux_v)
      if(allocated(zone4tempflux_glb))write(3,*)'zone4tempflux_glb',1+ubound(zone4tempflux_glb)-lbound(zone4tempflux_glb)
      if(allocated(zone4tempflux_u))write(3,*)'zone4tempflux_u',1+ubound(zone4tempflux_u)-lbound(zone4tempflux_u)
      if(allocated(zone4tempflux_v))write(3,*)'zone4tempflux_v',1+ubound(zone4tempflux_v)-lbound(zone4tempflux_v)
      if(allocated(zone4waterflux_glb))write(3,*)'zone4waterflux_glb',1+ubound(zone4waterflux_glb)-lbound(zone4waterflux_glb)
      if(allocated(zone4waterflux_u))write(3,*)'zone4waterflux_u',1+ubound(zone4waterflux_u)-lbound(zone4waterflux_u)
      if(allocated(zone4waterflux_v))write(3,*)'zone4waterflux_v',1+ubound(zone4waterflux_v)-lbound(zone4waterflux_v)
      if(allocated(zone4saltflux_glb_in))write(3,*)'zone4saltflux_glb_in',1+ubound(zone4saltflux_glb_in)-lbound(zone4saltflux_glb_in)
      if(allocated(zone4saltflux_u_in))write(3,*)'zone4saltflux_u_in',1+ubound(zone4saltflux_u_in)-lbound(zone4saltflux_u_in)
      if(allocated(zone4saltflux_v_in))write(3,*)'zone4saltflux_v_in',1+ubound(zone4saltflux_v_in)-lbound(zone4saltflux_v_in)
      if(allocated(zone4tempflux_glb_in))write(3,*)'zone4tempflux_glb_in',1+ubound(zone4tempflux_glb_in)-lbound(zone4tempflux_glb_in)
      if(allocated(zone4tempflux_u_in))write(3,*)'zone4tempflux_u_in',1+ubound(zone4tempflux_u_in)-lbound(zone4tempflux_u_in)
      if(allocated(zone4tempflux_v_in))write(3,*)'zone4tempflux_v_in',1+ubound(zone4tempflux_v_in)-lbound(zone4tempflux_v_in)
      if(allocated(zone4waterflux_glb_in))write(3,*)'zone4waterflux_glb_in',1+ubound(zone4waterflux_glb_in)-lbound(zone4waterflux_glb_in)
      if(allocated(zone4waterflux_u_in))write(3,*)'zone4waterflux_u_in',1+ubound(zone4waterflux_u_in)-lbound(zone4waterflux_u_in)
      if(allocated(zone4waterflux_v_in))write(3,*)'zone4waterflux_v_in',1+ubound(zone4waterflux_v_in)-lbound(zone4waterflux_v_in)
      if(allocated(zone4saltflux_glb_out))write(3,*)'zone4saltflux_glb_out',1+ubound(zone4saltflux_glb_out)-lbound(zone4saltflux_glb_out)
      if(allocated(zone4saltflux_u_out))write(3,*)'zone4saltflux_u_out',1+ubound(zone4saltflux_u_out)-lbound(zone4saltflux_u_out)
      if(allocated(zone4saltflux_v_out))write(3,*)'zone4saltflux_v_out',1+ubound(zone4saltflux_v_out)-lbound(zone4saltflux_v_out)
      if(allocated(zone4tempflux_glb_out))write(3,*)'zone4tempflux_glb_out',1+ubound(zone4tempflux_glb_out)-lbound(zone4tempflux_glb_out)
      if(allocated(zone4tempflux_u_out))write(3,*)'zone4tempflux_u_out',1+ubound(zone4tempflux_u_out)-lbound(zone4tempflux_u_out)
      if(allocated(zone4tempflux_v_out))write(3,*)'zone4tempflux_v_out',1+ubound(zone4tempflux_v_out)-lbound(zone4tempflux_v_out)
      if(allocated(zone4waterflux_glb_out))write(3,*)'zone4waterflux_glb_out',1+ubound(zone4waterflux_glb_out)-lbound(zone4waterflux_glb_out)
      if(allocated(zone4waterflux_u_out))write(3,*)'zone4waterflux_u_out',1+ubound(zone4waterflux_u_out)-lbound(zone4waterflux_u_out)
      if(allocated(zone4waterflux_v_out))write(3,*)'zone4waterflux_v_out',1+ubound(zone4waterflux_v_out)-lbound(zone4waterflux_v_out)
      if(allocated(zone4_mask))write(3,*)'zone4_mask',1+ubound(zone4_mask)-lbound(zone4_mask)
      if(allocated(zone4_flux_u_node))write(3,*)'zone4_flux_u_node',1+ubound(zone4_flux_u_node)-lbound(zone4_flux_u_node)
      if(allocated(zone4_flux_v_node))write(3,*)'zone4_flux_v_node',1+ubound(zone4_flux_v_node)-lbound(zone4_flux_v_node)
      write(3,*)'zone5saltflux_w'
      write(3,*)'zone5tempflux_w'
      write(3,*)'zone5waterflux_w'
      write(3,*)'zone5_nlayer'
      write(3,*)'zone5_max'
      write(3,*)'zone5_u_max'
      write(3,*)'zone5_v_max'
      write(3,*)'zone5_inv_dz'
      write(3,*)'zone5_stretch_dz'
      if(allocated(zone5saltflux_glb))write(3,*)'zone5saltflux_glb',1+ubound(zone5saltflux_glb)-lbound(zone5saltflux_glb)
      if(allocated(zone5saltflux_u))write(3,*)'zone5saltflux_u',1+ubound(zone5saltflux_u)-lbound(zone5saltflux_u)
      if(allocated(zone5saltflux_v))write(3,*)'zone5saltflux_v',1+ubound(zone5saltflux_v)-lbound(zone5saltflux_v)
      if(allocated(zone5tempflux_glb))write(3,*)'zone5tempflux_glb',1+ubound(zone5tempflux_glb)-lbound(zone5tempflux_glb)
      if(allocated(zone5tempflux_u))write(3,*)'zone5tempflux_u',1+ubound(zone5tempflux_u)-lbound(zone5tempflux_u)
      if(allocated(zone5tempflux_v))write(3,*)'zone5tempflux_v',1+ubound(zone5tempflux_v)-lbound(zone5tempflux_v)
      if(allocated(zone5waterflux_glb))write(3,*)'zone5waterflux_glb',1+ubound(zone5waterflux_glb)-lbound(zone5waterflux_glb)
      if(allocated(zone5waterflux_u))write(3,*)'zone5waterflux_u',1+ubound(zone5waterflux_u)-lbound(zone5waterflux_u)
      if(allocated(zone5waterflux_v))write(3,*)'zone5waterflux_v',1+ubound(zone5waterflux_v)-lbound(zone5waterflux_v)
      if(allocated(zone5saltflux_glb_in))write(3,*)'zone5saltflux_glb_in',1+ubound(zone5saltflux_glb_in)-lbound(zone5saltflux_glb_in)
      if(allocated(zone5saltflux_u_in))write(3,*)'zone5saltflux_u_in',1+ubound(zone5saltflux_u_in)-lbound(zone5saltflux_u_in)
      if(allocated(zone5saltflux_v_in))write(3,*)'zone5saltflux_v_in',1+ubound(zone5saltflux_v_in)-lbound(zone5saltflux_v_in)
      if(allocated(zone5tempflux_glb_in))write(3,*)'zone5tempflux_glb_in',1+ubound(zone5tempflux_glb_in)-lbound(zone5tempflux_glb_in)
      if(allocated(zone5tempflux_u_in))write(3,*)'zone5tempflux_u_in',1+ubound(zone5tempflux_u_in)-lbound(zone5tempflux_u_in)
      if(allocated(zone5tempflux_v_in))write(3,*)'zone5tempflux_v_in',1+ubound(zone5tempflux_v_in)-lbound(zone5tempflux_v_in)
      if(allocated(zone5waterflux_glb_in))write(3,*)'zone5waterflux_glb_in',1+ubound(zone5waterflux_glb_in)-lbound(zone5waterflux_glb_in)
      if(allocated(zone5waterflux_u_in))write(3,*)'zone5waterflux_u_in',1+ubound(zone5waterflux_u_in)-lbound(zone5waterflux_u_in)
      if(allocated(zone5waterflux_v_in))write(3,*)'zone5waterflux_v_in',1+ubound(zone5waterflux_v_in)-lbound(zone5waterflux_v_in)
      if(allocated(zone5saltflux_glb_out))write(3,*)'zone5saltflux_glb_out',1+ubound(zone5saltflux_glb_out)-lbound(zone5saltflux_glb_out)
      if(allocated(zone5saltflux_u_out))write(3,*)'zone5saltflux_u_out',1+ubound(zone5saltflux_u_out)-lbound(zone5saltflux_u_out)
      if(allocated(zone5saltflux_v_out))write(3,*)'zone5saltflux_v_out',1+ubound(zone5saltflux_v_out)-lbound(zone5saltflux_v_out)
      if(allocated(zone5tempflux_glb_out))write(3,*)'zone5tempflux_glb_out',1+ubound(zone5tempflux_glb_out)-lbound(zone5tempflux_glb_out)
      if(allocated(zone5tempflux_u_out))write(3,*)'zone5tempflux_u_out',1+ubound(zone5tempflux_u_out)-lbound(zone5tempflux_u_out)
      if(allocated(zone5tempflux_v_out))write(3,*)'zone5tempflux_v_out',1+ubound(zone5tempflux_v_out)-lbound(zone5tempflux_v_out)
      if(allocated(zone5waterflux_glb_out))write(3,*)'zone5waterflux_glb_out',1+ubound(zone5waterflux_glb_out)-lbound(zone5waterflux_glb_out)
      if(allocated(zone5waterflux_u_out))write(3,*)'zone5waterflux_u_out',1+ubound(zone5waterflux_u_out)-lbound(zone5waterflux_u_out)
      if(allocated(zone5waterflux_v_out))write(3,*)'zone5waterflux_v_out',1+ubound(zone5waterflux_v_out)-lbound(zone5waterflux_v_out)
      if(allocated(zone5_mask))write(3,*)'zone5_mask',1+ubound(zone5_mask)-lbound(zone5_mask)
      if(allocated(zone5_flux_u_node))write(3,*)'zone5_flux_u_node',1+ubound(zone5_flux_u_node)-lbound(zone5_flux_u_node)
      if(allocated(zone5_flux_v_node))write(3,*)'zone5_flux_v_node',1+ubound(zone5_flux_v_node)-lbound(zone5_flux_v_node)
      write(3,*)'zone6saltflux_w'
      write(3,*)'zone6tempflux_w'
      write(3,*)'zone6waterflux_w'
      write(3,*)'zone6_nlayer'
      write(3,*)'zone6_max'
      write(3,*)'zone6_u_max'
      write(3,*)'zone6_v_max'
      write(3,*)'zone6_inv_dz'
      write(3,*)'zone6_stretch_dz'
      if(allocated(zone6saltflux_glb))write(3,*)'zone6saltflux_glb',1+ubound(zone6saltflux_glb)-lbound(zone6saltflux_glb)
      if(allocated(zone6saltflux_u))write(3,*)'zone6saltflux_u',1+ubound(zone6saltflux_u)-lbound(zone6saltflux_u)
      if(allocated(zone6saltflux_v))write(3,*)'zone6saltflux_v',1+ubound(zone6saltflux_v)-lbound(zone6saltflux_v)
      if(allocated(zone6tempflux_glb))write(3,*)'zone6tempflux_glb',1+ubound(zone6tempflux_glb)-lbound(zone6tempflux_glb)
      if(allocated(zone6tempflux_u))write(3,*)'zone6tempflux_u',1+ubound(zone6tempflux_u)-lbound(zone6tempflux_u)
      if(allocated(zone6tempflux_v))write(3,*)'zone6tempflux_v',1+ubound(zone6tempflux_v)-lbound(zone6tempflux_v)
      if(allocated(zone6waterflux_glb))write(3,*)'zone6waterflux_glb',1+ubound(zone6waterflux_glb)-lbound(zone6waterflux_glb)
      if(allocated(zone6waterflux_u))write(3,*)'zone6waterflux_u',1+ubound(zone6waterflux_u)-lbound(zone6waterflux_u)
      if(allocated(zone6waterflux_v))write(3,*)'zone6waterflux_v',1+ubound(zone6waterflux_v)-lbound(zone6waterflux_v)
      if(allocated(zone6saltflux_glb_in))write(3,*)'zone6saltflux_glb_in',1+ubound(zone6saltflux_glb_in)-lbound(zone6saltflux_glb_in)
      if(allocated(zone6saltflux_u_in))write(3,*)'zone6saltflux_u_in',1+ubound(zone6saltflux_u_in)-lbound(zone6saltflux_u_in)
      if(allocated(zone6saltflux_v_in))write(3,*)'zone6saltflux_v_in',1+ubound(zone6saltflux_v_in)-lbound(zone6saltflux_v_in)
      if(allocated(zone6tempflux_glb_in))write(3,*)'zone6tempflux_glb_in',1+ubound(zone6tempflux_glb_in)-lbound(zone6tempflux_glb_in)
      if(allocated(zone6tempflux_u_in))write(3,*)'zone6tempflux_u_in',1+ubound(zone6tempflux_u_in)-lbound(zone6tempflux_u_in)
      if(allocated(zone6tempflux_v_in))write(3,*)'zone6tempflux_v_in',1+ubound(zone6tempflux_v_in)-lbound(zone6tempflux_v_in)
      if(allocated(zone6waterflux_glb_in))write(3,*)'zone6waterflux_glb_in',1+ubound(zone6waterflux_glb_in)-lbound(zone6waterflux_glb_in)
      if(allocated(zone6waterflux_u_in))write(3,*)'zone6waterflux_u_in',1+ubound(zone6waterflux_u_in)-lbound(zone6waterflux_u_in)
      if(allocated(zone6waterflux_v_in))write(3,*)'zone6waterflux_v_in',1+ubound(zone6waterflux_v_in)-lbound(zone6waterflux_v_in)
      if(allocated(zone6saltflux_glb_out))write(3,*)'zone6saltflux_glb_out',1+ubound(zone6saltflux_glb_out)-lbound(zone6saltflux_glb_out)
      if(allocated(zone6saltflux_u_out))write(3,*)'zone6saltflux_u_out',1+ubound(zone6saltflux_u_out)-lbound(zone6saltflux_u_out)
      if(allocated(zone6saltflux_v_out))write(3,*)'zone6saltflux_v_out',1+ubound(zone6saltflux_v_out)-lbound(zone6saltflux_v_out)
      if(allocated(zone6tempflux_glb_out))write(3,*)'zone6tempflux_glb_out',1+ubound(zone6tempflux_glb_out)-lbound(zone6tempflux_glb_out)
      if(allocated(zone6tempflux_u_out))write(3,*)'zone6tempflux_u_out',1+ubound(zone6tempflux_u_out)-lbound(zone6tempflux_u_out)
      if(allocated(zone6tempflux_v_out))write(3,*)'zone6tempflux_v_out',1+ubound(zone6tempflux_v_out)-lbound(zone6tempflux_v_out)
      if(allocated(zone6waterflux_glb_out))write(3,*)'zone6waterflux_glb_out',1+ubound(zone6waterflux_glb_out)-lbound(zone6waterflux_glb_out)
      if(allocated(zone6waterflux_u_out))write(3,*)'zone6waterflux_u_out',1+ubound(zone6waterflux_u_out)-lbound(zone6waterflux_u_out)
      if(allocated(zone6waterflux_v_out))write(3,*)'zone6waterflux_v_out',1+ubound(zone6waterflux_v_out)-lbound(zone6waterflux_v_out)
      if(allocated(zone6_mask))write(3,*)'zone6_mask',1+ubound(zone6_mask)-lbound(zone6_mask)
      if(allocated(zone6_flux_u_node))write(3,*)'zone6_flux_u_node',1+ubound(zone6_flux_u_node)-lbound(zone6_flux_u_node)
      if(allocated(zone6_flux_v_node))write(3,*)'zone6_flux_v_node',1+ubound(zone6_flux_v_node)-lbound(zone6_flux_v_node)
      if(allocated(anyvar3d))write(3,*)'anyvar3d',1+ubound(anyvar3d)-lbound(anyvar3d)
      if(allocated(sigma_fric_wu))write(3,*)'sigma_fric_wu',1+ubound(sigma_fric_wu)-lbound(sigma_fric_wu)
      if(allocated(sigma_fric_wv))write(3,*)'sigma_fric_wv',1+ubound(sigma_fric_wv)-lbound(sigma_fric_wv)
      if(allocated(dsig_t))write(3,*)'dsig_t',1+ubound(dsig_t)-lbound(dsig_t)
      if(allocated(anyv3dint))write(3,*)'anyv3dint',1+ubound(anyv3dint)-lbound(anyv3dint)
      if(allocated(sshobc_w))write(3,*)'sshobc_w',1+ubound(sshobc_w)-lbound(sshobc_w)
      if(allocated(velbarobc_u))write(3,*)'velbarobc_u',1+ubound(velbarobc_u)-lbound(velbarobc_u)
      if(allocated(velbarobc_v))write(3,*)'velbarobc_v',1+ubound(velbarobc_v)-lbound(velbarobc_v)
      if(ioffline_prv>=1.and.allocated(sshofl_w))write(3,*)'sshofl_w',1+ubound(sshofl_w)-lbound(sshofl_w)
      if(ioffline_prv>=1.and.allocated(velbarofl_u))write(3,*)'velbarofl_u',1+ubound(velbarofl_u)-lbound(velbarofl_u)
      if(ioffline_prv>=1.and.allocated(velbarofl_v))write(3,*)'velbarofl_v',1+ubound(velbarofl_v)-lbound(velbarofl_v)
      if(allocated(tem_delta_t))write(3,*)'tem_delta_t',1+ubound(tem_delta_t)-lbound(tem_delta_t)
      if(allocated(sal_delta_t))write(3,*)'sal_delta_t',1+ubound(sal_delta_t)-lbound(sal_delta_t)
      if(allocated(anyvar2d))write(3,*)'anyvar2d',1+ubound(anyvar2d)-lbound(anyvar2d)
      if(allocated(iriver))write(3,*)'iriver',1+ubound(iriver)-lbound(iriver)
      if(allocated(jriver))write(3,*)'jriver',1+ubound(jriver)-lbound(jriver)
      if(allocated(rankcoords))write(3,*)'rankcoords',1+ubound(rankcoords)-lbound(rankcoords)
      if(allocated(l2ij_out_u))write(3,*)'l2ij_out_u',1+ubound(l2ij_out_u)-lbound(l2ij_out_u)
      if(allocated(l2ij_out_v))write(3,*)'l2ij_out_v',1+ubound(l2ij_out_v)-lbound(l2ij_out_v)
      if(allocated(l2ij_out_w))write(3,*)'l2ij_out_w',1+ubound(l2ij_out_w)-lbound(l2ij_out_w)
      if(allocated(l2ij_in_u))write(3,*)'l2ij_in_u',1+ubound(l2ij_in_u)-lbound(l2ij_in_u)
      if(allocated(l2ij_in_v))write(3,*)'l2ij_in_v',1+ubound(l2ij_in_v)-lbound(l2ij_in_v)
      if(allocated(l2ij_in_w))write(3,*)'l2ij_in_w',1+ubound(l2ij_in_w)-lbound(l2ij_in_w)
      if(allocated(fillmask_t))write(3,*)'fillmask_t',1+ubound(fillmask_t)-lbound(fillmask_t)
      if(allocated(datesim))write(3,*)'datesim',1+ubound(datesim)-lbound(datesim)
      if(allocated(dateobc))write(3,*)'dateobc',1+ubound(dateobc)-lbound(dateobc)
      if(allocated(dateairsea))write(3,*)'dateairsea',1+ubound(dateairsea)-lbound(dateairsea)
      if(allocated(dateriver))write(3,*)'dateriver',1+ubound(dateriver)-lbound(dateriver)
      if(allocated(runoff_w))write(3,*)'runoff_w',1+ubound(runoff_w)-lbound(runoff_w)
      write(3,*)'drifter_out_sampling'
      write(3,*)'obc2dtype'
      write(3,*)'coef_diss_mangrove'
      write(3,*)'expnum'
      write(3,*)'cst_c0cub'
      write(3,*)'relativewind'
      write(3,*)'offset_sshobc'
      write(3,*)'biharm_2dfactor'
      write(3,*)'checkr0'
      write(3,*)'checkr1'
      write(3,*)'checkr2'
      write(3,*)'checkr3'
      write(3,*)'checkr4'
      write(3,*)'sponge_l'
      write(3,*)'sponge_dx_critic'
      write(3,*)'sponge_dx_width'
      write(3,*)'dzsurfmin'
      if(allocated(drifter_send_order_canal))write(3,*)'drifter_send_order_canal',1+ubound(drifter_send_order_canal)-lbound(drifter_send_order_canal)
      write(3,*)'drifter_onoff'
      if(allocated(alongresolriver))write(3,*)'alongresolriver',1+ubound(alongresolriver)-lbound(alongresolriver)
      if(allocated(crossresolriver))write(3,*)'crossresolriver',1+ubound(crossresolriver)-lbound(crossresolriver)
      write(3,*)'ww3_varmax'
      write(3,*)'ww3_type_grid'
      write(3,*)'type_unstructured'
      write(3,*)'type_structured'
      write(3,*)'vststep'
      write(3,*)'zprofile2d3d'
      write(3,*)'timestep_type'
      write(3,*)'timestep_leapfrog'
      write(3,*)'timestep_forwbckw'
      write(3,*)'bulk_core'
      write(3,*)'bulk_moon'
      write(3,*)'bulk_coare'
      write(3,*)'bulk_ecume'
      write(3,*)'bulk_scheme'
      write(3,*)'tideana_spinup'
      write(3,*)'tideana_delta'
      write(3,*)'ogcm_time_lag'
      write(3,*)'albedo_val'
      write(3,*)'tfilterfb'
      write(3,*)'tfilterlf'
      write(3,*)'ktide'
      write(3,*)'tideforces'
      write(3,*)'tideana_yesno'
      write(3,*)'ioffline'
      write(3,*)'tideanalysis_count'
      write(3,*)'ibl1_advbio'
      write(3,*)'ibl2_advbio'
      write(3,*)'jbl1_advbio'
      write(3,*)'jbl2_advbio'
      write(3,*)'ogcm_time_shift'
      write(3,*)'ieq1'
      write(3,*)'ieqimax'
      write(3,*)'jeq1'
      write(3,*)'jeqjmax'
      write(3,*)'ieq1_jeq1'
      write(3,*)'ieq1_jeqjmax'
      write(3,*)'ieqimax_jeq1'
      write(3,*)'ieqimax_jeqjmax'
      write(3,*)'dim_varid'
      write(3,*)'albedo_constant'
      write(3,*)'albedo_apel1987'
      write(3,*)'albedo_br1982'
      write(3,*)'albedo_br1986'
      write(3,*)'albedo_case'
      write(3,*)'il_oasis_time'
      write(3,*)'flag_steric_effect'
      write(3,*)'flag_remove_secondary_bassins'
      write(3,*)'one_kind1'
      write(3,*)'signe'
      write(3,*)'flag_maxbotstress'
      write(3,*)'flag_offline_binary'
      write(3,*)'flag_wstressbulk'
      write(3,*)'flag_write'
      write(3,*)'mangrove_scheme'
      write(3,*)'flag_meteo_land_plug'
      write(3,*)'flag_meteo_land_plug_wind'
      write(3,*)'flag_upwind_obc'
      write(3,*)'discard_lonlat_periodicity'
      write(3,*)'fplan2_grid'
      write(3,*)'flag_ts_effectivedensity'
      write(3,*)'flag_ts_quicklim'
      write(3,*)'flag_z2dv_outputs'
      write(3,*)'flag_ogcmname2date'
      write(3,*)'flag_meteoname2date'
      write(3,*)'ofl_bio'
      write(3,*)'drifter_output_files'
      write(3,*)'flag_tide3d_analysis'
      write(3,*)'flag_bathy_update'
      write(3,*)'loop1'
      write(3,*)'loop2'
      write(3,*)'loop3'
      write(3,*)'loopmaxtke'
      write(3,*)'loopmaxbio'
      write(3,*)'loopmaxts'
      write(3,*)'nairsea'
      write(3,*)'bi_onoff'
      write(3,*)'nc_or_bin_airsea'
      write(3,*)'loop_netcdf'
      write(3,*)'count_netcdfvar'
      write(3,*)'wavefile_prvtrec'
      write(3,*)'wavefile_nextrec'
      write(3,*)'wave_cpl_nextrec'
      write(3,*)'wave_cpl_period_iter'
      write(3,*)'wave_cpl_ww3_sdir'
      write(3,*)'wave_cpl_ww3_cdir'
      write(3,*)'wave_cpl_ww3_hs'
      write(3,*)'wave_cpl_ww3_hsw'
      write(3,*)'wave_cpl_ww3_foc'
      write(3,*)'wave_cpl_ww3_tw'
      write(3,*)'wave_cpl_ww3_tawx'
      write(3,*)'wave_cpl_ww3_tawy'
      write(3,*)'wave_cpl_ww3_twox'
      write(3,*)'wave_cpl_ww3_twoy'
      write(3,*)'wave_cpl_ww3_uss'
      write(3,*)'wave_cpl_ww3_vss'
      write(3,*)'wave_cpl_ww3_msk'
      write(3,*)'ofl_rec_now'
      write(3,*)'ofl_rec_max'
      write(3,*)'flag_kz_enhanced'
      write(3,*)'flag_net_ir'
      write(3,*)'flag_dt_adjust'
      write(3,*)'tide_interpolation'
      write(3,*)'tide_flagrotation'
      write(3,*)'flag_ssr24avr'
      write(3,*)'flag_abl'
      write(3,*)'flag_abl2'
      write(3,*)'flag_sequoia'
      write(3,*)'dimssr24prv'
      write(3,*)'flag_meteo_average'
      write(3,*)'flag_nemoffline'
      write(3,*)'trc_id'
      write(3,*)'vel_id'
      write(3,*)'ssh_id'
      write(3,*)'var_num'
      write(3,*)'flag_p0m_filter'
      write(3,*)'flag_refstate'
      write(3,*)'flag_linearfric'
      write(3,*)'linear_coef_mangrove'
      write(3,*)'flag_1dv'
      write(3,*)'flag_merged_levels'
      write(3,*)'flag1_smooth_h_mask'
      write(3,*)'flag3_smooth_h_mask'
      write(3,*)'flag_0status_option'
      write(3,*)'flag_rmnegval'
      write(3,*)'timemax'
      write(3,*)'dirmax'
      write(3,*)'freqmax'
      write(3,*)'var_misval'
      write(3,*)'un_r8'
      write(3,*)'heure'
      write(3,*)'suma'
      write(3,*)'sumb'
      write(3,*)'sum_mi'
      write(3,*)'grid_area'
      write(3,*)'grid_areaglb'
      write(3,*)'grid_volumeglb'
      write(3,*)'sumarchive'
      write(3,*)'lonmin'
      write(3,*)'lonmax'
      write(3,*)'latmin'
      write(3,*)'latmax'
      write(3,*)'wavefile_prvtime'
      write(3,*)'wavefile_nextime'
      write(3,*)'ofl_period_prev'
      write(3,*)'ofl_period_now'
      write(3,*)'ofl_period_next'
      write(3,*)'ofl_nextrec_time'
      write(3,*)'ofl_writime'
      write(3,*)'ofl_readtime_next'
      write(3,*)'ofl_readtime_prev'
      write(3,*)'biobc_nextfiletime'
      write(3,*)'biobc_prevfiletime'
      write(3,*)'stability_index'
      write(3,*)'iteration2d_max_r8'
      write(3,*)'iteration2d_upbound'
      write(3,*)'coef_linearfric'
      write(3,*)'gh'
      write(3,*)'gm'
      write(3,*)'nn'
      write(3,*)'tke2overeps'
      write(3,*)'tkeovereps'
      write(3,*)'gravoverrho'
      write(3,*)'sh'
      write(3,*)'sm'
      write(3,*)'heatfluxbias'
      write(3,*)'check0'
      write(3,*)'check1'
      write(3,*)'check2'
      write(3,*)'check3'
      write(3,*)'check4'
      write(3,*)'check5'
      write(3,*)'check6'
      write(3,*)'check7'
      write(3,*)'check8'
      write(3,*)'zero'
      write(3,*)'un'
      write(3,*)'cdseuil'
      write(3,*)'cdb_2dh'
      write(3,*)'small1'
      write(3,*)'deci'
      write(3,*)'decj'
      write(3,*)'deck'
      write(3,*)'rap1'
      write(3,*)'rap2'
      write(3,*)'rapi'
      write(3,*)'rapj'
      write(3,*)'rapk'
      write(3,*)'rap'
      write(3,*)'const0'
      write(3,*)'const1'
      write(3,*)'const2'
      write(3,*)'const3'
      write(3,*)'const4'
      write(3,*)'const5'
      write(3,*)'const6'
      write(3,*)'const7'
      write(3,*)'const8'
      write(3,*)'const9'
      write(3,*)'uv_10'
      write(3,*)'karman'
      write(3,*)'stefan'
      write(3,*)'z2m'
      write(3,*)'z10m'
      write(3,*)'pss0'
      write(3,*)'boltz'
      write(3,*)'planck'
      write(3,*)'avogadro'
      write(3,*)'ce'
      write(3,*)'cen'
      write(3,*)'ch'
      write(3,*)'chn'
      write(3,*)'cd'
      write(3,*)'cdn'
      write(3,*)'z_2'
      write(3,*)'z_10'
      write(3,*)'q_0'
      write(3,*)'r_0'
      write(3,*)'pvs_0'
      write(3,*)'psih_10'
      write(3,*)'psih_2'
      write(3,*)'phim'
      write(3,*)'zl_10'
      write(3,*)'zl_2'
      write(3,*)'falpha'
      write(3,*)'fbeta'
      write(3,*)'ro'
      write(3,*)'cp_air'
      write(3,*)'lv'
      write(3,*)'psim_10'
      write(3,*)'psim_2'
      write(3,*)'dte_lp'
      write(3,*)'inv_dte_lp'
      write(3,*)'inv_dti_lp'
      write(3,*)'inv_dti_fw'
      write(3,*)'inv_dti_fw_p2'
      write(3,*)'dte_fw'
      write(3,*)'dti_lpbef'
      write(3,*)'dti_lp'
      write(3,*)'dti_lpmax'
      write(3,*)'dti_lpsub'
      write(3,*)'dti_fwsub'
      write(3,*)'dti_fwsubio'
      write(3,*)'dti_fw'
      write(3,*)'dti_bef'
      write(3,*)'dti_now'
      write(3,*)'dtiratio'
      write(3,*)'dt_drf'
      write(3,*)'fbtfiltercoef'
      write(3,*)'assel0'
      write(3,*)'assel1'
      write(3,*)'assel2'
      write(3,*)'assel3'
      write(3,*)'wetdry_cst1'
      write(3,*)'wetdry_cst2'
      write(3,*)'wetdry_cst3'
      write(3,*)'h_inf'
      write(3,*)'h_inf_obc'
      write(3,*)'h_sup'
      write(3,*)'dist'
      write(3,*)'dist0'
      write(3,*)'dist1'
      write(3,*)'dist2'
      write(3,*)'dist3'
      write(3,*)'t0surf'
      write(3,*)'s0surf'
      write(3,*)'deg2rad'
      write(3,*)'rad2deg'
      write(3,*)'zmin'
      write(3,*)'zmax'
      write(3,*)'coa'
      write(3,*)'cob'
      write(3,*)'coc'
      write(3,*)'sol1'
      write(3,*)'sol2'
      write(3,*)'discri'
      write(3,*)'profm1'
      write(3,*)'profp1'
      write(3,*)'hmax'
      write(3,*)'hstepmin'
      write(3,*)'hstepmax'
      write(3,*)'grav'
      write(3,*)'cfl_sshmax'
      write(3,*)'cfl_hsshmax'
      write(3,*)'cfl_umax'
      write(3,*)'cfl_reduce'
      write(3,*)'relax_es'
      write(3,*)'relax_ext'
      write(3,*)'relax_int'
      write(3,*)'relax_ts'
      write(3,*)'relax_bpc'
      write(3,*)'momentum_input_depth'
      write(3,*)'z0s'
      write(3,*)'z1'
      write(3,*)'z2'
      write(3,*)'z3'
      write(3,*)'z4'
      write(3,*)'z5'
      write(3,*)'y0'
      write(3,*)'y2'
      write(3,*)'y3'
      write(3,*)'y4'
      write(3,*)'y5'
      write(3,*)'x0'
      write(3,*)'x1'
      write(3,*)'x2'
      write(3,*)'x3'
      write(3,*)'x4'
      write(3,*)'x5'
      write(3,*)'x6'
      write(3,*)'x7'
      write(3,*)'x8'
      write(3,*)'x9'
      write(3,*)'x10'
      write(3,*)'x11'
      write(3,*)'x12'
      write(3,*)'x13'
      write(3,*)'x14'
      write(3,*)'x20'
      write(3,*)'x21'
      write(3,*)'x22'
      write(3,*)'x33'
      write(3,*)'x44'
      write(3,*)'area'
      write(3,*)'tem_validmin'
      write(3,*)'tem_validmax'
      write(3,*)'sal_validmin'
      write(3,*)'sal_validmax'
      write(3,*)'xmd'
      write(3,*)'xmv'
      write(3,*)'xrd'
      write(3,*)'xrv'
      write(3,*)'xcpd'
      write(3,*)'xcpv'
      write(3,*)'xcl'
      write(3,*)'celsius2kelvin'
      write(3,*)'xlvtt'
      write(3,*)'xestt'
      write(3,*)'xgamw'
      write(3,*)'xbetaw'
      write(3,*)'xalpw'
      write(3,*)'zrvsrdm1'
      write(3,*)'z0_u'
      write(3,*)'qsat_sea_z'
      write(3,*)'airdensity'
      write(3,*)'sst_kelvin'
      write(3,*)'prs_atm_z'
      write(3,*)'tem_atm_z'
      write(3,*)'exner_atm_z'
      write(3,*)'delta_u'
      write(3,*)'delta_t'
      write(3,*)'delta_q'
      write(3,*)'delta_u_n'
      write(3,*)'exner_sea_z'
      write(3,*)'psifunctt'
      write(3,*)'psifunctu'
      write(3,*)'z0_q'
      write(3,*)'z0_t'
      write(3,*)'sst1000hpa_kelvin'
      write(3,*)'visa'
      write(3,*)'charnock'
      write(3,*)'qsat_atm_z'
      write(3,*)'ustar_bef'
      write(3,*)'qstar_bef'
      write(3,*)'tetastar_bef'
      write(3,*)'rayonterre'
      write(3,*)'northpole_lon'
      write(3,*)'northpole_lat'
      write(3,*)'southpole_lon'
      write(3,*)'southpole_lat'
      write(3,*)'phi0'
      write(3,*)'longi'
      write(3,*)'latit'
      write(3,*)'longi0'
      write(3,*)'latit0'
      write(3,*)'longi1'
      write(3,*)'latit1'
      write(3,*)'angle0'
      write(3,*)'alp_t'
      write(3,*)'alp_s'
      write(3,*)'pi'
      write(3,*)'rho'
      write(3,*)'inv_rho'
      write(3,*)'rhoair'
      write(3,*)'valmax'
      write(3,*)'vis'
      write(3,*)'tkee1'
      write(3,*)'tkee2'
      write(3,*)'tkee3'
      write(3,*)'tkeg2'
      write(3,*)'tkeg3'
      write(3,*)'tkeg4'
      write(3,*)'tkeg5'
      write(3,*)'tkeg6'
      write(3,*)'tkeb1'
      write(3,*)'tkeb2'
      write(3,*)'ctke1'
      write(3,*)'ctke2'
      write(3,*)'t0'
      write(3,*)'s0'
      write(3,*)'cp'
      write(3,*)'light_kpar1'
      write(3,*)'light_att1'
      write(3,*)'light_att2'
      write(3,*)'light_rat1'
      write(3,*)'light_rat2'
      write(3,*)'light_att2_val1'
      write(3,*)'light_att2_h1'
      write(3,*)'light_att2_val2'
      write(3,*)'light_att2_h2'
      write(3,*)'t0_base'
      write(3,*)'s0_base'
      write(3,*)'rho_base'
      write(3,*)'alp_t_base'
      write(3,*)'alp_s_base'
      write(3,*)'tfb0'
      write(3,*)'meteo_lonmin'
      write(3,*)'meteo_latmin'
      write(3,*)'meteo_lonmax'
      write(3,*)'meteo_latmax'
      write(3,*)'meteo_resol'
      write(3,*)'meteo_resol_u'
      write(3,*)'meteo_resol_v'
      write(3,*)'meteo_lonstr'
      write(3,*)'meteo_lonend'
      write(3,*)'meteo_londlt'
      write(3,*)'meteo_latstr'
      write(3,*)'meteo_latend'
      write(3,*)'meteo_latdlt'
      write(3,*)'var_lonmin'
      write(3,*)'var_latmin'
      write(3,*)'var_lonmax'
      write(3,*)'var_latmax'
      write(3,*)'ww3_lonmin'
      write(3,*)'ww3_latmin'
      write(3,*)'ww3_lonmax'
      write(3,*)'ww3_latmax'
      write(3,*)'ww3_dlon'
      write(3,*)'ww3_dlat'
      write(3,*)'tide_lonmin'
      write(3,*)'tide_latmin'
      write(3,*)'tide_lonmax'
      write(3,*)'tide_latmax'
      write(3,*)'tide_dlon'
      write(3,*)'tide_dlat'
      write(3,*)'dxb'
      write(3,*)'dyb'
      write(3,*)'dxa'
      write(3,*)'dya'
      write(3,*)'hmin'
      write(3,*)'h1d'
      write(3,*)'dlon'
      write(3,*)'dlat'
      write(3,*)'epsi'
      write(3,*)'lagrange_ssh'
      write(3,*)'diffu'
      write(3,*)'difnorm'
      write(3,*)'lup'
      write(3,*)'ldown'
      write(3,*)'zup'
      write(3,*)'zdown'
      write(3,*)'rbase'
      write(3,*)'small'
      write(3,*)'small3'
      write(3,*)'hgesig'
      write(3,*)'pgesig'
      write(3,*)'windfactor'
      write(3,*)'xdtk_out'
      write(3,*)'tfond'
      write(3,*)'sfond'
      write(3,*)'rfond'
      write(3,*)'c1streamf'
      write(3,*)'c2streamf'
      write(3,*)'rampe'
      write(3,*)'rampe_wind'
      write(3,*)'y1'
      write(3,*)'obc_hf_reset'
      write(3,*)'rho_0d'
      write(3,*)'tem_0d'
      write(3,*)'sal_0d'
      write(3,*)'rho_tmp'
      write(3,*)'cst_adv_hor'
      write(3,*)'cst_adv_ver'
      write(3,*)'cst_adv_vel'
      write(3,*)'ssh_avr_nest_out'
      write(3,*)'graph_nextime'
      write(3,*)'tidenodal_prev_rdv'
      write(3,*)'tidenodal_next_rdv'
      write(3,*)'tideana_modulo'
      write(3,*)'tideana_nextime'
      write(3,*)'cellboxfactor1'
      write(3,*)'cellboxfactor2'
      write(3,*)'rap_wave'
      write(3,*)'ihmax'
      write(3,*)'jhmax'
      write(3,*)'ihmin'
      write(3,*)'jhmin'
      write(3,*)'convect_yn'
      write(3,*)'nbvstepmin'
      write(3,*)'bio_relax_size'
      write(3,*)'unit_r4'
      write(3,*)'x0_r4'
      write(3,*)'x1_r4'
      write(3,*)'x2_r4'
      write(3,*)'x3_r4'
      write(3,*)'x4_r4'
      write(3,*)'time_r4'
      write(3,*)'discharge'
      write(3,*)'filval'
      write(3,*)'var_validmin'
      write(3,*)'var_validmax'
      write(3,*)'vdw_loc'
      write(3,*)'vup_loc'
      write(3,*)'zdw_loc'
      write(3,*)'var_scalefactor'
      write(3,*)'inv_scalefactor'
      write(3,*)'var_addoffset'
      write(3,*)'zup_loc'
      write(3,*)'hrmax'
      write(3,*)'relax_bio'
      write(3,*)'kmol_m'
      write(3,*)'kmol_h'
      write(3,*)'kmol_s'
      write(3,*)'upw_hrange1'
      write(3,*)'upw_hrange2'
      write(3,*)'inv_ekman_depth'
      write(3,*)'constant_km'
      write(3,*)'constant_kh'
      write(3,*)'z0b'
      write(3,*)'z0b_land'
      write(3,*)'z0b_rivers'
      write(3,*)'zlevel_land'
      write(3,*)'invloopmaxts'
      write(3,*)'invloopmaxu'
      write(3,*)'invloopmaxv'
      write(3,*)'sqrtgrav'
      write(3,*)'invgrav'
      write(3,*)'spinup_forcing'
      write(3,*)'relax_lwf'
      write(3,*)'dz_vertical_incr_fact'
      write(3,*)'ratio_negdif_ver'
      write(3,*)'ratio_negdif_hor'
      write(3,*)'ratio_bionegdif'
      write(3,*)'vqs_cst1'
      write(3,*)'vqs_cst2'
      write(3,*)'vqs_cst3'
      write(3,*)'ema_mu'
      write(3,*)'dz_over_z0_min'
      write(3,*)'coastal_viscosity'
      write(3,*)'quick_coef'
      write(3,*)'i'
      write(3,*)'j'
      write(3,*)'k'
      write(3,*)'flag3d'
      write(3,*)'lrec'
      write(3,*)'iteration3d'
      write(3,*)'compt1'
      write(3,*)'compt2'
      write(3,*)'compt3'
      write(3,*)'compt4'
      write(3,*)'kount0'
      write(3,*)'kount1'
      write(3,*)'kount2'
      write(3,*)'kount3'
      write(3,*)'kount4'
      write(3,*)'kount5'
      write(3,*)'kount6'
      write(3,*)'kount7'
      write(3,*)'kount8'
      write(3,*)'kount9'
      write(3,*)'kountmod'
      write(3,*)'kountrdv1'
      write(3,*)'kountrdv2'
      write(3,*)'kountrdv3'
      write(3,*)'kountrdv4'
      write(3,*)'substep_advbio'
      write(3,*)'subcycle_exchange'
      write(3,*)'subcycle_onoff'
      write(3,*)'subcycle_synchro'
      write(3,*)'subcycle_modulo'
      write(3,*)'dt_drf_over_dti_fw'
      write(3,*)'iteration3d_restart'
      write(3,*)'filvalshort'
      write(3,*)'quick_filter_points'
      write(3,*)'status'
      write(3,*)'forcedstatus'
      write(3,*)'decision'
      write(3,*)'ncid1'
      write(3,*)'ncid2'
      write(3,*)'dim_x_id'
      write(3,*)'dim_y_id'
      write(3,*)'dim_z_id'
      write(3,*)'dim_t_id'
      write(3,*)'dim_b_id'
      write(3,*)'max_x'
      write(3,*)'max_y'
      write(3,*)'max_z'
      write(3,*)'max_time_counter'
      write(3,*)'max_meteo_time_counter'
      write(3,*)'var_id'
      write(3,*)'var_nftype'
      write(3,*)'meteo_imax'
      write(3,*)'meteo_jmax'
      write(3,*)'meteo_kmax'
      write(3,*)'meteozoom_istr'
      write(3,*)'meteozoom_iend'
      write(3,*)'meteozoom_jstr'
      write(3,*)'meteozoom_jend'
      write(3,*)'meteofull_imax'
      write(3,*)'meteofull_jmax'
      write(3,*)'tide_imax'
      write(3,*)'tide_jmax'
      write(3,*)'tide_kmax'
      write(3,*)'tidezoom_istr'
      write(3,*)'tidezoom_iend'
      write(3,*)'tidezoom_jstr'
      write(3,*)'tidezoom_jend'
      write(3,*)'tidezoom_istr_t'
      write(3,*)'tidezoom_iend_t'
      write(3,*)'tidezoom_jstr_t'
      write(3,*)'tidezoom_jend_t'
      write(3,*)'tidezoom_istr_u'
      write(3,*)'tidezoom_iend_u'
      write(3,*)'tidezoom_jstr_u'
      write(3,*)'tidezoom_jend_u'
      write(3,*)'tidezoom_istr_v'
      write(3,*)'tidezoom_iend_v'
      write(3,*)'tidezoom_jstr_v'
      write(3,*)'tidezoom_jend_v'
      write(3,*)'tidefull_imax'
      write(3,*)'tidefull_jmax'
      write(3,*)'ww3_imax'
      write(3,*)'ww3_jmax'
      write(3,*)'ww3_kmax'
      write(3,*)'ww3_fmax'
      write(3,*)'ww3zoom_istr'
      write(3,*)'ww3zoom_iend'
      write(3,*)'ww3zoom_jstr'
      write(3,*)'ww3zoom_jend'
      write(3,*)'ww3full_imax'
      write(3,*)'ww3full_jmax'
      write(3,*)'jour'
      write(3,*)'kstop'
      write(3,*)'istr'
      write(3,*)'jstr'
      write(3,*)'kstr'
      write(3,*)'tstr'
      write(3,*)'bstr'
      write(3,*)'iend'
      write(3,*)'jend'
      write(3,*)'kend'
      write(3,*)'tend'
      write(3,*)'bend'
      write(3,*)'dimend'
      write(3,*)'kpvwave'
      write(3,*)'give_chanel9'
      write(3,*)'i0'
      write(3,*)'j0'
      write(3,*)'iteration2d'
      write(3,*)'iteration2d_begin'
      write(3,*)'iteration2d_max_now'
      write(3,*)'iteration2d_max_bef'
      write(3,*)'i2dh'
      write(3,*)'l1'
      write(3,*)'l1_sca'
      write(3,*)'l1_vec'
      write(3,*)'l2'
      write(3,*)'l2_sca'
      write(3,*)'l2_vec'
      write(3,*)'l3'
      write(3,*)'l3_sca'
      write(3,*)'l3_vec'
      write(3,*)'len1'
      write(3,*)'nc'
      write(3,*)'nc1'
      write(3,*)'itimets'
      write(3,*)'itimebio'
      write(3,*)'iadvec_ts_hor'
      write(3,*)'iadvec_ts_hor_upwind'
      write(3,*)'iadvec_ts_hor_quickest'
      write(3,*)'iadvec_ts_hor_quickest2'
      write(3,*)'iadvec_ts_ver'
      write(3,*)'iadvec_ts_ver_quickest'
      write(3,*)'iadvec_ts_ver_c2'
      write(3,*)'iadvec_ts_ver_quickest2'
      write(3,*)'iturbulence'
      write(3,*)'istreamf'
      write(3,*)'itime'
      write(3,*)'kts'
      write(3,*)'kuv'
      write(3,*)'ko'
      write(3,*)'jm1'
      write(3,*)'jm2'
      write(3,*)'jp1'
      write(3,*)'jp2'
      write(3,*)'im1'
      write(3,*)'im2'
      write(3,*)'ip1'
      write(3,*)'ip2'
      write(3,*)'iwind'
      write(3,*)'ip'
      write(3,*)'im'
      write(3,*)'jp'
      write(3,*)'jm'
      write(3,*)'kp'
      write(3,*)'km'
      write(3,*)'nbcanal'
      write(3,*)'gridtype1'
      write(3,*)'gridtype2'
      write(3,*)'point1'
      write(3,*)'point2'
      write(3,*)'sender'
      write(3,*)'receiver'
      write(3,*)'ipnoc'
      write(3,*)'i1'
      write(3,*)'i2'
      write(3,*)'i3'
      write(3,*)'i4'
      write(3,*)'i5'
      write(3,*)'i6'
      write(3,*)'i7'
      write(3,*)'i8'
      write(3,*)'i9'
      write(3,*)'i10'
      write(3,*)'i11'
      write(3,*)'j1'
      write(3,*)'j2'
      write(3,*)'j3'
      write(3,*)'j4'
      write(3,*)'j5'
      write(3,*)'j6'
      write(3,*)'j7'
      write(3,*)'j8'
      write(3,*)'j9'
      write(3,*)'j10'
      write(3,*)'j11'
      write(3,*)'k0'
      write(3,*)'k1'
      write(3,*)'k2'
      write(3,*)'k3'
      write(3,*)'k4'
      write(3,*)'k5'
      write(3,*)'k6'
      write(3,*)'k7'
      write(3,*)'k8'
      write(3,*)'k9'
      write(3,*)'kr'
      write(3,*)'kp1'
      write(3,*)'kp2'
      write(3,*)'km1'
      write(3,*)'km2'
      write(3,*)'key'
      write(3,*)'flag'
      write(3,*)'nbomax'
      write(3,*)'nbobuffermax'
      write(3,*)'nbobuffermax_c'
      write(3,*)'flag_stop'
      write(3,*)'itest'
      write(3,*)'itest1'
      write(3,*)'itest2'
      write(3,*)'itest3'
      write(3,*)'ioption'
      write(3,*)'nsmooth'
      write(3,*)'nriver'
      write(3,*)'ian0'
      write(3,*)'imois0'
      write(3,*)'ijour0'
      write(3,*)'iheure0'
      write(3,*)'iminute0'
      write(3,*)'iseconde0'
      write(3,*)'iref'
      write(3,*)'jref'
      write(3,*)'lname1'
      write(3,*)'lname2'
      write(3,*)'lname3'
      write(3,*)'lname4'
      write(3,*)'kmode'
      write(3,*)'kmodemax'
      write(3,*)'fgrid_or_wgrid'
      write(3,*)'fgrid_case'
      write(3,*)'wgrid_case'
      write(3,*)'typegrid'
      write(3,*)'typegrid_monopole'
      write(3,*)'typegrid_file'
      write(3,*)'typegrid_bipole'
      write(3,*)'istar'
      write(3,*)'jstar'
      write(3,*)'istop'
      write(3,*)'jstop'
      write(3,*)'isend'
      write(3,*)'jsend'
      write(3,*)'irecv'
      write(3,*)'jrecv'
      write(3,*)'izoomin'
      write(3,*)'izoomax'
      write(3,*)'jzoomin'
      write(3,*)'jzoomax'
      write(3,*)'iairsea'
      write(3,*)'ialbedo'
      write(3,*)'airseaoption'
      write(3,*)'iobc_f'
      write(3,*)'iobc_wv'
      write(3,*)'iobc_ogcm'
      write(3,*)'iobc_lr'
      write(3,*)'obc_option'
      write(3,*)'iobc_demo_wv'
      write(3,*)'iarchive'
      write(3,*)'imodeltrc'
      write(3,*)'imodelbio'
      write(3,*)'multiple'
      write(3,*)'kvarmax'
      write(3,*)'igesig'
      write(3,*)'isigfile'
      write(3,*)'ksecu'
      write(3,*)'kmaxtide'
      write(3,*)'kmaxtidep1'
      write(3,*)'kminserie'
      write(3,*)'nzctdmax'
      write(3,*)'nsctdmax'
      write(3,*)'ktctdmin'
      write(3,*)'ktctdmax'
      write(3,*)'kdtk_out'
      write(3,*)'k10'
      write(3,*)'k11'
      write(3,*)'nbinco'
      write(3,*)'nbequa'
      write(3,*)'nbsparse'
      write(3,*)'kland'
      write(3,*)'tidestep2'
      write(3,*)'tidestep3'
      write(3,*)'ncfilm_max'
      write(3,*)'in_out_tide'
      write(3,*)'kmax_dof'
      write(3,*)'k_in'
      write(3,*)'k_out'
      write(3,*)'used_unused_dom'
      write(3,*)'mergebathy_sponge'
      write(3,*)'rankhmax'
      write(3,*)'rankhmin'
      write(3,*)'id_tem'
      write(3,*)'id_tem2'
      write(3,*)'id_dtem'
      write(3,*)'id_sal'
      write(3,*)'id_sal2'
      write(3,*)'id_rhp'
      write(3,*)'id_rhop'
      write(3,*)'id_rhom'
      write(3,*)'id_rhf'
      write(3,*)'id_rhc'
      write(3,*)'id_rhpa'
      write(3,*)'id_rhpb'
      write(3,*)'id_z'
      write(3,*)'id_prs'
      write(3,*)'id_now'
      write(3,*)'id_aft'
      write(3,*)'id_eost'
      write(3,*)'id_eoss'
      write(3,*)'id_ssh'
      write(3,*)'id_breaker'
      write(3,*)'id_dz'
      write(3,*)'id_zt'
      write(3,*)'id_tdiv'
      write(3,*)'id_sdiv'
      write(3,*)'id_udiv'
      write(3,*)'id_vdiv'
      write(3,*)'id_bdiv'
      write(3,*)'id_bnegdif'
      write(3,*)'id_biobefo'
      write(3,*)'id_bioaftr'
      write(3,*)'id_u_now'
      write(3,*)'id_v_now'
      write(3,*)'id_u_rot'
      write(3,*)'id_v_rot'
      write(3,*)'id_ffreq'
      write(3,*)'id_coriolis1'
      write(3,*)'id_coriolis2'
      write(3,*)'id_flx'
      write(3,*)'id_fly'
      write(3,*)'id_uflx'
      write(3,*)'id_ufly'
      write(3,*)'id_vflx'
      write(3,*)'id_vfly'
      write(3,*)'id_tcn'
      write(3,*)'id_scn'
      write(3,*)'id_gradssh'
      write(3,*)'id_ofactort'
      write(3,*)'id_ofactors'
      write(3,*)'id_hybcoefu'
      write(3,*)'id_hybcoefv'
      write(3,*)'id_dxdydz'
      write(3,*)'id_wb'
      write(3,*)'id_prod'
      write(3,*)'id_buoy'
      write(3,*)'iDtOvRhCp'
      write(3,*)'id_ncu'
      write(3,*)'id_ncv'
      write(3,*)'id_webio'
      write(3,*)'id_varbef2'
      write(3,*)'id_kh_over_dz'
      write(3,*)'id_bihar_lim'
      write(3,*)'id_veltot'
      write(3,*)'id_velexp'
      write(3,*)'id_velimp'
      write(3,*)'id_wdrifter'
      write(3,*)'kreadgroup'
      write(3,*)'kreadgroupmax'
      write(3,*)'kread1'
      write(3,*)'kread2'
      write(3,*)'modulo_biotimestep'
      write(3,*)'obctime_bef2'
      write(3,*)'obctime_bef'
      write(3,*)'obctime_aft'
      write(3,*)'obctime_aft2'
      write(3,*)'obctime_order'
      write(3,*)'sch_imp_ts_u_loc'
      write(3,*)'sch_imp_ts_v_loc'
      write(3,*)'sch_imp_ts_u_glb'
      write(3,*)'sch_imp_ts_v_glb'
      write(3,*)'sch_imp_tke_u_loc'
      write(3,*)'sch_imp_tke_v_loc'
      write(3,*)'sch_imp_tke_u_glb'
      write(3,*)'sch_imp_tke_v_glb'
      write(3,*)'looplimit_hor'
      write(3,*)'ihybsig'
      write(3,*)'nhybsig'
      write(3,*)'tke_surf'
      write(3,*)'grh_out_mi'
      write(3,*)'nest_onoff_in'
      write(3,*)'nest_onoff_demo'
      write(3,*)'nest_full_in'
      write(3,*)'nest_full_out'
      write(3,*)'nest_onoff_out'
      write(3,*)'ncmin_airsea'
      write(3,*)'ncmin_river'
      write(3,*)'eos_author'
      write(3,*)'eos_comprs'
      write(3,*)'eos_linear'
      write(3,*)'obcfreeorfix'
      write(3,*)'iwve'
      write(3,*)'wave_obc_type'
      write(3,*)'dataperwavefile'
      write(3,*)'sp_or_db'
      write(3,*)'kbu'
      write(3,*)'kbumax'
      write(3,*)'kbumax_glb'
      write(3,*)'n_element'
      write(3,*)'relaxtype_ts'
      write(3,*)'obctype_ts'
      write(3,*)'obctype_p'
      write(3,*)'restart_file_y_or_n'
      write(3,*)'ncid'
      write(3,*)'x_xhl_dim'
      write(3,*)'y_xhl_dim'
      write(3,*)'z_xhl_dim'
      write(3,*)'x_yhl_dim'
      write(3,*)'y_yhl_dim'
      write(3,*)'z_yhl_dim'
      write(3,*)'x_zhl_dim'
      write(3,*)'y_zhl_dim'
      write(3,*)'z_zhl_dim'
      write(3,*)'x_zl_dim'
      write(3,*)'y_zl_dim'
      write(3,*)'z_zl_dim'
      write(3,*)'time_dim'
      write(3,*)'i_t_dim'
      write(3,*)'j_t_dim'
      write(3,*)'k_t_dim'
      write(3,*)'i_w_dim'
      write(3,*)'j_w_dim'
      write(3,*)'k_w_dim'
      write(3,*)'i_u_dim'
      write(3,*)'j_u_dim'
      write(3,*)'k_u_dim'
      write(3,*)'i_v_dim'
      write(3,*)'j_v_dim'
      write(3,*)'k_v_dim'
      write(3,*)'i_f_dim'
      write(3,*)'j_f_dim'
      write(3,*)'k_f_dim'
      write(3,*)'dayindex_size'
      write(3,*)'tide_year_min'
      write(3,*)'airsea_year_min'
      write(3,*)'tide_year_max'
      write(3,*)'airsea_year_max'
      write(3,*)'wave_year_min'
      write(3,*)'irelaxsst'
      write(3,*)'removetide'
      write(3,*)'ofl_rotation'
      write(3,*)'ofl_rhp'
      write(3,*)'ofl_tke'
      write(3,*)'ofl_surflux'
      write(3,*)'ale_selected'
      write(3,*)'flag_asselin'
      write(3,*)'freq'
      write(3,*)'dirw'
      write(3,*)'dirw_beg'
      write(3,*)'dirw_end'
      write(3,*)'freq_beg'
      write(3,*)'freq_end'
      write(3,*)'year_now'
      write(3,*)'month_now'
      write(3,*)'day_now'
      write(3,*)'hour_now'
      write(3,*)'minute_now'
      write(3,*)'second_now'
      write(3,*)'nd_send_est'
      write(3,*)'nd_send_ouest'
      write(3,*)'nd_send_nord'
      write(3,*)'nd_send_sud'
      write(3,*)'nd_send_out'
      write(3,*)'nd_recv_est'
      write(3,*)'nd_recv_ouest'
      write(3,*)'nd_recv_nord'
      write(3,*)'nd_recv_sud'
      write(3,*)'nd_send_sudouest'
      write(3,*)'nd_send_sudest'
      write(3,*)'nd_send_nordouest'
      write(3,*)'nd_send_nordest'
      write(3,*)'nd_recv_sudouest'
      write(3,*)'nd_recv_sudest'
      write(3,*)'nd_recv_nordouest'
      write(3,*)'nd_recv_nordest'
      write(3,*)'initial_main_status'
      write(3,*)'offline_init_status'
      write(3,*)'meteo_sealand_mask'
      write(3,*)'grid_i0'
      write(3,*)'grid_j0'
      write(3,*)'ifb'
      write(3,*)'dbefore'
      write(3,*)'dnow'
      write(3,*)'dafter'
      write(3,*)'dim_airsea'
      write(3,*)'rhp_zavr_xy'
      write(3,*)'ssr_id'
      write(3,*)'ir_id'
      write(3,*)'rain_id'
      write(3,*)'t2m_id'
      write(3,*)'t0m_id'
      write(3,*)'abl_id'
      write(3,*)'dp2m_id'
      write(3,*)'u10m_id'
      write(3,*)'v10m_id'
      write(3,*)'u100m_id'
      write(3,*)'v100m_id'
      write(3,*)'p0m_id'
      write(3,*)'ustrs_id'
      write(3,*)'vstrs_id'
      write(3,*)'slhf_id'
      write(3,*)'netir_id'
      write(3,*)'sshf_id'
      write(3,*)'t_flux_cumul'
      write(3,*)'s_flux_cumul'
      write(3,*)'cumuldeltaflux'
      write(3,*)'som0'
      write(3,*)'som2'
      write(3,*)'ssh_reservoir'
      write(3,*)'idealflux'
      write(3,*)'sum0'
      write(3,*)'sum1'
      write(3,*)'sum2'
      write(3,*)'sum3'
      write(3,*)'sum4'
      write(3,*)'sum5'
      write(3,*)'sum6'
      write(3,*)'sum7'
      write(3,*)'sum8'
      write(3,*)'sum9'
      write(3,*)'sum10'
      write(3,*)'sum11'
      write(3,*)'sum12'
      write(3,*)'sum13'
      write(3,*)'time0'
      write(3,*)'time1'
      write(3,*)'time2'
      write(3,*)'small2'
      write(3,*)'x1_r8'
      write(3,*)'x2_r8'
      write(3,*)'x3_r8'
      write(3,*)'x4_r8'
      write(3,*)'sum0glb'
      write(3,*)'sum1glb'
      write(3,*)'sum2glb'
      write(3,*)'sum3glb'
      write(3,*)'sum4glb'
      write(3,*)'sum5glb'
      write(3,*)'sum6glb'
      write(3,*)'emin'
      write(3,*)'epsmin'
      write(3,*)'elapsedtime_out'
      write(3,*)'elapsedtime_rst'
      write(3,*)'elapsedtime_now'
      write(3,*)'elapsedtime_last_writing'
      write(3,*)'elapsedtime_bef'
      write(3,*)'elapsedtime_aft'
      write(3,*)'cpu_seconds'
      write(3,*)'alpha'
      write(3,*)'eos_pgfzref'
      write(3,*)'eos_tkezref'
      write(3,*)'filvalr8'
      write(3,*)'dti_fw_i4'
      write(3,*)'elapsedtime_now_i4'
      write(3,*)'elapsedtime_now_i8'
      write(3,*)'elapsedtime_now_r16'
! jalon4 dyn_restart_check ne pas effacer
      close(3)
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !
      end subroutine dyn_restart_bounds_r
!.......................................................................

      subroutine dyn_restart_allo(txt_)
      use module_principal ; use module_offline
      implicit none
      character*1 txt_


      if(txt_=='a') then
       if(ioffline_prv/=0.and.ioffline==0)call offline_allocate('a')
      endif

      if(txt_=='d') then
       if(ioffline_prv/=0.and.ioffline==0)call offline_allocate('d')
      endif

      if(txt_=='a')return
      if(txt_=='d')return
      stop 'dyn_restart_allo Err on txt_'
      end subroutine dyn_restart_allo
