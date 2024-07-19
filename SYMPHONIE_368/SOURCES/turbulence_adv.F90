!______________________________________________________________________
! SYMPHONIE ocean model
! release 362 - last update: 04-01-23
!______________________________________________________________________
!...............................................................................
! Version date      Description des modifications
!         06-08-16: mise en service des nouvelles subroutines d'advection
!                   des champs turbulents
!         02-01-17  application des valeurs min A l'issue de la phase d'advection
!         13-01-17  modifs sur calcul du flux advectif
!         21-04-17  coordonnee s-z
!         15-05-17  Arret d'urgence si loopmax_ out of range
!         13-01-19  cumul des pivots evolutifs
! v247    01-03-19  commentaire devant une ligne test oubliEe
! v274    09-02-20  Message dans fichiers fort.xxx le 09-02-20
! v287    17-07-20  utiliser tmpdirname
!         14-08-20  obc_int_anyv3d renomme obc_mpi_anyv3d
! v349    27-06-22  i,j glb dans fichier fort.xxx
! v361    31-12-22  advection horizontale implicite
! v362    04-01-23  work_status_mpi_min_
!...............................................................................
!    _________                    .__                  .__                     ! 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................
      subroutine turbulence_adv_tkea  ! 05-01-16
      use module_principal
      use module_parallele
      implicit none
      integer ::                 &
!                id_tke_=3       &  
                 id_div_=5       & ! libre mais apres turbulence_adv_substep
                ,id_veldydz_u_=6 &
                ,id_veldxdz_v_=7 &
                ,loop_
#ifdef synopsis
       subroutinetitle='turbulence_adv_tkea'
       subroutinedescription=''
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pas de calcul de l'advection du modele 1DV
      if(flag_1dv==1) then
        tkea_w=tken_w
        return !14-02-16
      endif

! Nombre de sous-iterations advectives:
      call turbulence_adv_substep ! computes loopmaxtke, dti_fwsub (utilise id_dxdydz=5 mais ce dernier est libre ensuite)

! Attention id_veldydz_u_ ne peut etre egal A 1 ou 2 car ce sont les valeurs de id_prod et id_buoy
! qui doivent servir dans le cas du k-epsilon A deux équations
      call turbulence_adv_fluxlimit(id_veldydz_u_,id_veldxdz_v_) ! computes veldydz_u(i,j,k,id_velexp),veldxdz_v(i,j,k,id_velexp)


      do k=2,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3d(i,j,k,id_tke_)=tken_w(i,j,k) ! Note A ce stade tken_w=tkea_w
               tkea_w(i,j,k)=tken_w(i,j,k) ! Note A ce stade tken_w=tkea_w
      enddo ; enddo ; enddo


      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1
       anyv3d(i,j,k,id_div_)=0. !13-01-19
      enddo ; enddo ; enddo

      do loop_=1,loopmaxtke !*********>loopmaxtke


      do k=2,kmax 

! flux Oi
       do j=2,jmax-1 ; do i=2,imax
        xy_u(i,j,1)=0.5*( veldydz_u(i,j,k,id_velexp) *(tkea_w(i,j,k)+tkea_w(i-1,j,k)) & ! advflux_u(i  ,j)
                     -abs(veldydz_u(i,j,k,id_velexp))*(tkea_w(i,j,k)-tkea_w(i-1,j,k)) ) ! difflux_u(i  ,j)
       enddo         ; enddo


! Advection partielle Oi
       do j=2,jmax-1 ; do i=2,imax-1

! cumul pivot (placE avant le calcul de l'advection sinon tkea_w(i,j,k) est modifiE)
        anyv3d(i,j,k,id_div_)= & !13-01-19
        anyv3d(i,j,k,id_div_)+dti_fwsub*tkea_w(i,j,k)*(veldydz_u(i+1,j,k,id_velexp)-veldydz_u(i,j,k,id_velexp))

!       tkea_w(i,j,k)=                          &
        tkea_w(i,j,k)=(1-mask_vqs_tke_w(i,j,k))*tkea_w(i,j,k)+mask_vqs_tke_w(i,j,k)*( & !-CALCUL->
        tkea_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

           dti_fwsub*( &             !pmx>
     
                       xy_u(i+1,j,1)-xy_u(i,j,1) &
              -tkea_w(i,j,k)*(veldydz_u(i+1,j,k,id_velexp)-veldydz_u(i,j,k,id_velexp)) & !13-01-19

                     )/dxdy_t(i,j) & !pmx>


                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     ) &
                                                                                  )   !-CALCUL->

       enddo ; enddo

      enddo !k loop

! Lignes commentee depuis le tableau mask_vqs_tke_w !28-12-22
! A l'interieur de la couche merged pas d'advection:
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=2,kmerged_t(i,j)-1
!        tkea_w(i,j,k)=tken_w(i,j,k) 
!       enddo
!      enddo         ; enddo
!     endif                          !pmxpmx> !18-02-17

!     call turbulence_adv_surfavr(id_tke_) !13-01-17! Moyenne sur couches collEes
!     call obc_mpi_anyv3d(1,id_tke_,'z0')
      call obc_turbulence_tkea

      do k=2,kmax 

! flux Oj
       do j=2,jmax ; do i=2,imax-1
        xy_v(i,j,1)=0.5*( veldxdz_v(i,j,k,id_velexp) *(tkea_w(i,j,k)+tkea_w(i,j-1,k)) & ! advflux_v(i  ,j)
                     -abs(veldxdz_v(i,j,k,id_velexp))*(tkea_w(i,j,k)-tkea_w(i,j-1,k)) ) ! difflux_v(i  ,j)
       enddo         ; enddo

! Advection partielle Oj:
       do j=2,jmax-1 ; do i=2,imax-1

! cumul pivot (placE avant le calcul de l'advection sinon tkea_w(i,j,k) est modifiE)
        anyv3d(i,j,k,id_div_)= & !13-01-19
        anyv3d(i,j,k,id_div_)+dti_fwsub*tkea_w(i,j,k)*(veldxdz_v(i,j+1,k,id_velexp)-veldxdz_v(i,j,k,id_velexp))

        tkea_w(i,j,k)=(1-mask_vqs_tke_w(i,j,k))*tkea_w(i,j,k)+mask_vqs_tke_w(i,j,k)*( & !-CALCUL->
        tkea_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

           dti_fwsub*( &             !pmx>
     
                       xy_v(i,j+1,1)-xy_v(i,j,1) &
              -tkea_w(i,j,k)*(veldxdz_v(i,j+1,k,id_velexp)-veldxdz_v(i,j,k,id_velexp)) & !13-01-19

                     )/dxdy_t(i,j) & !pmx>

                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     )  &
                                                                                  )   !-CALCUL->

       enddo ; enddo

      enddo !k loop

! Lignes commentee depuis le tableau mask_vqs_tke_w !28-12-22
! A l'interieur de la couche merged pas d'advection:
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=2,kmerged_t(i,j)-1
!        tkea_w(i,j,k)=tken_w(i,j,k) 
!       enddo
!      enddo         ; enddo
!     endif                          !pmxpmx> !18-02-17


!     call turbulence_adv_surfavr(id_tke_) !13-01-17! Moyenne sur couches collEes
!     call obc_mpi_anyv3d(1,id_tke_,'z0')
      call obc_turbulence_tkea


      enddo                 !*********>loopmaxtke

      if(sch_imp_tke_u_glb==1.or.sch_imp_tke_v_glb==1)call turbulence_adv_tkea_impli(id_div_)

       do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1
        tkea_w(i,j,k)=max(emin,                 & !02-01-17
        tkea_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

               anyv3d(i,j,k,id_div_)/dxdy_t(i,j) & !13-01-19

                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     )     )
       enddo ; enddo ; enddo

      end subroutine turbulence_adv_tkea

!.........................................................................................

      subroutine turbulence_adv_tken  ! 05-01-16
      use module_principal
      use module_parallele
      implicit none
      integer ::                 &
!                id_tke_=3       &  
                 id_div_=5       & ! libre mais apres turbulence_adv_substep
                ,id_veldydz_u_=6 &
                ,id_veldxdz_v_=7 &
                ,loop_
#ifdef synopsis
       subroutinetitle='turbulence_adv_tken'
       subroutinedescription=''
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pas de calcul de l'advection du modele 1DV
      if(flag_1dv==1)return

! Nombre de sous-iterations advectives:
      call turbulence_adv_substep ! computes loopmaxtke, dti_fwsub (utilise id_dxdydz=5 mais ce dernier est libre ensuite)

! Attention id_veldydz_u_ ne peut etre egal A 1 ou 2 car ce sont les valeurs de id_prod et id_buoy
! qui doivent servir dans le cas du k-epsilon A deux équations
      call turbulence_adv_fluxlimit(id_veldydz_u_,id_veldxdz_v_) ! computes veldydz_u(i,j,k,id_velexp),veldxdz_v(i,j,k,id_velexp)

      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1
       anyv3d(i,j,k,id_div_)=0. !13-01-19
      enddo ; enddo ; enddo

      do loop_=1,loopmaxtke !*********>loopmaxtke


      do k=2,kmax 

! flux Oi
       do j=2,jmax-1 ; do i=2,imax
        xy_u(i,j,1)=0.5*( veldydz_u(i,j,k,id_velexp) *(tken_w(i,j,k)+tken_w(i-1,j,k)) & ! advflux_u(i  ,j)
                     -abs(veldydz_u(i,j,k,id_velexp))*(tken_w(i,j,k)-tken_w(i-1,j,k)) ) ! difflux_u(i  ,j)
       enddo         ; enddo


! Advection partielle Oi
       do j=2,jmax-1 ; do i=2,imax-1

! cumul pivot (placE avant le calcul de l'advection sinon tken_w(i,j,k) est modifiE)
        anyv3d(i,j,k,id_div_)= & !13-01-19
        anyv3d(i,j,k,id_div_)+dti_fwsub*tken_w(i,j,k)*(veldydz_u(i+1,j,k,id_velexp)-veldydz_u(i,j,k,id_velexp))

!       tken_w(i,j,k)=                          &
        tken_w(i,j,k)=(1-mask_vqs_tke_w(i,j,k))*tken_w(i,j,k)+mask_vqs_tke_w(i,j,k)*( & !-CALCUL->
        tken_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

           dti_fwsub*( &             !pmx>
     
                       xy_u(i+1,j,1)-xy_u(i,j,1) &
              -tken_w(i,j,k)*(veldydz_u(i+1,j,k,id_velexp)-veldydz_u(i,j,k,id_velexp)) & !13-01-19

                     )/dxdy_t(i,j) & !pmx>


                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     ) &
                                                                                  )   !-CALCUL->

       enddo ; enddo

      enddo !k loop

! Lignes commentee depuis le tableau mask_vqs_tke_w !28-12-22
! A l'interieur de la couche merged pas d'advection:
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=2,kmerged_t(i,j)-1
!        tken_w(i,j,k)=tken_w(i,j,k) 
!       enddo
!      enddo         ; enddo
!     endif                          !pmxpmx> !18-02-17

!     call turbulence_adv_surfavr(id_tke_) !13-01-17! Moyenne sur couches collEes
!     call obc_mpi_anyv3d(1,id_tke_,'z0')
      call obc_turbulence_tken

      do k=2,kmax 

! flux Oj
       do j=2,jmax ; do i=2,imax-1
        xy_v(i,j,1)=0.5*( veldxdz_v(i,j,k,id_velexp) *(tken_w(i,j,k)+tken_w(i,j-1,k)) & ! advflux_v(i  ,j)
                     -abs(veldxdz_v(i,j,k,id_velexp))*(tken_w(i,j,k)-tken_w(i,j-1,k)) ) ! difflux_v(i  ,j)
       enddo         ; enddo

! Advection partielle Oj:
       do j=2,jmax-1 ; do i=2,imax-1

! cumul pivot (placE avant le calcul de l'advection sinon tken_w(i,j,k) est modifiE)
        anyv3d(i,j,k,id_div_)= & !13-01-19
        anyv3d(i,j,k,id_div_)+dti_fwsub*tken_w(i,j,k)*(veldxdz_v(i,j+1,k,id_velexp)-veldxdz_v(i,j,k,id_velexp))

        tken_w(i,j,k)=(1-mask_vqs_tke_w(i,j,k))*tken_w(i,j,k)+mask_vqs_tke_w(i,j,k)*( & !-CALCUL->
        tken_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

           dti_fwsub*( &             !pmx>
     
                       xy_v(i,j+1,1)-xy_v(i,j,1) &
              -tken_w(i,j,k)*(veldxdz_v(i,j+1,k,id_velexp)-veldxdz_v(i,j,k,id_velexp)) & !13-01-19

                     )/dxdy_t(i,j) & !pmx>

                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     )  &
                                                                                  )   !-CALCUL->

       enddo ; enddo

      enddo !k loop

! Lignes commentee depuis le tableau mask_vqs_tke_w !28-12-22
! A l'interieur de la couche merged pas d'advection:
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=2,kmerged_t(i,j)-1
!        tken_w(i,j,k)=tken_w(i,j,k) 
!       enddo
!      enddo         ; enddo
!     endif                          !pmxpmx> !18-02-17


!     call turbulence_adv_surfavr(id_tke_) !13-01-17! Moyenne sur couches collEes
!     call obc_mpi_anyv3d(1,id_tke_,'z0')
      call obc_turbulence_tken


      enddo                 !*********>loopmaxtke

      if(sch_imp_tke_u_glb==1.or.sch_imp_tke_v_glb==1)call turbulence_adv_tken_impli(id_div_)

       do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1
        tken_w(i,j,k)=max(emin,                 & !02-01-17
        tken_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

               anyv3d(i,j,k,id_div_)/dxdy_t(i,j) & !13-01-19

                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     )     )
       enddo ; enddo ; enddo

      end subroutine turbulence_adv_tken

!.........................................................................................

      subroutine turbulence_adv_epsa  ! 05-01-16
      use module_principal
      use module_parallele
      implicit none
      integer ::                 &
!                id_eps_=3       &  
                 id_div_=5       & ! libre mais apres turbulence_adv_substep
                ,id_veldydz_u_=6 &
                ,id_veldxdz_v_=7 &
                ,loop_
#ifdef synopsis
       subroutinetitle='turbulence_adv_epsa'
       subroutinedescription=''
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pas de calcul de l'advection du modele 1DV
      if(flag_1dv==1) then
        epsa_w=epsn_w
        return !14-02-16
      endif

! Nombre de sous-iterations advectives:
!     call turbulence_adv_substep ! computes loopmaxtke, dti_fwsub (utilise id_dxdydz=5 mais ce dernier est libre ensuite)
! Attention id_veldydz_u_ ne peut etre egal A 1 ou 2 car ce sont les valeurs de id_prod et id_buoy
! qui doivent servir dans le cas du k-epsilon A deux équations
!     call turbulence_adv_fluxlimit(id_veldydz_u_,id_veldxdz_v_) ! computes veldydz_u(i,j,k,id_velexp),veldxdz_v(i,j,k,id_velexp)


      do k=2,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3d(i,j,k,id_eps_)=epsn_w(i,j,k) ! Note A ce stade epsn_w=epsa_w
               epsa_w(i,j,k)=epsn_w(i,j,k) ! Note A ce stade epsn_w=epsa_w
      enddo ; enddo ; enddo


      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1
       anyv3d(i,j,k,id_div_)=0. !13-01-19
      enddo ; enddo ; enddo

      do loop_=1,loopmaxtke !*********>loopmaxtke


      do k=2,kmax 

! flux Oi
       do j=2,jmax-1 ; do i=2,imax
        xy_u(i,j,1)=0.5*( veldydz_u(i,j,k,id_velexp) *(epsa_w(i,j,k)+epsa_w(i-1,j,k)) & ! advflux_u(i  ,j)
                     -abs(veldydz_u(i,j,k,id_velexp))*(epsa_w(i,j,k)-epsa_w(i-1,j,k)) ) ! difflux_u(i  ,j)
       enddo         ; enddo


! Advection partielle Oi
       do j=2,jmax-1 ; do i=2,imax-1

! cumul pivot (placE avant le calcul de l'advection sinon epsa_w(i,j,k) est modifiE)
        anyv3d(i,j,k,id_div_)= & !13-01-19
        anyv3d(i,j,k,id_div_)+dti_fwsub*epsa_w(i,j,k)*(veldydz_u(i+1,j,k,id_velexp)-veldydz_u(i,j,k,id_velexp))

!       epsa_w(i,j,k)=                          &
        epsa_w(i,j,k)=(1-mask_vqs_tke_w(i,j,k))*epsa_w(i,j,k)+mask_vqs_tke_w(i,j,k)*( & !-CALCUL->
        epsa_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

           dti_fwsub*( &             !pmx>
     
                       xy_u(i+1,j,1)-xy_u(i,j,1) &
              -epsa_w(i,j,k)*(veldydz_u(i+1,j,k,id_velexp)-veldydz_u(i,j,k,id_velexp)) & !13-01-19

                     )/dxdy_t(i,j) & !pmx>


                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     ) &
                                                                                  )   !-CALCUL->

       enddo ; enddo

      enddo !k loop

! Lignes commentee depuis le tableau mask_vqs_tke_w !28-12-22
! A l'interieur de la couche merged pas d'advection:
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=2,kmerged_t(i,j)-1
!        epsa_w(i,j,k)=epsn_w(i,j,k) 
!       enddo
!      enddo         ; enddo
!     endif                          !pmxpmx> !18-02-17

!     call turbulence_adv_surfavr(id_eps_) !13-01-17! Moyenne sur couches collEes
!     call obc_mpi_anyv3d(1,id_eps_,'z0')
      call obc_turbulence_epsa

      do k=2,kmax 

! flux Oj
       do j=2,jmax ; do i=2,imax-1
        xy_v(i,j,1)=0.5*( veldxdz_v(i,j,k,id_velexp) *(epsa_w(i,j,k)+epsa_w(i,j-1,k)) & ! advflux_v(i  ,j)
                     -abs(veldxdz_v(i,j,k,id_velexp))*(epsa_w(i,j,k)-epsa_w(i,j-1,k)) ) ! difflux_v(i  ,j)
       enddo         ; enddo

! Advection partielle Oj:
       do j=2,jmax-1 ; do i=2,imax-1

! cumul pivot (placE avant le calcul de l'advection sinon epsa_w(i,j,k) est modifiE)
        anyv3d(i,j,k,id_div_)= & !13-01-19
        anyv3d(i,j,k,id_div_)+dti_fwsub*epsa_w(i,j,k)*(veldxdz_v(i,j+1,k,id_velexp)-veldxdz_v(i,j,k,id_velexp))

        epsa_w(i,j,k)=(1-mask_vqs_tke_w(i,j,k))*epsa_w(i,j,k)+mask_vqs_tke_w(i,j,k)*( & !-CALCUL->
        epsa_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

           dti_fwsub*( &             !pmx>
     
                       xy_v(i,j+1,1)-xy_v(i,j,1) &
              -epsa_w(i,j,k)*(veldxdz_v(i,j+1,k,id_velexp)-veldxdz_v(i,j,k,id_velexp)) & !13-01-19

                     )/dxdy_t(i,j) & !pmx>

                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     )  &
                                                                                  )   !-CALCUL->

       enddo ; enddo

      enddo !k loop

! Lignes commentee depuis le tableau mask_vqs_tke_w !28-12-22
! A l'interieur de la couche merged pas d'advection:
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=2,kmerged_t(i,j)-1
!        epsa_w(i,j,k)=epsn_w(i,j,k) 
!       enddo
!      enddo         ; enddo
!     endif                          !pmxpmx> !18-02-17


!     call turbulence_adv_surfavr(id_eps_) !13-01-17! Moyenne sur couches collEes
!     call obc_mpi_anyv3d(1,id_eps_,'z0')
      call obc_turbulence_epsa


      enddo                 !*********>loopmaxtke

      if(sch_imp_tke_u_glb==1.or.sch_imp_tke_v_glb==1)call turbulence_adv_epsa_impli(id_div_)

       do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1
        epsa_w(i,j,k)=max(epsmin,                 & !02-01-17
        epsa_w(i,j,k)                           &

              -wetmask_t(i,j)*( & !ppp> ! A noter la multiplication par 0.25 qui disparait avec celle du denominateur

               anyv3d(i,j,k,id_div_)/dxdy_t(i,j) & !13-01-19

                              ) & !ppp>

                /(0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) & 
                        +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) & 
                 )     )     )
       enddo ; enddo ; enddo

      end subroutine turbulence_adv_epsa

!.........................................................................................

      subroutine turbulence_adv_substep
      use module_principal
      use module_parallele
      implicit none

! Note: j'ai verifie qu' anyv3d(:,:,:,id_dxdydz) est bien libre tout au long de
! la sequence de la fermeture turbulente.

! ETAPE 1: CALCULER ANYV3D(:,:,:,id_dxdydz)


! Standard case:
!     if(flag_merged_levels==0) then !standard case>

!      do k=2,kmax ; do j=1,jmax   ; do i=1,imax
       do k=2,kmax ; do j=0,jmax+1 ; do i=0,imax+1 ! plus grande pour l'advection implicite

                    anyv3d(i,j,k,id_dxdydz)=                         &
                    dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)  &
                                     +dz_t(i,j,k-1,1)+dz_t(i,j,k,1)) 

       enddo       ; enddo       ; enddo

! Couches merged:
! Le niveau k=kmin_w a une epaisseur egale A la demi-some de dz(kmin-1) et de
! dz(kmin), ce dernier etant A priori toujours suffisament epais pour que le calcul
! de l'advection soit viable:
!      do j=1,jmax   ; do i=1,imax
       do j=0,jmax+1 ; do i=0,imax+1 ! plus grande pour l'advection implicite
         do k=2,kmerged_t(i,j)-1
           anyv3d(i,j,k,id_dxdydz)=anyv3d(i,j,kmerged_t(i,j),id_dxdydz)
         enddo
       enddo       ; enddo


!     endif                          !standard case>
   

!..............................................
! ETAPE 2: Determiner le nombre de sous-iterations d'advection
      x1=0. ! reset loopmax decimal
      do k=2,kmax
       do j=1,jmax ; do i=1,imax

        xy_t(i,j,1)=0.5*dti_fw*wetmask_t(i,j)                        &
                   /anyv3d(i,j,k,id_dxdydz) !01-12-16

       enddo ; enddo
       do j=2,jmax-1 ; do i=2,imax
! loopmax decimal:
       x1=max(x1,                                               &
       abs(veldydz_u(i,j,k,1)+veldydz_u(i,j,k-1,1))*max(xy_t(i,j,1),xy_t(i-1,j,1)))
       enddo ; enddo
       do j=2,jmax ; do i=2,imax-1
! loopmax decimal:
       x1=max(x1,                                               &
       abs(veldxdz_v(i,j,k,1)+veldxdz_v(i,j,k-1,1))*max(xy_t(i,j,1),xy_t(i,j-1,1)))
       enddo ; enddo
      enddo ! k loop

      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)

! Plus tard on testera cette possibilite de limitation du flux:
!     loopmaxtke=min(int(x2)+1,4)
!     x2=2.1 ! bidouille pmx ! ligne commentee le 01-03-19
      loopmaxtke=    int(x2)+1   
      dti_fwsub=dti_fw/loopmaxtke

      if(par%rank==0)then
       open(unit=3,file=trim(tmpdirname)//'dti_tke',position='append') !17-07-20
        write(3,*)real(elapsedtime_now/86400.),loopmaxtke,real(x2),real(dti_fwsub)
       close(3)
      endif

      if(loopmaxtke>50.or.loopmaxtke<0) then !>>>

! Message dans fichiers fort.xxx le 09-02-20
       do k=2,kmax
        do j=1,jmax ; do i=1,imax
         xy_t(i,j,1)=0.5*dti_fw*wetmask_t(i,j)/anyv3d(i,j,k,id_dxdydz)
        enddo ; enddo
        do j=2,jmax-1 ; do i=2,imax
         x1=abs(veldydz_u(i,j,k,1)+veldydz_u(i,j,k-1,1))*max(xy_t(i,j,1),xy_t(i-1,j,1))
         if(x1>49) then !>>>
           write(10+par%rank,*)'.................'
           write(10+par%rank,*)'x1=',x1
           write(10+par%rank,*)'par%rank',par%rank
           write(10+par%rank,*)'i,j,k',i,j,k
           write(10+par%rank,*)'i,j glob',i+par%timax(1),j+par%tjmax(1) !27-06-22
           write(10+par%rank,*)'veldydz_u(i,j,k  ,1)',veldydz_u(i,j,k  ,1)
           write(10+par%rank,*)'vel_u(i,j,k  ,1:2)',vel_u(i,j,k  ,1:2)
           write(10+par%rank,*)'veldydz_u(i,j,k-1,1)',veldydz_u(i,j,k-1,1)
           write(10+par%rank,*)'vel_u(i,j,k-1,1:2)',vel_u(i,j,k-1,1:2)
           write(10+par%rank,*)'xy_t(i  ,j,1)',xy_t(i  ,j,1)
           write(10+par%rank,*)'xy_t(i-1,j,1)',xy_t(i-1,j,1)
           write(10+par%rank,*)'anyv3d(i  ,j,k,id_dxdydz)',anyv3d(i  ,j,k,id_dxdydz)
           write(10+par%rank,*)'dz_t(i  ,j,k,1:2)',dz_t(i  ,j,k,1:2)
           write(10+par%rank,*)'anyv3d(i-1,j,k,id_dxdydz)',anyv3d(i-1,j,k,id_dxdydz)
           write(10+par%rank,*)'dz_t(i-1,j,k,1:2)',dz_t(i-1,j,k,1:2)
         endif          !>>>
        enddo ; enddo
        do j=2,jmax ; do i=2,imax-1
         x1=abs(veldxdz_v(i,j,k,1)+veldxdz_v(i,j,k-1,1))*max(xy_t(i,j,1),xy_t(i,j-1,1))
         if(x1>49) then !>>>
           write(10+par%rank,*)'.................'
           write(10+par%rank,*)'x1=',x1
           write(10+par%rank,*)'par%rank',par%rank
           write(10+par%rank,*)'i,j,k',i,j,k
           write(10+par%rank,*)'i,j glob',i+par%timax(1),j+par%tjmax(1) !27-06-22
           write(10+par%rank,*)'veldxdz_v(i,j,k  ,1)',veldxdz_v(i,j,k  ,1)
           write(10+par%rank,*)'vel_v(i,j,k  ,1:2)',vel_v(i,j,k  ,1:2)
           write(10+par%rank,*)'veldxdz_v(i,j,k-1,1)',veldxdz_v(i,j,k-1,1)
           write(10+par%rank,*)'vel_v(i,j,k-1,1:2)',vel_v(i,j,k-1,1:2)
           write(10+par%rank,*)'xy_t(i,j  ,1)',xy_t(i,j  ,1)
           write(10+par%rank,*)'xy_t(i,j-1,1)',xy_t(i,j-1,1)
           write(10+par%rank,*)'anyv3d(i,j  ,k,id_dxdydz)',anyv3d(i,j  ,k,id_dxdydz)
           write(10+par%rank,*)'dz_t(i  ,j,k,1:2)',dz_t(i  ,j,k,1:2)
           write(10+par%rank,*)'anyv3d(i,j-1,k,id_dxdydz)',anyv3d(i,j-1,k,id_dxdydz)
           write(10+par%rank,*)'dz_t(i,j-1,k,1:2)',dz_t(i,j-1,k,1:2)
         endif          !>>>
        enddo ; enddo
       enddo

        call graph_out
        stop 'loopmax tke advection out of range' !15-05-17
      endif                              !>>>


!     write(6,*)'x1=',x1,loopmaxtke
!     if(iteration3d==2)stop 'coco'

!      write(6,*)'loopmaxtke',loopmaxtke
!      stop 'jojo'

      end subroutine turbulence_adv_substep

!.................................................................

      subroutine turbulence_adv_fluxlimit(id_veldydz_u_,id_veldxdz_v_) 
      use module_principal
      use module_parallele
      implicit none
      integer id_veldydz_u_,id_veldxdz_v_

      id_velexp=0 ; id_velimp=2
      if(loopmaxtke>looplimit_hor) then !pmx>

!..............................................
! Recalculer la moyenne de veldydz_u et de veldxdz_v sur k et k-1 en bornant de sorte 
! qu'ils ne depassent pas un nombre de courant de xxx

!     do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax
      do k=2,kmax ; do j=2,jmax-1 ; do i=1,imax+1 ! plus grande pour l'advection implicite

! NOTER QUE veldydz_u(i,j,k,id_velexp) correspond A la moyenne de veldydz_u(i,j,k,1) et de veldydz_u(i,j,k-1,1)
       veldydz_u(i,j,k,id_velexp)=max(min( &

         0.5*(veldydz_u(i,j,k,1)+veldydz_u(i,j,k-1,1))                 &

        ,+looplimit_hor*min(     anyv3d(i  ,j,k,id_dxdydz)             &
                         /max(wetmask_t(i  ,j),1.e-10)                 &
                                ,anyv3d(i-1,j,k,id_dxdydz)             &
                         /max(wetmask_t(i-1,j),1.e-10)                 &
                           )/dti_fw)                                   &
        ,-looplimit_hor*min(     anyv3d(i  ,j,k,id_dxdydz)             &
                         /max(wetmask_t(i  ,j),1.e-10)                 &
                                ,anyv3d(i-1,j,k,id_dxdydz)             &
                         /max(wetmask_t(i-1,j),1.e-10)                 &
                           )/dti_fw)


! NOTER QUE veldydz_u(i,j,k,id_velimp) correspond A la moyenne de veldydz_u(i,j,k,1) et de veldydz_u(i,j,k-1,1)
       veldydz_u(i,j,k,id_velimp)=0.5*(veldydz_u(i,j,k,1)+veldydz_u(i,j,k-1,1)) &
      -veldydz_u(i,j,k,id_velexp)

      enddo ; enddo ; enddo

!     do k=2,kmax ; do j=2,jmax   ; do i=2,imax-1
      do k=2,kmax ; do j=1,jmax+1 ; do i=2,imax-1 ! plus grande pour l'advection implicite

! NOTER QUE veldydz_u(i,j,k,id_velexp) correspond A la moyenne de veldydz_u(i,j,k,1) et de veldydz_u(i,j,k-1,1)
       veldxdz_v(i,j,k,id_velexp)=max(min( &

         0.5*(veldxdz_v(i,j,k,1)+veldxdz_v(i,j,k-1,1))                 &

        ,+looplimit_hor*min(     anyv3d(i,j  ,k,id_dxdydz)             &
                         /max(wetmask_t(i,j  ),1.e-10)                 &
                                ,anyv3d(i,j-1,k,id_dxdydz)             &
                         /max(wetmask_t(i,j-1),1.e-10)                 &
                           )/dti_fw)                                   &
        ,-looplimit_hor*min(     anyv3d(i,j  ,k,id_dxdydz)             &
                         /max(wetmask_t(i,j  ),1.e-10)                 &
                                ,anyv3d(i,j-1,k,id_dxdydz)             &
                         /max(wetmask_t(i,j-1),1.e-10)                 &
                           )/dti_fw)


! NOTER QUE veldxdz_v(i,j,k,id_velimp) correspond A la moyenne de veldxdz_v(i,j,k,1) et de veldxdz_v(i,j,k-1,1)
       veldxdz_v(i,j,k,id_velimp)=0.5*(veldxdz_v(i,j,k,1)+veldxdz_v(i,j,k-1,1)) &
      -veldxdz_v(i,j,k,id_velexp)

      enddo ; enddo ; enddo

      sch_imp_tke_u_loc=0
      if(maxval(abs(veldydz_u(:,:,2:kmax,id_velimp)))>0.)sch_imp_tke_u_loc=1
      call mpi_allreduce(sch_imp_tke_u_loc,sch_imp_tke_u_glb,1,mpi_integer,mpi_max,par%comm2d,ierr)

      sch_imp_tke_v_loc=0
      if(maxval(abs(veldxdz_v(:,:,2:kmax,id_velimp)))>0.)sch_imp_tke_v_loc=1
      call mpi_allreduce(sch_imp_tke_v_loc,sch_imp_tke_v_glb,1,mpi_integer,mpi_max,par%comm2d,ierr)

      else                              !pmx>

! NOTER QUE veldydz_u(i,j,k,id_velexp) correspond A la moyenne de veldydz_u(i,j,k,1) et de veldydz_u(i,j,k-1,1)
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax
       veldydz_u(i,j,k,id_velexp)=0.5*(veldydz_u(i,j,k,1)+veldydz_u(i,j,k-1,1)) 
      enddo ; enddo ; enddo
      do k=2,kmax ; do j=2,jmax ; do i=2,imax-1
       veldxdz_v(i,j,k,id_velexp)=0.5*(veldxdz_v(i,j,k,1)+veldxdz_v(i,j,k-1,1)) 
      enddo ; enddo ; enddo

      sch_imp_tke_u_loc=0
      sch_imp_tke_v_loc=0
      sch_imp_tke_u_glb=0
      sch_imp_tke_v_glb=0

      endif                             !pmx>

      end subroutine turbulence_adv_fluxlimit

!................................................................
#ifdef bidon
      subroutine turbulence_adv_surfavr(id_)
      use module_principal
      use module_parallele
      implicit none
      integer id_

! No merged levels if flag_merged_levels=0 !13-01-17
     if(flag_merged_levels==0)return

! Moyenner anyv3d sur les couches de surface "collees" 
! sur les niveaux k en coherence avec hypotheses de la
! subroutine turbulence_adv_substep

      do j=2,jmax-1
      do i=2,imax-1
       k1=int(truekmax_t(i,j))+1
       sum1=small2
       sum2=0.
       do k=k1+1,kmax
! On choisit 4eme arg=2 car a priori coherent avec critere de conservation des traceurs
        sum1=sum1+(0.25*( dz_t(i,j,k-1, 2)+dz_t(i,j,k, 2) &
                         +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) ) )    

        sum2=sum2+(0.25*( dz_t(i,j,k-1, 2)+dz_t(i,j,k, 2) &
                         +dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) ) )*anyv3d(i,j,k,id_)
       enddo
       x1=sum2/sum1
       do k=k1,kmax
        anyv3d(i,j,k,id_)=x1
       enddo
      enddo
      enddo

      end subroutine turbulence_adv_surfavr
#endif
!................................................................
!31-12-22
      subroutine turbulence_adv_tkea_impli(id_div_)
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down & ! 0= pas d'echange avec voisin aval , 1= echange
                ,id_div_            &
                ,work_status_mpi_min_ !04-01-23



      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_div_) se poursuit
!     do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1 
!      anyv3d(i,j,k,id_div_)=0. ! reset cumul temperature*divergence(courant) !26-01-17
!     enddo ; enddo ; enddo

!...........................
! Advection axe Oi
!...........................

! Passage boucle i croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0       ! 0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0       ! 0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0       ! 0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldydz_u(2,2:jmax-1,2:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=imax du present rank est vel_u en i=2 du voisin aval
      if(maxval(veldydz_u(imax,2:jmax-1,2:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldydz_u(2:imax-1,2:jmax-1,2:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 10 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(est)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(ouest)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call turbulence_adv_tkea_impobc(1,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

!      tkea_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tkea_w(i,j,k)                              &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
!       *tkea_w(i-1,j,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tkea_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tkea_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tkea_w(i,j,k) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
        *(tkea_w(i-1,j,k)-tkea_w(i,j,k))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tkea_w(i,j,k)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   10 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tkea_w(i,j,k)*dti_fw*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle i decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont: attention la boucle va jusqu'en i=imax-1 et utilise vel_u(i+1) donc vel_u(imax)
      if(minval(veldydz_u(imax,2:jmax-1,2:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=2 du present rank est vel_u en i=imax du voisin aval
      if(minval(veldydz_u(2,2:jmax-1,2:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldydz_u(3:imax,2:jmax-1,2:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 11 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(ouest)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(est)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call turbulence_adv_tkea_impobc(2,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=imax-1,2,-1

!      tkea_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tkea_w(i,j,k)                              &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
!       *tkea_w(i+1,j,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tkea_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tkea_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tkea_w(i,j,k) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
        *(tkea_w(i+1,j,k)-tkea_w(i,j,k))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tkea_w(i,j,k)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   11 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tkea_w(i,j,k)*dti_fw*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans turbulence_adv_tkea_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur tkea_w(i,j,k)
!     call obc_mpi_tkea_w(1,2,'zb')

!...........................
! Advection axe Oj
!...........................

! Passage boucle j croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldxdz_v(2:imax-1,2,2:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=jmax du present rank est vel_v en j=2 du voisin aval
      if(maxval(veldxdz_v(2:imax-1,jmax,2:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldxdz_v(2:imax-1,2:jmax-1,2:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 20 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(nord)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(sud)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call turbulence_adv_tkea_impobc(3,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

!      tkea_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tkea_w(i,j,k)                              &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
!       *tkea_w(i,j-1,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tkea_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tkea_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tkea_w(i,j,k) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
        *(tkea_w(i,j-1,k)-tkea_w(i,j,k))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tkea_w(i,j,k)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->
   20 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tkea_w(i,j,k)*dti_fw*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle j decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont: attention la boucle va jusqu'en jmax-1 et utilise vel_v(j+1) soit vel_v(jmax)
      if(minval(veldxdz_v(2:imax-1,jmax,2:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=2 du present rank est vel_v en j=jmax du voisin aval
      if(minval(veldxdz_v(2:imax-1,2,2:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldxdz_v(2:imax-1,3:jmax,2:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 21 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(sud)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(nord)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call turbulence_adv_tkea_impobc(4,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=jmax-1,2,-1 ; do i=2,imax-1

!      tkea_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tkea_w(i,j,k)                              &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
!       *tkea_w(i,j+1,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tkea_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tkea_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tkea_w(i,j,k) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
        *(tkea_w(i,j+1,k)-tkea_w(i,j,k))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tkea_w(i,j,k)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   21 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tkea_w(i,j,k)*dti_fw*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans turbulence_adv_tkea_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur tkea_w(i,j,k)
!     call obc_mpi_tkea_w(1,2,'zb')

      end subroutine turbulence_adv_tkea_impli

!..............................................................................
!31-12-22
      subroutine turbulence_adv_tkea_impobc(case_,idvar_,send_status_,recv_status_)
      use module_principal
      use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =2
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_,loop_,idvar_,case_,send_status_,recv_status_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_

      if(.not.allocated(adv_send_est)) then !>>>
               allocate(adv_send_est  (2*kmax*jmax)) ; adv_send_est=0.
               allocate(adv_send_ouest(2*kmax*jmax)) ; adv_send_ouest=0.
               allocate(adv_recv_est  (2*kmax*jmax)) ; adv_recv_est=0.
               allocate(adv_recv_ouest(2*kmax*jmax)) ; adv_recv_ouest=0.
      endif                                    !>>>
      if(.not.allocated(adv_send_nord)) then !>>>
               allocate(adv_send_nord  (2*kmax*imax)) ; adv_send_nord=0.
               allocate(adv_send_sud   (2*kmax*imax)) ; adv_send_sud=0.
               allocate(adv_recv_nord  (2*kmax*imax)) ; adv_recv_nord=0.
               allocate(adv_recv_sud   (2*kmax*imax)) ; adv_recv_sud=0.
      endif                                    !>>>


      if(case_==1) then !---echange-Oi-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        i=imax-1
        do k=2,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_est(k0)=tkea_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_est(1)     & ! envoyé au proc Est
                ,size(adv_send_est(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(est)            &
                ,tagest_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*jmax
       k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_ouest(1)    &
               ,size(adv_recv_ouest(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(ouest)                           &
               ,tagest_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=1
       k0=0
       do k=2,kmax ; do j=2,jmax-1
        k0=k0+1
        tkea_w(i,j,k)=adv_recv_ouest(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oi-boucle-croissante--->

      if(case_==2) then !---echange-Oi-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        i=2
        do k=2,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_ouest(k0)=tkea_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_ouest(1)     & ! envoyé au proc Est
                ,size(adv_send_ouest(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(ouest)         &
                ,tagouest_                  &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*jmax
       k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_est(1)     &
               ,size(adv_recv_est(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(est)         &
               ,tagouest_                &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=imax  
       k0=0
       do k=2,kmax ; do j=2,jmax-1
        k0=k0+1
        tkea_w(i,j,k)=adv_recv_est(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oi-boucle-decroissante--->

      if(case_==3) then !---echange-Oj-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        j=jmax-1
        do k=2,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_nord(k0)=tkea_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_nord(1)     & ! envoyé au proc Est
                ,size(adv_send_nord(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(nord)            &
                ,tagnord_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*imax
       k0=  (kmax-1)*(imax-2) ! puisque do k=2,kmax ; do i=2,imax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_sud(1)    &
               ,size(adv_recv_sud(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(sud)                           &
               ,tagnord_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=1
       k0=0
       do k=2,kmax ; do i=2,imax-1
        k0=k0+1
        tkea_w(i,j,k)=adv_recv_sud(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oj-boucle-croissante--->

      if(case_==4) then !---echange-Oj-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        j=2
        do k=2,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_sud(k0)=tkea_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_sud(1)     & ! envoyé au proc Est
                ,size(adv_send_sud(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(sud)         &
                ,tagsud_                    &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*imax
       k0=  (kmax-1)*(imax-2) ! puisque do k=2,kmax ; do i=2,imax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_nord(1)     &
               ,size(adv_recv_nord(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(nord)         &
               ,tagsud_                  &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=jmax
       k0=0
       do k=2,kmax ; do i=2,imax-1
        k0=k0+1
        tkea_w(i,j,k)=adv_recv_nord(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine turbulence_adv_tkea_impobc

!....................................................
!31-12-22
      subroutine turbulence_adv_epsa_impli(id_div_)
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down & ! 0= pas d'echange avec voisin aval , 1= echange
                ,id_div_            &
                ,work_status_mpi_min_



      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_div_) se poursuit
!     do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1 
!      anyv3d(i,j,k,id_div_)=0. ! reset cumul temperature*divergence(courant) !26-01-17
!     enddo ; enddo ; enddo

!...........................
! Advection axe Oi
!...........................

! Passage boucle i croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0       ! 0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0       ! 0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0       ! 0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldydz_u(2,2:jmax-1,2:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=imax du present rank est vel_u en i=2 du voisin aval
      if(maxval(veldydz_u(imax,2:jmax-1,2:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldydz_u(2:imax-1,2:jmax-1,2:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 10 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(est)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(ouest)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call turbulence_adv_epsa_impobc(1,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

!      epsa_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *epsa_w(i,j,k)                              &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
!       *epsa_w(i-1,j,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*epsa_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       epsa_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         epsa_w(i,j,k) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
        *(epsa_w(i-1,j,k)-epsa_w(i,j,k))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*epsa_w(i,j,k)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   10 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +epsa_w(i,j,k)*dti_fw*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle i decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldydz_u(imax,2:jmax-1,2:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=1 du present rank est vel_u en i=imax-1 du voisin aval
      if(minval(veldydz_u(2,2:jmax-1,2:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldydz_u(3:imax,2:jmax-1,2:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 11 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(ouest)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(est)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call turbulence_adv_epsa_impobc(2,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=imax-1,2,-1

!      epsa_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *epsa_w(i,j,k)                              &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
!       *epsa_w(i+1,j,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*epsa_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       epsa_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         epsa_w(i,j,k) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
        *(epsa_w(i+1,j,k)-epsa_w(i,j,k))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*epsa_w(i,j,k)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   11 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +epsa_w(i,j,k)*dti_fw*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans turbulence_adv_epsa_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur epsa_w(i,j,k)
!     call obc_mpi_epsa_w(1,2,'zb')

!...........................
! Advection axe Oj
!...........................

! Passage boucle j croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldxdz_v(2:imax-1,2,2:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=jmax du present rank est vel_v en j=2 du voisin aval
      if(maxval(veldxdz_v(2:imax-1,jmax,2:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldxdz_v(2:imax-1,2:jmax-1,2:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 20 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(nord)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(sud)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call turbulence_adv_epsa_impobc(3,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

!      epsa_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *epsa_w(i,j,k)                              &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
!       *epsa_w(i,j-1,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*epsa_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       epsa_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         epsa_w(i,j,k) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
        *(epsa_w(i,j-1,k)-epsa_w(i,j,k))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*epsa_w(i,j,k)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->
   20 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +epsa_w(i,j,k)*dti_fw*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle j decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldxdz_v(2:imax-1,jmax,2:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=1 du present rank est vel_v en j=jmax-1 du voisin aval
      if(minval(veldxdz_v(2:imax-1,2,2:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldxdz_v(2:imax-1,3:jmax,2:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->
       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 21 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(sud)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(nord)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call turbulence_adv_epsa_impobc(4,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=jmax-1,2,-1 ; do i=2,imax-1

!      epsa_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *epsa_w(i,j,k)                              &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
!       *epsa_w(i,j+1,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*epsa_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       epsa_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         epsa_w(i,j,k) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
        *(epsa_w(i,j+1,k)-epsa_w(i,j,k))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*epsa_w(i,j,k)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   21 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +epsa_w(i,j,k)*dti_fw*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans turbulence_adv_epsa_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur epsa_w(i,j,k)
!     call obc_mpi_epsa_w(1,2,'zb')

      end subroutine turbulence_adv_epsa_impli

!..............................................................................
!31-12-22
      subroutine turbulence_adv_epsa_impobc(case_,idvar_,send_status_,recv_status_)
      use module_principal
      use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =2
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_,loop_,idvar_,case_,send_status_,recv_status_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_

      if(.not.allocated(adv_send_est)) then !>>>
               allocate(adv_send_est  (2*kmax*jmax)) ; adv_send_est=0.
               allocate(adv_send_ouest(2*kmax*jmax)) ; adv_send_ouest=0.
               allocate(adv_recv_est  (2*kmax*jmax)) ; adv_recv_est=0.
               allocate(adv_recv_ouest(2*kmax*jmax)) ; adv_recv_ouest=0.
      endif                                    !>>>
      if(.not.allocated(adv_send_nord)) then !>>>
               allocate(adv_send_nord  (2*kmax*imax)) ; adv_send_nord=0.
               allocate(adv_send_sud   (2*kmax*imax)) ; adv_send_sud=0.
               allocate(adv_recv_nord  (2*kmax*imax)) ; adv_recv_nord=0.
               allocate(adv_recv_sud   (2*kmax*imax)) ; adv_recv_sud=0.
      endif                                    !>>>


      if(case_==1) then !---echange-Oi-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        i=imax-1
        do k=2,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_est(k0)=epsa_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_est(1)     & ! envoyé au proc Est
                ,size(adv_send_est(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(est)            &
                ,tagest_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*jmax
       k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_ouest(1)    &
               ,size(adv_recv_ouest(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(ouest)                           &
               ,tagest_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=1
       k0=0
       do k=2,kmax ; do j=2,jmax-1
        k0=k0+1
        epsa_w(i,j,k)=adv_recv_ouest(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oi-boucle-croissante--->

      if(case_==2) then !---echange-Oi-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        i=2
        do k=2,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_ouest(k0)=epsa_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_ouest(1)     & ! envoyé au proc Est
                ,size(adv_send_ouest(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(ouest)         &
                ,tagouest_                  &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*jmax
       k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_est(1)     &
               ,size(adv_recv_est(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(est)         &
               ,tagouest_                &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=imax  
       k0=0
       do k=2,kmax ; do j=2,jmax-1
        k0=k0+1
        epsa_w(i,j,k)=adv_recv_est(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oi-boucle-decroissante--->

      if(case_==3) then !---echange-Oj-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        j=jmax-1
        do k=2,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_nord(k0)=epsa_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_nord(1)     & ! envoyé au proc Est
                ,size(adv_send_nord(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(nord)            &
                ,tagnord_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*imax
       k0=  (kmax-1)*(imax-2) ! puisque do k=2,kmax ; do i=2,imax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_sud(1)    &
               ,size(adv_recv_sud(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(sud)                           &
               ,tagnord_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=1
       k0=0
       do k=2,kmax ; do i=2,imax-1
        k0=k0+1
        epsa_w(i,j,k)=adv_recv_sud(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oj-boucle-croissante--->

      if(case_==4) then !---echange-Oj-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        j=2
        do k=2,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_sud(k0)=epsa_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_sud(1)     & ! envoyé au proc Est
                ,size(adv_send_sud(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(sud)         &
                ,tagsud_                    &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*imax
       k0=  (kmax-1)*(imax-2) ! puisque do k=2,kmax ; do i=2,imax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_nord(1)     &
               ,size(adv_recv_nord(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(nord)         &
               ,tagsud_                  &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=jmax
       k0=0
       do k=2,kmax ; do i=2,imax-1
        k0=k0+1
        epsa_w(i,j,k)=adv_recv_nord(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine turbulence_adv_epsa_impobc

!....................................................
!31-12-22
      subroutine turbulence_adv_tken_impli(id_div_)
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down & ! 0= pas d'echange avec voisin aval , 1= echange
                ,id_div_            &
                ,work_status_mpi_min_



      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_div_) se poursuit
!     do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1 
!      anyv3d(i,j,k,id_div_)=0. ! reset cumul temperature*divergence(courant) !26-01-17
!     enddo ; enddo ; enddo

!...........................
! Advection axe Oi
!...........................

! Passage boucle i croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0       ! 0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0       ! 0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0       ! 0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldydz_u(2,2:jmax-1,2:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=imax du present rank est vel_u en i=2 du voisin aval
      if(maxval(veldydz_u(imax,2:jmax-1,2:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldydz_u(2:imax-1,2:jmax-1,2:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 10 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(est)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(ouest)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call turbulence_adv_tken_impobc(1,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

!      tken_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tken_w(i,j,k)                              &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
!       *tken_w(i-1,j,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tken_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tken_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tken_w(i,j,k) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &
        *(tken_w(i-1,j,k)-tken_w(i,j,k))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tken_w(i,j,k)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   10 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tken_w(i,j,k)*dti_fw*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle i decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldydz_u(imax,2:jmax-1,2:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=1 du present rank est vel_u en i=imax-1 du voisin aval
      if(minval(veldydz_u(2,2:jmax-1,2:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldydz_u(3:imax,2:jmax-1,2:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 11 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(ouest)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(est)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call turbulence_adv_tken_impobc(2,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=imax-1,2,-1

!      tken_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tken_w(i,j,k)                              &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
!       *tken_w(i+1,j,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tken_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tken_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tken_w(i,j,k) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &
        *(tken_w(i+1,j,k)-tken_w(i,j,k))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tken_w(i,j,k)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   11 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tken_w(i,j,k)*dti_fw*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans turbulence_adv_tken_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur tken_w(i,j,k)
!     call obc_mpi_tken_w(1,2,'zb')

!...........................
! Advection axe Oj
!...........................

! Passage boucle j croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldxdz_v(2:imax-1,2,2:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=jmax du present rank est vel_v en j=2 du voisin aval
      if(maxval(veldxdz_v(2:imax-1,jmax,2:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldxdz_v(2:imax-1,2:jmax-1,2:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 20 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(nord)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(sud)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call turbulence_adv_tken_impobc(3,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

!      tken_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tken_w(i,j,k)                              &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
!       *tken_w(i,j-1,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tken_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tken_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tken_w(i,j,k) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &
        *(tken_w(i,j-1,k)-tken_w(i,j,k))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tken_w(i,j,k)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->
   20 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tken_w(i,j,k)*dti_fw*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle j decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldxdz_v(2:imax-1,jmax,2:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=1 du present rank est vel_v en j=jmax-1 du voisin aval
      if(minval(veldxdz_v(2:imax-1,2,2:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldxdz_v(2:imax-1,3:jmax,2:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 21 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(sud)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(nord)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call turbulence_adv_tken_impobc(4,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=2,kmax ; do j=jmax-1,2,-1 ; do i=2,imax-1

!      tken_w(i,j,k)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       *tken_w(i,j,k)                              &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
!       *tken_w(i,j+1,k)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1)))                 &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tken_w(i,j,k)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tken_w(i,j,k)=mask_t(i,j,k)*( & !LSM>

         tken_w(i,j,k) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &
        *(tken_w(i,j+1,k)-tken_w(i,j,k))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*0.25*(dz_t(i,j,k-1,0)+dz_t(i,j,k,0)+dz_t(i,j,k-1,1)+dz_t(i,j,k,1))) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_fw &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tken_w(i,j,k)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   21 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       anyv3d(i,j,k,id_div_)= & 
       anyv3d(i,j,k,id_div_)  &
      +tken_w(i,j,k)*dti_fw*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans turbulence_adv_tken_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur tken_w(i,j,k)
!     call obc_mpi_tken_w(1,2,'zb')

      end subroutine turbulence_adv_tken_impli

!..............................................................................
!31-12-22
      subroutine turbulence_adv_tken_impobc(case_,idvar_,send_status_,recv_status_)
      use module_principal
      use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =2
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_,loop_,idvar_,case_,send_status_,recv_status_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_

      if(.not.allocated(adv_send_est)) then !>>>
               allocate(adv_send_est  (2*kmax*jmax)) ; adv_send_est=0.
               allocate(adv_send_ouest(2*kmax*jmax)) ; adv_send_ouest=0.
               allocate(adv_recv_est  (2*kmax*jmax)) ; adv_recv_est=0.
               allocate(adv_recv_ouest(2*kmax*jmax)) ; adv_recv_ouest=0.
      endif                                    !>>>
      if(.not.allocated(adv_send_nord)) then !>>>
               allocate(adv_send_nord  (2*kmax*imax)) ; adv_send_nord=0.
               allocate(adv_send_sud   (2*kmax*imax)) ; adv_send_sud=0.
               allocate(adv_recv_nord  (2*kmax*imax)) ; adv_recv_nord=0.
               allocate(adv_recv_sud   (2*kmax*imax)) ; adv_recv_sud=0.
      endif                                    !>>>


      if(case_==1) then !---echange-Oi-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        i=imax-1
        do k=2,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_est(k0)=tken_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_est(1)     & ! envoyé au proc Est
                ,size(adv_send_est(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(est)            &
                ,tagest_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*jmax
       k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_ouest(1)    &
               ,size(adv_recv_ouest(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(ouest)                           &
               ,tagest_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=1
       k0=0
       do k=2,kmax ; do j=2,jmax-1
        k0=k0+1
        tken_w(i,j,k)=adv_recv_ouest(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oi-boucle-croissante--->

      if(case_==2) then !---echange-Oi-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        i=2
        do k=2,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_ouest(k0)=tken_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_ouest(1)     & ! envoyé au proc Est
                ,size(adv_send_ouest(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(ouest)         &
                ,tagouest_                  &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*jmax
       k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_est(1)     &
               ,size(adv_recv_est(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(est)         &
               ,tagouest_                &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=imax  
       k0=0
       do k=2,kmax ; do j=2,jmax-1
        k0=k0+1
        tken_w(i,j,k)=adv_recv_est(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oi-boucle-decroissante--->

      if(case_==3) then !---echange-Oj-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        j=jmax-1
        do k=2,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_nord(k0)=tken_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_nord(1)     & ! envoyé au proc Est
                ,size(adv_send_nord(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(nord)            &
                ,tagnord_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*imax
       k0=  (kmax-1)*(imax-2) ! puisque do k=2,kmax ; do i=2,imax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_sud(1)    &
               ,size(adv_recv_sud(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(sud)                           &
               ,tagnord_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=1
       k0=0
       do k=2,kmax ; do i=2,imax-1
        k0=k0+1
        tken_w(i,j,k)=adv_recv_sud(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oj-boucle-croissante--->

      if(case_==4) then !---echange-Oj-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        j=2
        do k=2,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_sud(k0)=tken_w(i,j,k)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_sud(1)     & ! envoyé au proc Est
                ,size(adv_send_sud(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(sud)         &
                ,tagsud_                    &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
!      k0=2*kmax*imax
       k0=  (kmax-1)*(imax-2) ! puisque do k=2,kmax ; do i=2,imax-1
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_nord(1)     &
               ,size(adv_recv_nord(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(nord)         &
               ,tagsud_                  &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=jmax
       k0=0
       do k=2,kmax ; do i=2,imax-1
        k0=k0+1
        tken_w(i,j,k)=adv_recv_nord(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine turbulence_adv_tken_impobc

!....................................................
