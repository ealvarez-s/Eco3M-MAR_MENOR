










      module module_atmboundlayer
!______________________________________________________________________
! S model
! release S26  - last update: 19-06-16
!______________________________________________________________________
! version   date    description                                                !
! S.26    13-10-14  mise en service                                            !
!         02-11-14  pas de flux sur zones decouvertes                          !
!         09-12-14  subroutines renommees
!         19-11-15  hcla_ remplace par ablheight
!         24-11-15  ajout du vent dans la couche limite simplifiee
!         04-12-15  appel A checkmpi selon valeur du flag_abl2
!         22-12-15  ne pas passer par checkmpi3 si flag_abl2/=1
!         23-12-15  boucles reduites: 2:imax-1 2:jmax-1 
!         26-02-16  modif test sur flag_1dv
!         01-04-16  taille des boucles dans atmboundlayer_checkmpi_1dv
!         19-06-16  mises A jour suite au travail avec Leo et MNH
!______________________________________________________________________

! Description in:
! https://docs.google.com/document/d/1prt44S0lKX0OlF8YVUMFrFesJghx0D3i-UJ5SqtXW7s/edit
      use module_principal
      use module_parallele
      implicit none
      integer :: loopadv_,loopadvmax_ &
                ,kmax_abl=2 ! Nombre de couches sur la verticale
      integer ::        &
!........................
! Reservation memoire pour anyv3d
                 t0w_=1 &
                ,t2w_=2 &
               ,tbef_=3 &
               ,tnow_=4 &
               ,id_u_=5 &
               ,id_v_=6 &
              ,id_fu_=0 &
              ,id_fv_=7 &
!........................
              ,ibeg_    &
              ,istp_    &
              ,jbeg_    &
              ,jstp_  
      double precision :: w0_,w2_,w0bef_,w2bef_

contains

      subroutine atmboundlayer_initial
      implicit none
      end subroutine atmboundlayer_initial
!......................................................................
      subroutine atmboundlayer_sstmeteo(case_,t_)
      implicit none
      integer case_,t_ 

! A ce stade tem_t(:,:,:,2) est vacant (tem_t(1)=tem_t(2))
! On en profite pour charger tem_t(2) avec teta0_t la SST du modele meteo

      x2=(elapsedtime_now           -airseafile_prvtime(t2m_id))      &
        /(airseafile_nextime(t2m_id)-airseafile_prvtime(t2m_id))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x1=1.-x2
      do j=1,jmax
      do i=1,imax
! Bidouille transitoire pour test: on prend la sst de mercator....
!      tem_t(i,j,kmax,2)=timeweightobc(trc_id) *temobc_t(i,j,kmax,2) &
!                   +(1.-timeweightobc(trc_id))*temobc_t(i,j,kmax,0)
       tem_t(i,j,kmax,2)=x1*teta0_t(i,j,1)+x2*teta0_t(i,j,2)-273.15
      enddo
      enddo 

! Comptutes slhf and sshf with the SST of the meteorological model:
      call bulk_formulae(case_,2,2,0.)
! arg1: pahe initiale ou iterative
! arg2: dernier indice de tem_t
! arg3: dernier indice des variables meteo (pss_w(:,:,2) etc...)
! arg4: relativewind

! Retablir tem_t(2) dans sa valeur initiale
      do j=1,jmax
      do i=1,imax
       tem_t(i,j,kmax,2)=tem_t(i,j,kmax,1)
      enddo
      enddo

! La subroutine bulk_formulae renvoie slhf_w(:,:,1),sshf_w(:,:,1)
! On bascule sur slhf_w(:,:,2),sshf_w(:,:,2). L'interpolation temporelle 
! de slhf_w(i,j,0) et slhf_w(i,j,2) donnera le flux "dans ecmwf" a
! l'instant present pour balance avec slhf_w(:,:,1) qui, lui, sera
! calcule avec la vraie sst de S26.

! Archiver les flux bulk du modele meteo avec sa propre sst dans
! slhf(0) slhf(2) sshf(0) sshf(2) sachant que les formules bulk
! avec la SST de S26 donnent directement slhf(1) et sshf(1)
      do j=1,jmax
      do i=1,imax

       slhf_w(i,j,0)=slhf_w(i,j,2) ! t_=0 ou 2
       sshf_w(i,j,0)=sshf_w(i,j,2) ! t_=0 ou 2

       slhf_w(i,j,2)=slhf_w(i,j,1) ! t_=0 ou 2
       sshf_w(i,j,2)=sshf_w(i,j,1) ! t_=0 ou 2

! Pour pouvoir calculer un ecart du niveau de turbulence dans
! la couche de surface (10% de la couche limite) on memorise
! ustar calculE A partir des parametres du modele meteo (sa SST
! en particulier):
!      ustar_t(i,j,0)=ustar_t(i,j,2)
!      ustar_t(i,j,2)=ustar_t(i,j,1)

       kz_abl_w(i,j,0)=kz_abl_w(i,j,2)
       kz_abl_w(i,j,2)=kz_abl_w(i,j,1)

      enddo
      enddo 

      if(flag_abl2==1) then !ABL2ABL2ABL2>
       do j=1,jmax ; do i=1,imax+1
        wstress_u(i,j,0)=wstress_u(i,j,2)
        wstress_u(i,j,2)=wstress_u(i,j,1)
       enddo       ; enddo
       do j=1,jmax+1 ; do i=1,imax
        wstress_v(i,j,0)=wstress_v(i,j,2)
        wstress_v(i,j,2)=wstress_v(i,j,1)
       enddo       ; enddo
      endif                 !ABL2ABL2ABL2>

      end subroutine atmboundlayer_sstmeteo
!......................................................................
      subroutine atmboundlayer_windcfl
      implicit none



! Attention reset de x1=dti_fw deplace AVANT le debut de la boucle verticale
      x1=dti_fw
      do k=1,kmax_abl ! Boucle verticale 

!....................................
! EPAISSEUR DE LA COUCHE k:
      if(k< kmax_abl) then !pmxpmx>
        do j=0,jmax+1 ; do i=0,imax+1
          xy_t(i,j,1)=zw_abl(k+1)-zw_abl(k)
        enddo         ; enddo
      else                 !pmxpmx>
! La derniere couche a une epaisseur variant comme celle de la ABL
        do j=0,jmax+1 ; do i=0,imax+1
          xy_t(i,j,1)=ablheight_t(i,j,2)-zw_abl(k)
        enddo         ; enddo
      endif                !pmxpmx>

!....................................
! TRANSPORTS:
      do j=1,jmax ; do i=1,imax+1
! Oi Flux de masse:
        xy_u(i,j,1)=0.5*( & !oooo>
                xy_t(i  ,j,1)*uwindabl_t(i  ,j,1,2)*dy_t(i  ,j)  &
               +xy_t(i-1,j,1)*uwindabl_t(i-1,j,1,2)*dy_t(i-1,j)  &
                        )   !oooo>
      enddo ; enddo
      do j=1,jmax+1 ; do i=1,imax
! Oj Flux de masse:
        xy_v(i,j,1)=0.5*( & !oooo>
                xy_t(i,j  ,1)*vwindabl_t(i,j  ,1,2)*dx_t(i  ,j)  &
               +xy_t(i,j-1,1)*vwindabl_t(i,j-1,1,2)*dx_t(i,j-1)  &
                        )   !oooo>
      enddo ; enddo

! Attention reset de x1=dti_fw deplace AVANT le debut de la boucle verticale
      do j=1,jmax ; do i=1,imax

       x1=min(x1, & !****>

                -(2.*dxdy_t(i,j)*xy_t(i,j,1))              &
                /( &  !ooooo>
                  -1. & ! cste artificielle pour eviter singularite division par 0
                  -xy_u(i+1,j  ,1)-abs(xy_u(i+1,j  ,1))    &
                  +xy_u(i  ,j  ,1)-abs(xy_u(i  ,j  ,1))    &
                  -xy_v(i  ,j+1,1)-abs(xy_v(i  ,j+1,1))    &
                  +xy_v(i  ,j  ,1)-abs(xy_v(i  ,j  ,1))    &
                 ) &  !ooooo>
              )     !****>

      enddo ; enddo

      enddo ! boucle verticale k

      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_min,par%comm2d,ierr)
      dt_abl_max(0)=dt_abl_max(2)
      dt_abl_max(2)=x2

      if(par%rank==0)then !pmxpmx>
       open(unit=3,file='tmp/abl_dtimax',position='append')
        write(3,*)real(elapsedtime_now/86400.),real(dt_abl_max(2)),real(dti_fw)
       close(3)
      endif               !pmxpmx>

!#ifdef 1
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!       stop 'coucou2'

      end subroutine atmboundlayer_windcfl

!......................................................................

      subroutine atmboundlayer_balance
      implicit none
!     real :: hcla_=1000. ! hauteur de la CLA en metres
!     double precision :: w0_,w2_

!     ibeg_=1 ;  istp_=imax   ; jbeg_=1 ; jstp_=jmax    ! option grande boucle
      ibeg_=2 ;  istp_=imax-1 ; jbeg_=2 ; jstp_=jmax-1  ! option petite boucle

! Time interpolation of the 2 bounding best ABL time steps:
! loopadvmax_ = Number of advection sub-iterations = (Ocean model DT ) / ( ABL model best DT )
      x2=min(max(  &
       (elapsedtime_now            -airseafile_prvtime(u10m_id))      &
      /(airseafile_nextime(u10m_id)-airseafile_prvtime(u10m_id)),0.d00),1.d00) 


!      write(6,*)dt_abl_max(0) &
!               ,(1.-x2)*dt_abl_max(0)+x2*dt_abl_max(2) &
!               ,dt_abl_max(2)


       loopadvmax_= 1+int( dti_fw/( (1.-x2)*dt_abl_max(0)+x2*dt_abl_max(2) ) )
    
! dt_abl_max(1)= best ABL time step
       dt_abl_max(1)=dti_fw/loopadvmax_

! ADVECTION DE LA CLA:
      w2_=(elapsedtime_now            -airseafile_prvtime(u10m_id))   &
         /(airseafile_nextime(u10m_id)-airseafile_prvtime(u10m_id))
      w2_=min(max(w2_,0.d00),1.d00) 
      w0_=1.-w2_

! Idem mais au temps precedent:
      w2bef_=(elapsedtime_now-dti_fw     -airseafile_prvtime(u10m_id)) &
            /(airseafile_nextime(u10m_id)-airseafile_prvtime(u10m_id))
      w2bef_=min(max(w2bef_,0.d00),1.d00) 
      w0bef_=1.-w2bef_

      call atmboundlayer_massflux ! Transport anyv3d(:,:,:,id_u_)
                                  ! transport anyv3d(:,:,:,id_v_)
                                  ! Vitesse verticale wwindabl_w
                                  ! Epaisseur en 4 temps ! anyv3d(:,:,:,temps)
                                  ! avec temps: t0w_ t2m_ tnow_ tbef_


! Les echanges se font:
! 1. a l'initialisation (iteration3d==0) pour garantir la continuite mpi du champ initial
! 2. a la fin de chaque sous-iteration advective
! 3. apres la mise a jour par le delta de flux
      if(iteration3d==0) then !pmxpxm>
        call atmboundlayer_init_debug   ! Etat initial: verifier la memoire
        call atmboundlayer_mpi          ! mpi: 1: initialisation
      endif                   !pmxpxm>

! Pour tester les proprietes de conservation et pour cela
! voir au la C.L. au sommet de la ABL dans subroutine matrix_a
!     teta2delta_t(:,:,1)=1.
!        q2delta_t(:,:,1)=1.
!     teta2delta_t(:,:,2)=0.
!        q2delta_t(:,:,2)=0.
!     uwinddelta_t(:,:,1)=1.
!     uwinddelta_t(:,:,2)=0.
!     vwinddelta_t(:,:,1)=1.
!     vwinddelta_t(:,:,2)=0.
!     sshf_w=0.
!     slhf_w=0.
!     wstress_u=0.
!     wstress_v=0.

      do loopadv_=1,loopadvmax_ !-adv-adv-adv->

       rap2=real(loopadv_  )/real(loopadvmax_)
       rap1=real(loopadv_-1)/real(loopadvmax_)
       x1=(1.-rap1)
       x2=    rap1 
       x3=(1.-rap2)
       x4=    rap2 

! Champ teta2delta ...............................................
! Coef de la matrice principal ...................................
       call atmboundlayer_matrix_a('t')                           ! Advection verticale implicite
                                                                  ! Melange turbulent
       call atmboundlayer_flux_dteta                              ! Flux advectifs horizontaux
       call atmboundlayer_matrix_dteta                            ! Membre de droite de la matrice (termes explicites)
                                                                  ! Div. flux horizontaux
       call tridiagonalsolver(1,ibeg_,istp_,jbeg_,jstp_,kmax_abl) ! Solveur, arg1=1 car premiere variable
       call atmboundlayer_update_dteta                            ! Mise A jour avec la solution du solver

! Champ q2delta ..................................................
! Coef de la matrice principal ...................................
       call atmboundlayer_matrix_a('q')                           ! Advection verticale implicite
                                                                  ! Melange turbulent
       call atmboundlayer_flux_dq                                 ! Flux advectifs horizontaux
       call atmboundlayer_matrix_dq                               ! Membre de droite de la matrice (termes explicites)
                                                                  ! Div. flux horizontaux
       call tridiagonalsolver(1,ibeg_,istp_,jbeg_,jstp_,kmax_abl) ! Solveur, arg1=0 car variable suivante
       call atmboundlayer_update_dq                               ! Mise A jour avec la solution du solver

! Si option SABL2 alors ajouter les equations pour le vent:
       if(flag_abl2==1) then ! ABL2ABL2ABL2 > 

! Champ uwinddelta ...............................................
       call atmboundlayer_flux_du                                 ! Flux advectifs horizontaux
       call atmboundlayer_matrix_du                               ! Membre de droite de la matrice (termes explicites)
                                                                  ! Div. flux horizontaux
       call tridiagonalsolver(0,ibeg_,istp_,jbeg_,jstp_,kmax_abl) ! Solveur, arg1=0 car variable suivante
       call atmboundlayer_update_du                               ! Mise A jour avec la solution du solver

! Champ vwinddelta ...............................................
       call atmboundlayer_flux_dv                                 ! Flux advectifs horizontaux
       call atmboundlayer_matrix_dv                               ! Membre de droite de la matrice (termes explicites)
                                                                  ! Div. flux horizontaux
       call tridiagonalsolver(0,ibeg_,istp_,jbeg_,jstp_,kmax_abl) ! Solveur, arg1=0 car variable suivante
       call atmboundlayer_update_dv                               ! Mise A jour avec la solution du solver

       endif                 ! ABL2ABL2ABL2 > 

       call atmboundlayer_mpi ! mpi 2: fin de chaque sous-etape iteratie
      enddo                     !-adv-adv-adv->


!      if(par%rank==0)write(6,*)'ATTENTION GRAPHIQUES'
!      call graph_out
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'bibi2'

! Cas 1DV:
      if(flag_1dv==1)call atmboundlayer_1dv !27-11-15


      end subroutine atmboundlayer_balance

!......................................................................

      subroutine atmboundlayer_fieldcorrection
      implicit none

! On note que teta2_t(i,j,1) est toujours remis a jour par l'interpolation
! des echeances 0 et 2. Ces lignes ne viennent qu'ajouter la perturbation
! au champs de base avant d'entrer dans le calcul des formules bulk

! Update teta2
      do j=1,jmax ; do i=1,imax
      teta2_t(i,j,1)=teta2_t(i,j,1)+teta2delta_t(i,j,1)*mask_t(i,j,kmax)
      enddo       ; enddo

!     if(par%rank==0)write(6,*)'Q2 NON CORRIGE'
!#ifdef bidon
! Update q2
      do j=1,jmax ; do i=1,imax
! Par prudence, la condition n'est pas appliquee en zone continentale:
      if(mask_t(i,j,kmax)==1) then !1111111>

! Avant d'ajouter dq au champ q, appliquer des bornes garantissant 0<q<qsat
! Condition pour q>0:
       q2delta_t(i,j,1)=max(q2delta_t(i,j,1),-q2_t(i,j,1))
            q2_t(i,j,1)=q2_t(i,j,1)+q2delta_t(i,j,1)

!#ifdef bidon
! Condition pour q<qsat
! specific humidity at saturation at the atm. level 
         x1=exp(xalpw-xbetaw/teta2_t(i,j,1)-xgamw*log(teta2_t(i,j,1)))/(pss_w(i,j,1)-rhoair*grav*2.)
         qsat_atm_z=(xrd/xrv)*x1/(1.+((xrd/xrv)-1.)*x1)
         if(q2_t(i,j,1)>qsat_atm_z) then !phasechange>
          q2delta_t(i,j,1)=q2delta_t(i,j,1)+qsat_atm_z-q2_t(i,j,1)
! Prevoir un terme d'apport de chaleur A l'equation pour teta2delta
! il sera a ajouter aux termes sshf (signe A determiner en fonction convention de la SABL):
! phasechange_t(i,j)=(q2_t(i,j,1)-qsat_atm_z)/(rhoair*lv)/dt_abl_max(1)
! Puis borner q2:
          q2_t(i,j,1)=qsat_atm_z
!         write(10+par%rank,*)q2_t(i,j,1),qsat_atm_z
         endif                           !phasechange>
!     write(10+par%rank,*)qsat_atm_z*(pss_w(i,j,1)-rhoair*grav*2.)/0.622 &
!     ,teta2_t(i,j,1)-273.15 ! pour validation tableau livre Queney
!#endif

      endif                        !1111111>
      enddo       ; enddo
!#endif

!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'victory'


      if(flag_abl2==1) then !abl2abl2>
       do j=1,jmax   ; do i=1,imax  !01-04-16
        uwind_t(i,j,1)=uwind_t(i,j,1)+uwinddelta_t(i,j,1)*mask_t(i,j,kmax)
        vwind_t(i,j,1)=vwind_t(i,j,1)+vwinddelta_t(i,j,1)*mask_t(i,j,kmax)
       enddo       ; enddo
      endif                 !abl2abl2>

      end subroutine atmboundlayer_fieldcorrection

!......................................................................

      subroutine atmboundlayer_mpi
      implicit none
      integer loop_,id_teta2delta_,id_q2delta_ &
                   ,id_uwinddelta_,id_vwinddelta_

!     call get_type_echange('z0','teta2delta_t_z1',teta2delta_t,lbound(teta2delta_t),ubound(teta2delta_t),1,id_teta2delta_)
!     call get_type_echange('z0','q2delta_t_z1'   ,q2delta_t   ,lbound(q2delta_t)   ,ubound(q2delta_t)   ,1,id_q2delta_)
!     if(flag_abl2==1) then !>>>
!     call get_type_echange('z0','uwinddelta_t_z1',uwinddelta_t,lbound(uwinddelta_t),ubound(uwinddelta_t),1,id_uwinddelta_)
!     call get_type_echange('z0','vwinddelta_t_z1',vwinddelta_t,lbound(vwinddelta_t),ubound(vwinddelta_t),1,id_vwinddelta_)
!     endif                 !>>>
      call get_type_echange('z0','teta2delta_t_z1',teta2delta_t,lbound(teta2delta_t),ubound(teta2delta_t),id_teta2delta_)
      call get_type_echange('z0','q2delta_t_z1'   ,q2delta_t   ,lbound(q2delta_t)   ,ubound(q2delta_t)   ,id_q2delta_)
      if(flag_abl2==1) then !>>>
      call get_type_echange('z0','uwinddelta_t_z1',uwinddelta_t,lbound(uwinddelta_t),ubound(uwinddelta_t),id_uwinddelta_)
      call get_type_echange('z0','vwinddelta_t_z1',vwinddelta_t,lbound(vwinddelta_t),ubound(vwinddelta_t),id_vwinddelta_)
      endif                 !>>>

      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(teta2delta_t,id_teta2delta_,mpi_neighbor_list(loop_))
        call echange_voisin(q2delta_t   ,id_q2delta_   ,mpi_neighbor_list(loop_))
        if(flag_abl2==1) then !>>>
        call echange_voisin(uwinddelta_t,id_uwinddelta_,mpi_neighbor_list(loop_))
        call echange_voisin(vwinddelta_t,id_vwinddelta_,mpi_neighbor_list(loop_))
        endif                 !>>>
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine atmboundlayer_mpi
!..............................................................................

!................................................................................................
!................................................................................................

!................................................................................................

!................................................................................................

!................................................................................................

      subroutine atmboundlayer_1dv !27-11-15
      implicit none

      if(flag_abl==1) then !1111> !27-11-15
         teta2delta_t(:,:,1)=teta2delta_t(2,2,1)
            q2delta_t(:,:,1)=   q2delta_t(2,2,1)
         if(flag_abl2==1) then !2222>
           uwinddelta_t(:,:,1)=uwinddelta_t(2,2,1) 
           vwinddelta_t(:,:,1)=vwinddelta_t(2,2,1) 
         endif                 !2222>
      endif                !1111>

      end subroutine atmboundlayer_1dv

!................................................................................................

      subroutine atmboundlayer_massflux
      implicit none


! Boucle verticale
      do k=1,kmax_abl !VVVVVVVVVV>

!....................................
! EPAISSEUR DE LA COUCHE k aux temps des fichiers
      if(k< kmax_abl) then !pmxpmx>
        do j=0,jmax+1 ; do i=0,imax+1
          anyv3d(i,j,k,t0w_)=zw_abl(k+1)-zw_abl(k)
          anyv3d(i,j,k,t2w_)=zw_abl(k+1)-zw_abl(k)
        enddo         ; enddo
      else                 !pmxpmx>
! La derniere couche a une epaisseur variant comme celle de la ABL
        do j=0,jmax+1 ; do i=0,imax+1
          anyv3d(i,j,k,t0w_)=ablheight_t(i,j,t0w_)-zw_abl(k)
          anyv3d(i,j,k,t2w_)=ablheight_t(i,j,t2w_)-zw_abl(k)
        enddo         ; enddo
      endif                !pmxpmx>
!....................................
! EPAISSEUR DE LA COUCHE k aux temps des fichiers
!     if(k< kmax_abl) then !pmxpmx>
!       do j=0,jmax+1 ; do i=0,imax+1
!         anyv3d(i,j,k,t0w_)=0.1*ablheight_t(i,j,t0w_)
!         anyv3d(i,j,k,t2w_)=0.1*ablheight_t(i,j,t2w_)
!       enddo         ; enddo
!     else                 !pmxpmx>
! La derniere couche a une epaisseur variant comme celle de la ABL
!       do j=0,jmax+1 ; do i=0,imax+1
!         anyv3d(i,j,k,t0w_)=0.9*ablheight_t(i,j,t0w_)
!         anyv3d(i,j,k,t2w_)=0.9*ablheight_t(i,j,t2w_)
!       enddo         ; enddo
!     endif                !pmxpmx>

!....................................
! TRANSPORTS:
      do j=jbeg_,jstp_ ; do i=ibeg_,istp_+1
! Oi Flux de masse:
        anyv3d(i,j,k,id_u_)=0.5*( & !oooo>

         w0_*( anyv3d(i  ,j,k,t0w_)*uwindabl_t(i  ,j,1,t0w_)*dy_t(i  ,j)  &
              +anyv3d(i-1,j,k,t0w_)*uwindabl_t(i-1,j,1,t0w_)*dy_t(i-1,j)) &

        +w2_*( anyv3d(i  ,j,k,t2w_)*uwindabl_t(i  ,j,1,t2w_)*dy_t(i  ,j)  &
              +anyv3d(i-1,j,k,t2w_)*uwindabl_t(i-1,j,1,t2w_)*dy_t(i-1,j)) &

                                )   !oooo>
      enddo ; enddo

      do j=jbeg_,jstp_+1 ; do i=ibeg_,istp_
! Oj Flux de masse:
        anyv3d(i,j,k,id_v_)=0.5*( & !oooo>

         w0_*( anyv3d(i,j  ,k,t0w_)*vwindabl_t(i,j  ,1,t0w_)*dx_t(i  ,j)  &
              +anyv3d(i,j-1,k,t0w_)*vwindabl_t(i,j-1,1,t0w_)*dx_t(i,j-1)) &

        +w2_*( anyv3d(i,j  ,k,t2w_)*vwindabl_t(i,j  ,1,t2w_)*dx_t(i  ,j)  &
              +anyv3d(i,j-1,k,t2w_)*vwindabl_t(i,j-1,1,t2w_)*dx_t(i,j-1)) &

                        )   !oooo>
      enddo ; enddo

!...........................................................
! EPAISSEUR DE LA COUCHE k aux temps now et before du calcul
      do j=jbeg_,jstp_ ; do i=ibeg_,istp_
        anyv3d(i,j,k,tbef_)=w0bef_*anyv3d(i,j,k,t0w_)+w2bef_*anyv3d(i,j,k,t2w_)
        anyv3d(i,j,k,tnow_)=w0_   *anyv3d(i,j,k,t0w_)+w2_   *anyv3d(i,j,k,t2w_)
      enddo ; enddo

! Vitesse verticale au sommet de la CLA:
       do j=jbeg_,jstp_ ; do i=ibeg_,istp_
        wwindabl_w(i,j,k+1)=                                          &
        wwindabl_w(i,j,k  )                                           &
          -( anyv3d(i+1,j  ,k,id_u_)-anyv3d(i,j,k,id_u_)              &! -dU/dx
            +anyv3d(i  ,j+1,k,id_v_)-anyv3d(i,j,k,id_v_))/dxdy_t(i,j) &! -dV/dy
           -(anyv3d(i  ,j  ,k,tnow_)-anyv3d(i,j,k,tbef_))/dti_fw       ! -ddz/dt   
       enddo ; enddo

      enddo           !VVVVVVVVVV> ! Fin de boucle verticale

 
!      do k=1,kmax_abl
!      do j=jbeg_,jstp_ ; do i=ibeg_,istp_
!       write(10+par%rank,*) &
!       (anyv3d(i+1,j,k,id_u_)-anyv3d(i,j,k,id_u_)              &
!       +anyv3d(i,j+1,k,id_v_)-anyv3d(i,j,k,id_v_))/dxdy_t(i,j) &
!      ,wwindabl_w(i,j,k+1)-wwindabl_w(i,j,k)                     &
!      ,(anyv3d(i  ,j  ,k,tnow_)-anyv3d(i,j,k,tbef_))/dti_fw 
!      enddo ; enddo
!      enddo

!     call graph_out
!     stop 'coucou4'
      end subroutine atmboundlayer_massflux

!..............................................................................

      subroutine atmboundlayer_matrix_a(txt_)
      implicit none
      character*1 txt_

!.................................
! Pas de temps fois Coef melange vertical turbulent SUR epaisseur de la couche

! CAS TEMPERATURE
      if(txt_=='t') then !>>>>>>
       do k=2,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_
        xy_t(i,j,k)=dt_abl_max(1)                                 & !  dt
           /anyv3d(i,j,k-1,tnow_)                                 & 
!          /(0.5*(anyv3d(i,j,k-1,tnow_)+anyv3d(i,j,k,tnow_)))     & ! /dz risque = dz trop grand
!              *1.    & ! KZ trop grand
!              *0.5   & ! KZ trop grand
!              *0.2   & ! KZ trop grand
!              *0.1   & ! KZ trop grand
!              *0.01  & ! KZ un peu trop grand?
               *0.005 & ! ??
               *kz_abl_w(i,j,1)  
       enddo ; enddo ; enddo
      endif              !>>>>>>

! CAS HUMIDITE
      if(txt_=='q') then !>>>>>>
       do k=2,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_
        xy_t(i,j,k)=dt_abl_max(1)                                 & !  dt
           /anyv3d(i,j,k-1,tnow_)                                 & 
!              *0.2  & ! KZ trop petit
!              *0.5  & ! KZ un peu trop grand?
               *0.05 & ! ??
               *kz_abl_w(i,j,1)  
       enddo ; enddo ; enddo
      endif              !>>>>>>


! Flux aux limites
       do j=jbeg_,jstp_ ; do i=ibeg_,istp_
        xy_t(i,j,1         )=0.  ! Pas de flux A la Base de la ABL
!       xy_t(i,j,kmax_abl+1)=0.  ! Pas de flux au Sommet de la ABL
        xy_t(i,j,kmax_abl+1)=0.1*xy_t(i,j,2) ! flux au sommet
       enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_

!..........................
        tridia_in(i,j,k,3)=   &
! melange turbulent:
        -xy_t(i,j,k+1)        & 
! Advection verticale:
        +0.5*dt_abl_max(1)*(wwindabl_w(i,j,k+1)-abs(wwindabl_w(i,j,k+1)))

!..........................
        tridia_in(i,j,k,1)=   &
! melange turbulent:
        -xy_t(i,j,k)          & 
! Advection verticale:
        -0.5*dt_abl_max(1)*(wwindabl_w(i,j,k  )+abs(wwindabl_w(i,j,k  )))

!..........................
        tridia_in(i,j,k,2)=x4*anyv3d(i,j,k,tnow_)+x3*anyv3d(i,j,k,tbef_) & ! dz(t+1)
! melange turbulent:
        +xy_t(i,j,k+1)+xy_t(i,j,k)    &
! Advection verticale:
        +0.5*dt_abl_max(1)*( & !ooo>
                         wwindabl_w(i,j,k+1)+abs(wwindabl_w(i,j,k+1)) &
                        -wwindabl_w(i,j,k  )+abs(wwindabl_w(i,j,k  )) &
                           )   !ooo>

       enddo ; enddo ; enddo

! C.L. 
       do j=jbeg_,jstp_ ; do i=ibeg_,istp_
! Si gradient nul au sommet de la ABL:
!       tridia_in(i,j,kmax_abl,2)=tridia_in(i,j,kmax_abl,2)+tridia_in(i,j,kmax_abl,3)
!       tridia_in(i,j,kmax_abl,3)=0.
! Si Valeur nulle au dessus de la ABL:
        tridia_in(i,j,kmax_abl,3)=0.
       enddo            ; enddo 


      end subroutine atmboundlayer_matrix_a

!..............................................................................

      subroutine atmboundlayer_matrix_dteta
      implicit none



! Ajouter la contribution du flux au premier niveau de la variable
! avant de calculer tridia_in(:,:,:,4)
! A noter qu'en remplacant dt_abl_max(1) par dti_fw on pourrait ne
! faire cette operation qu'une fois (et non pas A tous les sous-pas de
! temps). Prevoir de faire un test.....
         k=1
         do j=jbeg_,jstp_ ; do i=ibeg_,istp_

          teta2delta_t(i,j,1)=                             &
          teta2delta_t(i,j,1)                              &

       -dt_abl_max(1)/(rhoair*cp_air*anyv3d(i,j,k,tnow_))  & !18-11-15
                     *(                                    &

                     !                    (slhf_w(i,j,1)   &
                     ! -w0_*(              slhf_w(i,j,0))  &
                     ! -w2_*(              slhf_w(i,j,2))) &

                                          +sshf_w(i,j,1)   &
                       -w0_*(              sshf_w(i,j,0))  &
                       -w2_*(              sshf_w(i,j,2))  &
                                                         ) &
                           *mask_t(i,j,kmax)               &
                        *wetmask_t(i,j)                     !02-11-14

          enddo ; enddo

! Calcul de la matrice RHS:

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_

         tridia_in(i,j,k,4)=                            & 

      teta2delta_t(i,j,k)                               &
       *(x2*anyv3d(i,j,k,tnow_)+x1*anyv3d(i,j,k,tbef_)) & !teta(t)*dz(t)

         -dt_abl_max(1)*( anyv3d(i+1,j,k,id_fu_)-anyv3d(i,j,k,id_fu_)  &
                         +anyv3d(i,j+1,k,id_fv_)-anyv3d(i,j,k,id_fv_))/dxdy_t(i,j)  

       enddo ; enddo ; enddo

      end subroutine atmboundlayer_matrix_dteta

!..............................................................................

      subroutine atmboundlayer_update_dteta
      implicit none

      do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_
       teta2delta_t(i,j,k)=tridia_out(i,j,k)  
      enddo ; enddo ; enddo

      end subroutine atmboundlayer_update_dteta

!..............................................................................

      subroutine atmboundlayer_flux_dteta
      implicit none

       x12=1./12.

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_+1

! Oi Flux de teta2delta:
        anyv3d(i,j,k,id_fu_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_u_)*( &       !11111>
          6.*(teta2delta_t(i  ,j,k)+teta2delta_t(i-1,j,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_u_))*( & !22222>
          6.*(teta2delta_t(i  ,j,k)-teta2delta_t(i-1,j,k))  & !up2
                                  ) & !22222>

            )     !ooo>

       enddo ; enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_+1 ; do i=ibeg_,istp_

! Oj Flux de teta2delta:
        anyv3d(i,j,k,id_fv_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_v_)*( &       !11111>
          6.*(teta2delta_t(i,j  ,k)+teta2delta_t(i,j-1,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_v_))*( & !22222>
          6.*(teta2delta_t(i,j  ,k)-teta2delta_t(i,j-1,k))  & !up2
                                  ) & !22222>

            )     !ooo>


       enddo ; enddo ; enddo

      end subroutine atmboundlayer_flux_dteta

!..............................................................................

      subroutine atmboundlayer_matrix_dq
      implicit none

! Ajouter la contribution du flux au premier niveau de la variable
! avant de calculer tridia_in(:,:,:,4)
! A noter qu'en remplacant dt_abl_max(1) par dti_fw on pourrait ne
! faire cette operation qu'une fois (et non pas A tous les sous-pas de
! temps). Prevoir de faire un test.....
         k=1
         do j=jbeg_,jstp_ ; do i=ibeg_,istp_

          q2delta_t(i,j,1)=                                &
          q2delta_t(i,j,1)                                 &

         -dt_abl_max(1)/(rhoair*lv*anyv3d(i,j,k,tnow_))    & !18-11-15
                       *(                                  &

                             +slhf_w(i,j,1)                &
                         -w0_*slhf_w(i,j,0)                & 
                         -w2_*slhf_w(i,j,2)                &

!      +(w0_*ablheight_t(i,j,t0w_)+w2_*ablheight_t(i,j,t2w_)) &
!     *1e-3*teta2delta_t(i,j,1)                               &

                                           )               &
                             *mask_t(i,j,kmax)             &
                            *wetmask_t(i,j)                     !02-11-14

          enddo ; enddo

! Calcul de la matrice RHS

      do k=1,kmax_abl 
      if(k==1)k0=1
      if(k==2)k0=-1
      if(k>2 )k0=0
!     k0=0
      do j=jbeg_,jstp_ ; do i=ibeg_,istp_

         tridia_in(i,j,k,4)=                            & 

         q2delta_t(i,j,k)                               &
       *(x2*anyv3d(i,j,k,tnow_)+x1*anyv3d(i,j,k,tbef_)) & !teta(t)*dz(t)

         -dt_abl_max(1)*( & !ooo>

           ( anyv3d(i+1,j,k,id_fu_)-anyv3d(i,j,k,id_fu_)              &
            +anyv3d(i,j+1,k,id_fv_)-anyv3d(i,j,k,id_fv_))/dxdy_t(i,j) & 

! DF/dz avec F=Delta_Kz*dq/dz avec Delta_Kz=Delta_ustar*karman*z avec
! z hauteur entre couches 1 et 2 donc anyv3d(i,j,1,tnow_), 
! dq/dz=(0-q2m)/H_couche_limite
! Quant au dz de dF/dz il disparait car l'equation est homogene A q*dz
        +k0*(     kz_abl_w(i,j,1)       &    ! Delta.Kz*Dq/Dz entre couches 1 et 2
             -w0_*kz_abl_w(i,j,0)       &
             -w2_*kz_abl_w(i,j,2) )     &
           *( w0_*q2_t(i,j,0)           &
             +w2_*q2_t(i,j,2) )*0.5d-4 &
!                    *q2_t(i,j,2)/ablheight_t(i,j,1)  &
!                    *1.e-6  &

                        )   !ooo>

!!     if(k==1.and.i+par%timax(1)==459.and.j+par%tjmax(1)==560) &
!      write(66,*) &
!       +k0*(     kz_abl_w(i,j,1)     &   ! Delta.Kz*Dq/Dz entre couches 1 et 2
!            -w0_*kz_abl_w(i,j,0)     &
!            -w2_*kz_abl_w(i,j,2) )   &
!          *( w0_*q2_t(i,j,0)         &
!            +w2_*q2_t(i,j,2) )*1.d-4 &
!       ,1./(rhoair*lv)    & !18-11-15
!                      *(                                  &
!                            +slhf_w(i,j,1)                &
!                        -w0_*slhf_w(i,j,0)                & 
!                        -w2_*slhf_w(i,j,2) )

      enddo ; enddo
      enddo

      end subroutine atmboundlayer_matrix_dq

!..............................................................................

      subroutine atmboundlayer_update_dq
      implicit none

      do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_
       q2delta_t(i,j,k)=tridia_out(i,j,k)  
!      q2delta_t(i,j,k)=0.
! Bidouille melanger dq sur la couche limite:
!      x1=w0_*anyv3d(i,j,1,t0w_)+w2_*anyv3d(i,j,1,t2w_)
!      x2=w0_*anyv3d(i,j,2,t0w_)+w2_*anyv3d(i,j,2,t2w_)
!      q2delta_t(i,j,k)=(x1*tridia_out(i,j,1)+tridia_out(i,j,2))/(x1+x2)
      enddo ; enddo ; enddo

! Au sommet de la couche rappel vers zEro avec une echelle de temps de
! 10 jours (Deremble et al 2013 Eq11)
!     kmax_abl ; x0=1./(1.+dt_abl_max(1)/(10.*86400.))
!     do j=jbeg_,jstp_ ; do i=ibeg_,istp_
!      q2delta_t(i,j,k)=q2delta_t(i,j,k)*x0
!     enddo ; enddo ; enddo
  


      end subroutine atmboundlayer_update_dq

!..............................................................................

      subroutine atmboundlayer_flux_dq
      implicit none

       x12=1./12.

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_+1

! Oi Flux de q2delta:
        anyv3d(i,j,k,id_fu_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_u_)*( &       !11111>
          6.*(q2delta_t(i  ,j,k)+q2delta_t(i-1,j,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_u_))*( & !22222>
          6.*(q2delta_t(i  ,j,k)-q2delta_t(i-1,j,k))  & !up2
                                  ) & !22222>

            )     !ooo>

       enddo ; enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_+1 ; do i=ibeg_,istp_

! Oj Flux de q2delta:
        anyv3d(i,j,k,id_fv_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_v_)*( &       !11111>
          6.*(q2delta_t(i,j  ,k)+q2delta_t(i,j-1,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_v_))*( & !22222>
          6.*(q2delta_t(i,j  ,k)-q2delta_t(i,j-1,k))  & !up2
                                  ) & !22222>

            )     !ooo>


       enddo ; enddo ; enddo

      end subroutine atmboundlayer_flux_dq

!..............................................................................

      subroutine atmboundlayer_matrix_dv
      implicit none

! Ajouter la contribution du flux au premier niveau de la variable
! avant de calculer tridia_in(:,:,:,4)
! A noter qu'en remplacant dt_abl_max(1) par dt_abl_max(1) on pourrait ne
! faire cette operation qu'une fois (et non pas A tous les sous-pas de
! temps). Prevoir de faire un test.....
        k=1
        do j=jbeg_,jstp_ ; do i=ibeg_,istp_
          vwinddelta_t(i,j,1)=                            & !24-11-15
          vwinddelta_t(i,j,1)                             &
          -dt_abl_max(1)/(rhoair*anyv3d(i,j,k,tnow_))     &          
            *0.5*( wstress_v(i,j,1)+wstress_v(i,j+1,1)    &
             -w0_*(wstress_v(i,j,0)+wstress_v(i,j+1,0))   &
             -w2_*(wstress_v(i,j,2)+wstress_v(i,j+1,2)))  &
                             *mask_t(i,j,kmax)            &
                          *wetmask_t(i,j)             
        enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_

         tridia_in(i,j,k,4)=                            & 

      vwinddelta_t(i,j,k)                               &
       *(x2*anyv3d(i,j,k,tnow_)+x1*anyv3d(i,j,k,tbef_)) & !teta(t)*dz(t)

         -dt_abl_max(1)*( anyv3d(i+1,j,k,id_fu_)-anyv3d(i,j,k,id_fu_)  &
                         +anyv3d(i,j+1,k,id_fv_)-anyv3d(i,j,k,id_fv_))/dxdy_t(i,j)  

       enddo ; enddo ; enddo

      end subroutine atmboundlayer_matrix_dv

!..............................................................................

      subroutine atmboundlayer_update_dv
      implicit none

      do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_
       vwinddelta_t(i,j,k)=tridia_out(i,j,k)  
      enddo ; enddo ; enddo

      end subroutine atmboundlayer_update_dv

!..............................................................................

      subroutine atmboundlayer_flux_dv
      implicit none

       x12=1./12.

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_+1

! Oi Flux de uwinddelta:
        anyv3d(i,j,k,id_fu_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_u_)*( &       !11111>
          6.*(vwinddelta_t(i  ,j,k)+vwinddelta_t(i-1,j,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_u_))*( & !22222>
          6.*(vwinddelta_t(i  ,j,k)-vwinddelta_t(i-1,j,k))  & !up2
                                  ) & !22222>

            )     !ooo>

       enddo ; enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_+1 ; do i=ibeg_,istp_

! Oj Flux de uwinddelta:
        anyv3d(i,j,k,id_fv_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_v_)*( &       !11111>
          6.*(vwinddelta_t(i,j  ,k)+vwinddelta_t(i,j-1,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_v_))*( & !22222>
          6.*(vwinddelta_t(i,j  ,k)-vwinddelta_t(i,j-1,k))  & !up2
                                  ) & !22222>

            )     !ooo>


       enddo ; enddo ; enddo

      end subroutine atmboundlayer_flux_dv

!..............................................................................

      subroutine atmboundlayer_matrix_du
      implicit none

! Ajouter la contribution du flux au premier niveau de la variable
! avant de calculer tridia_in(:,:,:,4)
! A noter qu'en remplacant dt_abl_max(1) par dt_abl_max(1) on pourrait ne
! faire cette operation qu'une fois (et non pas A tous les sous-pas de
! temps). Prevoir de faire un test.....
        k=1
        do j=jbeg_,jstp_ ; do i=ibeg_,istp_

          uwinddelta_t(i,j,1)=                            & !24-11-15
          uwinddelta_t(i,j,1)                             &
          -dt_abl_max(1)/(rhoair*anyv3d(i,j,k,tnow_))     &          
            *0.5*( wstress_u(i,j,1)+wstress_u(i+1,j,1)    &
             -w0_*(wstress_u(i,j,0)+wstress_u(i+1,j,0))   &
             -w2_*(wstress_u(i,j,2)+wstress_u(i+1,j,2)))  &
                             *mask_t(i,j,kmax)            &
                          *wetmask_t(i,j)             
        enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_

         tridia_in(i,j,k,4)=                            & 

      uwinddelta_t(i,j,k)                               &
       *(x2*anyv3d(i,j,k,tnow_)+x1*anyv3d(i,j,k,tbef_)) & !teta(t)*dz(t)

         -dt_abl_max(1)*( anyv3d(i+1,j,k,id_fu_)-anyv3d(i,j,k,id_fu_)  &
                         +anyv3d(i,j+1,k,id_fv_)-anyv3d(i,j,k,id_fv_))/dxdy_t(i,j)  

       enddo ; enddo ; enddo

      end subroutine atmboundlayer_matrix_du

!..............................................................................

      subroutine atmboundlayer_update_du
      implicit none

      do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_
       uwinddelta_t(i,j,k)=tridia_out(i,j,k)  
      enddo ; enddo ; enddo

      end subroutine atmboundlayer_update_du

!..............................................................................

      subroutine atmboundlayer_flux_du
      implicit none

       x12=1./12.

       do k=1,kmax_abl ; do j=jbeg_,jstp_ ; do i=ibeg_,istp_+1

! Oi Flux de uwinddelta:
        anyv3d(i,j,k,id_fu_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_u_)*( &       !11111>
          6.*(uwinddelta_t(i  ,j,k)+uwinddelta_t(i-1,j,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_u_))*( & !22222>
          6.*(uwinddelta_t(i  ,j,k)-uwinddelta_t(i-1,j,k))  & !up2
                                  ) & !22222>

            )     !ooo>

       enddo ; enddo ; enddo

       do k=1,kmax_abl ; do j=jbeg_,jstp_+1 ; do i=ibeg_,istp_

! Oj Flux de uwinddelta:
        anyv3d(i,j,k,id_fv_)=    &

        x12*(   & !ooo>

        anyv3d(i,j,k,id_v_)*( &       !11111>
          6.*(uwinddelta_t(i,j  ,k)+uwinddelta_t(i,j-1,k))  & !up2
                            ) &       !11111>

        -abs(anyv3d(i,j,k,id_v_))*( & !22222>
          6.*(uwinddelta_t(i,j  ,k)-uwinddelta_t(i,j-1,k))  & !up2
                                  ) & !22222>

            )     !ooo>


       enddo ; enddo ; enddo

      end subroutine atmboundlayer_flux_du

!..............................................................................

      subroutine atmboundlayer_init_debug   ! Etat initial: verifier la memoire
      implicit none

      ub4=ubound(anyv3d)
      if(ub4(3)<kmax_abl) stop 'Err 1548 atmboundlayer_init_debug'

      ub3=ubound(xy_t)
      if(ub3(3)<kmax_abl+1) stop 'Err 1549 atmboundlayer_init_debug'

      end subroutine atmboundlayer_init_debug   ! Etat initial: verifier la memoire

!..............................................................................

      end module module_atmboundlayer
     
