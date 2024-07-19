MODULE module_cpl_ww3
!______________________________________________________________________
! S model
! Couplage avec OASIS3-MCT2.XX
! release S26  - last update: 02-06-17
!______________________________________________________________________
!------------------------------------------------------------------------------!
! version   date      description                                              !
! S26       25-06-15  mise en service (Léo Seyfried)                           !
!           06-03-17  mise à jour pour la version 6.0 de WW3 + mise en forme   !
!                      (Fabien Rétif)                                          !
!           02-06-17  flag_2way_or_1way
!------------------------------------------------------------------------------!
#ifdef key_oasis_symp_ww3
   implicit none
   private

   !! * Accessibility
   public snd_fields_to_ww3 
   public rcv_fields_from_ww3

contains

!!======================================================================
  
  subroutine snd_fields_to_ww3(id_oasis_time)
   !==
   ! Envoi des champs à OASIS
   ! ==
   use module_parameter,   only : imax,jmax,kmax
   use module_cpl_oasis,   only : snd_fld, cpl_oasis_snd, il_nb_snd                           
   use module_principal,   only : tem_t,vel_u,vel_v,ssh_w,gridrotsin_t,gridrotcos_t,jeqjmax,jeq1,ieqimax,ieq1 &
                                 ,obcstatus,wetdry_cst2,small1,dz_t,depth_t,xy_t,k,x0,x1,x2,rho,id_u_rot      &
                                 ,id_v_rot,id_dz,id_ssh,sponge_t,velbar_u,velbar_v
   implicit none                       
   !! * Argument
   INTEGER, INTENT(IN) :: id_oasis_time
   
   !! * Local declaration
   REAL(kind=8), DIMENSION(imax-2,jmax-2) :: rla_oasis_snd
   INTEGER(kind=1) :: flag_2way_or_1way=1   ! 2 WAY (with retroaction symphonie >> WW3) !02-06-17
!  INTEGER(kind=1) :: flag_2way_or_1way=0   ! 1 WAY (no   retroaction symphonie >> WW3)
   INTEGER                                                   :: ib_do,i,j,ii,jj 
   LOGICAL                                                   :: ll_action   
!  REAL(kind=8), DIMENSION(imax,jmax) :: vel_u_rot,vel_v_rot 
   !!----------------------------------------------------------------------
   !! * Executable part

!   PRINT_DBG*, 'enter SND_FIELDS_FROM_ATM rioutine'
      id_dz=1 ; id_u_rot=2 ; id_v_rot=3

! CAS MODELE 3D (kmax>1)
      if(kmax>1) then  ! m[°v°]m >

       do j=1,jmax ; do i=1,imax
        xy_t(i,j,id_u_rot)=0.
        xy_t(i,j,id_v_rot)=0.
        xy_t(i,j,id_dz)=small1
       enddo ; enddo 
       x0=10. ! Epaisseur positive (metres) de la couche retroagissant avec les vagues 
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
          if(depth_t(i,j,k)>-x0) then !ooo>

            xy_t(i,j,id_dz)=xy_t(i,j,id_dz)+dz_t(i,j,k,1)

            xy_t(i,j,id_v_rot)=                          &
            xy_t(i,j,id_v_rot)+0.5*(                     &
             -( vel_u(i  ,j  ,k,1)                       &
               +vel_u(i+1,j  ,k,1))*gridrotsin_t(i,j)    &
             +( vel_v(i  ,j  ,k,1)                       &
               +vel_v(i  ,j+1,k,1))*gridrotcos_t(i,j))*dz_t(i,j,k,1)

            xy_t(i,j,id_u_rot)=                          &
            xy_t(i,j,id_u_rot)+0.5*(                     &
              ( vel_u(i  ,j  ,k,1)                       &
               +vel_u(i+1,j  ,k,1))*gridrotcos_t(i,j)    &
             +( vel_v(i  ,j  ,k,1)                       &
               +vel_v(i  ,j+1,k,1))*gridrotsin_t(i,j))*dz_t(i,j,k,1)

          endif                       !ooo>
       enddo ; enddo ; enddo 

      endif            ! m[°v°]m >

! CAS MODELE 2D (kmax=1)
      if(kmax==1) then ! m[0o0]m >

       do j=1,jmax ; do i=1,imax

            xy_t(i,j,id_dz)=1.

            xy_t(i,j,id_v_rot)=0.5*(                     &
             -( velbar_u(i  ,j  ,1)                      &
               +velbar_u(i+1,j  ,1))*gridrotsin_t(i,j)   &
             +( velbar_v(i  ,j  ,1)                      &
               +velbar_v(i  ,j+1,1))*gridrotcos_t(i,j))

            xy_t(i,j,id_u_rot)=0.5*(                     &
              ( velbar_u(i  ,j  ,1)                      &
               +velbar_u(i+1,j  ,1))*gridrotcos_t(i,j)   &
             +( velbar_v(i  ,j  ,1)                      &
               +velbar_v(i  ,j+1,1))*gridrotsin_t(i,j))

       enddo       ; enddo 

      endif            ! m[0o0]m >

      x1=3. ! Seuil max (m/s) appliquE au courant renvoyE au modele de vagues
      do j=1,jmax ; do i=1,imax
        xy_t(i,j,id_v_rot)=xy_t(i,j,id_v_rot)/xy_t(i,j,id_dz)
        xy_t(i,j,id_u_rot)=xy_t(i,j,id_u_rot)/xy_t(i,j,id_dz)
        x2=min(x1/sqrt(small1+xy_t(i,j,id_u_rot)**2+xy_t(i,j,id_v_rot)**2),1.)
        xy_t(i,j,id_u_rot)=x2*xy_t(i,j,id_u_rot)
        xy_t(i,j,id_v_rot)=x2*xy_t(i,j,id_v_rot)
      enddo ; enddo 

   DO ib_do = 1, il_nb_snd

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP_MSK') THEN
         rla_oasis_snd(1:imax-2,1:jmax-2) = 1.0          
      ENDIF

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP_SST') THEN
         ! On convertit en Kelvin
         rla_oasis_snd(1:imax-2,1:jmax-2) = DBLE(tem_t(2:imax-1,2:jmax-1,kmax,1)+273.15)      
      ENDIF 

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP__VZ') THEN
!        rla_oasis_snd(1:imax-2,1:jmax-2) = xy_t(2:imax-1,2:jmax-1,id_v_rot) ! DBLE(vel_v_rot(2:imax-1,2:jmax-1))
! La multiplication par (1.-sponge_t(i,j,1)) pour annuler la retroaction aux frontieres ouvertes
       do j=2,jmax-1 ; do i=2,imax-1
         rla_oasis_snd(i-1,j-1)=xy_t(i,j,id_v_rot)*(1.-sponge_t(i,j,1))*flag_2way_or_1way
       enddo         ; enddo
      ENDIF
 
      IF (snd_fld(ib_do)%cl_field_name == 'SYMP__UZ') THEN
!        rla_oasis_snd(1:imax-2,1:jmax-2) = xy_t(2:imax-1,2:jmax-1,id_u_rot) ! DBLE(vel_u_rot(2:imax-1,2:jmax-1))
! La multiplication par (1.-sponge_t(i,j,1)) pour annuler la retroaction aux frontieres ouvertes
       do j=2,jmax-1 ; do i=2,imax-1
         rla_oasis_snd(i-1,j-1)=xy_t(i,j,id_u_rot)*(1.-sponge_t(i,j,1))*flag_2way_or_1way
       enddo         ; enddo
      ENDIF

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP_SSH') THEN
!        rla_oasis_snd(1:imax-2,1:jmax-2) = DBLE(ssh_w(2:imax-1,2:jmax-1,1)-wetdry_cst2)
!        rla_oasis_snd(1:imax-2,1:jmax-2) = 0.0
! La multiplication par (1.-sponge_t(i,j,1)) pour annuler la retroaction aux frontieres ouvertes
       do j=2,jmax-1 ; do i=2,imax-1
         rla_oasis_snd(i-1,j-1)=(ssh_w(i,j,1)-wetdry_cst2)*(1.-sponge_t(i,j,1))*flag_2way_or_1way
       enddo         ; enddo
      ENDIF

      CALL cpl_oasis_snd(ib_do, id_oasis_time, rla_oasis_snd, ll_action)
      
   ENDDO

!   PRINT_DBG*, 'exit SND_FIELDS_FROM_ATM routine'

  end subroutine snd_fields_to_ww3


   !!======================================================================

  
  SUBROUTINE rcv_fields_from_ww3(id_oasis_time)

   USE module_parameter                  
   USE module_cpl_oasis, ONLY : rcv_fld, cpl_oasis_rcv, il_nb_rcv 
   USE module_principal
   USE module_parallele              
   !! * Argumentwstress_u
   INTEGER, INTENT(IN) :: id_oasis_time
   !! * Local declaration
   REAL(kind=8), DIMENSION(2:imax-1,2:jmax-1) :: rla_oasis_rcv
!   REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE :: field_rcv
   INTEGER                                                   :: ib_do, ii, jj
   LOGICAL                                                   :: ll_action   
   
   !!----------------------------------------------------------------------
   !! * Executable part

!   PRINT_DBG*, 'enter RCV_FIELDS_FROM_ATM routine'
!

   DO ib_do = 1, il_nb_rcv
      
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPSDIR') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            xy_t(2:imax-1,2:jmax-1,6) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_sdir=1 ! update receive flag
         ENDIF
      ENDIF
      
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPCDIR') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            xy_t(2:imax-1,2:jmax-1,5) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_cdir=1 ! update receive flag 
         ENDIF
      ENDIF
      
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP__HS') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            hs_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_hs=1 ! update receive flag
         ENDIF
      ENDIF
      
       IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_HSW') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            hsw_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_hsw=1 ! update receive flag
         ENDIF
      ENDIF
 
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_FOC') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            foc_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_foc=1 ! update receive flag
         ENDIF
      ENDIF
       
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP__TW') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            t_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_tw=1 ! update receive flag
         ENDIF
      ENDIF
 
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPTAWX') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            tawx_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_tawx=1 ! update receive flag
         ENDIF
      ENDIF
     
     IF (rcv_fld(ib_do)%cl_field_name == 'SYMPTAWY') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            tawy_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_tawy=1 ! update receive flag
         ENDIF
      ENDIF

 
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPTWOX') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            twox_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_twox=1 ! update receive flag
         ENDIF
      ENDIF
      
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPTWOY') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            twoy_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_twoy=1 ! update receive flag
         ENDIF
      ENDIF
 
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPUSSX') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            uss_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_uss=1 ! update receive flag
         ENDIF
      ENDIF
 
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPUSSY') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            vss_wave_t(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_vss=1 ! update receive flag
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_MSK') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            mask_wave_t(2:imax-1,2:jmax-1) = REAL(rla_oasis_rcv)
            wave_cpl_ww3_msk=1 ! update receive flag
         ENDIF
      ENDIF

   ENDDO
!   PRINT_DBG*, 'exit RCV_FIELDS_FROM_ATM routine'

  END SUBROUTINE rcv_fields_from_ww3

#else
   !! Dummy module 
#endif
END MODULE module_cpl_ww3
