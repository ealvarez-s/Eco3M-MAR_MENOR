MODULE module_cpl_surfex
#ifdef key_oasis_symp_surfex
   !!===========================================================================
   !!                   ***  MODULE  toolmarsatmcpl  ***
   !!
   !! Module used for coupling applications between MARS and WRF or Meso-NH with OASIS3-MCT 
   !!
   !!===========================================================================

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC snd_fields_to_atm 
   PUBLIC rcv_fields_from_atm

  CONTAINS


   !!======================================================================

  
  SUBROUTINE snd_fields_to_atm(id_oasis_time)

   !&E---------------------------------------------------------------------
   !&E                 *** ROUTINE snd_fields_to_atm  ***
   !&E
   !&E ** Purpose : Send coupling fields to WRF or Meso-NH
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : 
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E        2013-07 : A. THEVENIN (CERFACS) : Creation 
   !&E
   !&E---------------------------------------------------------------------

   USE module_parameter,   ONLY : imax,jmax,kmax
   USE module_cpl_oasis, ONLY : snd_fld, cpl_oasis_snd, il_nb_snd                           
   USE module_principal,    ONLY : tem_t,vel_u,vel_v,gridrotsin_t,gridrotcos_t          
                          
   !! * Argument
   INTEGER, INTENT(IN) :: id_oasis_time
   
   !! * Local declaration
   REAL(kind=8), DIMENSION(imax-2,jmax-2) :: rla_oasis_snd
   INTEGER                                                   :: ib_do,i,j 
   LOGICAL                                                   :: ll_action   
   REAL(kind=8), DIMENSION(imax,jmax) :: vel_u_rot,vel_v_rot   
   !!----------------------------------------------------------------------
   !! * Executable part

!   PRINT_DBG*, 'enter SND_FIELDS_FROM_ATM routine'

   DO ib_do = 1, il_nb_snd

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP_MSK') THEN
         rla_oasis_snd(:,:) = 1.0 
      ENDIF

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP_SST') THEN
         rla_oasis_snd(:,:) = DBLE(tem_t(2:imax-1,2:jmax-1,kmax,1)+273.15)      
      ENDIF 

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP__VZ') THEN
         do j=1,jmax ; do i=1,imax
          vel_v_rot(i,j)=0.5*(                            &
           -( vel_u(i  ,j  ,kmax,1)                       &
             +vel_u(i+1,j  ,kmax,1))*gridrotsin_t(i,j)    &
           +( vel_v(i  ,j  ,kmax,1)                       &
             +vel_v(i  ,j+1,kmax,1))*gridrotcos_t(i,j))
          enddo ; enddo
         rla_oasis_snd(:,:) = DBLE(vel_v_rot(2:imax-1,2:jmax-1))
      ENDIF

      IF (snd_fld(ib_do)%cl_field_name == 'SYMP__UZ') THEN
        do j=1,jmax ; do i=1,imax
          vel_u_rot(i,j)=0.5*(                            &
            ( vel_u(i  ,j  ,kmax,1)                       &
             +vel_u(i+1,j  ,kmax,1))*gridrotcos_t(i,j)    &
           +( vel_v(i  ,j  ,kmax,1)                       &
             +vel_v(i,j+1,kmax,1))*gridrotsin_t(i,j))
         enddo ; enddo
         rla_oasis_snd(:,:) = DBLE(vel_u_rot(2:imax-1,2:jmax-1))
      ENDIF

      ! ici on envoie a chaque pas de temps, du coup oasis peut cumuler faire
      ! des moyennes... puisque oasis recoit l info mais ne la redonne pas a
      ! WRF ou Meso-NH quand ce n est pas la date de couplage. 
      ! si on cherche a limiter les copies ci-dessus, on pert cette
      ! possibilite. La limitation peut se faire en precisant en dur la pas
      ! de couplage (MODt,120) par exemple.
      CALL cpl_oasis_snd(ib_do, id_oasis_time, rla_oasis_snd, ll_action)
      
   ENDDO

!   PRINT_DBG*, 'exit SND_FIELDS_FROM_ATM routine'

  END SUBROUTINE snd_fields_to_atm


   !!======================================================================

  
  SUBROUTINE rcv_fields_from_atm(id_oasis_time)

   !&E---------------------------------------------------------------------
   !&E                 *** ROUTINE rcv_fields_from_atm  ***
   !&E
   !&E ** Purpose : Receive coupling fields from WRF or Meso-NH
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : 
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E        2013-07 : A. THEVENIN (CERFACS) : Creation 
   !&E
   !&E---------------------------------------------------------------------
  
!#include "toolcpp.h"

   USE module_parameter                  
   USE module_cpl_oasis, ONLY : rcv_fld, cpl_oasis_rcv, il_nb_rcv 
   USE module_principal
   USE module_parallele              
   !! * Argumentwstress_u
   INTEGER, INTENT(IN) :: id_oasis_time
   
   !! * Local declaration
   REAL(kind=8), DIMENSION(2:imax-1,2:jmax-1) :: rla_oasis_rcv
!   REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE :: field_rcv
   INTEGER                                                   :: ib_do
   LOGICAL                                                   :: ll_action   
   
   !!----------------------------------------------------------------------
   !! * Executable part

!   PRINT_DBG*, 'enter RCV_FIELDS_FROM_ATM routine'

   DO ib_do = 1, il_nb_rcv

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPRAIN') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            precipi_w(1:imax-2,1:jmax-2,1) = REAL(rla_oasis_rcv)
         ENDIF
      ENDIF

!    IF (rcv_fld(ib_do)%cl_field_name == 'SYMPEVAP') THEN
!      CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
!      IF (ll_action) THEN
!         evap(1:imax,1:jmax,1) = REAL(rla_oasis_rcv)
!      ENDIF
!     ENDIF
      
      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_WAT') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
         precipi_w(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
      call get_type_echange('z0','precipi_w_z0_2',precipi_w,lbound(precipi_w),ubound(precipi_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(precipi_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      precipi_w(0:imax+1,0:jmax+1,1)=precipi_w(0:imax+1,0:jmax+1,2)
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_QSR') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
      ! (QSolaire)
            ssr_w(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
      call get_type_echange('z0','ssr_w_z0_2',ssr_w,lbound(ssr_w),ubound(ssr_w),2,i2)
      do loop3=1,subcycle_exchange
       call echange_voisin(ssr_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      ssr_w(0:imax+1,0:jmax+1,1)=ssr_w(0:imax+1,0:jmax+1,2)
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPHEAT') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            ! (IR + Qlat + Qsens in one of the corresponding arrays)
            slhf_w(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
      call get_type_echange('z0','slhf_w_z0_2',slhf_w,lbound(slhf_w),ubound(slhf_w),2,i2)
      do loop3=1,subcycle_exchange
       call echange_voisin(slhf_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      slhf_w(0:imax+1,0:jmax+1,1)=slhf_w(0:imax+1,0:jmax+1,2)
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP__IR') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            snsf_w(1:imax-2,1:jmax-2,1) = REAL(rla_oasis_rcv)
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_LAT') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            slhf_w(1:imax-2,1:jmax-2,1) = REAL(rla_oasis_rcv)
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP_SEN') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            sshf_w(1:imax-2,1:jmax-2,1) = REAL(rla_oasis_rcv)
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMP__PA') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
            pss_w(1:imax-2,1:jmax-2,1) = REAL(rla_oasis_rcv)
       call get_type_echange('za','pss_w_za_2',pss_w,lbound(pss_w),ubound(pss_w),2,i2)
      do loop3=1,subcycle_exchange
       call echange_voisin(pss_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
           pss_w(-1:imax+2,-1:jmax+2,1)=pss_w(-1:imax+2,-1:jmax+2,2)           
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPTAUX') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
       taux_w(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
      call get_type_echange('za','taux_w_za_2',taux_w,lbound(taux_w),ubound(taux_w),2,i2)
      do loop3=1,subcycle_exchange
       call echange_voisin(taux_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
           taux_w(-1:imax+2,-1:jmax+2,1)=taux_w(-1:imax+2,-1:jmax+2,2)
!           taux_w(-1:imax+2,-1:jmax+2,1)=0.0
         ENDIF
      ENDIF

      IF (rcv_fld(ib_do)%cl_field_name == 'SYMPTAUY') THEN
         CALL cpl_oasis_rcv(ib_do, id_oasis_time, rla_oasis_rcv, ll_action)
         IF (ll_action) THEN
          tauy_w(2:imax-1,2:jmax-1,2) = REAL(rla_oasis_rcv)
      call get_type_echange('za','tauy_w_za_2',tauy_w,lbound(tauy_w),ubound(tauy_w),2,i2)
      do loop3=1,subcycle_exchange
      call echange_voisin(tauy_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
      tauy_w(-1:imax+2,-1:jmax+2,1)=tauy_w(-1:imax+2,-1:jmax+2,2)
!       tauy_w(-1:imax+2,-1:jmax+2,1)=0.0
         ENDIF
      ENDIF

   ENDDO

!   PRINT_DBG*, 'exit RCV_FIELDS_FROM_ATM routine'

  END SUBROUTINE rcv_fields_from_atm

#else
   !! Dummy module 
#endif
END MODULE module_cpl_surfex
