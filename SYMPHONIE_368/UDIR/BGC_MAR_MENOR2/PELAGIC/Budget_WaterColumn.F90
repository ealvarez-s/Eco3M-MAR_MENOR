      SUBROUTINE BUDGET_WATERCOLUMN


!_____________________________________________________________________*
!                                                                     *
! LAST REVISION: 31 AOUT 2009                                         *
!                                                                     *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Calculates water column budget                                      *
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

      USE MODULE_PRINCIPAL
      USE ModuleDeclaration
      USE ModuleComBlock
      use module_parallele !#MPI
      use module_global
      IMPLICIT NONE

! Local variables

      INTEGER kountMOY

!---------------------------------------------------------------------*


! DEB BILAN PRODUCTION - RESPIRATION sur le plateau

! cumul
! claude tous les BIL_ pourraient etre apres le test sur kount qui permet d'ecrire puisque le cumul
! ne se fait pas ici
       BIL_PPB =  SUMPPBGDL_plat   / 86400. /12.                    &
                             *dti_fw        ! en mmolC

!      print*,'dti_fw Bud',dti_fw 


       BIL_NETPPB =  SUMNETPPBGDL_plat   / 86400. /12.              &
                             *dti_fw

      BIL_RESPPHYTO =  SUMRESPPHYTOGDL_plat   / 86400. /12.         &
                             *dti_fw 

      BIL_RESPZOO =  SUMRESPZOOGDL_plat   / 86400. /12.             &
                             *dti_fw

      BIL_RESPBACT =  (SUMRESPTOTGDL_plat                           &
                     - SUMRESPZOOGDL_plat                           &
                     - SUMRESPPHYTOGDL_plat )  / 86400. /12.        &
                             *dti_fw

      BIL_REMPOC =  SUMREMPOCGDL_plat   / 86400. /12.               &
                             *dti_fw

      BIL_GBAC =  SUMGBACGDL_plat   / 86400. /12.                   &
                             *dti_fw

      BIL_LOSTPOC =  SUMLOSTPOCGDL_plat   / 86400. /12.             &
                             *dti_fw


       BIL_NPB =  SUMNPBGDL_plat   / 86400. /14.                    &
                            *dti_fw

       BIL_RPB =  SUMRPBGDL_plat   / 86400. /14.                    &
                            *dti_fw

      BIL_REMPON =  SUMREMPONGDL_plat   / 86400. /14.               &
                            *dti_fw

      BIL_LOSTPON =  SUMLOSTPONGDL_plat   / 86400. /14.             &
                            *dti_fw

      BIL_UPTBACTDON =  SUMUPTBACTDONGDL_plat   / 86400. /14.       &
                            *dti_fw

      BIL_EXCHETERON =  SUMEXCHETERONGDL_plat   / 86400. /14.       &
                            *dti_fw


! instantane
      INS_PPB = INS_PPB +  INSPPBGDL_plat   / 86400. /12.           &
                            *dti_fw        ! en mmolC

      INS_NETPPB =  INS_NETPPB + INSNETPPBGDL_plat   / 86400. /12.  &
                            *dti_fw        ! en mmolC


      INS_RESPPHYTO =INS_RESPPHYTO + INSRESPPHYTOGDL_plat/ 86400./12.  &
                            *dti_fw

      INS_RESPZOO =  INS_RESPZOO + INSRESPZOOGDL_plat   / 86400. /12.  &
                            *dti_fw

      INS_RESPBACT =  INS_RESPBACT + (INSRESPTOTGDL_plat               &
                     - INSRESPZOOGDL_plat                              &
                     - INSRESPPHYTOGDL_plat )  / 86400. /12.           &
                            *dti_fw

      INS_GBAC =  INS_GBAC + INSGBACGDL_plat   / 86400. /12.           &
                            *dti_fw

      INS_LOSTPOC =  INS_LOSTPOC + INSLOSTPOCGDL_plat   / 86400. /12.  &
                            *dti_fw


       INS_NPB =  INS_NPB + INSNPBGDL_plat   / 86400. /14.             &
                            *dti_fw

       INS_RPB =  INS_RPB + INSRPBGDL_plat   / 86400. /14.             &
                            *dti_fw

      INS_LOSTPON =  INS_LOSTPON + INSLOSTPONGDL_plat   / 86400. /14.  &
                            *dti_fw

      INS_UPTBACTDON =INS_UPTBACTDON + INSUPTBACTDONGDL_plat/86400./14. &
                            *dti_fw

      INS_EXCHETERON =INS_EXCHETERON + INSEXCHETERONGDL_plat/86400./14. &
                            *dti_fw

! FIN BILAN PRODUCTION - RESPIRATION sur le plateau


! Bilan sediment sur le plateau du golfe du Lion
      SUM_HBIO_P=0.
      SUM_HBION_P=0.
      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
      SUM_HBIO_P=SUM_HBIO_P - SUM_EXPORTC_BOT_2D(I,J)  !  /86400.  Claude: me semble qu'il faut pas diviser par 86400
      SUM_HBION_P=SUM_HBION_P - SUM_EXPORTN_BOT_2D(I,J) ! /86400.   voir ds Budget_Export_Bottom.F90 versions diff Caroline et Pierre
       endif                                                             !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       endif
       endif
      enddo
      enddo

!     SUM_HBIO_P=SUM_HBIO_P*DXB*DYB
!     SUM_HBION_P=SUM_HBION_P*DXB*DYB
! Bilan sediment sur le plateau du golfe du Lion



! DEB BILAN COLONNE D'EAU
! Intégrale sur le domaine de la hauteur de sediment et concentration
! PAR CLASSE

      kountMOY = 12

      if(MOD(kount,kountMOY)==0)then
      

      sum_biopoc_p=0.
      sum_biodoc_p=0.
      sum_biopon_p=0.
      sum_biodon_p=0.
      sum_biono3_p=0.
      sum_bionh3_p=0.

      DO vb=1,vbmax

! start stock of POC
      if(vb==IZOONANOC.OR.               &
         vb==IZOOMICROC.OR.              &
         vb==IZOOMESOC.OR.               &
         vb==ISYNEC.OR.                  &
         vb==INANOC.OR.                  &
         vb==IDIAC.OR.                   &
         vb==IBACTC.OR.                  &
         vb==ILMOPC.OR.                  &
         vb==ISMOPC     )   then  !------>


      sumbio_p(vb)=0.

      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
           do k=1,kmax
       x1=dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
                     sumbio_p(vb)=sumbio_p(vb)                     &
                          +bio_t(i,j,k,vb)*x1
           enddo
       endif
       endif
       endif
      enddo  !J
      enddo  !I

      sum_biopoc_p=sum_biopoc_p+sumbio_p(vb)
      endif                         !------>
! end stock of POC

! start stock of DOC
      if(vb.EQ.IMODC      )   then  !------>


      sumbio_p(vb)=0.

      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
       do k=1,kmax
       x1=dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
                     sumbio_p(vb)=sumbio_p(vb)                     &
                          +bio_t(i,j,k,vb)*x1
           enddo
       endif                                                              !vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       endif
       endif
      enddo  !J
      enddo  !I

      sum_biodoc_p=sum_biodoc_p+sumbio_p(vb)
      endif                         !------>
! end stock of DOC


! start stock of PON
      if(vb==IZOONANOC.OR.     &
         vb==IZOOMICROC.OR.    &
         vb==IZOOMESOC.OR.     &
         vb==ISYNEN.OR.        &
         vb==INANON.OR.        &
         vb==IDIAN.OR.         &
         vb==IBACTC.OR.        &
         vb==ILMOPN.OR.        &
         vb==ISMOPN     )   then  !------>


      sumbion_p(vb)=0.

      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
           do k=1,kmax
       x1=dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
                     sumbion_p(vb)=sumbion_p(vb)                     &
                          +bio_t(i,j,k,vb)*x1
           enddo
       endif                                                              !vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       endif
       endif
      enddo  !J
      enddo  !I


! bilan général du stock de PON
       if(vb==IZOONANOC.OR.    &
          vb==IZOOMICROC.OR.   &
          vb==IZOOMESOC)       &
        sumbion_p(vb)=sumbion_p(vb)*NCZooNano
       if(vb==IBACTC) sumbion_p(vb)=sumbion_p(vb)*NCBact
        sum_biopon_p=sum_biopon_p+sumbion_p(vb)

      endif 
! end stock of PON

! start stock of DON
      if(vb==IMODN      )   then  !------>

      sumbion_p(vb)=0.

      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
         do k=1,kmax
         x1=dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
                     sumbion_p(vb)=sumbion_p(vb)                     &
                          +bio_t(i,j,k,vb)*x1
         enddo
       endif                                                             !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       endif
       endif
      enddo  !J
      enddo  !I

! bilan général du stock de DON
      sum_biodon_p=sum_biodon_p+sumbion_p(vb)
      endif
! end stock of DON

! start stock of NO3
      if(vb==INITRATE) then

      sumbion_p(vb)=0.

      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
           do k=1,kmax
         x1=dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
                     sumbion_p(vb)=sumbion_p(vb)                     &
                          +bio_t(i,j,k,vb)*x1
           enddo
       endif                                                             !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       endif
       endif
      enddo  !J
      enddo  !I
! bilan général du stock de carbone
      sum_biono3_p=sum_biono3_p+sumbion_p(vb)

      endif
! end stock of NO3

! start stock of NH4
      if(vb==IAMMONIUM) then

      sumbion_p(vb)=0.

      do j=nbio1,nbio2
      do i=mbio1,mbio2
       i2=i+par%timax(1)      ! ds le repere global
       j2=j+par%tjmax(1)
       if(mask_t(i,j,kmax+1)==1) then
       if (j2>=j1gdl.and.j2<=j2gdl.and.i2>=i1gdl.and.i2<=i2gdl) then     !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       if(-depth_w(i,j,kmin_w(i,j))<=200) then
           do k=1,kmax
         x1=dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
                     sumbion_p(vb)=sumbion_p(vb)                     &
                          +bio_t(i,j,k,vb)*x1
           enddo
       endif                                                             !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
       endif
       endif
      enddo  !J
      enddo  !I
! bilan général du stock de carbone
      sum_bionh3_p=sum_bionh3_p+sumbion_p(vb)

      endif ! test sur vb          !------>
! end stock of NH4

      enddo  !vb

! FIN BILAN COLONNE D'EAU



! Ecriture du fichier de Bilan debut
! on somme sur tous les procs
#ifdef parallele
! on peut sommer car ces variables ne se cumulent pas sur toute la simu
      call mpi_allreduce(BIL_NETPPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_NETPPB=sum1glb
      call mpi_allreduce(BIL_PPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_PPB=sum1glb
      call mpi_allreduce(BIL_TOUTFLEUVEPOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_TOUTFLEUVEPOC=sum1glb
      call mpi_allreduce(BIL_TOUTFLEUVEDOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_TOUTFLEUVEDOC=sum1glb
      call mpi_allreduce(SUM_HBIO_P,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      SUM_HBIO_P=sum1glb
      call mpi_allreduce(sum_biopoc_p,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_biopoc_p=sum1glb
      call mpi_allreduce(sum_biodoc_p,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_biodoc_p=sum1glb
      call mpi_allreduce(BIL_RESPPHYTO,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_RESPPHYTO=sum1glb
      call mpi_allreduce(BIL_REMPOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_REMPOC=sum1glb
      call mpi_allreduce(BIL_GBAC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_GBAC=sum1glb
      call mpi_allreduce(BIL_LOSTPOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_LOSTPOC=sum1glb
      call mpi_allreduce(BIL_RESPZOO,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_RESPZOO=sum1glb
      call mpi_allreduce(BIL_RESPBACT,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_RESPBACT=sum1glb
#endif
       if(par%rank==0)then
       open(unit=3,file='BILAN_GENERAL',position='append')
       write(3,104)float(kount)*dti_fw/86400.  &
              ,BIL_NETPPB         &    ! gain par ppb - perte par respiration, mort, exsu cum
              ,BIL_PPB            &    ! gain par ppb en mmolC
              ,BIL_TOUTFLEUVEPOC  &    ! apport des fleuves cum en mmolC
              ,BIL_TOUTFLEUVEDOC  &    ! apport des fleuves cum en mmolC
              ,SUM_HBIO_P         &    ! depot dans le sediment cum
              ,sum_biopoc_p       &    ! stock dans la colonne d'eau de POC
              ,sum_biodoc_p       &    ! stock dans la colonne d'eau de DOC
              ,BIL_RESPPHYTO      &    ! respiration du phyto
              ,BIL_REMPOC         &    ! remineralisation du POC
              ,BIL_GBAC           &    ! croissance du phyto
              ,BIL_LOSTPOC        &    ! bacteria mortality + Messy feeding
              ,BIL_RESPZOO        &    ! respiration of zoo
              ,BIL_RESPBACT            ! respiration of bact 

      close(3)
      endif

#ifdef parallele
      call mpi_allreduce(BIL_NPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_NPB=sum1glb
      call mpi_allreduce(BIL_RPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_RPB=sum1glb
      call mpi_allreduce(BIL_TOUTFLEUVEPON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_TOUTFLEUVEPON=sum1glb
      call mpi_allreduce(BIL_TOUTFLEUVEDON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_TOUTFLEUVEDON=sum1glb
      call mpi_allreduce(BIL_TOUTFLEUVENO3,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_TOUTFLEUVENO3=sum1glb
      call mpi_allreduce(BIL_TOUTFLEUVENH3,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_TOUTFLEUVENH3=sum1glb
      call mpi_allreduce(SUM_HBION_P,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      SUM_HBION_P=sum1glb
      call mpi_allreduce(sum_biopon_p,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_biopon_p=sum1glb
      call mpi_allreduce(sum_biodon_p,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_biodon_p=sum1glb
      call mpi_allreduce(sum_biono3_p,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_biono3_p=sum1glb
      call mpi_allreduce(sum_bionh3_p,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_bionh3_p=sum1glb
      call mpi_allreduce(BIL_REMPON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_REMPON=sum1glb
      call mpi_allreduce(BIL_LOSTPON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_LOSTPON=sum1glb
      call mpi_allreduce(BIL_UPTBACTDON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_UPTBACTDON=sum1glb
      call mpi_allreduce(BIL_EXCHETERON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      BIL_EXCHETERON=sum1glb
#endif

      if(par%rank==0)then
      open(unit=3,file='BILAN_GENERALN',position='append')
      write(3,105)float(kount)*dti_fw/86400.           &
              ,BIL_NPB            &   ! uptake NO3
              ,BIL_RPB            &   ! uptake NH3
              ,BIL_TOUTFLEUVEPON  &   ! apport des fleuves cum PON
              ,BIL_TOUTFLEUVEDON  &   ! apport des fleuves cum DON
              ,BIL_TOUTFLEUVENO3  &   ! apport des fleuves cum NO3
              ,BIL_TOUTFLEUVENH3  &   ! apport des fleuves cum NH3
              ,SUM_HBION_P        &   ! depot dans le sediment cum
              ,sum_biopon_p       &   ! stock dans la colonne d'eau PON
              ,sum_biodon_p       &   ! stock dans la colonne d'eau DON
              ,sum_biono3_p       &   ! stock dans la colonne d'eau NO3
              ,sum_bionh3_p       &   ! stock dans la colonne d'eau NH3
              ,BIL_REMPON         &   ! remineralisation du PON
              ,BIL_LOSTPON        &   ! bacteria mortality + Messy feeding + Phyto Exsudation
              ,BIL_UPTBACTDON     &   ! uptake of DON by bacteria
              ,BIL_EXCHETERON         ! excretion of DIN by heterophs 
      close(3)
      endif


#ifdef parallele
      call mpi_allreduce(INS_NETPPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_NETPPB=sum1glb
      call mpi_allreduce(INS_PPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_PPB=sum1glb
      call mpi_allreduce(INS_TOUTFLEUVEPOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_TOUTFLEUVEPOC=sum1glb
      call mpi_allreduce(INS_TOUTFLEUVEDOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_TOUTFLEUVEDOC=sum1glb
      call mpi_allreduce(INS_RESPPHYTO,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_RESPPHYTO=sum1glb
      call mpi_allreduce(INS_GBAC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_GBAC=sum1glb
      call mpi_allreduce(INS_LOSTPOC,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_LOSTPOC=sum1glb
      call mpi_allreduce(INS_RESPZOO,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_RESPZOO=sum1glb
      call mpi_allreduce(INS_RESPBACT,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_RESPBACT=sum1glb
! les tableaux _p1 ne sont pas sommés sur les procs car ce sont des buffers des _p qui sont deja sommés
#endif

       if(par%rank==0)then
       open(unit=3,file='insBILAN_GENERAL',position='append')
       write(3,106)float(kount)*dti_fw/86400.   &
              ,INS_NETPPB / kountMOY            &  ! gain par ppb - perte par respiration, mort, exsu cum
              ,INS_PPB / kountMOY               &  ! gain par ppb en mmolC
              ,INS_TOUTFLEUVEPOC / kountMOY     &  ! apport des fleuves cum en mmolC
              ,INS_TOUTFLEUVEDOC / kountMOY     &  ! apport des fleuves cum en mmolC
              ,(SUM_HBIO_P-SUM_HBIO_P1) / kountMOY  &   ! depot dans le sediment cum
              ,(sum_biopoc_p-sum_biopoc_p1) / kountMOY &     ! stock dans la colonne d'eau de POC
              ,(sum_biodoc_p-sum_biodoc_p1) / kountMOY &     ! stock dans la colonne d'eau de DOC
              ,INS_RESPPHYTO / kountMOY         &  ! respiration du phyto
              ,INS_GBAC / kountMOY              &  ! croissance du phyto
              ,INS_LOSTPOC / kountMOY           &  ! bacteria mortality + Messy feeding
              ,INS_RESPZOO / kountMOY           &  ! respiration of zoo
              ,INS_RESPBACT / kountMOY             ! respiration of bact
      close(3)
      endif

#ifdef parallele
      call mpi_allreduce(INS_NPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_NPB=sum1glb
      call mpi_allreduce(INS_RPB,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_RPB=sum1glb
      call mpi_allreduce(INS_TOUTFLEUVEPON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_TOUTFLEUVEPON=sum1glb
      call mpi_allreduce(INS_TOUTFLEUVEDON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_TOUTFLEUVEDON=sum1glb
      call mpi_allreduce(INS_TOUTFLEUVENO3,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_TOUTFLEUVENO3=sum1glb
      call mpi_allreduce(INS_TOUTFLEUVENH3,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_TOUTFLEUVENH3=sum1glb
      call mpi_allreduce(INS_LOSTPON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_LOSTPON=sum1glb
      call mpi_allreduce(INS_UPTBACTDON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_UPTBACTDON=sum1glb
      call mpi_allreduce(INS_EXCHETERON,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      INS_EXCHETERON=sum1glb
#endif
       if(par%rank==0)then
       open(unit=3,file='insBILAN_GENERALN',position='append')
       write(3,107)float(kount)*dti_fw/86400.     &
              ,INS_NPB / kountMOY                        &  ! uptake NO3
              ,INS_RPB / kountMOY                        &  ! uptake NH3
              ,INS_TOUTFLEUVEPON / kountMOY              &  ! apport des fleuves cum PON
              ,INS_TOUTFLEUVEDON / kountMOY              &  ! apport des fleuves cum DON
              ,INS_TOUTFLEUVENO3 / kountMOY              &  ! apport des fleuves cum NO3
              ,INS_TOUTFLEUVENH3 / kountMOY              &  ! apport des fleuves cum NH3
              ,(SUM_HBION_P-SUM_HBION_P1) / kountMOY     &  ! depot dans le sediment cum
              ,(sum_biopon_p-sum_biopon_p1) / kountMOY   &  ! stock dans la colonne d'eau PON
              ,(sum_biodon_p-sum_biodon_p1) / kountMOY   &  ! stock dans la colonne d'eau DON
              ,(sum_biono3_p-sum_biono3_p1) / kountMOY   &  ! stock dans la colonne d'eau NO3
              ,(sum_bionh3_p-sum_bionh3_p1) / kountMOY   &  ! stock dans la colonne d'eau NH3
              ,INS_LOSTPON  / kountMOY           &  ! bacteria mortality + Messy feeding + Phyto Exsudation
              ,INS_UPTBACTDON  / kountMOY                &  ! uptake of DON by bacteria
              ,INS_EXCHETERON / kountMOY                    ! excretion of DIN by heterophs
        endif


       close(3)

! Remise a zero apres moyenne et ecriture
      INS_NETPPB = 0.
      INS_PPB    = 0.
      INS_RESPPHYTO = 0. 
      INS_GBAC = 0.
      INS_LOSTPOC = 0.
      INS_RESPZOO = 0.  
      INS_RESPBACT = 0. 

      INS_NPB = 0.
      INS_RPB = 0.  
      INS_LOSTPON  = 0.
      INS_UPTBACTDON = 0.
      INS_EXCHETERON = 0.

! Stock dans la colonne d'eau
      sum_biopoc_p1=sum_biopoc_p      ! claude ?? on ecrit plus haut la difference des deux??
      sum_biodoc_p1=sum_biodoc_p
      sum_biopon_p1=sum_biopon_p
      sum_biodon_p1=sum_biodon_p
      sum_biono3_p1=sum_biono3_p
      sum_bionh3_p1=sum_bionh3_p

! Depot dans le sediment
      SUM_HBIO_P1=SUM_HBIO_P
      SUM_HBION_P1=SUM_HBION_P

! Fleuves
      INS_TOUTFLEUVEPOC=0.
      INS_TOUTFLEUVEDOC=0.
      INS_TOUTFLEUVEPON=0.
      INS_TOUTFLEUVEDON=0.
      INS_TOUTFLEUVENO3=0.
      INS_TOUTFLEUVENH3=0.

      endif ! mod kount

! Remise a zero a tous les kount
      INSNETPPBGDL_plat = 0.
      INSPPBGDL_plat    = 0.
      INSRESPPHYTOGDL_plat = 0.
      INSGBACGDL_plat = 0.
      INSLOSTPOCGDL_plat = 0.
      INSRESPZOOGDL_plat = 0.
      INSRESPBACTGDL_plat = 0.

      INSNPBGDL_plat = 0.
      INSRPBGDL_plat = 0.
      INSLOSTPONGDL_plat  = 0.
      INSUPTBACTDONGDL_plat = 0.
      INSEXCHETERONGDL_plat = 0.



! Ecriture du fichier de Bilan fin

  104   format(14(E14.5,1X))
  105   format(16(E14.5,1X))
  106   format(13(E14.5,1X))
  107   format(15(E14.5,1X))

      RETURN
      END

