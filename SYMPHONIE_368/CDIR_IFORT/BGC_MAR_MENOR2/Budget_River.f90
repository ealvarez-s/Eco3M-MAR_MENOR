










      SUBROUTINE BUDGET_RIVER


!_____________________________________________________________________*
!                                                                     *
! LAST REVISION: 24 JANUARY 2008                                      *
!                                                                     *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Calculates river inputs                                             *
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

      use module_principal
      use ModuleDeclaration
      use module_parallele !#MPI
      use module_global
      implicit none

! Local variables
!      double precision bil_fleuves_area(vbmax)

!---------------------------------------------------------------------*


      const1=dti_fw*0.5
      const2=dti_fw*0.5

      do vb=1,vbmax
       ins_fleuves(vb)=0.
      enddo

      do kr=1,nriver

!     if(river_dom(kr) == par%rank) then 
      if(rivertrc_inout(kr)==1) then
      if(riverdir(kr)/=0)then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      i=iriver(kr,2)
      j=jriver(kr,2)
!     print*,'river kr',kr,i,j 

      if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then

      if(riverdir(kr)==2) then ! ----------------->

      x3=mask_i_w(i)*mask_j_w(j)
      x4=mask_i_w(i)*mask_j_w(j-1)

      do vb=1,vbmax
      do k=kmin_v(i,j),kmax

!      if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

      ins_fleuve(kr,vb)=-                                  &             ! unite : mmolX
       veldxdz_v(i,j,k,1)*const1*(           &
        sign(1.,veldxdz_v(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4           )

      bil_fleuve(kr,vb)=bil_fleuve(kr,vb)+                 &             ! unite : mmolX
        ins_fleuve(kr,vb)

      ins_fleuves(vb)=-                                    &
       veldxdz_v(i,j,k,1)*const1*(           &
        sign(1.,veldxdz_v(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4           )

      bil_fleuves(vb)=bil_fleuves(vb)+                     &
        ins_fleuves(vb)

      bil_fleuves_area(vb)=bil_fleuves_area(vb)+           &
        ins_fleuves(vb)

!      endif

      enddo
      enddo


      elseif(riverdir(kr)==4) then ! ----------------->

      x3=mask_i_w(i)*mask_j_w(j)
      x4=mask_i_w(i)*mask_j_w(j-1)

      do vb=1,vbmax
      do k=kmin_v(i,j),kmax

!      if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

      ins_fleuve(kr,vb)=                                   &             ! unite : mmolX
       veldxdz_v(i,j,k,1)*const1*(           &
        sign(1.,veldxdz_v(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4           )

      bil_fleuve(kr,vb)=bil_fleuve(kr,vb)+                 &
        ins_fleuve(kr,vb)


      ins_fleuves(vb)=                                     &
       veldxdz_v(i,j,k,1)*const1*(           &
        sign(1.,veldxdz_v(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4           )

      bil_fleuves(vb)=bil_fleuves(vb)+                     &
        ins_fleuves(vb)

      bil_fleuves_area(vb)=bil_fleuves_area(vb)+           &
        ins_fleuves(vb)

!     endif

      enddo
      enddo


      elseif(riverdir(kr)==1) then ! ----------------->

      x3=mask_i_w(i)*mask_j_w(j)
      x4=mask_i_w(i-1)*mask_j_w(j)

      do vb=1,vbmax
      do k=kmin_u(i,j),kmax

!      if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

      ins_fleuve(kr,vb)=-                                  &
       veldydz_u(i,j,k,1)*const2*(           &
        sign(1.,veldydz_u(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4           )

      bil_fleuve(kr,vb)=bil_fleuve(kr,vb)+                 &
        ins_fleuve(kr,vb)

      ins_fleuves(vb)=-                                    &
       veldydz_u(i,j,k,1)*const2*(           &
        sign(1.,veldydz_u(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4           )

      bil_fleuves(vb)=bil_fleuves(vb)+                     &
        ins_fleuves(vb)

      bil_fleuves_area(vb)=bil_fleuves_area(vb)+           &
        ins_fleuves(vb)

!     endif 

      enddo
      enddo

      elseif(riverdir(kr)==3) then ! ----------------->

      x3=mask_i_w(i)*mask_j_w(j)
      x4=mask_i_w(i-1)*mask_j_w(j)

      do vb=1,vbmax
      do k=kmin_u(i,j),kmax

!       if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

      ins_fleuve(kr,vb)=                                   &
       veldydz_u(i,j,k,1)*const2*(           &
        sign(1.,veldydz_u(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4           )

      bil_fleuve(kr,vb)=bil_fleuve(kr,vb)+                 &
        ins_fleuve(kr,vb)

      ins_fleuves(vb)=                                     &
       veldydz_u(i,j,k,1)*const2*(           &
        sign(1.,veldydz_u(i,j,k,1))                        &
       *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)          &
        -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4           )

      bil_fleuves(vb)=bil_fleuves(vb)+                     &
        ins_fleuves(vb)

      bil_fleuves_area(vb)=bil_fleuves_area(vb)+           &
        ins_fleuves(vb)

!       endif
  
      enddo
      enddo

      endif !---------------------------------------------->

      else     ! riverdir(kr)=0
      if(flag_nemoffline==1) then    ! cas depot riviere traité par flux de surface
      i=iriver(kr,1)
      j=jriver(kr,1)
      IF(lon_t(i,j)*rad2deg>-5.6) then 

      if ((lon_t(i,j)*rad2deg<10).or.    &
        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<15.and.lat_t(i,j)*rad2deg>37.and.lat_t(i,j)*rad2deg<42).or.    &
        (lon_t(i,j)*rad2deg>10.and.lon_t(i,j)*rad2deg<12.25.and.lat_t(i,j)*rad2deg>42.and.lat_t(i,j)*rad2deg<44.25).or.    &
        (lon_t(i,j)*rad2deg>15.and.lon_t(i,j)*rad2deg<16.25.and.lat_t(i,j)*rad2deg>38.and.lat_t(i,j)*rad2deg<40.25))then

      x3=mask_i_w(i)*mask_j_w(j)
      do vb=1,vbmax
      ins_fleuves(vb)=fluxbio_w(i,j,vb,2)*x3*dxdy_t(i,j) 
      bil_fleuves(vb)=bil_fleuves(vb)+                     &
        ins_fleuves(vb)

      bil_fleuves_area(vb)=bil_fleuves_area(vb)+           &
        ins_fleuves(vb)
        enddo

        else   

      x3=mask_i_w(i)*mask_j_w(j)
      do vb=1,vbmax
      ins_fleuves_est(vb)=fluxbio_w(i,j,vb,2)*x3*dxdy_t(i,j)
      bil_fleuves_est(vb)=bil_fleuves_est(vb)+                     &
        ins_fleuves_est(vb)

      bil_fleuves_area_est(vb)=bil_fleuves_area_est(vb)+           &
        ins_fleuves_est(vb)
        enddo
        endif   ! end est

        ENDIF ! 5.6
      endif    ! riverdir(kr)=0

      endif ! i,j 

      endif                     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!        write(*,*)'bil_fleuves_area_est_nit BOUCLE = ',bil_fleuves_area_est(30)
!        write(*,*)'ins_fleuves_est_nit BOUCLE = ', ins_fleuves_est(30)

!       write(*,*)'bil_fleuves_area_west_nit BOUCLE = ',bil_fleuves_area(30)
!        write(*,*)'ins_fleuves_west_nit BOUCLE = ', ins_fleuves(30)


! bilan general de tous les fleuves, toutes les classes:
      BIL_TOUTFLEUVEPOC=0.
      BIL_TOUTFLEUVEDOC=0.
      BIL_TOUTFLEUVEPON=0.
      BIL_TOUTFLEUVEDON=0.
      BIL_TOUTFLEUVENO3=0.
      BIL_TOUTFLEUVENH3=0.

!      do vb=1,vbmax
!
!      IF(vb==ILMOPC.or.                                     &
!         vb==ISMOPC     )   then !---->
!
!      BIL_TOUTFLEUVEPOC=BIL_TOUTFLEUVEPOC+BIL_FLEUVES(vb)
!      INS_TOUTFLEUVEPOC=INS_TOUTFLEUVEPOC+INS_FLEUVES(vb)
!
!!     elseif(vb==IMODCL.or.vb==IMODCSL      )   then !---->                 &
!      elseif(vb==IMODC) then !---->                 &
!
!      BIL_TOUTFLEUVEDOC=BIL_TOUTFLEUVEDOC+BIL_FLEUVES(vb)
!      INS_TOUTFLEUVEDOC=INS_TOUTFLEUVEDOC+INS_FLEUVES(vb)
!
!      elseif(vb==ILMOPN.or.                                 &
!         vb==ISMOPN     )   then !---->
!
!
!      BIL_TOUTFLEUVEPON=BIL_TOUTFLEUVEPON+BIL_FLEUVES(vb)
!      INS_TOUTFLEUVEPON=INS_TOUTFLEUVEPON+INS_FLEUVES(vb)
!
!!     elseif(vb==IMODNL.or.vb==IMODNSL      )   then !---->                 
!      elseif(vb==IMODN) then !---->
!
!      BIL_TOUTFLEUVEDON=BIL_TOUTFLEUVEDON+BIL_FLEUVES(vb)
!      INS_TOUTFLEUVEDON=INS_TOUTFLEUVEDON+INS_FLEUVES(vb)
!
!      elseif(vb==INITRATE) then !---->
!
!      BIL_TOUTFLEUVENO3=BIL_TOUTFLEUVENO3+BIL_FLEUVES(vb)
!      INS_TOUTFLEUVENO3=INS_TOUTFLEUVENO3+INS_FLEUVES(vb)
!
!      elseif(vb==IAMMONIUM) then !---->
!
!     BIL_TOUTFLEUVENH3=BIL_TOUTFLEUVENH3+BIL_FLEUVES(vb)
!      INS_TOUTFLEUVENH3=INS_TOUTFLEUVENH3+INS_FLEUVES(vb)
!
!      endif


!      enddo          !vb


      endif          ! river_dom

      enddo          ! kr


      if(mod(iteration3d,12)==0)then
!èèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèè

!     write(6,*)par%rank,elapsedtime_now,bil_fleuves_area(30)/12.
      do vb=1,vbmax
      call mpi_allreduce(bil_fleuves_area(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      bil_fleuves_area(vb)=sum1glb
!      enddo
!
!      do vb=1,vbmax
      call mpi_allreduce(bil_fleuves_area_est(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      bil_fleuves_area_est(vb)=sum1glb
      enddo


! seul un proc travaille pour les dernieres operations sur les quantités
! globales et l'ecriture

      if(par%rank==0)then
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'River_ouest.out'                    &
              ,position='append')
         write(34,34)elapsedtime_now/86400.,                 &
                     (bil_fleuves_area(vb)/12.,vb=1,vbmax)        !mol on divise par le nombre de kount
         close(34)

         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'River_est.out'                    &
              ,position='append')
         write(34,34)elapsedtime_now/86400.,                 &
                     (bil_fleuves_area_est(vb)/12.,vb=1,vbmax)        !mol on divise par le nombre de kount
         close(34)

!        write(*,*)'bil_fleuves_area_est_nit FINAL = ',bil_fleuves_area_est(30)

! FIN: DOMAINE ENTIER
      endif                   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvv fin proc 0

      do vb=1,vbmax
       bil_fleuves_area(vb)=0.
       bil_fleuves_area_est(vb)=0.
      enddo

      endif !iteration3d éééééééééééééééééééééééééééééééééééééééééééééééééééééééééééé


 34      format(1(F13.5,1X),40(E,1X))


      RETURN
      END

