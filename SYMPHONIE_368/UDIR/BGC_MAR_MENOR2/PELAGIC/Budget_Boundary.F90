      SUBROUTINE BUDGET_BOUNDARY

!_____________________________________________________________________*
!                                                                     *
! LAST REVISION: 24 JANUARY 2008                                      *
!                                                                     *
!_____________________________________________________________________*


!_____________________________________________________________________*
!                                                                     *
! Calculates the exchanges at the lateral boundaries                  *
!                                                                     *
!_____________________________________________________________________*

!_____________________________________________________________________*
!                                                                     *
! Modifications:                                                      *
! 04-07-10 : New version                                              *
!_____________________________________________________________________*



!---------------------------------------------------------------------*
! Declarations:


! Global variables

      use module_principal
      use ModuleDeclaration
      use module_parallele !#MPI
      use module_global
      implicit none
   

!     double precision x3,x4 

!---------------------------------------------------------------------*


      const1=dti_fw*0.5
      const2=86400.*0.5
!
! Western boundary: Convention flux >0 vers l'Est, entrant dans le domaine
!
      i1=iobwest          ! domaine global
      i=i1-par%timax(1)   ! domaine du proc
      if(i>=mbio1.and.i<=mbio2)then     !iiiiiiiiiiiiiiiiiiiiiiiiiiii

      do vb=1,vbmax
      do j=nbio1,nbio2  ! repere proc  
      j2=j+par%tjmax(1)  ! repere global
      if (j2>=jobsouth.and.j2<=jobnorth) then      !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

!     print*,'lat',lat_t(i,j)*180/pi

      x3=mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)
      x4=mask_i_w(i-1)*mask_j_w(j)*mask_t(i-1,j,k)

       do k=kmin_u(i,j),kmax

!      if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

       x1=-veldydz_u(i,j,k,1)* ( & !mask_u(i,j,k)*  (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )


       sum_ob_west_tot(vb)=sum_ob_west_tot(vb)+           &
         x1*const1

!       endif !mask 

       enddo !k


      if (-depth_w(i,j,kmin_w(i,j))>=1500) then 
  

       do k=kmin_u(i,j),kmax

       
       x1=-veldydz_u(i,j,k,1)*mask_u(i,j,k)*   (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )


!       print*,'x1',x1 

       flux_ob_west_vert(k,vb)=flux_ob_west_vert(k,vb)+   &
         x1*const2

       sum_ob_west_vert(k,vb)=sum_ob_west_vert(k,vb)+     &
         x1*const1

        if(depth_t(i,j,k)>-200)then
       flux_ob_west(vb)=flux_ob_west(vb)+                 &
         x1*const2

       sum_ob_west(vb)=sum_ob_west(vb)+                   &
         x1*const1

        else

       flux_ob_westdeep(vb)=flux_ob_westdeep(vb)+                 &
         x1*const2

       sum_ob_westdeep(vb)=sum_ob_westdeep(vb)+                   &
         x1*const1

        endif

      enddo !k

      elseif(-depth_w(i,j,kmin_w(i,j))<1500.and.lat_t(i,j)*180/pi<40.5)then

       do k=kmin_u(i,j),kmax

       x1=-veldydz_u(i,j,k,1)*mask_u(i,j,k)*   (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )


        if(depth_t(i,j,k)>=-200) then
       flux_ob_southwestshelf(vb)=flux_ob_southwestshelf(vb)+          &
         x1*const2

       sum_ob_southwestshelf(vb)=sum_ob_southwestshelf(vb)+            &
         x1*const1

        else

       flux_ob_southwestdeepshelf(vb)=flux_ob_southwestdeepshelf(vb)+  &
         x1*const2

       sum_ob_southwestdeepshelf(vb)=sum_ob_southwestdeepshelf(vb)+    &
         x1*const1

        endif !depth200

       enddo !k

      else ! depth

       do k=kmin_u(i,j),kmax


       x1=-veldydz_u(i,j,k,1)*mask_u(i,j,k)*   (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )


        if(depth_t(i,j,k)>=-200) then
       flux_ob_northwestshelf(vb)=flux_ob_northwestshelf(vb)+          &
         x1*const2

       sum_ob_northwestshelf(vb)=sum_ob_northwestshelf(vb)+            &
         x1*const1

        else

       flux_ob_northwestdeepshelf(vb)=flux_ob_northwestdeepshelf(vb)+  &
         x1*const2

       sum_ob_northwestdeepshelf(vb)=sum_ob_northwestdeepshelf(vb)+    &
         x1*const1


       endif                                !--------------------------->
      enddo !k


      endif ! depth
  
      endif !j2                         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      enddo !j
      enddo !vb 
      endif !i                          !iiiiiiiiiiiiiiiiiiiiiiiiiiii
      

! Eastern boundary: Convention >0 vers l'Est, sortant du domaine

      i1=iobeast          ! domaine global
      i=i1-par%timax(1)   ! domaine du proc
      if(i>=mbio1.and.i<=mbio2)then     !iiiiiiiiiiiiiiiiiiiiiiiiiiii

      do vb=1,vbmax
      do j=nbio1,nbio2  ! repere proc
      j2=j+par%tjmax(1)  ! repere global
      if (j2>=jobsouth.and.j2<=jobnorth) then      !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      x3=mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)
      x4=mask_i_w(i-1)*mask_j_w(j)*mask_t(i-1,j,k)

      do k=kmin_u(i,j),kmax

!       if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

       x1=-veldydz_u(i,j,k,1)* ( & !mask_u(i,j,k)*   (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )

       sum_ob_east_tot(vb)=sum_ob_east_tot(vb)+                   &
         x1*const1

!       endif !mask

       enddo

      if (-depth_w(i,j,kmin_w(i,j))>=1500) then 

      do k=kmin_u(i,j),kmax

!       if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

       x1=-veldydz_u(i,j,k,1)*                 (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )

       flux_ob_east_vert(k,vb)=flux_ob_east_vert(k,vb)+   &
         x1*const2

       sum_ob_east_vert(k,vb)=sum_ob_east_vert(k,vb)+     &
         x1*const1

        if(depth_t(i,j,k)>-200)then
       flux_ob_east(vb)=flux_ob_east(vb)+                 &
         x1*const2

       sum_ob_east(vb)=sum_ob_east(vb)+                   &
         x1*const1

        else

       flux_ob_eastdeep(vb)=flux_ob_eastdeep(vb)+                 &
         x1*const2

       sum_ob_eastdeep(vb)=sum_ob_eastdeep(vb)+                   &
         x1*const1

        endif

      enddo !k

      elseif (-depth_w(i,j,kmin_w(i,j))<1500.and.    &
          lat_t(i,j)*180/pi<43.5)then

      do k=kmin_u(i,j),kmax


       x1=-veldydz_u(i,j,k,1)*mask_u(i,j,k)*   (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )


        if(depth_t(i,j,k)>-200)then
       flux_ob_southeastshelf(vb)=flux_ob_southeastshelf(vb)+          &
         x1*const2

       sum_ob_southeastshelf(vb)=sum_ob_southeastshelf(vb)+            &
         x1*const1

        else

       flux_ob_southeastdeepshelf(vb)=flux_ob_southeastdeepshelf(vb)+  &
         x1*const2

       sum_ob_southeastdeepshelf(vb)=sum_ob_southeastdeepshelf(vb)+    &
         x1*const1

        endif

      enddo !k

      else  !depth

      do k=kmin_u(i,j),kmax

!       if(mask_t(i,j,kmax+1)==1) then       !--------------------------->

       x1=-veldydz_u(i,j,k,1)*mask_u(i,j,k)*   (          &
        sign(1.,veldydz_u(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i-1,j,k,vb)*x4 )


        if(depth_t(i,j,k)>-200)then
       flux_ob_northeastshelf(vb)=flux_ob_northeastshelf(vb)+          &
         x1*const2

       sum_ob_northeastshelf(vb)=sum_ob_northeastshelf(vb)+            &
         x1*const1

        else

       flux_ob_northeastdeepshelf(vb)=flux_ob_northeastdeepshelf(vb)+  &
         x1*const2

       sum_ob_northeastdeepshelf(vb)=sum_ob_northeastdeepshelf(vb)+    &
         x1*const1

        endif

      enddo !k

      endif !depth

      endif !j2                         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      enddo !j
      enddo !vb
      endif !i                          !iiiiiiiiiiiiiiiiiiiiiiiiiiii



! Southern boundary. Convention: >0 vers le Nord,  entrant dans le domaine

      j1=jobsouth
      j=j1-par%tjmax(1)   ! domaine du proc
      if(j>=nbio1.and.j<=nbio2)then     !iiiiiiiiiiiiiiiiiiiiiiiiiiii

      do vb=1,vbmax
      do i=mbio1,mbio2    !repere proc
      i2=i+par%timax(1)  ! repere global
      if (i2>=iobwest.and.i2<=iobeast) then      !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      x3=mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)
      x4=mask_i_w(i)*mask_j_w(j-1)*mask_t(i,j-1,k)

       do k=kmin_v(i,j),kmax


       x1=-veldxdz_v(i,j,k,1)* ( & !mask_v(i,j,k)*   (          &
        sign(1.,veldxdz_v(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4 )

       sum_ob_south_tot(vb)=sum_ob_south_tot(vb)+                   &
         x1*const1


       enddo !k


      if (-depth_w(i,j,kmin_w(i,j))>=1500)then   !commenté le 12/04/12 pour prendre tous les points cotiers au sud

       do k=kmin_v(i,j),kmax

!       if(mask_t(i,j,kmax+1)==1) then     !--------------------------->

       x1=-veldxdz_v(i,j,k,1)*mask_v(i,j,k)*   (          &
        sign(1.,veldxdz_v(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4 )

       flux_ob_south_vert(k,vb)=flux_ob_south_vert(k,vb)+   &
         x1*const2

       sum_ob_south_vert(k,vb)=sum_ob_south_vert(k,vb)+     &
         x1*const1

        if(depth_t(i,j,k)>-200)then
       flux_ob_south(vb)=flux_ob_south(vb)+                 &
         x1*const2

       sum_ob_south(vb)=sum_ob_south(vb)+                   &
         x1*const1

       else

       flux_ob_southdeep(vb)=flux_ob_southdeep(vb)+                 &
         x1*const2

       sum_ob_southdeep(vb)=sum_ob_southdeep(vb)+                   &
         x1*const1

        endif

!       endif                                !--------------------------->
      enddo !k

      else !depth

       do k=kmin_v(i,j),kmax


       x1=-veldxdz_v(i,j,k,1)*mask_v(i,j,k)*   (          &
        sign(1.,veldxdz_v(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4 )

        if(depth_t(i,j,k)>-200)then
       flux_ob_southshelf(vb)=flux_ob_southshelf(vb)+                 &
         x1*const2

       sum_ob_southshelf(vb)=sum_ob_southshelf(vb)+                   &
         x1*const1

       else

       flux_ob_southdeepshelf(vb)=flux_ob_southdeepshelf(vb)+          &
         x1*const2

       sum_ob_southdeepshelf(vb)=sum_ob_southdeepshelf(vb)+            &
         x1*const1

        endif

      enddo !k


      endif !depth

      endif !i2                         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      enddo !i
      enddo !vb
      endif !j                            !iiiiiiiiiiiiiiiiiiiiiiiiiiii



! Northern boundary. Convention: >0 vers le Nord,  sortant du domaine

      j1=jobnorth
      j=j1-par%tjmax(1)   ! domaine du proc
      if(j>=nbio1.and.j<=nbio2)then     !iiiiiiiiiiiiiiiiiiiiiiiiiiii

      do vb=1,vbmax
      do i=mbio1,mbio2    !repere proc
      i2=i+par%timax(1)  ! repere global
      if (i2>=iobwest.and.i2<=iobeast) then      !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      x3=mask_i_w(i)*mask_j_w(j)*mask_t(i,j,k)
      x4=mask_i_w(i)*mask_j_w(j-1)*mask_t(i,j-1,k)

      if (-depth_w(i,j,kmin_w(i,j))>=1500)then

       do k=kmin_v(i,j),kmax


       x1=-veldxdz_v(i,j,k,1)*mask_v(i,j,k)*   (          &
        sign(1.,veldxdz_v(i,j,k,1))                       &
           *(bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4)           &
            -bio_t(i,j,k,vb)*x3-bio_t(i,j-1,k,vb)*x4 )

       flux_ob_north_vert(k,vb)=flux_ob_north_vert(k,vb)+   &
         x1*const2

       sum_ob_north_vert(k,vb)=sum_ob_north_vert(k,vb)+     &
         x1*const1

        if(depth_t(i,j,k)>-200)then
       flux_ob_north(vb)=flux_ob_north(vb)+                 &
         x1*const2

       sum_ob_north(vb)=sum_ob_north(vb)+                   &
         x1*const1

        endif

!       endif                                !--------------------------->
      enddo !k

      endif !depth

      endif !i2                         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      enddo !i
      enddo !vb
      endif !j                            !iiiiiiiiiiiiiiiiiiiiiiiiiiii


      if(mod(kount,12)==0)then !èèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèèè


#ifdef parallele
      do vb=1,vbmax
      call mpi_allreduce(sum_ob_south(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_south(vb)=sum1glb
      call mpi_allreduce(sum_ob_southdeep(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southdeep(vb)=sum1glb
      call mpi_allreduce(sum_ob_southshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_southdeepshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southdeepshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_east(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_east(vb)=sum1glb
      call mpi_allreduce(sum_ob_eastdeep(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_eastdeep(vb)=sum1glb
      call mpi_allreduce(sum_ob_southeastshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southeastshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_southeastdeepshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southeastdeepshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_northeastshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_northeastshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_northeastdeepshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_northeastdeepshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_west(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_west(vb)=sum1glb
      call mpi_allreduce(sum_ob_westdeep(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_westdeep(vb)=sum1glb
      call mpi_allreduce(sum_ob_southwestshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southwestshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_southwestdeepshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_southwestdeepshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_northwestshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_northwestshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_northwestdeepshelf(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_northwestdeepshelf(vb)=sum1glb
      call mpi_allreduce(sum_ob_west_tot(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_west_tot(vb)=sum1glb
      call mpi_allreduce(sum_ob_east_tot(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_east_tot(vb)=sum1glb
      call mpi_allreduce(sum_ob_south_tot(vb),sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum_ob_south_tot(vb)=sum1glb
      enddo
#endif

! seul un proc travaille pour les dernieres operations sur les quantités globales et l'ecriture

      

      if(par%rank==0)then                    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'frontierecol_shelfarea.out'               &
              ,position='append')
         write(34,34)float(kount)*dti_fw/86400.,                 &
                (-sum_ob_southwestshelf(vb)/12.  &
                -sum_ob_northwestshelf(vb)/12.  &
                +sum_ob_southeastshelf(vb)/12.  &
                +sum_ob_northeastshelf(vb)/12.  &
                -sum_ob_southshelf(vb)/12.  &
                ,vb=1,vbmax),                        & 
               (-sum_ob_southwestdeepshelf(vb)/12.  &
                -sum_ob_northwestdeepshelf(vb)/12.  &
                +sum_ob_southeastdeepshelf(vb)/12.  &
                +sum_ob_northeastdeepshelf(vb)/12.  &
                -sum_ob_southdeepshelf(vb)/12.  &
                ,vb=1,vbmax)
         close(34)


         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'frontierecol_deeparea.out'               &
              ,position='append')
         write(34,34)float(kount)*dti_fw/86400.,                 &
                 (-sum_ob_west(vb)/12.  &
                  +sum_ob_east(vb)/12.  &
                  -sum_ob_south(vb)/12.  &
                ,vb=1,vbmax),                        &
               (-sum_ob_westdeep(vb)/12.  &
                 +sum_ob_eastdeep(vb)/12.  &
                 -sum_ob_southdeep(vb)/12.  &
                 ,vb=1,vbmax)
         close(34)


         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'frontierecol_south.out'               &
              ,position='append')
         write(34,34)float(kount)*dti_fw/86400.,                 &
                (sum_ob_southshelf(vb)/12.     &
                +sum_ob_southdeepshelf(vb)/12. &
                ,vb=1,vbmax),              & 
                (sum_ob_south(vb)/12.  &
                +sum_ob_southdeep(vb)/12.  &
                ,vb=1,vbmax)
         close(34)

         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'frontierecol_east.out'               &
              ,position='append')
         write(34,34)float(kount)*dti_fw/86400.,                 &
                (sum_ob_southeastshelf(vb)/12.  &
                +sum_ob_northeastshelf(vb)/12.  &
                +sum_ob_southeastdeepshelf(vb)/12.  &
                +sum_ob_northeastdeepshelf(vb)/12.  &
                ,vb=1,vbmax),                   &
                (sum_ob_east(vb)/12.  &
                +sum_ob_eastdeep(vb)/12.  &
                ,vb=1,vbmax)
         close(34)
  

         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'frontierecol_west.out'               &
              ,position='append')
         write(34,34)float(kount)*dti_fw/86400.,                 &
                (sum_ob_southwestshelf(vb)/12.  &
                +sum_ob_northwestshelf(vb)/12.  &
                +sum_ob_southwestdeepshelf(vb)/12.  &
                +sum_ob_northwestdeepshelf(vb)/12.  &
                ,vb=1,vbmax),                   &
                (sum_ob_west(vb)/12.  &
                +sum_ob_westdeep(vb)/12.  &
                ,vb=1,vbmax)
         close(34)


         open(unit=34,file=                                      &
              dirgraph(1:lname4)//'frontierecol_wes.out'               &
              ,position='append')
         write(34,35)float(kount)*dti_fw/86400.,                 &
                (sum_ob_west_tot(vb)/12.,vb=1,vbmax),  &
                (sum_ob_east_tot(vb)/12.,vb=1,vbmax),  &
                (sum_ob_south_tot(vb)/12.,vb=1,vbmax)
         close(34)


! FIN: DOMAINE ENTIER
      endif                   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvv fin proc 0


! Remise a zero

        do vb=1,vbmax
           sum_ob_southwestshelf(vb)=0.
           sum_ob_northwestshelf(vb)=0.
           sum_ob_southeastshelf(vb)=0.
           sum_ob_northeastshelf(vb)=0.
           sum_ob_southshelf(vb)=0.
           sum_ob_southwestdeepshelf(vb)=0.
           sum_ob_northwestdeepshelf(vb)=0.
           sum_ob_southeastdeepshelf(vb)=0.
           sum_ob_northeastdeepshelf(vb)=0.
           sum_ob_southdeepshelf(vb)=0.
           sum_ob_west(vb)=0.
           sum_ob_east(vb)=0.
           sum_ob_south(vb)=0.
           sum_ob_westdeep(vb)=0.
           sum_ob_eastdeep(vb)=0.
           sum_ob_southdeep(vb)=0.
           sum_ob_west_tot(vb)=0.
           sum_ob_east_tot(vb)=0.
           sum_ob_south_tot(vb)=0.
        enddo

         endif !kount ééééééééééééééééééééééééééééééééééééééééééééééééééééééé


 34      format(1(F13.5,1X),66(E,1X))
 35      format(1(F13.5,1X),99(E,1X))



      return
      end
