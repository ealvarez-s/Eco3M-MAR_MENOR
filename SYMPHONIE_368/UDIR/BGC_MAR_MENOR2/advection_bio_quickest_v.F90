!.................................................................

      subroutine advection_bio_quickest_v
      use module_principal
      use module_parallele
      implicit none
!     integer :: id_biobef_=0   & ! identifiant variable bio "before"
      integer ::  loopmaxbio_v_ &
                 ,loop_


!..............................................
! PARTIE ITERATIVE:
!..............................................
      if(checkxyt(1)/=xy_t(imax/2,jmax/2,1)) & ! securite debug
      stop 'checkxyt(1)/=xy_t(imax/2,jmax/2,1) advection_bio_quickest_v'

! Calculer anyv3d(:,:,:,id_bnegdif) le champ filtrE des petites longueurs d'ondes verticales:
!     if(ratio_bionegdif>0.)call advection_bio_negdif_fields !21-07-21

      do j=1,jmax ; do i=1,imax ! boucles i j no 3

      loopmaxbio_v_=ceiling(xy_t(i,j,1))

!     if(i==imax/2.and.j==jmax/2) then
!      write(6,*)'nombre de courant ',xy_t(i,j,1)
!      write(6,*)'loopmaxbio_v_',loopmaxbio_v_
!     endif

       do loop_=1,loopmaxbio_v_ !iterative loop>

! Advection partielle direction Ok:
!            do k=3,kmax-2 ! boucle k no 3
!            do k=kmin_w(i,j)+2,kmax-2 ! boucle k no 3
             do k=kmerged_t(i,j)+3,kmax-3 ! boucle k no 3

                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k+3,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k-3,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17

             enddo       ! boucle k no 3

! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
! Cas particulier:
                  k=min(kmerged_t(i,j)+2,kmax-2)
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k+3,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
! Cas particulier:
                  k=min(kmerged_t(i,j)+1,kmax-2)
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k+3,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
! Cas particulier:
                  k=min(kmerged_t(i,j)  ,kmax-2)
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k+3,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
! Cas particulier:
                  k=kmax-2
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k-3,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
! Cas particulier:
                  k=kmax-1
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k-3,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
! Cas particulier:
                  k=kmax  
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,19) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,18) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k-3,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17

! Verification de la conservation entre avant et apres:
!                 if(par%rank==85.and.i==imax/2.and.j==jmax/2) then
!                  sum0=0.
!                  sum1=0.
!                  sum2=0.
!                  do k=kmerged_t(i,j),kmax ! anciennement 1,kmax ! boucle k no 5
!                   sum0=sum0+   (0.5*(dz_t(i,j,k,0)   &
!                                     +dz_t(i,j,k,1)))  
!                   sum1=sum1+   (0.5*(dz_t(i,j,k,0)   &
!                                     +dz_t(i,j,k,1))) &
!                     *anyv3d(i  ,j  ,k  ,id_biobefo)
!                   sum2=sum2+   (0.5*(dz_t(i,j,k,0)   &
!                                     +dz_t(i,j,k,1))) &
!                     *( anyv3d(i  ,j  ,k  ,id_bioaftr)  &
!                       -anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12) ) ! Correction de divergence dw/dz !18-01-17
!                  enddo
!                  write(par%rank+100,*)sum1,sum2
!                  write(par%rank+200,*)anyv3d(i,j,kmerged_t(i,j),id_webio),anyv3d(i,j,kmax+1,id_webio)
!                 endif

! Cumul du terme de correction de divergence: (A oter A la toute fin de l'advection)
!            do k=kmin_w(i,j),kmax ! anciennement 1,kmax ! boucle k no 4
             do k=kmerged_t(i,j),kmax ! anciennement 1,kmax ! boucle k no 4
!                   bio_t(i  ,j  ,k,vb)= &
!                   bio_t(i  ,j  ,k,vb)  &
                   anyv3d(i  ,j  ,k,id_bdiv)= & !21-07-21
                   anyv3d(i  ,j  ,k,id_bdiv)  &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12) & ! Correction de divergence dw/dz !18-01-17
                  *(0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))) !01-05-19
             enddo       ! boucle k no 4

!            do k=kmin_w(i,j),kmax ! anciennement 1,kmax ! boucle k no 5
             do k=kmerged_t(i,j),kmax ! anciennement 1,kmax ! boucle k no 5

               anyv3d(i,j,k,id_biobefo)=anyv3d(i,j,k,id_bioaftr)    ! &

! Contrer les valeurs negatives !10-02-17 (inactif si flag_rmnegval=0):
!               +( anyv1d(k+1,1)*(anyv1d(k+1,2)-anyv1d(k,2))    &
!                 +anyv1d(k  ,1)*(anyv1d(k-1,2)-anyv1d(k,2)) )/(dz_t(i,j,k,0)+dz_t(i,j,k,1))

             enddo       ! boucle k no 5

       enddo                 !iterative loop>

! Phase de finalistion expliquee dans:
! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit
!      do k=kmin_w(i,j),kmax ! anciennement 1,kmax ! boucle k no 6
!      do k=kmerged_t(i,j),kmax ! anciennement 1,kmax ! boucle k no 6
       do k=1,kmax ! boucle k no 6

!        bio_t(i,j,k,vb)=                                  &
!       anyv3d(i,j,k,id_biobefo)                           &
!       -bio_t(i,j,k,vb) ! A la toute fin oter le cumul du terme de correction de divergence !18-01-17

! Attention ici bio_t est homogene A bio_t*dz_t, par consequent dans l'etape suivante dans mixsed_bio
! bio_t n'est plus multipliE par dz_t:
         bio_t(i,j,k,vb)=                                           &
        anyv3d(i,j,k,id_biobefo)*0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))  & !01-05-19
!       -bio_t(i,j,k,vb) ! A la toute fin oter le cumul du terme de correction de divergence !18-01-17
       -anyv3d(i,j,k,id_bdiv)    !21-07-21


       enddo       ! boucle k no 6

      enddo       ; enddo ! boucles i j no 3

! Petit test de conservation mono-proc:
!     sum1=0.
!     sum2=0.
!     sum3=0.
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      sum1=sum1+(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dxdy_t(i,j)*mask_t(i,j,kmax)
!      sum2=sum2+(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dxdy_t(i,j)*mask_t(i,j,kmax)*bio_t(i,j,k,1)
!      sum3=sum3+(dz_t(i,j,k,2)+dz_t(i,j,k,1))*dxdy_t(i,j)*mask_t(i,j,kmax)*bio_t(i,j,k,2)
!     enddo       ; enddo       ; enddo
!     write(6,*)'sum2/sum1=',sum2/sum1
!     write(6,*)'sum3/sum1=',sum3/sum1

!     stop 'coco4'

      end subroutine advection_bio_quickest_v

!.................................................................
