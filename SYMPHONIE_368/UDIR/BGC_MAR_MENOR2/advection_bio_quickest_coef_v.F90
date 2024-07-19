!.................................................................

      subroutine advection_bio_quickest_coef_v
      use module_principal
      use module_parallele
      implicit none
      real :: alpha_,beta_
!     integer :: id_webio=20 ! identifiant vitesse omega+wsed explicite

       beta_=   ratio_bionegdif*0.25
      alpha_=1.-ratio_bionegdif*0.5

!..............................................
! Determiner le nombre de sous-iterations d'advection
! Celui ci sera: int(xy_t(i,j,1))+1
!     substep_advbio est lavVitesse (en nbre de courant) verticale limite pour schema explicite
      xy_t(:,:,1)=small2             ! et non pas zero pour eviter loopmaxbio=0 si courant nul
      do j=1,jmax ; do i=1,imax

! Couche merged (dans laquelle w est 100% implicite, donc w explicite = ! 0)
       do k=1,kmerged_t(i,j) ! kmerged=kmin ou kmerged=kmin+1 si kmin et kmin+1 sont fusionnees
          anyv3d(i,j,k,id_webio)=0.
       enddo
!...........
! Dans la colonne d'eau:
       do k=kmerged_t(i,j)+1,kmax

          x0=0.5*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1)            &
                    ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))           &
                    /dti_fw

! omega+wsed explicite bornee
          anyv3d(i,j,k,id_webio)=max(min(                           &
         omega_w(i,j,k,1)+wsed(k,vb)*wsed_explicit(vb)            &
                              , substep_advbio*x0)                &
                              ,-substep_advbio*x0)*wetmask_t(i,j) !13-11-16

! bidouille patrick
!         anyv3d(i,j,k,id_webio)=0.

             xy_t(i,j,1)                                   &
        =max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_webio))/x0)

       enddo ! boucle k

!...........
! commentees le 07-03-19
! Cas de la surface (dz en k-1...)
!      k=kmax+1
!         x0=0.5*min(dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1)            &
!                   ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))           &
!                   /dti_fw

!         anyv3d(i,j,k,id_webio)=max(min(                           &
!        omega_w(i,j,k,1)+wsed(k,vb)*wsed_explicit(vb)            &
!                             , substep_advbio*x0)                &
!                             ,-substep_advbio*x0)*wetmask_t(i,j)

!            xy_t(i,j,1)                                   &
!       =min(max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_webio))/x0),substep_advbio-small2)


! Flux de surface imposEs dans fluxbio. Flux "advectif nuls" => vitesse omega explicite et implicite toutes 2 nulles
          anyv3d(i,j,kmax+1,id_webio)=0. !07-03-19
! details dans https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit


      enddo ; enddo ! Boucle i,j

! BIDOUILLE PATRICK

      checkxyt(1)=xy_t(imax/2,jmax/2,1)

! AU FOND: - omega_w=0
!          - la vitesse de sedimention EXPLICITE est nulle
! Il ne reste donc plus qu'a etablir anyv3d(i,j,k,id_webio)=0:
!     do j=1,jmax ; do i=1,imax
!      do k=1,kmin_w(i,j)-1
!         anyv3d(i,j,k,id_webio)=0.
!      enddo
!     enddo ; enddo ! Boucle i,j


!      stop 'jiji'

!..............................................
      do j=1,jmax ; do i=1,imax ! boucles i j no 2

       dti_fwsubio=dti_fw/ceiling(xy_t(i,j,1))


       do k=kmerged_t(i,j)+1,kmax !anciennement 2,kmax

! Nombre de courant vis a vis de omega
       x0=min(1.,   abs(anyv3d(i,j,k,id_webio))*dti_fwsubio      &
                 /(0.5*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1)     &
                          ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))))  &
                 *upwindriver_t(i,j)                           & !x0 standard si upwindriver=1
             +(1.-upwindriver_t(i,j))                            !x0=1 si upwindriver=0   !18-01-17

! Vitesses positivEes negativEes:
       x1=0.5*(anyv3d(i,j,k,id_webio)+abs(anyv3d(i,j,k,id_webio)))
       x2=0.5*(anyv3d(i,j,k,id_webio)-abs(anyv3d(i,j,k,id_webio)))

!...........................
! En facteur de B(k):
        anyv1d(k,3)=x1*(0.5-0.5*x0-(1.-x0**2)/6.) &
                   +x2*(0.5+0.5*x0+(1.-x0**2)/3.)

! En facteur de B(k-1):
        anyv1d(k,2)=x1*(0.5+0.5*x0+(1.-x0**2)/3.) &
                   +x2*(0.5-0.5*x0-(1.-x0**2)/6.)

! En facteur de B(k-2):
        anyv1d(k,1)=-x1*(1.-x0**2)/6.

! En facteur de B(k+1):
        anyv1d(k,4)=-x2*(1.-x0**2)/6.

       enddo         ! boucle k no 1

! En k=1 et k=kmax+1 (si omega/=0) schEma upwind
! anciennement do k=1,kmax+1,kmax
! details dans https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
!      k=kmax+1
!       anyv1d(k,3)=0. ! En facteur de B(k)  !23-02-17
!       anyv1d(k,2)=0. ! En facteur de B(k-1)!23-02-17
!       anyv1d(k,1)=0. ! En facteur de B(k-2)
!       anyv1d(k,4)=0. ! En facteur de B(k+1)
        anyv1d(kmax+1,1:4)=0.
      do k=1,kmerged_t(i,j) ! kmin_w(i,j)
!       anyv1d(k,3)=0. ! En facteur de B(k)
!       anyv1d(k,2)=0. ! En facteur de B(k-1)
!       anyv1d(k,1)=0. ! En facteur de B(k-2)
!       anyv1d(k,4)=0. ! En facteur de B(k+1)
        anyv1d(k,1:4)=0.
      enddo


! RIVIERES: Ne pas ecraser la valeur imposee au point derriere la source (oU mask_t=0) !01-03-19
! en multipliant les coefs en facteur de bio_t(i,j,k,vb) par 1-mask et les coefs en facteurs
! des points voisins par mask

!      do k=kmin_w(i,j),kmax ! anciennement 1,kmax  ! boucle k no 2
       do k=kmerged_t(i,j),kmax

       x3=dti_fwsubio/(0.5*(dz_t(i,j,k,0)   &
                           +dz_t(i,j,k,1)))

! Constante:
!      anyv3d(i,j,k,12)=x3*(omega_w(i,j,k+1,1)-omega_w(i,j,k,1))
       anyv3d(i,j,k,12)=x3*(anyv3d(i,j,k+1,id_webio)-anyv3d(i,j,k,id_webio)) &
                          *mask_t(i,j,kmax) !01-03-19


!...........................
! En facteur de B(k+3):
        anyv3d(i,j,k,19)=   x3*mask_t(i,j,kmax)*( & !ooo>
                    -(                                  &
                                  -anyv1d(k+1,4))*beta_  & 
                                                ) & !ooo>

! En facteur de B(k+2):
        anyv3d(i,j,k,18)=   x3*mask_t(i,j,kmax)*( & !ooo>
                    -(                                  &
                      +anyv1d(k,4)-anyv1d(k+1,3))*beta_  & 
                     +(-anyv1d(k+1,4)           )*alpha_ & 
                                                ) & !ooo>

! En facteur de B(k+1):
        anyv3d(i,j,k,17)=   x3*mask_t(i,j,kmax)*( & !ooo>
                    -(            -anyv1d(k+1,4)        &
                      +anyv1d(k,3)-anyv1d(k+1,2))*beta_  & 
                     +(anyv1d(k,4)-anyv1d(k+1,3))*alpha_ & 
                                                ) & !ooo>

! En facteur de B(k  ):
        anyv3d(i,j,k,16)=1.+x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k,4)-anyv1d(k+1,3)        &
                      +anyv1d(k,2)-anyv1d(k+1,1))*beta_  & 
                     +(anyv1d(k,3)-anyv1d(k+1,2))*alpha_ & 
                                                ) & !ooo>
       
! En facteur de B(k-1):
        anyv3d(i,j,k,15)=  +x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k,3)-anyv1d(k+1,2)        &
                      +anyv1d(k,1)              )*beta_  & 
                     +(anyv1d(k,2)-anyv1d(k+1,1))*alpha_ & 
                                                ) & !ooo>

! En facteur de B(k-2):
        anyv3d(i,j,k,14)=  +x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k,2)-anyv1d(k+1,1)        &
                                                )*beta_  & 
                     +(anyv1d(k,1)              )*alpha_ & 
                                                ) & !ooo>
                        
! En facteur de B(k-3):
        anyv3d(i,j,k,13)=  +x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k, 1)                     &
                                                )*beta_  & 
                                                ) & !ooo>

       enddo         ! boucle k no 2

      enddo       ; enddo ! boucles i j no 2

! Optionnel: nombre de courant maximum
      x1=0.
      do j=1,jmax ; do i=1,imax
       x1=max(x1,xy_t(i,j,1))
      enddo       ; enddo
      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d ,ierr)
      if(par%rank==0)then
       open(unit=3,file='tmp/dti_bio_w',position='append')
          write(3,'(i3,1x,6(e14.7,1x))')vb,elapsedtime_now/86400,wsed(1,vb),x2 !25-01-17
       close(3)
      endif


      end subroutine advection_bio_quickest_coef_v

!.................................................................
