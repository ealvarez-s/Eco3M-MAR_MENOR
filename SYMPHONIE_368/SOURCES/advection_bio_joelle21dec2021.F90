      subroutine advection_bio
!______________________________________________________________________
! SYMPHONIE ocean model
! release 307 - last update: 27-08-21
!______________________________________________________________________
!    _________                    .__                  .__     m[°v°]m !
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!......................................................................
      use module_principal
      use module_parallele
      implicit none
!..............................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         14/12/01: bienvenue à ITIMEBIO
!         17/12/01: amenagements pour cas forward (economie de RAM)
!         26/12/02: amenagements pour remise en service schema forward
!         14/02/03: variable xy_z(i,j,1) était non initialisée. Ajout boucle 104
!                   de mise à zéro.
!         19/03/03: date de debut et de fin de calcul pour traceur
!         04/05/06: amenagements pour compiler en double precision
!         04/07/06: boucle 103: j avant i. Regroupement de constantes.
!         20/04/07: passage à la coordonnées curviligne
! 2009.3  01-10-09: suppression kount_bio
!         02-10-09: dti remplacé par dti_fw
! 2010.8  03-05-10  supression biohz et nouveau schema forward
! 2010.10 21-06-10  advection pour traceur avec sous-pas de temps
! 2010.22 01-05-11  Reecriture du cas sous-pas-de-temps advection tracer
! 2010.24 23-11-11  Ajouter une verification optionnelle (non active par defaut)
!                   du critere cfl pour l'advection
! S.26    06-10-13  float replaced by real
!         22-01-14  commentaire devant test sur masque
!         06-02-14  obc parallelisation
!         29-04-14  nloop devient substep_advbio defini dans namelist "notebook_advection"
!         16-07-14  Ajout des schemas advection_bioup2up3 advection_bio_tvd_pmx
!         20-07-14  Suite point precedent. Aspects C.L. fond, surface, fleuves
!         25-03-15  ibl1_advbio,ibl2_advbio,jbl1_advbio,jbl2_advbio definisse
!                   une zone buffer sans advection pour la bio
!         19-07-15  - utilisation de upwindriver_t au schema bioup2up3 pour schema
!                   100% up2 dans le voisinage des embouchures
!                   - bioup2up3: schema up2 dans zone intertidales
!         05-08-15  schema vst implique utilisation du masque de surface
!         17-03-16  verifier les dimensions A l'etat initial et aprEs lecture
!                   des fichiers restart
!         25-04-16  pas d'advection bio si flag_1dv=1
!         03-07-16  Ajout d'une subroutine calculant le nombre de sous-iterations advectives
!         04-07-16  appel A routine obc mpi deplacE pour conservation mpi
!         11-07-16  ajout routine advection_bio_rmnegval supprimant les valeurs negatives
!         12-11-16  Ajout schema QUICKEST
!         13-11-16  QUICKEST : deplacement wetmask
!         06-12-16  Suppression d'un return dans subroutine rmnegval et ajout d'un commentaire
!         13-12-16  Choix wsed implicite ou explicite selon les traceurs
!         18-01-17  - Prise en compte de upwindriver_t dans l'advection (voir aussi notebook_advection)
!                   - Correction d'un bug
!                   - Meilleure definition du terme de compensation de divergence
!         25-01-17  Modif ecriture fichier dti_bio_w
!         08-02-17  ajout advection_bio_min_max
!         12-02-17  Mises A jour cas NEMO offline
!         13-02-17  Mix Max modulo 100
!         23-02-17  Flux air/mer portE par fluxbio_w et non pas par omega
!         03-11-17  ajout subroutine advection_bio_rmnegval3d_plus
!         24-03-18  implementation grille verticale fusionnee
! v247    01-03-19  ne pas modifier la valeur du point en amont du point source riviere
! v249    07-03-19  - suite du point 23-02-17: meme si ca ne change rien il semble
!                   plus clair de prescrire vitesse nulle en kmax+1
!                   - vertmix_merged_levels_2t_any(1,id_bioaftr) tient compte du
!                   fait que le volume de controle de bio_t s'etablit sur 2 pas
!                   de temps
! v253    01-05-19  Evitement de la division par dz
!         06-05-19  fichiers fort.xxx si loopmaxbio>50
! v285    05-06-20  suite point precedent: suppresion dsig_u et v
! v287    14-08-20  obc_int_anyv3d renomme obc_mpi_anyv3d
! v303    21-07-21  - anyv3d(i,j,k,id_bdiv) en double precision plutot que bio_t temporaire real*4
!                   pour plus de precision sur l'advection
!                   - diffusion negatide
!                   - bilan bio
!         31-07-21  Advection sous forme flux
!         04-08-21  anyv3d(:,:,:,id_bnegdif)=0.      !04-08-21
! v304    07-08-21  pas de dif neg si vitesse de chute !07-08-21
!         11-08-21  if(vb<=vbmax_eco3ms.and.wsed(1,vb)>wsed(1,vb-1))call advection_bio_erreur_stop !11-08-21
! v307    27-08-21  id_bioaft_ devient id_bioaftr
!                   cond. limite anyv3d(i,j,k,id_bnegdif)=anyv3d(i,j,kmax-1,id_bnegdif) 
!..............................................................................
#ifdef synopsis
       subroutinetitle='advection_bio'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       if(flag_1dv==1)return !25-04-16 pas d'advection bio si flag_1dv=1

! Decommenter cette ligne pour diag min max sur bio_t
!      call advection_bio_min_max !08-02-17

! Nombre de sous-pas de temps advectifs
       if(iadvec_bio==3) then !pmxpmx> !12-11-16

! Verification des dimensions pour anyv3d
       if(iteration3d==0.or.iteration3d==iteration3d_restart) then !>>>
           call advection_quickest_checkdim
           call advection_checkdim_anyv1d
       endif                                                       !>>>

! advection horizontale pretraitement:
        call advection_bio_substep
        call advection_bio_quickest_coef_h_flux !31-07-21

      do vb=1,vbmax ! vbvbvbvb>

!     if(ratio_bionegdif>0. &
!       .and.wsed(1,vb)==0. & ! pas de dif neg si vitesse de chute !07-08-21
!                          ) then !>>>  
!       call advection_bio_negdif_fields !21-07-21
!     else                        !>>>
!       anyv3d(:,:,:,id_bnegdif)=0.      !04-08-21
!     endif                       !>>>

! advection horizontale:
        call advection_bio_quickest_h_flux !31-07-21

! advection verticale:
        k0=0
        if(vb==1)k0=1
        if(vb>1) then !ooo>
         if(wsed_explicit(vb)==1.and.wsed(1,vb)/=wsed(1,vb-1))k0=1 !13-12-16
!        if(vb<=vbmax_eco3ms.and.wsed(1,vb)>wsed(1,vb-1))call advection_bio_erreur_stop !11-08-21
        endif         !ooo>
        if(k0==1)call advection_bio_quickest_coef_v

        call advection_bio_quickest_v

      enddo         ! vbvbvbvb>

      return

       else                   !pmxpmx>

       stop 'Dans notebook_advection choisir iadvec_bio=3'

       call advection_bio_get_substep !03-07-16

       if(iadvec_bio==0) then !000000>
                              call advection_bioup
                              return
       endif                  !000000>

!      if(iadvec_bio==2) then !111111>
!                             call advection_bioup2up3
!                             return
!      endif                  !111111>

       endif                  !pmxpmx>

      stop ' No scheme to compute bio advection'
      end subroutine advection_bio

!..............................................................................

      subroutine advection_bioup
      use module_principal
      use module_parallele
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='advection_bioup'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...................................................................!
! NOTE SUR LE SCHEMA D'ADVECTION                                    !
! Schema upwind
!...................................................................!

!*********************************************************************
! ADVECTION VARIABLES DU MODELE BIOLOGIQUE ou SEDIMENTAIRE
! DEBUT:
!*********************************************************************
!     substep_advbio est defini dans namelist "notebook_advection"
      const1=0.5*dti_fw/real(substep_advbio) !06-10-13

      do j=1,jmax ; do i=1,imax
       xy_t(i,j,0)=1. ! Sert a annuler omega(kmax+1) si omega(kmax+1)>0
      enddo ; enddo

      do  k=1,kmax    ! Obligatoirement croissant a cause de xy_t(:,:,0)

       if(k==kmax) then ! -- Surface case --> (note: useless if omega(kmax+1)==0)
        do j=1,jmax ; do i=1,imax
         if(omega_w(i,j,kmax+1,1)>0.)xy_t(i,j,0)=0. ! No upward flux across surface
        enddo       ; enddo
       endif            ! -- Surface case -->

      do  j=1,jmax
      do  i=1,imax

      x1=   sign(un,veldydz_u(i+1,j,k,1))
      x2=   sign(un,veldydz_u(i,  j,k,1))
      x3=   sign(un,veldxdz_v(i,j+1,k,1))
      x4=   sign(un,veldxdz_v(i,j  ,k,1))
      x5=   sign(un,omega_w(i,j,k+1,1))
      x6=   sign(un,omega_w(i,j,k  ,1))

! Coef en facteur de bio_z(i,j,k):
      anyv3d(i,j,k,0)=                                               &
       (  veldydz_u(i+1,j  ,k  ,1)*const1*(-1.-x1)                   &
         +veldydz_u(i  ,j  ,k  ,1)*const1*( 1.-x2)                   &
         +veldxdz_v(i  ,j+1,k  ,1)*const1*(-1.-x3)                   &
         +veldxdz_v(i  ,j  ,k  ,1)*const1*( 1.-x4)  )/dxdy_t(i,j)    &
           +omega_w(i  ,j  ,k+1,1)*const1*(-1.-x5)*xy_t(i,j,0)       &
           +omega_w(i  ,j  ,k  ,1)*const1*( 1.-x6)

! Coef en facteur de bio_t(i+1,j,k):
      anyv3d(i,j,k,1)=veldydz_u(i+1,j  ,k  ,1)*const1*(-1.+x1)/dxdy_t(i,j)

! Coef en facteur de bio_t(i-1,j,k):
      anyv3d(i,j,k,2)=veldydz_u(i  ,j  ,k  ,1)*const1*( 1.+x2)/dxdy_t(i,j)

! Coef en facteur de bio_t(i,j+1,k):
      anyv3d(i,j,k,3)=veldxdz_v(i  ,j+1,k  ,1)*const1*(-1.+x3)/dxdy_t(i,j)

! Coef en facteur de bio_t(i,j-1,k):
      anyv3d(i,j,k,4)=veldxdz_v(i  ,j  ,k  ,1)*const1*( 1.+x4)/dxdy_t(i,j)

! Coef en facteur de bio_t(i,j,k+1):
      anyv3d(i,j,k,5)=omega_w(i  ,j  ,k+1,1)*const1*(-1.+x5)*xy_t(i,j,0)

! Coef en facteur de bio_t(i,j,k-1):
      anyv3d(i,j,k,6)=omega_w(i  ,j  ,k  ,1)*const1*( 1.+x6)

      enddo
      enddo
      enddo

      do loop_=1,substep_advbio

      call obc_bio_botsurf   ! obligatoirement avant obc_bio(1)
      call obc_bio(1)        ! rivieres
      call obc_bio_mpi('z1') !04-07-16

      rap2=real(loop_  )/real(substep_advbio)          !01-05-11
      rap1=real(loop_-1)/real(substep_advbio)          !01-05-11

      x1=(1.-rap1)*0.5
      x2=    rap1 *0.5
      x3=(1.-rap2)*0.5
      x4=    rap2 *0.5
! Uncomment this line to check the CFL criteria:
!     call advection_bioup_check_cfl    ! 23-11-11

! Note: si substep_advbio=1 rap2=1 rap1=0 x4=0.5 x3=0 x2=0 x1=0.5 -> dz(t+1/2)=0.5*(dz(t+1)+dz(t))
!                                                           dz(t-1/2)=0.5*(dz(t)+dz(t-1))

      do vb=1,vbmax
!     do vb=1,1
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      anyv3d(i,j,k,7)=(bio_t(i  ,j  ,k  ,vb)                       &
                     *( x2*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))  &   !01-05-11
                       +x1*(dz_t(i,j,k,now  )+dz_t(i,j,k,before))) &

           +anyv3d(i,j,k,0)*bio_t(i  ,j  ,k  ,vb)                  &
           +anyv3d(i,j,k,1)*bio_t(i+1,j  ,k  ,vb)                  &
           +anyv3d(i,j,k,2)*bio_t(i-1,j  ,k  ,vb)                  &
           +anyv3d(i,j,k,3)*bio_t(i  ,j+1,k  ,vb)                  &
           +anyv3d(i,j,k,4)*bio_t(i  ,j-1,k  ,vb)                  &
           +anyv3d(i,j,k,5)*bio_t(i  ,j  ,k+1,vb)                  &
           +anyv3d(i,j,k,6)*bio_t(i  ,j  ,k-1,vb)                  &
                   )/( x4*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))   &  !01-05-11
                      +x3*(dz_t(i,j,k,now  )+dz_t(i,j,k,before)))

      enddo    !i
      enddo    !j
      enddo    !k

!     call advection_bio_test_pmx1

! update bio_t (move anyv3d on bio_t) and bottom boundary condition
       call advection_biofinalize(7) !25-03-15

      enddo    !vb

!     call obc_bio_mpi('z1') !06-02-14

      enddo    !loop_

      return

!*********************************************************************
! ADVECTION VARIABLES DU MODELE BIOLOGIQUE ou SEDIMENTAIRE
! FIN
!*********************************************************************

      end subroutine advection_bioup

!.....................................................................

      subroutine advection_bioup_check_cfl    ! 23-11-11
      use module_principal
      implicit none
      double precision cfl_extrem1_,cfl_extrem2_,x0_
#ifdef synopsis
       subroutinetitle='advection_bioup_check_cfl'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     cfl_extrem1_=0.   ; cfl_extrem2_=1.   ! Tolerance zero
      cfl_extrem1_=-0.2 ; cfl_extrem2_=1.2  ! Tolerance + ou - 20% de la cfl stricte

      do k=1,kmax
      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmax)==1) then !111111111111111>

      x0_=1./( x4*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))   &  !01-05-11
                 +x3*(dz_t(i,j,k,now  )+dz_t(i,j,k,before)))

       do k0=1,6
        if(    anyv3d(i,j,k,k0)*x0_ > cfl_extrem2_        &
           .or.anyv3d(i,j,k,k0)*x0_ < cfl_extrem1_ ) then


        write(*,*)'CFL criteria not respected in routine adection_bioup'
        write(*,*)'Grid point i j k=',i,j,k
        stop ' STOP dans subroutine advection_bioup_check_cfl'
        endif

       enddo ! loop_ k0

      endif                     !111111111111111>
      enddo  ! loop_ i
      enddo  ! loop_ j
      enddo  ! loop_ k


      end subroutine advection_bioup_check_cfl

!.....................................................................

      subroutine advection_bioup3_checkdim
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='advection_bioup3_checkdim'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Verifier dimensions 2DH de bio_t
       ub4=ubound(bio_t) ; lb4=lbound(bio_t)
       if(ub4(2)<jmax+2)stop 'Erreur max dimension 2 bio_t'
       if(lb4(2)>-1)    stop 'Erreur min dimension 2 bio_t'
       if(ub4(1)<imax+2)stop 'Erreur max dimension 1 bio_t'
       if(lb4(1)>-1)    stop 'Erreur min dimension 1 bio_t'

! Verifier 4eme dimension de anyv3d
       ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
!      if(lb4(4)>0 .or.ub4(4)<13.or. &
!         lb4(3)>-1.or.ub4(3)<kmax+2) then !>>>>>>
!       deallocate(anyv3d)
!         allocate(anyv3d(-1:imax+2,-1:jmax+2,-1:kmax+2,0:13)) ; anyv3d=0
!      endif                               !>>>>>>
       if(lb4(4)>0 .or.ub4(4)<13) then !>>>>>>
        deallocate(anyv3d)
          allocate(anyv3d(-1:imax+2,-1:jmax+2,0:kmax+1,0:13)) ; anyv3d=0
       endif                           !>>>>>>

      end subroutine advection_bioup3_checkdim

!.....................................................................
#ifdef bidon
#endif
!.....................................................................
#ifdef bidon
#endif
!.....................................................................
#ifdef bidon
#endif
!________________________________________________________________________________
#ifdef bidon
#endif
!.....................................................................
#ifdef bidon
      subroutine advection_bio_test_pmx
      use module_principal
      use module_parallele
      implicit none
      sum1=0.
      sum3=0.
      do k=1,kmax
       kp2=min(k+2,kmax+1) ! zero z-gradient if z>0
       kp1=min(k+1,kmax+1) ! zero z-gradient if z>0
       km2=max(k-2,1)      ! zero z-gradient if z<-h
       km1=max(k-1,1)      ! zero z-gradient if z<-h
       do j=1,jmax
       do i=1,imax
        sum1=sum1+(                                                   &
        anyv3d(i,j,k,13)*( x4*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))  &
                          +x3*(dz_t(i,j,k,now  )+dz_t(i,j,k,before))) &
                        *dxdy_t(i,j)                                  &
       -bio_t(i,j,k,vb) *( x2*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))  &
                          +x1*(dz_t(i,j,k,now  )+dz_t(i,j,k,before))) &
                        *dxdy_t(i,j)   )                              &
                 *mask_i_w(i)*mask_j_w(j)                             &
                 *mask_t(i,j,kmax)
        sum3=sum3+(                                                   &
            -anyv3d(i,j,k,0 )*bio_t(i+2,j  ,k  ,vb)                  &
            -anyv3d(i,j,k,1 )*bio_t(i+1,j  ,k  ,vb)                  &
            -anyv3d(i,j,k,2 )*bio_t(i  ,j  ,k  ,vb)                  &
            -anyv3d(i,j,k,3 )*bio_t(i-1,j  ,k  ,vb)                  &
            -anyv3d(i,j,k,4 )*bio_t(i-2,j  ,k  ,vb)                  &
            -anyv3d(i,j,k,5 )*bio_t(i  ,j+2,k  ,vb)                  &
            -anyv3d(i,j,k,6 )*bio_t(i  ,j+1,k  ,vb)                  &
            -anyv3d(i,j,k,7 )*bio_t(i  ,j-1,k  ,vb)                  &
            -anyv3d(i,j,k,8 )*bio_t(i  ,j-2,k  ,vb)                  &
            -anyv3d(i,j,k,9 )*bio_t(i  ,j  ,kp2,vb)                  &
            -anyv3d(i,j,k,10)*bio_t(i  ,j  ,kp1,vb)                  &
            -anyv3d(i,j,k,11)*bio_t(i  ,j  ,km1,vb)                  &
            -anyv3d(i,j,k,12)*bio_t(i  ,j  ,km2,vb) )                &
                        *dxdy_t(i,j)                                 &
                 *mask_i_w(i)*mask_j_w(j)                            &
                 *mask_t(i,j,kmax)
#ifdef synopsis
       subroutinetitle='advection_bio_test_pmx'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!       if(mask_t(i,j,kmax)==1)write(6,*)mask_t(i,j,kmax),bio_t(i  ,j  ,k  ,vb)

       enddo
       enddo
      enddo
      sum2=0.
      k=kmax+1
      do j=1,jmax
      do i=1,imax
       if(omega_w(i,j,k,1)>0.) then
        k1=kmax
       else
        k1=kmax+1
       endif
       sum2=                                         &
       sum2-min(omega_w(i,j,k,1),zero)               &
            *dxdy_t(i,j)*dti_fw/real(substep_advbio) &
                    *bio_t(i,j,k1,vb)                &
                         *mask_t(i,j,kmax)              &
                       *mask_i_w(i)                  &
                       *mask_j_w(j)
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      if(par%rank==0) then
       write(6,*)'------'
       write(6,*)'sum1  ',sum1glb
       write(6,*)'sum2  ',sum2glb
       write(6,*)'sum3  ',sum3glb
      endif
      end subroutine advection_bio_test_pmx
#endif
!.....................................................................
#ifdef bidon
      subroutine advection_bio_test_pmx1
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='advection_bio_test_pmx1'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      sum1=0.
      do k=1,kmax
       do j=1,jmax
       do i=1,imax
        sum1=sum1+(                                                   &
        anyv3d(i,j,k, 7)*( x4*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))  &
                          +x3*(dz_t(i,j,k,now  )+dz_t(i,j,k,before))) &
                        *dxdy_t(i,j)                                  &
       -bio_t(i,j,k,vb) *( x2*(dz_t(i,j,k,after)+dz_t(i,j,k,now   ))  &
                          +x1*(dz_t(i,j,k,now  )+dz_t(i,j,k,before))) &
                        *dxdy_t(i,j)   )                              &
                 *mask_i_w(i)*mask_j_w(j)                             &
                 *mask_t(i,j,kmax)
       enddo
       enddo
      enddo
      sum2=0.
      k=kmax+1
      do j=1,jmax
      do i=1,imax
       if(omega_w(i,j,k,1)>0.) then
        k1=kmax
       else
        k1=kmax+1
       endif
       sum2=                                         &
       sum2-min(omega_w(i,j,k,1),zero)               &
            *dxdy_t(i,j)*dti_fw/real(substep_advbio) &
                    *bio_t(i,j,k1,vb)                &
                         *mask_t(i,j,kmax)              &
                       *mask_i_w(i)                  &
                       *mask_j_w(j)
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI
                        mpi_sum,par%comm2d,ierr)
      if(par%rank==0) then
       write(6,*)'------'
       write(6,*)'sum1  ',sum1glb
       write(6,*)'sum2  ',sum2glb
      endif
      end subroutine advection_bio_test_pmx1
#endif
!________________________________________________________________________________
      subroutine advection_biofinalize(ind4_)
      use module_principal ; use module_parallele
      implicit none
      integer ind4_
#ifdef synopsis
       subroutinetitle='advection_biofinalize'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do j=1,jmax
       if(j+par%tjmax(1)>=   1+jbl1_advbio.and. &
          j+par%tjmax(1)<=jglb-jbl2_advbio) then !jjjjj>

        do i=1,imax
         if(i+par%timax(1)>=   1+ibl1_advbio.and. &
            i+par%timax(1)<=iglb-ibl2_advbio) then !iiiii>

          do k=kmin_w(i,j),kmax
           bio_t(i,j,k,vb)=anyv3d(i,j,k,ind4_)
          enddo
          do k=1,kmin_w(i,j)-1
           bio_t(i,j,k,vb)=bio_t(i,j,kmin_w(i,j),vb)
          enddo

         endif                                     !iiiii>
        enddo

       endif                                      !jjjjj>
      enddo

      end subroutine advection_biofinalize

!.........................................................................

      subroutine advection_bio_get_substep !03-07-16
      use module_principal ; use module_parallele
      implicit none


      substep_advbio=1 ! Default value

      x1=0.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!#ifdef bidon
       x1=max(x1, &

                 wetmask_t(i,j)*dti_fw*(                               &

        ((veldydz_u(i+1,j,k,1)+abs(veldydz_u(i+1,j,k,1)))              &
        -(veldydz_u(i  ,j,k,1)-abs(veldydz_u(i  ,j,k,1)))              &
        +(veldxdz_v(i,j+1,k,1)+abs(veldxdz_v(i,j+1,k,1)))              &
        -(veldxdz_v(i,j  ,k,1)-abs(veldxdz_v(i,j  ,k,1))))/dxdy_t(i,j) &
          +(omega_w(i,j,k+1,1)+abs(  omega_w(i,j,k+1,1)))              &
          -(omega_w(i,j,k  ,1)-abs(  omega_w(i,j,k  ,1)))              &

                                       )/(2.*dz_t(i,j,k,2))            &

              )
!#endif
#ifdef bidon
       x2=       wetmask_t(i,j)*dti_fw*(                               &

        ((veldydz_u(i+1,j,k,1)+abs(veldydz_u(i+1,j,k,1)))              &
        -(veldydz_u(i  ,j,k,1)-abs(veldydz_u(i  ,j,k,1)))              &
        +(veldxdz_v(i,j+1,k,1)+abs(veldxdz_v(i,j+1,k,1)))              &
        -(veldxdz_v(i,j  ,k,1)-abs(veldxdz_v(i,j  ,k,1))))/dxdy_t(i,j) &
          +(omega_w(i,j,k+1,1)+abs(  omega_w(i,j,k+1,1)))              &
          -(omega_w(i,j,k  ,1)-abs(  omega_w(i,j,k  ,1)))              &

                                       )/(2.*dz_t(i,j,k,2))

        if(x2>x1) then
         x1=x2 ; i1=i ; j1=j ; k1=k
        endif
#endif




!      if(mask_t(i,j,kmax)==1)write(20+par%rank,*)i,j,k &

!               ,wetmask_t(i,j)*dti_fw*(                               &

!       ((veldydz_u(i+1,j,k,1)+abs(veldydz_u(i+1,j,k,1)))              &
!       -(veldydz_u(i  ,j,k,1)-abs(veldydz_u(i  ,j,k,1)))              &
!       +(veldxdz_v(i,j+1,k,1)+abs(veldxdz_v(i,j+1,k,1)))              &
!       -(veldxdz_v(i,j  ,k,1)-abs(veldxdz_v(i,j  ,k,1))))/dxdy_t(i,j) &
!         +(omega_w(i,j,k+1,1)+abs(  omega_w(i,j,k+1,1)))              &
!         -(omega_w(i,j,k  ,1)-abs(  omega_w(i,j,k  ,1)))              &

!                                      )/(2.*dz_t(i,j,k,2))


      enddo ; enddo ; enddo

!     i=i1 ; j=j1 ; k=k1
!     write(20+par%rank,*)x1,i1,j1,k1                          &
!     ,veldydz_u(i,j,k,1)/dy_u(i,j)/(0.5*(dz_t(i  ,j,k,1)      &
!                                        +dz_t(i-1,j,k,1)))    &
!     ,veldxdz_v(i,j,k,1)/dx_v(i,j)/(0.5*(dz_t(i,j  ,k,1)      &
!                                        +dz_t(i,j-1,k,1)))    &
!     ,omega_w(i,j,k,1)                                        &
!     ,dy_u(i,j),dx_v(i,j),dz_t(i  ,j,k,1)


      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d ,ierr)
      if(x2>1.)substep_advbio=int(x2)+1

      if(par%rank==0) then !>>>>>
       open(unit=3,file='tmp/substep_advbio',position='append')
       write(3,*)elapsedtime_now/86400.,substep_advbio
       close(3)
      endif                !>>>>>
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!#endif
!      stop 'bibi'



      end subroutine advection_bio_get_substep

!.........................................................................

      subroutine advection_bio_rmnegval !11-07-16
      use module_principal ; use module_parallele
      implicit none

      stop 'UTILISER advection_bio_rmnegval3d_plus'

! Détails et evolutions possibles dans: !06-12-16
! https://docs.google.com/document/d/15QsKJhPFiOjxU3os7omfhrhhxNYZBsRTprBmsyuNguc/edit

! Remove Negative Values

      call obc_bio_mpi('zb')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculer un bilan pour verification des proprietes de conservation
!     sum1=0.
!     sum2=0.
!     vb=1
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      sum1=sum1+dz_t(i,j,k,2)*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)
!      sum2=sum2+dz_t(i,j,k,2)*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)*bio_t(i,j,k,vb)
!     enddo
!     enddo
!     enddo
!     call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     if(par%rank==0)write(66,*)iteration3d,sum2glb/sum1glb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
        anyv3d(i,j,k,1)=0.25*min(dz_t(i,j,k,2),dz_t(i-1,j,k,2))*mask_u(i,j,kmax)
      enddo
      enddo
      enddo
! Ne pas diffuser A travers les C.L. que nous ne maitrisons pas suffisament:
      if(obcstatus(ieq1)==1)   anyv3d(1     ,:,:,1)=0.
      if(obcstatus(ieqimax)==1)anyv3d(imax+1,:,:,1)=0.


      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
        anyv3d(i,j,k,2)=0.25*min(dz_t(i,j,k,2),dz_t(i,j-1,k,2))*mask_v(i,j,kmax)
      enddo
      enddo
      enddo
! Ne pas diffuser A travers les C.L. que nous ne maitrisons pas suffisament:
      if(obcstatus(jeq1)==1)   anyv3d(:,1     ,:,2)=0.
      if(obcstatus(jeqjmax)==1)anyv3d(:,jmax+1,:,2)=0.

      do vb=1,vbmax

      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       anyv3d(i,j,k,0)=min(bio_t(i,j,k,vb),0.)
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

        bio_t(i,j,k,vb)=   &
        bio_t(i,j,k,vb)+(  & !ooooo>

        anyv3d(i+1,j,k,1)*(anyv3d(i+1,j,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i  ,j,k,1)*(anyv3d(i-1,j,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i,j+1,k,2)*(anyv3d(i,j+1,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i,j  ,k,2)*(anyv3d(i,j-1,k,0)-anyv3d(i,j,k,0)) &

                        )  & !ooooo>
                         /dz_t(i,j,k,2)


      enddo
      enddo
      enddo

      enddo ! vb

!     if(iteration3d==10) then
!     j=jmax/2 ; vb=1 ; k=kmax
!     do i=1,imax
!      write(20+par%rank,*)i+par%timax(1),bio_t(i,j,k,vb)
!     enddo
!     call graph_out
!     stop 'coucou'
!     endif


      end subroutine advection_bio_rmnegval

!.........................................................................

      subroutine advection_bio_substep
      use module_principal
      use module_parallele
      implicit none

!..............................................
! Determiner le nombre de sous-iterations d'advection
! Cet algo suppose que la couche merged est 100% melangee avec les
! couches sous-jacente
! et donc que l'on a affaire A une couche dont l'epaisseur peut etre
! consideree comme
! z_w(kmerged+1)-z_w(1) et non pas dz(kmerged) (methode 2). Si on garde
! dz(kmerged) (methode 1) l'algo considere
! un volume plus petit pour plus de securite (mais evntuellement un dt
! inutilement trop petit)
! D'autre part les flux qui entrent dans ce volume sont ceux du niveau
! k=kmerged mais egalement
! ceux sous kmerged, kundermin_t etant le niveau du premier flux lateral
! non nul
      x1=small2 ! et non pas zero pour eviter loopmaxbio=0
                ! si courant nul
! Niveaux standard:
      do j=1,jmax ; do i=1,imax
        do k=kmerged_t(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(veldydz_u(i  ,j  ,k,1))   &
                              ,abs(veldydz_u(i+1,j  ,k,1)))  &
                              ,abs(veldxdz_v(i  ,j  ,k,1)))  &
                              ,abs(veldxdz_v(i  ,j+1,k,1)))  &
                    *dti_fw                                  &
                    *wetmask_t(i,j)                          &
                       /(dxdy_t(i,j)*                        &
                         0.5*( dz_t(i,j,k,0)+dz_t(i,j,k,1))) &
               ) ! max
        enddo
      enddo       ; enddo
! Niveaux <= kmerged
      do j=1,jmax ; do i=1,imax
        sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
        do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
         sum1=sum1+veldydz_u(i  ,j  ,k,1)
         sum2=sum2+veldydz_u(i+1,j  ,k,1)
         sum3=sum3+veldxdz_v(i  ,j  ,k,1)
         sum4=sum4+veldxdz_v(i  ,j+1,k,1)
         sum5=sum5     +dz_t(i  ,j  ,k,0)
         sum6=sum6     +dz_t(i  ,j  ,k,2)
        enddo
         x1=max(x1,max(max(max(abs(sum1)  &
                              ,abs(sum2)) &
                              ,abs(sum3)) &
                              ,abs(sum4)) &
                    *dti_fw                                  &
                    *wetmask_t(i,j)                          &
                      /(dxdy_t(i,j)*min(sum5,sum6))          &
! methode 2
               ) ! fermeture max
      enddo       ; enddo


      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     loopmaxbio=int(x2)+1
      loopmaxbio=ceiling(x2)
      dti_fwsubio=dti_fw/loopmaxbio

      if(par%rank==0)then
       open(unit=3,file='tmp/dti_bio',position='append')
        write(3,*)real(elapsedtime_now/86400.),loopmaxbio,real(x2),real(dti_fwsubio)
       close(3)
      endif

!     write(6,*)'loopmaxbio horizontal',loopmaxbio
      if(loopmaxbio>50) then !m°v°m> !06-05-19
      do j=1,jmax ; do i=1,imax
        do k=kmerged_t(i,j)+1,kmax
         x2=       max(max(max(abs(veldydz_u(i  ,j  ,k,1))   &
                              ,abs(veldydz_u(i+1,j  ,k,1)))  &
                              ,abs(veldxdz_v(i  ,j  ,k,1)))  &
                              ,abs(veldxdz_v(i  ,j+1,k,1)))  &
                    *dti_fw                                  &
                    *wetmask_t(i,j)                          &
                       /(dxdy_t(i,j)*                        &
                         0.5*( dz_t(i,j,k,0)+dz_t(i,j,k,1)))  

             if(int(x2+1)>50) then !(0v0)>
               write(10+par%rank,*)'---------------'
               write(10+par%rank,*)'x2=            ',x2
               write(10+par%rank,*)'i,j,k loc      ',i,j,k
               write(10+par%rank,*)'i,j,k glb      ',i+par%timax(1),j+par%tjmax(1),k
               write(10+par%rank,*)'mask_t         ',mask_t(i,j,k)
               write(10+par%rank,*)'h_w            ',h_w(i,j)
               write(10+par%rank,*)'ssh_int_w(0:1) ',ssh_int_w(i,j,0:1)
               write(10+par%rank,*)'dz_t(i,j,k,0:1)',dz_t(i,j,k,0:1)
               write(10+par%rank,*)'dsig_t         ',dsig_t(i,j,k)
               write(10+par%rank,*)'kmin_w kmerge  ',kmin_w(i,j),kmerged_t(i,j)
               write(10+par%rank,*)'wetmask_t(i,j) ',wetmask_t(i,j)
               write(10+par%rank,*)'veldydz_u(i  ,j  ,k,1)',veldydz_u(i  ,j  ,k,1)
               write(10+par%rank,*)'veldydz_u(i+1,j  ,k,1)',veldydz_u(i+1,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j  ,k,1)',veldxdz_v(i ,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j+1,k,1)',veldxdz_v(i ,j+1,k,1)
               write(10+par%rank,*)'---------------'
             endif                 !(0v0)>

        enddo
      enddo       ; enddo
! Niveaux <= kmerged
      do j=1,jmax ; do i=1,imax
        sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
        do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
         sum1=sum1+veldydz_u(i  ,j  ,k,1)
         sum2=sum2+veldydz_u(i+1,j  ,k,1)
         sum3=sum3+veldxdz_v(i  ,j  ,k,1)
         sum4=sum4+veldxdz_v(i  ,j+1,k,1)
         sum5=sum5     +dz_t(i  ,j  ,k,0)
         sum6=sum6     +dz_t(i  ,j  ,k,2)
        enddo
         x2=       max(max(max(abs(sum1)  &
                              ,abs(sum2)) &
                              ,abs(sum3)) &
                              ,abs(sum4)) &
                    *dti_fw                                  &
                    *wetmask_t(i,j)                          &
                      /(dxdy_t(i,j)*min(sum5,sum6))           

             if(int(x2+1)>50) then !(0v0)>
               write(10+par%rank,*)'---------------'
               write(10+par%rank,*)'x2=            ',x2
               write(10+par%rank,*)'i,j,k loc      ',i,j,k
               write(10+par%rank,*)'i,j,k glb      ',i+par%timax(1),j+par%tjmax(1),k
               write(10+par%rank,*)'mask_t         ',mask_t(i,j,k)
               write(10+par%rank,*)'h_w            ',h_w(i,j)
               write(10+par%rank,*)'sum1           ',sum1
               write(10+par%rank,*)'sum2           ',sum1
               write(10+par%rank,*)'sum3           ',sum1
               write(10+par%rank,*)'sum4           ',sum1
               write(10+par%rank,*)'sum5           ',sum1
               write(10+par%rank,*)'sum6           ',sum1
               write(10+par%rank,*)'---------------'
             endif                 !(0v0)>

      enddo       ; enddo
      call graph_out
      stop 'loopmaxbio>50 see fortxxx error files'
      endif                  !m°v°m>

      end subroutine advection_bio_substep

!.................................................................

      subroutine advection_bio_quickest_coef_v
      use module_principal
      use module_parallele
      implicit none
!     integer :: id_webio=20 ! identifiant vitesse omega+wsed explicite


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

!...........................
! Section diffusion negative:
! En facteur de <B(k)>:
!       anyv1d(k,7)=-( x1*(   -0.5*x0*(1.-x0)-(1.-x0**2)/6.) &
!                     +x2*(   +0.5*x0*(1.-x0)+(1.-x0**2)/3.))

! En facteur de <B(k-1)>:
!       anyv1d(k,6)=-( x1*(   +0.5*x0*(1.-x0)+(1.-x0**2)/3.) &
!                     +x2*(   -0.5*x0*(1.-x0)-(1.-x0**2)/6.))

! En facteur de <B(k-2)>:
!       anyv1d(k,5)=x1*(1.-x0**2)/6.

! En facteur de <B(k+1)>:
!       anyv1d(k,8)=x2*(1.-x0**2)/6.

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

!...........................
! En facteur de B(k+2):
        anyv3d(i,j,k,15)=   x3*(           -anyv1d(k+1,4)) &
                              *mask_t(i,j,kmax) !01-03-19

! En facteur de B(k+1):
        anyv3d(i,j,k,14)=   x3*(anyv1d(k,4)-anyv1d(k+1,3)) &
                              *mask_t(i,j,kmax) !01-03-19

! En facteur de B(k  ):
        anyv3d(i,j,k,13)=1.+x3*(anyv1d(k,3)-anyv1d(k+1,2)) &
                              *mask_t(i,j,kmax) !01-03-19

! En facteur de B(k-1):
        anyv3d(i,j,k,16)=   x3*(anyv1d(k,2)-anyv1d(k+1,1)) &
                              *mask_t(i,j,kmax) !01-03-19

! En facteur de B(k-2):
        anyv3d(i,j,k,17)=   x3*(anyv1d(k,1)              ) &
                              *mask_t(i,j,kmax) !01-03-19

! Constante:
!      anyv3d(i,j,k,12)=x3*(omega_w(i,j,k+1,1)-omega_w(i,j,k,1))
       anyv3d(i,j,k,12)=x3*(anyv3d(i,j,k+1,id_webio)-anyv3d(i,j,k,id_webio)) &
                          *mask_t(i,j,kmax) !01-03-19

!...........................
! Section diffusion negative:
! En facteur de <B(k+2)>:
!       anyv3d(i,j,k,20)=   x3*(           -anyv1d(k+1,8)) &
!                             *mask_t(i,j,kmax) !01-03-19

! En facteur de <B(k+1)>:
!       anyv3d(i,j,k,19)=   x3*(anyv1d(k,8)-anyv1d(k+1,7)) &
!                             *mask_t(i,j,kmax) !01-03-19

! En facteur de <B(k  )>:
!       anyv3d(i,j,k,18)=  +x3*(anyv1d(k,7)-anyv1d(k+1,6)) &
!                             *mask_t(i,j,kmax) !01-03-19

! En facteur de <B(k-1)>:
!       anyv3d(i,j,k,21)=   x3*(anyv1d(k,6)-anyv1d(k+1,5)) &
!                             *mask_t(i,j,kmax) !01-03-19

! En facteur de <B(k-2)>:
!       anyv3d(i,j,k,22)=   x3*(anyv1d(k,5)              ) &
!                             *mask_t(i,j,kmax) !01-03-19


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
             do k=kmerged_t(i,j)+2,kmax-2 ! boucle k no 3

                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,18) &
!                 +anyv3d(i  ,j  ,k+1,id_bnegdif)*anyv3d(i,j,k,19) &
!                 +anyv3d(i  ,j  ,k+2,id_bnegdif)*anyv3d(i,j,k,20) &
!                 +anyv3d(i  ,j  ,k-1,id_bnegdif)*anyv3d(i,j,k,21) &
!                 +anyv3d(i  ,j  ,k-2,id_bnegdif)*anyv3d(i,j,k,22)  

             enddo       ! boucle k no 3

! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
! Cas particulier:
!                 k=2
!                 k=min(kmin_w(i,j)+1,kmax-1)
                  k=min(kmerged_t(i,j)+1,kmax-1)
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,18) &
!                 +anyv3d(i  ,j  ,k+1,id_bnegdif)*anyv3d(i,j,k,19) &
!                 +anyv3d(i  ,j  ,k+2,id_bnegdif)*anyv3d(i,j,k,20) &
!                 +anyv3d(i  ,j  ,k-1,id_bnegdif)*anyv3d(i,j,k,21) &
!                 +anyv3d(i  ,j  ,k-1,id_bnegdif)*anyv3d(i,j,k,22)  
! Cas particulier:
!                 k=1
!                 k=min(kmin_w(i,j),kmax-1)
                  k=min(kmerged_t(i,j),kmax-1)
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k+2,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,18) &
!                 +anyv3d(i  ,j  ,k+1,id_bnegdif)*anyv3d(i,j,k,19) &
!                 +anyv3d(i  ,j  ,k+2,id_bnegdif)*anyv3d(i,j,k,20) &
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,21) &
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,22)  

! Cas particulier:
                  k=kmax-1
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k+1,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,18) &
!                 +anyv3d(i  ,j  ,k+1,id_bnegdif)*anyv3d(i,j,k,19) &
!                 +anyv3d(i  ,j  ,k+1,id_bnegdif)*anyv3d(i,j,k,20) &
!                 +anyv3d(i  ,j  ,k-1,id_bnegdif)*anyv3d(i,j,k,21) &
!                 +anyv3d(i  ,j  ,k-2,id_bnegdif)*anyv3d(i,j,k,22)  
! Cas particulier:
                  k=kmax
                   anyv3d(i  ,j  ,k  ,id_bioaftr)=                 &
                   anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,13) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,14) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,15) &
                  +anyv3d(i  ,j  ,k-1,id_biobefo)*anyv3d(i,j,k,16) &
                  +anyv3d(i  ,j  ,k-2,id_biobefo)*anyv3d(i,j,k,17) &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*anyv3d(i,j,k,12)   ! Correction de divergence dw/dz !18-01-17
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,18) &
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,19) &
!                 +anyv3d(i  ,j  ,k  ,id_bnegdif)*anyv3d(i,j,k,20) &
!                 +anyv3d(i  ,j  ,k-1,id_bnegdif)*anyv3d(i,j,k,21) &
!                 +anyv3d(i  ,j  ,k-2,id_bnegdif)*anyv3d(i,j,k,22)  

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

      subroutine advection_checkdim_anyv1d
      use module_principal
      implicit none

! Verifier 2eme dimension de anyv1d
       ub2=ubound(anyv1d) ; lb2=lbound(anyv1d)
       if(lb2(2)>1.or.ub2(2)<4) then !>>>>>>
        deallocate(anyv1d)
          allocate(anyv1d(lb2(1):ub2(1),1:4)) ; anyv1d=0
       endif

      end subroutine advection_checkdim_anyv1d

!.................................................................

      subroutine advection_quickest_checkdim
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='advection_quickest_checkdim'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Verifier dimensions 2DH de bio_t
       ub4=ubound(bio_t) ; lb4=lbound(bio_t)
       if(ub4(2)<jmax+2)stop 'Erreur max dimension 2 bio_t'
       if(lb4(2)>-1)    stop 'Erreur min dimension 2 bio_t'
       if(ub4(1)<imax+2)stop 'Erreur max dimension 1 bio_t'
       if(lb4(1)>-1)    stop 'Erreur min dimension 1 bio_t'

! Verifier 4eme dimension de anyv3d
       ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
       if(lb4(4)>0 .or.ub4(4)<23) then !>>>>>> !21-07-21
        deallocate(anyv3d)
          allocate(anyv3d(-1:imax+2,-1:jmax+2,0:kmax+1,0:23)) ; anyv3d=0 !21-07-21
       endif                           !>>>>>>

      end subroutine advection_quickest_checkdim

!.....................................................................

      subroutine advection_bio_min_max !08-02-17
      use module_principal ; use module_parallele ; use module_s
      implicit none

      if(mod(iteration3d,100)/=0)return !13-02-17

! Calcule min et max et indices globaux correspondants des tableaux bio_t et archive ces
! info dans des fichiers tmp/biominVB et tmp/biomaxVB

      istr=2 ; iend=imax-1
      jstr=2 ; jend=jmax-1
      if(obcstatus(ieq1)==1)   istr=1
      if(obcstatus(ieqimax)==1)iend=imax
      if(obcstatus(jeq1)==1)   jstr=1
      if(obcstatus(jeqjmax)==1)jend=jmax

      do vb=1,vbmax

      x1= 1.e10
      x2=-1.e10
      do k=1,kmax ; do j=jstr,jend ; do i=istr,iend

      if(mask_t(i,j,kmax)==1) then !pmxpmx>

       if(bio_t(i,j,k,vb)>x2) then !maxmax>
        x2=bio_t(i,j,k,vb) ; i2=i ; j2=j ; k2=k
       endif                       !maxmax>

       if(bio_t(i,j,k,vb)<x1) then !minmin>
        x1=bio_t(i,j,k,vb) ; i1=i ; j1=j ; k1=k
       endif                       !minmin>

       if(isnan(bio_t(i,j,k,vb))) then !nan>
        write(6,*)'bio_t is NAN par%rank,i,j,k,vb',par%rank,i+par%timax(1),j+par%tjmax(1),k,vb
        stop 'bio_t is NAN'
       endif                           !nan>


      endif                     !pmxpmx>


      enddo ; enddo ; enddo

      call mpi_allreduce(x2,x22,1,mpi_double_precision,mpi_max,par%comm2d ,ierr)
      call mpi_allreduce(x1,x11,1,mpi_double_precision,mpi_min,par%comm2d ,ierr)

      if(x2==x22) then !max>
! Si cette condition est satisfaite on est dans le rank du maximum
       write(texte30,'(a,i0)')'tmp/biomax',vb
       open(unit=3,file=texte30,position='append')
        write(3,'(2(e14.7,1x),4i5)')elapsedtime_now/86400.,x22,par%rank &
                                   ,i2+par%timax(1),j2+par%tjmax(1),k2
       close(3)
      endif            !max>

       if(x1==x11.and.x11<0.) then !min>
! Si cette condition est satisfaite on est dans le rank du minimum
        write(texte30,'(a,i0)')'tmp/biomin',vb
        open(unit=3,file=texte30,position='append')
         write(3,'(2(e14.7,1x),4i5)')elapsedtime_now/86400.,x11,par%rank &
                                    ,i1+par%timax(1),j1+par%tjmax(1),k1
        close(3)
       endif                       !min>

      enddo ! vb loop
      call s_cpu('advection_bio_minmax',0)


      end subroutine advection_bio_min_max

!.........................................................................
      subroutine advection_bio_rmnegval3d !10-02-17
      use module_principal ; use module_parallele ; use module_s
      implicit none


! Détails et evolutions possibles dans: !06-12-16
! https://docs.google.com/document/d/15QsKJhPFiOjxU3os7omfhrhhxNYZBsRTprBmsyuNguc/edit

! Remove Negative Values

      call obc_bio_mpi('zb')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculer un bilan pour verification des proprietes de conservation
!     sum1=0.
!     sum2=0.
!     vb=1
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      sum1=sum1+dz_t(i,j,k,0)+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)
!      sum2=sum2+dz_t(i,j,k,0)+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)*bio_t(i,j,k,vb)
!     enddo
!     enddo
!     enddo
!     call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     if(par%rank==0)write(66,*)iteration3d,sum2glb/sum1glb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Coef diff Oi
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
        anyv3d(i,j,k,1)=(1./6.)*min(dz_t(i  ,j,k,0)+dz_t(i  ,j,k,1) &
                                   ,dz_t(i-1,j,k,0)+dz_t(i-1,j,k,1))*mask_u(i,j,kmax)
      enddo
      enddo
      enddo
! Ne pas diffuser A travers les C.L. que nous ne maitrisons pas suffisament:
      if(obcstatus(ieq1)==1)   anyv3d(1     ,:,:,1)=0.
      if(obcstatus(ieqimax)==1)anyv3d(imax+1,:,:,1)=0.


! Coef diff Oj
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
        anyv3d(i,j,k,2)=(1./6.)*min(dz_t(i,j  ,k,0)+dz_t(i,j  ,k,1)  &
                                   ,dz_t(i,j-1,k,0)+dz_t(i,j-1,k,1))*mask_v(i,j,kmax)
      enddo
      enddo
      enddo
! Ne pas diffuser A travers les C.L. que nous ne maitrisons pas suffisament:
      if(obcstatus(jeq1)==1)   anyv3d(:,1     ,:,2)=0.
      if(obcstatus(jeqjmax)==1)anyv3d(:,jmax+1,:,2)=0.

! Coef diff Ok
      do j=1,jmax
      do i=1,imax

       do k=kmin_w(i,j)+1,kmax
        anyv3d(i,j,k,3)=(1./6.)*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1)  &
                                   ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))
       enddo
       do k=1,kmin_w(i,j)
        anyv3d(i,j,k,3)=0.
       enddo
       anyv3d(i,j,kmax+1,3)=0.

      enddo
      enddo

! Diffusion

      do vb=1,vbmax

      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       anyv3d(i,j,k,0)=min(bio_t(i,j,k,vb),0.)
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

        bio_t(i,j,k,vb)=   &
        bio_t(i,j,k,vb)+(  & !ooooo>

        anyv3d(i+1,j,k,1)*(anyv3d(i+1,j,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i  ,j,k,1)*(anyv3d(i-1,j,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i,j+1,k,2)*(anyv3d(i,j+1,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i,j  ,k,2)*(anyv3d(i,j-1,k,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i,j,k+1,3)*(anyv3d(i,j,k+1,0)-anyv3d(i,j,k,0)) &
       +anyv3d(i,j,k  ,3)*(anyv3d(i,j,k-1,0)-anyv3d(i,j,k,0)) &

                        )  & !ooooo>
                         /(dz_t(i,j,k,0)+dz_t(i,j,k,1))


      enddo
      enddo
      enddo

      enddo ! vb
      call s_cpu('advection_bio_rmnegval3d',0)

      end subroutine advection_bio_rmnegval3d
!.........................................................................
      subroutine advection_bio_nemobottom(id_var_) !11-02-17
      use module_principal
      implicit none
      integer id_var_
       do j=1,jmax
       do i=1,imax
        do k=1,kmin_w(i,j)-1
         anyv3d(i,j,k,id_var_)=anyv3d(i,j,kmin_w(i,j),id_var_)
        enddo
       enddo
       enddo
      end subroutine advection_bio_nemobottom
!.........................................................................
      subroutine advection_bio_rmnegval3d_plus !03-11-17
      use module_principal ; use module_parallele ; use module_s
      implicit none

! Détails et evolutions possibles dans: !06-12-16
! https://docs.google.com/document/d/15QsKJhPFiOjxU3os7omfhrhhxNYZBsRTprBmsyuNguc/edit

! Remove Negative Values

      call obc_bio_mpi('zb')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculer un bilan pour verification des proprietes de conservation
!     sum1=0.
!     sum2=0.
!     vb=1
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      sum1=sum1+dz_t(i,j,k,0)+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)
!      sum2=sum2+dz_t(i,j,k,0)+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)*bio_t(i,j,k,vb)
!     enddo
!     enddo
!     enddo
!     call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     if(par%rank==0)write(66,*)iteration3d,sum2glb/sum1glb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Coef diff Oi
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
        anyv3d(i,j,k,1)=min(dz_t(i  ,j,k,0)+dz_t(i  ,j,k,1) &
                           ,dz_t(i-1,j,k,0)+dz_t(i-1,j,k,1))*mask_u(i,j,kmax)
      enddo
      enddo
      enddo
! Ne pas diffuser A travers les C.L. que nous ne maitrisons pas suffisament:
      if(obcstatus(ieq1)==1)   anyv3d(1     ,:,:,1)=0.
      if(obcstatus(ieqimax)==1)anyv3d(imax+1,:,:,1)=0.


! Coef diff Oj
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
        anyv3d(i,j,k,2)=min(dz_t(i,j  ,k,0)+dz_t(i,j  ,k,1)  &
                           ,dz_t(i,j-1,k,0)+dz_t(i,j-1,k,1))*mask_v(i,j,kmax)
      enddo
      enddo
      enddo
! Ne pas diffuser A travers les C.L. que nous ne maitrisons pas suffisament:
      if(obcstatus(jeq1)==1)   anyv3d(:,1     ,:,2)=0.
      if(obcstatus(jeqjmax)==1)anyv3d(:,jmax+1,:,2)=0.

! Coef diff Ok
      do j=1,jmax
      do i=1,imax

       do k=kmin_w(i,j)+1,kmax
        anyv3d(i,j,k,3)=min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1)  &
                           ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))
       enddo
       do k=1,kmin_w(i,j)
        anyv3d(i,j,k,3)=0.
       enddo
       anyv3d(i,j,kmax+1,3)=0.

      enddo
      enddo

! Diffusion

      do vb=1,vbmax

      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
!       anyv3d(i,j,k,0)=min(bio_t(i,j,k,vb),0.)
        anyv3d(i,j,k,0)=min(bio_t(i,j,k,vb),small1)
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

        bio_t(i,j,k,vb)=   &
        bio_t(i,j,k,vb)+(  & !ooooo>

        anyv3d(i+1,j,k,1)*(anyv3d(i+1,j,k,0)-anyv3d(i,j,k,0))   &
          *max(0.,-sign(1.d0,anyv3d(i+1,j,k,0)*anyv3d(i,j,k,0)))  & !=1 si l'un des deux seulement est <0. O sinon. !03-11-17
       +anyv3d(i  ,j,k,1)*(anyv3d(i-1,j,k,0)-anyv3d(i,j,k,0))   &
          *max(0.,-sign(1.d0,anyv3d(i-1,j,k,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j+1,k,2)*(anyv3d(i,j+1,k,0)-anyv3d(i,j,k,0))   &
          *max(0.,-sign(1.d0,anyv3d(i,j+1,k,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j  ,k,2)*(anyv3d(i,j-1,k,0)-anyv3d(i,j,k,0))   &
          *max(0.,-sign(1.d0,anyv3d(i,j-1,k,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j,k+1,3)*(anyv3d(i,j,k+1,0)-anyv3d(i,j,k,0))   &
          *max(0.,-sign(1.d0,anyv3d(i,j,k+1,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j,k  ,3)*(anyv3d(i,j,k-1,0)-anyv3d(i,j,k,0))   &
          *max(0.,-sign(1.d0,anyv3d(i,j,k-1,0)*anyv3d(i,j,k,0)))  &

                        )*0.16  & !ooooo> ! Note: 0.16=troncature inferieure de 1/6
                         /(dz_t(i,j,k,0)+dz_t(i,j,k,1))


      enddo
      enddo
      enddo

      enddo ! vb
      call s_cpu('advection_bio_rmnegval3d',0)

      end subroutine advection_bio_rmnegval3d_plus
!.........................................................................

      subroutine advection_bio_negdif_fields !21-07-21
      use module_principal ; use module_parallele
      implicit none
      integer loop_

!#ifdef bidon
! Filtre 3 points: !09-05-21
!      do k=2,kmax-1 ; do j=-1,jmax+2 ; do i=-1,imax+2
       do k=2,kmax-1 ; do j= 1,jmax   ; do i= 1,imax  
! Note: je prends le temps 1 et non pas le temps 0 car dans le cas de
! l'antidiffusion, il me semble que c'est le temps "1" qui va eviter
! l'amplification du mode numerique. On garde quand meme une EMA coef 0.5

        anyv3d(i,j,k,id_bnegdif)=ratio_bionegdif*( & !m°v°m>
!                               0.5*anyv3d(i,j,k  ,id_biobefo)   &
!                            +0.25*(anyv3d(i,j,k+1,id_biobefo)   &
!                                  +anyv3d(i,j,k-1,id_biobefo))  &
                                 0.5*bio_t(i,j,k  ,vb)   &
                              +0.25*(bio_t(i,j,k+1,vb)   &
                                    +bio_t(i,j,k-1,vb))  &
                                                 )   !m°v°m>
       enddo ; enddo ; enddo

! A priori, tel qu'est construit le shema d'advection vertical, ces conditions
! aux limites ne sont pas utiles, et donc elles sont commentees:
! CL surface et fond:
       do k=kmax,kmax+1
!       do j=-1,jmax+2 ; do i=-1,imax+2
        do j= 1,jmax   ; do i= 1,imax  
           anyv3d(i,j,k,id_bnegdif)=anyv3d(i,j,kmax-1,id_bnegdif) !27-08-21 
!        +(anyv3d(i,j,k,id_biobefo)-anyv3d(i,j,kmax-1,id_biobefo))*ratio_bionegdif !31-07-21
        enddo ; enddo
       enddo
       do k=0,1
!       do j=-1,jmax+2 ; do i=-1,imax+2
        do j= 1,jmax   ; do i= 1,imax  
           anyv3d(i,j,k,id_bnegdif)=anyv3d(i,j,2,id_bnegdif)   
!        +(anyv3d(i,j,k,id_biobefo)-anyv3d(i,j,2,id_biobefo))*ratio_bionegdif !31-07-21
        enddo ; enddo
       enddo
!#endif

#ifdef bidon
! Filtre 5 points:
       x0=1./16.
!      do k=3,kmax-2 ; do j=-1,jmax+2 ; do i=-1,imax+2
       do k=3,kmax-2 ; do j= 1,jmax   ; do i= 1,imax  
! Note: je prends le temps 1 et non pas le temps 0 car dans le cas de
! l'antidiffusion, il me semble que c'est le temps "1" qui va eviter
! l'amplification du mode numerique. On garde quand meme une EMA coef 0.5

        anyv3d(i,j,k,id_bnegdif)=              &
!                          ratio_bionegdif*(anyv3d(i,j,k  ,id_biobefo)  &
!                                   +x0*(  -anyv3d(i,j,k-2,id_biobefo)  &
!                                        +4*anyv3d(i,j,k-1,id_biobefo)  &
!                                        -6*anyv3d(i,j,k  ,id_biobefo)  &
!                                        +4*anyv3d(i,j,k+1,id_biobefo)  &
!                                          -anyv3d(i,j,k+2,id_biobefo)))
                            ratio_bionegdif*(bio_t(i,j,k  ,vb)  &
                                     +x0*(  -bio_t(i,j,k-2,vb)  &
                                          +4*bio_t(i,j,k-1,vb)  &
                                          -6*bio_t(i,j,k  ,vb)  &
                                          +4*bio_t(i,j,k+1,vb)  &
                                            -bio_t(i,j,k+2,vb)))
       enddo ; enddo ; enddo
! CL surface et fond:
       do k=kmax-1,kmax+1
!       do j=-1,jmax+2 ; do i=-1,imax+2
        do j= 1,jmax   ; do i= 1,imax  
         anyv3d(i,j,k,id_bnegdif)=anyv3d(i,j,kmax-2,id_bnegdif)+anyv3d(i,j,k,id_biobefo)-anyv3d(i,j,kmax-2,id_biobefo)
        enddo ; enddo
       enddo
       do k=0,2
!       do j=-1,jmax+2 ; do i=-1,imax+2
        do j= 1,jmax   ; do i= 1,imax  
         anyv3d(i,j,k,id_bnegdif)=anyv3d(i,j,3,id_bnegdif)+anyv3d(i,j,k,id_biobefo)-anyv3d(i,j,3,id_biobefo)
        enddo ; enddo
       enddo
#endif

      end subroutine advection_bio_negdif_fields

!..............................................................................

      subroutine advection_bio_quickest_coef_h_flux
      use module_principal
      use module_parallele
      implicit none

! Cette routine etant entierement independante de vb elle
! est appellee en dehors d'une boucle sur vb

! Details en:
!   https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit

      do k=1,kmax  !pmpmpmpmpmpmpmpmpmpmppmpmpmpmpmpmpm> !DEBUT DE BOUCLE SUR k

       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,1)=1.                                  &
                  /(dxdy_t(i,j)*0.5*( dz_t(i,j,k,0)     &
                                     +dz_t(i,j,k,1)))
       enddo ; enddo

! Note sur le role de upwindriver
! Si upwindriver=0 c'est que le schema d'advection doit etre 100% upwind. Cette condition
! s'obtient si le nombre de courant est artificiellement imposE A 1

! Nombre de courant vis a vis de u
       do j=1,jmax
       do i=1,imax+1
        xy_u(i,j,id_ncu)=min(1.,                             & !
        dti_fwsubio*abs(veldydz_u(i,j,k,1))                  & !
                        *max(xy_t(i,j,1),xy_t(i-1,j,1)) )    & ! Cas standard
              *0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j)) & ! si upwindriver=1

          +(1.-0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j)))  ! Cas particulier upwindriver=0 !18-01-17

       enddo
       enddo

! Nombre de courant vis a vis de v
       do j=1,jmax+1
       do i=1,imax
        xy_v(i,j,id_ncv)=min(1.,                             &
        dti_fwsubio*abs(veldxdz_v(i,j,k,1))                  &
                        *max(xy_t(i,j,1),xy_t(i,j-1,1)) )    & ! Cas standard
              *0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1)) & ! si upwindriver=1

          +(1.-0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1)))  ! Cas particulier upwindriver=0 !18-01-17
       enddo
       enddo

! Pre-coef etape intermediaire sur points flux u:
       do j=1,jmax
       do i=1,imax+1

! Rappel du flux quickest cas u>0
! F=u*( 0.5*(B(i)+B(i-1))-0.5*ncu*(B(i)-B(i-1))-(1-ncu**2)/6.*DRV2 )
! avec DRV2 la derivee seconde "masquee": B(i)-B(i-1)+(B(i-2)-B(i-1))*msk(i-2)
! Rappel du flux quickest cas u<0
! F=u*( 0.5*(B(i)+B(i-1))-0.5*ncu*(B(i-1)-B(i))-(1-ncu**2)/6.*DRV2 )
! avec DRV2 la derivee seconde "masquee": B(i-1)-B(i)+(B(i+1)-B(i))*msk(i+1)

        x1=0.5*(veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))
        x2=0.5*(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))

! noter que somme anyv3d(i,j,k,1:4)=veldydz_u(i,j,k,1)

! En facteur de B(i):
        anyv3d(i,j,k,3)=( & ! facteur B(i  )
        x1*(0.5-0.5*xy_u(i,j,id_ncu)-(1.-xy_u(i,j,id_ncu)**2)/6.) &
       +x2*(0.5+0.5*xy_u(i,j,id_ncu)+(1.-xy_u(i,j,id_ncu)**2)*(1.+mask_t(i+1,j,kmax))/6.) &
                                 )   ! facteur B(i  )

! En facteur de B(i-1):
        anyv3d(i,j,k,2)=( & ! facteur B(i-1)
        x1*(0.5+0.5*xy_u(i,j,id_ncu)+(1.-xy_u(i,j,id_ncu)**2)*(1.+mask_t(i-2,j,kmax))/6.) &
       +x2*(0.5-0.5*xy_u(i,j,id_ncu)-(1.-xy_u(i,j,id_ncu)**2)/6.)                         &
                                 )   ! facteur B(i-1)

! En facteur de B(i-2):
        anyv3d(i,j,k,1)=-x1*(1.-xy_u(i,j,id_ncu)**2)*mask_t(i-2,j,kmax)/6. ! facteur B(i-2)

! En facteur de B(i+1):
        anyv3d(i,j,k,4)=-x2*(1.-xy_u(i,j,id_ncu)**2)*mask_t(i+1,j,kmax)/6. ! facteur B(i+1)

       enddo
       enddo

! Pre-coef etape intermediaire sur points flux v:
       do j=1,jmax+1
       do i=1,imax

        x1=0.5*(veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))
        x2=0.5*(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))

! noter que somme anyv3d(i,j,k,5:8)=veldxdz_v(i,j,k,1)

! En facteur de B(j):
        anyv3d(i,j,k,7)= &
        x1*(0.5-0.5*xy_v(i,j,id_ncv)-(1.-xy_v(i,j,id_ncv)**2)/6.) &
       +x2*(0.5+0.5*xy_v(i,j,id_ncv)+(1.-xy_v(i,j,id_ncv)**2)*(1.+mask_t(i,j+1,kmax))/6.)

! En facteur de B(j-1):
        anyv3d(i,j,k,6)= &
        x1*(0.5+0.5*xy_v(i,j,id_ncv)+(1.-xy_v(i,j,id_ncv)**2)*(1.+mask_t(i,j-2,kmax))/6.) &
       +x2*(0.5-0.5*xy_v(i,j,id_ncv)-(1.-xy_v(i,j,id_ncv)**2)/6.)

! En facteur de B(j-2):
        anyv3d(i,j,k,5)=-x1*(1.-xy_v(i,j,id_ncv)**2)*mask_t(i,j-2,kmax)/6.

! En facteur de B(j+1):
        anyv3d(i,j,k,8)=-x2*(1.-xy_v(i,j,id_ncv)**2)*mask_t(i,j+1,kmax)/6.

       enddo
       enddo

      enddo       !pmpmpmpmpmpmpmpmpmpmppmpmpmpmpmpmpm> ! FIN DE BOUCLE SUR k


      end subroutine advection_bio_quickest_coef_h_flux

!.................................................................

      subroutine advection_bio_quickest_h_flux
      use module_principal
      use module_parallele
      use module_my_outputs !24-07-21
      implicit none
!     integer :: id_biobef_=0  & ! identifiant variable bio "before"
!     integer :: id_bioaft_=19 & ! identifiant variable bio "after"
      integer :: loop_


!     do vb=1,vbmax ! vbvbvbvb> boucle placEe dans le driver

      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_biobefo)=bio_t(i,j,k,vb)
! A partir de cette ligne le tableau bio_t(i,j,k,vb) est dupliquE dans anyv3d(i,j,k,id_biobefo). Il
! devient disponible calculer le terme supplementaire assurant la propriete de conservation de
! l'advection. Ici le reset du cumul (bio_t=0):
!      bio_t(i,j,k,vb)=0. !18-01-17
       anyv3d(i,j,k,id_bdiv)=0. !21-07-21
      enddo ; enddo ; enddo

      do loop_=1,loopmaxbio !iterative loop>

! RIVIERES: Ne pas ecraser la valeur imposee au point derriere la source (oU mask_t=0) !01-03-19
! en multipliant les coefs en facteur de bio_t(i,j,k,vb) par 1-mask et les coefs
! en facteurs des points voisins par mask

! Advection partielle direction Oi:
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

         anyv3d(i  ,j  ,k,id_bioaftr)=                                  &
         anyv3d(i  ,j  ,k,id_biobefo)*0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1)) &
                
         +dti_fwsubio*invdxdy_t(i,j)*mask_t(i,j,kmax)*( & !m°v°m>

! expression homogene A: +dti_fwbio*veldydz_u(i  ,j,k,1)*<bio(i  ,j,k)>:
!                                           anyv3d(i  ,j,k,3)*anyv3d(i  ,j,k,id_biobefo) &
!                                          +anyv3d(i  ,j,k,2)*anyv3d(i-1,j,k,id_biobefo) &
!                                          +anyv3d(i  ,j,k,1)*anyv3d(i-2,j,k,id_biobefo) &
!                                          +anyv3d(i  ,j,k,4)*anyv3d(i+1,j,k,id_biobefo) &
! expression homogene A: -dti_fwbio*veldydz_u(i+1,j,k,1)*<bio(i+1,j,k)>:
!                                          -anyv3d(i+1,j,k,3)*anyv3d(i+1,j,k,id_biobefo) &
!                                          -anyv3d(i+1,j,k,2)*anyv3d(i  ,j,k,id_biobefo) &
!                                          -anyv3d(i+1,j,k,1)*anyv3d(i-1,j,k,id_biobefo) &
!                                          -anyv3d(i+1,j,k,4)*anyv3d(i+2,j,k,id_biobefo) &
! Correction de divergence du/dx !18-01-17
!                          +(veldydz_u(i+1,j,k,1)-veldydz_u(i,j,k,1))*anyv3d(i  ,j  ,k,id_biobefo) &
! Equivalent mais plus compact
                        +(veldydz_u(i+1,j,k,1)-veldydz_u(i,j,k,1) &
                            + anyv3d(i  ,j,k,3)-anyv3d(i+1,j,k,2))*anyv3d(i  ,j,k,id_biobefo) &
                            +(anyv3d(i  ,j,k,2)-anyv3d(i+1,j,k,1))*anyv3d(i-1,j,k,id_biobefo) &
                            +(anyv3d(i  ,j,k,4)-anyv3d(i+1,j,k,3))*anyv3d(i+1,j,k,id_biobefo) &
                            + anyv3d(i  ,j,k,1)                   *anyv3d(i-2,j,k,id_biobefo) &
                                               -anyv3d(i+1,j,k,4) *anyv3d(i+2,j,k,id_biobefo) & 

                                                      )   !m°v°m>

! Cumul du terme de correction de divergence: (A oter A la toute fin de l'advection)
                         anyv3d(i  ,j  ,k,id_bdiv)= & !21-07-21
                         anyv3d(i  ,j  ,k,id_bdiv)  &
                        +anyv3d(i  ,j  ,k,id_biobefo)*(veldydz_u(i+1,j,k,1)-veldydz_u(i,j,k,1)) &
                                                   *dti_fwsubio*invdxdy_t(i,j)*mask_t(i,j,kmax)

       enddo ; enddo ; enddo

#ifdef bilanbio
       call my_outputs_zone1bioflux('i_bio',0)
#endif

       if(flag_nemoffline==1)call advection_bio_nemobottom(id_bioaftr)
       call vertmix_merged_levels_2t_any(1,id_bioaftr) !07-03-19

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_biobefo)=anyv3d(i,j,k,id_bioaftr)
       enddo ; enddo ; enddo

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
       call obc_mpi_anyv3d(1,id_biobefo,'zb')

! RIVIERES: Ne pas ecraser la valeur imposee au point derriere la source (oU mask_t=0) !01-03-19
! en multipliant les coefs en facteur de bio_t(i,j,k,vb) par 1-mask et les coefs
! en facteurs des points voisins par mask

! Advection partielle direction Oj:
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

         anyv3d(i  ,j  ,k,id_bioaftr)=                                  &
         anyv3d(i  ,j  ,k,id_biobefo)*0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1)) &
                
         +dti_fwsubio*invdxdy_t(i,j)*mask_t(i,j,kmax)*( & !m°v°m>

! expression homogene A: +dti_fwsubio*veldxdz_v(i,j  ,k,1)*<bio(i,j  ,k)>:
!                                           anyv3d(i,j  ,k,7)*anyv3d(i,j  ,k,id_biobefo) &
!                                          +anyv3d(i,j  ,k,6)*anyv3d(i,j-1,k,id_biobefo) &
!                                          +anyv3d(i,j  ,k,5)*anyv3d(i,j-2,k,id_biobefo) &
!                                          +anyv3d(i,j  ,k,8)*anyv3d(i,j+1,k,id_biobefo) &
! expression homogene A: +dti_fwsubio*veldxdz_v(i,j+1,k,1)*<bio(i,j+1,k)>:
!                                          -anyv3d(i,j+1,k,7)*anyv3d(i,j+1,k,id_biobefo) &
!                                          -anyv3d(i,j+1,k,6)*anyv3d(i,j  ,k,id_biobefo) &
!                                          -anyv3d(i,j+1,k,5)*anyv3d(i,j-1,k,id_biobefo) &
!                                          -anyv3d(i,j+1,k,8)*anyv3d(i,j+2,k,id_biobefo) &
! Correction de divergence dv/dy !18-01-17
!                         +(veldxdz_v(i,j+1,k,1)-veldxdz_v(i,j,k,1))*anyv3d(i  ,j  ,k,id_biobefo) &
! Equivalent mais plus compact
                        +(veldxdz_v(i,j+1,k,1)-veldxdz_v(i,j,k,1) &
                            + anyv3d(i,j  ,k,7)-anyv3d(i,j+1,k,6))*anyv3d(i,j  ,k,id_biobefo) &
                            +(anyv3d(i,j  ,k,6)-anyv3d(i,j+1,k,5))*anyv3d(i,j-1,k,id_biobefo) &
                            +(anyv3d(i,j  ,k,8)-anyv3d(i,j+1,k,7))*anyv3d(i,j+1,k,id_biobefo) &
                            + anyv3d(i,j  ,k,5)                   *anyv3d(i,j-2,k,id_biobefo) &
                                               -anyv3d(i,j+1,k,8) *anyv3d(i,j+2,k,id_biobefo) &

                                                      )   !m°v°m>

! Cumul du terme de correction de divergence: (A oter A la toute fin de l'advection)
                         anyv3d(i  ,j  ,k,id_bdiv)= & !21-07-21
                         anyv3d(i  ,j  ,k,id_bdiv)  &
                        +anyv3d(i  ,j  ,k,id_biobefo)*(veldxdz_v(i,j+1,k,1)-veldxdz_v(i,j,k,1)) &
                                                   *dti_fwsubio*invdxdy_t(i,j)*mask_t(i,j,kmax)

       enddo ; enddo ; enddo

#ifdef bilanbio
       call my_outputs_zone1bioflux('j_bio',0)
#endif

       if(flag_nemoffline==1)call advection_bio_nemobottom(id_bioaftr)
       call vertmix_merged_levels_2t_any(1,id_bioaftr) !07-03-19

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_biobefo)=anyv3d(i,j,k,id_bioaftr)
       enddo ; enddo ; enddo

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
       call obc_mpi_anyv3d(1,id_biobefo,'zb')

      enddo                 !iterative loop>

!     enddo         ! vbvbvbvb>


      end subroutine advection_bio_quickest_h_flux

!.................................................................
      subroutine advection_bio_erreur_stop
      use module_principal
      use module_parallele
      implicit none

      if(par%rank==0) then !000>

       write(6,*)'The "bio_t" tracers should be numbered in descending' &
       ,' order of sedimentation speed, to avoid unnecessary' &
       ,' recalculation of advection and vertical mixing coefficients.'
       write(6,*)'vb  ,wsed(1,vb)  :',vb  ,wsed(1,vb)
       write(6,*)'vb-1,wsed(1,vb-1):',vb-1,wsed(1,vb-1)

      endif                !000>
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
#endif
      stop 'advection_bio_erreur_stop'

      end subroutine advection_bio_erreur_stop
!.........................................................................
      subroutine advection_bio_rmnegval_vert !03-11-17
      use module_principal ; use module_parallele ; use module_s
      implicit none

! Détails et evolutions possibles dans: !06-12-16
! https://docs.google.com/document/d/15QsKJhPFiOjxU3os7omfhrhhxNYZBsRTprBmsyuNguc/edit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculer un bilan pour verification des proprietes de conservation
!     sum1=0.
!     sum2=0.
!     vb=1
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!      sum1=sum1+(dz_t(i,j,k,0)+dz_t(i,j,k,1))*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)
!      sum2=sum2+(dz_t(i,j,k,0)+dz_t(i,j,k,1))*dxdy_t(i,j)*mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)*bio_t(i,j,k,vb)
!     enddo
!     enddo
!     enddo
!     call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
!     if(par%rank==0)write(66,*)iteration3d,sum2glb/sum1glb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! En dehors de la boucle sur vb (puisque ne depend pas de bio_t) la mise
! A zero des C.L. de surface dt de fond de anyv3d(:,:,:,3)
      do j=1,jmax ; do i=1,imax
       do k=1,kmin_w(i,j)
        anyv3d(i,j,k,3)=0.
       enddo
      enddo ; enddo
      do j=1,jmax ; do i=1,imax
       anyv3d(i,j,kmax+1,3)=0.
      enddo ; enddo


! Boucle sur vb
      do vb=1,vbmax

      do j=1,jmax ; do i=1,imax

       do k=kmin_w(i,j)+1,kmax
        anyv3d(i,j,k,3)=                                                           &
               min(  abs(( dz_t(i,j,k-1,0)  +dz_t(i,j,k-1,1))*bio_t(i,j,k-1,vb))   &
                    ,abs(( dz_t(i,j,k  ,0)  +dz_t(i,j,k  ,1))*bio_t(i,j,k  ,vb)) ) & !fermer min
                   /(abs((bio_t(i,j,k  ,vb)-bio_t(i,j,k-1,vb)))+1.E-20)            &
       *max(0.,-sign(1.d0,bio_t(i,j,k  ,vb)*bio_t(i,j,k-1,vb)))                    &
                        *(bio_t(i,j,k  ,vb)-bio_t(i,j,k-1,vb))
       enddo

      enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax

        bio_t(i,j,k,vb)=       &
        bio_t(i,j,k,vb)+(      & !ooooo>

        anyv3d(i,j,k+1,3)-anyv3d(i,j,k,3) &

                        )      & !ooooo> 
                         /(dz_t(i,j,k,0)+dz_t(i,j,k,1))


      enddo ; enddo ; enddo

      enddo ! vb

      call s_cpu('advection_bio_rmnegval_vert',0)
      end subroutine advection_bio_rmnegval_vert
!.........................................................................
