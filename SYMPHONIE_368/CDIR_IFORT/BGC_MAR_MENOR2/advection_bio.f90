










      subroutine advection_bio
!______________________________________________________________________
! SYMPHONIE ocean model
! release 384 - last update: 09-02-24
!______________________________________________________________________
!    _________                    .__                  .__       m�v�m !
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
!         14/12/01: bienvenue � ITIMEBIO
!         17/12/01: amenagements pour cas forward (economie de RAM)
!         26/12/02: amenagements pour remise en service schema forward
!         14/02/03: variable xy_z(i,j,1) �tait non initialis�e. Ajout boucle 104
!                   de mise � z�ro.
!         19/03/03: date de debut et de fin de calcul pour traceur
!         04/05/06: amenagements pour compiler en double precision
!         04/07/06: boucle 103: j avant i. Regroupement de constantes.
!         20/04/07: passage � la coordonn�es curviligne
! 2009.3  01-10-09: suppression kount_bio
!         02-10-09: dti remplac� par dti_fw
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
! v310    05-11-21  modifs pour VQS 2
! v316    21-12-21  La partie diffusion negative (voir id_bnegdif) est commentee
!                   en raison de son cout et parce que le modele eco3m
!                   n'a pas ete encore calibre par rapport A cela. Il
!                   suffit d'enlever les commentaires pour reactualiser
!                   l'option.
! v320    03-01-22  dti_fwsubio renommE dti_fwsub_ dans l'advection verticale car les pas 
!                   de temps horizontal et vertical peuvent etre differents
! v374    18-11-23  modifs advection: C2 limited sur la verticale et nouveau limiteur
! v379    02-01-24  - suite du point precedent: debug anyv3d(i,j,k,id_cnk)
!                   - C2 limited sur l'horizontal
!                   - limiteur d'oscillations remplace limiteur de valeur negatives
!         03-01-24  Si w de surface non nul et negatif alors w explicite et OBC associee
! v384    09-02-24  limiteur diffErE
!..............................................................................

       if(flag_1dv==1)return !25-04-16 pas d'advection bio si flag_1dv=1

! Decommenter cette ligne pour diag min max sur bio_t
!      call advection_bio_min_max !08-02-17


! Pour obc en kmax+1 dans mixsed_bio !03-01-24
       bio_t(:,:,kmax+1,:)=bio_t(:,:,kmax,:) ; check_obcbio=bio_t(imax/2,jmax/2,kmax+1,1)

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
!.....................................................................
!.....................................................................
!________________________________________________________________________________
!.....................................................................
!.....................................................................
!________________________________________________________________________________
      subroutine advection_biofinalize(ind4_)
      use module_principal ; use module_parallele
      implicit none
      integer ind4_

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
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!#endif
!      stop 'bibi'



      end subroutine advection_bio_get_substep

!.........................................................................

      subroutine advection_bio_rmnegval !11-07-16
      use module_principal ; use module_parallele
      implicit none

      stop 'UTILISER advection_bio_rmnegval3d_plus'

! D�tails et evolutions possibles dans: !06-12-16
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
      do k=1,kmax !05-11-21
      do j=1,jmax ; do i=1,imax
!       do k=kmerged_t(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(veldydz_u(i  ,j  ,k,1))   &
                              ,abs(veldydz_u(i+1,j  ,k,1)))  &
                              ,abs(veldxdz_v(i  ,j  ,k,1)))  &
                              ,abs(veldxdz_v(i  ,j+1,k,1)))  &
                    *dti_fw                                  &
                    *wetmask_t(i,j)                          &
                       /(dxdy_t(i,j)*                        &
                         0.5*( dz_t(i,j,k,0)+dz_t(i,j,k,1))) &
               ) ! max
      enddo       ; enddo
      enddo
! Niveaux <= kmerged
!     do j=1,jmax ; do i=1,imax
!       sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
!       do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
!        sum1=sum1+veldydz_u(i  ,j  ,k,1)
!        sum2=sum2+veldydz_u(i+1,j  ,k,1)
!        sum3=sum3+veldxdz_v(i  ,j  ,k,1)
!        sum4=sum4+veldxdz_v(i  ,j+1,k,1)
!        sum5=sum5     +dz_t(i  ,j  ,k,0)
!        sum6=sum6     +dz_t(i  ,j  ,k,2)
!       enddo
!        x1=max(x1,max(max(max(abs(sum1)  &
!                             ,abs(sum2)) &
!                             ,abs(sum3)) &
!                             ,abs(sum4)) &
!                   *dti_fw                                  &
!                   *wetmask_t(i,j)                          &
!                     /(dxdy_t(i,j)*min(sum5,sum6))          &
! methode 2
!              ) ! fermeture max
!     enddo       ; enddo


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
      if(loopmaxbio>50) then !m�v�m> !06-05-19
      do k=1,kmax !05-11-21
      do j=1,jmax ; do i=1,imax
!       do k=kmerged_t(i,j)+1,kmax
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

      enddo       ; enddo
      enddo
! Niveaux <= kmerged
!     do j=1,jmax ; do i=1,imax
!       sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
!       do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
!        sum1=sum1+veldydz_u(i  ,j  ,k,1)
!        sum2=sum2+veldydz_u(i+1,j  ,k,1)
!        sum3=sum3+veldxdz_v(i  ,j  ,k,1)
!        sum4=sum4+veldxdz_v(i  ,j+1,k,1)
!        sum5=sum5     +dz_t(i  ,j  ,k,0)
!        sum6=sum6     +dz_t(i  ,j  ,k,2)
!       enddo
!        x2=       max(max(max(abs(sum1)  &
!                             ,abs(sum2)) &
!                             ,abs(sum3)) &
!                             ,abs(sum4)) &
!                   *dti_fw                                  &
!                   *wetmask_t(i,j)                          &
!                     /(dxdy_t(i,j)*min(sum5,sum6))           

!            if(int(x2+1)>50) then !(0v0)>
!              write(10+par%rank,*)'---------------'
!              write(10+par%rank,*)'x2=            ',x2
!              write(10+par%rank,*)'i,j,k loc      ',i,j,k
!              write(10+par%rank,*)'i,j,k glb      ',i+par%timax(1),j+par%tjmax(1),k
!              write(10+par%rank,*)'mask_t         ',mask_t(i,j,k)
!              write(10+par%rank,*)'h_w            ',h_w(i,j)
!              write(10+par%rank,*)'sum1           ',sum1
!              write(10+par%rank,*)'sum2           ',sum1
!              write(10+par%rank,*)'sum3           ',sum1
!              write(10+par%rank,*)'sum4           ',sum1
!              write(10+par%rank,*)'sum5           ',sum1
!              write(10+par%rank,*)'sum6           ',sum1
!              write(10+par%rank,*)'---------------'
!            endif                 !(0v0)>

!     enddo       ; enddo
      call graph_out
      stop 'loopmaxbio>50 see fortxxx error files'
      endif                  !m�v�m>

      end subroutine advection_bio_substep

!.................................................................

      subroutine advection_bio_quickest_coef_v
      use module_principal
      use module_parallele
      implicit none



!-----------------------------------------------------------------------------------------    
! CALCUL DE LA VITESSE VERTICALE EXPLICITE COHERENT AVEC LE CALCUL DE LA VITESSE IMPLICITE
! dans mixsed_bio le calcul de la vitesse implicite est l'exacte complementaire du calcul 
! suivant. Si l'un change, penser A corriger l'autre
       do k=2,kmax ; do j=1,jmax ; do i=1,imax

          anyv3d(i,j,k,id_webio)=max(min(                        &

          omega_w(i,j,k,1)+wsed(k,vb)*wsed_explicit(vb)          &
        , 0.5*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1),dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))*inv_dti_fw) &
        ,-0.5*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1),dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))*inv_dti_fw) &

                                   *wetmask_t(i,j)                 !13-11-16

       enddo       ; enddo       ; enddo
! Fond:
       k=1
       do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,id_webio)=0.
       enddo ; enddo
! Surface:
       k=kmax+1
       do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,id_webio)=0.
       enddo ; enddo
!-----------------------------------------------------------------------------------------    

! Nombre de courant vis a vis de omega explicite
       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_cnk)=                                   & !02-01-24
                    abs(anyv3d(i,j,k,id_webio))*dti_fw         & 
                 /(0.5*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1)     &
                          ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1)))   &
                 *upwindriver_t(i,j)                           & !standard si upwindriver=1
             +(1.-upwindriver_t(i,j))                            !      =1 si upwindriver=0 
       enddo ; enddo ; enddo

      end subroutine advection_bio_quickest_coef_v

!.................................................................

      subroutine advection_bio_quickest_v
      use module_principal
      use module_parallele
      implicit none
!     integer :: id_biobef_=0   & ! identifiant variable bio "before"
!     integer ::  loopmaxbio_v_ &
!                ,loop_

! Calculer le flux vertical: !18-11-23
       do k1=1,kmax+1
! Attention les lignes suivantes ne bornent que id_biobefo et pas id_flk ni id_webio qui utilisent k1
       km2=max(k1-2,1)
       km1=max(k1-1,1)
       kp1=min(k1+1,kmax)
       kp2=min(k1+2,kmax)
       k  =min(k1,kmax)
        do j=1,jmax ; do i=1,imax

        anyv3d(i,j,k1,id_flk)=  &
!#ifdef bidon
!SCHEMA C2+(cn**2)*GRAD + limiteur d'oscillations (mieux qu'un simple limiteur de valeurs negatives)
        0.5*    anyv3d(i,j,k1,id_webio) *(anyv3d(i,j,k,id_biobefo)+anyv3d(i,j,km1,id_biobefo)) &
       -0.5*abs(anyv3d(i,j,k1,id_webio))*(anyv3d(i,j,k,id_biobefo)-anyv3d(i,j,km1,id_biobefo)) &
           *( & !ooo>
             anyv3d(i,j,k,id_cnk)**2 &
! Limiteur d'oscillation (mieux que simple limiteur de valeurs negatives)
! si stable   alors =0 de sorte que cn**2 en facteur de la diffusion
! si instable alors =1-cn**2 de sorte que 1 en facteur de la diffusion
!     +max(0.d0,sign(1.d0,                         &
!      (abs(anyv3d(i,j,k1,id_webio))+anyv3d(i,j,k1,id_webio))            &
!     *( abs(anyv3d(i,j,k  ,id_biobefo)-anyv3d(i,j,km1,id_biobefo))*anyv3d(i,j,k,id_cnk)  &
!       -abs(anyv3d(i,j,km1,id_biobefo)-anyv3d(i,j,km2,id_biobefo)) )      &

!     +(abs(anyv3d(i,j,k1,id_webio))-anyv3d(i,j,k1,id_webio))            &
!     *( abs(anyv3d(i,j,k  ,id_biobefo)-anyv3d(i,j,km1,id_biobefo))*anyv3d(i,j,k,id_cnk)  &
!       -abs(anyv3d(i,j,kp1,id_biobefo)-anyv3d(i,j,k  ,id_biobefo)) )      &
!          ))*(1-anyv3d(i,j,k,id_cnk)**2)                     &

! Le limiteur ci-dessus ne convient pas toujours (par ex si courant non constant)
! on prefere donc le limiteur "diffErE" suivant: !09-02-24
      +max(0.d0,sign(1.d0,                         &
       (abs(anyv3d(i,j,k1,id_webio))+anyv3d(i,j,k1,id_webio))            &
      *(((anyv3d(i,j,k  ,id_biobefo)-anyv3d(i,j,km1,id_biobefo))**2)*anyv3d(i,j,k,id_cnk)  &
        -(anyv3d(i,j,k  ,id_biobefo)-anyv3d(i,j,km1,id_biobefo))    &
        *(anyv3d(i,j,km1,id_biobefo)-anyv3d(i,j,km2,id_biobefo)) )  &

      +(abs(anyv3d(i,j,k1,id_webio))-anyv3d(i,j,k1,id_webio))            &
      *(((anyv3d(i,j,k  ,id_biobefo)-anyv3d(i,j,km1,id_biobefo))**2)*anyv3d(i,j,k,id_cnk)  &
        -(anyv3d(i,j,k  ,id_biobefo)-anyv3d(i,j,km1,id_biobefo))    &
        *(anyv3d(i,j,kp1,id_biobefo)-anyv3d(i,j,k  ,id_biobefo)) )  &
           ))*(1-anyv3d(i,j,k,id_cnk)**2)                     &


            )   !ooo>
!#endif

       enddo ; enddo
      enddo



      do k=1,kmax ; do j=1,jmax ; do i=1,imax

                   anyv3d(i  ,j  ,k  ,id_bioaftr) &
                  =anyv3d(i  ,j  ,k  ,id_biobefo) &
                 -anyv3d(i,j,k,id_dtovedz)*( & !ooo>
             anyv3d(i,j,k+1,id_flk)-anyv3d(i,j,k,id_flk) &
           -(anyv3d(i,j,k+1,id_webio )-anyv3d(i,j,k,id_webio))*anyv3d(i,j,k,id_biobefo) &
                                           )   !ooo>

                   anyv3d(i  ,j  ,k,id_bdiv)= & !21-07-21
                   anyv3d(i  ,j  ,k,id_bdiv)  &
                  +anyv3d(i  ,j  ,k  ,id_biobefo)*(anyv3d(i,j,k+1,id_webio)-anyv3d(i,j,k,id_webio)) &
                                                 *dti_fw*mask_t(i,j,kmax)

      enddo ; enddo ; enddo

! Les lignes suivantes sont inutiles si pas de subcycling sur la verticale:
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!       anyv3d(i,j,k,id_biobefo)=anyv3d(i,j,k,id_bioaftr)
!     enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax

! Attention ici bio_t est homogene A bio_t*dz_t, par consequent dans l'etape suivante dans mixsed_bio
! bio_t n'est plus multipliE par dz_t:
         bio_t(i,j,k,vb)=                                           &
        anyv3d(i,j,k,id_bioaftr)*0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))  & !01-05-19
       -anyv3d(i,j,k,id_bdiv)    !21-07-21

      enddo ; enddo ; enddo

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

! Verifier dimensions 2DH de bio_t
       ub4=ubound(bio_t) ; lb4=lbound(bio_t)
       if(ub4(2)<jmax+2)stop 'Erreur max dimension 2 bio_t'
       if(lb4(2)>-1)    stop 'Erreur min dimension 2 bio_t'
       if(ub4(1)<imax+2)stop 'Erreur max dimension 1 bio_t'
       if(lb4(1)>-1)    stop 'Erreur min dimension 1 bio_t'

! Verifier 4eme dimension de anyv3d
       ub4=ubound(anyv3d) ; lb4=lbound(anyv3d)
       if(lb4(4)>0 .or.ub4(4)<12) then !>>>>>> !21-07-21
        deallocate(anyv3d)
          allocate(anyv3d(-1:imax+2,-1:jmax+2,0:kmax+1,0:12)) ; anyv3d=0 !21-07-21
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

      stop 'ne plus utilise advection_bio_rmnegval3d'

! D�tails et evolutions possibles dans: !06-12-16
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

       stop 'ne plus utiliser advection_bio_rmnegval3d_plus'

! D�tails et evolutions possibles dans: !06-12-16
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
          *max(0.d0,-sign(1.d0,anyv3d(i+1,j,k,0)*anyv3d(i,j,k,0)))  & !=1 si l'un des deux seulement est <0. O sinon. !03-11-17
       +anyv3d(i  ,j,k,1)*(anyv3d(i-1,j,k,0)-anyv3d(i,j,k,0))   &
          *max(0.d0,-sign(1.d0,anyv3d(i-1,j,k,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j+1,k,2)*(anyv3d(i,j+1,k,0)-anyv3d(i,j,k,0))   &
          *max(0.d0,-sign(1.d0,anyv3d(i,j+1,k,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j  ,k,2)*(anyv3d(i,j-1,k,0)-anyv3d(i,j,k,0))   &
          *max(0.d0,-sign(1.d0,anyv3d(i,j-1,k,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j,k+1,3)*(anyv3d(i,j,k+1,0)-anyv3d(i,j,k,0))   &
          *max(0.d0,-sign(1.d0,anyv3d(i,j,k+1,0)*anyv3d(i,j,k,0)))  &
       +anyv3d(i,j,k  ,3)*(anyv3d(i,j,k-1,0)-anyv3d(i,j,k,0))   &
          *max(0.d0,-sign(1.d0,anyv3d(i,j,k-1,0)*anyv3d(i,j,k,0)))  &

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

        anyv3d(i,j,k,id_bnegdif)=ratio_bionegdif*( & !m�v�m>
!                               0.5*anyv3d(i,j,k  ,id_biobefo)   &
!                            +0.25*(anyv3d(i,j,k+1,id_biobefo)   &
!                                  +anyv3d(i,j,k-1,id_biobefo))  &
                                 0.5*bio_t(i,j,k  ,vb)   &
                              +0.25*(bio_t(i,j,k+1,vb)   &
                                    +bio_t(i,j,k-1,vb))  &
                                                 )   !m�v�m>
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


      end subroutine advection_bio_negdif_fields

!..............................................................................

      subroutine advection_bio_quickest_coef_h_flux
      use module_principal
      use module_parallele
      implicit none
!     integer :: id_biobefo=0  & ! identifiant variable bio "before"
!               ,id_bioaftr=19   ! identifiant variable bio "after"

! Cette routine etant entierement independante de vb elle
! est appellee en dehors d'une boucle sur vb

! Details en:
!   https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit

! Nombre de courant vis a vis de u
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1

        anyv3d(i,j,k,id_nci)=min(1.,                                    &
        dti_fwsubio*abs(veldydz_u(i,j,k,1))                              &
            /min( dxdy_t(i  ,j)*0.5*( dz_t(i  ,j,k,0)+dz_t(i  ,j,k,1))   &
                 ,dxdy_t(i-1,j)*0.5*( dz_t(i-1,j,k,0)+dz_t(i-1,j,k,1)))) &
              *0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j))             & ! si upwindriver=1

          +(1.-0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j)))  ! Cas particulier upwindriver=0 !18-01-17

       enddo ; enddo ; enddo

! Nombre de courant vis a vis de v
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax

        anyv3d(i,j,k,id_ncj)=min(1.,                                    &
        dti_fwsubio*abs(veldxdz_v(i,j,k,1))                              & 
            /min( dxdy_t(i,j  )*0.5*( dz_t(i,j  ,k,0)+dz_t(i,j  ,k,1))   &
                 ,dxdy_t(i,j-1)*0.5*( dz_t(i,j-1,k,0)+dz_t(i,j-1,k,1)))) &
              *0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1))             & ! si upwindriver=1

          +(1.-0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1)))  ! Cas particulier upwindriver=0 !18-01-17

       enddo ; enddo ; enddo

! En facteur de l'advection horizontale: !18-11-23
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,id_dtovedxdydz)=                          &
                  dti_fwsubio/( dxdy_t(i,j)*0.5*(dz_t(i,j,k,0)   &
                                                +dz_t(i,j,k,1)))*mask_t(i,j,kmax)
       enddo ; enddo ; enddo

! En facteur de l'advection verticale: !18-11-23
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,id_dtovedz)=                          &
                  dti_fw/( 0.5*(dz_t(i,j,k,0)                &
                              +dz_t(i,j,k,1)))*mask_t(i,j,kmax)
       enddo ; enddo ; enddo

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

! EXPRIMER LE Flux Oi  -DEBUT !18-11-23
      do k=1,kmax  ; do j=1,jmax ; do i=1,imax+1


!#ifdef bidon
!....................................................
! SCHEMA 2
! C2 limited et son limiteur d'oscillation (mieux qu'un simple limiteur de valeurs negatives)
        anyv3d(i,j,k,id_fli)=  & 
! Terme C2 et sa diffusion avec cn(i)**2 en facteur
            0.5*veldydz_u(i,j,k,1) *(anyv3d(i,j,k,id_biobefo)+anyv3d(i-1,j,k,id_biobefo)) &
       -0.5*abs(veldydz_u(i,j,k,1))*(anyv3d(i,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo)) &
           *(anyv3d(i,j,k,id_nci)**2 &

! Limiteur d'oscillation (mieux que simple limiteur de valeurs negatives)
! si stable   alors =0 de sorte que cn**2 en facteur de la diffusion
! si instable alors =1-cn**2 de sorte que 1 en facteur de la diffusion
!     +max(0.d0,sign(1.d0,                         &
!      (abs(veldydz_u(i,j,k,1))+veldydz_u(i,j,k,1))            &
!     *( abs(anyv3d(i  ,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo))*anyv3d(i,j,k,id_nci)  &
!       -abs(anyv3d(i-1,j,k,id_biobefo)-anyv3d(i-2,j,k,id_biobefo)) )      &

!     +(abs(veldydz_u(i,j,k,1))-veldydz_u(i,j,k,1))            &
!     *( abs(anyv3d(i  ,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo))*anyv3d(i,j,k,id_nci)  &
!       -abs(anyv3d(i+1,j,k,id_biobefo)-anyv3d(i  ,j,k,id_biobefo)) )      &
!          ))*(1-anyv3d(i,j,k,id_nci)**2)                     &

! Le schema ci-dessus convient tout le temps si le champ de courant est constant
! mais pas dans les autres cas, aussi on prefere un schema qui, en plus
! de la contrainte precedente, agit auss en diff�r� sur des oscillations existante !09-02-24
      +max(0.d0,sign(1.d0,                         &
       (abs(veldydz_u(i,j,k,1))+veldydz_u(i,j,k,1))            &
      *(((anyv3d(i  ,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo))**2)*anyv3d(i,j,k,id_nci)  &
        -(anyv3d(i  ,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo))      &
        *(anyv3d(i-1,j,k,id_biobefo)-anyv3d(i-2,j,k,id_biobefo))   )  &

      +(abs(veldydz_u(i,j,k,1))-veldydz_u(i,j,k,1))            &
      *(((anyv3d(i  ,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo))**2)*anyv3d(i,j,k,id_nci)  &
        -(anyv3d(i  ,j,k,id_biobefo)-anyv3d(i-1,j,k,id_biobefo))      &
        *(anyv3d(i+1,j,k,id_biobefo)-anyv3d(i  ,j,k,id_biobefo)) )    &
           ))*(1-anyv3d(i,j,k,id_nci)**2)                     &

                              )   !-C2-UP**2
!....................................................
!#endif

       enddo ; enddo ; enddo ! boucle k
! EXPRIMER LE Flux Oi - FIN

! RIVIERES: Ne pas ecraser la valeur imposee au point derriere la source (oU mask_t=0) !01-03-19
! en multipliant les coefs en facteur de bio_t(i,j,k,vb) par 1-mask et les coefs
! en facteurs des points voisins par mask

! Advection partielle direction Oi:
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

! Alternative divergence de flux
        anyv3d(i,j,k,id_bioaftr)=  &
        anyv3d(i,j,k,id_biobefo)-( &
               anyv3d(i+1,j,k,id_fli)-anyv3d(i,j,k,id_fli) &
          -(veldydz_u(i+1,j,k,1)  -veldydz_u(i,j,k,1))*anyv3d(i  ,j ,k,id_biobefo) &
                                 )                           &
          *anyv3d(i,j,k,id_dtovedxdydz) 

! Cumul du terme de correction de divergence: (A oter A la toute fin de l'advection)
                         anyv3d(i  ,j  ,k,id_bdiv)= & !21-07-21
                         anyv3d(i  ,j  ,k,id_bdiv)  &
                        +anyv3d(i  ,j  ,k,id_biobefo)*(veldydz_u(i+1,j,k,1)-veldydz_u(i,j,k,1)) &
                                                   *dti_fwsubio*invdxdy_t(i,j)*mask_t(i,j,kmax)

       enddo ; enddo ; enddo


       if(flag_nemoffline==1)call advection_bio_nemobottom(id_bioaftr)
!      call vertmix_merged_levels_2t_any(1,id_bioaftr) !07-03-19

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_biobefo)=anyv3d(i,j,k,id_bioaftr)
       enddo ; enddo ; enddo

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
       call obc_mpi_anyv3d(1,id_biobefo,'zb')


! EXPRIMER LE Flux Oj  -DEBUT !18-11-23
      do k=1,kmax  ; do j=1,jmax+1 ; do i=1,imax


!#ifdef bidon
!....................................................
! SCHEMA 2
! C2 limited et son limiteur d'oscillation (mieux qu'un simple limiteur de valeurs negatives)
        anyv3d(i,j,k,id_flj)=  & 
! Terme C2 et sa diffusion avec cn(i)**2 en facteur
            0.5*veldxdz_v(i,j,k,1) *(anyv3d(i,j,k,id_biobefo)+anyv3d(i,j-1,k,id_biobefo)) &
       -0.5*abs(veldxdz_v(i,j,k,1))*(anyv3d(i,j,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo)) &
           *(anyv3d(i,j,k,id_ncj)**2 &

! Limiteur d'oscillation (mieux que simple limiteur de valeurs negatives)
! si stable   alors =0 de sorte que cn**2 en facteur de la diffusion
! si instable alors =1-cn**2 de sorte que 1 en facteur de la diffusion
!     +max(0.d0,sign(1.d0,                         &
!      (abs(veldxdz_v(i,j,k,1))+veldxdz_v(i,j,k,1))            &
!     *( abs(anyv3d(i,j  ,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo))*anyv3d(i,j,k,id_ncj)  &
!       -abs(anyv3d(i,j-1,k,id_biobefo)-anyv3d(i,j-2,k,id_biobefo)) )      &

!     +(abs(veldxdz_v(i,j,k,1))-veldxdz_v(i,j,k,1))            &
!     *( abs(anyv3d(i,j  ,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo))*anyv3d(i,j,k,id_ncj)  &
!       -abs(anyv3d(i,j+1,k,id_biobefo)-anyv3d(i,j  ,k,id_biobefo)) )      &
!          ))*(1-anyv3d(i,j,k,id_ncj)**2)                     &

! Le schema ci-dessus convient tout le temps si le champ de courant est constant
! mais pas dans les autres cas, aussi on prefere un schema qui, en plus
! de la contrainte precedente, agit auss en diff�r� sur des oscillations existante !09-02-24
      +max(0.d0,sign(1.d0,                         &
       (abs(veldxdz_v(i,j,k,1))+veldxdz_v(i,j,k,1))            &
      *(((anyv3d(i,j  ,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo))**2)*anyv3d(i,j,k,id_ncj)  &
        -(anyv3d(i,j  ,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo))   &
        *(anyv3d(i,j-1,k,id_biobefo)-anyv3d(i,j-2,k,id_biobefo)) ) &

      +(abs(veldxdz_v(i,j,k,1))-veldxdz_v(i,j,k,1))            &
      *(((anyv3d(i,j  ,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo))**2)*anyv3d(i,j,k,id_ncj)  &
        -(anyv3d(i,j  ,k,id_biobefo)-anyv3d(i,j-1,k,id_biobefo))   &
        *(anyv3d(i,j+1,k,id_biobefo)-anyv3d(i,j  ,k,id_biobefo)) ) &
           ))*(1-anyv3d(i,j,k,id_ncj)**2)                     &

                              )   !-C2-UP**2
!....................................................
!#endif

       enddo ; enddo ; enddo ! boucle k
! EXPRIMER LE Flux Oj - FIN


! RIVIERES: Ne pas ecraser la valeur imposee au point derriere la source (oU mask_t=0) !01-03-19
! en multipliant les coefs en facteur de bio_t(i,j,k,vb) par 1-mask et les coefs
! en facteurs des points voisins par mask

! Advection partielle direction Oj:
       do k=1,kmax ; do j=1,jmax ; do i=1,imax

! Alternative divergence de flux
        anyv3d(i,j,k,id_bioaftr)=  &
        anyv3d(i,j,k,id_biobefo)-( &
               anyv3d(i,j+1,k,id_flj)-anyv3d(i,j,k,id_flj) &
          -(veldxdz_v(i,j+1,k,1)  -veldxdz_v(i,j,k,1))*anyv3d(i  ,j ,k,id_biobefo) &
                                 ) &
          *anyv3d(i,j,k,id_dtovedxdydz) 

! Cumul du terme de correction de divergence: (A oter A la toute fin de l'advection)
                         anyv3d(i  ,j  ,k,id_bdiv)= & !21-07-21
                         anyv3d(i  ,j  ,k,id_bdiv)  &
                        +anyv3d(i  ,j  ,k,id_biobefo)*(veldxdz_v(i,j+1,k,1)-veldxdz_v(i,j,k,1)) &
                                                   *dti_fwsubio*invdxdy_t(i,j)*mask_t(i,j,kmax)

       enddo ; enddo ; enddo


       if(flag_nemoffline==1)call advection_bio_nemobottom(id_bioaftr)
!      call vertmix_merged_levels_2t_any(1,id_bioaftr) !07-03-19

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
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
      stop 'advection_bio_erreur_stop'

      end subroutine advection_bio_erreur_stop
!.........................................................................
      subroutine advection_bio_rmnegval_vert !03-11-17
      use module_principal ; use module_parallele ; use module_s
      implicit none

      stop 'ne plus utiliser advection_bio_rmnegval_vert'

! D�tails et evolutions possibles dans: !06-12-16
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
     *max(0.d0,-sign(1.d0,bio_t(i,j,k  ,vb)*bio_t(i,j,k-1,vb)))                    &
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
