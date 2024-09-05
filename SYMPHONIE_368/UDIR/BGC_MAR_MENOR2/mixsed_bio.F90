      subroutine mixsed_bio
!______________________________________________________________________
! SYMPHONIE ocean model
! release 303 - last update: 03-08-21
!______________________________________________________________________
      use module_principal ; use module_parallele ; use module_my_outputs
      implicit none
      integer :: i1_,i2_,j1_,j2_ &
                ,id_wi_=0        ! identifiant vitesse verticale residuelle implicite
!                id_kh_over_dz=1 dans module_principal
#ifdef synopsis
       subroutinetitle='mixsed_bio'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!......................................................................
! Version date      Description des modifications
!         01/11/01: inversion boucles i j
!         13/12/01: les coefs de la matrice principale sont séparés du
!                   terme de sedimentation qui est maintenant calculé
!                   à la fin et de maniere optionnelle (si wsed=0 pas
!                   de calculs). Ainsi séparé, on peut éviter de
!                   recalculer les coefs de la matrice principale à partir
!                   du 2eme traceurs puisqu'ils sont identiques à ceux
!                   du premier. (En toute rigueur ils sont les mêmes que
!                   ceux de T et S mais on ne va pas prendre de risques et tout
!                   de suite envisager le cas offline où la bio tourne sans
!                   la dynamique!
!         14/12/01: Bienvenue à ITMEBIO et à TENDANCEBIO_Z
!         25/10/02: Fluxbio_t homogene à grandeur primaire, n'est plus
!                   a etre multiplie par DTI
!         28/02/03: Amenagements pour coordonnee verticale hybride
!         19/03/03: date de debut et de fin de calcul pour traceur
!         02/12/04: suite du point precedent... le cas particulier
!                   KB.EQ.1 occasionnait un bug si le traceur n°1
!                   sort de la fenêtre de calcul alors que les autres
!                   traceurs continuent d'être calculés. J'ai jugé
!                   plus prudent (en attendant de mettre une selection
!                   plus sioux) de restreindre le test au seul update
!                   de BIOHZ et de calculer systematiquement les systemes
!         21/06/06: adaptation pour version 1D
!         05/03/07: si une variable bio devient negative, tendancebio est seuille
!         20/03/08: suppression de la modification precedente
! 2009.3  01-10-09: utilisation des nouveaux facteurs d'echelle verticale
!         02-10-09: dti_fw remplace dti
! 2010.2  22-12-09  facteur d'echelle pris au temps t+1
! 2010.4  22-01-10  prise en compte des bancs decouvrants
! 2010.8  03-05-10  Terminologie et biohz supprimé
! 2010.9  06-06-10  tridiagonalsolver renommé tridiagonalsolver
! 2010.25 08-02-12  kh_w renomme kh_w
! S.26    26-06-14  modif arg appel au solveur tridiag pour numero variables
!                   bio >=2 pour eviter de recalculer la bi-diagonalisation de la
!                   matrice principale
!         20-12-14  C.L. sous le fond
!         21-03-15  call subroutine obc_bio_mpi_fluxbio('bottom','z0')
!         27-03-15  La sedimentation est calculee en meme temps que le melange
!                   pour eviter l'expansion excessive du nepheloide dans le cas
!                   d'un calcul separe
!         23-04-15  suite 27-03-15, modif sur c.l. fond et surface sur tridia_in(:,:,:,2)
!         03-07-15  suite du point precedent qui empechait la sedimentation depuis la
!                   couche de surface. Ajout par consequent de flag_wsed_
!         04-07-15  oter la diffusivite associer a omega ne semble pas etre une
!                   bonne idee. Commentee en attendant de mieux comprendre
!         05-08-15  methode vst continue
!         09-11-16  - wsed peut etre traitee soit de maniere implicite soit explicite,
!                   dans ce dernier cas voir wsed dans l'advection verticale
!                   - dz en facteur de l'equation
!         13-12-16  Choix wsed implicite ou explicite selon les traceurs
!         25-01-17  utilisation plus coherente de wetmask
!         12-02-17  mises A jour NEMO offline
!         23-02-17  Flux air/mer portE par fluxbio_w et non pas omega
! v249    07-03-19  - suite du point precedent: il faut egalement garantir que
!                   la vitesse verticale implicite soit totalement nulle meme
!                   s'il etait peu probable qu'elle puisse etre non nulle
!                   - if(flag_merged_levels==1)call vertmix_merged_levels_bio
! v253    01-05-19  Evitement de la division par dz
! v254    11-05-19  Avec schema wetdry No6 pas besoin de multiplier par
!                   wetmask apres call tridiagonalsolver
! v256    28-05-19  cas modele 1DV multiplier bio par dz
!                   cas couches fusionnees distinguE du cas standard pour calcul du coef 
!                   de melange
!         29-05-19  suite du point precedent: calcul des C.L. fond et surface apres calcul
!                   du cas general
! v303    03-08-21  call my_outputs_zone1bioflux('botsurf',0) !03-08-21
!...............................................................................
!    _________                    .__                  .__                     !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................


!...................................................................!
! NOTE SUR LE SCHEMA IMPLICITE SUR LA VERTICALE:                    !
! On doit résoudre:                                                 !
!         TRIDIA_IN(I,J,K,1)*F(I,J,K-1)+                            !
!         TRIDIA_IN(I,J,K,2)*F(I,J,K  )+                            !
!         TRIDIA_IN(I,J,K,3)*F(I,J,K+1) = TRIDIA_IN(I,J,K,4)        !
! avec  F(K) représentant par exemple par vitesse                   !
! ou bien la température la salinité ... n'importe quelle variable  !
! comportant un schéma de diffusion verticale.                      !
!...................................................................!

      i1_=1
      i2_=imax
      j1_=1
      j2_=jmax

      if(i1dbio.eq.1)then !>>>>>>>>
       i1_=iptbio
       i2_=iptbio
       j1_=jptbio
       j2_=jptbio
      if(ioffline.eq.2)then !------------------------------------------>
       do k=1,kmax+1
           kh_w(iptbio,jptbio,k)=max(    kh_w(iptbio,jptbio,k),un*5.e-5)
       enddo
      else                  !------------------------------------------>
       write(6,*)'mixsed_bio je ne modifie pas difver '
       write(6,*)'car je suis en online.'
       write(6,*)'dans ce cas intervenir directement ds les coeffs de'
       write(6,*)'la matrice - voir tridia_in'
       stop 'mixsedbio'
      endif                 !------------------------------------------>
      endif               !>>>>>>>>

      if(flag3d==0) then !-1d-case-> !28-05-19 v256
! Dans le cas 1dV on ne passe pas par l'advection au sortir de laquelle le tableau
! bio_t est devenu (temporairement) homogene A bio fois dz. On convertit donc bio en bio*dz
! par les lignes suivantes pour l'homogeneitE du calcul de tridia_in(:,:,:,4) A venir
       do vb=1,vbmax ; do k=1,kmax ; do j=j1_,j2_ ; do i=i1_,i2_
         bio_t(i,j,k,vb)=                    &
         bio_t(i,j,k,vb) *0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))
       enddo         ; enddo       ; enddo        ; enddo
      endif              !-1d-case-> !28-05-19


! mpi obc on fluxbio_w(:,:,1)
      call obc_bio_mpi_fluxbio('bottom','z0') ! mpi obc on fluxbio_w(:,:,1) 21-03-15

!     call mixsed_bio_validation(0) ! pour verifier l'equilibre w.dB/dz = Kh.d/dz(dB/dz)

!......................................................
! Kh/dz
      if(flag_merged_levels==1) then !oooo> !28-05-19>

! Cas modele avec couches fusionnEes: (flag_merged_levels=1)
       do j=1,jmax ; do i=1,imax
        do k=2,kmerged_t(i,j)-1
! CommentE le 07-03-19
!       anyv3d(i,j,k,id_kh_over_dz)=100.  & ! melange intense dans la couche fusionnee
!                              /(depth_t(i,j,k)-depth_t(i,j,k-1))
! Adimensionnement de Kh pour que tridia_in (1 et 3) evitant la singularitE du schema precedent
! exposE A une division par quasi-zero si couche fusionnee trop fine. De maniere similaire
! A ce qui est fait pour T,S tridia_in(1 et 3) sont constants avec:
         anyv3d(i,j,k,id_kh_over_dz)=0.5/dti_fw !07-03-19
        enddo


! note sur cas "k=kmerged_t(i,j)" 
! Le cas general est kmerged_t=kmin_w+1 (kmerged est la couche au dessus de kmin_w)
! si cas special kmerged_t=kmin_w=1 la valeur sera ecrasee par la C.L. en k=1 A suivre
! si cas special kmerged_t=kmin_w>2 la ligne suivante augmente le melange entre kmin et la
! couche en dessous (fusionnEe). Le coef de melange est grand car dz(kmerged) est "normal"
        k=kmerged_t(i,j) 
        anyv3d(i,j,k,id_kh_over_dz)=10*dz_t(i,j,k,2)/dti_fw

        do k=kmerged_t(i,j)+1,kmax
         anyv3d(i,j,k,id_kh_over_dz)=kh_w(i,j,k)  &
                                /(depth_t(i,j,k)-depth_t(i,j,k-1))
        enddo
       enddo       ; enddo

      else                           !oooo> !28-05-19>

! Cas modele avec que des couches standards (flag_merged_levels=0)
       do k=2,kmax ; do j=1,jmax ; do i=1,imax
         anyv3d(i,j,k,id_kh_over_dz)=kh_w(i,j,k)  &
                                /(depth_t(i,j,k)-depth_t(i,j,k-1))
       enddo       ; enddo       ; enddo

      endif                          !oooo> !28-05-19>

! Enfin les C.L. fond et surface pour pas de melange au travers du fond et de la surface 
      do j=1,jmax ; do i=1,imax !29-05-19
       anyv3d(i,j,1     ,id_kh_over_dz)=0.
       anyv3d(i,j,kmax+1,id_kh_over_dz)=0.
      enddo       ; enddo
!......................................................
! TEST TENDECO
!     anyv3d(:,:,:,id_kh_over_dz)=0.
!     wsed(1,:)=0.

      do 200 vb=1,vbmax

! The principal matrix does not change if the sedimentation velocity does not change and
! consequently does not need to be recomputed for each variable.
      flag=0 ! reset flag=0 : the matrix does not need to be computed is the default hypothesis

! The principal matrix of the tridiagonal system is computed under the following conditions:
      if(vb==1)flag=1                             ! Condition 1: (initial condition) i.e. compute the matrix
                                                  !              in case of the first variable
      if( wsed(1,vb)/=wsed(1,max(vb-1,1)) )flag=1 ! Condition 2  The sedim velocity is different  !21-04-15
                                                  !              from the previous variable and so the matrix
                                                  !              needs to be computed

      if(flag==1) then   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@>        !13/12/01

! Determiner la vitesse verticale residuelle implicite
      k0=1
      if(vb>vbmax_eco3ms)k0=0 ! Dans MUSTANG le flux est externalise dans fluxbio_w
      do j=1,jmax ; do i=1,imax

! Au fond vitesse verticale nulle sauf sedimentation 100% implicite (hors cas mustang (voir k0))
        k=1
        anyv3d(i,j,k,id_wi_)=k0*wsed(k,vb) !vitesse de sedimentation partie implicite

! Dans la couche fusionnee vitesse verticale 100% implicite
       do k=2,kmerged_t(i,j)
          anyv3d(i,j,k,id_wi_)=omega_w(i,j,k,1)+wsed(k,vb)
       enddo

! Au dessus de la couche fusionnee
       do k=kmerged_t(i,j)+1,kmax


! Cette formule est coherente avec advection_bio
          x0=0.5*min(dz_t(i,j,k  ,0)+dz_t(i,j,k  ,1)            &
                    ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))           &
                    /dti_fw

! omega+wsed residuelle implicite
                   anyv3d(i,j,k,id_wi_)                         & !  Vitesse implicite=
                 =omega_w(i,j,k,1)+wsed(k,vb)                   & !  Vitesse totale
         -max(min(omega_w(i,j,k,1)+wsed(k,vb)*wsed_explicit(vb) & !  moins la vitesse verticale explicite
                 , substep_advbio*x0)                           & !  telle que definie dans advection_bio
                 ,-substep_advbio*x0)*wetmask_t(i,j)

! bidouille patrick
!                  anyv3d(i,j,k,id_wi_)=omega_w(i,j,k,1)+wsed(k,vb)

       enddo     ! boucle k

! Surface:
!     k=kmax+1
!         x0=0.5*min(dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1)   & ! coherent avec advection_bio
!                   ,dz_t(i,j,k-1,0)+dz_t(i,j,k-1,1))  &
!                   /dti_fw

! omega+wsed residuelle implicite
!                  anyv3d(i,j,k,id_wi_)                         & !  Vitesse implicite=
!                =omega_w(i,j,k,1)+wsed(k,vb)                   & !  Vitesse totale
!        -max(min(omega_w(i,j,k,1)+wsed(k,vb)*wsed_explicit(vb) & !  moins la vitesse verticale explicite
!                , substep_advbio*x0)                           & !  telle que definie dans advection_bio
!                ,-substep_advbio*x0)*wetmask_t(i,j)

! Flux de surface imposEs dans fluxbio. Flux "advectif nuls" => vitesse omega explicite et implicite toutes 2 nulles
          anyv3d(i,j,kmax+1,id_wi_)=0. !07-03-19


      enddo ; enddo    ! boucles i,j

      do k=1,kmax
      do j=j1_,j2_
      do i=i1_,i2_

! coef1 --->
          tridia_in(i,j,k,1)=dti_fw*mask_t(i,j,kmax)*( & !--->
!               -kh_w(i,j,k    )/(depth_t(i,j,k  )- depth_t(i,j,k-1)) &
                   -anyv3d(i,j,k,id_kh_over_dz)                       &
              -0.5*(anyv3d(i,j,k,id_wi_)+abs(anyv3d(i,j,k,id_wi_))) &

                                                  )   !--->

! coef3 --->
          tridia_in(i,j,k,3)=dti_fw*mask_t(i,j,kmax)*( & !--->
!               -kh_w(i,j,k+1  )/(depth_t(i,j,k+1)- depth_t(i,j,k  )) &
                   -anyv3d(i,j,k+1,id_kh_over_dz)                     &
              +0.5*(anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))&
                                                  )   !--->

       tridia_in(i,j,k,2)=0.5*(dz_t(i,j,k,1)+dz_t(i,j,k,2))           &
                            +dti_fw*mask_t(i,j,kmax)*( & !--->
!                kh_w(i,j,k    )/(depth_t(i,j,k  )- depth_t(i,j,k-1)) &
!               +kh_w(i,j,k+1  )/(depth_t(i,j,k+1)- depth_t(i,j,k  )) &
                    anyv3d(i,j,k  ,id_kh_over_dz)                     &
                   +anyv3d(i,j,k+1,id_kh_over_dz)                     &
              +0.5*( & !)))>
                     (anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
                    -(anyv3d(i,j,k  ,id_wi_)-abs(anyv3d(i,j,k  ,id_wi_))) &
                   ) & !)))>
                                                   )   !--->

! coef4 --->
! tendancebio est multiplie par max(bio,0)/max(bio,small) qui l'annule si bio<0
!     tridia_in(i,j,k,4)=(bio_t(i,j,k,vb)                     &
!         +dti_fw*tendancebio_t(i,j,k,vb)                     &
!                      )*mask_t(i,j,kmax)                        &    !28/02/03
!                    *0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))
! Avec la modif du schema d'advection, A ce stade bio_t est homogene A
! bio_t*dz_t de sorte que la multiplication par dz_t est restreinte A
! tendancebio_t: !01-05-19
!      tendancebio_t(i,j,k,vb)=1./86400.
      tridia_in(i,j,k,4)= bio_t(i,j,k,vb)                     &
          +dti_fw*tendancebio_t(i,j,k,vb)                     &
                        *mask_t(i,j,kmax)                        &    !28/02/03
                     *0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))

      enddo
      enddo
      enddo

! fluxbio_w(:,:,vb,:)=0. ! removing fluxes from 1 (bottom) and 2 (surface)


!......................................................
! Bottom & surface boundary conditions:
      do j=j1_,j2_
      do i=i1_,i2_

! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
! Principe de la modification du coef 2:
! En surface b(k+1)=b(k) --> a2=a2+a3 suivi de a3=0
! Au fond b(k-1)=b(k)    --> a2=a2+a1 suivi de a1=0
! CommentE le !23-02-17
! fond:
!      tridia_in(i,j,kmin_w(i,j),2)=tridia_in(i,j,kmin_w(i,j),2)  & !23-04-15
!                                  +tridia_in(i,j,kmin_w(i,j),1)
!      tridia_in(i,j,kmin_w(i,j),1)=0.

! CommentE le !23-02-17
! surface:
!      tridia_in(i,j,kmax,2)=tridia_in(i,j,kmax,2)+tridia_in(i,j,kmax,3)
!      tridia_in(i,j,kmax,3)=0.

       tridia_in(i,j,kmin_w(i,j),4)=tridia_in(i,j,kmin_w(i,j),4)  &
                            +dti_fw*fluxbio_w(i,j,vb,1)

       tridia_in(i,j,kmax,4)=                                  &
       tridia_in(i,j,kmax,4)+dti_fw*fluxbio_w(i,j,vb,2)

      enddo
      enddo

      call tridiagonalsolver(1,i1_,i2_,j1_,j2_,kmax)           !22/06/06


      else             !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@>        !13/12/01


      do k=1,kmax
      do j=j1_,j2_
      do i=i1_,i2_

! coef4 --->
!     tridia_in(i,j,k,4)=(bio_t(i,j,k,vb)                     &
!         +dti_fw*tendancebio_t(i,j,k,vb)                     &
!                      )*mask_t(i,j,kmax)                        &    !28/02/03
!                    *0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))
! Avec la modif du schema d'advection, A ce stade bio_t est homogene A
! bio_t*dz_t de sorte que la multiplication par dz_t est restreinte A
! tendancebio_t:
!      tendancebio_t(i,j,k,vb)=1./86400.
      tridia_in(i,j,k,4)= bio_t(i,j,k,vb)                     &
          +dti_fw*tendancebio_t(i,j,k,vb)                     &
                        *mask_t(i,j,kmax)                     &    !28/02/03
                     *0.5*(dz_t(i,j,k,0)+dz_t(i,j,k,1))

      enddo
      enddo
      enddo

!......................................................                !28/02/03
! LES CONDITIONS AUX LIMITES SUR Coef4:
      do j=j1_,j2_                                                      !28/02/03
      do i=i1_,i2_                                                      !28/02/03

! fond:
       tridia_in(i,j,kmin_w(i,j),4)=                                   &
       tridia_in(i,j,kmin_w(i,j),4)+dti_fw*fluxbio_w(i,j,vb,1)
! surface:
       tridia_in(i,j,kmax       ,4)=                                  & !28/02/03
       tridia_in(i,j,kmax       ,4)+dti_fw*fluxbio_w(i,j,vb,2)
      enddo                                                            !28/02/03
      enddo                                                            !28/02/03
!......................................................

!     call tridiagonalsolver(1,i1_,i2_,j1_,j2_,kmax)          !22/06/06
      call tridiagonalsolver(0,i1_,i2_,j1_,j2_,kmax)          !26-06-14

      endif            !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@>        !13/12/01

!     if(abs(wsed(vb)).lt.small1) then!*******************************>!13/12/01
! si vitesse de sedimentation est nulle on conclut ici le calcul sinon...

!     do k=1,kmax ; do j=j1_,j2_ ; do i=i1_,i2_
      do k=kmax,1,-1 ; do j=j1_,j2_ ; do i=i1_,i2_

! peut etre..      bio_t(i,j,k,vb)=tridia_out(i,j,k)*mask_t(i,j,kmax) !25-01-17
!       bio_t(i,j,k,vb)=tridia_out(i,j,k)*wetmask_t(i,j)*mask_t(i,j,kmax)
        bio_t(i,j,k,vb)=tridia_out(i,j,k)*mask_t(i,j,kmax) !11-05-19

      enddo       ; enddo        ; enddo

      if(flag_nemoffline==1) then !-nemocase-> !11-02-17

! Under the bottom:
       do j=j1_,j2_
       do i=i1_,i2_
       do k=1,kmin_w(i,j)-1
        bio_t(i,j,k,vb)=bio_t(i,j,kmin_w(i,j),vb) !20-12-14
       enddo
       enddo
       enddo

      endif                      !-nemocase->


! sedimentation implicite:
!     else                            !*******************************>!13/12/01
! sinon...si la vitesse de sedimentation est non nulle on poursuit par
! un calcul de type implicite:
!     const1=dti_fw*wsed(vb)
!     k=kmax
!     do j=j1_,j2_
!     do i=i1_,i2_
!      bio_t(i,j,k,vb)=wetmask_t(i,j)*mask_t(i,j,kmax)           &      !22-01-10
!                                *tridia_out(i,j,k)           &      !attention piege: ne pas prendre BIOHZ(....2)
!                           /(1.-const1/dz_t(i,j,k,1))
!     enddo
!     enddo
!     do j=j1_,j2_
!     do i=i1_,i2_
!      do k=kmax-1,kmin_w(i,j),-1 !20-12-14
!       bio_t(i,j,k  ,vb)=wetmask_t(i,j)*mask_t(i,j,kmax)               &
!       *(tridia_out(i,j,k)-bio_t(i,j,k+1,vb)*const1/dz_t(i,j,k,1))  &
!       /(1.                                 -const1/dz_t(i,j,k,1))
!      enddo
!      do k=1,kmin_w(i,j)-1
!       bio_t(i,j,k,vb)=bio_t(i,j,kmin_w(i,j),vb) !20-12-14
!      enddo
!     enddo
!     enddo
!     endif                           !*******************************>!13/12/01

  200 continue

!     call mixsed_bio_validation(1)     ! pour verifier l'equilibre w.dB/dz = Kh.d/dz(dB/dz)

#ifdef bilanbio
! appel placE ici pour pouvoir utiliser la valeur de bio_t de la sedimentation
! pour le calcul du flux correspondant dans le bilan:
     call my_outputs_zone1bioflux('botsurf',0) !03-08-21
     call my_outputs_zone2bioflux('botsurf',0) !03-08-21
     call my_outputs_zone3bioflux('botsurf',0) !03-08-21
     call my_outputs_zone4bioflux('botsurf',0) !03-08-21
#endif

! Homogeneisation forcEe de la couche fusionnee
      if(flag_merged_levels==1)call vertmix_merged_levels_bio !07-03-19

      end subroutine mixsed_bio

!......................................................................
#ifdef bidon
      subroutine mixsed_bio_validation(case_)
      use module_principal
      use module_parallele
      implicit none
      integer case_

! Cette routine sert a verifier le codage du melange et de la chute.
! On verifie que, dans le cas ou Kz et wsed sont constant, la solution B=B0.exp( z.wsed/Kz )
! stationnaire satisfait l'equilibre w.dB/dz=K.d/dz(dB/dz).
! Pour le cas particulier du premier niveau, contourner la difficultée en fixant le fluxbio du fond à
! fluxbio=K.dB/dz=K.(w/K).B

      if(case_==0) then !0000000>

! imposer un coef constant
       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
           kh_w(i,j,k)=0.01
        omega_w(i,j,k,1)=0.
       enddo ; enddo ; enddo

! imposer un profil analytique stationnaire
       do vb=1,vbmax
       if(wsed(vb)/=0.) then
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          bio_t(i,j,k,vb)=exp(depth_t(i,j,k)*wsed(vb)/kh_w(i,j,k))
         enddo ; enddo ; enddo

       endif
       enddo

      endif             !0000000>

      if(case_==1) then !1111111>

! Re-imposer le profil analytique stationnaire
! puisque celui-ci a ete modifie par le calcul...
       do vb=1,vbmax
       if(wsed(vb)/=0.) then
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          bio_t(i,j,k,vb)=exp(depth_t(i,j,k)*wsed(vb)/kh_w(i,j,k))
         enddo ; enddo ; enddo

       endif
       enddo

!     i=imax/2 ; j=jmax/2 ; vb=1
      i=74 ; j=36 ; vb=1
      write(6,*)'mask=',mask_t(i,j,kmax),h_w(i,j)
      do k=kmax-1,2,-1
       write(6,*)                                      &
               k,real(depth_t(i,j,k))                  &
                ,real(bio_t(i,j,k,vb))                 &
           ,real(tridia_in(i,j,k,1)*bio_t(i,j,k-1,vb)  &
                +tridia_in(i,j,k,2)*bio_t(i,j,k  ,vb)  &
                +tridia_in(i,j,k,3)*bio_t(i,j,k+1,vb)) &
                ,real(tridia_in(i,j,k,4))
      enddo
      stop 'jtbac'
      endif             !1111111>


      end subroutine mixsed_bio_validation
#endif
