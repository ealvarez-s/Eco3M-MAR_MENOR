      subroutine mix_sig_step(ichoix)
!______________________________________________________________________
!
! S model
! release 2010.25  - last update: 28-03-12
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      use module_global
      implicit none
      integer ichoix

!...............................................................................
!Version  Date      Description des modifications
!         14/11/02: mise en service
!         08/11/05: suite
!         10/07/03: suite (SMALL3)
!         02/08/03: bienvenue à NHYBSIG
!         08/10/03: debug i/o fichiers FORM='UNFORMATTED'
!         10/10/03: debug i/o fichiers FORM='UNFORMATTED' suite
!         24/03/04: lissage des niveaux z
!         24/03/04: sigma_z_hyb archivé de K=0,NR
!         26/03/04: c'est finalement  depth_z qui est archivé à la place de sigma
!                   car ca ne sert à rien d'avoir le nouveau sigma sans la
!                   nouvelle bathy ajustee à la nouvelle profondeur du fond
!                   lissée.
!         10/05/04: debug grille pour calcul gradient possible dans obc
!         17/06/04: suite du point precedent:
!         14/04/05: z(fond) adapté à la bathy et non l'inverse
!         21/04/05: dans la mise à jour precedente, seul le premier niveau au
!                   dessus du fond s'adapte. Maintenant on a la possibilité
!                   d'adapter plusieur niveau même si par défaut on reste sur
!                   l'option "1" niveau
!         25/07/05: suite du point 10/05/04 concernant les obc pour pouvoir
!                   calculer des conditions de gradient. Avant on adaptait la
!                   grille à la bathy car on pensait au respect du tranport des
!                   forcage. Or si on y reflechit bien, le transport indicident
!                   qui compte dans la condition de Flather est en i=2 et
!                   i=1 ce qui veut dire que l'on peut adapter la bathy sans
!                   modifier h_y(i=2) et donc c'est mieux car modifier la grille
!                   peut avoir de facheuses consequence comme genrer de forts
!                   gradient horizontaux de densité artificiels.
!         29/07/05: x6 remplace const1 pour deconflictualiser real et double
!                   precision
!         30/09/05: modif pour versions r8 possibles
!         19/11/05: la condition aux limites n'allait encore pas....
!         26/01/06: ANYV3D remplace ANYVAR3D
!         22/03/06: modif sur ecriture lecture fichier profondeur pour usage
!                   indifferent simple ou double precision. En fin de ce present
!                   fichier figure quelques lignes de fortran pour convertir les
!                   anciens fichiers au nouveau mode de lecture
!         10/05/06: fonctions compatibles avec double precision
!         23/10/03: Modif pour garantir que sigma(kmax+1)=1 en toutes
!                   circonstances....
!         21/03/07: Ajout de garde-fous et comentaires...
!         31/03/08: possibilite d'une grille hybride sans procedure iterative
!                   et sans fichier
!         16/10/08: C.L. sur   mask_x et   mask_y pour advection_ts.F
!         11/03/09: Amenagements pour parallelisation
! 2009.2  23-09-09: - mise à jour de la fonction gather
!                   - des boucles à la place de "array syntax"
!                   - contrainte dz(k) > dz(k+1)
!         29-09-09: - contrainte x1*dz(k+1) < dz(k) < x2*dz(k+1)
!                   - z s'ajuste à la vraie bathy sur un seul niveau
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
!         08-10-09: compatibilité avec f95: adequation des longueurs de
!                   chaines de charactere passées en argument de allocate_global
!         19-10-09: Afin d'eviter des incoherence entre la grille précalculée et
!                   d'eventuelle fausses manip sur le fichier de bathy en entree
!                   le niveau du fond est ajusté à la bathy dès la phase iterative
!                   et non pas apres la phase de lecture
!         21-10-09: routine z_to_xyr renommée maskz_to_maskxyr
! 2010.12 20-09-10  Possibilité de calcul en simple precision
! 2010.15 31-12-10  Evolution de la grille sigma-step (transition continue des marches)
!                   par initialisation de kbbc_
!         06-01-11  - Les 10 dernieres iterations ne servent qu'à lisser le maillage
!                   - Un nouveau parametrage de x1 x2 favorise l'epaississement des
!                   couches avec la profondeur
!                   avec la profondeur
! 2010.16 12-01-11  Ajuster la bathy pour que la premiere couche ne soit pas trop fine
! 2010.17 13-01-11  Ajout d'une C.L. de type "z1" sur kbbc_w
! 2010.18 25-01-11  Modif C.L. laterales & fond, possibilité de rappel
!                   vers grille de reference (par defaut rappel nul)
!         29-01-11  Bornage modifié apres experience de dz trop petit
!         30-01-11  Reset k2
!         19-02-11  Initialisation de kbbc_ depuis routine bottomboundary
!         23-02-11  Retour en arriere par rapport à kbbc autrement commentaire devant
!                   call bottomboundary('init kbbc  ')
! 2010.25 28-03-12  echange compatible avec pgf Stelios
!...............................................................................

      if(ichoix.eq.0)return
      if(i2dh.eq.0) then                                               !21/03/07
      write(6,*)
      write(6,*)'i2dh=0 ne pas passer par mix_sig_step'
      write(6,*)                                                        &
       'dans notebook_vertcoord oter l''option syst hybride'
      stop ' donc stop en attendant! '
      endif
#ifdef synopsis
       subroutinetitle='mix_sig_step'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     I=132
!     X1=1. ! facteur multiplicatif
!     OPEN(UNIT=3,FILE='avant.out',RECL=90000)
!     WRITE(3,*)'Profondeurs sens indice j croissant'
!     WRITE(3,*)'i est constant et égal à ',I
!     WRITE(3,'(300(I4,1X))')(J,J=1,NECO)
!     WRITE(3,*)
!     DO K=NR,1,-1
!     WRITE(3,'(300(I4,1X))')(
!    1  INT(ABS(PROFWK_Z(I,J,K)*X1             )),J=1,NECO)
!     ENDDO
!     CLOSE(3)

! Possibilite d'une grille hybride sans procedure iterative
! et sans fichier:
      if(nhybsig.eq.0) goto 2222                                       !31/03/08

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE ITERATIVE
! DEBUT:
      if(ichoix.eq.1) then
!_______________________________________________________________________________


!     i=imax/2 ; j=jmax/2
!     do k=kmaxp1,1,-1
!     write(*,*)depth_w(i,j,k),(sigma_w(i,j,k)-1.)*h_w(i,j)
!     enddo
!     stop 'coco'

      k2=2                                                              !23-09-09
      x5=0.125                                                          !23-09-09

      do 100 k1=1,nhybsig                                               !02/08/03

! condition aux limites:
      do j=1,jmax
      do i=1,imax
!       depth_w(i,j,0)=-hmax
        depth_w(i,j,0)=2.*depth_w(i,j,1)-depth_w(i,j,2)    !25-01-11
      enddo
      enddo

      x7=1.    ! incoherence hydro max
!     x8=0.01 ! rappel vers grille de reference
!     x8=0.001 ! rappel vers grille de reference
      x8=0.0 ! rappel vers grille de reference

!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax

!     anyv3d(i,j,k,2)=                                               &
!           abs(  depth_w(i  ,j,k+1)+ depth_w(i  ,j,k)               &
!               - depth_w(i,j-1,k+1)- depth_w(i,j-1,k) )             &
!         /(abs(  depth_w(i  ,j,k+1)- depth_w(i  ,j,k)               &
!               + depth_w(i,j-1,k+1)- depth_w(i,j-1,k))              &
!               +small1 )*mask_t(i,j-1,kmaxp1)*mask_t(i,j,kmaxp1)

!     anyv3d(i,j,k,3)=                                               &
!           abs(  depth_w(i  ,j,k+1)+ depth_w(i  ,j,k)               &
!               - depth_w(i-1,j,k+1)- depth_w(i-1,j,k) )             &
!         /(abs(  depth_w(i  ,j,k+1)- depth_w(i  ,j,k)               &
!               + depth_w(i-1,j,k+1)- depth_w(i-1,j,k))              &
!               +small1 )*mask_t(i-1,j,kmaxp1)*mask_t(i,j,kmaxp1)

!     enddo
!     enddo
!     enddo

      do 10 j=1,jmax
      do 10 i=1,imax
      if(mask_t(i,j,kmaxp1)==1)then !mskmskmsk>

      do k=1,kmax
!     km1=max(k-1,1)
      km1=    k-1

! Ces nouvelles lignes supposent que l'on calcule des C.L. en imax+1,0,jmax+1,0 !11/03/09
      jm1=j-1
      jp1=j+1
      im1=i-1
      ip1=i+1

      if( depth_w(i,j,k).ge. depth_w(i,j,k+1)) then
      write(6,*)i,j,k, depth_w(i,j,k), depth_w(i,j,k+1)
      stop 'erreur croisement de niveaux'
      endif


      x0=(0.5*h_w(i,j)/real(kmax+1))**2
      x10=h_w(i,j)

      x1=                                                               &
         min(                                                           &
         max(                                                           &
         max(                                                           &
         ( depth_w(i,j,k)- depth_w(im1,j,k))   & ! la pente des z ne doit pas etre
        *( depth_w(i,j,k)- depth_w(im1,j,k)+h_w(i,j)-h_w(im1,j))/x0 & ! plus
       ,                                        & ! grande que la pente des H
         ( depth_w(i,j,k)- depth_w(im1,j,k+1)) & ! inconsistence hydrostatique
        *( depth_w(i,j,k)- depth_w(im1,j,km1))/x0                       &
!        0.5*(anyv3d(i,j,k,3)+anyv3d(i,j,km1,3))-x7           &
            )                                                 &
            ,zero)                                            &
!           ,x5)
            ,x5)*mask_t(im1,j,kmaxp1)     !25-01-11

!     x1=min(max(0.5*(anyv3d(i,j,k,3)+anyv3d(i,j,km1,3))-2.,zero),x5)

      x2=                                                               &
         min(                                                           &
         max(                                                           &
         max(                                                           &
         ( depth_w(i,j,k)- depth_w(ip1,j,k))                            &
        *( depth_w(i,j,k)- depth_w(ip1,j,k)+h_w(i,j)-h_w(ip1,j))/x0     &
       ,                                                                &
         ( depth_w(i,j,k)- depth_w(ip1,j,k+1))                          &
        *( depth_w(i,j,k)- depth_w(ip1,j,km1))/x0                       &
!        0.5*(anyv3d(i+1,j,k,3)+anyv3d(i+1,j,km1,3))-x7      &
            )                                                           &
            ,zero)                                                      &
!           ,x5)
            ,x5)*mask_t(ip1,j,kmaxp1)

!     x2=min(max(0.5*(anyv3d(i+1,j,k,3)+anyv3d(i+1,j,km1,3))-2.,zero),x5)

      x3=                                                               &
         min(                                                           &
         max(                                                           &
         max(                                                           &
         ( depth_w(i,j,k)- depth_w(i,jm1,k))                            &
        *( depth_w(i,j,k)- depth_w(i,jm1,k)+h_w(i,j)-h_w(i,jm1))/x0     &
       ,                                                                &
         ( depth_w(i,j,k)- depth_w(i,jm1,k+1))                          &
        *( depth_w(i,j,k)- depth_w(i,jm1,km1))/x0                       &
!        0.5*(anyv3d(i,j,k,2)+anyv3d(i,j,km1,2))-x7                  &
            )                                                           &
            ,zero)                                                      &
!           ,x5)
            ,x5)*mask_t(i,jm1,kmaxp1)

!     x3=min(max(0.5*(anyv3d(i,j,k,2)+anyv3d(i,j,km1,2))-2.,zero),x5)

      x4=                                                               &
         min(                                                           &
         max(                                                           &
         max(                                                           &
         ( depth_w(i,j,k)- depth_w(i,jp1,k))                            &
        *( depth_w(i,j,k)- depth_w(i,jp1,k)+h_w(i,j)-h_w(i,jp1))/x0     &
       ,                                                                &
         ( depth_w(i,j,k)- depth_w(i,jp1,k+1))                          &
        *( depth_w(i,j,k)- depth_w(i,jp1,km1))/x0                       &
!        0.5*(anyv3d(i,j+1,k,2)+anyv3d(i,j+1,km1,2))-x7         &
            )                                                           &
            ,zero)                                                      &
!           ,x5)
            ,x5)*mask_t(i,jp1,kmaxp1)

!     x4=min(max(0.5*(anyv3d(i,j+1,k,2)+anyv3d(i,j+1,km1,2))-2.,zero),x5)


! Les 10 dernieres iterations ne servent qu'à lisser le maillage !06-01-11
        if(k1>nhybsig-4) then
         x1=0.125*mask_t(i,j,kmaxp1)*mask_t(im1,j  ,kmaxp1)   !25-01-11
         x2=0.125*mask_t(i,j,kmaxp1)*mask_t(ip1,j  ,kmaxp1)
         x3=0.125*mask_t(i,j,kmaxp1)*mask_t(i  ,jm1,kmaxp1)
         x4=0.125*mask_t(i,j,kmaxp1)*mask_t(i  ,jp1,kmaxp1)
         x8=0.
        endif

        anyv3d(i,j,k,1)=                                                &
         depth_w(i  ,j  ,k)*(1.-x1-x2-x3-x4)                            &
       + depth_w(im1,j  ,k)*x1                                          &
       + depth_w(ip1,j  ,k)*x2                                          &
       + depth_w(i  ,jm1,k)*x3                                          &
       + depth_w(i  ,jp1,k)*x4

        anyv3d(i,j,k,1)=(1.-x8)*anyv3d(i,j,k,1)     &        !25-01-11
                           +x8*(sigma_w(i,j,k)-1.)*h_w(i,j)

      enddo
      endif                         !mskmskmsk>
   10 continue

!  12 continue

      do 11 k=1,kmax
      do 11 j=1,jmax
      do 11 i=1,imax
       depth_w(i,j,k)=anyv3d(i,j,k,1)
   11 continue

! cas particulier du fond:
!     do j=1,jmax
!     do i=1,imax
!      depth_w(i,j,1)=max( min( depth_w(i,j,1),-h_w(i,j)) , -hmax )    !26/03/04
!     enddo
!     enddo

! cas particulier des niveaux trop proches les uns ds autres
      x6=0.01        !29/07/05 !29-01-11
      x1=-2. ; x2=-0.2
!     k=kmax
      do k=kmax,1,-1  !29-01-11
      do j=1,jmax
      do i=1,imax

       depth_w(i,j,k)= depth_w(i,j,k+1)                                 &
                -max(  depth_w(i,j,k+1)                                 &
                     - depth_w(i,j,k  ),x6 ) ! separation mini: 20cm
!      depth_w(i,j,k)=min(max(depth_w(i,j,k),x1),x2)

      enddo
      enddo
      enddo

! Assure que x1*dz(k+1) < dz(k) < x2*dz(k+1)
!     x1=0.5 ; x2=2.
      x1=0.1   ; x2=2000.                                                !06-01-11
!     x1=0.1    ; x2=2.                                                   !06-01-11
      do j=1,jmax
      do i=1,imax

      do k=kmax-1,1,-1                                                 !29-09-09

         depth_w(i,j,k  )=max(min(                                      &
         depth_w(i,j,k  )                                               &
       , depth_w(i,j,k+1)-x1*( depth_w(i,j,k+2)- depth_w(i,j,k+1)))     &
       , depth_w(i,j,k+1)-x2*( depth_w(i,j,k+2)- depth_w(i,j,k+1)))

      enddo

      if(depth_w(i,j,1)>-h_w(i,j)) then
      rap=-h_w(i,j)/depth_w(i,j,1)
        do k=kmax,1,-1
         depth_w(i,j,k)=rap*depth_w(i,j,k)
        enddo
      endif

      enddo
      enddo


!................................................!                     !11/03/09
! Conditions aux limites laterales pour PROFWK_Z
      call obc_depth(0)
!................................................!

      if(mod(k1,10)==0) then !xxxxxxxxxxxxxxx>
!     if(mod(k1,1)==0) then !xxxxxxxxxxxxxxx>
! bilan d'incoherence hydrostatique:

      x21=0.
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=2,jmax
      do i=2,imax

      if( depth_w(i,j  ,k)>-h_w(i,j  ).and.                             &
          depth_w(i,j-1,k)>-h_w(i,j-1))then !jjjjjjjjjjj>

      x22=                                                              &
            abs(  depth_w(i  ,j,k+1)+ depth_w(i  ,j,k)                  &
                - depth_w(i,j-1,k+1)- depth_w(i,j-1,k) )                &
          /(abs(  depth_w(i  ,j,k+1)- depth_w(i  ,j,k)                  &
                + depth_w(i,j-1,k+1)- depth_w(i,j-1,k))                 &
                +small1 )*mask_t(i,j-1,kmaxp1)*mask_t(i,j,kmaxp1)

       sum1=sum1+1.
       sum2=sum2+x22

       if(x22>x21)THen
        x21=x22
        i3=i
        j3=j
        k3=k
       endif

      endif                                 !jjjjjjjjjjj>

      if( depth_w(i  ,j,k)>-h_w(i  ,j).and.                             &
          depth_w(i-1,j,k)>-h_w(i-1,j))then !iiiiiiiiiii>

      x22=                                                              &
            abs(  depth_w(i  ,j,k+1)+ depth_w(i  ,j,k)                  &
                - depth_w(i-1,j,k+1)- depth_w(i-1,j,k) )                &
          /(abs(  depth_w(i  ,j,k+1)- depth_w(i  ,j,k)                  &
                + depth_w(i-1,j,k+1)- depth_w(i-1,j,k))                 &
                +small1 )*mask_t(i-1,j,kmaxp1)*mask_t(i,j,kmaxp1)

       sum1=sum1+1.
       sum2=sum2+x22

       if(x22>x21)THen
        x21=x22
        i3=i
        j3=j
        k3=k
       endif

      endif                                 !iiiiiiiiiii>

      enddo
      enddo
      enddo

      write(*,*)'Iteration:',k1
      write(*,*)'hydrostatic inconsistency. Mean value=',sum2/sum1
      write(*,*)'hydrostatic inconsistency.  Max value=',x21
      write(*,*)'and related cell box=',i3,j3,k3
       open(unit=3,file='iteration_grille.out')
        write(3,*)k1
       close(3)


      endif                  !xxxxxxxxxxxxxxx>

      if(mod(k1,50)==0) then
!     call sigma_levels_figure
!     pause
      endif

  100 continue
  101 continue


! Un dessin pour expliquer....

!     -----------------------------      Niveau entier sigma  8 ->Surface Ocean
!     - - - - - - - - - - - - - - -      1/2 niveau            7
!     -----------------------------      Niveau entier sigma  7
!     - - - - - - - - - - - - - - -      1/2 niveau            6
!     -----------------------------      Niveau entier sigma  6
!     - - - - - - - - - - - - - - -      1/2 niveau            5
!     -----------------------------      Niveau entier sigma  5
!     - - - - - - - - - - - - - - -      1/2 niveau            4
!     -----------------------------      Niveau entier sigma  4
!     - - - - - - - - - - - - - - -      1/2 niveau            3
!     -----------------------------      Niveau entier sigma  3  --> Fond Ocean
!     - - - - - - - - - - - - - - -      1/2 niveau            2 --> inutilisé
!     -----------------------------      Niveau entier sigma  2  --> inutilisé
!     - - - - - - - - - - - - - - -      1/2 niveau            1 --> inutilisé
!     -----------------------------      Niveau entier sigma  1  --> inutilisé

! exemple: les 2 premiers niveaux verticaux sont "sous le fond" de l'ocean.
! Cela implique:
!
! mask(1)=0
! mask(2)=0
! mask(3)=1
! ...
! mask(7)=1
! mask(8)=1

! Les niveaux sigma entier sont recalculés (normalisation) tels que
! sigma(1) et sigma(2) ne représentent plus rien (<0), sigma(3)=0 ... sigma(8)=1
! le fond et la surface coincideront respectivement avec sigma(3) & sigma(8).
! Les boucles verticales pour T,S,u,v iront de 3 à 7 alors que pour omega, tke
! etc... on ira de 3 à 8 avec omega(3)=omega(8)=0.
! Bref, le debut de boucle verticale est 3 pour toutes les variables.


!..............................................
! masquage:
!..............................................
! On note que si PROFWK_Z et H_Z ont fait l'objet d'un echange
! pour la parallelisation, alors kmin_z en beneficie de facto
!     do j=0,jmax+1
!     do i=0,imax+1
!      dist1=1.e10
!      k2=1
!      do k=1,kmax
!       dist2=abs(  depth_w(i,j,k) + h_w(i,j) )
!       if(dist2.lt.dist1) then
!         dist1=dist2
!         k2=k
!       endif
!      enddo
!     kmin_w(i,j)=min0(k2,kmax)
!     enddo
!     enddo
      do j=0,jmax+1
      do i=0,imax+1
      k2=1             !30-01-11
       do k=1,kmax
        if(depth_w(i,j,k)<-h_w(i,j).and.depth_w(i,j,k+1)>=-h_w(i,j)) then
!        if(0.1*abs(depth_w(i,j,k+1)+h_w(i,j))    &
!              <abs(depth_w(i,j,k  )+h_w(i,j))) then
!          k2=k+1
!        else
!          k2=k
!        endif
        k2=k
        endif
       kmin_w(i,j)=min0(k2,kmax)
       enddo
      enddo
      enddo

!..............................................
! The lowest level is adjusted to the bathymetry:                      !19-10-09
      do j=0,jmax+1
      do i=0,imax+1
        depth_w(i,j,kmin_w(i,j))=-h_w(i,j)
      enddo
      enddo

!..............................................
! La condition limite au fond = niveaux "quasi collés"
!..............................................
!     small3=1.e-3
      small3=1.e-5
      do j=0,jmax+1
      do i=0,imax+1
      do k=0,kmin_w(i,j)-1
        depth_w(i,j,k)= depth_w(i,j,kmin_w(i,j))+(k-kmin_w(i,j))*small3
      enddo
      enddo
      enddo

!     i=imax/2 ; j=20
!     do k=kmax,kmin_w(i,j),-1
!     write(*,*)k,depth_w(i,j,k)
!     enddo

! Optionnel: le lissage des niveaux z
!.........................................................
! etape 2: lissage des niveaux z                                       !24/03/04
!     CALL LISSEBATHY(2              ! choix n°2
!    &               ,0.             ! inutile
!    &               ,1,MECO,1,NECO  ! zone horizontale de lissage
!    &               ,20)            ! nbre de lissages
!.........................................................

!..................................................................!
! ECRITURE DES SIGMA ET KLOW APRES CONVERGENCE DU CALCUL

! Ecriture par plan via tableau simple precision anyvar2D afin         !22/03/06
! de travailler avec un fichier simple precision quelque soit
! la version du model_

!ECRITURE DANS DOMAINE 0
      call allocate_global('a','glob_anyvar2d       ',iglb,jglb,0)      !08-10-09
      call allocate_global('a','glob_kmin_w         ',iglb,jglb,0)      !08-10-09

      if (par%rank .eq. 0) then !----------------->
      open(unit=3                                                       &
       ,file=directory(1:lname1)//txtslash//'z_sig_step.out'            &
       ,access='direct'                                                 &
       ,form='unformatted'                                              &
       ,recl=(iglb+2)*(jglb+2)*4)
      endif                     !----------------->

      do k=0,kmax+1  ! debut boucle K

       do j=0,jmax+1
       do i=0,imax+1
!            anyvar2d(i             ,j             )= depth_z(i,j,k)
        glob_anyvar2d(i+par%timax(1),j+par%tjmax(1))= depth_w(i,j,k)    !23-09-09
       enddo
       enddo

#ifdef parallele
! Assemblage. ind indique les indices min et max de l'assemblage        !#mpi
      ind(1)=par%timax(1)+0 ;  ind(2)=par%timax(1)+imax+1               !23-09-09
      ind(3)=par%tjmax(1)+0 ;  ind(4)=par%tjmax(1)+jmax+1
      lb2=lbound(glob_anyvar2d) ; ub2=ubound(glob_anyvar2d)
      call  par_gatherall_2d(glob_anyvar2d,lb2,ub2,ind,par%nbdom)
#endif

      if (par%rank==0) then                                             !#mpi
      write(3,rec=k+1)((glob_anyvar2d(i,j),i=0,iglb+1),j=0,jglb+1)
      endif

      enddo      ! fin boucle K

      if (par%rank==0) close(3)                                         !#mpi

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out) ! assure la synchro, permet au domaine 0 de rattraper les autres
#endif

      texte90=directory(1:lname1)//txtslash//'kmin_w_hyb.out' !

       do j=0,jmax+1
       do i=0,imax+1
        glob_kmin_w(i+par%timax(1),j+par%tjmax(1))=kmin_w(i,j)          !23-09-09
       enddo
       enddo

#ifdef parallele
! Assemblage. ind indique les indices min et max de l'assemblage        !#mpi
      ind(1)=par%timax(1)+0 ;  ind(2)=par%timax(1)+imax+1               !23-09-09
      ind(3)=par%tjmax(1)+0 ;  ind(4)=par%tjmax(1)+jmax+1
      lb2=lbound(glob_kmin_w) ; ub2=ubound(glob_kmin_w)
      call  par_gatherall_2d(glob_kmin_w,lb2,ub2,ind,par%nbdom)
#endif

      if(par%rank.eq.0) then !------------------->
! longueur des enregistrements du tableau 2D:
      lrec=(iglb+2)*(jglb+2)*4
      open(unit=3,file=texte90,access='direct',form='unformatted',      & !08/10/03
       recl=lrec)
      write(3,rec=1)((glob_kmin_w(i,j),i=0,iglb+1),j=0,jglb+1)
      close(3)
      endif                  !------------------->

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out) ! assure la synchro, permet au domaine 0 de rattraper les autres
#endif
      call allocate_global('d','glob_anyvar2d       ',0,0,0)      !08-10-09
      call allocate_global('d','glob_kmin_w         ',0,0,0)      !08-10-09

!     call sigma_levels_figure
      stop 'la procedure iterative est arrivée au bout'
!..................................................................!

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE ITERATIVE
! FIN.
      endif
!_______________________________________________________________________________

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE LECTURE FICHIER
! DEBUT:
      if(ichoix.eq.2) then
!_______________________________________________________________________________

!..................................................................!
! LECTURE DES SIGMA ET KLOW APRES CONVERGENCE DU CALCUL
! Lecture par plan via tableau simple precision anyvar2D afin          !22/03/06
! de travailler avec un fichier simple precision quelque soit
! la version du model_

!#MPI
      call allocate_global('a','glob_anyvar2d       ',iglb,jglb,0)      !08-10-09
      call allocate_global('a','glob_kmin_w         ',iglb,jglb,0)      !08-10-09

      write(6,'(a,a)')'sur le point de lire:'                           &
       ,directory(1:lname1)//txtslash//'z_sig_step.out'

      open(unit=3                                                       &
       ,file=directory(1:lname1)//txtslash//'z_sig_step.out'            &
       ,access='direct'                                                 &
       ,form='unformatted'                                              &
       ,recl=(iglb+2)*(jglb+2)*4)

      do k=0,kmax+1  ! debut boucle K

      read(3,rec=k+1)((glob_anyvar2d(i,j),i=0,iglb+1),j=0,jglb+1)

       do j=0,jmax+1
       do i=0,imax+1
         depth_w(i,j,k)=glob_anyvar2d(i+par%timax(1),j+par%tjmax(1))    !23-09-09
       enddo
       enddo

      enddo      ! fin boucle K
      close(3)
      write(6,*)'ok'

      write(6,'(a,a)')'sur le point de lire:'                           &
       ,directory(1:lname1)//txtslash//'kmin_w_hyb.out' !
! longueur des enregistrements du tableau 2D:
      open(unit=3                                                       &
       ,file=directory(1:lname1)//txtslash//'kmin_w_hyb.out' & !
       ,access='direct'                                                 &
       ,form='unformatted'                                              &
       ,recl=(iglb+2)*(jglb+2)*4)
      read(3,rec=1)((glob_kmin_w(i,j),i=0,iglb+1),j=0,jglb+1)
      close(3)
      write(6,*)'ok'


       do j=0,jmax+1
       do i=0,imax+1
        kmin_w(i,j)=glob_kmin_w(i+par%timax(1),j+par%tjmax(1))          !23-09-09
       enddo
       enddo

      call allocate_global('d','glob_anyvar2d       ',0,0,0)      !08-10-09
      call allocate_global('d','glob_kmin_w         ',0,0,0)      !08-10-09

!..................................................................!

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE LECTURE FICHIER
! FIN.
      endif
!_______________________________________________________________________________

! Lignes commentées le 19-10-09
! & z readapté à la vraie bathy:                                     !19/11/05
!     do j=0,jmax+1
!     do i=0,imax+1

! lignes commentées le !29-09-09
!      do k=kmin_z(i,j)+1,kmax+1
!        depth_z(i,j,k)=-h_z(i,j)* depth_z(i,j,k)
!    &                           / depth_z(i,j,kmin_z(i,j))
!      enddo

!      depth_z(i,j,kmin_z(i,j))=-h_z(i,j)
!     enddo
!     enddo

 2222 continue                                                       !31/03/08

! Faire en sorte que la couche du fond ne soit pas trop fine: !12-01-11
      x1=0.5
      do j=0,jmax+1
      do i=0,imax+1
      k=kmin_w(i,j)
       if(        depth_w(i,j,k+1)-depth_w(i,j,k)           &
             <x1*(depth_w(i,j,k+2)-depth_w(i,j,k+1))) then

! Cas bathy ajustée:
!        depth_w(i,j,k)=depth_w(i,j,k+1)     &
!                  -x1*(depth_w(i,j,k+2)-depth_w(i,j,k+1))
!        h_w(i,j)=-depth_w(i,j,k)
! Cas profondeur du premier niveau ajustée:
         depth_w(i,j,k+1)=(depth_w(i,j,k)+x1*depth_w(i,j,k+2))/(1.+x1)

       endif
      enddo
      enddo

#ifdef bidon
! Reintroduire des niveaux dans la couche de fond:
      x1=0.25   ! pourcentage de l'épaisseur totale de la colonne d'eau
      x2=0.1    ! pourcentage de kmax pour le nbre de niveaux correspondant à x1
      do j=0,jmax+1
      do i=0,imax+1
        k0=kmin_w(i,j)+nint( x2*(kmax-kmin_w(i,j)) )
!       k0=kmin_w(i,j)+nint( x2*kmax )
        if(i==imax/2.and.j==134)write(*,*)k0,kmin_w(i,j)
        if(i==imax/2.and.j==135)write(*,*)k0,kmin_w(i,j)
        if(depth_w(i,j,k0)-depth_w(i,j,kmin_w(i,j))>x1*h_w(i,j)) then

!       write(*,*)'Avant',k0,kmin_w(i,j),depth_w(i,j,k0),(x1-1.)*h_w(i,j)
!       write(*,*)'Avant',depth_w(i,j,k0)-depth_w(i,j,kmin_w(i,j)),x1*h_w(i,j)

        rap=x1*h_w(i,j)/(depth_w(i,j,k0)-depth_w(i,j,kmin_w(i,j)))

        do k=1,kmax
         anyv1d(k,1)=depth_w(i,j,k+1)-depth_w(i,j,k)
        enddo

        sum1=0.
        do k=kmin_w(i,j),k0-1
         sum1=sum1+anyv1d(k,1)
        enddo
!       write(*,*)'Avant',sum1,rap
        sum1=0.
        do k=kmin_w(i,j),k0-1
         anyv1d(k,1)=rap*anyv1d(k,1)
         sum1=sum1+anyv1d(k,1)
         depth_w(i,j,k+1)=depth_w(i,j,k)+anyv1d(k,1)
        enddo
!       write(*,*)'Apres',sum1,depth_w(i,j,k0)

        sum2=0.
        do k=k0,kmax
         sum2=sum2+anyv1d(k,1)
        enddo
        rap=(1.-x1)*h_w(i,j)/max(sum2,small1)
        sum2=0.
        do k=k0,kmax
         anyv1d(k,1)=rap*anyv1d(k,1)
         depth_w(i,j,k+1)=depth_w(i,j,k)+anyv1d(k,1)
         sum2=sum2+anyv1d(k,1)
        enddo

!       write(*,*)'total Apres',sum1+sum2,h_w(i,j),depth_w(i,j,kmax+1)


!       stop
        endif
      enddo
      enddo
      stop

#endif

!----------------------------------------------------------------------!10/05/04
! TRAITEMENT SPECIAL DES FRONTIERES POUR QUE LES CONDITIONS DE GRADIENT
! NUL NE SOIENT PAS GENEES PAR UNE MARCHE D'ESCALIER:
! DEBUT:
!----------------------------------------------------------------------!10/05/04
! Même kmin (le + grand) à la frontière
! On n'oublie pas de modifier la bathy en consequence sinon le dz du bas
! risque d'être énorme et creer des pb pour les obc

       call obc_mixsigstep(0)

!----------------------------------------------------------------------!10/05/04
! TRAITEMENT SPECIAL DES FRONTIERES POUR QUE LES CONDITIONS DE GRADIENT
! NUL NE SOIENT PAS GENEES PAR UNE MARCHE D'ESCALIER:
! FIN.
!----------------------------------------------------------------------!10/05/04

! Couronne extraperipherique:
       call obc_h(0)
! A partir de h_z, deduire h_x, h_y et h_r
       call hz_to_hxyr

!..............................................
! On refait la condition limite au fond = niveaux "quasi collés"
!..............................................
!     small3=1.e-3
      small3=1.e-5
      do j=0,jmax+1
      do i=0,imax+1
      do k=0,kmin_w(i,j)-1
        depth_w(i,j,k)= depth_w(i,j,kmin_w(i,j))+(k-kmin_w(i,j))*small3
      enddo
      enddo
      enddo

!.........................................................
! réajuster les sigma aux z lissés:                                    !26/03/04
      do k=0,kmax+1
      do j=0,jmax+1
      do i=0,imax+1
       sigma_w(i,j,k)=( depth_w(i,j,k )- depth_w(i,j,kmin_w(i,j)))     & !23/10/06
                     /( depth_w(i,j,kmax+1)- depth_w(i,j,kmin_w(i,j)))
      enddo
      enddo
      enddo
!.........................................................

! Premier passgage:
! on ne prend pas encore en compte MORPHO_Z(I,J,NR) sinon on se
! retrouve d'office avec KLOW_X et KLOW_Y = NR-1 dans la partie
! masquée ce qui ne va pas. La mise à zero des mask à tous les
! niveaux verticaux en zone continental doit se faire apres que
! que KLOW_X et KLOW_Y aient ete calcules:
! Déduire les tableaux mask:
      do j=0,jmax+1
      do i=0,imax+1
       do k=1,kmin_w(i,j)-1
                mask_t(i,j,k)=0
       enddo
       do k=kmin_w(i,j),kmax
                mask_t(i,j,k)=1 ! premier passages
       enddo
      enddo
      enddo
! cas particulier: si il ne reste qu'un niveau on masque tout!
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_t(i,j,kmax).eq.0)  mask_t(i,j,kmax+1)=0
      enddo
      enddo

! OBC mask_c & Pour obtenir les autres points de la grille C:
      call maskt_to_maskuvp                                  !21-10-09
!.........................................................
! calculer les bornes inferieures des boucles verticales:
! pour les points de la grille C:
!.........................................................
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_u(i,j,k).eq.0)kmin_u(i,j)=min0(k+1,kmax)
       if(  mask_v(i,j,k).eq.0)kmin_v(i,j)=min0(k+1,kmax)
      enddo
      enddo
      enddo

! Deuxieme passage: mise à zero d'office des zones continentales:
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
           mask_t(i,j,k)=  mask_t(i,j,k)*  mask_t(i,j,kmax+1)
      enddo
      enddo
      enddo

! obc mask_c & obtenir les autres points de la grille C:
      call maskt_to_maskuvp                                  !21-10-09

!...........................................................
! Les niveaux sigma intermediares (SIGHL) sont déduits de
! SIGMA, ainsi que l'espacement sigma
!...........................................................
      do k=0,kmax
      do j=0,jmax+1
      do i=0,imax+1
       dsig_t(i,j,k)=sigma_w(i,j,k+1)-sigma_w(i,j,k)
      enddo
      enddo
      enddo

!...........................................................
! conditions limites verticales
!...........................................................
      small3=1e-6
      do 255 i=0,imax+1
      do 255 j=0,jmax+1
       dsig_t(i,j,kmax+1) =small3
  255 continue

! ATTENTION PROF3D_X ET PROF3D_Y jouent temporairement le rôle des
! profondeurs des niveaux ENTIER aux points de vitesse. Bientôt ils
! retrouveront leur vrai rôle, à savoir les profondeurs des  1/2 niveaux
! aux points de vitesse.
      do k=1,kmax+1
        do j=1,jmax+1
          do i=1,imax+1
             depth_u(i,j,k)=0.5*( depth_w(i,j,k)+ depth_w(i-1,j,k))
             depth_v(i,j,k)=0.5*( depth_w(i,j,k)+ depth_w(i,j-1,k))
          enddo
        enddo
      enddo

! initialisation de l'epaisseur de la colonne d'eau
      call z_thickness(0,2)

! calcul des profondeurs (en mètres) des points de grille:
      call z_levels(1)

!     call bottomboundary('init kbbc  ') !19-02-11

!     call sigma_levels_figure

      end subroutine mix_sig_step
