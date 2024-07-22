










      subroutine initial_sponge(ichoix)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 357 - last update: 21-10-22
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      use module_global
      implicit none
      integer(kind=1) ichoix ! ,flagspo_i1,flagspo_j1,flagspo_i2,flagspo_j2 !09-01-19

!......................................................................
! Version   date    Description
!         12/11/01: ajout de relaxmin_ext et de relaxmin_int
!         25/02/02: la distance sponge_l a une nouvelle signification:
!                   au dela de cette distance l'eponge est rigoureusement
!                   nulle. sponge_l/3. est maintenant l'echelle de
!                   decroissance exponentielle.
!                   disparition du seuil à l'interieur du domaine
!         27/02/02: on enleve les multiplications par le masque. Elles
!                   n'etaient pas utiles. Voir même embetantes pour
!                   certaines applications de model_streamf où on
!                   demasque certaines iles.
!         09/12/02: amenagements pratiques mineures pour calcul de la l'eponge
!         05/06/03: bienvenue à ICHOIX
!                   creation d'une zone eponge pour imbrication d'une
!                   bathy "parent" (ichoix=1)
!         09/06/03: debug cas special SPONGE_L=0.
!         25/03/04: const1=6. plus grand dans cas ichoix=0
!         24/04/04: La zone de bouclage d'update des tableaux de forcage sur
!                   zone eponge est mise en memoire dans des tableaux
!         15/06/04: Const1 calculé de sorte le poid de la relaxation soit
!                   divisé par 100 à une distance des frontieres = sponge_l
!                   Afin d'eviter une "marche" à la sortie de l'eponge on
!                   soustrait 0.01 à la fonction exp de sorte à x=sponge_l
!                   la relaxation est exactement nulle.
!         19/07/05: usage des tableaux spo_i1_x et cie à un 5eme cas qui est
!                   celui du domaine tout entier....
!         15/09/05: sponge_z forcee à zero. Attention ce signifie que l'on
!                   s'achemine progressivement vers le retrait des eponges
!                   sur les traceurs. A terme ce tableau devrait disparaitre
!         16/01/06: Retour en arriere par rapport au point precedent
!         16/02/06: Finalement le rappel vers T & S est distingué de celui
!                   sur u et v grace à relax_ts spécifié dans
!                   notebook_spongelayre
!         18/04/06: fonctions compatibles avec double precision
!         10/05/06: fonctions compatibles avec double precision
!         20/04/09: Parallelisation
!         27/04/09: Parallelisation
! 2009.2  10-06-09: correction de boucles (+jp1, +ip1)
!         02-09-09: Modif Cyril: nouvelle version de la routine gather
!         04-09-09: debug: -   mask_c au lieu de   mask_y
!                          - max0 au lieu de max
!                          - jglb+jp1 au lieu de jglb+1
!                   Le calcul des bornes de calcul est revu
!         08-09-09: Mise à jour de la fonction gather
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
!         08-10-09: compatibilité avec f95: adequation des longueurs de
!                   chaines de charactere passées en argument de allocate_global
!         13-10-09: Il faut prevoir au moins un point de forcage _z même si
!                   l'eponge est nulle
!         06-11-09: manquait l'initialisation de sponge_z et cie...
! 2010.12 05-09-10  debug: relax_ts remis à la place de relax_int
! 2010.20 22-04-11  relax_ts est séparé de sponge_t
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S26.1   05-02-13  eponges grille periodique
!         10-06-13  choix forme exp, lineaire, quadra
!         05-12-13  cas grille periodique
!         02-12-14  amenagement C.L. de h_inf_obc 
!         16-01-15  debug de l'algo de non chevauchement des boucles de la
!                   sponge layer. Ajout d'un routine de verification: 
!                   initial_sponge_doublon
!         01-02-15  ajout call obc_h(0)
!         18-05-18  eponge exponentielles
!         25-05-18  modif calcul des bornes des boucles de la sponge layer
!         02-06-18  si flag_upwind_obc=1  schema T,S upwind aux frontieres ouvertes
!         17-06-18  retour sur forme eponge avant mofif du 18-05-18
!         21-06-18  possibilite d'eponge indexee sur valeur de la resolution horizontale
!         09-01-19  possibilite d'enlever les eponges sur une frontiere
!                   meme si cette derniere est ouverte et non periodique
! v278    17-04-20  const2 pour la forme exponentielle
! v292    05-11-20  flag_stop=0 ! restaurer flag_stop avant de sortir                !05-11-20
! v299    18-03-21  utiliser obcstatus
! v301    06-04-21  STOP supprimE sur calcul sponge_l
! v309    11-10-21  forme demi cosinus au carre pour tendre plus vers zero A l'interieur !11-10-21
! v340    26-03-22  sponge_u(:,:,0) sponge_v(:,:,0) sans dimensions
! v347    28-04-22  (dist-3.) pour garantir que h_inf_=h_inf_obc sur les points 
!                   impliquEs dans les OBC
! v349    25-06-22  - Comme la dependance A la resolution comporte le risque d'avoir sponge=0 
!                   aux frontieres ouvertes on empeche le rappel d'etre  plus petit que la 
!                   valeur obtenue avec la dependance A la distance
!                   - La dependance A la resolution ne s'applique que
!                   sur T et S (u,v dependants de la distance) 
! v357    21-10-22  if(flag_sponge_txt=='full') then   !full grid-> !21-10-22
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

!     if(jperiodicboundary)flagspo_j1=0
!     if(jperiodicboundary)flagspo_j2=0
!     if(iperiodicboundary)flagspo_i1=0
!     if(iperiodicboundary)flagspo_i2=0


!*******************************************************************************
! CALCUL EPONGE POUR VARIABLES
! DEBUT:
      if(ichoix.eq.0)then                                              !05/06/03
!*******************************************************************************

       call allocate_global('a','glob_mask           ',iglb,jglb,0)    !08-10-09

       do loop1=1,3 ! loop1=1 node x, loop1=2 node y, loop1=3 node z

       if(loop1==1) then !111111>
        do j=1,jmax
        do i=1,imax+1
         glob_mask  (i+par%timax(1),j+par%tjmax(1))=  mask_u(i,j,kmax+1)
        enddo
        enddo
        ip1=1
        jp1=0
        ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax+1
        ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
       endif             !111111>

       if(loop1==2) then !222222>
        do j=1,jmax+1
        do i=1,imax
         glob_mask  (i+par%timax(1),j+par%tjmax(1))=  mask_v(i,j,kmax+1)
        enddo
        enddo
        ip1=0
        jp1=1
        ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
        ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax+1
       endif             !222222>

       if(loop1==3) then !333333>
        do j=1,jmax
        do i=1,imax
         glob_mask  (i+par%timax(1),j+par%tjmax(1))=  mask_t(i,j,kmax+1)
        enddo
        enddo
        ip1=0
        jp1=0
        ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
        ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
       endif             !333333>

      lb2=lbound(glob_mask) ; ub2=ubound(glob_mask)
      call  par_gatherall_2d(glob_mask,lb2,ub2,ind,par%nbdom)


! Afin d'avoir une division par 1000 de la relaxation à une distance = SPONGE_L
      const2=0.001 !17-04-20
      const1=-log(const2)

!---->
! La force du rappel est une fonction de la distance A la frontiere ouverte
! test if(flag_sponge_txt=='dist') commentE le 25-06-22
!     if(flag_sponge_txt=='dist') then !distance-distance->
! Meme si flag_sponge_txt=='dxdy' on calcule la depndance A la distance,
! au moins pour etablir une valeur minimum en dessous de laquelle la
! dependance A la resolution ne doit pas descendre !25-06-22
      do j=1,jmax+jp1
      do i=1,imax+ip1

        dist=1.e10

         do j1=1,jglb+jp1                                              !10-06-09

         i1=1
           if(flagspo_i1==1)                         & !09-01-19
           dist=min(dist                           &
           ,sqrt(real(i+par%timax(1)-i1)**2+       &
                 real(j+par%tjmax(1)-j1)**2)       &
           /max(glob_mask  (i1,j1)*un,small1)) !09-01-19

         i1=iglb+ip1
           if(flagspo_i2==1)                         & !09-01-19
           dist=min(dist                           &
           ,sqrt(real(i+par%timax(1)-i1)**2+       &
                 real(j+par%tjmax(1)-j1)**2)       &
           /max(glob_mask  (i1,j1)*un,small1))


         enddo ! J1

         do i1=1,iglb+ip1                                              !10-06-09

         j1=1
           if(flagspo_j1==1)                       & !09-01-19
           dist=min(dist                         &
           ,sqrt(real(i+par%timax(1)-i1)**2+     &
                 real(j+par%tjmax(1)-j1)**2)     &
           /max(glob_mask  (i1,j1)*un,small1))

         j1=jglb+jp1                                                   !04-09-09
           if(flagspo_j2==1)                       & !09-01-19
           dist=min(dist                         &
           ,sqrt(real(i+par%timax(1)-i1)**2+     &
                 real(j+par%tjmax(1)-j1)**2)     &
           /max(glob_mask  (i1,j1)*un,small1))

         enddo ! I1

       if(sponge_l > small1) then !%%%%%%%%%%%>

         anyv3d(i,j,1,loop1)=min(un,max(zero,           &
!         exp(-const1*dist/sponge_l)-const2             & ! forme exponentielle
!         ((sponge_l-dist)/sponge_l)**2                 & ! forme quadratique
!          (sponge_l-dist)/sponge_l                     & ! forme lineaire
!         0.5*(1.-cos(pi*(sponge_l-dist)/sponge_l))     & ! forme demi cosinus
         (0.5*(1.-cos(pi*(sponge_l-dist)/sponge_l)))**2 & ! forme demi cosinus au carre pour tendre plus vers zero A l'interieur !11-10-21
                                )      )

         if(dist > sponge_l)anyv3d(i,j,1,loop1)=0.

       else                       !%%%%%%%%%%%>

         anyv3d(i,j,1,loop1)=0.

       endif                      !%%%%%%%%%%%>

      enddo ! I
      enddo ! J

!     endif                            !distance-distance->
!---->
! La force du rappel est une fonction de resolution horizontale (plus basse Aux frontieres ouvertes)
! Note: comme la dependance A la resolution comporte le risque d'avoir sponge=0 aux frontieres ouvertes 
! (si dx est localement petit) alors on garantit un minimum de rappel en empechant le rappel d'etre 
! plus petit que la valeur obtenue avec la dependance A la distance, soit max(....,anyv3d(i,j,1,loop1))
      if(flag_sponge_txt=='dxdy') then   !resolution-resolution-> !21-06-18
! On applique le rappel depedant de la resolution que sur T et S seulement !25-06-22
        if(loop1==3) then !-t->
          do j=1,jmax ; do i=1,imax
            anyv3d(i,j,1,loop1)=max(0.5*(1.+tanh( (max(dx_t(i,j),dy_t(i,j)) - sponge_dx_critic )/sponge_dx_width ) ) -0.001 ,anyv3d(i,j,1,loop1)) !25-06-22
          enddo       ; enddo
        endif             !-t->
      endif                              !resolution-resolution->
      if(flag_sponge_txt=='full') then   !full grid-> !21-10-22
! On applique le rappel 'full grid' sur T et S seulement 
        if(loop1==3) then !-t->
          do j=1,jmax ; do i=1,imax
            anyv3d(i,j,1,loop1)=1.
          enddo       ; enddo
        endif             !-t->
      endif                              !full grid-> !21-10-22


      enddo ! LOOP1
      call allocate_global('d','glob_mask           ',0,0,0)           !08-10-09


! Initialiser les tableaux sponge à partir de anyv3d:                  !06-11-09
      do j=1,jmax
      do i=1,imax+1
       sponge_u(i,j,0)=anyv3d(i,j,1,1)*flag_upwind_obc !26-03-22
       sponge_u(i,j,1)=anyv3d(i,j,1,1)*relax_int
       sponge_u(i,j,2)=anyv3d(i,j,1,1)*relax_ext
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       sponge_v(i,j,0)=anyv3d(i,j,1,2)*flag_upwind_obc !26-03-22
       sponge_v(i,j,1)=anyv3d(i,j,1,2)*relax_int
       sponge_v(i,j,2)=anyv3d(i,j,1,2)*relax_ext
      enddo
      enddo
      do j=1,jmax
      do i=1,imax
       sponge_t(i,j,1)=anyv3d(i,j,1,3)                                !22-04-11
      enddo
      enddo
      if(flag_upwind_obc==1) then !m°v°m> !02-06-18
       if(obcstatus(ieq1)==1)upwindriver_t(0     ,0:jmax+1)=0.
       if(obcstatus(ieqimax)==1)upwindriver_t(imax+1,0:jmax+1)=0.
       if(obcstatus(jeq1)==1)upwindriver_t(0:imax+1,0     )=0.
       if(obcstatus(jeqjmax)==1)upwindriver_t(0:imax+1,jmax+1)=0.
       do j=1,jmax ; do i=1,imax
        upwindriver_t(i,j)=max(min(1.-sponge_t(i,j,1),upwindriver_t(i,j)),0.) 
       enddo       ; enddo
       call obc_scal_mpi_upwindriver('z1') !02-06-18
      endif                       !m°v°m>


! AFIN DE REDUIRE LA TAILLE DES BOUCLES D'AVANCEE DANS LE TEMPS DES TABLEAUX
! DE FORCAGE LES BOUCLES SE REDUIRONT A LA ZONE EPONGE. IL FAUT POUR CELA
! DEFINIR LES I ET J CONCERNES. ON DETERMINE ICI LES INTERVALLES DE CALCUL
! SUR I ET J POUR LES POINTS _T _U ET _V
! Modification du codage le !25-05-18
!......................................................................
! POINTS T
      i1=     1+nint(sponge_l)
      i2=iglb  -nint(sponge_l)
      j1=     1+nint(sponge_l)
      j2=jglb  -nint(sponge_l)

      i1=min(max(i1,1),iglb  )
      i2=min(max(i2,1),iglb  )
      j1=min(max(j1,1),jglb  )
      j2=min(max(j2,1),jglb  )

! "Evite doublon": il ne faut pas faire deux fois les
! operations de sommation dans update_tide. Donc pas de chevauchement.
      i2=max(i2,i1+1) !  Evite doublon !16-01-15
      j2=max(j2,j1+1) !  Evite doublon !16-01-15

! Frontiere i=1
      spo_i1_t(1)=1 
      spo_i2_t(1)=i1 
      spo_j1_t(1)=1 
      spo_j2_t(1)=jglb 

! Frontiere i=iglb
      spo_i1_t(2)=i2 
      spo_i2_t(2)=iglb 
      spo_j1_t(2)=1 
      spo_j2_t(2)=jglb 

! Frontiere j=1
      spo_i1_t(3)=i1+1  !evite doublon
      spo_i2_t(3)=i2-1  !evite doublon
      spo_j1_t(3)=1 
      spo_j2_t(3)=j1 

! Frontiere j=jglb
      spo_i1_t(4)=i1+1 !evite doublon
      spo_i2_t(4)=i2-1 !evite doublon
      spo_j1_t(4)=j2
      spo_j2_t(4)=jglb

! Les 4 frontieres d'un coup:                                          !19/07/05
      spo_i1_t(5)=1
      spo_i2_t(5)=iglb
      spo_j1_t(5)=1
      spo_j2_t(5)=jglb

! Passer de global en local
      do loop1=1,5
       spo_i1_t(loop1)=spo_i1_t(loop1)-par%timax(1)
       spo_i2_t(loop1)=spo_i2_t(loop1)-par%timax(1)
       spo_j1_t(loop1)=spo_j1_t(loop1)-par%tjmax(1)
       spo_j2_t(loop1)=spo_j2_t(loop1)-par%tjmax(1)
      enddo

! 1- Empecher les boucles de se derouler dans les domaines interieurs hors
! de la sponge layer. Principe: rendre la borne inferieure plus grande
! que la borne superieure.
! 2- Limiter les boucles aux bornes locales
      do loop1=1,5     !--loop1-->

       if(spo_i1_t(loop1)>imax) then
          spo_i1_t(loop1)=imax+1
          spo_i2_t(loop1)=imax
       else
          spo_i1_t(loop1)=max(spo_i1_t(loop1),1)
       endif

       if(spo_i2_t(loop1)<1) then
          spo_i2_t(loop1)=0
          spo_i1_t(loop1)=1
       else
          spo_i2_t(loop1)=min(spo_i2_t(loop1),imax)
       endif

       if(spo_j1_t(loop1)>jmax) then
          spo_j1_t(loop1)=jmax+1
          spo_j2_t(loop1)=jmax
       else
          spo_j1_t(loop1)=max(spo_j1_t(loop1),1)
       endif

       if(spo_j2_t(loop1)<1) then
          spo_j2_t(loop1)=0
          spo_j1_t(loop1)=1
       else
          spo_j2_t(loop1)=min(spo_j2_t(loop1),jmax)
       endif

      enddo            !--loop1-->

!......................................................................
! POINTS U
      i1=     1+nint(sponge_l)
      i2=iglb+1-nint(sponge_l)
      j1=     1+nint(sponge_l)
      j2=jglb  -nint(sponge_l)

      i1=min(max(i1,1),iglb+1)
      i2=min(max(i2,1),iglb+1)
      j1=min(max(j1,1),jglb  )
      j2=min(max(j2,1),jglb  )

! "Evite doublon": il ne faut pas faire deux fois les
! operations de sommation dans update_tide. Donc pas de chevauchement.
      i2=max(i2,i1+1) !  Evite doublon !16-01-15
      j2=max(j2,j1+1) !  Evite doublon !16-01-15

! Frontiere i=1
      spo_i1_u(1)=1 
      spo_i2_u(1)=i1 
      spo_j1_u(1)=1 
      spo_j2_u(1)=jglb 

! Frontiere i=iglb+1
      spo_i1_u(2)=i2 
      spo_i2_u(2)=iglb+1 
      spo_j1_u(2)=1 
      spo_j2_u(2)=jglb 

! Frontiere j=1
      spo_i1_u(3)=i1+1  !evite doublon
      spo_i2_u(3)=i2-1  !evite doublon
      spo_j1_u(3)=1 
      spo_j2_u(3)=j1 

! Frontiere j=jglb
      spo_i1_u(4)=i1+1 !evite doublon
      spo_i2_u(4)=i2-1 !evite doublon
      spo_j1_u(4)=j2
      spo_j2_u(4)=jglb

! Les 4 frontieres d'un coup:                                          !19/07/05
      spo_i1_u(5)=1
      spo_i2_u(5)=iglb+1
      spo_j1_u(5)=1
      spo_j2_u(5)=jglb

! Passer de global en local
      do loop1=1,5
       spo_i1_u(loop1)=spo_i1_u(loop1)-par%timax(1)
       spo_i2_u(loop1)=spo_i2_u(loop1)-par%timax(1)
       spo_j1_u(loop1)=spo_j1_u(loop1)-par%tjmax(1)
       spo_j2_u(loop1)=spo_j2_u(loop1)-par%tjmax(1)

      enddo

! 1- Empecher les boucles de se derouler dans les domaines interieurs hors
! de la sponge layer. Principe: rendre la borne inferieure plus grande
! que la borne superieure.
! 2- Limiter les boucles aux bornes locales
      do loop1=1,5     !--loop1-->

       if(spo_i1_u(loop1)>imax+1) then
          spo_i1_u(loop1)=imax+1+1
          spo_i2_u(loop1)=imax+1
       else
          spo_i1_u(loop1)=max(spo_i1_u(loop1),1)
       endif

       if(spo_i2_u(loop1)<1) then
          spo_i2_u(loop1)=0
          spo_i1_u(loop1)=1
       else
          spo_i2_u(loop1)=min(spo_i2_u(loop1),imax+1)
       endif

       if(spo_j1_u(loop1)>jmax) then
          spo_j1_u(loop1)=jmax+1
          spo_j2_u(loop1)=jmax
       else
          spo_j1_u(loop1)=max(spo_j1_u(loop1),1)
       endif

       if(spo_j2_u(loop1)<1) then
          spo_j2_u(loop1)=0
          spo_j1_u(loop1)=1
       else
          spo_j2_u(loop1)=min(spo_j2_u(loop1),jmax)
       endif

      enddo            !--loop1-->

!......................................................................
! POINTS V
      i1=     1+nint(sponge_l)
      i2=iglb  -nint(sponge_l)
      j1=     1+nint(sponge_l)
      j2=jglb+1-nint(sponge_l)

      i1=min(max(i1,1),iglb  )
      i2=min(max(i2,1),iglb  )
      j1=min(max(j1,1),jglb+1)
      j2=min(max(j2,1),jglb+1)

! "Evite doublon": il ne faut pas faire deux fois les
! operations de sommation dans update_tide. Donc pas de chevauchement.
      i2=max(i2,i1+1) !  Evite doublon !16-01-15
      j2=max(j2,j1+1) !  Evite doublon !16-01-15

! Frontiere i=1
      spo_i1_v(1)=1 
      spo_i2_v(1)=i1 
      spo_j1_v(1)=1 
      spo_j2_v(1)=jglb+1 

! Frontiere i=iglb
      spo_i1_v(2)=i2 
      spo_i2_v(2)=iglb 
      spo_j1_v(2)=1 
      spo_j2_v(2)=jglb+1 

! Frontiere j=1
      spo_i1_v(3)=i1+1  !evite doublon
      spo_i2_v(3)=i2-1  !evite doublon
      spo_j1_v(3)=1 
      spo_j2_v(3)=j1 

! Frontiere j=jglb+1
      spo_i1_v(4)=i1+1 !evite doublon
      spo_i2_v(4)=i2-1 !evite doublon
      spo_j1_v(4)=j2
      spo_j2_v(4)=jglb+1

! Les 4 frontieres d'un coup:                                          !19/07/05
      spo_i1_v(5)=1
      spo_i2_v(5)=iglb
      spo_j1_v(5)=1
      spo_j2_v(5)=jglb+1

! Passer de global en local
      do loop1=1,5
       spo_i1_v(loop1)=spo_i1_v(loop1)-par%timax(1)
       spo_i2_v(loop1)=spo_i2_v(loop1)-par%timax(1)
       spo_j1_v(loop1)=spo_j1_v(loop1)-par%tjmax(1)
       spo_j2_v(loop1)=spo_j2_v(loop1)-par%tjmax(1)
      enddo

! 1- Empecher les boucles de se derouler dans les domaines interieurs hors
! de la sponge layer. Principe: rendre la borne inferieure plus grande
! que la borne superieure.
! 2- Limiter les boucles aux bornes locales
      do loop1=1,5     !--loop1-->

       if(spo_i1_v(loop1)>imax) then
          spo_i1_v(loop1)=imax+1
          spo_i2_v(loop1)=imax
       else
          spo_i1_v(loop1)=max(spo_i1_v(loop1),1)
       endif

       if(spo_i2_v(loop1)<1) then
          spo_i2_v(loop1)=0
          spo_i1_v(loop1)=1
       else
          spo_i2_v(loop1)=min(spo_i2_v(loop1),imax)
       endif

       if(spo_j1_v(loop1)>jmax+1) then
          spo_j1_v(loop1)=jmax+1+1
          spo_j2_v(loop1)=jmax+1
       else
          spo_j1_v(loop1)=max(spo_j1_v(loop1),1)
       endif

       if(spo_j2_v(loop1)<1) then
          spo_j2_v(loop1)=0
          spo_j1_v(loop1)=1
       else
          spo_j2_v(loop1)=min(spo_j2_v(loop1),jmax+1)
       endif

      enddo            !--loop1-->

! Verifier qu'il n'y a pas de chevauchement dans les bornes des
! boucles de la sponge layer:
      call initial_sponge_doublon !16-01-15

!*******************************************************************************
! CALCUL EPONGE POUR VARIABLES
! FIN.
      flag_stop=0 ! restaurer flag_stop avant de sortir                !05-11-20
      return                                                           !05/06/03
      endif                                                            !05/06/03
!*******************************************************************************



!*******************************************************************************
! CALCUL EPONGE BATHY "PARENT" IMBRIQUEE
! DEBUT:
      if(ichoix.eq.1)then                                              !05/06/03
!*******************************************************************************
      call allocate_global('a','glob_mask           ',iglb,jglb,0)    !08-10-09

! eponge large 0 < CONST1 < 1, eponge etroite CONST1 > 1
      const1=1.

      do 40 i=1,imax
      do 40 j=1,jmax
         dist=1.e10
         do 41 j1=1,jglb

         i1=1
           dist=min(dist                                                &
           ,sqrt(real(i+par%timax(1)-i1)**2+real(j+par%tjmax(1)-j1)**2) &
           /max(un*glob_mask  (i1,j1),small1))

         i1=iglb
           dist=min(dist                                                &
           ,sqrt(real(i+par%timax(1)-i1)**2+real(j+par%tjmax(1)-j1)**2) &
           /max(un*glob_mask  (i1,j1),small1))

   41 continue
         do 42 i1=1,iglb

         j1=1
           dist=min(dist                                                &
           ,sqrt(real(i+par%timax(1)-i1)**2+real(j+par%tjmax(1)-j1)**2) &
           /max(un*glob_mask  (i1,j1),small1))

         j1=jglb
           dist=min(dist                                                &
           ,sqrt(real(i+par%timax(1)-i1)**2+real(j+par%tjmax(1)-j1)**2) &
           /max(un*glob_mask  (i1,j1),small1))

   42 continue

      if(sponge_l < small1) then !%%%%%%%%%%%>                        !09/06/03

      sponge_t(i,j,1)=exp(-dist/sponge_l)

      else                        !%%%%%%%%%%%>                        !09/06/03
      sponge_t(i,j,1)=0.
      endif                       !%%%%%%%%%%%>                        !09/06/03


   40 continue

      call allocate_global('d','glob_mask           ',0,0,0)    !08-10-09
!*******************************************************************************
! CALCUL EPONGE BATHY "PARENT" IMBRIQUEE
! FIN.
      return                                                           !05/06/03
      endif                                                            !05/06/03
!*******************************************************************************

      end subroutine initial_sponge

!..............................................................................

      subroutine initial_sponge_hmin
      use module_principal
      use module_parallele !#MPI
      use module_global
      implicit none
      double precision h_inf_
!     integer(kind=1) flagspo_i1,flagspo_j1,flagspo_i2,flagspo_j2 !09-01-19

!     flagspo_i1=1 ; flagspo_j1=1 ; flagspo_i2=1 ; flagspo_j2=1   !05-02-13 !09-01-19
!     if(jperiodicboundary)flagspo_j1=0
!     if(jperiodicboundary)flagspo_j2=0
!     if(iperiodicboundary)flagspo_i1=0
!     if(iperiodicboundary)flagspo_i2=0

      do 40 i=1,imax
      do 40 j=1,jmax


        dist=1.e10

!       if(flag_i_==1) then !>>>>>>>>>>>>>
         do j1=1,jglb                                             !10-06-09

         i1=1
           if(flagspo_i1==1)                                 & !09-01-19
           dist=min(dist                                   &
           ,sqrt(        real(i+par%timax(1)-i1)**2+       &
                         real(j+par%tjmax(1)-j1)**2        &
                 +1.)  & !Ajout d'une cst pour eviter dist=0 en i=1 !03-12-14
           /max(glob_mask  (i1,j1)*un,small1))

         i1=iglb
           if(flagspo_i2==1)                                 & !09-01-19
           dist=min(dist                                   &
           ,sqrt(real(i+par%timax(1)-i1)**2+       &
                 real(j+par%tjmax(1)-j1)**2        &
                 +1.)                                      &
           /max(glob_mask  (i1,j1)*un,small1))


         enddo ! J1
!       endif               !>>>>>>>>>>>>>

!       if(flag_j_==1) then !>>>>>>>>>>>>>
         do i1=1,iglb                                              !10-06-09

         j1=1
           if(flagspo_j1==1)                               & !09-01-19
           dist=min(dist                                 &
           ,sqrt(        real(i+par%timax(1)-i1)**2+     &
                         real(j+par%tjmax(1)-j1)**2      &
                 +1.)                                    &
           /max(glob_mask  (i1,j1)*un,small1))

         j1=jglb                                                   !04-09-09
           if(flagspo_j2==1)                               & !09-01-19
           dist=min(dist                                 &
           ,sqrt(        real(i+par%timax(1)-i1)**2+     &
                         real(j+par%tjmax(1)-j1)**2      &
                 +1.)                                    &
           /max(glob_mask  (i1,j1)*un,small1))

         enddo ! I1
!       endif               !>>>>>>>>>>>>>

!     x1=min(dist/30.,un)                  ! transition sur 30 points
      x1=max(min((dist-3.)/30.,1.),0.)     !28-04-22 (dist-3.) pour garantir que h_inf_=h_inf_obc sur les points impliquEs dans les OBC
      h_inf_=h_inf*x1+h_inf_obc*(1.-x1)

      h_w(i,j)=min(max(h_inf_,h_w(i,j)),h_sup)

   40 continue

      call obc_h(0) ! 01-02-15

      end subroutine initial_sponge_hmin
!.............................................................
      subroutine initial_sponge_doublon
      use module_principal
      use module_parallele

! Verifie que les boucles spo_ ne font pas de doublon.
! C'est important par exemple pour s'assurer que la 
! sommation de la solution de maree (update_tide) ne
! se fait pas en double (ou plus)
      xy_u(:,:,1)=0.
      xy_v(:,:,1)=0.
      xy_t(:,:,1)=0.

       do loop1=1,4 

        do j=spo_j1_u(loop1),spo_j2_u(loop1)
        do i=spo_i1_u(loop1),spo_i2_u(loop1)
            xy_u(i,j,1)=xy_u(i,j,1)+1.
         if(xy_u(i,j,1)>1.1) then
           write(10+par%rank,*)'doublon sponlayer u'
           write(10+par%rank,*)'par%rank,i,j loc',par%rank,i,j
           write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
           do loop2=1,4
            write(10+par%rank,*)'loop2',loop2
            write(10+par%rank,*)'spo_i1_u(loop2),spo_i2_u(loop2)',spo_i1_u(loop2),spo_i2_u(loop2)
            write(10+par%rank,*)'spo_j1_u(loop2),spo_j2_u(loop2)',spo_j1_u(loop2),spo_j2_u(loop2)
           enddo 
           stop 'Stop doublon u dans initial_sponge voir fort files'
         endif
        enddo
        enddo

        do j=spo_j1_v(loop1),spo_j2_v(loop1)
        do i=spo_i1_v(loop1),spo_i2_v(loop1)
            xy_v(i,j,1)=xy_v(i,j,1)+1.
         if(xy_v(i,j,1)>1.1) then
           write(10+par%rank,*)'doublon sponlayer v'
           write(10+par%rank,*)'par%rank,i,j loc',par%rank,i,j
           write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
           do loop2=1,4
            write(10+par%rank,*)'loop2',loop2
            write(10+par%rank,*)'spo_i1_v(loop2),spo_i2_v(loop2)',spo_i1_v(loop2),spo_i2_v(loop2)
            write(10+par%rank,*)'spo_j1_v(loop2),spo_j2_v(loop2)',spo_j1_v(loop2),spo_j2_v(loop2)
           enddo 
          stop 'Stop doublon v dans initial_sponge voir fort files'
         endif
        enddo
        enddo

        do j=spo_j1_t(loop1),spo_j2_t(loop1)
        do i=spo_i1_t(loop1),spo_i2_t(loop1)
            xy_t(i,j,1)=xy_t(i,j,1)+1.
         if(xy_t(i,j,1)>1.1) then
           write(10+par%rank,*)'doublon sponlayer t'
           write(10+par%rank,*)'par%rank,i,j loc',par%rank,i,j
           write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
           do loop2=1,4
            write(10+par%rank,*)'loop2',loop2
            write(10+par%rank,*)'spo_i1_t(loop2),spo_i2_t(loop2)',spo_i1_t(loop2),spo_i2_t(loop2)
            write(10+par%rank,*)'spo_j1_t(loop2),spo_j2_t(loop2)',spo_j1_t(loop2),spo_j2_t(loop2)
           enddo 
          stop 'Stop doublon t dans initial_sponge voir fort files'
         endif
        enddo
        enddo

       enddo !loop1

      end subroutine initial_sponge_doublon
