










      subroutine sigma_levels
!______________________________________________________________________
! SYMPHONIE ocean model
! release 361 - last update: 28-12-22
!______________________________________________________________________

      use module_principal ; use module_parallele ; use module_grid
      implicit none

!...............................................................................
! Version date     Description des modifications:
!  modif  15/07/01: passage à la coordonnée sigma généralisée.
!         17/04/02: modif sur melange en contrainte et sigma simple
!         24/03/03: changement de notebook:
!                   notebook_technum remplacé par notebook_vertcoord
!                   de + ce dernier est maintenant lu dans set_parameters.f
!         16/03/04: application d'un limiteur de resolution verticale
!         24/03/04: une autre distribution par defaut
!         26/03/04: amelioration du limiteur de resolution verticale
!         26/01/06: anyv3D remplace anyvar3D
!         27/12/07: par defaut, distribution verticale homogene
!         26/03/08: on adopte un autre calcul de dsig_x et dsig_y pour qu'en
!                   coordonnee generalisee, on ait dz_x = 0.5 (dz_z(i)+dz_z(i-1))
!                                                  dz_y = 0.5 (dz_z(i)+dz_z(j-1))
!         31/03/08: possibilite d'une grille hybride sans procedure iterative
! 2009.2  10-06-09: Remplacer h(ihmin,jhmin) par hmin
!                   Remplacer h(ihmax,jhmax) par hmax
!         09-09-09: Essai d'une grille verticale à increment de densité
!                   constant (à noter que le notebook_vertcoord n'est pas encore
!                   adapté à cet essai de grille)
! 2010.8  03-05-10  Supression du point precedent
! 2010.9  31-05-10  Ajout d'une grille partial step
! 2010.11 19-07-10  Seul le proc 0 recrit dans messages
!         26-07-10  Ajout d'une grille combinant partial step (ocean profond)
!                   sigma (shelf)
!         26-09-10  suite du point precedent
!         07-10-10  sous couvert de commenter/decommenter des lignes, possibilité
!                   de retrouver un partial step dans les tres petites profondeurs
!                   afin d'empecher la concentration des niveaux
! 2010.13 02-11-10  - affinage du cas 3
!                   - des figures 2DV de la grille
! 2010.14 24-11-10  modification sur distribution des niveaux de l'option sigma
!                   generalisee: l'exposant x1 remplace 0.5 pour eviter que le
!                   premier niveau au dessus du fond ne soit trop important et
!                   l'exposant 1 remplace 0.5 pour le cas de la contrainte sur la
!                   couche limite de fond pour au contraire eviter qu'il ne soit
!                   trop petit
! 2010.16 28-12-10  coupe 2DV graphique de la grille en cas de parallelisation
! S.26    19-03-14  amenagements pour nemo offline
!         08-04-14  ne plus ecrire les fichier z2dv...
!         30-04-14  ajout cas nemo offline
!         02-05-14  En cas de procedure offline sigma_w est lu dans le fichier de
!                   offline grille.nc
!         27-06-14  legere mofif de l'aiguillage ioffline=2
!         25-11-14  ajout subroutine sigma_levels_f
!         31-01-15  - modif de la fonction de transition sigma/sigma_generalisee
!                   dans le but de garantir dz(k+1)<dz(k)
!                   - faire des fichiers ascii profils de z,dz par proc au point hmax
!                   dans repertoire tmp
!         07-02-15  xy_t(:,:,1) est une elevation de surface non nulle en cas de zone
!                   de recouvrement...
!         11-03-15  sigma_w est alloue dans reset.F90
!         16-05-15  un peu de clarification....
!         25-05-15  mieux garantir dz/dk>0
!         05-06-15  suite du point precedent
!         05-08-15  la subroutine sigma_levels_hydrocons permet de satisfaire le
!                   critere de coherence hydrostatique de Haney 1991
!         24-01-17  un lien vers un doc explicatif
!         16-11-17  ajout subroutine sigma_levels_merged
!         06-12-17  modifs sigma_levels_merged
!         07-12-17  modifs sigma_levels_merged
!         25-01-18  ajout d'un seuil preservant la resolution verticale pres de la surface
!         27-01-18  suite point 25-01-18
!         14-02-18  nbvstepmin nombre minimum de niveauw sur la verticale
!                   hstepmax=min(hstepmax,hmax)
!         26-03-18  - botlevmerged_w supprimE
!                   - ajout subroutine sigma_levels_merged_finalize
!         14-04-18  ne pas passer par sigma_levels_dzsurfmin si kmax=1
!         14-06-18  lignes commentEes
!         28-10-18  if(flag_z2dv_outputs==1)call sigma_levels_2dv_plot
!         13-12-18  ajout d'un critere de proportionalite sur l'epaisseur des couches fusionnees
! v256    06-06-19  retablissement de dsig_t=max(dsig_t,small1) 
! v269    05-12-19  kmergedr4_ 
! v271    14-12-19  ameliorations de la construction de la grille hybride
! v274    09-02-20  Dans les zones avec le "miminum" de niveaux desormais kmin=kmerged
! v278    15-04-20  modif formule standard de la bathy enveloppe
!                   if(flag_merged_levels==1)call sigma_levels_envelope_bathy 
!                   subroutine sigma_levels_k_2_zw 
!                   rankhmax
! v279    25-04-20  ne plus utiliser hstepmax
! v281    07-05-20  dz_vertical_incr_fact
! v284    21-05-20  upwzone0_t prend une fonction differente pour desormais definir des zones upwind en 3D
! v285    04-06-20  subroutine sigma_levels_kmergedr4 basee sur dz plutot que dsig pour 
!                   pouvoir supprimer les tableaux dsig
!         06-06-20  dsig_t deduit de sigma_w
! v292    07-11-20  Securite sur sigma_w !07-11-20
!         22-11-20  ajout depth_w(kmax+1) au dessus de -h sigma simple   
! v297    05-03-21  ajout x4 en facteur du terme de courbure de la bathy enveloppe.
! v298    07-03-21  indices globaux pour les fichiers z2dv_OE et _SN
! v300    22-03-21  parametres VQS introduits dans notebook_vertcoord
! v301    29-04-21  parametres VQS introduits dans notebook_vertcoord
! v310    25-10-21  coordonnee VQS: la couche kmin peut etre plus grande que la couche 
!                   kmerged (juste au dessus de kmin)
!         02-11-21  - ajout dsigmerged_u et dsigmerged_v pour reducteur de flux dans 
!                   couple_modes_scalars
!                   - un seuil minimum sur delta sigma pour eviter dz trop petit en zone assechee
!         03-11-21  - ajout pgfratio_u,pgfratio_v !03-11-21
! v312    17-11-21  une autre formule pour dsigmerged_u
! v313    26-11-21  suite point precedent: if(kmerged_ (i,j)>kmin_ (i,j)) then !pmx>
! v318    23-12-21  - retablir call lissebathy(1,0.,1,iglb,1,jglb,10) !23-12-21
!                   - Il n'y a pas de raison que la bathy envelope soit plus profonde que hmax
! v319    27-12-21  if(kmin_u(i,j)/=kmerged_u(i,j)) then !m°v°m> !27-12-21
! v361    28-12-22  ajout mask_vqs_tke_w
!...............................................................................
!    _________                    .__                  .__                     !m[°v°]m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

! Details dans: !24-01-17
! https://docs.google.com/document/d/1HyHCalPTAEThRfmvib_4fWA7dMacELrUdSKx2tz2OIw/edit

! Ce sous programme calcule la distribution verticale des
! niveaux verticaux sigma, et des masques sur la direction
! verticale.
! Les niveaux sigma "entiers" s'appellent SiGMA_Z. Sur ces niveaux
! on trouve les parametres de la turbulence et la vitesse verticale.
! Les niveaux sigma "half level" ou intermediaire s'appellent SiGHL_Z.
! On y trouve les vitesses horizontales et les traceurs.
! L'espacement entre les niveaux SiGMA s'appelle DSiG_Z DiSiG_X DSiG_Y.
! La version "simple sigma" ne requiere pas une dimension 3 pour les
! tableaux (elle est présente neanmoins pour generaliser la coordonnee
! si on le souhaite) ce qui explique que l'on determine les parametres
! pour une colonne (i=j=1) en premier lieu et qu'enfin on repercute sur
! toutes les colonnes.

!................................................................
! Particular case where the grid is loaded from an external file:
      if(ioffline==2) then !grdgrdgrdgrdgrdgrdgrdgrd> !27-06-14

! Case where the grid is loaded from a NEMO grid file
       if(initialgridfile_txt/='none') then ! nemo case nemo case >
        call initial_depthw_from_file
        return
       endif                                ! nemo case nemo case >

! Case where the grid is loaded from a S grid file
! Cas simulation offline (bio) le tableau sigma_w est lu dans le fichier grille
! du repertoire OffLiNE et dsig_t est deduit de sigma_w
       call grid_sigma_offlinecase
       return

      endif                !grdgrdgrdgrdgrdgrdgrdgrd>


      if(igesig>=2)stop 'igesig>=2 obsolete'

      call sigma_levels_depth_w        ! depth_w first guess et couches merged
      call sigma_levels_consistency    ! dsig_t sigma_w
      call z_thickness(1,1)            ! computes hz_w, hz_u, hz_v
      call cellbox_thickness_sigma(1)  ! computes dz_t dz_u dz_v
      call z_levels(1)                 ! computes depth_w depth_t depth_u depth_v
      call sigma_levels_kmergedr4      ! computes kmergedr4_u kmergedr4_v 04-06-20
      if(flag_z2dv_outputs==1)call sigma_levels_2dv_plot !28-10-18
      if(ihybsig==1) stop 'stop sigma_levels sigstepgrid file created!'

!     deallocate(sigma_w)

! Faire des fichiers profil de distribution de z, par proc, au point hmax !31-01-15
! lignes en suivant commentEes le 14-06-18
 !    write(texte30,'(a,i0)')trim(tmpdirname)//'zlevelprofile_',par%rank
 !    open(unit=3,file=texte30)
 !    write(3,*)'k depth_w'
 !    do k=kmax+1,1,-1
 !     write(3,*)k,depth_w(ihmax,jhmax,k)
 !    enddo
 !    write(3,*)'----------------------'
 !    write(3,*)'k dz dz(k)/dz(k+1)'
 !    do k=kmax  ,1,-1
 !     if(k==kmax) then
 !      write(3,*)k,depth_w(ihmax,jhmax,k+1)-depth_w(ihmax,jhmax,k)
 !     else
 !      write(3,*)k,depth_w(ihmax,jhmax,k+1)-depth_w(ihmax,jhmax,k)  &
 !                ,(depth_w(ihmax,jhmax,k+1)-depth_w(ihmax,jhmax,k)) &
 !                /(depth_w(ihmax,jhmax,k+2)-depth_w(ihmax,jhmax,k+1))
 !     endif
 !    enddo
 !    close(3)

      end subroutine sigma_levels

!.............................................................

      subroutine sigma_levels_depth_w
      use module_principal ; use module_parallele ; use module_grid
      implicit none

!...........................................................
! ETAPE 1:
! Un profil type avant normalisation ( Sum(dsig)=1 )
!...........................................................

!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!     stop 'kokokoko3'

! Cas First guess homogene:
      if(isigfile==0) then !>>>>
! Reset par defaut: dsig=1
       dsig_t(1,1,:)=1.
       if(igesig==0) then !----simple-----> !07-05-20
        dsig_t(1,1,kmax)=1.
        do k=kmax-1,1,-1
         dsig_t(1,1,k)=dsig_t(1,1,k+1)*dz_vertical_incr_fact !07-05-20
!        dsig_t(1,1,k)=dsig_t(1,1,k+1) &
!          *(1.+(dz_vertical_incr_fact-1.)*real(kmax-k)/real(kmax-1))
        enddo
       endif              !----simple----->
      endif                !>>>>


! Cas First guess lu dans fichier dsigma.in
      if(isigfile==1) then !ooo>
       texte80(1)=directory(1:lname1)//txtslash//'dsigma.in'
       open(unit=3,file=texte80(1))
        do k=1,kmax
          read(3,*)dsig_t(1,1,k)
        enddo
       close(3)
      endif                !ooo>


!...........................................................
! ETAPE 2:
! 2eme passage: sigma 1ere estimation ! SiGMA est deduite de DSiG
!...........................................................
      i=1 ; j=1 ; sigma_w(i,j,1)=0.
      do k=1,kmax
        sigma_w(i,j,k+1)=sigma_w(i,j,k)+dsig_t(i,j,k)
      enddo

!...........................................................
! ETAPE 3:
! 3eme passage: normalisation sigma: SiGMA(1)=0. SiGMA(kmax)=1.
!...........................................................
      i=1 ; j=1
      do k=1,kmax+1
        sigma_w(i,j,k)=sigma_w(i,j,k)/sigma_w(i,j,kmax+1)
      enddo

! Appel placE apres calcul de sigma(1,1,k) utilise ensuite
      if(flag_merged_levels==1)call sigma_levels_envelope_bathy(1) !15-04-20

      ksecu=0

!-----------------------------------------------------
!  simple (horizontal-independent) sigma distribution:
      if(igesig==0) then !----simple----->

      ksecu=1
!...........................................................
! ETAPE 1:
! Les niveaux sigma intermediares (SiGHL) sont déduits de
! SiGMA, ainsi que l'espacement sigma maintenant normalisé.
!...........................................................
      i=1 ; j=1
      do k=1,kmax
       dsig_t(i,j,k)=sigma_w(i,j,k+1)-sigma_w(i,j,k)
      enddo

!...........................................................
! ETAPE 2: ! conditions limites inferieures et superieures
!...........................................................
      i=1 ; j=1
      sigma_w(i,j,0     )=0.
       dsig_t(i,j,kmax+1)=dsig_t(i,j,kmax)
       dsig_t(i,j,0     )=dsig_t(i,j,1   )        ! 15/07/01

!...........................................................
! ETAPE 3: grille 3D et depth_w
!...........................................................
      do k=0,kmax+1 ; do i=0,imax+1 ; do j=0,jmax+1
        sigma_w(i,j,k)=sigma_w(1,1,k)
         dsig_t(i,j,k)= dsig_t(1,1,k)          ! 15/07/01
      enddo ; enddo ; enddo

!...........
! xy_t(:,:,1) est une elevation de surface non nulle en cas bathy negative
      if(allocated(hcopy_w)) then !>>>>>>> !22-11-20
! Cas flag_merged_levels==1 oU la vraie bathy est temporairement dans hcopy_w
       do j=0,jmax+1 ; do i=0,imax+1
         xy_t(i,j,1)=max(0.,-hcopy_w(i,j)+0.001*wetdry_cst3) !22-11-20
       enddo         ; enddo
      else                        !>>>>>>>
! Cas flag_merged_levels==0
       do j=0,jmax+1 ; do i=0,imax+1
         xy_t(i,j,1)=max(0.,    -h_w(i,j)+0.001*wetdry_cst3)
       enddo         ; enddo
      endif                       !>>>>>>>

!...........
      do j=0,jmax+1
      do i=0,imax+1
       k=1
       depth_w(i,j,k)=sigma_w(i,j,k)*(h_w(i,j)+xy_t(i,j,1))-h_w(i,j) !22-11-20
       do k=2,kmax+1  ! attention forcement croissante !25-05-15
        depth_w(i,j,k)=max(sigma_w(i,j,k)*(h_w(i,j)+xy_t(i,j,1))-h_w(i,j) & !22-11-20
                          ,depth_w(i,j,k-1)+0.001) ! 25-05-15 !dz/dk>0 !05-06-15
       enddo
      enddo
      enddo

      endif              !----simple----->


!--------------------------------------------------------
!  Generalized (horizontal dependent) sigma distribution:
      if(igesig==1) then !---generalized--->

      ksecu=1

!...........................................................
! ETAPE 1:
! Les profondeurs sont relaxées vers des profondeurs
! "parfaites" fournies par la distribution "à priori" des
! niveaux sigma appliquée à une colonne d'eau dont
! l'épaisseur est spécifiée dans le notbook_technum
! (parametre HGESiG)
!...........................................................
      do k=1,kmax+1
       anyv3d(1,1,k,1)=(sigma_w(1,1,k)-1.)*hgesig
       anyv3d(2,2,k,1)=(sigma_w(1,1,k)-1.)*hmax
       anyv3d(3,3,k,1)=hmin*real(kmax+1-k)/real(1-kmax+1)
      enddo

! xy_t(:,:,1) est une elevation de surface non nulle en cas bathy negative
      if(allocated(hcopy_w)) then !>>>>>>> !22-11-20
! Cas flag_merged_levels==1 oU la vraie bathy est temporairement dans hcopy_w
       do j=0,jmax+1 ; do i=0,imax+1
         xy_t(i,j,1)=max(0.,-hcopy_w(i,j)+0.001*wetdry_cst3) !22-11-20
       enddo         ; enddo
      else                        !>>>>>>>
! Cas flag_merged_levels==0
       do j=0,jmax+1 ; do i=0,imax+1
         xy_t(i,j,1)=max(0.,    -h_w(i,j)+0.001*wetdry_cst3)
       enddo         ; enddo
      endif                       !>>>>>>>

! Details dans: !24-01-17
! https://docs.google.com/document/d/1HyHCalPTAEThRfmvib_4fWA7dMacELrUdSKx2tz2OIw/edit


      rap=pgesig*0.01
      x1=rap*1.+(1.-rap)*real(kmax+1)
      k1=min0(max0( nint( x1 ) , 1 ) ,kmax+1)

      do 170 j=0,jmax+1
      do 170 i=0,imax+1

      do 173 k=k1+1,kmax

      rap=0.
      if( hgesig.le.h_w(i,j)  ) then ! §§§§§§§§§>

       x1=real(k-k1)/real(kmax+1-k1)
       rap=(1.-exp(-x1/0.4))/(1.-exp(-1./0.4)) !31-01-15

      endif                          ! §§§§§§§§§>

      depth_w(i,j,k)=rap *(sigma_w(1,1,k)-1.)*hgesig    &
                +(1.-rap)*(sigma_w(1,1,k)-1.)*h_w(i,j)  &
                          +sigma_w(1,1,k)*xy_t(i,j,1) !07-02-15

  173 continue

      do 190 k=2,k1-1

      rap=0.
      if( hgesig.le.h_w(i,j)  ) then ! §§§§§§§§§>

      rap=((sigma_w(1,1,k )-sigma_w(1,1,k1))                     &
          /(sigma_w(1,1,1 )-sigma_w(1,1,k1)))             !17/04/02 !24-11-10 exposant 1 remplace 0.5

      endif                          ! §§§§§§§§§>

       depth_w(i,j,k)=rap*sigma_w(1,1,k)*hgesig       &
                +(1.-rap)*sigma_w(1,1,k)*h_w(i,j)     &
                -h_w(i,j)                             &
                         +sigma_w(1,1,k)*xy_t(i,j,1)


  190 continue

      k=1
       depth_w(i,j,k)=(sigma_w(1,1,k)-1)*h_w(i,j) &
                      +sigma_w(1,1,k)*xy_t(i,j,1)
      k=k1
       depth_w(i,j,k)=(sigma_w(1,1,k)-1)*h_w(i,j) &
                      +sigma_w(1,1,k)*xy_t(i,j,1)
      k=kmax+1
       depth_w(i,j,k)=(sigma_w(1,1,k)-1)*h_w(i,j) &
                      +sigma_w(1,1,k)*xy_t(i,j,1)

  170 continue



      endif              !---generalized--->


! debug notebook_vertcoord:
      if(ksecu.ne.1) then
      write(6,*)'anonmalie sur valeur coef igesig'
      write(6,*)'on ne sait pas si la grille doit être sigma'
      write(6,*)'ou si la grille doit être sigma generalisee'
      write(6,*)'verifier le notebook_vertcoord'
      stop 'sigma_levels.f'
      endif

!...............................................................................
      if(flag_merged_levels==1) then !m[°v°]> !16-11-17
! Retablir h_w dans sa valeur vraie:
       do j=0,jmax+1 ; do i=0,imax+1
        h_w(i,j)=xy_t(i,j,2)
       enddo         ; enddo
      endif                          !m[°v°]> !16-11-17
      if(flag_merged_levels==1)call sigma_levels_envelope_bathy(2) !15-04-20

      if(kmax>1.and. & !14-04-18
         dzsurfmin>0)call sigma_levels_dzsurfmin ! dz_t(kmax)=min(dz_t(kmax),dzsurfmin) !25-01-18

! Ajuster la position des niveaux depth_w pour ne pas depasser le critere
! d'incoherence hydrostatique defini par Haney (1991)
!     call sigma_levels_hydrocons !05-08-15

      if(par%rank==rankhmax) then !00000000000000000000000000> !15-04-20
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine sigma_levels:'
      write(3,*)
      if(igesig.eq.0) then
      write(3,*)'espacement entre niveaux sigma:'
      write(3,*)'Attention dzsurfmin n''est pas encore pris en compte'
      write(3,*)' niveau | dsig_t | dsig*hmax | z_t(hmax)'
      sum1=0.
      sum2=0.
      do k=kmax,1,-1
      sum1=sum1+dsig_t(1,1,k)*hmax
      write(3,'(i5,4x,f8.5,1x,f9.3,2x,f9.3)') &
        k,dsig_t(1,1,k),dsig_t(1,1,k)*hmax &
         ,-0.5*(sum1+sum2)
      sum2=sum1
      enddo
      endif
      if(igesig.eq.1) then
      write(3,*)'parametres systeme sigma generalisée:'
      write(3,*)'epaisseur de reference:',hgesig
      write(3,*)'partition surface fond:',pgesig
      write(3,*)'exemple de repartition au point le plus profond:'
      write(3,*)'level | zref(m) | zsimplesigma(m) | zgeneralsigma(m)'
      do k=kmax+1,1,-1
      texte11='        '
      if(k.eq.k1)texte11=' <split>'
      write(3,'(i6,f10.3,8x,f10.3,a11,f10.3)')                          &
                                             k                          &
                                            ,anyv3d(1,1,k,1)            &
                                            ,anyv3d(2,2,k,1)            &
                                            ,texte11                    &
                                            ,depth_w(ihmax,jhmax,k)
      enddo
      endif
      close(3)
      endif                !00000000000000000000000000>     !19-07-10

       if(vststep/=0.and.flag_merged_levels==1) &
       stop 'vststep/=0.and.flag_merged_levels==1'

      if(flag_merged_levels==1)call sigma_levels_merged

      end subroutine sigma_levels_depth_w

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!#ifdef bidon
      subroutine sigma_levels_2dv_plot
      use module_principal
      use module_parallele
      implicit none

      if(fgrid_or_wgrid==wgrid_case) then !wwwwwwww>

      write(texte30,'(a,i0)')trim(tmpdirname)//'z2dv_OE_',par%rank
      open(unit=3,file=trim(texte30))
      write(texte30,'(a,i0)')trim(tmpdirname)//'bathy_OE_',par%rank
      open(unit=4,file=trim(texte30))

      j=jmax/2
      do k=1,kmax+1
       if(mod(k,2)==0) then
        i1=1 ; i2=imax ; i3=1
       else
        i1=imax ; i2=1 ; i3=-1
       endif
       do i=i1,i2,i3
        if(i3==1) then !>>>
         write(3,*)i+par%timax(1),depth_w(i,j,k) !07-03-21
         if(k==2)write(4,*)i+par%timax(1),-h_w(i,j)
        else           !>>>
         write(3,*)i+par%timax(1),depth_w(i,j,k) !07-03-21
        endif          !>>>
       enddo
      enddo

      do i=1,imax
       if(mod(i,2)==0) then
        do k=1,kmax+1
          write(3,*)i+par%timax(1),depth_w(i,j,k) !07-03-21
        enddo
       else
        do k=kmax+1,1,-1
          write(3,*)i+par%timax(1), depth_w(i,j,k)
        enddo
       endif
      enddo

      close(3)
      close(4)

      write(texte30,'(a,i0)')trim(tmpdirname)//'z2dv_SN_',par%rank
      open(unit=3,file=trim(texte30))
      write(texte30,'(a,i0)')trim(tmpdirname)//'bathy_SN_',par%rank
      open(unit=4,file=trim(texte30))

       i=imax/2
      do k=1,kmax+1
       if(mod(k,2).EQ.0) THEN
        j1=1
        j2=jmax
        j3=1
       else
        j1=jmax
        j2=1
        j3=-1
       endif
      do j=j1,j2,j3
       if(j3==1) then !>>>
        write(3,*)j+par%tjmax(1),depth_w(i,j,k) !07-03-21
        if(k==2)write(4,*)j+par%tjmax(1),-h_w(i,j)
       else           !>>>
        write(3,*)j+par%tjmax(1),depth_w(i,j,k) !07-03-21
       endif          !>>>
      enddo
      enddo
      do j=1,jmax
       if(mod(j,2)==0) then
        do k=1,kmax+1
          write(3,*)j+par%tjmax(1), depth_w(i,j,k)
        enddo
       else
        do k=kmax+1,1,-1
          write(3,*)j+par%tjmax(1), depth_w(i,j,k)
        enddo
       endif
      enddo

      close(3)

      endif                               !wwwwwwww>

      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09

      end subroutine sigma_levels_2dv_plot
!#endif
!........................................................................
      subroutine initial_depthw_from_file
      use module_principal
      use module_parallele !#MPi
      implicit none

! Cette routine s'appuie sur le fait que dz_t a déjà été initialise dans la routine
! initial_mask_and_bathy !19-03-14

      do j=0,jmax+1
      do i=0,imax+1
       depth_w(i,j,kmax+1)=0.
      enddo
      enddo
      do k=kmax,1,-1
      do j=0,jmax+1
      do i=0,imax+1
       depth_w(i,j,k)=                                          &
       depth_w(i,j,k+1)-max(mask_t(i,j,k)*dz_t(i,j,k,1),1.d-3)
      enddo
      enddo
      enddo

      kmin_u(:,:)=kmax-1
      kmin_v(:,:)=kmax-1
      kmin_w(:,:)=kmax-1
      do k=kmax,1,-1
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,k)==1)kmin_w(i,j)=k
       if(mask_u(i,j,k)==1)kmin_u(i,j)=k
       if(mask_v(i,j,k)==1)kmin_v(i,j)=k
      enddo
      enddo
      enddo

! sigma:
      do k=0,kmax+1
      do j=0,jmax+1
      do i=0,imax+1
       sigma_w(i,j,k)=( depth_w(i,j,k )- depth_w(i,j,kmin_w(i,j)))     & !23/10/06
                     /( depth_w(i,j,kmax+1)- depth_w(i,j,kmin_w(i,j)))
      enddo
      enddo
      enddo

! dsig_t:
      do k=0,kmax
      do j=0,jmax+1
      do i=0,imax+1
       dsig_t(i,j,k)=sigma_w(i,j,k+1)-sigma_w(i,j,k)
      enddo
      enddo
      enddo
      do j=0,jmax+1
      do i=0,imax+1
       dsig_t(i,j,kmax+1)=1.d-6
      enddo
      enddo

! Cas nemoffline:
      if(flag_nemoffline==1) then !11111111111> !30-04-14

! Le facteur d'echelle e3u n'est pas toujours present dans les fichiers nemo
! Dans ce cas dz_u et dz_v de la grille initiale seront utilises tout au long
! du run. On les calcule maintenant:
       do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1
         dz_u(i,j,k,1)=0.5*(dz_t(i,j,k,1)+dz_t(i-1,j,k,1)) !=e3u
       enddo ; enddo ; enddo
! Partial step e3u deduit du partial step de e3t
       do j=0,jmax+1 ; do i=1,imax+1
         dz_u(i,j,kmin_u(i,j),1)=min(dz_t(i  ,j,kmin_u(i,j),1)      &
                                    ,dz_t(i-1,j,kmin_u(i,j),1))
       enddo ; enddo
       do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1
         dz_v(i,j,k,1)=0.5*(dz_t(i,j,k,1)+dz_t(i,j-1,k,1)) !=e3v
       enddo ; enddo ; enddo
! Partial step e3v deduit du partial step de e3t
       do j=1,jmax+1 ; do i=0,imax+1
         dz_v(i,j,kmin_v(i,j),1)=min(dz_t(i,j  ,kmin_v(i,j),1)      &
                                    ,dz_t(i,j-1,kmin_v(i,j),1))
       enddo ; enddo

      endif                       !11111111111>

      end subroutine initial_depthw_from_file

!..................................................................

      subroutine sigma_levels_mpi(txt_)
      use module_principal ; use module_parallele
      implicit none
      character(len=*) :: txt_
      integer loop_

      if(txt_=='sigma_f') then !ssssss>
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange('r1','sigma_f_r1'     &
                                  ,sigma_f         &
                           ,lbound(sigma_f)        &
                           ,ubound(sigma_f)        &
                           ,k0)
        do loop_=1, subcycle_exchange
         call echange_voisin(sigma_f,k0,mpi_neighbor_list(loop_))
        enddo
        call loc_wait() ! ----> important sinon pas d'echanges

       return
      endif                    !ssssss>

      if(txt_=='dsig_f') then  !ssssss>
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange('r1','dsig_f_r1'     &
                                  ,dsig_f         &
                           ,lbound(dsig_f)        &
                           ,ubound(dsig_f)        &
                           ,k0)
        do loop_=1, subcycle_exchange
         call echange_voisin(dsig_f,k0,mpi_neighbor_list(loop_))
        enddo
        call loc_wait() ! ----> important sinon pas d'echanges

       return
      endif                    !ssssss>


      if(txt_=='depth_f') then !ssssss>
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange('r1','depth_f_r1'     &
                                  ,depth_f         &
                           ,lbound(depth_f)        &
                           ,ubound(depth_f)        &
                           ,k0)
        do loop_=1, subcycle_exchange
         call echange_voisin(depth_f,k0,mpi_neighbor_list(loop_))
        enddo
        call loc_wait() ! ----> important sinon pas d'echanges

       return
      endif                    !ssssss>


      stop 'sigma_levels_mpi_f txt_ not recognized'
      end subroutine sigma_levels_mpi

!.................................................................

      subroutine sigma_levels_consistency
      use module_principal ; use module_parallele
      implicit none

       if(par%rank==0)write(6,*)'subroutine sigma_levels_consistency'

       do k=1,kmax+1 ; do j=0,jmax+1 ; do i=0,imax+1
        depth_w(i,j,k)=max(depth_w(i,j,k),-h_w(i,j))
       enddo         ; enddo         ; enddo

! sigma_w !06-06-20
       do k=1,kmax+1   ; do j=0,jmax+1 ; do i=0,imax+1
        sigma_w(i,j,k)=1.-(depth_w(i,j,k)-depth_w(i,j,kmax+1)) &
                         /(depth_w(i,j,1)-depth_w(i,j,kmax+1))
       enddo           ; enddo         ; enddo

!...............
! Securite sur sigma_w !07-11-20
! Eviter les niveaux collEs (note: il est important que sigma_w soit double precision)
       flag_stop=0
       do j=0,jmax+1 ; do i=0,imax+1
        sigma_w(i,j,kmax+1)=1.
        do k=2,kmax+1
! ce seuil sophistique est introduit pour eviter que dz devienne trop petit en zone assechee
         x1=1.e-5/max(h_w(i,j),wetdry_cst3) !02-11-21
         if(sigma_w(i,j,k)-sigma_w(i,j,k-1)<x1)sigma_w(i,j,k)=sigma_w(i,j,k-1)+x1 !02-11-21
!        if(sigma_w(i,j,k)-sigma_w(i,j,k-1)<1.e-9)sigma_w(i,j,k)=sigma_w(i,j,k-1)+1.e-9
        enddo
        if(sigma_w(i,j,kmax+1)>1.) then !debug>
         flag_stop=1
         write(10+par%rank,*)'err 942 sigma_w(i,j,kmax+1)>1'
         write(10+par%rank,*)'par%rank,i,j',par%rank,i,j
         goto 951
        endif                           !debug>
       enddo         ; enddo
 951   call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ;
       if(k0/=0)stop 'err 942 sigma_w(i,j,kmax+1)>1'

! dsig_t
       do k=1,kmax   ; do j=0,jmax+1 ; do i=0,imax+1
        dsig_t(i,j,k)=sigma_w(i,j,k+1)-sigma_w(i,j,k) !06-06-20
       enddo         ; enddo         ; enddo
       dsig_t=max(dsig_t,small1) !06-06-19
! note: retablissement du seuil suite A division par dsig_t(k=0)=0 constatEe



!      call sigma_levels_2dv_plot
!      stop 'bob'

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop ' Stop sigma_levels_consistency see fort.xxx files'

!     stop 'fifi'
      end subroutine sigma_levels_consistency

!...................................................................

      subroutine sigma_levels_merged !16-11-17
      use module_principal ; use module_parallele ; use module_grid
      implicit none

!     allocate(botlevmerged_w(0:imax+1,0:jmax+1)) ; botlevmerged_w=0.5

      flag_stop=0 !22-11-20
      do j=0,jmax+1 ; do i=0,imax+1

! Memoriser depth_w avant fusion:
       do k=1,kmax+1
        anyv1d(k,1)=depth_w(i,j,k)
       enddo

! Eventuellement modifier depth_w pour garantir un minimum de niveau
! par exemple si depth_w(:,:,kmax+1-10)>=-h_w on est certains d'avoir
! au moins 10 niveaux
!      k10=kmax+1-nbvstepmin        !14-02-18
!      if(depth_w(i,j,k10)<-h_w(i,j)) then !m0v0m
!       do k=kmax+1,k10,-1
!        rap=real(kmax+1-k)/real(kmax+1-k10)
!        depth_w(i,j,k)=depth_w(i,j,k)+rap*(-h_w(i,j)-depth_w(i,j,k10))
!       enddo
!      endif                               !m0v0m
! L'algo ci-dessus pouvait dans certains cas conduire A des dz negatifs. Il a
! donc EtE remplacE par l'algo en suivant le 14-12-19
       k10=kmax+1-nbvstepmin  

       if(depth_w(i,j,k10)<-h_w(i,j)) then !m0v0m
        x0=(depth_w(i,j,kmax+1)+h_w(i,j))/(depth_w(i,j,kmax+1)-depth_w(i,j,k10))
        if(x0<=0) then
           flag_stop=1
           write(10+par%rank,*)'Err sigma_levels_merged x0<0 x0=',x0
           write(10+par%rank,*)'i,j loc',i,j
           write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
           write(10+par%rank,*)'depth_w(i,j,kmax+1)',depth_w(i,j,kmax+1)
           write(10+par%rank,*)'depth_w(i,j,k10)',depth_w(i,j,k10),k10
           write(10+par%rank,*)'h_w(i,j)',h_w(i,j)
        endif
        depth_w(i,j,k10)=-h_w(i,j) 
        do k=k10+1,kmax+1
         depth_w(i,j,k)=depth_w(i,j,k-1)+x0*(anyv1d(k,1)-anyv1d(k-1,1)) !14-12-19
        enddo
       endif                               !m0v0m

       kmin_w(i,j)=1
       do k=1,kmax
        if(depth_w(i,j,k)<=-h_w(i,j))kmin_w(i,j)=k !09-02-20
       enddo
       kmin_w(i,j)=min(max(kmin_w(i,j),1),kmax-1)

! Couches fusionnees:
! L'epaisseur minimum des couche fusionnees doit etre beaucoup plus
! petite que la plus fine des couches operationnelles (a priori celle de
! surface) !13-12-18
!      x1=max(min(0.01*(depth_w(i,j,kmax+1)-depth_w(i,j,kmax)),0.001),small1) !13-12-18
!      x1=max(min(0.01*(depth_w(i,j,kmax+1)-depth_w(i,j,kmax)),0.0001),small1) !13-12-18
!      x1=max(min(0.01*(depth_w(i,j,kmax+1)-depth_w(i,j,kmax)),0.0   ),small1) !13-12-18
!      x1=0.
       x1=max(small1*(depth_w(i,j,kmax+1)-depth_w(i,j,1)) & ! pour que dsig_t soit small1 dans les couches ecrasees
             ,small1) 
       depth_w(i,j,1)=-h_w(i,j)
       do k=2,kmin_w(i,j)
         depth_w(i,j,k)=depth_w(i,j,k-1)+x1
       enddo
! mask_t = 0 si dz=0
!      do k=1,kmin_w(i,j)-1 ! supprimE tant que mpi non fait
!       mask_t(i,j,k)=0
!      enddo

! Au dessus des couches fusionnees:eviter les couches trop minces
       x1=depth_w(i,j,kmin_w(i,j))-depth_w(i,j,1) !07-12-17
       do k=kmin_w(i,j)+1,kmax+1
        depth_w(i,j,k)=depth_w(i,j,k)+x1 ! reporter translation depth_w(i,j,kmin_w(i,j))!06-12-17
        if(depth_w(i,j,k)-depth_w(i,j,k-1)<0.001) then
           depth_w(i,j,k)=depth_w(i,j,k-1)+0.001
        endif
       enddo

! 09-02-20 nouvelle strategie: par defaut kmerged=kmin+1 partout sauf...
! les cas particuliers avec k kmerged_t(i,j)=kmin_w(i,j) sont:
! 1- kmin_w(i,j)=kmax+1-nbvstepmin
! 2- kmin_w(i,j)=1 et dz(kmerged)<dz(kmin)
       k=kmin_w(i,j)
       kmerged_t(i,j)=kmin_w(i,j)+1
       if(kmin_w(i,j)==kmax+1-nbvstepmin)kmerged_t(i,j)=kmin_w(i,j) 
       if(kmin_w(i,j)==1.and.  &
            (depth_w(i,j,k+2)-depth_w(i,j,k+1))       &  !09-02-20
            <depth_w(i,j,k+1)-depth_w(i,j,k  ))kmerged_t(i,j)=kmin_w(i,j)




! En principe lorsque kmin_w=k10 (le seuil minimum) on s'attend A ce que kmerged=kmin
! sauf si un pb de precision machine a fait echouer les aiguillages mis en
! place auparavant. On verifie ce point et lance une alerte si ce n'est pas le cas
       if(kmin_w(i,j)==kmax+1-nbvstepmin.and.kmerged_t(i,j)/=kmin_w(i,j)) then !ooo> !09-02-20
         write(10+par%rank,*)'PB dans sigma_levels.F90'
         write(10+par%rank,*)'ANORMAL: ',  &
      'kmin_w(i,j)=kmax+1-nbvstepmin mais est different de kmerged_t'
         write(10+par%rank,*)'kmax+1-nbvstepmin',kmax+1-nbvstepmin
         write(10+par%rank,*)'kmin_w(i,j)',kmin_w(i,j)
         write(10+par%rank,*)'kmerged_t(i,j)',kmerged_t(i,j)
         write(10+par%rank,*)'i,j,par%rank',i,j,par%rank
         write(10+par%rank,*)'i,j,glob',i+par%timax(1),j+par%tjmax(1)
         write(10+par%rank,*)'dz(kmin)   ',depth_w(i,j,kmin_w(i,j)+1)-depth_w(i,j,kmin_w(i,j))
         write(10+par%rank,*)'dz(kmerged)',depth_w(i,j,kmerged_t(i,j)+1)-depth_w(i,j,kmerged_t(i,j))
       endif                                                                   !ooo> !09-02-20

! Enfin on s'assure que la couche kmin n'est pas "infiniment" (=trop) petite 
! par rapport A la couche kmerged !09-02-20
      k=kmin_w(i,j)
      if(depth_w(i,j,k+1)-depth_w(i,j,k)<0.01*(depth_w(i,j,k+2)-depth_w(i,j,k+1))) then
       depth_w(i,j,k+1)=(depth_w(i,j,k)+0.01*depth_w(i,j,k+2))/(1.+0.01)
      endif
! Message d'alerte
      if(depth_w(i,j,k+1)-depth_w(i,j,k)<0.009*(depth_w(i,j,k+2)-depth_w(i,j,k+1))) then
       write(10+par%rank,*)'ATTENTION dz(kmin)<0.009*dz(kmerged)'
       write(10+par%rank,*)'par%rank',par%rank
       write(10+par%rank,*)'i,j loc',i,j
       write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
       write(10+par%rank,*)'dz(kmin   )=',depth_w(i,j,k+1)-depth_w(i,j,k)
       write(10+par%rank,*)'dz(kmerged)=',depth_w(i,j,k+2)-depth_w(i,j,k+1)
      endif


      enddo         ; enddo ! i,j

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop 'STOP sigma_levels_merged, see fort.xxx error files'

! kmin_u kmin_v kmerged_u kmerged_v kundermin_t kundermin_u kundermin_v upwzone0_t
      call sigma_levels_merged_finalize

      end subroutine sigma_levels_merged

!..................................................................

      subroutine sigma_levels_merged_finalize     ! 26-03-18
      use module_principal ; use module_parallele ! use module_grid
      implicit none

! Deduire kmin_u kmin_v kmerged_u kmerged_v kundermin_t kundermin_u kundermin_v upwzone0_t
! de kmin_w et kmerged_t

      do i=1,imax+1 ; do j=0,jmax+1
          kmin_u(i,j)=min(   kmin_w(i,j),   kmin_w(i-1,j))
!      kmerged_u(i,j)=min(kmerged_t(i,j),kmerged_t(i-1,j))        
       kmerged_u(i,j)=max(kmerged_t(i,j),kmerged_t(i-1,j))
      enddo         ; enddo

      do i=0,imax+1 ; do j=1,jmax+1
          kmin_v(i,j)=min(   kmin_w(i,j),   kmin_w(i,j-1))
!      kmerged_v(i,j)=min(kmerged_t(i,j),kmerged_t(i,j-1))        
       kmerged_v(i,j)=max(kmerged_t(i,j),kmerged_t(i,j-1))
      enddo         ; enddo

! kundermin_t est le niveau "traceur" le plus bas ou on trouve le premier
! flux lateral non nul. Sous ce niveau il n'y a plus d'echange. 
      if(.not.allocated(kundermin_t))allocate(kundermin_t(imax,jmax))
      if(.not.allocated(kundermin_u))allocate(kundermin_u(2:imax,2:jmax-1))
      if(.not.allocated(kundermin_v))allocate(kundermin_v(2:imax-1,2:jmax))
      kundermin_t=1
      kundermin_u=1
      kundermin_v=1
      do i=1,imax   ; do j=1,jmax
       kundermin_t(i,j)=min(min(min(kmin_u(i,j),kmin_u(i+1,j)),kmin_v(i,j)),kmin_v(i,j+1))
      enddo         ; enddo

      do j=2,jmax-1 ; do i=2,imax   ! u
! l'advection u de u considere u(i-1:i+1,j)
       kundermin_u(i,j)=min(min(kmin_u(i,j),kmin_u(i+1,j)),kmin_u(i-1,j))
! l'advection v de u considere v(i,j) v(i,j+1) v(i-1,j) v(i-1,j+1)
       kundermin_u(i,j)=min(min(min(kmin_v(i,j),kmin_v(i,j+1)),kmin_v(i-1,j)),kmin_v(i-1,j+1))
      enddo         ; enddo

      do j=2,jmax   ; do i=2,imax-1 ! v
! l'advection v de v considere v(i,j-1:j+1)
       kundermin_v(i,j)=min(min(kmin_v(i,j),kmin_v(i,j+1)),kmin_v(i,j-1))
! l'advection v de v considere u(i,j) u(i+1,j) u(i,j-1) u(i+1,j-1)
       kundermin_v(i,j)=min(min(min(kmin_u(i,j),kmin_u(i+1,j)),kmin_u(i,j-1)),kmin_u(i+1,j-1))
      enddo         ; enddo

! Ces lignes sont supprimees le 21-05-20 car on se sert desormais differemment de upwzone0_t et ce apres avoir remarquE que l'usage qui
! en ete fait dans l'advection verticale de T,S etait redondant du fait que la vitesse verticale est 100% implicite de 1 a kmerged
!     upwzone0_t=100 ! 100 est la valeur standard (schema advection normal) - 0 = upwind 
!     do j=0,jmax+1
!     do i=0,imax+1
!      upwzone0_t(i,j,1:kmin_w(i,j))=0 ! 0==>100%upwind (100=valeur normale standard)
!     enddo
!     enddo

! Voir turbulence_adv pour comprendre l'utilite des lignes suivantes:
      mask_vqs_tke_w=1 !28-12-22
      if(flag_merged_levels==1) then !pmxpmx> 
       do j=2,jmax-1 ; do i=2,imax-1
        do k=2,kmerged_t(i,j)-1
         mask_vqs_tke_w(i,j,k)=0
        enddo
       enddo         ; enddo
      endif                          !pmxpmx> 

      end subroutine sigma_levels_merged_finalize

!.............................................

      subroutine sigma_levels_dzsurfmin !25-01-18 !27-01-18
      use module_principal ; use module_parallele ! use module_grid
      implicit none

! Preserve sufficient vertical resolution near the surface: dz(kmax)=min(dz(kmax),dzsurmin)


!     i=imax/2 ; j=jmax/2
!     write(10+par%rank,*)'---------------------------'
!     do k=kmax,1,-1
!      write(10+par%rank,*)k,depth_w(i,j,k),depth_w(i,j,k+1)-depth_w(i,j,k)
!     enddo
      do j=0,jmax+1 ; do i=0,imax+1

! x0 correction A apporter A dz en surface pour que dz <= dzsurfmin
         k=kmax
         x0=max(depth_w(i,j,k+1)-depth_w(i,j,k)-dzsurfmin,0.)
! loi delta_dz(k)=x1+x2*exp( (k-kmax)/x3 ) telle que sum(delta_dz(:))=0 et delta_dz(kmax)=x0
         sum1=0. ; sum2=0.  ; x3=5. ! x3=5. => ajustement des niveaux inferieurs jusqu'a 5m sous la surface
         do k=1,kmax
          sum1=sum1+1.
          sum2=sum2+exp(real(k-kmax)/x3)
         enddo
         x2=-sum1*x0/(sum2-sum1) ; x1=x0-x2
! anyv1d(k,1)=dz(k)+delta_dz(k)
         sum3=0.
         do k=kmax,1,-1
          anyv1d(k,1)=depth_w(i,j,k+1)-depth_w(i,j,k)-x1-x2*exp(real(k-kmax)/x3)
          sum3=sum3+x1+x2*exp(real(k-kmax)/x3)
!         if(i==imax/2.and.j==jmax/2)write(10+par%rank,*)'sum3',k,sum3,anyv1d(k,1)
         enddo

         do k=kmax,1,-1
          depth_w(i,j,k)=depth_w(i,j,k+1)-anyv1d(k,1)
         enddo

      enddo         ; enddo

!     i=imax/2 ; j=jmax/2
!     write(10+par%rank,*)'---------------------------'
!     do k=kmax,1,-1
!      write(10+par%rank,*)k,depth_w(i,j,k),depth_w(i,j,k+1)-depth_w(i,j,k)
!     enddo
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'bibi'

      end subroutine sigma_levels_dzsurfmin

!.............................................

      subroutine sigma_levels_envelope_bathy(case_) !15-04-20
      use module_principal ; use module_parallele 
      implicit none
      real , dimension(:,:,:) , allocatable :: h2envelope
      real , dimension(:)   , allocatable :: lon_envelope &
                                            ,lat_envelope &
                                           ,dist_envelope &
                                         ,weight_envelope
      integer :: dim1_envelope=0 , dim2_envelope=0 , case_ &
                , k1_,k2_

      if(case_==1) then !111111111111111111111111111>

       allocate(hcopy_w(0:imax+1,0:jmax+1))

      
       do j=0,jmax+1 ; do i=0,imax+1
! Faire une copie h_w avant de le bidouiller
        hcopy_w(i,j)=h_w(i,j)
       enddo ; enddo

       if(hstepmax<0.)hstepmax=hmax           !27-01-18
       hstepmax=min(hstepmax,hmax)            !14-02-18

      if(vqs_file=='none') then !-Fortran->

! Ces constantes sont pour avoir une coordonnee 100% sigma pour 0m < h < vqs_cst2 (vqs_cst2 typiquement entre 20m et 50m). !08-11-20
! Pour l'utilisateur:
! vqs_cst1: bathy envelope quand la vraie bathy =0. Il faut l'augmenter pour diminuer le nbre de niveaux sur le plateau
! vqs_cst2: valeur de la vraie bathy en dessous de laquelle on est 100% sigma
! vqs_cst3: est en facteur du terme de courbure de la bathy enveloppe. On joue sur celui ci pour eviter par exemple 
!           que niveaux et bathy aient des pentes de signes opposes
! Pour le developpeur:
! x1 : est la vraie bathy quand celle-ci = vqs_cst2
! x2 : la valeur de la bathy envelope quand la vraie bathy est egale x1=vqs_cst2
! lignes commentees le 29-04-21
!      vqs_cst1=200  ! 300. !400. !22-03-21
!      vqs_cst2=100. ! 60.        !22-03-21
!      vqs_cst3=0.5  ! 05-03-21   !22-03-21
       x1=vqs_cst2  
       x2=max(x1+vqs_cst3*x1*(1.-(x1/hstepmax))+vqs_cst1*(hstepmax-x1)/hstepmax,x1) !22-03-21
       do j=0,jmax+1 ; do i=0,imax+1

        if(hcopy_w(i,j)>vqs_cst2) then !>>>>
! si H (vraie) > vqs_cst2 alors coord vqs
! Attention si on modifie cette formule alors modifier aussi x2!
         h_w(i,j)=max(hcopy_w(i,j)+vqs_cst3*hcopy_w(i,j)*(1.-(hcopy_w(i,j)/hstepmax))+vqs_cst1*(hstepmax-hcopy_w(i,j))/hstepmax,hcopy_w(i,j)) !15-04-20
        else                     !>>>>
! si H (vraie) < vqs_cst2 (faibles profondeurs) alors coord sigma !08-11-20
         h_w(i,j)=max(hcopy_w(i,j)*x2/x1,0.)
        endif                    !>>>>

       enddo ; enddo

! Lisser l'envelope pour lisser la discontinuite entre la zone vqs et la  zone sigma
       call lissebathy(1,0.,1,iglb,1,jglb,10)
      endif                    !-Fortran->

      if(vqs_file/='none') then !-Fichier->

      if(par%rank==0) then !00000>
        dim1_envelope=0
        open(unit=3,file=trim(vqs_file))
        read(3,*) ! ligne pour rien
!       read(3,*) ! ligne pour rien
!1329   read(3,*,end=1330)x1,x2
!       dim1_envelope=dim1_envelope+1
!       goto 1329
        read(3,*)dim1_envelope,dim2_envelope
 1330   close(3)
!       write(6,*)'dim1_envelope=',dim1_envelope
      endif                !00000>

      call mpi_bcast(dim1_envelope,1,mpi_integer,0,par%comm2d,ierr)
      call mpi_bcast(dim2_envelope,1,mpi_integer,0,par%comm2d,ierr)

      allocate(     h2envelope(dim1_envelope,dim2_envelope,2)) ; h2envelope=0.
      allocate(   lon_envelope(dim2_envelope))                 ; lon_envelope=0.
      allocate(   lat_envelope(dim2_envelope))                 ; lat_envelope=0.
      allocate(  dist_envelope(dim2_envelope))                 ; dist_envelope=0.
      allocate(weight_envelope(dim2_envelope))                 ; weight_envelope=0.

      if(par%rank==0) then !00000>
        open(unit=3,file=trim(vqs_file))
          read(3,*) ! ligne https....
          read(3,*) ! ligne contenant dim1_envelope dim2_envelope

          do k2=1,dim2_envelope
          read(3,*) ! ligne pointillee de separation
          read(3,*)lon_envelope(k2),lat_envelope(k2),dist_envelope(k2),weight_envelope(k2)
          dist_envelope(k2)=dist_envelope(k2)*1000. ! conversion km en m
           do k1=dim1_envelope,1,-1 
!           write(6,*)'k2,k1,dim1_envelope',k2,k1,dim1_envelope
            read(3,*)h2envelope(k1,k2,1),h2envelope(k1,k2,2)
!          write(6,*)'h2envelope(k1,k2,1),h2envelope(k1,k2,2)',h2envelope(k1,k2,1),h2envelope(k1,k2,2)
           enddo
          enddo
        close(3)
      endif                !00000>
      call mpi_bcast(h2envelope,dim1_envelope*dim2_envelope*2,mpi_real,0,par%comm2d,ierr)
      call mpi_bcast(lon_envelope,dim2_envelope,mpi_real,0,par%comm2d,ierr)
      call mpi_bcast(lat_envelope,dim2_envelope,mpi_real,0,par%comm2d,ierr)
      call mpi_bcast(dist_envelope,dim2_envelope,mpi_real,0,par%comm2d,ierr)
      call mpi_bcast(weight_envelope,dim2_envelope,mpi_real,0,par%comm2d,ierr)

!......................................................................
! Conversions des % en bathy enveloppe
! Si h2envelope(k1_,2)<0 simplement changer le signe
! Si h2envelope(k1_,2)>0 il s'agit d'un % A convertir en bathy enveloppe
      do k2_=1,dim2_envelope
      do k1_=1,dim1_envelope

! Si h2envelope(k1_,k2_,2)<0 simplement changer le signe
       if(h2envelope(k1_,k2_,2)<0.)  then !pmx>

          h2envelope(k1_,k2_,2)=-h2envelope(k1_,k2_,2)

! Il n'y a pas de raison que la bathy envelope soit plus profonde que hmax: !23-12-21
          h2envelope(k1_,k2_,2)=min(h2envelope(k1_,k2_,2),hmax) !23-12-21

! Si h2envelope(k1_,k2_,2)>0 il s'agit d'un % A convertir en bathy enveloppe
       else                          !pmx>

! Le nombre de niveau sera le % de kmax+1:
         x4=h2envelope(k1_,k2_,2)/100.*real(kmax+1)
! Le numero du niveau le plus bas (on arrondit au plus proche) sera
! kmax+1 moins le nombre de niveau:
         k=nint( real(kmax+1)-x4)

! On veut que la profondeur de ce niveau soit le fond. Quelle bathy
! envelope donne ce resultat?
! x3 first guess la bathy "vraie" (h2envelope(k1_,k2_,1)
         x3=h2envelope(k1_,k2_,1)
! Compte tenu de k et x3 la routine suivante calcule x2 la profondeur de la
! formule sigma generalisee:
         call sigma_levels_k_2_zw

! Ce premier resultat est archivE dans x22:
         x33=x3 ; x22=x2
! Apportons une perturbation A x3:
         x3=x33*1.1
         call sigma_levels_k_2_zw

         do loop1=1,4
! La derivee dz/dhe
          if(abs(x3-x33)<0.01)goto 1418
          x0=(x2-x22)/(x3-x33)
          x33=x3 ; x22=x2
! Prochain guess:
          x3=x3+(-h2envelope(k1_,k2_,1)-x2)/x0
          call sigma_levels_k_2_zw
         enddo

 1418    h2envelope(k1_,k2_,2)=x3

!..........................................
! Pour VERIFICATION VALIDATION (commenter sinon)
! x3 la bathy "vraie" sera h2envelope(k1_,k2_,1)
         x3=h2envelope(k1_,k2_,2)
         call sigma_levels_k_2_zw
       if(par%rank==0.) then
            write(666,*)'VALIDATION'
            write(666,*)'Si x3     =',x3
            write(666,*)'Alors -x2 =',-x2
            write(666,*)'Est ce  = ?',h2envelope(k1_,k2_,1)
       endif
!..........................................

       endif                         !pmx>


      enddo ! k1_=1,dim1_envelope
      enddo ! k2_=1,dim2_envelope


!..............................................
! Interpoller les valeurs de la bathy enveloppe
       do j=0,jmax+1 ; do i=0,imax+1

        sum1=small1
        sum2=0.


        do k2=1,dim2_envelope

          x1=lon_envelope(k2)*deg2rad
          x2=lat_envelope(k2)*deg2rad
          call lonlat2distance(lon_t(i,j),lat_t(i,j),x1,x2,dist)
! poid distance:
          x0=weight_envelope(k2)*exp(-(dist/dist_envelope(k2))**2)

!        if(i+par%timax(1)==210.and.j+par%tjmax(1)==155) &
!        write(10+par%rank,*)'Gibraltar',k2,real(x0)

!        if(i+par%timax(1)==518.and.j+par%tjmax(1)==1060) &
!        write(10+par%rank,*)'Adri',k2,real(x0)

!        if(i+par%timax(1)==464.and.j+par%tjmax(1)==575) &
!        write(10+par%rank,*)'Sicile',k2,real(x0)

!        if(i+par%timax(1)==1526.and.j+par%tjmax(1)==511) &
!        write(10+par%rank,*)'Egypt',k2,real(x0)

          do k1=2,dim1_envelope

           if(hcopy_w(i,j)>=h2envelope(k1-1,k2,1)  &
         .and.hcopy_w(i,j)<=h2envelope(k1,k2,1)) then !PMXPMX>

              rap=( hcopy_w(i,j)   -h2envelope(k1-1,k2,1)) &
              /(h2envelope(k1,k2,1)-h2envelope(k1-1,k2,1))
!            h_w(i,j)=(1.-rap)*h2envelope(k1-1,k2,2)+rap*h2envelope(k1,k2,2)
            sum1=sum1+x0
            sum2=sum2+x0*( (1.-rap)*h2envelope(k1-1,k2,2)+rap*h2envelope(k1,k2,2) )

!        if(i+par%timax(1)==778.and.j+par%tjmax(1)==779) then
!            write(280,*)'coucou',k2,k1,real(x0),real( (1.-rap)*h2envelope(k1-1,k2,2)+rap*h2envelope(k1,k2,2) )
!        endif

           endif                                      !PMXPMX>

          enddo ! k1
!        if(hcopy_w(i,j)>h2envelope(dim1_envelope,k2,1))h_w(i,j)=h2envelope(dim1_envelope,k2,2)
!        if(hcopy_w(i,j)<h2envelope(1            ,k2,1))h_w(i,j)=h2envelope(1            ,k2,2)
         if(hcopy_w(i,j)>h2envelope(dim1_envelope,k2,1)) then !>>>
            sum1=sum1+x0
            sum2=sum2+x0*h2envelope(dim1_envelope,k2,2)
         endif                                                !>>>
         if(hcopy_w(i,j)<h2envelope(1            ,k2,1)) then !>>>
            sum1=sum1+x0
            sum2=sum2+x0*h2envelope(1            ,k2,2)
         endif                                                !>>>

        enddo ! k2

!        h_w(i,j)=max(h_w(i,j),hcopy_w(i,j))
         h_w(i,j)=max(sum2/sum1,hcopy_w(i,j))

       enddo ; enddo

      deallocate(h2envelope)
      deallocate(lon_envelope)
      deallocate(lat_envelope)
      deallocate(dist_envelope)
      deallocate(weight_envelope)

! Lisser l'envelope (pour lisser les discontinuites de la procedure fichier)
      call lissebathy(1,0.,1,iglb,1,jglb,20) !23-12-21

!     if(par%rank==280)write(666,*)'diantre',h_w(imax/2,jmax/2),hcopy_w(imax/2,jmax/2)

      endif                    !-Fichier->

      endif             !111111111111111111111111111>

      if(case_==2) then !222222222222222222222222222>
! Retablir h_w dans sa valeur vraie:
       do j=0,jmax+1 ; do i=0,imax+1
        h_w(i,j)=hcopy_w(i,j)
!       if(i+par%timax(1)==901.and.j+par%tjmax(1)==514)write(10+par%rank,*)'B',hcopy_w(i,j),h_w(i,j)
       enddo         ; enddo
       deallocate(hcopy_w)
      endif             !222222222222222222222222222>

      end subroutine sigma_levels_envelope_bathy

!.............................................

      subroutine sigma_levels_k_2_zw !15-04-20
      use module_principal ; use module_parallele
      implicit none

! Entree:
! k: niveau vertical "w"
! x3: bathymetrie

! Sortie
! x2: la profondeur depth_w correspondant A k et x3


! Cette routine considere en entree k, un numero de niveau vertical,
! x3 une valeur de bathymetrie,
! et calcule depth_w (x2) la profondeur correspondant A k et x3 compte
! tenu de la formule sigma generalisee du modele.

      rap=pgesig*0.01
      x1=rap*1.+(1.-rap)*real(kmax+1)
      k1=min0(max0( nint( x1 ) , 1 ) ,kmax+1)

!     do 173 k=k1+1,kmax
      if(k>=k1+1) then !>>>>>

      rap=0.
      if( hgesig.le.x3  ) then ! §§§§§§§§§>

       x1=real(k-k1)/real(kmax+1-k1)
       rap=(1.-exp(-x1/0.4))/(1.-exp(-1./0.4)) !31-01-15

      endif                          ! §§§§§§§§§>

      x2=rap *(sigma_w(1,1,k)-1.)*hgesig    &
                +(1.-rap)*(sigma_w(1,1,k)-1.)*x3  &
                          +sigma_w(1,1,k)*max(0., -x3+0.001*wetdry_cst3)

      endif            !>>>>>


! 173 continue

!     do 190 k=2,k1-1
      if(k<=k1-1) then !>>>>>

      rap=0.
      if( hgesig.le.x3  ) then ! §§§§§§§§§>

      rap=((sigma_w(1,1,k )-sigma_w(1,1,k1))                     &
          /(sigma_w(1,1,1 )-sigma_w(1,1,k1)))             !17/04/02 !24-11-10 exposant 1 remplace 0.5

      endif                          ! §§§§§§§§§>

       x2=rap*sigma_w(1,1,k)*hgesig &
                +(1.-rap)*sigma_w(1,1,k)*x3     &
                -x3                             &
                         +sigma_w(1,1,k)*max(0., -x3+0.001*wetdry_cst3)

      endif            !>>>>>


! 190 continue

      if(k==1)                              & 
       x2=(sigma_w(1,1,k)-1)*x3 &
                      +sigma_w(1,1,k)*max(0., -x3+0.001*wetdry_cst3)
      if(k==k1)                             & 
       x2=(sigma_w(1,1,k)-1)*x3 &
                      +sigma_w(1,1,k)*max(0., -x3+0.001*wetdry_cst3)
      if(k==kmax+1)                         & 
       x2=(sigma_w(1,1,k)-1)*x3 &
                      +sigma_w(1,1,k)*max(0., -x3+0.001*wetdry_cst3)

      end subroutine sigma_levels_k_2_zw

!.................................................................

      subroutine sigma_levels_kmergedr4 !04-06-20
      use module_principal ; use module_parallele
      implicit none

! La partie entiere de kmergedr4 = kmerge (sauf cas particuliers avec kmerged-kmin>1)
! La partie decimale de kmergedr4 indique la fraction de la couche kmerge qui intersecte la couche de fond
! definie comme la couche entre le fond et une hauteur egale A dz(kmerged) cette derniere etant
! continue sur l'horizontale
! Details dans:
! https://docs.google.com/presentation/d/1aUGTYvHSOEZIDCHzs-EBPKCYGJ9cXM1ntRLuaw9Ra0w/edit?folder=0BxfDfpz8eP5VU2RKRGZjU01IMkU#slide=id.g6c10187b84_0_1


! kmergedr4_u !05-12-19
       do j=0,jmax+1 ; do i=1,imax+1
          k1=kmerged_u(i,j)
! Cette initialisation convient au cas particulier avec kmin=kmerged: 
! car on obtient dans ce cas kmergedr4=kmin+1 comme ca le calcul de velbot 
! prendra 100% de vel(kmin) + 0% de vel(kmin+1) 
          kmergedr4_u(i,j)=kmin_u(i,j)+1
          sum1=0.
          do k=kmin_u(i,j),kmerged_u(i,j)
           sum1=sum1+dz_u(i,j,k,1)
           if(sum1>dz_u(i,j,k1,1)) then !m°v°m>
            kmergedr4_u(i,j)=k+(dz_u(i,j,k,1)-(sum1-dz_u(i,j,k1,1)) )/dz_u(i,j,k,1)
            goto 773
           endif                        !m°v°m>
          enddo
  773  continue
!ligne suivante pour avoir kmergedr4_u(i,j)=kmin_u(i,j)+1 si jamais dz_u(kmin_u)>dz_u(kmerged_u) :
        kmergedr4_u(i,j)=max(kmergedr4_u(i,j),real(kmin_u(i,j)+1)) !25-10-21
!      write(10+par%rank,*)real(xy_u(i,j,1)),kmergedr4_u(i,j),' U'
       enddo         ; enddo

! kmergedr4_v !05-12-19
       do j=1,jmax+1 ; do i=0,imax+1
          k1=kmerged_v(i,j)
          kmergedr4_v(i,j)=kmin_v(i,j)+1
          sum1=0.
          do k=kmin_v(i,j),kmerged_v(i,j)
           sum1=sum1+dz_v(i,j,k,1)
           if(sum1>dz_v(i,j,k1,1)) then !m°v°m>
            kmergedr4_v(i,j)=k+(dz_v(i,j,k,1)-(sum1-dz_v(i,j,k1,1)))/dz_v(i,j,k,1)
            goto 774
           endif                        !m°v°m>
          enddo
  774  continue
!ligne suivante pour avoir kmergedr4_v(i,j)=kmin_v(i,j)+1 si jamais dz_v(kmin_v)>dz_v(kmerged_v) :
        kmergedr4_v(i,j)=max(kmergedr4_v(i,j),real(kmin_v(i,j)+1)) !25-10-21
!      write(10+par%rank,*)real(xy_v(i,j,1)),kmergedr4_v(i,j),' V'
       enddo         ; enddo
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'TOTO' 

      if(flag_merged_levels==1) then !m°v°m> !02-11-21

! dsigmerged_u
      do j=1,jmax ; do i=1,imax+1
       dsigmerged_u(i,j)=1. ! par defaut !26-11-21
       if(kmerged_u(i,j)>kmin_u(i,j)) then !pmx> !26-11-21
        k=kmerged_u(i,j)-1
!      if(dsig_t(i-1,j,k)<dsig_t(i,j,k)) then !ooo>
!        x1=0.5*(dsig_t(i-1,j,k)+dsig_t(i,j,k)+dsig_t(i-1,j,k+1)+dsig_t(i,j,k+1)) &
!              /(dsig_t(i-1,j,k)+dsig_t(i-1,j,k+1))  
!        dsigmerged_u(i,j)=x1*dsig_t(i-1,j,k)
!      else                                   !ooo>
!        x1=0.5*(dsig_t(i-1,j,k)+dsig_t(i,j,k)+dsig_t(i-1,j,k+1)+dsig_t(i,j,k+1)) &
!              /(dsig_t(i  ,j,k)+dsig_t(i  ,j,k+1))  
!        dsigmerged_u(i,j)=x1*dsig_t(i  ,j,k)
!      endif                                  !ooo>
! https://docs.google.com/document/d/1NyusPr-6GD47EBXe_m13sj_D2Ux1BpEPOQQ9XkdYz_I/edit?usp=sharing
       dsigmerged_u(i,j)=min(1., & !17-11-21
         dz_t(i-1,j,k,1)*dz_u(i,j,k+1,1)/(dz_t(i-1,j,k+1,1)*dz_u(i,j,k,1)) &
        ,dz_t(i  ,j,k,1)*dz_u(i,j,k+1,1)/(dz_t(i  ,j,k+1,1)*dz_u(i,j,k,1)))

       endif                               !pmx> !26-11-21
      enddo ; enddo

! pgfratio_u !03-11-21
      do j=1,jmax ; do i=1,imax+1
       if(kmin_u(i,j)/=kmerged_u(i,j)) then !m°v°m> !27-12-21
        k=kmerged_u(i,j)-1
        pgfratio_u(i,j)=min(1.,dz_u(i,j,k  ,1) &
                              /dz_u(i,j,k+1,1))
       else                                 !m°v°m> !27-12-21
        pgfratio_u(i,j)=1.
       endif                                !m°v°m> !27-12-21
      enddo ; enddo

! dsigmerged_v
      do j=1,jmax+1 ; do i=1,imax
       dsigmerged_v(i,j)=1. ! par defaut !26-11-21
       if(kmerged_v(i,j)>kmin_v(i,j)) then !pmx> !26-11-21
        k=kmerged_v(i,j)-1
!      if(dsig_t(i,j-1,k)<dsig_t(i,j,k)) then !ooo>
!        x1=0.5*(dsig_t(i,j-1,k)+dsig_t(i,j,k)+dsig_t(i,j-1,k+1)+dsig_t(i,j,k+1)) &
!              /(dsig_t(i,j-1,k)+dsig_t(i,j-1,k+1))  
!        dsigmerged_v(i,j)=x1*dsig_t(i,j-1,k)
!      else                                   !ooo>
!        x1=0.5*(dsig_t(i,j-1,k)+dsig_t(i,j,k)+dsig_t(i,j-1,k+1)+dsig_t(i,j,k+1)) &
!              /(dsig_t(i  ,j,k)+dsig_t(i  ,j,k+1))  
!        dsigmerged_v(i,j)=x1*dsig_t(i  ,j,k)
!      endif                                  !ooo>
       dsigmerged_v(i,j)=min(1., & !17-11-21
         dz_t(i,j-1,k,1)*dz_v(i,j,k+1,1)/(dz_t(i,j-1,k+1,1)*dz_v(i,j,k,1)) &
        ,dz_t(i,j  ,k,1)*dz_v(i,j,k+1,1)/(dz_t(i,j  ,k+1,1)*dz_v(i,j,k,1)))

       endif                               !pmx> !26-11-21
      enddo ; enddo

! pgfratio_v !03-11-21
      do j=1,jmax+1 ; do i=1,imax
       if(kmin_v(i,j)/=kmerged_v(i,j)) then !m°v°m> !27-12-21
        k=kmerged_v(i,j)-1
        pgfratio_v(i,j)=min(1.,dz_v(i,j,k  ,1) &
                              /dz_v(i,j,k+1,1))
       else                                 !m°v°m> !27-12-21
        pgfratio_v(i,j)=1.
       endif                                !m°v°m> !27-12-21
      enddo ; enddo

      else                           !m°v°m>

! cas flag_merged_levels=0 :
       do j=1,jmax ; do i=1,imax+1
         dsigmerged_u(i,j)=0.5*(dsig_t(i,j,1)+dsig_t(i-1,j,1))
       enddo       ; enddo
       do j=1,jmax+1 ; do i=1,imax
         dsigmerged_v(i,j)=0.5*(dsig_t(i,j,1)+dsig_t(i,j-1,1))
       enddo       ; enddo

      endif                          !m°v°m>

      end subroutine sigma_levels_kmergedr4 !04-06-20
