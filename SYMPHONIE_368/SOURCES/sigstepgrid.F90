      subroutine sigstepgrid_driver(case_)
!______________________________________________________________________
! S model
! release S26 - last update: 05-08-15
!______________________________________________________________________
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
! S26     28-11-14  Ne pas passer en cas de fgrid (jusqu'a nouvel ordre)
!         09-01-15  Affichage ecran
!         28-01-15  Amelioration de l'algo d'evitement du chevauchement
!         02-02-15  sigma/step construite à partir du lissage de H
!         17-03-15  nbvstepmin (lu dans notebook_vertcoord.f) nombre de niveaux verticaux
!                   minimum 
!         24-05-15  Ajout d'une verification des dimensions du fichier sigstepgrifile
!         25-05-15  La partie de la grille qui n'est pas VST (celle ou le nombre
!                   de niveau est impose) doit etre lisse vis a vis du critere rmax...
!         21-07-15  ecrire le masque dans le fichier sigstep
!         05-08-15  des subroutines (commentees) pour des alternatives a la methode
!                   vst habituelle
!...............................................................................
      use module_principal ; use module_parallele
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='sigstepgrid_driver'
       subroutinedescription= &
       'sigstepgrid subroutine driver'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! La subroutine est appelee sigma_levels

      if(case_.eq.0)return
      if(i2dh.eq.0) then                                               !21/03/07
      write(6,*)
      write(6,*)'i2dh=0 ne pas passer par sigstepgrid'
      write(6,*)                                                        &
       'dans notebook_vertcoord oter l''option syst hybride'
      stop ' donc stop en attendant! '
      endif
      if(ihybsig==0) return

 
      if(case_==1) then !1111111>
! La subroutine est appelee sigma_levels
!     if(par%rank==0) then
!      i=0 ; j=86
!      do k=kmax+1,1,-1
!        write(6,*)'zz1',k,mask_t(i,j,k),depth_w(i,j,k) 
!      enddo
!      write(6,*)'!!',h_w(i,j),h0_w(i,j)
!      stop 'kiko'
!     endif

       call sigstepgrid_wgrid_case(case_)
       call sigstepgrid_netcdf_w

      endif             !1111111>


      if(case_==2) then !2222222>
! La subroutine est appelee par sigma_levels

       call sigstepgrid_netcdf_r ! load depth_f dsig_f h_f kmin_w

         call obc_mixsigstep(0)        ! obc on kmin_w
         do j=0,jmax+1 ; do i=0,imax+1 ! Maintain consistency between h_w and deph_w(:,:,1)
           h_w(i,j)=-depth_w(i,j,kmin_w(i,j))
         enddo ; enddo
         call obc_h(0) !Border conditions on h_w
         call hz_to_hxyr ! h_u h_v h_f from h_w

! mask_t mask_u mask_v kmin_u kmin_v from kminw:
         call sigstepgrid_kminmask_w2uv

! Compute consistency between depth_w sigma_w dsig_t:
!        call sigma_levels_consistency ! if f grid: dsig, dz, depth
         call sigma_levels_2dv_plot

      endif             !2222222>

      end subroutine sigstepgrid_driver

!................................................................

      subroutine sigstepgrid_wgrid_case(case_)
      use module_principal ; use module_parallele ; use module_global
      implicit none
      integer case_,loop_
      real :: dt_=0.125
#ifdef synopsis
       subroutinetitle='sigstepgrid_wgrid_case'
       subroutinedescription= &
       'Computes depth_w and kmin_w'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE ITERATIVE
! DEBUT:
      if(case_==1) then
!_______________________________________________________________________________

!     if(par%rank==1) then
!     i=2 ; j=165
!     do k=kmax+1,1,-1
!      write(6,*)depth_w(i,j,k)
!     enddo
!     stop 'kiko'
!     endif

      k2=2                                                              !23-09-09

! On stoke l'état initial de depth_w dans anyv3d(:,:,:,0)
! dans l'eventualite d'un rappel vers l'etat initial (si x8/=0)
      do k=1,kmax+1 ; do j=0,jmax+1 ; do i=0,imax+1
       anyv3d(i,j,k,0)=depth_w(i,j,k)
      enddo         ; enddo         ; enddo

      if(nhybsig/=0)stop 'filiere sigma/step lisse H ==> nhybsig=0' !03-02-15
      do 100 loop_=1,nhybsig                                    !02/08/03
      if(par%rank==0.and.mod(loop_,10)==0)write(6,*)'% iter=',100.*real(loop_)/real(nhybsig)

! condition aux limites:
      do j=1,jmax
      do i=1,imax
        depth_w(i,j,0)=2.*depth_w(i,j,1)-depth_w(i,j,2)    !25-01-11
      enddo
      enddo

      x7=1.      ! incoherence hydro max
!     x8=0.01    ! rappel vers grille de reference
      x8=0.0    ! rappel vers grille de reference

      do 10 j=1,jmax
      do 10 i=1,imax
!     if(mask_t(i,j,kmax)==1)then !mskmskmsk>

      do k=1,kmax
      km1=    k-1

! Ces nouvelles lignes supposent que l'on calcule des C.L. en imax+1,0,jmax+1,0 !11/03/09
      jm1=j-1
      jp1=j+1
      im1=i-1
      ip1=i+1

      if( depth_w(i,j,k)>depth_w(i,j,k+1)) then
       write(6,*)'boucle iterative=',loop_
       write(6,*)'i,j,k=',i,j,k
       write(6,*)'depth_w(i,j,k),depth_w(i,j,k+1)' &
                 ,depth_w(i,j,k),depth_w(i,j,k+1)
       stop ' sigstepgrid_wgrid_case erreur croisement de niveaux'
      endif


      x0=(0.5*h_w(i,j)/real(kmax+1))**2
      x10=h_w(i,j)

      x1=                                                           &
         min(                                                       &
         max(                                                       &
         max(                                                       &
         ( depth_w(i,j,k)- depth_w(im1,j,k))                        & ! la pente des z ne doit pas etre
        *( depth_w(i,j,k)- depth_w(im1,j,k)+h_w(i,j)-h_w(im1,j))/x0 & ! plus
       ,                                                            & ! grande que la pente des H
         ( depth_w(i,j,k)- depth_w(im1,j,k+1))                      & ! inconsistence hydrostatique
        *( depth_w(i,j,k)- depth_w(im1,j,km1))/x0                   &
            )                                                       &
            ,zero)                                                  &
            ,dt_)!*mask_t(im1,j,kmax)     !25-01-11


      x2=                                                              &
         min(                                                          &
         max(                                                          &
         max(                                                          &
         ( depth_w(i,j,k)- depth_w(ip1,j,k))                           &
        *( depth_w(i,j,k)- depth_w(ip1,j,k)+h_w(i,j)-h_w(ip1,j))/x0    &
       ,                                                               &
         ( depth_w(i,j,k)- depth_w(ip1,j,k+1))                         &
        *( depth_w(i,j,k)- depth_w(ip1,j,km1))/x0                      &
            )                                                          &
            ,zero)                                                     &
            ,dt_)!*mask_t(ip1,j,kmax)

      x3=                                                              &
         min(                                                          &
         max(                                                          &
         max(                                                          &
         ( depth_w(i,j,k)- depth_w(i,jm1,k))                           &
        *( depth_w(i,j,k)- depth_w(i,jm1,k)+h_w(i,j)-h_w(i,jm1))/x0    &
       ,                                                               &
         ( depth_w(i,j,k)- depth_w(i,jm1,k+1))                         &
        *( depth_w(i,j,k)- depth_w(i,jm1,km1))/x0                      &
            )                                                          &
            ,zero)                                                     &
            ,dt_)!*mask_t(i,jm1,kmax)

      x4=                                                              &
         min(                                                          &
         max(                                                          &
         max(                                                          &
         ( depth_w(i,j,k)- depth_w(i,jp1,k))                           &
        *( depth_w(i,j,k)- depth_w(i,jp1,k)+h_w(i,j)-h_w(i,jp1))/x0    &
       ,                                                               &
         ( depth_w(i,j,k)- depth_w(i,jp1,k+1))                         &
        *( depth_w(i,j,k)- depth_w(i,jp1,km1))/x0                      &
            )                                                          &
            ,zero)                                                     &
            ,dt_)!*mask_t(i,jp1,kmax)



! Les 10 dernieres iterations ne servent qu'à lisser le maillage !06-01-11
!       if(loop_>nhybsig-4) then
!        x1=0.125!*mask_t(i,j,kmax)*mask_t(im1,j  ,kmax)   !25-01-11
!        x2=0.125!*mask_t(i,j,kmax)*mask_t(ip1,j  ,kmax)
!        x3=0.125!*mask_t(i,j,kmax)*mask_t(i  ,jm1,kmax)
!        x4=0.125!*mask_t(i,j,kmax)*mask_t(i  ,jp1,kmax)
!        x8=0.
!       endif

        anyv3d(i,j,k,1)=                                               &
         depth_w(i  ,j  ,k)*(1.-x1-x2-x3-x4)                           &
       + depth_w(im1,j  ,k)*x1                                         &
       + depth_w(ip1,j  ,k)*x2                                         &
       + depth_w(i  ,jm1,k)*x3                                         &
       + depth_w(i  ,jp1,k)*x4

! Ici on effectue un rappel vers la grille initiale (si x8/=0)
        anyv3d(i,j,k,1)=(1.-x8)*anyv3d(i,j,k,1)            &        !25-01-11
                           +x8 *anyv3d(i,j,k,0)


      enddo
!     endif                         !mskmskmsk>
   10 continue

!  12 continue

      do 11 k=1,kmax
      do 11 j=1,jmax
      do 11 i=1,imax
!      if(mask_t(i,j,kmax)==1)depth_w(i,j,k)=anyv3d(i,j,k,1)
       depth_w(i,j,k)=anyv3d(i,j,k,1)
   11 continue

! Eviter le chevauchement en assurant que x1*dz(k+1) < dz(k) < x2*dz(k+1)
! avec comme condition initiale de l'algo: dz(surf)>x6
      x6=0.1 ; x1=0.3   ; x2=3.     !28-01-15
      do j=1,jmax
      do i=1,imax

       k=kmax
       depth_w(i,j,k)= depth_w(i,j,k+1)          &
                -max(  depth_w(i,j,k+1)          &
                     - depth_w(i,j,k  ),x6 ) 

      do k=kmax-1,1,-1                                                 !29-09-09

         depth_w(i,j,k  )=max(min(                                      &
         depth_w(i,j,k  )                                               &
       , depth_w(i,j,k+1)-x1*( depth_w(i,j,k+2)- depth_w(i,j,k+1)))     &
       , depth_w(i,j,k+1)-x2*( depth_w(i,j,k+2)- depth_w(i,j,k+1)))

      enddo

      if(depth_w(i,j,1)>-h_w(i,j)) then !>>>>>
        rap=-h_w(i,j)/depth_w(i,j,1)
        do k=kmax,1,-1
         depth_w(i,j,k)=rap*depth_w(i,j,k)
        enddo
      endif                             !>>>>>

      enddo
      enddo

!................................................!                     !11/03/09
! Conditions aux limites laterales pour PROFWK_Z
      call obc_depth(0)
!................................................!

! bilan d'incoherence hydrostatique:
!     if(mod(loop_,10)==0)call sigstepgrid_hydro_incons

      if(mod(loop_,100)==0)call sigma_levels_2dv_plot

  100 continue ! Fin de boucle iterative


! CALCULER KMIN METHODE VST "FOND":
      call sigstepgrid_bottom
! CALCULER KMIN METHODE VST "SURFACE":
!     call sigstepgrid_surface
! DISTORSION DE LA GRILLE
!     call sigstepgrid_compress

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE ITERATIVE
! FIN.
      endif
!_______________________________________________________________________________

!_______________________________________________________________________________
! CONSTRUCTION DU MELANGE COORDONNEE SIGMA COORDONNEE MARCHES D'ESCALIERS
! PROCEDURE LECTURE FICHIER
      if(case_==2)call sigstepgrid_netcdf_r

 2222 continue                                                       !31/03/08

! Border conditions on kmin_w:
       call obc_mixsigstep(0) ! obc on kmin_w

! Maintain consistency between h_w and deph_w(:,:,1)
      do j=0,jmax+1
      do i=0,imax+1
       h_w(i,j)=-depth_w(i,j,kmin_w(i,j))
      enddo
      enddo
! Border conditions on h_w:
       call obc_h(0)
! extend to other locations of the C-grid
       call hz_to_hxyr

! Premier passage:
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
                mask_t(i,j,k)=1*mask_t(i,j,kmax)
       enddo
      enddo
      enddo
! cas particulier: si il ne reste qu'un niveau on masque tout!
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax)==0)mask_t(i,j,kmax+1)=0
      enddo
      enddo

! OBC mask_c & Pour obtenir les autres points de la grille C:
      call maskt_to_maskuvp                                  !21-10-09

! Compute kmin_u and kmin_v
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_u(i,j,k)==0)kmin_u(i,j)=min(k+1,kmax)
       if(mask_v(i,j,k)==0)kmin_v(i,j)=min(k+1,kmax)
      enddo
      enddo
      enddo

! Deuxieme passage: mise à zero d'office des zones continentales:
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       mask_t(i,j,k)=mask_t(i,j,k)*mask_t(i,j,kmax)
      enddo
      enddo
      enddo

! obc mask_c & obtenir les autres points de la grille C:
      call maskt_to_maskuvp                                  !21-10-09

! Compute consistency between depth_w sigma_w dsig_t:
!     call sigma_levels_consistency 

      end subroutine sigstepgrid_wgrid_case

!..........................................................

      subroutine sigstepgrid_hydro_incons
      use module_principal
      use module_parallele !#mpi
      implicit none

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

      if(par%rank==0) then !000000>
        write(6,*)'Iteration:'
        write(6,*)'hydrostatic inconsistency. Mean value=',sum2/sum1
        write(6,*)'hydrostatic inconsistency.  Max value=',x21
        write(6,*)'and related cell box=',i3,j3,k3
        open(unit=3,file='iteration_grille.out')
         write(3,*)k1
        close(3)
      endif                !000000>


      end subroutine sigstepgrid_hydro_incons

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sigstepgrid_netcdf_w
      use module_principal
      use module_parallele !#mpi
      use pnetcdf
      implicit none
#ifdef synopsis
       subroutinetitle='sigstepgrid_netcdf_w'
       subroutinedescription= &
       'Writes depth_w and kmin_w in a netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      character txt_type_*6
      integer var_dims,var_type,tabdim(5) ! Que represente tabdim?
                                          ! Soit 1,2,3, l'ordre dans lequel
                                          ! sont listees les dimensions du fichiers
                                          ! netcdf. Tabdim indique (via le numero
                                          ! d'ordre) les dimensions dont dépend la
                                          ! variable

! ECRIRE LE FICHIER DE GRILLE:
      txt_type_='double'
!     txt_type_='real'

      do loop_netcdf=0,1

      count_netcdfvar=0

      filval=-9999.
      texte80(3)='none'
      texte80(4)='none'

! Nom du fichier:

      if(loop_netcdf==0) then  !§§§§§§§>
       status=nfmpi_create(par%comm2d,trim(sigstepgridfile),nf_clobber + NF_64BIT_OFFSET, MPI_INFO_NULL,ncid) !10-09-12
      else                     !§§§§§§§>
       status=nfmpi_open(par%comm2d,trim(sigstepgridfile),nf_write + NF_64BIT_OFFSET, MPI_INFO_NULL,ncid) !10-09-12
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop module_offline erreur 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! write variables dimensions:
      call graph_out_trueaxis_nfmpi

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lon_t(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lon_t(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='longitude_t'
      texte80(2)='degrees_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lon_f(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lon_f(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='longitude_f'
      texte80(2)='degrees_east'                              ! units
      texte80(3:4)='longitude_at_f_location'                 ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      call netcdf_main('_f')

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lat_t(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lat_t(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='latitude_t'
      texte80(2)='degrees_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lat_f(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lat_f(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='latitude_f'
      texte80(2)='degrees_north'                              ! units
      texte80(3:4)='latitude_at_f_location'                   ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      call netcdf_main('_f')


      if(fgrid_or_wgrid==fgrid_case) then !*********>

       if(allocated(depth_f)) then !fffff>
        if(loop_netcdf==1) then !--------->
         do k=1,kmax+1 ; do j=1,jmax+1 ; do i=1,imax+1
           anyv3d(i,j,k,1)=depth_f(i,j,k)
         enddo         ; enddo         ; enddo
        endif                   !--------->
        texte80(1)='depth_f'  ; texte80(2)='m' 
        texte80(3:4)='depth_at_f_location'   
        texte80(5)='ZYX' ; texte80(7)='double'
        call netcdf_main('_f')
       endif                       !fffff>

!      if(allocated(dsig_f))  then !fffff>
!       if(loop_netcdf==1) then !--------->
!        do k=1,kmax+1 ; do j=1,jmax+1 ; do i=1,imax+1
!          anyv3d(i,j,k,1)=dsig_f(i,j,k)
!        enddo         ; enddo         ; enddo
!       endif                   !--------->
!       texte80(1)='dsig_f'  ; texte80(2)='none' 
!       texte80(3:4)='delta_sigma_at_f_location'   
!       texte80(5)='ZYX' ; texte80(7)='double'
!       call netcdf_main('_f')
!      endif                       !fffff>

       if(allocated(h_f))     then !fffff>
        if(loop_netcdf==1) then !--------->
         do j=1,jmax+1 ; do i=1,imax+1
           anyv3d(i,j,1,1)=h_f(i,j)
         enddo         ; enddo   
        endif                   !--------->
        texte80(1)='h_f'  ; texte80(2)='m' 
        texte80(3:4)='sea_floor_depth_at_f_location'   
        texte80(5)='YX' ; texte80(7)='double'
        call netcdf_main('_f')
       endif                       !fffff>

      endif                               !*********>


      if(fgrid_or_wgrid==wgrid_case) then !WWWWWWW>

      if(loop_netcdf==1) then !--------->
      do k=1,kmax+1
      do j=0,jmax+1
      do i=0,imax+1
!         anyvar3d(i,j,k)=depth_w(i,j,k)
          anyv3d(i,j,k,1)=depth_w(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_w'
      texte80(2)='m'                             ! units
      texte80(3:4)='depth_at_w_location'    ! long_name
      texte80(5)='ZYX' ; texte80(7)='double' !01-05-14
      call netcdf_main('_w')

!     if(allocated(dsig_t)) then !ooooo> !28-11-14
!     if(loop_netcdf==1) then !--------->
!     do k=1,kmax
!     do j=0,jmax+1
!     do i=0,imax+1
!          anyv3d(i,j,k,1)=dsig_t(i,j,k) ! ecrire en double
!     enddo
!     enddo
!     enddo
!     endif                  !--------->
!     texte80(1)='dsig_t'
!     texte80(2)='none'                             ! units
!     texte80(3:4)='delta_sigma_at_t_location'   ! long_name
!     texte80(5)='ZYX' ; texte80(7)='double' !01-05-14
!     call netcdf_main('_t')
!     endif                       !ooooo> !28-11-14

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
            anyvar2d(i,j)=h_w(i,j) ! ecrire en real
          anyv3d(i,j,1,1)=h_w(i,j) ! ecrire en double
       else
            anyvar2d(i,j)=-9999.
          anyv3d(i,j,1,1)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='hm_w' ; texte80(2)='m'                     ! variable ; units
      texte80(3:4)='model_sea_floor_depth_below_geoid'
      texte80(5)='YX' ; texte80(7)='double' !01-05-14
      call netcdf_main('_w')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
            anyvar2d(i,j)=h_w(i,j) ! ecrire en real
          anyv3d(i,j,1,1)=h_w(i,j) ! ecrire en double
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_w' ; texte80(2)='m'                     ! variable ; units
      texte80(3:4)='sea_floor_depth_at_w_location'
      texte80(5)='YX' ; texte80(7)='double' !01-05-14
      call netcdf_main('_w')

      endif                               !WWWWWWW>

      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=kmin_w(i,j)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='kmin_w'   ; texte80(2)='none'      
      texte80(3)='first_level_from_bottom'    
      texte80(4)='first_level_from_bottom'  
      texte80(5)='YX' ; texte80(7)='short' 
      call netcdf_main('_w')

!     var_validmin=1.  ; var_validmax=real(kmax) !03-12-14
!     var_addoffset=0. ; var_scalefactor=1.
!     if(loop_netcdf==1) then !--------->
!     do j=0,jmax+1
!     do i=0,imax+1
!         anyvar2d(i,j)=ksurf_t(i,j)
!     enddo
!     enddo
!     endif                  !--------->
!     texte80(1)='ksurf_t'   ; texte80(2)='none'      
!     texte80(3)='first_level_under_surface'    
!     texte80(4)='first_level_under_surface'    
!     texte80(5)='YX' ; texte80(7)='short' 
!     call netcdf_main('_t')

! Ecrire le masque !21-07-15
      var_validmin=0.  ; var_validmax=1.
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=mask_t(i,j,kmax)
       enddo ; enddo
      endif                   !--------->
      texte80(1)='mask2d_t'   ; texte80(2)='none'      
      texte80(3)='sea_land_mask' ; texte80(4)='sea_land_mask'    
      texte80(5)='YX' ; texte80(7)='short' 
      call netcdf_main('_w')

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

! Definition of variables: done.
      call netcdf_general_attributes(ncid) !03-04-12
      status=nfmpi_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      status=nfmpi_close(ncid)


!#ifdef parallele
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!      call barriere(looproc_   ,3)
!#endif

      enddo  ! fin de boucle sur loop_netcdf
      if(par%rank==0) then !>>>>>
        write(6,'(a)')'sigstepgridfile created:'
        write(6,'(a)')trim(sigstepgridfile) !09-01-15
      endif                !>>>>>

      end subroutine sigstepgrid_netcdf_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sigstepgrid_netcdf_r
      use module_principal
      use module_parallele !#mpi
      use pnetcdf
      implicit none
      integer(kind=MPI_OFFSET_KIND) start(3)
      integer(kind=MPI_OFFSET_KIND) edge(3)
      integer(kind=MPI_OFFSET_KIND) lendim_
      integer ncid_
#ifdef synopsis
       subroutinetitle='sigstepgrid_netcdf_r'
       subroutinedescription= &
       'Reads depth_w and kmin_w in a netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      status=nfmpi_open(par%comm2d,trim(sigstepgridfile),nf_nowrite, MPI_INFO_NULL,ncid_)
      if(status/=0) then
       write(6,'(a,a)')'nf_open on ',trim(sigstepgridfile)
       stop ' stop nfmpi_open sigstepgrid_netcdf_r'
      endif

!...........
! Verification des dimensions du fichiers !24-05-15
! definir/verifier varstart:
       status=nfmpi_inq_dimid(ncid_,'ni_t',k0)       ;if(status/=0)stop 'Err1890'
       status=nfmpi_inq_dimlen(ncid_,k0,lendim_)     ;if(status/=0)stop 'Err1891'
       if(lendim_/=iglb+2)stop 'Err dim1 sigstepgridfile'

       status=nfmpi_inq_dimid(ncid_,'nj_t',k0)       ;if(status/=0)stop 'Err1892'
       status=nfmpi_inq_dimlen(ncid_,k0,lendim_)          ;if(status/=0)stop 'Err1893'
       if(lendim_/=jglb+2)stop 'Err dim2 sigstepgridfile'

       status=nfmpi_inq_dimid(ncid_,'nk_t',k0)       ;if(status/=0)stop 'Err1894'
       status=nfmpi_inq_dimlen(ncid_,k0,lendim_)          ;if(status/=0)stop 'Err1895'
       if(lendim_/=kmax)  stop 'Err dim3 sigstepgridfile'

!............
! kmin_w
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i
      status=nfmpi_inq_varid(ncid_,'kmin_w',var_id) 
      if(status/=0)stop ' stop trouve pas kmin_w'

       status=nfmpi_get_vara_int_all(ncid_,var_id,start(1:2) &
                                                  ,edge(1:2) &
                                  ,kmin_w(0:imax+1,0:jmax+1))
      if(status/=0)stop ' stop get kmin_w'

!............
! ksurf_t
!     start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
!     start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i
!     status=nfmpi_inq_varid(ncid_,'ksurf_t',var_id) 
!     if(status/=0)stop ' stop trouve pas ksurf_t'

!      status=nfmpi_get_vara_int_all(ncid_,var_id,start(1:2) &
!                                                 ,edge(1:2) &
!                                ,ksurf_t(0:imax+1,0:jmax+1))
!     if(status/=0)stop ' stop get ksurf_t'

!     do j=0,jmax+1 ; do i=1,imax+1
!      ksurf_u(i,j)=min(ksurf_t(i,j),ksurf_t(i-1,j))
!     enddo ; enddo
!     do j=1,jmax+1 ; do i=0,imax+1
!      ksurf_v(i,j)=min(ksurf_t(i,j),ksurf_t(i,j-1))
!     enddo ; enddo

!............
! depth_w
      start(3)=1                 ; edge(3)=kmax+1
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i
      status=nfmpi_inq_varid(ncid_,'depth_w',var_id) 
      if(status/=0)stop ' stop trouve pas depth_w'

       status=nfmpi_get_vara_double_all(ncid_,var_id,start(1:3) &
                                                     ,edge(1:3) &
                           ,depth_w(0:imax+1,0:jmax+1,1:kmax+1))
      if(status/=0)stop ' stop get depth_w'

      status=nfmpi_close(ncid_)


      end subroutine sigstepgrid_netcdf_r

!........................................................

      subroutine sigstepgrid_kminmask_w2uv
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='sigstepgrid_kminmask_w2uv'
       subroutinedescription= &
       'computes mask_t _u _v and kmin_u _v from kmin_w'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Compute mask_t from kmin_w
         do j=0,jmax+1 ; do i=0,imax+1
          do k=1,kmin_w(i,j)-1
                mask_t(i,j,k)=0
          enddo
          do k=kmin_w(i,j),kmax+1
                mask_t(i,j,k)=1*mask_t(i,j,kmax)
          enddo
         enddo ; enddo
! Compute Land/Sea mask on other grid nodes
         call maskt_to_maskuvp                                  !21-10-09
! Compute kmin_u and kmin_v
         do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
           if(mask_u(i,j,k)==0)kmin_u(i,j)=min(k+1,kmax)
           if(mask_v(i,j,k)==0)kmin_v(i,j)=min(k+1,kmax)
         enddo       ; enddo         ; enddo
! Compute consistency between depth_w sigma_w dsig_t:
!        call sigma_levels_consistency 

      end subroutine sigstepgrid_kminmask_w2uv

!........................................................

      subroutine sigstepgrid_bottom
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='sigstepgrid_bottom'
       subroutinedescription= &
       'Removes levels from bottom'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      h_w(:,:)=h0_w(:,:)
      call sigma_levels_2dv_plot

! Determiner kmin_w
       do j=0,jmax+1 ; do i=0,imax+1
        k2=1             !30-01-11
        do k=1,kmax
         if(depth_w(i,j,k)<=-h_w(i,j).and.depth_w(i,j,k+1)>-h_w(i,j)) then
          k2=k
         endif
!        kmin_w(i,j)=min(k2,kmax-5) ! au 6 niveaux...
         kmin_w(i,j)=k2
        enddo
       enddo  ; enddo

! Dans les zones où la grille VST n'est pas appliquee (kmax+1-kmin_w(i,j)<nbvstepmin)
! il faut appliquer un critere rmax a la grille sigma d'ou l'appel a lissebathy case 2:
       call lissebathy(2,0.,0,0,0,0,0) !25-05-15

!..............................................
! The lowest level is adjusted to the bathymetry:                      !19-10-09
! Change kmin_w if the number of vertical levels < k0 and then adjust depth_w
      k0=min(kmax,nbvstepmin) ! k0 = minimum number of vertical levels !17-03-15
      do j=0,jmax+1
      do i=0,imax+1

! Change kmin_w if the number of vertical levels < k0 and then adjust depth_w:
        if(kmax+1-kmin_w(i,j)<k0) then !>>>>>>

           kmin_w(i,j)=kmax+1-k0

! Previous sigma distribution starting from kmin_w+1 before depth_w adjustement
! attention c'est bien kmin+1 et non pas kmin pour une meilleure continuite horizontale
          do k=kmin_w(i,j),kmax+1
           anyv1d(k,1)=(depth_w(i,j,k     )-depth_w(i,j,kmin_w(i,j)+1)) &
                      /(depth_w(i,j,kmax+1)-depth_w(i,j,kmin_w(i,j)+1))  
          enddo
! Compute depth_w using the same sigma distribution but now with true h
          do k=kmin_w(i,j),kmax+1
           depth_w(i,j,k)=anyv1d(k,1)*(depth_w(i,j,kmax+1)+h_w(i,j)) &
                                                          -h_w(i,j)
          enddo
     
        endif                          !>>>>>>

! The lowest level is adjusted to the bathymetry:                      !19-10-09
        depth_w(i,j,kmin_w(i,j))=-h_w(i,j)
        if(depth_w(i,j,kmin_w(i,j))==depth_w(i,j,kmin_w(i,j)+1))then !warning>
           depth_w(i,j,kmin_w(i,j)+1)=0.5*(depth_w(i,j,kmin_w(i,j)  )  &
                                          +depth_w(i,j,kmin_w(i,j)+2))  
        endif                                                        !warning>
      enddo
      enddo

!..............................................
! La condition limite au fond = niveaux "quasi collés"
! and Special case of uncovered areas
!..............................................
      small3=1.d-5
      do j=0,jmax+1
      do i=0,imax+1
! Under the sea bottom level:
       do k=0,kmin_w(i,j)-1
        depth_w(i,j,k)= depth_w(i,j,kmin_w(i,j))+(k-kmin_w(i,j))*small3
       enddo
! Special case of uncovered areas
       do k=kmin_w(i,j),kmax+1
        if(depth_w(i,j,k)<-h_w(i,j))depth_w(i,j,k)=-h_w(i,j)+(k-kmin_w(i,j))*small3
       enddo
      enddo
      enddo

      end subroutine sigstepgrid_bottom

!.......................................................................
#ifdef bidon
      subroutine sigstepgrid_surface
      use module_principal ; use module_parallele
      implicit none
      integer kmax_
      real kmaxdec_
#ifdef synopsis
       subroutinetitle='sigstepgrid_bottom'
       subroutinedescription= &
       'Removes levels from surface'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      h_w(:,:)=h0_w(:,:)
      call sigma_levels_2dv_plot

! Determiner le niveau vertical decimal qui correpond a la
! bathymetry:
       ksurf_t=kmax
       do j=0,jmax+1 ; do i=0,imax+1
        deck=1

        do k=1,kmax
         if(depth_w(i,j,k)<=-h_w(i,j).and.depth_w(i,j,k+1)>-h_w(i,j)) then
          deck=k+(-h_w(i,j)-depth_w(i,j,k))/(depth_w(i,j,k+1)-depth_w(i,j,k))
         endif
        enddo

        kmaxdec_=kmax-(deck-1.)
        kmax_=min(int(kmaxdec_)+1,kmax)
        ksurf_t(i,j)=kmax_

! Essai d'un decoupage sigma simple:
        x1=h_w(i,j)/kmaxdec_ ! Resolution "decimale"
        depth_w(i,j,1)=-h_w(i,j)
        do k=2,kmax_
         depth_w(i,j,k)=depth_w(i,j,k-1)+x1

!---debug---->
         if(depth_w(i,j,k)>0.) then
           write(6,*)'Err z>0 i,j,k z',i,j,k,kmaxdec_
           stop 'Err 1348'
         endif
!---debug---->

        enddo
        do k=kmax_+1,kmax+1
         depth_w(i,j,k)=max(0.,-h_w(i,j))
        enddo

! Eviter que la derniere couche active ne soit trop fine:
            depth_w(i,j,kmax_)=                        &
        min(depth_w(i,j,kmax_),                        &
           (depth_w(i,j,kmax_+1)+0.5*depth_w(i,j,kmax_-1))/1.5)

        kmin_w=1

!       if(par%rank==0.and.i==imax/2.and.j==4) then
!        write(6,*)'kmaxdec_',kmaxdec_
!        write(6,*)'h et deck',-h_w(i,j),deck
!        do k=kmax,1,-1
!         write(6,*)k,depth_w(i,j,k)
!        enddo 
!        stop 'coco'
!       endif

       enddo  ; enddo
      end subroutine sigstepgrid_surface
#endif
!...................................................................
      subroutine sigstepgrid_compress
      use module_principal ; use module_parallele
      implicit none
      integer kmax_
      real kmaxdec_
#ifdef synopsis
       subroutinetitle='sigstepgrid_bottom'
       subroutinedescription= &
       'Removes levels from surface'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      

      do j=0,jmax+1 ; do i=0,imax+1

      if(par%rank==0.and.i==imax/2.and.j==20) then
       write(6,*)'h0 h',h0_w(i,j),h_w(i,j)
      endif

       sum1=0.
       do k=1,kmax
        rap=(real(kmax-k)/real(kmax-1))
        if(k<kmax/2) then
         rap=1
        else
         rap=0
        endif
        dsig_t(i,j,k)=rap*h_w(i,j)+(1.-rap)*h0_w(i,j)
        sum1=sum1+dsig_t(i,j,k)
       enddo

       sum2=0.
       do k=1,kmax
        dsig_t(i,j,k)=dsig_t(i,j,k)/sum1
        sum2=sum2+dsig_t(i,j,k)
        if(par%rank==0.and.i==imax/2.and.j==10)write(6,*)k,dsig_t(i,j,k)
       enddo
       if(par%rank==0.and.i==imax/2.and.j==10)write(6,*)'sum2=',sum2

       depth_w(i,j,1)=-h0_w(i,j)
       do k=1,kmax
                   depth_w(i,j,k+1)=    &
                   depth_w(i,j,k  )     &
         +h0_w(i,j)*dsig_t(i,j,k)
       enddo

      enddo         ; enddo

      h_w(:,:)=h0_w(:,:)
      call sigma_levels_2dv_plot
!     stop 'toto'

      end subroutine sigstepgrid_compress
