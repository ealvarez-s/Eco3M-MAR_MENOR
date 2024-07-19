      subroutine maskt_to_maskuvp
!______________________________________________________________________
! S model
! release S.26 - last update: 21-05-17
!______________________________________________________________________

      use module_principal
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='maskt_to_maskuvp'
       subroutinedescription= &
       'Sea-land mask at u v f points deduced from t points'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! modifs: 21/02/02: extension du calcul à toute la grille. Inutile pour
!                   le model_ mais requis par le model_ inverse VICOM
!         14/03/02: verification rigoureuse avec Francis
!         05/06/02: amenagements pour model_ inverse de Francis
!         05/02/03: amenagement pour filiere hybride sigma_z:
!                   le schema original est conservé pour ne pas avoir
!                   à faire une modif suplementaire sur bathycote et cie
!                   Par contre on ajoute un nouveau choix (choix=2) qui
!                   execute le masquage en boucle 3D. En pratique le 3eme
!                   argument n'est pas codé en dur à NR mais fait
!                   l'objet d'une boucle que l'on deroule ou qu'on ne
!                   deroule pas selon la valeur de ichoix.
!         16/10/08: La C.L. sur   mask_x &   mask_y est maintenant dans routine
!                   obc_morpho.F
!         06/03/09: suite à obc sur   mask_c, les boucles de calcul sont aggrandies.
!                   L'idee est que si le calcul de   mask_c est etendu spatialement
!                   on se peut se passer de C.L. pour   mask_x et   mask_y. Les
!                   C.L. repose alors uniquement sur   mask_c et la parallelisation
!                   s'en trouve clarifiee
!         02/04/09: modif algo   mask_x   mask_y   mask_c
!         21-10-09: - routine renommé (ex z_to_xyr) maskz_to_maskxyr
!                   - Pour eviter bug de parallelisation, la condition aux limites
!                     sur mask_c est appliquée dans la presente routine
!         22-10-09: suite point precedent (debug argument)
!         21-03-14  Pas de C.L. si le masque est lu dans la grille nemo (cas nemo offline)
!         16-07-14  echange mpi indispensable pour cas nemo
!         21-05-17 call obc_mask(2) ! commentE le !21-05-17
!...............................................................................

      if(initialgridfile_txt=='none') then !---------> !21-03-14

! Modifier le masque pour permettre la compatibilite avec
! les conditions aux limites radiatives
!      call obc_mask(2) ! commentE le !21-05-17

! C.L.
       call obc_mask(1) ! C.L. sur rang 1 ( 0, imax+1,  0, jmax+1)
       call obc_mask(3) ! C.L. sur rang 2 (-1, imax+2, -1, jmax+2)

! MPI:
       call obc_mask_mpi ! echange zc   !16-07-14

      else                                 !---------> !21-03-14

! Dans le cas de la grille NEMO faire seulement la C.L. rang 2 (pour schema advection_bioup2up3)
       call obc_mask(3) ! C.L. sur rang 2 (-1, imax+2, -1, jmax+2)

! MPI:
       call obc_mask_mpi ! echange zc   !16-07-14

      endif                                !---------> !21-03-14

      do 1000 k=1,kmax+1                                                  !05/02/03

!   mask_x:
      do j=0,jmax+1
      do i=0,imax+2
        mask_u(i,j,k)=  mask_t(i,j,k)*  mask_t(i-1,j,k)                !02/04/09
      enddo
      enddo
!   mask_y:
      do j=0,jmax+2
      do i=0,imax+1
        mask_v(i,j,k)=  mask_t(i,j,k)*  mask_t(i,j-1,k)                !02/04/09
      enddo
      enddo
!   mask_r:
      do j=0,jmax+1
      do i=0,imax+1
        mask_f(i,j,k)=  mask_t(i  ,j,k)*  mask_t(i  ,j-1,k)            & !02/04/09
                     *  mask_t(i-1,j,k)*  mask_t(i-1,j-1,k)            !02/04/09
      enddo
      enddo

 1000 continue                                                         !05/02/03


      end subroutine maskt_to_maskuvp
