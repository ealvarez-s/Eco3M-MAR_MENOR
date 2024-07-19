      subroutine tridiagonalsolver(case_,i1_,i2_,j1_,j2_,k2_)
!______________________________________________________________________
! S model
! release S.26 - last update: 25-06-14
!______________________________________________________________________

      use module_principal
      implicit none
      integer case_,i1_,i2_,j1_,j2_,k2_
!...............................................................................
! Version Date     Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         09/07/03: le systeme ECT est résolu en K=1 et K=NR
!         02/12/04: tableaux gama et beta remplacés par anyvar3D
! 2010.9  06-06-10  tridiagsyst_dbp renommé tridiagonalsolver
! S.26    25-06-14  algo revu (1 division en moins par rapport au precedent)
!                   Le second traceur (S, bio(i,j,k,2:)) ne donne pas
!                   pas lieu à un recalcul de la bi-diagonalisation de
!                   la matrice principale. Attention le tableau tridia_in
!                   a de nouvelles dimension (4eme dimension commence à 0)
!                   Par securite (temporairement) une subroutine
!                   tridiagonalsolver_cb verifie verifie les dimensions
!...............................................................................
#ifdef synopsis
       subroutinetitle='tridiagonalsolver'
       subroutinedescription= &
      'Transforms the 3-diagonal matrix of the vertical mixing implicit'&
      //' scheme into a 2-diagonal matrix and computes the solution.'

       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Résolution du systeme d'equations lineaires tridiagonal

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

! case_=0 correspond au cas particulier ou la transformation de la
! matrice principale tridiagonale en une matrice bidiagonale (archivee
! dans le tableau tridia_in(i,j,k,0) a déjà été faite lors d'une étape
! precedente. A priori cela concerne S, les traceurs bio à compter de
! la deuxieme variable.

! Check bounds:
      if(iteration3d==iteration3d_restart)call tridiagonalsolver_cb

      k=1
      if(case_/=0) then !------->
       do j=j1_,j2_ ; do i=i1_,i2_
        tridia_in(i,j,k,0)=1./tridia_in(i,j,k,2)
       enddo ; enddo
      endif             !------->
      do j=j1_,j2_ ; do i=i1_,i2_
       tridia_in(i,j,k,-1)=tridia_in(i,j,k,4)               &
                          *tridia_in(i,j,k,0)
      enddo ; enddo


      if(case_/=0) then !------->
       do k=2,k2_ ; do j=j1_,j2_  ; do i=i1_,i2_

        tridia_in(i,j,k,0)=1./(tridia_in(i,j,k  ,2)      &
                              -tridia_in(i,j,k  ,1)      &
                              *tridia_in(i,j,k-1,3)      &
                              *tridia_in(i,j,k-1,0))

       enddo ; enddo ; enddo
      endif             !------->

      do k=2,k2_ ; do j=j1_,j2_  ; do i=i1_,i2_

        tridia_in(i,j,k,-1)=( tridia_in(i,j,k  ,4)           &
                             -tridia_in(i,j,k  ,1)           &
                             *tridia_in(i,j,k-1,-1) )        &
                             *tridia_in(i,j,k  ,0)

      enddo ; enddo ; enddo


      do j=j1_,j2_ ; do i=i1_,i2_
       tridia_out(i,j,k2_)=tridia_in(i,j,k2_,-1)
      enddo ; enddo

      do k=k2_-1,1,-1
      do j=j1_,j2_
      do i=i1_,i2_

        tridia_out(i,j,k)= tridia_in(i,j,k,-1)                       &
                          -tridia_in(i,j,k,3)                        &
                         *tridia_out(i,j,k+1)                        &
                          *tridia_in(i,j,k,0)
      enddo
      enddo
      enddo

      end subroutine tridiagonalsolver

!...........................................................

      subroutine tridiagonalsolver_cb
      use module_principal
      implicit none
      integer,dimension(:),allocatable :: dim4d_
#ifdef synopsis
       subroutinetitle='tridiagonalsolver_cb'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Check bounds

      allocate(dim4d_(4))

      dim4d_=lbound(tridia_in)

      if(dim4d_(4)/=-1) then !-------->

       write(6,*)' Allocation error on tridia_in:'  &
                ,' 4th lower bound must 0.'
       stop ' Stop in subroutine tridiagonalsolver_cb'

      endif                  !-------->

      deallocate(dim4d_)

      end subroutine tridiagonalsolver_cb
