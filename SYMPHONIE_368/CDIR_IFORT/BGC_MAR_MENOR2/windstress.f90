










      subroutine windstress
!______________________________________________________________________
! SYMPHONIE ocean model
! release 290 - last update: 23-10-20
!______________________________________________________________________

      use module_principal
      implicit none

!..............................................................................
! modifs 01/11/01: inversion de l'ordre des boucles: J passe avant I
!        21/03/02: par defaut le vent est nul. Ca evite les erreurs
!                  d'etourderie comme par ex lancer model_streamf en
!                  pensant (à tort) que le vent est nul.
!        08-10-14  possibilite de construire des flux a partir des formules bulk
!        25-01-15  Attention par defaut wstress_ =0 pour prévoir le cas des
!                  simulation forcee par les vagues seule (sans modele meteo)
!                  car l'effet taw-two est un delta qui s'ajoute a wstress_...
!        05-08-16  ajout d'un commentaire
!        23-11-15  3 arguments passent dans bulk_formulae
!        24-08-16  deplacement de l'etiquette bidon pour permettre calcul wstress_w
!..............................................................................

!...........................................................................
!        25-01-15  Attention par defaut wstress_ =0 pour prévoir le cas des
!                  simulation forcee par les vagues seule (sans modele meteo)
!                  car l'effet taw-two est un delta qui s'ajoute a wstress_...
!                  wstress=rho*(ustar**2)
!                  wstress_u(:,:,1)=0. ; wstress_v(:,:,1)=0. ! reset necessaire a module_wave
! v290   23-10-20  La contrainte enoncee ci-dessus disparait avec la multiplication par 
!                  iairsea, desormais, de wstress_u dans module_wave, ce
!                  qui evite ce reset toutes les iterations.
!...............................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!......................................................................! m[°o°]m

      RETURN !23-10-20




! Module de la tension de vent:

      do j=1,jmax
      do i=1,imax

      wstress_w(i,j)=sqrt(                                              &
       ((wstress_u(i,j,1)+wstress_u(i+1,j,1))/2.)**2+                   &
       ((wstress_v(i,j,1)+wstress_v(i,j+1,1))/2.)**2)

      enddo
      enddo

!     CALL GRAPH_OUT
!     WRITE(6,*)'DANS SBR WINSTRESS'
!     WRITE(6,*)'-----------------------'

      end subroutine windstress
