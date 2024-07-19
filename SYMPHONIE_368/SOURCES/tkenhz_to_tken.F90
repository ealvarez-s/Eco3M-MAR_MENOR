      subroutine tkenhz_to_tken(ki2)
!______________________________________________________________________
!
! S model
! release 2010.10  - last update: 24-06-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
      integer ki2
#ifdef synopsis
       subroutinetitle='tkenhz_to_tken'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! OUTPUT: TKEN_Z

! Calcul de TKEN_Z a partir de TKENHZ_Z
! modifs 19/07/01: passage à la coordonnée sigma généralisée
! modifs 01/11/01: inversion de l'ordre des boucles: J passe avant I
!        02/07/06: schema forward donc TKENHZ n'est plus utilisé

      return                                                           !02/07/06

!     DO 100 K=1,NR
!     DO 101 J=1,NECO+1
!     DO 101 I=1,MECO+1
!     TKEN_Z(I,J,K)=TKENHZ_Z(I,J,K,KI2)/
!    1 (0.5*HZ_Z(I,J,KI2)*(DSIG_Z(I,J,K)+DSIG_Z(I,J,K-1)))             !19/07/01
! 101 CONTINUE
! 100 CONTINUE

      return
      end
