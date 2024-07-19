      subroutine sediment_bio
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
#ifdef synopsis
       subroutinetitle='sediment_bio'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! OUTPUT: FLUXBIO_Z
!...............................................................................
! Version Date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         25/10/02: Fluxbio_c homogene à grandeur primaire, n'est plus
!                   a etre multiplie par DTI
!         17/03/03: date de debut et de fin de calcul pour traceur
! 2010.8  03-05-10  Terminologie
!...............................................................................


      do 33 vb=1,vbmax
!     if(      kount.ge.kount_bio(vb,1)
!    &    .and.kount.le.kount_bio(vb,2)) then !§§§§§§§§§§§§§§§§§§§>    !19/03/03

      do 34 j=1,jmax
      do 34 i=1,imax
! fluxbio_c homogene à:
!cc   FLUXBIO_Z(I,J,KB)=    WSED(KB)*BIO_Z(I,J,1,KB)                   !25/10/02

! convention de signe: si les 2 flux sont positifs c'est une
! source pour la colonne d'eau si ils sont negatifs c'est un
! puit.
! flux erosion remise en suspension
! Attention uniquement pour la remise en suspension!
! la sedimentation se fait autrement.
      fluxbio_w(i,j,vb,1)=0.

! flux à l'interface air/mer
      fluxbio_w(i,j,vb,2)=0.
!c    IF(KB.EQ.1)FLUXBIO_Z(I,J,KB,2)= 1.E-3                            !25/10/02
!c    IF(KB.EQ.2)FLUXBIO_Z(I,J,KB,2)= 0.

   34 continue

!     endif                                   !§§§§§§§§§§§§§§§§§§§>    !19/03/03

   33 continue

      end subroutine sediment_bio

! amelioration:
! tout d'abord ajouter une 4eme dimension à fluxbio_c:
! fluxbio_c(i,j, q ,KB)
! q=1 remise en suspension
! q=0 deposition implicite recalculee apres coup pour pouvoir faire
! des bilans: fluxbio_c(i,j,0,KB)=DTI*WSED(KB)*BIOHZ(2)/HZ_Z(2) où 2 fait
! reference au temps t+1
! q=2 echange interface ocean atmosphere: penser à introduire ce cas en surface
