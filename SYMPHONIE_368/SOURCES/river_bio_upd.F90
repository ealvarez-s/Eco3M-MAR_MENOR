      subroutine river_bio_upd
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 18-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='river_bio_upd'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! modifs: 22/05/06: correction bug sur K1 (KR au lieu de KR)
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes


!----------------------------------------------------------------------
!     BUT:
!     ----
!     Ce programme lit les concentrations des fleuves dan un fichier
!     et met a jour à chaque pas de temps
!
!
!     OUTPUT:
!     ------
!     RIVER_BIO
!
!
!     INPUT:
!     ------
!     ILD_TO_SD  ! liquid discharge to solid discharge
!
!
!     MODIFS:
!     -------
!
!-----------------------------------------------------------------------
! mise en service: 05/04/06
!    mises a jour: 06/04/06 argument ICHOIX inutile
!-----------------------------------------------------------------------


      if(nriver.eq.0)return ! Pas de riviere?  alors on quitte la routine.
      if(vbmax.eq.0) return ! Pas de biologie? alors on quitte la routine.


! LE CAS PARTICULIER DE L'ETAT INITIAL
! 1. Definir les parametres d'avancee dans le temps (date repere & echantillonage)
! Debut:
      if(iteration3d==kount0) then !.................!       !18-04-11

        do kr=1,nriver ! debut boucle sur les rivieres
         i1=bioriv_date(1,kr)
         i2=bioriv_date(2,kr)
         i3=bioriv_date(3,kr)
         i4=bioriv_date(4,kr)
         i5=bioriv_date(5,kr)
         i6=bioriv_date(6,kr)
         write(6,*)'appel à date_to_kount dans river_bio_upd'
!        call date_to_kount
         call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!        bioriv_dt(kr,2)=xdtk_out                     ! date repere NC=1 en kount
!        bioriv_dt(kr,1)=bioriv_info(kr)*3600./dti_fw ! duree, en kount, entre 2 infos
         bioriv_dt(kr,2)=elapsedtime_out              !18-04-11
         bioriv_dt(kr,1)=bioriv_info(kr)*3600.        !18-04-11

        enddo          ! fin boucle sur les rivieres

      endif                    !.................!
! Fin.



      do 1000 kr=1,nriver ! boucle sur les rivieres
      if(ild_to_sd(kr).eq.1) then ! ILD ILD ILD ILD >                  !22/05/06


      ksecu=0 ! A priori on n'ouvre pas le fichier de "donnee riviere"

! Est on de part et d'autre d'une nouvelle echeance?
      x2=(elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1)     !18-04-11
      x1=(elapsedtime_bef-bioriv_dt(kr,2))/bioriv_dt(kr,1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
! Si la reponse est oui alors il faut lire le fichier:
      if(j2-j1.eq.1)         ksecu=1  ! indique que nous lirons le fichier
      if(iteration3d==kount0)ksecu=1  ! L'initialisation impose de lire le fichier !18-04-11


      if(ksecu.eq.1) then !*************************************>

      if(iteration3d==kount0)then !---------------->        !18-04-11
          k4=0 ! à l'initialisation on lit 2 champs
      else                        !---------------->                       !27/04/04
! sinon l'ancienne echeance est sauvée avant de lire la suivante:
          k4=2
          do vb=1,vbmax
           river_bio(vb,kr,0)=river_bio(vb,kr,2)
          enddo
      endif                       !---------------->



! Ouverture du fichier:
      open(unit=3,file=bioriv_file(kr)                                  &
                 ,access='direct'                                       &
                 ,recl=4*vbmax                                          &
                 ,form='unformatted')


! lecture du fichier:
      do k3=k4,2,2
         x2=(elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1)      !18-04-11
         nc=int(1.+x2)+k3/2
         read(3,rec=nc)(river_bio(vb,kr,k3),vb=1,vbmax)
      enddo

      close(3)
      endif               !*************************************>



! Interpolation temporelle entre 2 echeances du fichier:
      rap=      (elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1)  & !18-04-11
          -int( (elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1) )

      do vb=1,vbmax
       river_bio(vb,kr,1)=(1.-rap)*river_bio(vb,kr,0)                  &
                         +    rap *river_bio(vb,kr,2)
      enddo


      endif                       ! ILD ILD ILD ILD >
 1000 continue ! fin de boucle sur les rivieres

      return
      end
