      subroutine atm_bio_upd
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 18-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http:// 
!
!______________________________________________________________________

      use ModuleComBlock
      use UserDeclaration
      use ModuleDeclaration
      use module_principal
      implicit none
 

!----------------------------------------------------------------------
!     BUT:
!     ----
!     Ce programme lit les pCO2 mesurees dans l'atmosphere dans un fichier
!     et met a jour à chaque pas de temps
!
!
!     OUTPUT:
!     ------
!     pCO2Air
!
!
!     MODIFS:
!     -------
!
!-----------------------------------------------------------------------
! mise en service: 05/04/06
!    mises a jour: 06/04/06 argument ICHOIX inutile
!-----------------------------------------------------------------------


      if(vbmax.eq.0) return ! Pas de biologie? alors on quitte la routine.


! LE CAS PARTICULIER DE L'ETAT INITIAL
! 1. Definir les parametres d'avancee dans le temps (date repere & echantillonage)
! Debut:
      if(iteration3d==kount0) then !.................!       !18-04-11

         i1=bioatm_date(1)
         i2=bioatm_date(2)
         i3=bioatm_date(3)
         i4=bioatm_date(4)
         i5=bioatm_date(5)
         i6=bioatm_date(6)
!         write(6,*)'appel à date_to_kount dans river_bio_upd'
!        call date_to_kount
         call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!        bioatm_dt(kr,2)=xdtk_out                     ! date repere NC=1 en kount
!        bioatm_dt(kr,1)=bioatm_info(kr)*3600./dti_fw ! duree, en kount, entre 2 infos
         bioatm_dt(2)=elapsedtime_out              !18-04-11
         bioatm_dt(1)=bioatm_info*3600.        !18-04-11


      endif                    !.................!
! Fin.

      ksecu=0 ! A priori on n'ouvre pas le fichier de "donnee riviere"

! Est on de part et d'autre d'une nouvelle echeance?
      x2=(elapsedtime_now-bioatm_dt(2))/bioatm_dt(1)     !18-04-11
      x1=(elapsedtime_bef-bioatm_dt(2))/bioatm_dt(1)
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
           atm_bio(0)=atm_bio(2)
          enddo
      endif                       !---------------->



! Ouverture du fichier:
      open(unit=3,file=bioatm_file                                    &
                 ,access='direct'                                      &
                 ,recl=4                                               &
!                ,recl=4*kmax                                          &
                 ,form='unformatted')


! lecture du fichier:
      do k3=k4,2,2
         x2=(elapsedtime_now-bioatm_dt(2))/bioatm_dt(1)      !18-04-11
         nc=int(1.+x2)+k3/2
         read(3,rec=nc)atm_bio(k3)
      enddo

      close(3)
      endif               !*************************************>



! Interpolation temporelle entre 2 echeances du fichier:
      rap=      (elapsedtime_now-bioatm_dt(2))/bioatm_dt(1)  & !18-04-11
          -int( (elapsedtime_now-bioatm_dt(2))/bioatm_dt(1) )

       atm_bio(1)=(1.-rap)*atm_bio(0)                  &
                 +    rap *atm_bio(2)


      return
      end
