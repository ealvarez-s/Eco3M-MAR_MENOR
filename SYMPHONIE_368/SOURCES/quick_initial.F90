      subroutine quick_initial
!______________________________________________________________________
! SYMPHONIE ocean model
! release 310  - last update: 03-11-21
!______________________________________________________________________
!...............................................................................
!    _________                    .__                  .__                     !m[°v°]m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

      use module_principal
      use module_parallele   ! 06-12-13
      use module_q
      implicit none
#ifdef synopsis
       subroutinetitle='quick_initial'
       subroutinedescription='Simple initial state'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date    Descriptions des moficiations:
!         22/07/01: passage a la coordonnae sigma ganaralisae
!         22/02/02: appel a density.F
!         22/02/03: re-ecriture de la partie: quantite volumique decoule
!                   de la quantite primaire. (Plus propre maintenant)
!         04/11/03: passage d'un argument dans CALL DENSITY
!         16/01/04: par dafaut T=T0 et S=S0 pour que initial_state_eq.F ne
!                   perturbe pas les paramatres de l'aquation d'atat.
!         14/11/04: l'initialisation de la densita se fait dans initial_state_eq.F
!         25/08/05: ajout initialisation de la sse
!         02/07/06: equation tke forward entraine elimination du tableau tkenhz_z
!         05/06/07: initialisation de l'integralite de TOBC et SOBC pour eviter
!                   bug aux OBC en cas de non-ogcm aux frontieres
!         13/05/08: sauvegarde TEM_Z dans THZ_Z(2) (et aussi S) pour schema
!                   d'advection
! 2009.2  03-09-09: L'initialisation se fait par les tableaux obc "au temps 1".
!                   ce qui suppose de revoir l'appel a initial_with_obc dans initial_main.F90
!         04-09-09: Suite du point precedent: initialiser egalement les temps 0 & 2
! 2010.11 16-07-10  temobc & salobc renommas temobc & salobc
! 2010.12 03-10-10  suppression obcmin et obcmax
! S.26    28-06-13  use module_parallele
!         06-12-13  fonction max sur sshobc (wetdriyng)
!         22-10-15  ajout temref et salref
!         14-11-15  ajout zone seche
!         16-11-15  prise en compte flag_refstate
!         11-12-15  reset salinitE A s0
!         14-01-16   ajout du cas test baroclinic_jet
!         20-01-16   baroclinic_jet: la nouvelle moyenne zonale dispense de calculer
!                   temobc "non perturbE"
!         09-02-16  amenagement initialisation T et S
!         16-02-16  correction double precision
!         12-12-16  interpolation "2 colonnes" de l'etat de reference.
!                   Modification du cas de la surface (l'algo precedent
!                   conduisant A une singularitE pour dz tendant vers zero)
!         03-01-18  modif sur quick_nh
!         05-04-19  modif sur temobc, salobc
! v260    17-10-19  Par defaut quick_nh est vide
! v276    27-03-20  ajout subroutine quick_rho1DV_from_file
! v287    20-07-20  if(iobc_ogcm==0)call quick_ts1DV_from_file 
! v310    03-11-21  call tidal_ogcm_mix(k1) !03-11-21
! v363    08-01-23  ajout subroutine quick_initial_marmenor_polygon
!...............................................................................



!     call quick_baroclinic_jet !14-0-16
!     return

!***********************************************************************
! INITIALISATION SIMPLE ET RAPIDE DES PRINCIPAUX TABLEAUX:
! DEBUT:
!***********************************************************************

!_______________________________________________________________________________
! pour modifier l'etat initial on intervient a partir de cette ligne:

!..........
! Initialisation des tableaux de la famille "obc"

      do k1=0,2             !03-10-10
! Elevation de surface + epaisseur colonne d'eau                      !25/08/05
      do j=0,jmax+1
      do i=0,imax+1
       sshobc_w(i,j,k1)=0.
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
! ON PREND UNE FAIBLE ELEVATION DE SURFACE EXPRES CAR LE BILAN D'ENERGIE EST
! LINEARISE ET SUPPOSE H+SSH=H
!      sshobc_w(i,j,k1)=0.1*sin(2.*pi*real(i1-2)/real(iglb-2)*1) &
!                          *sin(2.*pi*real(j1-2)/real(jglb-2)*1)  
!      sshobc_w(i,j,k1)=max(real(sshobc_w(i,j,k1),kind=8),1.d-3*wetdry_cst3-h_w(i,j))!06-12-13
!      sshobc_w(i,j,k1)=0.1*exp(-(real(i1-iglb/2)/10.)**2)
!      sshobc_w(i,j,k1)=0.1*exp(-(real(j1-jglb/2)/10.)**2)*exp(-(real(i1-iglb/2)/10.)**2)
      enddo
      enddo
      enddo


! T et S:
! Profil de reference construit sur une grille z
! On definit les niveaux z
      if(flag_refstate==1) then !fffffff> !16-11-15
#ifdef bidonref

       call quick_zref
! Un etat independant de x,y sur cette grille z:
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
        temref_z(i,j,k)=10.+10.*exp(zref_z(k)/500.)
        salref_z(i,j,k)=38.5-1.*exp(zref_z(k)/500.)
       enddo       ; enddo         ; enddo
! Interpolation sur la grille s:
       id_tem=0 ; id_sal=1
       call quick_tsref_interp(0) ! return anyv3d(i,j,k,id_tem) anyv3d(i,j,k,id_sal)
!..........
       if(allocated(temref_t)) then !>>> !22-10-15
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         temref_t(i,j,k)=anyv3d(i,j,k,id_tem)
         salref_t(i,j,k)=anyv3d(i,j,k,id_sal)
        enddo       ; enddo              ; enddo
       endif                        !>>>
! Tableau OBC
       do k1=0,2             !03-10-10
       do j=0,jmax+1
       do i=0,imax+1
        do k=kmax,1,-1
! temperature
!        temobc_t(i,j,k,k1)=i2dh* t0     
         temobc_t(i,j,k,k1)=anyv3d(i,j,k,id_tem)
! salinite:
!        salobc_t(i,j,k,k1)=i2dh* s0     
         salobc_t(i,j,k,k1)=anyv3d(i,j,k,id_sal)
         enddo
       enddo
       enddo
       enddo
#endif

      else                      !fffffff> !16-11-15

! Tableau OBC
       do k1=0,2             !09-02-16
       do j=0,jmax+1
       do i=0,imax+1
        do k=kmax,1,-1
! temperature
         temobc_t(i,j,k,k1)=10. !05-04-19
!        temobc_t(i,j,k,k1)=10.+5.0*exp(depth_t(i,j,k)/50.)
! salinite:
         salobc_t(i,j,k,k1)=35.
!        salobc_t(i,j,k,k1)=35.-0.5*exp(depth_t(i,j,k)/50.)
         enddo
       enddo
       enddo
       enddo
       call tidal_ogcm_mix(0) !03-11-21
       call tidal_ogcm_mix(1) !03-11-21
       call tidal_ogcm_mix(2) !03-11-21

! quick_rho1DV_from_file est une routine pour interpoler un profil 1DV:
!     if(iobc_ogcm==0)call quick_rho1DV_from_file
!     if(iobc_ogcm==0)call quick_ts1DV_from_file !20-07-20

      endif                     !fffffff> !16-11-15
      

! courant moyen:
      do k1=0,2             !03-10-10
      do j=1,jmax+1
      do i=1,imax+1
       velbarobc_u(i,j,k1)=0.                                             !03-09-09
!      velbarobc_u(i,j,k1)=0.5*(sshobc_w(i,j,k1)+sshobc_w(i-1,j,k1))*sqrt(h_u(i,j)/grav)
       velbarobc_v(i,j,k1)=0.
      enddo
      enddo
      enddo

! courant total:
      do k1=0,2             !03-10-10
      do k=1,kmax
       do j=1,jmax
       do i=1,imax+1
!       velobc_u(i,j,k,k1)=0.                                           !03-09-09
        velobc_u(i,j,k,k1)=velbarobc_u(i,j,k1)
       enddo
       enddo
       do j=1,jmax+1
       do i=1,imax
        velobc_v(i,j,k,k1)=0.                                           !03-09-09
       enddo
       enddo
      enddo ! fin de boucle sur k
      enddo ! fin de boucle temps sur k1                               !04-09-09

!     call q_poisson_ondes_diag
!     stop 'toto'

! pour modifier l'etat initial on intervient jusqu'a cette ligne.
!-------------------------------------------------------------------------------

      end subroutine quick_initial

!...............................................................................

      subroutine quick_zref
      use module_principal ; use module_parallele   ! 06-12-13
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='quick_zref'
       subroutinedescription='Simple initial state for zref'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Profil de reference construit sur une grille z
! On definit les niveaux z
      x0=0.1
      k=kmax ; zref_z(k)=0.
      do loop_=1,2
      do k=kmax-1,1,-1
       if(loop_==1) then !lolololo>
        if(k==kmax-1) then
         zref_z(k)=-1.
        else
         zref_z(k)=(2.+x0)*zref_z(k+1)-(1.+x0)*zref_z(k+2)
        endif
       else              !lolololo>
         zref_z(k)=-hmax*zref_z(k)/zref_z(1)
       endif             !lolololo>
      enddo
      enddo
! Pour eviter que l'operation de normalisation (avec sa division) fasse
! que zref_z(1) ne soit pas egal a hmax (l'ennui c'est que ce niveau
! soit legerement ce qui pose pb pour l'interpolation de l'ogcm) on
! ajoute un petit delta:
      zref_z(1)=zref_z(1)-0.01

      end subroutine quick_zref

!...........................................................................
#ifdef bidonref
      subroutine quick_tsref_interp(case_)
      use module_principal ; use module_parallele 
      implicit none
      integer case_,shift_
#ifdef synopsis
       subroutinetitle='quick_tsref_interp'
       subroutinedescription=''
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      if(case_==0) then !00000000000000000>

!      do k=kmax,1,-1
!      write(64,*)zref_z(k),temref_z(1,1,k),salref_z(1,1,k)
!      enddo

       do j=0,jmax+1
       do i=0,imax+1

! Cas zone seche: !14-11-15
       if(depth_w(i,j,1)>=zref_z(kmax)) then !zzzzzzzzzzzzz>

        do k=1,kmax
         anyv3d(i,j,k,id_tem)=temref_z(i,j,kmax)
         anyv3d(i,j,k,id_sal)=salref_z(i,j,kmax)
        enddo

       else                                  !zzzzzzzzzzzzz>

! On ajuste la surface du champ de reference a celle du modele en
! rajoutant/otant la petite couche d'ecart dans laquelle T et S sont donnees
! par les valeurs de surface
       sum1=0.  ; sum2=0. 
       zref_z(kmax)=depth_w(i,j,kmax+1) !12-12-16
!      sum1=(depth_w(i,j,kmax+1)-zref_z(kmax))*temref_z(i,j,kmax)
!      sum2=(depth_w(i,j,kmax+1)-zref_z(kmax))*salref_z(i,j,kmax)
       sum3=0.  ; sum4=0. 
       k0=kmax-1

       do k=kmax,1,-1
     
        do while (zref_z(k0)>depth_w(i,j,k))
! Integrale methode trapezes
         sum1=sum1+(zref_z(k0+1)-zref_z(k0))*0.5*(temref_z(i,j,k0+1)+temref_z(i,j,k0))
         sum2=sum2+(zref_z(k0+1)-zref_z(k0))*0.5*(salref_z(i,j,k0+1)+salref_z(i,j,k0))
         k0=k0-1
        enddo

! ici k0 est en dessous de depth_w(i,j,k)
! on veut le reste de l'integrale entre k0+1 et depth_w(i,j,k)
! rap sert a trouver la valeur du profil de ref a z=depth_w(i,j,k):
        rap=(depth_w(i,j,k)-zref_z(k0))/(zref_z(k0+1)-zref_z(k0))


        x1=(zref_z(k0+1)-depth_w(i,j,k))*0.5*( &
           (1+rap)*temref_z(i,j,k0+1)+(1.-rap)*temref_z(i,j,k0))
        x2=(zref_z(k0+1)-depth_w(i,j,k))*0.5*( &
           (1+rap)*salref_z(i,j,k0+1)+(1.-rap)*salref_z(i,j,k0))


         anyv3d(i,j,k,id_tem)=(sum1+x1-sum3)/dz_t(i,j,k,1)
         anyv3d(i,j,k,id_sal)=(sum2+x2-sum4)/dz_t(i,j,k,1)

         sum3=sum3+dz_t(i,j,k,1)*anyv3d(i,j,k,id_tem)
         sum4=sum4+dz_t(i,j,k,1)*anyv3d(i,j,k,id_sal)

!       if(i==22.and.j==22)write(66,*)depth_t(i,j,k),anyv3d(i,j,k,id_tem),anyv3d(i,j,k,id_sal)

        enddo ! k

        zref_z(kmax)=0. !remettre la valeur par defaut!12-12-16

       endif                                 !zzzzzzzzzzzzz>

       enddo
       enddo

      return
      endif             !00000000000000000>

      if(case_==1.or.case_==2) then !----------------->

      if(case_==1)shift_=-1
      if(case_==2)shift_= 0
      
      do i=2,imax
      i1=i+shift_
      do j=2,jmax-1

! Cas zone seche: !14-11-15
       if(depth_w(i1,j,1)>=zref_z(kmax)) then !zzzzzzzzzzzzz>

        do k=1,kmax
         anyv3d(i1,j,k,id_tem)=0.5*(temref_z(i-1,j,kmax)+temref_z(i,j,kmax))
         anyv3d(i1,j,k,id_sal)=0.5*(salref_z(i-1,j,kmax)+salref_z(i,j,kmax))
        enddo

       else                                   !zzzzzzzzzzzzz>

! On ajuste la surface du champ de reference a celle du modele en
! rajoutant/otant la petite couche d'ecart dans laquelle T et S sont donnees
! par les valeurs de surface
       sum1=0.  ; sum2=0. 
       zref_z(kmax)=depth_w(i1,j,kmax+1) !12-12-16
!      sum1=(depth_w(i1,j,kmax+1)-zref_z(kmax))*0.5*(temref_z(i-1,j,kmax)+temref_z(i,j,kmax))
!      sum2=(depth_w(i1,j,kmax+1)-zref_z(kmax))*0.5*(salref_z(i-1,j,kmax)+salref_z(i,j,kmax))
       sum3=0.  ; sum4=0. 

       k0=kmax-1

       do k=kmax,1,-1

        do while (zref_z(k0)>depth_w(i1,j,k))
! Integrale methode trapezes
         sum1=sum1+(zref_z(k0+1)-zref_z(k0))*0.25*(       &
                temref_z(i  ,j,k0+1)+temref_z(i  ,j,k0)   &
               +temref_z(i-1,j,k0+1)+temref_z(i-1,j,k0) )
         sum2=sum2+(zref_z(k0+1)-zref_z(k0))*0.25*(       &
                salref_z(i  ,j,k0+1)+salref_z(i  ,j,k0)   &
               +salref_z(i-1,j,k0+1)+salref_z(i-1,j,k0) )
         k0=k0-1
        enddo

! ici k0 est en dessous de depth_w(i1,j,k)
! on veut le reste de l'integrale entre k0+1 et depth_w(i1,j,k)
! rap sert a trouver la valeur du profil de ref a z=depth_w(i1,j,k):
        rap=(depth_w(i1,j,k)-zref_z(k0))/(zref_z(k0+1)-zref_z(k0))

        x1=(zref_z(k0+1)-depth_w(i1,j,k))*0.25*(                      &
           (1+rap)*temref_z(i  ,j,k0+1)+(1.-rap)*temref_z(i  ,j,k0)   &
          +(1+rap)*temref_z(i-1,j,k0+1)+(1.-rap)*temref_z(i-1,j,k0))

        x2=(zref_z(k0+1)-depth_w(i1,j,k))*0.25*(                      &
           (1+rap)*salref_z(i  ,j,k0+1)+(1.-rap)*salref_z(i  ,j,k0)   &
          +(1+rap)*salref_z(i-1,j,k0+1)+(1.-rap)*salref_z(i-1,j,k0))

         anyv3d(i1,j,k,id_tem)=(sum1+x1-sum3)/dz_t(i1,j,k,1)
         anyv3d(i1,j,k,id_sal)=(sum2+x2-sum4)/dz_t(i1,j,k,1)
                
         sum3=sum3+dz_t(i1,j,k,1)*anyv3d(i1,j,k,id_tem)
         sum4=sum4+dz_t(i1,j,k,1)*anyv3d(i1,j,k,id_sal)

!       if(i==22.and.j==22)write(66,*)depth_t(i1,j,k),anyv3d(i1,j,k,id_tem),anyv3d(i1,j,k,id_sal)

       enddo  ! k

        zref_z(kmax)=0. !remettre la valeur par defaut !12-12-16
       endif                                  !zzzzzzzzzzzzz>
      enddo 
      enddo 

      return
      endif                         !----------------->

      if(case_==3.or.case_==4) then !----------------->

      if(case_==3)shift_=-1
      if(case_==4)shift_= 0
      
      do j=2,jmax
      j1=j+shift_
      do i=2,imax-1

! Cas zone seche: !14-11-15
       if(depth_w(i,j1,1)>=zref_z(kmax)) then !zzzzzzzzzzzzz>

        do k=1,kmax
         anyv3d(i,j1,k,id_tem)=0.5*(temref_z(i,j-1,kmax)+temref_z(i,j,kmax))
         anyv3d(i,j1,k,id_sal)=0.5*(salref_z(i,j-1,kmax)+salref_z(i,j,kmax))
        enddo

       else                                   !zzzzzzzzzzzzz>

! On ajuste la surface du champ de reference a celle du modele en
! rajoutant/otant la petite couche d'ecart dans laquelle T et S sont donnees
! par les valeurs de surface
       sum1=0.  ; sum2=0. 
       zref_z(kmax)=depth_w(i,j1,kmax+1) !12-12-16
!      sum1=(depth_w(i,j1,kmax+1)-zref_z(kmax))*0.5*(temref_z(i,j-1,kmax)+temref_z(i,j,kmax))
!      sum2=(depth_w(i,j1,kmax+1)-zref_z(kmax))*0.5*(salref_z(i,j-1,kmax)+salref_z(i,j,kmax))
       sum3=0.  ; sum4=0. 
       k0=kmax-1

       do k=kmax,1,-1

        do while (zref_z(k0)>depth_w(i,j1,k))
! Integrale methode trapezes
         sum1=sum1+(zref_z(k0+1)-zref_z(k0))*0.25*(       &
                temref_z(i,j  ,k0+1)+temref_z(i  ,j,k0)   &
               +temref_z(i,j-1,k0+1)+temref_z(i,j-1,k0) )
         sum2=sum2+(zref_z(k0+1)-zref_z(k0))*0.25*(       &
                salref_z(i,j  ,k0+1)+salref_z(i  ,j,k0)   &
               +salref_z(i,j-1,k0+1)+salref_z(i,j-1,k0) )
         k0=k0-1
        enddo

! ici k0 est en dessous de depth_w(i,j1,k)
! on veut le reste de l'integrale entre k0+1 et depth_w(i,j1,k)
! rap sert a trouver la valeur du profil de ref a z=depth_w(i,j1,k):
        rap=(depth_w(i,j1,k)-zref_z(k0))/(zref_z(k0+1)-zref_z(k0))

        x1=(zref_z(k0+1)-depth_w(i,j1,k))*0.25*(                      &
           (1+rap)*temref_z(i,j  ,k0+1)+(1.-rap)*temref_z(i  ,j,k0)   &
          +(1+rap)*temref_z(i,j-1,k0+1)+(1.-rap)*temref_z(i,j-1,k0))

        x2=(zref_z(k0+1)-depth_w(i,j1,k))*0.25*(                      &
           (1+rap)*salref_z(i,j  ,k0+1)+(1.-rap)*salref_z(i  ,j,k0)   &
          +(1+rap)*salref_z(i,j-1,k0+1)+(1.-rap)*salref_z(i,j-1,k0))

         anyv3d(i,j1,k,id_tem)=(sum1+x1-sum3)/dz_t(i,j1,k,1)
         anyv3d(i,j1,k,id_sal)=(sum2+x2-sum4)/dz_t(i,j1,k,1)
                
         sum3=sum3+dz_t(i,j1,k,1)*anyv3d(i,j1,k,id_tem)
         sum4=sum4+dz_t(i,j1,k,1)*anyv3d(i,j1,k,id_sal)

!       if(i==22.and.j==22)write(65,*)depth_t(i,j1,k),anyv3d(i,j1,k,id_tem),anyv3d(i,j1,k,id_sal)

       enddo  ! k

        zref_z(kmax)=0. !remettre la valeur par defaut !12-12-16
       endif                                  !zzzzzzzzzzzzz>

      enddo 
      enddo 

      return
      endif                         !----------------->

      stop 'Err 204 quick_tsref_interp case_ not found'
      end subroutine quick_tsref_interp
#endif
!...................................................................
      subroutine quick_baroclinic_jet !14-01-16
      use module_principal
      use module_parallele   ! 06-12-13
      implicit none 
      double precision icentre
      double precision rhomaxs, rhomaxn, bgdrhos, bgdrhon, zs1, zn1
      double precision drhos, dzs, dzn, drhosfs, drhosfn, drhosf, z0_ 
      double precision z0p, z_k, drhon, cff_perturb, x_, y_, Fyz
      double precision, dimension(:,:), allocatable :: rhoprof
      double precision lambda_, lambda2_, mssh, gmssh
      double precision surface_, gsurface_, rhobar, totalrho
      double precision asymth,zz,zasym,dzasym !16-02-16
           asymth(zz,zasym,dzasym) =                                   &
                                    (zz-zasym)                         &
                  *sqrt(1+0.5*((zz-zasym)+abs(zz-zasym))**2/dzasym**2) &
                           +zasym

!..........
! Initialisation des tableaux de la famille "obc"
      print*,'verif temp',i2dh, t0
! Parametre de configuration
      allocate(rhoprof(0:kmax,1:2))
      rhoprof(:,:)=0.
! Southern profile      
      rhomaxs=27.75
      bgdrhos=9.8d-6
      zs1=-1000.
      dzs=700.
      drhos=1.4
! Northern profile
      rhomaxn=27.7573
      bgdrhon=9.8d-6
      zn1=-400.
      dzn=300.
! Surface Charney mode
      drhosfs=1.5
      drhosfn=0.0
      z0_=-300.
      drhosf=0.0
      z0p=-110.

      lambda_=1600.*1000.
      lambda2_=2000.*1000.
      cff_perturb=0.02


! 1st step background density  
      do k=1,kmax
        rhoprof(k,1)=rhomaxs-bgdrhos*(depth_t(3,3,k)+h_w(3,3))
        rhoprof(k,2)=rhomaxn-bgdrhon*(depth_t(3,3,k)+h_w(3,3))
      enddo

! 2nd step 
      do k=1,kmax
         z_k=asymth(depth_t(3,3,k),zs1,1.3*dzs)
         rhoprof(k,1)=rhoprof(k,1)-drhos*0.5*(1.+tanh((z_k-zs1)/dzs))
      enddo
      drhon=-(rhoprof(kmax,1)-rhoprof(kmax,2))/(0.5+0.5*tanh((z_k-zn1)/dzn))
      do k=1,kmax
         z_k=asymth(depth_t(3,3,k),zn1,1.3*dzn)
         rhoprof(k,2)=rhoprof(k,2)-drhon*0.5*(1.+tanh((z_k-zn1)/dzn))
      enddo

      do k=1,kmax
         z_k=depth_t(3,3,k)
         rhoprof(k,1)=rhoprof(k,1)-drhosf*(exp((z_k-z0p)/abs(z0p)))/exp(1.)
         rhoprof(k,1)=rhoprof(k,1)                                     &
              -drhosfs*0.5*(1.+tanh((z_k-z0_)/abs(z0_)))/tanh(1.)
         rhoprof(k,2)=rhoprof(k,2)-drhosf*(exp((z_k-z0p)/abs(z0p)))/exp(1.)
         rhoprof(k,1)=rhoprof(k,1)                                     &
              -drhosfn*0.5*(1.+tanh((z_k-z0_)/abs(z0_)))/tanh(1.)
      enddo

! 3rd step fit north and south profil
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

      i1=i+par%timax(1)
      j1=j+par%tjmax(1)
! Garantir periodicite:
      if(i1==iglb)i1=2
      if(i1==1   )i1=iglb-1

      y_=real(j1  )/real(jglb  )-0.5
      x_=real(i1-2)/real(iglb-2)

! add perturbation and attenuation
        y_=y_+cff_perturb*exp(depth_t(i,j,k)/1000.)  &
                        *(exp(-(x_-0.5)**2/0.05))    &
                   *( 0.5*sin(2.*pi*x_)              &
                     +0.5*sin(6.*pi*x_) )

        y_=y_*pi*lambda2_/lambda_ + pi/2.

        if (y_.lt.0) then                 !----->
           Fyz=1.
        elseif (y_.gt.pi) then            !----->
           Fyz=0.  
        else                              !----->
           Fyz=1.-(y_-sin(y_)*cos(y_))/pi
        endif                             !----->

        rhp_t(i,j,k)=1000.+Fyz*rhoprof(k,1)+(1.-Fyz)*rhoprof(k,2)

! x1 est le poid de la sommation special mpi et couplage mode externe:
        x1=mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)  & ! *0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))
        *(depth_t(i,j,k)+h_w(i,j))/hz_w(i,j,1)      !remplacer sigma_w par son equivalent (z+h)/(h+ssh)!28-11-14
        sum2=sum2+x1*rhp_t(i,j,k)
        sum1=sum1+x1

      enddo      
      enddo      
      enddo      
#ifdef parallele
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      sum1=sum1glb ; sum2=sum2glb
#endif
      rho=sum2/sum1

        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         temobc_t(i,j,k,:)=rhp_t(i,j,k)-1000.
         salobc_t(i,j,k,:)=i2dh*s0         
        enddo ; enddo ; enddo 
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         tem_t(i,j,k,:)=rhp_t(i,j,k)-1000.
         sal_t(i,j,k,:)=i2dh*s0         
        enddo ; enddo ; enddo 

!      j=jmax/2 ; k=kmax
!      do i=1,imax
!       write(6,*)i+par%timax(1),tem_t(i,j,k,1)
!      enddo
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!#endif

! construction de la ssh
      do j=1,jmax
      do i=1,imax
        xy_t(i,j,1)=0.
      enddo
      enddo

      do 10 k=kmax,1,-1
      do 10 j=1,jmax
      do 10 i=1,imax
        x2=grav*(rhp_t(i,j,k)-rho)*dz_t(i,j,k,1)
        xy_t(i,j,1)=xy_t(i,j,1)+x2
        anyv3d(i,j,k,2)=xy_t(i,j,1)-0.5*x2
   10 continue

      do j=1,jmax
      do i=1,imax
        sshobc_w(i,j,:)=-anyv3d(i,j,1,2)/grav/rho
      enddo
      enddo

      sum1=0. ; sum2=0.
      do j=1,jmax ; do i=1,imax
       x1=mask_t(i,j,kmax+1)*mask_i_w(i)*mask_j_w(j)
       sum1=sum1+x1  ; sum2=sum2+x1*sshobc_w(i,j,1)
      enddo ; enddo
#ifdef parallele
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      sum1=sum1glb ; sum2=sum2glb
#endif
      x1=sum2/sum1
      do k0=0,2 ; do j=1,jmax ; do i=1,imax
       sshobc_w(i,j,k0)=sshobc_w(i,j,1)-x1
      enddo ; enddo ; enddo

! Deduire le courant geostrophique:
! Step 1 : calculer gradient de pression (rho_)
! Alternative: construire le courant géostrophique:
      do j=1,jmax ; do i=1,imax
       ssh_int_w(i,j,1)=sshobc_w(i,j,1)
!      do k=1,kmax
!       tem_t(i,j,k,1)=temobc_t(i,j,k,1)
!      enddo
      enddo ; enddo
      call pressure_gradient
      call graph_out_geostrophic_current('standard')

!     sum1=0. ; sum2=0.
!     do k=kmax,1,-1 ; do j=1,jmax ; do i=1,imax
!       x1=mask_t(i,j,k)*0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))
!       sum2=sum2+x1*rhp_t(i,j,k)
!       sum1=sum1+x1
!     enddo ; enddo ; enddo
!     write(6,*)sum2/sum1

! Step 2: elements du courant geostrophiques

! Step 4: finalisation (tableaux velobc et masquage)
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
       velobc_u(i,j,k,:)=anyv3d(i,j,k,3)*mask_u(i,j,k)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
       velobc_v(i,j,k,:)=anyv3d(i,j,k,4)*mask_v(i,j,k)
      enddo ; enddo ; enddo

      velbarobc_u(:,:,:)=0.d0
      velbarobc_v(:,:,:)=0.d0

! courant moyen:
      do k=1,kmax
        do j=1,jmax
        do i=1,imax+1
        velbarobc_u(i,j,:)=velbarobc_u(i,j,:)                        &
                              +velobc_u(i,j,k,: )                     &
                           *0.5*(dsig_t(i,j,k)+dsig_t(i-1,j,k))
        enddo
        enddo
        do j=1,jmax+1
        do i=1,imax
        velbarobc_v(i,j,:)=velbarobc_v(i,j,:)                        &
                              +velobc_v(i,j,k,:)                      &
                           *0.5*(dsig_t(i,j,k)+dsig_t(i,j-1,k))
        enddo
        enddo
     
      enddo ! fin de boucle sur k

! Ces lignes ne sont plus utile car la derniere version calcule le rappel A partir
! de la moyenne zonale de temobc-tem (idem pour vel et velbar....)
!20-01-16
#ifdef bidon
! Pour finir on calcule le background (le même profil mais avec perturbation nulle)
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
        y_=(j+par%tjmax(1))/real(jglb)-0.5
        x_=(i+par%timax(1))/real(iglb-2)
        y_=y_*pi*lambda2_/lambda_ + pi/2.
        if (y_.lt.0) then
           Fyz=1.
        elseif (y_.gt.pi) then
           Fyz=0.  
        else
           Fyz=1.-(y_-sin(y_)*cos(y_))/pi
        endif
        temobc_t(i,j,k,:)=Fyz*rhoprof(k,1)+(1.-Fyz)*rhoprof(k,2)
      enddo ; enddo ; enddo      
#endif


!-------------------------------------------------------------------------------

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
#endif

      end subroutine quick_baroclinic_jet
!...................................................................

      subroutine quick_nh
      use module_principal ; use module_parallele ; use module_q
      implicit none 

      RETURN !17-10-19

! Cas du domaine periodique. On ne connait pas la periode, on fixe la longueur d'onde (periodique)
      open(unit=3,file='kvector_input')
       read(3,*)i0
      close(3)
!     x0=real(iglb-2)*dxb/20.      ! longueur d'onde
      x0=real(iglb-2)*dxb/real(i0) ! longueur d'onde
      kvector_=2.*pi/x0
      ssh_amplitude_=0.001

! Cas NH
       c_=sqrt( grav/kvector_ *tanh(kvector_*hmax) )
! Cas Hydro
       if(nhpgf_reduce==0.)c_=sqrt(grav*hmax)

       period_=2.*pi/(kvector_*c_)

       if(par%rank==0) then
       write(6,*)'c_=',c_
       write(6,*)'sqrt(gH)',sqrt(grav*hmax)
       write(6,*)'period_ NH et H',period_,2.*pi/(kvector_*sqrt(grav*hmax))
       write(6,*)'longueur d''onde en m en indices',x0,x0/dxb
       write(6,*)'ssh_amplitude_',ssh_amplitude_
       endif

      x0=2.*pi*(iteration3d-1)*dti_fw/period_  ! +0.5*pi ! (+0.5*pi pour courant nul au depart)
      x1=2.*pi* iteration3d   *dti_fw/period_  ! +0.5*pi    ! (+0.5*pi pour courant nul au depart)
      x2=2.*pi*(iteration3d+1)*dti_fw/period_  ! +0.5*pi    ! (+0.5*pi pour courant nul au depart)

      x3=kvector_*dxb
      do j=0,jmax+1 ; do i=0,imax+1
       i1=modulo(i+par%timax(1),iglb-2) ! pour continuite mpi
       ssh_int_w(i,j,0)=ssh_amplitude_*cos(x3*real(i1-1)-x0)
       ssh_int_w(i,j,1)=ssh_amplitude_*cos(x3*real(i1-1)-x1)
       ssh_int_w(i,j,2)=ssh_amplitude_*cos(x3*real(i1-1)-x2)
      enddo         ; enddo


      do k10=0,2 ! indice du temps
       do j=0,jmax+1 ; do i=0,imax+1
        do k=1,kmax

! A propos de depth_t: on garde la profondeur "au repos" car c'est elle qui a priori correspond
! le mieux a la theorie lineaire
         q_t(i,j,k,k10)=                                               &
        grav*ssh_int_w(i,j,k10)*(cosh(kvector_*(h_w(i,j)+depth_t(i,j,k))) &
                                /cosh(kvector_*(h_w(i,j)               )) &
                                -1.)
        enddo
       enddo         ; enddo
      enddo !k10

! Pour la vitesse on pense Au decallage dans le temps et dans l'espace par
! rapport A la ssh
! Le time stepping "forward-backward est actuellement:
! ssh(0)    u(0)    ssh(1)   u(1)   ssh(2)    u(2)

      do k10=0,2 ! indice du temps
       x1=2.*pi*( real(iteration3d+k10)-0.5 )*dti_fw/period_ 
       do j=0,jmax+1 ; do i=1,imax+1
       i1=modulo(i+par%timax(1),iglb-2) ! pour continuite mpi
! x10 la ssh au point "u"
        x10=ssh_amplitude_*cos(x3*real(i1-1-0.5)-x1)

!      if(j==jmax/2.and.k10==1)                 & 
!       write(10+par%rank,*)x10,0.25*(ssh_int_w(i  ,j,k10)    &
!                                    +ssh_int_w(i-1,j,k10)    &
!                                    +ssh_int_w(i  ,j,k10+1)  &
!                                    +ssh_int_w(i-1,j,k10+1)  )

         do k=1,kmax
! x11 la pression q au point "u"
         x11=grav*x10*(cosh(kvector_*(h_u(i,j)+depth_u(i,j,k))) &
                      /cosh(kvector_*(h_u(i,j)               )) &
                                -1.)                             

          if(nhpgf_reduce==0.)x11=0.

          vel_u(i,j,k,k10)=(x10*grav+x11)/c_

         enddo
       enddo         ; enddo
      enddo !k10

      vel_v=0.

      call obc_int_mpi(-1,12)
      call obc_int_mpi( 0,12)
      call obc_int_mpi( 1,12)
      call obc_int_mpi( 2,12)

      end subroutine quick_nh

!...................................................................
#ifdef bidon
      subroutine quick_rho1DV_from_file !27-03-20
      use module_principal ; use module_parallele
      implicit none 
      double precision , dimension(:,:) , allocatable :: rhp1dvfile

! On ne rentre pas dans cette routine si un forcage "notebook_obcforcing" est prevu
      if(iobc_ogcm/=0) return

! Equation d'etat lineaire sans effet de compression:
  eos_author=0 
  eos_comprs=0  
  eos_linear=1 
  t0=11.6                                     
  s0=35.1                                    
  alp_t= 0.000184280361                          
  alp_s= 0.000757332833                   
  rho=   1026.74524

      if(par%rank==0) then !>>>>

!      open(unit=3,file='CORA4.3-IBI-raw_cluster5_RHO.txt')
       open(unit=3,file='profil_cluster_RHO.txt')
! Passage 1: taille du tableau
        read(3,*) ! Une ligne d'entete....
        k10=0
 836    read(3,*,end=837)x1,x2 
        k10=k10+1
        goto 836
 837    close(3)

      endif                !>>>>
      call mpi_bcast(k10,1,mpi_integer,0,par%comm2d,ierr)

      allocate(rhp1dvfile   (k10,2))

      if(par%rank==0) then !>>>>
!      open(unit=3,file='CORA4.3-IBI-raw_cluster5_RHO.txt')
       open(unit=3,file='profil_cluster_RHO.txt')
        read(3,*) ! Une ligne d'entete....
        do k=k10,1,-1
         read(3,*)rhp1dvfile(k,1),rhp1dvfile(k,2)
        enddo
       close(3)
      endif                !>>>>

      call mpi_bcast(rhp1dvfile,2*k10,mpi_double_precision,0,par%comm2d,ierr)

! DEBUG:
! Si le modele stoppe sur la ligne suivante c'est que l'axe vertical
! n'est pas dans le bon sens:
      do k=1,k10-1
      if(rhp1dvfile(k,1)>rhp1dvfile(k+1,1)) &
      stop 'rhp1dvfile(k,1)>rhp1dvfile(k+1,1) quick_rhp1DV_from_file'
      enddo

! Interpolation de la densite potentielle sur la grille symphonie
      do j=0,jmax+1
      do i=0,imax+1
      do k=1,kmax

        do k1=1,k10-1

         if(depth_t(i,j,k)>rhp1dvfile(k1,1).and.depth_t(i,j,k)<=rhp1dvfile(k1+1,1)) then
          rap=(depth_t(i,j,k)-rhp1dvfile(k1,1))/(rhp1dvfile(k1+1,1)-rhp1dvfile(k1,1))
          rhp_t(i,j,k)=rap*rhp1dvfile(k1+1,2)+(1-rap)*rhp1dvfile(k1,2)
         endif

        enddo
! Cas particulier fond & surface
        if(depth_t(i,j,k)<rhp1dvfile(1,1)  )rhp_t(i,j,k)=rhp1dvfile(1,2)
        if(depth_t(i,j,k)>rhp1dvfile(k10,1))rhp_t(i,j,k)=rhp1dvfile(k10,2)

      enddo
      enddo
      enddo

! Enlever rho
      rhp_t=rhp_t-rho

! hypothese sur la temperature
      temobc_t=t0

! Inversion de la salinite en passant par une equation d'etat lineaire
      x1=-rho*alp_t
      x2= rho*alp_s
      x3=-x1*t0-x2*s0
      do k1=0,2
      do j=0,jmax+1
      do i=0,imax+1
      do k=1,kmax
       salobc_t(i,j,k,k1)=(rhp_t(i,j,k)-x1*temobc_t(i,j,k,k1)-x3)/x2
!      if(par%rank==0.and.k1==1.and.i==imax/2.and.j==jmax/2) &
!        write(10+par%rank,*)k,salobc_t(i,j,k,k1),rhp_t(i,j,k),temobc_t(i,j,k,k1)
      enddo
      enddo
      enddo
      enddo

!     write(10+par%rank,*)salobc_t(imax/2,jmax/2,kmax,1),rhp_t(i,j,k)
!     if(par%rank==256) then
!     do k=k10,1,-1
!     write(64,*)rhp1dvfile(k,1:2)
!     enddo
!     endif
!     if(par%rank==256) then
!     do k=kmax,1,-1
!     write(63,*)depth_t(i,j,k),rhp_t(i,j,k)
!     enddo
!     endif

      deallocate(rhp1dvfile)

      end subroutine quick_rho1DV_from_file
#endif
!...................................................................
#ifdef bidon
      subroutine quick_ts1DV_from_file !20-07-20
      use module_principal ; use module_parallele
      implicit none 
      double precision , dimension(:) , allocatable :: z1dvfile,t1dvfile,s1dvfile

!.................................................................................................
! Initialiser la temperature potentielle et la salinitE A partir d'un fichier 1DV 3 colonnes: z,T,S
! Penser A decommenter l'appel A quick_ts1DV_from_file
!.................................................................................................


! On ne rentre pas dans cette routine si un forcage "notebook_obcforcing" est prevu
      if(iobc_ogcm/=0) return

      if(par%rank==0) then !>>>>

       open(unit=3,file='profil_zts.txt')
! Passage 1: taille du tableau
        read(3,*) ! Une ligne d'entete....
        k10=0
 836    read(3,*,end=837)x1,x2,x3
        k10=k10+1
        goto 836
 837    close(3)

      endif                !>>>>
      call mpi_bcast(k10,1,mpi_integer,0,par%comm2d,ierr)

      allocate(z1dvfile   (k10)) ; z1dvfile=0.
      allocate(t1dvfile   (k10)) ; t1dvfile=0.
      allocate(s1dvfile   (k10)) ; s1dvfile=0.

      if(par%rank==0) then !>>>>
       open(unit=3,file='profil_zts.txt')
        read(3,*) ! Une ligne d'entete....
        do k=k10,1,-1
         read(3,*)z1dvfile(k) & ! profondeur  (m)
                 ,t1dvfile(k) & ! temperature potentielle (degres C)
                 ,s1dvfile(k)   ! salinite
        enddo
       close(3)
      endif                !>>>>
      call mpi_bcast(z1dvfile,k10,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(t1dvfile,k10,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(s1dvfile,k10,mpi_double_precision,0,par%comm2d,ierr)

! DEBUG:
! Si le modele stoppe sur la ligne suivante c'est que l'axe vertical
! n'est pas dans le bon sens:
      do k=1,k10-1
      if(z1dvfile(k)>z1dvfile(k+1)) &
      stop 'z1dvfile(k)>z1dvfile(k+1) quick_ts1DV_from_file'
      enddo


! Interpolation de T et S sur la grille symphonie
      do j=0,jmax+1
      do i=0,imax+1
      do k=1,kmax

        do k1=1,k10-1

         if(depth_t(i,j,k)>z1dvfile(k1).and.depth_t(i,j,k)<=z1dvfile(k1+1)) then
          rap=(depth_t(i,j,k)-z1dvfile(k1))/(z1dvfile(k1+1)-z1dvfile(k1))
          temobc_t(i,j,k,0:2)=rap*t1dvfile(k1+1)+(1-rap)*t1dvfile(k1)
          salobc_t(i,j,k,0:2)=rap*s1dvfile(k1+1)+(1-rap)*s1dvfile(k1)

         endif

        enddo
! Cas particulier fond & surface
        if(depth_t(i,j,k)<z1dvfile(1)  )temobc_t(i,j,k,0:2)=t1dvfile(1)
        if(depth_t(i,j,k)<z1dvfile(1)  )salobc_t(i,j,k,0:2)=s1dvfile(1)
        if(depth_t(i,j,k)>z1dvfile(k10))temobc_t(i,j,k,0:2)=t1dvfile(k10)
        if(depth_t(i,j,k)>z1dvfile(k10))salobc_t(i,j,k,0:2)=s1dvfile(k10)

      enddo
      enddo
      enddo

      deallocate(z1dvfile)
      deallocate(t1dvfile)
      deallocate(s1dvfile)

      end subroutine quick_ts1DV_from_file
#endif
!...................................................................
!#ifdef bidon
      subroutine quick_initial_marmenor_polygon !08-01-23
      use module_principal ; use module_parallele ; use module_polygon
      implicit none
      double precision xlon,xlat

       texte90='../../../MAR_MENOR/BATHYMASK/initial_ts_in_marmenor.kml'
!'
       if(par%rank==0)write(6,'(a,a)')'Polygon file ',trim(texte90)
       call polygon_in(0,irep,0.d0,0.d0,trim(texte90))

       do i=0,imax+1 ; do j=0,jmax+1
        call polygon_in(1,irep,lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg,trim(texte90))
        if(irep==1)salobc_t(i,j,:,:)=44.
       enddo ; enddo

      end subroutine quick_initial_marmenor_polygon
!#endif
