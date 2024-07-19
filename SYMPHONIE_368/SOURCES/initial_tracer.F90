      subroutine initial_tracer
!______________________________________________________________________
! SYMPHONIE ocean model
! release 303 - last update: 19-07-21
!______________________________________________________________________
!...............................................................................
!    _________                    .__                  .__             ! m[°v°]m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal
      use module_parallele !#mpi
      implicit none
      namelist/notebook_tracer1/imodeltrc,vbmax,ksomax
      namelist/notebook_tracer2/wsed,river_bio,socard,sodate,tauR,radionucleide
#ifdef synopsis
       subroutinetitle='initial_tracer'
       subroutinedescription= &
       'Reads notebook_tracer. Initial state for passive tracer fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!______________________________________________________________________
! Version Date      Description des modifications
!         29/08/01: si strada n'est pas utilisé on ne lit pas le
!                   notebook_tracer. Ce qui évite de "planter" le modèle
!                   si l'utilisateur n'a pas pris le soin de bien
!                   construire son notebook_tracer
!         05/09/01: modification de l'ordre de lecture des profondeurs
!                   min et max de l'immersion initial du traceur (l'odre
!                   croissant est plus logique)
!         25/09/01: amenagements sur CALL DATE_TO_KOUNT
!         14/12/01: Bienvenue à ITIMEBIO
!         17/12/01: amenagements pour cas forward pour economie de RAM
!         22/11/02: annulation du schema forward
!         25/12/02: mise en place filiere restart pour les traceurs
!         27/01/03: debugueur conflits taille memoire
!         19/03/03: date de debut et de fin de calcul pour traceur
!                   & reset du tableau TENDANCEBIO_Z
!         02/12/04: suite du point precedent: si aucune date n'est précisée
!                   dans le notebook_tracer alors les dates de debut et fin
!                   de calcul des traceurs sont pris à kount0 et Ks (voir
!                   time_step.F)
!         23/06/05: Par defaut le traceur est initialisé à zéro
!         06/04/06: RIVER_BIO: un 3eme indice
!         18/04/06: fonctions compatibles avec double precision
!         10/05/06: fonctions compatibles avec double precision
!         18-06-09: suppression KOUNT_BIO
!         23-06-09: suite du point precedent: la lecture des dates associees
!                   à kount_bio sont commentees
! 2010.8  03-05-10  Terminologie
!         09-05-10  seul le proc 0 ecrit messages
! 2010.10 23-06-10  hot_restart renommé dyn_restart et bio_restart
! 2010.13 21-10-10  latlon_to_ij(1) devient latlon_to_ij('glb')
!         03-11-10  des arguments passés dans date_to_kount
! 2010.20 16-04-11  Calculs sur la base d'un temps en secondes
! 2010.25 26-02-12  supression dim_bio et onoff_bio
!         18-03-12  sources traceurs, cas mpi
! S.26    04-04-14  debug lecture de notebook_tracer
!         27-06-14  Encore plus de debug lecture de notebook_tracer!
!         01-12-16  Ajout d'une composante "radioactif" dans les traceurs via: 
!                   radionucleide(vb), alpha_radio(vb) et TauR(vb)
!         12-02-17  wsed dimension verticale
!         08-05-18  affichage ecran
! v249    05-03-19  initialisation de tableaux etc...
! v296    09-02-21  biobc_type(:)=0 
! v303    19-07-21  alpha_radio=0. !19-07-21 ! reset
!______________________________________________________________________



!*******************************************************************************
! LECTURE du notebook_tracer
! DEBUT:
!*******************************************************************************
      itimebio=dim_timebio ! 0=leapfrog 1=forward                      !26/12/02
! section debug debut:                                                 !26/12/02
      if(dim_timebio.ne.0.and.dim_timebio.ne.1) then !%%%%%%>
      write(6,*)'erreur dans parameter,'
      write(6,*)'dim_timebio doit etre egal a 0 ou a 1'
      stop ' dans initial_bio.f!'
      endif                                          !%%%%%%>
      if(itimebio.ne.dim_timebio) then !§§§§§§§§§§§§§§§§§§§§>
      write(6,*)'erreur dans parameter,'
      write(6,*)'mettre dim_timebio égal à ',itimebio,' puis relancer.'
      stop ' dans initial_bio.f!'
      endif                            !§§§§§§§§§§§§§§§§§§§§>
! section debug fin.                                                   !26/12/02

! notebook_tracer (Part I)
      imodeltrc=0 ; vbmax=1 ; ksomax=0
      open(100,file=nomfichier(10)) ! notebook_tracer
      read(100,nml=notebook_tracer1)
      close(100)
      vbmax=vbmax*imodeltrc
      if(imodeltrc==0)return                                         !29/08/01
      call biology_allocate
      if(ksomax>0) then !>>>
        allocate(socard(-1:9,ksomax))  ; socard=0. !05-03-19
        allocate(sodate(6,2,ksomax))   ; sodate=0. !05-03-19
      endif             !>>>
!     biobc_type(:)=1
      biobc_type(:)=0 !09-02-21

! notebook_tracer (Part II)
      wsed=0.
      open(100,file=nomfichier(10)) ! notebook_tracer
      read(100,nml=notebook_tracer2)
      close(100)

! Si radionucleide on initialise alpha_radio sinon il reste à zero
      alpha_radio=0. !19-07-21 ! reset
      do vb=1,vbmax
        if (radionucleide(vb)==1) alpha_radio(vb) = log(2.)/TauR(vb)/86400.
      enddo


! securité au cas ou les étourdis se trompent de signe:
      wsed(1,1:vbmax)=-abs(wsed(1,1:vbmax)) ! on est donc bien sur que wsed < 0 !12-02-17

      do kso=1,ksomax ! kso-kso-kso->

      i1=nint(socard(0,kso))
      if(i1==1) then  !----->
        latit1=socard(1,kso)*pi/180.
        longi1=socard(2,kso)*pi/180.
        call latlon_to_ij('loc')                                        !21-10-10
        socard(1,kso)=deci ; socard(2,kso)=decj
      else            !----->
        socard(1,kso)=socard(1,kso)-par%timax(1)
        socard(2,kso)=socard(2,kso)-par%tjmax(1)  !18-03-12
      endif           !----->

      i2=nint(socard(1,kso)) ; j2=nint(socard(2,kso))
      if(i2>=1.and.i2<=imax.and.j2>=1.and.j2<=jmax) then!2222222>
       if(mask_t(nint(socard(1,kso)),nint(socard(2,kso)),kmax+1).eq.0)then
       write(6,*)'attention source n°',kso,' en terre'
       write(6,*)'ses indices i j sont'                                &
        ,nint(socard(1,kso)),nint(socard(2,kso))
       stop 'donc dans initial_tracer.f'
      endif
      endif                                             !2222222>

      i1=nint(socard(-1,kso))
           x1=socard( 3,kso)
           x2=socard( 4,kso)

      i2=nint(socard(1,kso)) ; j2=nint(socard(2,kso))
      if(i2>=1.and.i2<=imax.and.j2>=1.and.j2<=jmax) then!2222222>
      if(i1.eq.1) then  ! ---------->
        k1= 9999
        k2=-9999
        do k=1,kmax
         if(depth_t(i2,j2,k).ge.x1)k1=min0(k1,k)                      !05/09/01
         if(depth_t(i2,j2,k).le.x2)k2=max0(k2,k)                      !05/09/01
        enddo
        k1=min0(k1,k2)
        socard(3,kso)=min0(kmax,max0(1,k1))
        socard(4,kso)=min0(kmax,max0(1,k2))

      else              ! ---------->

        socard(3,kso)=x1
        socard(4,kso)=x2

      endif             ! ---------->
      endif                                             !2222222>


      if(sodate(1,1,kso)/=0)                                           & !05-03-19
      call datetokount(sodate(1,1,kso),sodate(2,1,kso),sodate(3,1,kso) &
                      ,sodate(4,1,kso),sodate(5,1,kso),sodate(6,1,kso))
      socard(8,kso)=elapsedtime_out


      if(sodate(1,2,kso)/=0)                                           & !05-03-19
      call datetokount(sodate(1,2,kso),sodate(2,2,kso),sodate(3,2,kso) &
                      ,sodate(4,2,kso),sodate(5,2,kso),sodate(6,2,kso))
      socard(9,kso)=elapsedtime_out


!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
!#endif
!      stop 'BP'
      enddo           ! kso-kso-kso->


      if(imodeltrc.eq.1) then !*************************************>

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=4,file=trim(tmpdirname)//'messages',position='append')
      write(4,*)'-----------------------------------------------------'
      write(4,*)'subroutine initial_tracer:'
      write(4,*)
      write(4,*)'model_ traceurs strada (bio, sed..) activé'
      write(4,*)'nbre de variables: ',vbmax
      do kso=1,ksomax
      write(4,*)'source test dispersion n°',kso
      write(4,*)'i & j                 : ',socard(1,kso),socard(2,kso)
      write(4,*)'kmin et max           : ',socard(3,kso),socard(4,kso)
      write(4,*)'variables min et max  : ',socard(5,kso),socard(6,kso)
      write(4,*)'concentration         : ',socard(7,kso)
      write(4,*)'iterations min et max : ',socard(8,kso),socard(9,kso)
      enddo
      write(4,*)                                                       !02/12/04
      write(4,*)'n° d''itération de début et fin de calcul des traceurs'
      do vb=1,vbmax
!     WRITE(4,*)'Traceur n°',KB,' 1ere et derniere iterations:'
!    & ,NINT(KOUNT_BIO(KB,1)),NINT(KOUNT_BIO(KB,2))
      enddo

      close(4)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

      endif                   !*************************************>

!*******************************************************************************
! LECTURE du notebook_tracer
! FIN.
!*******************************************************************************


!*******************************************************************************
! INITIALISATION DES TABLEAUX:
! DEBUT:
!*******************************************************************************

      do 173 vb=1,vbmax
      do 173  k=1,kmax
      do 173  i=1,imax
      do 173  j=1,jmax
      tendancebio_t(i,j,k,vb)=0.                                       !19/03/03
      bio_t(i,j,k,vb)=0.                                               !23/06/05
!     bio_t(i,j,k,vb)=1.
  173 continue

! Cas d'un départ d'une simulation précédente:

      if(initial.eq.1) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!     write(6,*)
!     write(6,*)'....................................................'
!     write(6,*)'dans initial_bio.f: êtes vous sur que le fichier'
!     write(6,*)'restart du model_ bio (chanelbio_in) existe?'
!     write(6,*)'sinon commenter l`appel a hot_restart(11)'

!     call hot_restart(11)                                             !25/12/02
      call bio_restart('r') !23-06-10

!     write(6,*)'la lecture du fichier restart bio s`est bien passée,'
!     write(6,*)'on continue...'
!     write(6,*)'....................................................'


      endif                 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!*******************************************************************************
! INITIALISATION DES TABLEAUX:
! FIN.
!*******************************************************************************

      end subroutine initial_tracer
