      module module_grid
      use module_principal ; use module_parallele ; use module_global ; use module_s
      use module_webcanals
      implicit none
      include 'netcdf.inc'
      double precision , dimension(:,:) , allocatable :: ir8,jr8
      real*16 :: pi16,rad2deg16,deg2rad16
      integer :: range_lon_max=20 &
                ,range_lat_max=20 &
                ,range_lon        &
                ,range_lat
      double precision , dimension(:) , allocatable :: glob_res1, glob_res2

!______________________________________________________________________
! SYMPHONIE ocean model
! release 359 - last update: 23-11-22
!______________________________________________________________________

!_____________________________________________________________
! Version date      Description des modifications
!         13/08/01: plus d'infos dans le fichier "messages"
!                 + creation d'une carte d'identite de la grille
!                 + appel au sous programme latlon_to_ij(0)
!         19/06/03: modif format d'ecriture des latitudes & longitude
!                   dans fichier "messages"
!         21/08/03: ecriture de la periode inertielle moyenne dans le
!                   fichier messages
!         18/04/06: fonctions compatibles avec double precision
!         17/04/07: Passage e coordonnees curvilignes (voir ajout dx_y et cie...)
!         30/04/07: G defini e partir de latitude.
!         31/05/07: Correction test de convergence pour compatibilite e l'equateur
!         16/08/07: debug ecriture fichier messages
!         03/01/08: amelioration de l'increment iteratif
!         05/06/08: debug condition limite sur DY_Z
!         23/02/09: ajout de DX_R et DY_R
!         10/03/09: G est relie e PHI0 pour la parallelisation
!         13/03/09: une C.L. pour lon et lat en cas de parallelisation
!         02/04/09: augmentation taille des boucles pour parallelisation
!         12/05/09: G est relie e 45. et non plus e PHI0 qui pourrait disparaitre
!                   si des grilles plus sophistiquees que la grille mercator sont
!                   introduites e l'avenir...
!         04-06-09  ajouter le cas de grille dlon dlat constant
!         19-06-09  un critere de convergence plus precis
! 2009.3  05-10-09  ajout d'un "ifdef parallele"
! 2010.2  20-12-09  subroutine coriolis renommee initial_lonlat_dxdy_coriolis
! 2010.3  08-01-10  le pole de la grille n'est pas forcement le pole nord
! 2010.5  30-01-10  calcul iteratif de lat_t remplace par calcul direct
! 2010.6  02-02-10  renomme lon_t lat_t
!         03-02-10  calcul min et max pour lon et lat
! 2010.6  05-02-10  blindage calcul dx dy
! 2010.8  15-03-10  - suite du point precedent + amelioration de la precision
!                   des calculs via fonctions trigo double precision
!                   - cas d'une grille pour une ile typegrid=11
!         09-05-10  seul le proc 0 ecrit interrupteur
! 2010.12 20-09-10  possibilite de calcul en simple precision
! 2010.25 06-06-12  possibilite de calcul en simple precision
!         03-07-12  simplification de calculs induite par angle0=0
! S26.2   02-11-12  Possibilite d'initialiser la grille a partir d'un fichier
!         07-03-13  gridrot apres obc_latlon
!         10-09-13  debug
!         06-11-13  affichages ecran
!         07-06-12  call grid_mpi_obc
!         01-01-14  routine grid_toolbox
!         25-02-14  routine grid_toolbox devient module curvgrdtoolbox
!         02-04-14  lonmin lonmax ... deplace apres la lecture de la grille
!                   nemo et la C.L. obc_lonlat
!         02-05-14  En cas de procedure offline lon lat sont lues dans le
!                   fichier offline grille.nc
!         05-05-14  alarme incoherence notebook
!         26-06-14  correction ouverture fichier netcdf
!         27-06-14  modif sur lecture grille dans le cas ioffline=2
!         02-07-14  nouveaux echanges
!         06-07-14  de meilleurs affichages a l'ecran
!         17-07-14  de meilleurs affichages a l'ecran
!         25-07-14  possibilite d'initialiser lon lat avec fichier txt 4 colonnes i j lon lat
!         21-08-14  dissociation des boucles de calculs sur dx_u, dy_u et dx_v, dy_v
!         11-12-14  quadruple precision pour transformation bipolaire (Thom)
!         17-12-14  exemple de creation de tableaux globlon_t globlat_t ...
!         20-12-14  dsig_t lu a la place de sigma_w
!         28-12-14  ajout subroutine grid_rotation_f_location
!         07-01-15  procedure offline: lire kmin_w dans fichier grille
!         09-01-15  procedure offline: lire kmin_u kmin_v dans fichier grille
!         16-02-15  cas nemoffline entraine typegrid=typegrid_file
!         22-03-15  correction varcount lecture h_w
!         17-04-15  appel obc_lonlat_z2 pour cas nemo offline
!         13-06-15  message debug ecran
!         03-07-15  ajout aiguillage pour eviter des tests "bloquants"
!                   non justifies
!         10-09-16  Possibilite de lire le fichier grille de roms
!         07-02-17  aide au debugage filiere offline
!         15-04-17  dx_t et dy_t dans les reservoirs des fleuves
!         17-11-17  lecture mask format kind=1, h_w format real ==> double
!         19-01-18  ajout fichier lonlat au format name.ijhlonlatmask
!         01-02-18  on distingue name.ijhlonlatmask et name.ijhlonlat
!         20-02-18  deplacement ordre "open" avant boucle k
!         26-03-18  adaptation procedure offline "lecture" A la grille verticale fusionnEe
!         04-05-18  grid_full.nc nom possible pour fichier d'entree offline
!         18-05-18  message ecran
!         22-05-18  modele offline demarre sur fichier "lonlatfile"
!         07-07-18  ajout d'un lien https
!         06-09-18  invdx et cie desormais dans module_grid
!         02-10-18  possibilite de fichier lonlat_5col.txt
!         03-10-18  ajout fplan1_grid fplan2_grid
!         18-10-18  lecture format short
!         12-12-18  schema alternatif equivalent pour calculer l'angle de la grille
!         15-01-19  stop procedure occigen et suppression end=372
!         22-01-19  ajout bornage sur argument de fonction acos
! v245    30-01-19  grid_angle_t(i,j)=atan2(gridrotsin_t(i,j),gridrotcos_t(i,j))
! v252    14-04-19  call grid_webcanals_readfiles('dxdy')
! v253    24-04-19  if(webcanals_list=='none')return 
!         03-05-19  ajout lon,lat dans les canaux
!         04-05-19  tetes de canaux peuvent etre dans la zone de chevauchement mpi
!         05-05-19  choix du format d'entree entre i,j,h,dx,dy et i,j,h,dx,dy,lon,lat
! v257    17-06-19  Attention on ne cherche pas le rank si pas de connexion
!         18-06-19  ajout subroutine grid_webcanals_consistency
!         23-06-19  - reecriture compacte
!                   - continuite pour dx, dy gridtype2 (mer) entre points sender et receiver
! v258    17-07-19  if(i_canalcoord(k1,k3,k4,k5)>0.and.j_canalcoord(k1,k3,k4,k5)>0) then 
! v259    13-09-19  lonmin et cie calculEs A la fin de grid_driver
!         24-09-19  sqrtgrav=sqrt(grav) !24-09-19
! v270    12-12-19  stop
! v272    01-01-20  invgrav=1./grav     
! v273    29-01-20  correction message ecran
!         04-02-20  call grid_lonlat2ij_initial desormais dans set_parameters
! v276    12-03-20  plus de noms reconnus dans les fichiers netcdf
! v290    29-10-20  dim augmentEe pour drifters
! v292    24-11-20  lire dsig_t r8
! v303    14-07-21  - Deduire mask_t de globmask !14-07-21
!                   - Deduire flag_merged_levels de l'attribut general Vertical_coordinate 
!                     du fichier grid.nc 14-07-21
!         28-07-21  - Lire upwindriver_t !28-07-21
! v311    16-11-21  - call webcanals_gt2_to_gt2_any2d devient call webcanals_gt1_to_gt2_any2d
!                   - Si point1 non connectE (forcE par notebook_rivers) il est masquE et dx 
!                     garde sa valeur par defaut
! v324    12-02-22  correction bug lecture texte60
! v341    27-03-22  debug faille latlon2ij
! v345    06-04-22  do j=0,jmax+1 ; do i=0,imax+1 augmentees pour obc_bio
! v349    24-06-22  ecrire lon min max etc... dans tmp/messages
! v350    09-07-22  write(3,*) dx dy min max
! v352    18-08-22  webcanaux: gridrotcos_t etc... deduit de lon_t,lat_t 
! v359    23-11-22  Modifs pour le cas offline=2:
!                   - dx_t(:,:)=max(dx_t(:,:),dxb) 
!                   - dy_t(:,:)=max(dy_t(:,:),dyb) 
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

!$ This routine computes the longitude & latidude of the grid nodes,
!$ the horizontal scale factors and the coriolis parameter

contains

!     subroutine grid_lonlat_dxdy_coriolis
      subroutine grid_driver
      implicit none
#ifdef synopsis
       subroutinetitle='grid_driver'
       subroutinedescription=                                     &
       'Driver of the subroutines computing the horizontal grid,' &
       //' Coriolis parameter, lat-dependent gravity'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      typegrid_monopole=1
      typegrid_bipole=3
      typegrid_file=4

      if(southpole_lon==9999.or.     &
         southpole_lat==9999) then !----->
        typegrid=typegrid_monopole
      else                         !----->
        typegrid=typegrid_bipole
! initialisation des constantes en real*16
      pi16=acos(-1.q0) ; deg2rad16=pi/180.q0 ; rad2deg16=180.q0/pi
      endif                        !----->
      if(lonlatfile/='nofile')typegrid=typegrid_file
! Cas nemo offline:
      if(flag_nemoffline==1)typegrid=typegrid_file !16-02-15

      if(flag_nemoffline==0.and.ioffline==2                   &     !05-05-14
                           .and.typegrid/=typegrid_file) then !>>>> !02-05-14
       write(6,'(a,a,a)')'ioffline=2 & lonlatfile=nofile is'  &
      ,' an inconsistent choice. In notebook_grid lonlatfile' &
      ,' should be the OFFLINE grid file (grille.nc).'
       stop ' Stop subroutine grid_driver'
      endif                                            !>>>>

! Etape 1: calcul des points de latitude et longitude comme
!          dans la version precedente, c.a.d. projection
!          mercator dont les parametres sont donnes dans
!          notebook_grid:
      call grid_lonlat ! included in the present file

! Etape 2:
      if(ioffline/=2)call grid_dxdy  ! dx_t, dy_t

! dx_t & dy_t in rivers reservoirs:
      call set_rivers_reservoir_grid   !15-04-17
! dx_t & dy_t du reseau de canaux:
      call grid_webcanals_readfiles('dxdy')   !14-04-19
      call grid_webcanals_readfiles('lonlat') !02-05-19
      if(nbcanal>0) then !>>>  !23-06-19
! continuite pour dx, dy gridtype2 (mer) entre points sender et receiver !23-06-19
         call webcanals_gt1_to_gt2_any2d('dx_t') !16-11-21
         call webcanals_gt1_to_gt2_any2d('dy_t') !16-11-21
      endif              !>>>  !23-06-19

      flag_stop=0
      do j=0,jmax+1 !21-08-14
      do i=1,imax+1

         dx_u(i,j)=0.5*(dx_t(i,j)+dx_t(i-1,j))
         dy_u(i,j)=0.5*(dy_t(i,j)+dy_t(i-1,j))
       dxdy_u(i,j)=dx_u(i,j)*dy_u(i,j) ! surface (m2) maille u

       if(dx_u(i,j)<=0) then !debug>
        flag_stop=1 
        write(10+par%rank,*)'Err dx_u<=0 dx_u,i,j glb',dx_u(i,j),i+par%timax(1),j+par%tjmax(1)
       endif                 !debug>
       if(dy_u(i,j)<=0) then !debug>
        flag_stop=1 
        write(10+par%rank,*)'Err dy_u<=0 dy_u,i,j glb',dy_u(i,j),i+par%timax(1),j+par%tjmax(1)
       endif                 !debug>

      enddo
      enddo

      do j=1,jmax+1 !21-08-14
      do i=0,imax+1

         dx_v(i,j)=0.5*(dx_t(i,j)+dx_t(i,j-1))
         dy_v(i,j)=0.5*(dy_t(i,j)+dy_t(i,j-1))
       dxdy_v(i,j)=dx_v(i,j)*dy_v(i,j) ! surface (m2) maille v

       if(dx_v(i,j)<=0) then !debug>
        flag_stop=1 
        write(10+par%rank,*)'Err dx_v<=0 dx_v,i,j glb',dx_v(i,j),i+par%timax(1),j+par%tjmax(1)
       endif                 !debug>
       if(dy_v(i,j)<=0) then !debug>
        flag_stop=1 
        write(10+par%rank,*)'Err dy_v<=0 dy_v,i,j glb',dy_v(i,j),i+par%timax(1),j+par%tjmax(1)
       endif                 !debug>

      enddo
      enddo

      do j=1,jmax+1 !21-08-14
      do i=1,imax+1

        dx_f(i  ,j  )=0.25*(                                         & !23/02/09
        dx_t(i  ,j  )                                                &
       +dx_t(i-1,j  )                                                &
       +dx_t(i  ,j-1)                                                &
       +dx_t(i-1,j-1))

        dy_f(i  ,j  )=0.25*(                                         & !23/02/09
        dy_t(i  ,j  )                                                &
       +dy_t(i-1,j  )                                                &
       +dy_t(i  ,j-1)                                                &
       +dy_t(i-1,j-1))

      enddo
      enddo

      do j=0,jmax+1
      do i=0,imax+1

! dxdy_c calcule dans une boucle commencant e 0:
      dxdy_t(i,j)=dx_t(i,j)*dy_t(i,j) ! surface (m2) maille T S ssh

       if(dx_t(i,j)<=0) then !debug>
        flag_stop=1 
        write(10+par%rank,*)'Err dx_t<=0 dx_t,i,j glb',dx_t(i,j),i+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'i,j,local',i,j
        write(10+par%rank,*)'imax,jmax',imax,jmax
       endif                 !debug>
       if(dy_t(i,j)<=0) then !debug>
        flag_stop=1 
        write(10+par%rank,*)'Err dy_t<=0 dy_t,i,j glb',dy_t(i,j),i+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'i,j,local',i,j
        write(10+par%rank,*)'imax,jmax',imax,jmax
       endif                 !debug>

! une fois connue la latitude on trouve le coef. de Coriolis:
      coriolis_t(i,j)=2.*                                               &
                    (2.*pi/24./3600.) & ! rotation de la terre
                    *sin( lat_t(i,j))

      enddo
      enddo

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0) stop 'Err negative mesh size. See fortxxx files'

      do j=0,jmax+1 ; do i=0,imax+1 !06-05-22 augmentees pour obc_bio
       invdxdy_t(i,j)=1./dxdy_t(i,j)
         invdx_t(i,j)=  1./dx_t(i,j)
         invdy_t(i,j)=  1./dy_t(i,j)
      enddo         ; enddo
      do j=0,jmax+1 ; do i=1,imax+1 !29-10-20 dim augmentEe pour drifters
       invdxdy_u(i,j)=1./dxdy_u(i,j)
         invdx_u(i,j)=  1./dx_u(i,j)
         invdy_u(i,j)=  1./dy_u(i,j)
      enddo         ; enddo
      do j=1,jmax+1 ; do i=0,imax+1 !29-10-20 dim augmentEe pour drifters
       invdxdy_v(i,j)=1./dxdy_v(i,j)
         invdy_v(i,j)=  1./dy_v(i,j)
         invdx_v(i,j)=  1./dx_v(i,j)
      enddo       ; enddo
      do j=2,jmax ; do i=2,imax
       invdx_f(i,j)=1./dx_f(i,j)
       invdy_f(i,j)=1./dy_f(i,j)
      enddo       ; enddo

      if(fplan1_grid) then !fffffffffff> !15-06-14!03-10-18
       coriolis_t(:,:)=2.*(2.*pi/24./3600.)*sin(phi0)
      endif                !fffffffffff> !15-06-14

      grav=9.780327*(1.+0.0053024*sin(   45.*pi/180. )**2                 & !12/05/09
                       -0.0000058*sin(2.*45.*pi/180. )**2)
      sqrtgrav=sqrt(grav) !24-09-19
      invgrav=1./grav     !01-01-20

! Min et max des latitudes et longitudes !24-06-22
      x1=+9999.  ; x2=-9999. ; x3=+9999. ; x4=-9999.
      do j=0,jmax+1 ; do i=0,imax+1
        x1=min(x1,lon_t(i,j)*rad2deg)
        x3=min(x3,lat_t(i,j)*rad2deg)
        x2=max(x2,lon_t(i,j)*rad2deg)
        x4=max(x4,lat_t(i,j)*rad2deg)
      enddo ; enddo
      call mpi_allreduce(x1,x5,1,mpi_double_precision,mpi_min,par%comm2d,ierr)
      call mpi_allreduce(x3,x7,1,mpi_double_precision,mpi_min,par%comm2d,ierr)
      call mpi_allreduce(x2,x6,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
      call mpi_allreduce(x4,x8,1,mpi_double_precision,mpi_max,par%comm2d,ierr)

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')                    !13/08/01
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine grid_driver:' !24-06-22
      write(3,*)'lon min max',x5,x6       !24-06-22
      write(3,*)'lat min max',x7,x8       !24-06-22
      close(3)                                                         !13/08/01
      endif                !#mpi-->>-->                       !09-05-10

! Min et max dx, dy, !09-07-22
      x1=+9999.  ; x2=-9999. ; x3=+9999. ; x4=-9999.
      do j=2,jmax-1 ; do i=2,imax-1
        x1=min(x1,dx_t(i,j))
        x3=min(x3,dy_t(i,j))
        x2=max(x2,dx_t(i,j))
        x4=max(x4,dy_t(i,j))
      enddo ; enddo
      call mpi_allreduce(x1,x5,1,mpi_double_precision,mpi_min,par%comm2d,ierr)
      call mpi_allreduce(x3,x7,1,mpi_double_precision,mpi_min,par%comm2d,ierr)
      call mpi_allreduce(x2,x6,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
      call mpi_allreduce(x4,x8,1,mpi_double_precision,mpi_max,par%comm2d,ierr)

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')                    !13/08/01
      write(3,*)                                                       !13/08/01
      write(3,*)'dx min max',x5,x6       !09-07-22
      write(3,*)'dy min max',x7,x8       !09-07-22
      write(3,*)                                                       !13/08/01
      const1=max(coriolis_t(imax/2,jmax/2),small1)                       !16/08/07
      write(3,*)'periode inertielle rank0:(h & j)'                   & !21/08/03
       ,2.*pi/const1/3600.                                           &
       ,2.*pi/const1/3600./24.
      write(3,*)'grav=',grav
      close(3)                                                         !13/08/01
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif


! Grille generalisee necessite de connaitre le masque donc l'appel
! a la subroutine grid_check_lonlat2ij se fait apres lecture du masque
!     call grid_check_lonlat2ij !29-12-13

!..... deplace le 13-09-19
!& min & max values for longitude & latitude: !02-04-14
      lonmin= 1.d10 ; latmin= 1.d10 ; lonmax=-1.d10 ; latmax=-1.d10 !03-02-10
      do j=-1,jmax+2
      do i=-1,imax+2
       lonmin=min(lonmin,lon_t(i,j))                    !03-02-10
       latmin=min(latmin,lat_t(i,j))
       lonmax=max(lonmax,lon_t(i,j))
       latmax=max(latmax,lat_t(i,j))
      enddo
      enddo
      lonmin=lonmin*rad2deg ; lonmax=lonmax*rad2deg     !03-02-10
      latmin=latmin*rad2deg ; latmax=latmax*rad2deg

!     lonmin=minval(lon_t)*rad2deg
!     latmin=minval(lat_t)*rad2deg
!     lonmax=maxval(lon_t)*rad2deg
!     latmax=maxval(lat_t)*rad2deg

      return
      end subroutine grid_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine grid_dxdy
      implicit none
#ifdef synopsis
       subroutinetitle='grid_dxdy'
       subroutinedescription= &
       'Computes horizontal size meshes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(initialgridfile_txt/='none')then !----->
       call grid_dxdy_from_file
       return
      endif                               !----->

      if(fplan2_grid==1) then !fffffffffff> !15-06-14!03-10-18

! Constant horizontal resolution:
       dx_t(:,:)=dxb
       dy_t(:,:)=dyb

      else                    !fffffffffff> !03-07-15

       do j=0,jmax+1   
       do i=0,imax+1   

       x1=lon_u(i+1,j)
       if(x1-lon_u(i,j)<-pi)x1=x1+2.*pi
       if(x1-lon_u(i,j)> pi)x1=x1-2.*pi

! e1=(a+z)*sqrt( (dlon/di x cos(lat))**2 +(dlat/di)**2 ) ! ref: doc opa page 16
        dx_t(i,j)=rayonterre*sqrt(                               & !15-03-10
         ((     x1      - lon_u(i,j))*cos( lat_t(i,j)))**2       & !05-02-10
        +(  lat_u(i+1,j)- lat_u(i,j)                  )**2 )


!       if(dx_t(i,j)==0.) then !>>>>> !13-06-15
!        write(6,*)'Err dx=0 proc i j=',par%rank,i,j                        &
!       ,'lon_u(i+1,j),lon_u(i,j)=',rad2deg*lon_u(i+1,j),rad2deg*lon_u(i,j) &
!       ,'lat_u(i+1,j),lat_u(i,j)=',rad2deg*lat_u(i+1,j),rad2deg*lat_u(i,j) &
!       ,'Suggestion: verifier valeur fplan_grid dans notebook_grid' !18-05-18
!        stop 'Err dx=0 module_grid'
!       endif                  !>>>>>


       x1=lon_v(i,j+1)
       if(x1-lon_v(i,j)<-pi)x1=x1+2.*pi
       if(x1-lon_v(i,j)> pi)x1=x1-2.*pi

! e2=(a+z)*sqrt( (dlon/dj x cos(lat))**2 +(dlat/dj)**2 ) ! ref: doc opa page 16
        dy_t(i,j)=rayonterre*sqrt(                                 & !15-03-10
         ((     x1      - lon_v(i,j))*cos( lat_t(i,j)))**2         & !05-02-10
        +(  lat_v(i,j+1)- lat_v(i,j)                  )**2 )

!       if(dy_t(i,j)==0.) then !>>>>> !13-06-15
!        write(6,*)'Err dy=0 proc i j=',par%rank,i,j                        &
!       ,'lon_v(i,j+1),lon_v(i,j)=',rad2deg*lon_v(i,j+1),rad2deg*lon_v(i,j) &
!       ,'lat_v(i,j+1),lat_v(i,j)=',rad2deg*lat_v(i,j+1),rad2deg*lat_v(i,j) &
!       ,'Suggestion: verifier valeur fplan_grid dans notebook_grid' !18-05-18
!        stop 'Err dy=0 module_grid'
!       endif                  !>>>>>

      enddo
      enddo

      if(fplan2_grid==2) then !m°v°m> !05-10-18

! Imposer d(dx_t)/di=0 et d(dj_t)/di=0
! Principe: on fait la moyenne mpi_sum de dx_t(imax/2,:) et dy_t(imax/2,:)

       allocate(glob_res1(0:jglb+1)) ; glob_res1=0.
       allocate(glob_res2(0:jglb+1)) ; glob_res2=0.
       i=imax/2
       do j=0,jglb+1
        j1=j-par%tjmax(1)
        if(j1>=0.and.j1<=jmax+1)glob_res1(j)=dx_t(i,j1)*mask_j_w(j1)
       enddo
       call mpi_allreduce(glob_res1,glob_res2,jglb+2,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       do j=0,jglb+1
        j1=j-par%tjmax(1)
        if(j1>=0.and.j1<=jmax+1) then
         do i=0,imax+1
          dx_t(i,j1)=glob_res2(j)/real(nbdom_imax)
         enddo
        endif
       enddo

       i=imax/2
       do j=0,jglb+1
        j1=j-par%tjmax(1)
        if(j1>=0.and.j1<=jmax+1)glob_res1(j)=dy_t(i,j1)*mask_j_w(j1)
       enddo
       call mpi_allreduce(glob_res1,glob_res2,jglb+2,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       do j=0,jglb+1
        j1=j-par%tjmax(1)
        if(j1>=0.and.j1<=jmax+1) then
         do i=0,imax+1
          dy_t(i,j1)=glob_res2(j)/real(nbdom_imax)
         enddo
        endif
       enddo

       deallocate(glob_res1)
       deallocate(glob_res2)
      endif                   !m°v°m> !05-10-18

      if(fplan2_grid==3)stop 'fplan2_grid==3 A FAIRE'

      endif               !fffffffffff> !03-07-15

      call grid_mpi_obc !07-12-13

      flag_stop=0
! Etape debug apres grid_mpi_obc !02-10-18
      do j=0,jmax+1 ;  do i=0,imax+1
        if(dx_t(i,j)==0.) then !>>>>> !13-06-15
         write(6,*)'Err 367 dx=0 proc i j loc et glob',par%rank,i,j,i+par%timax(1),j+par%tjmax(1)
         write(6,*)'lon_u(i+1,j)  =',rad2deg*lon_u(i+1,j)
         write(6,*)'lon_u(i  ,j)  =',rad2deg*lon_u(i  ,j)
         write(6,*)'lat_u(i+1,j)  =',rad2deg*lat_u(i+1,j)
         write(6,*)'lat_u(i  ,j)  =',rad2deg*lat_u(i  ,j)
         write(6,*)'Suggestion: verif valeur fplan_grid notebook_grid' !18-05-18
         flag_stop=1 ; goto 416    !15-01-19
!        stop 'Err dx=0 module_grid'
        endif                  !>>>>>
        if(dy_t(i,j)==0.) then !>>>>> !13-06-15
         write(6,*)'Err 374 dy=0 proc i j=',par%rank,i,j,i+par%timax(1),j+par%tjmax(1) &
        ,'lon_v(i,j+1),lon_v(i,j)=',rad2deg*lon_v(i,j+1),rad2deg*lon_v(i,j) &
        ,'lat_v(i,j+1),lat_v(i,j)=',rad2deg*lat_v(i,j+1),rad2deg*lat_v(i,j) &
        ,'Suggestion: verifier valeur fplan_grid dans notebook_grid' !18-05-18
         flag_stop=1 ; goto 416    !15-01-19
!        stop 'Err dy=0 module_grid'
        endif                  !>>>>>

      enddo         ;  enddo
 416  call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'Error 435 in subroutine module_grid'

      if(typegrid==2)then                                            !04-06-09

      sum1=0.
      sum2=0.
      sum3=0.
      do j=1,jmax
      do i=1,imax
       x1=un*mask_i_w(i)*mask_j_w(j) !#MPI ne pas sommer 2 fois la zone de recouvrement
       sum1=sum1+x1
       sum2=sum2+x1*dx_t(i,j)
       sum3=sum3+x1*dy_t(i,j)
      enddo
      enddo
#ifdef parallele
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
#else
      sum1glb=sum1
      sum2glb=sum2
      sum3glb=sum3
#endif

! attention on ne modifie pas dxb et dyb quand typegrid=1 car dans
! ce cas DXB et DYB jouent un rele (interpolation) dans la projection
! de mercator:
       dxb=sum2glb/sum1glb
       dyb=sum3glb/sum1glb
       stop 'la grille dlon dlat pas encore testee!'

      endif


      end subroutine grid_dxdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine grid_lonlat
      implicit none
      integer :: ncid_,istr_=-9999,jstr_=-9999,iend_=-9999,jend_=-9999
#ifdef synopsis
       subroutinetitle='grid_lonlat'
       subroutinedescription= &
       'Computes longitude and latitude (radians) of grid nodes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      i0=grid_i0
      j0=grid_j0

! Compute the Mercator longitudes and latitudes:

      x0=rayonterre*cos(phi0)*longi0                                  !15-03-10
      y0=rayonterre*cos(phi0)*log( (1.d0+sin(latit0))/cos(latit0))    !15-03-10

! Creation d'une "carte d'identite" de la grille:                      !13/08/01
      gridcard(1,1)=phi0
      gridcard(2,1)=latit0
      gridcard(3,1)=longi0
      gridcard(4,1)=angle0
      gridcard(5,1)=real(i0)
      gridcard(6,1)=real(j0)
      gridcard(7,1)=rayonterre
      gridcard(8,1)=x0
      gridcard(9,1)=y0

!#MPI place apres la sauvegarde de I0 et J0 dans GRIDCARD(5,1) et GRIDCARD(6,1)
      i0=i0-par%timax(1)
      j0=j0-par%tjmax(1)

      if(typegrid==typegrid_monopole.or.    &
         typegrid==typegrid_bipole) then !1111111>                   !04-06-09

      const1=sqrt(dxb*dyb)*1.d-12                                      !19-06-09

      do j=-1,jmax+2   ! 0,NECO+1
      do i=-1,imax+2   ! 0,MECO+1                                      !02/04/09

      lon_t(i,j)=(x0+dxb*(i-i0))/(rayonterre*cos(phi0))
      lat_t(i,j)=atan(sinh( (y0+dyb*(j-j0) )/(rayonterre*cos(phi0)) ) )

      enddo
      enddo

#ifdef bidon
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      write(6,*)'bidouille'
      i=imax/2 ; j=jmax/2
      x1=lon_t(i,j)-lon_t(i-1,j)
      write(6,*)'x1=',x1*rad2deg

      sum1=0.
      do j=-1,jmax+2   ! 0,NECO+1
       sum1=sum1-0.3*x1
       do i=-1,imax+2   ! 0,MECO+1                                      !02/04/09
         lon_t(i,j)=lon_t(i,j)+sum1
!        write(64,*)lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
       enddo
      enddo

      sum1=0.
      do i=-1,imax+2   ! 0,MECO+1                                      !02/04/09
       sum1=sum1-0.1*x1
       do j=-1,jmax+2   ! 0,NECO+1
         lat_t(i,j)=lat_t(i,j)+sum1
         write(64,*)lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
       enddo
      enddo

!    stop 'naz'
#endif

      endif                !111111111111111111111111>


      if(typegrid==2) then !222222222222222222222222>                  !04-06-09
      stop ' Option typegrid==2 not available'
!     do j=-1,jmax+2   ! 0,NECO+1
!     do i=-1,imax+2   ! 0,MECO+1
!        lon_t(i,j)=longi0+real(i-i0)*dlon
!        lat_t(i,j)=latit0+real(j-j0)*dlat
!     enddo
!     enddo
      endif                !222222222222222222222222>


      if(typegrid==typegrid_bipole) then !bbbbbb>
        call grid_bipole_transformation
      else                               !bbbbbb>
           if(typegrid==typegrid_monopole) then !mmmmmm>
              call grid_monopole_transformation
           else                                 !mmmmmm>
              if(typegrid/=typegrid_file)stop 'Grid type unknown'
           endif                                !mmmmmm>
      endif                              !bbbbbb>

      if(par%rank==0)write(6,'(a,a)')'Open lonlatfile ',trim(lonlatfile)
!     if(typegrid==typegrid_file) then !fffffff>
      if(typegrid==typegrid_file.and.flag_nemoffline==0) then !fffffff> !16-02-15

       if(index(lonlatfile,'4col.txt')/=0.or.  &      
          index(lonlatfile,'5col.txt')/=0) then !??????????????> !25-07-14!02-10-18

! initialise istr_, jstr_
        open(unit=3,file=trim(lonlatfile))
         read(3,*)istr_,jstr_,x1,x2
        close(3)
        if(istr_==1) then ; iend_=imax   ; jend_=jmax   ; k2=iglb*jglb ; endif
        if(istr_==0) then ; iend_=imax+1 ; jend_=jmax+1 ; k2=(iglb+2)*(jglb+2) ; endif
        open(unit=3,file=trim(lonlatfile)) !20-02-18
        do k=1,k2
!         if(index(lonlatfile,'4col.txt')/=0)read(3,*,end=372)i1,j1,x1,x2
!         if(index(lonlatfile,'5col.txt')/=0)read(3,*,end=372)i1,j1,x0,x1,x2 !02-10-18
          if(index(lonlatfile,'4col.txt')/=0)read(3,*)i1,j1,x1,x2              !15-01-19
          if(index(lonlatfile,'5col.txt')/=0)read(3,*)i1,j1,x0,x1,x2 !02-10-18 !15-01-19
          i=i1-par%timax(1) ; j=j1-par%tjmax(1)
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then
           lon_t(i,j)=x1 ; lat_t(i,j)=x2
          endif
        enddo
  372   close(3)

       else                                    !??????????????>

!      if(index(lonlatfile,'.ijhlonlatmask')/=0) then !m0v0m> !19-01-18
       if(index(lonlatfile,'.ijhlonlat')/=0)     then !m0v0m> !01-02-18

        istr_=1 ; jstr_=1 ; iend_=imax ; jend_=jmax 
        open(unit=3,file=lonlatfile)
        do i1=1,iglb
         read(3,*) ! 'partie carte mask'
        enddo

        k0=0 ; if(index(lonlatfile,'.ijhlonlatmask')/=0)k0=1 !01-02-18
        do i1=1,iglb ; do j1=1,jglb
          if(k0==0)read(3,*)i,j,x1,x2,x3    !01-02-18
          if(k0==1)read(3,*)i,j,x1,x2,x3,k1 !01-02-18
          i=i-par%timax(1) ; j=j-par%tjmax(1)
          if(i>=1.and.i<=imax.and.j>=1.and.j<=jmax) then
           lon_t(i,j)=x2
           lat_t(i,j)=x3
          endif
        enddo       ; enddo


        close(3)

       else                                           !m0v0m> !19-01-18

       if(index(lonlatfile,'.nc')/=0) then !ooo> !19-01-18

       if(par%rank==0)write(6,'(a,a)')'Open lonlatfile ',trim(lonlatfile)
       status=nf_open(trim(lonlatfile),nf_nowrite,ncid_) !26-06-14
       if(status/=0)stop ' Erreur nf_open lonlatfile'

                    status=nf_inq_dimid(ncid_,'ni_t',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_,'nx_t',dim_x_id)   !12-03-20
       if(status/=0)status=nf_inq_dimid(ncid_,'xi_rho',dim_x_id) !10-09-16
       if(status/=0)stop ' stop erreur nf_inq_dimid ni_t 2014'
                    status=nf_inq_dimid(ncid_,'nj_t',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_,'ny_t',dim_y_id)    !12-03-20
       if(status/=0)status=nf_inq_dimid(ncid_,'eta_rho',dim_y_id) !10-09-16
       if(status/=0)stop ' stop erreur nf_inq_dimid nj_t 2014'
                    status=nf_inq_dimlen(ncid_,dim_x_id,i0)
       if(status/=0) stop 'Err 735 nf_inq_dimlen dim_x_id'
                    status=nf_inq_dimlen(ncid_,dim_y_id,j0)
       if(status/=0) stop 'Err 737 nf_inq_dimlen dim_y_id'

       if(i0==iglb)   then ; istr_=1 ; iend_=imax   ; endif
       if(j0==jglb)   then ; jstr_=1 ; jend_=jmax   ; endif
       if(i0==iglb+2) then ; istr_=0 ; iend_=imax+1 ; endif
       if(j0==jglb+2) then ; jstr_=0 ; jend_=jmax+1 ; endif
       if(istr_==-9999) then
        write(10+par%rank,*)'i0 iglb iglb+2',i0,iglb,iglb+2
        write(10+par%rank,*)'j0 jglb jglb+2',j0,jglb,jglb+2
        write(10+par%rank,*)'istr_,iend_',istr_,iend_
        write(10+par%rank,*)'jstr_,jend_',jstr_,jend_
        stop ' grid_lonlat istr iend jstr jend undefined. See fort.xxx'
       endif

       varstart(1)=1+par%timax(1) ; varcount(1)=iend_-istr_+1 ! (imax ou imax+2)
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jend_-jstr_+1 ! (jmax ou jmax+2)

                    status=nf_inq_varid(ncid_,'longitude_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_rho',var_id) !10-09-16
       if(status/=0)stop ' stop nf_inq_varid lon lonlatfile'
       status=nf_get_vara_double(ncid_,var_id,varstart(1:2)   &
                                             ,varcount(1:2)   &
                            ,lon_t(istr_:iend_,jstr_:jend_))
       if(status/=0)stop ' stop nf_get_var_double lon lonlatfile'

                    status=nf_inq_varid(ncid_,'latitude_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_rho',var_id) !10-09-16
       if(status/=0)stop ' stop nf_inq_varid lat lonlatfile'
       status=nf_get_vara_double(ncid_,var_id,varstart(1:2)   &
                                             ,varcount(1:2)   &
                            ,lat_t(istr_:iend_,jstr_:jend_))
       if(status/=0)stop ' stop nf_get_var_double lat lonlatfile'

!#ifdef bidon
       if(ioffline==2) then !ooo> !18-08-22

                      status=nf_inq_varid(ncid_,'dx_t',var_id)
         if(status/=0)stop ' stop nf_inq_varid dx_t lonlatfile'
                      status=nf_inq_vartype(ncid_,var_id,k0)
         if(status/=0)stop ' stop nf_inq_vartype dx_t lonlatfile'
         if(k0/=nf_double) stop ' Error type dx_t lonlatfile'
         status=nf_get_vara_double(ncid_,var_id,varstart(1:2)   &
                                               ,varcount(1:2)   &
                              ,dx_t(istr_:iend_,jstr_:jend_))
         if(status/=0)stop ' stop nf_get_var_double dx_t lonlatfile'

                      status=nf_inq_varid(ncid_,'dy_t',var_id)
         if(status/=0)stop ' stop nf_inq_varid dy_t lonlatfile'
                      status=nf_inq_vartype(ncid_,var_id,k0)
         if(status/=0)stop ' stop nf_inq_vartype dy_t lonlatfile'
         if(k0/=nf_double) stop ' Error type dy_t lonlatfile'
         status=nf_get_vara_double(ncid_,var_id,varstart(1:2)   &
                                               ,varcount(1:2)   &
                              ,dy_t(istr_:iend_,jstr_:jend_))
         if(status/=0)stop ' stop nf_get_var_double dy_t lonlatfile'

         if(istr_/=0.or.jstr_/=0.or.iend_/=imax+1.or.jend_/=jmax+1) then !debug>
!         write(6,*)'istr_,iend_,jstr_,jend',istr_,iend_,jstr_,jend_,imax+1,jmax+1
          stop 'Err: (istr_,iend_,jstr_,jend)=(0,imax+1,0,jmax+1) is expected when ioffline=2 for dx_t & dy_t in grid.nc' 
!'
         endif                                                          !debug>

! Si, dans le cas ioffline==2, dx et dy sont negatifs, c'est qu'on lit
! un fichier de grille excluant les zones continentales inutilisees. On
! impose une valeur positive pour eviter de stopper le modele:
         do j=jstr_,jend_ ; do i=istr_,iend_
          if(dx_t(i,j)<0.)dx_t(i,j)=dxb
          if(dy_t(i,j)<0.)dy_t(i,j)=dyb
         enddo            ; enddo

       endif                !ooo> !18-08-22
!#endif

       status=nf_close(ncid_)

       else                                !ooo> !19-01-18
        stop 'ERR 540 No file for lon,lat !'
       endif                               !ooo> !19-01-18

       endif                                          !m0v0m> !19-01-18

       endif                                   !??????????????>

! Conversion degres radian:
       do j=jstr_,jend_ ; do i=istr_,iend_
        lon_t(i,j)=lon_t(i,j)*deg2rad ; lat_t(i,j)=lat_t(i,j)*deg2rad
       enddo ; enddo

! C.L.
! Selon les bornes des boucles istr_ etc... il faut calculer des conditions aux
! limites en imax+2 seulement, ou imax+1 et imax+2...
       do i=istr_,iend_

        x1=lon_t(i,1)-lon_t(i,2)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
        if(jstr_/=0)lon_t(i,0 )=lon_t(i,1)+x1
                    lon_t(i,-1)=lon_t(i,1)+x1*2.

        x1=lon_t(i,jmax)-lon_t(i,jmax-1)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
        if(jend_/=jmax+1)lon_t(i,jmax+1)=lon_t(i,jmax)+x1
                         lon_t(i,jmax+2)=lon_t(i,jmax)+x1*2.

        x1=lat_t(i,1)-lat_t(i,2)
        if(jstr_/=0)lat_t(i,0 )=lat_t(i,1)+x1
                    lat_t(i,-1)=lat_t(i,1)+x1*2.

        x1=lat_t(i,jmax)-lat_t(i,jmax-1)
        if(jend_/=jmax+1)lat_t(i,jmax+1)=lat_t(i,jmax)+x1
                         lat_t(i,jmax+2)=lat_t(i,jmax)+x1*2.

       enddo

       do j=-1,jmax+2 ! Cette boucle est sur les bornes max (et non jstr_ et jend_).

        x1=lon_t(1,j)-lon_t(2,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
        if(istr_/=0)lon_t(0 ,j)=lon_t(1,j)+x1
                    lon_t(-1,j)=lon_t(1,j)+x1*2.

        x1=lon_t(imax,j)-lon_t(imax-1,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
        if(iend_/=imax+1)lon_t(imax+1,j)=lon_t(imax,j)+x1
                         lon_t(imax+2,j)=lon_t(imax,j)+x1*2.

        x1=lat_t(1,j)-lat_t(2,j)
        if(istr_/=0)lat_t(0 ,j)=lat_t(1,j)+x1
                    lat_t(-1,j)=lat_t(1,j)+x1*2.

        x1=lat_t(imax,j)-lat_t(imax-1,j)
        if(iend_/=imax+1)lat_t(imax+1,j)=lat_t(imax,j)+x1
                         lat_t(imax+2,j)=lat_t(imax,j)+x1*2.

       enddo


! commentE le 02-10-18 car applique 20 lignes plus bas....
!     call obc_lonlat(0)    !29-06-14

      endif                            !fffffff>

!.....

! nemo offline case: flag_nemoffline==1
!     if(initialgridfile_txt/='none')then !----->
      if(flag_nemoffline==1)         then !----->!05-05-14
       call grid_dim_from_file
       call grid_lonlat_from_file
      endif                               !----->

! Dans le cas de la parallelisation, appliquer une condition aux limites
! en i=1, i=imax, j=1, j=jmax et
! en i=0, i=imax+1, j=0, j=jmax+1
      call obc_lonlat(0)                                               !13/03/09

!..... deplace le 02-04-14
!& min & max values for longitude & latitude: !02-04-14
!     lonmin= 1.d10 ; latmin= 1.d10 ; lonmax=-1.d10 ; latmax=-1.d10 !03-02-10
!     do j=-1,jmax+2
!     do i=-1,imax+2
!      lonmin=min(lonmin,lon_t(i,j))                    !03-02-10
!      latmin=min(latmin,lat_t(i,j))
!      lonmax=max(lonmax,lon_t(i,j))
!      latmax=max(latmax,lat_t(i,j))
!     enddo
!     enddo
!     lonmin=lonmin*rad2deg ; lonmax=lonmax*rad2deg     !03-02-10
!     latmin=latmin*rad2deg ; latmax=latmax*rad2deg

#ifdef bidon
      do j=-1,jmax+2
      if(mod(j,2)==0) then
       do i=-1,imax+2                                       !02/04/09
       write(68,'(2(e13.6,1x))') lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
       enddo
      else
       do i=imax+2,-1,-1
       write(68,'(2(e13.6,1x))') lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
       enddo
      endif
      enddo
      do i=-1,imax+2                                       !02/04/09
      if(mod(i,2)==0) then
       do j=-1,jmax+2
       write(68,'(2(e13.6,1x))') lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
       enddo
      else
       do j=jmax+2,-1,-1
       write(68,'(2(e13.6,1x))') lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
       enddo
      endif
      enddo
      stop 'coco3'
#endif

!.............

! Tableaux de rotation si axes non paralleles lignes lon et lat   !07-03-13
      do j=0,jmax+1
      do i=0,imax+1
       x1=lon_t(i+1,j)-lon_t(i-1,j)
       if(x1<-pi)x1=x1+2.*pi
       if(x1> pi)x1=x1-2.*pi
! 07-07-18
! https://docs.google.com/presentation/d/1FmAXNCdY_vL5KCUvkY0AxOQ6u145xQc1q-3cfXL1shU/edit#slide=id.p
       x0=-atan2(lat_t(i+1,j)-lat_t(i-1,j),x1*cos(lat_t(i,j)))     !13-02-13
!      x0=-atan2(  lat_t(i+1,j)-lat_t(i-1,j)                     & ! rotation symphonie
!                 ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j)))    !15-03-10
       gridrotcos_t(i,j)=cos(x0)
       gridrotsin_t(i,j)=sin(x0)
       grid_angle_t(i,j)=atan2(gridrotsin_t(i,j),gridrotcos_t(i,j)) !30-01-19

! Meme algo mais calcul le long de j:
!      x1=lon_t(i,j+1)-lon_t(i,j-1)
!      if(x1<-pi)x1=x1+2.*pi
!      if(x1> pi)x1=x1-2.*pi
!      x0=-atan2(lat_t(i,j+1)-lat_t(i,j-1),x1*cos(lat_t(i,j)))+0.5*pi  !12-12-18

      enddo
      enddo

      do j=0,jmax+2   ! 1,NECO+1
      do i=0,imax+2   ! 1,MECO+1                                      !02/04/09

      x1=lon_t(i-1,j)                                                 !15-03-10
      if(x1-lon_t(i,j)<-pi)x1=x1+2.*pi
      if(x1-lon_t(i,j)> pi)x1=x1-2.*pi

        lat_u(i,j)=0.5*( lat_t(i  ,j  )+ lat_t(i-1,j  ))
!       lon_u(i,j)=0.5*( lon_t(i  ,j  )+ lon_t(i-1,j  ))
        lon_u(i,j)=0.5*( lon_t(i  ,j  )+ x1            )

      x2=lon_t(i,j-1)                                                 !15-03-10
      if(x2-lon_t(i,j)<-pi)x2=x2+2.*pi
      if(x2-lon_t(i,j)> pi)x2=x2-2.*pi

        lat_v(i,j)=0.5*( lat_t(i  ,j  )+ lat_t(i  ,j-1))
!       lon_v(i,j)=0.5*( lon_t(i  ,j  )+ lon_t(i  ,j-1))
        lon_v(i,j)=0.5*( lon_t(i  ,j  )+ x2            )

      x3=lon_t(i-1,j-1)                                               !15-03-10
      if(x3-lon_t(i,j)<-pi)x3=x3+2.*pi
      if(x3-lon_t(i,j)> pi)x3=x3-2.*pi

        lat_f(i,j)=                                                     &
       0.25*( lat_t(i,j)+ lat_t(i-1,j)+ lat_t(i,j-1)+ lat_t(i-1,j-1))
        lon_f(i,j)=                                                     &
!      0.25*( lon_t(i,j)+ lon_t(i-1,j)+ lon_t(i,j-1)+ lon_t(i-1,j-1))
       0.25*( lon_t(i,j)+ x1 + x2 +x3 )

      enddo
      enddo

#ifdef bidon
! Exemple de creation de tableaux globlon_t globlat_t ... !17-12-14
      if(.not.allocated(globlon_t))allocate(globlon_t(iglb  ,jglb  ))
      if(.not.allocated(globlat_t))allocate(globlat_t(iglb  ,jglb  ))
      if(.not.allocated(globlon_u))allocate(globlon_u(iglb+1,jglb  ))
      if(.not.allocated(globlat_u))allocate(globlat_u(iglb+1,jglb  ))
      if(.not.allocated(globlon_v))allocate(globlon_v(iglb  ,jglb+1))
      if(.not.allocated(globlat_v))allocate(globlat_v(iglb  ,jglb+1))
      do j=1,jmax ; do i=1,imax
       globlon_t(i+par%timax(1),j+par%tjmax(1))=lon_t(i,j)
       globlat_t(i+par%timax(1),j+par%tjmax(1))=lat_t(i,j)
      enddo       ; enddo
      ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
      ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
      call  par_gatherall_2d(globlon_t,lbound(globlon_t),ubound(globlon_t),ind,par%nbdom)
      call  par_gatherall_2d(globlat_t,lbound(globlat_t),ubound(globlat_t),ind,par%nbdom)

      do j=1,jmax ; do i=1,imax+1
       globlon_u(i+par%timax(1),j+par%tjmax(1))=lon_u(i,j)
       globlat_u(i+par%timax(1),j+par%tjmax(1))=lat_u(i,j)
      enddo       ; enddo
      ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax+1
      ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax
      call  par_gatherall_2d(globlon_u,lbound(globlon_u),ubound(globlon_u),ind,par%nbdom)
      call  par_gatherall_2d(globlat_u,lbound(globlat_u),ubound(globlat_u),ind,par%nbdom)

      do j=1,jmax+1 ; do i=1,imax
       globlon_v(i+par%timax(1),j+par%tjmax(1))=lon_v(i,j)
       globlat_v(i+par%timax(1),j+par%tjmax(1))=lat_v(i,j)
      enddo       ; enddo
      ind(1)=par%timax(1)+1 ;  ind(2)=par%timax(1)+imax
      ind(3)=par%tjmax(1)+1 ;  ind(4)=par%tjmax(1)+jmax+1
      call  par_gatherall_2d(globlon_v,lbound(globlon_v),ubound(globlon_v),ind,par%nbdom)
      call  par_gatherall_2d(globlat_v,lbound(globlat_v),ubound(globlat_v),ind,par%nbdom)

! Verification: tous les fichiers doivent etre identiques:
      do j=1,jglb+1 ; do i=1,iglb
       write(10+par%rank,*)i,j,globlon_v(i,j)*rad2deg &
                              ,globlat_v(i,j)*rad2deg
      enddo ; enddo
! Penser a deallocate!
      deallocate(globlon_t)
      deallocate(globlat_t)
      deallocate(globlon_u)
      deallocate(globlat_u)
      deallocate(globlon_v)
      deallocate(globlat_v)
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif
      stop 'Test longitudes latitudes globales termine'
#endif


      end subroutine grid_lonlat
!........................................................................
      subroutine grid_lonlat_from_file
      use pnetcdf
      implicit none
      integer ncid_
      integer(kind=MPI_OFFSET_KIND) start(4)
      integer(kind=MPI_OFFSET_KIND) edge(4)
#ifdef synopsis
       subroutinetitle='grid_lonlat_from_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      status=nfmpi_open(par%comm2d,initialgridfile_txt,nf_nowrite,MPI_INFO_NULL,ncid_)
      if(status/=0)stop 'erreur ouverture 2 fichier initialgridfile_txt'

! Lire les longitues:
      varstart(3)=1              ; varcount(3)=1        ! time
      varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2   ! j
      varstart(1)=1+par%timax(1) ; varcount(1)=imax+2   ! i
      start(1:4) = varstart(1:4)
      edge(1:4)  = varcount(1:4)

                   status=nfmpi_inq_varid(ncid_,'longitude_t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'glamt',var_id)
      if(status/=0)stop ' stop grid_lonlat var_id longitude'

      status=nfmpi_inq_vartype(ncid_,var_id,k0)
      if(k0==5) then  !----------------------------->
      status=nfmpi_get_vara_real_all(ncid_,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,0:jmax+1))
       lon_t(0:imax+1,0:jmax+1)=anyvar2d(0:imax+1,0:jmax+1)*deg2rad
      else            !----------------------------->
        stop 'grid_lonlat: cas longitude non real pas prevu'
      endif           !----------------------------->

! Lire les latitudes:
                   status=nfmpi_inq_varid(ncid_,'latitude_t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'gphit',var_id)
      if(status/=0)stop ' stop grid_lonlat var_id latitude'

      status=nfmpi_inq_vartype(ncid_,var_id,k0)
      if(k0==5) then  !----------------------------->
      status=nfmpi_get_vara_real_all(ncid_,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,0:jmax+1))
       lat_t(0:imax+1,0:jmax+1)=anyvar2d(0:imax+1,0:jmax+1)*deg2rad
      else            !----------------------------->
        stop 'grid_lonlat: cas latitude non real pas prevu'
      endif           !----------------------------->

      status=nfmpi_close(ncid_)

! Definir lont lat_ sur frontiere "z2"
      call obc_lonlat_z2 !17-04-15

      end subroutine grid_lonlat_from_file
!........................................................................
      subroutine grid_dim_from_file
      use pnetcdf
      implicit none
      integer ncid_
      integer(kind=MPI_OFFSET_KIND) start(4)
      integer(kind=MPI_OFFSET_KIND) edge(4)
      integer(kind=MPI_OFFSET_KIND) :: imax_,jmax_,kmax_
#ifdef synopsis
       subroutinetitle='grid_dim_from_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,'(a,a)')'Ouverture fichier ',trim(initialgridfile_txt)
      status=nfmpi_open(par%comm2d,trim(initialgridfile_txt),nf_nowrite,MPI_INFO_NULL,ncid_)
      if(status/=0) then !17-07-14
        write(6,'(a,a)')'Error 616 nfmpi_open file:',trim(initialgridfile_txt)
        stop 'erreur ouverture 3 fichier initialgridfile_txt'
      endif

                    status=nfmpi_inq_dimid(ncid_   ,'x_zhl',dim_x_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'x_ZHL',dim_x_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'imax_t',dim_x_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'ni_t',dim_x_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'ni',dim_x_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'x',dim_x_id)
       if(status/=0)stop 'grid_dxdy erreur dim_x_id'

                    status=nfmpi_inq_dimid(ncid_   ,'y_zhl',dim_y_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'y_ZHL',dim_y_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'jmax_t',dim_y_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'nj_t',dim_y_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'nj',dim_y_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'y',dim_y_id)
       if(status/=0)stop 'grid_dxdy erreur dim_y_id'

                    status=nfmpi_inq_dimid(ncid_   ,'z_zhl',dim_z_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'z_ZHL',dim_z_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'kmax_t',dim_z_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'nk_t',dim_z_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'nk',dim_z_id)
       if(status/=0)status=nfmpi_inq_dimid(ncid_   ,'z',dim_z_id)
       if(status/=0)stop 'grid_dxdy erreur dim_z_id'

      status=nfmpi_inq_dimlen(ncid_,dim_x_id,imax_)
      status=nfmpi_inq_dimlen(ncid_,dim_y_id,jmax_)
      status=nfmpi_inq_dimlen(ncid_,dim_z_id,kmax_)

      if(imax_/=iglb+2.or.jmax_/=jglb+2.or.kmax_/=kmax)then
      write(6,*)'Le fichier de grille ne correspond pas aux dimensions'
      write(6,*)'donnees dans notebook_grid'
      write(6,*)'imax_/=imax+2.or.jmax_/=jmax+2.or.kmax_/=kmax'
      write(6,*)imax_,imax+2
      write(6,*)jmax_,jmax+2
      write(6,*)kmax_,kmax
      stop 'grid_dxdy_from_file'
      endif

      status=nfmpi_close(ncid_)
      end subroutine grid_dim_from_file
!........................................................................
      subroutine grid_dxdy_from_file
      use pnetcdf
      implicit none
      integer ncid_
      integer(kind=MPI_OFFSET_KIND) start(4)
      integer(kind=MPI_OFFSET_KIND) edge(4)
#ifdef synopsis
       subroutinetitle='grid_dxdy_from_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,'(a,a)')'Ouverture fichier ',trim(initialgridfile_txt)
      status=nfmpi_open(par%comm2d,initialgridfile_txt,nf_nowrite,MPI_INFO_NULL,ncid_)
      if(status/=0)stop 'erreur ouverture 1 fichier initialgridfile_txt'

! Lire dx_t:
      varstart(3)=1              ; varcount(3)=1        ! time
      varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2   ! j
      varstart(1)=1+par%timax(1) ; varcount(1)=imax+2   ! i
      start(1:4) = varstart(1:4)
      edge(1:4)  = varcount(1:4)

                   status=nfmpi_inq_varid(ncid_,'dx_t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'e1t',var_id)
      if(status/=0)stop ' stop grid_dxdy var_id e1t'

      status=nfmpi_inq_vartype(ncid_,var_id,k0)
      if(k0==6) then  !----------------------------->
      status=nfmpi_get_vara_double_all(ncid_,var_id,start(1:3),edge(1:3) &
                                          ,dx_t(0:imax+1,0:jmax+1))
      else            !----------------------------->
        write(6,*)'k0=',k0
        stop 'grid_dxdy: cas e1t real pas prevu'
      endif           !----------------------------->

! Lire dy_t:
                   status=nfmpi_inq_varid(ncid_,'dy_t',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid_,'e2t',var_id)
      if(status/=0)stop ' stop grid_dxdy var_id latitude'

      status=nfmpi_inq_vartype(ncid_,var_id,k0)
      if(k0==6) then  !----------------------------->
      status=nfmpi_get_vara_double_all(ncid_,var_id,start(1:3),edge(1:3) &
                                          ,dy_t(0:imax+1,0:jmax+1))
      else            !----------------------------->
        write(6,*)'k0=',k0
        stop 'grid_dxdy: cas e2t real pas prevu'
      endif           !----------------------------->

!     do j=0,jmax+1
!     do i=0,imax+1
!      write(67,*)dx_t(i,j),dy_t(i,j)
!     enddo
!     enddo

      status=nfmpi_close(ncid_)
      end subroutine grid_dxdy_from_file

!..............................................................................................

      subroutine grid_bipole_transformation
! Transformation conforme "2 poles" de Bentsen et al, 1999.
! Ref:
! "Coordinate Transformation on a Sphere Using Conformal Mapping"
!  M. BENTSEN, G. EVENSEN, H. DRANGE, AND A. D. JENKINS
!  MWR, 1999, 2733-2740
      implicit none
!      double precision xa_,ya_,xb_,yb_,xc_,yc_,mu_,psi_,teta_,phi_   &
      real*16 ::       xa_,ya_,xb_,yb_,xc_,yc_,mu_,psi_,teta_,phi_   &
                      ,x_,y_,modzz_                                  &
                      ,xlon_pn_,ylat_pn_                             &
                      ,xlon_ps_,ylat_ps_                             &
                      ,cx_,cy_,cz_                                   &
                      ,xlon_c_,ylat_c_
!      complex ww_,zz_,za_,zb_,zc_
      complex*32 ::  ww_,zz_,za_,zb_,zc_
#ifdef synopsis
       subroutinetitle='grid_bipole_transformation'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Faire la formu_le 10 de l'article de Bentsen et al, 1999 (B99)
! Autrement faire les positions des points a b c des figures 1 2 3 4
! Attention dans B99, phi_=0 au pole Nord et phi_=pi au pole Sud
! Par consequent, quand on entre les positions des deux poles,
! on pense e faire lat=(90-lat)*pi/180 puisque logiquement
! on va specifier les positions en lon lat classiques...

! BOUCLE BIDON
!     i10=50
!     do i1=-i10,i10
!     write(6,*)i1
!     do j1=0,i10
!     northpole_lon=180.*real(i1)/real(i10)
!     northpole_lat=180.*real(j1)/real(i10)
!     do i2=-i10,i10
!     do j2=0,i10
!     ksecu=1
!     if(j1==0.and.j2==0)ksecu=0
!     if(j1==i10.and.j2==i10)ksecu=0
!     if(i2==i1.and.j2==j1)ksecu=0
!     if(ksecu==1) then !--------->
!     southpole_lon=180.*real(i2)/real(i10)
!     southpole_lat=180.*real(j2)/real(i10)

! Cas du pole Sud aux antipodes (c.a.d. reset a 9999 dans notebook_grid)
      if(southpole_lon==9999.or.southpole_lat==9999.) then !---->
         stop 'NE PAS PASSER PAR ROUTINE grid_bipole_transformation'
!      southpole_lon= northpole_lon+180.
!      if(southpole_lon>180.)southpole_lon=southpole_lon-360.
!      southpole_lat=-northpole_lat
      endif                                                !---->

! Point "a" position du nouveau "Pole Nord":
!     xlon_pn_=-5.6 ; ylat_pn_=36.5
!     xlon_pn_=30. ; ylat_pn_=45.
      xlon_pn_=northpole_lon ; ylat_pn_=northpole_lat
! Point "b" position nouveau "Pole Sud":
!     xlon_ps_=-5.6 ; ylat_ps_=35.5
!     xlon_ps_=-100 ; ylat_ps_=40.
      xlon_ps_=southpole_lon ; ylat_ps_=southpole_lat

      xlon_pn_=xlon_pn_*deg2rad16
      xlon_ps_=xlon_ps_*deg2rad16
      ylat_pn_=(90.-ylat_pn_)*deg2rad16
      ylat_ps_=(90.-ylat_ps_)*deg2rad16

! Le point "c" est un point median sur la courbe geodesique "ab" donne par les formu_les (11,12) de B99
! Eq 11:
      cx_=cos(xlon_pn_)*sin(ylat_pn_) + cos(xlon_ps_)*sin(ylat_ps_)
      cy_=sin(xlon_pn_)*sin(ylat_pn_) + sin(xlon_ps_)*sin(ylat_ps_)
      cz_=cos(ylat_pn_) + cos(ylat_ps_)
! Eq 12:
!     xlon_c_=atan(cy_/cx_)
!     ylat_c_=atan(sqrt(cx_*cx_+cy_*cy_)/cz_)
      xlon_c_=atan2(cy_,cx_)
      ylat_c_=atan2(sqrt(cx_*cx_+cy_*cy_),cz_)

!     write(6,*)'coucou',xlon_c_*rad2deg,90.-ylat_c_*rad2deg
!     write(6,*)'cxcy cz=',sqrt(cx_*cx_+cy_*cy_),cz_
!     call lonlat2distance(xlon_c_,0.5*pi-ylat_c_,xlon_pn_,0.5*pi-ylat_pn_,x1)
!     write(6,*)'dist1=',x1/1000.
!     call lonlat2distance(xlon_c_,0.5*pi-ylat_c_,xlon_ps_,0.5*pi-ylat_ps_,x2)
!     write(6,*)'dist2=',x2/1000.
!     write(65,'(6(e13.6,1x))')                          &
!       northpole_lon,northpole_lat                      &
!      ,southpole_lon,southpole_lat,x1,x2

!     if(abs(x2-x1)>1.)write(66,'(6(e13.6,1x))')        &
!       northpole_lon,northpole_lat                      &
!      ,southpole_lon,southpole_lat,x1,x2

!     endif !---------->

!     enddo
!     enddo
!     enddo
!     enddo
!     stop


! On applique e ces points la transformation z=f(lon,lat) de B99, z etant un nombre complexe
      xa_=tan(0.5*ylat_pn_)*cos(xlon_pn_)
      ya_=tan(0.5*ylat_pn_)*sin(xlon_pn_)
      xb_=tan(0.5*ylat_ps_)*cos(xlon_ps_)
      yb_=tan(0.5*ylat_ps_)*sin(xlon_ps_)
      xc_=tan(0.5*ylat_c_)*cos(xlon_c_)
      yc_=tan(0.5*ylat_c_)*sin(xlon_c_)
      za_=cmplx(xa_,ya_)
      zb_=cmplx(xb_,yb_)
      zc_=cmplx(xc_,yc_)
! Et a ce stade on a donc les coef a,b,c qui servent aux transformations Eq6 et Eq7 de B99

! Exemple de sequence oe l'on part des "fausses" lon lat mesurees dans l'univers tordu
! pour arriver aux "vraies" lon lat sur la terre
! Pour le moment on fait simple, mais e terme ces lon lat seront regulierement distribuees faeon "mercator"
! On va prendre un "cercle" de "fausse latitude" constante dans l'univers tordu et voir ce que ce cerle devient
! sur la terre

! Coordonnes (mu_,psi_) de ce cercle (attention aux conventions de latitude: 0=PN, pi=PS):
!     do j=1,100
!     do i=- 90, 800,10
      do j=-1,jmax+2
      do i=-1,imax+2                                       !02/04/09

! Point de depart: lon lat de l'univers tordu
!      psi_=(0.1+0.8*real(j-1)/99.)*pi
!      mu_=pi*real(i)/1000.   ! longitude univers tordu
       mu_=lon_t(i,j)
       psi_=0.5*pi16-lat_t(i,j)

! Passage de (mu_,psi_) au nombre complexe ww_ de B99 (Eq3):
       x_=tan(0.5*psi_)*cos(mu_)
       y_=tan(0.5*psi_)*sin(mu_)
       ww_=cmplx(x_,y_)

! Passage de ww_ e zz_ (Eq7 de B99):
       zz_= (-zb_*ww_*(zc_-za_) + za_*(zc_-zb_))/(-ww_*(zc_-za_)+zc_-zb_)

! Passage de zz_ e (teta_,phi_) vraies lon lat sur la terre (Eq2)
       modzz_=sqrt(real(zz_)**2.+aimag(zz_)**2)
       teta_=2.*atan(aimag(zz_)/(real(zz_)+modzz_))
       phi_=2.*atan(modzz_)

       lon_t(i,j)=teta_
       lat_t(i,j)=0.5*pi16-phi_

      enddo
      enddo

      end subroutine grid_bipole_transformation

!...................................................................

      subroutine grid_monopole_transformation
      implicit none
#ifdef synopsis
       subroutinetitle='grid_monopole_transformation'
       subroutinedescription= &
       'Grid pole location possibly different from Earth North Pole'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      const1=(90.-northpole_lat)*deg2rad
      const2=( 0.-northpole_lon)*deg2rad
      do j=-1,jmax+2
      do i=-1,imax+2

        x2=rayonterre*cos(lon_t(i,j))*cos(lat_t(i,j)) !15-03-10
        y2=rayonterre*sin(lon_t(i,j))*cos(lat_t(i,j)) !15-03-10
        z2=rayonterre*sin(lat_t(i,j))

        x1= cos(const1)*x2+sin(const1)*z2
        z1=-sin(const1)*x2+cos(const1)*z2
        y1=y2

        x0= cos(const2)*x1+sin(const2)*y1
        y0=-sin(const2)*x1+cos(const2)*y1


       lon_t(i,j)=atan2(y0,x0)                         !15-03-10
       lat_t(i,j)=asin(z1/rayonterre)

      enddo
      enddo

      end subroutine grid_monopole_transformation

!..............................................................................................

      subroutine grid_lonlat2ij(lon1_,lat1_,txt_)
      implicit none
      double precision lon1_,lat1_,dlon_di_,dlon_dj_,dlon_dm_
      character*3 txt_
      integer :: zone_
#ifdef synopsis
       subroutinetitle='grid_lonlat2ij'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! CommentE le 04-02-20 car il faut construire le tableau bien avant pour eviter les
! accident de lecture de fichier restart (voir dans set_parameters)
!     if(.not.allocated(lonlat2ij_t))call grid_lonlat2ij_initial !04-02-20

       range_lon=nint(1.+(rad2deg*lon1_-lonmin)*(range_lon_max-1)/(lonmax-lonmin))
       range_lat=nint(1.+(rad2deg*lat1_-latmin)*(range_lat_max-1)/(latmax-latmin))

       range_lon=min(max(range_lon,1),range_lon_max)
       range_lat=min(max(range_lat,1),range_lat_max)

!      write(6,*)'lon lat=',lon1_*rad2deg,lat1_*rad2deg
!      write(6,*)'range_lon=',range_lon
!      write(6,*)'range_lat=',range_lat
!      write(6,*)'lon lat',                        &
!                 (lonmin+real(range_lon    -1)    &
!                        /real(range_lon_max-1)    &
!                        *(lonmax-lonmin))         &
!                ,(latmin+real(range_lat    -1)    &
!                        /real(range_lat_max-1)    &
!                        *(latmax-latmin))

      deci=lonlat2ij_t(range_lon,range_lat,1)
      decj=lonlat2ij_t(range_lon,range_lat,2)

      const1=sin(lat1_) ; const2=cos(lat1_)

      zone_=0
      loop2=40
      x10=20.
      do loop1=1,loop2

      if(deci<1.    .and.par%tvoisin(ouest)==par%rank)     deci=imax-1.
      if(deci>imax  .and.par%tvoisin(est  )==par%rank)     deci=2.
      if(decj>jmax+1) then !jjjjjjjjj>
!      if(par%tvoisin(nord )==mpi_proc_null) then !mmmm>
        deci=0.125*imax+mod(zone_,4)*0.25*imax
        decj=0.75*jmax
        zone_=zone_+1
!       write(6,*)'zone_=',zone_
!      endif                                      !mmmm>
      endif                !jjjjjjjjj>

      deci=min(max(deci,0.001d00),imax+0.999d00)
      decj=min(max(decj,0.001d00),jmax+0.999d00)

!     write(6,*)'deci decj=',loop1,deci,decj

      i1=int(deci) ; j1=int(decj)
      rapi=deci-i1 ; rapj=decj-j1

      do j2=0,1 ; do i2=0,1

      i0=i1+i2 ; j0=j1+j2

      dlon_di_=lon_t(i0+1,j0  )-lon_t(i0-1,j0  )
      dlon_dj_=lon_t(i0  ,j0+1)-lon_t(i0  ,j0-1)
      dlon_dm_=lon1_           -lon_t(i0  ,j0  )

      if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
      if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
      if(dlon_dj_<-pi)dlon_dj_=dlon_dj_+2.*pi
      if(dlon_dj_> pi)dlon_dj_=dlon_dj_-2.*pi
      if(dlon_dm_<-pi)dlon_dm_=dlon_dm_+2.*pi
      if(dlon_dm_> pi)dlon_dm_=dlon_dm_-2.*pi

! Determinant principal:
      x1=( dlon_di_                                       &
          *( lat_t(i0  ,j0+1)- lat_t(i0  ,j0-1))          &
          -( lat_t(i0+1,j0  )- lat_t(i0-1,j0  ))          &
          *dlon_dj_ )*0.25

      ir8(i2,j2)=                                         &
       i0+max(min(( dlon_dm_                              &
          *( lat_t(i0  ,j0+1)- lat_t(i0,j0-1))            &
          -( lat1_           - lat_t(i0,j0))              &
          *dlon_dj_)/x1*0.5,x10),-x10)

      jr8(i2,j2)=                                         &
       j0+max(min(( dlon_di_                              &
          *( lat1_           - lat_t(i0  ,j0))            &
          -( lat_t(i0+1,j0  )- lat_t(i0-1,j0))            &
          *dlon_dm_   )/x1*0.5,x10),-x10)

!     write(6,*)'dx dy=',dx_t(i0,j0),dy_t(i0,j0)

      enddo ; enddo

      x2=deci ; x3=decj

      deci=(1.-rapi)*(1.-rapj)*ir8(0,0)   &
          +(1.-rapi)*    rapj *ir8(0,1)   &
          +    rapi *    rapj *ir8(1,1)   &
          +    rapi *(1.-rapj)*ir8(1,0)
      decj=(1.-rapi)*(1.-rapj)*jr8(0,0)   &
          +(1.-rapi)*    rapj *jr8(0,1)   &
          +    rapi *    rapj *jr8(1,1)   &
          +    rapi *(1.-rapj)*jr8(1,0)

!      write(6,*)'deci decj=',deci,decj
!      write(66,*)deci,decj
!      write(6,*)'lon lat correspondant='     &
!               ,lon_t(nint(deci),nint(decj))*rad2deg &
!               ,lat_t(nint(deci),nint(decj))*rad2deg

!     i1=nint(deci) ; j1=nint(decj)


       if(abs(x2-deci)+abs(x3-decj)<0.01)goto 963
      enddo

  963 continue

!     write(6,*)'AVANT',par%rank,deci,decj
!     if(par%rank==7)write(6,*)'deci decj=',deci,decj
!     if(par%rank==7)write(6,*)'loop1 loop2=',loop1,loop2
!     write(6,*)'deci decj=',deci,decj

      if(deci<0.or.deci>imax+1.or. &
         decj<0.or.decj>jmax+1.or.loop1>loop2)then             !----->
! Hors proc, deci=decj=0 et le flag sum0=0.
       deci=0. ; decj=0. ; sum0=0.
      else                                                       !----->
       deci=deci+par%timax(1) ; decj=decj+par%tjmax(1) ; sum0=1.
      endif                                                      !----->

!............................................................................
#ifdef parallele
      if(txt_=='glb') then !>>>>>>>>>>>>>>>>>
! On envoie la position e tous les autres proc. On utilise la fonction somme
! de sorte que si un point est e cheval sur plusieurs proc, on prendra la position
! moyenne. En dehors du proc deci=decj=0 ne modifie pas la somme....
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(deci,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(decj,sum2glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum0=sum0glb
      if(sum0/=0.) then !----->
        deci=sum1glb/sum0 ; decj=sum2glb/sum0
      endif             !----->
      endif                !>>>>>>>>>>>>>>>>>
#endif
!............................................................................

      if(sum0==0.) then !>>>>>
       deci=-999.  ; decj=-999.
      endif             !>>>>>

!     write(6,*)'APRES',par%rank,deci,decj

      end subroutine grid_lonlat2ij

!..............................................................................................

      subroutine grid_lonlat2ij_initial
      implicit none
#ifdef synopsis
       subroutinetitle='grid_lonlat2ij_initial'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Attribuer des valeurs d'indices de grilles e quelques valeurs de
! longitude et de latitude

!     write(6,*)'lon min et max=',lonmin,lonmax
!     write(6,*)'lat min et max=',latmin,latmax
      allocate(lonlat2ij_t(range_lon_max,range_lat_max,2))
      allocate(ir8(0:1,0:1))
      allocate(jr8(0:1,0:1))

      do range_lat=1,range_lat_max
      do range_lon=1,range_lon_max

! x1= gamme de longitude en radians
! x2= gamme de latitude en radians
       x1=deg2rad*(lonmin+real(range_lon    -1)    &
                         /real(range_lon_max-1)    &
                         *(lonmax-lonmin))
       x2=deg2rad*(latmin+real(range_lat    -1)    &
                         /real(range_lat_max-1)    &
                         *(latmax-latmin))

!      write(6,*)'x1 x2=',x1*rad2deg,x2*rad2deg

       lonlat2ij_t(range_lon,range_lat,1)=-999
       lonlat2ij_t(range_lon,range_lat,2)=-999
       i1=-999 ; j1=-999 ; dist1=1.d10
       do j=1,jmax
       do i=1,imax

        dist2=acos( min(max( sin(x2)*sin(lat_t(i,j))                   &
                            +cos(x2)*cos(lat_t(i,j))*cos(lon_t(i,j)-x1),-1.),1.)) !22-01-19

        if(dist2<dist1) then !>>>>
         i1=i ; j1=j ; dist1=dist2
        endif                !>>>>

       enddo
       enddo
       if(i1==-999)stop ' STOP i1=-999 dans grid_lonlat2ij_initial'
!      if(dist1<sqrt(dx_t(i1,j1)**2+dy_t(i1,j1)**2)) then !---->
        lonlat2ij_t(range_lon,range_lat,1)=i1
        lonlat2ij_t(range_lon,range_lat,2)=j1
!      endif                                              !---->

!     write(66,*)x1*rad2deg,x2*rad2deg
!     write(67,*)lon_t(i1,j1)*rad2deg,lat_t(i1,j1)*rad2deg


      enddo
      enddo

      end subroutine grid_lonlat2ij_initial

!..............................................................................................

      subroutine grid_bipole_lonlat2ij(lon1_8,lat1_8)
! Transformation conforme "2 poles" de Bentsen et al, 1999.
! Ref:
! "Coordinate Transformation on a Sphere Using Conformal Mapping"
!  M. BENTSEN, G. EVENSEN, H. DRANGE, AND A. D. JENKINS
!  MWR, 1999, 2733-2740
      implicit none
      double precision :: lon1_8,lat1_8
!      double precision xa_,ya_,xb_,yb_,xc_,yc_,mu_,psi_,teta_,phi_   &
      real*16 ::       xa_,ya_,xb_,yb_,xc_,yc_,mu_,psi_,teta_,phi_   &
                      ,x_,y_,modww_                                  &
                      ,xlon_pn_,ylat_pn_                             &
                      ,xlon_ps_,ylat_ps_                             &
                      ,cx_,cy_,cz_                                   &
                      ,xlon_c_,ylat_c_                               &
                      ,lon1_,lat1_                                   &
                      ,lon2_,lat2_
!      complex ww_,zz_,za_,zb_,zc_
      complex*32 ww_,zz_,za_,zb_,zc_

#ifdef synopsis
       subroutinetitle='grid_bipole_lonlat2ij'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      lon1_=lon1_8
      lat1_=lat1_8

! Faire la formu_le 10 de l'article de Bentsen et al, 1999 (B99)
! Autrement faire les positions des points a b c des figures 1 2 3 4
! Attention dans B99, phi_=0 au pole Nord et phi_=pi au pole Sud
! Par consequent, quand on entre les positions des deux poles,
! on pense e faire lat=(90-lat)*pi/180 puisque logiquement
! on va specifier les positions en lon lat classiques...

! Point "a" position du nouveau "Pole Nord":
      xlon_pn_=northpole_lon ; ylat_pn_=northpole_lat
! Point "b" position nouveau "Pole Sud":
      xlon_ps_=southpole_lon ; ylat_ps_=southpole_lat

      xlon_pn_=xlon_pn_*deg2rad16
      xlon_ps_=xlon_ps_*deg2rad16
      ylat_pn_=(90.-ylat_pn_)*deg2rad16
      ylat_ps_=(90.-ylat_ps_)*deg2rad16

! Le point "c" est un point median sur la courbe geodesique "ab" donne par les formu_les (11,12) de B99
! Eq 11:
      cx_=cos(xlon_pn_)*sin(ylat_pn_) + cos(xlon_ps_)*sin(ylat_ps_)
      cy_=sin(xlon_pn_)*sin(ylat_pn_) + sin(xlon_ps_)*sin(ylat_ps_)
      cz_=cos(ylat_pn_) + cos(ylat_ps_)
! Eq 12:
      xlon_c_=atan2(cy_,cx_)
      ylat_c_=atan2(sqrt(cx_*cx_+cy_*cy_),cz_)

! On applique e ces points la transformation z=f(lon,lat) de B99, z etant un nombre complexe
      xa_=tan(0.5*ylat_pn_)*cos(xlon_pn_)
      ya_=tan(0.5*ylat_pn_)*sin(xlon_pn_)
      xb_=tan(0.5*ylat_ps_)*cos(xlon_ps_)
      yb_=tan(0.5*ylat_ps_)*sin(xlon_ps_)
      xc_=tan(0.5*ylat_c_)*cos(xlon_c_)
      yc_=tan(0.5*ylat_c_)*sin(xlon_c_)
      za_=cmplx(xa_,ya_)
      zb_=cmplx(xb_,yb_)
      zc_=cmplx(xc_,yc_)
! Et a ce stade on a donc les coef a,b,c qui servent aux transformations Eq6 et Eq7 de B99

! Point de depart: lon lat de l'univers tordu
       teta_=lon1_
       phi_=0.5*pi16-lat1_

! Passage de (teta_,phi_) au nombre complexe zz_ de B99:
       x_=tan(0.5*phi_)*cos(teta_)
       y_=tan(0.5*phi_)*sin(teta_)
       zz_=cmplx(x_,y_)


! Passage de zz_ e ww_ :
       ww_= (zz_-za_)*(zc_-zb_)/((zz_-zb_)*(zc_-za_))

! Passage de ww_ e (mu_,psi_):
       modww_=sqrt(real(ww_)**2.+aimag(ww_)**2)
       mu_=2.*atan(aimag(ww_)/(real(ww_)+modww_))
       psi_=2.*atan(modww_)

       lon2_=mu_
       lat2_=0.5*pi16-psi_

       deci=                                                 &
         (rayonterre*cos(phi0)*lon2_                         &
         -gridcard(8,1) )/dxb                                &
         +gridcard(5,1)

       decj=                                                          &
         (rayonterre*cos(phi0)*log( (1.+sin(lat2_  ))/cos(lat2_  ) )  &
         -gridcard(9,1) )/dyb                                         &
         +gridcard(6,1)

! La fonction atan2 renvoie un angle lon2_ compris entre -pi et +pi ce qui !10-09-13
! induit un possibilite d'erreur sur deci qu'on leve en calculant
! la valeur de deci correspondant a lon2+2pi et en choississant la
! plus probable des 2 valeurs:
! Le test if(lon2_<0) a EtE commentE le !27-03-22 suite A faille constateE dans 
! bathymaker (voir notebook_grid_ile_de_re_50m.f)
!      if(lon2_<0.) then !nnnnnnnnnnnnn>
         x3=deci+rayonterre*cos(phi0)*2.*pi/dxb
         if(abs(x3-0.5*iglb)<abs(deci-0.5*iglb))deci=x3
!      endif             !nnnnnnnnnnnnn>

      end subroutine grid_bipole_lonlat2ij

!...................................................................

      subroutine grid_check_lonlat2ij
      implicit none
      integer i_,j_,flag_
#ifdef synopsis
       subroutinetitle='grid_check_lonlat2ij'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do j_=3,jmax-2   ! boucles evitant les obc pour que le test reste
      do i_=3,imax-2   ! compatible avec conditions periodiques

!Grille horizontale generalisee eventuellement chaotique dans le masque terre
! par consequent verification sur points de mer uniquement
      if(mask_t(i_,j_,kmaxp1)==1) then !11111111>

       call latlontoij(lon_t(i_,j_),lat_t(i_,j_),'loc')

       flag_=1

       if(abs(deci-i_)>0.01) then
        flag_=0
!       if(iperiodicboundary) then
!        if(nint(abs(deci-i_))==imax-2)flag_=1
!       endif
       endif

       if(abs(decj-j_)>0.01) then
         flag_=0
!       if(jperiodicboundary) then
!        if(nint(abs(decj-j_))==jmax-2)flag_=1
!       endif
       endif

       if(flag_==0) then !------>
        write(6,*)'Fonction latlontoij ne renvoie pas bons deci decj'
        write(6,*)'local imax jmax par%rank',imax,jmax,par%rank
        write(6,*)'local i j par%rank      ',i_,j_,par%rank
        write(6,*)'local deci decj par%rank',deci,decj,par%rank
        write(6,*)'iglb jglb               ',iglb,jglb
        write(6,*)'Global i j              ',i_+par%timax(1),j_+par%tjmax(1)
!       call mpi_finalize(ierr)
        stop 'STOP in subroutine grid_check_lonlat2ij'
       endif             !------>

      endif                          !11111111>

      enddo
      enddo

      end subroutine grid_check_lonlat2ij

!...................................................................

      subroutine grid_mpi_obc !07-12-13
      implicit none
      integer loop_,id_dx_t_za,id_dy_t_za
#ifdef synopsis
       subroutinetitle='grid_mpi_obc'
       subroutinedescription= &
       'Lateral boundary condition on dx_t dy_t ensuring mpi continuity'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
      call get_type_echange('za','dx_t_za',dx_t    &
                                   ,lbound(dx_t)   &
                                   ,ubound(dx_t)   &
                                       ,id_dx_t_za)
      call get_type_echange('za','dy_t_za',dy_t    &
                                   ,lbound(dy_t)   &
                                   ,ubound(dy_t)   &
                                       ,id_dy_t_za)
! Exchanges: !22-06-14
      do loop_=1,subcycle_exchange
        call echange_voisin(dx_t,id_dx_t_za,mpi_neighbor_list(loop_))
        call echange_voisin(dy_t,id_dy_t_za,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine grid_mpi_obc

!...................................................................

      subroutine grid_mask_h_offlinecase
      implicit none
      integer ncid_
#ifdef synopsis
       subroutinetitle='grid_mask_h_offlinecase'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'Initial mask_t h_w kmin_w'

      texte90=trim(lonlatfile) !07-02-17
      if(par%rank==0)write(6,'(a,a)')'OPEN ',trim(texte90)
      status=nf_open(trim(texte90),nf_nowrite,ncid_)
      if(status/=0)stop 'Err 1553 grid_mask_h_offlinecase nf_open'

!.................................................
! Verifier les dimensions du fichier d'entree ! 07-02-16
                   status=nf_inq_dimid(ncid_   ,'ni_t',dim_x_id)
      if(status/=0)stop 'Err 558 grid_mask_h_offlinecase nf_inq_dimid'
                   status=nf_inq_dimid(ncid_   ,'nj_t',dim_y_id)
      if(status/=0)stop 'Err 559 grid_mask_h_offlinecase nf_inq_dimid'
                   status=nf_inq_dimid(ncid_   ,'nk_t',dim_z_id)
      if(status/=0)stop 'Err 560 grid_mask_h_offlinecase nf_inq_dimid'

                   status=nf_inq_dimlen(ncid_,dim_x_id,i0)
      if(status/=0)stop 'Err 1 grid_mask_h_offlinecase nf_inq_dimlen'
                   status=nf_inq_dimlen(ncid_,dim_y_id,j0)
      if(status/=0)stop 'Err 2 grid_mask_h_offlinecase nf_inq_dimlen'
                   status=nf_inq_dimlen(ncid_,dim_z_id,k0)
      if(status/=0)stop 'Err 3 grid_mask_h_offlinecase nf_inq_dimlen'

      if(i0/=iglb+2.or.j0/=jglb+2.or.k0/=kmax) then !pmxpmx>
        if(par%rank==0) then !ooo>
         write(6,*)'--------------'
         write(6,'(a)')'Grid dimensions are inconsistent'
         write(6,'(a,a)')'In ',trim(nomfichier(2))
         write(6,*)iglb,jglb,kmax
         write(6,'(a,a)')'In ',trim(lonlatfile)
         write(6,*)i0-2,j0-2,k0
         write(6,*)'--------------'
        endif                !ooo>
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif
        stop ' Grid size insconsistency in grid_mask_h_offlinecase'
      endif                                         !pmxpmx>
!.................................................

! Deduire flag_merged_levels de l'attribut general Vertical_coordinate du fichier grid.nc 14-07-21
       status=nf_get_att_text(ncid_,nf_global,'Vertical_coordinate',texte30)
       if(status/=0)stop ' stop nf_get_att_text Vertical_coordinate'
       if(trim(texte30)=='s-z') then !ooo> !14-07-21
        flag_merged_levels=1
       else                          !ooo>
        flag_merged_levels=0
       endif                         !ooo>



       status=nf_inq_varid(ncid_,'mask_t',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid mask_t grid_mask_h_offlinecase'

       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_short)stop &
       ' type mask non prevu dans grid_mask_h_offlinecase'

! On lit mask_t(3D) mais il faut aussi glob_mask (2D) dans la subroutine
! set_rivers(case_=1) pour testet les positions des points d'entree des
! fleuves
! Lire glob_mask(2D):
       varstart(1)=1    ; varcount(1)=iglb+2
       varstart(2)=1    ; varcount(2)=jglb+2
!      varstart(3)=kmax ; varcount(3)=1
       status=nf_get_vara_int1(ncid_,var_id,varstart(1:2)   & !23-11-22
                                           ,varcount(1:2)   & !23-11-22
                            , glob_mask(0:iglb+1,0:jglb+1))
       if(status/=0)stop ' Erreur 1315 nf_get_vara_int'

! Deduire mask_t de globmask !14-07-21
      do j=0,jmax+1 ; do i=0,imax+1
        mask_t(i,j,:)=glob_mask(i+par%timax(1),j+par%tjmax(1))
      enddo ; enddo

! Lire kmin_w: !07-01-15
       status=nf_inq_varid(ncid_,'kmin_w',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid kmin_w grid_mask_h_offlinecase'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_short)stop &
       ' type kmin_w non prevu dans grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2
       status=nf_get_vara_int(ncid_,var_id,varstart(1:2)   &
                                          ,varcount(1:2)   &
                              ,kmin_w(0:imax+1,0:jmax+1))
       if(status/=0)stop ' Erreur 1545 nf_get_vara_int'

! Lire kmin_u: !09-01-15
       status=nf_inq_varid(ncid_,'kmin_u',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid kmin_u grid_mask_h_offlinecase'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_short)stop &
       ' type kmin_u non prevu dans grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+1
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2
       status=nf_get_vara_int(ncid_,var_id,varstart(1:2)   &
                                          ,varcount(1:2)   &
                              ,kmin_u(1:imax+1,0:jmax+1))
       if(status/=0)stop ' Erreur nf_get_vara_int kmin_u'


       if(flag_merged_levels==1) then !m°v°m> !14-07-21
! Lire kmergedr4_u
       status=nf_inq_varid(ncid_,'kmergedr4_u',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid kmergedr4_u grid_mask_h_offlinecase'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_real)stop &
       ' type kmergedr4_u non prevu dans grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+1
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2
       status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                           ,varcount(1:2)   &
                              ,kmergedr4_u(1:imax+1,0:jmax+1))
       if(status/=0)stop ' Erreur nf_get_vara_int kmergedr4_u'
       endif                          !m°v°m> !14-07-21

! Lire kmin_v: !09-01-15
       status=nf_inq_varid(ncid_,'kmin_v',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid kmin_v grid_mask_h_offlinecase'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_short)stop &
       ' type kmin_v non prevu dans grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+1
       status=nf_get_vara_int(ncid_,var_id,varstart(1:2)   &
                                          ,varcount(1:2)   &
                              ,kmin_v(0:imax+1,1:jmax+1))
       if(status/=0)stop ' Erreur nf_get_vara_int kmin_v'

       if(flag_merged_levels==1) then !m°v°m> !14-07-21
! Lire kmergedr4_v
       status=nf_inq_varid(ncid_,'kmergedr4_v',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid kmergedr4_v grid_mask_h_offlinecase'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_real)stop &
       ' type kmergedr4_v non prevu dans grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+1
       status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                           ,varcount(1:2)   &
                              ,kmergedr4_v(0:imax+1,1:jmax+1))
       if(status/=0)stop ' Erreur nf_get_vara_int kmergedr4_v'
       endif                          !m°v°m> !14-07-21

! Lire h_w
       status=nf_inq_varid(ncid_,'h_w',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid h_w grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2 !22-03-15
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2 !22-03-15

       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0==nf_double)  &
       status=nf_get_vara_double(ncid_,var_id,varstart(1:2)   &
                                             ,varcount(1:2)   &
                               ,anyv3d(0:imax+1,0:jmax+1,1,1))
!                                   ,h_w(0:imax+1,0:jmax+1))
           h_w(0:imax+1,0:jmax+1)=anyv3d(0:imax+1,0:jmax+1,1,1) !17-11-17
       if(k0==nf_real) then !...>
         status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                             ,varcount(1:2)   &
                               ,anyvar2d(0:imax+1,0:jmax+1))
         h_w(0:imax+1,0:jmax+1)=anyvar2d(0:imax+1,0:jmax+1)
       endif                !...>

! Lire upwindriver_t !28-07-21
! Note: on lit upwindriver_t de la simulation physique offline sans l'operation "min" de
! la valeur precedente
       status=nf_inq_varid(ncid_,'upwindriver_t',var_id)
       if(status/=0)stop &
       ' stop nf_inq_varid upwindriver_t grid_mask_h_offlinecase'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2 
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2 
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_real)stop 'k0/=nf_real upwindriver_t grid_mask_h_offlinecase'
       status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                           ,varcount(1:2)   &
                        ,upwindriver_t(0:imax+1,0:jmax+1))

       if(status/=0)stop &
       ' Erreur nf_get_vara_double h_w grid_mask_h_offlinecase'

      status=nf_close(ncid_)
      if(status/=0)stop ' Erreur nf_close grid_mask_h_offlinecase'


      end subroutine grid_mask_h_offlinecase

!...................................................................

      subroutine grid_sigma_offlinecase
      implicit none
      integer ncid_
#ifdef synopsis
       subroutinetitle='grid_sigma_offlinecase'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'Initial dsig_t'
!     status=1
!     k0=0
!     do while (status/=0) !---while status/=0---> !04-05-18
!      if(k0==0)texte250=trim(directory_offline)//'grille.nc'
!      if(k0==1)texte250=trim(directory_offline)//'grid.nc'
!      if(k0==2)texte250=trim(directory_offline)//'grid_full.nc' !04-05-18
!      if(k0==3)stop 'No grid file in grid_sigma_offlinecase'
!      if(par%rank==0)write(6,'(a,a)')'try nf_open ',trim(texte250)
!      status=nf_open(trim(texte250),nf_nowrite,ncid_)
!      k0=k0+1 
!      if(par%rank==0.and.status==0)write(6,'(a,a)')'OK  nf_open ',trim(texte250)
!     enddo                !---while status/=0--->
       status=nf_open(trim(lonlatfile),nf_nowrite,ncid_) !22-05-18
       if(status/=0) then !>>>
        write(6,'(a,a)')'status/=0 nf_open(trim(lonlatfile) for file=',trim(lonlatfile)
        stop 'Err 1785 nf_open(trim(lonlatfile)'
       else               !>>>
        if(par%rank==0)write(6,'(a,a)')'try nf_open ',trim(lonlatfile) !24-11-20
       endif              !>>>

       status=nf_inq_varid(ncid_,'sigma_w',var_id)

       if(status==0) then !oooooo>

       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0==nf_double) then !ddddd>
       status=nf_get_vara_double(ncid_,var_id,varstart(1:3)   &
                                             ,varcount(1:3)   &
                     ,anyv3d(0:imax+1,0:jmax+1,1:kmax+1,1))
         do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          dsig_t(i,j,k)=anyv3d(i,j,k+1,1)-anyv3d(i,j,k,1)
         enddo ; enddo ; enddo
       endif                  !ddddd>
       if(k0==nf_real) then !...>
         status=nf_get_vara_real(ncid_,var_id,varstart(1:3)   &
                                             ,varcount(1:3)   &
                      ,anyvar3d(0:imax+1,0:jmax+1,1:kmax+1))
         do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          dsig_t(i,j,k)=anyvar3d(i,j,k+1)-anyvar3d(i,j,k)
         enddo ; enddo ; enddo
       endif                !...>

       else               !oooooo> !20-12-14

        varstart(3)=1              ; varcount(3)=kmax
        status=nf_inq_varid(ncid_,'dsig_t',var_id)
        if(status/=0)stop &
       ' stop nf_inq_varid dsig_ tgrid_sigma_offlinecase'
        status=nf_inq_vartype(ncid_,var_id,k0)
        if(k0/=nf_double)stop 'vartyp dsig_t/=nf_double'
        status=nf_get_vara_double(ncid_,var_id,varstart(1:3)   &
                                              ,varcount(1:3)   &
     !                     ,dsig_t(0:imax+1,0:jmax+1,1:kmax))
                           ,anyv3d(0:imax+1,0:jmax+1,1:kmax,1)) !24-11-20
         do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          dsig_t(i,j,k)=anyv3d(i,j,k,1) !24-11-20
         enddo ; enddo ; enddo

       endif              !oooooo>


! LECTURE kmin_w et kmerge_t !26-03-18
       status=nf_inq_varid(ncid_,'kmin_w',var_id)
       if(status/=0)stop 'Err 1810 nf_inq_varid kmin_w'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_short)stop 'vartyp kmin_w/=nf_short'
       status=nf_get_vara_int(ncid_,var_id,varstart(1:2)   &
                                          ,varcount(1:2)   &
                              ,kmin_w(0:imax+1,0:jmax+1))
       if(status/=0)stop 'Err 1816 nf_get_vara_int kmin_w'

       status=nf_inq_varid(ncid_,'kmerged_t',var_id)
       if(status/=0)stop 'Err 1819 nf_inq_varid kmerged_t'
       status=nf_inq_vartype(ncid_,var_id,k0)
       if(k0/=nf_short)stop 'vartyp kmerged_t/=nf_short'
       status=nf_get_vara_int(ncid_,var_id,varstart(1:2)   &
                                          ,varcount(1:2)   &
                           ,kmerged_t(0:imax+1,0:jmax+1))
       if(status/=0)stop 'Err 1825 nf_get_vara_int kmerged_t'

      status=nf_close(ncid_)
      if(status/=0)stop ' Erreur nf_close grid_sigma_offlinecase'

! Deduire sigma_w de dsig_t
       sigma_w(:,:,1)=0.
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
       sigma_w(i,j,k+1)=sigma_w(i,j,k)+dsig_t(i,j,k) !25-11-20
       enddo       ; enddo         ; enddo
!      i=imax/2; j=jmax/2 ; write(10+par%Rank,*)'sigma',sigma_w(i,j,kmax+1)
! seuils et conditions limites inferieures et superieures
       dsig_t=max(dsig_t,small1)            !25-11-20
       dsig_t(:,:,kmax+1)=dsig_t(:,:,kmax)
       dsig_t(:,:,0     )=dsig_t(:,:,1   )

! kmin_u kmin_v kmerged_u kmerged_v kundermin_t kundermin_u kundermin_v
! upwzone0_t
      call sigma_levels_merged_finalize !26-03-18

      end subroutine grid_sigma_offlinecase

!...................................................................

#ifdef egrid
      subroutine grid_rotation_f_location !28-12-14
      implicit none
      double precision angle_,direction_
#ifdef synopsis
       subroutinetitle='grid_rotation_f_location'
       subroutinedescription=                                     &
       'Computes gridrotcos_f & gridrotsin_f'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(gridrotcos_f))allocate(gridrotcos_f(1:imax+1,1:jmax+1,8)) ; gridrotcos_f=0
      if(.not.allocated(gridrotsin_f))allocate(gridrotsin_f(1:imax+1,1:jmax+1,8)) ; gridrotsin_f=0
      if(.not.allocated(cocosisi_f)  )allocate(  cocosisi_f(1:imax+1,1:jmax+1,1)) ; cocosisi_f=0

! Tableaux de rotation si axes non paralleles lignes lon et lat   !07-03-13
      do j=1,jmax+1
      do i=1,imax+1


!..........................................................
! Partie 1: Axes diagonaux

       call lonlat2distance(lon_t(i  ,j  ),lat_t(i  ,j  )  &
                           ,lon_t(i-1,j-1),lat_t(i-1,j-1),dist1)
       dx_e_f(i,j)=dist1
!      call lonlat2distance(lon_t(i  ,j  ),lat_t(i-1,j-1)  &
!                          ,lon_t(i-1,j-1),lat_t(i-1,j-1),dist2)
!      call lonlat2distance(lon_t(i  ,j  ),lat_t(i-1,j-1)  &
!                          ,lon_t(i  ,j  ),lat_t(i  ,j  ),dist3)
!      dist1=dist1/rayonterre ; dist2=dist2/rayonterre ; dist3=dist3/rayonterre

! Angle Axe Oi Diagonal par rapport a axe Ouest-Est mesure depuis l'axe O-E autrement
! dit dans une grille standard a maille carree l'angle doit etre  45°
       x1=lon_t(i,j)-lon_t(i-1,j-1)
       if(x1<-pi)x1=x1+2.*pi
       if(x1> pi)x1=x1-2.*pi
       x0= atan2(lat_t(i,j)-lat_t(i-1,j-1)    &
             ,x1*cos(0.25*( lat_t(i  ,j  )    &
                           +lat_t(i-1,j  )    &
                           +lat_t(i  ,j-1)    &
                           +lat_t(i-1,j-1))))
!      x0=acos( (cos(dist3)-cos(dist1)*cos(dist2))/(sin(dist1)*sin(dist2)))
       gridrotcos_f(i,j,1)=cos(x0)
       gridrotsin_f(i,j,1)=sin(x0)
!      write(66,*)i,j,x0*rad2deg,rad2deg*acos( (cos(dist3)-cos(dist1)*cos(dist2))/(sin(dist1)*sin(dist2)))
!      write(66,*)dist1*rayonterre,dist2*rayonterre,dist3*rayonterre

!.............
       call lonlat2distance(lon_t(i-1,j  ),lat_t(i-1,j  )  &
                           ,lon_t(i  ,j-1),lat_t(i  ,j-1),dist1)
       dy_e_f(i,j)=dist1
!      call lonlat2distance(lon_t(i  ,j-1),lat_t(i-1,j  )  &
!                          ,lon_t(i  ,j-1),lat_t(i  ,j-1),dist2)
!      call lonlat2distance(lon_t(i  ,j-1),lat_t(i-1,j  )  &
!                          ,lon_t(i-1,j  ),lat_t(i-1,j  ),dist3)
!      dist1=dist1/rayonterre ; dist2=dist2/rayonterre ; dist3=dist3/rayonterre

! Angle Axe Oj Diagonal par rapport a axe Sud-Nord mesure depuis l'axe S-N autrement
! dit dans une grille standard a maille carree l'angle doit etre +45°
       x1=lon_t(i-1,j)-lon_t(i,j-1)
       if(x1<-pi)x1=x1+2.*pi
       if(x1> pi)x1=x1-2.*pi
       x2= atan2(lat_t(i-1,j)-lat_t(i,j-1)    &
             ,x1*cos(0.25*( lat_t(i  ,j  )    &
                           +lat_t(i-1,j  )    &
                           +lat_t(i  ,j-1)    &
                           +lat_t(i-1,j-1))))-0.5*pi
       gridrotcos_f(i,j,2)=cos(x2)
       gridrotsin_f(i,j,2)=sin(x2)
!      write(67,*)i,j,x2*rad2deg
!      write(67,*)i,j,x2*rad2deg,rad2deg*acos( (cos(dist3)-cos(dist1)*cos(dist2))/(sin(dist1)*sin(dist2)))

       gridrotcos_f(i,j,3)=     gridrotcos_f(i,j,1)     &
         /( gridrotcos_f(i,j,1)*gridrotcos_f(i,j,2)     &
           +gridrotsin_f(i,j,1)*gridrotsin_f(i,j,2) )

       gridrotsin_f(i,j,3)=     gridrotsin_f(i,j,1)     &
         /( gridrotcos_f(i,j,1)*gridrotcos_f(i,j,2)     &
           +gridrotsin_f(i,j,1)*gridrotsin_f(i,j,2) )

       gridrotcos_f(i,j,4)=     gridrotcos_f(i,j,2)     &
         /( gridrotcos_f(i,j,1)*gridrotcos_f(i,j,2)     &
           +gridrotsin_f(i,j,1)*gridrotsin_f(i,j,2) )

       gridrotsin_f(i,j,4)=     gridrotsin_f(i,j,2)     &
         /( gridrotcos_f(i,j,1)*gridrotcos_f(i,j,2)     &
           +gridrotsin_f(i,j,1)*gridrotsin_f(i,j,2) )

       cocosisi_f(i,j,1)=                               &
            gridrotcos_f(i,j,1)*gridrotcos_f(i,j,2)     &
           +gridrotsin_f(i,j,1)*gridrotsin_f(i,j,2)

!..........................................................
! Partie 2: Axe horizontal Axe vertical
! Angle Axe Oi horizontal par rapport a axe Ouest-Est mesure depuis l'axe O-E autrement
! dit dans une grille standard a maille carree l'angle doit etre 0°
       x1=lon_t(i,j)-lon_t(i-1,j)
       if(x1<-pi)x1=x1+2.*pi
       if(x1> pi)x1=x1-2.*pi

! Algo pour centrer (mais qu'on laisse de cote puisque le PGF n'est pas centre...)
!      x11=x1
!      x1= lon_t(i,j-1)-lon_t(i-1,j-1)
!      if(x1<-pi)x1=x1+2.*pi
!      if(x1> pi)x1=x1-2.*pi
!      x1=0.5*(x1+x11)
!      x0= atan2( 0.5*(lat_t(i,j  )-lat_t(i-1,j  )+lat_t(i,j-1)-lat_t(i-1,j-1))   &
!            ,x1*cos(0.25*( lat_t(i  ,j  )+lat_t(i-1,j  )+lat_t(i  ,j-1)+lat_t(i-1,j-1))))

       x0=atan2(lat_t(i,j)-lat_t(i-1,j),x1*cos(0.5*(lat_t(i,j)+lat_t(i-1,j))))
       gridrotcos_f(i,j,5)=cos(x0)
       gridrotsin_f(i,j,5)=sin(x0)

! Angle Axe Oj vertical par rapport a axe Sud-Nord mesure depuis l'axe S-N autrement
! dit dans une grille standard a maille carree l'angle doit etre +45°
       x1=lon_t(i,j)-lon_t(i,j-1)
       if(x1<-pi)x1=x1+2.*pi
       if(x1> pi)x1=x1-2.*pi
! Algo pour centrer (mais qu'on laisse de cote puisque le PGF n'est pas centre...)
!      x11=x1
!      x1=lon_t(i-1,j)-lon_t(i-1,j-1)
!      if(x1<-pi)x1=x1+2.*pi
!      if(x1> pi)x1=x1-2.*pi
!      x1=0.5*(x1+x11)
!      x2= atan2(0.5*(lat_t(i  ,j)-lat_t(i  ,j-1)+lat_t(i-1,j)-lat_t(i-1,j-1))   &
!            ,x1*cos(0.25*( lat_t(i  ,j  )+lat_t(i-1,j  )+lat_t(i  ,j-1)+lat_t(i-1,j-1))))-0.5*pi
       x2=atan2(lat_t(i,j)-lat_t(i,j-1),x1*cos(0.5*(lat_t(i,j)+lat_t(i,j-1))))-0.5*pi
       gridrotcos_f(i,j,6)=cos(x2)
       gridrotsin_f(i,j,6)=sin(x2)

       write(67,*)i,j,x0*rad2deg,x2*rad2deg

! Normalisation
        x5=1./( gridrotcos_f(i,j,5)*gridrotcos_f(i,j,6)  &
               +gridrotsin_f(i,j,5)*gridrotsin_f(i,j,6) )

       gridrotcos_f(i,j,7)=x5*gridrotcos_f(i,j,5)
       gridrotsin_f(i,j,7)=x5*gridrotsin_f(i,j,5)
       gridrotcos_f(i,j,8)=x5*gridrotcos_f(i,j,6)
       gridrotsin_f(i,j,8)=x5*gridrotsin_f(i,j,6)

      enddo
      enddo

      end subroutine grid_rotation_f_location
#endif

!...................................................................

#ifdef egrid
      subroutine grid_quadrangle_area
      implicit none
#ifdef synopsis
       subroutinetitle='grid_quadrangle_area'
       subroutinedescription=                                          &
       'Computes the area of a triangle'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      real*8 p_,a_,b_,c_

! a,b,c = length of the sides of the triangle
! p=0.5*(a+b+c)
! R=earth radius
! Area=4*R**2*arctan( sqrt( tan(0.5*p)*tan(0.5*(p-a))*tan(0.5*(p-b))*tan(0.5*(p-c)) ) )


! dxdy_t is the area of the small quadrangle defining the "c" box nested into the "e" box
! The 4 corners of this box are the 4 "f" points around the "t" point.
      do j=0,jmax+1 ; do i=0,imax+1
! Area of the triangle defined by the 3 points _f(i+1,j+1), _f(i+1,j), _f(i,j+1)
       call lonlat2distance(lon_f(i+1,j+1),lat_f(i+1,j+1)    &
                           ,lon_f(i+1,j  ),lat_f(i+1,j  ),a_)
       call lonlat2distance(lon_f(i+1,j+1),lat_f(i+1,j+1)    &
                           ,lon_f(i  ,j+1),lat_f(i  ,j+1),b_)
       call lonlat2distance(lon_f(i+1,j  ),lat_f(i+1,j  )    &
                           ,lon_f(i  ,j+1),lat_f(i  ,j+1),c_)
! Huilier method:
       a_=a_/rayonterre
       b_=b_/rayonterre
       c_=c_/rayonterre
       p_=0.5*(a_+b_+c_)
       x1=4*rayonterre**2*atan( sqrt( tan(0.5*p_)*tan(0.5*(p_-a_))*tan(0.5*(p_-b_))*tan(0.5*(p_-c_)) ) )
! Area of the triangle defined by the 3 points _f(i  ,j  ), _f(i+1,j), _f(i,j+1)
       call lonlat2distance(lon_f(i  ,j  ),lat_f(i  ,j  )    &
                           ,lon_f(i+1,j  ),lat_f(i+1,j  ),a_)
       call lonlat2distance(lon_f(i  ,j  ),lat_f(i  ,j  )    &
                           ,lon_f(i  ,j+1),lat_f(i  ,j+1),b_)
! quant a c_ il est deja calcule....
       a_=a_/rayonterre
       b_=b_/rayonterre
       p_=0.5*(a_+b_+c_)
       x2=4*rayonterre**2*atan( sqrt( tan(0.5*p_)*tan(0.5*(p_-a_))*tan(0.5*(p_-b_))*tan(0.5*(p_-c_)) ) )
! Area of the quadrangle given by the sum of the 2 triangle areas:
       dxdy_t(i,j)=x1+x2
       write(78,*)sqrt(dxdy_t(i,j)) &
!        ,sqrt(x1+x2) &
        ,dy_t(i,j),dx_t(i,j)
      enddo ; enddo

! dxdy_e_t is the area of the big quadrangle defining the "e" box
! Area of the quadrangle (i-1,j) (i+1,j) (i,j+1) (i,j-1) given by the sum of the 2 triangle areas
! The 4 corners of this box are the 4 "t" neighbor points around the "t" central point.
      do j=1,jmax  ; do i=1,imax
! Area of the triangle defined by the 3 points (i,j-1), (i+1,j), (i-1,j):
       call lonlat2distance(lon_t(i  ,j-1),lat_t(i  ,j-1)     &
                           ,lon_t(i+1,j  ),lat_t(i+1,j  ),a_)
       call lonlat2distance(lon_t(i  ,j-1),lat_t(i  ,j-1)     &
                           ,lon_t(i-1,j  ),lat_t(i-1,j  ),b_)
       call lonlat2distance(lon_t(i+1,j  ),lat_t(i+1,j  )  &
                           ,lon_t(i-1,j  ),lat_t(i-1,j  ),c_)
! Heron method (neglecting sphericity):
!     p_=0.5*(a_+b_+c_)
!     x1=sqrt( p_*(p_-a_)*(p_-b_)*(p_-c_))
!     write(6,*)'Triangle area Heron method  ',x1,sqrt(x1)
! Huilier method:
       a_=a_/rayonterre
       b_=b_/rayonterre
       c_=c_/rayonterre
       p_=0.5*(a_+b_+c_)
       x1=4*rayonterre**2*atan( sqrt( tan(0.5*p_)*tan(0.5*(p_-a_))*tan(0.5*(p_-b_))*tan(0.5*(p_-c_)) ) )
!     write(6,*)'Triangle area Huilier method',x1,sqrt(x1)
! Area of the triangle defined by the 3 points (i,j+1), (i+1,j), (i-1,j):
       call lonlat2distance(lon_t(i  ,j+1),lat_t(i  ,j+1)     &
                           ,lon_t(i+1,j  ),lat_t(i+1,j  ),a_)
       call lonlat2distance(lon_t(i  ,j+1),lat_t(i  ,j+1)     &
                           ,lon_t(i-1,j  ),lat_t(i-1,j  ),b_)
! quant a c_ il est deja calcule....
       a_=a_/rayonterre
       b_=b_/rayonterre
       p_=0.5*(a_+b_+c_)
       x2=4*rayonterre**2*atan( sqrt( tan(0.5*p_)*tan(0.5*(p_-a_))*tan(0.5*(p_-b_))*tan(0.5*(p_-c_)) ) )
!     write(6,*)'Triangle area Huilier method',x2,sqrt(x2)
! Area of the quadrangle (i-1,j) (i+1,j) (i,j+1) (i,j-1) given by the sum of the 2 triangle areas
       dxdy_e_t(i,j)=x1+x2
       write(79,*)sqrt(dxdy_e_t(i,j)) &
!        ,sqrt(x1+x2) &
        ,dy_f(i,j),dx_f(i,j)
!       ,sqrt(dx_t(i,j)**2+dy_t(i,j)**2)
      enddo ; enddo

      end subroutine grid_quadrangle_area
#endif

!...................................................................

      subroutine grid_webcanals_readfiles(field_)
      implicit none
      character(len=*) field_
      integer :: unit3_,unit4_,icnl_,point_

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_3
! coordonnees i,j des differents points:
! i_canalcoord(numero canal,point,gridtype,sender_or_receiver)
! j_canalcoord(numero canal,point,gridtype,sender_or_receiver)

! Point 1
! i_canalcoord(k,1,1,2),j_canalcoord(k,1,1,2) = i,j of the first point of the canal number k
! i_canalcoord(k,1,2,1),j_canalcoord(k,1,2,1) = i,j of the first point of the canal number k at its corresponding sea location

! Point 2
! i_canalcoord(k,2,1,2),j_canalcoord(k,2,1,2) = i,j of the last  point of the canal number k
! i_canalcoord(k,2,2,1),j_canalcoord(k,2,2,1) = i,j of the last  point of the canal number k at its corresponding sea location

! Direction du flux entrant en mer avec convention identique A celle de notebook_rivers
! canaldir( numero canal , numero du point 1:2 ) = 1,2,3,4
! le point 1 ou 2 entre en mer dans le sens:
!                                           i   croissant => convention 1
!                                           j   croissant => convention 2
!                                           i decroissant => convention 3
!                                           j decroissant => convention 4

      if(webcanals_list=='none')return !24-04-19

                unit3_=s_unit(7)
      open(unit=unit3_,file=trim(webcanals_list))
      read(unit3_,*)nbcanal  
      if(nbcanal<=0) then ! if no webcanal close & return
       close(unit3_) ; return
      endif               ! if no webcanal close & return

      read(unit3_,*)
      if(.not.allocated(i_canalcoord)) then  !allo-allo>
!               allocate(canalcoord(8,nbcanal)) ; canalcoord=0 
               allocate(i_canalcoord(nbcanal,2,2,2)) ; i_canalcoord=-999
               allocate(j_canalcoord(nbcanal,2,2,2)) ; j_canalcoord=-999
               allocate(    canaldir(nbcanal,2))     ; canaldir=-99 ! -99 nofield value for integer(kind=1)
      endif                                !allo-allo>

      point1=1    ; point2=2
      gridtype1=1 ; gridtype2=2
      sender=1    ; receiver=2
      do k=1,nbcanal ! (numero du canal)
       read(unit3_,'(a)')texte90
! ouvrir le fichier du canal numero "k" et charger les coordonnees puis h,dx,dy
                 unit4_=s_unit(8)
       open(unit=unit4_,file=trim(texte90))

         read(unit4_,*) & ! coords du point 1 dans le canal et en mer et direction
         i_canalcoord(k,point1,gridtype1,receiver),j_canalcoord(k,point1,gridtype1,receiver) &
        ,i_canalcoord(k,point1,gridtype2,sender)  ,j_canalcoord(k,point1,gridtype2,sender)   &
            ,canaldir(k,point1)
         read(unit4_,*) & ! coords du point 2 dans le canal et en mer et direction
         i_canalcoord(k,point2,gridtype1,receiver),j_canalcoord(k,point2,gridtype1,receiver) &
        ,i_canalcoord(k,point2,gridtype2,sender)  ,j_canalcoord(k,point2,gridtype2,sender)   &
            ,canaldir(k,point2)

! Les coordonnees du point "sender" dans la grille 1 sont deduites de celles du point "receiver"
         i_canalcoord(k,point1,gridtype1,sender)=i_canalcoord(k,point1,gridtype1,receiver)+1
         i_canalcoord(k,point2,gridtype1,sender)=i_canalcoord(k,point2,gridtype1,receiver)-1
         j_canalcoord(k,point1,gridtype1,sender)=j_canalcoord(k,point1,gridtype1,receiver)
         j_canalcoord(k,point2,gridtype1,sender)=j_canalcoord(k,point2,gridtype1,receiver)

! Les coordonnees du point "receiver" dans la grille 2 sont deduites de celles du point "sender"
! et de l'orientation du flux au point de connexion dans la grille 2
         i_canalcoord(k,point1,gridtype2,receiver)=i_canalcoord(k,point1,gridtype2,sender)
         j_canalcoord(k,point1,gridtype2,receiver)=j_canalcoord(k,point1,gridtype2,sender)
         i_canalcoord(k,point2,gridtype2,receiver)=i_canalcoord(k,point2,gridtype2,sender)
         j_canalcoord(k,point2,gridtype2,receiver)=j_canalcoord(k,point2,gridtype2,sender)
         if(canaldir(k,point1)==1)i_canalcoord(k,point1,gridtype2,receiver)=i_canalcoord(k,point1,gridtype2,sender)-1
         if(canaldir(k,point2)==1)i_canalcoord(k,point2,gridtype2,receiver)=i_canalcoord(k,point2,gridtype2,sender)-1
         if(canaldir(k,point1)==3)i_canalcoord(k,point1,gridtype2,receiver)=i_canalcoord(k,point1,gridtype2,sender)+1
         if(canaldir(k,point2)==3)i_canalcoord(k,point2,gridtype2,receiver)=i_canalcoord(k,point2,gridtype2,sender)+1
         if(canaldir(k,point1)==2)j_canalcoord(k,point1,gridtype2,receiver)=j_canalcoord(k,point1,gridtype2,sender)-1
         if(canaldir(k,point2)==2)j_canalcoord(k,point2,gridtype2,receiver)=j_canalcoord(k,point2,gridtype2,sender)-1
         if(canaldir(k,point1)==4)j_canalcoord(k,point1,gridtype2,receiver)=j_canalcoord(k,point1,gridtype2,sender)+1
         if(canaldir(k,point2)==4)j_canalcoord(k,point2,gridtype2,receiver)=j_canalcoord(k,point2,gridtype2,sender)+1



        if(field_=='hm') then !pmx> !23-06-19
! A l'etape 'hm' (c.a.d. h et mask, appel depuis initial_mask_and_bathy)
!  - Demasquer le receiver de la grille 2:
! Afin que les flux transmis au point receveir en mer ne soient pas
! annulEs par le mask initial, on demasque ce point !23-06-19
        do point_=point1,point2
         i=i_canalcoord(k,point_,gridtype2,receiver) ; j=j_canalcoord(k,point_,gridtype2,receiver)
         if(i-par%timax(1)>=0.and.i-par%timax(1)<=imax+1.and. &
            j-par%tjmax(1)>=0.and.j-par%tjmax(1)<=jmax+1) then !mask>

              mask_t(i-par%timax(1),j-par%tjmax(1),:)=1

! Cependant on ne peut pas ouvrir une frontiere fermee. Une connexion n'est donc
! pas possible si elle entraine de demasquer en i=1, i=iglb, i=1, j=jglb
           if(i==1.or.i==iglb.or.j==1.or.j==jglb) then !stop-debug>
             flag_stop=1
             write(10+par%rank,*)'Can not change border land/sea for' & 
            ,' channel',k, 'point 1'
           endif                                       !stop-debug>
         endif                                                 !mask>
        enddo ! boucle sur point_

        endif                 !pmx>

        read(unit4_,*)
        read(unit4_,'(a)')texte60 ! format....!12-02-22
        read(unit4_,*)
        do j=j_canalcoord(k,1,1,2),j_canalcoord(k,2,1,2)
        do i=i_canalcoord(k,1,1,2),i_canalcoord(k,2,1,2)

! Lire i,j,h_w(i,j),dx_t(i,j),dy_t(i,j):
         if(index(texte60,'lon')/=0.or. & !05-05-19
            index(texte60,'lat')/=0) then !ooo>
             read(unit4_,*)i1,j1,x1,x2,x3,x4,x5 !x1,x2,x3,x4,x5= h,dx,dy,lon,lat
         else                             !ooo>
             read(unit4_,*)i1,j1,x1,x2,x3       !x1,x2,x3= h,dx,dy
         endif                            !ooo>

! Verification de la coherence des indices i,j:
         if(i/=i1.or.j/=j1) then !>>> debug
          flag_stop=1
          write(10+par%rank,'(2a)')trim(texte90),' Err: i/=i1.or.j/=j1'
          write(10+par%rank,*)'i,j,i1,j1',i,j,i1,j1
         endif                   !>>> debug
         if(field_=='hm')glob_h(i,j)=x1 ! charger aussi glob_h est necessaire pour l'etape C a suivre

! Si le point i,j est dans le rank local charger les champs:
         if(i-par%timax(1)>=0.and.i-par%timax(1)<=imax+1.and. &
            j-par%tjmax(1)>=0.and.j-par%tjmax(1)<=jmax+1) then !ooo> inside local rank
           if(field_=='hm') then   !bathy & mask>
                            h_w(i-par%timax(1),j-par%tjmax(1))=x1
           endif                   !bathy & mask>
           if(field_=='dxdy') then !dx & dy>
           if(x2<0.and.canaldir(k,1)<0.and.i==i_canalcoord(k,1,1,2)) then !--->
! Attention ce cas particulier ne peut pas deduire dx de lon,lat si canaldir(k,1)<0 et i=i_canalcoord(k,1,1,2), 
! il garde sa valeur par defaut dans ce cas donc:                !16-11-21
                          x2=dx_t(i-par%timax(1),j-par%tjmax(1)) !16-11-21
           endif                                                          !--->
                             dx_t(i-par%timax(1),j-par%tjmax(1))=x2
                             dy_t(i-par%timax(1),j-par%tjmax(1))=x3
           endif                   !dx & dy>
         endif                                               !ooo> inside local rank

         if(field_=='lonlat') then !lonlat>
         if(i-par%timax(1)>=-1.and.i-par%timax(1)<=imax+2.and. &
            j-par%tjmax(1)>=-1.and.j-par%tjmax(1)<=jmax+2) then !ooo> inside local rank lon lat

         if(index(texte60,'lon')/=0.or. & !05-05-19
            index(texte60,'lat')/=0) then !ooo>
            lon_t(i-par%timax(1),j-par%tjmax(1))=x4*deg2rad
            lat_t(i-par%timax(1),j-par%tjmax(1))=x5*deg2rad
          endif                           !ooo>

         endif                                                  !ooo> inside local rank lon lat
         endif                     !lonlat>

        enddo ! i
        enddo ! j
       close(unit4_)
      enddo ! k (numero du canal)
      close(unit3_)

! flag_stop ?
      k0=0 
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop &
             ' Stop grid_webcanals_readfiles see fort.xxx err files' !29-01-20


! A l'etape mask and bathy faire glob_mask avant de tester les points d'entrees des rivieres
! afin de pouvoir introduire une riviere (via notebook_river) au bout d'un canal
! sans declencher les alarmes d'incoherence de positionnement des sources
      if(field_=='hm') then   !bathy & mask>

! Dans le canal "partie grille canal":
! A: demasquer les points
       do icnl_=1,nbcanal 
       do j=j_canalcoord(icnl_,point1,gridtype1,receiver),j_canalcoord(icnl_,point2,gridtype1,receiver)
       do i=i_canalcoord(icnl_,point1,gridtype1,receiver),i_canalcoord(icnl_,point2,gridtype1,receiver)
         glob_mask(i,j)=1
         if(i-par%timax(1)>=0.and.i-par%timax(1)<=imax+1.and. &
            j-par%tjmax(1)>=0.and.j-par%tjmax(1)<=jmax+1)mask_t(i-par%timax(1),j-par%tjmax(1),:)=1
       enddo ; enddo ; enddo

! Cas partculier de l'extremite point 1 du canal si non connectee A un
! point de mer (donc A priori intro d'une riviere via notebook_river):
! le point receiver est dans ce cas masquE et ses coordonnees sont
! celles qu'on indique dans notebook_rivers: !16-11-21
       do icnl_=1,nbcanal !16-11-21
        if(canaldir(icnl_,point1)<=0) then !-point-non-connecte->
         i=i_canalcoord(icnl_,point1,gridtype1,receiver)
         j=j_canalcoord(icnl_,point1,gridtype1,receiver)
         glob_mask(i,j)=0
         if(i-par%timax(1)>=0.and.i-par%timax(1)<=imax+1.and. &
            j-par%tjmax(1)>=0.and.j-par%tjmax(1)<=jmax+1)mask_t(i-par%timax(1),j-par%tjmax(1),:)=0
        endif                              !-point-non-connecte->
       enddo

       do icnl_=1,nbcanal

! B: a l'extremitE des canaux la bathy est imposee par la bathy des points en mer correspondant
!    si point non connectE ne pas faire !17-06-19
       do point_=1,2
       if(     i_canalcoord(icnl_,point_,gridtype2,sender)>0   &
          .and.j_canalcoord(icnl_,point_,gridtype2,sender)>0) then !m°v°m> !17-06-19
         i=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1) 
         j=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1) ! point dans la grille canal
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !>>>
          h_w(i,j)=glob_h(i_canalcoord(icnl_,point_,gridtype2,sender) &
                         ,j_canalcoord(icnl_,point_,gridtype2,sender))  !glob_h point 1 dans la grille mer
         endif                                              !>>>
       endif                                                       !m°v°m> !17-06-19
       enddo ! point_

! C: juste avant l'extremitE des canaux, les canaux imposent la bathy au point de mer correspondant
!    si point non connectE ne pas faire !17-06-19
       do point_=1,2
       if(     i_canalcoord(icnl_,point_,gridtype2,sender)>0   &
          .and.j_canalcoord(icnl_,point_,gridtype2,sender)>0) then !m°v°m> !17-06-19
         i=i_canalcoord(icnl_,point_,gridtype2,receiver) 
         j=j_canalcoord(icnl_,point_,gridtype2,receiver) ! point 1 dans la grille mer
         i=i-par%timax(1) ; j=j-par%tjmax(1)
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !>>>
          h_w(i,j)=glob_h(i_canalcoord(icnl_,point_,gridtype1,sender)  &
                         ,j_canalcoord(icnl_,point_,gridtype1,sender))  !glob_h juste avant l'extremite du point 1 dans la grille canal
         endif                                              !>>>
       endif                                                       !m°v°m> !17-06-19
       enddo ! point_
         
       enddo ! icnl_ canal

      endif                   !bathy & mask>

! Verifier que les points ont tous des coordonnees distinctes:
      flag_stop=0
      do k4=gridtype1,gridtype2
       if(k4==gridtype1)k5=receiver
       if(k4==gridtype2)k5=sender
       do k3=point1,point2
       do k1=1   ,nbcanal
       do k2=k1+1,nbcanal

        if(i_canalcoord(k1,k3,k4,k5)>0.and.j_canalcoord(k1,k3,k4,k5)>0) then !ooo> !17-07-19

        if(i_canalcoord(k1,k3,k4,k5)==i_canalcoord(k2,point1,k4,k5).and. &
           j_canalcoord(k1,k3,k4,k5)==j_canalcoord(k2,point1,k4,k5)) then !>>>>
         flag_stop=1 
         write(10+par%rank,'(3(a,i4))')'    channel',k1,' point',k3 &
         ,' gridtype=',k4
         write(10+par%rank,'(3(a,i4))')'and channel',k2,' point',point1 &
         ,' gridtype=',k4
         write(10+par%rank,'(a,2i5)')'have same coordinates ',i_canalcoord(k1,k3,k4,k5),j_canalcoord(k1,k3,k4,k5)
        endif                                                         !>>>>

        if(i_canalcoord(k1,k3,k4,k5)==i_canalcoord(k2,point2,k4,k5).and. &
           j_canalcoord(k1,k3,k4,k5)==j_canalcoord(k2,point2,k4,k5)) then !>>>>
         flag_stop=1 
         write(10+par%rank,'(3(a,i4))')'    channel',k1,' point',k3 &
         ,' gridtype=',k4
         write(10+par%rank,'(3(a,i4))')'and channel',k2,' point',point2 &
         ,' gridtype=',k4
         write(10+par%rank,'(a,2i5)')'have same coordinates ',i_canalcoord(k1,k3,k4,k5),j_canalcoord(k1,k3,k4,k5)
        endif                                                         !>>>>

        endif                                                                !ooo> !17-07-19

       enddo  !k2
       enddo  !k1
       enddo !k3
      enddo !k4
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0>0) stop &
      'Several canals have identical coordinates. See fort.xx err files'

      
! A la premiere etape (a priori celle qui calcule dx et dy) calculer les ranks des extremites
! des canaux et de leur connexion en mer:
      if(.not.allocated(canalrank))  then !rank>
              allocate(canalrank(nbcanal,2,2,2)) ; canalrank=-999 
              call grid_webcanals_rank

! Mettre les points de connexion en schema upwind
       do k1=1,nbcanal ; do k3=point1,point2 ; do k4=gridtype1,gridtype2 ; do k5=sender,receiver
         i=i_canalcoord(k1,k3,k4,k5)-par%timax(1)
         j=j_canalcoord(k1,k3,k4,k5)-par%tjmax(1)
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)upwindriver_t(i,j)=0.
       enddo ; enddo ; enddo ; enddo


      endif                               !rank>

! Longitudes, latitudes dans les canaux
! Connecter lon_t,lat_t aux points receveurs (des canaux)
! Deduire lon_u, lat_u dans les canaux  
! Connecter dx_t aux points receveurs des canaux

      if(field_=='lonlat') then !lonlat>
        do icnl_=1,nbcanal
         do point_=point1,point2
          i=i_canalcoord(icnl_,point_,gridtype1,2)-par%timax(1)
          j=j_canalcoord(icnl_,point_,gridtype1,2)-par%tjmax(1)
         enddo
        enddo
! Connecter lon_t,lat_t aux points receveurs (des canaux)
        call webcanals_gt2_to_gt1_any2d('lon_t') !03-05-19
        call webcanals_gt2_to_gt1_any2d('lat_t') !03-05-19

! continuite mpi sur lon_t, lat_
        call obc_lonlat !03-05-19

! Deduire lon_u, lat_u dans les canaux
       do icnl_=1,nbcanal
        j1=j_canalcoord(icnl_,point1,gridtype1,2)
        j=j1-par%tjmax(1)
        do i1=i_canalcoord(icnl_,point1,gridtype1,2)+1,i_canalcoord(icnl_,point2,gridtype1,2)
         i=i1-par%timax(1) 
         if(j>=0.and.j<=jmax+2.and.i>=0.and.i<=imax+2) then !>>>
           lat_u(i,j)=0.5*(lat_t(i-1,j)+lat_t(i,j))
           x1=lon_t(i-1,j)
           if(x1-lon_t(i,j)<-pi)x1=x1+2.*pi
           if(x1-lon_t(i,j)> pi)x1=x1-2.*pi
           lon_u(i,j)=0.5*(x1+lon_t(i,j))
         endif                                              !>>>
        enddo
       enddo

! Deduire gridrotcos_t et cie de lon_t,lat_t !18-08-22
       do icnl_=1,nbcanal
        j1=j_canalcoord(icnl_,point1,gridtype1,2)
        j=j1-par%tjmax(1)
        do i1=i_canalcoord(icnl_,point1,gridtype1,2),i_canalcoord(icnl_,point2,gridtype1,2)  
         i=i1-par%timax(1) 
! Note: ne pas disposer des valeurs en i+1 de l'extremite 2 et en i-1 de l'extremite 1 conduit au bornage de ip1 et im1
         ip1=min(i1+1,i_canalcoord(icnl_,point2,gridtype1,2))-par%timax(1)
         im1=max(i1-1,i_canalcoord(icnl_,point1,gridtype1,2))-par%timax(1)
         if(j>=0.and.j<=jmax+1.and.i>=0.and.i<=imax+1) then !>>>
! note: l'echange zc sur lon,lat permet d'avoir im1 en im1=-1 et ip1 en ip1=imax+2
            x1=lon_t(ip1,j)-lon_t(im1,j)
            if(x1<-pi)x1=x1+2.*pi
            if(x1> pi)x1=x1-2.*pi
            x0=-atan2(lat_t(ip1,j)-lat_t(im1,j),x1*cos(lat_t(i,j)))     !13-02-13
            gridrotcos_t(i,j)=cos(x0)
            gridrotsin_t(i,j)=sin(x0)
            grid_angle_t(i,j)=atan2(gridrotsin_t(i,j),gridrotcos_t(i,j))
         endif                                              !>>>
        enddo
       enddo

! Deduire la resolution dx des valeurs de longitude, latitude dans le
! cas oU les valeurs de dx donnEes dans le fichier sont negatives. Si
! les valeurs du fichiers sont positives elles sont considerees valides
! et par consequent elles ne sont pas recalculees
! La mer envoie aux points extremes des canaux (receveurs) dx_t ou dy_t selon l'orientation de la
! connexion:
        do icnl_=1,nbcanal
         do point_=point1,point2
          i=i_canalcoord(icnl_,point_,gridtype2,1)-par%timax(1)
          j=j_canalcoord(icnl_,point_,gridtype2,1)-par%tjmax(1)
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !ooo>
           if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
             xy_t(i,j,1)=dx_t(i,j)
           else                                                            !>>> 
             xy_t(i,j,1)=dy_t(i,j)
           endif                                                           !>>> 
          endif                                      !ooo>
         enddo
        enddo
        call webcanals_gt2_to_gt1_any2d('xy_t')
        call obc_ext_xy_t('za',1) !04-05-19
        do icnl_=1,nbcanal
         do point_=point1,point2
          i=i_canalcoord(icnl_,point_,gridtype1,2)-par%timax(1)
          j=j_canalcoord(icnl_,point_,gridtype1,2)-par%tjmax(1)
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !--->
              if(dx_t(i,j)<0.) & ! dx deduit de lon,lat seulement si dx<0 dans le fichier channel !03-05-19
              dx_t(i,j)=xy_t(i,j,1)
          endif                                      !--->
         enddo
        enddo

! Aux autres points des canaux dx_t est deduit de lon_u, lat_u
       do icnl_=1,nbcanal
        j1=j_canalcoord(icnl_,point1,gridtype1,2)
        j=j1-par%tjmax(1)
        do i1=i_canalcoord(icnl_,point1,gridtype1,1),i_canalcoord(icnl_,point2,gridtype1,1)
         i=i1-par%timax(1) 
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !ooo>
         x1=lon_u(i+1,j)
         if(x1-lon_u(i,j)<-pi)x1=x1+2.*pi
         if(x1-lon_u(i,j)> pi)x1=x1-2.*pi
      if(dx_t(i,j)<0.) & ! dx deduit de lon,lat seulement si dx<0 dans le fichier channel !03-05-19
         dx_t(i,j)=rayonterre*sqrt( ((     x1      - lon_u(i,j))*cos( lat_t(i,j)))**2       &
                                   +(  lat_u(i+1,j)- lat_u(i,j)                  )**2 )
         endif                                      !ooo>
        enddo ! i1
       enddo ! icnl_


      endif                     !lonlat>

      end subroutine grid_webcanals_readfiles

!......................................................................

      subroutine grid_webcanals_rank
      implicit none
      integer :: icnl_=0, point_=0, gridtype_=0, send_or_recv_=0, loop_=0

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_3
! calcule canalrank(index1,index2,index3,index4) 
! rank des points aux 2 bouts du canal et de leur position associee en mer.
! index1 = numero du canal
! index2 = extremite: 1=point1, 2=point2 (debut,fin)
! index3 = status: 1=canal, 2=mer
! index4 = status: 1=sender, 2=receveir

      do icnl_=1,nbcanal   ! numero de canal
      do point_=1,2        ! 1= point du debut, 2=point de fin
      do gridtype_=1,2     ! 1=extremite canal, 2=mer ou autre qu'extremite
      do send_or_recv_=1,2 ! 1=extremite canal, 2=mer ou autre qu'extremite

! Attention on ne cherche pas le rank si pas de connexion !17-06-19
! Si pas de connexion les coord du point en mer (c.a.d. gridtype2 iet sender) sont negatives (genre -99)
       if(     i_canalcoord(icnl_,point_,gridtype2,sender)>0   &
          .and.j_canalcoord(icnl_,point_,gridtype2,sender)>0) then !m°v°m> !17-06-19

      k0=0
      k1=0

       if(point_==1) then !11111>

        if(gridtype_==1) then !gt1> ! grille canal
         if(send_or_recv_==1) then !sender>
           i=i_canalcoord(icnl_,1,1,2)-par%timax(1)+1 ! note: il n'y a qu'une orientation 
           j=j_canalcoord(icnl_,1,1,2)-par%tjmax(1)   ! possible pour le canal: icroissant
         endif                     !sender>
         if(send_or_recv_==2) then !receiver>
           i=i_canalcoord(icnl_,1,1,2)-par%timax(1)   ! note: le receiver est comme une C.L.
           j=j_canalcoord(icnl_,1,1,2)-par%tjmax(1)   ! donc tout au bout du canal
         endif                     !receiver>
        endif                 !gt1> 

        if(gridtype_==2) then !gt2> ! grille mer
         if(send_or_recv_==1) then !sender>
           i=i_canalcoord(icnl_,1,2,1)-par%timax(1)   ! note: les coord en mer donnees en fichier 
           j=j_canalcoord(icnl_,1,2,1)-par%tjmax(1)   ! d'entree sont celles du sender
         endif                     !sender>
         if(send_or_recv_==2) then !receiver> ! note: le decallage en indice du point
           i=i_canalcoord(icnl_,1,2,1)-par%timax(1) ! receiver depend de l'orientation du
           j=j_canalcoord(icnl_,1,2,1)-par%tjmax(1) ! flux entrant en mer:
           if(canaldir(icnl_,point_)==1)i=i-1 ! 1=i croissant
           if(canaldir(icnl_,point_)==3)i=i+1 ! 3=i decroissant
           if(canaldir(icnl_,point_)==2)j=j-1 ! 2=j croissant
           if(canaldir(icnl_,point_)==4)j=j+1 ! 4=j decroissant
         endif                     !receiver>
        endif                 !gt2> 

       endif              !11111>

       if(point_==2) then !22222>

        if(gridtype_==1) then !gt1>
         if(send_or_recv_==1) then !sender>
           i=i_canalcoord(icnl_,2,1,2)-par%timax(1)-1 ! note: il n'y a qu'une orientation 
           j=j_canalcoord(icnl_,2,1,2)-par%tjmax(1)   ! possible pour le canal: icroissant
         endif                     !sender>
         if(send_or_recv_==2) then !receiver>
           i=i_canalcoord(icnl_,2,1,2)-par%timax(1)   ! note: le receiver est comme une C.L.
           j=j_canalcoord(icnl_,2,1,2)-par%tjmax(1)   ! donc tout au bout du canal
         endif                     !receiver>
        endif                 !gt1> 

        if(gridtype_==2) then !gt2>
         if(send_or_recv_==1) then !sender>
           i=i_canalcoord(icnl_,2,2,1)-par%timax(1)   ! note: les coord en mer donnees en fichier 
           j=j_canalcoord(icnl_,2,2,1)-par%tjmax(1)   ! d'entree sont celles du sender
         endif                     !sender>
         if(send_or_recv_==2) then !receiver> ! note: le decallage en indice du point
           i=i_canalcoord(icnl_,2,2,1)-par%timax(1) ! receiver depend de l'orientation du
           j=j_canalcoord(icnl_,2,2,1)-par%tjmax(1) ! flux entrant en mer:
           if(canaldir(icnl_,point_)==1)i=i-1 ! 1=i croissant
           if(canaldir(icnl_,point_)==3)i=i+1 ! 3=i decroissant
           if(canaldir(icnl_,point_)==2)j=j-1 ! 2=j croissant
           if(canaldir(icnl_,point_)==4)j=j+1 ! 4=j decroissant
         endif                     !receiver>
        endif                 !gt2> 

       endif              !22222>

        if(i>1.and.i<imax.and.j>1.and.j<jmax) then !ooo>
         k0=1 ; k1=par%rank
        endif                                      !ooo>
        call mpi_allreduce(k0,k2,1,mpi_integer,mpi_sum,par%comm2d,ierr)

        if(k2/=1) then !debug>
        write(10+par%rank,'(a)')'no rank found for'
        write(10+par%rank,'(8(a,i4))')'canal',icnl_,' point',point_ &
         ,' gridtype',gridtype_,' send_or_recv',send_or_recv_   &
        ,' rank=',canalrank(icnl_,point_,gridtype_,send_or_recv_) &
         ,' k2=',k2,' i=',i+par%timax(1),' j=',j+par%tjmax(1)
         if(i+par%timax(1)==1.or.i+par%timax(1)==iglb.or. &
            j+par%tjmax(1)==1.or.j+par%tjmax(1)==jglb)    &
         write(10+par%rank,'(3a)')'Can not change border land/sea mask ' &
         ,' => canal entrance shoud not be located at i=2 or iglb-1 '   &
         ,'or j=2 or jglb-1'

        flag_stop=1
        endif          !debug>

        call mpi_allreduce(k1,k2,1,mpi_integer,mpi_sum,par%comm2d,ierr)
        canalrank(icnl_,point_,gridtype_,send_or_recv_)=k2

        if(par%rank==0)write(6,'(5(a,i4))')'canal',icnl_,' point',point_ &
         ,' gridtype',gridtype_,' send_or_recv',send_or_recv_   &
        ,' rank',canalrank(icnl_,point_,gridtype_,send_or_recv_)

       endif                                                       !m°v°m>

      enddo ! send_or_recv_
      enddo ! gridtype_
      enddo ! point_ (1 ou 2)
      enddo ! icnl_


      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop 'stop grid_webcanals_rank. See fort.xxx err files'


! Les ranks receiver supplémentaires associes A la zone d'ombre mpi z0
! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
      if(.not.allocated(canalrankbis))  then !rank>
               allocate(canalrankbis(nbcanal,2,gridtype1:gridtype2,3)) ; canalrankbis=-999
      endif  

      do icnl_=1,nbcanal  
      do gridtype_=gridtype1,gridtype2
      do point_=1,2

! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype_,receiver)>=0) then !test de connexion a un de point de mer>

       irecv=i_canalcoord(icnl_,point_,gridtype_,receiver)-par%timax(1)
       jrecv=j_canalcoord(icnl_,point_,gridtype_,receiver)-par%tjmax(1) 

       loop_=0
       k0=-999 ! valeur par defaut pour pas de chevauchement
       if(jrecv==2.or.jrecv==jmax-1) then !>>>
         if(jrecv==2)       k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),4) ! voisin 'sud' du receiver
         if(jrecv==jmax-1)  k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),3) ! voisin 'nord' du receiver
       endif                              !>>>
       call mpi_allreduce(k0,k10,1,mpi_integer,mpi_max,par%comm2d ,ierr)
! k10 rank pour un echange au sud ou au nord
       if(k10>=0) then !>>>
         loop_=loop_+1
         canalrankbis(icnl_,point_,gridtype_,loop_)=k10
       endif           !>>>

       
       k0=-999 ! valeur par defaut pour pas de chevauchement
       if(irecv==2.or.irecv==imax-1) then !>>>
         if(irecv==2)       k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),1) ! voisin 'ouest' du receiver
         if(irecv==imax-1)  k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),2) ! voisin 'est' du receiver
       endif                              !>>>
       call mpi_allreduce(k0,k10,1,mpi_integer,mpi_max,par%comm2d ,ierr)
! k10 rank pour un echange a l'est ou a l'ouest
       if(k10>=0) then !>>>
         loop_=loop_+1
         canalrankbis(icnl_,point_,gridtype_,loop_)=k10
       endif           !>>>

! coins
       k0=-999
       if(irecv==2     .and.jrecv==2     )k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),5) ! voisin 'SO' du receiver
       if(irecv==2     .and.jrecv==jmax-1)k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),7) ! voisin 'NO' du receiver
       if(irecv==imax-1.and.jrecv==2     )k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),6) ! voisin 'SE' du receiver
       if(irecv==imax-1.and.jrecv==jmax-1)k0=par%gtvoisin(canalrank(icnl_,point_,gridtype_,receiver),8) ! voisin 'NE' du receiver
       call mpi_allreduce(k0,k10,1,mpi_integer,mpi_max,par%comm2d ,ierr)
       if(k10>=0) then !>>>
         loop_=loop_+1
         canalrankbis(icnl_,point_,gridtype_,loop_)=k10
       endif           !>>>

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! gridtype_
      enddo ! icnl_=1,nbcanal  

      end subroutine grid_webcanals_rank

!......................................................................

      subroutine grid_webcanals_consistency !18-06-19
      implicit none
      integer :: icnl_, point_

      if(ioffline==2) return !18-08-22
! cas ioffline 2 la grille lue dans grid.nc a deja recu des modifs qui ne passeraient pas les tests de verification de cette subroutine

      do icnl_=1,nbcanal   ! numero de canal
        do point_=point1,point2

         if(i_canalcoord(icnl_,point_,gridtype2,sender)>=1.and.i_canalcoord(icnl_,point_,gridtype2,sender)<=iglb.and. &
            j_canalcoord(icnl_,point_,gridtype2,sender)>=1.and.j_canalcoord(icnl_,point_,gridtype2,sender)<=jglb) then !>>>
            if(glob_mask(i_canalcoord(icnl_,point_,gridtype2,sender),j_canalcoord(icnl_,point_,gridtype2,sender))/=1) then !--->
             flag_stop=1
             write(10+par%rank,'(a)')' Err: mask/=1'
             write(10+par%rank,*)'icnl_=',icnl_
             write(10+par%rank,*)'point',point_  
             write(10+par%rank,*)'i,j',i_canalcoord(icnl_,point_,gridtype2,sender) &
                                      ,j_canalcoord(icnl_,point_,gridtype2,sender)  
             write(10+par%rank,*)'globmask ' &
                ,glob_mask(i_canalcoord(icnl_,point_,gridtype2,sender),j_canalcoord(icnl_,point_,gridtype2,sender))
            endif                                                  !--->
         endif                                                 !>>>

! On verifie d'autre part que ce point d'entree en mer est bien situE devant un point masque 
! et (optionnel) eventuellement on demasque ce point si interet algorithmique
! Point 1
         i=i_canalcoord(icnl_,point_,gridtype2,receiver) 
         j=j_canalcoord(icnl_,point_,gridtype2,receiver) ! point 1
         if(i>=1.and.i<=iglb.and.j>=1.and.j<=jglb) then !>>>
          if(glob_mask(i,j)/=0) then !pmx>
            flag_stop=1
            write(10+par%rank,*)'Err i,j,glob_mask/=0',i,j,glob_mask(i,j) &
            ,'. Channel entrance must be in front of a mask=0 point.' & 
            ,' Modify point',point_,' of channel ',icnl_
          endif                      !pmx>
         endif                                          !>>>

         enddo ! point_
        enddo  ! icnl_

        call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) !12-12-19
        if(k0/=0)stop 'Err grid_webcanals_consistency 2893'

      end subroutine grid_webcanals_consistency

!......................................................................

      end module module_grid
