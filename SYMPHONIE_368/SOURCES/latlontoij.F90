      subroutine latlontoij(lon1_,lat1_,txt_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 341 - last update: 13-06-18
!______________________________________________________________________
      use module_principal
      use module_parallele !#mpi
      use module_grid
      implicit none
      double precision,intent(in) :: lon1_,lat1_
      double precision lon2_,lat2_
      character*3 txt_
#ifdef synopsis
       subroutinetitle='latlontoij'
       subroutinedescription='Converts lon lat in grid indexes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!__________________________________________________________________________
! Version date      Description des modifications
!         18/04/06: fonctions compatibles avec double precision
!         07/05/09: cas grille dlat dlon constant
! 2010.3  08-01-10: le pole de la grille n'est pas forcement le pole nord
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.13 21-10-10  choix entre indices globaux ou locaux
! 2010.25 04-03-12  separation (lon1,lat1) et (lon2,lat2)
! S26     31-01-13  lonlat2distance is a routine that computes the geodesic
!                   distance between two points
!         02-06-13  pole_lon pole_lat renommes
!         05-06-13  cas de la grille bipole
!         10-09-13  debug
!         29-12-13  cas grille horizontale generalisee
!         12-12-14  message d'erreur
!         13-06-18  blinder la fonction acos
! v341    27-03-22  debug faille latlon2ij
!...............................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................



!___________________________________________________________________________
!     DONNER LES INDICES DE GRILLE CORRESPONDANT A UN COUPLE
!     LATITUDE, LONGITUDE  (LATIT1, LONGI1) (en radians)
!
!     Si txt_ (en argument) = 'glb' alors position sur la grille globale
!     Si txt_ (en argument) = 'loc' alors position sur la grille locale (celle du process)
!___________________________________________________________________________


! PROJECTION MERCARTOR :
! la projection mercator provient d'Haltiner et Williams:
      if(typegrid==typegrid_monopole) then !111111111111111111>       !07/05/09

!.....
! Rotation inverse de la grille !08-01-10
      const1=-(90.-northpole_lat)*deg2rad
      const2=-( 0.-northpole_lon)*deg2rad
      x0=rayonterre*cos(lon1_  )*cos(lat1_  )
      y0=rayonterre*sin(lon1_  )*cos(lat1_  )
      z1=rayonterre*sin(lat1_  )
      x1= cos(const2)*x0+sin(const2)*y0
      y2=-sin(const2)*x0+cos(const2)*y0
      x2= cos(const1)*x1+sin(const1)*z1
      z2=-sin(const1)*x1+cos(const1)*z1
      lon2_  =atan2(y2,x2)
      lat2_  =asin(z2/rayonterre)
!.....

      decj=(                                                          &
          rayonterre*cos(phi0)*log( (1.+sin(lat2_  ))/cos(lat2_  ) )  &
         -gridcard(9,1) )/dyb + gridcard(6,1)

      deci=(                                                 &
          rayonterre*cos(phi0)*lon2_                         &
         -gridcard(8,1) )/dxb + gridcard(5,1)

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

      endif                  !1111111111111111111>

!.........

! GRILLE DLON DLAT CONSTANTS:
      if(typegrid.eq.2) then !222222222222222222>                     !07/05/09

       deci=1.+real(imax-1)*(lon1_         - lon_t(1,1))                &
                           /( lon_t(imax,1)- lon_t(1,1))
       decj=1.+real(jmax-1)*(lat1_         - lat_t(1,1))                &
                           /( lat_t(1,jmax)- lat_t(1,1))


      endif                  !222222222222222222>       !07/05/09

      if(typegrid==typegrid_bipole)call grid_bipole_lonlat2ij(lon1_,lat1_)

      if(typegrid==typegrid_file)call grid_lonlat2ij(lon1_,lat1_,txt_) !29-12-13

      if(txt_=='glb')return                   !21-10-10
      if(txt_=='loc')then                     !21-10-10
       deci=deci-par%timax(1)
       decj=decj-par%tjmax(1)
       return
      endif
      write(6,*)'bad choice in latlontoij!'     !12-12-14
      stop ' STOP in subroutine latlontoij'     !21-10-10

      end subroutine latlontoij

!......................................................................

      subroutine lonlat2distance(lon1_,lat1_,lon2_,lat2_,dist_) !31-01-13
      use module_principal
      implicit none
      double precision lon1_,lat1_,lon2_,lat2_,dist_
#ifdef synopsis
       subroutinetitle='lonlat2distance'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Given the longitude and latitude in radians of two points this routine
! return dist_ the geodesic distance

! De petites erreurs liEe A la prEcision du codage des nombres pouvant conduire A des
! valeurs hors gamme pour la fonction ACOS (en pratique de part et d'autre de -1, +1) l'algo
! est en consequence blindE

      dist_=rayonterre*acos( &
             min(max(sin(lat1_)*sin(lat2_)+cos(lat1_)*cos(lat2_)*cos(lon2_-lon1_),-1.),1.) & !13-06-18
                           )


      end subroutine lonlat2distance
