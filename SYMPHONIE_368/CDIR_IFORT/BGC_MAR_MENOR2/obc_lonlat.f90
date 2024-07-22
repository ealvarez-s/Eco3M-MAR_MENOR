










      subroutine obc_lonlat(ichoix)
!______________________________________________________________________
! SYMPHONIE ocean model
! release S.26 - last update: 03-10-18
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      implicit none
      integer ichoix,id_lon_t_zc,id_lat_t_zc,loop_
!...............................................................................
! Version date      Description des modifications
!         02/04/09: ajout echange de type Z2
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
! 2010.6  05-02-10: renomme lon_t lat_t
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S.26    22-06-14  nouveaux echanges
!         17-04-15  ajout obc_lonlat_z2 (pour cas nemo offline)
!         03-10-18  ajout obc_lonlat_extrapol 
!...............................................................................
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

! Pour maintenir la compatibilite avec les version precedentes,
! par defaut, tout appel a obc_lonlat appelle obc_lonlat_mpi_zc
      call obc_lonlat_mpi_zc

! cas particulier des domaines academique periodique pour lesquels on n'applique pas 
! la periodicite aux champs lon,lat
      if(discard_lonlat_periodicity==1)call obc_lonlat_extrapol !03-10-18

      end subroutine obc_lonlat

!...............................................................................

      subroutine obc_lonlat_mpi_zc !17-04-15
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer ichoix,id_lon_t_zc,id_lat_t_zc,loop_

      call get_type_echange('zc','lon_t_zc',lon_t     &
                                    ,lbound(lon_t)    &
                                    ,ubound(lon_t)    &
                                    ,id_lon_t_zc)
      call get_type_echange('zc','lat_t_zc',lat_t     &
                                    ,lbound(lat_t)    &
                                    ,ubound(lat_t)    &
                                    ,id_lat_t_zc)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(lon_t,id_lon_t_zc,mpi_neighbor_list(loop_))
        call echange_voisin(lat_t,id_lat_t_zc,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_lonlat_mpi_zc

!...............................................................................

      subroutine obc_lonlat_z2 !17-04-15
      use module_principal
      implicit none

      lon_t(-1    ,:     )=lon_t(0     ,:     )
      lon_t(imax+2,:     )=lon_t(imax+1,:     )
      lon_t(:     ,-1    )=lon_t(:     ,0     )
      lon_t(:     ,jmax+2)=lon_t(:     ,jmax+1)

      lat_t(-1    ,:     )=lat_t(0     ,:     )
      lat_t(imax+2,:     )=lat_t(imax+1,:     )
      lat_t(:     ,-1    )=lat_t(:     ,0     )
      lat_t(:     ,jmax+2)=lat_t(:     ,jmax+1)

! mpi continuity:
      call obc_lonlat_mpi_z2 !17-04-15

      end subroutine obc_lonlat_z2

!...............................................................................

      subroutine obc_lonlat_mpi_z2 !17-04-15
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer ichoix,id_lon_t_z2,id_lat_t_z2,loop_

      call get_type_echange('z2','lon_t_z2',lon_t     &
                                    ,lbound(lon_t)    &
                                    ,ubound(lon_t)    &
                                    ,id_lon_t_z2)
      call get_type_echange('z2','lat_t_z2',lat_t     &
                                    ,lbound(lat_t)    &
                                    ,ubound(lat_t)    &
                                    ,id_lat_t_z2)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(lon_t,id_lon_t_z2,mpi_neighbor_list(loop_))
        call echange_voisin(lat_t,id_lat_t_z2,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_lonlat_mpi_z2

!...............................................................................

      subroutine obc_lonlat_extrapol !03-10-18
      use module_principal
      use module_parallele !#MPI
      implicit none

! Extrapoler les valeurs de lon,lat en dehors du domaine 
! Cette condition au limite peut etre appelee dans le cas d'un
! domaine academique periodique pour lequel la periodicitE de lon et lat
! est illogique car la periodicite n'est pas atteinte sur la sphere terrestrre
! Dans ce cas la condition de periodicite sur lon lat ne doit pas s'appliquer 
! ou (si on l'a applique pour garantir la continuite mpi entre domaines) on
! peut appeler cette routine juste apres pour retablir des valeurs coherentes
! sur les bords du domaine. Ceci suppose de retablir les valeurs en i=1, i=imax etc....
! detruite par la condition mpi qui a precedE

       do i=-1,imax+2

       if(1+par%tjmax(1)==1) then !ooo>
        x1=lon_t(i,2)-lon_t(i,3)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
                    lon_t(i,1 )=lon_t(i,2)+x1
                    lon_t(i,0 )=lon_t(i,2)+x1*2
                    lon_t(i,-1)=lon_t(i,2)+x1*3

        x1=lat_t(i,2)-lat_t(i,3)
                    lat_t(i,1 )=lat_t(i,2)+x1
                    lat_t(i,0 )=lat_t(i,2)+x1*2
                    lat_t(i,-1)=lat_t(i,2)+x1*3
       endif                      !ooo>

       if(jmax+par%tjmax(1)==jglb) then !ppp>
        x1=lon_t(i,jmax-1)-lon_t(i,jmax-2)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
                         lon_t(i,jmax  )=lon_t(i,jmax-1)+x1  
                         lon_t(i,jmax+1)=lon_t(i,jmax-1)+x1*2
                         lon_t(i,jmax+2)=lon_t(i,jmax-1)+x1*3  

        x1=lat_t(i,jmax-1)-lat_t(i,jmax-2)
                         lat_t(i,jmax  )=lat_t(i,jmax-1)+x1
                         lat_t(i,jmax+1)=lat_t(i,jmax-1)+x1*2
                         lat_t(i,jmax+2)=lat_t(i,jmax-1)+x1*3
       endif                            !ppp>

       enddo

       do j=-1,jmax+2 ! Cette boucle est sur les bornes max (et non jstr_ et jend_).

       if(1+par%timax(1)==1) then !mmm>
        x1=lon_t(2,j)-lon_t(3,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
                    lon_t(1 ,j)=lon_t(2,j)+x1  
                    lon_t(0 ,j)=lon_t(2,j)+x1*2
                    lon_t(-1,j)=lon_t(2,j)+x1*3 

        x1=lat_t(2,j)-lat_t(3,j)
                    lat_t(1 ,j)=lat_t(2,j)+x1  
                    lat_t(0 ,j)=lat_t(2,j)+x1*2
                    lat_t(-1,j)=lat_t(2,j)+x1*3
       endif                      !mmm>

       if(imax+par%timax(1)==iglb) then !xxx>
        x1=lon_t(imax-1,j)-lon_t(imax-2,j)
        if(x1<-pi)x1=x1+2.*pi
        if(x1> pi)x1=x1-2.*pi
                         lon_t(imax  ,j)=lon_t(imax-1,j)+x1  
                         lon_t(imax+1,j)=lon_t(imax-1,j)+x1*2
                         lon_t(imax+2,j)=lon_t(imax-1,j)+x1*3

        x1=lat_t(imax-1,j)-lat_t(imax-2,j)
                         lat_t(imax  ,j)=lat_t(imax-1,j)+x1  
                         lat_t(imax+1,j)=lat_t(imax-1,j)+x1*2
                         lat_t(imax+2,j)=lat_t(imax-1,j)+x1*3
       endif                            !xxx>

       enddo

      end subroutine obc_lonlat_extrapol
