










      subroutine moveforward
!______________________________________________________________________
! SYMPHONIE ocean model
! release 290 - last update: 01-10-20
!______________________________________________________________________
!    _________                    .__                  .__             !
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 

      use module_principal
      implicit none

!.......................................................................
! Version date      Description des modifications
!         15/09/06: mis en service
!         03/10/03: on rétablit le filtre temporelle sur T et S,
!                   qui avait temporairement été retiré le temps des tests
!                   academiques
!         08/06/07: bienvenue à hzdy_x et hzdx_y
!         13/05/08: TEM_Z et SAL_Z sont sauvés dans THZ_Z(2) et SHZ_Z(2)
!                   pour schema d'advection
!                   Taille des boucles ramenée à 1 MECO, 1,NECO
!         06-06-09  Boucles de calcul de hz_z etendues à i=0 j=0 pour parallelisation
! 2009.3  07-10-09  utilisation des nouveaux facteurs d'echelle verticale
!         10-10-09  tem et sal remplacent thz et shz
!         15-10-09  pour une bonne parallelisation il est preferable de n'updater
!                   que t et s "after" et ne pas toucher à t et s "before"
! 2010.7  15-02-10  "mise à jour en avant" des vitesses
!         16-02-10  move_ts.F90 renommé moveforward.F90
! 2010.8  10-03-10  sauvegarde de l'echeance "-1" pour vel, tem & sal
! 2010.9  01-06-10  update dz
! 2010.16 12-01-11  Filtre ordre élevé sur dz: ajout dz_t(-1)
! 2010.18 28-01-11  Possibilité de teste le vrai filtre d'Asselin
!         15-02-11  Asselin s'appelle avant la mise à jour
! S26     27-01-13  boucles etendues pour advection qdm o4
!         30-08-14  Lignes relatives au filtre FD sont commentees
!         28-11-14  suprression ssh_int_u ssh_int_v
!         16-12-15  boucles "f90"
!         28-01-16  suppression d'un stop
! v290    01-10-20  ajout de l'echeance -1 dans ssh_int_w
!.......................................................................

      vel_u(:,:,:,-1)=vel_u(:,:,:,0)
      vel_u(:,:,:, 0)=vel_u(:,:,:,1)
      vel_u(:,:,:, 1)=vel_u(:,:,:,2)

      velavr_u(:,:,0)=velavr_u(:,:,1)

      vel_v(:,:,:,-1)=vel_v(:,:,:,0)
      vel_v(:,:,:, 0)=vel_v(:,:,:,1)
      vel_v(:,:,:, 1)=vel_v(:,:,:,2)

      velavr_v(:,:,0)=velavr_v(:,:,1)

      tem_t(:,:,:,-1)=tem_t(:,:,:,0)
      tem_t(:,:,:, 0)=tem_t(:,:,:,1)
      tem_t(:,:,:, 1)=tem_t(:,:,:,2)

      sal_t(:,:,:,-1)=sal_t(:,:,:,0)
      sal_t(:,:,:, 0)=sal_t(:,:,:,1)
      sal_t(:,:,:, 1)=sal_t(:,:,:,2)

      ssh_int_w(:,:,-1)=ssh_int_w(:,:,0) !01-10-20
      ssh_int_w(:,:,0 )=ssh_int_w(:,:,1)
      ssh_int_w(:,:,1 )=ssh_int_w(:,:,2)

      hz_u(:,:,0)=hz_u(:,:,1)
      hz_u(:,:,1)=hz_u(:,:,2)

      hz_v(:,:,0)=hz_v(:,:,1)
      hz_v(:,:,1)=hz_v(:,:,2)

      hz_w(:,:,0)=hz_w(:,:,1)
      hz_w(:,:,1)=hz_w(:,:,2)

      dz_t(:,:,:,-1)=dz_t(:,:,:,0 ) ! Note: temps -1 necessaire A advection turbulence 
      dz_t(:,:,:,0 )=dz_t(:,:,:,1 )
      dz_t(:,:,:,1 )=dz_t(:,:,:,2 )

      dz_u(:,:,:,-1)=dz_u(:,:,:,0)
      dz_u(:,:,:,0 )=dz_u(:,:,:,1)
      dz_u(:,:,:,1 )=dz_u(:,:,:,2)

      dz_v(:,:,:,-1)=dz_v(:,:,:,0)
      dz_v(:,:,:,0 )=dz_v(:,:,:,1)
      dz_v(:,:,:,1 )=dz_v(:,:,:,2)

      end subroutine moveforward
