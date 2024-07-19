      subroutine add_bi(ki1,ki2,ki3)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 345 - last update: 07-04-22
!______________________________________________________________________
      use module_principal
      implicit none
      integer ki1,ki2,ki3
!...............................................................................
! modifs 26/03/07: mise en service
!        03/12/07: pas possible de mettre BI sur ztaobc(0) et zaobc(2)
!        01-06-09: Ajouter BI à ztaobc_z(1) puis rappeler initial_with_obc(1)
! v345   07-04-22  Si on initialise avec la SSH instantanee de symphonie en mode simple
!                  restart alors le BI est dEjA dans le signal...
!...............................................................................
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


! Ce programme ajoute aux champs de forcages et aux champs initiaux
! le barometre inverse

!_________________________________________________________________
! ATTENTION ON N'AJOUTE PAS LE BAROMETRE INVERSE SI ON EST DANS LE
! CAS D'UNE IMBRICATION DE SYMPHONIE DANS SYMPHONIE car l'effet
! est déjà inclus dans le model_ en amont
      if(nest_onoff_in.eq.1)return
      if(bi_onoff.eq.0)return

! Si on initialise avec la SSH instantanee de symphonie en mode simple
! restart alors le BI est dEjA dans les tableaux ssh_w et sshobc_w.
! Comme il ne faut pas que sshobc(0) et sshobc(2) contiennent le BI
! (car update_obcforcing ajoutera ensuite le BI A sshobc) on enleve
! le BI de sshobc(0) et sshobc(2)  
      if(simple_restart_file_txt/='none') then !pmx> !07-04-22
       const1=-1./grav/rho
       do j=1,jmax ; do i=1,imax

       sshobc_w(i,j,0)=(   sshobc_w(i,j,0)                              &
                     -const1*(pss_w(i,j,1)-pss_mean(1))                 &
                         )*  mask_t(i,j,kmax+1)
       sshobc_w(i,j,2)=(   sshobc_w(i,j,2)                              &
                     -const1*(pss_w(i,j,1)-pss_mean(1))                 &
                         )*  mask_t(i,j,kmax+1)

       enddo ; enddo
! Ensuite il ne faut pas ajouter le BI A sshobc(1) puisqu'il y est
! dEjA, donc on ne va pas A l'Etape suivante, donc return:
       return 
      endif                                    !pmx>

! Si on arrive ici c'est que simple_restart_file_txt=='none'

      const1=-1./grav/rho
      do j=1,jmax
      do i=1,imax

       sshobc_w(i,j,1)=(   sshobc_w(i,j,1)                              &
                     +const1*(pss_w(i,j,1)-pss_mean(1))                 &
                         )*  mask_t(i,j,kmax+1)

      enddo
      enddo

! Repercuter sur les autres variables:
      call initial_with_obc(1)

      end subroutine add_bi
