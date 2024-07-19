      subroutine strada(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 303 - last update: 24-07-21
!______________________________________________________________________
!    _________                    .__                  .__            !m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| _         !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \     ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/     !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >    !
!          \/\/          \/|__|        \/            \/        \/     !
!......................................................................

      use module_principal ; use module_biobc ; use module_biobalance
      use module_biology   ; use module_s ; use module_my_outputs
      use module_webcanals
      use module_parameter_sedim, only : l_sedim_s
      use module_sediment,        only : run_sediment , initial_sediment

      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='strada'
       subroutinedescription= &
       'Driver of the subroutines computing ' &
       //' the transport of passive and/or biogeochemical tracers'

       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Simulation of TRAcer Diffusion Advection

!......................................................................
! Version date      Description des modifications
!         14/12/01: bienvenue à ITIMEBIO
!         17/12/01: amenagements sur cas forward (economie de RAM)
!         26/12/02: on ne rappelle pas couple_modes et omega_upd. Les vitesses
!                   horizontales et verticale forward et leapfrog sont calculees
!                   une fois pour toute au premier appel de ces 2 routines
!         19/03/03: bienvenue à CALL SOURCE_TRACER (nouveau calcul des points
!                   d'émission des polluants)
!         25/03/03: choix entre 3 schemas d'advection pour la bio et la mes
!         11/07/05: nouveau schema d'advection tiré du schéma pour T et S mais
!                   version forward
!         06/04/06: choix entre fleuves "bio" ou "traceurs"
!                   Le seuil programme par Emilie est deplace en biohz_to_bio.F
!         07/04/06: Nouveaux appels a des conditions aux limites
!         22/06/06: cas model_ 1D bio
!         04/07/06: advection upstream
! 2009.3  02-10-09: commentaire devant sediment_bio & test sur I1D
!                   suppression dti
! 2010.8  03-05-10  suppression appels biohz_to_bio & asselin_bio
! S.26    28-03-13  call biobc_driver
!         27-04-13  call biobc_nudging_area
!         23-11-13  des appels a des diagnostiques
!         06-02-14  - parallelisation etat initial
!                   - check mpi
!         07-04-14  amelioration de la routine check mpi
!         01-07-14  amelioration de la routine check mpi
!         11-04-14  call nest_bio(1) supprime
!         14-04-14  advection_bioup remplace par advection_bio
!         16-07-14  gestion du stop dans check mpi
!         04-10-14  ajout flag_ dans checkmpi
!         25-04-16  call check_mpi_1dv_bio('A')
!         11-07-16  call advection_bio_rmnegval supprime les valeurs negatives
!         05-12-16  si flag_rmnegval=1 activer le retrait des valeurs negatives
!         13-12-16  Choix wsed implicite ou explicite selon les traceurs
!         18-01-17  lignes test supprimees
!         25-01-17  if(flag_nemoffline==1)vststep=1
!         12-02-17  Mises A jour cas NEMO offline
!         14-03-17  La concentration au point d'entree des fleuves est mise A jour
!                   juste avant le calcul de l'advection
! v249    04-03-19  subroutine strada_bio_minmax !04-03-19
!         05-03-19  warning si kmin_w/=1 et flag_merged_levels==0
! v252    22-04-19  ajout reseau de canaux
! v253    01-05-19  pas de porte entre advection et mixsedbio
!         07-05-19  subroutine strada_bio_minmax !04-03-19 !07-05-19
! v254    12-05-19  Echange canaux avant (puis apres) advection_bio_rmnegval3d_plus
! v303    19-07-21  if(modulo(iteration3d,modulo_biotimestep)==0)call biology !19-07-21        
!         24-07-21  call my_outputs_zone1bioflux('init',0) etc etc...
!...............................................................................

!......................................................................
! faire un initial_bio  (qui lira aussi fichier retsart)
! faire un close_bio (qui ecrira aussi fichier retsart)

!*******************************************************************************
!  EQUATIONS D'ADVECTION DIFFUSION DES TRACEURS (Biologie Sediments ...)
!  DEBUT:

      if(case_==1) then

!*******************************************************************************

      itime=itimebio                                                   !14/12/01

      call biobalance_gateA

!.............................................
! LES FORCAGES:
! Debut:
!.............................................

! calcul des flux eau / sediment
!     call sediment_bio                                                !02-10-09
! Mise a jour des forcages par les rivieres:
      if(imodelbio.eq.1)call river_bio_upd                             !06/04/06
      call s_cpu('rivers',0) !11-02-17

! calcul des sources de traceurs:
      if(ksomax.ge.1)call source_tracer                                !19/03/03
! Mise a jour des forcages aux limites ouvertes:
!     if(nest_onoff_in.eq.1)call nest_bio(2)
      if(biobc_extforcing==1)call biobc_driver                         !28-03-13
      call s_cpu('obc bio',0) !11-02-17

!.............................................
! LES FORCAGES:
! Fin.
!.............................................

! Faire des sorties pour des model_s imbriques:
!     if(nest_onoff_out.ge.1)call nest_bio(1) !11-07-14

      if(imodelbio==1) then !ECO3M-S>
           if(modulo(iteration3d,modulo_biotimestep)==0)call biology !19-07-21        
      endif                 !ECO3M-S>
      if(l_sedim_s)call run_sediment(1) !MUSTANG>
      call biobalance_gateB                                            !23-11-13
      call s_cpu('eco3m-mustang-gateB',0)

! Decommenter cette ligne pour verifier la parallelisation#
#ifdef checkmpi
      call strada_check_mpi_conservation                               !06-02-14
#endif

! Concentration at river input grid point (just before advection)
      if(nriver.ge.1)call obc_bio(1) ! river input  ! appel deplacE le !14-03-17
                     call biobalance_gateF                             !23-11-13

! advection:
      if(flag3d==1.and.flag3dbio==0)call advection_bio                       !14-07-14
! Pas de porte ici car bio_t est bio_t*dz_t !01-05-19
! flux turbulents & sédimentation
      call mixsed_bio
      call biobalance_gateD                                            !23-11-13
      call s_cpu('mixedbio+gateD',0)

! Echange canaux avant (puis apres) advection_bio_rmnegval3d_plus !12-05-19
      if(nbcanal>0) then !m°v°m> !12-05-19 
         call webcanals_gt1_to_gt2_bio
         call webcanals_gt2_to_gt1_bio
      endif              !m°v°m>

! Remove negative values (apres l'advection y compris la possible advection
! implicite dans mixsed_bio) !11-07-16
      if(flag_rmnegval==1)call advection_bio_rmnegval_vert
      if(flag_rmnegval==2)call advection_bio_rmnegval3d_plus !03-11-17

! Apply radioactive decay coefficient
      call halflife_radio

! conditions aux limites
      if(biobc_extforcing==1)call biobc_nudging_area ! Zone de rappel  !27-04-13
                             call biobalance_gateE                     !23-11-13

      if(nbcanal>0) then !m°v°m> !22-04-19
         call webcanals_gt1_to_gt2_bio
         call webcanals_gt2_to_gt1_bio
      endif              !m°v°m>
      call obc_bio(2)                ! Les conditions limites ouvertes !07/04/06

      if(l_sedim_s)call run_sediment(2) !MUSTANG !14-12-16

!     call strada_bio_minmax !04-03-19 !07-05-19
!     call my_outputs_bio_sum

#ifdef bilanbio
!!!!! call my_outputs_zone1bioflux('botsurf',0) ! deplace dans mixsedbio pour pouvoir deduire le flux de sedimentation implicite
      call my_outputs_zone1bioflux('tendancebio',0)
      if(modulo(iteration3d,10)==0) call my_outputs_zone1bioflux('mpi',0)
#endif


#ifdef checkmpi
      if(flag_1dv==1)call check_mpi_1dv_bio('A') !25-04-16
#endif


      call s_cpu('obc bio',0)

!*******************************************************************************
!  EQUATIONS D'ADVECTION DIFFUSION DES TRACEURS (Biologie Sediments ...)
!  FIN.
      return
      endif
!*******************************************************************************


!                                  /   /   /


!*******************************************************************************
!  INITIALISATION DES TRACEURS
!  DEBUT:
      if(case_==0) then
!*******************************************************************************

      call initial_tracer
      call initial_bio
      vbmax_eco3ms=vbmax
!     write(6,*)'vbmax_eco3ms=',vbmax_eco3ms
      call initial_sediment !13-12-16
!     write(6,*)'vbmax=',vbmax
!     stop 'coucou0'


!..............................................................
! Vitesse de sedimentation:
      if(vbmax>0)allocate(wsed_explicit(vbmax))        !13-12-16

      do vb=1,vbmax

! Explicite ou implicite ?:
       wsed_explicit(vb)=1 ! Par defaut schema explicite !13-12-16
!      wsed_explicit(vb)=0 ! implicite
       if(vb>vbmax_eco3ms)wsed_explicit(vb)=0 ! Eventuellement traceurs sedim implicit...

! Profil z:
       do k=2,kmax !12-02-17
        wsed(k,vb)=wsed(1,vb)
       enddo
       wsed(kmax+1,vb)=0.! Vitesse nulle A l'interface eau/air:

      enddo ! boucle sur vb
!..............................................................

      call obc_bio(0)  ! Reset pour OBC                                !07/04/06

!     call obc_bio_mpi('z0') !06-02-14
!     call obc_bio_mpi('z1') !06-02-14
      call obc_bio_mpi('zc') !14-07-14

!..............................
! Les requis du cas nemoffline:
      if(flag_nemoffline==1) then !pmxpmx> !12-02-17
         vststep=1 !25-01-17
         dz_t=max(dz_t,small1)
      endif
!..............................

      if(par%rank==0)write(6,*)'flag_nemoffline,vststep=',flag_nemoffline,vststep
      if(par%rank==0)write(6,*)'vbmax,vbmax_eco3ms=',vbmax,vbmax_eco3ms

      if(imodeltrc.eq.1.and.imodelbio.eq.1)then
       write(6,*)'les 2 1er boutons de notebook_tracer et notebook_bio'
       write(6,*)'ne peuvent être tous les 2 égaux à 1, il faut choisir'
       stop ' STOP donc dans strada.f!'
      endif

! Verifier que flag_merged_levels=1 si kmin_w/=1 dans le fichier de grille
! d'entre d'une simulation offline:
      if(flag_merged_levels==0) then !debug> !05-03-19
       flag_stop=0
       do j=1,jmax ; do i=1,imax
        if(kmin_w(i,j)/=1)flag_stop=1
       enddo       ; enddo
       call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       if(k0/=0) stop &
       'Err 238: use flag_merged_levels=1 (notebook_vertcoord)'
      endif                          !debug> !05-03-19

#ifdef bilanbio
      call my_outputs_zone1bioflux('init',0) !24-07-21
      if(flag_rmnegval==2) then !ooo> !01-08-21
! Attention le schema 3D supprimant les valeurs negatives n'est pas comptabilisE
! dans le bian, donc soit prendre le schema vertical (flag_rmnegval=1) ou aucun
! schema (flag_rmnegval=0)
       stop 'Err 280: choose flag_rmnegval=1 or flag_rmnegval=0'
      endif                     !ooo>
#endif

      if(imodelbio==1.and.ioffline==2) then !pmx> !01-08-21
       if(dti_fw>600.) then
        stop &
        'Warning: physical time step too large (dti_fw>600s detected)'
       endif
      endif                                 !pmx>

!*******************************************************************************
!  INITIALISATION DES TRACEURS
!  FIN.
      return
      endif
!*******************************************************************************


!                                  /   /   /


!*******************************************************************************
!  FERMER PROPREMENT LE MODELE DES TRACEURS
!  DEBUT:
      if(case_==2) then
!*******************************************************************************

      call close_bio

!*******************************************************************************
!  FERMER PROPREMENT LE MODELE DES TRACEURS
!  FIN.
      return
      endif
!*******************************************************************************

      end subroutine strada

!...............................................................................

      subroutine strada_check_mpi_conservation                    !06-02-14
      use module_principal
      use module_parallele
      implicit none
      integer idi_anyv3d_z0,loop_,flag_
#ifdef synopsis
       subroutinetitle='strada_check_mpi_conservation'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      flag_=0 !04-10-14

      do vb=1,vbmax

! Pour verifier la conservation de la parallelisation:
! ENCORE PLUS STRICTE!! Pour etre absolument certains que la procedure
! de verification n'intervient pas dans la solution, on n'echange pas bio_t
! mais la copie double de bio_t, a savoir anyv3d(i,j,k,2)....
      do k1=1,2
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
         anyv3d(i,j,k,k1)=bio_t(i,j,k,vb)
! Petit test pour verifier qu'une erreur est bien detectee: le code doit s'arreter si ligne
! suivante decommentee:
!        anyv3d(i,j,k,k1)=bio_t(i,j,k,vb)+0.00001*real(par%rank)
      enddo
      enddo
      enddo
      enddo

#ifdef parallele
      call get_type_echange('z0','anyv3d_z0',anyv3d,lbound(anyv3d),ubound(anyv3d),2,idi_anyv3d_z0) ! 'z0' = 'z0' et 'z1'
      ! Echanges
      do loop_=1,subcycle_exchange
      call echange_voisin(anyv3d,idi_anyv3d_z0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      if(par%tvoisin(sud)/=mpi_proc_null) then   !ssssssss> !01-07-14
      j1=1
      do k=1,kmax
      do i=1,imax
       if(anyv3d(i,j1,k,1)/=anyv3d(i,j1,k,2))then
        write(10+par%rank,*)'----'
!       write(10+par%rank,*)'mpi SUD non conservé en vb k=',vb,k
        write(10+par%rank,*)'mpi SUD non conservé en'
        write(10+par%rank,*)'vb =',vb
        write(10+par%rank,*)'i,j,k',i,j1,k
        write(10+par%rank,*)'du proc',par%rank
        write(10+par%rank,*)'proc voisin',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*)'anyv3d',anyv3d(i,j1,k,1),anyv3d(i,j1,k,2)
        write(10+par%rank,*)'delta ',anyv3d(i,j1,k,1)-anyv3d(i,j1,k,2)
        write(10+par%rank,*)'mask  ',mask_t(i,j1,kmax)
        write(10+par%rank,*)'kmin_w,kmerged_t',kmin_w(i,j1),kmerged_t(i,j1)
        if(k>kmerged_t(i,j1)) then !ooo>
        write(10+par%rank,*)'k >kmerged'
        else                       !ooo>
        write(10+par%rank,*)'k<=kmerged'
        endif                      !ooo>
        flag_=1
       endif
      enddo
      enddo
      endif                                      !ssssssss>

      if(par%tvoisin(nord)/=mpi_proc_null) then  !nnnnnnnn>
      j2=jmax
      do k=1,kmax
      do i=1,imax
       if(anyv3d(i,j2,k,1)/=anyv3d(i,j2,k,2))then
        write(10+par%rank,*)'----'
        write(10+par%rank,*)'mpi NORD non conservé en'
        write(10+par%rank,*)'vb   ',vb
        write(10+par%rank,*)'i,j,k',i,j2,k
        write(10+par%rank,*)'proc',par%rank
        write(10+par%rank,*)'proc voisin',par%tvoisin(nord)
        write(10+par%rank,*)'coordonnées',i,2
        write(10+par%rank,*)'anyv3d',anyv3d(i,j2,k,1),anyv3d(i,j2,k,2)     
        write(10+par%rank,*)'delta ',anyv3d(i,j2,k,1)-anyv3d(i,j2,k,2)       
        write(10+par%rank,*)'mask  ',mask_t(i,j2,kmax)
        write(10+par%rank,*)'kmin_w,kmerged_t',kmin_w(i,j2),kmerged_t(i,j2)
        if(k>kmerged_t(i,j2)) then !ooo>
        write(10+par%rank,*)'k >kmerged'
        else                       !ooo>
        write(10+par%rank,*)'k<=kmerged'
        endif                      !ooo>
        flag_=1
       endif
      enddo
      enddo
      endif                                      !nnnnnnnn>


      if(par%tvoisin(ouest)/=mpi_proc_null) then !oooooooo>
      i1=1
      do k=1,kmax
      do j=1,jmax
       if(anyv3d(i1,j,k,1)/=anyv3d(i1,j,k,2)) then
        write(10+par%rank,*)'----'
        write(10+par%rank,*)'mpi OUEST non conservé en'
        write(10+par%rank,*)'vb  ',vb
        write(10+par%rank,*)'i,j,k',i1,j,k
        write(10+par%rank,*)'proc',par%rank
        write(10+par%rank,*)'proc voisin',par%tvoisin(ouest)

        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*)'anyv3d',anyv3d(i1,j,k,1),anyv3d(i1,j,k,2)
        write(10+par%rank,*)'delta ',anyv3d(i1,j,k,1)-anyv3d(i1,j,k,2)  
        write(10+par%rank,*)'mask  ',mask_t(i1,j,kmax)
        write(10+par%rank,*)'kmin_w,kmerged_t',kmin_w(i1,j),kmerged_t(i1,j)
        if(k>kmerged_t(i1,j)) then !ooo>
        write(10+par%rank,*)'k >kmerged'
        else                       !ooo>
        write(10+par%rank,*)'k<=kmerged'
        endif                      !ooo>
        flag_=1
       endif
      enddo
      enddo
      endif                                      !oooooooo>


      if(par%tvoisin(est)/=mpi_proc_null) then   !eeeeeeee>
      i2=imax
      do k=1,kmax
      do j=1,jmax
       if(anyv3d(i2,j,k,1)/=anyv3d(i2,j,k,2)) then
        write(10+par%rank,*)'----'
        write(10+par%rank,*)'mpi EST non conservé en'
        write(10+par%rank,*)'vb   ',vb
        write(10+par%rank,*)'i,j,k',i2,j,k
        write(10+par%rank,*)'proc',par%rank
        write(10+par%rank,*)'proc voisin',par%tvoisin(est)
        write(10+par%rank,*)'coordonnées',2,j
        write(10+par%rank,*)'anyv3d',anyv3d(i2,j,k,1),anyv3d(i2,j,k,2)
        write(10+par%rank,*)'delta ',anyv3d(i2,j,k,1)-anyv3d(i2,j,k,2)
        write(10+par%rank,*)'mask  ',mask_t(i2,j,kmax)
        write(10+par%rank,*)'kmin_w,kmerged_t',kmin_w(i2,j),kmerged_t(i2,j)
        if(k>kmerged_t(i2,j)) then !ooo>
        write(10+par%rank,*)'k >kmerged'
        else                       !ooo>
        write(10+par%rank,*)'k<=kmerged'
        endif                      !ooo>
        flag_=1
       endif
      enddo
      enddo
      endif                                      !eeeeeeee>

!     k0=0                                       !16-07-14
!     call mpi_allreduce(flag_,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
!     if(k0>0) then
!       write(10+par%rank,*)'par%rank & nb erreur mpi=',par%rank,ksecu
!       stop 'dans strada check mpi'
!     endif

      enddo ! vb

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(flag_==1)stop ' Stop strada_check_mpi_conservation'

      end subroutine strada_check_mpi_conservation

!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine strada_bio_minmax !04-03-19 !07-05-19
      use module_principal ; use module_parallele
      implicit none

       do vb=1,vbmax

       x1=9999. ; x2=-x1
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
       if(mask_t(i,j,kmax)==1) then

        if(bio_t(i,j,k,vb)<x1) then
         x1=bio_t(i,j,k,vb) ; i1=i ; j1=j ; k1=k
        endif
        if(bio_t(i,j,k,vb)>x2) then
         x2=bio_t(i,j,k,vb) ; i2=i ; j2=j ; k2=k
        endif

       endif
       enddo ; enddo ; enddo

! Min et Max de T et S par proc ecrits dans des fichiers:
!      write(10+par%rank,*)iteration3d
!      write(10+par%rank,*)i1,j1,k1,x1
!      write(10+par%rank,*)i2,j2,k2,x2
!      write(10+par%rank,*)i3,j3,k3,x3
!      write(10+par%rank,*)i4,j4,k4,x4

      call mpi_allreduce(x1,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN BIO=',x0
      if(x1==x0) then !>>>
           write(6,*)'MIN BIO iteration ------------',iteration3d
           write(6,*)'MIN BIO vb',vb
           write(6,*)'MIN BIO val',x0 !,bio_t(i1,j1,k1,vb)
           write(6,*)'MIN BIO rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN BIO glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN BIO dz_t',dz_t(i1,j1,k1,2)
           write(6,*)'MIN BIO h_w',h_w(i1,j1)
      endif           !>>>

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX BIO=',x0
      if(x2==x0) then !>>>
           write(6,*)'MAX BIO iteration ------------',iteration3d
           write(6,*)'MAX BIO vb',vb
           write(6,*)'MAX BIO val',x0 !,bio_t(i2,j2,k2,vb)
           write(6,*)'MAX BIO rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX BIO glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX BIO dz_t',dz_t(i2,j2,k2,2)
           write(6,*)'MAX BIO h_w',h_w(i2,j2)
      endif           !>>>

      enddo ! boucle vb

      end subroutine strada_bio_minmax
#endif

!-----------------------------------------------------------------------------
