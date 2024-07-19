      subroutine obc_scal(case_)                                 !01-06-09
!______________________________________________________________________
! SYMPHONIE ocean model
! release 365 - last update: 26-01-23
!______________________________________________________________________
!Version date      Description des modifications
!        20/07/01: passage à la coordonnée sigma généralisée
!        25/11/01: prise en compte des petits fleuves
!        16/12/02: choix entre 3 schemas de conditions aux limites:
!        18/12/02: rappel vers climato sur choix 2
!        28/01/03: par defaut pas de rappel (const3=1.)
!        23/03/04: par defaut KI1=3
!        05/05/04: temperature des fleuves evolutive
!        11/05/04: debug thz vers tem et shz vers sal dans choix 3
!        18/05/04: passage des termes de gradient au temps t+1
!        28/05/04: choix entre leapfrog et forward pour ki1=3 (bienvenue à KI0)
!        29/05/04: passage de KI0, KI1, KI2 en arguments
!        11/06/04: modif sur condition choix1 c.a.d. condition upwind. T et S
!                  au point voisin dehors est évalué à partir du principe du
!                  "chemin parcouru"
!        29/06/04: finalisation cas precedent passe en KI1=4
!        01/08/05: Methode uptream enrichie de l'option masse d'eau repoussée
!                  quand le courant forcant est sortant. Cas 5.
!        15/09/05: Cas 8 pour le moment retenu. Mais ca bouge.....
!        30/09/05: Cas 5 et 7 commentés pour compiler avec -r8
!        31/05/06: A priori pas la peine de calculer 2 fois les coins....
!                  ... enfin c'est pas  encore sur...
!        21/04/07: seuls cas 1 et cas 8 retenus
!        17/01/08: nouvelles OBC: ajout d'un limiteur de contraste au cas
!                  numero 1
!        21/01/08: ajout d'un cas numero 3 couplé au cas numero 2
!        30/01/08: ajout d'un cas numero 4 couplé au cas numero 2
!        01/02/08: cas 3 modifié: perturbation proportionnelle à la
!                  varabilité de l'ogcm
!        04/03/08: cas 5 modifié: gradient symphonie = gradient ogcm à condition
!                  de ne pas trop s'ecarter de l'ogcm. L'ecart max est indexé sur
!                  la variabilité de l'ogcm.
!        13/05/08: Preparer des supplements de C.L. pour les couronnes en i=-1
!                  i=imax+2, j=-1, j=jmax+2 pour advection ordre élevé.
!        10/11/08: Amenagement de la C.L. pour tenir compte de ce que le
!                  sous-programme advection_ts.F ne fait pas de distinction
!                  des points interieurs et limites. C'est maintenant la C.L.
!                  qui fait la distinction de regime (entrant ou sortant).
!                  D'autre part la C.L. sur la deuxieme couronne peripherique,
!                  requise par le schema d'advection, est maintenant traitée dans
!                  la presente routine.
!        14/11/08: Ajout de la condition de Neumann (sur anomalies)
!        19/12/08: Par defaut X3=0 dans la condition ichoix=0
!        01-06-09: - Passage en arg de obctype_ts pour faire un appel avec un arg
!                  different dans initial_with_obc
!                  - Cas 1D: boucles etendues de -1 à imax+2, C.L. portant sur
!                  tem_c et non sur thz_c
!                  - Parallelisation
!        05-06-09: Attention à bien passer dans le cas 1DV...
! 2009.3 05-10-09: ajout d'un "ifdef parallele"
!        12-10-09: apres les amenagements sur thz et tem, il faut specifier ici
!                  l'update de tem_c(0) ...
!        16-10-09: echange T et S à t="now"
! 2010.11 16-07-10  temobc & salobc renommés temobc & salobc
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S26     27-01-13  C.L. au temps -1
!         12-05-14  nouveaux echanges
!         02-07-14  nouveaux echanges: subroutine obc_scal_mpi_z0_now
!         15-07-14  correction bug obc_scal_temobc & salobc
!         17-12-14  ajout obc_scal_rhp
!         14-10-15  ajout obc_scal_salref et obc_scal_temref
!         16-12-15  Calcul selon status obc
!         17-12-15  call obc_scal_moveforward commentE depuis extension boucle dans moveforward
!         19-03-17  subroutine obc_scal_t_jflux_cumul et cie calcule des bilan de flux
!                   "exacts" aux frontieres ouvertes (appels dans advection_scal)
!         21-05-18  ajout subroutine obc_scal_bottom_rhp(txt_) 
!         02-06-18  subroutine obc_scal_mpi_upwindriver(txt_)
! v299    18-03-21  utiliser obcstatus
! v309    01-10-21  mise A jour subroutine obc_scal_rhp pour echange zb
! v365    26-01-23  - obc_scal: Passage au temps "after" 
!                   - echange z2 sur rhp_t remplace echange zb 
!...............................................................................
!    _________                    .__                  .__                     !m[°v°]m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................
      use module_principal
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='obc_scal'
       subroutinedescription= &
       'Driver of the subroutines computing the lateral boundary ' &
       //'conditions for T and S'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Driver subroutine for T & S Open Boundary Conditions

      if(flag3d==0) then !------->

! cas du model_ 1D (I1D=0)
       call obc_scal_1dv

      else            !------->

!      call obc_scal_moveforward ! commentE depuis !17-12-15
       call obc_scal_wavecondition
!      call obc_scal_mpi_now

      endif           !------->

      end subroutine obc_scal

!...........................................................................

      subroutine obc_scal_moveforward
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='obc_scal_moveforward'
       subroutinedescription='Time update of T and S at OBC grid points'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(obcstatus(jeq1)==1) then !16-12-15>
       do k=1,kmax
       do i=-1,imax+2

! Update des obc au temps "beforebefore" avant calcul de celles au temps "before"
         sal_t(i,0     ,k,-1)=sal_t(i,0     ,k,0)
         sal_t(i,-1    ,k,-1)=sal_t(i,-1    ,k,0)
         tem_t(i,0     ,k,-1)=tem_t(i,0     ,k,0)
         tem_t(i,-1    ,k,-1)=tem_t(i,-1    ,k,0)

! Update des obc au temps "before" avant calcul de celles au temps "now": !12-10-09
         sal_t(i,0     ,k,0)=sal_t(i,0     ,k,1)
         sal_t(i,-1    ,k,0)=sal_t(i,-1    ,k,1)
         tem_t(i,0     ,k,0)=tem_t(i,0     ,k,1)
         tem_t(i,-1    ,k,0)=tem_t(i,-1    ,k,1)

       enddo ! i
       enddo  ! k
      endif                       !16-12-15>

      if(obcstatus(ieq1)==1) then !16-12-15>
       do k=1,kmax
       do j=1,jmax

! Update des obc au temps "before" avant calcul de celles au temps "now": !12-10-09
         sal_t(0     ,j,k,-1)=sal_t(0     ,j,k,0)
         sal_t(-1    ,j,k,-1)=sal_t(-1    ,j,k,0)
         tem_t(0     ,j,k,-1)=tem_t(0     ,j,k,0)
         tem_t(-1    ,j,k,-1)=tem_t(-1    ,j,k,0)

! Update des obc au temps "before" avant calcul de celles au temps "now": !12-10-09
         sal_t(0     ,j,k,0)=sal_t(0     ,j,k,1)
         sal_t(-1    ,j,k,0)=sal_t(-1    ,j,k,1)
         tem_t(0     ,j,k,0)=tem_t(0     ,j,k,1)
         tem_t(-1    ,j,k,0)=tem_t(-1    ,j,k,1)

       enddo ! j
       enddo  ! k
      endif                       !16-12-15>

      if(obcstatus(jeqjmax)==1) then !16-12-15>
       do k=1,kmax
       do i=-1,imax+2

! Update des obc au temps "beforebefore" avant calcul de celles au temps "before"
         sal_t(i,jmax+1,k,-1)=sal_t(i,jmax+1,k,0)
         sal_t(i,jmax+2,k,-1)=sal_t(i,jmax+2,k,0)
         tem_t(i,jmax+1,k,-1)=tem_t(i,jmax+1,k,0)
         tem_t(i,jmax+2,k,-1)=tem_t(i,jmax+2,k,0)

! Update des obc au temps "before" avant calcul de celles au temps "now": !12-10-09
         sal_t(i,jmax+1,k,0)=sal_t(i,jmax+1,k,1)
         sal_t(i,jmax+2,k,0)=sal_t(i,jmax+2,k,1)
         tem_t(i,jmax+1,k,0)=tem_t(i,jmax+1,k,1)
         tem_t(i,jmax+2,k,0)=tem_t(i,jmax+2,k,1)

       enddo ! i
       enddo  ! k
      endif                          !16-12-15>

      if(obcstatus(ieqimax)==1) then !16-12-15>
      do k=1,kmax
       do j=1,jmax

! Update des obc au temps "before" avant calcul de celles au temps "now": !12-10-09
         sal_t(imax+1,j,k,-1)=sal_t(imax+1,j,k,0)
         sal_t(imax+2,j,k,-1)=sal_t(imax+2,j,k,0)
         tem_t(imax+1,j,k,-1)=tem_t(imax+1,j,k,0)
         tem_t(imax+2,j,k,-1)=tem_t(imax+2,j,k,0)

! Update des obc au temps "before" avant calcul de celles au temps "now": !12-10-09
         sal_t(imax+1,j,k,0)=sal_t(imax+1,j,k,1)
         sal_t(imax+2,j,k,0)=sal_t(imax+2,j,k,1)
         tem_t(imax+1,j,k,0)=tem_t(imax+1,j,k,1)
         tem_t(imax+2,j,k,0)=tem_t(imax+2,j,k,1)

       enddo ! j
      enddo  ! k
      endif                          !16-12-15>


      end subroutine obc_scal_moveforward

!.....................................................................................

      subroutine obc_scal_wavecondition
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='obc_scal_wavecondition'
       subroutinedescription= &
       'Radiative open boundary conditions on T and S'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Passage au temps "after" le 26-01-23

! En situation de courant sortant, on applique un gradient nul, ce qui
! equivaut à faire fonctionner le schema d'advection comme un schema
! upwind.

! En situation de courant entrant:
! On applique: gradient model_ = gradient OGCM avec une contrainte qui
! empeche la valeur limite de s'eloigner de l'ogcm.
! La condition de gradient nul pour la deuxieme couronne peripherique, requise
! par la fonction "limiteur" du schéma d'advection, conduit le schema
! d'advection a fonctionné en mode upwind en cas de courant entrant.

! Ecart max autorisé par rapport à l'ogcm:                           !19/12/08
      x3=0. ! fois la variabilite de l'ogcm sur 2 echeances successives...
            ! si X3=0 la limite est l'ogcm.

! j=0 South open boundary (if any)
      if(obcstatus(jeq1)==1) then !-----------> !18-03-21

      do k=1,kmax
       do i=1,imax

       x0=0.5*sign(un, veldxdz_v(i,1     ,k,1))

!---- sal_t(0     ):

       x2=x3*abs(salobc_t(i,1  ,k,2)-salobc_t(i,1  ,k,0))

         sal_t(i,0     ,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        salobc_t(i,1     ,k,1)                                           & !18/12/02
       +max(-x2,min(x2,sal_t(i,2     ,k,2)-salobc_t(i,2     ,k,1)  ))  )  &

       +(0.5-x0)*sal_t(i,1   ,k,2)     ! condition courant sortant

!---- sal_t(-1    ):
         sal_t(i,-1    ,k,2)=sal_t(i,0     ,k,2)


!---- tem_t(0     ):

       x2=x3*abs(temobc_t(i,1  ,k,2)-temobc_t(i,1  ,k,0))

         tem_t(i,0     ,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        temobc_t(i,1     ,k,1)                                           & !18/12/02
       +max(-x2,min(x2,tem_t(i,2     ,k,2)-temobc_t(i,2     ,k,1)  ))  )  &

       +(0.5-x0)*tem_t(i,1   ,k,2)     ! condition courant sortant

!---- tem_t(-1    ):
         tem_t(i,-1    ,k,2)=tem_t(i,0     ,k,2)

       enddo ! i
      enddo  ! k

      endif                                      !----------->

! i=0 West open boundary (if any)
      if(obcstatus(ieq1)==1) then !-----------> !18-03-21

      do k=1,kmax
       do j=1,jmax

       x0=0.5*sign(un, veldydz_u(1     ,j,k,1))

!---- sal_t(0     ):

       x2=x3*abs(salobc_t(1  ,j,k,2)-salobc_t(1  ,j,k,0))

         sal_t(0     ,j,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        salobc_t(1     ,j,k,1)                                            &
       +max(-x2,min(x2,sal_t(2     ,j,k,2)-salobc_t(2     ,j,k,1)  )) )   &

       +(0.5-x0)*sal_t(1   ,j,k,2)     ! condition courant sortant

!---- sal_t(-1    ):
         sal_t(-1    ,j,k,2)=sal_t(0     ,j,k,2)

!---- tem_t(0     ):

       x2=x3*abs(temobc_t(1  ,j,k,2)-temobc_t(1  ,j,k,0))

         tem_t(0     ,j,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        temobc_t(1     ,j,k,1)                                            &
       +max(-x2,min(x2,tem_t(2     ,j,k,2)-temobc_t(2     ,j,k,1)  )) )   &

       +(0.5-x0)*tem_t(1   ,j,k,2)     ! condition courant sortant

!---- tem_t(-1    ):
         tem_t(-1    ,j,k,2)=tem_t(0     ,j,k,2)

       enddo ! j
      enddo  ! k

      endif                                      !----------->


! j=jmax+1 North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !-----------> !18-03-21

      do k=1,kmax
       do i=1,imax

       x0=0.5*sign(un,-veldxdz_v(i,jmax+1,k,1))

!---- sal_t(jmax+1):

       x2=x3*abs(salobc_t(i,jmax  ,k,2)-salobc_t(i,jmax  ,k,0))            !04/03/08

         sal_t(i,jmax+1,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        salobc_t(i,jmax  ,k,1)                                            &
       +max(-x2,min(x2,sal_t(i,jmax-1,k,2)-salobc_t(i,jmax-1,k,1) )) )    &

       +(0.5-x0)*sal_t(i,jmax,k,2)     ! condition courant sortant

!---- sal_t(jmax+2):
         sal_t(i,jmax+2,k,2)=sal_t(i,jmax+1,k,2)

!---- tem_t(jmax+1):

       x2=x3*abs(temobc_t(i,jmax  ,k,2)-temobc_t(i,jmax  ,k,0))            !04/03/08

         tem_t(i,jmax+1,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        temobc_t(i,jmax  ,k,1)                                            &
       +max(-x2,min(x2,tem_t(i,jmax-1,k,2)-temobc_t(i,jmax-1,k,1) )) )    &

       +(0.5-x0)*tem_t(i,jmax,k,2)     ! condition courant sortant

!---- tem_t(jmax+2):
         tem_t(i,jmax+2,k,2)=tem_t(i,jmax+1,k,2)


       enddo ! i
      enddo  ! k

      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> !18-03-21

      do k=1,kmax
       do j=1,jmax

       x0=0.5*sign(un,-veldydz_u(imax+1,j,k,1))

!---- sal_t(imax+1):

       x2=x3*abs(salobc_t(imax  ,j,k,2)-salobc_t(imax  ,j,k,0))

         sal_t(imax+1,j,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        salobc_t(imax  ,j,k,1)                                            &
       +max(-x2,min(x2,sal_t(imax-1,j,k,2)-salobc_t(imax-1,j,k,1)  )) )   &

       +(0.5-x0)*sal_t(imax,j,k,2)     ! condition courant sortant

!---- sal_t(imax+2):
         sal_t(imax+2,j,k,2)=sal_t(imax+1,j,k,2)


!---- tem_t(imax+1):

       x2=x3*abs(temobc_t(imax  ,j,k,2)-temobc_t(imax  ,j,k,0))

         tem_t(imax+1,j,k,2)=                                           &

        (0.5+x0)*(                   & ! condition courant entrant
        temobc_t(imax  ,j,k,1)                                            &
       +max(-x2,min(x2,tem_t(imax-1,j,k,2)-temobc_t(imax-1,j,k,1)  )) )   &

       +(0.5-x0)*tem_t(imax,j,k,2)     ! condition courant sortant

!---- tem_t(imax+2):
         tem_t(imax+2,j,k,2)=tem_t(imax+1,j,k,2)


       enddo ! j
      enddo  ! k

      endif                                      !----------->

      end subroutine obc_scal_wavecondition

!.....................................................................................

      subroutine obc_scal_mpi_now
      use module_principal
      use module_parallele
      implicit none
      integer loop_,id_tem_now_,id_sal_now_
#ifdef synopsis
       subroutinetitle='obc_scal_mpi_now'
       subroutinedescription= &
       'Lateral boundary conditions ensuring mpi continuity of the T ' &
       //'and S fields at time=now'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      if(subcycle_onoff==0) then !------->
       call get_type_echange('zb','tem_zb_n',tem_t,lbound(tem_t),ubound(tem_t),1,id_tem_now_)
       call get_type_echange('zb','sal_zb_n',sal_t,lbound(sal_t),ubound(sal_t),1,id_sal_now_)
      else                       !------->
       call get_type_echange('zc','tem_zc_n',tem_t,lbound(tem_t),ubound(tem_t),1,id_tem_now_)
       call get_type_echange('zc','sal_zc_n',sal_t,lbound(sal_t),ubound(sal_t),1,id_sal_now_)
      endif                      !------->

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(tem_t,id_tem_now_,mpi_neighbor_list(loop_))
         call echange_voisin(sal_t,id_sal_now_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_mpi_now

!.....................................................................................

      subroutine obc_scal_1dv
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='obc_scal_1dv'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! cas du model_ 1D (I1D=0)
       do k=1,kmax
        do j=-1,jmax+2                        !01-06-09
         do i=-1,imax+2
          tem_t(i,j,k,1)=tem_t(2,2,k,1)
          sal_t(i,j,k,1)=sal_t(2,2,k,1)
         enddo
        enddo
       enddo
      end subroutine obc_scal_1dv

!.....................................................................................

      subroutine obc_scal_mpi_after
      use module_principal
      use module_parallele
      implicit none
      integer loop_,id_tem_after_,id_sal_after_
#ifdef synopsis
       subroutinetitle='obc_scal_mpi_after'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      if(subcycle_onoff==0) then !------->
       call get_type_echange('zb','tem_zb_a',tem_t,lbound(tem_t),ubound(tem_t),2,id_tem_after_)
       call get_type_echange('zb','sal_zb_a',sal_t,lbound(sal_t),ubound(sal_t),2,id_sal_after_)
      else                       !------->
       call get_type_echange('zc','tem_zc_a',tem_t,lbound(tem_t),ubound(tem_t),2,id_tem_after_)
       call get_type_echange('zc','sal_zc_a',sal_t,lbound(sal_t),ubound(sal_t),2,id_sal_after_)
      endif                      !------->

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(tem_t,id_tem_after_,mpi_neighbor_list(loop_))
         call echange_voisin(sal_t,id_sal_after_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_mpi_after

!.....................................................................................

      subroutine obc_scal_checkmpi !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer loop_,id_tem_now_,id_sal_now_
#ifdef synopsis
       subroutinetitle='obc_scal_checkmpi'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Called by subroutine scalar_check_mpi_conservation
! z0 mpi exchange on tem_t and sal_t

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange('z0','tem_z0_n',tem_t,lbound(tem_t),ubound(tem_t),1,id_tem_now_)
       call get_type_echange('z0','sal_z0_n',sal_t,lbound(sal_t),ubound(sal_t),1,id_sal_now_)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(tem_t,id_tem_now_,mpi_neighbor_list(loop_))
         call echange_voisin(sal_t,id_sal_now_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_checkmpi

!.....................................................................................

      subroutine obc_scal_temobc(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,id_temobc_
#ifdef synopsis
       subroutinetitle='obc_scal_temobc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Called by subroutine module_ogcm
! z0 mpi exchange on temobc_t

#ifdef parallele
      write(texte30,'(a12,i0)')'temobc_t_z0_',t_     !15-07-14
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange('z0',trim(texte30)     &
                                  ,temobc_t         &
                           ,lbound(temobc_t)        &
                           ,ubound(temobc_t)        &
                           ,t_                      &
                           ,id_temobc_)

      ! Echanges
      do loop_=1, subcycle_exchange
       call echange_voisin(temobc_t,id_temobc_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_temobc

!.....................................................................................

      subroutine obc_scal_salobc(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,id_salobc_
#ifdef synopsis
       subroutinetitle='obc_scal_salobc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Called by subroutine module_ogcm
! z0 mpi exchange on salobc_t

#ifdef parallele
      write(texte30,'(a12,i0)')'salobc_t_z0_',t_     !15-07-14
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange('z0',trim(texte30)     &
                                  ,salobc_t         &
                           ,lbound(salobc_t)        &
                           ,ubound(salobc_t)        &
                           ,t_                      &
                           ,id_salobc_)

      ! Echanges
      do loop_=1, subcycle_exchange
       call echange_voisin(salobc_t,id_salobc_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_salobc

!.....................................................................................

      subroutine obc_scal_temref(txt_)
      use module_principal
      use module_parallele
      implicit none
      integer loop_
      character*2 txt_

#ifdef bidonref

! Called by subroutine module_ogcm
! z0 mpi exchange on temref_t

#ifdef parallele
      texte30='temref_t_'//txt_
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange(txt_,trim(texte30)     &
                                  ,temref_t         &
                           ,lbound(temref_t)        &
                           ,ubound(temref_t)        &
                           ,k0)

      ! Echanges
      do loop_=1, subcycle_exchange
       call echange_voisin(temref_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif
#endif

      end subroutine obc_scal_temref

!.....................................................................................

      subroutine obc_scal_salref(txt_)
      use module_principal
      use module_parallele
      implicit none
      integer loop_
      character*2 txt_

#ifdef bidonref

! Called by subroutine module_ogcm
! z0 mpi exchange on salref_t

#ifdef parallele
      texte30='salref_t_'//txt_
!!$! Nouvelle methode avec choix des voisins
       call get_type_echange(txt_,trim(texte30)     &
                                  ,salref_t         &
                           ,lbound(salref_t)        &
                           ,ubound(salref_t)        &
                           ,k0)

      ! Echanges
      do loop_=1, subcycle_exchange
       call echange_voisin(salref_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif
#endif

      end subroutine obc_scal_salref

!.....................................................................................
#ifdef bidon
      subroutine obc_scal_s_iflux_cumul(id_flx_,loop_)
      use module_principal ; use module_parallele 
      implicit none
      integer id_flx_,loop_ 

! Cette routine calcule le flux de sel et de temperature A travers la frontiere ouverte
! En l'absence d'autres sources que les frontieres (pas de flux A la surface) ce cumul
! equilibre la variation du contenu total de S et T tel qu'il est calculE par la subroutine
! my_outputs_sal_sum. Exemple de script gnuplot (les 2 courbes doivent se superpposer):
! plot "tmp/salmean" u 0:($2+$3) w l, "tmp/cumul_flux_ts" u 0:(-$3+459445309464.294) w l
 


      sum1=0.
      if(obcstatus(ieq1)==1) then    !pmxpmx>
        i=1
        do k=1,kmax ; do j=1,jmax
         sum1=sum1+anyv3d(i,j,k,id_flx_)*mask_j_u(j)*wetmask_t(i,j)
        enddo       ; enddo
      endif                          !pmxpmx>
      if(obcstatus(ieqimax)==1) then !pmxpmx>
        i=imax+1
        do k=1,kmax ; do j=1,jmax
         sum1=sum1-anyv3d(i,j,k,id_flx_)*mask_j_u(j)*wetmask_t(i-1,j)
        enddo       ; enddo
      endif                          !pmxpmx>
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      s_flux_cumul=s_flux_cumul+sum1glb*dti_lpsub*0.5

      end subroutine obc_scal_s_iflux_cumul
!...
      subroutine obc_scal_s_jflux_cumul(id_fly_,loop_)
      use module_principal ; use module_parallele 
      implicit none
      integer id_fly_,loop_ 

      sum1=0.
      if(obcstatus(jeq1)==1) then    !pmxpmx>
        j=1
        do k=1,kmax ; do i=1,imax
         sum1=sum1+anyv3d(i,j,k,id_fly_)*mask_i_v(i)*wetmask_t(i,j)
        enddo       ; enddo
      endif                          !pmxpmx>
      if(obcstatus(jeqjmax)==1) then !pmxpmx>
        j=jmax+1
        do k=1,kmax ; do i=1,imax
         sum1=sum1-anyv3d(i,j,k,id_fly_)*mask_i_v(i)*wetmask_t(i,j-1)
        enddo       ; enddo
      endif                          !pmxpmx>
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      s_flux_cumul=s_flux_cumul+sum1glb*dti_lpsub*0.5

      if(loop_==loopmaxts) then !styx>
      if(par%rank==0) then !cid>
        open(unit=3,file='tmp/cumul_flux_ts',position='append')
         write(3,*)real(elapsedtime_now/86400.),real(t_flux_cumul),real(s_flux_cumul)
        close(3)
      endif                !cid>
      endif                     !styx>
     

      end subroutine obc_scal_s_jflux_cumul
!...
      subroutine obc_scal_t_iflux_cumul(id_flx_,loop_)
      use module_principal ; use module_parallele 
      implicit none
      integer id_flx_,loop_ 

      sum1=0.
      if(obcstatus(ieq1)==1) then    !pmxpmx>
        i=1
        do k=1,kmax ; do j=1,jmax
         sum1=sum1+anyv3d(i,j,k,id_flx_)*mask_j_u(j)*wetmask_t(i,j)
        enddo       ; enddo
      endif                          !pmxpmx>
      if(obcstatus(ieqimax)==1) then !pmxpmx>
        i=imax+1
        do k=1,kmax ; do j=1,jmax
         sum1=sum1-anyv3d(i,j,k,id_flx_)*mask_j_u(j)*wetmask_t(i-1,j)
        enddo       ; enddo
      endif                          !pmxpmx>
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      t_flux_cumul=t_flux_cumul+sum1glb*dti_lpsub*0.5

      end subroutine obc_scal_t_iflux_cumul
!...
      subroutine obc_scal_t_jflux_cumul(id_fly_,loop_)
      use module_principal ; use module_parallele 
      implicit none
      integer id_fly_,loop_

      sum1=0.
      if(obcstatus(jeq1)==1) then    !pmxpmx>
        j=1
        do k=1,kmax ; do i=1,imax
         sum1=sum1+anyv3d(i,j,k,id_fly_)*mask_i_v(i)*wetmask_t(i,j)
        enddo       ; enddo
      endif                          !pmxpmx>
      if(obcstatus(jeqjmax)==1) then !pmxpmx>
        j=jmax+1
        do k=1,kmax ; do i=1,imax
         sum1=sum1-anyv3d(i,j,k,id_fly_)*mask_i_v(i)*wetmask_t(i,j-1)
        enddo       ; enddo
      endif                          !pmxpmx>
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      t_flux_cumul=t_flux_cumul+sum1glb*dti_lpsub*0.5

      end subroutine obc_scal_t_jflux_cumul
#endif
!.....................................................................................
      subroutine obc_scal_bottom_rhp(txt_)  !21-05-18
      use module_principal
      use module_parallele
      implicit none
      character txt_*2
      integer loop_,id_botrhp_

! Cette routine sert A la continuite mpi du champ rhp pour k=1 utilisE dans le diag
! d'instabilite de pente

! West open boundary (if any) !26-11-14
      if(obcstatus(ieq1)==1) then !-----------> !18-03-21
       do j=1,jmax
         rhp_t(-1:0,j,1)=rhp_t(1,j,1)
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> !18-03-21
       do j=1,jmax
         rhp_t(imax+1:imax+2,j,1)=rhp_t(imax,j,1)
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !-----------> !18-03-21
       do i=1,imax
         rhp_t(i,-1:0,1)=rhp_t(i,1,1)      
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !-----------> !18-03-21
       do i=1,imax
         rhp_t(i,jmax+1:jmax+2,1)=rhp_t(i,jmax,1)
       enddo
      endif                                      !----------->


#ifdef parallele
      write(texte30,'(a7,a2)')'botrhp_',txt_   
       call get_type_echange(txt_,trim(texte30)   & !15-07-14
                                     ,rhp_t       &
                              ,lbound(rhp_t)      &
                              ,ubound(rhp_t)      &
                              ,1                  & ! Seul le niveau du fond est echange
                              ,id_botrhp_)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(rhp_t,id_botrhp_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_bottom_rhp

!.................................................................

      subroutine obc_scal_mpi_upwindriver(txt_) !02-06-18
      use module_principal
      use module_parallele
      implicit none
      character*2 txt_
      integer loop_

      write(texte30,'(a9,a2)')'upwindriver_t_',txt_ !16-03-15

#ifdef parallele
      call get_type_echange(txt_,trim(texte30),upwindriver_t        &
                                       ,lbound(upwindriver_t)       &
                                       ,ubound(upwindriver_t)       &
                                       ,k0)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(upwindriver_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_mpi_upwindriver

!.................................................................

      subroutine obc_scal_rhp !01-10-21
      use module_principal
      use module_parallele
      implicit none
      integer loop_

! Cette routine sert A la continuite mpi du champ rhp pour k=1 utilisE dans le limiteur de l'advection horizontale

! West open boundary (if any) !26-11-14
      if(obcstatus(ieq1)==1) then !-----------> !18-03-21
       do k=1,kmax ; do j=1,jmax
         rhp_t(-1:0,j,k)=rhp_t(1,j,k)
       enddo ; enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> !18-03-21
       do k=1,kmax ; do j=1,jmax
         rhp_t(imax+1:imax+2,j,k)=rhp_t(imax,j,k)
       enddo ; enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !-----------> !18-03-21
       do k=1,kmax ; do i=1,imax
         rhp_t(i,-1:0,k)=rhp_t(i,1,k)      
       enddo ; enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !-----------> !18-03-21
       do k=1,kmax ; do i=1,imax
         rhp_t(i,jmax+1:jmax+2,k)=rhp_t(i,jmax,k)
       enddo ; enddo
      endif                                      !----------->


#ifdef parallele
       call get_type_echange('z2','rhp_t_z2'      & !26-01-23
                                     ,rhp_t       &
                              ,lbound(rhp_t)      &
                              ,ubound(rhp_t)      &
                              ,k0)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(rhp_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_scal_rhp

!.................................................................
