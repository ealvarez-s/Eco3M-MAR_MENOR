      subroutine obc_int
!______________________________________________________________________
! SYMPHONIE ocean model
! release 290 - last update: 30-10-20
!______________________________________________________________________
!    _________                    .__                  .__             !
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int'
       subroutinedescription='Radiative OBC on 3D velocity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         20/07/01: passage à la coordonnée sigma généralisée
!                  + amelioration sur SUM1 (n'aurait pas marcher pour SUM1<0)
!         25/11/01: prise en compte des petits fleuves (introduits en surface)
!         26/08/02: VHZOBC remplacé par VELOBC
!         24/11/02: gradient nul remplacé par une equation d'onde
!         28/11/02: gradient nul remplacé par une equation d'onde: suite debug
!         02/12/02: correction bug. Erreur de sens de propagation sur frontieres
!                   i=1 et j=1
!         15/12/02: les tableaux de vitesse des ondes sont maintenant calcules
!                   dans cwave_int.F
!         10/07/03: amenagement pour grille hybride
!         25/05/04: reecriture de l'equation d'onde. A priori c'est plus propre
!                   et on s'apercoit que pour un nombre d'onde = 1 on retombe
!                   bien sur une condition de gradient nul
!         25/06/04: un bon toiletage general....
!         15/09/05: des nouveaux cas ont été ajouté dans l'objectif de faire
!                   des test de comparaisons à l'avenir. Dans l'immediat c'est
!                   toujours le cas n°2 qui tient la preference.
!         09/01/07: multiplier les coins par mask et mettre ce calcul dans une
!                   subroutine
!         20/04/07: passage à la coordonnée curviligne
!         29/08/07: un seul choix, le schéma classique:
!                   v'(i,t+1)=( (1-a)*v'(i,t-1)+2*a*v'(i+-1,t) )/(1+a)
!         06/11/08: Intervention J. Floor: correction d'un bug sur la partie
!                   forcage de la condition radiative
!         05-06-09  Parallelisation
! 2009.3  01-10-09  utilisation des nouveaux facteurs d'echelle verticale
!         05-10-09  ajout d'un "ifdef parallele"
! 2010.8  10-03-10  beforebefore remplacé par before2
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S26     27-01-13  echange advection qdm o4
!         11-11-13  modif pour subcycling
!         08-05-14  call obc_int_mpi
!         02-07-14  subroutine obc_int_mpi4d
!         11-07-14  dt_obc devient dtobc(:)
!         21-08-14  ajout subroutine obc_int_cn_o4
!         16-03-15  ajout subroutine obc_int_anyvar2d
!         21-03-15  'za' dans la check list de obc_int_anyvar2d
!         15-07-15  modification des arguments passes dans obc_int_mpi
!         09-08-16  cas echange 'r' ajoutE
!         13-12-18  obc coins couches funionnees
! v252    22-04-19  ajout reseau canaux
! v287    14-08-20  subroutine obc_int_anyv3d(arg4,txt_) changee pour 
!                              obc_mpi_anyv3d(option_,arg4_,txt_) le
!                   nouveau parametre "option_" permettant d'inclure ou
!                   non des actions sur les canaux
! v290    30-10-20  subroutine obc_int_mpi_anyvar3d(txt_) !30-10-20
!...............................................................................

!...............................................................................
! Description:
!
! 1.      Calcule les conditions aux limites pour
! la partie cisaillee du courant. La composante moyenne
! est elle calculee par une condition de radiation calculee
! par obc_ext.F
!
! 2.     D'autre part la condition porte sur la perturbation.
! Soit u=u'+ug, où u est le courant, ug le courant geostrophique et u' la
! perturbation. C'est sur u' qu'on applique la condition, c.a.d.
! du'/dt + c*du'/dx = 0, soit du/dt+c*du/dx=dug/dt+c*dug/dx
! Le forçage au frontiere endosse le rôle du
! courant à variabilité lente ug (c'est pas choquant si ce dernier est filé
! par des sorties Mercator, ca le deviendra avec du nesting haute frequence).
!...............................................................................

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       i1=imax ; i3=i1-1 ; i2=imax+1 ; i4=i2-1 ; loop1=2 ; call obc_int_eastwest
      endif                                      !----------->

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       i1=1 ; i3=i1+1 ; i2=1 ; i4=i2+1 ; loop1=1 ; call obc_int_eastwest
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       j1=jmax ; j3=j1-1 ; j2=jmax+1 ; j4=j2-1 ; loop1=2 ; call obc_int_northsouth
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       j1=1 ; j3=j1+1 ; j2=1 ; j4=j2+1 ; loop1=1 ; call obc_int_northsouth
      endif                                      !----------->

! South-West open boundary (if any)
      if(obcstatus(ieq1_jeq1)==1)call obc_int_sudouest

! South-East open boundary (if any)
      if(obcstatus(ieqimax_jeq1)==1)call obc_int_sudest

! North-East open boundary (if any)
      if(obcstatus(ieqimax_jeqjmax)==1)call obc_int_nordest

! North-West open boundary (if any)
      if(obcstatus(ieq1_jeqjmax)==1)call obc_int_nordouest

! Lignes colonnes supplementaires pour advection Ordre 4 (gradient nul)
      call obc_int_o4

! cas du model_ 1D (I1D=0)
      if(flag3d==0) then   !1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D>
      do k=1,kmax ; do i=1,imax+1 ; do j=1,jmax+1
       vel_u(i,j,k,2)=vel_u(2,2,k,2)
       vel_v(i,j,k,2)=vel_v(2,2,k,2)
      enddo ; enddo ; enddo
! cas du model_ 1D: fin.
      endif             !1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D1D>

      end subroutine obc_int

!...............................................................................

      subroutine obc_int_mpi(t_,obctype_) !13-07-15
      use module_principal ; use module_parallele
      implicit none
      integer loop_,idi_u_u1_,idi_u_u2_,idi_v_v1_,idi_v_v2_,t_ &
             ,obctype_
             
#ifdef synopsis
       subroutinetitle='obc_int_mpi'
       subroutinedescription= &
       'mpi continuity of the 3D velocity field at time=after'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      if(obctype_==1.or.obctype_==12) then !u1v1u1v1>

       write(texte30,'(a,i0)')'vel_u1_',t_ 
            call get_type_echange('u1',trim(texte30),vel_u        &
                                             ,lbound(vel_u)       &
                                             ,ubound(vel_u),t_    &
                                                    ,idi_u_u1_)

       write(texte30,'(a,i0)')'vel_v1_',t_ 
            call get_type_echange('v1',trim(texte30),vel_v        &
                                             ,lbound(vel_v)       &
                                             ,ubound(vel_v),t_    &
                                                    ,idi_v_v1_)

      endif                                !u1v1u1v1>


      if(obctype_==2.or.obctype_==12) then !u2v2u2v2>

       write(texte30,'(a,i0)')'vel_u2_',t_ 
            call get_type_echange('u2',trim(texte30),vel_u        &
                                             ,lbound(vel_u)       &
                                             ,ubound(vel_u),t_    &
                                                    ,idi_u_u2_)

       write(texte30,'(a,i0)')'vel_v2_',t_ 
            call get_type_echange('v2',trim(texte30),vel_v        &
                                             ,lbound(vel_v)       &
                                             ,ubound(vel_v),t_    &
                                                    ,idi_v_v2_)

      endif                                !u2v2u2v2>

      ! Echanges
      do loop_=1, subcycle_exchange

        if(obctype_==1.or.obctype_==12) then !u1v1u1v1>
         call echange_voisin(vel_u,idi_u_u1_,mpi_neighbor_list(loop_))
         call echange_voisin(vel_v,idi_v_v1_,mpi_neighbor_list(loop_))
        endif                                !u1v1u1v1>

        if(obctype_==2.or.obctype_==12) then !u2v2u2v2>
         call echange_voisin(vel_u,idi_u_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(vel_v,idi_v_v2_,mpi_neighbor_list(loop_))
        endif                                !u2v2u2v2>

      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_int_mpi

!..............................................................

      subroutine obc_int_o4
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_o4'
       subroutinedescription= &
       '"+2 range" lateral boundary conditions on 3D velocity field' &
       //' required by high-order advection schemes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Lignes colonnes supplementaires pour advection O4

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> !11-11-13
       do k=1,kmax
       do j=1,jmax
        vel_u(imax+2,j,k,2)=vel_u(imax+1,j,k,2)
       enddo
       enddo
       do k=1,kmax
       do j=0,jmax+2
        vel_v(imax+1,j,k,2)=vel_v(imax,j,k,2)
       enddo
       enddo
      endif                                      !----------->

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       do k=1,kmax
       do j=1,jmax
        vel_u(0     ,j,k,2)=vel_u(1     ,j,k,2)
       enddo
       enddo
       do k=1,kmax
       do j=0,jmax+2
        vel_v(0     ,j,k,2)=vel_v(1   ,j,k,2)
       enddo
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do k=1,kmax
       do i=0,imax+2
        vel_u(i,jmax+1,k,2)=vel_u(i,jmax,k,2)
       enddo
       enddo
       do k=1,kmax
       do i=1,imax
        vel_v(i,jmax+2,k,2)=vel_v(i,jmax+1,k,2)
       enddo
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do k=1,kmax
       do i=0,imax+2
        vel_u(i,0     ,k,2)=vel_u(i,1   ,k,2)
       enddo
       enddo
       do k=1,kmax
       do i=1,imax
        vel_v(i,0     ,k,2)=vel_v(i,1     ,k,2)
       enddo
       enddo
      endif                                      !----------->

      end subroutine obc_int_o4
!............................................................................
      subroutine obc_int_northsouth
      use module_principal
      use module_parallele !#MPI
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_northsouth'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      const1=2.*dti_fw/dtobc(vel_id)
      const3=0.25*assel0

      do 880 k=1,kmax

! Composante _X:
      do 881 i=2,imax

       x0=dti_fw*cwi_int_u(i,loop1)/dy_u(i,j1)

             vel_u(i,j1,k,2)=(                     &
             vel_u(i,j1,k,0)*(1.-x0)               &

       +( velobc_u(i,j1,k,2)                           &
         -velobc_u(i,j1,k,0))*const1                   &
       +( velobc_u(i,j1,k,1)                           &
            +vel_u(i,j3,k,1)                           &
         -velobc_u(i,j3,k,1))*2.*x0                    & !06/11/08

!$ Laplacian time filter:
        +const3*(  (dz_u(i,j1,k,after)  +dz_u(i,j1,k,now))            &
                 *(vel_u(i,j1,k,now)   -vel_u(i,j1,k,before))         &
                  -(dz_u(i,j1,k,now)    +dz_u(i,j1,k,before))         &
                 *(vel_u(i,j1,k,before)-vel_u(i,j1,k,before2     )))  &
                   /dz_u(i,j1,k,after)                               &

                             )/(1.+x0)*mask_u(i,j1,k)

  881 continue

! Composante _Y:
      do 883 i=2,imax-1

      x0=dti_fw*cwi_int_v(i,loop1)/dy_v(i,j2)

             vel_v(i,j2,k,2)=(                       &
             vel_v(i,j2,k,0)*(1.-x0)                 &

       +( velobc_v(i,j2,k,2)                         &
         -velobc_v(i,j2,k,0))*const1                 &
       +( velobc_v(i,j2,k,1)                         &
            +vel_v(i,j4,k,1)                         &
         -velobc_v(i,j4,k,1))*2.*x0                  & !06/11/08

!$ Laplacian time filter:
        +const3*(  (dz_v(i,j2,k,after)  +dz_v(i,j2,k,now))            &
                 *(vel_v(i,j2,k,now)   -vel_v(i,j2,k,before))         &
                  -(dz_v(i,j2,k,now)    +dz_v(i,j2,k,before))         &
                 *(vel_v(i,j2,k,before)-vel_v(i,j2,k,before2     )))  &
                   /dz_v(i,j2,k,after)                       &

                             )/(1.+x0)*mask_v(i,j2,k)

  883 continue

  880 continue

      end subroutine obc_int_northsouth

!............................................................................

      subroutine obc_int_eastwest
      use module_principal
      use module_parallele !#MPI
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_eastwest'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      const1=2.*dti_fw/dtobc(vel_id) !11-07-14
      const3=0.25*assel0

      do 880 k=1,kmax

! Composante _X:
      do 882 j=2,jmax-1

       x0=dti_fw*cwj_int_u(j,loop1)/dx_u(i2,j)

             vel_u(i2,j,k,2)=(                                      &
             vel_u(i2,j,k,0)*(1.-x0)                                &

       +( velobc_u(i2,j,k,2)                                        &
         -velobc_u(i2,j,k,0))*const1                                &
       +( velobc_u(i2,j,k,1)                                        &
            +vel_u(i4,j,k,1)                                        &
         -velobc_u(i4,j,k,1))*2.*x0                                 & !06/11/08

!$ Laplacian time filter:
        +const3*(  (dz_u(i2,j,k,after)  +dz_u(i2,j,k,now))            &
                 *(vel_u(i2,j,k,now)   -vel_u(i2,j,k,before))         &
                  -(dz_u(i2,j,k,now)    +dz_u(i2,j,k,before))         &
                 *(vel_u(i2,j,k,before)-vel_u(i2,j,k,before2     )))  &
                   /dz_u(i2,j,k,after)                                &

                             )/(1.+x0)*mask_u(i2,j,k)

  882 continue

! Composante _Y:
      do 884 j=2,jmax

      x0=dti_fw*cwj_int_v(j,loop1)/dx_v(i1,j)

             vel_v(i1,j,k,2)=(                      &
             vel_v(i1,j,k,0)*(1.-x0)                &

       +( velobc_v(i1,j,k,2)                        &
         -velobc_v(i1,j,k,0))*const1                &
       +( velobc_v(i1,j,k,1)                        &
            +vel_v(i3,j,k,1)                        &
         -velobc_v(i3,j,k,1))*2.*x0                 & !06/11/08

!$ Laplacian time filter:
        +const3*(  (dz_v(i1,j,k,after)  +dz_v(i1,j,k,now))            &
                 *(vel_v(i1,j,k,now)   -vel_v(i1,j,k,before))         &
                  -(dz_v(i1,j,k,now)    +dz_v(i1,j,k,before))         &
                 *(vel_v(i1,j,k,before)-vel_v(i1,j,k,before2     )))  &
                   /dz_v(i1,j,k,after)                       &

                            )/(1.+x0)*mask_v(i1,j,k)

  884 continue

  880 continue

      end subroutine obc_int_eastwest

!............................................................................

      subroutine obc_int_sudouest
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_sudouest'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Les coins de la composante X:
! COIN n°1:
       i1=1 ; j1=1 ; i2=2 ; j2=2
!      do k=kmin_u(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_u(i1,j1,k,2)=                                    &
          (vel_u(i1,j2,k,2)+vel_u(i2,j1,k,2))/2.                &
         *mask_u(i1,j1,k)
       enddo

! Les coins de la composante Y:
! COIN n°1:
       i1=1 ; j1=1 ; i2=2 ; j2=2
!      do k=kmin_v(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_v(i1,j1,k,2)=                                    &
          (vel_v(i1,j2,k,2)+vel_v(i2,j1,k,2))/2.                &
         *mask_v(i1,j1,k)
       enddo

       end subroutine obc_int_sudouest

!..............................................................

      subroutine obc_int_sudest
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_sudest'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Les coins de la composante X:
! COIN n°2:
       i1=imax+1 ; j1=1 ; i2=imax ; j2=2
!      do k=kmin_u(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_u(i1,j1,k,2)=                                    &
          (vel_u(i1,j2,k,2)+vel_u(i2,j1,k,2))/2.                &
         *mask_u(i1,j1,k)
       enddo
! COIN n°2:
       i1=imax ; j1=1 ; i2=imax-1 ; j2=2
!      do k=kmin_v(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_v(i1,j1,k,2)=                                    &
          (vel_v(i1,j2,k,2)+vel_v(i2,j1,k,2))/2.                &
         *mask_v(i1,j1,k)
       enddo

      end subroutine obc_int_sudest

!..............................................................

      subroutine obc_int_nordest
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_nordest'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! COIN n°3:
       i1=imax+1 ; j1=jmax ; i2=imax ; j2=jmax-1
!      do k=kmin_u(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_u(i1,j1,k,2)=                                    &
          (vel_u(i1,j2,k,2)+vel_u(i2,j1,k,2))/2.                &
         *mask_u(i1,j1,k)
       enddo
! COIN n°3:
       i1=imax ; j1=jmax+1 ; i2=imax-1 ; j2=jmax
!      do k=kmin_v(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_v(i1,j1,k,2)=                                    &
          (vel_v(i1,j2,k,2)+vel_v(i2,j1,k,2))/2.                &
         *mask_v(i1,j1,k)
       enddo

       end subroutine obc_int_nordest

!..............................................................

      subroutine obc_int_nordouest
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='obc_int_nordouest'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! COIN n°4:
       i1=1 ; j1=jmax ; i2=2 ; j2=jmax-1
!      do k=kmin_u(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_u(i1,j1,k,2)=                                    &
          (vel_u(i1,j2,k,2)+vel_u(i2,j1,k,2))/2.                &
         *mask_u(i1,j1,k)
       enddo

! COIN n°4:
       i1=1 ; j1=jmax+1 ; i2=2 ; j2=jmax
!      do k=kmin_v(i1,j1),kmax
       do k=1,kmax !13-12-18
           vel_v(i1,j1,k,2)=                                    &
          (vel_v(i1,j2,k,2)+vel_v(i2,j1,k,2))/2.                &
         *mask_v(i1,j1,k)
       enddo

      end subroutine obc_int_nordouest

!...............................................................................

      subroutine obc_int_mpi4d !02-07-14
      use module_principal ; use module_parallele
      implicit none
      integer loop_,idi_u_u1_,idi_u_u2_,idi_v_v1_,idi_v_v2_
#ifdef synopsis
       subroutinetitle='obc_int_mpi4d'
       subroutinedescription= &
       'Lateral boundary condition ensuring mpi continuity of the ' &
       //'velocity 3D field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      call get_type_echange('u1','vel_u1_4d',vel_u                    &
                            ,lbound(vel_u),ubound(vel_u),idi_u_u1_)
      call get_type_echange('u2','vel_u2_4d',vel_u                    &
                            ,lbound(vel_u),ubound(vel_u),idi_u_u2_)

      call get_type_echange('v1','vel_v1_4d',vel_v                    &
                            ,lbound(vel_v),ubound(vel_v),idi_v_v1_)
      call get_type_echange('v2','vel_v2_4d',vel_v                    &
                            ,lbound(vel_v),ubound(vel_v),idi_v_v2_)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(vel_u,idi_u_u1_,mpi_neighbor_list(loop_))
         call echange_voisin(vel_u,idi_u_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(vel_v,idi_v_v1_,mpi_neighbor_list(loop_))
         call echange_voisin(vel_v,idi_v_v2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_int_mpi4d

!..............................................................

      subroutine obc_mpi_anyv3d(option_,arg4_,txt_)
      use module_principal ; use module_parallele ; use module_webcanals
      implicit none
      integer loop_,arg4_,option_
      character(len=*) txt_

! Note sur parametre "option_":
!    1- Tous les appels A la presente routine n'entrainent pas forcement un
!       traitement des canaux. Il faut pour cela que option_=1
!    2- La subroutine a ete renommee suite A l'ajout de ce parametre
!       afin de bien reconsiderer chaque appel vis A vis de cette
!       contingence

      if(option_==1.and. &       !14-08-20
         nbcanal>0) then !m°v°m> !22-04-19
! Comme les canaux sont along i, l'envoi de vitesse de gt1 A gt2
! ne concerne que u, d'oU le test sur u1 seulement (et pas aussi u2)
        if(txt_=='u1')call webcanals_gt1_to_gt2_anyvel(arg4_)
! Le test sur zb cible les traceurs tem, sal, bio
! Le test sur z0 cible la turbulence
        if(txt_=='zb'.or.txt_=='z0') then !>>>
                      call webcanals_gt1_to_gt2_anytrc(arg4_)
                      call webcanals_gt2_to_gt1_anytrc(arg4_)
        endif                             !>>>
      endif              !m°v°m>

      if(len(txt_)==1)write(texte30,'(a7,a1,a1,i0)')'anyv3d_',txt_,'_',arg4_
      if(len(txt_)==2)write(texte30,'(a7,a2,a1,i0)')'anyv3d_',txt_,'_',arg4_
      if(len(txt_)>2)stop 'obc_int_anyv3d len(txt_)>2' !09-08-16

!     if(txt_/='u1'.and.txt_/='v1'.and.txt_/='z1') &
!     stop 'Unrecognised argument in obc_int_anyv3d'

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      call get_type_echange(txt_,trim(texte30),anyv3d        &
                                       ,lbound(anyv3d)       &
                                       ,ubound(anyv3d)       &
                                       ,arg4_                &
                                       ,k0)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(anyv3d,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_mpi_anyv3d

!..............................................................

      subroutine obc_int_cn_o4 !21-08-14
      use module_principal ; use module_parallele
      implicit none
      integer idi_u_u2_,idi_v_v2_,loop_
#ifdef synopsis
       subroutinetitle='obc_int_cn_o4'
       subroutinedescription=                                          &
          '"+2 range" lateral boundary conditions on 3D dimensionless' &
       //' velocity veldtodx_u veldtody_o required by "upg" advection' &
       //' scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      stop 'je ne dois plus passer par obc_int_cn_o4'
#ifdef bidon

! Lignes colonnes supplementaires pour advection O4

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> !11-11-13
       do k=1,kmax
       do j=1,jmax
        veldtodx_u(imax+2,j,k,1)=veldtodx_u(imax+1,j,k,1)
       enddo
       enddo
       do k=1,kmax
       do j=0,jmax+2
        veldtody_v(imax+1,j,k,1)=veldtody_v(imax,j,k,1)
       enddo
       enddo
      endif                                      !----------->

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       do k=1,kmax
       do j=1,jmax
        veldtodx_u(0     ,j,k,1)=veldtodx_u(1     ,j,k,1)
       enddo
       enddo
       do k=1,kmax
       do j=0,jmax+2
        veldtody_v(0     ,j,k,1)=veldtody_v(1   ,j,k,1)
       enddo
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do k=1,kmax
       do i=0,imax+2
        veldtodx_u(i,jmax+1,k,1)=veldtodx_u(i,jmax,k,1)
       enddo
       enddo
       do k=1,kmax
       do i=1,imax
        veldtody_v(i,jmax+2,k,1)=veldtody_v(i,jmax+1,k,1)
       enddo
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do k=1,kmax
       do i=0,imax+2
        veldtodx_u(i,0     ,k,1)=veldtodx_u(i,1   ,k,1)
       enddo
       enddo
       do k=1,kmax
       do i=1,imax
        veldtody_v(i,0     ,k,1)=veldtody_v(i,1     ,k,1)
       enddo
       enddo
      endif                                      !----------->

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      call get_type_echange('u2'    &
                  ,'veldtodx_u2_t1' &
                   ,veldtodx_u      &
            ,lbound(veldtodx_u)     &
            ,ubound(veldtodx_u),1,idi_u_u2_)

      call get_type_echange('v2'    &
                  ,'veldtody_v2_t1' &
                   ,veldtody_v      &
            ,lbound(veldtody_v)     &
            ,ubound(veldtody_v),1,idi_v_v2_)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(veldtodx_u,idi_u_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(veldtody_v,idi_v_v2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

#endif
      end subroutine obc_int_cn_o4

!............................................................................

      subroutine obc_int_anyvar2d(txt_)
      use module_principal
      use module_parallele
      implicit none
      character*2 txt_
      integer loop_
#ifdef synopsis
       subroutinetitle='obc_int_anyvar2d'
       subroutinedescription= &
       'Lateral boundary conditions ensuring mpi continuity'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      write(texte30,'(a9,a2)')'anyvar2d_',txt_ !16-03-15

      if(txt_/='u1'.and.txt_/='v1'.and.txt_/='za') & !21-03-15
      stop 'Unrecognised argument in obc_int_anyvar2d'

#ifdef parallele
      call get_type_echange(txt_,trim(texte30),anyvar2d        &
                                       ,lbound(anyvar2d)       &
                                       ,ubound(anyvar2d)       &
                                       ,k0)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(anyvar2d,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_int_anyvar2d

!------------------------------------------------------------------------------

      subroutine obc_int_mpi_anyvar3d(txt_) !30-10-20
      use module_principal
      use module_parallele
      implicit none
      character*2 txt_
      integer loop_

      write(texte30,'(a9,a2)')'anyvar3d_',txt_ !16-03-15

#ifdef parallele
      call get_type_echange(txt_,trim(texte30),anyvar3d        &
                                       ,lbound(anyvar3d)       &
                                       ,ubound(anyvar3d)       &
                                       ,k0)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(anyvar3d,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_int_mpi_anyvar3d
