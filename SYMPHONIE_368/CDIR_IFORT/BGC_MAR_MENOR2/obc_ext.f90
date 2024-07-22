










      subroutine obc_ext(txt_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release S26 - last update: 25-05-17
!______________________________________________________________________
!    _________                    .__                  .__     m[�v�]m !
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!......................................................................
      use module_principal
      implicit none
      character txt_*5

!$ Compute the (physical) Open Boundary Conditions of the external mode momentum equations
!$ according to:  Marsaleix P., Auclair F., Estournel C., 2006, Considerations on Open
!$ Boundary Conditions for Regional and Coastal Ocean Models. Journal of Atmospheric and Oceanic
!$ Technology, 23,1604-1613, http://dx.doi.org/10.1175/JTECH1930.1

!...............................................................................
! Version Date      Description des modifications
!         25/11/01: introduction du cas riverdir=0 (fleuve introduit � la
!                   surface et non plus lateralement)
!         26/08/02: VBAROBC remplac� par VMEAOBC
!         11/06/04: procedure des coins revue. On calcule d'abord une ebauche
!                   des coins au moment de l'etape des composantes tgte � la
!                   frontiere. Ensuite un coef de lagrange est utilis� pour
!                   corriger le courant des coins afin d'etre compatible avec
!                   la sse. La perturbation du courant est proportionnelle au
!                   courant lui-m�me.
!         24/06/04: projet (les lignes sont en commentaire) de faire
!                   une correction du transport (correction bilan volume)
!                   qui soit proportionnelle � l'erreur sur le courant
!                   autrement dit l'ecart du transport � H*VMEAOBC.
!         12/07/04: Si relax_es < 0 pas de correction de la moyenne
!                   de l'elevation de la surface
!         03/08/05: Cas "Flather forcage libre"
!         10/11/05: debug du point precedent + extension de la methode
!                   aux composantes du transport tangeantes aux OBC
!         12/01/06: les variables termes contenant obc sont regroup�s
!                   et calcul�s une fois par iteration barocline
!                   La sse basse frequence voit sa moyenne control�e par
!                   ztaobc
!         20/01/06: amenagements sur conservation du toit
!         09/03/06: condition vitesses tangentes remis � la version originale
!                   car bug sinon. Appel � updateforcing rationnalise.
!         31/05/06: Correction d'un bug dans la condition au limite de T et S
!                   du au fait que les coins de la solution
!                   filtree (basse frequence) n'avait pas ete codes...
!         05/04/07: La condition aux limites integre maintenant la condition
!                   adaptative. Une version generalisee est employee, c'est �
!                   dire que la m�me condition de Flather est utilisee, il n'y
!                   a que la fonction de construire la condition de reference
!                   qui change.
!                   D'autre part, dans la cas standard, on a prevu une transition
!                   lineaire du forcage pour eviter les discontinuit�s li�es
!                   au time spliting, discontinuit�s qui se propageaient dans le
!                   domaine sous forme d'ondes de gravit�s.
!          17/04/07: Passage � coordonnees curvilignes (voir ajout dx_y et cie...)
!          13/12/07: Le point de variablilit� lente pour la vitesse est deplac�
!                    au point voisin.
!          12-05-09  le suivi de la moyenne de la ssh est introduit dans la filiere
!                    nesting
!          16-05-09  - Ajout du suivi de la moyenne de la ssh de la maree
!                    - Amortisseur de la correction moyenne revisit�
!          04-06-09  Parallelisation
! 2009.3   02-10-09  dte_lp remplace dte
!          05-10-09  ajout d'un "ifdef parallele"
!          09-10-09  parallelisation de la conservation du toit
! 2010.7   22-02-10  velbar devient la variable d'etat
!          01-03-10  correction du bug introduit avec le point precedent
! 2010.8   22-03-10  hssh_w remplace ssh_ext_w+h_w
!          24-03-10  La condition riviere sur fluxbar et la CL ouverte sur velbar
!                    sont disssossi�es car pas appliquees aux memes moments
! 2010.10  13-06-10  suppression ssh_ext_w
! 2010.13  01-11-10  river_inout remplace river_dom
! 2010.20  19-04-11  k2dite,k2dfin renomm�s iteration2d,iteration2d_max_now
! 2010.25  28-03-12  echange compatible avec pgf Stelios
! S26.1    27-01-13  echanges advection qdm o4
!          01-05-13  passage au time stepping fb
!          01-07-13  subcycle
!          12-05-14  nouveaux echanges
!          02-07-14  ajout subroutine obc_ext_mpi3d
!          15-07-14  ajout subroutine obc_ext_velbar_mpi
!          13-11-14  reorganisation pour plus de clarte
!          22-01-15  C.L. ssh dans les fleuves
!          08-02-15  ajout Flather sur courant
!          02-04-15  ajout echange r1 sur vlxbar_f et cie
!          18-08-15  ajout d'une conditionas aux limites pour les flux de qdm
!                    prise en compte de velbar_stokes et de omega_w dans
!                    l'inversion de l'equation de continuite
!          20-08-15  appel a call obc_ext_o4 inutile si appel a obc_ext_qdmflux
!          13-10-15  ajout obc_ext_xy_t
!          20-03-17  Ajout de subroutines pour controler les flux de volume aux detroits
!                    de Gibraltar et de Sicile
!          21-03-17  suite point precedent
!          25-05-17  Appel A obc_ext_o4
!........ .......................................................................


!$ Compute the (physical) Open Boundary Conditions of the external mode momentum equations
!$ according to:  Marsaleix P., Auclair F., Estournel C., 2006, Considerations on Open
!$ Boundary Conditions for Regional and Coastal Ocean Models. Journal of Atmospheric and Oceanic
!$ Technology, 23,1604-1613, http://dx.doi.org/10.1175/JTECH1930.1

      if(extmodtype.eq.'double')stop 'obc twin ext mode pas prevue'     !05/04/07

      if(txt_=='river') then !riverriverriverriver>
        call obc_ext_river 
        return !13-11-14
      endif                  !riverriverriverriver>

      if(flag_externalmode_lf_or_fb=='lf')call obc_ext_lf
      if(flag_externalmode_lf_or_fb=='fb')call obc_ext_fb

      stop 'obc_ext argument not recognized'
      end subroutine obc_ext

!.............................................................................

      subroutine obc_ext_lf
      use module_principal
      use module_parallele !#MPI
      implicit none
      stop 'obc_ext_lf not yet available'
      end subroutine obc_ext_lf

!.............................................................................

      subroutine obc_ext_fb
      use module_principal
      use module_parallele !#MPI
      use module_q
      implicit none
      integer t_,loop_
      integer :: idi_u_x_,idi_u_u1_,idi_u_u2_
      integer :: idi_v_y_,idi_v_v1_,idi_v_v2_

      if(extmodtype.eq.'double')stop 'obc twin ext mode pas prevue'     !05/04/07

! Reference solution:
      if(obcfreeorfix==0) then !rrrrrr>
       call obc_ext_reference_standard
      else                     !rrrrrr>
       call obc_ext_reference_adaptative
      endif                    !rrrrrr>



! Border conditions for tangential velocities:

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=2,imax
! note: fluxbar_u=Hdy*( stokes + u_euler )
           fluxbar_u(i,jmax  ,1)=(                            &
      velbarstokes_u(i,jmax  ,1)+                             &
! obc u_euler:
         velbarobc_u(i,jmax  ,1)                              &
           +velbar_u(i,jmax-1,2)                              &
        -velbarobc_u(i,jmax-1,1)                              &
                                )*mask_u(i,jmax,kmax)         &
       *(h_u(i,jmax)+0.5*(ssh_w(i,jmax,1)+ssh_w(i-1,jmax,1))) &
                                   *dy_u(i,jmax)

      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->
      do i=2,imax

! note: fluxbar_u=Hdy*( stokes + u_euler )
            fluxbar_u(i,1,1)=(                         &
       velbarstokes_u(i,1,1)+                          &
! obc u_euler:
          velbarobc_u(i,1,1)                           &
            +velbar_u(i,2,2)                           &
         -velbarobc_u(i,2,1)                           &
                             )*mask_u(i,1,kmax)        &
       *(h_u(i,1)+0.5*(ssh_w(i,1,1)+ssh_w(i-1,1,1)))   &
                                *dy_u(i,1)

      enddo
      endif                                          !----------->

! corners:
      if(obcstatus(ieqimax_jeqjmax)==1)  fluxbar_u(imax+1,jmax,1)=fluxbar_u(imax,jmax,1)
      if(obcstatus(ieq1_jeqjmax)==1)fluxbar_u(1     ,jmax,1)=fluxbar_u(2   ,jmax,1)
      if(obcstatus(ieqimax_jeq1)==1)   fluxbar_u(imax+1,1   ,1)=fluxbar_u(imax,1   ,1)
      if(obcstatus(ieq1_jeq1)==1) fluxbar_u(1     ,1   ,1)=fluxbar_u(2   ,1   ,1)


! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=2,jmax
! note: fluxbar_v=Hdx*( stokes + v_euler )
           fluxbar_v(imax  ,j,1)=(                             &
      velbarstokes_v(imax  ,j,1)+                              &
! obc v_euler:
         velbarobc_v(imax  ,j,1)                               &
           +velbar_v(imax-1,j,2)                               &
        -velbarobc_v(imax-1,j,1)                               &
                                 )*mask_v(imax,j,kmax)         &
       *(h_v(imax,j)+0.5*(ssh_w(imax,j,1)+ssh_w(imax,j-1,1)))  &
                                    *dx_v(imax  ,j)
      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->
      do j=2,jmax
! note: fluxbar_v=Hdx*( stokes + v_euler )
           fluxbar_v(1,j,1)=(                         &
      velbarstokes_v(1,j,1)+                          &
! obc v_euler:
         velbarobc_v(1,j,1)                           &
           +velbar_v(2,j,2)                           &
        -velbarobc_v(2,j,1)                           &
                           )*mask_v(1,j,kmax)         &
       *(h_v(1,j)+0.5*(ssh_w(1,j,1)+ssh_w(1,j-1,1)))  &
                              *dx_v(1,j)
      enddo
      endif                                          !----------->

! corners:
      if(obcstatus(ieqimax_jeqjmax)==1)  fluxbar_v(imax,jmax+1,1)=fluxbar_v(imax,jmax  ,1)
      if(obcstatus(ieqimax_jeq1)==1)   fluxbar_v(imax,1     ,1)=fluxbar_v(imax,2     ,1)
      if(obcstatus(ieq1_jeqjmax)==1)fluxbar_v(1   ,jmax+1,1)=fluxbar_v(1   ,jmax  ,1)
      if(obcstatus(ieq1_jeq1)==1) fluxbar_v(1   ,1     ,1)=fluxbar_v(1   ,2     ,1)


! Border conditions for sea surface height:

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=2,imax-1
                     ssh_w(i,jmax,2)=                           &
                                      sshrefobc_i(i,2)          &
                 +(                                             &
                    0.5*(velbar_v(i,jmax,1)+velbar_v(i,jmax,2)) &
                                     -vbrrefobc_i(i,2)          &
                                   )*sqr_hoverg_v(i,2)          !17-07-13

      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->

      do i=2,imax-1
                   ssh_w(i,1,2)=                               &
                               sshrefobc_i(i,1)                &
            -(                                                 &
               0.5*(velbar_v(i,2,1)+velbar_v(i,2,2))           &
                              -vbrrefobc_i(i,1)                &
                            )*sqr_hoverg_v(i,1)                !17-07-13

      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=2,jmax-1
                ssh_w(imax,j,2)=                               &
                                 sshrefobc_j(j,2)              &
           +(                                                  &
             0.5*(velbar_u(imax,j,1)+velbar_u(imax,j,2))       &
                                -vbrrefobc_j(j,2)              &
                              )*sqr_hoverg_u(j,2)              !17-07-13
      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->


      do j=2,jmax-1

                  ssh_w(1,j,2)=                                &
                             sshrefobc_j(j,1)                  &
           -(                                                  &
              0.5*(velbar_u(2,j,1)+velbar_u(2,j,2))            &
                            -vbrrefobc_j(j,1)                  &
                          )*sqr_hoverg_u(j,1)                  !17-07-13

      enddo
      endif                                          !----------->

! Corners:
      if(obcstatus(ieq1_jeq1)==1)  then !----------->
       i1=1    ; j1=1 ; i2=1 ; j2=2 ; i3=2 ; j3=1
       ssh_w(i1,j1,2)=0.5*(ssh_w(i2,j2,2)+ssh_w(i3,j3,2))
      endif                                          !----------->

      if(obcstatus(ieq1_jeqjmax)==1) then !----------->
       i1=1    ; j1=jmax ; i2=1 ; j2=jmax-1 ; i3=2 ; j3=jmax
       ssh_w(i1,j1,2)=0.5*(ssh_w(i2,j2,2)+ssh_w(i3,j3,2))
      endif                                          !----------->

      if(obcstatus(ieqimax_jeq1)==1)    then !----------->
       i1=imax ; j1=1 ; i2=imax ; j2=2 ; i3=imax-1 ; j3=1
       ssh_w(i1,j1,2)=0.5*(ssh_w(i2,j2,2)+ssh_w(i3,j3,2))
      endif                                          !----------->

      if(obcstatus(ieqimax_jeqjmax)==1)   then !----------->
       i1=imax ; j1=jmax ; i2=imax ; j2=jmax-1 ; i3=imax-1 ; j3=jmax
       ssh_w(i1,j1,2)=0.5*(ssh_w(i2,j2,2)+ssh_w(i3,j3,2))
      endif                                          !----------->

! Border conditions for normal velocities by inversing the continuity
! equation and using ssh and tangential velocities:

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=2,imax-1
      fluxbar_v(i,jmax+1,1)= (                                        &

      (velbar_v(i,jmax,2)+velbarstokes_v(i,jmax,1))                   & 
          *dx_v(i,jmax)                                               &
          *(h_v(i,jmax)+0.5*(ssh_w(i,jmax,1)+ssh_w(i,jmax-1,1)))      &

       -( fluxbar_u(i+1,jmax,1)-fluxbar_u(i,jmax,1) +                 &

        ((ssh_w(i,jmax,2)-ssh_w(i,jmax,1))/dte_fw                     &
                       +omega_w(i,jmax,kmax+1,1))*dxdy_t(i,jmax))     &

                             )*mask_v(i,jmax+1,kmax)

      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->
      do i=2,imax-1
      fluxbar_v(i,1,1) = (                                           &

      (velbar_v(i,2,2)+velbarstokes_v(i,2,1))                        &
          *dx_v(i,2)                                                 &
          *(h_v(i,2)+0.5*(ssh_w(i,2,1)+ssh_w(i,1,1)))                &

        +(fluxbar_u(i+1,1,1)-fluxbar_u(i,1,1) +                      &

        ((ssh_w(i,1,2)-ssh_w(i,1,1))/dte_fw                          &
                    +omega_w(i,1,kmax+1,1))*dxdy_t(i,1))             &

                         )*mask_v(i,1,kmax)

      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=2,jmax-1
      fluxbar_u(imax+1,j,1) = (                                       &

      (velbar_u(imax,j,2)+velbarstokes_u(imax,j,1))                   &
          *dy_u(imax,j)                                               &
          *(h_u(imax,j)+0.5*(ssh_w(imax,j,1)+ssh_w(imax-1,j,1)))      &

        -(fluxbar_v(imax,j+1,1)-fluxbar_v(imax,j,1) +                 &

        ((ssh_w(imax,j,2)-ssh_w(imax,j,1))/dte_fw                     &
                       +omega_w(imax,j,kmax+1,1))*dxdy_t(imax,j))     &

                              )*mask_u(imax+1,j,kmax)

      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->
      do j=2,jmax-1
      fluxbar_u(1,j,1) = (                                          &

      (velbar_u(2,j,2)+velbarstokes_u(2,j,1))                       &
          *dy_u(2,j)                                                &
          *(h_u(2,j)+0.5*(ssh_w(2,j,1)+ssh_w(1,j,1)))               &

        +(fluxbar_v(1,j+1,1)-fluxbar_v(1,j,1)+                      &

        ((ssh_w(1,j,2)-ssh_w(1,j,1))/dte_fw                         &
                    +omega_w(1,j,kmax+1,1))*dxdy_t(1,j))            &

                         )*mask_u(1,j,kmax)

      enddo
      endif                                          !----------->

! corners:
                                                                       !11/06/04
      if(obcstatus(ieqimax_jeqjmax)==1)   then !----------->
      i=imax
      j=jmax
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_w(i,j,2)- ssh_w(i,j,1))/dte_fw*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i+1,j  ,1))                                &
       +abs(fluxbar_v(i  ,j+1,1))                                &
                )

      fluxbar_u(i+1,j  ,1)=fluxbar_u(i+1,j  ,1)                  &
                   +x1*abs(fluxbar_u(i+1,j  ,1))                 &
                             *mask_u(i+1,j  ,kmax)
                                                                       !11/06/04
      fluxbar_v(i  ,j+1,1)=fluxbar_v(i  ,j+1,1)                  &
                   +x1*abs(fluxbar_v(i  ,j+1,1))                 &
                             *mask_v(i  ,j+1,kmax)

      endif                                          !----------->

      if(obcstatus(ieqimax_jeq1)==1)    then !----------->
      i=imax
      j=1
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_w(i,j,2)- ssh_w(i,j,1))/dte_fw*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i+1,j  ,1))                                &
       +abs(fluxbar_v(i  ,j  ,1))                                &
                )

      fluxbar_u(i+1,j  ,1)=fluxbar_u(i+1,j  ,1)                  &
                   +x1*abs(fluxbar_u(i+1,j  ,1))                 &
                             *mask_u(i+1,j  ,kmax)
                                                                       !11/06/04
      fluxbar_v(i  ,j  ,1)=fluxbar_v(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_v(i  ,j  ,1))                 &
                             *mask_v(i  ,j  ,kmax)

      endif                                          !----------->

      if(obcstatus(ieq1_jeqjmax)==1) then !----------->
      i=1
      j=jmax
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_w(i,j,2)- ssh_w(i,j,1))/dte_fw*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i  ,j  ,1))                                &
       +abs(fluxbar_v(i  ,j+1,1))                                &
                )

      fluxbar_u(i  ,j  ,1)=fluxbar_u(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_u(i  ,j  ,1))                 &
                             *mask_u(i  ,j  ,kmax)
                                                                       !11/06/04
      fluxbar_v(i  ,j+1,1)=fluxbar_v(i  ,j+1,1)                  &
                   +x1*abs(fluxbar_v(i  ,j+1,1))                 &
                             *mask_v(i  ,j+1,kmax)

      endif                                          !----------->

      if(obcstatus(ieq1_jeq1)==1)  then !----------->
      i=1
      j=1
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_w(i,j,2)- ssh_w(i,j,1))/dte_fw*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i  ,j  ,1))                                &
       +abs(fluxbar_v(i  ,j  ,1))                                &
                )

      fluxbar_u(i  ,j  ,1)=fluxbar_u(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_u(i  ,j  ,1))                 &
                             *mask_u(i  ,j  ,kmax)
                                                                      !11/06/04
      fluxbar_v(i  ,j  ,1)=fluxbar_v(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_v(i  ,j  ,1))                 &
                             *mask_v(i  ,j  ,kmax)


      endif                                          !----------->

!                         / / /


! ADJUST fluxbar to possible other requirements !20-03-17
!     call obc_ext_gibraltarflux
!     call obc_ext_sicileflux_i1 ! avec reservoir !


! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       i=1
        do j=1,jmax+1
             velbar_v(i,j,2)=(                                &
            fluxbar_v(i,j,1)                                  &
                /(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))) &
                /dx_v(i,j)                                    &
      -velbarstokes_v(i,j,1)                                  &
             )*mask_v(i,j,kmax)                          !01-03-10
        enddo

       i=1
        do j=1,jmax
             velbar_u(i,j,2)=(                                &
            fluxbar_u(i,j,1)                                  &
                /(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))) &
                /dy_u(i,j)                                    &
      -velbarstokes_u(i,j,1)                                  &
             )*mask_u(i,j,kmax)                          !01-03-10
        enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       j=1
       do i=1,imax+1
             velbar_u(i,j,2)=(                                &
            fluxbar_u(i,j,1)                                  &
                /(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))) &
                /dy_u(i,j)                                    &
      -velbarstokes_u(i,j,1)                                  &
             )*mask_u(i,j,kmax)                          !01-03-10
       enddo

       j=1
       do i=1,imax
             velbar_v(i,j,2)=(                                &
            fluxbar_v(i,j,1)                                  &
                /(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))) &
                /dx_v(i,j)                                    &
      -velbarstokes_v(i,j,1)                                  &
             )*mask_v(i,j,kmax)                          !01-03-10
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       i=imax
        do j=1,jmax+1
             velbar_v(i,j,2)=(                                &
            fluxbar_v(i,j,1)                                  &
                /(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))) &
                /dx_v(i,j)                                    &
      -velbarstokes_v(i,j,1)                                  &
             )*mask_v(i,j,kmax)                          !01-03-10
        enddo

       i=imax+1
        do j=1,jmax
             velbar_u(i,j,2)=(                                &
            fluxbar_u(i,j,1)                                  &
                /(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))) &
                /dy_u(i,j)                                    &
      -velbarstokes_u(i,j,1)                                  &
             )*mask_u(i,j,kmax)                          !01-03-10
        enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       j=jmax
       do i=1,imax+1
             velbar_u(i,j,2)=(                                &
            fluxbar_u(i,j,1)                                  &
                /(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))) &
                /dy_u(i,j)                                    &
      -velbarstokes_u(i,j,1)                                  &
             )*mask_u(i,j,kmax)                          !01-03-10
       enddo

       j=jmax+1
       do i=1,imax
             velbar_v(i,j,2)=(                                &
            fluxbar_v(i,j,1)                                  &
                /(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))) &
                /dx_v(i,j)                                    &
      -velbarstokes_v(i,j,1)                                  &
             )*mask_v(i,j,kmax)                          !01-03-10
       enddo
      endif                                      !----------->

! La condition aux limites sur les rivieres vient obligatoirement apres
! la conversion des flux en vitesses sinon les rivieres aux bords du domaine
! sont ecras�s par la conversion...

! Lignes et colonnes supplementaires pour advection o4:
      call obc_ext_o4 !25-05-17 

! Echange mpi u1 u2 v1 v2 sur velbar_u(:,:,2) et velbar_v(:,:,2)
      call obc_ext_mpi(2) !02-07-14

      end subroutine obc_ext_fb

!..............................................................

      subroutine obc_ext_o4
      use module_principal
      use module_parallele
      implicit none

! Lignes colonnes supplementaires pour advection O4

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       do j=1,jmax
        velbar_u(imax+2,j,2)=velbar_u(imax+1,j,2)
       enddo
       do j=0,jmax+2
        velbar_v(imax+1,j,2)=velbar_v(imax,j,2)
       enddo
      endif                                      !----------->

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       do j=1,jmax
        velbar_u(0     ,j,2)=velbar_u(1     ,j,2)
       enddo
       do j=0,jmax+2
        velbar_v(0     ,j,2)=velbar_v(1   ,j,2)
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do i=0,imax+2
        velbar_u(i,jmax+1,2)=velbar_u(i,jmax,2)
       enddo
       do i=1,imax
        velbar_v(i,jmax+2,2)=velbar_v(i,jmax+1,2)
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do i=0,imax+2
        velbar_u(i,0     ,2)=velbar_u(i,1   ,2)
       enddo
       do i=1,imax
        velbar_v(i,0     ,2)=velbar_v(i,1     ,2)
       enddo
      endif                                      !----------->


      end subroutine obc_ext_o4

!..............................................................

      subroutine obc_ext_mpi3d !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer loop_,idi_u_u1_,idi_u_u2_,idi_v_v1_,idi_v_v2_

!!$! Nouvelle methode avec choix des voisins
      call get_type_echange('u1','velbar_u1_3d'  &
                                 ,velbar_u       &
                          ,lbound(velbar_u)      &
                          ,ubound(velbar_u),idi_u_u1_)

      call get_type_echange('u2','velbar_u2_3d'  &
                                 ,velbar_u       &
                          ,lbound(velbar_u)      &
                          ,ubound(velbar_u),idi_u_u2_)

      call get_type_echange('v1','velbar_v1_3d'  &
                                 ,velbar_v       &
                          ,lbound(velbar_v)      &
                          ,ubound(velbar_v),idi_v_v1_)

      call get_type_echange('v2','velbar_v2_3d'  &
                                 ,velbar_v       &
                          ,lbound(velbar_v)      &
                          ,ubound(velbar_v),idi_v_v2_)

      ! Echanges
      !print *,"obc_ext : ",par%rank,subcycle_exchange
      do loop_=1, subcycle_exchange
         call echange_voisin(velbar_u,idi_u_u1_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_u,idi_u_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,idi_v_v1_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,idi_v_v2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ext_mpi3d

!..............................................................

      subroutine obc_ext_mpi(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_u1_,idi_u2_,idi_v1_,idi_v2_


!!$! Nouvelle methode avec choix des voisins
      write(texte30,'(a10,i0)')'velbar_u1_',t_
      call get_type_echange('u1',trim(texte30),velbar_u      &
                                       ,lbound(velbar_u)     &
                                       ,ubound(velbar_u)     &
                                       ,t_                   &
                                       ,idi_u1_)

      write(texte30,'(a10,i0)')'velbar_u2_',t_
      call get_type_echange('u2',trim(texte30),velbar_u      &
                                       ,lbound(velbar_u)     &
                                       ,ubound(velbar_u)     &
                                       ,t_                   &
                                       ,idi_u2_)

      write(texte30,'(a10,i0)')'velbar_v1_',t_
      call get_type_echange('v1',trim(texte30),velbar_v      &
                                       ,lbound(velbar_v)     &
                                       ,ubound(velbar_v)     &
                                       ,t_                   &
                                       ,idi_v1_)

      write(texte30,'(a10,i0)')'velbar_v2_',t_
      call get_type_echange('v2',trim(texte30),velbar_v      &
                                       ,lbound(velbar_v)     &
                                       ,ubound(velbar_v)     &
                                       ,t_                   &
                                       ,idi_v2_)

      ! Echanges
      !print *,"obc_ext : ",par%rank,subcycle_exchange
      do loop_=1,subcycle_exchange
         call echange_voisin(velbar_u,idi_u1_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_u,idi_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,idi_v1_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,idi_v2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges


      end subroutine obc_ext_mpi

!..............................................................

      subroutine obc_ext_xy_t(txt_,t_) !13-10-15
      use module_principal ; use module_parallele
      implicit none ; integer t_,loop_ ; character*2 txt_

      write(texte30,'(a5,a2,a1,i0)')'xy_t_',txt_,'_',t_

      call get_type_echange(txt_,trim(texte30),xy_t      &
                                       ,lbound(xy_t)     &
                                       ,ubound(xy_t)     &
                                       ,t_               &
                                       ,k0)
      do loop_=1,subcycle_exchange
         call echange_voisin(xy_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ext_xy_t

!..............................................................

      subroutine obc_ext_mpi_xy(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_u1_,idi_u2_,idi_v1_,idi_v2_


!!$! Nouvelle methode avec choix des voisins
      write(texte30,'(a6,i0)')'xy_u1_',t_
      call get_type_echange('u1',trim(texte30),xy_u      &
                                       ,lbound(xy_u)     &
                                       ,ubound(xy_u)     &
                                       ,t_               &
                                       ,idi_u1_)

      write(texte30,'(a6,i0)')'xy_u2_',t_
      call get_type_echange('u2',trim(texte30),xy_u      &
                                       ,lbound(xy_u)     &
                                       ,ubound(xy_u)     &
                                       ,t_               &
                                       ,idi_u2_)

      write(texte30,'(a6,i0)')'xy_v1_',t_
      call get_type_echange('v1',trim(texte30),xy_v      &
                                       ,lbound(xy_v)     &
                                       ,ubound(xy_v)     &
                                       ,t_               &
                                       ,idi_v1_)

      write(texte30,'(a6,i0)')'xy_v2_',t_
      call get_type_echange('v2',trim(texte30),xy_v      &
                                       ,lbound(xy_v)     &
                                       ,ubound(xy_v)     &
                                       ,t_               &
                                       ,idi_v2_)

      ! Echanges
      !print *,"obc_ext : ",par%rank,subcycle_exchange
      do loop_=1,subcycle_exchange
         call echange_voisin(xy_u,idi_u1_,mpi_neighbor_list(loop_))
         call echange_voisin(xy_u,idi_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(xy_v,idi_v1_,mpi_neighbor_list(loop_))
         call echange_voisin(xy_v,idi_v2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ext_mpi_xy

!..............................................................

      subroutine obc_ext_velbar_mpi(t_) !15-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_x_,idi_y_


!!$! Nouvelle methode avec choix des voisins
      write(texte30,'(a11,i0)')'velbar_u_x_',t_
      call get_type_echange('u1',trim(texte30),velbar_u      &
                                       ,lbound(velbar_u)     &
                                       ,ubound(velbar_u)     &
                                       ,t_                   &
                                       ,idi_x_)

      write(texte30,'(a11,i0)')'velbar_v_y_',t_
      call get_type_echange('v1',trim(texte30),velbar_v      &
                                       ,lbound(velbar_v)     &
                                       ,ubound(velbar_v)     &
                                       ,t_                   &
                                       ,idi_y_)

      ! Echanges
      !print *,"obc_ext : ",par%rank,subcycle_exchange
      do loop_=1,subcycle_exchange
         call echange_voisin(velbar_u,idi_x_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,idi_y_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges


      end subroutine obc_ext_velbar_mpi

!..............................................................


      subroutine obc_ext_river
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_x_,idi_y_

!$ Begin the Boundary conditions at the river grid nodes:

      if(flag_externalmode_lf_or_fb=='lf')t_=0
      if(flag_externalmode_lf_or_fb=='fb')t_=1 !01-05-13

      if(iteration3d==2) then !DEBUG>
! Verification du reset du cumul de fluxbar !02-10-19 suite A la
! possibilite de cumuler plusieurs fleuves sur un mEme point d'entrEe
! Voir aussi: https://docs.google.com/document/d/1d0NAgA5j5L7KfwBB41FQ8fsd4E1Z9gQS0kU44np7Vbo/edit
      flag_stop=0
      do 200 kr=1,nriver
      if(riverdir(kr)>0)then!>>>>>>>>>>>>>>>>>>>>>>>>>>>   
      if(rivervel_inout(kr)==1) then !rrrrrrrrrrrrr>
      if(mod(riverdir(kr),2).eq.0) then ! ----------------->
         if(fluxbar_v(iriver(kr,2),jriver(kr,2),1)/=0)flag_stop=1
      else
         if(fluxbar_u(iriver(kr,2),jriver(kr,2),1)/=0)flag_stop=1
      endif !---------------------------------------------->
      endif                          !rrrrrrrrrrrrr>

      endif                    !>>>>>>>>>>>>>>>>>>>>>>>>>>>   
  200 continue
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop 'Valeur nulle inattendue pour fluxbar obc_ext_river'
      endif                   !DEBUG

      do 100 kr=1,nriver

      if(riverdir(kr)>0)then!>>>>>>>>>>>>>>>>>>>>>>>>>>>            !25/11/01 !09-04-19

      if(rivertrc_inout(kr)==1)ssh_w(iriver(kr,1),jriver(kr,1),t_)= & !22-01-15
                 max(wetdry_cst2-h_w(iriver(kr,1),jriver(kr,1)),zero)

      if(rivervel_inout(kr)==1) then !rrrrrrrrrrrrr>

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->
         fluxbar_v(iriver(kr,2),jriver(kr,2),1)= &
         fluxbar_v(iriver(kr,2),jriver(kr,2),1)+riverflux(kr,1) !02-10-19
! Note la sommation fluxbar_v=fluxbar_v+... permet de cumuler les flux de plusieurs
! rivieres partageant un meme point d'entree. Les resctrictions concernant cette option
! sont discutees dans https://docs.google.com/document/d/1d0NAgA5j5L7KfwBB41FQ8fsd4E1Z9gQS0kU44np7Vbo/edit
! Un test de verification que le cumul demarre bien de 0 A chaque iteration a EtE ajoutE avant la
! boucle "do 100 kr=1,nriver" (voir plus haut)

! Pour le moment on assume velbar=0
!         velbar_v(iriver(kr,2),jriver(kr,2),1)=riverflux(kr,1)        &
!      /       h_v(iriver(kr,2),jriver(kr,2))                          &
!            /dx_v(iriver(kr,2),jriver(kr,2))

      else
         fluxbar_u(iriver(kr,2),jriver(kr,2),1)= &
         fluxbar_u(iriver(kr,2),jriver(kr,2),1)+riverflux(kr,1) !02-10-19
! Pour le moment on assume velbar=0
!         velbar_u(iriver(kr,2),jriver(kr,2),1)=riverflux(kr,1)        &
!      /       h_u(iriver(kr,2),jriver(kr,2))                          &
!            /dy_u(iriver(kr,2),jriver(kr,2))
      endif !---------------------------------------------->
      endif                          !rrrrrrrrrrrrr>

      endif                    !>>>>>>>>>>>>>>>>>>>>>>>>>>>            !25/11/01

  100 continue

!$ End of the Boundary conditions at the river grid nodes:

      end subroutine obc_ext_river
!..................................................................................

      subroutine obc_ext_reference_standard
      use module_principal
      use module_parallele
      implicit none
!     integer t_,loop_,idi_x_,idi_y_

! const1 et const2 assurent une variation lineaire du forcage          !05/04/07

      const1=real(iteration2d_max_now-iteration2d)/real(iteration2d_max_now-iteration2d+1)
      const2=1.-const1

      do i=2,imax-1

! FRONTIERE J=NECO:
      sshrefobc_i(i,2)=const1*sshrefobc_i(i,2)                          &
                      +const2*   sshobc_w(i,jmax,1)
      vbrrefobc_i(i,2)=const1*vbrrefobc_i(i,2)                          &
                      +const2*velbarobc_v(i,jmax,1)

! FRONTIERE J=1:
      sshrefobc_i(i,1)=const1*sshrefobc_i(i,1)                          &
                      +const2*   sshobc_w(i,1,1)
      vbrrefobc_i(i,1)=const1*vbrrefobc_i(i,1)                          &
                      +const2*velbarobc_v(i,2,1)

      enddo

      do j=2,jmax-1

! FRONTIERE I=MECO:
      sshrefobc_j(j,2)=const1*    sshrefobc_j(j,2)                      &
                      +const2*  sshobc_w(imax,j,1)
      vbrrefobc_j(j,2)=const1*    vbrrefobc_j(j,2)                      &
                      +const2*velbarobc_u(imax,j,1)

! FRONTIERE I=1:
      sshrefobc_j(j,1)=const1*sshrefobc_j(j,1)                          &
                      +const2* sshobc_w(1,j,1)
      vbrrefobc_j(j,1)=const1*vbrrefobc_j(j,1)                          &
                      +const2*velbarobc_u(2,j,1)

      enddo

      end subroutine obc_ext_reference_standard

!..................................................................................

      subroutine obc_ext_reference_adaptative
      use module_principal
      use module_parallele
      implicit none
!     integer t_,loop_,idi_x_,idi_y_
! On ne passe qu'une fois par iteration barocline donc on prend DTI_FW au lieu
! de DTE/2
      if(iteration2d.eq.1) then !>>>>>>>>>>>>>>>>>>>>>

      x1=2.*pi/coriolis_t(imax/2,jmax/2)*1.  ! duree lissage (ici 1 periodes inertielles)
      const2=dti_fw/x1
      const1=1.-const2

      do i=2,imax-1

! FRONTIERE J=NECO:
      sshrefobc_i(i,2)=const1*sshrefobc_i(i,2)        &
                      +const2*ssh_w(i,jmax-1,2)  
      vbrrefobc_i(i,2)=const1*vbrrefobc_i(i,2)        &
                      +const2*   velbar_v(i,jmax-1,2)                  !13/12/07

! FRONTIERE J=1:
      sshrefobc_i(i,1)=const1*sshrefobc_i(i,1)        &
                      +const2*ssh_w(i,2,2)
      vbrrefobc_i(i,1)=const1*vbrrefobc_i(i,1)        &
                      +const2*   velbar_v(i,3,2)                       !13/12/07

      enddo

      do j=2,jmax-1

! FRONTIERE I=MECO:
      sshrefobc_j(j,2)=const1*sshrefobc_j(j,2)   &
                      +const2*ssh_w(imax-1,j,2)    !10/11/05
      vbrrefobc_j(j,2)=const1*     vbrrefobc_j(j,2)   &
                      +const2* velbar_u(imax-1,j,2)                    !13/12/07

! FRONTIERE I=1:
      sshrefobc_j(j,1)=const1*sshrefobc_j(j,1)        &
                      +const2*ssh_w(2,j,2)
      vbrrefobc_j(j,1)=const1*vbrrefobc_j(j,1)        &
                      +const2* velbar_u(3,j,2)                         !13/12/07

      enddo

! Les 4 coins de sse pour pour pouvoir traiter la OBC pour T et S
! coin i=1 j=1
      sshrefobc_i(1,1)=0.5*(sshrefobc_i(2,1)+sshrefobc_j(2,1))
      sshrefobc_j(1,1)=sshrefobc_i(1,1)
! coin i=1 j=jmax
      sshrefobc_i(1   ,2)=0.5*(sshrefobc_i(2,2)+sshrefobc_j(jmax-1,1))
      sshrefobc_j(jmax,1)=sshrefobc_i(1,2)
! coin i=imax j=jmax
      sshrefobc_i(imax,2)=0.5*(sshrefobc_i(imax-1,2)                    &
                              +sshrefobc_j(jmax-1,2))
      sshrefobc_j(jmax,2)=sshrefobc_i(imax,2)
! coin i=imax j=1
      sshrefobc_i(imax,1)=0.5*(sshrefobc_i(imax-1,1)+sshrefobc_j(2,2))
      sshrefobc_j(1   ,2)=sshrefobc_i(imax,1)

      endif                !>>>>>>>>>>>>>>>>>>>>>

      end subroutine obc_ext_reference_adaptative

!..................................................................................

      subroutine obc_ext_sshmean_control
      use module_principal
      use module_parallele
      implicit none
      STOP'la conservation de la ssh estparallelisee maispas verifiee'
      end subroutine obc_ext_sshmean_control
!..................................................................................

!..............................................................

      subroutine obc_ext_egrid_xyf(exchtype_,indx1_,indx2_) !08-02-15
      use module_principal ; use module_parallele
      implicit none
      integer loop_,indx1_,indx2_,exchid1_,exchid2_
      character*2 exchtype_ ! 'r ' ou 'r1'

      write(texte30,'(a,i0)')'xy_f_'//exchtype_//'_',indx1_
      call get_type_echange('r1'              &
                            ,trim(texte30)    &
                                   ,xy_f      &
                            ,lbound(xy_f)     &
                            ,ubound(xy_f)     &
                            ,indx1_           &
                            ,exchid1_)


      write(texte30,'(a,i0)')'xy_f_'//exchtype_//'_',indx2_
      call get_type_echange( exchtype_        &
                            ,trim(texte30)    &
                                   ,xy_f      &
                            ,lbound(xy_f)     &
                            ,ubound(xy_f)     &
                            ,indx2_           &
                            ,exchid2_)

      do loop_=1,subcycle_exchange
         call echange_voisin(xy_f,exchid1_,mpi_neighbor_list(loop_))
         call echange_voisin(xy_f,exchid2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ext_egrid_xyf

!..............................................................

      subroutine obc_ext_qdmflux !18-08-15
      use module_principal ; use module_parallele
      implicit none

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)                  then !----------->
      do i=2,imax-1
! Note: on zerograde sur jmax-2 car flux(jmax-1) sous influence obc o4 (donc danger)
       yflux_t(i,jmax  )=yflux_t(i,jmax-2)
       yflux_t(i,jmax-1)=yflux_t(i,jmax-2)
      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)                    then !----------->
      do i=2,imax-1
! Note: on zerograde sur j=3 car flux(j=2) sous influence obc o4 (donc danger)
       yflux_t(i,1)=yflux_t(i,3)
       yflux_t(i,2)=yflux_t(i,3)
      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)                 then !----------->
      do j=2,jmax-1
! Note: on zerograde sur imax-2 car flux(imax-1) sous influence obc o4 (donc danger)
       xflux_t(imax  ,j)=xflux_t(imax-2,j)
       xflux_t(imax-1,j)=xflux_t(imax-2,j)

!     xflux_t(imax,j)=-0.25*( & !oooo>
!         (fluxbar_u(imax,j,0)+fluxbar_u(imax+1,j,0))   &
!     *( velbarobc_u(imax,j,1)+velbar_u(imax,j,before)) &
!     +abs(fluxbar_u(imax,j,0)+fluxbar_u(imax+1,j,0))   &
!     *(-velbarobc_u(imax,j,1)+velbar_u(imax,j,before)) &
!                           )   !oooo>
      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)                    then !----------->
      do j=2,jmax-1
! Note: on zerograde sur i=3 car flux(i=2) sous influence obc o4 (donc danger)
       xflux_t(1,j)=xflux_t(3,j)
       xflux_t(2,j)=xflux_t(3,j)

!     xflux_t(1,j)=-0.25*( & !oooo>
!         (fluxbar_u(1,j,0)+fluxbar_u(2,j,0))      &
!      *(velbarobc_u(2,j,1) +velbar_u(2,j,before)) &
!     +abs(fluxbar_u(1,j,0)+fluxbar_u(2,j,0))      &
!      *(velbarobc_u(2,j,1) -velbar_u(2,j,before))  &
!                        )   !oooo>
           
      enddo
      endif                                          !----------->


      end subroutine obc_ext_qdmflux

!..............................................................

      subroutine obc_ext_mpi_xyt(t_,txt_) !30-07-15
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_
      character*2 txt_


!!$! Nouvelle methode avec choix des voisins
      write(texte30,'(a6,a2,i0)')'xy_t_',txt_,t_
      call get_type_echange(txt_,trim(texte30),xy_t      &
                                       ,lbound(xy_t)     &
                                       ,ubound(xy_t)     &
                                       ,t_               &
                                       ,k0)

      ! Echanges
      !print *,"obc_ext : ",par%rank,subcycle_exchange
      do loop_=1,subcycle_exchange
         call echange_voisin(xy_t,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ext_mpi_xyt

!..............................................................

      subroutine obc_ext_o2o4 !20-08-15
      use module_principal
      use module_parallele
      implicit none

! Lignes colonnes supplementaires pour advection O4

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       do j=1,jmax
        velbar_u(imax+2,j,1)=velbar_u(imax  ,j,1)
        velbar_u(imax+1,j,1)=velbar_u(imax  ,j,1)
       enddo
       do j=0,jmax+2
        velbar_v(imax+1,j,1)=velbar_v(imax-1,j,1)
        velbar_v(imax  ,j,1)=velbar_v(imax-1,j,1)
       enddo
      endif                                      !----------->

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       do j=1,jmax
        velbar_u(0,j,1)=velbar_u(2,j,1)
        velbar_u(1,j,1)=velbar_u(2,j,1)
       enddo
       do j=0,jmax+2
        velbar_v(0,j,1)=velbar_v(2,j,1)
        velbar_v(1,j,1)=velbar_v(2,j,1)
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do i=0,imax+2
        velbar_u(i,jmax+1,1)=velbar_u(i,jmax-1,1)
        velbar_u(i,jmax  ,1)=velbar_u(i,jmax-1,1)
       enddo
       do i=1,imax
        velbar_v(i,jmax+2,1)=velbar_v(i,jmax  ,1)
        velbar_v(i,jmax+1,1)=velbar_v(i,jmax  ,1)
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do i=0,imax+2
        velbar_u(i,0,1)=velbar_u(i,2,1)
        velbar_u(i,1,1)=velbar_u(i,2,1)
       enddo
       do i=1,imax
        velbar_v(i,0,1)=velbar_v(i,2,1)
        velbar_v(i,1,1)=velbar_v(i,2,1)
       enddo
      endif                                      !----------->


      end subroutine obc_ext_o2o4

!..............................................................

      subroutine obc_ext_gibraltarflux !20-03-17
      use module_principal
      use module_parallele
      implicit none
! D�tails dans:
! https://docs.google.com/document/d/1_jIrcNnit1pSpq-8FyrIa3EDSf3WWP6XMidmXQsFTGo/edit

      if(kmaxtide>0)stop 'Cas Maree pas encore prevu'
! Quand la maree sera introduite, il faudra calculer le transport
! qui lui est associE afin de ne pas le compter dans la correction 
! du bilan pour que la maree ne soit pas impactEe par la contrainte
! sur le transport

      if(obcstatus(ieq1)==1) then !----->              
       sum0=0. ; sum1=0.
       do j=1,jmax
        sum0=sum0+mask_j_u(j)*mask_u(1,j,kmax)      *h_u(1,j)
        sum1=sum1+mask_j_u(j)*mask_u(1,j,kmax)*fluxbar_u(1,j,1)
       enddo
      else                        !----->
       sum0=0. ; sum1=0.
      endif                       !----->
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)

! ENTRER LE BILAN DE FLUX IDEAL:
! idealflux est le flux total (m3/s) 
!                                          idealflux=0.05d6.
                                           idealflux=0.
!     if(iteration3d==0.and.iteration2d==1)idealflux=sum1glb

! COEF DE LAGRANGE: ENTRER L'ECHELLE DE TEMPS POUR UNE CORRECTION EN TEMPS DIFFERE
      x0=(idealflux-sum1glb+cumuldeltaflux)/sum0glb & ! Correction en une iteration
         *dte_fw/86400.                               ! fractionnement sur temps long
!        *1.                                          ! fractionnement sur temps long

      if(obcstatus(ieq1)==1) then !----->              
       sum1=0.
       do j=1,jmax
        fluxbar_u(1,j,1)=fluxbar_u(1,j,1)+x0*mask_u(1,j,kmax)*h_u(1,j)
        sum1=sum1+mask_j_u(j)*mask_u(1,j,kmax)*fluxbar_u(1,j,1)
       enddo
      else                        !----->
       sum1=0.
      endif                       !----->
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      cumuldeltaflux=cumuldeltaflux+(idealflux-sum1glb)

      if(par%rank==0) then
      if(iteration2d==1) then
       open(unit=3,file='tmp/flux_gibraltar_obc',position='append')
          write(3,*)real(elapsedtime_now/86400.),real(sum1glb),real(cumuldeltaflux)
       close(3)
      endif
      endif

      end subroutine obc_ext_gibraltarflux

!..................................................................................

      subroutine obc_ext_sicileflux_imax !20-03-17
      use module_principal
      use module_parallele
      implicit none

! D�tails dans:
! https://docs.google.com/document/d/1_jIrcNnit1pSpq-8FyrIa3EDSf3WWP6XMidmXQsFTGo/edit

      if(kmaxtide>0)stop 'Cas Maree pas encore prevu'
! Quand la maree sera introduite, il faudra calculer le transport
! qui lui est associE afin de ne pas le compter dans la correction 
! du bilan pour que la maree ne soit pas impactEe par la contrainte
! sur le transport

      if(obcstatus(ieq1)==1) then !----->              
       sum0=0. ; sum1=0.
       do j=1,jmax
        sum0=sum0+mask_j_u(j)*mask_u(imax+1,j,kmax)      *h_u(imax+1,j)
        sum1=sum1+mask_j_u(j)*mask_u(imax+1,j,kmax)*fluxbar_u(imax+1,j,1)
       enddo
      else                        !----->
       sum0=0. ; sum1=0.
      endif                       !----->
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)

! ENTRER LE BILAN DE FLUX IDEAL:
! idealflux est le flux total (m3/s) (positif vers l'est)
!                                          idealflux=0.05d6.
                                           idealflux=0.
!     if(iteration3d==0.and.iteration2d==1)idealflux=sum1glb

! COEF DE LAGRANGE: ENTRER L'ECHELLE DE TEMPS POUR UNE CORRECTION EN TEMPS DIFFERE
      x0=(idealflux-sum1glb+cumuldeltaflux)/sum0glb & ! Correction en une iteration
         *dte_fw/86400.                               ! fractionnement sur temps long
!        *0.                                          ! fractionnement sur temps long

      if(obcstatus(ieq1)==1) then !----->              
       sum1=0.
       do j=1,jmax
        fluxbar_u(imax+1,j,1)=fluxbar_u(imax+1,j,1)+x0*mask_u(imax+1,j,kmax)*h_u(imax+1,j)
        sum1=sum1+mask_j_u(j)*mask_u(imax+1,j,kmax)*fluxbar_u(imax+1,j,1)
       enddo
      else                        !----->
       sum1=0.
      endif                       !----->
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      cumuldeltaflux=cumuldeltaflux+(idealflux-sum1glb)

      if(par%rank==0) then
      if(iteration2d==1) then
       open(unit=3,file='tmp/flux_sicile_obc',position='append')
          write(3,*)real(elapsedtime_now/86400.),real(sum1glb),real(cumuldeltaflux)
       close(3)
      endif
      endif

      end subroutine obc_ext_sicileflux_imax

!..................................................................................

      subroutine obc_ext_sicileflux_i1 !20-03-17
      use module_principal
      use module_parallele
      implicit none
! D�tails dans:
! https://docs.google.com/document/d/1_jIrcNnit1pSpq-8FyrIa3EDSf3WWP6XMidmXQsFTGo/edit

      signe=-1      ! Convention de signe pour obc en i ou j =1
       area=1.67d12 ! Superficie reservoir Oriental

      if(kmaxtide>0)stop 'Cas Maree pas encore prevu'
! Quand la maree sera introduite, il faudra calculer le transport
! qui lui est associE afin de ne pas le compter dans la correction 
! du bilan pour que la maree ne soit pas impactEe par la contrainte
! sur le transport

      if(obcstatus(jeq1)==1) then !----->              
       sum0=0. ; sum1=0.
       do i=1,imax
        sum0=sum0+mask_i_v(i)*mask_v(i,1,kmax)      *h_v(i,1)
        sum1=sum1+mask_i_v(i)*mask_v(i,1,kmax)*fluxbar_v(i,1,1)
       enddo
      else                        !----->
       sum0=0. ; sum1=0.
      endif                       !----->
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)

! ENTRER LE BILAN DE FLUX IDEAL:
! idealflux est le flux total (m3/s) 
!                                          idealflux=0.05d6.
                                           idealflux=0.
!     if(iteration3d==0.and.iteration2d==1)idealflux=sum1glb
! Methode des flux:
! Idealflux equilibre le bilan du flux ecmwf + rivieres et variation ssh moyenne du bassin oriental
!     x1=1.e-5/86400. ! x1=d<ssh>/dt du bassin oriental
!     idealflux=signe*area*x1 ! +signe*area*SurfaceRiverFlux !21-03-17

! COEF DE LAGRANGE: ENTRER L'ECHELLE DE TEMPS POUR UNE CORRECTION EN TEMPS DIFFERE
      x0=(idealflux-sum1glb+cumuldeltaflux)/sum0glb & ! Correction en une iteration
         *dte_fw/864000.  !10jours                    ! correction differee 
!        *0.01

      if(obcstatus(jeq1)==1) then !----->              
       sum1=0.
       do i=1,imax
        fluxbar_v(i,1,1)=fluxbar_v(i,1,1)+x0*mask_v(i,1,kmax)*h_v(i,1)
        sum1=sum1+mask_i_v(i)*mask_v(i,1,kmax)*fluxbar_v(i,1,1)
       enddo
      else                        !----->
       sum1=0.
      endif                       !----->
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      cumuldeltaflux=cumuldeltaflux+(idealflux-sum1glb)

! Evolution de la ssh moyenne dans le reservoir derriere la frontiere ouverte:
       ssh_reservoir=                         &
       ssh_reservoir+dte_fw*(                 &
                        signe*sum1glb         & ! Flux lateral en m3/s, >0 si entrant dans le reservoir
                           ! -surflux_        & ! Flux surface en m3/s, >0 si dirigE vers le haut c.a.d. sortant du reservoir
                            )/area              ! area = surperficie du reservoire

      if(par%rank==0) then
      if(iteration2d==1) then
       open(unit=3,file='tmp/flux_sicile_obc',position='append')
          write(3,'(4(1x,e14.7))') elapsedtime_now/86400. &
                                  ,sum1glb                & ! Flux instantannE
                                  ,cumuldeltaflux         & ! Flux cumulE
                                  ,ssh_reservoir            ! SSH du reservoir
       close(3)
      endif
      endif

      end subroutine obc_ext_sicileflux_i1

!.............................................................................

      subroutine obc_ext_nosplitting
      use module_principal
      use module_parallele !#MPI
      use module_q
      implicit none
      integer t_,loop_
      integer :: idi_u_x_,idi_u_u1_,idi_u_u2_
      integer :: idi_v_y_,idi_v_v1_,idi_v_v2_

      if(extmodtype.eq.'double')stop 'obc twin ext mode pas prevue'     !05/04/07

! Reference solution:
!     if(obcfreeorfix==0) then !rrrrrr>
!      call obc_ext_reference_standard
!     else                     !rrrrrr>
!      call obc_ext_reference_adaptative
!     endif                    !rrrrrr>



! Border conditions for tangential velocities:

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
       do i=1,imax+1 !2,imax
        fluxbar_u(i,jmax  ,1)= &
        fluxbar_u(i,jmax-1,1)  &
          *mask_u(i,jmax  ,kmax)
       enddo
      endif                              !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then    !----------->
       do i=1,imax+1 !2,imax
        fluxbar_u(i,1,1)=   &
        fluxbar_u(i,2,1)    &
          *mask_u(i,1,kmax) 
       enddo
      endif                              !----------->

! corners:
      if(obcstatus(ieqimax_jeqjmax)==1)fluxbar_u(imax+1,jmax,1)=fluxbar_u(imax,jmax,1)
      if(obcstatus(ieq1_jeqjmax)==1)   fluxbar_u(1     ,jmax,1)=fluxbar_u(2   ,jmax,1)
      if(obcstatus(ieqimax_jeq1)==1)   fluxbar_u(imax+1,1   ,1)=fluxbar_u(imax,1   ,1)
      if(obcstatus(ieq1_jeq1)==1)      fluxbar_u(1     ,1   ,1)=fluxbar_u(2   ,1   ,1)


! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
       do j=1,jmax+1 !2,jmax
        fluxbar_v(imax  ,j,1)=      &
        fluxbar_v(imax-1,j,1)       &
          *mask_v(imax  ,j,kmax)
       enddo
      endif                              !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->
       do j=1,jmax+1 !2,jmax
        fluxbar_v(1,j,1)=        &
        fluxbar_v(2,j,1)         &
          *mask_v(1,j,kmax) 
       enddo
      endif                           !----------->

! corners:
      if(obcstatus(ieqimax_jeqjmax)==1)fluxbar_v(imax,jmax+1,1)=fluxbar_v(imax,jmax  ,1)
      if(obcstatus(ieqimax_jeq1)==1)   fluxbar_v(imax,1     ,1)=fluxbar_v(imax,2     ,1)
      if(obcstatus(ieq1_jeqjmax)==1)   fluxbar_v(1   ,jmax+1,1)=fluxbar_v(1   ,jmax  ,1)
      if(obcstatus(ieq1_jeq1)==1)      fluxbar_v(1   ,1     ,1)=fluxbar_v(1   ,2     ,1)


! Border conditions for sea surface height:

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=1,imax !2,imax-1
                 ssh_int_w(i,jmax,2)=                              &
               sshrefobc_i(i,2)                                    &
            +(   fluxbar_v(i,jmax,1)                               &
                     -hz_v(i,jmax,1)*dx_v(i,jmax)*vbrrefobc_i(i,2) &
         )/(sqrt(grav*hz_v(i,jmax,1))                              &
                     *dx_v(i,jmax)   ) 

      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->

      do i=1,imax !2,imax-1
               ssh_int_w(i,1,2)=                           &
             sshrefobc_i(i,1)                              &
            -( fluxbar_v(i,2,1)                            &
                   -hz_v(i,1,1)*dx_v(i,1)*vbrrefobc_i(i,1) &
       )/(sqrt(grav*hz_v(i,1,1))                           &
                   *dx_v(i,1)   ) 

      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=1,jmax !2,jmax-1
                 ssh_int_w(imax,j,2)=                              &
                    sshrefobc_j(j,2)                               &
               +(fluxbar_u(imax,j,1)                               &
                     -hz_u(imax,j,1)*dy_u(imax,j)*vbrrefobc_j(j,2) &
!        )/(sqrt(grav*hz_u(imax,j,1))  &
         )/(c_                         & ! cas onde incidente
                     *dy_u(imax,j)   ) 
      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->

!     write(6,*)'Longueur',2.*pi/kvector_
!     write(6,*)'period_',period_
!     write(6,*)'c_',c_
!     write(6,*)'ssh_amplitude_',ssh_amplitude_
!     write(6,*)'cfl_reduce',cfl_reduce
!     write(6,*)'cfl_nh',cfl_nh
!     write(6,*)'nhpgf_reduce',nhpgf_reduce
!     stop 'coucou'

!#ifdef bidon
!..............................
! Cas onde incidente i croissant
      x0=min(1.,elapsedtime_now/(0.5*elapsedtime_end))*ssh_amplitude_

!     x1=sin(2.*pi*iteration3d*dti_fw/period_) &
!       *min(1.,iteration3d/(                  &
!                            1000.             &
!                            *0.5/cfl_reduce)) &  ! pour garder la meme rampe si on change le pas de temps
!       *min(1.,elapsedtime_now/period_)
!       *min(1.,elapsedtime_now/(0.5*elapsedtime_end))

      x1=x0*sin(2.*pi*iteration3d*dti_fw/period_)
      x2=x0*sin(2.*pi*iteration3d*dti_fw/period_-kvector_*0.5*dxb) ! ssh au point u

      do j=1,jmax 
! CAS AVEC UN ANGLE
! l'operation modulo sert pour assurer la continuite mpi:
!     x1=      &
!        -0.02 & ! dans le cas d'hamilton la frontiere est abaissEe de 2 cm
!    +x0*sin(2.*pi*iteration3d*dti_fw/period_-2.*pi*real(modulo(j+par%tjmax(1),jglb-2))/real(jglb-2))
!      sshrefobc_j(j,1)=ssh_amplitude_*x1

       sshrefobc_j(j,1)=x1
       vbrrefobc_j(j,1)=x2*c_/hz_u(1,j,1)

!      sshrefobc_j(j,1)*c_/(hz_u(1,j,1)          ) ! faute corrigee le 18/01/2018

      enddo

!#endif
!..............................

      do j=1,jmax !2,jmax-1

!#ifdef bidon
                 ssh_int_w(1,j,2)=                 &
                 sshrefobc_j(j,1)                  &
              -( fluxbar_u(2,j,1)                  &
                     -hz_u(1,j,1)*dy_u(1,j)*vbrrefobc_j(j,1)  &
!        )/(sqrt(grav*hz_u(1,j,1))    &
         )/(c_                        & ! cas onde incidente
                     *dy_u(1,j)   ) 
!#endif
!                 ssh_int_w(1,j,2)=sshrefobc_j(j,1)  ! c.l. clampee


      
!      if(j==jmax/2)write(64,*)real(elapsedtime_now)  &
!              ,real(ssh_int_w(1,j,2))       &
!              ,real(vbrrefobc_j(j,1))       &
!              ,real(fluxbar_u(2,j,1)/hz_u(1,j,1)/dy_u(1,j))

      enddo
      endif                                          !----------->

! Corners:
      if(obcstatus(ieq1_jeq1)==1)  then !----------->
       i1=1    ; j1=1 ; i2=1 ; j2=2 ; i3=2 ; j3=1
       ssh_int_w(i1,j1,2)=0.5*(ssh_int_w(i2,j2,2)+ssh_int_w(i3,j3,2))
      endif                                          !----------->

      if(obcstatus(ieq1_jeqjmax)==1) then !----------->
       i1=1    ; j1=jmax ; i2=1 ; j2=jmax-1 ; i3=2 ; j3=jmax
       ssh_int_w(i1,j1,2)=0.5*(ssh_int_w(i2,j2,2)+ssh_int_w(i3,j3,2))
      endif                                          !----------->

      if(obcstatus(ieqimax_jeq1)==1)    then !----------->
       i1=imax ; j1=1 ; i2=imax ; j2=2 ; i3=imax-1 ; j3=1
       ssh_int_w(i1,j1,2)=0.5*(ssh_int_w(i2,j2,2)+ssh_int_w(i3,j3,2))
      endif                                          !----------->

      if(obcstatus(ieqimax_jeqjmax)==1)   then !----------->
       i1=imax ; j1=jmax ; i2=imax ; j2=jmax-1 ; i3=imax-1 ; j3=jmax
       ssh_int_w(i1,j1,2)=0.5*(ssh_int_w(i2,j2,2)+ssh_int_w(i3,j3,2))
      endif                                          !----------->

! Border conditions for normal velocities by inversing the continuity
! equation and using ssh and tangential velocities:

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=1,imax !2,imax-1
      fluxbar_v(i,jmax+1,1)= (                                        &

      fluxbar_v(i,jmax  ,1)                                           &

       -( fluxbar_u(i+1,jmax,1)-fluxbar_u(i,jmax,1) +                 &

        ((ssh_int_w(i,jmax,2)-ssh_int_w(i,jmax,0))/dti_lp             &
!                              +omega_w(i,jmax,kmax+1,1)              &
                               )*dxdy_t(i,jmax))     &

                             )*mask_v(i,jmax+1,kmax)

      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->
      do i=1,imax !2,imax-1
      fluxbar_v(i,1,1) = (                             &

      fluxbar_v(i,2,1)                                 &

        +(fluxbar_u(i+1,1,1)-fluxbar_u(i,1,1) +        &

        ((ssh_int_w(i,1,2)-ssh_int_w(i,1,0))/dti_lp    &
!                           +omega_w(i,1,kmax+1,1)     &
                            )*dxdy_t(i,1))             &

                         )*mask_v(i,1,kmax)

      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=1,jmax !2,jmax-1

      fluxbar_u(imax+1,j,1) = (                              &

      fluxbar_u(imax  ,j,1)                                  &

        -(fluxbar_v(imax,j+1,1)-fluxbar_v(imax,j,1) +        &

        ((ssh_int_w(imax,j,2)-ssh_int_w(imax,j,0))/dti_lp    &
                               +omega_w(imax,j,kmax+1,1)     &
                               )*dxdy_t(imax,j))             &

                              )*mask_u(imax+1,j,kmax)

      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->
      do j=1,jmax !2,jmax-1
      fluxbar_u(1,j,1) = (                             &

      fluxbar_u(2,j,1)                                 &

        +(fluxbar_v(1,j+1,1)-fluxbar_v(1,j,1)+         &

        ((ssh_int_w(1,j,2)-ssh_int_w(1,j,0))/dti_lp    &
                            +omega_w(1,j,kmax+1,1)     &
                            )*dxdy_t(1,j))             &

                         )*mask_u(1,j,kmax)

      enddo
      endif                                          !----------->

! corners:
                                                                       !11/06/04
      if(obcstatus(ieqimax_jeqjmax)==1)   then !----------->
      i=imax
      j=jmax
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_int_w(i,j,2)- ssh_int_w(i,j,0))/dti_lp*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i+1,j  ,1))                                &
       +abs(fluxbar_v(i  ,j+1,1))                                &
                )

      fluxbar_u(i+1,j  ,1)=fluxbar_u(i+1,j  ,1)                  &
                   +x1*abs(fluxbar_u(i+1,j  ,1))                 &
                             *mask_u(i+1,j  ,kmax)
                                                                       !11/06/04
      fluxbar_v(i  ,j+1,1)=fluxbar_v(i  ,j+1,1)                  &
                   +x1*abs(fluxbar_v(i  ,j+1,1))                 &
                             *mask_v(i  ,j+1,kmax)

      endif                                          !----------->

      if(obcstatus(ieqimax_jeq1)==1)    then !----------->
      i=imax
      j=1
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_int_w(i,j,2)- ssh_int_w(i,j,0))/dti_lp*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i+1,j  ,1))                                &
       +abs(fluxbar_v(i  ,j  ,1))                                &
                )

      fluxbar_u(i+1,j  ,1)=fluxbar_u(i+1,j  ,1)                  &
                   +x1*abs(fluxbar_u(i+1,j  ,1))                 &
                             *mask_u(i+1,j  ,kmax)
                                                                       !11/06/04
      fluxbar_v(i  ,j  ,1)=fluxbar_v(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_v(i  ,j  ,1))                 &
                             *mask_v(i  ,j  ,kmax)

      endif                                          !----------->

      if(obcstatus(ieq1_jeqjmax)==1) then !----------->
      i=1
      j=jmax
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_int_w(i,j,2)- ssh_int_w(i,j,0))/dti_lp*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i  ,j  ,1))                                &
       +abs(fluxbar_v(i  ,j+1,1))                                &
                )

      fluxbar_u(i  ,j  ,1)=fluxbar_u(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_u(i  ,j  ,1))                 &
                             *mask_u(i  ,j  ,kmax)
                                                                       !11/06/04
      fluxbar_v(i  ,j+1,1)=fluxbar_v(i  ,j+1,1)                  &
                   +x1*abs(fluxbar_v(i  ,j+1,1))                 &
                             *mask_v(i  ,j+1,kmax)

      endif                                          !----------->

      if(obcstatus(ieq1_jeq1)==1)  then !----------->
      i=1
      j=1
      x1 = & ! rhs
        -(                                                       &
           (fluxbar_u(i+1,j  ,1)-   fluxbar_u(i,j,1))            &
          +(fluxbar_v(i  ,j+1,1)-   fluxbar_v(i,j,1))            &
          +( ssh_int_w(i,j,2)- ssh_int_w(i,j,0))/dti_lp*dxdy_t(i,j)      &
         ) / max(small1,                                         &
        abs(fluxbar_u(i  ,j  ,1))                                &
       +abs(fluxbar_v(i  ,j  ,1))                                &
                )

      fluxbar_u(i  ,j  ,1)=fluxbar_u(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_u(i  ,j  ,1))                 &
                             *mask_u(i  ,j  ,kmax)
                                                                      !11/06/04
      fluxbar_v(i  ,j  ,1)=fluxbar_v(i  ,j  ,1)                  &
                   -x1*abs(fluxbar_v(i  ,j  ,1))                 &
                             *mask_v(i  ,j  ,kmax)


      endif                                          !----------->


!                         / / /


! Ajuster les flux Normaux 3D A fluxbar:
! correction du type veldxdz_v=veldxdz_v+alpha*dz_v
      if(obcstatus(jeqjmax)==1)     then !----------->
! courant normal:
       k=1 ; j=jmax+1
       do i=1,imax
        xy_v(i,j,1)=veldxdz_v(i,j,k,1)
       enddo
       do k=2,kmax ; do i=1,imax
        xy_v(i,j,1)= &
        xy_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo       ; enddo
       do i=1,imax
        xy_v(i,j,1)=(fluxbar_v(i,j,1)-xy_v(i,j,1))/hz_v(i,j,1)
       enddo
       do k=1,kmax ; do i=1,imax
        veldxdz_v(i,j,k,1)= &
        veldxdz_v(i,j,k,1)+xy_v(i,j,1)*dz_v(i,j,k,1)
       enddo       ; enddo
! courant tangent
       k=1 ; j=jmax
       do i=1,imax+1
        xy_u(i,j,1)=veldydz_u(i,j,k,1)
       enddo
       do k=2,kmax ; do i=1,imax+1
        xy_u(i,j,1)= &
        xy_u(i,j,1)+veldydz_u(i,j,k,1)
       enddo       ; enddo
       do i=1,imax+1
        xy_u(i,j,1)=(fluxbar_u(i,j,1)-xy_u(i,j,1))/hz_u(i,j,1)
       enddo
       do k=1,kmax ; do i=1,imax+1
        veldydz_u(i,j,k,1)= &
        veldydz_u(i,j,k,1)+xy_u(i,j,1)*dz_u(i,j,k,1)
       enddo       ; enddo
      endif                              !----------->

      if(obcstatus(jeq1)==1)        then !----------->
! courant normal:
       k=1 ; j=1
       do i=1,imax
        xy_v(i,j,1)=veldxdz_v(i,j,k,1)
       enddo
       do k=2,kmax ; do i=1,imax
        xy_v(i,j,1)= &
        xy_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo       ; enddo
       do i=1,imax
        xy_v(i,j,1)=(fluxbar_v(i,j,1)-xy_v(i,j,1))/hz_v(i,j,1)
       enddo
       do k=1,kmax ; do i=1,imax
        veldxdz_v(i,j,k,1)= &
        veldxdz_v(i,j,k,1)+xy_v(i,j,1)*dz_v(i,j,k,1)
       enddo       ; enddo
! courant tangent
       k=1 ; j=1
       do i=1,imax+1
        xy_u(i,j,1)=veldydz_u(i,j,k,1)
       enddo
       do k=2,kmax ; do i=1,imax+1
        xy_u(i,j,1)= &
        xy_u(i,j,1)+veldydz_u(i,j,k,1)
       enddo       ; enddo
       do i=1,imax+1
        xy_u(i,j,1)=(fluxbar_u(i,j,1)-xy_u(i,j,1))/hz_u(i,j,1)
       enddo
       do k=1,kmax ; do i=1,imax+1
        veldydz_u(i,j,k,1)= &
        veldydz_u(i,j,k,1)+xy_u(i,j,1)*dz_u(i,j,k,1)
       enddo       ; enddo
      endif                              !----------->

      if(obcstatus(ieqimax)==1)     then !----------->
! courant normal
       k=1 ; i=imax+1
       do j=1,jmax
        xy_u(i,j,1)=veldydz_u(i,j,k,1)
       enddo
       do k=2,kmax ; do j=1,jmax
        xy_u(i,j,1)= &
        xy_u(i,j,1)+veldydz_u(i,j,k,1)
       enddo       ; enddo
       do j=1,jmax
        xy_u(i,j,1)=(fluxbar_u(i,j,1)-xy_u(i,j,1))/hz_u(i,j,1)
       enddo
       do k=1,kmax ; do j=1,jmax
        veldydz_u(i,j,k,1)= &
        veldydz_u(i,j,k,1)+xy_u(i,j,1)*dz_u(i,j,k,1)
       enddo       ; enddo
! courant tangent
       k=1 ; i=imax
       do j=1,jmax+1
        xy_v(i,j,1)=veldxdz_v(i,j,k,1)
       enddo
       do k=2,kmax ; do j=1,jmax+1
        xy_v(i,j,1)=xy_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo       ; enddo
       do j=1,jmax+1
        xy_v(i,j,1)=(fluxbar_v(i,j,1)-xy_v(i,j,1))/hz_v(i,j,1)
       enddo
       do k=1,kmax ; do j=1,jmax+1
        veldxdz_v(i,j,k,1)= &
        veldxdz_v(i,j,k,1)+xy_v(i,j,1)*dz_v(i,j,k,1)
       enddo       ; enddo
      endif                              !----------->

      if(obcstatus(ieq1   )==1)     then !----------->
! courant normal
       k=1 ; i=1
       do j=1,jmax
        xy_u(i,j,1)=veldydz_u(i,j,k,1)
       enddo
       do k=2,kmax ; do j=1,jmax
        xy_u(i,j,1)= &
        xy_u(i,j,1)+veldydz_u(i,j,k,1)
       enddo       ; enddo
       do j=1,jmax
        xy_u(i,j,1)=(fluxbar_u(i,j,1)-xy_u(i,j,1))/hz_u(i,j,1)
       enddo
       do k=1,kmax ; do j=1,jmax
        veldydz_u(i,j,k,1)= &
        veldydz_u(i,j,k,1)+xy_u(i,j,1)*dz_u(i,j,k,1)
       enddo       ; enddo
! courant tangent
       k=1 ; i=1
       do j=1,jmax+1 
        xy_v(i,j,1)=veldxdz_v(i,j,k,1)
       enddo
       do k=2,kmax ; do j=1,jmax+1
        xy_v(i,j,1)=xy_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo       ; enddo
       do j=1,jmax+1
        xy_v(i,j,1)=(fluxbar_v(i,j,1)-xy_v(i,j,1))/hz_v(i,j,1)
       enddo
       do k=1,kmax ; do j=1,jmax+1
        veldxdz_v(i,j,k,1)= &
        veldxdz_v(i,j,k,1)+xy_v(i,j,1)*dz_v(i,j,k,1)
       enddo       ; enddo

! Verif:
!      j=jmax/2
!      sum1=0.
!      do k=1,kmax
!       sum1=sum1+veldydz_u(i,j,k,1)
!       sum1=sum1+veldxdz_v(i,j,k,1)
!      enddo
!      write(557,*)sum1,fluxbar_v(i,j,1)
    

      endif                              !----------->

      end subroutine obc_ext_nosplitting

!.............................................................................
