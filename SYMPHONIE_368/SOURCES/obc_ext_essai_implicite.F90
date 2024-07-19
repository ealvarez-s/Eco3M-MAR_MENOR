      subroutine obc_ext(txt_)
!______________________________________________________________________
! S model
! release S26 - last update: 12-08-15
!______________________________________________________________________
      use module_principal
      implicit none
      character txt_*5
#ifdef synopsis
       subroutinetitle='obc_ext'
       subroutinedescription=                                &
        'driver routine for lateral boundary conditions for' &
       ' for barotropic fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ Compute the (physical) Open Boundary Conditions of the external mode momentum equations
!$ according to:  Marsaleix P., Auclair F., Estournel C., 2006, Considerations on Open
!$ Boundary Conditions for Regional and Coastal Ocean Models. Journal of Atmospheric and Oceanic
!$ Technology, 23,1604-1613, http://dx.doi.org/10.1175/JTECH1930.1

!...............................................................................
! Version Date      Description des modifications
!         25/11/01: introduction du cas riverdir=0 (fleuve introduit à la
!                   surface et non plus lateralement)
!         26/08/02: VBAROBC remplacé par VMEAOBC
!         11/06/04: procedure des coins revue. On calcule d'abord une ebauche
!                   des coins au moment de l'etape des composantes tgte à la
!                   frontiere. Ensuite un coef de lagrange est utilisé pour
!                   corriger le courant des coins afin d'etre compatible avec
!                   la sse. La perturbation du courant est proportionnelle au
!                   courant lui-même.
!         24/06/04: projet (les lignes sont en commentaire) de faire
!                   une correction du transport (correction bilan volume)
!                   qui soit proportionnelle à l'erreur sur le courant
!                   autrement dit l'ecart du transport à H*VMEAOBC.
!         12/07/04: Si relax_es < 0 pas de correction de la moyenne
!                   de l'elevation de la surface
!         03/08/05: Cas "Flather forcage libre"
!         10/11/05: debug du point precedent + extension de la methode
!                   aux composantes du transport tangeantes aux OBC
!         12/01/06: les variables termes contenant obc sont regroupés
!                   et calculés une fois par iteration barocline
!                   La sse basse frequence voit sa moyenne controlée par
!                   ztaobc
!         20/01/06: amenagements sur conservation du toit
!         09/03/06: condition vitesses tangentes remis à la version originale
!                   car bug sinon. Appel à updateforcing rationnalise.
!         31/05/06: Correction d'un bug dans la condition au limite de T et S
!                   du au fait que les coins de la solution
!                   filtree (basse frequence) n'avait pas ete codes...
!         05/04/07: La condition aux limites integre maintenant la condition
!                   adaptative. Une version generalisee est employee, c'est à
!                   dire que la même condition de Flather est utilisee, il n'y
!                   a que la fonction de construire la condition de reference
!                   qui change.
!                   D'autre part, dans la cas standard, on a prevu une transition
!                   lineaire du forcage pour eviter les discontinuités liées
!                   au time spliting, discontinuités qui se propageaient dans le
!                   domaine sous forme d'ondes de gravités.
!          17/04/07: Passage à coordonnees curvilignes (voir ajout dx_y et cie...)
!          13/12/07: Le point de variablilité lente pour la vitesse est deplacé
!                    au point voisin.
!          12-05-09  le suivi de la moyenne de la ssh est introduit dans la filiere
!                    nesting
!          16-05-09  - Ajout du suivi de la moyenne de la ssh de la maree
!                    - Amortisseur de la correction moyenne revisité
!          04-06-09  Parallelisation
! 2009.3   02-10-09  dte_lp remplace dte
!          05-10-09  ajout d'un "ifdef parallele"
!          09-10-09  parallelisation de la conservation du toit
! 2010.7   22-02-10  velbar devient la variable d'etat
!          01-03-10  correction du bug introduit avec le point precedent
! 2010.8   22-03-10  hssh_w remplace ssh_ext_w+h_w
!          24-03-10  La condition riviere sur fluxbar et la CL ouverte sur velbar
!                    sont disssossiées car pas appliquees aux memes moments
! 2010.10  13-06-10  suppression ssh_ext_w
! 2010.13  01-11-10  river_inout remplace river_dom
! 2010.20  19-04-11  k2dite,k2dfin renommés iteration2d,iteration2d_max_now
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
!          30-07-15  ajout obc_ext_mpi_xyt(t_,txt_) pour echanger xy_t
!          12-08-15  condition de flather par methode implicite
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
      implicit none
      integer t_,loop_
      integer :: idi_u_x_,idi_u_u1_,idi_u_u2_
      integer :: idi_v_y_,idi_v_v1_,idi_v_v2_
#ifdef synopsis
       subroutinetitle='obc_ext_fb'
       subroutinedescription= &
       'Lateral boundary conditions for barotropic fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
! Mes notes sont dans:
! https://docs.google.com/document/d/1bjpXJW7aLWJm_KqzvArKtg7hglD19mjUjC6cVchA9IM/edit

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=2,imax-1

      velbar_v(i,jmax+1,2)= &
        mask_v(i,jmax+1,kmax)*( & !oooo>

      vbrrefobc_i(i,2)-(1./sqr_hoverg_v(i,2))*( & !pppp>

      sshrefobc_i(i,2)-ssh_w(i,jmax,1)+ dte_fw*omega_w(i,jmax,kmax+1,1)&
                                      +(dte_fw/ dxdy_t(i,jmax))*( & !mmm>

!  fluxbar_v(i,jmax+1) except (since implicited) v_euler
       (                   velbarstokes_v(i,jmax+1,1))                 &
           *dy_v(i,jmax+1)                                             &
           *(h_v(i,jmax+1)+0.5*(ssh_w(i,jmax+1,1)+ssh_w(i,jmax,1)))    &

! -fluxbar_v(i,jmax):
      -(velbar_v(i,jmax,2)+velbarstokes_v(i,jmax  ,1))                 &
           *dy_v(i,jmax)                                               &
           *(h_v(i,jmax)+0.5*(ssh_w(i,jmax,1)+ssh_w(i,jmax-1,1)))      &

          +fluxbar_u(i+1,jmax,1)-fluxbar_u(i,jmax,1)                   &

                                                         ) & !mmm>

                                            ) & !pppp>

                        ) & !oooo>
      /(1.+(h_v(i,jmax+1)+0.5*(ssh_w(i,jmax+1,1)+ssh_w(i,jmax,1))) &
          *dy_v(i,jmax+1)                                  &
          *dte_fw/(dxdy_t(i,jmax)*sqr_hoverg_v(i,2)))

      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->
      do i=2,imax-1

      velbar_v(i,1,2)= &
        mask_v(i,1,kmax)*( & !oooo>

      vbrrefobc_i(i,1)+(1./sqr_hoverg_v(i,1))*( & !pppp>

      sshrefobc_i(i,1)-ssh_w(i,1,1)+ dte_fw*omega_w(i,1,kmax+1,1)&
                                   +(dte_fw/ dxdy_t(i,1))*( & !mmm>

!  fluxbar_v(i,2)
      (velbar_v(i,2,2)+velbarstokes_v(i,2,1))                       &
          *dy_v(i,2)                                                &
          *(h_v(i,2)+0.5*(ssh_w(i,2,1)+ssh_w(i,1,1)))               &

! -fluxbar_v(i,1) except (since implicited) v_euler
      -(                velbarstokes_v(i,1,1))                      &
          *dy_v(i,1)                                                &
          *(h_v(i,1)+0.5*(ssh_w(i,1,1)+ssh_w(i,0,1)))               &

          +fluxbar_u(i+1,1,1)-fluxbar_u(i,1,1)                      &

                                                         ) & !mmm>

                                            ) & !pppp>

                        ) & !oooo>
      /(1.+(h_v(i,1)+0.5*(ssh_w(i,1,1)+ssh_w(i,0,1))) &
          *dy_v(i,1)                                  &
          *dte_fw/(dxdy_t(i,1)*sqr_hoverg_v(i,1)))

      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=2,jmax-1 ! 1,jmax ! 2,jmax-1

       velbar_u(imax+1,j,2)=  &
         mask_u(imax+1,j,kmax)*( & !oooo>

       vbrrefobc_j(j,2)-(1./sqr_hoverg_u(j,2))*( & !pppp>

       sshrefobc_j(j,2)-ssh_w(imax,j,1)+dte_fw*omega_w(imax,j,kmax+1,1)&
                                      +(dte_fw/ dxdy_t(imax,j))*( & !mmm>

!  fluxbar_u(imax+1,j) except (since implicited) u_euler
       (                   velbarstokes_u(imax+1,j,1))                 &
           *dy_u(imax+1,j)                                             &
           *(h_u(imax+1,j)+0.5*(ssh_w(imax+1,j,1)+ssh_w(imax,j,1)))    &


! -fluxbar_u(imax  ,j)
      -(velbar_u(imax,j,2)+velbarstokes_u(imax,j,1))                   &
           *dy_u(imax,j)                                               &
           *(h_u(imax,j)+0.5*(ssh_w(imax,j,1)+ssh_w(imax-1,j,1)))      &

          +fluxbar_v(imax,j+1,1)-fluxbar_v(imax,j,1)                   &

                                                                ) & !mmm>
                                               ) & !pppp>
                               ) & !oooo>
      /(1.+(h_u(imax+1,j)+0.5*(ssh_w(imax+1,j,1)+ssh_w(imax,j,1))) &
          *dy_u(imax+1,j)                                          &
          *dte_fw/(dxdy_t(imax,j)*sqr_hoverg_u(j,2)))

      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->
      do j=2,jmax-1 ! 1,jmax !2,jmax-1

       velbar_u(1,j,2)= &
         mask_u(1,j,kmax)*( & !oooo>

      vbrrefobc_j(j,1)+(1./sqr_hoverg_u(j,1))*( & !pppp>

      sshrefobc_j(j,1)-ssh_w(1,j,1)+ dte_fw*omega_w(1,j,kmax+1,1)&
                                   +(dte_fw/ dxdy_t(1,j))*( & !mmm>

!  fluxbar_u(2,j)
      (velbar_u(2,j,2)+velbarstokes_u(2,j,1))                       &
          *dy_u(2,j)                                                &
          *(h_u(2,j)+0.5*(ssh_w(2,j,1)+ssh_w(1,j,1)))               &

! -fluxbar_u(1,j) except (since implicited) u_euler
      -(                velbarstokes_u(1,j,1))                      &
          *dy_u(1,j)                                                &
          *(h_u(1,j)+0.5*(ssh_w(1,j,1)+ssh_w(0,j,1)))               &

          +fluxbar_v(1,j+1,1)-fluxbar_v(1,j,1)                      &

                                                         ) & !mmm>

                                            ) & !pppp>

                        ) & !oooo>
      /(1.+(h_u(1,j)+0.5*(ssh_w(1,j,1)+ssh_w(0,j,1))) &
          *dy_u(1,j)                                  &
          *dte_fw/(dxdy_t(1,j)*sqr_hoverg_u(j,1)))

      enddo
      endif                                          !----------->

! Corners step1:
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

! corners step2:
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


      if(relax_es>=0.)call obc_ext_sshmean_control

!                         / / /

! tangential depth-averaged current

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       i=1
        do j=1,jmax+1
             velbar_v(i,j,2)=                                 &
               mask_v(i,j,kmax)*( & !ooo>

            fluxbar_v(i,j,1)                                  &
                /(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))) &
                /dx_v(i,j)                                    &

      -velbarstokes_v(i,j,1)                                  &

                                )   !ooo>
        enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       j=1
       do i=1,imax+1
             velbar_u(i,j,2)=                                 &
               mask_u(i,j,kmax)*( & !pmx>

            fluxbar_u(i,j,1)                                  &
                /(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))) &
                /dy_u(i,j)                                    &

      -velbarstokes_u(i,j,1)                                  &
                                )   !pmx>
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       i=imax
        do j=1,jmax+1
             velbar_v(i,j,2)=                                 &
               mask_v(i,j,kmax)*( & !pmx>

            fluxbar_v(i,j,1)                                  &
                /(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))) &
                /dx_v(i,j)                                    &

      -velbarstokes_v(i,j,1)                                  &
                                )   !pmx>
        enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       j=jmax
       do i=1,imax+1
             velbar_u(i,j,2)=                                 &
               mask_u(i,j,kmax)*( & !pmx>

            fluxbar_u(i,j,1)                                  &
                /(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))) &
                /dy_u(i,j)                                    &

      -velbarstokes_u(i,j,1)                                  &
                                )   !pmx>
       enddo
      endif                                      !----------->

! La condition aux limites sur les rivieres vient obligatoirement apres
! la conversion des flux en vitesses sinon les rivieres aux bords du domaine
! sont ecrasés par la conversion...

! Lignes et colonnes supplementaires pour advection o4:
      call obc_ext_o4

! Echange mpi u1 u2 v1 v2 sur velbar_u(:,:,2) et velbar_v(:,:,2)
      call obc_ext_mpi(2) !02-07-14

      end subroutine obc_ext_fb

!..............................................................

      subroutine obc_ext_o4
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='obc_ext_o4'
       subroutinedescription= &
       '"+2 range" lateral boundary conditions on barotropic velocity' &
       //' required by high-order advection schemes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
#ifdef synopsis
       subroutinetitle='obc_ext_mpi3d'
       subroutinedescription=&
       'Lateral boundary conditions ensuring mpi continuity of the' &
       //' barotropic velocity field (velbar_u velbar_v)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
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
#endif

      end subroutine obc_ext_mpi3d

!..............................................................

      subroutine obc_ext_mpi(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_u1_,idi_u2_,idi_v1_,idi_v2_
#ifdef synopsis
       subroutinetitle='obc_ext_mpi'
       subroutinedescription= &
       'mpi continuity of the barotropic velocity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


#ifdef parallele
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
#endif


      end subroutine obc_ext_mpi

!..............................................................

      subroutine obc_ext_mpi_xy(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_u1_,idi_u2_,idi_v1_,idi_v2_
#ifdef synopsis
       subroutinetitle='obc_ext_mpi_xy'
       subroutinedescription= &
       'mpi continuity of the temporary arrays xy_u and xy_v'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


#ifdef parallele
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
#endif

      end subroutine obc_ext_mpi_xy

!..............................................................

      subroutine obc_ext_velbar_mpi(t_) !15-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_x_,idi_y_
#ifdef synopsis
       subroutinetitle='obc_ext_velbar_mpi'
       subroutinedescription= &
       'mpi continuity of the barotropic velocity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


#ifdef parallele
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
#endif


      end subroutine obc_ext_velbar_mpi

!..............................................................

      subroutine obc_ext_river
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idi_x_,idi_y_
#ifdef synopsis
       subroutinetitle='obc_ext_river'
       subroutinedescription= &
       'river barotropic boundary condition'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ Begin the Boundary conditions at the river grid nodes:

      if(flag_externalmode_lf_or_fb=='lf')t_=0
      if(flag_externalmode_lf_or_fb=='fb')t_=1 !01-05-13

      do 100 kr=1,nriver

      if(riverdir(kr).ne.0)then!>>>>>>>>>>>>>>>>>>>>>>>>>>>            !25/11/01

      if(rivertrc_inout(kr)==1)ssh_w(iriver(kr,1),jriver(kr,1),t_)= & !22-01-15
                 max(wetdry_cst2-h_w(iriver(kr,1),jriver(kr,1)),zero)

      if(rivervel_inout(kr)==1) then !rrrrrrrrrrrrr>

      if(mod(riverdir(kr),2).eq.0) then ! ----------------->
         fluxbar_v(iriver(kr,2),jriver(kr,2),1)=riverflux(kr,1)
! Pour le moment on assume velbar=0
!         velbar_v(iriver(kr,2),jriver(kr,2),1)=riverflux(kr,1)        &
!      /       h_v(iriver(kr,2),jriver(kr,2))                          &
!            /dx_v(iriver(kr,2),jriver(kr,2))

      else
         fluxbar_u(iriver(kr,2),jriver(kr,2),1)=riverflux(kr,1)
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
#ifdef synopsis
       subroutinetitle='obc_ext_reference_standard'
       subroutinedescription=                                   &
       'computes sshrefobc vbrrefobc the reference solution of' &
       'the radiative barotropic open boundary condition'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
#ifdef synopsis
       subroutinetitle='obc_ext_reference_adaptative'
       subroutinedescription=                                   &
       'computes sshrefobc vbrrefobc the reference solution of' &
       'the radiative barotropic open boundary condition'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
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
!     integer t_,loop_,idi_x_,idi_y_
#ifdef synopsis
       subroutinetitle='obc_ext_sshmean_control'
       subroutinedescription=                                   &
        'SSH domain-average controled by adaptation of border'  &
       ' barotropic lateral fluxes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!_______________________________________________________________________________
! OPTION CONSERVATION DE LA HAUTEUR D'EAU
! Pour eviter que le bassin ne se vide ou se ne remplisse
! artificiellement, on controle le bilan de l'intégrale du
! transport normal à la frontière ouverte. Les composantes
! du transport à la frontiere sont perturbées pour que le
! bilan de l'intégrale soit équilibré
! Pour que cette option soit active il faut que RELAX_ES soit
! non nul.
! RELAX_ES est initialisé dans notebook_time.
! DEBUT:
!_______________________________________________________________________________
      pause 'la conservation de la ssh estparallelisee maispas verifiee'

      sum0=0.
      sum2=0.
! Les tests sur mpi_proc_null verifie que la frontiere n'est pas ouverte
! sur un autre sous-domaine:
      if(obcstatus(jeqjmax)==1) then !----->              !09-10-09
      do i=1,imax
       sum0=sum0+mask_v(i,jmax+1,kmax)*      h_v(i,jmax+1  )              &
              *mask_i_v(i)                                              !09-10-09
       sum2=sum2+mask_v(i,jmax+1,kmax)*xy_v(i,jmax+1,1)              &
              *mask_i_v(i)                                              !09-10-09

      enddo
      endif                                        !----->

      if(obcstatus(jeq1)==1) then !----->              !09-10-09
      do i=1,imax
       sum0=sum0+mask_v(i,1     ,kmax)*      h_v(i,1       )              &
              *mask_i_v(i)                                              !09-10-09
       sum2=sum2-mask_v(i,1     ,kmax)*xy_v(i,1     ,1)              &
              *mask_i_v(i)                                              !09-10-09
      enddo
      endif                                        !----->

      if(obcstatus(ieqimax)==1) then !----->              !09-10-09
      do j=1,jmax
       sum0=sum0+mask_u(imax+1,j,kmax)*      h_u(imax+1,j  )              &
              *mask_j_u(       j)                                       !09-10-09
       sum2=sum2+mask_u(imax+1,j,kmax)*xy_u(imax+1,j,1)              &
              *mask_j_u(       j)                                       !09-10-09
      enddo
      endif                                        !----->

      if(obcstatus(ieq1)==1) then !----->              !09-10-09
      do j=1,jmax
       sum0=sum0+mask_u(1     ,j,kmax)*      h_u(1     ,j  )              &
              *mask_j_u(       j)                                       !09-10-09
       sum2=sum2-mask_u(1     ,j,kmax)*xy_u(1     ,j,1)              &
              *mask_j_u(       j)                                       !09-10-09
      enddo
      endif                                        !----->

      if(nriver.ge.1) then !------------------------>
       do kr=1,nriver
       if(riverdir(kr).ne.0)sum2=sum2-abs(riverflux(kr,1))             !25/11/01
       enddo
      endif                !------------------------>

!     SOM2=SOM2+IOBC_AF*ZTAOBC_CUM(1) ! pour suivre les variations de sse de l'ogcm
      sum2=sum2+        sshobc_cum(1) & !                                !12-05-09
                       +sshtide_cum                                    !16-05-09

#ifdef parallele
      call mpi_allreduce(b2delev,b2delev_glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,         & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum0=sum0glb ;  sum2=sum2glb
#else
      b2delev_glb=b2delev ! attention à ne pas faire dans l'autre sens, c.a.d. b2delev=b2delev_glb
#endif

      lagrange_ssh=                                                     &
         b2delev_glb**2/(b2delev_glb**2+(100.*sum2)**2+small2)*     & ! Amortisseur !16-05-09
             (1.-relax_es)*lagrange_ssh                             & ! Accelerateur
            +    relax_es *(-b2delev_glb-sum2)/max(sum0,small2)     ! rappel
! L'amortisseur permet d'annuler l'accelerateur quand l'equilibre est retrouvé,
! evitant ainsi de nourrir une resonnance à la frequence du rappel.
! On considere que l'equilire est retrouvé quand B2DELEV est petit par rapport à SOM2
! , l'amortisseur tendant alors vers zero. Inversement un grand desequilbre est
! caracterisé par B2DELEV beaucoup plus grand que SOM2, l'amortisseur tendant vers 1.

! La ligne c1 permet d'accelerer la correction mais il est instable
! ce qui explique la presence d'un amortisseur (100*relax_es).
! Commenté par default. A utiliser avec precaution. A ne pas
! utiliser quand le protocole de correction n'est pas continue
! comme c'est le cas par exemple quand la correction depend de
! l'ecart au courant de forcage. On peut l'employer quand le protocole
! de correction est fixe comme par exemple proportionnel à la bathy car
! dans ce cas la "memoire" de lagrange_zta a un sens...

      if(obcstatus(jeqjmax)==1) then !----->              !09-10-09
       do i=1,imax
                     xy_v(i,jmax+1,1)=                             &
                    (xy_v(i,jmax+1,1)                              &
             +lagrange_ssh*h_v(i,jmax+1  ))                             &
                       *mask_v(i,jmax+1,kmax)

        b2delev=b2delev+mask_v(i,jmax+1,kmax)                             &
                    *xy_v(i,jmax+1,1)                              &
                     *mask_i_v(i)                                 !09-10-09
       enddo
      endif                                        !----->

      if(obcstatus(jeq1)==1) then !----->              !09-10-09
       do i=1,imax
                     xy_v(i,1,1)=                                  &
                    (xy_v(i,1,1)                                   &
             -lagrange_ssh*h_v(i,1  ))                                  &
                       *mask_v(i,1,kmax)

        b2delev=b2delev-mask_v(i,1,kmax)                                  &
                    *xy_v(i,1,1)                                   &
                     *mask_i_v(i)                                 !09-10-09
       enddo
      endif                                        !----->

      if(obcstatus(ieqimax)==1) then !----->              !09-10-09
       do j=1,jmax
                    xy_u(imax+1,j,1)=                              &
                   (xy_u(imax+1,j,1)                               &
            +lagrange_ssh*h_u(imax+1,j  ))                              &
                      *mask_u(imax+1,j,kmax)

       b2delev=b2delev+mask_u(imax+1,j,kmax)                              &
                   *xy_u(imax+1,j,1)                               &
                    *mask_j_u(       j)                                 !09-10-09
       enddo
      endif                                        !----->

      if(obcstatus(ieq1)==1) then !----->              !09-10-09
       do j=1,jmax
                     xy_u(1,j,1)=                                  &
                    (xy_u(1,j,1)                                   &
             -lagrange_ssh*h_u(1,j  ))                                  &
                       *mask_u(1,j,kmax)

        b2delev=b2delev-mask_u(1,j,kmax)                                  &
                    *xy_u(1,j,1)                                   &
                     *mask_j_u(  j)                                 !09-10-09
       enddo
      endif                                        !----->


      if(nriver.ge.1) then !------------------------>
       do kr=1,nriver
       if(riverdir(kr).ne.0)b2delev=b2delev-abs(riverflux(kr,1))       !25/11/01
       enddo
      endif                !------------------------>

!     B2DELEV=B2DELEV+IOBC_AF*ZTAOBC_CUM(1) ! pour suivre les variations de sse de l'ogcm
      b2delev=b2delev+sshobc_cum(1)                            & !12-05-09
                    +sshtide_cum                               !16-05-09

      end subroutine obc_ext_sshmean_control
!..................................................................................
#ifdef egrid
      subroutine obc_ext_egrid_fla
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='obc_ext_egrid_fla'
       subroutinedescription= &
        'Border conditions for barotropic fluxes' 
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Border conditions for barotropic fluxes

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)                 then !----------->
       j=jmax
       do i=1,imax 
! Composante 1:
        uflux2d_f(i+1,j+1,1)=      &
             0.5*mask_t(i,j,kmax)  & ! 50% partage avec la seconde composante du courant
           *sqrt(dxdy_t(i,j)*grav*(h_w(i,j)+ssh_w(i,j,1)))*ssh_w(i,j,0)
! Composante 2:
        vflux2d_f(i,j+1,1)=uflux2d_f(i+1,j+1,1)  
       enddo
      endif                                          !----------->
! j=1 south border:
      if(obcstatus(jeq1)==1)                    then !----------->
       j=1
       do i=1,imax 
! Composante 1:
        uflux2d_f(i,j,1)=      &
            -0.5*mask_t(i,j,kmax)  & ! 50% partage avec la seconde composante du courant
           *sqrt(dxdy_t(i,j)*grav*(h_w(i,j)+ssh_w(i,j,1)))*ssh_w(i,j,0)
! Composante 2:
        vflux2d_f(i+1,j,1)=uflux2d_f(i,j,1)  
       enddo
      endif                                          !----------->
! i=imax east border:
      if(obcstatus(ieqimax)==1)                 then !----------->
       i=imax
       do j=1,jmax ! 2,jmax-1 
! Composante 1:
        uflux2d_f(i+1,j+1,1)=      &
             0.5*mask_t(i,j,kmax)  & ! 50% partage avec la seconde composante du courant
           *sqrt(dxdy_t(i,j)*grav*(h_w(i,j)+ssh_w(i,j,1)))*ssh_w(i,j,0)
! Composante 2:
        vflux2d_f(i+1,j,1)=-uflux2d_f(i+1,j+1,1)  
       enddo
      endif                                          !----------->
! i=1 west border:
      if(obcstatus(ieq1)==1)                    then !----------->
       i=1 
       do j=1,jmax ! 2,jmax-1
! Composante 1:
        uflux2d_f(i,j,1)=      &
            -0.5*mask_t(i,j,kmax)  & ! 50% partage avec la seconde composante du courant
           *sqrt(dxdy_t(i,j)*grav*(h_w(i,j)+ssh_w(i,j,1)))*ssh_w(i,j,0)
! Composante 2:
        vflux2d_f(i,j+1,1)=-uflux2d_f(i,j,1)  
       enddo
      endif                                          !----------->
      end subroutine obc_ext_egrid_fla

!..........................................................................

      subroutine obc_ext_egrid_mpi(t_) !08-02-15
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_,idu1_,idv1_,idu2_,idv2_ &
             ,idexr_,ideyr_,idexr1_,ideyr1_
#ifdef synopsis
       subroutinetitle='obc_ext_egrid_mpi'
       subroutinedescription= &
        'mpi continuity of the barotropic velocity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
      write(texte30,'(a10,i0)')'flux2d_u1_',t_
      call get_type_echange('u1',trim(texte30),flux2d_u      &
                                       ,lbound(flux2d_u)     &
                                       ,ubound(flux2d_u)     &
                                       ,t_                   &
                                       ,idu1_)

      write(texte30,'(a10,i0)')'flux2d_v1_',t_
      call get_type_echange('v1',trim(texte30),flux2d_v      &
                                       ,lbound(flux2d_v)     &
                                       ,ubound(flux2d_v)     &
                                       ,t_                   &
                                       ,idv1_)

      write(texte30,'(a10,i0)')'flux2d_u2_',t_
      call get_type_echange('u2',trim(texte30),flux2d_u      &
                                       ,lbound(flux2d_u)     &
                                       ,ubound(flux2d_u)     &
                                       ,t_                   &
                                       ,idu2_)

      write(texte30,'(a10,i0)')'flux2d_v2_',t_
      call get_type_echange('v2',trim(texte30),flux2d_v      &
                                       ,lbound(flux2d_v)     &
                                       ,ubound(flux2d_v)     &
                                       ,t_                   &
                                       ,idv2_)


!!$! Nouvelle methode avec choix des voisins
      write(texte30,'(a11,i0)')'vlxbar_f_r_',t_
      call get_type_echange('r ',trim(texte30),vlxbar_f      &
                                       ,lbound(vlxbar_f)     &
                                       ,ubound(vlxbar_f)     &
                                       ,t_                   &
                                       ,idexr_)

      write(texte30,'(a11,i0)')'vlybar_f_r_',t_
      call get_type_echange('r ',trim(texte30),vlybar_f      &
                                       ,lbound(vlybar_f)     &
                                       ,ubound(vlybar_f)     &
                                       ,t_                   &
                                       ,ideyr_)

!     write(texte30,'(a12,i0)')'vlxbar_f_r1_',t_              !02-04-15
!     call get_type_echange('r1',trim(texte30),vlxbar_f      &
!                                      ,lbound(vlxbar_f)     &
!                                      ,ubound(vlxbar_f)     &
!                                      ,t_                   &
!                                      ,idexr1_)

!     write(texte30,'(a12,i0)')'vlybar_f_r1_',t_
!     call get_type_echange('r1',trim(texte30),vlybar_f      &
!                                      ,lbound(vlybar_f)     &
!                                      ,ubound(vlybar_f)     &
!                                      ,t_                   &
!                                      ,ideyr1_)

      do loop_=1,subcycle_exchange
         call echange_voisin(vlxbar_f,idexr_,mpi_neighbor_list(loop_))
         call echange_voisin(vlybar_f,ideyr_,mpi_neighbor_list(loop_))
!        call echange_voisin(vlxbar_f,idexr1_,mpi_neighbor_list(loop_))
!        call echange_voisin(vlybar_f,ideyr1_,mpi_neighbor_list(loop_))
         call echange_voisin(flux2d_u,idu1_,mpi_neighbor_list(loop_))
         call echange_voisin(flux2d_v,idv1_,mpi_neighbor_list(loop_))
         call echange_voisin(flux2d_u,idu2_,mpi_neighbor_list(loop_))
         call echange_voisin(flux2d_v,idv2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine obc_ext_egrid_mpi
#endif

!..............................................................

      subroutine obc_ext_egrid_xyf(exchtype_,indx1_,indx2_) !08-02-15
      use module_principal ; use module_parallele
      implicit none
      integer loop_,indx1_,indx2_,exchid1_,exchid2_
      character*2 exchtype_ ! 'r ' ou 'r1'
#ifdef synopsis
       subroutinetitle='obc_ext_egrid_xyf'
       subroutinedescription= &
        'mpi continuity of the barotropic velocity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
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
#endif

      end subroutine obc_ext_egrid_xyf
!..............................................................

      subroutine obc_ext_mpi_xyt(t_,txt_) !30-07-15
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_
      character*2 txt_
#ifdef synopsis
       subroutinetitle='obc_ext_mpi_xyt'
       subroutinedescription= &
       'mpi continuity of the temporary arrays xy_t'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


#ifdef parallele
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
#endif

      end subroutine obc_ext_mpi_xyt

!..............................................................

      subroutine obc_ext_qdmflux
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='obc_ext_qdmflux'
       subroutinedescription= &
        'Border condition for barotropic QDM fluxes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! j=jmax north border:
      if(obcstatus(jeqjmax)==1)     then !----------->
      do i=2,imax-1
       yflux_t(i,jmax)=yflux_t(i,jmax-1)
      enddo
      endif                                          !----------->

! j=1 south border:
      if(obcstatus(jeq1)==1)     then !----------->
      do i=2,imax-1
       yflux_t(i,1)=yflux_t(i,2)
      enddo
      endif                                          !----------->

! i=imax east border:
      if(obcstatus(ieqimax)==1)     then !----------->
      do j=2,jmax-1
       xflux_t(imax,j)=xflux_t(imax-1,j)
      enddo
      endif                                          !----------->

! i=1 west border:
      if(obcstatus(ieq1)==1)     then !----------->
      do j=2,jmax-1
       xflux_t(1,j)=xflux_t(2,j)
      enddo
      endif                                          !----------->


      end subroutine obc_ext_qdmflux
