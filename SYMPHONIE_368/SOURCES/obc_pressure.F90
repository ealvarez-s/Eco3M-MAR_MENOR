      subroutine obc_pressure(ichoix)
!______________________________________________________________________
!
! S model
! release 2010.10  - last update: 24-06-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________
      use module_principal
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='obc_pressure'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
!Version date      Description des modifications
!        24/02/08: mise en service
! 2009.3 01-10-09: utilisation des nouveaux facteurs d'echelle verticale
!...............................................................................

      stop 'obc_pressure non activé'

!----------------------------------------------------------------------
! Methode basique avec reference donnée par la basse frequence du model_
! Debut:
      if(ichoix.eq.1) then
!----------------------------------------------------------------------


! Phase speed:
      const5=1. ! phase speed (m/s)

! Coef equation d'etat linearisee:
      const1=-rho*alp_t
      const2= rho*alp_s
      const3=-const1*t0-const2*s0

      const8=min( dti_fw/86400. , un)
!     CONST8=0. ! reference = etat initial
!     CONST8=1.

!     WRITE(6,*)CONST8
!     STOP 'coco'
      if(kount.eq.kount0)const8=1.
      const7=1.-const8

!     STOP 'Ne passe Pas par obc_pressure'

      do loop1=1,2 ! debut boucle loop1

       if(loop1.eq.1)then !>>>>>
        i1=1
        i2=2
        i3=i1+1
        i4=i2+1
        i5=i1+2
        j1=1
        j2=2
        j3=j1+1
        j4=j2+1
        j5=j1+2
        const4=-rho*const5
       else
        i1=imax
        i2=imax
        i3=i1-1
        i4=i2-1
        i5=i1-2
        j1=jmax
        j2=jmax
        j3=j1-1
        j4=j2-1
        j5=j1-2
        const4=+rho*const5
       endif              !>>>>>

       do j=2,jmax-1
        xy_t(i1,j,1)=0.
        xy_t(i1,j,2)=0.
       enddo
       do i=2,imax-1
        xy_t(i,j1,1)=0.
        xy_t(i,j1,2)=0.
       enddo

      do k=kmax,1,-1

       do j=2,jmax-1

! L'anomalie de pression avant correction barotrope:
           anyv3d(i1,j,k,1)=                                            &
            veldydz_u(i2,j,k,0)                                         &
                /dz_u(i2,j,k,0)                                         &
                /dy_u(i2,j,0)                                          & !01-10-09
!    &    /dsig_x(i2,j,k)
!    &    /hzdy_x(i2,j,0)
          *const4

! Sa correction barotrope pour être strictement barocline:
            xy_t(i1,j,2)=                                               &
            xy_t(i1,j,2)                                                &
         +anyv3d(i1,j,k,1)                                              &
!    &   *dsig_x(i2,j,k)
           *dz_u(i2,j,k,1)                                              &
           /hz_u(i2,j,1)                                                & !01-10-09
       *  mask_u(i2,j,k)

       enddo

       do i=2,imax-1
! L'anomalie de pression avant correction barotrope:
           anyv3d(i,j1,k,1)=                                            &
            veldxdz_v(i,j2,k,0)                                         &
                /dz_v(i,j2,k,0)                                         &
                /dx_v(i,j2)                                            & !01-10-09
!    &    /dsig_y(i,j2,k)
!    &    /hzdx_y(i,j2,0)
          *const4

! Sa correction barotrope pour être strictement barocline:
            xy_t(i,j1,2)=                                               &
            xy_t(i,j1,2)                                                &
         +anyv3d(i,j1,k,1)                                              &
!    &   *dsig_y(i,j2,k)
         *dz_v(i,j2,k,1)                                                &
         /hz_v(i,j2,1)                                                  & !01-10-09
       *  mask_v(i,j2,k)
       enddo

      enddo

      do k=kmax,1,-1

       do j=2,jmax-1

! Ici calcul de la pression au point voisin I3. Mais attention cette pression
! doit servir de pression de reference et donc nous la voulons à la bonne
! profondeur de reference, c'est à dire z(I1). En z(I3) un effet de gradient
! pression serait induit par la pente des niveaux de calcul. Nous corrigeons
! l'effet de pente par un developpement de Taylor:
! densite(i3,z(i1))=densite(i3,z(i3))+deltaZ*(Ddensite/Dz) avec deltaZ=z(I1)-z(I3)
! et Ddensite/Dz donnee par le background: ( densiteobc(z(i1))-densiteobc(z(i3)) )
! / (z(i1)-z(i3)):

       x2=        g*dz_t(i1,j,k,1)                                      & !01-10-09
        *(const1*(tobc_t(i1,j,k,1)+tem_t(i3,j,k)-tobc_t(i3,j,k,1))      &
         +const2*(sobc_t(i1,j,k,1)+sal_t(i3,j,k)-sobc_t(i3,j,k,1))      &
         +const3)

       xy_t(i1,j,1)=xy_t(i1,j,1)+x2

! Dans GRDOBC3D il y a la reference. Celle ci est donnée par
! la valeur de pression et de courant aux points voisins des
! points de pression et courant impliqués dans la condition à la
! limite, filtrée de la haute fréquence. De sorte que pour un temps
! de filitrage tres court, cela revient à prendre un gradient nul pour
! p et u...
                grdobc3d_j(j,k,loop1)=                                  &
        const7* grdobc3d_j(j,k,loop1)                                   &
       +const8*(   xy_t(i1,j,1)-0.5*x2                                  &
                -(vel_u(i4,j,k,1)-velavr_u(i4,j,1))*const4  )

! pression totale = pression de reference + anomalie de pression:
       pressure_z(i1,j,k)=                                              &
          grdobc3d_j(j,k,loop1)                                         &
          +anyv3d(i1,j,k,1)                                             &
            -xy_t(i1,j,2)

       enddo

       do i=2,imax-1

       x2=        g*dz_t(i,j1,k,1)                                      & !01-10-09
        *(const1*(tobc_t(i,j1,k,1)+tem_t(i,j3,k)-tobc_t(i,j3,k,1))      &
         +const2*(sobc_t(i,j1,k,1)+sal_t(i,j3,k)-sobc_t(i,j3,k,1))      &
         +const3)

       xy_t(i,j1,1)=xy_t(i,j1,1)+x2

               grdobc3d_i(i,k,loop1)=                                   &
        const7*grdobc3d_i(i,k,loop1)                                    &
       +const8*(     xy_t(i,j1,1)-0.5*x2                                &
                  -(vel_v(i,j4,k,1)-velavr_v(i,j4,1))*const4  )

! pression totale = pression de reference + anomalie de pression:
       pressure_z(i,j1,k)=                                              &
       grdobc3d_i(i   ,k,loop1)                                         &
          +anyv3d(i,j1,k,1)                                             &
            -xy_t(i,j1,2)

       enddo

      enddo

      enddo ! fin boucle loop1
!----------------------------------------------------------------------
! Methode basique avec reference donnée par la basse frequence du model_
! Fin.
      return
      endif
!----------------------------------------------------------------------

      return
      end
