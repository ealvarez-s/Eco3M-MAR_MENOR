










      subroutine obc_ssh
!______________________________________________________________________
! SYMPHONIE ocean model
! release 299 - last update: 18-03-21
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer :: loop_,idi_ssh_z1_,idi_ssh_z0_,idi_ssh_za

!$ Description:
!$ ssh boundary conditions for far external grid nodes. This is not a momentum
!$ condition (the Flather is) but a simple procedure required to compute
!$ half/half average values on _X nodes and _Y nodes.


!...............................................................................
! Version date      Description des modifications
!         31-05-09  mis en service
! 2009.3  05-10-09  ajout d'un "ifdef parallele"
!         08-10-09  suppression de l'argument ichoix
! 2010.8  22-03-10  hssh_w remplace ssh_ext_w
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S.26    01-07-13  subcycle
!         25-03-14  Pour s'autoriser la possibilite de ne pas utiliser
!                   un gradient nul sur h_w mais seulement un gradient
!                   nul sur la ssh...
!         12-05-14  nouveaux echanges
!         02-07-14  nouveaux echanges
!         15-07-14  ajout subroutine obc_ssh_int_mpi
!                   + debug clef d'echange
!         25-09-14  Ajout subroutine obc_ssh_int_clamped 
!                   ajout subroutine obc_ssh_int_0grd
!         13-11-14  attention a la chronologie qui doit etre coherente
!                   avec les bornes des boucles
!         11-04-15  ajout subroutine obc_ssh_mpi
!         24-10-17  2 arguments dans la routine obc_ssh_int_mpi(txt_,t_)
! v259    01-10-19  call obc_ssh_flather !01-10-19
! v261    21-10-19  wetmask_u, wetmask_v
! v299    18-03-21  utiliser obcstatus
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................




! West open boundary (if any) !26-11-14
      if(obcstatus(ieq1)==1) then !----------->
       do j=1,jmax
         ssh_w(0     ,j,2)=ssh_w(1   ,j,2)
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> 
       do j=1,jmax
         ssh_w(imax+1,j,2)=ssh_w(imax,j,2)
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do i=0,imax+1
         ssh_w(i,0     ,2)=ssh_w(i,1   ,2)      
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do i=0,imax+1
         ssh_w(i,jmax+1,2)=ssh_w(i,jmax,2)
       enddo
      endif                                      !----------->

! mpi boundaries:
      call obc_ssh_mpi('z1',2)  !11-04-15

      end subroutine obc_ssh

!.................................................................

      subroutine obc_ssh_checkmpi !02-07-14
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer idi_ssh_z0,loop_

       call get_type_echange('z0','ssh_w_z0',ssh_w,lbound(ssh_w),ubound(ssh_w),2,idi_ssh_z0) ! 'za' = 'z0' et 'z1'
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(ssh_w,idi_ssh_z0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ssh_checkmpi

!.................................................................

      subroutine obc_sshobc_mpi(t_) !02-07-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,id_sshobc_w_z0_,loop_

      write(texte30,'(a12,i0)')'sshobc_w_z0_',t_    !15-07-14
       call get_type_echange('z0',trim(texte30)   & !15-07-14
                                  ,sshobc_w       &
                           ,lbound(sshobc_w)      &
                           ,ubound(sshobc_w)      &
                           ,t_                    &
                           ,id_sshobc_w_z0_)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(sshobc_w,id_sshobc_w_z0_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_sshobc_mpi

!.................................................................

      subroutine obc_ssh_int_mpi(txt_,t_) !24-10-17
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_
      character*2 txt_

      write(texte30,'(a,a,a,i0)')'ssh_int_w_',trim(txt_),'_',t_

       call get_type_echange(txt_,trim(texte30),ssh_int_w        &
                                        ,lbound(ssh_int_w)       &
                                        ,ubound(ssh_int_w),t_,k0)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(ssh_int_w,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ssh_int_mpi

!.................................................................

      subroutine obc_ssh_int_clamped !25-09-14
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_

! West open boundary (if any)
      if(obcstatus(ieq1)==1) then !----------->
       do j=1,jmax
         ssh_int_w(0,j,2)=sshobc_w(0,j,1)
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do i=0,imax+1
         ssh_int_w(i,0,2)=sshobc_w(i,0,1)
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !----------->
       do j=1,jmax
         ssh_int_w(imax+1,j,2)=sshobc_w(imax+1,j,1)
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do i=0,imax+1
         ssh_int_w(i,jmax+1,2)=sshobc_w(i,jmax+1,1)
       enddo
      endif                                      !----------->

      end subroutine obc_ssh_int_clamped !25-09-14

!.................................................................
      subroutine obc_ssh_int
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_

! West open boundary (if any) !26-11-14
      if(obcstatus(ieq1)==1) then !----------->
       do j=1,jmax
         ssh_int_w(0     ,j,2)=ssh_int_w(1   ,j,2)
       enddo
      endif                                      !----------->

! East open boundary (if any)
      if(obcstatus(ieqimax)==1) then !-----------> 
       do j=1,jmax
         ssh_int_w(imax+1,j,2)=ssh_int_w(imax,j,2)
       enddo
      endif                                      !----------->

! South open boundary (if any)
      if(obcstatus(jeq1)==1) then !----------->
       do i=0,imax+1
         ssh_int_w(i,0     ,2)=ssh_int_w(i,1   ,2)
       enddo
      endif                                      !----------->

! North open boundary (if any)
      if(obcstatus(jeqjmax)==1) then !----------->
       do i=0,imax+1
         ssh_int_w(i,jmax+1,2)=ssh_int_w(i,jmax,2)
       enddo
      endif                                      !----------->

      end subroutine obc_ssh_int

!.................................................................

      subroutine obc_ssh_mpi(txt_,t_) !11-04-15
      use module_principal
      use module_parallele
      implicit none
      integer t_,loop_
      character*2 txt_

      write(texte30,'(a,a,a,i0)')'ssh_w_',trim(txt_),'_',t_

       call get_type_echange(txt_,trim(texte30),ssh_w        &
                                        ,lbound(ssh_w)       &
                                        ,ubound(ssh_w),t_,k0)
      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(ssh_w,k0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine obc_ssh_mpi

!.................................................................

      subroutine obc_ssh_flather !01-10-19
      use module_principal ; use module_parallele ; use module_wave
      implicit none

!      sshlfwindow=1.-dti_fw/100.
!      cumultimesshlf=sshlfwindow*cumultimesshlf+1. ! dti_fw on somme 1 (et pas dti si sshlf fait une simple addition)
      
      if(obcstatus(jeq1)==1) then                !----------->
       j1=1     ! ssh obc
       j2=j1+1  ! v voisin 1
       j3=j1+1  ! ssh et u voisins 1
       j4=j2+1  ! v voisin 2
       j5=j3+1  ! ssh et u voisins 2
       j6=j4+1  ! v voisin 3
       j7=j5+1  
       j8=j6+1  ! v voisin 4
       j9=j7+1  
       loop1=1
       call obc_ssh_flather_northsouth
      endif                                      !----------->

      if(obcstatus(jeqjmax)==1) then             !----------->
       j1=jmax  ! ssh obc
       j2=j1    ! v voisin 1
       j3=j1-1  ! ssh et u voisins 1
       j4=j2-1  ! v voisin 2
       j5=j3-1  ! ssh et u voisins 2
       j6=j4-1  ! v voisin 3
       j7=j5-1  
       j8=j6-1  ! v voisin 4
       j9=j7-1  
       loop1=2
       call obc_ssh_flather_northsouth
      endif                                      !----------->

      if(obcstatus(ieqimax)==1) then             !----------->
       i1=imax     ! ssh obc
       i2=i1       ! u voisin 1
       i3=i1-1     ! ssh et v voisins 1
       i4=i2-1     ! u voisin 2
       i5=i3-1     ! ssh et v voisins 2
       i6=i4-1     ! u voisin 3
       i7=i5-1    
       i8=i6-1     ! u voisin 4
       i9=i7-1    
       loop1=2
       call obc_ssh_flather_eastwest
      endif                                      !----------->

      if(obcstatus(ieq1)==1) then                !----------->
       i1=1     ! ssh obc
       i2=i1+1  ! u voisin 1
       i3=i1+1  ! ssh et v voisins 1
       i4=i2+1  ! u voisin 2
       i5=i3+1  ! ssh et v voisins 2
       i6=i4+1  ! u voisin 3
       i7=i5+1    
       i8=i6+1  ! u voisin 4
       i9=i7+1    
       loop1=1
       call obc_ssh_flather_eastwest
      endif                                      !----------->

      end subroutine obc_ssh_flather

!.................................................................

      subroutine obc_ssh_flather_eastwest
      use module_principal ; use module_parallele ; use module_wave
      implicit none

!     stop 'obc_ssh_flather_eastwest'
! ICI UNE BIDOUILLE PATRICK POUR FAIRE UNE SORTE DE CONDITION
! PERIODIQUE ALONG I, bidouille requerant un decoupage mpi
! en tranches horizontale d'ou le test sur nbdom_imax=1
!      if(nbdom_imax/=1)stop 'Erreur nbdom_imax/=1'
!      if(loop1==2) then
!       do j=1,jmax
!           sshobc_w(i1,j,1)=0.
!           sshobc_w(i1,j,2)=0.
!           sshobc_w(i3,j,1)=0.
!           sshobc_w(i3,j,2)=0.
!           sshobc_w(i5,j,2)=0.
!       enddo
!      else
!       do j=1,jmax
!           sshobc_w(i1,j,1)=ssh_int_w(imax-3,j,1) 
!           sshobc_w(i1,j,2)=ssh_int_w(imax-3,j,2)
!           sshobc_w(i3,j,1)=ssh_int_w(imax-2,j,1)  
!           sshobc_w(i3,j,2)=ssh_int_w(imax-2,j,2) 
!           sshobc_w(i5,j,2)=ssh_int_w(imax-1,j,2)
!       enddo
!      endif

      k=kmax
      obc_scheme=6
!     x0=1. ! amplification pour tenir compte angle
      if(obc_scheme==0) then !000>
! FERMEE:
       do j=1,jmax
        ssh_int_w(i1,j,2)=0.
       enddo
      endif                  !000>

      if(obc_scheme==1) then !111>
! FLA 
       do j=1,jmax
         ssh_int_w(i1,j,2)=     &
            mask_t(i1,j,kmax)*( & !ooo>
                           sshwave_j_w(j,1,2,loop1) &
        +(i1-i3)*freq2pipeak/(kvectorpeak_j(j,loop1)*grav)*(vel_u(i2,j,k,2)-velwave_j_u(j,k,1,2,loop1)) & 
                              )   !ooo>
       enddo
      endif                  !111>

      if(obc_scheme==2) then !222>
! GDR 1 FLA NORMAL
       do j=1,jmax
        ssh_int_w(i1,j,2)= &
        ssh_int_w(i3,j,2)  &
       +(i1-i3)*(vel_u(i2,j,k,2)-vel_u(i4,j,k,2)) &
            *sinh(kvectorpeak_j(j,loop1)*max(0.,                h_w(i1,j))) &
            /cosh(kvectorpeak_j(j,loop1)*max(1.,depth_t(i1,j,k)+h_w(i1,j))) &
            /freq2pipeak*x0 
       enddo
      endif                  !222>

      if(obc_scheme==3) then !333>
! GDR 2 FLA NORMAL
       do j=1,jmax
         ssh_int_w(i1,j,2)= &
       2*ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)  &
           +(i1-i3)*(vel_u(i2,j,k,2)-2*vel_u(i4,j,k,2)+vel_u(i6,j,k,2)) &
            *sinh(kvectorpeak_j(j,loop1)*max(0.,                h_w(i1,j))) &
            /cosh(kvectorpeak_j(j,loop1)*max(1.,depth_t(i1,j,k)+h_w(i1,j))) &
            /freq2pipeak*x0       
       enddo
      endif                  !333>

      if(obc_scheme==23) then !232323>
! MIX GDR1 GDR 2 FLA NORMAL
       do j=1,jmax
         ssh_int_w(i1,j,2)= &
        0.5*(ssh_int_w(i3,j,2)  &
       +(i1-i3)*(vel_u(i2,j,k,2)-vel_u(i4,j,k,2)) &
            *sinh(kvectorpeak_j(j,loop1)*max(0.,                h_w(i1,j))) &
            /cosh(kvectorpeak_j(j,loop1)*max(1.,depth_t(i1,j,k)+h_w(i1,j))) &
            /freq2pipeak ) &
       +0.5*(2*ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)  &
           +(i1-i3)*(vel_u(i2,j,k,2)-2*vel_u(i4,j,k,2)+vel_u(i6,j,k,2)) &
            *sinh(kvectorpeak_j(j,loop1)*max(0.,                h_w(i1,j))) &
            /cosh(kvectorpeak_j(j,loop1)*max(1.,depth_t(i1,j,k)+h_w(i1,j))) &
            /freq2pipeak )      
       enddo
      endif                  !232323>

!     if(obc_scheme==4) then !444>
! GDR 3 FLA NORMAL
!      do j=1,jmax
!        ssh_int_w(i1,j,2)=min(max( &
!      3*ssh_int_w(i3,j,2)-3*ssh_int_w(i5,j,2)+ssh_int_w(i7,j,2) &
!     +(i1-i3)*(vel_u(i2,j,k,2)-3*vel_u(i4,j,k,2)+3*vel_u(i6,j,k,2)-vel_u(i8,j,k,2)) &
!           *sinh(kvectorpeak_j(j,loop1)*max(0.,                h_w(i1,j))) &
!           /cosh(kvectorpeak_j(j,loop1)*max(1.,depth_t(i1,j,k)+h_w(i1,j))) &
!           /freq2pipeak      &
!        ,min(ssh_int_w(i1,j,1),ssh_int_w(i3,j,1))) &
!        ,max(ssh_int_w(i1,j,1),ssh_int_w(i3,j,1)))  
!       ssh_int_w(i1,j,2)= &
!       0. *(i1-i3)*vel_u(i2,j,k,2) &
!           *sinh(kvectorpeak_j(j,loop1)*max(0.,                h_w(i1,j))) &
!           /cosh(kvectorpeak_j(j,loop1)*max(1.,depth_t(i1,j,k)+h_w(i1,j))) &
!           /freq2pipeak &
!       +1. *(3*ssh_int_w(i3,j,2)-3*ssh_int_w(i5,j,2)+ssh_int_w(i7,j,2) )
!      enddo
!     endif                  !444>

      if(obc_scheme==4) then !444>
! Sommerfeld implicite
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1))*x0

         ssh_int_w(i1,j,2)=  &
          sshobc_w(i1,j,2)+( &        !ooo>
                            ssh_int_w(i1,j,1)-sshobc_w(i1,j,1)  &
                       +x1*(ssh_int_w(i3,j,2)-sshobc_w(i3,j,2)) &
                           )/(1.+x1)  !ooo> 
       enddo
      endif                  !444>

      if(obc_scheme==5) then !555>
! Sommerfeld Gradient
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1))*x0

! ecriture 1 
!        ssh_int_w(i1,j,2)=(ssh_int_w(i1,j,1)+x1*ssh_int_w(i3,j,2))/(1.+x1) &
!       +ssh_int_w(i3,j,2)-(ssh_int_w(i3,j,1)+x1*ssh_int_w(i5,j,2))/(1.+x1)  
! ecriture 2 (equivalente):
           ssh_int_w(i1,j,2)=                  &
            sshobc_w(i1,j,2)                   &
          +ssh_int_w(i3,j,2)-sshobc_w(i3,j,2)  &
         +(ssh_int_w(i1,j,1)-sshobc_w(i1,j,1)-ssh_int_w(i3,j,1)+sshobc_w(i3,j,1) &
      +x1*(ssh_int_w(i3,j,2)-sshobc_w(i3,j,2)-ssh_int_w(i5,j,2)+sshobc_w(i5,j,2)))/(1.+x1)

!       if(loop1==1) then
!           ssh_int_w(i1,j,2)=ssh_int_w(imax-1,j,2)
!       endif

       enddo
      endif                  !555>

      if(obc_scheme==6) then !666>
       k=kmax
       do j=1,jmax

!.................................
! FLA , SOMMERFELD , GDR0
!       x2=(i1-i3)*freq2pipeak/(kvectorpeak_j(j,loop1)*grav)   ! c/g=cwavepeak/g
        x2=(i1-i3)*sqr_hoverg_u(j,loop1)                       ! c/g=sqrt(gh)/g=sqrt(h/g)

! x1 est fixe:
!       x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1)) ! cpeak*dt/dx
!       x1=grav*sqr_hoverg_u(j,loop1)*invdx_t(i1,j)*dti_fw        ! sqrt(gh)*dt/dx=g*sqrt(h/g)*dt*inv_dx
! Methode Orl 1
         x1=min(  abs((ssh_int_w(i3,j,2)-sshwave_j_w(j,2,2,loop1))-(ssh_int_w(i3,j,1)-sshwave_j_w(j,2,1,loop1))) &
                /(abs((ssh_int_w(i3,j,2)-sshwave_j_w(j,2,2,loop1))-(ssh_int_w(i5,j,2)-sshwave_j_w(j,3,2,loop1)))+small1) , 100. )


              ssh_int_w(i1,j,2)=max(-h_w(i1,j), & !--->
                 mask_t(i1,j,kmax)*( & !ooo>

! Condition de Flather:
               sshwave_j_w(j,1,2,loop1) +x2*(vel_u(i2,j,k,2)-velwave_j_u(j,k,1,2,loop1))   &

      +wetmask_u(i2,j)*wetmask_u(i4,j)*( & ! conditions additionnelles > !21-10-19

! + condition Sommerfeld
       (         &
             (ssh_int_w(i1,j,1)-sshwave_j_w(j,1,1,loop1))-x2*(vel_u(i2,j,k,1)-velwave_j_u(j,k,1,1,loop1))           &
        +x1*((ssh_int_w(i3,j,2)-sshwave_j_w(j,2,2,loop1))-x2*(vel_u(i4,j,k,2)-velwave_j_u(j,k,2,2,loop1)))          &
                                                  )/(1.+x1)   & 

! + condition GDR1
       -(   -(ssh_int_w(i3,j,2)-sshwave_j_w(j,2,2,loop1))+x2*(vel_u(i4,j,k,2)-velwave_j_u(j,k,2,2,loop1))+(        &
             (ssh_int_w(i3,j,1)-sshwave_j_w(j,2,1,loop1))-x2*(vel_u(i4,j,k,1)-velwave_j_u(j,k,2,1,loop1))          &
        +x1*((ssh_int_w(i5,j,2)-sshwave_j_w(j,3,2,loop1))-x2*(vel_u(i6,j,k,2)-velwave_j_u(j,k,3,2,loop1)))         &
                                                  )/(1.+x1) )    & 

                                       ) & ! conditions additionnelles >
                                   ) & !ooo>
                                   )              !--->

! Force de rappel sur la ssh moyenne
! Methode 1:
!    -sshlf_j_w(j,loop1)*0.01*dti_fw/cumultimesshlf   & 
!     *max(0.,sign(1.,(ssh_int_w(i1,j,2)-sshwave_j_w(j,1,2,loop1))  &
!     *(sshlf_j_w(j,loop1)/cumultimesshlf) ))   
! Methode 2 (du courant dans la BOX):
!     -dti_fw*hmax*velbar_u(i1,j,2)*0.001 ! 0.001 = inverse d'une distance A la ssh de reference "lointaine"
! Simple rappel
!     +sshwave_j_w(j,1,2,loop1)*dti_fw*0.01 &
!       ) & !ooo>rappel>
!                          /(1.+dti_fw*0.01)

       enddo
      endif                  !666>

      if(obc_scheme==7) then !777>
! Sommerfeld Gradient directionnel
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
!       x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1)) &
!     *sqrt( (abs(vel_u(i2,j,k,2)-vel_u(i4,j,k,2))+abs(vel_v(i1,j+1,k,2)-vel_v(i1,j,k,2))) &
!           /(abs(vel_u(i2,j,k,2)-vel_u(i4,j,k,2))+0.00001)  )

! Methode Orl
         x1=min(  abs(ssh_int_w(i3,j,2)-ssh_int_w(i3,j,1)) &
                /(abs(ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2))+small1) , 100. )

         ssh_int_w(i1,j,2)=  &
         ssh_int_w(i3,j,2)   &
       +(ssh_int_w(i1,j,1)-ssh_int_w(i3,j,1)+x1*(ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)))/(1.+x1)


       enddo
      endif                  !777>

      if(obc_scheme==8) then !888>
! Sommerfeld directionnel "explicite+borne"
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=min(1.,freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1)) &
      *sqrt( (abs(vel_u(i2,j,k,2)-vel_u(i4,j,k,2))+abs(vel_v(i1,j+1,k,2)-vel_v(i1,j,k,2))) &
            /(abs(vel_u(i2,j,k,2)-vel_u(i4,j,k,2))+0.00001)  ))

         ssh_int_w(i1,j,2)=(1-x1)*ssh_int_w(i1,j,1)+x1*ssh_int_w(i3,j,2)

       enddo
      endif                  !888>

      if(obc_scheme==9) then !999>
! Sommerfeld explicite
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1))

         ssh_int_w(i1,j,2)=(1-x1)*ssh_int_w(i1,j,1)+x1*ssh_int_w(i3,j,2)

       enddo
      endif                  !999>

      if(obc_scheme==10) then !1010>
! Sommerfeld gradient explicite
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1))

         ssh_int_w(i1,j,2)=   &
         ssh_int_w(i3,j,2)    &
         +(1-x1)*(ssh_int_w(i1,j,1)-ssh_int_w(i3,j,1))+x1*(ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2))

       enddo
      endif                  !1010>

      if(obc_scheme==11) then !11-11>
! Sommerfeld gradient directionnel "explicite+borne"
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=min(1.,freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1)) &
      *sqrt( (abs(vel_u(i2,j,k,2)-vel_u(i4,j,k,2))+abs(vel_v(i1,j+1,k,2)-vel_v(i1,j,k,2))) &
            /(abs(vel_u(i2,j,k,2)-vel_u(i4,j,k,2))+0.00001)  ))

         ssh_int_w(i1,j,2)=   &
         ssh_int_w(i3,j,2)    &
         +(1-x1)*(ssh_int_w(i1,j,1)-ssh_int_w(i3,j,1))+x1*(ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2))

       enddo
      endif                  !11-11>

      if(obc_scheme==12) then !1212>
! Sommerfeld derivee seconde implicite
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1))*x0

         ssh_int_w(i1,j,2)=                  &
       2*ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)+( & !pmx>

            ssh_int_w(i1,j,1)-2*ssh_int_w(i3,j,1)+ssh_int_w(i5,j,1)  &
       +x1*(ssh_int_w(i3,j,2)-2*ssh_int_w(i5,j,2)+ssh_int_w(i7,j,2)) &

                                             ) & !pmx>
                                              /(1.+x1)

       enddo
      endif                  !1212>

      if(obc_scheme==13) then !1313>
! Gradient nul
!      stop 'schema 13'
       do j=1,jmax
!        ssh_int_w(i1,j,2)=ssh_int_w(i3,j,2)   ! Gradient normal  nul
!        ssh_int_w(i1,j,2)=ssh_int_w(i3,j-1,1) ! Gradient oblique nul
         ssh_int_w(i1,j,2)=ssh_int_w(i3,j,2)+ssh_int_w(i3,j-1,1)-ssh_int_w(i5,j-1,1) ! mix des 2
       enddo
      endif                   !1313>
      if(obc_scheme==14) then !1414>
! Derivee seconde nulle
       do j=1,jmax
         ssh_int_w(i1,j,2)=2*ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)
       enddo
      endif                   !1414>

      if(obc_scheme==15) then !1515>
! Sommerfeld derivee seconde implicite "moyenne sur 2 points"
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1))

         ssh_int_w(i1,j,2)=0.5*(                  &
       2*ssh_int_w(i3,j,2)+ssh_int_w(i5,j,2)-ssh_int_w(i7,j,2)+( & !pmx>
            ssh_int_w(i1,j,1)-ssh_int_w(i3,j,1)-ssh_int_w(i5,j,1)+ssh_int_w(i7,j,1)  &
       +x1*(ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)-ssh_int_w(i7,j,2)+ssh_int_w(i9,j,1)) &

                                                               ) & !pmx>
                                              /(1.+x1) )

       enddo
      endif                  !1515>

      if(obc_scheme==16) then !1616>
! Sommerfeld derivee seconde implicite C moderee
       do j=1,jmax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dx_t(i1,j)*kvectorpeak_j(j,loop1)) &
          *0.5

! ecriture 1 
!        ssh_int_w(i1,j,2)=  &
!        ssh_int_w(i3,j,2)   &
!      +(ssh_int_w(i1,j,1)-ssh_int_w(i3,j,1)+x1*(ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)))/(1.+x1) &

!     -(-ssh_int_w(i3,j,2)   &
!       +ssh_int_w(i5,j,2)   &
!      +(ssh_int_w(i3,j,1)-ssh_int_w(i5,j,1)+x1*(ssh_int_w(i5,j,2)-ssh_int_w(i7,j,2)))/(1.+x1))
! ecriture 2 equivalente
         ssh_int_w(i1,j,2)=                  &
       2*ssh_int_w(i3,j,2)-ssh_int_w(i5,j,2)+( & !pmx>

            ssh_int_w(i1,j,1)-2*ssh_int_w(i3,j,1)+ssh_int_w(i5,j,1)  &
       +x1*(ssh_int_w(i3,j,2)-2*ssh_int_w(i5,j,2)+ssh_int_w(i7,j,2)) &

                                             ) & !pmx>
                                              /(1.+x1)

       enddo
      endif                  !1616>


      end subroutine obc_ssh_flather_eastwest

!.................................................................

      subroutine obc_ssh_flather_northsouth
      use module_principal ; use module_parallele ; use module_wave
      implicit none

!     stop 'obc_ssh_flather_northsouth'



      obc_scheme=1
      k=kmax

      if(obc_scheme==0) then !000>
! FERMEE:
       do i=1,imax
        ssh_int_w(i,j1,2)=0.
       enddo
      endif                  !000>

      if(obc_scheme==1) then !111>

!     if(.not.allocated(kvectorpeak_i))stop 'allocated(kvectorpeak_i)'
!     if(.not.allocated(velwave_i_v)) stop 'allocate velwave_i_v'

! FLA 
       do i=1,imax
        ssh_int_w(i,j1,2)= &
           mask_t(i,j1,kmax)*( & !ooo>
                          sshwave_i_w(i,1,2,loop1) &
       +(j1-j3)*freq2pipeak/(kvectorpeak_i(i,loop1)*grav)*(vel_v(i,j2,k,2)-velwave_i_v(i,k,1,2,loop1)) &
                             )   !ooo>
!       if(i+par%timax(1)==iglb/2.and.loop1==2)write(10+par%rank,*)real(ssh_int_w(i,j1,2)),sshwave_i_w(i,1,2,loop1)
       enddo

      endif                  !111>

      if(obc_scheme==2) then !222>
! GDR 1 FLA NORMAL
       do i=1,imax
        ssh_int_w(i,j1,2)= &
        ssh_int_w(i,j3,2)  &
       +(j1-j3)*(vel_v(i,j2,k,2)-vel_v(i,j4,k,2)) &
            *sinh(kvectorpeak_i(i,loop1)*max(0.,                h_w(i,j1))) &
            /cosh(kvectorpeak_i(i,loop1)*max(1.,depth_t(i,j1,k)+h_w(i,j1))) &
            /freq2pipeak 
       enddo
      endif                  !222>

      if(obc_scheme==3) then !333>
! GDR 2 FLA NORMAL
       do i=1,imax
         ssh_int_w(i,j1,2)= &
       2*ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)  & 
           +(j1-j3)*(vel_v(i,j2,k,2)-2*vel_v(i,j4,k,2)+vel_v(i,j6,k,2)) &
            *sinh(kvectorpeak_i(i,loop1)*max(0.,                h_w(i,j1))) &
            /cosh(kvectorpeak_i(i,loop1)*max(1.,depth_t(i,j1,k)+h_w(i,j1))) &
            /freq2pipeak       
       enddo
      endif                  !333>

!     if(obc_scheme==4) then !444>
! GDR 3 FLA NORMAL
!      do i=1,imax

!        ssh_int_w(i,j1,2)= min(max( &
!      3*ssh_int_w(i,j3,2)-3*ssh_int_w(i,j5,2)+ssh_int_w(i,j7,2)  & 
!     +(j1-j3)*(vel_v(i,j2,k,2)-3*vel_v(i,j4,k,2)+3*vel_v(i,j6,k,2)-vel_v(i,j8,k,2)) &
!           *sinh(kvectorpeak_i(i,loop1)*max(0.,                h_w(i,j1))) &
!           /cosh(kvectorpeak_i(i,loop1)*max(1.,depth_t(i,j1,k)+h_w(i,j1))) &
!           /freq2pipeak    &  
!       ,min(ssh_int_w(i,j1,1),ssh_int_w(i,j3,1))) &
!       ,max(ssh_int_w(i,j1,1),ssh_int_w(i,j3,1)))  

!      enddo
!     endif                  !444>

      if(obc_scheme==23) then !232323>
! 50% GDR1 50% GDR2
       do i=1,imax
        ssh_int_w(i,j1,2)=  &

       0.5*(ssh_int_w(i,j3,2)  &
       +(j1-j3)*(vel_v(i,j2,k,2)-vel_v(i,j4,k,2)) &
            *sinh(kvectorpeak_i(i,loop1)*max(0.,                h_w(i,j1))) &
            /cosh(kvectorpeak_i(i,loop1)*max(1.,depth_t(i,j1,k)+h_w(i,j1))) &
            /freq2pipeak ) &

      +0.5*(2*ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)  &
           +(j1-j3)*(vel_v(i,j2,k,2)-2*vel_v(i,j4,k,2)+vel_v(i,j6,k,2)) &
            *sinh(kvectorpeak_i(i,loop1)*max(0.,                h_w(i,j1))) &
            /cosh(kvectorpeak_i(i,loop1)*max(1.,depth_t(i,j1,k)+h_w(i,j1))) &
            /freq2pipeak )   
                               
       enddo
      endif                  !232323>

      if(obc_scheme==4) then !444>
! Sommerfeld implicit
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1))

         ssh_int_w(i,j1,2)=(ssh_int_w(i,j1,1)+x1*ssh_int_w(i,j3,2))/(1.+x1)  

       enddo
      endif                  !444>

      if(obc_scheme==5) then !555>
! Sommerfeld Gradient implicite
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) !&
!          *0.75 ! reducteur de vitesse

! ecriture 1:
!        ssh_int_w(i,j1,2)=(ssh_int_w(i,j1,1)+x1*ssh_int_w(i,j3,2))/(1.+x1) &
!       +ssh_int_w(i,j3,2)-(ssh_int_w(i,j3,1)+x1*ssh_int_w(i,j5,2))/(1.+x1)  
! ecriture 2 (equivalente):
         ssh_int_w(i,j1,2)=  &
         ssh_int_w(i,j3,2)   &
       +(ssh_int_w(i,j1,1)-ssh_int_w(i,j3,1)+x1*(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)))/(1.+x1)

       enddo
      endif                  !555>

      if(obc_scheme==6) then !666>
       k=kmax
! Sommerfeld directionnel
       do i=1,imax

!.................................
! FLA , SOMMERFELD , GDR0
        x2=(j1-j3)*freq2pipeak/(kvectorpeak_i(i,loop1)*grav) ! c/g=cwavepeak/g
!       x2=(j1-j3)*sqr_hoverg_v(i,loop1)                     ! c/g=sqrt(gh)/g=sqrt(h/g)

! x1 est fixe
!       x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) ! cwavepeak*dt/dy
!       x1=grav*sqr_hoverg_v(i,loop1)*invdy_t(i,j1)*dti_fw        ! sqrt(gh)*dt/dy=g*sqrt(h/g)*dt*inv_dy
! Methode Orl 1
         x1=min(  abs((ssh_int_w(i,j3,2)-sshwave_i_w(i,2,2,loop1))-(ssh_int_w(i,j3,1)-sshwave_i_w(i,2,1,loop1))) &
                /(abs((ssh_int_w(i,j3,2)-sshwave_i_w(i,2,2,loop1))-(ssh_int_w(i,j5,2)-sshwave_i_w(i,3,2,loop1)))+small1) , 100. )

! Terme basse frequence:
!     sshlf_i_w(i,loop1)=sshlfwindow*sshlf_i_w(i,loop1)+ssh_int_w(i,j1,2)-sshwave_i_w(i,1,2,loop1)

              ssh_int_w(i,j1,2)=max(-h_w(i,j1),  & !---> !17-10-19
                 mask_t(i,j1,kmax)*( & !msk>

! Condition de Flather:
                                sshwave_i_w(i,1,2,loop1) +x2*(vel_v(i,j2,k,2)-velwave_i_v(i,k,1,2,loop1))           &

      +wetmask_v(i,j2)*wetmask_v(i,j4)*( & ! conditions additionnelles > !21-10-19

! + condition Sommerfeld
       (          &
             (ssh_int_w(i,j1,1)-sshwave_i_w(i,1,1,loop1))-x2*(vel_v(i,j2,k,1)-velwave_i_v(i,k,1,1,loop1))           &
        +x1*((ssh_int_w(i,j3,2)-sshwave_i_w(i,2,2,loop1))-x2*(vel_v(i,j4,k,2)-velwave_i_v(i,k,2,2,loop1)))          &
                                                  )/(1.+x1)  & 

! + condition GDR1
       -(   -(ssh_int_w(i,j3,2)-sshwave_i_w(i,2,2,loop1))+x2*(vel_v(i,j4,k,2)-velwave_i_v(i,k,2,2,loop1))+(        &
             (ssh_int_w(i,j3,1)-sshwave_i_w(i,2,1,loop1))-x2*(vel_v(i,j4,k,1)-velwave_i_v(i,k,2,1,loop1))          &
        +x1*((ssh_int_w(i,j5,2)-sshwave_i_w(i,3,2,loop1))-x2*(vel_v(i,j6,k,2)-velwave_i_v(i,k,3,2,loop1)))         &
                                                  )/(1.+x1) )   &  

                                      ) & ! conditions additionnelles >
                                   ) & !msk>
                                   )               !--->

!    -sshlf_i_w(i,loop1)*0.01*dti_fw/cumultimesshlf   & 
!     *max(0.,sign(1.,(ssh_int_w(i,j1,2)-sshwave_i_w(i,1,2,loop1))  &
!     *(sshlf_i_w(i,loop1)/cumultimesshlf) ))   

!.................................
! Lignes supplementaires pour gradient oblique
!       -( &
!           -ssh_int_w(i,j1-1,1)+x2*vel_v(i,j2-1,k,1)+(         &
!            ssh_int_w(i,j1-1,0)-x2*vel_v(i,j2-1,k,0)           &
!       +x1*(ssh_int_w(i,j3-1,1)-x2*vel_v(i,j4-1,k,1))          &
!                                                 )/(1.+x1) &

! lignes supplementaires pour imbriquer en plus une condition GDR1
!      -(   -ssh_int_w(i,j3-1,1)+x2*vel_v(i,j4-1,k,1)+(        &
!            ssh_int_w(i,j3-1,0)-x2*vel_v(i,j4-1,k,0)          &
!       +x1*(ssh_int_w(i,j5-1,1)-x2*vel_v(i,j6-1,k,1))         &
!                                                 )/(1.+x1) )  &
!        )

!        if(i+par%timax(1)==iglb/2.and.loop1==2)write(10+par%rank,*) &
!        real(elapsedtime_now) &
!!      ,real(ssh_int_w(i,j1,2)),sshwave_i_w(i,1,2,loop1)   
!       ,real(sshlf_i_w(i,loop1)/cumultimesshlf) &
!     ,real(  max(0.,sign(1.,(ssh_int_w(i,j1,2)-sshwave_i_w(i,1,2,loop1))  &
!     *(sshlf_i_w(i,loop1)/cumultimesshlf) )) ) 

       enddo
      endif                  !666>

      if(obc_scheme==7) then !777>
! Sommerfeld Gradient directionnel
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
!       x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) &
!     *sqrt( (abs(vel_v(i,j2,k,2)-vel_v(i,j4,k,2))+abs(vel_u(i,j1,k,2)-vel_u(i+1,j1,k,2))) &
!           /(abs(vel_v(i,j2,k,2)-vel_v(i,j4,k,2))+0.00001)  )

! Methode Orl
         x1=min(  abs(ssh_int_w(i,j3,2)-ssh_int_w(i,j3,1)) &
                /(abs(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2))+small1) , 100. )

         ssh_int_w(i,j1,2)=  &
         ssh_int_w(i,j3,2)   &
       +(ssh_int_w(i,j1,1)-ssh_int_w(i,j3,1)+x1*(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)))/(1.+x1)

       enddo
      endif                  !777>

      if(obc_scheme==8) then !888>
! Sommerfeld directionnel "explicite+borne"
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=min(1.,freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) &
      *sqrt( (abs(vel_v(i,j2,k,2)-vel_v(i,j4,k,2))+abs(vel_u(i+1,j1,k,2)-vel_u(i,j1,k,2))) &
            /(abs(vel_v(i,j2,k,2)-vel_v(i,j4,k,2))+0.00001)  ))

         ssh_int_w(i,j1,2)=(1-x1)*ssh_int_w(i,j1,1)+x1*ssh_int_w(i,j3,2)

       enddo
      endif                  !888>

      if(obc_scheme==9) then !999>
! Sommerfeld explicit
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1))

         ssh_int_w(i,j1,2)=(1-x1)*ssh_int_w(i,j1,1)+x1*ssh_int_w(i,j3,2)

       enddo
      endif                  !999>

      if(obc_scheme==10) then !1010>
! Sommerfeld gradient explicite
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1))

         ssh_int_w(i,j1,2)=   &
         ssh_int_w(i,j3,2)    &
         +(1-x1)*(ssh_int_w(i,j1,1)-ssh_int_w(i,j3,1))+x1*(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2))

       enddo
      endif                   !1010>

      if(obc_scheme==11) then !1111>
! Sommerfeld gradient directionnel "explicite+borne"
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=min(1.,freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) &
      *sqrt( (abs(vel_v(i,j2,k,2)-vel_v(i,j4,k,2))+abs(vel_u(i+1,j1,k,2)-vel_u(i,j1,k,2))) &
            /(abs(vel_v(i,j2,k,2)-vel_v(i,j4,k,2))+0.00001)  ))

         ssh_int_w(i,j1,2)=   &
         ssh_int_w(i,j3,2)    &
         +(1-x1)*(ssh_int_w(i,j1,1)-ssh_int_w(i,j3,1))+x1*(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2))

       enddo
      endif                   !1111>

      if(obc_scheme==12) then !1212>
! Sommerfeld derivee seconde implicite
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) ! &
!           *0.75 ! reducteur de vitesse

! ecriture 1 
!        ssh_int_w(i,j1,2)=  &
!        ssh_int_w(i,j3,2)   &
!      +(ssh_int_w(i,j1,1)-ssh_int_w(i,j3,1)+x1*(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)))/(1.+x1) &

!     -(-ssh_int_w(i,j3,2)   &
!       +ssh_int_w(i,j5,2)   &
!      +(ssh_int_w(i,j3,1)-ssh_int_w(i,j5,1)+x1*(ssh_int_w(i,j5,2)-ssh_int_w(i,j7,2)))/(1.+x1))

! ecriture 2 equivalente
         ssh_int_w(i,j1,2)=                  &
       2*ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)+( & !pmx>

            ssh_int_w(i,j1,1)-2*ssh_int_w(i,j3,1)+ssh_int_w(i,j5,1)  &
       +x1*(ssh_int_w(i,j3,2)-2*ssh_int_w(i,j5,2)+ssh_int_w(i,j7,2)) &

                                             ) & !pmx>
                                              /(1.+x1)

       enddo
      endif                  !1212>

      if(obc_scheme==13) then !1313>
! Gradient nul
       do i=1,imax
         ssh_int_w(i,j1,2)=ssh_int_w(i,j3,2)
       enddo
      endif                   !1313>
      if(obc_scheme==14) then !1414>
! Derivee seconde nulle
       do i=1,imax
         ssh_int_w(i,j1,2)=2*ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)
       enddo
      endif                   !1414>

      if(obc_scheme==15) then !1515>
! Sommerfeld derivee seconde implicite "moyenne sur 2 points"
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1))  


         ssh_int_w(i,j1,2)=0.5*(                  &
       2*ssh_int_w(i,j3,2)+ssh_int_w(i,j5,2)-ssh_int_w(i,j7,2)+( & !pmx>
            ssh_int_w(i,j1,1)-ssh_int_w(i,j3,1)-ssh_int_w(i,j5,1)+ssh_int_w(i,j7,1)  &
       +x1*(ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)-ssh_int_w(i,j7,2)+ssh_int_w(i,j9,1)) &

                                                               ) & !pmx>
                                              /(1.+x1) )
       enddo
      endif                  !1515>

      if(obc_scheme==16) then !1616>
! Sommerfeld derivee seconde implicite C moderee
       do i=1,imax
! "r"= dt/dx*c avec c=w/k
        x1=freq2pipeak*dti_fw/(dy_t(i,j1)*kvectorpeak_i(i,loop1)) &
            *0.5 ! reducteur de vitesse

! ecriture 2 equivalente
         ssh_int_w(i,j1,2)=                  &
       2*ssh_int_w(i,j3,2)-ssh_int_w(i,j5,2)+( & !pmx>

            ssh_int_w(i,j1,1)-2*ssh_int_w(i,j3,1)+ssh_int_w(i,j5,1)  &
       +x1*(ssh_int_w(i,j3,2)-2*ssh_int_w(i,j5,2)+ssh_int_w(i,j7,2)) &

                                             ) & !pmx>
                                              /(1.+x1)

       enddo
      endif                  !1616>


!     endif             !>>>>>>>> ! test sur loop

      end subroutine obc_ssh_flather_northsouth
