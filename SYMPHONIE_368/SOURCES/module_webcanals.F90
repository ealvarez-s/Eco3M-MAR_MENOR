      module module_webcanals
!______________________________________________________________________
! SYMPHONIE ocean model
! release 311 - last update: 16-11-21
!______________________________________________________________________
! version   date    description                                                !
! v252    12-04-19  mise en service                                            !
! v253    02-05-19  subroutine webcanals_messages
!                   subroutine webcanals_gt2_to_gt1_any2d
!         04-05-19  permettre aux tetes de canaux d'etre dans la zone de chevauchement mpi
!         05-05-19  Puisque la vitesse "limite" n'est pas envoyee par gridtype2, gridtype1 
!                   construit lui meme sa condition aux limite "gradient nul" !05-05-19
!         15-05-19  suite du point precedent. La C.L. de gradient nul peut amplifier
!                   une resonance en situation de fort courant (situation frequente)
!                   on revient A la condition de courant nul (par defaut)
! v257    23-06-19  ajout subroutine webcanals_gt2_to_gt2_any2d
! v311    16-11-21  subroutine webcanals_gt2_to_gt2_any2d remplacEe par 
!                   subroutine webcanals_gt1_to_gt2_any2d
!...............................................................................
!    _________                    .__                  .__             !m°v°w
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal ; use module_parallele ; use module_parameter
      use module_global
      implicit none
!     integer :: nbcanal=0
!     integer,dimension(:,:),allocatable :: canalcoord
!     integer,dimension(:,:,:),allocatable :: canalrank

      integer,parameter :: nexchgmax_mpi=2           &
                          ,tag_gt1gt2_u2d     =1616  &  
                          ,tag_gt1gt2_ssh     =1617  &   
                          ,tag_gt2gt1_ssh     =1618  &
                          ,tag_gt2gt1_ts      =1620  &
                          ,tag_gt1gt2_ts      =1621  &
                          ,tag_gt1gt2_vel     =1622  &
                          ,tag_gt1gt2_bio     =1623  &
                          ,tag_gt2gt1_bio     =1624  &
                          ,tag_gt2gt1_anytrc  =1625  &
                          ,tag_gt1gt2_anytrc  =1626  &   
                          ,tag_gt1gt2_anyvel  =1627  &
                          ,tag_gt1gt2_sshint  =1628  &
                          ,tag_gt2gt1_sshint  =1629  &
                          ,tag_gt2gt1_any2d   =1630  &
                          ,tag_gt1gt2_any2d   =1631    !16-11-21
      integer :: nexchg_=0,icnl_=0,point_=0
      integer,dimension(nexchgmax_mpi   ) :: tabreq_canal1_mpi
!     integer,dimension(mpi_status_size,nexchgmax_mpi   ) :: status_canal1_mpi
      integer,dimension(mpi_status_size) ::                      &
                                            status_gt1gt2_u2d    &  
                                           ,status_gt1gt2_ssh    &
                                           ,status_gt2gt1_ssh    & 
                                       !   ,status_gt1gt2_ssh2   &
                                           ,status_gt2gt1_ts     &
                                           ,status_gt1gt2_ts     &
                                           ,status_gt1gt2_vel    &
                                           ,status_gt1gt2_bio    &
                                           ,status_gt2gt1_bio    &
                                           ,status_gt1gt2_anytrc &
                                           ,status_gt2gt1_anytrc &  
                                           ,status_gt1gt2_anyvel &
                                           ,status_gt2gt1_sshint &
                                           ,status_gt1gt2_sshint &
                                           ,status_gt2gt1_any2d  &
                                           ,status_gt1gt2_any2d !16-11-21
      double precision,dimension(:),allocatable ::     &
                     ar_send_canal1,ar_recv_canal1     & 
                    ,ar_send_canal2,ar_recv_canal2      
      real ,dimension(:,:),allocatable ::              &
                     ar_send_canal3,ar_recv_canal3
                 
!     integer(kind=1),dimension(:),allocatable :: canalinout

contains

!......................................................................

      subroutine webcanals_gt1_to_gt2_velbar! (param_)
      implicit none
      integer(kind=1) :: sign_=1
!     integer :: di_send_=0 & ! decallage d'indice i pour send ssh
!               ,di_recv_=0 & ! decallage d'indice i pour receive ssh
!               ,dj_recv_=0   ! decallage d'indice j pour receive ssh
      integer :: loop_=0

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

        if(point_==1) then !point>
! point 1
         sign_=1
         isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)     !canal
         jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)     !canal
         irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)     !mer
         jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)     !mer
         if(canaldir(icnl_,point_)==1) then ; irecv=irecv+1 ; sign_=-1 ; endif ! decallage _u et _t convention grille C
         if(canaldir(icnl_,point_)==2) then ; jrecv=jrecv+1 ; sign_=-1 ; endif ! decallage _v et _t convention grille C

! Puisque la vitesse "limite" n'est pas envoyee par gridtype2, gridtype1 construit lui meme sa condition aux limite "gradient nul" !05-05-19
!        if(jsend>=0.and.jsend<=jmax+1) then !>>>
!        if(isend-1>=0.and.isend<=imax+2)velbar_u(isend-1,jsend,2)=velbar_u(isend,jsend,2)
!        endif                               !>>>

        else               !point>
! point 2
         sign_=-1
         isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)+1   !canal
         jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)     !canal
         irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)     !mer
         jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)     !mer
         if(canaldir(icnl_,point_)==1) then ; irecv=irecv+1 ; sign_=1 ; endif ! decallage _u et _t convention grille C
         if(canaldir(icnl_,point_)==2) then ; jrecv=jrecv+1 ; sign_=1 ; endif ! decallage _v et _t convention grille C

! Puisque la vitesse "limite" n'est pas envoyee par gridtype2, gridtype1 construit lui meme sa condition aux limite "gradient nul" !05-05-19
!        if(jsend>=0.and.jsend<=jmax+1) then !>>>
!        if(isend>=0.and.isend+1<=imax+2)velbar_u(isend+1,jsend,2)=velbar_u(isend,jsend,2)
!        endif                               !>>>

        endif              !points 

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          velbar_u(irecv,jrecv,2)=velbar_u(isend,jsend,2)*sign_ &
                                     *dy_u(isend,jsend)         &
                                     /dy_u(irecv,jrecv)
    
         else                                                            !>>>
          velbar_v(irecv,jrecv,2)=velbar_u(isend,jsend,2)*sign_ &
                                     *dy_u(isend,jsend)         &
                                     /dx_v(irecv,jrecv)
         endif                                                           !>>>

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal1(1)=velbar_u(isend,jsend,2)*sign_ &
                              *dy_u(isend,jsend)          

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_u2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_u2d,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_u2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_u2d,ierr)

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          velbar_u(irecv,jrecv,2)=ar_recv_canal1(1) &
             /dy_u(irecv,jrecv)
         else                                                            !>>>
          velbar_v(irecv,jrecv,2)=ar_recv_canal1(1) &
             /dx_v(irecv,jrecv)
         endif                                                           !>>>

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>


! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
      do loop_=1,3
      if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          velbar_u(irecv,jrecv,2)=velbar_u(isend,jsend,2)*sign_ &
                                     *dy_u(isend,jsend)         &
                                     /dy_u(irecv,jrecv)
    
         else                                                            !>>>
          velbar_v(irecv,jrecv,2)=velbar_u(isend,jsend,2)*sign_ &
                                     *dy_u(isend,jsend)         &
                                     /dx_v(irecv,jrecv)
         endif                                                           !>>>

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal1(1)=velbar_u(isend,jsend,2)*sign_ &
                              *dy_u(isend,jsend)          

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_u2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_u2d,ierr)

        else if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_u2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_u2d,ierr)

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          velbar_u(irecv,jrecv,2)=ar_recv_canal1(1) &
             /dy_u(irecv,jrecv)
         else                                                            !>>>
          velbar_v(irecv,jrecv,2)=ar_recv_canal1(1) &
             /dx_v(irecv,jrecv)
         endif                                                           !>>>

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

      endif                                 !m°v°m>
      enddo ! loop_


      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'coucou'


      end subroutine webcanals_gt1_to_gt2_velbar

!......................................................................

      subroutine webcanals_gt1_to_gt2_ssh! (param_)
      implicit none
      integer loop_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         ssh_w(irecv,jrecv,2)=ssh_w(isend,jsend,2)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal1(1)=ssh_w(isend,jsend,2)

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_ssh,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_ssh,ierr)

        else  if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_ssh,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_ssh,ierr)
         ssh_w(irecv,jrecv,2)=ar_recv_canal1(1)

        endif                                                           !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype1,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

            ssh_w(irecv,jrecv,2)=ssh_w(isend,jsend,2)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

            ar_send_canal1(1)=ssh_w(isend,jsend,2)

            call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_ssh,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_ssh,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_ssh,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_ssh,ierr)
            ssh_w(irecv,jrecv,2)=ar_recv_canal1(1)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_



      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt1_to_gt2_ssh

!......................................................................

      subroutine webcanals_gt2_to_gt1_ssh! (param_)
      implicit none
      integer loop_

! Le receiver est le canal: gridtype=gridtype1
! le sender est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype1,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

         isend=i_canalcoord(icnl_,point_,gridtype2,sender)  -par%timax(1)     !gt2
         jsend=j_canalcoord(icnl_,point_,gridtype2,sender)  -par%tjmax(1)     !gt2
         irecv=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1)     !gt1 (canal)
         jrecv=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1)     !gt1 (canal)

       if(par%rank==canalrank(icnl_,point_,gridtype2,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype1,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         ssh_w(irecv,jrecv,2)=ssh_w(isend,jsend,2)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         ar_send_canal1(1)=ssh_w(isend,jsend,2)

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,receiver),tag_gt2gt1_ssh,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_ssh,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype1,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_ssh,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_ssh,ierr)

         ssh_w(irecv,jrecv,2)=ar_recv_canal1(1)

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0 !04-05-19
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype1,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype2,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype1,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

            ssh_w(irecv,jrecv,2)=ssh_w(isend,jsend,2)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

            ar_send_canal1(1)=ssh_w(isend,jsend,2)

            call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype1,loop_),tag_gt2gt1_ssh,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_ssh,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype1,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_ssh,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_ssh,ierr)
            ssh_w(irecv,jrecv,2)=ar_recv_canal1(1)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt2_to_gt1_ssh

!......................................................................

      subroutine webcanals_gt2_to_gt1_temsal! (param_)
      implicit none
      integer loop_

! Le receiver est le canal: gridtype=gridtype1
! le sender est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal2)) then !>>>
       allocate(ar_send_canal2(kmax*2)) ; ar_send_canal2=0.
       allocate(ar_recv_canal2(kmax*2)) ; ar_recv_canal2=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype1,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

         isend=i_canalcoord(icnl_,point_,gridtype2,sender)  -par%timax(1)     !gt2
         jsend=j_canalcoord(icnl_,point_,gridtype2,sender)  -par%tjmax(1)     !gt2
         irecv=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1)     !gt1 (canal)
         jrecv=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1)     !gt1 (canal)

       if(par%rank==canalrank(icnl_,point_,gridtype2,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype1,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             tem_t(irecv,jrecv,1:kmax,2)=    tem_t(isend,jsend,1:kmax,2)
             sal_t(irecv,jrecv,1:kmax,2)=    sal_t(isend,jsend,1:kmax,2)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)       =tem_t(isend,jsend,1:kmax,2)
         ar_send_canal2(kmax+1:kmax*2)=sal_t(isend,jsend,1:kmax,2)

         call mpi_issend(ar_send_canal2,kmax*2,mpi_double,canalrank(icnl_,point_,gridtype1,receiver),tag_gt2gt1_ts,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_ts,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype1,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax*2,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_ts,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_ts,ierr)

             tem_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax)
             sal_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(kmax+1:kmax*2)


        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0 !04-05-19
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype1,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype2,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype1,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             tem_t(irecv,jrecv,1:kmax,2)=    tem_t(isend,jsend,1:kmax,2)
             sal_t(irecv,jrecv,1:kmax,2)=    sal_t(isend,jsend,1:kmax,2)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

            ar_send_canal2(1:kmax)       =tem_t(isend,jsend,1:kmax,2)
            ar_send_canal2(kmax+1:kmax*2)=sal_t(isend,jsend,1:kmax,2)

            call mpi_issend(ar_send_canal2,kmax*2,mpi_double,canalrankbis(icnl_,point_,gridtype1,loop_),tag_gt2gt1_ts,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_ts,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype1,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal2,kmax*2,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_ts,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_ts,ierr)

             tem_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax)
             sal_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(kmax+1:kmax*2)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt2_to_gt1_temsal

!......................................................................

      subroutine webcanals_gt1_to_gt2_temsal! (param_)
      implicit none
      integer loop_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal2)) then !>>>
       allocate(ar_send_canal2(kmax*2)) ; ar_send_canal2=0.
       allocate(ar_recv_canal2(kmax*2)) ; ar_recv_canal2=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             tem_t(irecv,jrecv,1:kmax,2)=    tem_t(isend,jsend,1:kmax,2)
             sal_t(irecv,jrecv,1:kmax,2)=    sal_t(isend,jsend,1:kmax,2)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)       =tem_t(isend,jsend,1:kmax,2)
         ar_send_canal2(kmax+1:kmax*2)=sal_t(isend,jsend,1:kmax,2)

         call mpi_issend(ar_send_canal2,kmax*2,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_ts,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_ts,ierr)

        else  if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax*2,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_ts,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_ts,ierr)

             tem_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax)
             sal_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(kmax+1:kmax*2)

        endif                                                           !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype1,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             tem_t(irecv,jrecv,1:kmax,2)=    tem_t(isend,jsend,1:kmax,2)
             sal_t(irecv,jrecv,1:kmax,2)=    sal_t(isend,jsend,1:kmax,2)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

            ar_send_canal2(1:kmax)       =tem_t(isend,jsend,1:kmax,2)
            ar_send_canal2(kmax+1:kmax*2)=sal_t(isend,jsend,1:kmax,2)

            call mpi_issend(ar_send_canal2,kmax*2,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_ts,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_ts,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal2,kmax*2,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_ts,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_ts,ierr)

             tem_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax)
             sal_t(irecv,jrecv,1:kmax,2)=ar_recv_canal2(kmax+1:kmax*2)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_



      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt1_to_gt2_temsal

!......................................................................

      subroutine webcanals_gt1_to_gt2_vel! (param_)
      implicit none
      integer(kind=1) :: sign_=1
!     integer :: di_send_=0 & ! decallage d'indice i pour send ssh
!               ,di_recv_=0 & ! decallage d'indice i pour receive ssh
!               ,dj_recv_=0   ! decallage d'indice j pour receive ssh
      integer :: loop_=0

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal2)) then !>>>
! On alloue "au plus grand" A savoir pour l'envoie tem, sal soit 2*kmax
! mais pour les vitesses on n'envoie/recoit que kmax valeurs
       allocate(ar_send_canal2(2*kmax)) ; ar_send_canal2=0.
       allocate(ar_recv_canal2(2*kmax)) ; ar_recv_canal2=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

        if(point_==1) then !point>
! point 1
         sign_=1
         isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)     !canal
         jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)     !canal
         irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)     !mer
         jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)     !mer
         if(canaldir(icnl_,point_)==1) then ; irecv=irecv+1 ; sign_=-1 ; endif ! decallage _u et _t convention grille C
         if(canaldir(icnl_,point_)==2) then ; jrecv=jrecv+1 ; sign_=-1 ; endif ! decallage _v et _t convention grille C
! Puisque la vitesse "limite" n'est pas envoyee par gridtype2, gridtype1 construit lui meme sa condition aux limite "gradient nul" !05-05-19
!        if(jsend>=0.and.jsend<=jmax+1) then !>>>
!        if(isend-1>=0.and.isend<=imax+2)vel_u(isend-1,jsend,1:kmax,2)=vel_u(isend,jsend,1:kmax,2)
!        endif                               !>>>
        else               !point>
! point 2
         sign_=-1
         isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)+1   !canal
         jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)     !canal
         irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)     !mer
         jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)     !mer
         if(canaldir(icnl_,point_)==1) then ; irecv=irecv+1 ; sign_=1 ; endif ! decallage _u et _t convention grille C
         if(canaldir(icnl_,point_)==2) then ; jrecv=jrecv+1 ; sign_=1 ; endif ! decallage _v et _t convention grille C
! Puisque la vitesse "limite" n'est pas envoyee par gridtype2, gridtype1 construit lui meme sa condition aux limite "gradient nul" !05-05-19
!        if(jsend>=0.and.jsend<=jmax+1) then !>>>
!        if(isend>=0.and.isend+1<=imax+2)vel_u(isend+1,jsend,1:kmax,2)=vel_u(isend,jsend,1:kmax,2)
!        endif                               !>>>
        endif              !points 

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          vel_u(irecv,jrecv,1:kmax,2)=vel_u(isend,jsend,1:kmax,2)*sign_ &
                                      *dy_u(isend,jsend)         &
                                      /dy_u(irecv,jrecv)
    
         else                                                            !>>>
          vel_v(irecv,jrecv,1:kmax,2)=vel_u(isend,jsend,1:kmax,2)*sign_ &
                                      *dy_u(isend,jsend)         &
                                      /dx_v(irecv,jrecv)
         endif                                                           !>>>

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)=vel_u(isend,jsend,1:kmax,2)*sign_ &
                                *dy_u(isend,jsend)          

         call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_vel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_vel,ierr)


        else if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_vel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_vel,ierr)


         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          vel_u(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax) &
          /dy_u(irecv,jrecv)
         else                                                            !>>>
          vel_v(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax) &
          /dx_v(irecv,jrecv)
         endif                                                           !>>>

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>


! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
      do loop_=1,3
      if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          vel_u(irecv,jrecv,1:kmax,2)=vel_u(isend,jsend,1:kmax,2)*sign_ &
                                      *dy_u(isend,jsend)         &
                                      /dy_u(irecv,jrecv)
    
         else                                                            !>>>
          vel_v(irecv,jrecv,1:kmax,2)=vel_u(isend,jsend,1:kmax,2)*sign_ &
                                      *dy_u(isend,jsend)         &
                                      /dx_v(irecv,jrecv)
         endif                                                           !>>>

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)=vel_u(isend,jsend,1:kmax,2)*sign_ &
                                *dy_u(isend,jsend)          

         call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_vel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_vel,ierr)

        else if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_vel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_vel,ierr)

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          vel_u(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax) &
          /dy_u(irecv,jrecv)
         else                                                            !>>>
          vel_v(irecv,jrecv,1:kmax,2)=ar_recv_canal2(1:kmax) &
          /dx_v(irecv,jrecv)
         endif                                                           !>>>

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

      endif                                 !m°v°m>
      enddo ! loop_


      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'coucou'


      end subroutine webcanals_gt1_to_gt2_vel

!......................................................................

      subroutine webcanals_gt2_to_gt1_bio! (param_)
      implicit none
      integer loop_

! Le receiver est le canal: gridtype=gridtype1
! le sender est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal3)) then !>>>
       allocate(ar_send_canal3(kmax,vbmax)) ; ar_send_canal3=0.
       allocate(ar_recv_canal3(kmax,vbmax)) ; ar_recv_canal3=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype1,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

         isend=i_canalcoord(icnl_,point_,gridtype2,sender)  -par%timax(1)     !gt2
         jsend=j_canalcoord(icnl_,point_,gridtype2,sender)  -par%tjmax(1)     !gt2
         irecv=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1)     !gt1 (canal)
         jrecv=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1)     !gt1 (canal)

       if(par%rank==canalrank(icnl_,point_,gridtype2,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype1,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=    bio_t(isend,jsend,1:kmax,1:vbmax)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         ar_send_canal3(1:kmax,1:vbmax)       =bio_t(isend,jsend,1:kmax,1:vbmax)

         call mpi_issend(ar_send_canal3,kmax*vbmax,mpi_real,canalrank(icnl_,point_,gridtype1,receiver),tag_gt2gt1_bio,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_bio,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype1,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal3,kmax*vbmax,mpi_real,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_bio,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_bio,ierr)

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=ar_recv_canal3(1:kmax,1:vbmax)


        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0 !04-05-19
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype1,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype2,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype1,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=    bio_t(isend,jsend,1:kmax,1:vbmax)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

            ar_send_canal3(1:kmax,1:vbmax)       =bio_t(isend,jsend,1:kmax,1:vbmax)

            call mpi_issend(ar_send_canal3,kmax*vbmax,mpi_real,canalrankbis(icnl_,point_,gridtype1,loop_),tag_gt2gt1_bio,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_bio,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype1,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal3,kmax*vbmax,mpi_real,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_bio,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_bio,ierr)

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=ar_recv_canal3(1:kmax,1:vbmax)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt2_to_gt1_bio

!......................................................................

      subroutine webcanals_gt1_to_gt2_bio! (param_)
      implicit none
      integer loop_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal3)) then !>>>
       allocate(ar_send_canal3(kmax,vbmax)) ; ar_send_canal3=0.
       allocate(ar_recv_canal3(kmax,vbmax)) ; ar_recv_canal3=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=    bio_t(isend,jsend,1:kmax,1:vbmax)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal3(1:kmax,1:vbmax)       =bio_t(isend,jsend,1:kmax,1:vbmax)

         call mpi_issend(ar_send_canal3,kmax*vbmax,mpi_real,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_bio,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_bio,ierr)


        else  if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal3,kmax*vbmax,mpi_real,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_bio,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_bio,ierr)

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=ar_recv_canal3(1:kmax,1:vbmax)

        endif                                                           !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype1,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=    bio_t(isend,jsend,1:kmax,1:vbmax)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

            ar_send_canal3(1:kmax,1:vbmax)       =bio_t(isend,jsend,1:kmax,1:vbmax)

            call mpi_issend(ar_send_canal3,kmax*vbmax,mpi_real,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_bio,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_bio,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal3,kmax*vbmax,mpi_real,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_bio,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_bio,ierr)

             bio_t(irecv,jrecv,1:kmax,1:vbmax)=ar_recv_canal3(1:kmax,1:vbmax)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_



      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt1_to_gt2_bio

!......................................................................

      subroutine webcanals_gt2_to_gt1_anytrc(id_var_)
      implicit none
      integer id_var_,loop_

! Le receiver est le canal: gridtype=gridtype1
! le sender est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal2)) then !>>>
! On alloue 2 fois kmax A cause de l'envoi tem, sal mais ici on n'envoie que kmax valeurs
       allocate(ar_send_canal2(kmax*2)) ; ar_send_canal2=0.
       allocate(ar_recv_canal2(kmax*2)) ; ar_recv_canal2=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype1,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

         isend=i_canalcoord(icnl_,point_,gridtype2,sender)  -par%timax(1)     !gt2
         jsend=j_canalcoord(icnl_,point_,gridtype2,sender)  -par%tjmax(1)     !gt2
         irecv=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1)     !gt1 (canal)
         jrecv=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1)     !gt1 (canal)

       if(par%rank==canalrank(icnl_,point_,gridtype2,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype1,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             anyv3d(irecv,jrecv,1:kmax,id_var_)=    anyv3d(isend,jsend,1:kmax,id_var_)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)       =anyv3d(isend,jsend,1:kmax,id_var_)

         call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,receiver),tag_gt2gt1_anytrc,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_anytrc,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype1,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_anytrc,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_anytrc,ierr)

             anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax)


        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0 !04-05-19
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype1,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype2,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype1,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             anyv3d(irecv,jrecv,1:kmax,id_var_)=    anyv3d(isend,jsend,1:kmax,id_var_)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

            ar_send_canal2(1:kmax)       =anyv3d(isend,jsend,1:kmax,id_var_)

            call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrankbis(icnl_,point_,gridtype1,loop_),tag_gt2gt1_anytrc,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_anytrc,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype1,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_anytrc,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_anytrc,ierr)

             anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt2_to_gt1_anytrc

!......................................................................

      subroutine webcanals_gt1_to_gt2_anytrc(id_var_)
      implicit none
      integer loop_,id_var_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal2)) then !>>>
       allocate(ar_send_canal2(2*kmax)) ; ar_send_canal2=0.
       allocate(ar_recv_canal2(2*kmax)) ; ar_recv_canal2=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             anyv3d(irecv,jrecv,1:kmax,id_var_)=    anyv3d(isend,jsend,1:kmax,id_var_)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)       =anyv3d(isend,jsend,1:kmax,id_var_)

         call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_anytrc,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_anytrc,ierr)

        else  if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_anytrc,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_anytrc,ierr)

             anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax)

        endif                                                           !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype1,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

             anyv3d(irecv,jrecv,1:kmax,id_var_)=    anyv3d(isend,jsend,1:kmax,id_var_)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

            ar_send_canal2(1:kmax)       =anyv3d(isend,jsend,1:kmax,id_var_)

            call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_anytrc,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_anytrc,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_anytrc,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_anytrc,ierr)

             anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_



      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt1_to_gt2_anytrc

!......................................................................

      subroutine webcanals_gt1_to_gt2_anyvel(id_var_)
      implicit none
      integer(kind=1) :: sign_=1
!     integer :: di_send_=0 & ! decallage d'indice i pour send ssh
!               ,di_recv_=0 & ! decallage d'indice i pour receive ssh
!               ,dj_recv_=0   ! decallage d'indice j pour receive ssh
      integer :: loop_=0,id_var_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal2)) then !>>>
! On alloue "au plus grand" A savoir pour l'envoie tem, sal, ssh soit 1+2*kmax
! mais pour les vitesses on n'envoie/recoit que kmax valeurs
       allocate(ar_send_canal2(2*kmax)) ; ar_send_canal2=0.
       allocate(ar_recv_canal2(2*kmax)) ; ar_recv_canal2=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

        if(point_==1) then !point>
! point 1
         sign_=1
         isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)     !canal
         jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)     !canal
         irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)     !mer
         jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)     !mer
         if(canaldir(icnl_,point_)==1) then ; irecv=irecv+1 ; sign_=-1 ; endif ! decallage _u et _t convention grille C
         if(canaldir(icnl_,point_)==2) then ; jrecv=jrecv+1 ; sign_=-1 ; endif ! decallage _v et _t convention grille C
        else               !point>
! point 2
         sign_=-1
         isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)+1   !canal
         jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)     !canal
         irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)     !mer
         jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)     !mer
         if(canaldir(icnl_,point_)==1) then ; irecv=irecv+1 ; sign_=1 ; endif ! decallage _u et _t convention grille C
         if(canaldir(icnl_,point_)==2) then ; jrecv=jrecv+1 ; sign_=1 ; endif ! decallage _v et _t convention grille C
        endif              !points 

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=anyv3d(isend,jsend,1:kmax,id_var_)*sign_ &
                                              *dy_u(isend,jsend)         &
                                              /dy_u(irecv,jrecv)
    
         else                                                            !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=anyv3d(isend,jsend,1:kmax,id_var_)*sign_ &
                                              *dy_u(isend,jsend)         &
                                              /dx_v(irecv,jrecv)
         endif                                                           !>>>

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)=anyv3d(isend,jsend,1:kmax,id_var_)*sign_ &
                                 *dy_u(isend,jsend)          

         call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_anyvel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_anyvel,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_anyvel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_anyvel,ierr)

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax) &
           /dy_u(irecv,jrecv)
         else                                                            !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax) &
           /dx_v(irecv,jrecv)
         endif                                                           !>>>

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>


! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
      do loop_=1,3
      if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=anyv3d(isend,jsend,1:kmax,id_var_)*sign_ &
                                              *dy_u(isend,jsend)         &
                                              /dy_u(irecv,jrecv)
    
         else                                                            !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=anyv3d(isend,jsend,1:kmax,id_var_)*sign_ &
                                              *dy_u(isend,jsend)         &
                                              /dx_v(irecv,jrecv)
         endif                                                           !>>>

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal2(1:kmax)=anyv3d(isend,jsend,1:kmax,id_var_)*sign_ &
                                 *dy_u(isend,jsend)          

         call mpi_issend(ar_send_canal2,kmax,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_anyvel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_anyvel,ierr)

        else if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal2,kmax,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_anyvel,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_anyvel,ierr)

         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax) &
           /dy_u(irecv,jrecv)
         else                                                            !>>>
          anyv3d(irecv,jrecv,1:kmax,id_var_)=ar_recv_canal2(1:kmax) &
           /dx_v(irecv,jrecv)
         endif                                                           !>>>

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

      endif                                 !m°v°m>
      enddo ! loop_


      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'coucou'


      end subroutine webcanals_gt1_to_gt2_anyvel

!......................................................................

      subroutine webcanals_gt1_to_gt2_sshint! (param_)
      implicit none
      integer loop_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         ssh_int_w(irecv,jrecv,2)=ssh_int_w(isend,jsend,2)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         ar_send_canal1(1)=ssh_int_w(isend,jsend,2)

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_sshint,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_sshint,ierr)

        else  if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_sshint,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_sshint,ierr)
         ssh_int_w(irecv,jrecv,2)=ar_recv_canal1(1)

        endif                                                           !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype1,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

            ssh_int_w(irecv,jrecv,2)=ssh_int_w(isend,jsend,2)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

            ar_send_canal1(1)=ssh_int_w(isend,jsend,2)

            call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_sshint,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_sshint,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_sshint,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_sshint,ierr)
            ssh_int_w(irecv,jrecv,2)=ar_recv_canal1(1)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_



      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt1_to_gt2_sshint

!......................................................................

      subroutine webcanals_gt2_to_gt1_sshint! (param_)
      implicit none
      integer loop_

! Le receiver est le canal: gridtype=gridtype1
! le sender est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype1,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

         isend=i_canalcoord(icnl_,point_,gridtype2,sender)  -par%timax(1)     !gt2
         jsend=j_canalcoord(icnl_,point_,gridtype2,sender)  -par%tjmax(1)     !gt2
         irecv=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1)     !gt1 (canal)
         jrecv=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1)     !gt1 (canal)

       if(par%rank==canalrank(icnl_,point_,gridtype2,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype1,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         ssh_int_w(irecv,jrecv,2)=ssh_int_w(isend,jsend,2)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         ar_send_canal1(1)=ssh_int_w(isend,jsend,2)

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,receiver),tag_gt2gt1_sshint,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_sshint,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype1,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_sshint,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_sshint,ierr)

         ssh_int_w(irecv,jrecv,2)=ar_recv_canal1(1)

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0 !04-05-19
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype1,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype2,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype1,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

            ssh_int_w(irecv,jrecv,2)=ssh_int_w(isend,jsend,2)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

            ar_send_canal1(1)=ssh_int_w(isend,jsend,2)

            call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype1,loop_),tag_gt2gt1_sshint,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_sshint,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype1,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_sshint,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_sshint,ierr)
            ssh_int_w(irecv,jrecv,2)=ar_recv_canal1(1)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt2_to_gt1_sshint
!......................................................................

      subroutine webcanals_messages !02-05-19
      implicit none
      integer :: point_=0 , icnl_=1

! A toutes fins utiles on ecrit les valeurs de h,dx,dy,lon,lat
! pour eventuellement les reporter dans les fichiers
! de canaux

       do icnl_=1,nbcanal

        j=j_canalcoord(icnl_,point1,gridtype1,2)-par%tjmax(1)
        do i1=i_canalcoord(icnl_,point1,gridtype1,2),i_canalcoord(icnl_,point2,gridtype1,2)
        i=i1-par%timax(1)
        if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !m°v°m>

        write(texte30,'(a,i0,a,i0)')'tmp/messages_channel',icnl_,'_rank',par%rank
        open(3,file=trim(texte30),position='append')
        write(3,'(i5,1x,i5,3(1x,f11.4),1x,f12.7,1x,f12.7)')    &
        i+par%timax(1),j+par%tjmax(1),h_w(i,j)                 &
                                    ,dx_t(i,j),dy_t(i,j)       &
                                   ,lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg 
        close(3)

        endif                                      !m°v°m>

        enddo ! i1
       enddo  ! icnl_
     
      end subroutine webcanals_messages

!......................................................................

      subroutine webcanals_gt2_to_gt1_any2d(varname_)
      implicit none
      character(len=*)varname_
      integer loop_

! Le receiver est le canal: gridtype=gridtype1
! le sender est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype1,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148

         isend=i_canalcoord(icnl_,point_,gridtype2,sender)  -par%timax(1)     !gt2
         jsend=j_canalcoord(icnl_,point_,gridtype2,sender)  -par%tjmax(1)     !gt2
         irecv=i_canalcoord(icnl_,point_,gridtype1,receiver)-par%timax(1)     !gt1 (canal)
         jrecv=j_canalcoord(icnl_,point_,gridtype1,receiver)-par%tjmax(1)     !gt1 (canal)

       if(par%rank==canalrank(icnl_,point_,gridtype2,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype1,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =lon_t(isend,jsend)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dx_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dy_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)= xy_t(isend,jsend,1)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         if(varname_=='lon_t')ar_send_canal1(1)=lon_t(isend,jsend)
         if(varname_=='lat_t')ar_send_canal1(1)=lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t') ar_send_canal1(1)= dx_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t') ar_send_canal1(1)= dy_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t') ar_send_canal1(1)= xy_t(isend,jsend,1)

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,receiver),tag_gt2gt1_any2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_any2d,ierr)

        else if (par%rank==canalrank(icnl_,point_,gridtype1,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_any2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt2gt1_any2d,ierr)

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dx_t')  dx_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dy_t')  dy_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)=ar_recv_canal1(1)

        endif                                                      !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype1,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype2,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype1,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =lon_t(isend,jsend)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dx_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dy_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)= xy_t(isend,jsend,1)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype2,sender))  then         !MPI-MPI>

         if(varname_=='lon_t')ar_send_canal1(1)=lon_t(isend,jsend)
         if(varname_=='lat_t')ar_send_canal1(1)=lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t') ar_send_canal1(1)= dx_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t') ar_send_canal1(1)= dy_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t') ar_send_canal1(1)= xy_t(isend,jsend,1)

            call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype1,loop_),tag_gt2gt1_any2d,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_any2d,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype1,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,sender),tag_gt2gt1_any2d,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt2gt1_any2d,ierr)

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dx_t')  dx_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dy_t')  dy_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)=ar_recv_canal1(1)


           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_

      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt2_to_gt1_any2d

!......................................................................

      subroutine webcanals_gt1_to_gt2_any2d(varname_) !23-06-19!16-11-21
      implicit none
      character(len=*)varname_
      integer loop_

! Le sender est le canal: gridtype=gridtype1
! le receiver est la mer: gridtype=gridtype2

                 
      gridtype1=1
      gridtype2=2
      sender=1  
      receiver=2 

      if(.not.allocated(ar_send_canal1)) then !>>>
       allocate(ar_send_canal1(1)) ; ar_send_canal1=0.
       allocate(ar_recv_canal1(1)) ; ar_recv_canal1=0.
      endif                                   !>>>

      do icnl_=1,nbcanal  

      do point_=1,2
! Une des 2 extremites peut ne pas etre connectee a un point de mer 
! d'oU le test if suivant qui est la mer soit gridtype=2
      if(canalrank(icnl_,point_,gridtype2,receiver)>=0) then !test de connexion a un de point de mer>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype1,sender)  -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype1,sender)  -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype2,receiver)-par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype2,receiver)-par%tjmax(1)

       if(par%rank==canalrank(icnl_,point_,gridtype1,sender).and.&
          par%rank==canalrank(icnl_,point_,gridtype2,receiver)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =lon_t(isend,jsend)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dx_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dy_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)= xy_t(isend,jsend,1)

       else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

        if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         if(varname_=='lon_t')ar_send_canal1(1)=lon_t(isend,jsend)
         if(varname_=='lat_t')ar_send_canal1(1)=lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t') ar_send_canal1(1)= dx_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t') ar_send_canal1(1)= dy_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t') ar_send_canal1(1)= xy_t(isend,jsend,1)

         call mpi_issend(ar_send_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype2,receiver),tag_gt1gt2_any2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_any2d,ierr)

        else  if (par%rank==canalrank(icnl_,point_,gridtype2,receiver))  then !MPI-MPI>

         call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_any2d,par%comm2d,k0,ierr)
         call mpi_wait(k0,status_gt1gt2_any2d,ierr)

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dx_t')  dx_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dy_t')  dy_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)=ar_recv_canal1(1)

        endif                                                           !MPI-MPI>

       endif                                               !AAAAA>

! Ranks supplementaire si irec,jrecv est dans la zone d'ombre mpi z0
        do loop_=1,3 
         if(canalrankbis(icnl_,point_,gridtype2,loop_)>=0) then !m°v°m>

          if(par%rank==canalrank   (icnl_,point_,gridtype1,sender).and.&
             par%rank==canalrankbis(icnl_,point_,gridtype2,loop_)) then !AAAAA>
! sender et receiver ont un RANK IDENTIQUE = PAS DE MPI

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =lon_t(isend,jsend)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dx_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t')  dx_t(irecv,jrecv)  = dy_t(isend,jsend)
          if(varname_=='dy_t')  dy_t(irecv,jrecv)  = dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)= xy_t(isend,jsend,1)

          else                                                !AAAAA>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           if(par%rank==canalrank(icnl_,point_,gridtype1,sender))  then         !MPI-MPI>

         if(varname_=='lon_t')ar_send_canal1(1)=lon_t(isend,jsend)
         if(varname_=='lat_t')ar_send_canal1(1)=lat_t(isend,jsend)
         if(canaldir(icnl_,point_)==1.or.canaldir(icnl_,point_)==3) then !ooo> !16-11-21>
          if(varname_=='dx_t') ar_send_canal1(1)= dx_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dy_t(isend,jsend)
         else                                                            !ooo> !16-11-21>
! Si direction differente le sender envoie dy a la place de dx, ou dx a la place de dy
          if(varname_=='dx_t') ar_send_canal1(1)= dy_t(isend,jsend)
          if(varname_=='dy_t') ar_send_canal1(1)= dx_t(isend,jsend)
         endif                                                           !ooo> !16-11-21>
         if(varname_=='xy_t') ar_send_canal1(1)= xy_t(isend,jsend,1)

            call mpi_issend(ar_send_canal1,1,mpi_double,canalrankbis(icnl_,point_,gridtype2,loop_),tag_gt1gt2_any2d,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_any2d,ierr)

           else  if (par%rank==canalrankbis(icnl_,point_,gridtype2,loop_))  then !MPI-MPI>

            call mpi_irecv(ar_recv_canal1,1,mpi_double,canalrank(icnl_,point_,gridtype1,sender),tag_gt1gt2_any2d,par%comm2d,k0,ierr)
            call mpi_wait(k0,status_gt1gt2_any2d,ierr)

         if(varname_=='lon_t')lon_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='lat_t')lat_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dx_t')  dx_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='dy_t')  dy_t(irecv,jrecv)  =ar_recv_canal1(1)
         if(varname_=='xy_t')  xy_t(irecv,jrecv,1)=ar_recv_canal1(1)

           endif                                                           !MPI-MPI>

          endif                                               !AAAAA>
    
         endif                                                  !m°v°m>
        enddo ! loop_



      endif                                 !test de connexion a un de point de mer>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  


      end subroutine webcanals_gt1_to_gt2_any2d

!......................................................................

      end module module_webcanals
