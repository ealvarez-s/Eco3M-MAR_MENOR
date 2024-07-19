      module module_ruban_mpi
!______________________________________________________________________
! S model
! release S26 - last update: 13-10-15
!______________________________________________________________________
! Version  Date      Description des modifications
! S26      20-01-17  Mise en service
!______________________________________________________________________
      use module_principal
      use module_parallele
      implicit none

      integer,parameter :: tagouest_mpi     =5607  &
                          ,tagest_mpi       =5617  &
                          ,tagsud_mpi       =6607  &
                          ,tagnord_mpi      =6617  &
                          ,tagsudouest_mpi =15607  &
                          ,tagsudest_mpi   =15617  &
                          ,tagnordouest_mpi=16607  &
                          ,tagnordest_mpi  =16617  &
                          ,nexchgmax_mpi   =16

      integer,dimension(nexchgmax_mpi   ) :: tabreq_mpi
      integer,dimension(mpi_status_size,nexchgmax_mpi   ) :: tstatus_mpi

      integer   nordest_send_i(5) &
             ,nordouest_send_i(5) &
                ,sudest_send_i(5) &
              ,sudouest_send_i(5) &
               ,nordest_send_j(5) &
             ,nordouest_send_j(5) &
                ,sudest_send_j(5) &
              ,sudouest_send_j(5)

      integer   nordest_recv_i(5) &
             ,nordouest_recv_i(5) &
                ,sudest_recv_i(5) &
              ,sudouest_recv_i(5) &
               ,nordest_recv_j(5) &
             ,nordouest_recv_j(5) &
                ,sudest_recv_j(5) &
              ,sudouest_recv_j(5)

      integer                   &
              nb_send_nordest   &
             ,nb_send_nordouest &
             ,nb_send_sudest    &
             ,nb_send_sudouest  &
             ,nb_send_est       &
             ,nb_send_ouest     &
             ,nb_send_nord      &
             ,nb_send_sud       &
             ,nb_recv_nordest   &
             ,nb_recv_nordouest &
             ,nb_recv_sudest    &
             ,nb_recv_sudouest  &
             ,nb_recv_est       &
             ,nb_recv_ouest     &
             ,nb_recv_nord      &
             ,nb_recv_sud        

      double precision , dimension(:) , allocatable :: &
                                    ar_send_nordest    & 
                                   ,ar_send_nordouest  & 
                                   ,ar_send_sudest     & 
                                   ,ar_send_sudouest   & 
                                   ,ar_send_est        & 
                                   ,ar_send_ouest      & 
                                   ,ar_send_nord       & 
                                   ,ar_send_sud        &
                                   ,ar_recv_nordest    & 
                                   ,ar_recv_nordouest  & 
                                   ,ar_recv_sudest     & 
                                   ,ar_recv_sudouest   & 
                                   ,ar_recv_est        & 
                                   ,ar_recv_ouest      & 
                                   ,ar_recv_nord       & 
                                   ,ar_recv_sud    

      integer istr_mpi,iend_mpi,jstr_mpi,jend_mpi

      contains

!....................................................................

      subroutine ruban_mpi_allocate
      implicit none

       allocate(ar_send_nordest  (kmax*5))
       allocate(ar_send_nordouest(kmax*5))
       allocate(ar_send_sudest   (kmax*5)) 
       allocate(ar_send_sudouest (kmax*5))
       allocate(ar_send_est      (kmax*(jmax  )))
       allocate(ar_send_ouest    (kmax*(jmax  )))
       allocate(ar_send_nord     (kmax*(imax  )))
       allocate(ar_send_sud      (kmax*(imax  )))

       allocate(ar_recv_nordest  (kmax*5))
       allocate(ar_recv_nordouest(kmax*5))
       allocate(ar_recv_sudest   (kmax*5)) 
       allocate(ar_recv_sudouest (kmax*5))
       allocate(ar_recv_est      (kmax*(jmax  )))
       allocate(ar_recv_ouest    (kmax*(jmax  )))
       allocate(ar_recv_nord     (kmax*(imax  )))
       allocate(ar_recv_sud      (kmax*(imax  )))

      end subroutine ruban_mpi_allocate

!....................................................................

      subroutine ruban_mpi_beforesendrecv
      implicit none
      integer nexchg_  

!.............................................................................................
! Les domaines indiquent à leurs voisins le nombre d'element s'appretant à traverser les
! frontieres. 

      nexchg_   =0
! Coin Nord-Est:
      if (par%tvoisin(nordest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_nordest,1,mpi_integer,par%tvoisin(nordest)  &
              ,tagnordest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_nordest,1,mpi_integer,par%tvoisin(nordest)  &
              ,tagsudouest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
!     endif
      endif
! Coin Nord-Ouest:
      if (par%tvoisin(nordouest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_nordouest,1,mpi_integer,par%tvoisin(nordouest)  &
              ,tagnordouest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_nordouest,1,mpi_integer,par%tvoisin(nordouest)  &
              ,tagsudest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif
! Coin Sud-Ouest:
      if (par%tvoisin(sudouest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_sudouest,1,mpi_integer,par%tvoisin(sudouest)  &
              ,tagsudouest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_sudouest,1,mpi_integer,par%tvoisin(sudouest)  &
              ,tagnordest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif
! Coin Sud-Est:
      if (par%tvoisin(sudest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_sudest,1,mpi_integer,par%tvoisin(sudest)  &
              ,tagsudest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_sudest,1,mpi_integer,par%tvoisin(sudest)  &
              ,tagnordouest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_est,1,mpi_integer,par%tvoisin(est)  &
              ,tagest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_est,1,mpi_integer,par%tvoisin(est)  &
              ,tagouest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif
! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_ouest,1,mpi_integer,par%tvoisin(ouest)  &
              ,tagouest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_ouest,1,mpi_integer,par%tvoisin(ouest)  &
              ,tagest_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif
! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_nord,1,mpi_integer,par%tvoisin(nord)  &
              ,tagnord_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_nord,1,mpi_integer,par%tvoisin(nord)  &
              ,tagsud_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif
! Frontiere sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nb_send_sud,1,mpi_integer,par%tvoisin(sud)  &
              ,tagsud_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nb_recv_sud,1,mpi_integer,par%tvoisin(sud)  &
              ,tagnord_mpi   ,par%comm2d,tabreq_mpi   (nexchg_   ),ierr)
      endif

      if(nexchg_   >nexchgmax_mpi   ) &
      stop 'drifter_who_is_out nexchg_   >nexchgmax_mpi   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                  &
                      ,tabreq_mpi   (1:nexchg_   ) &
                      ,tstatus_mpi  (:,1:nexchg_ ) &
                      ,ierr)

      end subroutine ruban_mpi_beforesendrecv

!....................................................................

      subroutine ruban_mpi_sendrecv
      implicit none
      integer nexchg_  

! Preparer l'Echange:
      call ruban_mpi_beforesendrecv

!.............................................................................................
! ECHANGER LES RUBANS 

      nexchg_   =0

! Coin Nord-Est:
      if (par%tvoisin(nordest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_nordest)
      call mpi_issend(ar_send_nordest(1)    & ! envoyé au proc Est
                ,size(ar_send_nordest(1:k)) &
                ,mpi_double_precision                        &
                ,par%tvoisin(nordest)            &
                ,tagnordest_mpi                     &
                ,par%comm2d                      &
                ,tabreq_mpi   (nexchg_   )          &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_nordest)
      call mpi_irecv(ar_recv_nordest(1)     & ! tableau recu
               ,size(ar_recv_nordest(1:k))  &
               ,mpi_double_precision                         &
               ,par%tvoisin(nordest)             &
               ,tagsudouest_mpi                     &
               ,par%comm2d                       &
               ,tabreq_mpi   (nexchg_   )           &
               ,ierr)
      endif

! Coin Nord-Ouest:
      if (par%tvoisin(nordouest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_nordouest)
      call mpi_issend(ar_send_nordouest(1)    & ! envoyé au proc Est
                ,size(ar_send_nordouest(1:k)) &
                ,mpi_double_precision                          &
                ,par%tvoisin(nordouest)            &
                ,tagnordouest_mpi                     &
                ,par%comm2d                        &
                ,tabreq_mpi   (nexchg_   )            &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_nordouest)
      call mpi_irecv(ar_recv_nordouest(1)     & ! tableau recu
               ,size(ar_recv_nordouest(1:k))  &
               ,mpi_double_precision                           &
               ,par%tvoisin(nordouest)             &
               ,tagsudest_mpi                         &
               ,par%comm2d                         &
               ,tabreq_mpi   (nexchg_   )             &
               ,ierr)
      endif

! Coin Sud-Ouest:
      if (par%tvoisin(sudouest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_sudouest)
      call mpi_issend(ar_send_sudouest(1)    & ! envoyé au proc Est
                ,size(ar_send_sudouest(1:k)) &
                ,mpi_double_precision                         &
                ,par%tvoisin(sudouest)            &
                ,tagsudouest_mpi                     &
                ,par%comm2d                       &
                ,tabreq_mpi   (nexchg_   )           &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_sudouest)
      call mpi_irecv(ar_recv_sudouest(1)     & ! tableau recu
               ,size(ar_recv_sudouest(1:k))  &
               ,mpi_double_precision                          &
               ,par%tvoisin(sudouest)             &
               ,tagnordest_mpi                       &
               ,par%comm2d                        &
               ,tabreq_mpi   (nexchg_   )            &
               ,ierr)
      endif

! Coin Sud-Est:
      if (par%tvoisin(sudest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_sudest)
      call mpi_issend(ar_send_sudest(1)    & ! envoyé au proc Est
                ,size(ar_send_sudest(1:k)) &
                ,mpi_double_precision                       &
                ,par%tvoisin(sudest)            &
                ,tagsudest_mpi                     &
                ,par%comm2d                     &
                ,tabreq_mpi   (nexchg_   )         &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_sudest)
      call mpi_irecv(ar_recv_sudest(1)     & ! tableau recu
               ,size(ar_recv_sudest(1:k))  &
               ,mpi_double_precision                        &
               ,par%tvoisin(sudest)             &
               ,tagnordouest_mpi                   &
               ,par%comm2d                      &
               ,tabreq_mpi   (nexchg_   )          &
               ,ierr)
      endif


! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_est)
!     call mpi_issend(ar_send_est(1:k)  & ! envoyé au proc Est
      call mpi_issend(ar_send_est(1)    & ! envoyé au proc Est
                ,size(ar_send_est(1:k)) &
                ,mpi_double_precision                       &
                ,par%tvoisin(est)                           &
                ,tagest_mpi                                    &
                ,par%comm2d                                 &
                ,tabreq_mpi   (nexchg_   )                     &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_est)
!     call mpi_irecv(ar_recv_est(1:k)   & ! tableau recu du proc Ouest
      call mpi_irecv(ar_recv_est(1)     & ! tableau recu du proc Ouest
               ,size(ar_recv_est(1:k))  &
               ,mpi_double_precision                        &
               ,par%tvoisin(est)                            &
               ,tagouest_mpi                                   &
               ,par%comm2d                                  &
               ,tabreq_mpi   (nexchg_   )                      &
               ,ierr)
      endif

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_ouest)
!     call mpi_issend(ar_send_ouest(1:k)  &
      call mpi_issend(ar_send_ouest(1)    &
                ,size(ar_send_ouest(1:k)) &
                ,mpi_double_precision                         &
                ,par%tvoisin(ouest)                           &
                ,tagouest_mpi                                    &
                ,par%comm2d                                   &
                ,tabreq_mpi   (nexchg_   )                       &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_ouest)
!     call mpi_irecv(ar_recv_ouest(1:k)  &
      call mpi_irecv(ar_recv_ouest(1)    &
               ,size(ar_recv_ouest(1:k)) &
               ,mpi_double_precision                         &
               ,par%tvoisin(ouest)                           &
               ,tagest_mpi                                      &
               ,par%comm2d                                   &
               ,tabreq_mpi   (nexchg_   )                       &
               ,ierr)
      endif

! Frontiere nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_nord)
!     call mpi_issend(ar_send_nord(1:k)  &
      call mpi_issend(ar_send_nord(1)    &
                ,size(ar_send_nord(1:k)) &
                ,mpi_double_precision                        &
                ,par%tvoisin(nord)                           &
                ,tagnord_mpi                                    &
                ,par%comm2d                                  &
                ,tabreq_mpi   (nexchg_   )                      &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_nord)
!     call mpi_irecv(ar_recv_nord(1:k)  &
      call mpi_irecv(ar_recv_nord(1)    &
               ,size(ar_recv_nord(1:k)) &
               ,mpi_double_precision                        &
               ,par%tvoisin(nord)                           &
               ,tagsud_mpi                                     &
               ,par%comm2d                                  &
               ,tabreq_mpi   (nexchg_   )                      &
               ,ierr)
      endif

! Frontiere sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nb_send_sud)
      call mpi_issend(ar_send_sud(1)     &
                ,size(ar_send_sud(1:k))  &
                ,mpi_double_precision                        &
                ,par%tvoisin(sud)                            &
                ,tagsud_mpi                                     &
                ,par%comm2d                                  &
                ,tabreq_mpi   (nexchg_   )                      &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nb_recv_sud)
      call mpi_irecv(ar_recv_sud(1)     &
               ,size(ar_recv_sud(1:k))  &
               ,mpi_double_precision                        &
               ,par%tvoisin(sud)                            &
               ,tagnord_mpi                                    &
               ,par%comm2d                                  &
               ,tabreq_mpi   (nexchg_   )                      &
               ,ierr)
      endif

      if(nexchg_   >nexchgmax_mpi   ) &
      stop 'ar_obc nexchg_   >nexchgmax_mpi   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                  &
                      ,tabreq_mpi   (1:nexchg_   ) &
                      ,tstatus_mpi  (:,1:nexchg_ ) &
                      ,ierr)


      end subroutine ruban_mpi_sendrecv

!....................................................................

      subroutine ruban_mpi_cornerlist(type_,order_)
      implicit none
      character*1 type_
      integer order_


! LES POINTS TRACEURS:

      if(type_=='z') then !zzzzzzzzz>

       if(order_==0) then !00000>

          nordest_recv_i(1)=imax  
          nordest_recv_j(1)=jmax 
         sudouest_send_i(1)=2
         sudouest_send_j(1)=2

         sudouest_recv_i(1)=1
         sudouest_recv_j(1)=1
          nordest_send_i(1)=imax-1
          nordest_send_j(1)=jmax-1

        nordouest_recv_i(1)=1  
        nordouest_recv_j(1)=jmax 
           sudest_send_i(1)=imax-1
           sudest_send_j(1)=2

           sudest_recv_i(1)=imax  
           sudest_recv_j(1)=1
        nordouest_send_i(1)=2  
        nordouest_send_j(1)=jmax-1


       endif              !00000>

       if(order_==1) then !11111>

           nordest_recv_i(1:3)=(/ imax   , imax+1 , imax+1 /)
           nordest_recv_j(1:3)=(/ jmax+1 , jmax+1 , jmax   /)
          sudouest_send_i(1:3)=(/ 2      , 3      , 3      /)
          sudouest_send_j(1:3)=(/ 3      , 3      , 2      /)

        nordouest_recv_i(1:3)=(/ 1      , 0      , 0       /)
        nordouest_recv_j(1:3)=(/ jmax+1 , jmax+1 , jmax    /)
           sudest_send_i(1:3)=(/ imax-1 , imax-2 , imax-2  /)
           sudest_send_j(1:3)=(/ 3      , 3      , 2       /)

           sudest_recv_i(1:3)=(/ imax   , imax+1 , imax+1  /)
           sudest_recv_j(1:3)=(/ 0      , 0      , 1       /)
        nordouest_send_i(1:3)=(/ 2      , 3      , 3       /)
        nordouest_send_j(1:3)=(/ jmax-2 , jmax-2 , jmax-1  /)

         sudouest_recv_i(1:3)=(/ 1      , 0      , 0       /)
         sudouest_recv_j(1:3)=(/ 0      , 0      , 1       /)
          nordest_send_i(1:3)=(/ imax-1 , imax-2 , imax-2  /)
          nordest_send_j(1:3)=(/ jmax-2 , jmax-2 , jmax-1  /)

       endif              !11111>

       if(order_==2) then !22222>

          nordest_recv_i(1:5)=(/ imax   , imax+1 , imax+2  , imax+2 , imax+2 /)
          nordest_recv_j(1:5)=(/ jmax+2 , jmax+2 , jmax+2  , jmax+1 , jmax   /)
         sudouest_send_i(1:5)=(/ 2      , 3      , 4       , 4      , 4      /)
         sudouest_send_j(1:5)=(/ 4      , 4      , 4       , 3      , 2      /)

        nordouest_recv_i(1:5)=(/ 1      , 0      , -1      , -1     , -1     /)
        nordouest_recv_j(1:5)=(/ jmax+2 , jmax+2 , jmax+2  , jmax+1 , jmax   /)
           sudest_send_i(1:5)=(/ imax-1 , imax-2 , imax-3  , imax-3 , imax-3 /)
           sudest_send_j(1:5)=(/ 4      , 4      , 4       , 3      , 2      /)

           sudest_recv_i(1:5)=(/ imax   , imax+1 , imax+2  , imax+2 , imax+2 /)
           sudest_recv_j(1:5)=(/ -1     , -1     , -1      , 0      , 1      /)
        nordouest_send_i(1:5)=(/ 2      , 3      , 4       , 4      , 4      /)
        nordouest_send_j(1:5)=(/ jmax-3 , jmax-3 , jmax-3  , jmax-2 , jmax-1 /)

         sudouest_recv_i(1:5)=(/ 1      , 0      , -1      , -1     , -1     /)
         sudouest_recv_j(1:5)=(/ -1     , -1     , -1      , 0      , 1      /)
          nordest_send_i(1:5)=(/ imax-1 , imax-2 , imax-3  , imax-3 , imax-3 /)
          nordest_send_j(1:5)=(/ jmax-3 , jmax-3 , jmax-3  , jmax-2 , jmax-1 /)

       endif              !22222

      endif               !zzzzzzzzz>

! LES POINTS VITESSE U
! A poursuivre...
! LES POINTS VITESSE V
! A poursuivre...

      end subroutine ruban_mpi_cornerlist


!....................................................................

      subroutine ruban_mpi_indexes(type_,order_)
      implicit none
      character*1 type_
      integer        &
              order_   ! (z0=order_=0, z1=order_=1, etc...)

!..........
! les indices min max des bords:
      if(type_=='z') then !zzzzz>

       istr_mpi=2 ; iend_mpi=imax-1 ; jstr_mpi=2 ; jend_mpi=jmax-1

        if(par%tvoisin(nord) ==mpi_proc_null)jend_mpi=jmax
        if(par%tvoisin(sud)  ==mpi_proc_null)jstr_mpi=1
        if(par%tvoisin(est)  ==mpi_proc_null)iend_mpi=imax
        if(par%tvoisin(ouest)==mpi_proc_null)istr_mpi=1

      endif               !zzzzz>
!..........

!..........
! les indices (i,j) des coins:
      call ruban_mpi_cornerlist(type_,order_)
!..........

!..........
! Nombres d'elements echangEs
! Par les coins:
      nb_send_nordest=  kmax*(1+2*order_)
! Par les bords:
      nb_send_est= kmax*(jend_mpi-jstr_mpi+1)
      nb_send_nord=kmax*(iend_mpi-istr_mpi+1)
!..........
! Ici on a generalisE en envisageant la possibilitE que tous les bords et coins n'echangent pas
! pas le meme nombre d'elements en prevision d'une distribution mpi un jour plus complexe. Dans le
! cas present les echanges sont homogenes:
      nb_send_nordouest=nb_send_nordest
      nb_send_sudest   =nb_send_nordest
      nb_send_sudouest =nb_send_nordest
      nb_recv_nordest  =nb_send_nordest
      nb_recv_nordouest=nb_send_nordest
      nb_recv_sudest   =nb_send_nordest
      nb_recv_sudouest =nb_send_nordest

      nb_send_ouest=nb_send_est
      nb_recv_est  =nb_send_est
      nb_recv_ouest=nb_send_est
      nb_send_sud  =nb_send_nord
      nb_recv_nord =nb_send_nord
      nb_recv_sud  =nb_send_nord


      end subroutine ruban_mpi_indexes

!....................................................................

      subroutine ruban_mpi_tem(type_,order_,time_)
      implicit none
      character*1 type_
      integer        &
              order_ & ! (z0=order_=0, z1=order_=1, etc...)
             ,time_    ! Temps (now, before,...) de la variable


      call ruban_mpi_indexes(type_,order_) ! donne les indices des zones A Echanger (jstr_mpi,jend_mpi etc..)
                                           ! et compte les elements A Echanger

      if(.not.allocated(ar_send_est))call ruban_mpi_allocate
 
!..........
! Charger le ruban envoyeur avec les bords et coins de la variable A envoyer aux voisins:
      
      if (par%tvoisin(est) /= mpi_proc_null) then
! Etre sur le bord est et etre envoyeur 
       k1=0
       do k=1,kmax
       do j=jstr_mpi,jend_mpi
        k1=k1+1
        ar_send_est(k1)=tem_t(imax-1-order_,j,k,time_)  ! zO envoie imax-1, z1 imax-2, z2 imax-3
       enddo
       enddo
      endif

      if (par%tvoisin(ouest) /= mpi_proc_null) then
! Etre sur le bord ouest et etre envoyeur 
       k1=0
       do k=1,kmax
       do j=jstr_mpi,jend_mpi
        k1=k1+1
        ar_send_ouest(k1)=tem_t(2+order_,j,k,time_) ! z0 envoie 2, z1 3, z2 4
       enddo
       enddo
      endif

      if (par%tvoisin(nord) /= mpi_proc_null) then
! Etre sur le bord nord et etre envoyeur 
       k1=0
       do k=1,kmax
       do i=istr_mpi,iend_mpi
        k1=k1+1
        ar_send_nord(k1)=tem_t(i,jmax-1-order_,k,time_) 
       enddo
       enddo
      endif

      if (par%tvoisin(sud) /= mpi_proc_null) then
! Etre sur le bord sud et etre envoyeur 
       k1=0
       do k=1,kmax
       do i=istr_mpi,iend_mpi
        k1=k1+1
        ar_send_sud(k1)=tem_t(i,2+order_,k,time_) 
       enddo
       enddo
      endif

      if (par%tvoisin(nordest) /= mpi_proc_null) then
! Etre sur le coin nord-est et etre envoyeur
       k1=0
       do k2=1,(1+2*order_)
          i=nordest_send_i(k2)
          j=nordest_send_j(k2)
          do k=1,kmax
            k1=k1+1
            ar_send_nordest(k1)=tem_t(i,j,k,time_) 
          enddo
       enddo
      endif

      if (par%tvoisin(nordouest) /= mpi_proc_null) then
! Etre sur le coin nord-ouest et etre envoyeur
       k1=0
       do k2=1,(1+2*order_)
          i=nordouest_send_i(k2)
          j=nordouest_send_j(k2)
          do k=1,kmax
            k1=k1+1
            ar_send_nordouest(k1)=tem_t(i,j,k,time_) 
          enddo
       enddo
      endif

      if (par%tvoisin(sudest) /= mpi_proc_null) then
! Etre sur le coin sud-est et etre envoyeur
       k1=0
       do k2=1,(1+2*order_)
          i=sudest_send_i(k2)
          j=sudest_send_j(k2)
          do k=1,kmax
            k1=k1+1
            ar_send_sudest(k1)=tem_t(i,j,k,time_) 
          enddo
       enddo
      endif

      if (par%tvoisin(sudouest) /= mpi_proc_null) then
! Etre sur le coin sud-ouest et etre envoyeur
       k1=0
       do k2=1,(1+2*order_)
          i=sudouest_send_i(k2)
          j=sudouest_send_j(k2)
          do k=1,kmax
            k1=k1+1
            ar_send_sudouest(k1)=tem_t(i,j,k,time_) 
          enddo
       enddo
      endif
!..........



!..........
! FAIRE L'ECHANGE MPI:
      call ruban_mpi_sendrecv
!..........



!..........
! Charger les bords et coins de la variable tem_t avec le ruban recu

      if (par%tvoisin(est) /= mpi_proc_null) then
! receveur, bord est
       k1=0
       do k=1,kmax
       do j=jstr_mpi,jend_mpi
        k1=k1+1
        tem_t(imax+order_,j,k,time_)=ar_recv_est(k1)
       enddo
       enddo
      endif

      if (par%tvoisin(ouest) /= mpi_proc_null) then
! receveur, bord ouest
       k1=0
       do k=1,kmax
       do j=jstr_mpi,jend_mpi
        k1=k1+1
        tem_t(1-order_,j,k,time_)=ar_recv_ouest(k1)
       enddo
       enddo
      endif

      if (par%tvoisin(nord) /= mpi_proc_null) then
! receveur, bord nord
       k1=0
       do k=1,kmax
       do i=istr_mpi,iend_mpi
        k1=k1+1
        tem_t(i,jmax+order_,k,time_)=ar_recv_nord(k1)
       enddo
       enddo
      endif

      if (par%tvoisin(sud) /= mpi_proc_null) then
! receveur, bord sud
       k1=0
       do k=1,kmax
       do i=istr_mpi,iend_mpi
        k1=k1+1
        tem_t(i,1-order_,k,time_)=ar_recv_sud(k1)
       enddo
       enddo
      endif

      if (par%tvoisin(nordest) /= mpi_proc_null) then
! receveur, coin nord-est
       k1=0
       do k2=1,(1+2*order_)
          i=nordest_recv_i(k2)
          j=nordest_recv_j(k2)
          do k=1,kmax
            k1=k1+1
            tem_t(i,j,k,time_)=ar_recv_nordest(k1)
          enddo
       enddo
      endif

      if (par%tvoisin(nordouest) /= mpi_proc_null) then
! receveur, coin nord-ouest
       k1=0
       do k2=1,(1+2*order_)
          i=nordouest_recv_i(k2)
          j=nordouest_recv_j(k2)
          do k=1,kmax
            k1=k1+1
            tem_t(i,j,k,time_)=ar_recv_nordouest(k1)
          enddo
       enddo
      endif

      if (par%tvoisin(sudest) /= mpi_proc_null) then
! receveur, coin sud-est
       k1=0
       do k2=1,(1+2*order_)
          i=sudest_recv_i(k2)
          j=sudest_recv_j(k2)
          do k=1,kmax
            k1=k1+1
            tem_t(i,j,k,time_)=ar_recv_sudest(k1)
          enddo
       enddo
      endif

! receveur, coin sud-ouest
      if (par%tvoisin(sudouest) /= mpi_proc_null) then
       k1=0
       do k2=1,(1+2*order_)
          i=sudouest_recv_i(k2)
          j=sudouest_recv_j(k2)
          do k=1,kmax
            k1=k1+1
            tem_t(i,j,k,time_)=ar_recv_sudouest(k1)
          enddo
       enddo
      endif
!..........

      end subroutine ruban_mpi_tem


!....................................................................

      end module module_ruban_mpi
