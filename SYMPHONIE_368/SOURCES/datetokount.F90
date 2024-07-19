      subroutine datetokount(year_,month_,day_,hour_,minute_,second_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 326 - last update: 15-12-21
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
      double precision x0_,x1_,x2_,deltaday_
      integer kdtk0_,kdtk2_          &
             ,choix_,j0_,j1_,j2_,k1_  &
             ,k2_,ksecu_,loop_              &
             ,year_,month_,day_             &
             ,hour_,minute_,second_
#ifdef synopsis
       subroutinetitle='datetokount'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!......................................................................
! Version date     Description des modifications
!         25/09/01 creation de variable propre au sous programme
!         05/01/03 ajout d'un debuggeur si date non referencee dans dayindex
! 2010.2  19-12-09 trois methodes au choix pour trouver le numero du jour
!                  correspondant à la date entrante. La methode choix_=2
!                  est retenue car elle converge plus vite que les 2 autres.
! 2010.16 12-01-11 Contraindre j0_ à rester dans l'intervalle de definition
!                  de dayindex
! 2010.20 14-04-14 convertir la date en un temps ecoulé en secondes depuis
!                  le temps initial de la simulation
! S26     25-01-13 Trouver la date en dehors de la gamme de definition de dayindex
!         03-04-13 Debug du point precedent via introduction de routine datetokount_basic
!         12-10-17 Debug grands entiers
! v249    05-03-19 Ne pas calculer de conversion de date si cette derniere n'est pas
!                  definie
! v269    06-12-19 dti_now remplacE par dti_fw
! v326    15-12-21 write(10+par%rank,*)'Err 162 xdtk_out=',xdtk_out 
!.......................................................................
!    _________                    .__                  .__             !(°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!.......................................................................

! Connaissant I1,I2,I3,I4,I5,I6 correspondant respectivement à
! l'année, le mois, le jour, l'heure, la minute, la seconde
! date_to_kount.F retourne:

!   un numero d'iteration entier:  KDTK_OUT
!   et décimal: XDTK_OUT
!   elapsedtime_out: temps en secondes depuis le debut de la simulation

! attention, ce programme ne peut fonctionner que si le
! tableau DAYINDEX est initialisé, ce qui est fait dans time_step.F
! attention, ce programme ne peut fonctionner que si
! DAYSIM(0) et NUMBYEARS sont initialisés, ce qui est fait
! dans time_step.F
!.......................................................................

      if(year_==0.or.month_==0.or.day_==0)return !05-03-19

      choix_=2     ! selection de la methode
!     choix_=-999  ! selection de la methode
      ksecu_=0                                                          !05/01/03

      kdtk0_=year_*10000+month_*100+day_

      if(choix_==0) then !0000000000000000000000000000000>
! Par selection successive du demi-intervalle encadrant la date entrante

      j1_=1
      j2_=dayindex_size

      do loop_=1,dayindex_size

      j0_=(j1_+j2_)/2
      if ( kdtk0_ > dayindex(j0_) ) j1_=j0_
      if ( kdtk0_ < dayindex(j0_) ) j2_=j0_

      if ( kdtk0_ == dayindex(j0_)   ) then   !----->
           deltaday_=j0_
           ksecu_=1
           goto 97
      endif                                 !-----<

      if ( kdtk0_ == dayindex(j0_+1)   ) then !----->
           deltaday_=j0_+1
           ksecu_=1
           goto 97
      endif                                 !-----<

      enddo ! fin de boucle sur loop_
   97 continue

      endif              !0000000000000000000000000000000>


      if(choix_==1) then !1111111111111111111111111111111>
! méthode la plus simple. On compare un à un les valeurs de dayindex à la date entrante

      do 180 kdtk2_=1,366*numbyears
      if(kdtk0_.eq.dayindex(kdtk2_))then
      deltaday_=kdtk2_
      ksecu_=1                                                          !05/01/03
      goto 181
      endif
  180 continue
  181 continue

      endif              !1111111111111111111111111111111>

      if(choix_==2) then !2222222222222222222222222222222>
! méthode par regle de 3:

      x0_=real(year_)+(month_-1.+(day_-1.)/31.)/12.

      j0_=1

      do loop_=1,dayindex_size

      j1_=min(j0_,dayindex_size-1)                          !12-01-11
      k1_=dayindex(j1_)-int(dayindex(j1_)/10000)*10000
      k2_=k1_-int(k1_/100)*100
      x1_=real(dayindex(j1_)/10000)+(real(k1_/100)-1.+(k2_-1.)/31.)/12.

      j1_=j1_+1                                             !12-01-11
      k1_=dayindex(j1_)-int(dayindex(j1_)/10000)*10000
      k2_=k1_-int(k1_/100)*100
      x2_=real(dayindex(j1_)/10000)+(real(k1_/100)-1.+(k2_-1.)/31.)/12.

      j0_=max(min( j0_+nint((x0_-x1_)/(x2_-x1_)) &
                 ,dayindex_size),1)

      if(dayindex(j0_)==kdtk0_) then
      ksecu_=1
      deltaday_=j0_
      goto 156
      endif

      enddo ! fin de boucle sur loop_

  156 continue

      endif              !2222222222222222222222222222222>

      if(ksecu_/=0) then !--- debuggeur ----->                        !05/01/03
!---> DAYSIM est le numéro d'index (décimal) du jour référencé dans DAYINDEX
      elapsedtime_out=( deltaday_+(real(hour_)*3600.            & !14-04-11
                       +real(minute_)*60.+real(second_))/86400. &
                       -daysim(0)  )*3600.*24.
!     xdtk_out=elapsedtime_out/dti_now                                  !14-04-11
      xdtk_out=elapsedtime_out/dti_fw !06-12-19 (car dti_now peut etre nul)
      if(abs(xdtk_out)<huge(kdtk_out)) then !m°v°m> !15-12-21
       kdtk_out=nint(xdtk_out)
      else                                  !m°v°m>
! Si le nombre depasse la valeur max d'un entier, prendre 0 et faire une alerte dans fort.xxx
       kdtk_out=0
!      write(10+par%rank,*)'Err 162 xdtk_out=',xdtk_out & !15-12-21
!       ,' trop grand pour calculer kdtk_out'
      endif                                 !m°v°m>
   

      else               !--- debuggeur ----->                        !05/01/03

! Procedure de rattrapage si hors intervale de definition de dayindex   !25-01-13
       call datetokount_basic(year_,month_,day_,hour_,minute_,second_)

      endif              !--- debuggeur ----->                        !05/01/03

      end subroutine datetokount

!.............................................................

      subroutine datetokount_basic(year_,month_,day_,hour_,minute_,second_)
      use module_principal
      use module_parallele
      implicit none
      double precision day0deci_,day1deci_
      integer year_,month_,day_,hour_,minute_,second_,avantapres_    &
             ,d0_,mo0_,daycount_,a0_,stop_,a1_,mo1_,d1_                &
             ,h0_,mn0_,s0_     &
             ,h1_,mn1_,s1_     &
             ,a_,mo_,jpa_
#ifdef synopsis
       subroutinetitle='datetokount_basic'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      avantapres_=1

      if(second_>datesim(6,1))avantapres_=1
      if(second_<datesim(6,1))avantapres_=-1
      if(minute_>datesim(5,1))avantapres_=1
      if(minute_<datesim(5,1))avantapres_=-1
      if(hour_  >datesim(4,1))avantapres_=1
      if(hour_  <datesim(4,1))avantapres_=-1
      if(day_   >datesim(3,1))avantapres_=1
      if(day_   <datesim(3,1))avantapres_=-1
      if(month_ >datesim(2,1))avantapres_=1
      if(month_ <datesim(2,1))avantapres_=-1
      if(year_  >datesim(1,1))avantapres_=1
      if(year_  <datesim(1,1))avantapres_=-1

      if(avantapres_==-1) then !>>>>>>>>>>>>>>>>>>

       a0_=year_        ; mo0_=month_       ; d0_=day_
       h0_=hour_        ; mn0_=minute_      ; s0_=second_

       a1_=datesim(1,1) ; mo1_=datesim(2,1) ; d1_=datesim(3,1)
       h1_=datesim(4,1) ; mn1_=datesim(5,1) ; s1_=datesim(6,1)

      else                     !>>>>>>>>>>>>>>>>>

       a1_=year_        ; mo1_=month_       ; d1_=day_
       h1_=hour_        ; mn1_=minute_      ; s1_=second_

       a0_=datesim(1,1) ; mo0_=datesim(2,1) ; d0_=datesim(3,1)
       h0_=datesim(4,1) ; mn0_=datesim(5,1) ; s0_=datesim(6,1)

      endif                    !>>>>>>>>>>>>>>>>>

      daycount_=0


! Trouver le numero decimal du jour anterieur dans l'année a0_
! Annees bisextiles:
      a_=a0_
      if(mod(a_,4)==0) then
              jpm(2)=29
       else
              jpm(2)=28
      endif
      if(mod(a_,100)==0)jpm(2)=28 ! Exception1
      if(mod(a_,400)==0)jpm(2)=29 ! Exception2
      day0deci_=0.
      do mo_=1,mo0_-1
       day0deci_=day0deci_+jpm(mo_)
      enddo
      day0deci_=day0deci_+d0_-1.+real(h0_)/24.+real(mn0_)/1440.+real(s0_)/86400.

! Trouver le numero decimal du jour posterieur dans l'année a0_
      day1deci_=0.
      do a_=a0_,a1_-1
        if(mod(a_,4)==0) then
                jpa_=366
         else
                jpa_=365
        endif
        if(mod(a_,100)==0)jpa_=365  ! jpm(2)=28 ! Exception1
        if(mod(a_,400)==0)jpa_=366  ! jpm(2)=29 ! Exception2
       day1deci_=day1deci_+jpa_
      enddo ! fin de boucle sur a_
      a_=a1_
      if(mod(a_,4)==0) then
              jpm(2)=29
       else
              jpm(2)=28
      endif
      if(mod(a_,100)==0)jpm(2)=28 ! Exception1
      if(mod(a_,400)==0)jpm(2)=29 ! Exception2
      do mo_=1,mo1_-1
       day1deci_=day1deci_+jpm(mo_)
      enddo
      day1deci_=day1deci_+d1_-1.+real(h1_)/24.+real(mn1_)/1440.+real(s1_)/86400.

      elapsedtime_out=(day1deci_-day0deci_)*86400.*avantapres_


      xdtk_out=elapsedtime_out/dti_now                                  !14-04-11

      if(abs(xdtk_out)<214748364) then !12-10-17
         kdtk_out=nint(xdtk_out)
      else
         write(6,*)'WARNING abs(xdtk_out)>214748364 => kdtk_out=0'
         kdtk_out=0
      endif

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif

      end subroutine datetokount_basic
