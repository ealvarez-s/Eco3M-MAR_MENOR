










      subroutine chronos
!______________________________________________________________________
! S model
! release S26 - last update: 27-02-13
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      implicit none

!______________________________________________________________________________
! Version date      Description des modifications
!         26/04/05: interrupteur permet maintenant de creer un fichier restart
!                   sans arreter le model_
!         02/08/05: modif message à l'ecran
!         26/03/07: Calcul d'une rampe commune
!         04/03/08: Modif affichage (format) dans fichier kount
!         01/06/08: dans fichier "kount" KOUNT0_RST remplace KOUNT0
! 2009.3  02-10-09: dte_lp remplace dte
! 2010.8  09-05-10  seul le proc 0 peut ecrire fichier kount
! 2010.10 13-06-10  suite point precedent: tous les proc peuvent lire interrupteur
!         21-06-10  operations sur interrupteur regroupees dans io_switch
! 2010.14 20-11-10  ecriture de la date dans fichier kount
!                   kount=kount+1 apres ecritures fichier & ecran
! 2010.19 03-04-11  ajout variable elapsedtime
! 2010.20 16-04-11  elapsedtime_aft
!         19-04-11  Cas 2d obsolete
! 2010.22 27-04-11  Correction ecriture date de fin dans fichier "kount"
! S25.4   30-06-12  Seul le proc 0 ecrit à l'ecran
! S26     27-02-13  Messages dans kount
!______________________________________________________________________________

! OUTPUT: KOUNT

! modifs:
! 25/07/01: chronos.F ne marchait que dans le cas où l'heure, la minute
! la seconde du jour de reference etaient nuls (cas normal). Il a été
! généralisé au cas où ces 3 paramètres seraient non nuls.
! 29/08/01: le calcul de la date est fait par le sous programme kount_to_date.F

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE3D:
! DEBUT:
      if(i2dh.eq.1) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(mod(iteration3d,5).eq.0) then !*******************************>

! afficher la date à l'écran:
      if(par%rank==0)                                     &     !30-06-12
      write(6,'(1x,i2,1x,a,1x,i4,1x,i2,a1,i2,a1,i2,a1)')  &
         day_now,trim(month(month_now)),year_now          &     !20-11-10
        ,hour_now,'h'                                     &
        ,minute_now,'m',second_now,'s'


      endif                       !*******************************>

!____________________________________________________________
!////////// FICHIER KOUNT ///////////////////////////////////

      if(par%rank==0) then !#mpi>>>>>>>>>>>>>

      open(unit=3,file=trim(tmpdirname)//'kount')

      write(3,'(1x,i2,1x,a,1x,i4,1x,i2,a1,i2,a1,i2,a1)')  &
         day_now,trim(month(month_now)),year_now                 &     !20-11-10
        ,hour_now,'h'                                      &
        ,minute_now,'m',second_now,'s'
      write(3,*)'Iterative loop:                     ',iteration3d
      write(3,*)'Elapsedtime(days):                  '   &
                ,elapsedtime_now/86400.
      write(3,*)'Elapsedtime(days) from last restart:'   &
                ,(elapsedtime_now-elapsedtime_rst)/86400.
      write(3,*)'Elapsedtime(%)    from last restart:'   &
           ,100.*(elapsedtime_now-elapsedtime_rst)       &
                /(elapsedtime_end-elapsedtime_rst)

      write(3,'(a6,i2,1x,a,1x,i4,1x,i2,a1,i2,a1,i2,a1)')          & !27-02-13
         ' End: '                                                 &
        ,datesim(3,2),trim(month(datesim(2,2))),datesim(1,2)      &     !27-04-11
        ,datesim(4,2),'h'                                         &
        ,datesim(5,2),'m',datesim(6,2),'s'

      close(3)

      endif                !#mpi>>>>>>>>>>>>>

      call io_switch('r')  !(interrupteur)              !21-06-10

!////////// FICHIER KOUNT ///////////////////////////////////
!------------------------------------------------------------

!     kount=kount+1                                                     !20-11-10
      iteration3d=iteration3d+1

! Desormais calculés dans time_step_next:
!     elapsedtime_bef=elapsedtime_now                                   !03-04-11
!     elapsedtime_now=elapsedtime_now+dti_now                           !03-04-11
!     elapsedtime_aft=elapsedtime_now+dti_now                           !16-04-11
!     rampe=min(max(elapsedtime_now/86400.,zero),un)                    !03-04-11

!     rampe=min(max((kount -kount0)*dti_fw/86400.,zero),un)             !26/03/07

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE3D:
! FIN.
      return
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!                         /        /         /

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE2D:
! DEBUT:
      if(i2dh.eq.0) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      stop 'cas obsolete'     !19-04-11

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE2D:
! FIN.
      return
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine chronos
