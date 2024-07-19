      subroutine the_end(ichoix)
!______________________________________________________________________
! SYMPHONIE ocean model
! release S26 - last update: 09-06-18
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      use module_modeanalysis
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='the_end'
       subroutinedescription= &
       'Concludes the simulation. Writes the last outputs.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         14/01/02: bienvenue à l'analyse en harmoniques de zeta archivee dans
!                   film.data (on suppose que zeta est la premiere variable
!                   visualisee dans animation.F)
!         30/01/02: amenagement analyse harmonique
!         27/01/03: appels à CALL GRAPH_OUT
!         05/08/03: création d'un fichier de vérification de passage dans la
!                   routine (pratique quand on utilise le cluster)
!         29/09/03: bienvenue à CALL ANALYSEHARMONIQUE(2)
!         01/10/03: appel à graphique à la fin du programme
!         08/07/04: introduction routines ondelettes de Francis
!         21/03/07: Compter le nbre de lignes dans time.dat pour allouer
!                   la memoire.
!         22/01/08: debug allocation taille memoire
!         13/05/08: fermer proprement le fichier interrupteur
! 2009.3  08-11-09: optimobc_tide est renommé tide_analysis
! 2010.3  12-01-10: ecrire dans tmp
! 2010.8  09-05-10  seul le proc 0 ecrit interrupteur
! 2010.10 23-06-10  operations sur interrupteur regroupees dans io_switch
! 2010.11 27-06-10  time.dat devient dom_c//time.dat
! S26     07-02-13  Seul proc 0 ecrit a l'ecran
!         26-02-13  Message de fin de simulation
!         21-10-13  Message de fin de simulation
!         11-11-14  Ne pas appeler l'analyse harmonique si la duree d'analyse est
!                   insuffisante
!         19-04-15  Message ecran si tideana_yesno==1
!         02-07-15  fichier kount dans repertoire tmp
!         29-05-17  modif de message de fin de simulation
!         09-06-18 ajout return
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m°v°m
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


      if(ichoix.eq.0) then

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'the_end.dat')                              !12-01-10
      write(3,*)'entrée dans the_end.f ichoix.eq.0'                    !05/08/03
      close(3)
      endif                !#mpi-->>-->                       !09-05-10

!        kstop=0
!        kpvwave=0
!        give_chanel9=0                                                !13/05/08
!        call io_switch('w')                                           !23-06-10

! Analyse Harmonique:
      if(kmaxtide>0) then !------>
        if(flag_3dwaves_harmonics==1)call modeanalysis_harmonics(2)    !26-02-13
        if(tideana_yesno==1)then !AAAAA>               !19-04-15
          if(tideanalysis_count>2*kmaxtide) then !>>>> !11-11-14
           call tide_analysis(6)
          else                                   !>>>>
           if(par%rank==0) then                               !11-11-14
            write(6,'(a,i0,a)')'The number of analysed fields, '      &
                     ,tideanalysis_count                     &
                     ,', is too small to compute the tidal harmonics'
           endif
          endif                                  !>>>>
        endif                    !AAAAA>
      endif               !------>

      if(par%rank==0) then !#mpi-->>-->                                !09-05-10
      open(unit=3,file=trim(tmpdirname)//'the_end.dat',position='append')            !12-01-10
      write(3,*)'the_end.f ichoix.eq.0 avant graphique'                !05/08/03
      close(3)
      endif                !#mpi-->>-->                       !09-05-10

                call graph_out                          !              !01/10/03
                if(par%rank==0)write(6,*)'dans sbr the_end'            !

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'the_end.dat',position='append')            !12-01-10
      write(3,*)'the_end.f ichoix.eq.0 apres graphique'                !05/08/03
      close(3)

      open(unit=3,file=trim(tmpdirname)//'the_end.dat',position='append')            !12-01-10
      write(3,*)'the_end.f ichoix.eq.0 fin'                            !05/08/03
      close(3)
      endif                !#mpi-->>-->                       !09-05-10

      return
      endif

      if(ichoix==1) then
        if(par%rank==0) then !---------->

        if(kstop==0) then !ooo>
         open(unit=3,file=trim(tmpdirname)//'kount',position='append')  !O2-01-15
          write(3,*)
          write(3,*)'Congratulations for this successful run !'         !26-02-13
          write(3,*)
         close(3)
         write(6,*)'Congratulations for this successful run !'          !21-10-13
         write(6,*)
        endif             !ooo>

        if(kstop==1) then !pmx> !29-05-17
         open(unit=3,file=trim(tmpdirname)//'kount',position='append')  !O2-01-15
          write(3,*)
          write(3,*)'Run stopped by interrupteur file'         !26-02-13
          write(3,*)
         close(3)
         write(6,*)'Run stopped by interrupteur file'         !26-02-13
         write(6,*)
        endif             !pmx>

        return !09-06-18
        endif               !---------->

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
       stop
      endif

      end subroutine the_end
