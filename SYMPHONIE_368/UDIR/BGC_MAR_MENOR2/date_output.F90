      subroutine date_output(ichoix)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 257 - last update: 03-07-19
!______________________________________________________________________


      use module_principal
      use module_parallele !#mpi
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='date_output'
       subroutinedescription= &
       'Decides to call graph_out subroutine if the time conditions' &
       //' specified in notebook_grap are satisfied'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!______________________________________________________________________
! Version Date      Description des modifications
!         14/01/03: mise en service
!         16/01/03: afficher la date
!                   le calcul du nom du fichier est retabli dans graph_out.F
!         27/01/03: cas du model_1D pris en compte
!         10/07/03: ecriture d'un fichier chanel9 en même temps qu'un fichier
!                   graphique à échéance régulière.
!         25/04/05: appel à hot_restart commenté
!         19/01/07: compatible avec model_ 1D
!         30/01/08: appel à hot_restart reactivé
!         23/01/09: appel à hot_restart reactivé pour le cas "sorties aux 8 dates"
!         14-06-09: Pour etre sur ne de pas appeler hot_restart depuis un mauvais
!                   endroit, remplacer par give_chanel9=1
!         18-06-09: Cas 1D: appel à graph1D
! 2010.7  13-02-10  affichage ecran
! 2010.8  09-05-10  seul le proc 0 ecrit interrupteur
! 2010.10 21-06-10  operations sur interrupteur regroupees dans io_switch
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.20 15-04-11  Calculs sur la base d'un temps en secondes
! 2010.23 20-05-11  On ne demande plus l'ecriture d'un fichier restart quand
!                   on ecrit un fichier graphique
! S.26    05-11-13  seul par%rank=0 ecrit a l'ecran
!         02-11-14  appel graph_out_bio
!         08-11-14  graph_out_bio appele depuis graph_out
! v257    03-07-19  texte90=trim(directory)//txtslash//'notebook_dateoutput' 

!...............................................................................
!    _________                    .__                  .__             !m°v°w
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


      if(ichoix.eq.0)stop 'date_output: filière 8 Kounts obsolete'
!*******************************************************************************
! FAIRE DES SORTIES DU MODELE REGULIERES
! DEBUT:
      if(ichoix.eq.1) then
!*******************************************************************************

!     if(int(real(kount  )/kountgraph)                                  &
!       -int(real(kount-1)/kountgraph).eq.1) then !>>>>>>>>>>>>>>>>>>>

!      if( int(elapsedtime_now/graphperiod)       &
!         -int(elapsedtime_bef/graphperiod)==1) then !>>>>>>>>>>>>>>>>>>> !15-04-11

      if( int(elapsedtime_now/graphperiod)        &
         -int(elapsedtime_bef/graphperiod)==1.or. &
          int(elapsedtime_now/graphperiod_supp)        &
         -int(elapsedtime_bef/graphperiod_supp)==1) then !>>>>>>>>>>>>>>>>>>> !15-04-11



      if(flag3d==1) then !>>>>>
        ! if(imodelbio==1)call graph_out_bio                     !02-11-14
        ! if(imodeltrc==1)call graph_out_bio                     !08-11-14
        call graph_out                                         !18-06-09
      endif           !>>>>>
      if(flag3d==0)call graph1d_out                                       !18-06_09
      if(par%rank==0)write(6,*)'dans sbr date_output'
      if(par%rank==0)write(6,*)'-----------------------'

      endif                                        !>>>>>>>>>>>>>>>>>>>

!*******************************************************************************
! FAIRE DES SORTIES DU MODELE REGULIERES
! FIN.
      return
      endif
!*******************************************************************************


!                            /   /   /


!*******************************************************************************
! FAIRE DES SORTIES DU MODELE A DES DATES SPECIFIEES DANS notebook_output
! DEBUT:
      if(ichoix.eq.2) then
!*******************************************************************************

      ksecu=0
      do 114 k=1,kmax_dof
!      if(kount.eq.kdate_output(k)) then !§§§§§§§§§>
       if(     elapsedtime_now> tdate_output(k)     &
          .and.elapsedtime_bef<=tdate_output(k)) then !§§§§§§§§§>      !15-04-11
       ksecu=1
       goto 117
       endif                                          !§§§§§§§§§>
  114 continue
      return  ! si KSECU=0 pas de graphiques
  117 continue

      if(flag3d==1) then !>>>>>
        ! if(imodelbio==1)call graph_out_bio                     !02-11-14
        ! if(imodeltrc==1)call graph_out_bio                     !08-11-14
        call graph_out                                         !18-06-09
      endif           !>>>>>

      if(flag3d==0)call graph1d_out                                       !18-06_09
      if(par%rank==0)write(6,*)'dans sbr date_output'
      if(par%rank==0)write(6,*)'-----------------------'
!               stop 'date_output'

!*******************************************************************************
! FAIRE DES SORTIES DU MODELE A DES DATES SPECIFIEES DANS notebook_output
! FIN.
      return
      endif
!*******************************************************************************


!                        /   /   /


!*******************************************************************************
! FAIRE DES SORTIES DU MODELE A LA FIN DU RUN OU A LA DEMANDE INTERACTIVE
! DEBUT:
      if(ichoix.eq.3) then
!*******************************************************************************

      if(flag3d==1) then !>>>>>
         if(imodelbio==1)call graph_out_bio                     !02-11-14
        ! if(imodeltrc==1)call graph_out_bio                     !08-11-14
        call graph_out                                         !18-06-09
      endif           !>>>>>

      if(flag3d==0)call graph1d_out                                       !18-06_09
      if(par%rank==0)write(6,*)'dans sbr date_output'
      if(par%rank==0)write(6,*)'-----------------------'

!*******************************************************************************
! FAIRE DES SORTIES DU MODELE A LA FIN DU RUN OU A LA DEMANDE INTERACTIVE
! FIN.
      return
      endif
!*******************************************************************************


!                        /   /   /


!*******************************************************************************
! PHASE D'INITIALISATION: LECTURE DE notebook_dateoutput
! DEBUT:
      if(ichoix.eq.4) then
!*******************************************************************************

!     texte90=directory(1:lname1)//txtslash//'notebook_dateoutput'
      texte90=trim(directory)//txtslash//'notebook_dateoutput' !03-07-19
      open(unit=3,file=texte90)
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
       do 184 loop1=1,99999
        read(3,*,end=185)i1,i2,i3,i4,i5,i6
         if(loop1.gt.dim_dof) then !$$$$$$$$$$$$$$$$$$$$$$$$$>
          if(par%rank==0) then
           write(6,*)'nombre de dates dans notebook_dateoutput plus'
           write(6,*)'plus grand que dim_dof dans parameter.'
           write(6,*)'corriger puis relancer.'
          endif
          stop ' dans date_output.f'
         endif                     !$$$$$$$$$$$$$$$$$$$$$$$$$>
        if(par%rank==0)write(6,*) &
         loop1,'eme appel a date_to_kount dans date_output.f' !13-02-10
!       call date_to_kount
        call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!       kdate_output(loop1)=kdtk_out
        tdate_output(loop1)=elapsedtime_out    !15-04-11
  184   continue
  185  continue
      kmax_dof=loop1-1
      close(3)

!*******************************************************************************
! PHASE D'INITIALISATION: LECTURE DE notebook_dateoutput
! FIN.
      return
      endif
!*******************************************************************************

      end subroutine date_output
