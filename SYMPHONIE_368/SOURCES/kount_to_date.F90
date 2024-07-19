      subroutine kount_to_date(ki1)
!______________________________________________________________________
! S model
! release S.26  - last update: 24-12-13
!______________________________________________________________________

      use module_principal
      implicit none
      integer ki1
#ifdef synopsis
       subroutinetitle='kount_to_date'
       subroutinedescription= &
       'Converts time iteration into calendar date'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Converts iteration into calendar time

!..............................................................................
! Version date      Description des modifications
!         19/08/03: debug date: calcul sur modulo entier plutot que
!                   prendre partie entiere de modulo reel
!         20/08/03: debug: KI1 ne peut etre inferieur ou egal a zero
!         21/08/03: suite...
!         11/01/05: Calcul plus precis pour eviter bug de "minuit" (reste à
!                   à faire en 2D)
!         28/12/05: calcul forcé en double precision
! 2009.3  02-10-09: dte_lp remplace dte
! S.26    24-12-13  Ne pas passer par cette routine si ki1/=0
!..............................................................................

!.......................................................................
! Connaissant le numéro d'itération passé en argument KI1
! ce sous programme donne la date calendaire.
! I1 secondes
! I2 minutes
! I3 heures
! I7 jours
! I6 mois
! I5 années
!.......................................................................

      if(ki1/=0) then !----------------------------->
       write(6,*)'---------------------------------------------'
       write(6,*)'Argument dans kount_to_date=',ki1
       write(6,*)'Ne pas appeler routine kount_to_date' &
                ,' avec en argument une iteration non nulle.' &
                ,' A la place utiliser elapsedtimetodate'
       stop ' STOP dans routine kount_to_date'
      endif           !----------------------------->

!     IF(KI1.LE.0)THEN                                                 !20/08/03
      if(ki1.lt.0)then                                                 !21/08/03
      write(6,*)
      write(6,*)'anomalie reperée dans kount_to_date.f'
      write(6,*)'recherche de date pour iteration negative.'
      write(6,*)'en effet on passe ',ki1,' en argument.'
      write(6,*)'suggestions:'
      write(6,*)'verifier notebook_time:'
      write(6,*)'  date de debut anterieure à date reference ?'
      write(6,*)'verifier notebook_obcforcing:'
      write(6,*)'  date de premiere sortie (menu imbriqué) prématurée'
      write(6,*)' (notament à cause du décallage d''une 1/2 période)?'
      stop ' en attendant dans kount_to_date.f'
      endif                                                            !20/08/03


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE3D:
! DEBUT:
      if(i2dh.eq.1) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      x4_r8=ki1*dti_fw/86400.+daysim(0)
      i4=int(x4_r8)
      i5=nint(real(dayindex(i4))/10000.)  ! an
      i6=nint(real(dayindex(i4)-i5*10000)/100.) ! mois
      i7=dayindex(i4)-i5*10000-i6*100

                                                                       !11/01/05
! Exprimer le residus journalier en heures:
      x3_r8=(x4_r8-real(i4))*24.
      i3=int(x3_r8)

                                                                       !11/01/05
! Exprimer le residus horaire en minutes:
      x2_r8=(x3_r8-real(i3))*60.
      i2=int(x2_r8)

                                                                       !11/01/05
! Exprimer le residus minute en seconde:
      x1_r8=(x2_r8-real(i2))*60.
      i1=int(x1_r8)

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

!     I1=INT( AMOD(KI1*0.5*DTE,60.)       )
!     I2=INT( AMOD(KI1*0.5*DTE/60.,60.)   )
!     I3=INT( AMOD(KI1*0.5*DTE/3600.,24.) )

      i1=      mod(int(ki1*0.5*dte_lp),60 )                            !19/08/03
      i2=      mod(int(ki1*0.5*dte_lp/60.),60 )
      i3=      mod(int(ki1*0.5*dte_lp/3600.),24 )

       x4_r8=ki1*0.5*dte_lp/3600./24.+daysim(0)
      i4=int(x4_r8)
      i5=nint(real(dayindex(i4))/10000.)  ! an
      i6=nint(real(dayindex(i4)-i5*10000)/100.) ! mois
      i7=dayindex(i4)-i5*10000-i6*100

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE2D:
! FIN.
      return
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine kount_to_date
