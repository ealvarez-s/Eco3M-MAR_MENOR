










! ATTENTION NOUVELLE ROUTINE - NOUVEAU NOM !
      subroutine elapsedtime2date(time_,year_,month_,day_,hour_,minute_,second_)
!______________________________________________________________________
! S model
! release S.26 - last update: 01-04-13
!______________________________________________________________________

      use module_principal
      implicit none
      double precision time_
      integer year_,month_,day_,hour_,minute_,second_

! Converts iteration into calendar time

!..............................................................................
! Version date      Description des modifications
! 2010.19 03-04-10  Mise en service
! S.26    01-04-13  modif des noms des variables et passage en argument
!..............................................................................

      if(time_<0.)then
      write(6,*)
      write(6,*)'anomalie reperée dans elapsedtimetodate.f'
      write(6,*)'recherche de date pour iteration negative.'
      write(6,*)'en effet on passe ',time_,' en argument.'
      write(6,*)'suggestions:'
      write(6,*)'verifier notebook_time:'
      write(6,*)'  date de debut anterieure à date reference ?'
      write(6,*)'verifier notebook_obcforcing:'
      write(6,*)'  date de premiere sortie (menu imbriqué) prématurée'
      write(6,*)' (notament à cause du décallage d''une 1/2 période)?'
      stop ' en attendant dans elapsedtimetodate.f'
      endif                                                            !20/08/03


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE3D:
! DEBUT:
      if(i2dh.eq.1) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     x4_r8=ksecond_*dti_fw/86400.+daysim(0)
      x4_r8=time_/86400.+daysim(0)
      i4=int(x4_r8)
      year_=nint(real(dayindex(i4))/10000.)  ! an
      month_=nint(real(dayindex(i4)-year_*10000)/100.) ! mois
      day_=dayindex(i4)-year_*10000-month_*100

                                                                       !11/01/05
! Exprimer le residus journalier en heures:
      x3_r8=(x4_r8-real(i4))*24.
      hour_=int(x3_r8)

                                                                       !11/01/05
! Exprimer le residus horaire en minutes:
      x2_r8=(x3_r8-real(hour_))*60.
      minute_=int(x2_r8)

                                                                       !11/01/05
! Exprimer le residus minute en seconde:
      x1_r8=(x2_r8-real(minute_))*60.
      second_=int(x1_r8)

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

      stop 'elapsedtimetodate minute_dh==0 invalide'

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE2D:
! FIN.
      return
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine elapsedtime2date

!........................................................................
! ANCIENNE ROUTINE MAINTENU POUR COMPATIBILITE AVEC CLAUDE ET CAROLINE
      subroutine elapsedtimetodate(time_loc)
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 03-04-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
      double precision time_loc

! Converts iteration into calendar time

!..............................................................................
! Version date      Description des modifications
! 2010.19 03-04-10  Mise en service
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

      if(time_loc<0.)then
      write(6,*)
      write(6,*)'anomalie reperée dans elapsedtimetodate.f'
      write(6,*)'recherche de date pour iteration negative.'
      write(6,*)'en effet on passe ',time_loc,' en argument.'
      write(6,*)'suggestions:'
      write(6,*)'verifier notebook_time:'
      write(6,*)'  date de debut anterieure à date reference ?'
      write(6,*)'verifier notebook_obcforcing:'
      write(6,*)'  date de premiere sortie (menu imbriqué) prématurée'
      write(6,*)' (notament à cause du décallage d''une 1/2 période)?'
      stop ' en attendant dans elapsedtimetodate.f'
      endif                                                            !20/08/03


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE3D:
! DEBUT:
      if(i2dh.eq.1) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     x4_r8=ki1*dti_fw/86400.+daysim(0)
      x4_r8=time_loc/86400.+daysim(0)
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

      stop 'elapsedtimetodate i2dh==0 invalide'

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     CAS DU MODELE2D:
! FIN.
      return
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine elapsedtimetodate
