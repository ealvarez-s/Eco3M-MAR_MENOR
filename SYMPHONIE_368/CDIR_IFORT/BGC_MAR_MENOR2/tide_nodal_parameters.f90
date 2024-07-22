










      subroutine tide_nodal_parameters
!______________________________________________________________________
! S model
! release S26.1 - last update: 07-02-13
!______________________________________________________________________

      use module_principal
      use module_parallele
      use module_systeme
      implicit none
      double precision v0_loc                                          !24-12-09
      integer debug1_loc,debug2_loc


!...............................................................................
! version date      description des modifications
! 2010.2  22-12-09  mise en service
!         24-12-09  variable locale: v0_loc
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.20 15-04-11  Calculs sur la base d'un temps en secondes
! 2010.24 19-12-11  Ne pas permettre d'utiliser un fichier nodal dont
!                   la premiere date est posterieure à la date de demarrage
! S26     07-02-13  seul proc zero écrit à l'écran
!...............................................................................

!..........
! This routine reads the files containing the tidal nodal parameters
      do ktide=1,kmaxtide

       if(par%rank==0)write(6,'(A,A)')'Lecture fichier ',trim(tidenodalfile(ktide))
       open(unit=3,file=tidenodalfile(ktide))

! Debug section: check that all the nodal files come from the same directory
       debug2_loc=0
       if(index(tidenodalfile(ktide),'1950')/=0)debug2_loc=1950
       if(index(tidenodalfile(ktide),'2000')/=0)debug2_loc=2000
       if(debug2_loc==0)stop     &
       'unidentified nodal file detected in tide_nodal_parameters.F90'
       if(ktide==1)debug1_loc=debug2_loc
       if(debug1_loc/=debug2_loc) then
       write(*,*)   &
      'All the nodal files do not come from the same directory'
       stop ' Stop in subroutine tide_nodal_parameters'
       endif
       debug1_loc=debug2_loc

        read(3,*)
        read(3,*)
        read(3,*)

        do loop1=1,99999

          read(3,*,end=130)i1,i2,i3,i4,i5,i6,v0_loc     &
                          ,utide(ktide,2),ftide(ktide,2)

          utide(ktide,2)=utide(ktide,2)*deg2rad

          call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10

          if(elapsedtime_out<=elapsedtime_now) then  !-----------> !16-04-11

! Pas la peine d'interpoler ti0tide et v0tide car il s'agit de 2 parametres
! de synchonisation: en changer brutalement une fois par an ne gêne pas
! la continuité de la solution
          if(iteration3d==kount0) then !\initial\initial>
                v0tide(ktide)=v0_loc*deg2rad
                ti0tide=elapsedtime_out
          endif                        !\initial\initial<

           tidenodal_prev_rdv=elapsedtime_out
           utide(ktide,0)=utide(ktide,2)
           ftide(ktide,0)=ftide(ktide,2)

          else                            !----------->

! Cet algo suppose que la premiere date du fichier nodal soit avant la date de demarrage. !19-12-11
! On arrete le run si cette condition n'est pas remplie:
          if(loop1==1)then                                            !19-12-11
           write(*,*)
           write(*,*)'-------------------------------------------------'
           write(*,*)'ERREUR: first date in file'
           write(*,'(A)')trim(tidenodalfile(ktide))
           write(*,*)'should not be after the departure time of'
           write(*,*)'the simulation. Change for a more suitable file'
           stop ' STOP in  subroutine tide_nodal_parameters!'
          endif

             if(ktide>1)x1=tidenodal_next_rdv
             tidenodal_next_rdv=elapsedtime_out ! time of the next rendez-vous (read a new value)

! pour ne pas interpoler 2 valeurs de part et d'autre de 0° ( 0.5*(0+360)=180 ):
            if(utide(ktide,0) > utide(ktide,2)+pi)utide(ktide,0)=utide(ktide,0)-2.*pi
            if(utide(ktide,0) < utide(ktide,2)-pi)utide(ktide,0)=utide(ktide,0)+2.*pi

            goto 131 ! sortir de la boucle sur loop1

          endif                            !----------->

        enddo ! fin de boucle sur loop1

 130    stop 'fichier de parametres nodaux insuffisant'

 131    close(3)

      enddo ! fin de boucle sur ktide

      end subroutine tide_nodal_parameters
