










      subroutine graph1d_out
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 15-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none

!...............................................................................
! Version Date      Description des modifications
!         27/01/03: amenagements pour nouvelle facon de faire les sorties
!                   graphique. La construction du nom du fichier graphique
!                   est changée
!         19/01/07: nom du repertoire d'ecriture lu dans notebook_graph
! 2010.6  11-02-10  notebook_graph = nomfichier(21)
! 2010.20 15-04-11  Calculs sur la base d'un temps en secondes
!...............................................................................

!..............................................................................
! Lecture de notebook_graph
! Début:
!..............................................................................
      open(unit=3,file=nomfichier(21)) ! lecture de notebook_graph !11-02-10
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,'(a)')dirgraph                                            !13/04/06
      close(3)
!..............................................................................
! Lecture de notebook_graph
! Fin.
!..............................................................................

!                            /   /   /

!..............................................................................
! Nom du repertoire des sorties graphiques
! Debut:
!..............................................................................
      if(dirgraph(1:1).eq.'*')then !>>>>>>>>>>>>>>
! Si ancien notebook alors par default c'est ../GRAPHIQUES/
       lname4=14
       dirgraph(1:lname4)='../graphiques/'
      else                         !>>>>>>>>>>>>>>
       do k=1,90
         if(dirgraph(k:k).eq.' ') then !-------->
          lname4=k-1
          goto 21
         endif                         !-------->
       enddo
   21  continue
       if(dirgraph(lname4:lname4).eq.'/')lname4=lname4-1
       lname4=lname4+1
       dirgraph(lname4:lname4)='/'
      endif                        !>>>>>>>>>>>>>>
!..............................................................................
! Nom du repertoire des sorties graphiques
! Fin.
!..............................................................................

! calcul de la date                                                    !27/01/03
!     if(i2dh.eq.1)call kount_to_date(kount)
      if(i2dh.eq.1)call kount_to_date(iteration3d)

      write(texte80(1),                                                 &
       '(1x,i2,1x,a9,1x,i4,1x,a,i2,a1,i2,a1,i2)')                       &
       i7,month(i6),i5,'h:m:s ',i3,':',i2,':',i1

! puis afficher la date à l'écran:                                     !27/01/03
      write(6,*)'--------------------------------'
      write(6,*)'fichier graphique date:'
      write(6,'(a33)')texte80(1)(1:33)

! écrire année année mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

! écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
! puis on efface le "1" avec un caractere de séparation, "_"
! écrire ":" dans TEXTE90:
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04

      lrec=(imax+2)*(jmax+2)*4

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CONSTRUIRE LE NOM DU FICHIER DE SORTIE GRAPHIQUE                     !27/01/03
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! NOM DE FICHIER SORTIE PAR DEFAUT:                                    !27/01/03
      texte250=dirgraph(1:lname4)//'3dgraph.out'                       !13/04/06

!*******************************************************************************
! CAS DES SORTIES DU MODELE AUX 8 DATES SPECIFIEES DANS notebook_time
! DEBUT:
      if(idate_output.eq.0) then                                       !27/01/03
!*******************************************************************************

      if(iteration3d.eq.kount1)then                                          !12/05/05
      write(texte30(1:6),'(a6)')'_1.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount2)then
      write(texte30(1:6),'(a6)')'_2.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount3)then
      write(texte30(1:6),'(a6)')'_3.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount4)then
      write(texte30(1:6),'(a6)')'_4.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount5)then
      write(texte30(1:6),'(a6)')'_5.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount6)then
      write(texte30(1:6),'(a6)')'_6.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount7)then
      write(texte30(1:6),'(a6)')'_7.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif
      if(iteration3d.eq.kount8)then
      write(texte30(1:6),'(a6)')'_8.out'
      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30(1:6)         !13/04/06
      endif

!*******************************************************************************
! CAS DES SORTIES DU MODELE AUX 8 DATES SPECIFIEES DANS notebook_time
! FIN.
      endif
!*******************************************************************************


!                            /   /   /


!*******************************************************************************
! CAS DES SORTIES DU MODELE REGULIERES
! DEBUT:
      if(idate_output.eq.1) then                                       !27/01/03
!*******************************************************************************

!     if(int(real(kount  )/kountgraph)                                  &
!       -int(real(kount-1)/kountgraph).eq.1) then !>>>>>>>>>>>>>>>>>>>

      if(int(elapsedtime_now/graphperiod)                              &
        -int(elapsedtime_bef/graphperiod)==1) then !>>>>>>>>>>>>>>>>>>>   !15-04-11

!     k=int( real(kount) / kountgraph )
      k=int( elapsedtime_now/ graphperiod )                    !15-04-11

        if(k.le.9)then
           write(texte30,'(a,i1,a)')'_',k,'.out'
        else
         if(k.le.99)then
           write(texte30,'(a,i2,a)')'_',k,'.out'
         else
          if(k.le.999)then
           write(texte30,'(a,i3,a)')'_',k,'.out'
          else
           write(texte30,'(a,i4,a)')'_',k,'.out'
          endif
         endif
        endif

      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30              !13/04/06

      endif                                        !>>>>>>>>>>>>>>>>>>>

!*******************************************************************************
! FAIRE DES SORTIES DU MODELE REGULIERES
! FIN.
      endif
!*******************************************************************************


!                            /   /   /


!*******************************************************************************
! FAIRE DES SORTIES DU MODELE A DES DATES SPECIFIEES DANS notebook_output
! DEBUT:
      if(idate_output.eq.2) then                                       !27/01/03
!*******************************************************************************

      do k=1,kmax_dof
       if(     elapsedtime_now> tdate_output(k)     &
          .and.elapsedtime_bef<=tdate_output(k)) then !§§§§§§§§§>      !15-04-11

        if(k.le.9)then
           write(texte30,'(a,i1,a)')'_',k,'.out'
        else
         if(k.le.99)then
           write(texte30,'(a,i2,a)')'_',k,'.out'
         else
          if(k.le.999)then
           write(texte30,'(a,i3,a)')'_',k,'.out'
          else
           write(texte30,'(a,i4,a)')'_',k,'.out'
          endif
         endif
        endif

      texte250=dirgraph(1:lname4)//texte90(1:15)//texte30              !13/04/06

       endif                               !>>>>>>>>>>>>>>>>>>>>
      enddo

!*******************************************************************************
! FAIRE DES SORTIES DU MODELE A DES DATES SPECIFIEES DANS notebook_output
! FIN.
      endif
!*******************************************************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CONSTRUIRE LE NOM DU FICHIER DE SORTIE GRAPHIQUE                     !27/01/03
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      write(6,'(a)')texte250

      open(unit=12,file=texte250)

      i=2
      j=2
      do k=1,kmax+1
      k1=min0(k,kmax)
      write(12,'(10(e12.5,1x))')                                        &
                     depth_t(i,j,k1)                                    &
                   , depth_w(i,j,k)                                     &
                   ,   tem_t(i,j,k1,1)                                  &
                   ,   sal_t(i,j,k1,1)                                  &
                   ,   vel_u(i,j,k1,1)                                  &
                   ,   vel_v(i,j,k1,1)                                  &
                   ,    km_w(i,j,k1)                                    &
                   ,  tken_w(i,j,k)                                     &
                   ,  tkll_w(i,j,k)                                     &
                   ,  tkle_w(i,j,k)
      enddo

      close(12)

      end subroutine graph1d_out
