      subroutine river_upd(ichoix)                                     !09/08/01
!______________________________________________________________________
! SYMPHONIE ocean model
! release 301 - last update: 06-04-21
!______________________________________________________________________

      use module_principal ; use module_parallele ; use module_s
      implicit none
      integer ichoix,loop_          !09/08/01
#ifdef synopsis
       subroutinetitle='river_upd'
       subroutinedescription='Updates the river inputs'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         09/08/01: amenagements sur l'initialisation des fleuves
!         25/11/01: prise en compte des petits fleuves (introduits en surface)
!         22/04/02: debug initialisation
!         10/02/04: debug ouverture fichier
!         06/05/04: ajout d'une option permettant de commencer avant nc=1
!         27/07/04: debug temperature des fleuves "par d�faut"
!                   & on prefere calculer RIVER_T par la formule analytique
!                   � chaque iteration plutot que l'interpolation lin�aire
!                   entre 2 �ch�ances (ceci les anciennes lignes sont simplement
!                   comment�es pour pouvoir reservir un jour)
!         07/03/05: Des ecritures � l'ecran pour montrer que les fichiers sont
!                   lus correctement
!         08/03/05: produire un fichier ascii des debit des fleuves pour verif
!         18/03/05: debug relatif au point precedent....
!         28/04/05: aide au debugage en ajoutant de nouveaux messages � l'ecran
!         30/09/05: modif pour compiler avec -r8
!         17/04/07: Passage � coordonnees curvilignes (voir ajout dx_y et cie...)
!         23/01/08: renouvellement des decharges par modulo d'entier supprim�
!         08/04/08: Si zero riviere alors on sort tout de suite du cas
!                   "initialisation" pour eviter les problemes de date (9 fevrier)
!                   rencontr�s par Francis.
! 2009.3  04-11-09  ecrire les fichiers temporaire dans repertoire tmp
! 2010.4  26-01-10  possibilit� de lire des fichiers de rivieres ascii
! 2010.7  23-02-10  prise en compte de friver
! 2010.8  11-03-10  une rampe sur le debit des fleuves
!         05-05-10  Si plusieurs rivieres partagent le m�me point de grille
!                   alors la riviere portant le numero le plus grand cumule
!                   le debit des autres rivieres
!         09-05-10  - Ecriture d'un fichier binaire si par%rank=0
!                   - Prevoir le cas de valeurs flaguees (-999.)
! 2010.13 13-10-10  Seul le proc 0 ecrit
!         03-11-10  des arguments pass�s dans date_to_kount
! 2010.19 13-04-11  Ne pas faire les listes de fichiers binaires dans le
!                   cas d'un d�bit constant
! 2010.20 15-04-11  Renouvellement des forcages sur la base d'un temps
!                   exprim� en secondes
! 2010.24 03-12-11  Le fichier riviere au format ascii devient la norme
! 2010.25 02-03-12  Ajout d'une barriere mpi
! S26.1   20-09-12  Seul le rank 0 ecrit a l'ecran
!         08-04-14  liste binaire: ecriture partagee entre procs
!         08-09-16  Choisir le mois du minimum de temperature des rivieres
!         13-06-17  autoriser les debits negatifs (cas maree)
!         01-02-18  modif message A l'ecran
! v251    09-04-19  ajout resurgences
! v272    14-01-20  Seul par%rank=0 lit le fichier puis partage avec les autres proc 
! v296    18-02-21  if(riverdir(k1)==0)  riverflux(k1,k3)=-discharge !18-02-21 ! par la surface
!         27-02-21  Meme si le debit est academique constant, on applique la rampe transitoire
! v301    06-04-21  trim(tmpdirname)//riverfile(loop_)(1:len_)//'.binrec' !06-04-21
!...............................................................................
!    _________                    .__                  .__             ! (�o�) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! INITIALISATION DES FLEUVES. PARTIE AJOUTEE LE 09/08/01
! DEBUT
      if(ichoix.eq.1) then !************>                              !09/08/01
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(nriver.eq.0)return                                            !08/04/08

      call make_a_binary_file                                          !26-01-10

! Determiner l'iteration (decimale) correspondant � la temperature     !05/05/04
! minimale des fleuves: 
      do loop_=1,nriver
! Avant cette ligne river_timeref est le numero du mois donnE dans ! notebook_river
       call datetokount(datesim(1,1),nint(river_timeref(loop_)),1,0,0,0) ! a/m/j/h/m/s !08-09-16
       river_timeref(loop_)=elapsedtime_out
! Apres cette ligne river_timeref est un temps en seconde depuis l'annee
! initiale, le mois donnE dans notebook_river
      enddo

!     const1=dti_fw/60./60./24. !expression du pas temps en jours      !05/05/04
!     const1=2.*pi/86400./365.       ! Pour passer de secondes en jours
      const1=2.*pi/86400./365.25     ! Pour passer de secondes en jours!08-09-16
      do 999 k1=1,nriver

! Temp�rature instann�e:
      river_t(k1,1)=                                                    &
        (river_tmin(k1)+river_tmax(k1))/2.                              &
       -(river_tmax(k1)-river_tmin(k1))/2.                              &
       *cos(const1*(elapsedtime_now-river_timeref(k1)))            !15-04-11

      if(realriver(k1).eq.1) then !>>>>>>>>>>>>>> ! test d�plac� le    !27/07/04

! ajout le:                                                            !23/01/08
      i1=dateriver(1,k1)
      i2=dateriver(2,k1)
      i3=dateriver(3,k1)
      i4=dateriver(4,k1)
      i5=dateriver(5,k1)
      i6=dateriver(6,k1)
!     call date_to_kount
      call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
      riverdt(k1,2)=elapsedtime_out          !15-04-11
      riverdt(k1,1)=riverinfo(k1)*3600.      !15-04-11


! comment�es le                                                        !23/01/08
!     KOUNTMOD=RIVERDT(K1,1) ! echeance
!     K2      =RIVERDT(K1,2) ! repere NC=1

        if(par%rank==0)write(6,*)'sur le point de lire:'                              !07/03/05
        if(par%rank==0)write(6,'(a60)')riverfile(k1)
        lrec=4
        if(par%rank==0)                 & !14-01-20
        open(unit=3,file=riverfile(k1)  &
                   ,access='direct'     &
                   ,form='unformatted'  & !10/02/04
                   ,recl=lrec)
! NC   enregistrement present

        do 888 k3=0,2,2

! comment�es le:                                                       !23/01/08
!       NC=INT(1.+REAL(KOUNT-K2)/REAL(KOUNTMOD)+K3/2)
!       IF(MOD(KOUNT-K2,KOUNTMOD).EQ.0.AND.K3.EQ.2)NC=MAX0(NC-1,1)     !22/04/02
! Remplac�es par: le:                                                  !23/01/08
        x2=(elapsedtime_now-riverdt(k1,2))/riverdt(k1,1)               !18-04-11
        x1=(elapsedtime_bef-riverdt(k1,2))/riverdt(k1,1)
        j2=int(x2)
        j1=int(x1)
        if(x2.lt.0.)j2=j2-1
        if(x1.lt.0.)j1=j1-1
        nc=int(1.+x2)+k3/2
! cas particulier qui anticipe sur la procedure hors initialisation
! ichoix=2: (attention il faut prevoir que NC=0 n'existe pas)
        if(k3.eq.2.and.(j2-j1).eq.1)nc=nc-1
        nc=max0(nc,ncmin_river)                                        !06/05/04

        if(nc.le.0)then !---- debug ----->                             !06/05/04
         write(6,*)'je suis bloque dans river_upd choix 1'             !01-02-18
         write(6,*)'fleuve No',k1
         write(6,*)'car nc=',nc,' est <= 0 ce qui signifie que je suis'
         write(6,*)'en avance sur la premiere echeance disponible.'
         write(6,*)'pour continuer je fixe arbitrairement nc � 1'
         ncmin_river=1                                                 !06/05/04
         nc=max0(nc,ncmin_river)                                       !06/05/04
         stop 'stop river_upd' !09-08-14
        endif           !---- debug ----->


        if(par%rank==0)write(6,*)'.....................................'!22/04/02
        if(par%rank==0)write(6,*)'nc entier et decimal: ',nc,1.+x2+k3/2

! Seul par%rank=0 lit le fichier puis partage avec les autres proc !14-01-20
        if(par%rank==0)read(3,rec=nc)discharge                     !14-01-20
        s_typevar=0                                                !14-01-20
        call donne_le_type(discharge)                              !14-01-20
        if(s_typevar==4)call mpi_bcast(discharge,1,mpi_real,0,par%comm2d,ierr)             !14-01-20
        if(s_typevar==8)call mpi_bcast(discharge,1,mpi_double_precision,0,par%comm2d,ierr) !14-01-20
        if(s_typevar==0)stop 's_typevar non prevu. Completer interface'                    !14-01-20 

        discharge=discharge/friver(k1)                                 !23-02-10


        if(par%rank==0)write(6,*)'debit en m3/s en nc: ',discharge*friver(k1)

! riverflux converti en transport (avec le bon signe)
      if(riverdir(k1).eq.1) riverflux(k1,k3)= discharge
      if(riverdir(k1).eq.2) riverflux(k1,k3)= discharge
      if(riverdir(k1).eq.3) riverflux(k1,k3)=-discharge
      if(riverdir(k1).eq.4) riverflux(k1,k3)=-discharge
      if(riverdir(k1)==0)   riverflux(k1,k3)= discharge !18-02-21 ! par la surface
      if(riverdir(k1)==-1)  riverflux(k1,k3)= discharge !09-04-19 ! par le fond

! Temp�rature �ch�ance:
!cccccK10=(NC-1)*KOUNTMOD+K2 ! iteration correspondant � l'echeance    !23/01/08
      x10=(nc-1)*riverdt(k1,1)+riverdt(k1,2)                   !15-04-11
      river_t(k1,k3)=                                                 &
        (river_tmin(k1)+river_tmax(k1))/2.                            &
       -(river_tmax(k1)-river_tmin(k1))/2.                            &
       *cos(const1*(x10-river_timeref(k1)))                        !15-04-11

 888    continue

        if(par%rank==0)close(3) !14-01-20
        if(par%rank==0)write(6,*)'lecture ok'                                          !07/03/05

      endif                       !>>>>>>>>>>>>>>
  999 continue
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)     !13-10-10
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! INITIALISATION DES FLEUVES. PARTIE AJOUTEE LE 09/08/01
! FIN
      return
      endif                !************>                              !09/08/01
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



!                          / / /



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise � jour des �ch�ances:
! DEBUT:
      if(ichoix.eq.2) then !%%%%%%%%%%%%>                              !09/08/01
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do 1000 k1=1,nriver
      if(realriver(k1).eq.1) then

      call time_to_update_forcing_file(riverdt(k1,2),riverdt(k1,1))  !15-04-11
      x2=(elapsedtime_now-riverdt(k1,2))/riverdt(k1,1)               !15-04-11

      if(decision==1)  then                           !*******************> !15-04-11

! NC est le record au temps present ! NC+1 record pour echeance suivante:
      nc=max0(int(1.+x2),ncmin_river-1)

        lrec=4
        if(par%rank==0)                & !14-01-20
        open(unit=3,file=riverfile(k1) &
                   ,access='direct'    &
                   ,form='unformatted' & !10/02/04
                   ,recl=lrec)

        if(nc+1.le.0)then !---- debug ----->                           !06/05/04
        write(6,*)'je suis bloqu� dans river_upd choix 2'
        write(6,*)'fleuve n�',k1
        write(6,*)'car nc+1=',nc+1,' est <= 0 qui signifie que je suis'
        write(6,*)'en avance sur la premiere echeance disponible.'
        ncmin_river=1                                                 !06/05/04
        nc=max0(nc,ncmin_river-1)                                     !23/01/08
        stop !09-08-14
        endif             !---- debug ----->

! NC+1 enregistrement pour la prochaine echeance
        if(par%rank==0)read(3,rec=nc+1)discharge                   !14-01-20
        s_typevar=0                                                !14-01-20
        call donne_le_type(discharge)                              !14-01-20
        if(s_typevar==4)call mpi_bcast(discharge,1,mpi_real,0,par%comm2d,ierr)             !14-01-20
        if(s_typevar==8)call mpi_bcast(discharge,1,mpi_double_precision,0,par%comm2d,ierr) !14-01-20
        if(s_typevar==0)stop 's_typevar non prevu. Completer interface'                    !14-01-20 
        discharge=discharge/friver(k1)                                 !23-02-10
        if(par%rank==0)close(3) !14-01-20

        if(par%rank==0)write(6,*)'.....................................'!22/04/02
        if(par%rank==0)write(6,'(a60)')riverfile(k1)
        if(par%rank==0)write(6,*)'nc entier et decimal: ',nc,1.+x2                    !23/01/08
        if(par%rank==0)write(6,*)'debit en m3/s en nc+1: ',discharge*friver(k1)

      riverflux(k1,0)=riverflux(k1,2)
        river_t(k1,0)=  river_t(k1,2)

! riverflux converti en transport (avec le bon signe)
      if(riverdir(k1).eq.1) riverflux(k1,2)= discharge
      if(riverdir(k1).eq.2) riverflux(k1,2)= discharge
      if(riverdir(k1).eq.3) riverflux(k1,2)=-discharge
      if(riverdir(k1).eq.4) riverflux(k1,2)=-discharge
      if(riverdir(k1)==0)  riverflux(k1,2)= discharge                 !25/11/01
      if(riverdir(k1)==-1) riverflux(k1,2)= discharge  !09-04-19

      endif                                             !*******************>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise � jour des �ch�ances:
! FIN.
      endif
 1000 continue
      endif                !%%%%%%%%%%%%>                              !09/08/01
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



!                 ///         ///         ///



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 �cheances
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      const2=2.*pi/86400./365.         !15-04-11
      do 2000 k1=1,nriver

      const3=cos(const2*(elapsedtime_now-river_timeref(k1)))

      if(realriver(k1).eq.1) then !--------------------------------->

      rap=   (elapsedtime_now-riverdt(k1,2))/riverdt(k1,1)            & !23/01/08
       -int( (elapsedtime_now-riverdt(k1,2))/riverdt(k1,1) )


      riverflux(k1,1)=((1.-rap)*riverflux(k1,0)+rap*riverflux(k1,2))*rampe !11-03-10

      else                        !--------------------------------->

! Meme si le debit est academique constant, on applique la rampe transitoire !27-02-21
       riverflux(k1,1)=riverflux(k1,2)*rampe !27-02-21

      endif                       !--------------------------------->

! Temp�rature instann�e:                                               !27/07/04
      river_t(k1,1)=                                                    &
        (river_tmin(k1)+river_tmax(k1))/2.                              &
       -(river_tmax(k1)-river_tmin(k1))/2.                              &
       *const3


 2000 continue

!05-05-10  Si plusieurs rivieres partagent le m�me point de grille
!          alors la riviere portant le numero le plus grand cumule
!          le debit des autres rivieres:

      do k1=1,nriver
       if(river_no(k1)/=k1) then !ififififif>
          riverflux(river_no(k1),1)       &
         =riverflux(river_no(k1),1)       &
         +riverflux(k1,1)
          riverflux(k1,1)=0.
       endif                     !ififififif>
      enddo


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 �cheances
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      return
      end subroutine river_upd



!______________________________________________________________________________


      subroutine make_a_binary_file                                     !26-01-10
      use module_principal
      use module_parallele !#mpi
      implicit none
      real*4             &
        discharge_
      integer            &
        loop_            &
       ,len_             &
       ,secu_            &
       ,k_               &
       ,i_               &
       ,nc_
#ifdef synopsis
       subroutinetitle='make_a_binary_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!Cette routine lit le fichier ascii de riviere et ecrit un fichier binaire
!temporaire � acces direct.

      do loop_   =1,nriver
      if(nint(realriver(loop_   ))==1) then !111111111111111111>   !13-04-11

      secu_   =0
      len_   =len_trim( riverfile(loop_   ) )

      if(riverfile(loop_   )(len_   -5:len_   )=='binrec')  then !binbinbin>
      secu_   =1
      if(par%rank==0)write(6,*)'Les fichiers de riviere sont normalement ascii' !03-12-11
      if(par%rank==0)write(6,*)'Pour passer des fichiers binaires enlever'
      if(par%rank==0)write(6,*)'le stop dans river_upd.F90'
      stop 'STOP dans river_upd.F90' ! enlever cette ligne pour passer des fichiers rivieres format binaire
      endif                                                      !binbinbin>


      if(secu_   ==0) then !----------------------->

      texte250=riverfile(loop_   )
      do k_   =len_   ,1,-1
       if(texte250(k_   :k_   )=='/') then !:::::::::>
              riverfile(loop_   )=  &
!             'tmp'//riverfile(loop_   )(k_   :len_   )//'.binrec'
        trim(tmpdirname)//riverfile(loop_   )(k_   :len_   )//'.binrec' !06-04-21
        goto 488
       endif                               !:::::::::>
      enddo
      riverfile(loop_)= & ! Si on passe par ici c'est que le nom ne contient pas '/'
      trim(tmpdirname)//riverfile(loop_)(1:len_)//'.binrec' !06-04-21
      if(par%rank==0)write(6,*)nriver,loop_   ,secu_   ,trim(texte250),realriver(loop_   )
      goto 488
      stop 'river_upd erreur 471'
  488 continue


      k10=mod(loop_,nbdom) ! numero du rank qui va faire le fichier binaire !08-04-14
      if(par%rank==k10) then !#mpi>>>>>>>>>                                 !08-04-14
!     if(par%rank==0) i then !#mpi>>>>>>>>>      !09-05-10
        nc_   =0
        open(unit=4,file=texte250)
        open(unit=3,file=riverfile(loop_   )      &
                   ,access='direct'               &
                   ,form='unformatted'            &
                   ,recl=4)

!492    read(4,*,end=489)i_   ,discharge_
 492    read(4,*,end=489)      discharge_
!       discharge_   =max(discharge_   ,0.)!09-05-10 commentE le !13-06-17
        nc_   =nc_   +1
        write(3,rec=nc_   )discharge_
        goto 492

 489    close(4)
        close(3)
      endif                !#mpi>>>>>>>>>


      endif                !----------------------->
      endif                                 !111111111111111111>   !13-04-11
      enddo

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)     !02-03-12
#endif

      end subroutine make_a_binary_file
