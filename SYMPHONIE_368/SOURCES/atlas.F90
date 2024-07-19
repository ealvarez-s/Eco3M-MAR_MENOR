      subroutine atlas
!______________________________________________________________________
! SYMPHONIE ocean model
! release 287 - last update: 17-07-20
!______________________________________________________________________

      use module_principal ; use module_s ; use module_parallele ; use module_global
      implicit none
      integer loop_,rivercount_
#ifdef synopsis
       subroutinetitle='atlas'
       subroutinedescription= &
      'Reads notebook_atlas. Computes grid indexes from lon lat ' &
      //'values listed in notebook_atlas'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!...............................................................................
! Version date      Description des modifications
!         06/11/02: mise en service.
!         20/04/09: notebook_atlas devient NOMFICHIER(19)
! 2010.8  09-05-10  seul le proc 0 ecrit interrupteur
! 2010.13 01-11-10  messages à l'ecran
! S.26    28-12-13  grille horizontale generalisee: tous les
!                   proc calculent!
!         04-01-17  Aide au positionnement des fleuves
!         11-06-17  Aide au positionnement des fleuves
!         16-06-17  Routine atlas inutilisable sous certaines conditions
!         14-05-18  algo basE sur glob_mask A la place de mask_t
! v259    26-09-19  correction bug
! v260    13-10-19  plus d'infos dans le fichier tmp/notebook_rivers_missing
!         19-10-19  suite point precedent
! v287    17-07-20  utiliser tmpdirname
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


! Ce programme permet de localiser sur la grille horizontale du model_
! des sites repérés en latitudes et longitudes (décimales) listés dans
! le fichier notebook_atlas

! cette condition vient de l'impossibilite d'utiliser la routine
! latlontoij si la grille de symphonie n'est pas une grille analytique
! et que les domaines inutiles ont ete supprimEs (lon, lat ne sont pas
! definies sur toute la grille):
!     if(lonlatfile/='nofile'.and.nbdom/=nbdom_imax*nbdom_jmax)return !16-06-17

      if(.not.allocated(glob_mask)) &
      stop 'Err 54 atlas .not.allocated(glob_mask)'

      texte90=nomfichier(19) ! 'notebook_atlas'                        !29/04/09
      if(par%rank==0)write(6,*)'Subroutine atlas lecture de ',trim(texte90)                            !01-11-10

      do loop_=1,2 !--loop-->
      rivercount_=0 ! compteur de fleuve
! Avec loop_ on fait 2 tours de boucles. 
! Le premier tour serat A compter le nombre de fleuve afin de faire l'entete de tmp/notebook_rivers
! Le deuxieme tour sert A completer le tmp/notebook_rivers

      open(unit=4,file=texte90) ! texte90 = notebook_atlas
      read(4,*,end=100)

      if(par%rank==0.and.loop_==2) then !00000>
       open(unit=3,file=trim(tmpdirname)//'messages',position='append')
       write(3,*)'-------------------------------------'
       write(3,*)'subroutine atlas:'
      endif                             !00000>

      do 99 k=1,99999

      read(4,'(a)',end=101)texte250
      read(4,*    ,end=101)x1,x2

      latit1=x1*pi/180.
      longi1=x2*pi/180.
      call latlontoij(longi1,latit1,'glb')    !15-03-12

! Chercher le point viable le plus proche pour les fleuves !04-01-17
      i0=nint(deci) ! -par%timax(1)
      j0=nint(decj) ! -par%tjmax(1)
      k0=-1
      ksecu=0

  62  k0=k0+1
      do j=j0-k0,j0+k0
      do i=i0-k0,i0+k0

      if(i>=1.and.i<=iglb.and.j>=1.and.j<=jglb) then !pmxpmx>
       if(glob_mask(i,j)==0) then !>>>
        if(i/=1) then
           if(glob_mask(i-1,j)==1)ksecu=3
        endif
        if(i/=iglb) then
           if(glob_mask(i+1,j)==1)ksecu=1
        endif
        if(j/=jglb) then
           if(glob_mask(i,j+1)==1)ksecu=2
        endif
        if(j/=1) then
           if(glob_mask(i,j-1)==1)ksecu=4
        endif
        if(ksecu/=0)goto 63
       endif                      !>>>
      endif                                          !pmxpmx>

      enddo
      enddo
!     if(ksecu==0.and.k0<5)goto 62 ! si k0>=5 on renonce A trouver le point le plus proche
      if(k0<5)goto 62 ! si k0>=5 on renonce A trouver le point le plus proche

!..................................
  63  continue

      if(ksecu/=0.and.loop_==1)rivercount_=rivercount_+1 ! compteur de fleuves

      if(ksecu==0.and.par%rank==0.and.loop_==2) then ! m(OoO)m >
       k10=s_unit(7)
       open(unit=k10 &
        ,file=trim(tmpdirname)//'notebook_rivers_missing',position='append') !17-07-20
        write(k10,'(a)')'---------------------------'
        write(k10,'(a)')trim(texte250)  
        write(k10,'(a,1x,f10.5,1x,f10.5)')'River coordinates' &
       ,longi1*180/pi,latit1*180/pi !13-10-19
        write(k10,'(a)')'First corresponding grid point:'
        write(k10,*)'deci,decj',deci,decj
        write(k10,*)'i1=nint(deci),j1=nint(decj)',nint(deci),nint(decj)
        if(nint(deci)==-999.or.nint(decj)==-999) then !>>>
         write(k10,*)'-999 values for points outside de model domain'
        endif                                         !>>>
        write(k10,*)'Iterative process stopped after',k0,' loops'
       close(k10)
      endif                                         ! m(OoO)m >

      if(ksecu/=0.and.par%rank==0.and.loop_==2) then ! m[°v°]m >

       k10=s_unit(7)
       open(unit=k10,file='tmp/notebook_rivers',position='append')
        write(k10,'(a)')trim(texte250)
        write(k10,'(i6,a)')i & ! i+par%timax(1) & !26-09-19
        ,' ! river mouth location along i axis'
        write(k10,'(i6,a)')j & ! j+par%tjmax(1) & !26-09-19
        ,' ! river mouth location along j axis'
        write(k10,'(a)') &
        '0      ! length (grid nodes) of an additionnal water way'
        write(k10,'(i6,a)') &
        ksecu,' ! Flux Direction 1=incr i, 2=incr j, 3=decr i, 4=decr j'

! Blabla dans "messages"
       write(3,*)'-----------------------------------------------------'
       write(3,'(a5,a)')'site:',trim(texte250)
       write(3,'(a24,2(f8.3,2x))')'position (lat lon):     '             &
        ,latit1*180./pi                                                  &
        ,longi1*180./pi
       write(3,'(a24,2(f8.3,2x))')'position sur la grille: ',deci,decj
       write(3,*)'Parametrage notebook_rivers dans messages_fleuves'

      endif                             ! m[°v°]m >

      do k9=1,20
       read(4,'(a)',end=55)texte250 
       if(ksecu/=0.and.par%rank==0.and.loop_==2)write(k10,'(a)')trim(texte250)  !11-06-17
       if(texte250(1:5)=='*****')goto 55
      enddo
   55 continue

      if(ksecu/=0.and.par%rank==0.and.loop_==2)close(k10)

!..................................
      

!     if(par%rank==0.and.loop_==2) then !00000>
!      write(3,*)'-----------------------------------------------------'
!      write(3,'(a5,a)')'site:',trim(texte250)
!      write(3,'(a24,2(f8.3,2x))')'position (lat lon):     '             &
!       ,latit1*180./pi                                                  &
!       ,longi1*180./pi
!      write(3,'(a24,2(f8.3,2x))')'position sur la grille: ',deci,decj
!      write(3,*)'Parametrage notebook_rivers dans messages_fleuves'
!     endif                             !00000>

!     do k=1,20
!      read(4,'(a)',end=101)texte250 ; if(par%rank==0)write(3,'(a)')trim(texte250) ; if(texte250(1:5)=='*****')goto 55
!     enddo


  99  continue
  101 continue

      if(par%rank==0.and.loop_==2)close(3)

  100 continue
      close(4)


      if(par%rank==0.and.loop_==1) then !111>
       k10=s_unit(7)
       open(unit=k10,file= &
            trim(tmpdirname)//'notebook_rivers',position='append') !17-07-20
        write(k10,'(i0,a)')rivercount_,'     ! Number of river inputs'
        write(k10,'(a)')'_____________________________________________'
       close(k10)
      endif             !111>

      enddo !--loop-->

      if(par%rank==0)write(6,*)'ok'

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
#endif

      end subroutine atlas
