      subroutine set_rivers(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 363 - last update: 11-01-23
!______________________________________________________________________

      use module_principal ; use module_parallele ; use module_global ; use module_s
      implicit none
      integer case_,unit_
#ifdef synopsis
       subroutinetitle='set_rivers'
       subroutinedescription='Initial procedure for river inputs'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         03/09/01: Essai d'un canal en forme d'embouchure de "trompette"
!         30/10/01: apparition case_=3 initialisation T et S aux sources
!         25/11/01: organisation plus claire du notebook
!                   + apparition d'un nouveau type de fleuve: les petits fleuves
!                   qui ne sont pas introduits lateralement par un canal mais
!                   * la surface d'une colonne. Leur debit est si faible qu'ils
!                   ne generent pas de dynamique de panache (la maille est trop
!                   grande pour cela). Leur influence sur le champ de salinite
!                   intervient via un terme de flux de surface qui va legerement
!                   dessaller la couche de surface. RIVERDIR (=0) permet de les
!                   distinguer des autres.
!         28/03/02: Filiere fleuves introduits * la surface supprim*e
!         21/08/03: Initialisation de T S et ZTA * la cr*ation d'un canal (permet
!                   d'ajouter des canaux apres une analyse sans canaux)
!                   + securit*: on ne modifie pas la bathy lorsque l'on
!                     ouvre un canal si jamais on est dans le cas IHYBSIG=2
!                     c.a.d. que la grille hybride est lue dans un fichier,
!                     autrement dit on se contente de la bathy "sous le masque",
!                     donc attention * ce que celle ci soit convenable!
!         03/09/03: amelioration du message "creuser un canal"
!         20/02/04: amelioration initialisation de T & S dans le canal: on prend
!                   en compte la profondeur
!         05/05/04: temperature des fleuves evolutive
!         28/06/04: correction message * l'ecran
!         27/07/04: debugueur pour l'initialisation de T et S dans les canaux:
!         04/02/05: un message indique que la lecture du notebook_rivers se
!                   passe bien.
!         24/04/05: on autorise finalement la bathy d'etre modifiee dans les
!                   canaux m*me quand on lit un fichier de grille hybride.
!         12/11/05: Dans le cas o* les fleuves sont * priori d*j* creus*s
!                   on cherche tout de m*me * savoir o* se situe l'embouchure
!                   afin plus tard d'appliquer un schema d'advection particulier
!                   dans les cas. L'exercice consiste donc * retrouver des
!                   coordonn*es iriver(kr,3) jriver(kr,3)
!         04/08/06: modif concernant l'initialisatio de T et S dans les canaux
!                   des fleuves creus*s en debut de run. Le choix 3 (mis en
!                   place d'un coin de sel) est remplac* par une interpolation
!                   (verticale) du profil se situant juste devant l'embouchure.
!                   L'ex choix 3 passe en choix 4
!         17/04/07: Passage * coordonnees curvilignes (voir ajout dx_y et cie...)
!         29-04-09  - La longueur des canaux peut etre limitee par le decoupage
!                   en sous domaines en calcul parallele
!                   - Parallelisation
!         29-05-09  Parallelisation secteur ICHOIX=3
! 2009.3  05-10-09  ajout d'un "ifdef parallele"
!         11-10-09  tem_c et sal_c remplacent thz_c et shz_c
! 2010.7  22-02-10  - reparation d'un nom: h_fiver -> hriver
!                   - Si un fleuve porte X fois le m*me nom son d*bit est divis*
!                     par X
! 2010.8  21-03-10  ajout hssh_w
!         23-03-10  les canaux doivent *tre creus*s en 3D
!         26-03-10  Ajout de la parallelisation pour la construction de H dans les
!                   canaux
!         05-05-10  - Si plusieurs fleuves partagent un m*me point d'embouchure alors
!                   on cumule les flux -> creation river_no
!                   - Parallelisation: ne pas permettre qu'une embouchure soit
!                   attribuee * un domaine si elle doit se trouver sur une C.L.
!                   car forcement il existe un domaine o* elle appartient en
!                   tant que point interieur. (cette mesure permet d'envisager
!                   une C.L. de type gradient aux embouchure).
!         19-05-10  Passer outre certains tests quand la longueur du canal
!                   est nulle
! 2010.10 13-06-10  suppression ssh_ext_w
! 2010.13 29-10-10  mpi: upgrade decision d'attribution d'une riviere * son proc
!         01-11-10  river_inout remplace river_dom
! 2010.14 09-11-10  upwindriver_t est un tableau qui permet d'appliquer un
!                   schema 100% upwind * proximit* de l'embouchure des fleuves
!         24-11-10  riverupwdist lu dans notebook_river indique la distance (en
!                   points de grille depuis l'embouchure) d'application du
!                   schema d'advection upwind
! 2010.22 10-05-11  Application du facteur friver au cas des d*bits constants
! 2010.25 08-06-12  Elargissement de la zone couverte par river_inout
! S25.4   22-06-12  debug mpi sur l_river
! S26.1   07-04-13  Puisque dans certains cas (nemo offline par ex) on ne passe
!                   pas par set_rivers(1) on cree le cas 0 qui lit le notebook_rivers
!                   appele depuis set_parameters
!         02-04-14  debug mpi
!         19-04-14  rivertrc_inout rivervel_inout remplacent river_inout
!         27-06-14  ichoix devient case_ & ajout d'un commentaire pour case_=1
!         31-10-14  plus de message a l'ecran
!         21-01-15  debug du test de debug (oui je sais c'est c...)
!         08-09-16  choisir le mois de la temperature min des fleuves
!         25-02-17  upwindriver_t=0 <==> schema d'advection upwind 
!         07-03-17  mise A jour canaux
!         16-03-17  river_obctype choix du scheme de C.L. au point d'entree des rivieres
!         26-04-17  ajout de fonctionnalitEs (etat initial, z0, dx,dy) pour les reservoirs
!         01-05-17  debug point h amont reservoir
!         19-06-17  stop par flag_stop
!         02-05-18  Appliquer une condition (de gradient nul) sur la bathy du point source
!                   meme dans le cas oU l_river=0
!         24-05-18  debug du point 02-05-18
!         02-06-18  suppression ligne relative A iadvect_ts
!         11-11-18  amenagement canaux
!         06-12-18  faire des fichiers de canaux
!         07-12-18  relire un fichier canal
!         12-12-18  relire les lon, lat du fichier canal et en deduire gridrotcos gridrotsin
! v251    09-04-19  ajout resurgences
!         10-04-19  flag_groundwater=1 indique la presence de sources sous marine
! v252    14-04-19  amelioration des comptes rendus d'erreur 
! v253    01-05-19  Pas de schema upwind pour riverdir=-1 !01-05-19
!         06-05-19  if(riverdir(k)==-1)flag_groundwater=1 !06-05-19
! v271    13-12-19  z0b_rivers dans canaux des fleuves
! v296    18-02-21  rivieres de surface: !18-02-21
! v301    06-04-21  if(riverdir(k)==3)riverflux(k,1)=-riverflux(k,1) !06-04-21
!                   if(riverdir(k)==4)riverflux(k,1)=-riverflux(k,1) !06-04-21
! v363    11-01-23  debug longueur de boucle (+1) canaux
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      flag_stop=0

! lecture notebook_rivers
      open(unit=3,file=nomfichier(4))
      read(3,*)nriver
      if(nriver.gt.dim_river) then
      write(6,'(a,a20,a,i4)')'dans ',nomfichier(4),' nriver=',nriver
      write(6,*)'dans parameter dim_river=',dim_river
      write(6,*)'ce qui est insuffisant. il faudrait au moins ',nriver
      flag_stop=1 !19-06-17
      endif
      close(3)

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ; if(k0/=0)stop ' Stop 117 in set_rivers subroutine' !19-06-17
!     WRITE(6,*)DOM_C,PAR%TMECO(1),PAR%TNECO(1)

! rdv au bord du fleuve:
      if(nriver.ge.1) then

!**********************************************************************
      if(case_==0) then    !07-04-13
!......................................................................
! lecture notebook_rivers
      flag_groundwater=0 ! si 1 indique la presence de sources sous marine. 0 sinon. !06-05-19
      flag_surfriver=0   ! si 1 indique la presence de riviere de surface.  0 sinon. !18-02-21
      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(4)
      open(unit=3,file=nomfichier(4))
      read(3,*)nriver

      do 115 k=1,nriver
      realriver(k)=0
      riverinfo(k)=0.

      read(3,*) ! ligne pour rien
      read(3,'(a30)')rivername(k)
      read(3,*)iriver(k,3)                                             !25/11/01
      read(3,*)jriver(k,3)
      read(3,*)l_river(k)
      read(3,*)riverdir(k)
      k0=0
      if(riverdir(k)==-1)flag_groundwater=1 !06-05-19
      if(riverdir(k)==0)flag_surfriver=1 !18-02-21

      if(riverdir(k)==1)iriver(k,3)=iriver(k,3)-l_river(k) !07-03-17
      if(riverdir(k)==3)iriver(k,3)=iriver(k,3)+l_river(k)
      if(riverdir(k)==2)jriver(k,3)=jriver(k,3)-l_river(k)
      if(riverdir(k)==4)jriver(k,3)=jriver(k,3)+l_river(k)

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iriver(k,3)=iriver(k,3)-par%timax(1)                    !27-04-09
      jriver(k,3)=jriver(k,3)-par%tjmax(1)                    !27-04-09

      rivertrc_inout(k)=0                                              !01-11-10
      if(iriver(k,3)>=0.and.iriver(k,3)<=imax+1.and.             &     !08-06-12
         jriver(k,3)>=0.and.jriver(k,3)<=jmax+1)rivertrc_inout(k)=1
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      read(3,*)riverupwdist(k)        !24-11-10
      read(3,*)riverflux(k,1)
      if(riverdir(k)==3)riverflux(k,1)=-riverflux(k,1) !06-04-21
      if(riverdir(k)==4)riverflux(k,1)=-riverflux(k,1) !06-04-21


!     read(3,*)hriver(k)
      read(3,'(a)')texte90 !11-11-18
! compter le nmbre de vigules dans texte90
      k0=0
      do k1=1,40
       if(texte90(k1:k1)==',')k0=k0+1
      enddo

      if(index(texte90,'canal')==0) then !m°v°m> 07-12-18
! relire la ligne sachant le nombre d'info qu'elle contient avec la valeur de k0
       backspace 3 !11-11-18
       if(k0==0)then
        read(3,*)hriver(k)
       else
        read(3,*)hriver(k),crossresolriver(k,1),alongresolriver(k,1) !11-11-18
       endif
      else                               !m°v°m> 07-12-18
       if(l_river(k)/=0)riverchannelfile(k)=trim(texte90)
      endif                              !m°v°m> 07-12-18

      read(3,*)river_s(k)
      read(3,*)river_tmin(k)                                           !05/05/04
      read(3,*)river_timeref(k)                                        !08-09-16
      read(3,*)river_tmax(k)                                           !05/05/04
      read(3,*)river_obctype(k)                                        !16-03-17
      read(3,*)realriver(k)
      read(3,*)riverinfo(k)
      read(3,'(a)')riverfile(k)
      read(3,*)dateriver(1,k)                                           &
              ,dateriver(2,k)                                           &
              ,dateriver(3,k)                                           &
              ,dateriver(4,k)                                           &
              ,dateriver(5,k)                                           &
              ,dateriver(6,k)

      iriver(k,1)=iriver(k,3)
      iriver(k,2)=iriver(k,3)
      jriver(k,1)=jriver(k,3)
      jriver(k,2)=jriver(k,3)
      riverflux(k,0)=riverflux(k,1)
      riverflux(k,2)=riverflux(k,1)

  115 continue ! fin de boucle sur K
      close(3)
      if(par%rank==0)write(6,*)'ok!'

!Si un fleuve porte X fois le m*me nom son d*bit est divis* par X: !22-02-10
      do k=1,nriver
       texte80(1)=trim(rivername(k))
       texte80(2)=trim(riverfile(k))
       friver(k)=0.
       do k1=1,nriver
        if(trim(texte80(1))==trim(rivername(k1))) then !aaaa>
         friver(k)=friver(k)+1.
         if(trim(texte80(2))/=trim(riverfile(k1))) then !---->
          write(*,*)'-----------------------------------'
          write(*,*)'Warning: notebook_rivers is not correct.'
          write(*,*)'If rivers ',k,' and ',k1,' have the same name'
          write(*,*)'they should have the same data file.'
!         stop ' Stop in set_rivers.F90'
          flag_stop=1 !19-06-17
         endif                                          !---->
        else                                           !aaaa>
         if(trim(texte80(2))==trim(riverfile(k1))) then !....>
          if(realriver(k)==1.and.realriver(k1)==1) then !....> !08-04-13
          write(*,*)'-----------------------------------'
          write(*,*)'Warning: notebook_rivers is not correct.'
          write(*,*)'If rivers ',k,' and ',k1,' have the same data file'
          write(*,*)'they should have the same name.'
!         stop ' Stop in set_rivers.F90'
          flag_stop=1 !19-06-17
         endif                                          !....>
         endif                                          !....>
        endif                                          !aaaa>
       enddo
      if(realriver(k)==0)riverflux(k,1)=riverflux(k,1)/friver(k) !10-05-11
      enddo ! end of k loop

! Si plusieurs fleuves partagent un m*me point d'embouchure alors                  !05-05-10
! on cumule les flux -> creation river_no
! Exemple: 4 rivieres dont riviere no3 et riviere no4 au m*me point alors river_no(1)=1
! river_no(2)=2 river_no(3)=4 river_no(4)=4
      do k1=1,nriver
        river_no(k1)=k1
        do k2=k1+1,nriver
          if(iriver(k1,3)==iriver(k2,3).and.    &
             jriver(k1,3)==jriver(k2,3))then    !>>>>>>>>>>

             river_no(k1)=k2

          endif                                 !>>>>>>>>>>
        enddo
      enddo

!...........
! Definition d'une zone proche de l'embouchure o* on appliquera un
! sch*ma d'advection 100% upwind
      do j=0,jmax+1
      do i=0,imax+1
!      upwindriver_t(i,j)=1. ! ligne commentee le 25-02-17
! La fonction min garantit que si upwindriver_t a ete rabaissE par une contrainte
! additionnelle anterieure, cette derniere soit toujours prise en compte !25-02-17
       upwindriver_t(i,j)=min(1.,upwindriver_t(i,j)) !25-02-17 ! Le "min"
      enddo
      enddo
      do k=1,nriver
! Pas de schema upwind pour riverdir=-1 !01-05-19
      if(riverdir(k)>=0.and.riverdir(k)<=4) then !m°v°m> !01-05-19
       i1=iriver(k,3) ; j1=jriver(k,3)
       i0=riverupwdist(k)                          !24-11-10
       x0=real(i0)/2.

       do i=i1-i0,i1+i0
       do j=j1-i0,j1+i0
        if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !------->
                upwindriver_t(i,j)=                              &
        max(min(upwindriver_t(i,j)*un                            &
               ,sqrt(real(i-i1)**2+real(j-j1)**2)/x0-1.),zero)
        endif                                              !------->
       enddo
       enddo
      endif                                      !m°v°m> !01-05-19
      enddo

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out) !#MPI Section ! synchro processes
#endif
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ; if(k0/=0)stop ' Stop 271 in set_rivers subroutine' !19-06-17


!......................................................................
      endif ! case_==0
!**********************************************************************


!**********************************************************************
! Check the compatibility of the grid points river sources with the land-sea mask
      if(case_==1) then
!......................................................................

      do 100 k=1,nriver
!     l_river(k)=0

      i=iriver(k,3)
      j=jriver(k,3)


      if(rivertrc_inout(k)==1)     then !rrrrrrrrrrrrrrrrrr>

       if(mask_t(i,j,kmax+1).eq.1.and.riverdir(k)>0) then            !25/11/01 !09-04-19
        write(6,'(a,i4,a,a)')'Mask error for river n°',k,' name: ',trim(rivername(k)) !31-10-14
        write(6,'(a,i4,1x,i4)')'at (i,j) given in notebook_rivers=',i+par%timax(1),j+par%tjmax(1)
        write(6,*)'the sea land mask=1 instead of 0'
!       stop ' donc dans set_rivers.f'
        flag_stop=1 !19-06-17
       endif
! resurgences:
      if(riverdir(k)==-1.and.mask_t(i,j,kmax+1)==0) then !ooo>
         write(6,*)'Err on bottom fresh water source location.' &
                  ,' See fort.xx files'
         write(10+par%rank)'Err on bottom fresh water source location.'
         write(10+par%rank)'i,j (glob)',i+par%timax(1),j+par%tjmax(1)
         write(10+par%rank)'mask_t(kmax)=',mask_t(i,j,kmax) &
                          ,', but should be 1'
         flag_stop=1 
       endif                                             !ooo>
! rivieres de surface: !18-02-21
      if(riverdir(k)==0 .and.mask_t(i,j,kmax+1)==0) then !ooo>
         write(6,*)'Err on surface fresh water source location.' &
                  ,' See fort.xx files'
         write(10+par%rank)'Err on surface fresh water source location.'
         write(10+par%rank)'i,j (glob)',i+par%timax(1),j+par%tjmax(1)
         write(10+par%rank)'mask_t(kmax)=',mask_t(i,j,kmax) &
                          ,', but should be 1'
         flag_stop=1 
       endif                                             !ooo>
      endif                             !rrrrrrrrrrrrrrrrrr>

! EXPLICATIONS:
! IRIVER(K,3) et JRIVER(K,3) sont les coordonn*es i & j du point _z de
! l'embouchure du keme fleuve.
! L_RIVER(K) est la longueur du canal du keme fleuve.
! IRIVER(K,1) et JRIVER(K,1) sont les coordonn*es i & j du point _z de
! la source (du canal) du keme fleuve.
! IRIVER(K,2) et JRIVER(K,2) sont les coordonn*es i & j du point de courant
! (x ou y) de la source (du canal) du keme fleuve.

!-------------------------------------------------------------------
!     if(riverdir(k).eq.0) then                                        !25/11/01
! cas des petits fleuves introduits comme flux de surface
!      write(6,*)'filiere fleuve condition de surface non disponible'   !28/03/02
!     flag_stop=1 !19-06-17
!      if(  mask_t(i,j,kmax+1).ne.1) then
!       write(6,'(a,a)')'Mask error for river ',rivername(k)
!       write(6,*)'dans le cas des fleuves introduits comme condition'
!       write(6,*)'de surface, l''embouchure doit etre en mer!'
!     flag_stop=1 !19-06-17
!      endif
!      iriver(k,1)=iriver(k,3)
!      iriver(k,2)=iriver(k,1)
!      jriver(k,1)=jriver(k,3)
!      jriver(k,2)=jriver(k,3)
!     endif

!-------------------------------------------------------------------


      if(.not.allocated(glob_mask)) then
      write(6,*)'glob_mask must be usable but is not available!'
      stop ' STOP in subroutine set_rivers'
      endif

      rivervel_inout(k)=0

!-------------------------------------------------------------------

      if(riverdir(k)==1) then

      iriver(k,1)=iriver(k,3)
      iriver(k,2)=iriver(k,1)+1
      jriver(k,1)=jriver(k,3)
      jriver(k,2)=jriver(k,3)

      if(rivertrc_inout(k)==1)     then !rrrrrrrrrrrrrrrrrr>
       if(i+1<=imax+2) then
! On teste glob_mask et non pas mask_t car le point de grille voisin que l'on veut
! verifier est peut etre dans un sous domaine mpi voisin
!      if(   mask_t(i+1             ,j,kmax+1).ne.1) then
       if(glob_mask(i+1+par%timax(1),j+par%tjmax(1))/=1) then
        write(6,*)'.............................................'
        write(6,'(a,a)')'Mask error A for river ',rivername(k)
        write(6,*)'l''embouchure doit etre devant un point de mer (1)'
        write(6,*)'Local  indexes (i,j)=',i,j
        write(6,*)'global indexes (i,j)=',i+par%timax(1),j+par%tjmax(1)
        write(6,*)'imax jmax=           ',imax,jmax
        write(6,*)'iglb jglb=          =',iglb,jglb
        write(6,*)'par%rank=            ',par%rank
!       write(6,*)'mask_t(i  ,j,kmax+1)=',mask_t(i  ,j,kmax+1)
!       write(6,*)'mask_t(i+1,j,kmax+1)=',mask_t(i+1,j,kmax+1)
!       write(6,*)'mask_t(i  ,j,kmax+1)=',mask_t(i  ,j,kmax+1)
        write(6,*)'glob_mask(i+1,j)',glob_mask(i+1+par%timax(1),j+par%tjmax(1))
!       stop ' STOP in subroutine set_rivers'
        flag_stop=1 !19-06-17
       endif
       endif
      endif                             !rrrrrrrrrrrrrrrrrr>

      if(iriver(k,2)>=1.and.iriver(k,2)<=imax+1.and.               &
         jriver(k,2)>=1.and.jriver(k,2)<=jmax  )rivervel_inout(k)=1  !08-05-14

      endif
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      if(riverdir(k)==2) then

      jriver(k,1)=jriver(k,3)
      jriver(k,2)=jriver(k,1)+1
      iriver(k,1)=iriver(k,3)
      iriver(k,2)=iriver(k,3)

      if(rivertrc_inout(k)==1)     then !rrrrrrrrrrrrrrrrrr>
       if(j+1<=jmax+2) then
! On teste glob_mask et non pas mask_t car le point de grille voisin que l'on veut
! verifier est peut etre dans un sous domaine mpi voisin
!      if(  mask_t(i              ,j+1,kmax+1).ne.1) then
       if(glob_mask(i+par%timax(1),j+1+par%tjmax(1))/=1) then
        write(6,*)'.............................................'
        write(6,'(a,a)')'Mask error B for river ',rivername(k)
        write(6,*)'l''embouchure doit etre devant un point de mer (1)'
        write(6,*)'Local  indexes (i,j)=',i,j
        write(6,*)'global indexes (i,j)=',i+par%timax(1),j+par%tjmax(1)
        write(6,*)'imax jmax=           ',imax,jmax
        write(6,*)'iglb jglb=          =',iglb,jglb
        write(6,*)'par%rank=            ',par%rank
!       write(6,*)'mask_t(i,j  ,kmax+1)=',mask_t(i,j  ,kmax+1)
!       write(6,*)'mask_t(i,j+1,kmax+1)=',mask_t(i,j+1,kmax+1)
        write(6,*)'glob_mask(i,j+1)',glob_mask(i+par%timax(1),j+1+par%tjmax(1))
!       stop ' STOP in subroutine set_rivers'
        flag_stop=1 !19-06-17
       endif
       endif
      endif                             !rrrrrrrrrrrrrrrrrr>

      if(iriver(k,2)>=1.and.iriver(k,2)<=imax  .and.               &
         jriver(k,2)>=1.and.jriver(k,2)<=jmax+1)rivervel_inout(k)=1

      endif
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      if(riverdir(k)==3) then

!     riverflux(k,1)=-riverflux(k,1) !deplacE en amont le 06-04-21 
      iriver(k,1)=iriver(k,3)
      iriver(k,2)=iriver(k,1)
      jriver(k,1)=jriver(k,3)
      jriver(k,2)=jriver(k,3)

      if(rivertrc_inout(k)==1)     then !rrrrrrrrrrrrrrrrrr>
       if(i-1>=-1) then !21-01-15
! On teste glob_mask et non pas mask_t car le point de grille voisin que l'on veut
! verifier est peut etre dans un sous domaine mpi voisin
!      if(   mask_t(i-1             ,j,kmax+1).ne.1) then
       if(glob_mask(i-1+par%timax(1),j+par%tjmax(1))/=1) then
        write(6,*)'.............................................'
        write(6,'(a,a)')'Mask error C for river ',rivername(k)
        write(6,*)'l''embouchure doit etre devant un point de mer (1)'
        write(6,*)'Local  indexes (i,j)=',i,j
        write(6,*)'global indexes (i,j)=',i+par%timax(1),j+par%tjmax(1)
        write(6,*)'imax jmax=           ',imax,jmax
        write(6,*)'iglb jglb=          =',iglb,jglb
        write(6,*)'par%rank=            ',par%rank
        write(6,*)'glob_mask(i-1,j)',glob_mask(i-1+par%timax(1),j+par%tjmax(1))
!       write(6,*)'mask_t(i  ,j,kmax+1)=',mask_t(i  ,j,kmax+1)
!       write(6,*)'mask_t(i-1,j,kmax+1)=',mask_t(i-1,j,kmax+1)
!       stop ' STOP in subroutine set_rivers'
        flag_stop=1 !19-06-17
       endif
       endif
      endif                             !rrrrrrrrrrrrrrrrrr>

      if(iriver(k,2)>=1.and.iriver(k,2)<=imax+1.and.               &
         jriver(k,2)>=1.and.jriver(k,2)<=jmax  )rivervel_inout(k)=1

      endif
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      if(riverdir(k)==4) then

!     riverflux(k,1)=-riverflux(k,1) ! deplacE en amont le 06-04-21
      jriver(k,1)=jriver(k,3)
      jriver(k,2)=jriver(k,1)
      iriver(k,1)=iriver(k,3)
      iriver(k,2)=iriver(k,3)

      if(rivertrc_inout(k)==1)     then !rrrrrrrrrrrrrrrrrr>
       if(j-1>=-1) then !21-01-15
! On teste glob_mask et non pas mask_t car le point de grille voisin que l'on veut
! verifier est peut etre dans un sous domaine mpi voisin
!      if(   mask_t(i             ,j-1,kmax+1).ne.1) then
       if(glob_mask(i+par%timax(1),j-1+par%tjmax(1))/=1) then
        write(6,*)'.............................................'
        write(6,'(a,a)')'Mask error D for river ',rivername(k)
        write(6,*)'l''embouchure doit etre devant un point de mer (1)'
        write(6,*)'Local  indexes (i,j)=',i,j
        write(6,*)'global indexes (i,j)=',i+par%timax(1),j+par%tjmax(1)
        write(6,*)'imax jmax=           ',imax,jmax
        write(6,*)'iglb jglb=          =',iglb,jglb
        write(6,*)'par%rank=            ',par%rank
        write(6,*)'glob_mask(i,j-1)',glob_mask(i+par%timax(1),j-1+par%tjmax(1)) !14-04-19
!       write(6,*)'mask_t(i,j  ,kmax+1)=',mask_t(i,j  ,kmax+1)
!       write(6,*)'mask_t(i,j-1,kmax+1)=',mask_t(i,j-1,kmax+1)
!       stop ' STOP in subroutine set_rivers'
        flag_stop=1 !19-06-17
       endif
       endif
      endif                             !rrrrrrrrrrrrrrrrrr>

      if(iriver(k,2)>=1.and.iriver(k,2)<=imax  .and.               &
         jriver(k,2)>=1.and.jriver(k,2)<=jmax+1)rivervel_inout(k)=1

      endif
!-------------------------------------------------------------------

!     if(iriver(k,2)>=0.and.iriver(k,2)<=imax+1.and.               &
!        jriver(k,2)>=0.and.jriver(k,2)<=jmax+1)rivervel_inout(k)=1


  100 continue ! fin de boucle sur K

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out) !#MPI Section ! synchro processes
#endif
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ; if(k0/=0)stop ' Stop 493 in set_rivers subroutine' !19-06-17


!......................................................................
      endif    ! case_=1
!**********************************************************************



!**********************************************************************
      if(case_==2) then

! Creuser des canaux pour les fleuves
! Attention A ce stade les sources ont Ete deplacees au point amont des canaux
! A construire:
! Appelee apres les lissages de H, cette routine donne la valeur
! definitive de H dans les canaux, ce qui suppose de traiter les
! tableaux mask_t et h_w et non pas s'arreter A glob_mask et glob_h

! ATTENTION MEME SI L_RIVER=0 L'algo suivant agit au moins comme une condition qui
! assure que la bathy du point traceur source est identique au point oceanique adjacent (et pas n'importe quoi...)


      do k=1,nriver !RRRRRRRR> !07-03-17
       if(l_river(k)/=0) then !lllll> !02-05-2018 ce commentaire permet d'appliquer une condition sur H meme si l_river=0

       if(riverchannelfile(k)/='none') then !moLom> !07-12-18

! CAS PARTICULIER oU les parametres du canal sont lus dans un fichier
       unit_=s_unit(7)
       open(unit=unit_,file=trim(riverchannelfile(k))) !06-12-18
  553   read(unit_,*,end=554)i1,j1,k1,x1,x2,x3,x4,x5 !i,j,mask,h,dx,dy,lon,lat 
        i3=i1-par%timax(1) ; j3=j1-par%tjmax(1) ; glob_h(i1,j1)=x1 ; glob_mask(i1,j1)=k1
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then !>>>
           mask_t(i3,j3,:)=glob_mask(i1,j1) ; h_w(i3,j3)=glob_h(i1,j1)
          endif                                                  !>>>
        goto 553
  554  close(unit_)

       else                                 !moLom> !07-12-18

! CAS STANDARD utilisant les parametres du notebook_river: DEBUT
       write(texte30,'(a,i0,a,i0,a,a)')'tmp/canal',k,'_',par%rank,'_',trim(rivername(k)) !06-12-18
       unit_=s_unit(7)
  

!       x10=(alongresolriver(k,2)/alongresolriver(k,1))**(1./l_river(k)) ! alongriver loi de proportionalite
!       dx_t(i3,j3)=alongresolriver(k,1)*x10**(i-i1)     ! alongriver loi de proportionalite !11-11-18

        i1=iriver(k,3)+par%timax(1) ; j1=jriver(k,3)+par%tjmax(1) !Coord point source en amont du canal   

        if(riverdir(k)==1) then !-cas i croissant->
         j=j1 ; i2=i1+l_river(k) !Coord point embouchure
         do i=i1+1,i2
! Note: pas d'interpolation lineraire de i si grille irreguliere (si 
! sinon la pente n'est pas constante dans l'espace physique
! Plus loin les poids sont donnEs par le rapport des mailles afin 
! de garantir une pente constante dans l'espace physique
          glob_mask(i,j)=1
          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
           if(alongresolriver(k,1)-alongresolriver(k,2)/=0.) then !ooo>
            rap=(dx_t(i3,j3)-alongresolriver(k,2))/(alongresolriver(k,1)-alongresolriver(k,2))
           else                                                   !ooo>
            rap=real(i2-i)/real(l_river(k))
           endif                                                  !ooo>
           mask_t(i3,j3,:)=1 
              h_w(i3,j3)=(1.-rap)*glob_h(i2+1,j1)+rap*hriver(k)
           open(unit=unit_,file=trim(texte30),position='append') !06-12-18
!          write(unit_,'(2i6,5(1x,e17.7))')i,j,h_w(i3,j3),dy_t(i3,j3),dx_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg !06-12-18
!          write(unit_,'(2i6,5(1x,e17.7),a)')i,j,h_w(i3,j3),dx_t(i3,j3),dy_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg,'   ! i,j,h_w,dx_t,dy_t,lon_t,lat_t'' !06-12-18
           write(unit_,'(3i6,3(1x,f10.4),2(1x,f11.6),a)')i,j,mask_t(i3,j3,kmax),h_w(i3,j3),dx_t(i3,j3),dy_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg,'   ! i,j,mask,h,dx,dy,lon,lat' !07-12-18
           close(unit_)
          endif

         enddo
        endif                   !-cas i croissant->

        if(riverdir(k)==3) then !-cas i decroissant->
         j=j1 ; i2=i1-l_river(k)
         do i=i2,i1-1
          glob_mask(i,j)=1
          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
           if(alongresolriver(k,1)-alongresolriver(k,2)/=0.) then !ooo>
            rap=(dx_t(i3,j3)-alongresolriver(k,2))/(alongresolriver(k,1)-alongresolriver(k,2))
           else                                                   !ooo>
            rap=real(i-i2)/real(l_river(k))
           endif                                                  !ooo>
           mask_t(i3,j3,:)=1
              h_w(i3,j3)=(1.-rap)*glob_h(i2-1,j1)+rap*hriver(k)
           open(unit=unit_,file=trim(texte30),position='append') !06-12-18
           write(unit_,'(3i6,3(1x,f10.4),2(1x,f11.6),a)')i,j,mask_t(i3,j3,kmax),h_w(i3,j3),dx_t(i3,j3),dy_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg,'   ! i,j,mask,h,dx,dy,lon,lat' !07-12-18
           close(unit_)
          endif

         enddo
        endif                   !-cas i decroissant->

        if(riverdir(k)==2) then !-cas j croissant->
         i=i1 ; j2=j1+l_river(k) !Coord point embouchure
         do j=j1+1,j2
          glob_mask(i,j)=1
          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
           if(alongresolriver(k,1)-alongresolriver(k,2)/=0.) then !ooo>
            rap=(dy_t(i3,j3)-alongresolriver(k,2))/(alongresolriver(k,1)-alongresolriver(k,2))
           else                                                   !ooo>
            rap=real(j2-j)/real(l_river(k))
           endif                                                  !ooo>
           mask_t(i3,j3,:)=1
              h_w(i3,j3)=(1.-rap)*glob_h(i1,j2+1)+rap*hriver(k)
           open(unit=unit_,file=trim(texte30),position='append') !06-12-18
           write(unit_,'(3i6,3(1x,f10.4),2(1x,f11.6),a)')i,j,mask_t(i3,j3,kmax),h_w(i3,j3),dx_t(i3,j3),dy_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg,'   ! i,j,mask,h,dx,dy,lon,lat' !07-12-18
           close(unit_)
          endif

         enddo
        endif                   !-cas j croissant->

        if(riverdir(k)==4) then !-cas j decroissant->
         i=i1 ; j2=j1-l_river(k) !Coord point embouchure
         do j=j2,j1-1
          glob_mask(i,j)=1
          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
           if(alongresolriver(k,1)-alongresolriver(k,2)/=0.) then !ooo>
            rap=(dy_t(i3,j3)-alongresolriver(k,2))/(alongresolriver(k,1)-alongresolriver(k,2))
           else                                                   !ooo>
            rap=real(j-j2)/real(l_river(k))
           endif                                                  !ooo>
           mask_t(i3,j3,:)=1
              h_w(i3,j3)=(1.-rap)*glob_h(i1,j2-1)+rap*hriver(k)
           open(unit=unit_,file=trim(texte30),position='append') !06-12-18
           write(unit_,'(3i6,3(1x,f10.4),2(1x,f11.6),a)')i,j,mask_t(i3,j3,kmax),h_w(i3,j3),dx_t(i3,j3),dy_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg,'   ! i,j,mask,h,dx,dy,lon,lat' !07-12-18
           close(unit_)
          endif

         enddo
        endif                   !-cas j decroissant->

! Enfin le point i1,j1:
        glob_h(i1,j1)=hriver(k)
        i3=i1-par%timax(1) ; j3=j1-par%tjmax(1) 
        if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then !m°v°m> !07-12-18
              h_w(i3,j3)=glob_h(i1,j1) !01-05-17
        write(texte30,'(a,i0,a,i0,a,a)')'tmp/canal',k,'_',par%rank,'_s_',trim(rivername(k)) !06-12-18
        open(unit=unit_,file=trim(texte30),position='append') !06-12-18
        write(unit_,'(3i6,3(1x,f10.4),2(1x,f11.6),a)')i1,j1,mask_t(i3,j3,kmax),h_w(i3,j3),dx_t(i3,j3),dy_t(i3,j3),lon_t(i3,j3)*rad2deg,lat_t(i3,j3)*rad2deg,'   ! i,j,mask,h,dx,dy,lon,lat' !07-12-18
        close(unit_)
        endif                                                  !m°v°m> !07-12-18

! CAS STANDARD utilisant les parametres du notebook_river: FIN
       endif                                !moLom> !07-12-18

       else                   !lllll>

! Si pas de canal additionnel, la bathy du point source = bathy du point oceanique adjacent !02-05-18
        i=iriver(k,1)     ; j=jriver(k,1) 
        if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)then !m°v°w> !24-05-18
         i3=i+par%timax(1) ; j3=j+par%tjmax(1)
         if(riverdir(k)==1)i3=i3+1
         if(riverdir(k)==3)i3=i3-1
         if(riverdir(k)==2)j3=j3+1
         if(riverdir(k)==4)j3=j3-1
         i3=min(max(i3,1),iglb)
         j3=min(max(j3,1),jglb)
         h_w(i,j)=glob_h(i3,j3)
!        write(10+par%rank,'(a,a)')'River ',trim(rivername(k))
!        write(10+par%rank,*)'No ',k
!        write(10+par%rank,*)'i,j glob',i+par%timax(1),j+par%tjmax(1)
!        write(10+par%rank,*)'i3,j3,H',i3,j3,glob_h(i3,j3)
        endif                                             !m°v°w> !24-05-18
       

       endif                  !lllll>
      enddo         !RRRRRRRR>


      endif
!**********************************************************************


!***************************************************************************
! INITIALISATION DE LA TEMPERATURE ET DE LA SALINITE AUX POINTS SOURCES
! ET DANS LES CANAUX: OPTION REPORTER DANS LE CANAL LE PROFIL DEVANT
! L'EMBOUCHURE
! DEBUT
      if(case_==3) then
!***************************************************************************

      return

      do kr=1,nriver
       if(l_river(kr)/=0) then
        write(6,*)' Erreur parametrage notebookrivers!' &
       ,' Ne pas demander de creuser de canal pour la riviere ' &
       ,trim(rivername(kr))
!       stop ' STOP in subroutine set_rivers'
        flag_stop=1 !19-06-17
       endif
      enddo
      return !08-05-14

      ksecu=0
      do kr=1,nriver ! debut de boucle sur KR

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 29-05-09
! Si le sous-domaine correspond alors faire les tests de debug:
!     if(river_dom(kr) == par%rank) then !rrrrrrrrrrrrrrrrrr>
      if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>          !01-11-10
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 29-05-09

! S*curit* d*bugueur:                                                  !27/07/04
      if(river_s(kr).lt.-9999.)ksecu=1
      if(river_t(kr,1).lt.-9999.)ksecu=1
      if(ksecu.eq.1)then ! debug debug >

      write(6,*)'tentative d''initialisation de t & s dans les canaux'
      write(6,*)'des fleuves alors que les tableaux river_s et'
      write(6,*)'river_t n''ont pas encore *t* initialis*s'

!     stop ' dans set_rivers.f case_.eq.3'
      flag_stop=1 !19-06-17
      endif              ! debug debug >

      if(riverdir(kr)>0) then !0000000000000> !08-04-13 !09-04-19

      if(riverdir(kr).eq.1) then !11111111111111>

      i=iriver(kr,3)+1
      j=jriver(kr,3)

      j2=j
      j3=j
      j4=1

      i2=i
!     i3=i-(l_river(kr)+1)
      i3=i-l_river(kr)    !02-04-14
      i4=-1

      endif                      !11111111111111>

      if(riverdir(kr).eq.2) then !22222222222222>
      i=iriver(kr,3)
      j=jriver(kr,3)+1

      i2=i
      i3=i
      i4=1

      j2=j
!     j3=j-(l_river(kr)+1)
      j3=j-l_river(kr)
      j4=-1

      endif                      !22222222222222>

      if(riverdir(kr).eq.3) then !33333333333333>
      i=iriver(kr,3)-1
      j=jriver(kr,3)

      j2=j
      j3=j
      j4=1

      i2=i
!     i3=i+l_river(kr)+1
      i3=i+l_river(kr)
      i4=1

      endif                      !3333333333333>

      if(riverdir(kr).eq.4) then !4444444444444>
      i=iriver(kr,3)
      j=jriver(kr,3)-1

      i2=i
      i3=i
      i4=1

      j2=j
!     j3=j+l_river(kr)+1
      j3=j+l_river(kr)
      j4=1

      endif                      !4444444444444>

      do j1=j2,j3,j4 ! debut de boucle sur J1
      do i1=i2,i3,i4 ! debut de boucle sur I1

         ssh_int_w(i1,j1,0)=ssh_int_w(i,j,0)
         ssh_int_w(i1,j1,1)=ssh_int_w(i,j,1)
         ssh_int_w(i1,j1,2)=ssh_int_w(i,j,2)
             ssh_w(i1,j1,:)=ssh_int_w(i,j,2)



         do k1=kmin_w(i1,j1),kmax ! debut de boucle sur K1

          do k=kmin_w(i ,j ),kmax-1 ! debut de boucle sur K

           if( depth_t(i1,j1,k1).ge. depth_t(i,j,k).and.                &
               depth_t(i1,j1,k1).le. depth_t(i,j,k+1))then !>>>>>>

              rap=( depth_t(i1,j1,k1 )- depth_t(i,j,k))                 &
                 /( depth_t(i ,j ,k+1)- depth_t(i,j,k))

              tem_t(i1,j1,k1,1)=(1.-rap)*tem_t(i,j,k  ,1)               &
                                   +rap *tem_t(i,j,k+1,1)
              sal_t(i1,j1,k1,1)=(1.-rap)*sal_t(i,j,k  ,1)               &
                                   +rap *sal_t(i,j,k+1,1)

           endif                                          !>>>>>>

          enddo                   ! fin de boucke sur K

          if( depth_t(i1,j1,k1).le. depth_t(i,j,kmin_w(i,j)))then !------->
              tem_t(i1,j1,k1,1)=tem_t(i,j,kmin_w(i,j),1)
              sal_t(i1,j1,k1,1)=sal_t(i,j,kmin_w(i,j),1)
          endif                                                   !------->

          if( depth_t(i1,j1,k1).ge. depth_t(i,j,kmax))       then !------->
              tem_t(i1,j1,k1,1)=tem_t(i,j,kmax,1)
              sal_t(i1,j1,k1,1)=sal_t(i,j,kmax,1)
          endif                                                   !------->

          tem_t(i1,j1,k1,0)=tem_t(i1,j1,k1,1)                           !11-10-09
          tem_t(i1,j1,k1,2)=tem_t(i1,j1,k1,1)                           !11-10-09
          sal_t(i1,j1,k1,0)=sal_t(i1,j1,k1,1)                           !11-10-09
          sal_t(i1,j1,k1,2)=sal_t(i1,j1,k1,1)                           !11-10-09

         enddo                    ! fin de boucle sur K1

      enddo ! fin de boucle sur I1
      enddo ! fin de boucle sur J1

      endif                    !0000000000000>


!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Si le sous-domaine correspond alors faire les tests de debug:
      endif                          !rrrrrrrrrrrrrrrrrr>
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo ! fin de boucle sur KR
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ; if(k0/=0)stop ' Stop 756 in set_rivers subroutine' !19-06-17


!***************************************************************************
! INITIALISATION DE LA TEMPERATURE ET DE LA SALINITE AUX POINTS SOURCES
! ET DANS LES CANAUX: OPTION REPORTER DANS LE CANAL LE PROFIL DEVANT
! L'EMBOUCHURE
! FIN
      return
      endif
!***************************************************************************

!***************************************************************************
! INITIALISATION DE LA TEMPERATURE ET DE LA SALINITE AUX POINTS SOURCES
! ET DANS LES CANAUX                                                   !21/08/03
! DEBUT
      if(case_==4) then
!***************************************************************************
      ksecu=0
      do 201 kr=1,nriver                                               !21/08/03

! S*curit* d*bugueur:                                                  !27/07/04
      if(river_s(kr).lt.-9999.)ksecu=1
      if(river_t(kr,1).lt.-9999.)ksecu=1
      if(ksecu.eq.1)then ! debug debug >

      write(6,*)'tentative d''initialisation de t & s dans les canaux'
      write(6,*)'des fleuves alors que les tableaux river_s et'
      write(6,*)'river_t n''ont pas encore *t* initialis*s'

!     stop ' dans set_rivers.f case_.eq.3'
      flag_stop=1 !19-06-17
      endif              ! debug debug >

!-------------------------------------------------------------------
      if(riverdir(kr).eq.1) then
      i=iriver(kr,3)+1
      j=jriver(kr,3)

! Calcul analogue * celui de la construction du canal:
         do i1=i,i-(l_river(kr)+1),-1
         x1=real(i-i1)/amax1(real(l_river(kr)+1),1.e-10)
         rap=(exp(-5.*x1)-1.)/(exp(-5.)-1.)

         ssh_int_w(i1,j,0)=ssh_int_w(i,j,0)
         ssh_int_w(i1,j,1)=ssh_int_w(i,j,1)
         ssh_int_w(i1,j,2)=ssh_int_w(i,j,2)
             ssh_w(i1,j,:)=ssh_int_w(i,j,2)

          x2=rap                                                       !20/02/04
          do k=1,kmax !/////////////////////>
          rap=x2*exp( depth_t(i1,j ,k)/2.)                             !20/02/04

           tem_t(i1,j,k,1)=rap*river_t(kr,1)+(1.-rap)*tem_t(i,j,k,1)   !05/05/04
           sal_t(i1,j,k,1)=rap*river_s(kr)  +(1.-rap)*sal_t(i,j,k,1)

!           do k1=0,2 !----------------------------->
!            thz_c(i1,j,k,k1)=tem_c(i1,j,k,1)
!    &                 *     dsig_c(i1,j,k)
!    &                 *(       h_z(i1,j)
!    &                   +zta_int_z(i1,j,k1))
!            shz_c(i1,j,k,k1)=sal_c(i1,j,k,1)
!    &                 *     dsig_c(i1,j,k)
!    &                 *(       h_z(i1,j)
!    &                   +zta_int_z(i1,j,k1))
!           enddo     !----------------------------->
            tem_t(i1,j,k,0)=tem_t(i1,j,k,1)                            !11-10-09
            tem_t(i1,j,k,2)=tem_t(i1,j,k,1)                            !11-10-09
            sal_t(i1,j,k,0)=sal_t(i1,j,k,1)                            !11-10-09
            sal_t(i1,j,k,2)=sal_t(i1,j,k,1)                            !11-10-09



          enddo       !/////////////////////>

         enddo

      endif
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      if(riverdir(kr).eq.2) then
      i=iriver(kr,3)
      j=jriver(kr,3)+1

! Calcul analogue * celui de la construction du canal:
         do j1=j,j-(l_river(kr)+1),-1
         x1=real(j-j1)/amax1(real(l_river(kr)+1),1.e-10)
         rap=(exp(-5.*x1)-1.)/(exp(-5.)-1.)                            !03/09/01

         ssh_int_w(i,j1,0)=ssh_int_w(i,j,0)
         ssh_int_w(i,j1,1)=ssh_int_w(i,j,1)
         ssh_int_w(i,j1,2)=ssh_int_w(i,j,2)
             ssh_w(i,j1,:)=ssh_int_w(i,j,2)

          x2=rap                                                       !20/02/04
          do k=1,kmax !/////////////////////>
          rap=x2*exp( depth_t(i ,j1,k)/2.)                             !20/02/04

           tem_t(i,j1,k,1)=rap*river_t(kr,1)+(1.-rap)*tem_t(i,j,k,1)    !05/05/04
           sal_t(i,j1,k,1)=rap*river_s(kr  )+(1.-rap)*sal_t(i,j,k,1)

!           do k1=0,2 !----------------------------->
!            thz_c(i,j1,k,k1)=tem_c(i,j1,k,1)
!    &                 *     dsig_c(i,j1,k)
!    &                 *(       h_z(i,j1)
!    &                   +zta_int_z(i,j1,k1))
!            shz_c(i,j1,k,k1)=sal_c(i,j1,k,1)
!    &                 *     dsig_c(i,j1,k)
!    &                 *(       h_z(i,j1)
!    &                   +zta_int_z(i,j1,k1))
!           enddo     !----------------------------->
            tem_t(i,j1,k,0)=tem_t(i,j1,k,1)                            !11-10-09
            tem_t(i,j1,k,2)=tem_t(i,j1,k,1)                            !11-10-09
            sal_t(i,j1,k,0)=sal_t(i,j1,k,1)                            !11-10-09
            sal_t(i,j1,k,2)=sal_t(i,j1,k,1)                            !11-10-09

          enddo       !/////////////////////>


         enddo

      endif
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      if(riverdir(kr).eq.3) then
      i=iriver(kr,3)-1
      j=jriver(kr,3)

! Calcul analogue * celui de la construction du canal:
         do i1=i,i+l_river(kr)+1
         x1=real(i1-i)/amax1(real(l_river(kr)+1),1.e-10)
         rap=(exp(-5.*x1)-1.)/(exp(-5.)-1.)                            !03/09/01

         ssh_int_w(i1,j,0)=ssh_int_w(i,j,0)
         ssh_int_w(i1,j,1)=ssh_int_w(i,j,1)
         ssh_int_w(i1,j,2)=ssh_int_w(i,j,2)
             ssh_w(i1,j,:)=ssh_int_w(i,j,2)

          x2=rap                                                        !20/02/04
          do k=1,kmax !/////////////////////>
          rap=x2*exp( depth_t(i1,j ,k)/2.)                              !20/02/04

           tem_t(i1,j,k,1)=rap*river_t(kr,1)+(1.-rap)*tem_t(i,j,k,1)    !05/05/04
           sal_t(i1,j,k,1)=rap*river_s(kr)  +(1.-rap)*sal_t(i,j,k,1)

!           do k1=0,2 !----------------------------->
!            thz_c(i1,j,k,k1)=tem_c(i1,j,k,1)
!    &                 *     dsig_c(i1,j,k)
!    &                 *(       h_z(i1,j)
!    &                   +zta_int_z(i1,j,k1))
!            shz_c(i1,j,k,k1)=sal_c(i1,j,k,1)
!    &                 *     dsig_c(i1,j,k)
!    &                 *(       h_z(i1,j)
!    &                   +zta_int_z(i1,j,k1))
!           enddo     !----------------------------->
            tem_t(i1,j,k,0)=tem_t(i1,j,k,1)                            !11-10-09
            tem_t(i1,j,k,2)=tem_t(i1,j,k,1)                            !11-10-09
            sal_t(i1,j,k,0)=sal_t(i1,j,k,1)                            !11-10-09
            sal_t(i1,j,k,2)=sal_t(i1,j,k,1)                            !11-10-09
          enddo       !/////////////////////>

         enddo

      endif
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      if(riverdir(kr).eq.4) then
      i=iriver(kr,3)
      j=jriver(kr,3)-1

! Calcul analogue * celui de la construction du canal:
         do j1=j,j+l_river(kr)+1
         x1=real(j1-j)/amax1(real(l_river(kr)+1),1.e-10)
         rap=(exp(-5.*x1)-1.)/(exp(-5.)-1.)                            !03/09/01

         ssh_int_w(i,j1,0)=ssh_int_w(i,j,0)
         ssh_int_w(i,j1,1)=ssh_int_w(i,j,1)
         ssh_int_w(i,j1,2)=ssh_int_w(i,j,2)
            ssh_w(i,j1,0)=ssh_int_w(i,j,0)
            ssh_w(i,j1,1)=ssh_int_w(i,j,1)
            ssh_w(i,j1,2)=ssh_int_w(i,j,2)

          x2=rap                                                        !20/02/04
          do k=1,kmax !/////////////////////>
          rap=x2*exp( depth_t(i ,j1,k)/2.)                              !20/02/04

           tem_t(i,j1,k,1)=rap*river_t(kr,1)+(1.-rap)*tem_t(i,j,k,1)    !05/05/04
           sal_t(i,j1,k,1)=rap*river_s(kr)  +(1.-rap)*sal_t(i,j,k,1)

!           do k1=0,2 !----------------------------->
!            thz_c(i,j1,k,k1)=tem_c(i,j1,k,1)
!    &                 *     dsig_c(i,j1,k)
!    &                 *(       h_z(i,j1)
!    &                   +zta_int_z(i,j1,k1))
!            shz_c(i,j1,k,k1)=sal_c(i,j1,k,1)
!    &                 *     dsig_c(i,j1,k)
!    &                 *(       h_z(i,j1)
!    &                   +zta_int_z(i,j1,k1))
!           enddo     !----------------------------->
            tem_t(i,j1,k,0)=tem_t(i,j1,k,1)                            !11-10-09
            tem_t(i,j1,k,2)=tem_t(i,j1,k,1)                            !11-10-09
            sal_t(i,j1,k,0)=sal_t(i,j1,k,1)                            !11-10-09
            sal_t(i,j1,k,2)=sal_t(i,j1,k,1)                            !11-10-09
          enddo       !/////////////////////>

         enddo

      endif
!-------------------------------------------------------------------

  201 continue
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ; if(k0/=0)stop ' Stop 970 in set_rivers subroutine' !19-06-17


!***************************************************************************
! INITIALISATION DE LA TEMPERATURE ET DE LA SALINITE AUX POINTS SOURCES
! ET DANS LES CANAUX                                                   !21/08/03
! FIN
      return
      endif
!***************************************************************************



! rdv au bord du fleuve:
      endif

      end subroutine set_rivers

!................................................................................

      subroutine set_rivers_reservoir_sal(t_) !14-07-17
      use module_principal ; use module_parallele
      implicit none
      integer t_

! Salinity within rivers reservoirs

      do k=1,nriver !RRRRRRRR> !07-03-17

!      write(10+par%rank,*)k,l_river(k)

       if(l_river(k)/=0) then !lllll>
        i1=iriver(k,3)+par%timax(1) ; j1=jriver(k,3)+par%tjmax(1) !Coord point source en amont du canal   

        if(riverdir(k)==1) then !-cas i croissant->
         j=j1 ; i2=i1+l_river(k) !Coord point embouchure
         do i=i1+1,i2

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             salobc_t(i3,j3,:,t_)=river_s(k)
          endif

         enddo
        endif                   !-cas i croissant->

        if(riverdir(k)==3) then !-cas i decroissant->
         j=j1 ; i2=i1-l_river(k)
         do i=i2,i1-1

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             salobc_t(i3,j3,:,t_)=river_s(k)
          endif

         enddo
        endif                   !-cas i decroissant->

        if(riverdir(k)==2) then !-cas j croissant->
         i=i1 ; j2=j1+l_river(k) !Coord point embouchure
         do j=j1+1,j2

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             salobc_t(i3,j3,:,t_)=river_s(k)
          endif

         enddo
        endif                   !-cas j croissant->

        if(riverdir(k)==4) then !-cas j decroissant->
         i=i1 ; j2=j1-l_river(k) !Coord point embouchure
         do j=j2,j1-1

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             salobc_t(i3,j3,:,t_)=river_s(k)
          endif

         enddo
        endif                   !-cas j decroissant->

       endif                  !lllll>
      enddo         !RRRRRRRR>


      end subroutine set_rivers_reservoir_sal

!................................................................................

      subroutine set_rivers_reservoir_z0 !14-07-17
      use module_principal ; use module_parallele
      implicit none

! z0 bottom lenght in rivers reservoir
      x0=z0b_rivers  ! z0 default value given in notebook_visco !13-12-19
!     x0=0.01 ! User's choice for z0 (in meters)

      do k=1,nriver !RRRRRRRR> !07-03-17

!      write(10+par%rank,*)k,l_river(k)

       if(l_river(k)/=0) then !lllll>
        i1=iriver(k,3)+par%timax(1) ; j1=jriver(k,3)+par%tjmax(1) !Coord point source en amont du canal   

        if(riverdir(k)==1) then !-cas i croissant->
         j=j1 ; i2=i1+l_river(k) !Coord point embouchure
         do i=i1+1,i2

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             z0_w(i3,j3)=x0
          endif

         enddo
        endif                   !-cas i croissant->

        if(riverdir(k)==3) then !-cas i decroissant->
         j=j1 ; i2=i1-l_river(k)
         do i=i2,i1-1

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             z0_w(i3,j3)=x0
          endif

         enddo
        endif                   !-cas i decroissant->

        if(riverdir(k)==2) then !-cas j croissant->
         i=i1 ; j2=j1+l_river(k) !Coord point embouchure
         do j=j1+1,j2

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             z0_w(i3,j3)=x0
          endif

         enddo
        endif                   !-cas j croissant->

        if(riverdir(k)==4) then !-cas j decroissant->
         i=i1 ; j2=j1-l_river(k) !Coord point embouchure
         do j=j2,j1-1

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
             z0_w(i3,j3)=x0
          endif

         enddo
        endif                   !-cas j decroissant->

       endif                  !lllll>
      enddo         !RRRRRRRR>

! z0b des canaux du nouveau reseau de canaux: !13-12-19
! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_3
      do k=1,nbcanal ! (numero du canal)

! Les canaux eux-meme (c.a.d. gridtype1):
       do i1=i_canalcoord(k,point1,gridtype1,receiver),i_canalcoord(k,point2,gridtype1,receiver)  
       do j1=j_canalcoord(k,point1,gridtype1,receiver),j_canalcoord(k,point2,gridtype1,receiver)  
         i=i1-par%timax(1)
         j=j1-par%tjmax(1)
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)z0_w(i,j)=z0b_rivers !13-12-19
       enddo
       enddo

! Les points connectEs (c.a.d. gridtype2):
! Point 1:
         if(i_canalcoord(k,point1,gridtype2,sender)>0.and. &
            j_canalcoord(k,point1,gridtype2,sender)>0) then !--Si point1 connecte ->

          i1=i_canalcoord(k,point1,gridtype2,sender)
          j1=j_canalcoord(k,point1,gridtype2,sender)
          i=i1-par%timax(1)
          j=j1-par%tjmax(1)
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)z0_w(i,j)=z0b_rivers !13-12-19

          i1=i_canalcoord(k,point1,gridtype2,receiver)
          j1=j_canalcoord(k,point1,gridtype2,receiver)
          i=i1-par%timax(1)
          j=j1-par%tjmax(1)
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)z0_w(i,j)=z0b_rivers !13-12-19

         endif                                              !--Si point1 connecte ->

! Point 2:
         if(i_canalcoord(k,point2,gridtype2,sender)>0.and. &
            j_canalcoord(k,point2,gridtype2,sender)>0) then !--Si point2 connecte ->

         i1=i_canalcoord(k,point2,gridtype2,sender)
         j1=j_canalcoord(k,point2,gridtype2,sender)
         i=i1-par%timax(1)
         j=j1-par%tjmax(1)
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)z0_w(i,j)=z0b_rivers !13-12-19

         i1=i_canalcoord(k,point2,gridtype2,receiver)
         j1=j_canalcoord(k,point2,gridtype2,receiver)
         i=i1-par%timax(1)
         j=j1-par%tjmax(1)
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)z0_w(i,j)=z0b_rivers !13-12-19

         endif                                              !--Si point2 connecte ->

      enddo

      end subroutine set_rivers_reservoir_z0

!................................................................................

      subroutine set_rivers_reservoir_grid
      use module_principal ; use module_parallele ; use module_s
      implicit none
      integer unit_


!     RETURN ! If return dx, dy keep their default values

! dx_t & dy_t in rivers reservoirs:

      do k=1,nriver !RRRRRRRR> !07-03-17

       if((crossresolriver(k,1)>0.and.alongresolriver(k,1)>0) &
                              .or.riverchannelfile(k)/='none') then !m°v°m> !11-11-18

       if(par%rank==0)open(unit=3,file='tmp/messages',position='append')

!     alongresolriver(k,1)=1000.  ! user's choice for along river resolution
!     crossresolriver(k,1)=100.   ! user's choice for cross river resolution

       if(l_river(k)/=0) then !lllll>

       if(riverchannelfile(k)/='none') then !moLom> !07-12-18

! CAS PARTICULIER oU les parametres du canal sont lus dans un fichier
       unit_=s_unit(7)
! Premiere lecture pour charger dx,dy,lon,lat
       open(unit=unit_,file=trim(riverchannelfile(k))) !06-12-18
 1272   read(unit_,*,end=1274)i1,j1,k1,x1,x2,x3,x4,x5 !i,j,mask,h,dx,dy,lon,lat 
        i3=i1-par%timax(1) ; j3=j1-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then !>>>
           dx_t(i3,j3)=x2 ; dy_t(i3,j3)=x3 ! lon_t(i3,j3)=x4*deg2rad ; lat_t(i3,j3)=x5*deg2rad
          endif                                                  !>>>
        goto 1272
 1274  close(unit_)
! Ces lignes sont commentees (voir aussi lon_t(i3,j3)=x4*deg2rad ; lat_t(i3,j3)=x5*deg2rad au dessus)
! en attendant d'etre testees et qu'on ait besoin de les utiliser
! Une seconde lecture permet, une fois lon_t et lat_t lues et connues A la premiere lecture, !12-12-18
! d'en deduire l'orientation du canal et d'obtenir les tableaux de rotation
! https://docs.google.com/presentation/d/1FmAXNCdY_vL5KCUvkY0AxOQ6u145xQc1q-3cfXL1shU/edit#slide=id.p
!      open(unit=unit_,file=trim(riverchannelfile(k))) !06-12-18
!1290   read(unit_,*,end=1286)i1,j1
!       i=i1-par%timax(1) ; j=j1-par%tjmax(1) 
!         if(i>=1.and.i<=imax.and.j>=1.and.j<=jmax) then !>>>

!           if(riverdir(k)==1.or.riverdir(k)==3) then !along i>
!            x1=lon_t(i+1,j)-lon_t(i-1,j) ; if(x1<-pi)x1=x1+2.*pi ; if(x1> pi)x1=x1-2.*pi
!            x0=-atan2(lat_t(i+1,j)-lat_t(i-1,j),x1*cos(lat_t(i,j)))
!           endif                                     !along i>
!           if(riverdir(k)==2.or.riverdir(k)==4) then !along j>
!            x1=lon_t(i,j+1)-lon_t(i,j-1) ; if(x1<-pi)x1=x1+2.*pi ; if(x1> pi)x1=x1-2.*pi
!            x0=-atan2(lat_t(i,j+1)-lat_t(i,j-1),x1*cos(lat_t(i,j)))+0.5*pi 
!           endif                                     !along j>
!           gridrotcos_t(i,j)=cos(x0) ; gridrotsin_t(i,j)=sin(x0)

!         endif                                          !>>>
!       goto 1290
!1286  close(unit_)

       else                                 !moLom> !07-12-18

! CAS STANDARD utilisant les parametres du notebook_river: DEBUT

! i2,j2 est le point en aval du canal (premier point en mer)
         i2=iriver(k,3) ; j2=jriver(k,3)
         if(riverdir(k)==1)i2=i2+l_river(k)+1
         if(riverdir(k)==3)i2=i2-l_river(k)-1
         if(riverdir(k)==2)j2=j2+l_river(k)+1
         if(riverdir(k)==4)j2=j2-l_river(k)-1

        i1=iriver(k,3)+par%timax(1) ; j1=jriver(k,3)+par%tjmax(1) !Coord point source en amont du canal   

!       write(10+par%rank,*)'trouvemoi',i2,j2,i2+par%timax(1),j2+par%tjmax(1)

        k0=0 
        if(i2>1.and.i2<imax.and.j2>1.and.j2<jmax)k0=1
        call mpi_allreduce(k0,k1,1,mpi_integer,mpi_sum,par%comm2d,ierr) !Faire connaitre k0 A tous les dom
        if(k1>1) stop 'Err 1186 set_rivers_reservoir_grid' !Debug si k1>1 echec du test de non chevauchement sur 2 dom
        if(k1<1) stop 'Err 1187 set_rivers_reservoir_grid' 

! Faire connaitre dx_t(i2,j2) A tous les dom:
        x1=0.
        if(i2>1.and.i2<imax.and.j2>1.and.j2<jmax)x1=dx_t(i2,j2)
        call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        if(riverdir(k)==1.or.riverdir(k)==3)alongresolriver(k,2)=x2
        if(riverdir(k)==2.or.riverdir(k)==4)crossresolriver(k,2)=x2

! Faire connaitre dy_t(i2,j2) A tous les dom:
        y1=0.
        if(i2>1.and.i2<imax.and.j2>1.and.j2<jmax)y1=dy_t(i2,j2)
        call mpi_allreduce(y1,y2,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        if(riverdir(k)==1.or.riverdir(k)==3)crossresolriver(k,2)=y2
        if(riverdir(k)==2.or.riverdir(k)==4)alongresolriver(k,2)=y2

!       write(10+par%rank,*)'trouvemoi crossresolriver',y2

! note: dx(i)=x10*dx(i-1) soit dx(i)=dx(1)*x10**(i-1) !11-11-18
         x10=(alongresolriver(k,2)/alongresolriver(k,1))**(1./l_river(k)) ! alongriver loi de proportionalite

! Estimation de la longueur du canal
! Note: dx(i)=dx(i-1)*x10 ==>  SUM(dx)=dx(1)*(1+x10+x10**2+....+x10**(l_river-1))
         sum1=0.
         do k10=1,l_river(k)
          sum1=sum1+x10**(k10-1)
         enddo
         sum1=sum1*alongresolriver(k,1)
         if(par%rank==0)write(3,*)'longueur en km du canal de ',trim(rivername(k)),'=',sum1/1000.

! Attention les dx,dy du point de jonction avec l'ocean restent inchangEs
! La modification ne commence qu'au premier point A l'interieur du reservoir


        if(riverdir(k)==1) then !-cas i croissant->
         j=j1 ; i2=i1+l_river(k)+1  !11-01-23

         do i=i1,i2-1
          
          rap=real(i2-i)/real(l_river(k)+1) !11-01-23

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then

! Methode 1: progressive lineaire (de la resolution par default A la
! resolution souhaitee au dernier point du reservoir)
!          dx_t(i3,j3)=(1-rap)*x2+rap*alongresolriver(k,1)
           dx_t(i3,j3)=alongresolriver(k,1)*x10**(i-i1)     ! alongriver loi de proportionalite !11-11-18
           dy_t(i3,j3)=(1-rap)*y2+rap*crossresolriver(k,1)  ! crossriver loi lineaire
! Methode 2: tout le reservoir (sauf premier point) A la rEsolution souhaitEe
!          dx_t(i3,j3)=alongresolriver(k,1) ; dy_t(i3,j3)=x1
!       write(10+par%rank,*)'trouvebibi',i3+par%timax(1),j3+par%tjmax(1),dy_t(i3,j3)

          endif

         enddo
        endif                   !-cas i croissant->

        if(riverdir(k)==3) then !-cas i decroissant->
         j=j1 ; i2=i1-l_river(k)-1  !11-01-23

         do i=i2+1,i1
          
          rap=real(i-i2)/real(l_river(k)+1) !11-01-23

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then

! Methode 1: progressive lineaire (de la resolution par default A la
! resolution souhaitee au dernier point du reservoir)
!          dx_t(i3,j3)=(1-rap)*x2+rap*alongresolriver(k,1)
           dx_t(i3,j3)=alongresolriver(k,1)*x10**(i1-i)     ! alongriver loi de proportionalite !11-11-18
           dy_t(i3,j3)=(1-rap)*y2+rap*crossresolriver(k,1)
! Methode 2: tout le reservoir (sauf premier point) A la rEsolution souhaitEe
!          dx_t(i3,j3)=alongresolriver(k,1) ; dy_t(i3,j3)=x1

          endif

         enddo
        endif                   !-cas i decroissant->

        if(riverdir(k)==2) then !-cas j croissant->
         i=i1 ; j2=j1+l_river(k)+1 !11-01-23

         do j=j1,j2-1
          
          rap=real(j2-j)/real(l_river(k)+1) !11-01-23

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then

! Methode 1: progressive lineaire (de la resolution par default A la
! resolution souhaitee au dernier point du reservoir)
           dx_t(i3,j3)=(1-rap)*x2+rap*crossresolriver(k,1)
!          dy_t(i3,j3)=(1-rap)*y2+rap*alongresolriver(k,1)
           dy_t(i3,j3)=alongresolriver(k,1)*x10**(j-j1)     ! alongriver loi de proportionalite !11-11-18
! Methode 2: tout le reservoir (sauf premier point) A la rEsolution souhaitEe
!          dx_t(i3,j3)=alongresolriver(k,1) ; dy_t(i3,j3)=x1

          endif

         enddo
        endif                   !-cas j croissant->

        if(riverdir(k)==4) then !-cas j decroissant->
         i=i1 ; j2=j1-l_river(k)-1 !11-01-23

         do j=j2+1,j1
          
          rap=real(j-j2)/real(l_river(k)+1) !11-01-23

          i3=i-par%timax(1) ; j3=j-par%tjmax(1) 
          if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then

! Methode 1: progressive lineaire (de la resolution par default A la
! resolution souhaitee au dernier point du reservoir)
           dx_t(i3,j3)=(1-rap)*x2+rap*crossresolriver(k,1)
!          dy_t(i3,j3)=(1-rap)*y2+rap*alongresolriver(k,1)
           dy_t(i3,j3)=alongresolriver(k,1)*x10**(j1-j)     ! alongriver loi de proportionalite !11-11-18
! Methode 2: tout le reservoir (sauf premier point) A la rEsolution souhaitEe
!          dx_t(i3,j3)=alongresolriver(k,1) ; dy_t(i3,j3)=x1

          endif

         enddo
        endif                   !-cas j decroissant->


       endif                                !moLom> !07-12-18

       endif                  !lllll>
     
       if(par%rank==0)close(3) ! tmp/messages
       endif                                                       !m°v°m> !11-11-18
      enddo         !RRRRRRRR>

      end subroutine set_rivers_reservoir_grid
