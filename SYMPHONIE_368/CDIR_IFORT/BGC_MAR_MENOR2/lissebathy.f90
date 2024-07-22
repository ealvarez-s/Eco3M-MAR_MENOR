










      subroutine lissebathy(case_,rmax_,ideb_,ifin_,jdeb_,jfin_,key_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release S.26  - last update: 27-03-18
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
      real :: rmax_            &
             ,overslope_      &! taux de depassement de la pente par la pseudo pente
             ,rmx_
      integer case_,ideb_,ifin_,jdeb_,jfin_,key_
      integer i1_,j1_,i2_,j2_
      integer :: ipnocglob,loopmax_


!..............................................................................
! Version date      Description des modifications
!         08/10/01: KO part de 0 et non de 1
!         28/05/03: lissage de la bathy appliqué egalement dans le masque.
!                   raison: la methode MPV demasque les iles, il est donc
!                   important que la bathy soit lisse là aussi.
!         06/08/03: blindage pour eviter des divisions par zéro
!         24/03/04: ajout d'un choix n°2 pour lissage des niveaux z
!         18/04/06: fonctions compatibles avec double precision
!         10/05/06: fonctions compatibles avec double precision
!         24/03/09: parallelisation
!         30/03/09: Ne pas modifier la valeur des indices passés en arguments
!         01/04/09: modification pour parallelisation. creation variable
!                   globale IPNOCGLOB
! 2009.3  05-10-09: ajout d'un "ifdef parallele"
!         24-11-09: - suppression case_=2
!                   - pour parallelisation parfaite on ne lisse pas dans le masque
!         06-12-09: 0.125 -> 0.125d0
! 2010.8  09-05-10  seul le proc 0 ecrit messages
! 2010.12 20-09-10  Possibilité de calcul en simple precision
! 2010.13 25-10-10  Apres s'etre rendu compte que diffuser H sur la côte pouvait
!                   induire de tres fortes pentes (par ex si le talus touche la côte)
!                   on applique une C.L. de diffusion nulle à la côte
! S.26    17-07-13  lissage: attention à ponderation par 1/h quand h proche de zero
!                   voir h<0
!         14-04-14  cas case_=1, convertir les indices globlaux en indices locaux
!         14-07-14  key devient key_
!         22-10-14  Si rmax_ n'est pas donne dans notebook_bathy (val<0) rmax_=1./(kmax)
!         22-11-14  evite test sur mask
!         10-01-15  suppression mask dans lissage cas 1
!         02-02-15  sigma/step construite a partir du lissage de H
!         25-05-15  ajout case_2
!         02-07-15  Generation de la grille hybride: possibilite de reprendre
!                   le processus de convergence depuis la lecture d'un fichier
!         07-07-15  Suite du point precedent
!         12-07-15  sum0 evalue la convergence de l'algo de base egalement
!         09-11-15  modif de l'algo de lissage case 3 dans les zones ou
!                   h tend vers zero
!         25-11-15  modif lissage cas 3
!         12-02-16  h2min_ devient hrmax lu dans notebook_bathy
!         03-04-16  routine re-structurEe en subroutines
!         29-06-16  possibilite lissage cas 3 sur fentre i,j reduite
!...............................................................................
!    _________                    .__                  .__             ! m[0v0]m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!...............................................................................

      i1_=ideb_ ; i2_=ifin_ ; j1_=jdeb_ ; j2_=jfin_

! Convertir les indices globaux en indices locaux !29-06-16
      i1_=i1_-par%timax(1)  
      i2_=i2_-par%timax(1) 
      j1_=j1_-par%tjmax(1)  
      j2_=j2_-par%tjmax(1)  

      i1_=max(i1_,1   )        
      i2_=min(i2_,imax)
      j1_=max(j1_,1   )
      j2_=min(j2_,jmax)

      if(case_==1) then !11111>
       if(flag_smooth_h_mask==0)call lissebathy_regular(i1_,i2_,j1_,j2_,key_)
       if(flag_smooth_h_mask==1)call lissebathy_regular_mask(i1_,i2_,j1_,j2_,key_)
      endif             !11111>

      if(case_==2)call lissebathy_vst

      if(case_==3) then !33333>
       if(flag_smooth_h_mask==0)call lissebathy_rmax(i1_,i2_,j1_,j2_,rmax_)
       if(flag_smooth_h_mask==1)call lissebathy_rmax_mask(i1_,i2_,j1_,j2_,rmax_)
      endif             !33333>

      end subroutine lissebathy

!.................................................................

      subroutine lissebathy_rmax(i1_,i2_,j1_,j2_,rmax_)
      use module_principal
      use module_parallele !#MPI
      implicit none
      real rmax_           &
          ,overslope_      &! taux de depassement de la pente par la pseudo pente
          ,rmx_
      integer :: ipnocglob,loopmax_
      integer i1_,j1_,i2_,j2_

      if(ihybsig==1)then !>>>>>>>>
        ub2=ubound(h_w) ; lb2=lbound(h_w)
        if(.not.allocated(h0_w))allocate(h0_w(lb2(1):ub2(1),lb2(2):ub2(2)))
        h0_w(:,:)=h_w(:,:)
      endif              !>>>>>>>>

!     OPEN(UNIT=3,FILE=trim(tmpdirname)//'messages',ACCESS='APPEND')

! APPLICATION DU CRITERE dH/2H < RMAX
! (voir article Beckman and Haidvogel JPO 1993 pp 1736-1753)
! http://dx.doi.org/10.1175/1520-0485(1993)023<1736:NSOFAA>2.0.CO;2
! (voir article Mellor & cie JAOT 1994 pp 1126-1134)
! RMAX est compris entre 0 et 1

! Explications: la couchle la plus sensible a l'erreur sigma est
! a priori la couche du fond. Soientt zw(i,k=1) et zw(i,k=2) les profondeurs
! aux 2 facettes de la couche du fond. Il y a incoherence hydrostatique
! si zw(i-1,k=1) > zw(i,k=2)
! On a d'une part zw(i-1,1)=-h(i-1) et d'autre part
! dans le cas d'une coordonnee sigma uniforme et N couches z(i,2)=z(i,1)+h(i)/N=-h(i)+h(i)/N
! La condition zw(i-1,k=1) > zw(i,k=2) est equivalente a -h(i-1)>-h(i)+h(i)/N
! autrement dit (h(i)-h(i-1))/h(i) < 1/N soit rmax_=1/N
! Si rmax_ n'est pas donne dans notebook_bathy (valeur<0) rmax_=1/kmax
      if(rmax_<0)rmax_=1./real(kmax)                                     !22-10-14

      write(texte30,'(a,i0)')'ihybsig=',ihybsig
      if(ihybsig==1) then !>>>>>> !07-07-15
       texte30='Cold start'
       write(texte60,'(a,i0)') &
       trim(restartdirname)//'restart_input/lissebathy',par%rank
       open(unit=9,file=texte60,access='sequential',form='unformatted',err=155) !02-07-15
       if(par%rank==0)write(6,'(a,a)')'UTILISATION DU FICHIER RESTART ',trim(texte60)
       read(9,end=155)h_w
       texte30='Hot restart'
  155  close(9)
      endif               !>>>>>>
     
      if(ihybsig==1) then
       loopmax_=10000
      else
       loopmax_=10000
      endif

      do 111 ko=0,loopmax_
      overslope_=1.  ! 1.1 ! >=1 taux de depassement de la pente par la pseudo pente
      ipnoc=0

      do 101 j=j1_,j2_ ! 1,jmax !29-06-16
      do 101 i=i1_,i2_ ! 1,imax

      x1=min(0.125*un,max(zero,                                       & !06-12-09
            (abs(h_w(i,j)-h_w(i-1,j))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i-1,j) ,hrmax)                         & !25-11-15
!       /max(abs(h_w(i,j)+h_w(i-1,j)+hrmax),small1)                  & !09-11-15
         -rmax_)/(rmax_+small1))) !*mask_t(i-1,j,kmax+1) !10-01-15

      x2=min(0.125*un,max(zero,                                       &
            (abs(h_w(i,j)-h_w(i+1,j))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i+1,j) ,hrmax)                        &
!       /max(abs(h_w(i,j)+h_w(i+1,j)+hrmax),small1)                  &
         -rmax_)/(rmax_+small1))) !*mask_t(i+1,j,kmax+1)                                           !24-11-09

      x3=min(0.125*un,max(zero,                                       &
            (abs(h_w(i,j)-h_w(i,j-1))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i,j-1) ,hrmax)                  &
!       /max(abs(h_w(i,j)+h_w(i,j-1)+hrmax),small1)                  &
         -rmax_)/(rmax_+small1))) !*mask_t(i,j-1,kmax+1)                                            !24-11-09

      x4=min(0.125*un,max(zero,                                       &
            (abs(h_w(i,j)-h_w(i,j+1))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i,j+1),hrmax)                  &
!       /max(abs(h_w(i,j)+h_w(i,j+1)+hrmax),small1)                  &
         -rmax_)/(rmax_+small1))) !*mask_t(i,j+1,kmax+1)                                            !24-11-09

      if(ihybsig==1) then !>>>>>

! x11 doit etre negatif sinon on lisse:
         x11= (h_w(i,j) -h_w(i-1,j))             &
            *( h_w(i,j) -h_w(i-1,j)              &
            -(h0_w(i,j)-h0_w(i-1,j))*overslope_) &
      /max(abs(h_w(i,j) +h_w(i-1,j)),hrmax)

       x1=max(x1,x11)
       x1=min(0.125d0,max(0d0,x1))

! x12 doit etre negatif sinon on lisse:
           x12= (h_w(i,j) -h_w(i+1,j))             &
              *( h_w(i,j) -h_w(i+1,j)              &
              -(h0_w(i,j)-h0_w(i+1,j))*overslope_) &
        /max(abs(h_w(i,j) +h_w(i+1,j)),hrmax)
       x2=max(x2,x12)
       x2=min(0.125d0,max(0d0,x2))

! x13 doit etre negatif sinon on lisse:
           x13= (h_w(i,j) -h_w(i,j-1))             &
              *( h_w(i,j) -h_w(i,j-1)              &
              -(h0_w(i,j)-h0_w(i,j-1))*overslope_) &
        /max(abs(h_w(i,j) +h_w(i,j-1)),hrmax)

       x3=max(x3,x13)
       x3=min(0.125d0,max(0d0,x3))

! x14 doit etre negatif sinon on lisse:
           x14= (h_w(i,j) -h_w(i,j+1))             &
              *( h_w(i,j) -h_w(i,j+1)              &
              -(h0_w(i,j)-h0_w(i,j+1))*overslope_) &
        /max(abs(h_w(i,j) +h_w(i,j+1)),hrmax)

       x4=max(x4,x14)
       x4=min(0.125d0,max(0d0,x4))

      endif               !>>>>>

      xy_t(i,j,1)=(1.-x1-x2-x3-x4)*h_w(i  ,j  )                         &
                     +x1          *h_w(i-1,j  )                         &
                        +x2       *h_w(i+1,j  )                         &
                           +x3    *h_w(i  ,j-1)                         &
                              +x4 *h_w(i  ,j+1)

      if(x1+x2+x3+x4>small1)ipnoc=ipnoc+1

!     endif                      !-----------------<
  101 continue

      if(ihybsig/=1) then !>>>>>>

! Algo base
       sum0=0.
!      do j=1,jmax ; do i=1,imax
       do j=j1_,j2_  ; do i=i1_,i2_
          sum0=sum0+abs(h_w(i,j)-xy_t(i,j,1))*mask_i_w(i) & !12-07-15
                                             *mask_j_w(j)
                        h_w(i,j)=xy_t(i,j,1)  !10-01-15
       enddo ; enddo

      else                !>>>>>>

! Algo VST
!      x1=4./loopmax_ ! coefficient du terme de rappel vers h0
       x1=0.
       sum0=0.
!      do j=1,jmax ; do i=1,imax
       do j=j1_,j2_  ; do i=i1_,i2_
          sum0=sum0+abs(h_w(i,j)-xy_t(i,j,1))*mask_i_w(i) &
                                             *mask_j_w(j)
!       h_w(i,j)=max(xy_t(i,j,1),h_w(i,j))*(1.-x1)+x1*h0_w(i,j)
        h_w(i,j)=max( xy_t(i,j,1)*(1.-x1)   &
                     +h0_w(i,j)      *x1    & ! rappel vers h0
                    , h0_w(i,j) )


       enddo ; enddo

      endif               !>>>>>>

      call obc_h(0)                                                    !24/03/09


      call mpi_allreduce(sum0,sum0glb,1,mpi_double,                & !01/04/09
           mpi_sum,par%comm2d ,ierr)
      call mpi_allreduce(ipnoc,ipnocglob,1,mpi_integer,                & !01/04/09
           mpi_sum,par%comm2d ,ierr)
      ipnoc=ipnocglob
      sum0=sum0glb
      if(par%rank==0.and.mod(ko,100)==0) &
!     write(6,'(a,i0,1x,i0,1x,a)')'ko ipnoc=',ko,ipnoc,texte30
      write(6,*)'ko sum0=',ko,sum0,texte30
      if(ipnoc==0)goto 120


  111 continue

  120 continue

! Faire un fichier restart pour reprendre la convergence
      if(ihybsig==1) then !>>>>>> !07-07-15
       write(texte60,'(a,i0)') &
       trim(restartdirname)//'restart_input/lissebathy',par%rank
       open(unit=9,file=texte60,access='sequential',form='unformatted')
       write(9)h_w
       close(9)
      endif               !>>>>>>

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine lissebathy:'
      write(3,*)
      write(3,*)'rmax_ et nombre de lissages correspondants:',rmax_,ko
      close(3)
      endif                !#mpi-->>-->                       !09-05-10

!     j=jmax/2
!     do i=1,imax
!     write(20+par%rank,*)i,h0_w(i,j),h_w(i,j)
!     enddo
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'poulpe'

      end subroutine lissebathy_rmax

!.......................................................................

      subroutine lissebathy_vst
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer :: ipnocglob,loopmax_
      real :: rmx_

! Cet algo lisse la bathy des grille V.S.T dans la partie
! ou l'algo V.S.T est desactive, c'est a dire dans les zones
! ou le nombre de niveaux verticaux est borne par nbvstepmin
! rmax_ est indexe sur nbvstepmin

      k0=min(kmax,nbvstepmin)
      rmx_=1./real(k0)

      hrmax=1.
      do ko=0,1000 

      ipnoc=0

      do j=1,jmax ; do i=1,imax

      x1=min(0.125*un,max(zero,                                       & !06-12-09
        (abs(h_w(i,j)-h_w(i-1,j))                                     & !28/05/03
        /max(h_w(i,j)+h_w(i-1,j),hrmax)                               &
         -rmx_)/(rmx_+small1))) !*mask_t(i-1,j,kmax+1) !10-01-15

      x2=min(0.125*un,max(zero,                                       &
        (abs(h_w(i,j)-h_w(i+1,j))                                     & !28/05/03
        /max(h_w(i,j)+h_w(i+1,j),hrmax)                                  &
         -rmx_)/(rmx_+small1))) !*mask_t(i+1,j,kmax+1)                                           !24-11-09

      x3=min(0.125*un,max(zero,                                       &
        (abs(h_w(i,j)-h_w(i,j-1))                                     & !28/05/03
        /max(h_w(i,j)+h_w(i,j-1),hrmax)                                  &
         -rmx_)/(rmx_+small1))) !*mask_t(i,j-1,kmax+1)                                            !24-11-09

      x4=min(0.125*un,max(zero,                                       &
        (abs(h_w(i,j)-h_w(i,j+1))                                     & !28/05/03
        /max(h_w(i,j)+h_w(i,j+1),hrmax)                                  &
         -rmx_)/(rmx_+small1))) !*mask_t(i,j+1,kmax+1)                                            !24-11-09


      xy_t(i,j,1)=(1.-x1-x2-x3-x4)*h_w(i  ,j  )                         &
                     +x1          *h_w(i-1,j  )                         &
                        +x2       *h_w(i+1,j  )                         &
                           +x3    *h_w(i  ,j-1)                         &
                              +x4 *h_w(i  ,j+1)

      if(x1+x2+x3+x4>small1)ipnoc=ipnoc+1

      enddo ; enddo ! i j

       k0=min(kmax,nbvstepmin)
       sum0=0.
       do j=1,jmax ; do i=1,imax
        if(kmax+1-kmin_w(i,j)<k0) then
          sum0=sum0+abs(h_w(i,j)-xy_t(i,j,1))*mask_i_w(i) &
                                             *mask_j_w(j)
          h_w(i,j)=xy_t(i,j,1)
        endif
       enddo ; enddo

      call obc_h(0)                                                    !24/03/09

      call mpi_allreduce(ipnoc,ipnocglob,1,mpi_integer,                & !01/04/09
           mpi_sum,par%comm2d ,ierr)
      call mpi_allreduce(sum0,sum0glb,1,mpi_double,                & !01/04/09
           mpi_sum,par%comm2d ,ierr)
      ipnoc=ipnocglob
      sum0=sum0glb
      if(par%rank==0.and.mod(ko,100)==0) then
         write(6,*)'lissebathy case2 ',ko,' sum0=',sum0
      endif
      if(ipnoc==0)goto 57


      enddo ! ko loop

   57 continue

      end subroutine lissebathy_vst

!.................................................................

      subroutine lissebathy_regular(i1_,i2_,j1_,j2_,key_)
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer key_
      integer i1_,j1_,i2_,j2_

!*****************************************************************
! Five-Point Laplacian smoother
! Début:
!*****************************************************************

      call obc_h(0)                                                    !24/03/09

      do 12 k=1,key_
      do 10 j=j1_,j2_
      do 10 i=i1_,i2_

        xy_t(i,j,2)=0.9* h_w(i  ,j  )          & !10-01-15
                   +0.1*(h_w(i  ,j-1)+         &
                         h_w(i  ,j+1)+         &
                         h_w(i-1,j  )+         &
                         h_w(i+1,j  ) )*0.25


!     if(mask_t(i,j,kmax+1)==1) then !----------------->                   !24-11-09
!     x1=mask_t(i-1,j  ,kmax+1)*0.025                      !25-10-10
!     x2=mask_t(i+1,j  ,kmax+1)*0.025
!     x3=mask_t(i  ,j-1,kmax+1)*0.025
!     x4=mask_t(i  ,j+1,kmax+1)*0.025
!     xy_t(i,j,2)=(1.-x1-x2-x3-x4)*h_w(i  ,j  )          & !25-10-10
!                    +x1          *h_w(i-1,j  )          &
!                       +x2       *h_w(i+1,j  )          &
!                          +x3    *h_w(i  ,j-1)          &
!                             +x4 *h_w(i  ,j+1)
!     endif                      !-----------------<

   10 continue

      do 11 j=j1_,j2_
      do 11 i=i1_,i2_
       h_w(i,j)=xy_t(i,j,2) !10-01-15
!      h_w(i,j)=xy_t(i,j,2) *mask_t(i,j,kmax) & !22-11-14
!               +h_w(i,j)*(1-mask_t(i,j,kmax))
   11 continue

      call obc_h(0)                                                    !24/03/09

   12 continue

!*****************************************************************
! Five-Point Laplacian smoother
! Fin.
!*****************************************************************

      end subroutine lissebathy_regular

!.................................................................

      subroutine lissebathy_regular_mask(i1_,i2_,j1_,j2_,key_)
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer key_
      integer i1_,j1_,i2_,j2_
      real rmax_           &
          ,overslope_      &! taux de depassement de la pente par la pseudo pente
          ,rmx_

! Comme lissebathy_regular mais la bathy dans le masque ne compte pas

!*****************************************************************
! Five-Point Laplacian smoother
! Début:
!*****************************************************************

      call obc_h(0)                                                    !24/03/09

      do 12 k=1,key_
      do 10 j=j1_,j2_
      do 10 i=i1_,i2_

      x1=mask_t(i-1,j  ,kmax)*0.025                      !25-10-10
      x2=mask_t(i+1,j  ,kmax)*0.025
      x3=mask_t(i  ,j-1,kmax)*0.025
      x4=mask_t(i  ,j+1,kmax)*0.025
      xy_t(i,j,2)=(1.-x1-x2-x3-x4)*h_w(i  ,j  )          & !25-10-10
                     +x1          *h_w(i-1,j  )          &
                        +x2       *h_w(i+1,j  )          &
                           +x3    *h_w(i  ,j-1)          &
                              +x4 *h_w(i  ,j+1)

   10 continue

      do 11 j=j1_,j2_
      do 11 i=i1_,i2_
       h_w(i,j)=xy_t(i,j,2) *mask_t(i,j,kmax) & !22-11-14
                +h_w(i,j)*(1-mask_t(i,j,kmax))
   11 continue

      call obc_h(0)                                                    !24/03/09

   12 continue

!*****************************************************************
! Five-Point Laplacian smoother
! Fin.
!*****************************************************************

      end subroutine lissebathy_regular_mask

!.................................................................

      subroutine lissebathy_rmax_mask(i1_,i2_,j1_,j2_,rmax_)
      use module_principal
      use module_parallele !#MPI
      implicit none
      real rmax_           &
          ,overslope_      &! taux de depassement de la pente par la pseudo pente
          ,rmx_
      integer :: ipnocglob,loopmax_
      integer i1_,j1_,i2_,j2_

      if(ihybsig==1)then !>>>>>>>>
        ub2=ubound(h_w) ; lb2=lbound(h_w)
        if(.not.allocated(h0_w))allocate(h0_w(lb2(1):ub2(1),lb2(2):ub2(2)))
        h0_w(:,:)=h_w(:,:)
      endif              !>>>>>>>>

!     OPEN(UNIT=3,FILE=trim(tmpdirname)//'messages',ACCESS='APPEND')

! APPLICATION DU CRITERE dH/2H < RMAX
! (voir article Beckman and Haidvogel JPO 1993 pp 1736-1753)
! http://dx.doi.org/10.1175/1520-0485(1993)023<1736:NSOFAA>2.0.CO;2
! (voir article Mellor & cie JAOT 1994 pp 1126-1134)
! RMAX est compris entre 0 et 1

! Explications: la couchle la plus sensible a l'erreur sigma est
! a priori la couche du fond. Soientt zw(i,k=1) et zw(i,k=2) les profondeurs
! aux 2 facettes de la couche du fond. Il y a incoherence hydrostatique
! si zw(i-1,k=1) > zw(i,k=2)
! On a d'une part zw(i-1,1)=-h(i-1) et d'autre part
! dans le cas d'une coordonnee sigma uniforme et N couches z(i,2)=z(i,1)+h(i)/N=-h(i)+h(i)/N
! La condition zw(i-1,k=1) > zw(i,k=2) est equivalente a -h(i-1)>-h(i)+h(i)/N
! autrement dit (h(i)-h(i-1))/h(i) < 1/N soit rmax_=1/N
! Si rmax_ n'est pas donne dans notebook_bathy (valeur<0) rmax_=1/kmax
      if(rmax_<0)rmax_=1./real(kmax)                                     !22-10-14

      write(texte30,'(a,i0)')'ihybsig=',ihybsig
      if(ihybsig==1) then !>>>>>> !07-07-15
       texte30='Cold start'
       write(texte60,'(a,i0)') &
       trim(restartdirname)//'restart_input/lissebathy',par%rank
       open(unit=9,file=texte60,access='sequential',form='unformatted',err=155) !02-07-15
       if(par%rank==0)write(6,'(a,a)')'UTILISATION DU FICHIER RESTART ',trim(texte60)
       read(9,end=155)h_w
       texte30='Hot restart'
  155  close(9)
      endif               !>>>>>>
     
      if(ihybsig==1) then
       loopmax_=10000
      else
       loopmax_=10000
      endif

      do 111 ko=0,loopmax_
      overslope_=1.  ! 1.1 ! >=1 taux de depassement de la pente par la pseudo pente
      ipnoc=0

      do 101 j=j1_,j2_ ! 1,jmax
      do 101 i=i1_,i2_ ! 1,imax
      if(mask_t(i,j,kmax)==1) then !>>>>>>

      x1=min(0.125*un,max(zero,                                       & !06-12-09
            (abs(h_w(i,j)-h_w(i-1,j))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i-1,j) ,hrmax)                         & !25-11-15
!       /max(abs(h_w(i,j)+h_w(i-1,j)+hrmax),small1)                  & !09-11-15
         -rmax_)/(rmax_+small1))) !*mask_t(i-1,j,kmax+1) !10-01-15

      x2=min(0.125*un,max(zero,                                       &
            (abs(h_w(i,j)-h_w(i+1,j))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i+1,j) ,hrmax)                        &
!       /max(abs(h_w(i,j)+h_w(i+1,j)+hrmax),small1)                  &
         -rmax_)/(rmax_+small1))) !*mask_t(i+1,j,kmax+1)                                           !24-11-09

      x3=min(0.125*un,max(zero,                                       &
            (abs(h_w(i,j)-h_w(i,j-1))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i,j-1) ,hrmax)                  &
!       /max(abs(h_w(i,j)+h_w(i,j-1)+hrmax),small1)                  &
         -rmax_)/(rmax_+small1))) !*mask_t(i,j-1,kmax+1)                                            !24-11-09

      x4=min(0.125*un,max(zero,                                       &
            (abs(h_w(i,j)-h_w(i,j+1))                                 & !28/05/03
        /max(    h_w(i,j)+h_w(i,j+1),hrmax)                  &
!       /max(abs(h_w(i,j)+h_w(i,j+1)+hrmax),small1)                  &
         -rmax_)/(rmax_+small1))) !*mask_t(i,j+1,kmax+1)                                            !24-11-09

      if(ihybsig==1) then !>>>>>

! x11 doit etre negatif sinon on lisse:
         x11= (h_w(i,j) -h_w(i-1,j))             &
            *( h_w(i,j) -h_w(i-1,j)              &
            -(h0_w(i,j)-h0_w(i-1,j))*overslope_) &
      /max(abs(h_w(i,j) +h_w(i-1,j)),hrmax)

       x1=max(x1,x11)
       x1=min(0.125d0,max(0d0,x1))

! x12 doit etre negatif sinon on lisse:
           x12= (h_w(i,j) -h_w(i+1,j))             &
              *( h_w(i,j) -h_w(i+1,j)              &
              -(h0_w(i,j)-h0_w(i+1,j))*overslope_) &
        /max(abs(h_w(i,j) +h_w(i+1,j)),hrmax)
       x2=max(x2,x12)
       x2=min(0.125d0,max(0d0,x2))

! x13 doit etre negatif sinon on lisse:
           x13= (h_w(i,j) -h_w(i,j-1))             &
              *( h_w(i,j) -h_w(i,j-1)              &
              -(h0_w(i,j)-h0_w(i,j-1))*overslope_) &
        /max(abs(h_w(i,j) +h_w(i,j-1)),hrmax)

       x3=max(x3,x13)
       x3=min(0.125d0,max(0d0,x3))

! x14 doit etre negatif sinon on lisse:
           x14= (h_w(i,j) -h_w(i,j+1))             &
              *( h_w(i,j) -h_w(i,j+1)              &
              -(h0_w(i,j)-h0_w(i,j+1))*overslope_) &
        /max(abs(h_w(i,j) +h_w(i,j+1)),hrmax)

       x4=max(x4,x14)
       x4=min(0.125d0,max(0d0,x4))

      endif               !>>>>>

      x1=x1*mask_t(i-1,j  ,kmax)
      x2=x2*mask_t(i+1,j  ,kmax)
      x3=x3*mask_t(i  ,j-1,kmax)
      x4=x4*mask_t(i  ,j+1,kmax)

      xy_t(i,j,1)=(1.-x1-x2-x3-x4)*h_w(i  ,j  )                         &
                     +x1          *h_w(i-1,j  )                         &
                        +x2       *h_w(i+1,j  )                         &
                           +x3    *h_w(i  ,j-1)                         &
                              +x4 *h_w(i  ,j+1)

      if(x1+x2+x3+x4>small1)ipnoc=ipnoc+1

      endif                        !>>>>>>
  101 continue

      if(ihybsig/=1) then !>>>>>>

! Algo base
       sum0=0.
!      do j=1,jmax ; do i=1,imax
       do j=j1_,j2_ ; do i=i1_,i2_
        if(mask_t(i,j,kmax)==1) then !>>>>>>
          sum0=sum0+abs(h_w(i,j)-xy_t(i,j,1))*mask_i_w(i) & !12-07-15
                                             *mask_j_w(j)
                        h_w(i,j)=xy_t(i,j,1)  !10-01-15
        endif                        !>>>>>>
       enddo ; enddo

      else                !>>>>>>

! Algo VST
!      x1=4./loopmax_ ! coefficient du terme de rappel vers h0
       x1=0.
       sum0=0.
!      do j=1,jmax ; do i=1,imax
       do j=j1_,j2_ ; do i=i1_,i2_
       if(mask_t(i,j,kmax)==1) then !>>>>>>
          sum0=sum0+abs(h_w(i,j)-xy_t(i,j,1))*mask_i_w(i) &
                                             *mask_j_w(j)
!       h_w(i,j)=max(xy_t(i,j,1),h_w(i,j))*(1.-x1)+x1*h0_w(i,j)
        h_w(i,j)=max( xy_t(i,j,1)*(1.-x1)   &
                     +h0_w(i,j)      *x1    & ! rappel vers h0
                    , h0_w(i,j) )


       endif                        !>>>>>>
       enddo ; enddo

      endif               !>>>>>>

      call obc_h(0)                                                    !24/03/09


      call mpi_allreduce(sum0,sum0glb,1,mpi_double,                & !01/04/09
           mpi_sum,par%comm2d ,ierr)
      call mpi_allreduce(ipnoc,ipnocglob,1,mpi_integer,                & !01/04/09
           mpi_sum,par%comm2d ,ierr)
      ipnoc=ipnocglob
      sum0=sum0glb
      if(par%rank==0.and.mod(ko,100)==0) &
!     write(6,'(a,i0,1x,i0,1x,a)')'ko ipnoc=',ko,ipnoc,texte30
      write(6,*)'ko sum0=',ko,sum0,texte30
      if(ipnoc==0)goto 120


  111 continue

  120 continue

! Faire un fichier restart pour reprendre la convergence
      if(ihybsig==1) then !>>>>>> !07-07-15
       write(texte60,'(a,i0)') &
       trim(restartdirname)//'restart_input/lissebathy',par%rank
       open(unit=9,file=texte60,access='sequential',form='unformatted')
       write(9)h_w
       close(9)
      endif               !>>>>>>

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine lissebathy:'
      write(3,*)
      write(3,*)'rmax_ et nombre de lissages correspondants:',rmax_,ko
      close(3)
      endif                !#mpi-->>-->                       !09-05-10

!     j=50
!     do i=1,imax
!     write(20+par%rank,*)i,h0_w(i,j),h_w(i,j)
!     enddo
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'poulpe'

      end subroutine lissebathy_rmax_mask

!.......................................................................

      subroutine lissebathy_peaks_friendly(key_)
      use module_principal
      use module_parallele !#MPI
      implicit none
      integer key_

! Etape 1: identification des "sommets" A preserver. Cette etape
! s'obtient en identifiant les zones de faible derivees apres plusieurs
! lissages ayant supprimE les sommets secondaires 
      call obc_h(0)                                                    !24/03/09
      do j=1,jmax ; do i=1,imax
       xy_t(i,j,1)=h_w(i,j)
      enddo       ; enddo


! 50 lissages
      do 12 k=1,50

      do j=2,jmax-1 ; do i=2,imax-1

      x1=mask_t(i-1,j  ,kmax)*0.125                      !25-10-10
      x2=mask_t(i+1,j  ,kmax)*0.125
      x3=mask_t(i  ,j-1,kmax)*0.125
      x4=mask_t(i  ,j+1,kmax)*0.125
      xy_t(i,j,2)=(1.-x1-x2-x3-x4)*xy_t(i  ,j  ,1)          & !25-10-10
                     +x1          *xy_t(i-1,j  ,1)          &
                        +x2       *xy_t(i+1,j  ,1)          &
                           +x3    *xy_t(i  ,j-1,1)          &
                              +x4 *xy_t(i  ,j+1,1)

      enddo ; enddo

      do j=2,jmax-1 ; do i=2,imax-1

       xy_t(i,j,1)=xy_t(i,j,2) *mask_t(i,j,kmax) & !22-11-14
                   +h_w(i,j)*(1-mask_t(i,j,kmax))

      enddo ; enddo

      call obc_h_xy_t_z0(1)

   12 continue

! Les sommets sont associEs aux changements de signe de la derivee de
! cette bathymetrie lissee:
      do j=2,jmax-1 ; do i=2,imax-1
       if( (xy_t(i  ,j,1)-xy_t(i-1,j,1))    &
          *(xy_t(i+1,j,1)-xy_t(i  ,j,1))<0. &
       .or.(xy_t(i,j  ,1)-xy_t(i,j-1,1))    &
          *(xy_t(i,j+1,1)-xy_t(i,j  ,1))<0.) then !>>>
        xy_t(i,j,0)=0.
       else                                       !>>>
        xy_t(i,j,0)=1.
       endif                                      !>>>
      enddo         ; enddo
      call obc_h_xy_t_z0(0)

! Cette routine n'est pas terminee. Reste a refaire une etape de lissage
! utilisant xy_t(i,j,0) pour preserver les "sommets"
!..............................
! Bidouille pour visu:
!     do j=1,jmax ; do i=1,imax
!      z0_w(i,j)=xy_t(i,j,0)
!      cdb_t(i,j)=xy_t(i,j,1)
!     enddo       ; enddo
!     call graph_out
!     stop 'titi'
!..............................

!*****************************************************************
! Five-Point Laplacian smoother
! Fin.
!*****************************************************************

      end subroutine lissebathy_peaks_friendly

!.................................................................
