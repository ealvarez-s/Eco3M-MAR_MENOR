










      subroutine airseaflux_upd(ichoix)
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 17-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
      integer  ichoix

!...............................................................................
!Version  date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         19/03/02: amenagement etat initial similaire à obc_af.F et model_in.F
!         11/07/02: amenagements pour filiere formules bulk. Le 3eme argument
!                   des tableaux (0 1 ou 2) n'est plus codé en dur mais
!                   paramétré avec I0 et I2 pour eviter les warning à la
!                   compilation (avec formules bulk le 3eme argument est
!                   reduit à 1)
!         22/01/08: le renouvellement des echeances par modulo d'entier est
!                   supprime
! 2010.3  08-01-10: rotation calculee avec gridrotcos et gridrotsin
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.20 17-04-11  Calculs sur la base d'un temps en secondes
!...............................................................................


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PARTIE I
!     INITIALISATION DES FORCAGES EN DEBUT DE SIMULATION:
!     DEBUT:
      if (ichoix.eq.1) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(iairsea.ne.1)stop 'airseaflux_upd.f erreur1'
      if(ifb.ne.0)    stop 'airseaflux_upd.f erreur2'

      do k=1,nairsea                                                   !22/01/08
       i1=dateairsea(1,k)
       i2=dateairsea(2,k)
       i3=dateairsea(3,k)
       i4=dateairsea(4,k)
       i5=dateairsea(5,k)
       i6=dateairsea(6,k)
       call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!      airseadt(k,2)=xdtk_out
!      airseadt(k,1)=airseainfo(k,1)*3600./dti_fw
       airseadt(k,2)=elapsedtime_out          !17-04-11
       airseadt(k,1)=airseainfo(k,1)*3600.    !17-04-11
      enddo

      write(6,*)'le 22 jan 2008 j ai modifie la methode de'
      write(6,*)'renouvellement des echeances dans airseaflux_upd.f'
      write(6,*)'fini la methode des modulo d entier. le pb est'
      write(6,*)'que la validation de cette grosse modif n a pas'
      write(6,*)'été faite.'
      stop 'donc dans airseaflux_upd.f'

      i0=0                                                             !11/07/02
      i2=2                                                             !11/07/02

      do 1100 k1=1,6
       do 2000 k3=0,2,2

        x2=(elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)
        x1=(elapsedtime_bef-airseadt(k1,2))/airseadt(k1,1)
        j2=int(x2)
        j1=int(x1)
        if(x2.lt.0.)j2=j2-1
        if(x1.lt.0.)j1=j1-1
        nc=int(1.+x2)+k3/2
        if(k3.eq.2.and.(j2-j1).eq.1)nc=nc-1

        lrec=(imax+2)*(jmax+2)*4
        open(unit=3,file=airseafile(k1),access='direct',recl=lrec)
        read(3,rec=nc)anyvar2d
        close(3)

! Commentees le                                                        !22/01/08
!       X1= AIRSEAINFO(K1,2)*REAL(KOUNTMOD)*DTI_FW
! Remplacees par: le:                                                  !22/01/08
!       x1= airseainfo(k1,2)*dti_bd*airseadt(k1,1)                    &
        x1= airseainfo(k1,2)*       airseadt(k1,1)                    & !17-04-11
       +(1.-airseainfo(k1,2))*1.


!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.1)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,1)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.2)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,2)=anyvar2d(i,j)/x1
              enddo
              enddo

! Rotation de grille:   !08-01-10
      stop ' dans airseaflux_upd. Nouvelle rotation non verifiee'
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,1)=gridrotcos_t(i,j)*xy_t(i,j,1)-gridrotsin_t(i,j)*xy_t(i,j,2)
               xy_t(i,j,2)=gridrotsin_t(i,j)*xy_t(i,j,1)+gridrotcos_t(i,j)*xy_t(i,j,2)
              enddo
              enddo

! Regroupement aux points de vitesse:
              x0=windfactor*0.5
              do j=1,jmax+1
              do i=1,imax+1
               wstress_u(i,j,k3)=x0*(xy_t(i,j,1)+xy_t(i-1,j,1))
               wstress_v(i,j,k3)=x0*(xy_t(i,j,2)+xy_t(i,j-1,2))
              enddo
              enddo

            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.3)then
              do j=0,jmax+1
              do i=0,imax+1
               snsf_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.4)then
              do j=0,jmax+1
              do i=0,imax+1
               ssr_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.5)then
              do j=0,jmax+1
              do i=0,imax+1
               sshf_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.6)then
              do j=0,jmax+1
              do i=0,imax+1
               slhf_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !


 2000 continue
 1100 continue

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 écheances
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      k1=1                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 132 j=1,jmax+1
      do 132 i=1,imax+1
      wstress_u(i,j,1)=(1.-rap)*wstress_u(i,j,i0)+rap*wstress_u(i,j,i2)
      wstress_v(i,j,1)=(1.-rap)*wstress_v(i,j,i0)+rap*wstress_v(i,j,i2)
  132 continue

! Module de la tension de vent:
      do 133 j=1,jmax
      do 133 i=1,imax
      wstress_w(i,j)=sqrt(                                              &
       ((wstress_u(i,j,1)+wstress_u(i+1,j,1))/2.)**2+                   &
       ((wstress_v(i,j,1)+wstress_v(i,j+1,1))/2.)**2)
  133 continue

      k1=3                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 134 j=1,jmax
      do 134 i=1,imax
      snsf_w(i,j,1)=(1.-rap)*snsf_w(i,j,i0)+rap*snsf_w(i,j,i2)
  134 continue

      k1=4                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 135 j=1,jmax
      do 135 i=1,imax
      ssr_w(i,j,1)=(1.-rap)*ssr_w(i,j,i0)+rap*ssr_w(i,j,i2)
  135 continue

      k1=5                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 136 j=1,jmax
      do 136 i=1,imax
      sshf_w(i,j,1)=(1.-rap)*sshf_w(i,j,i0)+rap*sshf_w(i,j,i2)
  136 continue

      k1=6                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 137 j=1,jmax
      do 137 i=1,imax
      slhf_w(i,j,1)=(1.-rap)*slhf_w(i,j,i0)+rap*slhf_w(i,j,i2)
  137 continue

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 écheances
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PARTIE I
!     INITIALISATION DES FORCAGES EN DEBUT DE SIMULATION:
!     FIN.
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!                           /    /    /

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!            PARTIE II:
!     EVOLUTION DES FORCAGES ATMOSPHERIQUES AU COURS DE LA SIMULATION
! DEBUT:
      if (ichoix.eq.2) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      i0=0                                                             !11/07/02
      i2=2                                                             !11/07/02

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise à jour des échéances:
! DEBUT:
      do 1000 k1=1,6
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! commentees le                                                        !22/01/08
!     KOUNTMOD=AIRSEADT(K1,1) ! echeance
!     K2      =AIRSEADT(K1,2) ! repere NC=1
!     IF( MOD(KOUNT-K2,KOUNTMOD).EQ.0 ) THEN
! NC   enregistrement present NC+1 enregistrement pour dans 6 heures
!       NC=1+(KOUNT-K2)/KOUNTMOD
! Remplacees par: le:                                                  !22/01/08
      x2=(elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)
      x1=(elapsedtime_bef-airseadt(k1,2))/airseadt(k1,1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
      if(j2-j1.eq.1)  then                            !*******************>
! NC est le record au temps present ! NC+1 record pour echeance suivante:
      nc=int(1.+x2)

        lrec=(imax+2)*(jmax+2)*4
        open(unit=3,file=airseafile(k1),access='direct',recl=lrec)
        read(3,rec=nc+1)anyvar2d
        close(3)

! Commentees le                                                        !22/01/08
!       X1= AIRSEAINFO(K1,2)*REAL(KOUNTMOD)*DTI_FW
! Remplacees par: le:                                                  !22/01/08
!       x1= airseainfo(k1,2)*dti_bd*airseadt(k1,1)                      &
        x1= airseainfo(k1,2)*       airseadt(k1,1)             &    !17-04-11
       +(1.-airseainfo(k1,2))*1.


!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.1)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,1)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.2)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,2)=anyvar2d(i,j)/x1
              enddo
              enddo


! Rotation de grille:
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,1)=gridrotcos_t(i,j)*xy_t(i,j,1)-gridrotsin_t(i,j)*xy_t(i,j,2)
               xy_t(i,j,2)=gridrotsin_t(i,j)*xy_t(i,j,1)+gridrotcos_t(i,j)*xy_t(i,j,2)
              enddo
              enddo

! Regroupement aux points de vitesse:
              x0=windfactor*0.5
              do j=1,jmax+1
              do i=1,imax+1
               wstress_u(i,j,i0)=wstress_u(i,j,i2)
               wstress_v(i,j,i0)=wstress_v(i,j,i2)
               wstress_u(i,j,i2)=x0*(xy_t(i,j,1)+xy_t(i-1,j,1))
               wstress_v(i,j,i2)=x0*(xy_t(i,j,2)+xy_t(i,j-1,2))
              enddo
              enddo

            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.3)then
              do j=0,jmax+1
              do i=0,imax+1
               snsf_w(i,j,i0)=snsf_w(i,j,i2)
               snsf_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.4)then
              do j=0,jmax+1
              do i=0,imax+1
               ssr_w(i,j,i0)=ssr_w(i,j,i2)
               ssr_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.5)then
              do j=0,jmax+1
              do i=0,imax+1
               sshf_w(i,j,i0)=sshf_w(i,j,i2)
               sshf_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.6)then
              do j=0,jmax+1
              do i=0,imax+1
               slhf_w(i,j,i0)=slhf_w(i,j,i2)
               slhf_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !


      endif                                             !*******************>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise à jour des échéances:
! FIN.
 1000 continue
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



!                 ///         ///         ///



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 écheances
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      k1=1                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 122 j=1,jmax+1
      do 122 i=1,imax+1
      wstress_u(i,j,1)=(1.-rap)*wstress_u(i,j,i0)+rap*wstress_u(i,j,i2)
      wstress_v(i,j,1)=(1.-rap)*wstress_v(i,j,i0)+rap*wstress_v(i,j,i2)
  122 continue

! Module de la tension de vent:
      do 123 j=1,jmax
      do 123 i=1,imax
      wstress_w(i,j)=sqrt(                                              &
       ((wstress_u(i,j,1)+wstress_u(i+1,j,1))/2.)**2+                   &
       ((wstress_v(i,j,1)+wstress_v(i,j+1,1))/2.)**2)
  123 continue

      k1=3                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 124 j=1,jmax
      do 124 i=1,imax
      snsf_w(i,j,1)=(1.-rap)*snsf_w(i,j,i0)+rap*snsf_w(i,j,i2)
  124 continue

      k1=4                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 125 j=1,jmax
      do 125 i=1,imax
      ssr_w(i,j,1)=(1.-rap)*ssr_w(i,j,i0)+rap*ssr_w(i,j,i2)
  125 continue

      k1=5                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 126 j=1,jmax
      do 126 i=1,imax
      sshf_w(i,j,1)=(1.-rap)*sshf_w(i,j,i0)+rap*sshf_w(i,j,i2)
  126 continue

      k1=6                                                             !22/01/08
      rap=   (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)            & !22/01/08
       -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
      do 127 j=1,jmax
      do 127 i=1,imax
      slhf_w(i,j,1)=(1.-rap)*slhf_w(i,j,i0)+rap*slhf_w(i,j,i2)
  127 continue

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 écheances
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     OPEN(UNIT=3,FILE='mouchard1.dat',ACCESS='APPEND')
!     I=30
!     J=13
!     IF(MORPHO_Z(I,J,NR).EQ.1)WRITE(3,'(7(E11.5,1X))')
!    1  REAL(KOUNT)*DTI_FW/86400.
!    1 ,SNSF_Z(I,J,1)
!    1 ,SLHF_Z(I,J,1)
!    1 ,SSHF_Z(I,J,1)
!    1 , SSR_Z(I,J,1)
!    1 ,WSTRESS_X(I,J,1)
!    1 ,WSTRESS_Y(I,J,1)
!     CLOSE(3)

!     OPEN(UNIT=3,FILE='mouchard2.dat')
!     OPEN(UNIT=3,FILE='mouchard2.dat',ACCESS='APPEND')
!     I=20
!     J=13
!     K=2
!     DO I=1,MECO
!     WRITE(3,*)I,TKLL_Z(I,J,K),TKLE_Z(I,J,K)
!     ENDDO
!     IF(MORPHO_Z(I,J,NR).EQ.1)WRITE(3,'(7(E11.5,1X))')
!    1  REAL(KOUNT)*DTI_FW/86400.
!    1 ,VEL_X(I,J,NR-1,1)
!    1 ,VEL_Y(I,J,NR-1,1)
!    1 ,VEL_X(I,J,K,1)
!    1 ,VEL_Y(I,J,K,1)
!    1 ,TEM_Z(I,J,K)
!    1 ,SAL_Z(I,J,K)
!     CLOSE(3)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!            PARTIE II:
!     EVOLUTION DES FORCAGES ATMOSPHERIQUES AU COURS DE LA SIMULATION
!     FIN.
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      return
      end

