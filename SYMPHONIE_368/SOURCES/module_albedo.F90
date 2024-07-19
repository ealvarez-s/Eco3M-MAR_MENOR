      module module_albedo
!______________________________________________________________________
! S model
! release S26.1 - last update: 30-01-13
!______________________________________________________________________
      use module_principal
      implicit none
      double precision solardeclin
!--------------------------------------------------------------------------
! Version date Description des modifications
!         25/07/01: mise en service
!         01/12/05: On a enlev� la partie entiere devant kount dans X0
!                   qui nous a semble illogique puisque que daysim(0) est un reel
!                   portant un numero de jour decimal...
!         28/12/05: seuil sur l'absorption pour eviter division par zero
!                   dans temperature.F
!         20/01/06: Le calcul de l'absorption a �t� comment� car il me semblait
!                   bizarre....
!         18/04/06: fonctions compatibles avec double precision
!         04/05/06: fonctions compatibles avec double precision
! 2009.3  08-10-09: albedo devient albedo_z en 2D
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.12 20-09-10  Possibilit� de calcul en simple precision
! 2010.20 15-04-11  Calculs sur la base d'un temps exprim� en secondes
! 2010.25 27-06-12  Calculs remplac�s par une consante
! S26.1   30-01-13  Ajout de 2 autres methodes: br1982 et br1986
!--------------------------------------------------------------------------
contains

      subroutine albedo_upd
      implicit none
#ifdef synopsis
       subroutinetitle='albedo_upd'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! driver routine for albedo computation:
      call albedo_constant
!     call albedo_apel1987
!     call albedo_br1982
!     call albedo_br1986
      end subroutine albedo_upd

!......................................................................

      subroutine albedo_constant                 !27-06-12
      implicit none
#ifdef synopsis
       subroutinetitle='albedo_constant'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(iteration3d==0)albedo_w(:,:)=0.08
      end subroutine albedo_constant

!......................................................................

      subroutine albedo_br1982          !30-01-13
      implicit none
#ifdef synopsis
       subroutinetitle='albedo_br1982'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Briegleb & Ramanathan - Journal of Applied Meteorology - 1982

!  solardeclin est la d�clinaison du soleil : formule donn�e par pielke.
!  La formule peut etre verifiee en comparant a cette table:
!  http://www.wsanford.com/~wsanford/exo/sundials/DEC_Sun.html

!  Note DAYSIM(0)-1 est le numero du jour de reference (celui pour
!  lequel KOUNT=0).
!  Exemples: daysim(0)-1 du 1 janvier 0 h = 0.
!            daysim(0)-1 du 1 janvier 12h = 0.5
!            daysim(0)-1 du 2 janvier 12h = 1.5 ... etc...

      x0=2.*pi*(daysim(0)-1.+elapsedtime_now/86400.)/365. ! 2pi*(day number)/365. !15-04-11

      solardeclin=0.006918-0.39991*cos(x0)+0.070257*sin(x0)  &
                 -0.006758*cos(2*x0)+0.000907*sin(2*x0)      &
                 -0.002697*cos(3*x0)+0.00148*sin(3*x0)

      heure=mod(   elapsedtime_now/3600.                         &  !15-04-11
       +datesim(4,0)+datesim(5,0)/60.+datesim(6,0)/3600.,un*24.) & ! Heure � KOUNT=0
       -12.                                                      ! HEURE=0 � midi

!     do heure=-12.,12.,1.

      const1=180./pi/15.
      const2=pi/12.
      const5=180./pi
      do j=1,jmax
      do i=1,imax

      albedo_w(i,j)=0.05/(1.1*                                         &
                         (max(  sin(lat_t(i,j))*sin(solardeclin)       &
                               +cos(lat_t(i,j))*cos(solardeclin)       &
                               *cos( (heure+lon_t(i,j)*const1)*const2) &
                               , 0.0d0 ) )**1.4+0.15  )


!      if(i==imax.and.j==jmax) then
!       write(6,*)heure,albedo_w(i,j)
!      endif

      enddo
      enddo

!     enddo
!     stop 'bibi'
      end subroutine albedo_br1982

!......................................................................

      subroutine albedo_br1986          !30-01-13
      implicit none
      double precision mu_
#ifdef synopsis
       subroutinetitle='albedo_br1986'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Briegleb & Ramanathan - Journal of Applied Meteorology - 1986

!  solardeclin est la d�clinaison du soleil : formule donn�e par pielke.
!  La formule peut etre verifiee en comparant a cette table:
!  http://www.wsanford.com/~wsanford/exo/sundials/DEC_Sun.html

!  Note DAYSIM(0)-1 est le numero du jour de reference (celui pour
!  lequel KOUNT=0).
!  Exemples: daysim(0)-1 du 1 janvier 0 h = 0.
!            daysim(0)-1 du 1 janvier 12h = 0.5
!            daysim(0)-1 du 2 janvier 12h = 1.5 ... etc...

      x0=2.*pi*(daysim(0)-1.+elapsedtime_now/86400.)/365. ! 2pi*(day number)/365. !15-04-11

      solardeclin=0.006918-0.39991*cos(x0)+0.070257*sin(x0)  &
                 -0.006758*cos(2*x0)+0.000907*sin(2*x0)      &
                 -0.002697*cos(3*x0)+0.00148*sin(3*x0)

      heure=mod(   elapsedtime_now/3600.                         &  !15-04-11
       +datesim(4,0)+datesim(5,0)/60.+datesim(6,0)/3600.,un*24.) & ! Heure � KOUNT=0
       -12.                                                      ! HEURE=0 � midi

!     do heure=-12.,12.,1.

      const1=180./pi/15.
      const2=pi/12.
      const5=180./pi
      do j=1,jmax
      do i=1,imax

! mu_ is the cosine of the solar zenith angle:
              mu_=sin(lat_t(i,j))*sin(solardeclin)        &
                 +cos(lat_t(i,j))*cos(solardeclin)        &
         *cos( (heure+lon_t(i,j)*const1)*const2)

      albedo_w(i,j)=0.01*(  (2.6/( max(mu_,0.d00)**1.7d00+0.065d00  ))   &
                             +15.d00*(mu_-0.1d00)*(mu_-0.5d00)*(mu_-1.d00)    )

!      if(i==imax.and.j==jmax) then
!       write(6,*)heure,albedo_w(i,j)
!      endif

      enddo
      enddo

!     enddo
!     stop 'bibi'
      end subroutine albedo_br1986

!......................................................................

      subroutine albedo_apel1987
      implicit none
#ifdef synopsis
       subroutinetitle='albedo_apel1987'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!  solardeclin est la d�clinaison du soleil : formule donn�e par pielke.
!  La formule peut etre verifiee en comparant a cette table:
!  http://www.wsanford.com/~wsanford/exo/sundials/DEC_Sun.html

!  Note DAYSIM(0)-1 est le numero du jour de reference (celui pour
!  lequel KOUNT=0).
!  Exemples: daysim(0)-1 du 1 janvier 0 h = 0.
!            daysim(0)-1 du 1 janvier 12h = 0.5
!            daysim(0)-1 du 2 janvier 12h = 1.5 ... etc...

      x0=2.*pi*(daysim(0)-1.+elapsedtime_now/86400.)/365. ! 2pi*(day number)/365. !15-04-11

      solardeclin=0.006918-0.39991*cos(x0)+0.070257*sin(x0)  &
                 -0.006758*cos(2*x0)+0.000907*sin(2*x0)      &
                 -0.002697*cos(3*x0)+0.00148*sin(3*x0)

      heure=mod(   elapsedtime_now/3600.                         &  !15-04-11
       +datesim(4,0)+datesim(5,0)/60.+datesim(6,0)/3600.,un*24.) & ! Heure � KOUNT=0
       -12.                                                      ! HEURE=0 � midi

      const1=180./pi/15.
      const2=pi/12.
!     const3=0.001
      const4=90.
      const5=180./pi

!     do heure=-12.,12.,1.
      do j=1,jmax
      do i=1,imax

      albedo_w(i,j)=0.02                                               &
       +exp(4.e-6                                                      &
         *min(90.d0,                                                   &
              const5*acos(max(  sin(lat_t(i,j))*sin(solardeclin)       &
                               +cos(lat_t(i,j))*cos(solardeclin)       &
                               *cos( (heure+lon_t(i,j)*const1)*const2) &
                               ,0.001d0 ) ) )**3.1  )/100.

!      if(i==imax.and.j==jmax) then
!       write(6,*)heure,albedo_w(i,j)
!      endif
      enddo
      enddo

!     enddo
!     stop 'bibi'

! HEURE = heure locale par rapport � midi.
! convention pour le calcul de l'heure (HEURE):
! 1. il s'agit de l'heure locale: il faut donc ajouter le decalage
! horaire par rapport au temps T.U. (le temps du mod�le). Celui ci
! est donn� par la longitude (1heure pour 15 degr�s)
! 2. par convention HEURE est nul � midi.
!     heure=mod(real(kount)*dti_fw/3600.
!    1 +datesim(4,0)+datesim(5,0)/60.+datesim(6,0)/3600.,un*24.) ! Heure � KOUNT=0
!    2 -12.                                               ! HEURE=0 � midi
!    3 + lon_z(imax/2,jmax/2)*180./pi/15.      ! + decalage horaire local

! X5 est le cosinus de l'angle au zenith
!     x5=max(x1+x2*cos(heure*pi/12.),un*0.001)
!     WRITE(6,*)X5,'X5'

! X6 est l'angle au zenith
!     x6=min(un*90.,acos(x5)*180./pi)
!     WRITE(6,*)'zenith: ',X6
! formule Albedo tir�e d'Apel p.541 (figure 9-16)
!     albedo=0.02+exp(4.e-6*x6**3.1)/100.
! ci dessous 12 correspond au coefficient d absorption des eaux du cas II
! Figure 9-44 Apel
!                                                         comment�e le !20/01/06
!     ABSORPTION=MAX(12.*X5,0.1)                                       !28/12/05
!     WRITE(6,*)'ALBEDO',ALBEDO
!     WRITE(66,*)KOUNT*DTI_FW/86400.,ALBEDO,ABSORPTION

!c    IF(MOD(KOUNT,10).EQ.0) THEN
!c    OPEN(UNIT=3,FILE='helios.dat',ACCESS='APPEND')
!c    WRITE(3,'(I6,4(1X,E11.5))')KOUNT,HEURE,ALBEDO,X6,ABSORPTION
!c    WRITE(6,*)KOUNT,HEURE,ALBEDO,X6,ABSORPTION
!c    WRITE(6,*)'LATMOY:',SIN(PLAT_Z(MECO/2,NECO/2))
!c   1 ,PLAT_Z(MECO/2,NECO/2)*180./PI
!c    CLOSE(3)
!c    ENDIF
      end subroutine albedo_apel1987

      end module module_albedo
