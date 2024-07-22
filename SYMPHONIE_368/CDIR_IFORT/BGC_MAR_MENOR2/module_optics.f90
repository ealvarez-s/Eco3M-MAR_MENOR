










      module module_optics
!______________________________________________________________________
! SYMPHONIE ocean model
! release 272 - last update: 05-01-20
!______________________________________________________________________
      use module_principal ; use module_parallele
      implicit none
      double precision :: solardeclin,filvalr8_=-9999.
!...............................................................................
! Version date Description des modifications
!         25/07/01: mise en service
!         01/12/05: On a enlevE la partie entiere devant kount dans X0
!                   qui nous a semble illogique puisque que daysim(0) est un reel
!                   portant un numero de jour decimal...
!         28/12/05: seuil sur l'absorption pour eviter division par zero
!                   dans temperature.F
!         20/01/06: Le calcul de l'absorption a EtE commentE car il me semblait
!                   bizarre....
!         18/04/06: fonctions compatibles avec double precision
!         04/05/06: fonctions compatibles avec double precision
! 2009.3  08-10-09: albedo devient optics_albedo_z en 2D
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.12 20-09-10  PossibilitE de calcul en simple precision
! 2010.20 15-04-11  Calculs sur la base d'un temps exprimE en secondes
! 2010.25 27-06-12  Calculs remplacEs par une consante
! S26.1   30-01-13  Ajout de 2 autres methodes: br1982 et br1986
!         24-04-16  Ajout cas du modele 1DV
! v272    05-01-20  borne sup 30m
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
contains

      subroutine optics_notebook_optical
      implicit none
! notebook_optical
      namelist/notebook_optical/light_att1,light_att2,light_rat1,light_rat2,texte250 &
       ,light_att2_val1,light_att2_h1,light_att2_val2,light_att2_h2 &
       ,albedo_val,texte30


!...................................................
! Initial optical properties - Reads notebook_optical
! References:
! Paulson C. A., and J. J. Simpson, 1977: Irradiance !measurements in the upper ocean. J. Phys. Oceanogr., 7, !952-956.
! Maraldi et al, Ocean Science 2013


! Default values:
      light_att1=0.35       ! par1 (m)
      light_rat1=0.58       ! A1

      light_att2=23.        ! par2 (m) 
      light_rat2=0.42

      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(15)

      texte250='none' ; light_att1=filvalr8_ ; light_att2=filvalr8_
      light_att2_val1=filvalr8_ ; light_att2_h1=filvalr8_
      light_att2_val2=filvalr8_ ; light_att2_h2=filvalr8_ ; albedo_case=-9999
       open(100,file=nomfichier(15)) ! 'notebook_optical !24-05-15
       read(100,nml=notebook_optical)
      close(100)

! Select "albedo" case
      if(texte30=='constant')albedo_case=albedo_constant
      if(texte30=='apel1987')albedo_case=albedo_apel1987
      if(texte30=='br1982')  albedo_case=albedo_br1982
      if(texte30=='br1986')  albedo_case=albedo_br1986
      if(albedo_case==-9999)stop 'Err 69 albedo_case=-9999'

! Select "PAR2" case
      if(light_att1==filvalr8_)stop 'Err 62 light_att1==filvalr8_'
      light_kpar1  =1./light_att1
      if(texte250=='none') then !>>>>>
    
       if(light_att2_val1/=filvalr8_.and.light_att2_h1/=filvalr8_.and. &
          light_att2_val2/=filvalr8_.and.light_att2_h2/=filvalr8_) then !oooooo>

        call optics_initial_kpar2_hparam

       else                                                             !oooooo>

        if(light_att2==filvalr8_)stop 'Err 63 light_att2==filvalr8_'
        light_kpar2_w(:,:)=1./light_att2

       endif                                                            !oooooo>

      else                      !>>>>>

       call optics_initial_kpar2_file

      endif                     !>>>>>

! Ensure that light_rat1+light_rat2=1
      x1=(1.-light_rat1-light_rat2)/(light_rat1+light_rat2)
      light_rat1=light_rat1*(1.+x1)
      light_rat2=light_rat2*(1.+x1)

      end subroutine optics_notebook_optical

!..........................................................................

      subroutine optics_initial_kpar2_hparam
      implicit none
      integer ncid_
      include 'netcdf.inc'

      if(light_att2_h1>=light_att2_h2) then
        stop 'Err light_att2_h1>light_att2_h2'
      else
        x0=1./(light_att2_h2-light_att2_h1)
      endif

      do j=1,jmax ; do i=1,imax
        x1=min(max(x0*(h_w(i,j)-light_att2_h1),0.d0),1.d0)
        light_kpar2_w(i,j)=x1/light_att2_val2+(1.-x1)/light_att2_val1
      enddo ; enddo

      end subroutine optics_initial_kpar2_hparam

!..........................................................................

      subroutine optics_initial_kpar2_file
      implicit none
      integer ncid_
      include 'netcdf.inc'

! A l'origine, ce fichier etait produit par xscan sym-tools :

      if(par%rank==0)write(6,'(a,a)')'nf_open ',trim(texte250)

      status=nf_open(trim(texte250),nf_nowrite,ncid_)
      if(status/=0)stop 'Err 179 module_optics nf_open'

! check dimensions:
                    status=nf_inq_dimid(ncid_,'nx_t',k0)          
       if(status/=0)status=nf_inq_dimid(ncid_,'ni_t',k0) !05-01-20
       if(status/=0)stop 'par Err1890'
       status=nf_inq_dimlen(ncid_,k0,i1)             ;if(status/=0)stop 'par Err1891'
                    status=nf_inq_dimid(ncid_,'ny_t',k0)          
       if(status/=0)status=nf_inq_dimid(ncid_,'nj_t',k0) !05-01-20
       if(status/=0)stop 'par Err1892'
       status=nf_inq_dimlen(ncid_,k0,j1)             ;if(status/=0)stop 'par Err1893'
       if(i1/=iglb+2)stop 'Err dim1 netcdf par File'
       if(j1/=jglb+2)stop 'Err dim2 netcdf par File'
       varstart(1)=1+par%timax(1) ; varcount(1)=imax+2
       varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2


       status=nf_inq_varid(ncid_,'par_diffused_t',var_id)       ; if(status/=0)stop 'Err nf_inq_varid par'
       status=nf_inq_vartype(ncid_,var_id,k0)                   ; if(k0/=nf_real)stop 'Err uncorrect vartype par'
       status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                           ,varcount(1:2)   &
                             ,anyvar2d(0:imax+1,0:jmax+1))      ;  if(status/=0)stop 'Err 198 nf_get_vara_real par'


      status=nf_close(ncid_)
      if(status/=0)stop 'Err 2003 nf_close module_optics'

      do j=1,jmax ; do i=1,imax
!       light_kpar2_w(i,j)=1./min(max(anyvar2d(i,j),0.35),23.)
        light_kpar2_w(i,j)=1./min(max(anyvar2d(i,j),0.35),30.) !05-01-20
      enddo ; enddo

      end subroutine optics_initial_kpar2_file
!..........................................................................

      subroutine optics_albedo_upd
      implicit none

! driver routine for albedo computation:
      if(albedo_case==albedo_constant)then ; call optics_albedo_constant ; return ; endif
      if(albedo_case==albedo_apel1987)then ; call optics_albedo_apel1987 ; return ; endif
      if(albedo_case==albedo_br1982  )then ; call optics_albedo_br1982   ; return ; endif
      if(albedo_case==albedo_br1986  )then ; call optics_albedo_br1986   ; return ; endif

      stop 'Err 189 albedo_case not defined'
      end subroutine optics_albedo_upd

!......................................................................

      subroutine optics_albedo_constant                 !27-06-12
      implicit none

      if(iteration3d==iteration3d_restart)albedo_w(:,:)=albedo_val

      end subroutine optics_albedo_constant

!......................................................................

      subroutine optics_albedo_br1982          !30-01-13
      implicit none

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
      if(flag_1dv==1)albedo_w(:,:)=albedo_w(2,2) !24-04-16

!     enddo
!     stop 'bibi'
      end subroutine optics_albedo_br1982

!......................................................................

      subroutine optics_albedo_br1986          !30-01-13
      implicit none
      double precision mu_

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
      if(flag_1dv==1)albedo_w(:,:)=albedo_w(2,2) !24-04-16

!     enddo
!     stop 'bibi'
      end subroutine optics_albedo_br1986

!......................................................................

      subroutine optics_albedo_apel1987
      implicit none

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

      if(flag_1dv==1)albedo_w(:,:)=albedo_w(2,2) !24-04-16

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
      end subroutine optics_albedo_apel1987

!......................................................................

      end module module_optics
