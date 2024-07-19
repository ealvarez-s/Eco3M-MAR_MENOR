      subroutine usersbiobc
      use module_principal
      use module_parallele
      use ModuleDeclaration
      implicit none
      double precision InitValue(20,3000),Prof(20,3000)
      integer ichoix,nmax,K_pr,nmax1,nmax2,nmax3
      double precision salz,salmu,salsigma,p1,p2,p3,p4,p5

! coeff pour la relation oxygène/salinité polynôme 4
      salmu = 38.49 ! moyenne des salinité des mesures ctd
      salsigma = 0.085254 ! ecart type des salinité
      p1 = -8.9604d-5
      p2 = -0.016762
      p3 = -0.9261
      p4 = -15.234
      p5 = 203.9


! Modele 3D entouré de modèles 1D : on ne fait rien                       !24/03/08
!     if(biobc_type(vb)==4) CONTINUE

! Nutrients deduced from density for boundary points and nul gradient
       if(biobc_type(vb)==5) then
!!print*,'Dans obc_bio, biobc_type=6)'
!!       if(eqs_state1.eq.0)call density(0) ! linéaire
!!       if(eqs_state1.eq.1)call density(3) !Non line sans pr
         call equation_of_state('potential density',1)

! Southern and northern boundaries
       do k1=1,2

       if(k1==1) then
           if(par%tvoisin(sud  )/=mpi_proc_null)goto 595 ! si interieur domaine pas de C.L.
           j=1 ; j1=0
       endif
       if(k1==2) then
           if(par%tvoisin(nord )/=mpi_proc_null)goto 595 ! si interieur domaine pas de C.L.
           j=jmax ; j1=jmax+1
       endif

       do i=1,imax
       do k=kmin_w(i,j),kmax ! ajout caroline
       
        if (mask_t(i,j,k)==1) then
        if (mask_t(i,j1,k)==1) then

!! deduction of nutrients
!!       bio_t(i,j,k,iNITRATE)=26.5625*(rhp_t(i,j,k)+rho-1000.)-764.47
!!       if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!!       if(bio_t(i,j,k,iNITRATE)<=0.2) bio_t(i,j,k,iNITRATE)=0.2
!!
!!       bio_t(i,j,k,iPHOSPHATE)=1.6667*(rhp_t(i,j,k)+rho-1000.)-48.1
!!       if(bio_t(i,j,k,iPHOSPHATE)>=0.4 ) bio_t(i,j,k,iPHOSPHATE)=0.4
!!       if(bio_t(i,j,k,iPHOSPHATE)<=0.02) bio_t(i,j,k,iPHOSPHATE)=0.02
!!
!!       bio_t(i,j,k,iSILICE)=bio_t(i,j,k,iNITRATE)/1.4


!         bio_t(i,j,k,iNITRATE)=12.475*(rhp_t(i,j,k)+rho-1000.)-355.17
!         if(bio_t(i,j,k,iNITRATE)>=9) bio_t(i,j,k,iNITRATE)=9
!         if(bio_t(i,j,k,iNITRATE)<=0.032) bio_t(i,j,k,iNITRATE)=0.032
!
!         bio_t(i,j,k,iPHOSPHATE)=0.63649*(rhp_t(i,j,k)+rho-1000.)-18.176
!         if(bio_t(i,j,k,iPHOSPHATE)>=0.44) bio_t(i,j,k,iPHOSPHATE)=0.44
!         if(bio_t(i,j,k,iPHOSPHATE)<=0.02) bio_t(i,j,k,iPHOSPHATE)=0.02
!
!         bio_t(i,j,k,iSILICE)=12.037*(rhp_t(i,j,k)+rho-1000.)-343.27
!         if(bio_t(i,j,k,iSILICE)>=9) bio_t(i,j,k,iSILICE)=9
!         if(bio_t(i,j,k,iSILICE)<=0.03) bio_t(i,j,k,iSILICE)=0.03  !30/11/12





         if (                              &
         lon_t(i,j)>=9.5*deg2rad.and.                            & 
         lon_t(i,j)<13.*deg2rad.and.                            &
         lat_t(i,j)>41.5*deg2rad)then  !.and.lat_t(i,j)<44.*deg2rad)then


         if(vb==iNITRATE) then
         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513
         if(bio_t(i,j,k,iNITRATE)>=7.5) bio_t(i,j,k,iNITRATE)=7.5
         if(bio_t(i,j,k,iNITRATE)<=0.35) bio_t(i,j,k,iNITRATE)=0.35

         elseif(vb==iPHOSPHATE) then
         write(20+par%rank,*),i,j,k,rho,rhp_t(i,j,k)
         print*,'obc_bio',i,j,k,rho,rhp_t(i,j,k)
         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
         if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
         if(bio_t(i,j,k,iPHOSPHATE)<=0.0317) bio_t(i,j,k,iPHOSPHATE)=0.0317

         elseif(vb==iSILICE) then
         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
         if(bio_t(i,j,k,iSILICE)>=7.) bio_t(i,j,k,iSILICE)=7.
         if(bio_t(i,j,k,iSILICE)<=0.30) bio_t(i,j,k,iSILICE)=0.30

         endif


        elseif (                               &
        lon_t(i,j)>=9.5*deg2rad.and.                                 &
        lon_t(i,j)<16.*deg2rad.and.                                 &
        lat_t(i,j)<=41.5*deg2rad.and.lat_t(i,j)>38.5*deg2rad)then

         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513
         if(bio_t(i,j,k,iNITRATE)>=7.5) bio_t(i,j,k,iNITRATE)=7.5
         if(bio_t(i,j,k,iNITRATE)<=0.35) bio_t(i,j,k,iNITRATE)=0.35

         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
         if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
         if(bio_t(i,j,k,iPHOSPHATE)<=0.0317) bio_t(i,j,k,iPHOSPHATE)=0.0317

         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
         if(bio_t(i,j,k,iSILICE)>=7.) bio_t(i,j,k,iSILICE)=7.
         if(bio_t(i,j,k,iSILICE)<=0.30) bio_t(i,j,k,iSILICE)=0.30




! bioregion 12  subbassin algero-provençal

         elseif (                               &
         lon_t(i,j)>-2.*deg2rad.and.                                 &
         lon_t(i,j)<9.5*deg2rad.and.                                 &
         lat_t(i,j)<38.*deg2rad)then



         bio_t(i,j,k,iNITRATE)=8.38357*(rhp_t(i,j,k)+rho-1000.)-234.988
         if(bio_t(i,j,k,iNITRATE)>=10.) bio_t(i,j,k,iNITRATE)=10.
!         if(bio_t(i,j,k,iNITRATE)<=0.7015) bio_t(i,j,k,iNITRATE)=0.7015
         if(bio_t(i,j,k,iNITRATE)<=0.35) bio_t(i,j,k,iNITRATE)=0.35

         bio_t(i,j,k,iPHOSPHATE)=0.366129*(rhp_t(i,j,k)+rho-1000.)-10.261
         if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
         if(bio_t(i,j,k,iPHOSPHATE)<=0.0337) bio_t(i,j,k,iPHOSPHATE)=0.0337

         bio_t(i,j,k,iSILICE)=9.11274*(rhp_t(i,j,k)+rho-1000.)-256.249
         if(bio_t(i,j,k,iSILICE)>=9.5) bio_t(i,j,k,iSILICE)=9.5
         if(bio_t(i,j,k,iSILICE)<=0.9832) bio_t(i,j,k,iSILICE)=0.9832



! sicile
        elseif(                                  &
        lon_t(i,j)>=9.5*deg2rad.and.lon_t(i,j)<11.*deg2rad.and.         &
        lat_t(i,j)>36.5*deg2rad.and.lat_t(i,j)<=38.5*deg2rad)then

        bio_t(i,j,k,iNITRATE)=12.3021*(rhp_t(i,j,k)+rho-1000.)-349.366
        if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
        if(bio_t(i,j,k,iNITRATE)<=0.5475) bio_t(i,j,k,iNITRATE)=0.5475
 
        bio_t(i,j,k,iPHOSPHATE)=0.337703*(rhp_t(i,j,k)+rho-1000.)-9.56411
        if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0573) bio_t(i,j,k,iPHOSPHATE)=0.0573

        bio_t(i,j,k,iSILICE)=11.4883*(rhp_t(i,j,k)+rho-1000.)-326.085
        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
        if(bio_t(i,j,k,iSILICE)<=0.6763) bio_t(i,j,k,iSILICE)=0.6763
! malte
        elseif (                                &
        lon_t(i,j)>10.*deg2rad.and.lon_t(i,j)<=15.*deg2rad.and.        &
        lat_t(i,j)>34.*deg2rad.and.lat_t(i,j)<=36.5*deg2rad)then
       
        bio_t(i,j,k,iNITRATE)=4.64766*(rhp_t(i,j,k)+rho-1000.)-130.081
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.1432) bio_t(i,j,k,iNITRATE)=0.1432

        bio_t(i,j,k,iPHOSPHATE)=0.149974*(rhp_t(i,j,k)+rho-1000.)-4.1628
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0519) bio_t(i,j,k,iPHOSPHATE)=0.0519

        bio_t(i,j,k,iSILICE)=4.88587*(rhp_t(i,j,k)+rho-1000.)-136.875
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.1707) bio_t(i,j,k,iSILICE)=0.1707
! libye
        elseif(                                  &
        lon_t(i,j)>15.*deg2rad.and.lon_t(i,j)<=23.*deg2rad.and.        &
        lat_t(i,j)<=36.5*deg2rad)then
       
        bio_t(i,j,k,iNITRATE)=10.2623*(rhp_t(i,j,k)+rho-1000.)-239.769
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.1967) bio_t(i,j,k,iNITRATE)=0.1967

        bio_t(i,j,k,iPHOSPHATE)=0.345349*(rhp_t(i,j,k)+rho-1000.)-9.87983
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0328) bio_t(i,j,k,iPHOSPHATE)=0.0328

        bio_t(i,j,k,iSILICE)=10.548*(rhp_t(i,j,k)+rho-1000.)-301.996
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.3305) bio_t(i,j,k,iSILICE)=0.3305

! tun
        elseif(                                  &
        lon_t(i,j)>10.*deg2rad.and.lon_t(i,j)<=15.*deg2rad.and.        &
!        lat_t(i,j)>36.*deg2rad.and.lat_t(i,j)<38.*deg2rad)then
        lat_t(i,j)<=34.*deg2rad)then
       
        bio_t(i,j,k,iNITRATE)=4.8857*(rhp_t(i,j,k)+rho-1000.)-136.453
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.2118) bio_t(i,j,k,iNITRATE)=0.2118

        bio_t(i,j,k,iPHOSPHATE)=0.148005*(rhp_t(i,j,k)+rho-1000.)-4.11025
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0341) bio_t(i,j,k,iPHOSPHATE)=0.0341

        bio_t(i,j,k,iSILICE)=4.752*(rhp_t(i,j,k)+rho-1000.)-132.725
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.2682) bio_t(i,j,k,iSILICE)=0.2682


! ionnienne
        elseif(                                  &
        lon_t(i,j)>15.*deg2rad.and.lon_t(i,j)<22.*deg2rad.and.        &
        lat_t(i,j)>36.5*deg2rad.and.lat_t(i,j)<39.*deg2rad)then
!        lat_t(i,j)<39.*deg2rad)then

        bio_t(i,j,k,iNITRATE)=11.9114*(rhp_t(i,j,k)+rho-1000.)-341.931
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.096) bio_t(i,j,k,iNITRATE)=0.096

        bio_t(i,j,k,iPHOSPHATE)=0.332844*(rhp_t(i,j,k)+rho-1000.)-9.52928
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0348) bio_t(i,j,k,iPHOSPHATE)=0.0348

        bio_t(i,j,k,iSILICE)=11.6992*(rhp_t(i,j,k)+rho-1000.)-335.683
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.2943) bio_t(i,j,k,iSILICE)=0.2943

! cypre
        elseif(lon_t(i,j)>23.*deg2rad.and.        &
!        lat_t(i,j)<36.*deg2rad)then
         lat_t(i,j)<37.5*deg2rad)then
      
        bio_t(i,j,k,iNITRATE)=14.4713*(rhp_t(i,j,k)+rho-1000.)-417.015
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.1656) bio_t(i,j,k,iNITRATE)=0.1656

        bio_t(i,j,k,iPHOSPHATE)=0.55*(rhp_t(i,j,k)+rho-1000.)-16.09
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0232) bio_t(i,j,k,iPHOSPHATE)=0.0232

        bio_t(i,j,k,iSILICE)=14.0977*(rhp_t(i,j,k)+rho-1000.)-406.195
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.2344) bio_t(i,j,k,iSILICE)=0.2344

! NAD
        elseif(                                  &
        lon_t(i,j)>12.*deg2rad.and.lon_t(i,j)<20.*deg2rad.and.        &
        lat_t(i,j)>40.*deg2rad )then

        bio_t(i,j,k,iNITRATE)=16.9461*(rhp_t(i,j,k)+rho-1000.)-488.995
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.0623) bio_t(i,j,k,iNITRATE)=0.0623

        bio_t(i,j,k,iPHOSPHATE)=0.814105*(rhp_t(i,j,k)+rho-1000.)-23.5077
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0371) bio_t(i,j,k,iPHOSPHATE)=0.0371

        bio_t(i,j,k,iSILICE)=15.555*(rhp_t(i,j,k)+rho-1000.)-448.71
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.1188) bio_t(i,j,k,iSILICE)=0.1188

! SAD
        elseif(                                  &
        lon_t(i,j)>16.*deg2rad.and.lon_t(i,j)<20.*deg2rad.and.        &
        lat_t(i,j)>39.*deg2rad.and.lat_t(i,j)<40.*deg2rad)then

        bio_t(i,j,k,iNITRATE)=16.0433*(rhp_t(i,j,k)+rho-1000.)-462.552
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.0942) bio_t(i,j,k,iNITRATE)=0.0942

        bio_t(i,j,k,iPHOSPHATE)=0.5645*(rhp_t(i,j,k)+rho-1000.)-16.2724
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0422) bio_t(i,j,k,iPHOSPHATE)=0.0422

        bio_t(i,j,k,iSILICE)=16.0979*(rhp_t(i,j,k)+rho-1000.)-464.141
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.122) bio_t(i,j,k,iSILICE)=0.122


! alboran
        elseif(                                  &
        lon_t(i,j)<=-2.*deg2rad.and.lon_t(i,j)>-5.*deg2rad)then

!        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!        if(bio_t(i,j,k,iNITRATE)<=0.2128) bio_t(i,j,k,iNITRATE)=0.2128
!
!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.0753) bio_t(i,j,k,iPHOSPHATE)=0.0753
!
!        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
!        if(bio_t(i,j,k,iSILICE)<=0.2431) bio_t(i,j,k,iSILICE)=0.2431

! new 28032013 
        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924                 ! original
!         bio_t(i,j,k,iNITRATE)=8.38357*(rhp_t(i,j,k)+rho-1000.)-234.988     ! alg
        if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
!        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08981*(rhp_t(i,j,k)+rho-1000.)-30.6768)   ! alg
        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(2.00125*(rhp_t(i,j,k)+rho-1000.)-53.8274)  ! original

        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102    ! original
!         bio_t(i,j,k,iPHOSPHATE)=0.366129*(rhp_t(i,j,k)+rho-1000.)-10.261    ! alg
        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.890801*(rhp_t(i,j,k)+rho-1000.)-27.9247)  ! alg
        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(1.86852*(rhp_t(i,j,k)+rho-1000.)-52.7586)     ! original

        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489                        ! original
!         bio_t(i,j,k,iSILICE)=9.11274*(rhp_t(i,j,k)+rho-1000.)-256.249   ! alg
       if(bio_t(i,j,k,iSILICE)>=7.5) bio_t(i,j,k,iSILICE)=7.5
        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.78793*(rhp_t(i,j,k)+rho-1000.)-48.4092)    !original
!        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.05871*(rhp_t(i,j,k)+rho-1000.)-29.887)  ! alg





! ege
        elseif(                                   &
        lon_t(i,j)>23.*deg2rad.and.                                   &
        lat_t(i,j)>=37.5*deg2rad)then

        bio_t(i,j,k,iNITRATE)=9.27188*(rhp_t(i,j,k)+rho-1000.)-267.101
        if(bio_t(i,j,k,iNITRATE)>=4.) bio_t(i,j,k,iNITRATE)=4.
        if(bio_t(i,j,k,iNITRATE)<=0.2032) bio_t(i,j,k,iNITRATE)=0.2032

        bio_t(i,j,k,iPHOSPHATE)=0.366129*(rhp_t(i,j,k)+rho-1000.)-9.21097
        if(bio_t(i,j,k,iPHOSPHATE)>=0.2) bio_t(i,j,k,iPHOSPHATE)=0.2
        if(bio_t(i,j,k,iPHOSPHATE)<=0.043) bio_t(i,j,k,iPHOSPHATE)=0.043

        bio_t(i,j,k,iSILICE)=8.92964*(rhp_t(i,j,k)+rho-1000.)-257.004
        if(bio_t(i,j,k,iSILICE)>=4.) bio_t(i,j,k,iSILICE)=4.
        if(bio_t(i,j,k,iSILICE)<=0.3416) bio_t(i,j,k,iSILICE)=0.3416



!! gibraltar
        elseif(                                  &
        lon_t(i,j)<=-5.*deg2rad)then



!        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08981*(rhp_t(i,j,k)+rho-1000.)-30.6768)   ! alg
!        
!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.890801*(rhp_t(i,j,k)+rho-1000.)-27.9247)  ! alg
!
!        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
!        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.05871*(rhp_t(i,j,k)+rho-1000.)-29.887)  ! alg
!
!
        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
       if(bio_t(i,j,k,iNITRATE)>=10.) bio_t(i,j,k,iNITRATE)=10.
!       if(bio_t(i,j,k,iNITRATE)<=0.032) bio_t(i,j,k,iNITRATE)=0.032
      if(bio_t(i,j,k,iNITRATE)<=0.14) bio_t(i,j,k,iNITRATE)=0.14


        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
       if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
!       if(bio_t(i,j,k,iPHOSPHATE)<=0.02) bio_t(i,j,k,iPHOSPHATE)=0.02
       if(bio_t(i,j,k,iPHOSPHATE)<=0.03) bio_t(i,j,k,iPHOSPHATE)=0.03



        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
       if(bio_t(i,j,k,iSILICE)>=9.) bio_t(i,j,k,iSILICE)=9.
!       if(bio_t(i,j,k,iSILICE)<=0.03) bio_t(i,j,k,iSILICE)=0.03
       if(bio_t(i,j,k,iSILICE)<=0.12) bio_t(i,j,k,iSILICE)=0.12



!
!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!!!inflow
!        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iNITRATE)=0.4
!        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iPHOSPHATE)=0.12
!        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iSILICE)=0.79
!!!outflow
!        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iNITRATE)=8.
!        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iSILICE)=6.47
!!!!!!!!!!!!

!!        if(bio_t(i,j,k,iNITRATE)<=0.2128) bio_t(i,j,k,iNITRATE)=0.2128
!!
!!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!!        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!!        if(bio_t(i,j,k,iPHOSPHATE)<=0.0753) bio_t(i,j,k,iPHOSPHATE)=0.0753
!!       bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!!        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
!!        if(bio_t(i,j,k,iSILICE)<=0.2431) bio_t(i,j,k,iSILICE)=0.2431





! bioregion nord occidental (provençale + ligure)

       else
         
!       bio_t(i,j,k,iNITRATE)=12.475*(rhp_t(i,j,k)+rho-1000.)-355.17
       bio_t(i,j,k,iNITRATE)=16.0493*(rhp_t(i,j,k)+rho-1000.)-458.667
       if(bio_t(i,j,k,iNITRATE)>=9) bio_t(i,j,k,iNITRATE)=9
!       if(bio_t(i,j,k,iNITRATE)<=0.032) bio_t(i,j,k,iNITRATE)=0.032
      if(bio_t(i,j,k,iNITRATE)<=0.14) bio_t(i,j,k,iNITRATE)=0.14

!       bio_t(i,j,k,iPHOSPHATE)=0.63649*(rhp_t(i,j,k)+rho-1000.)-18.176
       bio_t(i,j,k,iPHOSPHATE)=0.842999*(rhp_t(i,j,k)+rho-1000.)-24.1636

       if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
!       if(bio_t(i,j,k,iPHOSPHATE)<=0.02) bio_t(i,j,k,iPHOSPHATE)=0.02
       if(bio_t(i,j,k,iPHOSPHATE)<=0.03) bio_t(i,j,k,iPHOSPHATE)=0.03

!       bio_t(i,j,k,iSILICE)=12.037*(rhp_t(i,j,k)+rho-1000.)-343.27
       bio_t(i,j,k,iSILICE)=15.7532*(rhp_t(i,j,k)+rho-1000.)-450.96

       if(bio_t(i,j,k,iSILICE)>=8) bio_t(i,j,k,iSILICE)=8
!       if(bio_t(i,j,k,iSILICE)<=0.03) bio_t(i,j,k,iSILICE)=0.03
       if(bio_t(i,j,k,iSILICE)<=0.12) bio_t(i,j,k,iSILICE)=0.12





       endif

         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5

         if(bio_t(i,j,k,iOXYGEN)>=260) bio_t(i,j,k,iOXYGEN)=260
         if(bio_t(i,j,k,iOXYGEN)<=170) bio_t(i,j,k,iOXYGEN)=170


! Gradient nul "z1"
         bio_t(i,j1,k,vb)=bio_t(i,j,k,vb)

         endif ! mask
         endif ! mask
         enddo ! k ajout caroline
         enddo ! i

  595    continue
         enddo ! k1



! Western and eastern boundaries
       do k1=1,2

       if(k1==1) then
          if(par%tvoisin(ouest)/=mpi_proc_null)goto 987 ! si interieur domaine pas de C.L.
          i=1 ; i1=0
       endif
       if(k1==2) then
          if(par%tvoisin(est  )/=mpi_proc_null)goto 987 ! si interieur domaine pas de C.L.
          i=imax ; i1=imax+1
       endif

       do j=1,jmax
       do k=kmin_w(i,j),kmax ! ajout caroline


       if (mask_t(i,j,k)==1) then
       if (mask_t(i1,j,k)==1) then

         if (                              &
         lon_t(i,j)>=9.5*deg2rad.and.                            & 
         lon_t(i,j)<13.*deg2rad.and.                            &
         lat_t(i,j)>41.5*deg2rad)then  !.and.lat_t(i,j)<44.*deg2rad)then



         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513
         if(bio_t(i,j,k,iNITRATE)>=7.5) bio_t(i,j,k,iNITRATE)=7.5
         if(bio_t(i,j,k,iNITRATE)<=0.35) bio_t(i,j,k,iNITRATE)=0.35

         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
         if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
         if(bio_t(i,j,k,iPHOSPHATE)<=0.0317) bio_t(i,j,k,iPHOSPHATE)=0.0317

         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
         if(bio_t(i,j,k,iSILICE)>=7.) bio_t(i,j,k,iSILICE)=7.
         if(bio_t(i,j,k,iSILICE)<=0.30) bio_t(i,j,k,iSILICE)=0.30



        elseif (                               &
        lon_t(i,j)>=9.5*deg2rad.and.                                 &
        lon_t(i,j)<16.*deg2rad.and.                                 &
        lat_t(i,j)<=41.5*deg2rad.and.lat_t(i,j)>38.5*deg2rad)then

         bio_t(i,j,k,iNITRATE)=15.278*(rhp_t(i,j,k)+rho-1000.)-435.513
         if(bio_t(i,j,k,iNITRATE)>=7.5) bio_t(i,j,k,iNITRATE)=7.5
         if(bio_t(i,j,k,iNITRATE)<=0.35) bio_t(i,j,k,iNITRATE)=0.35

         bio_t(i,j,k,iPHOSPHATE)=0.371067*(rhp_t(i,j,k)+rho-1000.)-10.4974
         if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
         if(bio_t(i,j,k,iPHOSPHATE)<=0.0317) bio_t(i,j,k,iPHOSPHATE)=0.0317

         bio_t(i,j,k,iSILICE)=12.8943*(rhp_t(i,j,k)+rho-1000.)-367.212
         if(bio_t(i,j,k,iSILICE)>=7.) bio_t(i,j,k,iSILICE)=7.
         if(bio_t(i,j,k,iSILICE)<=0.30) bio_t(i,j,k,iSILICE)=0.30



! bioregion 12  subbassin algero-provençal

         elseif (                               &
         lon_t(i,j)>-2.*deg2rad.and.                                 &
         lon_t(i,j)<9.5*deg2rad.and.                                 &
         lat_t(i,j)<38.*deg2rad)then



         bio_t(i,j,k,iNITRATE)=8.38357*(rhp_t(i,j,k)+rho-1000.)-234.988
         if(bio_t(i,j,k,iNITRATE)>=10.) bio_t(i,j,k,iNITRATE)=10.
!         if(bio_t(i,j,k,iNITRATE)<=0.7015) bio_t(i,j,k,iNITRATE)=0.7015
         if(bio_t(i,j,k,iNITRATE)<=0.35) bio_t(i,j,k,iNITRATE)=0.35

         bio_t(i,j,k,iPHOSPHATE)=0.366129*(rhp_t(i,j,k)+rho-1000.)-10.261
         if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
         if(bio_t(i,j,k,iPHOSPHATE)<=0.0337) bio_t(i,j,k,iPHOSPHATE)=0.0337

         bio_t(i,j,k,iSILICE)=9.11274*(rhp_t(i,j,k)+rho-1000.)-256.249
         if(bio_t(i,j,k,iSILICE)>=9.5) bio_t(i,j,k,iSILICE)=9.5
         if(bio_t(i,j,k,iSILICE)<=0.9832) bio_t(i,j,k,iSILICE)=0.9832



! sicile
        elseif(                                  &
        lon_t(i,j)>=9.5*deg2rad.and.lon_t(i,j)<11.*deg2rad.and.         &
        lat_t(i,j)>36.5*deg2rad.and.lat_t(i,j)<=38.5*deg2rad)then

        bio_t(i,j,k,iNITRATE)=12.3021*(rhp_t(i,j,k)+rho-1000.)-349.366
        if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
        if(bio_t(i,j,k,iNITRATE)<=0.5475) bio_t(i,j,k,iNITRATE)=0.5475
 
        bio_t(i,j,k,iPHOSPHATE)=0.337703*(rhp_t(i,j,k)+rho-1000.)-9.56411
        if(bio_t(i,j,k,iPHOSPHATE)>=0.3) bio_t(i,j,k,iPHOSPHATE)=0.3
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0573) bio_t(i,j,k,iPHOSPHATE)=0.0573

        bio_t(i,j,k,iSILICE)=11.4883*(rhp_t(i,j,k)+rho-1000.)-326.085
        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
        if(bio_t(i,j,k,iSILICE)<=0.6763) bio_t(i,j,k,iSILICE)=0.6763
! malte
        elseif (                                 &
        lon_t(i,j)>10.*deg2rad.and.lon_t(i,j)<=15.*deg2rad.and.        &
        lat_t(i,j)>34.*deg2rad.and.lat_t(i,j)<=36.5*deg2rad)then
       
        bio_t(i,j,k,iNITRATE)=4.64766*(rhp_t(i,j,k)+rho-1000.)-130.081
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.1432) bio_t(i,j,k,iNITRATE)=0.1432

        bio_t(i,j,k,iPHOSPHATE)=0.149974*(rhp_t(i,j,k)+rho-1000.)-4.1628
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0519) bio_t(i,j,k,iPHOSPHATE)=0.0519

        bio_t(i,j,k,iSILICE)=4.88587*(rhp_t(i,j,k)+rho-1000.)-136.875
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.1707) bio_t(i,j,k,iSILICE)=0.1707
! libye
        elseif(                                  &
        lon_t(i,j)>15.*deg2rad.and.lon_t(i,j)<=23.*deg2rad.and.        &
        lat_t(i,j)<=36.5*deg2rad)then
       
        bio_t(i,j,k,iNITRATE)=10.2623*(rhp_t(i,j,k)+rho-1000.)-239.769
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.1967) bio_t(i,j,k,iNITRATE)=0.1967

        bio_t(i,j,k,iPHOSPHATE)=0.345349*(rhp_t(i,j,k)+rho-1000.)-9.87983
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0328) bio_t(i,j,k,iPHOSPHATE)=0.0328

        bio_t(i,j,k,iSILICE)=10.548*(rhp_t(i,j,k)+rho-1000.)-301.996
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.3305) bio_t(i,j,k,iSILICE)=0.3305

! tun
        elseif(                                  &
        lon_t(i,j)>10.*deg2rad.and.lon_t(i,j)<=15.*deg2rad.and.        &
!        lat_t(i,j)>36.*deg2rad.and.lat_t(i,j)<38.*deg2rad)then
        lat_t(i,j)<=34.*deg2rad)then
       
        bio_t(i,j,k,iNITRATE)=4.8857*(rhp_t(i,j,k)+rho-1000.)-136.453
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.2118) bio_t(i,j,k,iNITRATE)=0.2118

        bio_t(i,j,k,iPHOSPHATE)=0.148005*(rhp_t(i,j,k)+rho-1000.)-4.11025
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0341) bio_t(i,j,k,iPHOSPHATE)=0.0341

        bio_t(i,j,k,iSILICE)=4.752*(rhp_t(i,j,k)+rho-1000.)-132.725
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.2682) bio_t(i,j,k,iSILICE)=0.2682


! ionnienne
        elseif(                                  &
        lon_t(i,j)>15.*deg2rad.and.lon_t(i,j)<22.*deg2rad.and.        &
        lat_t(i,j)>36.5*deg2rad.and.lat_t(i,j)<39.*deg2rad)then
!        lat_t(i,j)<39.*deg2rad)then

        bio_t(i,j,k,iNITRATE)=11.9114*(rhp_t(i,j,k)+rho-1000.)-341.931
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.096) bio_t(i,j,k,iNITRATE)=0.096

        bio_t(i,j,k,iPHOSPHATE)=0.332844*(rhp_t(i,j,k)+rho-1000.)-9.52928
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0348) bio_t(i,j,k,iPHOSPHATE)=0.0348

        bio_t(i,j,k,iSILICE)=11.6992*(rhp_t(i,j,k)+rho-1000.)-335.683
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.2943) bio_t(i,j,k,iSILICE)=0.2943

! cypre
        elseif(lon_t(i,j)>23.*deg2rad.and.        &
!        lat_t(i,j)<36.*deg2rad)then
         lat_t(i,j)<37.5*deg2rad)then
      
        bio_t(i,j,k,iNITRATE)=14.4713*(rhp_t(i,j,k)+rho-1000.)-417.015
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.1656) bio_t(i,j,k,iNITRATE)=0.1656

        bio_t(i,j,k,iPHOSPHATE)=0.55*(rhp_t(i,j,k)+rho-1000.)-16.09
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0232) bio_t(i,j,k,iPHOSPHATE)=0.0232

        bio_t(i,j,k,iSILICE)=14.0977*(rhp_t(i,j,k)+rho-1000.)-406.195
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.2344) bio_t(i,j,k,iSILICE)=0.2344

! NAD
        elseif(                                  &
        lon_t(i,j)>12.*deg2rad.and.lon_t(i,j)<20.*deg2rad.and.        &
        lat_t(i,j)>40.*deg2rad )then

        bio_t(i,j,k,iNITRATE)=16.9461*(rhp_t(i,j,k)+rho-1000.)-488.995
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.0623) bio_t(i,j,k,iNITRATE)=0.0623

        bio_t(i,j,k,iPHOSPHATE)=0.814105*(rhp_t(i,j,k)+rho-1000.)-23.5077
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0371) bio_t(i,j,k,iPHOSPHATE)=0.0371

        bio_t(i,j,k,iSILICE)=15.555*(rhp_t(i,j,k)+rho-1000.)-448.71
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.1188) bio_t(i,j,k,iSILICE)=0.1188

        

! SAD
        elseif(                                 &
        lon_t(i,j)>16.*deg2rad.and.lon_t(i,j)<20.*deg2rad.and.        &
        lat_t(i,j)>39.*deg2rad.and.lat_t(i,j)<40.*deg2rad)then

        bio_t(i,j,k,iNITRATE)=16.0433*(rhp_t(i,j,k)+rho-1000.)-462.552
        if(bio_t(i,j,k,iNITRATE)>=6.5) bio_t(i,j,k,iNITRATE)=6.5
        if(bio_t(i,j,k,iNITRATE)<=0.0942) bio_t(i,j,k,iNITRATE)=0.0942

        bio_t(i,j,k,iPHOSPHATE)=0.5645*(rhp_t(i,j,k)+rho-1000.)-16.2724
        if(bio_t(i,j,k,iPHOSPHATE)>=0.25) bio_t(i,j,k,iPHOSPHATE)=0.25
        if(bio_t(i,j,k,iPHOSPHATE)<=0.0422) bio_t(i,j,k,iPHOSPHATE)=0.0422

        bio_t(i,j,k,iSILICE)=16.0979*(rhp_t(i,j,k)+rho-1000.)-464.141
        if(bio_t(i,j,k,iSILICE)>=6.) bio_t(i,j,k,iSILICE)=6.
        if(bio_t(i,j,k,iSILICE)<=0.122) bio_t(i,j,k,iSILICE)=0.122
! ege
        elseif(                                   &
        lon_t(i,j)>23.*deg2rad.and.                                   &
        lat_t(i,j)>=37.5*deg2rad)then

        bio_t(i,j,k,iNITRATE)=9.27188*(rhp_t(i,j,k)+rho-1000.)-267.101
        if(bio_t(i,j,k,iNITRATE)>=4.) bio_t(i,j,k,iNITRATE)=4.
        if(bio_t(i,j,k,iNITRATE)<=0.2032) bio_t(i,j,k,iNITRATE)=0.2032

        bio_t(i,j,k,iPHOSPHATE)=0.366129*(rhp_t(i,j,k)+rho-1000.)-9.21097
        if(bio_t(i,j,k,iPHOSPHATE)>=0.2) bio_t(i,j,k,iPHOSPHATE)=0.2
        if(bio_t(i,j,k,iPHOSPHATE)<=0.043) bio_t(i,j,k,iPHOSPHATE)=0.043

        bio_t(i,j,k,iSILICE)=8.92964*(rhp_t(i,j,k)+rho-1000.)-257.004
        if(bio_t(i,j,k,iSILICE)>=4.) bio_t(i,j,k,iSILICE)=4.
        if(bio_t(i,j,k,iSILICE)<=0.3416) bio_t(i,j,k,iSILICE)=0.3416


! alboran
        elseif(                                  &
        lon_t(i,j)<=-2.*deg2rad.and.lon_t(i,j)>-5.*deg2rad)then

!        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!        if(bio_t(i,j,k,iNITRATE)<=0.2128) bio_t(i,j,k,iNITRATE)=0.2128
!
!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.0753) bio_t(i,j,k,iPHOSPHATE)=0.0753
!
!        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
!        if(bio_t(i,j,k,iSILICE)<=0.2431) bio_t(i,j,k,iSILICE)=0.2431


! new 28032013 
        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924                 ! original
!         bio_t(i,j,k,iNITRATE)=8.38357*(rhp_t(i,j,k)+rho-1000.)-234.988     ! alg
        if(bio_t(i,j,k,iNITRATE)>=8.) bio_t(i,j,k,iNITRATE)=8.
!        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08981*(rhp_t(i,j,k)+rho-1000.)-30.6768)   ! alg
        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(2.00125*(rhp_t(i,j,k)+rho-1000.)-53.8274)  ! original

        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102    ! original
!         bio_t(i,j,k,iPHOSPHATE)=0.366129*(rhp_t(i,j,k)+rho-1000.)-10.261    ! alg
        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.890801*(rhp_t(i,j,k)+rho-1000.)-27.9247)  ! alg
        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(1.86852*(rhp_t(i,j,k)+rho-1000.)-52.7586)     ! original

        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489                        ! original
!         bio_t(i,j,k,iSILICE)=9.11274*(rhp_t(i,j,k)+rho-1000.)-256.249   ! alg
        if(bio_t(i,j,k,iSILICE)>=7.5) bio_t(i,j,k,iSILICE)=7.5
        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.78793*(rhp_t(i,j,k)+rho-1000.)-48.4092)    !original
!        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.05871*(rhp_t(i,j,k)+rho-1000.)-29.887)  ! alg






!! gibraltar
        elseif(                                      &
        lon_t(i,j)<=-5.*deg2rad)then


!        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!        if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(1.08981*(rhp_t(i,j,k)+rho-1000.)-30.6768)   ! alg

        
        
!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.890801*(rhp_t(i,j,k)+rho-1000.)-27.9247)  ! alg


!        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
!        if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(1.05871*(rhp_t(i,j,k)+rho-1000.)-29.887)  ! alg



!
!        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!!
!!!        if(bio_t(i,j,k,iNITRATE)>=8.5) bio_t(i,j,k,iNITRATE)=8.5
!!!inflow
!        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iNITRATE)=0.4
!        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iPHOSPHATE)=0.12
!        if(rhp_t(i,j,k)+rho-1000.<27.9) bio_t(i,j,k,iSILICE)=0.79
!!!outflow
!        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iNITRATE)=8.
!        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(rhp_t(i,j,k)+rho-1000.>28.6) bio_t(i,j,k,iSILICE)=6.47
        bio_t(i,j,k,iNITRATE)=3.19339*(rhp_t(i,j,k)+rho-1000.)-84.8924
       if(bio_t(i,j,k,iNITRATE)>=10.) bio_t(i,j,k,iNITRATE)=10.
!       if(bio_t(i,j,k,iNITRATE)<=0.032) bio_t(i,j,k,iNITRATE)=0.032
      if(bio_t(i,j,k,iNITRATE)<=0.14) bio_t(i,j,k,iNITRATE)=0.14


        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
       if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
!       if(bio_t(i,j,k,iPHOSPHATE)<=0.02) bio_t(i,j,k,iPHOSPHATE)=0.02
       if(bio_t(i,j,k,iPHOSPHATE)<=0.03) bio_t(i,j,k,iPHOSPHATE)=0.03



        bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
       if(bio_t(i,j,k,iSILICE)>=9.) bio_t(i,j,k,iSILICE)=9.
!       if(bio_t(i,j,k,iSILICE)<=0.03) bio_t(i,j,k,iSILICE)=0.03
       if(bio_t(i,j,k,iSILICE)<=0.12) bio_t(i,j,k,iSILICE)=0.12










!        if(bio_t(i,j,k,iNITRATE)<=0.2128) bio_t(i,j,k,iNITRATE)=0.2128
!
!        bio_t(i,j,k,iPHOSPHATE)=0.135093*(rhp_t(i,j,k)+rho-1000.)-3.55102
!        if(bio_t(i,j,k,iPHOSPHATE)>=0.4) bio_t(i,j,k,iPHOSPHATE)=0.4
!        if(bio_t(i,j,k,iPHOSPHATE)<=0.0753) bio_t(i,j,k,iPHOSPHATE)=0.0753
!       bio_t(i,j,k,iSILICE)=3.48983*(rhp_t(i,j,k)+rho-1000.)-93.489
!        if(bio_t(i,j,k,iSILICE)>=8.) bio_t(i,j,k,iSILICE)=8.
!        if(bio_t(i,j,k,iSILICE)<=0.2431) bio_t(i,j,k,iSILICE)=0.2431





! bioregion nord occidental (provençale + ligure)

       else
         


         bio_t(i,j,k,iNITRATE)=1./(36.0319 - (1.23423*(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iNITRATE)>=9.2) bio_t(i,j,k,iNITRATE)=9.2
      if(bio_t(i,j,k,iNITRATE)<=1.) bio_t(i,j,k,iNITRATE)=exp(0.690112*(rhp_t(i,j,k)+rho-1000.)-19.7654)

!         write(20+par%rank,*),i,j,k,rho,rhp_t(i,j,k)
!         print*,'obc_bio',i,j,k,rho,rhp_t(i,j,k)
         bio_t(i,j,k,iPHOSPHATE)=exp(164.933 - (4828.84/(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iPHOSPHATE)>=0.45) bio_t(i,j,k,iPHOSPHATE)=0.45
       if(bio_t(i,j,k,iPHOSPHATE)<=0.06) bio_t(i,j,k,iPHOSPHATE)=exp(0.560291*(rhp_t(i,j,k)+rho-1000.)-18.9134)

         bio_t(i,j,k,iSILICE)= 1./(-37.9453 + (1108.41/(rhp_t(i,j,k)+rho-1000.)))
         if(bio_t(i,j,k,iSILICE)>=9) bio_t(i,j,k,iSILICE)=9
       if(bio_t(i,j,k,iSILICE)<=1.) bio_t(i,j,k,iSILICE)=exp(0.769009*(rhp_t(i,j,k)+rho-1000.)-21.6826)


       endif
         salz=(sal_t(i,j,k,1) - salmu )/salsigma
         bio_t(i,j,k,iOXYGEN) = p1*salz**4   &
                              + p2*salz**3   &
                              + p3*salz**2   &
                              + p4*salz      &
                              + p5

         if(bio_t(i,j,k,iOXYGEN)>=260) bio_t(i,j,k,iOXYGEN)=260
         if(bio_t(i,j,k,iOXYGEN)<=170) bio_t(i,j,k,iOXYGEN)=170

 
! Gradient nul "z1"
         bio_t(i1,j,k,vb)=bio_t(i,j,k,vb)

         endif ! mask
         endif ! mask
         enddo ! k ajout caroline
         enddo !j

  987    continue
         enddo ! k1

      endif ! biobc_type

!! Nutrients deduced from density in the sponge layers
       if(biobc_type(vb)==6) then

      endif ! biobc_type

!! Nutrients deduced observed profils for lon<-8
       if(biobc_type(vb)==7) then

      call obc_bio_upd

      endif ! biobc_type


      end subroutine usersbiobc
