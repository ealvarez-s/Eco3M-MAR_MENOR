










      module module_sedw
      implicit none

      real*4,dimension(:,:,:),allocatable ::  &
       rhos                                   &
      ,fr0
      real*4,dimension(:,:),allocatable ::  &
       db             &
      ,tauw           &
      ,tauc           &
      ,taucw          &
      ,taucwmax       &
      ,ustarc         &
      ,ustarw         &
      ,ustarcw        &
      ,ustarcwmax     &
      ,tauws          &
      ,taucs          &
      ,taucws         &
      ,taucwmaxs      &
      ,ustarcs        &
      ,ustarws        &
      ,ustarcws       &
      ,ustarcwmaxs    &
      ,ksf            &
      ,dm             &
      ,zta_rip        &
      ,lam_rip        &
      ,fr_s           &
      ,fr_v           &
      ,z0_sf          &
      ,poro           &
      ,frc            &
      ,tetacr         &
      ,tetas          &
      ,kb_z
      real*4,dimension(:),allocatable ::  &
       hfdmax         &
      ,waveinfo       &
      ,phi_sed
      integer,dimension(:),allocatable ::  &
       datewave
      real*4      &
       zta0       &
      ,tstar      &
      ,taucr      &
      ,pm         &
      ,dmc        &
      ,rhosc      &
      ,poroc      &
      ,lam0       &
      ,hwc        &
      ,twc        &
      ,dwc        &
      ,z0_bfc

      integer     &
       iwave_demo &
      ,idia       &
      ,iz0        &
      ,ni         &   !=1 si forcage par modele de vague
      ,nfd        &   !nombre de fond avec fraction
      ,kfd

      character         &
       wavefile(3)*60   &
      ,diafichier*44    &
      ,dmfichier*43     &
      ,frfichier*48     &
      ,rhosedfichier*48 &
      ,porofichier*45

      real*4 ::          &
       soulsby_a1=-0.06  &
      ,soulsby_a2= 1.7   &
      ,soulsby_a3=-0.29  &
      ,soulsby_a4= 0.29  &
      ,soulsby_m1= 0.67  &
      ,soulsby_m2=-0.29  &
      ,soulsby_m3= 0.09  &
      ,soulsby_m4= 0.42  &
      ,soulsby_n1= 0.75  &
      ,soulsby_n2=-0.27  &
      ,soulsby_n3= 0.11  &
      ,soulsby_n4=-0.02  &
      ,soulsby_ci= 0.80

contains

!...............................................................

      subroutine sedw_allocate
      use module_principal
      implicit none

      allocate(       rhos(0:imax+1,0:jmax+1,vbmax)) ; rhos=0
      allocate(        fr0(0:imax+1,0:jmax+1,vbmax)) ; fr0=0

      allocate(         db(0:imax+1,0:jmax+1)) ; db=0
      allocate(       tauw(0:imax+1,0:jmax+1)) ; tauw=0
      allocate(       tauc(0:imax+1,0:jmax+1)) ; tauc=0
      allocate(      taucw(0:imax+1,0:jmax+1)) ; taucw=0
      allocate(   taucwmax(0:imax+1,0:jmax+1)) ; taucwmax=0
      allocate(     ustarc(0:imax+1,0:jmax+1)) ; ustarc=0
      allocate(    ustarcw(0:imax+1,0:jmax+1)) ; ustarcw=0
      allocate( ustarcwmax(0:imax+1,0:jmax+1)) ; ustarcwmax=0
      allocate(      tauws(0:imax+1,0:jmax+1)) ; tauws=0
      allocate(      taucs(0:imax+1,0:jmax+1)) ; taucs=0
      allocate(     taucws(0:imax+1,0:jmax+1)) ; taucws=0
      allocate(  taucwmaxs(0:imax+1,0:jmax+1)) ; taucwmaxs=0
      allocate(    ustarcs(0:imax+1,0:jmax+1)) ; ustarcs=0
      allocate(    ustarws(0:imax+1,0:jmax+1)) ; ustarws=0
      allocate(   ustarcws(0:imax+1,0:jmax+1)) ; ustarcws=0
      allocate(ustarcwmaxs(0:imax+1,0:jmax+1)) ; ustarcwmaxs=0
      allocate(        ksf(0:imax+1,0:jmax+1)) ; ksf=0
      allocate(         dm(0:imax+1,0:jmax+1)) ; dm=0
      allocate(    zta_rip(0:imax+1,0:jmax+1)) ; zta_rip=0
      allocate(    lam_rip(0:imax+1,0:jmax+1)) ; lam_rip=0
      allocate(       fr_s(0:imax+1,0:jmax+1)) ; fr_s=0
      allocate(       fr_v(0:imax+1,0:jmax+1)) ; fr_v=0
      allocate(      z0_sf(0:imax+1,0:jmax+1)) ; z0_sf=0
      allocate(       poro(0:imax+1,0:jmax+1)) ; poro=0
      allocate(     tetacr(0:imax+1,0:jmax+1)) ; tetacr=0
      allocate(      tetas(0:imax+1,0:jmax+1)) ; tetas=0
      allocate(       kb_z(0:imax+1,0:jmax+1)) ; kb_z=0

      allocate(        frc(vbmax,3)) ; frc=0
      allocate(    phi_sed(vbmax))   ; phi_sed=0

      allocate(     hfdmax(0:4)) ; hfdmax=0
      allocate(   datewave(6))   ; datewave=0
      allocate(   waveinfo(2))   ; waveinfo=0

      end subroutine sedw_allocate

!...............................................................

      subroutine sedw_tauc
      use module_principal
      implicit none
      real  cdb_

      do j=1,jmax
      do i=1,imax
       cdb_=(0.4/log((h_w(i,j)+depth_t(i,j,kmin_w(i,j)))/z0_w(i,j)))**2
        tauc(i,j)=rho*cdb_*(0.5*( vel_u(i+1,j  ,kmin_u(i+1,j  ),1)**2  &!tension de cisaillement due aux courants
                                 +vel_u(i  ,j  ,kmin_u(i  ,j  ),1)**2  &
                                 +vel_v(i  ,j+1,kmin_v(i  ,j+1),1)**2  &
                                 +vel_v(i  ,j  ,kmin_v(i  ,j  ),1)**2 ))
      enddo
      enddo

      end subroutine sedw_tauc

!...............................................................

      subroutine sedw_tauc_tauw_taucw
      use module_principal
      implicit none
      real  a_, m_, n_, x_, gy_, phicw, phic,cdb_,log10fwovercd_

      do j=1,jmax
      do i=1,imax

       cdb_=(0.4/log((h_w(i,j)+depth_t(i,j,kmin_w(i,j)))/z0_w(i,j)))**2
        tauc(i,j)=rho*cdb_*(0.5*( vel_u(i+1,j  ,kmin_u(i+1,j  ),1)**2  &!tension de cisaillement due aux courants
                                 +vel_u(i  ,j  ,kmin_u(i  ,j  ),1)**2  &
                                 +vel_v(i  ,j+1,kmin_v(i  ,j+1),1)**2  &
                                 +vel_v(i  ,j  ,kmin_v(i  ,j  ),1)**2 ))

       log10fwovercd_=log10(fw(i,j)/cdb_)

       tauw(i,j)=max(0.5*rho*fw(i,j)*ubw(i,j)**2,0.0001d00)        ! tension de cisaillement due aux vagues

       taucw(i,j)=tauc(i,j)*(1.+ 1.2*(( tauw(i,j)/ &
                        (tauw(i,j)+tauc(i,j)+small2) + small2)**3.2)) ! soulsby 1995

! calcul de l'angle courant/vague
          phic=atan2( vel_v(i  ,j  ,kmin_v(i  ,j  ),1)     &
                     +vel_v(i  ,j+1,kmin_v(i  ,j+1),1)     &
                    , vel_u(i  ,j  ,kmin_u(i  ,j  ),1)     &
                     +vel_u(i+1,j  ,kmin_u(i+1,j  ),1) )

          phicw=abs(dir_wave_t(i,j,1)-phic)                     !01-02-12


      a_=(soulsby_a1+soulsby_a2*abs((cos(phicw)))**soulsby_ci)      &
        +(soulsby_a3+soulsby_a4*abs((cos(phicw)))**soulsby_ci)      &
         *log10fwovercd_
      m_=(soulsby_m1+soulsby_m2*abs((cos(phicw)))**soulsby_ci)      &
        +(soulsby_m3+soulsby_m4*abs((cos(phicw)))**soulsby_ci)      &
         *log10fwovercd_
      n_=(soulsby_n1+soulsby_n2*abs((cos(phicw)))**soulsby_ci)      &
        +(soulsby_n3+soulsby_n4*abs((cos(phicw)))**soulsby_ci)      &
         *log10fwovercd_

! paramètres adimensionnels
       x_= tauc(i,j)/(tauc(i,j)+tauw(i,j))
      gy_=max(zero,un+a_*x_**m_*(un-x_)**n_)

      taucwmax(i,j)=gy_*(tauc(i,j)+tauw(i,j))

      enddo
      enddo

      end subroutine sedw_tauc_tauw_taucw

!...............................................................

      subroutine sedw_taucs_tauws_taucws
      use module_principal
      implicit none
      real  a_, m_, n_, x_, gy_, phicw, phic,cdb_,log10fwovercd_

      do j=1,jmax
      do i=1,imax

       cdb_=(0.4/log((h_w(i,j)+depth_t(i,j,kmin_w(i,j)))/z0_sf(i,j)))**2
       taucs(i,j)=rho*cdb_*(0.5*( vel_u(i+1,j  ,kmin_u(i+1,j  ),1)**2  &!tension de cisaillement due aux courants
                                 +vel_u(i  ,j  ,kmin_u(i  ,j  ),1)**2  &
                                 +vel_v(i  ,j+1,kmin_v(i  ,j+1),1)**2  &
                                 +vel_v(i  ,j  ,kmin_v(i  ,j  ),1)**2 ))

       log10fwovercd_=log10(fw(i,j)/cdb_)

       tauws(i,j)=max(0.5*rho*fw(i,j)*ubw(i,j)**2,0.0001d00)        ! tension de cisaillement due aux vagues

       taucws(i,j)=taucs(i,j)*(1.+ 1.2*(( tauws(i,j)/ &
                        (tauws(i,j)+taucs(i,j)+small2) + small2)**3.2)) ! soulsby 1995

! calcul de l'angle courant/vague
          phic=atan2( vel_v(i  ,j  ,kmin_v(i  ,j  ),1)     &
                     +vel_v(i  ,j+1,kmin_v(i  ,j+1),1)     &
                    , vel_u(i  ,j  ,kmin_u(i  ,j  ),1)     &
                     +vel_u(i+1,j  ,kmin_u(i+1,j  ),1) )

          phicw=abs(dir_wave_t(i,j,1)-phic)                     !01-02-12


      a_=(soulsby_a1+soulsby_a2*abs((cos(phicw)))**soulsby_ci)      &
        +(soulsby_a3+soulsby_a4*abs((cos(phicw)))**soulsby_ci)      &
         *log10fwovercd_
      m_=(soulsby_m1+soulsby_m2*abs((cos(phicw)))**soulsby_ci)      &
        +(soulsby_m3+soulsby_m4*abs((cos(phicw)))**soulsby_ci)      &
         *log10fwovercd_
      n_=(soulsby_n1+soulsby_n2*abs((cos(phicw)))**soulsby_ci)      &
        +(soulsby_n3+soulsby_n4*abs((cos(phicw)))**soulsby_ci)      &
         *log10fwovercd_

! paramètres adimensionnels
       x_= taucs(i,j)/(taucs(i,j)+tauws(i,j))
      gy_=max(zero,un+a_*x_**m_*(un-x_)**n_)

      taucwmaxs(i,j)=gy_*(taucs(i,j)+tauws(i,j))

      enddo
      enddo

      end subroutine sedw_taucs_tauws_taucws

!..................................................

      end module module_sedw
