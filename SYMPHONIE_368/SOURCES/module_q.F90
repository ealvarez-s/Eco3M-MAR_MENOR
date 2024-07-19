      module module_q
      use module_principal ; use module_parallele ; use module_wave
      implicit none
!______________________________________________________________________
! SYMPHONIE ocean model
! release 347 - last update: 01-05-22
!______________________________________________________________________

!...............................................................................
! Version date      Description des modifications
!         xx-07-17  mise en service
!         10-11-17  mc_out remplcE par q_t(2)
!         23-11-17  correction sig1dpom_t
!         05-12-17  temps entre dans unites du coef de wave breaking
!         09-12-17  lecture notebook_nh deplacE
!         17-08-18  NH avec T et S
!         18-08-18  mise A jour du terme de diag d'equilibre NH
!         20-08-18  nouvelle alternance des champs: q(0)  u(1)  q(1)  u(2)  q(2)
!         26-08-18  ne pas utiliser le tableau "wetmask" pour le wetdrying de q pour
!                   ne pas etre en conflit avec le wetmask des traceurs
!         03-09-18  dbefore et cie definis comme des "parameter"
!         10-10-18  limite superieure sur le coef de wave breaking3d
!         21-10-18  ajout obc pour q en j=jmax, j=1
!         14-02-19  ajout acm_speed
!         11-03-19  ajout obc_q_i
! v259    01-10-19  Mise A jour des conditions aux limites
! v261    21-10-19  wetmask_u, wetmask_v
! v280    02-05-20  - schema calcul nhpgf_u nhpgf_v
!                   - nouvelle formule du breaker_t(i,j)
! v326    26-12-21  cas coordonnee VQS !26-12-21
! v347    01-05-22  ContinuitE mpi deplacEe en debut de boucle
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


! DESCRIPTION: https://docs.google.com/document/d/1wU5DVBww6YuYZ2VK0mAmTfH0rNPI97AOF46PuEIfFso/edit

      integer (kind=1) :: &
!                          flag_nh2d      & ! 0=NH off, 1=NH on , 2D case
!                         ,flag_nh3d_uv   & ! 0=NH off, 1=NH on 3D barotropic case, 2=NH on with T,S
                flag_nhdiag_timeaverage=0 &
               ,flag_surf_obc_q=0         & ! 0 pour dirichet, Neuman sinon
!                         ,nh_scheme=1      ! "leplussimple"
!                         ,nh_scheme=2      ! "sigma"
                          ,nh_scheme=3      ! "remapping" correction dq/dz evaluee au point central
!                         ,nh_scheme=4      ! "remapping" correction dq/dz evaluee  en i+1 (resp i-1,j+1,j-1)

!     real , dimension(:,:), allocatable ::   &
!      breaker_t

      integer , dimension(:), allocatable ::   &
             nhconvergence_histo1              &
            ,nhconvergence_histo2
      double precision, dimension(:), allocatable ::   &
             stokes_transport_sum1                     &
            ,stokes_transport_sum2                     &
            ,stokes_transport_sum3                      
      double precision :: check_sshwave=0.,dtmultiple=-999.
      real, dimension(:,:), allocatable ::   &
           vqsratio_t

      integer :: nh_freqrestart=10000 &
                ,nh_period_modulo=0   &
                ,retarded_pss_nb=4
      integer , private :: loop,loopmax=1,ktop &
!                         ,id_dzpgfu=2         & ! identifiant DZ*(dp/dx-dz/dx*dp/dk)
!                         ,id_dzpgfv=3         & ! identifiant DZ*(dp/dy-dz/dy*dp/dk)
                          ,id_dzodxpgfu=2      & ! identifiant dz/dx*(dp/dx-dz/dx*dp/dk)
                          ,id_dzodypgfv=3      & ! identifiant dz/dy*(dp/dy-dz/dy*dp/dk)
                          ,stokes_count=0      &
                          ,id_sshwave=3        &    
                          ,id_hydro_du=1       &
                          ,id_hydro_dv=2       &
                          ,computefullsolver=1    
      real :: dt2=0.  , damp=0.0    , c1_ , c2_ ,dtconvert,wavebreakfactor=2. &
             ,wavebreakcoef                   &
             ,w0_=0.5,w1_=0.,w2_=0.5,ww_            ! 50% explicite sur t-1, 50% implicite sur t+1
!            ,w0_=0.25,w1_=0.5,w2_=0.25,ww_      &  ! 25% explicite sur t-1, 50% explicite sur t, 25% implicite sur t+1
!            ,w0_=0.,w1_=0.,w2_=1.,ww_           &  ! 100% implicite sur t+1
!            ,w0_=0.4 ,w1_=0.,w2_=0.6 ,ww_       &  
!            ,nhpgf_reduce=0.98

      integer(kind=1) :: computesigma=1       &
                        ,computematrix=1       

#ifdef bidon
      real*16 :: sum_t_qdq=0.      &
                ,sum_t_udssh=0.    &
                ,sum_t_udu=0.      &
                ,sum_t_udq_h=0.    &
                ,sum_t_udq_v1=0.   &
                ,sum_t_udq_v2=0.   &
                ,sum_t_udq_err=0.  &
                ,sum_t_sshdssh=0.  &
                ,sum_ssh2_init=0.  &
                ,sum_ssh2_now=0.   &
                ,sum_q2_init=0.    &
                ,sum_q2_now=0.     &
                ,sum_u2_init=0.    &
                ,sum_u2_now=0.
#endif
!#ifdef bidon
! declaration pour subroutine q_poisson_ondes_diag
      integer iter_loop_
      double precision ::          &
                period_=-999.      & ! periode                   secondes
               ,pulsation_         & ! 2pi/periode               rad/s
               ,kvector_=-999.     & ! vecteur d'onde horizontal m**-1
               ,c_=0.              & ! vitesse de phase          m/s      pulsation_/kvector_
               ,beta_,alpha_,grav_over_2pi_,delta_c_,u_amplitude_,ssh_amplitude_
!     real , dimension(:) , allocatable :: kvector_oe, kvector_sn

!     real , dimension(:)   , allocatable :: kvector_oi, kvector_oj
!     double precision , dimension(:,:) , allocatable :: phase_wave
!     double precision valobc1,valobc2

!#endif
      
contains

!...............................

      subroutine q_initial
      implicit none

!     if(flag_nh3d>=1) then !ooo> !17-08-18
!      if(igesig/=0) &
!      stop 'Err 129 si NH alors sigma simple dans notebook_vertcoord'
!     endif                 !ooo>

      ema_mu=dti_fw/2. ! dti_fw/2s

      if(flag_nh3d==0)return

!     call wave_ww3_spec_sete_read
      call cwave_int_shortwaves

! Estimation de c_ et kvector_ pour h=hmax (oU l'on suppose que se situe les OBC)
    
!     period_=10. ! 9. ! 10. ! 250*dte_fw 
!     open(unit=3,file='periode_input')
!     read(3,*)period_
!     close(3)

      period_=10. ! Hamilton 
      periodpeak=period_
      freq2pipeak=2.*pi/periodpeak
       pulsation_=2.*pi/periodpeak
      grav_over_2pi_=grav/(2.*pi)

! wavebreakcoef=(3/8)*grav**3/2*dH/dx*(T/pi)**2 / 0.1
!     wavebreakcoef=(3./8.)*(grav**1.5) &
!                  *(1./30.)            & ! pente moyenne de la plage  hamilton
!                  *(periodpeak/pi)**2
      wavebreakcoef=1.

!     RETURN

!     write(6,*)'wavebreakcoef=',wavebreakcoef
!     write(6,*)'period_=',period_
!     stop 'coco'


! c_ is the phase speed
! On veut résoudre x=f(x)=alpha*tanh(beta/x) la vitesse de phase c_  etant l'inconnue x
! Le raisonnement est le suivant: f(x+dx)=x+dx
!                           soit: f(x)+dx f'(x)=x+dx
!                           soit: dx=(x-f(x))/(f'(x)-1)
! On estime une première valeur de x=alpha_ l'approximation pour c_ aux grandes profondeurs
      beta_=2.*pi*max(hmax,0.01d0)/period_
      alpha_=grav_over_2pi_*period_                ! c=g*t/2pi first guess for large depth

      if(beta_/alpha_>200.) then !fgfgfgfgfgfgfg>      !07-05-11

! Cas c=first guess valide:
       c_=alpha_ ! phase speed for large depth

      else                       !fgfgfgfgfgfgfg>

! Cas calcul itératif:
       delta_c_=0.01                  ! dx, on estime une première valeur de dx
       c_=alpha_+delta_c_             ! alpha_+dx
         do iter_loop_=1,200

          delta_c_=(c_-alpha_*tanh(beta_/c_))                & !  (x-f(x))
                  /(-alpha_*beta_/(c_*cosh(beta_/c_))**2-1.)   !  /(f'(x)-1)
          c_=max(c_+delta_c_,1.d-4)
          if(abs(delta_c_)<0.01)goto 25
         enddo
         write(6,*)c_,delta_c_
         stop 'waveforcing: le calcul de kn ne converge pas'
   25  continue
      endif                      !fgfgfgfgfgfgfg>

      kvector_=2.*pi/(c_*period_)
      kvectorpeak=kvector_
      
      if(par%rank==0) then !OOO>
       open(unit=3,file='tmp/messages',position='append')
        write(3,*)'..............................'
        write(3,*)'subroutine q_initial'
        write(3,*)'period_        ',period_
        write(3,*)'Parametres en suivant pour h=',hmax
        write(3,*)'  kvector_       ',kvector_
        write(3,*)'  Longueur d onde',2.*pi/kvector_
        write(3,*)'  ssh_amplitude_ ',ssh_amplitude_
        write(3,*)'  c_             ',c_
        write(3,*)'wavebreakcoef=',wavebreakcoef
       close(3)
      endif                !OOO>

      do i=1,imax 
       sqr_hoverg_v(i,2)=h_v(i,jmax)/max(c_,1.)
       sqr_hoverg_v(i,1)=h_v(i,2   )/max(c_,1.)
      enddo
      do j=1,jmax 
       sqr_hoverg_u(j,2)=h_u(imax,j)/max(c_,1.)
       sqr_hoverg_u(j,1)=h_u(2   ,j)/max(c_,1.)
      enddo
      cwi_int_u=c_
      cwi_int_v=c_
      cwj_int_v=c_
      cwj_int_u=c_

!     stop 'coco dingo'
      RETURN

#ifdef bidon
! Calculer les composantes du vecteur d'onde pour une condition limite en j=jmax
! Forcage de la forme sin(wt+k1*x+k2*y) 
! Note: la solution est entrante si elle se propage A contresens de Oj,
! autrement dit verifier qu'on a bien k2>0
      x0=90.                ! Angle de provenance
      x0=(270.-x0)*deg2rad   ! changement de convention repere OE,SN
      x1=kvector_*cos(x0)    ! kvector_oe
      x2=kvector_*sin(x0)    ! kvector_sn

! COPMPOSANTES VECTEUR D'ONDE
      j=jmax
      do i=1,imax
       kvector_oi(i)=x1*gridrotcos_t(i,j)-x2*gridrotsin_t(i,j)
       kvector_oj(i)=x1*gridrotsin_t(i,j)+x2*gridrotcos_t(i,j)
!      if(jmax+par%tjmax(1)==jglb)write(10+par%rank,*)i+par%timax(1),kvector_oi(i),kvector_oj(i)
      enddo
! Le domaine periodise a pour consequence que les valeurs de kvector_oi
! et kvector_oj ne sont pas correctes en i=1,2,iglb-1,iglb
! reparation etape 1: extrapolation
      if(1   +par%timax(1)==1)   then  !>>>
        kvector_oi(2)=2*kvector_oi(3)-kvector_oi(4)
        kvector_oi(1)=2*kvector_oi(2)-kvector_oi(3)
        kvector_oj(2)=2*kvector_oj(3)-kvector_oj(4)
        kvector_oj(1)=2*kvector_oj(2)-kvector_oj(3)
      endif                            !>>>
      if(imax+par%timax(1)==iglb) then !oo>
        kvector_oi(imax-1)=2*kvector_oi(imax-2)-kvector_oi(imax-3)
        kvector_oi(imax  )=2*kvector_oi(imax-1)-kvector_oi(imax-2)
        kvector_oj(imax-1)=2*kvector_oj(imax-2)-kvector_oj(imax-3)
        kvector_oj(imax  )=2*kvector_oj(imax-1)-kvector_oj(imax-2)
      endif                            !oo>

! Assurer la periodicite de kvector_oj  (celle de kvector_oi n etant a priori pas utile)
       valobc1=0.
       valobc2=0.
       if(1   +par%timax(1)==1   .and.jmax+par%tjmax(1)==jglb)valobc1=kvector_oj(2   )
       if(imax+par%timax(1)==iglb.and.jmax+par%tjmax(1)==jglb)valobc2=kvector_oj(imax)
       call mpi_allreduce(valobc1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(valobc2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       valobc1=sum1glb ; valobc2=sum2glb
! On corrige progressivement toute la frontiere:
       x10=valobc1-valobc2
       do i=1,imax
        i1=i+par%timax(1)
        rap=real(i1-2)/real(iglb-2)
        kvector_oj(i)=kvector_oj(i)+rap*x10 ! correction
       enddo
! Pour faire les choses strictement (au dernier chiffre apres la virgule) en i=iglb:
       if(imax+par%timax(1)==iglb)kvector_oj(imax)=valobc1
! Pour faire les choses strictement, envoyer la phase en i=iglb-1 A i=1:
       valobc2=0.
       if(imax+par%timax(1)==iglb.and.jmax+par%tjmax(1)==jglb) valobc2=kvector_oj(imax-1)
       call mpi_allreduce(valobc2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       valobc2=sum2glb
       if(1   +par%timax(1)==1)kvector_oj(1)=valobc2
!      do i=1,imax
!       if(jmax+par%tjmax(1)==jglb)write(10+par%rank,*)i+par%timax(1),kvector_oi(i),kvector_oj(i)
!      enddo

! PHASE
      do i=1,imax
       j=jmax
       phase_wave(i,1)=x1*rayonterre*lon_t(i,j)*cos(phi0)+x2*rayonterre*lat_t(i,j)
       j=jmax+1
       phase_wave(i,2)=x1*rayonterre*lon_v(i,j)*cos(phi0)+x2*rayonterre*lat_v(i,j)
!      if(jmax+par%tjmax(1)==jglb) then
!       write(10+par%rank,*)i+par%timax(1),phase_wave(i,1),phase_wave(i,2)
!      endif

      enddo
! extrapoler en i=1 et i=iglb (car avec la periodicite lon,lat en i=1 et i=iglb sont discontinus)
      if(1   +par%timax(1)==1)   phase_wave(1   ,:)=2*phase_wave(2     ,:)-phase_wave(3     ,:)
      if(imax+par%timax(1)==iglb)phase_wave(imax,:)=2*phase_wave(imax-1,:)-phase_wave(imax-2,:)
!     do i=1,imax
!     if(jmax+par%tjmax(1)==jglb)write(10+par%rank,*)i+par%timax(1),phase_wave(i,:)*rad2deg
!     enddo

!#ifdef bidon

! Apres extrapolation passer en modulo (en effet on n'extrapole pas un modulo)
      do i=1,imax
       phase_wave(i,:)=modulo(phase_wave(i,:),2*pi)
      enddo

! Assurer la periodicite de phase_wave
      do j=1,2  !>>>>
       valobc1=0.
       valobc2=0.
       if(1   +par%timax(1)==1   .and.jmax+par%tjmax(1)==jglb)valobc1=phase_wave(2   ,j)
       if(imax+par%timax(1)==iglb.and.jmax+par%tjmax(1)==jglb)valobc2=phase_wave(imax,j)
       call mpi_allreduce(valobc1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(valobc2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       valobc1=sum1glb ; valobc2=sum2glb
! periodicite si phase(i=2)=phase(iglb) soit valobc1=valobc2, l'erreur valobc2-valobc1 est maintenant connus de tous les domaines
! Corriger au plus court (ex: 4 degres plutot que -356 degres)
       x10=valobc1-valobc2
       if(x10> pi)x10=x10-2.*pi
       if(x10<-pi)x10=x10+2.*pi
! On corrige progressivement toute la frontiere:
       do i=1,imax
        i1=i+par%timax(1)
        rap=real(i1-2)/real(iglb-2)
        phase_wave(i,j)=phase_wave(i,j)+rap*x10 ! correction
       enddo

! Pour faire les choses strictement (au dernier chiffre apres la virgule) en i=iglb:
       if(imax+par%timax(1)==iglb)phase_wave(imax,j)=valobc1

! Pour faire les choses strictement, envoyer la phase en i=iglb-1 A i=1:
       valobc2=0.
       if(imax+par%timax(1)==iglb.and.jmax+par%tjmax(1)==jglb) valobc2=phase_wave(imax-1,j)
       call mpi_allreduce(valobc2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
       valobc2=sum2glb
       if(1   +par%timax(1)==1)phase_wave(1,j)=valobc2

      enddo     !>>>> fin de boucle sur j

#endif

#ifdef bidon
! Avec ces lignes on verifie que l'angle ajuste n'est pas loin de
! l'angle idel. Attention que ce calcul embrionaire n'est pas toujours
! valide, quand la phase (passee au modulo) est de part et d'autre de
! 360.
      j=jmax
      do i=2,imax-1
       if(jmax+par%tjmax(1)==jglb) then
! L'angle reel une fois l'ajustement fait: 
! deduire kvector_oi kvector_oj de la derivee de la phase
       x3=(phase_wave(i,1)-phase_wave(i-1,1))/dx_t(i,j)       ! kvector_oi
       x4=(phase_wave(i,2)-phase_wave(i  ,1))/(0.5*dy_t(i,j)) ! kvector_oi
       x5= x3*gridrotcos_t(i,j)+x4*gridrotsin_t(i,j)          !  kvector_oe
       x6=-x3*gridrotsin_t(i,j)+x4*gridrotcos_t(i,j)          !  kvector_sn
       write(10+par%rank,*)i+par%timax(1),real(mod(x0*rad2deg,360.)),real( mod(atan2(x6,x5)*rad2deg,360.) )
       endif
      enddo
#endif

!     do j=1,jmax
!     do i=1,imax
!      vel_u(i,j,:,:)=kvector_oi(i)
!      vel_v(i,j,:,:)=kvector_oj(i)
!     enddo
!     enddo
!     do j=1,jmax
!     do i=1,imax
!      ssh_int_w(i,j,:)=cos(kvector_oe(i)*real(i+par%timax(1))*dxb+kvector_sn(i)*real(j+par%tjmax(1))*dxb)
!      ssh_int_w(i,j,:)=cos(kvector_oe(i)*rayonterre*lon_t(i,j)*cos(lat_t(i,j)) &
!                          +kvector_sn(i)*rayonterre*lat_t(i,j) )
!      ssh_int_w(i,j,:)=cos(x1*rayonterre*lon_t(i,j)*cos(lat_t(i,j)) &
!                          +x2*rayonterre*lat_t(i,j) )
!      ssh_int_w(i,j,:)=cos(phase_wave(i,1))
!     enddo
!     enddo
!     call graph_out
!     stop 'coco'

      end subroutine q_initial

!...............................

      subroutine q_initial_beforetimestep
      implicit none

! notebook_nh part 2
      namelist/notebook_nh2/dtmultiple

! notebook_nh !09-12-17
      namelist/notebook_nh/flag_nh3d  & ! 0=NH off, 1=NH on
                          ,asselin_nh        &
                          ,cfl_nh            &
                          ,nh_freqrestart    &
                          ,wetdry_cstnh      &
                          ,nh2d_graph_period &
                          ,nh2d_graph_spinup &
                          ,nh_frozensigma    &
                          ,wavebreakfactor   &
                          ,flag_nhdiag_timeaverage &
                          ,nhpgf_reduce            &
                          ,nh_wavebreakscheme      &
                          ,ssh_amplitude_          &
                          ,brk_crit_slope          & 
                          ,brk_crit_h              & 
                          ,brk_crit_r              &       
                          ,brk_crit_var
            

      dtmultiple=-999.
      open(3,file=nomfichier(34))    ! lecture du namelist "notebook_nh" 
      read(3,nml=notebook_nh2)
      close(3)

      ssh_amplitude_=1.
      flag_nh3d=0
      asselin_nh=0.01
      cfl_nh=0.9
      wetdry_cstnh=0.01
      nh_freqrestart=10000
      nh2d_graph_period=1000
      nh2d_graph_spinup=0
      nh_frozensigma=0
      flag_nhdiag_timeaverage=0
      wavebreakfactor=2.
      nhpgf_reduce=1.
      nh_wavebreakscheme=0
      brk_crit_slope=0.28 
      brk_crit_var=2.8 
      brk_crit_h=0.8   
      brk_crit_r=0.75

      open(3,file=nomfichier(34))    ! lecture du namelist "notebook_nh" 
      read(3,nml=notebook_nh)
      close(3)

      inv_brkslope2=(1./brk_crit_slope)**2
      inv_brkvar=    1./brk_crit_var    
      inv_brkh=1./(0.5*brk_crit_h)

      if(flag_nh3d==flag_nh3d_nosplit_uv.or. &
         flag_nh3d==flag_nh3d_nosplit_tsuv) then !m°v°m>

! Desarmer le time-splitting
       flag_timesplitting=0

! Cas SANS time-splitting :
       dbefore=-1 
       dnow=0    
       dafter=1    
! Note: dans le cas SANS time splitting l'association des facteurs et des vitesses se fait de maniere analogue
! au module_external, A savoir dz(-1) avec vel(0) etc...
      endif                                      !m°v°m>

      if(flag_nh3d==flag_nh3d_timesplit_tsuv) then !w°v°w> !06-09-18
! confirmer le time-splitting (quoiqu'en principe rEglE par dEfaut donc en principe inutile)
       flag_timesplitting=1
       flag_surf_obc_q=1     ! Surface libre hydrostatique entraine condition de surface gradient nul sur q
!      flag_surf_obc_q=0     ! BIDOUILLE PATRICK Bidouille Patrick
      endif                                        !w°v°w>

!     wavebreakfactor=wavebreakfactor*dxb*dyb

!     open(unit=3,file='nhreduce_input')
!      read(3,*)nhpgf_reduce
!     close(3)

      end subroutine q_initial_beforetimestep

!...............................

      subroutine q_update_nhp
      implicit none
!     integer loop_

!     do loop_=1,1000
!     write(6,*)'iteration3d,loop_',iteration3d,loop_

      id_dz=1 ; id_zt=2

      if(nh_scheme==1)call q_moteur_leplussimple
      if(nh_scheme==2)call q_moteur_sigma
      if(nh_scheme==3) then !m[0v0]m>
!       if(flag_nh2d==1)   call q_moteur_remapping
!       if(flag_nh3d==flag_nh3d_nosplit_uv.or. &
!          flag_nh3d==flag_nh3d_nosplit_tsuv)call q_moteur3d_remapping !23-08-18
        if(flag_nh3d/=flag_nh3d_none)call q_moteur3d_remapping !07-09-18
      endif                 !m[0v0]m>
      if(nh_scheme==4)call q_moteur_remapping_b2

! Le cacul de breaker est placE ici car c'est le seul moment du modele
! ou on a les 3 valeurs distinctes de q(0), q(1), q(2)
      call q_trigger_wavebreaking

! La consequence du moveforward immediat est que q(1)=q(2) aussitot...
      call q_moveforward

! Diag, metrix...
      if(flag_nhdiag_timeaverage==1)call q_time_averaged

!     enddo ! loop_
!     if(iteration3d==1)stop 'hello'

      end subroutine q_update_nhp

!...............................

      subroutine q_allocate
      implicit none

! Que ktop soit kmax ou kmax+1 ne change rien au fait qu'il y a kmax vrais niveaux q
!     ktop=kmax   ! q_t est nul en k=kmax
      ktop=kmax+1 ! q_t est nul en k=kmax+1

!     allocate(kvector_oi(imax))   ; kvector_oi=0.
!     allocate(kvector_oj(imax))   ; kvector_oj=0.
!     allocate(phase_wave(imax,2)) ; phase_wave=0.


!     deallocate(omega_w)
!     allocate(omega_w(imax,jmax,kmax+1,0:2)) ; omega_w=0.

      allocate(anyv3dr4(0:imax+1,0:jmax+1,0:kmax+1,2))   ; anyv3dr4=0.
      allocate(      mc(0:imax+1,0:jmax+1,kmax+1,-1:4))  ; mc=0.
      allocate(     q_t(0:imax+1,0:jmax+1,0:kmax+1,-1:2)) ; q_t=0.
      allocate( sshr4_w(1:imax  ,1:jmax  ))              ; sshr4_w=0.
!     allocate(obc_q_i(0:imax+1,kmax,2))                 ; obc_q_i=0.  !11-03-19
!     allocate(obc_ub_i(0:imax+1,kmax,2))                ; obc_ub_i=0.  !11-03-19

!     allocate( invdx_t(1:imax  ,1:jmax))                ; invdx_t=0.
!     allocate( invdx_f(2:imax  ,2:jmax))                ; invdx_f=0.
!     allocate( invdx_u(1:imax+1,1:jmax))                ; invdx_u=0.
!     allocate( invdy_u(1:imax+1,1:jmax))                ; invdy_u=0.
!     allocate( invdy_t(1:imax  ,1:jmax  ))              ; invdy_t=0.
!     allocate( invdy_f(2:imax  ,2:jmax  ))              ; invdy_f=0.
!     allocate( invdy_v(1:imax  ,1:jmax+1))              ; invdy_v=0.
!     allocate( invdx_v(1:imax  ,1:jmax+1))              ; invdx_v=0.

      allocate( nhpgf2d_u(1:imax+1,1:jmax  ))            ; nhpgf2d_u=0.   
      allocate( nhpgf2d_v(1:imax  ,1:jmax+1))            ; nhpgf2d_v=0.   
!     allocate( velbarnh_u(1:imax+1,1:jmax  ))           ; velbarnh_u=0.   
!     allocate( velbarnh_v(1:imax  ,1:jmax+1))           ; velbarnh_v=0.   
      allocate( vqsratio_t(2:imax-1,2:jmax-1))           ; vqsratio_t=1.   
      allocate( nhpgf_u(1:imax+1,1:jmax  ,kmax))         ; nhpgf_u=0.   
      allocate( nhpgf_v(1:imax  ,1:jmax+1,kmax))         ; nhpgf_v=0.   
!     allocate(  dsigr4_t(kmax  ) ) ; dsigr4_t=1./real(kmax)
!     allocate(  dsig_w(kmax+1) ) ; dsig_w=1./real(kmax)
!     allocate(sig1dpom_t(kmax  ) ) ;  sig1dpom_t=0.
      allocate(  dsigr4_t(0:imax+1,0:jmax+1,kmax  ) ) ; dsigr4_t=1./real(kmax)
      allocate(    dsig_w(0:imax+1,0:jmax+1,kmax+1) ) ;   dsig_w=1./real(kmax)
      allocate(sig1dpom_t(0:imax+1,0:jmax+1,kmax  ) ) ;  sig1dpom_t=0.


      allocate(   breaker2d_t(0:imax+1,0:jmax+1))   ; breaker2d_t=0.
      allocate(  sshmin_tmp_w(0:imax+1,0:jmax+1))   ; sshmin_tmp_w=0.
      allocate(      sshmin_w(0:imax+1,0:jmax+1))   ; sshmin_w=0.
      allocate(      sshmax_w(0:imax+1,0:jmax+1))   ; sshmax_w=0.
      allocate(ema1_s_i(imax,2))          ; ema1_s_i=1.d-20
      allocate(ema2_s_i(imax,2))          ; ema2_s_i=1.d-20
! NOte: 2eme argument normalement =2. Si plus c'est qu'on calcule aussi sur des points interieurs pour comparaison
      allocate(ema1_s_j(jmax,4))          ; ema1_s_j=1.d-20
      allocate(ema2_s_j(jmax,4))          ; ema2_s_j=1.d-20
      allocate(ema1_q_j(jmax,kmax,2))          ; ema1_q_j=1.d-20
      allocate(ema2_q_j(jmax,kmax,2))          ; ema2_q_j=1.d-20

!         allocate(sshwave_i_w(0:imax+1     ,3,2,2)) ; sshwave_i_w=0. ! sshwave_i_w(i,pseudo j,temps,frontiere j=1 ou j=jmax)
!         allocate(qwave_i_w  (0:imax+1,kmax,3,2,2)) ; qwave_i_w=0.   ! qwave_i_w(i,k,pseudo j,temps,frontiere j=1 ou j=jmax)
!        allocate(velwave_i_v (0:imax+1,kmax,2,2,2)) ; velwave_i_v=0.
!        allocate(velwave_i_u (0:imax+1,kmax,2,2,2)) ; velwave_i_u=0.
          allocate(sshwave_i_w(0:imax+1     ,3,-1:4,2)) ; sshwave_i_w=0. ! sshwave_i_w(i,pseudo j,temps,frontiere j=1 ou j=jmax)
          allocate(qwave_i_w  (0:imax+1,kmax,3,-1:4,2)) ; qwave_i_w=0.   ! qwave_i_w(i,k,pseudo j,temps,frontiere j=1 ou j=jmax)
         allocate(velwave_i_v (0:imax+1,kmax,2,-1:4,2)) ; velwave_i_v=0.
         allocate(velwave_i_u (0:imax+1,kmax,2,-1:4,2)) ; velwave_i_u=0.

          allocate(sshwave_j_w(0:jmax+1     ,3,2,2)) ; sshwave_j_w=0. ! sshwave_i_w(j,pseudo i,temps,frontiere i=1 ou i=imax)
          allocate(qwave_j_w  (0:jmax+1,kmax,3,2,2)) ; qwave_j_w=0.   ! qwave_i_w(j,k,pseudo i,temps,frontiere i=1 ou i=imax)
         allocate(velwave_j_v (0:jmax+1,kmax,2,2,2)) ; velwave_j_v=0.
         allocate(velwave_j_u (0:jmax+1,kmax,2,2,2)) ; velwave_j_u=0.

!     allocate(     breaker_t(0:imax+1,0:jmax+1      ))   ; breaker_t=0.
!     allocate(   breaker2d_t(0:imax+1,0:jmax+1      ))   ; breaker2d_t=0.
!     allocate(  dt_over_retartedtime(retarded_pss_nb))   ; dt_over_retartedtime=0.
!     allocate(  retarded_pss_amp    (retarded_pss_nb))   ; retarded_pss_amp=0.
!     allocate(retarded_pss_w(imax,jmax,retarded_pss_nb)) ; retarded_pss_w=0.
      allocate(sshlwf_w(imax,jmax,1)) ; sshlwf_w=0.

      if(nh_wavebreakscheme==1) then !ovo>
        allocate(rwavebreak_t(imax,jmax)) ; rwavebreak_t=0.
      endif                          !ovo>

      if(flag_nhdiag_timeaverage==1) then !>>>>>>>>>>>
       allocate(variance_timeaveraged(imax  ,jmax  ))        ; variance_timeaveraged=0.
       allocate(     ssh_timeaveraged(0:imax+1,0:jmax+1))        ; ssh_timeaveraged=0.
       allocate(     flx2d_timeaveraged(imax+1,jmax  ))      ; flx2d_timeaveraged=0. 
       allocate(     flx3d_timeaveraged(imax+1,jmax  ,kmax)) ; flx3d_timeaveraged=0. 
       allocate(   u_euler_timeaveraged(imax+1,jmax  ,kmax)) ; u_euler_timeaveraged=0.
       allocate(     fly2d_timeaveraged(imax  ,jmax+1))      ; fly2d_timeaveraged=0. 
       allocate(     fly3d_timeaveraged(imax  ,jmax+1,kmax)) ; fly3d_timeaveraged=0. 
       allocate(   v_euler_timeaveraged(imax  ,jmax+1,kmax)) ; v_euler_timeaveraged=0.
       allocate(             sshmax_w(imax  ,jmax))          ; sshmax_w=-999.
       allocate(             sshmin_w(imax  ,jmax))          ; sshmin_w=+999.
       allocate(                 hs_w(imax  ,jmax))          ; hs_w=0. ; sum_avr_hs=0
      endif                               !>>>>>>>>>>>

      do k=1,kmax
        sig1dpom_t(:,:,k)=0.5*(-1.+real(k-1)/real(kmax  ) & !23-11-17
                               -1.+real(k  )/real(kmax  ))
      enddo

! Cas particulier de la localisation de la condition de surface (q(kmax+1)=0 en z=ssh)
      dsig_w(:,:,kmax+1)=0.5/real(kmax)

!     do k=kmax+1,1,-1
!      write(6,*)'dsig_w',k,dsig_w(i,j,k)
!     enddo
!     do k=kmax,1,-1
!      write(6,*)'dsigr4_t',k,dsigr4_t(i,j,k),dsig_t(1,1,k)
!     enddo
!     anyv1d(kmax+1,1)=0.+0.5*dsigr4_t(i,j,kmax)
!     do k=kmax,1,-1
!     anyv1d(k,1)=anyv1d(k+1,1)-dsigr4_t(i,j,k)
!      write(6,*)'dsig_t3D',dsig_t(1,1,k),sig1dpom_t(i,j,k)
!      write(6,*)'sig1dpom_t(i,j,k)',real(sig1dpom_t(i,j,k)),real(anyv1d(k,1))
!     enddo

! retrouver dsigr4_t dsig_w sig1dpom_t A partir de dsig_t:
      do j=0,jmax+1 ; do i=0,imax+1
       anyv1d(kmax+1,1)=0 ! sigma "pom" niveau w
       do k=kmax,1,-1
            anyv1d(k,1)=anyv1d(k+1,1)-dsig_t(i,j,k)
          dsigr4_t(i,j,k)=dsig_t(i,j,k)
        sig1dpom_t(i,j,k)=0.5*(anyv1d(k,1)+anyv1d(k+1,1))
!       write(6,*)'sigma pom w t',k,anyv1d(k,1),sig1dpom_t(i,j,k)
       enddo
      enddo       ; enddo

! Ici j'ai verifie que sigma_w pouvait etre utilisE:
!     i=imax/2 ; j=jmax/2
!     do k=kmax,1,-1
!      write(10+par%rank,*)k,dsig_t(i,j,k),sigma_w(i,j,k+1)-sigma_w(i,j,k)
!     enddo
!     stop 'ALBERT'

      do j=0,jmax+1 ; do i=0,imax+1
       do k=2,kmax
        dsig_w(i,j,k)=sig1dpom_t(i,j,k)-sig1dpom_t(i,j,k-1)
!       if(k==kmax/2.and.j==jmax/2)write(10+par%rank,*)i+par%timax(1) &
!         ,dsig_w(i,j,k),0.5*(dsig_t(i,j,k)+dsig_t(i,j,k-1)) &
!         ,' trouve coco'
       enddo
      enddo       ; enddo

      do j=0,jmax+1 ; do i=0,imax+1
! Surface:
       dsig_w(i,j,kmax+1)=0.-sig1dpom_t(i,j,kmax)
! fond: multiplication par 2 pour c.l. de gradient nulle "centree" sur z=-h
       dsig_w(i,j,1)=2.*(sig1dpom_t(i,j,1)-(-1.))
! Note du 26 dec 2021: cette condition ci-dessus n'est peut etre pas bien coherente avec la
! valeur de depth_t en k=0 utilisEe partout et notamment dans le calcul de NHPGF
! Du coup peut etre que l'aternative simple (q(0) correspond a q(z=-h))
! serait finalement preferable:
!      dsig_w(i,j,1)=sig1dpom_t(i,j,1)-(-1.)
      enddo       ; enddo
      if(flag_merged_levels==1) then !m°v°m>
       do j=2,jmax-1 ; do i=2,imax-1
!       if(kmerged_t(i,j)/=kmin_w(i,j)) then !ooo>
!        k=kmin_w(i,j)
!        vqsratio_t(i,j)=min(1.,dsig_t(i,j,k  ) &
!                              /dsig_t(i,j,k+1))
!       else                                 !ooo>
!        vqsratio_t(i,j)=1.
!       endif                                !ooo>
         vqsratio_t(i,j)=0.
       enddo ; enddo
      endif                          !m°v°m>

!     do k=kmax+1,1,-1
!      write(6,*)'dsig_w',k,dsig_w(i,j,k)
!     enddo
!     stop 'coco4'
      
!     do j=1,jmax ; do i=1,imax
!      invdx_t(i,j)=1./dx_t(i,j)
!      invdy_t(i,j)=1./dy_t(i,j)
!     enddo       ; enddo
!     do j=1,jmax ; do i=1,imax+1
!      invdx_u(i,j)=1./dx_u(i,j)
!      invdy_u(i,j)=1./dy_u(i,j)
!     enddo       ; enddo
!     do j=1,jmax+1 ; do i=1,imax
!      invdy_v(i,j)=1./dy_v(i,j)
!      invdx_v(i,j)=1./dx_v(i,j)
!     enddo       ; enddo
!     do j=2,jmax ; do i=2,imax
!      invdx_f(i,j)=1./dx_f(i,j)
!      invdy_f(i,j)=1./dy_f(i,j)
!     enddo       ; enddo

! gradient nul sur q aux bords du domaines
!     if(obcstatus(ieq1)==1)   invdx_u(1,:)=0.
!     if(obcstatus(ieqimax)==1)invdx_u(imax+1,:)=0.
!     if(obcstatus(jeq1)==1)   invdy_v(:,1)=0.
!     if(obcstatus(jeqjmax)==1)invdy_v(:,jmax+1)=0.

      x1_r4=99999.
      do j=1,jmax ; do i=1,imax
       if(mask_t(i,j,kmax)==1) then !m0v0m>

        x1_r4=min(x1_r4,2./( invdx_t(i,j)*invdx_u(i+1,j  )    &
                            +invdx_t(i,j)*invdx_u(i  ,j  )    &
                            +invdy_t(i,j)*invdy_v(i  ,j+1)    &
                            +invdy_t(i,j)*invdy_v(i  ,j  ) ) ) 

       endif                        !m0v0m>
      enddo       ; enddo
      call mpi_allreduce(x1_r4,dt2,1,mpi_real,mpi_min,par%comm2d ,ierr)
      dt2=cfl_nh*dt2
!     write(6,*)'dt2',dt2,0.45*dxb*dxb ; stop 'coco'
      dtconvert=dte_fw/sqrt(dt2)

!      write(6,*)'dte_fw/sqrt(dt2)',sqrt(dt2)/dte_fw,sqrt(grav*hmax)
!      stop 'alors?'

! Initialisation du vent
       uwind_t(:,:,1)=25. !(2*c)
       vwind_t(:,:,1)=0.

! acm_speed est la constante "alpha" des equations (6) et (9) dans Marsaleix et al, 2019 
      acm_speed=sqrt(dt2)/dte_fw !14-02-19

      if(par%rank==0) then
      open(unit=3,file='tmp/messages',position='append')
       write(3,*)'------------------'
       write(3,'(a,a)')'subroutine q_allocate'
! attention sqrt(dt2) n'est pas en secondes!!!
!      write(3,*)'NH ACM phase speed sqrt(dt2)/dte_fw=',sqrt(dt2)/dte_fw
       write(3,*)'NH ACM phase speed acm_speed',acm_speed !14-02-19
       write(3,*)'SQRT(grav*hmax)                    =',sqrt(grav*hmax)
       write(3,*)'Ratio:',sqrt(dt2)/dte_fw/sqrt(grav*hmax)
      close(3)
      endif

      end subroutine q_allocate

!...............................
#ifdef bidon
      subroutine q_equations
      implicit none
     
      if(.not.allocated(q_t))call q_allocate

! Definir dz_t, l'epaisseur au point q:
      x0=1./real(kmax) ; id_dz=1
      do k=0,kmax+1 ; do j=0,jmax+1 ; do i=0,imax+1
       anyv3dr4(i,j,k,id_dz)=x0*(h_w(i,j)+ssh_w(i,j,1))
!      anyv3dr4(i,j,k,id_dz)=1.
      enddo         ; enddo         ; enddo
! Definit z_t la profondeur de q
! Rappel: z_t(k)=z_t(k+1)-0.5*( dz_t(k+1)+dz_t(k) )
! arbitrairement z_t(kmax+1)=surface+0.5*dz_t(kmax+1)
      do j=0,jmax+1 ; do i=0,imax+1
       anyv3dr4(i,j,kmax+1,id_zt)=ssh_w(i,j,1)+0.5*anyv3dr4(i,j,kmax+1,id_dz)
      enddo         ; enddo
      do k=kmax,0,-1 ; do j=0,jmax+1 ; do i=0,imax+1
       anyv3dr4(i,j,k,id_zt)=anyv3dr4(i,j,k+1,id_zt)         &
               -0.5*(  anyv3dr4(i,j,k+1,id_dz)   &
                      +anyv3dr4(i,j,k  ,id_dz))
      enddo          ; enddo         ; enddo


      ww_=w2_ ; if(ww_==0.)ww_=1.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Si diffusion verticale 50% implicite 50% explicite alors diviser par 2:

       mc(i,j,k,1)=ww_*dt2/( (anyv3dr4(i,j,k  ,id_zt)-anyv3dr4(i,j,k-1,id_zt))   &
                            *(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) ) &
             *(-2.)
!      *min(0.,-2.+0.5*(invdx_t(i,j)*(anyv3dr4(i+1,j,k)-anyv3dr4(i-1,j,k,id_zt)))**2  &
!                 +0.5*(invdy_t(i,j)*(anyv3dr4(i,j+1,k)-anyv3dr4(i,j-1,k,id_zt)))**2)

       mc(i,j,k,3)=ww_*dt2/( (anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k  ,id_zt))   &
                            *(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) ) &
             *(-2.)
!      *min(0.,-2.+0.5*(invdx_t(i,j)*(anyv3dr4(i+1,j,k)-anyv3dr4(i-1,j,k,id_zt)))**2  &
!                 +0.5*(invdy_t(i,j)*(anyv3dr4(i,j+1,k)-anyv3dr4(i,j-1,k,id_zt)))**2)

      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       mc(i,j,k,2)=1.-mc(i,j,k,1)-mc(i,j,k,3)+damp     ! (damp=friction)
      enddo
      enddo
      enddo

      k=kmax
      do j=1,jmax
      do i=1,imax
!......
! Cas q_t(kmax+1)=q_t(kmax)
!      mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
!      mc(i,j,k,3)=0.
!......
! Cas q_t(kmax+1)=0
!      mc(i,j,k,3)=0.
! Cas q_t(kmax)=0
       mc(i,j,k,1)=0.
       mc(i,j,k,2)=1.
       mc(i,j,k,3)=0.
       mc(i,j,k,4)=0. !A priori inutile..... car initialisE A zero....

      enddo
      enddo

!.......................................................
! C.L. Fond: q_t(0)=q_t(1)
      k=1
      do j=1,jmax
      do i=1,imax
! C.L. Fond: q_t(0)=q_t(1):
       mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
!......
       mc(i,j,k,1)=0.
      enddo
      enddo
!.......................................................

     
!     i=imax/2 ; j=jmax/2 ; k=kmax/2
!     write(6,*)mc(i,j,k,1:4)
!     return
!     stop 'coucou'
    
!................
!     do loop=1,loopmax ! boucle iterative debut
      loop=1

      x2=real(loop)/real(loopmax) 
!     x2=1.

      x1=1.-x2

      if(loopmax/=1)qavr_t=0. ! reset

      ww_=1.
      if(w2_/=0.)ww_=1./w2_

      call q_moteur_sigma
      i=imax/2 ; j=jmax/2 ; k=kmax/2
      write(6,*)mc(i,j,k,1:4)
      return
      stop 'coucou'


      if(w2_/=0.) then !m0v0m>
       if(loop==1) then
        call q_solver(1,1,imax,1,jmax,1,kmax) ! donne q_t(2)
       else
        call q_solver(0,1,imax,1,jmax,1,kmax) ! donne q_t(2)
       endif
      else             !m0v0m>
       do k=1,kmax-1 ; do j=1,jmax ; do i=1,imax
        q_t(i,j,k,2)=mc(i,j,k,4)
       enddo         ; enddo       ; enddo
      endif            !m0v0m>

! En k=0 reporter la C.L sur q_t(t+1)
      do j=1,jmax ; do i=1,imax
! C.L.: q(0)=q(1):
       q_t(i,j,0,2)=q_t(i,j,1,2)*mask_t(i,j,kmax)
! C.L.: q(kmax+1)=q(kmax):
!      q_t(t+1)(i,j,kmax+1)=q_t(t+1)(i,j,kmax)*mask_t(i,j,kmax)
! C.L.: q(kmax+1)=0.
       q_t(i,j,kmax:kmax+1,2)=0.
      enddo ; enddo

! ZONE EPONGE: (Ne pourrait elle pas etre implicitee?)
      do k=1,kmax-1 ; do j=1,jmax ; do i=1,imax
       q_t(i,j,k,2)=q_t(i,j,k,2)*(1.-sponge_t(i,j,1))
      enddo         ; enddo       ; enddo

      call q_obc_mpi('z1') ! mpi 'z1' sur q_t(t+1)

!....................................
! PAS OBLIGATOIRE, POUR VERIFICATION:
      if(iteration2d==iteration2d_max_now)call q_solver_balance
!....................................

      if(loopmax/=1) then !>>>
       do k=1,kmax ; do j=1,jmax  ; do i=1,imax  
         qavr_t(i,j,k)=  q_t(i,j,k,2)+qavr_t(i,j,k)
       enddo       ; enddo        ; enddo
      endif               !>>>

!     enddo ! boucle iterative fin
!............

      if(loopmax/=1) then !>>>
       do k=1,kmax ; do j=1,jmax  ; do i=1,imax  
         q_t(i,j,k,2)=qavr_t(i,j,k)
       enddo       ; enddo        ; enddo
      endif               !>>>

! ATTENTION CETTE MOYENNE SUPPOSE UNE DISTRIBUTION REGULIERE DES NIVEAUX:
      if(par%rank==0.and.mod(iteration2d,100)==0) &
      write(6,*)'bidouille moyenne verticale suppose dz constant'

      do j=1,jmax ; do i=1,imax  
       q2davr_w(i,j)=0.
      enddo       ; enddo
      do k=kmax,1,-1 ! la boucle commence A kmax+1 mais comme q(kmax+1)=0 alors....
      do j=1,jmax  
      do i=1,imax  
       q2davr_w(i,j)=q2davr_w(i,j)+q_t(i,j,k,2)
      enddo
      enddo
      enddo

      x0=1./real(kmax*loopmax)
! Pour compenser le fait que les C.L du fond et de surface sont appliquees A une
! demi-maille de distance on applique un facteur de correction: 1+1/kmax
!     x0=x0*(1.+1./real(kmax))
      do j=1,jmax  
      do i=1,imax  
       q2davr_w(i,j)=x0*q2davr_w(i,j)
      enddo
      enddo

!     call q_poisson_local_1dv_diag
!     call q_pseudo3d_A

      end subroutine q_equations
#endif
!............................

      subroutine q_obc_mpi(txt_)
      implicit none
      character*2 txt_
      integer(kind=1) loop_

#ifdef parallele
!!$! Nouvelle methode avec choix des voisins
!      call get_type_echange(txt_,'q_t(t+1)_'//txt_,q_t(t+1),lbound(q_t(t+1)),ubound(q_t(t+1)),k10)
       call get_type_echange(txt_,'q_t_'//txt_,q_t,lbound(q_t),ubound(q_t),2,k10)

      ! Echanges
      do loop_=1, subcycle_exchange
         call echange_voisin(q_t,k10,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine q_obc_mpi

!............................

      subroutine q_solver(case_,i1_,i2_,j1_,j2_,k1_,k2_)
      implicit none
      integer case_,i1_,i2_,j1_,j2_,k1_,k2_

      k=k1_
      if(case_/=0) then !------->
       do j=j1_,j2_ ; do i=i1_,i2_
        mc(i,j,k,0)=1./mc(i,j,k,2)
       enddo ; enddo
      endif             !------->
      do j=j1_,j2_ ; do i=i1_,i2_
       mc(i,j,k,-1)=mc(i,j,k,4)               &
                   *mc(i,j,k,0)
      enddo ; enddo


      if(case_/=0) then !------->
       do k=k1_+1,k2_ ; do j=j1_,j2_  ; do i=i1_,i2_

        mc(i,j,k,0)=1./(mc(i,j,k  ,2)      &
                       -mc(i,j,k  ,1)      &
                       *mc(i,j,k-1,3)      &
                       *mc(i,j,k-1,0))

       enddo ; enddo ; enddo
      endif             !------->

      do k=k1_+1,k2_ ; do j=j1_,j2_  ; do i=i1_,i2_

        mc(i,j,k,-1)=( mc(i,j,k  ,4)           &
                      -mc(i,j,k  ,1)           &
                      *mc(i,j,k-1,-1) )        &
                      *mc(i,j,k  ,0)

      enddo ; enddo ; enddo


      do j=j1_,j2_ ; do i=i1_,i2_
       q_t(i,j,k2_,2)=mc(i,j,k2_,-1)
      enddo ; enddo

      do k=k2_-1,k1_,-1
      do j=j1_,j2_
      do i=i1_,i2_

        q_t(i,j,k,2)= mc(i,j,k,-1)                       &
                      -mc(i,j,k,3)                        &
                  *q_t(i,j,k+1,2)                        &
                      *mc(i,j,k,0)
      enddo
      enddo
      enddo

      end subroutine q_solver
      
!............................

      subroutine q_solver_balance
      implicit none

      ww_=1.
      if(w2_/=0.)ww_=1./w2_


!     i=imax/2 ; j=jmax/2 ; k=kmax/2
      i=10 ; j=jmax/2 ; k=kmax/2

      if(iteration2d==iteration2d_max_now) then !>>>>
      do k=ktop-1,1,-1

      write(10+par%rank,'(6(1x,e14.7))') &
                    q_t(i,j,k  ,2)-2*q_t(i,j,k,1)+q_t(i,j,k,0)  &

            ,   +anyvar2d(i,j)   &
                                              +dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
                ( q_t(i+1,j,k,1)-q_t(i  ,j,k,1) )*invdx_u(i+1,j)    &
               -( q_t(i  ,j,k,1)-q_t(i-1,j,k,1) )*invdx_u(i  ,j)    &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
                ( q_t(i,j+1,k,1)-q_t(i  ,j,k,1) )*invdy_v(i,j+1)    &
               -( q_t(i,j  ,k,1)-q_t(i,j-1,k,1) )*invdy_v(i  ,j)    &
                    ) & !ooo>

                                                   ) &  !ttt>


       -mc(i,j,k,1)   *ww_*(w0_*q_t(i,j,k-1,0)+w1_*q_t(i,j,k-1,1)+w2_*q_t(i,j,k-1,2)) &
       -mc(i,j,k,3)   *ww_*(w0_*q_t(i,j,k+1,0)+w1_*q_t(i,j,k+1,1)+w2_*q_t(i,j,k+1,2)) &
      -(mc(i,j,k,2)-1)*ww_*(w0_*q_t(i,j,k  ,0)+w1_*q_t(i,j,k  ,1)+w2_*q_t(i,j,k  ,2))

! Note 1: la mutiplication par ww_ correspond A une division par ww2_ qui s'explique par le fait
! que les coef mc contiennent ww2_ dont il faut donc les debarasser

! Note 2: on soustrait 1 A mc(2) pour ramener mc(2) A sa partie purement melange verticale car le 1
! correspond A l'operateur de time-stepping

! La commande gnuplot pour verifier
! plot "fort.x" u 0:1 title "unbalanced term" w l, "fort.x" u 0:2 title "2DH Laplacian Hydrostatic term" w l, "fort.x" u 0:3 title "3D Laplacian Non-Hydrostatic Pressure" w l
      enddo
      endif                                     !>>>>

      end subroutine q_solver_balance

!............................

      subroutine q_pseudo3d_moveforward
      implicit none
! AVANCER En AVANT:
       do k=1,kmax
       do j=2,jmax-1
       do i=2,imax
         vel_u(i,j,k,0)=vel_u(i,j,k,1)
         vel_u(i,j,k,1)=vel_u(i,j,k,2)
       enddo
       enddo
       enddo
       do k=1,kmax
       do j=2,jmax
       do i=2,imax-1
         vel_v(i,j,k,0)=vel_v(i,j,k,1)
         vel_v(i,j,k,1)=vel_v(i,j,k,2)
       enddo
       enddo
       enddo
       do j=1,jmax
       do i=1,imax
       do k=2,kmax+1
         omega_w(i,j,k,0)=omega_w(i,j,k,1)
         omega_w(i,j,k,1)=omega_w(i,j,k,2)
       enddo
       enddo
       enddo

      end subroutine q_pseudo3d_moveforward

!............................
#ifdef bidon
      subroutine q_pseudo3d
      implicit none

! Definir dz_t, l'epaisseur au point q:
      x0=1./real(kmax) ; id_dz=1
      do k=0,kmax+1 ; do j=0,jmax+1 ; do i=0,imax+1
       anyv3dr4(i,j,k,id_dz)=x0*(h_w(i,j)+ssh_w(i,j,1))
!      anyv3dr4(i,j,k,id_dz)=1.
      enddo         ; enddo         ; enddo
! Definit z_t la profondeur de q
! Rappel: z_t(k)=z_t(k+1)-0.5*( dz_t(k+1)+dz_t(k) )
! arbitrairement z_t(kmax+1)=surface+0.5*dz_t(kmax+1)
      do j=0,jmax+1 ; do i=0,imax+1
       anyv3dr4(i,j,kmax+1,id_zt)=ssh_w(i,j,1)+0.5*anyv3dr4(i,j,kmax+1,id_dz)
      enddo         ; enddo
      do k=kmax,0,-1 ; do j=0,jmax+1 ; do i=0,imax+1
       anyv3dr4(i,j,k,id_zt)=anyv3dr4(i,j,k+1,id_zt)         &
               -0.5*(  anyv3dr4(i,j,k+1,id_dz)   &
                      +anyv3dr4(i,j,k  ,id_dz))
      enddo          ; enddo         ; enddo

! COMPOSANTE U
       x0=1./real(kmax)
       sum1=0.
       sum2=0.
       sum6=0.
       sum3=0.
       sum5=0.
       do k=1,kmax
       do j=2,jmax-1
       do i=2,imax

!        x1=x0*(h_u(i,j)+0.5*(ssh_w(i,j,0)+ssh_w(i-1,j,0)))
!        x2=x0*(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1)))
!        x1=x0*h_u(i,j)
!        x2=x0*h_u(i,j)
         x1=0.5*(anyv3dr4(i-1,j,k,id_dz)+anyv3dr4(i,j,k,id_dz))
         x2=0.5*(anyv3dr4(i-1,j,k,id_dz)+anyv3dr4(i,j,k,id_dz))

! x3 pgf de ssh
         x3=x2*grav*(ssh_w(i  ,j,1)-ssh_w(i-1,j,1))/dx_u(i,j)

         if(i==imax/4.and.j==jmax/2.and.k==kmax/2) then
          write(69,*)x2*(q_t(i,j,k,1)   -q_t(i-1,j,k,1))/dx_u(i,j) &
        ,-(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i-1,j,k,id_zt)  )/dx_u(i,j)   &
        *0.25*( q_t(i,j,k+1,1) +q_t(i-1,j,k+1,1)            &
               -q_t(i,j,k-1,1) -q_t(i-1,j,k-1,1))
         endif

         vel_u(i,j,k,2)=mask_u(i,j,k)*(vel_u(i,j,k,1)*x1 &

           +dte_fw*( & !ooo>

              -x3-nhpgf_u(i,j,k) &
!             -x3-anyv3dr4(i,j,k,id_dzpgfu) &

                   ) & !ooo>

                         )/x2

        sum1=sum1+0.5*(vel_u(i,j,k,2)+vel_u(i,j,k,1))*x3 ! bilan ssh
        sum6=sum6+0.5*(vel_u(i,j,k,2)+vel_u(i,j,k,1))*nhpgf_u(i,j,k) ! anyv3dr4(i,j,k,id_dzpgfu) ! bilan q


        sum2=sum2+0.5*x2*(vel_u(i,j,k,2)+vel_u(i,j,k,1))*(vel_u(i,j,k,2)-vel_u(i,j,k,1))/dte_fw
        sum3=sum3+0.5*x2*(vel_u(i,j,k,1))*(vel_u(i,j,k,1))
        sum5=sum5+0.5*x2*(vel_u(i,j,k,2))*(vel_u(i,j,k,2))

       enddo
       enddo
       enddo

       stop 'caravane'


       sum7=0.
       do k=1,kmax
       do j=1,jmax
       do i=1,imax
        sum7=sum7-mask_t(i,j,k)*q_t(i,j,k,1)/dx_u(i,j) &
                 *0.5*( x0*h_u(i+1,j)*vel_u(i+1,j,k,2)-x0*h_u(i,j)*vel_u(i,j,k,2) &
                       +x0*h_u(i+1,j)*vel_u(i+1,j,k,1)-x0*h_u(i,j)*vel_u(i,j,k,1) )
       enddo
       enddo
       enddo

! COMPOSANTE V
!      sum2=0.
       do k=1,kmax
       do j=2,jmax
       do i=2,imax-1


         x1=x0*(h_v(i,j)+0.5*(ssh_w(i,j,0)+ssh_w(i,j-1,0)))
         x2=x0*(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1)))
!        x1=x0*h_v(i,j)
!        x2=x0*h_v(i,j)

         vel_v(i,j,k,2)=mask_v(i,j,k)*(vel_v(i,j,k,1)*x1 &

           +dte_fw*( & !ooo>

               -x2*(grav*(ssh_w(i  ,j,1)-ssh_w(i,j-1,1))            & 
                           +q_t(i  ,j,k,1)-q_t(i,j-1,k,1)  )/dy_v(i,j) &

                   ) & !ooo>

                         )/x2

       sum1=sum1+0.5*(vel_v(i,j,k,2)+vel_v(i,j,k,1)) &
                    *x2*grav*(ssh_w(i  ,j,1)-ssh_w(i,j-1,1))/dy_v(i,j)
       sum6=sum6+0.5*(vel_v(i,j,k,2)+vel_v(i,j,k,1)) &
                         *x2*(q_t(i  ,j,k,1)-q_t(i,j-1,k,1))/dy_v(i,j)

        sum2=sum2+0.5*x2*(vel_v(i,j,k,2)+vel_v(i,j,k,1))*(vel_v(i,j,k,2)-vel_v(i,j,k,1))/dte_fw

        sum3=sum3+0.5*x2*(vel_v(i,j,k,1))*(vel_v(i,j,k,1))
        sum5=sum5+0.5*x2*(vel_v(i,j,k,2))*(vel_v(i,j,k,2))

       enddo
       enddo
       enddo
       sum_t_udssh=sum_t_udssh+dte_fw*sum1

       do k=1,kmax
       do j=1,jmax
       do i=1,imax
        sum7=sum7-mask_t(i,j,k)*q_t(i,j,k,1)/dy_v(i,j) &
                 *0.5*( x0*h_v(i,j+1)*vel_v(i,j+1,k,2)-x0*h_v(i,j)*vel_v(i,j,k,2) &
                       +x0*h_v(i,j+1)*vel_v(i,j+1,k,1)-x0*h_v(i,j)*vel_v(i,j,k,1) )
       enddo
       enddo
       enddo
       write(6,*)'riri',sum6,sum7

! COMPOSANTE W
       sum4=0.
       do j=1,jmax
       do i=1,imax
       do k=2,kmax+1

         omega_w(i,j,k,2)=mask_t(i,j,k)*(omega_w(i,j,k,1)    &

           +dte_fw*( & !ooo>

            ( & 
               -w0_*(    q_t(i,j,k,0)   -q_t(i,j,k-1,0) ) &
               -w1_*(    q_t(i,j,k,1)   -q_t(i,j,k-1,1) ) &
               -w2_*( q_t(i,j,k,2)  -q_t(i,j,k-1,2)   ) &

            ) /( anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)   ) &

                   ) & !ooo>
                                        )
        
        sum4=sum4+0.5*(omega_w(i,j,k,2)+omega_w(i,j,k,1)) &
!                        *(q_t(i,j,k,1)   -q_t(i,j,k-1,1))
!                    *( q_t(i,j,k,2)  -q_t(i,j,k-1,2)  )  
              *(w0_*(    q_t(i,j,k,0)   -q_t(i,j,k-1,0) ) &
               +w1_*(    q_t(i,j,k,1)   -q_t(i,j,k-1,1) ) &
               +w2_*( q_t(i,j,k,2)  -q_t(i,j,k-1,2)   ) )

        sum2=sum2+0.5*(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)) &
                     *(omega_w(i,j,k,2)+omega_w(i,j,k,1)) &
                     *(omega_w(i,j,k,2)-omega_w(i,j,k,1))/dte_fw

        sum3=sum3+0.5*(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)) &
                     *(omega_w(i,j,k,1)) &
                     *(omega_w(i,j,k,1))

        sum5=sum5+0.5*(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)) &
                     *(omega_w(i,j,k,2)) &
                     *(omega_w(i,j,k,2))

       enddo
       enddo
       enddo
       sum_t_udq_h=sum_t_udq_h+dte_fw*sum6 ! sum u.dq/dx
       sum_t_udq_v2=sum_t_udq_v2+dte_fw*sum4 ! sum w.dq/dz effectif
       sum_t_udu=sum_t_udu+dte_fw*sum2     ! sum u.du/dt+w.dw/dt

       if(iteration3d==0.and.iteration2d==1)sum_u2_init=sum3
       sum_u2_now=sum5

       sum8=0.
       sum9=0.
       do j=1,jmax
       do i=1,imax
       do k=1,kmax
       sum8=sum8-q_t(i,j,k,1)*mask_t(i,j,k)*0.5*( &
        omega_w(i,j,k+1,2)-omega_w(i,j,k,2)        &
       +omega_w(i,j,k+1,1)-omega_w(i,j,k,1)      )
!      sum9=sum9-(q_t(i,j,k,2)-q_t(i,j,k,1))*mask_t(i,j,k)*0.5*( &
       sum9=sum9-(w2_*q_t(i,j,k,2)+w1_*q_t(i,j,k,1)+w0_*q_t(i,j,k,0)-q_t(i,j,k,1))*mask_t(i,j,k)*0.5*( &
        omega_w(i,j,k+1,2)-omega_w(i,j,k,2)        &
       +omega_w(i,j,k+1,1)-omega_w(i,j,k,1)      )
       enddo
       enddo
       enddo
       sum_t_udq_err=sum_t_udq_err+dte_fw*sum9 ! sum w.dq/dz erreur effectif-ideal
       sum_t_udq_v1=sum_t_udq_v1+dte_fw*sum8   ! sum w.dq/dz ideal
!      write(6,*)'loulou',real(sum4),real(sum8)+real(sum9)
       write(6,*)'loulou',real(sum_t_udq_v2),real(sum_t_udq_v1)+real(sum_t_udq_err)

       sum3=0.
       do j=1,jmax
       do i=1,imax

         sum3=sum3+mask_t(i,j,kmax)*ssh_w(i,j,1)*grav*0.5*( &
                     +(ssh_w(i,j,2)-ssh_w(i,j,0))/dte_fw    &
                                                          )

       enddo
       enddo
       sum_t_sshdssh=sum_t_sshdssh+dte_fw*sum3
       
       if(iteration3d==0.and.iteration2d==1) then !ooo>
        sum_ssh2_init=0.
        do j=1,jmax ; do i=1,imax
         sum_ssh2_init=sum_ssh2_init+0.5*grav*ssh_w(i,j,0)*ssh_w(i,j,1)*mask_t(i,j,kmax)
        enddo       ; enddo
       endif                                      !ooo>

       sum_ssh2_now=0.
       do j=1,jmax ; do i=1,imax
        sum_ssh2_now=sum_ssh2_now+0.5*grav*ssh_w(i,j,2)*ssh_w(i,j,1)*mask_t(i,j,kmax)
       enddo       ; enddo

       sum5=0.
       sum0=0.
       sum_q2_now=0.
       do j=1,jmax
       do i=1,imax
       do k=1,kmax
         sum5=sum5+mask_t(i,j,k)*0.5*q_t(i,j,k,1) &
      *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
                     *(q_t(i,j,k,2)-q_t(i,j,k,0))*(dtconvert*dtconvert)/dte_fw
!     sum1=sum1+0.5*(q_t(i,j,k,2)+q_t(i,j,k,1)) &
!     *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
!         *(q_t(i,j,k,2)-q_t(i,j,k,1))*(dtconvert*dtconvert)/dte_fw

        sum0=sum0+mask_t(i,j,k)*0.5*q_t(i,j,k,1) &
      *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
                     *(q_t(i,j,k,0))*(dtconvert*dtconvert)
        sum_q2_now=sum_q2_now+mask_t(i,j,k)*0.5*q_t(i,j,k,1) &
      *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
                     *(q_t(i,j,k,2))*(dtconvert*dtconvert)
       enddo
       enddo
       enddo
       sum_t_qdq=sum_t_qdq+sum5*dte_fw
       if(iteration3d==0.and.iteration2d==1)sum_q2_init=sum0

       write(6,*)'bobi',sum5,sum7+sum8
!      write(6,*)'bobi',real(sum_t_qdq),real(sum_t_udq_h+sum_t_udq_v2)
!      write(6,*)'bobi',real(sum_t_udq_v2),real(sum_t_udq_v1)+real(sum_t_udq_err)

!      write(66,*)sum1,sum3
       write(66,*)real(sum_t_udssh),real(sum_t_sshdssh),real(sum_ssh2_now-sum_ssh2_init)

!      write(67,*)sum2,sum1+sum6+sum4
!      write(67,*)real(sum_t_udu),real(sum_t_udssh+sum_t_udq_h)
!      write(67,*)real(sum_t_udu+sum_t_sshdssh),-real(sum_t_udq_h)-real(sum_t_udq_v2)
!      write(67,*)real(sum_t_udu+sum_t_sshdssh) &
!               ,-real(sum_t_udq_h+sum_t_udq_v1)-real(sum_t_udq_err)
       write(67,'(9(e14.7,1x))')            &
                ,real(sum_t_udu)            & ! C1
                ,real(sum_t_sshdssh)        & ! C2
                ,real(sum_t_qdq)            & ! C3
                ,-real(sum_t_udq_err)       & ! C4=C1+C2+C4 devrait etre zero idealement
                ,sum_u2_now-sum_u2_init     & ! C5=C1
                ,sum_ssh2_now-sum_ssh2_init & ! C6=C2
                ,sum_q2_now-sum_q2_init       ! C7=C3


      end subroutine q_pseudo3d
#endif
!............................
!#ifdef bidon
! REFERENCE: Kinsman Blair, 1984. WIND WAVES Their Generation and Propagation on the Ocean Surface. 
! Dover Publications, inc., New York. 676 pp. 
! Note: la premiere edition date de 1965
! La page 133 est plus particuliererment concernee par cette subroutine
      subroutine q_poisson_ondes_diag
      implicit none

      stop 'passe pas par q_poisson_ondes_diag'

#ifdef bidon
!     x0=real(iglb-2)*dxb/86. ! longueur d'onde
      x0=real(iglb-2)*dxb/20. ! longueur d'onde
      kvector_=2.*pi/x0
      ssh_amplitude_=0.001
      if(period_==-999.) then !aaa>

! Cas NH
       c_=sqrt( grav/kvector_ *tanh(kvector_*hmax) )
! Cas Hydro
!      c_=sqrt(grav*hmax)

       period_=2.*pi/(kvector_*c_)

       if(par%rank==0) then
       write(6,*)'c_=',c_
       write(6,*)'period_',period_
       write(6,*)'longueur d''onde',x0
       write(6,*)'ssh_amplitude_',ssh_amplitude_
       endif

      x1=2.*pi*iteration2d*dte_fw/period_+0.5*pi ! (+0.5*pi pour courant nul au depart)
      x2=kvector_*dxb
      do j=0,jmax+1 ; do i=0,imax+1
       sshobc_w(i,j,:)=ssh_amplitude_*cos(x2*real(i+par%timax(1)-1)-x1)
       do k=1,kmax
         q_t(i,j,k,0:1)=                                               &
        grav*sshobc_w(i,j,1)*(cosh(kvector_*(h_w(i,j)+depth_t(i,j,k))) &
                             /cosh(kvector_*(h_w(i,j)               )) &
                             -1.)
        enddo
      enddo         ; enddo
      do j=0,jmax+1 ; do i=1,imax+1
       velbarobc_u(i,j,1)=0.
       do k=1,kmax
        velobc_u(i,j,k,1)=( &  !ooo>
           0.5*(sshobc_w(i-1,j,1)+sshobc_w(i,j,1))*grav &
          +0.5*(     q_t(i-1,j,k,1)   +q_t(i,j,k,1))    &
                          )/c_ !ooo>

        velobc_u(i,j,k,:)=velobc_u(i,j,k,1)

        velbarobc_u(i,j,1)=velbarobc_u(i,j,1)   &
                                 +dz_u(i,j,k,1) &
                             *velobc_u(i,j,k,1)
       enddo
       velbarobc_u(i,j,1)=velbarobc_u(i,j,1) &
                                 /h_u(i,j)
       velbarobc_u(i,j,:)=velbarobc_u(i,j,1)
!      velbarobc_u(i,j,:)=0.5*(sshobc_w(i-1,j,1)+sshobc_w(i,j,1))*grav/sqrt(h_u(i,j)*grav)
!      velbarobc_u(i,j,:)=ssh_amplitude_*cos(x2*( real(i+par%timax(1)-1)-0.5 )-x1)*grav/sqrt(h_u(i,j)*grav)
      enddo         ; enddo
      velbarobc_v=0.
         velobc_v=0.

!       do k=1,100
!        write(6,*)'VEL Q INIT=0'
!       enddo
!       q_t=0.
!       velbarobc_u=0.
!       velobc_u=0.

      endif                   !aaa>
      return
#endif

! Si kvector n'a pas encore ete calculE alors calculer kvector_:
      if(kvector_==-999.) then !oooo>

!     period_=2.5 
      period_=1.  
      pulsation_=2.*pi/period_

      grav_over_2pi_=grav/(2.*pi)

!     i=imax/2 ; j=jmax/2 ; k=kmax/2
      i=1 ; j=jmax/2 ; k=kmax/2

! c_ is the phase speed
! On veut résoudre x=f(x)=alpha*tanh(beta/x) la vitesse de phase c_  etant l'inconnue x
! Le raisonnement est le suivant: f(x+dx)=x+dx
!                           soit: f(x)+dx f'(x)=x+dx
!                           soit: dx=(x-f(x))/(f'(x)-1)
! On estime une première valeur de x=alpha_ l'approximation pour c_ aux grandes profondeurs
      beta_=2.*pi*max(h_w(i,j),0.01d0)/period_
!     beta_=2.*pi*max(hmax,0.01d0)/period_
      alpha_=grav_over_2pi_*period_                ! c=g*t/2pi first guess for large depth

      if(beta_/alpha_>200.) then !fgfgfgfgfgfgfg>      !07-05-11

! Cas c=first guess valide:
       c_=alpha_ ! phase speed for large depth

      else                       !fgfgfgfgfgfgfg>

! Cas calcul itératif:
       delta_c_=0.01                  ! dx, on estime une première valeur de dx
       c_=alpha_+delta_c_             ! alpha_+dx
         do iter_loop_=1,200

          delta_c_=(c_-alpha_*tanh(beta_/c_))                & !  (x-f(x))
                  /(-alpha_*beta_/(c_*cosh(beta_/c_))**2-1.)   !  /(f'(x)-1)
          c_=max(c_+delta_c_,1.d-4)
          if(abs(delta_c_)<0.01)goto 25
         enddo
         write(6,*)c_,delta_c_
         stop 'waveforcing: le calcul de kn ne converge pas'
   25  continue

      endif                      !fgfgfgfgfgfgfg>

      kvector_=2.*pi/(c_*period_)

!     if(par%rank==0) then
!      write(6,*)'period_=   ',period_
!      write(6,*)'pulsation_=',pulsation_
!      write(6,*)'sqrt(gH) & phase speed:',sqrt(grav*h_w(i,j)),c_
!      write(6,*)'kvector:',kvector_
!      write(6,*)'longueur d onde:',c_*period_
!     endif
!     stop 'q_poisson_ondes_diag'

      endif                    !oooo>


!     ssh_amplitude_=0.01   ! amplitude de l'onde de surface
!     ssh_amplitude_=0.091  ! amplitude de l'onde de surface
!     ssh_amplitude_=0.08    ! amplitude de l'onde de surface
      ssh_amplitude_=0.09    ! amplitude de l'onde de surface
!     ssh_amplitude_=0.001  ! amplitude de l'onde de surface
      x1=2.*pi*iteration2d*dte_fw/period_+0.5*pi ! (+0.5*pi pour courant nul au depart)
      x2=kvector_*dxb
      do j=0,jmax+1 ; do i=0,imax+1
       sshobc_w(i,j,1)=ssh_amplitude_*cos(x2*real(i+par%timax(1)-1)-x1)
!      sshobc_w(i,j,1)=ssh_amplitude_ !*cos(x2*real(i+par%timax(1)-1)-x1)
      enddo         ; enddo
      do j=0,jmax+1 ; do i=0,imax+1
       sshobc_w(i,j,:)=sshobc_w(i,j,1)
      enddo         ; enddo
      do j=0,jmax+1 ; do i=1,imax+1
       velbarobc_u(i,j,:)=0.5*(sshobc_w(i-1,j,1)+sshobc_w(i,j,1))*grav/sqrt(h_u(i,j)*grav)
      enddo         ; enddo

!     if(par%rank==0) then
!      write(6,*)'ssh amplitude         :',ssh_amplitude_
!      write(6,*)'Stokes Transport      :',kvector_*(ssh_amplitude_**2)*c_/2.
!     endif

! Fonction de courant du point horizontal "u" et vertical "w"
      do k=1,kmax+1
      do j=1,jmax ; do i=1,imax+1
      x3=(real(k-1)/real(kmax) - 1.)*h_w(i,j) ! z
       anyv3d(i,j,k,1)=-c_*sinh(kvector_*(x3+h_w(i,j))) &
                          /sinh(kvector_*(  +h_w(i,j))) &
                *0.5*(sshobc_w(i-1,j,1)+sshobc_w(i,j,1))
      enddo       ; enddo
      enddo
! Composante u du courant:
      do j=1,jmax ; do i=1,imax+1
       velbarobc_u(i,j,1)=0.
      enddo       ; enddo
      do k=1,kmax
      do j=1,jmax ; do i=1,imax+1
       velobc_u(i,j,k,1)=-(anyv3d(i,j,k+1,1)-anyv3d(i,j,k,1)) & ! -dphi/dk
                  *real(kmax)/h_w(i,j)                          !      /dz
       velbarobc_u(i,j,1)=velbarobc_u(i,j,1)+velobc_u(i,j,k,1)
      enddo       ; enddo
      enddo
      do j=1,jmax ; do i=1,imax+1
       velbarobc_u(i,j,1)=velbarobc_u(i,j,1)/real(kmax)
      enddo       ; enddo
      do j=1,jmax ; do i=1,imax+1
       velbarobc_u(i,j,:)=velbarobc_u(i,j,1)
       do k=1,kmax
        velobc_u(i,j,k,:)=velobc_u(i,j,k,1)
       enddo
      enddo       ; enddo

! Vitesse verticale:
      do k=1,kmax+1
      do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,2)=(anyv3d(i+1,j,k,1)-anyv3d(i,j,k,1))/dxb ! dphi/dx
      enddo       ; enddo
      enddo

      u_amplitude_=velbarobc_u(1,jmax/2,1)

! Ici on modifie la vitesse de phase prise en compte dans obc_ext:
! Pour la C.L. en i=1
      do j=1,jmax
       sqr_hoverg_u(j,1)=h_u(2,j)/max(c_,1.)
      enddo

      if(par%rank==0) then
       i=1 ; j=jmax/2
       write(6,*)'--------------------------'
       write(6,*)'i,j',i,j
       write(6,*)'period_=   ',period_
       write(6,*)'pulsation_=',pulsation_
       write(6,*)'sqrt(gH) & phase speed:',sqrt(grav*h_w(i,j)),c_
       write(6,*)'kvector:',kvector_
       write(6,*)'longueur d onde:',c_*period_
       write(6,*)'ssh amplitude         :',ssh_amplitude_
       write(6,*)'Stokes Transport      :',kvector_*(ssh_amplitude_**2)*c_/2.
       write(6,*)'u_amplitude_=',u_amplitude_  
       write(6,*)'sqrt(gH) bis',h_u(2,j)/sqr_hoverg_u(j,1),sqrt(grav*h_u(2   ,j))
       write(6,*)'test',real(sqr_hoverg_u(j,1)),real(h_u(2,j)/sqrt(grav*h_u(2   ,j))),real(h_u(2,j)/c_)
!      write(6,*)'Angle (si obc_ext modifie direction):',rad2deg*(asin(2.*pi/kvector_ /(dyb*(jglb-2))))
       write(6,*)'--------------------------'
      endif

! UNE FOIS LES AMPLITUDES ET AUTRES PARAMETRES TROUVES ON REMET LES TABLEAUX OBC A ZERO
! POUR REALISER LA SIMULATION FORCEE PAR LA FRONTIERE OUVERTE:
       sshobc_w=0.
       velobc_u=0.
       velbarobc_u=0.
 
!     write(6,*)'BIDOUILLE INITIALISATION SSH'
!     do j=0,jmax+1
!     do i=0,imax+1
!      i1=i+par%timax(1)
!      j1=j+par%tjmax(1)
!      if(i1>iglb-50)sshobc_w(i,j,:)=-h_w(i,j)
!     enddo
!     enddo
      

#ifdef bidon
      x1=2.*pi*iteration2d*dte_fw/period_+0.5*pi ! (+0.5*pi pour courant nul au depart)
      x2=kvector_*dxb
      do j=1,jmax 
      do i=0,0
       sshobc_w(i,j,1)=ssh_amplitude_*cos(x2*real(i+par%timax(1)-1)-x1)

!      anyv3d(i,j,kmax,0)=0. ! Pression NH k_t=kmax
       q_t(i,j,kmax,1)=0. ! Pression NH k_t=kmax
       do k=kmax,2,-1
        q_t(i,j,k-1,1)=q_t(i,j,k,1)-h_w(i,j)/real(kmax) &
      *pulsation_**2*sinh(kvector_*(depth_w(i,j,k)+h_w(i,j))) &
                    /sinh(kvector_*(              +h_w(i,j)))*sshobc_w(i,j,1)
       enddo

      enddo       
      enddo
#endif

      end subroutine q_poisson_ondes_diag
!#endif
!............................

      subroutine q_moveforward
      implicit none

! FILTRE D'ASSELIN coef constant
        x1_r4=1.-asselin_nh
        x2_r4=0.5*asselin_nh

        do k=0,kmax+1
        do j=0,jmax+1
        do i=0,imax+1

#ifdef bidon
!.............
! AUGMENTER LE COEF D'ASSELIN SI INSTABILITE SUR Q
! x0_r4 est le coef d'asselin instantannE
      x0_r4=        &
! l'indice de "variation" varie entre 0 et 1 quand d2q/dt2=brkvar
! Quand ce dernier=1, alors le coef d'asselin = 1*(0.5-asselin_nh)+asselin_nh = 0.5 c.a.d. la valeur max du coef d'asselin
            min(1., & !MIN>
                  ( ( ( &
                  0.5*abs( & !ooo>
!            ssh_int_w(i,j  ,2)-ssh_int_w(i,j  ,1)-ssh_int_w(i,j  ,0)+ssh_int_w(i,j  ,-1)  &
         +invgrav*(q_t(i,j,k,2)      -q_t(i,j,k,1)      -q_t(i,j,k,0)      +q_t(i,j,k,-1)) & 
                        ) & !ooo>
                    )*inv_dti_fw_p2 )**2     &
                  )*inv_brkvar               &              
               )    & !MIN>
                *(0.5-asselin_nh)            &
                     +asselin_nh
!.............
#endif

#ifdef bidon
!.............
! AUGMENTER LE COEF D'ASSELIN SI INSTABILITE SUR Q
! x0_r4 est le coef d'asselin instantannE
      x0_r4=        &
! l'indice de "variation" varie entre 0 et 1 quand d2q/dt2=brkvar
! Quand ce dernier=1, alors le coef d'asselin = 1*(0.5-asselin_nh)+asselin_nh = 0.5 c.a.d. la valeur max du coef d'asselin
            min(1., & !MIN>
                  ( ( ( &
                         ( & !ooo>
         +invgrav*( & !m°v°m>
                   max( abs(q_t(i,j,k,2)-2*q_t(i,j,k,1)+q_t(i,j,k,0 ))  & ! max1
                       ,abs(q_t(i,j,k,1)-2*q_t(i,j,k,0)+q_t(i,j,k,-1))) & ! max2
                  ) & !m°v°m>
                        ) & !ooo>
                    )*inv_dti_fw_p2 )**2     &
                  )*inv_brkvar               &              
               )    & !MIN>
                *(0.5-asselin_nh)            &
                     +asselin_nh
!.............
#endif

!      if(i==imax/2.and.j==jmax/2.and.k==kmax/2.and.par%rank==0) then
!       write(10+par%rank,'(f12.3,1x,2(1x,e14.7),a)') &
!            elapsedtime_now  &
!           ,x0_r4            &
!           ,breaker2d_t(i,j) &
!           ,' trouvemoi87' 
!      endif

! FILTRE D'ASSELIN coef constant
            q_t(i,j,k,1)=x1_r4*q_t(i,j,k,1)+x2_r4*(q_t(i,j,k,0)+q_t(i,j,k,2))
! FILTRE D'ASSELIN coef ajustE
!           q_t(i,j,k,1)=(1.-x0_r4)*q_t(i,j,k,1)+0.5*x0_r4*(q_t(i,j,k,0)+q_t(i,j,k,2))
       enddo
       enddo
       enddo

! MOVE FORWARD

        do k=0,kmax+1
        do j=0,jmax+1
        do i=0,imax+1
            q_t(i,j,k,-1)=q_t(i,j,k,0)
            q_t(i,j,k,0 )=q_t(i,j,k,1)
            q_t(i,j,k,1 )=q_t(i,j,k,2)
       enddo
       enddo
       enddo

      end subroutine q_moveforward

!............................

      subroutine q_add_nhpgf2d
      implicit none

!     if(nh_scheme==1)call q_nhpgf2d_leplussimple
!     if(nh_scheme==2)call q_nhpgf2d
!     if(nh_scheme==3)call q_nhpgf2d_leplussimple
! Finalement tous les schemas fonctionnent mieux avec q_nhpgf2d:
      call q_nhpgf2d
!     call q_u3d_from_u2d

       do j=2,jmax-1
       do i=2,imax
                pgf_u(i,j,1)=  &
                pgf_u(i,j,1)   &
         +0.5*nhpgf2d_u(i,j) ! 0.5 car on ne prend que la moitie A cause de la partie pgf_u(0)
       enddo
       enddo

       do j=2,jmax
       do i=2,imax-1
                pgf_v(i,j,1)=  &
                pgf_v(i,j,1)   &
         +0.5*nhpgf2d_v(i,j) ! 0.5 car on ne prend que la moitie A cause de la partie pgf_u(0)
       enddo
       enddo

! Ajout d'un terme d'amplification des vagues:
!      call q_wind_induced_waves

      end subroutine q_add_nhpgf2d

!............................

      subroutine q_nhpgf2d
      implicit none

      id_dz=1 ; id_zt=2

! Calcul de la grille sigma si:
! - nh_frozensigma=0
! - Etat initial (iteration2d==1)
!     if(nh_frozensigma==0.or.iteration2d==1) then !m[°v°]m>
      if(computesigma==1) then                     !m[°v°]m>

      if(nh_frozensigma==0) then ! grid updated >
         do j=1,jmax ; do i=1,imax
!          sshr4_w(i,j)=0.999*sshr4_w(i,j)+0.001*max(ssh_w(i,j,1),0.1-h_w(i,j))
           sshr4_w(i,j)=max(ssh_w(i,j,1),0.001-h_w(i,j))
         enddo       ; enddo
      else                       ! frozen grid >
         do j=1,jmax ; do i=1,imax
            sshr4_w(i,j)=max(0.,0.001-h_w(i,j))
         enddo       ; enddo
      endif                      ! frozen grid >

! Definir dz_t, l'epaisseur au point q:
      x0=1./real(kmax) ; id_dz=1
      do k=0,kmax   ; do j=1,jmax   ; do i=1,imax
       anyv3dr4(i,j,k,id_dz)=x0*(h_w(i,j)+sshr4_w(i,j))
      enddo         ; enddo         ; enddo
! Definit z_t la profondeur de q
! Rappel: z_t(k)=z_t(k+1)-0.5*( dz_t(k+1)+dz_t(k) )
! arbitrairement z_t(kmax+1)=surface+0.5*dz_t(kmax+1)
      do j=1,jmax   ; do i=1,imax
       anyv3dr4(i,j,kmax+1,id_zt)=sshr4_w(i,j) ! +0.5*anyv3dr4(i,j,kmax+1,id_dz)
      enddo         ; enddo
      do j=1,jmax   ; do i=1,imax
       anyv3dr4(i,j,kmax,id_zt)=anyv3dr4(i,j,kmax+1,id_zt)-0.5*anyv3dr4(i,j,kmax,id_dz)
      enddo         ; enddo
      do k=kmax-1,0,-1 ; do j=1,jmax   ; do i=1,imax
       anyv3dr4(i,j,k,id_zt)=anyv3dr4(i,j,k+1,id_zt)         &
               -0.5*(  anyv3dr4(i,j,k+1,id_dz)   &
                      +anyv3dr4(i,j,k  ,id_dz))
      enddo          ; enddo         ; enddo

      checkr0=anyv3dr4(imax/2,jmax/2,kmax/2,id_dz)
      checkr1=anyv3dr4(imax/2,jmax/2,kmax/2,id_zt)

      if(nh_frozensigma==1)computesigma=0 ! A la prochaine iteration on ne passera plus par ces lignes

      endif                                       !m[°v°]m>

! Etape 1: calculer dz*nhgpf avec q au temps t+1 sachant qu'a l'iteration suivante
! cette info sera reutilisee pour calculer la membre de droite du laplacien non-hydrostatique

! Calculer dz*dq/dx-dz/dx*dq/dk:
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        nhpgf_u(i,j,k)=( &
        0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))*(q_t(i,j,k,1)-q_t(i-1,j,k,1))     &
       -(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i-1,j,k,id_zt))*0.25*( q_t(i,j,k+1,1) +q_t(i-1,j,k+1,1)           &
                                                  -q_t(i,j,k-1,1) -q_t(i-1,j,k-1,1))          &
                      )*invdx_u(i,j)*nhpgf_reduce
      enddo ; enddo ; enddo

! Calculer dz*dq/dy-dz/dy*dq/dk:
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        nhpgf_v(i,j,k)=( &
        0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))*(q_t(i,j,k,1)-q_t(i,j-1,k,1))     &
       -(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j-1,k,id_zt))*0.25*( q_t(i,j,k+1,1) +q_t(i,j-1,k+1,1)           &
                                                  -q_t(i,j,k-1,1) -q_t(i,j-1,k-1,1))          &
                      )*invdy_v(i,j)*nhpgf_reduce
      enddo ; enddo ; enddo

! Etape 2: calculer la moyenne verticale du nhpgf pour calcul du mode externe:

! Cette boucle est grande pour des imperatifs de continuite mpi
       k=1
       do j=2,jmax-1 ! 1,jmax   ! 2,jmax-1
       do i=2,imax   ! 1,imax+1 ! 2,imax
! note: si la mutiplication par dz est faite en amont alors commenter la fin de la ligne:
        nhpgf2d_u(i,j)=nhpgf_u(i,j,k) ! *(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))
!       nhpgf2d_u(i,j)=anyv3dr4(i,j,k,id_dzpgfu) ! *(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))
       enddo
       enddo
       do k=2,kmax
       do j=2,jmax-1 ! 1,jmax   ! 2,jmax-1
       do i=2,imax   ! 1,imax+1 ! 2,imax
        nhpgf2d_u(i,j)= &
        nhpgf2d_u(i,j)+nhpgf_u(i,j,k) ! anyv3dr4(i,j,k,id_dzpgfu) ! *(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))
! note: si la mutiplication par dz est faite en amont alors commenter la fin de la ligne
       enddo
       enddo
       enddo

! Cette boucle est grande pour des imperatifs de continuite mpi
       k=1
       do j=2,jmax   ! 1,jmax+1 ! 2,jmax
       do i=2,imax-1 ! 1,imax   ! 2,imax-1
! note: si la mutiplication par dz est faite en amont alors commenter la fin de la ligne:
        nhpgf2d_v(i,j)=nhpgf_v(i,j,k) ! anyv3dr4(i,j,k,id_dzpgfv) ! *(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))
       enddo
       enddo
       do k=2,kmax
       do j=2,jmax   ! 1,jmax+1 ! 2,jmax
       do i=2,imax-1 ! 1,imax   ! 2,imax-1
        nhpgf2d_v(i,j)=  &
! note: si la mutiplication par dz est faite en amont alors commenter la fin de la ligne:
        nhpgf2d_v(i,j)+nhpgf_v(i,j,k) ! anyv3dr4(i,j,k,id_dzpgfv) ! *(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))
       enddo
       enddo
       enddo

      end subroutine q_nhpgf2d

!............................

      subroutine q_nhpgf2d_leplussimple
      implicit none

! Calculer la pression qnh moyenne (sur la verticale).
! Que ktop=kmax ou ktop=kmax+1 cela ne change rien au fait qu'il y
! a kmax niveaux (simplement dans un cas q(kmax)=0 et dans l'autre non
! car c'est q(kmax+1) qui l'est (mais ne compte pas dans l'intergrale)
! L'integrale commence donc toujours A k=kmax
      k=kmax
      do j=1,jmax ; do i=1,imax
       q2davr_w(i,j)=q_t(i,j,k,1)*dsigr4_t(i,j,k)
      enddo       ; enddo
      do k=kmax-1,1,-1
      do j=1,jmax ; do i=1,imax
       q2davr_w(i,j)=  &
       q2davr_w(i,j)+q_t(i,j,k,1)*dsigr4_t(i,j,k)
      enddo       ; enddo
      enddo

! Cette boucle est grande pour des imperatifs de continuite mpi
       do j=2,jmax-1
       do i=2,imax

        nhpgf2d_u(i,j)= &
         mask_u(i,j,kmax)*( & !m0v0m>

! Methode simplifiee ignorant les variations de H et h:
!           (h_u(i,j)+0.5*(ssh_w(i,j,1 )+ssh_w(i-1,j,1)))*(           & ! (h+ssh)*[
!                       q2davr_w(i,j)-q2davr_w(i-1,j)                 & ! +dq/di
!                                                         )/dx_u(i,j) & ! ]/dx 

        (    (h_w(i  ,j)+ssh_w(i  ,j,1 ))*q2davr_w(i  ,j)        &
            -(h_w(i-1,j)+ssh_w(i-1,j,1 ))*q2davr_w(i-1,j)        &

       +0.5*(q_t(i,j,1,1)+q_t(i-1,j,1,1))*(h_w(i,j)-h_w(i-1,j))  &

        )/dx_u(i,j)                                              &
                          )   !m0v0m>

       enddo
       enddo

! Cette boucle est grande pour des imperatifs de continuite mpi
       do j=2,jmax
       do i=2,imax-1

        nhpgf2d_v(i,j)= &
         mask_v(i,j,kmax)*( & !m0v0m>

! Methode simplifiee ignorant les variations de H et h:
!           (h_v(i,j)+0.5*(ssh_w(i,j,1) +ssh_w(i,j-1,1)))*(           & ! (h+ssh)*[
!                       q2davr_w(i,j)-q2davr_w(i,j-1)                 & ! +dq/dj
!                                                         )/dy_v(i,j) & ! ]/dy

        (    (h_w(i,j  )+ssh_w(i,j  ,1 ))*q2davr_w(i,j  )        &
            -(h_w(i,j-1)+ssh_w(i,j-1,1 ))*q2davr_w(i,j-1)        &

       +0.5*(q_t(i,j,1,1)+q_t(i,j-1,1,1))*(h_w(i,j)-h_w(i,j-1))  &

        )/dy_v(i,j)                                              & 
                          )   !m0v0m>

       enddo
       enddo

      end subroutine q_nhpgf2d_leplussimple

!............................

      subroutine q_moteur_leplussimple
      implicit none

! rappel: z_t(k+1)-z_t(k)=H*dsig_w(k+1)
! rappel: z_w(k+1)-z_w(k)=H*dsig_t(k)

      ww_=w2_ ; if(ww_==0.)ww_=1.
      do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
       mc(i,j,k,1)=-ww_*dt2/( dsigr4_t(i,j,k)*dsig_w(i,j,k)*(h_w(i,j)+ssh_w(i,j,1))**2 )
      enddo       ; enddo       ;enddo

      do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
       mc(i,j,k,3)=mc(i,j,k,1)*dsig_w(i,j,k)/dsig_w(i,j,k+1)
      enddo       ; enddo       ;enddo

      do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
       mc(i,j,k,2)=1.-mc(i,j,k,1)-mc(i,j,k,3)
      enddo       ; enddo       ;enddo

!.......................................................
! C.L. Surface: q_t(ktop)=0.
      k=ktop
      do j=2,jmax-1 ; do i=2,imax-1
       mc(i,j,k,1)=0.
       mc(i,j,k,2)=1.
       mc(i,j,k,3)=0.
       mc(i,j,k,4)=0.
      enddo       ; enddo

!.......................................................
! C.L. Fond: q_t(0)=q_t(1)
      k=1
      do j=2,jmax-1 ; do i=2,imax-1
       mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
       mc(i,j,k,1)=0.
      enddo       ; enddo

      ww_=1.
      if(w2_/=0.)ww_=1./w2_

      do j=2,jmax-1 ; do i=2,imax-1

!      anyvar2d(i,j)=dt2*( & !ttt>

!      invdx_t(i,j)*( & !ooo>
!      ((velbar_u(i+1,j,1)-velbar_u(i+1,j,2))/dte_fw-nhpgf2d_u(i+1,j)/(h_u(i+1,j)+0.5*(ssh_w(i+1,j,1)+ssh_w(i,j,1)))) &
!     -((velbar_u(i  ,j,1)-velbar_u(i  ,j,2))/dte_fw-nhpgf2d_u(i  ,j)/(h_u(i  ,j)+0.5*(ssh_w(i-1,j,1)+ssh_w(i,j,1)))) &
!                   ) & !ooo>                                           

!     +invdy_t(i,j)*( & !ooo>
!      ((velbar_v(i,j+1,1)-velbar_v(i,j+1,2))/dte_fw-nhpgf2d_v(i,j+1)/(h_v(i,j+1)+0.5*(ssh_w(i,j+1,1)+ssh_w(i,j,1)))) &
!     -((velbar_v(i,j  ,1)-velbar_v(i,j  ,2))/dte_fw-nhpgf2d_v(i,j  )/(h_v(i,j  )+0.5*(ssh_w(i,j-1,1)+ssh_w(i,j,1)))) &
!                   ) & !ooo>

!                        )    !ttt>

!      anyvar2d(i,j)=dt2*( & !ttt>

!      invdx_t(i,j)*( & !ooo>
!            (grav*(ssh_w(i+1,j,1)-ssh_w(i  ,j,1))*invdx_u(i+1,j)) &
!           -(grav*(ssh_w(i  ,j,1)-ssh_w(i-1,j,1))*invdx_u(i  ,j)) &
!                   ) & !ooo>                                           

!     +invdy_t(i,j)*( & !ooo>
!            (grav*(ssh_w(i,j+1,1)-ssh_w(i,j  ,1))*invdy_v(i,j+1)) &
!           -(grav*(ssh_w(i,j  ,1)-ssh_w(i,j-1,1))*invdy_v(i,j  )) &
!                   ) & !ooo>

!                        )    !ttt>

       anyvar2d(i,j)=dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
       ((fluxbar_u(i+1,j,0)-fluxbar_u(i+1,j,1))/(dte_fw*dy_u(i+1,j))-nhpgf2d_u(i+1,j))/(h_u(i+1,j)+0.5*(ssh_w(i+1,j,1)+ssh_w(i,j,1))) &
      -((fluxbar_u(i  ,j,0)-fluxbar_u(i  ,j,1))/(dte_fw*dy_u(i  ,j))-nhpgf2d_u(i  ,j))/(h_u(i  ,j)+0.5*(ssh_w(i-1,j,1)+ssh_w(i,j,1))) &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
       ((fluxbar_v(i,j+1,0)-fluxbar_v(i,j+1,1))/(dte_fw*dx_v(i,j+1))-nhpgf2d_v(i,j+1))/(h_v(i,j+1)+0.5*(ssh_w(i,j+1,1)+ssh_w(i,j,1))) &
      -((fluxbar_v(i,j  ,0)-fluxbar_v(i,j  ,1))/(dte_fw*dx_v(i,j  ))-nhpgf2d_v(i,j  ))/(h_v(i,j  )+0.5*(ssh_w(i,j-1,1)+ssh_w(i,j,1))) &
                    ) & !ooo>

                         )    !ttt>

      enddo ; enddo

      do k=1,ktop-1 
      do j=2,jmax-1 ; do i=2,imax-1

        mc(i,j,k,4)=2*q_t(i,j,k,1)-q_t(i,j,k,0)                        &

                +anyvar2d(i,j)   &
                                              +dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
                ( q_t(i+1,j,k,1)-q_t(i  ,j,k,1) )*invdx_u(i+1,j)    &
               -( q_t(i  ,j,k,1)-q_t(i-1,j,k,1) )*invdx_u(i  ,j)    &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
                ( q_t(i,j+1,k,1)-q_t(i  ,j,k,1) )*invdy_v(i,j+1)    &
               -( q_t(i,j  ,k,1)-q_t(i,j-1,k,1) )*invdy_v(i  ,j)    &
                    ) & !ooo>

                                                   ) &  !ttt>


       -mc(i,j,k,1)    *ww_*(w0_*q_t(i,j,k-1,0)+w1_*q_t(i,j,k-1,1) )&
       -mc(i,j,k,3)    *ww_*(w0_*q_t(i,j,k+1,0)+w1_*q_t(i,j,k+1,1) )&
      -(mc(i,j,k,2)-1.)*ww_*(w0_*q_t(i,j,k  ,0)+w1_*q_t(i,j,k  ,1) ) 

! Note sur sponge_t: varie de 0 a 1 (sur la frontiere) ce qui equivaut en termes
! d'echelle de temps associee a un temps de rappel equivalent au pas de temps. C'est
! donc un rappel fort qui a pour but de supprimer la pression NH aux frontieres ouvertes

      enddo ; enddo
      enddo

! Solveur vertical tridiagonal:
      call q_solver(1,2,imax-1,2,jmax-1,1,ktop) ! donne q_t(t+1)

! Reporter les hypotheses de conditions aux limites au fond et en surface:
      do j=2,jmax-1 ; do i=2,imax-1
! C.L.: q(0)=q(1):
       q_t(i,j,0,2)=q_t(i,j,1,2)*mask_t(i,j,kmax)
! C.L.: q(kmax+1)=0.
       q_t(i,j,ktop:kmax+1,2)=0.
      enddo ; enddo

! ZONE EPONGE: (Ne pourrait elle pas etre implicitee?)
      do k=1,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1
       q_t(i,j,k,2)=q_t(i,j,k,2)*(1.-sponge_t(i,j,1))
      enddo         ; enddo       ; enddo

      call q_obc_mpi('z0') ! mpi 'z1' sur q_t(t+1)

      end subroutine q_moteur_leplussimple

!............................

      subroutine q_moteur_sigma
      implicit none

      if(checkr0/=anyv3dr4(imax/2,jmax/2,kmax/2,id_dz)) &
         stop 'err nh 1282'

      if(checkr1/=anyv3dr4(imax/2,jmax/2,kmax/2,id_zt))       &
         stop 'err nh 1286'

! Calculer (dz/dx)*[dq/dx-dz/dx*dq/dz]:
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       anyv3dr4(i,j,k,id_dzodxpgfu)=                          &
        nhpgf_u(i,j,k)                                        &
       /(0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))) &
       *(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i-1,j,k,id_zt))*invdx_u(i,j)

      enddo ; enddo ; enddo

! Calculer (dz/dy)*[dq/dy-dz/dy*dq/dz]:
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
       anyv3dr4(i,j,k,id_dzodypgfv)=                          &
        nhpgf_v(i,j,k)                                        &
       /(0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))) &
       *(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j-1,k,id_zt))*invdy_v(i,j)
      enddo ; enddo ; enddo

      ww_=w2_ ; if(ww_==0.)ww_=1.
      do k=1,ktop-1
      do j=2,jmax-1 ! 1,jmax
      do i=2,imax-1 ! 1,imax
       mc(i,j,k,1)=-ww_*dt2/((anyv3dr4(i,j,k  ,id_zt)-anyv3dr4(i,j,k-1,id_zt))*anyv3dr4(i,j,k,id_dz))
       mc(i,j,k,3)=-ww_*dt2/((anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k  ,id_zt))*anyv3dr4(i,j,k,id_dz))
      enddo
      enddo
      enddo

      do k=1,ktop-1
      do j=2,jmax-1 ! 1,jmax
      do i=2,imax-1 ! 1,imax
       mc(i,j,k,2)=1.-mc(i,j,k,1)-mc(i,j,k,3)
      enddo
      enddo
      enddo

!.......................................................
! C.L. Surface: q_t(ktop)=0.
      k=ktop
      do j=2,jmax-1 ; do i=2,imax-1 !     do j=1,jmax ; do i=1,imax
       mc(i,j,k,1)=0.
       mc(i,j,k,2)=1.
       mc(i,j,k,3)=0.
       mc(i,j,k,4)=0.
      enddo       ; enddo

!.......................................................
! C.L. Fond: q_t(0)=q_t(1)
      k=1
      do j=2,jmax-1 ; do i=2,imax-1 !     do j=1,jmax ; do i=1,imax
       mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
       mc(i,j,k,1)=0.
      enddo       ; enddo

#ifdef stokes
      stop 'ici fluxbar ne peut pas etre le courant lagrangien'
#endif

!.......................................................
! Tendance hydrostatique:
      do j=2,jmax-1 ; do i=2,imax-1 !     do j=1,jmax ; do i=1,imax

!      anyvar2d(i,j)=dt2*( & !ttt>

!      invdx_t(i,j)*( & !ooo>
!      ((velbar_u(i+1,j,1)-velbar_u(i+1,j,2))/dte_fw-nhpgf2d_u(i+1,j)/(h_u(i+1,j)+0.5*(ssh_w(i+1,j,1)+ssh_w(i,j,1)))) &
!     -((velbar_u(i  ,j,1)-velbar_u(i  ,j,2))/dte_fw-nhpgf2d_u(i  ,j)/(h_u(i  ,j)+0.5*(ssh_w(i-1,j,1)+ssh_w(i,j,1)))) &
!                   ) & !ooo>                                           

!     +invdy_t(i,j)*( & !ooo>
!      ((velbar_v(i,j+1,1)-velbar_v(i,j+1,2))/dte_fw-nhpgf2d_v(i,j+1)/(h_v(i,j+1)+0.5*(ssh_w(i,j+1,1)+ssh_w(i,j,1)))) &
!     -((velbar_v(i,j  ,1)-velbar_v(i,j  ,2))/dte_fw-nhpgf2d_v(i,j  )/(h_v(i,j  )+0.5*(ssh_w(i,j-1,1)+ssh_w(i,j,1)))) &
!                   ) & !ooo>

!                        )    !ttt>

        anyvar2d(i,j)=dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
       ((fluxbar_u(i+1,j,0)-fluxbar_u(i+1,j,1))/(dte_fw*dy_u(i+1,j))-nhpgf2d_u(i+1,j))/(h_u(i+1,j)+0.5*(ssh_w(i+1,j,1)+ssh_w(i,j,1))) &
      -((fluxbar_u(i  ,j,0)-fluxbar_u(i  ,j,1))/(dte_fw*dy_u(i  ,j))-nhpgf2d_u(i  ,j))/(h_u(i  ,j)+0.5*(ssh_w(i-1,j,1)+ssh_w(i,j,1))) &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
       ((fluxbar_v(i,j+1,0)-fluxbar_v(i,j+1,1))/(dte_fw*dx_v(i,j+1))-nhpgf2d_v(i,j+1))/(h_v(i,j+1)+0.5*(ssh_w(i,j+1,1)+ssh_w(i,j,1))) &
      -((fluxbar_v(i,j  ,0)-fluxbar_v(i,j  ,1))/(dte_fw*dx_v(i,j  ))-nhpgf2d_v(i,j  ))/(h_v(i,j  )+0.5*(ssh_w(i,j-1,1)+ssh_w(i,j,1))) &
                    ) & !ooo>

                         )    !ttt>
      enddo ; enddo

!.......................................................
! Laplacien non-hydrostatique:
      ww_=1.
      if(w2_/=0.)ww_=1./w2_
      do k=1,ktop-1 
      do j=2,jmax-1 ; do i=2,imax-1 !     do j=1,jmax ; do i=1,imax

        mc(i,j,k,4)=2*q_t(i,j,k,1)-q_t(i,j,k,0)                        &

                +anyvar2d(i,j)   &

      +dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
                 nhpgf_u(i+1,j,k)      &
                -nhpgf_u(i  ,j,k)      &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
                 nhpgf_v(i,j+1,k)      &
                -nhpgf_v(i,j  ,k)      &
                    ) & !ooo>

      -0.25*( & !(0L0)>
              anyv3dr4(i+1,j  ,k+1,id_dzodxpgfu)  &
             +anyv3dr4(i  ,j  ,k+1,id_dzodxpgfu)  &
             -anyv3dr4(i+1,j  ,k-1,id_dzodxpgfu)  &
             -anyv3dr4(i  ,j  ,k-1,id_dzodxpgfu)  &
             +anyv3dr4(i  ,j+1,k+1,id_dzodypgfv)  &
             +anyv3dr4(i  ,j  ,k+1,id_dzodypgfv)  &
             -anyv3dr4(i  ,j+1,k-1,id_dzodypgfv)  &
             -anyv3dr4(i  ,j  ,k-1,id_dzodypgfv)  &
            ) & !(0L0)>


           ) &  !ttt>
            /anyv3dr4(i,j,k,id_dz)                &
                                                    

       -mc(i,j,k,1)    *ww_*(w0_*q_t(i,j,k-1,0)+w1_*q_t(i,j,k-1,1) )&
       -mc(i,j,k,3)    *ww_*(w0_*q_t(i,j,k+1,0)+w1_*q_t(i,j,k+1,1) )&
      -(mc(i,j,k,2)-1.)*ww_*(w0_*q_t(i,j,k  ,0)+w1_*q_t(i,j,k  ,1) ) 

! Note sur sponge_t: varie de 0 a 1 (sur la frontiere) ce qui equivaut en termes
! d'echelle de temps associee a un temps de rappel equivalent au pas de temps. C'est
! donc un rappel fort qui a pour but de supprimer la pression NH aux frontieres ouvertes

#ifdef bidon
       if(i==40.and.j==40) then
        write(6,*)'toto' &
      ,invdx_t(i,j)*( & !ooo>
                ( q_t(i+1,j,k,1)-q_t(i  ,j,k,1) )*invdx_u(i+1,j)    &
               -( q_t(i  ,j,k,1)-q_t(i-1,j,k,1) )*invdx_u(i  ,j)    &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
                ( q_t(i,j+1,k,1)-q_t(i  ,j,k,1) )*invdy_v(i,j+1)    &
               -( q_t(i,j  ,k,1)-q_t(i,j-1,k,1) )*invdy_v(i  ,j)    &
                    ) & !ooo>
      ,    ( & !ttt>

       invdx_t(i,j)*( & !ooo>
                 nhpgf_u(i+1,j,k)      &
                -nhpgf_u(i  ,j,k)      &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
                 nhpgf_v(i,j+1,k)      &
                -nhpgf_v(i,j  ,k)      &
                    ) & !ooo>

      -0.25*( & !(0L0)>
              anyv3dr4(i+1,j  ,k+1,id_dzodxpgfu)  &
             +anyv3dr4(i  ,j  ,k+1,id_dzodxpgfu)  &
             -anyv3dr4(i+1,j  ,k-1,id_dzodxpgfu)  &
             -anyv3dr4(i  ,j  ,k-1,id_dzodxpgfu)  &
             +anyv3dr4(i  ,j+1,k+1,id_dzodypgfv)  &
             +anyv3dr4(i  ,j  ,k+1,id_dzodypgfv)  &
             -anyv3dr4(i  ,j+1,k-1,id_dzodypgfv)  &
             -anyv3dr4(i  ,j  ,k-1,id_dzodypgfv)  &
            ) & !(0L0)>


           ) &  !ttt>
            /anyv3dr4(i,j,k,id_dz)                 
       endif
#endif

      enddo ; enddo
      enddo

! Solveur vertical tridiagonal:
      call q_solver(1,2,imax-1,2,jmax-1,1,ktop) ! donne q_t(t+1)

! Reporter les hypotheses de conditions aux limites au fond et en surface:
      do j=2,jmax-1 ; do i=2,imax-1 !     do j=1,jmax ; do i=1,imax
! C.L.: q(0)=q(1):
       q_t(i,j,0,2)=q_t(i,j,1,2)*mask_t(i,j,kmax)
! C.L.: q(kmax+1)=0.
       q_t(i,j,ktop:kmax+1,2)=0.
      enddo ; enddo

! ZONE EPONGE: (Ne pourrait elle pas etre implicitee?)
      do k=1,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1 !     do j=1,jmax ; do i=1,imax
       q_t(i,j,k,2)=q_t(i,j,k,2)*(1.-sponge_t(i,j,1))
      enddo         ; enddo       ; enddo

      call q_obc_mpi('z0') ! mpi 'z1' sur q_t(t+1)

      end subroutine q_moteur_sigma

!-----------------------------------------------------------------------------

      subroutine q_diag_stokes_transport(i_)
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
      integer i_

!     if(iteration2d<iteration2d_max_now-1000)return

      i=i_ ; j=jmax/2

      if(.not.allocated(stokes_transport_sum1)) then !***>
        allocate(stokes_transport_sum1(1)) ; stokes_transport_sum1=0.
        allocate(stokes_transport_sum2(1)) ; stokes_transport_sum2=0.
        allocate(stokes_transport_sum3(1)) ; stokes_transport_sum3=0.
      endif                                          !***>


       stokes_transport_sum1(1)=stokes_transport_sum1(1)+dte_fw
       stokes_transport_sum2(1)=stokes_transport_sum2(1)+dte_fw*fluxbar_u(i,j,1)
       stokes_transport_sum3(1)=stokes_transport_sum3(1)+dte_fw*fluxbar_v(i,j,1)

!      if(  int((iteration2d  )*dte_fw/period_) &
!          -int((iteration2d-1)*dte_fw/period_)==1) then
!      if(mod(iteration2d,200)==0) then
       if(velbar_u(i,j,1)>0.and.velbar_u(i,j,0)<0)stokes_count=stokes_count+1
       if(stokes_count==100) then
        write(texte30,'(a,i0)')'tmp/nhdiag_stokestransport',i
        open(unit=3,file=texte30,position='append')
        write(3,*)iteration2d*dte_fw/3600. &
                 ,1.e5*stokes_transport_sum2(1)/dy_u(i,j)/stokes_transport_sum1(1) &
                 ,1.e5*stokes_transport_sum3(1)/dx_v(i,j)/stokes_transport_sum1(1)  
        close(3)
        stokes_transport_sum1=0.
        stokes_transport_sum2=0.
        stokes_transport_sum3=0.
        stokes_count=0
       endif

!      write(66,*)ssh_w(i,j,1),velbar_u(i,j,1),q_t(i,j,kmax,1)

      end subroutine q_diag_stokes_transport

!-----------------------------------------------------------------------------

#ifdef checkmpi
      subroutine q_check_mpi_conservation
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
      integer loop_

! Pour verifier la conservation de la parallelisation:
      do k=1,kmax
        j1=1 ; j2=jmax
        do i=1,imax
         anyv3dr4(i,j1,k,1)=q_t(i,j1,k,1)
         anyv3dr4(i,j2,k,1)=q_t(i,j2,k,1)
        enddo
        i1=1 ; i2=imax
        do j=1,jmax
         anyv3dr4(i1,j,k,1)=q_t(i1,j,k,1)
         anyv3dr4(i2,j,k,1)=q_t(i2,j,k,1)
        enddo
      enddo

      call get_type_echange('z0','q_z0_n',q_t,lbound(q_t),ubound(q_t),1,i0)
      do loop_=1, subcycle_exchange
         call echange_voisin(q_t,i0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait()

      ksecu=0
      do k=1,kmax

      j1=1 ; j2=jmax
      do i=1,imax

      if(mask_t(i,j1,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3dr4(i,j1,k,1)/=q_t(i,j1,k,1))then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i,j1,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i+par%timax(1),j1+par%tjmax(1),k
        write(10+par%rank,*)'q_t ',anyv3dr4(i,j1,k,1),q_t(i,j1,k,1)
        write(10+par%rank,*)'delta ',anyv3dr4(i,j1,k,1)-q_t(i,j1,k,1)
        write(10+par%rank,*)'mask_t',j1,j1+1,  mask_t(i,j1:j1+1,k)
       endif
      endif                      !mmmmmmmmmmmmm>

      if(mask_t(i,j2,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3dr4(i,j2,k,1)/=q_t(i,j2,k,1))then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i,j2,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(nord)
        write(10+par%rank,*)'coordonnées ',i,2
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i+par%timax(1),j2+par%tjmax(1),k
        write(10+par%rank,*)'q_t ',anyv3dr4(i,j2,k,1),q_t(i,j2,k,1)
        write(10+par%rank,*)'delta ',anyv3dr4(i,j2,k,1)-q_t(i,j2,k,1)
        write(10+par%rank,*)'mask_t',j2-1,j2,  mask_t(i,j2-1:j2,k)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! i

      i1=1 ; i2=imax
      do j=1,jmax

      if(mask_t(i1,j,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3dr4(i1,j,k,1)/=q_t(i1,j,k,1)) then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i1,j,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(ouest)
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i1+par%timax(1),j+par%tjmax(1),k
        write(10+par%rank,*)'q_t ',anyv3dr4(i1,j,k,1),q_t(i1,j,k,1)
        write(10+par%rank,*)'delta ',anyv3dr4(i1,j,k,1),q_t(i1,j,k,1)
        write(10+par%rank,*)'mask_t',i1,i1+1,  mask_t(i1:i1+1,j,k)
       endif
      endif                      !mmmmmmmmmmmmm>

      if(mask_t(i2,j,k)==1) then !mmmmmmmmmmmmm>
       if(anyv3dr4(i2,j,k,1)/=q_t(i2,j,k,1)) then
        write(10+par%rank,*)'----------------------------' ; ksecu=1
        write(10+par%rank,*)'local i j k',i2,j,k
        write(10+par%rank,*)'du proc     ',par%rank
        write(10+par%rank,*)'proc voisin ',par%tvoisin(est)
        write(10+par%rank,*)'coordonnées',2,j
        write(10+par%rank,*)'imax jmax ',imax,jmax
        write(10+par%rank,*)'global i j k',i2+par%timax(1),j+par%tjmax(1),k
        write(10+par%rank,*)'q_t ',anyv3dr4(i2,j,k,1),q_t(i2,j,k,1)
        write(10+par%rank,*)'delta ',anyv3dr4(i2,j,k,1)-q_t(i2,j,k,1)
        write(10+par%rank,*)'mask_t',i2-1,i2,  mask_t(i2-1:i2,j,k)
       endif
      endif                      !mmmmmmmmmmmmm>

      enddo ! j

      enddo ! kmax

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(ksecu==1) then
        write(6,'(a,i0,a)')   &
        ' mpi conservation error for q_t Details in fort.' &
         ,10+par%rank,' available in RDIR/CONFIG directory'
        stop ' STOP in routine q_check_mpi_conservation!'
      endif

      end subroutine q_check_mpi_conservation
#endif

!-----------------------------------------------------------------------------
#ifdef bidon
      subroutine q_pseudo3d_sigma
      implicit none

      ub4=ubound(omega_w)
      if(ub4(4)==1) then !>>>
       deallocate(omega_w)
       allocate(omega_w(imax,jmax,kmax+1,0:2)) ; omega_w=0.
      endif              !>>>

! COMPOSANTE U
!      sum1=0.
!      sum2=0.
!      sum6=0.
!      sum3=0.
!      sum5=0.
       do k=1,kmax
       do j=2,jmax-1
       do i=2,imax

         x1=dsigr4_t(i,j,k)*(h_u(i,j)+0.5*(ssh_w(i,j,0)+ssh_w(i-1,j,0)))
         x2=dsigr4_t(i,j,k)*(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1)))

! x3 pgf de ssh
         x3=x2*grav*(ssh_w(i  ,j,1)-ssh_w(i-1,j,1))/dx_u(i,j)

         vel_u(i,j,k,2)=mask_u(i,j,k)*(vel_u(i,j,k,1)*x1 &

           +dte_fw*( & !ooo>

              -x3-nhpgf_u(i,j,k) &

                   ) & !ooo>

                         )/x2

!       sum1=sum1+0.5*(vel_u(i,j,k,2)+vel_u(i,j,k,1))*x3 ! bilan ssh
!       sum6=sum6+0.5*(vel_u(i,j,k,2)+vel_u(i,j,k,1))*nhpgf_u(i,j,k) ! anyv3dr4(i,j,k,id_dzpgfu) ! bilan q


!       sum2=sum2+0.5*x2*(vel_u(i,j,k,2)+vel_u(i,j,k,1))*(vel_u(i,j,k,2)-vel_u(i,j,k,1))/dte_fw
!       sum3=sum3+0.5*x2*(vel_u(i,j,k,1))*(vel_u(i,j,k,1))
!       sum5=sum5+0.5*x2*(vel_u(i,j,k,2))*(vel_u(i,j,k,2))

       enddo
       enddo
       enddo

! Verifier adequation avec equation pour le courant moyen:
!     i=imax/2 ; j=jmax/2
!     sum0=0.
!     do k=1,kmax
!      sum0=sum0+dsigr4_t(i,j,k)*(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1)))*vel_u(i,j,k,2)
!     enddo
!     write(6,*)sum0/(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))),velbar_u(i,j,2)

!      stop 'caravane'
!#ifdef bidon


!      sum7=0.
!      do k=1,kmax
!      do j=1,jmax
!      do i=1,imax
!       sum7=sum7-mask_t(i,j,k)*q_t(i,j,k,1)/dx_u(i,j) &
!                *0.5*( x0*h_u(i+1,j)*vel_u(i+1,j,k,2)-x0*h_u(i,j)*vel_u(i,j,k,2) &
!                      +x0*h_u(i+1,j)*vel_u(i+1,j,k,1)-x0*h_u(i,j)*vel_u(i,j,k,1) )
!      enddo
!      enddo
!      enddo

! COMPOSANTE V
!      sum2=0.
       do k=1,kmax
       do j=2,jmax
       do i=2,imax-1


         x1=dsigr4_t(i,j,k)*(h_v(i,j)+0.5*(ssh_w(i,j,0)+ssh_w(i,j-1,0)))
         x2=dsigr4_t(i,j,k)*(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1)))

! x3 pgf de ssh
         x3=x2*grav*(ssh_w(i  ,j,1)-ssh_w(i,j-1,1))/dy_v(i,j)

         vel_v(i,j,k,2)=mask_v(i,j,k)*(vel_v(i,j,k,1)*x1 &

           +dte_fw*( & !ooo>

              -x3-nhpgf_v(i,j,k) &

                   ) & !ooo>

                         )/x2

!      sum1=sum1+0.5*(vel_v(i,j,k,2)+vel_v(i,j,k,1)) &
!                   *x2*grav*(ssh_w(i  ,j,1)-ssh_w(i,j-1,1))/dy_v(i,j)
!      sum6=sum6+0.5*(vel_v(i,j,k,2)+vel_v(i,j,k,1)) &
!                        *x2*(q_t(i  ,j,k,1)-q_t(i,j-1,k,1))/dy_v(i,j)

!       sum2=sum2+0.5*x2*(vel_v(i,j,k,2)+vel_v(i,j,k,1))*(vel_v(i,j,k,2)-vel_v(i,j,k,1))/dte_fw

!       sum3=sum3+0.5*x2*(vel_v(i,j,k,1))*(vel_v(i,j,k,1))
!       sum5=sum5+0.5*x2*(vel_v(i,j,k,2))*(vel_v(i,j,k,2))

       enddo
       enddo
       enddo
!      sum_t_udssh=sum_t_udssh+dte_fw*sum1

!      do k=1,kmax
!      do j=1,jmax
!      do i=1,imax
!       sum7=sum7-mask_t(i,j,k)*q_t(i,j,k,1)/dy_v(i,j) &
!                *0.5*( x0*h_v(i,j+1)*vel_v(i,j+1,k,2)-x0*h_v(i,j)*vel_v(i,j,k,2) &
!                      +x0*h_v(i,j+1)*vel_v(i,j+1,k,1)-x0*h_v(i,j)*vel_v(i,j,k,1) )
!      enddo
!      enddo
!      enddo
!      write(6,*)'riri',sum6,sum7

! COMPOSANTE W
!      sum4=0.
       do j=1,jmax
       do i=1,imax
       do k=2,kmax+1

         omega_w(i,j,k,2)=mask_t(i,j,k)*(omega_w(i,j,k,1)    &

           +dte_fw*( & !ooo>

            ( & 
               -w0_*(    q_t(i,j,k,0)   -q_t(i,j,k-1,0) ) &
               -w1_*(    q_t(i,j,k,1)   -q_t(i,j,k-1,1) ) &
               -w2_*( q_t(i,j,k,2)  -q_t(i,j,k-1,2)   ) &

            ) /( anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)   ) &

                   ) & !ooo>
                                        )
        
!       sum4=sum4+0.5*(omega_w(i,j,k,2)+omega_w(i,j,k,1)) &
!             *(w0_*(    q_t(i,j,k,0)   -q_t(i,j,k-1,0) ) &
!              +w1_*(    q_t(i,j,k,1)   -q_t(i,j,k-1,1) ) &
!              +w2_*( q_t(i,j,k,2)  -q_t(i,j,k-1,2)   ) )

!       sum2=sum2+0.5*(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)) &
!                    *(omega_w(i,j,k,2)+omega_w(i,j,k,1)) &
!                    *(omega_w(i,j,k,2)-omega_w(i,j,k,1))/dte_fw

!       sum3=sum3+0.5*(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)) &
!                    *(omega_w(i,j,k,1)) &
!                    *(omega_w(i,j,k,1))

!       sum5=sum5+0.5*(anyv3dr4(i,j,k,id_zt)-anyv3dr4(i,j,k-1,id_zt)) &
!                    *(omega_w(i,j,k,2)) &
!                    *(omega_w(i,j,k,2))

       enddo
       enddo
       enddo

!#ifdef bidon

!      sum_t_udq_h=sum_t_udq_h+dte_fw*sum6 ! sum u.dq/dx
!      sum_t_udq_v2=sum_t_udq_v2+dte_fw*sum4 ! sum w.dq/dz effectif
!      sum_t_udu=sum_t_udu+dte_fw*sum2     ! sum u.du/dt+w.dw/dt

       if(iteration3d==0.and.iteration2d==1)sum_u2_init=sum3
       sum_u2_now=sum5

       sum8=0.
       sum9=0.
       do j=1,jmax
       do i=1,imax
       do k=1,kmax
       sum8=sum8-q_t(i,j,k,1)*mask_t(i,j,k)*0.5*( &
        omega_w(i,j,k+1,2)-omega_w(i,j,k,2)        &
       +omega_w(i,j,k+1,1)-omega_w(i,j,k,1)      )
!      sum9=sum9-(q_t(i,j,k,2)-q_t(i,j,k,1))*mask_t(i,j,k)*0.5*( &
       sum9=sum9-(w2_*q_t(i,j,k,2)+w1_*q_t(i,j,k,1)+w0_*q_t(i,j,k,0)-q_t(i,j,k,1))*mask_t(i,j,k)*0.5*( &
        omega_w(i,j,k+1,2)-omega_w(i,j,k,2)        &
       +omega_w(i,j,k+1,1)-omega_w(i,j,k,1)      )
       enddo
       enddo
       enddo
       sum_t_udq_err=sum_t_udq_err+dte_fw*sum9 ! sum w.dq/dz erreur effectif-ideal
       sum_t_udq_v1=sum_t_udq_v1+dte_fw*sum8   ! sum w.dq/dz ideal
!      write(6,*)'loulou',real(sum4),real(sum8)+real(sum9)
       write(6,*)'loulou',real(sum_t_udq_v2),real(sum_t_udq_v1)+real(sum_t_udq_err)

       sum3=0.
       do j=1,jmax
       do i=1,imax

         sum3=sum3+mask_t(i,j,kmax)*ssh_w(i,j,1)*grav*0.5*( &
                     +(ssh_w(i,j,2)-ssh_w(i,j,0))/dte_fw    &
                                                          )

       enddo
       enddo
       sum_t_sshdssh=sum_t_sshdssh+dte_fw*sum3
       
       if(iteration3d==0.and.iteration2d==1) then !ooo>
        sum_ssh2_init=0.
        do j=1,jmax ; do i=1,imax
         sum_ssh2_init=sum_ssh2_init+0.5*grav*ssh_w(i,j,0)*ssh_w(i,j,1)*mask_t(i,j,kmax)
        enddo       ; enddo
       endif                                      !ooo>

       sum_ssh2_now=0.
       do j=1,jmax ; do i=1,imax
        sum_ssh2_now=sum_ssh2_now+0.5*grav*ssh_w(i,j,2)*ssh_w(i,j,1)*mask_t(i,j,kmax)
       enddo       ; enddo

       sum5=0.
       sum0=0.
       sum_q2_now=0.
       do j=1,jmax
       do i=1,imax
       do k=1,kmax
         sum5=sum5+mask_t(i,j,k)*0.5*q_t(i,j,k,1) &
      *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
                     *(q_t(i,j,k,2)-q_t(i,j,k,0))*(dtconvert*dtconvert)/dte_fw
!     sum1=sum1+0.5*(q_t(i,j,k,2)+q_t(i,j,k,1)) &
!     *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
!         *(q_t(i,j,k,2)-q_t(i,j,k,1))*(dtconvert*dtconvert)/dte_fw

        sum0=sum0+mask_t(i,j,k)*0.5*q_t(i,j,k,1) &
      *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
                     *(q_t(i,j,k,0))*(dtconvert*dtconvert)
        sum_q2_now=sum_q2_now+mask_t(i,j,k)*0.5*q_t(i,j,k,1) &
      *0.5*(anyv3dr4(i,j,k+1,id_zt)-anyv3dr4(i,j,k-1,id_zt)) & ! dz(k)=(z(k+1)-z(k-1))/2
                     *(q_t(i,j,k,2))*(dtconvert*dtconvert)
       enddo
       enddo
       enddo
       sum_t_qdq=sum_t_qdq+sum5*dte_fw
       if(iteration3d==0.and.iteration2d==1)sum_q2_init=sum0

       write(6,*)'bobi',sum5,sum7+sum8
!      write(6,*)'bobi',real(sum_t_qdq),real(sum_t_udq_h+sum_t_udq_v2)
!      write(6,*)'bobi',real(sum_t_udq_v2),real(sum_t_udq_v1)+real(sum_t_udq_err)

!      write(66,*)sum1,sum3
       write(66,*)real(sum_t_udssh),real(sum_t_sshdssh),real(sum_ssh2_now-sum_ssh2_init)

!      write(67,*)sum2,sum1+sum6+sum4
!      write(67,*)real(sum_t_udu),real(sum_t_udssh+sum_t_udq_h)
!      write(67,*)real(sum_t_udu+sum_t_sshdssh),-real(sum_t_udq_h)-real(sum_t_udq_v2)
!      write(67,*)real(sum_t_udu+sum_t_sshdssh) &
!               ,-real(sum_t_udq_h+sum_t_udq_v1)-real(sum_t_udq_err)
       write(67,'(9(e14.7,1x))')            &
                ,real(sum_t_udu)            & ! C1
                ,real(sum_t_sshdssh)        & ! C2
                ,real(sum_t_qdq)            & ! C3
                ,-real(sum_t_udq_err)       & ! C4=C1+C2+C4 devrait etre zero idealement
                ,sum_u2_now-sum_u2_init     & ! C5=C1
                ,sum_ssh2_now-sum_ssh2_init & ! C6=C2
                ,sum_q2_now-sum_q2_init       ! C7=C3


      end subroutine q_pseudo3d_sigma
#endif

!............................

      subroutine q_diag_ssh2_u2_q2
      implicit none

      if(mod(iteration2d,100)/=0)return
    
! Somme grav*0.5*ssh**2
      sum0=0.
      do j=1,jmax
      do i=1,imax
      sum0=sum0+mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j)*0.5*grav*ssh_w(i,j,1)*ssh_w(i,j,0)
      enddo
      enddo
      call mpi_allreduce(sum0,sum0glb,1,mpi_double,mpi_sum,par%comm2d,ierr)

      sum1=0.
      do j=2,jmax-1
      do i=2,imax
      sum1=sum1+mask_u(i,j,kmax)*mask_i_u(i)*mask_j_u(j)*dxdy_u(i,j) &
      *(h_u(i,j)+0.5*(ssh_w(i-1,j,0)+ssh_w(i,j,0)))*0.5*velbar_u(i,j,1)*velbar_u(i,j,0)
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
      sum1=sum1+mask_v(i,j,kmax)*mask_i_v(i)*mask_j_v(j)*dxdy_v(i,j) &
      *(h_u(i,j)+0.5*(ssh_w(i,j-1,0)+ssh_w(i,j,0)))*0.5*velbar_v(i,j,1)*velbar_v(i,j,0)
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double,mpi_sum,par%comm2d,ierr)

!     if(flag_nh2d==1) then !>>>
!      sum2=0.
!      do k=1,kmax
!      do j=1,jmax
!      do i=1,imax
!      sum2=sum2+mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)*dxdy_t(i,j) &
!        *(h_w(i,j)+ssh_w(i,j,1))*dsigr4_t(i,j,k)                         &
!        *0.5*q_t(i,j,k,1)*q_t(i,j,k,0)
!      enddo
!      enddo
!      enddo
!      call mpi_allreduce(sum2,sum2glb,1,mpi_double,mpi_sum,par%comm2d,ierr)
!     else                  !>>>
       sum2glb=0.
!     endif                 !>>>

      if(par%rank==0) then !>>>
       open(unit=3,file='tmp/nhdiag_ssh2_u2_q2',position='append')
       write(3,*)real(iteration2d*dte_fw/3600.),real(sum0glb),real(sum1glb),real(sum2glb*(dtconvert**2))
       close(3)
      endif                !>>>
!     if(par%rank==0)write(66,*)sum0glb,sum1glb,sum2glb*(dtconvert**2)

! QDM DANS LA DIRECTION OX
      sum0=0.
      do j=2,jmax-1
      do i=2,imax
      sum0=sum0+mask_u(i,j,kmax)*mask_i_u(i)*mask_j_u(j)*dx_u(i,j)*fluxbar_u(i,j,1)
      enddo
      enddo
      call mpi_allreduce(sum0,sum0glb,1,mpi_double,mpi_sum,par%comm2d,ierr)
      if(par%rank==0) then !>>>
       open(unit=3,file='tmp/nhdiag_qdm_oi',position='append')
       write(3,*)real(iteration2d*dte_fw/3600.),real(sum0glb)
       close(3)
      endif                !>>>

      end subroutine q_diag_ssh2_u2_q2

!............................

      subroutine q_moteur_remapping
      implicit none 

      if(checkr0/=anyv3dr4(imax/2,jmax/2,kmax/2,id_dz)) &
         stop 'err nh 1282'

      if(checkr1/=anyv3dr4(imax/2,jmax/2,kmax/2,id_zt))       &
         stop 'err nh 1286'

! Je commente ces lignes car elles sont dans add_nhpgf....
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3dr4(i,j,k,id_zt)=sig1d_t(k)*(h_w(i,j)+ssh_w(i,j,1))+ssh_w(i,j,1)
!     enddo       ; enddo  ; enddo
! Surface et fond
!     do j=1,jmax ; do i=1,imax
!      anyv3dr4(i,j,kmax+1,id_zt)=0.
!      anyv3dr4(i,j,0)=-h_w(i,j)
!     enddo  ; enddo

! rappel: z_t(k+1)-z_t(k)=H*dsig_w(k+1)
! rappel: z_w(k+1)-z_w(k)=H*dsig_t(k)
      
!----------------------------------------------------------------------------------------------
      if(computematrix==1) then !m[°v°]w>

! sshr4_w a EtE calculEe dans subroutine q_nhpgf2d
       ww_=w2_ ; if(ww_==0.)ww_=1.
       do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
        mc(i,j,k,1)=-ww_*dt2/( dsigr4_t(i,j,k)*dsig_w(i,j,k)*(h_w(i,j)+sshr4_w(i,j))**2 )
       enddo       ; enddo       ;enddo

       do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
        mc(i,j,k,3)=mc(i,j,k,1)*dsig_w(i,j,k)/dsig_w(i,j,k+1)
       enddo       ; enddo       ;enddo

       do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
        mc(i,j,k,2)=1.-mc(i,j,k,1)-mc(i,j,k,3)
       enddo       ; enddo       ;enddo

!.......................................................
! C.L. Surface: q_t(ktop)=0.
       k=ktop
       do j=2,jmax-1 ; do i=2,imax-1
        mc(i,j,k,1)=0.
        mc(i,j,k,2)=1.
        mc(i,j,k,3)=0.
       enddo       ; enddo

!.......................................................
! C.L. Fond: q_t(0)=q_t(1)
       k=1
       do j=2,jmax-1 ; do i=2,imax-1
        mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
        mc(i,j,k,1)=0.
       enddo       ; enddo

       checkr2=mc(imax/2,jmax/2,kmax/2,1)
       checkr3=mc(imax/2,jmax/2,kmax/2,2)
       checkr4=mc(imax/2,jmax/2,kmax/2,3)

       if(nh_frozensigma==1)computematrix=0 ! A la prochaine iteration on ne passera plus par ces lignes

      endif                     !m[°v°]w>
!----------------------------------------------------------------------------------------------

      if(checkr2/=mc(imax/2,jmax/2,kmax/2,1)) &
      stop 'checkr2/=mc(imax/2,jmax/2,kmax/2,1)'
      if(checkr3/=mc(imax/2,jmax/2,kmax/2,2)) &
      stop 'checkr3/=mc(imax/2,jmax/2,kmax/2,2)'
      if(checkr4/=mc(imax/2,jmax/2,kmax/2,3)) &
      stop 'checkr4/=mc(imax/2,jmax/2,kmax/2,3)'

!.......................................................
! C.L. Surface mc4
       k=ktop
       do j=2,jmax-1 ; do i=2,imax-1
        mc(i,j,k,4)=0.
       enddo       ; enddo

      ww_=1.
      if(w2_/=0.)ww_=1./w2_

      do j=2,jmax-1 ; do i=2,imax-1

       anyvar2d(i,j)=dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
       ((fluxbar_u(i+1,j,0)-fluxbar_u(i+1,j,1))/(dte_fw*dy_u(i+1,j))-nhpgf2d_u(i+1,j))/(h_u(i+1,j)+0.5*(sshr4_w(i+1,j)+sshr4_w(i,j))) &
      -((fluxbar_u(i  ,j,0)-fluxbar_u(i  ,j,1))/(dte_fw*dy_u(i  ,j))-nhpgf2d_u(i  ,j))/(h_u(i  ,j)+0.5*(sshr4_w(i-1,j)+sshr4_w(i,j))) &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
       ((fluxbar_v(i,j+1,0)-fluxbar_v(i,j+1,1))/(dte_fw*dx_v(i,j+1))-nhpgf2d_v(i,j+1))/(h_v(i,j+1)+0.5*(sshr4_w(i,j+1)+sshr4_w(i,j))) &
      -((fluxbar_v(i,j  ,0)-fluxbar_v(i,j  ,1))/(dte_fw*dx_v(i,j  ))-nhpgf2d_v(i,j  ))/(h_v(i,j  )+0.5*(sshr4_w(i,j-1)+sshr4_w(i,j))) &
                    ) & !ooo>

                         )    !ttt>

      enddo ; enddo

      do k=1,ktop-1 
      do j=2,jmax-1 ; do i=2,imax-1

!      if(i==imax/2.and.j==jmax/2)write(6,*)'2D',k,anyvar2d(i,j),velbar_u(i+1,j,1:2)

        mc(i,j,k,4)=2*q_t(i,j,k,1)-q_t(i,j,k,0)                        &

             -(q_t(i,j,k+1,1)-     q_t(i,j,k-1,1))                     &!b3-time
        /(anyv3dr4(i,j,k+1  ,id_zt)-anyv3dr4(i,j,k-1  ,id_zt))                     &!b3-time
        *(1.+sig1dpom_t(i,j,k))*(2*ssh_w(i,j,1)-ssh_w(i,j,0)-ssh_w(i,j,2)) &!b3-time

                +anyvar2d(i,j)   &
                                              +dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>

                ( q_t(i+1,j,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!               +(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b1
!          /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b1
!      *0.5*(anyv3dr4(i-1,j,k    ,id_zt)-anyv3dr4(i+1,j,k    ,id_zt)) & !b1

!               +(q_t(i+1,j,k+1,1)-     q_t(i+1,j,k-1,1)) & !b2
!          /(anyv3dr4(i+1,j,k+1  ,id_zt)-anyv3dr4(i+1,j,k-1  ,id_zt)) & !b2
!          *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i+1,j,k    ,id_zt)) & !b2

                +(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b3
           /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b3
           *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i+1,j,k    ,id_zt)) & !b3

                                -q_t(i  ,j,k,1) )*invdx_u(i+1,j)    &

               -( q_t(i  ,j,k,1)-q_t(i-1,j,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!                              -(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b1
!                         /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b1
!                     *0.5*(anyv3dr4(i+1,j,k    ,id_zt)-anyv3dr4(i-1,j,k    ,id_zt)) & !b1

!                              -(q_t(i-1,j,k+1,1)-     q_t(i-1,j,k-1,1)) & !b2
!                         /(anyv3dr4(i-1,j,k+1  ,id_zt)-anyv3dr4(i-1,j,k-1  ,id_zt)) & !b2
!                         *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i-1,j,k    ,id_zt)) & !b2

                               -(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b3
                          /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b3
                          *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i-1,j,k    ,id_zt)) & !b3

                                                )*invdx_u(i  ,j)    &

                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>

                ( q_t(i,j+1,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!               +(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b1
!          /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b1
!      *0.5*(anyv3dr4(i,j-1,k    ,id_zt)-anyv3dr4(i,j+1,k    ,id_zt)) & !b1

!               +(q_t(i,j+1,k+1,1)-     q_t(i,j+1,k-1,1)) & !b2
!          /(anyv3dr4(i,j+1,k+1  ,id_zt)-anyv3dr4(i,j+1,k-1  ,id_zt)) & !b2
!          *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j+1,k    ,id_zt)) & !b2

                +(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b3
           /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b3
           *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j+1,k    ,id_zt)) & !b3

                                -q_t(i  ,j,k,1) )*invdy_v(i,j+1)    &

               -( q_t(i  ,j,k,1)-q_t(i,j-1,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!                              -(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b1
!                         /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b1
!                     *0.5*(anyv3dr4(i,j+1,k    ,id_zt)-anyv3dr4(i,j-1,k    ,id_zt)) & !b1

!                              -(q_t(i,j-1,k+1,1)-     q_t(i,j-1,k-1,1)) & !b2
!                         /(anyv3dr4(i,j-1,k+1  ,id_zt)-anyv3dr4(i,j-1,k-1  ,id_zt)) & !b2
!                         *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j-1,k    ,id_zt)) & !b2

                               -(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b3
                          /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b3
                          *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j-1,k    ,id_zt)) & !b3

                                                )*invdy_v(i  ,j)    &

                    ) & !ooo>

                                                   ) &  !ttt>


       -mc(i,j,k,1)    *ww_*(w0_*q_t(i,j,k-1,0)+w1_*q_t(i,j,k-1,1) )&
       -mc(i,j,k,3)    *ww_*(w0_*q_t(i,j,k+1,0)+w1_*q_t(i,j,k+1,1) )&
      -(mc(i,j,k,2)-1.)*ww_*(w0_*q_t(i,j,k  ,0)+w1_*q_t(i,j,k  ,1) ) 

! Note sur sponge_t: varie de 0 a 1 (sur la frontiere) ce qui equivaut en termes
! d'echelle de temps associee a un temps de rappel equivalent au pas de temps. C'est
! donc un rappel fort qui a pour but de supprimer la pression NH aux frontieres ouvertes

      enddo ; enddo
      enddo

! Solveur vertical tridiagonal:
      call q_solver(computefullsolver,2,imax-1,2,jmax-1,1,ktop) ! donne q_t(t+1)
      if(nh_frozensigma==1)computefullsolver=0 ! A la prochaine iteration une partie seulement du solver sera calculEe

!........................
! ZONES EPONGE ET WETMASK: (Ne pourrait elle pas etre implicitee ou simplifiee?)
      do j=2,jmax-1 ; do i=2,imax-1
       wetmask_t(i,j)=min(1.,(h_w(i,j)+sshr4_w(i,j))/wetdry_cstnh)
      enddo         ; enddo
!     do k=0,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1
      do k=1,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! Boucle k commence A 1 si report des C.L. en suivant....
!      q_t(i,j,k,2)=q_t(i,j,k,2)*(1.-sponge_t(i,j,1))*wetmask_t(i,j)
       q_t(i,j,k,2)=q_t(i,j,k,2)                     *wetmask_t(i,j)
       if(k==kmax/2.and.j==jmax/2.and.i==imax/2)write(66,*)iteration2d,q_t(i,j,k,2)
      enddo         ; enddo       ; enddo


!........................
! Conditions ouvertes radiatives
      call q_obc_radiative

!........................
! Reporter les hypotheses de conditions aux limites au fond et en surface:
      do j=2,jmax-1 ; do i=2,imax-1
! C.L.: q(0)=q(1):
       q_t(i,j,0,2)=q_t(i,j,1,2)*mask_t(i,j,kmax)
! C.L.: q(kmax+1)=0.
       q_t(i,j,ktop:kmax+1,2)=0.
      enddo ; enddo

!........................
! Continuite MPI:
      call q_obc_mpi('z0') ! mpi 'z1' sur q_t(t+1)

!#ifdef bidon
! VERIFICATION DE L'EQUILIBRE:
      if(modulo(iteration2d,100)==0) then !>>>
! Ici verification possible de l'equilibre
! Attention que le filtre d'asselin se fait ulterieurement
! et par consequent n'est pas evalue ici
      i=imax/2 ; j=jmax/2 ; k=kmax/2
      write(texte30,'(a,i0)')'tmp/nhdiag_balance',par%rank
! opposE du Terme de time stepping:
       x1= -q_t(i,j,k,2)+2*q_t(i,j,k,1)-q_t(i,j,k,0)                  &
             -(q_t(i,j,k+1,1)-     q_t(i,j,k-1,1))                     &!b3-time
        /(anyv3dr4(i,j,k+1  ,id_zt)-anyv3dr4(i,j,k-1  ,id_zt))                     &!b3-time
        *(1.+sig1dpom_t(i,j,k))*(2*ssh_w(i,j,1)-ssh_w(i,j,0)-ssh_w(i,j,2))  !b3-time
! Forcage hydrostatique:
       x2=anyvar2d(i,j)
! Laplacien 3D de la pression NH deduit de l'equilibre elementaire
! m1*q(k-1,t+1)+m2*q(k  ,t+1)+m3*q(k+1,t+1)=m4 que l'on a pris
! soin prealablement de verifier
       x3=-q_t(i,j,k,2)+mc(i,j,k,4)-x1-x2       &
                 -mc(i,j,k,1)   *q_t(i,j,k-1,2) &
                -(mc(i,j,k,2)-1)*q_t(i,j,k  ,2) &
                 -mc(i,j,k,3)   *q_t(i,j,k+1,2)  

      open(unit=3,file=texte30,position='append')
! Equilibre elementaire qui doit etre verifiE:
!      write(3,*)mc(i,j,k,1)*q_t(i,j,k-1,2) &
!               +mc(i,j,k,2)*q_t(i,j,k  ,2) &
!               +mc(i,j,k,3)*q_t(i,j,k+1,2) &
!               ,mc(i,j,k,4)
! l'addition de 2+3+4 egale zero:
       write(3,*)real(iteration2d*dte_fw/3600.) & ! 1 temps en heure
                ,real(x1)                       & ! 2 terme de time-stepping (ou desequilibre)
                ,real(x2)                       & ! 3 forcage hydrostatique
                ,real(x3)                         ! 4 Laplacien 3D non-hydrostatique
! Si le solveur a bien convergE, 2 doit etre petit par rapport a 3 et 4
      close(3)
      endif                               !>>>
!#endif

      end subroutine q_moteur_remapping

!..............................

      subroutine q_wind_induced_waves
      implicit none

!     if(iteration2d/=1.and.                            &
!        check_sshwave/=xy_t(imax/2,jmax/2,id_sshwave)) &
!        stop 'Err check_sshwave'

! Details dans:
! https://docs.google.com/document/d/1IcV7uNcz4JW-jLP3c2WJp0Yx3_E7XBdQ-y3uzhYUEdM/edit

! Donner une valeur du vent:
! Equation 7 de Rolf Deigaard , Peter Nielsen, Coastal Engineering 139 (2018) 36–46
! https://doi.org/10.1016/j.coastaleng.2018.04.024
! P+=rho*grav*(-2)*beta/w*A*sin(kx-wt)=rho*grav*2*beta/(w*k)*Dssh/Dx
! Avec beta<<w

! COEF POUR methode 2
!      x1_r4=0.27*rhoair/rho
! COEF POUR methode 1:
!      x1_r4=0.27*rhoair/rho &
!           *0.5 !Pour equivalence avec methode 2, multiplier par V/Vrelatif...
! COEF POUR methode 3: 
! (Snyder et al, 1981) voir (60) dans
!https://docs.google.com/document/d/1IcV7uNcz4JW-jLP3c2WJp0Yx3_E7XBdQ-y3uzhYUEdM/edit?usp=sharing

       x2_r4=1./dti_lp

#ifdef bidon
       x1_r4=grav*rhoair/rho &
                 *8.8        & ! periode peak
                 /pi         &
                 *0.25       &
                 *40.    ! facteur pour contrer Asselin
!      x1_r4=0.

! Ajouter un dephasage et une dependance en 1/freq
!      x3_r4=1.
!      x3_r4=dti_fw/1. ; x1_r4=x1_r4/0.66 ! supprimer les periodes plus basses que 1    secondes
       x3_r4=dti_fw/8. ! supprimer les periodes plus basses que 3 secondes
#endif


#ifdef bidon
! Retards et facteurs d'amplification des pressions composites:
       if(iteration3d==0.or.iteration3d==iteration3d_restart) then !m°v°m>
        do loop1=1,retarded_pss_nb
! Retards (dt_over_retartedtime):
         if(loop1==1              )dt_over_retartedtime(loop1)=dti_fw/1.
         if(loop1==2              )dt_over_retartedtime(loop1)=dti_fw/2.
         if(loop1==3              )dt_over_retartedtime(loop1)=dti_fw/4.
         if(loop1==retarded_pss_nb)dt_over_retartedtime(loop1)=dti_fw/8.
! Amplifications (retarded_pss_amp):
         retarded_pss_amp(loop1)=grav*rhoair/rho*8.8/pi*0.25 &
                                *0.3 &
!                               /(dt_over_retartedtime(loop1)/dti_fw)**2
                                *5.

        enddo !loop1 
       endif                                                       !m°v°m>
#endif

       do j=1,jmax ; do i=1,imax

! Perturbation de pression de surface

#ifdef bidon
! Methode 1:
            mc(i,j,ktop,4)=( & !ooo>
            sqrt( uwind_t(i,j,1)**2 +vwind_t(i,j,1)**2)*             &
                           x1_r4*( & !pmx>

       0.5*(uwind_t(i,j,1)*( ssh_int_w(i+1,j,1)                      &
                            -ssh_int_w(i-1,j,1))*invdx_t(i,j)        &
           +vwind_t(i,j,1)*( ssh_int_w(i,j+1,1)                      &
                            -ssh_int_w(i,j-1,1))*invdy_t(i,j))       &

           +(ssh_int_w(i,j,2)-ssh_int_w(i,j,0))*x2_r4                &
!          +(ssh_int_w(i+1,j,2)-ssh_int_w(i+1,j,0) &
!           +ssh_int_w(i-1,j,2)-ssh_int_w(i-1,j,0))*x2_r4*0.5        &


                                 ) & !pmx>                             
                           ) & !ooo>
          *x3_r4+(1.-x3_r4)*mc(i,j,ktop,4) ! filtre BF
#endif


#ifdef bidon
! Methode 2:
            mc(i,j,ktop,4)=( & !ooo>

            sqrt( & !SQRT>
                  (     & !--->(Uwind-Cx)**2
                   uwind_t(i,j,1)   &
                 +(ssh_int_w(i,j,2)-ssh_int_w(i,j,0))*x2_r4                &
             *0.5*(ssh_int_w(i+1,j,1)-ssh_int_w(i-1,j,1))*invdx_t(i,j)     &
        /(   (0.5*(ssh_int_w(i+1,j,1)-ssh_int_w(i-1,j,1))*invdx_t(i,j))**2 &
            +(0.5*(ssh_int_w(i,j+1,1)-ssh_int_w(i,j-1,1))*invdy_t(i,j))**2 &
            +small1) &
                  )**2 &  !--->(Uwind-Cx)**2

                 +(     & !--->(Vwind-Cy)**2
                   vwind_t(i,j,1)   &
                 +(ssh_int_w(i,j,2)-ssh_int_w(i,j,0))*x2_r4                &
             *0.5*(ssh_int_w(i,j+1,1)-ssh_int_w(i,j-1,1))*invdy_t(i,j)     &
        /(   (0.5*(ssh_int_w(i+1,j,1)-ssh_int_w(i-1,j,1))*invdx_t(i,j))**2 &
            +(0.5*(ssh_int_w(i,j+1,1)-ssh_int_w(i,j-1,1))*invdy_t(i,j))**2 &
            +small1) &
                  )**2 &  !--->(Vwind-Cy)**2

                ) & !SQRT>


                          *x1_r4*( & !pmx>

       0.5*(uwind_t(i,j,1)*( ssh_int_w(i+1,j,1)                      &
                            -ssh_int_w(i-1,j,1))*invdx_t(i,j)        &
           +vwind_t(i,j,1)*( ssh_int_w(i,j+1,1)                      &
                            -ssh_int_w(i,j-1,1))*invdy_t(i,j))       &

           +(ssh_int_w(i,j,2)-ssh_int_w(i,j,0))*x2_r4                &


                                 ) & !pmx>
                           ) & !ooo>
          *x3_r4+(1.-x3_r4)*mc(i,j,ktop,4) ! filtre BF
#endif

#ifdef bidon
! Methode 3 "composite": ! Snyder et al 1981. Attention le vent est supposE A 5m!
                  anyvar2d(i,j)=rho*( & !ooo>
!!       retarded_pss_w(i,j,1) =rho*( & !ooo>
                                 ( & !pmx>

       0.5*(uwind_t(i,j,1)*( ssh_int_w(i+1,j,1)                      &
                            -ssh_int_w(i-1,j,1))*invdx_t(i,j)        &
           +vwind_t(i,j,1)*( ssh_int_w(i,j+1,1)                      &
                            -ssh_int_w(i,j-1,1))*invdy_t(i,j))       &

           +(ssh_int_w(i+1,j,2)-ssh_int_w(i+1,j,0) &
            +ssh_int_w(i-1,j,2)-ssh_int_w(i-1,j,0))*x2_r4*0.5        &


                                 ) & !pmx>                             
                           )   !ooo>
!!        *dt_over_retartedtime(1)+(1.-dt_over_retartedtime(1))*retarded_pss_w(i,j,1) ! filtre BF
#endif

       enddo ; enddo

#ifdef bidon
! Pressions retardees
!!     do loop1=2,retarded_pss_nb
       do loop1=1,retarded_pss_nb
! dt_over_retartedtime(loop1) est dti_fw/t1 oU t1 est le retard de retarded_pss_w(i,j,loop1)
! par rapport A retarded_pss_w(i,j,1), sachant que retarded_pss_w(i,j,1)
! a dEjA un retard par rapport A la pente instantanEe de la SSH 
! (voir valeur de x3_r4 appliquee A retarded_pss_w(i,j,1))
        do j=1,jmax ; do i=1,imax
         retarded_pss_w(i,j,loop1)=                                  &
                         anyvar2d(i,j)*dt_over_retartedtime(loop1)   &
!!       retarded_pss_w(i,j,1)        *dt_over_retartedtime(loop1)   &
        +retarded_pss_w(i,j,loop1)*(1.-dt_over_retartedtime(loop1))
        enddo       ; enddo        ! i,j
       enddo ! loop1
! La pression pss_w totale (celle en pratique utilisee par le modele)
! est un composite des differentes pressions retardEes
       k0=0
       do loop1=1,retarded_pss_nb
        do j=1,jmax ; do i=1,imax
                  pss_w(i,j,1)=     &
                  pss_w(i,j,1)*k0   &
        +retarded_pss_amp(loop1)*retarded_pss_w(i,j,loop1)
        enddo       ; enddo        ! i,j
        k0=1
       enddo ! loop1

       if(par%rank==0) then
        i=imax/2 ; j=jmax/2
        write(10+par%rank,'(f12.3,1x,5(1x,e14.7),a)') &
             elapsedtime_now  &
            ,anyvar2d(i,j) &
            ,retarded_pss_amp(1)*retarded_pss_w(i,j,1) &
            ,retarded_pss_amp(2)*retarded_pss_w(i,j,2) &
            ,retarded_pss_amp(3)*retarded_pss_w(i,j,3) &
            ,retarded_pss_amp(4)*retarded_pss_w(i,j,4) &
            ,' trouvemoi64' 
       endif
#endif

! Partie basee sur le windsress

! SSH "basse frequence"
!     x0_r4=dti_fw/60. ! Pour la methode "pente ssh"
      x0_r4=dti_fw/8.  ! Pour la methode "ssh-sshlwf"
      do j=1,jmax ; do i=1,imax
       sshlwf_w(i,j,1)=(1.-x0_r4)*sshlwf_w(i,j,1)+x0_r4*ssh_int_w(i,j,1)
      enddo       ; enddo

! Pour la methode "pente ssh"
!     x0_r4=2.*grav*rho*hmax/sqrt(grav*hmax)*3.5e-4
!     x1_r4=dti_fw/0.5 ! filtrage BF wstress_u

! Pour la methode "ssh-sshlwf"
      x0_r4=2.*grav*rho*hmax/sqrt(grav*hmax)*0.55e-5

      do j=1,jmax ; do i=2,imax

#ifdef bidon
! stress=x0*( Uwind-c)*|dsshdx|*sign(ssh)

! methode "pente ssh"
       wstress_u(i,j,1)=  & 

! Parce que la fonction sign est discontinue le filtrage via x1_r4 permet une transition progressive de wstress
       wstress_u(i,j,1)*(1.-x1_r4)+x1_r4* &
                        ( & 

         0.5*(uwind_t(i-1,j,1)+uwind_t(i,j,1)) & ! WIND
         *abs((ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1))*invdx_u(i  ,j))    &


           +( ssh_int_w(i  ,j,2)-ssh_int_w(i  ,j,0)                      &
             +ssh_int_w(i-1,j,2)-ssh_int_w(i-1,j,0))*0.5*inv_dti_lp      &
           *( ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1))                     &
       /(abs((ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1)))+1.e-4 ) &   

                        ) & 
! Methode 1: depend de ssh qui, via dependance A |dsshdx|, entraine une dependance quadratique A la ssh 
! qui entraine un risque d'emballement (A eviter par consequent)
!                        *x0_r4*0.5*(ssh_int_w(i  ,j,1)+ssh_int_w(i-1,j,1) &  
!                                    -sshlwf_w(i  ,j,1) -sshlwf_w(i-1,j,1))   
! Methode 2: depend de sign(1,ssh) pour casser la dependance quadratique A la ssh
         *x0_r4*sign(1., ssh_int_w(i  ,j,1)+ssh_int_w(i-1,j,1)  & 
                         -sshlwf_w(i  ,j,1) -sshlwf_w(i-1,j,1))   

#endif

!#ifdef bidon
! methode "ssh-sshlwf"

       wstress_u(i,j,1)=( & ! VENT moins VITESSE DE PHASE  >

         0.5*(uwind_t(i-1,j,1)+uwind_t(i,j,1)) & ! WIND

       +( & !ooo>  ! Vitesse phase DEBUT >

         ( ssh_int_w(i  ,j,2)-ssh_int_w(i  ,j,0)                      &
          +ssh_int_w(i-1,j,2)-ssh_int_w(i-1,j,0))*0.5*inv_dti_lp      &
        *( ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1))*invdx_u(i,j)        &

        +( ssh_int_w(i+1,j,2)-ssh_int_w(i+1,j,0)                      &
          +ssh_int_w(i  ,j,2)-ssh_int_w(i  ,j,0))*0.5*inv_dti_lp      &
        *( ssh_int_w(i+1,j,1)-ssh_int_w(i  ,j,1))*invdx_u(i+1,j)      &

        +( ssh_int_w(i-1,j,2)-ssh_int_w(i-1,j,0)                      &
          +ssh_int_w(i-2,j,2)-ssh_int_w(i-2,j,0))*0.5*inv_dti_lp      &
        *( ssh_int_w(i-1,j,1)-ssh_int_w(i-2,j,1))*invdx_u(i-1,j)      &

        ) & !ooo>
       /( & !pmx>
          ((ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1))*invdx_u(i  ,j))**2 &
         +((ssh_int_w(i+1,j,1)-ssh_int_w(i  ,j,1))*invdx_u(i+1,j))**2 &
         +((ssh_int_w(i-1,j,1)-ssh_int_w(i-2,j,1))*invdx_u(i-1,j))**2 &
         +1.e-8                                                       &
        ) & !pmx>  ! Vitesse phase FIN >

                        ) & ! VENT moins VITESSE DE PHASE  >
                         *x0_r4*0.5*(ssh_int_w(i  ,j,1)+ssh_int_w(i-1,j,1) &  
                                     -sshlwf_w(i  ,j,1) -sshlwf_w(i-1,j,1))   
!#endif

!      if(i==imax/2.and.j==jmax/2.and.par%rank==0) &
!       write(10+par%rank,'(f12.3,1x,3(1x,e14.7),a)') &
!      elapsedtime_now     &
!      ,wstress_u(i,j,1)   &
!      , 0.5*(uwind_t(i-1,j,1)+uwind_t(i,j,1)) & ! WIND
!        *abs((ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1))*invdx_u(i  ,j))    &


!      ,   +( ssh_int_w(i  ,j,2)-ssh_int_w(i  ,j,0)                      &
!            +ssh_int_w(i-1,j,2)-ssh_int_w(i-1,j,0))*0.5*inv_dti_lp      &
!          *( ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1))                     &
!      /(abs((ssh_int_w(i  ,j,1)-ssh_int_w(i-1,j,1)))+1.e-4 ) &   
!      ,ssh_int_w(i,j,1)   &
!      ,sshlwf_w(i,j,1)    &   
!      ,breaker2d_t(i,j)   &
!      ,'trouvemi'
      
      enddo       ; enddo

!......................................
! BRUIT A 5s sur le point i=imax/2:
!     if(modulo(iteration3d,196)==0) then 
!      i=imax/2 ; j=jmax/2
!      wstress_u(i,:,1)= &
!      wstress_u(i,:,1)  &
!       +0.5*(uwind_t(i-1,j,1)+uwind_t(i,j,1)) & ! WIND
!                        *x0_r4*10.
!     endif
!     if(par%rank==0) then
!       i=imax/2 ; j=jmax/2
!       write(10+par%rank,'(f12.3,1x,2(1x,e14.7),a)') &
!       elapsedtime_now     &
!      ,wstress_u(i,j,1)   &
!      ,breaker2d_t(i,j)   &
!      ,'trouvemoa'
!     endif
!......................................
 
      if(jmax>10) then
! Attention, dans le cas directionnel, faut il prendre la projection
! de C sur les axes Oi et Oj ou la vitesse de phase apparente (la
! premiere tend vers zero quand la seconde tend vers l'infini) ?
! A ce stade la reponse n'est pas si claire bien que je pencherais
! comme dans SWACH pour la seconde
           stop 'direction Oj reste A faire'
      endif

      stop 'PASSE PAS PAR LA!!!!'

      end subroutine q_wind_induced_waves

!..............................

      subroutine q_time_averaged
      implicit none
     
      stop 'PASSE PAS PAR q_time_averaged stp!'

!     if(.not.allocated(ssh_timeaveraged)) then
!      allocate(variance_timeaveraged(imax  ,jmax  ))        ; variance_timeaveraged=0.
!      allocate(     ssh_timeaveraged(imax  ,jmax  ))        ; ssh_timeaveraged=0.
!      allocate(     flx2d_timeaveraged(imax+1,jmax  ))      ; flx2d_timeaveraged=0. 
!      allocate(     flx3d_timeaveraged(imax+1,jmax  ,kmax)) ; flx3d_timeaveraged=0. 
!      allocate(   u_euler_timeaveraged(imax+1,jmax  ,kmax)) ; u_euler_timeaveraged=0.
!      allocate(     fly2d_timeaveraged(imax  ,jmax+1))      ; fly2d_timeaveraged=0. 
!      allocate(     fly3d_timeaveraged(imax  ,jmax+1,kmax)) ; fly3d_timeaveraged=0. 
!      allocate(   v_euler_timeaveraged(imax  ,jmax+1,kmax)) ; v_euler_timeaveraged=0.
!      allocate(             sshmax_w(imax  ,jmax))   ; sshmax_w=-999.
!      allocate(             sshmin_w(imax  ,jmax))   ; sshmin_w=+999.
!     endif

!      if(flag_nh2d==1)   k10=iteration2d
!      if(flag_nh3d_uv>=1)k10=iteration3d !17_08-18 (flag_nh3d_uv>1 envisagE)
       k10=iteration3d !23-08-18
       if( modulo(k10,nh2d_graph_period)==1 ) then !PMX>
              ssh_timeaveraged=0.
              flx2d_timeaveraged=0.
              fly2d_timeaveraged=0.
              flx3d_timeaveraged=0.
              fly3d_timeaveraged=0.
            u_euler_timeaveraged=0.
            v_euler_timeaveraged=0.
           variance_timeaveraged=0.
                            hs_w=0.
                      sum_avr_hs=0
                      sshmax_w=-999.
                      sshmin_w= 999.
       endif                                       !PMX>

       do j=1,jmax     ;  do i=1,imax
        ssh_timeaveraged(i,j)=ssh_timeaveraged(i,j)+ssh_int_w(i,j,1)
       enddo           ;  enddo
       do j=1,jmax     ;  do i=1,imax
        variance_timeaveraged(i,j)=variance_timeaveraged(i,j)+ssh_int_w(i,j,1)*ssh_int_w(i,j,1)
       enddo           ;  enddo
! Min Max
       x0_r4=1.-dti_fw/10000. ! pour analyser seulement la periode recente
       do j=1,jmax     ;  do i=1,imax
        sshmax_w(i,j)=max(sshmax_w(i,j)*x0_r4,ssh_int_w(i,j,1))  
        sshmin_w(i,j)=min(sshmin_w(i,j)*x0_r4           &
! Ces lignes ont pour but que le minimum puisse etre au dessus du niveau z=0m (en cas de surcote par ex)
!                         +(1.-x0_r4)*max(0.,-h_w(i,j)) &
!                         +(1.-x0_r4)*sshmax_w(i,j)     &
                          +(1.-x0_r4)*ssh_int_w(i,j,1)  &
                          ,ssh_int_w(i,j,1))  
       enddo           ; enddo

       do j=1,jmax     ;  do i=1,imax+1
        flx2d_timeaveraged(i,j)=flx2d_timeaveraged(i,j)+fluxbar_u(i,j,1)
       enddo           ;  enddo
       do k=1,kmax ; do j=1,jmax     ;  do i=1,imax+1
          flx3d_timeaveraged(i,j,k)=  flx3d_timeaveraged(i,j,k)+veldydz_u(i,j,k,1)
       enddo       ; enddo           ;  enddo
       do j=1,jmax+1   ;  do i=1,imax
        fly2d_timeaveraged(i,j)=fly2d_timeaveraged(i,j)+fluxbar_v(i,j,1)
       enddo           ;  enddo
       do k=1,kmax ; do j=1,jmax+1   ;  do i=1,imax
          fly3d_timeaveraged(i,j,k)=  fly3d_timeaveraged(i,j,k)+veldxdz_v(i,j,k,1)
       enddo       ; enddo           ;  enddo

! Moyenne du courant vel_u sur des niveaux fixes
       do j=1,jmax     ;  do i=1,imax+1
        k0=1
        do k=1,kmax
! Profondeur du niveau fixe (c.a.d. pour ssh=0):
         x1=min(max(sig1dpom_t(i,j,k)*h_u(i,j),depth_u(i,j,1)),depth_u(i,j,kmax)-0.00001)
! Chercher les niveaux qui encadrent
         do k1=k0,kmax-1
          if(depth_u(i,j,k1)<=x1.and.depth_u(i,j,k1+1)> x1) then !pmx>
            rap=(               x1-depth_u(i,j,k1)) &
               /(depth_u(i,j,k1+1)-depth_u(i,j,k1))  
            u_euler_timeaveraged(i,j,k)= &
            u_euler_timeaveraged(i,j,k)+rap *vel_u(i,j,k1+1,1) &
                                   +(1.-rap)*vel_u(i,j,k1  ,1)  
            k0=k1 ! comme la boucle est croissante le prochain point sera
                  ! forcement >=k1
          endif                                                  !pmx>
         enddo ! k1
        enddo ! k
       enddo           ;  enddo

! Moyenne du courant vel_v sur des niveaux fixes
       do j=1,jmax+1    ;  do i=1,imax
        k0=1
        do k=1,kmax
! Profondeur du niveau fixe (c.a.d. pour ssh=0):
         x1=min(max(sig1dpom_t(i,j,k)*h_v(i,j),depth_v(i,j,1)),depth_v(i,j,kmax)-0.00001)
! Chercher les niveaux qui encadrent
         do k1=k0,kmax-1
          if(depth_v(i,j,k1)<=x1.and.depth_v(i,j,k1+1)> x1) then !pmx>
            rap=(               x1-depth_v(i,j,k1)) &
               /(depth_v(i,j,k1+1)-depth_v(i,j,k1))
            v_euler_timeaveraged(i,j,k)= &
            v_euler_timeaveraged(i,j,k)+rap *vel_v(i,j,k1+1,1) &
                                   +(1.-rap)*vel_v(i,j,k1  ,1)
            k0=k1 ! comme la boucle est croissante le prochain point sera
                  ! forcement >=k1
          endif                                                  !pmx>
         enddo ! k1
        enddo ! k
       enddo           ;  enddo


       k10=iteration3d !23_08-18

       if( modulo(k10,nh2d_graph_period)==0 ) then !PMX>
! Fin de cumul (qui sera suivi d'archivage dans fichier netcdf)
        x0=1./real(nh2d_graph_period)
        do j=1,jmax     ;  do i=1,imax
         ssh_timeaveraged(i,j)=ssh_timeaveraged(i,j)*x0
        enddo           ;  enddo
        do j=1,jmax     ;  do i=1,imax+1
         flx2d_timeaveraged(i,j)=flx2d_timeaveraged(i,j)*x0/dy_u(i,j)
        enddo           ;  enddo
        do k=1,kmax ; do j=1,jmax     ;  do i=1,imax+1
          flx3d_timeaveraged(i,j,k)=  flx3d_timeaveraged(i,j,k)*x0/dy_u(i,j)
        u_euler_timeaveraged(i,j,k)=u_euler_timeaveraged(i,j,k)*x0
        enddo       ; enddo           ;  enddo
        do j=1,jmax+1   ;  do i=1,imax
         fly2d_timeaveraged(i,j)=fly2d_timeaveraged(i,j)*x0/dx_v(i,j)
        enddo           ;  enddo
        do k=1,kmax ; do j=1,jmax+1   ;  do i=1,imax
          fly3d_timeaveraged(i,j,k)=  fly3d_timeaveraged(i,j,k)*x0/dx_v(i,j)
        v_euler_timeaveraged(i,j,k)=v_euler_timeaveraged(i,j,k)*x0
        enddo       ; enddo           ;  enddo
        do j=1,jmax   ;  do i=1,imax
         variance_timeaveraged(i,j)=sqrt(x0*variance_timeaveraged(i,j)-ssh_timeaveraged(i,j)*ssh_timeaveraged(i,j))*4 !attention il s'agit de Hs
        enddo           ;  enddo
       endif                                       !PMX>

      end subroutine q_time_averaged

!..............................


      subroutine q_obc_radiative
      implicit none

! Frontiere en i=1 ou i=imax
      do loop1=1,2 ! Boucles sur 2 frontieres: (i=1 puis i=imax)

      if( (obcstatus(ieq1)   ==1.and.loop1==1).or.          &
          (obcstatus(ieqimax)==1.and.loop1==2)     )   then ! m[0v0]m

       obcnh_scheme=17


       if(loop1==1) then !obc i=1>
         i1=1    ; i3=i1+1 ; i5=i3+1 ; i7=i5+1 ! pour ssh et q
         i2=2    ; i4=i2+1 ; i6=i4+1 ; i8=i6+1 ! pour v
       endif           !obc i=1>
       if(loop1==2) then !obc i=imax>
         i1=imax ; i3=i1-1 ; i5=i3-1 ; i7=i5-1 ! pour ssh et q
         i2=imax ; i4=i2-1 ; i6=i4-1 ; i8=i6-1 ! pour v
       endif           !obc i=imax>

#ifdef bidon
      do k=1,ktop-1 ; do j=2,jmax-1
                  ema1_q_j(j,k,loop1)=    &
      (1.-ema_mu)*ema1_q_j(j,k,loop1) &
         +ema_mu *abs( (+q_t(i3,j,k,2)-qwave_j_w(j,k,2,2,loop1)  &
                        -q_t(i3,j,k,1)+qwave_j_w(j,k,2,1,loop1)  &
                        -q_t(i5,j,k,1)+qwave_j_w(j,k,3,1,loop1)  &
                        +q_t(i5,j,k,2)-qwave_j_w(j,k,3,2,loop1)) &
                      *(-q_t(i3,j,k,2)+qwave_j_w(j,k,2,2,loop1)  &
                        -q_t(i3,j,k,1)+qwave_j_w(j,k,2,1,loop1)  &
                        +q_t(i5,j,k,1)-qwave_j_w(j,k,3,1,loop1)  &
                        +q_t(i5,j,k,2)-qwave_j_w(j,k,3,2,loop1)))
                  ema2_q_j(j,k,loop1)=    &
      (1.-ema_mu)*ema2_q_j(j,k,loop1) &
         +ema_mu *(-q_t(i3,j,k,2)+qwave_j_w(j,k,2,2,loop1) &
                   -q_t(i3,j,k,1)+qwave_j_w(j,k,2,1,loop1) &
                   +q_t(i5,j,k,1)-qwave_j_w(j,k,3,1,loop1) &
                   +q_t(i5,j,k,2)-qwave_j_w(j,k,3,2,loop1))**2
      enddo         ; enddo
#endif

! Condition dirichlet analytique basEe sur ssh calculee (utilise kvectorpeak_j(j,loop1))
       if(obcnh_scheme==1) then !111>
        do k=1,ktop-1 ; do j=2,jmax-1
              q_t(i1,j,k,2)= &
           mask_t(i1,j,kmax)*( & !ooo>             
                   qwave_j_w(j,k,1,2,loop1)-grav*(ssh_int_w(i1,j,2)-sshwave_j_w(j,1,2,loop1))   &
       +freq2pipeak/kvectorpeak_j(j,loop1)*(i1-i3)*(vel_u(i2,j,k,2)-velwave_j_u(j,k,1,2,loop1)) &
                             )   !ooo>
        enddo         ; enddo
       endif                    !111>

      if(obcnh_scheme==17) then !1717>
        do k=1,ktop-1 ; do j=2,jmax-1
! Utilisant la vitesse de phase 2D de la SSH:
           q_t(i1,j,k,2)=(                                                &
                          ema1_s_j(j,loop1)*(-q_t(i1,j,k,1)+q_t(i3,j,k,1)+q_t(i3,j,k,2)) & 
                         +ema2_s_j(j,loop1)*( q_t(i1,j,k,1)+q_t(i3,j,k,1)-q_t(i3,j,k,2)) & 
                       )/(ema1_s_j(j,loop1)+ema2_s_j(j,loop1))

! Utilisant la vitesse de phase 3D de Q
!          q_t(i1,j,k,2)=(                                                &
!                         ema1_q_j(j,k,loop1)*(-q_t(i1,j,k,1)+q_t(i3,j,k,1)+q_t(i3,j,k,2)) & 
!                        +ema2_q_j(j,k,loop1)*( q_t(i1,j,k,1)+q_t(i3,j,k,1)-q_t(i3,j,k,2)) & 
!                      )/(ema1_q_j(j,k,loop1)+ema2_q_j(j,k,loop1))
        enddo         ; enddo
      endif                     !1717>

       endif                                                ! m[0v0]m  
      enddo ! boucle loop1

! Frontiere en j=1 ou j=jmax
      do loop1=1,2 ! Boucles sur 2 frontieres (j=1 puis j=jmax)

      if( (obcstatus(jeq1)   ==1.and.loop1==1).or.          &
          (obcstatus(jeqjmax)==1.and.loop1==2)     )   then ! m[0v0]m

       if(loop1==1)obcnh_scheme=1
       if(loop1==2)obcnh_scheme=17

       if(loop1==1) then !obc j=1>
         j1=1    ; j3=j1+1 ; j5=j3+1 ; j7=j5+1  ! pour ssh et q
         j2=2    ; j4=j2+1 ; j6=j4+1 ; j8=j6+1  ! pour v
       endif           !obc j=1>
       if(loop1==2) then !obc j=jmax>
         j1=jmax ; j3=j1-1 ; j5=j3-1 ; j7=j5-1 ! pour ssh et q
         j2=jmax ; j4=j2-1 ; j6=j4-1 ; j8=j6-1 ! pour v
       endif           !obc j=jmax>

! Condition dirichlet analytique basEe sur ssh calculee (utilise kvectorpeak_i(i,loop1))
       if(obcnh_scheme==1) then !111>
        do k=1,ktop-1 ; do i=2,imax-1
              q_t(i,j1,k,2)=                                     &
           mask_t(i,j1,kmax)*( & !ooo>             
                      qwave_i_w(i,k,1,2,loop1)-grav*(ssh_int_w(i,j1,2)-sshwave_i_w(i,1,2,loop1)) &
!       +freq2pipeak/kvectorpeak_i(i,loop1)*(j1-j3)*(vel_v(i,j2,k,2)-velwave_i_v(i,k,1,2,loop1)) &
        +freq2pipeak/kvectorpeak_i(i,loop1)*(j1-j3)*(vel_v(i,j2,k,2)-velwave_i_v(i,k,1,2,loop1)  &
                                                    +vel_v(i,j2,k,1)-velwave_i_v(i,k,1,1,loop1)) &
                             )   !ooo>

!       if(k==kmax/2.and.i==imax/2.and.loop1==1) then
!        write(300+par%rank,*)real(q_t(i,j1,k,2)),real(qwave_i_w(i,k,1,2,loop1)) &
!       ,real(vel_v(i,j2,k,2)),real(velwave_i_v(i,k,1,2,loop1)) &
!       ,real(ssh_int_w(i,j1,2)),real(sshwave_i_w(i,1,2,loop1))
!       endif
!       if(i+par%timax(1)==2476.and.k==kmax/2) write(20+par%rank,*)    &
!        real(elapsedtime_now)                                         &
!       ,real(q_t(i,j1,k,2)),real(qwave_i_w(i,k,1,2,loop1))            &
!       ,real(vel_v(i,j2,k,2)),real(velwave_i_v(i,k,1,2,loop1)) 

        enddo         ; enddo
       endif                    !111>

      if(obcnh_scheme==17) then !1717>
        do k=1,ktop-1 ; do i=2,imax-1
!          q_t(i,j1,k,2)=(                                                &
!                         ema1_s_i(i,loop1)*(-q_t(i,j1,k,1)+q_t(i,j3,k,1)+q_t(i,j3,k,2)) & 
!                        +ema2_s_i(i,loop1)*( q_t(i,j1,k,1)+q_t(i,j3,k,1)-q_t(i,j3,k,2)) & 
!                      )/(ema1_s_i(i,loop1)+ema2_s_i(i,loop1))
               q_t(i,j1,k,2) &
        =qwave_i_w(i,k,1,2,loop1)+( &
          ema1_s_i(i,loop1)*(-q_t(i,j1,k,1)+qwave_i_w(i,k,1,1,loop1)  &
                             +q_t(i,j3,k,1)-qwave_i_w(i,k,2,1,loop1)  &
                             +q_t(i,j3,k,2)-qwave_i_w(i,k,2,2,loop1)) & 
         +ema2_s_i(i,loop1)*(+q_t(i,j1,k,1)-qwave_i_w(i,k,1,1,loop1)  &
                             +q_t(i,j3,k,1)-qwave_i_w(i,k,2,1,loop1)  &
                             -q_t(i,j3,k,2)+qwave_i_w(i,k,2,2,loop1)) & 
                                  )/(ema1_s_i(i,loop1)+ema2_s_i(i,loop1))
        enddo         ; enddo
      endif                     !1717>

       endif                                                ! m[0v0]m  
      enddo ! boucle loop1

      end subroutine q_obc_radiative

!..............................

      subroutine q_u3d_from_u2d
      implicit none

! Cas ou on ajuste la moyenne seulement au moment des graphiques
      if(mod(iteration2d,nh2d_graph_period)==0  &
             .or.iteration2d==iteration2d_max_now) then !ooo> 

       do j=1,jmax ; do i=1,imax
        anyvar2d(i,j)=0. ! small1
       enddo       ; enddo
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
        anyvar2d(i,j)= &
        anyvar2d(i,j)+dsigr4_t(i,j,k)*abs(q_t(i,j,k,1)+grav*ssh_w(i,j,1))
       enddo    ; enddo       ; enddo
! u:
      do j=2,jmax-1 ; do i=2,imax
        x1=(anyvar2d(i-1,j)+anyvar2d(i,j))*velbar_u(i,j,1) &
         /((anyvar2d(i-1,j)+anyvar2d(i,j))**2+small1)
        do k=1,kmax
         vel_u(i,j,k,2)=x1*( abs(q_t(i  ,j  ,k,1)+grav*ssh_w(i  ,j  ,1)) &
                            +abs(q_t(i-1,j  ,k,1)+grav*ssh_w(i-1,j  ,1)))
        enddo
       enddo    ; enddo
! v
      do j=2,jmax ; do i=2,imax-1
        x1=(anyvar2d(i,j-1)+anyvar2d(i,j))*velbar_v(i,j,1) &
         /((anyvar2d(i,j-1)+anyvar2d(i,j))**2+small1)
        do k=1,kmax
         vel_v(i,j,k,2)=x1*( abs(q_t(i  ,j  ,k,1)+grav*ssh_w(i  ,j  ,1)) &
                            +abs(q_t(i  ,j-1,k,1)+grav*ssh_w(i  ,j-1,1)))
        enddo
       enddo    ; enddo

       call obc_int_mpi(2,1)

! Verification:
!      i=imax/2 ; j=jmax/2
!      i=5 ; j=jmax/2
!      sum1=0.
!      do k=1,kmax
!       sum1=sum1+dsigr4_t(i,j,k)*vel_u(i,j,k,2)
!      enddo
!      write(10+par%rank,*)sum1,velbar_u(i,j,1)

      endif                                             !ooo>

#ifdef bidon
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

! Version classique
!      vel_u(i,j,k,2)=                                  &
!      vel_u(i,j,k,2)-dte_fw*( & !>>>

!                     grav*( ssh_w(i  ,j,1)             &
!                           -ssh_w(i-1,j,1) )/dx_u(i,j) &

!      +nhpgf_u(i,j,k)/(0.5*(anyv3dr4(i,j,k,id_dz)      &
!                           +anyv3dr4(i-1,j,k,id_dz)))  &

!                            )*wetmask_u(i,j)   !>>>

! Version "transport" reduite A QPGF (suivie de l'ajustement au courant moyen)
       veldydz_u(i,j,k,1)=                                  &
       veldydz_u(i,j,k,1)-dte_fw*nhpgf_u(i,j,k)

      enddo       ; enddo         ; enddo
      do k=1,kmax ; do j=2,jmax   ; do i=2,imax-1
! Version "transport" reduite A QPGF (suivie de l'ajustement au courant moyen)
       veldxdz_v(i,j,k,1)=                                  &
       veldxdz_v(i,j,k,1)-dte_fw*nhpgf_v(i,j,k)
      enddo       ; enddo         ; enddo

! Cas ou on ajuste la moyenne seulement au moment des graphiques
      if(mod(iteration2d,nh2d_graph_period)==0  &
             .or.iteration2d==iteration2d_max_now) then !ooo> 
    
       k=1
       do j=2,jmax-1 ; do i=2,imax
        veldydz_u(i,j,0,1)=veldydz_u(i,j,k,1)
       enddo         ; enddo
       do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax
        veldydz_u(i,j,0,1)= &
        veldydz_u(i,j,0,1)+veldydz_u(i,j,k,1) ! ne pas utiliser xy_u peut etre occupE
       enddo       ; enddo         ; enddo
! Correction de la moyenne . Attention la division par kmax suppose couches homogenes
       do j=2,jmax-1 ; do i=2,imax
        veldydz_u(i,j,0,1)=(fluxbar_u(i,j,1)/dy_u(i,j)-veldydz_u(i,j,0,1))/real(kmax)
       enddo         ; enddo
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
         vel_u(i,j,k,2)= & ! Tableau "Courant dans graph_out"
        (veldydz_u(i,j,k,1)+veldydz_u(i,j,0,1)) &
       /(0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))) ! conversion transport courant
    
!       if(par%rank==0.and.i==imax/2.and.j==jmax/2.and.k==kmax/2) then
!        write(666,*)vel_u(i,j,k,2)
!       endif

       enddo       ; enddo         ; enddo

       k=1
       do j=2,jmax   ; do i=2,imax-1
        veldxdz_v(i,j,0,1)=veldxdz_v(i,j,k,1)
       enddo         ; enddo
       do k=2,kmax ; do j=2,jmax   ; do i=2,imax-1
        veldxdz_v(i,j,0,1)= &
        veldxdz_v(i,j,0,1)+veldxdz_v(i,j,k,1) 
       enddo       ; enddo         ; enddo
! Correction de la moyenne . Attention la division par kmax suppose couches homogenes
       do j=2,jmax   ; do i=2,imax-1
        veldxdz_v(i,j,0,1)=(fluxbar_v(i,j,1)/dx_v(i,j)-veldxdz_v(i,j,0,1))/real(kmax)
       enddo         ; enddo
       do k=1,kmax ; do j=2,jmax   ; do i=2,imax-1
         vel_v(i,j,k,2)= & ! Tableau "Courant dans graph_out"
        (veldxdz_v(i,j,k,1)+veldxdz_v(i,j,0,1)) &
       /(0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz)))
       enddo       ; enddo         ; enddo

       call obc_int_mpi(2,1)
       
      endif                                             !ooo>
#endif

#ifdef bidon
! Ne pas utiliser cette routine tant que modifier vel_u a un
! feedback sur le  velbot3d2d.....     
! D'autre part quand le gradient de pression passe par zero cette
! param ne marche plus...
      stop 'ne pas utiliser q_u3d_from_u2d'

! Le courant 3D a le cisaillement du gradient de pression et 
! sa moyenne egale le courant moyen:

      do j=2,jmax-1 ; do i=2,imax
       xy_u(i,j,0:1)=small1
      enddo         ; enddo

      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
! First guess:
       vel_u(i,j,k,1)=grav*( ssh_w(i  ,j,1)             &
                            -ssh_w(i-1,j,1) )/dx_u(i,j) &
       +nhpgf_u(i,j,k)/(0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz)))

       xy_u(i,j,0)=xy_u(i,j,0)+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))
       xy_u(i,j,1)=xy_u(i,j,1)+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))*vel_u(i,j,k,1)
      enddo ; enddo ; enddo

      do j=2,jmax-1 ; do i=2,imax
       xy_u(i,j,1)=velbar_u(i,j,1)*xy_u(i,j,0)/xy_u(i,j,1)
      enddo         ; enddo

      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
! Normalise:
       vel_u(i,j,k,2)=xy_u(i,j,1)*vel_u(i,j,k,1)
      enddo ; enddo ; enddo

!..................

      do j=2,jmax ; do i=2,imax-1
       xy_v(i,j,0:1)=small1
      enddo ; enddo
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
! First guess:
       vel_v(i,j,k,1)=grav*( ssh_w(i,j  ,1)             &
                            -ssh_w(i,j-1,1) )/dy_v(i,j) &
       +nhpgf_v(i,j,k)/(0.5*(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz)))

       xy_v(i,j,0)=xy_v(i,j,0)+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))
       xy_v(i,j,1)=xy_v(i,j,1)+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))*vel_v(i,j,k,1)

      enddo ; enddo ; enddo

      do j=2,jmax ; do i=2,imax-1
       xy_v(i,j,1)=velbar_v(i,j,1)*xy_v(i,j,0)/xy_v(i,j,1)
      enddo ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
! Normalise:
       vel_v(i,j,k,2)=xy_v(i,j,1)*vel_v(i,j,k,1)
      enddo       ; enddo       ; enddo

! VERIFICATION
      sum1=0. 
      sum2=0. 
      sum3=0. 
      sum4=0. 
      i=imax/2 ; j=jmax/2
      do k=1,kmax
       sum1=sum1+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))
       sum2=sum2+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))
       sum3=sum3+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i-1,j,k,id_dz))*vel_u(i,j,k,2)
       sum4=sum4+(anyv3dr4(i,j,k,id_dz)+anyv3dr4(i,j-1,k,id_dz))*vel_v(i,j,k,2)
      enddo 
      write(10+par%rank,*)'---------'
      write(10+par%rank,*)sum3/sum1,velbar_u(i,j,1)
      write(10+par%rank,*)sum4/sum2,velbar_v(i,j,1)
#endif


      end subroutine q_u3d_from_u2d

!..............................

      subroutine q_moteur_remapping_b2
      implicit none

! Je commente ces lignes car elles sont dans add_nhpgf....
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3dr4(i,j,k,id_zt)=sig1d_t(k)*(h_w(i,j)+ssh_w(i,j,1))+ssh_w(i,j,1)
!     enddo       ; enddo  ; enddo
! Surface et fond
!     do j=1,jmax ; do i=1,imax
!      anyv3dr4(i,j,kmax+1,id_zt)=0.
!      anyv3dr4(i,j,0)=-h_w(i,j)
!     enddo  ; enddo

! rappel: z_t(k+1)-z_t(k)=H*dsig_w(k+1)
! rappel: z_w(k+1)-z_w(k)=H*dsig_t(k)
      
! sshr4_w a EtE calculEe dans subroutine q_nhpgf2d
      ww_=w2_ ; if(ww_==0.)ww_=1.
      do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
       mc(i,j,k,1)=-ww_*dt2/( dsigr4_t(i,j,k)*dsig_w(i,j,k)*(h_w(i,j)+sshr4_w(i,j))**2 )
      enddo       ; enddo       ;enddo

      do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
       mc(i,j,k,3)=mc(i,j,k,1)*dsig_w(i,j,k)/dsig_w(i,j,k+1)
      enddo       ; enddo       ;enddo

      do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
       mc(i,j,k,2)=1.-mc(i,j,k,1)-mc(i,j,k,3)
      enddo       ; enddo       ;enddo

!.......................................................
! C.L. Surface: q_t(ktop)=0.
      k=ktop
      do j=2,jmax-1 ; do i=2,imax-1
       mc(i,j,k,1)=0.
       mc(i,j,k,2)=1.
       mc(i,j,k,3)=0.
       mc(i,j,k,4)=0.
      enddo       ; enddo

!.......................................................
! C.L. Fond: q_t(0)=q_t(1)
      k=1
      do j=2,jmax-1 ; do i=2,imax-1
       mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
       mc(i,j,k,1)=0.
      enddo       ; enddo

      ww_=1.
      if(w2_/=0.)ww_=1./w2_

      do j=2,jmax-1 ; do i=2,imax-1

       anyvar2d(i,j)=dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>
       ((fluxbar_u(i+1,j,0)-fluxbar_u(i+1,j,1))/(dte_fw*dy_u(i+1,j))-nhpgf2d_u(i+1,j))/(h_u(i+1,j)+0.5*(sshr4_w(i+1,j)+sshr4_w(i,j))) &
      -((fluxbar_u(i  ,j,0)-fluxbar_u(i  ,j,1))/(dte_fw*dy_u(i  ,j))-nhpgf2d_u(i  ,j))/(h_u(i  ,j)+0.5*(sshr4_w(i-1,j)+sshr4_w(i,j))) &
                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>
       ((fluxbar_v(i,j+1,0)-fluxbar_v(i,j+1,1))/(dte_fw*dx_v(i,j+1))-nhpgf2d_v(i,j+1))/(h_v(i,j+1)+0.5*(sshr4_w(i,j+1)+sshr4_w(i,j))) &
      -((fluxbar_v(i,j  ,0)-fluxbar_v(i,j  ,1))/(dte_fw*dx_v(i,j  ))-nhpgf2d_v(i,j  ))/(h_v(i,j  )+0.5*(sshr4_w(i,j-1)+sshr4_w(i,j))) &
                    ) & !ooo>

                         )    !ttt>

      enddo ; enddo

      do k=1,ktop-1 
      do j=2,jmax-1 ; do i=2,imax-1

        mc(i,j,k,4)=2*q_t(i,j,k,1)-q_t(i,j,k,0)                        &

             -(q_t(i,j,k+1,1)-     q_t(i,j,k-1,1))                     &!b3-time
        /(anyv3dr4(i,j,k+1  ,id_zt)-anyv3dr4(i,j,k-1  ,id_zt))                     &!b3-time
        *(1.+sig1dpom_t(i,j,k))*(2*ssh_w(i,j,1)-ssh_w(i,j,0)-ssh_w(i,j,2)) &!b3-time

                +anyvar2d(i,j)   &
                                              +dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>

                ( q_t(i+1,j,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!               +(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b1
!          /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b1
!      *0.5*(anyv3dr4(i-1,j,k    ,id_zt)-anyv3dr4(i+1,j,k    ,id_zt)) & !b1

                +(q_t(i+1,j,k+1,1)-     q_t(i+1,j,k-1,1)) & !b2
           /(anyv3dr4(i+1,j,k+1  ,id_zt)-anyv3dr4(i+1,j,k-1  ,id_zt)) & !b2
           *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i+1,j,k    ,id_zt)) & !b2

!               +(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b3
!          /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b3
!          *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i+1,j,k    ,id_zt)) & !b3

                                -q_t(i  ,j,k,1) )*invdx_u(i+1,j)    &

               -( q_t(i  ,j,k,1)-q_t(i-1,j,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!                              -(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b1
!                         /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b1
!                     *0.5*(anyv3dr4(i+1,j,k    ,id_zt)-anyv3dr4(i-1,j,k    ,id_zt)) & !b1

                               -(q_t(i-1,j,k+1,1)-     q_t(i-1,j,k-1,1)) & !b2
                          /(anyv3dr4(i-1,j,k+1  ,id_zt)-anyv3dr4(i-1,j,k-1  ,id_zt)) & !b2
                          *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i-1,j,k    ,id_zt)) & !b2

!                              -(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b3
!                         /(anyv3dr4(i  ,j,k+1  ,id_zt)-anyv3dr4(i  ,j,k-1  ,id_zt)) & !b3
!                         *(anyv3dr4(i  ,j,k    ,id_zt)-anyv3dr4(i-1,j,k    ,id_zt)) & !b3

                                                )*invdx_u(i  ,j)    &

                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>

                ( q_t(i,j+1,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!               +(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b1
!          /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b1
!      *0.5*(anyv3dr4(i,j-1,k    ,id_zt)-anyv3dr4(i,j+1,k    ,id_zt)) & !b1

                +(q_t(i,j+1,k+1,1)-     q_t(i,j+1,k-1,1)) & !b2
           /(anyv3dr4(i,j+1,k+1  ,id_zt)-anyv3dr4(i,j+1,k-1  ,id_zt)) & !b2
           *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j+1,k    ,id_zt)) & !b2

!               +(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b3
!          /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b3
!          *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j+1,k    ,id_zt)) & !b3

                                -q_t(i  ,j,k,1) )*invdy_v(i,j+1)    &

               -( q_t(i  ,j,k,1)-q_t(i,j-1,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!                              -(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b1
!                         /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b1
!                     *0.5*(anyv3dr4(i,j+1,k    ,id_zt)-anyv3dr4(i,j-1,k    ,id_zt)) & !b1

                               -(q_t(i,j-1,k+1,1)-     q_t(i,j-1,k-1,1)) & !b2
                          /(anyv3dr4(i,j-1,k+1  ,id_zt)-anyv3dr4(i,j-1,k-1  ,id_zt)) & !b2
                          *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j-1,k    ,id_zt)) & !b2

!                              -(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b3
!                         /(anyv3dr4(i,j  ,k+1  ,id_zt)-anyv3dr4(i,j  ,k-1  ,id_zt)) & !b3
!                         *(anyv3dr4(i,j  ,k    ,id_zt)-anyv3dr4(i,j-1,k    ,id_zt)) & !b3

                                                )*invdy_v(i  ,j)    &

                    ) & !ooo>

                                                   ) &  !ttt>


       -mc(i,j,k,1)    *ww_*(w0_*q_t(i,j,k-1,0)+w1_*q_t(i,j,k-1,1) )&
       -mc(i,j,k,3)    *ww_*(w0_*q_t(i,j,k+1,0)+w1_*q_t(i,j,k+1,1) )&
      -(mc(i,j,k,2)-1.)*ww_*(w0_*q_t(i,j,k  ,0)+w1_*q_t(i,j,k  ,1) ) 

! Note sur sponge_t: varie de 0 a 1 (sur la frontiere) ce qui equivaut en termes
! d'echelle de temps associee a un temps de rappel equivalent au pas de temps. C'est
! donc un rappel fort qui a pour but de supprimer la pression NH aux frontieres ouvertes

      enddo ; enddo
      enddo

! Solveur vertical tridiagonal:
      call q_solver(1,2,imax-1,2,jmax-1,1,ktop) ! donne q_t(t+1)

!........................
! ZONES EPONGE ET WETMASK: (Ne pourrait elle pas etre implicitee ou simplifiee?)
      do j=2,jmax-1 ; do i=2,imax-1
       wetmask_t(i,j)=min(1.,(h_w(i,j)+sshr4_w(i,j))/wetdry_cstnh)
      enddo         ; enddo
!     do k=0,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1
      do k=1,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! Boucle k commence A 1 si report des C.L. en suivant....
!      q_t(i,j,k,2)=q_t(i,j,k,2)*(1.-sponge_t(i,j,1))*wetmask_t(i,j)
       q_t(i,j,k,2)=q_t(i,j,k,2)                     *wetmask_t(i,j)
      enddo         ; enddo       ; enddo


!........................
! Conditions ouvertes radiatives
      call q_obc_radiative

!........................
! Reporter les hypotheses de conditions aux limites au fond et en surface:
      do j=2,jmax-1 ; do i=2,imax-1
! C.L.: q(0)=q(1):
       q_t(i,j,0,2)=q_t(i,j,1,2)*mask_t(i,j,kmax)
! C.L.: q(kmax+1)=0.
       q_t(i,j,ktop:kmax+1,2)=0.
      enddo ; enddo

!........................
! Continuite MPI:
      call q_obc_mpi('z0') ! mpi 'z1' sur q_t(t+1)

!#ifdef bidon
! VERIFICATION DE L'EQUILIBRE:
      if(modulo(iteration2d,100)==0) then !>>>
! Ici verification possible de l'equilibre
! Attention que le filtre d'asselin se fait ulterieurement
! et par consequent n'est pas evalue ici
      i=imax/2 ; j=jmax/2 ; k=kmax/2
      write(texte30,'(a,i0)')'tmp/nhdiag_balance',par%rank
! opposE du Terme de time stepping:
       x1= -q_t(i,j,k,2)+2*q_t(i,j,k,1)-q_t(i,j,k,0)                  &
             -(q_t(i,j,k+1,1)-     q_t(i,j,k-1,1))                     &!b3-time
        /(anyv3dr4(i,j,k+1  ,id_zt)-anyv3dr4(i,j,k-1  ,id_zt))                     &!b3-time
        *(1.+sig1dpom_t(i,j,k))*(2*ssh_w(i,j,1)-ssh_w(i,j,0)-ssh_w(i,j,2))  !b3-time
! Forcage hydrostatique:
       x2=anyvar2d(i,j)
! Laplacien 3D de la pression NH deduit de l'equilibre elementaire
! m1*q(k-1,t+1)+m2*q(k  ,t+1)+m3*q(k+1,t+1)=m4 que l'on a pris
! soin prealablement de verifier
       x3=-q_t(i,j,k,2)+mc(i,j,k,4)-x1-x2       &
                 -mc(i,j,k,1)   *q_t(i,j,k-1,2) &
                -(mc(i,j,k,2)-1)*q_t(i,j,k  ,2) &
                 -mc(i,j,k,3)   *q_t(i,j,k+1,2)  

      open(unit=3,file=texte30,position='append')
! Equilibre elementaire qui doit etre verifiE:
!      write(3,*)mc(i,j,k,1)*q_t(i,j,k-1,2) &
!               +mc(i,j,k,2)*q_t(i,j,k  ,2) &
!               +mc(i,j,k,3)*q_t(i,j,k+1,2) &
!               ,mc(i,j,k,4)
! l'addition de 2+3+4 egale zero:
       write(3,*)real(iteration2d*dte_fw/3600.) & ! 1 temps en heure
                ,real(x1)                       & ! 2 terme de time-stepping (ou desequilibre)
                ,real(x2)                       & ! 3 forcage hydrostatique
                ,real(x3)                         ! 4 Laplacien 3D non-hydrostatique
! Si le solveur a bien convergE, 2 doit etre petit par rapport a 3 et 4
      close(3)
      endif                               !>>>
!#endif

      end subroutine q_moteur_remapping_b2

!..............................

      subroutine q_wavebreaking
      implicit none
      integer loop_,id_u1_,id_u2_,id_v1_,id_v2_

!     RETURN

      id_flx=3
      id_fly=4
      id_breaker=5

! DEFERLEMENT PAR DIFFUSION HORIZONTALE MAXIMALE
      if(iteration2d/=iteration2d_begin) then
        if(checkxyt(2)/=xy_t(imax/2,jmax/2,id_flx)) &
        stop 'err checkxyt(2) q_wavebreaking'
        if(checkxyt(3)/=xy_t(imax/2,jmax/2,id_fly)) &
        stop 'err checkxyt(3) q_wavebreaking'
        if(checkxyt(4)/=xy_t(imax/2,jmax/2,id_breaker)) &
        stop 'err checkxyt(4) q_wavebreaking'
      endif

!     do j=0,jmax+1 ; do i=0,imax+1
!      ssh_w(i,j,1)=j*dxb*0.01
!     enddo         ; enddo
      
      x1=0.5*(1./0.14)**2
      do j=1,jmax ; do i=1,imax
! xy_t(i,j,id_breaker) varie de 0 (ssh plate) a 1 (pente d ssh = 0.14)
!      xy_t(i,j,id_breaker)=min(1.,                                      &
!        x1*( ( (ssh_w(i+1,j  ,1)-ssh_w(i  ,j  ,1))*invdx_u(i+1,j  ) )**2 &
!            +( (ssh_w(i  ,j  ,1)-ssh_w(i-1,j  ,1))*invdx_u(i  ,j  ) )**2 &
!            +( (ssh_w(i  ,j  ,1)-ssh_w(i  ,j-1,1))*invdy_v(i  ,j  ) )**2 &
!            +( (ssh_w(i  ,j+1,1)-ssh_w(i  ,j  ,1))*invdy_v(i  ,j+1) )**2 ))
!      xy_t(i,j,id_breaker)=0.5*abs(velbar_u(i+1,j,2)-velbar_u(i,j,2))

! u**2/c**2
       xy_t(i,j,id_breaker)=min(1.,             &
              dte_fw*wavebreakfactor            & !05-12-17 
                *0.5*( velbar_u(i+1,j  ,2)**2   &
                      +velbar_u(i  ,j  ,2)**2   &
                      +velbar_v(i  ,j  ,2)**2   &
                      +velbar_v(i  ,j+1,2)**2 ) &
             /(    grav*max(h_w(i  ,j  )    ,0.0001)))

! u/c
!      xy_t(i,j,id_breaker)=sqrt( min(1., &
!!               0.2* & ! Ca c'est une facteur d'amplification
!                0.5*( velbar_u(i+1,j  ,2)**2     &
!                     +velbar_u(i  ,j  ,2)**2     &
!                     +velbar_v(i  ,j  ,2)**2     &
!                     +velbar_v(i  ,j+1,2)**2 )   &
!               /(grav*max( h_w(i  ,j  ) ,0.0001))) )

!      write(666,*)i,j,xy_t(i,j,id_breaker)
!      if(iteration2d==1000.and.j==jmax/2)write(10+par%rank,*)i+par%timax(1),xy_t(i,j,id_breaker),h_w(i,j)
      enddo       ; enddo
      checkxyt(4)=xy_t(imax/2,jmax/2,id_breaker)
      
!     if(iteration2d==1001)stop 'cocote'

      do j=2,jmax-1 ; do i=1,imax

        xy_t(i,j,id_flx)=                                        &
!                        0.25 & ! Valeur theorique maximale: 0.25
                         0.2  & ! <0.25 permet de resorber les ultimes oscillations 2dx sur ssh
                            *xy_t(i,j,id_breaker)               &
              *min(h_u(i  ,j)+0.5*(ssh_w(i-1,j,1)+ssh_w(i,j,1))  &
                  ,h_u(i+1,j)+0.5*(ssh_w(i+1,j,1)+ssh_w(i,j,1))) &
            *(velbar_u(i+1,j,2)-velbar_u(i,j,2))
      enddo        ; enddo
      checkxyt(2)=xy_t(imax/2,jmax/2,id_flx)

      do j=2,jmax-1 ; do i=2,imax
           velbar_u(i,j,2)=mask_u(i,j,kmax)*( & !pmx>
           velbar_u(i,j,2)+( xy_t(i  ,j,id_flx)                   &
                            -xy_t(i-1,j,id_flx) )                 &
           /max(h_u(i,j)+0.5*(ssh_w(i-1,j,1)+ssh_w(i,j,1)),0.001) &
                                              ) !pmx>
      enddo         ; enddo

      do j=1,jmax ; do i=2,imax-1
        xy_t(i,j,id_fly)=                                        &
!                        0.25 & ! Valeur theorique maximale: 0.25
                         0.2  & ! Valeur theorique maximale: 0.25
                            *xy_t(i,j,id_breaker)               &
              *min(h_v(i,j  )+0.5*(ssh_w(i,j-1,1)+ssh_w(i,j,1))  &
                  ,h_v(i,j+1)+0.5*(ssh_w(i,j+1,1)+ssh_w(i,j,1))) &
            *(velbar_v(i,j+1,2)-velbar_v(i,j,2))
      enddo       ; enddo
      checkxyt(3)=xy_t(imax/2,jmax/2,id_fly)
      do j=2,jmax ; do i=2,imax-1
           velbar_v(i,j,2)=mask_v(i,j,kmax)*( & !pmx>
           velbar_v(i,j,2)+( xy_t(i,j  ,id_fly)                   &
                            -xy_t(i,j-1,id_fly) )                 &
           /max(h_v(i,j)+0.5*(ssh_w(i,j-1,1)+ssh_w(i,j,1)),0.001) &
                                              ) !pmx>
      enddo       ; enddo

!...........
      call get_type_echange('u1','ubar_u1_2',velbar_u  &
                                     ,lbound(velbar_u) &
                                     ,ubound(velbar_u) &
                                     ,2                &
                                     ,id_u1_)
      call get_type_echange('u2','ubar_u2_2',velbar_u  &
                                     ,lbound(velbar_u) &
                                     ,ubound(velbar_u) &
                                     ,2                &
                                     ,id_u2_)


      call get_type_echange('v1','vbar_u1_2',velbar_v  &
                                     ,lbound(velbar_v) &
                                     ,ubound(velbar_v) &
                                     ,2                &
                                     ,id_v1_)

      call get_type_echange('v2','vbar_u2_2',velbar_v  &
                                     ,lbound(velbar_v) &
                                     ,ubound(velbar_v) &
                                     ,2                &
                                     ,id_v2_)

      do loop_=1,subcycle_exchange
         call echange_voisin(velbar_u,id_u1_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_u,id_u2_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,id_v1_,mpi_neighbor_list(loop_))
         call echange_voisin(velbar_v,id_v2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
!...........

      RETURN

!      do j=1,jmax ; do i=1,imax
!           (0.5*(ssh_w(i+1,j,2)-ssh_w(i-1,j,2))/dx_t(i,j))**2 &
!          +(0.5*(ssh_w(i,j+1,2)-ssh_w(i,j-1,2))/dy_t(i,j))**2 &
!          -0.01
!      enddo       ; enddo

      end subroutine q_wavebreaking

!..............................
#ifdef bidon
      subroutine q_wavebreaking3d_old
      implicit none
      integer :: loop_,id_u1_,id_u2_,id_v1_,id_v2_ &
                ,looptime_,id_tken_                &
                ,looptimemax_=1
      real dti_wb_

! Sous pas de temps pour le wave breaking
      dti_wb_=real(dti_lp)/looptimemax_



!     RETURN

      id_uflx=3    ! xy_t
      id_ufly=3    ! xy_f
      id_vflx=4    ! xy_f
      id_vfly=4    ! xy_t
      id_breaker=5 ! xy_t xyf

! DEFERLEMENT PAR DIFFUSION HORIZONTALE MAXIMALE
      if(iteration2d/=iteration2d_begin) then

        flag_stop=0
        if(checkxyt(2)/=xy_t(imax/2,jmax/2,id_uflx)) &
          flag_stop=1

        if(checkxyt(5)/=xy_t(imax/2,jmax/2,id_vfly)) &
          flag_stop=1

        if(checkxyt(4)/=xy_t(imax/2,jmax/2,id_breaker)) &
          flag_stop=1

        if(checkxyf(1)/=xy_f(imax/2,jmax/2,id_ufly)) &
          flag_stop=1

        if(checkxyf(2)/=xy_f(imax/2,jmax/2,id_vflx)) &
          flag_stop=1

        if(checkxyf(3)/=xy_f(imax/2,jmax/2,id_breaker)) &
          flag_stop=1

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'Error 3705 module_q'

      endif

      x1_r4=        1./(1.-brk_crit_r)
      x2_r4=brk_crit_r/(1.-brk_crit_r)

      do j=1,jmax ; do i=1,imax

      breaker_t(i,j)=  & !02-05-20

      (1.-dti_fw/periodpeak)*max( & !MAX> 
       
              breaker_t(i,j)  &

        ,   x1_r4*(ssh_int_w(i,j,1)/max(h_w(i,j),0.03)*inv_brkh)**2-x2_r4         &

        ,   x1_r4*(                                            &  

              ( 0.5*( abs(ssh_int_w(i+1,j,1)-ssh_int_w(i,j,1)) &
                     +abs(ssh_int_w(i-1,j,1)-ssh_int_w(i,j,1)))*invdx_t(i,j) )**2  &

             +( 0.5*( abs(ssh_int_w(i,j+1,1)-ssh_int_w(i,j,1)) &
                     +abs(ssh_int_w(i,j-1,1)-ssh_int_w(i,j,1)))*invdy_t(i,j) )**2  &

                  )*inv_brkslope2-x2_r4                        & 

        ,   x1_r4*( ( &
                         abs( & !ooo>
                        !     ssh_int_w(i,j,2)-2*ssh_int_w(i,j,1)+ssh_int_w(i,j,0)   &
                             +ssh_int_w(i,j,1)-2*ssh_int_w(i,j,0)+ssh_int_w(i,j,-1)  &
                            ) & !ooo>
                    )*inv_dti_fw_p2     &
                  )*inv_brkvar-x2_r4    & 

                                )   !MAX>


          xy_t(i,j,id_breaker)=min( & !MIN> !02-05-20
                                   dx_t(i,j)**2                       & !1 
                                  ,dy_t(i,j)**2                       & !2
                                  ,dti_wb_*wavebreakfactor            & !3
                                          *wavebreakcoef              & !3
                                          *sqrt(max(h_w(i,j),0.0001)) & !3
                                          *breaker_t(i,j)             & !3
                                  )   !MIN>

#ifdef bidon
       if(i==imax/2.and.j==jmax/2) &
       write(10+par%rank,*)real(elapsedtime_now) &
        ,real(breaker_t(i,j)) &

        ,   real(x1_r4*(                                            &  

              ( 0.5*( abs(ssh_int_w(i+1,j,1)-ssh_int_w(i,j,1)) &
                     +abs(ssh_int_w(i-1,j,1)-ssh_int_w(i,j,1)))*invdx_t(i,j) )**2  &

             +( 0.5*( abs(ssh_int_w(i,j+1,1)-ssh_int_w(i,j,1)) &
                     +abs(ssh_int_w(i,j-1,1)-ssh_int_w(i,j,1)))*invdy_t(i,j) )**2  &

                  )*inv_brkslope2-x2_r4  )                      & 

        , real(  x1_r4*( ( &
                         abs( & !ooo>
                             +ssh_int_w(i,j,1)-2*ssh_int_w(i,j,0)+ssh_int_w(i,j,-1)  &
                            ) & !ooo>
                    )*inv_dti_fw_p2     &
                  )*inv_brkvar-x2_r4 )  & 
        ,' trouvemoi2'


#endif
     

      enddo       ; enddo
      checkxyt(4)=xy_t(imax/2,jmax/2,id_breaker)

      do j=2,jmax ; do i=2,imax

! pas de diffusion de u en y ni de v en x:
       xy_f(i,j,id_breaker)=0.25*( xy_t(i  ,j  ,id_breaker) &
                                  +xy_t(i-1,j  ,id_breaker) &
                                  +xy_t(i-1,j-1,id_breaker) &
                                  +xy_t(i  ,j-1,id_breaker) )

!      xy_f(i,j,id_breaker)=0.
!      if(h_f(i,j)<0.3) then
!         xy_f(i,j,id_breaker)=min(1.,dti_fw*wavebreakfactor*sqrt(max(h_f(i,j),0.03)/0.25))
!      endif

      enddo       ; enddo
      checkxyf(3)=xy_f(imax/2,jmax/2,id_breaker)



      do looptime_=1,looptimemax_    !PMXPMX>

      do k=1,kmax

! Flux X pour U:
      do j=2,jmax-1 ; do i=1,imax

        xy_t(i,j,id_uflx)=                         &
        xy_t(i,j,id_breaker)*                      &
                         0.1  & ! diffusion dans les 2 directions
                *mask_t(i,j,kmax)                  & ! slip condition on breaking
              *min(dz_u(i+1,j,k,1) ,dz_u(i,j,k,1)) &
                *(vel_u(i+1,j,k,2)-vel_u(i,j,k,2))*dy_t(i,j)*invdx_t(i,j)
      enddo        ; enddo
      checkxyt(2)=xy_t(imax/2,jmax/2,id_uflx)

! Flux Y pour U:
      do j=2,jmax   ; do i=2,imax
        xy_f(i,j,id_ufly)=                         &
        xy_f(i,j,id_breaker)*                      &
                         0.05  & ! diffusion dans les 2 directions + 1/2 somme des 2 gradients
                *mask_f(i,j,kmax)                    & ! slip condition on breaking
              *min(dz_u(i,j,k,1) ,dz_u(i,j-1,k,1))*( & !ooo>
                 (vel_u(i,j,k,2)-vel_u(i  ,j-1,k,2))*dx_f(i,j)*invdy_f(i,j) &
                 +vel_v(i,j,k,2)-vel_v(i-1,j  ,k,2)  &
                                                   )   !ooo>
      enddo        ; enddo
      checkxyf(1)=xy_f(imax/2,jmax/2,id_ufly)

! Flux Y pour V:
      do j=1,jmax ; do i=2,imax-1
        xy_t(i,j,id_vfly)=                         &
        xy_t(i,j,id_breaker)*                      &
                         0.1  & ! diffusion dans les 2 directions
                *mask_t(i,j,kmax)                    & ! slip condition on breaking
              *min(dz_v(i,j+1,k,1) ,dz_v(i,j,k,1)) &
                *(vel_v(i,j+1,k,2)-vel_v(i,j,k,2))*dx_t(i,j)*invdy_t(i,j)
      enddo       ; enddo
      checkxyt(5)=xy_t(imax/2,jmax/2,id_vfly)

! Flux X pour V:
      do j=2,jmax ; do i=2,imax  
        xy_f(i,j,id_vflx)=                         &
        xy_f(i,j,id_breaker)*                      &
                         0.05  & ! diffusion dans les 2 directions + 1/2 somme des 2 gradients
                *mask_f(i,j,kmax)                    & ! slip condition on breaking
              *min(dz_v(i,j,k,1) ,dz_v(i-1,j,k,1))*( & !ooo>
                 (vel_v(i,j,k,2)-vel_v(i-1,j  ,k,2))*dy_f(i,j)*invdx_f(i,j) &
                 +vel_u(i,j,k,2)-vel_u(i  ,j-1,k,2)  &
                                                   )   !ooo>
      enddo       ; enddo
      checkxyf(2)=xy_f(imax/2,jmax/2,id_vflx)

! Composante U
      do j=2,jmax-1 ; do i=2,imax
           vel_u(i,j,k,2)=mask_u(i,j,kmax)*( & !pmx>
           vel_u(i,j,k,2)+(                           &
                            xy_t(i  ,j  ,id_uflx)     &
                           -xy_t(i-1,j  ,id_uflx)     &
                           +xy_f(i  ,j+1,id_ufly)     &
                           -xy_f(i  ,j  ,id_ufly)     &

                          )/( & !ovo>
                               dxdy_u(i,j)*           &
                             max(dz_u(i,j,k,1),0.001) &
                            ) & !ovo>
                                             ) !pmx>

      enddo         ; enddo

! Composante V
      do j=2,jmax ; do i=2,imax-1
           vel_v(i,j,k,2)=mask_v(i,j,kmax)*( & !pmx>
           vel_v(i,j,k,2)+(                           &
                            xy_t(i  ,j  ,id_vfly)      &
                           -xy_t(i  ,j-1,id_vfly)      &
                           +xy_f(i+1,j  ,id_vflx)      &
                           -xy_f(i  ,j  ,id_vflx)      &

                           )/( & !ovo>
                               dxdy_v(i,j)*           &
                             max(dz_v(i,j,k,1),0.001) &
                             ) & !ovo>
                                              ) !pmx>
      enddo       ; enddo

      enddo ! k loop

#ifdef bidon
! ATTENTION A ACTIVER AUSSI LES ECHANGES MPI PLUS LOIN
! Apport de tke par breaking:
      if(iturbulence==2.or.iturbulence==3) &
      stop 'breaking tke input not possible with this tke scheme'
      x1_r4=0.1        & ! par ce que 0.1 est coef des flux
           *0.5*0.5    & ! pour la demi somme des vitesses sur level tke elevee au carre
           *0.5        & ! parce que xy_t(i,j,id_breaker) contient dti_lp alors que tke utilise dti_fw
!          *0.004        ! facteur reduction
           *0.01         ! facteur reduction du run 50 occigen
!          *0.001        ! facteur reduction du run 51 occigen
      do k=2,kmax
       do j=2,jmax-1
        do i=2,imax-1

         tken_w(i,j,k)=&
!        tken_w(i,j,k)+0.5*dti_fw*( & !pmx>
         tken_w(i,j,k)+x1_r4     *( & !pmx>

       xy_t(i,j,id_breaker)*( & !ooo>

       (( vel_u(i+1,j,k,2)+vel_u(i+1,j,k+1,2)                   &
         -vel_u(i  ,j,k,2)-vel_u(i  ,j,k+1,2))*invdx_t(i,j))**2 &

      +(( vel_v(i,j+1,k,2)+vel_v(i,j+1,k+1,2)                   &
         -vel_v(i,j  ,k,2)-vel_v(i,j  ,k+1,2))*invdy_t(i,j))**2 &

                            ) & !ooo>
      
                                  )   !pmx>
        enddo
       enddo
      enddo

!...........

      call get_type_echange('z0','tken_z0',tken_w   &
                                    ,lbound(tken_w) &
                                    ,ubound(tken_w) &
                                    ,id_tken_)
#endif

      call get_type_echange('u1','vel_u1_2',vel_u  &
                                    ,lbound(vel_u) &
                                    ,ubound(vel_u) &
                                    ,2             &
                                    ,id_u1_)

      call get_type_echange('u2','vel_u2_2',vel_u  &
                                    ,lbound(vel_u) &
                                    ,ubound(vel_u) &
                                    ,2             &
                                    ,id_u2_)

      call get_type_echange('v1','vel_v1_2',vel_v  &
                                    ,lbound(vel_v) &
                                    ,ubound(vel_v) &
                                    ,2             &
                                    ,id_v1_)

      call get_type_echange('v2','vel_v2_2',vel_v  &
                                    ,lbound(vel_v) &
                                    ,ubound(vel_v) &
                                    ,2             &
                                    ,id_v2_)

      do loop_=1,subcycle_exchange
#ifdef bidon
! ATTENTION A ACTIVER CET ECHANGE SI ALGO DE TRANSFERT DE BREAKING VERS TURBULENCE
         call echange_voisin(tken_w,id_tken_,mpi_neighbor_list(loop_))
#endif
         call echange_voisin(vel_u ,id_u1_  ,mpi_neighbor_list(loop_))
         call echange_voisin(vel_u ,id_u2_  ,mpi_neighbor_list(loop_))
         call echange_voisin(vel_v ,id_v1_  ,mpi_neighbor_list(loop_))
         call echange_voisin(vel_v ,id_v2_  ,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
!...........

      enddo                          !PMXPMX>

      RETURN
      end subroutine q_wavebreaking3d_old
#endif

!..............................

      subroutine q_moteur3d_remapping
      implicit none 

! Je commente ces lignes car elles sont dans add_nhpgf....
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3dr4(i,j,k,id_zt)=sig1d_t(k)*(h_w(i,j)+ssh_int_w(i,j,1))+ssh_int_w(i,j,1)
!     enddo       ; enddo  ; enddo
! Surface et fond
!     do j=1,jmax ; do i=1,imax
!      anyv3dr4(i,j,kmax+1,id_zt)=0.
!      anyv3dr4(i,j,0)=-h_w(i,j)
!     enddo  ; enddo

! rappel: z_t(k+1)-z_t(k)=H*dsig_w(k+1)
! rappel: z_w(k+1)-z_w(k)=H*dsig_t(k)
      
!----------------------------------------------------------------------------------------------
!     if(computematrix==1) then !m[°v°]w>

! sshr4_w a EtE calculEe dans subroutine q_nhpgf2d
       ww_=w2_ ; if(ww_==0.)ww_=1.
       do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
        mc(i,j,k,1)=-ww_*dt2/( dsigr4_t(i,j,k)*dsig_w(i,j,k)*(h_w(i,j)+ssh_int_w(i,j,1))**2 )
       enddo       ; enddo       ;enddo

       do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
        mc(i,j,k,3)=mc(i,j,k,1)*dsig_w(i,j,k)/dsig_w(i,j,k+1)
       enddo       ; enddo       ;enddo

       do k=1,ktop-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! boucle k stoppee en ktop-1 car q(ktop)=0
        mc(i,j,k,2)=1.-mc(i,j,k,1)-mc(i,j,k,3)+damp
       enddo       ; enddo       ;enddo

!.......................................................
! C.L. Surface: q_t(ktop)=0.
       k=ktop
       if(flag_surf_obc_q==0) then !m°v°m> !17-08-18
! CAS NORMAL
! C.L. q_t(ktop)=0.
        do j=2,jmax-1 ; do i=2,imax-1
         mc(i,j,k,1)=0.
         mc(i,j,k,2)=1.
         mc(i,j,k,3)=0.
        enddo       ; enddo
       else                        !m°v°m> !17-08-18
! CAS PARTICULIER surface libre hydrostatique
! C.L. q_t(ktop)=q_t(ktop-1)
        do j=2,jmax-1 ; do i=2,imax-1
         mc(i,j,k,1)=-1.
         mc(i,j,k,2)=1.
         mc(i,j,k,3)=0.
        enddo       ; enddo
       endif                       !m°v°m> !17-08-18

!.......................................................
! C.L. Fond etape 1 (voir aussi plus loin etape 2)
! condition q_t(0)=q_t(1)
       if(flag_merged_levels==0) then !m°v°m>
! cas coordonnee sigma
        k=1
        do j=2,jmax-1 ; do i=2,imax-1
         mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
         mc(i,j,k,1)=0.
        enddo       ; enddo
       else                          !m°v°m>
! cas coordonnee VQS !26-12-21
! Alors attention, il est important de calculer mc(1:3) ici car ca sert
! dans la boucle à suivre, mais pas mc(4)
        do j=2,jmax-1 ; do i=2,imax-1
         k=kmin_w(i,j)
         mc(i,j,k,2)=mc(i,j,k,2)+mc(i,j,k,1)
         mc(i,j,k,1)=0.
         mc(i,j,k,3)=mc(i,j,k,3)*vqsratio_t(i,j)-mc(i,j,k,2)*(1-vqsratio_t(i,j))
! note: vqsratio_t=0 conduit à q(k)=q(k+1)
         do k=1,kmin_w(i,j)-1
          mc(i,j,k,1)=0.
          mc(i,j,k,2)=1.
          mc(i,j,k,3)=-1.
         enddo
        enddo       ; enddo
       endif                         !m°v°m>
! La condition suivante a ete testee pour voir si le schema avait une
! influence sur le courant residuel (pas probant donc je garde l'autre
! methode, plus simple)
! condition oU le gradient nul s'etablit sur un z fixe:
! q_t(k=0)=q(zfixe) ou zfixe est z du niveau 1 au repos
!         =q(k=1)+dq/dz*(zfixe-z(k=1))
!         =q(k=1)+dq/dz*(-sigma*ssh)
!      k=1
!      do j=2,jmax-1 ; do i=2,imax-1
!       mc(i,j,k,2)=mc(i,j,k,2) &
!      +mc(i,j,k,1)*(1.+0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))*ssh_int_w(i,j,1) &
!                          /(depth_t(i,j,k+1)-depth_t(i,j,k)) )
!       mc(i,j,k,3)=mc(i,j,k,3) &
!      +mc(i,j,k,1)*(  -0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))*ssh_int_w(i,j,1) &
!                          /(depth_t(i,j,k+1)-depth_t(i,j,k)) )
!       mc(i,j,k,1)=0.
!      enddo       ; enddo

       checkr2=mc(imax/2,jmax/2,kmax/2,1)
       checkr3=mc(imax/2,jmax/2,kmax/2,2)
       checkr4=mc(imax/2,jmax/2,kmax/2,3)

!      if(nh_frozensigma==1)computematrix=0 ! A la prochaine iteration on ne passera plus par ces lignes

!     endif                     !m[°v°]w>
!----------------------------------------------------------------------------------------------

      flag_stop=0
      if(checkr2/=mc(imax/2,jmax/2,kmax/2,1)) &
                                            flag_stop=1
      if(checkr3/=mc(imax/2,jmax/2,kmax/2,2)) &
                                            flag_stop=1
      if(checkr4/=mc(imax/2,jmax/2,kmax/2,3)) &
                                            flag_stop=1
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'Error 4048 module_q'

!.......................................................
! C.L. Surface mc4 inutile car initialisE A zero
! Cas q(surf)=0
!      k=ktop
!      do j=2,jmax-1 ; do i=2,imax-1
!       mc(i,j,k,4)=0. ! A priori inutile car initialisE A zero
!       mc(i,j,k,4)=grav*ssh_int_w(i,j,2)*(exp(0.04*ssh_int_w(i,j,2))-1.)  
!      enddo       ; enddo
! Cas q(surf)=terme croissance Eq 7 Equation 7 de Rolf Deigaard , Peter Nielsen, Coastal Engineering 139 (2018) 36–46
!      call q_wind_induced_waves

      ww_=1.
      if(w2_/=0.)ww_=1./w2_

! Je calcule q(2) juste apres ssh(2) et apres vel(2) de sorte que je suppose
! l'ordre chrono suivant:
! q(0)   u(1)   q(1)   u(2)  q(2)
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        anyv3dr4(i,j,k,id_hydro_du)=                                         &
       ((vel_u(i,j,k,now)*dz_u(i,j,k,dnow)-vel_u(i,j,k,after)*dz_u(i,j,k,dafter) ) & !20-08-18
            /(dti_fw*0.5*(dz_u(i,j,k,dnow)                   +dz_u(i,j,k,dafter))) &
                      -nhpgf_u(i  ,j,k)*invdx_u(i  ,j))*mask_u(i,j,k)
       enddo       ; enddo         ; enddo
       do k=1,kmax ; do j=2,jmax   ; do i=2,imax-1
        anyv3dr4(i,j,k,id_hydro_dv)=                                         &
       ((vel_v(i,j,k,now)*dz_v(i,j,k,dnow)-vel_v(i,j,k,after)*dz_v(i,j,k,dafter) ) &
            /(dti_fw*0.5*(dz_v(i,j,k,dnow)                   +dz_v(i,j,k,dafter))) &
                      -nhpgf_v(i,j  ,k)*invdy_v(i,j  ))*mask_v(i,j,k)     
       enddo       ; enddo         ; enddo

! id_hydro_du,id_hydro_dv

      do k=1,ktop-1 
      kp1=min(k+1,kmax)
      km1=max(k-1,1   )
      do j=2,jmax-1 ; do i=2,imax-1

        mc(i,j,k,4)=2*q_t(i,j,k,1)-q_t(i,j,k,0)      &

            -(q_t(i,j,k+1,1)-  q_t(i,j,k-1,1))       &!b3-time
        /(depth_t(i,j,k+1)-depth_t(i,j,k-1))         &!b3-time
        *(1.+sig1dpom_t(i,j,k))*(2*ssh_int_w(i,j,1)-ssh_int_w(i,j,0)-ssh_int_w(i,j,2)) &!b3-time

! HYDROSTATIC:
      +dt2*( & !hydro>
       invdx_t(i,j)*( & !xxx>
         anyv3dr4(i+1,j,k,id_hydro_du)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        +(depth_t(i,j,k)-depth_u(i+1,j,k))*(anyv3dr4(i+1,j,kp1,id_hydro_du)-anyv3dr4(i+1,j,km1,id_hydro_du)) &
                                          /( depth_u(i+1,j,kp1)             -depth_u(i+1,j,km1)            ) &
        -anyv3dr4(  i,j,k,id_hydro_du)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        -(depth_t(i,j,k)-depth_u(i  ,j,k))*(anyv3dr4(i  ,j,kp1,id_hydro_du)-anyv3dr4(i  ,j,km1,id_hydro_du)) &
                                          /( depth_u(i  ,j,kp1)             -depth_u(i  ,j,km1)            ) &
                    ) & !xxx>
      +invdy_t(i,j)*( & !yyy>
         anyv3dr4(i,j+1,k,id_hydro_dv)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        +(depth_t(i,j,k)-depth_v(i,j+1,k))*(anyv3dr4(i,j+1,kp1,id_hydro_dv)-anyv3dr4(i,j+1,km1,id_hydro_dv)) &
                                          /( depth_v(i,j+1,kp1)             -depth_v(i,j+1,km1)            ) &
        -anyv3dr4(i,j  ,k,id_hydro_dv)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        -(depth_t(i,j,k)-depth_v(i,j  ,k))*(anyv3dr4(i,j  ,kp1,id_hydro_dv)-anyv3dr4(i,j  ,km1,id_hydro_dv)) &
                                          /( depth_v(i,j  ,kp1)             -depth_v(i,j  ,km1)            ) &
                    ) & !yyy>
            ) & !hydro>


                                              +dt2*( & !ttt>

       invdx_t(i,j)*( & !ooo>

                ( q_t(i+1,j,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!               +(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b1
!          /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k-1)) & !b1
!      *0.5*(depth_t(i-1,j,k  )-depth_t(i+1,j,k  )) & !b1

!               +(q_t(i+1,j,k+1,1)-     q_t(i+1,j,k-1,1)) & !b2
!          /(anyv3dr4(i+1,j,k+1  ,id_zt)-anyv3dr4(i+1,j,k-1  ,id_zt)) & !b2
!          *(depth_t(i  ,j,k  )-depth_t(i+1,j,k  )) & !b2

               +(q_t(i  ,j,k+1,1)-  q_t(i  ,j,k-1,1)) & !b3
           /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k-1))   & !b3
           *(depth_t(i  ,j,k  )-depth_t(i+1,j,k  ))   & !b3

                                -q_t(i  ,j,k,1) )*invdx_u(i+1,j)    &
                                                  *mask_u(i+1,j,k)  &

               -( q_t(i  ,j,k,1)-q_t(i-1,j,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!                              -(q_t(i  ,j,k+1,1)-     q_t(i  ,j,k-1,1)) & !b1
!                         /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k-1)) & !b1
!                     *0.5*(depth_t(i+1,j,k  )-depth_t(i-1,j,k  )) & !b1

!                              -(q_t(i-1,j,k+1,1)-     q_t(i-1,j,k-1,1)) & !b2
!                         /(anyv3dr4(i-1,j,k+1  ,id_zt)-anyv3dr4(i-1,j,k-1  ,id_zt)) & !b2
!                         *(depth_t(i  ,j,k  )-depth_t(i-1,j,k  )) & !b2

                              -(q_t(i  ,j,k+1,1)-  q_t(i  ,j,k-1,1)) & !b3
                          /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k-1))   & !b3
                          *(depth_t(i  ,j,k  )-depth_t(i-1,j,k  ))   & !b3

                                                )*invdx_u(i  ,j)    &
                                                  *mask_u(i  ,j,k)  &

                    ) & !ooo>                                           

      +invdy_t(i,j)*( & !ooo>

                ( q_t(i,j+1,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!               +(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b1
!          /(depth_t(i,j  ,k+1)-depth_t(i,j  ,k-1)) & !b1
!      *0.5*(depth_t(i,j-1,k  )-depth_t(i,j+1,k  )) & !b1

!               +(q_t(i,j+1,k+1,1)-     q_t(i,j+1,k-1,1)) & !b2
!          /(anyv3dr4(i,j+1,k+1  ,id_zt)-anyv3dr4(i,j+1,k-1  ,id_zt)) & !b2
!          *(depth_t(i,j  ,k  )-depth_t(i,j+1,k  )) & !b2

               +(q_t(i,j  ,k+1,1)-  q_t(i,j  ,k-1,1)) & !b3
           /(depth_t(i,j  ,k+1)-depth_t(i,j  ,k-1))   & !b3
           *(depth_t(i,j  ,k  )-depth_t(i,j+1,k  ))   & !b3

                                -q_t(i  ,j,k,1) )*invdy_v(i,j+1)    &
                                                  *mask_v(i,j+1,k)  &

               -( q_t(i  ,j,k,1)-q_t(i,j-1,k,1) &
! Le schema b1 ne fonctionne pas si point continental voisin car z(i+1)-z(i-1) indefini
!                              -(q_t(i,j  ,k+1,1)-     q_t(i,j  ,k-1,1)) & !b1
!                         /(depth_t(i,j  ,k+1)-depth_t(i,j  ,k-1)) & !b1
!                     *0.5*(depth_t(i,j+1,k  )-depth_t(i,j-1,k  )) & !b1

!                              -(q_t(i,j-1,k+1,1)-     q_t(i,j-1,k-1,1)) & !b2
!                         /(anyv3dr4(i,j-1,k+1  ,id_zt)-anyv3dr4(i,j-1,k-1  ,id_zt)) & !b2
!                         *(depth_t(i,j  ,k  )-depth_t(i,j-1,k  )) & !b2

                              -(q_t(i,j  ,k+1,1)-  q_t(i,j  ,k-1,1)) & !b3
                          /(depth_t(i,j  ,k+1)-depth_t(i,j  ,k-1))   & !b3
                          *(depth_t(i,j  ,k  )-depth_t(i,j-1,k  ))   & !b3

                                                )*invdy_v(i  ,j)    &
                                                  *mask_v(i  ,j,k)  &

                    ) & !ooo>

                                                   ) &  !ttt>


       -mc(i,j,k,1)    *ww_*(w0_*q_t(i,j,k-1,0)+w1_*q_t(i,j,k-1,1) )&
       -mc(i,j,k,3)    *ww_*(w0_*q_t(i,j,k+1,0)+w1_*q_t(i,j,k+1,1) )&
      -(mc(i,j,k,2)-1.)*ww_*(w0_*q_t(i,j,k  ,0)+w1_*q_t(i,j,k  ,1) ) 

! Note sur sponge_t: varie de 0 a 1 (sur la frontiere) ce qui equivaut en termes
! d'echelle de temps associee a un temps de rappel equivalent au pas de temps. C'est
! donc un rappel fort qui a pour but de supprimer la pression NH aux frontieres ouvertes

      enddo ; enddo
      enddo
      if(flag_merged_levels==1) then !m°v°m>
! Condition sur mc4 cas coordonnee VQS !26-12-21
! Alors attention, il est important de calculer mc(1:3) ici car ca sert
! dans la boucle à suivre, mais pas mc(4)
        do j=2,jmax-1 ; do i=2,imax-1
         k=kmin_w(i,j)
         mc(i,j,k,4)=mc(i,j,k,4)*vqsratio_t(i,j)
         do k=1,kmin_w(i,j)-1
          mc(i,j,k,4)=0.
         enddo
        enddo       ; enddo
       endif                         !m°v°m>

! Solveur vertical tridiagonal:
      call q_solver(computefullsolver,2,imax-1,2,jmax-1,1,ktop) ! donne q_t(t+1)
!     i=0.66*imax
!     j=jmax/2
!     do k=ktop,1,-1
!      write(6,*)k,q_t(i,j,k,2)
!     enddo
!     stop 'coco1'
!     if(nh_frozensigma==1)computefullsolver=0 ! A la prochaine iteration une partie seulement du solver sera calculEe

!........................
! ZONES EPONGE ET WETMASK: (Ne pourrait elle pas etre implicitee ou simplifiee?)
      do j=2,jmax-1 ; do i=2,imax-1
       anyvar2d(i,j)=min(1.,(h_w(i,j)+ssh_int_w(i,j,1))/wetdry_cstnh) & !26-08-18
                    *mask_t(i  ,j  ,kmax) &
                    *mask_t(i+1,j  ,kmax) &
                    *mask_t(i-1,j  ,kmax) &
                    *mask_t(i  ,j-1,kmax) &
                    *mask_t(i  ,j+1,kmax)
      enddo         ; enddo
!     do k=0,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1
      do k=1,kmax-1 ; do j=2,jmax-1 ; do i=2,imax-1 ! Boucle k commence A 1 si report des C.L. en suivant....
!      q_t(i,j,k,2)=q_t(i,j,k,2)*(1.-sponge_t(i,j,1))*wetmask_t(i,j)
       q_t(i,j,k,2)=q_t(i,j,k,2)                     *anyvar2d(i,j) !26-08-18
!      if(k==kmax/2.and.j==jmax/2.and.i==imax/2)write(10+par%rank,*)iteration3d,q_t(i,j,k,2)
      enddo         ; enddo       ; enddo

! BIDOUILLE PATRICK POUR RETROUVER UN CODE HYDROSTATIQUE
!!!!!!q_t=0.


!........................
! Conditions ouvertes radiatives
      call q_obc_radiative

!........................
! C.L. Fond etape 2
! Reporter les hypotheses de conditions aux limites au fond et en surface:
!     if(flag_merged_levels==0) then !m°v°m>
! Cas coordonnee sigma:
       do j=2,jmax-1 ; do i=2,imax-1
! C.L.: q(0)=q(1):
        q_t(i,j,0,2)=q_t(i,j,1,2)*mask_t(i,j,kmax)
! C.L.: q(0)=q(zfixe):
!       q_t(i,j,0,2)=mask_t(i,j,kmax)*( &
!       q_t(i,j,1,2)-(    q_t(i,j,2,2)  -q_t(i,j,1,2)) &
!                   /(depth_t(i,j,2)-depth_t(i,j,1)  ) &
!               *0.5*(sigma_w(i,j,2)+sigma_w(i,j,1))*ssh_int_w(i,j,1) )

! C.L.: q(kmax+1)=0.
!       q_t(i,j,ktop:kmax+1,2)=0. J'ai commentE cette ligne car elle me semble redondante
       enddo ; enddo
!     else                          !m°v°m>
! Cas coordonnee VQS: !26-12-21
!      do j=2,jmax-1 ; do i=2,imax-1
!       do k=kmerged_t(i,j)-1,0,-1
!        q_t(i,j,k,2)=q_t(i,j,k+1,2)*mask_t(i,j,kmax)
!       enddo
!      enddo ; enddo
!     endif                         !m°v°m>


!........................
! Continuite MPI:
      call q_obc_mpi('z0') ! mpi 'z1' sur q_t(t+1)


#ifdef bidon
! VERIFICATION DE L'EQUILIBRE:
      if(modulo(iteration3d,1000)==0) then !>>>
! Ici verification possible de l'equilibre
! Attention que le filtre d'asselin se fait ulterieurement
! et par consequent n'est pas evalue ici
!     i=imax/2 ; j=jmax/2 ; k=kmax/2
      write(texte30,'(a,i0)')'tmp/nhdiag_balance',par%rank

      if(.not.allocated(nhconvergence_histo1)) then
         allocate(nhconvergence_histo1(900:1000)) ; nhconvergence_histo1=0
         allocate(nhconvergence_histo2(900:1000)) ; nhconvergence_histo2=0
      endif

! reset A chaque iteration, ou pas, au choix....
!     nhconvergence_histo1=0
!     nhconvergence_histo2=0
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax-1

! opposE du Terme de time stepping:
       x1= -q_t(i,j,k,2)+2*q_t(i,j,k,1)-q_t(i,j,k,0)             &
             -(q_t(i,j,k+1,1)-  q_t(i,j,k-1,1))                   &!b3-time
         /(depth_t(i,j,k+1)-depth_t(i,j,k-1))                     &!b3-time
        *(1.+sig1dpom_t(i,j,k))*(2*ssh_int_w(i,j,1)-ssh_int_w(i,j,0)-ssh_int_w(i,j,2))  !b3-time
! Forcage hydrostatique:
       x2= & !anyvar2d(i,j)
!      dt2*( invdx_t(i,j)*(anyv3dr4(i+1,j,k,id_hydro_du)-anyv3dr4(i,j,k,id_hydro_du))   &
!           +invdy_t(i,j)*(anyv3dr4(i,j+1,k,id_hydro_dv)-anyv3dr4(i,j,k,id_hydro_dv)) )  
      +dt2*( & !hydro> !18-08-18
       invdx_t(i,j)*( & !xxx>
         anyv3dr4(i+1,j,k,id_hydro_du)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        +(depth_t(i,j,k)-depth_u(i+1,j,k))*(anyv3dr4(i+1,j,kp1,id_hydro_du)-anyv3dr4(i+1,j,km1,id_hydro_du)) &
                                          /( depth_u(i+1,j,kp1)             -depth_u(i+1,j,km1)            ) &
        -anyv3dr4(  i,j,k,id_hydro_du)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        -(depth_t(i,j,k)-depth_u(i  ,j,k))*(anyv3dr4(i  ,j,kp1,id_hydro_du)-anyv3dr4(i  ,j,km1,id_hydro_du)) &
                                          /( depth_u(i  ,j,kp1)             -depth_u(i  ,j,km1)            ) &
                    ) & !xxx>
      +invdy_t(i,j)*( & !yyy>
         anyv3dr4(i,j+1,k,id_hydro_dv)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        +(depth_t(i,j,k)-depth_v(i,j+1,k))*(anyv3dr4(i,j+1,kp1,id_hydro_dv)-anyv3dr4(i,j+1,km1,id_hydro_dv)) &
                                          /( depth_v(i,j+1,kp1)             -depth_v(i,j+1,km1)            ) &
        -anyv3dr4(i,j  ,k,id_hydro_dv)  & ! suivi de 2 lignes de remapping sur niveau depth_t:
        -(depth_t(i,j,k)-depth_v(i,j  ,k))*(anyv3dr4(i,j  ,kp1,id_hydro_dv)-anyv3dr4(i,j  ,km1,id_hydro_dv)) &
                                          /( depth_v(i,j  ,kp1)             -depth_v(i,j  ,km1)            ) &
                    ) & !yyy>
            )   !hydro>
! Laplacien 3D de la pression NH deduit de l'equilibre elementaire
! m1*q(k-1,t+1)+m2*q(k  ,t+1)+m3*q(k+1,t+1)=m4 que l'on a pris
! soin prealablement de verifier
       x3=-q_t(i,j,k,2)+mc(i,j,k,4)-x1-x2       &
                 -mc(i,j,k,1)   *q_t(i,j,k-1,2) &
                -(mc(i,j,k,2)-1)*q_t(i,j,k  ,2) &
                 -mc(i,j,k,3)   *q_t(i,j,k+1,2)  

! x0=poid proportionnel au courant (si pas de courant on se fiche de savoir
! le % de nonhydrostaticite)
       x0=vel_u(i  ,j  ,k,1)**2 &
         +vel_u(i+1,j  ,k,1)**2 &
         +vel_v(i  ,j  ,k,1)**2 &
         +vel_v(i  ,j+1,k,1)**2  

       sum1=sum1+mask_t(i,j,k)*x0
       sum2=sum2+mask_t(i,j,k)*x0*(1.-abs(x2+x3)/(max(abs(x2)+abs(x3),small1)))

!      k10=nint( 100.*(1.-abs(x2+x3)/(max(abs(x2)+abs(x3),small1))) )
!      if(k10<0.or.k10>100)stop 'k10<0.or.k10>100'

       k10=nint( 1000.*(1.-abs(x2+x3)/(max(abs(x2)+abs(x3),small1))) )
       if(k10<0.or.k10>1000)stop 'k10<0.or.k10>1000'
       k10=max(k10,900)
      
       nhconvergence_histo1(k10)=nhconvergence_histo1(k10)+1

      enddo
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(nhconvergence_histo1,nhconvergence_histo2,101,mpi_integer,mpi_sum,par%comm2d,ierr)

#ifdef bidon
      open(unit=3,file=texte30,position='append')
! Equilibre elementaire qui doit etre verifiE:
!      write(3,*)mc(i,j,k,1)*q_t(i,j,k-1,2) &
!               +mc(i,j,k,2)*q_t(i,j,k  ,2) &
!               +mc(i,j,k,3)*q_t(i,j,k+1,2) &
!               ,mc(i,j,k,4)
! l'addition de 2+3+4 egale zero:
       write(3,*)real(iteration3d*dti_fw/3600.) & ! 1 temps en heure
                ,real(x1)                       & ! 2 terme de time-stepping (ou desequilibre)
                ,real(x2)                       & ! 3 forcage hydrostatique
                ,real(x3)                         ! 4 Laplacien 3D non-hydrostatique
! Si le solveur a bien convergE, 2 doit etre petit par rapport a 3 et 4
      close(3)
#endif

       if(par%rank==0) then
        open(unit=3,file='tmp/nhdiag_convergence',position='append')
         write(3,*)real(iteration3d*dti_fw/3600.) & ! 1 temps en heure
                   ,sum2glb/max(sum1glb,small1)     ! Convergence du solveur: 1=parfait 0=echec
        close(3)                                    
        open(unit=3,file='tmp/nhdiag_convhisto',position='append')
         write(3,*)'---- iteration3: ',iteration3d
         do k=900,1000
          write(3,*)real(k)*0.1,nhconvergence_histo2(k)
         enddo
        close(3)                                    
       endif

      endif                               !>>>
#endif

      end subroutine q_moteur3d_remapping

!..............................

      subroutine q_nhpgf3d
      implicit none

! Calculer    dq/di-dz/di*dq/dz:
!     do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
!       nhpgf_u(i,j,k)= q_t(i,j,k,1)  -q_t(i-1,j,k,1)                  &
!                 -(depth_t(i,j,k)-depth_t(i-1,j,k))                   &
!               *0.25*( q_t(i,j,k+1,1)+q_t(i-1,j,k+1,1)                &
!                      -q_t(i,j,k-1,1)-q_t(i-1,j,k-1,1))/dz_u(i,j,k,1)  
!     enddo ; enddo ; enddo

! Calculer    dq/dj-dz/dj*dq/dz:
!     do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
!       nhpgf_v(i,j,k)= q_t(i,j,k,1)  -q_t(i,j-1,k,1)                  &
!                 -(depth_t(i,j,k)-depth_t(i,j-1,k))                   &
!               *0.25*( q_t(i,j,k+1,1)+q_t(i,j-1,k+1,1)                &
!                      -q_t(i,j,k-1,1)-q_t(i,j-1,k-1,1))/dz_v(i,j,k,1)  
!     enddo ; enddo ; enddo

! Meme si c'est identique (car le moveforward faisait que q(1)=q(2)
! je prefere cette ecriture qui rappelle que l'ordre 3d est le suivant:
! q(0)  u(0)  q(1)  u(1)  q(2) autrement dit que l'on passe de u(1)
! a u(2) en utilisant pgf(q(2))

! Calculer    dq/di-dz/di*dq/dz:
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
#ifdef bidon
! methode 1
        nhpgf_u(i,j,k)=( & ! -)>--
                        q_t(i,j,k,1)  -q_t(i-1,j,k,1)    & !20-08-18
                  -(depth_t(i,j,k)-depth_t(i-1,j,k))     &
                *0.25*( q_t(i,j,k+1,1)+q_t(i-1,j,k+1,1)  &
                       -q_t(i,j,k-1,1)-q_t(i-1,j,k-1,1))/dz_u(i,j,k,dnow) &
                       ) & ! -)>--
                        *nhpgf_reduce !20-11-17
#endif
!#ifdef bidon
! methode 2 !02-05-20
        nhpgf_u(i,j,k)=( & ! -)>--
                        q_t(i,j,k,1)  -q_t(i-1,j,k,1)    & !20-08-18
                  -(depth_t(i,j,k)-depth_t(i-1,j,k))     &
      *(    q_t(i,j,k+1,1)  +q_t(i-1,j,k+1,1)  -q_t(i,j,k-1,1)  -q_t(i-1,j,k-1,1)) &
      /(depth_t(i,j,k+1)+depth_t(i-1,j,k+1)-depth_t(i,j,k-1)-depth_t(i-1,j,k-1)) &
                       ) & ! -)>--
                        *nhpgf_reduce !20-11-17
!#endif
#ifdef bidon
! methode 3 !25-05-20
        nhpgf_u(i,j,k)=( & ! -)>--
                        q_t(i,j,k,1)  -q_t(i-1,j,k,1)    & !20-08-18
                  -(depth_t(i,j,k)-depth_t(i-1,j,k))     &
      *0.5*( & !ooo>
                  (    q_t(i  ,j,k+1,1)  -q_t(i  ,j,k-1,1)) &
                 /(depth_t(i  ,j,k+1)-depth_t(i  ,j,k-1)  ) &
                 +(    q_t(i-1,j,k+1,1)  -q_t(i-1,j,k-1,1)) &
                 /(depth_t(i-1,j,k+1)-depth_t(i-1,j,k-1)  ) &
           ) & !ooo>
                       ) & ! -)>--
                        *nhpgf_reduce !20-11-17
#endif
      enddo ; enddo ; enddo

! Calculer    dq/dj-dz/dj*dq/dz:
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
#ifdef bidon
! methode 1
        nhpgf_v(i,j,k)=( & ! -)>--
                        q_t(i,j,k,1)  -q_t(i,j-1,k,1)    & !20-08-18
                  -(depth_t(i,j,k)-depth_t(i,j-1,k))     &
                *0.25*( q_t(i,j,k+1,1)+q_t(i,j-1,k+1,1)  &
                       -q_t(i,j,k-1,1)-q_t(i,j-1,k-1,1))/dz_v(i,j,k,dnow) &
                       ) & ! -)>--
                        *nhpgf_reduce !20-11-17
#endif
!#ifdef bidon
! methode 2 !02-05-20
        nhpgf_v(i,j,k)=( & ! -)>--
                        q_t(i,j,k,1)  -q_t(i,j-1,k,1)    & !20-08-18
                  -(depth_t(i,j,k)-depth_t(i,j-1,k))     &
      *(    q_t(i,j,k+1,1)  +q_t(i,j-1,k+1,1)  -q_t(i,j,k-1,1)  -q_t(i,j-1,k-1,1)) &
      /(depth_t(i,j,k+1)+depth_t(i,j-1,k+1)-depth_t(i,j,k-1)-depth_t(i,j-1,k-1)) &
                       ) & ! -)>--
                        *nhpgf_reduce !20-11-17
!#endif
#ifdef bidon
! methode 3 !25-05-20
        nhpgf_v(i,j,k)=( & ! -)>--
                        q_t(i,j,k,1)  -q_t(i,j-1,k,1)    & !20-08-18
                  -(depth_t(i,j,k)-depth_t(i,j-1,k))     &
      *0.5*( & !ooo>
                  (    q_t(i,j  ,k+1,1)  -q_t(i,j  ,k-1,1)) &
                 /(depth_t(i,j  ,k+1)-depth_t(i,j  ,k-1)  ) &
                 +(    q_t(i,j-1,k+1,1)  -q_t(i,j-1,k-1,1)) &
                 /(depth_t(i,j-1,k+1)-depth_t(i,j-1,k-1)  ) &
           ) & !ooo>

                       ) & ! -)>--
                        *nhpgf_reduce !20-11-17
#endif
      enddo ; enddo ; enddo

      if(flag_merged_levels==1) then !pmx> !26-12-21

!     if(iteration3d==2) then
!     i=imax/2 ; j=jmax/2
!     do k=kmax,kmerged_u(i,j),-1
!      write(10+par%rank,*)depth_u(i,j,k),nhpgf_u(i,j,k)
!     enddo
!     do k=kmerged_u(i,j)-1,1,-1
!      write(110+par%rank,*)depth_u(i,j,k),nhpgf_u(i,j,k),' **'
!     enddo
!     endif

! Cas des coordonnees VQS
       do j=2,jmax-1 ; do i=2,imax
        do k=1,kmerged_u(i,j)-1
         nhpgf_u(i,j,k)= &
!        nhpgf_u(i,j,k)*pgfratio_u(i,j)+(1.-pgfratio_u(i,j))* &
         nhpgf_u(i,j,kmerged_u(i,j))
        enddo
       enddo         ; enddo
       do j=2,jmax ; do i=2,imax-1
        do k=1,kmerged_v(i,j)-1
         nhpgf_v(i,j,k)= &
!        nhpgf_v(i,j,k)*pgfratio_v(i,j)+(1.-pgfratio_v(i,j))* &
         nhpgf_v(i,j,kmerged_v(i,j))
        enddo
       enddo         ; enddo

!     if(iteration3d==2) then
!     i=imax/2 ; j=jmax/2
!     do k=kmax,kmerged_u(i,j),-1
!      write(210+par%rank,*)depth_u(i,j,k),nhpgf_u(i,j,k)
!     enddo
!     do k=kmerged_u(i,j)-1,1,-1
!      write(310+par%rank,*)depth_u(i,j,k),nhpgf_u(i,j,k),' **'
!     enddo
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'bibi'
!      endif

      endif                          !pmx> !26-12-21

      if(flag_timesplitting==1) then !w°v°w>

      if(flag_surf_obc_q==1) then !-dq/dz=0->
! Avec l'option AVEC time splitting et la condition dq/dz=0
! la surface libre est hydrostatique
! On ote par consequent la moyenne verticale du PGF NH:

       do j=2,jmax-1 ; do i=2,imax
        anyvar2d(i,j)=0.
       enddo         ; enddo 
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        anyvar2d(i,j)= &
        anyvar2d(i,j)+nhpgf_u(i,j,k)*dz_u(i,j,k,1) !05-06-20
       enddo ; enddo ; enddo
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        nhpgf_u(i,j,k)= &
        nhpgf_u(i,j,k)-anyvar2d(i,j)     &
                          /hz_u(i,j,1)     !05-06-20
       enddo ; enddo ; enddo

       do j=2,jmax ; do i=2,imax-1
        anyvar2d(i,j)=0.
       enddo       ; enddo 
       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        anyvar2d(i,j)= &
        anyvar2d(i,j)+nhpgf_v(i,j,k)*dz_v(i,j,k,1) !05-06-20
       enddo ; enddo ; enddo
       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        nhpgf_v(i,j,k)= &
        nhpgf_v(i,j,k)-anyvar2d(i,j)  &
                          /hz_v(i,j,1)     !05-06-20
       enddo ; enddo ; enddo

      endif                       !-dq/dz=0->>

      if(flag_surf_obc_q==0) then !-q=0->
! Avec l'option AVEC time splitting et la condition q=0, la surface
! libre est non-hydrostatique mais NHPGFbar est filtrE par EMA

! coef du filtre EMA:
!      x1_r4=0.1  ; x2_r4=1-x1_r4
       x1_r4=1.   ; x2_r4=1-x1_r4 ! x1_r4=1 --> CAS PARTICULIER SANS EMA

! nhpgf2d_u:
       do j=2,jmax-1 ; do i=2,imax
        anyvar2d(i,j)=0.
       enddo         ; enddo 
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        anyvar2d(i,j)= &
        anyvar2d(i,j)+nhpgf_u(i,j,k)*dz_u(i,j,k,1) 
       enddo ; enddo ; enddo
       do j=2,jmax-1 ; do i=2,imax
        nhpgf2d_u(i,j)=                               &
        nhpgf2d_u(i,j)*x2_r4                          &
                      +x1_r4*anyvar2d(i,j)/hz_u(i,j,1)
       enddo         ; enddo 
! pres3d2d_u:
       do j=2,jmax-1 ; do i=2,imax
        pres3d2d_u(i,j)= &
        pres3d2d_u(i,j)  &
        +nhpgf2d_u(i,j)  !*0.1 ! bidouille patrick
       enddo         ; enddo 

! nhpgf2d_v:
       do j=2,jmax ; do i=2,imax-1
        anyvar2d(i,j)=0.
       enddo       ; enddo 
       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        anyvar2d(i,j)= &
        anyvar2d(i,j)+nhpgf_v(i,j,k)*dz_v(i,j,k,1) 
       enddo ; enddo ; enddo
       do j=2,jmax ; do i=2,imax-1
        nhpgf2d_v(i,j)=                               &
        nhpgf2d_v(i,j)*x2_r4                          &
                      +x1_r4*anyvar2d(i,j)/hz_v(i,j,1)
       enddo       ; enddo 
! pres3d2d_v:
       do j=2,jmax ; do i=2,imax-1
        pres3d2d_v(i,j)= &
        pres3d2d_v(i,j)  &
        +nhpgf2d_v(i,j)  !*0.1 ! bidouille patrick
       enddo         ; enddo 

      endif                       !-q=0->

      endif                       !w°v°w>

      end subroutine q_nhpgf3d

!............................

      subroutine q_periode_0diag
      implicit none
      integer :: iteration_first0=-999  &
                ,zero_counter_=-3 ! Commence A -3 pour laisser passer spin-up
      double precision :: deduced_period=-999.

! period_=2.*pi/(kvector_*c_)
      i=imax/2 ; j=jmax/2

      if(ssh_int_w(i,j,0)*ssh_int_w(i,j,1)<0.) then !PMX>

       write(texte30,'(a,i0)')'tmp/nhdiag_period',par%rank
       zero_counter_=zero_counter_+1

       if(zero_counter_==0) then !reset
          if(nhpgf_reduce/=0.) then !>>>
           x2=sqrt( grav/kvector_ *tanh(kvector_*hmax) )                    ! vitesse de phase theorique
          else                      !>>>
           x2=sqrt( grav*hmax )                                             ! vitesse de phase theorique
          endif                     !>>>
          x3=2.*pi/(kvector_*x2)                                            ! periode theorique
          iteration_first0=iteration3d
          open(unit=3,file=texte30,access='append')
           write(3,*)cfl_reduce        ,'  cfl_reduce'
           write(3,*)nhpgf_reduce      ,' nhpgf_reduce'
           write(3,*)x2                ,' c theorique'
           write(3,*)x3                ,'periode theorique'
           write(3,*)2.*pi/kvector_    ,' longueur d onde'
           write(3,*)2.*pi/kvector_/dxb,' longueur d onde/dxb' 
           write(3,*)cfl_nh            ,' cfl_nh'
           write(3,*)sqrt(dt2)/dte_fw/sqrt(grav*hmax) &
           ,' Cq/sqrt(grav*hmax)'
          close(3)
       endif                     !reset

       if(zero_counter_>0.and.modulo(zero_counter_,2)==0) then !m0v0m>
        open(unit=3,file=texte30,access='append')
         x0=real(iteration3d-iteration_first0)*dti_fw*2/real(zero_counter_) ! periode
         x1=2.*pi/(kvector_*x0)                                             ! vitesse de phase
!        write(3,*)real(elapsedtime_now),real(x0),real(x3),real(x1),real(x2) 
         write(3,*)real(elapsedtime_now),real(x0),         real(x1)
!        if(abs(x0-deduced_period)<0.005) then !>>>
!         write(10+par%rank,*)zero_counter_,x0,deduced_period,x0-deduced_period    
!         stop ' deduced_period a convergE voir fichier fort'
!        endif                                 !>>>
         deduced_period=x0
        close(3)
       endif                                                    !m0v0m>

      endif                                             !PMX>

      end subroutine q_periode_0diag
!............................

      subroutine q_phasespeeddiag
      implicit none


      if(flag_nh3d/=0) then !3Dcase>

! Cas simple monochromatique, unidirectionnel Oi, resolution constante dxb
! pulsation=(d2(ssh)dt2)/ssh
! kvector  =(d2(ssh)dx2)/ssh
       sum1=0.
       sum2=0.
       sum3=0.
       sum4=0.
       j=jmax/2
       do i=2,imax-1

        x1=(ssh_int_w(i+1,j,0)-ssh_int_w(i  ,j,0))      &
          +(ssh_int_w(i+1,j,1)-ssh_int_w(i  ,j,1))         
        x2=(ssh_int_w(i+1,j,1)-ssh_int_w(i+1,j,0))      & 
          +(ssh_int_w(i  ,j,1)-ssh_int_w(i  ,j,0))         
       x11=(ssh_int_w(i+1,j,1)-ssh_int_w(i  ,j,1))      &
          +(ssh_int_w(i+1,j,2)-ssh_int_w(i  ,j,2))         
       x22=(ssh_int_w(i+1,j,2)-ssh_int_w(i+1,j,1))      & 
          +(ssh_int_w(i  ,j,2)-ssh_int_w(i  ,j,1))         

        sum1=sum1+x1**2+x11**2
        sum2=sum2+x2*x1+x22*x11

        x1=ssh_int_w(i,j,2)-2.*ssh_int_w(i  ,j,1)+ssh_int_w(i,j,0)
        x2=ssh_int_w(i,j,1) ! +ssh_int_w(i  ,j,2)+ssh_int_w(i,j,0)
        sum3=sum3+x1**2
        sum4=sum4+x2*x1

       enddo

       sum1=sum1/dxb
       sum2=sum2/dti_fw

       x1=iteration3d*dti_fw/3600
      endif                    !3Dcase>
!     if(flag_nh2d   ==1) then !2Dcase>
!       stop ' 2DCASE A FAIRE'
!     endif                    !2Dcase>

      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum4,sum4glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)

      if(par%rank==0) then !>>>
       open(unit=3,file='tmp/phasespeed',position='append')
        write(3,*)real(x1)                   &                      ! temps (minutes)
                 ,-sum2glb/sum1glb           &                      ! vitesse de phase
!                ,2.*pi/sqrt(abs(sum3glb/sum4glb*3./dti_fw/dti_fw)) ! periode moyenne
                 ,2.*pi/sqrt(abs(sum3glb/sum4glb   /dti_fw/dti_fw)) ! periode moyenne
       close(3)
      endif                !>>>

      end subroutine q_phasespeeddiag

!............................

      subroutine q_rwavebreak
      implicit none

      stop 'Passe pas par q_rwavebreak'

! https://docs.google.com/document/d/1wU5DVBww6YuYZ2VK0mAmTfH0rNPI97AOF46PuEIfFso/edit
! <Rb>=1.5*sqrt(grav/H))*sqrt( (dH/dx)**2 +(dH/dy)**2 )

! Calcul du coef de freinage HF
      do j=1,jmax ; do i=1,imax

       rwavebreak_t(i,j)=1.5*sqrt( & !pmx>
                   grav/max(h_w(i,j),0.03)*( & !0v0
        (0.5*(h_w(i+1,j)-h_w(i-1,j))*invdx_t(i,j))**2 &
       +(0.5*(h_w(i,j+1)-h_w(i,j-1))*invdy_t(i,j))**2 &
                                           ) & !0v0
                                 )   !pmx>

      enddo       ; enddo

      end subroutine q_rwavebreak

!..............................

      subroutine q_wavebreaking3d
      implicit none
      integer :: loop_,id_u1_,id_u2_,id_v1_,id_v2_ &
                ,looptime_,id_tken_                &
                ,looptimemax_=1
      real dti_wb_


! Sous pas de temps pour le wave breaking
      dti_wb_=real(dti_lp)/looptimemax_



!     RETURN

      id_uflx=3    ! xy_t
      id_ufly=3    ! xy_f
      id_vflx=4    ! xy_f
      id_vfly=4    ! xy_t
      id_breaker=5 ! xy_t xyf

! DEFERLEMENT PAR DIFFUSION HORIZONTALE MAXIMALE
      if(iteration2d/=iteration2d_begin) then

        stop 'JE SUIS LA!!!'

        flag_stop=0
        if(checkxyt(2)/=xy_t(imax/2,jmax/2,id_uflx)) &
          flag_stop=1

        if(checkxyt(5)/=xy_t(imax/2,jmax/2,id_vfly)) &
          flag_stop=1

        if(checkxyt(4)/=xy_t(imax/2,jmax/2,id_breaker)) &
          flag_stop=1

        if(checkxyf(1)/=xy_f(imax/2,jmax/2,id_ufly)) &
          flag_stop=1

        if(checkxyf(2)/=xy_f(imax/2,jmax/2,id_vflx)) &
          flag_stop=1

        if(checkxyf(3)/=xy_f(imax/2,jmax/2,id_breaker)) &
          flag_stop=1

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'Error 3705 module_q'

      endif

      do j=1,jmax ; do i=1,imax

          xy_t(i,j,id_breaker)=min( & !MIN> !02-05-20
                                   dx_t(i,j)**2                       & !1 
                                  ,dy_t(i,j)**2                       & !2
                                  ,dti_wb_*wavebreakfactor            & !3
!                                         *wavebreakcoef              & !3
!                                         *sqrt(max(h_w(i,j)*0.6,0.0001)) & !3
!                                         *sqrt(max(h_w(i,j),0.0001)) & !3
!                                         *sqrt(0.1) & !3
                                          *breaker2d_t(i,j)         & !3
                                  )   !MIN>

      enddo       ; enddo
      checkxyt(4)=xy_t(imax/2,jmax/2,id_breaker)

      do j=2,jmax ; do i=2,imax

! pas de diffusion de u en y ni de v en x:
       xy_f(i,j,id_breaker)=0.25*( xy_t(i  ,j  ,id_breaker) &
                                  +xy_t(i-1,j  ,id_breaker) &
                                  +xy_t(i-1,j-1,id_breaker) &
                                  +xy_t(i  ,j-1,id_breaker) )

      enddo       ; enddo
      checkxyf(3)=xy_f(imax/2,jmax/2,id_breaker)

      do looptime_=1,looptimemax_    !PMXPMX>

! ContinuitE mpi deplacEe en debut de boucle le !01-05-22 suite au
! deplacement de l'appel A q_wavebreaking3d dans internal_mode.F90 le !01-05-22
      call obc_int_mpi(2,12)

      do k=1,kmax

! Flux X pour U:
      do j=2,jmax-1 ; do i=1,imax

        xy_t(i,j,id_uflx)=                         &
        xy_t(i,j,id_breaker)*                      &
                         0.1  & ! diffusion dans les 2 directions
                *mask_t(i,j,kmax)                  & ! slip condition on breaking
              *min(dz_u(i+1,j,k,1) ,dz_u(i,j,k,1)) &
                *(vel_u(i+1,j,k,2)-vel_u(i,j,k,2))*dy_t(i,j)*invdx_t(i,j)
      enddo        ; enddo
      checkxyt(2)=xy_t(imax/2,jmax/2,id_uflx)

! Flux Y pour U:
      do j=2,jmax   ; do i=2,imax
        xy_f(i,j,id_ufly)=                         &
        xy_f(i,j,id_breaker)*                      &
                         0.05  & ! diffusion dans les 2 directions + 1/2 somme des 2 gradients
                *mask_f(i,j,kmax)                    & ! slip condition on breaking
              *min(dz_u(i,j,k,1) ,dz_u(i,j-1,k,1))*( & !ooo>
                 (vel_u(i,j,k,2)-vel_u(i  ,j-1,k,2))*dx_f(i,j)*invdy_f(i,j) &
                 +vel_v(i,j,k,2)-vel_v(i-1,j  ,k,2)  &
                                                   )   !ooo>
      enddo        ; enddo
      checkxyf(1)=xy_f(imax/2,jmax/2,id_ufly)

! Flux Y pour V:
      do j=1,jmax ; do i=2,imax-1
        xy_t(i,j,id_vfly)=                         &
        xy_t(i,j,id_breaker)*                      &
                         0.1  & ! diffusion dans les 2 directions
                *mask_t(i,j,kmax)                    & ! slip condition on breaking
              *min(dz_v(i,j+1,k,1) ,dz_v(i,j,k,1)) &
                *(vel_v(i,j+1,k,2)-vel_v(i,j,k,2))*dx_t(i,j)*invdy_t(i,j)
      enddo       ; enddo
      checkxyt(5)=xy_t(imax/2,jmax/2,id_vfly)

! Flux X pour V:
      do j=2,jmax ; do i=2,imax  
        xy_f(i,j,id_vflx)=                         &
        xy_f(i,j,id_breaker)*                      &
                         0.05  & ! diffusion dans les 2 directions + 1/2 somme des 2 gradients
                *mask_f(i,j,kmax)                    & ! slip condition on breaking
              *min(dz_v(i,j,k,1) ,dz_v(i-1,j,k,1))*( & !ooo>
                 (vel_v(i,j,k,2)-vel_v(i-1,j  ,k,2))*dy_f(i,j)*invdx_f(i,j) &
                 +vel_u(i,j,k,2)-vel_u(i  ,j-1,k,2)  &
                                                   )   !ooo>
      enddo       ; enddo
      checkxyf(2)=xy_f(imax/2,jmax/2,id_vflx)

! Composante U
      do j=2,jmax-1 ; do i=2,imax
           vel_u(i,j,k,2)=mask_u(i,j,kmax)*( & !pmx>
           vel_u(i,j,k,2)+(                           &
                            xy_t(i  ,j  ,id_uflx)     &
                           -xy_t(i-1,j  ,id_uflx)     &
                           +xy_f(i  ,j+1,id_ufly)     &
                           -xy_f(i  ,j  ,id_ufly)     &

                          )/( & !ovo>
                               dxdy_u(i,j)*           &
                             max(dz_u(i,j,k,1),0.001) &
                            ) & !ovo>
                                             ) !pmx>

      enddo         ; enddo

! Composante V
      do j=2,jmax ; do i=2,imax-1
           vel_v(i,j,k,2)=mask_v(i,j,kmax)*( & !pmx>
           vel_v(i,j,k,2)+(                           &
                            xy_t(i  ,j  ,id_vfly)      &
                           -xy_t(i  ,j-1,id_vfly)      &
                           +xy_f(i+1,j  ,id_vflx)      &
                           -xy_f(i  ,j  ,id_vflx)      &

                           )/( & !ovo>
                               dxdy_v(i,j)*           &
                             max(dz_v(i,j,k,1),0.001) &
                             ) & !ovo>
                                              ) !pmx>
      enddo       ; enddo

      enddo ! k loop

      enddo                          !PMXPMX>

      RETURN

      end subroutine q_wavebreaking3d

!..............................

      subroutine q_trigger_wavebreaking
      implicit none

!     x0_r4=dti_fw/600.  ! Pour la methode "ssh-sshlwf"
!     do j=1,jmax ; do i=1,imax
!      sshlwf_w(i,j,1)=(1.-x0_r4)*sshlwf_w(i,j,1)+x0_r4*ssh_int_w(i,j,1)
!     enddo       ; enddo

!     x0_r4=dti_fw/100. ! pour analyser seulement la periode recente
!     do j=1,jmax ; do i=1,imax
!      sshlwf_w(i,j,1)=(1.-x0_r4)*sshlwf_w(i,j,1)+x0_r4*ssh_int_w(i,j,1)
!     enddo       ; enddo

! Min Max
       x0_r4=1.-dti_fw/1000. ! pour analyser seulement la periode recente
       do j=1,jmax     ;  do i=1,imax
        sshmax_w(i,j)=max(sshmax_w(i,j)*x0_r4,ssh_int_w(i,j,1))  
        sshmin_w(i,j)=min(sshmin_w(i,j)*x0_r4           &
! Ces lignes ont pour but que le minimum puisse etre au dessus du niveau z=0m (en cas de surcote par ex)
!                         +(1.-x0_r4)*max(0.,-h_w(i,j)) &
!                         +(1.-x0_r4)*sshmax_w(i,j)     &
                          +(1.-x0_r4)*ssh_int_w(i,j,1)  &
                          ,ssh_int_w(i,j,1))  
       enddo           ; enddo

      x0_r4=1.-dti_fw/100. ! pour analyser seulement la periode recente
      do j=1,jmax     ;  do i=1,imax
       sshmin_tmp_w(i,j)=min(sshmin_tmp_w(i,j)*x0_r4 &
! Ces lignes ont pour but que le minimum puisse etre au dessus du niveau z=0m (en cas de surcote par ex)
!                                         +(1.-x0_r4)*sshlwf_w(i,j,1)  &
!                                         +(1.-x0_r4)*max(0.,-h_w(i,j)) &
!                                         +(1.-x0_r4)*sshmax_w(i,j) &
                                          +(1.-x0_r4)*ssh_int_w(i,j,2) &
                            ,ssh_int_w(i,j,2))  
      enddo           ; enddo
      if(iteration3d==0) then
       do j=1,jmax     ;  do i=1,imax
        sshmin_tmp_w(i,j)=max(0.,-h_w(i,j))
        sshmax_w(i,j)=sshmin_tmp_w(i,j)
        sshmin_w(i,j)=sshmin_tmp_w(i,j)
       enddo           ; enddo
      endif
       

      x1_r4=        1./(1.-brk_crit_r)
      x2_r4=brk_crit_r/(1.-brk_crit_r)

!     x0_r4=1.-dti_fw/(0.5*periodpeak)
!     x0_r4=   dti_fw/10. ! decroissance 10s
      x0_r4=   dti_fw/20. ! decroissance 20s

      do j=1,jmax ; do i=1,imax

      breaker2d_t(i,j)=   &
                       min(1.5, & !MIN>
         
            max( & !MAX> 
       
         0., &
!             breaker2d_t(i,j)*(1.-x0_r4) & ! methode exponentielle
              breaker2d_t(i,j)    -x0_r4  & ! methode lineaire

!!!!    ,   x1_r4*(ssh_int_w(i,j,2)/max(h_w(i,j),0.03)*inv_brkh)**2-x2_r4         &
        ,   x1_r4*(0.5*(ssh_int_w(i,j,2)-sshmin_tmp_w(i,j))/max(h_w(i,j),0.03)*inv_brkh)**2-x2_r4         &

!       ,   x1_r4*(                                            &  
!             ( 0.5*( abs(ssh_int_w(i+1,j,2)-ssh_int_w(i,j,2)) &
!                    +abs(ssh_int_w(i-1,j,2)-ssh_int_w(i,j,2)))*invdx_t(i,j) )**2  &
!            +( 0.5*( abs(ssh_int_w(i,j+1,2)-ssh_int_w(i,j,2)) &
!                    +abs(ssh_int_w(i,j-1,2)-ssh_int_w(i,j,2)))*invdy_t(i,j) )**2  &
!                 )*inv_brkslope2-x2_r4                        & 

!       ,   x1_r4*( ( ( &
!                 0.5*abs( & !ooo>
!            ssh_int_w(i,j  ,2)-ssh_int_w(i,j  ,1)-ssh_int_w(i,j  ,0)+ssh_int_w(i,j  ,-1)  &
!        +invgrav*(q_t(i,j,k,2)      -q_t(i,j,k,1)      -q_t(i,j,k,0)      +q_t(i,j,k,-1)) & 
!                       ) & !ooo>
!                   )*inv_dti_fw_p2 )**2     &
!                 )*inv_brkvar-x2_r4         & 

               ) & !MAX>
                            )    !MIN>

#ifdef bidon
       if(i==imax/2.and.j+par%tjmax(1)==jglb-60) then
!       write(10+par%rank,'(f12.3,1x,5(1x,e14.7),a)') &
        if(iteration3d>0)         &
        write(10+par%rank,*)      &
         real(elapsedtime_now)    &
!       ,real(breaker2d_t(i,j))   &
        ,real(ssh_int_w(i,j,2))   &
        ,real(sshmin_tmp_w(i,j))  &
        ,real(sshmin_w(i,j))      &
        ,real(sshmax_w(i,j))       
!       ,real(sshlwf_w(i,j,1))  
!       ,real(x1_r4*min((ssh_int_w(i,j,2)/max(h_w(i,j),0.03)*inv_brkh)**2,1.)-x2_r4)        &
!       ,real(x1_r4*min((                                            &  
!             ( 0.5*( abs(ssh_int_w(i+1,j,2)-ssh_int_w(i,j,2)) &
!                    +abs(ssh_int_w(i-1,j,2)-ssh_int_w(i,j,2)))*invdx_t(i,j) )**2  &
!            +( 0.5*( abs(ssh_int_w(i,j+1,2)-ssh_int_w(i,j,2)) &
!                    +abs(ssh_int_w(i,j-1,2)-ssh_int_w(i,j,2)))*invdy_t(i,j) )**2  &
!                 )*inv_brkslope2,1.)-x2_r4)                         
       endif
#endif


      enddo       ; enddo

      end subroutine q_trigger_wavebreaking

!............................
      end module module_q
