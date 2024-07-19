&notebook_nh

flag_nh3d=0               ! 0=NH off
                          ! 1=NH on 3D barotropic case (T,S constant) - NO TIME SPLITTING
                          ! 2=NH on 3D baroclinic case                - NO TIME SPLITTING
                          ! 3=NH on 3D baroclinic case WITH TIME SPLITTING and hydrostatic external mode

asselin_nh=0.02           ! q(t)=q(t)+0.5*asselin_nh*[q(t+1)-2*q(t)+q(t-1)]
cfl_nh=0.5 ! 0.7          ! dt_nh=dt_nh_therory*sqrt(cfl_nh)
nh_freqrestart=100000     ! restart file every "nh_freqrestart" 2d iterations
wetdry_cstnh=0.01         ! wetting/drying case: h+ssh (m) ubelow which nh pressure = 0

nh2d_graph_period=14200   ! 1 periode=710    ! periodicity (in number of 2d iterations) of output files
nh2d_graph_spinup=0       ! spinup period (in number of 2d iterations) before first output file

nh_frozensigma=0          ! 1= sigma levels are frozen. 0= sigma levels updated every iterations
flag_nhdiag_timeaverage=1 ! 1= compute mean ssh, hs, mean fluxes
wavebreakfactor=424.      ! calibration coef for wave breaking computation
!nhpgf_reduce=1.          ! NH PGF cancelled if nhpgf_reduce=0.
 nhpgf_reduce=0.98        ! 0.98 possibly improves accuracy

! Some parameters of academic test cases:
 ssh_amplitude_=0.001

/

&notebook_nh2

!dtmultiple=-999.         
 dtmultiple=1.  ! secondes  !  time step adjusted in order to have dtmultiple/dt=integer

/
