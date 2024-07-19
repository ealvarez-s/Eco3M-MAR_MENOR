&notebook_vertcoord

! Vertical distribution of the s levels
! https://docs.google.com/document/d/1lCiiOzhpXSxY0fMnc9iKGCnMSZWYbfQBy1ypdYAXQG0/edit

igesig=0       ! 0=sigma coordinate  1=generalized coordinate      
isigfile=0     ! sigma coodinate is: 0=computed or 1=read from "sigma.in" input file

!..............................................
! Relevant if igesig=1 (generalized coordinate)
hgesig=100.    ! H0. Regular sigma coordinate if h<H0. "Attracted" gen. coord. if h>H0
!..............................................


!..............................................
! Relevant if igesig=0 (sigma coordinate)
dz_vertical_incr_fact=1.08  ! Depth Resolution Increasing Factor
!..............................................

!.........  s-z hybrid grid section  .......................
flag_merged_levels=1 ! merged levels if flag_merged_levels=1
vqs_cst1=100.        ! (m) Bathymetry envelope when true bathymetry=0
vqs_cst2=50.         ! (m) True bathymetry below which the grid is 100% sigma
vqs_cst3=0.5         ! Curvature factor of the envelope bathymetry (avoid levels and bathymetry having opposite sign slopes)

! Relevant if flag_merged_levels=1:
nbvstepmin=5         ! Minimum number of levels

! Building the hybrid grid (VQS) from a file of "envelope bathymetry" characteristics:
! https://docs.google.com/document/d/1lCiiOzhpXSxY0fMnc9iKGCnMSZWYbfQBy1ypdYAXQG0/edit
!vqs_file='none'
 vqs_file='../../../BGC_MAR_MENOR/BATHYMASK/vqs_file'
!.........  s-z hybrid grid section  .......................

dzsurfmin=2. ! if>0 dzsurfmin(m) preserves surface resolution: surface delta z = MIN (surface delta z, dzsurfmin )

flag_z2dv_outputs=1 ! If=1 produce ascii files (named z2dv_OE_rank, z2dv_SN_rank) to visualize the position of grid levels in vertical planes


! Other parameters
vststep=0
isigfile=0     ! sigma coodinate is: 0=computed or 1=read from "sigma.in" input file

!Obsolete parameters
! Sigma / Step grid:
!ihybsig=0      ! 0 unused. 1=compute the Hybrid s-step grid. 2=read it in a file  
!nbvstepmin=5   ! Case IHYBSIG=1. Minimum number of levels (in the shallow areas)
!sigstepgridfile="../../../MAR_MENOR/BATHYMASK/sigstepgridfile.nc"

!nhybsig=0    ! Case IHYBSIG=1. Number of loops of the iterative process              
!hstepmin=3.    ! Case IHYBSIG=1. Define "the shallow area" with a H criteria  
!hstepmax=150.  ! Case IHYBSIG=1. H1. Use the kmax "declared" levels when h>H1   
!verticalgridfile="../../../MAR_MENOR/z_level.txt"
!ale_selected=0 ! "Arbitray Lagrangian Eulerian" grid (1=yes, 0=no)  

/
