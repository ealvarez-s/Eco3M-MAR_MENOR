&notebook_vertcoord

! Vertical distribution of the s levels

igesig=1       ! 0=sigma coordinate  1=generalized coordinate      
isigfile=0     ! sigma coodinate is: 0=computed or 1=read from "sigma.in" input file
hgesig=200.    ! H0. Regular sigma coordinate if h<H0. "Attracted" gen. coord.  if h>H0

!.........  s-z hybrid grid section  .......................
flag_merged_levels=0 ! merged levels if flag_merged_levels=1
! Relevant if flag_merged_levels=1:
hstepmax=3000.       ! Value of the bathymetry for which the grid counts kmax levels (for instance hstepmax=3000meters)
                     ! If hstepmax<0 hstepmax will be fixed at the maximum
                     ! bathymetry on the domain (i.e. hstepmax=hmax).
!.........  s-z hybrid grid section  .......................

dzsurfmin=2. ! if>0 dzsurfmin(m) preserves surface resolution: surface delta z = MIN (surface delta z, dzsurfmin )

vststep=0

/
