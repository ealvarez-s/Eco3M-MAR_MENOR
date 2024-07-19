&notebook_tracer1

 imodeltrc=0  ! Passive tracers model: 1=on  0=off
 vbmax=1      ! Number of independent 3D tracers
 ksomax=1     ! Ponctual sources number

/

&notebook_tracer2

! Settling velocity for each tracer in m/s (sign convention:  <0):
wsed(1,1)=-0.   
!wsed(1,2)=-0.003
!wsed(1,3)=-0.003

! Is the tracer a radionucleide or not (default =0)
Radionucleide(1)=1

! Half-life radioactive decay for each tracer if they are radioactive in days (default is 1000000.):
TauR(1)=1000000.   
!TauR(2)=8.0 ! for exemple Iode_131

! Concentrations at river mouth 
! Must be specified for all rivers (same order as in notebook_river) and for all tracers:
! river_bio(tracer counter,river counter,1)=100.
! river_bio(1,1,1)=100. ! variable n°1 river n°1
! river_bio(2,1,1)=100. ! variable n°2 river n°1

! river_bio(1,2,1)=100. ! variable n°1 river n°2
! river_bio(2,2,1)=100. ! variable n°2 river n°2

! river_bio(1,3,1)=100. ! variable n°1 river n°3
! river_bio(2,3,1)=100. ! variable n°2 river n°3

! river_bio(1,4,1)=100. ! variable n°1 river n°4
! river_bio(2,4,1)=100. ! variable n°2 river n°4

! Additional source points:
! Source n°1
! Horizontal coordinates:
  socard( 0,1)=0     ! if a=0 then b=i c=j / if a=1 then b=lat(°) c=lon(°)
  socard( 1,1)=709   ! b
  socard( 2,1)=334   ! c
! Depth:
  socard(-1,1)=1.    ! (a,b,c) if a=0 then b=kmin c=kmax / if a=1 then b=zmin c=zmax
  socard( 3,1)=-3.   ! b
  socard( 4,1)=-0.   ! c
! tracers concerned (min & max number)
  socard( 5,1)=1     ! tracers min number
  socard( 6,1)=1     ! tracers max number
  socard( 7,1)=100.0 ! Production term (per second)
! Dates:
  sodate(1:6,1,1)= 2013 , 01 , 01 , 01  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)
  sodate(1:6,2,1)= 2013 , 01 , 01 , 02  , 00  , 00  ! End   time (yyyy mm dd hh mm ss)
  
! Source n°2
! Horizontal coordinates:
! socard( 0,2)=0.   ! if a=0 then b=i c=j / if a=1 then b=lat(°) c=lon(°)
! socard( 1,2)=110. ! b
! socard( 2,2)=21.  ! c
! Depth:
! socard(-1,2)=1.   ! (a,b,c) if a=0 then b=kmin c=kmax / if a=1 then b=zmin c=zmax
! socard( 3,2)=-3.  ! b
! socard( 4,2)=-0.  ! c
! tracers concerned (min & max number)
! socard( 5,2)=1     ! tracers min number
! socard( 6,2)=2     ! tracers max number
! socard( 7,2)=0.001 ! Production term (per second)
! Dates:
! sodate(1:6,1,2)= 2009 , 03 , 01 , 00  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)
! sodate(1:6,2,2)= 2009 , 04 , 01 , 00  , 00  , 00  ! End   time (yyyy mm dd hh mm ss)

/
