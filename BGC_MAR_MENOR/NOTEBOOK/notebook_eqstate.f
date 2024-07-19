&notebook_eqstate

! EOS author: 
  eos_author=3 ! McDougall et al. (2003) EOS with updated coefficients by Jackett et al. (2006)
! eos_author=0 ! linear equation of state

! PGF:
! Compression effect in the EOS
  eos_comprs=1  ! Full EOS (that is including compression effect)
! eos_comprs=0  ! potential density only

! PGF:
! IF EOS_COMPRS=0 ONLY define the potential density in the PGF
! Define the level of reference (in meters, convention: positive) of the potential density
  eos_pgfzref=0.    ! Level of reference at z=0m
! eos_pgfzref=1000. ! Level of reference at z=-1000m

! PGF:
! Remove a (x,y) dependant z-weighted average of rhp before computing the PGF
! rhp_zavr_xy=0     ! NO
  rhp_zavr_xy=1     ! YES

! TKE SCHEME:
! Define depth of reference of the potential density used in the turbulence closure scheme
! If eos_tkezref < 0. the reference is the current depth (i.e. computed twice and derived at p constant)
! eos_tkezref=0.    ! Level of reference at z=0m
! eos_tkezref=1000. ! Level of reference at z=1000m
  eos_tkezref=-999. ! Level of reference at current depth (equivalent to MacDougall JPO 1987)

! ACADEMIC CASES:
! IF EOS_AUTHOR=0 ONLY define the linear EOS
  eos_linear=0 ! The coefficients of the linear EOS are deduced from EOS80 complete formula using the (T,S) initial field
! eos_linear=1 ! The coefficients are set to their default value (set_parameters.F90) or by the following lines:
! t0=    11.6352062                                     
! s0=    35.1110878                                         
! alp_t= 0.000184280361                          
! alp_s= 0.000757332833                   
! rho=   1026.74524                 

! STERIC EFFECT https://docs.google.com/document/d/1x9lWebKZzgpWGAlKypjulwtNvN-YiufIKb3EKN__yk8/edit
flag_steric_effect=0 ! Parametrization of the steric effect if flag_steric_effect=1 

! Min Max validity for temperature and salinity
tem_validmin=-999.   ! -3.
tem_validmax=999.    ! 40. 
sal_validmin=-999.   !  0.
sal_validmax=999.    ! 40.

/

Details sur l'EOS et le PGF
https://docs.google.com/document/d/1DPvKOZ3rV-FQfEaYsX-bTDwVmDg3wB4kbPcdS1TJqoc/edit
