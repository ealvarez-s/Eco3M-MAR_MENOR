&notebook_airseaflux

! ECMWF case
! https://docs.google.com/document/d/1rCOZKUdX_w_FgPzEeSS7g9ps3qKTRc00_CTjie_DK5A/edit

iairsea=1                        ! 1 = use the meteo fluxes , 0 otherwise
flag_meteodata='ecmwf'           ! Meteorological model key (ecmwf glorys)
airseaoption=1                   ! 1 = online interpolation  
irelaxsst=0                      ! 1 = relaxed sst (toward external ogcm sst) 0 otherwise (default value)
flag_meteo_average=0             ! 1 = non cumulated variables were averaged over the sampling period, 0 otherwise.
flag_ssr24avr=0                  ! 1 = 24h averaged solar flux (warning: average centre 12h before the current time)

! Land sea mask:
flag_lsm=0                  ! 1 = land-sea mask is not needed, 0 otherwise (run stops if lsm not found)
meteolandconvention=1       ! Value of the mask in land nodes (1 if ecmwf convention, 0 if surfex convention)
flag_meteo_land_plug=1      ! =0 land values are used, =1 land values replaced by nearest sea values (all fields except the wind).
flag_meteo_land_plug_wind=0 ! =0 land values are used, =1 land values replaced by nearest sea values (only the wind)
meteo_grid_file='none'      ! if LSM is not in variable files, give a grid file with LSM

! How I decumulate variables cumulated over time (SSR etc...)
! t0 convention for cumulated var:
 meteo_t0='file'        ! t0=0 first field of the 'file' 
!meteo_t0='time'        ! t0=0 when 'time'=0
!meteo_t0='user'        ! t0=0 when modulo(time,meteo_cumul_modulo)=0
!meteo_cumul_modulo=24  ! (if meteo_t0='user' only) time (in hours) on which the variables are cumulated
!meteo_t0='auto'        ! automatic detection of the solar flux reset
!meteo_t0='unli'        ! unlimited (the accumulation is never reset)

! Remove noise in ecmwf files (if any)
flag_p0m_filter=0                ! If 1 (0 otherwise) apply the 2d noise killer to ecmwf atmospheric pressure field

! Relative wind:
 relativewind=0.   ! Fraction of the surface current substracted to the wind speed
!relativewind=0.77 ! this case (although not recommended) correponds to (1-Sw)=(1-0.23) from Renault et al 2016, JPO, 46, 1685-1704

! Meteo grid overflow
flag_gridoverflow=0              ! 0=reject 1=accept meteo grid overflow (ocean grid larger than meteo grid)

! Simple Atmospheric Boundary Layer
flag_abl=0    ! 1= retro-active atmospheric boundary layer (temperature and humidity)     SABL1
flag_abl2=0   ! 1= retro-active atmospheric boundary layer (temperature, humidity, wind)  SABL2
! note: if flag_abl2=1 requires flag_abl=1 and relativewind=1.

! Bulk scheme:
bulk_scheme=3       ! 1=CORE, 2=MOON ,3=COARE
flag_wstressbulk=1  ! If flag_wstressbulk=0 wtress is given by the meteo model (instead bulk formulae)

! Files lists:
texte90='../../../MAR_MENOR/LIST/'  ! Path access to the following file lists
airseafile(1)='liste_ecmwf_ssr'
airseafile(2)='liste_ecmwf_ir'
airseafile(3)='liste_ecmwf_u10m'
airseafile(4)='liste_ecmwf_v10m'
airseafile(5)='liste_ecmwf_p0m'
airseafile(6)='liste_ecmwf_t2m'
airseafile(7)='liste_ecmwf_dp2m'
airseafile(8)='liste_ecmwf_rain'
!  airseafile(9)= 'liste_ecmwf_ustrs'  ! relevant if flag_wstressbulk=0
!  airseafile(10)='liste_ecmwf_vstrs'  ! relevant if flag_wstressbulk=0
/



ICI un exemple de reglage pour le modele AROME:

! AROME case

iairsea=1                        ! 1 = use the meteo fluxes , 0 otherwise
flag_meteodata='arome'           ! Meteorological model key (ecmwf glorys)
airseaoption=1                   ! 1 = online interpolation  
irelaxsst=0                      ! 1 = relaxed sst (toward external ogcm sst) 0 otherwise (default value)
flag_meteo_average=0             ! 1 = non cumulated variables were averaged over the sampling period, 0 otherwise.
flag_ssr24avr=0                  ! 1 = 24h averaged solar flux (warning: average centre 12h before the current time)

! Land sea mask:
flag_lsm=1                       ! 1 = land-sea mask is not needed, 0 otherwise (run stops if lsm not found)
meteolandconvention=1            ! Value of the mask in land nodes (1 if ecmwf convention, 0 if surfex convention)

meteo_t0='file' !='time'         ! t0 convention for cumulated var. t0=0 first field of the 'file' or t0=0 when 'time'=0
flag_p0m_filter=0                ! If 1 (0 otherwise) apply the 2d noise killer to ecmwf atmospheric pressure field
relativewind=0. !0.5             ! Fraction of the surface current substracted to the wind speed

! Simple Atmospheric Boundary Layer
flag_abl=0    ! 1= retro-active atmospheric boundary layer (temperature and humidity)     SABL1
flag_abl2=0   ! 1= retro-active atmospheric boundary layer (temperature, humidity, wind)  SABL2
! note: if flag_abl2=1 requires flag_abl=1 and relativewind=1.

! Bulk scheme:
bulk_scheme=3  ! 1=CORE, 2=MOON ,3=COARE

! Files lists:
texte90='../../../MAR_MENOR/LIST/'  ! Path access to the following file lists
airseafile(1)='liste_arome_ssr'
airseafile(2)='liste_arome_ir'
airseafile(3)='liste_arome_u10m'
airseafile(4)='liste_arome_v10m'
airseafile(5)='liste_arome_p0m'
airseafile(6)='liste_arome_t2m'
airseafile(7)='liste_arome_dp2m'
airseafile(8)='liste_arome_rain'
