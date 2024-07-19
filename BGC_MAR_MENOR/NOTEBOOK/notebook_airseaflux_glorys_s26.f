&notebook_airseaflux

! GLORYS case

iairsea=0                       ! 1 = use the meteo fluxes , 0 otherwise
flag_meteodata='glorys'         ! Meteorological model key (ecmwf glorys)
airseaoption=1                  ! 1 = online interpolation  
irelaxsst=0                     ! 1 = relaxed sst (toward external ogcm sst) 0 otherwise (default value)
flag_meteo_average=1            ! 1 = non cumulated variables were averaged over the sampling period, 0 otherwise.
flag_ssr24avr=0                 ! 1 = 24h averaged solar flux (warning: average centre 12h before the current time)
meteo_t0='file' !='time'        ! t0 convention for cumulated var. t0=0 first field of the 'file' or t0=0 when 'time'=0
texte90='../../../BGC_MAR_MENOR/LIST/' ! Path access to the following file lists

! Files lists:
airseafile(1)='list_glorys_ustrs'
airseafile(2)='list_glorys_vstrs'
airseafile(3)='list_glorys_ssr'
airseafile(4)='list_glorys_slhf'
airseafile(5)='list_glorys_netir'
airseafile(6)='list_glorys_sshf'
airseafile(7)='list_glorys_rain'
/
