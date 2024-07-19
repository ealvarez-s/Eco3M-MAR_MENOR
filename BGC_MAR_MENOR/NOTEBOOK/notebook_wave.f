
! https://docs.google.com/document/d/1KRe-gs3nHyn_N8ka3dtlzjc43aYWG9Fpr2Au58aKoBs/edit

! Taw, Two, Tao: https://docs.google.com/presentation/d/11aOwLq-XSytzzCTpDZkwHg2QHecDKI2jjT-7KltXE2I/edit#slide=id.p

&notebook_wave1
 iwve=0 ! 0: off 
!iwve=1 ! 1: Michaud et al, 2012       ! http://dx.doi.org/10.5194/os-8-657-2012
!iwve=2 ! 2: McWilliams 1999  (Eq35)   ! https://doi.org/10.1175/1520-0485(1999)029%3C2523:TWDOC%3E2.0.CO;2 
!iwve=3 ! 3: As case iwve=2 but with HS and orbital speed (Gentil et al 2022, https://doi.org/10.1016/j.csr.2022.104721)
/

&notebook_wave2

  txt_wavemodeltype='ww3_on_sgrid'    ! ww3 input field at "w" location of S grid (no interpolation)
! txt_wavemodeltype='ww3_native_grid' ! ww3 input field on native structured or unstructured ww3 grid (interpolation requiered)
wave_obc_type=1                       ! wave obc type
ww3_flagcumul=0                       ! wave model variables are 0=instantaneous 1=time-averaged

! flag_wavelength relevant if iwve=2:
flag_wavelength=1 ! =1: use the wavelength directly otherwise (=0) deduce it from the period and the dispersion relation

! v268: a single, multi-variable, list:
texte80(1)='../../../MAR_MENOR/LIST/list_ww3'

!texte80(1)='../../../MAR_MENOR/LIST/list_dir'
!texte80(2)='../../../MAR_MENOR/LIST/list_foc'
!texte80(3)='../../../MAR_MENOR/LIST/list_hs'
!texte80(4)='../../../MAR_MENOR/LIST/list_t'
!texte80(5)='../../../MAR_MENOR/LIST/list_taw'
!texte80(6)='../../../MAR_MENOR/LIST/list_two'
!texte80(7)='../../../MAR_MENOR/LIST/list_uss'
!texte80(8)='../../../MAR_MENOR/LIST/list_hsw'
!texte80(9)='../../../MAR_MENOR/LIST/list_bhd' ! Bernouilli Head
/

wave obc type:

WAVE_OBC_TYPE=1:
The ogcm forcing fields did not take the wave effect into account (ex: MERCATOR,COPERNICUS fields)

WAVE_OBC_TYPE=2:
The ogcm forcing fields took the wave effect into account and provide Lagrangian currents (ex: SYMPHONIE Offline files)

WAVE_OBC_TYPE=3:
The ogcm forcing fields took the wave effect into account and provide eulerian currents (never used so far)

Exemple de script pour faire les listes:
#!/bin/bash
DIR=$1
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_dir
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_foc
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_hs
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_t
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_taw
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_two
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_uss
ls $DIR/*nc > ../../../MAR_MENOR/LIST/list_hsw
