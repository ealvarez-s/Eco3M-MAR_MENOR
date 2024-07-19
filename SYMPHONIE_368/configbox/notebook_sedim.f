&notebook_sedim

! this namelist is read by the s26 routine module_parameter_sedim.F90

!--- parameter for the module sedim extract from MARS-SEDIM
!l_sedim_s=.true.                     ! switch .true. = ON, .false. = OFF
 l_sedim_s=.false.                     ! switch .true. = ON, .false. = OFF
!
!nb_var=3                           ! number of sediment variables
 nb_var=0                           ! number of sediment variables
                                   ! something as to be done about the obc and river discharge

!ksdmax=110                         ! number of maximum sediment layers stored into the sediment 
 ksdmax=0                           ! number of maximum sediment layers stored into the sediment 

file_parasedim='parasedim.txt'     ! file containing the parameter for the sedimental module

file_varsedim='parasubs.txt'       ! file containing the derscription of nature of the nb_var sediment variables
!...............................................................................

/
