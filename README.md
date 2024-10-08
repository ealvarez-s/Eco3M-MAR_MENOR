In Olympe: source ∼/.bashrc\
In IEO-cluster: follow instructions in *HT_compile_Eco3M_IEOcluster.pdf*

## CODE   ---------------------------------------------
Code in: SYMPHONIE_368/
 - original in SOURCES/
 - modifications in UDIR/BGC_MAR_MENOR2/

Depending on compiler, create:
  SYMPHONIE_368/CDIR_IFORT or SYMPHONIE_368/CDIR_GFORTRAN

Compile in SYMPHONIE_368/UDIR/BGC_MAR_MENOR2/ with *make*

Run in: SYMPHONIE_368/RDIR/BGC_MAR_MENOR2/

Directories not included, create them:\
SYMPHONIE_368/RDIR/BGC_MAR_MENOR2/restart_input/\
SYMPHONIE_368/RDIR/BGC_MAR_MENOR2/restart_outbis/\
SYMPHONIE_368/RDIR/BGC_MAR_MENOR2/restart_output/\
SYMPHONIE_368/RDIR/BGC_MAR_MENOR2/tmp/

## REGION ---------------------------------------------
Setup in: BGC_MAR_MENOR

Offline files with the circulation in: BGC_MAR_MENOR/OFFLINE/

Boundary conditions and initialization in: GLOBMED2/

Atmospheric forcing in: meteo/
