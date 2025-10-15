#__________________________________________________________
#__________________________________________________________
# Options de compilations: parametrer F77FLAGS
#__________________________________________________________
# Intel version
# EOS version : 
#######################################################################################################################################################################
# Currently Loaded Modulefiles:
#  1) intel/16.1.1             2) intelmpi/5.1.2.150       3) hdf5/1.8.14_intelmpi     4) intel/mkl.11.0.5         5) netcdf/4.3.2             6) pnetcdf/1.5.0_intelmpi
#######################################################################################################################################################################


# Compillateur, MPI wrapper -------------------------------------------------------------------
F90=mpif90
F90=mpiifort
#F90=mpiifx

#  option specifies where to put .mod files
CONF=BGC_MAR_MENOR2
COMP=IFORT
COMPMOD=-module  
#COMPMOD=-Mmodules

# Optimisatinos, debug, ... flags
# DO NOT USE -O3  optimization
#F77FLAGS= -g -traceback -O0 -ftrapuv -fpe0 -check all 
#F77FLAGS= -g -traceback -O0 -ftrapuv -fpe0 -check bounds
# F77FLAGS= -g -traceback -O0 -check bounds -fp-model precise -diag-disable=10448
F77FLAGS= -g -O2 -traceback -fp-model precise -diag-disable=10448

# F77FLAGS= -O2 -fallow-argument-mismatch -ffree-line-length-none
# F77FLAGS= -O2 -fbacktrace

## parallel netcdf and netcdf libraries

## con conda
NETINC=-I/home/alvarez/miniconda3/envs/belich-legos-r/include
LIB=-L/home/alvarez/miniconda3/envs/belich-legos-r/lib -lnetcdff -lnetcdf -lpnetcdf ###/home/alvarez/miniconda3/envs/belich-all/lib/libpnetcdf.a

## librerias por separado
##NETINC=-I/home/alvarez/local/pnetcdf-1.12.3/include -I/home/alvarez/local/netcdf-c-4.9.2/include -I/home/alvarez/local/netcdf-fortran-4.6.1/include
##LIB=-L/home/alvarez/local/pnetcdf-1.12.3/lib -lpnetcdf -L/home/alvarez/local/netcdf-c-4.9.2/lib -lnetcdf -L/home/alvarez/local/netcdf-fortran-4.6.1/lib -lnetcdff