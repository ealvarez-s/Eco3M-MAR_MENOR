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
#F90=mpif90
F90=mpiifort

#  option specifies where to put .mod files
CONF=BGC_MAR_MENOR2
COMP=IFORT
COMPMOD=-module  

# Optimisatinos, debug, ... flags
# DO NOT USE -O3  optimization
#F77FLAGS= -g -traceback -O0 -ftrapuv -fpe0 -check all 
#F77FLAGS= -g -traceback -O0 -ftrapuv -fpe0 -check bounds
# F77FLAGS= -g -traceback -O0 -check bounds -fp-model precise
 F77FLAGS= -O2 -traceback -fp-model precise

# parallel netcdf and netcdf librairie
NETINC=-I/usr/local/netcdf/4.6.1-intelmpi/include -I/usr/local/pnetcdf/1.9.0/include 

LIB=-L/usr/local/netcdf/4.6.1-intelmpi/lib -lnetcdff -Wl,-rpath,/usr/local/intel/2018.2.046/compilers_and_libraries/linux/lib/intel64 -Wl,-rpath,/usr/local/hdf5/1.10.2/intel_mpi/lib -lnetcdf -lnetcdf /usr/local/pnetcdf/1.9.0/lib/libpnetcdf.a


