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
CONF=REGION
COMP=IFORT
COMPMOD=-module  

# Optimisatinos, debug, ... flags
# DO NOT USE -O3  optimization
#F77FLAGS= -g -traceback -O0 -ftrapuv -fpe0 -check bounds
 F77FLAGS= -O2 -traceback -fp-model precise
#F77FLAGS= -traceback -fp-model precise -check bounds


# OASIS
#OASISDIR=/users/p09115/cestour/softs/oasis3-mct/compile_intelmpi_olympe
OASISDIR=/users/p1325/marsale/softs/oasis3-mct/compile_intelmpi_olympe
OASISLIB=$(OASISDIR)/lib/libpsmile.MPI1.a $(OASISDIR)/lib/libmct.a $(OASISDIR)/lib/libmpeu.a $(OASISDIR)/lib/libscrip.a

LIBBUILD        = $(OASISDIR)/build/lib
INCPSMILE       = -I$(LIBBUILD)/psmile.MPI1 -I$(LIBBUILD)/scrip

# parallel netcdf and netcdf librairie
NETINC=-I/usr/local/netcdf/4.6.1-intelmpi/include -I/usr/local/pnetcdf/1.9.0/include $(INCPSMILE)
#NETINC=-I/usr/local/netcdf/4.6.1/include -I/usr/local/pnetcdf/1.9.0/include $(INCPSMILE)

LIB=-L/usr/local/netcdf/4.6.1-intelmpi/lib -lnetcdff -Wl,-rpath,/usr/local/intel/2018.2.046/compilers_and_libraries/linux/lib/intel64 -Wl,-rpath,/usr/local/hdf5/1.10.2/intel_mpi/lib -lnetcdf -lnetcdf /usr/local/pnetcdf/1.9.0/lib/libpnetcdf.a $(OASISLIB)
#LIB=-L/usr/local/netcdf/4.6.1/lib -lnetcdff -Wl,-rpath,/usr/local/intel/2018.2.046/compilers_and_libraries/linux/lib/intel64 -Wl,-rpath,/usr/local/hdf5/1.10.2/intel_mpi/lib -lnetcdf -lnetcdf /usr/local/pnetcdf/1.9.0/lib/libpnetcdf.a $(OASISLIB)



###---------------------------------------------------------------------------------------------------------
#LIB=-lnetcdff -lnetcdf /panfs/panasas/cnt0026/lat0088/SHARED/intel/.17.0/intelmpi/5.1.3.258/lib/libpnetcdf.a

# Include repertoire pour netcdf.h, ...
#NETINC=-I/usr/local/netcdf.4.3.2/include -I/panfs/panasas/cnt0026/lat0088/SHARED/intel/.17.0/intelmpi/5.1.3.258/include/

#OLYMPE
#NETINC=-I/tmpdir/nguyen/CHALLENGE/softs/intel/18.2/openmpi/icc/mt/2.0.2.10/include/
#LIB=/tmpdir/nguyen/CHALLENGE/softs/intel/18.2/openmpi/icc/mt/2.0.2.10/lib/libnetcdff.a /tmpdir/nguyen/CHALLENGE/softs/intel/18.2/openmpi/icc/mt/2.0.2.10/lib/libnetcdf.a /tmpdir/nguyen/CHALLENGE/softs/intel/18.2/openmpi/icc/mt/2.0.2.10/lib/libpnetcdf.a /tmpdir/nguyen/CHALLENGE/softs/intel/18.2/openmpi/icc/mt/2.0.2.10/lib/libhdf5_hl.a /tmpdir/nguyen/CHALLENGE/softs/intel/18.2/openmpi/icc/mt/2.0.2.10/lib/libhdf5.a -lz -lcurl
