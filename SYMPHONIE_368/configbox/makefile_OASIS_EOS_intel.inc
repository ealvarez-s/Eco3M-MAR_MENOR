#__________________________________________________________
#__________________________________________________________
# Options de compilations: parametrer F77FLAGS
#__________________________________________________________
# Intel version
# EOS version : 
#######################################################################################################################################################################
# Currently Loaded Modulefiles:
#  1) intel/16.1.1             2) intelmpi/4.1.3.049      3) hdf5/1.8.14_intelmpi     4) intel/mkl.11.0.5         5) netcdf/4.3.2             6) pnetcdf/1.5.0_intelmpi
#######################################################################################################################################################################

#module load intel/16.1.1 intelmpi/4.1.3.049 hdf5/1.8.14_intelmpi intel/mkl.11.0.5 netcdf/4.3.2 pnetcdf/1.5.0_intelmpi
# Compillateur 
F90=mpiifort

# Config
CONF=GDL500MLA
COMP=INTEL
COMPMOD=-module


# Optimisatinos, debug, ... flags
# DO NOT USE -O3  optimization
F77FLAGS=  -O2 -traceback 

# OASIS
OASISDIR=/users/p1247/seyfried/oasis3-mct/oasis3-mct_comp
OASISLIB=$(OASISDIR)/lib/libpsmile.MPI1.a $(OASISDIR)/lib/libmct.a $(OASISDIR)/lib/libmpeu.a $(OASISDIR)/lib/libscrip.a

LIBBUILD        = $(OASISDIR)/build/lib
INCPSMILE       = -I$(LIBBUILD)/psmile.MPI1 -I$(LIBBUILD)/scrip

# parallel netcdf and netcdf librairie
LIB=-L/usr/local/netcdf.4.3.2/lib -lnetcdff -Wl,-rpath,/usr/local/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64 -Wl,-rpath,/usr/local/hdf5-1.8.13/lib -lnetcdf -lnetcdf -lpnetcdf /usr/local/parallel-netcdf/1.5.0_intelmpi/lib/libpnetcdf.a $(OASISLIB)

# Include repertoire pour netcdf.h, ...
NETINC=-I/usr/local/netcdf.4.3.2/include -I/usr/local/parallel-netcdf/1.5.0_intelmpi/include $(INCPSMILE)

#__________________________________________________________
#__________________________________________________________
