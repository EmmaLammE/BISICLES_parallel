# -*- Mode: Makefile -*- 

#LIBGFORTRAN = /usr/local/lib/gcc/11
LIBPFASST = /home/groups/jsuckale/liuwj/PFASST_icesheet/
USE_PFASST = TRUE
MPI = TRUE
HDF5 = TRUE
USE_PETSC=FALSE
FC_COMPILE = gfortran -O3 #mpifort #
LIBGFORTRAN = /share/software/user/open/gcc/12.2.0/lib64 -lmpi -lgfortran#/opt/homebrew/lib/gcc/13
OBJ_DIR = ../o/2d.Linux.64.mpicxx.mpifort.DEBUG.MPI# LIBOMP = /opt/homebrew/opt/libomp
OBJ_DIR = o/2d.Linux.64.mpicxx.gfortran.DEBUG.MPI

-include $(OBJ_DIR)/*.d
include .depend
FFLAGS  += -I$(LIBPFASST)/include
FFLAGS  += -I$(LIBPFASST)/build 
FFLAGS  += -I/share/software/user/open/openmpi/4.1.2/lib
FFLAGS  += -MD -MP -cpp
LDFLAGS += -lstdc++ -lquadmath -lmpi_mpifh #
FFLAGS  += -ffree-line-length-none
# LDFLAGS += -L/opt/homebrew/Cellar/gcc/13.1.0/lib/gcc/13/gcc/aarch64-apple-darwin22/13 -lgcc
# LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst -lmpi
# LDFLAGS += "-Wl,--copy-dt-needed-entries"
# LDFLAGS += -L/share/software/user/open/openmpi/4.1.2/lib -lmpi
# LDFLAGS += -L/home/groups/jsuckale/liuwj/Toolkits/hdf5-1.12.2/lib -lhdf5 -lz
# FFLAGS += "-Wl,--fno-underscore"
# FFLAGS += "-Wl,-no_compact_unwind"

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := driver

#all:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
#driver:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
all: clean libbike driver
driver: clean libbike

.NOTPARALLEL: all

include ../mk/Make.defs


tttt: tttt.txt
	echo "a" > tttt.txt

reallyclean: realclean libclean

libbike: 
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS) -L$(LIBGFORTRAN) -lgfortran

driver: clean libbike ../lib/libBisicles$(config).a 

alltest: libbike ../lib/libBisicles$(config).a 


# makes all of the fortran headers needed by the build (in case make fails to do it automagically)
fortranheaders: ../src/amrIceF_F.H ../src/ExtrapBCF_F.H ../src/AVCAMRPoissonOpF_F.H ../src/FIBCF_F.H ../src/BasalFrictionRelationF_F.H ../src/FourthOrderAverageF_F.H ../src/BisiclesF_F.H ../src/IceThermodynamicsF_F.H ../src/CalvingF_F.H ../src/IntInterpF_F.H ../src/ConstitutiveRelationF_F.H ../src/L1L2ConstitutiveRelationF_F.H ../src/ControlF_F.H ../src/NyeCrevasseF_F.H ../src/CrevasseF_F.H ../src/SigmaCSF_F.H ../src/DerivativesF_F.H ../util/DivergenceF_F.H


.PHONY: all lib libbike libclean alltest


