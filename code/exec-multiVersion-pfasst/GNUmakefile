# -*- Mode: Makefile -*- 
LIBPFASST = /home/groups/jsuckale/liuwj/PFASST_icesheet/
USE_PFASST = TRUE
MPI = TRUE
HDF5 = TRUE
USE_PETSC=FALSE
FC_COMPILE = gfortran -O3 #mpifort #
LIBGFORTRAN = /share/software/user/open/gcc/12.2.0/lib64 -lmpi -lgfortran#/opt/homebrew/lib/gcc/13
OBJ_DIR = ../o/2d.Linux.64.mpicxx.mpifort.DEBUG.MPI# LIBOMP = /opt/homebrew/opt/libomp
# LIBLLVM = /opt/homebrew/opt/llvm

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := driver
src_dirs = fsrc
src_dirs += csrc

#all:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
#driver:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
all:  libbike driver $(if $(USE_PFASST),  libbike driver pfasst)
driver:  libbike

.NOTPARALLEL: all

include ../mk/Make.defs
OBJ_DIR = o/2d.Linux.64.mpicxx.gfortran.DEBUG.MPI
-include $(OBJ_DIR)/*.d
include .depend
FFLAGS  += -I$(LIBPFASST)/include
FFLAGS  += -I$(LIBPFASST)/build 
FFLAGS  += -I/share/software/user/open/openmpi/4.1.2/lib
FFLAGS  += -MD -MP -cpp
LDFLAGS += -lstdc++ -lquadmath -lmpi_mpifh #
# LDFLAGS += -L/opt/homebrew/Cellar/gcc/13.1.0/lib/gcc/13/gcc/aarch64-apple-darwin22/13 -lgcc
LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst -lmpi
LDFLAGS += "-Wl,--copy-dt-needed-entries"
LDFLAGS += -L/share/software/user/open/openmpi/4.1.2/lib -lmpi
LDFLAGS += -L/home/groups/jsuckale/liuwj/Toolkits/hdf5-1.12.2/lib -lhdf5 -lz
FFLAGS += "-Wl,--fno-underscore"
FFLAGS += "-Wl,-no_compact_unwind"
# LDFLAGS += "-Wl,--fno-underscore"
# LDFLAGS += "-Wl,-no_compact_unwind"
# LDFLAGS += -L$(LIBLLVM)/lib
# CPPFLAGS+= -I$(LIBLLVM)/include -fopenmp
# LDFLAGS += -L$(LIBOMP)/lib -fopenmp
# CPPFLAGS+= -I$(LIBOMP)/include

FFLAGS  += -ffree-line-length-none
include .depend
include $(LIBPFASST)/.depend
# LIBFLAGS  +=  -L$(LIBBISICLE)/d/2d.Darwin.64.mpicxx.mpifort.DEBUG.MPI/
# include d/2d.Darwin.64.mpicxx.mpifort.DEBUG.MPI/bisicles_vector.d
# CFLAGS += -Wall -g



tttt: tttt.txt
	echo "a" > tttt.txt

reallyclean: realclean libclean

libbike: 
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

#LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)
LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS) -L$(LIBPFASST)/lib -lpfasst -L$(LIBGFORTRAN) -lgfortran
LIBFLAGS += -I$(LIBPFASST)/build -I$(LIBPFASST)/include -I/include
LIBFLAGS += -lstdc++ $(LIBPFASST)/lib/libpfasst.a



driver:  libbike ../lib/libBisicles$(config).a 

alltest: libbike ../lib/libBisicles$(config).a 

#include Makefile.fortran
pfasst: 
	mkdir -p include
	touch phony.mod
	mv *.mod include/


.PHONY: all lib libbike libclean alltest
