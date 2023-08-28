# -*- Mode: Makefile -*- 
LIBPFASST = /Users/emmaliu/packages/PFASST_icesheet
USE_PFASST = TRUE
FC_COMPILE = mpifort #gfortran -O3
LIBGFORTRAN = /opt/homebrew/lib/gcc/13
# LIBOMP = /opt/homebrew/opt/libomp
# LIBLLVM = /opt/homebrew/opt/llvm

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := driver
src_dirs = fsrc
src_dirs += csrc

#all:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
#driver:  libBisicles2d.Linux.64.g++.gfortran.DEBUG.a libbike 
all: clean libbike driver $(if $(USE_PFASST), clean libbike driver pfasst)
driver: clean libbike

.NOTPARALLEL: all

include ../mk/Make.defs
FFLAGS  += -I$(LIBPFASST)/include
FFLAGS  += -I$(LIBPFASST)/build 
FFLAGS  += -I/include
LDFLAGS += -lstdc++ -lmpi_mpifh -lquadmath 
LDFLAGS += -L/opt/homebrew/Cellar/gcc/13.1.0/lib/gcc/13/gcc/aarch64-apple-darwin22/13 -lgcc
LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst
LDFLAGS += "-Wl,-no_compact_unwind"
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



driver: clean libbike ../lib/libBisicles$(config).a 

alltest: libbike ../lib/libBisicles$(config).a 

#include Makefile.fortran
pfasst: 
	mkdir -p include
	mv *.mod include/


.PHONY: all lib libbike libclean alltest