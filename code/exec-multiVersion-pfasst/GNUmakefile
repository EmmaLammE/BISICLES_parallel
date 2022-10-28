# -*- Mode: Makefile -*- 

####################
# User-set variables
####################
# PFASST include directory
LIBPFASST = /Users/emmaliu/internship2022/installs/PFASST_design
USE_PFASST = TRUE
# fortran compiler 
FC_COMPILE = mpifort #gfortran -O3
FC_SRC_DIR = fsrc

LIBGFORTRAN = /usr/local/lib/gcc/11
# LIBGFORTRAN += /usr/local/lib/gcc/11/gcc/x86_64-apple-darwin21/11 

# include necessary mk in pfasst
#include Makefile.pfasstsrc
# include $(LIBPFASST)/.depend
# include Makefile.fortran


### fortran codes
FC_FILES = \
        $(FC_SRC_DIR)/pfasst_main.f90 \
        $(FC_SRC_DIR)/probin.f90 \
        $(FC_SRC_DIR)/sweeper.f90 \
        $(FC_SRC_DIR)/level.f90  
CPP_FILES = \
        driver.cpp
OBJ  = $(subst fsrc, build,$(FC_FILES:.f90=.o))
# OBJ += build/driver.o #$(build/driver.o,$(CPPSRC:.cpp=.o))


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
LDFLAGS += -lstdc++ -lmpi_mpifh -lgcc -lquadmath
LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst
FFLAGS  += -ffree-line-length-none
# FFLAGS += -fno-leading-underscore
# FFLAGS += -fno-underscoring
# FFLAGS += -fopenmp
include .depend
include $(LIBPFASST)/.depend
# include $(LIBPFASST)/Makefile.rules
EXE = $(ebase)$(config).ex


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
# XTRALDFLAGS = -L$(LIBPFASST)/lib -lpfasst -lstdc++
# LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst 
# LDFLAGS += -lstdc++ $(LIBPFASST)/lib/libpfasst.a
# LIBFLAGS += -I$(LIBPFASST)/.depend
# LIBFLAGS += -I/.depend

driver: clean libbike ../lib/libBisicles$(config).a


alltest: libbike ../lib/libBisicles$(config).a 


#include Makefile.fortran
pfasst: 
# 	mkdir -p build
	mkdir -p include
# 	@echo "Building pfasst related fortran codes ..."
# 	$(FC_COMPILE) $(FFLAGS) -c  $(FC_FILES)
# 	mv $(subst fsrc, .,$(FC_FILES:.f90=.o)) build/
# 	mv $(subst fsrc,.,$(FC_FILES:.f90=.mod)) include/
	mv *.mod include/

# linking: 
# 	# linking ...
# 	OBJECTS=./o/2d.Darwin.64.mpicxx.gfortran.DEBUG.MPI/driver.o 
# 	OBJECTS+=./o/2d.Darwin.64.mpicxx.gfortran.DEBUG.MPI/encap.o
# 	mpifort $(FFLAGS) ./o/2d.Darwin.64.mpicxx.gfortran.DEBUG.MPI/driver.o ./o/2d.Darwin.64.mpicxx.gfortran.DEBUG.MPI/encap.o -o $(EXE) $(LDFLAGS)



.PHONY: all lib libbike libclean alltest


