# -*- Mode: Makefile -*- 

####################
# User-set variables
####################
# PFASST include directory
LIBPFASST = /Users/emmaliu/internship2022/installs/PFASST_design
USE_PFASST = TRUE
FC_COMPILE = mpifort #gfortran -O3
LIBGFORTRAN = /usr/local/lib/gcc/11


### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := driver
src_dirs = fsrc
src_dirs += csrc

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
include .depend
include $(LIBPFASST)/.depend


tttt: tttt.txt
	echo "a" > tttt.txt

reallyclean: realclean libclean

libbike: 
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

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


