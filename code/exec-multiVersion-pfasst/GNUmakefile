# -*- Mode: Makefile -*- 

####################
# User-set variables
####################
# PFASST include directory
#LIBPFASST = /Users/emmaliu/internship2022/installs/PFASST_design
USE_PFASST = TRUE
FC_COMPILE = mpifort #gfortran -O3
#LIBGFORTRAN = /usr/local/lib/gcc/11


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
FFLAGS  += -I./include
#LDFLAGS += -lstdc++ -lmpi_mpifh -lgcc -lquadmath
LDFLAGS += -lstdc++  -lgcc 
LDFLAGS +=  -L$(LIBPFASST)/lib -lpfasst 
#XTRALDFLAGS += -lquadmath
FFLAGS  += -ffree-line-length-none
#include .depend
#include $(LIBPFASST)/.depend


tttt: tttt.txt
	echo "a" > tttt.txt

reallyclean: realclean libclean

libbike: 
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS) -L$(LIBPFASST)/lib -lpfasst -L$(LIBGFORTRAN) -lgfortran
LIBFLAGS += -I$(LIBPFASST)/build -I$(LIBPFASST)/include -I/include
LIBFLAGS += -lstdc++ $(LIBPFASST)/lib/libpfasst.a -lquadmath -ldl -lmpichfort -lmpich

# add specific dependencies for the stuff in fsrc, because it needs to be built in order
driver: clean libbike ../lib/libBisicles$(config).a myStuff


#build stuff in fsrc specially because it needs to be built in order
myStuff:
	$(MAKE) o/$(config)/a_probin.o
	$(MAKE) o/$(config)/b_encap.o
	$(MAKE) o/$(config)/c_level.o
	$(MAKE) o/$(config)/comm.o
	$(MAKE) o/$(config)/d_sweeper.o
	$(MAKE) o/$(config)/e_hooks.o
	$(MAKE) o/$(config)/pfasst_bisicles.o
	$(MAKE) o/$(config)/pfasst_main.o
















alltest: libbike ../lib/libBisicles$(config).a 


#include Makefile.fortran
pfasst: 
	mkdir -p include
	mv *.mod include/


.PHONY: all lib libbike libclean alltest


