# -*- Mode: Makefile -*- 

LIBGFORTRAN = /usr/local/lib/gcc/11

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


