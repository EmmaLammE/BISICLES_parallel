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
include $(LIBPFASST)/.depend
FFLAGS  += -I$(LIBPFASST)/include
FFLAGS  += -I$(LIBPFASST)/build 
FFLAGS  += -I/share/software/user/open/openmpi/4.1.2/lib
FFLAGS  += -MD -MP -cpp
FFLAGS += -I./include
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

# fixed fsrc dependency
PFASST_APP_UNITS := x_probin y_encap z_sweeper v_level u_hooks w_comm s_pfasst_bisicles t_pfasst_main

# Use the object paths selected by Chombo for linking.
pfasst_app_object = $(if $(strip $(objects)),$(or $(firstword $(filter %/$(1).o $(1).o,$(objects))),$(error No object for $(1) in Chombo objects)),$(OBJ_DIR)/$(1).o)

PFASST_APP_OBJECTS := $(foreach unit,$(PFASST_APP_UNITS),$(call pfasst_app_object,$(unit)))

# Store all application .mod files in one dedicated directory.
PFASST_APP_MODDIR := $(abspath $(dir $(firstword $(PFASST_APP_OBJECTS)))/pfasst-modules)

PFASST_APP_MODULES := probin encap pf_my_sweeper pf_my_level hooks pf_space_comm pfasst_bisicles pfasst_main
PFASST_APP_MODFILES := $(addprefix $(PFASST_APP_MODDIR)/,$(addsuffix .mod,$(PFASST_APP_MODULES)))

# Find application modules before similarly named external modules.
$(PFASST_APP_OBJECTS): FFLAGS := -I$(PFASST_APP_MODDIR) $(FFLAGS) #-J$(PFASST_APP_MODDIR)
$(PFASST_APP_OBJECTS): | $(PFASST_APP_MODDIR)

$(PFASST_APP_MODDIR): ; mkdir -p "$@"

# Rebuild if any generated application module is missing.
ifneq ($(words $(wildcard $(PFASST_APP_MODFILES))),$(words $(PFASST_APP_MODFILES)))
.PHONY: pfasst-missing-module
$(PFASST_APP_OBJECTS): pfasst-missing-module
endif

# Dependencies between the current Fortran source files.
$(call pfasst_app_object,x_probin): fsrc/x_probin.f90

$(call pfasst_app_object,y_encap): fsrc/y_encap.f90 \
    $(call pfasst_app_object,x_probin)

$(call pfasst_app_object,z_sweeper): fsrc/z_sweeper.f90 \
    $(call pfasst_app_object,y_encap)

$(call pfasst_app_object,v_level): fsrc/v_level.f90 \
    $(call pfasst_app_object,y_encap)

$(call pfasst_app_object,u_hooks): fsrc/u_hooks.f90 \
    $(call pfasst_app_object,y_encap) \
    $(call pfasst_app_object,z_sweeper)

$(call pfasst_app_object,w_comm): fsrc/w_comm.f90

$(call pfasst_app_object,s_pfasst_bisicles): fsrc/s_pfasst_bisicles.f90 \
    $(call pfasst_app_object,x_probin) \
    $(call pfasst_app_object,y_encap) \
    $(call pfasst_app_object,z_sweeper) \
    $(call pfasst_app_object,v_level)

$(call pfasst_app_object,t_pfasst_main): fsrc/t_pfasst_main.f90 \
    $(call pfasst_app_object,x_probin) \
    $(call pfasst_app_object,y_encap) \
    $(call pfasst_app_object,z_sweeper) \
    $(call pfasst_app_object,v_level) \
    $(call pfasst_app_object,u_hooks) \
    $(call pfasst_app_object,w_comm) \
    $(call pfasst_app_object,s_pfasst_bisicles)

.PHONY: pfasst
pfasst: $(PFASST_APP_OBJECTS)


.PHONY: all lib libbike libclean alltest

# Ordered modules; also rebuild after LibPFASST interface changes.
include Makefile.fortran-modules
