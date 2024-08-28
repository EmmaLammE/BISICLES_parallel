# Global definitions
FC = mpifort
FFLAGS = -g -Wno-line-truncation -fno-second-underscore -Wno-unused-parameter -Wno-unused-variable -m64 -DCH_MPI -I/home/groups/jsuckale/liuwj/PFASST_icesheet//include -I/home/groups/jsuckale/liuwj/PFASST_icesheet//build -cpp -MD -MP "-Wl,--fno-underscore" "-Wl,-no_compact_unwind" -ffree-line-length-none

OBJ_DIR = ../o/2d.Linux.64.mpicxx.mpifort.DEBUG.MPI
SOURCES = probin.f90 \
		  comm.f90 \
		  encap.f90 \
		  sweeper.f90 \
		  hooks.f90 \
		  levels.f90 \
		  pfasst_bisicles.f90 \
		  pfasst_main.f90
OBJECTS = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SOURCES))


all: clean $(OBJ_DIR)

$(OBJ_DIR): $(OBJECTS)
	$(FC) $(FFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: %.f90
	$(info Compiling: $@ from $<)
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -c $< -o $@

clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.d $(OBJ_DIR)/application

.PHONY: clean all

# Include dependencies
include .depend
