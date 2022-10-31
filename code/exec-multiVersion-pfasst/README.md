# BISICLES-PFASST example #
This is an example of simulating a simple ice stream using BISICLES combined with PFASST parallel in time scheme.

It contains mainly 2 additional src codes folders:

 - [`/fsrc`](#fsrc): The Fortran source codes;
 - [`/csrc`](#csrc): The Cpp source codes;
 
with some other folders and files:
- [`chombo_pfasst.hpp`](#chombo_pfasst.hpp): Binding file for the `driver.cpp` and the Fortran world;
- [`inputs.test`](#inputs.test): Input files for BISICLES (added some params for PFASST);
- [`pfasst_params.nml`](#pfasst_params.nml): Input files for PFASST;
- [`/data_convergence_test`](#data_convergence_test): Post processing ChomboCompare codes, plotting codes, and results data.

## Quick Run ##
Before compiling, please make sure the followings are installed and set:
- [PFASST](#PFASST): [Download PFASST](https://github.com/EmmaLammE/PFASST_icesheet). The PFASST library should be installed in the same directory as BISICLES.
- [`GNUMakefile`](#GNUMakefile): the only thing needs to be reset is the path of PFASST `LIBPFASST` (I think).

## Some other notes ##
If nothing is changed, the current codes will run 2 test: parallel and serial. Parallel one uses the objects and states of the `crse` settings in `inputs.test`, and Serial one uses that of the `fine` settings. The results are stored in the same folder as `pf_....` and `ref_....`.
