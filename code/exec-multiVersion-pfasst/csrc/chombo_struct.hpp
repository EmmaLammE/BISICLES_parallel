#ifndef CHOMBO_STRUCT_HPP
#define CHOMBO_STRUCT_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
// #include "_hypre_struct_mv.h"
// #include "_hypre_struct_ls.h"
// #include "HYPRE_struct_mv.h"
// #include "HYPRE_struct_ls.h"
#include <mpi.h>
#include "AmrIce.H"
#include "SundialsUtil.H"

using namespace std;

// ice sheet state, pointer to the amrice object
class ChomboStruct {
   protected:
      // chombo struct 
   	  // Vector<LevelData<FArrayBox>* > ice_thickness;
      //   AmrIce amrObjectCrse;

   public:

   	  // mpi stuffs
   	  MPI_Comm comm;
   	  int myid, num_procs;

   	  // chombo struct 
   	  int num_of_rows;



      ChomboStruct() {
         ;
      }
      ~ChomboStruct() {
         ;
      }

      void SetComm(MPI_Comm target_comm);

      void InitGrid(int num_grid_points);
      int GetNumRows(void); // get the num of rows of A matric, in 2d->grid_size^2
      //void PackPF2CH(double *ice_thickness_ptr_unwrapped); // wrap ptr into chombo struct
      // void SetDim(int in_dim);
      // int GetNumRows(void);
      // void SetInitCond(double val);
      // void InitGrid(int num_grid_points, int in_nrows = -1, int *extents = NULL);
      // double *HeatEquTrueSol(double t, int P, int Q, double init_cond);
};

#endif