#ifndef CHOMBO_STRUCT_HPP
#define CHOMBO_STRUCT_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <mpi.h>
#include "AmrIce.H"
#include "SundialsUtil.H"

using namespace std;

// ice sheet state, pointer to the amrice object
class ChomboStruct {
   protected:

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
};

#endif