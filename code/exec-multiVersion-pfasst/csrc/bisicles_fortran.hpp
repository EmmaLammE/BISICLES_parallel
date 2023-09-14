#ifndef BISICLES_FORTRAN_HPP
#define BISICLES_FORTRAN_HPP

#include "bisicles_vector.hpp"
#include "bisicles_dHdt.hpp"
#include "bisicles_holder.hpp"
#include <mpi.h>


extern "C"
{
   // bisicles Vector functions
   void BisiclesVectorCreate(BisiclesVector **bisicles_vector,int num_grid_points, \
                              AmrIceHolderClass *c_AmrIceHolderPtr);
   void BisiclesVectorDestroy(BisiclesVector *bisicles_vector);
   void BisiclesVectorSetVal(BisiclesVector *bisicles_vector, double y, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);
   void BisiclesVectorCopy(BisiclesVector *dest, BisiclesVector *src, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);
   int BisiclesCurrentVectorSize(BisiclesVector *bisicles_vector);
   Real* BisiclesVectorPack(BisiclesVector *bisicles_vector,\
                          AmrIceHolderClass *c_AmrIceHolderPtr, int level_id);
   void BisiclesVectorUnpack(BisiclesVector *bisicles_vector, Real* flattened_array, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);
   double BisiclesVectorNorm(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);
   void BisiclesVectorAxpy(BisiclesVector *y, double a, BisiclesVector *x, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);
   void BisiclesVectorL2Print(BisiclesVector *bisicles_vector);
   void BisiclesVectorLevelDataBox(BisiclesVector *bisicles_vector, BisiclesVector *x);
   double BisiclesVectorGetVal(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);


   // bisicles solver functions
   void BisiclesSolverInit(BisiclesdHdtSolver **BisiclesdHdtSolver,\
                        int pfasst_level_index,\
                        int num_grid_points,\
                        AmrIceHolderClass *c_AmrIceHolderPtr);
   void  BisiclesVectorSetHIC(BisiclesVector *y_0,\
                        AmrIceHolderClass *c_AmrIceHolderPtr);

   // void BisiclesSolverFEval(BisiclesdHdtSolver *bisicles_dHdt, HypreVector *y, double t, int pfasst_level_index, HypreVector *f);
   void BisiclesSolverFEval(BisiclesdHdtSolver *BisiclesdHdtSolver, BisiclesVector *y, double t, \
      int pfasst_level_index, BisiclesVector *f, double dt,int maxStep,bool evolve_velocity, AmrIceHolderClass *c_AmrIceHolderPtr);

   void PfasstBisiclesSaveResults(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr);
   
   // print out for debugging
   void PfasstPrintAmr(BisiclesVector *y, AmrIceHolderClass *c_AmrIceHolderPtr);

   // MPI comm transfer
   MPI_Fint PfasstBisicles_MPI_Comm_c2f(MPI_Comm *comm);\
}

   #endif