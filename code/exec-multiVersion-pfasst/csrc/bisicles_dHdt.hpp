#ifndef BISICLES_DHDT_SOLVER_HPP
#define BISICLES_DHDT_SOLVER_HPP

#include <iostream>
#include <cstdlib>
#include "chombo_struct.hpp"
#include "bisicles_holder.hpp"

using namespace std;

class BisiclesdHdtSolver : public ChomboStruct, public AmrIceHolderClass{
   protected:
      ;

   public:
      
      // AmrIce amrObjectCrse;
      // Vector<Vector<LevelData<FArrayBox>* > > crsedHdtVect;
      // Vector<LevelData<FArrayBox>* > crseH;


      BisiclesdHdtSolver(MPI_Comm in_comm = MPI_COMM_WORLD, int num_grid_points=16,\
                         AmrIceHolderClass *c_AmrIceHolderPtr=nullptr);
      ~BisiclesdHdtSolver(void);
      void Cleanup(void);


      void FEval(double t, int level_index, \
                 Vector<LevelData<FArrayBox>* > *rhs, \
                 Vector<LevelData<FArrayBox>* > *dHdtVect);
      // void ComputedHdt(AmrIce *amrObjectCrse,Vector<Vector<LevelData<FArrayBox>* > > crsedHdtVect,\
      //    Vector<LevelData<FArrayBox>* > crseH);
      // AmrIce getBisiclesObjPtr();
      // void AmrIceHolder(AmrIce *amrObjectCrsePtr);

   };


#endif