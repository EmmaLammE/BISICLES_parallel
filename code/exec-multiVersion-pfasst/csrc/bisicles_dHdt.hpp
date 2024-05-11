#ifndef BISICLES_DHDT_SOLVER_HPP
#define BISICLES_DHDT_SOLVER_HPP

#include <iostream>
#include <cstdlib>
// #include "chombo_struct.hpp"
#include "bisicles_holder.hpp"

using namespace std;

class BisiclesdHdtSolver :  public AmrIceHolderClass{
   protected:

   public:
      

      BisiclesdHdtSolver(int num_grid_points=16,\
                         AmrIceHolderClass *c_AmrIceHolderPtr=nullptr);
      ~BisiclesdHdtSolver(void);
      void Cleanup(void);


      void FEval(double t, int level_index, \
                 Vector<LevelData<FArrayBox>* > *rhs, \
                 Vector<LevelData<FArrayBox>* > *dHdtVect);
   };


#endif

