#include "bisicles_dHdt.hpp"
#include "bisicles_holder.hpp"



BisiclesdHdtSolver::BisiclesdHdtSolver(MPI_Comm in_comm,int num_grid_points,\
						  AmrIceHolderClass *c_AmrIceHolderPtr)
{
	InitGrid(num_grid_points);
	/* Emma TODO
   		add more functions
   		not sure what to add yet
   	   Emma END TODO
    */
}

BisiclesdHdtSolver::~BisiclesdHdtSolver(void)
{
   Cleanup();
}

void BisiclesdHdtSolver::Cleanup(void)
{
   /* Emma TODO
   		add clean up functions
   	  Emma END TOD
   */
}


void BisiclesdHdtSolver::FEval(double t, int level_index, \
										Vector<LevelData<FArrayBox>* > *rhs, \
										Vector<LevelData<FArrayBox>* > *dHdtVect)
{
	*rhs=*dHdtVect;
	
}

