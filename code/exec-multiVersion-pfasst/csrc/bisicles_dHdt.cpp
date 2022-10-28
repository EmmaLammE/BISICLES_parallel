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
	// 1. wrap pfasst data type of y(H) to chombo data type
	// 2. pass chombo data to dHdt
	// 3. compute dHdt by calling dHdt
	// 4. unwrap chombo data dHdt(f) to pfasst data type
	// 5. f->f_values
	// y->crseH

	
	//1. type y: pd; type f: ppd
	// cout << "bisicles_dHdt.cpp 0000 type of f " << typeid(f).name() << std::endl;
	// cout << "bisicles_dHdt.cpp 0000 length of f " << sizeof(*y) << std::endl;

	// y->PackPf2CH(double *y)

	*rhs=*dHdtVect;
	
}


// void BisiclesdHdtSolver::ComputedHdt(AmrIce amrObjectCrse,Vector<Vector<LevelData<FArrayBox>* > > crsedHdtVect,\
//          Vector<LevelData<FArrayBox>* > crseH)
// {
	
// 	amrObjectCrse.setParmParsePrefix("crse.");
// 	cout << "holding amr ice objects " << typeid(amrObjectCrse).name() << std::endl;
	
// }


// AmrIce BisiclesdHdtSolver::getBisiclesObjPtr()
// {
// 	call amrObject
// 	return amrObjectCrse	
// }


// void AmrIceHolder(AmrIce *amrObjectCrsePtr)
// {
// 	//cout << "holding amr ice objects " << typeid(amrObjectCrsePtr).name() << std::endl;
// 	// amrObjectCrse.setParmParsePrefix("crse.");
// 	// amrObjectCrsePtr->setParmParsePrefix("fine");
// }