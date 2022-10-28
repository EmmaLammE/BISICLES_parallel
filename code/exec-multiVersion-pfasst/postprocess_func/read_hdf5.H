using namespace std;

#include "AmrIce.H"
#include "SundialsUtil.H"

using namespace std;

void read_hdf5_to_leveldata(AmrIce* amrObjectPtr,\
							AmrIceHolderClass* AmrIceHolderPtr, \
							Vector<LevelData<FArrayBox>* > H,\
							Vector<LevelData<FArrayBox>* > Vel);