#ifndef PFASST_BISICLES_HPP
#define PFASST_BISICLES_HPP

#include "bisicles_holder.hpp"
#include "AmrIce.H"
#include "SundialsUtil.H"

using namespace std;

void Pf_Bisicles_setHolders(AmrIce* amrObjectPtr,\
							AmrIceHolderClass* AmrIceHolderPtr, \
							Vector<LevelData<FArrayBox>* > H,\
							Vector<LevelData<FArrayBox>* > Vel);

#endif