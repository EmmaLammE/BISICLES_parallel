#ifndef PFASST_BISICLES_HPP
#define PFASST_BISICLES_HPP

#include "bisicles_holder.hpp"
#include "AmrIce.H"
#include "SundialsUtil.H"

using namespace std;

void Pf_Bisicles_setHolders(AmrIce* amrObjectPtr,
                            AmrIceHolderClass* AmrIceHolderPtr,         
                            const Vector<LevelData<FArrayBox>* >& H,           
                            const Vector<LevelData<FArrayBox>* >& Vel);

#endif
