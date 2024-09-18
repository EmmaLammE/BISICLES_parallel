#include "pfasst_bisicles_setup.hpp"
#include "bisicles_holder.hpp"

// H passed in here is only for the coarsest level, not a vector of level data
void Pf_Bisicles_setHolders(AmrIce* amrObjectPtr,
                            AmrIceHolderClass* AmrIceHolderPtr,         
                            const Vector<LevelData<FArrayBox>* >& H,          
                            const Vector<LevelData<FArrayBox>* >& Vel)
{
    AmrIceHolderPtr->SetAmrIceObjPtr(amrObjectPtr);

    IceSheetState iceState;
    reshapeAndFill(iceState.ice_thickness, H);
    reshapeAndFill(iceState.ice_velocity,Vel);
    AmrIceHolderPtr->SetAmrIceStatePtr(&iceState);
    AmrIceHolderPtr->SetAmrIceState(iceState);

    // crsedHdtVect
    Vector<LevelData<FArrayBox>* > dHdtVect; // double check with Dan
    reshape(dHdtVect,H);
    AmrIceHolderPtr->SetAmrdHdtPtr(&dHdtVect);
    AmrIceHolderPtr->SetAmrdHdt(dHdtVect);

    // crseH
    Vector<LevelData<FArrayBox>* > nonConstH = H;
    AmrIceHolderPtr->SetAmrHPtr(&nonConstH);
    AmrIceHolderPtr->SetAmrH(nonConstH);
    // 9/17/24 -- Emma says it's safe to comment this out
    // AmrIceHolderPtr->SetAmrHBackup(nonConstH);
    // Vector<LevelData<FArrayBox>* > temp = AmrIceHolderPtr->GetAmrH(); // amr H is set corretely in all levels

    // set up max num of spatial level
    AmrIceHolderPtr->SetAmrNumLvl(H.size());
}
