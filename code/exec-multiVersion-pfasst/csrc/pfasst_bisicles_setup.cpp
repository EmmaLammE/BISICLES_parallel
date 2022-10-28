#include "pfasst_bisicles_setup.hpp"
#include "bisicles_holder.hpp"

void Pf_Bisicles_setHolders(AmrIce* amrObjectPtr,\
							AmrIceHolderClass* AmrIceHolderPtr, \
							Vector<LevelData<FArrayBox>* > H,\
							Vector<LevelData<FArrayBox>* > Vel)
{
	// AmrIce *amrObjectCrsePtr=amrObjectPtr;
 //    AmrIceHolderClass AmrIceHolderPtr;
    AmrIceHolderPtr->SetAmrIceObjPtr(amrObjectPtr);
    // AmrIceHolderPtr->SetAmrIceObj(*amrObject);

    // crseStateVect: 
    // double check with Dan: if anything else needs to be setup before passing to holder
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
    // cout<< "0000  pfasst_bisicles_setup.cpp size of reshaped dHdtVector "<<dHdtVect.size()<< std::endl;

    // crseH
    AmrIceHolderPtr->SetAmrHPtr(&H);
    AmrIceHolderPtr->SetAmrH(H);
    // cout<< "0000  pfasst_bisicles_setup.cpp size of reshaped H "<<H.size()<< std::endl;


    // set up max num of spatial level
    AmrIceHolderPtr->SetAmrNumLvl(H.size());
}
