#include "bisicles_holder.hpp"
#include "bisicles_dHdt.hpp"


// --------------- Amr Ice Object ---------------------//
void AmrIceHolderClass::SetAmrIceObjPtr(AmrIce *amrObjectCrsePtr)
{
	amrObjHolderPtr=amrObjectCrsePtr;
}

AmrIce *AmrIceHolderClass::GetAmrIceObjPtr(void)
{
	return amrObjHolderPtr;
}



// --------------- Amr Ice state ---------------------//
void AmrIceHolderClass::SetAmrIceStatePtr(IceSheetState *amrStateCrsePtr)
{
	amrStateHolderPtr=amrStateCrsePtr;
}

IceSheetState *AmrIceHolderClass::GetAmrIceStatePtr(void)
{
	return amrStateHolderPtr;
}

void AmrIceHolderClass::SetAmrIceState(IceSheetState amrStateCrse)
{
	amrStateHolder=amrStateCrse;
}

IceSheetState AmrIceHolderClass::GetAmrIceState(void)
{
	return amrStateHolder;
}


// --------------- Amr Ice dHdt ---------------------//
void AmrIceHolderClass::SetAmrdHdtPtr(Vector<LevelData<FArrayBox>* > *dHdtCrsePtr)
{
	dHdtHolderPtr=dHdtCrsePtr;
}

Vector<LevelData<FArrayBox>* > *AmrIceHolderClass::GetAmrdHdtPtr(void)
{
	return dHdtHolderPtr;
}

void AmrIceHolderClass::SetAmrdHdt(Vector<LevelData<FArrayBox>* > dHdtCrse)
{
	dHdtHolder=dHdtCrse;
}

Vector<LevelData<FArrayBox>* > AmrIceHolderClass::GetAmrdHdt(void)
{
	return dHdtHolder;
}

// --------------- Amr Ice H vector ---------------------//
void AmrIceHolderClass::SetAmrHPtr(Vector<LevelData<FArrayBox>* > *HCrsePtr)
{
	HHolderPtr=HCrsePtr;
}

Vector<LevelData<FArrayBox>* > *AmrIceHolderClass::GetAmrHPtr(void)
{
	return HHolderPtr;
}

void AmrIceHolderClass::SetAmrH(Vector<LevelData<FArrayBox>* > HCrse)
{
	HHolder=HCrse;
}

Vector<LevelData<FArrayBox>* > AmrIceHolderClass::GetAmrH(void)
{
	return HHolder;
}


// --------------- Amr Ice H vector ---------------------//
void AmrIceHolderClass::SetAmrNumLvl(int max_num_lvl)
{
	num_of_levels=max_num_lvl;
}


int AmrIceHolderClass::GetAmrNumLvl(void)
{
	return num_of_levels;
}

