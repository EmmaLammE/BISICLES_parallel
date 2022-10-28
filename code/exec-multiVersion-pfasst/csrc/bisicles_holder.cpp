#include "bisicles_holder.hpp"
#include "bisicles_dHdt.hpp"


// --------------- Amr Ice Object ---------------------//
void AmrIceHolderClass::SetAmrIceObjPtr(AmrIce *amrObjectCrsePtr)
{
	// cout << "bisicles_pfasst.cpp type of amrObjectCrsePtr " << amrObjectCrsePtr << std::endl;
	// cout << "bisicles_pfasst.cpp amrHolderPtr " << amrHolderPtr << std::endl;
	// amrObjectCrse.setParmParsePrefix("crse.");
	// amrObjectCrsePtr->setParmParsePrefix("fine");
	amrObjHolderPtr=amrObjectCrsePtr;

	cout << "bisicles_pfasst.cpp 0000 amrHolderPtr " << amrObjHolderPtr << std::endl;
}

AmrIce *AmrIceHolderClass::GetAmrIceObjPtr(void)
{
	// cout << "bisicles_pfasst.cpp  0000"<<std::endl;
	// cout << "bisicles_pfasst.cpp 1111 amrHolderPtr " << amrHolderPtr << std::endl;
	return amrObjHolderPtr;
}

// void AmrIceHolderClass::SetAmrIceObj(AmrIce amrObjectCrse)
// {
// 	amrObjHolder=amrObjectCrse;
// }

// AmrIce AmrIceHolderClass::GetAmrIceObj(void)
// {
// 	return amrObjHolder;
// }


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
	// cout << "bisicles_holder.cpp HHolder " << HHolder.size() << std::endl;
	// HHolder=0x10404f028
	return HHolder;
	// cout << "bisicles_holder.cpp done HHolder assign " << HHolder.size() << std::endl;
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


// // --------------- Amr Ice num of cells vector ---------------------//
// void AmrIceHolderClass::SetAmrNumCells(int max_num_lvl)
// {
// 	num_of_levels=max_num_lvl;
// }


// int AmrIceHolderClass::GetAmrNumCells(void)
// {
// 	return num_of_levels;
// }


