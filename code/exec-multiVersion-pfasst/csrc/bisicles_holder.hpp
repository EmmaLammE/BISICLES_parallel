#ifndef BISICLES_HOLDER_HPP
#define BISICLES_HOLDER_HPP


#include "AmrIce.H"
#include "SundialsUtil.H"

using namespace std;


class AmrIceHolderClass{
protected:
	AmrIce *amrObjHolderPtr;
	AmrIce amrObjHolder;
	IceSheetState *amrStateHolderPtr;
	IceSheetState amrStateHolder;
	Vector<LevelData<FArrayBox>* > *dHdtHolderPtr;
	Vector<LevelData<FArrayBox>* > dHdtHolder;
	Vector<LevelData<FArrayBox>* > *HHolderPtr;
	Vector<LevelData<FArrayBox>* > HHolder;
	Vector<LevelData<FArrayBox>* > HHolderBackup;
	int num_of_levels;
	// Vector<int> m_num_cells;

public:

	
	// AmrIce object;
	void SetAmrIceObjPtr(AmrIce *amrObjectCrsePtr);
	AmrIce *GetAmrIceObjPtr(void);

	// AmrIce ice state
	void SetAmrIceStatePtr(IceSheetState *amrStateHolderPtr);
	IceSheetState *GetAmrIceStatePtr(void);
	void SetAmrIceState(IceSheetState amrStateHolder);
	IceSheetState GetAmrIceState(void);

	// AmrIce dHdt vector
	void SetAmrdHdtPtr(Vector<LevelData<FArrayBox>* > *dHdtHolderPtr);
	Vector<LevelData<FArrayBox>* > *GetAmrdHdtPtr(void);
	void SetAmrdHdt(Vector<LevelData<FArrayBox>* > dHdtHolder);
	Vector<LevelData<FArrayBox>* > GetAmrdHdt(void);

	// AmrIce H vector
	void SetAmrHPtr(Vector<LevelData<FArrayBox>* > *HHolderPtr);
	Vector<LevelData<FArrayBox>* > *GetAmrHPtr(void);
	void SetAmrH(Vector<LevelData<FArrayBox>* > HHolder);
	void SetAmrHBackup(Vector<LevelData<FArrayBox>* > HHolder); // copy by value not pointer for H just in case
	Vector<LevelData<FArrayBox>* > GetAmrH(void);
	Vector<LevelData<FArrayBox>* > GetAmrHBackup(void);

	// AmrIce set up max num of levels
	void SetAmrNumLvl(int max_num_lvl);
	int GetAmrNumLvl(void);

};

#endif
