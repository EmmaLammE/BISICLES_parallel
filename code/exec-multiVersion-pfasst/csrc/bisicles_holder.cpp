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
	HHolder=HCrse; // HCrse is corretly intialized in different levels, passing in as pointer wil be updated in time
}

// deep copy of H holder into H holder backup
void AmrIceHolderClass::SetAmrHBackup(Vector<LevelData<FArrayBox>* > HCrse)
{
	Vector<LevelData<FArrayBox>* > Htemp(HCrse.size(), NULL);
	HHolderBackup = Htemp;

	for (int lvl=0; lvl<HCrse.size();lvl++)
    {
      const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl);
      HHolderBackup[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero); // FIXED the num of component=1, i.e. one d?? only Vx no Vy??
        
     DataIterator dit = (*HHolderBackup[lvl]).dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
         (*HHolderBackup[lvl])[dit].setVal(0.0);
      }
    }

	// Real test_min=0;
	// int nlvl=HCrse.size();
	// for (int lvl=0; lvl < nlvl; lvl++)
	// {
	// 	LevelData<FArrayBox>& ldf = *HCrse[lvl];
	// 	DisjointBoxLayout dbl = ldf.disjointBoxLayout();
	// 	DataIterator dit = ldf.dataIterator();
	// 	for (dit.reset(); dit.ok(); ++dit) 
	// 	{
	// 	const Box& box = dbl[dit()];
	// 	FArrayBox& fab = ldf[dit()];
	// 	test_min=fab.norm(box,1);
	// 	 cout<<"H backup after initialized "<<test_min<<endl;
	// 	} 
	// }
	
	// H Crse should be deep copied to HHolderBackup correctly
	for (int lvl=0; lvl<HCrse.size();lvl++)
	{
		const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl);
		// HHolderBackup[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero); // FIXED the num of component=1, i.e. one d?? only Vx no Vy??
		
		LevelData<FArrayBox>& ldf = *HHolderBackup[lvl]; // dest
		LevelData<FArrayBox>& xldf = *HCrse[lvl]; // src
		DisjointBoxLayout dbl = ldf.disjointBoxLayout();

		DataIterator dit = ldf.dataIterator();
		for (dit.begin(); dit.ok(); ++dit)
		{
			const Box& box = dbl[dit()];
			FArrayBox& fab = ldf[dit()];
			const FArrayBox& xfab = xldf[dit()];
			fab.copy(xfab,box);
		}
	}

	// for (int lvl=0; lvl < nlvl; lvl++)
	// {
	// 	LevelData<FArrayBox>& ldf = *HHolderBackup[lvl];
	// 	DisjointBoxLayout dbl = ldf.disjointBoxLayout();
	// 	DataIterator dit = ldf.dataIterator();
	// 	for (dit.reset(); dit.ok(); ++dit) 
	// 	{
	// 	const Box& box = dbl[dit()];
	// 	FArrayBox& fab = ldf[dit()];
	// 	test_min=fab.norm(box,1);
	// 	 cout<<"H backup after cpoied "<<test_min<<endl;
	// 	} 
	// }

    // int nlvl=HCrse.size();
    // for (int lvl=0; lvl < nlvl; lvl++)
    // {
    //   LevelData<FArrayBox>& ldf = *HCrse[lvl];
    //   LevelData<FArrayBox>& xldf = *HHolderBackup[lvl];
    //   DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    //   DataIterator dit = ldf.dataIterator();
    //   for (dit.reset(); dit.ok(); ++dit) 
    //      {
    //      const Box& box = dbl[dit()];
    //      FArrayBox& fab = ldf[dit()];
    //      const FArrayBox& xfab = xldf[dit()];
    //      fab.copy(xfab,box);
    //      } 
    // }
}

Vector<LevelData<FArrayBox>* > AmrIceHolderClass::GetAmrH(void)
{
	return HHolder;
}

Vector<LevelData<FArrayBox>* > AmrIceHolderClass::GetAmrHBackup(void)
{
	return HHolderBackup;
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

