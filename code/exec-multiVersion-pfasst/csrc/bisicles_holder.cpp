#include "bisicles_holder.hpp"
#include "bisicles_dHdt.hpp"
#include <iomanip>
#include <iostream>



// initialize count and run time
int AmrIceHolderClass::call_count_pack = 0;  // Initialize the static int
std::chrono::duration<double> AmrIceHolderClass::total_duration_pack(0);  // Initialize the static duration
int AmrIceHolderClass::call_count_unpack = 0;  // Initialize the static int
std::chrono::duration<double> AmrIceHolderClass::total_duration_unpack(0);  // Initialize the static duration
int AmrIceHolderClass::call_count_setstate = 0;  // Initialize the static int
std::chrono::duration<double> AmrIceHolderClass::total_duration_setstate(0);  // Initialize the static duration
int AmrIceHolderClass::call_count_dHdt = 0;  // Initialize the static int
std::chrono::duration<double> AmrIceHolderClass::total_duration_dHdt(0);  // Initialize the static duration

//  --------------- update call time and call counts ---------------------//
void AmrIceHolderClass::updateCallTimePack(std::chrono::duration<double> duration) {
    call_count_pack++;
    total_duration_pack += duration;
}
void AmrIceHolderClass::updateCallTimeUnpack(std::chrono::duration<double> duration) {
    call_count_unpack++;
    total_duration_unpack += duration;
}
void AmrIceHolderClass::updateCallTimeFEval(std::chrono::duration<double> duration_setstate,std::chrono::duration<double> duration_dHdt) {
    call_count_setstate++;
	call_count_dHdt++;
    total_duration_setstate += duration_setstate;
	total_duration_dHdt += duration_dHdt;
}
void AmrIceHolderClass::printStatistics(int number_procs,int new_time_size,int new_space_size,int new_time_rank,\
                                        int new_space_rank, std::chrono::duration<double> duration_pf,\
										bool pf_evolve_velocity,int numCrseIntervals,Real maxTime) {
	cout << "\n#################################################################################\n";
	cout << "Summary of performance statistics\n-------------------------------"<<endl;
	cout << "  Total num of procs       = "<< number_procs <<endl;
	cout << "  Total num of time procs  = "<< new_time_size <<endl;
	cout << "  Total num of space procs = "<< new_space_size <<endl;
	cout << "  Proc rank for this instance:  time = "<< new_time_rank << "  space =  "<< new_space_rank <<endl;
	cout << "Parameters that influence the run speeds (a lot)"<<endl;
	cout << "  Compute velocity? = "<< pf_evolve_velocity << endl;
	cout << "  Num of time steps = "<< numCrseIntervals << endl;
	cout << "  Final time        = "<< maxTime<<" years\n";
	cout << "  \nTotal simulation time using PFASST: ! "<< duration_pf.count() <<" ! sec"<<endl;
	cout << setfill('_') << setw(104) << "_" << endl << setfill(' ');
    cout << "| " << left << setw(22) << "Function name (short)" << "| " 
         << right << setw(18) << "call counts" << " | " 
         << right << setw(24) << "Total duration (sec)" << " | " 
         << right << setw(25) << "Ratio of total time (%)" << " |" << endl;
    cout << setfill('-') << setw(104) << "-" << endl << setfill(' ');

    cout << "| " << left << setw(22) << "Pack" << "| " 
         << right << setw(18) << call_count_pack << " | " 
         << right << setw(24) << fixed << setprecision(8) << total_duration_pack.count() << " | " 
         << right << setw(25) << fixed << setprecision(8) << (total_duration_pack.count()/duration_pf.count()*100) << " |" << endl;

    cout << "| " << left << setw(22) << "Unpack" << "| " 
         << right << setw(18) << call_count_unpack << " | " 
         << right << setw(24) << fixed << setprecision(8) << total_duration_unpack.count() << " | " 
         << right << setw(25) << fixed << setprecision(8) << (total_duration_unpack.count()/duration_pf.count()*100) << " |" << endl;

    cout << "| " << left << setw(22) << "setState(FEval)" << "| " 
         << right << setw(18) << call_count_setstate << " | " 
         << right << setw(24) << fixed << setprecision(8) << total_duration_setstate.count() << " | " 
         << right << setw(25) << fixed << setprecision(8) << (total_duration_setstate.count()/duration_pf.count()*100) << " |" << endl;

    cout << "| " << left << setw(22) << "dHdt(FEval)" << "| " 
         << right << setw(18) << call_count_dHdt << " | " 
         << right << setw(24) << fixed << setprecision(8) << total_duration_dHdt.count() << " | " 
         << right << setw(25) << fixed << setprecision(8) << (total_duration_dHdt.count()/duration_pf.count()*100) << " |" << endl;

    cout << setfill('_') << setw(104) << "_" << endl << setfill(' ');
    cout << setfill('-') << setw(104) << "-" << endl << setfill(' ');
}

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

	
	// H Crse should be deep copied to HHolderBackup correctly
	for (int lvl=0; lvl<HCrse.size();lvl++)
	{
		const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl);
		
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

