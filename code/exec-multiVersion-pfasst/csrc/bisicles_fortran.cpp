#include "bisicles_fortran.hpp"
#include "bisicles_vector.hpp"
#include "bisicles_holder.hpp"


#include "Misc.H"
#include "MayDay.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include "SPACE.H"
#include "BaseFabMacros.H"
#include "Box.H"
#include "BaseFab.H"
#include "REAL.H"
#include "SPACE.H"

#include "AmrIce.H"
#include "AmrIceBase.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelSigmaCS.H"
#include "IceVelocitySolver.H"
#include "RealVect.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "SurfaceFlux.H"
#include "BasalFriction.H"
#include "IceThicknessIBC.H"
#include "IceInternalEnergyIBC.H"
#include "CalvingModel.H"
#include "MuCoefficient.H"
#include "CH_HDF5.H"
#include "DomainDiagnosticData.H"
#include "ParmParse.H"
#include "SundialsUtil.H"
#include <sstream>
#ifdef SUNDIALS
#include "ChomboNVector.H"
#include <arkode/arkode_arkstep.h>    /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_parallel.h>   /* parallel N_Vector types, macros */
#include <sunlinsol/sunlinsol_pcg.h>  /* access to PCG SUNLinearSolver        */
#include <sundials/sundials_types.h>  /* defs. of realtype, sunindextype, etc */
#endif

#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <chrono>  // For high-resolution clock
using namespace std;


MPI_Comm glob_space_comm = MPI_COMM_NULL;

extern "C"
{

// bisicles vector functions
void BisiclesVectorCreate(BisiclesVector **bisicles_vector,int num_grid_points, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      *bisicles_vector = new BisiclesVector(num_grid_points,c_AmrIceHolderPtr);
   }

   void BisiclesVectorDestroy(BisiclesVector *bisicles_vector)
   {
      if (bisicles_vector != nullptr){
         delete bisicles_vector;
      }
   }

   void BisiclesVectorSetVal(BisiclesVector *bisicles_vector, double y, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      bisicles_vector->HSetVal2All(y,c_AmrIceHolderPtr);
   }

   void BisiclesVectorCopy(BisiclesVector *dest, BisiclesVector *src, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      dest->SetHVector(src);
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
   }

int BisiclesCurrentVectorSize(BisiclesVector *bisicles_vector)
   {
      const Vector<LevelData<FArrayBox>* >& test=bisicles_vector->GetHVector();
      
      int num_lvl = test.size();
      int num_total_cells = 0;
      int num_cells_per_lvl[num_lvl];
      for (int lvl=0; lvl < num_lvl; lvl++)
      {
         LevelData<FArrayBox>& ldf = *test[lvl];
         DisjointBoxLayout dbl = ldf.disjointBoxLayout();
         // some sizes
         num_cells_per_lvl[lvl] = dbl.numCells();
         num_total_cells += num_cells_per_lvl[lvl];
      }
      return num_total_cells;
   }

Real* BisiclesVectorPack(BisiclesVector *bisicles_vector, AmrIceHolderClass *c_AmrIceHolderPtr, int level_id)
{
    auto start_time = std::chrono::high_resolution_clock::now();  // Start time

    const Vector<LevelData<FArrayBox>*>& test = bisicles_vector->GetHVector();
    int num_lvl = test.size();
    vector<int> num_boxes_per_lvl_per_rank(num_lvl);
    vector<int> num_cells_per_lvl(num_lvl);
    vector<int> num_cells_per_lvl_per_rank(num_lvl);
    vector<int> num_cells_per_box(num_lvl);
    int num_total_cells_per_rank = 0;

    // Pre-calculate some numbers of cells and boxes
    vector<int> num_boxes_per_lvl(num_lvl);
    for (int lvl = 0; lvl < num_lvl; ++lvl) {
        LevelData<FArrayBox>& ldf = *test[lvl];
        DisjointBoxLayout dbl = ldf.disjointBoxLayout();
        num_cells_per_lvl[lvl] = dbl.numCells();
        num_boxes_per_lvl[lvl] = dbl.size();
        num_cells_per_box[lvl] = num_cells_per_lvl[lvl] / num_boxes_per_lvl[lvl];
        num_boxes_per_lvl_per_rank[lvl] = 0;
        num_cells_per_lvl_per_rank[lvl] = 0;
    }

    // Main loop
    for (int lvl = 0; lvl < num_lvl; ++lvl) {
        LevelData<FArrayBox>& ldf = *test[lvl];
        DataIterator dit = ldf.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            num_boxes_per_lvl_per_rank[lvl]++;
            num_cells_per_lvl_per_rank[lvl] += num_cells_per_box[lvl];
        }
        num_total_cells_per_rank += num_cells_per_lvl_per_rank[lvl];
    }

    Real** flattened_box;
    vector<Real*> flattened_array(num_lvl);
    for (int lvl = 0; lvl < num_lvl; ++lvl) {
        int num_cells = num_boxes_per_lvl_per_rank[lvl] * num_cells_per_box[lvl];
        flattened_array[lvl] = (Real*)malloc(num_cells * sizeof(Real));

        LevelData<FArrayBox>& ldf = *test[lvl];
        DisjointBoxLayout dbl = ldf.disjointBoxLayout();
        DataIterator dit = ldf.dataIterator();
        int box_index = 0;
        for (dit.reset(); dit.ok(); ++dit) {
            const Box& box = dbl[dit()];
            FArrayBox& fab = ldf[dit()];
            Real* flattened_data = fab.flattenAll(box);
            memcpy(flattened_array[lvl] + box_index * num_cells_per_box[lvl], flattened_data, num_cells_per_box[lvl] * sizeof(Real));
            free(flattened_data);
            box_index++;
        }
    }

    // Concatenate all level arrays
    Real* flattened_level = (Real*)malloc(num_total_cells_per_rank * sizeof(Real));
    int offset = 0;
    for (int lvl = 0; lvl < num_lvl; ++lvl) {
        memcpy(flattened_level + offset, flattened_array[lvl], num_boxes_per_lvl_per_rank[lvl] * num_cells_per_box[lvl] * sizeof(Real));
        offset += num_boxes_per_lvl_per_rank[lvl] * num_cells_per_box[lvl];
        free(flattened_array[lvl]); // Free each level array after copying
    }
    // check if flattened array corretcly copied all flattened box values
      //   cout<<"packing flattened_level "<<flattened_level<<endl;
      //   printf("flattened_level of size %lld:",num_total_cells_per_rank,"\n");
      //   for (int i = 0; i < num_total_cells_per_rank; i++) {
      //         printf("%lf ", flattened_level[i]);
      //   }

    // Update the timing and call count in AmrIceHolderClass
    auto end_time = std::chrono::high_resolution_clock::now();  // End time
    std::chrono::duration<double> duration = end_time - start_time;
    c_AmrIceHolderPtr->updateCallTimePack(duration);

   //  std::cout << "BisiclesVector Pack called " << call_count << " times. "
   //            << "Total time: " << total_duration.count() << " seconds." << std::endl;

    return flattened_level;
}





void BisiclesVectorUnpack(BisiclesVector *bisicles_vector, Real* flattened_level,
                          AmrIceHolderClass *c_AmrIceHolderPtr)
{
    auto start_time = std::chrono::high_resolution_clock::now();  // Start time

    Vector<LevelData<FArrayBox>* > HVector = bisicles_vector->GetHVector();
    int num_lvl = HVector.size();
    vector<int> num_boxes_per_lvl_per_rank(num_lvl);
    vector<int> num_cells_per_lvl(num_lvl);
    vector<int> num_cells_per_lvl_per_rank(num_lvl);
    vector<int> num_cells_per_box(num_lvl);
    int num_total_cells_per_rank = 0;

    // Calculate the number of cells and boxes per level
    for (int lvl=0; lvl < num_lvl; lvl++)
    {
        LevelData<FArrayBox>& ldf = *HVector[lvl];
        DisjointBoxLayout dbl = ldf.disjointBoxLayout();
        DataIterator dit = ldf.dataIterator();
        num_cells_per_lvl[lvl] = dbl.numCells();
        int num_boxes_per_lvl = dbl.size();
        num_cells_per_box[lvl] = num_cells_per_lvl[lvl] / num_boxes_per_lvl;
        num_boxes_per_lvl_per_rank[lvl] = 0;
        num_cells_per_lvl_per_rank[lvl] = 0;

        for (dit.reset(); dit.ok(); ++dit)
        {
            num_boxes_per_lvl_per_rank[lvl]++;
            num_cells_per_lvl_per_rank[lvl] += num_cells_per_box[lvl];
        }
        num_total_cells_per_rank += num_cells_per_lvl_per_rank[lvl];
    }

    int offset = 0;
    for (int lvl = 0; lvl < num_lvl; lvl++)
    {
        LevelData<FArrayBox>& ldf = *HVector[lvl];
        DataIterator dit = ldf.dataIterator();
        DisjointBoxLayout dbl = ldf.disjointBoxLayout();
        // Initialize storage for the boxes' data
        vector<Real*> flattened_box(num_boxes_per_lvl_per_rank[lvl]);
        int box_index=0;

        for (dit.reset(); dit.ok(); ++dit)
        {
            // Allocate and copy data from the flattened level to each box
            flattened_box[box_index] = (Real*)malloc(num_cells_per_box[lvl] * sizeof(Real));
            memcpy(flattened_box[box_index], flattened_level + offset, num_cells_per_box[lvl] * sizeof(Real));
            offset += num_cells_per_box[lvl];

            // Unpack the data back into the original FArrayBox
            const Box& box = dbl[dit()];
            FArrayBox& fab = ldf[dit()];
            fab.unpackAll(box, flattened_box[box_index]);
            box_index++;

            // Free the memory for this box after unpacking
            // free(flattened_box[box_index]);
        }
      //   free(flattened_box);
    }
    auto end_time = std::chrono::high_resolution_clock::now();  // End time
    std::chrono::duration<double> duration = end_time - start_time;
    c_AmrIceHolderPtr->updateCallTimeUnpack(duration);

   //  std::cout << "BisiclesVector Unpack called " << call_count << " times. "
   //            << "Total time: " << total_duration.count() << " seconds." << std::endl;
}




   // needs to be FIXED!!!!
   double BisiclesVectorNorm(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      return bisicles_vector->HL2Norm();
   }

   void BisiclesVectorAxpy(BisiclesVector *y, double a, BisiclesVector *x, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      y->HAxpy(a,x,c_AmrIceHolderPtr);
   }

   void BisiclesVectorL2Print(BisiclesVector *bisicles_vector)
   {
      bisicles_vector->PrintHL2norm();
   }

   void BisiclesVectorLevelDataBox(BisiclesVector *bisicles_vector, BisiclesVector *x)
   {
      const Vector<LevelData<FArrayBox>* >& test=bisicles_vector->GetHVector();
      pout()<<"\n--------------- EL - begin grid info ------------------\n";
      bisicles_vector->PFPrintLevelData(x);
      pout()<<"\n--------------- EL - end grid info ------------------\n";
   }

   double BisiclesVectorGetVal(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      return bisicles_vector->DataGet();
   }


   void PfasstBisiclesSaveResults(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      Vector<LevelData<FArrayBox>* > HVector_print=bisicles_vector->GetHVector();
      bisicles_vector->SaveSnapshot(c_AmrIceHolderPtr);
   }



   // bisicles solver functions
   void BisiclesSolverInit(BisiclesdHdtSolver **bisicles_dHdt, int pfasst_level_index, \
                           int num_grid_points,AmrIceHolderClass *c_AmrIceHolderPtr)
    {
            *bisicles_dHdt = new BisiclesdHdtSolver(num_grid_points,c_AmrIceHolderPtr);
    }

    void BisiclesVectorSetHIC(BisiclesVector *bisicles_vector,AmrIceHolderClass *c_AmrIceHolderPtr)
    {
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      cout<< "bisicles_fortran.cpp done assign H IC size of H "<<constH.size()<< std::endl;
    }


void BisiclesSolverFEval(BisiclesdHdtSolver *bisicles_dHdt, BisiclesVector *y, double t,\
   int pfasst_level_index, BisiclesVector *f, double dt,int maxStep,bool evolve_velocity, AmrIceHolderClass *c_AmrIceHolderPtr)

   {
      AmrIce *amrObjHolderPtr;
      amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
      IceSheetState iceState;
      iceState=y->GetIceState();

      Vector<LevelData<FArrayBox>* > dHdtVect;
      Vector<LevelData<FArrayBox>* > H=y->GetHVector(); // this H is the updating one from pfasst
      reshapeAndFill(iceState.ice_thickness, H); // H is from y vector in pfasst, should not change this line 
      
      // time setState function
      auto start_setState = std::chrono::high_resolution_clock::now();
      amrObjHolderPtr->setState(iceState,evolve_velocity); // constH is updated after setState
      auto end_setState = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration_setState = end_setState - start_setState;
      
      reshapeAndFill(dHdtVect,H);

      // time computedHdt function
      auto start_compute_dHdt = std::chrono::high_resolution_clock::now();
      amrObjHolderPtr->compute_dHdt(dHdtVect, H, dt, evolve_velocity);
      auto end_compute_dHdt = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration_compute_dHdt = end_compute_dHdt - start_compute_dHdt;// cout<<"dHdtVect norm 2: ";
      
      f->SetdHdtVector(dHdtVect,c_AmrIceHolderPtr);
      
      c_AmrIceHolderPtr->updateCallTimeFEval(duration_setState,duration_compute_dHdt);
   }


   void PfasstPrintAmr(BisiclesVector *y, AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      
      int nlvl=constH.size();
      for (int lvl=0; lvl < nlvl; lvl++)
      {
      LevelData<FArrayBox>& ldf = *constH[lvl];
      DisjointBoxLayout dbl = ldf.disjointBoxLayout();
      DataIterator dit = ldf.dataIterator();
      cout<<"in pfasst level "<<lvl<<endl;
      for (dit.reset(); dit.ok(); ++dit) 
         {
         const Box& box = dbl[dit()];
         FArrayBox& fab = ldf[dit()];
         cout<<"  box "<<box<<endl;
         } 
      }
   }

   MPI_Fint PfasstBisicles_MPI_Comm_c2f(MPI_Comm *comm) {
      return MPI_Comm_c2f(*comm);
   } 



}