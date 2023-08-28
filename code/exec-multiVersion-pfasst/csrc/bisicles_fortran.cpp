#include "bisicles_fortran.hpp"
#include "bisicles_vector.hpp"
#include "bisicles_holder.hpp"
#include "AmrIce.H"


#include "LevelData.H"
#include "FArrayBox.H"
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
#include "AmrIceBase.H"
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
using namespace std;


MPI_Comm glob_space_comm = MPI_COMM_NULL;

extern "C"
{

// bisicles vector functions
void BisiclesVectorCreate(BisiclesVector **bisicles_vector,int num_grid_points, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      *bisicles_vector = new BisiclesVector(num_grid_points,c_AmrIceHolderPtr);
      // Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
   //    Vector<LevelData<FArrayBox>* > temp=bisicles_vector->GetHVector();
   //    int nlvl=temp.size();
   // for (int lvl=0; lvl < nlvl; lvl++)
   // {
   //   LevelData<FArrayBox>& ldf = *temp[lvl];
   //   DisjointBoxLayout dbl = ldf.disjointBoxLayout();
   //   DataIterator dit = ldf.dataIterator();
   //   cout<<"bisicle vector create level "<<lvl<<endl;
   //   for (dit.reset(); dit.ok(); ++dit) 
   //    {
   //     const Box& box = dbl[dit()];
   //     FArrayBox& fab = ldf[dit()];
   //    //  test_min=fab.norm(box,1);
   //    cout<<"  box "<<box<<endl;
   //    } 
   // }

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
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      // cout<<"-- 2 bis_for:\n";
      // bisicles_vector->PrintLevelData(constH);

   }

   void BisiclesVectorCopy(BisiclesVector *dest, BisiclesVector *src, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      dest->SetHVector(src);
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      // cout<<"-- 1 bis_for:\n";
      // src->PrintLevelData(constH);
   }

Vector<LevelData<FArrayBox>* >* BisiclesVectorPack(BisiclesVector *bisicles_vector,int num_grid_points, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      cout<< "bisicles_fortran.cpp packing............................"<< std::endl;
      // return bisicles_vector->GetHDataPtr(c_AmrIceHolderPtr);
      return bisicles_vector->GetHDataPtr(c_AmrIceHolderPtr);
   }

   void BisiclesVectorUnpack(BisiclesVector *bisicles_vector, double y, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      cout<< "bisicles_fortran.cpp unpacking............................"<< std::endl;
      bisicles_vector->HSetVal2All(y,c_AmrIceHolderPtr);
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
      bisicles_vector->PFPrintLevelData(x);
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
      cout<< "-------------------------------------- done saving pfasst_bisicles results ----------------------------------------- "<< std::endl;
   }



   // bisicles solver functions
   void BisiclesSolverInit(BisiclesdHdtSolver **bisicles_dHdt, int pfasst_level_index, \
                           int num_grid_points,AmrIceHolderClass *c_AmrIceHolderPtr)
    {
      MPI_Comm newcomm;
      newcomm = glob_space_comm;
      // carefull about the pfasst_level_index, check if level_index=num_of_total_levels-pfasst_level_index
      int level_index = pfasst_level_index;
      *bisicles_dHdt = new BisiclesdHdtSolver(newcomm,num_grid_points,c_AmrIceHolderPtr);
      
//       Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
//       int nlvl=constH.size();
//    for (int lvl=0; lvl < nlvl; lvl++)
//    {
//      LevelData<FArrayBox>& ldf = *constH[lvl];
//      DisjointBoxLayout dbl = ldf.disjointBoxLayout();
//      DataIterator dit = ldf.dataIterator();
//      cout<<"bisicle init level "<<lvl<<endl;
//      for (dit.reset(); dit.ok(); ++dit) 
//       {
//        const Box& box = dbl[dit()];
//        FArrayBox& fab = ldf[dit()];
//       //  test_min=fab.norm(box,1);
//       cout<<"  box "<<box<<endl;
//       } 
//   }
    }

    void BisiclesVectorSetHIC(BisiclesVector *bisicles_vector,AmrIceHolderClass *c_AmrIceHolderPtr)
    {
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      cout<< "bisicles_fortran.cpp done assign H IC size of H "<<constH.size()<< std::endl;
    }


void BisiclesSolverFEval(BisiclesdHdtSolver *bisicles_dHdt, BisiclesVector *y, double t,\
   int pfasst_level_index, BisiclesVector *f, double dt,int maxStep, AmrIceHolderClass *c_AmrIceHolderPtr)

   {
      // AmrIce *amrObjHolderPtr;
      // amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
      // IceSheetState *iceStatePtr;
      // IceSheetState iceState;
      // iceStatePtr=y->GetIceStatePtr();
      // iceState=y->GetIceState();
   

      // // Vector<LevelData<FArrayBox>*> H_old = iceState.ice_thickness;
      // Vector<LevelData<FArrayBox>* > dHdtVect = c_AmrIceHolderPtr->GetAmrdHdt();
      // Vector<LevelData<FArrayBox>* > *dHdtVectPtr=c_AmrIceHolderPtr->GetAmrdHdtPtr();
      // Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH(); // double check with Dan
      // // cout<<"get amr H right from amr holder pointer\n";
      // // y->PrintLevelData(constH);

      // Vector<LevelData<FArrayBox>* > H=y->GetHVector(); // this H is the same as H_old
      // Vector<LevelData<FArrayBox>* > velo = amrObjHolderPtr->amrVelocity();

      // reshapeAndFill(iceState.ice_thickness, H); // H
      // // cout<<"H after reshape and fill ice state\n";
      // // y->PrintLevelData(constH);

      // bool recalculateVelocity = true;
      // // double check with Dan: dHdtVect & H new every time step?
      // // setState(IceSheetState& a_iceState... => 
      // amrObjHolderPtr->setState(iceState,recalculateVelocity); 
      // // cout<<"dHdt after set state\n";
      // // y->PrintLevelData(dHdtVect);
      // // cout<<"H after set state\n";
      // // y->PrintLevelData(constH);
      // reshape(dHdtVect,constH); // reshape the shape of dHdtVect into constH, function defined in SundialsUtil.cpp

      // amrObjHolderPtr->compute_dHdt(dHdtVect,constH,dt, recalculateVelocity); // dHdtVect is updated after this

      // f->SetdHdtVector(dHdtVect,c_AmrIceHolderPtr);




      AmrIce *amrObjHolderPtr;
      amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
      IceSheetState *iceStatePtr;
      IceSheetState iceState;
      iceStatePtr=y->GetIceStatePtr();
      iceState=y->GetIceState();

      Vector<LevelData<FArrayBox>*> H_old = iceState.ice_thickness;
      Vector<LevelData<FArrayBox>* > dHdtVect;
      Vector<LevelData<FArrayBox>* > *dHdtVectPtr=c_AmrIceHolderPtr->GetAmrdHdtPtr();
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH(); // double check with Dan
      Vector<LevelData<FArrayBox>* > constH_copy=c_AmrIceHolderPtr->GetAmrHBackup(); // 
      Vector<LevelData<FArrayBox>* > H=y->GetHVector(); // this H is the updating one from pfasst
      Vector<LevelData<FArrayBox>* > velo = amrObjHolderPtr->amrVelocity();

      // cout<<"H norm 2: ";
      // y->PrintL2norm(H);
      // cout<<"constH norm 2: ";
      // y->PrintL2norm(constH);

      reshapeAndFill(iceState.ice_thickness, H); // H is from y vector in pfasst, should not change this line 
      // reshapeAndFill(iceState.ice_velocity, velo); 
      // cout<<constH<<endl;
      bool recalculateVelocity = false;
      // cout<<"printing all entries in H \n";
      // y->PrintAllEntries(H);
      amrObjHolderPtr->setState(iceState,recalculateVelocity); // constH is updated after setState
      // cout<<"constH norm 2 after set state: ";
      // y->PrintL2norm(constH);
      
      // reshape(dHdtVect,constH); // this shouldn't matter, whether constH or constH_copy, whether reshape or reshapefill
      reshapeAndFill(dHdtVect,H);

      amrObjHolderPtr->compute_dHdt(dHdtVect,H,dt,recalculateVelocity); // dHdtVect is updated after this, H is used in the function be careful
      // cout<<"dHdtVect norm 2: ";
      // y->PrintL2norm(dHdtVect);
      f->SetdHdtVector(dHdtVect,c_AmrIceHolderPtr);
      // MayDay::Error("....exit here");
   }


   void PfasstPrintAmr(BisiclesVector *y, AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      // y->PrintLevelData(constH);
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
         //  test_min=fab.norm(box,1);
         cout<<"  box "<<box<<endl;
         } 
      }
   }

   MPI_Fint PfasstBisicles_MPI_Comm_c2f(MPI_Comm *comm) {
      cout<<"in mpi c2f: "<<comm<<endl;
      return MPI_Comm_c2f(*comm);
   } 



}