#include "bisicles_fortran.hpp"
// #include "cencap.hpp"
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
using namespace std;


MPI_Comm glob_space_comm = MPI_COMM_NULL;

extern "C"
{

// bisicles vector functions
void BisiclesVectorCreate(BisiclesVector **bisicles_vector,int num_grid_points, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      *bisicles_vector = new BisiclesVector(num_grid_points,c_AmrIceHolderPtr);
      // double intialization_factor = -1;
      // bisicles_vector->HAxpy(intialization_factor,bisicles_vector,c_AmrIceHolderPtr);
      // c_AmrIceHolderPtr=0x7ff7bfefb7c8
      // cout<< "bisicles_fortran.cpp vector init c_AmrIceHolderPtr "<<c_AmrIceHolderPtr<< std::endl;
      // cout<< "5555 bisicles_fortran.cpp 000 created bisicles_vector "<<bisicles_vector<< std::endl;
      // Vector<LevelData<FArrayBox>* > *crsedHdt=(*bisicles_vector)->GetVectorPtr();
      // cout<< "bisicles_fortran.cpp 666 get vector ptr crsedHdt "<<crsedHdt<< std::endl;

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
      // bisicles_vector->DataSetVal(y); double y
      // Vector<LevelData<FArrayBox>* > *AmrdHdtPtr=c_AmrIceHolderPtr->GetAmrdHdtPtr();
      // bisicles_vector->dHdtSetVal2All(y,c_AmrIceHolderPtr);
      bisicles_vector->HSetVal2All(y,c_AmrIceHolderPtr);
      // Vector<LevelData<FArrayBox>* > updated_vector=bisicles_vector->dHdtVector;

   }

   void BisiclesVectorCopy(BisiclesVector *dest, BisiclesVector *src, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      // dest->DataSetVal(src->DataGet());
      // dest->SetdHdtVector(src->GetdHdtVector());
      dest->SetHVector(src);
   }

   double *BisiclesVectorPack(BisiclesVector *bisicles_vector,int num_grid_points, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      // return bisicles_vector->DataGet();
      // return bisicles_vector->GetdHdtDataPtr();
      cout<< "bisicles_fortran.cpp packing............................"<< std::endl;
      return bisicles_vector->GetHDataPtr(c_AmrIceHolderPtr);
   }

   void BisiclesVectorUnpack(BisiclesVector *bisicles_vector, double y, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      // bisicles_vector->DataSetVal(y);
      // double check with Dan
      // bisicles_vector->dHdtSetVal2All(y,c_AmrIceHolderPtr);
      cout<< "bisicles_fortran.cpp unpacking............................"<< std::endl;
      bisicles_vector->HSetVal2All(y,c_AmrIceHolderPtr);
   }

   // needs to be FIXED!!!!
   double BisiclesVectorNorm(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      // return bisicles_vector->DataNorm();
      // return bisicles_vector->dHdtL2Norm();
      return bisicles_vector->HL2Norm();
   }

   void BisiclesVectorAxpy(BisiclesVector *y, double a, BisiclesVector *x, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      // y->DataAxpy(a, x->DataGet());
      // y->dHdtAxpy(a,x,c_AmrIceHolderPtr);
      // cout<< "bisicles_fortran.cpp  before axpy y "<<y<<" x "<<x<< std::endl;
      y->HAxpy(a,x,c_AmrIceHolderPtr);
   }

   void BisiclesVectorPrint(BisiclesVector *bisicles_vector)
   {
      // bisicles_vector->DataPrint();
      bisicles_vector->PrintHL2norm();
   }

   double BisiclesVectorGetVal(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      return bisicles_vector->DataGet();
   }


   void PfasstBisiclesSaveResults(BisiclesVector *bisicles_vector, \
                          AmrIceHolderClass *c_AmrIceHolderPtr)
   {
      cout<< "bisicles_fortran.cpp  before save snapshot bisicles_vectors ID "<<bisicles_vector<< std::endl;
      Vector<LevelData<FArrayBox>* > HVector_print=bisicles_vector->GetHVector();
      cout<< "bisicles_fortran.cpp  before save snapshot HVector size "<<HVector_print.size()<< std::endl;
      bisicles_vector->SaveSnapshot(c_AmrIceHolderPtr);
      // AmrIce *amrObjHolderPtr;
      // IceSheetState iceState;

      
      // amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
      // iceState=bisicles_vector->GetIceState();
      // Vector<LevelData<FArrayBox>* > H=bisicles_vector->GetHVector(); 

      // reshapeAndFill(iceState.ice_thickness, H);

      // c_AmrIceHolderPtr->SetAmrIceState(iceState);

      // // amrObjHolder=c_AmrIceHolderPtr->GetAmrIceObj();

      // amrObjHolderPtr->writePlotFile();
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
      // cout<< "bisicles_fortran.cpp solver init c_AmrIceHolderPtr "<<c_AmrIceHolderPtr<< std::endl;
    }

    void BisiclesVectorSetHIC(BisiclesVector *bisicles_vector,AmrIceHolderClass *c_AmrIceHolderPtr)
    {
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      cout<< "bisicles_fortran.cpp done assign H IC size of H "<<constH.size()<< std::endl;
      // bisicles_vector->SetHVector(constH);
    }













void BisiclesSolverFEval(BisiclesdHdtSolver *bisicles_dHdt, BisiclesVector *y, double t,\
   int pfasst_level_index, BisiclesVector *f, double dt,int maxStep, AmrIceHolderClass *c_AmrIceHolderPtr)
// void BisiclesSolverFEval(BisiclesdHdtSolver *bisicles_dHdt, HypreVector *y, double t, int pfasst_level_index, HypreVector *f)

   {
      // int level_index = PfasstToHypreLevelIndex(pfasst_level_index, hypre_solver->GetNumLevels());

      // int nrows = hypre_solver->GetNumRows();
      // double *f_values = (double *)malloc(nrows * sizeof(double));
      // if (piece == 2){
      //    hypre_solver->FEval(y->GetBoxValues(), t, level_index, &f_values);
      // }
      // else {
      //    for (int i = 0; i < nrows; i++){
      //       f_values[i] = 0.0;
      //    }
      // }
      // f->SetBoxValues(f_values);


    // f -> rhs, y -> lhs, dy/dt=f
    // size of *f=72, f=8, *bisicles_dHdt=88, *y=72
    // type of f=y=P14BisiclesVector, bisicles_dHdt=P18BisiclesdHdtSolver

   	 //int nrows = sizeof(*y);
       // cout << "bisicles_fortran.cpp 0000 row of bisicles_dHdt " << sizeof(*y) << std::endl;
       // cout << "bisicles_fortran.cpp 0000 col of y " << sizeof(*y[0])/sizeof(*y[0][0]) << std::endl;
       // int num_of_rows = bisicles_dHdt->GetNumRows();
       // Vector<LevelData<FArrayBox>* > crsedHdt=y->GetdHdtVector();
       // type of crsedHdt=6VectorIP9LevelDataI9FArrayBoxEE
       // cout<< "bisicles_fortran.cpp 666 &crsedHdt "<<&crsedHdt<< std::endl;
       // int nrows = sizeof(*bisicles_dHdt);
   	 // double *f_values = (double *)malloc(nrows * sizeof(double));
   	 // double *dHdt_values = (double *)malloc(nrows * sizeof(double));

       // call amrice obj, state vec from bisicles
       // AmrIceHolderClass amrHolder;
       // amrPtr=bisicles_dHdt->GetAmrIceState();
       // try create a new amrice obj and call functions
       // AmrIce *amrPtr;
       // BisiclesdHdtSolver amrPtr;
       // amrPtr->setParmParsePrefix("crse.");
       // cout<< "bisicles_fortran.cpp 000 c_AmrIceHolderPtr "<<c_AmrIceHolderPtr<< std::endl;
       // cout<< "bisicles_fortran.cpp 000 type of c_AmrIceHolderPtr "<<typeid(c_AmrIceHolderPtr).name()<< std::endl;
       // cout<< "bisicles_fortran.cpp 000 c_AmrIceHolderPtr "<<amrHolderPtr<< std::endl;
       AmrIce *amrObjHolderPtr;
       amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
       // double amr_dt=amrObjHolderPtr->dt();
       // const Vector<Real>& amr_dx=amrObjHolderPtr->amrDx();
       // int amr_finest_level=amrObjHolderPtr->finestLevel();
       // const Vector<DisjointBoxLayout>& grid_size=amrObjHolderPtr->grids();
       // cout<< "from f_eval: amr_dt "<<amr_dt<< std::endl; // dt=1e20
       // cout<< "from f_eval: amr_finest_level "<<amr_finest_level<< std::endl; // finest level=0
       // cout<< "from f_eval: grid_size "<<grid_size[0]<< std::endl;
      // typedef std::vector<Real> Amr_dx;
      // for (Amr_dx::size_type i=0; i<amr_dx.size(); ++i){
      // for (int i=0;i<amr_dx.size(); ++i){
      //   std::cout<<"i "<<i<<": " << *&amr_dx[i] << std::endl;
      // }

      IceSheetState *iceStatePtr;
      IceSheetState iceState;
      iceStatePtr=y->GetIceStatePtr();
      iceState=y->GetIceState();


      Vector<LevelData<FArrayBox>*> H_old = iceState.ice_thickness;
      // iceState=*amrStateHolderPtr;
      // IceSheetState icestatetemp=*amrStateHolderPtr;
      // Vector<LevelData<FArrayBox>* > ice_thick=iceState.ice_thickness;
      // cout<< "from f_eval: ice_thick "<<ice_thick<< std::endl; 

      Vector<LevelData<FArrayBox>* > dHdtVect;
      // cout<< "from f_eval: H vector "<<H<< std::endl; 

      Vector<LevelData<FArrayBox>* > *dHdtVectPtr=c_AmrIceHolderPtr->GetAmrdHdtPtr();
      Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH(); // double check with Dan
      Vector<LevelData<FArrayBox>* > H=y->GetHVector(); // this H is the same as H_old
      Vector<LevelData<FArrayBox>* > rhs;

      reshapeAndFill(iceState.ice_thickness, H); // H is read in correctly. verified
      
      // bool recalculateVelocity = false;
      amrObjHolderPtr->setState(iceState,false); // double check with Dan: dHdtVect & H new every time step?
      reshape(dHdtVect,constH);

      // Real test_min=0;
      // int nlvl=H.size();
      // for (int lvl=0; lvl < nlvl; lvl++)
      // {
      //   LevelData<FArrayBox>& ldf = *H[lvl];
      //   DisjointBoxLayout dbl = ldf.disjointBoxLayout();
      //   DataIterator dit = ldf.dataIterator();
      //   for (dit.reset(); dit.ok(); ++dit) 
      //   {
      //     const Box& box = dbl[dit()];
      //     FArrayBox& fab = ldf[dit()];
      //     test_min=fab.norm(box,1);
      //     cout<< "from bisicles_fortran: before f_eval norm of H "<<test_min << std::endl;
      //   }
      // }

      amrObjHolderPtr->compute_dHdt(dHdtVect,constH,dt, true); // dHdtVect is updated after this

      // for (int lvl=0; lvl < nlvl; lvl++)
      // {
      //   LevelData<FArrayBox>& ldf = *dHdtVect[lvl];
      //   DisjointBoxLayout dbl = ldf.disjointBoxLayout();
      //   DataIterator dit = ldf.dataIterator();
      //   for (dit.reset(); dit.ok(); ++dit) 
      //   {
      //     const Box& box = dbl[dit()];
      //     FArrayBox& fab = ldf[dit()];
      //     test_min=fab.norm(box,1);
      //     cout<< "from bisicles_fortran: after f_eval norm of dHdt "<<test_min << std::endl;
      //   }
      // }

      // LevelData<FArrayBox>& ldf = *dHdtVect[0];

      // cout<< "from f_eval: after compute H min of array box "<<test_min<< std::endl;
      // cout<< "from f_eval: after compute dHdt "<<dHdtVect[0]<< std::endl;
      // amrObjHolderPtr->run(iceState.time, maxStep);
      // amrObjHolderPtr->getState(iceState);


       // bisicles_dHdt->FEval(t, pfasst_level_index, &rhs, &dHdtVect);
       f->SetdHdtVector(dHdtVect,c_AmrIceHolderPtr);

      // MayDay::Abort("------------------------------- stop here --------------------------------------");

       // free(dHdtVect);
       // f->SetIceState();
       // bisicles_dHdt->FEval(y->DataPtrGet(), t, pfasst_level_index, &f_values);
       // f->DataPtrSetVal(f_values);
       // AmrIceHolderPtr.SetAmrIceState(amrObjectCrsePtr);
       // bisicles_dHdt->computedHdt(amrObjectCrse->DataPtrGet(),crsedHdtVect,crseH);

       // cout<< "from f_eval: dt "<<dt<< std::endl;

   	 // for (int i=0; i<nrows; i++){
   	 // 	//(*dHdt_values)[i] = crsedHdtVect[i];
   	 // 	// f_values[i]=amrObjectCrse.compute_dHdt(crsedHdtVect[i],crseH,crseDt, false);
   	 // }


   }



}