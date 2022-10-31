#include "bisicles_vector.hpp"
//#include <iostream>
#include <cstdlib>


#include "BISICLES_VERSION.H"
#include "Box.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "parstream.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "FineInterp.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "AmrIce.H"
#include "computeNorm.H" 
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "PiecewiseLinearFillPatch.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "DerivativesF_F.H"
#include "DivergenceF_F.H"
#include "computeSum.H"
#include "CONSTANTS.H"
#include "IceConstants.H"
#include "ExtrapBCF_F.H"
#include "amrIceF_F.H"
#include "BisiclesF_F.H"
#include "IceThermodynamics.H"
#include "JFNKSolver.H"
#include "InverseVerticallyIntegratedVelocitySolver.H"
#include "PetscIceSolver.H"
#include "RelaxSolver.H"
#ifdef CH_USE_FAS
#include "FASIceSolver.H"
#endif
#include "KnownVelocitySolver.H"
#include "VCAMRPoissonOp2.H"
#include "AMRPoissonOpF_F.H"
#include "CH_HDF5.H"
#include "IceUtility.H"
#include "LevelMappedDerivatives.H"


// ------- vector create should be good now ------- //
BisiclesVector::BisiclesVector(int num_grid_points, AmrIceHolderClass *c_AmrIceHolderPtr)
{
   AmrIce *amrObjHolderPtr;
   amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();

   // dHdtVector=c_AmrIceHolderPtr->GetAmrdHdt();
   // dHdtVectorPtr=c_AmrIceHolderPtr->GetAmrdHdtPtr();

   Vector<LevelData<FArrayBox>* > H_ref=c_AmrIceHolderPtr->GetAmrH();
   // HVectorPtr=c_AmrIceHolderPtr->GetAmrHPtr();


   // iceState=c_AmrIceHolderPtr->GetAmrIceState();
   // iceStatePtr=c_AmrIceHolderPtr->GetAmrIceStatePtr();
   
   // LevelData<FArrayBox>& ldf = *HVector[0];
      // DataIterator dit = ldf.dataIterator();
      // DisjointBoxLayout dbl = ldf.disjointBoxLayout();
      // Real test_min=0;
      // for (dit.reset(); dit.ok(); ++dit) {
      //    FArrayBox& fab = ldf[dit()];
      //    const Box& box = dbl[dit()];
      //    test_min=fab.norm(box);
      // }
   // cout<< "!!!!!!!!!!!!!!!!:0000 norm of H "<<test_min<< std::endl;
   
   int num_levels_bisicles_vec=c_AmrIceHolderPtr->GetAmrNumLvl();
   refineRatio=amrObjHolderPtr->refRatios();
   amrDx=amrObjHolderPtr->amrDx();
   Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();

   Vector<LevelData<FArrayBox>* > HVector_ref(num_levels_bisicles_vec, NULL);
   HVector=HVector_ref;
   Vector<LevelData<FArrayBox>* > dHdtVector_ref(num_levels_bisicles_vec, NULL);
   dHdtVector=dHdtVector_ref;
   for (int lvl=0; lvl<num_levels_bisicles_vec;lvl++)
   {
      // Vector<DisjointBoxLayout> current_grid = H_ref
      const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl);
      HVector[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero); // FIXED the num of component=1, i.e. one d?? only Vx no Vy??
      dHdtVector[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero);
      refineRatio[lvl]=2; // double check with Dan
      // cout<< "!!!!!!!!!!!!!!!!:0000 refine ratio "<<refineRatio[lvl]<<" dx "<<amrDx[lvl]<< std::endl;
      // default value for C0 is, fittingly enough, 0
        
     DataIterator dit = (*HVector[lvl]).dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
         (*HVector[lvl])[dit].setVal(0.0);
         (*dHdtVector[lvl])[dit].setVal(0.0);
      }
        
   }

   reshapeAndFill(iceState.ice_thickness, HVector);
   // reshapeAndFill(HVector,constH);
   // reshapeAndFill(dHdtVector,constH);

   // cout<< "!!!!!!!!!!!!!!!!:0000 num of levels "<<HVector.size()<< std::endl;


   InitGrid(num_grid_points);

}

BisiclesVector::~BisiclesVector(void)
{

}

Vector<LevelData<FArrayBox>* > BisiclesVector::GetdHdtVector(void)
{
   return dHdtVector;
}

Vector<LevelData<FArrayBox>* > *BisiclesVector::GetdHdtVectorPtr(void)
{
   return dHdtVectorPtr;
}

Vector<LevelData<FArrayBox>* > BisiclesVector::GetHVector(void)
{
   return HVector;
}

Vector<LevelData<FArrayBox>* > *BisiclesVector::GetHVectorPtr(void)
{
   return HVectorPtr;
}

IceSheetState BisiclesVector::GetIceState(void)
{
   return iceState;
}

IceSheetState *BisiclesVector::GetIceStatePtr(void)
{
   return iceStatePtr;
}

// double check with Hans/Dan
void BisiclesVector::SetdHdtVector(Vector<LevelData<FArrayBox>* > src,AmrIceHolderClass *c_AmrIceHolderPtr)
{
  
   Real test_min=0;
   int nlvl=HVector.size();
   double val=0;
   double a=1;
   Vector<LevelData<FArrayBox>* > xdHdtVector = src; // 

   for (int lvl=0; lvl < nlvl; lvl++)
   {
      // dHdtVector[lvl]=xdHdtVector[lvl];
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     LevelData<FArrayBox>& xldf = *xdHdtVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       const FArrayBox& xfab = xldf[dit()];
       test_min=xfab.norm(box,1);
       // cout<< "from bisicles_vector: before setH norm of H "<<test_min <<" num of components fab "<<fab.nComp()<<" num of components xfab "<<xfab.nComp()<< std::endl;
 
       // fab.setVal(val);      
       // fab.plus(xfab,a);
       fab.copy(xfab,box);
       test_min=fab.norm(box,1);
      } 
  
   // cout<< "from bisicles_vector: before setH norm of H<-dHdt  "<<test_min<< std::endl;
  }
}

// // currently passing in only dHdt vector, and assign that to the ice state
// void BisiclesVector::SetIceState(Vector<LevelData<FArrayBox>* > vector_value)
// {
//    // cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ setvector dHdtVector size "<<dHdtVector.size()<< std::endl;

//    for(int i=0;i<vector_value.size();i++)
//    {
//       dHdtVector[i]=vector_value[i];
//       // cout<< "bisicles_vector.cpp i "<<i<<" vector_value "<<dHdtVector[i]<< std::endl;
//    } 

void BisiclesVector::SetHVector(BisiclesVector* src)
{
   Real test_min=0;
   int nlvl=HVector.size();
   double val=0;
   double a=1;
   Vector<LevelData<FArrayBox>* > xHVector = src->HVector; // 
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     LevelData<FArrayBox>& xldf = *xHVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       const FArrayBox& xfab = xldf[dit()];
       test_min=fab.norm(box,2);
       // cout<< "from bisicles_vector: before setH norm of fab "<<test_min << " a " << a<< std::endl;
 
       fab.setVal(val);      
       fab.plus(xfab,a);
       test_min=fab.norm(box,2);
      } 
  
   // cout<< "from bisicles_vector: before setH norm of fab  "<<test_min<< std::endl;
  }
}


// int BisiclesVector::GetNumRows(void)
// {
//    return nrows;
// }

double *BisiclesVector::DataPtrGet(void)
{
   double *data = (double *)calloc(num_grid_points, sizeof(double));
   return data;
}

double BisiclesVector::DataGet(void)
{
   return data;
}

void BisiclesVector::DataPtrSetVal(double *y)
{
   data = *y;
}

void BisiclesVector::DataSetVal(double y)
{
   data = y;
}

void BisiclesVector::dHdtSetVal2All(double val,AmrIceHolderClass *c_AmrIceHolderPtr)
{

   // currently the num of level is only 1, before assigning, norm of dHdt=144.889
   for (int lvl=0; lvl < dHdtVector.size(); lvl++)
   {
     DataIterator dit = (*dHdtVector[lvl]).dataIterator();
     // currently dit.ok is only 1
     for (dit.reset(); dit.ok(); ++dit) {
       (*dHdtVector[lvl])[dit()].setVal(val);
     }
   }
   // this function is only called for initilization purposes %setval(0.0_pfdp)?? double check with Jordi
   // there should other places use setVal
   // c_AmrIceHolderPtr->SetAmrdHdtPtr(&dHdtVector);
   // c_AmrIceHolderPtr->SetAmrdHdt(dHdtVector); // after assigning, norm dHdt=0
}

void BisiclesVector::HSetVal2All(double val,AmrIceHolderClass *c_AmrIceHolderPtr)
{
   // // double check with Dan: do we need to assign values of all zero??

   // LevelData<FArrayBox>& ldf = *HVector[0];
   //    DataIterator dit = ldf.dataIterator();
   //    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
   //    Real test_min=0;
   //    for (dit.reset(); dit.ok(); ++dit) {
   //       FArrayBox& fab = ldf[dit()];
   //       const Box& box = dbl[dit()];
   //       test_min=fab.norm(box);
   //    }
   for (int lvl=0; lvl < HVector.size(); lvl++)
   {
     DataIterator dit = (*HVector[lvl]).dataIterator();
     // currently dit.ok is only 1
     for (dit.reset(); dit.ok(); ++dit) {
       (*HVector[lvl])[dit()].setVal(val);
     }
   }
   // this function is only called for initialization purposes %setval(0.0_pfdp)?? double check with Jordi
   // there should other places use setVal
   // c_AmrIceHolderPtr->SetAmrHPtr(&HVector);
   // c_AmrIceHolderPtr->SetAmrH(HVector); // after assigning, norm dHdt=0
}

double *BisiclesVector::GetdHdtDataPtr(void)
{
   // double check with Dan to find the total num of entries
   double *data = (double *)calloc(num_grid_points, sizeof(double));
   return data;
}

double *BisiclesVector::GetHDataPtr(AmrIceHolderClass *c_AmrIceHolderPtr)
{
   // count the grid cells in each level and store into vecotr
   AmrIce *amrObjHolderPtr;
   amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
   int nlvl=HVector.size();
   Vector<int> num_of_grid_points_vec(nlvl);
   int total_num_grid_points=0;

   for (int lvl=0; lvl < nlvl; lvl++)
   {
     const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl);
     LayoutIterator lit = current_grid_size.layoutIterator();
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();

     // int count=0;
     for (lit.begin(); lit.ok(); ++lit)
     {
       const Box& thisBox = current_grid_size.get(lit());
       num_of_grid_points_vec[lvl] += thisBox.numPts();
       // cout<< "from bisicles_vector: num of grid cells "<<num_of_grid_points_vec[lvl] << std::endl;
       // count++;

       // FArrayBox& fab = ldf[dit()];
      }
      total_num_grid_points +=num_of_grid_points_vec[lvl];
   }
   // cout<< "from bisicles_vector: num of grid cells "<<total_num_grid_points << std::endl;
   
   // allocate the vector for H and store the pointers in value
   double *values = (double *)calloc(total_num_grid_points, sizeof(double));

}

double BisiclesVector::DataNorm(void)
{
   cout<< "bisicles_vector.cpp abs data "<<abs(data)<< std::endl;
   return abs(data);
}

double BisiclesVector::dHdtL2Norm(void)
{
//    double dotProd = 0;
//    int nlvl = dHdtVector.size(); // nlvl=1
//    Vector<LevelData<FArrayBox>* > x=dHdtVector;
//    m_amrMask=dHdtVector;
//    Vector<LevelData<FArrayBox>* > m_amrData=dHdtVector;
//    BisiclesVector *xmask;
//    // xmask->m_amrMask = Vector<LevelData<FArrayBox>*>(nlvl,NULL);
//    refRatio=Vector<int>(nlvl);
//    for (int lvl=0; lvl < nlvl; lvl++)
//    {
//      refRatio[lvl]=1; // double check with Hans
//      m_amrMask[lvl]=NULL;
//    }
//    // xmask->defineMask(dHdtVector,refRatio);

//    // construct mask
//   long masksum = 0;
//   refRatio=Vector<int>(nlvl);

//   for (int lvl=0; lvl < nlvl; lvl++)
//   {
//     LevelData<FArrayBox>& ldf = *m_amrData[lvl];
//     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
//     m_amrMask[lvl] = new LevelData<FArrayBox>(dbl, 1, IntVect::Zero);
//     LevelData<FArrayBox>& mask = *m_amrMask[lvl];

//     // Loop over coarse dbl to set masks
//     for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit) {
//       // initialize this level mask to 1
//       const Box& box = dbl[dit()];
//       FArrayBox& fab = mask[dit()];
//       fab.setVal(1.0);
//       long boxsum = box.numPts();

//       // For each box in this dbl, mask=0 on each coarsened fine dbl box
//       if (lvl < nlvl-1) { // if there is a finer level
//         // loop over ALL finer level boxes, processes -> layout iterator
//         DisjointBoxLayout fdbl = m_amrData[lvl+1]->disjointBoxLayout();
//         for (LayoutIterator flit= fdbl.layoutIterator(); flit.ok(); ++flit) {
//           // set mask to 0 in intersection of any of coarse fine boxes
//           Box maskb(fdbl[flit()]);
//           maskb.coarsen(refRatio[lvl]);
//           maskb &= box;
//           if (!maskb.isEmpty())
//           {
//             fab.setVal(0.0, maskb, 0);
//             boxsum -= maskb.numPts();
//           }
//         }
//       }
//       // Check with mask sum vs. length
//       CH_assert(boxsum == (long) fab.sum(0));
//       pout() << "  Level " << lvl << ", box " << dbl[dit()] <<
//         ", sum of mask: " << boxsum << endl;
//       masksum += boxsum;
//     }
//   }

   


//    // dot product
//    for (int lvl=0; lvl < nlvl; lvl++)
//    {
//       LevelData<FArrayBox>& ldf = *dHdtVector[lvl];
//       LevelData<FArrayBox>& xldf = *(x[lvl]);
//       LevelData<FArrayBox>& mldf = *m_amrMask[lvl]; // check with Hans, mask? must?
//       DisjointBoxLayout dbl = ldf.disjointBoxLayout();
//       DataIterator dit = ldf.dataIterator();
//       // currently dit.ok is only 1
//       for (dit.reset(); dit.ok(); ++dit) {
//          const Box& box = dbl[dit()];
//          const FArrayBox& fab = ldf[dit()];
//          const FArrayBox& xfab = xldf[dit()];
//          const FArrayBox& mask = mldf[dit()];
//          FORT_MASKDOTPROD(CHF_CONST_FRA(fab), CHF_CONST_FRA(xfab), 
//              CHF_CONST_FRA(mask), CHF_BOX(box), CHF_REAL(dotProd));
//       }
//    }
//    Real mpidot = 0;
//    MPI_Allreduce(&dotProd, &mpidot, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
//    // cout<< "bisicles_vector.cpp mpidot "<<dotProd<< std::endl;
//    return mpidot;
}


double BisiclesVector::HL2Norm(void)
{
   double L2norm;
   //                  a_phi,             a_nRefFine,  a_dxCrse,Interval a_comps, int a_p  ,int a_lBase
   // Vector<LevelData<FArrayBox>* >& , Vector<int>& ,            ,(0,0),            2,         0)
   L2norm=computeNorm(HVector, refineRatio , amrDx[0], Interval(0,0),2, 0);
   // cout<< "from bisicles_vector: after norm "<<L2norm<< std::endl;
   // MayDay::Abort("------------------------------- stop here --------------------------------------");
   return L2norm;
}

void BisiclesVector::PrintHL2norm(void)
{
   Real test_min=0;
   int nlvl=HVector.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       test_min=fab.norm(box,1);
       // cout<< " ------------------- eprint() ....................... "<< std::endl;
       // cout<< "                       eprint sum of H vector  "<<test_min<< std::endl;
      } 
  }
// MayDay::Abort("------------------------------- stop here --------------------------------------");
}

void BisiclesVector::defineMask(Vector<LevelData<FArrayBox>* > m_amrData, Vector<int> refRatio)
{
   int nlvl=m_amrData.size(); // nlvl=1
   m_amrMask=dHdtVector;
   cout<< "bisicles_vector.cpp 2 length of m_amrMask "<<dHdtVector.size()<< std::endl;
   // m_amrMask = Vector<LevelData<FArrayBox>*>(nlvl);
  long masksum = 0;
  refRatio=Vector<int>(nlvl);

  for (int lvl=0; lvl < nlvl; lvl++)
  {
    refRatio[lvl]=1; // double check with Hans
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    m_amrMask[lvl] = new LevelData<FArrayBox>(dbl, 1, IntVect::Zero);
    LevelData<FArrayBox>& mask = *m_amrMask[lvl];

    // Loop over coarse dbl to set masks
    for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit) {
      // initialize this level mask to 1
      const Box& box = dbl[dit()];
      FArrayBox& fab = mask[dit()];
      fab.setVal(1.0);
      long boxsum = box.numPts();

      // For each box in this dbl, mask=0 on each coarsened fine dbl box
      if (lvl < nlvl-1) { // if there is a finer level
        // loop over ALL finer level boxes, processes -> layout iterator
        DisjointBoxLayout fdbl = m_amrData[lvl+1]->disjointBoxLayout();
        for (LayoutIterator flit= fdbl.layoutIterator(); flit.ok(); ++flit) {
          // set mask to 0 in intersection of any of coarse fine boxes
          Box maskb(fdbl[flit()]);
          maskb.coarsen(refRatio[lvl]);
          maskb &= box;
          if (!maskb.isEmpty())
          {
            fab.setVal(0.0, maskb, 0);
            boxsum -= maskb.numPts();
          }
        }
      }
      // Check with mask sum vs. length
      CH_assert(boxsum == (long) fab.sum(0));
      pout() << "  Level " << lvl << ", box " << dbl[dit()] <<
        ", sum of mask: " << boxsum << endl;
      masksum += boxsum;
    }
  }
}

void BisiclesVector::DataAxpy(double a, double x)
{
   data += a * x;
}

void BisiclesVector::dHdtAxpy(double a, BisiclesVector *x,AmrIceHolderClass *c_AmrIceHolderPtr)
{
   // this is scale
   Vector<LevelData<FArrayBox>* > m_amrData = x->dHdtVector; // norm of dHdt=0 before computation
   LevelData<FArrayBox>& ldf = *m_amrData[0];
      DataIterator dit = ldf.dataIterator();
      DisjointBoxLayout dbl = ldf.disjointBoxLayout();
      Real test_min=0;
      for (dit.reset(); dit.ok(); ++dit) {
         FArrayBox& fab = ldf[dit()];
         const Box& box = dbl[dit()];
         test_min=fab.norm(box);
      } // test_min is always 0, need to debug
      // cout<< "from bisicles_vector: before compute dHdt norm of test_min "<<test_min<< std::endl;

   int nlvl=dHdtVector.size();
   int m_comp = dHdtVector[0]->nComp();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *dHdtVector[lvl];
     LevelData<FArrayBox>& xldf = *m_amrData[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
     {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       const FArrayBox& xfab = xldf[dit()];
       // if (!inplace) // double check with Hans what it means
         fab.copy(xfab, box); 
       fab.mult(a, box, 0, m_comp); 
     }
   }
   
      dit = ldf.dataIterator();
      dbl = ldf.disjointBoxLayout();
      test_min=0;
      for (dit.reset(); dit.ok(); ++dit) {
         FArrayBox& fab = ldf[dit()];
         const Box& box = dbl[dit()];
         test_min=fab.norm(box);
      } // test_min is always 0, need to debug
      // cout<< "from bisicles_vector: after compute dHdt norm of test_min "<<test_min<< std::endl;
      // cout<< "from bisicles_vector: after compute dHdt axpy double a "<<a<< std::endl; // a is changing

   c_AmrIceHolderPtr->SetAmrdHdtPtr(&dHdtVector);
   c_AmrIceHolderPtr->SetAmrdHdt(dHdtVector);
}


void BisiclesVector::HAxpy(double a, BisiclesVector *x,AmrIceHolderClass *c_AmrIceHolderPtr)
{
   Real test_min=0;
   Vector<LevelData<FArrayBox>* > xHVector = x->HVector; // 

   int nlvl=HVector.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     LevelData<FArrayBox>& xldf = *xHVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DisjointBoxLayout xdbl = xldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
     {
       const Box& box = dbl[dit()];
       const Box& xbox = xdbl[dit()];
       FArrayBox& fab = ldf[dit()];
       const FArrayBox& xfab = xldf[dit()];
       test_min=fab.norm(box,1);
       // cout<< "from bisicles_vector: before axpy norm of y "<<test_min << " a " << a<< std::endl;
       // test_min=xfab.norm(xbox,1);
       // cout<< "from bisicles_vector: before axpy norm of x "<<test_min<< std::endl;
       
       fab.plus(xfab,a);

       // test_min=fab.norm(box,1);
       // cout<< "from bisicles_vector: after axpy norm of y "<<test_min<< std::endl;

     }
   }
   
   // MayDay::Abort("------------------------------- stop here --------------------------------------");
   // c_AmrIceHolderPtr->SetAmrHPtr(&HVector);
   // c_AmrIceHolderPtr->SetAmrH(HVector);
}

void BisiclesVector::SaveSnapshot(AmrIceHolderClass *c_AmrIceHolderPtr)
{
   AmrIce *amrObjHolderPtr;
   reshapeAndFill(iceState.ice_thickness, HVector);
   c_AmrIceHolderPtr->SetAmrIceState(iceState);
   amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
   amrObjHolderPtr->writePlotFile();
}

void BisiclesVector::DataPrint()
{
//   cout << data << '\n';
}
