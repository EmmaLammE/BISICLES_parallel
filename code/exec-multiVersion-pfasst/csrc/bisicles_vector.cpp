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


   Vector<LevelData<FArrayBox>* > H_ref=c_AmrIceHolderPtr->GetAmrH();
   
   int num_levels_bisicles_vec=c_AmrIceHolderPtr->GetAmrNumLvl();
   refineRatio=amrObjHolderPtr->refRatios();
   amrDx=amrObjHolderPtr->amrDx();

   Vector<LevelData<FArrayBox>* > HVector_ref(num_levels_bisicles_vec, NULL);
   HVector=HVector_ref;
   Vector<LevelData<FArrayBox>* > dHdtVector_ref(num_levels_bisicles_vec, NULL);
   dHdtVector=dHdtVector_ref;
   for (int lvl=0; lvl<num_levels_bisicles_vec;lvl++)
   {
      // Vector<DisjointBoxLayout> current_grid = H_ref
      const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl); // use: grid(), or getBoxes(), or surroundingNode? m_vect_coordSys()?
      // HVector[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero); // FIXED the num of component=1, i.e. one d?? only Vx no Vy??
      // dHdtVector[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero);
      // refineRatio[lvl]=2; // double check with Dan

      // get SigmaCS mapping
      RefCountedPtr<LevelSigmaCS> vect_coordSys = amrObjHolderPtr->geometry(lvl);
      const DisjointBoxLayout& phy_coord = vect_coordSys->grids();
      // LevelData<FArrayBox>* my_tempH = &(vect_coordSys->getH());
      // const DisjointBoxLayout& phy_coord = my_tempH->getBoxes();
      HVector[lvl] = new LevelData<FArrayBox>(phy_coord,1,IntVect::Zero); 
      dHdtVector[lvl] = new LevelData<FArrayBox>(phy_coord,1,IntVect::Zero); 
      IntVect sigmaCSGhost = vect_coordSys->ghostVect();
      HVector[lvl]->define(phy_coord, 1, sigmaCSGhost);
      dHdtVector[lvl]->define(phy_coord, 1, sigmaCSGhost);
      // cout<<"...bisicles vector creating sigmaCSGhost "<<sigmaCSGhost<<", level "<<lvl<<endl;
      DataIterator dit = (*HVector[lvl]).dataIterator();
      int iter=0;
      for (dit.reset(); dit.ok(); ++dit)
      {
         // cout<<" dit iter "<<iter<<endl;
         iter ++;
         // const FArrayBox& curHArray = (*HVector[lvl])[dit()];
         // Box edgeBox = surroundingNodes(curHArray.box(), 0); // sec var is direction
         // cout<<"...bisicles vector creating edgeBox "<<edgeBox<<endl;
      }

      // DataIterator dit = (*HVector[lvl]).dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
         (*HVector[lvl])[dit].setVal(0.0);
         (*dHdtVector[lvl])[dit].setVal(0.0);
      }
        
   }

   // bisicles vectors are intialized in pfasst corrrectly in different spatial levels
   // int nlvl=HVector.size();
   // for (int lvl=0; lvl < nlvl; lvl++)
   // {
   //   LevelData<FArrayBox>& ldf = *HVector[lvl];
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
   

   reshapeAndFill(iceState.ice_thickness, HVector);


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

const Vector<LevelData<FArrayBox>* >& BisiclesVector::GetHVector(void) const
{
   
      // Real test_min=0;
      // for (int lvl=0; lvl < HVector.size(); lvl++)
      // {
      // LevelData<FArrayBox>& ldf = *HVector[lvl];
      // DisjointBoxLayout dbl = ldf.disjointBoxLayout();
      // DataIterator dit = ldf.dataIterator();
      // for (dit.reset(); dit.ok(); ++dit) 
      //    {
      //    const Box& box = dbl[dit()];
      //    FArrayBox& fab = ldf[dit()];
      //    test_min=fab.norm(box,1);
      //    //  cout<<"      bisicles... box "<<box<<endl;
      //     cout<<"      bisicles... print L2 norm: "<<test_min<<endl;
      //    } 
      // }
      // cout<<"Hvector "<<HVector<<endl;
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

// double check with Hans/Dan, c_AmrIceHolderPtr does not seem to be used in this function
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
       fab.copy(xfab,box); // copy xfab to box
       test_min=fab.norm(box,1);
      } 
  
  }
//   Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      // cout<<"-- 5 bis_vec:\n";
      // PrintLevelData(constH);
}


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
 
       fab.setVal(val);      
       fab.plus(xfab,a);
       test_min=fab.norm(box,2);
      } 
  
  }
}


// double *BisiclesVector::DataPtrGet(void)
// {
//    double *data = (double *)calloc(num_grid_points, sizeof(double));
//    return data;
// }

double BisiclesVector::DataGet(void)
{
   return data;
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
}

void BisiclesVector::HSetVal2All(double val,AmrIceHolderClass *c_AmrIceHolderPtr)
{
   // // double check with Dan: do we need to assign values of all zero??
   int i=0;
   for (int lvl=0; lvl < HVector.size(); lvl++)
   {
   //   cout<<"lvl "<<lvl<<endl;
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) {
      const Box& box = dbl[dit()];
      (*HVector[lvl])[dit()].setVal(val);
     }
   }

   // setting some specific entries for packing/unpacking testing
   // int count=0;
   // for (int lvl=0; lvl < HVector.size(); lvl++)
   // {
   //   LevelData<FArrayBox>& ldf = *HVector[lvl];
   //   DisjointBoxLayout dbl = ldf.disjointBoxLayout();
   //   DataIterator dit = ldf.dataIterator();
   //   for (dit.reset(); dit.ok(); ++dit) {
   //    const Box& box = dbl[dit()];
   //    FArrayBox& fab = ldf[dit()];
   //    if(count==0){ // only setting specific values in first box
   //       cout<<"setting specific values in box "<<box<<endl;
   //       fab.setSpecificValTest(box);
   //    }
   //    count++;
   //   }
   // }
}

double *BisiclesVector::GetdHdtDataPtr(void)
{
   // double check with Dan to find the total num of entries
   double *data = (double *)calloc(num_grid_points, sizeof(double));
   return data;
}

Vector<LevelData<FArrayBox>* >* BisiclesVector::GetHDataPtr()
{
   // count the grid cells in each level and store into vecotr
   // AmrIce *amrObjHolderPtr;
   // amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
   int nlvl=HVector.size();
   Vector<int> num_of_grid_points_vec(nlvl);
   // int total_num_grid_points=0;
   // Real test_min=0;
   // for (int lvl=0; lvl < nlvl; lvl++)
   // {
   //   const DisjointBoxLayout& current_grid_size=amrObjHolderPtr->grids(lvl);
   //   LayoutIterator lit = current_grid_size.layoutIterator();
   //   LevelData<FArrayBox>& ldf = *HVector[lvl];
   //   DisjointBoxLayout dbl = ldf.disjointBoxLayout();
   //   DataIterator dit = ldf.dataIterator();
   //   for (dit.begin(); dit.ok(); ++dit)
   //   {
   //     const Box& thisBox = current_grid_size.get(dit());
   //     num_of_grid_points_vec[lvl] += thisBox.numPts();
   //     const Box& box = dbl[dit()];
   //     FArrayBox& fab = ldf[dit()];
   //     test_min=fab.norm(box,1);
   //    //  cout<<"      in packing H box "<<box<<endl;
   //    //  cout<<"      in packing H L2 norm: "<<test_min<<endl;
   //    }
   //    total_num_grid_points +=num_of_grid_points_vec[lvl];
   // }
   // allocate the vector for H and store the pointers in value
   // double *values = (double *)calloc(total_num_grid_points, sizeof(double));
   // cout<<"       in packing HVector "<<HVector<<endl;
   return &HVector;
   // Vector<LevelData<FArrayBox>* > constH=c_AmrIceHolderPtr->GetAmrH();
      // cout<<"-- 3 bis_vec:\n";
      // PrintLevelData(constH);
}

// double BisiclesVector::DataNorm(void)
// {
//    cout<< "bisicles_vector.cpp abs data "<<abs(data)<< std::endl;
//    return abs(data);
// }

double BisiclesVector::dHdtL2Norm(void)
{
}


double BisiclesVector::HL2Norm(void)
{
   double L2norm;
   // cout<<"      bisicles... in compute norm, HVector "<<HVector<<", size "<<HVector.size()<<endl;
   Real test_min=0;
   for (int lvl=0; lvl < HVector.size(); lvl++)
   {
   //   cout<<"      bisicles... level "<<lvl<<endl;
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       test_min=fab.norm(box,1);
      //  cout<<"      bisicles... box "<<box<<endl;
      //  cout<<"      bisicles... print L2 norm: "<<test_min<<endl;
      } 
   }
   // input: vector to compute norm, refine ratio, dx, interval?, a_p(1-L1,2-L2),
   // hanging up here with multiple procs
   L2norm=computeNorm(HVector, refineRatio , amrDx[0], Interval(0,0),2, 0); 
   // cout<<"      ... done L2norm computation\n";
   return L2norm;
}

void BisiclesVector::PrintHL2norm(void)
{
   Real test_min=0;
   int box_index=0;
   int nlvl=HVector.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *HVector[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     box_index=0;
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       test_min=fab.norm(box,1);
      //  if(box_index==1){
      //       cout<<" in printing L2 norm, the box to be packed box "<<box<<endl;
      //       fab.printAll(box);
      //       cout<<"\n";
      //   }
        box_index++;
       cout<<"   ... box "<<box<<", print L2 norm: "<<test_min<<endl;
      } 
  }
}

void BisiclesVector::PrintL2norm(Vector<LevelData<FArrayBox>* > Vect2Print)
{
   Real test_min=0;
   int nlvl=Vect2Print.size();
   // cout<<" # of levels "<<nlvl<<endl;
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *Vect2Print[lvl];
   //   cout<<"  level data "<<ldf.getBoxes();
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
   //   cout<<"  box "<<dbl;
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       test_min=fab.norm(box,1);
       cout<<"   "<<test_min<<endl;
      } 
   }
}

void BisiclesVector::defineMask(Vector<LevelData<FArrayBox>* > m_amrData, Vector<int> refRatio)
{
   int nlvl=m_amrData.size(); // nlvl=1
   m_amrMask=dHdtVector;
   cout<< "bisicles_vector.cpp 2 length of m_amrMask "<<dHdtVector.size()<< std::endl;
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
         // test_min=fab.norm(box);
      } // test_min is always 0, need to debug

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
         // test_min=fab.norm(box);
      } // test_min is always 0, need to debug

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
       
       fab.plus(xfab,a);


     }
   }
}

void BisiclesVector::SaveSnapshot(AmrIceHolderClass *c_AmrIceHolderPtr)
{
   AmrIce *amrObjHolderPtr;
   reshapeAndFill(iceState.ice_thickness, HVector);
   c_AmrIceHolderPtr->SetAmrIceState(iceState);
   amrObjHolderPtr=c_AmrIceHolderPtr->GetAmrIceObjPtr();
   amrObjHolderPtr->writePlotFile(); 
   // AmrIO.cpp 192, the prefix is read in at the beginning of the AmrIO as plot_prefix, later edited to filename
   // the actual write is in 1015
}

void BisiclesVector::PrintLevelData(Vector<LevelData<FArrayBox>* > level_data_to_print)
{
   int nlvl=level_data_to_print.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *level_data_to_print[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     cout<<"  level "<<lvl<<endl;
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
      //  test_min=fab.norm(box,1);
      cout<<"    box "<<box<<endl;
      } 
   }
}

void BisiclesVector::PFPrintLevelData(BisiclesVector *x)
{
   Vector<LevelData<FArrayBox>* > level_data_to_print = x->HVector;
   int nlvl=level_data_to_print.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *level_data_to_print[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     cout<<"  level "<<lvl<<endl;
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
      //  test_min=fab.norm(box,1);
      cout<<"    box "<<box<<endl;
      } 
   }
}

void BisiclesVector::PrintAllFArray(Vector<LevelData<FArrayBox>* > level_data_to_print)
{
   int nlvl=level_data_to_print.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FArrayBox>& ldf = *level_data_to_print[lvl];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     cout<<"level "<<lvl<<endl;
     for (dit.reset(); dit.ok(); ++dit) 
      {
       const Box& box = dbl[dit()];
       FArrayBox& fab = ldf[dit()];
       fab.printAll(box);
      } 
   }
}

void BisiclesVector::PrintAllFluxBox(Vector<LevelData<FluxBox>* > level_data_to_print)
{
   int nlvl=level_data_to_print.size();
   for (int lvl=0; lvl < nlvl; lvl++)
   {
     LevelData<FluxBox>& ldf = *level_data_to_print[lvl];
   //   LevelData<FArrayBox>& ldf1 = *a_H[lev];
     DisjointBoxLayout dbl = ldf.disjointBoxLayout();
     DataIterator dit = ldf.dataIterator();
     for (dit.reset(); dit.ok(); ++dit) 
      {
        FluxBox& thisFlux = ldf[dit];
        FArrayBox& thisFluxArrayBox = thisFlux.getFlux(0); // for now only check 0 direction
        const Box& thisFluxBox = thisFlux.box();
        thisFluxArrayBox.printAll(thisFluxBox);
      } 
   }
}
