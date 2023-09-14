#ifndef BISICLES_VECTOR_HPP
#define BISICLES_VECTOR_HPP

#include <iostream>
#include <cstdlib>
#include "chombo_struct.hpp"
#include "bisicles_holder.hpp"
#include "LevelData.H"
#include "FArrayBox.H"
#include "N_VectorOps_F.H"
// #include "N_VectorOps.ChF"


using namespace std;

class BisiclesVector: public ChomboStruct, public AmrIceHolderClass{
   protected:
      double data;
      Vector<LevelData<FArrayBox>* > *dHdtVectorPtr;
      Vector<LevelData<FArrayBox>* > dHdtVector;
      Vector<LevelData<FArrayBox>* > *HVectorPtr;
      Vector<LevelData<FArrayBox>* > HVector;
      IceSheetState iceState;
      IceSheetState *iceStatePtr;

      Vector<LevelData<FArrayBox>*> m_amrMask;

      long m_length = -1;
      Vector<int> refineRatio;
      Vector<Real> amrDx;
      // int num_of_grid_points; // num_of_grid_points=n^2;
      // Vector<Vector<LevelData<FArrayBox>* > > crsedHdtVect;
   public:
      int num_grid_points;

      BisiclesVector(int num_grid_points, AmrIceHolderClass *c_AmrIceHolderPtr);
      ~BisiclesVector(void);

      //int GetNumRows(void);
      Vector<LevelData<FArrayBox>* > GetdHdtVector(void);
      Vector<LevelData<FArrayBox>* > *GetdHdtVectorPtr(void);
      const Vector<LevelData<FArrayBox>* >& GetHVector(void) const;
      Vector<LevelData<FArrayBox>* > *GetHVectorPtr(void);
      IceSheetState GetIceState(void);
      IceSheetState *GetIceStatePtr(void);
      void SetdHdtVector(Vector<LevelData<FArrayBox>* > vector_value,AmrIceHolderClass *c_AmrIceHolderPtr);
      void SetHVector(BisiclesVector* vector_value);
      // double *DataPtrGet(void);
      double DataGet(void);

      void dHdtSetVal2All(double val,AmrIceHolderClass *c_AmrIceHolderPtr);
      void HSetVal2All(double val,AmrIceHolderClass *c_AmrIceHolderPtr);
      double *GetdHdtDataPtr(void);
      Vector<LevelData<FArrayBox>* >* GetHDataPtr(void);
      double dHdtL2Norm(void);
      double HL2Norm(void);
      void dHdtAxpy(double a, BisiclesVector *x,AmrIceHolderClass *c_AmrIceHolderPtr);
      void HAxpy(double a, BisiclesVector *x,AmrIceHolderClass *c_AmrIceHolderPtr);
      void defineMask(Vector<LevelData<FArrayBox>* > m_amrData, Vector<int> refRatio);
      void PrintHL2norm(void);
      void PrintL2norm(Vector<LevelData<FArrayBox>* > Vect2Print);
      void SaveSnapshot(AmrIceHolderClass *c_AmrIceHolderPtr);
      void PrintLevelData(Vector<LevelData<FArrayBox>* > level_data_to_print);
      void PFPrintLevelData(BisiclesVector *x);
      void PrintAllFArray(Vector<LevelData<FArrayBox>* > level_data_to_print);
      void PrintAllFluxBox(Vector<LevelData<FluxBox>* > level_data_to_print);
};

#endif
