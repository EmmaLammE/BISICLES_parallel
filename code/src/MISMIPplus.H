#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// collection of classes encapsulating the MISMIP+ problem

#ifndef _MISMIPplus_H_
#define _MISMIPplus_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"
#include "SurfaceFlux.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

/// melt forcing -- set specific experiment via parmparse
class MISMIPmelt: public SurfaceFlux
{
public:
  // default constructor
  MISMIPmelt()
  {
  }
  
  
  // constructor
  MISMIPmelt(ParmParse& a_pp);

  // virtual destructor
  virtual ~MISMIPmelt();

  /// factory method: return a pointer to a new SurfaceFlux object
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** 
      @param a_flux output flux data 
      @param     a_amrIce reference to the ice sheet state
      @param     a_level mesh level of a_flux
      @param     a_dt current timestep

      a_dt is included in case one needs integrals or averages over a
      timestep.  flux should be defined in meters per year in the current 
      implementation. 
  */  
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);
  
protected:
  string m_experiment;
  
};

// N.B, bedrock geometry is defined via the MISOMIPBedrockElevation class, defined in MarineIB.H

#include "NamespaceFooter.H"


#endif // multiple-include preventer
