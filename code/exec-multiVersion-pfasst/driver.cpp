#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// driver.cpp
//
//===========================================================================
#include <iostream>
#include <fstream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "AmrIce.H"
#include "ConstitutiveRelation.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFriction.H"
#include "BasalFrictionRelation.H"
#include "MuCoefficient.H"
#include "twistyStreamFriction.H"
#include "singularStreamFriction.H"
#include "GaussianBumpFriction.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "HumpIBC.H"
#include "LevelDataIBC.H"
#include "MultiLevelDataIBC.H"
#include "IceInternalEnergyIBC.H"
#include "LevelDataTemperatureIBC.H"
#include "VerticalConductionInternalEnergyIBC.H"
#include "LevelDataBasalFriction.H"
#include "PiecewiseLinearFlux.H"
#include "SurfaceFlux.H"
#include "IceConstants.H"
#include "AMRDamage.H"
#include "AMRMelange.H"
#include "DamageConstitutiveRelation.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
//#include "LevelDataSurfaceFlux.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "FineInterp.H"
#include "ReadLevelData.H"
#include "SundialsUtil.H"
#include "FABView.H"


#ifdef CH_USE_PETSC
#include "petsc.h"
#endif 
#include "Regression.H"

// include chombo-pfasst bindings
#include "chombo_pfasst.hpp"
#include "bisicles_vector.hpp"
#include "bisicles_fortran.hpp"
#include "bisicles_holder.hpp"
#include "pfasst_bisicles_setup.hpp"
#include <chrono>
using namespace std::chrono;
// #include <windows.h>
// #include "bisicles_vector.hpp"


/// types of basal friction (beta) distributions
/** SinusoidalBeta is the one for exp C in Pattyn et al (2008)
    guassianBump is used for the MISMIP3D perturbations tests.
 */
enum basalFrictionTypes {constantBeta = 0,
                         sinusoidalBeta,
                         sinusoidalBetay,
                         twistyStreamx,
       gaussianBump,
       singularStream,
                         NUM_BETA_TYPES};

/// main program for 2D ice sheet models
/**
   \callgraph
 */
int main(int argc, char* argv[]) {
  
int ierr = 0;
int rank, number_procs;
int new_space_rank, new_space_size, new_time_rank, new_time_size;
#ifdef CH_USE_PETSC
  ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
#else
  #ifdef CH_MPI
    // MPI_Comm sub_comm;
    MPI_Init(&argc, &argv);
  #endif 
#endif // end petsc conditional

FineInterp::s_default_boundary_limit_type = 0;

{ 

    if(argc < 2) 
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    char* in_file = argv[1];
    //std::string in_file_mangled(in_file);
    //in_file_mangled += "_mangled";
    //FileMangler(in_file,in_file_mangled.c_str());
    //ParmParse pp(argc-2,argv+2,NULL,in_file_mangled.c_str());

  
    ParmParse pp(argc-2,argv+2,NULL,in_file);
    ParmParse ppfasst("pf");
    bool USE_PF;
    int num_time_procs;
    ppfasst.get("USE_PF", USE_PF);  
    ppfasst.get("num_time_procs", num_time_procs); 
    // SYSTEM_INFO siSysInfo;
    // GetSystemInfo(&siSysInfo); 
    // int num_hardware_cores = siSysInfo.dwNumberOfProcessors;
    
    // Begin nested scope
    #ifdef CH_MPI
        MPI_Comm pf_comm; 
        MPI_Comm sub_comm; 
        if (USE_PF){
        // global mpi
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
        if(number_procs<num_time_procs){
          pout()<<"The number of time processor greater than total mpi processor. Number of time processor input: "
                << num_time_procs
                << ". Number of total processor: "
                << number_procs;
          MayDay::Error("Please check your input file and args for mpirun. Stop here.");
        }
        // split mpi into space mpi
        int space_color = round(rank%num_time_procs);
        MPI_Comm_split(MPI_COMM_WORLD, space_color, rank, &Chombo_MPI::comm); // make sure what key should put in here
        MPI_Comm_rank(Chombo_MPI::comm, &new_space_rank); // Get the rank in the new communicator
        MPI_Comm_size(Chombo_MPI::comm, &new_space_size); // Get the size of the new communicator
        // split mpi into time mpi
        int time_color = round(rank/num_time_procs);
        MPI_Comm_split(MPI_COMM_WORLD, time_color, rank, &pf_comm);
        MPI_Comm_rank(pf_comm, &new_time_rank); // Get the rank in the new communicator
        MPI_Comm_size(pf_comm, &new_time_size); // Get the size of the new communicator
        // cout<<"Chombo_MPI::comm "<<Chombo_MPI::comm<<endl;
        // should be correct, but double check if each time processor has the same num of bisicles processors
        if(new_time_size*new_space_size != number_procs){
            MayDay::Error("The number of time and space processors doesn't match the total number of processors. Stop here.");
        }
        }
        else{
        MPI_Barrier(Chombo_MPI::comm);
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
        MPI_Comm_size(Chombo_MPI::comm, &number_procs);
        }
    #else
          rank=0;
          number_procs=1;
    #endif 
    // MayDay::Error("Stop here.");


    ParmParse pp2("main");

    std::string poutBaseName = "pout";
    pp2.query("poutBaseName",poutBaseName);
    setPoutBaseName(poutBaseName);
    //make use of number_procs and rank for at least something.
    pout() << "number_procs = " << number_procs << ", rank = " << rank << std::endl;

    
    RealVect domainSize;
    Vector<Real> domSize(SpaceDim);
    pp2.getarr("domain_size", domSize, 0, SpaceDim);
    domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));

    // AmrIce amrPtr;
    // amrPtr.setParmParsePrefix("crse.");

    // create two AmrIce objects -- one coarse, one fine
    AmrIce amrObjectFine, amrObjectCrse;

    amrObjectFine.setParmParsePrefix("fine.");
    
    amrObjectCrse.setParmParsePrefix("crse.");

    
    // ---------------------------------------------
    // set constitutive relation & rate factor
    // For now, set the same for both AmrIce objects
    // ---------------------------------------------

    Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
    {
      ParmParse ppc("constants");
      ppc.query("seconds_per_unit_time",seconds_per_unit_time);
    }
    
    std::string rateFactorType = "constRate";
    pp2.query("rateFactor", rateFactorType);
    if (rateFactorType == "constRate")
      {
  ParmParse crPP("constRate");
  Real A = 9.2e-18 * seconds_per_unit_time/SECONDS_PER_TROPICAL_YEAR;
  crPP.query("A", A);
  ConstantRateFactor rateFactor(A);
  amrObjectFine.setRateFactor(&rateFactor);
        amrObjectCrse.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "arrheniusRate")
      {
  ArrheniusRateFactor rateFactor(seconds_per_unit_time);
  ParmParse arPP("ArrheniusRate");
  amrObjectFine.setRateFactor(&rateFactor);
        amrObjectCrse.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "patersonRate")
      {
  PatersonRateFactor rateFactor(seconds_per_unit_time);
  ParmParse arPP("PatersonRate");
  amrObjectFine.setRateFactor(&rateFactor);
        amrObjectCrse.setRateFactor(&rateFactor);
      }
    else if (rateFactorType == "zwingerRate")
      {
  ZwingerRateFactor rateFactor(seconds_per_unit_time);
  ParmParse arPP("ZwingerRate");
  amrObjectFine.setRateFactor(&rateFactor);
        amrObjectCrse.setRateFactor(&rateFactor);
      }

    ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");

    if (constRelPtr == NULL)
      {
  MayDay::Error("undefined constitutiveRelation in inputs");
      }

    amrObjectFine.setConstitutiveRelation(constRelPtr);
    amrObjectCrse.setConstitutiveRelation(constRelPtr);
 
    std::string basalRateFactorType = "";
    pp2.query("basalRateFactor", basalRateFactorType);
    
    if (basalRateFactorType == "patersonRate")
      {
  PatersonRateFactor rateFactor(seconds_per_unit_time);
  rateFactor.setA0(1.0);
  amrObjectFine.setBasalRateFactor(&rateFactor);
        amrObjectCrse.setBasalRateFactor(&rateFactor);
      }

    // ---------------------------------------------
    // set surface flux. 
    // ---------------------------------------------

    SurfaceFlux* surf_flux_ptr = SurfaceFlux::parse("surfaceFlux");
    if (surf_flux_ptr == NULL)
      {
  const std::string err("failed to parse surfaceFlux (maybe you have the old style surface_flux_type?");
  pout() << err << endl;
  MayDay::Error(err.c_str());
      }

    amrObjectFine.setSurfaceFlux(surf_flux_ptr);
    amrObjectCrse.setSurfaceFlux(surf_flux_ptr);    

    //delete surf_flux_ptr;
    // ---------------------------------------------
    // set basal (lower surface) flux. 
    // ---------------------------------------------
    
    SurfaceFlux* basal_flux_ptr = SurfaceFlux::parse("basalFlux");
    if (basal_flux_ptr == NULL)
      {
  const std::string err("failed to parse basalFlux (maybe you have the old style basal_flux_type?");
  pout() << err << endl;
  MayDay::Error(err.c_str());
      }

    amrObjectFine.setBasalFlux(basal_flux_ptr);
    amrObjectCrse.setBasalFlux(basal_flux_ptr); 
    //delete basal_flux_ptr;
     // ---------------------------------------------
    // set topography (bedrock) flux. 
    // ---------------------------------------------
    
    SurfaceFlux* topg_flux_ptr = SurfaceFlux::parse("topographyFlux");
    if (topg_flux_ptr == NULL)
      {
  topg_flux_ptr = new zeroFlux();
      }
    amrObjectFine.setTopographyFlux(topg_flux_ptr);
    amrObjectCrse.setTopographyFlux(topg_flux_ptr);     

    // ---------------------------------------------
    // set mu coefficient
    // ---------------------------------------------
    {
      MuCoefficient* muCoefPtr =  MuCoefficient::parseMuCoefficient("muCoefficient");

      if (muCoefPtr == NULL)
        {
          const std::string err("failed to parse muCoefficient");
          pout() << err << endl;
          MayDay::Error(err.c_str());
        }
      
      amrObjectFine.setMuCoefficient(muCoefPtr);
      amrObjectCrse.setMuCoefficient(muCoefPtr);      
      delete muCoefPtr;
    }

    // ---------------------------------------------
    // set basal friction coefficient and relation
    // ---------------------------------------------

    ParmParse geomPP("geometry");
    
    BasalFriction* basalFrictionPtr 
      = BasalFriction::parse("geometry", domainSize);
    
    if (basalFrictionPtr == NULL)
      {
  MayDay::Error("undefined  geometry.beta_type in inputs");
      }
    
    amrObjectFine.setBasalFriction(basalFrictionPtr);
    amrObjectCrse.setBasalFriction(basalFrictionPtr);
    
    BasalFrictionRelation* basalFrictionRelationPtr 
      = BasalFrictionRelation::parse("main",0);

    amrObjectFine.setBasalFrictionRelation(basalFrictionRelationPtr);
    amrObjectCrse.setBasalFrictionRelation(basalFrictionRelationPtr);    

    // ---------------------------------------------
    // set IBC -- this includes initial ice thickness, 
    // and basal geometry
    // ---------------------------------------------

    
    IceThicknessIBC* thicknessIBC = NULL;

    std::string problem_type;
    geomPP.get("problem_type", problem_type);
    if (problem_type == "basic")
      {
        thicknessIBC = new BasicThicknessIBC;
      }
    else if (problem_type == "VieliPayne")
      {
        VieliPayneIBC* ibcPtr = new VieliPayneIBC;
        ParmParse pvPP("vieliPayne");

        Real thickness, seaLevel, originElevation;
        RealVect basalSlope;
        pvPP.get("thickness", thickness);
        seaLevel = 0.0;
        pvPP.query("seaLevel", seaLevel);
        
        Vector<Real> vect(SpaceDim);
        pvPP.getarr("basal_slope", vect, 0, SpaceDim);
        basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));

        pvPP.get("originElevation", originElevation);

        ibcPtr->setParameters(thickness, basalSlope, 
                              originElevation, seaLevel);       

        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "marineIceSheet")
      {
        MarineIBC* ibcPtr = new MarineIBC;
        ParmParse mPP("marineIceSheet");
  
        Real  seaLevel;
        seaLevel = 0.0;
        mPP.query("seaLevel", seaLevel);
        
        std::string thicknessType = "constant";
        mPP.query("thickness_type",thicknessType);
        RefCountedPtr<RealFunction<RealVect> > thicknessFunction;
        Vector<RefCountedPtr<RealFunction<RealVect> > > bedrockFunction(1);
        if (thicknessType == "constant")
          {
            Real thickness;
            mPP.get("thickness", thickness);
            RefCountedPtr<RealFunction<RealVect> > ptr(new ConstantRealFunction<RealVect>(thickness));
            thicknessFunction =ptr;
          }
        else if (thicknessType == "compactSupportConstant")
          {
            Real thickness;
            Vector<Real> tmpIntVect(SpaceDim,0); 
            mPP.get("thickness", thickness);
            mPP.getarr("loBound", tmpIntVect, 0, SpaceDim);
            RealVect loBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            mPP.getarr("hiBound", tmpIntVect, 0, SpaceDim);
            RealVect hiBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            
            RefCountedPtr<RealFunction<RealVect> > ptr(new CompactSupportConstantRealFunction(thickness,loBound, hiBound));
            thicknessFunction =ptr;
            
          }
        else if (thicknessType == "circularSupportConstant")
          {
            Real thickness, supportRadius;
            Vector<Real> tmpRealVect(SpaceDim,0); 
            mPP.get("thickness", thickness);
            mPP.get("supportRadius", supportRadius);
            mPP.getarr("center", tmpRealVect, 0, SpaceDim);
            RealVect center(D_DECL(tmpRealVect[0],tmpRealVect[1],tmpRealVect[2]));

            RefCountedPtr<RealFunction<RealVect> > ptr(new CircularSupportConstantRealFunction(thickness,center, supportRadius));
            thicknessFunction =ptr;            
          }
        else if (thicknessType == "compactSupportInclinedPlane")
          {
            Real originThickness;
            Vector<Real> tmpIntVect(SpaceDim,0); 
            mPP.get("origin_thickness", originThickness);
            mPP.getarr("loBound", tmpIntVect, 0, SpaceDim);
            RealVect loBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            mPP.getarr("hiBound", tmpIntVect, 0, SpaceDim);
            RealVect hiBound(D_DECL(tmpIntVect[0], tmpIntVect[1], tmpIntVect[2]));
            RealVect thicknessSlope;
            Vector<Real> vect(SpaceDim);
            mPP.getarr("thickness_slope", vect, 0, SpaceDim);
            thicknessSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            RefCountedPtr<RealFunction<RealVect> > ptr(new CompactSupportInclinedPlaneFunction(originThickness, thicknessSlope,loBound, hiBound));
            thicknessFunction =ptr;
            
          }
        else if (thicknessType == "step")
          {
            int dir=0;
            Real leftThickness, rightThickness, cutoff;
            mPP.get("left_thickness", leftThickness);
            mPP.get("right_thickness", rightThickness);
            mPP.get("x_cutoff", cutoff);
            mPP.query("dir", dir);
            
            RefCountedPtr<RealFunction<RealVect> > ptr(new StepRealFunction(leftThickness, rightThickness, cutoff, dir));
            thicknessFunction =ptr;
          }
        else if (thicknessType == "flowline")
          {
            Real dx; 
            std::string file, set;
            mPP.get("thickness_flowline_dx", dx);
            mPP.get("thickness_flowline_file", file);
            mPP.get("thickness_flowline_set", set);
            RefCountedPtr<RealFunction<RealVect> > ptr(new ExtrudedPieceWiseLinearFlowline(file,set,dx));
            thicknessFunction = ptr;
          }
        else 
          {
            MayDay::Error("bad marineIceSheet.thicknessType");
          }
        
        
        std::string geometry = "plane";
        mPP.query("geometry",geometry);
        
        if (geometry == "plane")
          {
            //inclined plane geometry
            RealVect basalSlope;
            Vector<Real> vect(SpaceDim);
            mPP.getarr("basal_slope", vect, 0, SpaceDim);
            basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            Real originElevation;
            mPP.get("originElevation", originElevation);
            
            
            RefCountedPtr<RealFunction<RealVect> > ptr(new InclinedPlaneFunction(originElevation, basalSlope));
            bedrockFunction[0] =  ptr;
          }
        else if (geometry == "symmetricPlane")
          {
            //inclined plane geometry, symmetric about origin
            RealVect basalSlope;
            Vector<Real> vect(SpaceDim);
            mPP.getarr("basal_slope", vect, 0, SpaceDim);
            basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            
            Real originElevation;
            mPP.get("originElevation", originElevation);
            
            RealVect symmetryPoint(RealVect::Zero);
            mPP.getarr("symmetryPoint", vect, 0, SpaceDim);
            symmetryPoint = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            
            RefCountedPtr<RealFunction<RealVect> > ptr(new SymmetricInclinedPlaneFunction(originElevation, basalSlope, symmetryPoint));
            bedrockFunction[0] =  ptr;
          }
        else if (geometry == "gaussianHump")
          {
            //gaussian bedrock geometry, symmetric about origin
            RealVect center;
            Vector<Real> vect(SpaceDim);
            mPP.getarr("center", vect, 0, SpaceDim);
            center = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            
            RealVect radius;
            mPP.getarr("radius", vect, 0, SpaceDim);
            radius = RealVect(D_DECL(vect[0], vect[1], vect[2]));            
            
            Real magnitude;
            mPP.get("magnitude", magnitude);
            
            Real offset;
            mPP.get("offset", offset);
            
            RefCountedPtr<RealFunction<RealVect> > ptr(new GaussianFunction(center,
                                                                            radius,
                                                                            magnitude,
                                                                            offset));
            bedrockFunction[0] =  ptr;
          }        
        else if (geometry == "regroundingTest")
          {
            //inclined plane geometry with a Gaussian bump
            bedrockFunction.resize(2);
            
            RealVect basalSlope;
            Vector<Real> vect(SpaceDim);
            mPP.getarr("basal_slope", vect, 0, SpaceDim);
            basalSlope = RealVect(D_DECL(vect[0], vect[1], vect[2]));
            Real originElevation;
            mPP.get("originElevation", originElevation);
            
            // compose flat plane with Gaussian hump
            RefCountedPtr<RealFunction<RealVect> > ptr1(new InclinedPlaneFunction(originElevation, basalSlope));
            bedrockFunction[0] =  ptr1;
            
            Real bumpCenter;
            Real bumpRad;
            Real bumpMag;
            
            mPP.get("bumpCenter", bumpCenter);
            mPP.get("bumpRad", bumpRad);
            mPP.get("bumpMag",bumpMag);
            RefCountedPtr<RealFunction<RealVect> > ptr2(new GaussianFunctionX(bumpCenter, bumpRad, bumpMag));

            bedrockFunction[1] =  ptr2;
          }
        else if (geometry == "Schoof")
          {
            //geometry of Schoof, 2007
            Real originElevation;
            mPP.get("originElevation", originElevation);
            
            Real lengthScaleFactor = 1.0;
            mPP.query("schoofLengthScaleFactor",lengthScaleFactor);
            
            
            Real schoofCoeff2, schoofCoeff4, schoofCoeff6;
            mPP.get("schoofCoeff2", schoofCoeff2);
            mPP.get("schoofCoeff4", schoofCoeff4);
            mPP.get("schoofCoeff6", schoofCoeff6);
            
            //RefCountedPtr<RealFunction<RealVect> > schoofBedrock
            RefCountedPtr<RealFunction<RealVect> > ptr(new SchoofBedrockElevation(domainSize[SpaceDim-2] * lengthScaleFactor,
                                                                                  originElevation,
                                                                                  schoofCoeff2, schoofCoeff4, 
                                                                                  schoofCoeff6));
            
            bedrockFunction[0] = ptr;
            
            
          }
        
        else if (geometry == "Katz")
          {
            //geometry of Katz and Worster, 2010
            Real originElevation;
            mPP.get("originElevation", originElevation);
            
            Real lengthScaleFactor = 1.0;
            mPP.query("schoofLengthScaleFactor",lengthScaleFactor);
            
            Real katzAlpha, katzSigma;
            mPP.get("katzAlpha", katzAlpha);
            mPP.get("katzSigma", katzSigma);
            
            Real schoofCoeff2, schoofCoeff4, schoofCoeff6;
            mPP.get("schoofCoeff2", schoofCoeff2);
            mPP.get("schoofCoeff4", schoofCoeff4);
            mPP.get("schoofCoeff6", schoofCoeff6);
            
            //RefCountedPtr<RealFunction<RealVect> > katzBedrock
            RefCountedPtr<RealFunction<RealVect> > ptr(new KatzBedrockElevation(domainSize[SpaceDim-2],
                                                                                domainSize[SpaceDim-1],
                                                                                originElevation,
                                                                                katzAlpha, katzSigma,
                                                                                lengthScaleFactor,
                                                                                schoofCoeff2, schoofCoeff4, 
                                                                                schoofCoeff6));
            
            bedrockFunction[0] = ptr;
            
            
          }
        else if (geometry == "MISMIPplus")
          {
	    RefCountedPtr<RealFunction<RealVect> > ptr(new MISOMIPBedrockElevation());
            bedrockFunction[0] = ptr;
          }
        else
          {
            MayDay::Error("bad marineIceSheet.geometry");
          }
        
        ibcPtr->setParameters(thicknessFunction, bedrockFunction ,  seaLevel);
        thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
      }
     else if (problem_type == "hump")
       {
         HumpIBC* ibcPtr = new HumpIBC;
         ParmParse humpPP("hump");

         Real maxThickness, radSqr, baseElevation, minThickness, seaLevel;
         RealVect center, widthScale;
         
         // default values to be equivalent to hump in Glimmer-CISM
         radSqr = 0.125*domainSize[0]*domainSize[1];
         maxThickness = 2000.0*pow(radSqr,0.5);
         baseElevation = 0.0;
         minThickness = 0.0;
         widthScale = RealVect::Unit;
         // this just lowers the sea level so that it's not relevant...
         seaLevel = -10.0;
         center = 0.5*domainSize;

         humpPP.query("radSqr", radSqr);
         humpPP.query("maxThickness", maxThickness);
         humpPP.query("baseElevation", baseElevation);
         humpPP.query("minThickness", minThickness);
         if (humpPP.contains("center"))
           {
             Vector<Real> centerArr(SpaceDim);
             humpPP.getarr("center", centerArr, 0, SpaceDim);
             center = RealVect(D_DECL(centerArr[0], centerArr[1], 
                                      centerArr[2]));
           }

         if (humpPP.contains("widthScale"))
           {
             Vector<Real> factorArr(SpaceDim);
             humpPP.getarr("widthScale", factorArr, 0, SpaceDim);
             widthScale = RealVect(D_DECL(factorArr[0], factorArr[1], 
                                           factorArr[2]));
           }

         ibcPtr->setParameters(maxThickness,
                               radSqr,
                               baseElevation,
                               minThickness,
                               center,
                               seaLevel, 
                               widthScale);
         
         thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
       }
     else if (problem_type == "LevelData")
       {
         //read geometry from an AMR Hierarchy, store in LevelDataIBC
         ParmParse ildPP("inputLevelData");
         std::string infile;
         ildPP.get("geometryFile",infile);
         std::string thicknessName = "thck";
         ildPP.query("thicknessName",thicknessName);
         std::string topographyName = "topg";
         ildPP.query("topographyName",topographyName);
         
         // default values (only relevant when LevelData input doesn't
         // cover the domain
         Real defaultThickness = 0.0;
         ildPP.query("defaultThickness", defaultThickness);
         Real defaultTopography = -10000.0;
         ildPP.query("defaultTopography", defaultTopography);
         bool setDefaultValues = false;
         ildPP.query("setDefaultValues", setDefaultValues);
         
         RefCountedPtr<LevelData<FArrayBox> > levelThck
           (new LevelData<FArrayBox>());
         RefCountedPtr<LevelData<FArrayBox> > levelTopg
           (new LevelData<FArrayBox>());
         
         Real dx;
         
         Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
         vectData.push_back(levelThck);
         vectData.push_back(levelTopg);
         
         Vector<std::string> names(2);
         names[0] = thicknessName;
         names[1] = topographyName;
         readLevelData(vectData,dx,infile,names,1);
         
         // this is about removing ice from regions which
         // don't affect the dynamics of the region, but which 
         // can cause our solvers problems. 
         // for now, store these regions in a separate file. 
         // There's probably a better way to do this.
         
         // this will contain the boxes in the index space of the 
         // original LevelData in which the thickness will be cleared.
         
         if (ildPP.contains("clearThicknessRegionsFile"))
           {
             Vector<Box> clearBoxes;
             std::string clearFile;
             ildPP.get("clearThicknessRegionsFile", clearFile);
             
             if (procID() == uniqueProc(SerialTask::compute))
               {
                 ifstream is(clearFile.c_str(), ios::in);
                 if (is.fail())
                   {
                     MayDay::Error("Cannot open file with regions for thickness clearing");
                   }
                 // format of file: number of boxes, then list of boxes.
                 int numRegions;
                 is >> numRegions;
                 
                 // advance pointer in file
                 while (is.get() != '\n');
                 
                 clearBoxes.resize(numRegions);
                 
                 for (int i=0; i<numRegions; i++)
                   {
                     Box bx;
                     is >> bx;
                     while (is.get() != '\n');
                     
                     clearBoxes[i] = bx;
                   }
                 
               } // end if serial proc
             // broadcast results
             broadcast(clearBoxes, uniqueProc(SerialTask::compute));
             
             // now loop over the thickness levelData and set intersections
             // with boxes to zero
             
             DataIterator dit = levelThck->dataIterator();
             for (dit.begin(); dit.ok(); ++dit)
               {
                 FArrayBox& thickFab = levelThck->operator[](dit);
                 const Box& fabBox = thickFab.box();
                 for (int boxno=0; boxno<clearBoxes.size(); boxno++)
                   {
                     Box intersectBox(fabBox);
                     intersectBox &= clearBoxes[boxno];
                     if (!intersectBox.isEmpty())
                       {
                         thickFab.setVal(0.0,intersectBox,0);
                       } // end if there's an intersection
                   } // end loop over clearboxes
               } // end loop over grids in thickness levelData

           } // end if we're setting thickness to zero
       
         RealVect levelDx = RealVect::Unit * dx;
         LevelDataIBC* ptr = new LevelDataIBC(levelThck,levelTopg,levelDx,
                                              defaultThickness,
                                              defaultTopography,
                                              setDefaultValues
                                              );
         thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
     else if (problem_type == "MultiLevelData")
       {
         //read geometry from an AMR Hierarchy, store in MultiLevelDataIBC
         ParmParse ildPP("inputLevelData");
         std::string infile;
         ildPP.get("geometryFile",infile);
         std::string thicknessName = "thck";
         ildPP.query("thicknessName",thicknessName);
         std::string topographyName = "topg";
         ildPP.query("topographyName",topographyName);
         
         
         Real dx;
         Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectData;
         
         
         Vector<std::string> names(2);
         names[0] = thicknessName;
         names[1] = topographyName;
         Vector<int> refRatio;
         readMultiLevelData(vectData,dx,refRatio,infile,names,1);
         
         RealVect crseDx = RealVect::Unit * dx;
         MultiLevelDataIBC* ptr = new MultiLevelDataIBC
           (vectData[0],vectData[1],crseDx,refRatio);
         thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
         
       }
#ifdef HAVE_PYTHON
     else if (problem_type == "Python")
       {
         
         ParmParse pyPP("PythonIBC");
         std::string module;
         pyPP.get("module",module);
         std::string thckFuncName = "thickness";
         pyPP.query("thicknessFunction",thckFuncName);
         std::string topgFuncName = "topography";
         pyPP.query("topographyFunction",topgFuncName);
         std::string rhsFuncName = "";
         pyPP.query("RHSFunction",rhsFuncName);
         std::string faceVelFuncName = "";
         pyPP.query("faceVelFunction",faceVelFuncName);
         PythonInterface::PythonIBC* ptr = new PythonInterface::PythonIBC
           (module, thckFuncName, topgFuncName, rhsFuncName,faceVelFuncName);
         thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
       }
#endif
     else 
       {
         MayDay::Error("bad problem type");
       }
    
    amrObjectFine.setThicknessBC(thicknessIBC);
    amrObjectCrse.setThicknessBC(thicknessIBC);    
    
    {
      // ---------------------------------------------
      // set surface heat boundary data 
      // ---------------------------------------------

      SurfaceFlux* surf_heat_boundary_data_ptr = SurfaceFlux::parse("surfaceHeatBoundaryData");
      ParmParse pps("surfaceHeatBoundaryData");
      bool diri = false; // flux boundary data by default
      pps.query("Dirichlett",diri);
      bool temp = true; //temperature boundary data by default
      pps.query("Temperature",temp);
      if (surf_heat_boundary_data_ptr == NULL)
  {
    if (diri)
      {
        const std::string err("If surfaceHeatBoundaryData.Dirichlett = true, surfaceHeatBoundaryData.type must be set");
        pout() << err << endl;
        MayDay::Error(err.c_str());
      }
    else
      {
        const std::string warn("No surfaceHeatBoundaryData.type specified, so zero flux set. Only relevant for amr.isothermal = false");
        pout() << warn << endl;
        // only warn if we're on processor 0, otherwise we get 
        // nproc copies of this warning to stderr
        if (procID() == uniqueProc(SerialTask::compute))
    {
      MayDay::Warning(warn.c_str());
    }
        surf_heat_boundary_data_ptr = new zeroFlux();
      }
  }

      amrObjectFine.setSurfaceHeatBoundaryData(surf_heat_boundary_data_ptr, diri, temp);
      amrObjectCrse.setSurfaceHeatBoundaryData(surf_heat_boundary_data_ptr, diri, temp);      
      if (surf_heat_boundary_data_ptr != NULL)
        {
          delete surf_heat_boundary_data_ptr;
          surf_heat_boundary_data_ptr=NULL;
        }
    
      // ---------------------------------------------
      // set basal (lower surface) heat boundary data. 
      // ---------------------------------------------
      
      SurfaceFlux* basal_heat_boundary_data_ptr = SurfaceFlux::parse("basalHeatBoundaryData");
      if (basal_heat_boundary_data_ptr == NULL)
        {
          basal_heat_boundary_data_ptr = new zeroFlux();
        }
      
      amrObjectFine.setBasalHeatBoundaryData(basal_heat_boundary_data_ptr);
      amrObjectCrse.setBasalHeatBoundaryData(basal_heat_boundary_data_ptr);

      if (basal_heat_boundary_data_ptr != NULL)
  {
    delete basal_heat_boundary_data_ptr;
    basal_heat_boundary_data_ptr=NULL;
  }
      
    }
    
    {
      IceInternalEnergyIBC* internalEnergyIBC = NULL;
      ParmParse tempPP("temperature");
      std::string tempType("constant");
      tempPP.query("type",tempType);
      if (tempType == "constant")
  {
    Real T = 258.0;
    tempPP.query("value",T);
    ConstantIceTemperatureIBC* ptr = new ConstantIceTemperatureIBC(T);
    internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
  }
      else if (tempType == "LevelData")
  {
    ParmParse ildPP("inputLevelData");
    LevelDataTemperatureIBC* ptr = NULL;
    ptr = LevelDataTemperatureIBC::parse(ildPP); CH_assert(ptr != NULL);
    internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
  }
      else if (tempType == "VerticalConduction")
  {
    VerticalConductionInternalEnergyIBC* ptr = NULL;
    ptr = VerticalConductionInternalEnergyIBC::parse(tempPP); 
    CH_assert(ptr != NULL);
    internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>(ptr);
  }
#ifdef HAVE_PYTHON
      else if (tempType == "Python")
  {
    ParmParse pyPP("PythonIceTemperatureIBC");
    std::string module;
    pyPP.get("module",module);
    std::string funcName = "temperature";
    pyPP.query("function",funcName);
    internalEnergyIBC  = static_cast<IceInternalEnergyIBC*>
      (new PythonInterface::PythonIceTemperatureIBC(module, funcName));
  }
#endif
      else 
  {
    MayDay::Error("bad temperature/internal energy type");
  } 
      amrObjectFine.setInternalEnergyBC(internalEnergyIBC);
      amrObjectCrse.setInternalEnergyBC(internalEnergyIBC);
      
      if (internalEnergyIBC != NULL)
  {
    delete internalEnergyIBC;
  }
    }

    amrObjectFine.setDomainSize(domainSize);
    amrObjectCrse.setDomainSize(domainSize);    


    CalvingModel* calving_model_ptr = CalvingModel::parseCalvingModel("CalvingModel");
    if (calving_model_ptr == NULL)
      {
  calving_model_ptr = new NoCalvingModel;
      }
    amrObjectFine.setCalvingModel(calving_model_ptr);
    amrObjectCrse.setCalvingModel(calving_model_ptr);    
    
    { 
      //// TODO this is just a temporary means of initializing a 
      //// damage model observer etc. Once we have it working, it 
      //// probably needs to be bound up with MuCoefficient
      bool damage_model = false;
      pp2.query("damage_model",damage_model);
      if (damage_model)
  {
    ///currently not deleting this, need to decide
    ///how that will be done
    DamageIceObserver* ptr = new DamageIceObserver();
    amrObjectFine.addObserver(ptr);
    amrObjectCrse.addObserver(ptr);          

    //whatever the constitutive relation was, wrap
    //it up in a DamageConstitutiveRelation tied
    //to the DamageIceObserver components
    DamageConstitutiveRelation* dcrptr = 
      new DamageConstitutiveRelation(constRelPtr, &ptr->damage());
    amrObjectFine.setConstitutiveRelation(dcrptr);
    amrObjectCrse.setConstitutiveRelation(dcrptr);          

    CalvingModel* d_calving_model_ptr = new DamageCalvingModel(calving_model_ptr, &ptr->damage());
    amrObjectFine.setCalvingModel(d_calving_model_ptr);
    amrObjectCrse.setCalvingModel(d_calving_model_ptr);
          
    delete d_calving_model_ptr;
    
  }
    }

    
    {
      /// initialize the melange model
      bool melange_model = false;
      pp2.query("melange_model",melange_model);
      if (melange_model)
  {
    MelangeIceObserver* ptr = new MelangeIceObserver();
    amrObjectFine.addObserver(ptr);
    amrObjectCrse.addObserver(ptr);          
  }
    }
    
    // set up initial grids, initialize data, etc.
    amrObjectFine.initialize();
    amrObjectCrse.initialize();    


    // now set up to do time integration, maintaining discrete states
    // outside of BISICLES itself
    
    // how many crse and fine time-intervals there are
    int numCrseIntervals, numFineIntervals;
    pp2.get("numCrseIntervals", numCrseIntervals);
    pp2.get("numFineIntervals", numFineIntervals);

    
    // now create an array of dataholders for ice thickness, dH/dt, and velocity
    // Each single-time field is a Vector<LevelData<FArrayBox>* >, so we
    // maintain Vectors of those.
    // each Vector has numIntervals + 1 members because we need an extra one for the
    // final solution.

   
    // Vectors of coarse and fine ice sheet states
    Vector<IceSheetState> crseStateVect(numCrseIntervals+1);
    Vector<IceSheetState> fineStateVect(numFineIntervals+1);

    // DFM (7/7/22) -- we might want to make these refCountedPtrs    
    Vector<Vector<LevelData<FArrayBox>* > > finedHdtVect(numFineIntervals+1);
    Vector<Vector<LevelData<FArrayBox>* > > crsedHdtVect(numCrseIntervals+1);
    

    // initialize each time interval from initial time
    // Vector<LevelData<FArrayBox>* > crseH;
    // amrObjectCrse.amrThickness(crseH);
    ParmParse ppcrse("crse"); int crselevel=0;
    ppcrse.get("amr.maxLevel", crselevel);
    int finestLevel = 0;
    Vector<LevelData<FArrayBox>* > crseH(crselevel+1, NULL);
    for (int lvl=0; lvl<crselevel+1;lvl++)
      {
        const DisjointBoxLayout& current_grid_size=amrObjectCrse.grids(lvl);
        // only do this if the level is defined
        if (current_grid_size.isClosed())
          {
            // keep track of the finest defined level
            finestLevel = lvl;
            crseH[lvl] = new LevelData<FArrayBox>(current_grid_size,1,IntVect::Zero); 
            DataIterator dit = (*crseH[lvl]).dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                (*crseH[lvl])[dit].setVal(1000.0);
              }
          }
      }

    const Vector<LevelData<FArrayBox>* >& crseVel = amrObjectCrse.amrVelocity();
    for (int i=0; i<numCrseIntervals; i++)
      {
        reshapeAndFill(crseStateVect[i].ice_thickness, crseH);
        reshapeAndFill(crseStateVect[i].ice_velocity,crseVel);
      }


    // initialize each time interval from initial time
    Vector<LevelData<FArrayBox>* > fineH;
    amrObjectFine.amrThickness(fineH);
    const Vector<LevelData<FArrayBox>* >& fineVel = amrObjectFine.amrVelocity();    
    for (int i=0; i<numFineIntervals; i++)
      {
        reshapeAndFill(fineStateVect[i].ice_thickness, fineH);
        reshapeAndFill(fineStateVect[i].ice_velocity,fineVel);
      }
    
    
    int maxStep;
    Real maxTime;
    //Real startTime;
    pp2.get("maxTime", maxTime);
    pp2.get("maxStep", maxStep);
    Real crseDt = maxTime/numCrseIntervals;
    Real fineDt = maxTime/numFineIntervals;


    // ----------------------------------------------- //
    // ------------ create pfasst object ------------- //
    // ----------------------------------------------- //
    // pfasst setup
    bool PF_VERBOSE;
    string pf_plot_prefix;
    int pf_num_procs_per_time, pf_num_repeats;
    ppfasst.get("PF_VERBOSE", PF_VERBOSE);
    ppfasst.get("pf_plot_prefix", pf_plot_prefix);
    ppfasst.get("pf_num_procs_per_time", pf_num_procs_per_time);
    ppfasst.get("pf_num_repeats", pf_num_repeats);

    ParmParse pcrse("crse.amr");
    string crse_plot_prefix;
    bool pf_evolve_velocity = true;
    pcrse.get("plot_prefix",crse_plot_prefix);
    pcrse.query("evolve_velocity",pf_evolve_velocity);

    if (USE_PF){
      pout() <<"###########################################################################"<<endl;
      pout() <<"##############        EL - Start of PFASST simulation     ###################"<<endl;
      pout() <<"###########################################################################"<<endl;
      // pout()<< "  dt passed in: "<<crseDt<<", max T passed in: "<<maxTime<<", max steps passed in:"<<maxStep<< endl;

      AmrIceHolderClass AmrIceHolderPtr;
      MPI_Fint pf_comm_world = MPI_Comm_c2f(pf_comm);
      
      Pf_Bisicles_setHolders(&amrObjectCrse,&AmrIceHolderPtr,crseH,crseVel); // pass amrObjectCrse to AmrIceHolderPt, pretty sure is right
      
      
      // set up vector size (including all amr levels), i.e. total num of cells in all levels
      ParmParse ppcrseamr("crse.amr");
      Vector<int> ancells(3); 
      ppcrseamr.getarr("num_cells", ancells, 0, ancells.size());\
      int num_of_grids=ancells[0]*ancells[1];
      int num_total_cells=0;
      for (int lvl=0; lvl < crseH.size(); lvl++)
      {
        if (crseH[lvl] != NULL)
          {
            LevelData<FArrayBox>& ldf = *crseH[lvl];
            DisjointBoxLayout dbl = ldf.disjointBoxLayout();
            if (dbl.isClosed())
              {
                DataIterator dit = ldf.dataIterator();
                int num_cells_per_lvl = dbl.numCells();
                num_total_cells += num_cells_per_lvl;
              }
          }
      }
      reshape(crsedHdtVect[0],crseH);
      

      // call pfasst to solve
      Vector<LevelData<FArrayBox>* > ice_thick=crseStateVect[1].ice_thickness;
      auto start_pf = std::chrono::high_resolution_clock::now();
      Pf_Main(&AmrIceHolderPtr,pf_comm_world,numCrseIntervals,crseDt,maxTime,maxStep,pf_evolve_velocity,\
              num_total_cells,pf_num_repeats,PF_VERBOSE);
      auto end_pf = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration_pf = end_pf - start_pf;
      (&AmrIceHolderPtr)->printStatistics(number_procs,new_time_size, new_space_size,new_time_rank,\
                                          new_space_rank,duration_pf,pf_evolve_velocity,numCrseIntervals,maxTime);
      cout <<"  EL - results saved as: "<<crse_plot_prefix<<"...\n";
      cout <<"###########################################################################"<<endl;
      cout <<"##############        EL - End of PFASST simulation     ###################"<<endl;
      cout <<"###########################################################################"<<endl;
    } else {
    cout<<"\nUpdating crse-grained using serial................\n";
    cout<<"  results saved as: "<<crse_plot_prefix<<"...\n";
    cout<<"Not running coarse sequential\n";
    }


    // ------- EL - commenting out the serial run for now -----------
    // ParmParse pfine("fine.amr");
    // string fine_plot_prefix;
    // pfine.get("plot_prefix",fine_plot_prefix);
    // cout<<"\n..............Updating fine-grained using serial................\n";
    // cout<<"  results saved as: "<<fine_plot_prefix<<"...\n\n\n";

    // now do each fine-grained timestep
    // auto start = high_resolution_clock::now();
    // amrObjectFine.run(maxTime, maxStep);
    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);
    // cout << "Time taken by serial run: "
    //      << duration.count()/1e6<< " seconds" << endl;
    // don's really need the time interval loop, the loop is for illustration
    // for (int i=0; i<numFineIntervals; i++) 
    //   {
    //     // first, set state
    //     bool recalculateVelocity = true;
    //     //amrObjectFine.setState(fineHVect[i], fineTimeVect[i], recalculateVelocity);
    //     recalculateVelocity = false;
    //     amrObjectFine.setState(fineStateVect[i],recalculateVelocity); 

    //     // reshape dH/dt and then call computeDhDt
    //     reshape(finedHdtVect[i],fineH); // fineH
    //     amrObjectFine.compute_dHdt(finedHdtVect[i],fineH,fineDt, false);    // fineH     

    //     // now advance ice sheet
    //     cout<<"time "<<fineStateVect[i+1].time<<endl;
    //     // amrObjectFine.run(maxTime, maxStep);

    //     // now retrieve state and store
    //     amrObjectFine.getState(fineStateVect[i+1]);
    //     cout << "loop # " << i << std::endl;

      // }    

        

    
    // clean up
    if (constRelPtr != NULL)
      {
        delete constRelPtr;
        constRelPtr = NULL;
      }

    if (surf_flux_ptr != NULL)
      {
        delete surf_flux_ptr;
        surf_flux_ptr = NULL;
      }
    
    if (basal_flux_ptr != NULL)
      {
        delete basal_flux_ptr;
        basal_flux_ptr = NULL;
      }

    if (topg_flux_ptr != NULL)
      {
        delete topg_flux_ptr;
        topg_flux_ptr = NULL;
      }

    if (calving_model_ptr != NULL)
      {
  delete calving_model_ptr;
  calving_model_ptr = NULL;
      }
    
    if (basalFrictionPtr != NULL)
      {
  delete basalFrictionPtr;
  basalFrictionPtr = NULL;
      }

    if (basalFrictionRelationPtr != NULL)
      {
  delete basalFrictionRelationPtr;
  basalFrictionRelationPtr = NULL;
      }

    if (thicknessIBC != NULL)
      {
        delete thicknessIBC;
        thicknessIBC=NULL;
      }

     
    
#ifdef CH_USE_HDF5
    {
      // finally, carry out an optional regression test
      
      ParmParse ppr("regression");
      
      std::string result_hdf5("");
      ppr.query("result_hdf5", result_hdf5);

      if (result_hdf5 != "")
  {
    std::string reference_hdf5;
    ppr.get("reference_hdf5", reference_hdf5);

    Real tol = 1.0e-10;
    ppr.query("tol",tol);
          if (tol < HDF5NormTest(result_hdf5, reference_hdf5))
      {
        ierr = 1;
      }
  }

    }
#endif

  }  // end nested scope
  

  CH_TIMER_REPORT();

#ifdef HAVE_PYTHON
  Py_Finalize();
#endif

#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional

  return ierr;
}


