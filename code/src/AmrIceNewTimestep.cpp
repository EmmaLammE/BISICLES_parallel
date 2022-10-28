#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

using std::ifstream; 
using std::ios;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::string;
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
#include "CellToEdge.H"
#include "CH_HDF5.H"
#include "IceUtility.H"
#include "LevelMappedDerivatives.H"
#include "SundialsUtil.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

#include "NamespaceHeader.H"

void
AmrIce::newTimeStep(Real a_dt)
{

  CH_TIME("AmrIce::newTimestep");
    
  if (s_verbosity >=2) 
    {
      pout() << "New Timestep " << m_cur_step 
             << " Advancing solution from time " 
             << m_time << " ( " << time() << ")" " with dt = " << a_dt << endl;
    }

  // set member dt to be dt from argument 
  m_dt = a_dt;

  
  // first copy thickness into old thickness   
  for (int lev=0; lev <= m_finest_level ; lev++)
    {      
      LevelData<FArrayBox>& oldThickness = *m_old_thickness[lev];
      LevelData<FArrayBox>& currentThickness = m_vect_coordSys[lev]->getH();

      // this way we avoid communication and maintain ghost cells...
      DataIterator dit = oldThickness.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          oldThickness[dit].copy(currentThickness[dit],0, 0, 1);
        }

    }
        
  // assumption here is that we've already computed the current velocity 
  // field, most likely at initialization or at the end of the last timestep...
  // so, we don't need to recompute the velocity at the start.

  // use PatchGodunov hyperbolic solver
  

  // make a copy of m_vect_coordSys before it is overwritten
  // DFM -- toDo: move into a separate function to make this cleaner
  Vector<RefCountedPtr<LevelSigmaCS> > vectCoords_old (m_finest_level+1);
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      IntVect sigmaCSGhost = m_vect_coordSys[lev]->ghostVect();
      RealVect dx = m_amrDx[lev]*RealVect::Unit;
      vectCoords_old[lev] = RefCountedPtr<LevelSigmaCS> 
        (new LevelSigmaCS(m_amrGrids[lev], dx, sigmaCSGhost));
      LevelSigmaCS& levelCoords_old = *vectCoords_old[lev];
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      
      
      levelCoords_old.setIceDensity(levelCoords.iceDensity());
      levelCoords_old.setWaterDensity(levelCoords.waterDensity());
      levelCoords_old.setGravity(levelCoords.gravity());
      // todo replace the copies below with a deepCopy of levelCoords
      for (DataIterator dit( m_amrGrids[lev]); dit.ok(); ++dit)
        {
          FArrayBox& oldH = levelCoords_old.getH()[dit];
          const FArrayBox& H = levelCoords.getH()[dit];
          oldH.copy(H);
        }
      levelCoords_old.setTopography(levelCoords.getTopography());
      {
        LevelSigmaCS* crseCoords = (lev > 0)?&(*vectCoords_old[lev-1]):NULL;
        int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
        levelCoords_old.recomputeGeometry( crseCoords, refRatio);
      }
#if BISICLES_Z == BISICLES_LAYERED
      levelCoords_old.setFaceSigma(levelCoords.getFaceSigma());
#endif
    }



  
  // holder for dH/dt
  Vector<LevelData<FArrayBox>* > dHdt(m_finest_level+1, NULL);

  // local scratch space for new thickness
  Vector<LevelData<FArrayBox>* > newThickness(m_finest_level+1, NULL);

  
  // allocate storage
  for (int lev = finestTimestepLevel() ; lev>=0 ; lev--)
    {
      
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

      IntVect ghostVect = IntVect::Unit;      
      dHdt[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, 
                                           ghostVect);


      newThickness[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, 
                                           ghostVect);      


      LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];
            
      // ensure that ghost cells for thickness  are filled in
      if (lev > 0)
        {          
          int nGhost = levelOldThickness.ghostVect()[0];
          PiecewiseLinearFillPatch thicknessFiller(levelGrids, 
                                                   m_amrGrids[lev-1],
                                                   1, 
                                                   m_amrDomains[lev-1],
                                                   m_refinement_ratios[lev-1],
                                                   nGhost);
          
          // since we're not subcycling, don't need to interpolate in time
          Real time_interp_coeff = 0.0;
          thicknessFiller.fillInterp(levelOldThickness,
                                     *m_old_thickness[lev-1],
                                     *m_old_thickness[lev-1],
                                     time_interp_coeff,
                                     0, 0, 1);
          
          
          
        }
      // these are probably unnecessary...
      levelOldThickness.exchange();
      
    }
  
  // compute estimate of dH/dt (H) at current time
  bool recomputeVelocity = true;

  compute_dHdt(dHdt, m_old_thickness,
               m_dt/2, recomputeVelocity);
  
  
  // ignore temperarature, iceFraction updates and diagnostic output for now

  // for now, simple forward-euler update
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& levelGrids = grids(lev);
    LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];
    LevelData<FArrayBox>& levelNewThickness = *newThickness[lev];
    LevelData<FArrayBox>& leveldHdt = *dHdt[lev];

    leveldHdt.copyTo(levelNewThickness);
    
    DataIterator dit = levelGrids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        // scale dH/dt by dt
        levelNewThickness[dit].mult(m_dt / 2);
        // add in old H
        levelNewThickness[dit].plus(levelOldThickness[dit]);
      }
  } // end loop over levels      

  //  update geometry
  // dt to use in update
  Real updateDt = 0.5*m_dt;
  updateGeometryFromThickness(m_vect_coordSys, vectCoords_old,
                              newThickness, updateDt);       

  compute_dHdt(dHdt, newThickness,
               m_dt/2, recomputeVelocity);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& levelGrids = grids(lev);
    LevelData<FArrayBox>& levelOldThickness = *m_old_thickness[lev];
    LevelData<FArrayBox>& levelNewThickness = *newThickness[lev];
    LevelData<FArrayBox>& leveldHdt = *dHdt[lev];

    leveldHdt.copyTo(levelNewThickness);
    
    DataIterator dit = levelGrids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        // scale dH/dt by dt
        levelNewThickness[dit].mult(m_dt);
        // add in old H
        levelNewThickness[dit].plus(levelOldThickness[dit]);
      }
  } // end loop over levels      


  // update ice fraction through advection
  //advectIceFrac(m_iceFrac, m_volumeThicknessSource, m_faceVelAdvection,vectFluxes, a_dt);
  if (m_evolve_ice_frac)
    {
      advectIceFrac(m_iceFrac, m_faceVelAdvection, a_dt);
    }
  
  

  // note that we use "updateDt" here
  // since we are updating the geometry from t+dt/2 -> t+dt
  updateGeometryFromThickness(m_vect_coordSys, vectCoords_old,
                              newThickness, updateDt);       

  
  notifyObservers(Observer::PostGeometryUpdate);
  
  // clean up temporary storage
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      if (dHdt[lev] != NULL)
        {          
          delete dHdt[lev];
          dHdt[lev] = NULL;
        }
                       
      if (newThickness[lev] != NULL)
        {
          delete newThickness[lev];
          newThickness[lev] = NULL;
        }
 
    }
  
 
  //update time (velocity is to be computed at the step end)
  m_time += a_dt;
  m_cur_step += 1;

  // compute new ice velocity field
  /*if (m_evolve_velocity )
    {
      if (s_verbosity > 3) 
	{
	  pout() << "AmrIce::timeStep solveVelocityField() (step end) " << endl;
	}
      solveVelocityField();
    }
  
  if (s_verbosity > 0) 
    {
      pout () << "AmrIce::timestep " << m_cur_step
              << " --     end time = " 
	      << setiosflags(ios::fixed) << setprecision(6) << setw(12)
              << m_time  << " ( " << time() << " )"
        //<< " (" << m_time/secondsperyear << " yr)"
              << ", dt = " 
        //<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
              << a_dt
        //<< " ( " << a_dt/secondsperyear << " yr )"
	      << resetiosflags(ios::fixed)
              << endl;
    }*/

  
  int totalCellsAdvanced = 0;
  for (int lev=0; lev<m_num_cells.size(); lev++) 
    {
      totalCellsAdvanced += m_num_cells[lev];
    }
     
  if (s_verbosity > 0) 
    {
      pout() << "Time = " << m_time  
             << " cells advanced = " 
             << totalCellsAdvanced << endl;

      for (int lev=0; lev<m_num_cells.size(); lev++) 
        {
          pout () << "Time = " << m_time 
                  << "  level " << lev << " cells advanced = " 
                  << m_num_cells[lev] << endl;
        }
    }
}


// compute estimate of dH/dt at the current time
void
AmrIce::compute_dHdt(Vector<LevelData<FArrayBox>* >& a_dHdt,
                     Vector<LevelData<FArrayBox>* >& a_H, 
                     Real a_dt, bool a_recomputeVelocity)
{
  if (a_recomputeVelocity)
    {
      solveVelocityField();
    }

  // volumeThicknessSource is used in the semi-implicit scheme
  setToZero(m_volumeThicknessSource);
  
  Vector<LevelData<FluxBox>* > faceH(a_dHdt.size(), NULL);
  Vector<LevelData<FluxBox>* > vectFluxes(a_dHdt.size(), NULL);  
  for (int lev=0; lev<faceH.size(); lev++)
    {

      // allocate storage

      const LevelData<FArrayBox>& currentH = *a_H[lev];
      const DisjointBoxLayout& levelGrids = currentH.getBoxes();      
      faceH[lev] = new LevelData<FluxBox>(levelGrids, 1, IntVect::Zero);
      // vectFluxes needs to have a ghost cell to ensure that
      // the Copier works correctly (it's a bit complicated)
      vectFluxes[lev] = new LevelData<FluxBox>(levelGrids, 1, IntVect::Unit);

      
      LevelData<FArrayBox>& levelSTS = *m_surfaceThicknessSource[lev];
      LevelData<FArrayBox>& levelBTS = *m_basalThicknessSource[lev];
      LevelData<FArrayBox>& levelVTS = *m_volumeThicknessSource[lev];
      CH_assert(m_surfaceFluxPtr != NULL);
      
      // set surface thickness source
      m_surfaceFluxPtr->surfaceThicknessFlux(levelSTS, *this, lev, a_dt);
      
      // set basal thickness source
      m_basalFluxPtr->surfaceThicknessFlux(levelBTS, *this, lev, a_dt);
           
      // for now, just do simple cell-to-face averaging
      // (may want to go back to PPM at some point)
      CellToEdge(currentH,*faceH[lev]);
    }
  // dhDt is just div(faceH*faceVel) + thicknessSources

  // we already have a function to compute thicknessFluxes (faceH*faceVel)
  computeThicknessFluxes(vectFluxes, faceH, m_faceVelAdvection);

  // now compute divergence
  for (int lev=0; lev <= finestTimestepLevel() ; lev++)
    {
      DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FluxBox>& levelFlux = *vectFluxes[lev];
      LevelData<FArrayBox>& levelDhDt = *a_dHdt[lev];
      LevelData<FArrayBox>& levelDivThckFlux = *m_divThicknessFlux[lev];
      const LevelSigmaCS& levelCoords = *m_vect_coordSys[lev];
      const RealVect& dx = levelCoords.dx();              
      
      DataIterator dit = levelGrids.dataIterator();          

      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          FArrayBox& dHdt = levelDhDt[dit];
          FluxBox& thisFlux = levelFlux[dit];
          dHdt.setVal(0.0);
          
          // loop over directions and increment with div(F)
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // use the divergence from 
              // Chombo/example/fourthOrderMappedGrids/util/DivergenceF.ChF
              FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                              CHF_FRA(dHdt),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(dx[dir]),
                              CHF_INT(dir));
              
              
            } // end loop over directions
          
	  levelDivThckFlux[dit].copy(dHdt);


          dHdt *= -1.0;

          // // modify sources as required to balance div(uh)
          setStableSources((*m_surfaceThicknessSource[lev])[dit],
                           (*m_basalThicknessSource[lev])[dit],
                           (*m_volumeThicknessSource[lev])[dit],
                           levelDivThckFlux[dit],
                           levelCoords.getFloatingMask()[dit],
                           gridBox);


          // add in thickness source
          // if there are still diffusive fluxes to deal
          // with, the source term will be included then
          
	  //if (m_diffusionTreatment != IMPLICIT)
          
          if (m_frac_sources)
            {
              // scale surface fluxes by mask values
              const FArrayBox& thisFrac = (*m_iceFrac[lev])[dit];
              FArrayBox sources(gridBox,1);
              sources.setVal(0.0);
              sources.plus((*m_surfaceThicknessSource[lev])[dit], gridBox,
                           0, 0, 1);
              sources.plus((*m_basalThicknessSource[lev])[dit], gridBox, 
                           0, 0, 1);
              sources.plus((*m_volumeThicknessSource[lev])[dit], gridBox, 
                           0, 0, 1);
              
              sources.mult(thisFrac, gridBox, 0, 0, 1);
              dHdt.minus(sources, gridBox, 0, 0, 1);
              
            }
          else 
            {
              // just add in sources directly
              dHdt.plus((*m_surfaceThicknessSource[lev])[dit], gridBox,0,0,1);
              dHdt.plus((*m_basalThicknessSource[lev])[dit], gridBox,0,0,1);
              dHdt.plus((*m_volumeThicknessSource[lev])[dit], gridBox,0,0,1);
            }
      // pfasst test
      // cout << "DataIterator dit   " << dit << ", dHdt" << dHdt.getBoxes() << std::endl;

        } // end loop over boxes
    } // end loop over levels
      

  // clean up local storage
  for (int lev=0; lev<vectFluxes.size(); lev++)
    {
      if (vectFluxes[lev] != NULL)
        {
          delete vectFluxes[lev];
          vectFluxes[lev] = NULL;
        }
      if (faceH[lev] != NULL)
        {
          delete faceH[lev];
          faceH[lev] = NULL;
        }
    } // end loop over levels

  
}



/// set ice sheet state based on inputs
/**  New interface function for SUNDIALS and Parallel-in-time
     Take provided ice thickness and current time, reset ice sheet to match
     Recalculate velocity if required.
     TODO: What about cur_step?
     TODO: What about iceFrac?
     TODO: What about thermodynamics?
*/
void
AmrIce::setState(Vector<LevelData<FArrayBox>* >& a_thicknessVect,
                 Real a_cur_time, bool a_recalculateVelocity)
{
  // to maximize code reuse, pass an empty velocity vector in to the
  // oher setState 
  Vector<LevelData<FArrayBox>* > velocity;
  setState(a_thicknessVect, velocity, a_cur_time, a_recalculateVelocity);
}


/// set ice sheet state based on inputs
/**  New interface function for SUNDIALS and Parallel-in-time
     Take provided ice thickness and current time, reset ice sheet to match
     This version takes a provided velocity field as well. 
     TODO: What about cur_step?
*/
void
AmrIce::setState(Vector<LevelData<FArrayBox>* >& a_thicknessVect,
                 Vector<LevelData<FArrayBox>* >& a_velocityVect,
                 Real a_cur_time, bool a_recalculateVelocity)
{
  
  // this looks an awful lot like AmrIce::regrid, minus the part where actually tag for regrid
  // TODO:  Make a function which shares the common bits.

  CH_TIME("AmrIce::setState");

  if (s_verbosity > 3)
    {
      pout() << "AmrIce::setState" << endl;
    }


  // first, reset time
  if (s_verbosity > 3)
    {
      pout() << "resetting solution time from " << m_time
             << "  to " << a_cur_time << endl;
    }
  m_time = a_cur_time;
  
  // it's convenient to pull out the grids here
  Vector<DisjointBoxLayout> new_amrGrids(a_thicknessVect.size());
  for (int lev=0; lev<new_amrGrids.size(); lev++)
    {
      new_amrGrids[lev] = a_thicknessVect[lev]->getBoxes();
    }


  /// begin code lifted from AmrIce::regrid
  //test to see if grids have changed
  bool gridsSame = true;
  // assume that we'll always have the same level 0
  int lbase = 0;
  for (int lev=lbase+1; lev< new_amrGrids.size(); ++lev)
    {
      const DisjointBoxLayout oldDBL = m_amrGrids[lev];
      gridsSame &= oldDBL.sameBoxes(new_amrGrids[lev]);
    }


  // first do level 0 because it's the simplest (and because
  // we may wind up interpolating from it)
  
  /// overwrite thickness...
  LevelData<FArrayBox>& thisLevelH = m_vect_coordSys[0]->getH();
  
  DataIterator dit=new_amrGrids[0].dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      thisLevelH[dit].copy((*a_thicknessVect[0])[dit]);
    }

  /// recompute geometry
  {
    LevelSigmaCS* crseCoords = NULL;
    int refRatio = -1;
    m_vect_coordSys[0]->recomputeGeometry(crseCoords,refRatio);
  }
  
  /// overwrite velocity if needed
  if (a_velocityVect.size() > 0)
    {
      CH_assert(m_velocity[0]->getBoxes() == a_velocityVect[0]->getBoxes());

      
      DataIterator dit = m_velocity[0]->dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*m_velocity[0])[dit].copy((*a_velocityVect[0])[dit]);
        }      
    }
  
  // now loop through levels and redefine if necessary
  // this is more complicated than level 0 because we might have changed
  // grids
  
  for (int lev=lbase+1; lev< new_amrGrids.size(); ++lev)
    {
      const DisjointBoxLayout oldDBL = m_amrGrids[lev];
      
      m_amrGrids[lev] = new_amrGrids[lev];

      const DisjointBoxLayout newDBL = m_amrGrids[lev];      
      
      // first we need to regrid m_deltaTopography, it will be needed to 
      // regrid the bedrock topography & hence LevelSigmaCS
      m_deltaTopography[lev] = destructiveRegrid(m_deltaTopography[lev], newDBL,
                                                 m_deltaTopography[lev-1], m_refinement_ratios[lev-1]) ;
      
      // LevelSigmaCS regrid 
      {
        IntVect sigmaCSGhost = m_vect_coordSys[0]->ghostVect();
        RealVect dx = m_amrDx[lev]*RealVect::Unit;
        RefCountedPtr<LevelSigmaCS > oldCoordSys = m_vect_coordSys[lev];
        RefCountedPtr<LevelSigmaCS > auxCoordSys = (lev > 0)?m_vect_coordSys[lev-1]:oldCoordSys;
        
        m_vect_coordSys[lev] = RefCountedPtr<LevelSigmaCS >
          (new LevelSigmaCS(new_amrGrids[lev], dx, sigmaCSGhost));
        m_vect_coordSys[lev]->setIceDensity(auxCoordSys->iceDensity());
        m_vect_coordSys[lev]->setWaterDensity(auxCoordSys->waterDensity());
        m_vect_coordSys[lev]->setGravity(auxCoordSys->gravity());
        m_vect_coordSys[lev]->setBackgroundSlope(auxCoordSys->getBackgroundSlope());
#if BISICLES_Z == BISICLES_LAYERED
        m_vect_coordSys[lev]->setFaceSigma(auxCoordSys->getFaceSigma());
#endif		
        LevelSigmaCS* crsePtr = &(*m_vect_coordSys[lev-1]);
        int refRatio = m_refinement_ratios[lev-1];
        
        bool interpolate_zb = (m_interpolate_zb ||
                               !m_thicknessIBCPtr->regridIceGeometry
                               (*m_vect_coordSys[lev],dx,  m_domainSize, 
                                m_time,  crsePtr,refRatio ) );
        
        if (!interpolate_zb)
          {
            // need to re-apply accumulated bedrock (GIA). Could be optional?
            for (DataIterator dit(newDBL); dit.ok(); ++dit)
              {
                m_vect_coordSys[lev]->getTopography()[dit] += (*m_deltaTopography[lev])[dit];
              }
          }
        
        {
          //interpolate thickness & (maybe) topography
          // (will throw this thickness field away in the next step)
          
          bool interpolateThickness(true);
          bool preserveMask(true);
          bool interpolateTopographyGhost(true); 
          bool interpolateThicknessGhost(true); 
          bool preserveMaskGhost(true);
          m_vect_coordSys[lev]->interpFromCoarse(*m_vect_coordSys[lev-1],
                                                 m_refinement_ratios[lev-1],
                                                 interpolate_zb,
                                                 interpolateThickness, 
                                                 preserveMask,
                                                 interpolateTopographyGhost, 
                                                 interpolateThicknessGhost, 
                                                 preserveMaskGhost, 
                                                 m_regrid_thickness_interpolation_method);
        }


        LevelData<FArrayBox>& thisLevelH = m_vect_coordSys[lev]->getH();
        LevelData<FArrayBox>& thisLevelB = m_vect_coordSys[lev]->getTopography();
		
        // overwrite interpolated fields in valid regions with such valid old data as there is
        // ** this is where we replace interpolated old thickness with one that
        // we've passed in.  Do this fab-by-fab to maintain ghost cells
        DataIterator dit=new_amrGrids[lev].dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            thisLevelH[dit].copy((*a_thicknessVect[lev])[dit]);
          }
        
        if (oldDBL.isClosed())
          {	              
            const LevelData<FArrayBox>& oldLevelB = oldCoordSys->getTopography();
            oldLevelB.copyTo(thisLevelB);
          }
        
        //Defer to m_thicknessIBCPtr for boundary values - 
        //interpolation won't cut the mustard because it only fills
        //ghost cells overlying the valid regions.
        RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
        m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
                                          m_amrDomains[lev],levelDx, m_time, m_dt);
        
        
        m_vect_coordSys[lev]->exchangeTopography();
        
        {
          LevelSigmaCS* crseCoords = (lev > 0)?&(*m_vect_coordSys[lev-1]):NULL;
          int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:-1;
          m_vect_coordSys[lev]->recomputeGeometry(crseCoords,refRatio);
        }
      }

      // regrid other prognostic fields
      // (DFM 7/15/22) -- we may want to pass these in
      m_old_thickness[lev] = destructiveRegrid(m_old_thickness[lev], newDBL, m_old_thickness[lev-1], m_refinement_ratios[lev-1]) ;

      ///  DFM (7/15/22) -- will want to pass this in once we're doing marine ice sheets
      m_iceFrac[lev] = destructiveRegrid( m_iceFrac[lev], newDBL, m_iceFrac[lev-1],	m_refinement_ratios[lev-1]);
      // if we're not passing in a velocity, treat this like a regular regrid
      if (a_velocityVect.size() == 0)
        {
          m_velocity[lev] = destructiveRegrid( m_velocity[lev], newDBL,  m_velocity[lev-1], m_refinement_ratios[lev-1]);
        }
      else
        {
          // reshape and copy velocity. It is an error here if the velocity and thickness
          // aren't defined on the same grids
          CH_assert(a_velocityVect[lev]->getBoxes() == a_thicknessVect[lev]->getBoxes());

          LevelData<FArrayBox>& newVel = *a_velocityVect[lev];
          delete m_velocity[lev];
          m_velocity[lev] = new LevelData<FArrayBox>(new_amrGrids[lev], newVel.nComp(),
                                                     newVel.ghostVect());
          DataIterator dit = m_velocity[lev]->dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*m_velocity[lev])[dit].copy(newVel[dit]);
            }
        } // end if we're passing in a velocity
          
      // DFM (7/15/22) -- I don't think this is relevant in this context...
      // if we're passing in a velocity field, they should be properly set, and
      // if not, the velocity field doesn't make much sense before it's recomputed,
      // so no need to set ghost cells
#if 0
      {
        //handle ghost cells on the coarse-fine interface
        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                          m_amrDx[lev], m_refinement_ratios[lev-1],
                          2, m_amrDomains[lev]);
        qcfi.coarseFineInterp(*m_velocity[lev], *m_velocity[lev-1]);
        
        //boundary ghost cells
        m_thicknessIBCPtr->velocityGhostBC
          (*m_velocity[lev], *m_vect_coordSys[lev],m_amrDomains[lev],m_time);
      }
#endif // end vestigial velocity ghost-cell stuff

      
      //calved ice  regrid
      m_calvedIceThickness[lev] = destructiveRegrid(m_calvedIceThickness[lev], newDBL, m_calvedIceThickness[lev-1],  m_refinement_ratios[lev-1]);
      m_removedIceThickness[lev] = destructiveRegrid(m_removedIceThickness[lev], newDBL, m_removedIceThickness[lev-1],  m_refinement_ratios[lev-1]);
      m_addedIceThickness[lev] = destructiveRegrid(m_addedIceThickness[lev], newDBL, m_addedIceThickness[lev-1],  m_refinement_ratios[lev-1]);
      
      //internal energy regrid
      m_internalEnergy[lev] = destructiveRegrid( m_internalEnergy[lev], newDBL, m_internalEnergy[lev-1], m_refinement_ratios[lev-1]);
      m_tillWaterDepth[lev] = destructiveRegrid(m_tillWaterDepth[lev], newDBL, m_tillWaterDepth[lev-1], m_refinement_ratios[lev-1]);
      m_sInternalEnergy[lev] = destructiveRegrid( m_sInternalEnergy[lev], newDBL, m_sInternalEnergy[lev-1], m_refinement_ratios[lev-1]);
      m_bInternalEnergy[lev] = destructiveRegrid( m_bInternalEnergy[lev], newDBL, m_bInternalEnergy[lev-1], m_refinement_ratios[lev-1]);
      //internal energy boundary ghost cells
      m_internalEnergyIBCPtr->setIceInternalEnergyBC
        (* m_internalEnergy[lev], *m_tillWaterDepth[lev], *m_sInternalEnergy[lev], *m_bInternalEnergy[lev],
         *m_vect_coordSys[lev] );
      //since the internalEnergy data has changed
      m_A_valid = false;
      
      // no need to regrid, just reallocate
      if (!gridsSame)
        {
          if (m_velBasalC[lev] != NULL)
            {
              delete m_velBasalC[lev];
            }
          m_velBasalC[lev] = new LevelData<FArrayBox>(newDBL, 1, IntVect::Unit);
          
          if (m_cellMuCoef[lev] != NULL)
            {
              delete m_cellMuCoef[lev];
            }
          m_cellMuCoef[lev] = new LevelData<FArrayBox>(newDBL, 1, IntVect::Unit);
          
          if (m_velRHS[lev] != NULL)
            {
              delete m_velRHS[lev];
            }
          m_velRHS[lev] = new LevelData<FArrayBox>(newDBL, SpaceDim, 
                                                   IntVect::Zero);
          
          if (m_faceVelAdvection[lev] != NULL)
            {
              delete m_faceVelAdvection[lev];
            }
          m_faceVelAdvection[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);
          
          if (m_faceVelTotal[lev] != NULL)
            {
              delete m_faceVelTotal[lev];
            }
          m_faceVelTotal[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);
          
          
          if (m_diffusivity[lev] != NULL)
            {
              delete m_diffusivity[lev];
            }
          m_diffusivity[lev] = new LevelData<FluxBox>(newDBL, 1, IntVect::Unit);
          
          
          if (m_surfaceThicknessSource[lev] != NULL)
            {
              delete m_surfaceThicknessSource[lev];
            }
          m_surfaceThicknessSource[lev] = 
            new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit) ;
          
          if (m_basalThicknessSource[lev] != NULL)
            {
              delete m_basalThicknessSource[lev];
            }
          m_basalThicknessSource[lev] = 
            new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit) ;
          
          if (m_volumeThicknessSource[lev] != NULL)
            {
              delete m_volumeThicknessSource[lev];
            }
          m_volumeThicknessSource[lev] = 
            new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit) ;
          
          if (m_divThicknessFlux[lev] != NULL)
            {
              delete m_divThicknessFlux[lev];
            }
          m_divThicknessFlux[lev] = 
            new LevelData<FArrayBox>(newDBL,   1, IntVect::Zero) ;
          
          
          if (m_bHeatFlux[lev] != NULL)
            {
              delete m_bHeatFlux[lev];
            }
          m_bHeatFlux[lev] = 
            new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit);
          
          if (m_sHeatFlux[lev] != NULL)
            {
              delete m_sHeatFlux[lev];
            }
          m_sHeatFlux[lev] = 
            new LevelData<FArrayBox>(newDBL,   1, IntVect::Unit);
          
          if (m_layerXYFaceXYVel[lev] != NULL)
            {
              delete m_layerXYFaceXYVel[lev];
            }
          
          m_layerXYFaceXYVel[lev] = new LevelData<FluxBox>
            (newDBL, m_nLayers, IntVect::Unit);
          
          if (m_layerSFaceXYVel[lev] != NULL)
            {
              delete m_layerSFaceXYVel[lev];
            }
          
          m_layerSFaceXYVel[lev] = new LevelData<FArrayBox>
            (newDBL, SpaceDim*(m_nLayers + 1), IntVect::Unit);
          
          if (m_layerSFaceSVel[lev] != NULL)
            {
              delete m_layerSFaceSVel[lev];
            }
          
          m_layerSFaceSVel[lev] = new LevelData<FArrayBox>
            (newDBL, SpaceDim*(m_nLayers + 1), IntVect::Unit);
          
	} // end if !gridsSame  
      
    } // end loop over currently defined levels

  
  // now ensure that any remaining levels are null pointers
  // (in case of de-refinement)
  // unnecessary if gridsSame
  if (!gridsSame)
    {
      for (int lev = a_thicknessVect.size(); lev < m_old_thickness.size(); lev++)
        {
          if (m_old_thickness[lev] != NULL) 
            {
              delete m_old_thickness[lev];
              m_old_thickness[lev] = NULL;
            }
          
          if (m_velocity[lev] != NULL) 
            {
              delete m_velocity[lev];
              m_velocity[lev] = NULL;
            }
          
          if (m_internalEnergy[lev] != NULL) 
            {
              delete m_internalEnergy[lev];
              m_internalEnergy[lev] = NULL;
            }
          
          if (m_tillWaterDepth[lev] != NULL) 
            {
              delete m_internalEnergy[lev];
              m_internalEnergy[lev] = NULL;
            }
          
          if (m_iceFrac[lev] != NULL) 
            {
              delete m_iceFrac[lev];
              m_iceFrac[lev] = NULL;
            }
          
#if BISICLES_Z == BISICLES_LAYERED
          if (m_sInternalEnergy[lev] != NULL) 
            {
              delete m_sInternalEnergy[lev];
              m_sInternalEnergy[lev] = NULL;
            }
          if (m_bInternalEnergy[lev] != NULL) 
            {
              delete m_bInternalEnergy[lev];
              m_bInternalEnergy[lev] = NULL;
            }
#endif	      	      
          
          if (m_velRHS[lev] != NULL)
            {
              delete m_velRHS[lev];
              m_velRHS[lev] = NULL;
            }
          
          if (m_velBasalC[lev] != NULL)
            {
              delete m_velBasalC[lev];
              m_velBasalC[lev] = NULL;
            }
          
	  
          DisjointBoxLayout emptyDBL;
          m_amrGrids[lev] = emptyDBL;
        }
      
      m_finest_level = a_thicknessVect.size()-1;
    
  
  
      // set up counter of number of cells
      for (int lev=0; lev<=m_max_level; lev++)
        {
          m_num_cells[lev] = 0;
          if (lev <= m_finest_level) 
            {
              const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
              LayoutIterator lit = levelGrids.layoutIterator();
              for (lit.begin(); lit.ok(); ++lit)
                {
                  const Box& thisBox = levelGrids.get(lit());
                  m_num_cells[lev] += thisBox.numPts();
                }
            } 
        }
      
      
      // finally, set up covered_level flags
      m_covered_level.resize(m_max_level+1, 0);
      // note that finest level can't be covered.
      for (int lev=m_finest_level-1; lev>=0; lev--)
        {
          
          // if the next finer level is covered, then this one is too.
          if (m_covered_level[lev+1] == 1)
            {
              m_covered_level[lev] = 1;
            }
          else
            {
              // see if the grids finer than this level completely cover it
              IntVectSet fineUncovered(m_amrDomains[lev+1].domainBox());
              const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];
              
              LayoutIterator lit = fineGrids.layoutIterator();
              for (lit.begin(); lit.ok(); ++lit)
                {
                  const Box& thisBox = fineGrids.get(lit());
                  fineUncovered.minus_box(thisBox);
                }
              
              if (fineUncovered.isEmpty()) 
                {
                  m_covered_level[lev] = 1;
                }
            }
        } // end loop over levels to determine covered levels
    }  // end stuff we don't have to do if gridsSame
  
  // tempearture depends on internal energy
  updateTemperature();

  // this is a good time to check for remote ice
  // this shouldn't be necessary
#if 0
  if ((m_eliminate_remote_ice_after_regrid) 
      && !(m_eliminate_remote_ice))
    eliminateRemoteIce();
#endif
  
  //applyCalvingCriterion(CalvingModel::PostRegrid);

  if (m_evolve_velocity)
    {
      //velocity solver needs to be re-defined
      if (!gridsSame)
        {
          defineSolver();
        }
      //solve velocity field, but use the previous initial residual norm in place of this one
      //and force a solve even if other conditions (e.g the timestep interval condition) are not met
      if (a_recalculateVelocity)
        {
          solveVelocityField(true, m_velocitySolveInitialResidualNorm);
        }
    }
  else
    {
      CH_assert(m_evolve_velocity);
      MayDay::Error("AmrIce::regrid() not implemented for !m_evolve_velocity");
    }
  
  // end code lifted from AmrIce::regrid

  
}


/// set ice sheet state based on inputs -- uses the 
/**  New interface function for SUNDIALS and Parallel-in-time
     using an IceSheetState structure, reset ice sheet to match
     TODO: What about cur_step?
     TODO: What about iceFrac?
     TODO: What about thermodynamics?
*/
void
AmrIce::setState(IceSheetState& a_iceState,
                 bool a_recalculateVelocity)
{
  // for now, just use existing function
  setState(a_iceState.ice_thickness,
           a_iceState.ice_velocity,
           a_iceState.time,
           a_recalculateVelocity);  

}


  // companion function to setState -- loads current ice sheet state into structure
  /**  New interface function for SUNDIALS and Parallel-in-time
       loads current state into an IceSheetState structure
       TODO: What about cur_step?
       TODO: What about iceFrac?
       TODO: What about thermodynamics?
  */  
void
AmrIce::getState(IceSheetState& a_iceState)
{

  // thickness
  Vector<LevelData<FArrayBox>* > thickness(m_velocity.size());
  for (int i=0; i<thickness.size(); i++)
    {
      thickness[i] =  &(m_vect_coordSys[i]->getH());
    }
  reshapeAndFill(a_iceState.ice_thickness, thickness);

  // velocity
  reshapeAndFill(a_iceState.ice_velocity, m_velocity);

  // for now, don't worry about enthalpy or area fraction
 
  // current time
  a_iceState.time = m_time;
  

}


#include "NamespaceFooter.H"
