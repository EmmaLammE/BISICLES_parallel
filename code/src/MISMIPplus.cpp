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
#include "MISMIPplus.H"
#include "AmrIce.H"
#include "NamespaceHeader.H"

MISMIPmelt::MISMIPmelt(ParmParse& a_pp)
  : m_experiment("undefined")
{
  // require experiment to be defined in the inputs fie
  a_pp.get("experiment",m_experiment);
}

// virtual destructor
MISMIPmelt:: ~MISMIPmelt()
{
}

/// factory method: return a pointer to a new SurfaceFlux object
SurfaceFlux*
MISMIPmelt::new_surfaceFlux()
{
  MISMIPmelt* newPtr = new MISMIPmelt();
  newPtr->m_experiment = m_experiment;
  return static_cast<SurfaceFlux*>(newPtr);
  
}


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
void
MISMIPmelt::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
                                 const AmrIceBase& a_amrIce, 
                                 int a_level, Real a_dt)
{

  const RefCountedPtr<LevelSigmaCS> iceGeom = a_amrIce.geometry(a_level);
  const RealVect dx = iceGeom->dx();
  const LevelData<FArrayBox>& thickness = iceGeom->getH();
  const LevelData<FArrayBox>& topography = iceGeom->getTopography();
  const LevelData<FArrayBox>& surfaceHeight = iceGeom->getSurfaceHeight();
  
  DataIterator dit = a_flux.dataIterator();
  
  // experiment 0 is control and spinup (no melting other than calving-front)
  if ((m_experiment == "spin") || (m_experiment == "melt0"))
    {
      DataIterator dit = a_flux.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          a_flux[dit].setVal(0.0);
        }
    }
  else if (m_experiment == "melt4")
    {
      // experiment 4 is melt between 10 <= t <= 110
      Real time = a_amrIce.time();
      Real z0 = -100.0;
      //Real rhoi = a_amrIce.iceDensity();
      //Real rhow = a_amrIce.seawaterDensity();
      Real Omega = -0.2;
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisMelt = a_flux[dit];
          BoxIterator bit(thisMelt.box());
          const FArrayBox& thisThickness = thickness[dit];
          const FArrayBox& thisSurface = surfaceHeight[dit];
          const FArrayBox& thisTopo = topography[dit];
          // not melting above z_0
          // same as melt 3, but turn off after 100 years
          // also, massive melt if x>640km in order to keep calving front fixed.
          if ((time > 10.0) && (time < 110.0))
            {
              // loop over cells in this box (eventually move to fortran?)
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real zb = thisSurface(iv,0) - thisThickness(iv,0);
                  if  (zb < z0)
                    {
                      Real wct = zb - thisTopo(iv,0);
                      thisMelt(iv,0) = -1*Omega*(zb-z0)*tanh(wct/75.0);
                    }                   
                }
            }
          else
            {
              // no melt
              thisMelt.setVal(0.0);
            }
          
          // m = Omega * (z_bottom - z_0) * tanh( (z_bottom - z_bathymetry)/ 75)
        } // end loop over grids

      
    } // end experiment 4
  else if (m_experiment == "melt42")
    {
      // experiment 42 is same as 4, but don't turn off 
      Real time = a_amrIce.time();
      Real z0 = -100.0;
      //Real rhoi = a_amrIce.iceDensity();
      //Real rhow = a_amrIce.seawaterDensity();
      Real Omega = -0.2;
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisMelt = a_flux[dit];
          BoxIterator bit(thisMelt.box());
          const FArrayBox& thisThickness = thickness[dit];
          const FArrayBox& thisSurface = surfaceHeight[dit];
          const FArrayBox& thisTopo = topography[dit];
          // not melting above z_0
          // same as melt 3, but turn off after 100 years
          // also, massive melt if x>640km in order to keep calving front fixed.
          if (time > 10.0)
            {
              // loop over cells in this box (eventually move to fortran?)
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real zb = thisSurface(iv,0) - thisThickness(iv,0);
                  if  (zb < z0)
                    {
                      Real wct = zb - thisTopo(iv,0);
                      thisMelt(iv,0) = -1*Omega*(zb-z0)*tanh(wct/75.0);
                    }                   
                }
            }
          else
            {
              // no melt
              thisMelt.setVal(0.0);
            }
          
          // m = Omega * (z_bottom - z_0) * tanh( (z_bottom - z_bathymetry)/ 75)
        } // end loop over grids

      
    } // end experiment 42    
           
  else
    {
      MayDay::Error("MISMIPmelt::surfaceThicknessFlux == unknown experiment");
    }

  // add in extreme melt to preserve calving front
  bool applyCFsource = true;

  if (applyCFsource)
    {
      DataIterator dit = a_flux.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisMelt = a_flux[dit];
          BoxIterator bit(thisMelt.box());
          
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();              
              Real x = (Real(iv[0]) + 0.5)*dx[0];
              if (x > 640000)
                {
                  thisMelt(iv,0) -= 1.0e4;
                }
            } // end loop over cells
        } // end loop over boxes for calving-front melt      
    } // end if we're applying cf sourcex
}




#include "NamespaceFooter.H"
