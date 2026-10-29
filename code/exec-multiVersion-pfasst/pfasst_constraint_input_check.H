#ifndef PFASST_CONSTRAINT_INPUT_CHECK_H
#define PFASST_CONSTRAINT_INPUT_CHECK_H

#include "ParmParse.H"
#include "MayDay.H"
#include <string>

// Call once, before MPI communicator splitting and before AmrIce::initialize().
// This is a deliberately restricted thickness-only PFASST configuration check.
// It does not change any input option or any serial model silently.
static void checkPFASSTThicknessOnlyInput()
{
  const char* amrPrefixes[] = {"crse.amr", "fine.amr"};
  for (int k = 0; k < 2; ++k)
    {
      ParmParse pp(amrPrefixes[k]);
      int maxLevel = -1;
      pp.query("maxLevel", maxLevel);
      if (maxLevel < 0)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".maxLevel must be a nonnegative maximum AMR level.";
          MayDay::Error(message.c_str());
        }
      bool evolveIceFrac = true;
      pp.query("evolve_ice_frac", evolveIceFrac);
      if (evolveIceFrac)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".evolve_ice_frac must be explicitly false for thickness-only PFASST. "
            "PFASST does not evolve or communicate ice fraction.";
          MayDay::Error(message.c_str());
        }
      bool maskSources = false;
      pp.query("mask_sources", maskSources);
      if (maskSources)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".mask_sources must be false for this thickness-only PFASST patch.";
          MayDay::Error(message.c_str());
        }
      bool evolveThickness = true;
      pp.query("evolve_thickness", evolveThickness);
      if (!evolveThickness)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".evolve_thickness must be true for this PFASST patch.";
          MayDay::Error(message.c_str());
        }
      bool isothermal = true;
      pp.query("isothermal", isothermal);
      if (!isothermal)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".isothermal must be true: PFASST currently evolves thickness only.";
          MayDay::Error(message.c_str());
        }
      bool changeBed = false;
      pp.query("evolve_topography_fix_surface", changeBed);
      if (changeBed)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".evolve_topography_fix_surface must be false: bed elevation is not a PFASST state.";
          MayDay::Error(message.c_str());
        }
      bool remoteIce = false;
      pp.query("eliminate_remote_ice", remoteIce);
      if (remoteIce)
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".eliminate_remote_ice must be false for the supported DomainEdge projection.";
          MayDay::Error(message.c_str());
        }
      std::string diffusion = "none";
      pp.query("diffusion_treatment", diffusion);
      if (diffusion != "none")
        {
          const std::string message = std::string(amrPrefixes[k]) +
            ".diffusion_treatment must be none for this explicit thickness-only PFASST patch.";
          MayDay::Error(message.c_str());
        }
    }

  ParmParse calving("CalvingModel");
  std::string calvingType = "NoCalvingModel";
  calving.query("type", calvingType);
  if (calvingType != "DomainEdgeCalvingModel" && calvingType != "NoCalvingModel")
    MayDay::Error("This PFASST patch supports only DomainEdgeCalvingModel or NoCalvingModel.");
  if (calvingType == "DomainEdgeCalvingModel")
    {
      bool preserveSea = false;
      bool preserveLand = false;
      calving.query("preserveSea", preserveSea);
      calving.query("preserveLand", preserveLand);
      if (preserveSea || preserveLand)
        MayDay::Error("PFASST DomainEdge projection requires preserveSea=false and preserveLand=false.");
    }

  ParmParse topography("topographyFlux");
  std::string topographyType = "zeroFlux";
  topography.query("type", topographyType);
  if (topographyType != "zeroFlux")
    MayDay::Error("PFASST thickness-only integration requires absent or zeroFlux topographyFlux.");

  ParmParse mainOptions("main");
  bool damage = false;
  bool melange = false;
  mainOptions.query("damage_model", damage);
  mainOptions.query("melange_model", melange);
  if (damage || melange)
    MayDay::Error("This PFASST thickness-only patch does not evolve damage or melange observer states.");
}

#endif