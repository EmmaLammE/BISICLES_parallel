#ifndef _ICETHERMODYNAMICSF_F_H_
#define _ICETHERMODYNAMICSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPOSEINTERNALENERGYICE 
#define GUARDCOMPOSEINTERNALENERGYICE 
// Prototype for Fortran procedure composeinternalenergyice ...
//
void FORTRAN_NAME( COMPOSEINTERNALENERGYICE ,composeinternalenergyice )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(T)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(C)
      ,CHFp_CONST_REAL(L)
      ,CHFp_BOX(box) );

#define FORT_COMPOSEINTERNALENERGYICE FORTRAN_NAME( inlineCOMPOSEINTERNALENERGYICE, inlineCOMPOSEINTERNALENERGYICE)
#define FORTNT_COMPOSEINTERNALENERGYICE FORTRAN_NAME( COMPOSEINTERNALENERGYICE, composeinternalenergyice)

inline void FORTRAN_NAME(inlineCOMPOSEINTERNALENERGYICE, inlineCOMPOSEINTERNALENERGYICE)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(T)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(C)
      ,CHFp_CONST_REAL(L)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COMPOSEINTERNALENERGYICE");
 FORTRAN_NAME( COMPOSEINTERNALENERGYICE ,composeinternalenergyice )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(T)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(C)
      ,CHFt_CONST_REAL(L)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOMPOSEINTERNALENERGYICE 

#ifndef GUARDDECOMPOSEINTERNALENERGYICE 
#define GUARDDECOMPOSEINTERNALENERGYICE 
// Prototype for Fortran procedure decomposeinternalenergyice ...
//
void FORTRAN_NAME( DECOMPOSEINTERNALENERGYICE ,decomposeinternalenergyice )(
      CHFp_FRA(T)
      ,CHFp_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(P)
      ,CHFp_CONST_REAL(C)
      ,CHFp_CONST_REAL(L)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(Tr)
      ,CHFp_CONST_REAL(Tmin)
      ,CHFp_BOX(box) );

#define FORT_DECOMPOSEINTERNALENERGYICE FORTRAN_NAME( inlineDECOMPOSEINTERNALENERGYICE, inlineDECOMPOSEINTERNALENERGYICE)
#define FORTNT_DECOMPOSEINTERNALENERGYICE FORTRAN_NAME( DECOMPOSEINTERNALENERGYICE, decomposeinternalenergyice)

inline void FORTRAN_NAME(inlineDECOMPOSEINTERNALENERGYICE, inlineDECOMPOSEINTERNALENERGYICE)(
      CHFp_FRA(T)
      ,CHFp_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(P)
      ,CHFp_CONST_REAL(C)
      ,CHFp_CONST_REAL(L)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(Tr)
      ,CHFp_CONST_REAL(Tmin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DECOMPOSEINTERNALENERGYICE");
 FORTRAN_NAME( DECOMPOSEINTERNALENERGYICE ,decomposeinternalenergyice )(
      CHFt_FRA(T)
      ,CHFt_FRA(W)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA(P)
      ,CHFt_CONST_REAL(C)
      ,CHFt_CONST_REAL(L)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(Tr)
      ,CHFt_CONST_REAL(Tmin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDECOMPOSEINTERNALENERGYICE 

#ifndef GUARDCOLUMNTHERMODYAMICSSETCONSTANTS 
#define GUARDCOLUMNTHERMODYAMICSSETCONSTANTS 
// Prototype for Fortran procedure COLUMNTHERMODYAMICSSETCONSTANTS ...
//
void FORTRAN_NAME( COLUMNTHERMODYAMICSSETCONSTANTS ,columnthermodyamicssetconstants )(
      CHFp_CONST_REAL(a_seconds_per_unit_time)
      ,CHFp_CONST_REAL(a_ice_density)
      ,CHFp_CONST_REAL(a_water_density)
      ,CHFp_CONST_REAL(a_gravity)
      ,CHFp_CONST_REAL(a_ice_heat_capacity)
      ,CHFp_CONST_REAL(a_ice_latent_heat)
      ,CHFp_CONST_REAL(a_ice_conductivity)
      ,CHFp_CONST_REAL(a_moisture_conductivity)
      ,CHFp_CONST_REAL(a_ice_pressure_melt_factor)
      ,CHFp_CONST_REAL(a_triple_point)
      ,CHFp_CONST_REAL(a_water_fraction_drain)
      ,CHFp_CONST_REAL(a_water_fraction_max)
      ,CHFp_CONST_REAL(a_water_drain_factor)
      ,CHFp_CONST_REAL(a_till_water_drain_factor)
      ,CHFp_CONST_REAL(a_till_water_max) );

#define FORT_COLUMNTHERMODYAMICSSETCONSTANTS FORTRAN_NAME( inlineCOLUMNTHERMODYAMICSSETCONSTANTS, inlineCOLUMNTHERMODYAMICSSETCONSTANTS)
#define FORTNT_COLUMNTHERMODYAMICSSETCONSTANTS FORTRAN_NAME( COLUMNTHERMODYAMICSSETCONSTANTS, columnthermodyamicssetconstants)

inline void FORTRAN_NAME(inlineCOLUMNTHERMODYAMICSSETCONSTANTS, inlineCOLUMNTHERMODYAMICSSETCONSTANTS)(
      CHFp_CONST_REAL(a_seconds_per_unit_time)
      ,CHFp_CONST_REAL(a_ice_density)
      ,CHFp_CONST_REAL(a_water_density)
      ,CHFp_CONST_REAL(a_gravity)
      ,CHFp_CONST_REAL(a_ice_heat_capacity)
      ,CHFp_CONST_REAL(a_ice_latent_heat)
      ,CHFp_CONST_REAL(a_ice_conductivity)
      ,CHFp_CONST_REAL(a_moisture_conductivity)
      ,CHFp_CONST_REAL(a_ice_pressure_melt_factor)
      ,CHFp_CONST_REAL(a_triple_point)
      ,CHFp_CONST_REAL(a_water_fraction_drain)
      ,CHFp_CONST_REAL(a_water_fraction_max)
      ,CHFp_CONST_REAL(a_water_drain_factor)
      ,CHFp_CONST_REAL(a_till_water_drain_factor)
      ,CHFp_CONST_REAL(a_till_water_max) )
{
 CH_TIMELEAF("FORT_COLUMNTHERMODYAMICSSETCONSTANTS");
 FORTRAN_NAME( COLUMNTHERMODYAMICSSETCONSTANTS ,columnthermodyamicssetconstants )(
      CHFt_CONST_REAL(a_seconds_per_unit_time)
      ,CHFt_CONST_REAL(a_ice_density)
      ,CHFt_CONST_REAL(a_water_density)
      ,CHFt_CONST_REAL(a_gravity)
      ,CHFt_CONST_REAL(a_ice_heat_capacity)
      ,CHFt_CONST_REAL(a_ice_latent_heat)
      ,CHFt_CONST_REAL(a_ice_conductivity)
      ,CHFt_CONST_REAL(a_moisture_conductivity)
      ,CHFt_CONST_REAL(a_ice_pressure_melt_factor)
      ,CHFt_CONST_REAL(a_triple_point)
      ,CHFt_CONST_REAL(a_water_fraction_drain)
      ,CHFt_CONST_REAL(a_water_fraction_max)
      ,CHFt_CONST_REAL(a_water_drain_factor)
      ,CHFt_CONST_REAL(a_till_water_drain_factor)
      ,CHFt_CONST_REAL(a_till_water_max) );
}
#endif  // GUARDCOLUMNTHERMODYAMICSSETCONSTANTS 

#ifndef GUARDUPDATEINTERNALENERGY 
#define GUARDUPDATEINTERNALENERGY 
// Prototype for Fortran procedure UPDATEINTERNALENERGY ...
//
void FORTRAN_NAME( UPDATEINTERNALENERGY ,updateinternalenergy )(
      CHFp_FRA(energy)
      ,CHFp_FRA1(tillwaterdepth)
      ,CHFp_FRA1(senergy)
      ,CHFp_FRA1(benergy)
      ,CHFp_CONST_FRA1(sflux)
      ,CHFp_CONST_FRA1(bflux)
      ,CHFp_CONST_FRA1(tillWaterDrainFactor)
      ,CHFp_CONST_FIA1(floatingMaskOld)
      ,CHFp_CONST_FIA1(floatingMaskNew)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(thckold)
      ,CHFp_CONST_FRA1(thcknew)
      ,CHFp_CONST_FRA(usig)
      ,CHFp_CONST_VR(sigma)
      ,CHFp_CONST_VR(dsigma)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(nlay)
      ,CHFp_CONST_INT(sdiric)
      ,CHFp_BOX(box) );

#define FORT_UPDATEINTERNALENERGY FORTRAN_NAME( inlineUPDATEINTERNALENERGY, inlineUPDATEINTERNALENERGY)
#define FORTNT_UPDATEINTERNALENERGY FORTRAN_NAME( UPDATEINTERNALENERGY, updateinternalenergy)

inline void FORTRAN_NAME(inlineUPDATEINTERNALENERGY, inlineUPDATEINTERNALENERGY)(
      CHFp_FRA(energy)
      ,CHFp_FRA1(tillwaterdepth)
      ,CHFp_FRA1(senergy)
      ,CHFp_FRA1(benergy)
      ,CHFp_CONST_FRA1(sflux)
      ,CHFp_CONST_FRA1(bflux)
      ,CHFp_CONST_FRA1(tillWaterDrainFactor)
      ,CHFp_CONST_FIA1(floatingMaskOld)
      ,CHFp_CONST_FIA1(floatingMaskNew)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(thckold)
      ,CHFp_CONST_FRA1(thcknew)
      ,CHFp_CONST_FRA(usig)
      ,CHFp_CONST_VR(sigma)
      ,CHFp_CONST_VR(dsigma)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(nlay)
      ,CHFp_CONST_INT(sdiric)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATEINTERNALENERGY");
 FORTRAN_NAME( UPDATEINTERNALENERGY ,updateinternalenergy )(
      CHFt_FRA(energy)
      ,CHFt_FRA1(tillwaterdepth)
      ,CHFt_FRA1(senergy)
      ,CHFt_FRA1(benergy)
      ,CHFt_CONST_FRA1(sflux)
      ,CHFt_CONST_FRA1(bflux)
      ,CHFt_CONST_FRA1(tillWaterDrainFactor)
      ,CHFt_CONST_FIA1(floatingMaskOld)
      ,CHFt_CONST_FIA1(floatingMaskNew)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(thckold)
      ,CHFt_CONST_FRA1(thcknew)
      ,CHFt_CONST_FRA(usig)
      ,CHFt_CONST_VR(sigma)
      ,CHFt_CONST_VR(dsigma)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(nlay)
      ,CHFt_CONST_INT(sdiric)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATEINTERNALENERGY 

}

#endif
