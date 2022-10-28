#ifndef _SIGMACSF_F_H_
#define _SIGMACSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCENTEREDDIFF 
#define GUARDCENTEREDDIFF 
// Prototype for Fortran procedure CENTEREDDIFF ...
//
void FORTRAN_NAME( CENTEREDDIFF ,centereddiff )(
      CHFp_FRA1(deltaPhi)
      ,CHFp_FRA1(phi)
      ,CHFp_BOX(derivBox)
      ,CHFp_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_CENTEREDDIFF FORTRAN_NAME( inlineCENTEREDDIFF, inlineCENTEREDDIFF)
#define FORTNT_CENTEREDDIFF FORTRAN_NAME( CENTEREDDIFF, centereddiff)

inline void FORTRAN_NAME(inlineCENTEREDDIFF, inlineCENTEREDDIFF)(
      CHFp_FRA1(deltaPhi)
      ,CHFp_FRA1(phi)
      ,CHFp_BOX(derivBox)
      ,CHFp_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_CENTEREDDIFF");
 FORTRAN_NAME( CENTEREDDIFF ,centereddiff )(
      CHFt_FRA1(deltaPhi)
      ,CHFt_FRA1(phi)
      ,CHFt_BOX(derivBox)
      ,CHFt_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDCENTEREDDIFF 

#ifndef GUARDDEFINECELLGEOM 
#define GUARDDEFINECELLGEOM 
// Prototype for Fortran procedure DEFINECELLGEOM ...
//
void FORTRAN_NAME( DEFINECELLGEOM ,definecellgeom )(
      CHFp_FRA1(deltaFactors)
      ,CHFp_FRA1(deltaH)
      ,CHFp_FRA1(deltaZb)
      ,CHFp_REALVECT(dx)
      ,CHFp_BOX(box)
      ,CHFp_INT(derivDir) );

#define FORT_DEFINECELLGEOM FORTRAN_NAME( inlineDEFINECELLGEOM, inlineDEFINECELLGEOM)
#define FORTNT_DEFINECELLGEOM FORTRAN_NAME( DEFINECELLGEOM, definecellgeom)

inline void FORTRAN_NAME(inlineDEFINECELLGEOM, inlineDEFINECELLGEOM)(
      CHFp_FRA1(deltaFactors)
      ,CHFp_FRA1(deltaH)
      ,CHFp_FRA1(deltaZb)
      ,CHFp_REALVECT(dx)
      ,CHFp_BOX(box)
      ,CHFp_INT(derivDir) )
{
 CH_TIMELEAF("FORT_DEFINECELLGEOM");
 FORTRAN_NAME( DEFINECELLGEOM ,definecellgeom )(
      CHFt_FRA1(deltaFactors)
      ,CHFt_FRA1(deltaH)
      ,CHFt_FRA1(deltaZb)
      ,CHFt_REALVECT(dx)
      ,CHFt_BOX(box)
      ,CHFt_INT(derivDir) );
}
#endif  // GUARDDEFINECELLGEOM 

#ifndef GUARDFIXFACEH 
#define GUARDFIXFACEH 
// Prototype for Fortran procedure FIXFACEH ...
//
void FORTRAN_NAME( FIXFACEH ,fixfaceh )(
      CHFp_FRA1(xFaceH)
      ,CHFp_FRA1(yFaceH)
      ,CHFp_CONST_FRA1(cellH)
      ,CHFp_BOX(gridbox)
      ,CHFp_INT(xDir)
      ,CHFp_INT(yDir) );

#define FORT_FIXFACEH FORTRAN_NAME( inlineFIXFACEH, inlineFIXFACEH)
#define FORTNT_FIXFACEH FORTRAN_NAME( FIXFACEH, fixfaceh)

inline void FORTRAN_NAME(inlineFIXFACEH, inlineFIXFACEH)(
      CHFp_FRA1(xFaceH)
      ,CHFp_FRA1(yFaceH)
      ,CHFp_CONST_FRA1(cellH)
      ,CHFp_BOX(gridbox)
      ,CHFp_INT(xDir)
      ,CHFp_INT(yDir) )
{
 CH_TIMELEAF("FORT_FIXFACEH");
 FORTRAN_NAME( FIXFACEH ,fixfaceh )(
      CHFt_FRA1(xFaceH)
      ,CHFt_FRA1(yFaceH)
      ,CHFt_CONST_FRA1(cellH)
      ,CHFt_BOX(gridbox)
      ,CHFt_INT(xDir)
      ,CHFt_INT(yDir) );
}
#endif  // GUARDFIXFACEH 

#ifndef GUARDSURFACEHEIGHT 
#define GUARDSURFACEHEIGHT 
// Prototype for Fortran procedure surfaceheight ...
//
void FORTRAN_NAME( SURFACEHEIGHT ,surfaceheight )(
      CHFp_FRA1(zSurface)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(baseHeight)
      ,CHFp_CONST_REAL(iceDensity)
      ,CHFp_CONST_REAL(waterDensity)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_SURFACEHEIGHT FORTRAN_NAME( inlineSURFACEHEIGHT, inlineSURFACEHEIGHT)
#define FORTNT_SURFACEHEIGHT FORTRAN_NAME( SURFACEHEIGHT, surfaceheight)

inline void FORTRAN_NAME(inlineSURFACEHEIGHT, inlineSURFACEHEIGHT)(
      CHFp_FRA1(zSurface)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(baseHeight)
      ,CHFp_CONST_REAL(iceDensity)
      ,CHFp_CONST_REAL(waterDensity)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SURFACEHEIGHT");
 FORTRAN_NAME( SURFACEHEIGHT ,surfaceheight )(
      CHFt_FRA1(zSurface)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA1(baseHeight)
      ,CHFt_CONST_REAL(iceDensity)
      ,CHFt_CONST_REAL(waterDensity)
      ,CHFt_CONST_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSURFACEHEIGHT 

#ifndef GUARDSETFLOATINGMASK 
#define GUARDSETFLOATINGMASK 
// Prototype for Fortran procedure SETFLOATINGMASK ...
//
void FORTRAN_NAME( SETFLOATINGMASK ,setfloatingmask )(
      CHFp_FIA1(floatingMask)
      ,CHFp_CONST_FRA1(Zsurf)
      ,CHFp_CONST_FRA1(zBase)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_INT(anyFloating)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_SETFLOATINGMASK FORTRAN_NAME( inlineSETFLOATINGMASK, inlineSETFLOATINGMASK)
#define FORTNT_SETFLOATINGMASK FORTRAN_NAME( SETFLOATINGMASK, setfloatingmask)

inline void FORTRAN_NAME(inlineSETFLOATINGMASK, inlineSETFLOATINGMASK)(
      CHFp_FIA1(floatingMask)
      ,CHFp_CONST_FRA1(Zsurf)
      ,CHFp_CONST_FRA1(zBase)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_INT(anyFloating)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETFLOATINGMASK");
 FORTRAN_NAME( SETFLOATINGMASK ,setfloatingmask )(
      CHFt_FIA1(floatingMask)
      ,CHFt_CONST_FRA1(Zsurf)
      ,CHFt_CONST_FRA1(zBase)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_INT(anyFloating)
      ,CHFt_REAL(iceDensity)
      ,CHFt_REAL(waterDensity)
      ,CHFt_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETFLOATINGMASK 

#ifndef GUARDPRESERVEMASK 
#define GUARDPRESERVEMASK 
// Prototype for Fortran procedure preservemask ...
//
void FORTRAN_NAME( PRESERVEMASK ,preservemask )(
      CHFp_FRA1(lsrf)
      ,CHFp_FRA1(usrf)
      ,CHFp_FRA1(topg)
      ,CHFp_FRA1(thck)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_PRESERVEMASK FORTRAN_NAME( inlinePRESERVEMASK, inlinePRESERVEMASK)
#define FORTNT_PRESERVEMASK FORTRAN_NAME( PRESERVEMASK, preservemask)

inline void FORTRAN_NAME(inlinePRESERVEMASK, inlinePRESERVEMASK)(
      CHFp_FRA1(lsrf)
      ,CHFp_FRA1(usrf)
      ,CHFp_FRA1(topg)
      ,CHFp_FRA1(thck)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRESERVEMASK");
 FORTRAN_NAME( PRESERVEMASK ,preservemask )(
      CHFt_FRA1(lsrf)
      ,CHFt_FRA1(usrf)
      ,CHFt_FRA1(topg)
      ,CHFt_FRA1(thck)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_REAL(iceDensity)
      ,CHFt_REAL(waterDensity)
      ,CHFt_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRESERVEMASK 

#ifndef GUARDPRESERVEMASKR 
#define GUARDPRESERVEMASKR 
// Prototype for Fortran procedure preservemaskr ...
//
void FORTRAN_NAME( PRESERVEMASKR ,preservemaskr )(
      CHFp_FRA1(lsrf)
      ,CHFp_FRA1(usrf)
      ,CHFp_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_PRESERVEMASKR FORTRAN_NAME( inlinePRESERVEMASKR, inlinePRESERVEMASKR)
#define FORTNT_PRESERVEMASKR FORTRAN_NAME( PRESERVEMASKR, preservemaskr)

inline void FORTRAN_NAME(inlinePRESERVEMASKR, inlinePRESERVEMASKR)(
      CHFp_FRA1(lsrf)
      ,CHFp_FRA1(usrf)
      ,CHFp_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRESERVEMASKR");
 FORTRAN_NAME( PRESERVEMASKR ,preservemaskr )(
      CHFt_FRA1(lsrf)
      ,CHFt_FRA1(usrf)
      ,CHFt_FRA1(topg)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_REAL(iceDensity)
      ,CHFt_REAL(waterDensity)
      ,CHFt_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRESERVEMASKR 

#ifndef GUARDPRESERVEMASKO 
#define GUARDPRESERVEMASKO 
// Prototype for Fortran procedure preservemasko ...
//
void FORTRAN_NAME( PRESERVEMASKO ,preservemasko )(
      CHFp_FRA1(lsrf)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_PRESERVEMASKO FORTRAN_NAME( inlinePRESERVEMASKO, inlinePRESERVEMASKO)
#define FORTNT_PRESERVEMASKO FORTRAN_NAME( PRESERVEMASKO, preservemasko)

inline void FORTRAN_NAME(inlinePRESERVEMASKO, inlinePRESERVEMASKO)(
      CHFp_FRA1(lsrf)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_REAL(iceDensity)
      ,CHFp_REAL(waterDensity)
      ,CHFp_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRESERVEMASKO");
 FORTRAN_NAME( PRESERVEMASKO ,preservemasko )(
      CHFt_FRA1(lsrf)
      ,CHFt_CONST_FRA1(usrf)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_REAL(iceDensity)
      ,CHFt_REAL(waterDensity)
      ,CHFt_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRESERVEMASKO 

#ifndef GUARDODESTROYSCFI 
#define GUARDODESTROYSCFI 
// Prototype for Fortran procedure odestroyscfi ...
//
void FORTRAN_NAME( ODESTROYSCFI ,odestroyscfi )(
      CHFp_FRA1(thck)
      ,CHFp_FIA1(mask)
      ,CHFp_INT(count)
      ,CHFp_CONST_INT(xdir)
      ,CHFp_CONST_INT(ydir)
      ,CHFp_CONST_REAL(thresh)
      ,CHFp_BOX(box) );

#define FORT_ODESTROYSCFI FORTRAN_NAME( inlineODESTROYSCFI, inlineODESTROYSCFI)
#define FORTNT_ODESTROYSCFI FORTRAN_NAME( ODESTROYSCFI, odestroyscfi)

inline void FORTRAN_NAME(inlineODESTROYSCFI, inlineODESTROYSCFI)(
      CHFp_FRA1(thck)
      ,CHFp_FIA1(mask)
      ,CHFp_INT(count)
      ,CHFp_CONST_INT(xdir)
      ,CHFp_CONST_INT(ydir)
      ,CHFp_CONST_REAL(thresh)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ODESTROYSCFI");
 FORTRAN_NAME( ODESTROYSCFI ,odestroyscfi )(
      CHFt_FRA1(thck)
      ,CHFt_FIA1(mask)
      ,CHFt_INT(count)
      ,CHFt_CONST_INT(xdir)
      ,CHFt_CONST_INT(ydir)
      ,CHFt_CONST_REAL(thresh)
      ,CHFt_BOX(box) );
}
#endif  // GUARDODESTROYSCFI 

#ifndef GUARDDESTROYSCFI 
#define GUARDDESTROYSCFI 
// Prototype for Fortran procedure destroyscfi ...
//
void FORTRAN_NAME( DESTROYSCFI ,destroyscfi )(
      CHFp_FRA1(thck)
      ,CHFp_FIA1(mask)
      ,CHFp_INT(count)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(thresh)
      ,CHFp_BOX(box) );

#define FORT_DESTROYSCFI FORTRAN_NAME( inlineDESTROYSCFI, inlineDESTROYSCFI)
#define FORTNT_DESTROYSCFI FORTRAN_NAME( DESTROYSCFI, destroyscfi)

inline void FORTRAN_NAME(inlineDESTROYSCFI, inlineDESTROYSCFI)(
      CHFp_FRA1(thck)
      ,CHFp_FIA1(mask)
      ,CHFp_INT(count)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(thresh)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DESTROYSCFI");
 FORTRAN_NAME( DESTROYSCFI ,destroyscfi )(
      CHFt_FRA1(thck)
      ,CHFt_FIA1(mask)
      ,CHFt_INT(count)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(thresh)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDESTROYSCFI 

#ifndef GUARDSSETOPENSURFACE 
#define GUARDSSETOPENSURFACE 
// Prototype for Fortran procedure ssetopensurface ...
//
void FORTRAN_NAME( SSETOPENSURFACE ,ssetopensurface )(
      CHFp_FRA1(surf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_SSETOPENSURFACE FORTRAN_NAME( inlineSSETOPENSURFACE, inlineSSETOPENSURFACE)
#define FORTNT_SSETOPENSURFACE FORTRAN_NAME( SSETOPENSURFACE, ssetopensurface)

inline void FORTRAN_NAME(inlineSSETOPENSURFACE, inlineSSETOPENSURFACE)(
      CHFp_FRA1(surf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SSETOPENSURFACE");
 FORTRAN_NAME( SSETOPENSURFACE ,ssetopensurface )(
      CHFt_FRA1(surf)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSSETOPENSURFACE 

#ifndef GUARDSGLGRADS 
#define GUARDSGLGRADS 
// Prototype for Fortran procedure sglgrads ...
//
void FORTRAN_NAME( SGLGRADS ,sglgrads )(
      CHFp_FRA1(grads)
      ,CHFp_FRA1(surff)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(surf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_REAL(ratio)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_BOX(faceBox) );

#define FORT_SGLGRADS FORTRAN_NAME( inlineSGLGRADS, inlineSGLGRADS)
#define FORTNT_SGLGRADS FORTRAN_NAME( SGLGRADS, sglgrads)

inline void FORTRAN_NAME(inlineSGLGRADS, inlineSGLGRADS)(
      CHFp_FRA1(grads)
      ,CHFp_FRA1(surff)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(surf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_REAL(ratio)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box)
      ,CHFp_BOX(faceBox) )
{
 CH_TIMELEAF("FORT_SGLGRADS");
 FORTRAN_NAME( SGLGRADS ,sglgrads )(
      CHFt_FRA1(grads)
      ,CHFt_FRA1(surff)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_FRA1(surf)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_REAL(ratio)
      ,CHFt_CONST_REAL(seaLevel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box)
      ,CHFt_BOX(faceBox) );
}
#endif  // GUARDSGLGRADS 

#ifndef GUARDTHICKNESSOVERFLOTATION 
#define GUARDTHICKNESSOVERFLOTATION 
// Prototype for Fortran procedure thicknessoverflotation ...
//
void FORTRAN_NAME( THICKNESSOVERFLOTATION ,thicknessoverflotation )(
      CHFp_FRA1(p)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(sealevel)
      ,CHFp_BOX(box) );

#define FORT_THICKNESSOVERFLOTATION FORTRAN_NAME( inlineTHICKNESSOVERFLOTATION, inlineTHICKNESSOVERFLOTATION)
#define FORTNT_THICKNESSOVERFLOTATION FORTRAN_NAME( THICKNESSOVERFLOTATION, thicknessoverflotation)

inline void FORTRAN_NAME(inlineTHICKNESSOVERFLOTATION, inlineTHICKNESSOVERFLOTATION)(
      CHFp_FRA1(p)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(sealevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_THICKNESSOVERFLOTATION");
 FORTRAN_NAME( THICKNESSOVERFLOTATION ,thicknessoverflotation )(
      CHFt_FRA1(p)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_REAL(rhoi)
      ,CHFt_CONST_REAL(rhoo)
      ,CHFt_CONST_REAL(sealevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDTHICKNESSOVERFLOTATION 

}

#endif
