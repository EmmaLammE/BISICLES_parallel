#ifndef _AMRICEF_F_H_
#define _AMRICEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUNDIVIDEDGRAD 
#define GUARDUNDIVIDEDGRAD 
// Prototype for Fortran procedure UNDIVIDEDGRAD ...
//
void FORTRAN_NAME( UNDIVIDEDGRAD ,undividedgrad )(
      CHFp_FRA1(dq)
      ,CHFp_CONST_FRA(q)
      ,CHFp_BOX(dInterior)
      ,CHFp_BOX(loedge)
      ,CHFp_BOX(hiedge)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi) );

#define FORT_UNDIVIDEDGRAD FORTRAN_NAME( inlineUNDIVIDEDGRAD, inlineUNDIVIDEDGRAD)
#define FORTNT_UNDIVIDEDGRAD FORTRAN_NAME( UNDIVIDEDGRAD, undividedgrad)

inline void FORTRAN_NAME(inlineUNDIVIDEDGRAD, inlineUNDIVIDEDGRAD)(
      CHFp_FRA1(dq)
      ,CHFp_CONST_FRA(q)
      ,CHFp_BOX(dInterior)
      ,CHFp_BOX(loedge)
      ,CHFp_BOX(hiedge)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi) )
{
 CH_TIMELEAF("FORT_UNDIVIDEDGRAD");
 FORTRAN_NAME( UNDIVIDEDGRAD ,undividedgrad )(
      CHFt_FRA1(dq)
      ,CHFt_CONST_FRA(q)
      ,CHFt_BOX(dInterior)
      ,CHFt_BOX(loedge)
      ,CHFt_BOX(hiedge)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi) );
}
#endif  // GUARDUNDIVIDEDGRAD 

#ifndef GUARDSETONMASK 
#define GUARDSETONMASK 
// Prototype for Fortran procedure SETONMASK ...
//
void FORTRAN_NAME( SETONMASK ,setonmask )(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_INT(maskval)
      ,CHFp_CONST_REAL(fabval)
      ,CHFp_BOX(box) );

#define FORT_SETONMASK FORTRAN_NAME( inlineSETONMASK, inlineSETONMASK)
#define FORTNT_SETONMASK FORTRAN_NAME( SETONMASK, setonmask)

inline void FORTRAN_NAME(inlineSETONMASK, inlineSETONMASK)(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_INT(maskval)
      ,CHFp_CONST_REAL(fabval)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETONMASK");
 FORTRAN_NAME( SETONMASK ,setonmask )(
      CHFt_FRA1(fab)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_INT(maskval)
      ,CHFt_CONST_REAL(fabval)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETONMASK 

#ifndef GUARDSWEEPCONNECTED2D 
#define GUARDSWEEPCONNECTED2D 
// Prototype for Fortran procedure sweepconnected2d ...
//
void FORTRAN_NAME( SWEEPCONNECTED2D ,sweepconnected2d )(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(conn)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_BOX(box) );

#define FORT_SWEEPCONNECTED2D FORTRAN_NAME( inlineSWEEPCONNECTED2D, inlineSWEEPCONNECTED2D)
#define FORTNT_SWEEPCONNECTED2D FORTRAN_NAME( SWEEPCONNECTED2D, sweepconnected2d)

inline void FORTRAN_NAME(inlineSWEEPCONNECTED2D, inlineSWEEPCONNECTED2D)(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(conn)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWEEPCONNECTED2D");
 FORTRAN_NAME( SWEEPCONNECTED2D ,sweepconnected2d )(
      CHFt_FRA1(fab)
      ,CHFt_CONST_FRA1(conn)
      ,CHFt_CONST_REAL(tol)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWEEPCONNECTED2D 

#ifndef GUARDSETFLOATINGBETA 
#define GUARDSETFLOATINGBETA 
// Prototype for Fortran procedure SETFLOATINGBETA ...
//
void FORTRAN_NAME( SETFLOATINGBETA ,setfloatingbeta )(
      CHFp_FRA1(beta)
      ,CHFp_CONST_FIA1(floatingMask)
      ,CHFp_BOX(gridBox) );

#define FORT_SETFLOATINGBETA FORTRAN_NAME( inlineSETFLOATINGBETA, inlineSETFLOATINGBETA)
#define FORTNT_SETFLOATINGBETA FORTRAN_NAME( SETFLOATINGBETA, setfloatingbeta)

inline void FORTRAN_NAME(inlineSETFLOATINGBETA, inlineSETFLOATINGBETA)(
      CHFp_FRA1(beta)
      ,CHFp_CONST_FIA1(floatingMask)
      ,CHFp_BOX(gridBox) )
{
 CH_TIMELEAF("FORT_SETFLOATINGBETA");
 FORTRAN_NAME( SETFLOATINGBETA ,setfloatingbeta )(
      CHFt_FRA1(beta)
      ,CHFt_CONST_FIA1(floatingMask)
      ,CHFt_BOX(gridBox) );
}
#endif  // GUARDSETFLOATINGBETA 

#ifndef GUARDCORNFORDCORRECTION 
#define GUARDCORNFORDCORRECTION 
// Prototype for Fortran procedure CORNFORDCORRECTION ...
//
void FORTRAN_NAME( CORNFORDCORRECTION ,cornfordcorrection )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(Zsurf)
      ,CHFp_CONST_FIA1(floatingMask)
      ,CHFp_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(gridBox) );

#define FORT_CORNFORDCORRECTION FORTRAN_NAME( inlineCORNFORDCORRECTION, inlineCORNFORDCORRECTION)
#define FORTNT_CORNFORDCORRECTION FORTRAN_NAME( CORNFORDCORRECTION, cornfordcorrection)

inline void FORTRAN_NAME(inlineCORNFORDCORRECTION, inlineCORNFORDCORRECTION)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(Zsurf)
      ,CHFp_CONST_FIA1(floatingMask)
      ,CHFp_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(gridBox) )
{
 CH_TIMELEAF("FORT_CORNFORDCORRECTION");
 FORTRAN_NAME( CORNFORDCORRECTION ,cornfordcorrection )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(Zsurf)
      ,CHFt_CONST_FIA1(floatingMask)
      ,CHFt_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(gridBox) );
}
#endif  // GUARDCORNFORDCORRECTION 

#ifndef GUARDSETOPENSURFACE 
#define GUARDSETOPENSURFACE 
// Prototype for Fortran procedure SETOPENSURFACE ...
//
void FORTRAN_NAME( SETOPENSURFACE ,setopensurface )(
      CHFp_FRA1(surf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_BOX(box) );

#define FORT_SETOPENSURFACE FORTRAN_NAME( inlineSETOPENSURFACE, inlineSETOPENSURFACE)
#define FORTNT_SETOPENSURFACE FORTRAN_NAME( SETOPENSURFACE, setopensurface)

inline void FORTRAN_NAME(inlineSETOPENSURFACE, inlineSETOPENSURFACE)(
      CHFp_FRA1(surf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(seaLevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETOPENSURFACE");
 FORTRAN_NAME( SETOPENSURFACE ,setopensurface )(
      CHFt_FRA1(surf)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_REAL(seaLevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETOPENSURFACE 

#ifndef GUARDSETICEFREEVAL 
#define GUARDSETICEFREEVAL 
// Prototype for Fortran procedure seticefreeval ...
//
void FORTRAN_NAME( SETICEFREEVAL ,seticefreeval )(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_REAL(val)
      ,CHFp_BOX(cellBox) );

#define FORT_SETICEFREEVAL FORTRAN_NAME( inlineSETICEFREEVAL, inlineSETICEFREEVAL)
#define FORTNT_SETICEFREEVAL FORTRAN_NAME( SETICEFREEVAL, seticefreeval)

inline void FORTRAN_NAME(inlineSETICEFREEVAL, inlineSETICEFREEVAL)(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_REAL(val)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_SETICEFREEVAL");
 FORTRAN_NAME( SETICEFREEVAL ,seticefreeval )(
      CHFt_FRA1(fab)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_REAL(val)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDSETICEFREEVAL 

#ifndef GUARDSETFRONTFACEVT 
#define GUARDSETFRONTFACEVT 
// Prototype for Fortran procedure setfrontfacevt ...
//
void FORTRAN_NAME( SETFRONTFACEVT ,setfrontfacevt )(
      CHFp_FRA1(facevt)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_INT(k)
      ,CHFp_CONST_REAL(factor)
      ,CHFp_BOX(cellBox) );

#define FORT_SETFRONTFACEVT FORTRAN_NAME( inlineSETFRONTFACEVT, inlineSETFRONTFACEVT)
#define FORTNT_SETFRONTFACEVT FORTRAN_NAME( SETFRONTFACEVT, setfrontfacevt)

inline void FORTRAN_NAME(inlineSETFRONTFACEVT, inlineSETFRONTFACEVT)(
      CHFp_FRA1(facevt)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_INT(k)
      ,CHFp_CONST_REAL(factor)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_SETFRONTFACEVT");
 FORTRAN_NAME( SETFRONTFACEVT ,setfrontfacevt )(
      CHFt_FRA1(facevt)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_FRA1(usrf)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_INT(k)
      ,CHFt_CONST_REAL(factor)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDSETFRONTFACEVT 

#ifndef GUARDEXTRAPTOMARGIN 
#define GUARDEXTRAPTOMARGIN 
// Prototype for Fortran procedure extraptomargin ...
//
void FORTRAN_NAME( EXTRAPTOMARGIN ,extraptomargin )(
      CHFp_FRA1(uf)
      ,CHFp_FRA1(vface)
      ,CHFp_CONST_FRA1(ufin)
      ,CHFp_CONST_FRA1(uc)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FRA1(thk)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(facebox) );

#define FORT_EXTRAPTOMARGIN FORTRAN_NAME( inlineEXTRAPTOMARGIN, inlineEXTRAPTOMARGIN)
#define FORTNT_EXTRAPTOMARGIN FORTRAN_NAME( EXTRAPTOMARGIN, extraptomargin)

inline void FORTRAN_NAME(inlineEXTRAPTOMARGIN, inlineEXTRAPTOMARGIN)(
      CHFp_FRA1(uf)
      ,CHFp_FRA1(vface)
      ,CHFp_CONST_FRA1(ufin)
      ,CHFp_CONST_FRA1(uc)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FRA1(thk)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_EXTRAPTOMARGIN");
 FORTRAN_NAME( EXTRAPTOMARGIN ,extraptomargin )(
      CHFt_FRA1(uf)
      ,CHFt_FRA1(vface)
      ,CHFt_CONST_FRA1(ufin)
      ,CHFt_CONST_FRA1(uc)
      ,CHFt_CONST_FRA1(usrf)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_FRA1(thk)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDEXTRAPTOMARGIN 

#ifndef GUARDGLCORRECTION 
#define GUARDGLCORRECTION 
// Prototype for Fortran procedure GLCORRECTION ...
//
void FORTRAN_NAME( GLCORRECTION ,glcorrection )(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(Zsurf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(rhog)
      ,CHFp_BOX(gridBox) );

#define FORT_GLCORRECTION FORTRAN_NAME( inlineGLCORRECTION, inlineGLCORRECTION)
#define FORTNT_GLCORRECTION FORTRAN_NAME( GLCORRECTION, glcorrection)

inline void FORTRAN_NAME(inlineGLCORRECTION, inlineGLCORRECTION)(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(Zsurf)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(rhog)
      ,CHFp_BOX(gridBox) )
{
 CH_TIMELEAF("FORT_GLCORRECTION");
 FORTRAN_NAME( GLCORRECTION ,glcorrection )(
      CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA1(Zsurf)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(rhog)
      ,CHFt_BOX(gridBox) );
}
#endif  // GUARDGLCORRECTION 

#ifndef GUARDSETTHICKDIFF 
#define GUARDSETTHICKDIFF 
// Prototype for Fortran procedure SETTHICKDIFF ...
//
void FORTRAN_NAME( SETTHICKDIFF ,setthickdiff )(
      CHFp_FRA1(D)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_REAL(rg)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(faceBox) );

#define FORT_SETTHICKDIFF FORTRAN_NAME( inlineSETTHICKDIFF, inlineSETTHICKDIFF)
#define FORTNT_SETTHICKDIFF FORTRAN_NAME( SETTHICKDIFF, setthickdiff)

inline void FORTRAN_NAME(inlineSETTHICKDIFF, inlineSETTHICKDIFF)(
      CHFp_FRA1(D)
      ,CHFp_CONST_FRA1(C)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_REAL(rg)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(faceBox) )
{
 CH_TIMELEAF("FORT_SETTHICKDIFF");
 FORTRAN_NAME( SETTHICKDIFF ,setthickdiff )(
      CHFt_FRA1(D)
      ,CHFt_CONST_FRA1(C)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_REAL(rg)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(faceBox) );
}
#endif  // GUARDSETTHICKDIFF 

#ifndef GUARDSUBTRACTDVEL 
#define GUARDSUBTRACTDVEL 
// Prototype for Fortran procedure SUBTRACTDVEL ...
//
void FORTRAN_NAME( SUBTRACTDVEL ,subtractdvel )(
      CHFp_FRA1(faceVel)
      ,CHFp_CONST_FRA1(cellH)
      ,CHFp_CONST_FRA1(faceH)
      ,CHFp_CONST_FRA1(faceD)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(faceBox) );

#define FORT_SUBTRACTDVEL FORTRAN_NAME( inlineSUBTRACTDVEL, inlineSUBTRACTDVEL)
#define FORTNT_SUBTRACTDVEL FORTRAN_NAME( SUBTRACTDVEL, subtractdvel)

inline void FORTRAN_NAME(inlineSUBTRACTDVEL, inlineSUBTRACTDVEL)(
      CHFp_FRA1(faceVel)
      ,CHFp_CONST_FRA1(cellH)
      ,CHFp_CONST_FRA1(faceH)
      ,CHFp_CONST_FRA1(faceD)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(faceBox) )
{
 CH_TIMELEAF("FORT_SUBTRACTDVEL");
 FORTRAN_NAME( SUBTRACTDVEL ,subtractdvel )(
      CHFt_FRA1(faceVel)
      ,CHFt_CONST_FRA1(cellH)
      ,CHFt_CONST_FRA1(faceH)
      ,CHFt_CONST_FRA1(faceD)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(faceBox) );
}
#endif  // GUARDSUBTRACTDVEL 

#ifndef GUARDSUBTRACTDFLUX 
#define GUARDSUBTRACTDFLUX 
// Prototype for Fortran procedure SUBTRACTDFLUX ...
//
void FORTRAN_NAME( SUBTRACTDFLUX ,subtractdflux )(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(cellH)
      ,CHFp_CONST_FRA1(faceD)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(faceBox) );

#define FORT_SUBTRACTDFLUX FORTRAN_NAME( inlineSUBTRACTDFLUX, inlineSUBTRACTDFLUX)
#define FORTNT_SUBTRACTDFLUX FORTRAN_NAME( SUBTRACTDFLUX, subtractdflux)

inline void FORTRAN_NAME(inlineSUBTRACTDFLUX, inlineSUBTRACTDFLUX)(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(cellH)
      ,CHFp_CONST_FRA1(faceD)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(faceBox) )
{
 CH_TIMELEAF("FORT_SUBTRACTDFLUX");
 FORTRAN_NAME( SUBTRACTDFLUX ,subtractdflux )(
      CHFt_FRA1(faceFlux)
      ,CHFt_CONST_FRA1(cellH)
      ,CHFt_CONST_FRA1(faceD)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(faceBox) );
}
#endif  // GUARDSUBTRACTDFLUX 

#ifndef GUARDVELINITIALGUESS 
#define GUARDVELINITIALGUESS 
// Prototype for Fortran procedure VELINITIALGUESS ...
//
void FORTRAN_NAME( VELINITIALGUESS ,velinitialguess )(
      CHFp_FRA(vel)
      ,CHFp_FRA(rhs)
      ,CHFp_FRA1(beta)
      ,CHFp_BOX(box) );

#define FORT_VELINITIALGUESS FORTRAN_NAME( inlineVELINITIALGUESS, inlineVELINITIALGUESS)
#define FORTNT_VELINITIALGUESS FORTRAN_NAME( VELINITIALGUESS, velinitialguess)

inline void FORTRAN_NAME(inlineVELINITIALGUESS, inlineVELINITIALGUESS)(
      CHFp_FRA(vel)
      ,CHFp_FRA(rhs)
      ,CHFp_FRA1(beta)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VELINITIALGUESS");
 FORTRAN_NAME( VELINITIALGUESS ,velinitialguess )(
      CHFt_FRA(vel)
      ,CHFt_FRA(rhs)
      ,CHFt_FRA1(beta)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVELINITIALGUESS 

#ifndef GUARDUPWINDLAYERFLUX 
#define GUARDUPWINDLAYERFLUX 
// Prototype for Fortran procedure UPWINDLAYERFLUX ...
//
void FORTRAN_NAME( UPWINDLAYERFLUX ,upwindlayerflux )(
      CHFp_FRA1(f)
      ,CHFp_CONST_FRA1(phiD)
      ,CHFp_CONST_FRA1(phiU)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_BOX(box) );

#define FORT_UPWINDLAYERFLUX FORTRAN_NAME( inlineUPWINDLAYERFLUX, inlineUPWINDLAYERFLUX)
#define FORTNT_UPWINDLAYERFLUX FORTRAN_NAME( UPWINDLAYERFLUX, upwindlayerflux)

inline void FORTRAN_NAME(inlineUPWINDLAYERFLUX, inlineUPWINDLAYERFLUX)(
      CHFp_FRA1(f)
      ,CHFp_CONST_FRA1(phiD)
      ,CHFp_CONST_FRA1(phiU)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPWINDLAYERFLUX");
 FORTRAN_NAME( UPWINDLAYERFLUX ,upwindlayerflux )(
      CHFt_FRA1(f)
      ,CHFt_CONST_FRA1(phiD)
      ,CHFt_CONST_FRA1(phiU)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPWINDLAYERFLUX 

#ifndef GUARDFABMINPLUS 
#define GUARDFABMINPLUS 
// Prototype for Fortran procedure FABMINPLUS ...
//
void FORTRAN_NAME( FABMINPLUS ,fabminplus )(
      CHFp_FRA1(fs)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(p)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_REAL(fmax)
      ,CHFp_BOX(box) );

#define FORT_FABMINPLUS FORTRAN_NAME( inlineFABMINPLUS, inlineFABMINPLUS)
#define FORTNT_FABMINPLUS FORTRAN_NAME( FABMINPLUS, fabminplus)

inline void FORTRAN_NAME(inlineFABMINPLUS, inlineFABMINPLUS)(
      CHFp_FRA1(fs)
      ,CHFp_CONST_FRA1(f)
      ,CHFp_CONST_FRA1(p)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_REAL(fmax)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FABMINPLUS");
 FORTRAN_NAME( FABMINPLUS ,fabminplus )(
      CHFt_FRA1(fs)
      ,CHFt_CONST_FRA1(f)
      ,CHFt_CONST_FRA1(p)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_REAL(fmax)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFABMINPLUS 

#ifndef GUARDADVECTFRAC_SLC_SUBOPTIMAL 
#define GUARDADVECTFRAC_SLC_SUBOPTIMAL 
// Prototype for Fortran procedure ADVECTFRAC_SLC_SUBOPTIMAL ...
//
void FORTRAN_NAME( ADVECTFRAC_SLC_SUBOPTIMAL ,advectfrac_slc_suboptimal )(
      CHFp_FRA1(frac)
      ,CHFp_FRA1(dh)
      ,CHFp_CONST_FRA1(oldfrac)
      ,CHFp_CONST_FRA1(faceIceVel)
      ,CHFp_CONST_FRA1(faceCalvVel)
      ,CHFp_CONST_FRA1(faceIceFlux)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_REAL(dx)
      ,CHFp_REAL(dt)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir) );

#define FORT_ADVECTFRAC_SLC_SUBOPTIMAL FORTRAN_NAME( inlineADVECTFRAC_SLC_SUBOPTIMAL, inlineADVECTFRAC_SLC_SUBOPTIMAL)
#define FORTNT_ADVECTFRAC_SLC_SUBOPTIMAL FORTRAN_NAME( ADVECTFRAC_SLC_SUBOPTIMAL, advectfrac_slc_suboptimal)

inline void FORTRAN_NAME(inlineADVECTFRAC_SLC_SUBOPTIMAL, inlineADVECTFRAC_SLC_SUBOPTIMAL)(
      CHFp_FRA1(frac)
      ,CHFp_FRA1(dh)
      ,CHFp_CONST_FRA1(oldfrac)
      ,CHFp_CONST_FRA1(faceIceVel)
      ,CHFp_CONST_FRA1(faceCalvVel)
      ,CHFp_CONST_FRA1(faceIceFlux)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_REAL(dx)
      ,CHFp_REAL(dt)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_ADVECTFRAC_SLC_SUBOPTIMAL");
 FORTRAN_NAME( ADVECTFRAC_SLC_SUBOPTIMAL ,advectfrac_slc_suboptimal )(
      CHFt_FRA1(frac)
      ,CHFt_FRA1(dh)
      ,CHFt_CONST_FRA1(oldfrac)
      ,CHFt_CONST_FRA1(faceIceVel)
      ,CHFt_CONST_FRA1(faceCalvVel)
      ,CHFt_CONST_FRA1(faceIceFlux)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_REAL(dx)
      ,CHFt_REAL(dt)
      ,CHFt_REAL(eps)
      ,CHFt_BOX(box)
      ,CHFt_INT(idir) );
}
#endif  // GUARDADVECTFRAC_SLC_SUBOPTIMAL 

#ifndef GUARDABLATERATE_SLC_SUBOPTIMAL 
#define GUARDABLATERATE_SLC_SUBOPTIMAL 
// Prototype for Fortran procedure ablaterate_slc_suboptimal ...
//
void FORTRAN_NAME( ABLATERATE_SLC_SUBOPTIMAL ,ablaterate_slc_suboptimal )(
      CHFp_FRA1(facerate)
      ,CHFp_CONST_FRA1(frac)
      ,CHFp_CONST_FRA1(rate)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir) );

#define FORT_ABLATERATE_SLC_SUBOPTIMAL FORTRAN_NAME( inlineABLATERATE_SLC_SUBOPTIMAL, inlineABLATERATE_SLC_SUBOPTIMAL)
#define FORTNT_ABLATERATE_SLC_SUBOPTIMAL FORTRAN_NAME( ABLATERATE_SLC_SUBOPTIMAL, ablaterate_slc_suboptimal)

inline void FORTRAN_NAME(inlineABLATERATE_SLC_SUBOPTIMAL, inlineABLATERATE_SLC_SUBOPTIMAL)(
      CHFp_FRA1(facerate)
      ,CHFp_CONST_FRA1(frac)
      ,CHFp_CONST_FRA1(rate)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_ABLATERATE_SLC_SUBOPTIMAL");
 FORTRAN_NAME( ABLATERATE_SLC_SUBOPTIMAL ,ablaterate_slc_suboptimal )(
      CHFt_FRA1(facerate)
      ,CHFt_CONST_FRA1(frac)
      ,CHFt_CONST_FRA1(rate)
      ,CHFt_REAL(eps)
      ,CHFt_BOX(box)
      ,CHFt_INT(idir) );
}
#endif  // GUARDABLATERATE_SLC_SUBOPTIMAL 

#ifndef GUARDABLATERATECC_SLC_SUBOPTIMAL 
#define GUARDABLATERATECC_SLC_SUBOPTIMAL 
// Prototype for Fortran procedure ablateratecc_slc_suboptimal ...
//
void FORTRAN_NAME( ABLATERATECC_SLC_SUBOPTIMAL ,ablateratecc_slc_suboptimal )(
      CHFp_FRA1(facerate)
      ,CHFp_CONST_FRA1(frac)
      ,CHFp_CONST_FRA1(rate)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir) );

#define FORT_ABLATERATECC_SLC_SUBOPTIMAL FORTRAN_NAME( inlineABLATERATECC_SLC_SUBOPTIMAL, inlineABLATERATECC_SLC_SUBOPTIMAL)
#define FORTNT_ABLATERATECC_SLC_SUBOPTIMAL FORTRAN_NAME( ABLATERATECC_SLC_SUBOPTIMAL, ablateratecc_slc_suboptimal)

inline void FORTRAN_NAME(inlineABLATERATECC_SLC_SUBOPTIMAL, inlineABLATERATECC_SLC_SUBOPTIMAL)(
      CHFp_FRA1(facerate)
      ,CHFp_CONST_FRA1(frac)
      ,CHFp_CONST_FRA1(rate)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_ABLATERATECC_SLC_SUBOPTIMAL");
 FORTRAN_NAME( ABLATERATECC_SLC_SUBOPTIMAL ,ablateratecc_slc_suboptimal )(
      CHFt_FRA1(facerate)
      ,CHFt_CONST_FRA1(frac)
      ,CHFt_CONST_FRA1(rate)
      ,CHFt_REAL(eps)
      ,CHFt_BOX(box)
      ,CHFt_INT(idir) );
}
#endif  // GUARDABLATERATECC_SLC_SUBOPTIMAL 

#ifndef GUARDFINDFRONTVEL 
#define GUARDFINDFRONTVEL 
// Prototype for Fortran procedure FINDFRONTVEL ...
//
void FORTRAN_NAME( FINDFRONTVEL ,findfrontvel )(
      CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(uFaceVel)
      ,CHFp_CONST_FRA1(vFaceVel)
      ,CHFp_CONST_FRA1(cRate)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA(ccVel)
      ,CHFp_FRA(tmpVel)
      ,CHFp_FRA(relativeVel)
      ,CHFp_CONST_INT(normalCalving)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iStep)
      ,CHFp_BOX(box) );

#define FORT_FINDFRONTVEL FORTRAN_NAME( inlineFINDFRONTVEL, inlineFINDFRONTVEL)
#define FORTNT_FINDFRONTVEL FORTRAN_NAME( FINDFRONTVEL, findfrontvel)

inline void FORTRAN_NAME(inlineFINDFRONTVEL, inlineFINDFRONTVEL)(
      CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(uFaceVel)
      ,CHFp_CONST_FRA1(vFaceVel)
      ,CHFp_CONST_FRA1(cRate)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA(ccVel)
      ,CHFp_FRA(tmpVel)
      ,CHFp_FRA(relativeVel)
      ,CHFp_CONST_INT(normalCalving)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iStep)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FINDFRONTVEL");
 FORTRAN_NAME( FINDFRONTVEL ,findfrontvel )(
      CHFt_CONST_FRA1(oldFrac)
      ,CHFt_CONST_FRA1(uFaceVel)
      ,CHFt_CONST_FRA1(vFaceVel)
      ,CHFt_CONST_FRA1(cRate)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_FRA(ccVel)
      ,CHFt_FRA(tmpVel)
      ,CHFt_FRA(relativeVel)
      ,CHFt_CONST_INT(normalCalving)
      ,CHFt_CONST_REAL(epsFrac)
      ,CHFt_CONST_REAL(epsVel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iStep)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFINDFRONTVEL 

#ifndef GUARDADVECTFRAC 
#define GUARDADVECTFRAC 
// Prototype for Fortran procedure ADVECTFRAC ...
//
void FORTRAN_NAME( ADVECTFRAC ,advectfrac )(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(relativeFrontVel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(idir) );

#define FORT_ADVECTFRAC FORTRAN_NAME( inlineADVECTFRAC, inlineADVECTFRAC)
#define FORTNT_ADVECTFRAC FORTRAN_NAME( ADVECTFRAC, advectfrac)

inline void FORTRAN_NAME(inlineADVECTFRAC, inlineADVECTFRAC)(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(relativeFrontVel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_ADVECTFRAC");
 FORTRAN_NAME( ADVECTFRAC ,advectfrac )(
      CHFt_FRA1(frac)
      ,CHFt_CONST_FRA1(oldFrac)
      ,CHFt_CONST_FRA1(relativeFrontVel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(epsFrac)
      ,CHFt_CONST_REAL(epsVel)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDADVECTFRAC 

#ifndef GUARDADVECTFRACSIMPLE 
#define GUARDADVECTFRACSIMPLE 
// Prototype for Fortran procedure ADVECTFRACSIMPLE ...
//
void FORTRAN_NAME( ADVECTFRACSIMPLE ,advectfracsimple )(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(relativeFrontVel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(idir) );

#define FORT_ADVECTFRACSIMPLE FORTRAN_NAME( inlineADVECTFRACSIMPLE, inlineADVECTFRACSIMPLE)
#define FORTNT_ADVECTFRACSIMPLE FORTRAN_NAME( ADVECTFRACSIMPLE, advectfracsimple)

inline void FORTRAN_NAME(inlineADVECTFRACSIMPLE, inlineADVECTFRACSIMPLE)(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(relativeFrontVel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_ADVECTFRACSIMPLE");
 FORTRAN_NAME( ADVECTFRACSIMPLE ,advectfracsimple )(
      CHFt_FRA1(frac)
      ,CHFt_CONST_FRA1(oldFrac)
      ,CHFt_CONST_FRA1(relativeFrontVel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(epsFrac)
      ,CHFt_CONST_REAL(epsVel)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDADVECTFRACSIMPLE 

#ifndef GUARDGETFRAC 
#define GUARDGETFRAC 
// Prototype for Fortran procedure GETFRAC ...
//
void FORTRAN_NAME( GETFRAC ,getfrac )(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(copyFrac)
      ,CHFp_CONST_FRA(relativeVel)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_CONST_INT(iStep)
      ,CHFp_CONST_INT(icount)
      ,CHFp_INT(iOutofRange)
      ,CHFp_CONST_INT(lev)
      ,CHFp_CONST_INT(MaxIter)
      ,CHFp_BOX(box) );

#define FORT_GETFRAC FORTRAN_NAME( inlineGETFRAC, inlineGETFRAC)
#define FORTNT_GETFRAC FORTRAN_NAME( GETFRAC, getfrac)

inline void FORTRAN_NAME(inlineGETFRAC, inlineGETFRAC)(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA1(copyFrac)
      ,CHFp_CONST_FRA(relativeVel)
      ,CHFp_CONST_REAL(epsFrac)
      ,CHFp_CONST_REAL(epsVel)
      ,CHFp_CONST_INT(iStep)
      ,CHFp_CONST_INT(icount)
      ,CHFp_INT(iOutofRange)
      ,CHFp_CONST_INT(lev)
      ,CHFp_CONST_INT(MaxIter)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETFRAC");
 FORTRAN_NAME( GETFRAC ,getfrac )(
      CHFt_FRA1(frac)
      ,CHFt_CONST_FRA1(oldFrac)
      ,CHFt_CONST_FRA1(copyFrac)
      ,CHFt_CONST_FRA(relativeVel)
      ,CHFt_CONST_REAL(epsFrac)
      ,CHFt_CONST_REAL(epsVel)
      ,CHFt_CONST_INT(iStep)
      ,CHFt_CONST_INT(icount)
      ,CHFt_INT(iOutofRange)
      ,CHFt_CONST_INT(lev)
      ,CHFt_CONST_INT(MaxIter)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETFRAC 

#ifndef GUARDISOLATEDFRAC 
#define GUARDISOLATEDFRAC 
// Prototype for Fortran procedure ISOLATEDFRAC ...
//
void FORTRAN_NAME( ISOLATEDFRAC ,isolatedfrac )(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA(relativeVel)
      ,CHFp_CONST_INT(iStep)
      ,CHFp_CONST_INT(lev)
      ,CHFp_BOX(box) );

#define FORT_ISOLATEDFRAC FORTRAN_NAME( inlineISOLATEDFRAC, inlineISOLATEDFRAC)
#define FORTNT_ISOLATEDFRAC FORTRAN_NAME( ISOLATEDFRAC, isolatedfrac)

inline void FORTRAN_NAME(inlineISOLATEDFRAC, inlineISOLATEDFRAC)(
      CHFp_FRA1(frac)
      ,CHFp_CONST_FRA1(oldFrac)
      ,CHFp_CONST_FRA(relativeVel)
      ,CHFp_CONST_INT(iStep)
      ,CHFp_CONST_INT(lev)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ISOLATEDFRAC");
 FORTRAN_NAME( ISOLATEDFRAC ,isolatedfrac )(
      CHFt_FRA1(frac)
      ,CHFt_CONST_FRA1(oldFrac)
      ,CHFt_CONST_FRA(relativeVel)
      ,CHFt_CONST_INT(iStep)
      ,CHFt_CONST_INT(lev)
      ,CHFt_BOX(box) );
}
#endif  // GUARDISOLATEDFRAC 

#ifndef GUARDCOMPUTEZVEL 
#define GUARDCOMPUTEZVEL 
// Prototype for Fortran procedure COMPUTEZVEL ...
//
void FORTRAN_NAME( COMPUTEZVEL ,computezvel )(
      CHFp_FRA(uz)
      ,CHFp_FRA1(uzs)
      ,CHFp_CONST_FRA(ux)
      ,CHFp_CONST_FRA(uy)
      ,CHFp_CONST_FRA(divuhxy)
      ,CHFp_CONST_VR(fsig)
      ,CHFp_CONST_VR(csig)
      ,CHFp_CONST_VR(dsig)
      ,CHFp_CONST_FRA1(dsx)
      ,CHFp_CONST_FRA1(dhx)
      ,CHFp_CONST_FRA1(dsy)
      ,CHFp_CONST_FRA1(dhy)
      ,CHFp_CONST_FRA1(dst)
      ,CHFp_CONST_FRA1(dht)
      ,CHFp_CONST_FRA1(smb)
      ,CHFp_CONST_FRA1(bmb)
      ,CHFp_CONST_INT(nlay)
      ,CHFp_BOX(box) );

#define FORT_COMPUTEZVEL FORTRAN_NAME( inlineCOMPUTEZVEL, inlineCOMPUTEZVEL)
#define FORTNT_COMPUTEZVEL FORTRAN_NAME( COMPUTEZVEL, computezvel)

inline void FORTRAN_NAME(inlineCOMPUTEZVEL, inlineCOMPUTEZVEL)(
      CHFp_FRA(uz)
      ,CHFp_FRA1(uzs)
      ,CHFp_CONST_FRA(ux)
      ,CHFp_CONST_FRA(uy)
      ,CHFp_CONST_FRA(divuhxy)
      ,CHFp_CONST_VR(fsig)
      ,CHFp_CONST_VR(csig)
      ,CHFp_CONST_VR(dsig)
      ,CHFp_CONST_FRA1(dsx)
      ,CHFp_CONST_FRA1(dhx)
      ,CHFp_CONST_FRA1(dsy)
      ,CHFp_CONST_FRA1(dhy)
      ,CHFp_CONST_FRA1(dst)
      ,CHFp_CONST_FRA1(dht)
      ,CHFp_CONST_FRA1(smb)
      ,CHFp_CONST_FRA1(bmb)
      ,CHFp_CONST_INT(nlay)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COMPUTEZVEL");
 FORTRAN_NAME( COMPUTEZVEL ,computezvel )(
      CHFt_FRA(uz)
      ,CHFt_FRA1(uzs)
      ,CHFt_CONST_FRA(ux)
      ,CHFt_CONST_FRA(uy)
      ,CHFt_CONST_FRA(divuhxy)
      ,CHFt_CONST_VR(fsig)
      ,CHFt_CONST_VR(csig)
      ,CHFt_CONST_VR(dsig)
      ,CHFt_CONST_FRA1(dsx)
      ,CHFt_CONST_FRA1(dhx)
      ,CHFt_CONST_FRA1(dsy)
      ,CHFt_CONST_FRA1(dhy)
      ,CHFt_CONST_FRA1(dst)
      ,CHFt_CONST_FRA1(dht)
      ,CHFt_CONST_FRA1(smb)
      ,CHFt_CONST_FRA1(bmb)
      ,CHFt_CONST_INT(nlay)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOMPUTEZVEL 

#ifndef GUARDCOMPUTESIGMAVEL 
#define GUARDCOMPUTESIGMAVEL 
// Prototype for Fortran procedure COMPUTESIGMAVEL ...
//
void FORTRAN_NAME( COMPUTESIGMAVEL ,computesigmavel )(
      CHFp_FRA(usig)
      ,CHFp_CONST_FRA(ux)
      ,CHFp_CONST_FRA(uy)
      ,CHFp_CONST_FRA(divuhxy)
      ,CHFp_CONST_VR(dsig)
      ,CHFp_CONST_FRA1(dht)
      ,CHFp_CONST_FRA1(smb)
      ,CHFp_CONST_FRA1(bmb)
      ,CHFp_CONST_INT(nlay)
      ,CHFp_BOX(box) );

#define FORT_COMPUTESIGMAVEL FORTRAN_NAME( inlineCOMPUTESIGMAVEL, inlineCOMPUTESIGMAVEL)
#define FORTNT_COMPUTESIGMAVEL FORTRAN_NAME( COMPUTESIGMAVEL, computesigmavel)

inline void FORTRAN_NAME(inlineCOMPUTESIGMAVEL, inlineCOMPUTESIGMAVEL)(
      CHFp_FRA(usig)
      ,CHFp_CONST_FRA(ux)
      ,CHFp_CONST_FRA(uy)
      ,CHFp_CONST_FRA(divuhxy)
      ,CHFp_CONST_VR(dsig)
      ,CHFp_CONST_FRA1(dht)
      ,CHFp_CONST_FRA1(smb)
      ,CHFp_CONST_FRA1(bmb)
      ,CHFp_CONST_INT(nlay)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COMPUTESIGMAVEL");
 FORTRAN_NAME( COMPUTESIGMAVEL ,computesigmavel )(
      CHFt_FRA(usig)
      ,CHFt_CONST_FRA(ux)
      ,CHFt_CONST_FRA(uy)
      ,CHFt_CONST_FRA(divuhxy)
      ,CHFt_CONST_VR(dsig)
      ,CHFt_CONST_FRA1(dht)
      ,CHFt_CONST_FRA1(smb)
      ,CHFt_CONST_FRA1(bmb)
      ,CHFt_CONST_INT(nlay)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOMPUTESIGMAVEL 

#ifndef GUARDEVOLVEGROUNDEDBED 
#define GUARDEVOLVEGROUNDEDBED 
// Prototype for Fortran procedure evolvegroundedbed ...
//
void FORTRAN_NAME( EVOLVEGROUNDEDBED ,evolvegroundedbed )(
      CHFp_FRA1(newh)
      ,CHFp_FRA1(oldh)
      ,CHFp_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(box) );

#define FORT_EVOLVEGROUNDEDBED FORTRAN_NAME( inlineEVOLVEGROUNDEDBED, inlineEVOLVEGROUNDEDBED)
#define FORTNT_EVOLVEGROUNDEDBED FORTRAN_NAME( EVOLVEGROUNDEDBED, evolvegroundedbed)

inline void FORTRAN_NAME(inlineEVOLVEGROUNDEDBED, inlineEVOLVEGROUNDEDBED)(
      CHFp_FRA1(newh)
      ,CHFp_FRA1(oldh)
      ,CHFp_FRA1(topg)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EVOLVEGROUNDEDBED");
 FORTRAN_NAME( EVOLVEGROUNDEDBED ,evolvegroundedbed )(
      CHFt_FRA1(newh)
      ,CHFt_FRA1(oldh)
      ,CHFt_FRA1(topg)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEVOLVEGROUNDEDBED 

}

#endif
