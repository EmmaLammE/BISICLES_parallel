#ifndef _BISICLESF_F_H_
#define _BISICLESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDENFORCEWELLPOSEDCELL 
#define GUARDENFORCEWELLPOSEDCELL 
// Prototype for Fortran procedure enforcewellposedcell ...
//
void FORTRAN_NAME( ENFORCEWELLPOSEDCELL ,enforcewellposedcell )(
      CHFp_FRA1(C)
      ,CHFp_CONST_FRA1(mux)
      ,CHFp_CONST_FRA1(muy)
      ,CHFp_CONST_REAL(mu0)
      ,CHFp_CONST_REAL(C0)
      ,CHFp_BOX(box) );

#define FORT_ENFORCEWELLPOSEDCELL FORTRAN_NAME( inlineENFORCEWELLPOSEDCELL, inlineENFORCEWELLPOSEDCELL)
#define FORTNT_ENFORCEWELLPOSEDCELL FORTRAN_NAME( ENFORCEWELLPOSEDCELL, enforcewellposedcell)

inline void FORTRAN_NAME(inlineENFORCEWELLPOSEDCELL, inlineENFORCEWELLPOSEDCELL)(
      CHFp_FRA1(C)
      ,CHFp_CONST_FRA1(mux)
      ,CHFp_CONST_FRA1(muy)
      ,CHFp_CONST_REAL(mu0)
      ,CHFp_CONST_REAL(C0)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ENFORCEWELLPOSEDCELL");
 FORTRAN_NAME( ENFORCEWELLPOSEDCELL ,enforcewellposedcell )(
      CHFt_FRA1(C)
      ,CHFt_CONST_FRA1(mux)
      ,CHFt_CONST_FRA1(muy)
      ,CHFt_CONST_REAL(mu0)
      ,CHFt_CONST_REAL(C0)
      ,CHFt_BOX(box) );
}
#endif  // GUARDENFORCEWELLPOSEDCELL 

#ifndef GUARDMASKEDREPLACE 
#define GUARDMASKEDREPLACE 
// Prototype for Fortran procedure maskedreplace ...
//
void FORTRAN_NAME( MASKEDREPLACE ,maskedreplace )(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_INT(m)
      ,CHFp_BOX(box) );

#define FORT_MASKEDREPLACE FORTRAN_NAME( inlineMASKEDREPLACE, inlineMASKEDREPLACE)
#define FORTNT_MASKEDREPLACE FORTRAN_NAME( MASKEDREPLACE, maskedreplace)

inline void FORTRAN_NAME(inlineMASKEDREPLACE, inlineMASKEDREPLACE)(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_INT(m)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MASKEDREPLACE");
 FORTRAN_NAME( MASKEDREPLACE ,maskedreplace )(
      CHFt_FRA1(a)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_INT(m)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMASKEDREPLACE 

#ifndef GUARDPROXIMITYFILL 
#define GUARDPROXIMITYFILL 
// Prototype for Fortran procedure proximityfill ...
//
void FORTRAN_NAME( PROXIMITYFILL ,proximityfill )(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_FRA1(p)
      ,CHFp_CONST_REAL(n)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_REAL(m)
      ,CHFp_BOX(box) );

#define FORT_PROXIMITYFILL FORTRAN_NAME( inlinePROXIMITYFILL, inlinePROXIMITYFILL)
#define FORTNT_PROXIMITYFILL FORTRAN_NAME( PROXIMITYFILL, proximityfill)

inline void FORTRAN_NAME(inlinePROXIMITYFILL, inlinePROXIMITYFILL)(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_CONST_FRA1(p)
      ,CHFp_CONST_REAL(n)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_REAL(m)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PROXIMITYFILL");
 FORTRAN_NAME( PROXIMITYFILL ,proximityfill )(
      CHFt_FRA1(a)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_CONST_FRA1(p)
      ,CHFt_CONST_REAL(n)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_REAL(m)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPROXIMITYFILL 

#ifndef GUARDWATERDEPTH 
#define GUARDWATERDEPTH 
// Prototype for Fortran procedure waterdepth ...
//
void FORTRAN_NAME( WATERDEPTH ,waterdepth )(
      CHFp_FRA1(depth)
      ,CHFp_CONST_FRA1(thk)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_BOX(box) );

#define FORT_WATERDEPTH FORTRAN_NAME( inlineWATERDEPTH, inlineWATERDEPTH)
#define FORTNT_WATERDEPTH FORTRAN_NAME( WATERDEPTH, waterdepth)

inline void FORTRAN_NAME(inlineWATERDEPTH, inlineWATERDEPTH)(
      CHFp_FRA1(depth)
      ,CHFp_CONST_FRA1(thk)
      ,CHFp_CONST_FRA1(usrf)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_WATERDEPTH");
 FORTRAN_NAME( WATERDEPTH ,waterdepth )(
      CHFt_FRA1(depth)
      ,CHFt_CONST_FRA1(thk)
      ,CHFt_CONST_FRA1(usrf)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_BOX(box) );
}
#endif  // GUARDWATERDEPTH 

#ifndef GUARDZEROIFLESS 
#define GUARDZEROIFLESS 
// Prototype for Fortran procedure zeroifless ...
//
void FORTRAN_NAME( ZEROIFLESS ,zeroifless )(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(infab)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_BOX(box) );

#define FORT_ZEROIFLESS FORTRAN_NAME( inlineZEROIFLESS, inlineZEROIFLESS)
#define FORTNT_ZEROIFLESS FORTRAN_NAME( ZEROIFLESS, zeroifless)

inline void FORTRAN_NAME(inlineZEROIFLESS, inlineZEROIFLESS)(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(infab)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ZEROIFLESS");
 FORTRAN_NAME( ZEROIFLESS ,zeroifless )(
      CHFt_FRA1(fab)
      ,CHFt_CONST_FRA1(infab)
      ,CHFt_CONST_REAL(tol)
      ,CHFt_BOX(box) );
}
#endif  // GUARDZEROIFLESS 

#ifndef GUARDPWLFILL 
#define GUARDPWLFILL 
// Prototype for Fortran procedure pwlfill ...
//
void FORTRAN_NAME( PWLFILL ,pwlfill )(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_VR(xn)
      ,CHFp_CONST_VR(bn)
      ,CHFp_VR(dx)
      ,CHFp_VR(db)
      ,CHFp_BOX(box) );

#define FORT_PWLFILL FORTRAN_NAME( inlinePWLFILL, inlinePWLFILL)
#define FORTNT_PWLFILL FORTRAN_NAME( PWLFILL, pwlfill)

inline void FORTRAN_NAME(inlinePWLFILL, inlinePWLFILL)(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_VR(xn)
      ,CHFp_CONST_VR(bn)
      ,CHFp_VR(dx)
      ,CHFp_VR(db)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PWLFILL");
 FORTRAN_NAME( PWLFILL ,pwlfill )(
      CHFt_FRA1(a)
      ,CHFt_CONST_FRA1(x)
      ,CHFt_CONST_VR(xn)
      ,CHFt_CONST_VR(bn)
      ,CHFt_VR(dx)
      ,CHFt_VR(db)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPWLFILL 

#ifndef GUARDABSLIMITFAB 
#define GUARDABSLIMITFAB 
// Prototype for Fortran procedure abslimitfab ...
//
void FORTRAN_NAME( ABSLIMITFAB ,abslimitfab )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) );

#define FORT_ABSLIMITFAB FORTRAN_NAME( inlineABSLIMITFAB, inlineABSLIMITFAB)
#define FORTNT_ABSLIMITFAB FORTRAN_NAME( ABSLIMITFAB, abslimitfab)

inline void FORTRAN_NAME(inlineABSLIMITFAB, inlineABSLIMITFAB)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ABSLIMITFAB");
 FORTRAN_NAME( ABSLIMITFAB ,abslimitfab )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(limit)
      ,CHFt_BOX(box) );
}
#endif  // GUARDABSLIMITFAB 

#ifndef GUARDMODLIMITFAB 
#define GUARDMODLIMITFAB 
// Prototype for Fortran procedure modlimitfab ...
//
void FORTRAN_NAME( MODLIMITFAB ,modlimitfab )(
      CHFp_FRA(u)
      ,CHFp_INT(nlimit)
      ,CHFp_REAL(umax)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) );

#define FORT_MODLIMITFAB FORTRAN_NAME( inlineMODLIMITFAB, inlineMODLIMITFAB)
#define FORTNT_MODLIMITFAB FORTRAN_NAME( MODLIMITFAB, modlimitfab)

inline void FORTRAN_NAME(inlineMODLIMITFAB, inlineMODLIMITFAB)(
      CHFp_FRA(u)
      ,CHFp_INT(nlimit)
      ,CHFp_REAL(umax)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MODLIMITFAB");
 FORTRAN_NAME( MODLIMITFAB ,modlimitfab )(
      CHFt_FRA(u)
      ,CHFt_INT(nlimit)
      ,CHFt_REAL(umax)
      ,CHFt_CONST_REAL(limit)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMODLIMITFAB 

#ifndef GUARDMINFAB1 
#define GUARDMINFAB1 
// Prototype for Fortran procedure minfab1 ...
//
void FORTRAN_NAME( MINFAB1 ,minfab1 )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) );

#define FORT_MINFAB1 FORTRAN_NAME( inlineMINFAB1, inlineMINFAB1)
#define FORTNT_MINFAB1 FORTRAN_NAME( MINFAB1, minfab1)

inline void FORTRAN_NAME(inlineMINFAB1, inlineMINFAB1)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINFAB1");
 FORTRAN_NAME( MINFAB1 ,minfab1 )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(limit)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINFAB1 

#ifndef GUARDMAXFAB1 
#define GUARDMAXFAB1 
// Prototype for Fortran procedure maxfab1 ...
//
void FORTRAN_NAME( MAXFAB1 ,maxfab1 )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) );

#define FORT_MAXFAB1 FORTRAN_NAME( inlineMAXFAB1, inlineMAXFAB1)
#define FORTNT_MAXFAB1 FORTRAN_NAME( MAXFAB1, maxfab1)

inline void FORTRAN_NAME(inlineMAXFAB1, inlineMAXFAB1)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXFAB1");
 FORTRAN_NAME( MAXFAB1 ,maxfab1 )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(limit)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXFAB1 

#ifndef GUARDCHECKCOEF 
#define GUARDCHECKCOEF 
// Prototype for Fortran procedure checkcoef ...
//
void FORTRAN_NAME( CHECKCOEF ,checkcoef )(
      CHFp_FIA1(ok)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_CONST_FRA1(alpha)
      ,CHFp_CONST_REAL(mumin)
      ,CHFp_CONST_REAL(alphamin)
      ,CHFp_BOX(box) );

#define FORT_CHECKCOEF FORTRAN_NAME( inlineCHECKCOEF, inlineCHECKCOEF)
#define FORTNT_CHECKCOEF FORTRAN_NAME( CHECKCOEF, checkcoef)

inline void FORTRAN_NAME(inlineCHECKCOEF, inlineCHECKCOEF)(
      CHFp_FIA1(ok)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_CONST_FRA1(alpha)
      ,CHFp_CONST_REAL(mumin)
      ,CHFp_CONST_REAL(alphamin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHECKCOEF");
 FORTRAN_NAME( CHECKCOEF ,checkcoef )(
      CHFt_FIA1(ok)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(mu)
      ,CHFt_CONST_FRA1(alpha)
      ,CHFt_CONST_REAL(mumin)
      ,CHFt_CONST_REAL(alphamin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHECKCOEF 

#ifndef GUARDAPPLYMINVISCOSITY 
#define GUARDAPPLYMINVISCOSITY 
// Prototype for Fortran procedure applyminviscosity ...
//
void FORTRAN_NAME( APPLYMINVISCOSITY ,applyminviscosity )(
      CHFp_FRA1(mu)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(faceH)
      ,CHFp_CONST_FIA1(ok)
      ,CHFp_CONST_REAL(mumin)
      ,CHFp_BOX(box) );

#define FORT_APPLYMINVISCOSITY FORTRAN_NAME( inlineAPPLYMINVISCOSITY, inlineAPPLYMINVISCOSITY)
#define FORTNT_APPLYMINVISCOSITY FORTRAN_NAME( APPLYMINVISCOSITY, applyminviscosity)

inline void FORTRAN_NAME(inlineAPPLYMINVISCOSITY, inlineAPPLYMINVISCOSITY)(
      CHFp_FRA1(mu)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(faceH)
      ,CHFp_CONST_FIA1(ok)
      ,CHFp_CONST_REAL(mumin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_APPLYMINVISCOSITY");
 FORTRAN_NAME( APPLYMINVISCOSITY ,applyminviscosity )(
      CHFt_FRA1(mu)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(faceH)
      ,CHFt_CONST_FIA1(ok)
      ,CHFt_CONST_REAL(mumin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAPPLYMINVISCOSITY 

#ifndef GUARDUPWINDFLUXB 
#define GUARDUPWINDFLUXB 
// Prototype for Fortran procedure upwindfluxb ...
//
void FORTRAN_NAME( UPWINDFLUXB ,upwindfluxb )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_FRA1(s)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(fbox) );

#define FORT_UPWINDFLUXB FORTRAN_NAME( inlineUPWINDFLUXB, inlineUPWINDFLUXB)
#define FORTNT_UPWINDFLUXB FORTRAN_NAME( UPWINDFLUXB, upwindfluxb)

inline void FORTRAN_NAME(inlineUPWINDFLUXB, inlineUPWINDFLUXB)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_FRA1(s)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(fbox) )
{
 CH_TIMELEAF("FORT_UPWINDFLUXB");
 FORTRAN_NAME( UPWINDFLUXB ,upwindfluxb )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_FRA1(s)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(fbox) );
}
#endif  // GUARDUPWINDFLUXB 

#ifndef GUARDHOVERUNDERFLOTATION 
#define GUARDHOVERUNDERFLOTATION 
// Prototype for Fortran procedure hoverunderflotation ...
//
void FORTRAN_NAME( HOVERUNDERFLOTATION ,hoverunderflotation )(
      CHFp_FRA1(hab)
      ,CHFp_REAL(hmin)
      ,CHFp_REAL(hmax)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(sealevel)
      ,CHFp_BOX(box) );

#define FORT_HOVERUNDERFLOTATION FORTRAN_NAME( inlineHOVERUNDERFLOTATION, inlineHOVERUNDERFLOTATION)
#define FORTNT_HOVERUNDERFLOTATION FORTRAN_NAME( HOVERUNDERFLOTATION, hoverunderflotation)

inline void FORTRAN_NAME(inlineHOVERUNDERFLOTATION, inlineHOVERUNDERFLOTATION)(
      CHFp_FRA1(hab)
      ,CHFp_REAL(hmin)
      ,CHFp_REAL(hmax)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(topg)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(sealevel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HOVERUNDERFLOTATION");
 FORTRAN_NAME( HOVERUNDERFLOTATION ,hoverunderflotation )(
      CHFt_FRA1(hab)
      ,CHFt_REAL(hmin)
      ,CHFt_REAL(hmax)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_FRA1(topg)
      ,CHFt_CONST_REAL(rhoi)
      ,CHFt_CONST_REAL(rhoo)
      ,CHFt_CONST_REAL(sealevel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHOVERUNDERFLOTATION 

#ifndef GUARDINTEGRATEHEAVISIDE2D 
#define GUARDINTEGRATEHEAVISIDE2D 
// Prototype for Fortran procedure integrateheaviside2D ...
//
void FORTRAN_NAME( INTEGRATEHEAVISIDE2D ,integrateheaviside2d )(
      CHFp_FRA1(r)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_INT(n)
      ,CHFp_BOX(box) );

#define FORT_INTEGRATEHEAVISIDE2D FORTRAN_NAME( inlineINTEGRATEHEAVISIDE2D, inlineINTEGRATEHEAVISIDE2D)
#define FORTNT_INTEGRATEHEAVISIDE2D FORTRAN_NAME( INTEGRATEHEAVISIDE2D, integrateheaviside2d)

inline void FORTRAN_NAME(inlineINTEGRATEHEAVISIDE2D, inlineINTEGRATEHEAVISIDE2D)(
      CHFp_FRA1(r)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_INT(n)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INTEGRATEHEAVISIDE2D");
 FORTRAN_NAME( INTEGRATEHEAVISIDE2D ,integrateheaviside2d )(
      CHFt_FRA1(r)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_INT(n)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINTEGRATEHEAVISIDE2D 

}

#endif
