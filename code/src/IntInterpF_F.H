#ifndef _INTINTERPF_F_H_
#define _INTINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINTINTERPCONSTANT 
#define GUARDINTINTERPCONSTANT 
// Prototype for Fortran procedure intinterpconstant ...
//
void FORTRAN_NAME( INTINTERPCONSTANT ,intinterpconstant )(
      CHFp_FIA(fine)
      ,CHFp_CONST_FIA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) );

#define FORT_INTINTERPCONSTANT FORTRAN_NAME( inlineINTINTERPCONSTANT, inlineINTINTERPCONSTANT)
#define FORTNT_INTINTERPCONSTANT FORTRAN_NAME( INTINTERPCONSTANT, intinterpconstant)

inline void FORTRAN_NAME(inlineINTINTERPCONSTANT, inlineINTINTERPCONSTANT)(
      CHFp_FIA(fine)
      ,CHFp_CONST_FIA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_INTINTERPCONSTANT");
 FORTRAN_NAME( INTINTERPCONSTANT ,intinterpconstant )(
      CHFt_FIA(fine)
      ,CHFt_CONST_FIA(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDINTINTERPCONSTANT 

#ifndef GUARDPROLONGQUAD_ICE 
#define GUARDPROLONGQUAD_ICE 
// Prototype for Fortran procedure PROLONGQUAD_ICE ...
//
void FORTRAN_NAME( PROLONGQUAD_ICE ,prolongquad_ice )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INT(refRatio) );

#define FORT_PROLONGQUAD_ICE FORTRAN_NAME( inlinePROLONGQUAD_ICE, inlinePROLONGQUAD_ICE)
#define FORTNT_PROLONGQUAD_ICE FORTRAN_NAME( PROLONGQUAD_ICE, prolongquad_ice)

inline void FORTRAN_NAME(inlinePROLONGQUAD_ICE, inlinePROLONGQUAD_ICE)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INT(refRatio) )
{
 CH_TIMELEAF("FORT_PROLONGQUAD_ICE");
 FORTRAN_NAME( PROLONGQUAD_ICE ,prolongquad_ice )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(fineBox)
      ,CHFt_CONST_INT(refRatio) );
}
#endif  // GUARDPROLONGQUAD_ICE 

#ifndef GUARDPROLONGLINEAR_ICE 
#define GUARDPROLONGLINEAR_ICE 
// Prototype for Fortran procedure PROLONGLINEAR_ICE ...
//
void FORTRAN_NAME( PROLONGLINEAR_ICE ,prolonglinear_ice )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INT(refRatio) );

#define FORT_PROLONGLINEAR_ICE FORTRAN_NAME( inlinePROLONGLINEAR_ICE, inlinePROLONGLINEAR_ICE)
#define FORTNT_PROLONGLINEAR_ICE FORTRAN_NAME( PROLONGLINEAR_ICE, prolonglinear_ice)

inline void FORTRAN_NAME(inlinePROLONGLINEAR_ICE, inlinePROLONGLINEAR_ICE)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INT(refRatio) )
{
 CH_TIMELEAF("FORT_PROLONGLINEAR_ICE");
 FORTRAN_NAME( PROLONGLINEAR_ICE ,prolonglinear_ice )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(fineBox)
      ,CHFt_CONST_INT(refRatio) );
}
#endif  // GUARDPROLONGLINEAR_ICE 

}

#endif
