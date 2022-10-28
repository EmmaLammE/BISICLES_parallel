#ifndef _CALVINGF_F_H_
#define _CALVINGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEFFECTIVETHICKNESS 
#define GUARDEFFECTIVETHICKNESS 
// Prototype for Fortran procedure effectivethickness ...
//
void FORTRAN_NAME( EFFECTIVETHICKNESS ,effectivethickness )(
      CHFp_FRA1(he)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(m)
      ,CHFp_BOX(box) );

#define FORT_EFFECTIVETHICKNESS FORTRAN_NAME( inlineEFFECTIVETHICKNESS, inlineEFFECTIVETHICKNESS)
#define FORTNT_EFFECTIVETHICKNESS FORTRAN_NAME( EFFECTIVETHICKNESS, effectivethickness)

inline void FORTRAN_NAME(inlineEFFECTIVETHICKNESS, inlineEFFECTIVETHICKNESS)(
      CHFp_FRA1(he)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(m)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EFFECTIVETHICKNESS");
 FORTRAN_NAME( EFFECTIVETHICKNESS ,effectivethickness )(
      CHFt_FRA1(he)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_FRA1(m)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEFFECTIVETHICKNESS 

}

#endif
