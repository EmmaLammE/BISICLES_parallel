#ifndef _CREVASSEF_F_H_
#define _CREVASSEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVDVSTRESSB 
#define GUARDVDVSTRESSB 
// Prototype for Fortran procedure VDVSTRESSB ...
//
void FORTRAN_NAME( VDVSTRESSB ,vdvstressb )(
      CHFp_FRA1(k)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(hp)
      ,CHFp_CONST_FRA1(d)
      ,CHFp_CONST_FRA1(rxx)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhow)
      ,CHFp_CONST_REAL(gravity)
      ,CHFp_BOX(box) );

#define FORT_VDVSTRESSB FORTRAN_NAME( inlineVDVSTRESSB, inlineVDVSTRESSB)
#define FORTNT_VDVSTRESSB FORTRAN_NAME( VDVSTRESSB, vdvstressb)

inline void FORTRAN_NAME(inlineVDVSTRESSB, inlineVDVSTRESSB)(
      CHFp_FRA1(k)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(hp)
      ,CHFp_CONST_FRA1(d)
      ,CHFp_CONST_FRA1(rxx)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhow)
      ,CHFp_CONST_REAL(gravity)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VDVSTRESSB");
 FORTRAN_NAME( VDVSTRESSB ,vdvstressb )(
      CHFt_FRA1(k)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_FRA1(hp)
      ,CHFt_CONST_FRA1(d)
      ,CHFt_CONST_FRA1(rxx)
      ,CHFt_CONST_REAL(rhoi)
      ,CHFt_CONST_REAL(rhow)
      ,CHFt_CONST_REAL(gravity)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVDVSTRESSB 

#ifndef GUARDVDVSTRESSS 
#define GUARDVDVSTRESSS 
// Prototype for Fortran procedure VDVSTRESSS ...
//
void FORTRAN_NAME( VDVSTRESSS ,vdvstresss )(
      CHFp_FRA1(k)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(dw)
      ,CHFp_CONST_FRA1(d)
      ,CHFp_CONST_FRA1(rxx)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhow)
      ,CHFp_CONST_REAL(gravity)
      ,CHFp_BOX(box) );

#define FORT_VDVSTRESSS FORTRAN_NAME( inlineVDVSTRESSS, inlineVDVSTRESSS)
#define FORTNT_VDVSTRESSS FORTRAN_NAME( VDVSTRESSS, vdvstresss)

inline void FORTRAN_NAME(inlineVDVSTRESSS, inlineVDVSTRESSS)(
      CHFp_FRA1(k)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(dw)
      ,CHFp_CONST_FRA1(d)
      ,CHFp_CONST_FRA1(rxx)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhow)
      ,CHFp_CONST_REAL(gravity)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VDVSTRESSS");
 FORTRAN_NAME( VDVSTRESSS ,vdvstresss )(
      CHFt_FRA1(k)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_FRA1(dw)
      ,CHFt_CONST_FRA1(d)
      ,CHFt_CONST_FRA1(rxx)
      ,CHFt_CONST_REAL(rhoi)
      ,CHFt_CONST_REAL(rhow)
      ,CHFt_CONST_REAL(gravity)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVDVSTRESSS 

}

#endif
