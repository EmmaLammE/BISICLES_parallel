#ifndef _DERIVATIVESF_F_H_
#define _DERIVATIVESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCCDERIV 
#define GUARDCCDERIV 
// Prototype for Fortran procedure CCDERIV ...
//
void FORTRAN_NAME( CCDERIV ,ccderiv )(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_BOX(derivBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(derivDir) );

#define FORT_CCDERIV FORTRAN_NAME( inlineCCDERIV, inlineCCDERIV)
#define FORTNT_CCDERIV FORTRAN_NAME( CCDERIV, ccderiv)

inline void FORTRAN_NAME(inlineCCDERIV, inlineCCDERIV)(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_BOX(derivBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(derivDir) )
{
 CH_TIMELEAF("FORT_CCDERIV");
 FORTRAN_NAME( CCDERIV ,ccderiv )(
      CHFt_FRA1(deriv)
      ,CHFt_CONST_FRA1(data)
      ,CHFt_BOX(derivBox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(derivDir) );
}
#endif  // GUARDCCDERIV 

#ifndef GUARDCCDERIVMASK 
#define GUARDCCDERIVMASK 
// Prototype for Fortran procedure CCDERIVMASK ...
//
void FORTRAN_NAME( CCDERIVMASK ,ccderivmask )(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_FRA1(thk)
      ,CHFp_BOX(derivBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(derivDir) );

#define FORT_CCDERIVMASK FORTRAN_NAME( inlineCCDERIVMASK, inlineCCDERIVMASK)
#define FORTNT_CCDERIVMASK FORTRAN_NAME( CCDERIVMASK, ccderivmask)

inline void FORTRAN_NAME(inlineCCDERIVMASK, inlineCCDERIVMASK)(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_FRA1(thk)
      ,CHFp_BOX(derivBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(derivDir) )
{
 CH_TIMELEAF("FORT_CCDERIVMASK");
 FORTRAN_NAME( CCDERIVMASK ,ccderivmask )(
      CHFt_FRA1(deriv)
      ,CHFt_CONST_FRA1(data)
      ,CHFt_CONST_FRA1(thk)
      ,CHFt_BOX(derivBox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(derivDir) );
}
#endif  // GUARDCCDERIVMASK 

#ifndef GUARDFACEDERIV 
#define GUARDFACEDERIV 
// Prototype for Fortran procedure FACEDERIV ...
//
void FORTRAN_NAME( FACEDERIV ,facederiv )(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(ccData)
      ,CHFp_BOX(derivBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(derivDir)
      ,CHFp_INT(faceDir) );

#define FORT_FACEDERIV FORTRAN_NAME( inlineFACEDERIV, inlineFACEDERIV)
#define FORTNT_FACEDERIV FORTRAN_NAME( FACEDERIV, facederiv)

inline void FORTRAN_NAME(inlineFACEDERIV, inlineFACEDERIV)(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(ccData)
      ,CHFp_BOX(derivBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(derivDir)
      ,CHFp_INT(faceDir) )
{
 CH_TIMELEAF("FORT_FACEDERIV");
 FORTRAN_NAME( FACEDERIV ,facederiv )(
      CHFt_FRA1(deriv)
      ,CHFt_CONST_FRA1(ccData)
      ,CHFt_BOX(derivBox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(derivDir)
      ,CHFt_INT(faceDir) );
}
#endif  // GUARDFACEDERIV 

}

#endif
