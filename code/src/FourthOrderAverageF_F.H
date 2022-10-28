#ifndef _FOURTHORDERAVERAGEF_F_H_
#define _FOURTHORDERAVERAGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINCREMENTLAPLACIAN 
#define GUARDINCREMENTLAPLACIAN 
// Prototype for Fortran procedure INCREMENTLAPLACIAN ...
//
void FORTRAN_NAME( INCREMENTLAPLACIAN ,incrementlaplacian )(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(factor) );

#define FORT_INCREMENTLAPLACIAN FORTRAN_NAME( inlineINCREMENTLAPLACIAN, inlineINCREMENTLAPLACIAN)
#define FORTNT_INCREMENTLAPLACIAN FORTRAN_NAME( INCREMENTLAPLACIAN, incrementlaplacian)

inline void FORTRAN_NAME(inlineINCREMENTLAPLACIAN, inlineINCREMENTLAPLACIAN)(
      CHFp_FRA(lapPhi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(factor) )
{
 CH_TIMELEAF("FORT_INCREMENTLAPLACIAN");
 FORTRAN_NAME( INCREMENTLAPLACIAN ,incrementlaplacian )(
      CHFt_FRA(lapPhi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(gridBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(factor) );
}
#endif  // GUARDINCREMENTLAPLACIAN 

#ifndef GUARDINCREMENTLOSIDELAPLACIAN 
#define GUARDINCREMENTLOSIDELAPLACIAN 
// Prototype for Fortran procedure INCREMENTLOSIDELAPLACIAN ...
//
void FORTRAN_NAME( INCREMENTLOSIDELAPLACIAN ,incrementlosidelaplacian )(
      CHFp_FRA(lapPhi)
      ,CHFp_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_INT(dir)
      ,CHFp_REAL(factor) );

#define FORT_INCREMENTLOSIDELAPLACIAN FORTRAN_NAME( inlineINCREMENTLOSIDELAPLACIAN, inlineINCREMENTLOSIDELAPLACIAN)
#define FORTNT_INCREMENTLOSIDELAPLACIAN FORTRAN_NAME( INCREMENTLOSIDELAPLACIAN, incrementlosidelaplacian)

inline void FORTRAN_NAME(inlineINCREMENTLOSIDELAPLACIAN, inlineINCREMENTLOSIDELAPLACIAN)(
      CHFp_FRA(lapPhi)
      ,CHFp_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_INT(dir)
      ,CHFp_REAL(factor) )
{
 CH_TIMELEAF("FORT_INCREMENTLOSIDELAPLACIAN");
 FORTRAN_NAME( INCREMENTLOSIDELAPLACIAN ,incrementlosidelaplacian )(
      CHFt_FRA(lapPhi)
      ,CHFt_FRA(phi)
      ,CHFt_BOX(gridBox)
      ,CHFt_INT(dir)
      ,CHFt_REAL(factor) );
}
#endif  // GUARDINCREMENTLOSIDELAPLACIAN 

#ifndef GUARDINCREMENTHISIDELAPLACIAN 
#define GUARDINCREMENTHISIDELAPLACIAN 
// Prototype for Fortran procedure INCREMENTHISIDELAPLACIAN ...
//
void FORTRAN_NAME( INCREMENTHISIDELAPLACIAN ,incrementhisidelaplacian )(
      CHFp_FRA(lapPhi)
      ,CHFp_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_INT(dir)
      ,CHFp_REAL(factor) );

#define FORT_INCREMENTHISIDELAPLACIAN FORTRAN_NAME( inlineINCREMENTHISIDELAPLACIAN, inlineINCREMENTHISIDELAPLACIAN)
#define FORTNT_INCREMENTHISIDELAPLACIAN FORTRAN_NAME( INCREMENTHISIDELAPLACIAN, incrementhisidelaplacian)

inline void FORTRAN_NAME(inlineINCREMENTHISIDELAPLACIAN, inlineINCREMENTHISIDELAPLACIAN)(
      CHFp_FRA(lapPhi)
      ,CHFp_FRA(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_INT(dir)
      ,CHFp_REAL(factor) )
{
 CH_TIMELEAF("FORT_INCREMENTHISIDELAPLACIAN");
 FORTRAN_NAME( INCREMENTHISIDELAPLACIAN ,incrementhisidelaplacian )(
      CHFt_FRA(lapPhi)
      ,CHFt_FRA(phi)
      ,CHFt_BOX(gridBox)
      ,CHFt_INT(dir)
      ,CHFt_REAL(factor) );
}
#endif  // GUARDINCREMENTHISIDELAPLACIAN 

}

#endif
