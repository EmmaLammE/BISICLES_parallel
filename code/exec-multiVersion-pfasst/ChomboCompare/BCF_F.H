#ifndef _BCF_F_H_
#define _BCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEXTRAPGHOSTBC 
#define GUARDEXTRAPGHOSTBC 
// Prototype for Fortran procedure EXTRAPGHOSTBC ...
//
void FORTRAN_NAME( EXTRAPGHOSTBC ,extrapghostbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_EXTRAPGHOSTBC FORTRAN_NAME( inlineEXTRAPGHOSTBC, inlineEXTRAPGHOSTBC)
#define FORTNT_EXTRAPGHOSTBC FORTRAN_NAME( EXTRAPGHOSTBC, extrapghostbc)

inline void FORTRAN_NAME(inlineEXTRAPGHOSTBC, inlineEXTRAPGHOSTBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_EXTRAPGHOSTBC");
 FORTRAN_NAME( EXTRAPGHOSTBC ,extrapghostbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(bcBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDEXTRAPGHOSTBC 

#ifndef GUARDHOEXTRAPGHOSTBC 
#define GUARDHOEXTRAPGHOSTBC 
// Prototype for Fortran procedure HOEXTRAPGHOSTBC ...
//
void FORTRAN_NAME( HOEXTRAPGHOSTBC ,hoextrapghostbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_HOEXTRAPGHOSTBC FORTRAN_NAME( inlineHOEXTRAPGHOSTBC, inlineHOEXTRAPGHOSTBC)
#define FORTNT_HOEXTRAPGHOSTBC FORTRAN_NAME( HOEXTRAPGHOSTBC, hoextrapghostbc)

inline void FORTRAN_NAME(inlineHOEXTRAPGHOSTBC, inlineHOEXTRAPGHOSTBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_HOEXTRAPGHOSTBC");
 FORTRAN_NAME( HOEXTRAPGHOSTBC ,hoextrapghostbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(bcBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDHOEXTRAPGHOSTBC 

#ifndef GUARDREFLECTGHOSTBC 
#define GUARDREFLECTGHOSTBC 
// Prototype for Fortran procedure REFLECTGHOSTBC ...
//
void FORTRAN_NAME( REFLECTGHOSTBC ,reflectghostbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp)
      ,CHFp_CONST_REAL(scale) );

#define FORT_REFLECTGHOSTBC FORTRAN_NAME( inlineREFLECTGHOSTBC, inlineREFLECTGHOSTBC)
#define FORTNT_REFLECTGHOSTBC FORTRAN_NAME( REFLECTGHOSTBC, reflectghostbc)

inline void FORTRAN_NAME(inlineREFLECTGHOSTBC, inlineREFLECTGHOSTBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_REFLECTGHOSTBC");
 FORTRAN_NAME( REFLECTGHOSTBC ,reflectghostbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(bcBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDREFLECTGHOSTBC 

}

#endif
