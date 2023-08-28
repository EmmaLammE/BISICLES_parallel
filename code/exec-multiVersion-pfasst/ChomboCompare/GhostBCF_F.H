#ifndef _GHOSTBCF_F_H_
#define _GHOSTBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDBOXGHOSTBC 
#define GUARDBOXGHOSTBC 
// Prototype for Fortran procedure BOXGHOSTBC ...
//
void FORTRAN_NAME( BOXGHOSTBC ,boxghostbc )(
      CHFp_FRA(state)
      ,CHFp_CONST_FRA(neumfac)
      ,CHFp_CONST_FRA(dircfac)
      ,CHFp_CONST_FRA(inhmval)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_BOXGHOSTBC FORTRAN_NAME( inlineBOXGHOSTBC, inlineBOXGHOSTBC)
#define FORTNT_BOXGHOSTBC FORTRAN_NAME( BOXGHOSTBC, boxghostbc)

inline void FORTRAN_NAME(inlineBOXGHOSTBC, inlineBOXGHOSTBC)(
      CHFp_FRA(state)
      ,CHFp_CONST_FRA(neumfac)
      ,CHFp_CONST_FRA(dircfac)
      ,CHFp_CONST_FRA(inhmval)
      ,CHFp_BOX(facebox)
      ,CHFp_BOX(interiorBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_BOXGHOSTBC");
 FORTRAN_NAME( BOXGHOSTBC ,boxghostbc )(
      CHFt_FRA(state)
      ,CHFt_CONST_FRA(neumfac)
      ,CHFt_CONST_FRA(dircfac)
      ,CHFt_CONST_FRA(inhmval)
      ,CHFt_BOX(facebox)
      ,CHFt_BOX(interiorBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDBOXGHOSTBC 

}

#endif
