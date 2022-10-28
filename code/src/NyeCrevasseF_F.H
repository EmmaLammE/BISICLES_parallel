#ifndef _NYECREVASSEF_F_H_
#define _NYECREVASSEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUPLUSUT 
#define GUARDUPLUSUT 
// Prototype for Fortran procedure uplusut ...
//
void FORTRAN_NAME( UPLUSUT ,uplusut )(
      CHFp_FRA(e)
      ,CHFp_CONST_FRA(u)
      ,CHFp_INT(xxComp)
      ,CHFp_INT(xyComp)
      ,CHFp_INT(yxComp)
      ,CHFp_INT(yyComp)
      ,CHFp_BOX(box) );

#define FORT_UPLUSUT FORTRAN_NAME( inlineUPLUSUT, inlineUPLUSUT)
#define FORTNT_UPLUSUT FORTRAN_NAME( UPLUSUT, uplusut)

inline void FORTRAN_NAME(inlineUPLUSUT, inlineUPLUSUT)(
      CHFp_FRA(e)
      ,CHFp_CONST_FRA(u)
      ,CHFp_INT(xxComp)
      ,CHFp_INT(xyComp)
      ,CHFp_INT(yxComp)
      ,CHFp_INT(yyComp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPLUSUT");
 FORTRAN_NAME( UPLUSUT ,uplusut )(
      CHFt_FRA(e)
      ,CHFt_CONST_FRA(u)
      ,CHFt_INT(xxComp)
      ,CHFt_INT(xyComp)
      ,CHFt_INT(yxComp)
      ,CHFt_INT(yyComp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPLUSUT 

#ifndef GUARDEPLUSTRE 
#define GUARDEPLUSTRE 
// Prototype for Fortran procedure eplustre ...
//
void FORTRAN_NAME( EPLUSTRE ,eplustre )(
      CHFp_FRA(e)
      ,CHFp_CONST_FRA(u)
      ,CHFp_INT(xxComp)
      ,CHFp_INT(xyComp)
      ,CHFp_INT(yxComp)
      ,CHFp_INT(yyComp)
      ,CHFp_BOX(box) );

#define FORT_EPLUSTRE FORTRAN_NAME( inlineEPLUSTRE, inlineEPLUSTRE)
#define FORTNT_EPLUSTRE FORTRAN_NAME( EPLUSTRE, eplustre)

inline void FORTRAN_NAME(inlineEPLUSTRE, inlineEPLUSTRE)(
      CHFp_FRA(e)
      ,CHFp_CONST_FRA(u)
      ,CHFp_INT(xxComp)
      ,CHFp_INT(xyComp)
      ,CHFp_INT(yxComp)
      ,CHFp_INT(yyComp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EPLUSTRE");
 FORTRAN_NAME( EPLUSTRE ,eplustre )(
      CHFt_FRA(e)
      ,CHFt_CONST_FRA(u)
      ,CHFt_INT(xxComp)
      ,CHFt_INT(xyComp)
      ,CHFt_INT(yxComp)
      ,CHFt_INT(yyComp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEPLUSTRE 

#ifndef GUARDSYMTEIGEN 
#define GUARDSYMTEIGEN 
// Prototype for Fortran procedure symteigen ...
//
void FORTRAN_NAME( SYMTEIGEN ,symteigen )(
      CHFp_FRA(lambda)
      ,CHFp_CONST_FRA(T)
      ,CHFp_INT(xxComp)
      ,CHFp_INT(xyComp)
      ,CHFp_INT(yxComp)
      ,CHFp_INT(yyComp)
      ,CHFp_BOX(box) );

#define FORT_SYMTEIGEN FORTRAN_NAME( inlineSYMTEIGEN, inlineSYMTEIGEN)
#define FORTNT_SYMTEIGEN FORTRAN_NAME( SYMTEIGEN, symteigen)

inline void FORTRAN_NAME(inlineSYMTEIGEN, inlineSYMTEIGEN)(
      CHFp_FRA(lambda)
      ,CHFp_CONST_FRA(T)
      ,CHFp_INT(xxComp)
      ,CHFp_INT(xyComp)
      ,CHFp_INT(yxComp)
      ,CHFp_INT(yyComp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SYMTEIGEN");
 FORTRAN_NAME( SYMTEIGEN ,symteigen )(
      CHFt_FRA(lambda)
      ,CHFt_CONST_FRA(T)
      ,CHFt_INT(xxComp)
      ,CHFt_INT(xyComp)
      ,CHFt_INT(yxComp)
      ,CHFt_INT(yyComp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSYMTEIGEN 

#ifndef GUARDNYECREVASSEDEPTHTP 
#define GUARDNYECREVASSEDEPTHTP 
// Prototype for Fortran procedure nyecrevassedepthtp ...
//
void FORTRAN_NAME( NYECREVASSEDEPTHTP ,nyecrevassedepthtp )(
      CHFp_FRA1(dfab)
      ,CHFp_CONST_FRA1(dwfab)
      ,CHFp_CONST_FRA1(hfab)
      ,CHFp_CONST_FRA1(tpfab)
      ,CHFp_CONST_FRA1(habfab)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(g)
      ,CHFp_CONST_REAL(eps)
      ,CHFp_CONST_REAL(dmax)
      ,CHFp_BOX(box) );

#define FORT_NYECREVASSEDEPTHTP FORTRAN_NAME( inlineNYECREVASSEDEPTHTP, inlineNYECREVASSEDEPTHTP)
#define FORTNT_NYECREVASSEDEPTHTP FORTRAN_NAME( NYECREVASSEDEPTHTP, nyecrevassedepthtp)

inline void FORTRAN_NAME(inlineNYECREVASSEDEPTHTP, inlineNYECREVASSEDEPTHTP)(
      CHFp_FRA1(dfab)
      ,CHFp_CONST_FRA1(dwfab)
      ,CHFp_CONST_FRA1(hfab)
      ,CHFp_CONST_FRA1(tpfab)
      ,CHFp_CONST_FRA1(habfab)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(g)
      ,CHFp_CONST_REAL(eps)
      ,CHFp_CONST_REAL(dmax)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_NYECREVASSEDEPTHTP");
 FORTRAN_NAME( NYECREVASSEDEPTHTP ,nyecrevassedepthtp )(
      CHFt_FRA1(dfab)
      ,CHFt_CONST_FRA1(dwfab)
      ,CHFt_CONST_FRA1(hfab)
      ,CHFt_CONST_FRA1(tpfab)
      ,CHFt_CONST_FRA1(habfab)
      ,CHFt_CONST_REAL(rhoi)
      ,CHFt_CONST_REAL(rhoo)
      ,CHFt_CONST_REAL(g)
      ,CHFt_CONST_REAL(eps)
      ,CHFt_CONST_REAL(dmax)
      ,CHFt_BOX(box) );
}
#endif  // GUARDNYECREVASSEDEPTHTP 

#ifndef GUARDNYECREVASSEDEPTHT 
#define GUARDNYECREVASSEDEPTHT 
// Prototype for Fortran procedure nyecrevassedeptht ...
//
void FORTRAN_NAME( NYECREVASSEDEPTHT ,nyecrevassedeptht )(
      CHFp_FRA1(dfab)
      ,CHFp_CONST_FRA1(dwfab)
      ,CHFp_CONST_FRA1(hfab)
      ,CHFp_CONST_FRA1(tfab)
      ,CHFp_CONST_FRA1(habfab)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(g)
      ,CHFp_CONST_REAL(eps)
      ,CHFp_CONST_REAL(dlim)
      ,CHFp_BOX(box) );

#define FORT_NYECREVASSEDEPTHT FORTRAN_NAME( inlineNYECREVASSEDEPTHT, inlineNYECREVASSEDEPTHT)
#define FORTNT_NYECREVASSEDEPTHT FORTRAN_NAME( NYECREVASSEDEPTHT, nyecrevassedeptht)

inline void FORTRAN_NAME(inlineNYECREVASSEDEPTHT, inlineNYECREVASSEDEPTHT)(
      CHFp_FRA1(dfab)
      ,CHFp_CONST_FRA1(dwfab)
      ,CHFp_CONST_FRA1(hfab)
      ,CHFp_CONST_FRA1(tfab)
      ,CHFp_CONST_FRA1(habfab)
      ,CHFp_CONST_REAL(rhoi)
      ,CHFp_CONST_REAL(rhoo)
      ,CHFp_CONST_REAL(g)
      ,CHFp_CONST_REAL(eps)
      ,CHFp_CONST_REAL(dlim)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_NYECREVASSEDEPTHT");
 FORTRAN_NAME( NYECREVASSEDEPTHT ,nyecrevassedeptht )(
      CHFt_FRA1(dfab)
      ,CHFt_CONST_FRA1(dwfab)
      ,CHFt_CONST_FRA1(hfab)
      ,CHFt_CONST_FRA1(tfab)
      ,CHFt_CONST_FRA1(habfab)
      ,CHFt_CONST_REAL(rhoi)
      ,CHFt_CONST_REAL(rhoo)
      ,CHFt_CONST_REAL(g)
      ,CHFt_CONST_REAL(eps)
      ,CHFt_CONST_REAL(dlim)
      ,CHFt_BOX(box) );
}
#endif  // GUARDNYECREVASSEDEPTHT 

#ifndef GUARDFABMAX 
#define GUARDFABMAX 
// Prototype for Fortran procedure fabmax ...
//
void FORTRAN_NAME( FABMAX ,fabmax )(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_BOX(box) );

#define FORT_FABMAX FORTRAN_NAME( inlineFABMAX, inlineFABMAX)
#define FORTNT_FABMAX FORTRAN_NAME( FABMAX, fabmax)

inline void FORTRAN_NAME(inlineFABMAX, inlineFABMAX)(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FABMAX");
 FORTRAN_NAME( FABMAX ,fabmax )(
      CHFt_FRA1(a)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFABMAX 

#ifndef GUARDFABMIN 
#define GUARDFABMIN 
// Prototype for Fortran procedure fabmin ...
//
void FORTRAN_NAME( FABMIN ,fabmin )(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_BOX(box) );

#define FORT_FABMIN FORTRAN_NAME( inlineFABMIN, inlineFABMIN)
#define FORTNT_FABMIN FORTRAN_NAME( FABMIN, fabmin)

inline void FORTRAN_NAME(inlineFABMIN, inlineFABMIN)(
      CHFp_FRA1(a)
      ,CHFp_CONST_FRA1(b)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FABMIN");
 FORTRAN_NAME( FABMIN ,fabmin )(
      CHFt_FRA1(a)
      ,CHFt_CONST_FRA1(b)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFABMIN 

#ifndef GUARDCREVASSEMU 
#define GUARDCREVASSEMU 
// Prototype for Fortran procedure crevassemu ...
//
void FORTRAN_NAME( CREVASSEMU ,crevassemu )(
      CHFp_FRA1(mu)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REAL(eps)
      ,CHFp_BOX(box) );

#define FORT_CREVASSEMU FORTRAN_NAME( inlineCREVASSEMU, inlineCREVASSEMU)
#define FORTNT_CREVASSEMU FORTRAN_NAME( CREVASSEMU, crevassemu)

inline void FORTRAN_NAME(inlineCREVASSEMU, inlineCREVASSEMU)(
      CHFp_FRA1(mu)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REAL(eps)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CREVASSEMU");
 FORTRAN_NAME( CREVASSEMU ,crevassemu )(
      CHFt_FRA1(mu)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_FRA1(depth)
      ,CHFt_CONST_REAL(eps)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCREVASSEMU 

}

#endif
