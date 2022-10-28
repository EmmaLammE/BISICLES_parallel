#ifndef _N_VECTOROPS_F_H_
#define _N_VECTOROPS_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMASKDOTPROD 
#define GUARDMASKDOTPROD 
// Prototype for Fortran procedure maskdotprod ...
//
void FORTRAN_NAME( MASKDOTPROD ,maskdotprod )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(dotProd) );

#define FORT_MASKDOTPROD FORTRAN_NAME( inlineMASKDOTPROD, inlineMASKDOTPROD)
#define FORTNT_MASKDOTPROD FORTRAN_NAME( MASKDOTPROD, maskdotprod)

inline void FORTRAN_NAME(inlineMASKDOTPROD, inlineMASKDOTPROD)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(dotProd) )
{
 CH_TIMELEAF("FORT_MASKDOTPROD");
 FORTRAN_NAME( MASKDOTPROD ,maskdotprod )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_CONST_FRA(mask)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(dotProd) );
}
#endif  // GUARDMASKDOTPROD 

#ifndef GUARDARRAYPROD 
#define GUARDARRAYPROD 
// Prototype for Fortran procedure arrayprod ...
//
void FORTRAN_NAME( ARRAYPROD ,arrayprod )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg) );

#define FORT_ARRAYPROD FORTRAN_NAME( inlineARRAYPROD, inlineARRAYPROD)
#define FORTNT_ARRAYPROD FORTRAN_NAME( ARRAYPROD, arrayprod)

inline void FORTRAN_NAME(inlineARRAYPROD, inlineARRAYPROD)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg) )
{
 CH_TIMELEAF("FORT_ARRAYPROD");
 FORTRAN_NAME( ARRAYPROD ,arrayprod )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg) );
}
#endif  // GUARDARRAYPROD 

#ifndef GUARDARRAYDIV 
#define GUARDARRAYDIV 
// Prototype for Fortran procedure arraydiv ...
//
void FORTRAN_NAME( ARRAYDIV ,arraydiv )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg) );

#define FORT_ARRAYDIV FORTRAN_NAME( inlineARRAYDIV, inlineARRAYDIV)
#define FORTNT_ARRAYDIV FORTRAN_NAME( ARRAYDIV, arraydiv)

inline void FORTRAN_NAME(inlineARRAYDIV, inlineARRAYDIV)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg) )
{
 CH_TIMELEAF("FORT_ARRAYDIV");
 FORTRAN_NAME( ARRAYDIV ,arraydiv )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg) );
}
#endif  // GUARDARRAYDIV 

#ifndef GUARDARRAYSCL 
#define GUARDARRAYSCL 
// Prototype for Fortran procedure arrayscl ...
//
void FORTRAN_NAME( ARRAYSCL ,arrayscl )(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_CONST_REAL(c) );

#define FORT_ARRAYSCL FORTRAN_NAME( inlineARRAYSCL, inlineARRAYSCL)
#define FORTNT_ARRAYSCL FORTRAN_NAME( ARRAYSCL, arrayscl)

inline void FORTRAN_NAME(inlineARRAYSCL, inlineARRAYSCL)(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_CONST_REAL(c) )
{
 CH_TIMELEAF("FORT_ARRAYSCL");
 FORTRAN_NAME( ARRAYSCL ,arrayscl )(
      CHFt_CONST_FRA(x)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg)
      ,CHFt_CONST_REAL(c) );
}
#endif  // GUARDARRAYSCL 

#ifndef GUARDARRAYABS 
#define GUARDARRAYABS 
// Prototype for Fortran procedure arrayabs ...
//
void FORTRAN_NAME( ARRAYABS ,arrayabs )(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg) );

#define FORT_ARRAYABS FORTRAN_NAME( inlineARRAYABS, inlineARRAYABS)
#define FORTNT_ARRAYABS FORTRAN_NAME( ARRAYABS, arrayabs)

inline void FORTRAN_NAME(inlineARRAYABS, inlineARRAYABS)(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg) )
{
 CH_TIMELEAF("FORT_ARRAYABS");
 FORTRAN_NAME( ARRAYABS ,arrayabs )(
      CHFt_CONST_FRA(x)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg) );
}
#endif  // GUARDARRAYABS 

#ifndef GUARDADDCONST 
#define GUARDADDCONST 
// Prototype for Fortran procedure addconst ...
//
void FORTRAN_NAME( ADDCONST ,addconst )(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_CONST_REAL(b) );

#define FORT_ADDCONST FORTRAN_NAME( inlineADDCONST, inlineADDCONST)
#define FORTNT_ADDCONST FORTRAN_NAME( ADDCONST, addconst)

inline void FORTRAN_NAME(inlineADDCONST, inlineADDCONST)(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_CONST_REAL(b) )
{
 CH_TIMELEAF("FORT_ADDCONST");
 FORTRAN_NAME( ADDCONST ,addconst )(
      CHFt_CONST_FRA(x)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg)
      ,CHFt_CONST_REAL(b) );
}
#endif  // GUARDADDCONST 

#ifndef GUARDMASKSUMXP 
#define GUARDMASKSUMXP 
// Prototype for Fortran procedure masksumxp ...
//
void FORTRAN_NAME( MASKSUMXP ,masksumxp )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(sumxp)
      ,CHFp_INT(p) );

#define FORT_MASKSUMXP FORTRAN_NAME( inlineMASKSUMXP, inlineMASKSUMXP)
#define FORTNT_MASKSUMXP FORTRAN_NAME( MASKSUMXP, masksumxp)

inline void FORTRAN_NAME(inlineMASKSUMXP, inlineMASKSUMXP)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(sumxp)
      ,CHFp_INT(p) )
{
 CH_TIMELEAF("FORT_MASKSUMXP");
 FORTRAN_NAME( MASKSUMXP ,masksumxp )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(mask)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(sumxp)
      ,CHFt_INT(p) );
}
#endif  // GUARDMASKSUMXP 

#ifndef GUARDMASKMAXNORM 
#define GUARDMASKMAXNORM 
// Prototype for Fortran procedure maskmaxnorm ...
//
void FORTRAN_NAME( MASKMAXNORM ,maskmaxnorm )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(m) );

#define FORT_MASKMAXNORM FORTRAN_NAME( inlineMASKMAXNORM, inlineMASKMAXNORM)
#define FORTNT_MASKMAXNORM FORTRAN_NAME( MASKMAXNORM, maskmaxnorm)

inline void FORTRAN_NAME(inlineMASKMAXNORM, inlineMASKMAXNORM)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(m) )
{
 CH_TIMELEAF("FORT_MASKMAXNORM");
 FORTRAN_NAME( MASKMAXNORM ,maskmaxnorm )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(mask)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(m) );
}
#endif  // GUARDMASKMAXNORM 

#ifndef GUARDMASKMIN 
#define GUARDMASKMIN 
// Prototype for Fortran procedure maskmin ...
//
void FORTRAN_NAME( MASKMIN ,maskmin )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(m) );

#define FORT_MASKMIN FORTRAN_NAME( inlineMASKMIN, inlineMASKMIN)
#define FORTNT_MASKMIN FORTRAN_NAME( MASKMIN, maskmin)

inline void FORTRAN_NAME(inlineMASKMIN, inlineMASKMIN)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(m) )
{
 CH_TIMELEAF("FORT_MASKMIN");
 FORTRAN_NAME( MASKMIN ,maskmin )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(mask)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(m) );
}
#endif  // GUARDMASKMIN 

#ifndef GUARDMASKWTDSQ 
#define GUARDMASKWTDSQ 
// Prototype for Fortran procedure maskwtdsq ...
//
void FORTRAN_NAME( MASKWTDSQ ,maskwtdsq )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(norm) );

#define FORT_MASKWTDSQ FORTRAN_NAME( inlineMASKWTDSQ, inlineMASKWTDSQ)
#define FORTNT_MASKWTDSQ FORTRAN_NAME( MASKWTDSQ, maskwtdsq)

inline void FORTRAN_NAME(inlineMASKWTDSQ, inlineMASKWTDSQ)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(mask)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(norm) )
{
 CH_TIMELEAF("FORT_MASKWTDSQ");
 FORTRAN_NAME( MASKWTDSQ ,maskwtdsq )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_FRA(mask)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(norm) );
}
#endif  // GUARDMASKWTDSQ 

#ifndef GUARDWTDSIGNSQ 
#define GUARDWTDSIGNSQ 
// Prototype for Fortran procedure wtdsignsq ...
//
void FORTRAN_NAME( WTDSIGNSQ ,wtdsignsq )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(id)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(norm) );

#define FORT_WTDSIGNSQ FORTRAN_NAME( inlineWTDSIGNSQ, inlineWTDSIGNSQ)
#define FORTNT_WTDSIGNSQ FORTRAN_NAME( WTDSIGNSQ, wtdsignsq)

inline void FORTRAN_NAME(inlineWTDSIGNSQ, inlineWTDSIGNSQ)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(id)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(norm) )
{
 CH_TIMELEAF("FORT_WTDSIGNSQ");
 FORTRAN_NAME( WTDSIGNSQ ,wtdsignsq )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_FRA(id)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(norm) );
}
#endif  // GUARDWTDSIGNSQ 

#ifndef GUARDARRAYCOMP 
#define GUARDARRAYCOMP 
// Prototype for Fortran procedure arraycomp ...
//
void FORTRAN_NAME( ARRAYCOMP ,arraycomp )(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_CONST_REAL(c) );

#define FORT_ARRAYCOMP FORTRAN_NAME( inlineARRAYCOMP, inlineARRAYCOMP)
#define FORTNT_ARRAYCOMP FORTRAN_NAME( ARRAYCOMP, arraycomp)

inline void FORTRAN_NAME(inlineARRAYCOMP, inlineARRAYCOMP)(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_CONST_REAL(c) )
{
 CH_TIMELEAF("FORT_ARRAYCOMP");
 FORTRAN_NAME( ARRAYCOMP ,arraycomp )(
      CHFt_CONST_FRA(x)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg)
      ,CHFt_CONST_REAL(c) );
}
#endif  // GUARDARRAYCOMP 

#ifndef GUARDINVWCHK 
#define GUARDINVWCHK 
// Prototype for Fortran procedure invwchk ...
//
void FORTRAN_NAME( INVWCHK ,invwchk )(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_INT(nonzero) );

#define FORT_INVWCHK FORTRAN_NAME( inlineINVWCHK, inlineINVWCHK)
#define FORTNT_INVWCHK FORTRAN_NAME( INVWCHK, invwchk)

inline void FORTRAN_NAME(inlineINVWCHK, inlineINVWCHK)(
      CHFp_CONST_FRA(x)
      ,CHFp_FRA(z)
      ,CHFp_BOX(reg)
      ,CHFp_INT(nonzero) )
{
 CH_TIMELEAF("FORT_INVWCHK");
 FORTRAN_NAME( INVWCHK ,invwchk )(
      CHFt_CONST_FRA(x)
      ,CHFt_FRA(z)
      ,CHFt_BOX(reg)
      ,CHFt_INT(nonzero) );
}
#endif  // GUARDINVWCHK 

#ifndef GUARDCONSTRCHK 
#define GUARDCONSTRCHK 
// Prototype for Fortran procedure constrchk ...
//
void FORTRAN_NAME( CONSTRCHK ,constrchk )(
      CHFp_CONST_FRA(c)
      ,CHFp_CONST_FRA(x)
      ,CHFp_FRA(m)
      ,CHFp_BOX(reg)
      ,CHFp_INT(allpassed) );

#define FORT_CONSTRCHK FORTRAN_NAME( inlineCONSTRCHK, inlineCONSTRCHK)
#define FORTNT_CONSTRCHK FORTRAN_NAME( CONSTRCHK, constrchk)

inline void FORTRAN_NAME(inlineCONSTRCHK, inlineCONSTRCHK)(
      CHFp_CONST_FRA(c)
      ,CHFp_CONST_FRA(x)
      ,CHFp_FRA(m)
      ,CHFp_BOX(reg)
      ,CHFp_INT(allpassed) )
{
 CH_TIMELEAF("FORT_CONSTRCHK");
 FORTRAN_NAME( CONSTRCHK ,constrchk )(
      CHFt_CONST_FRA(c)
      ,CHFt_CONST_FRA(x)
      ,CHFt_FRA(m)
      ,CHFt_BOX(reg)
      ,CHFt_INT(allpassed) );
}
#endif  // GUARDCONSTRCHK 

#ifndef GUARDMINQUOT 
#define GUARDMINQUOT 
// Prototype for Fortran procedure minquot ...
//
void FORTRAN_NAME( MINQUOT ,minquot )(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(q)
      ,CHFp_INT(nmin) );

#define FORT_MINQUOT FORTRAN_NAME( inlineMINQUOT, inlineMINQUOT)
#define FORTNT_MINQUOT FORTRAN_NAME( MINQUOT, minquot)

inline void FORTRAN_NAME(inlineMINQUOT, inlineMINQUOT)(
      CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_BOX(reg)
      ,CHFp_REAL(q)
      ,CHFp_INT(nmin) )
{
 CH_TIMELEAF("FORT_MINQUOT");
 FORTRAN_NAME( MINQUOT ,minquot )(
      CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_BOX(reg)
      ,CHFt_REAL(q)
      ,CHFt_INT(nmin) );
}
#endif  // GUARDMINQUOT 

}

#endif
