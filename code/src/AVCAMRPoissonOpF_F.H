#ifndef _AVCAMRPOISSONOPF_F_H_
#define _AVCAMRPOISSONOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGSRBHELMHOLTZAVC1D 
#define GUARDGSRBHELMHOLTZAVC1D 
// Prototype for Fortran procedure GSRBHELMHOLTZAVC1D ...
//
void FORTRAN_NAME( GSRBHELMHOLTZAVC1D ,gsrbhelmholtzavc1d )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_CONST_INT(icol) );

#define FORT_GSRBHELMHOLTZAVC1D FORTRAN_NAME( inlineGSRBHELMHOLTZAVC1D, inlineGSRBHELMHOLTZAVC1D)
#define FORTNT_GSRBHELMHOLTZAVC1D FORTRAN_NAME( GSRBHELMHOLTZAVC1D, gsrbhelmholtzavc1d)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZAVC1D, inlineGSRBHELMHOLTZAVC1D)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_CONST_INT(icol) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZAVC1D");
 FORTRAN_NAME( GSRBHELMHOLTZAVC1D ,gsrbhelmholtzavc1d )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA1(lambda)
      ,CHFt_CONST_INT(icol) );
}
#endif  // GUARDGSRBHELMHOLTZAVC1D 

#ifndef GUARDGSRBHELMHOLTZAVC2D 
#define GUARDGSRBHELMHOLTZAVC2D 
// Prototype for Fortran procedure GSRBHELMHOLTZAVC2D ...
//
void FORTRAN_NAME( GSRBHELMHOLTZAVC2D ,gsrbhelmholtzavc2d )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_CONST_INT(icol)
      ,CHFp_CONST_INT(jcol) );

#define FORT_GSRBHELMHOLTZAVC2D FORTRAN_NAME( inlineGSRBHELMHOLTZAVC2D, inlineGSRBHELMHOLTZAVC2D)
#define FORTNT_GSRBHELMHOLTZAVC2D FORTRAN_NAME( GSRBHELMHOLTZAVC2D, gsrbhelmholtzavc2d)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZAVC2D, inlineGSRBHELMHOLTZAVC2D)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_CONST_INT(icol)
      ,CHFp_CONST_INT(jcol) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZAVC2D");
 FORTRAN_NAME( GSRBHELMHOLTZAVC2D ,gsrbhelmholtzavc2d )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA1(lambda)
      ,CHFt_CONST_INT(icol)
      ,CHFt_CONST_INT(jcol) );
}
#endif  // GUARDGSRBHELMHOLTZAVC2D 

#ifndef GUARDGSRBHELMHOLTZAVC3D 
#define GUARDGSRBHELMHOLTZAVC3D 
// Prototype for Fortran procedure GSRBHELMHOLTZAVC3D ...
//
void FORTRAN_NAME( GSRBHELMHOLTZAVC3D ,gsrbhelmholtzavc3d )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_CONST_INT(icol)
      ,CHFp_CONST_INT(jcol)
      ,CHFp_CONST_INT(kcol) );

#define FORT_GSRBHELMHOLTZAVC3D FORTRAN_NAME( inlineGSRBHELMHOLTZAVC3D, inlineGSRBHELMHOLTZAVC3D)
#define FORTNT_GSRBHELMHOLTZAVC3D FORTRAN_NAME( GSRBHELMHOLTZAVC3D, gsrbhelmholtzavc3d)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZAVC3D, inlineGSRBHELMHOLTZAVC3D)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_CONST_INT(icol)
      ,CHFp_CONST_INT(jcol)
      ,CHFp_CONST_INT(kcol) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZAVC3D");
 FORTRAN_NAME( GSRBHELMHOLTZAVC3D ,gsrbhelmholtzavc3d )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_CONST_FRA1(lambda)
      ,CHFt_CONST_INT(icol)
      ,CHFt_CONST_INT(jcol)
      ,CHFt_CONST_INT(kcol) );
}
#endif  // GUARDGSRBHELMHOLTZAVC3D 

#ifndef GUARDAVCCOMPUTEOP1D 
#define GUARDAVCCOMPUTEOP1D 
// Prototype for Fortran procedure AVCCOMPUTEOP1D ...
//
void FORTRAN_NAME( AVCCOMPUTEOP1D ,avccomputeop1d )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_AVCCOMPUTEOP1D FORTRAN_NAME( inlineAVCCOMPUTEOP1D, inlineAVCCOMPUTEOP1D)
#define FORTNT_AVCCOMPUTEOP1D FORTRAN_NAME( AVCCOMPUTEOP1D, avccomputeop1d)

inline void FORTRAN_NAME(inlineAVCCOMPUTEOP1D, inlineAVCCOMPUTEOP1D)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_AVCCOMPUTEOP1D");
 FORTRAN_NAME( AVCCOMPUTEOP1D ,avccomputeop1d )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDAVCCOMPUTEOP1D 

#ifndef GUARDAVCCOMPUTEOP2D 
#define GUARDAVCCOMPUTEOP2D 
// Prototype for Fortran procedure AVCCOMPUTEOP2D ...
//
void FORTRAN_NAME( AVCCOMPUTEOP2D ,avccomputeop2d )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_AVCCOMPUTEOP2D FORTRAN_NAME( inlineAVCCOMPUTEOP2D, inlineAVCCOMPUTEOP2D)
#define FORTNT_AVCCOMPUTEOP2D FORTRAN_NAME( AVCCOMPUTEOP2D, avccomputeop2d)

inline void FORTRAN_NAME(inlineAVCCOMPUTEOP2D, inlineAVCCOMPUTEOP2D)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_AVCCOMPUTEOP2D");
 FORTRAN_NAME( AVCCOMPUTEOP2D ,avccomputeop2d )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDAVCCOMPUTEOP2D 

#ifndef GUARDAVCCOMPUTEOP3D 
#define GUARDAVCCOMPUTEOP3D 
// Prototype for Fortran procedure AVCCOMPUTEOP3D ...
//
void FORTRAN_NAME( AVCCOMPUTEOP3D ,avccomputeop3d )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_AVCCOMPUTEOP3D FORTRAN_NAME( inlineAVCCOMPUTEOP3D, inlineAVCCOMPUTEOP3D)
#define FORTNT_AVCCOMPUTEOP3D FORTRAN_NAME( AVCCOMPUTEOP3D, avccomputeop3d)

inline void FORTRAN_NAME(inlineAVCCOMPUTEOP3D, inlineAVCCOMPUTEOP3D)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_AVCCOMPUTEOP3D");
 FORTRAN_NAME( AVCCOMPUTEOP3D ,avccomputeop3d )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDAVCCOMPUTEOP3D 

#ifndef GUARDAVCCOMPUTERES1D 
#define GUARDAVCCOMPUTERES1D 
// Prototype for Fortran procedure AVCCOMPUTERES1D ...
//
void FORTRAN_NAME( AVCCOMPUTERES1D ,avccomputeres1d )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_AVCCOMPUTERES1D FORTRAN_NAME( inlineAVCCOMPUTERES1D, inlineAVCCOMPUTERES1D)
#define FORTNT_AVCCOMPUTERES1D FORTRAN_NAME( AVCCOMPUTERES1D, avccomputeres1d)

inline void FORTRAN_NAME(inlineAVCCOMPUTERES1D, inlineAVCCOMPUTERES1D)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_AVCCOMPUTERES1D");
 FORTRAN_NAME( AVCCOMPUTERES1D ,avccomputeres1d )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDAVCCOMPUTERES1D 

#ifndef GUARDAVCCOMPUTERES2D 
#define GUARDAVCCOMPUTERES2D 
// Prototype for Fortran procedure AVCCOMPUTERES2D ...
//
void FORTRAN_NAME( AVCCOMPUTERES2D ,avccomputeres2d )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_AVCCOMPUTERES2D FORTRAN_NAME( inlineAVCCOMPUTERES2D, inlineAVCCOMPUTERES2D)
#define FORTNT_AVCCOMPUTERES2D FORTRAN_NAME( AVCCOMPUTERES2D, avccomputeres2d)

inline void FORTRAN_NAME(inlineAVCCOMPUTERES2D, inlineAVCCOMPUTERES2D)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_AVCCOMPUTERES2D");
 FORTRAN_NAME( AVCCOMPUTERES2D ,avccomputeres2d )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDAVCCOMPUTERES2D 

#ifndef GUARDAVCCOMPUTERES3D 
#define GUARDAVCCOMPUTERES3D 
// Prototype for Fortran procedure AVCCOMPUTERES3D ...
//
void FORTRAN_NAME( AVCCOMPUTERES3D ,avccomputeres3d )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_AVCCOMPUTERES3D FORTRAN_NAME( inlineAVCCOMPUTERES3D, inlineAVCCOMPUTERES3D)
#define FORTNT_AVCCOMPUTERES3D FORTRAN_NAME( AVCCOMPUTERES3D, avccomputeres3d)

inline void FORTRAN_NAME(inlineAVCCOMPUTERES3D, inlineAVCCOMPUTERES3D)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_AVCCOMPUTERES3D");
 FORTRAN_NAME( AVCCOMPUTERES3D ,avccomputeres3d )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDAVCCOMPUTERES3D 

#ifndef GUARDRESTRICTRESAVC1D 
#define GUARDRESTRICTRESAVC1D 
// Prototype for Fortran procedure RESTRICTRESAVC1D ...
//
void FORTRAN_NAME( RESTRICTRESAVC1D ,restrictresavc1d )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRESAVC1D FORTRAN_NAME( inlineRESTRICTRESAVC1D, inlineRESTRICTRESAVC1D)
#define FORTNT_RESTRICTRESAVC1D FORTRAN_NAME( RESTRICTRESAVC1D, restrictresavc1d)

inline void FORTRAN_NAME(inlineRESTRICTRESAVC1D, inlineRESTRICTRESAVC1D)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRESAVC1D");
 FORTRAN_NAME( RESTRICTRESAVC1D ,restrictresavc1d )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRESAVC1D 

#ifndef GUARDRESTRICTRESAVC2D 
#define GUARDRESTRICTRESAVC2D 
// Prototype for Fortran procedure RESTRICTRESAVC2D ...
//
void FORTRAN_NAME( RESTRICTRESAVC2D ,restrictresavc2d )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRESAVC2D FORTRAN_NAME( inlineRESTRICTRESAVC2D, inlineRESTRICTRESAVC2D)
#define FORTNT_RESTRICTRESAVC2D FORTRAN_NAME( RESTRICTRESAVC2D, restrictresavc2d)

inline void FORTRAN_NAME(inlineRESTRICTRESAVC2D, inlineRESTRICTRESAVC2D)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRESAVC2D");
 FORTRAN_NAME( RESTRICTRESAVC2D ,restrictresavc2d )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRESAVC2D 

#ifndef GUARDRESTRICTRESAVC3D 
#define GUARDRESTRICTRESAVC3D 
// Prototype for Fortran procedure RESTRICTRESAVC3D ...
//
void FORTRAN_NAME( RESTRICTRESAVC3D ,restrictresavc3d )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRESAVC3D FORTRAN_NAME( inlineRESTRICTRESAVC3D, inlineRESTRICTRESAVC3D)
#define FORTNT_RESTRICTRESAVC3D FORTRAN_NAME( RESTRICTRESAVC3D, restrictresavc3d)

inline void FORTRAN_NAME(inlineRESTRICTRESAVC3D, inlineRESTRICTRESAVC3D)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA1(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRESAVC3D");
 FORTRAN_NAME( RESTRICTRESAVC3D ,restrictresavc3d )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA1(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRESAVC3D 

#ifndef GUARDSUMAFACES 
#define GUARDSUMAFACES 
// Prototype for Fortran procedure SUMAFACES ...
//
void FORTRAN_NAME( SUMAFACES ,sumafaces )(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoefs)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(scale) );

#define FORT_SUMAFACES FORTRAN_NAME( inlineSUMAFACES, inlineSUMAFACES)
#define FORTNT_SUMAFACES FORTRAN_NAME( SUMAFACES, sumafaces)

inline void FORTRAN_NAME(inlineSUMAFACES, inlineSUMAFACES)(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoefs)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_SUMAFACES");
 FORTRAN_NAME( SUMAFACES ,sumafaces )(
      CHFt_FRA1(lhs)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoefs)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDSUMAFACES 

}

#endif
