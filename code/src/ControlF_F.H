#ifndef _CONTROLF_F_H_
#define _CONTROLF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADJRHSSPEEDCTRL 
#define GUARDADJRHSSPEEDCTRL 
// Prototype for Fortran procedure adjrhsspeedctrl ...
//
void FORTRAN_NAME( ADJRHSSPEEDCTRL ,adjrhsspeedctrl )(
      CHFp_FRA1(rhsx)
      ,CHFp_FRA1(rhsy)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(umx)
      ,CHFp_CONST_FRA1(umy)
      ,CHFp_CONST_FRA1(uox)
      ,CHFp_CONST_FRA1(uoy)
      ,CHFp_BOX(box) );

#define FORT_ADJRHSSPEEDCTRL FORTRAN_NAME( inlineADJRHSSPEEDCTRL, inlineADJRHSSPEEDCTRL)
#define FORTNT_ADJRHSSPEEDCTRL FORTRAN_NAME( ADJRHSSPEEDCTRL, adjrhsspeedctrl)

inline void FORTRAN_NAME(inlineADJRHSSPEEDCTRL, inlineADJRHSSPEEDCTRL)(
      CHFp_FRA1(rhsx)
      ,CHFp_FRA1(rhsy)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(umx)
      ,CHFp_CONST_FRA1(umy)
      ,CHFp_CONST_FRA1(uox)
      ,CHFp_CONST_FRA1(uoy)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADJRHSSPEEDCTRL");
 FORTRAN_NAME( ADJRHSSPEEDCTRL ,adjrhsspeedctrl )(
      CHFt_FRA1(rhsx)
      ,CHFt_FRA1(rhsy)
      ,CHFt_FRA1(misfit)
      ,CHFt_CONST_FRA1(umx)
      ,CHFt_CONST_FRA1(umy)
      ,CHFt_CONST_FRA1(uox)
      ,CHFt_CONST_FRA1(uoy)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADJRHSSPEEDCTRL 

#ifndef GUARDADJRHSVELCTRL 
#define GUARDADJRHSVELCTRL 
// Prototype for Fortran procedure adjrhsvelctrl ...
//
void FORTRAN_NAME( ADJRHSVELCTRL ,adjrhsvelctrl )(
      CHFp_FRA1(rhsx)
      ,CHFp_FRA1(rhsy)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(umx)
      ,CHFp_CONST_FRA1(umy)
      ,CHFp_CONST_FRA1(uox)
      ,CHFp_CONST_FRA1(uoy)
      ,CHFp_BOX(box) );

#define FORT_ADJRHSVELCTRL FORTRAN_NAME( inlineADJRHSVELCTRL, inlineADJRHSVELCTRL)
#define FORTNT_ADJRHSVELCTRL FORTRAN_NAME( ADJRHSVELCTRL, adjrhsvelctrl)

inline void FORTRAN_NAME(inlineADJRHSVELCTRL, inlineADJRHSVELCTRL)(
      CHFp_FRA1(rhsx)
      ,CHFp_FRA1(rhsy)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(umx)
      ,CHFp_CONST_FRA1(umy)
      ,CHFp_CONST_FRA1(uox)
      ,CHFp_CONST_FRA1(uoy)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADJRHSVELCTRL");
 FORTRAN_NAME( ADJRHSVELCTRL ,adjrhsvelctrl )(
      CHFt_FRA1(rhsx)
      ,CHFt_FRA1(rhsy)
      ,CHFt_FRA1(misfit)
      ,CHFt_CONST_FRA1(umx)
      ,CHFt_CONST_FRA1(umy)
      ,CHFt_CONST_FRA1(uox)
      ,CHFt_CONST_FRA1(uoy)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADJRHSVELCTRL 

#ifndef GUARDADJRHSLOGSPDCTRL 
#define GUARDADJRHSLOGSPDCTRL 
// Prototype for Fortran procedure adjrhslogspdctrl ...
//
void FORTRAN_NAME( ADJRHSLOGSPDCTRL ,adjrhslogspdctrl )(
      CHFp_FRA1(rhsx)
      ,CHFp_FRA1(rhsy)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(umx)
      ,CHFp_CONST_FRA1(umy)
      ,CHFp_CONST_FRA1(uox)
      ,CHFp_CONST_FRA1(uoy)
      ,CHFp_BOX(box) );

#define FORT_ADJRHSLOGSPDCTRL FORTRAN_NAME( inlineADJRHSLOGSPDCTRL, inlineADJRHSLOGSPDCTRL)
#define FORTNT_ADJRHSLOGSPDCTRL FORTRAN_NAME( ADJRHSLOGSPDCTRL, adjrhslogspdctrl)

inline void FORTRAN_NAME(inlineADJRHSLOGSPDCTRL, inlineADJRHSLOGSPDCTRL)(
      CHFp_FRA1(rhsx)
      ,CHFp_FRA1(rhsy)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(umx)
      ,CHFp_CONST_FRA1(umy)
      ,CHFp_CONST_FRA1(uox)
      ,CHFp_CONST_FRA1(uoy)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADJRHSLOGSPDCTRL");
 FORTRAN_NAME( ADJRHSLOGSPDCTRL ,adjrhslogspdctrl )(
      CHFt_FRA1(rhsx)
      ,CHFt_FRA1(rhsy)
      ,CHFt_FRA1(misfit)
      ,CHFt_CONST_FRA1(umx)
      ,CHFt_CONST_FRA1(umy)
      ,CHFt_CONST_FRA1(uox)
      ,CHFt_CONST_FRA1(uoy)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADJRHSLOGSPDCTRL 

#ifndef GUARDADJRHSMASSCTRL 
#define GUARDADJRHSMASSCTRL 
// Prototype for Fortran procedure adjrhsmassctrl ...
//
void FORTRAN_NAME( ADJRHSMASSCTRL ,adjrhsmassctrl )(
      CHFp_FRA(rhs)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_ADJRHSMASSCTRL FORTRAN_NAME( inlineADJRHSMASSCTRL, inlineADJRHSMASSCTRL)
#define FORTNT_ADJRHSMASSCTRL FORTRAN_NAME( ADJRHSMASSCTRL, adjrhsmassctrl)

inline void FORTRAN_NAME(inlineADJRHSMASSCTRL, inlineADJRHSMASSCTRL)(
      CHFp_FRA(rhs)
      ,CHFp_FRA1(misfit)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADJRHSMASSCTRL");
 FORTRAN_NAME( ADJRHSMASSCTRL ,adjrhsmassctrl )(
      CHFt_FRA(rhs)
      ,CHFt_FRA1(misfit)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADJRHSMASSCTRL 

#ifndef GUARDBOUNDEXPCTRL 
#define GUARDBOUNDEXPCTRL 
// Prototype for Fortran procedure boundexpctrl ...
//
void FORTRAN_NAME( BOUNDEXPCTRL ,boundexpctrl )(
      CHFp_FRA1(x)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) );

#define FORT_BOUNDEXPCTRL FORTRAN_NAME( inlineBOUNDEXPCTRL, inlineBOUNDEXPCTRL)
#define FORTNT_BOUNDEXPCTRL FORTRAN_NAME( BOUNDEXPCTRL, boundexpctrl)

inline void FORTRAN_NAME(inlineBOUNDEXPCTRL, inlineBOUNDEXPCTRL)(
      CHFp_FRA1(x)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_BOUNDEXPCTRL");
 FORTRAN_NAME( BOUNDEXPCTRL ,boundexpctrl )(
      CHFt_FRA1(x)
      ,CHFt_CONST_FRA1(a)
      ,CHFt_CONST_REAL(lb)
      ,CHFt_CONST_REAL(ub)
      ,CHFt_BOX(box) );
}
#endif  // GUARDBOUNDEXPCTRL 

#ifndef GUARDEXPCTRL 
#define GUARDEXPCTRL 
// Prototype for Fortran procedure expctrl ...
//
void FORTRAN_NAME( EXPCTRL ,expctrl )(
      CHFp_FRA1(x)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_BOX(box) );

#define FORT_EXPCTRL FORTRAN_NAME( inlineEXPCTRL, inlineEXPCTRL)
#define FORTNT_EXPCTRL FORTRAN_NAME( EXPCTRL, expctrl)

inline void FORTRAN_NAME(inlineEXPCTRL, inlineEXPCTRL)(
      CHFp_FRA1(x)
      ,CHFp_CONST_FRA1(a)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EXPCTRL");
 FORTRAN_NAME( EXPCTRL ,expctrl )(
      CHFt_FRA1(x)
      ,CHFt_CONST_FRA1(a)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEXPCTRL 

#ifndef GUARDBOUNDCTRL 
#define GUARDBOUNDCTRL 
// Prototype for Fortran procedure boundctrl ...
//
void FORTRAN_NAME( BOUNDCTRL ,boundctrl )(
      CHFp_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) );

#define FORT_BOUNDCTRL FORTRAN_NAME( inlineBOUNDCTRL, inlineBOUNDCTRL)
#define FORTNT_BOUNDCTRL FORTRAN_NAME( BOUNDCTRL, boundctrl)

inline void FORTRAN_NAME(inlineBOUNDCTRL, inlineBOUNDCTRL)(
      CHFp_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_BOUNDCTRL");
 FORTRAN_NAME( BOUNDCTRL ,boundctrl )(
      CHFt_FRA1(x)
      ,CHFt_CONST_REAL(lb)
      ,CHFt_CONST_REAL(ub)
      ,CHFt_BOX(box) );
}
#endif  // GUARDBOUNDCTRL 

#ifndef GUARDINCRBOUNDCTRL 
#define GUARDINCRBOUNDCTRL 
// Prototype for Fortran procedure incrboundctrl ...
//
void FORTRAN_NAME( INCRBOUNDCTRL ,incrboundctrl )(
      CHFp_FRA1(z)
      ,CHFp_CONST_FRA1(y)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) );

#define FORT_INCRBOUNDCTRL FORTRAN_NAME( inlineINCRBOUNDCTRL, inlineINCRBOUNDCTRL)
#define FORTNT_INCRBOUNDCTRL FORTRAN_NAME( INCRBOUNDCTRL, incrboundctrl)

inline void FORTRAN_NAME(inlineINCRBOUNDCTRL, inlineINCRBOUNDCTRL)(
      CHFp_FRA1(z)
      ,CHFp_CONST_FRA1(y)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INCRBOUNDCTRL");
 FORTRAN_NAME( INCRBOUNDCTRL ,incrboundctrl )(
      CHFt_FRA1(z)
      ,CHFt_CONST_FRA1(y)
      ,CHFt_CONST_FRA1(x)
      ,CHFt_CONST_REAL(lb)
      ,CHFt_CONST_REAL(ub)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINCRBOUNDCTRL 

#ifndef GUARDGRADBARRIERCTRL 
#define GUARDGRADBARRIERCTRL 
// Prototype for Fortran procedure gradbarrierctrl ...
//
void FORTRAN_NAME( GRADBARRIERCTRL ,gradbarrierctrl )(
      CHFp_FRA1(p)
      ,CHFp_FRA1(g)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(b)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_BOX(box) );

#define FORT_GRADBARRIERCTRL FORTRAN_NAME( inlineGRADBARRIERCTRL, inlineGRADBARRIERCTRL)
#define FORTNT_GRADBARRIERCTRL FORTRAN_NAME( GRADBARRIERCTRL, gradbarrierctrl)

inline void FORTRAN_NAME(inlineGRADBARRIERCTRL, inlineGRADBARRIERCTRL)(
      CHFp_FRA1(p)
      ,CHFp_FRA1(g)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(b)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GRADBARRIERCTRL");
 FORTRAN_NAME( GRADBARRIERCTRL ,gradbarrierctrl )(
      CHFt_FRA1(p)
      ,CHFt_FRA1(g)
      ,CHFt_CONST_FRA1(x)
      ,CHFt_CONST_REAL(b)
      ,CHFt_CONST_REAL(tol)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGRADBARRIERCTRL 

#ifndef GUARDMULTHATCTRL 
#define GUARDMULTHATCTRL 
// Prototype for Fortran procedure multhatctrl ...
//
void FORTRAN_NAME( MULTHATCTRL ,multhatctrl )(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) );

#define FORT_MULTHATCTRL FORTRAN_NAME( inlineMULTHATCTRL, inlineMULTHATCTRL)
#define FORTNT_MULTHATCTRL FORTRAN_NAME( MULTHATCTRL, multhatctrl)

inline void FORTRAN_NAME(inlineMULTHATCTRL, inlineMULTHATCTRL)(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MULTHATCTRL");
 FORTRAN_NAME( MULTHATCTRL ,multhatctrl )(
      CHFt_FRA1(fab)
      ,CHFt_CONST_FRA1(x)
      ,CHFt_CONST_REAL(lb)
      ,CHFt_CONST_REAL(ub)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMULTHATCTRL 

#ifndef GUARDHARDPOINTINCTRL 
#define GUARDHARDPOINTINCTRL 
// Prototype for Fortran procedure hardpointinctrl ...
//
void FORTRAN_NAME( HARDPOINTINCTRL ,hardpointinctrl )(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) );

#define FORT_HARDPOINTINCTRL FORTRAN_NAME( inlineHARDPOINTINCTRL, inlineHARDPOINTINCTRL)
#define FORTNT_HARDPOINTINCTRL FORTRAN_NAME( HARDPOINTINCTRL, hardpointinctrl)

inline void FORTRAN_NAME(inlineHARDPOINTINCTRL, inlineHARDPOINTINCTRL)(
      CHFp_FRA1(fab)
      ,CHFp_CONST_FRA1(x)
      ,CHFp_CONST_REAL(lb)
      ,CHFp_CONST_REAL(ub)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HARDPOINTINCTRL");
 FORTRAN_NAME( HARDPOINTINCTRL ,hardpointinctrl )(
      CHFt_FRA1(fab)
      ,CHFt_CONST_FRA1(x)
      ,CHFt_CONST_REAL(lb)
      ,CHFt_CONST_REAL(ub)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHARDPOINTINCTRL 

#ifndef GUARDCADOTBCTRL 
#define GUARDCADOTBCTRL 
// Prototype for Fortran procedure cadotbctrl ...
//
void FORTRAN_NAME( CADOTBCTRL ,cadotbctrl )(
      CHFp_FRA1(r)
      ,CHFp_CONST_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_FRA(b)
      ,CHFp_BOX(box) );

#define FORT_CADOTBCTRL FORTRAN_NAME( inlineCADOTBCTRL, inlineCADOTBCTRL)
#define FORTNT_CADOTBCTRL FORTRAN_NAME( CADOTBCTRL, cadotbctrl)

inline void FORTRAN_NAME(inlineCADOTBCTRL, inlineCADOTBCTRL)(
      CHFp_FRA1(r)
      ,CHFp_CONST_FRA1(c)
      ,CHFp_CONST_FRA(a)
      ,CHFp_CONST_FRA(b)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CADOTBCTRL");
 FORTRAN_NAME( CADOTBCTRL ,cadotbctrl )(
      CHFt_FRA1(r)
      ,CHFt_CONST_FRA1(c)
      ,CHFt_CONST_FRA(a)
      ,CHFt_CONST_FRA(b)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCADOTBCTRL 

#ifndef GUARDLIMITINCRCTRL 
#define GUARDLIMITINCRCTRL 
// Prototype for Fortran procedure limitincrctrl ...
//
void FORTRAN_NAME( LIMITINCRCTRL ,limitincrctrl )(
      CHFp_FRA1(h)
      ,CHFp_CONST_FRA1(dh)
      ,CHFp_CONST_FRA1(ho)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) );

#define FORT_LIMITINCRCTRL FORTRAN_NAME( inlineLIMITINCRCTRL, inlineLIMITINCRCTRL)
#define FORTNT_LIMITINCRCTRL FORTRAN_NAME( LIMITINCRCTRL, limitincrctrl)

inline void FORTRAN_NAME(inlineLIMITINCRCTRL, inlineLIMITINCRCTRL)(
      CHFp_FRA1(h)
      ,CHFp_CONST_FRA1(dh)
      ,CHFp_CONST_FRA1(ho)
      ,CHFp_CONST_REAL(limit)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LIMITINCRCTRL");
 FORTRAN_NAME( LIMITINCRCTRL ,limitincrctrl )(
      CHFt_FRA1(h)
      ,CHFt_CONST_FRA1(dh)
      ,CHFt_CONST_FRA1(ho)
      ,CHFt_CONST_REAL(limit)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLIMITINCRCTRL 

#ifndef GUARDUPDATEHCTRLB 
#define GUARDUPDATEHCTRLB 
// Prototype for Fortran procedure updatehctrlb ...
//
void FORTRAN_NAME( UPDATEHCTRLB ,updatehctrlb )(
      CHFp_FRA1(h)
      ,CHFp_FRA1(topg)
      ,CHFp_FRA1(dh)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(ho)
      ,CHFp_CONST_REAL(r)
      ,CHFp_CONST_REAL(dhmax)
      ,CHFp_CONST_REAL(hmin)
      ,CHFp_BOX(box) );

#define FORT_UPDATEHCTRLB FORTRAN_NAME( inlineUPDATEHCTRLB, inlineUPDATEHCTRLB)
#define FORTNT_UPDATEHCTRLB FORTRAN_NAME( UPDATEHCTRLB, updatehctrlb)

inline void FORTRAN_NAME(inlineUPDATEHCTRLB, inlineUPDATEHCTRLB)(
      CHFp_FRA1(h)
      ,CHFp_FRA1(topg)
      ,CHFp_FRA1(dh)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(ho)
      ,CHFp_CONST_REAL(r)
      ,CHFp_CONST_REAL(dhmax)
      ,CHFp_CONST_REAL(hmin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATEHCTRLB");
 FORTRAN_NAME( UPDATEHCTRLB ,updatehctrlb )(
      CHFt_FRA1(h)
      ,CHFt_FRA1(topg)
      ,CHFt_FRA1(dh)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_FRA1(ho)
      ,CHFt_CONST_REAL(r)
      ,CHFt_CONST_REAL(dhmax)
      ,CHFt_CONST_REAL(hmin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATEHCTRLB 

#ifndef GUARDUPDATEHCTRLC 
#define GUARDUPDATEHCTRLC 
// Prototype for Fortran procedure updatehctrlc ...
//
void FORTRAN_NAME( UPDATEHCTRLC ,updatehctrlc )(
      CHFp_FRA1(h)
      ,CHFp_FRA1(topg)
      ,CHFp_FRA1(dh)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(ho)
      ,CHFp_CONST_REAL(r)
      ,CHFp_CONST_REAL(dhmax)
      ,CHFp_CONST_REAL(hmin)
      ,CHFp_BOX(box) );

#define FORT_UPDATEHCTRLC FORTRAN_NAME( inlineUPDATEHCTRLC, inlineUPDATEHCTRLC)
#define FORTNT_UPDATEHCTRLC FORTRAN_NAME( UPDATEHCTRLC, updatehctrlc)

inline void FORTRAN_NAME(inlineUPDATEHCTRLC, inlineUPDATEHCTRLC)(
      CHFp_FRA1(h)
      ,CHFp_FRA1(topg)
      ,CHFp_FRA1(dh)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_CONST_FRA1(ho)
      ,CHFp_CONST_REAL(r)
      ,CHFp_CONST_REAL(dhmax)
      ,CHFp_CONST_REAL(hmin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATEHCTRLC");
 FORTRAN_NAME( UPDATEHCTRLC ,updatehctrlc )(
      CHFt_FRA1(h)
      ,CHFt_FRA1(topg)
      ,CHFt_FRA1(dh)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_CONST_FRA1(ho)
      ,CHFt_CONST_REAL(r)
      ,CHFt_CONST_REAL(dhmax)
      ,CHFt_CONST_REAL(hmin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATEHCTRLC 

#ifndef GUARDCONVOLVECTRLRB 
#define GUARDCONVOLVECTRLRB 
// Prototype for Fortran procedure convolvectrlrb ...
//
void FORTRAN_NAME( CONVOLVECTRLRB ,convolvectrlrb )(
      CHFp_FRA(u)
      ,CHFp_CONST_INT(rb)
      ,CHFp_BOX(box) );

#define FORT_CONVOLVECTRLRB FORTRAN_NAME( inlineCONVOLVECTRLRB, inlineCONVOLVECTRLRB)
#define FORTNT_CONVOLVECTRLRB FORTRAN_NAME( CONVOLVECTRLRB, convolvectrlrb)

inline void FORTRAN_NAME(inlineCONVOLVECTRLRB, inlineCONVOLVECTRLRB)(
      CHFp_FRA(u)
      ,CHFp_CONST_INT(rb)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONVOLVECTRLRB");
 FORTRAN_NAME( CONVOLVECTRLRB ,convolvectrlrb )(
      CHFt_FRA(u)
      ,CHFt_CONST_INT(rb)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONVOLVECTRLRB 

#ifndef GUARDCONVOLVECTRL 
#define GUARDCONVOLVECTRL 
// Prototype for Fortran procedure convolvectrl ...
//
void FORTRAN_NAME( CONVOLVECTRL ,convolvectrl )(
      CHFp_FRA(u)
      ,CHFp_FRA(r)
      ,CHFp_BOX(box) );

#define FORT_CONVOLVECTRL FORTRAN_NAME( inlineCONVOLVECTRL, inlineCONVOLVECTRL)
#define FORTNT_CONVOLVECTRL FORTRAN_NAME( CONVOLVECTRL, convolvectrl)

inline void FORTRAN_NAME(inlineCONVOLVECTRL, inlineCONVOLVECTRL)(
      CHFp_FRA(u)
      ,CHFp_FRA(r)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONVOLVECTRL");
 FORTRAN_NAME( CONVOLVECTRL ,convolvectrl )(
      CHFt_FRA(u)
      ,CHFt_FRA(r)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONVOLVECTRL 

}

#endif
