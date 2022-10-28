#ifndef _L1L2CONSTITUTIVERELATIONF_F_H_
#define _L1L2CONSTITUTIVERELATIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDL1L2LIMITFAB 
#define GUARDL1L2LIMITFAB 
// Prototype for Fortran procedure L1L2LIMITFAB ...
//
void FORTRAN_NAME( L1L2LIMITFAB ,l1l2limitfab )(
      CHFp_FRA1(fab)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(mv) );

#define FORT_L1L2LIMITFAB FORTRAN_NAME( inlineL1L2LIMITFAB, inlineL1L2LIMITFAB)
#define FORTNT_L1L2LIMITFAB FORTRAN_NAME( L1L2LIMITFAB, l1l2limitfab)

inline void FORTRAN_NAME(inlineL1L2LIMITFAB, inlineL1L2LIMITFAB)(
      CHFp_FRA1(fab)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(mv) )
{
 CH_TIMELEAF("FORT_L1L2LIMITFAB");
 FORTRAN_NAME( L1L2LIMITFAB ,l1l2limitfab )(
      CHFt_FRA1(fab)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(mv) );
}
#endif  // GUARDL1L2LIMITFAB 

#ifndef GUARDL1L2SUPRESSFAB 
#define GUARDL1L2SUPRESSFAB 
// Prototype for Fortran procedure L1L2SUPRESSFAB ...
//
void FORTRAN_NAME( L1L2SUPRESSFAB ,l1l2supressfab )(
      CHFp_FRA1(fab)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(mv) );

#define FORT_L1L2SUPRESSFAB FORTRAN_NAME( inlineL1L2SUPRESSFAB, inlineL1L2SUPRESSFAB)
#define FORTNT_L1L2SUPRESSFAB FORTRAN_NAME( L1L2SUPRESSFAB, l1l2supressfab)

inline void FORTRAN_NAME(inlineL1L2SUPRESSFAB, inlineL1L2SUPRESSFAB)(
      CHFp_FRA1(fab)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(mv) )
{
 CH_TIMELEAF("FORT_L1L2SUPRESSFAB");
 FORTRAN_NAME( L1L2SUPRESSFAB ,l1l2supressfab )(
      CHFt_FRA1(fab)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(mv) );
}
#endif  // GUARDL1L2SUPRESSFAB 

#ifndef GUARDANALYTICL1L2MU 
#define GUARDANALYTICL1L2MU 
// Prototype for Fortran procedure ANALYTICL1L2MU ...
//
void FORTRAN_NAME( ANALYTICL1L2MU ,analyticl1l2mu )(
      CHFp_FRA1(mu)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(epsSqr)
      ,CHFp_CONST_FRA1(torSqr)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(sigma)
      ,CHFp_CONST_REAL(nExponent)
      ,CHFp_CONST_REAL(epsSqr0) );

#define FORT_ANALYTICL1L2MU FORTRAN_NAME( inlineANALYTICL1L2MU, inlineANALYTICL1L2MU)
#define FORTNT_ANALYTICL1L2MU FORTRAN_NAME( ANALYTICL1L2MU, analyticl1l2mu)

inline void FORTRAN_NAME(inlineANALYTICL1L2MU, inlineANALYTICL1L2MU)(
      CHFp_FRA1(mu)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(epsSqr)
      ,CHFp_CONST_FRA1(torSqr)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(sigma)
      ,CHFp_CONST_REAL(nExponent)
      ,CHFp_CONST_REAL(epsSqr0) )
{
 CH_TIMELEAF("FORT_ANALYTICL1L2MU");
 FORTRAN_NAME( ANALYTICL1L2MU ,analyticl1l2mu )(
      CHFt_FRA1(mu)
      ,CHFt_CONST_FRA1(A)
      ,CHFt_CONST_FRA1(epsSqr)
      ,CHFt_CONST_FRA1(torSqr)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(sigma)
      ,CHFt_CONST_REAL(nExponent)
      ,CHFt_CONST_REAL(epsSqr0) );
}
#endif  // GUARDANALYTICL1L2MU 

#ifndef GUARDCOMPUTEL1L2RES 
#define GUARDCOMPUTEL1L2RES 
// Prototype for Fortran procedure COMPUTEL1L2RES ...
//
void FORTRAN_NAME( COMPUTEL1L2RES ,computel1l2res )(
      CHFp_REAL(res)
      ,CHFp_REAL(mu)
      ,CHFp_CONST_REAL(A)
      ,CHFp_CONST_REAL(esqr)
      ,CHFp_CONST_REAL(tsqr)
      ,CHFp_CONST_REAL(p) );

#define FORT_COMPUTEL1L2RES FORTRAN_NAME( inlineCOMPUTEL1L2RES, inlineCOMPUTEL1L2RES)
#define FORTNT_COMPUTEL1L2RES FORTRAN_NAME( COMPUTEL1L2RES, computel1l2res)

inline void FORTRAN_NAME(inlineCOMPUTEL1L2RES, inlineCOMPUTEL1L2RES)(
      CHFp_REAL(res)
      ,CHFp_REAL(mu)
      ,CHFp_CONST_REAL(A)
      ,CHFp_CONST_REAL(esqr)
      ,CHFp_CONST_REAL(tsqr)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_COMPUTEL1L2RES");
 FORTRAN_NAME( COMPUTEL1L2RES ,computel1l2res )(
      CHFt_REAL(res)
      ,CHFt_REAL(mu)
      ,CHFt_CONST_REAL(A)
      ,CHFt_CONST_REAL(esqr)
      ,CHFt_CONST_REAL(tsqr)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDCOMPUTEL1L2RES 

#ifndef GUARDL1L2SECANTSOLVE 
#define GUARDL1L2SECANTSOLVE 
// Prototype for Fortran procedure L1L2SECANTSOLVE ...
//
void FORTRAN_NAME( L1L2SECANTSOLVE ,l1l2secantsolve )(
      CHFp_REAL(mu)
      ,CHFp_REAL(mua)
      ,CHFp_CONST_REAL(A)
      ,CHFp_CONST_REAL(esqr)
      ,CHFp_CONST_REAL(tsqr)
      ,CHFp_CONST_REAL(p)
      ,CHFp_INT(niter)
      ,CHFp_CONST_REAL(tol) );

#define FORT_L1L2SECANTSOLVE FORTRAN_NAME( inlineL1L2SECANTSOLVE, inlineL1L2SECANTSOLVE)
#define FORTNT_L1L2SECANTSOLVE FORTRAN_NAME( L1L2SECANTSOLVE, l1l2secantsolve)

inline void FORTRAN_NAME(inlineL1L2SECANTSOLVE, inlineL1L2SECANTSOLVE)(
      CHFp_REAL(mu)
      ,CHFp_REAL(mua)
      ,CHFp_CONST_REAL(A)
      ,CHFp_CONST_REAL(esqr)
      ,CHFp_CONST_REAL(tsqr)
      ,CHFp_CONST_REAL(p)
      ,CHFp_INT(niter)
      ,CHFp_CONST_REAL(tol) )
{
 CH_TIMELEAF("FORT_L1L2SECANTSOLVE");
 FORTRAN_NAME( L1L2SECANTSOLVE ,l1l2secantsolve )(
      CHFt_REAL(mu)
      ,CHFt_REAL(mua)
      ,CHFt_CONST_REAL(A)
      ,CHFt_CONST_REAL(esqr)
      ,CHFt_CONST_REAL(tsqr)
      ,CHFt_CONST_REAL(p)
      ,CHFt_INT(niter)
      ,CHFt_CONST_REAL(tol) );
}
#endif  // GUARDL1L2SECANTSOLVE 

#ifndef GUARDCOMPUTEL1L2MU 
#define GUARDCOMPUTEL1L2MU 
// Prototype for Fortran procedure COMPUTEL1L2MU ...
//
void FORTRAN_NAME( COMPUTEL1L2MU ,computel1l2mu )(
      CHFp_FRA1(mu)
      ,CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(epsSqr)
      ,CHFp_CONST_FRA1(torSqr)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(sigma)
      ,CHFp_CONST_REAL(nExponent)
      ,CHFp_CONST_REAL(epsSqr0)
      ,CHFp_CONST_REAL(delta)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_CONST_INT(maxIter) );

#define FORT_COMPUTEL1L2MU FORTRAN_NAME( inlineCOMPUTEL1L2MU, inlineCOMPUTEL1L2MU)
#define FORTNT_COMPUTEL1L2MU FORTRAN_NAME( COMPUTEL1L2MU, computel1l2mu)

inline void FORTRAN_NAME(inlineCOMPUTEL1L2MU, inlineCOMPUTEL1L2MU)(
      CHFp_FRA1(mu)
      ,CHFp_FRA1(res)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(epsSqr)
      ,CHFp_CONST_FRA1(torSqr)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(sigma)
      ,CHFp_CONST_REAL(nExponent)
      ,CHFp_CONST_REAL(epsSqr0)
      ,CHFp_CONST_REAL(delta)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_CONST_INT(maxIter) )
{
 CH_TIMELEAF("FORT_COMPUTEL1L2MU");
 FORTRAN_NAME( COMPUTEL1L2MU ,computel1l2mu )(
      CHFt_FRA1(mu)
      ,CHFt_FRA1(res)
      ,CHFt_CONST_FRA1(A)
      ,CHFt_CONST_FRA1(epsSqr)
      ,CHFt_CONST_FRA1(torSqr)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(sigma)
      ,CHFt_CONST_REAL(nExponent)
      ,CHFt_CONST_REAL(epsSqr0)
      ,CHFt_CONST_REAL(delta)
      ,CHFt_CONST_REAL(tol)
      ,CHFt_CONST_INT(maxIter) );
}
#endif  // GUARDCOMPUTEL1L2MU 

#ifndef GUARDL1L2PHIINTEGRANDS 
#define GUARDL1L2PHIINTEGRANDS 
// Prototype for Fortran procedure L1L2PHIINTEGRANDS ...
//
void FORTRAN_NAME( L1L2PHIINTEGRANDS ,l1l2phiintegrands )(
      CHFp_FRA1(sixx)
      ,CHFp_FRA1(sixy)
      ,CHFp_FRA1(siyy)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_FRA(derivs)
      ,CHFp_INT(dudxComp)
      ,CHFp_INT(dudyComp)
      ,CHFp_INT(dvdxComp)
      ,CHFp_INT(dvdyComp)
      ,CHFp_BOX(box) );

#define FORT_L1L2PHIINTEGRANDS FORTRAN_NAME( inlineL1L2PHIINTEGRANDS, inlineL1L2PHIINTEGRANDS)
#define FORTNT_L1L2PHIINTEGRANDS FORTRAN_NAME( L1L2PHIINTEGRANDS, l1l2phiintegrands)

inline void FORTRAN_NAME(inlineL1L2PHIINTEGRANDS, inlineL1L2PHIINTEGRANDS)(
      CHFp_FRA1(sixx)
      ,CHFp_FRA1(sixy)
      ,CHFp_FRA1(siyy)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_FRA(derivs)
      ,CHFp_INT(dudxComp)
      ,CHFp_INT(dudyComp)
      ,CHFp_INT(dvdxComp)
      ,CHFp_INT(dvdyComp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_L1L2PHIINTEGRANDS");
 FORTRAN_NAME( L1L2PHIINTEGRANDS ,l1l2phiintegrands )(
      CHFt_FRA1(sixx)
      ,CHFt_FRA1(sixy)
      ,CHFt_FRA1(siyy)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA1(mu)
      ,CHFt_FRA(derivs)
      ,CHFt_INT(dudxComp)
      ,CHFt_INT(dudyComp)
      ,CHFt_INT(dvdxComp)
      ,CHFt_INT(dvdyComp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDL1L2PHIINTEGRANDS 

#ifndef GUARDL1L2UINTEGRANDS 
#define GUARDL1L2UINTEGRANDS 
// Prototype for Fortran procedure L1L2UINTEGRANDS ...
//
void FORTRAN_NAME( L1L2UINTEGRANDS ,l1l2uintegrands )(
      CHFp_FRA1(fzx)
      ,CHFp_FRA1(fzy)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_CONST_FRA1(epssqr)
      ,CHFp_CONST_FRA1(phizx)
      ,CHFp_CONST_FRA1(phizy)
      ,CHFp_CONST_REAL(n)
      ,CHFp_BOX(box) );

#define FORT_L1L2UINTEGRANDS FORTRAN_NAME( inlineL1L2UINTEGRANDS, inlineL1L2UINTEGRANDS)
#define FORTNT_L1L2UINTEGRANDS FORTRAN_NAME( L1L2UINTEGRANDS, l1l2uintegrands)

inline void FORTRAN_NAME(inlineL1L2UINTEGRANDS, inlineL1L2UINTEGRANDS)(
      CHFp_FRA1(fzx)
      ,CHFp_FRA1(fzy)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_CONST_FRA1(epssqr)
      ,CHFp_CONST_FRA1(phizx)
      ,CHFp_CONST_FRA1(phizy)
      ,CHFp_CONST_REAL(n)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_L1L2UINTEGRANDS");
 FORTRAN_NAME( L1L2UINTEGRANDS ,l1l2uintegrands )(
      CHFt_FRA1(fzx)
      ,CHFt_FRA1(fzy)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA1(A)
      ,CHFt_CONST_FRA1(mu)
      ,CHFt_CONST_FRA1(epssqr)
      ,CHFt_CONST_FRA1(phizx)
      ,CHFt_CONST_FRA1(phizy)
      ,CHFt_CONST_REAL(n)
      ,CHFt_BOX(box) );
}
#endif  // GUARDL1L2UINTEGRANDS 

#ifndef GUARDL1L2UIGRAND 
#define GUARDL1L2UIGRAND 
// Prototype for Fortran procedure L1L2UIGRAND ...
//
void FORTRAN_NAME( L1L2UIGRAND ,l1l2uigrand )(
      CHFp_FRA(fzi)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_CONST_FRA1(epssqr)
      ,CHFp_CONST_FRA(phizi)
      ,CHFp_CONST_REAL(n)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(box) );

#define FORT_L1L2UIGRAND FORTRAN_NAME( inlineL1L2UIGRAND, inlineL1L2UIGRAND)
#define FORTNT_L1L2UIGRAND FORTRAN_NAME( L1L2UIGRAND, l1l2uigrand)

inline void FORTRAN_NAME(inlineL1L2UIGRAND, inlineL1L2UIGRAND)(
      CHFp_FRA(fzi)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA1(A)
      ,CHFp_CONST_FRA1(mu)
      ,CHFp_CONST_FRA1(epssqr)
      ,CHFp_CONST_FRA(phizi)
      ,CHFp_CONST_REAL(n)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_L1L2UIGRAND");
 FORTRAN_NAME( L1L2UIGRAND ,l1l2uigrand )(
      CHFt_FRA(fzi)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA1(A)
      ,CHFt_CONST_FRA1(mu)
      ,CHFt_CONST_FRA1(epssqr)
      ,CHFt_CONST_FRA(phizi)
      ,CHFt_CONST_REAL(n)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDL1L2UIGRAND 

#ifndef GUARDL1L2MODLIMIT 
#define GUARDL1L2MODLIMIT 
// Prototype for Fortran procedure L1L2MODLIMIT ...
//
void FORTRAN_NAME( L1L2MODLIMIT ,l1l2modlimit )(
      CHFp_FRA(fab)
      ,CHFp_CONST_REAL(lim)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_L1L2MODLIMIT FORTRAN_NAME( inlineL1L2MODLIMIT, inlineL1L2MODLIMIT)
#define FORTNT_L1L2MODLIMIT FORTRAN_NAME( L1L2MODLIMIT, l1l2modlimit)

inline void FORTRAN_NAME(inlineL1L2MODLIMIT, inlineL1L2MODLIMIT)(
      CHFp_FRA(fab)
      ,CHFp_CONST_REAL(lim)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_L1L2MODLIMIT");
 FORTRAN_NAME( L1L2MODLIMIT ,l1l2modlimit )(
      CHFt_FRA(fab)
      ,CHFt_CONST_REAL(lim)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDL1L2MODLIMIT 

#ifndef GUARDL1L2ADDGSTRESS 
#define GUARDL1L2ADDGSTRESS 
// Prototype for Fortran procedure L1L2ADDGSTRESS ...
//
void FORTRAN_NAME( L1L2ADDGSTRESS ,l1l2addgstress )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_REAL(sigma)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_L1L2ADDGSTRESS FORTRAN_NAME( inlineL1L2ADDGSTRESS, inlineL1L2ADDGSTRESS)
#define FORTNT_L1L2ADDGSTRESS FORTRAN_NAME( L1L2ADDGSTRESS, l1l2addgstress)

inline void FORTRAN_NAME(inlineL1L2ADDGSTRESS, inlineL1L2ADDGSTRESS)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_REAL(sigma)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_L1L2ADDGSTRESS");
 FORTRAN_NAME( L1L2ADDGSTRESS ,l1l2addgstress )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(G)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_REAL(sigma)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDL1L2ADDGSTRESS 

#ifndef GUARDL1L2PHITILDESQR 
#define GUARDL1L2PHITILDESQR 
// Prototype for Fortran procedure L1L2PHITILDESQR ...
//
void FORTRAN_NAME( L1L2PHITILDESQR ,l1l2phitildesqr )(
      CHFp_FRA1(ptsqr)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA(grads)
      ,CHFp_CONST_REAL(rhog)
      ,CHFp_INT(ncomp)
      ,CHFp_BOX(box) );

#define FORT_L1L2PHITILDESQR FORTRAN_NAME( inlineL1L2PHITILDESQR, inlineL1L2PHITILDESQR)
#define FORTNT_L1L2PHITILDESQR FORTRAN_NAME( L1L2PHITILDESQR, l1l2phitildesqr)

inline void FORTRAN_NAME(inlineL1L2PHITILDESQR, inlineL1L2PHITILDESQR)(
      CHFp_FRA1(ptsqr)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA(grads)
      ,CHFp_CONST_REAL(rhog)
      ,CHFp_INT(ncomp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_L1L2PHITILDESQR");
 FORTRAN_NAME( L1L2PHITILDESQR ,l1l2phitildesqr )(
      CHFt_FRA1(ptsqr)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA(grads)
      ,CHFt_CONST_REAL(rhog)
      ,CHFt_INT(ncomp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDL1L2PHITILDESQR 

#ifndef GUARDL1L2DIFFFACTOR 
#define GUARDL1L2DIFFFACTOR 
// Prototype for Fortran procedure L1L2DIFFFACTOR ...
//
void FORTRAN_NAME( L1L2DIFFFACTOR ,l1l2difffactor )(
      CHFp_FRA1(D)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA(grads)
      ,CHFp_CONST_REAL(n)
      ,CHFp_CONST_REAL(rhog)
      ,CHFp_BOX(box) );

#define FORT_L1L2DIFFFACTOR FORTRAN_NAME( inlineL1L2DIFFFACTOR, inlineL1L2DIFFFACTOR)
#define FORTNT_L1L2DIFFFACTOR FORTRAN_NAME( L1L2DIFFFACTOR, l1l2difffactor)

inline void FORTRAN_NAME(inlineL1L2DIFFFACTOR, inlineL1L2DIFFFACTOR)(
      CHFp_FRA1(D)
      ,CHFp_CONST_FRA1(H)
      ,CHFp_CONST_FRA(grads)
      ,CHFp_CONST_REAL(n)
      ,CHFp_CONST_REAL(rhog)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_L1L2DIFFFACTOR");
 FORTRAN_NAME( L1L2DIFFFACTOR ,l1l2difffactor )(
      CHFt_FRA1(D)
      ,CHFt_CONST_FRA1(H)
      ,CHFt_CONST_FRA(grads)
      ,CHFt_CONST_REAL(n)
      ,CHFt_CONST_REAL(rhog)
      ,CHFt_BOX(box) );
}
#endif  // GUARDL1L2DIFFFACTOR 

#ifndef GUARDL1L2COMPUTEDIFFUSIVITY 
#define GUARDL1L2COMPUTEDIFFUSIVITY 
// Prototype for Fortran procedure L1L2COMPUTEDIFFUSIVITY ...
//
void FORTRAN_NAME( L1L2COMPUTEDIFFUSIVITY ,l1l2computediffusivity )(
      CHFp_FRA1(D)
      ,CHFp_CONST_FRA1(ux)
      ,CHFp_CONST_FRA1(uy)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(s)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dy)
      ,CHFp_BOX(box) );

#define FORT_L1L2COMPUTEDIFFUSIVITY FORTRAN_NAME( inlineL1L2COMPUTEDIFFUSIVITY, inlineL1L2COMPUTEDIFFUSIVITY)
#define FORTNT_L1L2COMPUTEDIFFUSIVITY FORTRAN_NAME( L1L2COMPUTEDIFFUSIVITY, l1l2computediffusivity)

inline void FORTRAN_NAME(inlineL1L2COMPUTEDIFFUSIVITY, inlineL1L2COMPUTEDIFFUSIVITY)(
      CHFp_FRA1(D)
      ,CHFp_CONST_FRA1(ux)
      ,CHFp_CONST_FRA1(uy)
      ,CHFp_CONST_FRA1(thck)
      ,CHFp_CONST_FRA1(s)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dy)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_L1L2COMPUTEDIFFUSIVITY");
 FORTRAN_NAME( L1L2COMPUTEDIFFUSIVITY ,l1l2computediffusivity )(
      CHFt_FRA1(D)
      ,CHFt_CONST_FRA1(ux)
      ,CHFt_CONST_FRA1(uy)
      ,CHFt_CONST_FRA1(thck)
      ,CHFt_CONST_FRA1(s)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dy)
      ,CHFt_BOX(box) );
}
#endif  // GUARDL1L2COMPUTEDIFFUSIVITY 

#ifndef GUARDL1L2CELLTOFACEHARMONIC 
#define GUARDL1L2CELLTOFACEHARMONIC 
// Prototype for Fortran procedure L1L2CELLTOFACEHARMONIC ...
//
void FORTRAN_NAME( L1L2CELLTOFACEHARMONIC ,l1l2celltofaceharmonic )(
      CHFp_FRA1(face)
      ,CHFp_CONST_FRA1(cell)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(tiny)
      ,CHFp_BOX(fbox) );

#define FORT_L1L2CELLTOFACEHARMONIC FORTRAN_NAME( inlineL1L2CELLTOFACEHARMONIC, inlineL1L2CELLTOFACEHARMONIC)
#define FORTNT_L1L2CELLTOFACEHARMONIC FORTRAN_NAME( L1L2CELLTOFACEHARMONIC, l1l2celltofaceharmonic)

inline void FORTRAN_NAME(inlineL1L2CELLTOFACEHARMONIC, inlineL1L2CELLTOFACEHARMONIC)(
      CHFp_FRA1(face)
      ,CHFp_CONST_FRA1(cell)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(tiny)
      ,CHFp_BOX(fbox) )
{
 CH_TIMELEAF("FORT_L1L2CELLTOFACEHARMONIC");
 FORTRAN_NAME( L1L2CELLTOFACEHARMONIC ,l1l2celltofaceharmonic )(
      CHFt_FRA1(face)
      ,CHFt_CONST_FRA1(cell)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(tiny)
      ,CHFt_BOX(fbox) );
}
#endif  // GUARDL1L2CELLTOFACEHARMONIC 

}

#endif
