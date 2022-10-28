#ifndef _CONSTITUTIVERELATIONF_F_H_
#define _CONSTITUTIVERELATIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSTRAININVARSSA 
#define GUARDSTRAININVARSSA 
// Prototype for Fortran procedure STRAININVARSSA ...
//
void FORTRAN_NAME( STRAININVARSSA ,straininvarssa )(
      CHFp_FRA1(epsSqr)
      ,CHFp_FRA(derivs)
      ,CHFp_INT(dudxComp)
      ,CHFp_INT(dudyComp)
      ,CHFp_INT(dvdxComp)
      ,CHFp_INT(dvdyComp)
      ,CHFp_BOX(box) );

#define FORT_STRAININVARSSA FORTRAN_NAME( inlineSTRAININVARSSA, inlineSTRAININVARSSA)
#define FORTNT_STRAININVARSSA FORTRAN_NAME( STRAININVARSSA, straininvarssa)

inline void FORTRAN_NAME(inlineSTRAININVARSSA, inlineSTRAININVARSSA)(
      CHFp_FRA1(epsSqr)
      ,CHFp_FRA(derivs)
      ,CHFp_INT(dudxComp)
      ,CHFp_INT(dudyComp)
      ,CHFp_INT(dvdxComp)
      ,CHFp_INT(dvdyComp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_STRAININVARSSA");
 FORTRAN_NAME( STRAININVARSSA ,straininvarssa )(
      CHFt_FRA1(epsSqr)
      ,CHFt_FRA(derivs)
      ,CHFt_INT(dudxComp)
      ,CHFt_INT(dudyComp)
      ,CHFt_INT(dvdxComp)
      ,CHFt_INT(dvdyComp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSTRAININVARSSA 

#ifndef GUARDCOMPUTEARRHENIUSA 
#define GUARDCOMPUTEARRHENIUSA 
// Prototype for Fortran procedure COMPUTEARRHENIUSA ...
//
void FORTRAN_NAME( COMPUTEARRHENIUSA ,computearrheniusa )(
      CHFp_FRA1(A)
      ,CHFp_CONST_FRA1(theta)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(nParameter)
      ,CHFp_CONST_REAL(mParameter)
      ,CHFp_CONST_REAL(B0)
      ,CHFp_CONST_REAL(thetaR)
      ,CHFp_CONST_REAL(Kexponent)
      ,CHFp_CONST_REAL(Cfactor)
      ,CHFp_CONST_REAL(RgasConst)
      ,CHFp_CONST_REAL(Q) );

#define FORT_COMPUTEARRHENIUSA FORTRAN_NAME( inlineCOMPUTEARRHENIUSA, inlineCOMPUTEARRHENIUSA)
#define FORTNT_COMPUTEARRHENIUSA FORTRAN_NAME( COMPUTEARRHENIUSA, computearrheniusa)

inline void FORTRAN_NAME(inlineCOMPUTEARRHENIUSA, inlineCOMPUTEARRHENIUSA)(
      CHFp_FRA1(A)
      ,CHFp_CONST_FRA1(theta)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(nParameter)
      ,CHFp_CONST_REAL(mParameter)
      ,CHFp_CONST_REAL(B0)
      ,CHFp_CONST_REAL(thetaR)
      ,CHFp_CONST_REAL(Kexponent)
      ,CHFp_CONST_REAL(Cfactor)
      ,CHFp_CONST_REAL(RgasConst)
      ,CHFp_CONST_REAL(Q) )
{
 CH_TIMELEAF("FORT_COMPUTEARRHENIUSA");
 FORTRAN_NAME( COMPUTEARRHENIUSA ,computearrheniusa )(
      CHFt_FRA1(A)
      ,CHFt_CONST_FRA1(theta)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(nParameter)
      ,CHFt_CONST_REAL(mParameter)
      ,CHFt_CONST_REAL(B0)
      ,CHFt_CONST_REAL(thetaR)
      ,CHFt_CONST_REAL(Kexponent)
      ,CHFt_CONST_REAL(Cfactor)
      ,CHFt_CONST_REAL(RgasConst)
      ,CHFt_CONST_REAL(Q) );
}
#endif  // GUARDCOMPUTEARRHENIUSA 

#ifndef GUARDCOMPUTEPATERSONA 
#define GUARDCOMPUTEPATERSONA 
// Prototype for Fortran procedure COMPUTEPATERSONA ...
//
void FORTRAN_NAME( COMPUTEPATERSONA ,computepatersona )(
      CHFp_FRA1(A)
      ,CHFp_CONST_FRA1(theta)
      ,CHFp_CONST_FRA1(theta0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(E)
      ,CHFp_CONST_REAL(A0)
      ,CHFp_CONST_REAL(R)
      ,CHFp_CONST_REAL(Qm)
      ,CHFp_CONST_REAL(Qp) );

#define FORT_COMPUTEPATERSONA FORTRAN_NAME( inlineCOMPUTEPATERSONA, inlineCOMPUTEPATERSONA)
#define FORTNT_COMPUTEPATERSONA FORTRAN_NAME( COMPUTEPATERSONA, computepatersona)

inline void FORTRAN_NAME(inlineCOMPUTEPATERSONA, inlineCOMPUTEPATERSONA)(
      CHFp_FRA1(A)
      ,CHFp_CONST_FRA1(theta)
      ,CHFp_CONST_FRA1(theta0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(E)
      ,CHFp_CONST_REAL(A0)
      ,CHFp_CONST_REAL(R)
      ,CHFp_CONST_REAL(Qm)
      ,CHFp_CONST_REAL(Qp) )
{
 CH_TIMELEAF("FORT_COMPUTEPATERSONA");
 FORTRAN_NAME( COMPUTEPATERSONA ,computepatersona )(
      CHFt_FRA1(A)
      ,CHFt_CONST_FRA1(theta)
      ,CHFt_CONST_FRA1(theta0)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(E)
      ,CHFt_CONST_REAL(A0)
      ,CHFt_CONST_REAL(R)
      ,CHFt_CONST_REAL(Qm)
      ,CHFt_CONST_REAL(Qp) );
}
#endif  // GUARDCOMPUTEPATERSONA 

#ifndef GUARDCOMPUTEZWINGERA 
#define GUARDCOMPUTEZWINGERA 
// Prototype for Fortran procedure COMPUTEZWINGERA ...
//
void FORTRAN_NAME( COMPUTEZWINGERA ,computezwingera )(
      CHFp_FRA1(A)
      ,CHFp_CONST_FRA1(theta)
      ,CHFp_CONST_FRA1(theta0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(E)
      ,CHFp_CONST_REAL(Ap)
      ,CHFp_CONST_REAL(R)
      ,CHFp_CONST_REAL(Qm)
      ,CHFp_CONST_REAL(Qp) );

#define FORT_COMPUTEZWINGERA FORTRAN_NAME( inlineCOMPUTEZWINGERA, inlineCOMPUTEZWINGERA)
#define FORTNT_COMPUTEZWINGERA FORTRAN_NAME( COMPUTEZWINGERA, computezwingera)

inline void FORTRAN_NAME(inlineCOMPUTEZWINGERA, inlineCOMPUTEZWINGERA)(
      CHFp_FRA1(A)
      ,CHFp_CONST_FRA1(theta)
      ,CHFp_CONST_FRA1(theta0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(E)
      ,CHFp_CONST_REAL(Ap)
      ,CHFp_CONST_REAL(R)
      ,CHFp_CONST_REAL(Qm)
      ,CHFp_CONST_REAL(Qp) )
{
 CH_TIMELEAF("FORT_COMPUTEZWINGERA");
 FORTRAN_NAME( COMPUTEZWINGERA ,computezwingera )(
      CHFt_FRA1(A)
      ,CHFt_CONST_FRA1(theta)
      ,CHFt_CONST_FRA1(theta0)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(E)
      ,CHFt_CONST_REAL(Ap)
      ,CHFt_CONST_REAL(R)
      ,CHFt_CONST_REAL(Qm)
      ,CHFt_CONST_REAL(Qp) );
}
#endif  // GUARDCOMPUTEZWINGERA 

#ifndef GUARDCOMPUTEGLENSMU0 
#define GUARDCOMPUTEGLENSMU0 
// Prototype for Fortran procedure COMPUTEGLENSMU0 ...
//
void FORTRAN_NAME( COMPUTEGLENSMU0 ,computeglensmu0 )(
      CHFp_FRA1(mu0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(nParameter) );

#define FORT_COMPUTEGLENSMU0 FORTRAN_NAME( inlineCOMPUTEGLENSMU0, inlineCOMPUTEGLENSMU0)
#define FORTNT_COMPUTEGLENSMU0 FORTRAN_NAME( COMPUTEGLENSMU0, computeglensmu0)

inline void FORTRAN_NAME(inlineCOMPUTEGLENSMU0, inlineCOMPUTEGLENSMU0)(
      CHFp_FRA1(mu0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(nParameter) )
{
 CH_TIMELEAF("FORT_COMPUTEGLENSMU0");
 FORTRAN_NAME( COMPUTEGLENSMU0 ,computeglensmu0 )(
      CHFt_FRA1(mu0)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(nParameter) );
}
#endif  // GUARDCOMPUTEGLENSMU0 

#ifndef GUARDCOMPUTEGLENSMU 
#define GUARDCOMPUTEGLENSMU 
// Prototype for Fortran procedure COMPUTEGLENSMU ...
//
void FORTRAN_NAME( COMPUTEGLENSMU ,computeglensmu )(
      CHFp_FRA1(mu)
      ,CHFp_CONST_FRA1(epsSqr)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(nExponent)
      ,CHFp_CONST_REAL(epsSqr0)
      ,CHFp_CONST_REAL(delta) );

#define FORT_COMPUTEGLENSMU FORTRAN_NAME( inlineCOMPUTEGLENSMU, inlineCOMPUTEGLENSMU)
#define FORTNT_COMPUTEGLENSMU FORTRAN_NAME( COMPUTEGLENSMU, computeglensmu)

inline void FORTRAN_NAME(inlineCOMPUTEGLENSMU, inlineCOMPUTEGLENSMU)(
      CHFp_FRA1(mu)
      ,CHFp_CONST_FRA1(epsSqr)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(nExponent)
      ,CHFp_CONST_REAL(epsSqr0)
      ,CHFp_CONST_REAL(delta) )
{
 CH_TIMELEAF("FORT_COMPUTEGLENSMU");
 FORTRAN_NAME( COMPUTEGLENSMU ,computeglensmu )(
      CHFt_FRA1(mu)
      ,CHFt_CONST_FRA1(epsSqr)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(nExponent)
      ,CHFt_CONST_REAL(epsSqr0)
      ,CHFt_CONST_REAL(delta) );
}
#endif  // GUARDCOMPUTEGLENSMU 

}

#endif
