#ifndef _AVERAGEHO_F_H_
#define _AVERAGEHO_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGEHO 
#define GUARDAVERAGEHO 
// Prototype for Fortran procedure AVERAGEHO ...
//
void FORTRAN_NAME( AVERAGEHO ,averageho )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_FRA1(lap)
      ,CHFp_BOX(b)
      ,CHFp_BOX(lapBox)
      ,CHFp_CONST_INT(nref)
      ,CHFp_BOX(fineRefBox)
      ,CHFp_INT(doHO)
      ,CHFp_INT(doAverage) );

#define FORT_AVERAGEHO FORTRAN_NAME( inlineAVERAGEHO, inlineAVERAGEHO)
#define FORTNT_AVERAGEHO FORTRAN_NAME( AVERAGEHO, averageho)

inline void FORTRAN_NAME(inlineAVERAGEHO, inlineAVERAGEHO)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_FRA1(lap)
      ,CHFp_BOX(b)
      ,CHFp_BOX(lapBox)
      ,CHFp_CONST_INT(nref)
      ,CHFp_BOX(fineRefBox)
      ,CHFp_INT(doHO)
      ,CHFp_INT(doAverage) )
{
 CH_TIMELEAF("FORT_AVERAGEHO");
 FORTRAN_NAME( AVERAGEHO ,averageho )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_FRA1(lap)
      ,CHFt_BOX(b)
      ,CHFt_BOX(lapBox)
      ,CHFt_CONST_INT(nref)
      ,CHFt_BOX(fineRefBox)
      ,CHFt_INT(doHO)
      ,CHFt_INT(doAverage) );
}
#endif  // GUARDAVERAGEHO 

}

#endif
