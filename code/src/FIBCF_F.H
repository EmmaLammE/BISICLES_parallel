#ifndef _FIBCF_F_H_
#define _FIBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFILLHOLES 
#define GUARDFILLHOLES 
// Prototype for Fortran procedure FILLHOLES ...
//
void FORTRAN_NAME( FILLHOLES ,fillholes )(
      CHFp_FRA1(data)
      ,CHFp_REAL(holeVal)
      ,CHFp_REAL(eps)
      ,CHFp_INT(numNeigbor)
      ,CHFp_BOX(box) );

#define FORT_FILLHOLES FORTRAN_NAME( inlineFILLHOLES, inlineFILLHOLES)
#define FORTNT_FILLHOLES FORTRAN_NAME( FILLHOLES, fillholes)

inline void FORTRAN_NAME(inlineFILLHOLES, inlineFILLHOLES)(
      CHFp_FRA1(data)
      ,CHFp_REAL(holeVal)
      ,CHFp_REAL(eps)
      ,CHFp_INT(numNeigbor)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FILLHOLES");
 FORTRAN_NAME( FILLHOLES ,fillholes )(
      CHFt_FRA1(data)
      ,CHFt_REAL(holeVal)
      ,CHFt_REAL(eps)
      ,CHFt_INT(numNeigbor)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFILLHOLES 

#ifndef GUARDNODETOCELL 
#define GUARDNODETOCELL 
// Prototype for Fortran procedure NODETOCELL ...
//
void FORTRAN_NAME( NODETOCELL ,nodetocell )(
      CHFp_CONST_FRA1(nodeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox) );

#define FORT_NODETOCELL FORTRAN_NAME( inlineNODETOCELL, inlineNODETOCELL)
#define FORTNT_NODETOCELL FORTRAN_NAME( NODETOCELL, nodetocell)

inline void FORTRAN_NAME(inlineNODETOCELL, inlineNODETOCELL)(
      CHFp_CONST_FRA1(nodeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_NODETOCELL");
 FORTRAN_NAME( NODETOCELL ,nodetocell )(
      CHFt_CONST_FRA1(nodeData)
      ,CHFt_FRA1(cellData)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDNODETOCELL 

#ifndef GUARDNODETOCELLCISMVEL 
#define GUARDNODETOCELLCISMVEL 
// Prototype for Fortran procedure NODETOCELLCISMVEL ...
//
void FORTRAN_NAME( NODETOCELLCISMVEL ,nodetocellcismvel )(
      CHFp_CONST_FRA(nodeData)
      ,CHFp_FRA(cellData)
      ,CHFp_CONST_INT(nGhost)
      ,CHFp_CONST_INT(offset)
      ,CHFp_BOX(cellBox) );

#define FORT_NODETOCELLCISMVEL FORTRAN_NAME( inlineNODETOCELLCISMVEL, inlineNODETOCELLCISMVEL)
#define FORTNT_NODETOCELLCISMVEL FORTRAN_NAME( NODETOCELLCISMVEL, nodetocellcismvel)

inline void FORTRAN_NAME(inlineNODETOCELLCISMVEL, inlineNODETOCELLCISMVEL)(
      CHFp_CONST_FRA(nodeData)
      ,CHFp_FRA(cellData)
      ,CHFp_CONST_INT(nGhost)
      ,CHFp_CONST_INT(offset)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_NODETOCELLCISMVEL");
 FORTRAN_NAME( NODETOCELLCISMVEL ,nodetocellcismvel )(
      CHFt_CONST_FRA(nodeData)
      ,CHFt_FRA(cellData)
      ,CHFt_CONST_INT(nGhost)
      ,CHFt_CONST_INT(offset)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDNODETOCELLCISMVEL 

#ifndef GUARDCELLTONODE 
#define GUARDCELLTONODE 
// Prototype for Fortran procedure CELLTONODE ...
//
void FORTRAN_NAME( CELLTONODE ,celltonode )(
      CHFp_FRA1(nodeData)
      ,CHFp_CONST_FRA1(cellData)
      ,CHFp_BOX(nodeBox) );

#define FORT_CELLTONODE FORTRAN_NAME( inlineCELLTONODE, inlineCELLTONODE)
#define FORTNT_CELLTONODE FORTRAN_NAME( CELLTONODE, celltonode)

inline void FORTRAN_NAME(inlineCELLTONODE, inlineCELLTONODE)(
      CHFp_FRA1(nodeData)
      ,CHFp_CONST_FRA1(cellData)
      ,CHFp_BOX(nodeBox) )
{
 CH_TIMELEAF("FORT_CELLTONODE");
 FORTRAN_NAME( CELLTONODE ,celltonode )(
      CHFt_FRA1(nodeData)
      ,CHFt_CONST_FRA1(cellData)
      ,CHFt_BOX(nodeBox) );
}
#endif  // GUARDCELLTONODE 

#ifndef GUARDCELLTONODECISMVELNOSHEAR 
#define GUARDCELLTONODECISMVELNOSHEAR 
// Prototype for Fortran procedure CELLTONODECISMVELNOSHEAR ...
//
void FORTRAN_NAME( CELLTONODECISMVELNOSHEAR ,celltonodecismvelnoshear )(
      CHFp_FRA(nodeData)
      ,CHFp_CONST_FRA1(cellData)
      ,CHFp_CONST_INT(nGhost)
      ,CHFp_CONST_INT(offset)
      ,CHFp_CONST_INTVECT(boxlovect)
      ,CHFp_BOX(nodeBox) );

#define FORT_CELLTONODECISMVELNOSHEAR FORTRAN_NAME( inlineCELLTONODECISMVELNOSHEAR, inlineCELLTONODECISMVELNOSHEAR)
#define FORTNT_CELLTONODECISMVELNOSHEAR FORTRAN_NAME( CELLTONODECISMVELNOSHEAR, celltonodecismvelnoshear)

inline void FORTRAN_NAME(inlineCELLTONODECISMVELNOSHEAR, inlineCELLTONODECISMVELNOSHEAR)(
      CHFp_FRA(nodeData)
      ,CHFp_CONST_FRA1(cellData)
      ,CHFp_CONST_INT(nGhost)
      ,CHFp_CONST_INT(offset)
      ,CHFp_CONST_INTVECT(boxlovect)
      ,CHFp_BOX(nodeBox) )
{
 CH_TIMELEAF("FORT_CELLTONODECISMVELNOSHEAR");
 FORTRAN_NAME( CELLTONODECISMVELNOSHEAR ,celltonodecismvelnoshear )(
      CHFt_FRA(nodeData)
      ,CHFt_CONST_FRA1(cellData)
      ,CHFt_CONST_INT(nGhost)
      ,CHFt_CONST_INT(offset)
      ,CHFt_CONST_INTVECT(boxlovect)
      ,CHFt_BOX(nodeBox) );
}
#endif  // GUARDCELLTONODECISMVELNOSHEAR 

}

#endif
