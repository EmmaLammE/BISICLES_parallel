#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef SUNDIALSUTIL_H
#define SUNDIALSUTIL_H

#include "FArrayBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

// helper functions for Sundials (and Sundials-like) activities


struct IceSheetState
{

Vector<LevelData<FArrayBox>* > ice_thickness;
Vector<LevelData<FArrayBox>* > iceAreaFraction;
Vector<LevelData<FArrayBox>* > ice_velocity;
Vector<LevelData<FArrayBox>* > ice_enthalpy;

Real time;

};


/// create a new LevelData<FArrayBox>*, interpolate from a_crseData and copy from a_oldData as needed, delete a_oldData
/** was a local-only function in AmrIceMesh, but we needed it for the setState function so 
    it got pulled out and put here.
*/
LevelData<FArrayBox>* destructiveRegrid(LevelData<FArrayBox>* a_oldData,
					const DisjointBoxLayout& a_newDBL,
				        const LevelData<FArrayBox>* a_crseData,
                                          int a_ratio);


/// reshape dest and fill with data from src
void reshapeAndFill(Vector<LevelData<FArrayBox>* >& a_dest,
                      const Vector<LevelData<FArrayBox>* >& a_src);


/// reshape dest and fill with data from src
void reshapeAndFill(LevelData<FArrayBox>& a_dest,
                    const LevelData<FArrayBox>& a_src);



/// reshape dest to match src
void reshape(Vector<LevelData<FArrayBox>* >& a_dest,
                      const Vector<LevelData<FArrayBox>* >& a_src);


/// reshape dest to match src
void reshape(LevelData<FArrayBox>& a_dest,
                    const LevelData<FArrayBox>& a_src);





#include "NamespaceFooter.H"

#endif
