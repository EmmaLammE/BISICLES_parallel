#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SundialsUtil.H"
#include "FineInterp.H"
#include "PiecewiseLinearFillPatch.H"

#include "NamespaceHeader.H"


/// create a new LevelData<FArrayBox>*, interpolate from a_crseData and copy from a_oldData as needed, delete a_oldData 
LevelData<FArrayBox>* destructiveRegrid(LevelData<FArrayBox>* a_oldData,
					const DisjointBoxLayout& a_newDBL,
				        const LevelData<FArrayBox>* a_crseData,
					int a_ratio)
{
  //CH_assert(a_oldData);
  CH_assert(a_crseData);
  LevelData<FArrayBox>* newData = new LevelData<FArrayBox>(a_newDBL, a_crseData->nComp(), a_crseData->ghostVect());
  CH_assert(newData);
  
  if (a_crseData)
    {
      CH_assert(a_crseData->nComp() == newData->nComp());
      CH_assert(a_crseData->ghostVect() == newData->ghostVect() );
  
      FineInterp interpolator(a_newDBL, newData->nComp(), a_ratio, newData->disjointBoxLayout().physDomain());
      interpolator.interpToFine(*newData, *a_crseData);
      
      PiecewiseLinearFillPatch ghostFiller
	(a_newDBL, a_crseData->disjointBoxLayout() ,  a_crseData->nComp(),
	 a_crseData->disjointBoxLayout().physDomain(), a_ratio, a_crseData->ghostVect()[0]);

      ghostFiller.fillInterp(*newData, *a_crseData, *a_crseData, 1.0, 0, 0,  a_crseData->nComp());
 
    }
  
  if (a_oldData)
    {
      if (a_oldData->isDefined())
	{
	  a_oldData->copyTo(*newData);
	}
      delete a_oldData;
    }
  
  newData->exchange();
  return newData;
  
}



/// reshape dest and fill with data from src
void reshapeAndFill(Vector<LevelData<FArrayBox>* >& a_dest,
                    const Vector<LevelData<FArrayBox>* >& a_src)
{

  // if not defined yet, set the number of levels here
  if (a_dest.size() == 0)
    {
      a_dest.resize(a_src.size());
    }

  // if a_dest was defined, it needs to be at least as big as src
  CH_assert(a_dest.size() >= a_src.size());
  // if a_dest size > a_src size, do we need to resize, or can we just
  // delete finer levels (assume the latter for now)
  
  for (int lev=0; lev<a_src.size(); lev++)
    {
      // eventually can be clever and not reshape level 0, just do a copy
      // (since level 0 grids never change)
      // if we need to, allocate a LevelData here
      if (a_dest[lev] == NULL)
        {
          a_dest[lev] = new LevelData<FArrayBox>;
        }
      
      reshapeAndFill(*a_dest[lev], *a_src[lev]);
    }

  // if there are any extra levels in dest, delete them
  for (int lev=a_src.size(); lev<a_dest.size(); lev++)
    {
      delete a_dest[lev];
      a_dest[lev] = NULL;
    }
}



/// reshape dest and fill with data from src
void reshapeAndFill(LevelData<FArrayBox>& a_dest,
                    const LevelData<FArrayBox>& a_src)
{
  const DisjointBoxLayout& grids = a_src.getBoxes();
  IntVect ghostVect = a_src.ghostVect();
  int nComp = a_src.nComp();

  a_dest.define(grids, nComp, ghostVect);
  // now do a fab-by-fab copy to get any ghost cells as well
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_dest[dit].copy(a_src[dit]);
    }

}



/// reshape dest to match src
void reshape(Vector<LevelData<FArrayBox>* >& a_dest,
             const Vector<LevelData<FArrayBox>* >& a_src)
{

  // if not defined yet, set the number of levels here
  if (a_dest.size() == 0)
    {
      a_dest.resize(a_src.size());
    }

  // if a_dest was defined, it needs to be at least as big as src
  CH_assert(a_dest.size() >= a_src.size()); // check if dest size >= src size
  // if a_dest size > a_src size, do we need to resize, or can we just
  // delete finer levels (assume the latter for now)
  // cout<<"in reshape in sundial, asrc size "<<a_src.size()<<", adest size "<<a_dest.size()<<endl;
  for (int lev=0; lev<a_src.size(); lev++)
    {
      // eventually can be clever and not reshape level 0, just do a copy
      // (since level 0 grids never change)
      // if we need to, allocate a LevelData here
      if (a_dest[lev] == NULL)
        {
          a_dest[lev] = new LevelData<FArrayBox>;
        }
      
      reshape(*a_dest[lev], *a_src[lev]);
    }

  // if there are any extra levels in dest, delete them
  for (int lev=a_src.size(); lev<a_dest.size(); lev++)
    {
      delete a_dest[lev];
      a_dest[lev] = NULL;
    }
}



/// reshape dest to match src
void reshape(LevelData<FArrayBox>& a_dest,
             const LevelData<FArrayBox>& a_src)
{
  const DisjointBoxLayout& grids = a_src.getBoxes();
  IntVect ghostVect = a_src.ghostVect();
  int nComp = a_src.nComp();
  // cout<<"  a_src.getBoxes "<<a_src.getBoxes()<<endl;
  a_dest.define(grids, nComp, ghostVect);
  // cout<<", a_dest box "<<a_dest.getBoxes()<<endl;

}


#include "NamespaceFooter.H"
