//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998 - 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__multigrid_templates_h
#define __deal2__multigrid_templates_h
#include <multigrid/multigrid.h>

#include <base/logstream.h>

#include "iostream"

using namespace std;









template <class VECTOR>
Multigrid<VECTOR>::Multigrid (const unsigned int minlevel,
			      const unsigned int maxlevel,
			      const MGMatrixBase<VECTOR>& matrix,
			      const MGCoarseGridBase<VECTOR>& coarse,
			      const MGTransferBase<VECTOR>& transfer,
			      const MGSmootherBase<VECTOR>& pre_smooth,
			      const MGSmootherBase<VECTOR>& post_smooth,
			      typename Multigrid<VECTOR>::Cycle cycle)
		:
		cycle_type(cycle),
		minlevel(minlevel),
		maxlevel(maxlevel),
		defect(minlevel,maxlevel),
		solution(minlevel,maxlevel),
		t(minlevel,maxlevel),
		matrix(&matrix, typeid(*this).name()),
		coarse(&coarse, typeid(*this).name()),
		transfer(&transfer, typeid(*this).name()),
		pre_smooth(&pre_smooth, typeid(*this).name()),
		post_smooth(&post_smooth, typeid(*this).name()),
		edge_down(0, typeid(*this).name()),
		edge_up(0, typeid(*this).name()),
		debug(0)
{}



//TODO: This function cannot work. At least maxlevel should be tested!
template <class VECTOR>
void
Multigrid<VECTOR>::reinit(const unsigned int min_level,
			  const unsigned int max_level)
{
  Assert(false, ExcNotImplemented());
  
  minlevel=min_level;
  maxlevel=max_level;
  defect.resize(minlevel, maxlevel);
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_maxlevel(unsigned int l)
{
  Assert (l <= maxlevel, ExcIndexRange(l,minlevel,maxlevel+1));
  Assert (l >= minlevel, ExcIndexRange(l,minlevel,maxlevel+1));
  maxlevel = l;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_minlevel(unsigned int l,
				bool relative)
{
  Assert (l <= maxlevel, ExcIndexRange(l,minlevel,maxlevel+1));
  minlevel = (relative)
	     ? (maxlevel-l)
	     : l;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_cycle(typename Multigrid<VECTOR>::Cycle c)
{
  cycle_type = c;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_debug(unsigned int d)
{
  debug = d;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_edge_matrices (const MGMatrixBase<VECTOR>& down,
				      const MGMatrixBase<VECTOR>& up)
{
  edge_down = &down;
  edge_up = &up;
}




template <class VECTOR>
void
Multigrid<VECTOR>::level_v_step(const unsigned int level)
{
  if (debug>0)
    deallog << "V-cycle entering level " << level << std::endl;
  
  if (level == minlevel)
    {
      (*coarse)(level, solution[level], defect[level]);
      return;
    }
  if (debug>1)
    deallog << "Smoothing on     level " << level << std::endl;
				   // smoothing of the residual by
				   // modifying s
  pre_smooth->smooth(level, solution[level], defect[level]);

  if (debug>1)
    deallog << "Residual on      level " << level << std::endl;
				   // t = A*solution[level]
  matrix->vmult(level, t[level], solution[level]);

				   // make t rhs of lower level The
				   // non-refined parts of the
				   // coarse-level defect already
				   // contain the global defect, the
				   // refined parts its restriction.
  for (unsigned int l = level;l>minlevel;--l)
    {
      t[l-1] = 0.;
      if (l==level && edge_down != 0)
	{
	  edge_down->vmult(level, t[level-1], solution[level]);
	}
      transfer->restrict_and_add (l, t[l-1], t[l]);
      defect[l-1] -= t[l-1];
    }

				   // do recursion
  solution[level-1] = 0.;
  level_v_step(level-1);

				   // reset size of the auxiliary
				   // vector, since it has been
				   // resized in the recursive call to
				   // level_v_step directly above
  t[level] = 0.;

				   // do coarse grid correction
  transfer->prolongate(level, t[level], solution[level-1]);
  
  solution[level] += t[level];

  if (edge_up != 0)
    {
      edge_up->Tvmult(level, t[level], solution[level-1]);
      defect[level] -= t[level];
    }
  
  if (debug>1)
    deallog << "Smoothing on     level " << level << std::endl;
				   // post-smoothing
  post_smooth->smooth(level, solution[level], defect[level]);

  if (debug>1)
    deallog << "V-cycle leaving  level " << level << std::endl;
}



template <class VECTOR>
void
Multigrid<VECTOR>::level_step(const unsigned int level,
			      Cycle cycle)
{
  char cychar = '?';
  switch (cycle)
    {
      case v_cycle: cychar = 'V'; break;
      case f_cycle: cychar = 'F'; break;
      case w_cycle: cychar = 'W'; break;
    }
  
				   // Not actually the defect yet, but
				   // we do not want to spend yet
				   // another vector.
  if (level>minlevel)
    {
      defect2[level-1] = 0.;
      transfer->restrict_and_add (level, defect2[level-1], defect2[level]);
    }
  
				   // We get an update of the defect
				   // from the previous level in t and
				   // from two levels above in
				   // defect2. This is subtracted from
				   // the original defect.
  t[level].equ(-1.,defect2[level],1.,defect[level]);
  
  if (level == minlevel)
    {
      if (debug>0)
	deallog << cychar << "-cycle solving  level " << level << std::endl;
  
      (*coarse)(level, solution[level], t[level]);
      return;
    }
  if (debug>0)
    deallog << cychar << "-cycle entering level " << level << std::endl;
  
  if (debug>1)
    deallog << "Smoothing on     level " << level << std::endl;
				   // smoothing of the residual by
				   // modifying s
  pre_smooth->smooth(level, solution[level], t[level]);

  if (debug>1)
    deallog << "Residual on      level " << level << std::endl;
				   // t = A*solution[level]
  matrix->vmult(level, t[level], solution[level]);
				   // make t rhs of lower level The
				   // non-refined parts of the
				   // coarse-level defect already
				   // contain the global defect, the
				   // refined parts its restriction.
  if (edge_down != 0)
    edge_down->vmult_add(level, defect2[level-1], solution[level]);
  transfer->restrict_and_add (level, defect2[level-1], t[level]);
				   // do recursion
  solution[level-1] = 0.;
				   // Every cycle need one recursion,
				   // the V-cycle, which is included
				   // here for the sake of the
				   // F-cycle, needs only one,
  level_step(level-1, cycle);
				   // If we solve exactly, then we do
				   // not need a second coarse grid
				   // step.
  if (level>minlevel+1)
    {
				       // while the W-cycle repeats itself
      if (cycle == w_cycle)
	level_step(level-1, cycle);
				       // and the F-cycle does a V-cycle
				       // after an F-cycle.
      else if (cycle == f_cycle)
	level_step(level-1, v_cycle);
    }
  
				   // reset size of the auxiliary
				   // vector, since it has been
				   // resized in the recursive call to
				   // level_v_step directly above
  t[level] = 0.;
				   // do coarse grid correction
  transfer->prolongate(level, t[level], solution[level-1]);
  
  solution[level] += t[level];

  
  if (edge_up != 0)
    {
      edge_up->Tvmult(level, t[level], solution[level-1]);
    }

  t[level].sadd(-1.,-1.,defect2[level],1.,defect[level]);

  if (debug>1)
    deallog << "Smoothing on     level " << level << std::endl;
				   // post-smoothing
  post_smooth->smooth(level, solution[level], t[level]);

  if (debug>1)
    deallog << cychar << "-cycle leaving  level " << level << std::endl;
}


template <class VECTOR>
void
Multigrid<VECTOR>::cycle()
{
				   // The defect vector has been
				   // initialized by copy_to_mg. Now
				   // adjust the other vectors.
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);
  if (cycle_type != v_cycle)
    defect2.resize(minlevel, maxlevel);
    
  for (unsigned int level=minlevel; level<=maxlevel;++level)
    {
      solution[level].reinit(defect[level]);
      t[level].reinit(defect[level]);
      if (cycle_type != v_cycle)
	defect2[level].reinit(defect[level]);
    }
  
  if (cycle_type == v_cycle)
    level_v_step (maxlevel);
  else
    level_step (maxlevel, cycle_type);
}


template <class VECTOR>
void
Multigrid<VECTOR>::vcycle()
{
				   // The defect vector has been
				   // initialized by copy_to_mg. Now
				   // adjust the other vectors.
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);
  
  for (unsigned int level=minlevel; level<=maxlevel;++level)
    {
      solution[level].reinit(defect[level]);
      t[level].reinit(defect[level]);
    }
  level_v_step (maxlevel);
}




#endif
