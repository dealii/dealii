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
  solution[level] = 0.;
  
  if (level == minlevel)
    {
      (*coarse)(level, solution[level], defect[level]);
      return;
    }

//  deallog << "Pre-smooth " << level << endl;
  
				   // smoothing of the residual by
				   // modifying s
  pre_smooth->smooth(level, solution[level], defect[level]);

//  deallog << "vmult " << level << endl;
				   // t = A*solution[level]
  matrix->vmult(level, t[level], solution[level]);

//  deallog << "restrict " << level << endl;
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

//    deallog << "recursion " << level << endl;
				   // do recursion
  level_v_step(level-1);

				   // reset size of the auxiliary
				   // vector, since it has been
				   // resized in the recursive call to
				   // level_v_step directly above
  t[level] = 0.;

				   // do coarse grid correction
//  deallog << "prolongate " << level << endl;
  transfer->prolongate(level, t[level], solution[level-1]);
  
  solution[level] += t[level];

  if (edge_up != 0)
    {
      edge_up->Tvmult(level, t[level], solution[level-1]);
      defect[level] -= t[level];
    }
//  deallog << "Post-smooth " << level << endl;
				   // post-smoothing
  post_smooth->smooth(level, solution[level], defect[level]);
//  deallog << "ready " << level << endl;
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
//  abort ();
}




#endif
