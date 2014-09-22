// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__multigrid_templates_h
#define __deal2__multigrid_templates_h
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/base/logstream.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN


template <class VECTOR>
Multigrid<VECTOR>::Multigrid (const unsigned int minlevel,
                              const unsigned int maxlevel,
                              const MGMatrixBase<VECTOR> &matrix,
                              const MGCoarseGridBase<VECTOR> &coarse,
                              const MGTransferBase<VECTOR> &transfer,
                              const MGSmootherBase<VECTOR> &pre_smooth,
                              const MGSmootherBase<VECTOR> &post_smooth,
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
  edge_out(0, typeid(*this).name()),
  edge_in(0, typeid(*this).name()),
  edge_down(0, typeid(*this).name()),
  edge_up(0, typeid(*this).name()),
  debug(0)
{}



template <class VECTOR>
void
Multigrid<VECTOR>::reinit(const unsigned int min_level,
                          const unsigned int max_level)
{
  minlevel=min_level;
  maxlevel=max_level;
  defect.resize(minlevel, maxlevel);
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_maxlevel (const unsigned int l)
{
  Assert (l <= maxlevel, ExcIndexRange(l,minlevel,maxlevel+1));
  Assert (l >= minlevel, ExcIndexRange(l,minlevel,maxlevel+1));
  maxlevel = l;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_minlevel (const unsigned int l,
                                 const bool relative)
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
Multigrid<VECTOR>::set_debug (const unsigned int d)
{
  debug = d;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_edge_matrices (const MGMatrixBase<VECTOR> &down,
                                      const MGMatrixBase<VECTOR> &up)
{
  edge_out = &down;
  edge_in = &up;
}


template <class VECTOR>
void
Multigrid<VECTOR>::set_edge_flux_matrices (const MGMatrixBase<VECTOR> &down,
                                           const MGMatrixBase<VECTOR> &up)
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
  if (debug>2)
    deallog << "V-cycle  Defect norm   " << defect[level].l2_norm()
            << std::endl;

  if (level == minlevel)
    {
      if (debug>0)
        deallog << "Coarse level           " << level << std::endl;
      (*coarse)(level, solution[level], defect[level]);
      return;
    }
  if (debug>1)
    deallog << "Smoothing on     level " << level << std::endl;
  // smoothing of the residual by
  // modifying s
//  defect[level].print(std::cout, 2,false);
//  std::cout<<std::endl;
  pre_smooth->smooth(level, solution[level], defect[level]);
//  solution[level].print(std::cout, 2,false);

  if (debug>2)
    deallog << "Solution norm          " << solution[level].l2_norm()
            << std::endl;

  if (debug>1)
    deallog << "Residual on      level " << level << std::endl;
  // t = A*solution[level]
  matrix->vmult(level, t[level], solution[level]);

  if (debug>2)
    deallog << "Residual norm          " << t[level].l2_norm()
            << std::endl;
//  std::cout<<std::endl;
//  t[level].print(std::cout, 2,false);

  // make t rhs of lower level The
  // non-refined parts of the
  // coarse-level defect already
  // contain the global defect, the
  // refined parts its restriction.
  for (unsigned int l = level; l>minlevel; --l)
    {
      t[l-1] = 0.;
      if (l==level && edge_out != 0)
        {
          edge_out->vmult_add(level, t[level], solution[level]);
          if (debug>2)
            deallog << "Norm     t[" << level << "] " << t[level].l2_norm() << std::endl;
        }

      if (l==level && edge_down != 0)
        edge_down->vmult(level, t[level-1], solution[level]);

      transfer->restrict_and_add (l, t[l-1], t[l]);
      if (debug>3)
        deallog << "restrict t[" << l-1 << "] " << t[l-1].l2_norm() << std::endl;
      defect[l-1] -= t[l-1];
      if (debug>3)
        deallog << "defect   d[" << l-1 << "] " << defect[l-1].l2_norm() << std::endl;
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
  if (debug>2)
    deallog << "Prolongate norm        " << t[level].l2_norm() << std::endl;

  solution[level] += t[level];

  if (edge_in != 0)
    {
      edge_in->Tvmult(level, t[level], solution[level]);
      defect[level] -= t[level];
    }

  if (edge_up != 0)
    {
      edge_up->Tvmult(level, t[level], solution[level-1]);
      defect[level] -= t[level];
    }

  if (debug>2)
    deallog << "V-cycle  Defect norm   " << defect[level].l2_norm()
            << std::endl;

  if (debug>1)
    deallog << "Smoothing on     level " << level << std::endl;
  // post-smoothing

//  std::cout<<std::endl;
//  defect[level].print(std::cout, 2,false);
  post_smooth->smooth(level, solution[level], defect[level]);
//  solution[level].print(std::cout, 2,false);
//  std::cout<<std::endl;

  if (debug>2)
    deallog << "Solution norm          " << solution[level].l2_norm()
            << std::endl;

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
    case v_cycle:
      cychar = 'V';
      break;
    case f_cycle:
      cychar = 'F';
      break;
    case w_cycle:
      cychar = 'W';
      break;
    }

  if (debug>0)
    deallog << cychar << "-cycle entering level " << level << std::endl;

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

  if (debug>2)
    deallog << cychar << "-cycle defect norm    " << t[level].l2_norm()
            << std::endl;

  if (level == minlevel)
    {
      if (debug>0)
        deallog << cychar << "-cycle coarse level   " << level << std::endl;

      (*coarse)(level, solution[level], t[level]);
      return;
    }
  if (debug>1)
    deallog << cychar << "-cycle smoothing level " << level << std::endl;
  // smoothing of the residual by
  // modifying s
  pre_smooth->smooth(level, solution[level], t[level]);

  if (debug>2)
    deallog << cychar << "-cycle solution norm    "
            << solution[level].l2_norm() << std::endl;

  if (debug>1)
    deallog << cychar << "-cycle residual level   " << level << std::endl;
  // t = A*solution[level]
  matrix->vmult(level, t[level], solution[level]);
  // make t rhs of lower level The
  // non-refined parts of the
  // coarse-level defect already
  // contain the global defect, the
  // refined parts its restriction.
  if (edge_out != 0)
    edge_out->vmult_add(level, t[level], solution[level]);

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


  if (edge_in != 0)
    edge_in->Tvmult(level, t[level], solution[level]);

  if (edge_up != 0)
    edge_up->Tvmult(level, t[level], solution[level-1]);

  t[level].sadd(-1.,-1.,defect2[level],1.,defect[level]);

  if (debug>2)
    deallog << cychar << "-cycle  Defect norm    " << t[level].l2_norm()
            << std::endl;

  if (debug>1)
    deallog << cychar << "-cycle smoothing level " << level << std::endl;
  // post-smoothing
  post_smooth->smooth(level, solution[level], t[level]);

  if (debug>1)
    deallog << cychar << "-cycle leaving level   " << level << std::endl;
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

  for (unsigned int level=minlevel; level<=maxlevel; ++level)
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

  for (unsigned int level=minlevel; level<=maxlevel; ++level)
    {
      solution[level].reinit(defect[level]);
      t[level].reinit(defect[level]);
    }
  level_v_step (maxlevel);
}


template <class VECTOR>
void
Multigrid<VECTOR>::vmult(VECTOR &dst, const VECTOR &src) const
{
  Multigrid<VECTOR> &mg = const_cast<Multigrid<VECTOR>&>(*this);
  mg.defect[maxlevel] = src;
  for (unsigned int level=maxlevel; level>minlevel; --level)
    {
      mg.defect[level-1] = 0.;
      mg.transfer->restrict_and_add (level,
                                     mg.defect[level-1],
                                     mg.defect[level]);
    }

  mg.cycle();
  dst = mg.solution[maxlevel];
}


template <class VECTOR>
void
Multigrid<VECTOR>::vmult_add(VECTOR &dst, const VECTOR &src) const
{
  Multigrid<VECTOR> &mg = const_cast<Multigrid<VECTOR>&>(*this);
  mg.defect[maxlevel] = src;
  mg.cycle();
  dst += mg.solution[maxlevel];
}


template <class VECTOR>
void
Multigrid<VECTOR>::Tvmult(VECTOR &, const VECTOR &) const
{
  Assert(false, ExcNotImplemented());
}


template <class VECTOR>
void
Multigrid<VECTOR>::Tvmult_add(VECTOR &, const VECTOR &) const
{
  Assert(false, ExcNotImplemented());
}

DEAL_II_NAMESPACE_CLOSE

#endif
