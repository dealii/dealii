// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_multigrid_templates_h
#define dealii_multigrid_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>

#include <deal.II/multigrid/multigrid.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN

template <typename VectorType>
void
Multigrid<VectorType>::reinit(const unsigned int min_level,
                              const unsigned int max_level)
{
  Assert(min_level >= matrix->get_minlevel(),
         ExcLowerRangeType<unsigned int>(min_level, matrix->get_minlevel()));
  Assert(max_level <= matrix->get_maxlevel(),
         ExcLowerRangeType<unsigned int>(matrix->get_maxlevel(), max_level));
  Assert(min_level <= max_level,
         ExcLowerRangeType<unsigned int>(max_level, min_level));
  minlevel = min_level;
  maxlevel = max_level;
  // solution, t and defect2 are resized in cycle()
  defect.resize(minlevel, maxlevel);
}



template <typename VectorType>
void
Multigrid<VectorType>::set_maxlevel(const unsigned int l)
{
  reinit(minlevel, l);
}



template <typename VectorType>
void
Multigrid<VectorType>::set_minlevel(const unsigned int l, const bool relative)
{
  const unsigned int new_minlevel = (relative) ? (maxlevel - l) : l;
  reinit(new_minlevel, maxlevel);
}



template <typename VectorType>
void
Multigrid<VectorType>::set_cycle(typename Multigrid<VectorType>::Cycle c)
{
  cycle_type = c;
}



template <typename VectorType>
void
Multigrid<VectorType>::set_edge_matrices(const MGMatrixBase<VectorType> &down,
                                         const MGMatrixBase<VectorType> &up)
{
  edge_out = &down;
  edge_in  = &up;
}



template <typename VectorType>
void
Multigrid<VectorType>::set_edge_flux_matrices(
  const MGMatrixBase<VectorType> &down,
  const MGMatrixBase<VectorType> &up)
{
  edge_down = &down;
  edge_up   = &up;
}



template <typename VectorType>
void
Multigrid<VectorType>::level_v_step(const unsigned int level)
{
  if (level == minlevel)
    {
      this->signals.coarse_solve(true, level);
      (*coarse)(level, solution[level], defect[level]);
      this->signals.coarse_solve(false, level);
      return;
    }

  // smoothing of the residual
  this->signals.pre_smoother_step(true, level);
  pre_smooth->apply(level, solution[level], defect[level]);
  this->signals.pre_smoother_step(false, level);

  // compute residual on level, which includes the (CG) edge matrix
  matrix->vmult(level, t[level], solution[level]);
  if (edge_out != nullptr)
    {
      edge_out->vmult_add(level, t[level], solution[level]);
    }
  t[level].sadd(-1.0, 1.0, defect[level]);

  // Get the defect on the next coarser level as part of the (DG) edge matrix
  // and then the main part by the restriction of the transfer
  if (edge_down != nullptr)
    {
      edge_down->vmult(level, t[level - 1], solution[level]);
      defect[level - 1] -= t[level - 1];
    }

  this->signals.restriction(true, level);
  transfer->restrict_and_add(level, defect[level - 1], t[level]);
  this->signals.restriction(false, level);

  // do recursion
  level_v_step(level - 1);

  // do coarse grid correction
  this->signals.prolongation(true, level);
  transfer->prolongate(level, t[level], solution[level - 1]);
  this->signals.prolongation(false, level);

  solution[level] += t[level];

  // get in contribution from edge matrices to the defect
  if (edge_in != nullptr)
    {
      edge_in->Tvmult(level, t[level], solution[level]);
      defect[level] -= t[level];
    }
  if (edge_up != nullptr)
    {
      edge_up->Tvmult(level, t[level], solution[level - 1]);
      defect[level] -= t[level];
    }

  // post-smoothing
  this->signals.post_smoother_step(true, level);
  post_smooth->smooth(level, solution[level], defect[level]);
  this->signals.post_smoother_step(false, level);
}



template <typename VectorType>
void
Multigrid<VectorType>::level_step(const unsigned int level, Cycle cycle)
{
  // Combine the defect from the initial copy_to_mg with the one that has come
  // from the finer level by the transfer
  defect2[level] += defect[level];
  defect[level] = typename VectorType::value_type(0.);

  if (level == minlevel)
    {
      this->signals.coarse_solve(true, level);
      (*coarse)(level, solution[level], defect2[level]);
      this->signals.coarse_solve(false, level);
      return;
    }

  // smoothing of the residual
  this->signals.pre_smoother_step(true, level);
  pre_smooth->apply(level, solution[level], defect2[level]);
  this->signals.pre_smoother_step(false, level);

  // compute residual on level, which includes the (CG) edge matrix
  matrix->vmult(level, t[level], solution[level]);
  if (edge_out != nullptr)
    edge_out->vmult_add(level, t[level], solution[level]);
  t[level].sadd(-1.0, 1.0, defect2[level]);

  // Get the defect on the next coarser level as part of the (DG) edge matrix
  // and then the main part by the restriction of the transfer
  if (edge_down != nullptr)
    edge_down->vmult(level, defect2[level - 1], solution[level]);
  else
    defect2[level - 1] = typename VectorType::value_type(0.);

  this->signals.restriction(true, level);
  transfer->restrict_and_add(level, defect2[level - 1], t[level]);
  this->signals.restriction(false, level);

  // Every cycle starts with a recursion of its type.
  level_step(level - 1, cycle);

  // For W and F-cycle, repeat the process on the next coarser level except
  // for the coarse solver which we invoke just once
  if (level > minlevel + 1)
    {
      // while the W-cycle repeats itself, ...
      if (cycle == w_cycle)
        level_step(level - 1, cycle);
      // ... the F-cycle does a V-cycle after an F-cycle, ...
      else if (cycle == f_cycle)
        level_step(level - 1, v_cycle);
      // ... and the V-cycle does nothing.
    }

  // do coarse grid correction
  this->signals.prolongation(true, level);
  transfer->prolongate(level, t[level], solution[level - 1]);
  this->signals.prolongation(false, level);
  solution[level] += t[level];

  // get in contribution from edge matrices to the defect
  if (edge_in != nullptr)
    {
      edge_in->Tvmult(level, t[level], solution[level]);
      defect2[level] -= t[level];
    }

  if (edge_up != nullptr)
    {
      edge_up->Tvmult(level, t[level], solution[level - 1]);
      defect2[level] -= t[level];
    }

  // post-smoothing
  this->signals.post_smoother_step(true, level);
  post_smooth->smooth(level, solution[level], defect2[level]);
  this->signals.post_smoother_step(false, level);
}



template <typename VectorType>
void
Multigrid<VectorType>::cycle()
{
  // The defect vector has been initialized by copy_to_mg. Now adjust the
  // other vectors.
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);
  if (cycle_type != v_cycle)
    defect2.resize(minlevel, maxlevel);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      // the vectors for level>minlevel will be overwritten by the apply()
      // method of the smoother -> do not force them to be zeroed out here
      solution[level].reinit(defect[level], level > minlevel);
      t[level].reinit(defect[level], level > minlevel);
      if (cycle_type != v_cycle)
        defect2[level].reinit(defect[level]);
    }

  if (cycle_type == v_cycle)
    level_v_step(maxlevel);
  else
    level_step(maxlevel, cycle_type);
}



template <typename VectorType>
void
Multigrid<VectorType>::vcycle()
{
  // The defect vector has been initialized by copy_to_mg. Now adjust the
  // other vectors.
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      solution[level].reinit(defect[level], level > minlevel);
      t[level].reinit(defect[level], level > minlevel);
    }
  level_v_step(maxlevel);
}



template <typename VectorType>
boost::signals2::connection
Multigrid<VectorType>::connect_coarse_solve(
  const std::function<void(const bool, const unsigned int)> &slot)
{
  return this->signals.coarse_solve.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
Multigrid<VectorType>::connect_restriction(
  const std::function<void(const bool, const unsigned int)> &slot)
{
  return this->signals.restriction.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
Multigrid<VectorType>::connect_prolongation(
  const std::function<void(const bool, const unsigned int)> &slot)
{
  return this->signals.prolongation.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
Multigrid<VectorType>::connect_pre_smoother_step(
  const std::function<void(const bool, const unsigned int)> &slot)
{
  return this->signals.pre_smoother_step.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
Multigrid<VectorType>::connect_post_smoother_step(
  const std::function<void(const bool, const unsigned int)> &slot)
{
  return this->signals.post_smoother_step.connect(slot);
}


DEAL_II_NAMESPACE_CLOSE

#endif
