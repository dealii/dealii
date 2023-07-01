// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2023 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_selector_h
#define dealii_matrix_free_evaluation_selector_h

#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/shape_info.h>

DEAL_II_NAMESPACE_OPEN



/**
 * This class chooses an appropriate evaluation strategy based on the template
 * parameters and the shape_info variable, providing a short-cut to some
 * internal functions.
 */
template <int dim, int fe_degree, int n_q_points_1d, typename Number>
struct SelectEvaluator
{
  /**
   * Chooses an appropriate evaluation strategy for the evaluate function, i.e.
   * this calls internal::FEEvaluationImpl::evaluate(),
   * internal::FEEvaluationImplCollocation::evaluate() or
   * internal::FEEvaluationImplTransformToCollocation::evaluate() with
   * appropriate template parameters.
   */
  static void
  evaluate(const unsigned int                     n_components,
           const EvaluationFlags::EvaluationFlags evaluation_flag,
           const Number *                         values_dofs,
           FEEvaluationData<dim, Number, false> & eval);

  /**
   * Chooses an appropriate evaluation strategy for the integrate function, i.e.
   * this calls internal::FEEvaluationImpl::integrate(),
   * internal::FEEvaluationImplCollocation::integrate() or
   * internal::FEEvaluationImplTransformToCollocation::integrate() with
   * appropriate template parameters.
   */
  static void
  integrate(const unsigned int                     n_components,
            const EvaluationFlags::EvaluationFlags integration_flag,
            Number *                               values_dofs,
            FEEvaluationData<dim, Number, false> & eval,
            const bool sum_into_values_array = false);
};

//----------------------Implementation for SelectEvaluator---------------------
#ifndef DOXYGEN

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
inline void
SelectEvaluator<dim, fe_degree, n_q_points_1d, Number>::evaluate(
  const unsigned int                     n_components,
  const EvaluationFlags::EvaluationFlags evaluation_flag,
  const Number *                         values_dofs,
  FEEvaluationData<dim, Number, false> & eval)
{
  Assert(fe_degree >= 0 && n_q_points_1d > 0, ExcInternalError());

  internal::FEEvaluationImplSelector<dim, Number, false>::template run<
    fe_degree,
    n_q_points_1d>(n_components, evaluation_flag, values_dofs, eval);
}



template <int dim, int fe_degree, int n_q_points_1d, typename Number>
inline void
SelectEvaluator<dim, fe_degree, n_q_points_1d, Number>::integrate(
  const unsigned int                     n_components,
  const EvaluationFlags::EvaluationFlags integration_flag,
  Number *                               values_dofs,
  FEEvaluationData<dim, Number, false> & eval,
  const bool                             sum_into_values_array)
{
  Assert(fe_degree >= 0 && n_q_points_1d > 0, ExcInternalError());

  internal::FEEvaluationImplSelector<dim, Number, true>::
    template run<fe_degree, n_q_points_1d>(
      n_components, integration_flag, values_dofs, eval, sum_into_values_array);
}
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
