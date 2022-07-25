// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_template_factory_templates_h
#define dealii_matrix_free_evaluation_template_factory_templates_h


#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/evaluation_selector.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/evaluation_template_factory_internal.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  struct FastEvaluationSupported
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run()
    {
      return fe_degree != -1;
    }
  };



  template <int dim, typename Number>
  void
  FEEvaluationFactory<dim, Number>::evaluate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number *                         values_dofs,
    FEEvaluationData<dim, Number, false> & fe_eval)
  {
    instantiation_helper_run<1, FEEvaluationImplEvaluateSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      evaluation_flag,
      values_dofs,
      fe_eval);
  }



  template <int dim, typename Number>
  void
  FEEvaluationFactory<dim, Number>::integrate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    Number *                               values_dofs,
    FEEvaluationData<dim, Number, false> & fe_eval,
    const bool                             sum_into_values_array)
  {
    instantiation_helper_run<1, FEEvaluationImplIntegrateSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      integration_flag,
      values_dofs,
      fe_eval,
      sum_into_values_array);
  }



  template <int dim, typename Number>
  bool
  FEEvaluationFactory<dim, Number>::fast_evaluation_supported(
    const unsigned int given_degree,
    const unsigned int n_q_points_1d)
  {
    return instantiation_helper_run<1, FastEvaluationSupported>(given_degree,
                                                                n_q_points_1d);
  }



  template <int dim, typename Number>
  void
  CellwiseInverseMassFactory<dim, Number>::apply(
    const unsigned int                          n_components,
    const FEEvaluationData<dim, Number, false> &fe_eval,
    const Number *                              in_array,
    Number *                                    out_array)
  {
    const unsigned int fe_degree = fe_eval.get_shape_info().data[0].fe_degree;
    instantiation_helper_run<1,
                             CellwiseInverseMassMatrixImplBasic<dim, Number>>(
      fe_degree, fe_degree + 1, n_components, fe_eval, in_array, out_array);
  }



  template <int dim, typename Number>
  void
  CellwiseInverseMassFactory<dim, Number>::apply(
    const unsigned int           n_components,
    const unsigned int           fe_degree,
    const AlignedVector<Number> &inverse_shape,
    const AlignedVector<Number> &inverse_coefficients,
    const Number *               in_array,
    Number *                     out_array)
  {
    instantiation_helper_run<
      1,
      CellwiseInverseMassMatrixImplFlexible<dim, Number>>(fe_degree,
                                                          fe_degree + 1,
                                                          n_components,
                                                          inverse_shape,
                                                          inverse_coefficients,
                                                          in_array,
                                                          out_array);
  }



  template <int dim, typename Number>
  void
  CellwiseInverseMassFactory<dim, Number>::transform_from_q_points_to_basis(
    const unsigned int                          n_components,
    const FEEvaluationData<dim, Number, false> &fe_eval,
    const Number *                              in_array,
    Number *                                    out_array)
  {
    const unsigned int fe_degree = fe_eval.get_shape_info().data[0].fe_degree;
    const unsigned int n_q_points_1d =
      fe_eval.get_shape_info().data[0].n_q_points_1d;
    instantiation_helper_run<
      1,
      CellwiseInverseMassMatrixImplTransformFromQPoints<dim, Number>>(
      fe_degree, n_q_points_1d, n_components, fe_eval, in_array, out_array);
  }

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
