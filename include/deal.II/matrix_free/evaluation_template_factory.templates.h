// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
  template <int dim, typename Number>
  void
  FEEvaluationFactory<dim, Number>::evaluate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number                          *values_dofs,
    FEEvaluationData<dim, Number, false>  &fe_eval)
  {
    instantiation_helper_run<1, FEEvaluationImplSelector<dim, Number, false>>(
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
    Number                                *values_dofs,
    FEEvaluationData<dim, Number, false>  &fe_eval,
    const bool                             sum_into_values_array)
  {
    instantiation_helper_run<1, FEEvaluationImplSelector<dim, Number, true>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      integration_flag,
      values_dofs,
      fe_eval,
      sum_into_values_array);
  }



  // It is important that this file sits in the same compilation unit as the
  // evaluate() and integrate() calls of this class, to ensure that all
  // options choose the same code path when compiling FEFaceEvaluationFactory
  // outside of deal.II.
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
    const Number                               *in_array,
    Number                                     *out_array)
  {
    const unsigned int fe_degree = fe_eval.get_shape_info().data[0].fe_degree;
    instantiation_helper_run<1,
                             CellwiseInverseMassMatrixImplBasic<dim, Number>>(
      fe_degree, fe_degree + 1, n_components, fe_eval, in_array, out_array);
  }



  template <int dim, typename Number>
  void
  CellwiseInverseMassFactory<dim, Number>::apply(
    const unsigned int                          n_components,
    const FEEvaluationData<dim, Number, false> &fe_eval,
    const ArrayView<const Number>              &inverse_coefficients,
    const bool                                  dyadic_coefficients,
    const Number                               *in_array,
    Number                                     *out_array)
  {
    const unsigned int fe_degree = fe_eval.get_shape_info().data[0].fe_degree;
    instantiation_helper_run<
      1,
      CellwiseInverseMassMatrixImplFlexible<dim, Number>>(fe_degree,
                                                          fe_degree + 1,
                                                          n_components,
                                                          fe_eval,
                                                          inverse_coefficients,
                                                          dyadic_coefficients,
                                                          in_array,
                                                          out_array);
  }



  template <int dim, typename Number>
  void
  CellwiseInverseMassFactory<dim, Number>::transform_from_q_points_to_basis(
    const unsigned int                          n_components,
    const FEEvaluationData<dim, Number, false> &fe_eval,
    const Number                               *in_array,
    Number                                     *out_array)
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
