// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
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

#include <deal.II/matrix_free/evaluation_kernels_face.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/evaluation_template_factory_internal.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::evaluate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number                          *values_dofs,
    FEEvaluationData<dim, Number, true>   &fe_eval)
  {
    instantiation_helper_run<1,
                             FEFaceEvaluationImplEvaluateSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      evaluation_flag,
      values_dofs,
      fe_eval);
  }



  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::project_to_face(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number                          *values_dofs,
    FEEvaluationData<dim, Number, true>   &fe_eval)
  {
    instantiation_helper_degree_run<
      1,
      FEFaceEvaluationImplProjectToFaceSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      n_components,
      evaluation_flag,
      values_dofs,
      fe_eval);
  }



  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::evaluate_in_face(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    FEEvaluationData<dim, Number, true>   &fe_eval)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplEvaluateInFaceSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      evaluation_flag,
      fe_eval);
  }



  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::integrate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    Number                                *values_dofs,
    FEEvaluationData<dim, Number, true>   &fe_eval,
    const bool                             sum_into_values)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplIntegrateSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      integration_flag,
      values_dofs,
      fe_eval,
      sum_into_values);
  }



  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::collect_from_face(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    Number                                *values_dofs,
    FEEvaluationData<dim, Number, true>   &fe_eval,
    const bool                             sum_into_values)
  {
    instantiation_helper_degree_run<
      1,
      FEFaceEvaluationImplCollectFromFaceSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      n_components,
      integration_flag,
      values_dofs,
      fe_eval,
      sum_into_values);
  }



  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::integrate_in_face(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    FEEvaluationData<dim, Number, true>   &fe_eval)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplIntegrateInFaceSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      integration_flag,
      fe_eval);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEFaceEvaluationGatherFactory<dim, Number, VectorizedArrayType>::evaluate(
    const unsigned int                                n_components,
    const EvaluationFlags::EvaluationFlags            evaluation_flag,
    const Number                                     *src_ptr,
    const std::vector<ArrayView<const Number>>       *sm_ptr,
    FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplGatherEvaluateSelector<dim,
                                                 Number,
                                                 VectorizedArrayType>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      evaluation_flag,
      src_ptr,
      sm_ptr,
      fe_eval);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEFaceEvaluationGatherFactory<dim, Number, VectorizedArrayType>::integrate(
    const unsigned int                                n_components,
    const EvaluationFlags::EvaluationFlags            integration_flag,
    Number                                           *dst_ptr,
    const std::vector<ArrayView<const Number>>       *sm_ptr,
    FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplIntegrateScatterSelector<dim,
                                                   Number,
                                                   VectorizedArrayType>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      integration_flag,
      dst_ptr,
      sm_ptr,
      fe_eval);
  }



  // It is important that this file sits in the same compilation unit as the
  // evaluate() and integrate() calls of this class, to ensure that all
  // options choose the same code path when compiling FEFaceEvaluationFactory
  // outside of deal.II.
  template <int dim, typename Number>
  bool
  FEFaceEvaluationFactory<dim, Number>::fast_evaluation_supported(
    const unsigned int given_degree,
    const unsigned int n_q_points_1d)
  {
    return instantiation_helper_run<1, FastEvaluationSupported>(given_degree,
                                                                n_q_points_1d);
  }
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
