// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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
  template <int dim, typename Number>
  void
  FEFaceEvaluationFactory<dim, Number>::evaluate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number *                         values_dofs,
    FEEvaluationData<dim, Number, true> &  fe_eval)
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
  FEFaceEvaluationFactory<dim, Number>::integrate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    Number *                               values_dofs,
    FEEvaluationData<dim, Number, true> &  fe_eval)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplIntegrateSelector<dim, Number>>(
      fe_eval.get_shape_info().data[0].fe_degree,
      fe_eval.get_shape_info().data[0].n_q_points_1d,
      n_components,
      integration_flag,
      values_dofs,
      fe_eval);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEFaceEvaluationGatherFactory<dim, Number, VectorizedArrayType>::evaluate(
    const unsigned int                                n_components,
    const EvaluationFlags::EvaluationFlags            evaluation_flag,
    const Number *                                    src_ptr,
    const std::vector<ArrayView<const Number>> *      sm_ptr,
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
    Number *                                          dst_ptr,
    const std::vector<ArrayView<const Number>> *      sm_ptr,
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

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
