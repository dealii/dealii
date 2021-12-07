// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_template_factory_hanging_nodes_templates_h
#define dealii_matrix_free_evaluation_template_factory_hanging_nodes_templates_h


#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_kernels_hanging_nodes.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/evaluation_template_factory_internal.h>
#include <deal.II/matrix_free/fe_evaluation_base_data.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEEvaluationHangingNodesFactory<dim, Number, VectorizedArrayType>::apply(
    const unsigned int n_components,
    const unsigned int fe_degree,
    const FEEvaluationBaseData<dim, Number, false, VectorizedArrayType>
      &                                            fe_eval,
    const bool                                     transpose,
    const std::array<MatrixFreeFunctions::ConstraintKinds,
                     VectorizedArrayType::size()> &c_mask,
    VectorizedArrayType *                          values)
  {
    instantiation_helper_degree_run<
      1,
      FEEvaluationImplHangingNodes<dim, VectorizedArrayType, false>>(
      fe_degree, n_components, fe_eval, transpose, c_mask, values);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEEvaluationHangingNodesFactory<dim, Number, VectorizedArrayType>::apply(
    const unsigned int n_components,
    const unsigned int fe_degree,
    const FEEvaluationBaseData<dim, Number, true, VectorizedArrayType> &fe_eval,
    const bool                                     transpose,
    const std::array<MatrixFreeFunctions::ConstraintKinds,
                     VectorizedArrayType::size()> &c_mask,
    VectorizedArrayType *                          values)
  {
    instantiation_helper_degree_run<
      1,
      FEEvaluationImplHangingNodes<dim, VectorizedArrayType, true>>(
      fe_degree, n_components, fe_eval, transpose, c_mask, values);
  }
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
