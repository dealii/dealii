// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_evaluation_template_factory_hanging_nodes_templates_h
#define dealii_matrix_free_evaluation_template_factory_hanging_nodes_templates_h


#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_kernels_hanging_nodes.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/evaluation_template_factory_internal.h>
#include <deal.II/matrix_free/shape_info.h>

#include <array>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEEvaluationHangingNodesFactory<dim, Number, VectorizedArrayType>::apply(
    const unsigned int                             n_components,
    const unsigned int                             fe_degree,
    const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
    const bool                                     transpose,
    const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                     VectorizedArrayType::size()> &c_mask,
    VectorizedArrayType                           *values)
  {
    instantiation_helper_degree_run<
      1,
      FEEvaluationImplHangingNodes<dim, VectorizedArrayType>>(
      fe_degree, n_components, shape_info, transpose, c_mask, values);
  }
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
