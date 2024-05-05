// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_evaluation_template_factory_h
#define dealii_matrix_free_evaluation_template_factory_h


#include <deal.II/base/config.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN


template <int, typename, bool>
class FEEvaluationData;


namespace internal
{
  template <int dim, typename Number>
  struct FEEvaluationFactory
  {
    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number                          *values_dofs,
             FEEvaluationData<dim, Number, false>  &fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number                                *values_dofs,
              FEEvaluationData<dim, Number, false>  &fe_eval,
              const bool                             sum_into_values_array);

    static bool
    fast_evaluation_supported(const unsigned int given_degree,
                              const unsigned int n_q_points_1d);
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationFactory
  {
    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number                          *values_dofs,
             FEEvaluationData<dim, Number, true>   &fe_eval);

    static void
    project_to_face(const unsigned int                     n_components,
                    const EvaluationFlags::EvaluationFlags evaluation_flag,
                    const Number                          *values_dofs,
                    FEEvaluationData<dim, Number, true>   &fe_eval);

    static void
    evaluate_in_face(const unsigned int                     n_components,
                     const EvaluationFlags::EvaluationFlags evaluation_flag,
                     FEEvaluationData<dim, Number, true>   &fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number                                *values_dofs,
              FEEvaluationData<dim, Number, true>   &fe_eval,
              const bool                             sum_into_values);

    static void
    collect_from_face(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags integration_flag,
                      Number                                *values_dofs,
                      FEEvaluationData<dim, Number, true>   &fe_eval,
                      const bool                             sum_into_values);

    static void
    integrate_in_face(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags integration_flag,
                      FEEvaluationData<dim, Number, true>   &fe_eval);

    static bool
    fast_evaluation_supported(const unsigned int given_degree,
                              const unsigned int n_q_points_1d);
  };



  template <int dim, typename Number, typename VectorizedArrayType>
  struct FEFaceEvaluationGatherFactory
  {
    static void
    evaluate(const unsigned int                                n_components,
             const EvaluationFlags::EvaluationFlags            evaluation_flag,
             const Number                                     *src_ptr,
             const std::vector<ArrayView<const Number>>       *sm_ptr,
             FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval);

    static void
    integrate(const unsigned int                          n_components,
              const EvaluationFlags::EvaluationFlags      integration_flag,
              Number                                     *dst_ptr,
              const std::vector<ArrayView<const Number>> *sm_ptr,
              FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval);
  };



  template <int dim, typename Number>
  struct CellwiseInverseMassFactory
  {
    static void
    apply(const unsigned int                          n_components,
          const FEEvaluationData<dim, Number, false> &fe_eval,
          const Number                               *in_array,
          Number                                     *out_array);

    static void
    apply(const unsigned int                          n_components,
          const FEEvaluationData<dim, Number, false> &fe_eval,
          const ArrayView<const Number>              &inverse_coefficients,
          const bool                                  dyadic_coefficients,
          const Number                               *in_array,
          Number                                     *out_array);

    static void
    transform_from_q_points_to_basis(
      const unsigned int                          n_components,
      const FEEvaluationData<dim, Number, false> &fe_eval,
      const Number                               *in_array,
      Number                                     *out_array);
  };

  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  struct FEEvaluationHangingNodesFactory
  {
    static void
    apply(const unsigned int                             n_components,
          const unsigned int                             fe_degree,
          const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
          const bool                                     transpose,
          const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                           VectorizedArrayType::size()> &c_mask,
          VectorizedArrayType                           *values);
  };

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
