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
             const Number *                         values_dofs,
             FEEvaluationData<dim, Number, false> & fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number *                               values_dofs,
              FEEvaluationData<dim, Number, false> & fe_eval,
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
             const Number *                         values_dofs,
             FEEvaluationData<dim, Number, true> &  fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number *                               values_dofs,
              FEEvaluationData<dim, Number, true> &  fe_eval);
  };



  template <int dim, typename Number, typename VectorizedArrayType>
  struct FEFaceEvaluationGatherFactory
  {
    static void
    evaluate(const unsigned int                                n_components,
             const EvaluationFlags::EvaluationFlags            evaluation_flag,
             const Number *                                    src_ptr,
             const std::vector<ArrayView<const Number>> *      sm_ptr,
             FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval);

    static void
    integrate(const unsigned int                          n_components,
              const EvaluationFlags::EvaluationFlags      integration_flag,
              Number *                                    dst_ptr,
              const std::vector<ArrayView<const Number>> *sm_ptr,
              FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval);
  };



  template <int dim, typename Number>
  struct CellwiseInverseMassFactory
  {
    static void
    apply(const unsigned int                          n_components,
          const FEEvaluationData<dim, Number, false> &fe_eval,
          const Number *                              in_array,
          Number *                                    out_array);

    static void
    apply(const unsigned int           n_components,
          const unsigned int           fe_degree,
          const AlignedVector<Number> &inverse_shape,
          const AlignedVector<Number> &inverse_coefficients,
          const Number *               in_array,
          Number *                     out_array);

    static void
    transform_from_q_points_to_basis(
      const unsigned int                          n_components,
      const FEEvaluationData<dim, Number, false> &fe_eval,
      const Number *                              in_array,
      Number *                                    out_array);
  };

  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  struct FEEvaluationHangingNodesFactory
  {
    static void
    apply(const unsigned int n_components,
          const unsigned int fe_degree,
          const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &shape_info,
          const bool                                                 transpose,
          const std::array<MatrixFreeFunctions::compressed_constraint_kind,
                           VectorizedArrayType::size()> &            c_mask,
          VectorizedArrayType *                                      values);
  };

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
