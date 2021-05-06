// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


template <int, typename, bool, typename>
class FEEvaluationBaseData;


namespace internal
{
  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  struct FEEvaluationFactory
  {
    static void
    evaluate(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &shape_info,
      VectorizedArrayType *values_dofs_actual,
      VectorizedArrayType *values_quad,
      VectorizedArrayType *gradients_quad,
      VectorizedArrayType *hessians_quad,
      VectorizedArrayType *scratch_data);

    static void
    integrate(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags integration_flag,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &shape_info,
      VectorizedArrayType *values_dofs_actual,
      VectorizedArrayType *values_quad,
      VectorizedArrayType *gradients_quad,
      VectorizedArrayType *scratch_data,
      const bool           sum_into_values_array);

    static bool
    fast_evaluation_supported(const unsigned int given_degree,
                              const unsigned int n_q_points_1d);
  };



  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  struct FEFaceEvaluationFactory
  {
    static void
    evaluate(const unsigned int n_components,
             const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
             const VectorizedArrayType *   values_array,
             VectorizedArrayType *         values_quad,
             VectorizedArrayType *         gradients_quad,
             VectorizedArrayType *         scratch_data,
             const bool                    evaluate_values,
             const bool                    evaluate_gradients,
             const unsigned int            face_no,
             const unsigned int            subface_index,
             const unsigned int            face_orientation,
             const Table<2, unsigned int> &orientation_map);

    static void
    integrate(const unsigned int n_components,
              const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
              VectorizedArrayType *         values_array,
              VectorizedArrayType *         values_quad,
              VectorizedArrayType *         gradients_quad,
              VectorizedArrayType *         scratch_data,
              const bool                    integrate_values,
              const bool                    integrate_gradients,
              const unsigned int            face_no,
              const unsigned int            subface_index,
              const unsigned int            face_orientation,
              const Table<2, unsigned int> &orientation_map);

    static bool
    gather_evaluate(
      const unsigned int                          n_components,
      const std::size_t                           n_face_orientations,
      const Number *                              src_ptr,
      const std::vector<ArrayView<const Number>> *sm_ptr,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
      const MatrixFreeFunctions::DoFInfo &                       dof_info,
      VectorizedArrayType *                                      values_quad,
      VectorizedArrayType *                                      gradients_quad,
      VectorizedArrayType *                                      scratch_data,
      const bool         evaluate_values,
      const bool         evaluate_gradients,
      const unsigned int active_fe_index,
      const unsigned int first_selected_component,
      const std::array<unsigned int, VectorizedArrayType::size()> cells,
      const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
      const unsigned int                                          subface_index,
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
      const std::array<unsigned int, VectorizedArrayType::size()>
                                    face_orientations,
      const Table<2, unsigned int> &orientation_map);

    static bool
    integrate_scatter(
      const unsigned int                          n_components,
      const std::size_t                           n_face_orientations,
      Number *                                    dst_ptr,
      const std::vector<ArrayView<const Number>> *sm_ptr,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
      const MatrixFreeFunctions::DoFInfo &                       dof_info,
      VectorizedArrayType *                                      values_array,
      VectorizedArrayType *                                      values_quad,
      VectorizedArrayType *                                      gradients_quad,
      VectorizedArrayType *                                      scratch_data,
      const bool         integrate_values,
      const bool         integrate_gradients,
      const unsigned int active_fe_index,
      const unsigned int first_selected_component,
      const std::array<unsigned int, VectorizedArrayType::size()> cells,
      const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
      const unsigned int                                          subface_index,
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
      const std::array<unsigned int, VectorizedArrayType::size()>
                                    face_orientations,
      const Table<2, unsigned int> &orientation_map);

    static bool
    fast_evaluation_supported(const unsigned int given_degree,
                              const unsigned int n_q_points_1d);
  };



  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  struct CellwiseInverseMassFactory
  {
    static void
    apply(const unsigned int n_components,
          const unsigned int fe_degree,
          const FEEvaluationBaseData<dim, Number, false, VectorizedArrayType>
            &                        fe_eval,
          const VectorizedArrayType *in_array,
          VectorizedArrayType *      out_array);

    static void
    apply(const unsigned int                        n_components,
          const unsigned int                        fe_degree,
          const AlignedVector<VectorizedArrayType> &inverse_shape,
          const AlignedVector<VectorizedArrayType> &inverse_coefficients,
          const VectorizedArrayType *               in_array,
          VectorizedArrayType *                     out_array);

    static void
    transform_from_q_points_to_basis(
      const unsigned int                        n_components,
      const unsigned int                        fe_degree,
      const AlignedVector<VectorizedArrayType> &inverse_shape,
      const VectorizedArrayType *               in_array,
      VectorizedArrayType *                     out_array);
  };

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
