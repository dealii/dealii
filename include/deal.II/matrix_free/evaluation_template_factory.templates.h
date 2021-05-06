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


#ifndef dealii_matrix_free_evaluation_template_factory_templates_h
#define dealii_matrix_free_evaluation_template_factory_templates_h


#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/evaluation_selector.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#ifndef FE_EVAL_FACTORY_DEGREE_MAX
#  define FE_EVAL_FACTORY_DEGREE_MAX 6
#endif

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <int degree, typename EvaluatorType, typename... Args>
  bool
  instantiation_helper_run(const unsigned int given_degree,
                           const unsigned int n_q_points_1d,
                           Args &... args)
  {
    if (given_degree == degree)
      {
        if (n_q_points_1d == degree + 1)
          return EvaluatorType::template run<degree, degree + 1>(args...);
        else if (n_q_points_1d == degree + 2)
          return EvaluatorType::template run<degree, degree + 2>(args...);
        else if (n_q_points_1d == degree)
          return EvaluatorType::template run<degree, degree>(args...);
        else if (n_q_points_1d == (3 * degree) / 2 + 1)
          return EvaluatorType::template run<degree, (3 * degree) / 2 + 1>(
            args...);
        else
          // slow path
          return EvaluatorType::template run<-1, 0>(args...);
      }
    else if (degree < FE_EVAL_FACTORY_DEGREE_MAX)
      return instantiation_helper_run<
        (degree < FE_EVAL_FACTORY_DEGREE_MAX ? degree + 1 : degree),
        EvaluatorType>(given_degree, n_q_points_1d, args...);
    else
      // slow path
      return EvaluatorType::template run<-1, 0>(args...);
  }

  struct FastEvaluationSupported
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run()
    {
      return fe_degree != -1;
    }
  };



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEEvaluationFactory<dim, Number, VectorizedArrayType>::evaluate(
    const unsigned int                                         n_components,
    const EvaluationFlags::EvaluationFlags                     evaluation_flag,
    const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &shape_info,
    VectorizedArrayType *values_dofs_actual,
    VectorizedArrayType *values_quad,
    VectorizedArrayType *gradients_quad,
    VectorizedArrayType *hessians_quad,
    VectorizedArrayType *scratch_data)
  {
    instantiation_helper_run<
      1,
      FEEvaluationImplEvaluateSelector<dim, VectorizedArrayType>>(
      shape_info.data[0].fe_degree,
      shape_info.data[0].n_q_points_1d,
      n_components,
      evaluation_flag,
      shape_info,
      values_dofs_actual,
      values_quad,
      gradients_quad,
      hessians_quad,
      scratch_data);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEEvaluationFactory<dim, Number, VectorizedArrayType>::integrate(
    const unsigned int                                         n_components,
    const EvaluationFlags::EvaluationFlags                     integration_flag,
    const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &shape_info,
    VectorizedArrayType *values_dofs_actual,
    VectorizedArrayType *values_quad,
    VectorizedArrayType *gradients_quad,
    VectorizedArrayType *scratch_data,
    const bool           sum_into_values_array)
  {
    instantiation_helper_run<
      1,
      FEEvaluationImplIntegrateSelector<dim, VectorizedArrayType>>(
      shape_info.data[0].fe_degree,
      shape_info.data[0].n_q_points_1d,
      n_components,
      integration_flag,
      shape_info,
      values_dofs_actual,
      values_quad,
      gradients_quad,
      scratch_data,
      sum_into_values_array);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  bool
  FEEvaluationFactory<dim, Number, VectorizedArrayType>::
    fast_evaluation_supported(const unsigned int given_degree,
                              const unsigned int n_q_points_1d)
  {
    return instantiation_helper_run<1, FastEvaluationSupported>(given_degree,
                                                                n_q_points_1d);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::evaluate(
    const unsigned int                                         n_components,
    const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
    const VectorizedArrayType *                                values_array,
    VectorizedArrayType *                                      values_quad,
    VectorizedArrayType *                                      gradients_quad,
    VectorizedArrayType *                                      scratch_data,
    const bool                                                 evaluate_values,
    const bool                    evaluate_gradients,
    const unsigned int            face_no,
    const unsigned int            subface_index,
    const unsigned int            face_orientation,
    const Table<2, unsigned int> &orientation_map)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplEvaluateSelector<dim, VectorizedArrayType>>(
      data.data[0].fe_degree,
      data.data[0].n_q_points_1d,
      n_components,
      data,
      values_array,
      values_quad,
      gradients_quad,
      scratch_data,
      evaluate_values,
      evaluate_gradients,
      face_no,
      subface_index,
      face_orientation,
      orientation_map);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::integrate(
    const unsigned int                                         n_components,
    const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
    VectorizedArrayType *                                      values_array,
    VectorizedArrayType *                                      values_quad,
    VectorizedArrayType *                                      gradients_quad,
    VectorizedArrayType *                                      scratch_data,
    const bool                                                 integrate_values,
    const bool                    integrate_gradients,
    const unsigned int            face_no,
    const unsigned int            subface_index,
    const unsigned int            face_orientation,
    const Table<2, unsigned int> &orientation_map)
  {
    instantiation_helper_run<
      1,
      FEFaceEvaluationImplIntegrateSelector<dim, VectorizedArrayType>>(
      data.data[0].fe_degree,
      data.data[0].n_q_points_1d,
      n_components,
      data,
      values_array,
      values_quad,
      gradients_quad,
      scratch_data,
      integrate_values,
      integrate_gradients,
      face_no,
      subface_index,
      face_orientation,
      orientation_map);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  bool
  FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::gather_evaluate(
    const unsigned int                          n_components,
    const std::size_t                           n_face_orientations,
    const Number *                              src_ptr,
    const std::vector<ArrayView<const Number>> *sm_ptr,
    const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
    const MatrixFreeFunctions::DoFInfo &                       dof_info,
    VectorizedArrayType *                                      values_quad,
    VectorizedArrayType *                                      gradients_quad,
    VectorizedArrayType *                                      scratch_data,
    const bool                                                 evaluate_values,
    const bool         evaluate_gradients,
    const unsigned int active_fe_index,
    const unsigned int first_selected_component,
    const std::array<unsigned int, VectorizedArrayType::size()> cells,
    const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
    const unsigned int                                          subface_index,
    const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
    const std::array<unsigned int, VectorizedArrayType::size()>
                                  face_orientations,
    const Table<2, unsigned int> &orientation_map)
  {
    return instantiation_helper_run<
      1,
      FEFaceEvaluationImplGatherEvaluateSelector<dim,
                                                 Number,
                                                 VectorizedArrayType>>(
      data.data[0].fe_degree,
      data.data[0].n_q_points_1d,
      n_components,
      n_face_orientations,
      src_ptr,
      sm_ptr,
      data,
      dof_info,
      values_quad,
      gradients_quad,
      scratch_data,
      evaluate_values,
      evaluate_gradients,
      active_fe_index,
      first_selected_component,
      cells,
      face_nos,
      subface_index,
      dof_access_index,
      face_orientations,
      orientation_map);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  bool
  FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::integrate_scatter(
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
    const bool                                                 integrate_values,
    const bool         integrate_gradients,
    const unsigned int active_fe_index,
    const unsigned int first_selected_component,
    const std::array<unsigned int, VectorizedArrayType::size()> cells,
    const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
    const unsigned int                                          subface_index,
    const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
    const std::array<unsigned int, VectorizedArrayType::size()>
                                  face_orientations,
    const Table<2, unsigned int> &orientation_map)
  {
    return instantiation_helper_run<
      1,
      FEFaceEvaluationImplIntegrateScatterSelector<dim,
                                                   Number,
                                                   VectorizedArrayType>>(
      data.data[0].fe_degree,
      data.data[0].n_q_points_1d,
      n_components,
      n_face_orientations,
      dst_ptr,
      sm_ptr,
      data,
      dof_info,
      values_array,
      values_quad,
      gradients_quad,
      scratch_data,
      integrate_values,
      integrate_gradients,
      active_fe_index,
      first_selected_component,
      cells,
      face_nos,
      subface_index,
      dof_access_index,
      face_orientations,
      orientation_map);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  bool
  FEFaceEvaluationFactory<dim, Number, VectorizedArrayType>::
    fast_evaluation_supported(const unsigned int given_degree,
                              const unsigned int n_q_points_1d)
  {
    return instantiation_helper_run<1, FastEvaluationSupported>(given_degree,
                                                                n_q_points_1d);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  CellwiseInverseMassFactory<dim, Number, VectorizedArrayType>::apply(
    const unsigned int n_components,
    const unsigned int fe_degree,
    const FEEvaluationBaseData<dim, Number, false, VectorizedArrayType>
      &                        fe_eval,
    const VectorizedArrayType *in_array,
    VectorizedArrayType *      out_array)
  {
    instantiation_helper_run<
      1,
      CellwiseInverseMassMatrixImplBasic<dim, VectorizedArrayType>>(
      fe_degree, fe_degree + 1, n_components, fe_eval, in_array, out_array);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  CellwiseInverseMassFactory<dim, Number, VectorizedArrayType>::apply(
    const unsigned int                        n_components,
    const unsigned int                        fe_degree,
    const AlignedVector<VectorizedArrayType> &inverse_shape,
    const AlignedVector<VectorizedArrayType> &inverse_coefficients,
    const VectorizedArrayType *               in_array,
    VectorizedArrayType *                     out_array)
  {
    instantiation_helper_run<
      1,
      CellwiseInverseMassMatrixImplFlexible<dim, VectorizedArrayType>>(
      fe_degree,
      fe_degree + 1,
      n_components,
      inverse_shape,
      inverse_coefficients,
      in_array,
      out_array);
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  void
  CellwiseInverseMassFactory<dim, Number, VectorizedArrayType>::
    transform_from_q_points_to_basis(
      const unsigned int                        n_components,
      const unsigned int                        fe_degree,
      const AlignedVector<VectorizedArrayType> &inverse_shape,
      const VectorizedArrayType *               in_array,
      VectorizedArrayType *                     out_array)
  {
    instantiation_helper_run<
      1,
      CellwiseInverseMassMatrixImplTransformFromQPoints<dim,
                                                        VectorizedArrayType>>(
      fe_degree,
      fe_degree + 1,
      n_components,
      inverse_shape,
      in_array,
      out_array);
  }

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
