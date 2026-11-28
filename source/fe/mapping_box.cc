// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_box.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>



DEAL_II_NAMESPACE_OPEN

DeclExceptionMsg(
  ExcCellNotAssociatedWithBox,
  "You are using MappingBox, but the incoming element is not associated with a"
  "Bounding Box.");



/**
 * Return whether the incoming element has a BoundingBox associated to it.
 * Simplicial and quad-hex meshes are supported.
 */
template <typename CellType>
bool
has_box(const CellType &cell,
        const std::map<types::global_cell_index, types::global_cell_index>
          &translator)
{
  Assert((cell->reference_cell().is_hyper_cube() ||
          cell->reference_cell().is_simplex()),
         ExcNotImplemented());
  Assert((translator.find(cell->active_cell_index()) != translator.cend()),
         ExcCellNotAssociatedWithBox());

  return true;
}



template <int dim, int spacedim>
MappingBox<dim, spacedim>::MappingBox(
  const std::vector<BoundingBox<dim>> &input_boxes,
  const std::map<types::global_cell_index, types::global_cell_index>
    &global_to_polytope)
{
  Assert(input_boxes.size() > 0,
         ExcMessage("Invalid number of bounding boxes."));

  // copy boxes and map
  boxes.resize(input_boxes.size());
  for (unsigned int i = 0; i < input_boxes.size(); ++i)
    boxes[i] = input_boxes[i];
  polytope_translator = global_to_polytope;
}



template <int dim, int spacedim>
MappingBox<dim, spacedim>::InternalData::InternalData(const Quadrature<dim> &q)
  : cell_extents(numbers::signaling_nan<Tensor<1, dim>>())
  , traslation(numbers::signaling_nan<Tensor<1, dim>>())
  , inverse_cell_extents(numbers::signaling_nan<Tensor<1, dim>>())
  , volume_element(numbers::signaling_nan<double>())
  , quadrature_points(q.get_points())
{}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::InternalData::reinit(const UpdateFlags update_flags,
                                                const Quadrature<dim> &)
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values(). use the transitive hull of the required
  // flags
  this->update_each = update_flags;
}



template <int dim, int spacedim>
std::size_t
MappingBox<dim, spacedim>::InternalData::memory_consumption() const
{
  return (Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption(cell_extents) +
          MemoryConsumption::memory_consumption(traslation) +
          MemoryConsumption::memory_consumption(inverse_cell_extents) +
          MemoryConsumption::memory_consumption(volume_element));
}



template <int dim, int spacedim>
bool
MappingBox<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



template <int dim, int spacedim>
bool
#if DEAL_II_VERSION_GTE(9, 8, 0)
MappingBox<dim, spacedim>::is_compatible_with(
  const ReferenceCell<dim> &reference_cell) const
#else
MappingBox<dim, spacedim>::is_compatible_with(
  const ReferenceCell &reference_cell) const
#endif
{
  Assert(dim == reference_cell.get_dimension(),
         ExcMessage("The dimension of your mapping (" +
                    Utilities::to_string(dim) +
                    ") and the reference cell cell_type (" +
                    Utilities::to_string(reference_cell.get_dimension()) +
                    " ) do not agree."));

  return reference_cell.is_hyper_cube() || reference_cell.is_simplex();
}



template <int dim, int spacedim>
UpdateFlags
MappingBox<dim, spacedim>::requires_update_flags(const UpdateFlags in) const
{
  // this mapping is pretty simple in that it can basically compute
  // every piece of information wanted by FEValues without requiring
  // computing any other quantities. boundary forms are one exception
  // since they can be computed from the normal vectors without much
  // further ado
  UpdateFlags out = in;
  if (out & update_boundary_forms)
    out |= update_normal_vectors;

  return out;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingBox<dim, spacedim>::get_data(const UpdateFlags      update_flags,
                                    const Quadrature<dim> &q) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>();
  data_ptr->reinit(requires_update_flags(update_flags), q);

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingBox<dim, spacedim>::get_subface_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  (void)update_flags;
  (void)quadrature;
  DEAL_II_NOT_IMPLEMENTED();
  return {};
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::update_cell_extents(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const InternalData                                         &data) const
{
  // Compute start point and sizes along axes. The vertices to be looked at
  // are 1, 2, 4 compared to the base vertex 0.
  if (cell_similarity != CellSimilarity::translation)
    {
      const BoundingBox<dim> &current_box =
        boxes[polytope_translator.at(cell->active_cell_index())];
      const std::pair<Point<dim>, Point<dim>> &bdary_points =
        current_box.get_boundary_points();

      for (unsigned int d = 0; d < dim; ++d)
        {
          const double cell_extent_d = current_box.side_length(d);
          data.cell_extents[d]       = cell_extent_d;

          data.traslation[d] =
            .5 * (bdary_points.first[d] +
                  bdary_points.second[d]); // midpoint of each interval

          Assert(cell_extent_d != 0.,
                 ExcMessage("Cell does not appear to be Cartesian!"));
          data.inverse_cell_extents[d] = 1. / cell_extent_d;
        }
    }
}



namespace
{
  template <int dim>
  void
  transform_quadrature_points(
    const BoundingBox<dim>            &box,
    const ArrayView<const Point<dim>> &unit_quadrature_points,
    std::vector<Point<dim>>           &quadrature_points)
  {
    for (unsigned int i = 0; i < quadrature_points.size(); ++i)
      quadrature_points[i] = box.unit_to_real(unit_quadrature_points[i]);
  }
} // namespace



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::maybe_update_cell_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const InternalData                                         &data,
  const ArrayView<const Point<dim>> &unit_quadrature_points,
  std::vector<Point<dim>>           &quadrature_points) const
{
  if (data.update_each & update_quadrature_points)
    transform_quadrature_points(
      boxes[polytope_translator.at(cell->active_cell_index())],
      unit_quadrature_points,
      quadrature_points);
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::maybe_update_normal_vectors(
  const unsigned int           face_no,
  const InternalData          &data,
  std::vector<Tensor<1, dim>> &normal_vectors) const
{
  // compute normal vectors. All normals on a face have the same value.
  if (data.update_each & update_normal_vectors)
    {
      Assert(face_no < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
      std::fill(normal_vectors.begin(),
                normal_vectors.end(),
                GeometryInfo<dim>::unit_normal_vector[face_no]);
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::maybe_update_jacobian_derivatives(
  const InternalData              &data,
  const CellSimilarity::Similarity cell_similarity,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  if (cell_similarity != CellSimilarity::translation)
    {
      if (data.update_each & update_jacobian_grads)
        for (unsigned int i = 0; i < output_data.jacobian_grads.size(); ++i)
          output_data.jacobian_grads[i] = DerivativeForm<2, dim, spacedim>();

      if (data.update_each & update_jacobian_pushed_forward_grads)
        for (unsigned int i = 0;
             i < output_data.jacobian_pushed_forward_grads.size();
             ++i)
          output_data.jacobian_pushed_forward_grads[i] = Tensor<3, spacedim>();

      if (data.update_each & update_jacobian_2nd_derivatives)
        for (unsigned int i = 0;
             i < output_data.jacobian_2nd_derivatives.size();
             ++i)
          output_data.jacobian_2nd_derivatives[i] =
            DerivativeForm<3, dim, spacedim>();

      if (data.update_each & update_jacobian_pushed_forward_2nd_derivatives)
        for (unsigned int i = 0;
             i < output_data.jacobian_pushed_forward_2nd_derivatives.size();
             ++i)
          output_data.jacobian_pushed_forward_2nd_derivatives[i] =
            Tensor<4, spacedim>();

      if (data.update_each & update_jacobian_3rd_derivatives)
        for (unsigned int i = 0;
             i < output_data.jacobian_3rd_derivatives.size();
             ++i)
          output_data.jacobian_3rd_derivatives[i] =
            DerivativeForm<4, dim, spacedim>();

      if (data.update_each & update_jacobian_pushed_forward_3rd_derivatives)
        for (unsigned int i = 0;
             i < output_data.jacobian_pushed_forward_3rd_derivatives.size();
             ++i)
          output_data.jacobian_pushed_forward_3rd_derivatives[i] =
            Tensor<5, spacedim>();
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::maybe_update_volume_elements(
  const InternalData &data) const
{
  if (data.update_each & update_volume_elements)
    {
      double volume = data.cell_extents[0];
      for (unsigned int d = 1; d < dim; ++d)
        volume *= data.cell_extents[d];
      data.volume_element = volume;
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::maybe_update_jacobians(
  const InternalData              &data,
  const CellSimilarity::Similarity cell_similarity,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // "compute" Jacobian at the quadrature points, which are all the
  // same
  if (data.update_each & update_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int i = 0; i < output_data.jacobians.size(); ++i)
        {
          output_data.jacobians[i] = DerivativeForm<1, dim, spacedim>();
          for (unsigned int j = 0; j < dim; ++j)
            output_data.jacobians[i][j][j] = data.cell_extents[j];
        }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::maybe_update_inverse_jacobians(
  const InternalData              &data,
  const CellSimilarity::Similarity cell_similarity,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // "compute" inverse Jacobian at the quadrature points, which are
  // all the same
  if (data.update_each & update_inverse_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int i = 0; i < output_data.inverse_jacobians.size(); ++i)
        {
          output_data.inverse_jacobians[i] = Tensor<2, dim>();
          for (unsigned int j = 0; j < dim; ++j)
            output_data.inverse_jacobians[i][j][j] =
              data.inverse_cell_extents[j];
        }
}



template <int dim, int spacedim>
CellSimilarity::Similarity
MappingBox<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim>                                      &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(has_box(cell, polytope_translator), ExcCellNotAssociatedWithBox());

  // convert data object to internal data for this class. fails with
  // an exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);


  update_cell_extents(cell, cell_similarity, data);

  maybe_update_cell_quadrature_points(cell,
                                      data,
                                      quadrature.get_points(),
                                      output_data.quadrature_points);

  // compute Jacobian determinant. all values are equal and are the
  // product of the local lengths in each coordinate direction
  if (data.update_each & (update_JxW_values | update_volume_elements))
    if (cell_similarity != CellSimilarity::translation)
      {
        double J = data.cell_extents[0];
        for (unsigned int d = 1; d < dim; ++d)
          J *= data.cell_extents[d];
        data.volume_element = J;
        if (data.update_each & update_JxW_values)
          for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
            output_data.JxW_values[i] = quadrature.weight(i);
      }


  maybe_update_jacobians(data, cell_similarity, output_data);
  maybe_update_jacobian_derivatives(data, cell_similarity, output_data);
  maybe_update_inverse_jacobians(data, cell_similarity, output_data);

  return cell_similarity;
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  (void)cell;
  (void)face_no;
  (void)subface_no;
  (void)quadrature;
  (void)internal_data;
  (void)output_data;
  DEAL_II_NOT_IMPLEMENTED();
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::fill_fe_immersed_surface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const NonMatching::ImmersedSurfaceQuadrature<dim>          &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  AssertDimension(dim, spacedim);
  Assert(has_box(cell, polytope_translator), ExcCellNotAssociatedWithBox());

  // Convert data object to internal data for this class. Fails with an
  // exception if that is not possible.
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);


  update_cell_extents(cell, CellSimilarity::none, data);

  maybe_update_cell_quadrature_points(cell,
                                      data,
                                      quadrature.get_points(),
                                      output_data.quadrature_points);

  if (data.update_each & update_normal_vectors)
    for (unsigned int i = 0; i < output_data.normal_vectors.size(); ++i)
      output_data.normal_vectors[i] = quadrature.normal_vector(i);

  if (data.update_each & update_JxW_values)
    for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
      output_data.JxW_values[i] = quadrature.weight(i);

  maybe_update_volume_elements(data);
  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::transform(
  const ArrayView<const Tensor<1, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<1, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_kind)
    {
      case mapping_covariant:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d = 0; d < dim; ++d)
              output[i][d] = input[i][d] * data.inverse_cell_extents[d];
          return;
        }

      case mapping_contravariant:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d = 0; d < dim; ++d)
              output[i][d] = input[i][d] * data.cell_extents[d];
          return;
        }
      case mapping_piola:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));
          Assert(data.update_each & update_volume_elements,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_volume_elements"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d = 0; d < dim; ++d)
              output[i][d] =
                input[i][d] * data.cell_extents[d] / data.volume_element;
          return;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_kind)
    {
      case mapping_covariant:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] =
                  input[i][d1][d2] * data.inverse_cell_extents[d2];
          return;
        }

      case mapping_contravariant:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2];
          return;
        }

      case mapping_covariant_gradient:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] *
                                    data.inverse_cell_extents[d2] *
                                    data.inverse_cell_extents[d1];
          return;
        }

      case mapping_contravariant_gradient:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] *
                                    data.inverse_cell_extents[d1];
          return;
        }

      case mapping_piola:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));
          Assert(data.update_each & update_volume_elements,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_volume_elements"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] /
                                    data.volume_element;
          return;
        }

      case mapping_piola_gradient:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));
          Assert(data.update_each & update_volume_elements,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_volume_elements"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] *
                                    data.inverse_cell_extents[d1] /
                                    data.volume_element;
          return;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::transform(
  const ArrayView<const Tensor<2, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_kind)
    {
      case mapping_covariant:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] =
                  input[i][d1][d2] * data.inverse_cell_extents[d2];
          return;
        }

      case mapping_contravariant:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2];
          return;
        }

      case mapping_covariant_gradient:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] *
                                    data.inverse_cell_extents[d2] *
                                    data.inverse_cell_extents[d1];
          return;
        }

      case mapping_contravariant_gradient:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] *
                                    data.inverse_cell_extents[d1];
          return;
        }

      case mapping_piola:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));
          Assert(data.update_each & update_volume_elements,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_volume_elements"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] /
                                    data.volume_element;
          return;
        }

      case mapping_piola_gradient:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));
          Assert(data.update_each & update_volume_elements,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_volume_elements"));

          for (unsigned int i = 0; i < output.size(); ++i)
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] *
                                    data.inverse_cell_extents[d1] /
                                    data.volume_element;
          return;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_kind)
    {
      case mapping_covariant_gradient:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int q = 0; q < output.size(); ++q)
            for (unsigned int i = 0; i < spacedim; ++i)
              for (unsigned int j = 0; j < spacedim; ++j)
                for (unsigned int k = 0; k < spacedim; ++k)
                  {
                    output[q][i][j][k] = input[q][i][j][k] *
                                         data.inverse_cell_extents[j] *
                                         data.inverse_cell_extents[k];
                  }
          return;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::transform(
  const ArrayView<const Tensor<3, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_kind)
    {
      case mapping_contravariant_hessian:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int q = 0; q < output.size(); ++q)
            for (unsigned int i = 0; i < spacedim; ++i)
              for (unsigned int j = 0; j < spacedim; ++j)
                for (unsigned int k = 0; k < spacedim; ++k)
                  {
                    output[q][i][j][k] = input[q][i][j][k] *
                                         data.cell_extents[i] *
                                         data.inverse_cell_extents[j] *
                                         data.inverse_cell_extents[k];
                  }
          return;
        }

      case mapping_covariant_hessian:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int q = 0; q < output.size(); ++q)
            for (unsigned int i = 0; i < spacedim; ++i)
              for (unsigned int j = 0; j < spacedim; ++j)
                for (unsigned int k = 0; k < spacedim; ++k)
                  {
                    output[q][i][j][k] = input[q][i][j][k] *
                                         (data.inverse_cell_extents[i] *
                                          data.inverse_cell_extents[j]) *
                                         data.inverse_cell_extents[k];
                  }

          return;
        }

      case mapping_piola_hessian:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));
          Assert(data.update_each & update_volume_elements,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_volume_elements"));

          for (unsigned int q = 0; q < output.size(); ++q)
            for (unsigned int i = 0; i < spacedim; ++i)
              for (unsigned int j = 0; j < spacedim; ++j)
                for (unsigned int k = 0; k < spacedim; ++k)
                  {
                    output[q][i][j][k] =
                      input[q][i][j][k] *
                      (data.cell_extents[i] / data.volume_element *
                       data.inverse_cell_extents[j]) *
                      data.inverse_cell_extents[k];
                  }

          return;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
Point<spacedim>
MappingBox<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &p) const
{
  Assert(has_box(cell, polytope_translator), ExcCellNotAssociatedWithBox());
  Assert(dim == spacedim, ExcNotImplemented());

  return boxes[polytope_translator.at(cell->active_cell_index())].unit_to_real(
    p);
}



template <int dim, int spacedim>
Point<dim>
MappingBox<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim>                                      &p) const
{
  Assert(has_box(cell, polytope_translator), ExcCellNotAssociatedWithBox());
  Assert(dim == spacedim, ExcNotImplemented());

  return boxes[polytope_translator.at(cell->active_cell_index())].real_to_unit(
    p);
}



template <int dim, int spacedim>
void
MappingBox<dim, spacedim>::transform_points_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<spacedim>>                     &real_points,
  const ArrayView<Point<dim>>                                &unit_points) const
{
  Assert(has_box(cell, polytope_translator), ExcCellNotAssociatedWithBox());
  AssertDimension(real_points.size(), unit_points.size());

  if (dim != spacedim)
    DEAL_II_NOT_IMPLEMENTED();
  for (unsigned int i = 0; i < real_points.size(); ++i)
    unit_points[i] =
      boxes[polytope_translator.at(cell->active_cell_index())].real_to_unit(
        real_points[i]);
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingBox<dim, spacedim>::clone() const
{
  return std::make_unique<MappingBox<dim, spacedim>>(*this);
}


//---------------------------------------------------------------------------
// explicit instantiations
#include "fe/mapping_box.inst"


DEAL_II_NAMESPACE_CLOSE
