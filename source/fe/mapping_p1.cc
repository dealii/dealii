// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
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
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_internal.h>
#include <deal.II/fe/mapping_p1.h>

#include <deal.II/grid/tria_iterator.h>

#include <algorithm>
#include <cmath>
#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace
{
  template <int dim, int spacedim>
  DerivativeForm<1, dim, spacedim>
  compute_linear_transformation(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  {
    DerivativeForm<1, dim, spacedim> result;
    for (unsigned int j = 0; j < spacedim; ++j)
      for (unsigned int i = 1; i < dim + 1; ++i)
        result[j][i - 1] = cell->vertex(i)[j] - cell->vertex(0)[j];
    return result;
  }
} // namespace



template <int dim, int spacedim>
MappingP1<dim, spacedim>::InternalData::InternalData(
  const ArrayView<const Point<dim>> &quadrature_points)
{
  quadrature.initialize(quadrature_points);
}



template <int dim, int spacedim>
MappingP1<dim, spacedim>::InternalData::InternalData(
  const Quadrature<dim> &quadrature)
  : quadrature(quadrature)
{}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::InternalData::reinit(
  const UpdateFlags      update_flags,
  const Quadrature<dim> &quadrature)
{
  this->quadrature  = quadrature;
  this->update_each = update_flags;
}



template <int dim, int spacedim>
std::size_t
MappingP1<dim, spacedim>::InternalData::memory_consumption() const
{
  return (Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption(affine_component) +
          MemoryConsumption::memory_consumption(linear_component) +
          MemoryConsumption::memory_consumption(covariant) +
          MemoryConsumption::memory_consumption(determinant) +
          MemoryConsumption::memory_consumption(quadrature));
}



template <int dim, int spacedim>
bool
MappingP1<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



template <int dim, int spacedim>
bool
MappingP1<dim, spacedim>::is_compatible_with(
  const ReferenceCell &reference_cell) const
{
  return reference_cell.is_simplex();
}



template <int dim, int spacedim>
UpdateFlags
MappingP1<dim, spacedim>::requires_update_flags(const UpdateFlags in) const
{
  // Like MappingCartesian, this mapping is simple and has minimal update
  // interdependencies
  UpdateFlags out = in;
  if (out & update_boundary_forms)
    out |= update_normal_vectors;

  return out;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingP1<dim, spacedim>::get_data(const UpdateFlags      update_flags,
                                   const Quadrature<dim> &quadrature) const
{
  auto data_ptr = std::make_unique<InternalData>(quadrature);

  // verify that we have computed the transitive hull of the required
  // flags and that FEValues has faithfully passed them on to us
  Assert(update_flags == requires_update_flags(update_flags),
         ExcInternalError());

  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values(). use the transitive hull of the required
  // flags
  data_ptr->update_each = requires_update_flags(update_flags);

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingP1<dim, spacedim>::get_face_data(
  const UpdateFlags               update_flags,
  const hp::QCollection<dim - 1> &quadrature) const
{
  auto data_ptr = std::make_unique<InternalData>(
    QProjector<dim>::project_to_all_faces(ReferenceCells::get_simplex<dim>(),
                                          quadrature));

  // verify that we have computed the transitive hull of the required
  // flags and that FEValues has faithfully passed them on to us
  Assert(update_flags == requires_update_flags(update_flags),
         ExcInternalError());

  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  data_ptr->update_each = update_flags;

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingP1<dim, spacedim>::get_subface_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  auto data_ptr = std::make_unique<InternalData>(
    QProjector<dim>::project_to_all_subfaces(ReferenceCells::get_simplex<dim>(),
                                             quadrature));

  // verify that we have computed the transitive hull of the required
  // flags and that FEValues has faithfully passed them on to us
  Assert(update_flags == requires_update_flags(update_flags),
         ExcInternalError());

  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  data_ptr->update_each = update_flags;

  return data_ptr;
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::update_transformation(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const InternalData                                         &data) const
{
  data.affine_component = cell->vertex(0);
  data.linear_component = compute_linear_transformation<dim, spacedim>(cell);
  data.covariant        = data.linear_component.covariant_form();
  data.determinant      = data.linear_component.determinant();
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const InternalData                                         &data,
  const typename QProjector<dim>::DataSetDescriptor          &offset,
  std::vector<Point<spacedim>> &quadrature_points) const
{
  Assert(cell->vertex(0) == data.affine_component, ExcInternalError());
  for (unsigned int i = 0; i < quadrature_points.size(); ++i)
    quadrature_points[i] =
      data.affine_component +
      apply_transformation(data.linear_component,
                           data.quadrature.point(offset + i));
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::maybe_update_normal_vectors(
  const unsigned int                face_no,
  const InternalData               &data,
  std::vector<Tensor<1, spacedim>> &normal_vectors) const
{
  const Tensor<1, dim> ref_normal_vector =
    ReferenceCells::get_simplex<dim>().template unit_normal_vectors<dim>(
      face_no);
  Tensor<1, spacedim> normal_vector =
    apply_transformation(data.covariant, ref_normal_vector);
  normal_vector /= normal_vector.norm();

  std::fill(normal_vectors.begin(), normal_vectors.end(), normal_vector);
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::maybe_update_jacobian_derivatives(
  const InternalData              &data,
  const CellSimilarity::Similarity cell_similarity,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // The Jacobian is constant so its derivatives are zero
  if (cell_similarity != CellSimilarity::translation)
    {
      if (data.update_each & update_jacobian_grads)
        std::fill(output_data.jacobian_grads.begin(),
                  output_data.jacobian_grads.end(),
                  DerivativeForm<2, dim, spacedim>());

      if (data.update_each & update_jacobian_pushed_forward_grads)
        std::fill(output_data.jacobian_pushed_forward_grads.begin(),
                  output_data.jacobian_pushed_forward_grads.end(),
                  Tensor<3, spacedim>());

      if (data.update_each & update_jacobian_2nd_derivatives)
        std::fill(output_data.jacobian_2nd_derivatives.begin(),
                  output_data.jacobian_2nd_derivatives.end(),
                  DerivativeForm<3, dim, spacedim>());

      if (data.update_each & update_jacobian_pushed_forward_2nd_derivatives)
        std::fill(output_data.jacobian_pushed_forward_2nd_derivatives.begin(),
                  output_data.jacobian_pushed_forward_2nd_derivatives.end(),
                  Tensor<4, spacedim>());

      if (data.update_each & update_jacobian_3rd_derivatives)
        std::fill(output_data.jacobian_3rd_derivatives.begin(),
                  output_data.jacobian_3rd_derivatives.end(),
                  DerivativeForm<4, dim, spacedim>());

      if (data.update_each & update_jacobian_pushed_forward_3rd_derivatives)
        std::fill(output_data.jacobian_pushed_forward_3rd_derivatives.begin(),
                  output_data.jacobian_pushed_forward_3rd_derivatives.end(),
                  Tensor<5, spacedim>());
    }
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::maybe_update_jacobians(
  const InternalData              &data,
  const CellSimilarity::Similarity cell_similarity,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // "compute" Jacobian at the quadrature points, which are all the
  // same
  if (data.update_each & update_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      std::fill(output_data.jacobians.begin(),
                output_data.jacobians.end(),
                data.linear_component);
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::maybe_update_inverse_jacobians(
  const InternalData              &data,
  const CellSimilarity::Similarity cell_similarity,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  if (data.update_each & update_inverse_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      {
        const auto inverse = data.covariant.transpose();
        std::fill(output_data.inverse_jacobians.begin(),
                  output_data.inverse_jacobians.end(),
                  inverse);
      }
}



template <int dim, int spacedim>
CellSimilarity::Similarity
MappingP1<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim>                                      &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  update_transformation(cell, data);
  if (data.update_each & update_quadrature_points)
    transform_quadrature_points(cell,
                                data,
                                QProjector<dim>::DataSetDescriptor::cell(),
                                output_data.quadrature_points);

  // A special property of this mapping is that the Jacobian is constant over
  // the cell
  if (data.update_each & update_JxW_values)
    if (cell_similarity != CellSimilarity::translation)
      {
        for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
          output_data.JxW_values[i] = data.determinant * quadrature.weight(i);
      }

  maybe_update_jacobians(data, cell_similarity, output_data);
  maybe_update_jacobian_derivatives(data, cell_similarity, output_data);
  maybe_update_inverse_jacobians(data, cell_similarity, output_data);

  if (data.update_each & update_normal_vectors)
    {
      Assert(spacedim == dim + 1,
             ExcMessage("There is no (unique) cell normal for " +
                        Utilities::int_to_string(dim) +
                        "-dimensional cells in " +
                        Utilities::int_to_string(spacedim) +
                        "-dimensional space. This only works if the "
                        "space dimension is one greater than the "
                        "dimensionality of the mesh cells."));

      Tensor<1, spacedim> normal;
      // avoid warnings by only computing cross products in supported dimensions
      if constexpr (dim == 1 && spacedim == 2)
        normal = cross_product_2d(-data.linear_component.transpose()[0]);
      else if constexpr (dim == 2 && spacedim == 3)
        {
          const auto transpose = data.linear_component.transpose();
          normal               = cross_product_3d(transpose[0], transpose[1]);
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();
      normal /= normal.norm();

      if (cell->direction_flag() == false)
        normal *= -1.0;

      std::fill(output_data.normal_vectors.begin(),
                output_data.normal_vectors.end(),
                normal);
    }

  return cell_similarity;
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::fill_mapping_data_for_generic_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>>                          &unit_points,
  const UpdateFlags                                           update_flags,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  if (update_flags == update_default)
    return;

  Assert(update_flags & update_inverse_jacobians ||
           update_flags & update_jacobians ||
           update_flags & update_quadrature_points,
         ExcNotImplemented());

  output_data.initialize(unit_points.size(), update_flags);

  InternalData data(unit_points);
  data.update_each = update_flags;
  update_transformation(cell, data);

  if (data.update_each & update_quadrature_points)
    transform_quadrature_points(cell,
                                data,
                                QProjector<dim>::DataSetDescriptor::cell(),
                                output_data.quadrature_points);

  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1> &quadrature_collection,
  const typename Mapping<dim, spacedim>::InternalDataBase &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  update_transformation(cell, data);
  const auto offset =
    QProjector<dim>::DataSetDescriptor::face(cell->reference_cell(),
                                             face_no,
                                             cell->combined_face_orientation(
                                               face_no),
                                             quadrature_collection);
  if (data.update_each & update_quadrature_points)
    transform_quadrature_points(cell,
                                data,
                                offset,
                                output_data.quadrature_points);

  maybe_update_normal_vectors(face_no, data, output_data.normal_vectors);

  if (data.update_each & (update_JxW_values | update_boundary_forms))
    {
      // Since the quadrature weights presently sum to
      // cell->reference_cell().face_measure(face_no), we have to rescale so
      // they sum to the area of the face
      const double J = cell->face(face_no)->measure() /
                       cell->reference_cell().face_measure(face_no);
      if (data.update_each & update_JxW_values)
        for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
          output_data.JxW_values[i] = J * data.quadrature.weight(i + offset);

      if (data.update_each & update_boundary_forms)
        for (unsigned int i = 0; i < output_data.boundary_forms.size(); ++i)
          output_data.boundary_forms[i] = J * output_data.normal_vectors[i];
    }

  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  update_transformation(cell, data);
  const auto offset =
    QProjector<dim>::DataSetDescriptor::subface(cell->reference_cell(),
                                                face_no,
                                                subface_no,
                                                cell->combined_face_orientation(
                                                  face_no),
                                                quadrature.size(),
                                                cell->subface_case(face_no));
  if (data.update_each & update_quadrature_points)
    transform_quadrature_points(cell,
                                data,
                                offset,
                                output_data.quadrature_points);

  maybe_update_normal_vectors(face_no, data, output_data.normal_vectors);

  if (data.update_each & (update_JxW_values | update_boundary_forms))
    {
      // Same as fill_fe_face_values()
      const double J =
        cell->face(face_no)->measure() /
        cell->face(face_no)->reference_cell().volume() /
        // TODO: once we support 3d refinement this should be updated to the
        // simplex version of GeometryInfo::subface_ratio()
        (dim == 2 ? 2 : 4);
      if (data.update_each & update_JxW_values)
        for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
          output_data.JxW_values[i] = J * quadrature.weight(i);

      if (data.update_each & update_boundary_forms)
        for (unsigned int i = 0; i < output_data.boundary_forms.size(); ++i)
          output_data.boundary_forms[i] = J * output_data.normal_vectors[i];
    }

  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform(
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
          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] = apply_transformation(data.covariant, input[i]);
          return;
        }

      case mapping_contravariant:
        {
          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] = apply_transformation(data.linear_component, input[i]);
          return;
        }
      case mapping_piola:
        {
          auto transformation = data.linear_component;
          Assert(data.determinant > 0.0, ExcDivideByZero());
          for (unsigned int d = 0; d < spacedim; ++d)
            transformation[d] *= 1.0 / data.determinant;
          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] = apply_transformation(transformation, input[i]);
          return;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform(
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
          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] = apply_transformation(data.covariant, input[i]);

          return;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform(
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
            output[i] = apply_transformation(data.covariant, input[i]);
          return;
        }

      case mapping_contravariant:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] = apply_transformation(data.linear_component, input[i]);
          return;
        }

      case mapping_covariant_gradient:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            {
              const DerivativeForm<1, spacedim, dim> A =
                apply_transformation(data.covariant, transpose(input[i]));
              output[i] = apply_transformation(data.covariant, A.transpose());
            }
          return;
        }

      case mapping_contravariant_gradient:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            {
              const DerivativeForm<1, spacedim, dim> A =
                apply_transformation(data.linear_component,
                                     transpose(input[i]));
              output[i] = apply_transformation(data.covariant, A.transpose());
            }

          return;
        }

      case mapping_piola_gradient:
        {
          Assert(data.update_each & update_contravariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_contravariant_transformation"));

          auto scaled_contravariant = data.linear_component;
          for (unsigned int d = 0; d < spacedim; ++d)
            scaled_contravariant[d] *= 1.0 / data.determinant;
          for (unsigned int i = 0; i < output.size(); ++i)
            {
              const DerivativeForm<1, spacedim, dim> A =
                apply_transformation(data.covariant, input[i]);
              const Tensor<2, spacedim> T =
                apply_transformation(scaled_contravariant, A.transpose());

              output[i] = transpose(T);
            }

          return;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform(
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
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] =
              internal::apply_covariant_gradient(data.covariant, input[i]);

          return;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform(
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

          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] =
              internal::apply_contravariant_hessian(data.covariant,
                                                    data.linear_component,
                                                    input[i]);

          return;
        }

      case mapping_covariant_hessian:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] =
              internal::apply_covariant_hessian(data.covariant, input[i]);

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

          for (unsigned int i = 0; i < output.size(); ++i)
            output[i] = internal::apply_piola_hessian(data.covariant,
                                                      data.linear_component,
                                                      data.determinant,
                                                      input[i]);

          return;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
Point<spacedim>
MappingP1<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &p) const
{
  const DerivativeForm<1, dim, spacedim> linear_component =
    compute_linear_transformation<dim, spacedim>(cell);
  const Tensor<1, spacedim> sheared = apply_transformation(linear_component, p);
  return cell->vertex(0) + sheared;
}



template <int dim, int spacedim>
Point<dim>
MappingP1<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim>                                      &p) const
{
  const DerivativeForm<1, spacedim, dim> linear_component =
    compute_linear_transformation<dim, spacedim>(cell)
      .covariant_form()
      .transpose();
  const Tensor<1, spacedim> offset = cell->vertex(0);
  return Point<dim>(apply_transformation(linear_component, p - offset));
}



template <int dim, int spacedim>
void
MappingP1<dim, spacedim>::transform_points_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<spacedim>>                     &real_points,
  const ArrayView<Point<dim>>                                &unit_points) const
{
  const DerivativeForm<1, spacedim, dim> linear_component =
    compute_linear_transformation<dim, spacedim>(cell)
      .covariant_form()
      .transpose();
  const Tensor<1, spacedim> offset = cell->vertex(0);
  for (unsigned int i = 0; i < real_points.size(); ++i)
    unit_points[i] = Point<dim>(
      apply_transformation(linear_component, real_points[i] - offset));
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingP1<dim, spacedim>::clone() const
{
  return std::make_unique<MappingP1<dim, spacedim>>(*this);
}


//---------------------------------------------------------------------------
// explicit instantiations
#include "fe/mapping_p1.inst"


DEAL_II_NAMESPACE_CLOSE
