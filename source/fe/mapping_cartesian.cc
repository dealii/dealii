// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2019 by the deal.II authors
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

#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <cmath>


DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
MappingCartesian<dim, spacedim>::InternalData::InternalData(
  const Quadrature<dim> &q)
  : cell_extents(numbers::signaling_nan<Tensor<1, dim>>())
  , volume_element(numbers::signaling_nan<double>())
  , quadrature_points(q.get_points())
{}



template <int dim, int spacedim>
std::size_t
MappingCartesian<dim, spacedim>::InternalData::memory_consumption() const
{
  return (Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption(cell_extents) +
          MemoryConsumption::memory_consumption(volume_element) +
          MemoryConsumption::memory_consumption(quadrature_points));
}



template <int dim, int spacedim>
bool
MappingCartesian<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



template <int dim, int spacedim>
UpdateFlags
MappingCartesian<dim, spacedim>::requires_update_flags(
  const UpdateFlags in) const
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
MappingCartesian<dim, spacedim>::get_data(const UpdateFlags      update_flags,
                                          const Quadrature<dim> &q) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std_cxx14::make_unique<InternalData>(q);
  auto &data = dynamic_cast<InternalData &>(*data_ptr);

  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values(). use the transitive hull of the required
  // flags
  data.update_each = requires_update_flags(update_flags);

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingCartesian<dim, spacedim>::get_face_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std_cxx14::make_unique<InternalData>(
      QProjector<dim>::project_to_all_faces(quadrature));
  auto &data = dynamic_cast<InternalData &>(*data_ptr);

  // verify that we have computed the transitive hull of the required
  // flags and that FEValues has faithfully passed them on to us
  Assert(update_flags == requires_update_flags(update_flags),
         ExcInternalError());

  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  data.update_each = update_flags;

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingCartesian<dim, spacedim>::get_subface_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std_cxx14::make_unique<InternalData>(
      QProjector<dim>::project_to_all_subfaces(quadrature));
  auto &data = dynamic_cast<InternalData &>(*data_ptr);

  // verify that we have computed the transitive hull of the required
  // flags and that FEValues has faithfully passed them on to us
  Assert(update_flags == requires_update_flags(update_flags),
         ExcInternalError());

  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  data.update_each = update_flags;

  return data_ptr;
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::update_cell_extents(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const InternalData &                                        data) const
{
  // Compute start point and sizes
  // along axes.  Strange vertex
  // numbering makes this complicated
  // again.
  if (cell_similarity != CellSimilarity::translation)
    {
      const Point<dim> &start = cell->vertex(0);
      switch (dim)
        {
          case 1:
            data.cell_extents[0] = cell->vertex(1)(0) - start(0);
            break;
          case 2:
            data.cell_extents[0] = cell->vertex(1)(0) - start(0);
            data.cell_extents[1] = cell->vertex(2)(1) - start(1);
            break;
          case 3:
            data.cell_extents[0] = cell->vertex(1)(0) - start(0);
            data.cell_extents[1] = cell->vertex(2)(1) - start(1);
            data.cell_extents[2] = cell->vertex(4)(2) - start(2);
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_cell_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const InternalData &                                        data,
  std::vector<Point<dim>> &quadrature_points) const
{
  if (data.update_each & update_quadrature_points)
    {
      const typename QProjector<dim>::DataSetDescriptor offset =
        QProjector<dim>::DataSetDescriptor::cell();

      transform_quadrature_points(cell, data, offset, quadrature_points);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_face_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const InternalData &                                        data,
  std::vector<Point<dim>> &quadrature_points) const
{
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  if (data.update_each & update_quadrature_points)
    {
      const typename QProjector<dim>::DataSetDescriptor offset =
        QProjector<dim>::DataSetDescriptor::face(face_no,
                                                 cell->face_orientation(
                                                   face_no),
                                                 cell->face_flip(face_no),
                                                 cell->face_rotation(face_no),
                                                 quadrature_points.size());


      transform_quadrature_points(cell, data, offset, quadrature_points);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_subface_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const InternalData &                                        data,
  std::vector<Point<dim>> &quadrature_points) const
{
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange(sub_no, GeometryInfo<dim>::max_children_per_face);
  if (cell->face(face_no)->has_children())
    {
      AssertIndexRange(sub_no, cell->face(face_no)->n_children());
    }

  if (data.update_each & update_quadrature_points)
    {
      const typename QProjector<dim>::DataSetDescriptor offset =
        QProjector<dim>::DataSetDescriptor::subface(
          face_no,
          sub_no,
          cell->face_orientation(face_no),
          cell->face_flip(face_no),
          cell->face_rotation(face_no),
          quadrature_points.size(),
          cell->subface_case(face_no));

      transform_quadrature_points(cell, data, offset, quadrature_points);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const InternalData &                                        data,
  const typename QProjector<dim>::DataSetDescriptor &         offset,
  std::vector<Point<dim>> &quadrature_points) const
{
  // let @p{start} be the origin of a local coordinate system. it is chosen
  // as the (lower) left vertex
  const Point<dim> &start = cell->vertex(0);

  for (unsigned int i = 0; i < quadrature_points.size(); ++i)
    {
      quadrature_points[i] = start;
      for (unsigned int d = 0; d < dim; ++d)
        quadrature_points[i](d) +=
          data.cell_extents[d] * data.quadrature_points[i + offset](d);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_normal_vectors(
  const unsigned int           face_no,
  const InternalData &         data,
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
MappingCartesian<dim, spacedim>::maybe_update_jacobian_derivatives(
  const InternalData &             data,
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
CellSimilarity::Similarity
MappingCartesian<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim> &                                     quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal data for this class. fails with
  // an exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);


  update_cell_extents(cell, cell_similarity, data);

  maybe_update_cell_quadrature_points(cell,
                                      data,
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
            output_data.JxW_values[i] = J * quadrature.weight(i);
      }
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

  maybe_update_jacobian_derivatives(data, cell_similarity, output_data);

  // "compute" inverse Jacobian at the quadrature points, which are
  // all the same
  if (data.update_each & update_inverse_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int i = 0; i < output_data.inverse_jacobians.size(); ++i)
        {
          output_data.inverse_jacobians[i] = Tensor<2, dim>();
          for (unsigned int j = 0; j < dim; ++j)
            output_data.inverse_jacobians[i][j][j] = 1. / data.cell_extents[j];
        }

  return cell_similarity;
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  update_cell_extents(cell, CellSimilarity::none, data);

  maybe_update_face_quadrature_points(cell,
                                      face_no,
                                      data,
                                      output_data.quadrature_points);

  maybe_update_normal_vectors(face_no, data, output_data.normal_vectors);

  // first compute Jacobian determinant, which is simply the product
  // of the local lengths since the jacobian is diagonal
  double J = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    if (d != GeometryInfo<dim>::unit_normal_direction[face_no])
      J *= data.cell_extents[d];

  if (data.update_each & update_JxW_values)
    for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
      output_data.JxW_values[i] = J * quadrature.weight(i);

  if (data.update_each & update_boundary_forms)
    for (unsigned int i = 0; i < output_data.boundary_forms.size(); ++i)
      output_data.boundary_forms[i] = J * output_data.normal_vectors[i];

  if (data.update_each & update_volume_elements)
    {
      J = data.cell_extents[0];
      for (unsigned int d = 1; d < dim; ++d)
        J *= data.cell_extents[d];
      data.volume_element = J;
    }

  if (data.update_each & update_jacobians)
    for (unsigned int i = 0; i < output_data.jacobians.size(); ++i)
      {
        output_data.jacobians[i] = DerivativeForm<1, dim, spacedim>();
        for (unsigned int d = 0; d < dim; ++d)
          output_data.jacobians[i][d][d] = data.cell_extents[d];
      }

  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);

  if (data.update_each & update_inverse_jacobians)
    for (unsigned int i = 0; i < output_data.inverse_jacobians.size(); ++i)
      {
        output_data.inverse_jacobians[i] = DerivativeForm<1, dim, spacedim>();
        for (unsigned int d = 0; d < dim; ++d)
          output_data.inverse_jacobians[i][d][d] = 1. / data.cell_extents[d];
      }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal data for this class. fails with
  // an exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  update_cell_extents(cell, CellSimilarity::none, data);

  maybe_update_subface_quadrature_points(
    cell, face_no, subface_no, data, output_data.quadrature_points);

  maybe_update_normal_vectors(face_no, data, output_data.normal_vectors);

  // first compute Jacobian determinant, which is simply the product
  // of the local lengths since the jacobian is diagonal
  double J = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    if (d != GeometryInfo<dim>::unit_normal_direction[face_no])
      J *= data.cell_extents[d];

  if (data.update_each & update_JxW_values)
    {
      // Here, cell->face(face_no)->n_children() would be the right
      // choice, but unfortunately the current function is also called
      // for faces without children (see tests/fe/mapping.cc). Add
      // following switch to avoid diffs in tests/fe/mapping.OK
      const unsigned int n_subfaces =
        cell->face(face_no)->has_children() ?
          cell->face(face_no)->n_children() :
          GeometryInfo<dim>::max_children_per_face;
      for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
        output_data.JxW_values[i] = J * quadrature.weight(i) / n_subfaces;
    }

  if (data.update_each & update_boundary_forms)
    for (unsigned int i = 0; i < output_data.boundary_forms.size(); ++i)
      output_data.boundary_forms[i] = J * output_data.normal_vectors[i];

  if (data.update_each & update_volume_elements)
    {
      J = data.cell_extents[0];
      for (unsigned int d = 1; d < dim; ++d)
        J *= data.cell_extents[d];
      data.volume_element = J;
    }

  if (data.update_each & update_jacobians)
    for (unsigned int i = 0; i < output_data.jacobians.size(); ++i)
      {
        output_data.jacobians[i] = DerivativeForm<1, dim, spacedim>();
        for (unsigned int d = 0; d < dim; ++d)
          output_data.jacobians[i][d][d] = data.cell_extents[d];
      }

  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);

  if (data.update_each & update_inverse_jacobians)
    for (unsigned int i = 0; i < output_data.inverse_jacobians.size(); ++i)
      {
        output_data.inverse_jacobians[i] = DerivativeForm<1, spacedim, dim>();
        for (unsigned int d = 0; d < dim; ++d)
          output_data.inverse_jacobians[i][d][d] = 1. / data.cell_extents[d];
      }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform(
  const ArrayView<const Tensor<1, dim>> &                  input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<1, spacedim>> &                   output) const
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
              output[i][d] = input[i][d] / data.cell_extents[d];
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
        Assert(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>> &                   output) const
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
                output[i][d1][d2] = input[i][d1][d2] / data.cell_extents[d2];
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
                output[i][d1][d2] = input[i][d1][d2] / data.cell_extents[d2] /
                                    data.cell_extents[d1];
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
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] /
                                    data.cell_extents[d1];
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
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] /
                                    data.cell_extents[d1] / data.volume_element;
          return;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform(
  const ArrayView<const Tensor<2, dim>> &                  input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>> &                   output) const
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
                output[i][d1][d2] = input[i][d1][d2] / data.cell_extents[d2];
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
                output[i][d1][d2] = input[i][d1][d2] / data.cell_extents[d2] /
                                    data.cell_extents[d1];
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
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] /
                                    data.cell_extents[d1];
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
                output[i][d1][d2] = input[i][d1][d2] * data.cell_extents[d2] /
                                    data.cell_extents[d1] / data.volume_element;
          return;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>> &                   output) const
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
                    output[q][i][j][k] = input[q][i][j][k] /
                                         data.cell_extents[j] /
                                         data.cell_extents[k];
                  }
          return;
        }
      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform(
  const ArrayView<const Tensor<3, dim>> &                  input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>> &                   output) const
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
                    output[q][i][j][k] =
                      input[q][i][j][k] * data.cell_extents[i] /
                      data.cell_extents[j] / data.cell_extents[k];
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
                    output[q][i][j][k] =
                      input[q][i][j][k] / data.cell_extents[i] /
                      data.cell_extents[j] / data.cell_extents[k];
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
                      input[q][i][j][k] * data.cell_extents[i] /
                      data.volume_element / data.cell_extents[j] /
                      data.cell_extents[k];
                  }

          return;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
Point<spacedim>
MappingCartesian<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim> &                                          p) const
{
  Tensor<1, dim>   length;
  const Point<dim> start = cell->vertex(0);
  switch (dim)
    {
      case 1:
        length[0] = cell->vertex(1)(0) - start(0);
        break;
      case 2:
        length[0] = cell->vertex(1)(0) - start(0);
        length[1] = cell->vertex(2)(1) - start(1);
        break;
      case 3:
        length[0] = cell->vertex(1)(0) - start(0);
        length[1] = cell->vertex(2)(1) - start(1);
        length[2] = cell->vertex(4)(2) - start(2);
        break;
      default:
        Assert(false, ExcNotImplemented());
    }

  Point<dim> p_real = cell->vertex(0);
  for (unsigned int d = 0; d < dim; ++d)
    p_real(d) += length[d] * p(d);

  return p_real;
}



template <int dim, int spacedim>
Point<dim>
MappingCartesian<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim> &                                     p) const
{
  if (dim != spacedim)
    Assert(false, ExcNotImplemented());
  const Point<dim> &start = cell->vertex(0);
  Point<dim>        real  = p;
  real -= start;

  switch (dim)
    {
      case 1:
        real(0) /= cell->vertex(1)(0) - start(0);
        break;
      case 2:
        real(0) /= cell->vertex(1)(0) - start(0);
        real(1) /= cell->vertex(2)(1) - start(1);
        break;
      case 3:
        real(0) /= cell->vertex(1)(0) - start(0);
        real(1) /= cell->vertex(2)(1) - start(1);
        real(2) /= cell->vertex(4)(2) - start(2);
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
  return real;
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingCartesian<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<MappingCartesian<dim, spacedim>>(*this);
}


//---------------------------------------------------------------------------
// explicit instantiations
#include "mapping_cartesian.inst"


DEAL_II_NAMESPACE_CLOSE
