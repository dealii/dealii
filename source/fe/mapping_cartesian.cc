// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
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
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <cmath>
#include <memory>


DEAL_II_NAMESPACE_OPEN

DeclExceptionMsg(
  ExcCellNotCartesian,
  "You are using MappingCartesian, but the incoming cell is not Cartesian.");



/**
 * Return whether the incoming cell is of Cartesian shape. This is determined by
 * checking if the smallest BoundingBox that encloses the cell has the same
 * vertices as the cell itself.
 */
template <typename CellType>
bool
is_cartesian(const CellType &cell)
{
  if (!cell->reference_cell().is_hyper_cube())
    return false;

  // The tolerances here are somewhat larger than the square of the machine
  // epsilon, because we are going to compare the square of distances (to
  // avoid computing square roots).
  const double abs_tol           = 1e-30;
  const double rel_tol           = 1e-28;
  const auto   bounding_box      = cell->bounding_box();
  const auto  &bounding_vertices = bounding_box.get_boundary_points();
  const auto   bb_diagonal_length_squared =
    bounding_vertices.first.distance_square(bounding_vertices.second);

  for (const unsigned int v : cell->vertex_indices())
    {
      // Choose a tolerance that takes into account both that vertices far
      // away from the origin have only a finite number of digits
      // that are considered correct (an "absolute tolerance"), as well as that
      // vertices are supposed to be close to the corresponding vertices of the
      // bounding box (a tolerance that is "relative" to the size of the cell).
      //
      // We need to do it this way because when a vertex is far away from
      // the origin, computing the difference between two vertices is subject
      // to cancellation.
      const double tolerance = std::max(abs_tol * cell->vertex(v).norm_square(),
                                        rel_tol * bb_diagonal_length_squared);

      if (cell->vertex(v).distance_square(bounding_box.vertex(v)) > tolerance)
        return false;
    }

  return true;
}



template <int dim, int spacedim>
MappingCartesian<dim, spacedim>::InternalData::InternalData(
  const Quadrature<dim> &q)
  : cell_extents(numbers::signaling_nan<Tensor<1, dim>>())
  , inverse_cell_extents(numbers::signaling_nan<Tensor<1, dim>>())
  , volume_element(numbers::signaling_nan<double>())
  , quadrature_points(q.get_points())
{}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::InternalData::reinit(
  const UpdateFlags update_flags,
  const Quadrature<dim> &)
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values(). use the transitive hull of the required
  // flags
  this->update_each = update_flags;
}



template <int dim, int spacedim>
std::size_t
MappingCartesian<dim, spacedim>::InternalData::memory_consumption() const
{
  return (Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption(cell_extents) +
          MemoryConsumption::memory_consumption(inverse_cell_extents) +
          MemoryConsumption::memory_consumption(volume_element));
}



template <int dim, int spacedim>
bool
MappingCartesian<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



template <int dim, int spacedim>
bool
MappingCartesian<dim, spacedim>::is_compatible_with(
  const ReferenceCell &reference_cell) const
{
  Assert(dim == reference_cell.get_dimension(),
         ExcMessage("The dimension of your mapping (" +
                    Utilities::to_string(dim) +
                    ") and the reference cell cell_type (" +
                    Utilities::to_string(reference_cell.get_dimension()) +
                    " ) do not agree."));

  return reference_cell.is_hyper_cube();
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
    std::make_unique<InternalData>();
  data_ptr->reinit(requires_update_flags(update_flags), q);

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingCartesian<dim, spacedim>::get_face_data(
  const UpdateFlags               update_flags,
  const hp::QCollection<dim - 1> &quadrature) const
{
  AssertDimension(quadrature.size(), 1);

  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(QProjector<dim>::project_to_all_faces(
      ReferenceCells::get_hypercube<dim>(), quadrature[0]));
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
    std::make_unique<InternalData>(QProjector<dim>::project_to_all_subfaces(
      ReferenceCells::get_hypercube<dim>(), quadrature));
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
  const InternalData                                         &data) const
{
  // Compute start point and sizes along axes. The vertices to be looked at
  // are 1, 2, 4 compared to the base vertex 0.
  if (cell_similarity != CellSimilarity::translation)
    {
      const Point<dim> start = cell->vertex(0);
      for (unsigned int d = 0; d < dim; ++d)
        {
          const double cell_extent_d = cell->vertex(1 << d)[d] - start[d];
          data.cell_extents[d]       = cell_extent_d;
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
    const Tensor<1, dim>                               first_vertex,
    const Tensor<1, dim>                               cell_extents,
    const ArrayView<const Point<dim>>                 &unit_quadrature_points,
    const typename QProjector<dim>::DataSetDescriptor &offset,
    std::vector<Point<dim>>                           &quadrature_points)
  {
    for (unsigned int i = 0; i < quadrature_points.size(); ++i)
      {
        quadrature_points[i] = first_vertex;
        for (unsigned int d = 0; d < dim; ++d)
          quadrature_points[i][d] +=
            cell_extents[d] * unit_quadrature_points[i + offset][d];
      }
  }
} // namespace



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_cell_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const InternalData                                         &data,
  const ArrayView<const Point<dim>> &unit_quadrature_points,
  std::vector<Point<dim>>           &quadrature_points) const
{
  if (data.update_each & update_quadrature_points)
    {
      const auto offset = QProjector<dim>::DataSetDescriptor::cell();

      transform_quadrature_points(cell->vertex(0),
                                  data.cell_extents,
                                  unit_quadrature_points,
                                  offset,
                                  quadrature_points);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_face_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const InternalData                                         &data,
  std::vector<Point<dim>> &quadrature_points) const
{
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  if (data.update_each & update_quadrature_points)
    {
      const auto offset = QProjector<dim>::DataSetDescriptor::face(
        ReferenceCells::get_hypercube<dim>(),
        face_no,
        cell->combined_face_orientation(face_no),
        quadrature_points.size());


      transform_quadrature_points(cell->vertex(0),
                                  data.cell_extents,
                                  make_array_view(data.quadrature_points),
                                  offset,
                                  quadrature_points);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_subface_quadrature_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const InternalData                                         &data,
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
      const auto offset = QProjector<dim>::DataSetDescriptor::subface(
        ReferenceCells::get_hypercube<dim>(),
        face_no,
        sub_no,
        cell->combined_face_orientation(face_no),
        quadrature_points.size(),
        cell->subface_case(face_no));

      transform_quadrature_points(cell->vertex(0),
                                  data.cell_extents,
                                  make_array_view(data.quadrature_points),
                                  offset,
                                  quadrature_points);
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_normal_vectors(
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
                ReferenceCells::get_hypercube<dim>()
                  .template face_normal_vector<dim>(face_no));
    }
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::maybe_update_jacobian_derivatives(
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
MappingCartesian<dim, spacedim>::maybe_update_volume_elements(
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
MappingCartesian<dim, spacedim>::maybe_update_jacobians(
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
MappingCartesian<dim, spacedim>::maybe_update_inverse_jacobians(
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
MappingCartesian<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim>                                      &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(is_cartesian(cell), ExcCellNotCartesian());

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
            output_data.JxW_values[i] = J * quadrature.weight(i);
      }


  maybe_update_jacobians(data, cell_similarity, output_data);
  maybe_update_jacobian_derivatives(data, cell_similarity, output_data);
  maybe_update_inverse_jacobians(data, cell_similarity, output_data);

  return cell_similarity;
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_mapping_data_for_generic_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>>                          &unit_points,
  const UpdateFlags                                           update_flags,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  if (update_flags == update_default)
    return;

  Assert(is_cartesian(cell), ExcCellNotCartesian());

  Assert(update_flags & update_inverse_jacobians ||
           update_flags & update_jacobians ||
           update_flags & update_quadrature_points,
         ExcNotImplemented());

  output_data.initialize(unit_points.size(), update_flags);

  InternalData data;
  data.update_each = update_flags;

  update_cell_extents(cell, CellSimilarity::none, data);

  maybe_update_cell_quadrature_points(cell,
                                      data,
                                      unit_points,
                                      output_data.quadrature_points);

  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(is_cartesian(cell), ExcCellNotCartesian());
  AssertDimension(quadrature.size(), 1);

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
      output_data.JxW_values[i] = J * quadrature[0].weight(i);

  if (data.update_each & update_boundary_forms)
    for (unsigned int i = 0; i < output_data.boundary_forms.size(); ++i)
      output_data.boundary_forms[i] = J * output_data.normal_vectors[i];

  maybe_update_volume_elements(data);
  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(is_cartesian(cell), ExcCellNotCartesian());

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

  maybe_update_volume_elements(data);
  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_immersed_surface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const NonMatching::ImmersedSurfaceQuadrature<dim>          &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  AssertDimension(dim, spacedim);
  Assert(is_cartesian(cell), ExcCellNotCartesian());

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
      {
        // The normals are n = J^{-T} * \hat{n} before normalizing.
        Tensor<1, dim>        normal;
        const Tensor<1, dim> &ref_space_normal = quadrature.normal_vector(i);
        for (unsigned int d = 0; d < dim; ++d)
          {
            normal[d] = ref_space_normal[d] * data.inverse_cell_extents[d];
          }
        normal /= normal.norm();
        output_data.normal_vectors[i] = normal;
      }

  if (data.update_each & update_JxW_values)
    for (unsigned int i = 0; i < output_data.JxW_values.size(); ++i)
      {
        const Tensor<1, dim> &ref_space_normal = quadrature.normal_vector(i);

        // J^{-T} \times \hat{n}
        Tensor<1, dim> invJTxNormal;
        double         det_jacobian = 1.;
        for (unsigned int d = 0; d < dim; ++d)
          {
            det_jacobian *= data.cell_extents[d];
            invJTxNormal[d] =
              ref_space_normal[d] * data.inverse_cell_extents[d];
          }
        output_data.JxW_values[i] =
          det_jacobian * invJTxNormal.norm() * quadrature.weight(i);
      }

  maybe_update_volume_elements(data);
  maybe_update_jacobians(data, CellSimilarity::none, output_data);
  maybe_update_jacobian_derivatives(data, CellSimilarity::none, output_data);
  maybe_update_inverse_jacobians(data, CellSimilarity::none, output_data);
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform(
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
MappingCartesian<dim, spacedim>::transform(
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
MappingCartesian<dim, spacedim>::transform(
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
MappingCartesian<dim, spacedim>::transform(
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
MappingCartesian<dim, spacedim>::transform(
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
MappingCartesian<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &p) const
{
  Assert(is_cartesian(cell), ExcCellNotCartesian());
  Assert(dim == spacedim, ExcNotImplemented());

  Point<dim> unit = cell->vertex(0);

  // Go through vertices with numbers 1, 2, 4
  for (unsigned int d = 0; d < dim; ++d)
    unit[d] += (cell->vertex(1 << d)[d] - unit[d]) * p[d];

  return unit;
}



template <int dim, int spacedim>
Point<dim>
MappingCartesian<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim>                                      &p) const
{
  Assert(is_cartesian(cell), ExcCellNotCartesian());
  Assert(dim == spacedim, ExcNotImplemented());

  const Point<dim> start = cell->vertex(0);
  Point<dim>       real  = p;

  // Go through vertices with numbers 1, 2, 4
  for (unsigned int d = 0; d < dim; ++d)
    real[d] = (real[d] - start[d]) / (cell->vertex(1 << d)[d] - start[d]);

  return real;
}



template <int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::transform_points_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<spacedim>>                     &real_points,
  const ArrayView<Point<dim>>                                &unit_points) const
{
  Assert(is_cartesian(cell), ExcCellNotCartesian());
  AssertDimension(real_points.size(), unit_points.size());

  if (dim != spacedim)
    DEAL_II_NOT_IMPLEMENTED();

  const Point<dim> start = cell->vertex(0);

  // Go through vertices with numbers 1, 2, 4
  std::array<double, dim> inverse_lengths;
  for (unsigned int d = 0; d < dim; ++d)
    inverse_lengths[d] = 1. / (cell->vertex(1 << d)[d] - start[d]);

  for (unsigned int i = 0; i < real_points.size(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      unit_points[i][d] = (real_points[i][d] - start[d]) * inverse_lengths[d];
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingCartesian<dim, spacedim>::clone() const
{
  return std::make_unique<MappingCartesian<dim, spacedim>>(*this);
}


//---------------------------------------------------------------------------
// explicit instantiations
#include "fe/mapping_cartesian.inst"


DEAL_II_NAMESPACE_CLOSE
