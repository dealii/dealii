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
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <memory>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
MappingQ<dim, spacedim>::InternalData::InternalData()
  : use_mapping_q1_on_current_cell(false)
{}



template <int dim, int spacedim>
std::size_t
MappingQ<dim, spacedim>::InternalData::memory_consumption() const
{
  return (
    Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
    MemoryConsumption::memory_consumption(use_mapping_q1_on_current_cell) +
    MemoryConsumption::memory_consumption(mapping_q1_data) +
    MemoryConsumption::memory_consumption(mapping_qp_data));
}



template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const unsigned int degree,
                                  const bool         use_mapping_q_on_all_cells)
  : polynomial_degree(degree)
  ,

  // see whether we want to use *this* mapping objects on *all* cells,
  // or defer to an explicit Q1 mapping on interior cells. if
  // degree==1, then we are already that Q1 mapping, so we don't need
  // it; if dim!=spacedim, there is also no need for anything because
  // we're most likely on a curved manifold
  use_mapping_q_on_all_cells(degree == 1 || use_mapping_q_on_all_cells ||
                             (dim != spacedim))
  ,
  // create a Q1 mapping for use on interior cells (if necessary)
  // or to create a good initial guess in transform_real_to_unit_cell()
  q1_mapping(std::make_shared<MappingQGeneric<dim, spacedim>>(1))
  ,

  // create a Q_p mapping; if p=1, simply share the Q_1 mapping already
  // created via the shared_ptr objects
  qp_mapping(this->polynomial_degree > 1 ?
               std::make_shared<const MappingQGeneric<dim, spacedim>>(degree) :
               q1_mapping)
{}



template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const MappingQ<dim, spacedim> &mapping)
  : polynomial_degree(mapping.polynomial_degree)
  , use_mapping_q_on_all_cells(mapping.use_mapping_q_on_all_cells)
{
  // Note that we really do have to use clone() here, since mapping.q1_mapping
  // may be MappingQ1Eulerian and mapping.qp_mapping may be MappingQEulerian.
  std::shared_ptr<const Mapping<dim, spacedim>> other_q1_map =
    mapping.q1_mapping->clone();
  q1_mapping = std::dynamic_pointer_cast<const MappingQGeneric<dim, spacedim>>(
    other_q1_map);
  Assert(q1_mapping != nullptr, ExcInternalError());
  Assert(q1_mapping->get_degree() == 1, ExcInternalError());

  // Same as the other constructor: if possible reuse the Q1 mapping
  if (this->polynomial_degree == 1)
    {
      qp_mapping = q1_mapping;
    }
  else
    {
      std::shared_ptr<const Mapping<dim, spacedim>> other_qp_map =
        mapping.qp_mapping->clone();
      qp_mapping =
        std::dynamic_pointer_cast<const MappingQGeneric<dim, spacedim>>(
          other_qp_map);
      Assert(qp_mapping != nullptr, ExcInternalError());
    }
}



template <int dim, int spacedim>
unsigned int
MappingQ<dim, spacedim>::get_degree() const
{
  return polynomial_degree;
}



template <int dim, int spacedim>
inline bool
MappingQ<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



template <int dim, int spacedim>
UpdateFlags
MappingQ<dim, spacedim>::requires_update_flags(const UpdateFlags in) const
{
  return (q1_mapping->requires_update_flags(in) |
          qp_mapping->requires_update_flags(in));
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingQ<dim, spacedim>::get_data(const UpdateFlags      update_flags,
                                  const Quadrature<dim> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>();
  auto &data = dynamic_cast<InternalData &>(*data_ptr);

  // build the Q1 and Qp internal data objects in parallel
  Threads::Task<
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>>
    do_get_data = Threads::new_task(&MappingQGeneric<dim, spacedim>::get_data,
                                    *qp_mapping,
                                    update_flags,
                                    quadrature);

  if (!use_mapping_q_on_all_cells)
    data.mapping_q1_data = Utilities::dynamic_unique_cast<
      typename MappingQGeneric<dim, spacedim>::InternalData>(
      std::move(q1_mapping->get_data(update_flags, quadrature)));

  // wait for the task above to finish and use returned value
  data.mapping_qp_data = Utilities::dynamic_unique_cast<
    typename MappingQGeneric<dim, spacedim>::InternalData>(
    std::move(do_get_data.return_value()));
  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingQ<dim, spacedim>::get_face_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>();
  auto &data = dynamic_cast<InternalData &>(*data_ptr);

  // build the Q1 and Qp internal data objects in parallel
  Threads::Task<
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>>
    do_get_data =
      Threads::new_task(&MappingQGeneric<dim, spacedim>::get_face_data,
                        *qp_mapping,
                        update_flags,
                        quadrature);

  if (!use_mapping_q_on_all_cells)
    data.mapping_q1_data = Utilities::dynamic_unique_cast<
      typename MappingQGeneric<dim, spacedim>::InternalData>(
      std::move(q1_mapping->get_face_data(update_flags, quadrature)));

  // wait for the task above to finish and use returned value
  data.mapping_qp_data = Utilities::dynamic_unique_cast<
    typename MappingQGeneric<dim, spacedim>::InternalData>(
    std::move(do_get_data.return_value()));
  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingQ<dim, spacedim>::get_subface_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>();
  auto &data = dynamic_cast<InternalData &>(*data_ptr);

  // build the Q1 and Qp internal data objects in parallel
  Threads::Task<
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>>
    do_get_data =
      Threads::new_task(&MappingQGeneric<dim, spacedim>::get_subface_data,
                        *qp_mapping,
                        update_flags,
                        quadrature);

  if (!use_mapping_q_on_all_cells)
    data.mapping_q1_data = Utilities::dynamic_unique_cast<
      typename MappingQGeneric<dim, spacedim>::InternalData>(
      std::move(q1_mapping->get_subface_data(update_flags, quadrature)));

  // wait for the task above to finish and use returned value
  data.mapping_qp_data = Utilities::dynamic_unique_cast<
    typename MappingQGeneric<dim, spacedim>::InternalData>(
    std::move(do_get_data.return_value()));
  return data_ptr;
}


// Note that the CellSimilarity flag is modifiable, since MappingQ can need to
// recalculate data even when cells are similar.
template <int dim, int spacedim>
CellSimilarity::Similarity
MappingQ<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim> &                                     quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is in the interior of the domain
  data.use_mapping_q1_on_current_cell =
    !(use_mapping_q_on_all_cells || cell->has_boundary_lines());


  // call the base class. we need to ensure that the flag indicating whether
  // we can use some similarity has to be modified - for a general MappingQ,
  // the data needs to be recomputed anyway since then the mapping changes the
  // data. this needs to be known also for later operations, so modify the
  // variable here. this also affects the calculation of the next cell -- if
  // we use Q1 data on the next cell, the data will still be invalid.
  const CellSimilarity::Similarity updated_cell_similarity =
    ((data.use_mapping_q1_on_current_cell == false) &&
         (this->polynomial_degree > 1) ?
       CellSimilarity::invalid_next_cell :
       cell_similarity);

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_values(cell,
                               updated_cell_similarity,
                               quadrature,
                               *data.mapping_q1_data,
                               output_data);
  else
    qp_mapping->fill_fe_values(cell,
                               updated_cell_similarity,
                               quadrature,
                               *data.mapping_qp_data,
                               output_data);

  return updated_cell_similarity;
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell =
    !(use_mapping_q_on_all_cells || cell->has_boundary_lines());

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_face_values(
      cell, face_no, quadrature, *data.mapping_q1_data, output_data);
  else
    qp_mapping->fill_fe_face_values(
      cell, face_no, quadrature, *data.mapping_qp_data, output_data);
}


template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell =
    !(use_mapping_q_on_all_cells || cell->has_boundary_lines());

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_subface_values(cell,
                                       face_no,
                                       subface_no,
                                       quadrature,
                                       *data.mapping_q1_data,
                                       output_data);
  else
    qp_mapping->fill_fe_subface_values(cell,
                                       face_no,
                                       subface_no,
                                       quadrature,
                                       *data.mapping_qp_data,
                                       output_data);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const Tensor<1, dim>> &                  input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<1, spacedim>> &                   output) const
{
  AssertDimension(input.size(), output.size());

  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);
  Assert(data != nullptr, ExcInternalError());

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    q1_mapping->transform(input, mapping_kind, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    qp_mapping->transform(input, mapping_kind, *data->mapping_qp_data, output);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>> &                   output) const
{
  AssertDimension(input.size(), output.size());
  Assert((dynamic_cast<const typename MappingQ<dim, spacedim>::InternalData *>(
            &mapping_data) != nullptr),
         ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    q1_mapping->transform(input, mapping_kind, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    qp_mapping->transform(input, mapping_kind, *data->mapping_qp_data, output);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const Tensor<2, dim>> &                  input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>> &                   output) const
{
  AssertDimension(input.size(), output.size());
  Assert((dynamic_cast<const typename MappingQ<dim, spacedim>::InternalData *>(
            &mapping_data) != nullptr),
         ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    q1_mapping->transform(input, mapping_kind, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    qp_mapping->transform(input, mapping_kind, *data->mapping_qp_data, output);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>> &                   output) const
{
  AssertDimension(input.size(), output.size());
  Assert((dynamic_cast<const typename MappingQ<dim, spacedim>::InternalData *>(
            &mapping_data) != nullptr),
         ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    q1_mapping->transform(input, mapping_kind, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    qp_mapping->transform(input, mapping_kind, *data->mapping_qp_data, output);
}



template <int dim, int spacedim>
void
MappingQ<dim, spacedim>::transform(
  const ArrayView<const Tensor<3, dim>> &                  input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>> &                   output) const
{
  AssertDimension(input.size(), output.size());
  Assert((dynamic_cast<const typename MappingQ<dim, spacedim>::InternalData *>(
            &mapping_data) != nullptr),
         ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    q1_mapping->transform(input, mapping_kind, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    qp_mapping->transform(input, mapping_kind, *data->mapping_qp_data, output);
}



template <int dim, int spacedim>
Point<spacedim>
MappingQ<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim> &                                          p) const
{
  // first see, whether we want to use a linear or a higher order
  // mapping, then either use our own facilities or that of the Q1
  // mapping we store
  if (use_mapping_q_on_all_cells || cell->has_boundary_lines())
    return qp_mapping->transform_unit_to_real_cell(cell, p);
  else
    return q1_mapping->transform_unit_to_real_cell(cell, p);
}



template <int dim, int spacedim>
Point<dim>
MappingQ<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim> &                                     p) const
{
  if (cell->has_boundary_lines() || use_mapping_q_on_all_cells ||
      (dim != spacedim))
    return qp_mapping->transform_real_to_unit_cell(cell, p);
  else
    return q1_mapping->transform_real_to_unit_cell(cell, p);
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ<dim, spacedim>::clone() const
{
  return std::make_unique<MappingQ<dim, spacedim>>(
    this->polynomial_degree, this->use_mapping_q_on_all_cells);
}



// explicit instantiations
#include "mapping_q.inst"


DEAL_II_NAMESPACE_CLOSE
