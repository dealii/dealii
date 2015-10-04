// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/utilities.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <numeric>

DEAL_II_NAMESPACE_OPEN


template<int dim, int spacedim>
MappingQ<dim,spacedim>::InternalData::InternalData (const unsigned int polynomial_degree)
  :
  MappingQGeneric<dim,spacedim>::InternalData(polynomial_degree),
  use_mapping_q1_on_current_cell(false)
{}



template<int dim, int spacedim>
std::size_t
MappingQ<dim,spacedim>::InternalData::memory_consumption () const
{
  return (MappingQGeneric<dim,spacedim>::InternalData::memory_consumption () +
          MemoryConsumption::memory_consumption (use_mapping_q1_on_current_cell) +
          MemoryConsumption::memory_consumption (mapping_q1_data));
}



template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const unsigned int degree,
                                  const bool use_mapping_q_on_all_cells)
  :
  MappingQGeneric<dim,spacedim>(degree),

  // see whether we want to use *this* mapping objects on *all* cells,
  // or defer to an explicit Q1 mapping on interior cells. if
  // degree==1, then we are already that Q1 mapping, so we don't need
  // it; if dim!=spacedim, there is also no need for anything because
  // we're most likely on a curved manifold
  use_mapping_q_on_all_cells (degree==1
                              ||
                              use_mapping_q_on_all_cells
                              ||
                              (dim != spacedim)),
  // create a Q1 mapping for use on interior cells (if necessary)
  // or to create a good initial guess in transform_real_to_unit_cell()
  q1_mapping (new MappingQ1<dim,spacedim>())
{}



template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const MappingQ<dim,spacedim> &mapping)
  :
  MappingQGeneric<dim,spacedim>(mapping.get_degree()),
  use_mapping_q_on_all_cells (mapping.use_mapping_q_on_all_cells),
  // clone the Q1 mapping for use on interior cells (if necessary)
  // or to create a good initial guess in transform_real_to_unit_cell()
  q1_mapping (mapping.q1_mapping->clone())
{}



template<int dim, int spacedim>
typename MappingQ<dim,spacedim>::InternalData *
MappingQ<dim,spacedim>::get_data (const UpdateFlags update_flags,
                                  const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(this->polynomial_degree);

  // fill the data of both the Q_p and the Q_1 objects in parallel
  Threads::Task<>
  task = Threads::new_task (std_cxx11::function<void()>
                            (std_cxx11::bind(&MappingQGeneric<dim,spacedim>::InternalData::initialize,
                                             data,
                                             update_flags,
                                             std_cxx11::cref(quadrature),
                                             quadrature.size())));
  if (!use_mapping_q_on_all_cells)
    data->mapping_q1_data.reset (q1_mapping->get_data (update_flags, quadrature));

  task.join ();
  return data;
}



template<int dim, int spacedim>
typename MappingQ<dim,spacedim>::InternalData *
MappingQ<dim,spacedim>::get_face_data (const UpdateFlags update_flags,
                                       const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(this->polynomial_degree);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));

  // fill the data of both the Q_p and the Q_1 objects in parallel
  Threads::Task<>
  task = Threads::new_task (std_cxx11::function<void()>
                            (std_cxx11::bind(&MappingQGeneric<dim,spacedim>::InternalData::initialize_face,
                                             data,
                                             update_flags,
                                             std_cxx11::cref(q),
                                             quadrature.size())));
  if (!use_mapping_q_on_all_cells)
    data->mapping_q1_data.reset (q1_mapping->get_face_data (update_flags, quadrature));

  task.join ();
  return data;
}



template<int dim, int spacedim>
typename MappingQ<dim,spacedim>::InternalData *
MappingQ<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
                                          const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(this->polynomial_degree);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));

  // fill the data of both the Q_p and the Q_1 objects in parallel
  Threads::Task<>
  task = Threads::new_task (std_cxx11::function<void()>
                            (std_cxx11::bind(&MappingQGeneric<dim,spacedim>::InternalData::initialize_face,
                                             data,
                                             update_flags,
                                             std_cxx11::cref(q),
                                             quadrature.size())));
  if (!use_mapping_q_on_all_cells)
    data->mapping_q1_data.reset (q1_mapping->get_subface_data (update_flags, quadrature));

  task.join ();
  return data;
}


// Note that the CellSimilarity flag is modifiable, since MappingQ can need to
// recalculate data even when cells are similar.
template<int dim, int spacedim>
CellSimilarity::Similarity
MappingQ<dim,spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const CellSimilarity::Similarity                           cell_similarity,
                const Quadrature<dim>                                     &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0, ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is in the interior of the domain
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());


  // call the base class. we need to ensure that the flag indicating whether
  // we can use some similarity has to be modified - for a general MappingQ,
  // the data needs to be recomputed anyway since then the mapping changes the
  // data. this needs to be known also for later operations, so modify the
  // variable here. this also affects the calculation of the next cell -- if
  // we use Q1 data on the next cell, the data will still be invalid.
  const CellSimilarity::Similarity updated_cell_similarity
    = ((data.use_mapping_q1_on_current_cell == false)
       &&
       (this->polynomial_degree > 1)
       ?
       CellSimilarity::invalid_next_cell
       :
       cell_similarity);

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_values (cell,
                                updated_cell_similarity,
                                quadrature,
                                *data.mapping_q1_data,
                                output_data);
  else
    MappingQGeneric<dim,spacedim>::fill_fe_values(cell,
                                                  updated_cell_similarity,
                                                  quadrature,
                                                  data,
                                                  output_data);

  return updated_cell_similarity;
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int                                         face_no,
                     const Quadrature<dim-1>                                   &quadrature,
                     const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                     internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_face_values (cell,
                                     face_no,
                                     quadrature,
                                     *data.mapping_q1_data,
                                     output_data);
  else
    MappingQGeneric<dim,spacedim>::fill_fe_face_values(cell,
                                                       face_no,
                                                       quadrature,
                                                       data,
                                                       output_data);
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                        const unsigned int                                         face_no,
                        const unsigned int                                         subface_no,
                        const Quadrature<dim-1>                                   &quadrature,
                        const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                        internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_subface_values (cell,
                                        face_no,
                                        subface_no,
                                        quadrature,
                                        *data.mapping_q1_data,
                                        output_data);
  else
    MappingQGeneric<dim,spacedim>::fill_fe_subface_values(cell,
                                                          face_no,
                                                          subface_no,
                                                          quadrature,
                                                          data,
                                                          output_data);
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<Tensor<1,dim> > >   input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           VectorSlice<std::vector<Tensor<1,spacedim> > >          output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQ<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<DerivativeForm<1, dim ,spacedim>  > >  input,
           const MappingType                                                          mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase                    &mapping_data,
           VectorSlice<std::vector<Tensor<2,spacedim> > >                             output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQ<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}


template<int dim, int spacedim>
void MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<Tensor<2, dim> > >  input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           VectorSlice<std::vector<Tensor<2,spacedim> > >          output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQ<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<DerivativeForm<2, dim ,spacedim>  > >  input,
           const MappingType                                                          mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase                    &mapping_data,
           VectorSlice<std::vector<Tensor<3,spacedim> > >                             output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQ<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}


template<int dim, int spacedim>
void MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<Tensor<3, dim> > >  input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           VectorSlice<std::vector<Tensor<3,spacedim> > >          output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQ<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());
  const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data);

  // check whether we should in fact work on the Q1 portion of it
  if (data->use_mapping_q1_on_current_cell)
    return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);
  else
    // otherwise use the full mapping
    MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}


template<int dim, int spacedim>
Point<spacedim>
MappingQ<dim,spacedim>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim>                                 &p) const
{
  // first see, whether we want to use a linear or a higher order
  // mapping, then either use our own facilities or that of the Q1
  // mapping we store
  if (use_mapping_q_on_all_cells || cell->has_boundary_lines())
    return this->MappingQGeneric<dim,spacedim>::transform_unit_to_real_cell (cell, p);
  else
    return q1_mapping->transform_unit_to_real_cell (cell, p);
}



template<int dim, int spacedim>
Point<dim>
MappingQ<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
{
  if (cell->has_boundary_lines()
      ||
      use_mapping_q_on_all_cells
      ||
      (dim!=spacedim) )
    return MappingQGeneric<dim,spacedim>::transform_real_to_unit_cell(cell, p);
  else
    return q1_mapping->transform_real_to_unit_cell(cell, p);
}



template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingQ<dim,spacedim>::clone () const
{
  return new MappingQ<dim,spacedim>(this->polynomial_degree);
}



// explicit instantiations
#include "mapping_q.inst"


DEAL_II_NAMESPACE_CLOSE
