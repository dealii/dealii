// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <array>
#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim, class VectorType, int spacedim>
MappingQ1Eulerian<dim, VectorType, spacedim>::MappingQ1Eulerian(
  const DoFHandler<dim, spacedim> &shiftmap_dof_handler,
  const VectorType &               euler_transform_vectors)
  : MappingQGeneric<dim, spacedim>(1)
  , euler_transform_vectors(&euler_transform_vectors)
  , shiftmap_dof_handler(&shiftmap_dof_handler)
{}



template <int dim, class VectorType, int spacedim>
std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
MappingQ1Eulerian<dim, VectorType, spacedim>::get_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell> vertices;
  // The assertions can not be in the constructor, since this would
  // require to call dof_handler.distribute_dofs(fe) *before* the mapping
  // object is constructed, which is not necessarily what we want.

  // TODO: Only one of these two assertions should be relevant
  AssertDimension(spacedim, shiftmap_dof_handler->get_fe().n_dofs_per_vertex());
  AssertDimension(shiftmap_dof_handler->get_fe(0).n_components(), spacedim);

  AssertDimension(shiftmap_dof_handler->n_dofs(),
                  euler_transform_vectors->size());

  // cast the Triangulation<dim>::cell_iterator into a
  // DoFHandler<dim>::cell_iterator which is necessary for access to
  // DoFCellAccessor::get_dof_values()
  typename DoFHandler<dim, spacedim>::cell_iterator dof_cell(
    *cell, shiftmap_dof_handler);

  // We require the cell to be active since we can only then get nodal
  // values for the shifts
  Assert(dof_cell->is_active() == true, ExcInactiveCell());

  // now get the values of the shift vectors at the vertices
  Vector<typename VectorType::value_type> mapping_values(
    shiftmap_dof_handler->get_fe().dofs_per_cell);
  dof_cell->get_dof_values(*euler_transform_vectors, mapping_values);

  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    {
      Point<spacedim> shift_vector;

      // pick out the value of the shift vector at the present
      // vertex. since vertex dofs are always numbered first, we can
      // access them easily
      for (unsigned int j = 0; j < spacedim; ++j)
        shift_vector[j] = mapping_values(i * spacedim + j);

      // compute new support point by old (reference) value and added
      // shift
      vertices[i] = cell->vertex(i) + shift_vector;
    }
  return vertices;
}



template <int dim, class VectorType, int spacedim>
std::vector<Point<spacedim>>
MappingQ1Eulerian<dim, VectorType, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  const std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
    vertices = this->get_vertices(cell);

  std::vector<Point<spacedim>> a(GeometryInfo<dim>::vertices_per_cell);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    a[i] = vertices[i];

  return a;
}



template <int dim, class VectorType, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ1Eulerian<dim, VectorType, spacedim>::clone() const
{
  return std::make_unique<MappingQ1Eulerian<dim, VectorType, spacedim>>(*this);
}



template <int dim, class VectorType, int spacedim>
CellSimilarity::Similarity
MappingQ1Eulerian<dim, VectorType, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity,
  const Quadrature<dim> &                                  quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // call the function of the base class, but ignoring
  // any potentially detected cell similarity between
  // the current and the previous cell
  MappingQGeneric<dim, spacedim>::fill_fe_values(
    cell,
    CellSimilarity::invalid_next_cell,
    quadrature,
    internal_data,
    output_data);
  // also return the updated flag since any detected
  // similarity wasn't based on the mapped field, but
  // the original vertices which are meaningless
  return CellSimilarity::invalid_next_cell;
}



// explicit instantiations
#include "mapping_q1_eulerian.inst"


DEAL_II_NAMESPACE_CLOSE
