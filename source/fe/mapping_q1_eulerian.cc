// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <boost/container/small_vector.hpp>

#include <array>
#include <memory>


DEAL_II_NAMESPACE_OPEN


template <int dim, typename VectorType, int spacedim>
MappingQ1Eulerian<dim, VectorType, spacedim>::MappingQ1Eulerian(
  const DoFHandler<dim, spacedim> &shiftmap_dof_handler,
  const VectorType                &euler_transform_vectors)
  : MappingQ1Eulerian(shiftmap_dof_handler,
                      euler_transform_vectors,
                      ComponentMask())
{}



template <int dim, typename VectorType, int spacedim>
MappingQ1Eulerian<dim, VectorType, spacedim>::MappingQ1Eulerian(
  const DoFHandler<dim, spacedim> &shiftmap_dof_handler,
  const VectorType                &euler_transform_vectors,
  const ComponentMask             &component_mask)
  : MappingQ<dim, spacedim>(1)
  , euler_transform_vectors(&euler_transform_vectors)
  , shiftmap_dof_handler(&shiftmap_dof_handler)
  , component_mask(component_mask)
{
  // The finite element may only be supplied later by distribute_dofs().
  if (shiftmap_dof_handler.get_fe_collection().size() != 0)
    get_vertex_dof_indices();
}



template <int dim, typename VectorType, int spacedim>
const std::array<unsigned int, spacedim> &
MappingQ1Eulerian<dim, VectorType, spacedim>::get_vertex_dof_indices() const
{
  std::scoped_lock lock(vertex_dof_indices_mutex);
  const auto      &fe = shiftmap_dof_handler->get_fe();
  if (cached_fe_collection.size() != 0 && &cached_fe_collection[0] == &fe)
    return this->vertex_dof_indices;

  const unsigned int                 n_vertex_dofs = fe.n_dofs_per_vertex();
  std::array<unsigned int, spacedim> vertex_dof_indices{};
  if (component_mask.represents_the_all_selected_mask())
    {
      AssertDimension(n_vertex_dofs, spacedim);
      AssertDimension(fe.n_components(), spacedim);
      for (unsigned int d = 0; d < spacedim; ++d)
        vertex_dof_indices[d] = d;
    }
  else
    {
      AssertThrow(component_mask.size() == fe.n_components(),
                  ExcDimensionMismatch(component_mask.size(),
                                       fe.n_components()));
      AssertThrow(component_mask.n_selected_components() == spacedim,
                  ExcMessage(
                    "Select exactly spacedim displacement components."));
      unsigned int d = 0;
      for (unsigned int c = 0; c < fe.n_components(); ++c)
        if (component_mask[c])
          {
            unsigned int index = numbers::invalid_unsigned_int;
            for (unsigned int j = 0; j < n_vertex_dofs; ++j)
              if (fe.get_nonzero_components(j)[c])
                {
                  AssertThrow(
                    fe.is_primitive(j) &&
                      index == numbers::invalid_unsigned_int,
                    ExcMessage(
                      "Each selected component must have one primitive vertex DoF."));
                  index = j;
                }
            AssertThrow(
              index != numbers::invalid_unsigned_int,
              ExcMessage(
                "Each selected component must have one primitive vertex DoF."));
            // Sampling vertex coefficients is valid only if other selected
            // shape functions do not contribute at vertices (FE_Q_DG0, for
            // example, also has a cell-constant mode).
            for (unsigned int j =
                   fe.reference_cell().n_vertices() * n_vertex_dofs;
                 j < fe.n_dofs_per_cell();
                 ++j)
              if (fe.get_nonzero_components(j)[c])
                for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                  AssertThrow(
                    fe.shape_value_component(
                      j, GeometryInfo<dim>::unit_cell_vertex(v), c) == 0.,
                    ExcMessage(
                      "Selected nonvertex shape functions must vanish at vertices."));
            vertex_dof_indices[d++] = index;
          }
    }

  cached_fe_collection =
    hp::FECollection<dim, spacedim>(shiftmap_dof_handler->get_fe_collection());
  this->vertex_dof_indices = vertex_dof_indices;
  return this->vertex_dof_indices;
}



template <int dim, typename VectorType, int spacedim>
boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                               ReferenceCells::max_n_vertices<dim>()
#else
                               GeometryInfo<dim>::vertices_per_cell
#endif
                               >
MappingQ1Eulerian<dim, VectorType, spacedim>::get_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                                 ReferenceCells::max_n_vertices<dim>()
#else
                                 GeometryInfo<dim>::vertices_per_cell
#endif
                                 >
                     vertices(cell->n_vertices());
  const auto        &vertex_dof_indices = get_vertex_dof_indices();
  const unsigned int n_vertex_dofs =
    shiftmap_dof_handler->get_fe().n_dofs_per_vertex();

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
    shiftmap_dof_handler->get_fe().n_dofs_per_cell());
  dof_cell->get_dof_values(*euler_transform_vectors, mapping_values);

  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    {
      Point<spacedim> shift_vector;

      // pick out the value of the shift vector at the present
      // vertex. since vertex dofs are always numbered first, we can
      // access them easily
      for (unsigned int j = 0; j < spacedim; ++j)
        shift_vector[j] =
          mapping_values(i * n_vertex_dofs + vertex_dof_indices[j]);

      // compute new support point by old (reference) value and added
      // shift
      vertices[i] = cell->vertex(i) + shift_vector;
    }
  return vertices;
}



template <int dim, typename VectorType, int spacedim>
void
MappingQ1Eulerian<dim, VectorType, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  boost::container::small_vector<Point<spacedim>, 200>       &points) const
{
  const auto vertices = this->get_vertices(cell);
  points.assign(vertices.begin(), vertices.end());
}



template <int dim, typename VectorType, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ1Eulerian<dim, VectorType, spacedim>::clone() const
{
  return std::make_unique<MappingQ1Eulerian<dim, VectorType, spacedim>>(
    *shiftmap_dof_handler, *euler_transform_vectors, component_mask);
}



template <int dim, typename VectorType, int spacedim>
CellSimilarity::Similarity
MappingQ1Eulerian<dim, VectorType, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity,
  const Quadrature<dim>                                   &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // call the function of the base class, but ignoring
  // any potentially detected cell similarity between
  // the current and the previous cell
  MappingQ<dim, spacedim>::fill_fe_values(cell,
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
#include "fe/mapping_q1_eulerian.inst"


DEAL_II_NAMESPACE_CLOSE
