// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, class EulerVectorType, int spacedim>
MappingQ1Eulerian<dim, EulerVectorType, spacedim>::
MappingQ1Eulerian (const EulerVectorType  &euler_transform_vectors,
                   const DoFHandler<dim,spacedim> &shiftmap_dof_handler)
  :
  euler_transform_vectors(&euler_transform_vectors),
  shiftmap_dof_handler(&shiftmap_dof_handler)
{}



template <int dim, class EulerVectorType, int spacedim>
void
MappingQ1Eulerian<dim, EulerVectorType, spacedim>::
compute_mapping_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               std::vector<Point<spacedim> > &a) const
{

  // The assertions can not be in the
  // constructor, since this would
  // require to call
  // dof_handler.distribute_dofs(fe)
  // *before* the mapping object is
  // constructed, which is not
  // necessarily what we want.

//TODO: Only one of these two assertions should be relevant
  AssertDimension (spacedim, shiftmap_dof_handler->get_fe().n_dofs_per_vertex());
  AssertDimension(shiftmap_dof_handler->get_fe().n_components(), spacedim);

  AssertDimension (shiftmap_dof_handler->n_dofs(), euler_transform_vectors->size());

  // cast the
  // Triangulation<dim>::cell_iterator
  // into a
  // DoFHandler<dim>::cell_iterator
  // which is necessary for access to
  // DoFCellAccessor::get_dof_values()
  typename DoFHandler<dim,spacedim>::cell_iterator dof_cell (*cell, shiftmap_dof_handler);

  // We require the cell to be active
  // since we can only then get nodal
  // values for the shifts
  Assert (dof_cell->active() == true, ExcInactiveCell());

  // for Q1 elements, the number of
  // support points should equal the
  // number of vertices
  a.resize(GeometryInfo<dim>::vertices_per_cell);

  // now get the values of the shift
  // vectors at the vertices
  Vector<double> mapping_values (shiftmap_dof_handler->get_fe().dofs_per_cell);
  dof_cell->get_dof_values (*euler_transform_vectors, mapping_values);


  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      Point<spacedim> shift_vector;

      // pick out the value of the
      // shift vector at the present
      // vertex. since vertex dofs
      // are always numbered first,
      // we can access them easily
      for (unsigned int j=0; j<spacedim; ++j)
        shift_vector[j] = mapping_values(i*spacedim+j);

      // compute new support point by
      // old (reference) value and
      // added shift
      a[i] = cell->vertex(i) + shift_vector;
    }
}



template <int dim, class EulerVectorType, int spacedim>
Mapping<dim,spacedim> *
MappingQ1Eulerian<dim, EulerVectorType, spacedim>::clone () const
{
  return new MappingQ1Eulerian<dim,EulerVectorType,spacedim>(*this);
}



template<int dim, class EulerVectorType, int spacedim>
void
MappingQ1Eulerian<dim,EulerVectorType,spacedim>::fill_fe_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Quadrature<dim>                                     &q,
  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
  std::vector<Point<spacedim> >                             &quadrature_points,
  std::vector<double>                                       &JxW_values,
  std::vector< DerivativeForm<1,dim,spacedim>   >              &jacobians,
  std::vector<DerivativeForm<2,dim,spacedim>    >   &jacobian_grads,
  std::vector<DerivativeForm<1,spacedim,dim>    >   &inverse_jacobians,
  std::vector<Point<spacedim> >                             &normal_vectors,
  CellSimilarity::Similarity                           &cell_similarity) const
{
  // disable any previously detected
  // similarity and then enter the
  // respective function of the base class.
  cell_similarity = CellSimilarity::invalid_next_cell;
  MappingQ1<dim,spacedim>::fill_fe_values (cell, q, mapping_data,
                                           quadrature_points, JxW_values, jacobians,
                                           jacobian_grads, inverse_jacobians,
                                           normal_vectors, cell_similarity);
}



// explicit instantiations
#include "mapping_q1_eulerian.inst"


DEAL_II_NAMESPACE_CLOSE
