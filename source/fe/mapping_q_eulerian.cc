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

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q_eulerian.h>

DEAL_II_NAMESPACE_OPEN



// .... MAPPING Q EULERIAN CONSTRUCTOR

template <int dim, class EulerVectorType, int spacedim>
MappingQEulerian<dim, EulerVectorType, spacedim>::
MappingQEulerian (const unsigned int degree,
                  const EulerVectorType &euler_vector,
                  const DoFHandler<dim,spacedim> &euler_dof_handler)
  :
  MappingQ<dim,spacedim>(degree, true),
  euler_vector(&euler_vector),
  euler_dof_handler(&euler_dof_handler),
  support_quadrature(degree),
  fe_values(euler_dof_handler.get_fe(),
            support_quadrature,
            update_values | update_q_points)
{ }



template <int dim, class EulerVectorType, int spacedim>
Mapping<dim,spacedim> *
MappingQEulerian<dim, EulerVectorType, spacedim>::clone () const
{
  return new MappingQEulerian<dim,EulerVectorType,spacedim>(this->get_degree(),
                                                            *euler_vector,
                                                            *euler_dof_handler);
}



// .... SUPPORT QUADRATURE CONSTRUCTOR

template <int dim, class EulerVectorType, int spacedim>
MappingQEulerian<dim,EulerVectorType,spacedim>::
SupportQuadrature::
SupportQuadrature (const unsigned int map_degree)
  :
  Quadrature<dim>(Utilities::fixed_power<dim>(map_degree+1))
{
  // first we determine the support points on the unit cell in lexicographic
  // order, which are (in accordance with MappingQ) the support points of
  // QGaussLobatto.
  const QGaussLobatto<dim> q_iterated(map_degree+1);
  const unsigned int n_q_points = q_iterated.size();

  // we then need to define a renumbering vector that allows us to go from a
  // lexicographic numbering scheme to a hierarchic one.  this fragment is
  // taking almost verbatim from the MappingQ class.
  std::vector<unsigned int> renumber(n_q_points);
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(map_degree-1);

  FETools::lexicographic_to_hierarchic_numbering (
    FiniteElementData<dim> (dpo, 1, map_degree), renumber);

  // finally we assign the quadrature points in the required order.
  for (unsigned int q=0; q<n_q_points; ++q)
    this->quadrature_points[renumber[q]] = q_iterated.point(q);
}



// .... COMPUTE MAPPING SUPPORT POINTS

template <int dim, class EulerVectorType, int spacedim>
void
MappingQEulerian<dim, EulerVectorType, spacedim>::
compute_mapping_support_points
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 std::vector<Point<spacedim> > &a) const
{

  // first, basic assertion with respect to vector size,

  const types::global_dof_index n_dofs  = euler_dof_handler->n_dofs();
  const types::global_dof_index vector_size = euler_vector->size();

  AssertDimension(vector_size,n_dofs);

  // we then transform our tria iterator into a dof iterator so we can access
  // data not associated with triangulations
  typename DoFHandler<dim,spacedim>::cell_iterator dof_cell(*cell, euler_dof_handler);

  Assert (dof_cell->active() == true, ExcInactiveCell());

  // our quadrature rule is chosen so that each quadrature point corresponds
  // to a support point in the undeformed configuration.  we can then query
  // the given displacement field at these points to determine the shift
  // vector that maps the support points to the deformed configuration.

  // we assume that the given field contains dim displacement components, but
  // that there may be other solution components as well (e.g. pressures).
  // this class therefore assumes that the first dim components represent the
  // actual shift vector we need, and simply ignores any components after
  // that.  this implies that the user should order components appropriately,
  // or create a separate dof handler for the displacements.

  const unsigned int n_support_pts = support_quadrature.size();
  const unsigned int n_components  = euler_dof_handler->get_fe().n_components();

  Assert (n_components >= spacedim, ExcDimensionMismatch(n_components, spacedim) );

  std::vector<Vector<double> > shift_vector(n_support_pts,Vector<double>(n_components));

  // fill shift vector for each support point using an fe_values object. make
  // sure that the fe_values variable isn't used simulatenously from different
  // threads
  Threads::Mutex::ScopedLock lock(fe_values_mutex);
  fe_values.reinit(dof_cell);
  fe_values.get_function_values(*euler_vector, shift_vector);

  // and finally compute the positions of the support points in the deformed
  // configuration.
  a.resize(n_support_pts);
  for (unsigned int q=0; q<n_support_pts; ++q)
    {
      a[q] = fe_values.quadrature_point(q);
      for (unsigned int d=0; d<spacedim; ++d)
        a[q](d) += shift_vector[q](d);
    }
}



template<int dim, class EulerVectorType, int spacedim>
void
MappingQEulerian<dim,EulerVectorType,spacedim>::fill_fe_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Quadrature<dim>                                     &q,
  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
  std::vector<Point<spacedim> >                             &quadrature_points,
  std::vector<double>                                       &JxW_values,
  std::vector<DerivativeForm<1,dim,spacedim> >       &jacobians,
  std::vector<DerivativeForm<2,dim,spacedim>  >     &jacobian_grads,
  std::vector<DerivativeForm<1,spacedim,dim>  >     &inverse_jacobians,
  std::vector<Point<spacedim> >                             &normal_vectors,
  CellSimilarity::Similarity                           &cell_similarity) const
{
  // disable any previously detected similarity and hand on to the respective
  // function of the base class.
  cell_similarity = CellSimilarity::invalid_next_cell;
  MappingQ<dim,spacedim>::fill_fe_values (cell, q, mapping_data,
                                          quadrature_points, JxW_values, jacobians,
                                          jacobian_grads, inverse_jacobians,
                                          normal_vectors, cell_similarity);
}



// explicit instantiations
#include "mapping_q_eulerian.inst"


DEAL_II_NAMESPACE_CLOSE
