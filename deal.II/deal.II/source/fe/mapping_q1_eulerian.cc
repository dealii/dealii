//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <fe/mapping_q1_eulerian.h>
#include <lac/vector.h>
#include <lac/petsc_vector.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, class EulerVectorType, int spacedim>
MappingQ1Eulerian<dim, EulerVectorType, spacedim>::
MappingQ1Eulerian (const EulerVectorType  &euler_transform_vectors,
		   const DoFHandler<dim,spacedim> &shiftmap_dof_handler)
                   :
		   euler_transform_vectors(euler_transform_vectors),
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
  Assert (spacedim == shiftmap_dof_handler->get_fe().n_dofs_per_vertex(),
          ExcWrongNoOfComponents());
  Assert (shiftmap_dof_handler->get_fe().n_components() == spacedim,
	  ExcWrongNoOfComponents());

  Assert (shiftmap_dof_handler->n_dofs() == euler_transform_vectors.size(),
          ExcWrongVectorSize(euler_transform_vectors.size(),
			     shiftmap_dof_handler->n_dofs()));

				   // cast the
				   // Triangulation<dim>::cell_iterator
				   // into a
				   // DoFHandler<dim>::cell_iterator
				   // which is necessary for access to
				   // DoFCellAccessor::get_dof_values()
  typename DoFHandler<dim,spacedim>::cell_iterator
    dof_cell (const_cast<Triangulation<dim,spacedim> *> (&(cell->get_triangulation())),
	      cell->level(),
	      cell->index(),
	      shiftmap_dof_handler);

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
  dof_cell->get_dof_values (euler_transform_vectors, mapping_values);


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
  std::vector<Tensor<2,spacedim> >                          &jacobians,
  std::vector<Tensor<3,spacedim> >                          &jacobian_grads,
  std::vector<Tensor<2,spacedim> >                          &inverse_jacobians,
  std::vector<Point<spacedim> >                             &cell_normal_vectors,
  enum CellSimilarity::Similarity                           &cell_similarity) const
{
				   // disable any previously detected
				   // similarity and then enter the
				   // respective function of the base class.
  cell_similarity = CellSimilarity::invalid_next_cell;
  MappingQ1<dim,spacedim>::fill_fe_values (cell, q, mapping_data,
					   quadrature_points, JxW_values, jacobians,
					   jacobian_grads, inverse_jacobians,
					   cell_normal_vectors, cell_similarity);
}



// explicit instantiation
template class MappingQ1Eulerian<deal_II_dimension, Vector<double> >;
#ifdef DEAL_II_USE_PETSC
template class MappingQ1Eulerian<deal_II_dimension, PETScWrappers::Vector>;
#endif

// Explicit instantiation for codimension one problems.
#if deal_II_dimension != 3
template class MappingQ1Eulerian<deal_II_dimension, Vector<double>, deal_II_dimension+1 >;
#	ifdef DEAL_II_USE_PETSC
template class MappingQ1Eulerian<deal_II_dimension, PETScWrappers::Vector, deal_II_dimension+1>;
#	endif
#endif

DEAL_II_NAMESPACE_CLOSE
