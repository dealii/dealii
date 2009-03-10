//---------------------------------------------------------------------------
//    $Id: mapping_q_eulerian.cc 14038 2006-10-23 02:46:34Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/utilities.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <lac/petsc_vector.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>
#include <fe/mapping_q_eulerian.h>

DEAL_II_NAMESPACE_OPEN



// .... MAPPING Q EULERIAN CONSTRUCTOR

template <int dim, class EulerVectorType, int spacedim>
MappingQEulerian<dim, EulerVectorType, spacedim>::
MappingQEulerian (const unsigned int degree,
		  const EulerVectorType &euler_vector,
		  const DoFHandler<dim> &euler_dof_handler)
		:
		MappingQ<dim,spacedim>(degree, true),
		euler_vector(euler_vector),
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
						   euler_vector,
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
				   // first we determine the support points
				   // on the unit cell in lexicographic order.
				   // for this purpose we can use an interated
				   // trapezoidal quadrature rule.
  const QTrapez<1> q1d;
  const QIterated<dim> q_iterated(q1d,map_degree);
  const unsigned int n_q_points = q_iterated.n_quadrature_points;

				   // we then need to define a renumbering 
				   // vector that allows us to go from a 
				   // lexicographic numbering scheme to a hierarchic
				   // one.  this fragment is taking almost verbatim
				   // from the MappingQ class.

  std::vector<unsigned int> renumber(n_q_points);
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(map_degree-1);

  FETools::lexicographic_to_hierarchic_numbering (
    FiniteElementData<dim> (dpo, 1, map_degree), renumber);

				   // finally we assign the quadrature points in the
				   // required order.

  for(unsigned int q=0; q<n_q_points; ++q)
    this->quadrature_points[renumber[q]] = q_iterated.point(q);
}



// .... COMPUTE MAPPING SUPPORT POINTS

template <int dim, class EulerVectorType, int spacedim>
void
MappingQEulerian<dim, EulerVectorType, spacedim>::
compute_mapping_support_points
(const typename Triangulation<dim>::cell_iterator &cell,
 std::vector<Point<dim> > &a) const
{

				   // first, basic assertion
				   // with respect to vector size,

  const unsigned int n_dofs      = euler_dof_handler->n_dofs();
  const unsigned int vector_size = euler_vector.size();

  Assert (n_dofs == vector_size,ExcWrongVectorSize(vector_size,n_dofs));

				   // we then transform our tria iterator
				   // into a dof iterator so we can 
				   // access data not associated with
				   // triangulations

  typename DoFHandler<dim>::cell_iterator dof_cell 
    (const_cast<Triangulation<dim> *> (&(cell->get_triangulation())),
     cell->level(),
     cell->index(),
     euler_dof_handler);

  Assert (dof_cell->active() == true, ExcInactiveCell());

				   // our quadrature rule is chosen
				   // so that each quadrature point
				   // corresponds to a support point
				   // in the undeformed configuration.
				   // we can then query the given
				   // displacement field at these points
				   // to determine the shift vector that
				   // maps the support points to the 
				   // deformed configuration.

				   // we assume that the given field contains
				   // dim displacement components, but
				   // that there may be other solution
				   // components as well (e.g. pressures).
				   // this class therefore assumes that the
				   // first dim components represent the
				   // actual shift vector we need, and simply
				   // ignores any components after that.  
				   // this implies that the user should order
				   // components appropriately, or create a
				   // separate dof handler for the displacements.

  const unsigned int n_support_pts = support_quadrature.n_quadrature_points;
  const unsigned int n_components  = euler_dof_handler->get_fe().n_components();

  Assert (n_components >= dim, ExcWrongNoOfComponents() );

  std::vector<Vector<double> > shift_vector(n_support_pts,Vector<double>(n_components));

				   // fill shift vector for each
				   // support point using an fe_values
				   // object. make sure that the
				   // fe_values variable isn't used
				   // simulatenously from different
				   // threads
  Threads::ThreadMutex::ScopedLock lock(fe_values_mutex);
  fe_values.reinit(dof_cell);
  fe_values.get_function_values(euler_vector,shift_vector);

				   // and finally compute the positions of the 
				   // support points in the deformed
				   // configuration.

  a.resize(n_support_pts);
  for(unsigned int q=0; q<n_support_pts; ++q)
    {
      a[q] = fe_values.quadrature_point(q);
      for(unsigned int d=0; d<dim; ++d) 
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
  std::vector<Tensor<2,spacedim> >                          &jacobians,
  std::vector<Tensor<3,spacedim> >                          &jacobian_grads,
  std::vector<Tensor<2,spacedim> >                          &inverse_jacobians,
  std::vector<Point<spacedim> >                             &cell_normal_vectors,
  enum CellSimilarity::Similarity                           &cell_similarity) const
{
				   // disable any previously detected
				   // similarity and hand on to the respective
				   // function of the base class.
  cell_similarity = CellSimilarity::invalid_next_cell;
  MappingQ<dim,spacedim>::fill_fe_values (cell, q, mapping_data,
					  quadrature_points, JxW_values, jacobians,
					  jacobian_grads, inverse_jacobians,
					  cell_normal_vectors, cell_similarity);
}



// explicit instantiation
template class MappingQEulerian<deal_II_dimension, Vector<double> >;
#ifdef DEAL_II_USE_PETSC
template class MappingQEulerian<deal_II_dimension, PETScWrappers::Vector>;
#endif

DEAL_II_NAMESPACE_CLOSE
