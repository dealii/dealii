//----------------------------  mapping_q1_eulerian.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors and Michael Stadler
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_q1_eulerian.cc  ---------------------------


#include <fe/mapping_q1_eulerian.h>
#include <lac/vector.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>



template <int dim>
MappingQ1Eulerian<dim>::MappingQ1Eulerian ( const Vector<double>  &euler_transform_vectors,
					    const DoFHandler<dim> &shiftmap_dof_handler)
                     :
		     euler_transform_vectors(euler_transform_vectors),
		     shiftmap_dof_handler(&shiftmap_dof_handler)
{}



template <int dim>
void
MappingQ1Eulerian<dim>::compute_mapping_support_points(
  const typename Triangulation<dim>::cell_iterator &cell,
  typename std::vector<Point<dim> > &a) const
{

				   // The assertions can not be in the
				   // constructor, since this would
				   // require to call
				   // dof_handler.distribute_dofs(fe)
				   // *before* the mapping object is
				   // constructed, which is not
				   // necessarily what we want.
  Assert (dim == shiftmap_dof_handler->get_fe().n_dofs_per_vertex(),
          ExcWrongNoOfComponents());
  Assert (shiftmap_dof_handler->get_fe().n_components() == dim,
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
  typename DoFHandler<dim>::cell_iterator
    dof_cell (const_cast<Triangulation<dim> *> (&(cell->get_triangulation())),
	      cell->level(),
	      cell->index(),
	      shiftmap_dof_handler);

				   // We require the cell to be
				   // active.  This is determined by
				   // the user when looping over all
				   // active cells in the problem
				   // code.
  Assert (dof_cell->active() == true,
          ExcInactiveCell());

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
      Point<dim> shift_vector;

				       // pick out the value of the
				       // shift vector at the present
				       // vertex. since vertex dofs
				       // are always numbered first,
				       // we can access them easily
      for (unsigned int j=0; j<dim; ++j)
	shift_vector[j] = mapping_values(vertex_mapping[i]*dim+j);

				       // compute new support point by
				       // old (reference) value and
				       // added shift
      a[i] = cell->vertex(vertex_mapping[i]) + shift_vector;
    }
}


// explicit instantiation
template class MappingQ1Eulerian<deal_II_dimension>;
