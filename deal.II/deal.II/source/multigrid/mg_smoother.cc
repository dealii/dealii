//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <grid/tria.h>
#include <dofs/dof_constraints.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_smoother.templates.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>

#include <algorithm>


//////////////////////////////////////////////////////////////////////



#if deal_II_dimension == 1

template <class VECTOR>
MGSmootherContinuous<VECTOR>::MGSmootherContinuous (
  const MGDoFHandler<1> &/*mg_dof*/,
  const unsigned int steps)
		:
		steps(steps)
{
  Assert (false, ExcNotImplemented());
}

#endif


#if deal_II_dimension > 1

template <class VECTOR>
template <int dim>
MGSmootherContinuous<VECTOR>::MGSmootherContinuous (
  const MGDoFHandler<dim> &mg_dof,
  const unsigned int steps)
		:
		steps(steps)
{
  const unsigned int n_levels = mg_dof.get_tria().n_levels();
  
				   // allocate the right number of
				   // elements
  interior_boundary_dofs.resize (n_levels-1);

				   // use a temporary to store the
				   // indices. this allows to not allocate
				   // to much memory by later copying the
				   // content of this vector to its final
				   // destination
  std::vector<unsigned int> boundary_dofs;

				   // temporary to hold the dof indices
				   // on a face between to levels
  std::vector<unsigned int> dofs_on_face (mg_dof.get_fe().dofs_per_face);

  for (unsigned int level=1; level<n_levels; ++level)
    {
      boundary_dofs.clear ();

				       // for each cell on this level:
				       // find out whether a face is
				       // at the boundary of this level's
				       // cells and if so add the dofs
				       // to the interior boundary dofs
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if ((cell->neighbor(face).state() == IteratorState::valid) &&
	      (static_cast<unsigned int>(cell->neighbor(face)->level())
	       == level-1))
	    {
					       // get indices of this face
	      cell->face(face)->get_mg_dof_indices (dofs_on_face);
					       // append them to the levelwise
					       // list
	      boundary_dofs.insert (boundary_dofs.end(),
				    dofs_on_face.begin(),
				    dofs_on_face.end());
	    };

				       // now sort the list of interior boundary
				       // dofs and eliminate duplicates
      std::sort (boundary_dofs.begin(), boundary_dofs.end());
      boundary_dofs.erase (std::unique (boundary_dofs.begin(),
					boundary_dofs.end()),
			   boundary_dofs.end());

				       // now finally copy the result
				       // for this level its destination
      interior_boundary_dofs[level-1] = boundary_dofs;
    };      
}

#endif


// explicit instantiations

// don't do the following instantiation in 1d, since there is a specialized
// function there
#if deal_II_dimension > 1
template MGSmootherContinuous<Vector<float> >::MGSmootherContinuous (
  const MGDoFHandler<deal_II_dimension>&,
  const unsigned int);
template MGSmootherContinuous<Vector<double> >::MGSmootherContinuous (
  const MGDoFHandler<deal_II_dimension>&,
  const unsigned int);
template MGSmootherContinuous<BlockVector<float> >::MGSmootherContinuous (
  const MGDoFHandler<deal_II_dimension>&,
  const unsigned int);
template MGSmootherContinuous<BlockVector<double> >::MGSmootherContinuous (
  const MGDoFHandler<deal_II_dimension>&,
  const unsigned int);
#endif

