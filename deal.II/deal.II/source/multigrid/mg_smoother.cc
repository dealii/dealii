//----------------------------  mg_smoother.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_smoother.cc  ---------------------------


#include <grid/tria.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <multigrid/multigrid.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_smoother.templates.h>
#include <lac/vector.h>

#include <algorithm>


//////////////////////////////////////////////////////////////////////


MGSmootherBase::~MGSmootherBase()
{};


//////////////////////////////////////////////////////////////////////


void
MGSmootherIdentity::smooth (const unsigned int,
			    Vector<double>       &,
			    const Vector<double> &) const
{};


#if deal_II_dimension == 1

MGSmoother::MGSmoother (const MGDoFHandler<1> &/*mg_dof*/, unsigned int steps)
		:
		steps(steps)
{
  Assert (false, ExcNotImplemented());
};

#endif


#if deal_II_dimension > 1

template <int dim>
MGSmoother::MGSmoother (const MGDoFHandler<dim> &mg_dof, unsigned int steps)
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
  vector<unsigned int> boundary_dofs;

				   // temporary to hold the dof indices
				   // on a face between to levels
  vector<unsigned int> dofs_on_face (mg_dof.get_fe().dofs_per_face);

  for (unsigned int level=1; level<n_levels; ++level)
    {
      boundary_dofs.clear ();

				       // for each cell on this level:
				       // find out whether a face is
				       // at the boundary of this level's
				       // cells and if so add the dofs
				       // to the interior boundary dofs
      for (MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if ((cell->neighbor(face).state() == valid) &&
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
      sort (boundary_dofs.begin(), boundary_dofs.end());
      boundary_dofs.erase (unique (boundary_dofs.begin(),
				   boundary_dofs.end()),
			   boundary_dofs.end());

				       // now finally copy the result
				       // for this level its destination
      interior_boundary_dofs[level-1] = boundary_dofs;
    };      
};

#endif

//////////////////////////////////////////////////////////////////////

void
MGSmoother::set_zero_interior_boundary (const unsigned int  level,
					Vector<double>      &u) const
{
  if (level==0)
    return;
  else
    for (vector<unsigned int>::const_iterator p=interior_boundary_dofs[level-1].begin();
	 p!=interior_boundary_dofs[level-1].end(); ++p)
      u(*p) = 0;
};

//////////////////////////////////////////////////////////////////////


// explicit instantiations
// don't do the following instantiation in 1d, since there is a specialized
// function there
#if deal_II_dimension > 1
template MGSmoother::MGSmoother (const MGDoFHandler<deal_II_dimension>&, unsigned int);
#endif

template
MGSmootherRelaxation<float>
::MGSmootherRelaxation(const MGDoFHandler<deal_II_dimension>        &mg_dof,
		       const MGLevelObject<SparseMatrix<float> > &matrix,
		       const function_ptr                            relaxation,
		       const unsigned int                            steps,
		       const double                                  omega);

template
MGSmootherRelaxation<double>
::MGSmootherRelaxation(const MGDoFHandler<deal_II_dimension>         &mg_dof,
		       const MGLevelObject<SparseMatrix<double> > &matrix,
		       const function_ptr                             relaxation,
		       const unsigned int                             steps,
		       const double                                   omega);

