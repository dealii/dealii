/* $Id$ */

#include <grid/tria.h>
#include <grid/mg_dof.h>
#include <grid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <numerics/mg_smoother.h>
#include <lac/vector.h>

#include <algorithm>


#if deal_II_dimension == 1

MGSmoother::MGSmoother (const MGDoFHandler<1> &/*mg_dof*/) 
{
  Assert (false, ExcNotImplemented());
};

#endif


#if deal_II_dimension > 1

template <int dim>
MGSmoother::MGSmoother (const MGDoFHandler<dim> &mg_dof) 
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
  vector<int> boundary_dofs;

				   // temporary to hold the dof indices
				   // on a face between to levels
  vector<int> dofs_on_face (mg_dof.get_fe().dofs_per_face);

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



MGSmoother::~MGSmoother () 
{};



void
MGSmoother::set_zero_interior_boundary (const unsigned int  level,
					Vector<float>      &u) const
{
  if (level==0)
    return;
  else
    for (vector<int>::const_iterator p=interior_boundary_dofs[level-1].begin();
	 p!=interior_boundary_dofs[level-1].end(); ++p)
      u(*p) = 0;
};




void
MGSmoother::post_smooth (const unsigned int  level,
			 Vector<float>      &u) const
{
  pre_smooth (level, u);
};




// explicit instantiations
// don't do the following instantiation in 1d, since there is a specialized
// function there
#if deal_II_dimension > 1
template MGSmoother::MGSmoother (const MGDoFHandler<deal_II_dimension>&);
#endif
