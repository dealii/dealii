/* $Id$ */

#include <grid/mg_dof.h>
#include <grid/dof_levels.h>
#include <grid/dof_accessor.h>
#include <grid/dof_constraints.h>
#include <grid/tria_levels.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>



/* ------------------------ MGVertexDoFs ----------------------------------- */

template <int dim>
MGDoFHandler<dim>::MGVertexDoFs::MGVertexDoFs (const unsigned int coarsest_level,
					       const unsigned int n_levels,
					       const unsigned int dofs_per_vertex) :
		coarsest_level (coarsest_level),
		indices (new int[n_levels * dofs_per_vertex])
{
  Assert (indices != 0, ExcNoMemory ());

  for (unsigned int i=0; i<n_levels*dofs_per_vertex; ++i)
    indices[i] = -1;
};



template <int dim>
MGDoFHandler<dim>::MGVertexDoFs::~MGVertexDoFs () {
  delete[] indices;
};







/* ------------------------ MGDoFHandler ------------------------------------- */

template <int dim>
MGDoFHandler<dim>::MGDoFHandler (Triangulation<dim> *tria) :
		DoFHandler<dim> (tria)
{};



template <int dim>
MGDoFHandler<dim>::~MGDoFHandler () {};




template <int dim>
void MGDoFHandler<dim>::distribute_dofs (const FiniteElementBase<dim> &fe) {
				   // first distribute global dofs
  DoFHandler<dim>::distribute_dofs (fe);


				   // reserve space for the MG dof numbers
  reserve_space ();
  
};




#if deal_II_dimension == 1

template <>
void MGDoFHandler<1>::reserve_space () {
};

#endif


#if deal_II_dimension == 2

template <>
void MGDoFHandler<2>::reserve_space () {
};

#endif


// explicite instantiations
template class MGDoFHandler<deal_II_dimension>;
