//------------------  get_functions_float.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  get_functions_float.cc  ------------------------


// this function tests the correctness of the implementation of matrix free
// operations in getting the function values, the function gradients, and the
// function Laplacians on a mesh with all different cell types (same as in
// matrix_vector_06) with hanging nodes, boundary conditions for float values.

#include "../tests.h"
#include <deal.II/base/function.h>
#include "create_mesh.h"

std::ofstream logfile("get_functions_float/output");

#include "get_functions_common.h"


template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  create_mesh (tria);
  tria.refine_global(4-dim);

				// refine a few cells
  for (unsigned int i=0; i<10-3*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
	cell = tria.begin_active (),
	endc = tria.end();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
	if (counter % (7-i) == 0)
	  cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  VectorTools::interpolate_boundary_values (dof, 1, ZeroFunction<dim>(),
					    constraints);
  constraints.close();

  do_test <dim, fe_degree, float> (dof, constraints);
}

