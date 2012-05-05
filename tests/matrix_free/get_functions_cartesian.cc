//------------------  get_functions_cartesian.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  get_functions_cartesian.cc  ------------------------


// this function tests the correctness of the implementation of matrix free
// operations in getting the function values, the function gradients, and the
// function Laplacians on a cartesian mesh (hyper cube). This tests whether
// cartesian meshes are treated correctly. The test case is without any
// constraints

#include "../tests.h"


std::ofstream logfile("get_functions_cartesian/output");

#include "get_functions_common.h"


template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  ConstraintMatrix constraints;
  constraints.close();
  do_test<dim, fe_degree, double> (dof, constraints);
}

