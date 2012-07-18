//------------------  matrix_vector_04.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  matrix_vector_04.cc  ------------------------


// this function tests the correctness of the implementation of matrix free
// matrix-vector products by comparing with the result of deal.II sparse
// matrix. The mesh uses a hypershell mesh without hanging nodes (only cell
// type: 2)

#include "../tests.h"

std::ofstream logfile("matrix_vector_04/output");

#include "matrix_vector_common.h"


template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_shell (tria, Point<dim>(),
                              0.5, 1., 96, true);
  static const HyperShellBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
  tria.set_boundary (1, boundary);
  if (dim == 2)
    tria.refine_global (2);

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  do_test<dim, fe_degree, double> (dof, constraints);
}
