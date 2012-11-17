//------------------------  compress_constraints.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------------  compress_constraints.cc  ------------------------

// this function tests whether the compression of constraint weights
// (constraint pool) works properly

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

#include "create_mesh.h"

std::ofstream logfile("compress_constraints/output");



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  create_mesh (tria);
  tria.begin_active ()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active ();
  endc = tria.end();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(1);
  for (unsigned int i=0; i<10-3*dim; ++i)
    {
      cell = tria.begin_active ();
      endc = tria.end();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
        if (counter % (7-i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim> fe (2);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values (dof, 0, ZeroFunction<dim>(),
                                            constraints);
  constraints.close();

  const QGauss<1> quad(2);
  MatrixFree<dim> mf;
  mf.reinit (dof, constraints, quad);

  deallog << "Number of hanging nodes: "
          << constraints.n_constraints() << std::endl;
  deallog << "Number of different constraint weights: "
          << mf.n_constraint_pool_entries() << std::endl;
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision (3);

  {
    deallog.threshold_double(5.e-11);
    deallog.push("2d");
    test<2>();
    deallog.pop();
    deallog.push("3d");
    test<3>();
    deallog.pop();
  }
}
