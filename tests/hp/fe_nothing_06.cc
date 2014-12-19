// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test that FE_Nothing works as intended


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/lac/constraint_matrix.h>


#include <fstream>

// Create a mesh with hanging
// nodes and FEQ/FENothing interfaces
// in several admissible configurations
// we'd like to check.  In 2D the mesh
// looks like the following.
//
// +---------+----+----+
// |         | 1  | 1  |
// |    0    +----+----+
// |         | 1  | 1  |
// +----+----+----+----+
// | 0  | 0  |         |
// +----+----+    1    |
// | 0  | 0  |         |
// +----+----+---------+
//
// We then attempt to make hanging node
// constraints on this mesh.

template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  GridGenerator :: hyper_cube (triangulation, -0.5, 0.5);
  triangulation.refine_global(1);

  {
    typename Triangulation<dim> :: active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();

    for (; cell != endc; cell++)
      {
        Point<dim> center = cell->center();

        if (center[0] < 0)
          {
            cell->set_subdomain_id(1);
          }

        double h=0;
        for (unsigned d=0; d<dim; ++d) h += center[d];

        if (std::fabs(h) + 1e-6 > 0.25*dim)
          cell->set_refine_flag();
      }

    triangulation.execute_coarsening_and_refinement();
  }

  // create fe_collection and
  // distribute dofs. in this test,
  // we are looking at FESystems.

  hp::FECollection<dim>    fe_collection;

  fe_collection.push_back (FESystem<dim>(FE_Q<dim>(2), dim,
                                         FE_Q<dim>(1), 1));

  fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(), dim,
                                         FE_Nothing<dim>(), 1));

  hp::DoFHandler<dim>      dof_handler (triangulation);

  // loop over cells, and set cells
  // within a circle to be of type
  // FE_Nothing, while outside the
  // circle to be of type FE_Q(1)
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell != endc; cell++)
      {
        if (cell->subdomain_id()==1 )
          cell->set_active_fe_index(1);
        else
          cell->set_active_fe_index(0);
      }

    dof_handler.distribute_dofs (fe_collection);
  }


  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells()
          << std::endl
          << "   Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;


  // .... test constraint handling

  ConstraintMatrix constraints;

  DoFTools::make_hanging_node_constraints (dof_handler, constraints);

  constraints.close();

  deallog << "   Number of constraints:        "
          << constraints.n_constraints()
          << std::endl;

}



int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "Try dim == 1" << std::flush << std::endl;
  test<1> ();

  deallog << "Try dim == 2" << std::flush << std::endl;
  test<2> ();

  deallog << "Try dim == 3" << std::flush << std::endl;
  test<3> ();

  deallog << "OK" << std::endl;
}
