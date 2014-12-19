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
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/lac/constraint_matrix.h>


#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  GridGenerator :: hyper_cube (triangulation, -0.5, 0.5);
  triangulation.refine_global(4);

  hp::FECollection<dim>    fe_collection;

  fe_collection.push_back (FE_Q<dim>(1));
  fe_collection.push_back (FE_Nothing<dim>());

  hp::DoFHandler<dim>      dof_handler (triangulation);

  // loop over cells, and set cells
  // within a circle to be of type
  // FE_Nothing, while outside the
  // circle to be of type FE_Q(1)

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell != endc; cell++)
    {
      Point<dim> center = cell->center();
      if (std::sqrt(center.square()) < 0.25 )
        cell->set_active_fe_index(1);
      else
        cell->set_active_fe_index(0);
    }

  dof_handler.distribute_dofs (fe_collection);

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

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
