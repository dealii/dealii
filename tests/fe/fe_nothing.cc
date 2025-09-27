// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <iostream>

#include "../tests.h"

const double eps = 1e-10;


template <int dim>
void
test2cells(const unsigned int p1 = 2, const unsigned int p2 = 1)
{
  deallog << "2cells: " << dim << "d, p1=" << p1 << ", p2=" << p2 << std::endl;
  Triangulation<dim> triangulation;
  {
    Triangulation<dim> triangulationL;
    Triangulation<dim> triangulationR;
    GridGenerator::hyper_cube(triangulationL,
                              -1,
                              0); // create a square [-1,0]^d domain
    GridGenerator::hyper_cube(triangulationR,
                              -1,
                              0); // create a square [-1,0]^d domain
    Point<dim> shift_vector;
    shift_vector[0] = 1.0;
    GridTools::shift(shift_vector, triangulationR);
    GridGenerator::merge_triangulations(triangulationL,
                                        triangulationR,
                                        triangulation);
  }

  DoFHandler<dim> dof_handler(triangulation);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(p1), 1, FE_Q<dim>(p2), 1));
  fe_collection.push_back(
    FESystem<dim>(FE_Q<dim>(p1), 1, FE_Nothing<dim>(1, true), 1));

  // default is 0
  dof_handler.begin_active()->set_active_fe_index(1);

  dof_handler.distribute_dofs(fe_collection);

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());

  dof_handler.clear();
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4) << std::fixed;
  deallog.depth_console(0);


  try
    {
      test2cells<2>();
      test2cells<3>();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
