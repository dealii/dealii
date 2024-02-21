// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// A test that checks that FE_Nothing does not dominate the
// neighboring Q1. This works fine for a single, scalar element as
// shown in the _01 test, but fails at the time the test was written
// for systems


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <iostream>
#include <vector>

#include "../tests.h"


const unsigned int dim = 2;

void
print_dofs(const DoFHandler<2>::active_cell_iterator &cell)
{
  deallog << "DoFs on cell=" << cell << ": ";

  std::vector<types::global_dof_index> dof_indices(
    cell->get_fe().dofs_per_cell);
  cell->get_dof_indices(dof_indices);
  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    deallog << dof_indices[i] << ' ';
  deallog << std::endl;
}


int
main()
{
  initlog();

  Triangulation<dim>        triangulation;
  std::vector<unsigned int> subdivisions(dim, 1U);
  subdivisions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            subdivisions,
                                            Point<dim>(0, 0),
                                            Point<dim>(2, 1));

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Nothing<dim>(), 1));
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));

  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.begin_active()->set_active_fe_index(1);

  dof_handler.distribute_dofs(fe_collection);

  print_dofs(dof_handler.begin_active());
  print_dofs(std::next(dof_handler.begin_active()));

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());
}
