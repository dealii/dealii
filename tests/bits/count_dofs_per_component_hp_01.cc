// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// check
//   DoFTools::
//   count_dofs_per_component (...);
// for the hp-case


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  // define DoFhandler and FEs
  FE_Q<dim> u(2);
  FE_Q<dim> p(1);

  FESystem<dim> fe_system(u, 2, p, 1);

  DoFHandler<dim> dof_handler(triangulation);

  // distribute dofs
  dof_handler.distribute_dofs(fe_system);

  // count dofs per component and show them on the screen
  const std::vector<types::global_dof_index> dofs_per_component =
    DoFTools::count_dofs_per_fe_component(dof_handler);

  for (unsigned int i = 0; i < 3; ++i)
    {
      deallog << "DoFs in the " << i
              << ". component: " << dofs_per_component.at(i) << std::endl;
    }
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
  return 0;
}
