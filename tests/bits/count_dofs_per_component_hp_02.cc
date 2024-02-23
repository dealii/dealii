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


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"

// check
//   DoFTools::
//   count_dofs_per_component (...);
// for the hp-case
//
// in contrast to the _01 test also check that this works if the element
// collection has more than one element


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  // define DoFhandler and FEs
  FE_Q<dim>     u1(2);
  FE_Q<dim>     p1(1);
  FESystem<dim> fe_system1(u1, 2, p1, 1);
  FE_Q<dim>     u2(3);
  FE_Q<dim>     p2(2);
  FESystem<dim> fe_system2(u2, 2, p2, 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe_system1);
  fe_collection.push_back(fe_system2);

  DoFHandler<dim> hp_dof_handler(triangulation);
  hp_dof_handler.begin_active()->set_active_fe_index(1);

  // distribute dofs
  hp_dof_handler.distribute_dofs(fe_collection);

  // count dofs per component and show them on the screen
  const std::vector<types::global_dof_index> dofs_per_component_hp =
    DoFTools::count_dofs_per_fe_component(hp_dof_handler);

  for (unsigned int i = 0; i < 3; ++i)
    {
      deallog << "DoFs in the " << i
              << ". component: " << dofs_per_component_hp.at(i) << std::endl;
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
