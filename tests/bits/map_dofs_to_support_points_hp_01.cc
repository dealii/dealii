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
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"

// check
//   DoFTools::
//   map_dofs_to_support_points(...);
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

  MappingQ<dim>              mapping(2);
  hp::MappingCollection<dim> mapping_collection(mapping);

  DoFHandler<dim> dof_handler(triangulation);

  // distribute dofs
  dof_handler.distribute_dofs(fe_system);

  // now map the dofs to the support points and show them on the screen
  std::vector<Point<dim>> map(dof_handler.n_dofs());
  std::vector<Point<dim>> hp_map(dof_handler.n_dofs());

  DoFTools::map_dofs_to_support_points(mapping, dof_handler, map);
  DoFTools::map_dofs_to_support_points(mapping_collection, dof_handler, hp_map);

  // output the elements
  for (unsigned int i = 0; i < hp_map.size(); ++i)
    {
      // both maps should contain the same
      Assert(hp_map[i] == map[i], ExcInternalError());
      deallog << hp_map[i] << ' ';
    }
  deallog << std::endl;
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
