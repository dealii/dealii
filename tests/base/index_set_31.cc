// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet::split_by_block ()

#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  FESystem<dim>   fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  const std::vector<IndexSet> partition =
    dof_handler.locally_owned_dofs().split_by_block(dofs_per_block);
  for (unsigned int i = 0; i < 2; ++i)
    {
      partition[i].print(deallog);
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test<2>();
  test<3>();
}
