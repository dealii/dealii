// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria[2];

  GridGenerator::hyper_cube(tria[0]);
  GridGenerator::hyper_cube(tria[1]);

  tria[0].refine_global(2);
  tria[1].refine_global(2);

  tria[0].begin_active()->set_refine_flag();
  tria[0].execute_coarsening_and_refinement();

  tria[1].last_active()->set_refine_flag();
  tria[1].execute_coarsening_and_refinement();

  tria[1].last_active()->set_refine_flag();
  tria[1].execute_coarsening_and_refinement();

  DoFHandler<dim> dh0(tria[0]);
  DoFHandler<dim> dh1(tria[1]);

  using CellList =
    std::list<std::pair<typename DoFHandler<dim>::cell_iterator,
                        typename DoFHandler<dim>::cell_iterator>>;

  const CellList cell_list = GridTools::get_finest_common_cells(dh0, dh1);
  for (typename CellList::const_iterator cell_pair = cell_list.begin();
       cell_pair != cell_list.end();
       ++cell_pair)
    deallog << cell_pair->first << ' ' << cell_pair->second << std::endl;
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
