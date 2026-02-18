// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

/*
 * Test that copying a FilteredIterator works
 */

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation);

  FilteredIterator<typename DoFHandler<dim>::level_cell_iterator> begin(
    IteratorFilters::LocallyOwnedLevelCell(), dof_handler.begin());
  FilteredIterator<typename DoFHandler<dim>::level_cell_iterator> end(
    IteratorFilters::LocallyOwnedLevelCell(), dof_handler.end());
  end = begin;

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<2>();

  return 0;
}
