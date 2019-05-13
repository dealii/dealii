// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
