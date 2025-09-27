// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test, if the grid is written correctly in dx format. the vertices have to be
// renumbered after coarsening, as the dx format uses an implicit vertex
// numbering.

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  GridOut grid_out;
  grid_out.write_dx(tria, deallog.get_file_stream());

  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active();
  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      cell->set_coarsen_flag();
      ++cell;
    }
  tria.execute_coarsening_and_refinement();
  grid_out.write_dx(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<2>();
  test<3>();
}
