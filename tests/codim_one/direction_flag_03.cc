// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// there are no direction flags if dim==spacedim. make sure we always
// get back true in that case

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  volume_mesh.refine_global(1);

  for (typename Triangulation<dim>::active_cell_iterator cell =
         volume_mesh.begin_active();
       cell != volume_mesh.end();
       ++cell)
    AssertThrow(cell->direction_flag() == true, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
