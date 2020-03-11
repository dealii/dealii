// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


// there are no direction flags if dim==spacedim. make sure we always
// get back true in that case

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace std;


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
