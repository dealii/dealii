// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like grid_in_out but test persistent triangulation

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/persistent_tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <string>

#include "../tests.h"

template <int dim, int spacedim>
void
test(std::string filename)
{
  deallog.get_file_stream()
    << " Tria<" << dim << ',' << spacedim << ">: " << std::endl;
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);
  gi.read_ucd(in);

  PersistentTriangulation<dim, spacedim>                      ptria(tria);
  typename Triangulation<dim, spacedim>::active_cell_iterator cell;

  ptria.restore();
  for (unsigned int i = 0; i < 2; ++i)
    {
      for (cell = ptria.begin_active(); cell != ptria.end(); ++cell)
        if (cell->center()[0] < .5)
          cell->set_refine_flag();
      ptria.execute_coarsening_and_refinement();
    }

  ptria.write_flags(deallog.get_file_stream());
}

int
main()
{
  initlog();

  test<2, 3>(SOURCE_DIR "/grids/square.inp");
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}
