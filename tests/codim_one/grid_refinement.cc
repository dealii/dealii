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



// see what happens when creating a surface mesh and then refining it

#include "../tests.h"

// all include files you need here

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <string>

template <int dim, int spacedim>
void
test(std::string filename)
{
  SphericalManifold<dim, spacedim> boundary;
  Triangulation<dim, spacedim>     tria;
  GridIn<dim, spacedim>            gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);
  gi.read_ucd(in);
  tria.set_all_manifold_ids(1);
  tria.set_manifold(1, boundary);

  GridOut grid_out;
  grid_out.set_flags(GridOutFlags::Ucd(true));
  for (unsigned int cycle = 0; cycle < 3; ++cycle)
    {
      tria.refine_global(1);
      grid_out.write_msh(tria, deallog.get_file_stream());
    }
}

int
main()
{
  initlog();

  deallog << "Test<1,2>" << std::endl;
  test<1, 2>(SOURCE_DIR "/grids/circle_1.inp");

  deallog << std::endl << "Test<2,3>" << std::endl;
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}
