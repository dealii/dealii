// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test output and input of meshes, to see if the meshes are different
// after a couple of global refinements

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  GridOut go;
  go.set_flags(GridOutFlags::Ucd(false, true, true));
  go.set_flags(GridOutFlags::Msh(true, true));
  Triangulation<dim> tr;

  GridGenerator::hyper_cube_with_cylindrical_hole(tr, .3, .4, 1, 1, false);
  tr.reset_manifold(0);
  GridTools::copy_boundary_to_manifold_id(tr);
  CylindricalManifold<dim> boundary(2);
  tr.set_manifold(1, boundary);
  {
    std::ofstream grid_file("coarse_grid.inp");
    go.write_ucd(tr, grid_file);
    grid_file.close();
  }
  tr.refine_global(1);
  deallog << "Writing refined from constructor" << std::endl;
  go.write_ucd(tr, out);
  tr.clear();

  GridIn<dim> gi;
  gi.attach_triangulation(tr);
  {
    deallog << "Read coarse grid" << std::endl;
    std::ifstream grid_file("coarse_grid.inp");
    gi.read_ucd(grid_file);
    grid_file.close();
  }

  GridTools::map_boundary_to_manifold_ids({1}, {0}, tr);
  tr.set_manifold(0, boundary);
  tr.refine_global(1);
  deallog << "Writing refined from file" << std::endl;
  go.write_ucd(tr, out);

  {
    std::ofstream grid_file("grid.msh");
    go.write_msh(tr, grid_file);
    grid_file.close();
  }
  tr.clear();
}

int
main()
{
  initlog();

  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
