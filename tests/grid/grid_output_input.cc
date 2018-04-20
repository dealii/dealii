// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Test output and input of meshes, to see if the meshes are different
// after a couple of global refinements

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>



template <int dim>
void test(std::ostream &out)
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

  GridTools::map_boundary_to_manifold_ids({1}, {0},tr);
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

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
