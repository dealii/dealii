// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test extract_boundary_mesh on a 3D simplex mesh

#include "deal.II/grid/grid_tools_geometry.h"
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
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> hex_tria, triangulation;
  GridGenerator::hyper_cube(hex_tria);
  hex_tria.refine_global(3);
  GridGenerator::convert_hypercube_to_simplex_mesh(hex_tria, triangulation);

  // now extract the surface mesh
  Triangulation<dim - 1, dim> triangulation_surface;

  GridGenerator::extract_boundary_mesh(triangulation, triangulation_surface);
  double measure_boundary_mesh = GridTools::volume(triangulation_surface);
  if constexpr (dim == 2)
    AssertThrow(
      std::fabs(measure_boundary_mesh - 4.0) < 1e-12,
      ExcMessage(
        "The measure of the boundary of the unit square should be 4, but is " +
        std::to_string(measure_boundary_mesh)));
  else if constexpr (dim == 3)
    AssertThrow(
      std::fabs(measure_boundary_mesh - 6.0) < 1e-12,
      ExcMessage(
        "The measure of the boundary of the unit cube should be 6, but is " +
        std::to_string(measure_boundary_mesh)));
  deallog << "Ok" << std::endl;
}


int
main()
{
  initlog();
  deallog.depth_console(0);

  test<2>(); // we know this already works...
  test<3>();
}
