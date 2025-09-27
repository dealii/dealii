// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that GridGenerator::merge_triangulations correctly copies
// manifold ids if requested.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  SphericalManifold<dim>                spherical_manifold;
  TransfiniteInterpolationManifold<dim> inner_manifold;
  Triangulation<dim>                    ball_triangulation;
  GridGenerator::hyper_ball(ball_triangulation);
  ball_triangulation.set_all_manifold_ids(1);
  ball_triangulation.set_all_manifold_ids_on_boundary(0);
  ball_triangulation.set_manifold(0, spherical_manifold);
  inner_manifold.initialize(ball_triangulation);
  ball_triangulation.set_manifold(1, inner_manifold);

  Triangulation<dim> cube_triangulation;
  GridGenerator::hyper_cube(cube_triangulation);

  Triangulation<dim> merged_triangulation;
  GridGenerator::merge_triangulations(
    ball_triangulation, cube_triangulation, merged_triangulation, 1.e-12, true);
  TransfiniteInterpolationManifold<dim> inner_manifold_copied;
  merged_triangulation.set_manifold(0, spherical_manifold);
  inner_manifold_copied.initialize(merged_triangulation);
  merged_triangulation.set_manifold(1, inner_manifold_copied);

  ball_triangulation.refine_global(1);
  cube_triangulation.refine_global(1);
  merged_triangulation.refine_global(1);

  std::stringstream separate;
  GridOut().write_gnuplot(ball_triangulation, separate);
  GridOut().write_gnuplot(cube_triangulation, separate);

  std::stringstream merged;
  GridOut().write_gnuplot(merged_triangulation, merged);
  deallog << merged.str() << std::endl;
  AssertThrow(separate.str() == merged.str(), ExcInternalError());
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();

  return 0;
}
