// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check GridTools::transform with the mesh used in the "Possibilities for
// extensions" section of step-38. The test exists because the function
// originally did not allow application to meshes that were already refined.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
Point<dim>
warp(const Point<dim> &p)
{
  Point<dim> q = p;
  q[dim - 1] *= 10;

  if (dim >= 2)
    q[0] += 2 * std::sin(q[dim - 1]);
  if (dim >= 3)
    q[1] += 2 * std::cos(q[dim - 1]);

  return q;
}


template <int dim, int spacedim>
void
save_mesh(const Triangulation<dim, spacedim> &tria)
{
  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();

  Triangulation<2, 3> triangulation;

  SphericalManifold<3> boundary_description;
  Triangulation<3>     volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);

  volume_mesh.set_manifold(0, boundary_description);
  volume_mesh.refine_global(3);

  static SphericalManifold<3 - 1, 3> surface_description;
  triangulation.set_manifold(0, surface_description);

  const std::set<types::boundary_id> boundary_ids = {0};
  GridGenerator::extract_boundary_mesh(volume_mesh,
                                       triangulation,
                                       boundary_ids);
  triangulation.reset_manifold(0);
  GridTools::transform(&warp<3>, triangulation);

  deallog << "Surface mesh has " << triangulation.n_active_cells() << " cells."
          << std::endl;
  save_mesh(triangulation);

  return 0;
}
