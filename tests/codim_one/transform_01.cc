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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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

  std::set<types::boundary_id> boundary_ids;
  boundary_ids.insert(0);

  GridGenerator::extract_boundary_mesh(
    volume_mesh, triangulation, boundary_ids);
  triangulation.reset_manifold(0);
  GridTools::transform(&warp<3>, triangulation);

  deallog << "Surface mesh has " << triangulation.n_active_cells() << " cells."
          << std::endl;
  save_mesh(triangulation);

  return 0;
}
