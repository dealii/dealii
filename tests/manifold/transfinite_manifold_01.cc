// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

// Test that the transfinite manifold creates good points for a circle, a
// cylinder, and a sphere

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"


template <int dim, int spacedim>
void
do_test(const Triangulation<dim, spacedim> &tria)
{
  for (typename Triangulation<dim, spacedim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    {
      deallog << "Lines on cell with center: " << cell->center() << std::endl;
      for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell;
           ++line)
        deallog << cell->line(line)->center(/*respect_manifold=*/true)
                << std::endl;
      deallog << "Faces on cell with center: " << cell->center() << std::endl;
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        deallog << cell->face(face)->center(/*respect_manifold=*/true)
                << std::endl;
      deallog << "Center with manifold: " << cell->center(true) << std::endl;
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        if (cell->at_boundary(face))
          {
            std::vector<Point<spacedim>> points;
            points.push_back(cell->face(face)->vertex(0));
            points.push_back(cell->face(face)->vertex(1));
            std::vector<double> weights(2);
            weights[0] = 0.1;
            weights[1] = 0.9;
            Point<spacedim> p =
              cell->get_manifold().get_new_point(make_array_view(points),
                                                 make_array_view(weights));
            Point<spacedim> pref =
              cell->face(face)->get_manifold().get_new_point(
                make_array_view(points), make_array_view(weights));
            deallog << "Distance between cell manifold and face manifold: "
                    << (pref - p) << std::endl;
            weights[0] = 0.55;
            weights[1] = 0.45;
            p    = cell->get_manifold().get_new_point(make_array_view(points),
                                                   make_array_view(weights));
            pref = cell->face(face)->get_manifold().get_new_point(
              make_array_view(points), make_array_view(weights));
            deallog << "Distance between cell manifold and face manifold: "
                    << (pref - p) << std::endl;
          }
    }
  deallog << std::endl;
}

template <int dim, int spacedim>
void
test_polar()
{
  deallog << "Testing with PolarManifold dim=" << dim
          << ", spacedim=" << spacedim << std::endl;

  PolarManifold<dim, spacedim>                    polar_manifold;
  TransfiniteInterpolationManifold<dim, spacedim> manifold;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_ball(tria);

  // set all entities to the transfinite manifold except for the boundary where
  // we put the polar manifold
  tria.set_all_manifold_ids(1);
  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, polar_manifold);
  manifold.initialize(tria);
  tria.set_manifold(1, manifold);

  do_test(tria);
}

template <int dim, int spacedim>
void
test_spherical()
{
  deallog << "Testing with SphericalManifold dim=" << dim
          << ", spacedim=" << spacedim << std::endl;

  SphericalManifold<dim, spacedim>                polar_manifold;
  TransfiniteInterpolationManifold<dim, spacedim> manifold;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_ball(tria);

  // set all entities to the transfinite manifold except for the boundary where
  // we put the polar manifold
  tria.set_all_manifold_ids(1);
  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, polar_manifold);
  manifold.initialize(tria);
  tria.set_manifold(1, manifold);

  do_test(tria);
}

void
test_cylinder(unsigned int ref = 1)
{
  const unsigned int dim = 3, spacedim = 3;
  deallog << "Testing with CylindricalManifold in 3d" << std::endl;

  CylindricalManifold<dim, spacedim>              cylinder_manifold;
  TransfiniteInterpolationManifold<dim, spacedim> manifold;
  Triangulation<dim, spacedim>                    tria;
  GridGenerator::cylinder(tria);
  tria.set_all_manifold_ids(1);

  for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        if (cell->at_boundary(face))
          {
            bool cell_at_surfaces = true;
            for (unsigned int i = 1; i < GeometryInfo<dim>::vertices_per_face;
                 ++i)
              if (std::abs(cell->face(face)->vertex(i)[0] -
                           cell->face(face)->vertex(0)[0]) > 1e-10)
                cell_at_surfaces = false;
            if (cell_at_surfaces == false)
              cell->face(face)->set_all_manifold_ids(0);
          }
    }
  tria.set_manifold(0, cylinder_manifold);
  manifold.initialize(tria);
  tria.set_manifold(1, manifold);

  do_test(tria);
}

int
main()
{
  initlog();

  test_polar<2, 2>();
  test_spherical<2, 2>();
  test_spherical<3, 3>();
  test_cylinder();

  return 0;
}
