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


// Test that clearing a triangulation is behaving correctly for the transfinite
// interpolation on a test that is similar to transfinite_manifold_01.cc

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

  {
    Triangulation<dim, spacedim> tria;
    GridGenerator::hyper_ball(tria);

    // set all entities to the transfinite manifold except for the boundary
    // where we put the polar manifold
    tria.set_all_manifold_ids(1);
    tria.set_all_manifold_ids_on_boundary(0);
    tria.set_manifold(0, polar_manifold);
    manifold.initialize(tria);
    tria.set_manifold(1, manifold);

    do_test(tria);

    // clear the triangulation, set up a similar problem again and check that it
    // is sane.
    tria.clear();
    GridGenerator::hyper_ball(tria);
    tria.set_all_manifold_ids(1);
    tria.set_all_manifold_ids_on_boundary(0);
    tria.set_manifold(0, polar_manifold);
    manifold.initialize(tria);
    tria.set_manifold(1, manifold);

    do_test(tria);
  }
  {
    Triangulation<dim, spacedim> tria;
    GridGenerator::hyper_ball(tria);

    // set all entities to the transfinite manifold except for the boundary
    // where we put the polar manifold
    tria.set_all_manifold_ids(1);
    tria.set_all_manifold_ids_on_boundary(0);
    tria.set_manifold(0, polar_manifold);
    manifold.initialize(tria);
    tria.set_manifold(1, manifold);

    do_test(tria);
  }
}


int
main()
{
  initlog();

  test_polar<2, 2>();

  return 0;
}
