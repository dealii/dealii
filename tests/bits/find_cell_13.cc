// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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

// take a 2d mesh, provide a marked vertex such that the point is
// outside of active cells adjacent to the marked vertex
// and check that the function throws an exception
// that it couldn't find an arbitrary point's cell

// while not providing such marked vertex the function
// executes normally as expected

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

void check(Triangulation<2>& tria)
{
  // generate some random points bounded by [0., 0.2)^2 in R^2 space
  // any point in this domain should be inside one of the cells
  Point<2> p(random_value<double>(), random_value<double>() / 8.);

  try
    {
      // don't mark any vertex at first
      std::vector<bool> marked_vertices(tria.n_vertices(), false);

      // find the closes vertex to (1.,1.)
      // (default call to find_closest_vertex() without passing marked_vertices)
      unsigned int closest_vertex
        = GridTools::find_closest_vertex(tria, Point<2>(1., 1.));

      // mark the vertex closest to (1.,1.)
      marked_vertices[closest_vertex] = true;

      auto vertices = tria.get_vertices();
      deallog << vertices[closest_vertex] << " is the closest vertex"
              << std::endl;

      GridTools::find_active_cell_around_point(tria, p, marked_vertices);

      // The test fails if the control reaches here.
      deallog << "Garbage text to make the diff fail if the control "
              << "reaches here:"
              << " Point: " << p;
      deallog << std::endl; // Flush deallog buffer
    }
  catch(GridTools::ExcPointNotFound<2>)
    {
      deallog
        << "The first call to the function find_closest_vertex() has thrown "
        << "an exception. It is supposed to throw.";
      deallog << std::endl;

      // The default function call doesn't throw exceptions
      // for this use case
      Triangulation<2>::active_cell_iterator cell
        = GridTools::find_active_cell_around_point(tria, p);

      // check if the below function call actually finds appropriate cell
      Assert(p.distance(cell->center()) < cell->diameter() / 2,
             ExcInternalError());
      deallog << "Test passed!";
      deallog << std::endl;
    }
  catch(...)
    {
      Assert(false, ExcInternalError());
    }
}

int
main()
{
  initlog();

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube(coarse_grid);
    coarse_grid.refine_global(2);
    check(coarse_grid);
  }

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    static const SphericalManifold<2> boundary;
    coarse_grid.set_manifold(0, boundary);
    coarse_grid.refine_global(2);
    check(coarse_grid);
  }
}
