// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that GridTools::laplace_transform can deal with interior
// nodes being pinned to a new location as well. The test itself
// doesn't make much sense since it leads to a few inverted cells, but
// it allows for easy visual inspection that the desired result
// happens.
//
// (Testcase adapted from one by Denis Davydov.)


#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



int
main()
{
  const int dim = 2;

  Triangulation<dim>                 tria;
  std::map<unsigned int, Point<dim>> new_points;
  const unsigned int                 N = 8;
  GridGenerator::subdivided_hyper_cube(tria, N, -5, 5);

  // find the vertex at the origin
  Triangulation<dim>::active_cell_iterator cell =
    GridTools::find_active_cell_around_point(tria, Point<dim>());

  unsigned int best_vertex =
    cell->vertex_index(0); // vertex number on local triangulation
  Point<dim> best_pos  = cell->vertex(0);
  double     best_dist = Point<dim>().distance(best_pos);

  for (unsigned int vertex_no = 1;
       vertex_no < GeometryInfo<dim>::vertices_per_cell;
       vertex_no++)
    {
      const double dist = Point<dim>().distance(cell->vertex(vertex_no));
      if (dist < best_dist)
        {
          best_pos    = cell->vertex(vertex_no);
          best_vertex = cell->vertex_index(vertex_no);
          best_dist   = dist;
        }
    }
  // move the point at the origin by 1 unit to the right
  new_points[best_vertex] = Point<dim>();
  new_points[best_vertex][0] += 1.;

  // now pin all of the points on the boundary
  cell                                          = tria.begin_active();
  Triangulation<dim>::active_cell_iterator endc = tria.end();

  for (; cell != endc; ++cell)
    if (cell->at_boundary() == true)
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        if (cell->face(face)->at_boundary() == true)
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
               ++v)
            {
              unsigned int vertex_number = cell->face(face)->vertex_index(v);
              new_points[vertex_number]  = cell->face(face)->vertex(v);
            }

  // then compute new point locations and output the result
  GridTools::laplace_transform<dim>(new_points, tria, nullptr, true);
  std::ofstream out("output");
  GridOut       grid_out;
  grid_out.write_eps(tria, out);
}
