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


// Similar to grid_transform_02.cc but use coefficient function and solve
// for displacement field.


#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"


template <int dim>
class L2_inverse : public Function<dim>
{
public:
  L2_inverse(const std::vector<Point<dim>> &distance_source)
    : Function<dim>()
    , distance_source(distance_source)
  {
    Assert(distance_source.size() > 0, ExcNotImplemented());
  }

  virtual double
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const
  {
    double l2_inverse = std::numeric_limits<double>::max();

    for (unsigned int d = 0; d < distance_source.size(); d++)
      l2_inverse = std::min((p - distance_source[d]).norm_square(), l2_inverse);

    l2_inverse = std::max(l2_inverse, 1.e-5);

    return 1.0 / l2_inverse;
  }

private:
  const std::vector<Point<dim>> distance_source;
};


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

  // store the current location to be used in coefficient function
  std::vector<Point<dim>> metric;
  metric.push_back(best_pos);

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
  L2_inverse<dim> coefficient(metric);
  GridTools::laplace_transform(new_points, tria, &coefficient);
  std::ofstream out("output");
  GridOut       grid_out;
  grid_out.write_eps(tria, out);
}
