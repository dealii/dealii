// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2017 by the deal.II authors
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

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>


using namespace dealii;

Triangulation<2,2> create_grid()
{
  Triangulation<2,2> triangulation;

  double inner_radius = .5;
  double outer_radius = 1.;

  Point<2> center(0., 0.);

  GridGenerator::hyper_cube_with_cylindrical_hole(triangulation, inner_radius, outer_radius);
  triangulation.reset_manifold(0);
  triangulation.refine_global(1);

  Triangulation<2>::active_cell_iterator

  cell = triangulation.begin_active(),
  endc = triangulation.end();

  for (; cell!=endc; ++cell)
    {
      for (unsigned int v=0; v < GeometryInfo<2>::vertices_per_cell; ++v)
        {
          const double distance_from_center = center.distance(cell->vertex(v));

          if (std::fabs(distance_from_center - inner_radius) < .25)
            {
              cell->set_refine_flag();
              break;
            }
        }
    }

  triangulation.execute_coarsening_and_refinement();

  return triangulation;
}

int main()
{
  initlog();

  GridOut grid_out;
  GridOutFlags::Svg svg_flags;

  svg_flags.coloring = GridOutFlags::Svg::level_number;
  svg_flags.label_material_id = true;
  svg_flags.background = GridOutFlags::Svg::transparent;

  grid_out.set_flags(svg_flags);
  grid_out.write_svg(create_grid(), deallog.get_file_stream());

  return 0;
}
