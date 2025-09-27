// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



Triangulation<2, 2>
create_grid()
{
  Triangulation<2, 2> triangulation;

  double inner_radius = .5;
  double outer_radius = 1.;

  Point<2> center(0., 0.);

  GridGenerator::half_hyper_shell(
    triangulation, Point<2>(), inner_radius, outer_radius, 0, true);

  Triangulation<2>::active_cell_iterator

    cell = triangulation.begin_active(),
    endc = triangulation.end();

  return triangulation;
}

int
main()
{
  initlog();

  GridOut           grid_out;
  GridOutFlags::Svg svg_flags;

  svg_flags.coloring           = GridOutFlags::Svg::level_number;
  svg_flags.label_boundary_id  = true;
  svg_flags.background         = GridOutFlags::Svg::transparent;
  svg_flags.label_level_number = true;
  svg_flags.label_cell_index   = true;
  svg_flags.draw_legend        = true;
  svg_flags.draw_colorbar      = true;

  grid_out.set_flags(svg_flags);
  grid_out.write_svg(create_grid(), deallog.get_file_stream());

  return 0;
}
