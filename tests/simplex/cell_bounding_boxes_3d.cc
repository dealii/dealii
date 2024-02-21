// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test bounding box creation for simplex 3D tria.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
create_simplex_grid(Triangulation<dim> &tria_simplex)
{
  tria_simplex.clear();
  Triangulation<dim> tria_standard;
  GridGenerator::hyper_cube(tria_standard);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_standard, tria_simplex);
}


template <int dim>
void
compute_bounding_boxes(Triangulation<dim> &tria)
{
  for (const auto &cell : tria.active_cell_iterators())
    {
      const auto cell_bb = cell->bounding_box();
      deallog << "Center : " << cell_bb.center()
              << " , volume : " << cell_bb.volume() << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(5);
  Triangulation<3> tria;
  create_simplex_grid(tria);
  compute_bounding_boxes(tria);
}
