// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

// test GridGenerator::replicate_triangulation

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

void
test()
{
  {
    Triangulation<1, 1> tr1;
    GridGenerator::hyper_cube(tr1);
    std::vector<unsigned int> reps1 = {2};
    Triangulation<1>          res1;
    GridGenerator::replicate_triangulation(tr1, reps1, res1);
    deallog << "1d triangulation has the following cell centers:" << std::endl;
    for (const auto &cell : res1.active_cell_iterators())
      deallog << '(' << cell->center() << ')' << std::endl;

    deallog << "Number of global vertices: " << res1.n_vertices() << std::endl;
  }

  {
    Triangulation<2, 2> tr2;
    GridGenerator::hyper_cube_with_cylindrical_hole(tr2, 0.3, 0.5);
    const std::vector<unsigned int> reps2 = {{2, 2}};
    Triangulation<2, 2>             res2;
    GridGenerator::replicate_triangulation(tr2, reps2, res2);
    deallog << "2d triangulation has the following cell centers:" << std::endl;
    for (const auto &cell : res2.active_cell_iterators())
      deallog << '(' << cell->center() << ')' << std::endl;

    deallog << "Number of global vertices: " << res2.n_vertices() << std::endl;
  }

  {
    Triangulation<3, 3>       tr3;
    std::vector<unsigned int> reps3 = {2, 2, 2};
    GridGenerator::subdivided_hyper_rectangle(tr3,
                                              reps3,
                                              Point<3>(-0.7, -0.5, -1.5),
                                              Point<3>(0.5, 0.5, 2.));
    Triangulation<3, 3> res3;
    GridGenerator::replicate_triangulation(tr3, reps3, res3);
    deallog << "3d triangulation has the following cell centers:" << std::endl;
    for (const auto &cell : res3.active_cell_iterators())
      deallog << '(' << cell->center() << ')' << std::endl;

    deallog << "Number of global vertices: " << res3.n_vertices() << std::endl;
  }
}

int
main()
{
  initlog();
  test();
}
