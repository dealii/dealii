// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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
  Triangulation<1, 1> tr1;
  GridGenerator::hyper_cube(tr1);
  std::array<Point<1>, 2>     diags1 = {{Point<1>(), Point<1>(1.)}};
  std::array<unsigned int, 1> reps1  = {{2}};
  Triangulation<1>            res1;
  GridGenerator::replicate_triangulation(tr1, diags1, reps1, res1);
  deallog << "1d triangulation has the following cell centers:" << std::endl;
  for (typename Triangulation<1>::active_cell_iterator cell =
         res1.begin_active();
       cell != res1.end();
       ++cell)
    deallog << "(" << cell->center()[0] << ")" << std::endl;

  deallog << "Number of global vertices: " << res1.n_vertices() << std::endl;

  Triangulation<2, 2> tr2;
  GridGenerator::hyper_cube_with_cylindrical_hole(tr2, 0.3, 0.5);
  const std::array<Point<2>, 2> diags2 = {
    {Point<2>(-0.5, -0.5), Point<2>(0.5, 0.5)}};
  const std::array<unsigned int, 2> reps2 = {{2, 2}};
  Triangulation<2, 2>               res2;
  GridGenerator::replicate_triangulation(tr2, diags2, reps2, res2);
  deallog << "2d triangulation has the following cell centers:" << std::endl;
  for (typename Triangulation<2>::active_cell_iterator cell =
         res2.begin_active();
       cell != res2.end();
       ++cell)
    deallog << "(" << cell->center()[0] << "," << cell->center()[1] << ")"
            << std::endl;

  deallog << "Number of global vertices: " << res2.n_vertices() << std::endl;

  Triangulation<3, 3>     tr3;
  std::array<Point<3>, 2> diags3 = {
    {Point<3>(-0.7, -0.5, -1.5), Point<3>(0.5, 0.5, 2.)}};
  std::array<unsigned int, 3> reps3 = {{2, 2, 2}};
  std::vector<unsigned int>   rp(reps3.begin(), reps3.end());
  GridGenerator::subdivided_hyper_rectangle(tr3,
                                            rp,
                                            diags3.front(),
                                            diags3.back());
  Triangulation<3, 3> res3;
  GridGenerator::replicate_triangulation(tr3, diags3, reps3, res3);
  deallog << "3d triangulation has the following cell centers:" << std::endl;
  for (typename Triangulation<3>::active_cell_iterator cell =
         res3.begin_active();
       cell != res3.end();
       ++cell)
    deallog << "(" << cell->center()[0] << "," << cell->center()[1] << ","
            << cell->center()[2] << ")" << std::endl;

  deallog << "Number of global vertices: " << res3.n_vertices() << std::endl;
}

int
main()
{
  initlog();

  test();

  return 0;
}
