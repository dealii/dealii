// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
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
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
Point<dim>
trans_func(const Point<dim> &p);

template <>
Point<2>
trans_func(const Point<2> &p)
{
  Point<2> r(p[0] + p[1] * p[1], p[1]);
  return r;
}
template <>
Point<3>
trans_func(const Point<3> &p)
{
  Point<3> r(p[0] + p[1] * p[1], p[1], p[2]);
  return r;
}



template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();


  deallog << "Unchanged grid:" << std::endl;
  GridOut().write_gnuplot(tria, deallog.get_file_stream());

  {
    std::ofstream f("grid1");
    GridOut().write_gnuplot(tria, f);
  }

  GridTools::transform(trans_func<dim>, tria);

  deallog << "transformed grid:" << std::endl;
  GridOut().write_gnuplot(tria, deallog.get_file_stream());
  {
    std::ofstream f("grid2");
    GridOut().write_gnuplot(tria, f);
  }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
