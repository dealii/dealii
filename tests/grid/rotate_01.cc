// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test GridTools::rotate(), output is checked for correctness visually

#include <deal.II/base/numbers.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim, int spacedim>
void
test();

template <>
void
test<1, 2>()
{
  const int                    dim      = 1;
  const int                    spacedim = 2;
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_sphere<spacedim>(tria);

  // GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(numbers::PI / 3.0, tria);
  GridOut().write_gnuplot(tria, deallog.get_file_stream());
}

template <>
void
test<2, 2>()
{
  const int                            dim      = 2;
  const int                            spacedim = 2;
  Triangulation<dim, spacedim>         tria;
  Point<spacedim>                      origin;
  std::array<Tensor<1, spacedim>, dim> edges;
  edges[0] = Point<spacedim>(2.0, 0.0) - Point<spacedim>();
  edges[1] = Point<spacedim>(0.2, 0.8) - Point<spacedim>();
  GridGenerator::subdivided_parallelepiped<dim, spacedim>(tria, origin, edges);

  // GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(numbers::PI / 3.0, tria);
  GridOut().write_gnuplot(tria, deallog.get_file_stream());
}

// version for <1,3>, <2,3>, and <3,3>
template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim>         tria;
  Point<spacedim>                      origin(0.1, 0.2, 0.3);
  std::array<Tensor<1, spacedim>, dim> edges;
  edges[0] = Point<spacedim>(2.0, 0.0, 0.1) - Point<spacedim>();
  if (dim >= 2)
    edges[1] = Point<spacedim>(0.2, 0.8, 0.15) - Point<spacedim>();
  if (dim >= 3)
    edges[2] = Point<spacedim>(0.0, 0.0, 0.1) - Point<spacedim>();

  GridGenerator::subdivided_parallelepiped<dim, spacedim>(tria, origin, edges);

  // GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(Tensor<1, 3>({1., 0., 0.}), numbers::PI / 3.0, tria);
  // GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(Tensor<1, 3>({0., 1., 0.}), numbers::PI / 10.0, tria);
  // GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(Tensor<1, 3>({0., 0., 1.}), -numbers::PI / 5.0, tria);
  GridOut().write_gnuplot(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();

  test<1, 2>();
  test<2, 2>();
  test<1, 3>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
