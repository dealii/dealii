// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// test GridTools::rotate(), output is checked for correctness visually

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


template <int dim, int spacedim>
void test ();

template <>
void test<2,2> ()
{
  const int dim = 2;
  const int spacedim = 2;
  Triangulation<dim, spacedim> tria;
  Point<spacedim> origin;
  std_cxx11::array<Tensor<1,spacedim>,dim> edges;
  edges[0] = Point<spacedim>(2.0, 0.0)-Point<spacedim>();
  edges[1] = Point<spacedim>(0.2, 0.8)-Point<spacedim>();
  GridGenerator::subdivided_parallelepiped<dim, spacedim> (tria,
                                                           origin,
                                                           edges);

  //GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(numbers::PI/3.0, tria);
  GridOut().write_gnuplot (tria, deallog.get_file_stream());
}

// version for <1,3>, <2,3>, and <3,3>
template <int dim, int spacedim>
void test ()
{
  Triangulation<dim, spacedim> tria;
  Point<spacedim> origin(0.1, 0.2, 0.3);
  std_cxx11::array<Tensor<1,spacedim>,dim> edges;
  edges[0] = Point<spacedim>(2.0, 0.0, 0.1)-Point<spacedim>();
  if (dim>=2)
    edges[1] = Point<spacedim>(0.2, 0.8, 0.15)-Point<spacedim>();
  if (dim>=3)
    edges[2] = Point<spacedim>(0.0, 0.0, 0.1)-Point<spacedim>();

  GridGenerator::subdivided_parallelepiped<dim, spacedim> (tria,
                                                           origin,
                                                           edges);

  //GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(numbers::PI/3.0, 0, tria);
  //GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(numbers::PI/10.0, 1, tria);
  //GridOut().write_gnuplot (tria, deallog.get_file_stream());
  GridTools::rotate(-numbers::PI/5.0, 2, tria);
  GridOut().write_gnuplot (tria, deallog.get_file_stream());
}


int main ()
{
  initlog ();

  test<2,2> ();
  test<1,3> ();
  test<2,3> ();
  test<3,3> ();

  return 0;
}
