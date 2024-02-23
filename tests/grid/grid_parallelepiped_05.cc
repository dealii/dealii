// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test subdivided_parallelepiped

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
Point<dim>
point(double x = 0, double y = 0, double z = 0)
{
  Point<dim> p;
  if (dim > 0)
    p[0] = x;
  if (dim > 1)
    p[1] = y;
  if (dim > 2)
    p[2] = z;
  return p;
}


// The simplest test case is to create a parallelepiped grid with a
// number of subdivisions and output the result.
template <int dim, int spacedim>
void
check(bool subdivide)
{
  deallog << "dim " << dim << " spacedim " << spacedim << std::endl;

  Point<spacedim> origin;
  for (unsigned int d = 0; d < spacedim; ++d)
    origin[d] = 0.1 + d * 1.0;

  std::array<Tensor<1, spacedim>, dim> edges;
  switch (dim)
    {
      case 1:
        edges[0] = point<spacedim>(0.5, 0.02, 0.03);
        break;

      case 2:
        edges[0] = point<spacedim>(0.70, 0.25, -0.01);
        edges[1] = point<spacedim>(0.15, 0.50, 0.03);
        break;

      case 3:
        edges[0] = point<spacedim>(0.10, 0.50, 0.70);
        edges[1] = point<spacedim>(1.50, 0.25, 0.70);
        edges[2] = point<spacedim>(0.25, 0.50, 0.3);
        break;

      default:
        Assert(false, ExcInternalError());
    }


  std::vector<unsigned int> subdivisions;
  if (subdivide)
    {
      subdivisions.resize(dim);
      for (unsigned int d = 0; d < dim; ++d)
        subdivisions[d] = d + 1;
    }


  bool colorize = false;

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::subdivided_parallelepiped<dim, spacedim>(
    triangulation, origin, edges, subdivisions, colorize);

  GridOut grid_out;

  grid_out.write_gnuplot(triangulation, deallog.get_file_stream());
}

int
main()
{
  initlog();

  check<1, 1>(true);
  check<1, 2>(true);
  check<1, 3>(true);

  check<2, 2>(true);
  check<2, 3>(true);

  check<3, 3>(true);
  check<3, 3>(false);


  deallog << "OK" << std::endl;
}
