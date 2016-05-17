// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2016 by the deal.II authors
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


// test subdivided_parallelepiped

#include "../tests.h"
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

template <int dim>
Point<dim> point(double x=0, double y=0, double z=0)
{
  Point<dim> p;
  if (dim>0) p[0] = x;
  if (dim>1) p[1] = y;
  if (dim>2) p[2] = z;
  return p;
}


// The simplest test case is to create a parallelepiped grid with a
// number of subdivisions and output the result.
template<int dim, int spacedim>
void check (bool subdivide)
{
  deallog << "dim " << dim << " spacedim " << spacedim << std::endl;

  Point<spacedim> origin;
  for (unsigned int d=0; d<spacedim; ++d)
    origin[d] = 0.1+d*1.0;

  std_cxx11::array<Tensor<1,spacedim>,dim> edges;
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
      Assert (false, ExcInternalError ());
    }


  std::vector<unsigned int> subdivisions;
  if (subdivide)
    {
      subdivisions.resize(dim);
      for (unsigned int d=0; d<dim; ++d)
        subdivisions[d]=d+1;
    }


  bool colorize = false;

  Triangulation<dim,spacedim> triangulation;
  GridGenerator::subdivided_parallelepiped<dim,spacedim> (triangulation,
                                                          origin,
                                                          edges,
                                                          subdivisions,
                                                          colorize);

  GridOut grid_out;

  grid_out.write_gnuplot (triangulation, deallog.get_file_stream());
}

int main ()
{
  initlog();

  check<1,1> (true);
  check<1,2> (true);
  check<1,3> (true);

  check<2,2> (true);
  check<2,3> (true);

  check<3,3> (true);
  check<3,3> (false);


  deallog<< "OK" << std::endl;
}
