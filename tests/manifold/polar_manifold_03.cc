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


// Test the push_forward and pull_back mechanisms

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>


// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim " << dim
          << ", spacedim " << spacedim << std::endl;

  PolarManifold<dim,spacedim> manifold;

  Triangulation<dim,spacedim> tria;
  Point<spacedim> p0;
  Point<spacedim> p1;
  p0[0] = .2;
  p1[0] = 1;
  p0[1] = .1;

  if (spacedim == 2)
    {
      p1[1] = 2*numbers::PI-.1; // theta
    }
  else if (spacedim == 3)
    {
      p1[1] = numbers::PI-.1;
      p1[2] = 2*numbers::PI-.1;
    }

  GridGenerator::hyper_rectangle (tria, p0, p1);
  tria.refine_global(3);

  const std::vector<Point<spacedim> > &vertices = tria.get_vertices();

  for (unsigned int i=0; i<vertices.size(); ++i)
    {
      Point<spacedim> p0 = manifold.push_forward(vertices[i]);
      Point<spacedim> p1 = manifold.pull_back(p0);

      if (p1.distance(vertices[i]) > 1e-10)
        deallog << "ERROR! d: " << p1.distance(vertices[i])
                << " - " << p1 << " != " << vertices[i] << std::endl;
    }



}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<2,2>();
  test<3,3>();

  return 0;
}

