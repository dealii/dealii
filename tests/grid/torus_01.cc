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

// test GridTools::torus() and TorusManifold

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <fstream>
#include <iomanip>


template <int dim, int spacedim>
void test ();

template <>
void test<3,3> ()
{
  const int dim = 3;
  const int spacedim = 3;
  Triangulation<dim, spacedim> triangulation;

  GridGenerator::torus(triangulation, 1.0, 0.4);

  static const SphericalManifold<3> desc_sphere;
  static const TorusManifold<3> desc_torus(1.0, 0.4);
  triangulation.set_manifold (0, desc_torus);
  triangulation.set_manifold (1, desc_sphere);

  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement ();

  std::ofstream out ("grid-1.vtk");
  GridOut grid_out;
  grid_out.write_vtk (triangulation, out);


  GridOut().write_gnuplot (triangulation, deallog.get_file_stream());
}


int main ()
{
  initlog ();

  test<3,3> ();

  return 0;
}
