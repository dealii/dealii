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


// Test TensorProductManifold by refining and generating normals for
// a manually constructed cylinder hull.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>


#include <deal.II/grid/tensor_product_manifold.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>


void test()
{
  std::ostream &out = deallog.get_file_stream();

  FunctionManifold<1,1> F("x","x");
  SphericalManifold<2,2> G;

  TensorProductManifold<2, 1,1,1, 2,2,2> manifold(F, G);

  // make a hull of a cylinder
  Triangulation<2,3> tria;
  {
    Triangulation<3,3> volume_tria;
    GridGenerator::cylinder(volume_tria);
    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);
    GridGenerator::extract_boundary_mesh(volume_tria, tria, boundary_ids);
  }
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  tria.refine_global(1);

  out << "set view equal xyz" << std::endl
      << "splot '-' with lines, '-' with vectors " << std::endl;
  GridOut().write_gnuplot (tria, out);
  out << "e" << std::endl;

  Triangulation<2,3>::active_cell_iterator it = tria.begin_active();
  for (; it!=tria.end(); ++it)
    {
      Point<3> p = it->center(true);
      Tensor<1,3> t1 = manifold.get_tangent_vector(p, it->vertex(0));
      Tensor<1,3> t2 = manifold.get_tangent_vector(p, it->vertex(1));
      Tensor<1,3> n = cross_product_3d(t1, t2);
      n/=-n.norm();
      out << it->center() << " " << n << std::endl;
    }
  out << "e" << std::endl;
}



int main ()
{
  initlog();

  test();

  return 0;
}

