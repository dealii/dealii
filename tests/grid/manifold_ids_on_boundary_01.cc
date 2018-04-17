// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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

// Test set_manifold_ids_on_boundary(b_id,man_id)

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <iostream>


void dim_2(std::ostream &os)
{
  const unsigned int d=2;
  const Point<d> center(0,0);
  const double inner = 0.2,
               outer = 1.0;
  Triangulation<d> tr;

  GridGenerator::hyper_cube_with_cylindrical_hole(tr,inner,outer,outer,true);
  tr.set_all_manifold_ids_on_boundary(numbers::flat_manifold_id);
  tr.reset_manifold(0);
  static SphericalManifold<d> boundary;
  GridTools::copy_boundary_to_manifold_id(tr);
  tr.set_manifold(1,boundary);

  tr.refine_global(2);


  GridOut gout;
  gout.write_vtk(tr, os);
}

void dim_3(std::ostream &os)
{
  const unsigned int d=3;

  const Point<d> center(0,0,0);
  const double inner = 0.2,
               outer = 1.0;
  Triangulation<d> tr;

  GridGenerator::hyper_cube_with_cylindrical_hole(tr,inner,outer,outer,true);
  tr.set_all_manifold_ids_on_boundary(numbers::flat_manifold_id);
  tr.reset_manifold(0);
  static CylindricalManifold<d> boundary(2);
  GridTools::copy_boundary_to_manifold_id(tr);
  tr.set_manifold(1,boundary);

  tr.refine_global(1);
  GridOut gout;
  gout.write_vtk(tr, os);
}


int main()
{
  initlog(true);
  std::ostream &logfile = deallog.get_file_stream();
  dim_2(logfile);
  dim_3(logfile);
}
