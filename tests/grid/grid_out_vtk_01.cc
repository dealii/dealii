// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
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


// test GridOut::write_vtk

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim, int spacedim>
void
test(std::ostream &logfile)
{
  Triangulation<dim, spacedim> tria;
  std::vector<unsigned int>    legs(2 * dim, 1);
  if (dim > 1)
    GridGenerator::hyper_cross(tria, legs, true);
  else
    GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  GridOut grid_out;
  grid_out.write_vtk(tria, logfile);
}


int
main()
{
  initlog("output");
  test<1, 1>(deallog.get_file_stream());
  test<1, 2>(deallog.get_file_stream());
  test<2, 2>(deallog.get_file_stream());
  test<2, 3>(deallog.get_file_stream());
  test<3, 3>(deallog.get_file_stream());
}
