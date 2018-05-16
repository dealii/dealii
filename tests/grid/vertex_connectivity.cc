// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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



// check GridTools::get_vertex_connectivity_of_cells for two different meshes

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>


template <int dim>
void
test()
{
  // First check connectivity of plain mesh
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 3);
  DynamicSparsityPattern dsp;
  GridTools::get_vertex_connectivity_of_cells(tria, dsp);
  dsp.print(deallog.get_file_stream());

  // Refine a few cells and check again
  tria.begin_active()->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  GridTools::get_vertex_connectivity_of_cells(tria, dsp);
  dsp.print(deallog.get_file_stream());
}


int
main ()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
