// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test cell similarity over GridTools::transform (here: scale)

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <iostream>

#include "../tests.h"



int
main()
{
  initlog();

  // there used to be a bug in the cell similarity detection beyond the
  // GridTools::transform method , but cell similarity is only enabled without
  // threads. to make sure this test is effective, manually set the thread
  // limit to 1.
  MultithreadInfo::set_thread_limit(1);

  Triangulation<2> triangulation;
  FE_DGQ<2>        fe(0);
  QMidpoint<2>     qf_cell;

  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  FEValues<2> fe_values(fe, qf_cell, update_JxW_values);

  // compute the volume of the mesh
  fe_values.reinit(triangulation.begin_active());
  const double volume_before = fe_values.JxW(0);
  deallog << volume_before << std::endl;

  // shrink the mesh
  GridTools::scale(0.5, triangulation);

  // Now we measure the volume again:
  fe_values.reinit(triangulation.begin_active());
  const double volume_after = fe_values.JxW(0);
  deallog << volume_after << std::endl;
}
