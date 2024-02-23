// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
